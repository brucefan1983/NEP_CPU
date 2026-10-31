"""Numerical regressions and independent consistency checks for the public API."""
import copy
import hashlib
import itertools
import json
import math
import unittest

from fixtures import (ALL_FLAGS, FLEXIBLE_ZBL, ZBL_COEFFICIENTS, angular_axis_frame,
                      angular_cases, angular_descriptor_oracle, cluster,
                      asymmetric_model, model_text, pair, raw_charges, read_block,
                      small_cell, supercell, write_models)
from support import (HERE, ROOT, displaced, flatten, matvec, read_matrix, read_xyz,
                     strained, tensor_transform, transform)


def configure(directory, backends):
    global WORK, BACKENDS, MODELS, REFERENCE, ANGULAR_MODELS, ANGULAR_REFERENCE
    WORK, BACKENDS = directory, backends
    MODELS = write_models(directory)
    REFERENCE = json.loads((HERE / 'references.json').read_text())
    ANGULAR_MODELS = {}
    for name, options in angular_cases().items():
        path = directory / ('angular_' + name + '.txt')
        path.write_text(model_text(**options))
        ANGULAR_MODELS[name] = path
    ANGULAR_REFERENCE = json.loads((HERE / 'angular_references.json').read_text())


class APITests(unittest.TestCase):
    def close(self, actual, expected, label, atol=2e-10, rtol=2e-10):
        """Check shape, finiteness and every component, with useful diagnostics."""
        def shape(x):
            if not isinstance(x, list):
                return ()
            return (len(x), tuple(shape(y) for y in x))
        self.assertEqual(shape(actual), shape(expected), label + ': shape mismatch')
        a, b = list(flatten(actual)), list(flatten(expected))
        self.assertTrue(all(math.isfinite(x) for x in a + b), label + ': NaN/Inf')
        errors = [(abs(x - y) / (atol + rtol * abs(y)), i, x, y) for i, (x, y) in enumerate(zip(a, b))]
        if errors:
            ratio, i, x, y = max(errors)
            self.assertLessEqual(ratio, 1, '{}: component {}, actual={:.17g}, reference={:.17g}, '
                                 '|diff|={:.3g}, allowed={:.3g}'.format(
                                     label, i, x, y, abs(x - y), atol + rtol * abs(y)))

    def compare_result(self, actual, expected, label, atol=2e-10, rtol=2e-10):
        for field, values in expected.items():
            self.assertIn(field, actual, label)
            self.close(actual[field], values, label + '/' + field, atol, rtol)

    def test_reference_inputs(self):
        """Fixture edits must be intentional; production source is free to change."""
        for path, digest in REFERENCE['input_sha256'].items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT / path).read_bytes()).hexdigest(), digest)
        for name, digest in REFERENCE['synthetic_sha256'].items():
            with self.subTest(model=name):
                self.assertEqual(hashlib.sha256(MODELS[name].read_bytes()).hexdigest(), digest)

    def test_saved_cpu_references(self):
        """All atoms/components, including all nine per-atom virial entries."""
        for mode, backend in BACKENDS.items():
            for name in ('nep', 'qnep'):
                with self.subTest(mode=mode, fixture=name):
                    directory = HERE / 'data' / name
                    actual = backend.run(directory / 'nep.txt', [read_xyz(directory)])[0]
                    for field, filename in [('force', 'force_analytical_ref.out'),
                                            ('virial', 'virial_ref.out'), ('descriptor', 'descriptor_ref.out')]:
                        expected = read_matrix(directory / filename)
                        # Analytic double precision is the original reference.
                        # Table interpolation has a distinct, frozen baseline below.
                        tol = 2e-10 if mode == 'analytic' else 5e-4
                        self.close(actual[field], expected, name + '/' + mode + '/' + field, tol, tol)
                    self.compare_result(actual, REFERENCE[mode][name], name + '/' + mode)

    def test_saved_gpu_forces(self):
        """Existing float32 GPUMD results; no GPU or GPUMD executable needed."""
        for mode, backend in BACKENDS.items():
            for name, tol in [('nep', 1e-5), ('qnep', 1.2e-4)]:
                with self.subTest(mode=mode, fixture=name):
                    directory = HERE / 'data' / name
                    actual = backend.run(directory / 'nep.txt', [read_xyz(directory)])[0]['force']
                    gpu = [row[:3] for row in read_matrix(directory / 'force_gpu.out')]
                    self.close(actual, gpu, name + '/GPU force', tol, 0)

    def test_saved_finite_difference_references(self):
        if 'analytic' not in BACKENDS:
            self.skipTest('Finite differences target the analytic implementation')
        backend = BACKENDS['analytic']
        # Three atoms per species, spanning low/high coordination and pairs
        # closest to the radial (atom 8) and angular (atom 130) cutoffs.
        atoms = (5, 8, 66, 130, 153, 238)
        for name, step in [('nep', 2e-5), ('qnep', 2e-4)]:
            with self.subTest(fixture=name):
                directory = HERE / 'data' / name
                frame = read_xyz(directory)
                frames = [displaced(frame, i, d, sign * step)
                          for i in atoms for d in range(3) for sign in (-1, 1)]
                energies = [r['total_energy'] for r in backend.run(directory / 'nep.txt', frames, 'energy')]
                fd = [(energies[k] - energies[k + 1]) / (2 * step) for k in range(0, len(energies), 2)]
                fd = [fd[i:i + 3] for i in range(0, len(fd), 3)]
                saved = read_matrix(directory / 'force_finite_difference_ref.out')
                self.close(fd, [saved[i] for i in atoms],
                           name + '/saved FD', 2e-7, 2e-7)
                # qNEP's historical 2e-4 step crosses a cutoff for some pairs.
                # Preserve that regression, but test derivative consistency at
                # 2e-6 separately instead of weakening the force assertion.
                # The closest pair is only 8.85e-6 Angstrom below the cutoff.
                if name == 'qnep':
                    step = 2e-6
                    frames = [displaced(frame, i, d, sign * step)
                              for i in atoms for d in range(3) for sign in (-1, 1)]
                    energies = [r['total_energy'] for r in backend.run(directory / 'nep.txt', frames, 'energy')]
                    fd = [(energies[k] - energies[k + 1]) / (2 * step) for k in range(0, len(energies), 2)]
                    fd = [fd[i:i + 3] for i in range(0, len(fd), 3)]
                actual = backend.run(directory / 'nep.txt', [frame])[0]['force']
                self.close(fd, [actual[i] for i in atoms], name + '/force = -dE/dr', 2e-7, 2e-7)

    def test_synthetic_references(self):
        for mode, backend in BACKENDS.items():
            for name, model in MODELS.items():
                with self.subTest(mode=mode, model=name):
                    action = 'features' if name in ('ordinary', 'typewise') else 'full'
                    self.compare_result(backend.run(model, [cluster()], action)[0],
                                        REFERENCE[mode][name], mode + '/' + name)
            with self.subTest(mode=mode, fixture='small_cell'):
                self.compare_result(backend.run(MODELS['typewise'], [small_cell()])[0],
                                    REFERENCE[mode]['small_cell'], mode + '/small_cell')

    def check_derivatives(self, backend, model, frame, all_coordinates=True, action='energy'):
        base = backend.run(model, [frame], 'd3' if action == 'd3' else 'full')[0]
        coordinates = list(itertools.product(range(len(frame['type'])), range(3)))
        if not all_coordinates:
            coordinates = [(0, 0), (2, 1), (5, 2)]
        h = 1e-5
        frames = [displaced(frame, i, d, sign * h) for i, d in coordinates for sign in (-1, 1)]
        values = [r['total_energy'] for r in backend.run(model, frames, action)]
        fd = [(values[i] - values[i + 1]) / (2 * h) for i in range(0, len(values), 2)]
        self.close(fd, [base['force'][i][d] for i, d in coordinates], 'force = -dE/dr', 3e-7, 3e-7)
        components = list(itertools.product(range(3), repeat=2))
        frames = [strained(frame, a, b, sign * h) for a, b in components for sign in (-1, 1)]
        values = [r['total_energy'] for r in backend.run(model, frames, action)]
        fd = [(values[i] - values[i + 1]) / (2 * h) for i in range(0, len(values), 2)]
        # r'_a = r_a + eps_ab r_b; virial is stored as r_a f_b.
        virial = [sum(row[b * 3 + a] for row in base['virial']) for a, b in components]
        self.close(fd, virial, 'virial = -dE/dstrain (transposed indices)', 5e-7, 5e-7)

    def test_force_and_virial_derivatives(self):
        if 'analytic' not in BACKENDS:
            self.skipTest('Interpolated values and derivatives are independent tables')
        for name in ('ordinary', 'typewise', 'charge1', 'charge2', 'zbl_fixed',
                     'zbl_typewise', 'zbl_flexible', 'zbl_charge1', 'zbl_charge2'):
            with self.subTest(model=name):
                self.check_derivatives(BACKENDS['analytic'], MODELS[name], cluster())
        with self.subTest(fixture='small periodic cell'):
            self.check_derivatives(BACKENDS['analytic'], MODELS['typewise'], small_cell())

    def test_all_many_body_flag_combinations(self):
        """64 layouts, each with energy/force/virial checks for NEP and qNEP1/2/3."""
        if 'analytic' not in BACKENDS:
            self.skipTest('Finite differences target the analytic implementation')
        backend = BACKENDS['analytic']
        for flags in itertools.product((0, 1), repeat=6):
            for kind in ('nep', 'qnep1', 'qnep2', 'qnep3'):
                with self.subTest(flags=flags, kind=kind):
                    model = WORK / 'flags.txt'
                    model.write_text(model_text(kind=kind, flags=flags))
                    actual = backend.run(model, [cluster()], 'load')[0]
                    self.assertEqual(actual['dimension'], 3 + 3 * (4 + sum(flags)))
                    self.check_derivatives(backend, model, cluster(), all_coordinates=False)

    def test_charge_and_bec(self):
        for mode, backend in BACKENDS.items():
            for name in ('charge1', 'charge2'):
                with self.subTest(mode=mode, model=name):
                    self.check_charge_and_bec(backend, MODELS[name], cluster())

    def check_charge_and_bec(self, backend, model, frame):
        base = backend.run(model, [frame])[0]
        q = [row[0] for row in base['charge']]
        self.close(sum(q), 0.0, 'charge neutrality', 2e-12, 0)
        self.close([sum(row[d] for row in base['bec']) for d in range(9)],
                   [0.0] * 9, 'BEC acoustic sum', 2e-11, 0)
        raw = raw_charges(model, base['descriptor'], frame['type'])
        mean = sum(raw) / len(raw)
        self.close(q, [x - mean for x in raw], 'ANN charge / mean subtraction')
        if backend.mode != 'analytic':
            return
        # GPUMD's bond-based BEC, not d(sum(q_i*r_i))/dr.
        # Half bond moments use derivatives of *raw* ANN charges.
        h, n = 1e-5, len(q)
        shifted = [displaced(frame, j, b, s * h) for j in range(n)
                   for b in range(3) for s in (-1, 1)]
        outputs = backend.run(model, shifted)
        raw_shifted = [raw_charges(model, x['descriptor'], frame['type']) for x in outputs]
        numerical = [[q[i] if d in (0, 4, 8) else 0. for d in range(9)] for i in range(n)]
        for j in range(n):
            for b in range(3):
                minus, plus = raw_shifted[2 * (3 * j + b):2 * (3 * j + b) + 2]
                for i in range(n):
                    if i == j:
                        continue
                    derivative = (plus[i] - minus[i]) / (2 * h)
                    for a in range(3):
                        term = 0.5 * (frame['position'][j][a] - frame['position'][i][a]) * derivative
                        numerical[i][3 * a + b] += term
                        numerical[j][3 * a + b] -= term
        self.close(base['bec'], [[1.7 * x for x in row] for row in numerical],
                   'BEC / finite differences of raw charges', 2e-8, 2e-8)

    def test_translation_permutation_and_rotation(self):
        # A proper, non-axis-aligned rotation with rational entries.
        rotation = [[0.36, -0.48, 0.8], [0.8, 0.60, 0.0], [-0.48, 0.64, 0.6]]
        frame = cluster()
        translated = copy.deepcopy(frame)
        translated['position'] = [[r[d] + (0.73, -1.21, 0.42)[d] for d in range(3)] for r in frame['position']]
        order = [3, 0, 5, 1, 4, 2]
        permuted = copy.deepcopy(frame)
        for key in ('type', 'position'):
            permuted[key] = [frame[key][i] for i in order]
        rotated = transform(frame, rotation)
        for mode, backend in BACKENDS.items():
            for name, model in MODELS.items():
                with self.subTest(mode=mode, model=name):
                    base, moved, perm, rot = backend.run(model, [frame, translated, permuted, rotated])
                    self.compare_result(moved, base, 'translation', 3e-9, 3e-9)
                    for field in ('energy', 'force', 'virial', 'charge', 'bec', 'descriptor'):
                        if field in base:
                            self.close(perm[field], [base[field][i] for i in order], 'permutation/' + field, 3e-9, 3e-9)
                    for field in ('energy', 'charge', 'descriptor'):
                        if field in base:
                            self.close(rot[field], base[field], 'rotation invariant/' + field, 3e-9, 3e-9)
                    for field in ('virial', 'bec'):
                        if field in base:
                            self.close(rot[field], [tensor_transform(row, rotation) for row in base[field]],
                                       'tensor rotation/' + field, 3e-9, 3e-9)
                    if 'force' in base:
                        self.close(rot['force'], [matvec(rotation, row) for row in base['force']], 'force rotation', 3e-9, 3e-9)
                        self.close([sum(row[d] for row in base['force']) for d in range(3)],
                                   [0.] * 3, 'zero net force', 2e-10, 0)
                    if 'response' in base:
                        self.close(perm['response'], base['response'], 'response permutation', 3e-9, 3e-9)
                        response = base['response'][0]
                        if name == 'dipole':
                            expected = matvec(rotation, response)
                        else:
                            xx, yy, zz, xy, yz, zx = response
                            rr = tensor_transform([xx, xy, zx, xy, yy, yz, zx, yz, zz], rotation)
                            expected = [rr[d] for d in (0, 4, 8, 1, 5, 6)]
                        self.close(rot['response'], [expected], 'response rotation', 3e-9, 3e-9)

    def test_periodic_images_and_supercell(self):
        for mode, backend in BACKENDS.items():
            for name in ('ordinary', 'typewise', 'charge1', 'charge2'):
                with self.subTest(mode=mode, model=name):
                    frame = small_cell()
                    wrapped = copy.deepcopy(frame)
                    # Different lattice translations for different atoms.
                    for i, r in enumerate(wrapped['position']):
                        shift = matvec(frame['box'], [i % 3 - 1, (i + 1) % 3 - 1, (i + 2) % 3 - 1])
                        wrapped['position'][i] = [x + y for x, y in zip(r, shift)]
                    base, moved = backend.run(MODELS[name], [frame, wrapped])
                    self.compare_result(moved, base, 'periodic image', 3e-9, 3e-9)
                    if name in ('ordinary', 'typewise'):
                        # Small-box explicit images vs a supercell large enough
                        # to enter the cell-list neighbor search.
                        large = backend.run(MODELS[name], [supercell(frame)])[0]
                        for field in ('energy', 'force', 'virial', 'descriptor'):
                            self.close(large[field], base[field] * 27, 'supercell/' + field, 3e-9, 3e-9)

    def test_reuse_and_atom_count_changes(self):
        for mode, backend in BACKENDS.items():
            for name in ('ordinary', 'typewise', 'charge1', 'charge2'):
                with self.subTest(mode=mode, model=name):
                    a, b, c = pair((0, 1), 2.3), cluster(), small_cell()
                    sequence = [a, b, c, a, b]
                    shared = backend.run(MODELS[name], sequence)
                    fresh = [backend.run(MODELS[name], [frame])[0] for frame in sequence[:3]]
                    for result, i in zip(shared, (0, 1, 2, 0, 1)):
                        self.compare_result(result, fresh[i], 'same object / fresh object')

    def test_typewise_cutoffs(self):
        for mode, backend in BACKENDS.items():
            for types in itertools.product((0, 1), repeat=2):
                with self.subTest(mode=mode, types=types):
                    # A type-pair cutoff is the arithmetic mean of its two
                    # element cutoffs. With only this pair, a uniform model
                    # at that cutoff is an independent equivalent input.
                    radial = sum((5.0, 4.4)[t] for t in types) / 2
                    angular = sum((4.0, 3.6)[t] for t in types) / 2
                    uniform = WORK / 'pair_uniform.txt'
                    uniform.write_text(model_text().replace('cutoff 5 4 100 100',
                        'cutoff {} {} 100 100'.format(radial, angular)))
                    distances = [1.4, angular - 1e-5, angular, angular + 1e-5,
                                 radial - 1e-5, radial, radial + 1e-5, radial + 0.2]
                    frames = [pair(types, r) for r in distances]
                    actual = backend.run(MODELS['typewise'], frames)
                    expected = backend.run(uniform, frames)
                    for r, a, b in zip(distances, actual, expected):
                        self.compare_result(a, b, 'pair cutoff at r={}'.format(r))
                    for result in actual[-3:]:
                        self.close(result['force'], [[0.] * 3] * 2, 'outside radial cutoff', 2e-12, 0)

    def test_zbl_pair_oracle(self):
        # Zero ANN isolates the pair potential. This Python expression is
        # independent of NEP's force and virial kernels, including the switch.
        for mode, backend in BACKENDS.items():
            for kind in ('fixed', 'typewise', 'flexible'):
                zbl = {'fixed': 'zbl 1 3.5', 'typewise': 'zbl 2 4 0.65', 'flexible': 'zbl 0 0'}[kind]
                model = WORK / ('zero_zbl_' + kind + '.txt')
                model.write_text(model_text(zbl=zbl, zero=True))
                for types in itertools.product((0, 1), repeat=2):
                    inner, outer, coeff = 1., 3.5, ZBL_COEFFICIENTS
                    if kind == 'typewise':
                        inner, outer = 0., min(4., 0.65 * sum((1., 1.38667)[t] for t in types))
                    if kind == 'flexible':
                        inner, outer, *coeff = FLEXIBLE_ZBL[sum(types)]
                    distances = sorted(set([0.9, 2.7, inner + 1e-6, outer - 1e-6, outer, outer + 1e-6]))
                    distances = [r for r in distances if r > 0.1]
                    outputs = backend.run(model, [pair(types, r) for r in distances])
                    for r, result in zip(distances, outputs):
                        with self.subTest(mode=mode, zbl=kind, types=types, r=r):
                            energy, derivative = 0., 0.
                            if r < outer:
                                z1, z2 = [(6, 14)[t] for t in types]
                                screening = (z1 ** 0.23 + z2 ** 0.23) * 2.134563
                                terms = [coeff[k] * math.exp(-coeff[k + 1] * screening * r) for k in (0, 2, 4, 6)]
                                phi = sum(terms)
                                dp = sum(-coeff[2 * k + 1] * screening * x for k, x in enumerate(terms))
                                switch, ds = 1., 0.
                                if r >= inner:
                                    w = math.pi / (outer - inner)
                                    switch, ds = 0.5 * (1 + math.cos(w * (r - inner))), -0.5 * w * math.sin(w * (r - inner))
                                k = 14.399645 * z1 * z2
                                energy = k * phi * switch / r
                                derivative = k * ((dp / r - phi / r ** 2) * switch + phi * ds / r)
                            direction = (0.36, 0.48, 0.8)
                            force = [[sign * derivative * u for u in direction] for sign in (1, -1)]
                            virial = [[-0.5 * r * derivative * a * b for a in direction for b in direction]] * 2
                            self.close(result['energy'], [[0.5 * energy]] * 2, 'ZBL atomic energies', 3e-10, 3e-11)
                            self.close(result['force'], force, 'ZBL force', 3e-10, 3e-11)
                            self.close(result['virial'], virial, 'ZBL per-atom virial', 3e-10, 3e-11)

    def test_dftd3(self):
        model = ROOT / 'tests/api/data/dftd3/nep.txt'
        frame = cluster()
        frame['type'] = [0] * len(frame['type'])
        for mode, backend in BACKENDS.items():
            with self.subTest(mode=mode):
                base = backend.run(model, [frame])[0]
                d3 = backend.run(model, [frame], 'd3')[0]
                combined = backend.run(model, [frame], 'with_d3')[0]
                self.compare_result(d3, REFERENCE[mode]['d3'], 'DFTD3 snapshot')
                for field in ('energy', 'force', 'virial'):
                    expected = [[x + y for x, y in zip(a, b)] for a, b in zip(base[field], d3[field])]
                    self.close(combined[field], expected, 'NEP + DFTD3/' + field)
                if mode == 'analytic':
                    self.check_derivatives(backend, model, frame, action='d3')

    def test_latent_and_parameter_gradients(self):
        for mode, backend in BACKENDS.items():
            with self.subTest(mode=mode):
                model, frame = MODELS['ordinary'], cluster()
                base = backend.run(model, [frame], 'features')[0]
                lines = model.read_text().splitlines()
                dim, neurons, h = 33, 8, 1e-5
                start = next(i for i, line in enumerate(lines) if line.startswith('@ parameters ')) + 1
                common_bias = float(lines[start + 2 * neurons * (dim + 2)])
                self.close([[sum(row) - common_bias] for row in base['latent']], base['energy'], 'latent reconstruction')
                # Perturb input weight, hidden bias and output weight for both
                # species. B_projection is per-atom, neuron-major within atom.
                for t, j, component in itertools.product((0, 1), (0, 7), (0, dim, dim + 1)):
                    with self.subTest(type=t, neuron=j, parameter=component):
                        offset = j * dim if component == 0 else neurons * dim + j
                        if component == dim + 1:
                            offset += neurons
                        index = start + t * neurons * (dim + 2) + offset
                        energies = []
                        for sign in (-1, 1):
                            changed = lines[:]
                            changed[index] = format(float(lines[index]) + sign * h, '.17g')
                            path = WORK / 'perturbed_parameter.txt'
                            path.write_text('\n'.join(changed) + '\n')
                            energies.append(backend.run(path, [frame])[0]['energy'])
                        derivative = [(p[0] - m[0]) / (2 * h) for m, p in zip(*energies)]
                        expected = [row[j * (dim + 2) + component] if tt == t else 0.
                                    for tt, row in zip(frame['type'], base['projection'])]
                        self.close(derivative, expected, 'B_projection parameter derivative', 2e-8, 2e-8)

    def test_rejected_models_and_api_misuse(self):
        for mode, backend in BACKENDS.items():
            cases = []
            for tag in ('nep3', 'nep3_zbl', 'nep3_dipole', 'nep3_polarizability',
                        'nep4', 'nep4_zbl', 'nep4_charge1', 'nep4_charge2',
                        'nep4_zbl_charge1', 'nep4_zbl_charge2', 'nep4_dipole', 'nep4_polarizability',
                        'nep5', 'nep5_zbl', 'nep4_charge3', 'nep4_zbl_charge3', 'unknown'):
                path = WORK / (tag + '.txt')
                path.write_text(tag + ' 2 C Si\n')
                cases.append((path, 'load', 'expected first line'))
            for name, text, expected in [
                ('two_layers', model_text().replace('ann 8 0', 'ann 8 4'), 'one hidden layer'),
                ('short_ann', model_text().replace('ann 8 0', 'ann 8'), "incorrect number of values in 'ann'"),
                ('bad_count', 'nep 3 C Si\n', 'element count'),
                ('short_header', 'nep\n', 'expected first line')]:
                path = WORK / (name + '.txt')
                path.write_text(text)
                cases.append((path, 'load', expected))
            cases.append((WORK / 'does_not_exist.txt', 'load', 'Failed to open'))
            cases += [(MODELS['ordinary'], 'charge', 'qNEP model'),
                      (MODELS['charge1'], 'ordinary', 'qNEP model'),
                      (MODELS['ordinary'], 'dipole', 'non-dipole'),
                      (MODELS['ordinary'], 'polarizability', 'non-polarizability')]
            for name in ('dipole', 'polarizability', 'charge1', 'charge2'):
                cases.append((MODELS[name], 'bad_lammps', 'only supports non-charge'))
            for name in ('dipole', 'polarizability'):
                cases.append((MODELS[name], 'ordinary', 'non-potential'))
            for name in ('ordinary', 'charge1'):
                for field, message in [('position', 'position'), ('energy', 'potential'), ('force', 'force'), ('virial', 'virial')]:
                    cases.append((MODELS[name], 'bad_' + field, 'Type and ' + message + ' sizes'))
            cases += [(MODELS['charge1'], 'bad_charge', 'Type and charge sizes'),
                      (MODELS['charge1'], 'bad_bec', 'Type and BEC sizes')]
            for model, action, message in cases:
                with self.subTest(mode=mode, model=model.name, action=action):
                    result, output = backend.invoke(model, [cluster()], action)
                    self.assertEqual(result.returncode, 1, 'Expected controlled exit(1), not a crash: ' + result.stdout + result.stderr)
                    self.assertIn(message, result.stdout + result.stderr)
                    self.assertEqual(output, '', 'Rejected input produced numerical results')

    def test_q222_requires_zero_or_two(self):
        for mode, backend in BACKENDS.items():
            with self.subTest(mode=mode):
                path = WORK / 'q222_marker_1.txt'
                path.write_text(model_text().replace('l_max 4 2', 'l_max 4 1'))
                result, output = backend.invoke(path, [cluster()], 'load')
                self.assertEqual(result.returncode, 1)
                self.assertIn('q_222 value in l_max must be 0 or 2', result.stdout)
                self.assertEqual(output, '')

    def test_model_record_syntax(self):
        """Record order, case, comments and line grouping do not change a model."""
        for mode, backend in BACKENDS.items():
            for name, original in MODELS.items():
                with self.subTest(mode=mode, model=name):
                    lines = original.read_text().splitlines()
                    records, i = [], 1
                    while i < len(lines):
                        tokens = lines[i].split()
                        if tokens[0] == '@':
                            count = int(tokens[2])
                            values = lines[i + 1:i + 1 + count]
                            record = [lines[i].upper()] + [
                                '\t'.join(values[j:j + 7]) + ' # values'
                                for j in range(0, count, 7)]
                            i += count + 1
                        else:
                            record = [lines[i].upper() + ' # field']
                            i += 1
                        records.append('\r\n'.join(record))
                    text = lines[0].upper() + ' # elements\r\n\r\n'
                    text += '\r\n\r\n# next record\r\n'.join(reversed(records))
                    changed = WORK / 'reordered.txt'
                    changed.write_bytes(text.encode('ascii'))  # CRLF, no final newline
                    base = backend.run(original, [cluster()])[0]
                    self.compare_result(backend.run(changed, [cluster()])[0], base, 'record syntax')
                    stored = backend.run(changed, [cluster()], 'layout')[0]
                    self.assertEqual(stored['parameters'][0], read_block(original, 'parameters'))

    def test_invalid_model_records(self):
        text = model_text()
        cases = {
            'blank_first_line': '\n' + text,
            'bom': '\ufeff' + text,
            'duplicate_element': text.replace('C Si', 'C c'),
            'unknown_element': text.replace('C Si', 'C Xx'),
            'zero_elements': text.replace('nep 2 C Si', 'nep 0'),
            'version': text.replace('format_version 1', 'format_version 2'),
            'duplicate_field': text + 'ANN 8 0\n',
            'duplicate_block': text + '@ PARAMETERS 0\n',
            'field_block_collision': text + '@ ANN 0\n',
            'unknown_field': text + 'unknown 0\n',
            'unknown_block': text + '@ unknown 1\n0\n',
            'block_header': text.replace('@ parameters ', '@ parameters 1 '),
            'negative_count': text + '@ unknown -1\n',
            'huge_count': text + '@ unknown 999999999999999999999\n',
            'truncated_block': text.rsplit('\n', 2)[0] + '\n',
            'extra_block_value': text.rstrip() + ' 0\n',
            'wrong_scaler_count': text.replace('@ q_scaler 33\n', '@ q_scaler 34\n') + '1\n',
            'wrong_parameter_count': text.replace('@ parameters 633\n', '@ parameters 634\n0\n'),
            'wrong_kind': text.replace('energy_reference 0 0', '@ energy_reference 2\n0 0'),
            'cutoff_arity': text.replace('cutoff 5 4 100 100', 'cutoff 5 4 5 100 100'),
            'cutoff_order': text.replace('cutoff 5 4', 'cutoff 4 5'),
            'cutoff_low': text.replace('cutoff 5 4', 'cutoff 5 2.9'),
            'cutoff_high': text.replace('cutoff 5 4', 'cutoff 101 4'),
            'negative_neighbors': text.replace('100 100', '-1 100'),
            'n_max': text.replace('n_max 2 2', 'n_max 16 2'),
            'basis_size': text.replace('basis_size 2 2', 'basis_size 2 13'),
            'dimension': text.replace('n_max 2 2', 'n_max 2 12'),
            'l_max': text.replace('l_max 4 2', 'l_max 9 2'),
            'flag': text.replace('l_max 4 2 1', 'l_max 4 2 2'),
            'ann_zero': text.replace('ann 8 0', 'ann 0 0'),
            'ann_large': text.replace('ann 8 0', 'ann 121 0'),
            'integer_suffix': text.replace('ann 8 0', 'ann 8x 0'),
            'reference_count': text.replace('energy_reference 0 0', 'energy_reference 0'),
            'reference_underflow': text.replace('energy_reference 0 0', 'energy_reference 1e-999 0'),
            'reference_overflow': text.replace('energy_reference 0 0', 'energy_reference 1e999 0'),
            'zbl_order': text + 'zbl 2 1\n',
            'zbl_factor': text + 'zbl 1 3 0.4\n',
            'flexible_missing': text + 'zbl 0 0\n',
            'flexible_factor': model_text(zbl='zbl 0 0').replace('zbl 0 0\n', 'zbl 0 0 0.8\n'),
            'flexible_radii': model_text(zbl='zbl 0 0').replace('@ zbl_parameters 30\n0.6\n1.5',
                                                                            '@ zbl_parameters 30\n1.5\n0.6'),
            'charge_typewise': model_text(kind='qnep1', typewise=True),
            'charge_flexible': model_text(kind='qnep2', zbl='zbl 0 0'),
            'response_reference': model_text(kind='nep_dipole') + 'energy_reference 0 0\n',
            'response_zbl': model_text(kind='nep_polarizability', zbl='zbl 1 3'),
        }
        for field in ('format_version', 'model', 'cutoff', 'n_max', 'basis_size',
                      'l_max', 'ann', 'energy_reference'):
            cases['missing_' + field] = '\n'.join(
                line for line in text.splitlines() if not line.startswith(field + ' ')) + '\n'
        for block in ('parameters', 'q_scaler'):
            lines = text.splitlines()
            start = next(i for i, line in enumerate(lines) if line.startswith('@ ' + block + ' '))
            count = int(lines[start].split()[2])
            cases['missing_' + block] = '\n'.join(lines[:start] + lines[start + count + 1:]) + '\n'
        for name in ('nep_temperature', 'nep_spin', 'nep_vdw', 'unknown'):
            cases['model_' + name] = text.replace('model nep\n', 'model ' + name + '\n')
        for value in ('nan', 'inf', '-inf', '1e40', '1e-60', '1e-999', '0x1p0', '1x', '1,2'):
            scaler_start = text.index('@ q_scaler 33\n') + len('@ q_scaler 33\n')
            scaler_end = text.index('\n', scaler_start)
            cases['number_' + value] = text[:scaler_start] + value + text[scaler_end:]
        cases['control_character'] = text.replace('ann 8 0', 'ann\x008 0')
        for mode, backend in BACKENDS.items():
            for name, content in cases.items():
                with self.subTest(mode=mode, case=name):
                    path = WORK / 'invalid_record.txt'
                    path.write_text(content)
                    result, output = backend.invoke(path, [cluster()], 'load')
                    self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                    self.assertIn(str(path), result.stdout + result.stderr)
                    self.assertEqual(output, '')

    def test_energy_reference(self):
        references = [1.23456789012345, -2.34567890123456]
        frame = cluster()
        for mode, backend in BACKENDS.items():
            for name, original in MODELS.items():
                if name in ('dipole', 'polarizability'):
                    continue
                actions = ['features', 'd3', 'with_d3'] if name == 'ordinary' else ['full']
                for action in actions:
                    with self.subTest(mode=mode, model=name, action=action):
                        path = WORK / 'energy_reference.txt'
                        path.write_text(original.read_text().replace('energy_reference 0 0',
                            'energy_reference ' + ' '.join(map(repr, references))))
                        expected = backend.run(original, [frame], action)[0]
                        if action != 'd3':
                            shifts = [references[t] for t in frame['type']]
                            expected['total_energy'] += sum(shifts)
                            expected['energy'] = [[row[0] + shift]
                                                  for row, shift in zip(expected['energy'], shifts)]
                        self.compare_result(backend.run(path, [frame], action)[0], expected,
                                            'only atomic energies change')
            # Reference energies are float64, unlike float32-range model blocks.
            path = WORK / 'large_reference.txt'
            path.write_text(model_text().replace('energy_reference 0 0', 'energy_reference 1e40 1e-60'))
            self.assertEqual(backend.run(path, [frame], 'load')[0]['dimension'], 33)

    def check_lammps(self, backend, model, frame):
        cpu = backend.run(model, [frame])[0]
        expected = {k: cpu[k] for k in ('total_energy', 'energy', 'force', 'virial')}
        expected['total_virial'] = [[sum(row[d] for row in cpu['virial'])
                                     for d in (0, 4, 8, 1, 2, 5)]]
        for action, fields in [('lammps', ('energy', 'virial')),
                               ('lammps_energy', ('energy',)),
                               ('lammps_virial', ('virial',)), ('lammps_total', ())]:
            with self.subTest(action=action):
                selected = {k: expected[k] for k in ('total_energy', 'total_virial', 'force') + fields}
                actual = backend.run(model, [frame], action)[0]
                self.assertEqual(set(actual), set(selected))
                self.compare_result(actual, selected, 'LAMMPS reversed type map / optional buffers')

    def test_lammps_energy_reference(self):
        for mode, backend in BACKENDS.items():
            for name in ('ordinary', 'typewise', 'zbl_fixed', 'zbl_typewise', 'zbl_flexible'):
                for reference in ('0 0', '1.23456789012345 -2.34567890123456'):
                    with self.subTest(mode=mode, model=name, reference=reference):
                        path = WORK / 'lammps_reference.txt'
                        path.write_text(MODELS[name].read_text().replace('energy_reference 0 0',
                                                                        'energy_reference ' + reference))
                        self.check_lammps(backend, path, cluster())

    def test_rank_diagnostics(self):
        """Only rank zero prints; every failing rank must still exit with 1."""
        legacy = WORK / 'rank_legacy.txt'
        legacy.write_text('nep4 2 C Si\n')
        malformed = WORK / 'rank_malformed.txt'
        malformed.write_text(model_text().replace('ann 8 0', 'ann 8'))
        cases = [(WORK / 'rank_missing.txt', 'load', 'Failed to open'),
                 (legacy, 'load', 'expected first line'),
                 (malformed, 'load', "incorrect number of values in 'ann'"),
                 (MODELS['ordinary'], 'bad_type_map', 'There is no element H'),
                 (MODELS['charge1'], 'bad_lammps', 'only supports non-charge')]
        for mode, backend in BACKENDS.items():
            expected = backend.run(MODELS['ordinary'], [cluster()])[0]
            for rank in (0, 1, 3):
                for model, action, message in cases:
                    with self.subTest(mode=mode, rank=rank, model=model.name, action=action):
                        result, output = backend.invoke(model, [cluster()], action, rank=rank)
                        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                        self.assertEqual(output, '')
                        if rank == 0:
                            self.assertEqual((result.stdout + result.stderr).count(message), 1)
                        else:
                            self.assertEqual(result.stdout + result.stderr, '')
                with self.subTest(mode=mode, rank=rank, case='valid model'):
                    result, output = backend.invoke(MODELS['ordinary'], [cluster()], rank=rank)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    # The rank flag controls diagnostics, not the numerical output.
                    self.compare_result(json.loads(output), expected, 'rank flag')
                    if rank != 0:
                        self.assertEqual(result.stdout + result.stderr, '')

    def test_asymmetric_coefficient_layout(self):
        """Check memory order and pair descriptors independently of the C++ index helper."""
        frame = cluster()
        frame['type'] = [0, 1, 2, 0, 1, 2]

        def radial_sum(coefficients, distance, cutoff):
            if distance >= cutoff:
                return 0.0
            angle = math.acos(2 * (distance / cutoff - 1) ** 2 - 1)
            fc = (1 + math.cos(math.pi * distance / cutoff)) / 2
            return sum(c * (1 + math.cos(k * angle)) * fc / 2
                       for k, c in enumerate(coefficients))

        for n_max, basis_size in [((1, 2), (4, 1)), ((4, 0), (1, 3))]:
            for kind in ('nep', 'qnep1', 'qnep2', 'qnep3', 'nep_dipole', 'nep_polarizability'):
                text, coefficients = asymmetric_model(kind, n_max, basis_size)
                model = WORK / 'asymmetric.txt'
                model.write_text(text)
                scalers = read_block(model, 'q_scaler')
                parameters = read_block(model, 'parameters')
                c_count = 9 * sum((n + 1) * (b + 1) for n, b in zip(n_max, basis_size))
                for segment in coefficients:
                    self.assertNotEqual(segment[0][1], segment[1][0])
                # These distances lie on both table grids (cutoffs 5 and 4).
                pairs = [pair(types, distance) for types in itertools.product(range(3), repeat=2)
                         for distance in (2.0, 3.5, 4.5)]
                for mode, backend in BACKENDS.items():
                    with self.subTest(mode=mode, model=kind, n_max=n_max, basis_size=basis_size):
                        stored = backend.run(model, [frame], 'layout')[0]
                        self.assertEqual(stored['parameters'][0], parameters)
                        self.assertEqual(stored['coefficient_offset'], len(parameters) - c_count)
                        for geometry, result in zip(pairs, backend.run(model, pairs)):
                            distance = math.sqrt(sum((b - a) ** 2 for a, b in zip(*geometry['position'])))
                            for i, t1 in enumerate(geometry['type']):
                                t2 = geometry['type'][1 - i]
                                radial, angular = [
                                    [radial_sum(row, distance, cutoff) for row in segment[t1][t2]]
                                    for segment, cutoff in zip(coefficients, (5.0, 4.0))]
                                # Spherical-harmonic addition theorem for a single neighbor.
                                expected = radial + [(2 * l + 1) / (4 * math.pi) * g * g
                                                     for l in range(1, 5) for g in angular]
                                expected = [q * scale for q, scale in zip(expected, scalers)]
                                self.close(result['descriptor'][i][:len(expected)], expected,
                                           'ordered-pair radial and angular descriptors')
                        if kind == 'nep':
                            self.check_lammps(backend, model, frame)
                        if kind.startswith('qnep'):
                            self.check_charge_and_bec(backend, model, frame)
                        if mode == 'analytic' and kind in ('nep', 'qnep1', 'qnep2', 'qnep3'):
                            self.check_derivatives(backend, model, frame, all_coordinates=False)
                        if backend.openmp:
                            one = backend.run(model, [frame])[0]
                            many = backend.run(model, [frame], threads=3)[0]
                            self.compare_result(many, one, 'asymmetric layout / OpenMP', 3e-9, 3e-9)

    def test_angular_references(self):
        self.assertEqual(set(ANGULAR_REFERENCE['model_sha256']), set(ANGULAR_MODELS))
        frame = cluster()
        encoded = json.dumps(frame, sort_keys=True, separators=(',', ':')).encode()
        self.assertEqual(hashlib.sha256(encoded).hexdigest(), ANGULAR_REFERENCE['frame_sha256'])
        for name, model in ANGULAR_MODELS.items():
            with self.subTest(model=name):
                self.assertEqual(hashlib.sha256(model.read_bytes()).hexdigest(),
                                 ANGULAR_REFERENCE['model_sha256'][name])
        for mode, backend in BACKENDS.items():
            self.assertEqual(set(ANGULAR_REFERENCE[mode]), set(ANGULAR_MODELS))
            for name, model in ANGULAR_MODELS.items():
                with self.subTest(mode=mode, model=name):
                    self.compare_result(backend.run(model, [frame])[0],
                                        ANGULAR_REFERENCE[mode][name], 'angular reference/' + name)

    def test_angular_descriptor_oracle(self):
        if 'analytic' not in BACKENDS:
            self.skipTest('Exact radial functions target the analytic implementation')
        backend = BACKENDS['analytic']
        for name, model in ANGULAR_MODELS.items():
            if not name.startswith(('nep_l', 'typewise_l')):
                continue
            for frame in (cluster(), angular_axis_frame()):
                with self.subTest(model=name, positions=frame['position']):
                    expected = angular_descriptor_oracle(model, frame)
                    actual = backend.run(model, [frame])[0]['descriptor']
                    self.close([row[3:3 + len(expected[0])] for row in actual], expected,
                               'three-body descriptors / Legendre addition theorem')

    def test_angular_derivatives(self):
        if 'analytic' not in BACKENDS:
            self.skipTest('Finite differences target the analytic implementation')
        backend = BACKENDS['analytic']
        for name, model in ANGULAR_MODELS.items():
            if 'dipole' in name or 'polarizability' in name:
                continue
            with self.subTest(model=name):
                self.check_derivatives(backend, model, cluster(), all_coordinates=False)
        for l_max in range(2, 9):
            with self.subTest(l_max=l_max, flags='three-body only'):
                path = WORK / 'angular_three_body.txt'
                path.write_text(model_text(l_max=l_max, flags=(0,) * 6))
                self.check_derivatives(backend, path, cluster(), all_coordinates=False)
        for l_max in (2, 4, 5, 8):
            with self.subTest(l_max=l_max, geometry='axes and near pole'):
                self.check_derivatives(backend, ANGULAR_MODELS['nep_l{}'.format(l_max)],
                                       angular_axis_frame(), all_coordinates=False)

    def test_angular_symmetries(self):
        rotation = [[0.36, -0.48, 0.8], [0.8, 0.60, 0.0], [-0.48, 0.64, 0.6]]
        frame = cluster()
        rotated = transform(frame, rotation)
        order = [3, 0, 5, 1, 4, 2]
        permuted = copy.deepcopy(frame)
        for field in ('type', 'position'):
            permuted[field] = [frame[field][i] for i in order]
        for mode, backend in BACKENDS.items():
            for name, model in ANGULAR_MODELS.items():
                with self.subTest(mode=mode, model=name):
                    base, rot, perm = backend.run(model, [frame, rotated, permuted])
                    for field in ('energy', 'charge', 'descriptor'):
                        if field in base:
                            self.close(rot[field], base[field], 'angular rotation/' + field, 3e-9, 3e-9)
                    for field in ('virial', 'bec'):
                        if field in base:
                            self.close(rot[field], [tensor_transform(row, rotation) for row in base[field]],
                                       'angular tensor rotation/' + field, 3e-9, 3e-9)
                    if 'force' in base:
                        self.close(rot['force'], [matvec(rotation, row) for row in base['force']],
                                   'angular force rotation', 3e-9, 3e-9)
                        self.close([sum(row[d] for row in base['force']) for d in range(3)],
                                   [0.] * 3, 'angular zero net force', 2e-10, 0)
                    for field in ('energy', 'force', 'virial', 'charge', 'bec', 'descriptor'):
                        if field in base:
                            self.close(perm[field], [base[field][i] for i in order],
                                       'angular permutation/' + field, 3e-9, 3e-9)
                    if 'response' in base:
                        self.close(perm['response'], base['response'], 'angular response permutation', 3e-9, 3e-9)
                        response = base['response'][0]
                        if name.startswith('nep_dipole'):
                            expected = matvec(rotation, response)
                        else:
                            xx, yy, zz, xy, yz, zx = response
                            tensor = tensor_transform([xx, xy, zx, xy, yy, yz, zx, yz, zz], rotation)
                            expected = [tensor[d] for d in (0, 4, 8, 1, 5, 6)]
                        self.close(rot['response'], [expected], 'angular response rotation', 3e-9, 3e-9)

    def test_angular_lammps(self):
        for mode, backend in BACKENDS.items():
            for name, model in ANGULAR_MODELS.items():
                if not name.startswith(('nep_l', 'typewise_l', 'zbl_')):
                    continue
                with self.subTest(mode=mode, model=name):
                    self.check_lammps(backend, model, cluster())
            for l_max in (2, 5, 8):
                with self.subTest(mode=mode, l_max=l_max, geometry='axes and near pole'):
                    self.check_lammps(backend, ANGULAR_MODELS['nep_l{}'.format(l_max)], angular_axis_frame())

    def test_angular_charge_and_bec(self):
        for mode, backend in BACKENDS.items():
            for name, model in ANGULAR_MODELS.items():
                if name.startswith('qnep'):
                    with self.subTest(mode=mode, model=name):
                        self.check_charge_and_bec(backend, model, cluster())

    def test_angular_reuse_and_periodic_images(self):
        for mode, backend in BACKENDS.items():
            for name in ('typewise_l2', 'typewise_l5', 'typewise_l8', 'qnep3_l8'):
                with self.subTest(mode=mode, model=name):
                    model = ANGULAR_MODELS[name]
                    frames = [pair((1, 0), 2.3), cluster(), small_cell(), angular_axis_frame()]
                    fresh = [backend.run(model, [frame])[0] for frame in frames]
                    order = (0, 1, 2, 3, 0, 2, 1)
                    shared = backend.run(model, [frames[i] for i in order])
                    for value, i in zip(shared, order):
                        self.compare_result(value, fresh[i], 'angular reused object')
                    shifted = copy.deepcopy(frames[2])
                    for i, r in enumerate(shifted['position']):
                        image = (i % 3 - 1, 1 - i % 2, i % 2)
                        delta = matvec(shifted['box'], image)
                        shifted['position'][i] = [x + y for x, y in zip(r, delta)]
                    self.compare_result(backend.run(model, [shifted])[0], fresh[2],
                                        'angular periodic images', 3e-9, 3e-9)
                    if name.startswith('typewise'):
                        large, small = backend.run(model, [supercell(frames[2]), frames[2]])
                        self.compare_result(small, fresh[2], 'angular shrink after supercell', 3e-9, 3e-9)
                        for field in ('energy', 'force', 'virial', 'descriptor'):
                            self.close(large[field], fresh[2][field] * 27,
                                       'angular small cell / supercell/' + field, 3e-9, 3e-9)

    def test_angular_openmp(self):
        if not any(b.openmp for b in BACKENDS.values()):
            self.skipTest('OpenMP disabled with --no-openmp')
        for mode, backend in BACKENDS.items():
            for name, model in ANGULAR_MODELS.items():
                with self.subTest(mode=mode, model=name):
                    frames = [cluster(), small_cell()]
                    one = backend.run(model, frames, threads=1)
                    many = backend.run(model, frames, threads=3)
                    for a, b in zip(one, many):
                        self.compare_result(a, b, 'angular 1 vs 3 threads', 3e-9, 3e-9)

    def test_qnep3_self_energy(self):
        """Compare against qNEP2 plus the independent Ewald self-energy formula."""
        q2, q3 = WORK / 'self_qnep2.txt', WORK / 'self_qnep3.txt'
        for cutoff, zbl in [(5.0, None), (6.0, 'zbl 1 3.5')]:
            text = model_text(kind='qnep2', zbl=zbl).replace('cutoff 5 4', 'cutoff {} 4'.format(cutoff))
            text = text.replace('energy_reference 0 0', 'energy_reference 1.23456789012345 -2.34567890123456')
            q2.write_text(text)
            q3.write_text(text.replace('model qnep2', 'model qnep3'))
            coefficient = -14.399645 * math.sqrt(math.pi) / cutoff
            for mode, backend in BACKENDS.items():
                for frame in (cluster(), small_cell()):
                    with self.subTest(mode=mode, cutoff=cutoff, box=frame['box']):
                        before = backend.run(q2, [frame])[0]
                        after = backend.run(q3, [frame])[0]
                        for field in ('charge', 'bec', 'descriptor'):
                            self.close(after[field], before[field], 'qNEP3 unchanged ' + field)
                        shifts = [coefficient * q[0] ** 2 for q in before['charge']]
                        self.assertLess(sum(shifts), -1e-4)
                        self.close(after['energy'], [[e[0] + s] for e, s in zip(before['energy'], shifts)],
                                   'qNEP3 atomic self energy')
                        self.close(after['total_energy'], before['total_energy'] + sum(shifts),
                                   'qNEP3 total self energy')
                        if mode == 'analytic':
                            self.check_derivatives(backend, q3, frame)
                        if backend.openmp:
                            self.compare_result(backend.run(q3, [frame], threads=3)[0], after,
                                                'qNEP3 / OpenMP', 3e-9, 3e-9)
        # With zero ANN charge, the new term must vanish exactly.
        q2.write_text(model_text(kind='qnep2', zero=True))
        q3.write_text(model_text(kind='qnep3', zero=True))
        for mode, backend in BACKENDS.items():
            self.compare_result(backend.run(q3, [cluster()])[0], backend.run(q2, [cluster()])[0],
                                mode + '/zero charge')

    def test_qnep3_rejections(self):
        text = model_text(kind='qnep3')
        cases = [(model_text(kind='qnep3', typewise=True), 'load', 'uniform cutoffs'),
                 (model_text(kind='qnep3', zbl='zbl 0 0'), 'load', 'Flexible ZBL'),
                 (text.replace('ann 8 0', 'ann 8 4'), 'load', 'one hidden layer'),
                 (text, 'bad_lammps', 'only supports non-charge'),
                 (text, 'ordinary', 'qNEP')]
        model = WORK / 'rejected_qnep3.txt'
        for mode, backend in BACKENDS.items():
            for content, action, message in cases:
                with self.subTest(mode=mode, action=action, message=message):
                    model.write_text(content)
                    result, output = backend.invoke(model, [cluster()], action)
                    self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                    self.assertIn(message, result.stdout + result.stderr)
                    self.assertEqual(output, '')

    def test_openmp_consistency(self):
        if not any(b.openmp for b in BACKENDS.values()):
            self.skipTest('OpenMP disabled with --no-openmp')
        for mode, backend in BACKENDS.items():
            for name in ('ordinary', 'typewise', 'charge1', 'charge2', 'zbl_flexible'):
                with self.subTest(mode=mode, model=name):
                    one = backend.run(MODELS[name], [cluster(), small_cell()], threads=1)
                    many = backend.run(MODELS[name], [cluster(), small_cell()], threads=3)
                    for a, b in zip(one, many):
                        self.compare_result(a, b, '1 vs 3 OpenMP threads', 3e-9, 3e-9)
