"""Numerical regressions and independent consistency checks for the public API."""
import copy
import hashlib
import itertools
import json
import math
import unittest

from fixtures import (ALL_FLAGS, FLEXIBLE_ZBL, ZBL_COEFFICIENTS, cluster,
                      model_text, pair, raw_charges, small_cell, supercell, write_models)
from support import (HERE, ROOT, displaced, flatten, matvec, read_matrix, read_xyz,
                     strained, tensor_transform, transform)


def configure(directory, backends):
    global WORK, BACKENDS, MODELS, REFERENCE
    WORK, BACKENDS = directory, backends
    MODELS = write_models(directory)
    REFERENCE = json.loads((HERE / 'references.json').read_text())


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
        for name, step in [('nep', 2e-5), ('qnep', 2e-4)]:
            with self.subTest(fixture=name):
                directory = HERE / 'data' / name
                frame = read_xyz(directory)
                frames = [displaced(frame, i, d, sign * step)
                          for i in range(len(frame['type'])) for d in range(3) for sign in (-1, 1)]
                energies = [r['total_energy'] for r in backend.run(directory / 'nep.txt', frames, 'energy')]
                fd = [(energies[k] - energies[k + 1]) / (2 * step) for k in range(0, len(energies), 2)]
                fd = [fd[i:i + 3] for i in range(0, len(fd), 3)]
                self.close(fd, read_matrix(directory / 'force_finite_difference_ref.out'),
                           name + '/saved FD', 2e-7, 2e-7)
                # qNEP's historical 2e-4 step crosses a cutoff for some pairs.
                # Preserve that regression, but test derivative consistency at
                # 2e-6 separately instead of weakening the force assertion.
                # The closest pair is only 8.85e-6 Angstrom below the cutoff.
                if name == 'qnep':
                    step = 2e-6
                    frames = [displaced(frame, i, d, sign * step)
                              for i in range(len(frame['type'])) for d in range(3) for sign in (-1, 1)]
                    energies = [r['total_energy'] for r in backend.run(directory / 'nep.txt', frames, 'energy')]
                    fd = [(energies[k] - energies[k + 1]) / (2 * step) for k in range(0, len(energies), 2)]
                    fd = [fd[i:i + 3] for i in range(0, len(fd), 3)]
                actual = backend.run(directory / 'nep.txt', [frame])[0]['force']
                self.close(fd, actual, name + '/force = -dE/dr', 2e-7, 2e-7)

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
        """64 layouts, each with energy/force/virial checks for NEP and qNEP1/2."""
        if 'analytic' not in BACKENDS:
            self.skipTest('Finite differences target the analytic implementation')
        backend = BACKENDS['analytic']
        for flags in itertools.product((0, 1), repeat=6):
            for kind in ('nep4', 'nep4_charge1', 'nep4_charge2'):
                with self.subTest(flags=flags, kind=kind):
                    model = WORK / 'flags.txt'
                    model.write_text(model_text(kind=kind, flags=flags))
                    actual = backend.run(model, [cluster()], 'load')[0]
                    self.assertEqual(actual['dimension'], 3 + 3 * (4 + sum(flags)))
                    self.check_derivatives(backend, model, cluster(), all_coordinates=False)

    def test_charge_and_bec(self):
        frame = cluster()
        for mode, backend in BACKENDS.items():
            for name in ('charge1', 'charge2'):
                with self.subTest(mode=mode, model=name):
                    model = MODELS[name]
                    base = backend.run(model, [frame])[0]
                    q = [row[0] for row in base['charge']]
                    self.close(sum(q), 0.0, 'charge neutrality', 2e-12, 0)
                    self.close([sum(row[d] for row in base['bec']) for d in range(9)],
                               [0.0] * 9, 'BEC acoustic sum', 2e-11, 0)
                    raw = raw_charges(model, base['descriptor'], frame['type'])
                    mean = sum(raw) / len(raw)
                    self.close(q, [x - mean for x in raw], 'ANN charge / mean subtraction')
                    if mode != 'analytic':
                        continue
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
                model.write_text(model_text(kind='nep4_zbl', zbl=zbl, zero=True))
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
        model = ROOT / 'test_dftd3/C_2022_NEP4.txt'
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
                common_bias = float(lines[6 + 2 * neurons * (dim + 2)])
                self.close([[sum(row) - common_bias] for row in base['latent']], base['energy'], 'latent reconstruction')
                # Perturb input weight, hidden bias and output weight for both
                # species. B_projection is per-atom, neuron-major within atom.
                for t, j, component in itertools.product((0, 1), (0, 7), (0, dim, dim + 1)):
                    with self.subTest(type=t, neuron=j, parameter=component):
                        offset = j * dim if component == 0 else neurons * dim + j
                        if component == dim + 1:
                            offset += neurons
                        index = 6 + t * neurons * (dim + 2) + offset
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
                        'nep5', 'nep5_zbl', 'nep4_charge3', 'nep4_zbl_charge3', 'unknown'):
                path = WORK / (tag + '.txt')
                path.write_text(tag + ' 2 C Si\n')
                cases.append((path, 'load', 'unsupported NEP model'))
            for name, text, expected in [
                ('two_layers', model_text().replace('ANN 8 0', 'ANN 8 4'), 'one hidden layer'),
                ('short_ann', model_text().replace('ANN 8 0', 'ANN 8'), 'ANN num_neurons 0'),
                ('bad_count', 'nep4 3 C Si\n', 'atom symbols'),
                ('short_header', 'nep4\n', 'at least 3 items')]:
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

    def test_q222_legacy_marker(self):
        for mode, backend in BACKENDS.items():
            with self.subTest(mode=mode):
                old = WORK / 'q222_marker_1.txt'
                old.write_text(model_text().replace('l_max 4 2', 'l_max 4 1'))
                a, b = [backend.run(p, [cluster()])[0] for p in (MODELS['ordinary'], old)]
                self.compare_result(a, b, 'q222 markers 1 and 2')

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
