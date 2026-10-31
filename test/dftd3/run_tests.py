#!/usr/bin/env python3
"""Regression tests for the D3 dataset tool; uses only the Python standard library."""
import io
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'dftd3'))
sys.path.insert(0, str(ROOT / 'test' / 'api'))
import dftd3
from fixtures import cluster


def dataset(mixed=False):
    frame = cluster()
    lattice = [frame['box'][i][j] for j in range(3) for i in range(3)]
    header = ('Lattice="{}" energy = -12.3 virial="{}" stress="{}" weight=0.25 '
              'config_type="Mixed CASE" pbc="T T T" note="keep \\"this\\"" flag '
              'Properties=id:I:1:forces:R:3:species:S:1:pos:R:3:label:S:1\n').format(
                  ' '.join(map(str, lattice)), ' '.join(str(i * .1) for i in range(9)),
                  ' '.join(str(i * -.01) for i in range(9)))
    rows = []
    for i, (t, r) in enumerate(zip(frame['type'], frame['position'])):
        rows.append('{}  0.1 -0.2 0.3\t{} {}  "atom {}"\n'.format(
            i, ('C', 'Si')[t] if mixed else 'C', ' '.join(map(str, r)), i))
    return str(len(rows)) + '\n' + header + ''.join(rows)


class D3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build = tempfile.TemporaryDirectory(prefix='nep_d3_tests_')
        cls.addClassCleanup(cls.build.cleanup)
        print('Compiling D3 helper...', flush=True)
        cls.executable = dftd3.compile_backend(cls.build.name, os.environ.get('CXX', 'g++'),
                                             os.environ.get('CXXFLAGS', ''))
        cls.reference = json.loads((ROOT / 'test/api/references.json').read_text())['analytic']['d3']

    def setUp(self):
        self.work = tempfile.TemporaryDirectory(dir=self.build.name)
        self.addCleanup(self.work.cleanup)
        self.directory = Path(self.work.name)

    def run_conversion(self, text, mode='add', threads=1, functional='pbe'):
        source, target = self.directory / 'input.xyz', self.directory / 'output.xyz'
        source.write_text(text)
        if target.exists():
            target.unlink()
        count = dftd3.convert(source, target, self.executable, mode, functional, 8.0, 6.0, threads)
        self.assertEqual(source.read_text(), text)
        result = target.read_text()
        frames = list(dftd3.read_frames(io.StringIO(result)))
        self.assertEqual(count, len(frames))
        return result, frames

    def close(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for x, y in zip(actual, expected):
            self.assertTrue(math.isfinite(x) and math.isfinite(y))
            self.assertAlmostEqual(x, y, delta=2e-10 * (1 + abs(y)))

    def unchanged_fields(self, before, after):
        ignored = ('energy', 'virial', 'stress')
        def header_without_results(frame):
            return dftd3.replace_spans(frame.header, [(v[1], v[2], '_')
                                       for key, v in frame.fields.items() if key in ignored])
        self.assertEqual(header_without_results(before), header_without_results(after))
        self.assertEqual(before.atoms, after.atoms)
        for old, new, old_tokens, new_tokens in zip(before.rows, after.rows, before.tokens, after.tokens):
            def without_forces(row, tokens):
                return dftd3.replace_spans(row, [(tokens[j].start(), tokens[j].end(), '_')
                    for j in range(before.force_offset, before.force_offset + 3)])
            self.assertEqual(without_forces(old, old_tokens), without_forces(new, new_tokens))

    def check_reference(self, text, correction):
        before = next(dftd3.read_frames(io.StringIO(text)))
        total = [sum(row[d] for row in correction['virial']) for d in (0, 3, 6, 1, 4, 7, 2, 5, 8)]
        for mode, sign in [('add', 1), ('subtract', -1)]:
            with self.subTest(mode=mode):
                _, frames = self.run_conversion(text, mode)
                after = frames[0]
                self.close([after.energy], [before.energy + sign * correction['total_energy']])
                for old, new, force in zip(before.forces, after.forces, correction['force']):
                    self.close(new, [x + sign * y for x, y in zip(old, force)])
                self.close(after.tensors['virial'], [x + sign * y for x, y in zip(before.tensors['virial'], total)])
                self.close(after.tensors['stress'], [x - sign * y / before.volume
                                                    for x, y in zip(before.tensors['stress'], total)])
                self.unchanged_fields(before, after)

    def test_saved_d3_reference(self):
        self.check_reference(dataset(), self.reference)

    def test_mixed_species_reference(self):
        reference = json.loads((Path(__file__).parent / 'reference.json').read_text())['result']
        self.check_reference(dataset(mixed=True), reference)

    def test_roundtrip_and_changing_atom_count(self):
        example = (ROOT / 'tools/dftd3/example.xyz').read_text()
        text = dataset(mixed=True) + example + dataset()
        added, _ = self.run_conversion(text)
        _, frames = self.run_conversion(added, 'subtract')
        originals = list(dftd3.read_frames(io.StringIO(text)))
        self.assertEqual(len(frames), len(originals))
        for old, new in zip(originals, frames):
            self.unchanged_fields(old, new)
            self.close([old.energy], [new.energy])
            for a, b in zip(old.forces, new.forces):
                self.close(a, b)
            for name in old.tensors:
                self.close(old.tensors[name], new.tensors[name])

    def test_optional_tensors_and_singular_force(self):
        text = dataset()
        before = next(dftd3.read_frames(io.StringIO(text)))
        for tensor in ('virial', 'stress'):
            text = text.replace(tensor + '="' + before.fields[tensor][0] + '" ', '')
        text = text.replace('forces:R:3', 'force:R:3').replace('pbc="T T T" ', '')
        _, frames = self.run_conversion(text)
        self.assertEqual(frames[0].tensors, {})
        self.assertNotIn('pbc', frames[0].fields)
        self.close([frames[0].energy], [before.energy + self.reference['total_energy']])

    def test_openmp(self):
        one, _ = self.run_conversion(dataset(), threads=1)
        three, _ = self.run_conversion(dataset(), threads=3)
        self.assertEqual(one, three)

    def test_invalid_inputs_leave_no_output(self):
        text = dataset()
        cases = ['', text.replace('energy = -12.3', 'other_energy=-12.3'),
                 text.replace('energy = -12.3', 'energy=nan'),
                 text.replace('pbc="T T T"', 'pbc="T T F"'),
                 text.replace('pbc="T T T"', 'pbc'),
                 text.replace('species:S:1', 'other:S:1'),
                 text.replace('forces:R:3', 'other:R:3'),
                 text.replace('\tC ', '\tXx '), text.rsplit('\n', 2)[0] + '\n',
                 text.replace('energy = -12.3', 'energy=-12.3 energy=2'),
                 text.replace('Lattice="20.0', 'Lattice="0.0'),
                 text.replace('forces:R:3', 'forces:R:2')]
        for i, content in enumerate(cases):
            with self.subTest(case=i):
                with self.assertRaises(ValueError):
                    self.run_conversion(content)
                self.assertFalse((self.directory / 'output.xyz').exists())
        # A failure after a successfully processed frame must also discard the output.
        with self.assertRaises(ValueError):
            self.run_conversion(text + '3\ntruncated\n')
        self.assertFalse((self.directory / 'output.xyz').exists())

    def test_existing_output_and_invalid_settings(self):
        source, target = self.directory / 'input.xyz', self.directory / 'output.xyz'
        source.write_text(dataset())
        target.write_text('existing output\n')
        for output in (source, target):
            with self.assertRaises(ValueError):
                dftd3.convert(source, output, self.executable, 'add', 'pbe', 8., 6.)
        self.assertEqual(target.read_text(), 'existing output\n')
        target.unlink()
        for cutoffs in [(0., 6.), (8., -1.), (4., 6.), (math.nan, 6.)]:
            with self.assertRaises(ValueError):
                dftd3.convert(source, target, self.executable, 'add', 'pbe', *cutoffs)
        self.assertFalse(target.exists())

    def test_unsupported_functional(self):
        with self.assertRaisesRegex(RuntimeError, 'not supported'):
            self.run_conversion(dataset(), functional='unknown_functional')
        self.assertFalse((self.directory / 'output.xyz').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
