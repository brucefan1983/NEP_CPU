#!/usr/bin/env python3
"""Add or subtract NEP_CPU's D3(BJ) correction in an extended XYZ dataset."""
import argparse
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ELEMENTS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
            'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La '
            'Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi '
            'Po At Rn Fr Ra Ac Th Pa U Np Pu').split()
ATOMIC_NUMBER = {symbol: i + 1 for i, symbol in enumerate(ELEMENTS)}
VALUE = r'''"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\{[^}]*\}|\[[^\]]*\]|[^\s"'{}\[\]]+'''
FIELD = re.compile(r'([^\s=]+)(?:\s*=\s*(' + VALUE + r'))?(?=\s|$)')
TOKEN = re.compile(VALUE)


def unquote(value):
    return value[1:-1] if value and value[0] in ('"', "'", '[', '{') else value


def numbers(text, count):
    if text is None:
        raise ValueError('missing numerical value')
    result = [float(x) for x in text.split()]
    if len(result) != count or not all(math.isfinite(x) for x in result):
        raise ValueError('expected {} finite numbers'.format(count))
    return result


def replace_spans(text, edits):
    for start, end, value in sorted(edits, reverse=True):
        text = text[:start] + value + text[end:]
    return text


def parse_header(line):
    fields, position = {}, 0
    while position < len(line):
        if line[position].isspace():
            position += 1
            continue
        match = FIELD.match(line, position)
        if match is None:
            raise ValueError('malformed extended XYZ header near ' + repr(line[position:position + 40]))
        key = match[1].lower()
        if key in fields:
            raise ValueError('duplicate header field: ' + key)
        fields[key] = (unquote(match[2]), match.start(2), match.end(2))
        position = match.end()
    return fields


class Frame:
    def __init__(self, count_line, header, rows):
        self.count_line, self.header, self.rows = count_line, header, rows
        self.fields = parse_header(header)
        fields = self.fields
        for key in ('lattice', 'properties', 'energy'):
            if key not in fields or fields[key][0] is None:
                raise ValueError('missing ' + key)
        self.energy = numbers(fields['energy'][0], 1)[0]
        lattice = numbers(fields['lattice'][0], 9)
        self.box = [lattice[i] for i in (0, 3, 6, 1, 4, 7, 2, 5, 8)]
        a, b, c, d, e, f, g, h, i = self.box
        self.volume = abs(a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g))
        if not math.isfinite(self.volume) or self.volume <= 0:
            raise ValueError('Lattice must have nonzero volume')
        if 'pbc' in fields and (fields['pbc'][0] or '').lower().split() not in (
                ['t', 't', 't'], ['true', 'true', 'true'], ['1', '1', '1']):
            raise ValueError('the D3 core requires periodic boundary conditions in all three directions')
        self.tensors = {key: numbers(fields[key][0], 9) for key in ('virial', 'stress') if key in fields}
        properties = fields['properties'][0].split(':')
        if len(properties) % 3:
            raise ValueError('Properties must contain name:type:count triples')
        columns, offset = {}, 0
        for j in range(0, len(properties), 3):
            name, kind, width = properties[j:j + 3]
            name, width = name.lower(), int(width)
            if name in columns or width <= 0 or kind not in ('S', 'R', 'I', 'L'):
                raise ValueError('invalid Properties entry: ' + name)
            columns[name] = (offset, kind, width)
            offset += width
        forces = [name for name in ('force', 'forces') if name in columns]
        if len(forces) != 1:
            raise ValueError('Properties must contain exactly one force or forces field')
        for name, kind, width in [('species', 'S', 1), ('pos', 'R', 3), (forces[0], 'R', 3)]:
            if name not in columns or columns[name][1:] != (kind, width):
                raise ValueError('Properties requires {}:{}:{}'.format(name, kind, width))
        self.force_offset = columns[forces[0]][0]
        self.tokens, self.atoms, self.forces = [], [], []
        for row in rows:
            tokens = list(TOKEN.finditer(row))
            gaps = [row[a:b] for a, b in zip([0] + [t.end() for t in tokens],
                                             [t.start() for t in tokens] + [len(row)])]
            if len(tokens) != offset or any(gap.strip() for gap in gaps):
                raise ValueError('atom row does not match Properties')
            values = [unquote(token[0]) for token in tokens]
            symbol = values[columns['species'][0]]
            if symbol not in ATOMIC_NUMBER:
                raise ValueError('unsupported element: ' + symbol)
            start = columns['pos'][0]
            position = numbers(' '.join(values[start:start + 3]), 3)
            self.atoms.append((ATOMIC_NUMBER[symbol], position))
            start = self.force_offset
            self.forces.append(numbers(' '.join(values[start:start + 3]), 3))
            self.tokens.append(tokens)

    def corrected(self, energy, virial, forces, sign):
        changes = {'energy': format(self.energy + sign * energy, '.17g')}
        # Extended XYZ stores matrices column by column; the API uses row order.
        virial = [virial[i] for i in (0, 3, 6, 1, 4, 7, 2, 5, 8)]
        for key, original in self.tensors.items():
            factor = sign if key == 'virial' else -sign / self.volume
            changes[key] = '"' + ' '.join(format(x + factor * v, '.17g')
                                         for x, v in zip(original, virial)) + '"'
        edits = [(self.fields[key][1], self.fields[key][2], value) for key, value in changes.items()]
        yield self.count_line
        yield replace_spans(self.header, edits)
        for row, tokens, original, correction in zip(self.rows, self.tokens, self.forces, forces):
            edits = [(tokens[self.force_offset + d].start(), tokens[self.force_offset + d].end(),
                      format(original[d] + sign * correction[d], '.17g')) for d in range(3)]
            yield replace_spans(row, edits).rstrip('\r\n') + '\n'


def read_frames(stream):
    index = 0
    while True:
        count_line = stream.readline()
        if not count_line:
            return
        index += 1
        try:
            count = int(count_line)
            if count <= 0:
                raise ValueError('atom count must be positive')
            header = stream.readline()
            rows = [stream.readline() for _ in range(count)]
            if not header or not all(rows):
                raise ValueError('truncated frame')
            yield Frame(count_line, header, rows)
        except (ValueError, TypeError) as exc:
            raise ValueError('frame {}: {}'.format(index, exc)) from exc


def compile_backend(directory, cxx, flags='', openmp=True):
    executable = Path(directory) / 'dftd3'
    command = shlex.split(cxx) + ['-std=c++11', '-O3'] + shlex.split(flags)
    if openmp:
        command += ['-fopenmp']
    command += ['-I', str(ROOT / 'src'), str(HERE / 'dftd3.cpp')]
    command += [str(ROOT / 'src' / name)
                for name in ('nep.cpp', 'nep_model_file.cpp', 'ewald_nep.cpp', 'neighbor_nep.cpp')]
    command += ['-o', str(executable)]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError('compilation failed:\n' + result.stdout + result.stderr)
    return executable


def convert(source, target, executable, mode, functional, cutoff, cn_cutoff, threads=1):
    if not (math.isfinite(cutoff) and math.isfinite(cn_cutoff) and 0 < cn_cutoff <= cutoff):
        raise ValueError('cutoffs must satisfy 0 < cn-cutoff <= cutoff')
    if mode not in ('add', 'subtract') or threads < 1:
        raise ValueError('invalid mode or thread count')
    source, target = Path(source), Path(target)
    if source.resolve() == target.resolve() or (target.exists() and source.samefile(target)):
        raise ValueError('input and output must be different files')
    if target.exists():
        raise ValueError('output file already exists: ' + str(target))
    sign = 1 if mode == 'add' else -1
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE')
    temporary = None
    try:
        with source.open() as input_file, tempfile.TemporaryFile(mode='w+') as errors:
            command = [str(executable), functional, repr(cutoff), repr(cn_cutoff)]
            with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=errors, text=True, env=env) as process:
                try:
                    with tempfile.NamedTemporaryFile(mode='w', dir=target.parent, delete=False) as output:
                        temporary = Path(output.name)
                        count = 0
                        for count, frame in enumerate(read_frames(input_file), 1):
                            process.stdin.write(str(len(frame.atoms)) + '\n')
                            process.stdin.write(' '.join(format(x, '.17g') for x in frame.box) + '\n')
                            for z, position in frame.atoms:
                                process.stdin.write(str(z) + ' ' + ' '.join(format(x, '.17g') for x in position) + '\n')
                            process.stdin.flush()
                            values = numbers(process.stdout.readline(), 10)
                            forces = [numbers(process.stdout.readline(), 3) for _ in frame.atoms]
                            output.writelines(frame.corrected(values[0], values[1:], forces, sign))
                        if not count:
                            raise ValueError('input dataset is empty')
                    process.stdin.close()
                    if process.wait() != 0:
                        raise RuntimeError('D3 calculation failed')
                except BaseException as exc:
                    if process.poll() is None:
                        process.terminate()
                    process.wait()
                    errors.seek(0)
                    detail = errors.read().strip()
                    if detail:
                        raise RuntimeError(detail) from exc
                    raise
        # Publish only a complete dataset, without overwriting an existing file.
        os.link(temporary, target)
        return count
    finally:
        if temporary is not None:
            temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('add', 'subtract'))
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--functional', required=True, help='D3(BJ) functional, e.g. pbe')
    parser.add_argument('--cutoff', required=True, type=float, help='D3 energy cutoff in Angstrom')
    parser.add_argument('--cn-cutoff', required=True, type=float, help='coordination-number cutoff in Angstrom')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'g++'))
    parser.add_argument('--cxxflags', default=os.environ.get('CXXFLAGS', ''))
    parser.add_argument('--no-openmp', action='store_true')
    args = parser.parse_args()
    try:
        with tempfile.TemporaryDirectory(prefix='nep_dftd3_') as directory:
            executable = compile_backend(directory, args.cxx, args.cxxflags, not args.no_openmp)
            count = convert(args.input, args.output, executable, args.mode, args.functional,
                            args.cutoff, args.cn_cutoff, args.threads)
        print('Wrote {} frames to {}'.format(count, args.output))
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
