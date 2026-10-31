"""Standard-library-only helpers; test decisions and tolerances live in test_api.py."""
import copy
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
IDENTITY = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def flatten(value):
    if isinstance(value, list):
        for item in value:
            yield from flatten(item)
    else:
        yield value


def matvec(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def tensor_transform(flat, rotation):
    a = [flat[3 * d:3 * d + 3] for d in range(3)]
    return list(flatten(matmul(matmul(rotation, a), transpose(rotation))))


def transform(frame, matrix):
    result = copy.deepcopy(frame)
    result['box'] = matmul(matrix, frame['box'])
    result['position'] = [matvec(matrix, r) for r in frame['position']]
    return result


def displaced(frame, atom, axis, delta):
    result = copy.deepcopy(frame)
    result['position'][atom][axis] += delta
    return result


def strained(frame, a, b, delta):
    matrix = copy.deepcopy(IDENTITY)
    matrix[a][b] += delta
    return transform(frame, matrix)


def read_xyz(directory):
    symbols = (directory / 'nep.txt').read_text().splitlines()[0].split()[2:]
    lines = (directory / 'xyz.in').read_text().splitlines()
    n = int(lines[0])
    cell = list(map(float, lines[1].split()))
    rows = [line.split() for line in lines[2:]]
    if len(rows) != n:
        raise ValueError('Unexpected atom count in ' + str(directory))
    return dict(box=transpose([cell[i:i + 3] for i in (0, 3, 6)]),
                type=[symbols.index(row[0]) for row in rows],
                position=[list(map(float, row[1:4])) for row in rows])


def read_matrix(path):
    return [list(map(float, row.split())) for row in path.read_text().splitlines() if row.strip()]


class Backend:
    def __init__(self, directory, mode, cxx, flags, openmp=False, threads=1):
        self.directory = Path(directory)
        self.mode = mode
        self.openmp = openmp
        self.threads = threads
        self.executable = self.directory / ('driver_' + mode)
        command = shlex.split(cxx) + ['-std=c++11', '-O2'] + shlex.split(flags)
        if openmp:
            command += ['-fopenmp']
        if mode == 'table':
            command += ['-DUSE_TABLE_FOR_RADIAL_FUNCTIONS']
        command += ['-I', str(ROOT / 'src'), str(HERE / 'driver.cpp')]
        command += [str(ROOT / 'src' / f)
                    for f in ('nep.cpp', 'nep_model_file.cpp', 'ewald_nep.cpp', 'neighbor_nep.cpp')]
        command += ['-o', str(self.executable)]
        print('Compiling ' + mode + ': ' + shlex.join(command), flush=True)
        result = subprocess.run(command, text=True, capture_output=True, timeout=180)
        if result.returncode:
            raise RuntimeError('Compilation failed:\n' + result.stdout + result.stderr)

    def invoke(self, model, frames, action='full', threads=None, rank=0):
        # A fresh process isolates exit(1) and memory corruption; frames within
        # that process share one initialized NEP instance (MD-like reuse).
        with tempfile.TemporaryDirectory(prefix='case_', dir=self.directory) as tmp:
            source, target = Path(tmp) / 'input.txt', Path(tmp) / 'output.jsonl'
            with source.open('w') as stream:
                for frame in frames:
                    stream.write('{} {}\n'.format(action, len(frame['type'])))
                    stream.write(' '.join(format(x, '.17g') for x in flatten(frame['box'])) + '\n')
                    for t, r in zip(frame['type'], frame['position']):
                        stream.write(str(t) + ' ' + ' '.join(format(x, '.17g') for x in r) + '\n')
            env = dict(os.environ, OMP_NUM_THREADS=str(threads or self.threads), OMP_DYNAMIC='FALSE')
            result = subprocess.run([str(self.executable), str(model), str(source), str(target), str(rank)],
                                    text=True, capture_output=True, env=env, timeout=120)
            text = target.read_text() if target.exists() else ''
            return result, text

    def run(self, model, frames, action='full', threads=None):
        # Bound process time and temporary file size for long frame sequences.
        # Short MD-like reuse sequences stay intact.
        if len(frames) > 128:
            values = []
            for start in range(0, len(frames), 128):
                values.extend(self.run(model, frames[start:start + 128], action, threads))
            return values
        result, text = self.invoke(model, frames, action, threads)
        if result.returncode:
            raise AssertionError('{} / {} / {} exited {}:\n{}{}'.format(
                self.mode, Path(model).name, action, result.returncode, result.stdout, result.stderr))
        try:
            values = [json.loads(line) for line in text.splitlines()]
        except ValueError as exc:
            raise AssertionError('Non-JSON output (possibly NaN/Inf): ' + text[:500]) from exc
        if len(values) != len(frames):
            raise AssertionError('Missing result frames: {} != {}'.format(len(values), len(frames)))
        for frame in values:
            for key, value in frame.items():
                if not all(math.isfinite(x) for x in flatten(value)):
                    raise AssertionError('Nonfinite {} in {}'.format(key, model))
        return values
