#!/usr/bin/env python3
"""Measure standalone NEP_CPU API throughput with reproducible workloads."""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import shlex
import statistics
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FRAME_COUNT = 8


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_xyz(directory):
    symbols = (directory / 'nep.txt').read_text().splitlines()[0].split()[2:]
    lines = (directory / 'xyz.in').read_text().splitlines()
    cell = list(map(float, lines[1].split()))
    rows = [line.split() for line in lines[2:]]
    if len(rows) != int(lines[0]) or len(cell) != 9:
        raise ValueError('Invalid structure: ' + str(directory))
    return dict(box=[[cell[3 * j + i] for j in range(3)] for i in range(3)],
                type=[symbols.index(row[0]) for row in rows],
                position=[list(map(float, row[1:4])) for row in rows])


def replicate(frame, count):
    result = dict(box=[[count * x for x in row] for row in frame['box']], type=[], position=[])
    for image in itertools.product(range(count), repeat=3):
        shift = [sum(x * k for x, k in zip(row, image)) for row in frame['box']]
        result['type'].extend(frame['type'])
        result['position'].extend([[x + y for x, y in zip(r, shift)] for r in frame['position']])
    return result


def diamond(count, two_types=False):
    basis = [(0., 0., 0.), (0., .5, .5), (.5, 0., .5), (.5, .5, 0.)]
    basis += [tuple(x + .25 for x in r) for r in basis]
    cell = dict(box=[[3.567 if i == j else 0. for j in range(3)] for i in range(3)],
                type=[i % 2 if two_types else 0 for i in range(8)],
                position=[[3.567 * x for x in r] for r in basis])
    return replicate(cell, count)


def synthetic_model(l_max):
    """Two-species kernel workload, not a fitted physical potential."""
    dim, neurons = 5 + 5 * (l_max + 1), 30
    ann_count = 2 * neurons * (dim + 2) + 1
    count = ann_count + 2 * 4 * 5 * 9
    rng = random.Random(39017)
    lines = ['nep 2 C Si', 'format_version 1', 'model nep', 'cutoff 5 4 200 100',
             'n_max 4 4', 'basis_size 8 8', 'l_max {} 2 0'.format(l_max), 'ann 30 0',
             'energy_reference 0 0', '@ parameters ' + str(count)]
    lines += [format(rng.uniform(-.15, .15), '.17g') for _ in range(count)]
    lines += ['@ q_scaler ' + str(dim)] + ['1'] * dim
    return '\n'.join(lines) + '\n'


def write_frames(path, frame):
    with path.open('w') as out:
        out.write('{} {}\n'.format(len(frame['type']), FRAME_COUNT))
        out.write(' '.join(format(x, '.17g') for row in frame['box'] for x in row) + '\n')
        out.write(' '.join(map(str, frame['type'])) + '\n')
        for step in range(FRAME_COUNT):
            phase = 2 * math.pi * step / FRAME_COUNT
            for i, r in enumerate(frame['position']):
                out.write(' '.join(format(x + .01 * math.sin(phase + .71 * i + 1.13 * d), '.17g')
                                   for d, x in enumerate(r)) + '\n')


def workloads(directory, only):
    data = ROOT / 'tests/api/data'
    carbon = ROOT / 'tests/lammps/diamond/nep.txt'
    pbte = read_xyz(data / 'nep')
    cases = [('diamond_1728', carbon, diamond(6)),
             ('diamond_13824', carbon, diamond(12)),
             ('pbte_250', data / 'nep/nep.txt', pbte),
             ('pbte_2000', data / 'nep/nep.txt', replicate(pbte, 2)),
             ('qnep_250', data / 'qnep/nep.txt', read_xyz(data / 'qnep'))]
    for l_max in (6, 8):
        model = directory / 'synthetic_l{}.txt'.format(l_max)
        model.write_text(synthetic_model(l_max))
        cases.append(('synthetic_l{}_1728'.format(l_max), model, diamond(6, two_types=True)))
    result = []
    for name, model, frame in cases:
        if only and only not in name:
            continue
        source = directory / (name + '.frames')
        write_frames(source, frame)
        result.append(dict(case=name, atoms=len(frame['type']), model=model, source=source,
                           model_sha256=digest(model), frames_sha256=digest(source)))
    if not result:
        raise ValueError('No case matches --only ' + only)
    return result


def command_output(command, env=None, timeout=300):
    result = subprocess.run(command, text=True, capture_output=True, env=env, timeout=timeout)
    if result.returncode:
        raise RuntimeError('Command failed: {}\n{}{}'.format(
            shlex.join(map(str, command)), result.stdout, result.stderr))
    return result.stdout.strip()


def machine_info():
    cpu = platform.processor()
    cpuinfo = Path('/proc/cpuinfo')
    if cpuinfo.exists():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith('model name'):
                cpu = line.split(':', 1)[1].strip()
                break
    return dict(cpu=cpu, architecture=platform.machine(), system=platform.platform(),
                logical_cpus=os.cpu_count(),
                affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None)


def git_info():
    try:
        revision = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True)
        status = subprocess.run(['git', '-C', str(ROOT), 'status', '--porcelain'],
                                capture_output=True, text=True)
        return dict(commit=revision.stdout.strip() if revision.returncode == 0 else None,
                    dirty=bool(status.stdout) if status.returncode == 0 else None)
    except FileNotFoundError:
        return dict(commit=None, dirty=None)


def same_signatures(actual, expected):
    return (len(actual) == len(expected) and
            all(len(a) == len(b) and all(math.isfinite(x) and math.isfinite(y) and
                abs(x - y) <= 1e-8 * (1 + abs(y)) for x, y in zip(a, b))
                for a, b in zip(actual, expected)))


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def positive_float(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='new JSON result file (default: benchmarks/results/)')
    parser.add_argument('--compare', type=Path, help='compare against an earlier JSON result')
    parser.add_argument('--threads', type=positive_int, nargs='+', default=[1], help='thread counts (default: 1)')
    parser.add_argument('--mode', choices=('analytic', 'table'), default='analytic')
    parser.add_argument('--only', default='', help='run case names containing this text')
    parser.add_argument('--seconds', type=positive_float, default=1., help='target seconds per batch (default: 1)')
    parser.add_argument('--repeats', type=positive_int, default=5, help='timed batches per case (default: 5)')
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'g++'))
    parser.add_argument('--cxxflags', default=os.environ.get('CXXFLAGS', ''), help='additional compiler flags')
    args = parser.parse_args()
    try:
        timestamp = datetime.now(timezone.utc)
        output = args.output or HERE / 'results' / timestamp.strftime('nep_%Y%m%dT%H%M%S_%fZ.json')
        output = output.resolve()
        if output.exists():
            raise ValueError('Output already exists: ' + str(output))
        baseline = json.loads(args.compare.read_text()) if args.compare else None
        if baseline is not None and baseline.get('schema_version') != 1:
            raise ValueError('Unsupported baseline format')
        compiler = shlex.split(args.cxx)
        if not compiler:
            raise ValueError('Empty compiler command')
        flags = ['-std=c++11', '-O3', '-fopenmp'] + shlex.split(args.cxxflags)
        if args.mode == 'table':
            flags += ['-DUSE_TABLE_FOR_RADIAL_FUNCTIONS']
        env = dict(os.environ, OMP_DYNAMIC='FALSE', OMP_MAX_ACTIVE_LEVELS='1')
        env.setdefault('OMP_PROC_BIND', 'close')
        env.setdefault('OMP_PLACES', 'cores')
        omp_keys = ('OMP_DYNAMIC', 'OMP_MAX_ACTIVE_LEVELS', 'OMP_PROC_BIND', 'OMP_PLACES',
                    'OMP_WAIT_POLICY', 'OMP_THREAD_LIMIT', 'GOMP_CPU_AFFINITY', 'KMP_AFFINITY')
        report = dict(schema_version=1, created_utc=timestamp.isoformat(), mode=args.mode,
                      machine=machine_info(), git=git_info(),
                      build=dict(compiler=compiler, version=command_output(compiler + ['--version']).splitlines()[0],
                                 flags=flags, driver_sha256=digest(HERE / 'driver.cpp')),
                      omp_environment={key: env[key] for key in omp_keys if key in env},
                      source_sha256={p.relative_to(ROOT).as_posix(): digest(p)
                                     for p in sorted((ROOT / 'src').iterdir()) if p.is_file()},
                      script_sha256=digest(Path(__file__).resolve()),
                      target_seconds=args.seconds, repeats=args.repeats, results=[])
        previous = {}
        if baseline is not None:
            report['compare_to'] = dict(file=args.compare.name, sha256=digest(args.compare),
                                        created_utc=baseline['created_utc'], git=baseline['git'])
            for key in ('mode', 'machine', 'build', 'omp_environment'):
                if report[key] != baseline[key]:
                    raise ValueError('Cannot compare: {} differs from the baseline'.format(key))
            previous = {(row['case'], row['threads']): row for row in baseline['results']}
        print('CPU: ' + report['machine']['cpu'], flush=True)
        print('Compiler: ' + report['build']['version'], flush=True)
        print('Build: {}; mode: {}'.format(shlex.join(flags), args.mode), flush=True)
        problems = []
        with tempfile.TemporaryDirectory(prefix='nep_benchmark_') as tmp:
            directory = Path(tmp)
            cases = workloads(directory, args.only)
            threads = list(dict.fromkeys(args.threads))
            if baseline is not None:
                for case in cases:
                    for count in threads:
                        old = previous.get((case['case'], count))
                        if old is None:
                            raise ValueError('Baseline is missing {} / {} threads'.format(case['case'], count))
                        for key in ('model_sha256', 'frames_sha256'):
                            if case[key] != old[key]:
                                raise ValueError('Cannot compare {}: {} differs'.format(case['case'], key))
                        if not math.isfinite(old['median_ms']) or old['median_ms'] <= 0:
                            raise ValueError('Invalid baseline timing')
            executable = directory / 'driver'
            sources = [HERE / 'driver.cpp'] + sorted((ROOT / 'src').glob('*.cpp'))
            print('Compiling benchmark driver...', flush=True)
            command_output(compiler + flags + ['-I', str(ROOT / 'src')] + list(map(str, sources)) +
                           ['-o', str(executable)])
            for case in cases:
                for count in threads:
                    print('Running {}: {} atoms, {} thread(s)...'.format(case['case'], case['atoms'], count), flush=True)
                    target = directory / 'timing.json'
                    command_output([str(executable), str(case['model']), str(case['source']), str(target),
                                    str(args.repeats), str(args.seconds)],
                                   env=dict(env, OMP_NUM_THREADS=str(count)), timeout=3600)
                    row = json.loads(target.read_text())
                    if row['threads'] != count or row['atoms'] != case['atoms'] or row['frames'] != FRAME_COUNT:
                        raise RuntimeError('Driver workload or thread count differs from the request')
                    elapsed = row['elapsed_seconds']
                    if len(elapsed) != args.repeats or not all(math.isfinite(t) and t > 0 for t in elapsed):
                        raise RuntimeError('Invalid driver timings')
                    timings = [1000 * t / row['steps_per_repeat'] for t in elapsed]
                    median = statistics.median(timings)
                    row.update({key: case[key] for key in ('case', 'model_sha256', 'frames_sha256')})
                    row.update(median_ms=median, min_ms=min(timings), max_ms=max(timings),
                               spread_percent=100 * (max(timings) - min(timings)) / median,
                               katom_steps_per_second=case['atoms'] / median)
                    if baseline is not None:
                        old = previous[(case['case'], count)]
                        row['baseline_ms'] = old['median_ms']
                        row['signatures_match'] = same_signatures(row['signatures'], old['signatures'])
                        if row['signatures_match']:
                            row['speedup'] = old['median_ms'] / median
                        else:
                            problems.append('{} / {} threads: numerical signatures differ'.format(case['case'], count))
                    report['results'].append(row)
        print('\n{:<20} {:>7} {:>7} {:>12} {:>12} {:>9}{}'.format(
            'Case', 'Atoms', 'Threads', 'ms/call', 'katom-step/s', 'Spread %',
            '  Baseline ms   Speedup' if baseline else ''))
        for row in report['results']:
            comparison = ''
            if baseline is not None:
                comparison = ' {:>12.3f}   '.format(row['baseline_ms'])
                comparison += '{:.3f}x'.format(row['speedup']) if row['signatures_match'] else 'DIFF'
            print('{case:<20} {atoms:>7} {threads:>7} {median_ms:>12.3f} '
                  '{katom_steps_per_second:>12.2f} {spread_percent:>9.1f}'.format(**row) + comparison)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x') as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print('\nSaved: ' + str(output))
        for problem in problems:
            print('ERROR: ' + problem, file=sys.stderr)
        return 1 if problems else 0
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
