"""Instrument temporary benchmark sources; keep the normal timing build unchanged."""
import hashlib
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
STAGES = ('setup', 'neighbors', 'descriptor', 'radial', 'angular', 'bec', 'electrostatics', 'zbl')
LABELS = ('Setup', 'Neighbor', 'Desc/ANN', 'Radial', 'Angular', 'BEC', 'Coulomb', 'ZBL', 'Other')


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError('Profiling source layout changed; update benchmarks/stage_profile.py: ' + old.strip())
    return text.replace(old, new, 1)


def instrument_nep(text):
    starts = list(re.finditer(r'^void NEP::compute\(', text, re.MULTILINE))
    if len(starts) != 2:
        raise RuntimeError('Profiling requires the two NEP::compute overloads; update benchmarks/stage_profile.py')
    for match in reversed(starts):
        opening = text.index('\n{\n', match.start())
        closing = text.index('\n}\n', opening) + 2
        body = text[opening:closing]
        charged = 'std::vector<double>& charge' in text[match.start():opening]
        markers = [('  find_neighbor_list_small_box(', 'neighbors'),
                   ('  find_descriptor_small_box(', 'descriptor')]
        if charged:
            markers += [('  zero_total_charge(', 'bec'), ('  ewald.find_force(', 'electrostatics')]
        markers += [('  find_force_radial_small_box(', 'radial'),
                    ('  find_force_angular_small_box(', 'angular'),
                    ('  if (zbl.enabled) {', 'zbl')]
        positions = [body.find(marker) for marker, stage in markers]
        if any(pos < 0 for pos in positions) or positions != sorted(positions):
            raise RuntimeError('Profiling compute stages changed; update benchmarks/stage_profile.py')
        for marker, stage in markers:
            body = replace_once(body, marker, '  stage_timer.start(nep_benchmark::' + stage + ');\n' + marker)
        body = replace_once(body, '\n{\n', '\n{\n  nep_benchmark::Call stage_timer;\n')
        text = text[:opening] + body + text[closing:]
    return '#include "stage_timer.h"\n' + text


def instrument_driver(text):
    text = replace_once(text, '      elapsed.push_back(measure(steps));',
                        '      nep_benchmark::begin_batch();\n'
                        '      elapsed.push_back(measure(steps));\n'
                        '      nep_benchmark::end_batch();')
    text = replace_once(text, '    out << "]}\\n";',
                        '    out << \']\';\n'
                        '    nep_benchmark::write(out);\n'
                        '    out << "}\\n";')
    return '#include "stage_timer.h"\n' + text


def build(directory, root, compiler, flags, command_output):
    directory = directory / 'profile'
    directory.mkdir()
    nep = directory / 'nep.cpp'
    driver = directory / 'driver.cpp'
    nep.write_text(instrument_nep((root / 'src/nep.cpp').read_text()))
    driver.write_text(instrument_driver((HERE / 'driver.cpp').read_text()))
    executable = directory / 'driver'
    sources = [driver, nep] + [p for p in sorted((root / 'src').glob('*.cpp')) if p.name != 'nep.cpp']
    command_output(compiler + flags + ['-I', str(root / 'src'), '-I', str(HERE)] +
                   list(map(str, sources)) + ['-o', str(executable)])
    paths = {'instrumenter': Path(__file__), 'timer': HERE / 'stage_timer.h',
             'driver': driver, 'nep': nep}
    metadata = dict(stages=list(STAGES), sha256={name: hashlib.sha256(path.read_bytes()).hexdigest()
                                               for name, path in paths.items()})
    return executable, metadata


def summarize(raw, expected, repeats, same_signatures):
    for key in ('atoms', 'frames', 'threads'):
        if raw[key] != expected[key]:
            raise RuntimeError('Profile workload differs: ' + key)
    if not same_signatures(raw['signatures'], expected['signatures']):
        raise RuntimeError('Profile numerical signatures differ from the uninstrumented build')
    elapsed, stages, counts = raw['elapsed_seconds'], raw['stage_seconds'], raw['stage_calls']
    steps = raw['steps_per_repeat']
    if steps <= 0 or not len(elapsed) == len(stages) == len(counts) == repeats:
        raise RuntimeError('Invalid profile batch counts')
    totals = [0.] * (len(STAGES) + 1)
    for duration, values, calls in zip(elapsed, stages, counts):
        if calls != steps or len(values) != len(STAGES):
            raise RuntimeError('Profile did not cover every compute call')
        if not math.isfinite(duration) or duration <= 0 or any(not math.isfinite(x) or x < 0 for x in values):
            raise RuntimeError('Invalid profile timing')
        remainder = duration - sum(values)
        if remainder < -1e-9 * duration:
            raise RuntimeError('Profile stages exceed the enclosing batch time')
        for i, value in enumerate(values + [max(0., remainder)]):
            totals[i] += value
    total = sum(elapsed)
    names = STAGES + ('other',)
    raw.update(signatures_match=True,
               mean_ms={name: 1000 * value / (steps * repeats) for name, value in zip(names, totals)},
               share_percent={name: 100 * value / total for name, value in zip(names, totals)},
               total_mean_ms=1000 * total / (steps * repeats))
    return raw


def print_report(rows):
    names = STAGES + ('other',)
    for key, title, decimals in [('mean_ms', 'Profile stages (mean ms/call)', 3),
                                 ('share_percent', 'Profile stages (%)', 1)]:
        print('\n' + title)
        print('{:<20} {:>7} '.format('Case', 'Threads') + ''.join('{:>10}'.format(s) for s in LABELS) +
              ('{:>12}'.format('Total ms') if key == 'mean_ms' else ''))
        for row in rows:
            profile = row['profile']
            line = '{:<20} {:>7} '.format(row['case'], row['threads'])
            line += ''.join('{:>10.{precision}f}'.format(profile[key][name], precision=decimals) for name in names)
            if key == 'mean_ms':
                line += '{:>12.3f}'.format(profile['total_mean_ms'])
            print(line)
