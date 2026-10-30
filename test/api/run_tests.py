#!/usr/bin/env python3
"""Run from anywhere: python3 /path/to/NEP_CPU/test/api/run_tests.py."""
import argparse
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
from support import Backend
import test_api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('analytic', 'table', 'both'), default='both')
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'g++'), help='C++ compiler (default: CXX or g++)')
    parser.add_argument('--cxxflags', default=os.environ.get('CXXFLAGS', ''), help='additional compiler flags')
    parser.add_argument('--openmp', action='store_true', help='compile with OpenMP and compare 1 vs 3 threads')
    parser.add_argument('--only', default='', help='run test methods whose names contain this string')
    args = parser.parse_args()
    modes = ('analytic', 'table') if args.mode == 'both' else (args.mode,)
    tests = [name for name in unittest.defaultTestLoader.getTestCaseNames(test_api.APITests) if args.only in name]
    if not tests:
        parser.error('No matching tests')
    try:
        with tempfile.TemporaryDirectory(prefix='nep_api_tests_') as tmp:
            backends = {mode: Backend(tmp, mode, args.cxx, args.cxxflags, args.openmp) for mode in modes}
            test_api.configure(Path(tmp), backends)
            suite = unittest.TestSuite(test_api.APITests(name) for name in tests)
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            return 0 if result.wasSuccessful() else 1
    except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
