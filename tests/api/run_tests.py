#!/usr/bin/env python3
"""Run from anywhere: python3 /path/to/NEP_CPU/tests/api/run_tests.py."""
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
import test_dftd3


class TestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.passed = 0
        self.failed_tests = set()

    def addSuccess(self, test):
        super().addSuccess(test)
        self.passed += 1

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.failed_tests.add(test)

    def addError(self, test, err):
        super().addError(test, err)
        self.failed_tests.add(test)

    def addSubTest(self, test, subtest, err):
        super().addSubTest(test, subtest, err)
        if err is not None:
            self.failed_tests.add(test)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('analytic', 'table', 'both'), default='both')
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'g++'), help='C++ compiler (default: CXX or g++)')
    parser.add_argument('--cxxflags', default=os.environ.get('CXXFLAGS', ''), help='additional compiler flags')
    openmp = parser.add_mutually_exclusive_group()
    openmp.add_argument('--openmp', action='store_true', default=True,
                        help='compile with OpenMP and compare 1 vs 3 threads (default)')
    openmp.add_argument('--no-openmp', dest='openmp', action='store_false',
                        help='compile without OpenMP and skip the thread comparison')
    parser.add_argument('--only', default='', help='run test methods whose names contain this string')
    args = parser.parse_args()
    modes = ('analytic', 'table') if args.mode == 'both' else (args.mode,)
    tests = {cls: [name for name in unittest.defaultTestLoader.getTestCaseNames(cls) if args.only in name]
             for cls in (test_api.APITests, test_dftd3.D3ToolTests)}
    if not any(tests.values()):
        parser.error('No matching tests')
    try:
        with tempfile.TemporaryDirectory(prefix='nep_api_tests_') as tmp:
            if tests[test_api.APITests]:
                backends = {mode: Backend(tmp, mode, args.cxx, args.cxxflags, args.openmp) for mode in modes}
                test_api.configure(Path(tmp), backends)
            if tests[test_dftd3.D3ToolTests]:
                test_dftd3.configure(Path(tmp), args.cxx, args.cxxflags, args.openmp)
            suite = unittest.TestSuite(cls(name) for cls, names in tests.items() for name in names)
            result = unittest.TextTestRunner(verbosity=2, resultclass=TestResult).run(suite)
            print('\n{} passed, {} failed, {} skipped'.format(
                result.passed, len(result.failed_tests), len(result.skipped)), file=sys.stderr)
            return 0 if result.wasSuccessful() else 1
    except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
