# API performance benchmarks

From the repository root:

```sh
python3 benchmarks/run_benchmarks.py
```

Requires Python 3.8+ (standard library only) and a C++11 compiler with OpenMP
support (`g++` by default). No Python binding, LAMMPS, MPI, or GPU is needed.
The default build uses `-O3`, double precision, analytic radial functions,
and one OpenMP thread. Compilation and generated inputs use a temporary
directory that is removed on exit. Results are saved in `benchmarks/results/`.
Existing result files are never overwritten.

| Case | Atoms | Model and purpose |
| --- | ---: | --- |
| `diamond_1728`, `diamond_13824` | 1,728 / 13,824 | Carbon model from `tests/lammps/diamond`; two workload sizes related to issue #75 |
| `pbte_250`, `pbte_2000` | 250 / 2,000 | Real two-species model and nonorthogonal structure from `tests/api/data/nep`; original cell and 2x2x2 replication |
| `qnep_250` | 250 | Real model and structure from `tests/api/data/qnep`; includes electrostatics and BEC through the charge API |
| `synthetic_l6_1728`, `synthetic_l8_1728` | 1,728 | Deterministic two-species models, `n_max 4 4`, `basis_size 8 8`, 30 neurons; exercise high angular orders, not physical predictions |

Each case uses eight structures with deterministic, periodic displacements of
amplitude 0.01 Angstrom. One process initializes one NEP object and reuses its
output buffers across all calls. This follows repeated use from a Python
program, but measures the C++ API itself, excluding Python binding overhead.
It does not integrate an MD trajectory.

Timing uses `std::chrono::steady_clock` around batches of complete `compute`
calls, including the API's neighbor construction and output initialization.
File I/O, model initialization, first allocation, and an initial warmup cycle
are excluded. A separate calibration cycle selects a whole number of eight-frame
cycles targeting at least one second per batch. Five timed batches are reported
as median milliseconds per call and thousands of atom-steps per second.
`Spread %` is `(maximum - minimum) / median`; repeat noisy runs on an idle machine.

Save a baseline, then compare after a code change:

```sh
python3 benchmarks/run_benchmarks.py --output benchmarks/results/baseline.json
python3 benchmarks/run_benchmarks.py --compare benchmarks/results/baseline.json --output benchmarks/results/optimized.json
```

`Speedup` is baseline median time divided by current median time. Comparisons
require matching machine information, compiler and flags, thread count, OpenMP
settings, driver, model, and structure hashes. Source hashes and Git state are
recorded but may change between versions. Calibration can choose different
batch lengths; every batch samples the same eight frames equally.

Optional examples:

```sh
python3 benchmarks/run_benchmarks.py --threads 1 2 4 --output threads.json
python3 benchmarks/run_benchmarks.py --mode table --output table.json
python3 benchmarks/run_benchmarks.py --only pbte_250 --seconds 2 --repeats 5
```

`--cxx` / `CXX` select the compiler, and `--cxxflags` / `CXXFLAGS` add flags.
`OMP_DYNAMIC` is disabled and nested teams are disabled. Unless already set,
`OMP_PROC_BIND=close` and `OMP_PLACES=cores` bind threads to cores. The observed
thread count is checked. Use the same CPU allocation and settings for comparisons.

Untimed numerical signatures cover each frame: mean atomic energy, force RMS,
maximum absolute force component, nine mean virial components, and, for qNEP,
charge/BEC RMS. Nonfinite outputs or changes after repeated calls fail the run.
A baseline signature mismatch is shown as `DIFF` and returns a nonzero status.
These signatures are sanity checks; use `python3 tests/api/run_tests.py` for
correctness regression tests. Benchmarks are not part of that default test suite.

## Stage profiling

```sh
python3 benchmarks/run_benchmarks.py --profile --output benchmarks/results/profile.json
```

`--profile` adds a separate instrumented run for each case and thread count.
The usual throughput table and `--compare` use the original uninstrumented
driver, so existing baseline files remain usable. Profiling builds temporary
copies of the two `compute` overloads and the benchmark driver; production
sources and the API are unchanged. Profiling checks the expected stage call
sites and their order, and stops if these anchors are missing or ambiguous.

Two extra tables report mean milliseconds per call and percentages:

| Column | Work included |
| --- | --- |
| `Setup` | Input checks, buffer preparation, and output initialization |
| `Neighbor` | Radial/angular neighbor construction and relative displacements |
| `Desc/ANN` | Descriptors, neural networks, and descriptor derivatives |
| `Radial`, `Angular` | Radial and angular forces and virials |
| `BEC` | qNEP charge neutrality correction, BEC calculation, and BEC scaling |
| `Coulomb` | qNEP Ewald calculation, real-space or self-energy correction, and mean-charge-derivative correction |
| `ZBL` | Optional ZBL block |
| `Other` | Remainder of the enclosing batch, including calls, timer bookkeeping, and the benchmark loop |

Warmup and calibration are excluded from stage statistics. Times are wall times,
including completion of OpenMP work, not sums of worker CPU times. Stage times
and percentages use all calls in the profiled batches and add up to the profiled
total. This total can differ from the normal table's median because profiling
is a separate run and instrumentation can affect generated code and timing.
Use the normal table to judge speedups; use stage shares to locate bottlenecks.

Each profiled run must match the normal run's numerical signatures and account
for every timed `compute` call. Raw batch times, stage times, call counts, and
instrumentation hashes are saved in the JSON result.
