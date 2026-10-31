# C++ API regression tests

From the repository root:

```sh
python3 test/api/run_tests.py
```

Requires Python 3.8+ (standard library only) and a C++11 compiler with OpenMP
support (`g++` by default). No NumPy, MATLAB, GPU, GPUMD, LAMMPS, or Python
binding is required.
Python compiles a small C++ adapter, runs the API, and checks the results.
Both analytic and table builds, including the OpenMP consistency check, run by
default. Compilation and all generated files use a temporary directory, removed
on exit. Existing examples and their
reference files are never overwritten. A failed assertion, crash, or compilation
failure returns a nonzero exit status.

```sh
python3 test/api/run_tests.py --mode analytic
python3 test/api/run_tests.py --mode table
python3 test/api/run_tests.py --no-openmp
python3 test/api/run_tests.py --only typewise
```

`--cxx` / `CXX` select the compiler; `--cxxflags` / `CXXFLAGS` add compiler flags.
For example, `--cxxflags='-O0 -g'`. The default build adds `-fopenmp` and checks
one versus three threads. Other checks use one thread. For compilers without
OpenMP support, `--no-openmp` builds serially and skips only the thread comparison.
`--openmp` remains accepted explicitly. Paths are resolved from the script, so
it also runs from outside the repository.

## Coverage

| Check | Purpose |
| --- | --- |
| `data/nep`, `data/qnep` saved references | Every force component, all nine per-atom virial components, and every descriptor component |
| Saved finite-difference and GPU forces | Reproduce the original examples and compare with their GPUMD float32 results |
| Additional frozen results | Atomic energies, charges, BEC, typewise cutoffs, response models, latent space, B projection, ZBL, DFTD3, and the table build |
| Energy derivatives | Coordinate finite differences for forces and all nine homogeneous strains for total virial |
| All 64 many-body flag layouts | Selected force components and all strain components for ordinary NEP, qNEP1, and qNEP2 |
| Charge response | Neutrality, BEC sum rule, and finite differences of raw ANN charges for the GPUMD bond-based BEC |
| Symmetries | Translation, atom permutation, proper rotation of forces/tensors/dipole/polarizability |
| Periodic geometry | Triclinic cells, independent atom image shifts, small cell versus 3x3x3 supercell |
| Object reuse | Repeated computation, changing geometry and atom count on one initialized object |
| Typewise cutoffs | Each ordered type pair, equivalent uniform-cutoff model, and both sides of radial/angular cutoffs |
| ZBL | Independent Python pair-energy/force/virial formula; fixed, typewise, and flexible switching |
| DFTD3 | Frozen result, finite differences, and NEP + DFTD3 additivity |
| Feature API | Latent-space energy reconstruction and B-projection derivatives with respect to ANN parameters |
| Rejections | Removed/unknown models, second hidden layer, short ANN header, wrong compute overload, wrong output sizes, and the LAMMPS model-type guard |
| Version-1 model files | Record order, case, comments, numerical blocks, and rejection of malformed or unsupported records |
| Parameter layout | Three species, asymmetric ordered-pair coefficients, unequal radial/angular dimensions, exact parameter storage, independent pair descriptors, and force/virial/BEC consistency |
| Reference energies | Per-element energy shifts with unchanged forces, virials, charges, BEC, and features |
| LAMMPS API | Reversed element mapping, total and per-atom energies/virials, and optional per-atom output buffers |
| Rank diagnostics | Root-only output and nonzero exit, tested through the API rank flag |

The LAMMPS checks call the C++ API directly; this suite does not
run a LAMMPS executable or test MPI/ghost-atom integration. See `test/lammps` for
LAMMPS examples. Reinitializing an existing object with a second model is outside
this suite's intended usage.

## References and tolerances

The PbTe models, structures, historical `*_ref.out` files, and `force_gpu.out`
files are stored in `data/nep/` and `data/qnep/` next to the test scripts.
The carbon model used for D3 checks is in `data/dftd3/`.
`references.json` supplements them with results from
unmodified commit `9bc09df9fb96dc0bf3facad6b34ce02deeaf1ffa`. It records input
SHA-256 hashes and separate analytic/table results. Synthetic C/Si and C/Si/O fixtures have
nonuniform descriptor scalers and distinct ANN parameters for each species;
they exercise code paths, not physical accuracy.

Fixed CPU results use a componentwise `atol + rtol * abs(reference)` comparison,
normally with both tolerances `2e-10`. This permits floating-point roundoff
without treating changes as bitwise failures. Finite-difference tolerances are
larger because of subtraction and truncation errors. Values, shapes, and
finiteness are checked; comparisons never accept NaN/Inf.

The saved GPU force tolerances are absolute `1e-5 eV/Angstrom` for ordinary NEP
and `1.2e-4 eV/Angstrom` for qNEP. These are fixture-specific, based on the
preserved float32 results, not general CPU/GPU accuracy bounds. Table results
are checked against their own frozen results with the tight CPU tolerance;
the additional historical analytic comparison allows interpolation error.

The original qNEP finite-difference step, `2e-4 Angstrom`, crosses a cutoff for
some pairs. The suite reproduces those saved values and separately uses
`2e-6 Angstrom` to test force/energy consistency (the closest pair is only about
`8.85e-6 Angstrom` below the cutoff). The adapter sums atomic energies in
`long double` to reduce cancellation; API outputs remain doubles. Strict coordinate and strain
finite differences are applied to the analytic build: independently interpolated
function and derivative tables do not form an exact derivative pair.

There is no automatic reference-update option. Before intentionally changing
reference values, establish the cause, validate against an independent check
(e.g. GPUMD or finite differences), review the numerical changes, and update
the recorded source commit and input hashes as appropriate. A snapshot protects
existing behavior; the independent checks are needed to detect a consistently
wrong energy/force implementation.
