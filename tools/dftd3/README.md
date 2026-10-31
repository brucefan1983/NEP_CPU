# D3 dataset corrections

Add or subtract the D3(BJ) correction using NEP_CPU's existing D3 implementation:

```sh
python3 tools/dftd3/dftd3.py add input.xyz with_d3.xyz --functional pbe --cutoff 12 --cn-cutoff 6
python3 tools/dftd3/dftd3.py subtract with_d3.xyz without_d3.xyz --functional pbe --cutoff 12 --cn-cutoff 6
```

Requires Python 3.8+ (standard library only) and a C++11 compiler with OpenMP
support. The script compiles its C++ helper in a temporary directory and removes
it on exit. Use `--cxx` / `CXX` to select the compiler, `--threads` to select the
thread count (default: 1), or `--no-openmp` for a serial build.
No NEP model file is needed.

The input is extended XYZ with `Lattice`, total `energy`, and
`Properties=species:S:1:pos:R:3:forces:R:3`; reordered properties, extra columns,
and the singular name `force` are also accepted. Energy and forces use eV and
eV/Angstrom. Existing nine-component `virial` (eV) and `stress` (eV/Angstrom^3)
are corrected; absent tensors are not added. Other header fields and atom
columns are preserved. All three directions are periodic; omitted `pbc` is
treated as `T T T`.

Choose the functional and cutoffs used for the desired correction. The cutoff
arguments are in Angstrom; the coordination-number cutoff must not exceed the
energy cutoff. Adding and subtracting with the same settings are inverse
operations. Existing output files are not overwritten, and a failed conversion
does not leave a partial output dataset.

`example.xyz` is the original `test_dftd3/train.xyz` dataset. This tool replaces
the old executable's add/subtract modes (3 and 4).

Run all regression tests, including this tool's tests, with:

```sh
python3 tests/api/run_tests.py
```
