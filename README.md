# NEP_CPU

NEP_CPU is part of the [GPUMD](https://github.com/brucefan1983/GPUMD) project. It provides a standalone CPU implementation of GPUMD's NEP model evaluation routines for use in Python packages and LAMMPS.

* A **standalone** C++ class (`NEP`) with **no external dependencies**.

* An interface of the `NEP` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). **It can be run with MPI**.

## Model files

Only model files with `format_version 1` are accepted. See the [GPUMD format specification](https://github.com/brucefan1983/GPUMD/blob/master/potentials/nep/nep_model_format.md).

Convert an older model file using GPUMD's [nep_model_converter.py](https://github.com/brucefan1983/GPUMD/blob/master/potentials/nep/nep_model_converter.py):

```sh
python3 /path/to/GPUMD/potentials/nep/nep_model_converter.py nep_old.txt nep.txt
```

## C++ library

The `NEP` class and its interfaces are defined in [`src/nep.h`](src/nep.h).

### Build and install with CMake

With CMake 3.16 or later:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --parallel
cmake --install build
```

The default is a static library with analytic radial functions and no OpenMP.
Add `-DBUILD_SHARED_LIBS=ON` for a shared library, `-DNEP_CPU_OPENMP=ON` for
OpenMP, or `-DNEP_CPU_USE_TABLE=ON` for tabulated radial functions.

Downstream CMake projects can use:

```cmake
find_package(NEP_CPU CONFIG REQUIRED)
target_link_libraries(program PRIVATE nep::nep_cpu)
```

Configure the downstream project with `-DCMAKE_PREFIX_PATH=/path/to/install`
and include `nep.h` in its C++ code. The target supplies the required include
directories, compile definitions, and link dependencies. The same target is
also available when adding NEP_CPU with `add_subdirectory`.

### Compile directly with your program

Compile your program together with `src/nep.cpp`, `src/nep_model_file.cpp`, `src/ewald_nep.cpp`, and `src/neighbor_nep.cpp`, using C++11 or later.

### Python packages

The `NEP` class powers the following Python packages:

* `calorine`: https://gitlab.com/materials-modeling/calorine
* `PyNEP`: https://github.com/bigd4/PyNEP
* `somd`: https://github.com/initqp/somd
* `NepTrainKit`: https://github.com/aboys-cb/NepTrainKit

## LAMMPS interface

The NEP-LAMMPS interface supports NEP potential models, including their ZBL variants, with one hidden layer. qNEP, dipole, and polarizability models are not supported.

### Build

* Step 1: Copy the files in `src/` into `interface/lammps/USER-NEP/`.

Command:

```
cd ${software}/NEP_CPU
cp src/* interface/lammps/USER-NEP
```

* Step 2: Now you can copy the `USER-NEP/` folder into `YOUR_LAMMPS_PATH/src/` and start to compile LAMMPS in your favorite way.

  If you are compiling LAMMPS using `make`, run `make yes-USER-NEP` before the final compilation. Replace `machine` below with the desired LAMMPS build target.
  
```
cd ${software}/NEP_CPU
cp -r interface/lammps/USER-NEP ${software}/lammps/src

cd ${software}/lammps/src
make yes-USER-NEP
make machine
```

  If you are using `cmake`, copy `interface/lammps/USER-NEP.cmake` to the LAMMPS CMake package directory and add `USER-NEP` to `CMakeLists.txt`. Then enable the package with `-D PKG_USER-NEP=on`, along with any other CMake options.

```
cd ${software}/NEP_CPU
cp -r interface/lammps/USER-NEP ${software}/lammps/src
cp interface/lammps/USER-NEP.cmake ${software}/lammps/cmake/Modules/Packages/

cd ${software}/lammps/cmake
sed -i '/foreach(PKG_WITH_INCL / s/)/ USER-NEP)/' CMakeLists.txt
sed -i '/set(STANDARD_PACKAGES/,/)/ s/)/  \n  USER-NEP)/' CMakeLists.txt

cd ..; mkdir build; cd build
cmake -D PKG_USER-NEP=on ../cmake
cmake --build .
make install
```

### Usage

* `atom_style` can be `atomic` and `full`
* `units` must be `metal`
* Set `pair_style` to `nep`. The first two arguments of `pair_coeff` must be `* *`, followed by the model file and one element mapping for each LAMMPS atom type. For example, if you have a NEP model `NEP_HBCN.txt` and your data file contains only carbon, use:
  ```shell
  pair_style nep   
  pair_coeff * *  NEP_HBCN.txt C
  ```
  List the element symbols in LAMMPS atom-type order. Use `NULL` for types handled by other sub-styles in a hybrid setup. In this example, LAMMPS atom type `1` maps to element `C` in the NEP model.
* The interface also supports multi-element systems and hybrid potentials. For a model `NEP_PdCuNiP.txt` containing Pd, Cu, Ni, and P, an example is:
  ```shell
  pair_style hybrid/overlay nep nep lj/cut 1.0
  pair_coeff * * nep 1 NEP_PdCuNiP.txt Cu   Ni   NULL
  pair_coeff * * nep 2 NEP_PdCuNiP.txt NULL NULL Pd
  pair_coeff 1*2 3 lj/cut 1.0 1.0
  ```
  List the sub-styles after `hybrid/overlay`. In each `pair_coeff` command, specify the sub-style name and its instance number if that sub-style appears more than once.
  Here, the first NEP instance describes the subsystem of Cu and Ni atoms (types `1` and `2`), and the second describes the Pd atoms (type `3`). The `lj/cut` term describes interactions between these two subsystems.

* To calculate the potential (virial) contribution to the heat current, use the 9-component per-atom virial from `compute centroid/stress/atom`.

See [`tests/lammps/`](tests/lammps/) for interface tests.

## D3 dataset corrections

Use the [Python D3 tool](tools/dftd3/README.md) to add or subtract D3(BJ)
energy, force, and virial corrections in extended XYZ datasets.

## Tests and benchmarks

* [API regression tests](tests/api/README.md)
* [Performance benchmarks](benchmarks/README.md)

## Citation

* As NEP_CPU is part of GPUMD, please cite the following GPUMD 4.0 paper when using NEP_CPU:

  * Ke Xu, Hekai Bu, Shuning Pan, Eric Lindgren, Yongchao Wu, Yong Wang, Jiahui Liu, Keke Song, Bin Xu, Yifan Li, Tobias Hainer, Lucas Svensson, Julia Wiktor, Rui Zhao, Hongfu Huang, Cheng Qian, Shuo Zhang, Zezhu Zeng, Bohan Zhang, Benrui Tang, Yang Xiao, Zihan Yan, Jiuyang Shi, Zhixin Liang, Junjie Wang, Ting Liang, Shuo Cao, Yanzhou Wang, Penghua Ying, Nan Xu, Chengbing Chen, Yuwen Zhang, Zherui Chen, Xin Wu, Wenwu Jiang, Esme Berger, Yanlong Li, Shunda Chen, Alexander J. Gabourie, Haikuan Dong, Shiyun Xiong, Ning Wei, Yue Chen, Jianbin Xu, Feng Ding, Zhimei Sun, Tapio Ala-Nissila, Ari Harju, Jincheng Zheng, Pengfei Guan, Paul Erhart, Jian Sun, Wengen Ouyang, Yanjing Su, Zheyong Fan, [GPUMD 4.0: A high-performance molecular dynamics package for versatile materials simulations with machine-learned potentials]( https://doi.org/10.1002/mgea.70028), MGE Advances **3**, e70028 (2025).

* If you use the LAMMPS interface of the `NEP` class, a proper citation for LAMMPS is also suggested. 

## License

The NEP_CPU core is licensed under LGPL-3.0-or-later; see the [LGPL v3](LICENSES/LGPL-3.0.txt) and [GPL v3](LICENSES/GPL-3.0.txt) license texts. The LAMMPS interface retains its existing GPL license and copyright notices.
