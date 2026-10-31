# NEP_CPU

# What does this repository contain?

* A **standalone** C++ implementation (a class called `NEP`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD). We stress that **there is no external dependence**.

* An interface of the `NEP` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). **It can be run with MPI**.

# Model files

Only model files with `format_version 1` are accepted. See the [GPUMD format specification](https://github.com/brucefan1983/GPUMD/blob/master/potentials/nep/nep_model_format.md).

Convert an older model file using GPUMD's [nep_model_converter.py](https://github.com/brucefan1983/GPUMD/blob/master/potentials/nep/nep_model_converter.py):

```sh
python3 /path/to/GPUMD/potentials/nep/nep_model_converter.py nep_old.txt nep.txt
```

# The standalone C++ implementation of NEP

* The `NEP` C++ class and the interfaces are defined in `src/nep.h`.

* Compile your program together with `src/nep.cpp`, `src/nep_model_file.cpp`, `src/ewald_nep.cpp`, and `src/neighbor_nep.cpp`, using C++11 or later.
  
* The following folders contain some testing code and results:
  * `test/api/`
  * `test/lammps/`
  * `test_dftd3/`
  
* The `NEP` C++ class is used as an engine powering the following Python packages:
  * `calorine`: https://gitlab.com/materials-modeling/calorine
  * `PyNEP`: https://github.com/bigd4/PyNEP
  * `somd`: https://github.com/initqp/somd
  * `NepTrainKit`: https://github.com/aboys-cb/NepTrainKit
  
# The NEP-LAMMPS interface

The NEP-LAMMPS interface supports NEP potential models, including their ZBL variants, with one hidden layer. qNEP, dipole, and polarizability models are not supported.

## Build the NEP-LAMMPS interface

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

## Use the NEP-LAMMPS interface

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
  
# Citation

* If you directly or indirectly use the `NEP` class here, you are suggested to cite the following paper:

  * Ke Xu, Hekai Bu, Shuning Pan, Eric Lindgren, Yongchao Wu, Yong Wang, Jiahui Liu, Keke Song, Bin Xu, Yifan Li, Tobias Hainer, Lucas Svensson, Julia Wiktor, Rui Zhao, Hongfu Huang, Cheng Qian, Shuo Zhang, Zezhu Zeng, Bohan Zhang, Benrui Tang, Yang Xiao, Zihan Yan, Jiuyang Shi, Zhixin Liang, Junjie Wang, Ting Liang, Shuo Cao, Yanzhou Wang, Penghua Ying, Nan Xu, Chengbing Chen, Yuwen Zhang, Zherui Chen, Xin Wu, Wenwu Jiang, Esme Berger, Yanlong Li, Shunda Chen, Alexander J. Gabourie, Haikuan Dong, Shiyun Xiong, Ning Wei, Yue Chen, Jianbin Xu, Feng Ding, Zhimei Sun, Tapio Ala-Nissila, Ari Harju, Jincheng Zheng, Pengfei Guan, Paul Erhart, Jian Sun, Wengen Ouyang, Yanjing Su, Zheyong Fan, [GPUMD 4.0: A high-performance molecular dynamics package for versatile materials simulations with machine-learned potentials]( https://doi.org/10.1002/mgea.70028), MGE Advances **3**, e70028 (2025).

* If you use the LAMMPS interface of the `NEP` class, a proper citation for LAMMPS is also suggested. 

