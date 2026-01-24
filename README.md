# NEP_CPU

# What does this repository contain?

* A **standalone** C++ implementation (a class called `NEP`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD). We stress that **there is no external dependence**. The implementation works for all versions of NEP.

* An interface of the `NEP` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). **It can be run with MPI**.

# The standalone C++ implementation of NEP

* The `NEP` C++ class and the interfaces are defined in `src/nep.h`.
  
* The following folders contain some testing code and results:
  * `test/lammps/`
  * `test_nep/`
  * `test_qnep/`
  * `test_dftd3/`
  
* The `NEP` C++ class is used as an engine powering the following Python packages:
  * `calorine`: https://gitlab.com/materials-modeling/calorine
  * `PyNEP`: https://github.com/bigd4/PyNEP
  * `somd`: https://github.com/initqp/somd
  * `NepTrainKit`: https://github.com/aboys-cb/NepTrainKit
  
# The NEP-LAMMPS interface

Note: qNEP is not supported in the NEP-LAMMPS interface.

## Build the NEP-LAMMPS interface

* step 1: Copy the files in `src/` into `interface/lammps/USER-NEP/`.

Command:

```
cd ${software}/NEP_CPU
cp src/* interface/lammps/USER-NEP
```

* Step 2: Now you can copy the `USER-NEP/` folder into `YOUR_LAMMPS_PATH/src/` and start to compile LAMMPS in your favorite way.

  If you are compiling LAMMPS using `make`, ensure that you run "make yes-USER-NEP" before the final compilation.
  
```
cd ${software}/NEP_CPU
cp -r interface/lammps/USER-NEP ${software}/lammps/src

cd ${software}/lammps/src
make yes-USER-NEP
make machine
```

  If you are using `cmake`, you need to copy `/interface/lammps/USER-NEP.cmake` to LAMMPS cmake package directory and add `USER-NEP` to `CMakeLists.txt`. Then enable NEP by add "PKG_NEP=on", along with your other -D flags during the configuration step.

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
* Specify the `pair_style` in the same way as other potentials in LAMMPS (the first 2 arguments must be * * so as to span all atom types). For example, if you have a NEP model `NEP_HBCN.txt`, and your data file just have element carbon, you can set
  ```shell
  pair_style nep   
  pair_coeff * *  NEP_HBCN.txt C
  ```
  Firstly, we should set `pair_style` to `nep`, showed in the first line. Then we need set the NEP potential file and atom types by the command `pair_coeff`. Two asterisks `* *` mean every atom type will be set an element type or `NULL`. `NULL` means this potential doesn't consider the atom type. In this example, we set atom type `1` in LAMMPS data file to element `C` in NEP potential file. 
* The interface also supports multi-element system and hybrid potentials. Take a NEP model `NEP_PdCuNiP.txt` as an example. In this NEP model file, the first line is `nep3 4 Pd Cu Ni P`. Then in your LAMMPS input file, the next setting is allowed:
  ```shell
  pair_style hybrid/overlay nep nep lj/cut 1.0
  pair_coeff * * nep 1 NEP_PdCuNiP.txt Cu   Ni   NULL
  pair_coeff * * nep 2 NEP_PdCuNiP.txt NULL NULL Pd
  pair_ceoff 1*2 3 lj/cut 1.0 1.0
  ```
  The `pair_style` should be set `hybrid/overly` or other hybrid methods in LAMMPS. The hybrid potentials should be set after hybrid method. Then, in command `pair_coeff` we need set potential name again to identify which potential is setting for and the number of the potential if more than one. 
  Here, we set two NEP potentials. The first one just computes the NEP potential between atom type `1` `Cu` and `2` `Ni`, and of themselves. The second computes NEP potential of atom type `3` `Pd` itself.

* If you want to calculate the heat current correctly, use the following command to get the 9-component per-atom virial:
  ```shell
  compute 1 all centroid/stress/atom NULL
  ```
  
# Citation

* If you directly or indirectly use the `NEP` class here, you are suggested to cite the following paper:

  * Ke Xu, Hekai Bu, Shuning Pan, Eric Lindgren, Yongchao Wu, Yong Wang, Jiahui Liu, Keke Song, Bin Xu, Yifan Li, Tobias Hainer, Lucas Svensson, Julia Wiktor, Rui Zhao, Hongfu Huang, Cheng Qian, Shuo Zhang, Zezhu Zeng, Bohan Zhang, Benrui Tang, Yang Xiao, Zihan Yan, Jiuyang Shi, Zhixin Liang, Junjie Wang, Ting Liang, Shuo Cao, Yanzhou Wang, Penghua Ying, Nan Xu, Chengbing Chen, Yuwen Zhang, Zherui Chen, Xin Wu, Wenwu Jiang, Esme Berger, Yanlong Li, Shunda Chen, Alexander J. Gabourie, Haikuan Dong, Shiyun Xiong, Ning Wei, Yue Chen, Jianbin Xu, Feng Ding, Zhimei Sun, Tapio Ala-Nissila, Ari Harju, Jincheng Zheng, Pengfei Guan, Paul Erhart, Jian Sun, Wengen Ouyang, Yanjing Su, Zheyong Fan, [GPUMD 4.0: A high-performance molecular dynamics package for versatile materials simulations with machine-learned potentials]( https://doi.org/10.1002/mgea.70028), MGE Advances **3**, e70028 (2025).

* If you use the LAMMPS interface of the `NEP` class, a proper citation for LAMMPS is also suggested. 

