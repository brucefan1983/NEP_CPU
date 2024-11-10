# NEP_CPU

# What does this repository contain?

* A **standalone** C++ implementation (a class called `NEP3`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD). We stress that **there is no external dependence**. Even though the class name is `NEP3`, the implementation works for all versions of NEP.

* An interface of the `NEP3` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). **It can be run with MPI**.

# The standalone C++ implementation of NEP

* The `NEP3` C++ class is defined in the following three files:
  * `src/nep.h`
  * `src/nep.cpp`
  * `src/dftd3para.h`
  
* There is an option to use tables to speed up the calculations for the radial functions in NEP. To enable it, one can change line 20 of `src/nep.h`:

```
// #define USE_TABLE_FOR_RADIAL_FUNCTIONS
```
  
* The following folders contain some testing code and results:
  * `test/`
  * `test_dftd3/`
  
* The `NEP3` C++ class is used as an engine powering the following Python packages:
  * `calorine`: https://gitlab.com/materials-modeling/calorine
  * `PyNEP`: https://github.com/bigd4/PyNEP
  * `somd`: https://github.com/initqp/somd
  
# The NEP-LAMMPS interface

## Build the NEP-LAMMPS interface

* step 1: Copy the files in `src/` into `interface/lammps/USER-NEP/` such that you have the following files in `interface/lammps/USER-NEP/`:
  * `nep.h`
  * `nep.cpp`
  * `dftd3para.h`
  * `pair_NEP.h`
  * `pair_NEP.cpp`
  * `install.sh`
  
* Step 2: Now you can copy the `USER-NEP/` folder into `YOUR_LAMMPS_PATH/src/` and start to compile LAMMPS in your favorite way. Good luck!
  * Reminder: If you are compiling LAMMPS using `make`, ensure that you run "make yes-NEP" before the final compilation. For cmake, please add "PKG_NEP=on" to enable NEP, along with your other -D flags during the configuration step.
  
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
  pair_style hybrid/overlay nep nep ij/cut 1.0
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

* If you directly or indirectly use the `NEP3` class here, you are suggested to cite the following paper:

  * Zheyong Fan, Yanzhou Wang, Penghua Ying, Keke Song, Junjie Wang, Yong Wang, Zezhu Zeng, Ke Xu, Eric Lindgren, J. Magnus Rahm, Alexander J. Gabourie, Jiahui Liu, Haikuan Dong, Jianyang Wu, Yue Chen, Zheng Zhong, Jian Sun, Paul Erhart, Yanjing Su, Tapio Ala-Nissila,
[GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617), The Journal of Chemical Physics **157**, 114801 (2022).

* If you use the LAMMPS interface of the `NEP3` class, a proper citation for LAMMPS is also suggested. 

