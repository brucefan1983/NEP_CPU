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
  
## Use the NEP-LAMMPS interface

* `units` must be `metal`
  
* Specify the `pair_style` in the following way:
  
  ```shell
  A.for pure nep potential
  pair_style nep 
  pair_coeff * * YOUR_NEP_MODEL_FILE.txt Pd Cu Ni P  # YOUR_NEP_MODEL_FILE.txt is your NEP model file (with path)
  
  B.for hybrid potential （e.g. NEP +  LJ/CUT)
  pair_style 	hybrid  nep 	 lj/cut 12.0   
  pair_coeff  * *     nep nep.txt  Pd Cu Ni P  NULL NULL NULL NULL   NULL # NULL means the type that is not described by NEP
  ```
  
  The **./example** directory contains detailed usage examples.
  
* For multi-element system, the atom types must be carefully set. Take a NEP model `NEP_PdCuNiP.txt` as an example. In this NEP model file, the first line is `nep3 4 Pd Cu Ni P`. Then in your LAMMPS input file, you must set 

  * Pd atoms to type 1
  * Cu atoms to type 2
  * Ni atoms to type 3
  * P atoms to type 4

* Some atom types can be missing in the simulated system. For example you can use the above NEP model to simulate a CuNi binary alloy. It is important to make sure to still set Cu atoms to type 2 and Ni atoms to type 3. In this case, atom types 1 and 4 are missing.

* One can hybrid NEP with other potentials in LAMMPS.

* If you want to calculate the heat current correctly, use the following command to get the 9-component per-atom virial:
  ```shell
  compute 1 all centroid/stress/atom NULL
  ```


# Citation

* If you directly or indirectly use the `NEP3` class here, you are suggested to cite the following paper:

  * Zheyong Fan, Yanzhou Wang, Penghua Ying, Keke Song, Junjie Wang, Yong Wang, Zezhu Zeng, Ke Xu, Eric Lindgren, J. Magnus Rahm, Alexander J. Gabourie, Jiahui Liu, Haikuan Dong, Jianyang Wu, Yue Chen, Zheng Zhong, Jian Sun, Paul Erhart, Yanjing Su, Tapio Ala-Nissila,
[GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617), The Journal of Chemical Physics **157**, 114801 (2022).

* If you use the LAMMPS interface of the `NEP3` class, a proper citation for LAMMPS is also suggested. 

