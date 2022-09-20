# NEP_CPU

# What does this repository contain?

* A **standalone** C++ implementation (a class called `NEP3`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD). We stress that **there is no external dependence**.

* An interface of the `NEP3` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). **It can be run with MPI**.

# The standalone C++ implementation of NEP

* The `NEP3` C++ class is defined in the following two files:
  * `src/nep.h`
  * `src/nep.cpp`
  
* The following file is used to test the above class:
  * `src/main.cpp`
  
* The following folder contains some testing results:
  * `test/`
  
* The `NEP3` C++ class is used as an engine powering the following two Python packages:
  * `PyNEP`: https://github.com/bigd4/PyNEP
  * `calorine`: https://gitlab.com/materials-modeling/calorine
  
# The NEP-LAMMPS interface

## Build the NEP-LAMMPS interface

* step 1: Copy `src/nep.h` and `src/nep.cpp` into `interface/lammps/USER-NEP/` such that you have the following files in `interface/lammps/USER-NEP/`:
  * `nep.h`
  * `nep.cpp`
  * `pair_NEP.h`
  * `pair_NEP.cpp`

* Step 2: Check the version of your LAMMPS to be installed. Translate it into a number such as 20220324 (year, month, day). Then check the beginning of `pair_NEP.cpp`, where you can find the following line
```
#define LAMMPS_VERSION_NUMBER 20220324 // use the new neighbor list starting from this version
```
If your LAMMPS version is not the one written here, change this line to use your version number.

* Step 3: Now you can copy the `USER-NEP/` folder into `YOUR_LAMMPS_PATH/src/` and start to compile LAMMPS in your favorite way. Good luck!
  
## Use the NEP-LAMMPS interface

* `atom_style` can only be `atomic`
* `units` must be `metal`
* Specify the `pair_style` in the following way:
  ```shell
  pair_style nep YOUR_NEP_MODEL_FILE.txt  # YOUR_NEP_MODEL_FILE.txt is your NEP model file (with path)
  pair_coeff * *                          # This format is fixed
  ```
  
* For multi-element system, the atom types must be carefully set. Take a NEP model `NEP_PdCuNiP.txt` as an example. In this NEP model file, the first line is `nep3 4 Pd Cu Ni P`. Then in your LAMMPS input file, you must set 
  * Pb atoms to type 1
  * Cu atoms to type 2
  * Ni atoms to type 3
  * P atoms to type 4
  
* If you want to calculate the heat current correctly, use the following command to get the 9-component per-atom virial:
  ```shell
  compute 1 all centroid/stress/atom NULL
  ```
  
# Citation

* If you directly or indirectly use the `NEP3` class here, you are suggested to cite the following paper:

  * Zheyong Fan, Yanzhou Wang, Penghua Ying, Keke Song, Junjie Wang, Yong Wang, Zezhu Zeng, Ke Xu, Eric Lindgren, J. Magnus Rahm, Alexander J. Gabourie, Jiahui Liu, Haikuan Dong, Jianyang Wu, Yue Chen, Zheng Zhong, Jian Sun, Paul Erhart, Yanjing Su, Tapio Ala-Nissila,
[GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617), The Journal of Chemical Physics. (2022) In press.

* If you use the LAMMPS interface of the `NEP3` class, a proper citation for LAMMPS is also suggested. 

