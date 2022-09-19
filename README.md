# NEP_CPU

# What does this repository contain?

* A **standalone** C++ implementation (a class called `NEP3`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD). We stress that **there is no externel dependence**.

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
  
# Build the NEP-LAMMPS interface

* step 1: copy `src/nep.h` and `src/nep.cpp` into `interface/lammps/USER-NEP/` such that you have the following files in `interface/lammps/USER-NEP/`:
  * `nep.h`
  * `nep.cpp`
  * `pair_NEP.h`
  * `pair_NEP.cpp`

* Step 2: Check the version of your LAMMPS to be installed. Translate it into a number such as 20220324 (year, month, day). Then check the begining of `pair_NEP.cpp`, where you can find the following line
```
#define LAMMPS_VERSION_NUMBER 20220324 // use the new neighbor list starting from this version
```
If your LAMMPS version is not the one written here, change it to your version number.

* Step 3: Now you can copy the `USER-NEP/` folder into your LAMMPS package and start compiling LAMMPS in your favorite way. Good luck!
  
* Step 4: Start to use NEP in LAMMPS:
  ```
  atom_style atomic                       # Can only be atomic
  units metal                             # Can only be metal
  pair_style nep YOUR_NEP_MODEL_FILE.txt  # Put your NEP potential file in the current working directory
  pair_coeff * *                          # This format is fixed
  ```

* Speed information: NEP-LAMMPS with 1 CPU can achieve about 1/1000 of the speed of NEP-GPUMD with one A100 GPU.

# Citation

* If you directly or indirectly use the `NEP3` class here, you are suggested to cite the following paper:

  * Zheyong Fan, Yanzhou Wang, Penghua Ying, Keke Song, Junjie Wang, Yong Wang, Zezhu Zeng, Ke Xu, Eric Lindgren, J. Magnus Rahm, Alexander J. Gabourie, Jiahui Liu, Haikuan Dong, Jianyang Wu, Yue Chen, Zheng Zhong, Jian Sun, Paul Erhart, Yanjing Su, Tapio Ala-Nissila,
[GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617), The Journal of Chemical Physics. (2022) In press.

* If you use the LAMMPS interface of the `NEP3` class, a proper citation for LAMMPS is also suggested. 

