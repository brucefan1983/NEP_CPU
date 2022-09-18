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
  
# The NEP-LAMMPS interface

* Applies to LAMMPS that is no older than the [20220324 version](https://github.com/lammps/lammps/releases/tag/patch_24Mar2022).

* The NEP-LAMMPS interface consists of the following files:
  * `interface/lammps/install.sh`
  * `interface/lammps/USER-NEP/pair_NEP.h`
  * `interface/lammps/USER-NEP/pair_NEP.cpp`
  * `interface/lammps/USER-NEP/install.sh`
  
* How to install?
  * Modify `interface/lammps/install.sh`:
    * Change `$1` to the path of your LAMMPS package.
    * Change `serial` to `mpi`, `lmp_serial` to `lmp_mpi`, and `clean-serial` to `clean-mpi` if you want to build an MPI version.
  * Run the `interface/lammps/install.sh` file to install.
  
* Tips for using NEP in LAMMPS:
  ```
  atom_style atomic
  units metal
  pair_style nep YOUR_NEP_MODEL_FILE.txt
  pair_coeff * *
  ```

* Speed information: NEP-LAMMPS with 1 CPU can achieve about 1/1000 of the speed of NEP-GPUMD with one A100 GPU.

# Citation

* If you directly or indirectly use the `NEP3` class here, you are suggested to cite the following paper:

  * Zheyong Fan, Yanzhou Wang, Penghua Ying, Keke Song, Junjie Wang, Yong Wang, Zezhu Zeng, Ke Xu, Eric Lindgren, J. Magnus Rahm, Alexander J. Gabourie, Jiahui Liu, Haikuan Dong, Jianyang Wu, Yue Chen, Zheng Zhong, Jian Sun, Paul Erhart, Yanjing Su, Tapio Ala-Nissila,
[GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617), The Journal of Chemical Physics. (2022) In press.

* If you use the LAMMPS interface of the `NEP3` class, a proper citation for LAMMPS is also suggested. 

