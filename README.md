# NEP_CPU

# What does this repository contain?

* A **standalone** C++ implementation (a class called `NEP3`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD). We stress that **there is no externel dependence**.

* An interface of the `NEP3` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). **It can be run with MPI**.

# File organization

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
  
* The LAMMPS interface code can be found in the folder `interface/lammps/USER-NEP/`, which contains the following files:
  * `pair_NEP.h`
  * `pair_NEP.cpp`
  * `install.sh`

# Speed information

* The computational speed of NEP in LAMMPS with one CPU thread is about 1/1000 of that of NEP in GPUMD with one A100 GPU.

# Citation

* If you directly or indirectly use the `NEP3` class here, you are suggested to cite the following paper:

  * Zheyong Fan, Yanzhou Wang, Penghua Ying, Keke Song, Junjie Wang, Yong Wang, Zezhu Zeng, Ke Xu, Eric Lindgren, J. Magnus Rahm, Alexander J. Gabourie, Jiahui Liu, Haikuan Dong, Jianyang Wu, Yue Chen, Zheng Zhong, Jian Sun, Paul Erhart, Yanjing Su, Tapio Ala-Nissila,
[GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617), The Journal of Chemical Physics. (2022) In press.

* If you use the LAMMPS interface of NEP, a proper citation for LAMMPS is also suggested. 

