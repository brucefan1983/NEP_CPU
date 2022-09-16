# NEP_CPU

* This repository contains:
  * A C++ implementation (a class called `NEP3`) of the neuroevolution potential (NEP) as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD).
  * An interface of the `NEP3` class to the CPU version of LAMMPS (https://github.com/lammps/lammps). It can be run with MPI.

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

* The computational speed of NEP in LAMMPS with one CPU thread is about 1/1000 of that of NEP in GPUMD with one A100 GPU.
