# NEP_CPU

* This repository contains:
  * A C++ implementation (a class called `NEP3`) of the NEP (neuroevolution potential) machine-learned potential as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD).
  * An interface of NEP to the CPU version LAMMPS (https://github.com/lammps/lammps). Can run with MPI.

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
  
* The standalone C++ class is also used as an engine powering the interface to LAMMPS. The interface code can be found in the folder `interface/lammps/USER-NEP/`, which contains the following files:
  * `pair_NEP.h`
  * `pair_NEP.cpp`
  * `install.sh`
