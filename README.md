# NEP_CPU

* This is a serial CPU implementation (a standalone C++ class) of the NEP machine-learned potential as introduced in the GPUMD package (https://github.com/brucefan1983/GPUMD)

* `NEP_CPU` is mainly used as an engine powering the following two Python packages:
  * `PyNEP`: https://github.com/bigd4/PyNEP
  * `calorine`: https://gitlab.com/materials-modeling/calorine
  
* `NEP_CPU` is not designed for large-scale atomistic simulations. For this purpose, please `GPUMD` is the better choice.
