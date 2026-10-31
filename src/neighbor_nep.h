/*
    Copyright 2017 Zheyong Fan and GPUMD development team
    This file is part of GPUMD.
    GPUMD is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GPUMD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with GPUMD.  If not, see <https://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
A CPU implementation of the neuroevolution potential (NEP)
Ref: Zheyong Fan et al., Neuroevolution machine learning potentials:
Combining high accuracy and low cost in atomistic simulations and application to
heat transport, Phys. Rev. B. 104, 104309 (2021).
------------------------------------------------------------------------------*/

#pragma once
#include <vector>

void find_neighbor_list_small_box(
  const double rc_radial,
  const double rc_angular,
  const int N,
  const int MN, 
  const std::vector<double>& box,
  const std::vector<double>& position,
  int* num_cells,
  double* ebox,
  std::vector<int>& g_NN_radial,
  std::vector<int>& g_NL_radial,
  std::vector<int>& g_NN_angular,
  std::vector<int>& g_NL_angular,
  std::vector<double>& r12);
