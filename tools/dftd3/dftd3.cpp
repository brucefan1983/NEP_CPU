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

#include "nep.h"
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

int main(int argc, char* argv[])
{
  if (argc != 4) {
    return 1;
  }
  std::ostream output(std::cout.rdbuf());
  std::cout.rdbuf(std::cerr.rdbuf());
  output << std::setprecision(17);
  const double cutoff = std::stod(argv[2]);
  const double cutoff_cn = std::stod(argv[3]);
  NEP nep;
  for (int n = 0; n < 94; ++n) {
    nep.dftd3.atomic_number[n] = n;
  }
  nep.set_dftd3_para_all(argv[1], cutoff, cutoff_cn);

  int N;
  while (std::cin >> N) {
    if (N <= 0) {
      return 1;
    }
    std::vector<int> type(N);
    std::vector<double> box(9), position(N * 3);
    std::vector<double> potential(N), force(N * 3), virial(N * 9);
    for (double& value : box) {
      std::cin >> value;
    }
    for (int n = 0; n < N; ++n) {
      std::cin >> type[n];
      --type[n];
      if (type[n] < 0 || type[n] >= 94) {
        return 1;
      }
      for (int d = 0; d < 3; ++d) {
        std::cin >> position[n + N * d];
      }
    }
    if (!std::cin) {
      return 1;
    }
    nep.compute_dftd3(argv[1], cutoff, cutoff_cn, type, box, position, potential, force, virial);
    output << std::accumulate(potential.begin(), potential.end(), 0.0L);
    for (int d = 0; d < 9; ++d) {
      output << " " << std::accumulate(virial.begin() + N * d, virial.begin() + N * (d + 1), 0.0L);
    }
    output << "\n";
    for (int n = 0; n < N; ++n) {
      output << force[n] << " " << force[n + N] << " " << force[n + N * 2] << "\n";
    }
    output.flush();
  }
  return std::cin.eof() ? 0 : 1;
}
