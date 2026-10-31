// Timed standalone API calls; Python prepares workloads and reports comparisons.
#include "nep.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace
{
void check_finite(const std::vector<double>& values)
{
  for (double value : values) {
    if (!std::isfinite(value)) {
      throw std::runtime_error("Nonfinite API output.");
    }
  }
}

double rms(const std::vector<double>& values)
{
  long double sum = 0.0;
  for (double value : values) {
    sum += static_cast<long double>(value) * value;
  }
  return std::sqrt(static_cast<double>(sum / values.size()));
}

std::vector<double> signature(
  const std::vector<double>& energy,
  const std::vector<double>& force,
  const std::vector<double>& virial,
  const std::vector<double>& charge,
  const std::vector<double>& bec,
  bool charged)
{
  check_finite(energy);
  check_finite(force);
  check_finite(virial);
  const std::size_t n = energy.size();
  long double total = 0.0;
  for (double value : energy) {
    total += value;
  }
  double force_max = 0.0;
  for (double value : force) {
    force_max = std::max(force_max, std::abs(value));
  }
  std::vector<double> result = {static_cast<double>(total / n), rms(force), force_max};
  for (int d = 0; d < 9; ++d) {
    total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      total += virial[d * n + i];
    }
    result.push_back(static_cast<double>(total / n));
  }
  if (charged) {
    check_finite(charge);
    check_finite(bec);
    result.push_back(rms(charge));
    result.push_back(rms(bec));
  }
  return result;
}

void write_array(std::ostream& out, const std::vector<double>& values)
{
  out << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      out << ',';
    }
    out << values[i];
  }
  out << ']';
}
} // namespace

int main(int argc, char** argv)
{
  try {
    if (argc != 6) {
      throw std::runtime_error("Usage: driver MODEL FRAMES OUTPUT REPEATS SECONDS");
    }
    const int repeats = std::stoi(argv[4]);
    const double seconds = std::stod(argv[5]);
    if (repeats < 1 || !std::isfinite(seconds) || seconds <= 0.0) {
      throw std::runtime_error("REPEATS and SECONDS must be positive.");
    }
    std::ifstream input(argv[2]);
    int n = 0;
    int frame_count = 0;
    input >> n >> frame_count;
    if (!input || n <= 0 || frame_count <= 0) {
      throw std::runtime_error("Invalid frame header.");
    }
    std::vector<double> box(9);
    std::vector<int> type(n);
    std::vector<std::vector<double>> position(frame_count, std::vector<double>(n * 3));
    for (double& value : box) {
      input >> value;
    }
    for (int& value : type) {
      input >> value;
    }
    for (auto& frame : position) {
      for (int i = 0; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
          input >> frame[d * n + i];
        }
      }
    }
    if (!input) {
      throw std::runtime_error("Incomplete frame data.");
    }
    int threads = 1;
#if defined(_OPENMP)
#pragma omp parallel
    {
#pragma omp single
      threads = omp_get_num_threads();
    }
#endif
    NEP nep;
    nep.init_from_file(argv[1], false);
    const bool charged = nep.paramb.charge_mode != 0;
    std::vector<double> energy(n), force(n * 3), virial(n * 9), charge(n), bec(n * 9);
    auto compute = [&](int frame) {
      if (charged) {
        nep.compute(type, box, position[frame], energy, force, virial, charge, bec);
      } else {
        nep.compute(type, box, position[frame], energy, force, virial);
      }
    };
    std::vector<std::vector<double>> signatures;
    for (int frame = 0; frame < frame_count; ++frame) {
      compute(frame);
      signatures.push_back(signature(energy, force, virial, charge, bec, charged));
    }
    auto measure = [&](int steps) {
      const auto start = std::chrono::steady_clock::now();
      for (int step = 0; step < steps; ++step) {
        compute(step % frame_count);
      }
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    };
    const double calibration = measure(frame_count);
    if (!std::isfinite(calibration) || calibration <= 0.0) {
      throw std::runtime_error("Invalid calibration time.");
    }
    const double cycles = std::max(1.0, std::ceil(seconds / calibration));
    if (cycles > 10000000.0 / frame_count) {
      throw std::runtime_error("Requested timing interval is too long.");
    }
    const int steps = static_cast<int>(cycles) * frame_count;
    std::vector<double> elapsed;
    for (int repeat = 0; repeat < repeats; ++repeat) {
      elapsed.push_back(measure(steps));
    }
    const auto final_signature = signature(energy, force, virial, charge, bec, charged);
    for (std::size_t i = 0; i < final_signature.size(); ++i) {
      const double expected = signatures.back()[i];
      if (std::abs(final_signature[i] - expected) > 1.0e-8 * (1.0 + std::abs(expected))) {
        throw std::runtime_error("API outputs changed after repeated calls.");
      }
    }
    std::ofstream out(argv[3]);
    out << std::setprecision(17);
    out << "{\"atoms\":" << n << ",\"frames\":" << frame_count << ",\"threads\":" << threads
        << ",\"steps_per_repeat\":" << steps << ",\"elapsed_seconds\":";
    write_array(out, elapsed);
    out << ",\"signatures\":[";
    for (std::size_t i = 0; i < signatures.size(); ++i) {
      if (i) {
        out << ',';
      }
      write_array(out, signatures[i]);
    }
    out << "]}\n";
    out.close();
    if (!out) {
      throw std::runtime_error("Cannot write timing output.");
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
