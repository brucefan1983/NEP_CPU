// Thin executable adapter for the Python tests. No numerical assertions live here.
#include "nep.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

static void write_array(
  std::ostream& out, const char* name, const std::vector<double>& values,
  int rows, int cols, bool soa = true)
{
  out << ",\"" << name << "\":[";
  for (int i = 0; i < rows; ++i) {
    if (i) out << ',';
    out << '[';
    for (int j = 0; j < cols; ++j) {
      if (j) out << ',';
      out << values[soa ? j * rows + i : i * cols + j];
    }
    out << ']';
  }
  out << ']';
}

int main(int argc, char** argv)
{
  if (argc != 4) {
    std::cerr << "Usage: driver MODEL INPUT OUTPUT\n";
    return 2;
  }
  NEP nep;
  nep.init_from_file(argv[1], true);
  std::ifstream input(argv[2]);
  std::ofstream out(argv[3]);
  if (!input || !out) return 2;
  out << std::setprecision(17);
  std::string action;
  int n;
  while (input >> action >> n) {
    if (n <= 0) return 2;
    std::vector<int> type(n);
    std::vector<double> box(9), pos(n * 3);
    for (double& v : box) input >> v;
    for (int i = 0; i < n; ++i) {
      input >> type[i];
      for (int d = 0; d < 3; ++d) input >> pos[d * n + i];
    }
    if (!input) return 2;
    // Deliberately nonzero: outputs must be overwritten, not accumulated.
    std::vector<double> energy(n, 7.125), force(n * 3, 7.125), virial(n * 9, 7.125);
    std::vector<double> charge(n, 7.125), bec(n * 9, 7.125);
    std::vector<double> response;
    if (action == "bad_position") pos.pop_back();
    if (action == "bad_energy") energy.pop_back();
    if (action == "bad_force") force.pop_back();
    if (action == "bad_virial") virial.pop_back();
    if (action == "bad_charge") charge.pop_back();
    if (action == "bad_bec") bec.pop_back();

    if (action == "load") {
      out << "{\"dimension\":" << nep.annmb.dim << "}\n";
      continue;
    } else if (action == "bad_lammps") {
      // Only used with unsupported models, which must reject before any access.
      double total = 0, vv[6] = {};
      nep.compute_for_lammps(0, 0, nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, total, vv, nullptr, nullptr, nullptr);
    } else if (action == "dipole" || (action == "full" && nep.paramb.model_type == 1)) {
      response.assign(3, 7.125);
      nep.find_dipole(type, box, pos, response);
    } else if (action == "polarizability" || (action == "full" && nep.paramb.model_type == 2)) {
      response.assign(6, 7.125);
      nep.find_polarizability(type, box, pos, response);
    } else if (action == "d3") {
      nep.compute_dftd3("pbe", 8.0, 6.0, type, box, pos, energy, force, virial);
    } else if (action == "with_d3") {
      nep.compute_with_dftd3("pbe", 8.0, 6.0, type, box, pos, energy, force, virial);
    } else if (action == "ordinary" || (!nep.paramb.charge_mode && action != "charge")) {
      nep.compute(type, box, pos, energy, force, virial);
    } else {
      nep.compute(type, box, pos, energy, force, virial, charge, bec);
    }
    // Extra accumulation precision reduces cancellation in coordinate finite
    // differences. The API's per-atom energies remain unchanged doubles.
    out << "{\"total_energy\":" << (response.empty() ?
      std::accumulate(energy.begin(), energy.end(), 0.0L) : 0.0L);
    if (action != "energy") {
      if (response.empty()) {
        write_array(out, "energy", energy, n, 1);
        write_array(out, "force", force, n, 3);
        write_array(out, "virial", virial, n, 9);
        if (nep.paramb.charge_mode) {
          write_array(out, "charge", charge, n, 1);
          write_array(out, "bec", bec, n, 9);
        }
      } else {
        write_array(out, "response", response, 1, response.size());
      }
      if (action == "full" || action == "features") {
        std::vector<double> descriptor(n * nep.annmb.dim, 7.125);
        nep.find_descriptor(type, box, pos, descriptor);
        write_array(out, "descriptor", descriptor, n, nep.annmb.dim);
      }
      if (action == "features") {
        std::vector<double> latent(n * nep.annmb.num_neurons1, 7.125);
        nep.find_latent_space(type, box, pos, latent);
        write_array(out, "latent", latent, n, nep.annmb.num_neurons1);
        const int width = nep.annmb.num_neurons1 * (nep.annmb.dim + 2);
        std::vector<double> projection(n * width, 7.125);
        nep.find_B_projection(type, box, pos, projection);
        write_array(out, "projection", projection, n, width, false);
      }
    }
    out << "}\n";
  }
  return input.eof() && out ? 0 : 2;
}
