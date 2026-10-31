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
  if (argc != 4 && argc != 5) {
    std::cerr << "Usage: driver MODEL INPUT OUTPUT [RANK]\n";
    return 2;
  }
  NEP nep;
  nep.init_from_file(argv[1], argc == 4 || std::string(argv[4]) == "0");
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
    const bool lammps = action == "lammps" || action == "lammps_total" ||
                        action == "lammps_energy" || action == "lammps_virial";
    const bool atom_energy = action != "lammps_total" && action != "lammps_virial";
    const bool atom_virial = action != "lammps_total" && action != "lammps_energy";
    double lammps_energy = 0.0;
    std::vector<double> lammps_virial(6, 0.0);
    if (action == "bad_position") pos.pop_back();
    if (action == "bad_energy") energy.pop_back();
    if (action == "bad_force") force.pop_back();
    if (action == "bad_virial") virial.pop_back();
    if (action == "bad_charge") charge.pop_back();
    if (action == "bad_bec") bec.pop_back();

    if (action == "layout") {
      out << "{\"coefficient_offset\":" << (nep.annmb.c - nep.parameters.data());
      write_array(out, "parameters", nep.parameters, 1, nep.parameters.size());
      out << "}\n";
      continue;
    } else if (action == "load") {
      out << "{\"dimension\":" << nep.annmb.dim << "}\n";
      continue;
    } else if (action == "bad_type_map") {
      int mapping[2] = {-1, 0};
      char hydrogen[] = "H";
      char* elements[] = {hydrogen};
      nep.update_type_map(1, mapping, elements);
    } else if (action == "bad_lammps") {
      // Only used with unsupported models, which must reject before any access.
      double total = 0, vv[6] = {};
      nep.compute_for_lammps(0, 0, nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, total, vv, nullptr, nullptr, nullptr);
    } else if (lammps) {
      // Full neighbor list for a cluster with no periodic images. Reverse the
      // LAMMPS type numbering to exercise the element mapping independently.
      std::vector<int> ilist(n), nn(n, n - 1), atom_type(n), mapping(nep.paramb.num_types + 1);
      std::vector<std::vector<int>> neighbors(n);
      std::vector<std::vector<double>> positions(n, std::vector<double>(3));
      std::vector<std::vector<double>> forces(n, std::vector<double>(3, 0.0));
      std::vector<std::vector<double>> virials(n, std::vector<double>(9, 0.0));
      std::vector<int*> nl(n);
      std::vector<double*> x(n), f(n), v(n);
      std::vector<std::string> elements(nep.element_list.rbegin(), nep.element_list.rend());
      std::vector<char*> names(elements.size());
      mapping[0] = -1;
      for (std::size_t t = 0; t < elements.size(); ++t) {
        names[t] = &elements[t][0];
        mapping[t + 1] = t;
      }
      nep.update_type_map(elements.size(), mapping.data(), names.data());
      for (int i = 0; i < n; ++i) {
        ilist[i] = n - 1 - i;
        atom_type[i] = nep.paramb.num_types - type[i];
        for (int j = 0; j < n; ++j) {
          if (i != j) {
            neighbors[i].push_back(j);
          }
        }
        for (int d = 0; d < 3; ++d) {
          positions[i][d] = pos[d * n + i];
        }
        nl[i] = neighbors[i].data();
        x[i] = positions[i].data();
        f[i] = forces[i].data();
        v[i] = virials[i].data();
        energy[i] = 0.0;
      }
      nep.compute_for_lammps(
        n, n, ilist.data(), nn.data(), nl.data(), atom_type.data(), mapping.data(), x.data(),
        lammps_energy, lammps_virial.data(), atom_energy ? energy.data() : nullptr, f.data(),
        atom_virial ? v.data() : nullptr);
      const int order[9] = {0, 3, 4, 6, 1, 5, 7, 8, 2};
      for (int i = 0; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
          force[d * n + i] = forces[i][d];
        }
        for (int d = 0; d < 9; ++d) {
          virial[d * n + i] = virials[i][order[d]];
        }
      }
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
    out << "{\"total_energy\":"
        << (lammps
              ? lammps_energy
              : (response.empty() ? std::accumulate(energy.begin(), energy.end(), 0.0L) : 0.0L));
    if (lammps) {
      write_array(out, "total_virial", lammps_virial, 1, 6);
    }
    if (action != "energy") {
      if (response.empty()) {
        if (atom_energy) {
          write_array(out, "energy", energy, n, 1);
        }
        write_array(out, "force", force, n, 3);
        if (atom_virial) {
          write_array(out, "virial", virial, n, 9);
        }
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
