/*
    Copyright 2022 Zheyong Fan, Junjie Wang, Eric Lindgren, Hekai Bu
    This file is part of NEP_CPU.
    NEP_CPU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    NEP_CPU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with NEP_CPU.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
Usage:
    Compile:
        g++ -O3 main.cpp ../src/nep.cpp # Without openMP support
        g++ -O3 -fopenmp main.cpp ../src/nep.cpp # With openMP support
    run:
        export OMP_NUM_THREADS=6 # 6 is the number of the threads to be used
        ./a.out
------------------------------------------------------------------------------*/

#include "../include/nep_utilities.h"

const int num_repeats = 1000;

namespace NEP_UTILITIES
{
struct Atom {
  int N;
  std::vector<int> type;
  std::vector<double> box, position, potential, force, virial;
};


std::vector<std::string> get_atom_symbols()
{
  std::ifstream input_potential("nep.txt");
  if (!input_potential.is_open()) {
    std::cout << "Error: cannot open nep.txt.\n";
    exit(1);
  }

  std::string potential_name;
  input_potential >> potential_name;
  int number_of_types;
  input_potential >> number_of_types;
  std::vector<std::string> atom_symbols(number_of_types);
  for (int n = 0; n < number_of_types; ++n) {
    input_potential >> atom_symbols[n];
  }

  input_potential.close();
  return atom_symbols;
}

void readXYZ(Atom& atom)
{
  std::cout << "Reading xyz.in.\n";

  std::ifstream input_file("xyz.in");

  if (!input_file) {
    std::cout << "Cannot open xyz.in\n";
    exit(1);
  }

  input_file >> atom.N;
  std::cout << "    Number of atoms is " << atom.N << ".\n";

  atom.box.resize(9);
  input_file >> atom.box[0];
  input_file >> atom.box[3];
  input_file >> atom.box[6];
  input_file >> atom.box[1];
  input_file >> atom.box[4];
  input_file >> atom.box[7];
  input_file >> atom.box[2];
  input_file >> atom.box[5];
  input_file >> atom.box[8];

  std::cout << "    Box matrix h = [a, b, c] is\n";
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = 0; d2 < 3; ++d2) {
      std::cout << "\t" << atom.box[d1 * 3 + d2];
    }
    std::cout << "\n";
  }

  std::vector<std::string> atom_symbols = get_atom_symbols();

  atom.type.resize(atom.N);
  atom.position.resize(atom.N * 3);
  atom.potential.resize(atom.N);
  atom.force.resize(atom.N * 3);
  atom.virial.resize(atom.N * 9);

  for (int n = 0; n < atom.N; n++) {
    std::string atom_symbol_tmp;
    input_file >> atom_symbol_tmp >> atom.position[n] >> atom.position[n + atom.N] >>
      atom.position[n + atom.N * 2];
    bool is_allowed_element = false;
    for (int t = 0; t < atom_symbols.size(); ++t) {
      if (atom_symbol_tmp == atom_symbols[t]) {
        atom.type[n] = t;
        is_allowed_element = true;
      }
    }
    if (!is_allowed_element) {
      std::cout << "There is atom in xyz.in that is not allowed in the used NEP potential.\n";
      exit(1);
    }
  }
}

void timing(Atom& atom, NEP_NS::NEP3& nep3)
{
  std::cout << "Started timing.\n";

  clock_t time_begin = clock();

  for (int n = 0; n < num_repeats; ++n) {
    nep3.compute(atom.type, atom.box, atom.position, atom.potential, atom.force, atom.virial);
  }

  clock_t time_finish = clock();
  double time_used = (time_finish - time_begin) / double(CLOCKS_PER_SEC);
  std::cout << "    Number of atoms = " << atom.N << ".\n";
  std::cout << "    Number of steps = " << num_repeats << ".\n";
  std::cout << "    Time used = " << time_used << " s.\n";

  double speed = atom.N * num_repeats / time_used;
  double cost = 1000 / speed;
  std::cout << "    Computational speed = " << speed << " atom-step/second.\n";
  std::cout << "    Computational cost = " << cost << " mini-second/atom-step.\n";
}

void compare_analytical_and_finite_difference(Atom& atom, NEP_NS::NEP3& nep3)
{
  std::cout << "Started validating force.\n";

  std::vector<double> force_finite_difference(atom.force.size());
  std::vector<double> position_copy(atom.position.size());
  for (int n = 0; n < atom.position.size(); ++n) {
    position_copy[n] = atom.position[n];
  }

  const double delta = 2.0e-5;

  for (int n = 0; n < atom.N; ++n) {
    for (int d = 0; d < 3; ++d) {
      atom.position[n + d * atom.N] = position_copy[n + d * atom.N] - delta; // negative shift

      nep3.compute(atom.type, atom.box, atom.position, atom.potential, atom.force, atom.virial);

      double energy_negative_shift = 0.0;
      for (int n = 0; n < atom.N; ++n) {
        energy_negative_shift += atom.potential[n];
      }

      atom.position[n + d * atom.N] = position_copy[n + d * atom.N] + delta; // positive shift

      nep3.compute(atom.type, atom.box, atom.position, atom.potential, atom.force, atom.virial);

      double energy_positive_shift = 0.0;
      for (int n = 0; n < atom.N; ++n) {
        energy_positive_shift += atom.potential[n];
      }

      force_finite_difference[n + d * atom.N] =
        (energy_negative_shift - energy_positive_shift) / (2.0 * delta);

      atom.position[n + d * atom.N] = position_copy[n + d * atom.N]; // back to original position
    }
  }

  nep3.compute(atom.type, atom.box, atom.position, atom.potential, atom.force, atom.virial);

  std::ofstream output_file("force_analytical.out");

  if (!output_file.is_open()) {
    std::cout << "Cannot open force_analytical.out\n";
    exit(1);
  }
  output_file << std::setprecision(15);
  for (int n = 0; n < atom.N; ++n) {
    output_file << atom.force[n] << " " << atom.force[n + atom.N] << " "
                << atom.force[n + atom.N * 2] << "\n";
  }
  output_file.close();
  std::cout << "    analytical forces are written into force_analytical.out.\n";

  std::ofstream output_finite_difference("force_finite_difference.out");

  if (!output_finite_difference.is_open()) {
    std::cout << "Cannot open force_finite_difference.out\n";
    exit(1);
  }
  output_finite_difference << std::setprecision(15);
  for (int n = 0; n < atom.N; ++n) {
    output_finite_difference << force_finite_difference[n] << " "
                             << force_finite_difference[n + atom.N] << " "
                             << force_finite_difference[n + atom.N * 2] << "\n";
  }
  output_finite_difference.close();
  std::cout << "    finite-difference forces are written into force_finite_difference.out.\n";

  std::ofstream output_file_virial("virial.out");

  if (!output_file_virial.is_open()) {
    std::cout << "Cannot open virial.out\n";
    exit(1);
  }
  output_file_virial << std::setprecision(15);
  for (int n = 0; n < atom.N; ++n) {
    output_file_virial << atom.virial[n + atom.N * 0] << " " << atom.virial[n + atom.N * 1] << " "
                       << atom.virial[n + atom.N * 2] << " " << atom.virial[n + atom.N * 3] << " "
                       << atom.virial[n + atom.N * 4] << " " << atom.virial[n + atom.N * 5] << " "
                       << atom.virial[n + atom.N * 6] << " " << atom.virial[n + atom.N * 7] << " "
                       << atom.virial[n + atom.N * 8] << "\n";
  }
  output_file_virial.close();
  std::cout << "    virials are written into virial.out.\n";
}

void get_descriptor(Atom& atom, NEP_NS::NEP3& nep3)
{
  std::cout << "Getting descriptor.\n";

  std::vector<double> descriptor(atom.N * nep3.annmb.dim);

  nep3.find_descriptor(atom.type, atom.box, atom.position, descriptor);

  std::ofstream output_file("descriptor.out");

  if (!output_file.is_open()) {
    std::cout << "Cannot open descriptor.out\n";
    exit(1);
  }

  output_file << std::setprecision(15);

  for (int n = 0; n < atom.N; ++n) {
    for (int d = 0; d < nep3.annmb.dim; ++d) {
      output_file << descriptor[d * atom.N + n] << " ";
    }
    output_file << "\n";
  }
  output_file.close();
  std::cout << "    descriptors are written into descriptor.out.\n";
}

std::string remove_spaces_step1(const std::string& line)
{
  std::vector<int> indices_for_spaces(line.size(), 0);
  for (int n = 0; n < line.size(); ++n) {
    if (line[n] == '=') {
      for (int k = 1; n - k >= 0; ++k) {
        if (line[n - k] == ' ' || line[n - k] == '\t') {
          indices_for_spaces[n - k] = 1;
        } else {
          break;
        }
      }
      for (int k = 1; n + k < line.size(); ++k) {
        if (line[n + k] == ' ' || line[n + k] == '\t') {
          indices_for_spaces[n + k] = 1;
        } else {
          break;
        }
      }
    }
  }

  std::string new_line;
  for (int n = 0; n < line.size(); ++n) {
    if (!indices_for_spaces[n]) {
      new_line += line[n];
    }
  }

  return new_line;
}

std::string remove_spaces(const std::string& line_input)
{
  auto line = remove_spaces_step1(line_input);

  std::vector<int> indices_for_spaces(line.size(), 0);
  for (int n = 0; n < line.size(); ++n) {
    if (line[n] == '\"') {
      if (n == 0) {
        std::cout << "The second line of the .xyz file should not begin with \"." << std::endl;
        exit(1);
      } else {
        if (line[n - 1] == '=') {
          for (int k = 1; n + k < line.size(); ++k) {
            if (line[n + k] == ' ' || line[n + k] == '\t') {
              indices_for_spaces[n + k] = 1;
            } else {
              break;
            }
          }
        } else {
          for (int k = 1; n - k >= 0; ++k) {
            if (line[n - k] == ' ' || line[n - k] == '\t') {
              indices_for_spaces[n - k] = 1;
            } else {
              break;
            }
          }
        }
      }
    }
  }

  std::string new_line;
  for (int n = 0; n < line.size(); ++n) {
    if (!indices_for_spaces[n]) {
      new_line += line[n];
    }
  }

  return new_line;
}

std::vector<std::string> get_tokens(const std::string& line)
{
  std::istringstream iss(line);
  std::vector<std::string> tokens{
    std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
  return tokens;
}

std::vector<std::string> get_tokens_without_unwanted_spaces(std::ifstream& input)
{
  std::string line;
  std::getline(input, line);
  auto line_without_unwanted_spaces = remove_spaces(line);
  std::istringstream iss(line_without_unwanted_spaces);
  std::vector<std::string> tokens{
    std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
  return tokens;
}

std::vector<std::string> get_tokens(std::ifstream& input)
{
  std::string line;
  std::getline(input, line);
  std::istringstream iss(line);
  std::vector<std::string> tokens{
    std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
  return tokens;
}

int get_int_from_token(const std::string& token, const char* filename, const int line)
{
  int value = 0;
  try {
    value = std::stoi(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:\n";
    std::cout << "    File:          " << filename << std::endl;
    std::cout << "    Line:          " << line << std::endl;
    std::cout << "    Error message: " << e.what() << std::endl;
    exit(1);
  }
  return value;
}

double get_double_from_token(const std::string& token, const char* filename, const int line)
{
  double value = 0;
  try {
    value = std::stod(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:\n";
    std::cout << "    File:          " << filename << std::endl;
    std::cout << "    Line:          " << line << std::endl;
    std::cout << "    Error message: " << e.what() << std::endl;
    exit(1);
  }
  return value;
}



void read_force(
  const int num_columns,
  const int species_offset,
  const int pos_offset,
  const int force_offset,
  std::ifstream& input,
  Structure& structure)
{
  structure.atom_symbol.resize(structure.num_atom);
  structure.x.resize(structure.num_atom);
  structure.y.resize(structure.num_atom);
  structure.z.resize(structure.num_atom);
  structure.fx.resize(structure.num_atom);
  structure.fy.resize(structure.num_atom);
  structure.fz.resize(structure.num_atom);

  for (int na = 0; na < structure.num_atom; ++na) {
    std::vector<std::string> tokens = get_tokens(input);
    if (tokens.size() != num_columns) {
      std::cout << "Number of items for an atom line mismatches properties." << std::endl;
      exit(1);
    }
    structure.atom_symbol[na] = tokens[0 + species_offset];
    structure.x[na] = get_double_from_token(tokens[0 + pos_offset], __FILE__, __LINE__);
    structure.y[na] = get_double_from_token(tokens[1 + pos_offset], __FILE__, __LINE__);
    structure.z[na] = get_double_from_token(tokens[2 + pos_offset], __FILE__, __LINE__);
    if (num_columns > 4) {
      structure.fx[na] = get_double_from_token(tokens[0 + force_offset], __FILE__, __LINE__);
      structure.fy[na] = get_double_from_token(tokens[1 + force_offset], __FILE__, __LINE__);
      structure.fz[na] = get_double_from_token(tokens[2 + force_offset], __FILE__, __LINE__);
    }
  }
}

void read_one_structure(std::ifstream& input, Structure& structure)
{
  std::vector<std::string> tokens = get_tokens_without_unwanted_spaces(input);
  for (auto& token : tokens) {
    std::transform(
      token.begin(), token.end(), token.begin(), [](unsigned char c) { return std::tolower(c); });
  }

  if (tokens.size() == 0) {
    std::cout << "The second line for each frame should not be empty." << std::endl;
    exit(1);
  }

  bool has_energy_in_exyz = false;
  for (const auto& token : tokens) {
    const std::string energy_string = "energy=";
    if (token.substr(0, energy_string.length()) == energy_string) {
      has_energy_in_exyz = true;
      structure.energy = get_double_from_token(
        token.substr(energy_string.length(), token.length()), __FILE__, __LINE__);
    }
  }
  if (!has_energy_in_exyz) {
    std::cout << "'energy' is missing in the second line of a frame." << std::endl;
    exit(1);
  }

  structure.weight = 1.0f;
  for (const auto& token : tokens) {
    const std::string weight_string = "weight=";
    if (token.substr(0, weight_string.length()) == weight_string) {
      structure.weight = get_double_from_token(
        token.substr(weight_string.length(), token.length()), __FILE__, __LINE__);
      if (structure.weight <= 0.0f || structure.weight > 100.0f) {
        std::cout << "Configuration weight should > 0 and <= 100." << std::endl;
        exit(1);
      }
    }
  }

  bool has_lattice_in_exyz = false;
  for (int n = 0; n < tokens.size(); ++n) {
    const std::string lattice_string = "lattice=";
    if (tokens[n].substr(0, lattice_string.length()) == lattice_string) {
      has_lattice_in_exyz = true;
      for (int m = 0; m < 9; ++m) {
        structure.box[m] = get_double_from_token(
          tokens[n + m].substr(
            (m == 0) ? (lattice_string.length() + 1) : 0,
            (m == 8) ? (tokens[n + m].length() - 1) : tokens[n + m].length()),
          __FILE__, __LINE__);
      }
    }
  }
  if (!has_lattice_in_exyz) {
    std::cout << "'lattice' is missing in the second line of a frame." << std::endl;
    exit(1);
  }

  structure.has_virial = false;
  for (int n = 0; n < tokens.size(); ++n) {
    const std::string virial_string = "virial=";
    if (tokens[n].substr(0, virial_string.length()) == virial_string) {
      structure.has_virial = true;
      for (int m = 0; m < 9; ++m) {
        structure.virial[m] = get_double_from_token(
          tokens[n + m].substr(
            (m == 0) ? (virial_string.length() + 1) : 0,
            (m == 8) ? (tokens[n + m].length() - 1) : tokens[n + m].length()),
          __FILE__, __LINE__);
      }
    }
  }

  int species_offset = 0;
  int pos_offset = 0;
  int force_offset = 0;
  int num_columns = 0;
  for (int n = 0; n < tokens.size(); ++n) {
    const std::string properties_string = "properties=";
    if (tokens[n].substr(0, properties_string.length()) == properties_string) {
      std::string line = tokens[n].substr(properties_string.length(), tokens[n].length());
      for (auto& letter : line) {
        if (letter == ':') {
          letter = ' ';
        }
      }
      std::vector<std::string> sub_tokens = get_tokens(line);
      int species_position = -1;
      int pos_position = -1;
      int force_position = -1;
      for (int k = 0; k < sub_tokens.size() / 3; ++k) {
        if (sub_tokens[k * 3] == "species") {
          species_position = k;
        }
        if (sub_tokens[k * 3] == "pos") {
          pos_position = k;
        }
        if (sub_tokens[k * 3] == "force" || sub_tokens[k * 3] == "forces") {
          force_position = k;
        }
      }
      if (species_position < 0) {
        std::cout << "'species' is missing in properties." << std::endl;
        exit(1);
      }
      if (pos_position < 0) {
        std::cout << "'pos' is missing in properties." << std::endl;
        exit(1);
      }
      if (force_position < 0) {
        std::cout << "'force' or 'forces' is missing in properties." << std::endl;
        exit(1);
      }
      for (int k = 0; k < sub_tokens.size() / 3; ++k) {
        if (k < species_position) {
          species_offset += get_int_from_token(sub_tokens[k * 3 + 2], __FILE__, __LINE__);
        }
        if (k < pos_position) {
          pos_offset += get_int_from_token(sub_tokens[k * 3 + 2], __FILE__, __LINE__);
        }
        if (k < force_position) {
          force_offset += get_int_from_token(sub_tokens[k * 3 + 2], __FILE__, __LINE__);
        }
        num_columns += get_int_from_token(sub_tokens[k * 3 + 2], __FILE__, __LINE__);
      }
    }
  }

  read_force(num_columns, species_offset, pos_offset, force_offset, input, structure);
}

void read(const std::string& inputfile, std::vector<Structure>& structures)
{
  std::ifstream input(inputfile);
  if (!input.is_open()) {
    std::cout << "Failed to open " << inputfile << std::endl;
    exit(1);
  } else {
    while (true) {
      std::vector<std::string> tokens = get_tokens(input);
      if (tokens.size() == 0) {
        break;
      } else if (tokens.size() > 1) {
        std::cout << "The first line for each frame should have one value." << std::endl;
        exit(1);
      }
      Structure structure;
      structure.num_atom = get_int_from_token(tokens[0], __FILE__, __LINE__);
      if (structure.num_atom < 1) {
        std::cout << "Number of atoms for each frame should >= 1." << std::endl;
        exit(1);
      }
      read_one_structure(input, structure);
      structures.emplace_back(structure);
    }
    input.close();
  }
}

void write_one_structure(std::ofstream& output, const Structure& structure)
{
  output << structure.num_atom << "\n";
  output << "lattice=\"";
  for (int m = 0; m < 9; ++m) {
    output << structure.box[m];
    if (m != 8) {
      output << " ";
    }
  }
  output << "\" ";

  output << "energy=" << structure.energy << " ";

  if (structure.has_virial) {
    output << "virial=\"";
    for (int m = 0; m < 9; ++m) {
      output << structure.virial[m];
      if (m != 8) {
        output << " ";
      }
    }
    output << "\" ";
  }

  output << "Properties=species:S:1:pos:R:3:force:R:3\n";

  for (int n = 0; n < structure.num_atom; ++n) {
    output << structure.atom_symbol[n] << " " << structure.x[n] << " " << structure.y[n] << " "
           << structure.z[n] << " " << structure.fx[n] << " " << structure.fy[n] << " "
           << structure.fz[n] << "\n";
  }
}

void write(const std::string& outputfile, const std::vector<Structure>& structures)
{
  std::ofstream output(outputfile);
  if (!output.is_open()) {
    std::cout << "Failed to open " << outputfile << std::endl;
    exit(1);
  }
  for (int nc = 0; nc < structures.size(); ++nc) {
    write_one_structure(output, structures[nc]);
  }
  output.close();
}

std::vector<std::string> get_atom_symbols(std::string& nep_file)
{
  std::ifstream input_potential(nep_file);
  if (!input_potential.is_open()) {
    std::cout << "Failed to open " << nep_file << std::endl;
    exit(1);
  }

  std::string potential_name;
  input_potential >> potential_name;
  int number_of_types;
  input_potential >> number_of_types;
  std::vector<std::string> atom_symbols(number_of_types);
  for (int n = 0; n < number_of_types; ++n) {
    input_potential >> atom_symbols[n];
  }

  input_potential.close();
  return atom_symbols;
}

void calculate_one_structure(
  NEP_NS::NEP3& nep3,
  std::vector<std::string>& atom_symbols,
  Structure& structure,
  int mode,
  std::string& functional,
  double D3_cutoff,
  double D3_cutoff_cn)
{
  std::vector<double> box(9);
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = 0; d2 < 3; ++d2) {
      box[d1 * 3 + d2] = structure.box[d2 * 3 + d1];
    }
  }

  std::vector<int> type(structure.num_atom);
  std::vector<double> position(structure.num_atom * 3);
  std::vector<double> potential(structure.num_atom);
  std::vector<double> force(structure.num_atom * 3);
  std::vector<double> virial(structure.num_atom * 9);

  for (int n = 0; n < structure.num_atom; n++) {
    position[n] = structure.x[n];
    position[n + structure.num_atom] = structure.y[n];
    position[n + structure.num_atom * 2] = structure.z[n];

    bool is_allowed_element = false;
    for (int t = 0; t < atom_symbols.size(); ++t) {
      if (structure.atom_symbol[n] == atom_symbols[t]) {
        type[n] = t;
        is_allowed_element = true;
      }
    }
    if (!is_allowed_element) {
      std::cout << "There is atom not allowed in the used NEP potential.\n";
      exit(1);
    }
  }

  if (mode == 0) {
    nep3.compute(type, box, position, potential, force, virial);
  } else if (mode == 1 || mode == 3 || mode == 4) {
    nep3.compute_dftd3(
      functional, D3_cutoff, D3_cutoff_cn, type, box, position, potential, force, virial);
  } else if (mode == 2) {
    nep3.compute_with_dftd3(
      functional, D3_cutoff, D3_cutoff_cn, type, box, position, potential, force, virial);
  }

  if (mode < 3) {
    structure.energy = 0.0;
    for (int n = 0; n < structure.num_atom; n++) {
      structure.energy += potential[n];
      structure.fx[n] = force[0 * structure.num_atom + n];
      structure.fy[n] = force[1 * structure.num_atom + n];
      structure.fz[n] = force[2 * structure.num_atom + n];
    }
    for (int d = 0; d < 9; ++d) {
      structure.virial[d] = 0.0;
      for (int n = 0; n < structure.num_atom; n++) {
        structure.virial[d] += virial[d * structure.num_atom + n];
      }
    }
  } else {
    double factor = (mode == 3) ? 1.0 : -1.0;
    for (int n = 0; n < structure.num_atom; n++) {
      structure.energy += factor * potential[n];
      structure.fx[n] += factor * force[0 * structure.num_atom + n];
      structure.fy[n] += factor * force[1 * structure.num_atom + n];
      structure.fz[n] += factor * force[2 * structure.num_atom + n];
    }
    for (int d = 0; d < 9; ++d) {
      for (int n = 0; n < structure.num_atom; n++) {
        structure.virial[d] += factor * virial[d * structure.num_atom + n];
      }
    }
  }
}

void calculate(
  std::string& nep_file,
  std::vector<Structure>& structures,
  int mode,
  std::string& functional,
  double D3_cutoff,
  double D3_cutoff_cn)
{
  NEP_NS::NEP3 nep3(nep_file);
  std::vector<std::string> atom_symbols = get_atom_symbols(nep_file);
  for (int nc = 0; nc < structures.size(); ++nc) {
    calculate_one_structure(
      nep3, atom_symbols, structures[nc], mode, functional, D3_cutoff, D3_cutoff_cn);
  }
}

int main_test(int argc, char* argv[])
{
  if (argc < 5) {
    std::cout << "Usage:\n";
    std::cout
      << argv[0]
      << " nep_file input_xyz_file output_xyz_file mode [functional D3_cutoff D3_cutoff_cn]\n";
    exit(1);
  } else {
    std::string nep_file = argv[1];
    std::string input_file = argv[2];
    std::string output_file = argv[3];
    int mode = atoi(argv[4]);
    std::string functional;
    double D3_cutoff = 0.0;
    double D3_cutoff_cn = 0.0;
    if (mode < 0 || mode > 4) {
      std::cout << "mode can only be 0, 1, 2, 3, 4.\n";
      exit(1);
    }
    if (mode > 0) {
      if (argc != 8) {
        std::cout << "Usage:\n";
        std::cout
          << argv[0]
          << " nep_file input_xyz_file output_xyz_file mode functional D3_cutoff D3_cutoff_cn\n";
        exit(1);
      } else {
        functional = argv[5];
        D3_cutoff = atof(argv[6]);
        D3_cutoff_cn = atof(argv[7]);
      }
    }
    std::vector<Structure> structures;
    read(input_file, structures);
    std::cout << "Read " << structures.size() << " structures from " << input_file << std::endl;

    if (mode == 0) {
      std::cout << "Calculate NEP only." << std::endl;
    } else if (mode == 1) {
      std::cout << "Calculate DFT-D3 only." << std::endl;
    } else if (mode == 2) {
      std::cout << "Calculate NEP in combination with DFT-D3." << std::endl;
    } else if (mode == 3) {
      std::cout << "Add DFT-D3 to existing data." << std::endl;
    } else if (mode == 4) {
      std::cout << "Subtract DFT-D3 from existing data." << std::endl;
    }

    clock_t time_begin = clock();
    calculate(nep_file, structures, mode, functional, D3_cutoff, D3_cutoff_cn);
    clock_t time_finish = clock();
    double time_used = (time_finish - time_begin) / double(CLOCKS_PER_SEC);
    std::cout << "Time used = " << time_used << " s.\n";

    write(output_file, structures);
  }

  std::cout << "Done." << std::endl;
  return EXIT_SUCCESS;
}

}