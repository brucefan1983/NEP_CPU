/*
    Copyright 2025 Zheyong Fan, Junjie Wang, Eric Lindgren, Hekai Bu
    NEP_CPU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with NEP_CPU.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include "nep.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace NEP_UTILITIES
{

struct Atom {
  int N;
  std::vector<int> type;
  std::vector<double> box, position, potential, force, virial;
};

struct Structure {
  int num_atom;
  int has_virial;
  double energy;
  double weight;
  double virial[9];
  double box[9];
  std::vector<std::string> atom_symbol;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> fx;
  std::vector<double> fy;
  std::vector<double> fz;
};

void readXYZ(Atom& atom);
void timing(Atom& atom, NEP_NS::NEP3& nep3);
void compare_analytical_and_finite_difference(Atom& atom, NEP_NS::NEP3& nep3);
void get_descriptor(Atom& atom, NEP_NS::NEP3& nep3);
std::vector<std::string> get_atom_symbols(void);

std::string remove_spaces_step1(const std::string& line);
std::string remove_spaces(const std::string& line_input);
std::vector<std::string> get_tokens(const std::string& line);
std::vector<std::string> get_tokens_without_unwanted_spaces(std::ifstream& input);
std::vector<std::string> get_tokens(std::ifstream& input);
int get_int_from_token(const std::string& token, const char* filename, const int line);
double get_double_from_token(const std::string& token, const char* filename, const int line);

void read_force(
  const int num_columns,
  const int species_offset,
  const int pos_offset,
  const int force_offset,
  std::ifstream& input,
  Structure& structure);

void read_one_structure(std::ifstream& input, Structure& structure);
void read(const std::string& inputfile, std::vector<Structure>& structures);
void write_one_structure(std::ofstream& output, const Structure& structure);
void write(const std::string& outputfile, const std::vector<Structure>& structures);
std::vector<std::string> get_atom_symbols(std::string& nep_file);
void calculate_one_structure(
  NEP_NS::NEP3& nep3,
  std::vector<std::string>& atom_symbols,
  Structure& structure,
  int mode,
  std::string& functional,
  double D3_cutoff,
  double D3_cutoff_cn);
}