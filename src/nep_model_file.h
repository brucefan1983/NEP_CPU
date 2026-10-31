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

#pragma once
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

class NEP_File
{
public:
  struct Field {
    std::string name;
    std::vector<std::string> values;
    std::size_t line = 0;
  };

  struct Block {
    std::string name;
    std::vector<double> values;
    std::size_t line = 0;
  };

  std::vector<std::string> elements;
  std::vector<Field> fields;
  std::vector<Block> blocks;

  const Field* find_field(const std::string& name) const;
  const Block* find_block(const std::string& name) const;
  static bool parse_size(const std::string& token, std::size_t& value);
  static bool parse_double(const std::string& token, double& value);
  bool read(std::istream& input, std::string& error);
  bool check_unhandled(
    const std::vector<std::string>& handled_fields,
    const std::vector<std::string>& handled_blocks,
    std::string& error) const;

private:
  bool validate(std::string& error);
};

class NEP_Model_File
{
public:
  struct ModelLayout {
    int dim_radial = 0;
    int dim_angular = 0;
    int dim = 0;
    int ann_count = 0;
    int global_parameters = 0;
    int number_of_variables_ann_1 = 0;
    int number_of_variables_ann = 0;
    int number_of_variables_descriptor = 0;
    int number_of_variables = 0;
  };

  std::string name;
  std::vector<std::string> elements;
  int num_types = 0;
  int num_neurons1 = 0;
  int num_neurons2 = 0;
  int n_max_radial = 0;
  int n_max_angular = 0;
  int basis_size_radial = 0;
  int basis_size_angular = 0;
  int L_max = 0;
  int has_q_222 = 0;
  int has_q_1111 = 0;
  int has_q_112 = 0;
  int has_q_123 = 0;
  int has_q_233 = 0;
  int has_q_134 = 0;
  bool has_multiple_cutoffs = false;
  std::vector<double> rc_radial;
  std::vector<double> rc_angular;
  int max_radial = 0;
  int max_angular = 0;
  bool enable_zbl = false;
  bool flexible_zbl = false;
  bool use_typewise_cutoff_zbl = false;
  double zbl_rc_inner = 0.0;
  double zbl_rc_outer = 0.0;
  double typewise_cutoff_zbl_factor = 0.0;
  ModelLayout layout;
  std::vector<double> energy_reference;
  std::vector<double> parameters;
  std::vector<double> q_scaler;
  std::vector<double> zbl_parameters;

  bool read(std::istream& input, std::string& error);

private:
  bool parse(const NEP_File& document, std::string& error);
  static bool fail(std::string& error, const std::string& message, std::size_t line = 0);
  static bool read_integer(
    const NEP_File::Field& field,
    std::size_t index,
    int minimum,
    int maximum,
    int& value,
    std::string& error);
  static bool
  read_float(const NEP_File::Field& field, std::size_t index, double& value, std::string& error);
  static const NEP_File::Field* require_field(
    const NEP_File& document,
    const char* name,
    std::size_t minimum,
    std::size_t maximum,
    std::vector<std::string>& handled,
    std::string& error);
  static bool read_block(
    const NEP_File& document,
    const char* name,
    int count,
    std::vector<double>& values,
    std::vector<std::string>& handled,
    std::string& error);
  bool
  read_descriptor(const NEP_File& document, std::vector<std::string>& handled, std::string& error);
  bool read_ann(const NEP_File& document, std::vector<std::string>& handled, std::string& error);
  bool
  read_reference(const NEP_File& document, std::vector<std::string>& handled, std::string& error);
  bool read_zbl(
    const NEP_File& document,
    std::vector<std::string>& fields,
    std::vector<std::string>& blocks,
    std::string& error);
  void calculate_layout(int extra_inputs, int extra_outputs, int globals, int ann_count);
  bool
  read_parameters(const NEP_File& document, std::vector<std::string>& blocks, std::string& error);
  bool read_nep(
    const NEP_File& document,
    std::vector<std::string>& fields,
    std::vector<std::string>& blocks,
    std::string& error);
  bool read_charge(
    const NEP_File& document,
    std::vector<std::string>& fields,
    std::vector<std::string>& blocks,
    std::string& error);
  bool read_tnep(
    const NEP_File& document,
    std::vector<std::string>& fields,
    std::vector<std::string>& blocks,
    std::string& error);
};
