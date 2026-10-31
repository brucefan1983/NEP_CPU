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

#include "nep_model_file.h"
#include <cmath>
#include <istream>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <utility>

static std::string lower(std::string text)
{
  for (char& c : text) {
    if (c >= 'A' && c <= 'Z') {
      c += 'a' - 'A';
    }
  }
  return text;
}

static bool identifier(const std::string& text)
{
  if (text.empty()) {
    return false;
  }
  for (std::size_t i = 0; i < text.size(); ++i) {
    const char c = text[i];
    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_' ||
          (i > 0 && c >= '0' && c <= '9'))) {
      return false;
    }
  }
  return true;
}

static bool token_valid(const std::string& token)
{
  if (token.empty()) {
    return false;
  }
  for (const unsigned char c : token) {
    if (c <= 32 || c == 127 || c == '#') {
      return false;
    }
  }
  return true;
}

static bool fail(std::string& error, const std::string& message, std::size_t line = 0)
{
  error = "NEP file" + (line ? " line " + std::to_string(line) : "") + ": " + message;
  return false;
}

static std::string element_symbol(const std::string& token)
{
  static const char* const symbols[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu"};
  const std::string name = lower(token);
  for (const char* symbol : symbols) {
    if (lower(symbol) == name) {
      return symbol;
    }
  }
  return "";
}

static int next_line(std::istream& input, std::string& line)
{
  line.clear();
  try {
    std::getline(input, line);
  } catch (const std::ios_base::failure&) {
    if (!input.eof() || input.bad()) {
      return -1;
    }
  } catch (...) {
    if (!input.bad()) {
      throw;
    }
    return -1;
  }
  if (input.bad() || (input.fail() && !input.eof())) {
    return -1;
  }
  return input.fail() ? 0 : 1;
}

static bool tokenize(std::string line, std::vector<std::string>& tokens)
{
  tokens.clear();
  if (!line.empty() && line.back() == '\r') {
    line.pop_back();
  }
  const std::size_t comment = line.find('#');
  if (comment != std::string::npos) {
    line.resize(comment);
  }
  for (const unsigned char c : line) {
    if ((c < 32 && c != '\t') || c == 127) {
      return false;
    }
  }
  std::istringstream input(line);
  input.imbue(std::locale::classic());
  std::string token;
  while (input >> token) {
    tokens.push_back(token);
  }
  return true;
}

bool NEP_File::validate(std::string& error)
{
  if (elements.empty()) {
    return fail(error, "empty element list.", 1);
  }
  std::set<std::string> seen_elements;
  for (std::string& element : elements) {
    const std::string symbol = element_symbol(element);
    if (symbol.empty()) {
      return fail(error, "unknown element '" + element + "'.", 1);
    }
    if (!seen_elements.insert(symbol).second) {
      return fail(error, "duplicate element '" + symbol + "'.", 1);
    }
    element = symbol;
  }
  std::set<std::string> names;
  for (Field& field : fields) {
    field.name = lower(field.name);
    if (!identifier(field.name)) {
      return fail(error, "invalid field name.", field.line);
    }
    if (!names.insert(field.name).second) {
      return fail(error, "duplicate record '" + field.name + "'.", field.line);
    }
    if (field.values.empty()) {
      return fail(error, "empty field '" + field.name + "'.", field.line);
    }
    for (const std::string& value : field.values) {
      if (!token_valid(value)) {
        return fail(error, "invalid token in '" + field.name + "'.", field.line);
      }
    }
    if (field.name == "model") {
      if (field.values.size() != 1 || !identifier(field.values[0])) {
        return fail(error, "model requires one identifier.", field.line);
      }
      field.values[0] = lower(field.values[0]);
    }
  }
  for (Block& block : blocks) {
    block.name = lower(block.name);
    if (!identifier(block.name)) {
      return fail(error, "invalid block name.", block.line);
    }
    if (!names.insert(block.name).second) {
      return fail(error, "duplicate record '" + block.name + "'.", block.line);
    }
    for (const double value : block.values) {
      if (!std::isfinite(value)) {
        return fail(error, "non-finite value in '" + block.name + "'.", block.line);
      }
    }
  }
  const Field* version = find_field("format_version");
  if (version == nullptr) {
    return fail(error, "missing format_version.");
  }
  if (version->values.size() != 1 || version->values[0] != "1") {
    return fail(error, "unsupported format_version (expected 1).", version->line);
  }
  if (find_field("model") == nullptr) {
    return fail(error, "missing model.");
  }
  return true;
}

static bool contains(const std::vector<std::string>& names, const std::string& name)
{
  for (const std::string& entry : names) {
    if (lower(entry) == lower(name)) {
      return true;
    }
  }
  return false;
}

const NEP_File::Field* NEP_File::find_field(const std::string& name) const
{
  for (const Field& field : fields) {
    if (lower(field.name) == lower(name)) {
      return &field;
    }
  }
  return nullptr;
}

const NEP_File::Block* NEP_File::find_block(const std::string& name) const
{
  for (const Block& block : blocks) {
    if (lower(block.name) == lower(name)) {
      return &block;
    }
  }
  return nullptr;
}

bool NEP_File::parse_size(const std::string& token, std::size_t& value)
{
  if (token.empty()) {
    return false;
  }
  std::size_t parsed = 0;
  for (const char c : token) {
    if (c < '0' || c > '9') {
      return false;
    }
    const unsigned int digit = static_cast<unsigned int>(c - '0');
    if (parsed > (std::numeric_limits<std::size_t>::max() - digit) / 10) {
      return false;
    }
    parsed = parsed * 10 + digit;
  }
  value = parsed;
  return true;
}

bool NEP_File::parse_double(const std::string& token, double& value)
{
  std::size_t i = 0;
  if (i < token.size() && (token[i] == '+' || token[i] == '-')) {
    ++i;
  }
  bool digits = false;
  bool nonzero = false;
  bool point = false;
  for (; i < token.size(); ++i) {
    const char c = token[i];
    if (c >= '0' && c <= '9') {
      digits = true;
      nonzero = nonzero || c != '0';
    } else if (c == '.' && !point) {
      point = true;
    } else {
      break;
    }
  }
  if (!digits) {
    return false;
  }
  if (i < token.size() && (token[i] == 'e' || token[i] == 'E')) {
    ++i;
    if (i < token.size() && (token[i] == '+' || token[i] == '-')) {
      ++i;
    }
    const std::size_t start = i;
    while (i < token.size() && token[i] >= '0' && token[i] <= '9') {
      ++i;
    }
    if (i == start) {
      return false;
    }
  }
  if (i != token.size()) {
    return false;
  }
  std::istringstream input(token);
  input.imbue(std::locale::classic());
  double parsed;
  if (!(input >> parsed) || !std::isfinite(parsed) || (parsed == 0.0 && nonzero)) {
    return false;
  }
  value = parsed;
  return true;
}

bool NEP_File::read(std::istream& input, std::string& error)
{
  *this = NEP_File();
  error.clear();
  NEP_File parsed;
  std::string line;
  std::vector<std::string> tokens;
  int status = next_line(input, line);
  if (status < 0) {
    return fail(error, "input stream failure.", 1);
  }
  if (status == 0 || !tokenize(line, tokens) || tokens.size() < 3 || lower(tokens[0]) != "nep") {
    return fail(error, "expected first line: nep <num_types> <elements...>.", 1);
  }
  std::size_t num_types;
  if (!parse_size(tokens[1], num_types) || num_types == 0 || num_types != tokens.size() - 2) {
    return fail(error, "element count does not match first line.", 1);
  }
  parsed.elements.assign(tokens.begin() + 2, tokens.end());
  std::size_t line_number = 1;
  std::size_t remaining = 0;
  while ((status = next_line(input, line)) > 0) {
    ++line_number;
    if (!tokenize(line, tokens)) {
      return fail(error, "invalid control character.", line_number);
    }
    if (tokens.empty()) {
      continue;
    }
    if (remaining != 0) {
      Block& block = parsed.blocks.back();
      if (tokens.size() > remaining) {
        return fail(error, "too many values in block '" + block.name + "'.", line_number);
      }
      for (const std::string& token : tokens) {
        double value;
        if (!parse_double(token, value)) {
          return fail(
            error, "expected finite decimal value in block '" + block.name + "'.", line_number);
        }
        block.values.push_back(value);
      }
      remaining -= tokens.size();
    } else if (tokens[0] == "@") {
      if (tokens.size() != 3 || !identifier(tokens[1]) || !parse_size(tokens[2], remaining)) {
        return fail(error, "expected block header: @ <name> <count>.", line_number);
      }
      Block block;
      block.name = lower(tokens[1]);
      block.line = line_number;
      if (remaining > block.values.max_size()) {
        return fail(error, "block count too large.", line_number);
      }

      parsed.blocks.push_back(std::move(block));
    } else {
      if (tokens.size() < 2 || !identifier(tokens[0])) {
        return fail(error, "expected field: <name> <values...>.", line_number);
      }
      Field field;
      field.name = lower(tokens[0]);
      field.values.assign(tokens.begin() + 1, tokens.end());
      field.line = line_number;
      parsed.fields.push_back(std::move(field));
    }
  }
  if (status < 0) {
    return fail(error, "input stream failure.", line_number + 1);
  }
  if (remaining != 0) {
    return fail(
      error, "truncated block '" + parsed.blocks.back().name + "'.", parsed.blocks.back().line);
  }
  if (!parsed.validate(error)) {
    return false;
  }
  *this = std::move(parsed);
  return true;
}

bool NEP_File::check_unhandled(
  const std::vector<std::string>& handled_fields,
  const std::vector<std::string>& handled_blocks,
  std::string& error) const
{
  error.clear();
  for (const Field& field : fields) {
    if (!contains(handled_fields, field.name)) {
      return fail(error, "unhandled field '" + field.name + "'.", field.line);
    }
  }
  for (const Block& block : blocks) {
    if (!contains(handled_blocks, block.name)) {
      return fail(error, "unhandled block '" + block.name + "'.", block.line);
    }
  }
  return true;
}

bool NEP_Model_File::fail(std::string& error, const std::string& message, std::size_t line)
{
  error = "NEP model" + (line ? " line " + std::to_string(line) : "") + ": " + message;
  return false;
}

bool NEP_Model_File::read_integer(
  const NEP_File::Field& field,
  std::size_t index,
  int minimum,
  int maximum,
  int& value,
  std::string& error)
{
  std::size_t parsed;
  if (
    !NEP_File::parse_size(field.values[index], parsed) ||
    parsed < static_cast<std::size_t>(minimum) || parsed > static_cast<std::size_t>(maximum)) {
    return fail(
      error,
      "'" + field.name + "' requires an integer in [" + std::to_string(minimum) + ", " +
        std::to_string(maximum) + "].",
      field.line);
  }
  value = static_cast<int>(parsed);
  return true;
}

// Validate the float range without rounding CPU model values to float.
static bool valid_float(double value)
{
  return std::isfinite(value) && std::abs(value) <= std::numeric_limits<float>::max() &&
         (value == 0.0 || static_cast<float>(value) != 0.0);
}

bool NEP_Model_File::read_float(
  const NEP_File::Field& field, std::size_t index, double& value, std::string& error)
{
  double parsed;
  if (!NEP_File::parse_double(field.values[index], parsed) || !valid_float(parsed)) {
    return fail(error, "invalid float value in '" + field.name + "'.", field.line);
  }
  value = parsed;
  return true;
}

const NEP_File::Field* NEP_Model_File::require_field(
  const NEP_File& document,
  const char* name,
  std::size_t minimum,
  std::size_t maximum,
  std::vector<std::string>& handled,
  std::string& error)
{
  const NEP_File::Field* field = document.find_field(name);
  if (field == nullptr) {
    fail(error, "missing field '" + std::string(name) + "'.");
    return nullptr;
  }
  if (field->values.size() < minimum || field->values.size() > maximum) {
    fail(error, "incorrect number of values in '" + field->name + "'.", field->line);
    return nullptr;
  }
  handled.push_back(name);
  return field;
}

bool NEP_Model_File::read_block(
  const NEP_File& document,
  const char* name,
  int count,
  std::vector<double>& values,
  std::vector<std::string>& handled,
  std::string& error)
{
  const NEP_File::Block* block = document.find_block(name);
  if (block == nullptr) {
    return fail(error, "missing block '" + std::string(name) + "'.");
  }
  if (block->values.size() != static_cast<std::size_t>(count)) {
    return fail(
      error,
      "block '" + block->name + "' requires " + std::to_string(count) + " values, found " +
        std::to_string(block->values.size()) + ".",
      block->line);
  }
  values.resize(count);
  for (int i = 0; i < count; ++i) {
    if (!valid_float(block->values[i])) {
      return fail(error, "invalid float value in block '" + block->name + "'.", block->line);
    }
    values[i] = block->values[i];
  }
  handled.push_back(name);
  return true;
}

bool NEP_Model_File::read_descriptor(
  const NEP_File& document, std::vector<std::string>& handled, std::string& error)
{
  const NEP_File::Field* field =
    require_field(document, "cutoff", 4, 2 * num_types + 2, handled, error);
  if (field == nullptr) {
    return false;
  }
  has_multiple_cutoffs = field->values.size() != 4;
  const int count = has_multiple_cutoffs ? num_types : 1;
  if (field->values.size() != static_cast<std::size_t>(2 * count + 2)) {
    return fail(
      error, "cutoff requires two radii or one pair per element, followed by two neighbor counts.",
      field->line);
  }
  rc_radial.resize(num_types);
  rc_angular.resize(num_types);
  for (int t = 0; t < num_types; ++t) {
    const int offset = has_multiple_cutoffs ? 2 * t : 0;
    if (
      !read_float(*field, offset, rc_radial[t], error) ||
      !read_float(*field, offset + 1, rc_angular[t], error)) {
      return false;
    }
    if (rc_angular[t] < 3.0 || rc_angular[t] > rc_radial[t] || rc_radial[t] > 100.0) {
      return fail(error, "cutoff requires 3 <= angular <= radial <= 100 Angstrom.", field->line);
    }
  }
  const int max_count = std::numeric_limits<int>::max();
  if (
    !read_integer(*field, 2 * count, 0, max_count, max_radial, error) ||
    !read_integer(*field, 2 * count + 1, 0, max_count, max_angular, error)) {
    return false;
  }
  field = require_field(document, "n_max", 2, 2, handled, error);
  if (
    field == nullptr || !read_integer(*field, 0, 0, 15, n_max_radial, error) ||
    !read_integer(*field, 1, 0, 12, n_max_angular, error)) {
    return false;
  }
  field = require_field(document, "basis_size", 2, 2, handled, error);
  if (
    field == nullptr || !read_integer(*field, 0, 0, 16, basis_size_radial, error) ||
    !read_integer(*field, 1, 0, 12, basis_size_angular, error)) {
    return false;
  }
  field = require_field(document, "l_max", 3, 7, handled, error);
  if (field == nullptr || !read_integer(*field, 0, 2, 8, L_max, error)) {
    return false;
  }
  int flags[6] = {0, 0, 0, 0, 0, 0};
  for (std::size_t i = 1; i < field->values.size(); ++i) {
    if (!read_integer(*field, i, 0, i == 1 ? 2 : 1, flags[i - 1], error)) {
      return false;
    }
  }
  if (flags[0] == 1) {
    return fail(error, "the q_222 value in l_max must be 0 or 2.", field->line);
  }
  has_q_222 = flags[0] != 0;
  has_q_1111 = flags[1];
  has_q_112 = flags[2];
  has_q_123 = flags[3];
  has_q_233 = flags[4];
  has_q_134 = flags[5];
  layout.dim_radial = n_max_radial + 1;
  layout.dim_angular = (n_max_angular + 1) * (L_max + has_q_222 + has_q_1111 + has_q_112 +
                                              has_q_123 + has_q_233 + has_q_134);
  if (layout.dim_angular > 90) {
    return fail(error, "number of angular descriptors exceeds 90.", field->line);
  }
  layout.dim = layout.dim_radial + layout.dim_angular;
  return true;
}

bool NEP_Model_File::read_ann(
  const NEP_File& document, std::vector<std::string>& handled, std::string& error)
{
  const NEP_File::Field* field = require_field(document, "ann", 2, 2, handled, error);
  if (
    field == nullptr || !read_integer(*field, 0, 1, 120, num_neurons1, error) ||
    !read_integer(*field, 1, 0, 120, num_neurons2, error)) {
    return false;
  }
  if (num_neurons2 != 0) {
    return fail(error, "NEP_CPU supports only one hidden layer.", field->line);
  }
  return true;
}

bool NEP_Model_File::read_reference(
  const NEP_File& document, std::vector<std::string>& handled, std::string& error)
{
  const NEP_File::Field* field =
    require_field(document, "energy_reference", num_types, num_types, handled, error);
  if (field == nullptr) {
    return false;
  }
  energy_reference.resize(num_types);
  for (int t = 0; t < num_types; ++t) {
    if (!NEP_File::parse_double(field->values[t], energy_reference[t])) {
      return fail(error, "invalid double value in energy_reference.", field->line);
    }
  }
  return true;
}

bool NEP_Model_File::read_zbl(
  const NEP_File& document,
  std::vector<std::string>& fields,
  std::vector<std::string>& blocks,
  std::string& error)
{
  if (document.find_field("zbl") == nullptr) {
    return true;
  }
  const NEP_File::Field* field = require_field(document, "zbl", 2, 3, fields, error);
  if (
    field == nullptr || !read_float(*field, 0, zbl_rc_inner, error) ||
    !read_float(*field, 1, zbl_rc_outer, error)) {
    return false;
  }
  enable_zbl = true;
  flexible_zbl = zbl_rc_inner == 0.0 && zbl_rc_outer == 0.0;
  if (flexible_zbl) {
    if (name != "nep") {
      return fail(error, "Flexible ZBL is only supported for model nep.", field->line);
    }
    if (field->values.size() != 2) {
      return fail(error, "flexible ZBL cannot use a typewise cutoff factor.", field->line);
    }
    const int count = 10 * num_types * (num_types + 1) / 2;
    if (count > 550) {
      return fail(error, "NEP_CPU supports flexible ZBL for at most 10 elements.", field->line);
    }
    if (!read_block(document, "zbl_parameters", count, zbl_parameters, blocks, error)) {
      return false;
    }
    for (int i = 0; i < count; i += 10) {
      if (zbl_parameters[i] < 0.0 || zbl_parameters[i + 1] <= zbl_parameters[i]) {
        return fail(
          error, "each zbl_parameters pair requires 0 <= inner < outer.",
          document.find_block("zbl_parameters")->line);
      }
    }
  } else {
    if (zbl_rc_inner < 0.0 || zbl_rc_outer <= zbl_rc_inner) {
      return fail(error, "zbl requires 0 <= inner < outer, or 0 0 for flexible ZBL.", field->line);
    }
    if (field->values.size() == 3) {
      use_typewise_cutoff_zbl = true;
      if (!read_float(*field, 2, typewise_cutoff_zbl_factor, error)) {
        return false;
      }
      if (typewise_cutoff_zbl_factor < 0.5) {
        return fail(error, "typewise ZBL cutoff factor must be at least 0.5.", field->line);
      }
    }
  }
  return true;
}

void NEP_Model_File::calculate_layout(
  int extra_inputs, int extra_outputs, int globals, int ann_count)
{
  layout.dim += extra_inputs;
  layout.ann_count = ann_count;
  layout.global_parameters = globals;
  layout.number_of_variables_ann_1 = (layout.dim + 2 + extra_outputs) * num_neurons1;
  layout.number_of_variables_ann = layout.number_of_variables_ann_1 * num_types + globals;
  layout.number_of_variables_descriptor =
    num_types * num_types *
    ((n_max_radial + 1) * (basis_size_radial + 1) + (n_max_angular + 1) * (basis_size_angular + 1));
  layout.number_of_variables =
    layout.number_of_variables_ann * ann_count + layout.number_of_variables_descriptor;
}

bool NEP_Model_File::read_parameters(
  const NEP_File& document, std::vector<std::string>& blocks, std::string& error)
{
  if (layout.dim > 103) {
    return fail(error, "descriptor dimension exceeds 103.", document.find_field("l_max")->line);
  }
  return read_block(
           document, "parameters", layout.number_of_variables, parameters, blocks, error) &&
         read_block(document, "q_scaler", layout.dim, q_scaler, blocks, error);
}

bool NEP_Model_File::read(std::istream& input, std::string& error)
{
  *this = NEP_Model_File();
  NEP_File document;
  if (!document.read(input, error)) {
    return false;
  }
  return parse(document, error);
}

bool NEP_Model_File::parse(const NEP_File& document, std::string& error)
{
  NEP_Model_File parsed;
  parsed.name = document.find_field("model")->values[0];
  parsed.elements = document.elements;
  parsed.num_types = static_cast<int>(parsed.elements.size());
  std::vector<std::string> fields = {"format_version", "model"};
  std::vector<std::string> blocks;
  bool success = false;
  if (parsed.name == "nep") {
    success = parsed.read_nep(document, fields, blocks, error);
  } else if (parsed.name == "qnep1" || parsed.name == "qnep2") {
    success = parsed.read_charge(document, fields, blocks, error);
  } else if (parsed.name == "nep_dipole" || parsed.name == "nep_polarizability") {
    success = parsed.read_tnep(document, fields, blocks, error);
  } else {
    return fail(
      error, "unsupported model '" + parsed.name + "'.", document.find_field("model")->line);
  }
  if (!success || !document.check_unhandled(fields, blocks, error)) {
    return false;
  }
  *this = std::move(parsed);
  return true;
}

bool NEP_Model_File::read_nep(
  const NEP_File& document,
  std::vector<std::string>& fields,
  std::vector<std::string>& blocks,
  std::string& error)
{
  if (!read_descriptor(document, fields, error) || !read_ann(document, fields, error)) {
    return false;
  }
  calculate_layout(0, 0, 1, 1);
  return read_reference(document, fields, error) && read_zbl(document, fields, blocks, error) &&
         read_parameters(document, blocks, error);
}

bool NEP_Model_File::read_charge(
  const NEP_File& document,
  std::vector<std::string>& fields,
  std::vector<std::string>& blocks,
  std::string& error)
{
  if (!read_descriptor(document, fields, error) || !read_ann(document, fields, error)) {
    return false;
  }
  if (has_multiple_cutoffs) {
    return fail(
      error, name + " supports only uniform cutoffs.", document.find_field("cutoff")->line);
  }
  calculate_layout(0, 1, 2, 1);
  return read_reference(document, fields, error) && read_zbl(document, fields, blocks, error) &&
         read_parameters(document, blocks, error);
}

bool NEP_Model_File::read_tnep(
  const NEP_File& document,
  std::vector<std::string>& fields,
  std::vector<std::string>& blocks,
  std::string& error)
{
  if (!read_descriptor(document, fields, error) || !read_ann(document, fields, error)) {
    return false;
  }
  calculate_layout(0, 0, 1, name == "nep_polarizability" ? 2 : 1);
  return read_parameters(document, blocks, error);
}
