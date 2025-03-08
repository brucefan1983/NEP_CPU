/*
    Copyright 2022 Zheyong Fan, Junjie Wang, Eric Lindgren
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

#pragma once
#include <string>
#include <vector>

// #define USE_TABLE_FOR_RADIAL_FUNCTIONS

// NEP namespace for parameters, structs, classes and functions
namespace NEP_NS {

// ---------- PARAMETER ----------

const int MAX_NEURON = 200; // maximum number of neurons in the hidden layer
const int MN = 1000;        // maximum number of neighbors for one atom
const int NUM_OF_ABC = 80;  // 3 + 5 + 7 + 9 + 11 + 13 + 15 + 17 for L_max = 8
const int MAX_NUM_N = 20;   // n_max+1 = 19+1
const int MAX_DIM = MAX_NUM_N * 7;
const int MAX_DIM_ANGULAR = MAX_NUM_N * 6;
const double C3B[NUM_OF_ABC] = {
  0.238732414637843, 0.119366207318922, 0.119366207318922, 0.099471839432435, 0.596831036594608,
  0.596831036594608, 0.149207759148652, 0.149207759148652, 0.139260575205408, 0.104445431404056,
  0.104445431404056, 1.044454314040563, 1.044454314040563, 0.174075719006761, 0.174075719006761,
  0.011190581936149, 0.223811638722978, 0.223811638722978, 0.111905819361489, 0.111905819361489,
  1.566681471060845, 1.566681471060845, 0.195835183882606, 0.195835183882606, 0.013677377921960,
  0.102580334414698, 0.102580334414698, 2.872249363611549, 2.872249363611549, 0.119677056817148,
  0.119677056817148, 2.154187022708661, 2.154187022708661, 0.215418702270866, 0.215418702270866,
  0.004041043476943, 0.169723826031592, 0.169723826031592, 0.106077391269745, 0.106077391269745,
  0.424309565078979, 0.424309565078979, 0.127292869523694, 0.127292869523694, 2.800443129521260,
  2.800443129521260, 0.233370260793438, 0.233370260793438, 0.004662742473395, 0.004079899664221,
  0.004079899664221, 0.024479397985326, 0.024479397985326, 0.012239698992663, 0.012239698992663,
  0.538546755677165, 0.538546755677165, 0.134636688919291, 0.134636688919291, 3.500553911901575,
  3.500553911901575, 0.250039565135827, 0.250039565135827, 0.000082569397966, 0.005944996653579,
  0.005944996653579, 0.104037441437634, 0.104037441437634, 0.762941237209318, 0.762941237209318,
  0.114441185581398, 0.114441185581398, 5.950941650232678, 5.950941650232678, 0.141689086910302,
  0.141689086910302, 4.250672607309055, 4.250672607309055, 0.265667037956816, 0.265667037956816};
const double C4B[5] = {
  -0.007499480826664, -0.134990654879954, 0.067495327439977, 0.404971964639861, -0.809943929279723};
const double C5B[3] = {0.026596810706114, 0.053193621412227, 0.026596810706114};

const double Z_COEFFICIENT_1[2][2] = {{0.0, 1.0}, {1.0, 0.0}};

const double Z_COEFFICIENT_2[3][3] = {{-1.0, 0.0, 3.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}};

const double Z_COEFFICIENT_3[4][4] = {
  {0.0, -3.0, 0.0, 5.0}, {-1.0, 0.0, 5.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0}};

const double Z_COEFFICIENT_4[5][5] = {
  {3.0, 0.0, -30.0, 0.0, 35.0},
  {0.0, -3.0, 0.0, 7.0, 0.0},
  {-1.0, 0.0, 7.0, 0.0, 0.0},
  {0.0, 1.0, 0.0, 0.0, 0.0},
  {1.0, 0.0, 0.0, 0.0, 0.0}};

const double Z_COEFFICIENT_5[6][6] = {
  {0.0, 15.0, 0.0, -70.0, 0.0, 63.0}, {1.0, 0.0, -14.0, 0.0, 21.0, 0.0},
  {0.0, -1.0, 0.0, 3.0, 0.0, 0.0},    {-1.0, 0.0, 9.0, 0.0, 0.0, 0.0},
  {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

const double Z_COEFFICIENT_6[7][7] = {
  {-5.0, 0.0, 105.0, 0.0, -315.0, 0.0, 231.0}, {0.0, 5.0, 0.0, -30.0, 0.0, 33.0, 0.0},
  {1.0, 0.0, -18.0, 0.0, 33.0, 0.0, 0.0},      {0.0, -3.0, 0.0, 11.0, 0.0, 0.0, 0.0},
  {-1.0, 0.0, 11.0, 0.0, 0.0, 0.0, 0.0},       {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

const double Z_COEFFICIENT_7[8][8] = {{0.0, -35.0, 0.0, 315.0, 0.0, -693.0, 0.0, 429.0},
                                      {-5.0, 0.0, 135.0, 0.0, -495.0, 0.0, 429.0, 0.0},
                                      {0.0, 15.0, 0.0, -110.0, 0.0, 143.0, 0.0, 0.0},
                                      {3.0, 0.0, -66.0, 0.0, 143.0, 0.0, 0.0, 0.0},
                                      {0.0, -3.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0},
                                      {-1.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                      {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                      {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

const double Z_COEFFICIENT_8[9][9] = {
  {35.0, 0.0, -1260.0, 0.0, 6930.0, 0.0, -12012.0, 0.0, 6435.0},
  {0.0, -35.0, 0.0, 385.0, 0.0, -1001.0, 0.0, 715.0, 0.0},
  {-1.0, 0.0, 33.0, 0.0, -143.0, 0.0, 143.0, 0.0, 0.0},
  {0.0, 3.0, 0.0, -26.0, 0.0, 39.0, 0.0, 0.0, 0.0},
  {1.0, 0.0, -26.0, 0.0, 65.0, 0.0, 0.0, 0.0, 0.0},
  {0.0, -1.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {-1.0, 0.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

const double K_C_SP = 14.399645; // 1/(4*PI*epsilon_0)
const double PI = 3.141592653589793;
const double PI_HALF = 1.570796326794897;
const int NUM_ELEMENTS = 94;
const std::string ELEMENTS[NUM_ELEMENTS] = {
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
  "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
  "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu"};
double COVALENT_RADIUS[NUM_ELEMENTS] = {
  0.426667, 0.613333, 1.6,     1.25333, 1.02667, 1.0,     0.946667, 0.84,    0.853333, 0.893333,
  1.86667,  1.66667,  1.50667, 1.38667, 1.46667, 1.36,    1.32,     1.28,    2.34667,  2.05333,
  1.77333,  1.62667,  1.61333, 1.46667, 1.42667, 1.38667, 1.33333,  1.32,    1.34667,  1.45333,
  1.49333,  1.45333,  1.53333, 1.46667, 1.52,    1.56,    2.52,     2.22667, 1.96,     1.85333,
  1.76,     1.65333,  1.53333, 1.50667, 1.50667, 1.44,    1.53333,  1.64,    1.70667,  1.68,
  1.68,     1.64,     1.76,    1.74667, 2.78667, 2.34667, 2.16,     1.96,    2.10667,  2.09333,
  2.08,     2.06667,  2.01333, 2.02667, 2.01333, 2.0,     1.98667,  1.98667, 1.97333,  2.04,
  1.94667,  1.82667,  1.74667, 1.64,    1.57333, 1.54667, 1.48,     1.49333, 1.50667,  1.76,
  1.73333,  1.73333,  1.81333, 1.74667, 1.84,    1.89333, 2.68,     2.41333, 2.22667,  2.10667,
  2.02667,  2.04,     2.05333, 2.06667};

// ---------- STRUCT ----------

// struct for many body parameters and settings
struct ParaMB {
    bool use_typewise_cutoff = false;
    bool use_typewise_cutoff_zbl = false;
    double typewise_cutoff_radial_factor = 2.5;
    double typewise_cutoff_angular_factor = 2.0;
    double typewise_cutoff_zbl_factor = 0.65;

    int model_type = 0; // 0=potential, 1=dipole, 2=polarizability
    int version = 4;
    double rc_radial = 0.0;
    double rc_angular = 0.0;
    double rcinv_radial = 0.0;
    double rcinv_angular = 0.0;
    int n_max_radial = 0;
    int n_max_angular = 0;
    int L_max = 0;
    int dim_angular;
    int num_L;
    int basis_size_radial = 8;
    int basis_size_angular = 8;
    int num_types_sq = 0;
    int num_c_radial = 0;
    int num_types = 0;
    double q_scaler[140];
    int atomic_numbers[94];
  };

// struct for neural network parameters
struct ANN {
  int dim = 0;
  int num_neurons1 = 0;
  int num_para = 0;
  int num_para_ann = 0;
  const double* w0[94];
  const double* b0[94];
  const double* w1[94];
  const double* b1;
  const double* c;
  // for the scalar part of polarizability
  const double* w0_pol[94];
  const double* b0_pol[94];
  const double* w1_pol[94];
  const double* b1_pol;
};

// struct for ZBL parameters and settings
struct ZBL {
  bool enabled = false;
  bool flexibled = false;
  int num_types;
  double rc_inner = 1.0;
  double rc_outer = 2.0;
  double para[550];
};

// struct for D3 parameters
struct DFTD3 {
  double s6 = 0.0;
  double s8 = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double rc_radial = 20.0;
  double rc_angular = 10.0;
  int atomic_number[94]; // H to Pu
  std::vector<double> cn;
  std::vector<double> dc6_sum;
  std::vector<double> dc8_sum;
};


// ---------- CLASS ----------

// NEP class
class NEP3
{
public:
  NEP3();
  NEP3(const std::string& potential_filename);

  void init_from_file(const std::string& potential_filename, const bool is_rank_0);

  void update_type_map(const int ntype, int* type_map, char** elements);

  // type[num_atoms] should be integers 0, 1, ..., mapping to the atom types in nep.txt in order
  // box[9] is ordered as ax, bx, cx, ay, by, cy, az, bz, cz
  // position[num_atoms * 3] is ordered as x[num_atoms], y[num_atoms], z[num_atoms]
  // potential[num_atoms]
  // force[num_atoms * 3] is ordered as fx[num_atoms], fy[num_atoms], fz[num_atoms]
  // virial[num_atoms * 9] is ordered as v_xx[num_atoms], v_xy[num_atoms], v_xz[num_atoms],
  // v_yx[num_atoms], v_yy[num_atoms], v_yz[num_atoms], v_zx[num_atoms], v_zy[num_atoms],
  // v_zz[num_atoms]
  // descriptor[num_atoms * dim] is ordered as d0[num_atoms], d1[num_atoms], ...

  void compute(
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& potential,
    std::vector<double>& force,
    std::vector<double>& virial);

  void compute_with_dftd3(
    const std::string& xc,
    const double rc_potential,
    const double rc_coordination_number,
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& potential,
    std::vector<double>& force,
    std::vector<double>& virial);

  void compute_dftd3(
    const std::string& xc,
    const double rc_potential,
    const double rc_coordination_number,
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& potential,
    std::vector<double>& force,
    std::vector<double>& virial);

  void find_descriptor(
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& descriptor);

  void find_latent_space(
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& latent_space);

  void find_B_projection(
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& B_projection);

  void find_dipole(
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& dipole // 3 components, for the whole box
  );

  void find_polarizability(
    const std::vector<int>& type,
    const std::vector<double>& box,
    const std::vector<double>& position,
    std::vector<double>& polarizability // 6 components, for the whole box
  );

  void compute_for_lammps(
    int nlocal,              // atom->nlocal
    int inum,                // list->inum
    int* ilist,              // list->ilist
    int* numneigh,           // list->numneigh
    int** firstneigh,        // list->firstneigh
    int* type,               // atom->type
    int* type_map,           // map from atom type to element
    double** x,              // atom->x
    double& total_potential, // total potential energy for the current processor
    double total_virial[6],  // total virial for the current processor
    double* potential,       // eatom or nullptr
    double** f,              // atom->f
    double** virial          // cvatom or nullptr
  );

  int num_atoms = 0;
  int num_cells[3];
  double ebox[18];
  ParaMB paramb;
  ANN annmb;
  ZBL zbl;
  DFTD3 dftd3;
  std::vector<int> NN_radial, NL_radial, NN_angular, NL_angular;
  std::vector<double> r12;
  std::vector<double> Fp;
  std::vector<double> sum_fxyz;
  std::vector<double> parameters;
  std::vector<std::string> element_list;
  void update_potential(double* parameters, ANN& ann);
  void allocate_memory(const int N);

#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  std::vector<double> gn_radial;   // tabulated gn_radial functions
  std::vector<double> gnp_radial;  // tabulated gnp_radial functions
  std::vector<double> gn_angular;  // tabulated gn_angular functions
  std::vector<double> gnp_angular; // tabulated gnp_angular functions
  void construct_table(double* parameters);
#endif

  bool set_dftd3_para_one(
    const std::string& functional_input,
    const std::string& functional_library,
    const double s6,
    const double a1,
    const double s8,
    const double a2);
  void set_dftd3_para_all(
    const std::string& functional_input,
    const double rc_potential,
    const double rc_coordination_number);
};

// ---------- FUNCTION ----------

void complex_product(const double a, const double b, double& real_part, double& imag_part);

void apply_ann_one_layer(
  const int dim,
  const int num_neurons1,
  const double* w0,
  const double* b0,
  const double* w1,
  const double* b1,
  double* q,
  double& energy,
  double* energy_derivative,
  double* latent_space,
  bool need_B_projection,
  double* B_projection);

void apply_ann_one_layer_nep5(
  const int dim,
  const int num_neurons1,
  const double* w0,
  const double* b0,
  const double* w1,
  const double* b1,
  double* q,
  double& energy,
  double* energy_derivative,
  double* latent_space);


void find_fc(double rc, double rcinv, double d12, double& fc);


void find_fc_and_fcp(double rc, double rcinv, double d12, double& fc, double& fcp);


void find_fc_and_fcp_zbl(double r1, double r2, double d12, double& fc, double& fcp);


void find_phi_and_phip_zbl(double a, double b, double x, double& phi, double& phip);


void find_f_and_fp_zbl(
  double zizj,
  double a_inv,
  double rc_inner,
  double rc_outer,
  double d12,
  double d12inv,
  double& f,
  double& fp);


void find_f_and_fp_zbl(
  double* zbl_para, double zizj, double a_inv, double d12, double d12inv, double& f, double& fp);


void find_fn(const int n, const double rcinv, const double d12, const double fc12, double& fn);


void find_fn_and_fnp(
  const int n,
  const double rcinv,
  const double d12,
  const double fc12,
  const double fcp12,
  double& fn,
  double& fnp);


void find_fn(const int n_max, const double rcinv, const double d12, const double fc12, double* fn);


void find_fn_and_fnp(
  const int n_max,
  const double rcinv,
  const double d12,
  const double fc12,
  const double fcp12,
  double* fn,
  double* fnp);


void get_f12_4body(
  const double d12,
  const double d12inv,
  const double fn,
  const double fnp,
  const double Fp,
  const double* s,
  const double* r12,
  double* f12);


void get_f12_5body(
  const double d12,
  const double d12inv,
  const double fn,
  const double fnp,
  const double Fp,
  const double* s,
  const double* r12,
  double* f12);


template <int L>
void calculate_s_one(
  const int n, const int n_max_angular_plus_1, const double* Fp, const double* sum_fxyz, double* s);


template <int L>
void accumulate_f12_one(
  const double d12inv,
  const double fn,
  const double fnp,
  const double* s,
  const double* r12,
  double* f12);


void accumulate_f12(
  const int L_max,
  const int num_L,
  const int n,
  const int n_max_angular_plus_1,
  const double d12,
  const double* r12,
  double fn,
  double fnp,
  const double* Fp,
  const double* sum_fxyz,
  double* f12);


template <int L>
void accumulate_s_one(
  const double x12, const double y12, const double z12, const double fn, double* s);


void accumulate_s(
  const int L_max, const double d12, double x12, double y12, double z12, const double fn, double* s);


template <int L>
double find_q_one(const double* s);


void find_q(
  const int L_max,
  const int num_L,
  const int n_max_angular_plus_1,
  const int n,
  const double* s,
  double* q);


#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
const int table_length = 2001;
const int table_segments = table_length - 1;
const double table_resolution = 0.0005;

void find_index_and_weight(
  const double d12_reduced,
  int& index_left,
  int& index_right,
  double& weight_left,
  double& weight_right);


void construct_table_radial_or_angular(
  const int version,
  const int num_types,
  const int num_types_sq,
  const int n_max,
  const int basis_size,
  const double rc,
  const double rcinv,
  const double* c,
  double* gn,
  double* gnp);

#endif

void find_descriptor_small_box(
  const bool calculating_potential,
  const bool calculating_descriptor,
  const bool calculating_latent_space,
  const bool calculating_polarizability,
  ParaMB& paramb,
  ANN& annmb,
  const int N,
  const int* g_NN_radial,
  const int* g_NL_radial,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12_radial,
  const double* g_y12_radial,
  const double* g_z12_radial,
  const double* g_x12_angular,
  const double* g_y12_angular,
  const double* g_z12_angular,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_radial,
  const double* g_gn_angular,
#endif
  double* g_Fp,
  double* g_sum_fxyz,
  double* g_potential,
  double* g_descriptor,
  double* g_latent_space,
  double* g_virial,
  bool calculating_B_projection,
  double* g_B_projection);


void find_force_radial_small_box(
  const bool is_dipole,
  ParaMB& paramb,
  ANN& annmb,
  const int N,
  const int* g_NN,
  const int* g_NL,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_Fp,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gnp_radial,
#endif
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial);


void find_force_angular_small_box(
  const bool is_dipole,
  ParaMB& paramb,
  ANN& annmb,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_Fp,
  const double* g_sum_fxyz,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_angular,
  const double* g_gnp_angular,
#endif
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial);


void find_force_ZBL_small_box(
  const int N,
  ParaMB& paramb,
  const ZBL& zbl,
  const int* g_NN,
  const int* g_NL,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_pe);


void find_dftd3_coordination_number(
  DFTD3& dftd3,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12);


void add_dftd3_force(
  DFTD3& dftd3,
  const int N,
  const int* g_NN_radial,
  const int* g_NL_radial,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_potential,
  double* g_force,
  double* g_virial);


void add_dftd3_force_extra(
  const DFTD3& dftd3,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_force,
  double* g_virial);


void find_descriptor_for_lammps(
  ParaMB& paramb,
  ANN& annmb,
  int nlocal,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_radial,
  const double* g_gn_angular,
#endif
  double* g_Fp,
  double* g_sum_fxyz,
  double& g_total_potential,
  double* g_potential);


void find_force_radial_for_lammps(
  ParaMB& paramb,
  ANN& annmb,
  int nlocal,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
  double* g_Fp,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gnp_radial,
#endif
  double** g_force,
  double g_total_virial[6],
  double** g_virial);


void find_force_angular_for_lammps(
  ParaMB& paramb,
  ANN& annmb,
  int nlocal,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
  double* g_Fp,
  double* g_sum_fxyz,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_angular,
  const double* g_gnp_angular,
#endif
  double** g_force,
  double g_total_virial[6],
  double** g_virial);


void find_force_ZBL_for_lammps(
  ParaMB& paramb,
  const ZBL& zbl,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
  double** g_force,
  double g_total_virial[6],
  double** g_virial,
  double& g_total_potential,
  double* g_potential);


double get_area_one_direction(const double* a, const double* b);

double get_area(const int d, const double* cpu_h);

double get_det(const double* cpu_h);

double get_volume(const double* cpu_h);

void get_inverse(double* cpu_h);

bool get_expanded_box(const double rc, const double* box, int* num_cells, double* ebox);


void applyMicOne(double& x12);


void apply_mic_small_box(const double* ebox, double& x12, double& y12, double& z12);


void findCell(
  const double* box,
  const double* thickness,
  const double* r,
  double cutoffInverse,
  const int* numCells,
  int* cell);


void applyPbcOne(double& sx);


void applyPbc(const int N, const double* box, double* x, double* y, double* z);


void find_neighbor_list_large_box(
  const double rc_radial,
  const double rc_angular,
  const int N,
  const std::vector<double>& box,
  const std::vector<double>& position,
  int* num_cells,
  double* ebox,
  std::vector<int>& g_NN_radial,
  std::vector<int>& g_NL_radial,
  std::vector<int>& g_NN_angular,
  std::vector<int>& g_NL_angular,
  std::vector<double>& r12);


void find_neighbor_list_small_box(
  const double rc_radial,
  const double rc_angular,
  const int N,
  const std::vector<double>& box,
  const std::vector<double>& position,
  int* num_cells,
  double* ebox,
  std::vector<int>& g_NN_radial,
  std::vector<int>& g_NL_radial,
  std::vector<int>& g_NN_angular,
  std::vector<int>& g_NL_angular,
  std::vector<double>& r12);


}
