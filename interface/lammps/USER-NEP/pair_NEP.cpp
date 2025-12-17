/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Junjie Wang (Nanjing University)
------------------------------------------------------------------------- */

#include "pair_NEP.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define LAMMPS_VERSION_NUMBER 20220324 // use the new neighbor list starting from this version

using namespace LAMMPS_NS;

PairNEP::PairNEP(LAMMPS* lmp) : Pair(lmp)
{
#if LAMMPS_VERSION_NUMBER >= 20201130
  centroidstressflag = CENTROID_AVAIL;
#else
  centroidstressflag = 2;
#endif

  restartinfo = 0;
  manybody_flag = 1;

  single_enable = 0;
  one_coeff = 1;

  inited = false;
  allocated = 0;
}

PairNEP::~PairNEP()
{
  if (copymode)
    return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete[] type_map;
  }
}

void PairNEP::allocate()
{
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  map = new int[n + 1];
  type_map = new int[n + 1];
  allocated = 1;
}

void PairNEP::coeff(int narg, char** arg)
{
  if (!allocated)
    allocate();

  bool is_rank_0 = (comm->me == 0);
  model_filename = arg[2];

  // get atom types in lammps input file
  map_element2type(narg-3, arg+3);

  // copy the type map list
  int ntype = atom->ntypes;
  for (int i = 0; i <= ntype; ++i) {
    type_map[i] = map[i];
  }

  // read the potential file
  std::string model_path_from_lammps = utils::get_potential_file_path(model_filename);
  nep_model.init_from_file(model_path_from_lammps, is_rank_0);
  // update the type map to elements in potential file
  nep_model.update_type_map(atom->ntypes, type_map, elements);

  // get cutoff from NEP model
  cutoff = nep_model.paramb.rc_radial_max;
  cutoffsq = cutoff * cutoff;
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      cutsq[i][j] = cutoffsq;
}

void PairNEP::settings(int narg, char** arg)
{
  if (narg > 0) {
    error->all(FLERR, "Illegal pair_style command; nep doesn't require any parameter");
  }
}

void PairNEP::init_style()
{
#if LAMMPS_VERSION_NUMBER >= 20220324
  neighbor->add_request(this, NeighConst::REQ_FULL);
#else
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
#endif

}

double PairNEP::init_one(int i, int j) { return cutoff; }

void PairNEP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
  }
  double total_potential = 0.0;
  double total_virial[6] = {0.0};
  double* per_atom_potential = nullptr;
  double** per_atom_virial = nullptr;
  if (eflag_atom) {
    per_atom_potential = eatom;
  }
  if (cvflag_atom) {
    per_atom_virial = cvatom;
  }

  nep_model.compute_for_lammps(
    atom->nlocal, list->inum, list->ilist, list->numneigh, list->firstneigh, atom->type, type_map, 
    atom->x, total_potential, total_virial, per_atom_potential, atom->f, per_atom_virial);

  if (eflag) {
    eng_vdwl += total_potential;
  }
  if (vflag) {
    for (int component = 0; component < 6; ++component) {
      virial[component] += total_virial[component];
    }
  }
}
