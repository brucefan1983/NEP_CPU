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
#include "nep.h"
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace LAMMPS_NS;

PairNEP::PairNEP(LAMMPS* lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;

  single_enable = 0;

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

  allocated = 1;
}

void PairNEP::coeff(int narg, char** arg)
{
  if (!allocated)
    allocate();
}

void PairNEP::settings(int narg, char** arg)
{
  if (narg != 1)
    error->all(FLERR, "Illegal pair_style command");
  strcpy(model_filename, arg[0]);
}

void PairNEP::init_style()
{
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // if (inited) ~nep_model();
  nep_model.init_from_file(model_filename);
  inited = true;
  cutoff = nep_model.paramb.rc_radial;
  cutoffsq = cutoff * cutoff;
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      cutsq[i][j] = cutoffsq;

  // if (comm->nprocs != 1) error->all(FLERR, "no parallel plz");
}

double PairNEP::init_one(int i, int j) { return cutoff; }

void PairNEP::compute(int eflag, int vflag)
{
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
  double energy = 0;
  double* p_site_en = NULL;
  double** p_site_virial = NULL;
  if (eflag_atom)
    p_site_en = eatom;
  if (vflag_atom)
    p_site_virial = cvatom;
  nep_model.compute_for_lammps(
    list->inum, list->ilist, list->numneigh, list->firstneigh, atom->type, atom->x, energy,
    p_site_en, atom->f, p_site_virial);
  if (eflag_global) {
    // eng_vdwl += energy;    mix with other potential
    eng_vdwl = energy;
  }
  if (vflag_fdotr) {
    virial_fdotr_compute();
  }
}
