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
   Contributing authors: 
   Junjie Wang (Nanjing University)
   Wenhao Luo  (Sun Yat-sen University)
   Xi     Tan  (Huazhong University of Science and Technology)  
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
#include <iostream> 

//#define LAMMPS_VERSION_NUMBER 20240522 // use the new neighbor list starting from this version

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairNEP::PairNEP(LAMMPS* lmp) : Pair(lmp)
{

  single_enable = 0;    // 1 if single() routine exists
  restartinfo = 0;      // 1 if pair style writes restart info
  one_coeff = 1;        // 1 if allows only one coeff * * call
  manybody_flag = 1;    // 1 if a manybody potential
}


/* ---------------------------------------------------------------------- */

PairNEP::~PairNEP()
{
  if (copymode)
    return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNEP::allocate()
{
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 1;
    }
  
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");    
  allocated = 1;
}



/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNEP::settings(int narg, char** arg)
{

    // default values
    
    // process optional keywords  //now is none!
        
    if (narg != 0)
      error->all(FLERR, "Illegal pair_style command");	
}



/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNEP::coeff(int narg, char** arg)
{
  int n = atom->ntypes;

  if (!allocated) allocate();

  if (narg != 3 + n) error->all(FLERR, "Incorrect args for pair coefficients");

  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  int *map = new int[n + 1];
  for (int i = 0; i < n; i++) map[i] = 0;

  emap = "";
   for (int i = 3; i < narg; i++) {
     if (strcmp(arg[i], "NULL") != 0) {
       if (!emap.empty()) emap += ",";
       emap += std::to_string(i - 1) + ":" + arg[i];
       map[i - 2] = 1;
     } else {
       map[i - 2] = -1;
     }
   }

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] > 0 && map[j] > 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
  
   // read potential file 
   strcpy(model_filename, arg[2]);
   
  delete[] map;
}




/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairNEP::init_style()
{
 
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style NEP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style NEP requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this,NeighConst::REQ_FULL);
 

  bool is_rank_0 = (comm->me == 0);
  nep_model.init_from_file(model_filename, is_rank_0);
  inited = true;
  cutoff = nep_model.paramb.rc_radial;
  cutoffsq = cutoff * cutoff;
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      cutsq[i][j] = cutoffsq;
}



/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNEP::init_one(int i, int j) 
{ 

  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutoff;
}



/* ---------------------------------------------------------------------- */
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
   atom->nlocal, list->inum, list->ilist, list->numneigh, list->firstneigh, atom->type, atom->x, total_potential,
    total_virial, per_atom_potential, atom->f, per_atom_virial);

  if (eflag) {
    eng_vdwl += total_potential;
  }
  if (vflag) {
    for (int component = 0; component < 6; ++component) {
      virial[component] += total_virial[component];
    }
  }
}
