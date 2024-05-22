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

#ifdef PAIR_CLASS

PairStyle(nep, PairNEP)

#else

#ifndef LMP_PAIR_NEP_H
#define LMP_PAIR_NEP_H

#include "nep.h"
#include "pair.h"
#include <string>

namespace LAMMPS_NS
{
class PairNEP : public Pair
{
public:
  double cutoff;
  NEP3 nep_model;
  PairNEP(class LAMMPS*);
  virtual ~PairNEP();
  virtual void coeff(int, char**);
  virtual void settings(int, char**);
  virtual double init_one(int, int);
  virtual void init_style();
  virtual void compute(int, int);

protected:
  bool inited;
  char model_filename[1000];
  double cutoffsq;
  void allocate();
  
  std::string emap;   //map
};
} // namespace LAMMPS_NS

#endif
#endif

