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
};
} // namespace LAMMPS_NS

#endif
#endif
