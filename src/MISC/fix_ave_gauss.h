/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ave/gauss,FixAveGauss)

#else

#ifndef LMP_FIX_AVE_GAUSS_H
#define LMP_FIX_AVE_GAUSS_H

#include "fix.h"


/**

fix 42 all ave/gauss v_foo v_bar window 500 every 100 delay 0 delay 5

- window: int>0
  number of timesteps to use for statistics
- every: int>0
  update statistics every this many timesteps
  values are valid but not updated in between
- delay: int>0
  keep the last N updates ("every"*"delay" timesteps ago)
  can be given more than once, if none is given, "delay 0" is implied
  values can be queried by accessing the computed array global in order given

f_42[i]    => value of i at first delay
f_42[i][j] => value of i at jth delay

**/


namespace LAMMPS_NS {

class FixAveGauss : public Fix {
 public:
  FixAveGauss(class LAMMPS *, int, char **);
  ~FixAveGauss();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);

 private:
  int me;
  
  int nfull;
  bigint nvalid, nvalid_last, nfull_next;

  int nvalues;
  char **ids;
  int *varindex;

  int nwindow;
  int iwindow, window_limit;
  double **window_list;
  
  int ndelay,nresult;
  int iresult;
  int *delays;
  double **result_list;

  void append_values(bigint);
  void update_results(bigint);
  void options(int, int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No values in fix ave/gauss command

Self-explanatory.

E: Invalid fix ave/gauss off column

Self-explanatory.

E: Variable name for fix ave/gauss does not exist

Self-explanatory.

E: Fix ave/gauss variable is not equal-style variable

Self-explanatory.

E: Invalid timestep reset for fix ave/gauss

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

*/
