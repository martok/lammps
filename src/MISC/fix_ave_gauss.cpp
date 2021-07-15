/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_ave_gauss.h"

#include <cmath>
#include "compute.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

static const int max_delays = 10;

FixAveGauss::FixAveGauss(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0), varindex(nullptr), ids(nullptr), delays(nullptr),
  window_list(nullptr), result_list(nullptr)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/gauss command");

  MPI_Comm_rank(world,&me);

  nevery = 1;
  global_freq = 1;

  dynamic_group_allow = 1;

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"v_",2) == 0) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/gauss command");
  
  delays = new int[max_delays];

  options(iarg,narg,arg);

  // shift args
  arg = &arg[3];

  // parse values

  ids = new char*[nvalues];
  varindex = new int[nvalues];
  
  for (int i = 0; i < nvalues; i++) {
    int n = strlen(arg[i]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[i][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      error->all(FLERR,"Illegal fix ave/gauss command");
    };

    n = strlen(suffix) + 1;
    ids[i] = new char[n];
    strcpy(ids[i],suffix);
    delete [] suffix;
  }

  // setup and error check
  // set variable_length if any compute is variable length
  for (int i = 0; i < nvalues; i++) {
    int ivariable = input->variable->find(ids[i]);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for fix ave/gauss does not exist");
    if (input->variable->equalstyle(ivariable) == 0)
      error->all(FLERR,"Fix ave/gauss variable is not equal-style variable");
  }

  // allocate memory for averaging
 
  window_list = nullptr;
  result_list = nullptr;
  
  // one window of length nwindow  
  memory->create(window_list, nwindow, nvalues, "ave/gauss:window_list");
  for (int i = 0; i < nwindow; i++)
    for (int j = 0; j < nvalues; j++)
      window_list[i][j] = 0.0;
  
  
  // the last N results, taking the longest delay, produce 2 outputs per value
  int delaymax = 0;
  for (int i = 0; i < ndelay; i++) {
    if (delays[i] > delaymax)
      delaymax = delays[i];
  }
  nresult = delaymax + 1;
  memory->create(result_list, nresult, nvalues*2, "ave/gauss:result_list");
  for (int i = 0; i < nresult; i++)
    for (int j = 0; j < nvalues*2; j++)
      result_list[i][j] = 0.0;
  
  // this fix produces a global vector and array
  vector_flag = 1;
  size_vector = nvalues*2;
  extvector = 0;
  array_flag = 1;
  size_array_rows = nvalues*2;
  size_array_cols = ndelay;
  extarray = 0;

  // initializations

  iwindow = window_limit = 0;
  iresult = 0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = update->ntimestep;
  nvalid_last = -1;
  nfull_next = update->ntimestep;
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveGauss::~FixAveGauss()
{
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] varindex;
  
  delete [] delays;

  memory->destroy(window_list);
  memory->destroy(result_list);
}

/* ---------------------------------------------------------------------- */

int FixAveGauss::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    int ivariable = input->variable->find(ids[i]);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for fix ave/gauss does not exist");
    varindex[i] = ivariable;
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    nvalid = update->ntimestep;
    nfull_next = update->ntimestep;
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveGauss::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::end_of_step()
{
  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/gauss");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;
  
  append_values(ntimestep);
  
  if (ntimestep == nfull_next) {
    update_results(ntimestep);
    nfull_next = ntimestep + nfull;
  }
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::append_values(bigint ntimestep)
{
  int i,m;
  double scalar;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  double *vector(window_list[iwindow]);

  for (i = 0; i < nvalues; i++) {
    m = varindex[i];
    scalar = input->variable->compute_equal(m);
    vector[i] = scalar;
  }
  
  iwindow += 1;
  if (iwindow >= nwindow) {
    window_limit = 1;
    iwindow = 0;
  }

  nvalid += 1;
  modify->addstep_compute(nvalid);
}

/* ---------------------------------------------------------------------- */

void FixAveGauss::update_results(bigint ntimestep)
{
  int count = nwindow;
  if (window_limit == 0)
    count = iwindow;

  if (count<1) return;

  double invcount = 1.0 / count;
  double *result = result_list[iresult];
  for (int i = 0; i < nvalues; i++)
    result[i*2] = result[i*2+1] = 0.0;

  // first pass: mean
  for (int j = 0; j < count; j++) {
    double *row = window_list[j];
    for (int i = 0; i < nvalues; i++)
      result[i*2] += row[i] * invcount;
  }
  
  // second pass: de-biased variance
  for (int j = 0; j < count; j++) {
    double *row = window_list[j];
    for (int i = 0; i < nvalues; i++) {
      double x = row[i] - result[i*2];
      result[i*2+1] += x * x * invcount;
    }
  }

  // return as stddev
  for (int i = 0; i < nvalues; i++)
    result[i*2+1] = sqrt(result[i*2+1]);
    
  iresult += 1;
  if (iresult >= nresult)
    iresult = 0;
}


/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveGauss::compute_scalar()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveGauss::compute_vector(int i)
{
  return compute_array(i, 0);
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveGauss::compute_array(int i, int j)
{
  if (i >= nvalues*2) return 0.0;
  if (j >= ndelay) return 0.0;
  int row = (iresult - 1 - delays[j] + nresult) % nresult;
  if (row >= nresult) return 0.0;
  return result_list[row][i];
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveGauss::options(int iarg, int narg, char **arg)
{
  // option defaults

  nfull = 1;
  nwindow = 100;
  ndelay = 0;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"window") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/gauss command");
      nwindow = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/gauss command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/gauss command");
      nfull = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nfull <= 0) error->all(FLERR,"Illegal fix ave/gauss command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/gauss command");
      int delay = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (ndelay < 0) error->all(FLERR,"Illegal fix ave/gauss command");
      if (ndelay+1==max_delays) error->all(FLERR,"Illegal fix ave/gauss command");
      delays[ndelay] = delay;
      ndelay += 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/gauss command");
  }
  
  if (ndelay == 0) {
    ndelay = 1;
    delays[0] = 0;
  }
}
