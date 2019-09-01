#ifndef LMP_DEBUGMATH_H
#define LMP_DEBUGMATH_H

#ifndef _GNU_SOURCE
  #define LAMMPS__GNU_SOURCE
  #define _GNU_SOURCE
#endif
#include <fenv.h>

// https://stackoverflow.com/a/2949452
static inline void debug_math_enable()
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
}

#define LAMMPS_DEBUGMATH_ENABLE

class DBG_FP_Exceptions
{
public:
  static const int COMMON = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
public: 
  DBG_FP_Exceptions(int enableexcepts)
  {
#if defined(LAMMPS_DEBUGMATH_ENABLE)
    this->saved = fegetexcept();
    if (enableexcepts) {
      this->enable(enableexcepts);
    }
#endif
  }
  
  ~DBG_FP_Exceptions() {
#if defined(LAMMPS_DEBUGMATH_ENABLE)
    fedisableexcept(FE_ALL_EXCEPT);
    feenableexcept(this->saved);
#endif
  }
  
  void enable(int excepts) {
#if defined(LAMMPS_DEBUGMATH_ENABLE)
    feenableexcept(excepts);
#endif
  }
  
  void disable(int excepts) {
#if defined(LAMMPS_DEBUGMATH_ENABLE)
    fedisableexcept(excepts);
#endif
  }
private:
   int saved;
};


#ifdef LAMMPS__GNU_SOURCE
  #undef _GNU_SOURCE
  #undef LAMMPS__GNU_SOURCE
#endif

#endif
