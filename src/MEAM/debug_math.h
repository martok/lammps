#ifndef LMP_DEBUGMATH_H
#define LMP_DEBUGMATH_H

#ifndef _GNU_SOURCE
  #define LAMMPS__GNU_SOURCE
  #define _GNU_SOURCE
#endif
#include <fenv.h>

#if defined(__APPLE__) && defined(__MACH__)

// Public domain polyfill for feenableexcept on OS X
// http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c

inline int fegetexcept()
{
    static fenv_t fenv;

    if (fegetenv(&fenv)) {
        return -1;
    }
    return ~fenv.__control & FE_ALL_EXCEPT;
}

inline int feenableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // previous masks
    unsigned int old_excepts;

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr &= ~(new_excepts << 7);

    return fesetenv(&fenv) ? -1 : old_excepts;
}

inline int fedisableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // all previous masks
    unsigned int old_excepts;

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // mask
    fenv.__control |= new_excepts;
    fenv.__mxcsr |= new_excepts << 7;

    return fesetenv(&fenv) ? -1 : old_excepts;
}

#endif

#if defined(_MSC_VER)

// Polyfill for feenableexcept on MSVC

inline int fegetexcept()
{
    unsigned int cw;
    _controlfp_s(&cw, 0, 0);
    return ~cw & _MCW_EM;
}

inline int feenableexcept(unsigned int excepts)
{
    unsigned int new_excepts = excepts & _MCW_EM;

    unsigned int old_excepts, cw;
    _controlfp_s(&old_excepts, 0, 0);
    if (excepts)
        _controlfp_s(&cw, old_excepts & ~new_excepts, _MCW_EM);
    return ~old_excepts & _MCW_EM;
}

inline int fedisableexcept(unsigned int excepts)
{
    unsigned int new_excepts = excepts & _MCW_EM;
    unsigned int old_excepts, cw;
    _controlfp_s(&old_excepts, 0, 0);
    if (excepts)
        _controlfp_s(&cw, old_excepts | new_excepts, _MCW_EM);
    return ~old_excepts & _MCW_EM;
}

#endif


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
