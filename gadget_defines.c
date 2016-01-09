/* These functions will normally be defined by gadget. However, for the purposes of
 * the standalone tests, we want our own definitions*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gadget_defines.h"
#include "omega_nu_single.h"

/*These functions need bodies; normally this is provided by gadget*/
void * mymalloc_fullinfo(const char * string, size_t size, const char *func, const char *file, int line)
{
    return malloc(size);
}

void myfree_fullinfo(void * ptr, const char *func, const char *file, int line)
{
    free(ptr);
}

/* This is a bit goofy, because the hubble function in gadget uses global state.
 * So in the tests we need to set the global state with init_hubble_function before using it,
 * to keep the function definition the same.
 * We COULD just use an internal version, but then we would break compatibility if anyone has
 * used an odd hubble history for eg, some dark energy model. */
static _omega_nu * m_omnu;
static double Omega_nonu;
static double OmegaLambda;
static double m_Hubble;

void init_hubble_function(_omega_nu * omnu, const double Omega0, const double UnitTime_in_s)
{
    m_omnu = omnu;
    Omega_nonu = Omega0 - get_omega_nu(omnu, 1);
    OmegaLambda = 1 - Omega0;
    m_Hubble = HUBBLE * UnitTime_in_s;
}

double hubble_function(double a)
{
    if(!m_omnu) {
        terminate("init_hubble_function was not called in test suite before this!\n");
    }
    /* Matter + Lambda: neglect curvature*/
    double omega_tot = Omega_nonu/pow(a,3) + OmegaLambda;
    /*Neutrinos*/
    omega_tot += get_omega_nu(m_omnu, a);
    /*Radiation*/
    omega_tot += get_omegag(m_omnu, a);
    return m_Hubble * sqrt(omega_tot);
}
