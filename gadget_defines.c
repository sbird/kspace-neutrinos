/* These functions will normally be defined by gadget. However, for the purposes of
 * the standalone tests, we want our own definitions*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gadget_defines.h"
#include "omega_nu_single.h"

//Forward define terminate, because we'll need it.
void terminate(const char * string)
{
    fprintf(stderr, "Error: %s\n",string);
    exit(1);
}

void * mymalloc(const char * string, size_t size)
{
    return malloc(size);
}

void myfree(void * ptr)
{
    free(ptr);
}

/* This is a bit goofy, because the hubble function in gadget uses global state.
 * So in the tests we need to set the global state with init_hubble_function before using it,
 * to keep the function definition the same.
 * We COULD just use an internal version, but then we would break compatibility if anyone has
 * used an odd hubble history for eg, some dark energy model. */
static _omega_nu omnu;
static double m_HubbleParam;

void init_hubble_function(const double MNu[], const double Omega0, const double a0, const double HubbleParam)
{
    init_omega_nu(&omnu, MNu, Omega0, a0, HubbleParam);
    m_HubbleParam = HubbleParam;

}

#define STEFAN_BOLTZMANN 5.670373e-5
#define OMEGAR (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*LIGHTCGS*LIGHTCGS*LIGHTCGS*HUBBLE*HUBBLE*m_HubbleParam*m_HubbleParam)*pow(T_CMB0,4))

double hubble_function(double a)
{
    if(!omnu.RhoNuTab[0]) {
        terminate("init_hubble_function was not called in test suite before this!\n");
    }
    /* Matter + Lambda: neglect curvature*/
    double omega_tot = omnu.Omega0/pow(a,3) + (1-omnu.Omega0);
    /*Neutrinos*/
    omega_tot += get_omega_nu(&omnu, a) - get_omega_nu(&omnu, 1)/pow(a,3);
    /*Radiation*/
    omega_tot += OMEGAR/pow(a,4);
    return m_HubbleParam * HUBBLE * sqrt(omega_tot);
}
