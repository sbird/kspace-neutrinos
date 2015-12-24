#ifndef KSPACE_NEUTRINO_VARS
#define KSPACE_NEUTRINO_VARS

#include <stdlib.h>
#include "kspace_neutrino_const.h"

//Global variables that need to be set from a parameter file
struct __kspace_params {
  char	KspaceTransferFunction[500];
  double TimeTransfer;
  double OmegaBaryonCAMB;
  double InputSpectrum_UnitLength_in_cm;
  double MNu[NUSPECIES];
#if defined HYBRID_NEUTRINOS
    /*Two parameters for the hybrid neutrinos.
    If this is true, then we proceed using the analytic method for all neutrinos.
    If this is false, then we cut off the analytic method at q < qcrit (specified using vcrit, below) and use
    particles for the slower neutrinos.*/
    int slow_neutrinos_analytic;
    /*Critical velocity above which to treat neutrinos with particles.
    Note this is unperturbed velocity *TODAY*
    To get velocity at redshift z, multiply by (1+z)*/
    double vcrit;
    //Time at which to turn on the particle neutrinos.
    //Ultimately we want something better than this.
    double nu_crit_time;
#endif
} kspace_params;

//Function which sets the above variables
int set_kspace_vars(char * tag[], void *addr[], int id [], int nt);

#endif
