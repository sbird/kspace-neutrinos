#ifndef KSPACE_NEUTRINO_VARS
#define KSPACE_NEUTRINO_VARS

#include <stdlib.h>
  /* for three massive neutrino species:
   * Could be made configurable at some point
   * Neutrino masses are in eV*/
  #define NUSPECIES 3

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

//These are global variables that come from the All structure.
struct __kspace_vars {
  double OmegaNu;
  double HubbleParam;
  double UnitLength_in_cm;
  double UnitTime_in_s;
  double Hubble;
  double BoxSize;
  double TimeBegin;
  double Omega0;
  double TimeMax;
} kspace_vars;

// static int ThisTask=0;

/* Note Omega0, the total non-relativistic matter, includes neutrinos (and radiation). */
#define H0      (kspace_vars.HubbleParam)                /* H0 in units of H100*/

/*Light speed in internal units. C is defined in allvars.h to be lightspeed in cm/s*/
#define LIGHT (C*kspace_vars.UnitTime_in_s/kspace_vars.UnitLength_in_cm)

#define H100   kspace_vars.Hubble /* 100 km/s/Mpc in units of 1/UnitTime. */
/*The time at which we first start our integrator:
 * NOTE! This is not All.TimeBegin, but the time of the transfer function file,
 * so that we can support restarting from snapshots.*/
#define A0      (kspace_params.TimeTransfer)

//Function which sets the above variables
int set_kspace_vars(char * tag[], void *addr[], int id [], int nt);

//Forward define terminate, because we'll need it.
void terminate(const char *);
#ifndef mymalloc
#define mymalloc(x,y) malloc(y)
#endif

#ifndef myfree
#define myfree(x) free(x)
#endif

#endif
