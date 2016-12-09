#ifndef KSPACE_NEUTRINO_CONSTANTS
#define KSPACE_NEUTRINO_CONSTANTS
/*This file defines constants specific to the neutrinos.*/
/* for three massive neutrino species:
 * Could be made configurable at some point
 * Neutrino masses are in eV*/
#define NUSPECIES 3

/*Speed of light in cm/s: in allvars.h this is called 'C'*/
#define  LIGHTCGS           2.99792458e10
#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/

#define FLOAT   1e-6            /*Floating point accuracy*/
/*Number of bins in integrations*/
#define GSL_VAL 200

/*M_PI is technically a GNU extension, so define it here just in case*/
#ifndef M_PI
#define M_PI 3.14159265358979323846L
#endif


#endif
