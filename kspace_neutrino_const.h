#ifndef KSPACE_NEUTRINO_CONSTANTS
#define KSPACE_NEUTRINO_CONSTANTS
/*This file defines constants specific to the neutrinos.*/
/* for three massive neutrino species:
 * Could be made configurable at some point
 * Neutrino masses are in eV*/
#define NUSPECIES 3

/*Speed of light in cm/s: in allvars.h this is called 'C'*/
#define  LIGHTCGS           2.99792458e10
#define  GRAVITY     6.67408e-8 /*Newton's constant in cgs*/
#define  HUBBLE          3.24077929e-18	/* 100 km/s in h/sec */
#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/

#define FLOAT   1e-6            /*Floating point accuracy*/
/*Number of bins in integrations*/
#define GSL_VAL 200

/*M_PI is technically a GNU extension, so define it here just in case*/
#ifndef M_PI
#define M_PI 3.14159265358979323846L
#endif


#endif
