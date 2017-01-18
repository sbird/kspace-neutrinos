#ifndef KSPACE_NEUTRINO_CONSTANTS
#define KSPACE_NEUTRINO_CONSTANTS
/**\file
 * Defines constants specific to the neutrinos.*/

/** Number of massive neutrino species: 3
 * Could be made configurable at some point
 * Neutrino masses are in eV*/
#define NUSPECIES 3

/** Speed of light in cm/s: in allvars.h this is called 'C'*/
#define  LIGHTCGS           2.99792458e10
/** The Boltzmann constant in units of eV/K*/
#define BOLEVK 8.61734e-5

/** Floating point accuracy*/
#define FLOAT_ACC   1e-6
/** Number of bins in integrations*/
#define GSL_VAL 200

/** M_PI is technically a GNU extension, so define it here just in case*/
#ifndef M_PI
#define M_PI 3.14159265358979323846L
#endif


#endif
