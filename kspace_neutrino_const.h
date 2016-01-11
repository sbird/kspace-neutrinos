#ifndef KSPACE_NEUTRINO_CONSTANTS
#define KSPACE_NEUTRINO_CONSTANTS
/*This file defines constants specific to the neutrinos.
 * It is separate from gadget_defines.h in that those constants are also defined in gadget*/

/* for three massive neutrino species:
 * Could be made configurable at some point
 * Neutrino masses are in eV*/
#define NUSPECIES 3

#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/

#define FLOAT   1e-6            /*Floating point accuracy*/
/*Number of bins in integrations*/
#define GSL_VAL 200

#endif
