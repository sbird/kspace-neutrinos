/*This file contains routines for computing OmegaNu, in the full relativistic approximation.
This is needed for adding massive neutrinos to any simulation.

double OmegaNu(double) is the only externally accessible function.
*/

#include "kspace_neutrinos_func.h"
#include "kspace_neutrinos_vars.h"
#include "kspace_neutrino_const.h"
#include <math.h>
#include "omega_nu_single.h"
