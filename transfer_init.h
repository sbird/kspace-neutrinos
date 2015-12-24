#ifndef TRANSFER_INIT_H
#define TRANSFER_INIT_H

#include "omega_nu_single.h"
/* Structure to store the initial transfer functions from CAMB.
 * We store transfer functions because we want to use the
 * CDM + Baryon total matter power spectrum from the
 * first timestep of Gadget, so that possible Rayleigh scattering
 * in the initial conditions is included in the neutrino and radiation components. */
struct _transfer_init_table {
    int NPowerTable;
    double *logk;
    double *T_nu;
};
typedef struct _transfer_init_table _transfer_init_table;

void allocate_transfer_init_table(_transfer_init_table *t_init, int nk_in, const double BoxSize, const double UnitLength_in_cm, const double InputSpectrum_UnitLength_in_cm, const double OmegaBaryonCAMB, _omega_nu * omnu);

#endif //TRANSFER_INIT_H
