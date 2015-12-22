/*This file contains routines for computing OmegaNu, in the full relativistic approximation.
This is needed for adding massive neutrinos to any simulation.

double OmegaNu(double) is the only externally accessible function.
*/

#include "kspace_neutrinos_func.h"
#include "kspace_neutrinos_vars.h"
#include "kspace_neutrino_const.h"
#include <math.h>
#include "omega_nu_single.h"

/*Pointers to the array of structures we use to store rho_nu*/
static _rho_nu_single * RhoNuTab[NUSPECIES];

/* Return the total matter density in neutrinos.
 * rho_nu and friends are not externally callable*/
double OmegaNu(double a)
{
        double rhonu=0;
        int mi;
        double OmegaNu_one[NUSPECIES];
        for(mi=0; mi<NUSPECIES; mi++){
             int mmi;
             for(mmi=0; mmi<mi; mmi++){
                 if(fabs(kspace_params.MNu[mi] -kspace_params.MNu[mmi]) < FLOAT)
                   break;
             }
             if(mmi==mi) {
                 /*Allocate the table if we didn't do this already*/
                 if(!RhoNuTab[mmi]) {
                    RhoNuTab[mmi] = mymalloc("RhoNuTab", sizeof(_rho_nu_single));
                    rho_nu_init(RhoNuTab[mmi], kspace_params.TimeTransfer, kspace_params.MNu[mmi]);
                 }
                 OmegaNu_one[mi]=omega_nu_single(RhoNuTab[mmi], a);
             }
            rhonu+=OmegaNu_one[mmi];
        }
        return rhonu;
}
