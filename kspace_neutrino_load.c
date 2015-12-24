/* File to contain various functions for saving and resuming kspace neutrinos*/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "kspace_neutrinos_vars.h"

//These are defined in begrun.c
#define STRING 2
#define REAL 1

//Setup the config files to load the needed variables
int set_kspace_vars(char * tag[], void *addr[], int id [], int nt)
{
      strcpy(tag[nt], "KspaceTransferFunction");
      addr[nt] = kspace_params.KspaceTransferFunction;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimeTransfer");
      addr[nt] = &kspace_params.TimeTransfer;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryonCAMB");
      addr[nt] = &kspace_params.OmegaBaryonCAMB;
      id[nt++] = REAL;

      strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
      addr[nt] = &kspace_params.InputSpectrum_UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "MNue");
      addr[nt] = &(kspace_params.MNu[0]);
      id[nt++] = REAL;
      strcpy(tag[nt], "MNum");
      addr[nt] = &(kspace_params.MNu[1]);
      id[nt++] = REAL;
      strcpy(tag[nt], "MNut");
      addr[nt] = &(kspace_params.MNu[2]);
      id[nt++] = REAL;
#if defined HYBRID_NEUTRINOS
    strcpy(tag[nt], "VCRIT");
    addr[nt] = &(kspace_params.vcrit);
    id[nt++] = REAL;
    strcpy(tag[nt], "NuPartTime");
    addr[nt] = &(kspace_params.nu_crit_time);
    id[nt++] = REAL;
#endif
    return nt;
}
