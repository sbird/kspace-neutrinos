/*This file contains functions which need to be called from the PM code in Gadget.
 add_nu_power_to_rhogrid is the main public function. */

#ifdef KSPACE_NEUTRINOS_2
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kspace_neutrinos_2.h"
#include "kspace_neutrino_const.h"
#include "gadget_defines.h"
#include "transfer_init.h"
#include "delta_tot_table.h"
#include "delta_pow.h"
#include "powerspectrum.h"

//Global variables that need to be set from a parameter file
struct __kspace_params {
  char	KspaceTransferFunction[500];
  double TimeTransfer;
  double OmegaBaryonCAMB;
  double InputSpectrum_UnitLength_in_cm;
  double MNu[NUSPECIES];
#if defined HYBRID_NEUTRINOS
    /*Critical velocity above which to treat neutrinos with particles.
    Note this is unperturbed velocity *TODAY*
    To get velocity at redshift z, multiply by (1+z)*/
    double vcrit;
    //Time at which to turn on the particle neutrinos.
    //Ultimately we want something better than this.
    double nu_crit_time;
#endif
} kspace_params;

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

static _transfer_init_table transfer_init;

static _delta_tot_table delta_tot_table;

static _omega_nu omeganu_table;

void broadcast_transfer_table(_transfer_init_table *t_init, int ThisTask)
{
  MPI_Bcast(&(t_init->NPowerTable), 1,MPI_INT,0,MYMPI_COMM_WORLD);
  /*Allocate the memory unless we are on task 0, in which case it is already allocated*/
  if(ThisTask!=0)
    t_init->logk = (double *) mymalloc("Transfer_functions", 2*t_init->NPowerTable* sizeof(double));
  t_init->T_nu=t_init->logk+t_init->NPowerTable;
  /*Broadcast the arrays*/
  MPI_Bcast(t_init->logk,2*(t_init->NPowerTable),MPI_DOUBLE,0,MYMPI_COMM_WORLD);
}

/** This function loads the initial transfer functions from CAMB transfer files.
 * One processor 0 it reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu_init and friends and has length NPowerTable is then broadcast to all processors.
 * Then, on all processors, it allocates memory for delta_tot_table.*/
void allocate_kspace_memory(const int nk_in, const int ThisTask,const double BoxSize, const double UnitLength_in_cm, const double Omega0, const double HubbleParam)
{
  if(omeganu_table.RhoNuTab[0] == 0)
      init_omega_nu(&omeganu_table, kspace_params.MNu, Omega0, kspace_params.TimeTransfer, HubbleParam);
  /*We only need this for initialising delta_tot later.
   * ThisTask is needed so we only read the transfer functions on task 0, serialising disc access.*/
  if(ThisTask==0)
    allocate_transfer_init_table(&transfer_init, nk_in, BoxSize, UnitLength_in_cm, kspace_params.InputSpectrum_UnitLength_in_cm, kspace_params.OmegaBaryonCAMB, kspace_params.KspaceTransferFunction, &omeganu_table);
  /*Broadcast data to other processors*/
  broadcast_transfer_table(&transfer_init, ThisTask);
  /*Set the private copy of the task in delta_tot_table*/
  delta_tot_table.ThisTask = ThisTask;
  allocate_delta_tot_table(&delta_tot_table, nk_in, kspace_params.TimeTransfer, ThisTask);
  /*Check that if we are restarting from a snapshot, we successfully read a table*/
/*   if(fabs(kspace_vars.TimeBegin - d_tot->TimeTransfer) >1e-4 && (!d_tot->ia)) */
/*      terminate("Transfer function not at the same time as simulation start (are you restarting from a snapshot?) and could not read delta_tot table\n"); */
}

#define TARGETBINS 300              /* Number of bins in the smoothed power spectrum*/

/* This function adds the neutrino power spectrum to the
 * density grid. It calls the gadget power spectrum routines, which output their
 * results in SumPower (the total power in all modes in a bin),
 * SumPowerUncorrected (the total power, minus a correction for the bin width)
 * and CountModes (the number of modes per bin).
 * SumPower[1] is from box scale to grid cell scale,
 * SumPower[0] is the folded power on smaller scales.
 * It also touches fft_of_rhogrid, which is the fourier transformed density grid.
 */
void add_nu_power_to_rhogrid(int save, const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int PMGRID, int ThisTask, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir, const double total_mass)
{
  /*Some of the neutrinos will be relativistic at early times. However, the transfer function for the massless neutrinos
   * is very similar to the transfer function for the massive neutrinos, so treat them the same*/
  const double OmegaNua3 = get_omega_nu(&omeganu_table, Time)*pow(Time,3);
  /*kspace_prefac = M_nu / M_cdm */
  const double kspace_prefac = OmegaNua3/(omeganu_table.Omega0-get_omega_nu(&omeganu_table, 1));
  int i,x,y,z;
  /*Calculate the power for kspace neutrinos*/
  double delta_cdm_curr[TARGETBINS];  /* (square root of) the power spectrum.*/
  double keff[TARGETBINS];      /* (binned) k values for the power spectrum*/
  double delta_nu_curr[TARGETBINS];  /*The square root of the neutrino power spectrum*/
  long long int count[TARGETBINS];
  /*We calculate the power spectrum at every timestep
   * because we need it as input to the neutrino power spectrum.
   * This function stores the total power*no. modes.*/
  total_powerspectrum(PMGRID, fft_of_rhogrid, TARGETBINS, slabstart_y, nslab_y, delta_cdm_curr, count, keff, total_mass);
  /*Get delta_cdm_curr , which is P(k)^1/2, and convert P(k) to physical units. */
  const double scale=pow(2*M_PI/BoxSize,3);
  for(i=0;i<TARGETBINS;i++){
      delta_cdm_curr[i] = sqrt(delta_cdm_curr[i]/scale);
  }
  /*This sets up P_nu_curr.*/
  get_delta_nu_update(&delta_tot_table, Time, TARGETBINS, keff, delta_cdm_curr,  delta_nu_curr);
  for(i=0;i<TARGETBINS;i++){
          if(isnan(delta_nu_curr[i]) || delta_nu_curr[i] < 0){
                  char err[300];
                  snprintf(err,300,"delta_nu_curr=%g z=%d SmoothPow=%g kk=%g\n",delta_nu_curr[i],i,delta_cdm_curr[i],keff[i]);
                  terminate(err);
          }
  }
  /*Sets up the interpolation for get_neutrino_powerspec*/
  _delta_pow d_pow;
  /*We want to interpolate in log space*/
  for(i=0;i<TARGETBINS;i++){
      keff[i] = log(keff[i]);
  }
  init_delta_pow(&d_pow, keff, delta_nu_curr, delta_cdm_curr, TARGETBINS);
  /*Add P_nu to fft_of_rhgrid*/
  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
        {
          double kx,ky,kz,k2,smth;
          int ip;
          kx = x > PMGRID/2 ? x-PMGRID : x;
          ky = y > PMGRID/2 ? y-PMGRID : y;
          kz = z > PMGRID/2 ? z-PMGRID : z;

          k2 = kx*kx + ky*ky + kz*kz;
          if(k2 <= 0)
              continue;
          /*Change the units of k to match those of logkk*/
          k2=log(sqrt(k2)*2*M_PI/BoxSize);
          /* Note get_neutrino_powerspec returns delta_nu / P_cdm^1/2, which is dimensionless.
           * We have delta_t = (M_cdm+M_nu)*delta_cdm (1-f_nu + f_nu (delta_nu / delta_cdm)^1/2)
           * which gives the right power spectrum, once we divide by
           * M_cdm +M_nu in powerspec*/
          smth=(1+kspace_prefac * get_dnudcdm_powerspec(&d_pow, k2));
          ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
          fft_of_rhogrid[ip].re *= smth;
          fft_of_rhogrid[ip].im *= smth;
        }
  /*If this is being called to save all particle types, save a file with the neutrino power spectrum as well.*/
  if(save){
            save_nu_power(&d_pow, Time, snapnum, OutputDir);
  }
  free_d_pow(&d_pow);
  return;
}

#endif //KSPACE_NEUTRINOS_2
