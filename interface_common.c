/* This file contains functions which need to be called from the PM code in Gadget.
 * Everything required to interface with the main code should be in this file.
 * All MPI communication is also done here.
 * add_nu_power_to_rhogrid is the main public function. */
#include "interface_common.h"

#include <stdio.h>
/* #include <string.h> */
#include <math.h>
#include "kspace_neutrino_const.h"
#include "gadget_defines.h"
#include "omega_nu_single.h"
#include "transfer_init.h"
#include "delta_tot_table.h"
#include "delta_pow.h"

/*Global neutrino module parameters*/
struct __kspace_params kspace_params;

static _transfer_init_table transfer_init;

static _delta_tot_table delta_tot_table;

static _omega_nu omeganu_table;

/*Compute the matter density in neutrinos*/
double OmegaNu(double a)
{
    return get_omega_nu(&omeganu_table, a);
}

/* Compute the matter density in neutrinos, 
 * excluding density in particles.*/
double OmegaNu_nopart(double a)
{
    return get_omega_nu_nopart(&omeganu_table, a);
}

void save_nu_state(char * savedir)
{
    if(delta_tot_table.ThisTask == 0)
        save_all_nu_state(&delta_tot_table, savedir);
}

void broadcast_transfer_table(_transfer_init_table *t_init, int ThisTask, MPI_Comm MYMPI_COMM_WORLD)
{
  MPI_Bcast(&(t_init->NPowerTable), 1,MPI_INT,0,MYMPI_COMM_WORLD);
  /*Allocate the memory unless we are on task 0, in which case it is already allocated*/
  if(ThisTask != 0)
    t_init->logk = (double *) mymalloc("Transfer_functions", 2*t_init->NPowerTable* sizeof(double));
  t_init->T_nu=t_init->logk+t_init->NPowerTable;
  /*Broadcast the arrays*/
  MPI_Bcast(t_init->logk,2*(t_init->NPowerTable),MPI_DOUBLE,0,MYMPI_COMM_WORLD);
}

void broadcast_delta_tot_table(_delta_tot_table *d_tot, const int nk_in, MPI_Comm MYMPI_COMM_WORLD)
{
  /*Broadcast array sizes*/
  MPI_Bcast(&(d_tot->ia), 1,MPI_INT,0,MYMPI_COMM_WORLD);
  if(d_tot->ia > 0) {
      MPI_Bcast(&(d_tot->nk), 1,MPI_INT,0,MYMPI_COMM_WORLD);
      /*Broadcast data for scalefact and delta_tot, Delta_tot is allocated as the same block of memory as scalefact.
        Not all this memory will actually have been used, but it is easiest to bcast all of it.*/
      MPI_Bcast(d_tot->scalefact,d_tot->namax*(nk_in+1),MPI_DOUBLE,0,MYMPI_COMM_WORLD);
  }
}

/** This function loads the initial transfer functions from CAMB transfer files.
 * One processor 0 it reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu_init and friends and has length NPowerTable is then broadcast to all processors.
 * Then, on all processors, it allocates memory for delta_tot_table.
 * This must be called *EARLY*, before OmegaNu, just after the parameters are read.
 * Arguments:
 * nk_in - number of bins desired in the neutrino power spectrum
 * ThisTask - MPI rank
 * BoxSize - size of box in internal units
 * UnitTime_in_s, UnitLength_in_cm - conversion factors from internal units to cgs.
 * Omega0 - total matter density (including massive neutrinos and baryons but not including radiation)
 * tcmb0 - present-day CMB temperature
 * snapdir - snapshot directory to try to read state and resume from
 * TimeMax - Final redshift desired, sets number of output redshift bins
 * MYMPI_COMM_WORLD - MPI  communicator
 * Global state used: omeganu_table, delta_tot_table, transfer_init */
void allocate_kspace_memory(const int nk_in, const int ThisTask, const double BoxSize, const double UnitTime_in_s, const double UnitLength_in_cm, const double Omega0, const double HubbleParam, const double tcmb0, const char * snapdir, const double TimeMax, MPI_Comm MYMPI_COMM_WORLD)
{
  /*First make sure kspace_params is propagated to all processors*/
  MPI_Bcast(&kspace_params,sizeof(kspace_params),MPI_BYTE,0,MYMPI_COMM_WORLD);
  /*Now initialise the background*/
  init_omega_nu(&omeganu_table, kspace_params.MNu, kspace_params.TimeTransfer, HubbleParam, tcmb0);
  if(kspace_params.hybrid_neutrinos_on)
    init_hybrid_nu(&omeganu_table.hybnu, kspace_params.MNu, kspace_params.vcrit, LIGHTCGS * UnitTime_in_s/UnitLength_in_cm, kspace_params.nu_crit_time, omeganu_table.kBtnu);
  /*We only need this for initialising delta_tot later.
   * ThisTask is needed so we only read the transfer functions on task 0, serialising disc access.*/
  if(ThisTask==0) {
    allocate_transfer_init_table(&transfer_init, BoxSize, UnitLength_in_cm, kspace_params.InputSpectrum_UnitLength_in_cm, kspace_params.OmegaBaryonCAMB, get_omega_nu(&omeganu_table, 1), Omega0, kspace_params.KspaceTransferFunction);
  }
  /*Broadcast data to other processors*/
  broadcast_transfer_table(&transfer_init, ThisTask, MYMPI_COMM_WORLD);
  /*Set the private copy of the task in delta_tot_table*/
  delta_tot_table.ThisTask = ThisTask;
  allocate_delta_tot_table(&delta_tot_table, nk_in, kspace_params.TimeTransfer, TimeMax, Omega0, &omeganu_table, UnitTime_in_s, UnitLength_in_cm, 1);
  /*Read the saved data from a snapshot if present*/
  if(ThisTask==0) {
  	read_all_nu_state(&delta_tot_table, snapdir);
  }
  /*Broadcast save-data to other processors*/
  broadcast_delta_tot_table(&delta_tot_table, nk_in, MYMPI_COMM_WORLD);
}

/* This function calls the integrator to compute the neutrino power spectrum,
 * taking as input a pre-computed matter power spectrum, assumed to have the same units as stored in transfer_init.
 * neutrino power spectrum is stored in _delta_pow and returned.
 * Memory allocated here must be freed later.
 * Arguments:
 * Time - scale factor, a.
 * keff_in - k values for each power bin. Has units of UnitLength_in_cm passed to transfer_init
 * P_cdm - Normalised matter power spectrum. Has units of UnitLength_in_cm passed to transfer_init
 * Nmodes - number of modes in each bin. Used only to see if bin is nonempty.
 * MYMPI_COMM_WORLD - MPI communicator to use
 * Global state used: delta_tot_table, transfer_init, omeganu_table
 * Returns: _delta_pow, containing delta_nu/delta_cdm*/
_delta_pow compute_neutrino_power_from_cdm(const double Time, const double keff_in[], const double P_cdm[], const long long int Nmodes[], const int nk_in, MPI_Comm MYMPI_COMM_WORLD)
{
  int i;
  /*Allocate memory to copy power spectrum to. The caller may free the input before we are done.*/
  double * delta_cdm_curr = mymalloc("temp_power_spectrum", 3*nk_in*sizeof(double));
  /*The square root of the neutrino power spectrum*/
  double * delta_nu_curr = delta_cdm_curr+nk_in;
  /* (binned) k values for the power spectrum*/
  double * keff = delta_cdm_curr+2*nk_in;
  if(!delta_cdm_curr || !delta_nu_curr || !keff)
      terminate("Could not allocate temporary memory for power spectra\n");
  /*Get delta_cdm_curr , which is P(k)^1/2, and skip bins with zero modes. */
  int nk_nonzero = 0;
  for(i=0;i<nk_in;i++){
      if (Nmodes[i] == 0)
          continue;
      delta_cdm_curr[nk_nonzero] = sqrt(P_cdm[i]);
      keff[nk_nonzero] = keff_in[i];
      nk_nonzero++;
  }
  /*This sets up P_nu_curr.*/
  get_delta_nu_update(&delta_tot_table, Time, nk_nonzero, keff, delta_cdm_curr,  delta_nu_curr, &transfer_init);

  MPI_Barrier(MYMPI_COMM_WORLD);
  if(delta_tot_table.ThisTask==0)
	printf("Done get_delta_nu_update on all processors\n");
  /*Sets up the interpolation for get_neutrino_powerspec*/
  _delta_pow d_pow;
  /*We want to interpolate in log space*/
  for(i=0;i<nk_in;i++){
      keff[i] = log(keff[i]);
  }
  /*kspace_prefac = M_nu (analytic) / M_particles */
  const double OmegaNu_nop = get_omega_nu_nopart(&omeganu_table, Time);
  /* Note if (hybrid) neutrino particles are off, this is zero.
   * We cannot just use OmegaNu(1) as we need to know
   * whether hybrid neutrinos are on at this redshift.*/
  const double omega_hybrid = get_omega_nu(&omeganu_table, Time) - OmegaNu_nop;
  /* Omega0 - Omega in neutrinos + Omega in particle neutrinos = Omega in particles*/
  const double kspace_prefac = OmegaNu_nop/(delta_tot_table.Omeganonu/pow(Time,3) + omega_hybrid);
  init_delta_pow(&d_pow, keff, delta_nu_curr, delta_cdm_curr, nk_in, kspace_prefac);
  return d_pow;
}
