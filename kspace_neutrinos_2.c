/* This file contains functions which need to be called from the PM code in Gadget.
 * Everything required to interface with the main code should be in this file.
 * All MPI communication is also done here.
 * add_nu_power_to_rhogrid is the main public function. */
#include "kspace_neutrinos_2.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kspace_neutrino_const.h"
#include "gadget_defines.h"
#include "omega_nu_single.h"
#include "transfer_init.h"
#include "delta_tot_table.h"
#include "delta_pow.h"
#include "powerspectrum.h"

/*Global neutrino module parameters*/
struct __kspace_params kspace_params;

/*Setup the config files to load the needed variables.
 * This is an example config file reader specific to P-Gadget3.*/
int set_kspace_vars(char tag[][50], void *addr[], int id [], int nt)
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
      strcpy(tag[nt], "HybridNeutrinosOn");
      addr[nt] = &(kspace_params.hybrid_neutrinos_on);
      id[nt++] = INT;
      strcpy(tag[nt], "Vcrit");
      addr[nt] = &(kspace_params.vcrit);
      id[nt++] = REAL;
      strcpy(tag[nt], "NuPartTime");
      addr[nt] = &(kspace_params.nu_crit_time);
      id[nt++] = REAL;
      return nt;
}

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

/* This function calculates the matter power spectrum, then calls the integrator to compute the neutrino power spectrum,
 * which is stored in _delta_pow and returned.
 * Arguments:
 * Time - scale factor, a.
 * BoxSize - size of the box in internal units.
 * fft_of_rhogrid - Fourier transformed density grid.
 * pmgrid - size of one dimension of the density grid.
 * slabstart_y - for slab parallelized FFT routines, this is the start index of the FFT on this rank.
 * nslab_y - number of elements of the FFT on this rank.
 * MYMPI_COMM_WORLD - MPI communicator to use
 * Global state used: delta_tot_table, transfer_init, omeganu_table
 * Returns: _delta_pow, containing delta_nu/delta_cdm*/
_delta_pow compute_neutrino_power_spectrum(const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int pmgrid, int slabstart_y, int nslab_y, MPI_Comm MYMPI_COMM_WORLD)
{
  int i, nk_in;
  const int nk_allocated = delta_tot_table.nk_allocated;
  /*Calculate the power for kspace neutrinos*/
  /* (square root of) the power spectrum.*/
  double * delta_cdm_curr = mymalloc("temp_power_spectrum", 3*nk_allocated*sizeof(double));
  /*The square root of the neutrino power spectrum*/
  double * delta_nu_curr = delta_cdm_curr+nk_allocated;
  /* (binned) k values for the power spectrum*/
  double * keff = delta_cdm_curr+2*nk_allocated;
  long long int * count = mymalloc("temp_modecount", nk_allocated*sizeof(long long int));
  const double scale=pow(2*M_PI/BoxSize,3);
  if(!delta_cdm_curr || !delta_nu_curr || !keff || !count)
      terminate("Could not allocate temporary memory for power spectra\n");
  /*We calculate the power spectrum at every timestep
   * because we need it as input to the neutrino power spectrum.
   * This function stores the total power*no. modes.*/
  nk_in = total_powerspectrum(pmgrid, fft_of_rhogrid, nk_allocated, slabstart_y, nslab_y, delta_cdm_curr, count, keff, MYMPI_COMM_WORLD);
  /*Don't need count memory any more*/
  myfree(count);
  /*Get delta_cdm_curr , which is P(k)^1/2, and convert P(k) to physical units. */
  for(i=0;i<nk_in;i++){
      delta_cdm_curr[i] = sqrt(delta_cdm_curr[i]/scale);
      keff[i] *= (2*M_PI/BoxSize);
  }
  /*This sets up P_nu_curr.*/
  get_delta_nu_update(&delta_tot_table, Time, nk_in, keff, delta_cdm_curr,  delta_nu_curr, &transfer_init);

  MPI_Barrier(MYMPI_COMM_WORLD);
  if(delta_tot_table.ThisTask==0)
	printf("Done get_delta_nu_update on all processors\n");
  /*Sets up the interpolation for get_neutrino_powerspec*/
  _delta_pow d_pow;
  /*We want to interpolate in log space*/
  for(i=0;i<nk_in;i++){
      keff[i] = log(keff[i]);
  }
  /*Some of the neutrinos will be relativistic at early times. However, the transfer function for the massless neutrinos
   * is very similar to the transfer function for the massive neutrinos, so treat them the same*/
  const double OmegaNua3 = get_omega_nu_nopart(&omeganu_table, Time)*pow(Time,3);
  /*kspace_prefac = M_nu (analytic) / M_particles */
  const double OmegaNu_nop = get_omega_nu_nopart(&omeganu_table, Time);
  const double kspace_prefac = OmegaNu_nop*pow(Time,3)/(delta_tot_table.Omeganonu + get_omega_nu(&omeganu_table, Time) - OmegaNu_nop);
  init_delta_pow(&d_pow, keff, delta_nu_curr, delta_cdm_curr, nk_in, kspace_prefac);
  return d_pow;
}

/* This function adds the neutrino power spectrum to the
 * density grid. It calls the internal power spectrum routine and the neutrino integrator.
 * It then adds the neutrino power to fft_of_rhogrid, which is the fourier transformed density grid from the PM code.
 * Arguments:
 * Time - scale factor, a.
 * BoxSize - size of the box in internal units.
 * fft_of_rhogrid - Fourier transformed density grid.
 * pmgrid - size of one dimension of the density grid.
 * slabstart_y - for slab parallelized FFT routines, this is the start index of the FFT on this rank.
 * nslab_y - number of elements of the FFT on this rank.
 * snapnum - number of snapshot to save neutrino power spectrum as powerspec_nu_$(snapnum).txt
 * OutputDir - output directory for neutrino power spectrum.
 * MYMPI_COMM_WORLD - MPI communicator to use
 */
void add_nu_power_to_rhogrid(const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int pmgrid, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir, MPI_Comm MYMPI_COMM_WORLD)
{
  int x,y,z;
  _delta_pow d_pow = compute_neutrino_power_spectrum(Time, BoxSize, fft_of_rhogrid, pmgrid, slabstart_y, nslab_y, MYMPI_COMM_WORLD);
  /*Add P_nu to fft_of_rhgrid*/
  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < pmgrid; x++)
      for(z = 0; z < pmgrid / 2 + 1; z++)
        {
          double kx,ky,kz,k2,smth;
          int ip;
          kx = x > pmgrid/2 ? x-pmgrid : x;
          ky = y > pmgrid/2 ? y-pmgrid : y;
          kz = z > pmgrid/2 ? z-pmgrid : z;

          k2 = kx*kx + ky*ky + kz*kz;
          if(k2 <= 0)
              continue;
          /*Change the units of k to match those of logkk*/
          k2=log(sqrt(k2)*2*M_PI/BoxSize);
          /* Note get_neutrino_powerspec returns delta_nu / P_cdm^1/2, which is dimensionless.
           * We have delta_t = (M_cdm+M_nu)*delta_cdm (1-f_nu + f_nu (delta_nu / delta_cdm)^1/2)
           * which gives the right power spectrum, once we divide by
           * M_cdm +M_nu in powerspec*/
          smth=(1+get_dnudcdm_powerspec(&d_pow, k2));
          if(isnan(smth))
                terminate("delta_nu or delta_cdm is nan\n");
          ip = pmgrid * (pmgrid / 2 + 1) * (y - slabstart_y) + (pmgrid / 2 + 1) * x + z;
          fft_of_rhogrid[ip].re *= smth;
          fft_of_rhogrid[ip].im *= smth;
        }
  MPI_Barrier(MYMPI_COMM_WORLD);
  if(delta_tot_table.ThisTask==0)
	printf("Done adding neutrinos to grid on all processors\n");
  /*If this is being called to save all particle types, save a file with the neutrino power spectrum as well.*/
  if(snapnum >= 0 && delta_tot_table.ThisTask == 0){
            if(save_nu_power(&d_pow, Time, snapnum, OutputDir))
                terminate("Could not save neutrino power\n");
  }
  /*Free memory*/
  free_d_pow(&d_pow);
  myfree(d_pow.delta_cdm_curr);
  return;
}
