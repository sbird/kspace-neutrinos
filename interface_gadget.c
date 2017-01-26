/* This file contains functions which need to be called from the PM code in Gadget-3.
 * Interface functions which assume FFTW2 are here, others are in interface_common.c
 * All MPI communication is also done here.
 * add_nu_power_to_rhogrid is the main public function. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "interface_gadget.h"
#include "powerspectrum.h"
#include "kspace_neutrino_const.h"
#include "gadget_defines.h"
#include "omega_nu_single.h"
#include "transfer_init.h"
#include "delta_tot_table.h"
#include "delta_pow.h"

/*Global neutrino module parameters*/
extern _transfer_init_table transfer_init;

extern _delta_tot_table delta_tot_table;

extern _omega_nu omeganu_table;

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
  /*Free memory*/
  free_d_pow(&d_pow);
  myfree(d_pow.delta_cdm_curr);
  return;
}
