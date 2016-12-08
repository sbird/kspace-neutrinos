/*Global header, to be included in gadget's pm code. This defines the external interface to the kspace neutrino code*/
#ifndef KSPACE_NEUTRINOS_GLOBAL
#define KSPACE_NEUTRINOS_GLOBAL

#ifdef KSPACE_NEUTRINOS
#error "KSPACE_NEUTRINOS_2 is incompatible with KSPACE_NEUTRINOS"
#endif

/*We only need this for fftw_complex*/
#ifdef NOTYPEPREFIX_FFTW
#include        <fftw.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <dfftw.h>	/* double precision FFTW */
#else
#include     <sfftw.h>
#endif
#endif

#include <mpi.h>

/* Return the total matter density in all neutrino species.
 * This is not just OmegaNu(1)/a^3 because at early times neutrinos are relativistic.
 * The density in neutrino particles is included even if hybrid neutrinos are enabled.
 * Should be called from within the hubble function.
 * Arguments: a - scale factor. */
double OmegaNu(double a);

/** This function allocates memory for the neutrino tables, and loads the initial transfer
 * functions from CAMB transfer files.
 * One processor 0 it reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu_init and friends and has length NPowerTable is then broadcast to all processors.
 * Then, on all processors, it allocates memory for delta_tot_table.
 * This must be called *EARLY*, before OmegaNu is called for the first time (as that function
 * uses state set up here), just after the parameters are read.
 * Arguments:
 * nk_in - number of bins desired in the neutrino power spectrum
 * ThisTask - MPI rank
 * BoxSize - size of box in internal units
 * UnitTime_in_s, UnitLength_in_cm - conversion factors from internal units to cgs.
 * Omega0 - total matter density (including massive neutrinos and baryons but not including radiation)
 * tcmb0 - present-day CMB temperature
 * snapdir - snapshot directory to try to read state and resume from
 * TimeMax - Final redshift desired, sets number of output redshift bins
 * MYMPI_COMM_WORLD - MPI  communicator*/
void allocate_kspace_memory(const int nk_in, const int ThisTask,const double BoxSize, const double UnitTime_in_s, const double UnitLength_in_cm, const double Omega0, const double HubbleParam, const double tcmb0, const char * snapdir, const double TimeMax, MPI_Comm MYMPI_COMM_WORLD);

/* Main function to include neutrino power in the code. Call it from pm_periodic.c
 * This function computes the neutrino power spectrum and adds it to the
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
 * OutputDir - output directory for neutrino power spectrum. Neutrino power will be saved if this is non-null.
 * MYMPI_COMM_WORLD - MPI communicator to use
 */
void add_nu_power_to_rhogrid(const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int pmgrid, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir, MPI_Comm MYMPI_COMM_WORLD);

/* Function which sets up the parameter reader to read kspace neutrino parameters from the parameter file. 
 * It will store them in a static variable, kspace_params, in the translation unit where the function is defined
 * (which is the same as the above functions). */
int set_kspace_vars(char tag[][50], void *addr[], int id [], int nt);

/*Save the internal state of the integrator to disc.
 * Arguments: savedir - Output file is savedir/delta_tot_nu.txt.
 * Each row of the output file contains a scale factor and
 * the total matter power spectrum at that scale factor.*/
void save_nu_state(char * savedir);
/*KSPACE_NEUTRINOS_GLOBAL*/
#endif
