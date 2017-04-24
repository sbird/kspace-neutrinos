/** \file
 * Global header, to be included in gadget's pm code. This defines the external interface to the kspace neutrino code.
 * These routines are specific to FFTW2 using codes.*/
#ifndef KSPACE_NEUTRINOS_GADGET
#define KSPACE_NEUTRINOS_GADGET

#ifdef KSPACE_NEUTRINOS
#error "KSPACE_NEUTRINOS_2 is incompatible with KSPACE_NEUTRINOS"
#endif

#include "interface_common.h"

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

/** Main function to include neutrino power in the code. Call it from pm_periodic.c
 * This function computes the neutrino power spectrum and adds it to the
 * density grid. It calls the internal power spectrum routine and the neutrino integrator.
 * It then adds the neutrino power to fft_of_rhogrid, which is the fourier transformed density grid from the PM code.
 * @param Time scale factor, a.
 * @param BoxSize size of the box in internal units.
 * @param fft_of_rhogrid Fourier transformed density grid.
 * @param pmgrid size of one dimension of the density grid.
 * @param slabstart_y for slab parallelized FFT routines, this is the start index of the FFT on this rank.
 * @param nslab_y number of elements of the FFT on this rank.
 * @param snapnum number of snapshot to save neutrino power spectrum as powerspec_nu_$(snapnum).txt
 * @param OutputDir output directory for neutrino power spectrum. Neutrino power will be saved if this is non-null.
 * @param MYMPI_COMM_WORLD MPI communicator to use
 */
void add_nu_power_to_rhogrid(const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int pmgrid, int slabstart_y, int nslab_y, MPI_Comm MYMPI_COMM_WORLD);

/** Function which sets up the parameter reader to read kspace neutrino parameters from the parameter file. 
 * It will store them in a static variable, kspace_params, in the translation unit where the function is defined
 * (which is the same as the above functions).
 * This is an example for Gadget-3, using that codes parameter reading scheme.
 * It may need changing for your code!*/
int set_kspace_vars(char tag[][50], void *addr[], int id [], int nt);

/*Compute the total matter power spectrum: will be saved somewhere where save_total_power can read it out.*/
void compute_total_power_spectrum(const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int pmgrid, int slabstart_y, int nslab_y, MPI_Comm MYMPI_COMM_WORLD);

/** Save a file containing the total power spectrum.
 * Output to OutputDir/powerspec_nu_$(snapnum).txt
 * File format is:
 * Time
 * Nbins
 * k   P(k)   (repeated Nbins times)
 * .*/
int save_total_power(const double Time, const int snapnum, const char * OutputDir);
#endif
