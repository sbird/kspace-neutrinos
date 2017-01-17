#ifndef POWERSPEC_H
#define POWERSPEC_H
#include <mpi.h>

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

/** Compute the total powerspectrum from a Fourier-transformed density field in outfield, and store it in power.
 * Before use you may wish to normalise by dividing by count^2 and 2*MPI/BoxSize. This assumes an FFTW2, slab decomposed, FFT.
 * Arguments:
 * @param dims size of one dimension of the density grid.
 * @param outfield Fourier transformed density grid.
 * @param nrbins Number of bins allocated for the output power spectrum.
 * @param startslab for slab parallelized FFT routines, this is the start index of the FFT on this rank.
 * @param nslab number of elements of the FFT on this rank.
 * @param power array for the output power spectrum. Returns power, in units of the box, UNNORMALISED.
 * Only bins with > 0 modes are included.
 * @param count Number of modes in each bin.
 * @param keffs effective k, averaged over all modes in a bin
 * @param MYMPI_COMM_WORLD MPI communicator.
 * @returns the number of bins in the output power spectrum, all of which have a non-zero number of modes.*/
int total_powerspectrum(const int dims, fftw_complex *outfield, const int nrbins, const int startslab, const int nslab, double *power, long long int *count, double *keffs, const MPI_Comm MYMPI_COMM_WORLD);

#endif
