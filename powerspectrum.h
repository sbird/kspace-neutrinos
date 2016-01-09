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

/* Compute the total powerspectrum from a Fourier-transformed density field in outfield, and store it in power.
 * Before use you may wish to normalise by dividing by count*count*/
void total_powerspectrum(const int dims, fftw_complex *outfield, const int nrbins, const int startslab, const int nslab, double *power, long long int *count, double *keffs, const double total_mass, const MPI_Comm MYMPI_COMM_WORLD);

#endif
