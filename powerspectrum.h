#ifndef POWERSPEC_H
#define POWERSPEC_H

/*So that fftw_complex is defined*/
#include "kspace_neutrinos_2.h"

/* Compute the total powerspectrum from a Fourier-transformed density field in outfield, and store it in power.
 * Before use you may wish to normalise by dividing by count*count*/
void total_powerspectrum(const int dims, fftw_complex *outfield, const int nrbins, const int startslab, const int nslab, double *power, long long int *count, double *keffs, const double total_mass);

/*Rebin the power computed by total_powerspectrum by rebinning it to have even bins in logk. This reduces mode-by-mode dispersion*/
void rebin_power(double logk[],double delta_cdm_curr[],const int nbins, double Kbin[], double SumPower[], long long CountModes[], const int bins_ps, const int pmgrid, const double BoxSize);

#endif
