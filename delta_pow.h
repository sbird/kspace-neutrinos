#ifndef DELTA_POW_H
#define DELTA_POW_H

#include <gsl/gsl_interp.h>

struct _delta_pow{
    double *logkk;
    double *delta_nu_curr;
    gsl_interp *spline_nu;
    gsl_interp_accel * acc_nu;
    double *delta_cdm_curr;
    gsl_interp *spline_cdm;
    gsl_interp_accel * acc;
    int nbins;
    double norm;
};
typedef struct _delta_pow _delta_pow;

/* Initialise the structure for given power spectra. Note only pointers are stored; no copy of the actual array is made.
 * The thing interpolated is delta_nu_curr/delta_cdm_curr.*/
void init_delta_pow(_delta_pow *d_pow, double logkk[], double delta_nu_curr[], double delta_cdm_curr[], int nbins, double norm);

/*Get P_nu(k)/P_cdm(k) for arbitrary k*/
double get_dnudcdm_powerspec(_delta_pow *d_pow, double kk);

/*Save P_nu(k) to disc. Returns 0 on success.*/
int save_nu_power(_delta_pow *d_pow, const double Time, const int snapnum, const char * OutputDir);

/*Free memory for the GSL structure*/
void free_d_pow(_delta_pow * d_pow);
#endif
