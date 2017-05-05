#ifndef DELTA_POW_H
#define DELTA_POW_H
/**\file
 * Functions for manipulating _delta_pow, which stores and interpolates a power spectrum.
 */
#include <gsl/gsl_interp.h>

/** Opaque structure storing the interpolation functions and pointers to the memory
 * holding the power spectra*/
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

/** Initialise the structure for given power spectra. This initialises various interpolation routines 
 * from precomputed delta_cdm and delta_nu.
 * Note only pointers are stored; no copy of the actual array is made, so you must not free memory while d_pow is active.
 * The thing interpolated is delta_nu_curr/delta_cdm_curr.
 * @param d_pow (opaque) structure being initialised.
 * @param logkk Array containing log-k
 * @param delta_nu_curr array containing neutrino power spectrum at bins specified by logkk
 * @param delta_cdm_curr array containing CDM power spectrum at bins specified by logkk
 * @param nbins number of bins in earlier arrays.
 * @param norm constant which multiplies the value returned by get_dnudcdm_powerspec. 
 * Default call in interface_common.c sets it to Omega_nu/(Omega_0-Omega_nu)*/
void init_delta_pow(_delta_pow *d_pow, double logkk[], double delta_nu_curr[], double delta_cdm_curr[], int nbins, double norm);

/**Get P_nu(k)/P_cdm(k) for arbitrary k. This will become: 
 * \f$\delta_\nu = \delta_{CDM}(\vec{k}) \left(\delta_\nu(k) / \delta_{CDM} (k)\right)\f$,
 * thus we get the right powerspectrum.
 * @param d_pow (opaque) structure containing stored power spectrum and GSL interpolators.
 * @param kk log(k) value to get delta_nu at
 * @returns Omega_nu/(Omega_0-Omega_nu) * delta_nu / delta_CDM
 * */
double get_dnudcdm_powerspec(_delta_pow *d_pow, double kk);

/** Free memory for the GSL structure*/
void free_d_pow(_delta_pow * d_pow);

double get_delta_tot(const double delta_nu_curr, const double delta_cdm_curr, const double OmegaNua3, const double Omeganonu, const double Omeganu1);
#endif
