#ifndef DELTA_POW_H
#define DELTA_POW_H
struct _delta_pow{
    double *logkk;
    double *delta_nu_curr;
    gsl_interp *spline_nu;
    gsl_interp_accel * acc_nu;
    double *delta_cdm_curr;
    gsl_interp *spline_cdm;
    gsl_interp_accel * acc;
    int nbins;
};
typedef struct _delta_pow _delta_pow;

void init_delta_pow(_delta_pow *d_pow, double logkk[], double delta_nu_curr[], double delta_cdm_curr[], int nbins);

double get_dnudcdm_powerspec(_delta_pow *d_pow, double kk);

void save_nu_power(_delta_pow *d_pow, const double Time, const int snapnum, const char * OutputDir);

void free_d_pow(_delta_pow * d_pow);
#endif
