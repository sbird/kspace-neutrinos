#include <stdio.h>
#include <gsl/gsl_interp.h>
#include <math.h>
#include "delta_pow.h"

void init_delta_pow(_delta_pow *d_pow, double logkk[], double delta_ratio[], int nbins, double norm)
{
  /*Set up array pointers*/
  d_pow->logkk = logkk;
  d_pow->delta_ratio = delta_ratio;
  d_pow->nbins = nbins;
  d_pow->norm = norm;
  /*Set up interpolation structures*/
  d_pow->acc = gsl_interp_accel_alloc();
  d_pow->spline=gsl_interp_alloc(gsl_interp_linear,nbins);
  gsl_interp_init(d_pow->spline,d_pow->logkk,d_pow->delta_ratio,nbins);
}

double get_dnudcdm_powerspec(_delta_pow *d_pow, double kk)
{
        double dnudcdm;
        /* Floating point roundoff and the binning means there may be a mode just beyond the box size.
         * For now we assume P(k) is constant on these large scales.
         * At some point in the future, linear extrapolation could be used.*/
        if(kk < d_pow->logkk[0] && kk > d_pow->logkk[0]-log(2) ){
            kk = d_pow->logkk[0];
        }
        if(kk < d_pow->logkk[0]-log(2)){
            fprintf(stderr,"trying to extract a k= %g < min stored = %g \n",kk,d_pow->logkk[0]);
            kk = d_pow->logkk[0];
        }
        /*This is just to guard against floating point roundoff*/
        if( kk > d_pow->logkk[d_pow->nbins-1])
                kk=d_pow->logkk[d_pow->nbins-1];
        dnudcdm=gsl_interp_eval(d_pow->spline,d_pow->logkk, d_pow->delta_ratio,kk,d_pow->acc);
        return d_pow->norm * dnudcdm;
}

void free_d_pow(_delta_pow * d_pow)
{
  gsl_interp_free(d_pow->spline);
  gsl_interp_accel_free(d_pow->acc);
}
