#ifndef OMEGA_NU_SINGLE_H
#define OMEGA_NU_SINGLE_H
/*This file contains routines for computing the matter density in a single neutrino species*/
#include <gsl/gsl_interp.h>

#define NRHOTAB 500
/* Tables for rho_nu: stores precomputed values between
 * simulation start and a M_nu = 20 kT_nu for a single neutrino species.*/
struct _rho_nu_single {
    double loga[NRHOTAB];
    double rhonu[NRHOTAB];
    gsl_interp * interp;
    gsl_interp_accel * acc;
    /*Neutrino mass for this structure*/
    double mnu;
    /*Prefactor to turn density into matter density omega*/
    double omega_prefac;
};
typedef struct _rho_nu_single _rho_nu_single;

/* Initialise the tables for the structure, by doing numerical integration*/
void rho_nu_init(_rho_nu_single * rho_nu_tab, double a0, double mnu);

/* Compute the density, either by doing the */
double rho_nu(_rho_nu_single * rho_nu_tab, double a);

double omega_nu_single(_rho_nu_single * rho_nu_tab, double a);

#endif