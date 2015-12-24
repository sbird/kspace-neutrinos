#ifndef OMEGA_NU_SINGLE_H
#define OMEGA_NU_SINGLE_H
/*This file contains routines for computing the matter density in a single neutrino species*/
#include <gsl/gsl_interp.h>

/* for three massive neutrino species:
* Could be made configurable at some point
* Neutrino masses are in eV*/
#define NUSPECIES 3
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
void rho_nu_init(_rho_nu_single * rho_nu_tab, double a0, const double mnu, double HubbleParam);

/* Compute the density, either by doing the */
double rho_nu(_rho_nu_single * rho_nu_tab, double a);

double omega_nu_single(_rho_nu_single * rho_nu_tab, double a);

/*These are the structures you should call externally*/
struct _omega_nu {
    /*Pointers to the array of structures we use to store rho_nu*/
    _rho_nu_single * RhoNuTab[NUSPECIES];
    /* Which species have the same mass and can thus be counted together.
     */
    int nu_degeneracies[NUSPECIES];
    double MNu[NUSPECIES];
    /*Matter fraction*/
    double Omega0;
};
typedef struct _omega_nu _omega_nu;

/*Initialise the above structure, allocating memory for the subclass rho_nu_single*/
void init_omega_nu(_omega_nu * omnu, const double MNu[], const double Omega0);

/* Return the total matter density in neutrinos.*/
double OmegaNu(_omega_nu *omnu, double a);

#endif
