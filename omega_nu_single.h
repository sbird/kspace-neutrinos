#ifndef OMEGA_NU_SINGLE_H
#define OMEGA_NU_SINGLE_H
/*This file contains routines for computing the matter density in a single neutrino species*/
#include <gsl/gsl_interp.h>

/* for three massive neutrino species:
* Could be made configurable at some point
* Neutrino masses are in eV*/
#define NUSPECIES 3
/* Tables for rho_nu: stores precomputed values between
 * simulation start and a M_nu = 20 kT_nu for a single neutrino species.*/
struct _rho_nu_single {
    double * loga;
    double * rhonu;
    gsl_interp * interp;
    gsl_interp_accel * acc;
    /*Neutrino mass for this structure*/
    double mnu;
#ifdef HYBRID_NEUTRINOS
    /* If this is zero, then we proceed using the analytic method for all neutrinos.
    If this is nonzero, then we assume this fraction of neutrino mass is not followed by the analytic integrator.
    Instead cut off the analytic method at q < qcrit (specified using vcrit, below) and use
    particles for the slower neutrinos.*/
    double nufrac_low;
#endif
};
typedef struct _rho_nu_single _rho_nu_single;

/* Initialise the tables for the structure, by doing numerical integration*/
void rho_nu_init(_rho_nu_single * rho_nu_tab, double a0, const double mnu, double HubbleParam);

/* Compute the density, either by doing the */
double rho_nu(_rho_nu_single * rho_nu_tab, double a);

/*These are the structures you should call externally*/
struct _omega_nu {
    /*Pointers to the array of structures we use to store rho_nu*/
    _rho_nu_single * RhoNuTab[NUSPECIES];
    /* Which species have the same mass and can thus be counted together.*/
    int nu_degeneracies[NUSPECIES];
    /*Matter fraction*/
    double Omega0;
    /*Prefactor to turn density into matter density omega*/
    double rhocrit;
#ifdef HYBRID_NEUTRINOS
    /*Are the neutrinos still analytic?*/
    int neutrinos_not_analytic;
    /* Critical velocity above which to treat neutrinos with particles.
    Note this is unperturbed velocity *TODAY*
    To get velocity at redshift z, multiply by (1+z)*/
    double vcrit;
    /* Time at which to turn on the particle neutrinos.
     * Ultimately we want something better than this.*/
    double nu_crit_time;
#endif
};
typedef struct _omega_nu _omega_nu;

/*Initialise the above structure, allocating memory for the subclass rho_nu_single*/
void init_omega_nu(_omega_nu * omnu, const double MNu[], const double Omega0, const double a0, const double HubbleParam);

/* Return the total matter density in neutrinos.*/
double get_omega_nu(_omega_nu *omnu, double a);

/*Return the photon matter density*/
double get_omegag(_omega_nu * omnu, double a);

/*Get the matter density in a single neutrino species*/
double omega_nu_single(_omega_nu * rho_nu_tab, double a, int i);


#ifdef HYBRID_NEUTRINOS
/*Check whether the neutrinos are analytic or not*/
int slow_neutrinos_analytic(_omega_nu * omnu, const double a, const double light);
#endif

#endif
