#ifndef OMEGA_NU_SINGLE_H
#define OMEGA_NU_SINGLE_H
/*This file contains routines for computing the matter density in a single neutrino species*/
#include <gsl/gsl_interp.h>

#ifndef KSPACE_NEUTRINOS_TEST
#include "../gadgetconfig.h"
#endif

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
};
typedef struct _rho_nu_single _rho_nu_single;

/* Initialise the tables for the structure, by doing numerical integration*/
void rho_nu_init(_rho_nu_single * rho_nu_tab, double a0, const double mnu, double HubbleParam);

/* Compute the density, either by doing the */
double rho_nu(_rho_nu_single * rho_nu_tab, double a);

#ifdef HYBRID_NEUTRINOS
struct _hybrid_nu {
    /* This is the fraction of neutrino mass not followed by the analytic integrator.
    The analytic method is cutoff at q < qcrit (specified using vcrit, below) and use
    particles for the slower neutrinos.*/
    double nufrac_low[NUSPECIES];
    /* Time at which to turn on the particle neutrinos.
     * Ultimately we want something better than this.*/
    double nu_crit_time;
    /*Critical velocity as a fraction of lightspeed*/
    double vcrit;
};
typedef struct _hybrid_nu _hybrid_nu;

/*Set up parameters for the hybrid neutrinos
 * vcrit: Critical velocity above which to treat neutrinos with particles.
 *   Note this is unperturbed velocity *TODAY*
 *   To get velocity at redshift z, multiply by (1+z)
 * light: speed of light in internal units
 * nu_crit_time: critical time to make neutrino particles live
 */
void init_hybrid_nu(_hybrid_nu * hybnu, const double mnu[], const double vcrit, const double light, const double nu_crit_time);

/* Returns the fraction of neutrinos currently traced by particles.
 * When neutrinos are fully analytic at early times, returns 0.
 * Last argument: neutrino species to use.
 */
double particle_nu_fraction(_hybrid_nu * hybnu, const double a, int i);

/*Integrate the fermi-dirac kernel between 0 and qc to find the fraction of neutrinos that are particles*/
double nufrac_low(const double qc);
#endif

/*These are the structures you should call externally*/
struct _omega_nu {
    /*Pointers to the array of structures we use to store rho_nu*/
    _rho_nu_single * RhoNuTab[NUSPECIES];
    /* Which species have the same mass and can thus be counted together.*/
    int nu_degeneracies[NUSPECIES];
    /*Prefactor to turn density into matter density omega*/
    double rhocrit;
#ifdef HYBRID_NEUTRINOS
    _hybrid_nu hybnu;
#endif
};
typedef struct _omega_nu _omega_nu;

/*Initialise the above structure, allocating memory for the subclass rho_nu_single*/
void init_omega_nu(_omega_nu * omnu, const double MNu[], const double a0, const double HubbleParam);

/* Return the total matter density in neutrinos.*/
double get_omega_nu(_omega_nu *omnu, double a);

/*Return the photon matter density*/
double get_omegag(_omega_nu * omnu, double a);

/*Get the matter density in a single neutrino species*/
double omega_nu_single(_omega_nu * rho_nu_tab, double a, int i);

#endif
