#ifndef OMEGA_NU_SINGLE_H
#define OMEGA_NU_SINGLE_H
/*This file contains routines for computing the matter density in a single neutrino species*/
#include <gsl/gsl_interp.h>
#include "kspace_neutrino_const.h"

/* Ratio between the massless neutrino temperature and the CMB temperature.
 * Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 * We use the CLASS default value, chosen so that omega_nu = m_nu / 93.14 h^2
 * At time of writing this is T_nu / T_gamma = 0.71611.
 * See https://github.com/lesgourg/class_public/blob/master/explanatory.ini
 */
#define TNUCMB     (pow(4/11.,1/3.)*1.00328)              /* Neutrino + antineutrino background temperature in Kelvin */

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
void rho_nu_init(_rho_nu_single * rho_nu_tab, double a0, const double mnu, const double HubbleParam, const double kBtnu);

/* Compute the density, either by looking up in a table, or a simple calculation in the limits, or by direct integration.*/
double rho_nu(_rho_nu_single * rho_nu_tab, const double a, const double kT);

/*The following functions and structures are used for hybrid neutrinos only*/
struct _hybrid_nu {
    /*True if hybrid neutrinos are enabled*/
    int enabled;
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
void init_hybrid_nu(_hybrid_nu * const hybnu, const double mnu[], const double vcrit, const double light, const double nu_crit_time, const double kBtnu);

/* Returns the fraction of neutrinos currently traced by particles.
 * When neutrinos are fully analytic at early times, returns 0.
 * Last argument: neutrino species to use.
 */
double particle_nu_fraction(const _hybrid_nu * const hybnu, const double a, int i);

/*Integrate the fermi-dirac kernel between 0 and qc to find the fraction of neutrinos that are particles*/
double nufrac_low(const double qc);

/*End hybrid neutrino-only structures.*/

/*These are the structures you should call externally*/
struct _omega_nu {
    /*Pointers to the array of structures we use to store rho_nu*/
    _rho_nu_single * RhoNuTab[NUSPECIES];
    /* Which species have the same mass and can thus be counted together.*/
    int nu_degeneracies[NUSPECIES];
    /* Prefactor to turn density into matter density omega*/
    double rhocrit;
    /*neutrino temperature times Boltzmann constant*/
    double kBtnu;
    /*CMB temperature*/
    double tcmb0;
    /* Pointer to structure for hybrid neutrinos. */
    _hybrid_nu hybnu;
};
typedef struct _omega_nu _omega_nu;

/*Initialise the above structure, allocating memory for the subclass rho_nu_single*/
void init_omega_nu(_omega_nu * const omnu, const double MNu[], const double a0, const double HubbleParam, const double tcmb0);

/* Return the total matter density in neutrinos.*/
double get_omega_nu(const _omega_nu * const omnu, const double a);

/* Return the total matter density in neutrinos, excluding active particles.*/
double get_omega_nu_nopart(const _omega_nu * const omnu, const double a);

/*Return the photon matter density*/
double get_omegag(const _omega_nu * const omnu, const double a);

/*Get the matter density in a single neutrino species*/
double omega_nu_single(const _omega_nu * const rho_nu_tab, const double a, const int i);

#endif
