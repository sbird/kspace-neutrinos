#ifndef DELTA_TOT_TABLE_H
#define DELTA_TOT_TABLE_H

#include "transfer_init.h"
#include "omega_nu_single.h"

/*Now we want to define a static object to store all previous delta_tot.
 * This object needs a constructor, a few private data members, and a way to be read and written from disk.
 * nk is fixed, delta_tot, scalefact and ia are updated in get_delta_nu_update*/
struct _delta_tot_table {
    /* Number of actually non-zero k values stored in each power spectrum*/
    int nk;
    /* Size of arrays allocated to store power spectra*/
    int nk_allocated;
    /*Maximum number of redshifts to store. Redshifts are stored every delta a = 0.01 */
    int namax;
    /* Number of already "recorded" time steps, i.e. scalefact[0...ia-1] is recorded.
    * Current time corresponds to index ia (but is only recorded if sufficiently far from previous time).
    * Caution: ia here is different from Na in get_delta_nu (Na = ia+1).*/
    int ia;
    /* MPI rank of this processor*/
    int ThisTask;
    /*Prefactor for use in get_delta_nu. Should be 3/2 Omega_m H^2 /c */
    double delta_nu_prefac;
    /*Set to unity once the init routine has run.*/
    int delta_tot_init_done;
    /*If greater than 0, intermediate files will be saved and status output will be displayed*/
    int debug;
    /* Pointer to nk arrays of length namax containing the total power spectrum.*/
    double **delta_tot;
    /* Array of length namax containing scale factors at which the power spectrum is stored*/
    double * scalefact;
    /*Pointer to array of length nk storing initial neutrino power spectrum*/
    double * delta_nu_init;
    /*Pointer to array of length nk storing the last neutrino power spectrum we saw, for a first estimate
    * of the new delta_tot */
    double * delta_nu_last;
    /*Pointer to a structure for computing omega_nu*/
    const _omega_nu * omnu;
    /*Matter density excluding neutrinos*/
    double Omeganonu;
    /*Light speed in internal units. C is defined in allvars.h to be lightspeed in cm/s*/
    double light;
    /*The time at which we first start our integrator:
     * NOTE! This is not All.TimeBegin, but the time of the transfer function file,
     * so that we can support restarting from snapshots.*/
    double TimeTransfer;
};
typedef struct _delta_tot_table _delta_tot_table;

/*This function allocates memory for delta_tot_table*/
void allocate_delta_tot_table(_delta_tot_table *d_tot, const int nk_in, const double TimeTransfer, const double TimeMax, const double Omega0, const _omega_nu * const omnu, const double UnitTime_in_s, const double UnitLength_in_cm, int debug);

/*Frees the memory allocated above*/
void free_delta_tot_table(_delta_tot_table *d_tot);

/*Initialise the data in delta_tot_init from the transfer table data in transfer_init.
 This is separate from allocate_delta_tot_table because we need some information not available when the memory needs to be allocated*/
void delta_tot_init(_delta_tot_table * const d_tot, const int nk_in, const double wavenum[], const double delta_cdm_curr[], const _transfer_init_table * const t_init, const double Time);

/*Update the last value of delta_tot in the table with a new value computed
 from the given delta_cdm_curr and delta_nu_curr.
 If overwrite is true, overwrite the existing final entry.*/
void update_delta_tot(_delta_tot_table * const d_tot, const double a, const double delta_cdm_curr[], const double delta_nu_curr[], const int overwrite);

/*Function called by add_nu_power_to_rhogrid*/
void get_delta_nu_update(_delta_tot_table * const d_tot, const double a, const int nk_in, const double keff[], const double P_cdm_curr[], double delta_nu_curr[]);

/*This function does the work and updates delta_nu_curr*/
void get_delta_nu(const _delta_tot_table * const d_tot, const double a, const double wavenum[], double delta_nu_curr[], const double mnu);

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(const _delta_tot_table * const d_tot, const double a, const double wavenum[],  double delta_nu_curr[]);

/*Save a single line in the delta_tot table to a file*/
void save_delta_tot(const _delta_tot_table *const d_tot, const int iia, char * savedir);

/*Save a complete delta_tot table to disc*/
void save_all_nu_state(const _delta_tot_table * const d_tot, char * savedir);

/* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
 * Must be called before delta_tot_init, or resuming wont work*/
void read_all_nu_state(_delta_tot_table * const d_tot, const char * savedir);

/*Fit to the special function J(x) that is accurate to better than 3% relative and 0.07% absolute*/
double specialJ(const double x, const double vcmnubylight);

/* Free-streaming length for a non-relativistic particle of momentum q = T0, from scale factor ai to af.
 * Result is in Unit_Length.*/
double fslength(const double logai, const double logaf, const double mnu, const double light);

#endif
