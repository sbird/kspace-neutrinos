#ifndef DELTA_TOT_TABLE_H
#define DELTA_TOT_TABLE_H

#include "transfer_init.h"
#include "omega_nu_single.h"

/*Now we want to define a static object to store all previous delta_tot.
 * This object needs a constructor, a few private data members, and a way to be read and written from disk.
 * nk is fixed, delta_tot, scalefact and ia are updated in get_delta_nu_update*/
struct _delta_tot_table {
    /* Number of k values stored in each power spectrum*/
    int nk;
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
    /*Set to unity once the init routine has run*/
    int delta_tot_init_done;
    /* Pointer to nk arrays of length namax containing the total power spectrum.*/
    double **delta_tot;
    /* Array of length namax containing scale factors at which the power spectrum is stored*/
    double *scalefact;
    /*Pointer to array of length nk storing initial neutrino power spectrum*/
    double *delta_nu_init;
    /*Pointer to array of length nk storing the last neutrino power spectrum we saw, for a first estimate
    * of the new delta_tot */
    double *delta_nu_last;
    /*Light speed in internal units. C is defined in allvars.h to be lightspeed in cm/s*/
    double light;
    /*The time at which we first start our integrator:
     * NOTE! This is not All.TimeBegin, but the time of the transfer function file,
     * so that we can support restarting from snapshots.*/
    double TimeTransfer;
};
typedef struct _delta_tot_table _delta_tot_table;

/*This function allocates memory for delta_tot_table*/
void allocate_delta_tot_table(_delta_tot_table *d_tot, int nk_in, const double TimeTransfer);

/*Initialise the data in delta_tot_init from the transfer table data in transfer_init.
 This is separate from allocate_delta_tot_table because we need some information not available when the memory needs to be allocated*/
void delta_tot_init(_delta_tot_table *d_tot, int nk_in, double wavenum[], double delta_cdm_curr[], _transfer_init_table *t_init, const double Omega0, _omega_nu * omnu, const double UnitTime_in_s, const double UnitLength_in_cm);

/*Function called by add_nu_power_to_rhogrid*/
void get_delta_nu_update(_delta_tot_table *d_tot, double a, int nk_in, double wavenum[], double P_cdm_curr[], double delta_nu_curr[], _omega_nu * omnu);

/*This function does the work and updates delta_nu_curr*/
void get_delta_nu(_delta_tot_table *d_tot, double a, int Na, double wavenum[], double delta_nu_curr[],double mnu);

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(_delta_tot_table *d_tot, double a, int Na, double wavenum[],  double delta_nu_curr[], _omega_nu * omnu);

/*Save a single line in the delta_tot table to a file*/
void save_delta_tot(_delta_tot_table *d_tot, int iia, char * savedir);

/*Save a complete delta_tot table to disc*/
void save_all_nu_state(_delta_tot_table *d_tot, char * savedir);

/* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
 * Must be called before delta_tot_init, or resuming wont work*/
void read_all_nu_state(_delta_tot_table *d_tot, char * savedir, double Time);

#endif
