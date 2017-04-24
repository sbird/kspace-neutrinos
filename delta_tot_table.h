#ifndef DELTA_TOT_TABLE_H
#define DELTA_TOT_TABLE_H

/**\file
 * Contains calculations for the Fourier-space semi-linear neutrino method
 * described in Ali-Haimoud and Bird 2012.
 * delta_tot_table stores the state of the integrator, which includes the matter power spectrum over all past time.
 * This file contains routines for manipulating this structure; updating it by computing a new neutrino power spectrum,
 * from the non-linear CDM power, and saving and loading the structure to and from disc.
 */
#include "transfer_init.h"
#include "omega_nu_single.h"

/** Now we want to define a static object to store all previous delta_tot.
 * This object needs a constructor, a few private data members, and a way to be read and written from disk.
 * nk is fixed, delta_tot, scalefact and ia are updated in get_delta_nu_update*/
struct _delta_tot_table {
    /** Number of actually non-zero k values stored in each power spectrum*/
    int nk;
    /** Size of arrays allocated to store power spectra*/
    int nk_allocated;
    /** Maximum number of redshifts to store. Redshifts are stored every delta a = 0.01 */
    int namax;
    /** Number of already "recorded" time steps, i.e. scalefact[0...ia-1] is recorded.
    * Current time corresponds to index ia (but is only recorded if sufficiently far from previous time).
    * Caution: ia here is different from Na in get_delta_nu (Na = ia+1).*/
    int ia;
    /** MPI rank of this processor*/
    int ThisTask;
    /** Prefactor for use in get_delta_nu. Should be 3/2 Omega_m H^2 /c */
    double delta_nu_prefac;
    /** Set to unity once the init routine has run.*/
    int delta_tot_init_done;
    /** If greater than 0, intermediate files will be saved and status output will be displayed*/
    int debug;
    /** Pointer to nk arrays of length namax containing the total power spectrum.*/
    double **delta_tot;
    /** Array of length namax containing scale factors at which the power spectrum is stored*/
    double * scalefact;
    /** Pointer to array of length nk storing initial neutrino power spectrum*/
    double * delta_nu_init;
    /** Pointer to array of length nk storing the last neutrino power spectrum we saw, for a first estimate
    * of the new delta_tot */
    double * delta_nu_last;
    /**Pointer to array storing the effective wavenumbers for the above power spectra*/
    double * wavenum;
    /** Pointer to a structure for computing omega_nu*/
    const _omega_nu * omnu;
    /** Matter density excluding neutrinos*/
    double Omeganonu;
    /** Light speed in internal units. C is defined in allvars.h to be lightspeed in cm/s*/
    double light;
    /** The time at which we first start our integrator:
     * NOTE! This is not All.TimeBegin, but the time of the transfer function file,
     * so that we can support restarting from snapshots.*/
    double TimeTransfer;
};
typedef struct _delta_tot_table _delta_tot_table;

/** Allocates memory for delta_tot_table.
 * @param d_tot structure to initialise
 * @param nk_in Number of bins stored in each power spectrum.
 * @param TimeTransfer Scale factor of the transfer functions.
 * @param TimeMax Final scale factor up to which we will need memory.
 * @param Omega0 Matter density at z=0.
 * @param omnu Pointer to structure containing pre-computed tables for evaluating neutrino matter densities.
 * @param UnitTime_in_s Time unit of the simulation in s.
 * @param UnitLength_in_cm Length unit of the simulation in cm
 * @param debug If this is > zero, there will be extra output.*/
void allocate_delta_tot_table(_delta_tot_table *d_tot, const int nk_in, const double TimeTransfer, const double TimeMax, const double Omega0, const _omega_nu * const omnu, const double UnitTime_in_s, const double UnitLength_in_cm, int debug);

/** Frees the memory allocated above*/
void free_delta_tot_table(_delta_tot_table *d_tot);

/** Constructor. transfer_init_tabulate must be called before this function.
 * Initialises delta_tot (including from a file) and delta_nu_init from the transfer functions.
 * read_all_nu_state must be called before this if you want reloading from a snapshot to work
 * Note delta_cdm_curr includes baryons, and is only used if not resuming.
 * @param d_tot Structure allocated with enough memory to hold the power spectra.
 * @param nk_in number of k bins for power spectra.
 * @param wavenum Values of k (not log k!) for each power spectrum bin.
 * @param delta_cdm_curr Current (at this timestep) value of the CDM (+baryon) power spectrum. Source term for the neutrino power
 * @param t_init Table for the initial transfer function.
 * @param Time Current scale factor. */
void delta_tot_init(_delta_tot_table * const d_tot, const int nk_in, const double wavenum[], const double delta_cdm_curr[], const _transfer_init_table * const t_init, const double Time);

/** Update the last value of delta_tot in the table with a new value computed
 from the given delta_cdm_curr and delta_nu_curr.
 If overwrite is true, overwrite the existing final entry.*/
void update_delta_tot(_delta_tot_table * const d_tot, const double a, const double delta_cdm_curr[], const double delta_nu_curr[], const int overwrite);

/** Callable function to calculate the power spectra.
 * Calls the rest of the code internally.
 * In reality we will be given P_cdm(current) but not delta_tot.
 * Here is the full function that deals with this
 * @param d_tot contains the state of the integrator; samples of the total power spectrum at earlier times.
 * @param a is the current scale factor
 * @param nk_in is the number of k bins in delta_cdm_curr and keff.
 * @param keff is an array of length nk containing (natural) log k
 * @param P_cdm_curr array of length nk containing the square root of the current cdm power spectrum
 * @param delta_nu_curr is an array of length nk which stores the square root of the current neutrino power spectrum. Main output of the function.
 * @param transfer_init is a pointer to the structure containing transfer tables.
******************************************************************************************************/
void get_delta_nu_update(_delta_tot_table * const d_tot, const double a, const int nk_in, const double keff[], const double P_cdm_curr[], double delta_nu_curr[], _transfer_init_table * transfer_init);

/** Main function: given tables of wavenumbers, total delta at Na earlier times (< = a),
 * and initial conditions for neutrinos, computes the current delta_nu.
 * @param d_tot Initialised structure for storing total matter density.
 * @param a Current scale factor.
 * @param wavenum Values of k (not log k!) for each power spectrum bin.
 * @param delta_nu_curr Pointer to array to store square root of neutrino power spectrum. Main output.
 * @param mnu Neutrino mass in eV.*/
void get_delta_nu(const _delta_tot_table * const d_tot, const double a, const double wavenum[], double delta_nu_curr[], const double mnu);

/** Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(const _delta_tot_table * const d_tot, const double a, const double wavenum[],  double delta_nu_curr[]);

/** Save a single line in the delta_tot table to a file*/
void save_delta_tot(const _delta_tot_table *const d_tot, const int iia, char * savedir);

/** Save a complete delta_tot table to disc*/
void save_all_nu_state(const _delta_tot_table * const d_tot, char * savedir);

/** Save P_nu(k) to disc.
 * @param d_tot will save delta_nu_last from d_tot.
 * @param Time output scale factor.
 * @param snapnum Number of snapshot. Will save to powerspec_nu_$(snapnum).txt.
 * @param OutputDir Output directory to save data in.
 * @returns 0 on success*/
int save_nu_power(const _delta_tot_table * const d_tot, const double Time, const int snapnum, const char * OutputDir);

/** Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
 * Must be called before delta_tot_init, or resuming wont work*/
void read_all_nu_state(_delta_tot_table * const d_tot, char * savedir);

/** Fit to the special function J(x) that is accurate to better than 3% relative and 0.07% absolute*/
double specialJ(const double x, const double vcmnubylight);

/** Free-streaming length (times Mnu/k_BT_nu, which is dimensionless) for a non-relativistic
particle of momentum q = T0, from scale factor ai to af.
Arguments:
@param logai log of initial scale factor
@param logaf log of final scale factor
@param mnu Neutrino mass in eV
@param light speed of light in internal length units.
@returns free-streaming length in Unit_Length/Unit_Time (same units as light parameter).
*/
double fslength(const double logai, const double logaf, const double light);

#endif
