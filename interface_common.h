/** \file
 * Global header, to be included in gadget's pm code. This defines the external interface to the kspace neutrino code.
 * These routines are generically useful for all N-body codes.*/
#ifndef KSPACE_NEUTRINOS_GLOBAL
#define KSPACE_NEUTRINOS_GLOBAL

#include "kspace_neutrino_const.h"
#include "delta_pow.h"
#include <mpi.h>

/**Global variables that need to be set from a parameter file*/
extern struct __kspace_params {
  /*File containing CAMB-formatted transfer functions*/
  char	KspaceTransferFunction[500];
  /*Scale factor of CAMB output*/
  double TimeTransfer;
  /*Units for the CAMB transfer functions in cm. By default 1 Mpc*/
  double InputSpectrum_UnitLength_in_cm;
  /*Neutrino masses in eV*/
  double MNu[NUSPECIES];
  /*Flag to enable hybrid neutrinos*/
  int hybrid_neutrinos_on;
  /*These variables are only used if hybrid_neutrinos_on = 1*/
  /*Critical velocity above which to treat neutrinos with particles.
  Note this is unperturbed velocity *TODAY* in Gadget units.
  To get velocity at redshift z, multiply by (1+z)*/
  double vcrit;
  /*Scale factor at which to turn on the particle neutrinos.*/
  double nu_crit_time;
} kspace_params;

/** Return the total matter density in all neutrino species.
 * This is not just OmegaNu(1)/a^3 because at early times neutrinos are relativistic.
 * The density in neutrino particles is included even if hybrid neutrinos are enabled.
 * Should be called from within the hubble function.
 * @param a scale factor. */
double OmegaNu(double a);

/** Return the matter density in neutrinos, 
 * excluding density in (hybrid) particles.
 * Identical to OmegaNu if hybrid particles are off.
 * Mostly useful for Gadget's check_omega.
 * @param a scale factor. */
double OmegaNu_nopart(double a);

void InitOmegaNu(const double TimeBegin,const double HubbleParam, const double tcmb0);

/** This function allocates memory for the neutrino tables, and loads the initial transfer
 * functions from CAMB transfer files.
 * One processor 0 it reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu_init and friends and has length NPowerTable is then broadcast to all processors.
 * Then, on all processors, it allocates memory for delta_tot_table.
 * This must be called *EARLY*, before OmegaNu is called for the first time (as that function
 * uses state set up here), just after the parameters are read.
 * Global state used: omeganu_table, delta_tot_table, transfer_init
 * @param nk_in number of bins desired in the neutrino power spectrum
 * @param ThisTask MPI rank
 * @param BoxSize size of box in internal units
 * @param UnitTime_in_s, UnitLength_in_cm - conversion factors from internal units to cgs.
 * @param Omega0 total matter density (including massive neutrinos and baryons but not including radiation)
 * @param HubbleParam Hubble parameter, eg, 0.7.
 * @param tcmb0 present-day CMB temperature in K
 * @param snapdir snapshot directory to try to read state and resume from
 * @param TimeMax Final scale factor desired, sets number of output redshift bins
 * @param MYMPI_COMM_WORLD MPI  communicator*/
void allocate_kspace_memory(const int nk_in, const int ThisTask,const double BoxSize, const double UnitTime_in_s, const double UnitLength_in_cm, const double Omega0, const char * snapdir, const double TimeMax, MPI_Comm MYMPI_COMM_WORLD);

/** This function calls the integrator to compute the neutrino power spectrum,
 * taking as input a pre-computed matter power spectrum, assumed to have the same units as stored in transfer_init.
 * neutrino power spectrum is stored in _delta_pow and returned.
 * Memory allocated here must be freed later.
 * Global state used: delta_tot_table, transfer_init, omeganu_table
 * @param Time scale factor, a.
 * @param nk_in Size of keff_in and P_cdm
 * @param keff_in k values for each power bin. Has units of UnitLength_in_cm passed to transfer_init
 * @param P_cdm Normalised matter power spectrum. Has units of UnitLength_in_cm passed to transfer_init
 * @param Nmodes number of modes in each bin. Used only to see if bin is nonempty.
 * @param MYMPI_COMM_WORLD MPI communicator to use
 * @returns _delta_pow, containing delta_nu/delta_cdm*/
_delta_pow compute_neutrino_power_from_cdm(const double Time, const double keff_in[], const double P_cdm[], const long int Nmodes[], const int nk_in, MPI_Comm MYMPI_COMM_WORLD);

/** Save the internal state of the integrator to disc.
 * @param savedir Output file is savedir/delta_tot_nu.txt.
 * Each row of the output file contains a scale factor and
 * the total matter power spectrum at that scale factor.*/
void save_nu_state(char * savedir);

/** Allocate memory and copy integrator internal state to it.
 * This may change and should be used with caution.
 * Included to allow saving the integrator state
 * using native save/load routines of the N-body code.
 * @param scalefact pointer containing address of array containing scale factors (malloced).
 * @param delta_tot pointer containing address of 2D array containing
 * delta_t values at each scale factor (malloced).
 * @param nk pointer to number of k values in delta_tot.
 * @param ia number of scale factors stored.*/
void get_nu_state(double ** scalefact, double ** delta_tot, size_t* nk, size_t* ia);

/** Get a pointer to the internal state of the integrator.
 * This may change and should be used with caution.
 * Included to allow loading the integrator state
 * from native save files of the N-body code.
 * @param scalefact pointer to array containing scale factors.
 * @param delta_tot pointer to 2D array containing
 * delta_t values at each scale factor.
 * @param nk number of k values in delta_tot.
 * @param ia number of scale factors stored.*/
void set_nu_state(double * scalefact, double * delta_tot, const size_t nk, const size_t ia, MPI_Comm MYMPI_COMM_WORLD);

/** Save a file containing the neutrino power spectrum.
 * Output to OutputDir/powerspec_nu_$(snapnum).txt
 * File format is:
 * Time
 * Nbins
 * k   P(k)   (repeated Nbins times)
 * .*/
int save_neutrino_power(const double Time, const int snapnum, const char * OutputDir);
/*KSPACE_NEUTRINOS_GLOBAL*/
#endif
