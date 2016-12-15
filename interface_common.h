/*Global header, to be included in gadget's pm code. This defines the external interface to the kspace neutrino code*/
#ifndef KSPACE_NEUTRINOS_GLOBAL
#define KSPACE_NEUTRINOS_GLOBAL

#include "kspace_neutrino_const.h"
#include "delta_pow.h"
#include <mpi.h>

/*Global variables that need to be set from a parameter file*/
extern struct __kspace_params {
  /*File containing CAMB-formatted transfer functions*/
  char	KspaceTransferFunction[500];
  /*Scale factor of CAMB output*/
  double TimeTransfer;
  /*OmegaBaryon used for CAMB*/
  double OmegaBaryonCAMB;
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

/* Return the total matter density in all neutrino species.
 * This is not just OmegaNu(1)/a^3 because at early times neutrinos are relativistic.
 * The density in neutrino particles is included even if hybrid neutrinos are enabled.
 * Should be called from within the hubble function.
 * Arguments: a - scale factor. */
double OmegaNu(double a);

/* Compute the matter density in neutrinos, 
 * excluding density in particles.
 * Mostly useful for Gadget's check_omega.*/
double OmegaNu_nopart(double a);

/** This function allocates memory for the neutrino tables, and loads the initial transfer
 * functions from CAMB transfer files.
 * One processor 0 it reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu_init and friends and has length NPowerTable is then broadcast to all processors.
 * Then, on all processors, it allocates memory for delta_tot_table.
 * This must be called *EARLY*, before OmegaNu is called for the first time (as that function
 * uses state set up here), just after the parameters are read.
 * Arguments:
 * nk_in - number of bins desired in the neutrino power spectrum
 * ThisTask - MPI rank
 * BoxSize - size of box in internal units
 * UnitTime_in_s, UnitLength_in_cm - conversion factors from internal units to cgs.
 * Omega0 - total matter density (including massive neutrinos and baryons but not including radiation)
 * tcmb0 - present-day CMB temperature
 * snapdir - snapshot directory to try to read state and resume from
 * TimeMax - Final redshift desired, sets number of output redshift bins
 * MYMPI_COMM_WORLD - MPI  communicator*/
void allocate_kspace_memory(const int nk_in, const int ThisTask,const double BoxSize, const double UnitTime_in_s, const double UnitLength_in_cm, const double Omega0, const double HubbleParam, const double tcmb0, const char * snapdir, const double TimeMax, MPI_Comm MYMPI_COMM_WORLD);

/*Compute the neutrino power spectrum using an externally computed matter power spectrum*/
_delta_pow compute_neutrino_power_from_cdm(const double Time, const double keff_in[], const double P_cdm[], const long long int Nmodes[], const int nk_in, MPI_Comm MYMPI_COMM_WORLD);

/*Save the internal state of the integrator to disc.
 * Arguments: savedir - Output file is savedir/delta_tot_nu.txt.
 * Each row of the output file contains a scale factor and
 * the total matter power spectrum at that scale factor.*/
void save_nu_state(char * savedir);
/*KSPACE_NEUTRINOS_GLOBAL*/
#endif
