#ifndef TRANSFER_INIT_H
#define TRANSFER_INIT_H

/** \file 
 * Transfer Functions: this file contains routines to read CAMB transfer functions into memory*/

/** Structure to store the initial transfer functions from CAMB.
 * We store transfer functions because we want to use the
 * CDM + Baryon total matter power spectrum from the
 * first timestep of Gadget, so that possible Rayleigh scattering
 * in the initial conditions is included in the neutrino and radiation components. */
struct _transfer_init_table {
    int NPowerTable;
    double *logk;
    /*This is T_nu / (T_not-nu), where T_not-nu is a weighted average of T_cdm and T_baryon*/
    double *T_nu;
};
typedef struct _transfer_init_table _transfer_init_table;

/** This function loads the initial transfer functions from CAMB transfer files.
 * It reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu and logk with length NPowerTable.
 * @param t_init Structure to initialise with the transfer function table.
 * @param BoxSize Size of simualtion box. Used to set maximal transfer function value to store.
 * @param UnitLength_in_cm Units of the stored transfer function in cm/h. Should be those expected by the simulation code, usually kpc/h.
 * @param InputSpectrum_UnitLength_in_cm Units of the CAMB transfer function in cm/h. Usually Mpc/h.
 * @param KspaceTransferFunction pointer to string containing the transfer function filename.*/
void allocate_transfer_init_table(_transfer_init_table *t_init, const double BoxSize, const double UnitLength_in_cm, const double InputSpectrum_UnitLength_in_cm, const char * KspaceTransferFunction);

/**Free memory for a transfer function table*/
void free_transfer_init_table(_transfer_init_table *t_init);

/*TRANSFER_INIT_H*/
#endif
