/*Global header, to be included in gadget's pm code. This defines the external interface to the kspace neutrino code*/
#ifndef KSPACE_NEUTRINOS_GLOBAL
#define KSPACE_NEUTRINOS_GLOBAL

#ifdef KSPACE_NEUTRINOS_2

#ifdef NEUTRINOS
#error "Cannot define particle based and Fourier-space neutrinos at the same time!"
#endif

#ifdef KSPACE_NEUTRINOS
#error "KSPACE_NEUTRINOS_2 is incompatible with KSPACE_NEUTRINOS"
#endif

/*We only need this for fftw_complex*/
#ifdef NOTYPEPREFIX_FFTW
#include        <fftw.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <dfftw.h>	/* double precision FFTW */
#else
#include     <sfftw.h>
#endif
#endif

#include <mpi.h>

/*This needs to change to match gadget.*/
#define MYMPI_COMM_WORLD MPI_COMM_WORLD

/* Return the total matter density in all neutrino species.*/
double OmegaNu(double a);

/* This sets up various structures for the kspace neutrinos, allocates memory, and reads saved data from disc. */
void allocate_kspace_memory(const int nk_in, const int ThisTask,const double BoxSize, const double UnitTime_in_s, const double UnitLength_in_cm, const double Omega0, const double HubbleParam, const char * snapdir, const double Time);

/* Main function, called from pm_periodic.c. 
   Computes the neutrino power, then adds it to the Fourier grid.*/
void add_nu_power_to_rhogrid(int save, const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int PMGRID, int ThisTask, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir, const double total_mass);

/* Function which sets up the parameter reader to read kspace neutrino parameters from the parameter file. 
 * It will store them in a static variable, kspace_params, in the translation unit where the function is defined
 * (which is the same as the above functions). */
int set_kspace_vars(char * tag[], void *addr[], int id [], int nt);

#endif //KSPACE_NEUTRINOS_2

#endif //KSPACE_NEUTRINOS_GLOBAL
