#ifndef KSPACE_NEUTRINOS_FUNC
#define KSPACE_NEUTRINOS_FUNC
/*File to define globally accessible functions for massive neutrinos. This is the header that should be included in proto.h*/

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

/* Return the total matter density in all neutrino species.*/
double OmegaNu(double a);

#ifdef KSPACE_NEUTRINOS_2

#ifdef NEUTRINOS
#error "Cannot define particle based and Fourier-space neutrinos at the same time!"
#endif

#ifdef KSPACE_NEUTRINOS
#error "KSPACE_NEUTRINOS_2 is incompatible with KSPACE_NEUTRINOS"
#endif

/* These functions only need to be around if we actually have kspace neutrinos. They are not needed for particle neutrinos*/
/* Main function, called from pm_periodic.c. 
   Computes the neutrino power, then adds it to the Fourier grid.*/
void add_nu_power_to_rhogrid(int save, const double Time, const double Omega0, const double BoxSize, fftw_complex *fft_of_rhogrid, const int PMGRID, int ThisTask, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir, const double total_mass);

/*Functions to load data for the neutrino powerspectrum from the disc*/
void transfer_init_tabulate(int nk_in, int ThisTask);

/*Forward define the hubble function*/
double hubble_function(double a);

#endif //KSPACE_NEUTRINOS_2

#endif //KSPACE_NEUTRINOS_FUNC
