#ifndef KSPACE_NEUTRINOS_FUNC
#define KSPACE_NEUTRINOS_FUNC
/*File to define globally accessible functions for massive neutrinos. This is the header that should be included in proto.h*/

/* Return the total matter density in all neutrino species.*/
double OmegaNu(double a);
/*The matter density in a single neutrino species */
double OmegaNu_single(double a,double mnu, int sp);


//Forward define terminate, because we'll need it.
void terminate(const char *);
#ifndef mymalloc
#define mymalloc(x,y) malloc(y)
#endif

#ifndef myfree
#define myfree(x) free(x)
#endif

#ifdef KSPACE_NEUTRINOS_2
#include <gsl/gsl_interp.h>

#ifdef NEUTRINOS
#error "Cannot define particle based and Fourier-space neutrinos at the same time!"
#endif

#ifdef KSPACE_NEUTRINOS
#error "KSPACE_NEUTRINOS_2 is incompatible with KSPACE_NEUTRINOS"
#endif

/*These functions only need to be around if we actually have kspace neutrinos. They are not needed for particle neutrinos*/
void get_delta_nu_update(double a, int nk_in, double wavenum[], double P_cdm_curr[], double delta_nu_curr[], int ThisTask);
double get_neutrino_powerspec(double kk_in, double SmoothK[], double SmoothPowerNu[],
			      gsl_interp * SplinePowNu, gsl_interp_accel * acc_nu, double SmoothPower[],
			      gsl_interp * SplinePow, gsl_interp_accel * acc, int nbins);
void transfer_init_tabulate(int nk_in, int ThisTask);
void rebin_power(double SmoothK[], double SmoothPower[], int nbins, double Kbin[], double SumPower[],
		 long long CountModes[], int bins_ps, int pmgrid);
void save_all_nu_state(char * savedir);
void read_all_nu_state(char * savedir, double Time);
#endif //KSPACE_NEUTRINOS_2

#endif //KSPACE_NEUTRINOS_FUNC
