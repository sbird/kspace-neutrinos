#ifndef KSPACE_NEUTRINO_CONSTANTS
#define KSPACE_NEUTRINO_CONSTANTS
#ifndef T_CMB0 
#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */
#endif

#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/
#define HBAR    6.582119e-16  /*hbar in units of eV s*/
#define FLOAT   1e-6            /*Floating point accuracy*/
#ifndef C
#define  C           2.9979e10 /*Speed of light in cm/s*/
#endif
#ifndef GRAVITY
#define  GRAVITY     6.672e-8 /*Newton's constant in cgs*/
#endif
#ifndef HUBBLE
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#endif


#ifndef TNU 
/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 *The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)              /* Neutrino + antineutrino background temperature in Kelvin */
#endif
#ifndef MYMPI_COMM_WORLD
#define MYMPI_COMM_WORLD MPI_COMM_WORLD
#endif

/*Number of bins in integrations*/
#define GSL_VAL 200


#endif
