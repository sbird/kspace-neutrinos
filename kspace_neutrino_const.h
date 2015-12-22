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

/*With slightly relativistic massive neutrinos, for consistency we need to include radiation.
 * A note on normalisation (as of 08/02/2012):
 * CAMB appears to set Omega_Lambda + Omega_Matter+Omega_K = 1,
 * calculating Omega_K in the code and specifying Omega_Lambda and Omega_Matter in the paramfile.
 * This means that Omega_tot = 1+ Omega_r + Omega_g, effectively
 * making h0 (very) slightly larger than specified.
 */
#ifndef STEFAN_BOLTZMANN
#define STEFAN_BOLTZMANN 5.670373e-5
#endif

#ifdef INCLUDE_RADIATION
/*Stefan-Boltzmann constant in cgs units*/
/* Omega_g = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H^2)*/
#define OMEGAG (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*C*C*C*HUBBLE*HUBBLE)*pow(T_CMB0,4)/All.HubbleParam/All.HubbleParam)
#if (defined KSPACE_NEUTRINOS_2) || (defined NEUTRINOS)
    /*Neutrinos are included elsewhere*/
    #define OMEGAR OMEGAG
#else
    /*Neutrinos are included in the radiation*/
    /*For massless neutrinos, rho_nu/rho_g = 7/8 (T_nu/T_cmb)^4 *N_eff, but we absorbed N_eff into T_nu above*/
    #define OMEGANU (OMEGAG*7/8.*pow(TNU/T_CMB0,4)*3)
    /*With massless neutrinos only, add the neutrinos to the radiation*/
    #define OMEGAR (OMEGAG+OMEGANU)
#endif
#else
        /*Default is no radiation*/
        #define OMEGAR 0.
#endif


#ifndef MYMPI_COMM_WORLD
#define MYMPI_COMM_WORLD MPI_COMM_WORLD
#endif

/*Number of bins in integrations*/
#define GSL_VAL 200


#endif
