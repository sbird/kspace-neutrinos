#ifndef KSPACE_NEUTRINOS_PRIVATE
#define KSPACE_NEUTRINOS_PRIVATE
/* Header file for defines and structures which would normally be included in gadget's proto.h or allvars.c
 * When we use this code as part of gadget, these functions are defined in other translation units.
 * When used with the test suite, we define them in kspace_neutrinos_private.c.
 * Take care that the constants are in sync with the rest of gadget!*/

/*Speed of light in cm/s: in allvars.h this is called 'C'*/
#define  LIGHTCGS           2.99792458e10
#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */

/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 * We use the CLASS default value, chosen so that omega_nu = m_nu / 93.14 h^2
 * At time of writing this is T_nu / T_gamma = 0.71611.
 * See https://github.com/lesgourg/class_public/blob/master/explanatory.ini
 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00328)              /* Neutrino + antineutrino background temperature in Kelvin */

#define  GRAVITY     6.67408e-8 /*Newton's constant in cgs*/
#define  HUBBLE          3.24077929e-18	/* 100 km/s in h/sec */

/*Forward define the hubble function*/
double hubble_function(double a);

#include <stdlib.h>

#ifdef MPI_VERSION
#define  terminate(x) {fprintf(stderr,"code termination, function '%s()', file '%s', line %d: '%s'\n", __FUNCTION__, __FILE__, __LINE__, x); MPI_Abort(MPI_COMM_WORLD, 1); exit(1);}
#else
#define  terminate(x) {fprintf(stderr,"code termination, function '%s()', file '%s', line %d: '%s'\n", __FUNCTION__, __FILE__, __LINE__, x); exit(1);}
#endif

/*Definitions from gadget's allvars.h: these are macros, so we have to repeat them here or include allvars.h.
 If you do define DISABLE_MEMORY_MANAGER in Gadget's Config.sh you will need to define it again in our Makefile.*/
#ifndef DISABLE_MEMORY_MANAGER
    #define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
    #define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
    /*These functions need bodies; normally this is provided by gadget. We fake it in gadget_defines.c for the tests.*/
    void * mymalloc_fullinfo(const char * string, size_t size, const char *func, const char *file, int line);
    void myfree_fullinfo(void * ptr, const char *func, const char *file, int line);
#else
    #define  mymalloc(x, y)            malloc(y)
    #define  myfree(x)                 free(x)
#endif

//These are defined in begrun.c
#define STRING 2
#define REAL 1

#endif //KSPACE_NEUTRINOS_FUNC
