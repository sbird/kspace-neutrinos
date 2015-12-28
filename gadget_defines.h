#ifndef KSPACE_NEUTRINOS_PRIVATE
#define KSPACE_NEUTRINOS_PRIVATE
/* Header file for defines and structures which would normally be included in gadget's proto.h or allvars.c
 * When we use this code as part of gadget, these functions are defined in other translation units.
 * When used with the test suite, we define them in kspace_neutrinos_private.c.
 * Take care that the constants are in sync with the rest of gadget!*/

/*Speed of light in cm/s: in allvars.h this is called 'C'*/
#define  LIGHTCGS           2.997926e10

#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */


/* Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 *The correction is (3.046/3)^(1/4), for N_eff = 3.046 */
#define TNU     (T_CMB0*pow(4/11.,1/3.)*1.00381)              /* Neutrino + antineutrino background temperature in Kelvin */

#define  GRAVITY     6.672e-8 /*Newton's constant in cgs*/
#define  HUBBLE          3.2407789e-18	/* 100 km/s in h/sec */

/*Forward define the hubble function*/
double hubble_function(double a);

//Forward define terminate, because we'll need it.
void terminate(const char *);

void * mymalloc(const char *, size_t size);
void myfree(void * ptr);

//These are defined in begrun.c
#define STRING 2
#define REAL 1

#endif //KSPACE_NEUTRINOS_FUNC
