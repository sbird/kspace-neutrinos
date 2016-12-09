#ifndef KSPACE_NEUTRINOS_PRIVATE
#define KSPACE_NEUTRINOS_PRIVATE
/* Header file for defines and structures which would normally be included in gadget's proto.h or allvars.c
 * When we use this code as part of gadget, these functions are defined in other translation units.
 * When used with the test suite, we define them in kspace_neutrinos_private.c.
 * Take care that the constants are in sync with the rest of gadget!*/

#define  HUBBLE          3.24077929e-18	/* 100 km/s in h/sec */

/*Forward define the hubble function*/
double hubble_function(double a);

#include <stdlib.h>

#ifdef MPI_VERSION
#define  terminate(x) {fprintf(stderr,"code termination, function '%s()', file '%s', line %d: '%s'\n", __FUNCTION__, __FILE__, __LINE__, x); MPI_Abort(MPI_COMM_WORLD, 1); exit(1);}
#else
#define  terminate(x) {fprintf(stderr,"code termination, function '%s()', file '%s', line %d: '%s'\n", __FUNCTION__, __FILE__, __LINE__, x); exit(1);}
#endif

/*Definitions from gadget's allvars.h: these are macros, so we have to repeat them here or include allvars.h.*/
#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
/*These functions need bodies; normally this is provided by gadget. We fake it in gadget_defines.c for the tests.*/
void * mymalloc_fullinfo(const char * string, size_t size, const char *func, const char *file, int line);
void myfree_fullinfo(void * ptr, const char *func, const char *file, int line);

/*These are defined in begrun.c. They correspond to types and are very specific to P-Gadget3!
 * Only used in set_kspace_vars*/
#define INT 3
#define STRING 2
#define REAL 1

/*KSPACE_NEUTRINOS_PRIVATE*/
#endif
