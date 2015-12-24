#ifndef KSPACE_NEUTRINOS_PRIVATE
#define KSPACE_NEUTRINOS_PRIVATE
/* Header file for functions and structures which would normally be included in gadget's proto.h
 * When we use this code as part of gadget, these functions are defined in other translation units.
 * When used with the test suite, we define them in kspace_neutrinos_private.c.*/

/*Forward define the hubble function*/
double hubble_function(double a);

//Forward define terminate, because we'll need it.
void terminate(const char *);

void * mymalloc(const char *, size_t size);
void myfree(void * ptr);

#ifndef MYMPI_COMM_WORLD
#define MYMPI_COMM_WORLD MPI_COMM_WORLD
#endif

#endif //KSPACE_NEUTRINOS_FUNC
