/* These functions will normally be defined by gadget. However, for the purposes of
 * the standalone tests, we want our own definitions*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gadget_defines.h"
#include "omega_nu_single.h"

/*These functions need bodies; normally this is provided by gadget*/
void * mymalloc_fullinfo(const char * string, size_t size, const char *func, const char *file, int line)
{
    return malloc(size);
}

void myfree_fullinfo(void * ptr, const char *func, const char *file, int line)
{
    free(ptr);
}

