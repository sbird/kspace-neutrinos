/* These functions will normally be defined by gadget. However, for the purposes of
 * the standalone tests, we want our own definitions*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "gadget_defines.h"
#include "omega_nu_single.h"

/*These functions need bodies; normally this is provided by gadget*/
void terminate(int ierr, const char * fmt, ...)
{
    va_list va;
    va_start(va, fmt);
    printf(fmt, va);
    va_end(va);
#ifdef MPI_VERSION
    MPI_Abort(MPI_COMM_WORLD, ierr);
#endif
    exit(ierr);
}

/*  This function writes a message.
 *
 *  if ierr > 0, the message is uncollective.
 *  if ierr <= 0, the message is 'collective', only the root rank prints the message. A barrier is applied.
 */

void message(int ierr, const char * fmt, ...)
{
    int rank = 0;
#ifdef MPI_VERSION
    if(ierr == 0)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if(ierr > 0 || rank == 0) {
        va_list va;
        char buf[4096];
        va_start(va, fmt);
        vsprintf(buf, fmt, va);
        va_end(va);
        printf("%s", buf);
    }
}

void * mymalloc_fullinfo(const char * string, size_t size, const char *func, const char *file, int line)
{
    return malloc(size);
}

void myfree_fullinfo(void * ptr, const char *func, const char *file, int line)
{
    free(ptr);
}

