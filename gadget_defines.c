/* These functions will normally be defined by gadget. However, for the purposes of
 * the standalone tests, we want our own definitions*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "gadget_defines.h"
#include "omega_nu_single.h"
#include "mpi.h"

extern int ThisTask;

/*These functions need bodies; normally this is provided by gadget*/
void terminate(int ierr, const char * fmt, ...)
{
    va_list va;
    char buf[4096];
    va_start(va, fmt);
    vsprintf(buf, fmt, va);
    va_end(va);
    printf("%s",buf);
    MPI_Abort(MPI_COMM_WORLD, ierr);
    exit(ierr);
}

/*  This function writes a message.
 *
 *  if ierr > 0, the message is uncollective.
 *  if ierr <= 0, the message is 'collective', only the root rank prints the message. A barrier is applied.
 */

void message(int ierr, const char * fmt, ...)
{
    if(ierr > 0 || ThisTask == 0) {
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

