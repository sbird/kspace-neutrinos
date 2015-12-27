/* These functions will normally be defined by gadget. However, for the purposes of 
 * the standalone tests, we want our own definitions*/
#include <stdio.h>
#include <stdlib.h>
/*Forward define the hubble function*/
//double hubble_function(double a);

//Forward define terminate, because we'll need it.
void terminate(const char * string)
{
    fprintf(stderr, "Error: %s\n",string);
    exit(1);
}

void * mymalloc(const char * string, size_t size)
{
    return malloc(size);
}

void myfree(void * ptr)
{
    free(ptr);
}
