#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdio.h>
#include <math.h>
#include "powerspectrum.h"
#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw.h>	/* double precision FFTW */
#else
#include     <srfftw.h>
#endif
#endif

/*Test the total powerspectrum on one processor only*/
static void test_total_powerspectrum(void **state) {
    (void) state;
    fftw_real field[4*4*2*(4/2+1)]={0};
    double pow[10];
    long long int count[10];
    double keffs[10];
    fftw_complex * outfield;
    outfield = (fftw_complex *) &field[0];
    //Pad the input
    for(int i=0; i<32; i++) {
        int ii = i/4;
        field[6*ii+i%4]=1;
    }
    field[0]=2;
    rfftwnd_plan pl = rfftw3d_create_plan(4,4,4,FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    double total_mass = 4*4*4;
    rfftwnd_one_real_to_complex(pl, &field[0], outfield);
    /* Compute the total powerspectrum from a Fourier-transformed density field in outfield, and store it in power.
     * Before use you may wish to normalise by dividing by count*count*/
    total_powerspectrum(4,&outfield[0],10,0, 4, pow,count,keffs, total_mass);
    assert_true(fabs(keffs[2]-1.73205) < 1e-5);
    assert_true(count[1]==12);
    assert_true(count[0]==6);
    assert_true(fabs(pow[0]-0.0677526) < 1e-5*0.04);
    assert_true(fabs(pow[1]-0.000565561) < 1e-5*0.005);
    assert_true(fabs(pow[2]-0.000860793) < 1e-5*0.003);
    assert_true(fabs(pow[9]) < 1e-6);
    rfftwnd_destroy_plan(pl);
}

static int setup_mpi(void **state) {
    int ac=1;
    char * str = "powerspectrum_test";
    char **av = &str;
    MPI_Init(&ac, &av);
    return 0;
}

static int teardown_mpi(void **state) {
    MPI_Finalize();
    return 0;
}


int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_total_powerspectrum),
    };
    return cmocka_run_group_tests(tests, setup_mpi, teardown_mpi);
}
