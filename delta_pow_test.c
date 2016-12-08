#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdio.h>
#include <math.h>
#include "delta_pow.h"

static void test_init_delta_pow(void **state) {
    /*Load initialised structure*/
    _delta_pow * d_pow = (_delta_pow *) *state;
    /*Check it initialised correctly.*/
    assert_true(d_pow);
    assert_true(d_pow->logkk);
    assert_true(d_pow->delta_nu_curr);
    assert_true(d_pow->spline_nu);
    assert_true(d_pow->acc_nu);
    assert_true(d_pow->delta_cdm_curr);
    assert_true(d_pow->spline_cdm);
    assert_true(d_pow->acc);
}

static void test_get_dnudcdm_powerspec(void **state) {
    /*Load initialised structure*/
    _delta_pow * d_pow = (_delta_pow *) *state;
    /*Check it initialised correctly.*/
    assert_true(d_pow);
    /*Check that we get the desired answer at interpolation points*/
    for(int i=0; i<10; i++)
        assert_true(get_dnudcdm_powerspec(d_pow, d_pow->logkk[5*i]) == d_pow->delta_nu_curr[5*i]/d_pow->delta_cdm_curr[5*i]);
    /*Check that we get a reasonable answer in between*/
    for(int i=0; i<50; i++) {
        double logk = (d_pow->logkk[5*i+1] + d_pow->logkk[5*i+2])/2;
        double pk = (d_pow->delta_nu_curr[5*i+1]/d_pow->delta_cdm_curr[5*i+1] + d_pow->delta_nu_curr[5*i+2]/d_pow->delta_cdm_curr[5*i+2])/2;
        /*   if (!(fabs(get_dnudcdm_powerspec(d_pow, logk) - pk) < 1e-3*pk))
             printf("i=%d %g %g %g\n", i, logk, pk, get_dnudcdm_powerspec(d_pow, logk));*/
        /*Of course this is not as accurate as the cubic spline, but if linear and
         * spline interpolation give the same answer, the spline interpolation should be very accurate.*/
        assert_true(fabs(get_dnudcdm_powerspec(d_pow, logk) - pk) < 5e-3*pk);
    }
    /*Check that we get something if we are slightly past the edges of the interpolator*/
    double pksmall = d_pow->delta_nu_curr[0]/d_pow->delta_cdm_curr[0];
    assert_true(fabs(get_dnudcdm_powerspec(d_pow, d_pow->logkk[0]-0.01) - pksmall) < 5e-3*pksmall);
    double pklarge = d_pow->delta_nu_curr[d_pow->nbins-1]/d_pow->delta_cdm_curr[d_pow->nbins-1];
    assert_true(fabs(get_dnudcdm_powerspec(d_pow, d_pow->logkk[d_pow->nbins-1]+0.01) - pklarge) < 5e-3*pklarge);
}

static void test_save_nu_power(void **state) {
    /*Load initialised structure*/
    _delta_pow * d_pow = (_delta_pow *) *state;
    /*Check it initialised correctly.*/
    assert_true(d_pow);
    /*Save P_nu(k) to disc*/
    save_nu_power(d_pow, 1/3., 005, "./");
    /*Load it again and check it is the same.*/
    FILE * fd = fopen("testdata/powerspec_nu_004.txt", "r");
    assert_true(fd > 0);
    /*Read redshift*/
    double scale;
    assert_true(fscanf(fd, "%lg", &scale));
    assert_true(fabs(scale - 1./3) < 1e-5);
    /*Read number of bins*/
    int nbins;
    assert_true(fscanf(fd,"%d",&nbins));
    assert_true(nbins == d_pow->nbins);
    for(int i=0; i< nbins; i++){
        /*Read kvalues and P(k)*/
        double logkk, delta_nu_curr;
        assert_true(fscanf(fd, "%lg %lg", &logkk, &delta_nu_curr) == 2);
        assert_true(fabs(log(logkk) - d_pow->logkk[i]) < 1e-5*fabs(logkk));
        assert_true(fabs(sqrt(delta_nu_curr) - d_pow->delta_nu_curr[i]) < 1e-5*delta_nu_curr);
    }
    fclose(fd);
}


/*Test we can initialise a delta_pow structure from disc correctly.
 *Note if one of these assertions is false we won't actually get an error; just a message saying group setup failed.*/
static int setup_delta_pow(void **state) {
    /*Initialise the structure by reading something from the test file*/
    double * logkk;
    double * delta_nu_curr;
    double * delta_cdm_curr;
    _delta_pow * d_pow = NULL;
    /*Initialise the structure and use it as the state.*/
    d_pow = (_delta_pow *) malloc(sizeof(_delta_pow));
    FILE * fd;
    double scale;
    int nbins;
    /*First read the neutrino power*/
    if((fd = fopen("testdata/powerspec_nu_004.txt", "r"))){
        /*Read redshift*/
        assert_true(fscanf(fd, "%lg", &scale));
        /*Read number of bins*/
        assert_true(fscanf(fd,"%d",&nbins));
        logkk = (double *) malloc(nbins*sizeof(double));
        delta_nu_curr = (double *) malloc(nbins*sizeof(double));
        delta_cdm_curr = (double *) malloc(nbins*sizeof(double));
        if(!delta_nu_curr || !delta_cdm_curr){
            printf("Failed to allocate memory. nbins=%d\n", nbins);
            return 1;
        }
        for(int i=0; i< nbins; i++){
            /*Read kvalues and P(k)*/
            assert_true(fscanf(fd, "%lg %lg", logkk+i, delta_nu_curr+i) == 2);
        }
        fclose(fd);
    }
    else {
        printf("Could not open powerspec_nu_004.txt\n");
        return 1;
    }

    /*Convert the units*/
    for(int i=0; i < nbins; i++){
        logkk[i] = log(logkk[i]);
        delta_nu_curr[i] = sqrt(delta_nu_curr[i]);
    }
    /*Then read CDM power. We use the linear power spectrum, just because it's easier.*/
    if((fd = fopen("testdata/powerspec_cdm_004.txt", "r"))){
        /*Read redshift*/
        int nbins2;
        double dummy;
        assert_true(fscanf(fd, "%lg", &dummy));
        /*Read number of bins*/
        assert_true(fscanf(fd,"%d",&nbins2));
        assert_true(nbins2 == nbins);
        for(int i=0; i< nbins; i++){
            assert_true(fscanf(fd, "%lg", delta_cdm_curr+i));
        }
        fclose(fd);
    }
    else {
        printf("Could not open powerspec_cdm_004.txt\n");
        return 1;
    }
    /*Now initialise data structure*/
    init_delta_pow(d_pow, logkk, delta_nu_curr, delta_cdm_curr, nbins);
    *state = (void *) d_pow;
    return 0;
}


static int teardown_delta_pow(void **state) {
    _delta_pow * d_pow = (_delta_pow *) *state;
    if(!d_pow)
        return 0;
    free_d_pow(d_pow);
    free(d_pow->delta_nu_curr);
    free(d_pow->delta_cdm_curr);
    free(d_pow);
    return 0;
}

int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_init_delta_pow),
        cmocka_unit_test(test_get_dnudcdm_powerspec),
        cmocka_unit_test(test_save_nu_power),
    };
    return cmocka_run_group_tests(tests, setup_delta_pow, teardown_delta_pow);
}
