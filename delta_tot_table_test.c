#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdio.h>
#include <math.h>
#include "delta_tot_table.h"
#include "delta_pow.h"
#include "transfer_init.h"
#include "omega_nu_single.h"
#include "gadget_defines.h"

/*Initialise the hubble function: this is defined in gadget_defines.c and forward declared here,
 *to avoid accidentally using it anywhere except in tests.*/
void init_hubble_function(const double MNu[], const double Omega0, const double a0, const double HubbleParam, const double UnitTime_in_s);

#if 0
/*Update the last value of delta_tot in the table with a new value computed
 from the given delta_cdm_curr and delta_nu_curr.
 If overwrite is true, overwrite the existing final entry.*/
void update_delta_tot(_delta_tot_table *d_tot, _omega_nu * omnu, double a, double delta_cdm_curr[], double delta_nu_curr[], int overwrite);

/*Function called by add_nu_power_to_rhogrid*/
void get_delta_nu_update(_delta_tot_table *d_tot, _omega_nu * omnu, double a, int nk_in, double keff[], double P_cdm_curr[], double delta_nu_curr[]);

/*This function does the work and updates delta_nu_curr*/
void get_delta_nu(_delta_tot_table *d_tot, double a, int Na, double wavenum[], double delta_nu_curr[],double mnu, double vcrit);

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(_delta_tot_table *d_tot, double a, int Na, double wavenum[],  double delta_nu_curr[], _omega_nu * omnu);

#endif

/* Test that the allocations are done correctly.
 * delta_tot is still empty (but allocated) after this.*/
static void test_allocate_delta_tot_table(void **state)
{
    _delta_tot_table d_tot;
    allocate_delta_tot_table(&d_tot, 300, 0.01, 1);
    assert_true(d_tot.ia == 0);
    assert_true(d_tot.namax > 10);
    assert_true(d_tot.scalefact);
    assert_true(d_tot.delta_nu_init);
    assert_true(d_tot.delta_nu_last);
    assert_true(d_tot.delta_tot);
    for(int i=0; i<d_tot.nk; i++){
        assert_true(d_tot.delta_tot[i]);
    }
    free_delta_tot_table(&d_tot);
}

static void test_save_resume(void **state)
{
    _delta_pow * d_pow = (_delta_pow *) *state;
    _delta_tot_table d_tot;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1);
    /* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
     * Must be called before delta_tot_init, or resuming wont work*/
    read_all_nu_state(&d_tot, "testdata/", 0.33333333);
    assert_true(d_tot.ia == 25);
    assert_true(fabs(d_tot.scalefact[0]/log(0.01)-1) < 1e-5);
    for(int i=1; i < d_tot.ia; i++) {
        assert_true(d_tot.scalefact[i] > d_tot.scalefact[i-1]);
        for(int kk=0; kk < d_tot.nk; kk++)
            assert_true(d_tot.delta_tot[kk][i] > 0);
    }
    /*Now check that calling delta_tot_init after this works.*/
    /*Set up the background*/
    _omega_nu omnu;
    const double MNu[3] = {0.15, 0.15, 0.15};
    init_omega_nu(&omnu, MNu, 0.2793, 0.01, 0.7);
    /*Set up the transfer table*/
    _transfer_init_table transfer;
    const double UnitLength_in_cm = 3.085678e21;
    allocate_transfer_init_table(&transfer, 500, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", &omnu);
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Note that the delta_cdm_curr used is actually at z=2, so the transfer function isn't really right, but who cares.*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, &transfer, &omnu, UnitTime_in_s, UnitLength_in_cm);
    assert_true(d_tot.ia == 25);
    /*Check we initialised delta_nu_init and delta_nu_last*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(d_tot.delta_nu_init[ik] > 0);
        assert_true(d_tot.delta_nu_last[ik] > 0);
    }
    /*Check saving works: we should have saved a table in delta_tot_init, so try to load it again.*/
    _delta_tot_table d_tot2;
    allocate_delta_tot_table(&d_tot2, d_pow->nbins, 0.01, 1);
    read_all_nu_state(&d_tot2, NULL, 0.33333333);
    assert_true(d_tot.ia == d_tot2.ia);
    assert_true(d_tot.nk == d_tot2.nk);
    for(int i=1; i < d_tot.ia; i++) {
        assert_true(d_tot.scalefact[i] == d_tot2.scalefact[i]);
        for(int kk=0; kk < d_tot.nk; kk++)
            assert_true(d_tot.delta_tot[kk][i] == d_tot2.delta_tot[kk][i]);
    }
}

/* Test that we can initialise delta_tot without resuming.*/
static void test_delta_tot_init(void **state)
{
    _delta_pow * d_pow = (_delta_pow *) *state;
    _delta_tot_table d_tot;
    /*Set up the background*/
    _omega_nu omnu;
    const double MNu[3] = {0.15, 0.15, 0.15};
    init_omega_nu(&omnu, MNu, 0.2793, 0.01, 0.7);
    /*Set up the transfer table*/
    _transfer_init_table transfer;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1);
    const double UnitLength_in_cm = 3.085678e21;
    allocate_transfer_init_table(&transfer, 500, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", &omnu);
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Note that the delta_cdm_curr used is actually at z=2, so the transfer function isn't really right, but who cares.*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, &transfer, &omnu, UnitTime_in_s, UnitLength_in_cm);
    assert_true(d_tot.ia == 1);
    assert_true(d_tot.scalefact[0] == log(0.01));
    /*Check the initial power spectra were created properly*/
    const double OmegaNua3=get_omega_nu(&omnu, 0.01)*pow(0.01,3);
    const double OmegaMa = omnu.Omega0 - get_omega_nu(&omnu, 1) + OmegaNua3;
    const double fnu = OmegaNua3/OmegaMa;
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(d_tot.delta_nu_init[ik] > 0);
        assert_true(d_tot.delta_nu_last[ik] > 0);
        /*These two should be initially the same, although one is created by calling the integrator.*/
        assert_true(fabs(d_tot.delta_nu_last[ik]/ d_tot.delta_nu_init[ik] -1) < 1e-4);
        double delta_tot_before = (1.-fnu)* d_pow->delta_cdm_curr[ik] + fnu*d_tot.delta_nu_init[ik];
        assert_true(fabs(d_tot.delta_tot[ik][0]/delta_tot_before-1) < 1e-4);
    }
    free_delta_tot_table(&d_tot);
}

/*We are using delta_pow as a source for the current state of the integrator*/
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
        printf("Could not open testdata/powerspec_nu_004.txt\n");
        return 1;
    }

    /*Convert the units*/
    for(int i=0; i < nbins; i++){
        /*We want k, not logk for this*/
        /*logkk[i] = log(logkk[i]);*/
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
        printf("Could not open testdata/powerspec_cdm_004.txt\n");
        return 1;
    }
    /*We have actually read the total matter power; so for the CDM power we should subtract off the neutrino power*/
    _omega_nu omnu;
    const double MNu[3] = {0.15, 0.15, 0.15};
    init_omega_nu(&omnu, MNu, 0.2793, 0.01, 0.7);
    const double OmegaNua3=get_omega_nu(&omnu, 0.01)*pow(0.01,3);
    const double OmegaMa = omnu.Omega0 - get_omega_nu(&omnu, 1) + OmegaNua3;
    const double fnu = OmegaNua3/OmegaMa;
    for (int ik = 0; ik < nbins; ik++){
         delta_cdm_curr[ik] = (delta_cdm_curr[ik] - fnu*delta_nu_curr[ik])/(1.-fnu);
    }
    /*Now initialise data structure*/
    init_delta_pow(d_pow, logkk, delta_nu_curr, delta_cdm_curr, nbins);
    *state = (void *) d_pow;
    const double UnitLength_in_cm = 3.085678e21;
//     _transfer_init_table transfer;
//     allocate_transfer_init_table(&transfer, 500, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", &omnu);
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Set up the global variables for the hubble function before we do anything else!*/
    init_hubble_function(MNu, 0.2793, 0.01, 0.7, UnitTime_in_s);

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
        cmocka_unit_test(test_allocate_delta_tot_table),
        cmocka_unit_test(test_save_resume),
        cmocka_unit_test(test_delta_tot_init),
    };
    return cmocka_run_group_tests(tests, setup_delta_pow, teardown_delta_pow);
}
