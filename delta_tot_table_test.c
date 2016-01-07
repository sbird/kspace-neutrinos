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
void init_hubble_function(_omega_nu * omnu, const double UnitTime_in_s);

/*A struct to hold some useful pointers*/
struct _test_state {
    _delta_pow * d_pow;
    _transfer_init_table * transfer;
    _omega_nu * omnu;
};
typedef struct _test_state test_state;

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
    test_state * ts = (test_state *) *state;
    _delta_pow * d_pow = (_delta_pow *) ts->d_pow;
    _omega_nu * omnu = (_omega_nu *) ts->omnu;
    _transfer_init_table * transfer = (_transfer_init_table *) ts->transfer;
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
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Note that the delta_cdm_curr used is actually at z=2, so the transfer function isn't really right, but who cares.*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, transfer, omnu, UnitTime_in_s, UnitLength_in_cm);
    assert_true(d_tot.ia == 25);
    /*Check we initialised delta_nu_init and delta_nu_last*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(d_tot.delta_nu_init[ik] > 0);
        assert_true(d_tot.delta_nu_last[ik] > 0);
    }
    /*Check saving works: we should also have saved a table in delta_tot_init, but save one in the test data directory and try to load it again.*/
    save_all_nu_state(&d_tot, "testdata/");
    _delta_tot_table d_tot2;
    allocate_delta_tot_table(&d_tot2, d_pow->nbins, 0.01, 1);
    read_all_nu_state(&d_tot2, "testdata/", 0.33333333);
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
    test_state * ts = (test_state *) *state;
    _delta_pow * d_pow = (_delta_pow *) ts->d_pow;
    _omega_nu * omnu = (_omega_nu *) ts->omnu;
    _transfer_init_table * transfer = (_transfer_init_table *) ts->transfer;
    _delta_tot_table d_tot;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1);
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Note that the delta_cdm_curr used is actually at z=2, so the transfer function isn't really right, but who cares.*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, transfer, omnu, UnitTime_in_s, UnitLength_in_cm);
    assert_true(d_tot.ia == 1);
    assert_true(d_tot.scalefact[0] == log(0.01));
    /*Check the initial power spectra were created properly*/
    const double OmegaNua3=get_omega_nu(omnu, 0.01)*pow(0.01,3);
    const double OmegaMa = omnu->Omega0 - get_omega_nu(omnu, 1) + OmegaNua3;
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

/*Check that the fits to the special J function are accurate, by doing explicit integration.*/
static void test_specialJ(void **state)
{
    /*Check against mathematica computed values:
    Integrate[(Sinc[q*x])*(q^2/(Exp[q] + 1)), {q, 0, Infinity}]*/
    assert_true(specialJ(0,-1) == 1);
    assert_true(fabs(specialJ(1,-1) - 0.2117) < 1e-3);
    assert_true(fabs(specialJ(2,-1) - 0.0223807) < 1e-3);
    assert_true(fabs(specialJ(0.5,-1) - 0.614729) < 1e-3);
    assert_true(fabs(specialJ(0.3,-1) - 0.829763) < 1e-3);
}

/* Check that we accurately work out the free-streaming length.
 * Free-streaming length for a non-relativistic particle of momentum q = T0, from scale factor ai to af.
 * The 'light' argument defines the units.
 * Test values use the following mathematica snippet:
 kB = 8.61734*10^(-5);
 Tnu = 2.7255*(4/11.)^(1/3.)*1.00328;
 omegar = 5.04672*10^(-5);
 Hubble[a_] := 3.085678*10^(21)/10^5*3.24077929*10^(-18)*Sqrt[0.2793/a^3 + (1 - 0.2793) + omegar/a^4]
  fs[a_, Mnu_] := kB*Tnu/(a*Mnu)/(a*Hubble[a])
  fslength[ai_, af_, Mnu_] := 299792*NIntegrate[fs[Exp[loga], Mnu], {loga, Log[ai], Log[af]}]
 */
static void test_fslength(void **state)
{
    /*Note that MNu is the mass of a single neutrino species:
     *we use large masses so that we don't have to compute omega_nu in mathematica.*/
    assert_true(fabs(fslength(0.5, 1,0.45, 299792.)/ 1272.92 -1 ) < 1e-5);
    assert_true(fabs(fslength(0.1, 0.5,0.6, 299792.)/ 5427.8 -1 ) < 1e-5);
}

static void test_get_delta_nu_update(void **state)
{
    /*Initialise stuff*/
    test_state * ts = (test_state *) *state;
    _delta_pow * d_pow = (_delta_pow *) ts->d_pow;
    _omega_nu * omnu = (_omega_nu *) ts->omnu;
    _transfer_init_table * transfer = (_transfer_init_table *) ts->transfer;
    _delta_tot_table d_tot;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1);
    /* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
     * Must be called before delta_tot_init, or resuming wont work*/
    read_all_nu_state(&d_tot, "testdata/", 0.33333333);
    /*Then init delta_tot*/
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, transfer, omnu, UnitTime_in_s, UnitLength_in_cm);
    /*Check that we will actually do something*/
    assert_true(log(0.3333333)-d_tot.scalefact[d_tot.ia-1] > 1e-5);
    /*So now we have a fully initialised d_tot. Update it!*/
    double delta_nu_curr[d_pow->nbins];
    get_delta_nu_update(&d_tot, 0.33333333, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, delta_nu_curr);
    /*Check that we did indeed update*/
    assert_true(d_tot.ia == 25);
    /*Check that we get the same answer as the saved powerspectrum*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(fabs(delta_nu_curr[ik]/ d_pow->delta_nu_curr[ik] -1) < 1e-2);
    }
    /*Throw out the last stored power spectrum and try again, so that we will perform a save*/
    d_tot.ia--;
    get_delta_nu_update(&d_tot, 0.33333333, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, delta_nu_curr);
    assert_true(d_tot.ia == 25);
    /*Check that we get the same answer as the saved powerspectrum*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        /*Be a bit more generous with the error as we have fewer datapoints now*/
        assert_true(fabs(delta_nu_curr[ik]/ d_pow->delta_nu_curr[ik] -1) < 3e-2);
    }
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
    test_state * ts = (test_state *) malloc(sizeof(test_state));
    ts->omnu = malloc(sizeof(_omega_nu));
    ts->transfer = malloc(sizeof(_transfer_init_table));
    ts->d_pow = d_pow;
    const double MNu[3] = {0.15, 0.15, 0.15};
    init_omega_nu(ts->omnu, MNu, 0.2793, 0.01, 0.7);
    const double OmegaNua3=get_omega_nu(ts->omnu, 0.01)*pow(0.01,3);
    const double OmegaMa = ts->omnu->Omega0 - get_omega_nu(ts->omnu, 1) + OmegaNua3;
    const double fnu = OmegaNua3/OmegaMa;
    for (int ik = 0; ik < nbins; ik++){
         delta_cdm_curr[ik] = (delta_cdm_curr[ik] - fnu*delta_nu_curr[ik])/(1.-fnu);
    }
    /*Now initialise data structure*/
    init_delta_pow(d_pow, logkk, delta_nu_curr, delta_cdm_curr, nbins);
    const double UnitLength_in_cm = 3.085678e21;
    allocate_transfer_init_table(ts->transfer, 500, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", ts->omnu);
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Set up the global variables for the hubble function before we do anything else!*/
    init_hubble_function(ts->omnu, UnitTime_in_s);
    *state = (void *) ts;
    return 0;
}

static int teardown_delta_pow(void **state) {
    test_state * ts = (test_state *) *state;
    if(!ts)
        return 0;
    free(ts->omnu);
    free(ts->transfer);
    _delta_pow * d_pow = ts->d_pow;
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
        cmocka_unit_test(test_specialJ),
        cmocka_unit_test(test_fslength),
        cmocka_unit_test(test_get_delta_nu_update),
    };
    return cmocka_run_group_tests(tests, setup_delta_pow, teardown_delta_pow);
}
