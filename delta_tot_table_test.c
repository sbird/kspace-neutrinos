#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "delta_tot_table.h"
#include "delta_pow.h"
#include "transfer_init.h"
#include "omega_nu_single.h"
#include "kspace_neutrino_const.h"
#include "gadget_defines.h"

#define  T_CMB0      2.7255	/* present-day CMB temperature, from Fixsen 2009 */
/* This is a bit goofy, because the hubble function in gadget uses global state.
 * So in the tests we need to set the global state with init_hubble_function before using it,
 * to keep the function definition the same.
 * We COULD just use an internal version, but then we would break compatibility if anyone has
 * used an odd hubble history for eg, some dark energy model. */
static _omega_nu * m_omnu;
static double Omega_nonu;
static double OmegaLambda;
static double m_Hubble;

void init_hubble_function(_omega_nu * omnu, const double Omega0, const double UnitTime_in_s)
{
    m_omnu = omnu;
    Omega_nonu = Omega0 - get_omega_nu(omnu, 1);
    OmegaLambda = 1 - Omega0;
    m_Hubble = HUBBLE * UnitTime_in_s;
}

double hubble_function(double a)
{
    if(!m_omnu) {
        terminate(1,"init_hubble_function was not called in test suite before this!\n");
    }
    /* Matter + Lambda: neglect curvature*/
    double omega_tot = Omega_nonu/pow(a,3) + OmegaLambda;
    /*Neutrinos*/
    omega_tot += get_omega_nu(m_omnu, a);
    /*Radiation*/
    omega_tot += get_omegag(m_omnu, a);
    return m_Hubble * sqrt(omega_tot);
}

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
    test_state * ts = (test_state *) *state;
    _delta_pow * d_pow = (_delta_pow *) ts->d_pow;
    _transfer_init_table * transfer = (_transfer_init_table *) ts->transfer;
    _delta_tot_table d_tot;
    _omega_nu omnu;
    double MNu[3] = {0, 0, 0};
    init_omega_nu(&omnu, MNu, 0.01, 0.7,T_CMB0);
    allocate_delta_tot_table(&d_tot, 300, 0.01, 1, 0.2793, &omnu, 1, 1, 0);
    assert_true(d_tot.ia == 0);
    assert_true(d_tot.namax > 10);
    assert_true(d_tot.scalefact);
    assert_true(d_tot.delta_nu_init);
    assert_true(d_tot.delta_nu_last);
    assert_true(d_tot.delta_tot);
    for(int i=0; i<d_tot.nk_allocated; i++){
        assert_true(d_tot.delta_tot[i]);
    }
    /* Check that we do not crash when neutrino mass is zero*/
    double delta_nu_curr[d_pow->nbins];
    get_delta_nu_update(&d_tot, 0.02, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, delta_nu_curr, transfer);
    free_delta_tot_table(&d_tot);
}

static void test_save_resume(void **state)
{
    test_state * ts = (test_state *) *state;
    _delta_pow * d_pow = (_delta_pow *) ts->d_pow;
    _omega_nu * omnu = (_omega_nu *) ts->omnu;
    _transfer_init_table * transfer = (_transfer_init_table *) ts->transfer;
    _delta_tot_table d_tot;
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1, 0.2793, omnu, UnitTime_in_s, UnitLength_in_cm, 0);
    /* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
     * Must be called before delta_tot_init, or resuming wont work*/
    read_all_nu_state(&d_tot, "testdata/delta_tot_nu.txt");
    assert_true(d_tot.ia == 25);
    assert_true(fabs(d_tot.scalefact[0]/log(0.01)-1) < 1e-5);
    for(int i=1; i < d_tot.ia; i++) {
        assert_true(d_tot.scalefact[i] > d_tot.scalefact[i-1]);
        for(int kk=0; kk < d_tot.nk; kk++)
            assert_true(d_tot.delta_tot[kk][i] > 0);
    }
    /*Now check that calling delta_tot_init after this works.*/
    /*Note that the delta_cdm_curr used is actually at z=2, so the transfer function isn't really right, but who cares.*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, transfer, 0.3333333);
    assert_true(d_tot.ia == 25);
    /*Check we initialised delta_nu_init and delta_nu_last*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(d_tot.delta_nu_init[ik] > 0);
        assert_true(d_tot.delta_nu_last[ik] > 0);
    }
    /*Check saving works: we should also have saved a table in delta_tot_init, but save one in the test data directory and try to load it again.*/
    save_all_nu_state(&d_tot, "testdata/delta_tot_nu.txt");
    _delta_tot_table d_tot2;
    allocate_delta_tot_table(&d_tot2, d_pow->nbins, 0.01, 1, 0.2793, omnu, UnitTime_in_s, UnitLength_in_cm, 0);
    read_all_nu_state(&d_tot2, "testdata/delta_tot_nu.txt");
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
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1, 0.2793, omnu, UnitTime_in_s, UnitLength_in_cm, 0);
    /*Note that the delta_cdm_curr used is actually at z=2, so the transfer function isn't really right, but who cares.*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, transfer,0.01);
    assert_true(d_tot.ia == 1);
    assert_true(d_tot.scalefact[0] == log(0.01));
    /*Check the initial power spectra were created properly*/
    const double OmegaNua3=get_omega_nu(omnu, 0.01)*pow(0.01,3);
    const double OmegaNu1=get_omega_nu(omnu, 1);
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(d_tot.delta_nu_init[ik] > 0);
        assert_true(d_tot.delta_nu_last[ik] > 0);
        /*These two should be initially the same, although one is created by calling the integrator.*/
        assert_true(fabs(d_tot.delta_nu_last[ik]/ d_tot.delta_nu_init[ik] -1) < 1e-4);
        double delta_tot_before = get_delta_tot(d_tot.delta_nu_init[ik],d_pow->delta_cdm_curr[ik],OmegaNua3,d_tot.Omeganonu, OmegaNu1);
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
    /*Test that it is ok when truncated*/
    /*Mathematica: Jfrac[x_, qc_] := NIntegrate[(Sinc[q*x])*(q^2/(Exp[q] + 1)), {q, qc, Infinity}]/(3*Zeta[3]/2) */
    assert_true(fabs(specialJ(0,1) - 0.940437) < 1e-4);
    assert_true(fabs(specialJ(0.5,1) - 0.556557) < 1e-4);
    assert_true(fabs(specialJ(1,0.1) - 0.211611) < 1e-4);
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
    double kT = BOLEVK*TNUCMB*T_CMB0;
    /*fslength function returns fslength * (MNu / kT)*/
    assert_true(fabs(fslength(log(0.5), log(1), 299792.)/ 1272.92/(0.45/kT) -1 ) < 1e-5);
    assert_true(fabs(fslength(log(0.1), log(0.5),299792.)/ 5427.8/(0.6/kT) -1 ) < 1e-5);
}

static void test_get_delta_nu_update(void **state)
{
    /*Initialise stuff*/
    test_state * ts = (test_state *) *state;
    _delta_pow * d_pow = (_delta_pow *) ts->d_pow;
    _omega_nu * omnu = (_omega_nu *) ts->omnu;
    _transfer_init_table * transfer = (_transfer_init_table *) ts->transfer;
    _delta_tot_table d_tot;
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    allocate_delta_tot_table(&d_tot, d_pow->nbins, 0.01, 1, 0.2793, omnu, UnitTime_in_s, UnitLength_in_cm, 0);
    /* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
     * Must be called before delta_tot_init, or resuming wont work*/
    read_all_nu_state(&d_tot, "testdata/delta_tot_nu.txt");
    /*Then init delta_tot*/
    delta_tot_init(&d_tot, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, transfer,0.33333333);
    /*Check that we will actually do something*/
    assert_true(log(0.3333333)-d_tot.scalefact[d_tot.ia-1] > 1e-5);
    /*So now we have a fully initialised d_tot. Update it!*/
    double delta_nu_curr[d_pow->nbins];
    get_delta_nu_update(&d_tot, 0.33333333, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, delta_nu_curr, transfer);
    /*Check that we did indeed update*/
    assert_true(d_tot.ia == 25);
    /*Check that we get the same answer as the saved powerspectrum*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        assert_true(fabs(delta_nu_curr[ik]/ d_pow->delta_nu_curr[ik] -1) < 1e-2);
    }
    /*Throw out the last stored power spectrum and try again, so that we will perform a save*/
    d_tot.ia--;
    get_delta_nu_update(&d_tot, 0.33333333, d_pow->nbins, d_pow->logkk, d_pow->delta_cdm_curr, delta_nu_curr, transfer);
    assert_true(d_tot.ia == 25);
    /*Check that we get the same answer as the saved powerspectrum*/
    for(int ik=0; ik < d_tot.nk; ik++) {
        /*Be a bit more generous with the error as we have fewer datapoints now*/
        assert_true(fabs(delta_nu_curr[ik]/ d_pow->delta_nu_curr[ik] -1) < 3e-2);
    }
}

/*Load transfer functions from CAMB files.*/
void load_camb_transfer(char * transfer_file, char * matterpow_file, int nk_read, double *delta_cdm, double * delta_nu, double * keffs, double kmin)
{
    int count=0;
    char string[1000];
    /* We aren't interested in modes on scales larger than twice the boxsize*/
    const double scale=1000;
    FILE * fd = fopen(transfer_file, "r");
    if(!fd){
        printf("Could not open '%s' for read.\n", transfer_file);
        exit(1);
    }
    memset(delta_cdm, 0, nk_read*sizeof(double));
    memset(delta_nu, 0, nk_read*sizeof(double));
    while(count < nk_read) {
        /*Note for T_cdm we actually use the CDM+baryons function.*/
        double k, T_cdm,T_nu, dummy, T_tot;
        /* read line from transfer function file */
        char * ret=fgets(string,1000,fd);
        /*End of file*/
        if(!ret)
            break;
        /* Skip comments*/
        if(string[0] == '#')
            continue;
        /* read transfer function file from CAMB */
        if(sscanf(string, " %lg %lg %lg %lg %lg %lg %lg %lg", &k, &dummy, &dummy, &dummy, &dummy, &T_nu, &T_tot, &T_cdm) == 8){
            if(k > kmin){
                /*Combine the massive and massless neutrinos.*/
                /*Set up the total transfer for all the species with particles*/
                delta_cdm[count] = pow(T_cdm/T_tot,2);
                delta_nu[count] = pow(T_nu/T_tot,2);
                count++;
            }
        }
        else
            break;
    }
    fclose(fd);
    assert_true(count == nk_read);

    /*Now read the matter power spectrum*/
    fd = fopen(matterpow_file, "r");
    if(!fd){
        printf("Could not open '%s' for read.\n", matterpow_file);
        exit(1);
    }

    count=0;
    while(count < nk_read) {
        /*Note for T_cdm we actually use the CDM+baryons function.*/
        double k, Pk;
        /* read line from transfer function file */
        char * ret=fgets(string,1000,fd);
        /*End of file*/
        if(!ret)
            break;
        /* Skip comments*/
        if(string[0] == '#')
            continue;
        /* read transfer function file from CAMB */
        if(sscanf(string, "%lg %lg", &k, &Pk) == 2){
            if(k > kmin){
                /*Combine the massive and massless neutrinos.*/
                /*Set up the total transfer for all the species with particles*/
                delta_cdm[count] *= Pk * pow(scale, 3);
                delta_nu[count] *= Pk * pow(scale, 3);
                delta_cdm[count] = sqrt(delta_cdm[count]);
                delta_nu[count] = sqrt(delta_nu[count]);
                keffs[count] = k/scale;
                count++;
            }
        }
        else
            break;
    }
    fclose(fd);
    assert_true(count == nk_read);

    return;
}

#define NREAD 200

/*Check that we reproduce the results of linear theory, as given by CAMB, for neutrinos.*/
static void test_reproduce_linear(void **state)
{
    test_state * ts = (test_state *) *state;
    _omega_nu * omnu = (_omega_nu *) ts->omnu;
    _transfer_init_table transfer;
    const double UnitLength_in_cm = 3.085678e21;
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    allocate_transfer_init_table(&transfer, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, get_omega_nu(omnu, 1), 0.2793, "camb_linear/ics_transfer_0.01.dat");
    /* We will build state by loading the CAMB CDM transfer function at different redshifts,
     * and repeatedly using get_delta_nu_update to advance the internal state of the neutrino code.
     * Then we will compare the CAMB neutrino transfer function to the delta_nu from the neutrino code.*/
    double keffs[NREAD];
    double delta_nu_camb[NREAD];
    double delta_nu[NREAD];
    double delta_cdm[NREAD];
    /*Get the first transfer file*/
    load_camb_transfer("camb_linear/ics_transfer_0.01.dat", "camb_linear/ics_matterpow_0.01.dat", NREAD, delta_cdm, delta_nu, keffs, 2*M_PI/512.);
    _delta_tot_table d_tot;
    allocate_delta_tot_table(&d_tot, NREAD, 0.01, 1, 0.2793, omnu, UnitTime_in_s, UnitLength_in_cm, 0);
    /*Initialise*/
    delta_tot_init(&d_tot, NREAD, keffs, delta_cdm, &transfer,0.01);
    /* Desired accuracy. The first few integrations are less accurate.
     * For the first values we have fewer integration points,
     * and later we assume non-relativistic neutrinos.
     * This is not very important. */
    /* There is also a specific range around k=0.64 where it is slightly less accurate than 1%.
     * This is probably CAMB's fault; presumably it is switching integration method there.*/
    double acc = 0.05;
    for(int i=0; i< 99; i++) {
        double scalefact = 0.01 + i*0.01;
        char tfile[150], mfile[150];
        if(i == 8)
            acc = 2e-2;
        if(i == 22)
            acc = 1.2e-2;
        snprintf(tfile, 150, "camb_linear/ics_transfer_%2g.dat",scalefact);
        snprintf(mfile, 150, "camb_linear/ics_matterpow_%2g.dat",scalefact);
        load_camb_transfer(tfile, mfile, NREAD, delta_cdm, delta_nu_camb, keffs, 2*M_PI/512.);
        get_delta_nu_update(&d_tot, scalefact, NREAD, keffs, delta_cdm, delta_nu, &transfer);
        if(d_tot.ia != i+1)
            printf("%d %d \n", d_tot.ia, i+1);
        assert_true(d_tot.ia == i+1);
        for(int k = 0; k < NREAD; k++) {
//             if(fabs(delta_nu_camb[k] - delta_nu[k]) > 1.2e-2*delta_nu[k])
//                  printf("i = %d k=%g : dnu = %g %g diff %g\n", i, 1000*keffs[k], delta_nu[k], delta_nu_camb[k], fabs(delta_nu_camb[k]/ delta_nu[k]-1));
            assert_true(fabs(delta_nu_camb[k] - delta_nu[k]) < acc*delta_nu[k]);
        }
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
        if(!delta_nu_curr || !delta_cdm_curr || !logkk){
            printf("Failed to allocate memory. nbins=%d\n", nbins);
            free(d_pow);
            free(logkk);
            free(delta_nu_curr);
            free(delta_cdm_curr);
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
        free(d_pow);
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
        free(d_pow);
        free(logkk);
        free(delta_nu_curr);
        free(delta_cdm_curr);
        return 1;
    }
    /*We have actually read the total matter power; so for the CDM power we should subtract off the neutrino power*/
    test_state * ts = (test_state *) malloc(sizeof(test_state));
    ts->omnu = malloc(sizeof(_omega_nu));
    ts->transfer = malloc(sizeof(_transfer_init_table));
    ts->d_pow = d_pow;
    const double MNu[3] = {0.15, 0.15, 0.15};
    init_omega_nu(ts->omnu, MNu, 0.01, 0.7,T_CMB0);
    const double OmegaNua3=get_omega_nu(ts->omnu, 0.01)*pow(0.01,3);
    const double OmegaMa = 0.2793 - get_omega_nu(ts->omnu, 1) + OmegaNua3;
    const double fnu = OmegaNua3/OmegaMa;
    for (int ik = 0; ik < nbins; ik++){
         delta_cdm_curr[ik] = (delta_cdm_curr[ik] - fnu*delta_nu_curr[ik])/(1.-fnu);
    }
    /*Now initialise data structure*/
    init_delta_pow(d_pow, logkk, delta_nu_curr, delta_cdm_curr, nbins,1.);
    const double UnitLength_in_cm = 3.085678e21;
    allocate_transfer_init_table(ts->transfer, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, get_omega_nu(ts->omnu, 1), 0.2793, "testdata/ics_transfer_99.dat");
    const double UnitTime_in_s = UnitLength_in_cm / 1e5;
    /*Set up the global variables for the hubble function before we do anything else!*/
    init_hubble_function(ts->omnu, 0.2793, UnitTime_in_s);
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
        cmocka_unit_test(test_reproduce_linear),
    };
    return cmocka_run_group_tests(tests, setup_delta_pow, teardown_delta_pow);
}
