#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdio.h>
#include <math.h>
#include "omega_nu_single.h"
#include "gadget_defines.h"

/* A test case that checks initialisation. */
static void test_rho_nu_init(void **state) {
    (void) state; /* unused */
    double mnu = 0.06;
    _rho_nu_single rho_nu_tab;
    /*Initialise*/
    rho_nu_init(&rho_nu_tab, 0.01, mnu, 0.7);
    /*Check everything initialised ok*/
    assert_true(rho_nu_tab.mnu == mnu);
    assert_true(rho_nu_tab.acc);
    assert_true(rho_nu_tab.interp);
    assert_true(rho_nu_tab.loga);
    assert_true(rho_nu_tab.rhonu);
    /*Check that loga is correctly ordered (or interpolation won't work)*/
    for(int i=1; i<NRHOTAB; i++){
        assert_true(rho_nu_tab.loga[i] > rho_nu_tab.loga[i-1]);
    }
}

/* Check that the table gives the right answer. */
static void test_omega_nu_single(void **state) {
    (void) state; /* unused */
    double mnu = 0.5;
    double hubble = 0.7;
    _rho_nu_single rho_nu_tab;
    /*Initialise*/
    rho_nu_init(&rho_nu_tab, 0.01, mnu, hubble);
    assert_true(rho_nu_tab.mnu == mnu);
    /*This is the critical density at z=0:
     * we allow it to change by what is (at the time of writing) the uncertainty on G.*/
    assert_true(fabs(rho_nu_tab.rhocrit - 1.8784e-29*hubble*hubble) < 5e-5*rho_nu_tab.rhocrit);
    /*Check everything initialised ok*/
    double omnuz0 = omega_nu_single(&rho_nu_tab, 1);
    /*Check redshift scaling*/
    assert_true(fabs(omnuz0/pow(0.5,3) - omega_nu_single(&rho_nu_tab, 0.5)) < 5e-5*omnuz0);
    /*Check not just a^-3 scaling*/
    assert_true(omnuz0/pow(0.01,3) <  omega_nu_single(&rho_nu_tab, 0.01));
    /*This fails...is it the right value?*/
    printf("val: %g %g\n", omnuz0, mnu/93.14/hubble/hubble/0.999787);
    assert_true(fabs(omnuz0 - mnu/93.14/hubble/hubble) < 1e-6*omnuz0);
}

/*Check massless neutrinos work*/
#define STEFAN_BOLTZMANN 5.670373e-5
#define OMEGAR (4*STEFAN_BOLTZMANN*8*M_PI*GRAVITY/(3*LIGHTCGS*LIGHTCGS*LIGHTCGS*HUBBLE*HUBBLE*HubbleParam*HubbleParam)*pow(T_CMB0,4))
static void test_massless_omega_nu_single(void **state) {
    _rho_nu_single rho_nu_tab_nomass;
    /*Initialise*/
    double HubbleParam = 0.7;
    rho_nu_init(&rho_nu_tab_nomass, 0.01, 0., HubbleParam);
    /*Check gives the right value*/
    double omnunomassz0 = omega_nu_single(&rho_nu_tab_nomass, 1);
    assert_true(omnunomassz0 - OMEGAR*7./8.*pow(pow(4/11.,1/3.)*1.00381,4)< 1e-5*omnunomassz0);
    assert_true(omnunomassz0/pow(0.5,4) == omega_nu_single(&rho_nu_tab_nomass, 0.5));
}

static void test_omega_nu_init_degenerate(void **state) {
    /*Check we correctly initialise omega_nu with degenerate neutrinos*/
    _omega_nu omnu;
    /*Initialise*/
    double MNu[3] = {0.2,0.2,0.2};
    init_omega_nu(&omnu, MNu, 0.3, 0.01, 0.7);
    /*Check that we initialised the right number of arrays*/
    assert_int_equal(omnu.nu_degeneracies[0], 3);
    assert_int_equal(omnu.nu_degeneracies[1], 0);
    assert_true(omnu.RhoNuTab[0]);
    assert_false(omnu.RhoNuTab[1]);
    assert_true(omnu.Omega0 == 0.3);
}

static void test_omega_nu_init_nondeg(void **state) {
    /*Now check that it works with a non-degenerate set*/
    _omega_nu omnu;
    /*Initialise*/
    double MNu[3] = {0.2,0.1,0.3};
    init_omega_nu(&omnu, MNu, 0.3, 0.01, 0.7);
    /*Check that we initialised the right number of arrays*/
    for(int i=0; i<3; i++) {
        assert_int_equal(omnu.nu_degeneracies[i], 1);
        assert_true(omnu.RhoNuTab[i]);
    }
}

static void test_get_omega_nu(void **state) {
    /*Check that we get the right answer from get_omega_nu, in terms of rho_nu*/
    _omega_nu omnu;
    /*Initialise*/
    double MNu[3] = {0.2,0.1,0.3};
    init_omega_nu(&omnu, MNu, 0.3, 0.01, 0.7);
    double total =0;
    for(int i=0; i<3; i++) {
        total += omega_nu_single(omnu.RhoNuTab[i], 0.5);
    }
    assert_true(get_omega_nu(&omnu, 0.5) - total < 1e-6*total);
}

int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_rho_nu_init),
        cmocka_unit_test(test_omega_nu_single),
        cmocka_unit_test(test_massless_omega_nu_single),
        cmocka_unit_test(test_omega_nu_init_degenerate),
        cmocka_unit_test(test_omega_nu_init_nondeg),
        cmocka_unit_test(test_get_omega_nu),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
