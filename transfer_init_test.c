#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdio.h>
#include <math.h>
#include "transfer_init.h"
#include "gadget_defines.h"


// void allocate_transfer_init_table(_transfer_init_table *t_init, int nk_in, const double BoxSize, const double UnitLength_in_cm, const double InputSpectrum_UnitLength_in_cm, const double OmegaBaryonCAMB, const char * KspaceTransferFunction, _omega_nu * omnu);

static void test_transfer_init(void **state) {
    (void) state;
    /*First we need the background density*/
    _omega_nu omnu;
    double MNu[3] = {0.15,0.15,0.15};
    init_omega_nu(&omnu, MNu, 0.2793, 0.01, 0.7);
    /*Check we can correctly read a transfer table from disc.*/
    _transfer_init_table transfer;
    /*Initialise the table. Use kpc units and a very large boxsize, so we get the whole table.*/
    const double UnitLength_in_cm = 3.085678e21;
    allocate_transfer_init_table(&transfer, 500, 1000000000, UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", &omnu);
    /*Check that we read all the rows*/
    assert_true(transfer.NPowerTable == 336);
    free_transfer_init_table(&transfer);
    /*Check the boxsize bit is working*/
    allocate_transfer_init_table(&transfer, 500, 512000, UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", &omnu);
    assert_true(transfer.logk[0] > log(M_PI/512000));
    assert_true(transfer.NPowerTable == 271);
    /*Check we are scaling correctly with unit system*/
    _transfer_init_table transfer2;
    allocate_transfer_init_table(&transfer2, 500, 512, 1e3*UnitLength_in_cm, UnitLength_in_cm*1e3, 0.0463, "testdata/ics_transfer_99.dat", &omnu);
    assert_true(transfer2.logk[0] == transfer.logk[0]+log(1e3));
    assert_true(transfer2.T_nu[0] == transfer.T_nu[0]);
    assert_true(transfer2.NPowerTable == transfer.NPowerTable);
    /*Check the read value has not changed: these are just the values read the first time I wrote the test*/
    assert_true(fabs(transfer.T_nu[0] - 0.508478) < 1e-6);
    assert_true(fabs(transfer.T_nu[30] - 0.0122563) < 1e-6);
    free_transfer_init_table(&transfer);
    free_transfer_init_table(&transfer2);
}


int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_transfer_init),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
