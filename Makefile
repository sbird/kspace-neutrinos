OPT   +=  -DKSPACE_NEUTRINOS_2  # Enable kspace neutrinos
OPT   += -DHYBRID_NEUTRINOS

LFLAGS += -lgsl -lgslcblas -lpthread
CFLAGS +=-O2 -ffast-math -g -Wall -fopenmp ${OPT}
LFLAGS += -lm -lgomp

OBJS = transfer_init.o delta_tot_table.o powerspectrum.o delta_pow.o kspace_neutrinos_2.o omega_nu_single.o
INCL = kspace_neutrino_const.h kspace_neutrinos_2.h powerspectrum.h delta_pow.h omega_nu_single.h gadget_defines.h transfer_init.h delta_tot_table.h Makefile

.PHONY : clean all test

all: lib

lib: ${OBJS}
	ar rcs libkspace_neutrinos_2.a $^

test: run_omega_nu_single_test run_transfer_init_test run_powerspectrum_test run_delta_pow_test run_delta_tot_table_test

run_%_test: %_test
	./$^

%.o: %.c ${INCL}
	$(CC) -c $(CFLAGS) $< -o $@

kspace_neutrinos_2.o: kspace_neutrinos_2.c ${INCL}
	mpicc -c $(CFLAGS) $< -o $@

powerspectrum.o: powerspectrum.c ${INCL}
	mpicc -c $(CFLAGS) $< -o $@

%_test: %_test.c %.o omega_nu_single.o gadget_defines.o
	$(CC) $(CFLAGS) $^ -o $@ -lcmocka $(LFLAGS)

delta_tot_table_test: delta_tot_table_test.c delta_tot_table.o delta_pow.o transfer_init.o omega_nu_single.o gadget_defines.o
	$(CC) $(CFLAGS) $^ -o $@ -lcmocka $(LFLAGS)

#This needs MPI
#The fftw link must match the include in powerspectrum_test.c
powerspectrum_test: powerspectrum_test.c powerspectrum.o omega_nu_single.o gadget_defines.o
	mpicc $(CFLAGS) $^ -o $@ -lcmocka $(LFLAGS) -lsfftw -lsrfftw

clean:
	rm -f $(OBJS) gadget_defines.o
