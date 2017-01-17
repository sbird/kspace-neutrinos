LFLAGS += -lgsl -lgslcblas -lpthread
CFLAGS +=-O2 -ffast-math -g -Wall -fopenmp ${OPT}
LFLAGS += -lm -lgomp

OBJS = transfer_init.o delta_tot_table.o powerspectrum.o delta_pow.o interface_common.o omega_nu_single.o interface_gadget.o
INCL = kspace_neutrino_const.h interface_common.h interface_gadget.h powerspectrum.h delta_pow.h omega_nu_single.h gadget_defines.h transfer_init.h delta_tot_table.h Makefile

.PHONY : clean all test doc

all: lib

doc:
	doxygen

lib: ${OBJS}
	ar rcs libkspace_neutrinos_2.a $^

test: run_omega_nu_single_test run_transfer_init_test run_powerspectrum_test run_delta_pow_test run_delta_tot_table_test

run_%_test: %_test
	./$^

%.o: %.c ${INCL}
	$(CC) -c $(CFLAGS) $< -o $@

interface_common.o: interface_common.c ${INCL}
	mpicc -c $(CFLAGS) $< -o $@

interface_gadget.o: interface_gadget.c ${INCL}
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
	mpicc $(CFLAGS) $^ -o $@ -lcmocka $(LFLAGS) -lsrfftw -lsfftw

clean:
	rm -f $(OBJS) gadget_defines.o
