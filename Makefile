OPT   +=  -DKSPACE_NEUTRINOS_2  # Enable kspace neutrinos
OPT   += -DHYBRID_NEUTRINOS

LFLAGS += -lgsl -lgslcblas -lpthread
CFLAGS +=-O2 -ffast-math -g -Wall -fopenmp ${OPT}
LFLAGS += -lm -lgomp

OBJS = transfer_init.o delta_tot_table.o powerspectrum.o delta_pow.o kspace_neutrinos_2.o omega_nu_single.o
INCL = kspace_neutrino_const.h kspace_neutrinos_2.h powerspectrum.h delta_pow.h omega_nu_single.h gadget_defines.h transfer_init.h delta_tot_table.h Makefile

.PHONY : clean all test

all: ${OBJS}

test: omega_nu_single_test #transfer_init_test
	./$^

%.o: %.c ${INCL}
	$(CC) -c $(CFLAGS) $< -o $@

%_test: %_test.c %.o gadget_defines.o
	$(CC) $(CFLAGS) $^ -o $@ -lcmocka $(LFLAGS) 

clean:
	rm -f $(OBJS)
