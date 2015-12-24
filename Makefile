OPT   +=  -DKSPACE_NEUTRINOS_2  # Enable kspace neutrinos
OPT   += -DHYBRID_NEUTRINOS

LFLAGS += -lgsl -lgslcblas -lpthread
CFLAGS +=-O2 -ffast-math -g -c -Wall -fopenmp ${OPT}
LFLAGS += -lm -lgomp

OBJS = transfer_init.o delta_tot_table.o powerspectrum.o delta_pow.o kspace_neutrinos_add_pm_functions.o omega_nu_single.o
INCL = kspace_neutrino_const.h kspace_neutrinos_func.h powerspectrum.h delta_pow.h omega_nu_single.h kspace_neutrinos_private.h transfer_init.h delta_tot_table.h Makefile

.PHONY : clean all test

all: ${OBJS}

test: btest
	./$^

%.o: %.c ${INCL}
	$(CC) $(CFLAGS) $< -o $@

test.o: test.cpp ${INCL}
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@

btest: test.o ${OBJS}
	${LINK} ${LFLAGS} -lboost_unit_test_framework $^ -o  $@
clean:
	rm -f $(OBJS)
