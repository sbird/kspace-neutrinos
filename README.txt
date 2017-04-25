==About this code==
This code is an extension for Gadget-3 to allow the cheap inclusion of massive neutrinos.
It implements the linear response method detailed in Ali-Haimoud & Bird 2013 (arxiv:1209.0461)
and the hybrid neutrino method in Bird & Ali-Haimoud 2017 (in prep). It has a full test suite, 
which can be run with 'make test' and aims to be easily portable to a variety of different codes.
You are welcome to use it as desired, but please cite Ali-Haimoud & Bird 2013, and 
Bird & Ali-Haimoud 2017 if you use the hybrid neutrinos.

==Accuracy limits==

The equations implemented here are accurate for neutrinos
which are non-relativistic by z=10 (M_nu > 0.005 eV). They are
also accurate for lighter neutrinos within the horizon, because
the clustering there is roughly zero. However, they are *not*
accurate for relativistic neutrinos on superhorizon scales.
We also do not include radiation (photon) perturbations.
If you are interested in this limit you should use the output
of a Boltzmann code, such as CAMB, directly.

=Massless neutrinos=
Note that disabling the integrator is not exactly
equivalent to enabling it with MNu = 0. When
the integrator is enabled, OmegaNu is included in Omega0, the matter density.
When it is disabled, OmegaNu is included in OmegaR, the radiation density.
For best accuracy, disable the neutrino integrator when neutrinos are massless.

==Building==
'make all' will make a static library
'make test' will perform the runtime tests.
'make doc' will run doxygen (if installed) and generate html documentation,
which can be browsed by opening doc/html/index.html in a web browser.

==Dependencies==

GSL
The FFT library used in your PM routine.
cmocka for the test suite.

==Parameters of the code==

The neutrino linear response code takes the following parameters, which should 
be set in the parameter file of your N-body code. Parameter file name shows the default key 
as written in the parameter file, while internal variable gives the member of the 
kspace_params structure in interface_common.c . Default values are given where relevant:
Parameter file name         Internal Variable       Default     Description

STRINGS:
KspaceTransferFunction      KspaceTransferFunction    -         File containing CAMB formatted output transfer functions.
                                                                Used to set initial conditions for the neutrino integration.
FLOATS:
TimeTransfer                TimeTransfer              -         Scale factor from which the neutrino integration should start.
                                                                Must be equal to the scale factor of the simulation initial conditions, and should
                                                                generally be the scale factor at which the CAMB transfer functions were generated.
InputSpectrum_UnitLength_in_cm   ""                3.085678e24  Units of the CAMB transfer function in cm. By default 1 Mpc.
MNue                        MNu[0]                    -         Mass of the lightest neutrino in eV.
MNum                        MNu[1]                    -         Second neutrino mass in eV.
MNut                        MNu[2]                    -         Third neutrino mass. Note the observed mass splitting is not enforced.
Vcrit                       vcrit                    500        Critical velocity in the Fermi-Dirac distribution below which the neutrinos
                                                                are followed with particles, if hybrid neutrinos are on.
NuPartTime                  nu_crit_time              0.3333    Scale factor at which to 'turn on', ie, make active gravitators, 
                                                                the particle neutrinos, if hybrid neutrinos are on.
INTS:
HybridNeutrinosOn           hybrid_neutrinos_on       0         Whether hybrid neutrinos are enabled.

Note that total_powerspectrum returns a power spectrum which is in units of the box, and unnormalised, 
that is, P(k) * N^2, where N is the number of modes in each bin. After investigation, no attempt 
is made to smooth the power spectrum by averaging neighbouring bins.

==Output Files==

Output is saved with every snapshot, and with a restart. 
The main output is the neutrino power spectrum, which is saved to: 
$(All.OutpurDir)/powerspectrum_nu_$(SnapNum).txt
The format of this file is: ( k, P_nu(k) ).
Units are: 1/L, L^3, where L is Gadget internal length units.
A short python script for reading it is found in plot_nu_power.py

The code's internal state is saved to delta_tot_nu.txt.
This contains a table containing the total matter power spectrum, 
delta_tot, as a function of redshift.
Each row is formatted as : "# log(a) delta_tot(k)"
The wavenumbers are not stored, and so this file is not portable to other simulations.
I may change this format in future.

==Using kspace neutrinos with your version of Gadget.==

===MP-Gadget===
The easiest way to use the code is with the public MP-Gadget,
into which it is natively included, and which has many other features.
To enable, set the parameters:
RadiationOn = 1
MassiveNuLinRespOn = 1

The first enables radiation in the background, and the second enables
the neutrino integrator.

When using MP-Gadget, the neutrino state is saved inside the snapshot,
not in delta_tot_nu.txt.

===Gadget-2===
As a convenience, we include patches to the
public version of Gadget-2, which add the required function calls.
To use kspace-neutrinos with Gadget-2:
1) Extract this code to the "Gadget-2.0.7/Gadget2/kspace-neutrinos" subdirectory.
2) Execute the script Gadget-2.0.7/Gadget2/kspace-neutrinos/gadget-2/apply-patches"
which will patch the copy of Gadget-2.0.7 within which it finds itself.
3) Add the two Makefile options:
OPT   += -DINCLUDE_RADIATION
OPT   += -DKSPACE_NEUTRINOS_2
to your Makefile to enable kspace neutrinos.
The first enables radiation density in the background Hubble expansion,
the second enables massive neutrinos.

NOTE When comparing to massless neutrino simulations, you should
enable INCLUDE_RADIATION but not KSPACE_NEUTRINOS_2.

For convenience our patches also add code to output the total
matter powerspectrum (in the same units as powerspec_nu_***.txt)
on every Gadget-2 snapshot.

===Gadget-3===
We internally maintain patches to Gadget-3 which incorporate the neutrino code.
Since Gadget-3 is not public at this time, and many different versions exist,
these patches often require some work to apply. If you are interested in using
them, we encourage you to contact Simeon Bird (spb@ias.edu) directly.

==Porting the neutrino library to a new code==

This version of the neutrino integration library is written
to be easy to port to different versions of Gadget, or other 
non-Gadget-based N-body codes. However, due to the wide 
variety of codes in use, some manual adjustment may still be necessary.
This document details the steps to take.

All interfacing between the neutrino integrator and the rest of 
the code goes through the routines in interface_common.c 
and interface_common.h Ideally therefore, you should just 
call these routines at the correct points in your code 
and add the .c files to your Makefile. 
Interfaces specific to Gadget-3, in particular assuming FFTW2 and Gadget's 
parameter reading routines, are found in interface_gadget.[ch]. 
If you are using Gadget-3, you can just include these files.
Note we do not use global Gadget variables, nor Gadget configuration switches.

The main routines are:
0. InitOmegaNu(a_start, h0, tcmb): Initialises the table for OmegaNu. Should be called before all other functions.
1. OmegaNu(a): the matter density in neutrinos, should be added to the Hubble function
2. allocate_kspace_memory(): allocates and sets up the neutrino module. Do it before calling OmegaNu.

3. add_nu_power_to_rhogrid(): call this inside your PM routine to add the neutrino power to the grid,
Further documentation is provided inside interface_gadget.h
4. save_nu_state(): Saves the internal state of the neutrino integrator to disc, so that resuming from a snapshot works.
5. save_nu_power(): Call this to save the neutrino power spectrum whenever you make a snapshot, or otherwise save the DM power.

Note that add_nu_power_to_rhogrid assumes the (slab-decomposed) FFTW 2, with a type complex number type fftw_complex,
as this is used in almost all gadget versions. If this does not match your code, the routine 
compute_neutrino_power_from_cdm takes a pre-computed matter power spectrum, and you should adapt the for loop
in add_nu_power_to_rhogrid to your own FFT routines.

The .c files which need to be compiled in are:
delta_pow.c - GSL interpolation for neutrino power spectra
delta_tot_table.c - Core integrator that computes delta_nu given a matter power spectrum.
interface_common.c - Routines to do the messy business of interfacing with Gadget. 
                       Also stores the state for the neutrino code in the form of global variables.
interface_gadget.c - Interface routines which assume FFTW2 and are only suitable for Gadget-3.
                       Also stores the state for the neutrino code in the form of global variables.
powerspectrum.c - Routine to compute the power spectrum of a Fourier-transformed density field, 
                    divided up between processors as by FFTW2.
omega_nu_single.c -  Routines to compute OmegaNu and OmegaR 
transfer_init.c - Routine to read and parse CAMB formatter transfer functions.

Other c files are: 
*_test.c - cmocka tests for each module.
gadget_defines.c - support infrastructure normally in gadget for the tests.

You should also provide a routine to read the parameters required 
by the neutrino integrator from your code's parameter file. 
An example routine for Gadget-3 called set_kspace_vars is provided in
interface_gadget.h

Your code should provide the routines declared in gadget_defines.h.
These are:
mymalloc, myfree - free and allocate memory.
hubble_function - H(a)
endrun - ends the simulation.
message - Prints a message to the screen.
In case your code does not provide them, there are simple examples 
defined in gadget_defines.c, primarily for the test suite.
