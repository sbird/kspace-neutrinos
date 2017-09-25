"""Routine to plot the neutrino power spectrum, as output by this code."""

import math
import numpy as np
import scipy.interpolate

def load_genpk(path,box):
    """Load a GenPk format power spectum, plotting the DM and the neutrinos (if present)
    Does not plot baryons."""
    #Load DM P(k)
    matpow=np.loadtxt(path)
    scale=2*math.pi/box
    #Adjust Fourier convention to match CAMB.
    simk=matpow[1:,0]*scale
    Pk=matpow[1:,1]/scale**3*(2*math.pi)**3
    return (simk,Pk)

def get_nu_power(filename):
    """Reads the neutrino power spectrum.
    Format is: ( k, P_nu(k) ).
    Units are: 1/L, L^3, where L is
    Gadget internal length units for
    Gadget-2 and Mpc/h for MP-Gadget."""
    data = np.loadtxt(filename)
    k = data[:,0]
    #Convert fourier convention to CAMB.
    pnu = data[:,1]
    return (k, pnu)

def get_camb_nu_power(matpow, transfer):
    """Plot the neutrino power spectrum from CAMB.
    This is just the matter power multiplied
    by the neutrino transfer function.
    CAMB internal units are used.
    Assume they have the same k binning."""
    matter = np.loadtxt(matpow)
    trans = np.loadtxt(transfer)
    #Adjust Fourier convention to match CAMB.
    tnufac = (trans[:,5]/trans[:,6])**2
    return matter[:,0], matter[:,1]*tnufac

def get_hyb_nu_power(nu_filename, genpk_neutrino, box, npart, part_prop):
    """Get the total matter power spectrum when some of it is in particles, some analytic."""
    (k_part,pk_part)=load_genpk(genpk_neutrino,box)
    (k_sl, pk_sl) = get_nu_power(nu_filename)
    intp=scipy.interpolate.InterpolatedUnivariateSpline(np.log(k_part),pk_part)
    pk_part_r = intp(np.log(k_sl))
    shot=(box/npart)**3/(2*math.pi**2)*np.ones(np.size(pk_part_r))
    pk = (part_prop*np.sqrt(pk_part_r-shot)+(1-part_prop)*np.sqrt(pk_sl))**2
    return (k_sl, pk)
