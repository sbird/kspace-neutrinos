/*This file contains routines for computing OmegaNu, in the full relativistic approximation.
This is needed for adding massive neutrinos to any simulation.

double OmegaNu(double) is the only externally accessible function.
double OmegaNu_single(double a,double mnu, int sp) can be called from the main integrator.
*/

#include "kspace_neutrinos_func.h"
#include "kspace_neutrino_const.h"
#include "kspace_neutrinos_vars.h"
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

/*Note q carries units of eV/c. kT/c has units of eV/c.
 * M_nu has units of eV  Here c=1. */
double rho_nu_int(double q, void * params)
{
        double amnu = *((double *)params);
        double epsilon = sqrt(q*q+amnu*amnu);
        double f0 = 1./(exp(q/(BOLEVK*TNU))+1);
        return q*q*epsilon*f0;
}

/* Tables for rho_nu: stores precomputed values between
 * simulation start and a M_nu = 20 kT_nu*/
#define NRHOTAB 500
double RhoNuLogA[NUSPECIES][NRHOTAB];
double RhoNuTab[NUSPECIES][NRHOTAB];
gsl_interp * RhoNuTab_interp[NUSPECIES];
gsl_interp_accel * RhoNuTab_acc[NUSPECIES];
/*Get the conversion factor to go from (eV/c)^4 to g/cm^3
 * for a **single** neutrino species. */
double get_rho_nu_conversion()
{
        /*q has units of eV/c, so rho_nu now has units of (eV/c)^4*/
        double convert=4*M_PI*2; /* The factor of two is for antineutrinos*/
        /*rho_nu_val now has units of eV^4*/
        /*To get units of density, divide by (c*hbar)**3 in eV s and cm/s */
        const double chbar=1./(2*M_PI*C*HBAR);
        convert*=(chbar*chbar*chbar);
        /*Now has units of (eV)/(cm^3)*/
        /* 1 eV = 1.60217646 × 10-12 g cm^2 s^(-2) */
        /* So 1eV/c^2 = 1.7826909604927859e-33 g*/
        /*So this is rho_nu_val in g /cm^3*/
        convert*=(1.60217646e-12/C/C);
        return convert;
}

/*Seed a pre-computed table of rho_nu values for speed*/
void tabulate_rho_nu(double af,double mnu,int sp)
{
     int i;
     double abserr;
     const double logA0=log(A0);
     const double logaf=log(af);
     gsl_function F;
     gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
     F.function = &rho_nu_int;
     for(i=0; i< NRHOTAB; i++){
        double param;
        RhoNuLogA[sp][i]=logA0+i*(logaf-logA0)/(NRHOTAB-1);
        param=mnu*exp(RhoNuLogA[sp][i]);
        F.params = &param;
        gsl_integration_qag (&F, 0, 500*BOLEVK*TNU,0 , 1e-9,GSL_VAL,6,w,&(RhoNuTab[sp][i]), &abserr);
        RhoNuTab[sp][i]=RhoNuTab[sp][i]/pow(exp(RhoNuLogA[sp][i]),4)*get_rho_nu_conversion();
     }
     gsl_integration_workspace_free (w);
     RhoNuTab_acc[sp] = gsl_interp_accel_alloc();
     RhoNuTab_interp[sp]=gsl_interp_alloc(gsl_interp_cspline,NRHOTAB);
     if(!RhoNuTab_interp[sp] || !RhoNuTab_acc[sp] || gsl_interp_init(RhoNuTab_interp[sp],RhoNuLogA[sp],RhoNuTab[sp],NRHOTAB))
         terminate("Could not initialise tables for RhoNu\n");
     return;
}

/* Value of kT/aM_nu on which to switch from the
 * analytic expansion to the numerical integration*/
#define NU_SW 50
//1.878 82(24) x 10-29 h02 g/cm3 = 1.053 94(13) x 104 h02 eV/cm3
/*Finds the physical density in neutrinos for a single neutrino species*/
double rho_nu(double a,double mnu,int sp)
{
        double rho_nu_val;
        double amnu=a*mnu;
        const double kT=BOLEVK*TNU;
        const double kTamnu2=(kT*kT/amnu/amnu);
        /*Do it analytically if we are in a regime where we can
         * The next term is 141682 (kT/amnu)^8.
         * At kT/amnu = 8, higher terms are larger and the series stops converging.
         * Don't go lower than 50 here. */
        if(NU_SW*NU_SW*kTamnu2 < 1){
            /*Heavily non-relativistic*/
            /*The constants are Riemann zetas: 3,5,7 respectively*/
            rho_nu_val=mnu*pow(kT/a,3)*(1.5*1.202056903159594+kTamnu2*45./4.*1.0369277551433704+2835./32.*kTamnu2*kTamnu2*1.0083492773819229+80325/32.*kTamnu2*kTamnu2*kTamnu2*1.0020083928260826)*get_rho_nu_conversion();
        }
        else if(amnu < 1e-6*kT){
            /*Heavily relativistic: we could be more accurate here,
             * but in practice this will only be called for massless neutrinos, so don't bother.*/
            rho_nu_val=7*pow(M_PI*kT/a,4)/120.*get_rho_nu_conversion();
        }
        else{
            if(!(RhoNuTab_interp[sp]))
                tabulate_rho_nu(NU_SW*kT/mnu,mnu,sp);
            rho_nu_val=gsl_interp_eval(RhoNuTab_interp[sp],RhoNuLogA[sp],RhoNuTab[sp],log(a),RhoNuTab_acc[sp]);
        }
        return rho_nu_val;
}


/* Return the matter density in a single neutrino species.
 * Not externally callable*/
double OmegaNu_single(double a,double mnu, int sp)
{
        double rhonu;
        rhonu=rho_nu(a,mnu,sp);
        rhonu /= (3* HUBBLE* HUBBLE / (8 * M_PI * GRAVITY));
        rhonu /= kspace_vars.HubbleParam*kspace_vars.HubbleParam;
        //Remove neutrino density which is particle based, if necessary
#ifdef HYBRID_NEUTRINOS
        if(! All.slow_neutrinos_analytic){
            if(!nufrac_low)
                nufrac_low = _nufrac_low(mnu);
            rhonu*=(1-nufrac_low);
        }
#endif
        return rhonu;
}

/* Return the total matter density in neutrinos.
 * rho_nu and friends are not externally callable*/
double OmegaNu(double a)
{
        double rhonu=0;
        int mi;
        double OmegaNu_one[NUSPECIES];
        for(mi=0; mi<NUSPECIES; mi++){
             int mmi;
             for(mmi=0; mmi<mi; mmi++){
                 if(fabs(kspace_params.MNu[mi] -kspace_params.MNu[mmi]) < FLOAT)
                   break;
             }
             if(mmi==mi)
                 OmegaNu_one[mi]=OmegaNu_single(a,kspace_params.MNu[mi],mi);
            rhonu+=OmegaNu_one[mmi];
        }
        return rhonu;
}
