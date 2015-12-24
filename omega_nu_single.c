#include "omega_nu_single.h"

#include "kspace_neutrino_const.h"
#include "kspace_neutrinos_private.h"

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define HBAR    6.582119e-16  /*hbar in units of eV s*/

void init_omega_nu(_omega_nu * omnu, const double MNu[], const double Omega0)
{
    /*Store matter fraction*/
    omnu->Omega0 = Omega0;
    /*First compute which neutrinos are degenerate with each other*/
    for(int mi=0; mi<NUSPECIES; mi++){
        int mmi;
        omnu->nu_degeneracies[mi]=0;
        for(mmi=0; mmi<mi; mmi++){
            if(fabs(MNu[mi] -MNu[mmi]) < FLOAT){
                omnu->nu_degeneracies[mmi]+=1;
                break;
            }
        }
        if(mmi==mi) {
            omnu->nu_degeneracies[mi]=1;
        }
    }
    /*Now allocate a table for the species we want*/
    for(int mi=0; mi<NUSPECIES; mi++){
        if(omnu->nu_degeneracies[mi])
            omnu->RhoNuTab[mi] = (_rho_nu_single *) mymalloc("RhoNuTab", sizeof(_rho_nu_single));
    }
}
/* Return the total matter density in neutrinos.
 * rho_nu and friends are not externally callable*/
double get_omega_nu(_omega_nu *omnu, double a)
{
        double rhonu=0;
        for(int mi=0; mi<NUSPECIES; mi++) {
            if(omnu->nu_degeneracies[mi] > 0)
                 rhonu += omnu->nu_degeneracies[mi] * omega_nu_single(omnu->RhoNuTab[mi], a);
        }
        return rhonu;
}


/* Value of kT/aM_nu on which to switch from the
 * analytic expansion to the numerical integration*/
#define NU_SW 50

/*Note q carries units of eV/c. kT/c has units of eV/c.
 * M_nu has units of eV  Here c=1. */
double rho_nu_int(double q, void * params)
{
        double amnu = *((double *)params);
        double epsilon = sqrt(q*q+amnu*amnu);
        double f0 = 1./(exp(q/(BOLEVK*TNU))+1);
        return q*q*epsilon*f0;
}

/*Get the conversion factor to go from (eV/c)^4 to g/cm^3
 * for a **single** neutrino species. */
double get_rho_nu_conversion()
{
        /*q has units of eV/c, so rho_nu now has units of (eV/c)^4*/
        double convert=4*M_PI*2; /* The factor of two is for antineutrinos*/
        /*rho_nu_val now has units of eV^4*/
        /*To get units of density, divide by (c*hbar)**3 in eV s and cm/s */
        const double chbar=1./(2*M_PI*LIGHTCGS*HBAR);
        convert*=(chbar*chbar*chbar);
        /*Now has units of (eV)/(cm^3)*/
        /* 1 eV = 1.60217646 × 10-12 g cm^2 s^(-2) */
        /* So 1eV/c^2 = 1.7826909604927859e-33 g*/
        /*So this is rho_nu_val in g /cm^3*/
        convert*=(1.60217646e-12/LIGHTCGS/LIGHTCGS);
        return convert;
}

/*Seed a pre-computed table of rho_nu values for speed*/
void rho_nu_init(_rho_nu_single * rho_nu_tab, double a0, const double mnu, double HubbleParam)
{
     int i;
     double abserr;
     const double logA0=log(a0);
     const double logaf=log(NU_SW*BOLEVK*TNU/mnu);
     gsl_function F;
     gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
     F.function = &rho_nu_int;
     /*Initialise constants*/
     rho_nu_tab->mnu = mnu;
     rho_nu_tab->omega_prefac = (3* HUBBLE* HUBBLE / (8 * M_PI * GRAVITY))*HubbleParam*HubbleParam;
     for(i=0; i< NRHOTAB; i++){
        double param;
        rho_nu_tab->loga[i]=logA0+i*(logaf-logA0)/(NRHOTAB-1);
        param=mnu*exp(rho_nu_tab->loga[i]);
        F.params = &param;
        gsl_integration_qag (&F, 0, 500*BOLEVK*TNU,0 , 1e-9,GSL_VAL,6,w,&(rho_nu_tab->rhonu[i]), &abserr);
        rho_nu_tab->rhonu[i]=rho_nu_tab->rhonu[i]/pow(exp(rho_nu_tab->loga[i]),4)*get_rho_nu_conversion();
     }
     gsl_integration_workspace_free (w);
     rho_nu_tab->acc = gsl_interp_accel_alloc();
     rho_nu_tab->interp=gsl_interp_alloc(gsl_interp_cspline,NRHOTAB);
     if(!rho_nu_tab->interp || !rho_nu_tab->acc || gsl_interp_init(rho_nu_tab->interp,rho_nu_tab->loga,rho_nu_tab->rhonu,NRHOTAB))
         terminate("Could not initialise tables for RhoNu\n");
     return;
}


//1.878 82(24) x 10-29 h02 g/cm3 = 1.053 94(13) x 104 h02 eV/cm3
/*Finds the physical density in neutrinos for a single neutrino species*/
double rho_nu(_rho_nu_single * rho_nu_tab, double a)
{
        double rho_nu_val;
        double amnu=a*rho_nu_tab->mnu;
        const double kT=BOLEVK*TNU;
        const double kTamnu2=(kT*kT/amnu/amnu);
        /*Do it analytically if we are in a regime where we can
         * The next term is 141682 (kT/amnu)^8.
         * At kT/amnu = 8, higher terms are larger and the series stops converging.
         * Don't go lower than 50 here. */
        if(NU_SW*NU_SW*kTamnu2 < 1){
            /*Heavily non-relativistic*/
            /*The constants are Riemann zetas: 3,5,7 respectively*/
            rho_nu_val=rho_nu_tab->mnu*pow(kT/a,3)*(1.5*1.202056903159594+kTamnu2*45./4.*1.0369277551433704+2835./32.*kTamnu2*kTamnu2*1.0083492773819229+80325/32.*kTamnu2*kTamnu2*kTamnu2*1.0020083928260826)*get_rho_nu_conversion();
        }
        else if(amnu < 1e-6*kT){
            /*Heavily relativistic: we could be more accurate here,
             * but in practice this will only be called for massless neutrinos, so don't bother.*/
            rho_nu_val=7*pow(M_PI*kT/a,4)/120.*get_rho_nu_conversion();
        }
        else{
            rho_nu_val=gsl_interp_eval(rho_nu_tab->interp,rho_nu_tab->loga,rho_nu_tab->rhonu,log(a),rho_nu_tab->acc);
        }
        return rho_nu_val;
}

/* Return the matter density in a single neutrino species.
 * Not externally callable*/
double omega_nu_single(_rho_nu_single * rho_nu_tab, double a)
{
        double rhonu=rho_nu(rho_nu_tab, a);
        rhonu /= rho_nu_tab->omega_prefac; 
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
