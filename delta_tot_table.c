/*This file contains calculations for the Fourier-space semi-linear neutrino method
 * described in Ali-Haimoud and Bird 2012.
 * delta_tot_table stores the state of the integrator, which includes the matter power spectrum over all past time.
 * This file contains routines for manipulating this structure; updating it by computing a new neutrino power spectrum,
 * from the non-linear CDM power, and saving and loading the structure to and from disc.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_bessel.h>

#include "delta_tot_table.h"
#include "kspace_neutrinos_private.h"
#include "kspace_neutrinos_vars.h"
#include "kspace_neutrinos_func.h"
#include "kspace_neutrino_const.h"

/*Allocate memory for delta_tot_table. This is separate from delta_tot_init because we need to allocate memory
 * before we have the information needed to initialise it*/
void allocate_delta_tot_table(_delta_tot_table *d_tot, int nk_in)
{
   /*Memory allocations need to be done on all processors*/
   d_tot->nk=nk_in;
   /* Allocate memory for delta_tot here, so that we can have further memory allocated and freed
    * before delta_tot_init is called. The number nk here should be larger than the actual value needed.*/
   /*Allocate pointers to each k-vector*/
   d_tot->namax=ceil(100*(kspace_vars.TimeMax-kspace_params.TimeTransfer))+2;
   d_tot->ia=0;
   d_tot->delta_tot =(double **) mymalloc("kspace_delta_tot",nk_in*sizeof(double *));
   /*Allocate list of scale factors, and space for delta_tot, in one operation.*/
   d_tot->scalefact = (double *) mymalloc("kspace_scalefact",d_tot->namax*(nk_in+1)*sizeof(double));
   /*Allocate actual data. Note that this means data can be accessed either as:
    * delta_tot[k][a] OR as
    * delta_tot[0][a+k*namax] */
   d_tot->delta_tot[0] = d_tot->scalefact+d_tot->namax;
   for(int count=1; count< nk_in; count++)
        d_tot->delta_tot[count] = d_tot->delta_tot[0] + count*d_tot->namax;
   /*Allocate space for the initial neutrino power spectrum*/
   d_tot->delta_nu_init =(double *) mymalloc("kspace_delta_nu_init",2*nk_in*sizeof(double));
   d_tot->delta_nu_last=d_tot->delta_nu_init+nk_in;
}

void handler (const char * reason, const char * file, int line, int gsl_errno)
{
    printf("GSL_ERROR in file: %s, line %d, errno:%d\n",file, line, gsl_errno);
    terminate(reason);
}

/* Constructor. transfer_init_tabulate must be called before this function.
 * Initialises delta_tot (including from a file) and delta_nu_init from the transfer functions.
 * read_all_nu_state must be called before this if you want reloading from a snapshot to work
 * Note delta_cdm_curr includes baryons*/
void delta_tot_init(_delta_tot_table *d_tot, int nk_in, double wavenum[], double delta_cdm_curr[], _transfer_init_table *t_init, const double Omega0)
{
    if(nk_in > d_tot->nk){
           char err[500];
           sprintf(err,"input power of %d is longer than memory of %d\n",nk_in,d_tot->nk);
           terminate(err);
    }
    gsl_set_error_handler(handler);
    d_tot->nk=nk_in;
    /*Construct delta_nu_init from the transfer functions.*/
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp *spline;
    const double OmegaNua3=OmegaNu(t_init->TimeTransfer)*pow(t_init->TimeTransfer,3);
    const double OmegaNu_today = OmegaNu(1);
    /*Set the prefactor for delta_nu*/
    d_tot->delta_nu_prefac = 1.5 *Omega0 *H100*H100/LIGHT;
    int ik;
    /*Check that if we are restarting from a snapshot, we successfully read a table*/
    if(fabs(kspace_vars.TimeBegin - t_init->TimeTransfer) >1e-4 && (!d_tot->ia))
        terminate("Transfer function not at the same time as simulation start (are you restarting from a snapshot?) and could not read delta_tot table\n");
    /*Initialise the first delta_tot to use the first timestep's delta_cdm_curr
        * so that it includes potential Rayleigh scattering. */
    d_tot->scalefact[0]=log(t_init->TimeTransfer);
    spline=gsl_interp_alloc(gsl_interp_cspline,t_init->NPowerTable);
    gsl_interp_init(spline,t_init->logk,t_init->T_nu,t_init->NPowerTable);
    for(ik=0;ik<d_tot->nk;ik++){
            double T_nubyT_0 = gsl_interp_eval(spline,t_init->logk,t_init->T_nu,log(wavenum[ik]),acc);
            /*The total power spectrum using neutrinos and radiation from the CAMB transfer functions:
            * The CAMB transfer functions are defined such that
            * P_cdm ~ T_cdm^2 (and some other constant factors)
            * then P_t = (Omega_cdm P_cdm + Omega_nu P_nu)/(Omega_cdm + Omega_nu)
            *          = P_cdm (Omega_cdm+ Omega_nu (P_nu/P_cdm)) / (Omega_cdm +Omega_nu)
            *          = P_cdm (Omega_cdm+ Omega_nu (T_nu/T_cdm)^2) / (Omega_cdm+Omega_nu) */
            double CDMtoTot=((Omega0-OmegaNu_today)+pow(T_nubyT_0,2)*OmegaNua3)/(Omega0-OmegaNu_today+OmegaNua3);
            /*We only want to use delta_cdm_curr if we are not restarting*/
            if(!d_tot->ia)
                    d_tot->delta_tot[ik][0] = delta_cdm_curr[ik]*sqrt(CDMtoTot);
            /* Also initialise delta_nu_init here to save time later.
            * Use the first delta_tot, in case we are resuming.*/
            d_tot->delta_nu_init[ik] = d_tot->delta_tot[ik][0]/sqrt(CDMtoTot)*fabs(T_nubyT_0);
    }
    gsl_interp_accel_free(acc);
    gsl_interp_free(spline);
    if(!d_tot->ia)
            d_tot->ia=1;
    if(d_tot->ThisTask==0){
        /*Get a clean debug restart file*/
#ifndef NOCALLSOFSYSTEM
        system("mv delta_tot_nu.txt delta_tot_nu.txt.bak");
#else
        fd=fopen("delta_tot_nu.txt","w");
        fclose(fd);
#endif
        save_all_nu_state(d_tot, NULL);
    }
    /*Initialise delta_nu_last*/
    get_delta_nu_combined(d_tot, exp(d_tot->scalefact[d_tot->ia-1])-2*FLOAT, d_tot->ia, wavenum, d_tot->delta_nu_last);
    d_tot->delta_tot_init_done=1;
    return;
}

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(_delta_tot_table *d_tot, double a, int Na, double wavenum[],  double delta_nu_curr[])
{
    double Omega_nu_tot=OmegaNu(a);
    int mi;
    double delta_nu_single[NUSPECIES][d_tot->nk];
    /*Initialise delta_nu_curr*/
    for(mi=0;mi<d_tot->nk;mi++)
        delta_nu_curr[mi]=0;
    /*Get each neutrinos species and density separately and add them to the total.
     * Neglect perturbations in massless neutrinos.*/
    for(mi=0;mi< NUSPECIES; mi++){
         if(kspace_params.MNu[mi] > 0){
             int ik,mmi;
             double OmegaNu_frac=OmegaNu_single(a,kspace_params.MNu[mi],mi)/Omega_nu_tot;
             for(mmi=0; mmi<mi; mmi++){
                 if(fabs(kspace_params.MNu[mi] -kspace_params.MNu[mmi]) < FLOAT)
                   break;
             }
             if(mmi==mi)
                get_delta_nu(d_tot, a, Na, wavenum, delta_nu_single[mi],kspace_params.MNu[mi]);
             for(ik=0; ik<d_tot->nk; ik++)
                 delta_nu_curr[ik]+=delta_nu_single[mmi][ik]*OmegaNu_frac;
         }
    }
    return;
}

/** Callable function to calculate the power spectra.
 * Calls the rest of the code internally.
 * In reality we will be given P_cdm(current) but not delta_tot.
 * Here is the full function that deals with this
 * Static variables used:
 * nk is the number of k values stored in each power spectrum
 * double scalefact[] is an array of length namax containing log (scale factor) at which the power spectrum is stored
 * double **delta_tot is a pointer to nk arrays of length namax containing the total power spectrum.
 * int ia is the number of already "recorded" time steps, i.e. scalefact[0...ia-1] is recorded.
 * Current time corresponds to index ia (but is only recorded if sufficiently far from previous time).
 * Caution: ia here is different from Na in get_delta_nu (Na = ia+1).
 * @param a is the current scale factor
 * @param logk is an array of length nk containing (natural) log k
 * @param delta_cdm_curr is an array of length nk containing the square root of the current cdm power spectrum
 * @param delta_nu_curr is an array of length nk which stores the square root of the current neutrino power spectrum. Main output of the function.
******************************************************************************************************/
void get_delta_nu_update(_delta_tot_table *d_tot, double a, int nk_in, double logk[], double delta_cdm_curr[], double delta_nu_curr[])
{
  int ik;
  double wavenum[nk_in];
  const double OmegaNua3=OmegaNu(a)*pow(a,3);
  const double Omega0 = kspace_vars.Omega0 - kspace_vars.OmegaNu + OmegaNua3;
  const double fnu = OmegaNua3/Omega0;
  for (ik = 0; ik < nk_in; ik++)
           wavenum[ik]=exp(logk[ik]);
  /* Get a delta_nu_curr from CAMB.*/
  if(!d_tot->delta_tot_init_done){
      terminate("Should have called delta_tot_init first\n");
       /*Initialise delta_tot, setting ia = 1
        * and signifying that we are ready to leave the relativistic regime*/
//        delta_tot_init(d_tot, nk_in, wavenum,delta_cdm_curr, ThisTask, &transfer_init);
  }

  /*If we get called twice with the same scale factor, do nothing*/
  if(log(a)-d_tot->scalefact[d_tot->ia-1] < FLOAT){
       for (ik = 0; ik < d_tot->nk; ik++)
               delta_nu_curr[ik] = d_tot->delta_nu_last[ik];
       return;
  }

  /* We need some estimate for delta_tot(current time) to obtain delta_nu(current time).
     Even though delta_tot(current time) is not directly used (the integrand vanishes at a = a(current)),
     it is indeed needed for interpolation */
  /* It was checked that using delta_tot(current time) = delta_cdm(current time) leads to no more than 2%
     error on delta_nu (and moreover for large k). A simple estimate for delta_nu decreases the maximum
     relative error on delta_nu to ~1E-4. So we only need one step. */
   d_tot->scalefact[d_tot->ia] = log(a);
   /* First estimate for delta_tot(a) */
   for (ik = 0; ik < d_tot->nk; ik++) {
      /*Guard against floating point zeros*/
      d_tot->delta_tot[ik][d_tot->ia] = (1.-fnu)*delta_cdm_curr[ik]+fnu*d_tot->delta_nu_last[ik];
   }
   /*Get the new delta_nu_curr*/
   get_delta_nu_combined(d_tot, a, d_tot->ia+1, wavenum, delta_nu_curr);
   /*Update delta_nu_last*/
   for (ik = 0; ik < d_tot->nk; ik++)
       d_tot->delta_nu_last[ik]=delta_nu_curr[ik];
   /* Decide whether we save the current time or not */
   if (a > exp(d_tot->scalefact[d_tot->ia-1]) + 0.01) {
       /* If so update delta_tot(a) correctly */
       for (ik = 0; ik < d_tot->nk; ik++){
               d_tot->delta_tot[ik][d_tot->ia] = fnu*delta_nu_curr[ik]+(1.-fnu)*delta_cdm_curr[ik];
       }
       d_tot->ia++;
#ifdef HYBRID_NEUTRINOS
    //Check whether we want to stop the particle neutrinos from being tracers.
    set_slow_neutrinos_analytic();
#endif
       /*printf("Updating delta_tot: a=%f, Na=%d, last=%f\n",a,ia,exp(scalefact[ia-2]));*/
       if(d_tot->ThisTask==0)
        save_delta_tot(d_tot, d_tot->ia-1, NULL);
   }
   return;
}

/*Save a single delta_nu power spectrum into a file*/
void save_delta_tot(_delta_tot_table *d_tot, int iia, char * savedir)
{
        FILE *fd;
        int i;
        char * dfile;
        //NULL means use current directory
        if (savedir == NULL){
            dfile = "delta_tot_nu.txt";
        }
        else{
            int nbytes = sizeof(char)*(strlen(savedir)+25);
            dfile = mymalloc("filename", nbytes);
            if(!dfile){
                    char err[150];
                    snprintf(err,150,"Unable to allocate %d bytes for filename\n",nbytes);
                    terminate(err);
            }
            dfile = strncpy(dfile, savedir, nbytes);
            dfile = strncat(dfile, "delta_tot_nu.txt",25);
        }
        if(!(fd = fopen(dfile, "a")))
             terminate("Could not open delta_tot_nu.txt for writing!\n");
        /*Write log scale factor*/
        fprintf(fd, "# %le ", d_tot->scalefact[iia]);
        /*Write kvalues*/
        for(i=0;i<d_tot->nk; i++)
                fprintf(fd,"%le ",d_tot->delta_tot[i][iia]);
        fprintf(fd,"\n");
        fclose(fd);
        if(savedir != NULL)
            myfree(dfile);
   return;
}

/* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
 * Must be called before delta_tot_init, or resuming wont work*/
void read_all_nu_state(_delta_tot_table *d_tot, char * savedir, double Time)
{
    FILE* fd;
    char * dfile;
    if (savedir == NULL){
        dfile = "delta_tot_nu.txt";
    }
    else{
        int nbytes = sizeof(char)*(strlen(savedir)+25);
        dfile = mymalloc("filename", nbytes);
        if(!dfile){
                char err[150];
                snprintf(err,150,"Unable to allocate %d bytes for filename\n",nbytes);
                terminate(err);
        }
        dfile = strncpy(dfile, savedir, nbytes);
        dfile = strncat(dfile, "delta_tot_nu.txt",25);
    }
    /*Load delta_tot from a file, if such a file exists. Allows resuming.*/
    if((fd = fopen(dfile, "r"))){
            /*Read redshifts; Initial one is known already*/
            int iia;
            for(iia=0; iia< d_tot->namax;iia++){
                    double scale;
                    int ik;
                    if(fscanf(fd, "# %lg ", &scale) != 1)
                            break;
                    d_tot->scalefact[iia]=scale;
                    /*Only read until we reach the present day*/
                    if(log(Time) <= scale)
                            break;
                    /*Read kvalues*/
                    /*If we do not have a complete delta_tot for one redshift, we stop
                    * unless this is the first line, in which case we use it to set nk */
                    for(ik=0;ik<d_tot->nk; ik++){
                            if(fscanf(fd, "%lg ", &(d_tot->delta_tot[ik][iia])) != 1){
                                if(iia != 0){
                                    char err[150];
                                    snprintf(err,150,"Incomplete delta_tot in delta_tot_nu.txt; a=%g\n",exp(d_tot->scalefact[iia]));
                                    terminate(err);
                                }
                                else{
                                    d_tot->nk = ik+1;
                                    break;
                                }
                            }
                    }
            }
            /*If our table starts at a different time from the simulation, stop.*/
            if(fabs(d_tot->scalefact[0] - log(kspace_params.TimeTransfer)) > 1e-4){
                    char err[250];
                    snprintf(err,250," delta_tot_nu.txt starts wih a=%g, transfer function is at a=%g\n",exp(d_tot->scalefact[0]),kspace_params.TimeTransfer);
                    terminate(err);
            }

            if(iia > 0)
                    d_tot->ia=iia;
            printf("Read %d stored power spectra from delta_tot_nu.txt\n",iia);
            fclose(fd);
    }
    if (savedir != NULL)
        myfree(dfile);
}

/* Function to save all the internal state of the neutrino integrator to disc.
 * Must be called for resume to work*/
void save_all_nu_state(_delta_tot_table *d_tot, char * savedir)
{
    int ik;
    for(ik=0; ik< d_tot->ia; ik++)
         save_delta_tot(d_tot, ik, savedir);
}

/*What follows are private functions for the integration routine get_delta_nu*/

#ifdef HYBRID_NEUTRINOS
//Fraction of neutrino mass not followed by the analytic integrator.
static double nufrac_low=0;
double _nufrac_low(double mnu);
double Jfrac_high(double x, double mnu);
int set_slow_neutrinos_analytic();
#endif

double fslength(double ai, double af,double mnu);
double find_ai(double af, double k, double kfsl,double mnu);
double specialJ(double x, double mnu);
double specialJ_fit(double x);

/*Returns kT / a M_nu (which is dimensionless) in the relativistic limit
 * where it is kT / (a^2 m_nu^2 + (kT)^2)^(1/2)
 * Takes a* M_nu as argument. */
inline double kTbyaM(double amnu)
{
  const double kT=BOLEVK*TNU;
  return kT/amnu;
}

/*Kernel function for the fslength integration*/
double fslength_int(double loga, void *params)
{
    double mnu = *((double *)params);
    double a = exp(loga);
    return kTbyaM(a*mnu)/(a*hubble_function(a));
}

/******************************************************************************************************
Free-streaming length for a non-relativistic particle of momentum q = T0, from scale factor ai to af.
Result is in Unit_Length.
******************************************************************************************************/

double fslength(double ai, double af,double mnu)
{
  double abserr;
  double fslength_val;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
  F.function = &fslength_int;
  F.params = &mnu;
  gsl_integration_qag (&F, log(ai), log(af), 0, 1e-6,GSL_VAL,6,w,&(fslength_val), &abserr);
  gsl_integration_workspace_free (w);
  return LIGHT*fslength_val;
}

/**************************************************************************************************
Fit to the special function J(x) that is accurate to better than 3% relative and 0.07% absolute
    J(x) = Integrate[(Sin[q*x]/(q*x))*(q^2/(Exp[q] + 1)), {q, 0, Infinity}]
    and J(0) = 1.
    Mathematica gives this in terms of the PolyGamma function:
   (PolyGamma[1, 1/2 - i x/2] - PolyGamma[1, 1 - i x/2] -    PolyGamma[1, 1/2 + i x/2] +
   PolyGamma[1, 1 + i x/2])/(12 x Zeta[3]), which we could evaluate exactly if we wanted to.
***************************************************************************************************/
double specialJ_fit(double x)
{

  double x2, x4, x8;
  if (x <= 0.)
      return 1.;
  x2 = x*x;
  x4 = x2*x2;
  x8 = x4*x4;

  return (1.+ 0.0168 * x2 + 0.0407* x4)/(1. + 2.1734 * x2 + 1.6787 * exp(4.1811*log(x)) +  0.1467 * x8);
}

#ifdef HYBRID_NEUTRINOS
/*Function to decide whether slow neutrinos are treated analytically or not.
 * Sets kspace_vars.slow_neutrinos_analytic.
 * Should be called every few timesteps.*/
int set_slow_neutrinos_analytic()
{
    double val=1;
    //Just use a redshift cut for now. Really we want something more sophisticated,
    //based on the shot noise and average overdensity.
    if (kspace_vars.Time > kspace_vars.nu_crit_time){
        if(kspace_vars.slow_neutrinos_analytic && d_tot->ThisTask==0)
            printf("Particle neutrinos start to gravitate NOW: nufrac_low is: %g\n",nufrac_low);
        val = 0;
    }
    kspace_vars.slow_neutrinos_analytic = val;
    return val;
}


//Fermi-Dirac kernel for below
double fermi_dirac_kernel(double x, void * params)
{
  return x * x / (exp(x) + 1);
}

//This is integral f_0(q) q^2 dq between 0 and qc to compute the fraction of OmegaNu which is in particles.
double _nufrac_low(double mnu)
{
    double qc = mnu * kspace_vars.vcrit/LIGHT/(BOLEVK*TNU);
    /*These functions are so smooth that we don't need much space*/
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
    double abserr;
    gsl_function F;
    F.function = &fermi_dirac_kernel;
    F.params = NULL;
    double total_fd;
    gsl_integration_qag (&F, 0, qc, 0, 1e-6,100,GSL_INTEG_GAUSS61, w,&(total_fd), &abserr);
    //divided by the total F-D probability (which is 3 Zeta(3)/2 ~ 1.8 if MAX_FERMI_DIRAC is large enough
    total_fd /= 1.5*1.202056903159594;
    gsl_integration_workspace_free (w);
    return total_fd;
}

//Asymptotic series expansion from YAH. Not good when qc * x is small, but fine otherwise.
double II(double x, double qc, int n)
{
    return (n*n+n*n*n*qc+n*qc*x*x - x*x)* qc*gsl_sf_bessel_j0(qc*x) + (2*n+n*n*qc+qc*x*x)*cos(qc*x);
}

//This is an approximation to integral f_0(q) q^2 j_0(qX) dq between qc and infinity.
//It gives the fraction of the integral that is due to neutrinos above a certain threshold.
double Jfrac_high(double x, double mnu)
{
    double qc = mnu * kspace_vars.vcrit/LIGHT/(BOLEVK*TNU);
    double integ=0;
    for(int n=1; n<20; n++)
    {
        integ+= pow((-1),n)*exp(-n*qc)/(n*n+x*x)/(n*n+x*x)*II(x,qc,n);
    }
    //Normalise with integral(f_0(q)q^2 dq), same as I(X). So that as qc-> infinity, this -> specialJ_fit(x)
    integ /= 1.8031;
    return integ;
}

/* Fourier transform of truncated Fermi Dirac distribution, with support on q > qc only.
    qc is a dimensionless momentum (normalized to TNU),
    Uses an expansion of the FD distribution for qc << 1, very accurate for qc <= 1
    mnu is in eV. x has units of inverse dimensionless momentum
 */
double specialJ(double x, double mnu)
{
  if( !kspace_vars.slow_neutrinos_analytic ) {
   return Jfrac_high(x, mnu);
  }
  else {
    return specialJ_fit(x);
  }
}

#else //Now for single-component neutrinos

double specialJ(double x, double mnu)
{
    return specialJ_fit(x);
}

#endif //HYBRID_NEUTRINOS

/*A structure for the parameters for the below integration kernel*/
struct _delta_nu_int_params
{
    /*This is the current time, not the time at which
     * this contribution was seen by the neutrinos*/
    double a;
    /*Current wavenumber*/
    double k;
    double mnu;
    gsl_interp_accel *acc;
    gsl_interp *spline;
    /*Make sure this is at the same k as above*/
    double * delta_tot;
    double * scale;
};
typedef struct _delta_nu_int_params delta_nu_int_params;

/*Integration kernel for below*/
double get_delta_nu_int(double logai, void * params)
{
    delta_nu_int_params * p = (delta_nu_int_params *) params;
    double ai = exp(logai);
    double fsl_aia = fslength(ai, p->a,p->mnu);
    double delta_tot_at_a = gsl_interp_eval(p->spline,p->scale,p->delta_tot,logai,p->acc);
    return fsl_aia/(ai*hubble_function(ai)) *specialJ(p->k*fsl_aia, p->mnu) * delta_tot_at_a/(ai*kTbyaM(ai*p->mnu));
}


/*****************************************************************************************************
Main function: given tables of wavenumbers, total delta at Na earlier times (<= a),
and initial conditions for neutrinos, computes the current delta_nu.
Na is the number of currently stored time steps.
Requires transfer_init_tabulate to have been called prior to first call.
******************************************************************************************************/

void get_delta_nu(_delta_tot_table * d_tot, double a, int Na, double wavenum[],  double delta_nu_curr[],double mnu)
{
  double fsl_A0a,deriv_prefac;
  int ik;
  if(d_tot->ThisTask == 0)
          printf("Start get_delta_nu: a=%g Na =%d wavenum[0]=%g delta_tot[0]=%g m_nu=%g\n",a,Na,wavenum[0],d_tot->delta_tot[0][Na-1],mnu);

  fsl_A0a = fslength(A0, a,mnu);
  /*Precompute factor used to get delta_nu_init. This assumes that delta ~ a, so delta-dot is roughly 1.*/
  deriv_prefac = A0*(hubble_function(A0)/LIGHT)/kTbyaM(A0*mnu);

  for (ik = 0; ik < d_tot->nk; ik++) {
      /* Initial condition piece, assuming linear evolution of delta with a up to startup redshift */
      /* This assumes that delta ~ a, so delta-dot is roughly 1. */
      /* Also ignores any difference in the transfer functions between species.
       * This will be good if all species have similar masses, or
       * if two species are massless.
       * Also, since at early times the clustering is tiny, it is very unlikely to matter.*/
      delta_nu_curr[ik] = specialJ(wavenum[ik]*fsl_A0a, mnu)*d_tot->delta_nu_init[ik] *(1.+ deriv_prefac*fsl_A0a);
  }
  /*If only one time given, we are still at the initial time*/
  if(Na > 1){
        delta_nu_int_params params;
        params.acc = gsl_interp_accel_alloc();
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
        gsl_function F;
        F.function = &get_delta_nu_int;
        F.params=&params;
        /*Use cubic interpolation*/
        if(Na > 2)
                params.spline=gsl_interp_alloc(gsl_interp_cspline,Na);
        /*Unless we have only two points*/
        else
                params.spline=gsl_interp_alloc(gsl_interp_linear,Na);
        if(!params.spline || !params.acc || !w )
              terminate("Error initialising and allocating memory for gsl interpolator and integrator.\n");
        params.a=a;
        params.scale=d_tot->scalefact;
        params.mnu=mnu;
        for (ik = 0; ik < d_tot->nk; ik++) {
            double abserr,d_nu_tmp;
            params.k=wavenum[ik];
            params.delta_tot=d_tot->delta_tot[ik];
            gsl_interp_init(params.spline,params.scale,params.delta_tot,Na);
            gsl_integration_qag (&F, log(A0), log(a), 0, 1e-6,GSL_VAL,6,w,&d_nu_tmp, &abserr);
            delta_nu_curr[ik] += d_tot->delta_nu_prefac * d_nu_tmp;
         }
         gsl_integration_workspace_free (w);
         gsl_interp_free(params.spline);
         gsl_interp_accel_free(params.acc);
   }
   if(d_tot->ThisTask == 0){
          printf("delta_nu_curr[0] is %g\n",delta_nu_curr[0]);
          for(ik=0; ik< 5; ik++)
          printf("k %g d_nu %g\n",wavenum[40*ik], delta_nu_curr[40*ik]);
   }
   return;
}
