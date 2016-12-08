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
#include <unistd.h>

#include "delta_tot_table.h"
#include "gadget_defines.h"
#include "kspace_neutrino_const.h"

/*Allocate memory for delta_tot_table. This is separate from delta_tot_init because we need to allocate memory
 * before we have the information needed to initialise it*/
void allocate_delta_tot_table(_delta_tot_table *d_tot, const int nk_in, const double TimeTransfer, const double TimeMax, const double Omega0, const _omega_nu * const omnu, const double UnitTime_in_s, const double UnitLength_in_cm, int debug)
{
   int count;
   /*Memory allocations need to be done on all processors*/
   d_tot->nk_allocated=nk_in;
   d_tot->nk=nk_in;
   /*Store starting time*/
   d_tot->TimeTransfer = TimeTransfer;
   /* Allocate memory for delta_tot here, so that we can have further memory allocated and freed
    * before delta_tot_init is called. The number nk here should be larger than the actual value needed.*/
   /*Allocate pointers to each k-vector*/
   d_tot->namax=ceil(100*(TimeMax-TimeTransfer))+2;
   d_tot->ia=0;
   d_tot->delta_tot =(double **) mymalloc("kspace_delta_tot",nk_in*sizeof(double *));
   /*Allocate list of scale factors, and space for delta_tot, in one operation.*/
   d_tot->scalefact = (double *) mymalloc("kspace_scalefact",d_tot->namax*(nk_in+1)*sizeof(double));
   /*Allocate actual data. Note that this means data can be accessed either as:
    * delta_tot[k][a] OR as
    * delta_tot[0][a+k*namax] */
   d_tot->delta_tot[0] = d_tot->scalefact+d_tot->namax;
   for(count=1; count< nk_in; count++)
        d_tot->delta_tot[count] = d_tot->delta_tot[0] + count*d_tot->namax;
   /*Allocate space for the initial neutrino power spectrum*/
   d_tot->delta_nu_init =(double *) mymalloc("kspace_delta_nu_init",2*nk_in*sizeof(double));
   d_tot->delta_nu_last=d_tot->delta_nu_init+nk_in;
   /*Setup pointer to the matter density*/
   d_tot->omnu = omnu;
   /*Set the prefactor for delta_nu, and the units system*/
   d_tot->light = LIGHTCGS * UnitTime_in_s/UnitLength_in_cm;
   d_tot->delta_nu_prefac = 1.5 *Omega0 * HUBBLE * HUBBLE * pow(UnitTime_in_s,2)/d_tot->light;
   /*Matter fraction excluding neutrinos*/
   d_tot->Omeganonu = Omega0 - get_omega_nu(omnu, 1);
   /*Whether we save intermediate files and output diagnostics*/
   d_tot->debug = debug;
}

/*Free memory for delta_tot_table.*/
void free_delta_tot_table(_delta_tot_table *d_tot)
{
    myfree(d_tot->delta_tot);
    myfree(d_tot->scalefact);
    myfree(d_tot->delta_nu_init);
}

void handler (const char * reason, const char * file, int line, int gsl_errno)
{
    printf("GSL_ERROR in file: %s, line %d, errno:%d\n",file, line, gsl_errno);
    terminate(reason);
}

/* Constructor. transfer_init_tabulate must be called before this function.
 * Initialises delta_tot (including from a file) and delta_nu_init from the transfer functions.
 * read_all_nu_state must be called before this if you want reloading from a snapshot to work
 * Note delta_cdm_curr includes baryons, and is only used if not resuming.*/
void delta_tot_init(_delta_tot_table * const d_tot, const int nk_in, const double wavenum[], const double delta_cdm_curr[], const _transfer_init_table * const t_init)
{
    int ik;
    if(nk_in > d_tot->nk_allocated){
           char err[500];
           sprintf(err,"input power of %d is longer than memory of %d\n",nk_in,d_tot->nk_allocated);
           terminate(err);
    }
    gsl_set_error_handler(handler);
    d_tot->nk=nk_in;
    /*Construct delta_nu_init from the transfer functions.*/
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp *spline;
    const double OmegaNua3=get_omega_nu(d_tot->omnu, d_tot->TimeTransfer)*pow(d_tot->TimeTransfer,3);
    /*Initialise delta_nu_init to use the first timestep's delta_cdm_curr
     * so that it includes potential Rayleigh scattering. */
    spline=gsl_interp_alloc(gsl_interp_cspline,t_init->NPowerTable);
    gsl_interp_init(spline,t_init->logk,t_init->T_nu,t_init->NPowerTable);
    for(ik=0;ik<d_tot->nk;ik++){
            double T_nubyT_notnu = gsl_interp_eval(spline,t_init->logk,t_init->T_nu,log(wavenum[ik]),acc);
            /*The total power spectrum using neutrinos and radiation from the CAMB transfer functions:
             * The CAMB transfer functions are defined such that
             * delta_cdm ~ T_cdm (and some other constant factors)
             * then delta_t = (Omega_cdm delta_cdm + Omega_nu delta_nu)/(Omega_cdm + Omega_nu)
             *          = delta_cdm (Omega_cdm+ Omega_nu (delta_nu/delta_cdm)) / (Omega_cdm +Omega_nu)
             *          = delta_cdm (Omega_cdm+ Omega_nu (delta_nu/delta_cdm)) / (Omega_cdm+Omega_nu) */
            const double OmegaMa = (d_tot->Omeganonu+OmegaNua3);
            const double fnu = OmegaNua3/OmegaMa;
            const double CDMtoTot=1-fnu+T_nubyT_notnu*fnu;
            /*If we are not restarting, initialise the first delta_tot and delta_nu_init*/
            if(d_tot->ia == 0)
                d_tot->delta_tot[ik][0] = delta_cdm_curr[ik]*CDMtoTot;
            /* Also initialise delta_nu_init here to save time later.
             * Use the first delta_tot, in case we are resuming.*/
            d_tot->delta_nu_init[ik] = d_tot->delta_tot[ik][0]/CDMtoTot*fabs(T_nubyT_notnu);
    }
    gsl_interp_accel_free(acc);
    gsl_interp_free(spline);

    /*If we are not restarting, make sure we set the scale factor*/
    if(d_tot->ia == 0) {
        d_tot->scalefact[0]=log(d_tot->TimeTransfer);
        d_tot->ia=1;
    }
    if(d_tot->ThisTask==0 && d_tot->debug){
        save_all_nu_state(d_tot, NULL);
    }
    /*Initialise delta_nu_last*/
    get_delta_nu_combined(d_tot, exp(d_tot->scalefact[d_tot->ia-1]), wavenum, d_tot->delta_nu_last);
    d_tot->delta_tot_init_done=1;
    return;
}

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(const _delta_tot_table * const d_tot, const double a, const double wavenum[],  double delta_nu_curr[])
{
    const double Omega_nu_tot=get_omega_nu_nopart(d_tot->omnu, a);
    int mi;
    /*Initialise delta_nu_curr*/
    memset(delta_nu_curr, 0, d_tot->nk*sizeof(double));
    /*Get each neutrinos species and density separately and add them to the total.
     * Neglect perturbations in massless neutrinos.*/
    for(mi=0; mi<NUSPECIES; mi++) {
            if(d_tot->omnu->nu_degeneracies[mi] > 0) {
                 int ik;
                 double delta_nu_single[d_tot->nk];
                 const double omeganu = d_tot->omnu->nu_degeneracies[mi] * omega_nu_single(d_tot->omnu, a, mi);
                 get_delta_nu(d_tot, a, wavenum, delta_nu_single,d_tot->omnu->RhoNuTab[mi]->mnu);
                 for(ik=0; ik<d_tot->nk; ik++)
                    delta_nu_curr[ik]+=delta_nu_single[ik]*omeganu/Omega_nu_tot;
            }
    }
    return;
}

/*Update the last value of delta_tot in the table with a new value computed
 from the given delta_cdm_curr and delta_nu_curr.
 If overwrite is true, overwrite the existing final entry.*/
void update_delta_tot(_delta_tot_table * const d_tot, const double a, const double delta_cdm_curr[], const double delta_nu_curr[], const int overwrite)
{
  const double OmegaNua3=get_omega_nu_nopart(d_tot->omnu, a)*pow(a,3);
  const double OmegaMa = d_tot->Omeganonu + get_omega_nu(d_tot->omnu, a)*pow(a,3);
  const double fnu = OmegaNua3/OmegaMa;
  int ik;
  if(!overwrite)
    d_tot->ia++;
  /*Update the scale factor*/
  d_tot->scalefact[d_tot->ia-1] = log(a);
  /* Update delta_tot(a)*/
  for (ik = 0; ik < d_tot->nk; ik++){
    d_tot->delta_tot[ik][d_tot->ia-1] = fnu*delta_nu_curr[ik]+(1.-fnu)*delta_cdm_curr[ik];
  }
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
 * @param d_tot contains the state of the integrator; samples of the total power spectrum at earlier times.
 * @param a is the current scale factor
 * @param nk_in is the number of k bins in delta_cdm_curr and keff.
 * @param keff is an array of length nk containing (natural) log k
 * @param delta_cdm_curr is an array of length nk containing the square root of the current cdm power spectrum
 * @param delta_nu_curr is an array of length nk which stores the square root of the current neutrino power spectrum. Main output of the function.
******************************************************************************************************/
void get_delta_nu_update(_delta_tot_table * const d_tot, const double a, const int nk_in, const double keff[], const double delta_cdm_curr[], double delta_nu_curr[])
{
  int ik;
  /* Get a delta_nu_curr from CAMB.*/
  if(!d_tot->delta_tot_init_done)
      terminate("Should have called delta_tot_init first\n");
  if(nk_in != d_tot->nk)
      terminate("Number of kbins differs from stored delta_tot\n");
  if(d_tot->nk < 2){
      char err[150];
      snprintf(err,150,"Number of kbins is unreasonably small: %d\n",d_tot->nk);
      terminate(err);
  }
  /*If we get called twice with the same scale factor, do nothing*/
  if(log(a)-d_tot->scalefact[d_tot->ia-1] < FLOAT){
       for (ik = 0; ik < d_tot->nk; ik++)
               delta_nu_curr[ik] = d_tot->delta_nu_last[ik];
       return;
  }

   /*We need some estimate for delta_tot(current time) to obtain delta_nu(current time).
     Even though delta_tot(current time) is not directly used (the integrand vanishes at a = a(current)),
     it is indeed needed for interpolation */
   /*It was checked that using delta_tot(current time) = delta_cdm(current time) leads to no more than 2%
     error on delta_nu (and moreover for large k). A simple estimate for delta_nu decreases the maximum
     relative error on delta_nu to ~1E-4. So we only need one step. */
   /*This increments the number of stored spectra, although the last one is not yet final.*/
   update_delta_tot(d_tot, a, delta_cdm_curr, d_tot->delta_nu_last, 0);
   /*Get the new delta_nu_curr*/
   get_delta_nu_combined(d_tot, a, keff, delta_nu_curr);
   /*Update delta_nu_last*/
   for (ik = 0; ik < d_tot->nk; ik++)
       d_tot->delta_nu_last[ik]=delta_nu_curr[ik];
   /* Decide whether we save the current time or not */
   if (a >= exp(d_tot->scalefact[d_tot->ia-2]) + 0.009) {
       /* If so update delta_tot(a) correctly, overwriting current power spectrum */
       update_delta_tot(d_tot, a, delta_cdm_curr, delta_nu_curr, 1);
       /*printf("Updating delta_tot: a=%f, Na=%d, last=%f\n",a,ia,exp(scalefact[ia-2]));*/
       if(d_tot->ThisTask==0 && d_tot->debug)
          save_delta_tot(d_tot, d_tot->ia-1, NULL);
   }
   /*Otherwise discard the last powerspectrum*/
   else
       d_tot->ia--;
   return;
}

/* Reads data from snapdir / delta_tot_nu.txt into delta_tot, if present.
 * Must be called before delta_tot_init, or resuming wont work*/
void read_all_nu_state(_delta_tot_table * const d_tot, const char * savedir)
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
    fd = fopen(dfile, "r");
    if(!fd) {
        return;
    }
    /*Read redshifts; Initial one is known already*/
    int iia;
    for(iia=0; iia< d_tot->namax;iia++){
            double scale;
            int ik;
            if(fscanf(fd, "# %lg ", &scale) != 1)
                    break;
            d_tot->scalefact[iia]=scale;
            /*Read kvalues*/
            /*If we do not have a complete delta_tot for one redshift, we stop
            * unless this is the first line, in which case we use it to set nk */
            for(ik=0;ik<d_tot->nk; ik++){
                    if(fscanf(fd, "%lg ", &(d_tot->delta_tot[ik][iia])) != 1){
                        if(iia != 0){
                            char err[150];
                            snprintf(err,150,"Expected %d k values, got %d for delta_tot in %s; a=%g\n",d_tot->nk, ik, dfile, exp(d_tot->scalefact[iia]));
                            terminate(err);
                        }
                        else{
                            d_tot->nk = ik;
                            break;
                        }
                    }
            }
    }
    /*If our table starts at a different time from the simulation, stop.*/
    if(fabs(d_tot->scalefact[0] - log(d_tot->TimeTransfer)) > 1e-4){
            char err[250];
            snprintf(err,250,"%s starts wih a=%g, transfer function is at a=%g\n",dfile, exp(d_tot->scalefact[0]),d_tot->TimeTransfer);
            terminate(err);
    }

    if(iia > 0)
            d_tot->ia=iia;
    if(d_tot->debug)
        printf("Read %d stored power spectra from %s\n",iia, dfile);
    fclose(fd);
    if (savedir != NULL)
        myfree(dfile);
}

/*Save a single delta_nu power spectrum into a file*/
void save_delta_tot(const _delta_tot_table * const d_tot, const int iia, char * savefile)
{
    FILE *fd;
    int i;
    char * dfile;
    /*NULL means use current directory*/
    if (savefile == NULL){
        dfile = "delta_tot_nu.txt";
    }
    else {
        dfile = savefile;
    }
    if(!(fd = fopen(dfile, "a"))) {
            char err[300];
            snprintf(err,300,"Could not open %s for writing!\n",dfile);
            terminate(err);
    }
    /*Write log scale factor*/
    fprintf(fd, "# %le ", d_tot->scalefact[iia]);
    /*Write kvalues*/
    for(i=0;i<d_tot->nk; i++)
            fprintf(fd,"%le ",d_tot->delta_tot[i][iia]);
    fprintf(fd,"\n");
    fclose(fd);
    return;
}

/* Function to save all the internal state of the neutrino integrator to disc.
 * Must be called for resume to work*/
void save_all_nu_state(const _delta_tot_table * const d_tot, char * savedir)
{

    int ik;
    char * savefile;
    if(savedir) {
            int nbytes = sizeof(char)*(strlen(savedir)+25);
            savefile = mymalloc("filename", nbytes);
            if(!savefile){
                    char err[150];
                    snprintf(err,150,"Unable to allocate %d bytes for filename\n",nbytes);
                    terminate(err);
            }
            savefile = strncpy(savefile, savedir, nbytes);
            savefile = strncat(savefile, "delta_tot_nu.txt",25);
    }
    else {
        savefile = "delta_tot_nu.txt";
    }
    /*Get a clean debug restart file*/
    /*Check whether old file exists*/
    if(access( savefile, F_OK ) != -1 ) {
        /*If it does make the new file name and rename the file*/
        int nbytes = sizeof(char)*(strlen(savefile)+6);
        char * bak_savefile = mymalloc("filename2", nbytes);
        if(bak_savefile) {
            bak_savefile = strncpy(bak_savefile, savefile, strlen(savefile)+1);
            bak_savefile = strncat(bak_savefile, ".bak",6);
            rename(savefile, bak_savefile);
            myfree(bak_savefile);
        }
    }
    for(ik=0; ik< d_tot->ia; ik++)
         save_delta_tot(d_tot, ik, savefile);
    if(savedir)
        myfree(savefile);
}

/*What follows are private functions for the integration routine get_delta_nu*/

/*Kernel function for the fslength integration*/
double fslength_int(const double loga, void *params)
{
    /*This should be M_nu / k_B T_nu (which is dimensionless)*/
    double mnubykT = *((double *)params);
    const double a = exp(loga);
    return 1./a/mnubykT/(a*hubble_function(a));
}

/******************************************************************************************************
Free-streaming length for a non-relativistic particle of momentum q = T0, from scale factor ai to af.
Arguments:
logai - log of initial scale factor
logaf - log of final scale factor
mnubykT - M_nu / k_B T_nu (dimensionless)
light - speed of light in internal units.
Result is in Unit_Length/Unit_Time.
******************************************************************************************************/
double fslength(const double logai, const double logaf, const double mnubykT, const double light)
{
  double abserr;
  double fslength_val;
  /*This is to avoid a gcc warning: it is still const really.*/
  double mnu_ncst = mnubykT;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
  F.function = &fslength_int;
  F.params = &mnu_ncst;
  if(logai >= logaf)
      return 0;
  gsl_integration_qag (&F, logai, logaf, 0, 1e-6,GSL_VAL,6,w,&(fslength_val), &abserr);
  gsl_integration_workspace_free (w);
  return light*fslength_val;
}

/**************************************************************************************************
Fit to the special function J(x) that is accurate to better than 3% relative and 0.07% absolute
    J(x) = Integrate[(Sin[q*x]/(q*x))*(q^2/(Exp[q] + 1)), {q, 0, Infinity}]
    and J(0) = 1.
    Mathematica gives this in terms of the PolyGamma function:
   (PolyGamma[1, 1/2 - i x/2] - PolyGamma[1, 1 - i x/2] -    PolyGamma[1, 1/2 + i x/2] +
   PolyGamma[1, 1 + i x/2])/(12 x Zeta[3]), which we could evaluate exactly if we wanted to.
***************************************************************************************************/
inline double specialJ_fit(const double x)
{

  double x2, x4, x8;
  if (x <= 0.)
      return 1.;
  x2 = x*x;
  x4 = x2*x2;
  x8 = x4*x4;

  return (1.+ 0.0168 * x2 + 0.0407* x4)/(1. + 2.1734 * x2 + 1.6787 * exp(4.1811*log(x)) +  0.1467 * x8);
}

/*Asymptotic series expansion from YAH. Not good when qc * x is small, but fine otherwise.*/
inline double II(const double x, const double qc, const int n)
{
    return (n*n+n*n*n*qc+n*qc*x*x - x*x)* qc*gsl_sf_bessel_j0(qc*x) + (2*n+n*n*qc+qc*x*x)*cos(qc*x);
}

/* Fourier transform of truncated Fermi Dirac distribution, with support on q > qc only.
 * qc is a dimensionless momentum (normalized to TNU),
 * mnu is in eV. x has units of inverse dimensionless momentum
 * This is an approximation to integral f_0(q) q^2 j_0(qX) dq between qc and infinity.
 * It gives the fraction of the integral that is due to neutrinos above a certain threshold.
 * Arguments: vcmnu is vcrit*mnu/LIGHT */
inline double Jfrac_high(const double x, const double qc)
{
    double integ=0;
    int n;
    for(n=1; n<20; n++)
    {
        integ+= -1*pow((-1),n)*exp(-n*qc)/(n*n+x*x)/(n*n+x*x)*II(x,qc,n);
    }
    /*Normalise with integral(f_0(q)q^2 dq), same as I(X). So that as qc-> infinity, this -> specialJ_fit(x)*/
    integ /= 1.8031;
    return integ;
}

/*Function that picks whether to use the truncated integrator or not*/
double specialJ(const double x, const double qc)
{
  if( qc > 0 ) {
   return Jfrac_high(x, qc);
  }
  return specialJ_fit(x);
}

/*A structure for the parameters for the below integration kernel*/
struct _delta_nu_int_params
{
    /*Current wavenumber*/
    double k;
    /*Neutrino mass divided by k_B T_nu*/
    double mnubykT;
    gsl_interp_accel *acc;
    gsl_interp *spline;
    /*Precomputed free-streaming lengths*/
    gsl_interp_accel *fs_acc;
    gsl_interp *fs_spline;
    double * fslengths;
    double * fsscales;
    /*Make sure this is at the same k as above*/
    double * delta_tot;
    double * scale;
    /* qc is a dimensionless momentum (normalized to TNU): v_c * mnu / (k_B * T_nu).
     * This is the critical momentum for hybrid neutrinos: it is unused if
     * hybrid neutrinos are not defined, but left here to save ifdefs.*/
    double qc;
};
typedef struct _delta_nu_int_params delta_nu_int_params;

/*Integration kernel for below*/
double get_delta_nu_int(double logai, void * params)
{
    delta_nu_int_params * p = (delta_nu_int_params *) params;
    double fsl_aia = gsl_interp_eval(p->fs_spline,p->fsscales,p->fslengths,logai,p->fs_acc);
    double delta_tot_at_a = gsl_interp_eval(p->spline,p->scale,p->delta_tot,logai,p->acc);
    double specJ = specialJ(p->k*fsl_aia, p->qc);
    double ai = exp(logai);
    return fsl_aia/(ai*hubble_function(ai)) * specJ * delta_tot_at_a * p->mnubykT;
}

/*****************************************************************************************************
Main function: given tables of wavenumbers, total delta at Na earlier times (<= a),
and initial conditions for neutrinos, computes the current delta_nu.
Na is the number of currently stored time steps.
Requires transfer_init_tabulate to have been called prior to first call.
******************************************************************************************************/

void get_delta_nu(const _delta_tot_table * const d_tot, const double a, const double wavenum[], double delta_nu_curr[],const double mnu)
{
  double fsl_A0a,deriv_prefac;
  int ik;
  /* Variable is unused unless we have hybrid neutrinos,
   * but we define it anyway to save ifdeffing later.*/
  double qc = 0;
  /*Number of stored power spectra. This includes the initial guess for the next step*/
  const int Na = d_tot->ia;
  const double mnubykT = mnu /d_tot->omnu->kBtnu;
  if(d_tot->ThisTask == 0 && d_tot->debug)
          printf("Start get_delta_nu: a=%g Na =%d wavenum[0]=%g delta_tot[0]=%g m_nu=%g\n",a,Na,wavenum[0],d_tot->delta_tot[0][Na-1],mnu);

  fsl_A0a = fslength(log(d_tot->TimeTransfer), log(a),mnubykT, d_tot->light);
   /* Check whether the particle neutrinos are active at this point.
    * If they are we want to truncate our integration.
    * Only do this is hybrid neutrinos are activated in the param file.*/
   if(particle_nu_fraction(&d_tot->omnu->hybnu, a, 0) > 0) {
        qc = (d_tot->omnu->hybnu.vcrit / d_tot->light) * mnubykT;
/*         if(d_tot->omnu->neutrinos_not_analytic && d_tot->ThisTask==0) */
/*             printf("Particle neutrinos start to gravitate NOW: a=%g nufrac_low is: %g\n",a, d_tot->omnu->nufrac_low[0]); */
   }
  /*Precompute factor used to get delta_nu_init. This assumes that delta ~ a, so delta-dot is roughly 1.*/
  deriv_prefac = d_tot->TimeTransfer*(hubble_function(d_tot->TimeTransfer)/d_tot->light)* d_tot->TimeTransfer*mnubykT;
  for (ik = 0; ik < d_tot->nk; ik++) {
      /* Initial condition piece, assuming linear evolution of delta with a up to startup redshift */
      /* This assumes that delta ~ a, so delta-dot is roughly 1. */
      /* Also ignores any difference in the transfer functions between species.
       * This will be good if all species have similar masses, or
       * if two species are massless.
       * Also, since at early times the clustering is tiny, it is very unlikely to matter.*/
      delta_nu_curr[ik] = specialJ(wavenum[ik]*fsl_A0a, qc)*d_tot->delta_nu_init[ik] *(1.+ deriv_prefac*fsl_A0a);
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
        if(Na > 2) {
                params.spline=gsl_interp_alloc(gsl_interp_cspline,Na);
        }
        /*Unless we have only two points*/
        else {
                params.spline=gsl_interp_alloc(gsl_interp_linear,Na);
        }
        params.scale=d_tot->scalefact;
        params.mnubykT=mnubykT;
        params.qc = qc;
        /* Massively over-sample the free-streaming lengths.
         * Interpolation is least accurate where the free-streaming length -> 0,
         * which is exactly where it doesn't matter, but
         * we still want to be safe. */
        int Nfs = Na*16;
        params.fs_acc = gsl_interp_accel_alloc();
        params.fs_spline=gsl_interp_alloc(gsl_interp_cspline,Nfs);

        /*Pre-compute the free-streaming lengths, which are scale-independent*/
        double * fslengths = mymalloc("fslengths", Nfs* sizeof(double));
        double * fsscales = mymalloc("fsscales", Nfs* sizeof(double));
        for(ik=0; ik < Nfs; ik++) {
            fsscales[ik] = log(d_tot->TimeTransfer) + ik*(log(a) - log(d_tot->TimeTransfer))/(Nfs-1.);
            fslengths[ik] = fslength(fsscales[ik], log(a),mnubykT, d_tot->light);
        }
        params.fslengths = fslengths;
        params.fsscales = fsscales;

        if(!params.spline || !params.acc || !w || !params.fs_spline || !params.fs_acc || !fslengths || !fsscales)
              terminate("Error initialising and allocating memory for gsl interpolator and integrator.\n");

        gsl_interp_init(params.fs_spline,params.fsscales,params.fslengths,Nfs);
        for (ik = 0; ik < d_tot->nk; ik++) {
            double abserr,d_nu_tmp;
            params.k=wavenum[ik];
            params.delta_tot=d_tot->delta_tot[ik];
            gsl_interp_init(params.spline,params.scale,params.delta_tot,Na);
            gsl_integration_qag (&F, log(d_tot->TimeTransfer), log(a), 0, 1e-6,GSL_VAL,6,w,&d_nu_tmp, &abserr);
            delta_nu_curr[ik] += d_tot->delta_nu_prefac * d_nu_tmp;
         }
         gsl_integration_workspace_free (w);
         gsl_interp_free(params.spline);
         gsl_interp_accel_free(params.acc);
         myfree(fsscales);
         myfree(fslengths);
   }
   if(d_tot->ThisTask == 0 && d_tot->debug){
          printf("delta_nu_curr[0] is %g\n",delta_nu_curr[0]);
          for(ik=0; ik< 5; ik++)
            printf("k %g d_nu %g\n",wavenum[d_tot->nk/6*ik], delta_nu_curr[d_tot->nk/6*ik]);
   }
   return;
}
