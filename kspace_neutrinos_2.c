/*This file contains calculations for the Fourier-space semi-linear neutrino method
 * described in Ali-Haimoud and Bird 2012.
 * It has two parts: the first deals with calculating the background evolution of the neutrinos,
 * by integrating their energy density numerically.
 * The second contains routines for finding the neutrino power spectrum from a (non-linear) CDM power spectrum.
 */
#include "allvars.h"
#include "proto.h"

#ifdef KSPACE_NEUTRINOS_2
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <sys/types.h>
#include <glob.h>

#ifdef NEUTRINOS
#error "Cannot define particle based and Fourier-space neutrinos at the same time!"
#endif

#ifdef KSPACE_NEUTRINOS
#error "KSPACE_NEUTRINOS_2 is incompatible with KSPACE_NEUTRINOS"
#endif

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
/* Note Omega0, the total non-relativistic matter, includes neutrinos (and radiation). */
#define FLOAT   1e-6            /*Floating point accuracy*/
#define BOLEVK 8.61734e-5        /*The Boltzmann constant in units of eV/K*/
#define H0      (All.HubbleParam)                /* H0 in units of H100*/
/*The time at which we first start our integrator:
 * NOTE! This is not All.TimeBegin, but the time of the transfer function file,
 * so that we can support restarting from snapshots.*/
#define A0      (All.TimeTransfer)
/*Light speed in internal units. C is defined in allvars.h to be lightspeed in cm/s*/
#define LIGHT (C*All.UnitTime_in_s/All.UnitLength_in_cm)
#define HBAR    6.582119e-16  /*hbar in units of eV s*/

#define H100   All.Hubble /* 100 km/s/Mpc in units of 1/UnitTime. */

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
#define GSL_VAL 200
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
     double logA0=log(A0);
     double logaf=log(af);
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
            if(!(RhoNuTab[sp]))
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
        rhonu /= All.HubbleParam*All.HubbleParam;
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
                 if(fabs(All.MNu[mi] -All.MNu[mmi]) < FLOAT)
                   break;
             }
             if(mmi==mi)
                 OmegaNu_one[mi]=OmegaNu_single(a,All.MNu[mi],mi);
            rhonu+=OmegaNu_one[mmi];
        }
        return rhonu;
}

/*Arrays to store the initial transfer functions from CAMB.
 * We store transfer functions because we want to use the
 * CDM + Baryon total matter power spectrum from the
 * first timestep of Gadget, so that for possible Rayleigh scattering
 * in the initial conditions is included in the neutrino and radiation components. */
static int NPowerTable;

static double *logk_init;
static double *T_nu_init;

/*Now we want to define a static object to store all previous delta_tot.
 * This object needs a constructor, a few private data members, and a way to be read and written from disk.
 * nk is fixed, delta_tot, scalefact and ia are updated in get_delta_nu_update*/

/* Number of k values stored in each power spectrum*/
static int nk;

/*Maximum number of redshifts to store. Redshifts are stored every delta a = 0.01 */
static int namax;
/* Number of already "recorded" time steps, i.e. scalefact[0...ia-1] is recorded.
 * Current time corresponds to index ia (but is only recorded if sufficiently far from previous time).
 * Caution: ia here is different from Na in get_delta_nu (Na = ia+1).*/
static int ia;
/* Pointer to nk arrays of length namax containing the total power spectrum.*/
static double **delta_tot;
/* Array of length namax containing scale factors at which the power spectrum is stored*/
static double *scalefact;
/*Pointer to array of length nk storing initial neutrino power spectrum*/
static double *delta_nu_init;
/*Pointer to array of length nk storing the last neutrino power spectrum we saw, for a first estimate
 * of the new delta_tot */
static double *delta_nu_last;


/** This function loads the initial transfer functions from CAMB transfer files.
 * Note it uses the global parameter All.KspaceTransferFunction for the filename to read.
 * Output stored in T_nu_init and friends and has length NPowerTable. */
void transfer_init_tabulate(int nk_in)
{
  FILE *fd;
  int count;
  char string[1000];
  /* We aren't interested in modes on scales larger than twice the boxsize*/
  /*Normally 1000*/
  const double scale=(All.InputSpectrum_UnitLength_in_cm / All.UnitLength_in_cm);
  const double kmin=M_PI/All.BoxSize*scale;
  /*We only need this for initialising delta_tot later, which is only done on task 0.
   * So only read the transfer functions on that task*/
  if(ThisTask==0){
     /*Set up the table length with the first file found*/
     if(!(fd = fopen(All.KspaceTransferFunction, "r"))){
         sprintf(string, "Can't read input transfer function in file '%s' on task %d\n", All.KspaceTransferFunction, ThisTask);
         terminate(string);
     }
     if(All.TimeTransfer > All.TimeBegin + 1e-4){
         snprintf(string, 1000,"Transfer function is at a=%g but you tried to start the simulation earlier, at a=%g\n", All.TimeTransfer, All.TimeBegin);
         terminate(string);
     }

     NPowerTable = 0;
     while(1){
         double k, T_cdm, T_b, dummy, T_nu, T_tot;
         char * ret;
         /* read transfer function file */
         ret=fgets(string,1000,fd);
         /*End of file*/
         if(! ret)
             break;
         /* Skip comments*/
         if(string[0] == '#')
             continue;
         if(sscanf(string, " %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &dummy, &dummy, &T_nu, &T_tot) == 7){
                 if(k > kmin)
                           NPowerTable++;
         }
         else
           break;
     }

     fclose(fd);
     if(ThisTask == 0)
       printf("Found transfer function, using %d rows. Min k used is %g.\n", NPowerTable,kmin);
     logk_init = (double *) mymalloc("Transfer_functions", 2*NPowerTable* sizeof(double));
     T_nu_init=logk_init+NPowerTable;

     /*Now open the file*/
     if(!(fd = fopen(All.KspaceTransferFunction, "r"))){
         sprintf(string, "Can't read input transfer function in file '%s' on task %d\n", All.KspaceTransferFunction, ThisTask);
         terminate(string);
     }
     count=0;
     while(count < NPowerTable){
       /*T_g stores radiation, T_rnu stores massless/relativistic neutrinos*/
       double k, T_b, T_cdm,T_g, T_rnu, T_nu, T_tot;
       /*T_0tot stores all the species for which there are particles*/
       double T_0tot;
       /*We may have "faked" baryons by including them into the DM in Gadget.
        * In this case All.OmegaBaryon will be zero, but CAMB will have been fed a different value.
        * For this reason we have the variable All.OmegaBaryonCAMB*/
       char * ret;
       /* read transfer function file */
       ret=fgets(string,1000,fd);
       /*End of file*/
       if(!ret)
           break;
       /* Skip comments*/
       if(string[0] == '#')
           continue;
       /* read transfer function file from CAMB */
       if(sscanf(string, " %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &T_g, &T_rnu, &T_nu, &T_tot) == 7){
           if(k > kmin){

             /*Combine the massive and massless neutrinos.*/
             /*Set up the total transfer for all the species with particles*/
             T_0tot=((All.Omega0-All.OmegaBaryonCAMB-All.OmegaNu)*T_cdm+All.OmegaBaryonCAMB*T_b)/(All.Omega0-All.OmegaNu);
             T_nu_init[count]= T_nu/T_0tot;
             /*k has units of 1/Mpc, need 1/kpc */
             k /= scale; /* Convert to internal units*/
             logk_init[count] = log(k);
             count++;
           }
       }
       else
         break;
     }
     if(count < NPowerTable){
         sprintf(string, "Expected %d rows in  file '%s' but only found %d on task %d\n", NPowerTable, All.KspaceTransferFunction, count,ThisTask);
         terminate(string);
     }
     fclose(fd);
   }
   /*Memory allocations need to be done on all processors*/
   nk=nk_in;
   /*Allocate memory for delta_tot here, so that we can have further memory allocated and freed
    * before delta_tot_init is called. The number nk here should be larger than the actual value needed.*/
   /*Allocate pointers to each k-vector*/
   namax=ceil(100*(All.TimeMax-All.TimeTransfer))+2;
   ia=0;
   delta_tot =(double **) mymalloc("kspace_delta_tot",nk*sizeof(double *));
   /*Allocate list of scale factors, and space for delta_tot, in one operation.*/
   scalefact = (double *) mymalloc("kspace_scalefact",namax*(nk+1)*sizeof(double));
   /*Allocate actual data. Note that this means data can be accessed either as:
    * delta_tot[k][a] OR as
    * delta_tot[0][a+k*namax] */
   delta_tot[0] = scalefact+namax;
   for(count=1; count< nk; count++)
        delta_tot[count] = delta_tot[0] + count*namax;
   /*Allocate space for the initial neutrino power spectrum*/
   delta_nu_init =(double *) mymalloc("kspace_delta_nu_init",2*nk*sizeof(double));
   delta_nu_last=delta_nu_init+nk;
   return;
}

#define MINMODES 1
/*This function rebins Power[1] to have a constant lower number of bins
 *
 * @param logk stores the output wavenumbers
 * @param delta_cdm_curr is the square root of the rebinned power
 * @param nbins is the length of the output vectors
 * @param Kbin is the wavenumbers of the input arrays
 * @param SumPower is the input power, calculated in powerspec in pm_periodic
 * @param CountModes is the number of modes in each bin in SumPower
 * @param bins_ps the number of bins in SumPower. Defined as BINS_PS in pm_periodic.
 * */
void rebin_power(double logk[],double delta_cdm_curr[],int nbins, double Kbin[], double SumPower[], long long CountModes[], int bins_ps)
{
        int i;
        //Number of modes in a bin: CountModes[1]
        //Total power, without dividing by the number of modes: SumPower[1]
        int istart=0,iend=0;
        /*We will have issues if count overflows a 64 bit integer. */
        unsigned long long count=0;
        /*This is a Fourier conversion factor */
        const double scale=pow(2*M_PI/All.BoxSize,3);
        double dlogK = (log(Kbin[0]*PMGRID)-log(Kbin[0]))/nbins;
        double logK_A[bins_ps];
        int MaxIndex=0;
        double SlogK[nbins];
        double Spk[nbins];
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_interp *spline;
        for(i=0; i<bins_ps; i++)
                logK_A[i] = log(Kbin[i]);
        /*First smooth the array so that there are enough modes in each bin.
         * These bins will be uneven, and their number will be uncertain, but less than nbins*/
        while(iend < bins_ps){
                count+=CountModes[iend];
                iend++;
                if (count >= MINMODES && (logK_A[iend-1]-logK_A[istart] >= dlogK)){
                        double pk = 0,kk=0;
                        for(i=istart;i<iend;i++){
                                pk+=SumPower[i];
                                kk+=logK_A[i]* CountModes[i];
                        }
                        pk/=(count*scale);
                        kk/=count;
                        //This is a correction from what is done in powerspec
                        pk *= PowerSpec_Efstathiou(exp(kk));
                        SlogK[MaxIndex] = kk;
                        Spk[MaxIndex]=pk;
                        MaxIndex++;
                        istart=iend;
                        count=0;
                        if(MaxIndex >= nbins)
                                break;
                }
        }
        /*So now we have an array which is not zero, but which has an uncertain binning.
         * Rebin it to a regular grid by interpolation, and find the exponent. */
        for(i=0; i<nbins; i++)
                logK_A[i] = log(Kbin[0])+dlogK*i-0.01;
        /*Final term is just to guard against floating point roundoff*/
        /*Allocate an interpolating spline to rebin this onto a regular grid.
         * Use linear spline because the power spectrum will be jagged due to Rayleigh scattering.*/
        spline=gsl_interp_alloc(gsl_interp_linear,MaxIndex);
        if(spline == NULL || acc == NULL || gsl_interp_init(spline,SlogK,Spk,MaxIndex))
                terminate("Error initialising and allocating memory for gsl interpolator.\n");

        for(i=0; i< nbins; i++){
                double x=logK_A[i];
                if(x < SlogK[0])
                        x=SlogK[0];
                /*This will just be floating point roundoff*/
                if(x > SlogK[MaxIndex-1])
                        x=SlogK[MaxIndex-1];
                logk[i]=logK_A[i];
                delta_cdm_curr[i]=gsl_interp_eval(spline,SlogK,Spk,x,acc);
                /*Guard against floating point negatives and nans*/
                if (delta_cdm_curr[i] < 0)
                    delta_cdm_curr[i] = 0;
                delta_cdm_curr[i] = sqrt(delta_cdm_curr[i]);
        }
        gsl_interp_free(spline);
        gsl_interp_accel_free(acc);

        return;
}

#if defined(PMGRID)
double fslength(double ai, double af,double mnu);
double find_ai(double af, double k, double kfsl,double mnu);
double specialJ_fit(double x);
void get_delta_nu(double a, int Na, double wavenum[], double delta_nu_curr[],const double Omega0,double mnu);

void save_delta_tot(ia)
{
   if(ThisTask==0){
        FILE *fd;
        int i;
        if(!(fd = fopen("delta_tot_nu.txt", "a")))
             terminate("Could not open delta_tot_nu.txt for writing!\n");
        /*Write log scale factor*/
        fprintf(fd, "# %le ", scalefact[ia]);
        /*Write kvalues*/
        for(i=0;i<nk; i++)
                fprintf(fd,"%le ",delta_tot[i][ia]);
        fprintf(fd,"\n");
        fclose(fd);
   }
   return;
}

/*Returns kT / a M_nu (which is dimensionless) in the relativistic limit
 * where it is kT / (a^2 m_nu^2 + (kT)^2)^(1/2)
 * Takes a* M_nu as argument. */
inline double kTbyaM(double amnu)
{
  const double kT=BOLEVK*TNU;
  return kT/amnu;
}

/* Constructor. transfer_init_tabulate must be called before this function.
 * Initialises delta_tot (including from a file) and delta_nu_init from the transfer functions.
 * Note delta_cdm_curr includes baryons*/
void delta_tot_init(int nk_in, double wavenum[], double delta_cdm_curr[])
{
   if(nk_in > nk){
           char err[500];
           sprintf(err,"input power of %d is longer than memory of %d\n",nk_in,nk);
           terminate(err);
   }
   nk=nk_in;
    /*Construct delta_nu_init only if we are on the first processor.
     * This avoids races with the file access*/
   if(ThisTask==0){
        int ik;
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_interp *spline;
        FILE* fd;
        double OmegaNua3=OmegaNu(All.TimeTransfer)*pow(All.TimeTransfer,3);
        /*Load delta_tot from a file, if such a file exists. Allows resuming.*/
        if((fd = fopen("delta_tot_nu.txt", "r"))){
             /*Read redshifts; Initial one is known already*/
             int iia;
             for(iia=0; iia< namax;iia++){
                     double scale;
                     if(fscanf(fd, "# %lg ", &scale) != 1)
                             break;
                     scalefact[iia]=scale;
                     /*Only read until we reach the present day*/
                     if(log(All.Time) <= scale)
                             break;
                     /*Read kvalues*/
                     for(ik=0;ik<nk; ik++)
                             /*If we do not have a complete delta_tot for one redshift, we stop*/
                             if(fscanf(fd, "%lg ", &(delta_tot[ik][iia])) != 1){
                                     char err[150];
                                     snprintf(err,150,"Incomplete delta_tot in delta_tot_nu.txt; a=%g\n",exp(scalefact[iia]));
                                     terminate(err);
                             }
             }
             /*If our table starts at a different time from the simulation, stop.*/
             if(fabs(scalefact[0] - log(All.TimeTransfer)) > 1e-4){
                     char err[250];
                     snprintf(err,250," delta_tot_nu.txt starts wih a=%g, transfer function is at a=%g\n",exp(scalefact[0]),All.TimeTransfer);
                     terminate(err);
             }

             if(iia > 0)
                     ia=iia;
             printf("Read %d stored power spectra from delta_tot_nu.txt\n",iia);
             fclose(fd);
             /*Get a clean restart file*/
#ifndef NOCALLSOFSYSTEM
             system("mv delta_tot_nu.txt delta_tot_nu.txt.bak");
#else
             fd=fopen("delta_tot_nu.txt","w");
             fclose(fd);
#endif
        }
        /*Check that if we are restarting from a snapshot, we successfully read a table*/
        if(fabs(All.TimeBegin - All.TimeTransfer) >1e-4 && (!ia))
            terminate("Transfer function not at the same time as simulation start (are you restarting from a snapshot?) and could not read delta_tot table\n");
        /*Initialise the first delta_tot to use the first timestep's delta_cdm_curr
         * so that it includes potential Rayleigh scattering. */
        scalefact[0]=log(All.TimeTransfer);
        spline=gsl_interp_alloc(gsl_interp_cspline,NPowerTable);
        gsl_interp_init(spline,logk_init,T_nu_init,NPowerTable);
        for(ik=0;ik<nk;ik++){
               double T_nubyT_0 = gsl_interp_eval(spline,logk_init,T_nu_init,log(wavenum[ik]),acc);
               /*The total power spectrum using neutrinos and radiation from the CAMB transfer functions:
                * The CAMB transfer functions are defined such that
                * P_cdm ~ T_cdm^2 (and some other constant factors)
                * then P_t = (Omega_cdm P_cdm + Omega_nu P_nu)/(Omega_cdm + Omega_nu)
                *          = P_cdm (Omega_cdm+ Omega_nu (P_nu/P_cdm)) / (Omega_cdm +Omega_nu)
                *          = P_cdm (Omega_cdm+ Omega_nu (T_nu/T_cdm)^2) / (Omega_cdm+Omega_nu) */
               double CDMtoTot=((All.Omega0-All.OmegaNu)+pow(T_nubyT_0,2)*OmegaNua3)/(All.Omega0-All.OmegaNu+OmegaNua3);
               /*We only want to use delta_cdm_curr if we are not restarting*/
               if(!ia)
                     delta_tot[ik][0] = delta_cdm_curr[ik]*sqrt(CDMtoTot);
               /* Also initialise delta_nu_init here to save time later.
                * Use the first delta_tot, in case we are resuming.*/
               delta_nu_init[ik] = delta_tot[ik][0]/sqrt(CDMtoTot)*fabs(T_nubyT_0);
        }
        gsl_interp_accel_free(acc);
        gsl_interp_free(spline);
        if(!ia)
             ia=1;
        for(ik=0; ik< ia; ik++)
             save_delta_tot(ik);
   }
   /*Broadcast data to other processors*/
   MPI_Bcast(&ia,1,MPI_INT,0,MPI_COMM_WORLD);
   /*scalefact is a chunk of memory that also includes the delta_tot values*/
   MPI_Bcast(scalefact,namax*(1+nk),MPI_DOUBLE,0,MPI_COMM_WORLD);
   /*Finally delta_nu_init*/
   MPI_Bcast(delta_nu_init,nk,MPI_DOUBLE,0,MPI_COMM_WORLD);
   return;

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
***************************************************************************************************/

double specialJ_fit(double x)
{
  double x2, x4, x8;
  if (x <= 0.) return 1.;
  x2 = x*x;
  x4 = x2*x2;
  x8 = x4*x4;

  return (1.+ 0.0168 * x2 + 0.0407* x4)/(1. + 2.1734 * x2 + 1.6787 * exp(4.1811*log(x)) +  0.1467 * x8);
}

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
    return fsl_aia/(ai*hubble_function(ai)) *specialJ_fit(p->k*fsl_aia) * delta_tot_at_a/(ai*kTbyaM(ai*p->mnu));
}

/*****************************************************************************************************
Main function: given tables of wavenumbers, total delta at Na earlier times (<= a),
and initial conditions for neutrinos, computes the current delta_nu.
Na is the number of currently stored time steps.
Requires transfer_init_tabulate to have been called prior to first call.
******************************************************************************************************/

void get_delta_nu(double a, int Na, double wavenum[],  double delta_nu_curr[],const double Omega0,double mnu)
{
  double fsl_A0a,deriv_prefac;
  int ik;
  if(ThisTask == 0)
          printf("Start get_delta_nu: a=%g Na =%d wavenum[0]=%g delta_tot[0]=%g m_nu=%g\n",a,Na,wavenum[0],delta_tot[0][Na-1],mnu);

  fsl_A0a = fslength(A0, a,mnu);
  /*Precompute factor used to get delta_nu_init. This assumes that delta ~ a, so delta-dot is roughly 1.*/
  deriv_prefac = A0*(hubble_function(A0)/LIGHT)/kTbyaM(A0*mnu);

  for (ik = 0; ik < nk; ik++) {
      /* Initial condition piece, assuming linear evolution of delta with a up to startup redshift */
      /* This assumes that delta ~ a, so delta-dot is roughly 1. */
      /* Also ignores any difference in the transfer functions between species.
       * This will be good if all species have similar masses, or
       * if two species are massless.
       * Also, since at early times the clustering is tiny, it is very unlikely to matter.*/
      delta_nu_curr[ik] = specialJ_fit(wavenum[ik]*fsl_A0a)*delta_nu_init[ik] *(1.+ deriv_prefac*fsl_A0a);
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
        params.scale=scalefact;
        params.mnu=mnu;
        for (ik = 0; ik < nk; ik++) {
            double abserr,d_nu_tmp;
            params.k=wavenum[ik];
            params.delta_tot=delta_tot[ik];
            gsl_interp_init(params.spline,params.scale,params.delta_tot,Na);
            gsl_integration_qag (&F, log(A0), log(a), 0, 1e-6,GSL_VAL,6,w,&d_nu_tmp, &abserr);
            delta_nu_curr[ik] += 1.5 *Omega0 *H100*H100/LIGHT*d_nu_tmp;
         }
         gsl_integration_workspace_free (w);
         gsl_interp_free(params.spline);
         gsl_interp_accel_free(params.acc);
   }
  if(ThisTask == 0){
          printf("delta_nu_curr[0] is %g\n",delta_nu_curr[0]);
          for(ik=0; ik< 5; ik++)
          printf("k %g d_nu %g\n",wavenum[40*ik], delta_nu_curr[40*ik]);
  }
   return;
}

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value isfor all neutrino species*/
void get_delta_nu_combined(double a, int Na, double wavenum[],  double delta_nu_curr[],const double Omega0)
{
    double Omega_nu_tot=OmegaNu(a);
    int mi;
    double delta_nu_single[NUSPECIES][nk];
    /*Initialise delta_nu_curr*/
    for(mi=0;mi<nk;mi++)
        delta_nu_curr[mi]=0;
    /*Get each neutrinos species and density separately and add them to the total.
     * Neglect perturbations in massless neutrinos.*/
    for(mi=0;mi< NUSPECIES; mi++){
         if(All.MNu[mi] > 0){
             int ik,mmi;
             double OmegaNu_frac=OmegaNu_single(a,All.MNu[mi],mi)/Omega_nu_tot;
             for(mmi=0; mmi<mi; mmi++){
                 if(fabs(All.MNu[mi] -All.MNu[mmi]) < FLOAT)
                   break;
             }
             if(mmi==mi)
                get_delta_nu(a, Na, wavenum, delta_nu_single[mi],Omega0,All.MNu[mi]);
             for(ik=0; ik<nk; ik++)
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
void get_delta_nu_update(double a, int nk_in, double logk[], double delta_cdm_curr[], double delta_nu_curr[])
{
  int ik;
  double wavenum[nk_in];
  const double OmegaNua3=OmegaNu(a)*pow(a,3);
  const double Omega0 = All.Omega0 - All.OmegaNu + OmegaNua3;
  const double fnu = OmegaNua3/Omega0;
  for (ik = 0; ik < nk_in; ik++)
           wavenum[ik]=exp(logk[ik]);
  /* Get a delta_nu_curr from CAMB.*/
  if(!ia){
       /*Initialise delta_tot, setting ia = 1
        * and signifying that we are ready to leave the relativistic regime*/
       delta_tot_init(nk_in, wavenum,delta_cdm_curr);
       /*Initialise delta_nu_last*/
       get_delta_nu_combined(a, ia, wavenum, delta_nu_last,Omega0);
  }

  /*If we get called twice with the same scale factor, do nothing*/
  if(log(a)-scalefact[ia-1] < FLOAT){
       for (ik = 0; ik < nk; ik++)
               delta_nu_curr[ik] = delta_nu_last[ik];
       return;
  }

  /* We need some estimate for delta_tot(current time) to obtain delta_nu(current time).
     Even though delta_tot(current time) is not directly used (the integrand vanishes at a = a(current)),
     it is indeed needed for interpolation */
  /* It was checked that using delta_tot(current time) = delta_cdm(current time) leads to no more than 2%
     error on delta_nu (and moreover for large k). A simple estimate for delta_nu decreases the maximum
     relative error on delta_nu to ~1E-4. So we only need one step. */
   scalefact[ia] = log(a);
   /* First estimate for delta_tot(a) */
   for (ik = 0; ik < nk; ik++) {
      /*Guard against floating point zeros*/
      delta_tot[ik][ia] = (1.-fnu)*delta_cdm_curr[ik]+fnu*delta_nu_last[ik];
   }
   /*Get the new delta_nu_curr*/
   get_delta_nu_combined(a, ia+1, wavenum, delta_nu_curr,Omega0);
   /*Update delta_nu_last*/
   for (ik = 0; ik < nk; ik++)
       delta_nu_last[ik]=delta_nu_curr[ik];
   /* Decide whether we save the current time or not */
   if (a > exp(scalefact[ia-1]) + 0.01) {
       /* If so update delta_tot(a) correctly */
       for (ik = 0; ik < nk; ik++){
               delta_tot[ik][ia] = fnu*delta_nu_curr[ik]+(1.-fnu)*delta_cdm_curr[ik];
       }
       ia++;
       if (ThisTask==0){
              printf("Updating delta_tot: a=%f, Na=%d, last=%f\n",a,ia,exp(scalefact[ia-2]));
              save_delta_tot(ia-1);
       }
   }

   return;
}

/*Get the neutrino power spectrum. This will become:
 * \delta_\nu = \delta_{CDM}(\vec{k}) * delta_\nu(k) / delta_{CDM} (k),
 * thus we get the right powerspectrum.
 * @param kk log(k) value to get delta_nu at
 * @param logkk[] vector of k values corresponding to delta_nu_curr
 * @param delta_nu_curr vector of delta_nu
 * @param spline_nu interpolating spline for delta_nu_curr
 * @param acc_nu accelerator for spline_nu
 * @param delta_cdm_curr vector of delta_cdm
 * @param spline_cdm interpolating spline for delta_cdm_curr
 * @param acc accelerator for spline_cdm
 * @param nbins number of bins in delta_nu_curr and delta_cdm_curr
 * @returns delta_nu / delta_CDM
 * */
double get_neutrino_powerspec(double kk, double logkk[], double delta_nu_curr[], gsl_interp *spline_nu, gsl_interp_accel * acc_nu, double delta_cdm_curr[], gsl_interp *spline_cdm, gsl_interp_accel * acc,int nbins)
{
        double delta_cdm,delta_nu;
        int ret;
        if(kk < logkk[0]){
                char err[300];
                snprintf(err,300,"trying to extract a k= %g < min stored = %g \n",kk,logkk[0]);
                terminate(err);
        }
        /*This is just to guard against floating point roundoff*/
        if( kk > logkk[nbins-1])
                kk=logkk[nbins-1];
        delta_cdm=gsl_interp_eval(spline_cdm,logkk, delta_cdm_curr,kk,acc);
        delta_nu=gsl_interp_eval(spline_nu,logkk, delta_nu_curr, kk,acc_nu);
        if(isnan(delta_cdm) || isnan(delta_nu))
                terminate("delta_nu or delta_cdm is nan\n");
        return delta_nu/delta_cdm;
}


#endif //PMGRID

#endif //KSPACE_NEUTRINOS_2
