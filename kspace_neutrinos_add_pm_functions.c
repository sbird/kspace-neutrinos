/*This file contains functions which interface closely with the PM grid in gadget.
 add_nu_power_to_rhogrid is the main public function. */

#ifdef KSPACE_NEUTRINOS_2
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>
#include "kspace_neutrinos_func.h"
#include "kspace_neutrino_const.h"
#include "kspace_neutrinos_vars.h"
#include "kspace_neutrinos_private.h"
#include "transfer_init.h"
#include "delta_tot_table.h"

#ifndef MYMPI_COMM_WORLD
#define MYMPI_COMM_WORLD MPI_COMM_WORLD
#endif

void save_nu_power(const double Time, const double logkk[], const double delta_nu[],const int nbins, const int snapnum, const char * OutputDir)
{
    FILE *fd;
    int i;
    char nu_fname[1000];
    /*The last underscore in the filename will be just before the snapshot number.
    * This is daft, but works.*/
    snprintf(nu_fname, 1000,"%s/powerspec_nu_%d", OutputDir, snapnum);
    if(!(fd = fopen(nu_fname, "w"))){
        char buf[1000];
        snprintf(buf, 1000, "can't open file `%s` for writing\n", nu_fname);
        terminate(buf);
    }
    fprintf(fd, "%g\n", Time);
    fprintf(fd, "%d\n", nbins);
    for(i = 0; i < nbins; i++){
        fprintf(fd, "%g %g\n", exp(logkk[i]), delta_nu[i]*delta_nu[i]);
    }
    fclose(fd);
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
 * @param pmgrid is the number of bins in the fourier mesh, PMGRID in allvars.h
 * */
void rebin_power(double logk[],double delta_cdm_curr[],const int nbins, double Kbin[], double SumPower[], long long CountModes[], const int bins_ps, const int pmgrid, const double BoxSize)
{
        int i;
        //Number of modes in a bin: CountModes[1]
        //Total power, without dividing by the number of modes: SumPower[1]
        int istart=0,iend=0;
        /*We will have issues if count overflows a 64 bit integer. */
        unsigned long long count=0;
        /*This is a Fourier conversion factor */
        const double scale=pow(2*M_PI/BoxSize,3);
        double dlogK = (log(Kbin[0]*pmgrid)-log(Kbin[0]))/nbins;
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

/*Helper function for 1D window function.*/
inline fftw_real onedinvwindow(int kx, int n)
{
    //Return \pi x /(n sin(\pi x n)) unless x = 0, in which case return 1.
    return kx ? M_PI*kx/(n*sin(M_PI*kx/(fftw_real)n)) : 1.0;
}

/*The inverse window function of the CiC procedure above. Need to deconvolve this for the power spectrum.
  Only has an effect for k > Nyquist/4.*/
fftw_real invwindow(int kx, int ky, int kz, int n)
{
    if(n == 0)
        return 0;
    fftw_real iwx = onedinvwindow(kx, n);
    fftw_real iwy = onedinvwindow(ky, n);
    fftw_real iwz = onedinvwindow(kz, n);
    return pow(iwx*iwy*iwz,2);
}

/**Little macro to work the storage order of the FFT.*/
#define KVAL(n) ((n)<=dims/2 ? (n) : ((n)-dims))

/* This computes the power in an array on one processor, for an MPI transform.
 * Normalisation and reduction is done in the caller.*/
void total_powerspectrum(const int dims, fftw_complex *outfield, const int nrbins, const int startslab, const int nslab, double *power, long long int *count, double *keffs, const double total_mass)
{
    /*First we sum the power on this processor, then we do an MPI_allgather*/
    double powerpriv[nrbins];
    double keffspriv[nrbins];
    long long int countpriv[nrbins];
    /*How many bins per unit (log) interval in k?*/
    const int binsperunit=nrbins/ceil(log(sqrt(3)*dims/2.0));
    /* Now we compute the powerspectrum in each direction.
     * FFTW is unnormalised, so we need to scale by the length of the array
     * (we do this later). */
    memset(powerpriv, 0, nrbins*sizeof(double));
    memset(countpriv, 0, nrbins*sizeof(long long int));
    memset(keffspriv, 0, nrbins*sizeof(double));
    /* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
     * Use the symmetry of the real fourier transform to half the final dimension.*/
    for(int i=startslab; i<startslab+nslab;i++){
        int indx=(i-startslab)*dims*(dims/2+1);
        for(int j=0; j<dims; j++){
            int indy=j*(dims/2+1);
            /* The k=0 and N/2 mode need special treatment here, 
                * as they alone are not doubled.*/
            /*Do k=0 mode.*/
            int index=indx+indy;
            double kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2));
            //We don't want the 0,0,0 mode as that is just the mean of the field.
            if (kk > 0) {
                int psindex=floor(binsperunit*log(kk));
                powerpriv[psindex] += (outfield[index].re*outfield[index].re+outfield[index].im*outfield[index].im)*pow(invwindow(KVAL(i),KVAL(j),0,dims),2);
                keffspriv[psindex]+=kk;
                countpriv[psindex]++;
            }
            /*Now do the k=N/2 mode*/
            index=indx+indy+dims/2;
            kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(dims/2),2));
            int psindex=floor(binsperunit*log(kk));
            powerpriv[psindex] += (outfield[index].re*outfield[index].re+outfield[index].im*outfield[index].im)*pow(invwindow(KVAL(i),KVAL(j),KVAL(dims/2),dims),2);
            keffspriv[psindex]+=kk;
            countpriv[psindex]++;
            /*Now do the rest. Because of the symmetry, each mode counts twice.*/
            for(int k=1; k<dims/2; k++){
                    index=indx+indy+k;
                    kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(k),2));
                    int psindex=floor(binsperunit*log(kk));
                    powerpriv[psindex]+=2*(outfield[index].re*outfield[index].re+outfield[index].im*outfield[index].im)*pow(invwindow(KVAL(i),KVAL(j),KVAL(k),dims),2);
                    countpriv[psindex]+=2;
                    keffspriv[psindex]+=2*kk;
            }
        }
    }
    /*Now sum the different contributions*/
    MPI_Allgather(countpriv, nrbins * sizeof(long long), MPI_BYTE, count, nrbins * sizeof(long long), MPI_BYTE, MYMPI_COMM_WORLD);
    MPI_Allgather(powerpriv, nrbins * sizeof(double), MPI_BYTE, power, nrbins * sizeof(double), MPI_BYTE, MYMPI_COMM_WORLD);
    MPI_Allgather(keffspriv, nrbins * sizeof(double), MPI_BYTE, keffs, nrbins * sizeof(double), MPI_BYTE, MYMPI_COMM_WORLD);
    /*Normalise by the total mass in the array*/
    for(int i=0; i<nrbins;i++) {
        power[i]/=total_mass*total_mass;
        if(count[i])
            keffs[i]/=count[i];
    }
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


static _transfer_init_table transfer_init;

static _delta_tot_table delta_tot_table;

static _omega_nu omeganu_table;

void broadcast_transfer_table(_transfer_init_table *t_init, int ThisTask)
{
  MPI_Bcast(&(t_init->NPowerTable), 1,MPI_INT,0,MYMPI_COMM_WORLD);
  /*Allocate the memory unless we are on task 0, in which case it is already allocated*/
  if(ThisTask!=0)
    t_init->logk = (double *) mymalloc("Transfer_functions", 2*t_init->NPowerTable* sizeof(double));
  t_init->T_nu=t_init->logk+t_init->NPowerTable;
  /*Broadcast the arrays*/
  MPI_Bcast(t_init->logk,2*(t_init->NPowerTable),MPI_DOUBLE,0,MYMPI_COMM_WORLD);
}

/** This function loads the initial transfer functions from CAMB transfer files.
 * One processor 0 it reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu_init and friends and has length NPowerTable is then broadcast to all processors.
 * Then, on all processors, it allocates memory for delta_tot_table.*/
void transfer_init_tabulate(const int nk_in, const int ThisTask,const double BoxSize, const double UnitLength_in_cm, const double Omega0)
{
  if(omeganu_table.RhoNuTab[0] == 0)
      init_omega_nu(&omeganu_table, kspace_params.MNu, Omega0);
  /*We only need this for initialising delta_tot later.
   * ThisTask is needed so we only read the transfer functions on task 0, serialising disc access.*/
  if(ThisTask==0)
    allocate_transfer_init_table(&transfer_init, nk_in, BoxSize, UnitLength_in_cm, kspace_params.InputSpectrum_UnitLength_in_cm, kspace_params.OmegaBaryonCAMB, kspace_params.KspaceTransferFunction, &omeganu_table);
  /*Broadcast data to other processors*/
  broadcast_transfer_table(&transfer_init, ThisTask);
  /*Set the private copy of the task in delta_tot_table*/
  delta_tot_table.ThisTask = ThisTask;
  allocate_delta_tot_table(&delta_tot_table, nk_in, kspace_params.TimeTransfer, ThisTask);
  /*Check that if we are restarting from a snapshot, we successfully read a table*/
/*   if(fabs(kspace_vars.TimeBegin - d_tot->TimeTransfer) >1e-4 && (!d_tot->ia)) */
/*      terminate("Transfer function not at the same time as simulation start (are you restarting from a snapshot?) and could not read delta_tot table\n"); */
}


#define TARGETBINS 300              /* Number of bins in the smoothed power spectrum*/

/* This function adds the neutrino power spectrum to the
 * density grid. It calls the gadget power spectrum routines, which output their
 * results in SumPower (the total power in all modes in a bin),
 * SumPowerUncorrected (the total power, minus a correction for the bin width)
 * and CountModes (the number of modes per bin).
 * SumPower[1] is from box scale to grid cell scale,
 * SumPower[0] is the folded power on smaller scales.
 * It also touches fft_of_rhogrid, which is the fourier transformed density grid.
 */
void add_nu_power_to_rhogrid(int save, const double Time, const double BoxSize, fftw_complex *fft_of_rhogrid, const int PMGRID, int ThisTask, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir, const double total_mass)
{
  /*Some of the neutrinos will be relativistic at early times. However, the transfer function for the massless neutrinos 
   * is very similar to the transfer function for the massive neutrinos, so treat them the same*/
  const double OmegaNua3 = OmegaNu(&omeganu_table, Time)*pow(Time,3);
  /*kspace_prefac = M_nu / M_cdm */
  const double kspace_prefac = OmegaNua3/(omeganu_table.Omega0-OmegaNu(&omeganu_table, 1));
  int i,x,y,z;
  /*Calculate the power for kspace neutrinos*/
  /* Interpolation structures for the GSL*/
  gsl_interp_accel *acc_nu = gsl_interp_accel_alloc();
  gsl_interp *spline_nu;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spline_cdm;
  double delta_cdm_curr[TARGETBINS];  /*A rebinned power spectrum, smoothed over several modes.*/
  double logkk[TARGETBINS];      /*log k values for the rebinned power spectrum*/
  double delta_nu_curr[TARGETBINS];  /*The square root of the neutrino power spectrum*/
  /*Total power, before rebinning*/
  double sumpower[PMGRID];
  double keffs[PMGRID];
  long long int count[PMGRID];
  /*We calculate the power spectrum at every timestep
   * because we need it as input to the neutrino power spectrum.
   * This function stores the total power*no. modes.*/
  total_powerspectrum(PMGRID, fft_of_rhogrid, PMGRID, slabstart_y, nslab_y, sumpower, count, keffs, total_mass);
  /*Get delta_cdm_curr , which is P(k)^1/2, and logkk*/
  rebin_power(logkk,delta_cdm_curr,TARGETBINS,keffs,sumpower,count,PMGRID, PMGRID, BoxSize);
  /*This sets up P_nu_curr.*/
  get_delta_nu_update(&delta_tot_table, Time, TARGETBINS, logkk, delta_cdm_curr,  delta_nu_curr, &omeganu_table);
  for(i=0;i<TARGETBINS;i++){
          if(isnan(delta_nu_curr[i]) || delta_nu_curr[i] < 0){
                  char err[300];
                  snprintf(err,300,"delta_nu_curr=%g z=%d SmoothPow=%g logkk=%g\n",delta_nu_curr[i],i,delta_cdm_curr[i],exp(logkk[i]));
                  terminate(err);
          }
  }
  /*Sets up the interpolation for get_neutrino_powerspec*/
  spline_cdm=gsl_interp_alloc(gsl_interp_cspline,TARGETBINS);
  gsl_interp_init(spline_cdm,logkk,delta_cdm_curr,TARGETBINS);
  spline_nu=gsl_interp_alloc(gsl_interp_cspline,TARGETBINS);
  gsl_interp_init(spline_nu,logkk,delta_nu_curr,TARGETBINS);
  /*Add P_nu to fft_of_rhgrid*/
  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
        {
          double kx,ky,kz,k2,smth;
          int ip;
          kx = x > PMGRID/2 ? x-PMGRID : x;
          ky = y > PMGRID/2 ? y-PMGRID : y;
          kz = z > PMGRID/2 ? z-PMGRID : z;

          k2 = kx*kx + ky*ky + kz*kz;
          if(k2 <= 0)
              continue;
          /*Change the units of k to match those of logkk*/
          k2=log(sqrt(k2)*2*M_PI/BoxSize);
          /* Note get_neutrino_powerspec returns delta_nu / P_cdm^1/2, which is dimensionless.
           * We have delta_t = (M_cdm+M_nu)*delta_cdm (1-f_nu + f_nu (delta_nu / delta_cdm)^1/2)
           * which gives the right power spectrum, once we divide by
           * M_cdm +M_nu in powerspec*/
          smth=(1+kspace_prefac * get_neutrino_powerspec(k2,logkk, delta_nu_curr, spline_nu,acc_nu, delta_cdm_curr,spline_cdm,acc,TARGETBINS));
          ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
          fft_of_rhogrid[ip].re *= smth;
          fft_of_rhogrid[ip].im *= smth;
        }
  /*If this is being called to save all particle types, save a file with the neutrino power spectrum as well.*/
  if(save){
            save_nu_power(Time, logkk, delta_nu_curr,TARGETBINS, snapnum, OutputDir);
  }
  gsl_interp_free(spline_nu);
  gsl_interp_accel_free(acc_nu);
  gsl_interp_free(spline_cdm);
  gsl_interp_accel_free(acc);
  return;
}

#endif //KSPACE_NEUTRINOS_2
