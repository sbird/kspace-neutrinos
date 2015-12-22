/*This file contains functions which interface closely with the PM grid in gadget.
 add_nu_power_to_rhogrid is the main public function. */

#ifdef KSPACE_NEUTRINOS_2

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include "kspace_neutrinos_func.h"

/*We only need this for fftw_complex*/
#ifdef NOTYPEPREFIX_FFTW
#include        <fftw.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <dfftw.h>	/* double precision FFTW */
#else
#include     <sfftw.h>
#endif
#endif

//Forward define the Efstathiou power spectrum, but we want to remove this later
double PowerSpec_Efstathiou(double k);
/*Forward declaration of powerspec*/
void powerspec(int flag, int *typeflag);

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
void add_nu_power_to_rhogrid(int save, const double Time, const double Omega0, const double BoxSize, double Kbin[], double SumPowerUncorrected[], long long int CountModes[], const int BINS_PS, fftw_complex *fft_of_rhogrid, const int PMGRID, int ThisTask, int slabstart_y, int nslab_y, const int snapnum, const char * OutputDir)
{
  /*Some of the neutrinos will be relativistic at early times. However, the transfer function for the massless neutrinos 
   * is very similar to the transfer function for the massive neutrinos, so treat them the same*/
  const double OmegaNua3 = OmegaNu(Time)*pow(Time,3);
  /*kspace_prefac = M_nu / M_cdm */
  const double kspace_prefac = OmegaNua3/(Omega0-OmegaNu(1));
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
  /*We calculate the power spectrum at every timestep
   * because we need it as input to the neutrino power spectrum.
   * Use the upper bits of flag to not save the potential.*/
  int typelist2[6]={1,1,1,1,1,1};
  /*Zero the power spectra now, because it only gets zerod inside powerspec
   * during the foldonitself calculation, which is only called when mode=1*/
  memset(SumPowerUncorrected, 0, BINS_PS*sizeof(double));
  memset(CountModes, 0, BINS_PS*sizeof(long long int));
  /*Get SumPower, which is P(k)*/
  powerspec(1+ (2<<29),typelist2);
  /*Get delta_cdm_curr , which is P(k)^1/2, and logkk*/
  rebin_power(logkk,delta_cdm_curr,TARGETBINS,Kbin,SumPowerUncorrected,CountModes,BINS_PS, PMGRID, BoxSize);
  /*Now we don't need SumPower[1] anymore: zero it again so it doesn't mess up state.*/
  memset(SumPowerUncorrected, 0, BINS_PS*sizeof(double));
  memset(CountModes, 0, BINS_PS*sizeof(long int));
  /*This sets up P_nu_curr.*/
  get_delta_nu_update(Time, TARGETBINS, logkk, delta_cdm_curr,  delta_nu_curr, ThisTask);
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

          k2 = kx * kx + ky * ky + kz * kz;
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
