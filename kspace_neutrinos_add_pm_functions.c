#ifdef KSPACE_NEUTRINOS_2

#include <gsl/gsl_interp.h>

void save_nu_power(double logkk[], double delta_nu[],int nbins)
{
  if(ThisTask == 0)
    {
      FILE *fd;
      char buf[500];
      int i;
      char nu_fname[500];
      /*The last underscore in the filename will be just before the snapshot number.
       * This is daft, but works.*/
      char* under=strrchr(power_spec_fname,'_');
      if(under)
        snprintf(nu_fname, 500,"%s/powerspec_nu_%s", All.OutputDir, under+1);
      else
        snprintf(nu_fname, 500,"%s/powerspec_nu_XXX.txt", All.OutputDir);

      if(!(fd = fopen(nu_fname, "w"))){
	  sprintf(buf, "can't open file `%s`\n", power_spec_fname);
	  terminate(buf);
      }
      fprintf(fd, "%g\n", All.Time);
      fprintf(fd, "%d\n", nbins);
      for(i = 0; i < nbins; i++){
          fprintf(fd, "%g %g\n", exp(logkk[i]), delta_nu[i]*delta_nu[i]);
      }
      fclose(fd);

    }
    return;
}

/* This function adds the neutrino power spectrum to the
 * density grid. It calls the gadget power spectrum routines, which output their
 * results in SumPower (the total power in all modes in a bin),
 * SumPowerUncorrected (the total power, minus a correction for the bin width)
 * and CountModes (the number of modes per bin).
 * SumPower[1] is from box scale to grid cell scale,
 * SumPower[0] is the folded power on smaller scales.
 * It also touches fft_of_rhogrid, which is the fourier transformed density grid.
 */
void add_nu_power_to_rhogrid(int save)
{
  /*Some of the neutrinos will be relativistic at early times. However, the transfer function for the massless neutrinos 
   * is very similar to the transfer function for the massive neutrinos, so treat them the same*/
  const double OmegaNua3 = OmegaNu(All.Time)*pow(All.Time,3);
  /*kspace_prefac = M_nu / M_cdm */
  const double kspace_prefac = OmegaNua3/(All.Omega0-All.OmegaNu);
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
  for(i = 0; i < BINS_PS; i++)
  {
      SumPower[1][i] = 0;
      SumPowerUncorrected[1][i] = 0;
      CountModes[1][i] = 0;
  }
  /*Get SumPower, which is P(k)*/
  powerspec(1+ (2<<29),typelist2);
  /*Get delta_cdm_curr , which is P(k)^1/2, and logkk*/
  rebin_power(logkk,delta_cdm_curr,TARGETBINS,Kbin,SumPower[1],CountModes[1],BINS_PS);
  /*Now we don't need SumPower[1] anymore: zero it again so it doesn't mess up state.*/
  for(i = 0; i < BINS_PS; i++)
  {
      SumPower[1][i] = 0;
      SumPowerUncorrected[1][i] = 0;
      CountModes[1][i] = 0;
  }
  /*This sets up P_nu_curr.*/
  get_delta_nu_update(All.Time, TARGETBINS, logkk, delta_cdm_curr,  delta_nu_curr);
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
          k2=log(sqrt(k2)*2*M_PI/All.BoxSize);
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
            save_nu_power(logkk, delta_nu_curr,TARGETBINS);
  }
  gsl_interp_free(spline_nu);
  gsl_interp_accel_free(acc_nu);
  gsl_interp_free(spline_cdm);
  gsl_interp_accel_free(acc);
  return;
}

#endif //KSPACE_NEUTRINOS_2
