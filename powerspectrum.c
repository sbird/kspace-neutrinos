#include <math.h>
#include <string.h>
#include <gsl/gsl_interp.h>
#include "powerspectrum.h"
#include "gadget_defines.h"

/*Helper function for 1D window function.*/
static inline fftw_real onedinvwindow(int kx, int n)
{
    /*Return \pi x /(n sin(\pi x / n)) unless x = 0, in which case return 1.*/
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
 * Returns the number of bins with non-zero mode counts.
 * The power spectrum returned is normalised conventionally, and
 * reordered to omit zero bins.*/
int total_powerspectrum(const int dims, fftw_complex *outfield, const int nrbins, const int startslab, const int nslab, double *power, long long int *count, double *keffs, MPI_Comm MYMPI_COMM_WORLD)
{
    /*First we sum the power on this processor, then we do an MPI_allgather*/
    double powerpriv[nrbins];
    double keffspriv[nrbins];
    long long int countpriv[nrbins];
    /*How many bins per unit (log) interval in k?*/
    const double binsperunit=(nrbins-1)/log(sqrt(3)*dims/2.0);
    int i, nonzero;
    double total_mass2 = 0;
    /* First element of the FFT stores the total mass, on the processor with the first slab.
     * Note this may not be the rank 0 processor! */
    if(startslab == 0){
        total_mass2 = outfield[0].re*outfield[0].re + outfield[0].im*outfield[0].im;
    }
    /* Now we compute the powerspectrum in each direction.
     * FFTW is unnormalised, so we need to scale by the length of the array
     * (we do this later). */
    memset(powerpriv, 0, nrbins*sizeof(double));
    memset(countpriv, 0, nrbins*sizeof(long long int));
    memset(keffspriv, 0, nrbins*sizeof(double));
    /* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
     * Use the symmetry of the real fourier transform to half the final dimension.*/
    for(i=startslab; i<startslab+nslab;i++){
        int j, indx=(i-startslab)*dims*(dims/2+1);
        for(j=0; j<dims; j++){
            int k, indy=j*(dims/2+1);
            /* The k=0 and N/2 mode need special treatment here,
                * as they alone are not doubled.*/
            /*Do k=0 mode.*/
            int index=indx+indy;
            double kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2));
            /*We don't want the 0,0,0 mode as that is just the mean of the field.*/
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
            for(k=1; k<dims/2; k++){
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
    MPI_Allreduce(countpriv, count, nrbins, MPI_LONG_LONG_INT, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(powerpriv, power, nrbins, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    MPI_Allreduce(keffspriv, keffs, nrbins, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    /*Make sure total_mass is same on all processors.*/
    MPI_Allreduce(MPI_IN_PLACE, &total_mass2, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
    message(0,"Total powerspectrum mass: %g\n", sqrt(total_mass2));
    /*Normalise by the total mass in the array*/
    for(i=0; i<nrbins;i++) {
        power[i]/=total_mass2;
        if(count[i]) {
            keffs[i]/=count[i];
            power[i] /= count[i];
        }
    }
    /*Remove bins with zero modes*/
    for(i=0, nonzero=0; i<nrbins;i++) {
        if(count[i]) {
            if(nonzero < i) {
                power[nonzero] = power[i];
                keffs[nonzero] = keffs[i];
                count[nonzero] = count[i];
            }
            nonzero++;
        }
    }
    return nonzero;
}
