#include <stdio.h>
#include <math.h>
#include "transfer_init.h"
#include "kspace_neutrinos_vars.h"


/** This function loads the initial transfer functions from CAMB transfer files.
 * It reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu and logk with length NPowerTable.*/
void allocate_transfer_init_table(_transfer_init_table *t_init, int nk_in)
{
    FILE *fd;
    int count;
    char string[1000];
    /* We aren't interested in modes on scales larger than twice the boxsize*/
    /*Normally 1000*/
    const double scale=(kspace_params.InputSpectrum_UnitLength_in_cm / kspace_vars.UnitLength_in_cm);
    const double kmin=M_PI/kspace_vars.BoxSize*scale;
    /*Set up the table length with the first file found*/
    if(!(fd = fopen(kspace_params.KspaceTransferFunction, "r"))){
        snprintf(string, 1000, "Can't read input transfer function in file '%s'\n", kspace_params.KspaceTransferFunction);
        terminate(string);
    }
    if(kspace_params.TimeTransfer > kspace_vars.TimeBegin + 1e-4){
        snprintf(string, 1000,"Transfer function is at a=%g but you tried to start the simulation earlier, at a=%g\n", kspace_params.TimeTransfer, kspace_vars.TimeBegin);
        terminate(string);
    }

    t_init->TimeTransfer = kspace_params.TimeTransfer;
    t_init->NPowerTable = 0;
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
            t_init->NPowerTable++;
        }
        else
        break;
    }

    fclose(fd);
    printf("Found transfer function, using %d rows. Min k used is %g.\n", t_init->NPowerTable,kmin);
    t_init->logk = (double *) mymalloc("Transfer_functions", 2*t_init->NPowerTable* sizeof(double));
    t_init->T_nu=t_init->logk+t_init->NPowerTable;

    /*Now open the file*/
    if(!(fd = fopen(kspace_params.KspaceTransferFunction, "r"))){
        sprintf(string, "Can't read input transfer function in file '%s'\n", kspace_params.KspaceTransferFunction);
        terminate(string);
    }
    count=0;
    while(count < t_init->NPowerTable){
    /*T_g stores radiation, T_rnu stores massless/relativistic neutrinos*/
    double k, T_b, T_cdm,T_g, T_rnu, T_nu, T_tot;
    /*T_0tot stores all the species for which there are particles*/
    double T_0tot;
    /*We may have "faked" baryons by including them into the DM in Gadget.
    * In this case All.OmegaBaryon will be zero, but CAMB will have been fed a different value.
    * For this reason we have the variable kspace_params.OmegaBaryonCAMB*/
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
            T_0tot=((kspace_vars.Omega0-kspace_params.OmegaBaryonCAMB-kspace_vars.OmegaNu)*T_cdm+kspace_params.OmegaBaryonCAMB*T_b)/(kspace_vars.Omega0-kspace_vars.OmegaNu);
            t_init->T_nu[count]= T_nu/T_0tot;
            /*k has units of 1/Mpc, need 1/kpc */
            k /= scale; /* Convert to internal units*/
            t_init->logk[count] = log(k);
            count++;
        }
    }
    else
        break;
    }
    if(count < t_init->NPowerTable){
        sprintf(string, "Expected %d rows in  file '%s' but only found %d\n", t_init->NPowerTable, kspace_params.KspaceTransferFunction, count);
        terminate(string);
    }
    fclose(fd);
    return;
}
