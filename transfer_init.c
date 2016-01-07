#include <stdio.h>
#include <math.h>
#include "transfer_init.h"
#include "gadget_defines.h"


/** This function loads the initial transfer functions from CAMB transfer files.
 * It reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu and logk with length NPowerTable.*/
void allocate_transfer_init_table(_transfer_init_table *t_init, int nk_in, const double BoxSize, const double UnitLength_in_cm, const double InputSpectrum_UnitLength_in_cm, const double OmegaBaryonCAMB, const char * KspaceTransferFunction, _omega_nu * omnu)
{
    FILE *fd;
    int count;
    char string[1000];
    /* We aren't interested in modes on scales larger than twice the boxsize*/
    /*Normally 1000*/
    const double scale=(InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
    const double kmin=M_PI/BoxSize*scale;
    /*Set up the table length with the first file found*/
    if(!(fd = fopen(KspaceTransferFunction, "r"))){
        snprintf(string, 1000, "Can't read input transfer function in file '%s'\n", KspaceTransferFunction);
        terminate(string);
    }
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
    if(!(fd = fopen(KspaceTransferFunction, "r"))){
        sprintf(string, "Can't read input transfer function in file '%s'\n", KspaceTransferFunction);
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
    * For this reason we have the variable OmegaBaryonCAMB*/
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
            T_0tot=((omnu->Omega0-OmegaBaryonCAMB-get_omega_nu(omnu, 1))*T_cdm+OmegaBaryonCAMB*T_b)/(omnu->Omega0-get_omega_nu(omnu, 1));
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
        sprintf(string, "Expected %d rows in  file '%s' but only found %d\n", t_init->NPowerTable, KspaceTransferFunction, count);
        terminate(string);
    }
    fclose(fd);
    return;
}

void free_transfer_init_table(_transfer_init_table *t_init)
{
    myfree(t_init->logk);
}
