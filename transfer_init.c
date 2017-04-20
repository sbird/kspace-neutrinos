#include <stdio.h>
#include <math.h>
#include "transfer_init.h"
#include "gadget_defines.h"

/** This function loads the initial transfer functions from CAMB transfer files.
 * It reads the transfer tables from CAMB into the transfer_init structure.
 * Output stored in T_nu and logk with length NPowerTable.*/
void allocate_transfer_init_table(_transfer_init_table *t_init, const double BoxSize, const double UnitLength_in_cm, const double InputSpectrum_UnitLength_in_cm, const double OmegaNu, const double Omega0, const char * KspaceTransferFunction)
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
        endrun(2019,"Can't read input transfer function in file '%s'\n", KspaceTransferFunction);
    }
    t_init->NPowerTable = 0;
    while(1){
        double k, T_cdm, T_b, dummy, T_nu, T_tot, T_nonu;
        char * ret;
        /* read transfer function file */
        ret=fgets(string,1000,fd);
        /*End of file*/
        if(! ret)
        break;
        /* Skip comments*/
        if(string[0] == '#')
            continue;
        if(sscanf(string, " %lg %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &dummy, &dummy, &T_nu, &T_tot, &T_nonu) == 8){
        if(k > kmin)
            t_init->NPowerTable++;
        }
        else
        break;
    }

    fclose(fd);
    message(1,"Found transfer function, using %d rows. Min k used is %g.\n", t_init->NPowerTable,kmin);
    t_init->logk = (double *) mymalloc("Transfer_functions", 2*t_init->NPowerTable* sizeof(double));
    t_init->T_nu=t_init->logk+t_init->NPowerTable;

    /*Now open the file*/
    if(!(fd = fopen(KspaceTransferFunction, "r"))){
        endrun(2020,"Can't read input transfer function in file '%s'\n", KspaceTransferFunction);
    }
    count=0;
    while(count < t_init->NPowerTable){
    /*T_g stores radiation, T_rnu stores massless/relativistic neutrinos*/
    double k, T_b, T_cdm,T_g, T_rnu, T_nu, T_tot, T_nonu;
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
    if(sscanf(string, " %lg %lg %lg %lg %lg %lg %lg %lg", &k, &T_cdm, &T_b, &T_g, &T_rnu, &T_nu, &T_tot, &T_nonu) == 8){
        if(k > kmin){
            /*Set up the total transfer for all the particle species excluding neutrinos*/
            t_init->T_nu[count]= T_nu/T_nonu;
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
        endrun(2021,"Expected %d rows in  file '%s' but only found %d\n", t_init->NPowerTable, KspaceTransferFunction, count);
    }
    fclose(fd);
    return;
}

void free_transfer_init_table(_transfer_init_table *t_init)
{
    myfree(t_init->logk);
}
