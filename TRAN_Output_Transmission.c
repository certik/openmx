#include <stdio.h>
#include <stdlib.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"
#include "tran_variables.h"

#ifndef eV2Hartree
#define eV2Hartree    27.2113845                
#endif


/* 
 * write tran_transmission to the file "transmission"
 *   input tran_transmission_energyrange, 
 *         tran_transmission_energydiv
 *         tran_transmission 
 */

void TRAN_Output_Transmission(int SpinP_switch)
{
  int iw,k;
  dcomplex w;
  FILE *fp;

  if ( tran_transmission_on ) {

  char *fname="transmission";

  printf("TRAN_Output_Transmission file=%s\n",fname);

  if ( ( fp =fopen(fname,"w") )== NULL ) {
    printf("\ncan not open file to write transmission\n");
    printf("write transmission to stdout\n");
    fp = stdout;
  }


  fprintf(fp,"# SpinP_switch= %d\n",SpinP_switch);
  fprintf(fp,"# tran_transmission_energydiv= %d\n", tran_transmission_energydiv);

  for (k=0;k<=SpinP_switch; k++) {

    fprintf(fp,"# spin= %d\n",k);
    fprintf(fp,"# iw w.real(au) w.imag(au)  w.real(eV) w.imag(eV) trans.real trans.imag\n");

    for (iw=0;iw<tran_transmission_energydiv ; iw++) {

      w.r = tran_transmission_energyrange[0]+
        (tran_transmission_energyrange[1]-tran_transmission_energyrange[0])*
        (double)iw/(tran_transmission_energydiv-1);
      w.i = tran_transmission_energyrange[2];

      fprintf(fp,"%d %le %le %le %le %le %le\n", iw, w.r, w.i,
               w.r*eV2Hartree, w.i*eV2Hartree, 
              tran_transmission[k][iw].r,
                tran_transmission[k][iw].i );

    } /* iw */
    fprintf(fp,"\n\n");
  } /* k */

  if ( fp!=stdout && fp!=NULL )  fclose(fp);

  }



  /*********************************************************************
         write transmission_iv
   ********************************************************************/

  if ( tran_transmission_iv_on ) {

  char *fname="transmission_iv";

  printf("TRAN_Output_Transmission iv file=%s\n",fname);

  if ( ( fp =fopen(fname,"w") )== NULL ) {
    printf("\ncan not open file to write transmission\n");
    printf("write transmission to stdout\n");
    fp = stdout;
  }


  fprintf(fp,"# SpinP_switch= %d\n",SpinP_switch);
  fprintf(fp,"# tran_transmission_iv_energydiv= %d\n", tran_transmission_iv_energydiv);

  for (k=0;k<=SpinP_switch; k++) {

    fprintf(fp,"# spin= %d\n",k);
    fprintf(fp,"# iw w.real(au) w.imag(au)  w.real(eV) w.imag(eV) trans.real trans.imag\n");

    for (iw=0;iw<tran_transmission_iv_energydiv ; iw++) {

      w.r = tran_transmission_iv_energyrange[0]+
        (tran_transmission_iv_energyrange[1]-tran_transmission_iv_energyrange[0])*
        ((double)iw+0.5)/(tran_transmission_iv_energydiv-1);
      w.i = tran_transmission_iv_energyrange[2];

      fprintf(fp,"%d %le %le %le %le %le %le\n", iw, w.r, w.i,
               w.r*eV2Hartree, w.i*eV2Hartree,
              tran_transmission_iv[k][iw].r,
                tran_transmission_iv[k][iw].i );

    } /* iw */
    fprintf(fp,"\n\n");
  } /* k */

  if ( fp!=stdout && fp!=NULL )  fclose(fp);

  }


  printf("TRAN_Output_Transmission out\n");


}

