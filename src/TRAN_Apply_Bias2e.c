#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"

/*
 apply source-drain voltage to the electrode
*/


void  TRAN_Apply_Bias2e(
       MPI_Comm comm1,

       double voltage, /* applied bias */

       int SpinP_switch,
       int atomnum,

       int *WhatSpecies,
       int *Spe_Total_CNO,
       int *FNAN,
       int **natn,
       
       int Ngrid1,
       int Ngrid2,
       int Ngrid3,

        double ****OLP,

                  /* output: overwritten */
        double *ChemP, 
        double *****H, 
        double *dVHart_Grid
)
{
  if (fabs(voltage)<1.0e-50) return;

  printf("add voltage =%lf to the electrode\n",voltage);


  {
    int GA_AN,wanA,tnoA;
    int LB_AN,GB_AN,wanB,tnoB; 
    int i,j,k; 
    int myid;

    MPI_Comm_rank(comm1,&myid);

    /* ChemP */
    (*ChemP) +=  voltage; 

    /*   <i|H+V|j> = <i|H|j> + V<i|j>  */
    
    for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];

        for (i=0; i<tnoA; i++){
          for (j=0; j<tnoB; j++){
            for (k=0;k<=SpinP_switch; k++) {
	      H[k][GA_AN][LB_AN][i][j] += voltage*OLP[GA_AN][LB_AN][i][j];
            }
          }
        } /* i */
      } /* LB_AN */
    }   /* GA_AN */
  }

  { 
    int i;

    /* grid_value += voltage */

    for (i=0;i<Ngrid1*Ngrid2*Ngrid3; i++) {
      dVHart_Grid[i] += voltage ; 
    }
  }

}



