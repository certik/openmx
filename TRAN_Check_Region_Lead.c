#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tran_variables.h"


/*  output: none
    return value:  1 = OK
                   0 = NG

   purpose:
   to confirm that ...|L0|L1|L2|L3|... has overlapping PAO only between  L1 and L2.

*/
int TRAN_Check_Region_Lead(
  int atomnum,
  int *WhatSpecies, 
  double *Spe_Atom_Cut1,
  double **Gxyz,
  double tv[4][4]
)
{

   int ct_AN;
   int wanA;
   double rcutA;

   int ct_BN;
   int wanB;
   double rcutB;

   double rcutAB;

   int ix,iy,iz;
   int nz=3,ny=3;

   double A[4], B[4], diff[4];
   double len;

   int i;


   /* this routine only in the case of TRAN_output_hks!=0 */
   if ( TRAN_output_hks ==0 ) {
          return 1;
   }


   for (ct_AN=1;ct_AN<=atomnum;ct_AN++) {
     wanA = WhatSpecies[ct_AN];
     rcutA = Spe_Atom_Cut1[wanA];
     for (i=1;i<=3;i++) { A[i]= Gxyz[ct_AN][i]; }
     for (ct_BN=1;ct_BN<=atomnum;ct_BN++) {
       wanB = WhatSpecies[ct_BN];
       rcutB = Spe_Atom_Cut1[wanB];
       
       rcutAB = rcutA+rcutB;

       for (ix=-2;ix<=2; ix+=4) {
	 for (iy=-ny;iy<=ny;iy++) {
	   for (iz=-nz;iz<=nz;iz++) {
	     for (i=1;i<=3;i++) { B[i] = Gxyz[ct_BN][i]+ tv[1][i]*ix + tv[2][i]*iy+tv[3][i]*iz; } 
	     for (i=1;i<=3;i++) { diff[i] = A[i]-B[i]; }
	     len = sqrt( diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3] );
	     if ( len <= rcutAB ) { 
                printf("\n\nTRAN_Check_Region_Lead()\n");
	        printf("\nThe length between atomA=%d and atomB=%d is too short for the transport calculation.\n",ct_AN, ct_BN); 
		printf("distance=%lf rcutA=%lf rcutB=%lf\n",len,rcutA,rcutB);
	        return 0; 
	     }
	     
	   } /* ix */
	 } /* iy */
       } /* iz */

     } /* ct_BN */

   } /* ct_AN */

   return 1;
}



