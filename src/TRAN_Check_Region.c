#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tran_variables.h"


/*  output: none
    return value:  1 = OK
                   0 = NG

   purpose:
   to confirm that ...|L1|L|C|R|R1|...  has no overlapping PAO over the nearest neighboring region.

*/
int TRAN_Check_Region(
		      int atomnum,
		      int *WhatSpecies, 
		      double *Spe_Atom_Cut1,
		      double **Gxyz
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

  int iregion;

  double A[4], B[4], diff[4];
  double len;

  int i;

  char *regionstr[2]={"region L","region R"};

  int error;

  error=0;


  /* center <-> L1 or R1 */

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++) {

    if ( TRAN_region[ct_AN]%10 ==1 ) { /* center */

      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];
      for (i=1; i<=3; i++) { A[i]= Gxyz[ct_AN][i]; }

      for (ct_BN=1; ct_BN<=atomnum; ct_BN++) {

	if ( TRAN_region[ct_BN]%10 == 2  || TRAN_region[ct_BN]%10 == 3 ) { /* left or right */

	  if ( TRAN_region[ct_BN]%10 == 2 ) {ix=-1; iregion=0;}/* left */
	  if ( TRAN_region[ct_BN]%10 == 3 ) {ix=+1; iregion=1;}/* right */

	  wanB = WhatSpecies[ct_BN];
	  rcutB = Spe_Atom_Cut1[wanB];
       
	  rcutAB = rcutA + rcutB;

	  for (iy=-ny;iy<=ny;iy++) {
	    for (iz=-nz;iz<=nz;iz++) {

	      for (i=1;i<=3;i++) {
		B[i] = Gxyz[ct_BN][i] + tv_e[iregion][1][i]*ix
             		              + tv_e[iregion][2][i]*iy
		                      + tv_e[iregion][3][i]*iz;
	      } 

	      for (i=1;i<=3;i++) { diff[i] = A[i]-B[i]; }
	      len = sqrt( diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3] );
	      if ( len <= rcutAB ) { 
                printf("\n\nTRAN_Check_Region_Lead()\n");
                printf("\nThe length between atomA=%d(region C) and atomB=%d(%s) is too short for the transport calculation.\n",
		       ct_AN, ct_BN,regionstr[iregion]); 
                printf("distance=%lf rcutA=%lf rcutB=%lf\n",len,rcutA,rcutB);
                error=1;
		goto LRcheck;
              }
	     
	    } /* ix */
	  } /* iy */

	} /* left or right */
      } /* ct_BN */

    } /* center */
  } /* ct_AN */


 LRcheck: 

  /* L <-> R */
  for (ct_AN=1;ct_AN<=atomnum;ct_AN++) {
    if ( TRAN_region[ct_AN]%10 ==2 ) { /* left */
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];
      for (i=1;i<=3;i++) { A[i]= Gxyz[ct_AN][i]; }

      for (ct_BN=1;ct_BN<=atomnum;ct_BN++) { 
        if ( TRAN_region[ct_AN]%10 ==3 ) { /* right */
          wanB = WhatSpecies[ct_BN];
          rcutB = Spe_Atom_Cut1[wanB];

          rcutAB = rcutA+rcutB;

          ix =0;
          for (iy=-ny;iy<=ny;iy++) {
            for (iz=-nz;iz<=nz;iz++) {
              for (i=1;i<=3;i++) {

                B[i] = Gxyz[ct_BN][i] + tv_e[iregion][1][i]*ix
                                      + tv_e[iregion][2][i]*iy
                                      + tv_e[iregion][3][i]*iz;
              }
              for (i=1;i<=3;i++) { diff[i] = A[i]-B[i]; }
              len = sqrt( diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3] );
              if ( len <= rcutAB ) { 
                printf("\n\nTRAN_Check_Region_Lead()\n");
                printf("\nThe length between atomA=%d(region L) and atomB=%d(region R) is too short for the transport calculation.\n",ct_AN, ct_BN);
                printf("distance=%lf rcutA=%lf rcutB=%lf\n",len,rcutA,rcutB);
                error=1;
		goto lastproc;
              }

	    }/* iz */
	  } /* iy */
        } /* right */
      } /* ct_BN */
    } /* left */
  } /* ct_AN */

 lastproc:

  if (error) { return 0; }
  else { return 1; }
}
