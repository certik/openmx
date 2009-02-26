/**********************************************************************
  TRNA_Calc_GC_LorR.c:
 
    TRAN_Calc_GC_LorR.c is a subroutine to calculate G_{C_L} or G_{C_R}
    in the non-equilibrium case.  

  Log of TRAN_Calc_GC_LorR.c:

     17/Dec/2005  Released by T.Ozaki, Taisuke Ozaki Copyright (C) 

***********************************************************************/

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
#include "lapack_prototypes.h"


/*
 *  calculate GCLorR
 *     GCLorR (w) = GC(w) ( w^* SC - HC - \Sigma_LorR(w^*) ) GC(w)^* 
 * 
 *
 *  no implicit variables 
 */


void TRAN_Calc_GC_LorR(
		       int iw_method,       /*  input */
                       dcomplex w,          /*  input */
                       double ChemP_e[2],   /*  input */
                       int nc,              /*  input */
                       int ne[2],           /*  input */
                       dcomplex *SigmaLorR, /*  input */
                       dcomplex *GC,        /*  input */
                       dcomplex *HCC,       /*  input */ 
                       dcomplex *SCC,       /*  input */ 
                       dcomplex *v1,        /* work, nc*nc */
                       dcomplex *GCLorR     /*  output */
                      )
#define GC_ref(i,j)  GC[nc*((j)-1)+(i)-1]
#define SigmaLorR_ref(i,j) SigmaLorR[nc*((j)-1)+(i)-1]
#define v1_ref(i,j)  v1[nc*((j)-1)+(i)-1] 
#define GCLorR_ref(i,j)  GCLorR[nc*((j)-1)+(i)-1]
#define HCC_ref(i,j) HCC[nc*((j)-1)+(i)-1]
#define SCC_ref(i,j) SCC[nc*((j)-1)+(i)-1]

{
  int i,j;
  int full_flag;
  int side;
  dcomplex alpha,beta;

  alpha.r = 1.0;
  alpha.i = 0.0;
  beta.r  = 0.0;
  beta.i  = 0.0;

  /* find a flag for how the H and S are divided */

  if (iw_method==3){
    if (ChemP_e[0]<ChemP_e[1])
      full_flag =  1;
    else 
      full_flag =  0;
  }
  else if (iw_method==4){
    if (ChemP_e[1]<ChemP_e[0])
      full_flag =  1;
    else 
      full_flag =  0;
  }

  if ( fabs(ChemP_e[1]-ChemP_e[0])<1.0e-50 ) full_flag = 2;

  /* GCLorR = -(\Sigma_LorR(w))^* */

  for (j=1; j<=nc; j++) {
    for (i=1; i<=nc; i++) {
      GCLorR_ref(i,j).r = -SigmaLorR_ref(i,j).r;
      GCLorR_ref(i,j).i =  SigmaLorR_ref(i,j).i;
    }
  }

  /* in case of the left side */

  if (iw_method==3){
    for (j=1; j<=(nc-ne[1]); j++) {
      for (i=1; i<=(nc-ne[1]); i++) {
        GCLorR_ref(i,j).r +=  w.r*SCC_ref(i,j).r - w.i*SCC_ref(i,j).i - HCC_ref(i,j).r; 
        GCLorR_ref(i,j).i += -w.i*SCC_ref(i,j).r - w.r*SCC_ref(i,j).i + HCC_ref(i,j).i; 
      }
    }
  }

  /* in case of the right side */

  else if (iw_method==4){
    for (j=(ne[0]+1); j<=nc; j++) {
      for (i=(ne[0]+1); i<=nc; i++) {
        GCLorR_ref(i,j).r +=  w.r*SCC_ref(i,j).r - w.i*SCC_ref(i,j).i - HCC_ref(i,j).r;
        GCLorR_ref(i,j).i += -w.i*SCC_ref(i,j).r - w.r*SCC_ref(i,j).i + HCC_ref(i,j).i;
      }
    }
  }

  /* correction of the central region */

  if (full_flag==0){
    for (j=(ne[0]+1); j<=(nc-ne[1]); j++) {
      for (i=(ne[0]+1); i<=(nc-ne[1]); i++) {
	GCLorR_ref(i,j).r -=  w.r*SCC_ref(i,j).r - w.i*SCC_ref(i,j).i - HCC_ref(i,j).r;
	GCLorR_ref(i,j).i -= -w.i*SCC_ref(i,j).r - w.r*SCC_ref(i,j).i + HCC_ref(i,j).i; 
      }
    }
  } 

  else if (full_flag==2){
    for (j=(ne[0]+1); j<=(nc-ne[1]); j++) {
      for (i=(ne[0]+1); i<=(nc-ne[1]); i++) {
	GCLorR_ref(i,j).r -= 0.5*(  w.r*SCC_ref(i,j).r - w.i*SCC_ref(i,j).i - HCC_ref(i,j).r );
	GCLorR_ref(i,j).i -= 0.5*( -w.i*SCC_ref(i,j).r - w.r*SCC_ref(i,j).i + HCC_ref(i,j).i );
      }
    }
  } 

  /* GC(w) ( w^* SC - HC -2 \Sigma_LorR(w^*) ) */

  F77_NAME(zgemm,ZGEMM)("N","N",&nc,&nc,&nc,&alpha, GC, &nc, GCLorR, &nc, &beta, v1, &nc);

  /* do complex conjugate */

  for (j=1;j<=nc;j++) {
    for (i=1;i<=nc;i++) {
      GC_ref(i,j).i = -GC_ref(i,j).i;
    }
  }

  /* GC(w) ( w^* SC - HC -2 \Sigma_LorR(w^*) ) GC(w)^* */
  
  F77_NAME(zgemm,ZGEMM)("N","N",&nc,&nc,&nc,&alpha, v1, &nc, GC, &nc, &beta, GCLorR, &nc);
  
  /* again do complex conjugate */
  
  for (j=1;j<=nc;j++) {
    for (i=1;i<=nc;i++) {
      GC_ref(i,j).i = -GC_ref(i,j).i;
    }
  }
}



