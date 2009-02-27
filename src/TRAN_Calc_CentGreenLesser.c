#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "lapack_prototypes.h"



/*
 *  calculate G^<(w) 
 *     G^< (w) = G^R(w) ( Gamma_L(w) nf(w-ChemP_L) + Gamma_R(w) nf(w-ChemP_R) ) G^A(w) 
 * 
 *  ... w_weight includes -1 factor 
 *    rho(w) = -1/PI Im G(w) 
 *
 *  ... therefore, this function returns -G^< (w) 
 *
 *
 *  no implicit variables 
 */
void TRAN_Calc_CentGreenLesser(
                                /* input */
                      dcomplex w,
                      double ChemP_e[2],

                      int nc, 
                      dcomplex *SigmaL,
                      dcomplex *SigmaR, 
                      dcomplex *GC, /*  input */
                      dcomplex *v1, /* work, nc*nc */
                      dcomplex *Gless /*  output */
                      )
#define GC_ref(i,j)  GC[nc*((j)-1)+(i)-1]
#define SigmaL_ref(i,j) SigmaL[nc*((j)-1)+(i)-1]
#define SigmaR_ref(i,j) SigmaR[nc*((j)-1)+(i)-1]
#define v1_ref(i,j)  v1[nc*((j)-1)+(i)-1] 
#define Gless_ref(i,j)  Gless[nc*((j)-1)+(i)-1]

{
  int i,j;
  int side;
  dcomplex alpha,beta;

  alpha.r=1.0;
  alpha.i=0.0;
  beta.r=0.0;
  beta.i=0.0;


  for (j=1;j<=nc;j++) {
    for (i=1;i<=nc;i++) {
      Gless_ref(i,j).r = 0.0;
      Gless_ref(i,j).i = 0.0;
    }
  }

  side=0;
  if ( w.r <= ChemP_e[side] )  {

    for (j=1;j<=nc;j++) {
      for (i=1;i<=nc;i++) {
	Gless_ref(i,j).r += 0.0;   /* Gamma_L */
	Gless_ref(i,j).i += -SigmaL_ref(i,j).i*2.0; 
      }
    }

  }
#ifdef DEBUG
  TRAN_Print2_dcomplex("self+",nc,nc,Gless);
#endif


  side=1;
  if ( w.r <= ChemP_e[side] )  {

    for (j=1;j<=nc;j++) {
      for (i=1;i<=nc;i++) {
	Gless_ref(i,j).r +=  0.0; 
	Gless_ref(i,j).i +=  -SigmaR_ref(i,j).i*2.0;   /* Gamma_R */
      }
    }

  }

#ifdef DEBUG
  TRAN_Print2_dcomplex("self+self",nc,nc,Gless);
#endif

  if (w.r > ChemP_e[0] && w.r >  ChemP_e[1] )  {
    /* then Gless =0 as it is here */
    return;
  }


  /* now Gless = Gamma_L f(w-EF_L) + Gamma_R f(w-EF_R) */

  F77_NAME(zgemm,ZGEMM)("N","N",&nc,&nc,&nc,&alpha, GC, &nc, Gless, &nc, &beta, v1, &nc);
  /*   G^R(w) * ( Gamma_L f(w-EF_L) + Gamma_R f(w-EF_R)  ) */

#ifdef DEBUG
  TRAN_Print2_dcomplex("GxSelf",nc,nc,v1);
#endif


  for (j=1;j<=nc;j++) {
    for (i=1;i<=nc;i++) {
      GC_ref(i,j).i = - GC_ref(i,j).i;   /* G^A(w) */
    }
  }

  F77_NAME(zgemm,ZGEMM)("N","N",&nc,&nc,&nc,&alpha, v1, &nc, GC, &nc, &beta, Gless, &nc);
  /* v1=  G^R(w) * ( Gamma_L f(w-EF_L) + Gamma_R f(w-EF_R)  )  G^A(w) */

  for (j=1;j<=nc;j++) {
    for (i=1;i<=nc;i++) {
      GC_ref(i,j).i = - GC_ref(i,j).i;   /* G^R(w) */
    }
  }

}



