#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"

void TRAN_Calc_CentGreen(
                      dcomplex w,
                      int nc,
                      dcomplex *sigmaL,
                      dcomplex *sigmaR,
                      dcomplex *HCC,
                      dcomplex *SCC,
                      dcomplex *GC /* output */
                      )
;

/* interface from fortran */
void tran_calc_centgreen( dcomplex *w, int *nc, dcomplex *sigmaL, dcomplex *sigmaR,
                          dcomplex *HCC, dcomplex *SCC, dcomplex *GC)
{
  TRAN_Calc_CentGreen( *w, *nc, sigmaL, sigmaR, HCC, SCC, GC); 
}


/*
 *  calculate G(w) = ( w SCC - HCC - SigmaL -SigmaR ) ^-1 
 *
 * charge is calculated as 
 *     -1/PI int dw Im G(w+id) 
 *        -1/PI *dw is included in w_weight[]
 *
 *  no implicit variables 
 */
void TRAN_Calc_CentGreen(
			 /* input */
			 dcomplex w,
			 int nc, 
			 dcomplex *sigmaL,
			 dcomplex *sigmaR, 
			 dcomplex *HCC,
			 dcomplex *SCC,
			 dcomplex *GC /* output */
			 )
#define HCC_ref(i,j) HCC[nc*((j)-1)+(i)-1]
#define SCC_ref(i,j) SCC[nc*((j)-1)+(i)-1]

#define GC_ref(i,j)  GC[nc*((j)-1)+(i)-1]
#define sigmaL_ref(i,j) sigmaL[nc*((j)-1)+(i)-1]
#define sigmaR_ref(i,j) sigmaR[nc*((j)-1)+(i)-1]

{
  int i,j;
  int pos;

  /* w SCC - HCC - SigmaL -SigmaR */
    
  for (i=1;i<=nc;i++) {
    for (j=1;j<=nc;j++) {

#if 0
      GC_ref(i,j).r = w.r*SCC_ref(i,j).r - w.i*SCC_ref(i,j).i - HCC_ref(i,j).r ;
      GC_ref(i,j).i = w.r*SCC_ref(i,j).i + w.i*SCC_ref(i,j).r - HCC_ref(i,j).i ;
#else
      GC_ref(i,j).r = w.r*SCC_ref(i,j).r - w.i*SCC_ref(i,j).i - HCC_ref(i,j).r
                     - sigmaL_ref(i,j).r - sigmaR_ref(i,j).r;
      GC_ref(i,j).i = w.r*SCC_ref(i,j).i + w.i*SCC_ref(i,j).r - HCC_ref(i,j).i
                     - sigmaL_ref(i,j).i - sigmaR_ref(i,j).i;
#endif

    }
  }

  Lapack_LU_Zinverse(nc,GC);
}



