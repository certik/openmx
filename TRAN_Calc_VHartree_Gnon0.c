#include <stdio.h>
#include <math.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"

double Dot_Product(double a[4], double b[4]);

/*
 *  int d G_// exp(iG_// r_//) 
 *         { 
 *            dVH(G_//,x_l+1) ( exp(-Gx(x_(l+1) -x) -exp(-Gx( x_(l+1) -2x_0-x)) )
 *          + dVH(G_//,x_0) ( exp(-Gx(x-x_0)) - exp(-Gx(2x_(l+1)-x0-x) )
 *        }  /  ( 1-exp(-2Gx(x_(l+1)-x_0)) ) 
 *
 *   G_// .ne. 0 
 *
 *  dVH = VH_electrode-VH_center 
 *
 */


void TRAN_Calc_VHartree_Gnon0(
			      int k2, int k3,  /* k1.ne.0 and k2.ne.0 */
			      int ix0, int ixlp1, double gtv[4][4],
			      double   v1_r, double v1_i, double v2_r, double v2_i, 
			      dcomplex v1_c, dcomplex v2_c,
			      double *ReV2, double *ImV2) /* output, overwritten */
{

  int k1;
  double len;
  double G;
  double v[4];
  int i;
  double f1,f2,f3;
  double r1,r2,r3;
  int idim=1;

  len = Dot_Product(gtv[idim],gtv[idim]);
  len = sqrt(len);

  r1 = ix0*len;
  r2 = ixlp1*len;

  for (i=1;i<=3;i++) {
    v[i]=0.0;
  }
  for (i=1;i<=3;i++) {
    v[i]+=gtv[2][i]*k2+gtv[3][i]*k3;
  }
  G = sqrt(v[3]*v[3]+v[1]*v[1]+v[2]*v[2]);
     /* but, v[1] = 0 */


  for (k1=ix0+1; k1<ixlp1; k1++) {

    r3 = k1*len;

    f1 = exp(-G*(r2-r3)) - exp(-G*(r2-2.0*r1+r3));
    f2 = exp(-G*(r3-r1))-exp(-G*(2.0*r2-r1-r3));
    f3 = 1.0-exp(-2.0*G*(r2-r1));

    ReV2[k1] += ((v2_r-v2_c.r)*f1+(v1_r-v1_c.r)*f2)/f3;
    ImV2[k1] += ((v2_i-v2_c.i)*f1+(v1_i-v1_c.i)*f2)/f3;

  }

#if 0
  printf("Gnon0 Re %d %d\n",k1,k2);
  for (k3=iz0+1; k3<izlp1; k3++) {
      printf("%lf ",ReV2[k3]);
  }
  printf("Gnon0 Im %d %d\n",k1,k2);
  for (k3=iz0+1; k3<izlp1; k3++) {
      printf("%lf ",ImV2[k3]);
  }
#endif

}

