#include <stdio.h>
#include <math.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"

double  Dot_Product(double a[4], double b[4]); 

/*
 *    effect of the boundary condition
 * 
 *
 *
 *  VH = VH0 + dVH, VH0 = calculatef from rho_C(x,y,z) 
 *                  dVH = effect of the boundary condition 
 *
 *  calculate dVH from
 *   V_electrode(xi,ky,kz) - VH_central(xi,ky,kz),  xi=x_0 and x_(l+1) 
 *
 *----------------------------------
 *
 *  Solve d^2/dx^2 dVH(x) = 0, with the boundary condition,dVH(x0) and dVH(x_l+1)
 *
 *  d/dx dVH(x) = a
 *  dVH(x) = a(x-x0) + b 
 * 
 *  dVH(x0) = b
 *  dVH(x_l+1) = a(x_l+1 -x_0) + b 
 *
 */

void TRAN_Calc_VHartree_G0(dcomplex vl, dcomplex vr,   /* value of the electrode */
                      double vl_c_r,  double vl_c_i,   /* value of left edge (x_0) of the C region , calculated by FFT */
                      double vr_c_r,  double vr_c_i,   /* value of right edge the (x_(l+1)) C region ,calculated by FFT */ 
                      int ix0, int ixlp1  , double gtv1[4],
                      double *pot_r, double *pot_i  /* input,output  ---  overwrite */
)
{
   double r1,r2;
   double len;
   dcomplex vb,va;
   int i1;

/*debug
 *   printf("vl=%le %le, vr= %le %le\n", vl.r,vl.i, vr.r, vr.i);
 *  printf("vl_c_r=%le %le\n", vl_c_r,vl_c_i);
 *  printf("vr_c_r=%le %le\n", vr_c_r,vr_c_i);
 */


   len = Dot_Product(gtv1,gtv1);
   len = sqrt(len);


   r1 =  len*ix0;
   r2 = len*ixlp1;


   vb.r = vl.r - vl_c_r;
   vb.i = vl.i - vl_c_i;

   va.r =( vr.r-vr_c_r-vb.r )/ (r2-r1);
   va.i =( vr.i-vr_c_i-vb.i )/ (r2-r1);

   for (i1=ix0+1;i1<ixlp1;i1++) {
       r2 = len*i1;
       pot_r[i1] += va.r*(r2-r1)+vb.r;  /* add it */
       pot_i[i1] += va.i*(r2-r1)+vb.i;  /* add it */
   }

#if 0
  printf("G0 Re\n");
   for (i3=iz0+1;i3<izlp1;i3++) {
     printf("%lf ",pot_r[i3]);
   }
  printf("G0 Im\n");
   for (i3=iz0+1;i3<izlp1;i3++) {
     printf("%lf ",pot_i[i3]);
   }

#endif


}






