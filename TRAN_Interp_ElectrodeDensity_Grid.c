

/************************************************
  set density_grid on the electrode 

  interpolate Density_Grid_e 
  and setup Density_Grid on the electrodes  ( mesh points are different )

***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#include "tran_prototypes.h"

#define print_stdout 0



double Dot_Product(double a[4], double b[4]);
void Cross_Product(double a[4], double b[4], double c[4]);

#define PI            3.1415926535897932384626


static void    Convert_Cell_tv(double tv[4][4], double pos[4],
  /* output*/ double xyz[4])
{
  double a[4];
  double rtv[4][4]; 
  double tmpv[4], CellV;

  Cross_Product(tv[2],tv[3],tmpv);
  CellV = Dot_Product(tv[1],tmpv); 
  
  Cross_Product(tv[2],tv[3],tmpv);
  rtv[1][1] = 2.0*PI*tmpv[1]/CellV;
  rtv[1][2] = 2.0*PI*tmpv[2]/CellV;
  rtv[1][3] = 2.0*PI*tmpv[3]/CellV;

  Cross_Product(tv[3],tv[1],tmpv);
  rtv[2][1] = 2.0*PI*tmpv[1]/CellV;
  rtv[2][2] = 2.0*PI*tmpv[2]/CellV;
  rtv[2][3] = 2.0*PI*tmpv[3]/CellV;
  
  Cross_Product(tv[1],tv[2],tmpv);
  rtv[3][1] = 2.0*PI*tmpv[1]/CellV;
  rtv[3][2] = 2.0*PI*tmpv[2]/CellV;
  rtv[3][3] = 2.0*PI*tmpv[3]/CellV;

  if (print_stdout){
    printf("Convert_Cell_tv, pos= %lf %lf %lf ", pos[1], pos[2], pos[3]);
  }

  a[1] = Dot_Product(pos,rtv[1])*0.5/PI;
  a[2] = Dot_Product(pos,rtv[2])*0.5/PI;
  a[3] = Dot_Product(pos,rtv[3])*0.5/PI;

  if (print_stdout){
    printf("x= %lf %lf %lf ", a[1], a[2], a[3]);
  }

  a[1] *= sqrt(Dot_Product(tv[1],tv[1]));
  a[2] *= sqrt(Dot_Product(tv[2],tv[2]));
  a[3] *= sqrt(Dot_Product(tv[3],tv[3]));

  xyz[1]= a[1];
  xyz[2]= a[2];
  xyz[3]= a[3];

  if (print_stdout){
    printf("x*|tv|= %lf %lf %lf\n", xyz[1], xyz[2], xyz[3]);
  }
}



/**************************************************

   N_e = Ngrid1_e*  Ngrid2_e* Ngrid3_e
   Density_Grid_e[][N_e]
   mesh point is  tv_e[4][4]

   N = Ngrid1*  Ngrid2* Ngrid3
   Density_Grid[][N]
   mesh point is  tv[4][4]

 purpose:
   interpolate Density_Grid_e to fit the mesh points of Density_Grid


   This subroutine is to be used also so as to interpolate data of Vpot_Grid. 


 note: 
   Density_Grid_ref(i,j,k) = Density_Grid[MN] 
            MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
     i=0...Ngrid1-1
     j=
     k=

**************************************************/


void TRAN_Interp_ElectrodeDensity_Grid(
    MPI_Comm   comm1,
   double offset[4],  /* input. use offset[1..3],  = dx */
                               /*    left edge is gtv*offset */

   double *Density_Grid_e ,  /* input */
   int Ngrid1_e, int Ngrid2_e, int Ngrid3_e, /* input */
   double tv_e[4][4],  /* input, unit cell vector of electrode, from inputfile not from saved data */
   double gtv_e[4][4],  /* input,  cell vector of electrode */

   double *Density_Grid, /* <=== output */
   int Ngrid1, int Ngrid2, int Ngrid3,
    int l1boundary[2], /* input */
   double tv[4][4], /* input, unit cell vector */
   double gtv[4][4] /* input,  cell vector */
)
{

   int i,j,k;
   int id;
   int po; /* error flag */
   int n_in[3], n_out[4];
   double angle,sign,factor;
   double alpha[3],beta[3],gamma[3];
   double grids[4],grids_e[4];

   double Xgtv_e[4], Xgtv[4], Xtv_e[4], Xoffset[4];


   /* mesh on the electrode is gtv_e*(1/2+integer) */

   /* mesh on the CLR system is gtv*(1/2+integer) */

   /* first check tv[1] and tv[2] have the same direction as tv_e[1] and tv_e[2] */
   po=0; /* error flag */
   for (id=2;id<=3;id++) {
      for (i=1;i<=3;i++) {
         if ( 1.0e-20<fabs(tv_e[id][i]-tv[id][i]) ) po++;
       }
    }
    if (po) {
        /* tv and tv_electrode is different */
        printf("a and/or b axis of tv and tv_electrode is different\n");
        exit(0);
     }


  /* check whether the direction of tv[1] and tv_e[1] is 0 or pi */
  angle =  Dot_Product( tv[1], tv_e[1]);
  if (angle<0) {
     printf("c direction of this system is different from that of the electrode\n");
     exit(0);
  }

  /***************************************************************************
   *     positions 
   *     electrode mesh                    LCR mesh 
   *     (1/2+i_e)*gtv_e  => dx + (1/2+i0+i)*gtv,   integer: i_e,i,i0
   *    i_e = [0:Ngrid_e-1]              i=[0:N1-1]
   * 
   *     (1/2+i_e)*tv_e/Ngrid_e 
   *     
   *     [(1/2+i_e)/Ngrid_e] *tv_e      
   *                                     dx + (1/2+i0+i)*gtv
   *                                    = (dx + (1/2+i0+i)*gtv )/tv_e * tv_e
   *                                    = (dx/tv_e + (1/2+i0+i)*gtv/tv_e) * tv_e
   *                                    = (dx/tv_e + (1/2+i0+i)*gtv/gtv_e/Ngrid) *tv_e
   *                                    = {(dx/tv_e + (1/2+i0+i)*gtv/gtv_e/Ngrid)*N1 /N1 } *tv_e
   *
   *                 
   *                                 
   *      0<[(1/2+i_e)/Ngrid_e]<1       
   *                                  alpha =  dx/tv_e*Ngrid + (1/2+i0)*tv/tv_e, beta = tv/tv_e 
   *
   *   
   ***************************************************************************/


  for (i=1;i<=3;i++) {
      gamma[i-1]=0.5;
  }
  grids_e[1]=Ngrid1_e;
  grids_e[2]=Ngrid2_e;
  grids_e[3]=Ngrid3_e;

  grids[1] = l1boundary[1]-l1boundary[0]+1;
  grids[2] = Ngrid2;
  grids[3] = Ngrid3;
      
#if 1
  Convert_Cell_tv(tv,offset,Xoffset);
  for (i=1;i<=3;i++) {
       Xgtv[i]= sqrt(Dot_Product(gtv[i],gtv[i]));
       Xtv_e[i] = sqrt(Dot_Product(tv_e[i],tv_e[i]));
       Xgtv_e[i] = sqrt(Dot_Product(gtv_e[i],gtv_e[i]));
  }
  for (i=1;i<=3;i++) {
      factor = Xgtv[i]/Xgtv_e[i]/(double)grids_e[i]*(double)grids[i];
      if (i==1) {
      alpha[i-1]=Xoffset[i]/Xtv_e[i]*grids[i]+ (0.5+l1boundary[0])*factor;
      }
      else {
      alpha[i-1]=Xoffset[i]/Xtv_e[i]*grids[i]+ (0.5)*factor;
      }
                   /* start from l3boundary[0] */
      beta[i-1] = factor;
  }

#else
  for (i=1;i<=3;i++) {
      factor = sqrt(Dot_Product(gtv[i],gtv[i]))/sqrt(Dot_Product(gtv_e[i],gtv_e[i]))/grids_e[i]*grids[i];
      alpha[i-1]=offset[i]/sqrt(Dot_Product(tv_e[i],tv_e[i]))*grids[i]+(0.5+l3boundary[0])*factor; 
                   /* start from l3boundary[0] */
      beta[i-1] = factor;
  }
#endif

/* debug */
/*   alpha[2]=0.5*(l3boundary[1]-l3boundary[0]+1)/Ngrid3_e;
 *   beta[2]=1.0*(l3boundary[1]-l3boundary[0]+1)/Ngrid3_e;
 */


  if (print_stdout){
    for (i=0;i<3;i++) {
      printf("alpha,beta,gamma=%lf %lf %lf\n",alpha[i],beta[i],gamma[i]);
    }
  }

  /* interpolate it */


     /* interpolation */
  n_in[0]=Ngrid1_e;
  n_in[1]=Ngrid2_e;
  n_in[2]=Ngrid3_e;
  n_out[0]=l1boundary[1]-l1boundary[0]+1;
  n_out[1]=Ngrid2;
  n_out[2]=Ngrid3;
  n_out[3]=l1boundary[1]-l1boundary[0];

  if (print_stdout){
    printf("nin = %d %d %d\n", n_in[0],  n_in[1], n_in[2]);
    printf("nout = %d %d %d %d\n", n_out[0],  n_out[1], n_out[2],n_out[3]);
  }

  TRAN_FFT_Dinterpolation3D(comm1, n_in, Density_Grid_e, n_out, Density_Grid, 
                            alpha,beta,gamma);
}



