#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"
#include "tran_variables.h"
#include "lapack_prototypes.h"

#define PI            3.1415926535897932384626

typedef struct {
    double a,b;
} dlists;



/* Taisuke Ozaki Copyright (C) */
static void Eigen_DGGEVX( int n, double **a, double **s, double *eval,
                          double *resr, double *resi );
/* Taisuke Ozaki Copyright (C) */
static void zero_cfrac( int n, dcomplex *zp, dcomplex *Rp );
/* Taisuke Ozaki Copyright (C) */
static int dlists_cmp(const dlists *x, const dlists *y);
/* Taisuke Ozaki Copyright (C) */
static void dqsort(long n, double *a, double *b);




/*
*  (0,3) ->   (1,3)
*    A          |
*    |          V
*  (0,2)      (1,2) ->(1+|VR-VL|,2)
*/

void TRAN_Set_IntegPath_Square(void)
{
  static int method=1; /* log or linear for the (1,3)->(1,2) path  */

  int i,itotal;
  double x;

  double fac;
  dcomplex dx;


  tran_omega_n_scf=0;
  for (i=0;i<=3;i++) {
    tran_omega_n_scf += tran_square_path_div[i];
  }

  if (tran_bias_apply) {
    tran_omega_n_scf += tran_square_path_bias_div;
  }


     
  /* allocation */
  tran_omega_scf = (dcomplex*)malloc(sizeof(dcomplex)*tran_omega_n_scf);
  tran_omega_weight_scf=(dcomplex*)malloc(sizeof(dcomplex)*tran_omega_n_scf);
  tran_integ_method_scf=(int*)malloc(sizeof(int)*tran_omega_n_scf);


  /* calculate freq */  

  /* spec = -1/PI Im G^R(w) */
  fac = 1.0/PI;
  itotal=0;

  /* (0,2)->(0,3) divided by div[0]  */

  dx.r = 0.0;
  dx.i = ( tran_square_path_ene[3]- tran_square_path_ene[2] ) / tran_square_path_div[0];
   
  for (i=0;i<tran_square_path_div[0];i++) {
    x = (double)i+0.5;
    tran_omega_scf[itotal].r = tran_square_path_ene[0] + dx.r*x;
    tran_omega_scf[itotal].i = tran_square_path_ene[2] + dx.i*x;

    tran_omega_weight_scf[itotal].r = -dx.r*fac;
    tran_omega_weight_scf[itotal].i = -dx.i*fac;

    tran_integ_method_scf[itotal] = 1; /* equiv */

    itotal++;

  }

  /* (0,3) -> (1,3)  by div[1] */

  dx.r = ( tran_square_path_ene[1]- tran_square_path_ene[0] ) / tran_square_path_div[1];
  dx.i = 0.0;

  for (i=0;i<tran_square_path_div[1];i++) {
    x=(double)i+0.5;

    tran_omega_scf[itotal].r = tran_square_path_ene[0] + dx.r*x;
    tran_omega_scf[itotal].i = tran_square_path_ene[3] + dx.i*x;

    tran_omega_weight_scf[itotal].r = -dx.r*fac;
    tran_omega_weight_scf[itotal].i = -dx.i*fac;

    tran_integ_method_scf[itotal] = 1; /* equiv */

    itotal++;

  }

  /* (1,3) -> (1,2) by div[2] */ 
  if (method==1) {  /* linear */
    dx.r = 0.0;
    dx.i = ( tran_square_path_ene[2]- tran_square_path_ene[3] ) / tran_square_path_div[2];

    for (i=0;i<tran_square_path_div[2];i++) {
      x=i+0.5;
      tran_omega_scf[itotal].r = tran_square_path_ene[1] + dx.r*x;
      tran_omega_scf[itotal].i = tran_square_path_ene[3] + dx.i*x;

      tran_omega_weight_scf[itotal].r = -dx.r*fac;
      tran_omega_weight_scf[itotal].i = -dx.i*fac;

      tran_integ_method_scf[itotal] = 1; /* equiv */

      itotal++;

    }

  }
  else if (method==2) {  /* log */
    /*
      c                  assuming                         
      c                     theta(1/2+i) = exp(th1+delta(1/2+i))
      c                boundary condition
      c                     theta(0) = exp(th1) = exp( log(w4) ) = w4
      c                     theta(N3)= exp(th1+delta*N3)) = exp(th1+(th2-th1)/N3*N3) = exp(th2)= w3
      c                                   delta = (th2-th1)/N3
      c                integral weight 
      c                     int d theta = int  exp(th1+delta x) delta dx 
      c                             =>  weight = exp(th1+delta(1/2+i))delta  for ith 
      c
    */
    double th1, th2,delta;

    th1 = log( tran_square_path_ene[3] );
    th2 = log (tran_square_path_ene[2] );

    delta = (th2-th1)/tran_square_path_div[2]; 

    for (i=0;i<tran_square_path_div[2];i++) {
      x = 0.5+i;
      tran_omega_scf[itotal].r = tran_square_path_ene[1] ;
      tran_omega_scf[itotal].i = exp(th1+delta*x);

      tran_omega_weight_scf[itotal].r = 0.0;
      tran_omega_weight_scf[itotal].i = -delta*exp(th1+delta*x)*fac;

      tran_integ_method_scf[itotal] = 1; /* equiv */

      itotal++;
    }

  } /* method */


  if ( tran_bias_apply ) {
    /* bias case,  (NEGF) */
    /*        path (1,2)-> (| V_L - V_R | + biasV ,2) */
    dx.r = ( fabs(tran_biasvoltage_e[0]-tran_biasvoltage_e[1]) + 
             tran_square_path_bias_expandenergy ) / tran_square_path_bias_div;
    dx.i =  0;

    for (i=0;i<tran_square_path_bias_div;i++) {
      x=(double)i+0.5;

      tran_omega_scf[itotal].r = tran_square_path_ene[1] + dx.r*x;
      tran_omega_scf[itotal].i = tran_square_path_ene[2] + dx.i*x;

      /*  -i/pi  int  A dw 
       *     -i dw    =  (dw.r+ i dw.i)*(-i) = ( dw.i -i dw.r )
       */
      tran_omega_weight_scf[itotal].r =  dx.i*fac;
      tran_omega_weight_scf[itotal].i = -dx.r*fac;

      tran_integ_method_scf[itotal] = 2; /* non-equiv. */

      itotal++;

    }

  }


  /*debug*/
  printf("integral path\n");
  printf("id  freq                weight\n");
  for (i=0;i<tran_omega_n_scf;i++) {
    printf("%d %lf %lf %lf %lf\n",i, tran_omega_scf[i].r, tran_omega_scf[i].i, 
	   tran_omega_weight_scf[i].r, tran_omega_weight_scf[i].i);
  }
  printf("\n");
  /*end debug*/ 
  
}





/*****************************************************************************************
*****************************************************************************************
*****************************************************************************************/

/*
  spectra = -1/pi Im G^R(w)
*/

/*  
       ----------
     /            \ (1-gamma,3) 
    /              -------------(1,3)
    |                x
    |                x
    |                x
    ---+------------------------------>
    (0,2)

   n = int_C A(x) nF(x) dx 

  (0,2) -> (1-gamma,3)   divided by n0
  (1-gamma,3)->(1,3)     divided by n1

   if (applied_bias)     additional points

 use
double tran_thermalarc_path_ene[5];
int tran_thermalarc_path_ene_fix[5];
int tran_thermalarc_path_div[2];
double tran_thermalarc_path_bias_expandenergy;  for NEGF 
int tran_thermalarc_path_bias_div;     for NEGF 

*/
void      TRAN_Set_IntegPath_ThermalArc(void)
{

  double f,factor,beta;
  double chem;
  int n_max;
  int i, itotal;
  double gamma, w1,w2,w3,w4;
  double a,b;
  double th1,th2,dth;
   

  /* find chem=chemical potential */
  chem = ( ChemP_e[0] < ChemP_e[1] )? ChemP_e[0] : ChemP_e[1] ;
  if ( tran_bias_apply ) {
    chem -= tran_thermalarc_path_bias_expandenergy; 
  }


  n_max=100;

  f= PI*(2*n_max+1)*tran_temperature; 
  if ( f<  tran_thermalarc_path_ene[3] ) {
    printf("PI*(2*n_max+1)*tran_temperature (%lf) < tran_thermalarc_path_ene[3] (%lf)\n",
	   f,tran_thermalarc_path_ene[3]);
    exit(0);
  }

  /* count max thermal freq */
  for (i=0;i<=n_max;i++) {
    f = (2.0*(double)i+1.0)*PI*tran_temperature;
    if (f> tran_thermalarc_path_ene[3]) {
      break;
    }
  }
  n_max = i;


  /* storage size */
  tran_omega_n_scf=0;
  for (i=0;i<=3;i++) {
    tran_omega_n_scf += tran_square_path_div[i];
  }
  
  tran_omega_n_scf += n_max; /* contribution of thermal freq */

  

  if (tran_bias_apply) {
    tran_omega_n_scf += tran_square_path_bias_div;
  }
  

  /* allocation */
  tran_omega_scf = (dcomplex*)malloc(sizeof(dcomplex)*tran_omega_n_scf);
  tran_omega_weight_scf=(dcomplex*)malloc(sizeof(dcomplex)*tran_omega_n_scf);
  tran_integ_method_scf=(int*)malloc(sizeof(int)*tran_omega_n_scf);


  /* calculate freq */  

  /* spec = -1/PI Im G^R(w) */
  factor = 1.0/PI;
  itotal=0;


  /************
   1st contribution

   calculate arc 

        -----------
      /             \
     /                \(2,4)-(gamma,0) 
    |
    |
    ---------+---------------------------
    (1,3)   (b,0)

    radius = a 
    center = (b,0)

    a^2 = (w1-b)^2+w3^2
    a^2 = (w2-gamma-b)^2+w4^2

    0 = w1^2 -2 w1 b +  w3^2
    - (  (w2-gamma)^2 - 2(w2-gamma) b  + w4^2 )

    b = ( w1^2 + w3^2 - (w2-gamma)^2 - w4^2 ) / ( w1- (w2-gamma) )/2.0
    

  ******************/


  w1 = tran_thermalarc_path_ene[0];
  w2 = tran_thermalarc_path_ene[1];
  w3 = tran_thermalarc_path_ene[2];
  w4 = tran_thermalarc_path_ene[3];
  gamma = tran_thermalarc_path_ene[4];

  b = ( sqr( w1 ) +sqr(w3) - sqr( w2-gamma ) - sqr( w4 ) )/
    ( w1- w2+ gamma )*0.5;

  a = sqrt( sqr(w1-b) + sqr(w3) );


  /* comment 
       (1,0) is 0 degree
       (sqrt(2)/2,sqrt(2)/2) is 45 degree 
       (-1,0) is 180 degree */

  th1 = acos( (w1-b)/a );
  th2 = asin( w4/a );

  dth = (th2-th1)/tran_thermalarc_path_div[0];

  /*   int_C a(x) nF(x) dx  
      = int_C a(x(t)) (nF(x(t)) dx/dt) dt ,
      here calculate (nF(x(t)) dx/dt) */ 

  for (i=0;i<tran_thermalarc_path_div[0] ;i++) {
    double x,th,d;
    dcomplex z,z2,z1,fermi,dxdth;

    x=0.5+i;
    th = th1 + dth*x;
    z.r = a*cos(th)+b;
    z.i = a*sin(th);
    tran_omega_scf[itotal].r = z.r;
    tran_omega_scf[itotal].i = z.i;

    z2.r = beta*(z.r-chem); z2.i = beta*(z.i);
    z_exp_inline(z2,z1);
    z1.r += 1.0;   /* z1 = exp(beta(z-chem))+1  */
    d = sqr(z1.r)+ sqr(z1.i);
    fermi.r = z1.r/d; fermi.i = -z1.i/d;   /* fermi = 1/ (nf(beta(z-chem))+1 ) */
    dxdth.r= -a *sin(th)* dth; dxdth.i = a* cos(th)* dth; 
    z_mul_inline(fermi, dxdth, z1);  
     
    tran_omega_weight_scf[itotal].r =  z1.r*factor;
    tran_omega_weight_scf[itotal].i =  z1.i*factor;
    tran_integ_method_scf[itotal]=1;

    itotal++;
  }

  /*
    (2,4)-(gamma,0)  -> (2,4) 
  */

  { 
    dcomplex w_int; 
  
    w_int.r = gamma/w2; w_int.i =0.0;
    for (i=0;i<tran_thermalarc_path_div[1];i++) {
      double x,d;
      dcomplex z,z1,z2,fermi;
      
      x=0.5+i;
      z.r = w2-gamma+w_int.r*x; z.i = w4;
      tran_omega_scf[itotal].r = z.r;
      tran_omega_scf[itotal].i = z.i;

      z2.r = beta*(z.r-chem); z2.i = beta*(z.i);
      z_exp_inline(z2,z1);   /* exp(beta(z-chem)) */
      z1.r += 1.0;   /* z1 = exp(beta(z-chem))+1  */
      d = sqr(z1.r)+ sqr(z1.i);
      fermi.r = z1.r/d; fermi.i = -z1.i/d;

      z_mul_inline( w_int, fermi, z2);
      
      tran_omega_weight_scf[itotal].r = -z2.r*factor;
      tran_omega_weight_scf[itotal].i = -z2.i*factor;
      tran_integ_method_scf[itotal]=1;
      
      itotal++;
    }
  }

  /* contour integral */

  for (i=0;i<n_max;i++)  {
    double th;
    th = (2.0*(double)i+1.0)*PI*tran_temperature;
    tran_omega_scf[itotal].r = chem; 
    tran_omega_scf[itotal].i= th;
    tran_omega_weight_scf[itotal].r=0; 
    tran_omega_weight_scf[itotal].i=-factor*2.0*PI*tran_temperature;
    tran_integ_method_scf[itotal]=1;
    itotal++;
  }


  /* applied bias */

  if (tran_bias_apply) {
    dcomplex w_int;
    factor = 1.0/(2.0*PI);
    w_int.r = (fabs(ChemP_e[0]-ChemP_e[1])+ tran_thermalarc_path_bias_expandenergy )/ tran_thermalarc_path_bias_div;
    w_int.i=0;
    for (i=0;i<tran_thermalarc_path_bias_div;i++) {
      double x;
      x=0.5+i;
      tran_omega_scf[itotal].r = chem+ w_int.r*x;
      tran_omega_scf[itotal].i = w2;
      tran_omega_weight_scf[itotal].r=  w_int.i*factor;
      tran_omega_weight_scf[itotal].i= -w_int.r*factor;
      
      tran_integ_method_scf[itotal]=2;
      itotal++;
      
    }
  }

}






/* Taisuke Ozaki Copyright (C) */

void TRAN_Set_IntegPath_GaussHG( double kBvalue, double Electronic_Temperature )
{
  int p;
  double beta,R;

   R = 1.0e+12;

  if ( tran_bias_apply ){
    if (TRAN_Kspace_grid2==1 && TRAN_Kspace_grid3==1)
      tran_omega_n_scf = 2*(tran_num_poles + 1);
    else 
      tran_omega_n_scf = 4*(tran_num_poles + 1);
  }
  else{ 
    if (TRAN_Kspace_grid2==1 && TRAN_Kspace_grid3==1)
      tran_omega_n_scf = tran_num_poles + 1;
    else 
      tran_omega_n_scf = 2*(tran_num_poles + 1);
  }

  /* allocation */
  tran_omega_scf = (dcomplex*)malloc(sizeof(dcomplex)*tran_omega_n_scf);
  tran_omega_weight_scf=(dcomplex*)malloc(sizeof(dcomplex)*tran_omega_n_scf);
  tran_integ_method_scf=(int*)malloc(sizeof(int)*tran_omega_n_scf);

  zero_cfrac( tran_num_poles, tran_omega_scf, tran_omega_weight_scf ); 
  beta = 1.0/kBvalue/Electronic_Temperature;

  /*
  for (p=0; p<tran_num_poles; p++){
    printf("p=%3d zp.r=%15.12f zp.i=%15.12f Rp.r=%15.12f Rp.i=%15.12f\n",
            p, tran_zp[p].r, tran_zp[p].i, tran_Rp[p].r, tran_Rp[p].i ); 
  }
  */

  /* bias case */

  if ( tran_bias_apply ){

    /*************************************
     (1) finite bias and only Gamma point

           p=0,tran_num_poles-1  

              poles for mu_L

           p=tran_num_poles  

              iR for mu_L

           p=tran_num_poles+1,2*tran_num_poles  

              poles for mu_R
         
           p=2*tran_num_poles+1  

              iR for mu_R
    *************************************/

    if (TRAN_Kspace_grid2==1 && TRAN_Kspace_grid3==1){

      /* contribution of poles */

      for (p=0; p<tran_num_poles; p++){
	tran_omega_scf[p].r = ChemP_e[0];
	tran_omega_scf[p].i = tran_omega_scf[p].i/beta;
	tran_omega_weight_scf[p].r = -2.0*tran_omega_weight_scf[p].r/beta;
	tran_omega_weight_scf[p].i = 0.0;
	tran_integ_method_scf[p] = 3; 

	tran_omega_scf[tran_num_poles+p+1].r = ChemP_e[1];
	tran_omega_scf[tran_num_poles+p+1].i = tran_omega_scf[p].i;
	tran_omega_weight_scf[tran_num_poles+p+1].r = tran_omega_weight_scf[p].r;
	tran_omega_weight_scf[tran_num_poles+p+1].i =  0.0;
	tran_integ_method_scf[tran_num_poles+p+1] = 4; 
      }

      /* contribution of the half circle contour integral */
  
      tran_omega_scf[tran_num_poles].r = 0.0;
      tran_omega_scf[tran_num_poles].i = R;
      tran_omega_weight_scf[tran_num_poles].r = 0.0;
      tran_omega_weight_scf[tran_num_poles].i = 0.5*R; 
      tran_integ_method_scf[tran_num_poles] = 3;

      tran_omega_scf[2*tran_num_poles+1].r = 0.0;
      tran_omega_scf[2*tran_num_poles+1].i = R;
      tran_omega_weight_scf[2*tran_num_poles+1].r = 0.0;
      tran_omega_weight_scf[2*tran_num_poles+1].i = 0.5*R; 
      tran_integ_method_scf[2*tran_num_poles+1] = 4;

    }

    /*************************************
     (2) finite bias and lots of k-points

           p=0,tran_num_poles-1  

              poles for retarded with mu_L

           p=tran_num_poles  

              iR for retarded with mu_L

           p=tran_num_poles+1,2*tran_num_poles  

              poles for advanced with mu_L
         
           p=2*tran_num_poles+1  

              iR for advanced with mu_L


           p=2*tran_num_poles+2,3*tran_num_poles+1  

              poles for retarded with mu_R

           p=3*tran_num_poles + 2  

              iR for retarded with mu_L

           p=3*tran_num_poles+3,4*tran_num_poles+2  

              poles for advanced with mu_L
         
           p=4*tran_num_poles+3  

              iR for advanced with mu_L

    *************************************/

    else{

      /* contribution of poles */

      for (p=0; p<tran_num_poles; p++){

        /* retarded with mu_L */

	tran_omega_scf[p].r = ChemP_e[0];
	tran_omega_scf[p].i = tran_omega_scf[p].i/beta;
	tran_omega_weight_scf[p].r = -tran_omega_weight_scf[p].r/beta;
	tran_omega_weight_scf[p].i = 0.0;
	tran_integ_method_scf[p] = 3;

        /* advanced with mu_L */

	tran_omega_scf[tran_num_poles+p+1].r = ChemP_e[0];
	tran_omega_scf[tran_num_poles+p+1].i = -tran_omega_scf[p].i;
	tran_omega_weight_scf[tran_num_poles+p+1].r = tran_omega_weight_scf[p].r;
	tran_omega_weight_scf[tran_num_poles+p+1].i = 0.0;
	tran_integ_method_scf[tran_num_poles+p+1] = 3; 

        /* retarded with mu_R */

	tran_omega_scf[2*tran_num_poles+p+2].r = ChemP_e[1];
	tran_omega_scf[2*tran_num_poles+p+2].i = tran_omega_scf[p].i;
	tran_omega_weight_scf[2*tran_num_poles+p+2].r = tran_omega_weight_scf[p].r;
	tran_omega_weight_scf[2*tran_num_poles+p+2].i = 0.0;
	tran_integ_method_scf[2*tran_num_poles+p+2] = 4; 

        /* advanced with mu_R */

	tran_omega_scf[3*tran_num_poles+p+3].r = ChemP_e[1];
	tran_omega_scf[3*tran_num_poles+p+3].i = -tran_omega_scf[p].i;
	tran_omega_weight_scf[3*tran_num_poles+p+3].r = tran_omega_weight_scf[p].r;
	tran_omega_weight_scf[3*tran_num_poles+p+3].i = 0.0;
	tran_integ_method_scf[3*tran_num_poles+p+3] = 4; 

      }

      /* contribution of the half circle contour integral */
  
      tran_omega_scf[tran_num_poles].r = 0.0;
      tran_omega_scf[tran_num_poles].i = R;
      tran_omega_weight_scf[tran_num_poles].r = 0.0;
      tran_omega_weight_scf[tran_num_poles].i = 0.25*R; 
      tran_integ_method_scf[tran_num_poles] = 3;

      tran_omega_scf[2*tran_num_poles+1].r = 0.0;
      tran_omega_scf[2*tran_num_poles+1].i = R;
      tran_omega_weight_scf[2*tran_num_poles+1].r = 0.0;
      tran_omega_weight_scf[2*tran_num_poles+1].i = 0.25*R; 
      tran_integ_method_scf[2*tran_num_poles+1] = 3;

      tran_omega_scf[3*tran_num_poles+2].r = 0.0;
      tran_omega_scf[3*tran_num_poles+2].i = R;
      tran_omega_weight_scf[3*tran_num_poles+2].r = 0.0;
      tran_omega_weight_scf[3*tran_num_poles+2].i = 0.25*R; 
      tran_integ_method_scf[3*tran_num_poles+2] = 4;

      tran_omega_scf[4*tran_num_poles+3].r = 0.0;
      tran_omega_scf[4*tran_num_poles+3].i = R;
      tran_omega_weight_scf[4*tran_num_poles+3].r = 0.0;
      tran_omega_weight_scf[4*tran_num_poles+3].i = 0.25*R; 
      tran_integ_method_scf[4*tran_num_poles+3] = 4;
    } 
  }

  /* zero-bias case */

  else {

    /*************************************
     (3) zero-bias and only Gamma point
    *************************************/

    if (TRAN_Kspace_grid2==1 && TRAN_Kspace_grid3==1){

      /* contribution of poles */

      for (p=0; p<tran_num_poles; p++){
	tran_omega_scf[p].r = ChemP_e[0];
	tran_omega_scf[p].i = tran_omega_scf[p].i/beta;
	tran_omega_weight_scf[p].r = -2.0*tran_omega_weight_scf[p].r/beta;
	tran_omega_weight_scf[p].i =  0.0;
      }

      /* contribution of the half circle contour integral */

      tran_omega_scf[tran_num_poles].r = 0.0;
      tran_omega_scf[tran_num_poles].i = R;
      tran_omega_weight_scf[tran_num_poles].r = 0.0; 
      tran_omega_weight_scf[tran_num_poles].i = 0.5*R;

      /* set the integral method in 1 */

      for (p=0; p<=tran_num_poles; p++){
	tran_integ_method_scf[p] = 1; 
      }
    }

    /*************************************
     (4) zero-bias and lots of k-points

           p=0,tran_num_poles-1  

              poles for retarded

           p=tran_num_poles  

              iR for retarded

           p=tran_num_poles+1,2*tran_num_poles  

              poles for advanced
         
           p=2*tran_num_poles+1  

              iR for advanced
    *************************************/

    else{
      
      /* contribution of poles */

      for (p=0; p<tran_num_poles; p++){

	tran_omega_scf[p].r = ChemP_e[0];
	tran_omega_scf[p].i = tran_omega_scf[p].i/beta;
	tran_omega_weight_scf[p].r = -tran_omega_weight_scf[p].r/beta;
	tran_omega_weight_scf[p].i = 0.0;
	tran_integ_method_scf[p] = 1; 

	tran_omega_scf[tran_num_poles+p+1].r = ChemP_e[0];
	tran_omega_scf[tran_num_poles+p+1].i =-tran_omega_scf[p].i;
	tran_omega_weight_scf[tran_num_poles+p+1].r = tran_omega_weight_scf[p].r;
	tran_omega_weight_scf[tran_num_poles+p+1].i =  0.0;
	tran_integ_method_scf[tran_num_poles+p+1] = 1; 
      }

      /* contribution of the half circle contour integral */

      tran_omega_scf[tran_num_poles].r = 0.0;
      tran_omega_scf[tran_num_poles].i = R;
      tran_omega_weight_scf[tran_num_poles].r = 0.0;
      tran_omega_weight_scf[tran_num_poles].i = 0.25*R; 
      tran_integ_method_scf[tran_num_poles] = 1;

      tran_omega_scf[2*tran_num_poles+1].r = 0.0;
      tran_omega_scf[2*tran_num_poles+1].i = R;
      tran_omega_weight_scf[2*tran_num_poles+1].r = 0.0;
      tran_omega_weight_scf[2*tran_num_poles+1].i = 0.25*R; 
      tran_integ_method_scf[2*tran_num_poles+1] = 1;

    } 
  }  

  /*
  printf("ChemP_e[0]=%15.12f\n",ChemP_e[0]);
  printf("ChemP_e[1]=%15.12f\n",ChemP_e[1]);
  */

}










/* Taisuke Ozaki Copyright (C) */

void zero_cfrac( int n, dcomplex *zp, dcomplex *Rp ) 
{
  static int i,j,N;
  static double **a,**s,*eval,*resr,*resi;

  /* check input parameters */

  if (n<=0){
    printf("\ncould not find the number of zeros\n\n");
    MPI_Finalize();
    exit(0);
  }

  /* the total number of zeros including minus value */

  N = 2*n + 1;

  /* allocation of arrays */

  a = (double**)malloc(sizeof(double*)*(N+2));
  for (i=0; i<(N+2); i++){
    a[i] = (double*)malloc(sizeof(double)*(N+2));
  }

  s = (double**)malloc(sizeof(double*)*(N+2));
  for (i=0; i<(N+2); i++){
    s[i] = (double*)malloc(sizeof(double)*(N+2));
  }

  eval = (double*)malloc(sizeof(double)*(n+3));
  resr = (double*)malloc(sizeof(double)*(n+3));
  resi = (double*)malloc(sizeof(double)*(n+3));

  /* initialize arrays */

  for (i=0; i<(N+2); i++){
    for (j=0; j<(N+2); j++){
      a[i][j] = 0.0;
      s[i][j] = 0.0;
    }
  }

  /* set matrix elements */

  s[2][1] =  1.0;
  s[2][2] = -0.5;

  for (i=3; i<=N; i++){
     s[i][i-1] =  -0.5;
     s[i-1][i] =   0.5;
  }

  a[1][1] = -2.0;
  a[1][2] =  1.0;
  a[2][2] = -1.0;

  for (i=3; i<=N; i++){
    a[i][i] = -(2.0*(double)i - 3.0);
  }

  /* diagonalization */

  Eigen_DGGEVX( N, a, s, eval, resr, resi );

  for (i=0; i<n; i++){
    zp[i].r = 0.0;
    zp[i].i = eval[i+1];
    Rp[i].r = resr[i+1];
    Rp[i].i = 0.0;
  }

  /* print result */

  /*
  for (i=1; i<=n; i++){
    printf("i=%4d  eval=%18.14f   resr=%18.15f resi=%18.15f\n",i,eval[i],resr[i],resi[i]);
  }
  */

  /*
  for (i=1; i<=n; i++){
    printf("%4d 0.0 %18.14f\n",i,eval[i]);
  }

  MPI_Finalize();
  exit(0);
  */

  /* free of arrays */

  for (i=0; i<(N+2); i++){
    free(a[i]);
  }
  free(a);

  for (i=0; i<(N+2); i++){
    free(s[i]);
  }
  free(s);

  free(eval);
  free(resr);
  free(resi);
}


/* Taisuke Ozaki Copyright (C) */

void dqsort(long n, double *a, double *b)
{
  int i;
  dlists *AB;

  AB = (dlists*)malloc(sizeof(dlists)*n);

  for (i=0; i<n; i++){
    AB[i].a = a[i+1];     
    AB[i].b = b[i+1];
  }

  qsort(AB, n, sizeof(dlists), (int(*)(const void*, const void*))dlists_cmp);

  for (i=0; i<n; i++){
    a[i+1] = AB[i].a;
    b[i+1] = AB[i].b;
  }

  free(AB);
}

/* Taisuke Ozaki Copyright (C) */

int dlists_cmp(const dlists *x, const dlists *y)
{
  return (x->a < y->a ? -1 :
          y->a < x->a ?  1 : 0);
}


/* Taisuke Ozaki Copyright (C) */

void Eigen_DGGEVX( int n, double **a, double **s, double *eval, double *resr, double *resi )
{
  static int i,j,k,l,num;

  static char balanc = 'N';
  static char jobvl = 'V';
  static char jobvr = 'V';
  static char sense = 'B';
  static double *A;
  static double *b;
  static double *alphar;
  static double *alphai;
  static double *beta;
  static double *vl;
  static double *vr;
  static double *lscale;
  static double *rscale;
  static double abnrm;
  static double bbnrm;
  static double *rconde;
  static double *rcondv;
  static double *work;
  static double *tmpvecr,*tmpveci;
  static INTEGER *iwork;
  static INTEGER lda,ldb,ldvl,ldvr,ilo,ihi;
  static INTEGER lwork,info;
  static logical *bwork; 
  static double sumr,sumi,tmpr,tmpi;
  static double *sortnum;

  lda = n;
  ldb = n;
  ldvl = n;
  ldvr = n;

  A = (double*)malloc(sizeof(double)*n*n);
  b = (double*)malloc(sizeof(double)*n*n);
  alphar = (double*)malloc(sizeof(double)*n);
  alphai = (double*)malloc(sizeof(double)*n);
  beta = (double*)malloc(sizeof(double)*n);

  vl = (double*)malloc(sizeof(double)*n*n);
  vr = (double*)malloc(sizeof(double)*n*n);

  lscale = (double*)malloc(sizeof(double)*n);
  rscale = (double*)malloc(sizeof(double)*n);

  rconde = (double*)malloc(sizeof(double)*n);
  rcondv = (double*)malloc(sizeof(double)*n);

  lwork = 2*n*n + 12*n + 16;
  work = (double*)malloc(sizeof(double)*lwork);

  iwork = (INTEGER*)malloc(sizeof(INTEGER)*(n+6));
  bwork = (logical*)malloc(sizeof(logical)*n);

  tmpvecr = (double*)malloc(sizeof(double)*(n+2));
  tmpveci = (double*)malloc(sizeof(double)*(n+2));

  sortnum = (double*)malloc(sizeof(double)*(n+2));

  /* convert two dimensional arrays to one-dimensional arrays */

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
       A[j*n+i]= a[i+1][j+1];
       b[j*n+i]= s[i+1][j+1];
    }
  }

  /* call dggevx_() */

  F77_NAME(dggevx,DGGEVX)(
           &balanc, &jobvl, & jobvr, &sense, &n, A, &lda, b, &ldb,
           alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi,
           lscale, rscale, &abnrm, &bbnrm, rconde, rcondv, work,
           &lwork, iwork, bwork, &info );

  if (info!=0){
    printf("Errors in dggevx_() info=%2d\n",info);
  }

  /*
  for (i=0; i<n; i++){
    printf("i=%4d  %18.13f %18.13f %18.13f\n",i,alphar[i],alphai[i],beta[i]);
  }
  printf("\n");
  */

  num = 0;
  for (i=0; i<n; i++){

    if ( 1.0e-13<fabs(beta[i]) && 0.0<alphai[i]/beta[i] ){

      /* normalize the eigenvector */

      for (j=0; j<n; j++) {

        sumr = 0.0;
        sumi = 0.0;

        for (k=0; k<n; k++) {
          sumr += s[j+1][k+1]*vr[i*n    +k];
          sumi += s[j+1][k+1]*vr[(i+1)*n+k];
        }
        
        tmpvecr[j] = sumr;
        tmpveci[j] = sumi;
      }

      sumr = 0.0;
      sumi = 0.0;

      for (k=0; k<n; k++) {
        sumr += vl[i*n+k]*tmpvecr[k] + vl[(i+1)*n+k]*tmpveci[k];
        sumi += vl[i*n+k]*tmpveci[k] - vl[(i+1)*n+k]*tmpvecr[k];
      }

      /* calculate zero point and residue */

      eval[num+1] = alphai[i]/beta[i];
      tmpr =  vr[i*n]*vl[i*n] + vr[(i+1)*n]*vl[(i+1)*n];
      tmpi = -vr[i*n]*vl[(i+1)*n] + vr[(i+1)*n]*vl[i*n];
      resr[num+1] =  tmpi/sumi;
      resi[num+1] = -tmpr/sumi;

      num++;
    }
    else{
      /*
      printf("i=%4d  %18.13f %18.13f %18.13f\n",i+1,alphar[i],alphai[i],beta[i]);
      */
    }
  }

  /* check round-off error */

  for (i=1; i<=num; i++){
    if (1.0e-8<fabs(resi[i])){
      printf("Could not calculate zero points and residues due to round-off error\n");
      MPI_Finalize();
      exit(0);
    }
  }

  /* sorting */

  dqsort(num,eval,resr);

  /* free arraies */

  free(A);
  free(b);
  free(alphar);
  free(alphai);
  free(beta);

  free(vl);
  free(vr);

  free(lscale);
  free(rscale);

  free(rconde);
  free(rcondv);

  free(work);

  free(iwork);
  free(bwork);

  free(tmpvecr);
  free(tmpveci);
  free(sortnum);
}


















void TRAN_Set_IntegPath( double kBvalue, double Electronic_Temperature )
{

  if ( tran_integ_pathtype==1 ) {
    TRAN_Set_IntegPath_Square();
  }
  else if ( tran_integ_pathtype==10 ) {
    TRAN_Set_IntegPath_ThermalArc();
  }

  /* Taisuke Ozaki Copyright (C) */
  else if ( tran_integ_pathtype==3 ) {
    TRAN_Set_IntegPath_GaussHG( kBvalue, Electronic_Temperature );
  }
  else {
    printf("abort tran_integ_pathtype=%d, not supported\n",tran_integ_pathtype);
    exit(10);
  }
}


