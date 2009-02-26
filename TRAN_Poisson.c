/**********************************************************************
  Poisson.c:

     Poisson.c is a subrutine to solve Poisson's equation using
     fast Fourier transformation.

  Log of Poisson.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#define  measure_time   0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"


#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"


#ifdef fftw2 
#include <fftw.h>
#else
#include <fftw3.h> 
#endif

/* input: Density_Grid (-ADensity_Grid )
   or input: ReV1, ImV1  ( which are FFT of Density_Grid (-ADensity_Grid ) )
   work: ReV2, ImV2
   output: dVHart_Grid = rho(G)/G^2 

   fft_charge_flag=1:   density_matrix mixing
   fft_charge_flag=0:   k-transformed Density_Grid mixing
*/



double TRAN_Poisson(
                    int fft_charge_flag,
		    double ***ReV1, double ***ImV1,
		    double ***ReV2, double ***ImV2)
     /*#define grid_l_ref(i,j,k) ( (i)*Ngrid2*(l3l[1]-l3l[0]+1)+(j)*(l3l[1]-l3l[0]+1)+(k)-l3l[0] )
      *#define grid_r_ref(i,j,k) ( (i)*Ngrid2*(l3r[1]-l3r[0]+1)+(j)*(l3r[1]-l3r[0]+1)+(k)-l3r[0] )
      */
{ 
  static int ct_AN,n1,n2,n3,k1,k2,k3,kk2,N3[4];
  static int Cwan,NO0,NO1,Rn,N,Hwan,i,j;
  static int nn0,nn1,nnn1,MN,MN0,MN1,MN2;
  static double time0;
  static int h_AN,Gh_AN,Rnh,spin,Nc,GNc,GN,Nh,Nog;
  static double tmp0,tmp1,sk1,sk2,sk3,tot_den;
  static double x,y,z,Gx,Gy,Gz,Eden[2],DenA,G2;
  static double ReTmp,ImTmp,c_coe,p_coe;
  static double *tmp_array0;
  static double *tmp_array1;
  static double TStime,TEtime;
  static int numprocs,myid,tag=999,ID;

  /* TRAN extension */
  fftw_complex  *fftin, *fftout;
  fftw_plan p; 
  int ix0,ixlp1; 
  int l1l[2], l1r[2];
  double *rev, *imv;


  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<TRAN_Poisson>  Poisson's equation using FFT...\n");

  dtime(&TStime);

#if 0
  {
    double R[4];
    R[1]=0.0; R[2]=0.0; R[3]=0.0;
    TRAN_Print_Grid_Cell0("zCharge_Density", Grid_Origin, gtv,
			  Ngrid1,Ngrid2, 0,Ngrid3-1, R,
			  My_Cell0, 
			  Density_Grid[0]);

  }
  {
    double R[4];
    double *v;
    int i;
    R[1]=0.0; R[2]=0.0; R[3]=0.0;
    v= (double*)malloc(sizeof(double)*Ngrid1*Ngrid2*Ngrid3);
    for (i=0;i< Ngrid1*Ngrid2*Ngrid3; i++) {
      v[i]= Density_Grid[0][i]+Density_Grid[1][i]-ADensity_Grid[i]; 
    }
    TRAN_Print_Grid_Cell0("zCharge_DensitySum", Grid_Origin, gtv,
			  Ngrid1,Ngrid2, 0,Ngrid3-1, R,
			  My_Cell0,
			  v);
    free(v);

  }
#endif

  /****************************************************
                 FFT of charge density 
  ****************************************************/

  if (fft_charge_flag==1) FFT_Density(0,ReV1,ImV1,ReV2,ImV2);

#ifdef DEBUG 
  /*debug*/
  { 
    char name[100];
    sprintf(name,"ReV20.%d",myid);
    TRAN_Print_Grid_Startv(name,  Ngrid1,My_NGrid2_Poisson, Ngrid3,Start_Grid2[myid], ReV2);
    sprintf(name,"ImV20.%d",myid);
    TRAN_Print_Grid_Startv(name,  Ngrid1,My_NGrid2_Poisson,Ngrid3,Start_Grid2[myid],  ImV2);
  }
#endif




  /****************************************************
                       4*PI/G2/N^3
  ****************************************************/

  /************************
     x -> k,  factor=1/N 
  ************************/

  tmp0 = 4.0*PI/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k1=0; k1<Ngrid1; k1++){

      if (k1<Ngrid1/2) sk1 = (double)k1;
      else             sk1 = (double)(k1 - Ngrid1);

      for (k3=0; k3<Ngrid3; k3++){

        if (k3<Ngrid3/2) sk3 = (double)k3;
        else             sk3 = (double)(k3 - Ngrid3);

        Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
        Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
        Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
        G2 = Gx*Gx + Gy*Gy + Gz*Gz;

        if (k1==0 && kk2==0 && k3==0){
          ReV2[k2][k1][k3] = 0.0;
          ImV2[k2][k1][k3] = 0.0;
        }
        else{
          ReV2[k2][k1][k3] = tmp0*ReV2[k2][k1][k3]/G2;
          ImV2[k2][k1][k3] = tmp0*ImV2[k2][k1][k3]/G2; 
        }
      }
    }
  }

  /****************************************************
   *   TRAN extension 
   ****************************************************/

  /* now (ReV2,ImV2) = Hartree potential(kx,ky,kz) */

  /* factor=1, in kx->x */

#ifdef fftw2
  fftin =(fftw_complex*)malloc(sizeof(fftw_complex)*Ngrid1);
  fftout=(fftw_complex*)malloc(sizeof(fftw_complex)*Ngrid1);
  p = fftw_create_plan(Ngrid1,1,FFTW_ESTIMATE);
#else
  fftin =fftw_malloc(sizeof(fftw_complex)*Ngrid1);
  fftout=fftw_malloc(sizeof(fftw_complex)*Ngrid1);
  p = fftw_plan_dft_1d(Ngrid1,fftin,fftout,1,FFTW_ESTIMATE);
#endif

  
  /*parallel global: k2  0:Ngrid2-1 */
  /*parallel local: k2  0:My_NGrid2_Poisson-1 */
  /*parallel local_to_global   k2 = kk2 + Start_Grid2[myid]  */
  for (k2=0;k2<My_NGrid2_Poisson ;k2++) {
    for (k3=0;k3<Ngrid3;k3++) {
      
      for (k1=0;k1<Ngrid1;k1++) {

#ifdef fftw2
        c_re(fftin[k1]) = ReV2[k2][k1][k3];
        c_im(fftin[k1]) = ImV2[k2][k1][k3];
#else
        fftin[k1][0] = ReV2[k2][k1][k3];
        fftin[k1][1] = ImV2[k2][k1][k3];
#endif

      }

#ifdef fftw2
      fftw_one(p, fftin, fftout); 
#else
      fftw_execute(p); 
#endif

      for (k1=0;k1<Ngrid1;k1++) {

#ifdef fftw2
        ReV2[k2][k1][k3] = c_re(fftout[k1]);
        ImV2[k2][k1][k3] = c_im(fftout[k1]);
#else
        ReV2[k2][k1][k3] = fftout[k1][0];
        ImV2[k2][k1][k3] = fftout[k1][1];
#endif

      }
    }
  }

  fftw_destroy_plan(p);

  /* now, (ReV2,ImV2)=Hartree potential(x,ky,kz) */

  /* boundary */

  l1l[0]=0;
  l1l[1]=TRAN_grid_bound[0]; 
  l1r[0]=TRAN_grid_bound[1];
  l1r[1]=Ngrid1-1; 

  ix0 =   TRAN_grid_bound[0];
  ixlp1 = TRAN_grid_bound[1];

  /*
   * V_e(x) given
   *
   * V_c(x) <= solved by rho/ (G^2) 
   *
   *  d^2/dx^2 V(x) = rho(x), linear equation
   *  d^2/dx^2 (VH(x) + dVH(x)) = ( rho(x) +  0 )
   *
   *  d^2/dx^2 VH(x) = rho(x) , solved via FFT
   *  d^2/dx^2 dVH(x) =0 with boundary condition, V_e(x)-V_c(x) at x=x_0 and x_(l+1)   
   *
   *
   * add correction,dVH(x), to the region x=[ix0+1:ixlp1-1]  
   */

  rev = (double*)malloc(sizeof(double)*Ngrid1);
  imv = (double*)malloc(sizeof(double)*Ngrid1);

#define grid_l_ref(i,j,k) ( ((i)-l1l[0])*Ngrid2*Ngrid3 + (j)*Ngrid3+ (k) )
#define grid_r_ref(i,j,k) ( ((i)-l1r[0])*Ngrid2*Ngrid3 + (j)*Ngrid3+ (k) )


  /*parallel global: kk2  0:Ngrid2-1 */
  /*parallel local: k2  0:My_NGrid2_Poisson-1 */
  /*parallel local_to_global   kk2 = k2 + Start_Grid2[myid]  */
  for (k2=0;k2<My_NGrid2_Poisson;k2++) {
    kk2 = k2 + Start_Grid2[myid];
    for (k3=0;k3<Ngrid3;k3++) {
      if (k3==0 && kk2==0) {  /* G_para==0 */
        for (k1=0;k1<Ngrid1;k1++) { rev[k1]=ReV2[k2][k1][k3]; imv[k1]=ImV2[k2][k1][k3]; }
	TRAN_Calc_VHartree_G0(
			      ElectrodedVHart_Grid_c[0][grid_l_ref(ix0,kk2,k3)],     /* (x,ky,kz), left edge */
			      ElectrodedVHart_Grid_c[1][grid_r_ref(ixlp1,kk2,k3)],   /* right edge           */
			      ReV2[k2][ix0][k3], ImV2[k2][ix0][k3],                  /* x0,ky,kz             */
			      ReV2[k2][ixlp1][k3], ImV2[k2][ixlp1][k3],              /* x_(l+1),ky,kz}       */
			      ix0, ixlp1, gtv[1],
			      rev, imv  /* output */
			      );
        for (k1=0;k1<Ngrid1;k1++) { ReV2[k2][k1][k3]=rev[k1]; ImV2[k2][k1][k3]=imv[k1]; }

      }
      else {         /* G_para .ne.0 */
        for (k1=0;k1<Ngrid1;k1++) { rev[k1]=ReV2[k2][k1][k3]; imv[k1]=ImV2[k2][k1][k3]; }

	TRAN_Calc_VHartree_Gnon0(
				 kk2,k3, ix0, ixlp1,gtv,
				 ReV2[k2][ix0][k3],ImV2[k2][ix0][k3],   /* x0,ky,kz */
				 ReV2[k2][ixlp1][k3], ImV2[k2][ixlp1][k3], /* x_(l+1),ky,kz */
				 ElectrodedVHart_Grid_c[0][grid_l_ref(ix0,kk2,k3)], /* (x,ky,kz) */
				 ElectrodedVHart_Grid_c[1][grid_r_ref(ixlp1,kk2,k3) ],
				 rev, imv  /* output */
				 );
        for (k1=0;k1<Ngrid1;k1++) { ReV2[k2][k1][k3]=rev[k1]; ImV2[k2][k1][k3]=imv[k1]; }
            
      }
    }
  }
 
  free(imv);
  free(rev);

#ifdef DEBUG 
  /*debug*/
  {
    char name[100];
    sprintf(name,"ReV2b.%d",myid);
    TRAN_Print_Grid_Startv(name,  Ngrid1,My_NGrid2_Poisson, Ngrid3,Start_Grid2[myid], ReV2);
    sprintf(name,"ImV2b.%d",myid);
    TRAN_Print_Grid_Startv(name,  Ngrid1,My_NGrid2_Poisson,Ngrid3,Start_Grid2[myid],  ImV2);
  }
#endif


  /* overwrite  ElectrodedVHart_Grid_c */
  /* (ReV2,ImV2), z=[0:TRAN_grid_bound[0]], [TRAN_grid_bound[1]:Ngrid3-1]  <= ElectrodedVHart_Grid_c */
  /* global 0:Ngrid2-1 */
  /* local  0:My_NGrid2_Poisson-1 */
  /* local_to_global   k2 = kk2 + Start_Grid2[myid]  */ 
  TRAN_Overwrite_V2(Ngrid1,Ngrid2,Ngrid3, TRAN_grid_bound,My_NGrid2_Poisson, Start_Grid2[myid], 
		    ElectrodedVHart_Grid_c,
		    ReV2, ImV2 ); /* output */

#ifdef DEBUG
  /*debug*/
  {  char name[100];
  sprintf(name,"ReV2c.%d",myid);
  TRAN_Print_Grid_Startv(name,  Ngrid1,My_NGrid2_Poisson,Ngrid3,Start_Grid2[myid], ReV2);
  sprintf(name,"ImV2c.%d",myid);
  TRAN_Print_Grid_Startv(name,  Ngrid1,My_NGrid2_Poisson,Ngrid3,Start_Grid2[myid], ImV2);
  }
#endif


  /*  now (ReV2, ImV2) is a smooth function */
  /* (ReV2,ImV2)=Hartree potential(x,ky,kz) with effects of boundary condition */

  {
    /* (kx,ky,z) -> (kx,ky,kz) */

#ifdef fftw2
    p = fftw_create_plan(Ngrid1, -1, FFTW_ESTIMATE);
#else
    p = fftw_plan_dft_1d(Ngrid1,fftin,fftout,-1,FFTW_ESTIMATE);
#endif

    tmp0=1.0/(double)Ngrid1;  
    for (k3=0;k3<Ngrid3;k3++) {
      /* global k2 0:Ngrid2-1 */
      /* local  k2  0:My_NGrid2_Poisson-1 */
      /* local_to_global k2 = k2+ Start_Grid2[myid] */
      for (k2=0;k2<My_NGrid2_Poisson;k2++) {

	for (k1=0;k1<Ngrid1;k1++) {

#ifdef fftw2
	  c_re(fftin[k1]) = ReV2[k2][k1][k3];
	  c_im(fftin[k1]) = ImV2[k2][k1][k3];
#else
	  fftin[k1][0]= ReV2[k2][k1][k3];
	  fftin[k1][1]= ImV2[k2][k1][k3];
#endif

	}
#if 0
        {
	  printf("fftin(%d)->", Ngrid1); 
	  for (k1=0;k1<Ngrid2;k1++) { printf("%lf ",fftin[k1][0]); } printf("\n");
        }
#endif

#ifdef fftw2
	fftw_one(p, fftin, fftout);
#else
	fftw_execute(p);
#endif

#if 0
        {
	  printf("fftout(%d)->", Ngrid1); 
	  for (k1=0;k1<Ngrid2;k1++) { printf("%lf ",fftout[k1][0]); } printf("\n");
        }
#endif

	for (k1=0;k1<Ngrid1;k1++) {

#ifdef fftw2
	  ReV2[k2][k1][k3]=c_re(fftout[k1])*tmp0;
	  ImV2[k2][k1][k3]=c_im(fftout[k1])*tmp0;
#else
	  ReV2[k2][k1][k3]=fftout[k1][0]*tmp0;
	  ImV2[k2][k1][k3]=fftout[k1][1]*tmp0;
#endif

	}
      }
    }

    fftw_destroy_plan(p);

  }

  /* (ReV2,ImV2)=Hartree potential(kx,ky,kz) with effects of boundary condition */

#if 0
  /*debug*/
  { int i; double R[4];  for (i=1;i<=3;i++) R[i]=0.0;
  TRAN_Print_Grid_v("ReV2d", Grid_Origin,gtv, Ngrid1,Ngrid2, 0,Ngrid3-1, R, ReV2);
  TRAN_Print_Grid_v("ImV2d", Grid_Origin,gtv, Ngrid1,Ngrid2, 0,Ngrid3-1, R, ImV2);
  }
#endif

#ifdef fftw2
  free(fftout);
  free(fftin);
#else
  fftw_free(fftout);
  fftw_free(fftin);
#endif


  /****************************************************
        find the Hartree potential in real space
  ****************************************************/

  Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,dVHart_Grid,dVHart_Grid); 

#ifdef DEBUG
  {
    int i; double R[4];  for (i=1;i<=3;i++) R[i]=0.0;
    TRAN_Print_Grid_Cell0( "dVHart_Grid", Grid_Origin,gtv, Ngrid1,Ngrid2, 0,Ngrid3-1, R,
			   My_Cell0, dVHart_Grid);
  }

#endif

  /* for time */

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}

