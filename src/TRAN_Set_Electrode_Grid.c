#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef fftw2 
#include <fftw.h>
#else
#include <fftw3.h> 
#endif


#include "tran_variables.h"
#include "tran_prototypes.h" 

#define print_stdout 0

#ifdef MAX
#undef MAX
#endif

#define MAX(x,y)  ( (x)>(y) )? (x):(y)  

double Dot_Product(double a[4], double b[4]);


static void  TRAN_FFT_Electrode_Grid(MPI_Comm comm1, int isign);





/*
 *  interpolate values on grids 
 *
 *   input
 *             Density_Grid_e
 *             dVHart_Grid_e
 * 
 *   
 *   output 
 *               double **ElectrodeDensity_Grid, (interpolated values)
 *               double **ElectrodedVHart_Grid,  (interpolated values)
 */


/*
 *  This is also used for Vtot_Grid 
 */


/* implicit input  from trans_variables.h 
 *      double  **Density_Grid_e
 *      double  **dVHart_Grid_e
 */
void TRAN_Set_Electrode_Grid(
			     MPI_Comm comm1,
			     double *Grid_Origin,   /* origin of the grid */
			     double tv[4][4],       /* unit vector  of the cell*/
			     double Left_tv[4][4],  /* unit vector  left */
			     double Right_tv[4][4], /* unit vector right */
			     double gtv[4][4],      /* unit vector  of the grid point, which is gtv*integer */
			     int Ngrid1,
			     int Ngrid2,
			     int Ngrid3,            /* # of c grid points */
			     double *Density_Grid   /* work */
			     )
{
  double offset[4];
  double R[4];
  int    l1[2];
  int    i,j,k;
  int    side;
  int spin;
  int n1,n2,n3;
  int id;
  double dx[4],dx_e[4];
  int idim=1;

  double chargeeps=1.0e-12;
  int myid;

  MPI_Comm_rank(comm1, &myid);

  if (myid==Host_ID){
    printf("<TRAN_Set_Electrode_Grid>\n");
  }

  if (print_stdout) printf("interpolation start\n");

  /* allocation */
  if (print_stdout){
    printf("%d %d %d %d %d\n",Ngrid1,Ngrid2,Ngrid3,TRAN_grid_bound[0], TRAN_grid_bound[1]);
  }

  for (side=0;side<2;side++) {
    ElectrodeDensity_Grid[side] = (double**)malloc(sizeof(double*)*(SpinP_switch_e[side]+1));
  }

  side=0;
  l1[0]= 0;
  l1[1]= TRAN_grid_bound[0];
  for (spin=0;spin<=SpinP_switch_e[side];spin++) {
    ElectrodeDensity_Grid[side][spin]= 
      (double*)malloc(sizeof(double)*Ngrid3*Ngrid2*(l1[1]-l1[0]+1));
  }
  ElectrodedVHart_Grid[side]=(double*)malloc(sizeof(double)*Ngrid3*Ngrid2*(l1[1]-l1[0]+1));

  side=1;
  l1[0]= TRAN_grid_bound[1];
  l1[1]= Ngrid1-1;
  for (spin=0;spin<=SpinP_switch_e[side];spin++) {
    ElectrodeDensity_Grid[side][spin]=
      (double*)malloc(sizeof(double)*Ngrid3*Ngrid2*(l1[1]-l1[0]+1));
  }
  ElectrodedVHart_Grid[side]=(double*)malloc(sizeof(double)*Ngrid3*Ngrid2*(l1[1]-l1[0]+1));


  /**********************************************
                      left side 
  **********************************************/

  side=0;

  /* find the fist atom of L */

  for (i=1;i<=atomnum;i++) {
    if ( TRAN_region[i]%10==2 && TRAN_Original_Id[i]==1){ 
      id = i; 
      break;
    }
  }

  if (print_stdout){
    printf("the first atom of L = %d\n",id);
  }

#if 1
  for (i=1;i<=3;i++) {
    dx_e[i] = Gxyz_e[side][1][i] - Grid_Origin_e[side][i];
  }

  for (i=1;i<=3;i++) {
    dx[i] = Gxyz[id][i] - Grid_Origin[i];
  }

  if (print_stdout){
    printf("shift of the cell(L)=%lf %lf %lf\n", dx[1]-dx_e[1],dx[2]-dx_e[2],dx[3]-dx_e[3]);
  }

  for (i=1;i<=3;i++) {
    offset[i]= -dx[i] + dx_e[i];
  }

  if (print_stdout){
    printf("shift of the cell(L)=%lf %lf %lf\n", offset[1], offset[2], offset[3]);
  }

#else
  for (i=1;i<=3;i++) {
    offset[i]= 0.0;
  }

  if (print_stdout){
    printf("shift of the cell(L)=%lf %lf %lf\n", offset[1], offset[2], offset[3]);
  }

#endif

  l1[0]= 0;
  l1[1]= TRAN_grid_bound[0];

  if (print_stdout){
    printf("left boundary %d %d\n",l1[0],l1[1]);
    printf("interpolation %d %d %d -> %d:%d %d %d\n",
	   Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side],
	   l1[0], l1[1],Ngrid2,Ngrid3);
  }

  for (spin=0;spin<=SpinP_switch_e[side]; spin++) {

    TRAN_Interp_ElectrodeDensity_Grid( comm1,
				       offset, Density_Grid_e[side][spin],
				       Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side],
				       tv_e[side], gtv_e[side],
				       Density_Grid,Ngrid1,Ngrid2,Ngrid3, l1, tv,gtv);
#define grid_idx(i,j,k) ( ((i)-l1[0])*Ngrid2*Ngrid3+(j)*Ngrid3 + (k) ) 
      
    for (i=l1[0];i<=l1[1];i++){
      for (j=0;j<Ngrid2;j++) {
	for (k=0; k<Ngrid3;k++) {
	  ElectrodeDensity_Grid[side][spin][ grid_idx(i,j,k) ] =
	    Density_Grid[ grid_idx(i,j,k) ];
	}
      }
    }

    for (i=0;i<Ngrid3*Ngrid2*(l1[1]-l1[0]+1); i++) {
      if ( ElectrodeDensity_Grid[side][spin][i] < chargeeps ) {
	ElectrodeDensity_Grid[side][spin][i] = 0.0;
      }
    }


#ifdef DEBUG
#undef DEBUG
#endif


#ifdef DEBUG
    { char name[100];int i; double R[4]; for (i=1;i<=3;i++) R[i]=0.0; 
    sprintf(name,"Density_Grid_el.%d",myid);
    TRAN_Print_Grid(name, Grid_Origin_e[side], gtv_e[side], 
		    Ngrid1_e[side],  Ngrid2_e[side],0,  Ngrid3_e[side]-1 ,
		    R, Density_Grid_e[side][0] ); /* spin=0*/
    sprintf(name,"ElectrodeDensity_Gridl.%d",myid);
    TRAN_Print_Grid(name,Grid_Origin,gtv,
		    l1[1]-l1[0]+1, Ngrid2, 0,Ngrid3-1,
		    R, ElectrodeDensity_Grid[side][0]); /* spin=0 */
    }
       
#endif



  }

  TRAN_Interp_ElectrodeDensity_Grid(comm1,
				    offset, dVHart_Grid_e[side],
				    Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side],
				    tv_e[side], gtv_e[side],
				    Density_Grid,Ngrid1,Ngrid2,Ngrid3,l1, tv,gtv);

  for (i=l1[0];i<=l1[1];i++){
    for (j=0;j<Ngrid2;j++) {
      for (k=0; k<Ngrid3;k++) {
        ElectrodedVHart_Grid[side][ grid_idx(i,j,k) ] =
          Density_Grid[ grid_idx(i,j,k) ]; 
      }
    }
  }

#ifdef DEBUG
  { char name[100]; int i; double R[4]; for (i=1;i<=3;i++) R[i]=0.0;
  sprintf(name,"dVHart_Grid_el.%d", myid);
  TRAN_Print_Grid(name,Grid_Origin_e[side], gtv_e[side],
		  Ngrid1_e[side],Ngrid2_e[side], 0,Ngrid3_e[side]-1,
		  R, dVHart_Grid_e[side]);
  sprintf(name,"ElectrodedVHart_Gridl.%d", myid);
  TRAN_Print_Grid(name,Grid_Origin, gtv,
		  l1[1]-l1[0]+1,Ngrid2,0,Ngrid3-1,
		  R, ElectrodedVHart_Grid[side] );
  }
#endif


  /**********************************************
                    right side 
  **********************************************/

  side=1;

  /* find the fist atom of R */
  for (i=1;i<=atomnum;i++) {
    if ( TRAN_region[i]%10== 3  && TRAN_Original_Id[i]==1 ) {
      id = i;
      break;
    }
  }

  if (print_stdout){
    printf("the first atom of R is =%d\n",id);
  }

  l1[0]= TRAN_grid_bound[1];
  l1[1]= Ngrid1-1;

  if (print_stdout){
    printf("right boundary %d %d\n",l1[0],l1[1]);
  }


  i=1;
  for (j=1;j<=3;j++) {
    R[j]=tv[i][j]-1.0*Right_tv[i][j];
  }
  /* Right cell starts at Origin+R 
   * Density at Origin+R = Density_e at Origin_e
   */
#if 1
  for (i=1;i<=3;i++) {
    dx_e[i]= Gxyz_e[side][1][i] - Grid_Origin_e[side][i];
  }
  for (i=1;i<=3;i++) {
    dx[i] = Gxyz[id][i] - (Grid_Origin[i]+R[i]); 
  }

  if (print_stdout){ 
    printf("shift of the cell(R)=%lf %lf %lf\n", dx[1]-dx_e[1],dx[2]-dx_e[2],dx[3]-dx_e[3]);
  }

  for (i=1;i<=3;i++) {
    offset[i]= -dx[i]+dx_e[i]-R[i]+gtv[idim][i]*l1[0]; 
  }

  if (print_stdout){
    printf("shift of the cell(R)=%lf %lf %lf\n", offset[1], offset[2],offset[3]);
  }

#else
  for (i=1;i<=3;i++) {
    offset[i]=  gtv[3][i]*l3[0]-R[i]; 
  }
  offset[3]= sqrt( Dot_Product(offset,offset) );
  offset[1]= 0.0;
  offset[2]= 0.0;

  if (print_stdout){
    printf("shift of the cell(R)=%lf %lf %lf\n", offset[1], offset[2],offset[3]);
  }

#endif

  l1[1] = l1[1]- l1[0];
  l1[0] = 0;

  if (print_stdout){
    printf("right boundary (shift) %d %d\n",l1[0],l1[1]);
  }

  for (spin=0;spin<= SpinP_switch_e[side]; spin++) {
    TRAN_Interp_ElectrodeDensity_Grid(comm1,
				      offset, Density_Grid_e[side][spin],
				      Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side],
				      tv_e[side], gtv_e[side],
				      Density_Grid,Ngrid1,Ngrid2,Ngrid3,l1, tv,gtv);

    /* ElectrodeDensity_Grid[2][ Ngrid1*Ngrid2*(l3[1]-l3[0]+1)] */


    for (i=l1[0];i<=l1[1];i++){  /* l1[0]=0 */
      for (j=0;j<Ngrid2;j++) {
	for (k=0; k<Ngrid3;k++) {
	  ElectrodeDensity_Grid[side][spin][ grid_idx(i,j,k) ] =
	    Density_Grid[ grid_idx(i,j,k) ]; 
	}
      }
    }
    for (i=0;i<Ngrid3*Ngrid2*(l1[1]-l1[0]+1); i++) {
      if (ElectrodeDensity_Grid[side][spin][i] < chargeeps ) {
	ElectrodeDensity_Grid[side][spin][i] =0.0;
      }
    }


  }

#ifdef DEBUG
  { char name[100];
  sprintf(name,"Density_Grid_er.%d",myid);
  TRAN_Print_Grid(name, Grid_Origin_e[side], gtv_e[side],
		  Ngrid1_e[side],Ngrid2_e[side],0,Ngrid3_e[side]-1,
		  R, Density_Grid_e[side][0]);
  sprintf(name,"ElectrodeDensity_Gridr.%d",myid);
  TRAN_Print_Grid(name, Grid_Origin,gtv,
		  l1[1]-l1[0]+1,Ngrid2,0, Ngrid3-1,
		  R, ElectrodeDensity_Grid[side][0]);
  }
#endif


  TRAN_Interp_ElectrodeDensity_Grid(comm1,
				    offset, dVHart_Grid_e[side],
				    Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side],
				    tv_e[side], gtv_e[side],
				    Density_Grid,Ngrid1,Ngrid2,Ngrid3, l1,tv,gtv);

  /* ElectrodeDensity_Grid[2][ Ngrid1*Ngrid2*(l3[1]-l3[0]+1)] */


  for (i=l1[0];i<=l1[1];i++){
    for (j=0;j<Ngrid2;j++) {
      for (k=0; k<Ngrid3;k++ ) {
	ElectrodedVHart_Grid[side][ grid_idx(i,j,k)]=
	  Density_Grid[ grid_idx(i,j,k) ]; 
      }
    }
  }

#ifdef DEBUG
  {
    char name[100];
    sprintf(name,"dVHart_Grid_er.%d",myid);
    TRAN_Print_Grid(name, Grid_Origin_e[side], gtv_e[side],
		    Ngrid1_e[side],Ngrid2_e[side],0,Ngrid3_e[side]-1,
		    R, dVHart_Grid_e[side] );
    sprintf(name,"ElectrodedVHart_Gridr.%d",myid);
    TRAN_Print_Grid(name,Grid_Origin,gtv,
		    l1[1]-l1[0]+1,Ngrid2,0, Ngrid3-1,
		    R, ElectrodedVHart_Grid[side] );
  }
#endif


  /*********************************************************
         FFT:   ElectrodedVHart_Grid(x,y,z) -> ElectrodedVHart_Grid(kx,ky,z)
  *********************************************************/
  TRAN_FFT_Electrode_Grid(comm1, -1);

  if (print_stdout){
    printf("interpolation end\n");
  }

}



/*********************************************************
         allocate ElectrodedVHart_Grid_c
         and 
         FFT  ElectrodedVHart_Grid(x,y,z) -> ElectrodedVHart_Grid_c(kx,ky,z)
*********************************************************/

static void  TRAN_FFT_Electrode_Grid(MPI_Comm comm1, int isign)
     /* #define grid_e_ref(i,j,k) ( (i)*Ngrid2*(l3[1]-l3[0]+1)+ (j)*(l3[1]-l3[0]+1) + (k)-l3[0] ) 
      *#define fft2d_ref(i,j)  ( (i)*Ngrid2+(j) ) 
      */
{

  int side;
  int i,j,k;
  int l1[2];

#ifdef fftw2
  fftwnd_plan p;
#else
  fftw_plan p;
#endif
  fftw_complex *in,*out;

  double factor;
  int myid;

  MPI_Comm_rank(comm1,&myid);

  if (print_stdout){
    printf("TRAN_FFT_Electrode_Grid in\n");
  }

  /* allocation 
   *   ElectrodedVHart_Grid_c is a global variable 
   *  do not free these 
   */
  l1[0]= 0;
  l1[1]= TRAN_grid_bound[0];
  ElectrodedVHart_Grid_c[0]=(dcomplex*)malloc(sizeof(dcomplex)*Ngrid3*Ngrid2*(l1[1]-l1[0]+1) );
  /* ElectrodedVHart_Grid_c = Vh_electrode(kx,ky, z=[0:TRAN_grid_bound[0]]), TRAN_grid_bound: integer */

  l1[0]= TRAN_grid_bound[1];
  l1[1]= Ngrid1-1;
  ElectrodedVHart_Grid_c[1]=(dcomplex*)malloc(sizeof(dcomplex)*Ngrid3*Ngrid2*(l1[1]-l1[0]+1) ); 
  /* ElectrodedVHart_Grid_c = Vh_electrode(kx,ky, z=[TRAN_grid_bound[1]:Ngrid3-1]), TRAN_grid_bound: integer */


  /* allocation for fft ,
   *  free these at last 
   */

#ifdef fftw2
  in  = (fftw_complex*)malloc(sizeof(fftw_complex)*Ngrid3*Ngrid2);
  out = (fftw_complex*)malloc(sizeof(fftw_complex)*Ngrid3*Ngrid2);
#else
  in  = fftw_malloc(sizeof(fftw_complex)*Ngrid3*Ngrid2); 
  out = fftw_malloc(sizeof(fftw_complex)*Ngrid3*Ngrid2); 
#endif
 

#ifdef fftw2
  p=fftw2d_create_plan(Ngrid2,Ngrid3,isign,FFTW_ESTIMATE);
#else
  p=fftw_plan_dft_2d(Ngrid2,Ngrid3,in,out,isign,FFTW_ESTIMATE);
#endif


  /* left side */

  side=0;
  l1[0]= 0;
  l1[1]= TRAN_grid_bound[0];
  factor = 1.0/( (double)(Ngrid2*Ngrid3) ) ;

#define grid_e_ref(i,j,k) ( ( (i)-l1[0])*Ngrid2*Ngrid3+(j)*Ngrid3+(k) )
#define fft2d_ref(j,k)  ( (j)*Ngrid3+ (k) )

  for (i=l1[0];i<=l1[1];i++) {
    for (j=0;j<Ngrid2;j++){
      for (k=0;k<Ngrid3;k++) {

#ifdef fftw2
	c_re(in[fft2d_ref(j,k)]) = ElectrodedVHart_Grid[side][grid_e_ref(i,j,k)]; 
	c_im(in[fft2d_ref(j,k)]) = 0.0;
#else
	in[fft2d_ref(j,k)][0]= ElectrodedVHart_Grid[side][grid_e_ref(i,j,k)]; 
	in[fft2d_ref(j,k)][1]= 0.0;
#endif

      }
    }

#ifdef fftw2
    fftwnd_one(p, in, out);
#else
    fftw_execute(p);
#endif

    for (j=0;j<Ngrid2;j++) {
      for (k=0;k<Ngrid3;k++) {

#ifdef fftw2
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].r = c_re(out[fft2d_ref(j,k)])*factor;
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].i = c_im(out[fft2d_ref(j,k)])*factor;
#else
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].r = out[fft2d_ref(j,k)][0]*factor;
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].i = out[fft2d_ref(j,k)][1]*factor;
#endif

      }
    }
  } /* k */

#ifdef DEBUG
  /*debug*/
  { char name[100];int i; double R[4]; for(i=1;i<=3;i++) R[i]=0.0; 
  sprintf(name,"ElectrodedVHart_Grid_c_lr.%d",myid);
  TRAN_Print_Grid_c(name,"ElectrodedVHart_Grid_c_li",
		    Grid_Origin, gtv, l1[1]-l1[0]+1,Ngrid2, 0, Ngrid3-1, R, ElectrodedVHart_Grid_c[side]);
  }
#endif


  /* right side */

  side=1;
  l1[0]= TRAN_grid_bound[1];
  l1[1]= Ngrid1-1;
  factor = 1.0/( (double) Ngrid2*Ngrid3 );

  for (i=l1[0];i<=l1[1];i++) {
    for (j=0;j<Ngrid2;j++) {
      for (k=0;k<Ngrid3;k++) {

#ifdef fftw2
	c_re(in[fft2d_ref(j,k)]) = ElectrodedVHart_Grid[side][grid_e_ref(i,j,k)];
        c_im(in[fft2d_ref(j,k)]) = 0.0;
#else
	in[fft2d_ref(j,k)][0]= ElectrodedVHart_Grid[side][grid_e_ref(i,j,k)];
        in[fft2d_ref(j,k)][1]= 0.0;
#endif

      }
    }

#ifdef fftw2
    fftwnd_one(p, in, out);
#else
    fftw_execute(p);
#endif

    for (j=0;j<Ngrid2;j++) {
      for (k=0;k<Ngrid3;k++) {

#ifdef fftw2
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].r = c_re(out[fft2d_ref(j,k)])*factor;
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].i = c_im(out[fft2d_ref(j,k)])*factor;
#else
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].r = out[fft2d_ref(j,k)][0]*factor;
        ElectrodedVHart_Grid_c[side][grid_e_ref(i,j,k)].i = out[fft2d_ref(j,k)][1]*factor;
#endif

      }
    }
  } /* k */

#ifdef DEBUG 
  /*debug*/
  { char name[100];int i; double R[4]; for(i=1;i<=3;i++) R[i]=0.0; 
  sprintf(name,"ElectrodedVHart_Grid_c_rr.%d",myid);
  TRAN_Print_Grid_c(name,"ElectrodedVHart_Grid_c_ri",
		    Grid_Origin, gtv,l1[1]-l1[0]+1,Ngrid2, 0, Ngrid3-1, R, ElectrodedVHart_Grid_c[side]);
  }
#endif


#ifdef fftw2
  fftwnd_destroy_plan(p);
#else
  fftw_destroy_plan(p);
#endif



#ifdef fftw2
  free(out);
  free(in);
#else
  fftw_free(out);
  fftw_free(in);
#endif

  if (print_stdout){
    printf("TRAN_FFT_Electrode_Grid out\n");
  }
}
