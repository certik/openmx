#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_variables.h"
#include "tran_prototypes.h"

 
#ifdef MAX 
#undef MAX
#endif
#define MAX(a,b) ((a)>(b))?  (a):(b) 

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b))?  (a):(b)


/*
 *  output: TRAN_region[] is changed 
 *          TRAN_grid_bound[2];
 *
 * The region where potential change starts from TRAN_grid_bound[0]+1 to TRAN_grid_bound[1]-1.
 * The whole region is 0 to Ngrid1-1
 *
 * electrode is [0:TRAN_grid_bound[0]] and [TRAN_grid_bound[1]: Ngird1-1] 
 */
void TRAN_Calc_GridBound(
  MPI_Comm mpi_comm_level1,
  int atomnum,
  int *WhatSpecies,
  double *Spe_Atom_Cut1,
  int Ngrid1, 
  double *Grid_Origin, 
  double **Gxyz,
  double tv[4][4],
  double gtv[4][4],
  double Right_tv[4][4]
)
{
  int ct_AN, wanA;
  int i,j;
  int side;
  double rcutA;
  double tmp_l,tmp_r,tmp_e; 
  int idim;
  double R[4];
  int myid,numprocs;

  /* MPI */ 
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  idim = 1; /* a-axis */

  for (i=1; i<=3; i++) { R[i]= tv[idim][i]-Right_tv[idim][i]; }

  if (myid==Host_ID){
    printf("\n\n<TRAN_Calc_GridBound>\n\n");
  }

  /********************************************************** 
                          left side
  ***********************************************************/

  /* find max. position where no overlap to the left edge=0 */

  tmp_l = 0.0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    if ( TRAN_region[ct_AN]==2 ) {  /* left region */
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];
      if ( (Gxyz[ct_AN][1] - rcutA) < Grid_Origin[1] ) {
        TRAN_region[ct_AN] += 10;
        tmp_l = MAX(tmp_l, Gxyz[ct_AN][1] + rcutA);
      }
    }
  }


  if (myid==Host_ID){
    printf("Grid_Origin=%lf %lf %lf\n",Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]); 
    printf("tv1=%lf %lf %lf\n",tv[1][1], tv[1][2], tv[1][3]);
    printf("tv2=%lf %lf %lf\n",tv[2][1], tv[2][2], tv[2][3]);
    printf("tv3=%lf %lf %lf\n",tv[3][1], tv[3][2], tv[3][3]);
    printf("Gtv1, tmp_l, %lf %lf\n", gtv[idim][1],tmp_l);
  }


  /***************************************************
        grid point is Grid_Origin + gtv*integer 
         Grid_Origin + gtv*integer = tmp_l 
  ***************************************************/
 
  tmp_l = (tmp_l-Grid_Origin[1])/gtv[idim][1] ;
  TRAN_grid_bound[0]= (int)tmp_l;

  if (myid==Host_ID){
    printf("grid_bound (left) = %lf %d\n",tmp_l, TRAN_grid_bound[0]);
  }

  side=0;
  tmp_l = Grid_Origin[1] + gtv[idim][1]*TRAN_grid_bound[0];

  for (i=0; i<Ngrid1_e[side]; i++){
    tmp_e = Grid_Origin_e[side][1] + gtv_e[side][idim][1]*(double)i;
    if (tmp_e > tmp_l) {
      TRAN_grid_bound_e[side]= i;
      TRAN_grid_bound_diff[side]= tmp_e - tmp_l;
      break;
    }
  }

  if (myid==Host_ID){
    printf("grid_bound_e(left) ( > grid_bound(left) )=%d diff=%lf\n",
	   TRAN_grid_bound_e[side], TRAN_grid_bound_diff[side]);
  }

#if 0

  if (myid==Host_ID){

    side=0;
    printf("Grid_Origin_e[left] = %lf %lf %lf\n",
            Grid_Origin_e[side][1],Grid_Origin_e[side][2], Grid_Origin_e[side][3]);
    side=1;
    printf("Grid_Origin_e[right]= %lf %lf %lf\n",
            Grid_Origin_e[side][1],Grid_Origin_e[side][2], Grid_Origin_e[side][3]);

    side=0;
    printf("gtv_e[left][3] =%lf %lf %lf\n", gtv_e[side][3][1],  gtv_e[side][3][2],  gtv_e[side][3][3]) ;
    side=1;
    printf("gtv_e[right][3]=%lf %lf %lf\n", gtv_e[side][3][1],  gtv_e[side][3][2],  gtv_e[side][3][3]) ;

    i=1;
    printf("Gxyz[1]-Origin=%lf %lf %lf\n", Gxyz[i][1]-Grid_Origin[1], 
	   Gxyz[i][2]-Grid_Origin[2], Gxyz[i][3]-Grid_Origin[3]);

    printf("(Origin+tv[3])-Gxyz[atomnum]=%lf %lf %lf\n",
	   Grid_Origin[1]+tv[3][1]-Gxyz[atomnum][1],
	   Grid_Origin[2]+tv[3][2]-Gxyz[atomnum][2],
	   Grid_Origin[3]+tv[3][3]-Gxyz[atomnum][3]);
  }

  for (j=0;j<40;j++)  {
    side=0;
    tmp_l=  Grid_Origin[1]+ gtv[3][1]*j;

    for (i=0;i<Ngrid3_e[side];i++) {
      tmp_e = Grid_Origin_e[side][1] + gtv_e[side][3][1]*i;
      if (tmp_e >= tmp_l) {
	TRAN_grid_bound_e[side]= i;
	TRAN_grid_bound_diff[side]= tmp_e - tmp_l;
	break;
      }
    }

    if (myid==Host_ID){
      printf("grid_bound_e(left) ( > grid_bound(left) ) j=%d i=%d tmp_l=%lf tmp_e=%lf diff=%lf\n",
              j, TRAN_grid_bound_e[side], tmp_l,tmp_e, TRAN_grid_bound_diff[side] );
    }
  }

#endif

  /********************************************************** 
                        right side 
  ***********************************************************/

  /* find min. position where no overlap to the right edge=tv[3][1] */

  tmp_r=1.0e+50;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    if ( TRAN_region[ct_AN]==3 ) {  /* right region */
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];

      if ( (Gxyz[ct_AN][1] + rcutA)  > (Grid_Origin[1]+tv[idim][1]) ) {
	TRAN_region[ct_AN] += 10;
	tmp_r = MIN(tmp_r, Gxyz[ct_AN][1] - rcutA);
      }
    }
  }

  if (myid==Host_ID){
    printf("Gtv1, tmp_r, %lf %lf %lf\n", gtv[idim][1],tmp_r,Grid_Origin[1]+tv[idim][1]-tmp_r);
  }

  tmp_r = (tmp_r-Grid_Origin[1])/gtv[idim][1];
  TRAN_grid_bound[1]= (int) (tmp_r)+1;

  if (myid==Host_ID){
    printf("grid_bound (right) = %lf %d\n",tmp_r, TRAN_grid_bound[1]);
    printf("grid_bound(left,right) and  Ngrid1= %d %d %d\n",
           TRAN_grid_bound[0], TRAN_grid_bound[1], Ngrid1);

    printf("TRAN_region\n");
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      printf("%d ",TRAN_region[ct_AN]);
    }
    printf("\n\n");
  }

  side=1;
  tmp_r = Grid_Origin[1]+gtv[idim][1]*TRAN_grid_bound[1]; 

  for (i=Ngrid3_e[side]-1; i>=0;i--) {
#if 0
    tmp_e = Grid_Origin_e[side][1]+ (tv[3][1]- Right_tv[3][1]) + gtv_e[side][3][1]*i; 
#else
    tmp_e = Grid_Origin_e[side][1]+ R[1] + gtv_e[side][idim][1]*i;

#endif
    if ( tmp_e < tmp_r ) {
      TRAN_grid_bound_e[side]=i;
      TRAN_grid_bound_diff[side] = tmp_e - tmp_r;
      break;
    }
  }


  if (myid==Host_ID){
    printf("grid_bound_e(right) ( < grid_bound(right) )=%d diff=%lf\n",
	     TRAN_grid_bound_e[side], TRAN_grid_bound_diff[side]);
  }


#if 0

  if (myid==Host_ID){
    printf("tv[3][1]- Right_tv[3][1]=%lf %lf %lf\n",tv[3][1]- Right_tv[3][1],
	    tv[3][2]- Right_tv[3][2], tv[3][3]- Right_tv[3][3]);
  }

  for (j=Ngrid3-1;j>=Ngrid3-40;j--)  {
    side=1;
    tmp_r=  Grid_Origin[1]+ gtv[3][1]*j;

    for (i=Ngrid3_e[side]-1;i>=0;i--) {
      tmp_e = Grid_Origin_e[side][1]+ (tv[3][1]- Right_tv[3][1]) + gtv_e[side][3][1]*i;

      if (tmp_e <= tmp_r) {
	TRAN_grid_bound_e[side]= i;
	TRAN_grid_bound_diff[side]= tmp_e - tmp_r;
	break;
      }
    }

    if (myid==Host_ID){
      printf("grid_bound_e(right) ( > grid_bound(right) ) j=%d i=%d tmp_r=%lf tmp_e=%lf diff=%lf\n",
	     j, TRAN_grid_bound_e[side], tmp_r, tmp_e, TRAN_grid_bound_diff[side]);
    }

  }

#endif





#if 0

  if (myid==Host_ID){
    printf("Origin_e=%lf %lf %lf\n",Grid_Origin_e[side][1],Grid_Origin_e[side][2],Grid_Origin_e[side][3]);
    printf("tv[3]=%lf %lf %lf\n",tv[3][1],tv[3][2],tv[3][3]);
    printf("Origin_e+tv=%lf %lf %lf\n",Grid_Origin_e[side][1]+tv[3][1],
	   Grid_Origin_e[side][2]+tv[3][2],
	   Grid_Origin_e[side][3]+tv[3][3]);
    printf("gtv_e[side][3]=%lf %lf %lf\n",gtv_e[side][3][1], gtv_e[side][3][2], gtv_e[side][3][3]);
  }

  side=1;
  for (i=0;i<=3;i++) {
    R[i]= Grid_Origin_e[side][1]+ (tv[3][i]- Right_tv[3][i]) + gtv_e[side][3][i]*Ngrid3_e[side];
  }

  if (myid==Host_ID){
    printf("R=%lf %lf %lf\n",R[1],R[2],R[3]);
  }

  side=1;
  for (i=0;i<=3;i++) {
    R[i]=  (tv[3][i]- Right_tv[3][i]) + gtv_e[side][3][i]*Ngrid3_e[side];
  }

  if (myid==Host_ID){
    printf("dR=%lf %lf %lf\n",R[1],R[2],R[3]);
  }

#endif


}


