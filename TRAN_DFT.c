/**********************************************************************
  TRAN_DFT.c:

  TRAN_DFT.c is a subroutine to perform self-consistent calculations
  of a central region with left and right infinite leads based on
  a non-equilibrium Green's function method.

  Log of TRAN_DFT.c:

     11/Dec/2005   Released by H.Kino
     06/June/2005  Modified by T.Ozaki

***********************************************************************/

#define WRITE_DENSITY 0
#define MEASURE_TIME  0

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"
 
void dtime(double *);

void Make_Comm_Worlds(
   MPI_Comm MPI_Curret_Comm_WD,   
   int myid0,
   int numprocs0,
   int Num_Comm_World, 
   int *myworld1, 
   MPI_Comm *MPI_CommWD,     /* size: Num_Comm_World */
   int *NPROCS1_ID,          /* size: numprocs0 */
   int *Comm_World1,         /* size: numprocs0 */
   int *NPROCS1_WD,          /* size: Num_Comm_World */
   int *Comm_World_StartID   /* size: Num_Comm_World */
   );



int Get_OneD_HS_Col(int set_flag, double ****RH, double *H1, int *MP, 
                    int *order_GA, int *My_NZeros, int *is1, int *is2);

static void TRAN_Add_MAT(
    int mode,
    int NUM_c,
    dcomplex w_weight,
    dcomplex *v,  
    dcomplex *out
    );


static void TRAN_DFT_Kdependent(
			  /* input */
			  MPI_Comm comm1,
                          int parallel_mode,
                          int numprocs,
                          int myid,
			  int level_stdout,
			  int iter,
			  int SpinP_switch,
                          double k2,
                          double k3,
                          int k_op,
                          int *order_GA,
                          double **DM1,
                          double **H1,
                          double *S1,
			  double *****nh,  /* H */
			  double *****ImNL, /* not used, s-o coupling */
			  double ****CntOLP, 
			  int atomnum,
			  int Matomnum,
			  int *WhatSpecies,
			  int *Spe_Total_CNO,
			  int *FNAN,
			  int **natn, 
			  int **ncn,
			  int *M2G, 
			  int *G2ID, 
			  int **atv_ijk,
			  int *List_YOUSO,
			  /* output */
			  double *****CDM,  /* output, charge density */
			  double *****EDM,  /* not used */
			  double Eele0[2], double Eele1[2]); /* not used */






double TRAN_DFT(
		/* input */
                MPI_Comm comm1,
                int level_stdout,
		int iter, 
		int SpinP_switch,
		double *****nh,  /* H */
		double *****ImNL, /* not used, s-o coupling */
		double ****CntOLP, 
		int atomnum,
		int Matomnum,
		int *WhatSpecies,
		int *Spe_Total_CNO,
		int *FNAN,
		int **natn, 
		int **ncn,
		int *M2G, 
		int *G2ID,
                int *F_G2M, 
		int **atv_ijk,
		int *List_YOUSO,
		/* output */
		double *****CDM,  /* output, density matrix */
		double *****EDM,  /* not used */
                double ***TRAN_DecMulP, /* output, partial DecMulP */
		double Eele0[2], double Eele1[2], 
                double ChemP_e0[2]) 
{
  int numprocs0,myid0,ID;
  int i2,i3,k_op,l1,l2,l3,RnB;
  int k,E_knum,S_knum,T_knum,num_kloop0;
  int parallel_mode,kloop,kloop0;
  int i,j,spin,MA_AN,GA_AN,wanA,tnoA;
  int LB_AN,GB_AN,wanB,tnoB;
  int **op_flag,*T_op_flag,*T_k_ID;
  double *T_KGrids2,*T_KGrids3;
  double k2,k3,tmp;
  double TStime,TEtime;

  int *MP;
  int *order_GA;
  int *My_NZeros;
  int *SP_NZeros;
  int *SP_Atoms;  

  int size_H1;
  int myworld1; 
  int numprocs1,myid1;
  int Num_Comm_World1;
  int *NPROCS_ID1;
  int *Comm_World1;
  int *NPROCS_WD1;
  int *Comm_World_StartID1;
  MPI_Comm *MPI_CommWD1;

  double **DM1,*TDM1;
  double **H1,*S1;

  MPI_Comm_size(comm1,&numprocs0);
  MPI_Comm_rank(comm1,&myid0);

  dtime(&TStime);

  if (myid0==Host_ID){
    printf("<TRAN_DFT>\n");
  }

  /***********************************
  ChemP_e0 will be used in outputfile1
  ************************************/

  ChemP_e0[0] = ChemP_e[0];
  ChemP_e0[1] = ChemP_e[1];

  /***********************************
            initialize CDM
  ************************************/

  for (spin=0; spin<=SpinP_switch; spin++) {
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++) {
      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	GB_AN = natn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];
	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){
	    CDM[spin][MA_AN][LB_AN][i][j] = 0.0;
	  }
	}
      }
    }
  }

  /***********************************
        set up operation flag
  ************************************/

  op_flag = (int**)malloc(sizeof(int*)*TRAN_Kspace_grid2); 
  for (i2=0; i2<TRAN_Kspace_grid2; i2++){
    op_flag[i2] = (int*)malloc(sizeof(int)*TRAN_Kspace_grid3); 
    for (i3=0; i3<TRAN_Kspace_grid3; i3++){
      op_flag[i2][i3] = -999;
    }
  }

  for (i2=0; i2<TRAN_Kspace_grid2; i2++){
    for (i3=0; i3<TRAN_Kspace_grid3; i3++){

      if (op_flag[i2][i3]<0){ 

	if ( (TRAN_Kspace_grid2-1-i2)==i2 && (TRAN_Kspace_grid3-1-i3)==i3 ){
	  op_flag[i2][i3] = 1;
	}
	else{
	  op_flag[i2][i3] = 2;
	  op_flag[TRAN_Kspace_grid2-1-i2][TRAN_Kspace_grid3-1-i3] = 0;
	}
      }

    }
  }

  /***********************************
       one-dimentionalize for MPI
  ************************************/

  T_knum = 0;
  for (i2=0; i2<TRAN_Kspace_grid2; i2++){
    for (i3=0; i3<TRAN_Kspace_grid3; i3++){
      if (0<op_flag[i2][i3]) T_knum++;  
    }
  }         

  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);
  T_op_flag = (int*)malloc(sizeof(int)*T_knum);
  T_k_ID = (int*)malloc(sizeof(int)*T_knum);

  T_knum = 0;

  for (i2=0; i2<TRAN_Kspace_grid2; i2++){

    k2 = -0.5 + (2.0*(double)i2+1.0)/(2.0*(double)TRAN_Kspace_grid2) + Shift_K_Point;

    for (i3=0; i3<TRAN_Kspace_grid3; i3++){

      k3 = -0.5 + (2.0*(double)i3+1.0)/(2.0*(double)TRAN_Kspace_grid3) - Shift_K_Point;

      if (0<op_flag[i2][i3]){  

        T_KGrids2[T_knum] = k2;
        T_KGrids3[T_knum] = k3;
        T_op_flag[T_knum] = op_flag[i2][i3];

        T_knum++;        
      }
    }
  }

  /***************************************************
   allocate calculations of k-points into processors 
  ***************************************************/

  if (numprocs0<T_knum){

    /* set parallel_mode */
    parallel_mode = 0;

    /* allocation of kloop to ID */     

    for (ID=0; ID<numprocs0; ID++){
      tmp = (double)T_knum/(double)numprocs0;
      S_knum = (int)((double)ID*(tmp+1.0e-12)); 
      E_knum = (int)((double)(ID+1)*(tmp+1.0e-12)) - 1;
      if (ID==(numprocs0-1)) E_knum = T_knum - 1;
      if (E_knum<0)          E_knum = 0;

      for (k=S_knum; k<=E_knum; k++){
        /* ID in the first level world */
        T_k_ID[k] = ID;
      }
    }

    /* find own informations */

    tmp = (double)T_knum/(double)numprocs0; 
    S_knum = (int)((double)myid0*(tmp+1.0e-12)); 
    E_knum = (int)((double)(myid0+1)*(tmp+1.0e-12)) - 1;
    if (myid0==(numprocs0-1)) E_knum = T_knum - 1;
    if (E_knum<0)             E_knum = 0;

    num_kloop0 = E_knum - S_knum + 1;

  }

  else {

    /* set parallel_mode */
    parallel_mode = 1;
    num_kloop0 = 1;

    Num_Comm_World1 = T_knum;

    NPROCS_ID1 = (int*)malloc(sizeof(int)*numprocs0);
    Comm_World1 = (int*)malloc(sizeof(int)*numprocs0);
    NPROCS_WD1 = (int*)malloc(sizeof(int)*Num_Comm_World1);
    Comm_World_StartID1 = (int*)malloc(sizeof(int)*Num_Comm_World1);
    MPI_CommWD1 = (MPI_Comm*)malloc(sizeof(MPI_Comm)*Num_Comm_World1);

    Make_Comm_Worlds(comm1, myid0, numprocs0, Num_Comm_World1, &myworld1, MPI_CommWD1, 
		     NPROCS_ID1, Comm_World1, NPROCS_WD1, Comm_World_StartID1);

    MPI_Comm_size(MPI_CommWD1[myworld1],&numprocs1);
    MPI_Comm_rank(MPI_CommWD1[myworld1],&myid1);

    S_knum = myworld1;

    /* allocate k-points into processors */
    
    for (k=0; k<T_knum; k++){
      /* ID in the first level world */
      T_k_ID[k] = Comm_World_StartID1[k];
    }

  }

  /*************************************************************
   one-dimensitonalize H and S and store them in a compact form  
  *************************************************************/

  MP = (int*)malloc(sizeof(int)*(atomnum+1));
  order_GA = (int*)malloc(sizeof(int)*(atomnum+1));

  My_NZeros = (int*)malloc(sizeof(int)*numprocs0);
  SP_NZeros = (int*)malloc(sizeof(int)*numprocs0);
  SP_Atoms = (int*)malloc(sizeof(int)*numprocs0);

  size_H1 = Get_OneD_HS_Col(0, nh[0], S1, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

  DM1 = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (k=0; k<(SpinP_switch+1); k++){
    DM1[k] = (double*)malloc(sizeof(double)*size_H1);
    for (i=0; i<size_H1; i++) DM1[k][i] = 0.0;
  }

  TDM1 = (double*)malloc(sizeof(double)*size_H1);

  H1 = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (k=0; k<(SpinP_switch+1); k++){
    H1[k] = (double*)malloc(sizeof(double)*size_H1);
  }

  S1 = (double*)malloc(sizeof(double)*size_H1);

  for (k=0; k<(SpinP_switch+1); k++){
    size_H1 = Get_OneD_HS_Col(1, nh[k], H1[k], MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  }

  size_H1 = Get_OneD_HS_Col(1, CntOLP, S1, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

  /***********************************************************
   start "kloop0"
  ***********************************************************/

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = S_knum + kloop0;

    k2 = T_KGrids2[kloop];
    k3 = T_KGrids3[kloop];
    k_op = T_op_flag[kloop];

    if (parallel_mode){

      TRAN_DFT_Kdependent(MPI_CommWD1[myworld1],
			  parallel_mode, numprocs1, myid1,
			  level_stdout, iter, SpinP_switch, k2, k3, k_op, order_GA,
                          DM1,H1,S1,
                          nh, ImNL, CntOLP,
			  atomnum, Matomnum, WhatSpecies, Spe_Total_CNO, FNAN,
			  natn, ncn, M2G, G2ID, atv_ijk, List_YOUSO, CDM, EDM, Eele0, Eele1);
    }
    else{

      TRAN_DFT_Kdependent(comm1,
			  parallel_mode, 1, 0,
			  level_stdout, iter, SpinP_switch, k2, k3, k_op, order_GA,
                          DM1,H1,S1,
                          nh, ImNL, CntOLP,
			  atomnum, Matomnum, WhatSpecies, Spe_Total_CNO, FNAN,
			  natn, ncn, M2G, G2ID, atv_ijk, List_YOUSO, CDM, EDM, Eele0, Eele1);
    }

  } /* kloop0 */

  /* MPI communication of DM1 */

  tmp = 1.0/(double)(TRAN_Kspace_grid2*TRAN_Kspace_grid3);
  
  for (k=0; k<=SpinP_switch; k++) {

    int itot0,Anum,Bnum;  
    double co,si,kRn;  

    MPI_Allreduce( DM1[k], TDM1, size_H1, MPI_DOUBLE, MPI_SUM, comm1);

    itot0 = 0; 

    for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {

      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];
      MA_AN = F_G2M[GA_AN]; 

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){

	GB_AN = natn[GA_AN][LB_AN];
	RnB = ncn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];
	Bnum = MP[GB_AN];
	l1 = atv_ijk[RnB][1];
	l2 = atv_ijk[RnB][2];
	l3 = atv_ijk[RnB][3];

	kRn = k2*(double)l2 + k3*(double)l3;
	si = (double)k_op*sin(2.0*PI*kRn);
	co = (double)k_op*cos(2.0*PI*kRn);

	for (i=0;i<tnoA;i++) {
	  for (j=0;j<tnoB;j++) {

	    if (1<=MA_AN && MA_AN<=Matomnum){   
              CDM[k][MA_AN][LB_AN][i][j] = TDM1[itot0]*tmp;
	    }

            itot0++;        

	  }
	}
      }
    }

  } /* k */

  /***********************************************************
           overwrite CMD with l1=-1 or l1=1 by zero
  ***********************************************************/

  for (spin=0; spin<=SpinP_switch; spin++) {
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++) {

      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){

	GB_AN = natn[GA_AN][LB_AN];
        RnB = ncn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];
        l1 = atv_ijk[RnB][1];
        l2 = atv_ijk[RnB][2];
        l3 = atv_ijk[RnB][3];

        /* DM between C-L or C-R */
        if (l1==-1 || l1==1){
	  for (i=0; i<tnoA; i++){
	    for (j=0; j<tnoB; j++){
	      CDM[spin][MA_AN][LB_AN][i][j] = 0.0;
	    }
	  }
        }
      }
    }
  }

  /***********************************************************
     calculate (partial) decomposed Mulliken population by 
               density matrix between CL or CR

     This contribution will be added in Mulliken_Charge.c
  ***********************************************************/

  {
    int MA_AN,GA_AN,wanA,tnoA,GA_AN_e,LB_AN_e,iside;
    int direction,Rn_e,GB_AN_e;
    double sum;

    /* initialize */

    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){

      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (spin=0; spin<=SpinP_switch; spin++) {
	for (i=0; i<tnoA; i++){
          TRAN_DecMulP[spin][MA_AN][i] = 0.0;
	}
      }
    }

    /* Left lead */

    iside = 0;
    direction = -1;

    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){

      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      if (TRAN_region[GA_AN]%10==2){

        GA_AN_e =  TRAN_Original_Id[GA_AN];

        for (spin=0; spin<=SpinP_switch; spin++) {
	  for (i=0; i<tnoA; i++){

	    sum = 0.0;

	    for (LB_AN_e=0; LB_AN_e<=FNAN_e[iside][GA_AN_e]; LB_AN_e++){

	      GB_AN_e = natn_e[iside][GA_AN_e][LB_AN_e];
	      Rn_e = ncn_e[iside][GA_AN_e][LB_AN_e];
	      wanB = WhatSpecies_e[iside][GB_AN_e];
	      tnoB = Spe_Total_CNO_e[iside][wanB];
	      l1 = atv_ijk_e[iside][Rn_e][1];

	      if (l1==direction) {
		for (j=0; j<tnoB; j++){
		  sum += OLP_e[iside][0][GA_AN_e][LB_AN_e][i][j]*
                         DM_e[iside][0][spin][GA_AN_e][LB_AN_e][i][j]; 
		}
	      }
	    }

	    TRAN_DecMulP[spin][MA_AN][i] = sum;

	  } /* i */
	} /* spin */
      }

    } /* MA_AN */

    /* Right lead */

    iside = 1;
    direction = 1;

    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){

      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      if (TRAN_region[GA_AN]%10==3){

        GA_AN_e = TRAN_Original_Id[GA_AN];

        for (spin=0; spin<=SpinP_switch; spin++) {
	  for (i=0; i<tnoA; i++){

	    sum = 0.0;

	    for (LB_AN_e=0; LB_AN_e<=FNAN_e[iside][GA_AN_e]; LB_AN_e++){

	      GB_AN_e = natn_e[iside][GA_AN_e][LB_AN_e];
	      Rn_e = ncn_e[iside][GA_AN_e][LB_AN_e];
	      wanB = WhatSpecies_e[iside][GB_AN_e];
	      tnoB = Spe_Total_CNO_e[iside][wanB];
	      l1 = atv_ijk_e[iside][Rn_e][1];

	      if (l1==direction) {
		for (j=0; j<tnoB; j++){
		  sum += OLP_e[iside][0][GA_AN_e][LB_AN_e][i][j]*
                         DM_e[iside][0][spin][GA_AN_e][LB_AN_e][i][j]; 
		}
	      }
	    }

	    TRAN_DecMulP[spin][MA_AN][i] = sum;

	  } /* i */
	} /* spin */
      } 
    } /* MA_AN */

  }

  /* set Eele as 0 */

  Eele0[0] = 0.0;
  Eele0[1] = 0.0;
  Eele1[0] = 0.0;
  Eele1[1] = 0.0;

  /* free arrays */

  for (i2=0; i2<TRAN_Kspace_grid2; i2++){
    free(op_flag[i2]);
  }
  free(op_flag);

  free(T_KGrids2);
  free(T_KGrids3);
  free(T_op_flag);
  free(T_k_ID);

  if (T_knum<=numprocs0){

    if (Num_Comm_World1<=numprocs0){
      MPI_Comm_free(&MPI_CommWD1[myworld1]);
    }

    free(NPROCS_ID1);
    free(Comm_World1);
    free(NPROCS_WD1);
    free(Comm_World_StartID1);
    free(MPI_CommWD1);
  }

  free(MP);
  free(order_GA);

  free(My_NZeros);
  free(SP_NZeros);
  free(SP_Atoms);

  for (k=0; k<(SpinP_switch+1); k++){
    free(DM1[k]);
  }
  free(DM1);

  free(TDM1);

  for (k=0; k<(SpinP_switch+1); k++){
    free(H1[k]);
  }
  free(H1);

  free(S1);

  /* for elapsed time */
  dtime(&TEtime);
  return TEtime - TStime;
}




/*
 *   calculate CDM from nh, CntOLP, ... 
 *
 *       an alternative routine for Cluster_DFT or Band_DFT 
 *
 */




static void TRAN_DFT_Kdependent(
			  /* input */
			  MPI_Comm comm1,
                          int parallel_mode,
                          int numprocs,
                          int myid,
			  int level_stdout,
			  int iter,
			  int SpinP_switch,
                          double k2,
                          double k3,
                          int k_op,
                          int *order_GA,
                          double **DM1,
                          double **H1,
                          double *S1,
			  double *****nh,  /* H */
			  double *****ImNL, /* not used, SO-coupling */
			  double ****CntOLP, 
			  int atomnum,
			  int Matomnum,
			  int *WhatSpecies,
			  int *Spe_Total_CNO,
			  int *FNAN,
			  int **natn, 
			  int **ncn,
			  int *M2G, 
			  int *G2ID, 
			  int **atv_ijk,
			  int *List_YOUSO,
			  /* output */
			  double *****CDM,  /* output, charge density */
			  double *****EDM,  /* not used */
			  double Eele0[2], double Eele1[2]) /* not used */

#define GC_ref(i,j) GC[ NUM_c*((j)-1) + (i)-1 ] 
#define Gless_ref(i,j) Gless[ NUM_c*((j)-1) + (i)-1 ]

{
  int i,j,k,iside; 
  int *MP;
  int  iw,iw_method;
  dcomplex w, w_weight;
  dcomplex *GC,*GRL,*GRR,*SigmaL, *SigmaR; 
  dcomplex *v1,*Gless,*GCLorR;
  dcomplex **v2;
  double dum;
  double TStime,TEtime;
  int MA_AN, GA_AN, wanA, tnoA, Anum;
  int LB_AN, GB_AN, wanB, tnoB, Bnum; 

  /* debug */
  double **Density;
  double density_sum;
  /* end debug*/
  
  static int ID;
  int **iwIdx, Miwmax, Miw,iw0; 
  double time_a0, time_a1, time_a2; 
  
  /* parallel setup */

  iwIdx=(int**)malloc(sizeof(int*)*numprocs);
  Miwmax = (tran_omega_n_scf)/numprocs+1;
  for (i=0; i<numprocs; i++) {
    iwIdx[i]=(int*)malloc(sizeof(int)*Miwmax);
  }

  TRAN_Distribute_Node_Idx(0, tran_omega_n_scf-1, numprocs, Miwmax,
                           iwIdx); /* output */

  /* setup MP */
  TRAN_Set_MP(0, atomnum, WhatSpecies, Spe_Total_CNO, &NUM_c, MP);
  MP = (int*)malloc(sizeof(int)*(NUM_c+1));
  TRAN_Set_MP(1, atomnum, WhatSpecies, Spe_Total_CNO, &NUM_c, MP);
  
  /*debug */

  if (WRITE_DENSITY){
    Density = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
    for (k=0;k<=SpinP_switch;k++) {
      Density[k] = (double*)malloc(sizeof(double)*NUM_c);
      for (i=0;i<NUM_c;i++) { Density[k][i]=0.0; }
    }
  }
  /*end debug */
  
  /* initialize */
  TRAN_Set_Value_double(SCC,NUM_c*NUM_c,    0.0,0.0);
  TRAN_Set_Value_double(SCL,NUM_c*NUM_e[0], 0.0,0.0);
  TRAN_Set_Value_double(SCR,NUM_c*NUM_e[1], 0.0,0.0);
  for (k=0; k<=SpinP_switch; k++) {
    TRAN_Set_Value_double(HCC[k],NUM_c*NUM_c,    0.0,0.0);
    TRAN_Set_Value_double(HCL[k],NUM_c*NUM_e[0], 0.0,0.0);
    TRAN_Set_Value_double(HCR[k],NUM_c*NUM_e[1], 0.0,0.0);
  }

  /* set Hamiltonian and overlap matrices of left and right leads */

  TRAN_Set_SurfOverlap(comm1,"left", k2, k3);
  TRAN_Set_SurfOverlap(comm1,"right",k2, k3);

  /* set CC, CL and CR */

  TRAN_Set_CentOverlap(   comm1,
                          3,
                          SpinP_switch, 
                          k2,
                          k3,
                          order_GA,
                          H1,
                          S1,
                          nh,      /* input */
                          CntOLP,  /* input */
                          atomnum,
			  Matomnum,
			  M2G,
			  G2ID,
                          WhatSpecies,
                          Spe_Total_CNO,
                          FNAN,
                          natn,
                          ncn,
                          atv_ijk);


  if (MEASURE_TIME){
    dtime(&time_a0);
  }
  
  /* allocate */

  v2 = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1));

  for (k=0; k<=SpinP_switch; k++) {
    v2[k] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_c);
    TRAN_Set_Value_double( v2[k], NUM_c*NUM_c, 0.0, 0.0);
  }
  
  GC = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  GRL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[0]* NUM_e[0]);
  GRR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[1]* NUM_e[1]);
  SigmaL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  v1 = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  Gless = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c); 
  GCLorR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c); 
  
  if (2<=level_stdout){
    printf("NUM_c=%d, NUM_e= %d %d\n",NUM_c, NUM_e[0], NUM_e[1]);
    printf("# of freq. to calculate G =%d\n",tran_omega_n_scf);
  }

  /*parallel global iw 0:tran_omega_n_scf-1 */
  /*parallel local  Miw 0:Miwmax-1 */
  /*parllel variable iw=iwIdx[myid][Miw] */

  for (Miw=0; Miw<Miwmax; Miw++) {

    iw = iwIdx[myid][Miw];

    if (iw>=0) {
      w = tran_omega_scf[iw];
      w_weight  = tran_omega_weight_scf[iw];
      iw_method = tran_integ_method_scf[iw]; 
    }
    
    /*
    printf("Miwmax=%3d Miw=%3d iw=%d of %d w=%le %le weight=%le %le method=%d\n",
            Miwmax,Miw,iw,tran_omega_n_scf, w.r,w.i, w_weight.r, w_weight.i, iw_method);
    */

    for (k=0; k<=SpinP_switch; k++) {

      if (iw>=0) {

        iside=0;

        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max,
                                   tran_surfgreen_eps, GRL);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL);

        iside=1;

        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max,
                                   tran_surfgreen_eps, GRR);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR);

        TRAN_Calc_CentGreen(w, NUM_c, SigmaL,SigmaR, HCC[k], SCC, GC);

        /***********************************************
                             G_{C} 
        ***********************************************/

	if (iw_method==1) {

	  /* start debug */

          if (WRITE_DENSITY){
  	    for (i=0;i<NUM_c;i++) {
	      /* imag( GC * weight ) */
	      Density[k][i] += GC_ref(i+1,i+1).r*w_weight.i + GC_ref(i+1,i+1).i*w_weight.r;
	    }
	  }

	  /* end debug */

	}  /* iw_method */

        /***********************************************
              based on the lesser Green's function
        ***********************************************/

        else if (iw_method==2)  {

          if (2<=level_stdout){
            printf("G_Lesser, iw=%d (of %d) w=%le %le\n",iw,tran_omega_n_scf,w.r, w.i);
	  }

          TRAN_Calc_CentGreenLesser(w, ChemP_e, NUM_c, SigmaL,SigmaR,GC,v1, Gless);

#ifdef DEBUG
	  TRAN_Print2_dcomplex("Gless",NUM_c,NUM_c,Gless);
	  printf("exit after TRAN_Print2_dcomplex Gless\n");
	  exit(0);
#endif

          if (WRITE_DENSITY){
	    /* real( Gless * w_weight )  */
	    for (i=0;i<NUM_c;i++) {
	      Density[k][i] += Gless_ref(i+1,i+1).r*w_weight.r - Gless_ref(i+1,i+1).i*w_weight.i;
	    }

	    /*        printf("iw=%d w=%lf %lf weight=%lf %lf GC=%lf %lf, val=%lf\n",
	     *          iw,w.r , w.i, w_weight.r, w_weight.i, Density[k][0]);
	     */
	    /*end debug */
	  }
	}

        /***********************************************
              G_{C_L} in the non-equilibrium case
                  based on the GaussHG method
        ***********************************************/

        else if (iw_method==3){
	  TRAN_Calc_GC_LorR(iw_method, w, ChemP_e, NUM_c, NUM_e, SigmaL, GC, HCC[k], SCC, v1, GCLorR);
        }

        /***********************************************
              G_{C_R} in the non-equilibrium case
                  based on the GaussHG method
        ***********************************************/

        else if (iw_method==4){
	  TRAN_Calc_GC_LorR(iw_method, w, ChemP_e, NUM_c, NUM_e, SigmaR, GC, HCC[k], SCC, v1, GCLorR);
        }

	else {
	  printf("error, iw_method=%d",iw_method);
	  exit(10);
	}

        /***********************************************
            add it to construct the density matrix
        ***********************************************/

        if      (iw_method==1) {
          TRAN_Add_MAT( 1, NUM_c, w_weight, GC,     v2[k]);
        }
        else if (iw_method==2) {
          TRAN_Add_MAT( 2, NUM_c, w_weight, Gless,  v2[k]);
        }
        else if (iw_method==3) {
          TRAN_Add_MAT( 1, NUM_c, w_weight, GCLorR, v2[k]);
        }
        else if (iw_method==4) {
          TRAN_Add_MAT( 1, NUM_c, w_weight, GCLorR, v2[k]);
        }

      }  /* iw>=0 */
    } /* for k */
  } /* iw */

  if (MEASURE_TIME){
    dtime(&time_a1);
  }

  free(GCLorR);
  free(Gless);
  free(v1);
  free(SigmaR);
  free(SigmaL);
  free(GRR);
  free(GRL);
  free(GC);

  /***********************************************
          calculation of density matrix
  ***********************************************/

  {
    int l1,l2,l3,RnB;
    int size_v3,itot,itot0;
    double kRn,si,co,re,im;
    double *my_v3;
    double *v3;

    /* find the size of v3 */

    size_v3 = 0;

    for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {

      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){

	GB_AN = natn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];
	Bnum = MP[GB_AN];

	for (i=0;i<tnoA;i++) {
	  for (j=0;j<tnoB;j++) {
	    size_v3++;
	    size_v3++;
	  }
	}
      }
    }  

    /* allocate arrays */
  
    my_v3 = (double*)malloc(sizeof(double)*size_v3);
    v3 = (double*)malloc(sizeof(double)*size_v3);

    /* set up v3 */

#define v_idx(i,j)   ( ((j)-1)*NUM_c + (i)-1 ) 

    for (k=0; k<=SpinP_switch; k++) {

      itot = 0;

      for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {

	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

        for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){

	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  Bnum = MP[GB_AN];

	  for (i=0;i<tnoA;i++) {
	    for (j=0;j<tnoB;j++) {
	      my_v3[itot++] = v2[k][ v_idx( Anum+i, Bnum+j) ].r;
	      my_v3[itot++] = v2[k][ v_idx( Anum+i, Bnum+j) ].i;
	    }
	  }
	}
      }  

      if (parallel_mode){
        MPI_Allreduce( my_v3, v3, itot, MPI_DOUBLE, MPI_SUM, comm1);
      }
      else {
        for (i=0; i<itot; i++) {
          v3[i] = my_v3[i]; 
	}
      }

      /* v3 -> CDM */

      itot = 0;
      itot0 = 0; 

      for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {

	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){

	  GB_AN = natn[GA_AN][LB_AN];
	  RnB = ncn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  Bnum = MP[GB_AN];
	  l1 = atv_ijk[RnB][1];
	  l2 = atv_ijk[RnB][2];
	  l3 = atv_ijk[RnB][3];

	  kRn = k2*(double)l2 + k3*(double)l3;
	  si = (double)k_op*sin(2.0*PI*kRn);
	  co = (double)k_op*cos(2.0*PI*kRn);

	  for (i=0;i<tnoA;i++) {
	    for (j=0;j<tnoB;j++) {
	      re = v3[itot++];
	      im = v3[itot++];
	      DM1[k][itot0++] += (re*co + im*si)/(double)numprocs; 
	    }
	  }
	}
      }

    } /* k */

    /* free arrays */

    free(my_v3);
    free(v3);
  }

  if (MEASURE_TIME){
    MPI_Barrier(comm1);
    dtime(&time_a2);
    printf("TRAN_DFT(%d)> calculaiton (%le)\n",myid, time_a1-time_a0 ); 
  }


  /* free arrays */

  for (k=0; k<=SpinP_switch; k++) {
    free(v2[k]);
  }
  free(v2);

  if (WRITE_DENSITY){
    for (k=SpinP_switch; k>=0; k--) {
      free(Density[k]);
    }
    free(Density);
  }

  free(MP);

  for (i=0;i<numprocs;i++) {
    free(iwIdx[i]);
  }
  free(iwIdx);

}



static void TRAN_Add_MAT(
    int mode,
    int NUM_c,
    dcomplex w_weight,
    dcomplex *v,
    dcomplex *out
)
{
  int i;
  double dum; 

  switch (mode) {
  case 0:
    for (i=0;i<NUM_c*NUM_c;i++) {
      out[i].r =  0.0;
      out[i].i =  0.0;
    }
    break;
  case 1:
    for (i=0;i<NUM_c*NUM_c;i++) {
      out[i].r += ( v[i].r*w_weight.r - v[i].i*w_weight.i );
      out[i].i += ( v[i].r*w_weight.i + v[i].i*w_weight.r );
    }
    break;
  case 2:
    for (i=0;i<NUM_c*NUM_c;i++) {
      /* real(G^less*weight) */
      dum = v[i].r*w_weight.r - v[i].i*w_weight.i;
      out[i].r += dum;
    }
    break;

  }
}

