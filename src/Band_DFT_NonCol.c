/**********************************************************************
  Band_DFT_NonCol.c:

     Band_DFT_NonCol.c is a subroutine to perform band calculations
     based on a non-collinear DFT

  Log of Band_DFT_NonCol.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#define  measure_time   0

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif



static double Band_non_collinear(
             int knum_i, int knum_j, int knum_k,
             int SpinP_switch,
             double *****nh,
             double *****ImNL,
             double ****CntOLP,
             double *****CDM,
             double *****EDM,
             double Eele0[2], double Eele1[2]);

static double Band_non_collinear2(int SCF_iter,
				  double *koS,
				  dcomplex **S,
				  int knum_i, int knum_j, int knum_k,
				  int SpinP_switch,
				  double *****nh,
				  double *****ImNL,
				  double ****CntOLP,
				  double *****CDM,
				  double *****EDM,
				  double Eele0[2], double Eele1[2]);


double Band_DFT_NonCol(int SCF_iter,
                       double *koS,
                       dcomplex **S,
                       int knum_i, int knum_j, int knum_k,
		       int SpinP_switch,
		       double *****nh,
		       double *****ImNL,
		       double ****CntOLP,
		       double *****CDM,
		       double *****EDM,
		       double Eele0[2], double Eele1[2])
{
  double time0;

  time0 = Band_non_collinear2(SCF_iter,koS,S,knum_i, knum_j, knum_k, SpinP_switch,
                               nh, ImNL, CntOLP, CDM, EDM, Eele0, Eele1 );

  /*
  time0 = Band_non_collinear( knum_i, knum_j, knum_k, SpinP_switch,
                              nh, ImNL, CntOLP, CDM, EDM, Eele0, Eele1 );
  */

  return time0;
}





static double Band_non_collinear2(int SCF_iter,
				  double *koS,
				  dcomplex **S,
				  int knum_i, int knum_j, int knum_k,
				  int SpinP_switch,
				  double *****nh,
				  double *****ImNL,
				  double ****CntOLP,
				  double *****CDM,
				  double *****EDM,
				  double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  int i,j,k,l,n,n2,wan;
  int ii1,m,mn,lmn,lmax;
  int i1,j1,jj1,jj2,po,spin,n1;
  int ia,ib,ja,jb;
  int num2,RnB,l1,l2,l3,loop_num;
  int ct_AN,h_AN,wanA,tnoA,wanB,tnoB;
  int MA_AN,GA_AN,Anum;
  double time0,sum_weights;
  double Imdum2,Imsum,ImsumE;
  int LB_AN,GB_AN,Bnum;
  double EV_cut0;
  double snum_i,snum_j,snum_k,k1,k2,k3,sum,sumi,Num_State,FermiF;
  double sum_r0,sum_i0;
  double sum_r1,sum_i1;
  double sum_r00,sum_i00;
  double sum_r01,sum_i01;
  double sum_r10,sum_i10;
  double sum_r11,sum_i11;
  double av_num;
  double x,Dnum,Dnum2,AcP,ChemP_MAX,ChemP_MIN;
  double **EIGEN;
  double d1,d2,d3,d4,d5,d6,d7,d8,d9,d10;
  double *VecFkw;
  double *VecFkwE;
  dcomplex **H,**C;
  dcomplex Ctmp1,Ctmp2;
  int ***k_op,*T_k_op,*T_k_ID;
  int ii,ij,ik;
  int S_knum,E_knum,T_knum,e1,s1,kloop,kloop0,num_kloop0;
  double *KGrids1,*KGrids2,*KGrids3;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  double u2,v2,uv,vu,tmp;
  double My_Eele1[2]; 
  double kw,eig,Fkw;
  double TZ,dum,sumE,kRn,si,co;
  double Resum,ResumE,Redum,Redum2,Imdum;

  double ResumA,ResumEA,RedumA,Redum2A,ImdumA;
  double ResumB,ResumEB,RedumB,Redum2B,ImdumB;
  double ResumC,ResumEC,RedumC,Redum2C,ImdumC;
  double Imdum2D,ImsumED;
  double u2A,v2A,uvA,vuA;
  double u2B,v2B,uvB,vuB;
  double u2C,v2C,uvC,vuC;

  double TStime,TEtime,SiloopTime,EiloopTime;
  double Stime,Etime;
  double FermiEps=1.0e-14;
  double x_cut=30.0;
  char file_EV[YOUSO10];
  FILE *fp_EV;
  char buf[fp_bsize];          /* setvbuf */

  int *MP;
  int *order_GA;
  double *ko;
  double *S1;
  double *RH0;
  double *RH1;
  double *RH2;
  double *RH3;
  double *IH0;
  double *IH1;
  double *IH2;

  double *CDM0;
  double *CDM1;
  double *CDM2;
  double *CDM3;

  double *EDM0;
  double *EDM1;
  double *EDM2;
  double *EDM3;

  double *iDM0;
  double *iDM1;

  int *is1,*ie1;
  int *is12,*ie12;
  int *is2,*ie2;
  MPI_Status *stat_send;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  int *My_NZeros;
  int *SP_NZeros;
  int *SP_Atoms;
  MPI_Comm *MPI_CommWD_CDM1; 
  int *MPI_CDM1_flag;  

  double lumos;
  int Rn,AN;
  int size_H1;
  int all_knum;
  int MaxN,ks,i0;
  int ID,ID0,ID1;
  int numprocs0,myid0;
  int numprocs1,myid1;
  int myworld1;
  int parallel_mode;
  int Num_Comm_World1;
  int *NPROCS_ID1;
  int *NPROCS_WD1;
  int *Comm_World1;
  int *Comm_World_StartID1;
  MPI_Comm *MPI_CommWD1;

  double time1,time2,time3;
  double time4,time5,time6;
  double time7,time8,time9;
  double time10,time11,time12;
  double time13,time14,time15;
  double time16,time17,time18;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs0);
  MPI_Comm_rank(mpi_comm_level1,&myid0);

  dtime(&TStime);

  time1 = 0.0;
  time2 = 0.0;
  time3 = 0.0;
  time4 = 0.0;
  time5 = 0.0;
  time6 = 0.0;
  time7 = 0.0;
  time8 = 0.0;
  time9 = 0.0;
  time10 = 0.0;
  time11 = 0.0;
  time12 = 0.0;
  time13 = 0.0;
  time14 = 0.0;
  time15 = 0.0;
  time16 = 0.0;
  time17 = 0.0;
  time18 = 0.0;

  /****************************************************
             calculation of the array size
  ****************************************************/

  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n += Spe_Total_CNO[wanA];
  }
  n2 = 2*n + 2;

  /****************************************************
                        find TZ
  ****************************************************/

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ += Spe_Core_Charge[wan];
  }

  /***********************************************
     find the number of states to be solved 
  ***********************************************/

  lumos = (double)2*n*0.050;
  if (lumos<60.0) lumos = 60.0;
  MaxN = (TZ-system_charge) + (int)lumos;
  if ( (2*n)<MaxN ) MaxN = 2*n;

  /****************************************************
   Allocation

   int       MP[List_YOUSO[1]]
   double    ko[n2]
   double    koS[n+1];
   dcomplex  H[n2][n2]
   dcomplex  S[n2][n2]
   double    M1[n2]
   dcomplex  C[n2][n2]
   double    EIGEN[List_YOUSO[27]]
                  [List_YOUSO[28]]
                  [List_YOUSO[29]]
                  [n2] 
   double    KGrids1[List_YOUSO[27]]
   double    KGrids2[List_YOUSO[28]]
   double    KGrids3[List_YOUSO[29]]
   int       k_op[knum_i][knum_j][knum_k]
  ****************************************************/

  VecFkw = (double*)malloc(sizeof(double)*n2);
  VecFkwE = (double*)malloc(sizeof(double)*n2);

  H = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    H[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  C = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    C[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  KGrids1 = (double*)malloc(sizeof(double)*List_YOUSO[27]);
  KGrids2 = (double*)malloc(sizeof(double)*List_YOUSO[28]);
  KGrids3 = (double*)malloc(sizeof(double)*List_YOUSO[29]);

  k_op=(int***)malloc(sizeof(int**)*(knum_i));
  for (i=0;i<knum_i; i++) {
    k_op[i]=(int**)malloc(sizeof(int*)*(knum_j));
    for (j=0;j<knum_j; j++) {
      k_op[i][j]=(int*)malloc(sizeof(int)*(knum_k));
    }
  }
  /*  k_op[0:knum_i-1][0:knum_j-1][0:knum_k-1] */

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  order_GA = (int*)malloc(sizeof(int)*(List_YOUSO[1]+1));

  My_NZeros = (int*)malloc(sizeof(int)*numprocs0);
  SP_NZeros = (int*)malloc(sizeof(int)*numprocs0);
  SP_Atoms = (int*)malloc(sizeof(int)*numprocs0);

  MPI_CommWD_CDM1 = (MPI_Comm*)malloc(sizeof(MPI_Comm)*numprocs0);
  MPI_CDM1_flag = (int*)malloc(sizeof(int)*numprocs0);
  for (i=0; i<numprocs0; i++) MPI_CDM1_flag[i] = 0;

  ko = (double*)malloc(sizeof(double)*n2);

  /****************************************************
    PrintMemory
  ****************************************************/

  if (firsttime){

    PrintMemory("Band_DFT: ko",sizeof(double)*n2,NULL);
    PrintMemory("Band_DFT: H",sizeof(dcomplex)*n2*n2,NULL);
    PrintMemory("Band_DFT: S",sizeof(dcomplex)*n2*n2,NULL);
    PrintMemory("Band_DFT: C",sizeof(dcomplex)*n2*n2,NULL);
    PrintMemory("Band_DFT: M1",sizeof(double)*n2,NULL);

    PrintMemory("Band_DFT: EIGEN",
          sizeof(double)*List_YOUSO[27]*
                         List_YOUSO[28]*
                         List_YOUSO[29]*n2,NULL);
    PrintMemory("Band_DFT: ko_op",
          sizeof(int)*knum_i*knum_j*knum_k,NULL);
  }

  /* for PrintMemory */
  firsttime=0;

  /***********************************************
              k-points by regular mesh 
  ***********************************************/

  if (way_of_kpoint==1){

    snum_i = (double)knum_i;
    snum_j = (double)knum_j;
    snum_k = (double)knum_k;

    for (i=0; i<=(knum_i-1); i++){
      if (knum_i==1){
	k1 = 0.0;
      }
      else {
	k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_i) + Shift_K_Point;
      }
      KGrids1[i]=k1;
    }
    for (i=0; i<=(knum_j-1); i++){
      if (knum_j==1){
	k1 = 0.0;
      }
      else {
	k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_j) - Shift_K_Point;
      }
      KGrids2[i]=k1;
    }
    for (i=0; i<=(knum_k-1); i++){
      if (knum_k==1){
	k1 = 0.0;
      }
      else {
	k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_k) + 2.0*Shift_K_Point;
      }
      KGrids3[i]=k1;
    }

    if (myid0==Host_ID){
      printf(" KGrids1: ");
      for (i=0;i<=knum_i-1;i++) printf("%9.5f ",KGrids1[i]);
      printf("\n");
      printf(" KGrids2: ");
      for (i=0;i<=knum_j-1;i++) printf("%9.5f ",KGrids2[i]);
      printf("\n");
      printf(" KGrids3: ");
      for (i=0;i<=knum_k-1;i++) printf("%9.5f ",KGrids3[i]);
      printf("\n");
    }
    /**************************************************************
   k_op[i][j][k]: weight of DOS 
                 =0   no calc.
                 =1   G-point
             note that there is no simple inversion symmetry for
               the wave functions in noncollinear case 
                 =1  other  
    *************************************************************/

    for (i=0;i<=knum_i-1;i++) {
      for (j=0;j<=knum_j-1;j++) {
	for (k=0;k<=knum_k-1;k++) {
	  k_op[i][j][k]=-999;
	}
      }
    }
    for (i=0;i<=knum_i-1;i++) {
      for (j=0;j<=knum_j-1;j++) {
	for (k=0;k<=knum_k-1;k++) {
	  if ( k_op[i][j][k]==-999 ) {
	    k_inversion(i,j,k,knum_i,knum_j,knum_k,&ii,&ij,&ik);
	    if ( i==ii && j==ij && k==ik ) {
	      /*
		printf("* : %d %d %d (%f %f %f)\n",i,j,k,
		KGrids1[i], KGrids2[j], KGrids3[k]);
	      */
	      k_op[i][j][k]    = 1;
	    }

	    else {

	      /* note that there is no simple inversion symmetry for
		 the wave functions in noncollinear case */

	      k_op[i][j][k]    = 1;
	      k_op[ii][ij][ik] = 1;
	    }
	  }
	} /* k */
      } /* j */
    } /* i */

    /***********************************
       one-dimentionalize for MPI
    ************************************/

    T_knum = 0;
    for (i=0; i<=(knum_i-1); i++){
      for (j=0; j<=(knum_j-1); j++){
	for (k=0; k<=(knum_k-1); k++){
	  if (0<k_op[i][j][k]){
	    T_knum++;
	  }
	}
      }
    }

    /* allocation of arrays */

    T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);
    T_k_op    = (int*)malloc(sizeof(int)*T_knum);
    T_k_ID    = (int*)malloc(sizeof(int)*T_knum);

    EIGEN = (double**)malloc(sizeof(double*)*T_knum);
    for (j=0; j<T_knum; j++){
      EIGEN[j] = (double*)malloc(sizeof(double)*n2);
      for (i=0; i<n2; i++)  EIGEN[j][i] = 1.0e+5;
    }

    /* set T_KGrid1,2,3 and T_k_op */

    T_knum = 0;
    for (i=0; i<=(knum_i-1); i++){
      for (j=0; j<=(knum_j-1); j++){
	for (k=0; k<=(knum_k-1); k++){
	  if (0<k_op[i][j][k]){

	    T_KGrids1[T_knum] = KGrids1[i];
	    T_KGrids2[T_knum] = KGrids2[j];
	    T_KGrids3[T_knum] = KGrids3[k];

	    T_k_op[T_knum]    = k_op[i][j][k];

	    T_knum++;
	  }
	}
      }
    }
  }

  /***********************************************
                 Monkhorst-Pack k-points 
  ***********************************************/

  else if (way_of_kpoint==2){

    T_knum = num_non_eq_kpt; 

    T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);
    T_k_op    = (int*)malloc(sizeof(int)*T_knum);
    T_k_ID    = (int*)malloc(sizeof(int)*T_knum);

    EIGEN = (double**)malloc(sizeof(double*)*T_knum);
    for (j=0; j<T_knum; j++){
      EIGEN[j] = (double*)malloc(sizeof(double)*n2);
      for (i=0; i<n2; i++)  EIGEN[j][i] = 1.0e+5;
    }
   
    for (k=0; k<T_knum; k++){
      T_KGrids1[k] = NE_KGrids1[k];
      T_KGrids2[k] = NE_KGrids2[k];
      T_KGrids3[k] = NE_KGrids3[k];
      T_k_op[k]    = NE_T_k_op[k];
    }
  }

  /***********************************************
            calculate the sum of weights
  ***********************************************/

  sum_weights = 0.0;
  for (k=0; k<T_knum; k++){
    sum_weights += (double)T_k_op[k];
  }

  /***********************************************
         allocate k-points into proccessors 
  ***********************************************/

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
        /* ID in the zeroth level world */
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

    Make_Comm_Worlds(mpi_comm_level1, myid0, numprocs0, Num_Comm_World1, &myworld1, MPI_CommWD1, 
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

  /****************************************************
   find all_knum
   if (all_knum==1), all the calculation will be made 
   by the first diagonalization loop, and the second 
   diagonalization will be skipped. 
  ****************************************************/

  MPI_Allreduce(&num_kloop0, &all_knum, 1, MPI_INT, MPI_PROD, mpi_comm_level1);

  if (all_knum==1){

    /* allocation */ 

    stat_send = malloc(sizeof(MPI_Status)*numprocs1);
    request_send = malloc(sizeof(MPI_Request)*numprocs1);
    request_recv = malloc(sizeof(MPI_Request)*numprocs1);

    is1 = (int*)malloc(sizeof(int)*numprocs1);
    ie1 = (int*)malloc(sizeof(int)*numprocs1);

    is12 = (int*)malloc(sizeof(int)*numprocs1);
    ie12 = (int*)malloc(sizeof(int)*numprocs1);

    is2 = (int*)malloc(sizeof(int)*numprocs1);
    ie2 = (int*)malloc(sizeof(int)*numprocs1);

    /* make is1 and ie1, is12 and ie12 */ 

    if ( numprocs1<=n ){

      av_num = (double)n/(double)numprocs1;

      for (ID=0; ID<numprocs1; ID++){
	is1[ID] = (int)(av_num*(double)ID) + 1; 
	ie1[ID] = (int)(av_num*(double)(ID+1)); 
      }

      is1[0] = 1;
      ie1[numprocs1-1] = n; 

    }

    else{

      for (ID=0; ID<n; ID++){
	is1[ID] = ID + 1; 
	ie1[ID] = ID + 1;
      }
      for (ID=n; ID<numprocs1; ID++){
	is1[ID] =  1;
	ie1[ID] = -2;
      }
    }

    for (ID=0; ID<numprocs1; ID++){
      is12[ID] = 2*is1[ID] - 1;
      ie12[ID] = 2*ie1[ID];
    }

    /* make is2 and ie2 */ 

    if ( numprocs1<=MaxN ){

      av_num = (double)MaxN/(double)numprocs1;

      for (ID=0; ID<numprocs1; ID++){
	is2[ID] = (int)(av_num*(double)ID) + 1; 
	ie2[ID] = (int)(av_num*(double)(ID+1)); 
      }

      is2[0] = 1;
      ie2[numprocs1-1] = MaxN; 
    }

    else{
      for (ID=0; ID<MaxN; ID++){
	is2[ID] = ID + 1; 
	ie2[ID] = ID + 1;
      }
      for (ID=MaxN; ID<numprocs1; ID++){
	is2[ID] =  1;
	ie2[ID] = -2;
      }
    }

    /**********************************************
       for MPI communication of CDM1 and EDM1        
    **********************************************/

    {
    MPI_Group new_group,old_group; 
    int *new_ranks; 

    new_ranks = (int*)malloc(sizeof(int)*T_knum);

    /* ID: zeroth world */

    for (ID=0; ID<numprocs0; ID++){

      for (i=0; i<T_knum; i++){
        if (Comm_World_StartID1[i]<=ID){
          ks = i;
        }
      }

      new_ranks[0] = ID; 
      if (myid0==ID) MPI_CDM1_flag[ID] = 1;

      for (i=(ks+1); i<(T_knum+ks); i++){
        i0 = i % T_knum;
        /* id in the zeroth world */
        ID0 = Comm_World_StartID1[i0] + (ID % NPROCS_WD1[i0]);
        new_ranks[i-ks] = ID0;
        if (myid0==ID0) MPI_CDM1_flag[ID] = 1;
      }

      MPI_Comm_group(mpi_comm_level1, &old_group);

      /* define a new group */
      MPI_Group_incl(old_group,T_knum,new_ranks,&new_group);
      MPI_Comm_create(mpi_comm_level1, new_group, &MPI_CommWD_CDM1[ID]);
      MPI_Group_free(&new_group);
    }

    free(new_ranks); /* never forget cleaning! */
    }

  } /* if (all_knum==1) */

  /****************************************************
     store in each processor all the matrix elements
        for overlap and Hamiltonian matrices
  ****************************************************/

  dtime(&Stime);

  /* get size_H1 */

  size_H1 = Get_OneD_HS_Col(0, CntOLP,   &tmp, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

  S1 = (double*)malloc(sizeof(double)*size_H1);
  RH0 = (double*)malloc(sizeof(double)*size_H1);
  RH1 = (double*)malloc(sizeof(double)*size_H1);
  RH2 = (double*)malloc(sizeof(double)*size_H1);
  RH3 = (double*)malloc(sizeof(double)*size_H1);
  IH0 = (double*)malloc(sizeof(double)*size_H1);
  IH1 = (double*)malloc(sizeof(double)*size_H1);
  IH2 = (double*)malloc(sizeof(double)*size_H1);

  /* set RH0, RH1, RH2, RH3, IH0, IH1, IH2 */

  size_H1 = Get_OneD_HS_Col(1, CntOLP, S1,  MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  size_H1 = Get_OneD_HS_Col(1, nh[0],  RH0, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  size_H1 = Get_OneD_HS_Col(1, nh[1],  RH1, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  size_H1 = Get_OneD_HS_Col(1, nh[2],  RH2, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  size_H1 = Get_OneD_HS_Col(1, nh[3],  RH3, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
       && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){  
    
    /* nothing is done. */
  }
  else {
    size_H1 = Get_OneD_HS_Col(1, ImNL[0], IH0, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
    size_H1 = Get_OneD_HS_Col(1, ImNL[1], IH1, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
    size_H1 = Get_OneD_HS_Col(1, ImNL[2], IH2, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  }
 
  /* allocate
     CDM0, CDM1, CDM2, CDM3
     EDM0, EDM1, EDM2, EDM3
     iDM0, iDM1
  */

  CDM0 = (double*)malloc(sizeof(double)*size_H1);
  CDM1 = (double*)malloc(sizeof(double)*size_H1);
  CDM2 = (double*)malloc(sizeof(double)*size_H1);
  CDM3 = (double*)malloc(sizeof(double)*size_H1);

  EDM0 = (double*)malloc(sizeof(double)*size_H1);
  EDM1 = (double*)malloc(sizeof(double)*size_H1);
  EDM2 = (double*)malloc(sizeof(double)*size_H1);
  EDM3 = (double*)malloc(sizeof(double)*size_H1);

  iDM0 = (double*)malloc(sizeof(double)*size_H1);
  iDM1 = (double*)malloc(sizeof(double)*size_H1);

  dtime(&Etime);
  time1 += Etime - Stime; 

  /****************************************************
                      start kloop
  ****************************************************/

  dtime(&SiloopTime);

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    dtime(&Stime);

    kloop = S_knum + kloop0;

    k1 = T_KGrids1[kloop];
    k2 = T_KGrids2[kloop];
    k3 = T_KGrids3[kloop];

    /* make S and H */

    if (SCF_iter==1 || all_knum!=1){

      for (i=1; i<=n; i++){
        for (j=1; j<=n; j++){
	  S[i][j] = Complex(0.0,0.0);
        } 
      } 
    }

    for (i=1; i<=2*n; i++){
      for (j=1; j<=2*n; j++){
	H[i][j] = Complex(0.0,0.0);
      } 
    } 

    /* non-spin-orbit coupling and non-LDA+U */
    if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
	&& Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){  

      k = 0;
      for (AN=1; AN<=atomnum; AN++){
	GA_AN = order_GA[AN];
	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  Rn = ncn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  Bnum = MP[GB_AN];

	  l1 = atv_ijk[Rn][1];
	  l2 = atv_ijk[Rn][2];
	  l3 = atv_ijk[Rn][3];
	  kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

	  si = sin(2.0*PI*kRn);
	  co = cos(2.0*PI*kRn);

	  for (i=0; i<tnoA; i++){
	    for (j=0; j<tnoB; j++){

	      H[Anum+i  ][Bnum+j  ].r += co*RH0[k];
	      H[Anum+i  ][Bnum+j  ].i += si*RH0[k];

	      H[Anum+i+n][Bnum+j+n].r += co*RH1[k];
	      H[Anum+i+n][Bnum+j+n].i += si*RH1[k];
            
	      H[Anum+i  ][Bnum+j+n].r += co*RH2[k] - si*RH3[k];
	      H[Anum+i  ][Bnum+j+n].i += si*RH2[k] + co*RH3[k];

	      k++;

	    }

            if (SCF_iter==1 || all_knum!=1){
              
              k -= tnoB; 

	      for (j=0; j<tnoB; j++){

	        S[Anum+i  ][Bnum+j  ].r += co*S1[k];
	        S[Anum+i  ][Bnum+j  ].i += si*S1[k];

		k++;
	      }
	    }

	  }
	}
      }
    }

    /* spin-orbit coupling or LDA+U */
    else {  

      k = 0;
      for (AN=1; AN<=atomnum; AN++){
	GA_AN = order_GA[AN];
	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  Rn = ncn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  Bnum = MP[GB_AN];

	  l1 = atv_ijk[Rn][1];
	  l2 = atv_ijk[Rn][2];
	  l3 = atv_ijk[Rn][3];
	  kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

	  si = sin(2.0*PI*kRn);
	  co = cos(2.0*PI*kRn);

	  for (i=0; i<tnoA; i++){
	    for (j=0; j<tnoB; j++){

	      H[Anum+i  ][Bnum+j  ].r += co*RH0[k] - si*IH0[k];
	      H[Anum+i  ][Bnum+j  ].i += si*RH0[k] + co*IH0[k];

	      H[Anum+i+n][Bnum+j+n].r += co*RH1[k] - si*IH1[k];
	      H[Anum+i+n][Bnum+j+n].i += si*RH1[k] + co*IH1[k];
            
	      H[Anum+i  ][Bnum+j+n].r += co*RH2[k] - si*(RH3[k]+IH2[k]);
	      H[Anum+i  ][Bnum+j+n].i += si*RH2[k] + co*(RH3[k]+IH2[k]);

	      k++;

	    }

            if (SCF_iter==1 || all_knum!=1){

	      k -= tnoB; 

	      for (j=0; j<tnoB; j++){

		S[Anum+i  ][Bnum+j  ].r += co*S1[k];
		S[Anum+i  ][Bnum+j  ].i += si*S1[k];

		k++;
	      }
	    }

	  }
	}
      }
    }

    /* set off-diagonal part */

    for (i=1; i<=n; i++){
      for (j=1; j<=n; j++){
	H[j+n][i].r = H[i][j+n].r;
	H[j+n][i].i =-H[i][j+n].i;
      } 
    } 

    dtime(&Etime);
    time2 += Etime - Stime; 

    /* diagonalize S */

    dtime(&Stime);

    if (parallel_mode==0){
      EigenBand_lapack(S,ko,n,1);
    }
    else if (SCF_iter==1 || all_knum!=1){
      Eigen_PHH(MPI_CommWD1[myworld1],S,ko,n,n,1);
    }

    dtime(&Etime);
    time3 += Etime - Stime; 

    if (SCF_iter==1 || all_knum!=1){

      if (3<=level_stdout){
	printf(" myid0=%2d spin=%2d kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       myid0,spin,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
	for (i=1; i<=n; i++){
	  printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[i]);
	}
      }

      /* minus eigenvalues to 1.0e-14 */

      for (l=1; l<=n; l++){
	if (ko[l]<0.0) ko[l] = 1.0e-14;
	koS[l] = ko[l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) ko[l] = 1.0/sqrt(ko[l]);

      /* S * 1.0/sqrt(ko[l]) */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*ko[j1];
	  S[i1][j1].i = S[i1][j1].i*ko[j1];
	} 
      } 
    }

    /****************************************************
                  set H' and diagonalize it
    ****************************************************/

    /* U'^+ * H * U * M1 */

    /* transpose S */

    for (i1=1; i1<=n; i1++){
      for (j1=i1+1; j1<=n; j1++){
	Ctmp1 = S[i1][j1];
	Ctmp2 = S[j1][i1];
	S[i1][j1] = Ctmp2;
	S[j1][i1] = Ctmp1;
      }
    }

    /* for parallelization in the first world */

    if (all_knum==1){

      /* H * U' */
      /* C is distributed by row in each processor */

      dtime(&Stime);

#pragma omp parallel shared(C,S,H,n,is1,ie1,myid1) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=is1[myid1]+OMPID; j1<=ie1[myid1]; j1+=Nthrds){
	  for (i1=1; i1<=2*n; i1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    for (l=1; l<=n; l++){
	      sum_r0 += H[i1][l  ].r*S[j1][l].r - H[i1][l  ].i*S[j1][l].i;
	      sum_i0 += H[i1][l  ].r*S[j1][l].i + H[i1][l  ].i*S[j1][l].r;

	      sum_r1 += H[i1][n+l].r*S[j1][l].r - H[i1][n+l].i*S[j1][l].i;
	      sum_i1 += H[i1][n+l].r*S[j1][l].i + H[i1][n+l].i*S[j1][l].r;
	    }

	    C[2*j1-1][i1].r = sum_r0;
	    C[2*j1-1][i1].i = sum_i0;

	    C[2*j1  ][i1].r = sum_r1;
	    C[2*j1  ][i1].i = sum_i1;
	  }
	}     

      } /* #pragma omp parallel */

      dtime(&Etime);
      time4 += Etime - Stime; 

      /* U'^+ H * U' */
      /* H is distributed by row in each processor */

      dtime(&Stime);

#pragma omp parallel shared(C,S,H,n,is1,ie1,myid1) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r00,sum_i00,sum_r01,sum_i01,sum_r10,sum_i10,sum_r11,sum_i11,l,jj1,jj2)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=is1[myid1]+OMPID; j1<=ie1[myid1]; j1+=Nthrds){
	  for (i1=1; i1<=n; i1++){

	    sum_r00 = 0.0;
	    sum_i00 = 0.0;

	    sum_r01 = 0.0;
	    sum_i01 = 0.0;

	    sum_r10 = 0.0;
	    sum_i10 = 0.0;

	    sum_r11 = 0.0;
	    sum_i11 = 0.0;

	    jj1 = 2*j1 - 1;
	    jj2 = 2*j1;

	    for (l=1; l<=n; l++){

	      sum_r00 +=  S[i1][l].r*C[jj1][l  ].r + S[i1][l].i*C[jj1][l  ].i;
	      sum_i00 +=  S[i1][l].r*C[jj1][l  ].i - S[i1][l].i*C[jj1][l  ].r;

	      sum_r01 +=  S[i1][l].r*C[jj1][l+n].r + S[i1][l].i*C[jj1][l+n].i;
	      sum_i01 +=  S[i1][l].r*C[jj1][l+n].i - S[i1][l].i*C[jj1][l+n].r;

	      sum_r10 +=  S[i1][l].r*C[jj2][l  ].r + S[i1][l].i*C[jj2][l  ].i;
	      sum_i10 +=  S[i1][l].r*C[jj2][l  ].i - S[i1][l].i*C[jj2][l  ].r;

	      sum_r11 +=  S[i1][l].r*C[jj2][l+n].r + S[i1][l].i*C[jj2][l+n].i;
	      sum_i11 +=  S[i1][l].r*C[jj2][l+n].i - S[i1][l].i*C[jj2][l+n].r;
	    }

	    H[jj1][2*i1-1].r = sum_r00;
	    H[jj1][2*i1-1].i = sum_i00;

	    H[jj1][2*i1  ].r = sum_r01;
	    H[jj1][2*i1  ].i = sum_i01;

	    H[jj2][2*i1-1].r = sum_r10;
	    H[jj2][2*i1-1].i = sum_i10;

	    H[jj2][2*i1  ].r = sum_r11;
	    H[jj2][2*i1  ].i = sum_i11;

	  }
	}

      } /* #pragma omp parallel */

      /* broadcast H */

      BroadCast_ComplexMatrix(MPI_CommWD1[myworld1],H,2*n,is12,ie12,myid1,numprocs1,
			      stat_send,request_send,request_recv);

      dtime(&Etime);
      time5 += Etime - Stime; 

    }

    else{

      /* H * U' */

      dtime(&Stime);

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=n; j1+=Nthrds){
	  for (i1=1; i1<=2*n; i1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    for (l=1; l<=n; l++){
	      sum_r0 += H[i1][l  ].r*S[j1][l].r - H[i1][l  ].i*S[j1][l].i;
	      sum_i0 += H[i1][l  ].r*S[j1][l].i + H[i1][l  ].i*S[j1][l].r;

	      sum_r1 += H[i1][n+l].r*S[j1][l].r - H[i1][n+l].i*S[j1][l].i;
	      sum_i1 += H[i1][n+l].r*S[j1][l].i + H[i1][n+l].i*S[j1][l].r;
	    }

	    C[2*j1-1][i1].r = sum_r0;
	    C[2*j1-1][i1].i = sum_i0;

	    C[2*j1  ][i1].r = sum_r1;
	    C[2*j1  ][i1].i = sum_i1;
	  }
	}     

      } /* #pragma omp parallel */

      dtime(&Etime);
      time6 += Etime - Stime; 

      /* U'^+ H * U' */

      dtime(&Stime);

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r00,sum_i00,sum_r01,sum_i01,sum_r10,sum_i10,sum_r11,sum_i11,l,jj1,jj2)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=n; j1+=Nthrds){
	  for (i1=1; i1<=n; i1++){

	    sum_r00 = 0.0;
	    sum_i00 = 0.0;

	    sum_r01 = 0.0;
	    sum_i01 = 0.0;

	    sum_r10 = 0.0;
	    sum_i10 = 0.0;

	    sum_r11 = 0.0;
	    sum_i11 = 0.0;

	    jj1 = 2*j1 - 1;
	    jj2 = 2*j1;

	    for (l=1; l<=n; l++){

	      sum_r00 +=  S[i1][l].r*C[jj1][l  ].r + S[i1][l].i*C[jj1][l  ].i;
	      sum_i00 +=  S[i1][l].r*C[jj1][l  ].i - S[i1][l].i*C[jj1][l  ].r;

	      sum_r01 +=  S[i1][l].r*C[jj1][l+n].r + S[i1][l].i*C[jj1][l+n].i;
	      sum_i01 +=  S[i1][l].r*C[jj1][l+n].i - S[i1][l].i*C[jj1][l+n].r;

	      sum_r10 +=  S[i1][l].r*C[jj2][l  ].r + S[i1][l].i*C[jj2][l  ].i;
	      sum_i10 +=  S[i1][l].r*C[jj2][l  ].i - S[i1][l].i*C[jj2][l  ].r;

	      sum_r11 +=  S[i1][l].r*C[jj2][l+n].r + S[i1][l].i*C[jj2][l+n].i;
	      sum_i11 +=  S[i1][l].r*C[jj2][l+n].i - S[i1][l].i*C[jj2][l+n].r;
	    }

	    H[jj1][2*i1-1].r = sum_r00;
	    H[jj1][2*i1-1].i = sum_i00;

	    H[jj1][2*i1  ].r = sum_r01;
	    H[jj1][2*i1  ].i = sum_i01;

	    H[jj2][2*i1-1].r = sum_r10;
	    H[jj2][2*i1-1].i = sum_i10;

	    H[jj2][2*i1  ].r = sum_r11;
	    H[jj2][2*i1  ].i = sum_i11;

	  }
	}

      } /* #pragma omp parallel */

      dtime(&Etime);
      time7 += Etime - Stime; 

    }

    /* H to C (transposition) */

    for (i1=1; i1<=2*n; i1++){
      for (j1=1; j1<=2*n; j1++){
	C[j1][i1].r = H[i1][j1].r;
	C[j1][i1].i = H[i1][j1].i;
      }
    }

    /* penalty for ill-conditioning states */

    EV_cut0 = Threshold_OLP_Eigen;

    for (i1=1; i1<=n; i1++){

      if (koS[i1]<EV_cut0){
	C[2*i1-1][2*i1-1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	C[2*i1  ][2*i1  ].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
      }

      /* cutoff the interaction between the ill-conditioned state */

      if (1.0e+3<C[2*i1-1][2*i1-1].r){
	for (j1=1; j1<=2*n; j1++){
	  C[2*i1-1][j1    ] = Complex(0.0,0.0);
	  C[j1    ][2*i1-1] = Complex(0.0,0.0);
	  C[2*i1  ][j1    ] = Complex(0.0,0.0);
	  C[j1    ][2*i1  ] = Complex(0.0,0.0);
	}
	C[2*i1-1][2*i1-1] = Complex(1.0e+4,0.0);
	C[2*i1  ][2*i1  ] = Complex(1.0e+4,0.0);
      }
    }

    /* solve eigenvalue problem */

    dtime(&Stime);

    n1 = 2*n;

    if (parallel_mode==0){
      EigenBand_lapack(C, ko, 2*n, all_knum);
    }
    else{
      /*  The output C matrix is distributed by column. */
      Eigen_PHH(MPI_CommWD1[myworld1],C,ko,2*n,MaxN,0);
    }
 
    dtime(&Etime);
    time8 += Etime - Stime; 

    for (l=1; l<=MaxN; l++)  EIGEN[kloop][l] = ko[l];

    if (3<=level_stdout && 0<=kloop){
      printf(" myid0=%2d  kloop %i, k1 k2 k3 %10.6f %10.6f %10.6f\n",
              myid0,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
      for (i1=1; i1<=2*n; i1++){
	printf("  Eigenvalues of Kohn-Sham %2d  %15.12f\n", i1, ko[i1]);
      }
    }
  
    /* calculation of wave functions */

    if (all_knum==1){

      dtime(&Stime);

      /*  The H matrix is distributed by row */

      for (i1=1; i1<=2*n; i1++){
	for (j1=is2[myid1]; j1<=ie2[myid1]; j1++){
	  H[j1][i1] = C[i1][j1];
	}
      }

      /* transpose */

      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      /* C is distributed by row in each processor */

#pragma omp parallel shared(C,S,H,n,is2,ie2,myid1) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l,l1)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=is2[myid1]+OMPID; j1<=ie2[myid1]; j1+=Nthrds){
	  for (i1=1; i1<=n; i1++){

	    sum_r0 = 0.0; 
	    sum_i0 = 0.0;

	    sum_r1 = 0.0; 
	    sum_i1 = 0.0;

	    l1 = 0; 

	    for (l=1; l<=n; l++){

	      l1++; 

	      sum_r0 +=  S[i1][l].r*H[j1][l1].r - S[i1][l].i*H[j1][l1].i;
	      sum_i0 +=  S[i1][l].r*H[j1][l1].i + S[i1][l].i*H[j1][l1].r;

	      l1++; 

	      sum_r1 +=  S[i1][l].r*H[j1][l1].r - S[i1][l].i*H[j1][l1].i;
	      sum_i1 +=  S[i1][l].r*H[j1][l1].i + S[i1][l].i*H[j1][l1].r;
	    } 

	    C[j1][i1  ].r = sum_r0;
	    C[j1][i1  ].i = sum_i0;

	    C[j1][i1+n].r = sum_r1;
	    C[j1][i1+n].i = sum_i1;

	  }
	}

      } /* #pragma omp parallel */

      /* broadcast C:
       C is distributed by row in each processor
      */

      BroadCast_ComplexMatrix(MPI_CommWD1[myworld1],C,2*n,is2,ie2,myid1,numprocs1,
			      stat_send,request_send,request_recv);

      /* C to H (transposition)
         H consists of column vectors
      */ 

      for (i1=1; i1<=MaxN; i1++){
        for (j1=1; j1<=2*n; j1++){
	  H[j1][i1] = C[i1][j1];
	}
      }

      dtime(&Etime);
      time9 += Etime - Stime; 

    } /* if (all_knum==1) */
  } /* kloop0 */

  /****************************************************
     MPI:

     EIGEN
  ****************************************************/

  dtime(&Stime);

  for (kloop=0; kloop<T_knum; kloop++){
    ID = T_k_ID[kloop];
    MPI_Bcast(&EIGEN[kloop][0], MaxN+1, MPI_DOUBLE, ID, mpi_comm_level1);
  } 

  dtime(&Etime);
  time10 += Etime - Stime; 
  
  /****************************************************
           band energy in a finite temperature
  ****************************************************/

  dtime(&Stime);

  po = 0;
  loop_num = 0;
  ChemP_MAX = 10.0;  
  ChemP_MIN =-10.0;
    
  do {
    loop_num++;

    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);
    Num_State = 0.0;

    for (kloop=0; kloop<T_knum; kloop++){

      for (l=1; l<=2*n; l++){

	x = (EIGEN[kloop][l] - ChemP)*Beta;

	if      (x<=-x_cut)  FermiF = 1.0;
	else if (x>=x_cut)   FermiF = 0.0;
	else                 FermiF = 1.0/(1.0 + exp(x));

	Num_State += FermiF*(double)T_k_op[kloop];

      } /* l */
    } /* kloop */

    Num_State = Num_State/sum_weights;

    Dnum = TZ - Num_State - system_charge;
    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<10e-13) po = 1;
  }
  while (po==0 && loop_num<2000);

  Eele0[0] = 0.0;
  Eele0[1] = 0.0;

  for (kloop=0; kloop<T_knum; kloop++){
    for (l=1; l<=2*n; l++){
      x = (EIGEN[kloop][l] - ChemP)*Beta;
      if      (x<=-x_cut)  FermiF = 1.0;
      else if (x_cut<=x)   FermiF = 0.0;
      else                 FermiF = 1.0/(1.0 + exp(x));
      Eele0[0] += FermiF*EIGEN[kloop][l]*(double)T_k_op[kloop];
    } /* l */
  } /* kloop */

  Eele0[0] = Eele0[0]/sum_weights;
  Uele = Eele0[0];

  if (2<=level_stdout){
    printf("ChemP=%18.15f, Eele0[0]=%18.15f\n",ChemP,Eele0[0]);fflush(stdout);
  }

  dtime(&Etime);
  time11 += Etime - Stime; 

  /****************************************************
      if all_knum==1, calculate CDM, EDM, and iDM
  ****************************************************/

  if (all_knum==1){

    dtime(&Stime);

    kloop = S_knum;

    /* weight of k-point */

    kw = (double)T_k_op[kloop];

    /* initialize 
       CDM0, CDM1, CDM2, CDM3
       EDM0, EDM1, EDM2, EDM3
       iDM0, iDM1
    */

    for (i=0; i<size_H1; i++){
      CDM0[i] = 0.0;
      CDM1[i] = 0.0;
      CDM2[i] = 0.0;
      CDM3[i] = 0.0;

      EDM0[i] = 0.0;
      EDM1[i] = 0.0;
      EDM2[i] = 0.0;
      EDM3[i] = 0.0;

      iDM0[i] = 0.0;
      iDM1[i] = 0.0;
    }

    /* store Fermi function */

    po = 0;
    l = 1;
    do{

      eig = EIGEN[kloop][l];
      x = (eig - ChemP)*Beta;

      if      (x<-x_cut)  FermiF = 1.0;
      else if (x>x_cut)   FermiF = 0.0;
      else                FermiF = 1.0/(1.0 + exp(x));

      VecFkw[l]  = FermiF*kw; 
      VecFkwE[l] = VecFkw[l]*eig;

      if ( FermiF<=FermiEps && po==0 ) {
	lmax = l;  
	po = 1;
      }

      l++;

    } while(po==0 && l<=MaxN);

    if (po==0) lmax = MaxN;

    /* predetermination of k 
       H is used as temporal array 
    */ 

    k = 0;
    for (AN=1; AN<=atomnum; AN++){

      GA_AN = order_GA[AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      H[0][AN].r = (double)k;

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	GB_AN = natn[GA_AN][LB_AN];
	Rn = ncn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];
	k += tnoA*tnoB;
      }
    }

    /* calculate CDM, EDM, and iDM */

    k1 = T_KGrids1[kloop];
    k2 = T_KGrids2[kloop];
    k3 = T_KGrids3[kloop];

#pragma omp parallel shared(numprocs0,n,My_NZeros,SP_NZeros,MPI_CDM1_flag,CDM0,EDM0,iDM0,CDM1,EDM1,iDM1,CDM2,EDM2,CDM3,EDM3,VecFkwE,VecFkw,H,lmax,k1,k2,k3,atv_ijk,ncn,natn,FNAN,MP,Spe_Total_CNO,WhatSpecies,order_GA,atomnum) private(OMPID,Nthrds,Nprocs,AN,GA_AN,wanA,tnoA,Anum,k,LB_AN,GB_AN,RnB,wanB,tnoB,Bnum,l1,l2,l3,kRn,si,co,i,ia,ib,j,ja,jb,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,l,tmp,po,ID,RedumA,ImdumA,Redum2A,RedumB,ImdumB,Redum2B,RedumC,ImdumC,Redum2C,Imdum2D)
    { 

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (AN=1+OMPID; AN<=atomnum; AN+=Nthrds){
	GA_AN = order_GA[AN];
	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

	k = (int)H[0][AN].r;

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  RnB = ncn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  Bnum = MP[GB_AN];

	  l1 = atv_ijk[RnB][1];
	  l2 = atv_ijk[RnB][2];
	  l3 = atv_ijk[RnB][3];

	  kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

	  si = sin(2.0*PI*kRn);  
	  co = cos(2.0*PI*kRn);  

	  for (i=0; i<tnoA; i++){

	    ia = Anum + i;
	    ib = Anum + i + n;    
              
	    for (j=0; j<tnoB; j++){

	      ja = Bnum + j;
	      jb = Bnum + j + n;

	      po = 0;
	      ID = 0;
	      do { 
		if (MPI_CDM1_flag[ID] && SP_NZeros[ID]<=k && k<(SP_NZeros[ID]+My_NZeros[ID])) po = 1;
		ID++;
	      } while (po==0 && ID<numprocs0);

	      if (po==1){

		d1  = 0.0; /* CDM[0]    */
		d2  = 0.0; /* EDM[0]    */
		d3  = 0.0; /* iDM[0][0] */
		d4  = 0.0; /* CDM[1]    */
		d5  = 0.0; /* EDM[1]    */
		d6  = 0.0; /* iDM[0][1] */
		d7  = 0.0; /* CDM[2]    */
		d8  = 0.0; /* EDM[2]    */
		d9  = 0.0; /* CDM[3]    */
		d10 = 0.0; /* EDM[3]    */

		for (l=1; l<=lmax; l++){

		  /* Re11 */

		  RedumA = H[ia][l].r*H[ja][l].r + H[ia][l].i*H[ja][l].i;
		  ImdumA = H[ia][l].r*H[ja][l].i - H[ia][l].i*H[ja][l].r;
		  Redum2A = co*RedumA - si*ImdumA;

		  d1 += VecFkw[l]*Redum2A;
		  d2 += VecFkwE[l]*Redum2A;

		  /* Im11 */

		  d3 += VecFkw[l]*(si*RedumA + co*ImdumA);

		  /* Re22 */

		  RedumB = H[ib][l].r*H[jb][l].r + H[ib][l].i*H[jb][l].i;
		  ImdumB = H[ib][l].r*H[jb][l].i - H[ib][l].i*H[jb][l].r;
		  Redum2B = co*RedumB - si*ImdumB;

		  d4 += VecFkw[l]*Redum2B;
		  d5 += VecFkwE[l]*Redum2B;

		  /* Im22 */

		  d6 += VecFkw[l]*(si*RedumB + co*ImdumB);

		  /* Re12 */

		  RedumC = H[ia][l].r*H[jb][l].r + H[ia][l].i*H[jb][l].i;
		  ImdumC = H[ia][l].r*H[jb][l].i - H[ia][l].i*H[jb][l].r;
		  Redum2C = co*RedumC - si*ImdumC;

		  d7 += VecFkw[l]*Redum2C;
		  d8 += VecFkwE[l]*Redum2C;

		  /* Im12 */

		  Imdum2D = co*ImdumC + si*RedumC;

		  d9  += VecFkw[l]*Imdum2D;
		  d10 += VecFkwE[l]*Imdum2D;

		} /* l */

		CDM0[k] += d1;
		EDM0[k] += d2;
		iDM0[k] += d3;
		CDM1[k] += d4;
		EDM1[k] += d5;
		iDM1[k] += d6;
		CDM2[k] += d7;
		EDM2[k] += d8;
		CDM3[k] += d9;
		EDM3[k] += d10;

	      }
 
	      /* increment of k */

	      k++;

	    }
	  }
	} /* LB_AN */   
      } /* AN */     

    } /* #pragma omp parallel */

    /*******************************************************
       sum of
       CDM0, CDM1, CDM2, CDM3
       EDM0, EDM1, EDM2, EDM3
       iDM0, iDM1
       by Allreduce in MPI
    *******************************************************/

    dtime(&Etime);
    time12 += Etime - Stime; 

    /* CDM and iDM */

    dtime(&Stime);

    for (ID=0; ID<numprocs0; ID++){
      if (MPI_CDM1_flag[ID]){
        k = SP_NZeros[ID];

        MPI_Reduce(&CDM0[k], &RH0[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
        MPI_Reduce(&CDM1[k], &RH1[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
        MPI_Reduce(&CDM2[k], &RH2[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
        MPI_Reduce(&CDM3[k], &RH3[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
        MPI_Reduce(&iDM0[k], &IH0[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
        MPI_Reduce(&iDM1[k], &IH1[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
      }
    }

    dtime(&Etime);
    time13 += Etime - Stime; 

    dtime(&Stime);

    k = 0;
    for (AN=1; AN<=atomnum; AN++){
      GA_AN = order_GA[AN];
      MA_AN = F_G2M[GA_AN]; 

      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	GB_AN = natn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];

	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){

	    if (1<=MA_AN && MA_AN<=Matomnum){   
	      CDM[0][MA_AN][LB_AN][i][j] = RH0[k];
	      CDM[1][MA_AN][LB_AN][i][j] = RH1[k];
	      CDM[2][MA_AN][LB_AN][i][j] = RH2[k];
	      CDM[3][MA_AN][LB_AN][i][j] = RH3[k];
	      iDM[0][0][MA_AN][LB_AN][i][j] = IH0[k];
	      iDM[0][1][MA_AN][LB_AN][i][j] = IH1[k];
	    }

	    k++;

	  }
	}
      }
    }

    dtime(&Etime);
    time14 += Etime - Stime; 

    /* EDM */

    if (SCF_iter%3==1){ 

      dtime(&Stime);

      for (ID=0; ID<numprocs0; ID++){
	if (MPI_CDM1_flag[ID]){
	  k = SP_NZeros[ID];
	  MPI_Reduce(&EDM0[k], &RH0[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
	  MPI_Reduce(&EDM1[k], &RH1[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
	  MPI_Reduce(&EDM2[k], &RH2[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
	  MPI_Reduce(&EDM3[k], &RH3[k], My_NZeros[ID], MPI_DOUBLE, MPI_SUM, 0, MPI_CommWD_CDM1[ID]);
	}
      }

      dtime(&Etime);
      time15 += Etime - Stime; 

      dtime(&Stime);

      k = 0;
      for (AN=1; AN<=atomnum; AN++){
	GA_AN = order_GA[AN];
	MA_AN = F_G2M[GA_AN]; 

	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];

	  for (i=0; i<tnoA; i++){
	    for (j=0; j<tnoB; j++){

	      if (1<=MA_AN && MA_AN<=Matomnum){   
		EDM[0][MA_AN][LB_AN][i][j] = RH0[k];
		EDM[1][MA_AN][LB_AN][i][j] = RH1[k];
		EDM[2][MA_AN][LB_AN][i][j] = RH2[k];
		EDM[3][MA_AN][LB_AN][i][j] = RH3[k];
	      }

	      k++;

	    }
	  }
	}
      }

      dtime(&Etime);
      time16 += Etime - Stime; 
    }

  } /*if (all_knum==1) */

  dtime(&EiloopTime);
  if (myid0==Host_ID){
    printf("<Band_DFT>  Eigen, time=%lf\n", EiloopTime-SiloopTime);
  }
  SiloopTime=EiloopTime;
  fflush(stdout);

  /****************************************************
    if all_knum!=1, calculate 

        CDM[0]  Re alpha alpha density matrix
        CDM[1]  Re beta  beta  density matrix
        CDM[2]  Re alpha beta  density matrix
        CDM[3]  Im alpha beta  density matrix
        iDM[0][0]  Im alpha alpha density matrix
        iDM[0][1]  Im beta  beta  density matrix

        EDM[0]  Re alpha alpha energy density matrix
        EDM[1]  Re beta  beta  energy density matrix
        EDM[2]  Re alpha beta  energy density matrix
        EDM[3]  Im alpha beta  energy density matrix
  ****************************************************/

  dtime(&SiloopTime);

  if (all_knum!=1){

    /* initialize 
       CDM0, CDM1, CDM2, CDM3
       EDM0, EDM1, EDM2, EDM3
       iDM0, iDM1
    */

    for (i=0; i<size_H1; i++){
      CDM0[i] = 0.0;
      CDM1[i] = 0.0;
      CDM2[i] = 0.0;
      CDM3[i] = 0.0;

      EDM0[i] = 0.0;
      EDM1[i] = 0.0;
      EDM2[i] = 0.0;
      EDM3[i] = 0.0;

      iDM0[i] = 0.0;
      iDM1[i] = 0.0;
    }

    /* for kloop */

    for (kloop0=0; kloop0<num_kloop0; kloop0++){

      kloop = S_knum + kloop0;

      k1 = T_KGrids1[kloop];
      k2 = T_KGrids2[kloop];
      k3 = T_KGrids3[kloop];

      /* make S and H */

      for (i=1; i<=n; i++){
	for (j=1; j<=n; j++){
	  S[i][j] = Complex(0.0,0.0);
	} 
      } 

      for (i=1; i<=2*n; i++){
	for (j=1; j<=2*n; j++){
	  H[i][j] = Complex(0.0,0.0);
	} 
      } 

      /* non-spin-orbit coupling and non-LDA+U */

      if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
	  && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){  

	k = 0;
	for (AN=1; AN<=atomnum; AN++){
	  GA_AN = order_GA[AN];
	  wanA = WhatSpecies[GA_AN];
	  tnoA = Spe_Total_CNO[wanA];
	  Anum = MP[GA_AN];

	  for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	    GB_AN = natn[GA_AN][LB_AN];
	    Rn = ncn[GA_AN][LB_AN];
	    wanB = WhatSpecies[GB_AN];
	    tnoB = Spe_Total_CNO[wanB];
	    Bnum = MP[GB_AN];

	    l1 = atv_ijk[Rn][1];
	    l2 = atv_ijk[Rn][2];
	    l3 = atv_ijk[Rn][3];
	    kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

	    si = sin(2.0*PI*kRn);
	    co = cos(2.0*PI*kRn);

	    for (i=0; i<tnoA; i++){
	      for (j=0; j<tnoB; j++){

		S[Anum+i  ][Bnum+j  ].r += co*S1[k];
		S[Anum+i  ][Bnum+j  ].i += si*S1[k];

		H[Anum+i  ][Bnum+j  ].r += co*RH0[k];
		H[Anum+i  ][Bnum+j  ].i += si*RH0[k];

		H[Anum+i+n][Bnum+j+n].r += co*RH1[k];
		H[Anum+i+n][Bnum+j+n].i += si*RH1[k];
            
		H[Anum+i  ][Bnum+j+n].r += co*RH2[k] - si*RH3[k];
		H[Anum+i  ][Bnum+j+n].i += si*RH2[k] + co*RH3[k];

		k++;

	      }
	    }
	  }
	}
      }

      /* spin-orbit coupling or LDA+U */

      else {  

	k = 0;
	for (AN=1; AN<=atomnum; AN++){
	  GA_AN = order_GA[AN];
	  wanA = WhatSpecies[GA_AN];
	  tnoA = Spe_Total_CNO[wanA];
	  Anum = MP[GA_AN];

	  for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	    GB_AN = natn[GA_AN][LB_AN];
	    Rn = ncn[GA_AN][LB_AN];
	    wanB = WhatSpecies[GB_AN];
	    tnoB = Spe_Total_CNO[wanB];
	    Bnum = MP[GB_AN];

	    l1 = atv_ijk[Rn][1];
	    l2 = atv_ijk[Rn][2];
	    l3 = atv_ijk[Rn][3];
	    kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

	    si = sin(2.0*PI*kRn);
	    co = cos(2.0*PI*kRn);

	    for (i=0; i<tnoA; i++){
	      for (j=0; j<tnoB; j++){

		S[Anum+i  ][Bnum+j  ].r += co*S1[k];
		S[Anum+i  ][Bnum+j  ].i += si*S1[k];

		H[Anum+i  ][Bnum+j  ].r += co*RH0[k] - si*IH0[k];
		H[Anum+i  ][Bnum+j  ].i += si*RH0[k] + co*IH0[k];

		H[Anum+i+n][Bnum+j+n].r += co*RH1[k] - si*IH1[k];
		H[Anum+i+n][Bnum+j+n].i += si*RH1[k] + co*IH1[k];
            
		H[Anum+i  ][Bnum+j+n].r += co*RH2[k] - si*(RH3[k]+IH2[k]);
		H[Anum+i  ][Bnum+j+n].i += si*RH2[k] + co*(RH3[k]+IH2[k]);

		k++;

	      }
	    }
	  }
	}
      }

      /* set off-diagonal part */

      for (i=1; i<=n; i++){
	for (j=1; j<=n; j++){
	  H[j+n][i].r = H[i][j+n].r;
	  H[j+n][i].i =-H[i][j+n].i;
	} 
      } 

      /* diagonalize S */

      if (parallel_mode==0){
	EigenBand_lapack(S,ko,n,1);
      }
      else{
	Eigen_PHH(MPI_CommWD1[myworld1],S,ko,n,n,1);
      }

      if (3<=level_stdout){
	printf(" myid0=%2d spin=%2d kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       myid0,spin,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
	for (i=1; i<=n; i++){
	  printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[i]);
	}
      }

      /* minus eigenvalues to 1.0e-14 */

      for (l=1; l<=n; l++){
	if (ko[l]<0.0) ko[l] = 1.0e-14;
	koS[l] = ko[l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) ko[l] = 1.0/sqrt(ko[l]);

      /* S * 1.0/sqrt(ko[l]) */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*ko[j1];
	  S[i1][j1].i = S[i1][j1].i*ko[j1];
	} 
      } 

      /****************************************************
                  set H' and diagonalize it
      ****************************************************/

      /* U'^+ * H * U * M1 */

      /* transpose S */

      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      /* H * U' */

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=n; j1+=Nthrds){
	  for (i1=1; i1<=2*n; i1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    for (l=1; l<=n; l++){
	      sum_r0 += H[i1][l  ].r*S[j1][l].r - H[i1][l  ].i*S[j1][l].i;
	      sum_i0 += H[i1][l  ].r*S[j1][l].i + H[i1][l  ].i*S[j1][l].r;

	      sum_r1 += H[i1][n+l].r*S[j1][l].r - H[i1][n+l].i*S[j1][l].i;
	      sum_i1 += H[i1][n+l].r*S[j1][l].i + H[i1][n+l].i*S[j1][l].r;
	    }

	    C[2*j1-1][i1].r = sum_r0;
	    C[2*j1-1][i1].i = sum_i0;

	    C[2*j1  ][i1].r = sum_r1;
	    C[2*j1  ][i1].i = sum_i1;
	  }
	}     

      } /* #pragma omp parallel */

      /* U'^+ H * U' */

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (i1=1+OMPID; i1<=n; i1+=Nthrds){
	  for (j1=1; j1<=2*n; j1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    for (l=1; l<=n; l++){
	      sum_r0 +=  S[i1][l].r*C[j1][l  ].r + S[i1][l].i*C[j1][l  ].i;
	      sum_i0 +=  S[i1][l].r*C[j1][l  ].i - S[i1][l].i*C[j1][l  ].r;

	      sum_r1 +=  S[i1][l].r*C[j1][l+n].r + S[i1][l].i*C[j1][l+n].i;
	      sum_i1 +=  S[i1][l].r*C[j1][l+n].i - S[i1][l].i*C[j1][l+n].r;
	    }

	    H[2*i1-1][j1].r = sum_r0;
	    H[2*i1-1][j1].i = sum_i0;

	    H[2*i1  ][j1].r = sum_r1;
	    H[2*i1  ][j1].i = sum_i1;

	  }
	}

      } /* #pragma omp parallel */

      /* H to C */

      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
	  C[i1][j1].r = H[i1][j1].r;
	  C[i1][j1].i = H[i1][j1].i;
	}
      }

      /* penalty for ill-conditioning states */

      EV_cut0 = Threshold_OLP_Eigen;

      for (i1=1; i1<=n; i1++){

	if (koS[i1]<EV_cut0){
	  C[2*i1-1][2*i1-1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  C[2*i1  ][2*i1  ].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	}

	/* cutoff the interaction between the ill-conditioned state */

	if (1.0e+3<C[2*i1-1][2*i1-1].r){
	  for (j1=1; j1<=2*n; j1++){
	    C[2*i1-1][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1-1] = Complex(0.0,0.0);
	    C[2*i1  ][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1  ] = Complex(0.0,0.0);
	  }
	  C[2*i1-1][2*i1-1] = Complex(1.0e+4,0.0);
	  C[2*i1  ][2*i1  ] = Complex(1.0e+4,0.0);
	}
      }

      /* solve eigenvalue problem */

      dtime(&Stime);

      n1 = 2*n;

      if (parallel_mode==0){
	EigenBand_lapack(C, ko, 2*n, 1);
      }
      else{
        Eigen_PHH(MPI_CommWD1[myworld1],C,ko,2*n,MaxN,1);
      } 

      dtime(&Etime);
      time5 += Etime - Stime; 

      for (l=1; l<=MaxN; l++)  EIGEN[kloop][l] = ko[l];

      if (3<=level_stdout && 0<=kloop){
	printf(" myid0=%2d  kloop %i, k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       myid0,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
	for (i1=1; i1<=2*n; i1++){
	  printf("  Eigenvalues of Kohn-Sham %2d  %15.12f\n", i1, ko[i1]);
	}
      }

      /****************************************************
          transformation to the original eigenvectors
                 NOTE JRCAT-244p and JAIST-2122p 
                    C = U * lambda^{-1/2} * D
      ****************************************************/

      /* transpose */

      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      dtime(&Stime);

      /* transpose */

      for (i1=1; i1<=2*n; i1++){
	for (j1=i1+1; j1<=2*n; j1++){
	  Ctmp1 = C[i1][j1];
	  Ctmp2 = C[j1][i1];
	  C[i1][j1] = Ctmp2;
	  C[j1][i1] = Ctmp1;
	}
      }

#pragma omp parallel shared(C,S,H,n,n1) private(OMPID,Nthrds,Nprocs,i1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l,l1)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (i1=1+OMPID; i1<=n; i1+=Nthrds){
	  for (j1=1; j1<=n1; j1++){

	    sum_r0 = 0.0; 
	    sum_i0 = 0.0;

	    sum_r1 = 0.0; 
	    sum_i1 = 0.0;

	    l1 = 0; 

	    for (l=1; l<=n; l++){

	      l1++; 

	      sum_r0 +=  S[i1][l].r*C[j1][l1].r - S[i1][l].i*C[j1][l1].i;
	      sum_i0 +=  S[i1][l].r*C[j1][l1].i + S[i1][l].i*C[j1][l1].r;

	      l1++; 

	      sum_r1 +=  S[i1][l].r*C[j1][l1].r - S[i1][l].i*C[j1][l1].i;
	      sum_i1 +=  S[i1][l].r*C[j1][l1].i + S[i1][l].i*C[j1][l1].r;
	    } 

	    H[i1  ][j1].r = sum_r0;
	    H[i1  ][j1].i = sum_i0;

	    H[i1+n][j1].r = sum_r1;
	    H[i1+n][j1].i = sum_i1;
	  }
	}

      } /* #pragma omp parallel */

      dtime(&Etime);
      time6 += Etime - Stime; 

      /****************************************************
        calculate DM, EDM, and iDM
      ****************************************************/

      dtime(&Stime);

      /* weight of k-point */

      kw = (double)T_k_op[kloop];

      /* store Fermi function */

      po = 0;
      l = 1;
      do{

	eig = EIGEN[kloop][l];
	x = (eig - ChemP)*Beta;

	if      (x<-x_cut)  FermiF = 1.0;
	else if (x>x_cut)   FermiF = 0.0;
	else                FermiF = 1.0/(1.0 + exp(x));

        VecFkw[l]  = FermiF*kw; 
        VecFkwE[l] = VecFkw[l]*eig;

	if ( FermiF<=FermiEps && po==0 ) {
	  lmax = l;  
	  po = 1;
	}

        l++;

      } while(po==0 && l<=MaxN);

      if (po==0) lmax = MaxN;

      /* predetermination of k 
	 H is used as temporal array 
      */ 

      k = 0;
      for (AN=1; AN<=atomnum; AN++){

	GA_AN = order_GA[AN];
	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];

	H[0][AN].r = (double)k;

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  Rn = ncn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  k += tnoA*tnoB;
	}
      }

      /* DM, EDM, and iDM */

#pragma omp parallel shared(n,CDM0,EDM0,iDM0,CDM1,EDM1,iDM1,CDM2,EDM2,CDM3,EDM3,VecFkwE,VecFkw,H,lmax,k1,k2,k3,atv_ijk,ncn,natn,FNAN,MP,Spe_Total_CNO,WhatSpecies,order_GA,atomnum) private(OMPID,Nthrds,Nprocs,AN,GA_AN,wanA,tnoA,Anum,k,LB_AN,GB_AN,RnB,wanB,tnoB,Bnum,l1,l2,l3,kRn,si,co,i,ia,ib,j,ja,jb,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,l,RedumA,ImdumA,Redum2A,RedumB,ImdumB,Redum2B,RedumC,ImdumC,Redum2C,Imdum2D)
      { 

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (AN=1+OMPID; AN<=atomnum; AN+=Nthrds){

	  GA_AN = order_GA[AN];
	  wanA = WhatSpecies[GA_AN];
	  tnoA = Spe_Total_CNO[wanA];
	  Anum = MP[GA_AN];

	  k = (int)H[0][AN].r;

	  for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	    GB_AN = natn[GA_AN][LB_AN];
	    RnB = ncn[GA_AN][LB_AN];
	    wanB = WhatSpecies[GB_AN];
	    tnoB = Spe_Total_CNO[wanB];
	    Bnum = MP[GB_AN];

	    l1 = atv_ijk[RnB][1];
	    l2 = atv_ijk[RnB][2];
	    l3 = atv_ijk[RnB][3];

	    kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

	    si = sin(2.0*PI*kRn);  
	    co = cos(2.0*PI*kRn);  

	    for (i=0; i<tnoA; i++){

	      ia = Anum + i;
	      ib = Anum + i + n;    
              
	      for (j=0; j<tnoB; j++){

		ja = Bnum + j;
		jb = Bnum + j + n;

		d1  = 0.0; /* CDM[0]    */
		d2  = 0.0; /* EDM[0]    */
		d3  = 0.0; /* iDM[0][0] */
		d4  = 0.0; /* CDM[1]    */
		d5  = 0.0; /* EDM[1]    */
		d6  = 0.0; /* iDM[0][1] */
		d7  = 0.0; /* CDM[2]    */
		d8  = 0.0; /* EDM[2]    */
		d9  = 0.0; /* CDM[3]    */
		d10 = 0.0; /* EDM[3]    */

		for (l=1; l<=lmax; l++){

		  /* Re11 */

		  RedumA = H[ia][l].r*H[ja][l].r + H[ia][l].i*H[ja][l].i;
		  ImdumA = H[ia][l].r*H[ja][l].i - H[ia][l].i*H[ja][l].r;
		  Redum2A = co*RedumA - si*ImdumA;

		  d1 += VecFkw[l]*Redum2A;
		  d2 += VecFkwE[l]*Redum2A;

		  /* Im11 */

		  d3 += VecFkw[l]*(si*RedumA + co*ImdumA);

		  /* Re22 */

		  RedumB = H[ib][l].r*H[jb][l].r + H[ib][l].i*H[jb][l].i;
		  ImdumB = H[ib][l].r*H[jb][l].i - H[ib][l].i*H[jb][l].r;
		  Redum2B = co*RedumB - si*ImdumB;

		  d4 += VecFkw[l]*Redum2B;
		  d5 += VecFkwE[l]*Redum2B;

		  /* Im22 */

		  d6 += VecFkw[l]*(si*RedumB + co*ImdumB);

		  /* Re12 */

		  RedumC = H[ia][l].r*H[jb][l].r + H[ia][l].i*H[jb][l].i;
		  ImdumC = H[ia][l].r*H[jb][l].i - H[ia][l].i*H[jb][l].r;
		  Redum2C = co*RedumC - si*ImdumC;

		  d7 += VecFkw[l]*Redum2C;
		  d8 += VecFkwE[l]*Redum2C;

		  /* Im12 */

		  Imdum2D = co*ImdumC + si*RedumC;

		  d9  += VecFkw[l]*Imdum2D;
		  d10 += VecFkwE[l]*Imdum2D;

		}

		CDM0[k] += d1;
		EDM0[k] += d2;
		iDM0[k] += d3;
		CDM1[k] += d4;
		EDM1[k] += d5;
		iDM1[k] += d6;
		CDM2[k] += d7;
		EDM2[k] += d8;
		CDM3[k] += d9;
		EDM3[k] += d10;
 
		/* increment of k */

		k++;

	      }
	    }
	  }   
	}     

      } /* #pragma omp parallel */

      dtime(&Etime);
      time7 += Etime - Stime; 

    } /* kloop0 */

    /*******************************************************
       sum of
       CDM0, CDM1, CDM2, CDM3
       EDM0, EDM1, EDM2, EDM3
       iDM0, iDM1
       by Allreduce in MPI
    *******************************************************/

    if (parallel_mode){
      tmp = 1.0/(double)numprocs1;
      for (i=0; i<size_H1; i++){
	CDM0[i] *= tmp;
	CDM1[i] *= tmp;
	CDM2[i] *= tmp;
	CDM3[i] *= tmp;
	EDM0[i] *= tmp;
	EDM1[i] *= tmp;
	EDM2[i] *= tmp;
	EDM3[i] *= tmp;
	iDM0[i] *= tmp;
	iDM1[i] *= tmp;
      }
    }

    /* CDM and iDM */

    MPI_Allreduce(&CDM0[0], &RH0[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&CDM1[0], &RH1[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&CDM2[0], &RH2[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&CDM3[0], &RH3[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&iDM0[0], &IH0[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&iDM1[0], &IH1[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    k = 0;
    for (AN=1; AN<=atomnum; AN++){
      GA_AN = order_GA[AN];
      MA_AN = F_G2M[GA_AN]; 

      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	GB_AN = natn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];

	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){

	    if (1<=MA_AN && MA_AN<=Matomnum){   
	      CDM[0][MA_AN][LB_AN][i][j] = RH0[k];
	      CDM[1][MA_AN][LB_AN][i][j] = RH1[k];
	      CDM[2][MA_AN][LB_AN][i][j] = RH2[k];
	      CDM[3][MA_AN][LB_AN][i][j] = RH3[k];
	      iDM[0][0][MA_AN][LB_AN][i][j] = IH0[k];
	      iDM[0][1][MA_AN][LB_AN][i][j] = IH1[k];
	    }

	    k++;

	  }
	}
      }
    }

    /* EDM */

    MPI_Allreduce(&EDM0[0], &RH0[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&EDM1[0], &RH1[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&EDM2[0], &RH2[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    MPI_Allreduce(&EDM3[0], &RH3[0], size_H1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    k = 0;
    for (AN=1; AN<=atomnum; AN++){
      GA_AN = order_GA[AN];
      MA_AN = F_G2M[GA_AN]; 

      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	GB_AN = natn[GA_AN][LB_AN];
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];

	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){

	    if (1<=MA_AN && MA_AN<=Matomnum){   
	      EDM[0][MA_AN][LB_AN][i][j] = RH0[k];
	      EDM[1][MA_AN][LB_AN][i][j] = RH1[k];
	      EDM[2][MA_AN][LB_AN][i][j] = RH2[k];
	      EDM[3][MA_AN][LB_AN][i][j] = RH3[k];
	    }

	    k++;

	  }
	}
      }
    }

  } /* if (all_knum!=1) */

  dtime(&EiloopTime);
  if (myid0==Host_ID){
    printf("<Band_DFT>  DM, time=%lf\n", EiloopTime-SiloopTime);
  }
  SiloopTime=EiloopTime;
  fflush(stdout);

  /********************************************
      normalization of CDM, EDM, and iDM
  ********************************************/

  dum = 1.0/sum_weights;

  for (spin=0; spin<=SpinP_switch; spin++){
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];    
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];
        Bnum = MP[GB_AN];

        for (i=0; i<tnoA; i++){
          for (j=0; j<tnoB; j++){
            CDM[spin][MA_AN][LB_AN][i][j] = CDM[spin][MA_AN][LB_AN][i][j]*dum;
            EDM[spin][MA_AN][LB_AN][i][j] = EDM[spin][MA_AN][LB_AN][i][j]*dum;
          }
        }

        if (spin==0){
          for (i=0; i<tnoA; i++){
            for (j=0; j<tnoB; j++){
              iDM[0][0][MA_AN][LB_AN][i][j] = iDM[0][0][MA_AN][LB_AN][i][j]*dum;
              iDM[0][1][MA_AN][LB_AN][i][j] = iDM[0][1][MA_AN][LB_AN][i][j]*dum;
            }
          }
	}
      }
    }
  }

  /****************************************************
                     bond energies
  ****************************************************/

  My_Eele1[0] = 0.0;
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];    
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];

    for (j=0; j<=FNAN[GA_AN]; j++){
      GB_AN = natn[GA_AN][j];  
      wanB = WhatSpecies[GB_AN];
      tnoB = Spe_Total_CNO[wanB];

      /* non-spin-orbit coupling and non-LDA+U */  
      if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
          && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
        for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
            My_Eele1[0] = My_Eele1[0]
                         + CDM[0][MA_AN][j][k][l]*nh[0][MA_AN][j][k][l]
                         + CDM[1][MA_AN][j][k][l]*nh[1][MA_AN][j][k][l]
                         + 2.0*CDM[2][MA_AN][j][k][l]*nh[2][MA_AN][j][k][l]
                         - 2.0*CDM[3][MA_AN][j][k][l]*nh[3][MA_AN][j][k][l];
	  }
        }
      }

      /* spin-orbit coupling or LDA+U */  
      else {
        for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
            My_Eele1[0] = My_Eele1[0]
                         + CDM[0][MA_AN][j][k][l]*nh[0][MA_AN][j][k][l]
                         - iDM[0][0][MA_AN][j][k][l]*ImNL[0][MA_AN][j][k][l]
                         + CDM[1][MA_AN][j][k][l]*nh[1][MA_AN][j][k][l]
                         - iDM[0][1][MA_AN][j][k][l]*ImNL[1][MA_AN][j][k][l]
                         + 2.0*CDM[2][MA_AN][j][k][l]*nh[2][MA_AN][j][k][l]
                         - 2.0*CDM[3][MA_AN][j][k][l]*(nh[3][MA_AN][j][k][l]
                                                    +ImNL[2][MA_AN][j][k][l]);
	  }
        }
      }
  
    }
  }

  /* MPI, My_Eele1 */

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_Eele1[0], &Eele1[0], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  Eele1[1] = 0.0;

  /*
  printf("Eele00=%15.12f Eele01=%15.12f\n",Eele0[0],Eele0[1]);
  printf("Eele10=%15.12f Eele11=%15.12f\n",Eele1[0],Eele1[1]);
  MPI_Finalize();
  exit(0);
  */

  /****************************************************
                        Output
  ****************************************************/

  if (myid0==Host_ID){
  
    strcpy(file_EV,".EV");
    fnjoint(filepath,filename,file_EV);

    if ((fp_EV = fopen(file_EV,"w")) != NULL){

#ifdef xt3
      setvbuf(fp_EV,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp_EV,"\n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"           Eigenvalues (Hartree) of SCF KS-eq.           \n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"***********************************************************\n\n");
      fprintf(fp_EV,"   Chemical Potential (Hatree) = %18.14f\n",ChemP);
      fprintf(fp_EV,"   Number of States            = %18.14f\n",Num_State);
      fprintf(fp_EV,"   Eigenvalues\n\n");

      for (kloop=0; kloop<T_knum; kloop++){

	if (0<T_k_op[kloop]){

	  fprintf(fp_EV,"\n");
	  fprintf(fp_EV,"   kloop=%i\n",kloop);
	  fprintf(fp_EV,"   k1=%10.5f k2=%10.5f k3=%10.5f\n\n",
		  T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);

	  for (l=1; l<=MaxN; l++){
	    fprintf(fp_EV,"%5d  %18.14f\n",l,EIGEN[kloop][l]);
	  }
	}
      }
      fclose(fp_EV);
    }
    else{
      printf("Failure of saving the EV file.\n");
      fclose(fp_EV);
    }  
  }

  /****************************************************
                       free arrays
  ****************************************************/

  free(VecFkw);
  free(VecFkwE);

  for (j=0; j<n2; j++){
    free(H[j]);
  }
  free(H);

  for (j=0; j<n2; j++){
    free(C[j]);
  }
  free(C);

  free(KGrids1); free(KGrids2); free(KGrids3);

  for (i=0;i<knum_i; i++) {
    for (j=0;j<knum_j; j++) {
      free(k_op[i][j]);
    }
    free(k_op[i]);
  }
  free(k_op);

  free(T_KGrids1);
  free(T_KGrids2);
  free(T_KGrids3);
  free(T_k_op);
  free(T_k_ID);

  for (j=0; j<T_knum; j++){
    free(EIGEN[j]);
  }
  free(EIGEN);

  free(MP);
  free(order_GA);

  free(My_NZeros);
  free(SP_NZeros);
  free(SP_Atoms);

  free(ko);

  if (all_knum==1){
    for (ID=0; ID<numprocs0; ID++){
      if (MPI_CDM1_flag[ID]){
        MPI_Comm_free(&MPI_CommWD_CDM1[ID]);
      }
    }
  }

  free(MPI_CommWD_CDM1);
  free(MPI_CDM1_flag);

  free(S1);
  free(RH0);
  free(RH1);
  free(RH2);
  free(RH3);

  free(IH0);
  free(IH1);
  free(IH2);

  free(CDM0);
  free(CDM1);
  free(CDM2);
  free(CDM3);

  free(EDM0);
  free(EDM1);
  free(EDM2);
  free(EDM3);

  free(iDM0);
  free(iDM1);

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

  if (all_knum==1){

    free(stat_send);
    free(request_send);
    free(request_recv);

    free(is1);
    free(ie1);
    free(is12);
    free(ie12);
    free(is2);
    free(ie2);
  }

  /* for timing */

  if (measure_time){
    printf("myid0=%2d  time1 =%10.5f\n",myid0,time1);fflush(stdout);
    printf("myid0=%2d  time2 =%10.5f\n",myid0,time2);fflush(stdout);
    printf("myid0=%2d  time3 =%10.5f\n",myid0,time3);fflush(stdout);
    printf("myid0=%2d  time4 =%10.5f\n",myid0,time4);fflush(stdout);
    printf("myid0=%2d  time5 =%10.5f\n",myid0,time5);fflush(stdout);
    printf("myid0=%2d  time6 =%10.5f\n",myid0,time6);fflush(stdout);
    printf("myid0=%2d  time7 =%10.5f\n",myid0,time7);fflush(stdout);
    printf("myid0=%2d  time8 =%10.5f\n",myid0,time8);fflush(stdout);
    printf("myid0=%2d  time9 =%10.5f\n",myid0,time9);fflush(stdout);
    printf("myid0=%2d  time10=%10.5f\n",myid0,time10);fflush(stdout);
    printf("myid0=%2d  time11=%10.5f\n",myid0,time11);fflush(stdout);
    printf("myid0=%2d  time12=%10.5f\n",myid0,time12);fflush(stdout);
    printf("myid0=%2d  time13=%10.5f\n",myid0,time13);fflush(stdout);
    printf("myid0=%2d  time14=%10.5f\n",myid0,time14);fflush(stdout);
    printf("myid0=%2d  time15=%10.5f\n",myid0,time15);fflush(stdout);
    printf("myid0=%2d  time16=%10.5f\n",myid0,time16);fflush(stdout);
  }

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}




static double Band_non_collinear(
             int knum_i, int knum_j, int knum_k,
             int SpinP_switch,
             double *****nh,
             double *****ImNL,
             double ****CntOLP,
             double *****CDM,
             double *****EDM,
             double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  int i,j,k,l,n,n2,wan;
  int ii1,jj1,m,mn,lmn,lmax;
  int *MP,*arpo;
  int i1,j1,po,spin,n1;
  int ia,ib,ja,jb;
  int num2,RnB,l1,l2,l3,loop_num;
  int ct_AN,h_AN,wanA,tnoA,wanB,tnoB;
  int MA_AN,GA_AN,Anum;
  double time0;
  double Imdum2,Imsum,ImsumE;
  int LB_AN,GB_AN,Bnum;
  double EV_cut0;
  double snum_i,snum_j,snum_k,k1,k2,k3,sum,sumi,Num_State,FermiF;
  double x,Dnum,Dnum2,AcP,ChemP_MAX,ChemP_MIN;
  double *ko,*M1,**EIGEN,*tmp_array0;
  double *koS;
  double *****Dummy_ImNL;
  double d1,d2,d3,d4,d5,d6,d7,d8,d9,d10;
  double *VecFermiF;
  dcomplex **H,**S,**C,**H2,**TmpM;
  dcomplex Ctmp1,Ctmp2;
  int ***k_op,*T_k_op;
  int ii,ij,ik;
  int S_knum,E_knum,T_knum,e1,s1,kloop,kloop0,num_kloop0;
  double *KGrids1,*KGrids2,*KGrids3;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  double u2,v2,uv,vu,tmp;
  double My_Eele1[2]; 
  double TmpEigen,Tmpk_op;
  double TZ,dum,sumE,kRn,si,co;
  double Resum,ResumE,Redum,Redum2,Imdum;
  double TStime,TEtime,SiloopTime,EiloopTime;
  double FermiEps=1.0e-14;
  double x_cut=30.0;
  char file_EV[YOUSO10];
  FILE *fp_EV;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID,ID1;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /*
  printf(" Band_DFT start ID=%2d\n",myid);fflush(stdout);
  */

  dtime(&TStime);

  /****************************************************
             calculation of the array size
  ****************************************************/

  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n  = n + Spe_Total_CNO[wanA];
  }
  n2 = 2*n + 2;

  /****************************************************
   Allocation

   int       MP[List_YOUSO[1]]
   double    ko[n2]
   double    koS[n+1];
   dcomplex  H[n2][n2]
   dcomplex  S[n2][n2]
   double    M1[n2]
   dcomplex  C[n2][n2]
   double    EIGEN[List_YOUSO[27]]
                  [List_YOUSO[28]]
                  [List_YOUSO[29]]
                  [n2] 
   double    KGrids1[List_YOUSO[27]]
   double    KGrids2[List_YOUSO[28]]
   double    KGrids3[List_YOUSO[29]]
   int       k_op[knum_i][knum_j][knum_k]
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  ko = (double*)malloc(sizeof(double)*n2);
  koS = (double*)malloc(sizeof(double)*(n+1));
  arpo = (int*)malloc(sizeof(int)*numprocs);

  tmp_array0 = (double*)malloc(sizeof(double)*2*n2*n2);

  H = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    H[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  H2 = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    H2[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  S = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (i=0; i<n2; i++){
    S[i] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  M1 = (double*)malloc(sizeof(double)*n2);

  VecFermiF = (double*)malloc(sizeof(double)*n2);

  C = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    C[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  KGrids1 = (double*)malloc(sizeof(double)*List_YOUSO[27]);
  KGrids2 = (double*)malloc(sizeof(double)*List_YOUSO[28]);
  KGrids3 = (double*)malloc(sizeof(double)*List_YOUSO[29]);

  k_op=(int***)malloc(sizeof(int**)*(knum_i));
  for (i=0;i<knum_i; i++) {
    k_op[i]=(int**)malloc(sizeof(int*)*(knum_j));
    for (j=0;j<knum_j; j++) {
      k_op[i][j]=(int*)malloc(sizeof(int)*(knum_k));
    }
  }
  /*  k_op[0:knum_i-1][0:knum_j-1][0:knum_k-1] */

  /* non-spin-orbit coupling and non-LDA+U */
  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
      && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){

    Dummy_ImNL = (double*****)malloc(sizeof(double****)*1);
    Dummy_ImNL[0] = (double****)malloc(sizeof(double***)*1);
    Dummy_ImNL[0][0] = (double***)malloc(sizeof(double**)*1);
    Dummy_ImNL[0][0][0] = (double**)malloc(sizeof(double*)*1);
    Dummy_ImNL[0][0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  /****************************************************
    PrintMemory
  ****************************************************/

  if (firsttime){

    PrintMemory("Band_DFT: ko",sizeof(double)*n2,NULL);
    PrintMemory("Band_DFT: tmp_array0",sizeof(double)*n2*n2*2,NULL);
    PrintMemory("Band_DFT: H",sizeof(dcomplex)*n2*n2,NULL);
    PrintMemory("Band_DFT: S",sizeof(dcomplex)*n2*n2,NULL);
    PrintMemory("Band_DFT: C",sizeof(dcomplex)*n2*n2,NULL);
    PrintMemory("Band_DFT: M1",sizeof(double)*n2,NULL);

    PrintMemory("Band_DFT: EIGEN",
          sizeof(double)*List_YOUSO[27]*
                         List_YOUSO[28]*
                         List_YOUSO[29]*n2,NULL);
    PrintMemory("Band_DFT: ko_op",
          sizeof(int)*knum_i*knum_j*knum_k,NULL);
  }

  /* for PrintMemory */
  firsttime=0;

  snum_i = (double)knum_i;
  snum_j = (double)knum_j;
  snum_k = (double)knum_k;

  for (i=0; i<=(knum_i-1); i++){
    if (knum_i==1){
      k1 = 0.0;
    }
    else {
      k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_i) + Shift_K_Point;
    }
    KGrids1[i]=k1;
  }
  for (i=0; i<=(knum_j-1); i++){
    if (knum_j==1){
      k1 = 0.0;
    }
    else {
      k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_j) - Shift_K_Point;
    }
    KGrids2[i]=k1;
  }
  for (i=0; i<=(knum_k-1); i++){
    if (knum_k==1){
      k1 = 0.0;
    }
    else {
      k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_k) + 2.0*Shift_K_Point;
    }
    KGrids3[i]=k1;
  }

  if (myid==Host_ID){
    printf(" KGrids1: ");
    for (i=0;i<=knum_i-1;i++) printf("%9.5f ",KGrids1[i]);
    printf("\n");
    printf(" KGrids2: ");
    for (i=0;i<=knum_j-1;i++) printf("%9.5f ",KGrids2[i]);
    printf("\n");
    printf(" KGrids3: ");
    for (i=0;i<=knum_k-1;i++) printf("%9.5f ",KGrids3[i]);
    printf("\n");
  }
  /**************************************************************
   k_op[i][j][k]: weight of DOS 
                 =0   no calc.
                 =1   G-point
             note that there is no simple inversion symmetry for
               the wave functions in noncollinear case 
                 =1  other  
  *************************************************************/

  for (i=0;i<=knum_i-1;i++) {
    for (j=0;j<=knum_j-1;j++) {
      for (k=0;k<=knum_k-1;k++) {
	k_op[i][j][k]=-999;
      }
    }
  }
  for (i=0;i<=knum_i-1;i++) {
    for (j=0;j<=knum_j-1;j++) {
      for (k=0;k<=knum_k-1;k++) {
	if ( k_op[i][j][k]==-999 ) {
	  k_inversion(i,j,k,knum_i,knum_j,knum_k,&ii,&ij,&ik);
	  if ( i==ii && j==ij && k==ik ) {
	    /*
	      printf("* : %d %d %d (%f %f %f)\n",i,j,k,
	      KGrids1[i], KGrids2[j], KGrids3[k]);
	    */
	    k_op[i][j][k]    = 1;
	  }

          else {

            /* note that there is no simple inversion symmetry for
               the wave functions in noncollinear case */

	    k_op[i][j][k]    = 1;
	    k_op[ii][ij][ik] = 1;
	  }
	}
      } /* k */
    } /* j */
  } /* i */

  /***********************************
       one-dimentionalize for MPI
  ************************************/

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        if (0<k_op[i][j][k]){
          T_knum++;
        }
      }
    }
  }

  /* allocation of arrays */

  T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);
  T_k_op    = (int*)malloc(sizeof(int)*T_knum);

  EIGEN = (double**)malloc(sizeof(double*)*T_knum);
  for (j=0; j<T_knum; j++){
    EIGEN[j] = (double*)malloc(sizeof(double)*n2);
    for (i=0; i<n2; i++)  EIGEN[j][i] = 1.0e+5;
  }

  /* set T_KGrid1,2,3 and T_k_op */

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        if (0<k_op[i][j][k]){

	  T_KGrids1[T_knum] = KGrids1[i];
	  T_KGrids2[T_knum] = KGrids2[j];
	  T_KGrids3[T_knum] = KGrids3[k];

	  T_k_op[T_knum]    = k_op[i][j][k];

          T_knum++;
        }
      }
    }
  }

  /* allocate k-points into proccessors */

  if (T_knum<=myid){
    S_knum = -10;
    E_knum = -100;
    num_kloop0 = 1;
  }
  else if (T_knum<numprocs) {
    S_knum = myid;
    E_knum = myid;
    num_kloop0 = 1;
  }
  else {
    tmp = (double)T_knum/(double)numprocs; 
    num_kloop0 = (int)tmp + 1;

    S_knum = (int)((double)myid*(tmp+0.0001)); 
    E_knum = (int)((double)(myid+1)*(tmp+0.0001)) - 1;
    if (myid==(numprocs-1)) E_knum = T_knum - 1;
    if (E_knum<0)           E_knum = 0;
  }

  /* find TZ */

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ = TZ + Spe_Core_Charge[wan];
  }

  /****************************************************
                      start kloop
  ****************************************************/

  dtime(&SiloopTime);

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S and diagonalize it */

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

        k1 = T_KGrids1[kloop];
        k2 = T_KGrids2[kloop];
        k3 = T_KGrids3[kloop];

        Overlap_Band(ID,CntOLP,TmpM,MP,k1,k2,k3);
        n = TmpM[0][0].r;

        if (myid==ID){
	  for (i1=1; i1<=n; i1++){
	    for (j1=1; j1<=n; j1++){
	      S[i1][j1].r = TmpM[i1][j1].r;
	      S[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 

      }
    } /* ID */

    kloop = arpo[myid];

    if (0<=kloop){

      EigenBand_lapack(S,ko,n,1);

      if (3<=level_stdout){
	printf(" myid=%2d kloop %2d  k1 k2 k3 %9.5f %9.5f %9.5f\n",
	       myid,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
        for (i1=1; i1<=n; i1++){
          printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[i1]);
        }
      }

      /* minus eigenvalues to 1.0e-14 */
      for (l=1; l<=n; l++){
        if (ko[l]<0.0){
          ko[l] = 1.0e-14;

          if (2<=level_stdout){
  	    printf("found an eigenvalue smaller than %10.8f of OLP kloop=%2d\n",
                           Threshold_OLP_Eigen,kloop);
	  }
	}

        koS[l] = ko[l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[l]);

      /* S * M1  */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*M1[j1];
	  S[i1][j1].i = S[i1][j1].i*M1[j1];
	} 
      } 

    } /* if (0<=kloop) */

    /* set H' and diagonalize it */

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

        k1 = T_KGrids1[kloop];
        k2 = T_KGrids2[kloop];
        k3 = T_KGrids3[kloop];

        /****************************************************
                   make a full Hamiltonian matrix
        ****************************************************/

        /* non-spin-orbit coupling and non-LDA+U */
        if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
            && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0)
          Hamiltonian_Band_NC(ID, nh, Dummy_ImNL, TmpM, MP, k1, k2, k3);

        /* spin-orbit coupling or LDA+U */
        else
          Hamiltonian_Band_NC(ID, nh,       ImNL, TmpM, MP, k1, k2, k3);

	if (myid==ID){
	  for (i1=1; i1<=2*n; i1++){
	    for (j1=1; j1<=2*n; j1++){
	      H[i1][j1].r = TmpM[i1][j1].r;
	      H[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
	} 

      } /* if (0<=kloop) */
    } /* ID */

    kloop = arpo[myid];

    if (0<=kloop){

      /****************************************************
                     M1 * U^+ * H * U * M1
      ****************************************************/

      /* transpose S */

      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      /* H * U * M1 */

      for (j1=1; j1<=n; j1++){
        for (i1=1; i1<=2*n; i1++){
	  for (m=0; m<=1; m++){

	    sum  = 0.0;
	    sumi = 0.0;
            mn = m*n;
	    for (l=1; l<=n; l++){
	      sum  += H[i1][l+mn].r*S[j1][l].r - H[i1][l+mn].i*S[j1][l].i;
	      sumi += H[i1][l+mn].r*S[j1][l].i + H[i1][l+mn].i*S[j1][l].r;
	    }

	    jj1 = 2*j1 - 1 + m;

	    C[jj1][i1].r = sum;
	    C[jj1][i1].i = sumi;
	  }
	}
      }     

      /* M1 * U^+ H * U * M1 */

      for (i1=1; i1<=n; i1++){

	for (m=0; m<=1; m++){

	  ii1 = 2*i1 - 1 + m;

	  for (j1=1; j1<=2*n; j1++){
	    sum  = 0.0;
	    sumi = 0.0;
            mn = m*n;
	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*C[j1][l+mn].r + S[i1][l].i*C[j1][l+mn].i;
	      sumi +=  S[i1][l].r*C[j1][l+mn].i - S[i1][l].i*C[j1][l+mn].r;
	    }
	    H[ii1][j1].r = sum;
	    H[ii1][j1].i = sumi;
	  }
	}
      }     

      /* H to C */

      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
	  C[i1][j1].r = H[i1][j1].r;
	  C[i1][j1].i = H[i1][j1].i;
	}
      }

      /* penalty for ill-conditioning states */

      EV_cut0 = Threshold_OLP_Eigen;

      for (i1=1; i1<=n; i1++){

	if (koS[i1]<EV_cut0){
	  C[2*i1-1][2*i1-1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  C[2*i1  ][2*i1  ].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	}

	/* cutoff the interaction between the ill-conditioned state */

	if (1.0e+3<C[2*i1-1][2*i1-1].r){
	  for (j1=1; j1<=2*n; j1++){
	    C[2*i1-1][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1-1] = Complex(0.0,0.0);
	    C[2*i1  ][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1  ] = Complex(0.0,0.0);
	  }
	  C[2*i1-1][2*i1-1] = Complex(1.0e+4,0.0);
	  C[2*i1  ][2*i1  ] = Complex(1.0e+4,0.0);
	}
      }

      /* solve eigenvalue problem */

      n1 = 2*n;
      EigenBand_lapack(C, ko, n1, 1);

      for (l=1; l<=n1; l++)  EIGEN[kloop][l] = ko[l];

    } /* if (0<=kloop) */ 

    if (2<=level_stdout && 0<=kloop){
      printf(" myid=%2d  kloop %i, k1 k2 k3 %10.6f %10.6f %10.6f\n",
              myid,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
      for (i1=1; i1<=2*n; i1++){
	printf("  Eigenvalues of Kohn-Sham %2d  %15.12f\n", i1, ko[i1]);
      }
    }
  
  } /* kloop0 */

  /****************************************************
     MPI:

     EIGEN
  ****************************************************/

  tmp = (double)T_knum/(double)numprocs; 

  for (kloop=0; kloop<T_knum; kloop++){

    for (ID=0; ID<numprocs; ID++){

      if (T_knum<=ID){
	s1 = -10;
	e1 = -100;
      }
      else if (T_knum<numprocs) {
	s1 = ID;
	e1 = ID;
      }
      else {
	s1 = (int)((double)ID*(tmp+0.0001)); 
	e1 = (int)((double)(ID+1)*(tmp+0.0001)) - 1;
	if (ID==(numprocs-1)) e1 = T_knum - 1;
	if (e1<0)             e1 = 0;
      }

      if (s1<=kloop && kloop<=e1)  ID1 = ID;
    }

    MPI_Bcast(&EIGEN[kloop][0], 2*n+1, MPI_DOUBLE, ID1, mpi_comm_level1);

  } /* kloop */
  
  /****************************************************
           band energy in a finite temperature
  ****************************************************/

  po = 0;
  loop_num = 0;
  ChemP_MAX = 10.0;  
  ChemP_MIN =-10.0;
    
  do {
    loop_num++;

    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);
    Num_State = 0.0;

    for (kloop=0; kloop<T_knum; kloop++){

      for (l=1; l<=2*n; l++){

	x = (EIGEN[kloop][l] - ChemP)*Beta;

	if      (x<=-x_cut)  FermiF = 1.0;
	else if (x>=x_cut)   FermiF = 0.0;
	else                 FermiF = 1.0/(1.0 + exp(x));

	Num_State += FermiF*(double)T_k_op[kloop];

      } /* l */
    } /* kloop */

    Num_State = Num_State/snum_i/snum_j/snum_k;

    Dnum = TZ - Num_State - system_charge;
    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<10e-13) po = 1;
  }
  while (po==0 && loop_num<2000);

  Eele0[0] = 0.0;
  Eele0[1] = 0.0;

  for (kloop=0; kloop<T_knum; kloop++){
    for (l=1; l<=2*n; l++){
      x = (EIGEN[kloop][l] - ChemP)*Beta;
      if      (x<=-x_cut)  FermiF = 1.0;
      else if (x_cut<=x)   FermiF = 0.0;
      else                 FermiF = 1.0/(1.0 + exp(x));
      Eele0[0] += FermiF*EIGEN[kloop][l]*(double)T_k_op[kloop];
    } /* l */
  } /* kloop */

  Eele0[0] = Eele0[0]/snum_i/snum_j/snum_k;
  Uele = Eele0[0];

  if (2<=level_stdout){
    printf("ChemP=%18.15f, Eele0[0]=%18.15f\n",ChemP,Eele0[0]);fflush(stdout);
  }

  dtime(&EiloopTime);
  if (myid==Host_ID){
    printf("<Band_DFT>  Eigen, time=%lf\n", EiloopTime-SiloopTime);
  }
  SiloopTime=EiloopTime;
  fflush(stdout);

  /****************************************************
    diagonalization for calculating density matrix

        CDM[0]  Re alpha alpha density matrix
        CDM[1]  Re beta  beta  density matrix
        CDM[2]  Re alpha beta  density matrix
        CDM[3]  Im alpha beta  density matrix
        iDM[0][0]  Im alpha alpha density matrix
        iDM[0][1]  Im beta  beta  density matrix

        EDM[0]  Re alpha alpha energy density matrix
        EDM[1]  Re beta  beta  energy density matrix
        EDM[2]  Re alpha beta  energy density matrix
        EDM[3]  Im alpha beta  energy density matrix
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];    
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];
        Bnum = MP[GB_AN];
        for (i=0; i<tnoA; i++){
          for (j=0; j<tnoB; j++){
            CDM[spin][MA_AN][LB_AN][i][j] = 0.0;
            EDM[spin][MA_AN][LB_AN][i][j] = 0.0;
	    iDM[0][0][MA_AN][LB_AN][i][j] = 0.0;
	    iDM[0][1][MA_AN][LB_AN][i][j] = 0.0;
          }
        }
      }
    }
  }

  /* for kloop */

  dtime(&SiloopTime);

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S and diagonalize it */

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

        k1 = T_KGrids1[kloop];
        k2 = T_KGrids2[kloop];
        k3 = T_KGrids3[kloop];

        Overlap_Band(ID,CntOLP,TmpM,MP,k1,k2,k3);
        n = TmpM[0][0].r;

        if (myid==ID){
	  for (i1=1; i1<=n; i1++){
	    for (j1=1; j1<=n; j1++){
	      S[i1][j1].r = TmpM[i1][j1].r;
	      S[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 

      }
    }

    kloop = arpo[myid];

    if (0<=kloop){

      EigenBand_lapack(S,ko,n,1);

      if (3<=level_stdout){
	printf(" myid=%2d kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       myid,kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
        for (i1=1; i1<=n; i1++){
          printf("  Eigenvalues of OLP (DM)  %2d  %15.12f\n",i1,ko[i1]);
        }
      }

      /* minus eigenvalues to 1.0e-14 */
      for (l=1; l<=n; l++){
        if (ko[l]<0.0){
          ko[l] = 1.0e-14;

          if (2<=level_stdout){
  	    printf("found an eigenvalue smaller than %10.8f of OLP kloop=%2d\n",
                           Threshold_OLP_Eigen,kloop);
	  }
	}

        koS[l] = ko[l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[l]);

      /* S * M1  */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*M1[j1];
	  S[i1][j1].i = S[i1][j1].i*M1[j1];
	} 
      } 

    } /* if (0<=kloop) */

    /* set H' and diagonalize it */

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

	k1 = T_KGrids1[kloop];
	k2 = T_KGrids2[kloop];
	k3 = T_KGrids3[kloop];

        /****************************************************
           make a full Hamiltonian matrix
        ****************************************************/

        /* non-spin-orbit coupling and non-LDA+U */  
        if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
            && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0)
          Hamiltonian_Band_NC(ID, nh, Dummy_ImNL, TmpM, MP, k1, k2, k3);

        /* spin-orbit coupling or LDA+U */  
        else
          Hamiltonian_Band_NC(ID, nh,       ImNL, TmpM, MP, k1, k2, k3);

	if (myid==ID){
	  for (i1=1; i1<=2*n; i1++){
	    for (j1=1; j1<=2*n; j1++){
	      H[i1][j1].r = TmpM[i1][j1].r;
	      H[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
	} 

      } /* if (0<=kloop) */
    } /* ID */

    kloop = arpo[myid];

    if (0<=kloop){

      /****************************************************
                     M1 * U^t * H * U * M1
      ****************************************************/

      /* transpose S */

      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      /* H * U * M1 */

      for (j1=1; j1<=n; j1++){
        for (i1=1; i1<=2*n; i1++){
	  for (m=0; m<=1; m++){

	    sum  = 0.0;
	    sumi = 0.0;
            mn = m*n;

	    for (l=1; l<=n; l++){
	      sum  += H[i1][l+mn].r*S[j1][l].r - H[i1][l+mn].i*S[j1][l].i;
	      sumi += H[i1][l+mn].r*S[j1][l].i + H[i1][l+mn].i*S[j1][l].r;
	    }

	    jj1 = 2*j1 - 1 + m;

	    C[jj1][i1].r = sum;
	    C[jj1][i1].i = sumi;
	  }
	}
      }     

      /* M1 * U^+ H * U * M1 */

      for (i1=1; i1<=n; i1++){

	for (m=0; m<=1; m++){

	  ii1 = 2*i1 - 1 + m;

	  for (j1=1; j1<=2*n; j1++){
	    sum  = 0.0;
	    sumi = 0.0;
            mn = m*n;
	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*C[j1][l+mn].r + S[i1][l].i*C[j1][l+mn].i;
	      sumi +=  S[i1][l].r*C[j1][l+mn].i - S[i1][l].i*C[j1][l+mn].r;
	    }
	    H[ii1][j1].r = sum;
	    H[ii1][j1].i = sumi;
	  }
	}
      }     

      /* H to C */

      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
	  C[i1][j1].r = H[i1][j1].r;
	  C[i1][j1].i = H[i1][j1].i;
	}
      }

      /* penalty for ill-conditioning states */

      EV_cut0 = Threshold_OLP_Eigen;

      for (i1=1; i1<=n; i1++){

	if (koS[i1]<EV_cut0){
	  C[2*i1-1][2*i1-1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  C[2*i1  ][2*i1  ].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	}

	/* cutoff the interaction between the ill-conditioned state */

	if (1.0e+3<C[2*i1-1][2*i1-1].r){
	  for (j1=1; j1<=2*n; j1++){
	    C[2*i1-1][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1-1] = Complex(0.0,0.0);
	    C[2*i1  ][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1  ] = Complex(0.0,0.0);
	  }
	  C[2*i1-1][2*i1-1] = Complex(1.0e+4,0.0);
	  C[2*i1  ][2*i1  ] = Complex(1.0e+4,0.0);
	}
      }

      /* solve eigenvalue problem */

      n1 = 2*n ;
      EigenBand_lapack(C, ko, n1, 1);

      if (3<=level_stdout && 0<=kloop){
        printf("  kloop %i, k1 k2 k3 %10.6f %10.6f %10.6f\n",
                kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
	for (i1=1; i1<=2*n; i1++){
	  printf("  Eigenvalues of Kohn-Sham(DM) %2d %15.12f\n",i1,ko[i1]);
	}
      }

      /****************************************************
          Transformation to the original eigenvectors
                 NOTE JRCAT-244p and JAIST-2122p 
                    C = U * lambda^{-1/2} * D
      ****************************************************/

      /*
      printf("C.r myid=%2d kloop=%2d\n",myid,kloop);
      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
          printf("%6.3f ",C[i1][j1].r); 
	}
        printf("\n");
      }

      printf("C.i myid=%2d kloop=%2d\n",myid,kloop);
      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
          printf("%6.3f ",C[i1][j1].i); 
	}
        printf("\n");
      }
      */

      /* transpose */
      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      /*
      printf("S.r myid=%2d kloop=%2d\n",myid,kloop);
      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
          printf("%6.3f ",S[i1][j1].r); 
	}
        printf("\n");
      }

      printf("S.i myid=%2d kloop=%2d\n",myid,kloop);
      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
          printf("%6.3f ",S[i1][j1].i); 
	}
        printf("\n");
      }
      */



      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
	  H2[i1][j1].r = 0.0;
	  H2[i1][j1].i = 0.0;
	}
      }

      for (m=0; m<=1; m++){
	for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n1; j1++){

	    sum  = 0.0; 
	    sumi = 0.0;

	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*C[2*(l-1)+1+m][j1].r
	             - S[i1][l].i*C[2*(l-1)+1+m][j1].i;
	      sumi +=  S[i1][l].r*C[2*(l-1)+1+m][j1].i
		     + S[i1][l].i*C[2*(l-1)+1+m][j1].r;
	    } 
	    H2[i1+m*n][j1].r = sum;
	    H2[i1+m*n][j1].i = sumi;
	  }
	}
      }

    } /* if (0<=kloop )*/

    /****************************************************
        calculate DM and EDM
    ****************************************************/

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

	/* MPI: H2 */

	k = 0;

	for (i1=0; i1<=2*n; i1++){
	  for (j1=0; j1<=2*n; j1++){
	    tmp_array0[k] = H2[i1][j1].r;
	    k++;
	  }
	}
         
	for (i1=0; i1<=2*n; i1++){
	  for (j1=0; j1<=2*n; j1++){
	    tmp_array0[k] = H2[i1][j1].i;
	    k++;
	  }
	}

	MPI_Barrier(mpi_comm_level1);
	MPI_Bcast(&tmp_array0[0], 2*(2*n+1)*(2*n+1), MPI_DOUBLE, ID, mpi_comm_level1);

	k = 0;

	for (i1=0; i1<=2*n; i1++){
	  for (j1=0; j1<=2*n; j1++){
	    H[i1][j1].r = tmp_array0[k];
	    k++;
	  }
	}
         
	for (i1=0; i1<=2*n; i1++){
	  for (j1=0; j1<=2*n; j1++){
	    H[i1][j1].i = tmp_array0[k];
	    k++;
	  }
	}

	/****************************************************
         calc of DM and EDM
	****************************************************/

	k1 = T_KGrids1[kloop];
	k2 = T_KGrids2[kloop];
	k3 = T_KGrids3[kloop];

	/* store Fermi function */

        po = 0;
	for (l=1; l<=2*n; l++){
	  TmpEigen = EIGEN[kloop][l];

	  x = (TmpEigen - ChemP)*Beta;

	  if      (x<-x_cut)  FermiF = 1.0;
	  else if (x>x_cut)   FermiF = 0.0;
	  else                FermiF = 1.0/(1.0 + exp(x));

	  VecFermiF[l] = FermiF; 

          if ( FermiF<=FermiEps && po==0 ) {
            lmax = l;  
            po = 1;
	  }
	}

	/* loop for MA_AN */

        Tmpk_op = (double)T_k_op[kloop];

	for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
	  GA_AN = M2G[MA_AN];
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

	    kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;
	    si = sin(2.0*PI*kRn);  /* Imaginary part */
	    co = cos(2.0*PI*kRn);  /* Real part      */
              
	    for (i1=0; i1<tnoA; i1++){
	      for (j1=0; j1<tnoB; j1++){

                d1  = 0.0; /* CDM[0]    */
                d2  = 0.0; /* EDM[0]    */
                d3  = 0.0; /* iDM[0][0] */
                d4  = 0.0; /* CDM[1]    */
                d5  = 0.0; /* EDM[1]    */
                d6  = 0.0; /* iDM[0][1] */
                d7  = 0.0; /* CDM[2]    */
                d8  = 0.0; /* EDM[2]    */
                d9  = 0.0; /* CDM[3]    */
                d10 = 0.0; /* EDM[3]    */

                ia = Anum + i1;
                ib = Anum + i1 + n;    
                ja = Bnum + j1;
                jb = Bnum + j1 + n;

		for (l=1; l<=lmax; l++){

	          FermiF = VecFermiF[l]; 
                  TmpEigen = EIGEN[kloop][l];

		  /* Re11 */

		  u2 = H[ia][l].r*H[ja][l].r;
		  v2 = H[ia][l].i*H[ja][l].i;
		  uv = H[ia][l].r*H[ja][l].i;
		  vu = H[ia][l].i*H[ja][l].r;

		  Redum = FermiF*( u2 + v2 );
		  Imdum = FermiF*( uv - vu );
		  Redum2 = (co*Redum - si*Imdum)*Tmpk_op;
		  ResumE = TmpEigen*Redum2;

		  d1 += Redum2;
		  d2 += ResumE;

		  /*
		  CDM[0][MA_AN][LB_AN][i1][j1] += Redum2;
		  EDM[0][MA_AN][LB_AN][i1][j1] += ResumE;
		  */

		  /* Im11 */

		  Imdum2 = (si*Redum + co*Imdum)*Tmpk_op;

		  d3 += Imdum2;

		  /*
		  iDM[0][0][MA_AN][LB_AN][i1][j1] += Imdum2;
		  */

		  /* Re22 */

		  u2 = H[ib][l].r*H[jb][l].r;
		  v2 = H[ib][l].i*H[jb][l].i;
		  uv = H[ib][l].r*H[jb][l].i;
		  vu = H[ib][l].i*H[jb][l].r;

		  Redum = FermiF*( u2 + v2 );
		  Imdum = FermiF*( uv - vu );
		  Redum2 = (co*Redum - si*Imdum)*Tmpk_op;
		  ResumE = TmpEigen*Redum2;

		  d4 += Redum2;
		  d5 += ResumE;

		  /*
		  CDM[1][MA_AN][LB_AN][i1][j1] += Redum2;
		  EDM[1][MA_AN][LB_AN][i1][j1] += ResumE;
		  */

		  /* Im22 */

		  Imdum2 = (si*Redum + co*Imdum)*Tmpk_op;

		  d6 += Imdum2;

		  /*
	          iDM[0][1][MA_AN][LB_AN][i1][j1] += Imdum2;
		  */

		  /* Re12 */

		  u2 = H[ia][l].r*H[jb][l].r;
		  v2 = H[ia][l].i*H[jb][l].i;
		  uv = H[ia][l].r*H[jb][l].i;
		  vu = H[ia][l].i*H[jb][l].r;

		  Redum = FermiF*( u2 + v2);
		  Imdum = FermiF*( uv - vu);
		  Redum2 = (co*Redum - si*Imdum)*Tmpk_op;
		  ResumE = TmpEigen*Redum2;

		  d7 += Redum2;
		  d8 += ResumE;

		  /*
		  CDM[2][MA_AN][LB_AN][i1][j1] += Redum2;
		  EDM[2][MA_AN][LB_AN][i1][j1] += ResumE;
		  */

		  /* Im12 */

		  u2 = H[ia][l].r*H[jb][l].r;
		  v2 = H[ia][l].i*H[jb][l].i;
		  uv = H[ia][l].r*H[jb][l].i;
		  vu = H[ia][l].i*H[jb][l].r;

		  Redum = FermiF*( u2 + v2);
		  Imdum = FermiF*( uv - vu);
		  Imdum2 = (co*Imdum + si*Redum)*Tmpk_op;
		  ImsumE = TmpEigen*Imdum2;

		  d9  += Imdum2;
		  d10 += ImsumE;

		  /*
		  CDM[3][MA_AN][LB_AN][i1][j1] += Imdum2;
		  EDM[3][MA_AN][LB_AN][i1][j1] += ImsumE;
		  */

		} /* l    */

                /* d1: CDM[0]    */
                /* d2: EDM[0]    */
                /* d3: iDM[0][0] */
                /* d4: CDM[1]    */
                /* d5: EDM[1]    */
                /* d6: iDM[0][1] */
                /* d7: CDM[2]    */
                /* d8: EDM[2]    */
                /* d9: CDM[3]    */
                /* d10:EDM[3]    */

		CDM[0][MA_AN][LB_AN][i1][j1]    += d1;
		EDM[0][MA_AN][LB_AN][i1][j1]    += d2;
 	        iDM[0][0][MA_AN][LB_AN][i1][j1] += d3;
	        CDM[1][MA_AN][LB_AN][i1][j1]    += d4;
	        EDM[1][MA_AN][LB_AN][i1][j1]    += d5;
	        iDM[0][1][MA_AN][LB_AN][i1][j1] += d6;
	        CDM[2][MA_AN][LB_AN][i1][j1]    += d7;
		EDM[2][MA_AN][LB_AN][i1][j1]    += d8;
		CDM[3][MA_AN][LB_AN][i1][j1]    += d9;
		EDM[3][MA_AN][LB_AN][i1][j1]    += d10;

	      }   /* j1     */
	    }     /* i1     */
	  }       /* LB_AN  */
	}         /* GA_AN  */
      }           /* if (0<=kloop) */
    }             /* ID     */

  }               /* kloop0 */


  dtime(&EiloopTime);
  if (myid==Host_ID){
    printf("<Band_DFT>  DM, time=%lf\n", EiloopTime-SiloopTime);
  }
  SiloopTime=EiloopTime;
  fflush(stdout);

  dum = 1.0/(snum_i*snum_j*snum_k);
  for (spin=0; spin<=SpinP_switch; spin++){
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];    
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];
        Bnum = MP[GB_AN];

        for (i=0; i<tnoA; i++){
          for (j=0; j<tnoB; j++){
            CDM[spin][MA_AN][LB_AN][i][j] = CDM[spin][MA_AN][LB_AN][i][j]*dum;
            EDM[spin][MA_AN][LB_AN][i][j] = EDM[spin][MA_AN][LB_AN][i][j]*dum;
          }
        }

        if (spin==0){
          for (i=0; i<tnoA; i++){
            for (j=0; j<tnoB; j++){
              iDM[0][0][MA_AN][LB_AN][i][j] = iDM[0][0][MA_AN][LB_AN][i][j]*dum;
              iDM[0][1][MA_AN][LB_AN][i][j] = iDM[0][1][MA_AN][LB_AN][i][j]*dum;
            }
          }
	}
      }
    }
  }

  /****************************************************
                     bond energies
  ****************************************************/

  My_Eele1[0] = 0.0;
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];    
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];

    for (j=0; j<=FNAN[GA_AN]; j++){
      GB_AN = natn[GA_AN][j];  
      wanB = WhatSpecies[GB_AN];
      tnoB = Spe_Total_CNO[wanB];

      /* non-spin-orbit coupling and non-LDA+U */  
      if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
          && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
        for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
            My_Eele1[0] = My_Eele1[0]
                         + CDM[0][MA_AN][j][k][l]*nh[0][MA_AN][j][k][l]
                         + CDM[1][MA_AN][j][k][l]*nh[1][MA_AN][j][k][l]
                         + 2.0*CDM[2][MA_AN][j][k][l]*nh[2][MA_AN][j][k][l]
                         - 2.0*CDM[3][MA_AN][j][k][l]*nh[3][MA_AN][j][k][l];
	  }
        }
      }

      /* spin-orbit coupling or LDA+U */  
      else {
        for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
            My_Eele1[0] = My_Eele1[0]
                         + CDM[0][MA_AN][j][k][l]*nh[0][MA_AN][j][k][l]
                         - iDM[0][0][MA_AN][j][k][l]*ImNL[0][MA_AN][j][k][l]
                         + CDM[1][MA_AN][j][k][l]*nh[1][MA_AN][j][k][l]
                         - iDM[0][1][MA_AN][j][k][l]*ImNL[1][MA_AN][j][k][l]
                         + 2.0*CDM[2][MA_AN][j][k][l]*nh[2][MA_AN][j][k][l]
                         - 2.0*CDM[3][MA_AN][j][k][l]*(nh[3][MA_AN][j][k][l]
                                                    +ImNL[2][MA_AN][j][k][l]);
	  }
        }
      }
  
    }
  }

  /* MPI, My_Eele1 */

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_Eele1[0], &Eele1[0], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  Eele1[1] = 0.0;

  /*
  printf("Eele00=%15.12f Eele01=%15.12f\n",Eele0[0],Eele0[1]);
  printf("Eele10=%15.12f Eele11=%15.12f\n",Eele1[0],Eele1[1]);
  MPI_Finalize();
  exit(0);
  */

  /****************************************************
                        Output
  ****************************************************/

  if (myid==Host_ID){
  
    strcpy(file_EV,".EV");
    fnjoint(filepath,filename,file_EV);

    if ((fp_EV = fopen(file_EV,"w")) != NULL){

#ifdef xt3
      setvbuf(fp_EV,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp_EV,"\n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"           Eigenvalues (Hartree) of SCF KS-eq.           \n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"***********************************************************\n\n");
      fprintf(fp_EV,"   Chemical Potential (Hatree) = %18.14f\n",ChemP);
      fprintf(fp_EV,"   Number of States            = %18.14f\n",Num_State);
      fprintf(fp_EV,"   Eigenvalues\n\n");

      m = 0;
      for (i=0; i<=(knum_i-1); i++){
        for (j=0; j<=(knum_j-1); j++){
	  for (k=0; k<=(knum_k-1); k++){

            if (k_op[i][j][k]>0){

  	      fprintf(fp_EV,"\n");
	      fprintf(fp_EV,"   i=%i j=%i k=%i\n",i,j,k);
	      fprintf(fp_EV,"   k1=%10.5f k2=%10.5f k3=%10.5f\n\n",
                      KGrids1[i],KGrids2[j],KGrids3[k]);
              for (l=1; l<=n*2; l++){
	        fprintf(fp_EV,"%5d  %18.14f\n",l,EIGEN[m][l]);
	      }

              m++;

	    }
	  }
        }
      }
      fclose(fp_EV);
    }
    else{
      printf("Failure of saving the EV file.\n");
      fclose(fp_EV);
    }  
  }

  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);
  free(ko);
  free(koS);
  free(arpo);

  free(tmp_array0);

  for (j=0; j<n2; j++){
    free(H[j]);
  }
  free(H);

  for (j=0; j<n2; j++){
    free(TmpM[j]);
  }
  free(TmpM);

  for (j=0; j<n2; j++){
    free(H2[j]);
  }
  free(H2);

  for (i=0; i<n2; i++){
    free(S[i]);
  }
  free(S);

  free(M1);
  free(VecFermiF);

  for (j=0; j<n2; j++){
    free(C[j]);
  }
  free(C);

  free(KGrids1); free(KGrids2); free(KGrids3);

  for (i=0;i<knum_i; i++) {
    for (j=0;j<knum_j; j++) {
      free(k_op[i][j]);
    }
    free(k_op[i]);
  }
  free(k_op);

  /* non-spin-orbit coupling and non-LDA+U */  
  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
      && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){

    free(Dummy_ImNL[0][0][0][0]);
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  free(T_KGrids1);
  free(T_KGrids2);
  free(T_KGrids3);
  free(T_k_op);

  for (j=0; j<T_knum; j++){
    free(EIGEN[j]);
  }
  free(EIGEN);

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}

