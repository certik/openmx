/**********************************************************************
  Band_DFT_Dosout.c:

     Band_DFT_Dosout.c is a subroutine to calculate the density of
     states based on the band calculation.

  Log of Band_DFT_Dosout.c:

     12/May/2003  Released by H.Kino
     25/Dec/2003  a non-collinear part (added by T.Ozaki)

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




static double Band_DFT_Dosout_Col(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP);

static double Band_DFT_DosoutGauss_Col(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP);

static double Band_DFT_Dosout_NonCol(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP);

static double Band_DFT_DosoutGauss_NonCol(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP);




double Band_DFT_Dosout( int knum_i, int knum_j, int knum_k,
                        int SpinP_switch,
                        double *****nh,
                        double *****ImNL,
                        double ****CntOLP)
{
  static double time0;

  if ( (SpinP_switch==0 || SpinP_switch==1) && Dos_fileout){
    time0 = Band_DFT_Dosout_Col( knum_i, knum_j, knum_k, SpinP_switch, nh, ImNL, CntOLP);
  }

  else if ( (SpinP_switch==0 || SpinP_switch==1) && DosGauss_fileout){
    time0 = Band_DFT_DosoutGauss_Col( knum_i, knum_j, knum_k, SpinP_switch, nh, ImNL, CntOLP);
  }

  else if (SpinP_switch==3 && Dos_fileout){
    time0 = Band_DFT_Dosout_NonCol( knum_i, knum_j, knum_k, SpinP_switch, nh, ImNL, CntOLP);
  }

  else if (SpinP_switch==3 && DosGauss_fileout){
    time0 = Band_DFT_DosoutGauss_NonCol( knum_i, knum_j, knum_k, SpinP_switch, nh, ImNL, CntOLP);
  }

  return time0;
}








static double Band_DFT_DosoutGauss_Col(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP)
{
  int i,j,k,spin,l,i1,j1,n1;
  int n, wanA;
  int *MP;
  int MA_AN, GA_AN, tnoA,Anum, LB_AN;
  int GB_AN, wanB, tnoB, Bnum, RnB;
  int l1,l2,l3,ik;
  int kloop,kloop0,S_knum,E_knum,T_knum,num_kloop0; 
  int ie,iemin,iemax,n1min,iecenter,iewidth;
  int MaxL,e1,s1;
  int i_vec[10];

  double sum,sumi,u2,v2,uv,vu,tmp,pi2,xa,eg,x,de;
  double factor,EV_cut0;
  double kRn,si,co,Redum,Imdum,Redum2,Resum;
  double TStime,TEtime,time0;
  double OLP_eigen_cut=Threshold_OLP_Eigen;

  double *koS;
  double **ko; int N_ko, i_ko[10];
  dcomplex **H;  int N_H,  i_H[10];
  dcomplex **S;  int N_S,  i_S[10];
  double *M1,***EIGEN;
  dcomplex **C;  int N_C,  i_C[10];
  double *KGrids1, *KGrids2, *KGrids3;
  double *SD;
  dcomplex ***H2,**TmpM,**TmpS;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  int *Ti_KGrids1,*Tj_KGrids2,*Tk_KGrids3,*arpo;
  dcomplex Ctmp1,Ctmp2;

  double ****Dummy_ImNL;
  double snum_i, snum_j, snum_k; 
  double k1,k2,k3;
  double *r_energy;
  double time1,time2,time3;
  double time4,time5,time6;
  double time7,time8,time9,time10;
  double Stime1,Etime1;
  float ***Dos;
  float *DosVec;
  char file_ev[YOUSO10],file_eig[YOUSO10];
  FILE *fp_ev,*fp_eig;
  char buf1[fp_bsize];          /* setvbuf */
  char buf2[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID,ID1,ID2,tag;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){
    printf("\n<Band_DFT_Dosout>\n"); fflush(stdout);
  }

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

  /****************************************************
             calculation of the array size
  ****************************************************/

  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n  += Spe_Total_CNO[wanA];
  }

  /****************************************************
   Allocation of arrays
  ****************************************************/

  DosVec = (float*)malloc(sizeof(float)*(atomnum+1)*List_YOUSO[7]);

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  arpo = (int*)malloc(sizeof(int)*numprocs);

  N_ko=2; i_ko[0]=List_YOUSO[23]; i_ko[1]=n+1;
  ko=(double**)malloc_multidimarray("double",N_ko, i_ko); 

  koS = (double*)malloc(sizeof(double)*(n+1));

  N_H=2; i_H[0]=i_H[1]=n+1;
  H=(dcomplex**)malloc_multidimarray("dcomplex",N_H, i_H);

  N_S=2; i_S[0]=i_S[1]=n+1;
  S=(dcomplex**)malloc_multidimarray("dcomplex",N_S, i_S);

  M1 = (double*)malloc(sizeof(double)*(n+1));

  N_C=2; i_C[0]=i_C[1]=n+1;
  C=(dcomplex**)malloc_multidimarray("dcomplex",N_C, i_C);

  KGrids1 = (double*)malloc(sizeof(double)*knum_i);
  KGrids2 = (double*)malloc(sizeof(double)*knum_j);
  KGrids3 = (double*)malloc(sizeof(double)*knum_k);

  SD = (double*)malloc(sizeof(double)*(atomnum+1)*List_YOUSO[7]*2);

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  TmpS = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    TmpS[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  H2 = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    H2[i] = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (j=0; j<n+1; j++){
      H2[i][j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }
  }

  /* allocate Dos */

  Dos = (float***)malloc(sizeof(float**)*DosGauss_Num_Mesh);
  for (ik=0; ik<DosGauss_Num_Mesh; ik++) {
    Dos[ik] = (float**)malloc(sizeof(float*)*(SpinP_switch+1) );
    for (spin=0; spin<=SpinP_switch; spin++) {
      Dos[ik][spin] = (float*)malloc(sizeof(float)*(atomnum+1)*List_YOUSO[7]);
      for (i=0; i<(atomnum+1)*List_YOUSO[7]; i++)  Dos[ik][spin][i] = 0.0;
    }
  }

  /* allocate r_energy */

  r_energy = (double*)malloc(sizeof(double)*DosGauss_Num_Mesh);

  /* set up energies where DOS is calculated */

  Dos_Erange[0] += ChemP;
  Dos_Erange[1] += ChemP;

  de = (Dos_Erange[1]-Dos_Erange[0])/(double)DosGauss_Num_Mesh;
  for (i=0; i<DosGauss_Num_Mesh; i++) {
    r_energy[i] = Dos_Erange[0] + de*(double)i;
  }

  /* no spin-orbit coupling */
  if (SO_switch==0){
    Dummy_ImNL = (double****)malloc(sizeof(double***)*1);
    Dummy_ImNL[0] = (double***)malloc(sizeof(double**)*1);
    Dummy_ImNL[0][0] = (double**)malloc(sizeof(double*)*1);
    Dummy_ImNL[0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  snum_i = knum_i;
  snum_j = knum_j;
  snum_k = knum_k;

  /* set up  k-grids */
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
    for (i=0;i<=knum_i-1;i++) printf("%f ",KGrids1[i]);
    printf("\n");
    printf(" KGrids2: ");
    for (i=0;i<=knum_j-1;i++) printf("%f ",KGrids2[i]);
    printf("\n");
    printf(" KGrids3: ");
    for (i=0;i<=knum_k-1;i++) printf("%f ",KGrids3[i]);
    printf("\n");
  }

  /***********************************
       one-dimentionalize for MPI
  ************************************/

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        T_knum++;
      }
    }
  }

  T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);

  Ti_KGrids1 = (int*)malloc(sizeof(int)*T_knum);
  Tj_KGrids2 = (int*)malloc(sizeof(int)*T_knum);
  Tk_KGrids3 = (int*)malloc(sizeof(int)*T_knum);

  EIGEN  = (double***)malloc(sizeof(double**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    EIGEN[i] = (double**)malloc(sizeof(double*)*T_knum);
    for (j=0; j<T_knum; j++){
      EIGEN[i][j] = (double*)malloc(sizeof(double)*(n+1));
    }
  }

  /* set T_KGrid1,2,3 */

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){

	T_KGrids1[T_knum] = KGrids1[i];
	T_KGrids2[T_knum] = KGrids2[j];
	T_KGrids3[T_knum] = KGrids3[k];

	Ti_KGrids1[T_knum] = i;
	Tj_KGrids2[T_knum] = j;
	Tk_KGrids3[T_knum] = k;

        T_knum++;
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

  factor = 1.0/(double)(knum_i * knum_j * knum_k); /* normalization factor */
  pi2 = sqrt(PI);

  /* for kloop */

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S */

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (measure_time==1) dtime(&Stime1);
      
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

	      TmpS[i1][j1].r = TmpM[i1][j1].r;
	      TmpS[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 
       
      }

      if (measure_time==1){
        dtime(&Etime1);
        time1 += Etime1 - Stime1;
      }
    }

    kloop = arpo[myid];

    if (0<=kloop){

      if (measure_time==1) dtime(&Stime1);

      EigenBand_lapack(S,ko[0],n, 1);

      if (measure_time==1){
        dtime(&Etime1);
        time2 += Etime1 - Stime1;
      }

      if (3<=level_stdout && 0<=kloop){
	printf("  kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
	for (i1=1; i1<=n; i1++){
	  printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[0][i1]);
	}
      }

      /* minus eigenvalues to 1.0e-14 */

      for (l=1; l<=n; l++){
        if (ko[0][l]<0.0) ko[0][l] = 1.0e-14;
        koS[l] = ko[0][l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[0][l]);

      /* S * M1  */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*M1[j1];
	  S[i1][j1].i = S[i1][j1].i*M1[j1];
	} 
      } 

    } /* if (0<=kloop) */

    /* set H */

    for (spin=0; spin<=SpinP_switch; spin++){

      for (ID=0; ID<numprocs; ID++){

        kloop = arpo[ID];

        if (0<=kloop){
          k1 = T_KGrids1[kloop];
          k2 = T_KGrids2[kloop];
          k3 = T_KGrids3[kloop];

          Hamiltonian_Band(ID, nh[spin], TmpM, MP, k1, k2, k3);

          if (myid==ID){
	    for (i1=1; i1<=n; i1++){
	      for (j1=1; j1<=n; j1++){
	        H[i1][j1].r = TmpM[i1][j1].r;
	        H[i1][j1].i = TmpM[i1][j1].i;
	      } 
	    } 
          } 
	}
      }

      kloop = arpo[myid];

      if (0<=kloop){

        /****************************************************
 	                 M1 * U^t * H * U * M1
        ****************************************************/

        if (measure_time==1) dtime(&Stime1);

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

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum,sumi,l)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();
 
	  for (i1=1+OMPID; i1<=n; i1+=Nthrds){

	    for (j1=1; j1<=n; j1++){

	      sum  = 0.0;
	      sumi = 0.0;

	      for (l=1; l<=n; l++){
		sum  += H[i1][l].r*S[j1][l].r - H[i1][l].i*S[j1][l].i;
		sumi += H[i1][l].r*S[j1][l].i + H[i1][l].i*S[j1][l].r;
	      }

	      C[j1][i1].r = sum;
	      C[j1][i1].i = sumi;
	    }
	  }     
	} /* #pragma omp parallel */
 
        /* M1 * U^+ H * U * M1 */

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum,sumi,l)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (i1=1+OMPID; i1<=n; i1+=Nthrds){
	    for (j1=1; j1<=n; j1++){

	      sum  = 0.0;
	      sumi = 0.0;

	      for (l=1; l<=n; l++){
		sum  +=  S[i1][l].r*C[j1][l].r + S[i1][l].i*C[j1][l].i;
		sumi +=  S[i1][l].r*C[j1][l].i - S[i1][l].i*C[j1][l].r;
	      }
	      H[i1][j1].r = sum;
	      H[i1][j1].i = sumi;
	    }
	  }     
	} /* #pragma omp parallel */

        /* H to C */

	for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n; j1++){
	    C[i1][j1] = H[i1][j1];
	  }
	}

       /* penalty for ill-conditioning states */

	EV_cut0 = Threshold_OLP_Eigen;

	for (i1=1; i1<=n; i1++){

	  if (koS[i1]<EV_cut0){
	    C[i1][i1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  }
 
	  /* cutoff the interaction between the ill-conditioned state */
 
	  if (1.0e+3<C[i1][i1].r){
	    for (j1=1; j1<=n; j1++){
	      C[i1][j1] = Complex(0.0,0.0);
	      C[j1][i1] = Complex(0.0,0.0);
	    }
	    C[i1][i1].r = 1.0e+4;
	  }
	}

        if (measure_time==1){
          dtime(&Etime1);
          time3 += Etime1 - Stime1;
        }

        /* solve eigenvalue problem */

        if (measure_time==1) dtime(&Stime1);

	n1 = n;
        EigenBand_lapack(C,ko[spin],n1,1);

	for (i1=1; i1<=n; i1++) EIGEN[spin][kloop][i1] = ko[spin][i1];

        if (measure_time==1){
          dtime(&Etime1);
          time4 += Etime1 - Stime1;
        }

        /****************************************************
	     transformation to the original eigen vectors.
	     NOTE JRCAT-244p and JAIST-2122p 
	     C = U * lambda^{-1/2} * D
        ****************************************************/

        if (measure_time==1) dtime(&Stime1);

	/* transpose */

	for (i1=1; i1<=n; i1++){
	  for (j1=i1+1; j1<=n; j1++){
	    Ctmp1 = S[i1][j1];
	    Ctmp2 = S[j1][i1];
	    S[i1][j1] = Ctmp2;
	    S[j1][i1] = Ctmp1;
	  }
	}

	/* transpose */

	for (i1=1; i1<=n; i1++){
	  for (j1=i1+1; j1<=n; j1++){
	    Ctmp1 = C[i1][j1];
	    Ctmp2 = C[j1][i1];
	    C[i1][j1] = Ctmp2;
	    C[j1][i1] = Ctmp1;
	  }
	}

        /* shift */

	for (j1=1; j1<=n; j1++){
	  for (l=n; 1<=l; l--){
   	    C[j1][l] = C[j1][l];
	  }
	}

#pragma omp parallel shared(C,S,H2,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum,sumi,l)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (i1=1+OMPID; i1<=n; i1+=Nthrds){
	    for (j1=1; j1<=n1; j1++){

	      sum  = 0.0;
	      sumi = 0.0;

	      for (l=1; l<=n; l++){
		sum  +=  S[i1][l].r*C[j1][l].r - S[i1][l].i*C[j1][l].i;
		sumi +=  S[i1][l].r*C[j1][l].i + S[i1][l].i*C[j1][l].r;
	      }

	      H2[spin][j1][i1].r = sum;
	      H2[spin][j1][i1].i = sumi;

	    }
	  }
	} /* #pragma omp parallel */

        if (measure_time==1){
          dtime(&Etime1);
          time5 += Etime1 - Stime1;
        }

      } /* if (0<=kloop) */
    } /* spin */

    /****************************************************
        store LDOS
    ****************************************************/

    kloop = arpo[myid];

    if (0<=kloop){

      for (spin=0; spin<=SpinP_switch; spin++){

	k1 = T_KGrids1[kloop];
	k2 = T_KGrids2[kloop];
	k3 = T_KGrids3[kloop];

	i = Ti_KGrids1[kloop];
	j = Tj_KGrids2[kloop];
	k = Tk_KGrids3[kloop];

	for (l=1; l<=n; l++){

	  /* initialize */

	  for (i1=0; i1<(atomnum+1)*List_YOUSO[7]; i1++){
	    SD[i1] = 0.0;
	  }

	  /* calculate SD */

#pragma omp parallel shared(TmpS,SD,spin,l,H2,MP,Spe_Total_CNO,WhatSpecies,atomnum) private(OMPID,Nthrds,Nprocs,GA_AN,wanA,tnoA,Anum,GB_AN,wanB,tnoB,Bnum,i1,j1,u2,v2,uv,vu,Redum,Imdum)
	  { 

	    /* get info. on OpenMP */ 

	    OMPID = omp_get_thread_num();
	    Nthrds = omp_get_num_threads();
	    Nprocs = omp_get_num_procs();

	    for (GA_AN=1+OMPID; GA_AN<=atomnum; GA_AN+=Nthrds){

	      wanA = WhatSpecies[GA_AN];
	      tnoA = Spe_Total_CNO[wanA];
	      Anum = MP[GA_AN];

	      for (GB_AN=1; GB_AN<=atomnum; GB_AN++){

		wanB = WhatSpecies[GB_AN];
		tnoB = Spe_Total_CNO[wanB];
		Bnum = MP[GB_AN];

		for (i1=0; i1<tnoA; i1++){
		  for (j1=0; j1<tnoB; j1++){

		    u2 = H2[spin][l][Anum+i1].r*H2[spin][l][Bnum+j1].r;
		    v2 = H2[spin][l][Anum+i1].i*H2[spin][l][Bnum+j1].i;
		    uv = H2[spin][l][Anum+i1].r*H2[spin][l][Bnum+j1].i;
		    vu = H2[spin][l][Anum+i1].i*H2[spin][l][Bnum+j1].r;

		    Redum = (u2 + v2);
		    Imdum = (uv - vu);

		    SD[Anum+i1] += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  } /* j1 */
		} /* i1 */
	      } /* GB_AN */
	    } /* GA_AN */

	  } /* #pragma omp parallel */

          /*  calculate contribution to Gaussian DOS */

#pragma omp parallel shared(pi2,SD,factor,Dos,r_energy,DosGauss_Width,DosGauss_Num_Mesh,Dos_Erange,l,kloop,spin,EIGEN,atomnum) private(OMPID,Nthrds,Nprocs,GA_AN,wanA,tnoA,Anum,eg,x,iecenter,iewidth,iemin,iemax,ie,xa,i1)
	  { 

	    /* get info. on OpenMP */ 

	    OMPID = omp_get_thread_num();
	    Nthrds = omp_get_num_threads();
	    Nprocs = omp_get_num_procs();

	    for (GA_AN=1+OMPID; GA_AN<=atomnum; GA_AN+=Nthrds){

	      wanA = WhatSpecies[GA_AN];
	      tnoA = Spe_Total_CNO[wanA];
	      Anum = MP[GA_AN];

	      /*************************
                  exp( -(x/a)^2 ) 
                  a= DosGauss_Width
	          x=-3a : 3a is counted 
	      **************************/

	      eg = EIGEN[spin][kloop][l];
	      x = (eg-Dos_Erange[0])/(Dos_Erange[1]-Dos_Erange[0])*(DosGauss_Num_Mesh-1);
	      iecenter = (int)x ; 

	      iewidth = DosGauss_Width*3.0/(Dos_Erange[1]-Dos_Erange[0])*(DosGauss_Num_Mesh-1)+3;

	      iemin = iecenter - iewidth;
	      iemax = iecenter + iewidth;

	      if (iemin<0) iemin=0;
	      if (iemax>=DosGauss_Num_Mesh) iemax=DosGauss_Num_Mesh-1;

	      if ( 0<=iemin && iemin<DosGauss_Num_Mesh && 0<=iemax && iemax<DosGauss_Num_Mesh ) {

		for (ie=iemin;ie<=iemax;ie++) {
		  xa = (eg-r_energy[ie])/DosGauss_Width;

		  for (i1=0; i1<tnoA; i1++){
		    Dos[ie][spin][Anum+i1] += (float)(factor*SD[Anum+i1]* exp(-xa*xa)/(DosGauss_Width*pi2));
		  }
		}
	      }

	    }   /* GA_AN */

	  } /* #pragma omp parallel */

	} /* l */
      } /* spin */   
    } /* if (0<=kloop) */

  } /* kloop0 */

  /****************************************************************
                    write eigenvalues and eigenvectors
  ****************************************************************/

  if (measure_time==1) dtime(&Stime1);

  if (myid==Host_ID){

    sprintf(file_eig,"%s%s.Dos.val",filepath,filename);

    if ( (fp_eig=fopen(file_eig,"w"))==NULL ) {
      printf("can not open a file %s\n",file_eig);
    }
  
    else{

      /* write *.Dos.val */

      printf("write eigenvalues\n");

      fprintf(fp_eig,"mode        7\n");
      fprintf(fp_eig,"NonCol      0\n");
      fprintf(fp_eig,"N           %d\n",n);
      fprintf(fp_eig,"Nspin       %d\n",SpinP_switch);
      fprintf(fp_eig,"Erange      %lf %lf\n",Dos_Erange[0],Dos_Erange[1]);
      fprintf(fp_eig,"Kgrid       %d %d %d\n",1,1,1);
      fprintf(fp_eig,"atomnum     %d\n",atomnum);
      fprintf(fp_eig,"<WhatSpecies\n");
      for (i=1; i<=atomnum; i++) {
        fprintf(fp_eig,"%d ",WhatSpecies[i]);
      }
      fprintf(fp_eig,"\nWhatSpecies>\n");
      fprintf(fp_eig,"SpeciesNum     %d\n",SpeciesNum);
      fprintf(fp_eig,"<Spe_Total_CNO\n");
      for (i=0;i<SpeciesNum;i++) {
        fprintf(fp_eig,"%d ",Spe_Total_CNO[i]);
      }
      fprintf(fp_eig,"\nSpe_Total_CNO>\n");

      MaxL = Supported_MaxL; 

      fprintf(fp_eig,"MaxL           %d\n",Supported_MaxL);
      fprintf(fp_eig,"<Spe_Num_CBasis\n");
      for (i=0;i<SpeciesNum;i++) {
        for (l=0;l<=MaxL;l++) {
	  fprintf(fp_eig,"%d ",Spe_Num_CBasis[i][l]);
        }
        fprintf(fp_eig,"\n");
      }
      fprintf(fp_eig,"Spe_Num_CBasis>\n");
      fprintf(fp_eig,"ChemP       %lf\n",ChemP);

      fprintf(fp_eig,"irange      %d %d\n",0,DosGauss_Num_Mesh-1);
      fprintf(fp_eig,"<Eigenvalues\n");
      for (spin=0; spin<=SpinP_switch; spin++) {
        fprintf(fp_eig,"%d %d %d ",0,0,0);

        for (ik=0; ik<DosGauss_Num_Mesh; ik++) {
          fprintf(fp_eig,"%lf ",r_energy[ik]);
	}
        fprintf(fp_eig,"\n");
      }  
      fprintf(fp_eig,"Eigenvalues>\n");

      fclose(fp_eig);
    }

  } /* if (myid==Host_ID) */

   /* write *.Dos.vec */
  
  if (myid==Host_ID){
    printf("write eigenvectors\n");
    sprintf(file_ev,"%s%s.Dos.vec",filepath,filename);
    if ( (fp_ev=fopen(file_ev,"w"))==NULL ) {
      printf("cannot open a file %s\n",file_ev);
    }
  }

  for (spin=0; spin<=SpinP_switch; spin++) {
    for (ik=0; ik<DosGauss_Num_Mesh; ik++) {

      i_vec[0] = 0; 
      i_vec[1] = 0;
      i_vec[2] = 0;

      MPI_Reduce(&Dos[ik][spin][1], &DosVec[1], n, MPI_FLOAT, MPI_SUM, Host_ID, mpi_comm_level1);

      if (myid==Host_ID){
        fwrite(i_vec,sizeof(int),3,fp_ev);
        fwrite(&DosVec[1],sizeof(float),n,fp_ev);
      }
    }
  }

  if (myid==Host_ID){
    fclose(fp_ev);
  }

  if (measure_time==1){
    dtime(&Etime1);
    time10 += Etime1 - Stime1;
  }

  if (measure_time==1){
    printf("myid=%4d  time1=%15.12f\n",myid,time1);
    printf("myid=%4d  time2=%15.12f\n",myid,time2);
    printf("myid=%4d  time3=%15.12f\n",myid,time3);
    printf("myid=%4d  time4=%15.12f\n",myid,time4);
    printf("myid=%4d  time5=%15.12f\n",myid,time5);
    printf("myid=%4d  time6=%15.12f\n",myid,time6);
    printf("myid=%4d  time7=%15.12f\n",myid,time7);
    printf("myid=%4d  time8=%15.12f\n",myid,time8);
    printf("myid=%4d  time9=%15.12f\n",myid,time9);
    printf("myid=%4d time10=%15.12f\n",myid,time10);
  }

  /****************************************************
                       free arrays
  ****************************************************/

  free(DosVec);

  free(MP);
  free(arpo);

  for (i=0; i<i_ko[0]; i++){
    free(ko[i]);
  }
  free(ko);

  free(koS);

  for (i=0; i<i_H[0]; i++){
    free(H[i]);
  }
  free(H);

  for (i=0; i<i_S[0]; i++){
    free(S[i]);
  }
  free(S);

  free(M1);

  for (i=0; i<i_C[0]; i++){
    free(C[i]);
  }
  free(C);

  free(KGrids1); free(KGrids2);free(KGrids3);
  
  free(SD);

  for (j=0; j<n+1; j++){
    free(TmpM[j]);
  }
  free(TmpM);

  for (j=0; j<n+1; j++){
    free(TmpS[j]);
  }
  free(TmpS);

  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<n+1; j++){
      free(H2[i][j]);
    }
    free(H2[i]);
  }
  free(H2);

  /* no spin-orbit coupling */
  if (SO_switch==0){
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  free(T_KGrids1);
  free(T_KGrids2);
  free(T_KGrids3);

  free(Ti_KGrids1);
  free(Tj_KGrids2);
  free(Tk_KGrids3);

  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<T_knum; j++){
      free(EIGEN[i][j]);
    }
    free(EIGEN[i]);
  }
  free(EIGEN);

  for (ik=0; ik<DosGauss_Num_Mesh; ik++) {
    for (spin=0; spin<=SpinP_switch; spin++) {
      free(Dos[ik][spin]);
    }
    free(Dos[ik]);
  }
  free(Dos);

  free(r_energy);

  /* for elapsed time */

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;

}







static double Band_DFT_DosoutGauss_NonCol(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP)
{
  int i,j,k,spin,l,i1,j1,ik,n1;
  int n, wanA,ii1,jj1,m,n2;
  int *MP;
  int MA_AN, GA_AN, tnoA,Anum, LB_AN;
  int GB_AN, wanB, tnoB, Bnum, RnB;
  int l1,l2,l3;
 
  int ie,iemin,iemax,iemin0,iemax0,n1min,mn;
  int MaxL,e1,s1,T_knum,S_knum,E_knum;
  int kloop,kloop0,num_kloop0;
  int i_vec[10];

  double EV_cut0;
  double sum,sumi,u2,v2,uv,vu,tmp;
  double sum_r0,sum_i0,sum_r1,sum_i1;
  double d0,d1,d2,d3;
  double kRn,si,co,Redum,Imdum;
  double Redum2,Resum,Imdum2;
  double TStime,TEtime,time0;
  double OLP_eigen_cut=Threshold_OLP_Eigen;
  double theta,phi,sit,cot,sip,cop,de;
  double tmp1,tmp2,tmp3,eg,x,xa,pi2,factor;
  double *ko; int N_ko, i_ko[10];
  double *koS;
  double *r_energy;
  dcomplex **H;  int N_H,  i_H[10];
  dcomplex **S;  int N_S,  i_S[10];
  dcomplex Ctmp1,Ctmp2;
  double *M1,**EIGEN;
  dcomplex **C;  int N_C,  i_C[10];
  double *KGrids1, *KGrids2, *KGrids3;
  double *SD;
  dcomplex **TmpM,**TmpS;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  int *Ti_KGrids1,*Tj_KGrids2,*Tk_KGrids3,*arpo;
  int iecenter,iewidth;

  double *****Dummy_ImNL;
  float **Dos;
  float *DosVec;
  double snum_i, snum_j, snum_k; 
  double k1,k2,k3;

  char file_ev[YOUSO10],file_eig[YOUSO10];
  FILE *fp_ev,*fp_eig;
  char buf1[fp_bsize];          /* setvbuf */
  char buf2[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID,ID1,ID2,tag;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){
    printf("\n<Band_DFT_Dosout>\n"); fflush(stdout);
  }

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
   Allocation of arrays
  ****************************************************/

  DosVec = (float*)malloc(sizeof(float)*2*(atomnum+1)*List_YOUSO[7]);

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);

  arpo = (int*)malloc(sizeof(int)*numprocs);

  ko = (double*)malloc(sizeof(double)*n2);
  koS = (double*)malloc(sizeof(double)*(n+1));

  H = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    H[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  S = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (i=0; i<n2; i++){
    S[i] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  M1 = (double*)malloc(sizeof(double)*n2);

  C = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    C[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  KGrids1 = (double*)malloc(sizeof(double)*knum_i);
  KGrids2 = (double*)malloc(sizeof(double)*knum_j);
  KGrids3 = (double*)malloc(sizeof(double)*knum_k);

  SD = (double*)malloc(sizeof(double)*(atomnum+1)*List_YOUSO[7]*2);

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  TmpS = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    TmpS[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  /* allocate Dos */

  Dos = (float**)malloc(sizeof(float*)*DosGauss_Num_Mesh);
  for (ik=0; ik<DosGauss_Num_Mesh; ik++) {
    Dos[ik] = (float*)malloc(sizeof(float)*(atomnum+1)*List_YOUSO[7]*2);
    for (i=0; i<(atomnum+1)*List_YOUSO[7]*2; i++)  Dos[ik][i] = 0.0;
  }

  /* allocate r_energy */

  r_energy = (double*)malloc(sizeof(double)*DosGauss_Num_Mesh);

  /* set up energies where DOS is calculated */

  Dos_Erange[0] += ChemP;
  Dos_Erange[1] += ChemP;

  de = (Dos_Erange[1]-Dos_Erange[0])/(double)DosGauss_Num_Mesh;
  for (i=0; i<DosGauss_Num_Mesh; i++) {
    r_energy[i] = Dos_Erange[0] + de*(double)i;
  }

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
    set up k-grids
  ****************************************************/

  snum_i = knum_i;
  snum_j = knum_j;
  snum_k = knum_k;

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
    for (i=0;i<=knum_i-1;i++) printf("%f ",KGrids1[i]);
    printf("\n");
    printf(" KGrids2: ");
    for (i=0;i<=knum_j-1;i++) printf("%f ",KGrids2[i]);
    printf("\n");
    printf(" KGrids3: ");
    for (i=0;i<=knum_k-1;i++) printf("%f ",KGrids3[i]);
    printf("\n");
  }

  /***********************************
       one-dimentionalize for MPI
  ************************************/

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        T_knum++;
      }
    }
  }

  T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);

  Ti_KGrids1 = (int*)malloc(sizeof(int)*T_knum);
  Tj_KGrids2 = (int*)malloc(sizeof(int)*T_knum);
  Tk_KGrids3 = (int*)malloc(sizeof(int)*T_knum);

  EIGEN = (double**)malloc(sizeof(double*)*T_knum);
  for (j=0; j<T_knum; j++){
    EIGEN[j] = (double*)malloc(sizeof(double)*n2);
  }

  /* set T_KGrid1,2,3 */

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){

	T_KGrids1[T_knum] = KGrids1[i];
	T_KGrids2[T_knum] = KGrids2[j];
	T_KGrids3[T_knum] = KGrids3[k];

	Ti_KGrids1[T_knum] = i;
	Tj_KGrids2[T_knum] = j;
	Tk_KGrids3[T_knum] = k;

        T_knum++;
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

  factor = 1.0/(double)(knum_i * knum_j * knum_k); /* normalization factor */
  pi2 = sqrt(PI);

  /* for kloop */

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S */

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

	      TmpS[i1][j1].r = TmpM[i1][j1].r;
	      TmpS[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 

      }
    }

    kloop = arpo[myid];

    if (0<=kloop){

      EigenBand_lapack(S,ko,n,1);

      if (2<=level_stdout && 0<=kloop){
	printf("  kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
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

    }

    /****************************************************
       make a full Hamiltonian matrix
    ****************************************************/

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){
        k1 = T_KGrids1[kloop];
        k2 = T_KGrids2[kloop];
        k3 = T_KGrids3[kloop];

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
      }
    }

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

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,jj1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (i1=1+OMPID; i1<=2*n; i1+=Nthrds){

	  jj1 = 1;
	  for (j1=1; j1<=n; j1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    for (l=1; l<=n; l++){

	      sum_r0 += H[i1][l  ].r*S[j1][l].r - H[i1][l  ].i*S[j1][l].i;
	      sum_i0 += H[i1][l  ].r*S[j1][l].i + H[i1][l  ].i*S[j1][l].r;

	      sum_r1 += H[i1][l+n].r*S[j1][l].r - H[i1][l+n].i*S[j1][l].i;
	      sum_i1 += H[i1][l+n].r*S[j1][l].i + H[i1][l+n].i*S[j1][l].r;
	    }

	    C[jj1][i1].r = sum_r0;
	    C[jj1][i1].i = sum_i0;  jj1++;

	    C[jj1][i1].r = sum_r1;
	    C[jj1][i1].i = sum_i1;  jj1++;

	  }
	}     

      } /* #pragma omp parallel */

      /* M1 * U^+ H * U * M1 */

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,ii1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=2*n; j1+=Nthrds){

	  ii1 = 1;
	  for (i1=1; i1<=n; i1++){

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

	    H[ii1][j1].r = sum_r0;
	    H[ii1][j1].i = sum_i0; ii1++;

	    H[ii1][j1].r = sum_r1;
	    H[ii1][j1].i = sum_i1; ii1++;

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

      n1 = 2*n;
      EigenBand_lapack(C, ko, n1, 1);

      for (i1=1; i1<=2*n; i1++) EIGEN[kloop][i1] = ko[i1];

      /****************************************************
	   Transformation to the original eigenvectors.
	   NOTE JRCAT-244p and JAIST-2122p 
	   C = U * lambda^{-1/2} * D
      ****************************************************/

      /* transpose */
      for (i1=1; i1<=2*n; i1++){
	for (j1=i1+1; j1<=2*n; j1++){
	  Ctmp1 = C[i1][j1];
	  Ctmp2 = C[j1][i1];
	  C[i1][j1] = Ctmp2;
	  C[j1][i1] = Ctmp1;
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

#pragma omp parallel shared(C,S,H,n,n1) private(OMPID,Nthrds,Nprocs,i1,j1,l1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=n1; j1+=Nthrds){
	  for (i1=1; i1<=n; i1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    l1 = 1;
	    for (l=1; l<=n; l++){

	      sum_r0 += S[i1][l].r*C[j1][l1].r
		       -S[i1][l].i*C[j1][l1].i; 
	      sum_i0 += S[i1][l].r*C[j1][l1].i
		       +S[i1][l].i*C[j1][l1].r; l1++;

	      sum_r1 += S[i1][l].r*C[j1][l1].r
		       -S[i1][l].i*C[j1][l1].i;
	      sum_i1 += S[i1][l].r*C[j1][l1].i
	 	       +S[i1][l].i*C[j1][l1].r; l1++;
	    } 

	    H[j1][i1  ].r = sum_r0;
	    H[j1][i1  ].i = sum_i0;

	    H[j1][i1+n].r = sum_r1;
	    H[j1][i1+n].i = sum_i1;
	  }
	}

      } /* #pragma omp parallel */

    } /* if (0<=kloop) */

    /****************************************************
        store LDOS
    ****************************************************/

    kloop = arpo[myid];

    if (0<=kloop){

      k1 = T_KGrids1[kloop];
      k2 = T_KGrids2[kloop];
      k3 = T_KGrids3[kloop];

      i = Ti_KGrids1[kloop];
      j = Tj_KGrids2[kloop];
      k = Tk_KGrids3[kloop];

      for (l=1; l<=2*n; l++){

	/* calculate SD */

#pragma omp parallel shared(TmpS,SD,H,l,MP,Spe_Total_CNO,WhatSpecies,atomnum) private(OMPID,Nthrds,Nprocs,GA_AN,wanA,tnoA,Anum,theta,phi,sit,cot,sip,cop,d0,d1,d2,d3,GB_AN,wanB,tnoB,Bnum,i1,j1,u2,v2,uv,vu,Redum,Imdum,tmp1,tmp2,tmp3)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (GA_AN=1+OMPID; GA_AN<=atomnum; GA_AN+=Nthrds){

	    wanA = WhatSpecies[GA_AN];
	    tnoA = Spe_Total_CNO[wanA];
	    Anum = MP[GA_AN];
	    theta = Angle0_Spin[GA_AN];
	    phi   = Angle1_Spin[GA_AN];
	    sit = sin(theta);
	    cot = cos(theta);
	    sip = sin(phi);
	    cop = cos(phi);

	    for (i1=0; i1<tnoA; i1++){

	      d0 = 0.0;
	      d1 = 0.0;
	      d2 = 0.0;
	      d3 = 0.0;

	      for (GB_AN=1; GB_AN<=atomnum; GB_AN++){

		wanB = WhatSpecies[GB_AN];
		tnoB = Spe_Total_CNO[wanB];
		Bnum = MP[GB_AN];

		for (j1=0; j1<tnoB; j1++){

		  /* Re11 */
		  u2 = H[l][Anum+i1].r*H[l][Bnum+j1].r;
		  v2 = H[l][Anum+i1].i*H[l][Bnum+j1].i;
		  uv = H[l][Anum+i1].r*H[l][Bnum+j1].i;
		  vu = H[l][Anum+i1].i*H[l][Bnum+j1].r;
		  Redum = (u2 + v2);
		  Imdum = (uv - vu);
		  d0 += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  /* Re22 */
		  u2 = H[l][Anum+i1+n].r*H[l][Bnum+j1+n].r;
		  v2 = H[l][Anum+i1+n].i*H[l][Bnum+j1+n].i;
		  uv = H[l][Anum+i1+n].r*H[l][Bnum+j1+n].i;
		  vu = H[l][Anum+i1+n].i*H[l][Bnum+j1+n].r;
		  Redum = (u2 + v2);
		  Imdum = (uv - vu);
		  d1 += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  /* Re12 */
		  u2 = H[l][Anum+i1].r*H[l][Bnum+j1+n].r;
		  v2 = H[l][Anum+i1].i*H[l][Bnum+j1+n].i;
		  uv = H[l][Anum+i1].r*H[l][Bnum+j1+n].i;
		  vu = H[l][Anum+i1].i*H[l][Bnum+j1+n].r;
		  Redum = (u2 + v2);
		  Imdum = (uv - vu);
		  d2 += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  /* Im12
		     conjugate complex of Im12 due to difference in the definition
		     between density matrix and charge density
		  */
 
		  d3 += -(Imdum*TmpS[Anum+i1][Bnum+j1].r + Redum*TmpS[Anum+i1][Bnum+j1].i);
                    
		} /* j1 */
	      } /* GB_AN */

	      tmp1 = 0.5*(d0 + d1);
	      tmp2 = 0.5*cot*(d0 - d1);
	      tmp3 = (d2*cop - d3*sip)*sit;

	      SD[2*(Anum-1)+i1]      = tmp1 + tmp2 + tmp3;
	      SD[2*(Anum-1)+tnoA+i1] = tmp1 - tmp2 - tmp3;

	    } /* i1 */
	  } /* GA_AN */

	} /* #pragma omp parallel */

          /*  calculate contribution to Gaussian DOS */

	for (GA_AN=1; GA_AN<=atomnum; GA_AN++){

	  wanA = WhatSpecies[GA_AN];
	  tnoA = Spe_Total_CNO[wanA];
	  Anum = MP[GA_AN];

	  /*************************
                  exp( -(x/a)^2 ) 
                  a= DosGauss_Width
	          x=-3a : 3a is counted 
	  **************************/

	  eg = EIGEN[kloop][l];
	  x = (eg-Dos_Erange[0])/(Dos_Erange[1]-Dos_Erange[0])*(DosGauss_Num_Mesh-1);
	  iecenter = (int)x ; 

	  iewidth = DosGauss_Width*3.0/(Dos_Erange[1]-Dos_Erange[0])*(DosGauss_Num_Mesh-1)+3;

	  iemin = iecenter - iewidth;
	  iemax = iecenter + iewidth;

	  if (iemin<0) iemin=0;
	  if (iemax>=DosGauss_Num_Mesh) iemax=DosGauss_Num_Mesh-1;

	  if ( 0<=iemin && iemin<DosGauss_Num_Mesh && 0<=iemax && iemax<DosGauss_Num_Mesh ) {

	    for (ie=iemin;ie<=iemax;ie++) {
	      xa = (eg-r_energy[ie])/DosGauss_Width;

	      for (i1=0; i1<tnoA; i1++){
		Dos[ie][2*(Anum-1)+i1     ] += (float)(factor*SD[2*(Anum-1)+i1     ]*exp(-xa*xa)/(DosGauss_Width*pi2));
		Dos[ie][2*(Anum-1)+tnoA+i1] += (float)(factor*SD[2*(Anum-1)+tnoA+i1]*exp(-xa*xa)/(DosGauss_Width*pi2));
	      }
	    }
	  }

	} /* GA_AN */
      } /* l */
    } /* if (0<=kloop) */
  } /* kloop0 */

  /****************************************************************
                    write eigenvalues and eigenvectors
  ****************************************************************/

  if (myid==Host_ID){

    sprintf(file_eig,"%s%s.Dos.val",filepath,filename);

    if ( (fp_eig=fopen(file_eig,"w"))==NULL ) {
      printf("can not open a file %s\n",file_eig);
    }
  
    else{

      /* write *.Dos.val */

      printf("write eigenvalues\n");

      fprintf(fp_eig,"mode        7\n");
      fprintf(fp_eig,"NonCol      1\n");
      fprintf(fp_eig,"N           %d\n",n);
      fprintf(fp_eig,"Nspin       %d\n",1);  /* switch to 1 */
      fprintf(fp_eig,"Erange      %lf %lf\n",Dos_Erange[0],Dos_Erange[1]);
      fprintf(fp_eig,"Kgrid       %d %d %d\n",1,1,1);
      fprintf(fp_eig,"atomnum     %d\n",atomnum);
      fprintf(fp_eig,"<WhatSpecies\n");
      for (i=1; i<=atomnum; i++) {
        fprintf(fp_eig,"%d ",WhatSpecies[i]);
      }
      fprintf(fp_eig,"\nWhatSpecies>\n");
      fprintf(fp_eig,"SpeciesNum     %d\n",SpeciesNum);
      fprintf(fp_eig,"<Spe_Total_CNO\n");
      for (i=0;i<SpeciesNum;i++) {
        fprintf(fp_eig,"%d ",Spe_Total_CNO[i]);
      }
      fprintf(fp_eig,"\nSpe_Total_CNO>\n");

      MaxL = Supported_MaxL; 

      fprintf(fp_eig,"MaxL           %d\n",Supported_MaxL);
      fprintf(fp_eig,"<Spe_Num_CBasis\n");
      for (i=0;i<SpeciesNum;i++) {
        for (l=0;l<=MaxL;l++) {
	  fprintf(fp_eig,"%d ",Spe_Num_CBasis[i][l]);
        }
        fprintf(fp_eig,"\n");
      }
      fprintf(fp_eig,"Spe_Num_CBasis>\n");
      fprintf(fp_eig,"ChemP       %lf\n",ChemP);

      fprintf(fp_eig,"<SpinAngle\n");
      for (i=1; i<=atomnum; i++) {
        fprintf(fp_eig,"%lf %lf\n",Angle0_Spin[i],Angle1_Spin[i]);
      }
      fprintf(fp_eig,"SpinAngle>\n");

      fprintf(fp_eig,"irange      %d %d\n",0,DosGauss_Num_Mesh-1);
      fprintf(fp_eig,"<Eigenvalues\n");
      for (spin=0; spin<=1; spin++) {
        fprintf(fp_eig,"%d %d %d ",0,0,0);

        for (ik=0; ik<DosGauss_Num_Mesh; ik++) {
          fprintf(fp_eig,"%lf ",r_energy[ik]);
	}
        fprintf(fp_eig,"\n");
      }  
      fprintf(fp_eig,"Eigenvalues>\n");

      fclose(fp_eig);
    }

  } /* if (myid==Host_ID) */

  if (myid==Host_ID){

    /* write *.Dos.vec */
    printf("write eigenvectors\n");
    sprintf(file_ev,"%s%s.Dos.vec",filepath,filename);
    if ( (fp_ev=fopen(file_ev,"w"))==NULL ) {
      printf("can not open a file %s\n",file_ev);
    }
  }

  for (ik=0; ik<DosGauss_Num_Mesh; ik++) {

    i_vec[0] = 0; 
    i_vec[1] = 0;
    i_vec[2] = 0;

    MPI_Reduce(&Dos[ik][0], &DosVec[0], 2*n, MPI_FLOAT, MPI_SUM, Host_ID, mpi_comm_level1);

    if (myid==Host_ID){
      fwrite(i_vec,sizeof(int),3,fp_ev);
      fwrite(&DosVec[0],sizeof(float),2*n,fp_ev);
    }

  }

  if (myid==Host_ID){
    fclose(fp_ev);
  }

  /****************************************************
                       free arrays
  ****************************************************/

  free(DosVec);

  free(MP);

  free(arpo);

  free(ko);
  free(koS);

  for (j=0; j<n2; j++){
    free(H[j]);
  }
  free(H);

  for (i=0; i<n2; i++){
    free(S[i]);
  }
  free(S);

  free(M1);

  for (j=0; j<n2; j++){
    free(C[j]);
  }
  free(C);

  free(KGrids1); free(KGrids2);free(KGrids3);

  free(SD);

  for (j=0; j<n2; j++){
    free(TmpM[j]);
  }
  free(TmpM);

  for (j=0; j<n+1; j++){
    free(TmpS[j]);
  }
  free(TmpS);

  for (ik=0; ik<DosGauss_Num_Mesh; ik++) {
    free(Dos[ik]);
  }
  free(Dos);

  free(r_energy);

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

  free(Ti_KGrids1);
  free(Tj_KGrids2);
  free(Tk_KGrids3);

  for (j=0; j<T_knum; j++){
    free(EIGEN[j]);
  }
  free(EIGEN);

  /* for elapsed time */

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}







static double Band_DFT_Dosout_Col(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP)
{
  int i,j,k,spin,l,i1,j1,n1;
  int n, wanA;
  int *MP;
  int MA_AN, GA_AN, tnoA,Anum, LB_AN;
  int GB_AN, wanB, tnoB, Bnum, RnB;
  int l1,l2,l3;
  int kloop,kloop0,S_knum,E_knum,T_knum,num_kloop0; 
  int ie,iemin,iemax,iemin0,iemax0,n1min;
  int MaxL,e1,s1;
  int i_vec[10];

  double sum,sumi,u2,v2,uv,vu,tmp;
  double kRn,si,co,Redum,Imdum,Redum2,Resum;
  double TStime,TEtime,time0;
  double OLP_eigen_cut=Threshold_OLP_Eigen;

  double EV_cut0;
  double *koS;
  double **ko; int N_ko, i_ko[10];
  dcomplex **H;  int N_H,  i_H[10];
  dcomplex **S;  int N_S,  i_S[10];
  double *M1,***EIGEN;
  dcomplex **C;  int N_C,  i_C[10];
  double *KGrids1, *KGrids2, *KGrids3;
  float *SD; int N_SD, i_SD[10];
  dcomplex ***H2,**TmpM,**TmpS;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  int *Ti_KGrids1,*Tj_KGrids2,*Tk_KGrids3,*arpo;
  int *num_allocated_k;
  dcomplex Ctmp1,Ctmp2;

  double ****Dummy_ImNL;
  double snum_i, snum_j, snum_k; 
  double k1,k2,k3;
  char file_ev[YOUSO10],file_eig[YOUSO10];
  char file_ev0[YOUSO10];
  FILE *fp_ev0,*fp_ev,*fp_eig;
  char buf1[fp_bsize];          /* setvbuf */
  char buf2[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID,ID1,ID2,tag;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){
    printf("\n<Band_DFT_Dosout>\n"); fflush(stdout);
  }

  dtime(&TStime);

  /****************************************************
             calculation of the array size
  ****************************************************/

  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n  += Spe_Total_CNO[wanA];
  }

  /****************************************************
   Allocation of arrays
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  arpo = (int*)malloc(sizeof(int)*numprocs);

  num_allocated_k = (int*)malloc(sizeof(int)*numprocs);

  N_ko=2; i_ko[0]=List_YOUSO[23]; i_ko[1]=n+1;
  ko=(double**)malloc_multidimarray("double",N_ko, i_ko); 

  koS = (double*)malloc(sizeof(double)*(n+1));

  N_H=2; i_H[0]=i_H[1]=n+1;
  H=(dcomplex**)malloc_multidimarray("dcomplex",N_H, i_H);

  N_S=2; i_S[0]=i_S[1]=n+1;
  S=(dcomplex**)malloc_multidimarray("dcomplex",N_S, i_S);

  M1 = (double*)malloc(sizeof(double)*(n+1));

  N_C=2; i_C[0]=i_C[1]=n+1;
  C=(dcomplex**)malloc_multidimarray("dcomplex",N_C, i_C);

  KGrids1 = (double*)malloc(sizeof(double)*knum_i);
  KGrids2 = (double*)malloc(sizeof(double)*knum_j);
  KGrids3 = (double*)malloc(sizeof(double)*knum_k);

  SD=(float*)malloc(sizeof(float)*(atomnum+1)*List_YOUSO[7]);

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  TmpS = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    TmpS[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  H2 = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    H2[i] = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (j=0; j<n+1; j++){
      H2[i][j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }
  }

  /* no spin-orbit coupling */
  if (SO_switch==0){
    Dummy_ImNL = (double****)malloc(sizeof(double***)*1);
    Dummy_ImNL[0] = (double***)malloc(sizeof(double**)*1);
    Dummy_ImNL[0][0] = (double**)malloc(sizeof(double*)*1);
    Dummy_ImNL[0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  snum_i = knum_i;
  snum_j = knum_j;
  snum_k = knum_k;

  /* set up  k-grids */
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
    for (i=0;i<=knum_i-1;i++) printf("%f ",KGrids1[i]);
    printf("\n");
    printf(" KGrids2: ");
    for (i=0;i<=knum_j-1;i++) printf("%f ",KGrids2[i]);
    printf("\n");
    printf(" KGrids3: ");
    for (i=0;i<=knum_k-1;i++) printf("%f ",KGrids3[i]);
    printf("\n");
  }

  /***********************************
       one-dimentionalize for MPI
  ************************************/

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        T_knum++;
      }
    }
  }

  T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);

  Ti_KGrids1 = (int*)malloc(sizeof(int)*T_knum);
  Tj_KGrids2 = (int*)malloc(sizeof(int)*T_knum);
  Tk_KGrids3 = (int*)malloc(sizeof(int)*T_knum);

  EIGEN  = (double***)malloc(sizeof(double**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    EIGEN[i] = (double**)malloc(sizeof(double*)*T_knum);
    for (j=0; j<T_knum; j++){
      EIGEN[i][j] = (double*)malloc(sizeof(double)*(n+1));
    }
  }

  /* set T_KGrid1,2,3 */

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){

	T_KGrids1[T_knum] = KGrids1[i];
	T_KGrids2[T_knum] = KGrids2[j];
	T_KGrids3[T_knum] = KGrids3[k];

	Ti_KGrids1[T_knum] = i;
	Tj_KGrids2[T_knum] = j;
	Tk_KGrids3[T_knum] = k;

        T_knum++;
      }
    }
  }

  /* allocate k-points into proccessors */

  if (T_knum<=myid){
    S_knum = -10;
    E_knum = -100;
    num_kloop0 = 1;
    num_allocated_k[myid] = 0;
  }
  else if (T_knum<numprocs) {
    S_knum = myid;
    E_knum = myid;
    num_kloop0 = 1;
    num_allocated_k[myid] = 1;
  }
  else {
    tmp = (double)T_knum/(double)numprocs; 
    num_kloop0 = (int)tmp + 1;

    S_knum = (int)((double)myid*(tmp+0.0001)); 
    E_knum = (int)((double)(myid+1)*(tmp+0.0001)) - 1;
    if (myid==(numprocs-1)) E_knum = T_knum - 1;
    if (E_knum<0)           E_knum = 0;
    num_allocated_k[myid] = E_knum - S_knum + 1;
  }

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&num_allocated_k[ID], 1, MPI_INT, ID, mpi_comm_level1);
  }

  /****************************************************************
                      find iemin and iemax
  *****************************************************************/

  iemin=n; iemax=1; n1min=1;

  k1 = 0.0;
  k2 = 0.0;
  k3 = 0.0;

  Overlap_Band(Host_ID,CntOLP,S,MP,k1,k2,k3);

  if (myid==Host_ID){

    n = S[0][0].r;
    EigenBand_lapack(S, ko[0], n, 1);

    if (3<=level_stdout){
      printf("  k1 k2 k3 %10.6f %10.6f %10.6f\n",k1,k2,k3);
      for (i1=1; i1<=n; i1++){
	printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[0][i1]);
      }
    }

    /* minus eigenvalues to 1.0e-14 */

    for (l=1; l<=n; l++){
      if (ko[0][l]<0.0) ko[0][l] = 1.0e-14;
      koS[l] = ko[0][l];
    }

    /* calculate S*1/sqrt(ko) */

    for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[0][l]);

    /* S * M1  */

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
	S[i1][j1].r = S[i1][j1].r*M1[j1];
	S[i1][j1].i = S[i1][j1].i*M1[j1];
      } 
    } 

  } /* if (myid==Host_ID) */

  spin = 0;

  /****************************************************
               make a full Hamiltonian matrix
  ****************************************************/

  Hamiltonian_Band(Host_ID, nh[spin], H, MP, k1, k2, k3);

  if (myid==Host_ID){

    /****************************************************
                    M1 * U^t * H * U * M1
    ****************************************************/

    /* transpose S */

    if (spin==0){
      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }
    }

    /* H * U * M1 */
 
    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
	sum = 0.0; sumi=0.0;
	for (l=1; l<=n; l++){
	  sum  += H[i1][l].r*S[j1][l].r - H[i1][l].i*S[j1][l].i;
	  sumi += H[i1][l].r*S[j1][l].i + H[i1][l].i*S[j1][l].r;
	}

	C[j1][i1].r = sum;
	C[j1][i1].i = sumi;
      }
    }     

    /* M1 * U^+ H * U * M1 */

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
	sum = 0.0;
	sumi=0.0;
	for (l=1; l<=n; l++){
	  sum  +=  S[i1][l].r*C[j1][l].r + S[i1][l].i*C[j1][l].i;
	  sumi +=  S[i1][l].r*C[j1][l].i - S[i1][l].i*C[j1][l].r;
	}
	H[i1][j1].r = sum;
	H[i1][j1].i = sumi;
      }
    }     

    /* H to C */

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
	C[i1][j1] = H[i1][j1];
      }
    }

    /* penalty for ill-conditioning states */

    EV_cut0 = Threshold_OLP_Eigen;

    for (i1=1; i1<=n; i1++){

      if (koS[i1]<EV_cut0){
	C[i1][i1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
      }
 
      /* cutoff the interaction between the ill-conditioned state */
 
      if (1.0e+3<C[i1][i1].r){
	for (j1=1; j1<=n; j1++){
	  C[i1][j1] = Complex(0.0,0.0);
	  C[j1][i1] = Complex(0.0,0.0);
	}
	C[i1][i1].r = 1.0e+4;
      }
    }

    /* solve eigenvalue problem */

    n1 = n;
    EigenBand_lapack(C,ko[spin],n1, 1);

    if (n1min<n1) n1min=n1;

    iemin0=1;
    for (i1=1;i1<=n1;i1++) {
      if (ko[spin][i1]>(ChemP+Dos_Erange[0]) ) {
	iemin0=i1-1;
	break;
      }
    }
    if (iemin0<1) iemin0=1;

    iemax0=n1;
    for (i1=iemin0;i1<=n1;i1++) {
      if (ko[spin][i1]>(ChemP+Dos_Erange[1]) ) {
	iemax0=i1;
	break;
      }
    }
    if (iemax0>n1) iemax0=n1;

    if (iemin>iemin0) iemin=iemin0;
    if (iemax<iemax0) iemax=iemax0;

  } /* if (myid==Host_ID) */

  /* add a buffer to iemin and iemax */

  iemin -= 5;
  iemax += 5;

  if (iemin<1) iemin = 1;
  if (n<iemax) iemax = n;

  if (myid==Host_ID){
    if (n1min<iemax) iemax=n1min;
    printf(" iemin and iemax= %d %d\n",iemin,iemax);
  }

  /* MPI, iemin, iemax */
  MPI_Barrier(mpi_comm_level1);
  MPI_Bcast(&iemin, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Bcast(&iemax, 1, MPI_INT, Host_ID, mpi_comm_level1);

  /****************************************************************
                     eigenvalues and eigenvectors
   ****************************************************************/

  if (myid==Host_ID){

    sprintf(file_eig,"%s%s.Dos.val",filepath,filename);
    if ( (fp_eig=fopen(file_eig,"w"))==NULL ) {

#ifdef xt3
      setvbuf(fp_eig,buf1,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      printf("can not open a file %s\n",file_eig);
    }

    if ( fp_eig==NULL ) {
      goto Finishing;
    }
  }

  sprintf(file_ev,"%s%s.Dos.vec%d",filepath,filename,myid);
  if ( (fp_ev=fopen(file_ev,"w"))==NULL ) {
    printf("cannot open a file %s\n",file_ev);
  }
  if ( fp_ev==NULL ) {
    goto Finishing;
  }

  if (myid==Host_ID){

    if (fp_eig) {
      fprintf(fp_eig,"mode        1\n");
      fprintf(fp_eig,"NonCol      0\n");
      fprintf(fp_eig,"N           %d\n",n);
      fprintf(fp_eig,"Nspin       %d\n",SpinP_switch);
      fprintf(fp_eig,"Erange      %lf %lf\n",Dos_Erange[0],Dos_Erange[1]);
      fprintf(fp_eig,"irange      %d %d\n",iemin,iemax);
      fprintf(fp_eig,"Kgrid       %d %d %d\n",knum_i,knum_j,knum_k);
      fprintf(fp_eig,"atomnum     %d\n",atomnum);
      fprintf(fp_eig,"<WhatSpecies\n");
      for (i=1;i<=atomnum;i++) {
        fprintf(fp_eig,"%d ",WhatSpecies[i]);
      }
      fprintf(fp_eig,"\nWhatSpecies>\n");
      fprintf(fp_eig,"SpeciesNum     %d\n",SpeciesNum);
      fprintf(fp_eig,"<Spe_Total_CNO\n");
      for (i=0;i<SpeciesNum;i++) {
        fprintf(fp_eig,"%d ",Spe_Total_CNO[i]);
      }
      fprintf(fp_eig,"\nSpe_Total_CNO>\n");

      MaxL=Supported_MaxL; 

      fprintf(fp_eig,"MaxL           %d\n",Supported_MaxL);
      fprintf(fp_eig,"<Spe_Num_CBasis\n");
      for (i=0;i<SpeciesNum;i++) {
        for (l=0;l<=MaxL;l++) {
          fprintf(fp_eig,"%d ",Spe_Num_CBasis[i][l]);
        }
        fprintf(fp_eig,"\n");
      }
      fprintf(fp_eig,"Spe_Num_CBasis>\n");

      fprintf(fp_eig,"ChemP       %lf\n",ChemP);

      fprintf(fp_eig,"<Eigenvalues\n");

    }
    if (fp_eig) {
      printf(" write eigenvalues\n");
    }
    if (fp_ev) {
      printf(" write eigenvectors\n");
    }
  }

  /* for kloop */

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S */

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

	      TmpS[i1][j1].r = TmpM[i1][j1].r;
	      TmpS[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 
       
      }
    }

    kloop = arpo[myid];

    if (0<=kloop){

      EigenBand_lapack(S,ko[0],n, 1);

      if (3<=level_stdout && 0<=kloop){
	printf("  kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
	for (i1=1; i1<=n; i1++){
	  printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[0][i1]);
	}
      }

      /* minus eigenvalues to 1.0e-14 */

      for (l=1; l<=n; l++){
	if (ko[0][l]<0.0) ko[0][l] = 1.0e-14;
	koS[l] = ko[0][l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[0][l]);

      /* S * M1  */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*M1[j1];
	  S[i1][j1].i = S[i1][j1].i*M1[j1];
	} 
      } 

    } /* if (0<=kloop) */

    /* set H */

    for (spin=0; spin<=SpinP_switch; spin++){

      for (ID=0; ID<numprocs; ID++){

        kloop = arpo[ID];

        if (0<=kloop){
          k1 = T_KGrids1[kloop];
          k2 = T_KGrids2[kloop];
          k3 = T_KGrids3[kloop];

          Hamiltonian_Band(ID, nh[spin], TmpM, MP, k1, k2, k3);

          if (myid==ID){
	    for (i1=1; i1<=n; i1++){
	      for (j1=1; j1<=n; j1++){
	        H[i1][j1].r = TmpM[i1][j1].r;
	        H[i1][j1].i = TmpM[i1][j1].i;
	      } 
	    } 
          } 
	}
      }

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

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum,sumi,l)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();
 
	  for (i1=1+OMPID; i1<=n; i1+=Nthrds){

	    for (j1=1; j1<=n; j1++){

	      sum  = 0.0;
	      sumi = 0.0;

	      for (l=1; l<=n; l++){
		sum  += H[i1][l].r*S[j1][l].r - H[i1][l].i*S[j1][l].i;
		sumi += H[i1][l].r*S[j1][l].i + H[i1][l].i*S[j1][l].r;
	      }

	      C[j1][i1].r = sum;
	      C[j1][i1].i = sumi;
	    }
	  }     
	} /* #pragma omp parallel */

        /* M1 * U^+ H * U * M1 */

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum,sumi,l)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (i1=1+OMPID; i1<=n; i1+=Nthrds){
	    for (j1=1; j1<=n; j1++){

	      sum  = 0.0;
	      sumi = 0.0;

	      for (l=1; l<=n; l++){
		sum  +=  S[i1][l].r*C[j1][l].r + S[i1][l].i*C[j1][l].i;
		sumi +=  S[i1][l].r*C[j1][l].i - S[i1][l].i*C[j1][l].r;
	      }
	      H[i1][j1].r = sum;
	      H[i1][j1].i = sumi;
	    }
	  }     
	} /* #pragma omp parallel */

        /* H to C */

	for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n; j1++){
	    C[i1][j1] = H[i1][j1];
	  }
	}

	/* penalty for ill-conditioning states */

	EV_cut0 = Threshold_OLP_Eigen;

	for (i1=1; i1<=n; i1++){

	  if (koS[i1]<EV_cut0){
	    C[i1][i1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  }
 
	  /* cutoff the interaction between the ill-conditioned state */
 
	  if (1.0e+3<C[i1][i1].r){
	    for (j1=1; j1<=n; j1++){
	      C[i1][j1] = Complex(0.0,0.0);
	      C[j1][i1] = Complex(0.0,0.0);
	    }
	    C[i1][i1].r = 1.0e+4;
	  }
	}

        /* solve eigenvalue problem */

	n1 = n;
        EigenBand_lapack(C,ko[spin],n1, 1);

	for (i1=1; i1<=n; i1++) EIGEN[spin][kloop][i1] = ko[spin][i1];

        /****************************************************
	     transformation to the original eigen vectors.
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

	/* transpose */

	for (i1=1; i1<=n; i1++){
	  for (j1=i1+1; j1<=n; j1++){
	    Ctmp1 = C[i1][j1];
	    Ctmp2 = C[j1][i1];
	    C[i1][j1] = Ctmp2;
	    C[j1][i1] = Ctmp1;
	  }
	}

        /* shift */

	for (j1=1; j1<=n; j1++){
	  for (l=n; 1<=l; l--){
   	    C[j1][l] = C[j1][l];
	  }
	}

#pragma omp parallel shared(C,S,H2,n) private(OMPID,Nthrds,Nprocs,i1,j1,sum,sumi,l)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (i1=1+OMPID; i1<=n; i1+=Nthrds){
	    for (j1=1; j1<=n1; j1++){

	      sum  = 0.0;
	      sumi = 0.0;

	      for (l=1; l<=n; l++){
		sum  +=  S[i1][l].r*C[j1][l].r - S[i1][l].i*C[j1][l].i;
		sumi +=  S[i1][l].r*C[j1][l].i + S[i1][l].i*C[j1][l].r;
	      }

	      H2[spin][j1][i1].r = sum;
	      H2[spin][j1][i1].i = sumi;

	    }
	  }
	} /* #pragma omp parallel */

      } /* if (0<=kloop) */
    } /* spin */

    /****************************************************
        store LDOS
    ****************************************************/

    kloop = arpo[myid];

    if (0<=kloop){

      for (spin=0; spin<=SpinP_switch; spin++){

	k1 = T_KGrids1[kloop];
	k2 = T_KGrids2[kloop];
	k3 = T_KGrids3[kloop];

	i = Ti_KGrids1[kloop];
	j = Tj_KGrids2[kloop];
	k = Tk_KGrids3[kloop];

	for (l=iemin; l<=iemax; l++){

	  /* initialize */

	  for (i1=0; i1<(atomnum+1)*List_YOUSO[7]; i1++){
	    SD[i1] = 0.0;
	  }

	  /* calculate SD */

#pragma omp parallel shared(TmpS,SD,spin,l,H2,MP,Spe_Total_CNO,WhatSpecies,atomnum) private(OMPID,Nthrds,Nprocs,GA_AN,wanA,tnoA,Anum,GB_AN,wanB,tnoB,Bnum,i1,j1,u2,v2,uv,vu,Redum,Imdum)
	  { 

	    /* get info. on OpenMP */ 

	    OMPID = omp_get_thread_num();
	    Nthrds = omp_get_num_threads();
	    Nprocs = omp_get_num_procs();

	    for (GA_AN=1+OMPID; GA_AN<=atomnum; GA_AN+=Nthrds){

	      wanA = WhatSpecies[GA_AN];
	      tnoA = Spe_Total_CNO[wanA];
	      Anum = MP[GA_AN];

	      for (GB_AN=1; GB_AN<=atomnum; GB_AN++){

		wanB = WhatSpecies[GB_AN];
		tnoB = Spe_Total_CNO[wanB];
		Bnum = MP[GB_AN];

		for (i1=0; i1<tnoA; i1++){
		  for (j1=0; j1<tnoB; j1++){

		    u2 = H2[spin][l][Anum+i1].r*H2[spin][l][Bnum+j1].r;
		    v2 = H2[spin][l][Anum+i1].i*H2[spin][l][Bnum+j1].i;
		    uv = H2[spin][l][Anum+i1].r*H2[spin][l][Bnum+j1].i;
		    vu = H2[spin][l][Anum+i1].i*H2[spin][l][Bnum+j1].r;

		    Redum = (u2 + v2);
		    Imdum = (uv - vu);

		    SD[Anum+i1] += (float)(Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i);

		  } /* j1 */
		} /* i1 */
	      } /* GB_AN */
	    } /* GA_AN */

	  } /* #pragma omp parallel */

	  /*********************************************
                   writing a binary file 
	  *********************************************/

	  i_vec[0] = Ti_KGrids1[kloop];
	  i_vec[1] = Tj_KGrids2[kloop];
	  i_vec[2] = Tk_KGrids3[kloop];

	  fwrite(i_vec,sizeof(int),3,fp_ev);
	  fwrite(&SD[1],sizeof(float),n,fp_ev);

	} /* l */
      } /* spin */
    } /* if (kloop<=0) */ 

  } /* kloop0        */

  /****************************************************
     MPI:

     EIGEN
  ****************************************************/

  tmp = (double)T_knum/(double)numprocs; 

  for (spin=0; spin<=SpinP_switch; spin++){
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

      MPI_Bcast(&EIGEN[spin][kloop][0], n+1, MPI_DOUBLE, ID1, mpi_comm_level1);
    }
  }

  if (myid==Host_ID){
    if (fp_eig) {

      for (kloop=0; kloop<T_knum; kloop++){
	for (spin=0; spin<=SpinP_switch; spin++){

          i = Ti_KGrids1[kloop];
          j = Tj_KGrids2[kloop];
          k = Tk_KGrids3[kloop];

	  fprintf(fp_eig,"%d %d %d ",i,j,k);
	  for (ie=iemin;ie<=iemax;ie++) {
	    fprintf(fp_eig,"%lf ",EIGEN[spin][kloop][ie]);
	  }
	  fprintf(fp_eig,"\n");

	}
      }

      fprintf(fp_eig,"Eigenvalues>\n");
    }
  }

Finishing: 

  if (myid==Host_ID){
    if (fp_eig) { 
      fclose(fp_eig);
    }
  }

  if (fp_ev) {
    fclose(fp_ev);
  }

  /****************************************************
     merge *.Dos.vec#
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){

    sprintf(file_ev0,"%s%s.Dos.vec",filepath,filename);
    if ( (fp_ev0=fopen(file_ev0,"w"))==NULL ) {
      printf("cannot open a file %s\n",file_ev0);
    }

    for (ID=0; ID<numprocs; ID++){

      sprintf(file_ev,"%s%s.Dos.vec%d",filepath,filename,ID);

      if ( 1<=num_allocated_k[ID] ){ 

	if ( (fp_ev=fopen(file_ev,"r"))==NULL ) {
	  printf("cannot open a file %s\n",file_ev);
	}
 
        for (k=0; k<num_allocated_k[ID]; k++){
	  for (spin=0; spin<=SpinP_switch; spin++){
   	    for (l=iemin; l<=iemax; l++){
	      fread( i_vec,  sizeof(int),  3,fp_ev);
	      fwrite(i_vec,  sizeof(int),  3,fp_ev0);
	      fread( &SD[1], sizeof(float),n,fp_ev);
	      fwrite(&SD[1], sizeof(float),n,fp_ev0);
	    }
	  }
	}

	fclose(fp_ev); 
      }
    }

    fclose(fp_ev0); 
  }  

  MPI_Barrier(mpi_comm_level1);

  /* delete files */
  if (myid==Host_ID){
    for (ID=0; ID<numprocs; ID++){
      sprintf(file_ev,"%s%s.Dos.vec%d",filepath,filename,ID);
      remove(file_ev);
    }  
  }
  
  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);
  free(arpo);
  free(num_allocated_k);

  for (i=0; i<i_ko[0]; i++){
    free(ko[i]);
  }
  free(ko);

  free(koS);

  for (i=0; i<i_H[0]; i++){
    free(H[i]);
  }
  free(H);

  for (i=0; i<i_S[0]; i++){
    free(S[i]);
  }
  free(S);

  free(M1);

  for (i=0; i<i_C[0]; i++){
    free(C[i]);
  }
  free(C);

  free(KGrids1); free(KGrids2);free(KGrids3);
  
  free(SD);

  for (j=0; j<n+1; j++){
    free(TmpM[j]);
  }
  free(TmpM);

  for (j=0; j<n+1; j++){
    free(TmpS[j]);
  }
  free(TmpS);

  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<n+1; j++){
      free(H2[i][j]);
    }
    free(H2[i]);
  }
  free(H2);

  /* no spin-orbit coupling */
  if (SO_switch==0){
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  free(T_KGrids1);
  free(T_KGrids2);
  free(T_KGrids3);

  free(Ti_KGrids1);
  free(Tj_KGrids2);
  free(Tk_KGrids3);

  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<T_knum; j++){
      free(EIGEN[i][j]);
    }
    free(EIGEN[i]);
  }
  free(EIGEN);

  /* for elapsed time */

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;

}







static double Band_DFT_Dosout_NonCol(
                           int knum_i, int knum_j, int knum_k,
                           int SpinP_switch,
                           double *****nh,
                           double *****ImNL,
                           double ****CntOLP)
{
  int i,j,k,spin,l,i1,j1,n1;
  int n, wanA,ii1,jj1,m,n2;
  int *MP;
  int MA_AN, GA_AN, tnoA,Anum, LB_AN;
  int GB_AN, wanB, tnoB, Bnum, RnB;
  int l1,l2,l3;
 
  int ie,iemin,iemax,iemin0,iemax0,n1min,mn;
  int MaxL,e1,s1,T_knum,S_knum,E_knum;
  int kloop,kloop0,num_kloop0;
  int i_vec[10];

  double EV_cut0;
  double sum_r0,sum_i0,sum_r1,sum_i1;
  double sum,sumi,u2,v2,uv,vu,tmp;
  double kRn,si,co,Redum,Imdum,Redum2,Resum,Imdum2;
  double theta,phi,sit,cot,sip,cop,tmp1,tmp2,tmp3;
  double TStime,TEtime,time0;
  double OLP_eigen_cut=Threshold_OLP_Eigen;

  double *ko; int N_ko, i_ko[10];
  double *koS;
  double d0,d1,d2,d3;
  dcomplex **H;  int N_H,  i_H[10];
  dcomplex **S;  int N_S,  i_S[10];
  dcomplex **TmpS;
  dcomplex Ctmp1,Ctmp2;
  double *M1,**EIGEN;
  dcomplex **C;  int N_C,  i_C[10];
  double *KGrids1, *KGrids2, *KGrids3;
  float *SD; int N_SD, i_SD[10];
  dcomplex **TmpM;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  int *Ti_KGrids1,*Tj_KGrids2,*Tk_KGrids3,*arpo;
  int *num_allocated_k;

  double *****Dummy_ImNL;

  double snum_i, snum_j, snum_k; 
  double k1,k2,k3;

  char file_ev[YOUSO10],file_ev0[YOUSO10],file_eig[YOUSO10];
  FILE *fp_ev,*fp_ev0,*fp_eig;

  char buf1[fp_bsize];          /* setvbuf */
  char buf2[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID,ID1,ID2,tag;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;
 
  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){
    printf("\n<Band_DFT_Dosout>\n"); fflush(stdout);
  }

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
   Allocation of arrays
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);

  arpo = (int*)malloc(sizeof(int)*numprocs);
  num_allocated_k = (int*)malloc(sizeof(int)*numprocs);

  ko = (double*)malloc(sizeof(double)*n2);
  koS = (double*)malloc(sizeof(double)*(n+1));

  H = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    H[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  S = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (i=0; i<(n+1); i++){
    S[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  TmpS = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (i=0; i<(n+1); i++){
    TmpS[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  M1 = (double*)malloc(sizeof(double)*n2);

  C = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    C[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  KGrids1 = (double*)malloc(sizeof(double)*knum_i);
  KGrids2 = (double*)malloc(sizeof(double)*knum_j);
  KGrids3 = (double*)malloc(sizeof(double)*knum_k);

  SD = (float*)malloc(sizeof(float)*(atomnum+1)*List_YOUSO[7]*2);

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

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
    set up k-grids
  ****************************************************/

  snum_i = knum_i;
  snum_j = knum_j;
  snum_k = knum_k;

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
    for (i=0;i<=knum_i-1;i++) printf("%f ",KGrids1[i]);
    printf("\n");
    printf(" KGrids2: ");
    for (i=0;i<=knum_j-1;i++) printf("%f ",KGrids2[i]);
    printf("\n");
    printf(" KGrids3: ");
    for (i=0;i<=knum_k-1;i++) printf("%f ",KGrids3[i]);
    printf("\n");
  }

  /***********************************
       one-dimentionalize for MPI
  ************************************/

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        T_knum++;
      }
    }
  }

  T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);

  Ti_KGrids1 = (int*)malloc(sizeof(int)*T_knum);
  Tj_KGrids2 = (int*)malloc(sizeof(int)*T_knum);
  Tk_KGrids3 = (int*)malloc(sizeof(int)*T_knum);

  EIGEN = (double**)malloc(sizeof(double*)*T_knum);
  for (j=0; j<T_knum; j++){
    EIGEN[j] = (double*)malloc(sizeof(double)*n2);
  }

  /* set T_KGrid1,2,3 */

  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){

	T_KGrids1[T_knum] = KGrids1[i];
	T_KGrids2[T_knum] = KGrids2[j];
	T_KGrids3[T_knum] = KGrids3[k];

	Ti_KGrids1[T_knum] = i;
	Tj_KGrids2[T_knum] = j;
	Tk_KGrids3[T_knum] = k;

        T_knum++;
      }
    }
  }

  /* allocate k-points into proccessors */

  if (T_knum<=myid){
    S_knum = -10;
    E_knum = -100;
    num_kloop0 = 1;
    num_allocated_k[myid] = 0;
  }
  else if (T_knum<numprocs) {
    S_knum = myid;
    E_knum = myid;
    num_kloop0 = 1;
    num_allocated_k[myid] = 1;
  }
  else {
    tmp = (double)T_knum/(double)numprocs; 
    num_kloop0 = (int)tmp + 1;

    S_knum = (int)((double)myid*(tmp+0.0001)); 
    E_knum = (int)((double)(myid+1)*(tmp+0.0001)) - 1;
    if (myid==(numprocs-1)) E_knum = T_knum - 1;
    if (E_knum<0)           E_knum = 0;
    num_allocated_k[myid] = E_knum - S_knum + 1;
  }

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&num_allocated_k[ID], 1, MPI_INT, ID, mpi_comm_level1);
  }

  /****************************************************************
                      find iemin and iemax
  *****************************************************************/

  iemin=n2; iemax=1; n1min=1;

  k1 = 0.0;
  k2 = 0.0;
  k3 = 0.0;

  Overlap_Band(Host_ID,CntOLP,S,MP,k1,k2,k3);

  if (myid==Host_ID){

    n = S[0][0].r;
    EigenBand_lapack(S, ko, n, 1);

    if (2<=level_stdout){
      printf("  k1 k2 k3 %10.6f %10.6f %10.6f\n",k1,k2,k3);
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

  } /* if (myid==Host_ID) */

  /****************************************************
             make a full Hamiltonian matrix
  ****************************************************/

  /* non-spin-orbit coupling and non-LDA+U */
  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
      && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0)
    Hamiltonian_Band_NC(Host_ID, nh, Dummy_ImNL, H, MP, k1, k2, k3);

  /* spin-orbit coupling or LDA+U */
  else
    Hamiltonian_Band_NC(Host_ID, nh, ImNL, H, MP, k1, k2, k3);

  if (myid==Host_ID){

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

    for (i1=1; i1<=2*n; i1++){

      jj1 = 1;
      for (j1=1; j1<=n; j1++){

	sum_r0 = 0.0;
	sum_i0 = 0.0;

	sum_r1 = 0.0;
	sum_i1 = 0.0;

	for (l=1; l<=n; l++){

	  sum_r0 += H[i1][l  ].r*S[j1][l].r - H[i1][l  ].i*S[j1][l].i;
	  sum_i0 += H[i1][l  ].r*S[j1][l].i + H[i1][l  ].i*S[j1][l].r;

	  sum_r1 += H[i1][l+n].r*S[j1][l].r - H[i1][l+n].i*S[j1][l].i;
	  sum_i1 += H[i1][l+n].r*S[j1][l].i + H[i1][l+n].i*S[j1][l].r;
	}

	C[jj1][i1].r = sum_r0;
	C[jj1][i1].i = sum_i0;  jj1++;

	C[jj1][i1].r = sum_r1;
	C[jj1][i1].i = sum_i1;  jj1++;

      }
    }     

    /* M1 * U^+ H * U * M1 */

    for (j1=1; j1<=2*n; j1++){

      ii1 = 1;
      for (i1=1; i1<=n; i1++){

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

	H[ii1][j1].r = sum_r0;
	H[ii1][j1].i = sum_i0; ii1++;

	H[ii1][j1].r = sum_r1;
	H[ii1][j1].i = sum_i1; ii1++;

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

    if (n1min<n1) n1min=n1;

    iemin0=1;
    for (i1=1;i1<=n1;i1++) {
      if (ko[i1]> (ChemP+Dos_Erange[0])) {
	iemin0=i1-1;
	break;
      }
    }
    if (iemin0<1) iemin0=1;

    iemax0=n1;
    for (i1=iemin0;i1<=n1;i1++) {
      if (ko[i1]> (ChemP+Dos_Erange[1])) {
	iemax0=i1;
	break;
      }
    }
    if (iemax0>n1) iemax0=n1;

    if (iemin>iemin0) iemin=iemin0;
    if (iemax<iemax0) iemax=iemax0;

  }   /* if (myid==Host_ID) */

  /* add a buffer to iemin and iemax */

  iemin -= 5;
  iemax += 5;

  if (iemin<1)  iemin = 1;
  if (n2<iemax) iemax = n2;

  if (myid==Host_ID){
    if (n1min<iemax) iemax=n1min;
    printf(" iemin and iemax= %d %d\n",iemin,iemax);fflush(stdout);
  }

  /* MPI, iemin, iemax */
  MPI_Barrier(mpi_comm_level1);
  MPI_Bcast(&iemin, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Bcast(&iemax, 1, MPI_INT, Host_ID, mpi_comm_level1);

  /****************************************************************
                     eigenvalues and eigenvectors
   ****************************************************************/

  if (myid==Host_ID){

    sprintf(file_eig,"%s%s.Dos.val",filepath,filename);
    if ( (fp_eig=fopen(file_eig,"w"))==NULL ) {

#ifdef xt3
      setvbuf(fp_eig,buf1,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      printf("cannot open a file %s\n",file_eig);
    }

    if ( fp_eig==NULL ) {
      goto Finishing;
    }
  }

  sprintf(file_ev,"%s%s.Dos.vec%d",filepath,filename,myid);
  if ( (fp_ev=fopen(file_ev,"w"))==NULL ) {
    printf("cannot open a file %s\n",file_ev);
  }
  if ( fp_ev==NULL ) {
    goto Finishing;
  }

  if (myid==Host_ID){

    if (fp_eig) {
      fprintf(fp_eig,"mode        1\n");
      fprintf(fp_eig,"NonCol      1\n");
      fprintf(fp_eig,"N           %d\n",n);
      fprintf(fp_eig,"Nspin       %d\n",1); /* switch to 1 */
      fprintf(fp_eig,"Erange      %lf %lf\n",Dos_Erange[0],Dos_Erange[1]);
      fprintf(fp_eig,"Kgrid       %d %d %d\n",knum_i,knum_j,knum_k);
      fprintf(fp_eig,"atomnum     %d\n",atomnum);
      fprintf(fp_eig,"<WhatSpecies\n");
      for (i=1;i<=atomnum;i++) {
        fprintf(fp_eig,"%d ",WhatSpecies[i]);
      }
      fprintf(fp_eig,"\nWhatSpecies>\n");
      fprintf(fp_eig,"SpeciesNum     %d\n",SpeciesNum);
      fprintf(fp_eig,"<Spe_Total_CNO\n");
      for (i=0;i<SpeciesNum;i++) {
        fprintf(fp_eig,"%d ",Spe_Total_CNO[i]);
      }
      fprintf(fp_eig,"\nSpe_Total_CNO>\n");

      MaxL=Supported_MaxL; 

      fprintf(fp_eig,"MaxL           %d\n",Supported_MaxL);
      fprintf(fp_eig,"<Spe_Num_CBasis\n");
      for (i=0;i<SpeciesNum;i++) {
        for (l=0;l<=MaxL;l++) {
          fprintf(fp_eig,"%d ",Spe_Num_CBasis[i][l]);
        }
        fprintf(fp_eig,"\n");
      }
      fprintf(fp_eig,"Spe_Num_CBasis>\n");

      fprintf(fp_eig,"ChemP       %lf\n",ChemP);

      fprintf(fp_eig,"<SpinAngle\n");
      for (i=1; i<=atomnum; i++) {
        fprintf(fp_eig,"%lf %lf\n",Angle0_Spin[i],Angle1_Spin[i]);
      }
      fprintf(fp_eig,"SpinAngle>\n");

      printf(" write eigenvalues\n");
      printf(" write eigenvectors\n");
    }
  }

  /* for kloop */

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S */

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

	      TmpS[i1][j1].r = TmpM[i1][j1].r;
	      TmpS[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 

      }
    }

    kloop = arpo[myid];

    if (0<=kloop){

      EigenBand_lapack(S,ko,n,1);

      if (2<=level_stdout && 0<=kloop){
	printf("  kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       kloop,T_KGrids1[kloop],T_KGrids2[kloop],T_KGrids3[kloop]);
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

    }

    /****************************************************
       make a full Hamiltonian matrix
    ****************************************************/

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){
        k1 = T_KGrids1[kloop];
        k2 = T_KGrids2[kloop];
        k3 = T_KGrids3[kloop];

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
      }
    }

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

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,jj1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (i1=1+OMPID; i1<=2*n; i1+=Nthrds){

	  jj1 = 1;
	  for (j1=1; j1<=n; j1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    for (l=1; l<=n; l++){

	      sum_r0 += H[i1][l  ].r*S[j1][l].r - H[i1][l  ].i*S[j1][l].i;
	      sum_i0 += H[i1][l  ].r*S[j1][l].i + H[i1][l  ].i*S[j1][l].r;

	      sum_r1 += H[i1][l+n].r*S[j1][l].r - H[i1][l+n].i*S[j1][l].i;
	      sum_i1 += H[i1][l+n].r*S[j1][l].i + H[i1][l+n].i*S[j1][l].r;
	    }

	    C[jj1][i1].r = sum_r0;
	    C[jj1][i1].i = sum_i0;  jj1++;

	    C[jj1][i1].r = sum_r1;
	    C[jj1][i1].i = sum_i1;  jj1++;

	  }
	}     

      } /* #pragma omp parallel */

      /* M1 * U^+ H * U * M1 */

#pragma omp parallel shared(C,S,H,n) private(OMPID,Nthrds,Nprocs,i1,ii1,j1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=2*n; j1+=Nthrds){

	  ii1 = 1;
	  for (i1=1; i1<=n; i1++){

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

	    H[ii1][j1].r = sum_r0;
	    H[ii1][j1].i = sum_i0; ii1++;

	    H[ii1][j1].r = sum_r1;
	    H[ii1][j1].i = sum_i1; ii1++;

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

      n1 = 2*n;
      EigenBand_lapack(C, ko, n1, 1);

      for (i1=1; i1<=2*n; i1++) EIGEN[kloop][i1] = ko[i1];

      /****************************************************
	   Transformation to the original eigenvectors.
	   NOTE JRCAT-244p and JAIST-2122p 
	   C = U * lambda^{-1/2} * D
      ****************************************************/

      /* transpose */
      for (i1=1; i1<=2*n; i1++){
	for (j1=i1+1; j1<=2*n; j1++){
	  Ctmp1 = C[i1][j1];
	  Ctmp2 = C[j1][i1];
	  C[i1][j1] = Ctmp2;
	  C[j1][i1] = Ctmp1;
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

#pragma omp parallel shared(C,S,H,n,n1) private(OMPID,Nthrds,Nprocs,i1,j1,l1,sum_r0,sum_i0,sum_r1,sum_i1,l)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (j1=1+OMPID; j1<=n1; j1+=Nthrds){
	  for (i1=1; i1<=n; i1++){

	    sum_r0 = 0.0;
	    sum_i0 = 0.0;

	    sum_r1 = 0.0;
	    sum_i1 = 0.0;

	    l1 = 1;
	    for (l=1; l<=n; l++){

	      sum_r0 += S[i1][l].r*C[j1][l1].r
		       -S[i1][l].i*C[j1][l1].i; 
	      sum_i0 += S[i1][l].r*C[j1][l1].i
		       +S[i1][l].i*C[j1][l1].r; l1++;

	      sum_r1 += S[i1][l].r*C[j1][l1].r
		       -S[i1][l].i*C[j1][l1].i;
	      sum_i1 += S[i1][l].r*C[j1][l1].i
	 	       +S[i1][l].i*C[j1][l1].r; l1++;
	    } 

	    H[j1][i1  ].r = sum_r0;
	    H[j1][i1  ].i = sum_i0;

	    H[j1][i1+n].r = sum_r1;
	    H[j1][i1+n].i = sum_i1;
	  }
	}

      } /* #pragma omp parallel */

    } /* if (0<=kloop) */

    /****************************************************
        store LDOS
    ****************************************************/

    kloop = arpo[myid];

    if (0<=kloop){

      k1 = T_KGrids1[kloop];
      k2 = T_KGrids2[kloop];
      k3 = T_KGrids3[kloop];

      i = Ti_KGrids1[kloop];
      j = Tj_KGrids2[kloop];
      k = Tk_KGrids3[kloop];

      for (l=iemin; l<=iemax; l++){

	/* calculate SD */

#pragma omp parallel shared(TmpS,SD,H,l,MP,Spe_Total_CNO,WhatSpecies,atomnum) private(OMPID,Nthrds,Nprocs,GA_AN,wanA,tnoA,Anum,theta,phi,sit,cot,sip,cop,d0,d1,d2,d3,GB_AN,wanB,tnoB,Bnum,i1,j1,u2,v2,uv,vu,Redum,Imdum,tmp1,tmp2,tmp3)
	{ 

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (GA_AN=1+OMPID; GA_AN<=atomnum; GA_AN+=Nthrds){

	    wanA = WhatSpecies[GA_AN];
	    tnoA = Spe_Total_CNO[wanA];
	    Anum = MP[GA_AN];
	    theta = Angle0_Spin[GA_AN];
	    phi   = Angle1_Spin[GA_AN];
	    sit = sin(theta);
	    cot = cos(theta);
	    sip = sin(phi);
	    cop = cos(phi);

            for (i1=0; i1<tnoA; i1++){

	      d0 = 0.0;
	      d1 = 0.0;
	      d2 = 0.0;
	      d3 = 0.0;

	      for (GB_AN=1; GB_AN<=atomnum; GB_AN++){

		wanB = WhatSpecies[GB_AN];
		tnoB = Spe_Total_CNO[wanB];
		Bnum = MP[GB_AN];

		for (j1=0; j1<tnoB; j1++){

		  /* Re11 */
		  u2 = H[l][Anum+i1].r*H[l][Bnum+j1].r;
		  v2 = H[l][Anum+i1].i*H[l][Bnum+j1].i;
		  uv = H[l][Anum+i1].r*H[l][Bnum+j1].i;
		  vu = H[l][Anum+i1].i*H[l][Bnum+j1].r;
		  Redum = (u2 + v2);
		  Imdum = (uv - vu);
                  d0 += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  /* Re22 */
		  u2 = H[l][Anum+i1+n].r*H[l][Bnum+j1+n].r;
		  v2 = H[l][Anum+i1+n].i*H[l][Bnum+j1+n].i;
		  uv = H[l][Anum+i1+n].r*H[l][Bnum+j1+n].i;
		  vu = H[l][Anum+i1+n].i*H[l][Bnum+j1+n].r;
		  Redum = (u2 + v2);
		  Imdum = (uv - vu);
                  d1 += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  /* Re12 */
		  u2 = H[l][Anum+i1].r*H[l][Bnum+j1+n].r;
		  v2 = H[l][Anum+i1].i*H[l][Bnum+j1+n].i;
		  uv = H[l][Anum+i1].r*H[l][Bnum+j1+n].i;
		  vu = H[l][Anum+i1].i*H[l][Bnum+j1+n].r;
		  Redum = (u2 + v2);
		  Imdum = (uv - vu);
                  d2 += Redum*TmpS[Anum+i1][Bnum+j1].r - Imdum*TmpS[Anum+i1][Bnum+j1].i;

		  /* Im12
		     conjugate complex of Im12 due to difference in the definition
		     between density matrix and charge density
		  */
 
                  d3 += -(Imdum*TmpS[Anum+i1][Bnum+j1].r + Redum*TmpS[Anum+i1][Bnum+j1].i);
                    
		} /* j1 */
	      } /* GB_AN */

	      tmp1 = 0.5*(d0 + d1);
	      tmp2 = 0.5*cot*(d0 - d1);
	      tmp3 = (d2*cop - d3*sip)*sit;

	      SD[2*(Anum-1)+i1]      = (float)(tmp1 + tmp2 + tmp3);
	      SD[2*(Anum-1)+tnoA+i1] = (float)(tmp1 - tmp2 - tmp3);

	    } /* i1 */
	  } /* GA_AN */

	} /* #pragma omp parallel */

	/*********************************************
                 writing a binary file 
	*********************************************/

	i_vec[0] = i;
	i_vec[1] = j;
	i_vec[2] = k;

	fwrite(i_vec,sizeof(int),3,fp_ev);
	fwrite(&SD[0],sizeof(float),2*n,fp_ev);

      } /* l */
    } /* if (kloop<=0) */ 

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
  }

  if (myid==Host_ID){
    if (fp_eig) {

      fprintf(fp_eig,"irange      %d %d\n",iemin,iemax);
      fprintf(fp_eig,"<Eigenvalues\n");

      for (kloop=0; kloop<T_knum; kloop++){
        for (spin=0; spin<=1; spin++){

          i = Ti_KGrids1[kloop];
          j = Tj_KGrids2[kloop];
          k = Tk_KGrids3[kloop];

          fprintf(fp_eig,"%d %d %d ",i,j,k);
	  for (ie=iemin; ie<=iemax; ie++) {
	    fprintf(fp_eig,"%lf ",EIGEN[kloop][ie]);
	  }
	  fprintf(fp_eig,"\n");

	}
      }

      fprintf(fp_eig,"Eigenvalues>\n");

    } /* fp_eig */
  } /* if (myid==Host_ID) */

Finishing: 

  if (myid==Host_ID){
    if (fp_eig) { 
      fclose(fp_eig);
    }
  }

  if (fp_ev) {
    fclose(fp_ev);
  }

  /****************************************************
     merge *.Dos.vec#
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){

    sprintf(file_ev0,"%s%s.Dos.vec",filepath,filename);
    if ( (fp_ev0=fopen(file_ev0,"w"))==NULL ) {
      printf("cannot open a file %s\n",file_ev0);
    }

    for (ID=0; ID<numprocs; ID++){

      sprintf(file_ev,"%s%s.Dos.vec%d",filepath,filename,ID);

      if ( 1<=num_allocated_k[ID] ){ 

	if ( (fp_ev=fopen(file_ev,"r"))==NULL ) {
	  printf("cannot open a file %s\n",file_ev);
	}
 
        for (k=0; k<num_allocated_k[ID]; k++){
	  for (l=iemin; l<=iemax; l++){
	    fread( i_vec,  sizeof(int),  3,fp_ev);
	    fwrite(i_vec,  sizeof(int),  3,fp_ev0);
	    fread( &SD[0], sizeof(float),2*n,fp_ev);
	    fwrite(&SD[0], sizeof(float),2*n,fp_ev0);
	  }
	}

	fclose(fp_ev); 
      }
    }

    fclose(fp_ev0); 
  }  

  MPI_Barrier(mpi_comm_level1);
  
  /* delete files */
  if (myid==Host_ID){
    for (ID=0; ID<numprocs; ID++){
      sprintf(file_ev,"%s%s.Dos.vec%d",filepath,filename,ID);
      remove(file_ev);
    }  
  }

  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);
  free(arpo);
  free(num_allocated_k);

  free(ko);
  free(koS);

  for (j=0; j<n2; j++){
    free(H[j]);
  }
  free(H);

  free(SD);

  for (i=0; i<(n+1); i++){
    free(TmpS[i]);
  }
  free(TmpS);

  free(M1);

  for (j=0; j<n2; j++){
    free(C[j]);
  }
  free(C);

  free(KGrids1); free(KGrids2);free(KGrids3);
  
  for (j=0; j<n2; j++){
    free(TmpM[j]);
  }
  free(TmpM);

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

  free(Ti_KGrids1);
  free(Tj_KGrids2);
  free(Tk_KGrids3);

  for (j=0; j<T_knum; j++){
    free(EIGEN[j]);
  }
  free(EIGEN);

  /* for elapsed time */

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}
