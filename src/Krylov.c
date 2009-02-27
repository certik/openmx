/**********************************************************************
  Krylov.c:

     Krylov.c is a subroutine to perform a Krylov subspace
     method developed by T.Ozaki.

  Log of Krylov.c

     10/June/2005  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

#define  measure_time   0
#define  error_check    0
#define  cutoff_value  Threshold_OLP_Eigen



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

#ifdef TRAN
#include "tran_prototypes.h"
#endif


static void Generate_pMatrix( int spin, int Mc_AN, double *****Hks, double ****OLP0, double **invS,
                              double *****Krylov_U, double ***Krylov_U_OLP, double **inv_RS, int *MP,
                              int *Msize2, int *Msize3, int *Msize4, int Msize2_max, 
                              double **tmpvec0, double **tmpvec1, double **tmpvec2 ); 
                              
static void Generate_pMatrix2( int spin, int Mc_AN, double *****Hks, double ****OLP0, double *****Krylov_U, 
                               int *MP, int *Msize2, int *Msize3, 
                               double **tmpvec1 );
static void Krylov_IOLP( int Mc_AN, double ****OLP0, double ***Krylov_U_OLP, double **inv_RS, 
                         int *MP, int *Msize2, int *Msize4, int Msize2_max, 
                         double **tmpvec0, double **tmpvec1 );
static void S_orthonormalize_vec( int Mc_AN, int ct_on, double **vec,
                                  double **workvec, double ****OLP0,
                                  double **tmpmat0, double *ko, double *iko, int *MP, int *Msize2 );
static void Embedding_Matrix(int spin, int Mc_AN, double *****Hks,
                             double *****Krylov_U, double ****EC_matrix, 
			     int *MP, int *Msize, int *Msize2, int *Msize3, 
                             double **tmpvec1);

static void Inverse_S_by_Cholesky(int Mc_AN, double ****OLP0, double **invS, int *MP, int NUM, double *LoS);

static void Save_DOS_Col(double ******Residues, double ****OLP0, double ***EVal, int **LO_TC, int **HO_TC);

static double Krylov_Col(char *mode,
			 int SCF_iter,
			 double *****Hks,
			 double ****OLP0,
			 double *****CDM,
			 double *****EDM,
			 double Eele0[2],
			 double Eele1[2]);





double Krylov(char *mode,
              int SCF_iter,
              double *****Hks,
              double *****ImNL,
              double ****OLP0,
              double *****CDM,
              double *****EDM,
              double Eele0[2], double Eele1[2])
{
  double time0;

  /****************************************************
         collinear without spin-orbit coupling
  ****************************************************/

  if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==0 ){
    time0 = Krylov_Col(mode,SCF_iter, Hks, OLP0, CDM, EDM, Eele0, Eele1);
  }

  /****************************************************
         collinear with spin-orbit coupling
  ****************************************************/

  else if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==1 ){
    printf("Spin-orbit coupling is not supported for collinear DFT calculations.\n");
    MPI_Finalize();
    exit(1);
  }

  /****************************************************
   non-collinear with and without spin-orbit coupling
  ****************************************************/

  else if (SpinP_switch==3){
    printf("The O(N) Krylov subspace method is not supported for non-collinear DFT calculations.\n");
    MPI_Finalize();
    exit(1);
  }

  return time0;
}
















static double Krylov_Col(char *mode,
			 int SCF_iter,
			 double *****Hks,
			 double ****OLP0,
			 double *****CDM,
			 double *****EDM,
			 double Eele0[2],
			 double Eele1[2])
{
  static int firsttime=1,recalc_firsttime=1,recalc_flag;
  int Mc_AN,Gc_AN,i,is,js,Gi,wan,wanA,wanB,Anum;
  int size1,size2,num,NUM0,NUM,NUM1,n2,Cwan,Hwan,Rn2;
  int ih,ig,ian,j,kl,jg,jan,Bnum,m,n,spin,i2,ip;
  int k,l,i1,j1,P_min,m_size,q1,q2,csize,Residues_size;
  int h_AN1,Mh_AN1,h_AN2,Gh_AN1,Gh_AN2,wan1,wan2;
  int po,po1,loopN,tno1,tno2,h_AN,Gh_AN,rl1,rl2,rl;
  int MA_AN,GA_AN,tnoA,GB_AN,tnoB,ct_on;
  int Msize2_max;
  double My_TZ,TZ,sum,FermiF,time0,srt;
  double sum00,sum10,sum20,sum30;
  double sum01,sum11,sum21,sum31;
  double sum02,sum12,sum22,sum32;
  double sum03,sum13,sum23,sum33;
  double tmp0,tmp1,tmp2,tmp3,b2,co,x0,Sx,Dx,xmin,xmax;
  double My_Num_State,Num_State,x,Dnum;
  double emin,emax,de;
  double TStime,TEtime,Stime1,Etime1;
  double Stime2,Etime2;
  double time1,time2,time3,time4,time5;
  double time6,time7,time8,time9,time10;
  double time11,time12,time13,time14,time15,time16;
  double Erange;
  double My_Eele0[2],My_Eele1[2];
  double max_x=30.0;
  double ChemP_MAX,ChemP_MIN,spin_degeneracy;
  double spetrum_radius;
  double ***EVal;
  double ******Residues;
  double ***PDOS_DC;
  double *tmp_array;
  double *tmp_array2;

  int *MP;

  /*****************************************************
  Msize:   \sum_FNAN Spe_Total_CNO  
  Msize2:  \sum_FNAN+SNAN Spe_Total_CNO  
  Msize3:  rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]
  Msize4:  rlmax_EC2[Mc_AN]*EKC_core_size[Mc_AN]
  Msize5:  dimension for the last column of Residues
  *****************************************************/

  int *Msize;
  int *Msize2;
  int *Msize3;
  int *Msize4;
  int *Msize5;
  int **LO_TC;
  int **HO_TC;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  double Stime_atom, Etime_atom;

  MPI_Status stat;
  MPI_Request request;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  MPI_Barrier(mpi_comm_level1);
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);

  if (measure_time==1){ 
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
  }

  /****************************************************
                   allocation of arrays:
  ****************************************************/

  Msize = (int*)malloc(sizeof(int)*(Matomnum+1));
  Msize2 = (int*)malloc(sizeof(int)*(Matomnum+1));
  Msize3 = (int*)malloc(sizeof(int)*(Matomnum+1));
  Msize4 = (int*)malloc(sizeof(int)*(Matomnum+1));
  Msize2_max = 0;

  /* find Msize */

  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (Mc_AN==0){
      Msize[Mc_AN] = 1;
    }
    else{
      NUM = 0;
      for (i=0; i<=FNAN[Gc_AN]; i++){
	Gi = natn[Gc_AN][i];
	wanA = WhatSpecies[Gi];
	NUM += Spe_Total_CNO[wanA];
      }
      Msize[Mc_AN] = NUM;
    }
  }

  /* find Msize2 and Msize2_max */

  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (Mc_AN==0){
      Msize2[Mc_AN] = 1;
    }
    else{
      NUM = 0;
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	Gi = natn[Gc_AN][i];
	wanA = WhatSpecies[Gi];
	NUM += Spe_Total_CNO[wanA];
      }
      Msize2[Mc_AN] = NUM;
    }
 
    if (Msize2_max<Msize2[Mc_AN]) Msize2_max = Msize2[Mc_AN] + 4;
  }

  /* find Msize3 and Msize4 */

  Msize3[0] = 1;
  Msize4[0] = 1;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    wan = WhatSpecies[Gc_AN];
    ct_on = Spe_Total_CNO[wan];
    Msize3[Mc_AN] = rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN];
    Msize4[Mc_AN] = rlmax_EC2[Mc_AN]*EKC_core_size[Mc_AN];
  }

  m_size = 0;

  EVal = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    EVal[spin] = (double**)malloc(sizeof(double*)*(Matomnum+1));

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      n2 = Msize3[Mc_AN] + 2;
      m_size += n2;
      EVal[spin][Mc_AN] = (double*)malloc(sizeof(double)*n2);
    }
  }

  if (firsttime){
  PrintMemory("Krylov: EVal",  sizeof(double)*m_size,NULL);
  }

  if (2<=level_stdout){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        printf("<Krylov> myid=%4d Mc_AN=%4d Gc_AN=%4d Msize=%4d\n",
        myid,Mc_AN,M2G[Mc_AN],Msize[Mc_AN]);
    }
  }

  /****************************************************
    allocation of arrays:

    double PDOS[SpinP_switch+1]
               [Matomnum+1]
               [n2]
  ****************************************************/

  m_size = 0;
  PDOS_DC = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    PDOS_DC[spin] = (double**)malloc(sizeof(double*)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      n2 = Msize3[Mc_AN] + 4;
      m_size += n2;
      PDOS_DC[spin][Mc_AN] = (double*)malloc(sizeof(double)*n2);
    }
  }

  if (firsttime){
  PrintMemory("Krylov: PDOS_DC",sizeof(double)*m_size,NULL);
  }

  /****************************************************
    allocation of arrays:

    int LO_TC[SpinP_switch+1][Matomnum+1]
    int HO_TC[SpinP_switch+1][Matomnum+1]
  ****************************************************/

  LO_TC = (int**)malloc(sizeof(int*)*(SpinP_switch+1));
  for (spin=0; spin<(SpinP_switch+1); spin++){
    LO_TC[spin] = (int*)malloc(sizeof(int)*(Matomnum+1));
  }

  HO_TC = (int**)malloc(sizeof(int*)*(SpinP_switch+1));
  for (spin=0; spin<(SpinP_switch+1); spin++){
    HO_TC[spin] = (int*)malloc(sizeof(int)*(Matomnum+1));
  }

  /****************************************************
    allocation of array:

    double Residues[SpinP_switch+1]
                   [Matomnum+1]
                   [FNAN[Gc_AN]+1]
                   [Spe_Total_CNO[Gc_AN]] 
                   [Spe_Total_CNO[Gh_AN]] 
                   [HO_TC-LO_TC+3]
  ****************************************************/

  Residues_size = 0;
  Residues = (double******)malloc(sizeof(double*****)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    Residues[spin] = (double*****)malloc(sizeof(double****)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
        Gc_AN = 0;
        FNAN[0] = 0;
        tno1 = 1;
      }
      else{
        wanA = WhatSpecies[Gc_AN];
        tno1 = Spe_Total_CNO[wanA];
      }

      Residues[spin][Mc_AN] =
           (double****)malloc(sizeof(double***)*(FNAN[Gc_AN]+1));

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno2 = 1;
        }
        else {
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
        }

        Residues[spin][Mc_AN][h_AN] = (double***)malloc(sizeof(double**)*tno1);
        for (i=0; i<tno1; i++){
          Residues[spin][Mc_AN][h_AN][i] = (double**)malloc(sizeof(double*)*tno2);
          /* note that it is allocated once more in the loop */
        }
      }
    }
  }

  Residues_size = 1;
  for (spin=0; spin<=SpinP_switch; spin++){
    Residues[spin][0][0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  /****************************************************
     initialize density and energy density matrices
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	Gh_AN = natn[Gc_AN][h_AN];
	wanB = WhatSpecies[Gh_AN];
	tno2 = Spe_Total_CNO[wanB];
	for (i=0; i<tno1; i++){
	  for (j=0; j<tno2; j++){
	    CDM[spin][Mc_AN][h_AN][i][j] = 0.0;
	    EDM[spin][Mc_AN][h_AN][i][j] = 0.0;
	  }
	}
      }
    }
  }

  /****************************************************
   MPI

   Hks
  ****************************************************/

  if (measure_time==1) dtime(&Stime1);

  if (SCF_iter==1){

    /***********************************
               set data size
    ************************************/

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
	tag = 999;

	/* find data size to send block data */
	if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

	  size1 = 0;
	  for (spin=0; spin<=SpinP_switch; spin++){
	    for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
	      Mc_AN = Snd_MAN[IDS][n];
	      Gc_AN = Snd_GAN[IDS][n];
	      Cwan = WhatSpecies[Gc_AN]; 
	      tno1 = Spe_Total_CNO[Cwan];
	      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
		Gh_AN = natn[Gc_AN][h_AN];        
		Hwan = WhatSpecies[Gh_AN];
		tno2 = Spe_Total_CNO[Hwan];
		size1 += tno1*tno2;
	      }
	    }
	  }

	  Snd_HFS_Size[IDS] = size1;
	  MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
	}
	else{
	  Snd_HFS_Size[IDS] = 0;
	}

	/* receiving of size of data */

	if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){

	  MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
	  Rcv_HFS_Size[IDR] = size2;
	}
	else{
	  Rcv_HFS_Size[IDR] = 0;
	}

	if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);
      }
      else{
	Snd_HFS_Size[IDS] = 0;
	Rcv_HFS_Size[IDR] = 0;
      }
    }
  }

  /***********************************
             data transfer
  ************************************/

  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /*****************************
              sending of data 
      *****************************/

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

        size1 = Snd_HFS_Size[IDS];

        /* allocation of array */

        tmp_array = (double*)malloc(sizeof(double)*size1);

        /* multidimentional array to vector array */

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  tmp_array[num] = Hks[spin][Mc_AN][h_AN][i][j];
                  num++;
                } 
              } 
	    }
          }
	}

        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){

        size2 = Rcv_HFS_Size[IDR];
        
        /* allocation of array */
        tmp_array2 = (double*)malloc(sizeof(double)*size2);
        
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          Mc_AN = S_TopMAN[IDR] - 1;  /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];

            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  Hks[spin][Mc_AN][h_AN][i][j] = tmp_array2[num];
                  num++;
		}
	      }
	    }
	  }        
	}

        /* freeing of array */
        free(tmp_array2);
      }

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array); /* freeing of array */
      }
    }
  }

  /****************************************************
   MPI

   OLP0
  ****************************************************/

  /***********************************
             set data size
  ************************************/

  if (SCF_iter==1){

    /***********************************
               data transfer
    ************************************/

    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*****************************
                sending of data 
        *****************************/

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

          size1 = Snd_HFS_Size[IDS]/(1+SpinP_switch);

          /* allocation of array */

          tmp_array = (double*)malloc(sizeof(double)*size1);

          /* multidimentional array to vector array */

          num = 0;

          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  tmp_array[num] = OLP0[Mc_AN][h_AN][i][j];
                  num++;
                } 
              } 
	    }
          }

          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
        }

        /*****************************
           receiving of block data
        *****************************/

        if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
          
          size2 = Rcv_HFS_Size[IDR]/(1+SpinP_switch);
        
         /* allocation of array */
          tmp_array2 = (double*)malloc(sizeof(double)*size2);
         
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          num = 0;
          Mc_AN = S_TopMAN[IDR] - 1; /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];

            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  OLP0[Mc_AN][h_AN][i][j] = tmp_array2[num];
                  num++;
   	        }
	      }
	    }
	  }        

          /* freeing of array */
          free(tmp_array2);

        }

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
          MPI_Wait(&request,&stat);
          free(tmp_array); /* freeing of array */
	}
      }
    }
  }

  if (measure_time==1){ 
    dtime(&Etime1);
    time1 = Etime1 - Stime1;      
  }

  /***********************************************
    for regeneration of the buffer matrix
  ***********************************************/

  if (sqrt(NormRD[0])<(0.2+0.10*(double)atomnum) && recalc_firsttime){
    recalc_flag = 1;
    recalc_firsttime = 0;
  }
  else{
    recalc_flag = 0;
  }

  if (SCF_iter==1) recalc_firsttime = 1;

  if (error_check==1){
    printf("SCF_iter=%2d recalc_firsttime=%2d\n",SCF_iter,recalc_firsttime);
  }

  if (measure_time==1) dtime(&Stime2);

#pragma omp parallel shared(List_YOUSO,Residues,EDM,CDM,HO_TC,LO_TC,ChemP,EVal,RMI1,S_G2M,EC_matrix,recalc_flag,recalc_EM,Krylov_U,SpinP_switch,EKC_core_size,rlmax_EC,rlmax_EC2,time11,time10,time9,time8,time7,time6,time5,time4,time3,time2,Hks,OLP0,EKC_Exact_invS_flag,SCF_iter,Msize3,Msize2,Msize,natn,FNAN,SNAN,Spe_Total_CNO,WhatSpecies,M2G,Matomnum) 
  {
    int OMPID,Nthrds,Nprocs;
    int Mc_AN,Gc_AN,wan,ct_on,spin;
    int ig,ian,ih,kl,jg,jan,Bnum,m,n,rl;
    int Anum,i,j,k,Gi,wanA,NUM,n2,csize,is,i2;
    int i1,rl1,js,ip,po1,tno1,h_AN,Gh_AN,wanB,tno2;
    int *MP;

    double **invS;
    double *LoS;
    double **C;
    double **H_DC,*ko;
    double ***Krylov_U_OLP;
    double **inv_RS;
    double sum00,sum10,sum20,sum30,sum;
    double tmp0,tmp1,tmp2,tmp3;
    double Erange;
    double Stime_atom,Etime_atom;
    double Stime1,Etime1;
    double **tmpvec0;
    double **tmpvec1;
    double **tmpvec2;

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* allocation of arrays */

    MP = (int*)malloc(sizeof(int)*List_YOUSO[2]);

    tmpvec0 = (double**)malloc(sizeof(double*)*EKC_core_size_max);
    for (i=0; i<EKC_core_size_max; i++){
      tmpvec0[i] = (double*)malloc(sizeof(double)*Msize2_max);
    }

    tmpvec1 = (double**)malloc(sizeof(double*)*EKC_core_size_max);
    for (i=0; i<EKC_core_size_max; i++){
      tmpvec1[i] = (double*)malloc(sizeof(double)*Msize2_max);
    }

    tmpvec2 = (double**)malloc(sizeof(double*)*EKC_core_size_max);
    for (i=0; i<EKC_core_size_max; i++){
      tmpvec2[i] = (double*)malloc(sizeof(double)*Msize2_max);
    }

    if (firsttime && OMPID==0){
      PrintMemory("Krylov: tmpvec0",sizeof(double)*EKC_core_size_max*Msize2_max,NULL);
      PrintMemory("Krylov: tmpvec1",sizeof(double)*EKC_core_size_max*Msize2_max,NULL);
      PrintMemory("Krylov: tmpvec2",sizeof(double)*EKC_core_size_max*Msize2_max,NULL);
    }

    /***********************************************
     main loop of calculation 
    ***********************************************/

    for (Mc_AN=1+OMPID; Mc_AN<=Matomnum; Mc_AN+=Nthrds){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      wan = WhatSpecies[Gc_AN];
      ct_on = Spe_Total_CNO[wan];

      /* MP array */

      Anum = 1;
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	MP[i] = Anum;
	Gi = natn[Gc_AN][i];
	wanA = WhatSpecies[Gi];
	Anum += Spe_Total_CNO[wanA];
      }
      NUM = Anum - 1;
      n2 = NUM + 40;

      /***********************************************
                 allocation of arrays:
      ***********************************************/

      if (Msize[Mc_AN]<Msize3[Mc_AN])
	csize = Msize3[Mc_AN] + 40;
      else
	csize = Msize[Mc_AN] + 40;

      H_DC = (double**)malloc(sizeof(double*)*csize);
      for (i=0; i<csize; i++){
	H_DC[i] = (double*)malloc(sizeof(double)*csize);
      }

      ko = (double*)malloc(sizeof(double)*csize);

      C = (double**)malloc(sizeof(double*)*csize);
      for (i=0; i<csize; i++){
	C[i] = (double*)malloc(sizeof(double)*csize);
      }

      /***********************************************
        calculate the inverse of overlap matrix
      ***********************************************/

      if (SCF_iter==1 && Msize3[Mc_AN]<Msize2[Mc_AN] && EKC_Exact_invS_flag==1){

	LoS = (double*)malloc(sizeof(double)*(Msize2[Mc_AN]+3)*(Msize2[Mc_AN]+3));

	invS = (double**)malloc(sizeof(double*)*(Msize2[Mc_AN]+3));
	for (i=0; i<(Msize2[Mc_AN]+3); i++){
	  invS[i] = (double*)malloc(sizeof(double)*(Msize2[Mc_AN]+3));
	}

	if (measure_time==1) dtime(&Stime1);

	Inverse_S_by_Cholesky(Mc_AN, OLP0, invS, MP, NUM, LoS); 

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time2 += Etime1 - Stime1;      
	}
      }

      else if (SCF_iter==1 && Msize3[Mc_AN]<Msize2[Mc_AN] && EKC_invS_flag==1){

	Krylov_U_OLP = (double***)malloc(sizeof(double**)*rlmax_EC2[Mc_AN]);
	for (i=0; i<rlmax_EC2[Mc_AN]; i++){
	  Krylov_U_OLP[i] = (double**)malloc(sizeof(double*)*EKC_core_size[Mc_AN]);
	  for (j=0; j<EKC_core_size[Mc_AN]; j++){
	    Krylov_U_OLP[i][j] = (double*)malloc(sizeof(double)*(Msize2[Mc_AN]+3));
	    for (k=0; k<(Msize2[Mc_AN]+3); k++)  Krylov_U_OLP[i][j][k] = 0.0;
	  }
	}  

	inv_RS = (double**)malloc(sizeof(double*)*(rlmax_EC2[Mc_AN]+1)*EKC_core_size[Mc_AN]);
	for (i=0; i<(rlmax_EC2[Mc_AN]+1)*EKC_core_size[Mc_AN]; i++){
	  inv_RS[i] = (double*)malloc(sizeof(double)*(rlmax_EC2[Mc_AN]+1)*EKC_core_size[Mc_AN]);
	}      

	if (measure_time==1) dtime(&Stime1);

	Krylov_IOLP( Mc_AN, OLP0, Krylov_U_OLP, inv_RS, MP, Msize2, Msize4, Msize2_max, tmpvec0, tmpvec1 );

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time2 += Etime1 - Stime1;      
	}
      }

      for (spin=0; spin<=SpinP_switch; spin++){

	/****************************************************
                generate a preconditioning matrix
	****************************************************/

	if (measure_time==1) dtime(&Stime1);

	if (SCF_iter==1 && Msize3[Mc_AN]<Msize2[Mc_AN]){
	  Generate_pMatrix( spin, Mc_AN, Hks, OLP0, invS, Krylov_U, Krylov_U_OLP, inv_RS, MP, 
                            Msize2, Msize3, Msize4, Msize2_max, tmpvec0, tmpvec1, tmpvec2 );
	}
	else if (SCF_iter==1){
	  Generate_pMatrix2( spin, Mc_AN, Hks, OLP0, Krylov_U, MP, Msize2, Msize3, tmpvec1 );
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time3 += Etime1 - Stime1;      
	}

	if (measure_time==1) dtime(&Stime1);

	if (recalc_EM==1 || SCF_iter<=3 || recalc_flag==1){
	  Embedding_Matrix( spin, Mc_AN, Hks, Krylov_U, EC_matrix, MP, Msize, Msize2, Msize3, tmpvec1);
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time4 += Etime1 - Stime1;      
	}

	/****************************************************
                construct the Hamiltonian matrix
	****************************************************/

	if (measure_time==1) dtime(&Stime1);

	for (i=0; i<=FNAN[Gc_AN]; i++){
	  ig = natn[Gc_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  Anum = MP[i];
	  ih = S_G2M[ig]; /* S_G2M should be used */

	  for (j=0; j<=FNAN[Gc_AN]; j++){

	    kl = RMI1[Mc_AN][i][j];
	    jg = natn[Gc_AN][j];
	    jan = Spe_Total_CNO[WhatSpecies[jg]];
	    Bnum = MP[j];

	    if (0<=kl){
	      for (m=0; m<ian; m++){
		for (n=0; n<jan; n++){
		  H_DC[Anum+m][Bnum+n] = Hks[spin][ih][kl][m][n];
		}
	      }
	    }

	    else{
	      for (m=0; m<ian; m++){
		for (n=0; n<jan; n++){
		  H_DC[Anum+m][Bnum+n] = 0.0;
		}
	      }
	    }
	  }
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time5 += Etime1 - Stime1;      
	}

	/****************************************************
                   transform u1^+ * H_DC * u1
	****************************************************/

	/* H_DC * u1 */

	if (measure_time==1) dtime(&Stime1);

	/* original version */
	/*
	  for (i=1; i<=Msize[Mc_AN]; i++){
	  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
          for (n=0; n<EKC_core_size[Mc_AN]; n++){

	  sum = 0.0;
	  for (j=1; j<=Msize[Mc_AN]; j++){
	  sum += H_DC[i][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
	  }

	  C[rl*EKC_core_size[Mc_AN]+n+1][i] = sum;
	  }      
	  }      
	  } 
	*/
     
	/* unrolling version */

	for (i=1; i<=(Msize[Mc_AN]-3); i+=4){
	  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum00 = 0.0;
	      sum10 = 0.0;
	      sum20 = 0.0;
	      sum30 = 0.0;

	      for (j=1; j<=Msize[Mc_AN]; j++){
		sum00 += H_DC[i+0][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
		sum10 += H_DC[i+1][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
		sum20 += H_DC[i+2][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
		sum30 += H_DC[i+3][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
	      }

	      C[rl*EKC_core_size[Mc_AN]+n+1][i+0] = sum00;
	      C[rl*EKC_core_size[Mc_AN]+n+1][i+1] = sum10;
	      C[rl*EKC_core_size[Mc_AN]+n+1][i+2] = sum20;
	      C[rl*EKC_core_size[Mc_AN]+n+1][i+3] = sum30;

	    }      
	  }      
	} 

	is = Msize[Mc_AN] - Msize[Mc_AN]%4 + 1;

	for (i=is; i<=Msize[Mc_AN]; i++){
	  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (j=1; j<=Msize[Mc_AN]; j++){
		sum += H_DC[i][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
	      }

	      C[rl*EKC_core_size[Mc_AN]+n+1][i] = sum;
	    }      
	  }      
	} 

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time6 += Etime1 - Stime1;      
	}

	/* u1^+ * H_DC * u1 */

	if (measure_time==1) dtime(&Stime1);

	/* original version */

	/*
	  for (rl1=0; rl1<rlmax_EC[Mc_AN]; rl1++){
	  for (m=0; m<EKC_core_size[Mc_AN]; m++){
	  for (rl2=rl1; rl2<rlmax_EC[Mc_AN]; rl2++){
	  for (n=0; n<EKC_core_size[Mc_AN]; n++){

	  sum = 0.0;

	  i2 = rl2*EKC_core_size[Mc_AN] + n + 1;

	  for (i=1; i<=Msize[Mc_AN]; i++){
	  sum += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2][i];
	  }

	  H_DC[rl1*EKC_core_size[Mc_AN]+m+1][rl2*EKC_core_size[Mc_AN]+n+1] = sum;
	  H_DC[rl2*EKC_core_size[Mc_AN]+n+1][rl1*EKC_core_size[Mc_AN]+m+1] = sum;
	  }
	  }
	  }
	  }
	*/

	/* unrolling version */

	for (i2=1; i2<=(rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]-3); i2+=4){
	  for (i1=i2; i1<=rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]; i1++){

	    rl1 = (i1-1)/EKC_core_size[Mc_AN];
	    m = (i1-1) % EKC_core_size[Mc_AN];

	    sum00 = 0.0;
	    sum10 = 0.0;
	    sum20 = 0.0;
	    sum30 = 0.0;

	    for (i=1; i<=Msize[Mc_AN]; i++){
	      /* transpose */
	      sum00 += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2+0][i];
	      sum10 += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2+1][i];
	      sum20 += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2+2][i];
	      sum30 += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2+3][i];
	    }

	    H_DC[i1][i2+0] = sum00;
	    H_DC[i2+0][i1] = sum00;

	    H_DC[i1][i2+1] = sum10;
	    H_DC[i2+1][i1] = sum10;

	    H_DC[i1][i2+2] = sum20;
	    H_DC[i2+2][i1] = sum20;

	    H_DC[i1][i2+3] = sum30;
	    H_DC[i2+3][i1] = sum30;
	  }
	}

	is = rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN] - (rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN])%4 + 1;

	for (i2=is; i2<=rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]; i2++){
	  for (i1=i2; i1<=rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]; i1++){

	    rl1 = (i1-1)/EKC_core_size[Mc_AN];
	    m = (i1-1) % EKC_core_size[Mc_AN];

	    sum = 0.0;

	    for (i=1; i<=Msize[Mc_AN]; i++){
	      /* transpose */
	      sum += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2][i];
	    }

	    H_DC[i1][i2] = sum;
	    H_DC[i2][i1] = sum;
	  }
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time7 += Etime1 - Stime1;      
	}

	/* correction for ZeroNum */

	m = (int)Krylov_U[spin][Mc_AN][0][0][0];
	for (i=1; i<=m; i++) H_DC[i][i] = 1.0e+3;

	/****************************************************
            H0 = u1^+ * H_DC * u1 + D 
	****************************************************/

	if (measure_time==1) dtime(&Stime1);

	for (i=1; i<=Msize3[Mc_AN]; i++){
	  for (j=1; j<=Msize3[Mc_AN]; j++){
	    H_DC[i][j] += EC_matrix[spin][Mc_AN][i][j];
	  }
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time8 += Etime1 - Stime1;      
	}

	/****************************************************
           diagonalize
	****************************************************/

	if (measure_time==1) dtime(&Stime1);

	Eigen_lapack(H_DC,ko,Msize3[Mc_AN],Msize3[Mc_AN]);

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time9 += Etime1 - Stime1;      
	}

	/********************************************
           back transformation of eigenvectors
                      c = u1 * b
	*********************************************/

	if (measure_time==1) dtime(&Stime1);

	for (i=1; i<=Msize[Mc_AN]; i++){
	  for (j=1; j<=Msize3[Mc_AN]; j++){
	    C[i][j] = 0.0;       
	  }
	}

	/* original version */
	/*
	  for (i=1; i<=Msize[Mc_AN]; i++){
	  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
	  for (n=0; n<EKC_core_size[Mc_AN]; n++){

	  tmp1 = Krylov_U[spin][Mc_AN][rl][n][i];
	  i1 = rl*EKC_core_size[Mc_AN] + n + 1; 

	  for (j=1; j<=Msize3[Mc_AN]; j++){
	  C[i][j] += tmp1*H_DC[i1][j];
	  }

	  }
	  }     
	  }
	*/

	/* unrolling version */

	for (i=1; i<=(Msize[Mc_AN]-3); i+=4){

	  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      i1 = rl*EKC_core_size[Mc_AN] + n + 1; 

	      tmp0 = Krylov_U[spin][Mc_AN][rl][n][i+0];
	      tmp1 = Krylov_U[spin][Mc_AN][rl][n][i+1];
	      tmp2 = Krylov_U[spin][Mc_AN][rl][n][i+2];
	      tmp3 = Krylov_U[spin][Mc_AN][rl][n][i+3];

	      for (j=1; j<=(Msize3[Mc_AN]-3); j+=4){
		C[i+0][j+0] += tmp0*H_DC[i1][j+0];
		C[i+0][j+1] += tmp0*H_DC[i1][j+1];
		C[i+0][j+2] += tmp0*H_DC[i1][j+2];
		C[i+0][j+3] += tmp0*H_DC[i1][j+3];

		C[i+1][j+0] += tmp1*H_DC[i1][j+0];
		C[i+1][j+1] += tmp1*H_DC[i1][j+1];
		C[i+1][j+2] += tmp1*H_DC[i1][j+2];
		C[i+1][j+3] += tmp1*H_DC[i1][j+3];

		C[i+2][j+0] += tmp2*H_DC[i1][j+0];
		C[i+2][j+1] += tmp2*H_DC[i1][j+1];
		C[i+2][j+2] += tmp2*H_DC[i1][j+2];
		C[i+2][j+3] += tmp2*H_DC[i1][j+3];

		C[i+3][j+0] += tmp3*H_DC[i1][j+0];
		C[i+3][j+1] += tmp3*H_DC[i1][j+1];
		C[i+3][j+2] += tmp3*H_DC[i1][j+2];
		C[i+3][j+3] += tmp3*H_DC[i1][j+3];
	      }

	      js = Msize3[Mc_AN] - Msize3[Mc_AN]%4 + 1;

	      for (j=js; j<=Msize3[Mc_AN]; j++){
		C[i+0][j] += tmp0*H_DC[i1][j];
		C[i+1][j] += tmp1*H_DC[i1][j];
		C[i+2][j] += tmp2*H_DC[i1][j];
		C[i+3][j] += tmp3*H_DC[i1][j];
	      }

	    }
	  }     
	}

	is = Msize[Mc_AN] - Msize[Mc_AN]%4 + 1;

	for (i=is; i<=Msize[Mc_AN]; i++){
	  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      tmp1 = Krylov_U[spin][Mc_AN][rl][n][i];
	      i1 = rl*EKC_core_size[Mc_AN] + n + 1; 

	      for (j=1; j<=Msize3[Mc_AN]; j++){
		C[i][j] += tmp1*H_DC[i1][j];
	      }

	    }
	  }     
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time10 += Etime1 - Stime1;      
	}

	if (measure_time==1) dtime(&Stime1);

	/***********************************************
          store eigenvalues and residues of poles
	***********************************************/

	for (i=1; i<=Msize3[Mc_AN]; i++){
	  EVal[spin][Mc_AN][i-1] = ko[i];
	}

	/******************************************************
        set an energy range (-Erange+ChemP to Erange+ChemP)
        of eigenvalues used to store the Residues.
	******************************************************/

	Erange = 0.367493245;  /* in hartree, corresponds to 10 eV */

	/***********************************************
          find LO_TC and HO_TC
	***********************************************/

	/* LO_TC */ 
	i = 0;
	ip = 1;
	po1 = 0;
	do{
	  if ( (ChemP-Erange)<EVal[spin][Mc_AN][i]){
	    ip = i;
	    po1 = 1; 
	  }
	  i++;
	} while (po1==0 && i<Msize3[Mc_AN]);

	LO_TC[spin][Mc_AN] = ip;

	/* HO_TC */ 
	i = 0;
	ip = Msize3[Mc_AN];
	po1 = 0;
	do{
	  if ( (ChemP+Erange)<EVal[spin][Mc_AN][i]){
	    ip = i;
	    po1 = 1; 
	  }
	  i++;
	} while (po1==0 && i<Msize3[Mc_AN]);

	HO_TC[spin][Mc_AN] = ip;

	/***********************************************
          store residues of poles
	***********************************************/

	wanA = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_CNO[wanA];

	for (i=0; i<tno1; i++){
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    Gh_AN = natn[Gc_AN][h_AN];
	    wanB = WhatSpecies[Gh_AN];
	    tno2 = Spe_Total_CNO[wanB];
	    Bnum = MP[h_AN];
	    for (j=0; j<tno2; j++){

	      for (i1=1; i1<LO_TC[spin][Mc_AN]; i1++){
		tmp1 = C[1+i][i1]*C[Bnum+j][i1];
		CDM[spin][Mc_AN][h_AN][i][j] += tmp1;
		EDM[spin][Mc_AN][h_AN][i][j] += tmp1*EVal[spin][Mc_AN][i1-1];
	      }            

	      /* <allocation of Residues */
	      n2 = HO_TC[spin][Mc_AN] - LO_TC[spin][Mc_AN] + 3;
	      Residues[spin][Mc_AN][h_AN][i][j] = (double*)malloc(sizeof(double)*n2);
	      Residues_size += n2;
	      /* allocation of Residues> */

	      for (i1=LO_TC[spin][Mc_AN]; i1<=HO_TC[spin][Mc_AN]; i1++){
		Residues[spin][Mc_AN][h_AN][i][j][i1-LO_TC[spin][Mc_AN]] = C[1+i][i1]*C[Bnum+j][i1];
	      }
	    }
	  }
	}

	if (measure_time==1 && OMPID==0){ 
	  dtime(&Etime1);
	  time11 += Etime1 - Stime1;      
	}

      } /* spin */

      /***********************************************
                    freeing of arrays:
      ***********************************************/

      for (i=0; i<csize; i++){
	free(H_DC[i]);
      }
      free(H_DC);

      free(ko);

      for (i=0; i<csize; i++){
	free(C[i]);
      }
      free(C);

      if (SCF_iter==1 && Msize3[Mc_AN]<Msize2[Mc_AN] && EKC_Exact_invS_flag==1){

	free(LoS);

	for (i=0; i<(Msize2[Mc_AN]+3); i++){
	  free(invS[i]);
	}
	free(invS);
      }

      else if (SCF_iter==1 && Msize3[Mc_AN]<Msize2[Mc_AN] && EKC_invS_flag==1){

	for (i=0; i<rlmax_EC2[Mc_AN]; i++){
	  for (j=0; j<EKC_core_size[Mc_AN]; j++){
	    free(Krylov_U_OLP[i][j]);
	  }
	  free(Krylov_U_OLP[i]);
	}  
	free(Krylov_U_OLP);

	for (i=0; i<(rlmax_EC2[Mc_AN]+1)*EKC_core_size[Mc_AN]; i++){
	  free(inv_RS[i]);
	}      
	free(inv_RS);
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of array */

    free(MP);

    for (i=0; i<EKC_core_size_max; i++){
      free(tmpvec0[i]);
    }
    free(tmpvec0);

    for (i=0; i<EKC_core_size_max; i++){
      free(tmpvec1[i]);
    }
    free(tmpvec1);

    for (i=0; i<EKC_core_size_max; i++){
      free(tmpvec2[i]);
    }
    free(tmpvec2);

  } /* #pragma omp parallel */

  if (measure_time==1){ 
    dtime(&Etime2);
    time16 = Etime2 - Stime2;      
  }

  if (firsttime){
  PrintMemory("Krylov: Residues",sizeof(double)*Residues_size,NULL);
  }

  if ( strcasecmp(mode,"scf")==0 ){

  /****************************************************
               calculate the projected DOS
  ****************************************************/

  if (measure_time==1) dtime(&Stime1);

#pragma omp parallel shared(time_per_atom,Residues,LO_TC,HO_TC,EDM,CDM,OLP0,natn,FNAN,PDOS_DC,Msize3,Spe_Total_CNO,WhatSpecies,M2G,SpinP_switch,Matomnum) private(OMPID,Nthrds,Nprocs,Mc_AN,Stime_atom,spin,Gc_AN,wanA,tno1,i1,i,h_AN,Gh_AN,wanB,tno2,j,tmp1,Etime_atom)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (Mc_AN=1+OMPID; Mc_AN<=Matomnum; Mc_AN+=Nthrds){

      dtime(&Stime_atom);

      for (spin=0; spin<=SpinP_switch; spin++){

	Gc_AN = M2G[Mc_AN];
	wanA = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_CNO[wanA];

	for (i1=0; i1<=(Msize3[Mc_AN]+1); i1++){
	  PDOS_DC[spin][Mc_AN][i1] = 0.0;
	}

	for (i=0; i<tno1; i++){
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    Gh_AN = natn[Gc_AN][h_AN];
	    wanB = WhatSpecies[Gh_AN];
	    tno2 = Spe_Total_CNO[wanB];
	    for (j=0; j<tno2; j++){

	      tmp1 = OLP0[Mc_AN][h_AN][i][j];

	      PDOS_DC[spin][Mc_AN][0] += tmp1*CDM[spin][Mc_AN][h_AN][i][j];
	      PDOS_DC[spin][Mc_AN][1] += tmp1*EDM[spin][Mc_AN][h_AN][i][j];

	      for (i1=0; i1<(HO_TC[spin][Mc_AN]-LO_TC[spin][Mc_AN]+1); i1++){
		PDOS_DC[spin][Mc_AN][i1+2] += Residues[spin][Mc_AN][h_AN][i][j][i1]*tmp1;
	      }

	    }            
	  }        
	}

      } /* spin */

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

  } /* #pragma omp parallel */

  if (measure_time==1){ 
    dtime(&Etime1);
    time12 += Etime1 - Stime1;      
  }

  /****************************************************
           find the total number of electrons 
  ****************************************************/

  if (measure_time==1) dtime(&Stime1);

  My_TZ = 0.0;
  for (i=1; i<=Matomnum; i++){
    Gc_AN = M2G[i];
    wan = WhatSpecies[Gc_AN];
    My_TZ += Spe_Core_Charge[wan];
  }

  /* MPI, My_TZ */

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_TZ, &TZ, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /****************************************************
                find the chemical potential
  ****************************************************/

  po = 0;
  loopN = 0;

  ChemP_MAX = 15.0;  
  ChemP_MIN =-15.0;
  if      (SpinP_switch==0) spin_degeneracy = 2.0;
  else if (SpinP_switch==1) spin_degeneracy = 1.0;

  do {
    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);

    My_Num_State = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      for (spin=0; spin<=SpinP_switch; spin++){

        dtime(&Stime_atom);

        Gc_AN = M2G[Mc_AN];

        My_Num_State += spin_degeneracy*PDOS_DC[spin][Mc_AN][0];

        for (i=0; i<(HO_TC[spin][Mc_AN]-LO_TC[spin][Mc_AN]+1); i++){

          x = (EVal[spin][Mc_AN][i-1+LO_TC[spin][Mc_AN]] - ChemP)*Beta;
          if (x<=-max_x) x = -max_x;
          if (max_x<=x)  x = max_x;
          FermiF = 1.0/(1.0 + exp(x));
          My_Num_State += spin_degeneracy*FermiF*PDOS_DC[spin][Mc_AN][i+2];
	}

        dtime(&Etime_atom);
        time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }
    }

    /* MPI, My_Num_State */

    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Num_State, &Num_State, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    Dnum = (TZ - Num_State) - system_charge;
    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<1.0e-11) po = 1;


    if (myid==Host_ID && 3<=level_stdout){
      printf("  ChemP=%15.12f TZ=%15.12f Num_state=%15.12f\n",ChemP,TZ,Num_State); 
    }

    loopN++;
  }
  while (po==0 && loopN<1000); 

  if (measure_time==1){ 
    dtime(&Etime1);
    time13 += Etime1 - Stime1;      
  }

  /****************************************************
        eigenenergy by summing up eigenvalues
  ****************************************************/

  if (measure_time==1) dtime(&Stime1);

  My_Eele0[0] = 0.0;
  My_Eele0[1] = 0.0;
  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      My_Eele0[spin] += PDOS_DC[spin][Mc_AN][1];

      for (i=0; i<(HO_TC[spin][Mc_AN]-LO_TC[spin][Mc_AN]+1); i++){

        x = (EVal[spin][Mc_AN][i-1+LO_TC[spin][Mc_AN]] - ChemP)*Beta;

        if (x<=-max_x) x = -max_x;
        if (max_x<=x)  x = max_x;
        FermiF = 1.0/(1.0 + exp(x));
        My_Eele0[spin] += FermiF*EVal[spin][Mc_AN][i-1+LO_TC[spin][Mc_AN]]*PDOS_DC[spin][Mc_AN][i+2];
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }
  }

  /* MPI, My_Eele0 */
  for (spin=0; spin<=SpinP_switch; spin++){
    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Eele0[spin], &Eele0[spin], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  }

  if (SpinP_switch==0){
    Eele0[1] = Eele0[0];
  }

  if (measure_time==1){ 
    dtime(&Etime1);
    time14 += Etime1 - Stime1;      
  }

  if (measure_time==1) dtime(&Stime1);

#pragma omp parallel shared(time_per_atom,EDM,CDM,Residues,natn,max_x,Beta,ChemP,EVal,LO_TC,HO_TC,Spe_Total_CNO,WhatSpecies,M2G,SpinP_switch,Matomnum) private(OMPID,Nthrds,Nprocs,Mc_AN,spin,Stime_atom,Gc_AN,wanA,tno1,i1,x,FermiF,h_AN,wanB,tno2,i,j,tmp1,Etime_atom)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (Mc_AN=1+OMPID; Mc_AN<=Matomnum; Mc_AN+=Nthrds){
      for (spin=0; spin<=SpinP_switch; spin++){

	dtime(&Stime_atom);

	Gc_AN = M2G[Mc_AN];
	wanA = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_CNO[wanA];

	for (i1=0; i1<(HO_TC[spin][Mc_AN]-LO_TC[spin][Mc_AN]+1); i1++){

	  x = (EVal[spin][Mc_AN][i1-1+LO_TC[spin][Mc_AN]] - ChemP)*Beta;
	  if (x<=-max_x) x = -max_x;
	  if (max_x<=x)  x = max_x;
	  FermiF = 1.0/(1.0 + exp(x));

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    Gh_AN = natn[Gc_AN][h_AN];
	    wanB = WhatSpecies[Gh_AN];
	    tno2 = Spe_Total_CNO[wanB];
	    for (i=0; i<tno1; i++){
	      for (j=0; j<tno2; j++){
		tmp1 = FermiF*Residues[spin][Mc_AN][h_AN][i][j][i1];
		CDM[spin][Mc_AN][h_AN][i][j] += tmp1;
		EDM[spin][Mc_AN][h_AN][i][j] += tmp1*EVal[spin][Mc_AN][i1-1+LO_TC[spin][Mc_AN]];
	      }
	    }
	  }
	}

	dtime(&Etime_atom);
	time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }
    }

  } /* #pragma omp parallel */

  /****************************************************
                      bond energies
  ****************************************************/

  My_Eele1[0] = 0.0;
  My_Eele1[1] = 0.0;
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];    
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];

    for (j=0; j<=FNAN[GA_AN]; j++){
      GB_AN = natn[GA_AN][j];  
      wanB = WhatSpecies[GB_AN];
      tnoB = Spe_Total_CNO[wanB];

      for (k=0; k<tnoA; k++){
	for (l=0; l<tnoB; l++){
	  for (spin=0; spin<=SpinP_switch; spin++){
	    My_Eele1[spin] += CDM[spin][MA_AN][j][k][l]*Hks[spin][MA_AN][j][k][l];
	  }
	}
      }
  
    }
  }

  /* MPI, My_Eele1 */
  MPI_Barrier(mpi_comm_level1);
  for (spin=0; spin<=SpinP_switch; spin++){
    MPI_Allreduce(&My_Eele1[spin], &Eele1[spin], 1, MPI_DOUBLE,
                   MPI_SUM, mpi_comm_level1);
  }

  if (SpinP_switch==0){
    Eele1[1] = Eele1[0];
  }

  if (3<=level_stdout && myid==Host_ID){
    printf("  Eele00=%15.12f Eele01=%15.12f\n",Eele0[0],Eele0[1]);
    printf("  Eele10=%15.12f Eele11=%15.12f\n",Eele1[0],Eele1[1]);
  }

  if (measure_time==1){ 
    dtime(&Etime1);
    time15 += Etime1 - Stime1;      
  }

  } /* if ( strcasecmp(mode,"scf")==0 ) */

  else if ( strcasecmp(mode,"dos")==0 ){
    Save_DOS_Col(Residues,OLP0,EVal,LO_TC,HO_TC);
  }

  if (measure_time==1){ 
    printf("myid=%2d time1 =%5.3f time2 =%5.3f time3 =%5.3f time4 =%5.3f time5 =%5.3f\n",
            myid,time1,time2,time3,time4,time5);
    printf("myid=%2d time6 =%5.3f time7 =%5.3f time8 =%5.3f time9 =%5.3f time10=%5.3f\n",
            myid,time6,time7,time8,time9,time10);
    printf("myid=%2d time11=%5.3f time12=%5.3f time13=%5.3f time14=%5.3f time15=%5.3f\n",
            myid,time11,time12,time13,time14,time15);
    printf("myid=%2d time16=%5.3f\n",myid,time16);
  }

  /****************************************************
    freeing of arrays:

  ****************************************************/

  free(Msize);
  free(Msize2);
  free(Msize3);
  free(Msize4);

  for (spin=0; spin<(SpinP_switch+1); spin++){
    free(LO_TC[spin]);
  }
  free(LO_TC);

  for (spin=0; spin<(SpinP_switch+1); spin++){
    free(HO_TC[spin]);
  }
  free(HO_TC);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(EVal[spin][Mc_AN]);
    }
    free(EVal[spin]);
  }
  free(EVal);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
        Gc_AN = 0;
        FNAN[0] = 0;
        tno1 = 1;
      }
      else{
        wanA = WhatSpecies[Gc_AN];
        tno1 = Spe_Total_CNO[wanA];
      }

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno2 = 1;
        }
        else {
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
        }

        for (i=0; i<tno1; i++){
          for (j=0; j<tno2; j++){
            free(Residues[spin][Mc_AN][h_AN][i][j]);
	  }
          free(Residues[spin][Mc_AN][h_AN][i]);
        }
        free(Residues[spin][Mc_AN][h_AN]);
      }
      free(Residues[spin][Mc_AN]);
    }
    free(Residues[spin]);
  }
  free(Residues);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(PDOS_DC[spin][Mc_AN]);
    }
    free(PDOS_DC[spin]);
  }
  free(PDOS_DC);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  if (measure_time==1){ 
  printf("total time=%15.12f\n",time0);
  }

  /* for PrintMemory */
  firsttime=0;

  return time0;
}





















void Generate_pMatrix( int spin, int Mc_AN, double *****Hks, double ****OLP0, double **invS,
                       double *****Krylov_U, double ***Krylov_U_OLP, double **inv_RS, int *MP, 
                       int *Msize2, int *Msize3, int *Msize4, int Msize2_max, 
                       double **tmpvec0, double **tmpvec1, double **tmpvec2 ) 
{
  int rl,rl0,rl1,ct_AN,fan,san,can,wan,ct_on,i,j;
  int n,Anum,Bnum,k,ian,ih,kl,jg,ig,jan,m,m1,n1;
  int ZeroNum,Gh_AN,wanB,m1s,is;
  int rl00,rl01,rl02,rl03,rl04,rl05,rl06,rl07;
  int mm0,mm1,mm2,mm3,mm4,mm5,mm6,mm7;

  double time1,time2,time3,time4,time5;
  double time6,time7,time8,time9,time10;
  double Stime1,Etime1;
  double sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7;
  double sum,dum,tmp0,tmp1,tmp2,tmp3,rcutA,r0;
  double **Utmp,**matRS0,**matRS1;
  double **tmpmat0;
  double *ko,*iko;
  double **FS;
  double ***U0;
  int numprocs,myid,ID,tag=999;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  ct_AN = M2G[Mc_AN];
  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];
  rcutA = Spe_Atom_Cut1[wan]; 

  /* allocation of arrays */

  Utmp = (double**)malloc(sizeof(double*)*rlmax_EC[Mc_AN]);
  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    Utmp[i] = (double*)malloc(sizeof(double)*EKC_core_size[Mc_AN]);
  }

  U0 = (double***)malloc(sizeof(double**)*rlmax_EC[Mc_AN]);
  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    U0[i] = (double**)malloc(sizeof(double*)*EKC_core_size[Mc_AN]);
    for (j=0; j<EKC_core_size[Mc_AN]; j++){
      U0[i][j] = (double*)malloc(sizeof(double)*(Msize2[Mc_AN]+3));
      for (k=0; k<(Msize2[Mc_AN]+3); k++)  U0[i][j][k] = 0.0;
    }
  }  

  tmpmat0 = (double**)malloc(sizeof(double*)*(EKC_core_size[Mc_AN]+4));
  for (i=0; i<(EKC_core_size[Mc_AN]+4); i++){
    tmpmat0[i] = (double*)malloc(sizeof(double)*(EKC_core_size[Mc_AN]+4));
  }

  FS = (double**)malloc(sizeof(double*)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  for (i=0; i<(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]; i++){
    FS[i] = (double*)malloc(sizeof(double)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  }

  ko = (double*)malloc(sizeof(double)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  iko = (double*)malloc(sizeof(double)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);

  matRS0 = (double**)malloc(sizeof(double*)*(EKC_core_size[Mc_AN]+2));
  for (i=0; i<(EKC_core_size[Mc_AN]+2); i++){
    matRS0[i] = (double*)malloc(sizeof(double)*(Msize4[Mc_AN]+3));
  }

  matRS1 = (double**)malloc(sizeof(double*)*(Msize4[Mc_AN]+3));
  for (i=0; i<(Msize4[Mc_AN]+3); i++){
    matRS1[i] = (double*)malloc(sizeof(double)*(EKC_core_size[Mc_AN]+2));
  }

  /****************************************************
      initialize 
  ****************************************************/

  if (measure_time==1){ 
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
  }

  if (measure_time==1) dtime(&Stime1);

  for (i=0; i<EKC_core_size_max; i++){
    for (j=0; j<Msize2_max; j++){
      tmpvec0[i][j] = 0.0;
      tmpvec1[i][j] = 0.0;
    }
  }

  /* find the nearest atom with distance of r0 */   

  r0 = 1.0e+10;
  for (k=1; k<=FNAN[ct_AN]; k++){
    Gh_AN = natn[ct_AN][k];
    wanB = WhatSpecies[Gh_AN];
    if (Dis[ct_AN][k]<r0)  r0 = Dis[ct_AN][k]; 
  }

  /* starting vector */

  m = 0;  
  for (k=0; k<=FNAN[ct_AN]; k++){

    Gh_AN = natn[ct_AN][k];
    wanB = WhatSpecies[Gh_AN];

    if ( Dis[ct_AN][k]<(scale_rc_EKC[Mc_AN]*r0) ){

      Anum = MP[k] - 1; 

      for (i=0; i<Spe_Total_CNO[wanB]; i++){

        tmpvec0[m][Anum+i] = 1.0;

        m++;  
      }
    }
  }

  S_orthonormalize_vec( Mc_AN, ct_on, tmpvec0, tmpvec1, OLP0, tmpmat0, ko, iko, MP, Msize2 );

  for (n=0; n<EKC_core_size[Mc_AN]; n++){
    for (i=0; i<Msize2[Mc_AN]; i++){
      U0[0][n][i] = tmpvec0[n][i];
    }
  }

  if (measure_time==1){ 
    dtime(&Etime1);
    time1 = Etime1 - Stime1;      
  }

  /****************************************************
           generate Krylov subspace vectors
  ****************************************************/

  for (rl=0; rl<(rlmax_EC[Mc_AN]-1); rl++){

    if (measure_time==1) dtime(&Stime1);

    /*******************************************************
                            H * |Wn)
    *******************************************************/

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=0; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=0; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
        Bnum = MP[j] - 1;

        if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += Hks[spin][ih][kl][m][k]*tmpvec0[n][Bnum+k];
	      }

              tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    if (measure_time==1){ 
      dtime(&Etime1);
      time2 += Etime1 - Stime1;      
    }

    /*******************************************************
           S^{-1} * H * |Wn)
    *******************************************************/

    if (EKC_Exact_invS_flag==1){ 

      if (measure_time==1) dtime(&Stime1);

      /* original version */

      /*
      for (n=0; n<EKC_core_size[Mc_AN]; n++){
        for (i=0; i<Msize2[Mc_AN]; i++){
	  sum = 0.0;
	  for (j=0; j<Msize2[Mc_AN]; j++){
	    sum += invS[i][j]*tmpvec1[n][j];
	  }
	  tmpvec0[n][i] = sum;
        } 
      }
      */

      /* unrolling version */

      for (i=0; i<(Msize2[Mc_AN]-3); i+=4){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){

	  sum0 = 0.0;
	  sum1 = 0.0;
	  sum2 = 0.0;
	  sum3 = 0.0;

	  for (j=0; j<Msize2[Mc_AN]; j++){
	    sum0 += invS[i+0][j]*tmpvec1[n][j];
	    sum1 += invS[i+1][j]*tmpvec1[n][j];
	    sum2 += invS[i+2][j]*tmpvec1[n][j];
	    sum3 += invS[i+3][j]*tmpvec1[n][j];
	  }

	  tmpvec0[n][i+0] = sum0;
	  tmpvec0[n][i+1] = sum1;
	  tmpvec0[n][i+2] = sum2;
	  tmpvec0[n][i+3] = sum3;
        } 
      }

      is = Msize2[Mc_AN] - Msize2[Mc_AN]%4;

      for (i=is; i<Msize2[Mc_AN]; i++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  sum = 0.0;
	  for (j=0; j<Msize2[Mc_AN]; j++){
	    sum += invS[i][j]*tmpvec1[n][j];
	  }
	  tmpvec0[n][i] = sum;
        } 
      }

      if (measure_time==1){ 
	dtime(&Etime1);
	time3 += Etime1 - Stime1;      
      }
    }

    /*******************************************************
          U * RS^-1 * U^+ * H * |Wn)
    *******************************************************/

    else if (EKC_invS_flag==1){ 

      /*  U^+ * H * |Wn) */

      for (rl0=0; rl0<rlmax_EC2[Mc_AN]; rl0++){
	for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  for (m=0; m<EKC_core_size[Mc_AN]; m++){

	    sum = 0.0;
	    for (i=0; i<Msize2[Mc_AN]; i++){
	      sum += Krylov_U_OLP[rl0][n][i]*tmpvec1[m][i];            
	    }

	    /* transpose the later calcualtion */
	    matRS0[m][rl0*EKC_core_size[Mc_AN]+n] = sum;
	  }
	}  
      }

      /*  RS^-1 * U^+ * H * |Wn) */

      for (rl0=0; rl0<rlmax_EC2[Mc_AN]; rl0++){
	for (n=0; n<EKC_core_size[Mc_AN]; n++){

	  for (m=0; m<EKC_core_size[Mc_AN]; m++){

	    sum = 0.0;
	    for (i=0; i<Msize4[Mc_AN]; i++){
	      sum += inv_RS[rl0*EKC_core_size[Mc_AN]+n][i]*matRS0[m][i];
	    }

	    matRS1[rl0*EKC_core_size[Mc_AN]+n][m] = sum;
	  }
	}
      }         

      /*  U * RS^-1 * U^+ * H * |Wn) */

      for (n=0; n<EKC_core_size[Mc_AN]; n++){
	for (i=0; i<Msize2[Mc_AN]; i++){
	  tmpvec0[n][i] = 0.0;
	}
      }

      for (rl0=0; rl0<rlmax_EC2[Mc_AN]; rl0++){
	for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  for (m=0; m<EKC_core_size[Mc_AN]; m++){
	    tmp0 = matRS1[rl0*EKC_core_size[Mc_AN]+n][m];
	    for (i=0; i<Msize2[Mc_AN]; i++){
	      tmpvec0[m][i] += Krylov_U_OLP[rl0][n][i]*tmp0;
	    }
	  }
	}
      }
    }

    else {
      for (n=0; n<EKC_core_size[Mc_AN]; n++){
        for (i=0; i<Msize2[Mc_AN]; i++){
	  tmpvec0[n][i]  = tmpvec1[n][i];
        }
      }
    }

    if (measure_time==1) dtime(&Stime1);

    /*************************************************************
     S-orthogonalization by a classical block Gram-Schmidt method 
    *************************************************************/

    /* |tmpvec2) = S * |tmpvec0) */

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec2[n][i] = 0.0;
      }
    }

    for (i=0; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=0; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[j] - 1;

	if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += OLP0[ih][kl][m][k]*tmpvec0[n][Bnum+k];
	      }

	      tmpvec2[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    if (measure_time==1){ 
      dtime(&Etime1);
      time4 += Etime1 - Stime1;      
    }

    if (measure_time==1) dtime(&Stime1);

    for (rl0=0; rl0<=rl; rl0++){

      /* (U_rl0|tmpvec2) */

      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
          sum = 0.0;
  	  for (i=0; i<Msize2[Mc_AN]; i++){
            sum += U0[rl0][m][i]*tmpvec2[n][i];
	  }
          tmpmat0[m][n] = sum;
	}
      }

      /* |tmpvec0) - |U_rl0) * (U_rl0|tmpvec2) */

      for (n=0; n<EKC_core_size[Mc_AN]; n++){
	for (k=0; k<EKC_core_size[Mc_AN]; k++){
	  dum = tmpmat0[k][n];
	  for (i=0; i<Msize2[Mc_AN]; i++)  tmpvec0[n][i] -= U0[rl0][k][i]*dum;
        }
      }

    }

    if (measure_time==1){ 
      dtime(&Etime1);
      time5 += Etime1 - Stime1;      
    }

    /*************************************************************
                   S-orthonormalization of tmpvec0
    *************************************************************/

    if (measure_time==1) dtime(&Stime1);

    S_orthonormalize_vec( Mc_AN, ct_on, tmpvec0, tmpvec1, OLP0, tmpmat0, ko, iko, MP, Msize2 );

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
        U0[rl+1][n][i] = tmpvec0[n][i];
      }
    }

    if (measure_time==1){ 
      dtime(&Etime1);
      time6 += Etime1 - Stime1;      
    }

  } /* rl */

  /************************************************************
              orthogonalization by diagonalization
  ************************************************************/

  if (measure_time==1) dtime(&Stime1);

  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

    /*  S * |Vn) */

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=0; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=0; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[j] - 1;

	if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += OLP0[ih][kl][m][k]*U0[rl][n][Bnum+k];
	      }
	      tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    for (rl0=rl; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
  	  sum = 0.0;
	  for (i=0; i<Msize2[Mc_AN]; i++){
            sum += U0[rl0][m][i]*tmpvec1[n][i]; 
	  }
          FS[rl0*EKC_core_size[Mc_AN]+m+1][rl*EKC_core_size[Mc_AN]+n+1] = sum;
          FS[rl*EKC_core_size[Mc_AN]+n+1][rl0*EKC_core_size[Mc_AN]+m+1] = sum;
	}
      }
    }
  }

  if (measure_time==1){ 
    dtime(&Etime1);
    time7 += Etime1 - Stime1;      
  }

  if (measure_time==1) dtime(&Stime1);

  Eigen_lapack(FS,ko,Msize3[Mc_AN],Msize3[Mc_AN]);

  if (measure_time==1){ 
    dtime(&Etime1);
    time8 += Etime1 - Stime1;      
  }

  ZeroNum = 0;

  for (i=1; i<=Msize3[Mc_AN]; i++){

    if (error_check==1){
      printf("spin=%2d Mc_AN=%2d i=%3d ko[i]=%18.15f\n",spin,Mc_AN,i,ko[i]);
    }

    if (cutoff_value<ko[i]){
      ko[i]  = sqrt(ko[i]);
      iko[i] = 1.0/ko[i];
    }
    else{
      ZeroNum++;
      ko[i]  = 0.0;
      iko[i] = 0.0;
    }
  }  

  if (error_check==1){
    printf("spin=%2d Mc_AN=%2d ZeroNum=%2d\n",spin,Mc_AN,ZeroNum);
  }

  for (i=1; i<=Msize3[Mc_AN]; i++){
    for (j=1; j<=Msize3[Mc_AN]; j++){
      FS[i][j] = FS[i][j]*iko[j];
    }
  }

  /* transpose for later calculation */
  for (i=1; i<=Msize3[Mc_AN]; i++){
    for (j=i+1; j<=Msize3[Mc_AN]; j++){
      tmp1 = FS[i][j];
      tmp2 = FS[j][i];
      FS[i][j] = tmp2;
      FS[j][i] = tmp1;
    }
  }

  if (measure_time==1) dtime(&Stime1);

  /* U0 * U * lamda^{-1/2} */ 

  /* original version */

  /*
  for (i=0; i<Msize2[Mc_AN]; i++){
    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){

        m1 = rl0*EKC_core_size[Mc_AN] + m + 1;

	sum = 0.0; 
	for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

          n1 = rl*EKC_core_size[Mc_AN] + 1;

	  for (n=0; n<EKC_core_size[Mc_AN]; n++){
	    sum += U0[rl][n][i]*FS[m1][n1+n]; 
	  }
	}

        Utmp[rl0][m] = sum;        
      }
    }

    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        U0[rl0][m][i] = Utmp[rl0][m];
      }
    }
  }
  */

  /* unrolling version */

  for (i=0; i<Msize2[Mc_AN]; i++){

    for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
      for (n=0; n<EKC_core_size[Mc_AN]; n++){
        Utmp[rl][n] = U0[rl][n][i];
      }
    }

    for (m1=1; m1<=(rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]-3); m1+=4){

      rl00 = (m1+0-1)/EKC_core_size[Mc_AN];
      rl01 = (m1+1-1)/EKC_core_size[Mc_AN];
      rl02 = (m1+2-1)/EKC_core_size[Mc_AN];
      rl03 = (m1+3-1)/EKC_core_size[Mc_AN];

      mm0 = (m1+0-1)%EKC_core_size[Mc_AN];
      mm1 = (m1+1-1)%EKC_core_size[Mc_AN];
      mm2 = (m1+2-1)%EKC_core_size[Mc_AN];
      mm3 = (m1+3-1)%EKC_core_size[Mc_AN];

      sum0 = 0.0; 
      sum1 = 0.0; 
      sum2 = 0.0; 
      sum3 = 0.0; 

      for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

	n1 = rl*EKC_core_size[Mc_AN] + 1;

	for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  sum0 += Utmp[rl][n]*FS[m1+0][n1+n]; 
	  sum1 += Utmp[rl][n]*FS[m1+1][n1+n]; 
	  sum2 += Utmp[rl][n]*FS[m1+2][n1+n]; 
	  sum3 += Utmp[rl][n]*FS[m1+3][n1+n]; 
	}
      }

      U0[rl00][mm0][i] = sum0;        
      U0[rl01][mm1][i] = sum1;        
      U0[rl02][mm2][i] = sum2;        
      U0[rl03][mm3][i] = sum3;        
    }

    m1s = rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN] - (rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN])%4 + 1;

    for (m1=m1s; m1<=rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]; m1++){

      rl0 = (m1-1)/EKC_core_size[Mc_AN];
      m = (m1-1)%EKC_core_size[Mc_AN];

      sum = 0.0; 

      for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

	n1 = rl*EKC_core_size[Mc_AN] + 1;

	for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  sum += Utmp[rl][n]*FS[m1][n1+n]; 
	}
      }

      U0[rl0][m][i] = sum;        
    }
  } /* i */

  /*
  for (i=0; i<Msize2[Mc_AN]; i++){

    for (m1=1; m1<=(rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]-7); m1+=8){

      rl00 = (m1+0-1)/EKC_core_size[Mc_AN];
      rl01 = (m1+1-1)/EKC_core_size[Mc_AN];
      rl02 = (m1+2-1)/EKC_core_size[Mc_AN];
      rl03 = (m1+3-1)/EKC_core_size[Mc_AN];
      rl04 = (m1+4-1)/EKC_core_size[Mc_AN];
      rl05 = (m1+5-1)/EKC_core_size[Mc_AN];
      rl06 = (m1+6-1)/EKC_core_size[Mc_AN];
      rl07 = (m1+7-1)/EKC_core_size[Mc_AN];

      mm0 = (m1+0-1)%EKC_core_size[Mc_AN];
      mm1 = (m1+1-1)%EKC_core_size[Mc_AN];
      mm2 = (m1+2-1)%EKC_core_size[Mc_AN];
      mm3 = (m1+3-1)%EKC_core_size[Mc_AN];
      mm4 = (m1+4-1)%EKC_core_size[Mc_AN];
      mm5 = (m1+5-1)%EKC_core_size[Mc_AN];
      mm6 = (m1+6-1)%EKC_core_size[Mc_AN];
      mm7 = (m1+7-1)%EKC_core_size[Mc_AN];

      sum0 = 0.0; 
      sum1 = 0.0; 
      sum2 = 0.0; 
      sum3 = 0.0; 
      sum4 = 0.0; 
      sum5 = 0.0; 
      sum6 = 0.0; 
      sum7 = 0.0; 

      for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

	n1 = rl*EKC_core_size[Mc_AN] + 1;

	for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  sum0 += U0[rl][n][i]*FS[m1+0][n1+n]; 
	  sum1 += U0[rl][n][i]*FS[m1+1][n1+n]; 
	  sum2 += U0[rl][n][i]*FS[m1+2][n1+n]; 
	  sum3 += U0[rl][n][i]*FS[m1+3][n1+n]; 
	  sum4 += U0[rl][n][i]*FS[m1+4][n1+n]; 
	  sum5 += U0[rl][n][i]*FS[m1+5][n1+n]; 
	  sum6 += U0[rl][n][i]*FS[m1+6][n1+n]; 
	  sum7 += U0[rl][n][i]*FS[m1+7][n1+n]; 
	}
      }

      Utmp[rl00][mm0] = sum0;        
      Utmp[rl01][mm1] = sum1;        
      Utmp[rl02][mm2] = sum2;        
      Utmp[rl03][mm3] = sum3;        
      Utmp[rl04][mm4] = sum4;        
      Utmp[rl05][mm5] = sum5;        
      Utmp[rl06][mm6] = sum6;        
      Utmp[rl07][mm7] = sum7;        
    }

    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        U0[rl00][mm0][i] = Utmp[rl00][mm0];
        U0[rl01][mm1][i] = Utmp[rl01][mm1];
        U0[rl02][mm2][i] = Utmp[rl02][mm2];
        U0[rl03][mm3][i] = Utmp[rl03][mm3];
        U0[rl04][mm4][i] = Utmp[rl04][mm4];
        U0[rl05][mm5][i] = Utmp[rl05][mm5];
        U0[rl06][mm6][i] = Utmp[rl06][mm6];
        U0[rl07][mm7][i] = Utmp[rl07][mm7];
      }
    }

    m1s = rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN] - (rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN])%8 + 1;

    for (m1=m1s; m1<=rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]; m1++){

      rl0 = (m1-1)/EKC_core_size[Mc_AN];
      m = (m1-1)%EKC_core_size[Mc_AN];

      sum = 0.0; 

      for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

	n1 = rl*EKC_core_size[Mc_AN] + 1;

	for (n=0; n<EKC_core_size[Mc_AN]; n++){
	  sum += U0[rl][n][i]*FS[m1][n1+n]; 
	}
      }

      Utmp[rl0][m] = sum;        
    }

    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        U0[rl0][m][i] = Utmp[rl0][m];
      }
    }

  }
  */


  Krylov_U[spin][Mc_AN][0][0][0] = ZeroNum;

  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
        Krylov_U[spin][Mc_AN][rl][n][i+1] = U0[rl][n][i];
      }
    }
  }

  if (measure_time==1){ 
    dtime(&Etime1);
    time9 += Etime1 - Stime1;      
  }

  /************************************************************
         check the orthonormality of Krylov vectors
  ************************************************************/

  if (error_check==1){

    for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

      for (n=0; n<EKC_core_size[Mc_AN]; n++){
	for (i=0; i<Msize2[Mc_AN]; i++){
	  tmpvec1[n][i]  = 0.0;
	}
      }

      for (i=0; i<=can; i++){

	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];
	Anum = MP[i] - 1;
	ih = S_G2M[ig];

	for (j=0; j<=can; j++){

	  kl = RMI1[Mc_AN][i][j];
	  jg = natn[ct_AN][j];
	  jan = Spe_Total_CNO[WhatSpecies[jg]];
	  Bnum = MP[j] - 1;

	  if (0<=kl){

	    for (m=0; m<ian; m++){
	      for (n=0; n<EKC_core_size[Mc_AN]; n++){

		sum = 0.0;
		for (k=0; k<jan; k++){
		  sum += OLP0[ih][kl][m][k]*U0[rl][n][Bnum+k];
		}

		tmpvec1[n][Anum+m] += sum;
	      }
	    } 
	  }
	}
      }

      for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
	for (m=0; m<EKC_core_size[Mc_AN]; m++){
	  for (n=0; n<EKC_core_size[Mc_AN]; n++){
	    sum = 0.0;
	    for (i=0; i<Msize2[Mc_AN]; i++){
	      sum += U0[rl0][m][i]*tmpvec1[n][i]; 
	    }

	    if (rl==rl0 && m==n){
	      if ( 1.0e-10<fabs(sum-1.0) ) {
		printf("A spin=%2d Mc_AN=%2d rl=%2d rl0=%2d m=%2d n=%2d sum=%18.15f\n",
		       spin,Mc_AN,rl,rl0,m,n,sum);
	      }
	    }
	    else{
	      if ( 1.0e-10<fabs(sum) ) {
		printf("B spin=%2d Mc_AN=%2d rl=%2d rl0=%2d m=%2d n=%2d sum=%18.15f\n",
		       spin,Mc_AN,rl,rl0,m,n,sum);
	      }
	    }

	  }
	}
      }
    }
  }

  if (measure_time==1){ 
    printf("pMatrix myid=%2d time1 =%5.3f time2 =%5.3f time3 =%5.3f time4 =%5.3f\n",
            myid,time1,time2,time3,time4);
    printf("pMatrix myid=%2d time5 =%5.3f time6 =%5.3f time7 =%5.3f time8 =%5.3f\n",
            myid,time5,time6,time7,time8);
    printf("pMatrix myid=%2d time9 =%5.3f\n",myid,time9);

  }

  /* freeing of arrays */

  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    free(Utmp[i]);
  }
  free(Utmp);

  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    for (j=0; j<EKC_core_size[Mc_AN]; j++){
      free(U0[i][j]);
    }
    free(U0[i]);
  }  
  free(U0);

  for (i=0; i<(EKC_core_size[Mc_AN]+4); i++){
    free(tmpmat0[i]);
  }
  free(tmpmat0);

  for (i=0; i<(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]; i++){
    free(FS[i]);
  }
  free(FS);

  free(ko);
  free(iko);

  for (i=0; i<(EKC_core_size[Mc_AN]+2); i++){
    free(matRS0[i]);
  }
  free(matRS0);

  for (i=0; i<(Msize4[Mc_AN]+3); i++){
    free(matRS1[i]);
  }
  free(matRS1);
}








void Generate_pMatrix2( int spin, int Mc_AN, double *****Hks, double ****OLP0, 
                        double *****Krylov_U, int *MP, int *Msize2, int *Msize3, 
                        double **tmpvec1)
{
  int rl,rl0,rl1,ct_AN,fan,san,can,wan,ct_on,i,j;
  int n,Anum,Bnum,k,ian,ih,kl,jg,ig,jan,m,m1,n1;
  int ZeroNum,rl_half;
  double sum,dum,tmp0,tmp1,tmp2,tmp3;
  double **Utmp;
  double *ko,*iko;
  double **FS;
  double ***U0;

  ct_AN = M2G[Mc_AN];
  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];

  /* allocation of arrays */

  Utmp = (double**)malloc(sizeof(double*)*rlmax_EC[Mc_AN]);
  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    Utmp[i] = (double*)malloc(sizeof(double)*EKC_core_size[Mc_AN]);
  }

  U0 = (double***)malloc(sizeof(double**)*rlmax_EC[Mc_AN]);
  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    U0[i] = (double**)malloc(sizeof(double*)*EKC_core_size[Mc_AN]);
    for (j=0; j<EKC_core_size[Mc_AN]; j++){
      U0[i][j] = (double*)malloc(sizeof(double)*(Msize2[Mc_AN]+3));
      for (k=0; k<(Msize2[Mc_AN]+3); k++)  U0[i][j][k] = 0.0;
    }
  }  

  FS = (double**)malloc(sizeof(double*)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  for (i=0; i<(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]; i++){
    FS[i] = (double*)malloc(sizeof(double)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  }

  ko = (double*)malloc(sizeof(double)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  iko = (double*)malloc(sizeof(double)*(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]);

  /****************************************************
   initialize 
  ****************************************************/

  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
        U0[rl][n][i] = 0.0;
      }
    }
  }

  i = 0;
  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      if (i<Msize2[Mc_AN]) U0[rl][n][i] = 1.0;
      i++;
    }
  }

  /************************************************************
              orthogonalization by diagonalization
  ************************************************************/

  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

    /*  S * |Vn) */

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=0; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=0; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[j] - 1;

	if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += OLP0[ih][kl][m][k]*U0[rl][n][Bnum+k];
	      }
	      tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    for (rl0=rl; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
  	  sum = 0.0;
	  for (i=0; i<Msize2[Mc_AN]; i++){
            sum += U0[rl0][m][i]*tmpvec1[n][i]; 
	  }
          FS[rl0*EKC_core_size[Mc_AN]+m+1][rl*EKC_core_size[Mc_AN]+n+1] = sum;
          FS[rl*EKC_core_size[Mc_AN]+n+1][rl0*EKC_core_size[Mc_AN]+m+1] = sum;
	}
      }
    }
  }

  Eigen_lapack(FS,ko,Msize3[Mc_AN],Msize3[Mc_AN]);

  ZeroNum = 0;

  for (i=1; i<=Msize3[Mc_AN]; i++){

    if (error_check==1){
      printf("spin=%2d Mc_AN=%2d i=%3d ko[i]=%18.15f\n",spin,Mc_AN,i,ko[i]);
    }

    if (cutoff_value<ko[i]){
      ko[i]  = sqrt(ko[i]);
      iko[i] = 1.0/ko[i];
    }
    else{
      ZeroNum++;
      ko[i]  = 0.0;
      iko[i] = 0.0;
    }
  }  

  if (error_check==1){
    printf("spin=%2d Mc_AN=%2d ZeroNum=%2d\n",spin,Mc_AN,ZeroNum);
  }

  for (i=1; i<=Msize3[Mc_AN]; i++){
    for (j=1; j<=Msize3[Mc_AN]; j++){
      FS[i][j] = FS[i][j]*iko[j];
    }
  }

  /* transpose for later calculation */
  for (i=1; i<=Msize3[Mc_AN]; i++){
    for (j=i+1; j<=Msize3[Mc_AN]; j++){
      tmp1 = FS[i][j];
      tmp2 = FS[j][i];
      FS[i][j] = tmp2;
      FS[j][i] = tmp1;
    }
  }

  /* U0 * U * lamda^{-1/2} */ 

  for (i=0; i<Msize2[Mc_AN]; i++){
    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){

        m1 = rl0*EKC_core_size[Mc_AN] + m + 1;

	sum = 0.0; 
	for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

          n1 = rl*EKC_core_size[Mc_AN] + 1;

	  for (n=0; n<EKC_core_size[Mc_AN]; n++){
	                     /* transpose */
	    sum += U0[rl][n][i]*FS[m1][n1+n]; 
	  }
	}

        Utmp[rl0][m] = sum;        
      }
    }

    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        U0[rl0][m][i] = Utmp[rl0][m];
      }
    }
  }

  Krylov_U[spin][Mc_AN][0][0][0] = ZeroNum;

  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
        Krylov_U[spin][Mc_AN][rl][n][i+1] = U0[rl][n][i];
      }
    }
  }

  /************************************************************
         check the orthonormality of Krylov vectors
  ************************************************************/

  if (error_check==1){

    for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

      for (n=0; n<EKC_core_size[Mc_AN]; n++){
	for (i=0; i<Msize2[Mc_AN]; i++){
	  tmpvec1[n][i]  = 0.0;
	}
      }

      for (i=0; i<=can; i++){

	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];
	Anum = MP[i] - 1;
	ih = S_G2M[ig];

	for (j=0; j<=can; j++){

	  kl = RMI1[Mc_AN][i][j];
	  jg = natn[ct_AN][j];
	  jan = Spe_Total_CNO[WhatSpecies[jg]];
	  Bnum = MP[j] - 1;

	  if (0<=kl){

	    for (m=0; m<ian; m++){
	      for (n=0; n<EKC_core_size[Mc_AN]; n++){

		sum = 0.0;
		for (k=0; k<jan; k++){
		  sum += OLP0[ih][kl][m][k]*U0[rl][n][Bnum+k];
		}

		tmpvec1[n][Anum+m] += sum;
	      }
	    } 
	  }
	}
      }

      for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
	for (m=0; m<EKC_core_size[Mc_AN]; m++){
	  for (n=0; n<EKC_core_size[Mc_AN]; n++){
	    sum = 0.0;
	    for (i=0; i<Msize2[Mc_AN]; i++){
	      sum += U0[rl0][m][i]*tmpvec1[n][i]; 
	    }

	    if (rl==rl0 && m==n){
	      if ( 1.0e-10<fabs(sum-1.0) ) {
		printf("A spin=%2d Mc_AN=%2d rl=%2d rl0=%2d m=%2d n=%2d sum=%18.15f\n",
		       spin,Mc_AN,rl,rl0,m,n,sum);
	      }
	    }
	    else{
	      if ( 1.0e-10<fabs(sum) ) {
		printf("B spin=%2d Mc_AN=%2d rl=%2d rl0=%2d m=%2d n=%2d sum=%18.15f\n",
		       spin,Mc_AN,rl,rl0,m,n,sum);
	      }
	    }

	  }
	}
      }
    }
  }

  /* freeing of arrays */

  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    free(Utmp[i]);
  }
  free(Utmp);

  for (i=0; i<rlmax_EC[Mc_AN]; i++){
    for (j=0; j<EKC_core_size[Mc_AN]; j++){
      free(U0[i][j]);
    }
    free(U0[i]);
  }  
  free(U0);

  for (i=0; i<(rlmax_EC[Mc_AN]+2)*EKC_core_size[Mc_AN]; i++){
    free(FS[i]);
  }
  free(FS);

  free(ko);
  free(iko);

}







void Embedding_Matrix(int spin, int Mc_AN, double *****Hks,
                      double *****Krylov_U, double ****EC_matrix, 
                      int *MP, int *Msize, int *Msize2, int *Msize3, 
                      double **tmpvec1 )
{
  int ct_AN,fan,san,can,wan,ct_on;
  int rl,rl0,m,n,i,j,k,kl,jg,jan,ih,ian;
  int Anum,Bnum,ig;
  double sum,tmp1,tmp2,tmp3;

  ct_AN = M2G[Mc_AN];
  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];

  /*******************************
            u1^+ C^+ u2 
  *******************************/
   
  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

    /* C^+ u2 */

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=0; i<=fan; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=fan+1; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[j];

	if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += Hks[spin][ih][kl][m][k]*Krylov_U[spin][Mc_AN][rl][n][Bnum+k];
	      }

	      tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    /* u1^+ C^+ u2 */

    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
  	  sum = 0.0;
	  for (i=0; i<Msize[Mc_AN]; i++){
            sum += Krylov_U[spin][Mc_AN][rl0][m][i+1]*tmpvec1[n][i]; 
	  }
          EC_matrix[spin][Mc_AN][rl0*EKC_core_size[Mc_AN]+m+1][rl*EKC_core_size[Mc_AN]+n+1] = sum;
	}
      }
    }

  } /* rl */

  /*******************************
            u2^+ C u1 
  *******************************/
  
  for (i=1; i<=Msize3[Mc_AN]; i++){
    for (j=i; j<=Msize3[Mc_AN]; j++){

      tmp1 = EC_matrix[spin][Mc_AN][i][j];
      tmp2 = EC_matrix[spin][Mc_AN][j][i];
      tmp3 = tmp1 + tmp2;   

      EC_matrix[spin][Mc_AN][i][j] = tmp3; 
      EC_matrix[spin][Mc_AN][j][i] = tmp3;
    }
  }

  /*******************************
             u2^+ B u2 
  *******************************/

  for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){

    /* B u2 */

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=fan+1; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=fan+1; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[j];

	if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += Hks[spin][ih][kl][m][k]*Krylov_U[spin][Mc_AN][rl][n][Bnum+k];
	      }

	      tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    /* u2^+ B u2 */

    for (rl0=0; rl0<rlmax_EC[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
  	  sum = 0.0;
	  for (i=Msize[Mc_AN]; i<Msize2[Mc_AN]; i++){
            sum += Krylov_U[spin][Mc_AN][rl0][m][i+1]*tmpvec1[n][i]; 
	  }
          EC_matrix[spin][Mc_AN][rl0*EKC_core_size[Mc_AN]+m+1][rl*EKC_core_size[Mc_AN]+n+1] += sum;
	}
      }
    }

  } /* rl */

}










void Krylov_IOLP( int Mc_AN, double ****OLP0, double ***Krylov_U_OLP, double **inv_RS, int *MP, 
                  int *Msize2, int *Msize4, int Msize2_max, double **tmpvec0, double **tmpvec1 )
{
  int rl,ct_AN,fan,san,can,wan,ct_on,i,j;
  int n,Anum,Bnum,k,ian,ih,kl,jg,ig,jan,m,m1,n1;
  int rl0,ZeroNum,Neumann_series,ns,Gh_AN,wanB;
  double sum,dum,tmp0,tmp1,tmp2,tmp3,rcutA,r0;
  double **tmpmat0;
  double *ko,*iko;
  double **FS;

  ct_AN = M2G[Mc_AN];
  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];
  rcutA = Spe_Atom_Cut1[wan]; 

  /* allocation of arrays */

  tmpmat0 = (double**)malloc(sizeof(double*)*(EKC_core_size[Mc_AN]+4));
  for (i=0; i<(EKC_core_size[Mc_AN]+4); i++){
    tmpmat0[i] = (double*)malloc(sizeof(double)*(EKC_core_size[Mc_AN]+4));
  }

  FS = (double**)malloc(sizeof(double*)*(rlmax_EC2[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  for (i=0; i<(rlmax_EC2[Mc_AN]+2)*EKC_core_size[Mc_AN]; i++){
    FS[i] = (double*)malloc(sizeof(double)*(rlmax_EC2[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  }

  ko = (double*)malloc(sizeof(double)*(rlmax_EC2[Mc_AN]+2)*EKC_core_size[Mc_AN]);
  iko = (double*)malloc(sizeof(double)*(rlmax_EC2[Mc_AN]+2)*EKC_core_size[Mc_AN]);

  /****************************************************
      initialize 
  ****************************************************/

  for (i=0; i<EKC_core_size_max; i++){
    for (j=0; j<Msize2_max; j++){
      tmpvec0[i][j] = 0.0;
      tmpvec1[i][j] = 0.0;
    }
  }

  /* find the nearest atom with distance of r0 */   

  r0 = 1.0e+10;
  for (k=1; k<=FNAN[ct_AN]; k++){
    Gh_AN = natn[ct_AN][k];
    wanB = WhatSpecies[Gh_AN];
    if (Dis[ct_AN][k]<r0)  r0 = Dis[ct_AN][k]; 
  }

  /* starting vector */

  m = 0;  
  for (k=0; k<=FNAN[ct_AN]; k++){

    Gh_AN = natn[ct_AN][k];
    wanB = WhatSpecies[Gh_AN];

    if ( Dis[ct_AN][k]<(scale_rc_EKC[Mc_AN]*r0) ){

      Anum = MP[k] - 1; 
      
      for (i=0; i<Spe_Total_CNO[wanB]; i++){

        tmpvec0[m][Anum+i]         = 1.0;
        Krylov_U_OLP[0][m][Anum+i] = 1.0;

        m++;  
      }
    }
  }

  /*
  for (i=0; i<EKC_core_size[Mc_AN]; i++){
    tmpvec0[i][i]         = 1.0;
    Krylov_U_OLP[0][i][i] = 1.0;
  }
  */

  /****************************************************
           generate Krylov subspace vectors
  ****************************************************/

  for (rl=0; rl<(rlmax_EC2[Mc_AN]-1); rl++){

    /*******************************************************
                            S * |Wn)
    *******************************************************/

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=0; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=0; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
        Bnum = MP[j] - 1;

        if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += OLP0[ih][kl][m][k]*tmpvec0[n][Bnum+k];
	      }

              tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    /*************************************************************
      orthogonalization by a modified block Gram-Schmidt method 
    *************************************************************/

    /* |tmpvec1) := (I - \sum_{rl0} |U_rl0)(U_rl0|)|tmpvec1) */

    for (rl0=0; rl0<=rl; rl0++){

      /* (U_rl0|tmpvec1) */

      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
          sum = 0.0;
  	  for (i=0; i<Msize2[Mc_AN]; i++){
            sum += Krylov_U_OLP[rl0][m][i]*tmpvec1[n][i];
	  }
          tmpmat0[m][n] = sum;
	}
      }

      /* |tmpvec1) := |tmpvec1) - |U_rl0)(U_rl0|tmpvec1) */

      for (n=0; n<EKC_core_size[Mc_AN]; n++){
	for (k=0; k<EKC_core_size[Mc_AN]; k++){
	  dum = tmpmat0[k][n];
	  for (i=0; i<Msize2[Mc_AN]; i++)  tmpvec1[n][i] -= Krylov_U_OLP[rl0][k][i]*dum;
        }
      }
    }

    /*************************************************************
                       normalization of tmpvec1
    *************************************************************/

    for (m=0; m<EKC_core_size[Mc_AN]; m++){
      for (n=0; n<EKC_core_size[Mc_AN]; n++){

        sum = 0.0;
        for (i=0; i<Msize2[Mc_AN]; i++){
          sum += tmpvec1[m][i]*tmpvec1[n][i]; 
	}

        tmpmat0[m+1][n+1] = sum;
      }
    }

    /* diagonalize tmpmat0 */ 

    if ( EKC_core_size[Mc_AN]==1){
      ko[1] = tmpmat0[1][1];
      tmpmat0[1][1] = 1.0;
    }
    else{
      Eigen_lapack( tmpmat0, ko, EKC_core_size[Mc_AN], EKC_core_size[Mc_AN] );
    }

    ZeroNum = 0;

    for (n=1; n<=EKC_core_size[Mc_AN]; n++){
      if (cutoff_value<ko[n]){
	ko[n]  = sqrt(ko[n]);
	iko[n] = 1.0/ko[n];
      }
      else{
	ZeroNum++;
	ko[n]  = 0.0;
	iko[n] = 0.0;
      }

      if (error_check==1){
        printf("rl=%3d ko=%18.15f\n",rl,ko[n]);
      }
    }

    if (error_check==1){
      printf("rl=%3d ZeroNum=%3d\n",rl,ZeroNum);
    } 

    /* tmpvec0 = tmpvec1 * tmpmat0^{-1/2} */ 

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec0[n][i] = 0.0;
      }
    }

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (k=0; k<EKC_core_size[Mc_AN]; k++){
	dum = tmpmat0[k+1][n+1]*iko[n+1];
	for (i=0; i<Msize2[Mc_AN]; i++)  tmpvec0[n][i] += tmpvec1[k][i]*dum;
      }
    }

    /*************************************************************
                         store Krylov vectors
    *************************************************************/

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
        Krylov_U_OLP[rl+1][n][i] = tmpvec0[n][i];
      }
    }

  } /* rl */

  /************************************************************
       calculate the inverse of the reduced overlap matrix
  ************************************************************/

  /*  construct the reduced overlap matrix */

  for (rl=0; rl<rlmax_EC2[Mc_AN]; rl++){

    /*  S * |Vn) */

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
	tmpvec1[n][i]  = 0.0;
      }
    }

    for (i=0; i<=can; i++){

      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[i] - 1;
      ih = S_G2M[ig];

      for (j=0; j<=can; j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[ct_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[j] - 1;

	if (0<=kl){

	  for (m=0; m<ian; m++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

	      sum = 0.0;
	      for (k=0; k<jan; k++){
		sum += OLP0[ih][kl][m][k]*Krylov_U_OLP[rl][n][Bnum+k];
	      }
	      tmpvec1[n][Anum+m] += sum;
	    }
	  } 
	}
      }
    }

    for (rl0=rl; rl0<rlmax_EC2[Mc_AN]; rl0++){
      for (m=0; m<EKC_core_size[Mc_AN]; m++){
        for (n=0; n<EKC_core_size[Mc_AN]; n++){
  	  sum = 0.0;
	  for (i=0; i<Msize2[Mc_AN]; i++){
            sum += Krylov_U_OLP[rl0][m][i]*tmpvec1[n][i]; 
	  }
          FS[rl0*EKC_core_size[Mc_AN]+m+1][rl*EKC_core_size[Mc_AN]+n+1] = sum;
          FS[rl*EKC_core_size[Mc_AN]+n+1][rl0*EKC_core_size[Mc_AN]+m+1] = sum;
	}
      }
    }
  }

  /* diagonalize FS */

  Eigen_lapack(FS,ko,Msize4[Mc_AN],Msize4[Mc_AN]);

  /* find ill-conditioned eigenvalues */

  ZeroNum = 0;
  for (i=1; i<=Msize4[Mc_AN]; i++){

    if (error_check==1){
      printf("Mc_AN=%2d i=%3d ko[i]=%18.15f\n",Mc_AN,i,ko[i]);
    }

    if (cutoff_value<ko[i]){
      iko[i] = 1.0/ko[i];
    }
    else{
      ZeroNum++;
      ko[i]  = 0.0;
      iko[i] = 0.0;
    }
  }  

  if (error_check==1){
    printf("Mc_AN=%2d ZeroNum=%2d\n",Mc_AN,ZeroNum);
  }

  /* construct the inverse */

  for (i=1; i<=Msize4[Mc_AN]; i++){
    for (j=1; j<=Msize4[Mc_AN]; j++){
      sum = 0.0;
      for (k=1; k<=Msize4[Mc_AN]; k++){
        sum += FS[i][k]*iko[k]*FS[j][k];
      }
      inv_RS[i-1][j-1] = sum;
    }
  }

  /* symmetrization of inv_RS */

  for (i=1; i<=Msize4[Mc_AN]; i++){
    for (j=i+1; j<=Msize4[Mc_AN]; j++){
      tmp0 = inv_RS[i-1][j-1];
      tmp1 = inv_RS[j-1][i-1];
      tmp2 = 0.5*(tmp0 + tmp1);
      inv_RS[i-1][j-1] = tmp2;
      inv_RS[j-1][i-1] = tmp2;
    }
  }


  





  /*
  {


    double mat1[1000][1000];
    double tsum;
    int myid;

    MPI_Comm_rank(mpi_comm_level1,&myid);

    if (myid==0){

      printf("check normalization\n");

      for (rl=0; rl<rlmax_EC2[Mc_AN]; rl++){
	for (m=0; m<EKC_core_size[Mc_AN]; m++){
	  for (rl0=0; rl0<rlmax_EC2[Mc_AN]; rl0++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

              sum = 0.0; 
   	      for (i=0; i<Msize2[Mc_AN]; i++){
                sum += Krylov_U_OLP[rl][m][i]*Krylov_U_OLP[rl0][n][i]; 
	      }
              printf("rl=%3d rl0=%3d m=%3d n=%3d  <|>=%18.15f\n",rl,rl0,m,n,sum);
	    }
	  }
	}
      }


      printf("\n\ninvS\n");


      for (rl=0; rl<rlmax_EC2[Mc_AN]; rl++){
	for (m=0; m<EKC_core_size[Mc_AN]; m++){

	  for (i=0; i<Msize2[Mc_AN]; i++){

	    sum = 0.0;

	    for (rl0=0; rl0<rlmax_EC2[Mc_AN]; rl0++){
	      for (n=0; n<EKC_core_size[Mc_AN]; n++){
		sum += inv_RS[rl*EKC_core_size[Mc_AN]+m][rl0*EKC_core_size[Mc_AN]+n]*Krylov_U_OLP[rl0][n][i]; 
	      }
	    }

	    mat1[rl*EKC_core_size[Mc_AN]+m][i] = sum; 
	  }
	}
      }


      tsum = 0.0;

      for (i=0; i<Msize2[Mc_AN]; i++){
	for (j=0; j<Msize2[Mc_AN]; j++){

	  sum = 0.0;
	  for (rl=0; rl<rlmax_EC2[Mc_AN]; rl++){
	    for (m=0; m<EKC_core_size[Mc_AN]; m++){
	      sum += Krylov_U_OLP[rl][m][i]*mat1[rl*EKC_core_size[Mc_AN]+m][j];
	    }
	  }

	  printf("i=%4d j=%4d %18.15f\n",i,j,sum);

          tsum += fabs(sum);

	}    
      }    


      printf("\n\ntsum=%18.15f\n",tsum);

    }


    MPI_Finalize();
    exit(0);
  }
  */








  /* freeing of arrays */

  for (i=0; i<(EKC_core_size[Mc_AN]+4); i++){
    free(tmpmat0[i]);
  }
  free(tmpmat0);

  for (i=0; i<(rlmax_EC2[Mc_AN]+2)*EKC_core_size[Mc_AN]; i++){
    free(FS[i]);
  }
  free(FS);

  free(ko);
  free(iko);
}












void S_orthonormalize_vec( int Mc_AN, int ct_on, double **vec,
                           double **workvec, double ****OLP0,
                           double **tmpmat0, double *ko, double *iko, 
                           int *MP, int *Msize2 )
{
  int n,i,j,can,san,fan,ct_AN;
  int k,m,ZeroNum,Anum,Bnum,ih;
  int kl,jg,jan,ig,ian;
  double dum,sum;

  ct_AN = M2G[Mc_AN];
  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;

  /* S|Vn) */

  for (n=0; n<EKC_core_size[Mc_AN]; n++){
    for (i=0; i<Msize2[Mc_AN]; i++){
      workvec[n][i]  = 0.0;
    }
  }

  for (i=0; i<=can; i++){

    ig = natn[ct_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    Anum = MP[i] - 1;
    ih = S_G2M[ig];

    for (j=0; j<=can; j++){

      kl = RMI1[Mc_AN][i][j];
      jg = natn[ct_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];
      Bnum = MP[j] - 1;

      if (0<=kl){

	for (m=0; m<ian; m++){
	  for (n=0; n<EKC_core_size[Mc_AN]; n++){

	    sum = 0.0;
	    for (k=0; k<jan; k++){
	      sum += OLP0[ih][kl][m][k]*vec[n][Bnum+k];
	    }

	    workvec[n][Anum+m] += sum;
	  }
	} 
      }
    }
  }

  /*  (Vn|S|Vn) */

  for (m=0; m<EKC_core_size[Mc_AN]; m++){
    for (n=m; n<EKC_core_size[Mc_AN]; n++){
      sum = 0.0;
      for (i=0; i<Msize2[Mc_AN]; i++){
	sum += vec[m][i]*workvec[n][i];
      }
      tmpmat0[m+1][n+1] = sum;
      tmpmat0[n+1][m+1] = sum;
    }
  }

  /* diagonalize tmpmat0 */ 

  if ( EKC_core_size[Mc_AN]==1){
    ko[1] = tmpmat0[1][1];
    tmpmat0[1][1] = 1.0;
  }
  else{
    Eigen_lapack( tmpmat0, ko, EKC_core_size[Mc_AN], EKC_core_size[Mc_AN] );
  }

  ZeroNum = 0;

  for (n=1; n<=EKC_core_size[Mc_AN]; n++){
    if (cutoff_value<ko[n]){
      ko[n]  = sqrt(ko[n]);
      iko[n] = 1.0/ko[n];
    }
    else{
      ZeroNum++;
      ko[n]  = 0.0;
      iko[n] = 0.0;
    }
  }

  /* U0 = vec * tmpmat0^{-1/2} */ 

  for (n=0; n<EKC_core_size[Mc_AN]; n++){
    for (i=0; i<Msize2[Mc_AN]; i++){
      workvec[n][i] = 0.0;
    }
  }

  for (n=0; n<EKC_core_size[Mc_AN]; n++){
    for (k=0; k<EKC_core_size[Mc_AN]; k++){
      dum = tmpmat0[k+1][n+1]*iko[n+1];
      for (i=0; i<Msize2[Mc_AN]; i++)  workvec[n][i] += vec[k][i]*dum;
    }
  }

  for (n=0; n<EKC_core_size[Mc_AN]; n++){
    for (i=0; i<Msize2[Mc_AN]; i++){
      vec[n][i] = workvec[n][i];
    }
  }

}




void Inverse_S_by_Cholesky(int Mc_AN, double ****OLP0, double **invS, int *MP, int NUM, double *LoS)
{
  char  *UPLO="U";
  INTEGER N,lda,info,lwork;
  int Gc_AN,i,j,k,Anum,Gi,wanA,time1,time2;
  int ig,ian,ih,kl,jg,jan,Bnum,m,n;
  double tmp0,tmp1;

  N = NUM;
  Gc_AN = M2G[Mc_AN];

  /* OLP0 to invS */

  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){

    ig = natn[Gc_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    Anum = MP[i] - 1;
    ih = S_G2M[ig];

    for (j=0; j<=(FNAN[Gc_AN]+SNAN[Gc_AN]); j++){

      kl = RMI1[Mc_AN][i][j];
      jg = natn[Gc_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];
      Bnum = MP[j] - 1;

      if (0<=kl){
	for (m=0; m<ian; m++){
	  for (n=0; n<jan; n++){
	    invS[Anum+m][Bnum+n] = OLP0[ih][kl][m][n];
	  }
	}
      }

      else{
	for (m=0; m<ian; m++){
	  for (n=0; n<jan; n++){
	    invS[Anum+m][Bnum+n] = 0.0;
	  }
	}
      }
    }
  }

  /* invS -> LoS */

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
       LoS[j*N+i]= invS[i][j];
    }
  }

  /* call dpotrf_() in clapack */

  lda = N;
  F77_NAME(dpotrf,DPOTRF)(UPLO, &N, LoS, &lda, &info);

  if (info!=0){
    printf("Error in dpotrf_() which is called from Embedding_Cluster.c  info=%2d\n",info);
  }

  /* call dpotri_() in clapack */

  lwork = N;
  F77_NAME(dpotri,DPOTRI)(UPLO, &N, LoS, &lda, &info);

  if (info!=0){
    printf("Error in dpotri_() which is called from Embedding_Cluster.c  info=%2d\n",info);
  }

  /* LoS -> invS */

  for (j=0; j<N; j++) {
    m = j*N; 
    for (i=0; i<=j; i++) {
      invS[i][j] = LoS[m+i];
      invS[j][i] = LoS[m+i];
    }
  }

}


void Save_DOS_Col(double ******Residues, double ****OLP0, double ***EVal, int **LO_TC, int **HO_TC)
{
  int spin,Mc_AN,wanA,Gc_AN,tno1;
  int i1,i,j,MaxL,l,h_AN,Gh_AN,wanB,tno2;
  double Stime_atom,Etime_atom; 
  double sum;
  int i_vec[10];  
  char file_eig[YOUSO10],file_ev[YOUSO10];
  FILE *fp_eig, *fp_ev;
  int numprocs,myid,ID,tag;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID){
    printf("The DOS is supported for a range from -10 to 10 eV for the O(N) Krylov subspace method.\n");
  }

  /* open file pointers */

  if (myid==Host_ID){

    sprintf(file_eig,"%s%s.Dos.val",filepath,filename);
    if ( (fp_eig=fopen(file_eig,"w"))==NULL ) {
      printf("cannot open a file %s\n",file_eig);
    }
  }
  
  sprintf(file_ev, "%s%s.Dos.vec%i",filepath,filename,myid);
  if ( (fp_ev=fopen(file_ev,"w"))==NULL ) {
    printf("cannot open a file %s\n",file_ev);
  }

  /****************************************************
                   save *.Dos.vec
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];

      fprintf(fp_ev,"<AN%dAN%d\n",Gc_AN,spin);
      fprintf(fp_ev,"%d\n",(HO_TC[spin][Mc_AN]-LO_TC[spin][Mc_AN]+1));

      for (i1=0; i1<(HO_TC[spin][Mc_AN]-LO_TC[spin][Mc_AN]+1); i1++){

	fprintf(fp_ev,"%4d  %10.6f  ",i1,EVal[spin][Mc_AN][i1-1+LO_TC[spin][Mc_AN]]);

	for (i=0; i<tno1; i++){

	  sum = 0.0;
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    Gh_AN = natn[Gc_AN][h_AN];
	    wanB = WhatSpecies[Gh_AN];
	    tno2 = Spe_Total_CNO[wanB];
	    for (j=0; j<tno2; j++){
	      sum += Residues[spin][Mc_AN][h_AN][i][j][i1]*
		               OLP0[Mc_AN][h_AN][i][j];
	    }
	  }

	  fprintf(fp_ev,"%8.5f",sum);
	}
	fprintf(fp_ev,"\n");
      }

      fprintf(fp_ev,"AN%dAN%d>\n",Gc_AN,spin);

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }
  }

  /****************************************************
                   save *.Dos.val
  ****************************************************/

  if ( (fp_eig=fopen(file_eig,"w")) != NULL ) {

    fprintf(fp_eig,"mode        5\n");
    fprintf(fp_eig,"NonCol      0\n");
    /*      fprintf(fp_eig,"N           %d\n",n); */
    fprintf(fp_eig,"Nspin       %d\n",SpinP_switch);
    fprintf(fp_eig,"Erange      %lf %lf\n",Dos_Erange[0],Dos_Erange[1]);
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

  }

  /* close file pointers */

  if (myid==Host_ID){
    if (fp_eig) fclose(fp_eig);
  }

  if (fp_ev)  fclose(fp_ev);
}

