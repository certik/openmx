/**********************************************************************
  Krylov_Dosout.c:

     Krylov_Dosout.c is a subroutine to calculate density of states by 
     a Krylov subspace method.

  Log of Krylov_Dosout.c

     10/June/2005  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
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

#ifdef TRAN
#include "tran_prototypes.h"
#endif


static void Generate_pMatrix( int spin, int Mc_AN, double *****Hks, double ****OLP0, double **invS,
                              double *****Krylov_U, double ***Krylov_U_OLP, double **inv_RS );
static void Generate_pMatrix2( int spin, int Mc_AN, double *****Hks, double ****OLP0, double *****Krylov_U ); 
static void Krylov_IOLP( int Mc_AN, double ****OLP0, double ***Krylov_U_OLP, double **inv_RS );
static void S_orthonormalize_vec( int Mc_AN, int ct_on, double **vec,
                                  double **workvec, double ****OLP0,
                                  double **tmpmat0, double *ko, double *iko );
static void Embedding_Matrix(int spin, int Mc_AN, double *****Hks,
                             double *****Krylov_U, double ****EC_matrix);
static void Inverse_S_by_Cholesky(int Mc_AN, double ****OLP0, double **invS, int *MP, int NUM);
static double Krylov_Dosout_Col(int SCF_iter,  double *****Hks, double ****OLP0 );

/**************************************************
   Msize:   \sum_FNAN Spe_Total_CNO  
   Msize2:  \sum_FNAN+SNAN Spe_Total_CNO  
   Msize3:  rlmax_EC[Mc_AN]*EKC_core_size[Mc_AN]
   Msize4:  rlmax_EC2[Mc_AN]*EKC_core_size[Mc_AN]
   Msize5:  dimension for the last column of Residues
**************************************************/

static int *MP;
static int *Msize;
static int Msize2_max;
static int *Msize2;
static int *Msize3;
static int *Msize4;
static int *Msize5;
static double **tmpvec0;
static double **tmpvec1;
static double **tmpvec2;
static double ***U0;
static double ***Krylov_U_OLP;
static double **inv_RS;

/* for Inverse_S_by_Cholesky() */
static double **invS;
static double *LoS;






double Krylov_Dosout( double *****Hks,
		      double *****ImNL,
		      double ****OLP0 )
{
  double time0;

  /****************************************************
         collinear without spin-orbit coupling
  ****************************************************/

  if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==0 ){
    time0 = Krylov_Dosout_Col(2, Hks, OLP0 );
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
    printf("Now, Krylov method is not supported for non-collinear DFT calculations.\n");
    MPI_Finalize();
    exit(1);
  }

  return time0;
}
















static double Krylov_Dosout_Col(int SCF_iter,  double *****Hks, double ****OLP0 )
{
  static int firsttime=1,recalc_firsttime=1,recalc_flag;
  int Mc_AN,Gc_AN,i,Gi,wan,wanA,wanB,Anum;
  int size1,size2,num,NUM0,NUM,NUM1,n2,Cwan,Hwan,Rn2;
  int ih,ig,ian,j,kl,jg,jan,Bnum,m,n,spin,i2;
  int k,l,i1,j1,P_min,m_size,q1,q2,csize;
  int h_AN1,Mh_AN1,h_AN2,Gh_AN1,Gh_AN2,wan1,wan2;
  int po,loopN,tno1,tno2,h_AN,Gh_AN,rl1,rl2,rl;
  int MA_AN,GA_AN,tnoA,GB_AN,tnoB,ct_on,MaxL;
  double My_TZ,TZ,sum,FermiF,time0,srt;
  double tmp1,tmp2,b2,co,x0,Sx,Dx,xmin,xmax;
  double My_Num_State,Num_State,x,Dnum;
  double emin,emax,de;
  double TStime,TEtime;
  double My_Eele0[2],My_Eele1[2];
  double max_x=30.0;
  double ChemP_MAX,ChemP_MIN,spin_degeneracy;
  double spetrum_radius;
  double **H_DC,*ko;
  double **C;
  double ***EVal;
  double ******Residues;
  double ***PDOS_DC;
  double *tmp_array;
  double *tmp_array2;
  int *Snd_H_Size,*Rcv_H_Size;
  int *Snd_S_Size,*Rcv_S_Size;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  double Stime_atom, Etime_atom;
  char file_eig[YOUSO10],file_ev[YOUSO10];
  FILE *fp_eig, *fp_ev;
  char buf1[fp_bsize];          /* setvbuf */
  char buf2[fp_bsize];          /* setvbuf */

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);

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

  Snd_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Snd_S_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_S_Size = (int*)malloc(sizeof(int)*numprocs);

  m_size = 0;
  MP = (int*)malloc(sizeof(int)*List_YOUSO[2]);

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
  PrintMemory("Embedding_Cluster: EVal",  sizeof(double)*m_size,NULL);
  }

  if (2<=level_stdout){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        printf("<Krylov> myid=%4d Mc_AN=%4d Gc_AN=%4d Msize=%4d\n",
        myid,Mc_AN,M2G[Mc_AN],Msize[Mc_AN]);
    }
  }

  /****************************************************
    allocation of arrays:

    double Residues[SpinP_switch+1]
                   [Matomnum+1]
                   [FNAN[Gc_AN]+1]
                   [Spe_Total_CNO[Gc_AN]] 
                   [Spe_Total_CNO[Gh_AN]] 
                   [NUM2]
     To reduce the memory size, the size of NUM2 is
     needed to be found in the loop.  
  ****************************************************/

  m_size = 0;
  Residues = (double******)malloc(sizeof(double*****)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    Residues[spin] = (double*****)malloc(sizeof(double****)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
        Gc_AN = 0;
        FNAN[0] = 1;
        tno1 = 1;
        n2 = 1;
      }
      else{
        wanA = WhatSpecies[Gc_AN];
        tno1 = Spe_Total_CNO[wanA];
        n2 = Msize3[Mc_AN] + 2;
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
          for (j=0; j<tno2; j++){
            Residues[spin][Mc_AN][h_AN][i][j] = (double*)malloc(sizeof(double)*n2);
	  }
        }

        m_size += tno1*tno2*n2;
      }
    }
  }

  if (firsttime){
  PrintMemory("Embedding_Cluster: Residues",sizeof(double)*m_size,NULL);
  }

  /****************************************************
    allocation of arrays:

    double PDOS[SpinP_switch+1]
               [Matomnum+1]
               [NUM]
  ****************************************************/

  m_size = 0;
  PDOS_DC = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    PDOS_DC[spin] = (double**)malloc(sizeof(double*)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      n2 = Msize3[Mc_AN] + 2;
      m_size += n2;
      PDOS_DC[spin][Mc_AN] = (double*)malloc(sizeof(double)*n2);
    }
  }

  if (firsttime)
  PrintMemory("Embedding_Cluster: PDOS_DC",sizeof(double)*m_size,NULL);

  /****************************************************
    allocation of arrays:
  ****************************************************/

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

  /****************************************************
   MPI

   Hks
  ****************************************************/

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

        Snd_H_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_H_Size[IDS] = 0;
      }

      /* receiving of size of data */

      if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){

        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_H_Size[IDR] = size2;
      }
      else{
        Rcv_H_Size[IDR] = 0;
      }

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);
    }
    else{
      Snd_H_Size[IDS] = 0;
      Rcv_H_Size[IDR] = 0;
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

        size1 = Snd_H_Size[IDS];

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

        size2 = Rcv_H_Size[IDR];
        
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

  if (SCF_iter<=1){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
        tag = 999;
 
        /* find data size to send block data */
        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
          size1 = 0;
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
                  size1++; 
                } 
              } 
	    }
          }

          Snd_S_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_S_Size[IDS] = 0;
        }

        /* receiving of size of data */
 
        if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_S_Size[IDR] = size2;
        }
        else{
          Rcv_S_Size[IDR] = 0;
        }

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);
      }
      else{
        Snd_S_Size[IDS] = 0;
        Rcv_S_Size[IDR] = 0;
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

          size1 = Snd_S_Size[IDS];

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
          
          size2 = Rcv_S_Size[IDR];
        
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

  /***********************************************
     main loop of calculation 
  ***********************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

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

      Inverse_S_by_Cholesky(Mc_AN, OLP0, invS, MP, NUM); 
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

      Krylov_IOLP( Mc_AN, OLP0, Krylov_U_OLP, inv_RS );
    }

    for (spin=0; spin<=SpinP_switch; spin++){

      /****************************************************
                generate a preconditioning matrix
      ****************************************************/

      if (SCF_iter==1 && Msize3[Mc_AN]<Msize2[Mc_AN]){
        Generate_pMatrix( spin, Mc_AN, Hks, OLP0, invS, Krylov_U, Krylov_U_OLP, inv_RS );
      }
      else if (SCF_iter==1){
        Generate_pMatrix2( spin, Mc_AN, Hks, OLP0, Krylov_U );
      }

      if (recalc_EM==1 || SCF_iter<=3 || recalc_flag==1){
        Embedding_Matrix( spin, Mc_AN, Hks, Krylov_U, EC_matrix);
      }

      /****************************************************
                  construct Hamiltonian matrix
      ****************************************************/

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

      /****************************************************
                   transform u1^+ * H_DC * u1
      ****************************************************/

      /* H_DC * u1 */

      for (i=1; i<=Msize[Mc_AN]; i++){
        for (rl=0; rl<rlmax_EC[Mc_AN]; rl++){
          for (n=0; n<EKC_core_size[Mc_AN]; n++){
            sum = 0.0;
	    for (j=1; j<=Msize[Mc_AN]; j++){
              sum += H_DC[i][j]*Krylov_U[spin][Mc_AN][rl][n][j];  
	    }
            /* transpose */
            C[rl*EKC_core_size[Mc_AN]+n+1][i] = sum;
	  }      
	}      
      }      

      /* u1^+ * H_DC * u1 */

      for (rl1=0; rl1<rlmax_EC[Mc_AN]; rl1++){
        for (m=0; m<EKC_core_size[Mc_AN]; m++){
	  for (rl2=rl1; rl2<rlmax_EC[Mc_AN]; rl2++){
	    for (n=0; n<EKC_core_size[Mc_AN]; n++){

              sum = 0.0;

              i2 = rl2*EKC_core_size[Mc_AN] + n + 1;

	      for (i=1; i<=Msize[Mc_AN]; i++){
                                                  /* transpose */
                sum += Krylov_U[spin][Mc_AN][rl1][m][i]*C[i2][i];
	      }

              H_DC[rl1*EKC_core_size[Mc_AN]+m+1][rl2*EKC_core_size[Mc_AN]+n+1] = sum;
              H_DC[rl2*EKC_core_size[Mc_AN]+n+1][rl1*EKC_core_size[Mc_AN]+m+1] = sum;
	    }
	  }
	}
      }

      /* correction for ZeroNum */

      m = (int)Krylov_U[spin][Mc_AN][0][0][0];
      for (i=1; i<=m; i++) H_DC[i][i] = 1.0e+3;

      /****************************************************
            H0 = u1^+ * H_DC * u1 + D 
      ****************************************************/

      for (i=1; i<=Msize3[Mc_AN]; i++){
        for (j=1; j<=Msize3[Mc_AN]; j++){
          H_DC[i][j] += EC_matrix[spin][Mc_AN][i][j];
        }
      }

      /****************************************************
           diagonalize
      ****************************************************/

      Eigen_lapack(H_DC,ko,Msize3[Mc_AN],Msize3[Mc_AN]);

      /********************************************
           back transformation of eigenvectors
                      c = u1 * b
      *********************************************/

      for (i=1; i<=Msize[Mc_AN]; i++){
        for (j=1; j<=Msize3[Mc_AN]; j++){
          C[i][j] = 0.0;       
	}
      }

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

      /***********************************************
          store eigenvalues and residues of poles
      ***********************************************/

      for (i=1; i<=Msize3[Mc_AN]; i++){
        EVal[spin][Mc_AN][i-1] = ko[i];
      }

      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];

      for (i=0; i<tno1; i++){
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
          Bnum = MP[h_AN];
          for (j=0; j<tno2; j++){
            for (i1=1; i1<=Msize3[Mc_AN]; i1++){
              Residues[spin][Mc_AN][h_AN][i][j][i1-1] = C[1+i][i1]*C[Bnum+j][i1];
	    }
	  }
	}
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

  /****************************************************
                   fprintf *.Dos.val
  ****************************************************/

  sprintf(file_eig,"%s%s.Dos.val",filepath,filename);

  if (myid==Host_ID){
    if ( (fp_eig=fopen(file_eig,"w")) != NULL ) {

#ifdef xt3
      setvbuf(fp_eig,buf1,_IOFBF,fp_bsize);  /* setvbuf */
#endif

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

      fclose(fp_eig);
    }
    else{
      printf("failure of saving %s\n",file_eig);
    }
  }

  /****************************************************
              calculate projected DOS
                      and 
               fprintf *.Dos.vec
  ****************************************************/

  sprintf(file_ev, "%s%s.Dos.vec%i",filepath,filename,myid);

  if ( (fp_ev=fopen(file_ev,"w")) != NULL ) {

#ifdef xt3
    setvbuf(fp_ev,buf2,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    for (spin=0; spin<=SpinP_switch; spin++){
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	dtime(&Stime_atom);

	Gc_AN = M2G[Mc_AN];
	wanA = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_CNO[wanA];

        fprintf(fp_ev,"<AN%dAN%d\n",Gc_AN,spin);
        fprintf(fp_ev,"%d\n",Msize3[Mc_AN]);

        for (i1=0; i1<Msize3[Mc_AN]; i1++){

          fprintf(fp_ev,"%4d  %10.6f  ",i1,EVal[spin][Mc_AN][i1]);

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

    fclose(fp_ev);
  }
  else {
    printf("failure of saving %s\n",file_ev);
  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  free(Snd_H_Size);
  free(Rcv_H_Size);

  free(Snd_S_Size);
  free(Rcv_S_Size);

  free(MP);
  free(Msize);
  free(Msize2);
  free(Msize3);
  free(Msize4);

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
        FNAN[0] = 1;
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

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  /* for PrintMemory */
  firsttime=0;

  return time0;
}





















void Generate_pMatrix( int spin, int Mc_AN, double *****Hks, double ****OLP0, double **invS,
                       double *****Krylov_U, double ***Krylov_U_OLP, double **inv_RS ) 
{
  int rl,rl0,rl1,ct_AN,fan,san,can,wan,ct_on,i,j;
  int n,Anum,Bnum,k,ian,ih,kl,jg,ig,jan,m,m1,n1;
  int ZeroNum,Gh_AN,wanB;
  double sum,dum,tmp0,tmp1,tmp2,tmp3,rcutA,r0;
  double **Utmp,**matRS0,**matRS1;
  double **tmpmat0;
  double *ko,*iko;
  double **FS;
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

  S_orthonormalize_vec( Mc_AN, ct_on, tmpvec0, tmpvec1, OLP0, tmpmat0, ko, iko );

  for (n=0; n<EKC_core_size[Mc_AN]; n++){
    for (i=0; i<Msize2[Mc_AN]; i++){
      U0[0][n][i] = tmpvec0[n][i];
    }
  }

  /****************************************************
           generate Krylov subspace vectors
  ****************************************************/

  for (rl=0; rl<(rlmax_EC[Mc_AN]-1); rl++){

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

    /*******************************************************
           S^{-1} * H * |Wn)
    *******************************************************/

    if (EKC_Exact_invS_flag==1){ 

      for (n=0; n<EKC_core_size[Mc_AN]; n++){
        for (i=0; i<Msize2[Mc_AN]; i++){
	  sum = 0.0;
	  for (j=0; j<Msize2[Mc_AN]; j++){
	    sum += invS[i][j]*tmpvec1[n][j];
	  }
	  tmpvec0[n][i] = sum;
        } 
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

    /*************************************************************
                   S-orthonormalization of tmpvec0
    *************************************************************/

    S_orthonormalize_vec( Mc_AN, ct_on, tmpvec0, tmpvec1, OLP0, tmpmat0, ko, iko );

    for (n=0; n<EKC_core_size[Mc_AN]; n++){
      for (i=0; i<Msize2[Mc_AN]; i++){
        U0[rl+1][n][i] = tmpvec0[n][i];
      }
    }

  } /* rl */

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








void Generate_pMatrix2( int spin, int Mc_AN, double *****Hks, double ****OLP0, double *****Krylov_U )
{
  int rl,rl0,rl1,ct_AN,fan,san,can,wan,ct_on,i,j;
  int n,Anum,Bnum,k,ian,ih,kl,jg,ig,jan,m,m1,n1;
  int ZeroNum,rl_half;
  double sum,dum,tmp0,tmp1,tmp2,tmp3;
  double **Utmp;
  double *ko,*iko;
  double **FS;

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
                      double *****Krylov_U, double ****EC_matrix)
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










void Krylov_IOLP( int Mc_AN, double ****OLP0, double ***Krylov_U_OLP, double **inv_RS )
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
                           double **tmpmat0, double *ko, double *iko )
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




void Inverse_S_by_Cholesky(int Mc_AN, double ****OLP0, double **invS, int *MP, int NUM)
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


  /*
  printf("S\n");

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      printf("%12.8f ",invS[i][j]);
    }
    printf("\n");
  }
  */

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


  /*
  {
    int myid;
    double tsum;

    MPI_Comm_rank(mpi_comm_level1,&myid);
  
    if (myid==0){

      printf("invS\n");

      tsum = 0.0;

      for (i=0;i<N;i++) {
	for (j=0;j<N;j++) {
	  printf("i=%4d j=%4d %18.15f\n",i,j,invS[i][j]);

          tsum += fabs(invS[i][j]);

	}
      }

      printf("\n\ntsum=%18.15f\n",tsum);

    }


    MPI_Finalize();
    exit(0);

  }
  */

}

