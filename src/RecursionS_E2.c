/**********************************************************************
  RecursionS.c:

     RecursionS.c is a subroutine to solve a generalized eigenvalue
     problem with an overlap matrix for both spin and non-spin
     polarized systems using a generalized recursion method in O(N) operation.

  Log of RecursionS.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#define  measure_time   1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef TRAN
#include "tran_prototypes.h"
#endif





static void Inverse_S_by_LU(int Mc_AN, double ****OLP0, double ***S12);
static void TSB_Lanczos(int Mc_AN, int spin, double **HP, double **HR);
static void GBTD00(int T_switch, int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S,
                   dcomplex **cinv, double *be2, double *cterR, double *cterI,
                   double **cinvR, double **cinvI, double **ctp1R, double **ctp1I,
                   double **ctp2R, double **ctp2I,
                   INTEGER *ipiv, dcomplex *a1d, dcomplex *work);
static void RecurG5(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
                    INTEGER *ipiv, dcomplex *a1d, dcomplex *work,
                    dcomplex **Ctp, dcomplex **Ctp1, dcomplex **Ctp2, dcomplex **Ctp3 );
static void T_Conserve_MP(int T_switch);
static void Calc_DM_with_S(int T_switch, double *****CDM, double *****EDM, double ****IOLP);
static void Inverse_ComplexMat(int n, dcomplex **a, INTEGER *ipiv,
			       dcomplex *a1d, dcomplex *work);
static int SVD_fact_inv(int rl, int n, double **a,
                        double **B, double **IB,
                        double **C, double **IC);










static void RecurG(int ct_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S);
static void RecurG2(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S);
static void RecurG3(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
                    INTEGER *ipiv, dcomplex *a1d, dcomplex *work);
static void RecurG4(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
                    INTEGER *ipiv, dcomplex *a1d, dcomplex *work);
static void RecurG6(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
                    INTEGER *ipiv, dcomplex *a1d, dcomplex *work);







static void Save_Recursion();
static void Output_RcnCof(FILE *fp);
static void Output_LU(FILE *fp);

static int **RCN_Rank;
static double Uele_OS[4];
static double Uele_IS[4];

static double *****al;
static double *****be;
static double *****ibe;
static double *****ga;
static double *****tga;
static double *****iga;

static double **HR;
static double **HP;
static int *Msize,*MP;

static double **BeGa;
static double **UmRn;
static double **RnUm;

static double **TmpCoes;
static double **RCpre;
static double **RCcnt;
static double **RRcnt;
static double **LCpre;
static double **LCcnt;
static double **LRcnt;
static double ***RU;
static double *****LU;

/* for Inverse_S_by_LU() */
static INTEGER *piv;
static double *LoS,*work;

static int Msize_max;






double RecursionS_E(int SCF_iter,
                    double *****Hks,
                    double ****OLP0,
                    double *****CDM,
                    double *****EDM,
                    double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  static int i,spin;
  static double time0;
  static int m,n,j,Rn,k,l,rl;
  static int ig,ian,ih,kl,jg,jan;
  static int Mc_AN,Gc_AN,tno2,num; 
  static int size1,size2,n2,Gi,wanA,NUM; 
  static int Anum,Bnum;
  static double sum;
  static double TStime,TEtime;
  static int size_HPS;
  static int tno0,tno1,Cwan,h_AN; 
  static int Gh_AN,Hwan,GA_AN,GB_AN,tnoA,tnoB,wanB;
  static double *tmp_array;
  static double *tmp_array2;
  static int *Snd_H_Size,*Rcv_H_Size;
  static int *Snd_S_Size,*Rcv_S_Size;
  static int numprocs,myid,ID,tag;
  static double Stime_proc, Etime_proc;

  MPI_Status stat;
  MPI_Request request; 

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);
  Timetool_times("RecursionS_B","start");

  /****************************************************
    allocation of arrays:
  ****************************************************/

  /**********************
    find the matrix size
  ***********************/

  Msize = (int*)malloc(sizeof(int)*(Matomnum+1));

  Msize_max = 0;

  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (Mc_AN==0){
      Gc_AN = 0;
      FNAN[0] = 1;
      SNAN[0] = 0;
      n2 = 1;
      Msize[Mc_AN] = 1;
    }
    else{
      Anum = 1;
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	Gi = natn[Gc_AN][i];
	wanA = WhatSpecies[Gi];
	Anum += Spe_Total_CNO[wanA];
      }
      NUM = Anum - 1;
      Msize[Mc_AN] = NUM;
    }
 
    if (Msize_max<Msize[Mc_AN]) Msize_max = Msize[Mc_AN] + 4;
  }

  al = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    al[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      al[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        al[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          al[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
          for (m=0; m<List_YOUSO[7]; m++) al[i][j][k][l][m] = 0.0;
	}
      }
    }
  }

  be = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    be[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      be[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        be[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          be[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  ibe = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    ibe[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      ibe[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        ibe[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          ibe[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  ga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    ga[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      ga[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        ga[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          ga[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  tga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    tga[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      tga[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        tga[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          tga[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  iga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    iga[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      iga[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        iga[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          iga[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  LU = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    LU[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    for (j=0; j<=Matomnum; j++){
      LU[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        LU[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          LU[i][j][k][l] = (double*)malloc(sizeof(double)*Msize_max); 
          for (m=0; m<Msize_max; m++) LU[i][j][k][l][m] = 0.0; 
	}
      }
    }
  }

  RCN_Rank = (int**)malloc(sizeof(int*)*(Matomnum+1)); 
  for (i=0; i<=Matomnum; i++){
    RCN_Rank[i] = (int*)malloc(sizeof(int)*List_YOUSO[3]); 
  }

  Snd_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Snd_S_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_S_Size = (int*)malloc(sizeof(int)*numprocs);

  MP = (int*)malloc(sizeof(int)*List_YOUSO[2]);

  /****************************************************
     allocation of arrays
  ****************************************************/

  HR = (double**)malloc(sizeof(double*)*Msize_max); 
  for (i=0; i<Msize_max; i++){
    HR[i] = (double*)malloc(sizeof(double)*Msize_max); 
  }

  HP = (double**)malloc(sizeof(double*)*Msize_max); 
  for (i=0; i<Msize_max; i++){
    HP[i] = (double*)malloc(sizeof(double)*Msize_max); 
  }

  BeGa = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    BeGa[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  UmRn = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    UmRn[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  RnUm = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    RnUm[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  TmpCoes = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    TmpCoes[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  RCpre = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    RCpre[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  RCcnt = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    RCcnt[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  RRcnt = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    RRcnt[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  LCpre = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    LCpre[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  LCcnt = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    LCcnt[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  LRcnt = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    LRcnt[i] = (double*)malloc(sizeof(double)*Msize_max);
  }
  
  RU = (double***)malloc(sizeof(double**)*List_YOUSO[3]);
  for (i=0; i<List_YOUSO[3]; i++){
    RU[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      RU[i][j] = (double*)malloc(sizeof(double)*Msize_max);
    }
  }  

  piv = (INTEGER*)malloc(sizeof(INTEGER)*Msize_max);
  LoS = (double*)malloc(sizeof(double)*Msize_max*Msize_max);
  work = (double*)malloc(sizeof(double)*Msize_max);

  /****************************************************
   MPI

   Hks
  ****************************************************/

  /***********************************
             set data size
  ************************************/

  for (ID=0; ID<numprocs; ID++){

    if (ID!=myid){
      tag = 999;

      /* find data size to send block data */
      if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0){

        size1 = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[ID]+S_Snd_Num[ID]); n++){
            Mc_AN = Snd_MAN[ID][n];
            Gc_AN = Snd_GAN[ID][n];
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

        Snd_H_Size[ID] = size1;
        MPI_Isend(&size1, 1, MPI_INT, ID, tag, mpi_comm_level1, &request);

      }
      else{
        Snd_H_Size[ID] = 0;
      }

      /* receiving of size of data */

      if ((F_Rcv_Num[ID]+S_Rcv_Num[ID])!=0){


        MPI_Recv(&size2, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
        Rcv_H_Size[ID] = size2;
      }
      else{
        Rcv_H_Size[ID] = 0;
      }

      if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0) MPI_Wait(&request,&stat);
    }
    else{
      Snd_H_Size[ID] = 0;
      Rcv_H_Size[ID] = 0;
    }
  }

  /***********************************
             data transfer
  ************************************/

  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    if (ID!=myid){

      /*****************************
              sending of data 
      *****************************/

      if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0){

        size1 = Snd_H_Size[ID];

        /* allocation of array */

        tmp_array = (double*)malloc(sizeof(double)*size1);

        /* multidimentional array to vector array */

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[ID]+S_Snd_Num[ID]); n++){
            Mc_AN = Snd_MAN[ID][n];
            Gc_AN = Snd_GAN[ID][n];
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

        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, ID, tag, mpi_comm_level1, &request);

      }

      /*****************************
         receiving of block data
      *****************************/

      if ((F_Rcv_Num[ID]+S_Rcv_Num[ID])!=0){

        size2 = Rcv_H_Size[ID];
        
        /* allocation of array */
        tmp_array2 = (double*)malloc(sizeof(double)*size2);
        
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          Mc_AN = S_TopMAN[ID] - 1;  /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[ID]+S_Rcv_Num[ID]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[ID][n];
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

      if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0){
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

  if (SCF_iter<=2){

    for (ID=0; ID<numprocs; ID++){

      if (ID!=myid){
        tag = 999;

        /* find data size to send block data */
        if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0){
          size1 = 0;
          for (n=0; n<(F_Snd_Num[ID]+S_Snd_Num[ID]); n++){
            Mc_AN = Snd_MAN[ID][n];
            Gc_AN = Snd_GAN[ID][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              size1 += tno1*tno2;
	    }
          }

          Snd_S_Size[ID] = size1;
          MPI_Isend(&size1, 1, MPI_INT, ID, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_S_Size[ID] = 0;
        }

        /* receiving of size of data */
 
        if ((F_Rcv_Num[ID]+S_Rcv_Num[ID])!=0){
          MPI_Recv(&size2, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
          Rcv_S_Size[ID] = size2;
        }
        else{
          Rcv_S_Size[ID] = 0;
        }

        if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0) MPI_Wait(&request,&stat);
      }
      else{
        Snd_S_Size[ID] = 0;
        Rcv_S_Size[ID] = 0;
      }
    }

    /***********************************
               data transfer
    ************************************/

    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      if (ID!=myid){

        /*****************************
                sending of data 
        *****************************/

        if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0){

          size1 = Snd_S_Size[ID];

          /* allocation of array */

          tmp_array = (double*)malloc(sizeof(double)*size1);

          /* multidimentional array to vector array */

          num = 0;

          for (n=0; n<(F_Snd_Num[ID]+S_Snd_Num[ID]); n++){
            Mc_AN = Snd_MAN[ID][n];
            Gc_AN = Snd_GAN[ID][n];
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

          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, ID, tag, mpi_comm_level1, &request);
        }

        /*****************************
           receiving of block data
        *****************************/

        if ((F_Rcv_Num[ID]+S_Rcv_Num[ID])!=0){
          
          size2 = Rcv_S_Size[ID];
        
         /* allocation of array */
          tmp_array2 = (double*)malloc(sizeof(double)*size2);
         
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);

          num = 0;
          Mc_AN = S_TopMAN[ID] - 1; /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[ID]+S_Rcv_Num[ID]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[ID][n];
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

        if ((F_Snd_Num[ID]+S_Snd_Num[ID])!=0){
          MPI_Wait(&request,&stat);
          free(tmp_array); /* freeing of array */
	}
      }
    }
  }

  /****************************************************
     calculate the inverse matrix of overlap matrix
  ****************************************************/

  if (measure_time==1) dtime(&Stime_proc);

  if (SCF_iter==1){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Inverse_S_by_LU(Mc_AN, OLP0, S12); 
    }
  }

  if (measure_time==1){
    dtime(&Etime_proc);
    printf("myid=%2d  Time Inverse_S_by_LU = %15.12f\n",myid,Etime_proc-Stime_proc);
  }

  /****************************************************
          two-sided block Lanczos transformation
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    /***********************************************
     MP array

     Note:  MP indicates the starting position of
            atom i in arraies H and S
    ***********************************************/
    
    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      MP[i] = Anum;
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum += Spe_Total_CNO[wanA];
    }

    for (spin=0; spin<=SpinP_switch; spin++){

      /****************************************************
                          construct HR
      ****************************************************/

      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	ig = natn[Gc_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];
	Anum = MP[i];
	ih = S_G2M[ig];

	for (j=0; j<=(FNAN[Gc_AN]+SNAN[Gc_AN]); j++){

	  kl = RMI1[Mc_AN][i][j];
	  jg = natn[Gc_AN][j];
	  jan = Spe_Total_CNO[WhatSpecies[jg]];
	  Bnum = MP[j];

	  if (0<=kl){
	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
		HR[Anum+m][Bnum+n] = Hks[spin][ih][kl][m][n];
	      }
	    }
	  }

	  else{
	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
		HR[Anum+m][Bnum+n] = 0.0;
	      }
	    }
	  }
	}
      }
               
      /****************************************************
                       calculate S^{-1}H
      ****************************************************/

      if (measure_time==1) dtime(&Stime_proc);

      for (i=1; i<=Msize[Mc_AN]; i++){
        for (j=1; j<=Msize[Mc_AN]; j++){
          sum = 0.0;
          for (k=1; k<=Msize[Mc_AN]; k++){
            /* note: HR is real and symmetric */
            sum += S12[Mc_AN][i][k]*HR[j][k]; 
	  }
          HP[i][j] = sum;
	}
      }

      if (measure_time==1){
        dtime(&Etime_proc);
        printf("myid=%2d  Time S^{-1}H         = %15.12f\n",myid,Etime_proc-Stime_proc);
      }

      /****************************************************
                       copy HP^t to HR 
      ****************************************************/

      for (i=1; i<=Msize[Mc_AN]; i++){
        for (j=1; j<=Msize[Mc_AN]; j++){
          HR[i][j] = HP[j][i];           
	}
      }

      /****************************************************
             two-sided block Lanczos transformation
      ****************************************************/

      if (measure_time==1) dtime(&Stime_proc);

      TSB_Lanczos(Mc_AN,spin,HP,HR);

      if (measure_time==1){
        dtime(&Etime_proc);
        printf("myid=%2d  Time TSB_Lanczos     = %15.12f\n",myid,Etime_proc-Stime_proc);
      }

    }
  }

  /****************************************************
        search the chemical potential to maintain
                the number of electrons
  ****************************************************/

  if (measure_time==1) dtime(&Stime_proc);
  
  T_Conserve_MP(T_switch);

  if (measure_time==1){
    dtime(&Etime_proc);
    printf("myid=%2d  Time T_Conserve_MP   = %15.12f\n",myid,Etime_proc-Stime_proc);
  }
  
  /****************************************************
             Calculate the density matrix
  ****************************************************/

  if (measure_time==1) dtime(&Stime_proc);

  Calc_DM_with_S(T_switch,CDM,EDM,IOLP);

  if (measure_time==1){
    dtime(&Etime_proc);
    printf("myid=%2d  Time Calc_DM_with_S  = %15.12f\n",myid,Etime_proc-Stime_proc);
  }

  if (SpinP_switch==0){
    Eele0[0] = Uele_OS[0];
    Eele0[1] = Uele_OS[0];
    Eele1[0] = Uele_IS[0];   
    Eele1[1] = Uele_IS[0];
  }
  else if (SpinP_switch==1){
    Eele0[0] = Uele_OS[0];
    Eele0[1] = Uele_OS[1];
    Eele1[0] = Uele_IS[0];   
    Eele1[1] = Uele_IS[1];
  }
  
  /****************************************************
   Output several informations of recursion algorithm
  ****************************************************/

  /*
  if (2<=level_fileout) Save_Recursion();
  */

















  /****************************************************
                     Bond Energies
  ****************************************************/

  /*
  sum = 0.0;
  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      GA_AN = M2G[Mc_AN];    
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (j=0; j<=FNAN[GA_AN]; j++){
	GB_AN = natn[GA_AN][j];  
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];

	for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
	    sum += CDM[spin][Mc_AN][j][k][l]*Hks[spin][Mc_AN][j][k][l];
	  }
	}
      }
    }
  }



  printf("Bond Energy=%15.12f\n",sum);

  sum = 0.0;
  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      GA_AN = M2G[Mc_AN];    
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];

      for (j=0; j<=FNAN[GA_AN]; j++){
	GB_AN = natn[GA_AN][j];  
	wanB = WhatSpecies[GB_AN];
	tnoB = Spe_Total_CNO[wanB];

	for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
	    sum += fabs(CDM[spin][Mc_AN][j][k][l]);
	  }
	}
      }
    }
  }


  printf("sum of CDM=%15.12f\n",sum);


  MPI_Finalize();
  exit(0);
  */













  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(al[i][j][k][l]);
	}
        free(al[i][j][k]);
      }
      free(al[i][j]);
    }
    free(al[i]);
  }
  free(al);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(be[i][j][k][l]);
	}
        free(be[i][j][k]);
      }
      free(be[i][j]);
    }
    free(be[i]);
  }
  free(be);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(ibe[i][j][k][l]);
	}
        free(ibe[i][j][k]);
      }
      free(ibe[i][j]);
    }
    free(ibe[i]);
  }
  free(ibe);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(ga[i][j][k][l]);
	}
        free(ga[i][j][k]);
      }
      free(ga[i][j]);
    }
    free(ga[i]);
  }
  free(ga);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(tga[i][j][k][l]);
	}
        free(tga[i][j][k]);
      }
      free(tga[i][j]);
    }
    free(tga[i]);
  }
  free(tga);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(iga[i][j][k][l]);
	}
        free(iga[i][j][k]);
      }
      free(iga[i][j]);
    }
    free(iga[i]);
  }
  free(iga);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=Matomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(LU[i][j][k][l]);
	}
        free(LU[i][j][k]);
      }
      free(LU[i][j]);
    }
    free(LU[i]);
  }
  free(LU);

  for (i=0; i<=Matomnum; i++){
    free(RCN_Rank[i]);
  }
  free(RCN_Rank);

  free(Snd_H_Size);
  free(Rcv_H_Size);
  free(Snd_S_Size);
  free(Rcv_S_Size);

  free(Msize);
  free(MP);

  for (i=0; i<Msize_max; i++){
    free(HR[i]);
  }
  free(HR);

  for (i=0; i<Msize_max; i++){
    free(HP[i]);
  }
  free(HP);

  for (i=0; i<List_YOUSO[7]; i++){
    free(BeGa[i]);
  }
  free(BeGa);

  for (i=0; i<List_YOUSO[7]; i++){
    free(UmRn[i]);
  }
  free(UmRn);

  for (i=0; i<List_YOUSO[7]; i++){
    free(RnUm[i]);
  }
  free(RnUm);

  for (i=0; i<List_YOUSO[7]; i++){
    free(TmpCoes[i]);
  }
  free(TmpCoes);

  for (i=0; i<List_YOUSO[7]; i++){
    free(RCpre[i]);
  }
  free(RCpre);

  for (i=0; i<List_YOUSO[7]; i++){
    free(RCcnt[i]);
  }
  free(RCcnt);

  for (i=0; i<List_YOUSO[7]; i++){
    free(RRcnt[i]);
  }
  free(RRcnt);

  for (i=0; i<List_YOUSO[7]; i++){
    free(LCpre[i]);
  }
  free(LCpre);

  for (i=0; i<List_YOUSO[7]; i++){
    free(LCcnt[i]);
  }
  free(LCcnt);

  for (i=0; i<List_YOUSO[7]; i++){
    free(LRcnt[i]);
  }
  free(LRcnt);
  
  for (i=0; i<List_YOUSO[3]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RU[i][j]);
    }
    free(RU[i]);
  }  
  free(RU);

  free(piv);
  free(LoS);
  free(work);

  /* for elapsed time */
  
  dtime(&TEtime);
  Timetool_times("RecursionS_B","end");
  time0 = TEtime - TStime;
  return time0;
}






















void TSB_Lanczos(int Mc_AN, int spin, double **HP, double **HR)
{
  static int firsttime=1;
  static int i,j,k,l,m,n,ie,i1,j1,po,ig,jg,Rni,Rnj,Rn;
  static int kl,l1,l2,l3,ian,jan,ct_AN;
  static int m1,m2,m3,ZeroNum,Rloop_END;
  static int fan,san,can,rl,rl0,wan,ct_on;
  static double sum,sum1,sum2,dum,dum1,sumx,sumy,sumz,xn,xa;

  ct_AN = M2G[Mc_AN];
  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];

  for (i=0; i<List_YOUSO[7]; i++){
    for (j=0; j<=Msize[Mc_AN]; j++){
      RCpre[i][j] = 0.0;
      RCcnt[i][j] = 0.0;
      RRcnt[i][j] = 0.0;
      LCpre[i][j] = 0.0;
      LCcnt[i][j] = 0.0;
      LRcnt[i][j] = 0.0;
    }
  }

  for (i=0; i<List_YOUSO[3]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      for (k=0; k<Msize_max; k++){
        RU[i][j][k] = 0.0;
      }  
    }
  }  

  for (i=0; i<ct_on; i++){
    RCcnt[i][i+1] = 1.0; 
    LCcnt[i][i+1] = 1.0;
    LU[spin][Mc_AN][0][i][i+1] = 1.0;
    RU[0][i][i+1] = 1.0;
  }

  for (k=0; k<List_YOUSO[3]; k++){
    for (i=0; i<List_YOUSO[7]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        be[spin][Mc_AN][k][i][j] = 0.0;
        ga[spin][Mc_AN][k][i][j] = 0.0;
      }
    }
  }

  /****************************************************
     Recursion of two-sided block Lanczos algorithm
  ****************************************************/

  Rloop_END = 0;
  RNUM[ct_AN] = rlmax;
  RCN_Rank[Mc_AN][0] = ct_on;

  /****************************************************
    H_TD = A0  B1  
           C1  A1  B1
               C1  A2  B2  
                   C2  A3  B3
                       C3  A4 ..
                       .......
    A -> al
    B -> be
    C -> ga
  ****************************************************/


  rl = 0;

  do {

    /****************************************************
      H'|WRn}
    ****************************************************/

    for (n=0; n<ct_on; n++){
      for (i=1; i<=Msize[Mc_AN]; i++){
	sum = 0.0;
	for (j=1; j<=Msize[Mc_AN]; j++){
	  sum += HP[i][j]*RCcnt[n][j];
	}
	RRcnt[n][i] = sum;
      } 
    }

    /****************************************************
      {WLn|H'
    ****************************************************/

    for (n=0; n<ct_on; n++){
      for (i=1; i<=Msize[Mc_AN]; i++){
	sum = 0.0;
	for (j=1; j<=Msize[Mc_AN]; j++){
	  sum += LCcnt[n][j]*HR[i][j];
	}
	LRcnt[n][i] = sum;
      } 
    }

    /****************************************************
       Alpha_n = {WLn|H'|WRn}

                  Alpha_n is calculated as
                {WLn|H'|WRn} = <LCcnt|RRcnt>.
    ****************************************************/

    for (m=0; m<ct_on; m++){
      for (n=0; n<ct_on; n++){
        sum = 0.0; 
        for (i=1; i<=Msize[Mc_AN]; i++){
          sum += LCcnt[m][i]*RRcnt[n][i];
	}
        al[spin][Mc_AN][rl][m][n] = sum;
      }
    }

    if (rl!=RNUM[ct_AN]){

      /****************************************************
             |RRcnt} = H'|WRn} - |WRn-1}Bn - |WRn}An
      ****************************************************/

      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++)      TmpCoes[n][i]  = 0.0;
      }

      for (n=0; n<ct_on; n++){

        if (1<=rl){ 
          for (k=0; k<ct_on; k++){
            dum = be[spin][Mc_AN][rl][k][n];
            for (i=1; i<=Msize[Mc_AN]; i++)  TmpCoes[n][i] += RCpre[k][i]*dum;
	  }
	}

	for (k=0; k<ct_on; k++){
          dum = al[spin][Mc_AN][rl][k][n];
          for (i=1; i<=Msize[Mc_AN]; i++)    TmpCoes[n][i] += RCcnt[k][i]*dum;
	}

        for (i=1; i<=Msize[Mc_AN]; i++)      RRcnt[n][i]   -= TmpCoes[n][i]; 
      }         

      /****************************************************
             {LRcnt| = {WLn|H' - Cn{WLn-1| - An{WLn|
      ****************************************************/

      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++)      TmpCoes[n][i]  = 0.0;
      }

      for (n=0; n<ct_on; n++){

        if (1<=rl){ 
	  for (k=0; k<ct_on; k++){
            dum = ga[spin][Mc_AN][rl][n][k];
            for (i=1; i<=Msize[Mc_AN]; i++)  TmpCoes[n][i] += dum*LCpre[k][i];
	  }
	}

	for (k=0; k<ct_on; k++){
          dum = al[spin][Mc_AN][rl][n][k];
          for (i=1; i<=Msize[Mc_AN]; i++)    TmpCoes[n][i] += dum*LCcnt[k][i];
	}

        for (i=1; i<=Msize[Mc_AN]; i++)      LRcnt[n][i]   -= TmpCoes[n][i]; 
      }

      /****************************************************
                       Rebiorthogonalization
       by the two-sided modified block Gram-Schmidt method
      ****************************************************/

      /* |RRcnt) := (I - \sum_{rl0} |RU_rl0)(LU_rl0|)|RRcnt) */

      for (rl0=0; rl0<=rl; rl0++){

        /* (LU_rl0|RRcnt) */

	for (m=0; m<ct_on; m++){
	  for (n=0; n<ct_on; n++){

	    sum = 0.0; 
	    for (i=1; i<=Msize[Mc_AN]; i++){
	      sum += LU[spin][Mc_AN][rl0][m][i]*RRcnt[n][i];
	    }

	    UmRn[m][n] = sum;
	  }
	}

        /* |RRcnt) := |RRcnt) - |RU_rl0)(LU_rl0|RRcnt) */

        for (n=0; n<ct_on; n++){
          for (i=1; i<=Msize[Mc_AN]; i++)    TmpCoes[n][i]  = 0.0;
        }

        for (n=0; n<ct_on; n++){
          for (k=0; k<ct_on; k++){
            dum = UmRn[k][n];
            for (i=1; i<=Msize[Mc_AN]; i++)  TmpCoes[n][i] += RU[rl0][k][i]*dum;
	  }
          for (i=1; i<=Msize[Mc_AN]; i++)    RRcnt[n][i]   -= TmpCoes[n][i]; 
  	}
      }

      /* (LRcnt| := (LRcnt| (I - \sum_{rl0}) |RU_rl0)(LU_rl0|) */

      for (rl0=0; rl0<=rl; rl0++){

        /* (LRcnt|RU_rl0) */

        for (m=0; m<ct_on; m++){
	  for (n=0; n<ct_on; n++){

	    sum = 0.0; 
	    for (i=1; i<=Msize[Mc_AN]; i++){
	      sum += LRcnt[m][i]*RU[rl0][n][i];
	    }

            RnUm[m][n]  = sum;
	  }
	}

        /* (LRcnt| := (LRcnt| - (LRcnt|RU_rl0)(LU_rl0| */

        for (n=0; n<ct_on; n++){
          for (i=1; i<=Msize[Mc_AN]; i++)    TmpCoes[n][i]  = 0.0;
        }

        for (n=0; n<ct_on; n++){
          for (k=0; k<ct_on; k++){
            dum = RnUm[n][k];
            for (i=1; i<=Msize[Mc_AN]; i++)  TmpCoes[n][i] += dum*LU[spin][Mc_AN][rl0][k][i];
	  }
          for (i=1; i<=Msize[Mc_AN]; i++)    LRcnt[n][i]   -= TmpCoes[n][i]; 
	}
      }

      /****************************************************
                    Bn+1 * Cn+1 = {LRcnt|RRcnt}
      ****************************************************/

      for (m=0; m<ct_on; m++){
	for (n=0; n<ct_on; n++){

	  sum = 0.0; 
	  for (i=1; i<=Msize[Mc_AN]; i++)  sum += LRcnt[m][i]*RRcnt[n][i];

	  BeGa[m][n] = sum;
	}
      }

      /****************************************************
	1) Perform the singular value decomposition of 
           {LRcnt|RRcnt} as U*W*V^t.
        2) Then, B_{n+1} = UW^{1/2} and C_{n+1}=W^{1/2}V^t.
        3) It is easy to calculate the inverses of B and C.
      ****************************************************/

      printf("BeGa\n");
      for (m=0; m<ct_on; m++){
	for (n=0; n<ct_on; n++){
          printf("%8.4f ",BeGa[m][n]);
	}
        printf("\n");
      }

      /* singular value decomposition */

      RCN_Rank[Mc_AN][rl+1] = SVD_fact_inv(
                                rl,ct_on-1,BeGa,
                                be[spin][Mc_AN][rl+1],ibe[spin][Mc_AN][rl+1],
                                ga[spin][Mc_AN][rl+1],iga[spin][Mc_AN][rl+1]);

      /* The following termination condition couble be modified */

      if (RCN_Rank[Mc_AN][rl+1]==0){
        RNUM[ct_AN] = rl;
        Rloop_END = 1;
      }
      else if (RCN_Rank[Mc_AN][rl+1]==1 && ct_on!=1 && Rloop_END!=-1){
        Rloop_END = -1;
      }
      else if (Rloop_END==-1){
        RNUM[ct_AN] = rl;
        Rloop_END = 1;
      }

      /*
      printf("RCN_Rank[Mc_AN][rl+1]=%i  Rloop_END=%i  ct_on=%i\n",
              RCN_Rank[Mc_AN][rl+1],Rloop_END,ct_on);
      */

      /****************************************************
         store the transposed ga into tga
      ****************************************************/

      for (m=0; m<ct_on; m++){
        for (n=0; n<ct_on; n++){
          tga[spin][Mc_AN][rl+1][n][m] = ga[spin][Mc_AN][rl+1][m][n];
	}
      }

      /****************************************************
                     |WRn+1} = |RRcnt}ICn+1
      ****************************************************/

      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++){
          RCpre[n][i] = RCcnt[n][i];
	}
      }

      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++) TmpCoes[n][i]  = 0.0;
      }

      for (n=0; n<ct_on; n++){
        for (k=0; k<ct_on; k++){
          dum = iga[spin][Mc_AN][rl+1][k][n];
          for (i=1; i<=Msize[Mc_AN]; i++){
            TmpCoes[n][i] += RRcnt[k][i]*dum;
	  }
	}

        for (i=1; i<=Msize[Mc_AN]; i++){
	  RCcnt[n][i]    = TmpCoes[n][i];
	  RU[rl+1][n][i] = TmpCoes[n][i];   
	}
      } 

      /****************************************************
                      {WLn+1| = IBn+1{LRcnt|
      ****************************************************/

      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++){
          LCpre[n][i] = LCcnt[n][i];
	}
      }

      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++) TmpCoes[n][i]  = 0.0;
      }

      for (n=0; n<ct_on; n++){
        for (k=0; k<ct_on; k++){
          dum = ibe[spin][Mc_AN][rl+1][n][k];
          for (i=1; i<=Msize[Mc_AN]; i++){
            TmpCoes[n][i] += dum*LRcnt[k][i];
	  }
	}
        for (i=1; i<=Msize[Mc_AN]; i++){
	  LCcnt[n][i]                 = TmpCoes[n][i];
	  LU[spin][Mc_AN][rl+1][n][i] = TmpCoes[n][i];   
	}
      }	

    } /* end of if (rl!=RNUM[ct_AN]) */

    rl++;

    /****************************************************
                       End of rl loop
    ****************************************************/

  } while (rl<=RNUM[ct_AN] && Rloop_END<=0);

}











int SVD_fact_inv(int rl, int n, double **a,
                 double **B, double **IB,
                 double **C, double **IC)
{

  static int i,j,j1,k,rankN,po,num;
  static double sum,rank_criterion,rank_criterion2,scale;
  static double *a1,*vt,*u,*w,*work,*w2,*iw2;
  static double **vt1,**u1;
  static char jobu  = 'A';
  static char jobvt = 'A';
  static INTEGER ROW,COL,lda,ldu,ldvt,lwork,info;

  /****************************************************
    allocation of arrays:

    static double a1[YOUSO7*YOUSO7];
    static double vt[YOUSO7*YOUSO7];
    static double u[YOUSO7*YOUSO7];
    static double work[YOUSO7*YOUSO7];
    static double w[YOUSO7+2];
    static double w2[YOUSO7+2];
    static double iw2[YOUSO7+2];
    static double vt1[YOUSO7][YOUSO7];
    static double u1[YOUSO7][YOUSO7];
  ****************************************************/

  /* a1 */
  a1 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* vt */
  vt = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* u */
  u = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* work */
  work = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* w */
  w = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 

  /* w2 */
  w2 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 

  /* iw2 */
  iw2 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 

  /* vt1 */
  vt1 = (double**)malloc(sizeof(double*)*(List_YOUSO[7]+2)); 
  for (i=0; i<(List_YOUSO[7]+2); i++){
    vt1[i] = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 
  }

  /* u1 */
  u1 = (double**)malloc(sizeof(double*)*(List_YOUSO[7]+2)); 
  for (i=0; i<(List_YOUSO[7]+2); i++){
    u1[i] = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 
  }

  /****************************************************
                      From 0 to n
  ****************************************************/

  num = -1;  
  for (j=0; j<=n; j++){
    for (i=0; i<=n; i++){
      num++;
      a1[num] = a[i][j];
    }
  }

  /****************************************************
             singular value decomposition
  ****************************************************/

  ROW  = n + 1;
  COL  = n + 1;
  lda  = n + 1;
  ldu  = n + 1;
  ldvt = n + 1;
  lwork = 5*List_YOUSO[7];
  F77_NAME(dgesvd,DGESVD)(&jobu, &jobvt, &ROW, &COL, a1, &lda, w, u,
          &ldu, vt, &ldvt, work, &lwork, &info);

  num = -1;  
  for (j=0; j<=n; j++){
    for (i=0; i<=n; i++){
      num++;
      u1[i][j]  =  u[num];
      vt1[i][j] = vt[num];
    }
  }

  /****************************************************
           calculate rank of the BeGa matrix
  ****************************************************/

  rank_criterion  = 1.0e-11;
  rank_criterion2 = 1.0e-12;

  rankN = -1;
  po = 0;
  do{
    rankN++;
    if (  (fabs(w[rankN]/w[0])<rank_criterion && fabs(w[rankN])<100.0*rank_criterion)
        ||
           fabs(w[rankN])<rank_criterion2 && rankN==0) po = 1; 
  } while(rankN<n && po==0);
  if (po==0) rankN++;


  printf("rankN=%2d\n",rankN);
  for (i=0; i<=n; i++){
    printf("i=%2d w=%18.15f\n",i,w[i]);
  }

  for (i=0; i<=n; i++){
    if (i<rankN){
       w2[i] = sqrt(w[i]);
      iw2[i] = 1.0/w2[i];
    }
    else{
       w2[i] = 0.0;
      iw2[i] = 0.0;
    }
  }

  /****************************************************
       calculate beta, gamma, and its inverses
  ****************************************************/

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      B[i][j] = u1[i][j]*w2[j];
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      IB[i][j] = u1[j][i]*iw2[i];
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      C[i][j] = vt1[i][j]*w2[i];
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      IC[i][j] = vt1[j][i]*iw2[j];
    }
  }

  /****************************************************
    freeing of arrays:

    static double a1[YOUSO7*YOUSO7];
    static double vt[YOUSO7*YOUSO7];
    static double u[YOUSO7*YOUSO7];
    static double work[YOUSO7*YOUSO7];
    static double w[YOUSO7+2];
    static double w2[YOUSO7+2];
    static double iw2[YOUSO7+2];
    static double vt1[YOUSO7][YOUSO7];
    static double u1[YOUSO7][YOUSO7];
  ****************************************************/

  /* a1 */
  free(a1);

  /* vt */
  free(vt);

  /* u */
  free(u);

  /* work */
  free(work);

  /* w */
  free(w);

  /* w2 */
  free(w2);

  /* iw2 */
  free(iw2);

  /* vt1 */
  for (i=0; i<(List_YOUSO[7]+2); i++){
    free(vt1[i]);
  }
  free(vt1);

  /* u1 */
  for (i=0; i<(List_YOUSO[7]+2); i++){
    free(u1[i]);
  }
  free(u1);

  /* return */
  return rankN;
}










void T_Conserve_MP(int T_switch)
{

  /****************************************************
    T_switch is a parameter to select a terminator
    in the Green's function.

    T_switch = 1 is no terminator.
    T_switch = 2 is a square root terminator.
  ****************************************************/

  static int firsttime=1;
  static int spin,i,j,k,n,p,po,Gc_AN,rl_loop,doko;
  static int wan,tno,itnum,Num_loop,Mc_AN;
  static double TN,TX,My_TN,My_TX,TZ,TQ,pTQ,Delta,pDelta;
  static double E_Temp0,rc,ac,DChemP,dtnum;
  static double dp,dM,dum,dum1;
  static dcomplex EpC,EpP,CN,CX;
  static dcomplex ***G00;
  static dcomplex *CNd,*CXd;
  static dcomplex cdum,Csum3;
  static int numprocs,myid,ID,tag;

  /* for GBTD00() */
  static dcomplex **cinv;
  static double *be2;
  static double *cterR,*cterI;
  static double **cinvR,**cinvI;
  static double **ctp1R,**ctp1I,**ctp2R,**ctp2I;
  static INTEGER *ipiv;
  static dcomplex *a1d,*work;

  MPI_Status stat;
  MPI_Request request; 

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    allocation of arrays:
  ****************************************************/

  G00 = (dcomplex***)malloc(sizeof(dcomplex**)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    G00[k] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
    for (i=0; i<List_YOUSO[7]; i++){
      G00[k][i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
    }
  }  

  CNd = (dcomplex*)malloc(sizeof(dcomplex)*(Matomnum+1)); 
  CXd = (dcomplex*)malloc(sizeof(dcomplex)*(Matomnum+1)); 

  /* for GBTD00() */

  be2 = (double*)malloc(sizeof(double)*List_YOUSO[7]); 

  cterR = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  cterI = (double*)malloc(sizeof(double)*List_YOUSO[7]); 

  cinv = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinv[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  cinvR = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinvR[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  cinvI = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinvI[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp1R = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp1R[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp1I = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp1I[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp2R = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp2R[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp2I = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp2I[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ipiv = (INTEGER*)malloc(sizeof(INTEGER)*(List_YOUSO[7]+2));
  a1d = (dcomplex*)malloc(sizeof(dcomplex)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2));
  work = (dcomplex*)malloc(sizeof(dcomplex)*(List_YOUSO[7]+2));

  /****************************************************
    start calculation
  ****************************************************/

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ = TZ + Spe_Core_Charge[wan];
  }

  printf("<Recursion>  Ideal num.of.electrons = %15.10f\n",TZ);

  po = 150;
  Num_loop = 0;
  TQ = 10.0;
  ac = 1.0;
  rc = 1.0;
  DChemP = 0.0040*eV2Hartree;
  pTQ = 10000;
  itnum = Av_num;
  dtnum = Av_num;

  while(0<po){

    Num_loop++;

    /****************************************************
                     Total_N and Total_X
    ****************************************************/

    CN = Complex(0.0,0.0);
    CX = Complex(0.0,0.0);

    for (p=0; p<POLES; p++){

      EpC = RCadd(ChemP,Ep[p]);
      EpP.r = EpC.r;
      EpP.i = EpC.i;

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

        Gc_AN = M2G[Mc_AN];
	wan = WhatSpecies[Gc_AN];
	tno = Spe_Total_CNO[wan];

        /****************************************************
                          multiple inverse
        ****************************************************/

        for (spin=0; spin<=SpinP_switch; spin++){
          GBTD00( T_switch, Mc_AN, spin, EpP, G00[spin],
                  cinv, be2, cterR, cterI, cinvR, cinvI,
                  ctp1R, ctp1I, ctp2R, ctp2I,
                  ipiv, a1d, work );
        }

        /****************************************************
                    calculate Mulliken population
        ****************************************************/

	Csum3.r = 0.0;
	Csum3.i = 0.0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (i=0; i<tno; i++){
 	    Csum3.r += G00[spin][i][i].r;
	    Csum3.i += G00[spin][i][i].i;
	  }
        }

	CNd[Mc_AN].r = Csum3.r;
	CNd[Mc_AN].i = Csum3.i;
        
        /****************************************************
                          response functions
        ****************************************************/
        
	Csum3.r = 0.0;
	Csum3.i = 0.0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (i=0; i<tno; i++){
	    Csum3.r += G00[spin][i][i].r*G00[spin][i][i].r - G00[spin][i][i].i*G00[spin][i][i].i;
	    Csum3.i += 2.0*G00[spin][i][i].i*G00[spin][i][i].r;
          }
        }
	CXd[Mc_AN].r = Csum3.r;
	CXd[Mc_AN].i = Csum3.i;

      } /* end of Mc_AN loop */

      /****************************************************
         summation of Mulliken populations and
            responce functions over atoms
      ****************************************************/

      Csum3.r = 0.0;
      Csum3.i = 0.0;

      for (i=1; i<=Matomnum; i++){
	Csum3.r += zp[p].r*CNd[i].r - zp[p].i*CNd[i].i;
	Csum3.i += zp[p].i*CNd[i].r + zp[p].r*CNd[i].i;
      }

      CN.r += Csum3.r;
      CN.i += Csum3.i;
      Csum3.r = 0.0;
      Csum3.i = 0.0;

      for (i=1; i<=Matomnum; i++){
	Csum3.r += zp[p].r*CXd[i].r - zp[p].i*CXd[i].i;
	Csum3.i += zp[p].i*CXd[i].r + zp[p].r*CXd[i].i;
      }

      CX.r += Csum3.r;
      CX.i += Csum3.i;

    } /* end of POLE loop */ 

    if (SpinP_switch==0){
      My_TN = 4.0*CN.r/Beta;
      My_TX = 4.0*CX.r/Beta;
    }
    else if (SpinP_switch==1){
      My_TN = 2.0*CN.r/Beta;
      My_TX = 2.0*CX.r/Beta;
    }

    /****************************************************
      MPI
    ****************************************************/

     MPI_Allreduce(&My_TN, &TN, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
     MPI_Allreduce(&My_TX, &TX, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /****************************************************
               correction of chemical potential
    ****************************************************/

    TQ = TZ - TN - Given_Total_Charge;

    if (fabs(TQ)<(CN_Error*(double)atomnum)) po=0;
    else                                     po--;

    if (0<po){
      if (0.1<fabs(TQ/TZ) || 1.0<fabs(TQ)){
	if (0.0<TQ){
	  ChemP = ChemP + DChemP;
          doko = 1;
	}
	else{
	  ChemP = ChemP - DChemP;
          doko = 2;
	}
	if (fabs(pTQ)<fabs(TQ)) DChemP = 0.5*DChemP;
      }
      else{
	if (fabs(pTQ)<fabs(TQ) && pTQ*TQ<0.0){
	  rc = rc + 3.0;
	  Delta = -0.5*pDelta;
	  pDelta = Delta;
          doko = 3;
	}
	else{
	  Delta = -ac*(TZ-TN-Given_Total_Charge)/TX/rc;
	  pDelta = Delta;
          doko = 4;          

	  if ((0.5<fabs(TQ/pTQ)) && (fabs(TQ/TZ)<0.05)){
	    if (0.0<(TQ/pTQ)){
	      ac = ac + 0.5;
              doko = 5;          
	    }
	    else{
	      ac = ac/3.0;
              doko = 6;
	    }
	  }
	}
	ChemP = ChemP + Delta;
      }
      pTQ = TQ;
    }

    printf("<Recursion> %3d where=%2d  Num.of.eles= %15.10f   ChemP= %15.10f\n",
            Num_loop,doko,TN,eV2Hartree*ChemP);
  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[7]; i++){
      free(G00[k][i]);
    }
    free(G00[k]);
  }  
  free(G00);

  free(CNd);
  free(CXd);

  /* for GBTD00() */

  free(be2);

  free(cterR);
  free(cterI);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinv[i]);
  }
  free(cinv);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinvR[i]);
  }
  free(cinvR);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinvI[i]);
  }
  free(cinvI);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp1R[i]);
  }
  free(ctp1R);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp1I[i]);
  }
  free(ctp1I);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp2R[i]);
  }
  free(ctp2R);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp2I[i]);
  }
  free(ctp2I);

  free(ipiv);
  free(a1d);
  free(work);
}













void Calc_DM_with_S(int T_switch, double *****CDM, double *****EDM, double ****IOLP)
{
  /****************************************************
    T_switch is a parameter to select a terminator
    in the Green's function.
    T_switch = 1 is no terminator.
    T_switch = 2 is a square root terminator.
  ****************************************************/
  static int firsttime=1;
  static int Gc_AN,Mc_AN,spin,rl_loop,Anum,Bnum,Gi,wanA;
  static int wan,tno,itnum,ig,ino1,jg,jno,kg,kj,kno;
  static int i,j,k,p,l,m,n,q,rl,i1,j1,h_AN;
  static double sum,sum2,sumx,sumy,sumz,dtnum;

  static dcomplex EpC,EpP,CE,CRho,CRes,CERho;
  static dcomplex ***G00;
  static dcomplex ****G0;

  static double ****Rho;
  static double ****ERho;
  static double **E00;
  static double *Fx,*Fy,*Fz;
  static double **Uele_temp;
  static double **Dmat,**Dmat2;
  static double dum,dum1,dum2;
  static double my_Uele_IS[3];

  /* for GBTD00() */
  static dcomplex **cinv;
  static double *be2;
  static double *cterR,*cterI;
  static double **cinvR,**cinvI;
  static double **ctp1R,**ctp1I,**ctp2R,**ctp2I;
  static INTEGER *ipiv;
  static dcomplex *a1d,*work;

  /* for RecurG */
  static dcomplex **Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  G00 = (dcomplex***)malloc(sizeof(dcomplex**)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    G00[k] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
    for (i=0; i<List_YOUSO[7]; i++){
      G00[k][i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
    }
  }  

  G0 = (dcomplex****)malloc(sizeof(dcomplex***)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    G0[k] = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[3]);
    for (i=0; i<List_YOUSO[3]; i++){
      G0[k][i] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
        G0[k][i][j] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
      }
    }
  }  

  Rho = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    Rho[k] = (double***)malloc(sizeof(double**)*List_YOUSO[3]);
    for (i=0; i<List_YOUSO[3]; i++){
      Rho[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
        Rho[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }
  }  

  ERho = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    ERho[k] = (double***)malloc(sizeof(double**)*List_YOUSO[3]);
    for (i=0; i<List_YOUSO[3]; i++){
      ERho[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
        ERho[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }
  }  

  E00 = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    E00[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  Fx = (double*)malloc(sizeof(double)*(atomnum+1));
  Fy = (double*)malloc(sizeof(double)*(atomnum+1));
  Fz = (double*)malloc(sizeof(double)*(atomnum+1));

  Uele_temp = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    Uele_temp[k] = (double*)malloc(sizeof(double)*(Matomnum+1)); 
  }

  Dmat = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (j=0; j<List_YOUSO[7]; j++){
    Dmat[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  Dmat2 = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (j=0; j<List_YOUSO[7]; j++){
    Dmat2[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  /* for GBTD00 */

  be2 = (double*)malloc(sizeof(double)*List_YOUSO[7]); 

  cterR = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  cterI = (double*)malloc(sizeof(double)*List_YOUSO[7]); 

  cinv = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinv[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  cinvR = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinvR[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  cinvI = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinvI[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp1R = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp1R[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp1I = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp1I[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp2R = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp2R[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ctp2I = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp2I[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  ipiv = (INTEGER*)malloc(sizeof(INTEGER)*(List_YOUSO[7]+2));
  a1d = (dcomplex*)malloc(sizeof(dcomplex)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2));
  work = (dcomplex*)malloc(sizeof(dcomplex)*(List_YOUSO[7]+2));

  /* for RecurG */
  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  if (firsttime) {
  firsttime=0;
  }

  /****************************************************
                   calculation of DM
  ****************************************************/

  itnum = Av_num;
  dtnum = Av_num;

  for (spin=0; spin<=SpinP_switch; spin++){
    my_Uele_IS[spin] = 0.0;
  }

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan = WhatSpecies[Gc_AN];
    tno = Spe_Total_CNO[wan];

    /* calculation of MP */

    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      MP[i] = Anum;
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum += Spe_Total_CNO[wanA];
    }

    /* initialize E00 */

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<tno; i++) E00[spin][i] = 0.0;
    }

    /* initialize Rho and ERho */

    for (spin=0; spin<=SpinP_switch; spin++){
      for (rl=0; rl<=RNUM[Gc_AN]; rl++){
	for (i=0; i<tno; i++){
	  for (j=0; j<tno; j++){
             Rho[spin][rl][i][j] = 0.0;
            ERho[spin][rl][i][j] = 0.0;
	  }
	}
      }
    }

    /* a loop for p */

    for (p=0; p<POLES; p++){

      EpC.r = ChemP + Ep[p].r;
      EpC.i = Ep[p].i;
      EpP.r = EpC.r;
      EpP.i = EpC.i;

      /****************************************************
                        multiple inverse
      ****************************************************/

      for (spin=0; spin<=SpinP_switch; spin++){
        GBTD00( T_switch, Mc_AN, spin, EpP, G00[spin],
                cinv, be2, cterR, cterI, cinvR, cinvI,
                ctp1R, ctp1I, ctp2R, ctp2I,
                ipiv, a1d, work );
      }

      /****************************************************
                      recurrence relation
      ****************************************************/

      for (spin=0; spin<=SpinP_switch; spin++){
        RecurG5( Mc_AN, spin, EpP, G00[spin], G0[spin], ipiv, a1d, work, Ctp, Ctp1, Ctp2, Ctp3 );
      } 

      for (spin=0; spin<=SpinP_switch; spin++){
        for (rl=0; rl<=RNUM[Gc_AN]; rl++){
	  for (i=0; i<tno; i++){
	    for (j=0; j<tno; j++){

	      CRho.r = zp[p].r*G0[spin][rl][i][j].r
                     - zp[p].i*G0[spin][rl][i][j].i;
	      CRho.i = zp[p].i*G0[spin][rl][i][j].r
                     + zp[p].r*G0[spin][rl][i][j].i;

	      Rho[spin][rl][i][j] += CRho.r;
              CERho.r = CRho.r*EpC.r - CRho.i*EpC.i;
	      ERho[spin][rl][i][j] += CERho.r;

              if (rl==0 && i==j) E00[spin][i] += CERho.r;
            }
	  }
        }
      }

    } /* end of p loop */

    for (spin=0; spin<=SpinP_switch; spin++){
      for (rl=0; rl<=RNUM[Gc_AN]; rl++){
        for (i=0; i<tno; i++){
	  for (j=0; j<tno; j++){
            Rho[spin][rl][i][j]  = 2.0*Rho[spin][rl][i][j]/Beta;
            ERho[spin][rl][i][j] = 2.0*ERho[spin][rl][i][j]/Beta;
	  }
        }
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<tno; i++){
        E00[spin][i] = 2.0*E00[spin][i]/Beta;
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      sum = 0.0;
      for (i=0; i<tno; i++) sum += E00[spin][i];
      Uele_temp[spin][Mc_AN] = sum;
    }

    /****************************************************
                       intersite Uele
    ****************************************************/

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<tno; i++){
        for (k=0; k<tno; k++){
          my_Uele_IS[spin] += Rho[spin][0][i][k]*al[spin][Mc_AN][0][k][i];
        }
      }
      for (i=0; i<tno; i++){
        for (k=0; k<tno; k++){
          my_Uele_IS[spin] += Rho[spin][1][i][k]*ga[spin][Mc_AN][1][k][i];
        }
      }
    }
 



    for (rl=0; rl<=RNUM[Gc_AN]; rl++){

      /*
      printf("Rho  Gc_AN=%2d rl=%2d\n",Gc_AN,rl);
      */

      sum = 0.0;   
      for (i=0; i<tno; i++){
        for (k=0; k<tno; k++){
          sum = sum + fabs(Rho[0][rl][i][k]);
	  /*
          printf("%15.12f ",Rho[0][rl][i][k]);
	  */

	}
	/*
        printf("\n");
	*/
      }

      printf("Gc_AN=%2d rl=%2d  Sum of Rho = %18.15f\n",Gc_AN,rl,sum); 
    }




    /****************************************************
          transformation of DM and EDM from L basis
                  into the original basis
    ****************************************************/

    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<=FNAN[Gc_AN]; j++){
        jg = natn[Gc_AN][j];
        jno = Spe_Total_CNO[WhatSpecies[jg]];
        for (p=0; p<tno; p++){
	  for (q=0; q<jno; q++){                
	    CDM[spin][Mc_AN][j][p][q] = 0.0;
	    EDM[spin][Mc_AN][j][p][q] = 0.0;
	  }
	}
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<=FNAN[Gc_AN]; j++){

        jg = natn[Gc_AN][j];
        jno = Spe_Total_CNO[WhatSpecies[jg]];
        Anum = MP[j]; 

        for (rl=0; rl<=RNUM[Gc_AN]; rl++){
          for (k=0; k<=(FNAN[Gc_AN]+SNAN[Gc_AN]); k++){

            kg = natn[Gc_AN][k];
            kno = Spe_Total_CNO[WhatSpecies[kg]];
            Bnum = MP[k];

            /* Rho * LU^t */

            /* store LU temporary */
	    for (m=0; m<tno; m++){
              for (q=0; q<kno; q++){
                cinvR[q][m] = LU[spin][Mc_AN][rl][m][Bnum+q];
	      }        
	    }    

            for (p=0; p<tno; p++){
              for (q=0; q<kno; q++){

	        sum  = 0.0;
		sum2 = 0.0;
	        for (m=0; m<tno; m++){
                  sum  +=  Rho[spin][rl][p][m]*cinvR[q][m]; 
                  sum2 += ERho[spin][rl][p][m]*cinvR[q][m];
		}
 	        Dmat[p][q]  = sum;
		Dmat2[p][q] = sum2;

	      }
	    } 

            /* Rho * LU^t * S^{-1} */

            /* store S12 into cinvR temporary */
  	    for (m=0; m<kno; m++){
	      for (q=0; q<jno; q++){                
                cinvR[q][m] = S12[Mc_AN][Bnum+m][Anum+q];
	      }        
	    }    

            for (p=0; p<tno; p++){
	      for (q=0; q<jno; q++){                
		sum  = 0.0;
		sum2 = 0.0;
		for (m=0; m<kno; m++){
		  sum  +=  Dmat[p][m]*cinvR[q][m];
		  sum2 += Dmat2[p][m]*cinvR[q][m];
		}
		CDM[spin][Mc_AN][j][p][q] += sum;
		EDM[spin][Mc_AN][j][p][q] += sum2;
	      }
	    }
	  }
	}
      }
    }

  } /* end of Mc_AN loop */

  /*  summation of Uele_temp for Uele_OS  */

  for (spin=0; spin<=SpinP_switch; spin++){
    sum = 0.0;
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      sum += Uele_temp[spin][Mc_AN];
    }

    MPI_Allreduce(&sum, &Uele_OS[spin], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  }

  /*  call MPI to calculate Uele_IS by summing up my_Uele_IS */

  for (spin=0; spin<=SpinP_switch; spin++){
    MPI_Allreduce(&my_Uele_IS[spin], &Uele_IS[spin], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  }

  if (SpinP_switch==0){
    printf("Uele_OS[0]=%15.12f   Uele_OS[1]=%15.12f\n",Uele_OS[0],Uele_OS[0]);
    printf("Uele_IS[0]=%15.12f   Uele_IS[1]=%15.12f\n",Uele_IS[0],Uele_IS[0]);
  }
  else if (SpinP_switch==1){
    printf("Uele_OS[0]=%15.12f   Uele_OS[1]=%15.12f\n",Uele_OS[0],Uele_OS[1]);
    printf("Uele_IS[0]=%15.12f   Uele_IS[1]=%15.12f\n",Uele_IS[0],Uele_IS[1]);
  }


  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[7]; i++){
      free(G00[k][i]);
    }
    free(G00[k]);
  }  
  free(G00);

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[3]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(G0[k][i][j]);
      }
      free(G0[k][i]);
    }
    free(G0[k]);
  }  
  free(G0);

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[3]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(Rho[k][i][j]);
      }
      free(Rho[k][i]);
    }
    free(Rho[k]);
  }  
  free(Rho);

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[3]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(ERho[k][i][j]);
      }
      free(ERho[k][i]);
    }
    free(ERho[k]);
  }  
  free(ERho);

  for (k=0; k<=SpinP_switch; k++){
    free(E00[k]);
  }
  free(E00);

  free(Fx);
  free(Fy);
  free(Fz);

  for (k=0; k<=SpinP_switch; k++){
    free(Uele_temp[k]);
  }
  free(Uele_temp);

  for (j=0; j<List_YOUSO[7]; j++){
    free(Dmat[j]);
  }
  free(Dmat);

  for (j=0; j<List_YOUSO[7]; j++){
    free(Dmat2[j]);
  }
  free(Dmat2);

  /* for GBTD00() */

  free(be2);

  free(cterR);
  free(cterI);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinv[i]);
  }
  free(cinv);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinvR[i]);
  }
  free(cinvR);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinvI[i]);
  }
  free(cinvI);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp1R[i]);
  }
  free(ctp1R);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp1I[i]);
  }
  free(ctp1I);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp2R[i]);
  }
  free(ctp2R);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp2I[i]);
  }
  free(ctp2I);

  free(ipiv);
  free(a1d);
  free(work);

  /* for RecurG */
  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);

}



























void RecurG(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S)
{

  static int i,j,k,rl,wan;
  static int tno,tno1,ct_AN;
  static int rlm,rlm2;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2,**Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3;
  static dcomplex Cdum,Cdum2,Cdum3;

  ct_AN = M2G[Mc_AN];

  /****************************************************
    allocation of arrays:
  ****************************************************/

  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        G0S[rl][i][j].r = 0.0;
        G0S[rl][i][j].i = 0.0;
      }
    }
  }

  /****************************************************
    start recursion loop
  ****************************************************/

  for (rl=1; rl<=RNUM[ct_AN]; rl++){

    rlm = rl - 1;
    rlm2 = rl - 2;

    /****************************************************
                        G0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  Ctp[i][j].r = 0.0;
	  Ctp[i][j].i = 0.0;
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rlm2]; k++){
	    Cdum.r = G0S[rlm2][i][k].r*be[spin][Mc_AN][rlm][k][j];
	    Cdum.i = G0S[rlm2][i][k].i*be[spin][Mc_AN][rlm][k][j];
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
    }

    /****************************************************
                       G0(n-1)*(Z-A(n-1))
    ****************************************************/

    for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
        dum = -al[spin][Mc_AN][rlm][i][j];
        if (i==j){
          Ctp1[i][i].r = dum + EpP.r;
          Ctp1[i][i].i = EpP.i;
        } 
        else{
          Ctp1[i][j].r = dum;
          Ctp1[i][j].i = 0.0;
        } 
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
                 - G0S[rlm][i][k].i*Ctp1[k][j].i;
          Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
                 + G0S[rlm][i][k].r*Ctp1[k][j].i;
          Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
        }
        Ctp2[i][j].r = Csum.r;
        Ctp2[i][j].i = Csum.i;
      }
    }

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  if (i==j){
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r - 1.0;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	  else{
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }
    }

    /****************************************************
                        [ ]*(Cn)^{-1}
    ****************************************************/

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = Ctp3[i][k].r*iga[spin][Mc_AN][rl][k][j];
	  Cdum.i = Ctp3[i][k].i*iga[spin][Mc_AN][rl][k][j];
	  Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
	}
	G0S[rl][i][j].r = Csum.r;
	G0S[rl][i][j].i = Csum.i;
      }
    }
  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);
}











void RecurG2(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S)
{
  /****************************************************
    scaled recurrence relation
  ****************************************************/

  static int i,j,k,rl,wan;
  static int tno,tno1,ct_AN;
  static int rlm,rlm2;
  static double dum,dum1,dum2,scaleF,scaleF2;
  static dcomplex Csum,Csum2,**Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3;
  static dcomplex Cdum,Cdum2,Cdum3;

  scaleF  = 1.0e-7;
  scaleF2 = scaleF*scaleF;

  ct_AN = M2G[Mc_AN];

  /****************************************************
    allocation of arrays:
  ****************************************************/

  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }


  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }
  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    rlm = rl - 1;
    rlm2 = rl - 2;

    /****************************************************
                        g0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  Ctp[i][j].r = 0.0;
	  Ctp[i][j].i = 0.0;
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rlm2]; k++){
	    Cdum.r = G0S[rlm2][i][k].r*be[spin][Mc_AN][rlm][k][j];
	    Cdum.i = G0S[rlm2][i][k].i*be[spin][Mc_AN][rlm][k][j];
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
    }

    /****************************************************
                     scaleF*g0(n-1)*(Z-A(n-1))
    ****************************************************/

    for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
        dum = -al[spin][Mc_AN][rlm][i][j];
        if (i==j){
          Ctp1[i][i].r = dum + EpP.r;
          Ctp1[i][i].i = EpP.i;
        } 
        else{
          Ctp1[i][j].r = dum;
          Ctp1[i][j].i = 0.0;
        } 
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r - G0S[rlm][i][k].i*Ctp1[k][j].i;
          Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r + G0S[rlm][i][k].r*Ctp1[k][j].i;
          Csum.r += Cdum.r;
	  Csum.i += Cdum.i;
        }
        Ctp2[i][j].r = scaleF*Csum.r;
        Ctp2[i][j].i = scaleF*Csum.i;
      }
    }

    /****************************************************
         scaleF*g0(n-1)*(Z-A(n-1)) - g0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  if (i==j){
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r - scaleF;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	  else{
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }
    }

    /****************************************************
                        [ ]*(Cn)^{-1}
    ****************************************************/

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = Ctp3[i][k].r*iga[spin][Mc_AN][rl][k][j];
	  Cdum.i = Ctp3[i][k].i*iga[spin][Mc_AN][rl][k][j];
	  Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
	}
	G0S[rl][i][j].r = Csum.r/scaleF2;
	G0S[rl][i][j].i = Csum.i/scaleF2;
      }
    }
  }

  /****************************************************
    rescaling of G0S
  ****************************************************/

  dum = 1.0;
  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    dum *= scaleF; 
    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	G0S[rl][i][j].r *= dum;
	G0S[rl][i][j].i *= dum; 
      }
    }
  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);
}









void RecurG3(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
             INTEGER *ipiv, dcomplex *a1d, dcomplex *work)
{
  /****************************************************
    decreasing order
  ****************************************************/

  static int i,j,k,rl,wan;
  static int tno,tno1,ct_AN;
  static int rlm,rlm2,match_rl;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2,**Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3,**match_G;
  static dcomplex Cdum,Cdum2,Cdum3;

  ct_AN = M2G[Mc_AN];

  /****************************************************
    allocation of arrays:
  ****************************************************/

  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  match_G = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    match_G[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        G0S[rl][i][j].r = 0.0;
        G0S[rl][i][j].i = 0.0;
      }
    }
  }

  /****************************************************
    start recursion loop
  ****************************************************/

  for (rl=1; rl<=RNUM[ct_AN]; rl++){

    rlm = rl - 1;
    rlm2 = rl - 2;

    /****************************************************
                        G0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  Ctp[i][j].r = 0.0;
	  Ctp[i][j].i = 0.0;
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rlm2]; k++){
	    Cdum.r = G0S[rlm2][i][k].r*be[spin][Mc_AN][rlm][k][j];
	    Cdum.i = G0S[rlm2][i][k].i*be[spin][Mc_AN][rlm][k][j];
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
    }

    /****************************************************
                       G0(n-1)*(Z-A(n-1))
    ****************************************************/

    for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
        dum = -al[spin][Mc_AN][rlm][i][j];
        if (i==j){
          Ctp1[i][i].r = dum + EpP.r;
          Ctp1[i][i].i = EpP.i;
        } 
        else{
          Ctp1[i][j].r = dum;
          Ctp1[i][j].i = 0.0;
        } 
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
                 - G0S[rlm][i][k].i*Ctp1[k][j].i;
          Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
                 + G0S[rlm][i][k].r*Ctp1[k][j].i;
          Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
        }
        Ctp2[i][j].r = Csum.r;
        Ctp2[i][j].i = Csum.i;
      }
    }

    /****************************************************
           G0(n-1)*(Z-A(n-1)) - G0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  if (i==j){
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r - 1.0;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	  else{
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }
    }

    /****************************************************
                        [ ]*(Cn)^{-1}
    ****************************************************/

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = Ctp3[i][k].r*iga[spin][Mc_AN][rl][k][j];
	  Cdum.i = Ctp3[i][k].i*iga[spin][Mc_AN][rl][k][j];
	  Csum.r += Cdum.r;
	  Csum.i += Cdum.i;
	}
	G0S[rl][i][j].r = Csum.r;
	G0S[rl][i][j].i = Csum.i;
      }
    }

  } /* rl */




  /*
  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    printf("\n");
    printf("A G0%2d.r\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].r);
      }
      printf("\n");
    }
    printf("A G0%2d.i\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].i);
      }
      printf("\n");
    }
  }
  */
















  match_rl = 10;

  if (match_rl<(RNUM[ct_AN]-1) ){

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	match_G[i][j].r = G0S[match_rl][i][j].r;
	match_G[i][j].i = G0S[match_rl][i][j].i;
      }
    }


    dum = 0.0;

    rl = RNUM[ct_AN];
    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	if (dum<fabs(G0S[rl][i][j].r)) dum = fabs(G0S[rl][i][j].r);
	if (dum<fabs(G0S[rl][i][j].i)) dum = fabs(G0S[rl][i][j].i);
      }
    }

    for (rl=RNUM[ct_AN]-1; rl<=RNUM[ct_AN]; rl++){
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	  G0S[rl][i][j].r /= dum;
	  G0S[rl][i][j].i /= dum;
	}
      }
    }




    for (rl=(RNUM[ct_AN]-2); match_rl<=rl; rl--){

      /****************************************************
                        G0(n+2)*C(n+2)
      ****************************************************/

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+2]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl+2]; k++){
	    Cdum.r = G0S[rl+2][i][k].r*ga[spin][Mc_AN][rl+2][k][j];
	    Cdum.i = G0S[rl+2][i][k].i*ga[spin][Mc_AN][rl+2][k][j];
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
     
      /****************************************************
                       G0(n+1)*(Z-A(n+1))
      ****************************************************/

      for (i=0; i<RCN_Rank[Mc_AN][rl+1]; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  dum = -al[spin][Mc_AN][rl+1][i][j];
	  if (i==j){
	    Ctp1[i][i].r = dum + EpP.r;
	    Ctp1[i][i].i = EpP.i;
	  } 
	  else{
	    Ctp1[i][j].r = dum;
	    Ctp1[i][j].i = 0.0;
	  } 
	}
      }

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl+1]; k++){
	    Cdum.r = G0S[rl+1][i][k].r*Ctp1[k][j].r
	      - G0S[rl+1][i][k].i*Ctp1[k][j].i;
	    Cdum.i = G0S[rl+1][i][k].i*Ctp1[k][j].r
	      + G0S[rl+1][i][k].r*Ctp1[k][j].i;
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp2[i][j].r = Csum.r;
	  Ctp2[i][j].i = Csum.i;
	}
      }

      /****************************************************
           G0(n+1)*(Z-A(n+1)) - G0(n+2)*C(n+2)
      ****************************************************/

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }

      /****************************************************
                        [ ]*(Bn+1)^{-1}
      ****************************************************/

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl+1]; k++){
	    Cdum.r = Ctp3[i][k].r*ibe[spin][Mc_AN][rl+1][k][j];
	    Cdum.i = Ctp3[i][k].i*ibe[spin][Mc_AN][rl+1][k][j];
	    Csum.r += Cdum.r;
	    Csum.i += Cdum.i;
	  }
	  G0S[rl][i][j].r = Csum.r;
	  G0S[rl][i][j].i = Csum.i;
	}
      }

    } /* rl */












    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	Ctp3[i][j].r = G0S[match_rl][i][j].r;
	Ctp3[i][j].i = G0S[match_rl][i][j].i;
      }
    }


    Inverse_ComplexMat(RCN_Rank[Mc_AN][match_rl]-1, G0S[match_rl], ipiv, a1d, work);

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][match_rl]; k++){
	  Cdum.r = match_G[i][k].r*G0S[match_rl][k][j].r
	    - match_G[i][k].i*G0S[match_rl][k][j].i;
	  Cdum.i = match_G[i][k].i*G0S[match_rl][k][j].r
	    + match_G[i][k].r*G0S[match_rl][k][j].i;
	  Csum.r += Cdum.r;
	  Csum.i += Cdum.i;
	}
	Ctp[i][j].r = Csum.r;
	Ctp[i][j].i = Csum.i;
      }
    }


    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	G0S[match_rl][i][j].r = Ctp3[i][j].r;
	G0S[match_rl][i][j].i = Ctp3[i][j].i;
      }
    }



    for (rl=RNUM[ct_AN]; match_rl<=rl; rl--){

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl]; k++){
	    Cdum.r = Ctp[i][k].r*G0S[rl][k][j].r
	      - Ctp[i][k].i*G0S[rl][k][j].i;
	    Cdum.i = Ctp[i][k].i*G0S[rl][k][j].r
	      + Ctp[i][k].r*G0S[rl][k][j].i;
	    Csum.r += Cdum.r;
	    Csum.i += Cdum.i;
	  }
	  Ctp2[i][j].r = Csum.r;
	  Ctp2[i][j].i = Csum.i;
	}
      }

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	  G0S[rl][i][j].r = Ctp2[i][j].r;
	  G0S[rl][i][j].i = Ctp2[i][j].i;
	}
      }
    }








    /*
      if (match_rl<RNUM[ct_AN]){
      for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
      G0S[match_rl][i][j].r = match_G[i][j].r;
      G0S[match_rl][i][j].i = match_G[i][j].i;
      }
      }
      }
    */






    /*
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
      printf("\n");
      printf("B G0%2d.r\n",rl);
      for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
      printf("%15.12f ",G0S[rl][i][j].r);
      }
      printf("\n");
      }
      printf("B G0%2d.i\n",rl);
      for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
      printf("%15.12f ",G0S[rl][i][j].i);
      }
      printf("\n");
      }
      }

      MPI_Finalize();
      exit(0);
    */


  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);

  for (i=0; i<List_YOUSO[7]; i++){
    free(match_G[i]);
  }
  free(match_G);

}









void RecurG4(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
             INTEGER *ipiv, dcomplex *a1d, dcomplex *work)
{
  /****************************************************
    decreasing order
  ****************************************************/

  static int i,j,k,rl,wan;
  static int tno,tno1,ct_AN;
  static int rlm,rlm2,match_rl;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2,**Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3,**match_G;
  static dcomplex Cdum,Cdum2,Cdum3;

  ct_AN = M2G[Mc_AN];

  /****************************************************
    allocation of arrays:
  ****************************************************/

  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  match_G = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    match_G[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        G0S[rl][i][j].r = 0.0;
        G0S[rl][i][j].i = 0.0;
      }
    }
  }

  /****************************************************
    start recursion loop
  ****************************************************/

  for (rl=1; rl<=RNUM[ct_AN]; rl++){

    rlm = rl - 1;
    rlm2 = rl - 2;

    /****************************************************
                        G0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  Ctp[i][j].r = 0.0;
	  Ctp[i][j].i = 0.0;
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rlm2]; k++){
	    Cdum.r = G0S[rlm2][i][k].r*be[spin][Mc_AN][rlm][k][j];
	    Cdum.i = G0S[rlm2][i][k].i*be[spin][Mc_AN][rlm][k][j];
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
    }

    /****************************************************
                       G0(n-1)*(Z-A(n-1))
    ****************************************************/

    for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
        dum = -al[spin][Mc_AN][rlm][i][j];
        if (i==j){
          Ctp1[i][i].r = dum + EpP.r;
          Ctp1[i][i].i = EpP.i;
        } 
        else{
          Ctp1[i][j].r = dum;
          Ctp1[i][j].i = 0.0;
        } 
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
                 - G0S[rlm][i][k].i*Ctp1[k][j].i;
          Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
                 + G0S[rlm][i][k].r*Ctp1[k][j].i;
          Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
        }
        Ctp2[i][j].r = Csum.r;
        Ctp2[i][j].i = Csum.i;
      }
    }

    /****************************************************
           G0(n-1)*(Z-A(n-1)) - G0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  if (i==j){
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r - 1.0;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	  else{
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }
    }

    /****************************************************
                        [ ]*(Cn)^{-1}
    ****************************************************/

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	  Cdum.r = Ctp3[i][k].r*iga[spin][Mc_AN][rl][k][j];
	  Cdum.i = Ctp3[i][k].i*iga[spin][Mc_AN][rl][k][j];
	  Csum.r += Cdum.r;
	  Csum.i += Cdum.i;
	}
	G0S[rl][i][j].r = Csum.r;
	G0S[rl][i][j].i = Csum.i;
      }
    }

  } /* rl */




  /*
  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    printf("\n");
    printf("A G0%2d.r\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].r);
      }
      printf("\n");
    }
    printf("A G0%2d.i\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].i);
      }
      printf("\n");
    }
  }
  */





  /*
  rl = RNUM[ct_AN] + 1;
  rlm = rl - 1;
  rlm2 = rl - 2;

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<RCN_Rank[Mc_AN][rlm2]; k++){
	Cdum.r = G0S[rlm2][i][k].r*be[spin][Mc_AN][rlm][k][j];
	Cdum.i = G0S[rlm2][i][k].i*be[spin][Mc_AN][rlm][k][j];
	Csum.r = Csum.r + Cdum.r;
	Csum.i = Csum.i + Cdum.i;
      }
      Ctp[i][j].r = Csum.r;
      Ctp[i][j].i = Csum.i;
    }
  }



  for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      dum = -al[spin][Mc_AN][rlm][i][j];
      if (i==j){
	Ctp1[i][i].r = dum + EpP.r;
	Ctp1[i][i].i = EpP.i;
      } 
      else{
	Ctp1[i][j].r = dum;
	Ctp1[i][j].i = 0.0;
      } 
    }
  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
	  - G0S[rlm][i][k].i*Ctp1[k][j].i;
	Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
	  + G0S[rlm][i][k].r*Ctp1[k][j].i;
	Csum.r = Csum.r + Cdum.r;
	Csum.i = Csum.i + Cdum.i;
      }
      Ctp2[i][j].r = Csum.r;
      Ctp2[i][j].i = Csum.i;
    }
  }


  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
      Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
    }
  }


  printf("Real\n");

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      printf("%15.12f ",Ctp3[i][j].r); 
    }
    printf("\n");
  }

  printf("Im\n");

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      printf("%15.12f ",Ctp3[i][j].r); 
    }
    printf("\n");
  }
  */





  
















  match_rl = 4;

  if (match_rl<RNUM[ct_AN] ){

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	match_G[i][j].r = G0S[match_rl][i][j].r;
	match_G[i][j].i = G0S[match_rl][i][j].i;
      }
    }

    rl = RNUM[ct_AN];
    dum = 0.0;

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	if (dum<fabs(G0S[rl][i][j].r)) dum = fabs(G0S[rl][i][j].r);
	if (dum<fabs(G0S[rl][i][j].i)) dum = fabs(G0S[rl][i][j].i);
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	G0S[rl][i][j].r /= dum;
	G0S[rl][i][j].i /= dum;
      }
    }


    for (rl=(RNUM[ct_AN]-1); match_rl<=rl; rl--){

      /****************************************************
                        G0(n+2)*C(n+2)
      ****************************************************/
  
      if (rl==(RNUM[ct_AN]-1)){
        for (i=0; i<=tno1; i++){
	  for (j=0; j<=tno1; j++){
	    Ctp[i][j].r = 0.0;
	    Ctp[i][j].i = 0.0;
	  }
        }
      }

      else{
	for (i=0; i<=tno1; i++){
	  for (j=0; j<RCN_Rank[Mc_AN][rl+2]; j++){
	    Csum.r = 0.0;
	    Csum.i = 0.0;
	    for (k=0; k<RCN_Rank[Mc_AN][rl+2]; k++){
	      Cdum.r = G0S[rl+2][i][k].r*ga[spin][Mc_AN][rl+2][k][j];
	      Cdum.i = G0S[rl+2][i][k].i*ga[spin][Mc_AN][rl+2][k][j];
	      Csum.r = Csum.r + Cdum.r;
	      Csum.i = Csum.i + Cdum.i;
	    }
	    Ctp[i][j].r = Csum.r;
	    Ctp[i][j].i = Csum.i;
	  }
	}
      }
     
      /****************************************************
                       G0(n+1)*(Z-A(n+1))
      ****************************************************/

      for (i=0; i<RCN_Rank[Mc_AN][rl+1]; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  dum = -al[spin][Mc_AN][rl+1][i][j];
	  if (i==j){
	    Ctp1[i][i].r = dum + EpP.r;
	    Ctp1[i][i].i = EpP.i;
	  } 
	  else{
	    Ctp1[i][j].r = dum;
	    Ctp1[i][j].i = 0.0;
	  } 
	}
      }

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl+1]; k++){
	    Cdum.r = G0S[rl+1][i][k].r*Ctp1[k][j].r
	           - G0S[rl+1][i][k].i*Ctp1[k][j].i;
	    Cdum.i = G0S[rl+1][i][k].i*Ctp1[k][j].r
	           + G0S[rl+1][i][k].r*Ctp1[k][j].i;
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp2[i][j].r = Csum.r;
	  Ctp2[i][j].i = Csum.i;
	}
      }

      /****************************************************
           G0(n+1)*(Z-A(n+1)) - G0(n+2)*C(n+2)
      ****************************************************/

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }

      /****************************************************
                        [ ]*(Bn+1)^{-1}
      ****************************************************/

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl+1]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl+1]; k++){
	    Cdum.r = Ctp3[i][k].r*ibe[spin][Mc_AN][rl+1][k][j];
	    Cdum.i = Ctp3[i][k].i*ibe[spin][Mc_AN][rl+1][k][j];
	    Csum.r += Cdum.r;
	    Csum.i += Cdum.i;
	  }
	  G0S[rl][i][j].r = Csum.r;
	  G0S[rl][i][j].i = Csum.i;
	}
      }

    } /* rl */












    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	Ctp3[i][j].r = G0S[match_rl][i][j].r;
	Ctp3[i][j].i = G0S[match_rl][i][j].i;
      }
    }

    Inverse_ComplexMat(RCN_Rank[Mc_AN][match_rl]-1,G0S[match_rl],ipiv, a1d, work);

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][match_rl]; k++){
	  Cdum.r = match_G[i][k].r*G0S[match_rl][k][j].r
	         - match_G[i][k].i*G0S[match_rl][k][j].i;
	  Cdum.i = match_G[i][k].i*G0S[match_rl][k][j].r
	         + match_G[i][k].r*G0S[match_rl][k][j].i;
	  Csum.r += Cdum.r;
	  Csum.i += Cdum.i;
	}
	Ctp[i][j].r = Csum.r;
	Ctp[i][j].i = Csum.i;
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
	G0S[match_rl][i][j].r = Ctp3[i][j].r;
	G0S[match_rl][i][j].i = Ctp3[i][j].i;
      }
    }



    for (rl=RNUM[ct_AN]; match_rl<=rl; rl--){

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[Mc_AN][rl]; k++){
	    Cdum.r = Ctp[i][k].r*G0S[rl][k][j].r
	           - Ctp[i][k].i*G0S[rl][k][j].i;
	    Cdum.i = Ctp[i][k].i*G0S[rl][k][j].r
	           + Ctp[i][k].r*G0S[rl][k][j].i;
	    Csum.r += Cdum.r;
	    Csum.i += Cdum.i;
	  }
	  Ctp2[i][j].r = Csum.r;
	  Ctp2[i][j].i = Csum.i;
	}
      }

      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	  G0S[rl][i][j].r = Ctp2[i][j].r;
	  G0S[rl][i][j].i = Ctp2[i][j].i;
	}
      }
    }








    /*
      if (match_rl<RNUM[ct_AN]){
      for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][match_rl]; j++){
      G0S[match_rl][i][j].r = match_G[i][j].r;
      G0S[match_rl][i][j].i = match_G[i][j].i;
      }
      }
      }
    */





    /*
    for (rl=0; rl<=RNUM[ct_AN]; rl++){
      printf("\n");
      printf("B G0%2d.r\n",rl);
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  printf("%15.12f ",G0S[rl][i][j].r);
	}
	printf("\n");
      }
      printf("B G0%2d.i\n",rl);
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  printf("%15.12f ",G0S[rl][i][j].i);
	}
	printf("\n");
      }
    }

    MPI_Finalize();
    exit(0);
    */


  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);

  for (i=0; i<List_YOUSO[7]; i++){
    free(match_G[i]);
  }
  free(match_G);

}













void RecurG5(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
             INTEGER *ipiv, dcomplex *a1d, dcomplex *work,
             dcomplex **Ctp, dcomplex **Ctp1, dcomplex **Ctp2, dcomplex **Ctp3 )
{
  /****************************************************
    decreasing order
  ****************************************************/

  static int i,j,k,rl,wan;
  static int tno,tno1,ct_AN;
  static int rlm,rlm2;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2;
  static dcomplex Cdum,Cdum2,Cdum3;

  ct_AN = M2G[Mc_AN];
  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        G0S[rl][i][j].r = 0.0;
        G0S[rl][i][j].i = 0.0;
      }
    }
  }

  /****************************************************
    calculate G01
  ****************************************************/

  rl = 1;
  rlm = rl - 1;

  /****************************************************
                    G00*(Z-A0) - I
  ****************************************************/

  for (i=0; i<tno; i++){
    for (j=0; j<tno; j++){
      /* transpose fot the later calculation */ 
      Ctp1[j][i].r = -al[spin][Mc_AN][rlm][i][j];
      Ctp1[j][i].i = 0.0;
    }
    Ctp1[i][i].r += EpP.r;
    Ctp1[i][i].i += EpP.i;
  }

  for (i=0; i<tno; i++){
    for (j=0; j<tno; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<tno; k++){
        /* due to the previous transposition of Ctp1 */
	Cdum.r = G0S[rlm][i][k].r*Ctp1[j][k].r - G0S[rlm][i][k].i*Ctp1[j][k].i;
	Cdum.i = G0S[rlm][i][k].i*Ctp1[j][k].r + G0S[rlm][i][k].r*Ctp1[j][k].i;
	Csum.r += Cdum.r;
	Csum.i += Cdum.i;
      }
      Ctp3[i][j].r = Csum.r;
      Ctp3[i][j].i = Csum.i;
    }
    Ctp3[i][i].r -= 1.0;
  }

  /****************************************************
              G01 = [G00*(Z-A0) - I]*(C1)^{-1}
  ****************************************************/

  for (i=0; i<tno; i++){
    for (j=0; j<tno; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<tno; k++){
	Cdum.r = Ctp3[i][k].r*iga[spin][Mc_AN][rl][k][j];
	Cdum.i = Ctp3[i][k].i*iga[spin][Mc_AN][rl][k][j];
	Csum.r += Cdum.r;
	Csum.i += Cdum.i;
      }
      G0S[rl][i][j].r = Csum.r;
      G0S[rl][i][j].i = Csum.i;
    }
  }

  /****************************************************
        calculate a series of ratio functions 
        by the following recurrence relation

        R_n = Cn*[ Z - An - R_{n+1}*B_{n+1} ]^{-1}
  ****************************************************/

  for (rl=RNUM[ct_AN]; 2<=rl; rl--){

    /* use Ctp3 to store R_{n+1} */


    /* set R_{RNUM[ct_AN]+1} */
    /* This part can be modified by SRT. */

    if (rl==RNUM[ct_AN]){
      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){

          Ctp3[i][j].r = 0.0;
          Ctp3[i][j].i = 0.0;

          /* store the transpositin */
          G0S[rl+1][j][i].r = 0.0;
          G0S[rl+1][j][i].i = 0.0;
	}
      }     
    }

    /****************************************************
          R_n = Cn*[ Z - An - R_{n+1}*B_{n+1} ]^{-1}
    ****************************************************/

    /* R_{n+1}*B_{n+1} */    
    /* This part can be modified by SRT. */

    if (rl==RNUM[ct_AN]){
      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){
	  Ctp[i][j].r = 0.0;
	  Ctp[i][j].i = 0.0;
	}
      }
    }
    else{
      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<tno; k++){
	    Csum.r += Ctp3[i][k].r*tga[spin][Mc_AN][rl+1][j][k];
	    Csum.i += Ctp3[i][k].i*tga[spin][Mc_AN][rl+1][j][k];
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
    }

    /* Z - An - R_{n+1}*B_{n+1} */    

    for (i=0; i<tno; i++){
      for (j=0; j<tno; j++){
	Ctp1[i][j].r =  -Ctp[i][j].r - al[spin][Mc_AN][rl][i][j];
	Ctp1[i][j].i =  -Ctp[i][j].i;
      }
      Ctp1[i][i].r += EpP.r;
      Ctp1[i][i].i += EpP.i;
    }

    /* [Z - An - R_{n+1}*B_{n+1}]^{-1} */    
    
    Inverse_ComplexMat( tno-1, Ctp1, ipiv, a1d, work );

    /* transposition of Ctp1 */
    for (i=0; i<=tno1; i++){
      for (j=(i+1); j<=tno1; j++){
        Cdum  = Ctp1[i][j];
        Cdum2 = Ctp1[j][i]; 
        Ctp1[i][j] = Cdum2;
        Ctp1[j][i] = Cdum;
      }
    }

    /* Cn*[Z - An - R_{n+1}*B_{n+1}]^{-1} */    

    for (i=0; i<tno; i++){
      for (j=0; j<tno; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<tno; k++){
	  Csum.r += be[spin][Mc_AN][rl][i][k]*Ctp1[j][k].r;
	  Csum.i += be[spin][Mc_AN][rl][i][k]*Ctp1[j][k].i;
	}

	/* store R_n into G0S as the transposition */

	G0S[rl][j][i].r = Csum.r;
	G0S[rl][j][i].i = Csum.i;

	/* store R_n into Ctp3 */

	Ctp3[i][j].r = Csum.r;
	Ctp3[i][j].i = Csum.i;
      }
    }

  } /* rl */

  /****************************************************
        calculate off-diagonal Green's functions
            by G_{0,n} = G_{0,n-1} * R_{n}
  ****************************************************/

  for (rl=2; rl<=RNUM[ct_AN]; rl++){

    for (i=0; i<tno; i++){
      for (j=0; j<tno; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<tno; k++){
          /* G0S_{rl} stores the transposition of R_{rl} */       
	  Csum.r += G0S[rl-1][i][k].r*G0S[rl][j][k].r - G0S[rl-1][i][k].i*G0S[rl][j][k].i;
	  Csum.i += G0S[rl-1][i][k].r*G0S[rl][j][k].i + G0S[rl-1][i][k].i*G0S[rl][j][k].r;
	}
	Ctp[i][j].r = Csum.r;
	Ctp[i][j].i = Csum.i;
      }
    }

    for (i=0; i<tno; i++){
      for (j=0; j<tno; j++){
	G0S[rl][i][j].r = Ctp[i][j].r;
	G0S[rl][i][j].i = Ctp[i][j].i;
      }
    }
  }




  /*
  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    printf("\n");
    printf("A G0%2d.r\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].r);
      }
      printf("\n");
    }

    printf("A G0%2d.i\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].i);
      }
      printf("\n");
    }
  }

  MPI_Finalize();
  exit(0);
  */

}






void RecurG6(int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S,
             INTEGER *ipiv, dcomplex *a1d, dcomplex *work)
{
  /****************************************************
    increasing order
  ****************************************************/

  static int i,j,k,rl,wan;
  static int tno,tno1,ct_AN;
  static int rlm,rlm2;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2,**Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3;
  static dcomplex Cdum,Cdum2,Cdum3;

  ct_AN = M2G[Mc_AN];

  /****************************************************
    allocation of arrays:
  ****************************************************/

  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        G0S[rl][i][j].r = 0.0;
        G0S[rl][i][j].i = 0.0;
      }
    }
  }

  /****************************************************
    calculate G01
  ****************************************************/

  rl = 1;
  rlm = rl - 1;

  /****************************************************
                    G00*(Z-A0) - I
  ****************************************************/

  for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      dum = -al[spin][Mc_AN][rlm][i][j];
      if (i==j){
	Ctp1[i][i].r = EpP.r + dum;
	Ctp1[i][i].i = EpP.i;
      } 
      else{
	Ctp1[i][j].r = dum;
	Ctp1[i][j].i = 0.0;
      } 
    }
  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
	       - G0S[rlm][i][k].i*Ctp1[k][j].i;
	Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
	       + G0S[rlm][i][k].r*Ctp1[k][j].i;
	Csum.r = Csum.r + Cdum.r;
	Csum.i = Csum.i + Cdum.i;
      }

      if (i==j){
	Ctp3[i][j].r = Csum.r - 1.0;
	Ctp3[i][j].i = Csum.i;
      }
      else{
	Ctp3[i][j].r = Csum.r;
	Ctp3[i][j].i = Csum.i;
      }
    }
  }

  /****************************************************
              G01 = [G00*(Z-A0) - I]*(C1)^{-1}
  ****************************************************/

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	Cdum.r = Ctp3[i][k].r*iga[spin][Mc_AN][rl][k][j];
	Cdum.i = Ctp3[i][k].i*iga[spin][Mc_AN][rl][k][j];
	Csum.r += Cdum.r;
	Csum.i += Cdum.i;
      }
      G0S[rl][i][j].r = Csum.r;
      G0S[rl][i][j].i = Csum.i;
    }
  }

  /****************************************************
    G00^{-1}
  ****************************************************/

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      Ctp[i][j].r = G0S[0][i][j].r;
      Ctp[i][j].i = G0S[0][i][j].i;
    }
  }

  Inverse_ComplexMat(tno1,Ctp, ipiv, a1d, work);

  /****************************************************
    R1 = G00^{-1} * G01
  ****************************************************/

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][1]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<=tno1; k++){
	Csum.r += Ctp[i][k].r*G0S[1][k][j].r - Ctp[i][k].i*G0S[1][k][j].i;
	Csum.i += Ctp[i][k].r*G0S[1][k][j].i + Ctp[i][k].i*G0S[1][k][j].r;
      }
      Ctp1[i][j].r = Csum.r;
      Ctp1[i][j].i = Csum.i;
    }
  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][1]; j++){
      Ctp3[i][j].r = Ctp1[i][j].r;
      Ctp3[i][j].i = Ctp1[i][j].i;
    }
  }

  /****************************************************
        calculate a series of ratio functions 
        by the following recurrence relation

  Rn = [Z - A_{n-1} - R_{n-1}^{-1} * C_{n-1}] * Bn^{-1}
  ****************************************************/

  for (rl=2; rl<=RNUM[ct_AN]; rl++){

    /* R_{n-1}^{-1} */    

    if (rl==2){
      for (i=0; i<RCN_Rank[Mc_AN][rl-1]; i++){
        for (j=0; j<RCN_Rank[Mc_AN][rl-1]; j++){
          Ctp[i][j].r = Ctp3[i][j].r;
          Ctp[i][j].i = Ctp3[i][j].i;
        }
      }
    }
    else{
      for (i=0; i<RCN_Rank[Mc_AN][rl-1]; i++){
        for (j=0; j<RCN_Rank[Mc_AN][rl-1]; j++){
          Ctp[i][j].r = G0S[rl-1][i][j].r;
          Ctp[i][j].i = G0S[rl-1][i][j].i;
        }
      }
    }

    Inverse_ComplexMat(RCN_Rank[Mc_AN][rl-1]-1,Ctp, ipiv, a1d, work);

    /* R_{n-1}^{-1} * C_{n-1} */

    for (i=0; i<RCN_Rank[Mc_AN][rl]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl-1]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rl-1]; k++){
	  Csum.r += Ctp[i][k].r*be[spin][Mc_AN][rl-1][k][j];
	  Csum.i += Ctp[i][k].i*be[spin][Mc_AN][rl-1][k][j];
	}
	Ctp1[i][j].r = Csum.r;
	Ctp1[i][j].i = Csum.i;
      }
    }

    /* Z - A_{n-1} - R_{n-1}^{-1} * C_{n-1} */

    for (i=0; i<RCN_Rank[Mc_AN][rl]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
        if (i==j){
	  Ctp2[i][j].r = EpP.r - al[spin][Mc_AN][rl-1][i][j] - Ctp1[i][j].r;
	  Ctp2[i][j].i = EpP.i                               - Ctp1[i][j].i;
	}
        else{
	  Ctp2[i][j].r =        -al[spin][Mc_AN][rl-1][i][j] - Ctp1[i][j].r;
	  Ctp2[i][j].i =                                     - Ctp1[i][j].i;
	}
      }
    }

    /* Rn = [Z - A_{n-1} - R_{n-1}^{-1} * C_{n-1}] * Bn^{-1} */

    for (i=0; i<RCN_Rank[Mc_AN][rl]; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rl]; k++){
	  Csum.r += Ctp2[i][k].r*iga[spin][Mc_AN][rl][k][j];
	  Csum.i += Ctp2[i][k].i*iga[spin][Mc_AN][rl][k][j];
	}
	G0S[rl][i][j].r = Csum.r;
	G0S[rl][i][j].i = Csum.i;
      }
    }
  }

  /****************************************************
        calculate off-diagonal Green's functions
            by G_{0,n} = G_{0,n-1} * R_{n}
  ****************************************************/

  for (rl=2; rl<=RNUM[ct_AN]; rl++){

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[Mc_AN][rl]; k++){
	  Csum.r += G0S[rl-1][i][k].r*G0S[rl][k][j].r - G0S[rl-1][i][k].i*G0S[rl][k][j].i;
	  Csum.i += G0S[rl-1][i][k].r*G0S[rl][k][j].i + G0S[rl-1][i][k].i*G0S[rl][k][j].r;
	}
	Ctp[i][j].r = Csum.r;
	Ctp[i][j].i = Csum.i;
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[Mc_AN][rl]; j++){
	G0S[rl][i][j].r = Ctp[i][j].r;
	G0S[rl][i][j].i = Ctp[i][j].i;
      }
    }
  }




  /*
  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    printf("\n");
    printf("A G0%2d.r\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].r);
      }
      printf("\n");
    }

    printf("A G0%2d.i\n",rl);
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",G0S[rl][i][j].i);
      }
      printf("\n");
    }
  }

  MPI_Finalize();
  exit(0);
  */







  /*
  rl = RNUM[ct_AN] + 1;
  rlm = rl - 1;
  rlm2 = rl - 2;

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<RCN_Rank[Mc_AN][rlm2]; k++){
	Cdum.r = G0S[rlm2][i][k].r*be[spin][Mc_AN][rlm][k][j];
	Cdum.i = G0S[rlm2][i][k].i*be[spin][Mc_AN][rlm][k][j];
	Csum.r = Csum.r + Cdum.r;
	Csum.i = Csum.i + Cdum.i;
      }
      Ctp[i][j].r = Csum.r;
      Ctp[i][j].i = Csum.i;
    }
  }



  for (i=0; i<RCN_Rank[Mc_AN][rlm]; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      dum = -al[spin][Mc_AN][rlm][i][j];
      if (i==j){
	Ctp1[i][i].r = dum + EpP.r;
	Ctp1[i][i].i = EpP.i;
      } 
      else{
	Ctp1[i][j].r = dum;
	Ctp1[i][j].i = 0.0;
      } 
    }
  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Csum.r = 0.0;
      Csum.i = 0.0;
      for (k=0; k<RCN_Rank[Mc_AN][rlm]; k++){
	Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
	  - G0S[rlm][i][k].i*Ctp1[k][j].i;
	Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
	  + G0S[rlm][i][k].r*Ctp1[k][j].i;
	Csum.r = Csum.r + Cdum.r;
	Csum.i = Csum.i + Cdum.i;
      }
      Ctp2[i][j].r = Csum.r;
      Ctp2[i][j].i = Csum.i;
    }
  }


  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
      Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
    }
  }


  printf("Real\n");

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      printf("%15.12f ",Ctp3[i][j].r); 
    }
    printf("\n");
  }

  printf("Im\n");

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[Mc_AN][rlm]; j++){
      printf("%15.12f ",Ctp3[i][j].r); 
    }
    printf("\n");
  }
  */




  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);
}






void GBTD00(int T_switch, int Mc_AN, int spin, dcomplex EpP, dcomplex **G00S,
            dcomplex **cinv, double *be2, double *cterR, double *cterI,
            double **cinvR, double **cinvI, double **ctp1R, double **ctp1I,
            double **ctp2R, double **ctp2I,
            INTEGER *ipiv, dcomplex *a1d, dcomplex *work)
{
  static int i,j,k,l,m,n,rl,wan,Avnum,ct_AN;
  static int tno,itnum,rl_loop;
  static double dtnum,Av_al,d1,d2;
  static double dum,sum,xd,yd;
  static double c1,s1,c2,s2;
  static double tr,ti,ai,bi;
  static double Csum3R,Csum3I,cdumR,cdumI;

  ct_AN = M2G[Mc_AN];
  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  itnum = Av_num;
  dtnum = Av_num;

  /****************************************************
                       Non-Terminator 
  ****************************************************/

  if (T_switch==1){
    for (i=0; i<tno; i++){
      cterR[i] = 0.0;
      cterI[i] = 0.0;
    }
  }

  /****************************************************
                  square root terminator
  ****************************************************/

  else{

    /****************************************************
                        sum of Bn*Cn
    ****************************************************/

    Avnum = 0;    
    for (i=0; i<tno; i++) be2[i] = 0.0;

    for (rl=RNUM[ct_AN]; (RNUM[ct_AN]-itnum+1)<=rl; rl--){
      Avnum += RCN_Rank[Mc_AN][rl-1];
      for (i=0; i<tno; i++){
	sum = 0.0;
	for (k=0; k<tno; k++){
	  sum += be[spin][Mc_AN][rl][i][k]*tga[spin][Mc_AN][rl][i][k];
	}
	be2[i] += 2.0*sum;
      }
    }

    /****************************************************
                         averaging
    ****************************************************/

    sum = 0.0;
    for (i=0; i<tno; i++) sum += be2[i];

    sum = sum/(double)Avnum;
    for (i=0; i<tno; i++){
      be2[i] = sum;
    }

    Avnum = 0;    
    sum = 0.0;
    for (rl=RNUM[ct_AN]; (RNUM[ct_AN]-itnum+1)<=rl; rl--){
      Avnum += RCN_Rank[Mc_AN][rl];
      for (i=0; i<tno; i++){  
	sum += al[spin][Mc_AN][rl][i][i];
      }
    }
    Av_al = sum/(double)Avnum;

    /****************************************************
         calculation of the square root terminator
    ****************************************************/

    for (i=0; i<tno; i++){

      d1 = EpP.r - Av_al;
      d2 = EpP.i;
      xd = d1*d1 - d2*d2 - 2.0*be2[i];
      yd = 2.0*d1*d2;
      dum = sqrt(xd*xd+yd*yd);
      c1 = xd/dum;
      s1 = yd/dum;

      if (yd<0.0){
	c2 = -sqrt((1.0+c1)*0.5);
	s2 = sqrt(1.0-c2*c2);
      }
      else{
	c2 = sqrt((1.0+c1)*0.5);
	s2 = sqrt(1.0-c2*c2);
      }
          
      dum = sqrt(dum);
      d1 = d1 - dum*c2;
      d2 = d2 - dum*s2;
      cterR[i] = d1/be2[i];
      cterI[i] = d2/be2[i];
    }          
  } /* else */

  /****************************************************
            calculation of multiple inverse
  ****************************************************/

  for (i=0; i<tno; i++){
    for (j=0; j<tno; j++){
      cinvR[i][j] = 0.0;
      cinvI[i][j] = 0.0;
    }
    cinvR[i][i] = cterR[i];
    cinvI[i][i] = cterI[i];
  }

  if (T_switch==1)   rl_loop = RNUM[ct_AN];
  else               rl_loop = RNUM[ct_AN] - 1;

  /* start rl loop */

  for (rl=rl_loop; 0<=rl; rl--){

    if (rl==RNUM[ct_AN]){

      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){
	  /* transpose for calculation of Bn+1*[]^-1 */
	  cinvR[j][i] = -al[spin][Mc_AN][rl][i][j];
	  cinvI[j][i] = 0.0;
	}
	cinvR[i][i] += EpP.r;
	cinvI[i][i] += EpP.i;
      }
    }
    else{

      /****************************************************
                           Bn+1*[]^-1
      ****************************************************/

      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){
	  Csum3R = 0.0;
	  Csum3I = 0.0;
	  for (k=0; k<tno; k++){
	    Csum3R += be[spin][Mc_AN][rl+1][i][k]*cinvR[j][k];
	    Csum3I += be[spin][Mc_AN][rl+1][i][k]*cinvI[j][k];
	  } 
	  ctp1R[i][j] = Csum3R;
	  ctp1I[i][j] = Csum3I;
	}
      }

      /****************************************************
                        Bn+1 * []^-1 * Cn+1
      ****************************************************/

      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){
	  Csum3R = 0.0;
	  Csum3I = 0.0;
	  for (k=0; k<tno; k++){
            Csum3R += ctp1R[i][k]*tga[spin][Mc_AN][rl+1][j][k];
	    Csum3I += ctp1I[i][k]*tga[spin][Mc_AN][rl+1][j][k];
          }
	  ctp2R[i][j] = Csum3R;
	  ctp2I[i][j] = Csum3I;
	}
      }

      /****************************************************
                    Z - Alpha - Bn+1*[]^-1*Cn+1
      ****************************************************/

      for (i=0; i<tno; i++){
	for (j=0; j<tno; j++){
	  cinvR[i][j] = -ctp2R[i][j] - al[spin][Mc_AN][rl][i][j];
	  cinvI[i][j] = -ctp2I[i][j];
	}
	cinvR[i][i] += EpP.r;
	cinvI[i][i] += EpP.i;
      }
    }

    /****************************************************
                 [Z-Alpha-Beta*[]^-1*Beta]^-1
    ****************************************************/

    for (i=0; i<tno; i++){
      for (j=0; j<tno; j++){
        cinv[i][j].r = cinvR[i][j];
        cinv[i][j].i = cinvI[i][j];
      }
    }

    Inverse_ComplexMat( tno-1, cinv, ipiv, a1d, work );

    for (i=0; i<tno; i++){
      for (j=0; j<tno; j++){
        /* transpose for calculation of Bn+1*[]^-1 */
        cinvR[j][i] = cinv[i][j].r;
        cinvI[j][i] = cinv[i][j].i;
      }
    }

  } /* rl loop */

  for (i=0; i<tno; i++){
    for (j=0; j<tno; j++){
      G00S[i][j] = cinv[i][j];
    }
  }
}












void Save_Recursion()
{
  static char fileRN[YOUSO10] = ".rcn";
  static FILE *fp_Rcn;

  strcpy(fileRN,".rcn");
  fnjoint(filepath,filename,fileRN);
  if ((fp_Rcn = fopen(fileRN,"w")) != NULL){
    Output_RcnCof(fp_Rcn);
    Output_LU(fp_Rcn);
    fclose(fp_Rcn);
  }
  else
    printf("Failure of saving the recusion logfile.\n");
}








void Output_LU(FILE *fp_Rcn)
{
  static int ct_AN,rl,i,j,wan1,TNO1,wan2,TNO2;
  static int spin,Gh_AN,h_AN;

  for (spin=0; spin<=SpinP_switch; spin++){

    if (spin==0){
      fprintf(fp_Rcn,"******************************************\n");
      fprintf(fp_Rcn,"   Up spin  - Left side Lanczos vectors   \n");
      fprintf(fp_Rcn,"******************************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"*****************************************\n");
      fprintf(fp_Rcn,"  Down spin - Left side Lanczos vectors  \n");
      fprintf(fp_Rcn,"*****************************************\n");
    }

    /*
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        fprintf(fp_Rcn,"%i %i\n",ct_AN,rl);
        for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){
          Gh_AN = natn[ct_AN][h_AN];
          wan2 = WhatSpecies[Gh_AN];
          TNO2 = Spe_Total_CNO[wan2];

          for (j=0; j<TNO2; j++){
            for (i=0; i<TNO1; i++){
              fprintf(fp_Rcn,"%19.14f ",LU[spin][ct_AN][rl][h_AN][i][j]);
            }
            fprintf(fp_Rcn,"\n");
          }
        }
      }   
    }
    */

  }

}









void Output_RcnCof(FILE *fp_Rcn)
{
  static int spin,ct_AN,rl,i,j,wan1,TNO1;

  for (spin=0; spin<=SpinP_switch; spin++){

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************ Alpha UP ************\n");
      fprintf(fp_Rcn,"**********************************\n");
    } 
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"*********** Alpha DOWN ***********\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",al[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************ Beta UP *************\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"*********** Beta DOWN ************\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",be[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"******* Inverse of Beta UP *******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"****** Inverse of Beta DOWN ******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",ibe[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************* Gamma UP ***********\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************ Gamma DOWN **********\n");
      fprintf(fp_Rcn,"**********************************\n");
    }


    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",ga[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"****** Inverse of Gamma UP *******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"***** Inverse of Gamma DOWN ******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",iga[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

  }

}








void Inverse_S_by_LU(int Mc_AN, double ****OLP0, double ***S12)
{
  static INTEGER N,lda,info,lwork;
  static int Gc_AN,i,Anum,Gi,wanA;
  static int ig,ian,ih,j,kl,jg,jan,Bnum,m,n;

  Gc_AN = M2G[Mc_AN];

  /* making of MP */

  Anum = 0;
  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
    MP[i] = Anum;
    Gi = natn[Gc_AN][i];
    wanA = WhatSpecies[Gi];
    Anum += Spe_Total_CNO[wanA];
  }

  N = Anum;

  /* OLP0 to S12 */

  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){

    ig = natn[Gc_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    Anum = MP[i];
    ih = S_G2M[ig];

    for (j=0; j<=(FNAN[Gc_AN]+SNAN[Gc_AN]); j++){

      kl = RMI1[Mc_AN][i][j];
      jg = natn[Gc_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];
      Bnum = MP[j];

      if (0<=kl){
	for (m=0; m<ian; m++){
	  for (n=0; n<jan; n++){
	    S12[Mc_AN][Anum+m][Bnum+n] = OLP0[ih][kl][m][n];
	  }
	}
      }

      else{
	for (m=0; m<ian; m++){
	  for (n=0; n<jan; n++){
	    S12[Mc_AN][Anum+m][Bnum+n] = 0.0;
	  }
	}
      }
    }
  }

  /* S12 -> LoS */

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
       LoS[i*N+j]= S12[Mc_AN][i][j];
    }
  }

  /* call dgetrf_() in clapack */

  lda = N;
  F77_NAME(dgetrf,DGETRF)(&N, &N, LoS, &lda, piv, &info);

  if (info!=0){
    printf("Error in dgetrf_() which is called from RecursionS  info=%2d\n",info);
  }

  /* call dgetri_() in clapack */

  lwork = N;
  F77_NAME(dgetri,DGETRI)(&N, LoS, &lda, piv, work, &lwork, &info);
  if (info!=0){
    printf("Error in dgetri_() which is called from RecursionS  info=%2d\n",info);
  }

  /* LoS -> S12 */

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      S12[Mc_AN][i+1][j+1] = LoS[i*N+j];
    }
  }
}










void Inverse_ComplexMat(int n, dcomplex **a, INTEGER *ipiv,
                        dcomplex *a1d, dcomplex *work)
{
  static int i,j,N;
  static INTEGER lda,info,lwork;

  N = n + 1;

  lda = N;
  lwork = N;

  /*
  ipiv = (INTEGER*)malloc(sizeof(INTEGER)*N);
  a1d = (dcomplex*)malloc(sizeof(dcomplex)*N*N);
  work = (dcomplex*)malloc(sizeof(dcomplex)*N);
  */

  /****************************************************
      a -> a1d
  ****************************************************/

  for (i=0;i<=n;i++) {
    for (j=0;j<=n;j++) {
       a1d[i*N+j]= a[i][j];
    }
  }

  /****************************************************
                Call zgetrf_() in clapack
  ****************************************************/

  F77_NAME(zgetrf,ZGETRF)(&N, &N, a1d, &lda, ipiv, &info);

  if (info!=0){
    printf("Error in zgetrf_() which is called from RecursionS info=%2d\n",info);
  }

  /****************************************************
                Call zgetri_() in clapack
  ****************************************************/

  F77_NAME(zgetri,ZGETRI)(&N, a1d, &lda, ipiv, work, &lwork, &info);
  if (info!=0){
    printf("Error in zgetri_() which is called from RecursionS info=%2d\n",info);
  }

  /****************************************************
               a1d -> a
  ****************************************************/

  for (i=0; i<=n; i++) {
    for (j=0; j<=n; j++) {
      a[i][j] = a1d[i*N+j];
    }
  }
}
