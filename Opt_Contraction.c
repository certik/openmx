/**********************************************************************
  Opt_Contraction.c:

    Opt_Contraction.c is a subroutine to update the contraction
    coefficients using the gradient of the total energy with respect
    to the contraction coefficients for the orbital optimization method.

  Log of Opt_Contraction.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

static void D_UCntCoes(double ***CntCoes0, double ***D_CntCoes,
                       double *****H,      double *****OLP,
                       double *****CDM,    double *****EDM);
static void D_RCntCoes(double ***CntCoes0, double ***D_CntCoes,
                       double *****H,      double *****OLP,
                       double *****CDM,    double *****EDM);
static void D_ACntCoes(double ***CntCoes0,
                       double ***D_CntCoes,
                       double ***DCntCoes_Spe,
                       double *****H, double *****OLP,
                       double *****CDM,
                       double *****EDM);
static void Change_CntCoes(double alpha,
                           double *****OLP,
                           double ***CntCoes,
                           double ***D_CntCoes,
                           double ***CntCoes0);
static double NormD_UCnt(double ***D_CntCoes);
static double NormD_RCnt(double ***D_CntCoes);
static double NormD_ACnt(double ***DCntCoes_Spe);

double Opt_Contraction(
         double *****H, double *****OLP,
         double *****CDM,
         double *****EDM)
{
  static int firsttime=1;
  int i,j,Mc_AN,Gc_AN,Cwan,be,p;
  int n,num,wan,size1,size2;
  double nd[4],al[4],x01,x12,x20,y12,y01;
  double a,b,c,alpha,norm_by_alpha;
  double norm_deri;
  double ***CntCoes0;
  double ***D_CntCoes;
  double ***D_CntCoes0;
  double ***DCntCoes_Spe;
  double ***DCntCoes_Spe0;
  double *tmp_array;
  double *tmp_array2;
  int *Snd_CntCoes_Size;
  int *Rcv_CntCoes_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /********************************************************
    allocation of arrays:

    int Snd_CntCoes_Size[numprocs];
    int Rcv_CntCoes_Size[numprocs];

    double CntCoes0[Matomnum+MatonumF+1]
                          [List_YOUSO[7]]
                          [List_YOUSO[24]];

    double D_CntCoes[Matomnum+1]
                           [List_YOUSO[7]]
                           [List_YOUSO[24]];

    double D_CntCoes0[Matomnum+1]
                            [List_YOUSO[7]]
                            [List_YOUSO[24]];

    double DCntCoes_Spe[SpeciesNum+1]
                              [List_YOUSO[7]]
                              [List_YOUSO[24]];
  ********************************************************/

  Snd_CntCoes_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_CntCoes_Size = (int*)malloc(sizeof(int)*numprocs);

  CntCoes0 = (double***)malloc(sizeof(double**)*(Matomnum+MatomnumF+1));
  for (i=0; i<=(Matomnum+MatomnumF); i++){
    CntCoes0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      CntCoes0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  D_CntCoes = (double***)malloc(sizeof(double**)*(Matomnum+1));
  for (i=0; i<=Matomnum; i++){
    D_CntCoes[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      D_CntCoes[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  D_CntCoes0 = (double***)malloc(sizeof(double**)*(Matomnum+1));
  for (i=0; i<=Matomnum; i++){
    D_CntCoes0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      D_CntCoes0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  DCntCoes_Spe = (double***)malloc(sizeof(double**)*(SpeciesNum+1));
  for (i=0; i<=SpeciesNum; i++){
    DCntCoes_Spe[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      DCntCoes_Spe[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  DCntCoes_Spe0 = (double***)malloc(sizeof(double**)*(SpeciesNum+1));
  for (i=0; i<=SpeciesNum; i++){
    DCntCoes_Spe0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      DCntCoes_Spe0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  /********************************************************
    PrintMemory
  ********************************************************/

  if (firsttime) {
    PrintMemory("Opt_Contraction: CntCoes0",sizeof(double)*
           (Matomnum+1)*List_YOUSO[7]*List_YOUSO[24], NULL);
    PrintMemory("Opt_Contraction: D_CntCoes",sizeof(double)*
           (Matomnum+1)*List_YOUSO[7]*List_YOUSO[24], NULL);
    PrintMemory("Opt_Contraction: D_CntCoes0",sizeof(double)*
           (Matomnum+1)*List_YOUSO[7]*List_YOUSO[24], NULL);
    firsttime=0;
  }

  /********************************************************
    start calc.
  ********************************************************/

  norm_by_alpha = 0.01;

  /********************************************************
    unrestricted contraction 
  ********************************************************/

  if (RCnt_switch==0){

    D_UCntCoes(CntCoes,D_CntCoes,H,OLP,CDM,EDM);

    if (firsttime){
      norm_deri = NormD_UCnt(D_CntCoes);
      Cnt_scaling = norm_by_alpha/norm_deri;
    }

    al[0] = 0.5*Cnt_scaling;
    al[1] = 1.0*Cnt_scaling;
    al[2] = 3.0*Cnt_scaling;

    for(i=0; i<=2; i++){
      Change_CntCoes(al[i],OLP,CntCoes,D_CntCoes,CntCoes0);
      D_UCntCoes(CntCoes0,D_CntCoes0,H,OLP,CDM,EDM);
      nd[i] = NormD_UCnt(D_CntCoes0);
      if (2<=level_stdout){
        printf("<Opt_Contraction> i=%i al=%15.12f  nd=%15.12f\n",
               i,al[i],nd[i]);
      }
    }
  }

  /********************************************************
    restricted contraction 
  ********************************************************/

  else if (Cnt_switch==1 && ACnt_switch==0){
    D_RCntCoes(CntCoes,D_CntCoes,H,OLP,CDM,EDM);

    if (firsttime){
      norm_deri = NormD_ACnt(D_CntCoes);
      Cnt_scaling = norm_by_alpha/norm_deri;
    }

    al[0] = 0.5*Cnt_scaling;
    al[1] = 1.0*Cnt_scaling;
    al[2] = 3.0*Cnt_scaling;

    for(i=0; i<=2; i++){
      Change_CntCoes(al[i],OLP,CntCoes,D_CntCoes,CntCoes0);

      D_RCntCoes(CntCoes0,D_CntCoes0,H,OLP,CDM,EDM);

      nd[i] = NormD_RCnt(D_CntCoes0);
      if (2<=level_stdout){
        printf("<Opt_Contraction> i=%i al=%15.12f  nd=%15.12f\n",
               i,al[i],nd[i]);
      }
    }
  }

  /********************************************************
    restricted spcies contraction 
  ********************************************************/

  else if (Cnt_switch==1 && ACnt_switch==1){
    D_ACntCoes(CntCoes, D_CntCoes, DCntCoes_Spe, H, OLP, CDM, EDM);

    if (firsttime){
      norm_deri = NormD_ACnt(DCntCoes_Spe);
      Cnt_scaling = norm_by_alpha/norm_deri;
    }

    al[0] = 0.5*Cnt_scaling;
    al[1] = 1.0*Cnt_scaling;
    al[2] = 3.0*Cnt_scaling;

    for(i=0; i<=2; i++){
      Change_CntCoes(al[i],OLP,CntCoes,D_CntCoes,CntCoes0);

      D_ACntCoes(CntCoes0, D_CntCoes0, DCntCoes_Spe0, H, OLP, CDM, EDM);

      nd[i] = NormD_ACnt(DCntCoes_Spe0);

      if (2<=level_stdout){
        printf("<Opt_Contraction> i=%i al=%15.12f  nd=%15.12f\n",
               i,al[i],nd[i]);
      }
    }
  }

  x01 = al[0] - al[1];
  x12 = al[1] - al[2];
  x20 = al[2] - al[0];
  y12 = nd[1] - nd[2]; 
  y01 = nd[0] - nd[1];

  a = y12/(x12*x20) - y01/(x01*x20);
  b = y01/x01 - a*(al[0] + al[1]);
  c = nd[0] - a*al[0]*al[0] - b*al[0];
  alpha = -0.5*b/largest(a,10e-13);

  if (al[2]<alpha) alpha = al[2];
  if (alpha<0.0)   alpha = -0.1*alpha;

  Cnt_scaling = alpha;  

  Change_CntCoes(alpha,OLP,CntCoes,D_CntCoes,CntCoes0);

  if (RCnt_switch==0){
    D_UCntCoes(CntCoes0,D_CntCoes0,H,OLP,CDM,EDM);
    norm_deri = NormD_UCnt(D_CntCoes0);
  }
  else if (Cnt_switch==1 && ACnt_switch==0){
    D_RCntCoes(CntCoes0,D_CntCoes0,H,OLP,CDM,EDM);
    norm_deri = NormD_RCnt(D_CntCoes0);
  }
  else if (Cnt_switch==1 && ACnt_switch==1){
    D_ACntCoes(CntCoes0,D_CntCoes0,DCntCoes_Spe0,H,OLP,CDM,EDM);
    norm_deri = NormD_ACnt(DCntCoes_Spe0);
  }

  if (norm_deri<Oopt_NormD[1]){

    Change_CntCoes(alpha,OLP,CntCoes,D_CntCoes,CntCoes);
    Oopt_NormD[1] = norm_deri;  

    if (myid==Host_ID){
      printf("<Opt_Contraction> Scaling factor      = %15.12f\n",alpha);
      printf("<Opt_Contraction> Norm of derivatives = %15.12f \n",norm_deri);
    }
  }
  else{

    /*
    if (myid==Host_ID){
      printf("<Opt_Contraction> Could not optimize due to increasing the norm\n");
      printf("<Opt_Contraction> Scaling factor      = %15.12f\n",alpha);
      printf("<Opt_Contraction> Norm of derivatives = %15.12f \n",norm_deri);
    }
    */

    Change_CntCoes(alpha,OLP,CntCoes,D_CntCoes,CntCoes);
    Oopt_NormD[1] = norm_deri;  

    if (myid==Host_ID){
      printf("<Opt_Contraction> Scaling factor      = %15.12f\n",alpha);
      printf("<Opt_Contraction> Norm of derivatives = %15.12f \n",norm_deri);
    }

    /* reduce Cnt_scaling */

    Cnt_scaling = 0.7*Cnt_scaling;
  }

  if (2<=level_stdout){
    printf("<Opt_Contraction> Contraction coefficients\n");
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      for (be=0; be<Spe_Total_CNO[Cwan]; be++){
        for (p=0; p<Spe_Specified_Num[Cwan][be]; p++){
          printf("  Mc_AN=%2d Gc_AN=%2d be=%2d p=%2d  %15.12f\n",
                    Mc_AN,Gc_AN,be,p,CntCoes[Mc_AN][be][p]);
        }
      }
    }
  }

  /****************************************************
    MPI:

    CntCoes
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
      if (F_Snd_Num[IDS]!=0){

        size1 = 0;
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Mc_AN = Snd_MAN[IDS][n];
          Gc_AN = Snd_GAN[IDS][n];
          wan = WhatSpecies[Gc_AN]; 
          for (be=0; be<Spe_Total_CNO[wan]; be++){
            size1 += Spe_Specified_Num[wan][be];
	  }
	}

        Snd_CntCoes_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_CntCoes_Size[IDS] = 0;
      }

      /* receiving of size of data */

      if (F_Rcv_Num[IDR]!=0){
        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_CntCoes_Size[IDR] = size2;
      }
      else{
        Rcv_CntCoes_Size[IDR] = 0;
      }
    
      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);

    }
    else {
      Snd_CntCoes_Size[myid] = 0;
      Rcv_CntCoes_Size[myid] = 0;
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

      if (F_Snd_Num[IDS]!=0){

        size1 = Snd_CntCoes_Size[IDS];

        /* allocation of array */
        tmp_array = (double*)malloc(sizeof(double)*size1);

        /* multidimentional array to vector array */

        num = 0;
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Mc_AN = Snd_MAN[IDS][n];
          Gc_AN = Snd_GAN[IDS][n];
          wan = WhatSpecies[Gc_AN]; 
          for (be=0; be<Spe_Total_CNO[wan]; be++){
            for (p=0; p<Spe_Specified_Num[wan][be]; p++){
              tmp_array[num] = CntCoes[Mc_AN][be][p];
              num++;
  	    }
	  }
        }

        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if (F_Rcv_Num[IDR]!=0){

        size2 = Rcv_CntCoes_Size[IDR];

        /* allocation of array */
        tmp_array2 = (double*)malloc(sizeof(double)*size2);

        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        num = 0;
        Mc_AN = F_TopMAN[IDR] - 1;
        for (n=0; n<F_Rcv_Num[IDR]; n++){
          Mc_AN++;
          Gc_AN = Rcv_GAN[IDR][n];
          wan = WhatSpecies[Gc_AN];
          for (be=0; be<Spe_Total_CNO[wan]; be++){
            for (p=0; p<Spe_Specified_Num[wan][be]; p++){
              CntCoes[Mc_AN][be][p] = tmp_array2[num];
              num++;
	    }
  	  }
        }

        /* freeing of array */
        free(tmp_array2);
      }

      if (F_Snd_Num[IDS]!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array); /* freeing of array */
      }
    }
  }

  /********************************************************
    freeing of arrays:

    int Snd_CntCoes_Size[numprocs];
    int Rcv_CntCoes_Size[numprocs];

    double CntCoes0[Matomnum+MatomnumF+1]
                          [List_YOUSO[7]]
                          [List_YOUSO[24]];

    double D_CntCoes[Matomnum+1]
                           [List_YOUSO[7]]
                           [List_YOUSO[24]];

    double D_CntCoes0[Matomnum+1]
                            [List_YOUSO[7]]
                            [List_YOUSO[24]];
  ********************************************************/

  free(Snd_CntCoes_Size);
  free(Rcv_CntCoes_Size);

  for (i=0; i<=(Matomnum+MatomnumF); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(CntCoes0[i][j]);
    }
    free(CntCoes0[i]);
  }
  free(CntCoes0);

  for (i=0; i<=Matomnum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(D_CntCoes[i][j]);
    }
    free(D_CntCoes[i]);
  }
  free(D_CntCoes);

  for (i=0; i<=Matomnum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(D_CntCoes0[i][j]);
    }
    free(D_CntCoes0[i]);
  }
  free(D_CntCoes0);

  for (i=0; i<=SpeciesNum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(DCntCoes_Spe[i][j]);
    }
    free(DCntCoes_Spe[i]);
  }
  free(DCntCoes_Spe);

  for (i=0; i<=SpeciesNum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(DCntCoes_Spe0[i][j]);
    }
    free(DCntCoes_Spe0[i]);
  }
  free(DCntCoes_Spe0);

  if (firsttime) {
    firsttime=0;
  }

  /* return */
  return Oopt_NormD[1];
}

void D_UCntCoes(double ***CntCoes0,
                double ***D_CntCoes,
                double *****H, double *****OLP,
                double *****CDM,
                double *****EDM)
{
  int Cwan,Mc_AN,Gc_AN,Mh_AN,p,q,al,be,p0,q0;
  int h_AN,Gh_AN,Hwan,spin;
  double sum,sum0;

  /****************************************************
            Calculate D_CntCoes[Mc_AN][al][p]
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
      
    for (al=0; al<Spe_Total_CNO[Cwan]; al++){
      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	p0 = Spe_Trans_Orbital[Cwan][al][p];

	sum = 0.0;
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          Mh_AN = F_G2M[Gh_AN];
	  for (be=0; be<Spe_Total_CNO[Hwan]; be++){
	    for (q=0; q<Spe_Specified_Num[Hwan][be]; q++){
	      q0 = Spe_Trans_Orbital[Hwan][be][q];

	      sum0 = 0.0;
	      for (spin=0; spin<=SpinP_switch; spin++){
		sum0 = sum0 +
		   CDM[spin][Mc_AN][h_AN][al][be]*H[spin][Mc_AN][h_AN][p0][q0]
		  -EDM[spin][Mc_AN][h_AN][al][be]*OLP[0][Mc_AN][h_AN][p0][q0];
	      } 
	      if (SpinP_switch==0) sum0 = 2.0*sum0; 

	      sum = sum + 2.0*sum0*CntCoes0[Mh_AN][be][q];
        
	    }
	  }
	}

	D_CntCoes[Mc_AN][al][p] = sum;

	/*
        printf("step=%2d Mc_AN=%2d Gc_AN=%2d al=%2d p=%2d  %15.12f\n",
 	        step,Mc_AN,Gc_AN,al,p,sum);
	*/

      }
    }
  }

}




void Change_CntCoes(double alpha,
                    double *****OLP,
                    double ***CntCoes,
                    double ***D_CntCoes,
                    double ***CntCoes0)
{
  int Cwan,Mc_AN,Gc_AN;
  int p,q,al,p0,q0;
  int n,num,size1,size2,wan,be;
  double sum,tmp0;
  double *tmp_array;
  double *tmp_array2;
  int *Snd_CntCoes_Size;
  int *Rcv_CntCoes_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /********************************************************
    allocation of arrays:

    int Snd_CntCoes_Size[numprocs];
    int Rcv_CntCoes_Size[numprocs];
  ********************************************************/

  Snd_CntCoes_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_CntCoes_Size = (int*)malloc(sizeof(int)*numprocs);

  /****************************************************
                     Change CntCoes
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    for (al=0; al<Spe_Total_CNO[Cwan]; al++){
      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	 CntCoes0[Mc_AN][al][p] = CntCoes[Mc_AN][al][p]
                                  -alpha*D_CntCoes[Mc_AN][al][p];
      }
    }
  }

  /****************************************************
                     Normalization
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    for (al=0; al<Spe_Total_CNO[Cwan]; al++){

      sum = 0.0;
      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	p0 = Spe_Trans_Orbital[Cwan][al][p];

	for (q=0; q<Spe_Specified_Num[Cwan][al]; q++){
          q0 = Spe_Trans_Orbital[Cwan][al][q];

          tmp0 = CntCoes0[Mc_AN][al][p]*CntCoes0[Mc_AN][al][q];
          sum += tmp0*OLP[0][Mc_AN][0][p0][q0]; 
        }
      }

      tmp0 = 1.0/sqrt(sum);
      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
        CntCoes0[Mc_AN][al][p] = CntCoes0[Mc_AN][al][p]*tmp0;
      }        

    }
  }

  /****************************************************
    MPI:

    CntCoes0
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
      if (F_Snd_Num[IDS]!=0){

        size1 = 0;
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Mc_AN = Snd_MAN[IDS][n];
          Gc_AN = Snd_GAN[IDS][n];
          wan = WhatSpecies[Gc_AN]; 
          for (be=0; be<Spe_Total_CNO[wan]; be++){
            size1 += Spe_Specified_Num[wan][be];
	  }
	}

        Snd_CntCoes_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_CntCoes_Size[IDS] = 0;
      }

      /* receiving of size of data */

      if (F_Rcv_Num[IDR]!=0){
        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_CntCoes_Size[IDR] = size2;
      }
      else{
        Rcv_CntCoes_Size[IDR] = 0;
      }
    
      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);

    }
    else {
      Snd_CntCoes_Size[myid] = 0;
      Rcv_CntCoes_Size[myid] = 0;
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

      if (F_Snd_Num[IDS]!=0){

        size1 = Snd_CntCoes_Size[IDS];

        /* allocation of array */
        tmp_array = (double*)malloc(sizeof(double)*size1);

        /* multidimentional array to vector array */

        num = 0;
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Mc_AN = Snd_MAN[IDS][n];
          Gc_AN = Snd_GAN[IDS][n];
          wan = WhatSpecies[Gc_AN]; 
          for (be=0; be<Spe_Total_CNO[wan]; be++){
            for (p=0; p<Spe_Specified_Num[wan][be]; p++){
              tmp_array[num] = CntCoes0[Mc_AN][be][p];
              num++;
  	    }
	  }
        }

        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if (F_Rcv_Num[IDR]!=0){

        size2 = Rcv_CntCoes_Size[IDR];

        /* allocation of array */
        tmp_array2 = (double*)malloc(sizeof(double)*size2);

        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        num = 0;
        Mc_AN = F_TopMAN[IDR] - 1;
        for (n=0; n<F_Rcv_Num[IDR]; n++){
          Mc_AN++;
          Gc_AN = Rcv_GAN[IDR][n];
          wan = WhatSpecies[Gc_AN];
          for (be=0; be<Spe_Total_CNO[wan]; be++){
            for (p=0; p<Spe_Specified_Num[wan][be]; p++){
              CntCoes0[Mc_AN][be][p] = tmp_array2[num];
              num++;
	    }
  	  }
        }

        /* freeing of array */
        free(tmp_array2);
      }

      if (F_Snd_Num[IDS]!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array);  /* freeing of array */ 
      }
    }
  }

  /********************************************************
    allocation of arrays:

    int Snd_CntCoes_Size[numprocs];
    int Rcv_CntCoes_Size[numprocs];
  ********************************************************/

  free(Snd_CntCoes_Size);
  free(Rcv_CntCoes_Size);

}

void D_RCntCoes(double ***CntCoes0,
                double ***D_CntCoes,
                double *****H, double *****OLP,
                double *****CDM,
                double *****EDM)
{
  static int firsttime=1;
  int i,j,Cwan,Mc_AN,Gc_AN,Mh_AN;
  int p,q,al,al0,be,p0,q0;
  int h_AN,Gh_AN,Hwan,spin;
  int L0,Mul0,M0;
  double sum,sum0;
  double *TmpD;
  double ***D0_CntCoes;

  /********************************************************
    allocation of arrays:

    double TmpD[List_YOUSO[24]];

    double D0_CntCoes[Matomnum+1]
                            [List_YOUSO[7]]
                            [List_YOUSO[24]];

  ********************************************************/

  TmpD = (double*)malloc(sizeof(double)*List_YOUSO[24]);

  D0_CntCoes = (double***)malloc(sizeof(double**)*(Matomnum+1));
  for (i=0; i<=Matomnum; i++){
    D0_CntCoes[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      D0_CntCoes[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  /* PrintMemory */

  if (firsttime) {
    PrintMemory("D_RCntCoes: D0_CntCoes",sizeof(D0_CntCoes),NULL);
    firsttime=0;
  }

  /****************************************************
            Calculate D_CntCoes[Mc_AN][al][p]
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    for (al=0; al<Spe_Total_CNO[Cwan]; al++){
      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	p0 = Spe_Trans_Orbital[Cwan][al][p];

	sum = 0.0;
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          Mh_AN = F_G2M[Gh_AN];
	  for (be=0; be<Spe_Total_CNO[Hwan]; be++){
	    for (q=0; q<Spe_Specified_Num[Hwan][be]; q++){
	      q0 = Spe_Trans_Orbital[Hwan][be][q];

	      sum0 = 0.0;
	      for (spin=0; spin<=SpinP_switch; spin++){
		sum0 = sum0 +
		   CDM[spin][Mc_AN][h_AN][al][be]*H[spin][Mc_AN][h_AN][p0][q0]
		  -EDM[spin][Mc_AN][h_AN][al][be]*OLP[0][Mc_AN][h_AN][p0][q0];
	      } 
	      if (SpinP_switch==0) sum0 = 2.0*sum0; 

	      sum = sum + 2.0*sum0*CntCoes0[Mh_AN][be][q];
        
	    }
	  }
	}

	D0_CntCoes[Mc_AN][al][p] = sum;

	/*
          printf("step=%2d Mc_AN=%2d Gc_AN=%2d al=%2d p=%2d  %15.12f\n",
	  step,Mc_AN,Gc_AN,al,p,sum);
	*/

      }
    }

    /****************************************************
          taking into account of the restriction
    ****************************************************/
          
    al = -1;
    for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){

	for (p=0; p<Spe_Specified_Num[Cwan][al+1]; p++){
	  TmpD[p] = 0.0;
	}

	al0 = al;
	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    TmpD[p] = TmpD[p] + D0_CntCoes[Mc_AN][al][p];
	  }
	}
	al = al0;

	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    D_CntCoes[Mc_AN][al][p] = TmpD[p];
	  }
	}
      }
    }      
  }

  /********************************************************
    freeing of arrays:

    double TmpD[List_YOUSO[24]];

    double D0_CntCoes[Matomnum+1]
                            [List_YOUSO[7]]
                            [List_YOUSO[24]];
  ********************************************************/

  free(TmpD);

  for (i=0; i<=Matomnum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(D0_CntCoes[i][j]);
    }
    free(D0_CntCoes[i]);
  }
  free(D0_CntCoes);

}




void D_ACntCoes(double ***CntCoes0,
                double ***D_CntCoes,
                double ***DCntCoes_Spe,
                double *****H, double *****OLP,
                double *****CDM,
                double *****EDM)
{
  static int firsttime=1;
  int i,j,Cwan,Mc_AN,Gc_AN,Mh_AN;
  int p,q,al,al0,be,p0,q0;
  int h_AN,Gh_AN,Hwan,spin;
  int L0,Mul0,M0,wan;
  double sum,sum0;
  double *TmpD;
  double ***D0_CntCoes;
  double ***My_DCntCoes_Spe;

  /********************************************************
    allocation of arrays:

    double TmpD[List_YOUSO[24]];

    double D0_CntCoes[Matomnum+1]
                            [List_YOUSO[7]]
                            [List_YOUSO[24]];


    double My_DCntCoes_Spe[SpeciesNum+1]
                                 [List_YOUSO[7]]
                                 [List_YOUSO[24]];
  ********************************************************/

  TmpD = (double*)malloc(sizeof(double)*List_YOUSO[24]);

  D0_CntCoes = (double***)malloc(sizeof(double**)*(Matomnum+1));
  for (i=0; i<=Matomnum; i++){
    D0_CntCoes[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      D0_CntCoes[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  My_DCntCoes_Spe = (double***)malloc(sizeof(double**)*(SpeciesNum+1));
  for (i=0; i<=SpeciesNum; i++){
    My_DCntCoes_Spe[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      My_DCntCoes_Spe[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
    }
  }

  /* PrintMemory */

  if (firsttime) {
    PrintMemory("D_RCntCoes: D0_CntCoes",sizeof(D0_CntCoes),NULL);
    firsttime=0;
  }

  /****************************************************
            Calculate D_CntCoes[Mc_AN][al][p]
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    for (al=0; al<Spe_Total_CNO[Cwan]; al++){
      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	p0 = Spe_Trans_Orbital[Cwan][al][p];

	sum = 0.0;
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          Mh_AN = F_G2M[Gh_AN];
	  for (be=0; be<Spe_Total_CNO[Hwan]; be++){
	    for (q=0; q<Spe_Specified_Num[Hwan][be]; q++){
	      q0 = Spe_Trans_Orbital[Hwan][be][q];

	      sum0 = 0.0;
	      for (spin=0; spin<=SpinP_switch; spin++){
		sum0 = sum0 +
		   CDM[spin][Mc_AN][h_AN][al][be]*H[spin][Mc_AN][h_AN][p0][q0]
		  -EDM[spin][Mc_AN][h_AN][al][be]*OLP[0][Mc_AN][h_AN][p0][q0];
	      } 
	      if (SpinP_switch==0) sum0 = 2.0*sum0; 

	      sum = sum + 2.0*sum0*CntCoes0[Mh_AN][be][q];
        
	    }
	  }
	}

	D0_CntCoes[Mc_AN][al][p] = sum;

	/*
          printf("step=%2d Mc_AN=%2d Gc_AN=%2d al=%2d p=%2d  %15.12f\n",
	  step,Mc_AN,Gc_AN,al,p,sum);
	*/

      }
    }

    /****************************************************
          taking into account of the restriction
    ****************************************************/
          
    al = -1;
    for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){

	for (p=0; p<Spe_Specified_Num[Cwan][al+1]; p++){
	  TmpD[p] = 0.0;
	}

	al0 = al;
	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    TmpD[p] = TmpD[p] + D0_CntCoes[Mc_AN][al][p];
	  }
	}
	al = al0;

	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    D_CntCoes[Mc_AN][al][p] = TmpD[p];
	  }
	}
      }
    }      
  }

 /********************************************************
    calculate derivatives with respect to CntCoes_Spe
 ********************************************************/

  for (wan=0; wan<SpeciesNum; wan++){
    for (i=0; i<List_YOUSO[7]; i++){
      for (j=0; j<List_YOUSO[24]; j++){
        My_DCntCoes_Spe[wan][i][j] = 0.0;
      }
    }
  }

  /* local sum in a proccessor */
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    wan = WhatSpecies[Gc_AN];
    for (al=0; al<Spe_Total_CNO[wan]; al++){
      for (p=0; p<Spe_Specified_Num[wan][al]; p++){
        My_DCntCoes_Spe[wan][al][p] += D_CntCoes[Mc_AN][al][p];
      }        
    }
  }

  /* global sum by MPI */
  for (wan=0; wan<SpeciesNum; wan++){
    for (al=0; al<Spe_Total_CNO[wan]; al++){
      for (p=0; p<Spe_Specified_Num[wan][al]; p++){
        MPI_Allreduce(&My_DCntCoes_Spe[wan][al][p], &DCntCoes_Spe[wan][al][p],
                       1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
      }
    }
  }    

  /* DCntCoes_Spe to DCntCoes */
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    wan = WhatSpecies[Gc_AN];
    for (al=0; al<Spe_Total_CNO[wan]; al++){
      for (p=0; p<Spe_Specified_Num[wan][al]; p++){
        D_CntCoes[Mc_AN][al][p] = DCntCoes_Spe[wan][al][p];
      }        
    }
  }

  /********************************************************
    freeing of arrays:

    double TmpD[List_YOUSO[24]];

    double D0_CntCoes[Matomnum+1]
                            [List_YOUSO[7]]
                            [List_YOUSO[24]];
  ********************************************************/

  free(TmpD);

  for (i=0; i<=Matomnum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(D0_CntCoes[i][j]);
    }
    free(D0_CntCoes[i]);
  }
  free(D0_CntCoes);

  for (i=0; i<=SpeciesNum; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(My_DCntCoes_Spe[i][j]);
    }
    free(My_DCntCoes_Spe[i]);
  }
  free(My_DCntCoes_Spe);
}









double NormD_UCnt(double ***D_CntCoes)
{
  int Cwan,Mc_AN,Gc_AN,p,q,al;
  int L0,Mul0,M0;
  double My_NormD,NormD;

  /****************************************************
                 Calculate NormD_UCnt
  ****************************************************/

  My_NormD = 0.0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    al = -1;
    for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){
	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    My_NormD = My_NormD + D_CntCoes[Mc_AN][al][p]*
                                  D_CntCoes[Mc_AN][al][p];  
	  }
	}
      }
    }      
  }

  /* MPI My_NormD */
  MPI_Allreduce(&My_NormD, &NormD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  return NormD;
}


double NormD_RCnt(double ***D_CntCoes)
{
  int Cwan,Mc_AN,Gc_AN,p,q,al;
  int L0,Mul0,M0;
  double My_NormD,NormD;

  /****************************************************
                 Calculate NormD_RCnt
  ****************************************************/

  My_NormD = 0.0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    al = -1;
    for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){
	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    if (M0==0) My_NormD = My_NormD + D_CntCoes[Mc_AN][al][p]*
                                             D_CntCoes[Mc_AN][al][p];  
	  }
	}
      }
    }      
  }

  /* MPI My_NormD */
  MPI_Allreduce(&My_NormD, &NormD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  return NormD;
}


double NormD_ACnt(double ***DCntCoes_Spe)
{
  int Cwan,Mc_AN,Gc_AN,p,q,al;
  int L0,Mul0,M0;
  double My_NormD,NormD;

  /****************************************************
                 Calculate NormD_ACnt
  ****************************************************/

  My_NormD = 0.0;

  for (Cwan=0; Cwan<SpeciesNum; Cwan++){
    al = -1;
    for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){
	for (M0=0; M0<=2*L0; M0++){
	  al++;
	  for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	    if (M0==0) My_NormD = My_NormD + DCntCoes_Spe[Cwan][al][p]*
                                             DCntCoes_Spe[Cwan][al][p];  
	  }
	}
      }
    }      
  }

  NormD = My_NormD;

  return NormD;
}

