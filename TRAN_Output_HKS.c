#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_variables.h"
#include "tran_prototypes.h" 


/****************************************************************
   purpose: save and load data to restart without input file 

   from RestartFileDFT.c
------------------------------------------

"read" mode allocates 

int *WhatSpecies_l;
int *Spe_Total_CNO_l;
int *Spe_Total_NO_l;
int *FNAN_l;
int **natn_l;
int **ncn_l;
int **atv_ijk_l;

double *****OLP_l;
double *****H_l;
double ******DM_l;

or 

int *WhatSpecies_r;
int *Spe_Total_CNO_r;
int *Spe_Total_NO_r;
int *FNAN_r;
int **natn_r;
int **ncn_r;
int **atv_ijk_r;

double *****OLP_r;
double *****H_r;
double ******DM_r;


*****************************************************************/

/*
  e.g. Overlap_Band, which gives hints of data to be saved
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
  GA_AN = M2G[MA_AN];  
  wanA = WhatSpecies[GA_AN]; int* WhatSpecies[atomnum+1]
  tnoA = Spe_Total_CNO[wanA]; int* Spe_Total_CNO[SpeciesNum]
  Anum = MP[GA_AN];   int *MP, neglect!

  for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
  GB_AN = natn[GA_AN][LB_AN]; int** natn[atomnum+1][Max_FSNAN*ScaleSize+1]
  Rn = ncn[GA_AN][LB_AN];  int** ncn[atomnum+1][Max_FSNAN*ScaleSize+1]
  wanB = WhatSpecies[GB_AN]; 
  tnoB = Spe_Total_CNO[wanB];

  l1 = atv_ijk[Rn][1];  int** atv_ijk[TCpyCell+1][4];
  l2 = atv_ijk[Rn][2];
  l3 = atv_ijk[Rn][3];
  kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

  si = sin(2.0*PI*kRn);
  co = cos(2.0*PI*kRn);
  Bnum = MP[GB_AN];
  for (i=0; i<tnoA; i++){
  for (j=0; j<tnoB; j++){
  s = OLP[MA_AN][LB_AN][i][j]; double****
  size: OLP[4]
  [Matomnum+MatomnumF+MatomnumS+1]
  [FNAN[Gc_AN]+1]
  [Spe_Total_NO[Cwan]]
  [Spe_Total_NO[Hwan]]

  int *Spe_Total_NO  Spe_Total_NO[SpeciesNum]
*/


/***************************************************************************/

int TRAN_Output_HKS(char *fileHKS)
{
  FILE *fp;
  int i_vec[100],i2_vec[2];
  double d_vec[100];
  int i,id,j;
  int size1,size,vsize;
  int *ia_vec;
  int Gc_AN,Mc_AN, tno0, Cwan, h_AN, tno1, Gh_AN, Hwan,k,m,N;

  int myid,numprocs;
  int tag=99;
  double *v;

  MPI_Status status;
  MPI_Request request;

  MPI_Comm_rank(mpi_comm_level1, &myid);
  MPI_Comm_size(mpi_comm_level1, &numprocs);
   
  if (myid==Host_ID) {

    /* make a filename */
    if ( (fp=fopen(fileHKS,"w"))==NULL) {
      printf("can not open %s\n",fileHKS);
      exit(0);
    }

    /* save data to the file (*fp) */

    /* parameter to allocate memory */
    i=0;
    i_vec[i++]=1; /* major version */
    i_vec[i++]=0; /* minor version*/

    fwrite(i_vec,sizeof(int),i,fp);

    i=0;
    i_vec[i++]= SpinP_switch;
    i_vec[i++]= atomnum;
    i_vec[i++]= SpeciesNum;
    i_vec[i++]= Max_FSNAN;
    i_vec[i++]= TCpyCell;
    i_vec[i++]= Matomnum;
    i_vec[i++]= MatomnumF;
    i_vec[i++]= MatomnumS;
    i_vec[i++]= Ngrid1;
    i_vec[i++]= Ngrid2;
    i_vec[i++]= Ngrid3;
    i_vec[i++]= Num_Cells0; 


    fwrite(i_vec,sizeof(int),i,fp);


#ifdef DEBUG
    printf("Scale =%lf\n",ScaleSize);
    printf("%d %d %d %d %d %d\n",SpinP_switch,atomnum,SpeciesNum,
	   Max_FSNAN, TCpyCell, Matomnum );
    printf("%d %d %d %d %d %d\n",MatomnumF,MatomnumS,Ngrid1,
	   Ngrid2,Ngrid3,Num_Cells0);

#endif

    id=0;
    d_vec[id++]= ScaleSize;
    for(j=1;j<=3;j++) {
      for (k=1;k<=3;k++) {
	d_vec[id++]= tv[j][k];
      }
    }
    for(j=1;j<=3;j++) {
      for (k=1;k<=3;k++) {
	d_vec[id++]= gtv[j][k];
      }
    }
    for (i=1;i<=3;i++) {
      d_vec[id++] = Grid_Origin[i];
    }
    printf("Grid_Origin=%lf %lf %lf\n",Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);

    for (i=1;i<=3;i++) {
      d_vec[id++] = Gxyz[1][i];
    }
    d_vec[id++] = ChemP; 

    fwrite(d_vec,sizeof(double),id,fp);




    /*  data in arrays */
    fwrite(WhatSpecies, sizeof(int), atomnum+1, fp);
    fwrite(Spe_Total_CNO, sizeof(int), SpeciesNum, fp);
    fwrite(Spe_Total_NO,sizeof(int),SpeciesNum,fp);

#if 0
    printf("spe %d %d %d\n",WhatSpecies[1],Spe_Total_CNO[1],Spe_Total_NO[1]);
#endif

    fwrite(FNAN,sizeof(int),atomnum+1,fp);

    size1=(int)Max_FSNAN*ScaleSize+1;
    for (i=0;i<= atomnum; i++) {
      fwrite(natn[i],sizeof(int),size1,fp);
    }
    for (i=0;i<= atomnum; i++) {
      fwrite(ncn[i],sizeof(int),size1,fp);
    }
#if 0
    printf("ncn %d %d \n",natn[1][1], ncn[1][1]);
    printf("ncn\n");
    for (i=0;i<=atomnum;i++) {
      for (j=0;j<size1;j++) {
	printf("%d ",ncn[i][j]);
      }
      printf("\n");
    }

#endif

 
    /*  printf("atv_ijk\n"); */
    size1=(TCpyCell+1)*4;
    ia_vec=(int*)malloc(sizeof(int)*size1);
    id=0;
    for (i=0;i<TCpyCell+1;i++) {
      for (j=0;j<=3;j++) {
	ia_vec[id++]=atv_ijk[i][j];
	/*      printf("%d ",atv_ijk[i][j]); */
      }
      /*    printf("\n"); */

    }
    fwrite(ia_vec,sizeof(int),size1,fp);
    free(ia_vec);

  }  /* if (myid==Host_ID) */


  /* allocate v */

  v = (double*)malloc(sizeof(double)*List_YOUSO[8]*List_YOUSO[7]*List_YOUSO[7]);

  /* OLP,  this is complex */

  for (k=0;k<4;k++) {

    int m,ID;
    /*global  Gc_AN  1:atomnum */
    /*variable ID = G2ID[Gc_AN] */
    /*variable Mc_AN = G2M[Gc_AN] */
    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      Mc_AN = F_G2M[Gc_AN];
      ID = G2ID[Gc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_NO[Cwan];

      /* OLP into v */

      if (myid==ID) {

	vsize = 0;

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    for (j=0;j<tno1;j++) {
	      v[vsize] =  OLP[k][Mc_AN][h_AN][i][j];
	      vsize++;
	    }
	  }
	} 

        /* Isend */

        if (myid!=Host_ID){
          MPI_Isend(&vsize, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&status);
          MPI_Isend(&v[0], vsize, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&status);
	}
        else{
          fwrite(v, sizeof(double), vsize, fp);
        }
      }

      /* Recv */

      else if (ID!=myid && myid==Host_ID){
        MPI_Recv(&vsize, 1, MPI_INT, ID, tag, mpi_comm_level1, &status);
        MPI_Recv(&v[0], vsize, MPI_DOUBLE, ID, tag, mpi_comm_level1, &status);
        fwrite(v, sizeof(double), vsize, fp);
      }

    } /* Gc_AN */
  } /* k */



  /* H */
  for (k=0; k<=SpinP_switch; k++){
    int ID;

    /*global  Gc_AN  1:atomnum */
    /*variable ID = G2ID[Gc_AN] */
    /*variable Mc_AN = G2M[Gc_AN] */
    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      ID    = G2ID[Gc_AN]; 
      Mc_AN = F_G2M[Gc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_NO[Cwan];  

      /* H into v */

      if (myid==ID) {

	vsize = 0;

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

          for (i=0; i<tno0; i++){
	    for (j=0;j<tno1; j++) {
	      v[vsize] = H[k][Mc_AN][h_AN][i][j];
	      vsize++;
	    }
	  }
	}

        /* Isend */

        if (myid!=Host_ID){
          MPI_Isend(&vsize, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&status);
          MPI_Isend(&v[0], vsize, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&status);
	}
        else{
          fwrite(v, sizeof(double), vsize, fp);
        }
      }

      /* Recv */

      else if (ID!=myid && myid==Host_ID){
        MPI_Recv(&vsize, 1, MPI_INT, ID, tag, mpi_comm_level1, &status);
        MPI_Recv(&v[0], vsize, MPI_DOUBLE, ID, tag, mpi_comm_level1, &status);
        fwrite(v, sizeof(double), vsize, fp);
      }

    }
  } /* k */



  /* DM */

  for (m=0; m<1; m++){
    int ID;

    for (k=0; k<=SpinP_switch; k++){
      /*global  Gc_AN  1:atomnum */
      /*variable ID = G2ID[Gc_AN] */
      /*variable Mc_AN = G2M[Gc_AN] */
      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

        ID = G2ID[ Gc_AN]; 
	Mc_AN = F_G2M[Gc_AN];
	Cwan = WhatSpecies[Gc_AN];
	tno0 = Spe_Total_NO[Cwan];  

        /* DM into v */

        if (myid==ID) {

    	  vsize = 0;

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      for (j=0;j<tno1;j++) {
		v[vsize] = DM[m][k][Mc_AN][h_AN][i][j] ;
  	        vsize++;
	      }
	    }
	  }

          /* Isend */

	  if (myid!=Host_ID){
	    MPI_Isend(&vsize, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
	    MPI_Wait(&request,&status);
	    MPI_Isend(&v[0], vsize, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
	    MPI_Wait(&request,&status);
	  }
	  else{
	    fwrite(v, sizeof(double), vsize, fp);
	  }
	}

        /* Recv */

        else if (ID!=myid && myid==Host_ID){
          MPI_Recv(&vsize, 1, MPI_INT, ID, tag, mpi_comm_level1, &status);
          MPI_Recv(&v[0], vsize, MPI_DOUBLE, ID, tag, mpi_comm_level1, &status);
          fwrite(v, sizeof(double), vsize, fp);
        }

      } /* Gc_AN */
    } /* k */
  } /* m */

  /* free v */

  free(v);

  N=Ngrid1*Ngrid2*Ngrid3;

  /* Density_Grid */
  for (m=0;m<=SpinP_switch;m++) {
#if 1
    TRAN_Output_HKS_Write_Grid(mpi_comm_level1, 
			       My_NGrid1_Poisson,Ngrid2,Ngrid3,
			       Num_Snd_Grid1,Snd_Grid1,Num_Rcv_Grid1,Rcv_Grid1,My_Cell0,Start_Grid1,End_Grid1,
			       Density_Grid[m],NULL,fp);
#else 
    for (i=0;i<Ngrid1;i++) {
      /* fwrite(Density_Grid[m],sizeof(double), N,fp); */
      fwrite(&Density_Grid[m][My_Cell0[i]*Ngrid2*Ngrid3],sizeof(double),Ngrid2*Ngrid3,fp);
    }
#endif
  } /* m */







  /* Vpot_Grid */
  /* fwrite(dVHart_Grid,sizeof(double), N,fp); */
#if 1
  TRAN_Output_HKS_Write_Grid(mpi_comm_level1, 
			     My_NGrid1_Poisson,Ngrid2,Ngrid3,
			     Num_Snd_Grid1,Snd_Grid1,Num_Rcv_Grid1,Rcv_Grid1,My_Cell0,Start_Grid1,End_Grid1,
			     dVHart_Grid,NULL,fp);
#else
  for (i=0;i<Ngrid1;i++) {
    fwrite(&dVHart_Grid[My_Cell0[i]*Ngrid2*Ngrid3],sizeof(double),Ngrid2*Ngrid3,fp);
  }
    
  /* debug */
  i = My_Cell0[Ngrid1-1];
  j = Ngrid2-1;
  k = Ngrid3-1;
  printf("the last of dVHart_Grid = %20.10le\n",dVHart_Grid[ i*Ngrid2*Ngrid3+ j*Ngrid3+ k ]); 
#endif


  if (myid==Host_ID) 
    fclose(fp);

  return 1;

}



	  
  
