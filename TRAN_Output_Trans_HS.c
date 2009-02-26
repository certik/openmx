#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"

#define print_stdout 0

void  TRAN_Output_Trans_HS(
        MPI_Comm comm1,
        int SpinP_switch, 
        double ChemP ,
        double *****H,
        double *****OLP,
        int atomnum,
        int SpeciesNum,
        int *WhatSpecies,
        int *Spe_Total_CNO,
        int *FNAN,
        int **natn,
        int **ncn,
        int *G2ID,
        int **atv_ijk,
        int Max_FSNAN,
        double ScaleSize,
        int *F_G2M,      
        int TCpyCell,
        int *List_YOUSO,
        char *filepath,
        char *filename,
        char *fname  )
{
  FILE *fp;
  double v[20];
  int iv[20];
  int *ia_vec;
  int i,j,k,id,vsize;
  int iside,spin,size1;
  int Gc_AN,Mc_AN,GA_AN,GB_AN,LB_AN;
  int h_AN,Gh_AN,Hwan,tno1,tno0,Cwan;
  int wanA,wanB,tnoA,tnoB; 
  double *v1;
  char name[100];
  int numprocs,myid,tag=999,ID;
  MPI_Status status;
  MPI_Request request;

  MPI_Comm_size(comm1 ,&numprocs);
  MPI_Comm_rank(comm1, &myid);

  if (myid==Host_ID){
    printf("<TRAN_Output_Trans_HS>\n");
  }

  if (myid==Host_ID){

    /* work only myid==0 */

    sprintf(name,"%s%s.%s",filepath,filename,fname);

    if ( (fp=fopen(name,"w"))==NULL ) {
      printf("can not open %s\n",name);
      printf("in TRAN_Output_Trans_HS\n");
      exit(0);
    }
  }

  if (myid==Host_ID) {

    /* SpinP_switch, NUM_c, and NUM_e */
   
    i=0;
    iv[i++]= SpinP_switch;
    iv[i++]= NUM_c;
    iv[i++]= NUM_e[0];
    iv[i++]= NUM_e[1];
    fwrite(iv,sizeof(int),i,fp);

    /* chemical potential */

    i=0;
    v[i++]= ChemP;
    v[i++]= ChemP_e[0];
    v[i++]= ChemP_e[1];
    fwrite(v,sizeof(double),i,fp);

    /* tran_bias_apply */

    iv[0] = tran_bias_apply;
    fwrite(iv,sizeof(int),1,fp);

    /* the number of atoms */

    i=0;
    iv[i++]= atomnum;
    iv[i++]= atomnum_e[0];
    iv[i++]= atomnum_e[1];
    fwrite(iv,sizeof(int),i,fp);

    /* the number of species */

    i=0;
    iv[i++]= SpeciesNum;
    iv[i++]= SpeciesNum_e[0];
    iv[i++]= SpeciesNum_e[1];
    fwrite(iv,sizeof(int),i,fp);

    /* TCpyCell */

    i=0;
    iv[i++]= TCpyCell;
    iv[i++]= TCpyCell_e[0];
    iv[i++]= TCpyCell_e[1];
    fwrite(iv,sizeof(int),i,fp);

    /* TRAN_region */

    fwrite(TRAN_region,sizeof(int),atomnum+1,fp);

    /* TRAN_Original_Id */

    fwrite(TRAN_Original_Id,sizeof(int),atomnum+1,fp);

  }

  /**********************************************
       informations of the central region
  **********************************************/

  if (myid==Host_ID) {

    /* information of central region */

    fwrite(WhatSpecies,   sizeof(int), atomnum+1,  fp);
    fwrite(Spe_Total_CNO, sizeof(int), SpeciesNum, fp);

    fwrite(FNAN,sizeof(int),atomnum+1,fp);
    fwrite(&Max_FSNAN,sizeof(int),1,fp);
    fwrite(&ScaleSize,sizeof(double),1,fp);

    size1=(int)Max_FSNAN*ScaleSize+1;
    for (i=0; i<=atomnum; i++) {
      fwrite(natn[i],sizeof(int),size1,fp);
    }
    for (i=0; i<=atomnum; i++) {
      fwrite(ncn[i],sizeof(int),size1,fp);
    }

    size1=(TCpyCell+1)*4;
    ia_vec=(int*)malloc(sizeof(int)*size1);
    id=0;
    for (i=0; i<TCpyCell+1; i++) {
      for (j=0; j<=3; j++) {
	ia_vec[id++]=atv_ijk[i][j];
      }
    }
    fwrite(ia_vec,sizeof(int),size1,fp);
    free(ia_vec);
  }
   
  v1 = (double*)malloc(sizeof(double)*List_YOUSO[8]*List_YOUSO[7]*List_YOUSO[7]);

  /* OLP,  this is complex */

  for (k=0; k<4; k++) {

    int m,ID;
    /*global  Gc_AN  1:atomnum */
    /*variable ID = G2ID[Gc_AN] */
    /*variable Mc_AN = G2M[Gc_AN] */

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      Mc_AN = F_G2M[Gc_AN];
      ID = G2ID[Gc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];

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
	    tno1 = Spe_Total_CNO[Hwan];
	  } 

	  for (i=0;  i<tno0; i++){
	    for (j=0; j<tno1;j++) {
	      v1[vsize] =  OLP[k][Mc_AN][h_AN][i][j];
	      vsize++;
	    }
	  }
	} 

        /* Isend */

        if (myid!=Host_ID){
          tag = 999;
          MPI_Isend(&vsize, 1, MPI_INT, Host_ID, tag, comm1, &request);
          MPI_Wait(&request,&status);
          tag = 999;
          MPI_Isend(&v1[0], vsize, MPI_DOUBLE, Host_ID, tag, comm1, &request);
          MPI_Wait(&request,&status);
	}
        else{
          fwrite(v1, sizeof(double), vsize, fp);
        }
      }

      /* Recv */

      else if (ID!=myid && myid==Host_ID){
        tag = 999;
        MPI_Recv(&vsize, 1, MPI_INT, ID, tag, comm1, &status);
        tag = 999;
        MPI_Recv(&v1[0], vsize, MPI_DOUBLE, ID, tag, comm1, &status);
        fwrite(v1, sizeof(double), vsize, fp);
      }

      MPI_Barrier(comm1);

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
      tno0 = Spe_Total_CNO[Cwan];  

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
	    tno1 = Spe_Total_CNO[Hwan];
	  } 

          for (i=0; i<tno0; i++){
	    for (j=0;j<tno1; j++) {
	      v1[vsize] = H[k][Mc_AN][h_AN][i][j];
	      vsize++;
	    }
	  }
	}

        /* Isend */

        if (myid!=Host_ID){
          tag = 999;
          MPI_Isend(&vsize, 1, MPI_INT, Host_ID, tag, comm1, &request);
          MPI_Wait(&request,&status);
          tag = 999;
          MPI_Isend(&v1[0], vsize, MPI_DOUBLE, Host_ID, tag, comm1, &request);
          MPI_Wait(&request,&status);
	}
        else{
          fwrite(v1, sizeof(double), vsize, fp);
        }
      }

      /* Recv */

      else if (ID!=myid && myid==Host_ID){
        tag = 999;
        MPI_Recv(&vsize, 1, MPI_INT, ID, tag, comm1, &status);
        tag = 999;
        MPI_Recv(&v1[0], vsize, MPI_DOUBLE, ID, tag, comm1, &status);
        fwrite(v1, sizeof(double), vsize, fp);
      }

      MPI_Barrier(comm1);

    } /* Gc_AN */
  } /* k */

  free(v1);

  /**********************************************
              informations of leads
  **********************************************/

  if (myid==Host_ID){
    for (iside=0; iside<=1; iside++) {

      fwrite(WhatSpecies_e[iside],   sizeof(int), atomnum_e[iside]+1,  fp);
      fwrite(Spe_Total_CNO_e[iside], sizeof(int), SpeciesNum_e[iside], fp);
      fwrite(FNAN_e[iside],          sizeof(int), atomnum_e[iside]+1,  fp);
      fwrite(&Max_FSNAN_e[iside],    sizeof(int), 1,                   fp);
      fwrite(&ScaleSize_e[iside],    sizeof(double),1,                 fp);

      size1=(int)Max_FSNAN_e[iside]*ScaleSize_e[iside]+1;
      for (i=0; i<=atomnum_e[iside]; i++) {
	fwrite(natn_e[iside][i],sizeof(int),size1,fp);
      }
      for (i=0; i<=atomnum_e[iside]; i++) {
	fwrite(ncn_e[iside][i],sizeof(int),size1,fp);
      }

      size1=(TCpyCell_e[iside]+1)*4;
      ia_vec=(int*)malloc(sizeof(int)*size1);
      id=0;
      for (i=0; i<TCpyCell_e[iside]+1; i++) {
	for (j=0; j<=3; j++) {
	  ia_vec[id++]=atv_ijk_e[iside][i][j];
	}
      }
      fwrite(ia_vec,sizeof(int),size1,fp);
      free(ia_vec);

      /* overlap matrix */

      for (k=0; k<4 ; k++){
	for (GA_AN=1; GA_AN<=atomnum_e[iside]; GA_AN++){
	  wanA = WhatSpecies_e[iside][GA_AN];
	  tnoA = Spe_Total_CNO_e[iside][wanA];

	  for (LB_AN=0; LB_AN<=FNAN_e[iside][GA_AN]; LB_AN++){

	    GB_AN = natn_e[iside][GA_AN][LB_AN];
	    wanB = WhatSpecies_e[iside][GB_AN];
	    tnoB = Spe_Total_CNO_e[iside][wanB];

	    for (i=0; i<tnoA; i++){
	      fwrite(OLP_e[iside][k][GA_AN][LB_AN][i], sizeof(double), tnoB, fp);
	    }
	  }
	}   
      }

      /* Hamiltonian matrix */

      for (k=0; k<=SpinP_switch; k++){
	for (GA_AN=1; GA_AN<=atomnum_e[iside]; GA_AN++){
	  wanA = WhatSpecies_e[iside][GA_AN];
	  tnoA = Spe_Total_CNO_e[iside][wanA];

	  for (LB_AN=0; LB_AN<=FNAN_e[iside][GA_AN]; LB_AN++){

	    GB_AN = natn_e[iside][GA_AN][LB_AN];
	    wanB = WhatSpecies_e[iside][GB_AN];
	    tnoB = Spe_Total_CNO_e[iside][wanB];

	    for (i=0; i<tnoA; i++){
	      fwrite(H_e[iside][k][GA_AN][LB_AN][i], sizeof(double), tnoB, fp);
	    }
	  }
	}
      }
    }
  }

  /**********************************************
              close the file pointer 
  **********************************************/

  if (myid==Host_ID){
    fclose(fp);
  }

}



