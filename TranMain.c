/**********************************************************************
  TranMain.c:

  TranMain.c is a code to calculate electric transmission and current
  through the central region from the left lead to the right lead,
  in which the Hamiltonian and overlap matrices calculated by a
  non-equilibrium Green's function method are used  

  Log of TranMain.c:

     11/Dec/2005  released by H.Kino
     26/Feb/2006  modified by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "Inputtools.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#include "tran_prototypes.h"

#define eV2Hartree    27.2113845                
#define kB            0.00008617251324000000   /* eV/K  */          
#define PI            3.1415926535897932384626

#define Host_ID  0
#define PrintLevel  0

#define SCC_ref(i,j) ( ((j)-1)*NUM_c + (i)-1 )
#define SCL_ref(i,j) ( ((j)-1)*NUM_c + (i)-1 )
#define SCR_ref(i,j) ( ((j)-1)*NUM_c + (i)-1 )
#define S00l_ref(i,j) ( ((j)-1)*NUM_e[0]+(i)-1 )
#define S00r_ref(i,j) ( ((j)-1)*NUM_e[1]+(i)-1 )

int SpinP_switch,SpinP_switch2;
int NUM_c, NUM_e[2];
double E_Temp;

/* the center region */

double ChemP;
double *****OLP;
double *****H;
int atomnum;
int SpeciesNum;
int *WhatSpecies;
int *Spe_Total_CNO;
int *FNAN;
int **natn;
int **ncn;
int **atv_ijk; 
int Max_FSNAN;
double ScaleSize;
int TCpyCell;
int *TRAN_region;
int *TRAN_Original_Id;

/* the leads region */

double ChemP_e[2];
double *****OLP_e[2];
double *****H_e[2];
int atomnum_e[2];
int SpeciesNum_e[2];
int *WhatSpecies_e[2];
int *Spe_Total_CNO_e[2];
int *FNAN_e[2];
int **natn_e[2];
int **ncn_e[2];
int **atv_ijk_e[2]; 
int Max_FSNAN_e[2];
double ScaleSize_e[2];
int TCpyCell_e[2];

/* k-dependent matrices */

dcomplex *S00_e[2];
dcomplex *S01_e[2];
dcomplex **H00_e[2];
dcomplex **H01_e[2];

dcomplex *SCC;
dcomplex *SCL;
dcomplex *SCR;
dcomplex **HCC;
dcomplex **HCL;
dcomplex **HCR;

dcomplex ****tran_transmission;
dcomplex **tran_transmission_iv;

int tran_surfgreen_iteration_max;
double tran_surfgreen_eps;

int tran_bias_apply;
int SCF_tran_bias_apply;
double tran_biasvoltage_e[2];
int tran_transmission_on;
double tran_transmission_energyrange[3];
int tran_transmission_energydiv;
int tran_transmission_iv_on;
double tran_transmission_iv_energyrange[3];
int tran_transmission_iv_energydiv;
int TRAN_TKspace_grid2,TRAN_TKspace_grid3;
double ***current;
  
char filepath[100];
char filename[100];







void MTRAN_Read_Tran_HS( char *filepath, char *filename, char *ext ) 
{
  FILE *fp;
  int iv[20];
  double v[20];
  int i,j,k,id;
  int Gc_AN,Cwan,tno0,tno1;
  int h_AN,Gh_AN,Hwan;
  int iside,spin;
  int size1;
  int *ia_vec;
  char fname[300];

  sprintf(fname,"%s%s.%s",filepath,filename,ext);

  if (  (fp=fopen(fname,"r"))==NULL ) {
    printf("can not open %s\n",fname);
    printf("in MTRAN_Read_Tran_HS\n");
    exit(0);
  }

  /* SpinP_switch, NUM_c, and NUM_e */

  i=0;
  fread(iv, sizeof(int),4,fp);
  SpinP_switch = iv[i++];
  NUM_c    = iv[i++];
  NUM_e[0] = iv[i++];
  NUM_e[1] = iv[i++];

  if (PrintLevel){
    printf("spin=%d NUM_c=%d NUM_e=%d %d\n", SpinP_switch, NUM_c, NUM_e[0], NUM_e[1]);
  }

  /* chemical potential */

  i=0;
  fread(v,sizeof(double),3,fp);
  ChemP     = v[i++]; 
  ChemP_e[0]= v[i++];
  ChemP_e[1]= v[i++];

  if (PrintLevel){
    printf("ChemP= %10.5f  %10.5f  %10.5f\n", ChemP, ChemP_e[0], ChemP_e[1]);
  }

  /* tran_bias_apply */

  fread(iv, sizeof(int),1,fp);
  SCF_tran_bias_apply = iv[0];

  /* the number of atoms */

  i=0;
  fread(iv, sizeof(int),3,fp);
  atomnum      = iv[i++];
  atomnum_e[0] = iv[i++];
  atomnum_e[1] = iv[i++];

  if (PrintLevel){
    printf("atomnum=%d atomnum_e=%d %d\n", atomnum, atomnum_e[0], atomnum_e[1]);
  }

  /* the number of species */

  i=0;
  fread(iv, sizeof(int),3,fp);
  SpeciesNum      = iv[i++];
  SpeciesNum_e[0] = iv[i++];
  SpeciesNum_e[1] = iv[i++];

  if (PrintLevel){
    printf("SpeciesNum=%d SpeciesNum_e=%d %d\n",SpeciesNum, SpeciesNum_e[0], SpeciesNum_e[1]);
  }

  /* TCpyCell */

  i=0;
  fread(iv, sizeof(int),3,fp);
  TCpyCell      = iv[i++];
  TCpyCell_e[0] = iv[i++];
  TCpyCell_e[1] = iv[i++];

  if (PrintLevel){
    printf("TCpyCell=%d TCpyCell_e=%d %d\n",TCpyCell, TCpyCell_e[0], TCpyCell_e[1]);
  }

  /* TRAN_region */

  TRAN_region = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(TRAN_region,sizeof(int),atomnum+1,fp);

  /* TRAN_Original_Id */

  TRAN_Original_Id = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(TRAN_Original_Id,sizeof(int),atomnum+1,fp);

  /**********************************************
       informations of the central region
  **********************************************/

  WhatSpecies = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(WhatSpecies,   sizeof(int), atomnum+1,  fp);

  if (PrintLevel){
    for (i=0; i<=atomnum; i++) {
      printf("i=%2d WhatSpecies=%2d\n",i,WhatSpecies[i]);
    }
  }

  Spe_Total_CNO = (int*)malloc(sizeof(int)*(SpeciesNum));
  fread(Spe_Total_CNO, sizeof(int), SpeciesNum, fp);

  FNAN = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(FNAN,sizeof(int), atomnum+1,fp);

  if (PrintLevel){
    for (i=0; i<=atomnum; i++) {
      printf("i=%2d FNAN=%2d\n",i,FNAN[i]);
    }   
  }
  
  fread(&Max_FSNAN,sizeof(int),1,fp);
  fread(&ScaleSize,sizeof(double),1,fp);

  if (PrintLevel){
    printf("Max_FSNAN=%2d ScaleSize=%10.5f\n",Max_FSNAN,ScaleSize);
  }

  size1=(int)(Max_FSNAN*ScaleSize) + 1;

  natn = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (i=0; i<=(atomnum); i++) {
    natn[i] = (int*)malloc(sizeof(int)*size1);
  }
  for (i=0; i<=(atomnum); i++) {
    fread(natn[i],sizeof(int),size1,fp);
  }

  if (PrintLevel){
    for (i=0; i<=(atomnum); i++) {
      for (j=0; j<size1; j++) {
        printf("i=%3d j=%3d  natn=%2d\n",i,j,natn[i][j]);   
      }
    }
  }

  ncn = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (i=0; i<=(atomnum); i++) {
    ncn[i] = (int*)malloc(sizeof(int)*size1);
  }
  for (i=0; i<=(atomnum); i++) {
    fread(ncn[i],sizeof(int),size1,fp);
  }

  if (PrintLevel){
    for (i=0; i<=(atomnum); i++) {
      for (j=0; j<size1; j++) {
        printf("i=%3d j=%3d  ncn=%2d\n",i,j,natn[i][j]);   
      }
    }
  }

  size1=(TCpyCell+1)*4;
  ia_vec=(int*)malloc(sizeof(int)*size1);
  fread(ia_vec,sizeof(int),size1,fp);

  atv_ijk = (int**)malloc(sizeof(int*)*(TCpyCell+1));
  for (i=0; i<(TCpyCell+1); i++) {
    atv_ijk[i] = (int*)malloc(sizeof(int)*4);
  }

  id=0;
  for (i=0; i<(TCpyCell+1); i++) {
    for (j=0; j<=3; j++) {
      atv_ijk[i][j] = ia_vec[id++];
    }

    if (PrintLevel){
      printf("atv_ijk %3d   %2d %2d %2d\n",i,atv_ijk[i][1],atv_ijk[i][2],atv_ijk[i][3]);
    }

  }
  free(ia_vec);

  /* OLP */

  OLP = (double*****)malloc(sizeof(double****)*4);
  for (k=0; k<4; k++) {

    OLP[k] = (double****)malloc(sizeof(double***)*(atomnum+1));

    FNAN[0] = 0;
    for (Gc_AN=0; Gc_AN<=(atomnum); Gc_AN++){

      if (Gc_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_CNO[Cwan];
      }

      OLP[k][Gc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Gc_AN==0){
          tno1 = 1;
	}
        else {
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_CNO[Hwan];
	}

        OLP[k][Gc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0);

        for (i=0; i<tno0; i++){
          OLP[k][Gc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1);

          if (Gc_AN!=0){
	    fread(OLP[k][Gc_AN][h_AN][i],sizeof(double), tno1, fp); 

	    /*
            for (j=0; j<tno1; j++){
              printf("k=%2d Gc_AN=%2d h_AN=%2d i=%2d j=%2d OLP=%15.12f\n",
                      k,Gc_AN,h_AN,i,j,OLP[k][Gc_AN][h_AN][i][j]);
	    }        
	    */

	  }
        }
      }
    }
  }

  /* H */

  H = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
  for (k=0; k<(SpinP_switch+1); k++) {

    H[k] = (double****)malloc(sizeof(double***)*(atomnum+1));

    FNAN[0] = 0;
    for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){

      if (Gc_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_CNO[Cwan];
      }

      H[k][Gc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Gc_AN==0){
          tno1 = 1;
	}
        else {
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_CNO[Hwan];
	}

        H[k][Gc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0);

        for (i=0; i<tno0; i++){
          H[k][Gc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1);
          if (Gc_AN!=0){
            fread(H[k][Gc_AN][h_AN][i],sizeof(double),tno1,fp); 

	    /*
            for (j=0; j<tno1; j++){
              printf("k=%2d Gc_AN=%2d h_AN=%2d i=%2d j=%2d H=%15.12f\n",
                      k,Gc_AN,h_AN,i,j,H[k][Gc_AN][h_AN][i][j]);
	    } 
	    */       

	  }
        }
      }
    }
  }

  /**********************************************
              informations of leads
  **********************************************/

  for (iside=0; iside<=1; iside++) {

    WhatSpecies_e[iside] = (int*)malloc(sizeof(int)*(atomnum_e[iside]+1));
    fread(WhatSpecies_e[iside],   sizeof(int), atomnum_e[iside]+1,  fp);

    if (PrintLevel){
      for (i=0; i<=atomnum_e[iside]; i++) {
        printf("iside=%2d i=%2d WhatSpecies_e=%2d\n",iside,i,WhatSpecies_e[iside][i]);
      }
    }

    Spe_Total_CNO_e[iside] = (int*)malloc(sizeof(int)*SpeciesNum_e[iside]);
    fread(Spe_Total_CNO_e[iside], sizeof(int), SpeciesNum_e[iside], fp);

    FNAN_e[iside] = (int*)malloc(sizeof(int)*(atomnum_e[iside]+1));
    fread(FNAN_e[iside],          sizeof(int), atomnum_e[iside]+1,  fp);

    fread(&Max_FSNAN_e[iside],    sizeof(int), 1,                   fp);
    fread(&ScaleSize_e[iside],    sizeof(double),1,                 fp);

    size1=(int)(Max_FSNAN_e[iside]*ScaleSize_e[iside]) + 1;

    natn_e[iside] = (int**)malloc(sizeof(int*)*(atomnum_e[iside]+1));
    for (i=0; i<=atomnum_e[iside]; i++) {
      natn_e[iside][i] = (int*)malloc(sizeof(int)*size1);
    }
    for (i=0; i<=atomnum_e[iside]; i++) {
      fread(natn_e[iside][i],sizeof(int),size1,fp);
    }

    ncn_e[iside] = (int**)malloc(sizeof(int*)*(atomnum_e[iside]+1));
    for (i=0; i<=atomnum_e[iside]; i++) {
      ncn_e[iside][i] = (int*)malloc(sizeof(int)*size1);
    }
    for (i=0; i<=atomnum_e[iside]; i++) {
      fread(ncn_e[iside][i],sizeof(int),size1,fp);
    }

    size1=(TCpyCell_e[iside]+1)*4;
    ia_vec=(int*)malloc(sizeof(int)*size1);
    fread(ia_vec,sizeof(int),size1,fp);

    atv_ijk_e[iside] = (int**)malloc(sizeof(int*)*(TCpyCell_e[iside]+1));
    for (i=0; i<(TCpyCell_e[iside]+1); i++) {
      atv_ijk_e[iside][i] = (int*)malloc(sizeof(int)*4);
    }

    id=0;
    for (i=0; i<(TCpyCell_e[iside]+1); i++) {
      for (j=0; j<=3; j++) {
	atv_ijk_e[iside][i][j] = ia_vec[id++];
      }
      /*
      printf("atv_ijk_e %3d   %2d %2d %2d\n", 
              i,atv_ijk_e[iside][i][1],atv_ijk_e[iside][i][2],atv_ijk_e[iside][i][3]);
      */

    }
    free(ia_vec);

    /* overlap matrix */

    OLP_e[iside] = (double*****)malloc(sizeof(double****)*4);
    for (k=0; k<4; k++) {

      OLP_e[iside][k] = (double****)malloc(sizeof(double***)*(atomnum_e[iside]+1));

      FNAN_e[iside][0] = 0;
      for (Gc_AN=0; Gc_AN<=atomnum_e[iside]; Gc_AN++){

	if (Gc_AN==0){
          tno0 = 1;
	}
        else{
  	  Cwan = WhatSpecies_e[iside][Gc_AN];
  	  Cwan = WhatSpecies_e[iside][Gc_AN];
	  tno0 = Spe_Total_CNO_e[iside][Cwan];
	}

	OLP_e[iside][k][Gc_AN] = (double***)malloc(sizeof(double**)*(FNAN_e[iside][Gc_AN]+1));

	for (h_AN=0; h_AN<=FNAN_e[iside][Gc_AN]; h_AN++){
	  if (Gc_AN==0){
            tno1 = 1;
	  }
          else{
  	    Gh_AN = natn_e[iside][Gc_AN][h_AN];
	    Hwan = WhatSpecies_e[iside][Gh_AN];
	    tno1 = Spe_Total_CNO_e[iside][Hwan];
	  }

	  OLP_e[iside][k][Gc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0);

	  for (i=0; i<tno0; i++){
	    OLP_e[iside][k][Gc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1);
	    if (Gc_AN!=0){
	      fread(OLP_e[iside][k][Gc_AN][h_AN][i],sizeof(double),tno1,fp); 
	    }
	  }
	}
      }
    }

    /* Hamiltonian matrix */

    H_e[iside] = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (k=0; k<(SpinP_switch+1); k++) {

      H_e[iside][k] = (double****)malloc(sizeof(double***)*(atomnum_e[iside]+1));

      FNAN_e[iside][0] = 0;
      for (Gc_AN=0; Gc_AN<=atomnum_e[iside]; Gc_AN++){

	if (Gc_AN==0){
          tno0 = 1;
	}
        else{
  	  Cwan = WhatSpecies_e[iside][Gc_AN];
	  tno0 = Spe_Total_CNO_e[iside][Cwan];
	}

	H_e[iside][k][Gc_AN] = (double***)malloc(sizeof(double**)*(FNAN_e[iside][Gc_AN]+1));

	for (h_AN=0; h_AN<=FNAN_e[iside][Gc_AN]; h_AN++){

	  if (Gc_AN==0){
            tno1 = 1;
	  }
          else{ 
	    Gh_AN = natn_e[iside][Gc_AN][h_AN];
 	    Hwan = WhatSpecies_e[iside][Gh_AN];
	    tno1 = Spe_Total_CNO_e[iside][Hwan];
	  }

	  H_e[iside][k][Gc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0);

	  for (i=0; i<tno0; i++){
	    H_e[iside][k][Gc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1);
	    if (Gc_AN!=0){
	      fread(H_e[iside][k][Gc_AN][h_AN][i],sizeof(double),tno1,fp); 
	    }
	  }
	}
      }
    }

  }

  /**********************************************
              close the file pointer
  **********************************************/

  fclose(fp);
}














void MTRAN_Transmission(
                        MPI_Comm comm1,
			int SpinP_switch,
                        double ChemP_e[2],
			int NUM_c,
			int NUM_e[2],
			dcomplex **H00_e[2],
			dcomplex *S00_e[2],
			dcomplex **H01_e[2],
			dcomplex *S01_e[2],
			dcomplex **HCC,
			dcomplex **HCL,
			dcomplex **HCR,
			dcomplex *SCC,
			dcomplex *SCL,
			dcomplex *SCR, 
			double tran_surfgreen_iteration_max,
			double tran_surfgreen_eps, 
			double tran_transmission_energyrange[3],
			int tran_transmission_energydiv, 
			dcomplex **tran_transmission
			)
{
  dcomplex w;
  dcomplex *GRL,*GRR;
  dcomplex *GC_R,*GC_A;
  dcomplex *SigmaL_R,*SigmaL_A;
  dcomplex *SigmaR_R,*SigmaR_A;
  dcomplex *v1,*v2;
  dcomplex value;

  int iw,k,iside;
  int myid,numprocs,ID;
  int tag=99;

  int **iwIdx;
  int Miw,Miwmax ;
  int i,j;
  MPI_Status status;

  v1 = (dcomplex*) malloc(sizeof(dcomplex)*NUM_c*NUM_c);
  v2 = (dcomplex*) malloc(sizeof(dcomplex)*NUM_c*NUM_c);

  /* allocate */
  GRL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[0]* NUM_e[0]);
  GRR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[1]* NUM_e[1]);

  GC_R = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  GC_A = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaL_R = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaL_A = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaR_R = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaR_A = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);

  /* initialize */

  for (k=0; k<=1; k++) {
    for (iw=0; iw<tran_transmission_energydiv ; iw++) {
      tran_transmission[k][iw].r = 0.0;
      tran_transmission[k][iw].i = 0.0;
    }
  }

  /*parallel setup*/
  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  iwIdx=(int**)malloc(sizeof(int*)*numprocs);
  Miwmax = (tran_transmission_energydiv)/numprocs+1;
  for (i=0;i<numprocs;i++) {
    iwIdx[i]=(int*)malloc(sizeof(int)*Miwmax);
  }
  TRAN_Distribute_Node_Idx(0, tran_transmission_energydiv-1, numprocs, Miwmax,
                           iwIdx); /* output */

  /*parallel global iw 0:tran_transmission_energydiv-1 */
  /*parallel local  Miw 0:Miwmax-1 */
  /*parallel variable iw=iwIdx[myid][Miw] */

  for (Miw=0;Miw<Miwmax ; Miw++) {

    iw = iwIdx[myid][Miw];

    if ( iw>=0 ) {

      w.r = tran_transmission_energyrange[0] + ChemP_e[0]
            +
  	    (tran_transmission_energyrange[1]-tran_transmission_energyrange[0])*
	    (double)iw/(tran_transmission_energydiv-1);

      w.i = tran_transmission_energyrange[2];

      /*
      printf("iw=%d of %d  w= % 9.6e % 9.6e \n" ,iw, tran_transmission_energydiv,  w.r,w.i);
      */

      for (k=0; k<=SpinP_switch; k++) {

        /*****************************************************************
         Note that retarded and advanced Green functions and self energies
         are not conjugate comlex in case of the k-dependent case. 
        **************************************************************/ 

        /* in case of retarded ones */ 

	iside=0;
	TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside],
				   tran_surfgreen_iteration_max, tran_surfgreen_eps, GRL);

	TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL_R);
        
	iside=1;
	TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside],
				   tran_surfgreen_iteration_max, tran_surfgreen_eps, GRR);

	TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR_R);

	TRAN_Calc_CentGreen(w, NUM_c, SigmaL_R, SigmaR_R, HCC[k], SCC, GC_R);

        /* in case of advanced ones */ 

        w.i = -w.i;

	iside=0;
	TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside],
				   tran_surfgreen_iteration_max, tran_surfgreen_eps, GRL);

	TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL_A);
        
	iside=1;
	TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside],
				   tran_surfgreen_iteration_max, tran_surfgreen_eps, GRR);

	TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR_A);

	TRAN_Calc_CentGreen(w, NUM_c, SigmaL_A, SigmaR_A, HCC[k], SCC, GC_A);

        w.i = -w.i;

        /* calculation of transmission  */ 

	TRAN_Calc_OneTransmission(NUM_c, SigmaL_R, SigmaL_A, SigmaR_R, SigmaR_A, GC_R, GC_A, v1, v2 ,&value);

	tran_transmission[k][iw].r = value.r; 
	tran_transmission[k][iw].i = value.i;

        if (PrintLevel){
          printf("k=%2d w.r=%6.3f w.i=%6.3f value.r=%15.12f value.i=%15.12f\n",k,w.r,w.i,value.r,value.i);
	}

        if (SpinP_switch==0){
  	  tran_transmission[1][iw].r = value.r; 
	  tran_transmission[1][iw].i = value.i;
        }

      } /* for k */
    } /* if ( iw>=0 ) */
  } /* iw */


  /*parallel communication */
  for (k=0; k<=1; k++) {
    for (ID=0;ID<numprocs;ID++) {

      for (Miw=0;Miw<Miwmax ; Miw++) {
	double v[2];
	iw = iwIdx[ID][Miw];

	MPI_Barrier(comm1);

	v[0]= tran_transmission[k][iw].r;
	v[1]= tran_transmission[k][iw].i;

	if (iw>0) {
          MPI_Bcast(v, 2, MPI_DOUBLE, ID, comm1);
	}
	tran_transmission[k][iw].r=v[0];
	tran_transmission[k][iw].i=v[1];
      }
    }
  }

  free(GC_R);
  free(GC_A);
  free(SigmaL_R);
  free(SigmaL_A);
  free(SigmaR_R);
  free(SigmaR_A);

  free(GRR);
  free(GRL);
  free(v2);
  free(v1);

  for (i=0;i<numprocs;i++) {
    free(iwIdx[i]);
  }
  free(iwIdx);
}









void MTRAN_Current(MPI_Comm comm1,
		   int SpinP_switch,
                   double ChemP_e[2],
                   double E_Temp,
                   double tran_transmission_energyrange[3],
                   int tran_transmission_energydiv,
 		   dcomplex **tran_transmission,
                   double *current
		   )
{
  double dw,fL,fR,Beta,xL,xR;
  dcomplex w;
  dcomplex value;
  int iw,k,iside;
  int i,j;
  int myid,numprocs,ID;

  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  dw = (tran_transmission_energyrange[1]-tran_transmission_energyrange[0])/
       (double)(tran_transmission_energydiv-1);

  Beta = 1.0/kB/E_Temp;

  for (k=0; k<=SpinP_switch; k++) {

    current[k] = 0.0;

    for (iw=0; iw<tran_transmission_energydiv ; iw++) {

      w.r = tran_transmission_energyrange[0] + ChemP_e[0] + dw*(double)iw;
      w.i = tran_transmission_energyrange[2];

      xL = (w.r - ChemP_e[0])*Beta;
      fL = 1.0/(1.0 + exp(xL));

      xR = (w.r - ChemP_e[1])*Beta;
      fR = 1.0/(1.0 + exp(xR));
      current[k] -= (fL-fR)*tran_transmission[k][iw].r*dw;
    }

    /* atomic unit */
    current[k] /= (2.0*PI);

    /**********************************************************

     Current:

       convert the unit from a.u to ampere

       The unit of current is given by eEh/bar{h} in
       atomic unit, where e is the elementary charge,
       Eh Hartree, and bar{h} h/{2pi} given by

        e = 1.60217653 * 10^{-19} C
        Eh = 4.35974417 * 10^{-18} J
        bar{h} = 1.05457168 * 10^{-34} J s

      Therefore, 

      1 a.u.
         = 1.60217653 * 10^{-19} C * 4.35974417 * 10^{-18} J
           / (1.05457168 * 10^{-34} J s )
         = 6.6236178 * 10^{-3} [Cs^{-1}=ampere] 

     Electric potential:

       convert the unit from a.u to volt

       The unit of electric potential is given by Eh/e in
       atomic unit.

      Therefore, 

      1 a.u.
         = 4.35974417 * 10^{-18} J/ (1.60217653 * 10^{-19} C) 
         = 27.21138 [JC^{-1}=volt] 
    *********************************************************/
    
    current[k] *= 0.0066236178;

    if (PrintLevel){
      printf("myid=%2d k=%3d current= % 9.6e\n",myid,k,current[k]);fflush(stdout); 
    }
  }

  if (SpinP_switch==0) current[1] = current[0];
}





















void MTRAN_Output_Transmission(
        MPI_Comm comm1,
        char *fname,
        double k2,
        double k3,
        int SpinP_switch,
        int tran_transmission_energydiv,
        double tran_transmission_energyrange[3],
        dcomplex **tran_transmission
        )
{
  int iw,k;
  dcomplex w;
  int myid;
  FILE *fp;

  MPI_Comm_rank(comm1,&myid);
  if (myid!=Host_ID) { return; }

  printf("  %s\n",fname);

  if ( ( fp =fopen(fname,"w") )== NULL ) {
    printf("\ncan not open file to write transmission\n");
    printf("write transmission to stdout\n");
    exit(0);
  }

  fprintf(fp,"#/************************************************************/\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#  Current:\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    convert the unit from a.u to ampere\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    The unit of current is given by eEh/bar{h} in\n");
  fprintf(fp,"#    atomic unit, where e is the elementary charge,\n");
  fprintf(fp,"#    Eh Hartree, and bar{h} h/{2pi} given by\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#    e = 1.60217653 * 10^{-19} C\n");
  fprintf(fp,"#    Eh = 4.35974417 * 10^{-18} J\n");
  fprintf(fp,"#    bar{h} = 1.05457168 * 10^{-34} J s\n");    
  fprintf(fp,"#\n");
  fprintf(fp,"#    Therefore,\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#    1 a.u.\n");
  fprintf(fp,"#    = 1.60217653 * 10^{-19} C * 4.35974417 * 10^{-18} J\n");
  fprintf(fp,"#    / (1.05457168 * 10^{-34} J s )\n");
  fprintf(fp,"#    = 6.6236178 * 10^{-3} [Cs^{-1}=ampere]\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#  Electric potential:\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    convert the unit from a.u to volt\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    The unit of electric potential is given by Eh/e in\n");
  fprintf(fp,"#    atomic unit.\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    Therefore,\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    1 a.u.\n");
  fprintf(fp,"#    = 4.35974417 * 10^{-18} J/ (1.60217653 * 10^{-19} C)\n"); 
  fprintf(fp,"#    = 27.21138 [JC^{-1}=volt]\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#***********************************************************/\n");
  fprintf(fp,"#\n");
  fprintf(fp,"# k2= %8.4f  k3= %8.4f\n",k2,k3);
  fprintf(fp,"# SpinP_switch= %d\n",SpinP_switch);
  fprintf(fp,"# Chemical potential (Hartree) Left, Right Leads=% 9.6e % 9.6e\n",
            ChemP_e[0],ChemP_e[1]);
  fprintf(fp,"# diff Chemical potential (Hartree)=% 9.6e\n",
            ChemP_e[0]-ChemP_e[1]);
  fprintf(fp,"# tran_transmission_energydiv= %d\n", tran_transmission_energydiv);
  fprintf(fp,"#\n");
  fprintf(fp,"# iw w.real(au) w.imag(au) w.real(eV) w.imag(eV) trans.real(up) trans.imag(up) trans.real(down) trans.imag(down)\n");
  fprintf(fp,"#\n");

  for (iw=0; iw<tran_transmission_energydiv; iw++){

    w.r = tran_transmission_energyrange[0] + ChemP_e[0]
      + 
      (tran_transmission_energyrange[1]-tran_transmission_energyrange[0])*
      (double)iw/(tran_transmission_energydiv-1);

    w.i = tran_transmission_energyrange[2];

    fprintf(fp,"%3d % 9.6e % 9.6e % 9.6e % 9.6e % 9.6e % 9.6e % 9.6e % 9.6e\n",
	    iw,
	    w.r,
	    w.i,
	    w.r*eV2Hartree, w.i*eV2Hartree, 
	    tran_transmission[0][iw].r,
	    tran_transmission[0][iw].i,
	    tran_transmission[1][iw].r,
	    tran_transmission[1][iw].i );

  } /* iw */

  fclose(fp);
}


void MTRAN_Output_Current(
        MPI_Comm comm1,
        char *fname,
        int TRAN_TKspace_grid2,
        int TRAN_TKspace_grid3,
        int SpinP_switch,
        double ***current
        )
{
  int iw,i2,i3;
  double k2,k3;
  double crt0,crt1;
  dcomplex w;
  int myid;
  FILE *fp;

  MPI_Comm_rank(comm1,&myid);
  if (myid!=Host_ID) { return; }

  printf("  %s\n",fname);

  if ( ( fp =fopen(fname,"w") )== NULL ) {
    printf("\ncan not open file to write current\n");
    printf("write current to stdout\n");
    exit(0);
  }

  /* total current */

  crt0 = 0.0;
  crt1 = 0.0;
  for (i2=0; i2<TRAN_TKspace_grid2; i2++){
    for (i3=0; i3<TRAN_TKspace_grid3; i3++){
      crt0 += current[i2][i3][0];
      crt1 += current[i2][i3][1];
    }
  }

  crt0 /= (double)(TRAN_TKspace_grid2*TRAN_TKspace_grid3);
  crt1 /= (double)(TRAN_TKspace_grid2*TRAN_TKspace_grid3);

  fprintf(fp,"#/************************************************************/\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#  Current:\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    convert the unit from a.u to ampere\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    The unit of current is given by eEh/bar{h} in\n");
  fprintf(fp,"#    atomic unit, where e is the elementary charge,\n");
  fprintf(fp,"#    Eh Hartree, and bar{h} h/{2pi} given by\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#    e = 1.60217653 * 10^{-19} C\n");
  fprintf(fp,"#    Eh = 4.35974417 * 10^{-18} J\n");
  fprintf(fp,"#    bar{h} = 1.05457168 * 10^{-34} J s\n");    
  fprintf(fp,"#\n");
  fprintf(fp,"#    Therefore,\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#    1 a.u.\n");
  fprintf(fp,"#    = 1.60217653 * 10^{-19} C * 4.35974417 * 10^{-18} J\n");
  fprintf(fp,"#    / (1.05457168 * 10^{-34} J s )\n");
  fprintf(fp,"#    = 6.6236178 * 10^{-3} [Cs^{-1}=ampere]\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#  Electric potential:\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    convert the unit from a.u to volt\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    The unit of electric potential is given by Eh/e in\n");
  fprintf(fp,"#    atomic unit.\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    Therefore,\n");
  fprintf(fp,"#\n");
  fprintf(fp,"#    1 a.u.\n");
  fprintf(fp,"#    = 4.35974417 * 10^{-18} J/ (1.60217653 * 10^{-19} C)\n"); 
  fprintf(fp,"#    = 27.21138 [JC^{-1}=volt]\n"); 
  fprintf(fp,"#\n");
  fprintf(fp,"#***********************************************************/\n");
  fprintf(fp,"#\n");
  fprintf(fp,"# SpinP_switch= %d\n",SpinP_switch);
  fprintf(fp,"# Chemical potential (Hartree) Left, Right Leads=% 9.6e % 9.6e\n",
            ChemP_e[0],ChemP_e[1]);
  fprintf(fp,"# diff Chemical potential (Hartree)=% 9.6e\n",
            ChemP_e[0]-ChemP_e[1]);
  fprintf(fp,"# Corresponding bias voltage (Volt)=% 9.6e\n",
            (ChemP_e[0]-ChemP_e[1])*27.21138 );
  fprintf(fp,"# average current for up and down spins (ampere)=% 9.6e % 9.6e\n",crt0,crt1);
  fprintf(fp,"#\n");
  fprintf(fp,"# i2 i3 k2 k3  current (ampere, up spin)  current (ampere, down spin)\n");
  fprintf(fp,"#\n");

  for (i2=0; i2<TRAN_TKspace_grid2; i2++){
    k2 = -0.5 + (2.0*(double)i2+1.0)/(2.0*(double)TRAN_TKspace_grid2);
    for (i3=0; i3<TRAN_TKspace_grid3; i3++){
      k3 = -0.5 + (2.0*(double)i3+1.0)/(2.0*(double)TRAN_TKspace_grid3);

      fprintf(fp,"%3d %3d %8.4f %8.4f % 9.6e % 9.6e\n", 
              i2,i3,k2,k3,current[i2][i3][0],current[i2][i3][1]);
    }
  }

  fprintf(fp,"\n\n");

  fclose(fp);
}








void MTRAN_Input(int argc,
                 char *fname,
                 double ChemP_e[2],
                 int *TRAN_TKspace_grid2,
                 int *TRAN_TKspace_grid3,
		 int *SpinP_switch,
		 double *E_Temp,
		 int *tran_surfgreen_iteration_max,
		 double *tran_surfgreen_eps,
		 int *tran_bias_apply,
		 double tran_biasvoltage_e[2],
		 int  *tran_transmission_on,
		 double tran_transmission_energyrange[3],
		 int *tran_transmission_energydiv,
		 int  *tran_transmission_iv_on, 
		 double tran_transmission_iv_energyrange[3],
		 int *tran_transmission_iv_energydiv,
		 dcomplex ****(*tran_transmission),
		 dcomplex **(*tran_transmission_iv) 
		 )
{
  int i,i2,i3;
  double r_vec[20];
  int i_vec[20];
  int i_vec2[20];
  char *s_vec[20];

  input_open(fname);

  s_vec[0]="Off"; s_vec[1]="On"; s_vec[2]="NC";
  i_vec[0]=0    ; i_vec[1]=1   ; i_vec[2]=3;
  input_string2int("scf.SpinPolarization", SpinP_switch, 3, s_vec,i_vec);
  input_double("scf.ElectronicTemperature", E_Temp, 300);
  /* chage the unit from K to a.u. */
  *E_Temp = *E_Temp/eV2Hartree;

  input_int(   "Tran.Surfgreen.iterationmax", tran_surfgreen_iteration_max, 100);
  input_double("Tran.Surfgreen.convergeeps", tran_surfgreen_eps, 1.0e-5);

  input_logical("Tran.bias.apply",tran_bias_apply,0);
  if ( *tran_bias_apply ) {
    i=0;
    r_vec[i++]=0.0;
    r_vec[i++]=0.0;
    input_doublev("Tran.bias.voltage",i,tran_biasvoltage_e,r_vec);
  }
  else {
    tran_biasvoltage_e[0]=0.0;
    tran_biasvoltage_e[1]=0.0;
  }

  /****  k-points parallel to the layer, which are used for the transmission calc. ****/

  i_vec2[0]=1;
  i_vec2[1]=1;
  input_intv("Tran.tran.Kgrid",2,i_vec,i_vec2);
  *TRAN_TKspace_grid2 = i_vec[0];
  *TRAN_TKspace_grid3 = i_vec[1];

  if (*TRAN_TKspace_grid2<=0){
    printf("Tran.tran.Kgrid should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 

  if (*TRAN_TKspace_grid3<=0){
    printf("Tran.tran.Kgrid should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 

  input_logical("Tran.transmission.on",tran_transmission_on,1);
  if (tran_transmission_on) {
    i=0;
    r_vec[i++] = -2.0;
    r_vec[i++] = 2.0;
    r_vec[i++] = 1.0e-6;
    input_doublev("Tran.transmission.energyrange",i, tran_transmission_energyrange, r_vec);
    input_int("Tran.transmission.energydiv",tran_transmission_energydiv,100);

    *tran_transmission = (dcomplex****)malloc(sizeof(dcomplex***)*(*TRAN_TKspace_grid2));
    for (i2=0; i2<(*TRAN_TKspace_grid2); i2++){
      (*tran_transmission)[i2] = (dcomplex***)malloc(sizeof(dcomplex**)*(*TRAN_TKspace_grid3));
      for (i3=0; i3<(*TRAN_TKspace_grid3); i3++){
        (*tran_transmission)[i2][i3] = (dcomplex**)malloc(sizeof(dcomplex*)*3);
        for (i=0; i<3; i++) {
          (*tran_transmission)[i2][i3][i] = (dcomplex*)malloc(sizeof(dcomplex)*(*tran_transmission_energydiv));
	}
      }
    }
  }
  else {
    tran_transmission_energydiv=0;
    *tran_transmission=NULL;
  }
  input_logical("Tran.transmission.IV.on",tran_transmission_iv_on,1);
  if (*tran_transmission_iv_on) {
       
    i=0;
    r_vec[i++] = (ChemP_e[0]<ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
    r_vec[i++] = (ChemP_e[0]>ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
    r_vec[i++] = 1.0e-6;
    input_doublev("Tran.transmission.IV.energyrange",i, tran_transmission_iv_energyrange, r_vec);
    input_int("Tran.transmission.IV.energydiv",tran_transmission_iv_energydiv,15);
       
    (*tran_transmission_iv) = (dcomplex**)malloc(sizeof(dcomplex*)*(*SpinP_switch+1));
    for (i=0;i<= *SpinP_switch; i++) {
      (*tran_transmission_iv)[i]=(dcomplex*)malloc(sizeof(dcomplex)*(*tran_transmission_iv_energydiv));
    }

  }
  else {
    tran_transmission_iv_energydiv=0;
    (*tran_transmission_iv)=NULL;
  }

  input_close();
}


void   MTRAN_Input_Sys(int argc,char *file,
            char *filepath,
            char *filename)
{


  input_open(file);
  input_string("System.CurrrentDirectory",filepath,"./");
  input_string("System.Name",filename,"default");

  input_close();

}


void MTRAN_Set_MP(
        int job, 
        int anum, int *WhatSpecies, int *Spe_Total_CNO, 
        int *NUM,  /* output */
        int *MP    /* output */
)
{
  int Anum, i, wanA, tnoA;

 /* setup MP */
  Anum = 1;
  for (i=1; i<=anum; i++){
    if (job) MP[i]=Anum; 
    wanA = WhatSpecies[i];
    tnoA = Spe_Total_CNO[wanA];
    Anum += tnoA;
  }

  *NUM=Anum-1;
} 






void MTRAN_Set_SurfOverlap(char *position, 
                           double k2,
                           double k3,
                           int SpinP_switch,
                           int atomnum_e[2],
                           double *****OLP_e[2],
                           double *****H_e[2],
                           int SpeciesNum_e[2], 
                           int *WhatSpecies_e[2], 
                           int *Spe_Total_CNO_e[2], 
                           int *FNAN_e[2],
                           int **natn_e[2], 
                           int **ncn_e[2], 
                           int **atv_ijk_e[2],
			   dcomplex *S00_e[2],
			   dcomplex *S01_e[2],
                           dcomplex **H00_e[2],
			   dcomplex **H01_e[2]
                           )
#define S00_ref(i,j) ( ((j)-1)*NUM+(i)-1 ) 
{
  int NUM,n2;
  int Anum, Bnum, wanA, tnoA;
  int i,j,k;
  int GA_AN;
  int GB_AN, LB_AN,wanB, tnoB,Rn;
  int l1,l2,l3; 
  int spin,MN;
  int SpinP_switch_e[2]; 
  int direction,iside;
  double si,co,kRn;
  double s,h[10];
  static double epscutoff=1.0e-8;
  double epscutoff2;
  int *MP;

  /*debug */
  char msg[100];
  /*end debug */

  SpinP_switch_e[0] = SpinP_switch;
  SpinP_switch_e[1] = SpinP_switch;

  /* position -> direction */

  if      ( strcasecmp(position,"left")==0) {
    direction=-1;
    iside=0;
  }
  else if ( strcasecmp(position,"right")==0) {
    direction= 1;
    iside=1;
  } 

  /* set MP */
  MTRAN_Set_MP(0, atomnum_e[iside], WhatSpecies_e[iside], Spe_Total_CNO_e[iside], &NUM, MP);
  MP = (int*)malloc(sizeof(int)*(NUM+1));
  MTRAN_Set_MP(1, atomnum_e[iside], WhatSpecies_e[iside], Spe_Total_CNO_e[iside], &NUM, MP);

  n2 = NUM;   

  for (i=0; i<n2*n2; i++){
    S00_e[iside][i].r = 0.0;
    S00_e[iside][i].i = 0.0;
    S01_e[iside][i].r = 0.0;
    S01_e[iside][i].i = 0.0;
  }

  for (k=0; k<=SpinP_switch_e[iside]; k++) {
    for (i=0; i<n2*n2; i++){
      H00_e[iside][k][i].r = 0.0;
      H00_e[iside][k][i].i = 0.0;
      H01_e[iside][k][i].r = 0.0;
      H01_e[iside][k][i].i = 0.0;
    }
  }

  for (GA_AN=1; GA_AN<=atomnum_e[iside]; GA_AN++){

    wanA = WhatSpecies_e[iside][GA_AN];
    tnoA = Spe_Total_CNO_e[iside][wanA];
    Anum = MP[GA_AN];

    for (LB_AN=0; LB_AN<=FNAN_e[iside][GA_AN]; LB_AN++){

      GB_AN = natn_e[iside][GA_AN][LB_AN];
      Rn = ncn_e[iside][GA_AN][LB_AN];
      wanB = WhatSpecies_e[iside][GB_AN];
      tnoB = Spe_Total_CNO_e[iside][wanB];
      Bnum = MP[GB_AN];

      l1 = atv_ijk_e[iside][Rn][1];
      l2 = atv_ijk_e[iside][Rn][2];
      l3 = atv_ijk_e[iside][Rn][3];

      kRn = k2*(double)l2 + k3*(double)l3;
      si = sin(2.0*PI*kRn);
      co = cos(2.0*PI*kRn);

      if (l1==0) {
	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){

	    S00_e[iside][S00_ref(Anum+i,Bnum+j)].r += co*OLP_e[iside][0][GA_AN][LB_AN][i][j];
	    S00_e[iside][S00_ref(Anum+i,Bnum+j)].i += si*OLP_e[iside][0][GA_AN][LB_AN][i][j];

	    for (k=0;k<=SpinP_switch_e[iside]; k++ ){
	      H00_e[iside][k][S00_ref(Anum+i,Bnum+j)].r += co*H_e[iside][k][GA_AN][LB_AN][i][j];
	      H00_e[iside][k][S00_ref(Anum+i,Bnum+j)].i += si*H_e[iside][k][GA_AN][LB_AN][i][j];
	    }
	  }
	}
      }

      if (l1==direction) {

	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){

	    S01_e[iside][S00_ref(Anum+i,Bnum+j)].r += co*OLP_e[iside][0][GA_AN][LB_AN][i][j];
	    S01_e[iside][S00_ref(Anum+i,Bnum+j)].i += si*OLP_e[iside][0][GA_AN][LB_AN][i][j];

	    for (k=0; k<=SpinP_switch_e[iside]; k++ ){
	      H01_e[iside][k][S00_ref(Anum+i,Bnum+j)].r += co*H_e[iside][k][GA_AN][LB_AN][i][j];
	      H01_e[iside][k][S00_ref(Anum+i,Bnum+j)].i += si*H_e[iside][k][GA_AN][LB_AN][i][j];
	    }
	  }
	}
      }
    }
  }  /* GA_AN */
  
  for (GA_AN=1; GA_AN<=atomnum_e[iside]; GA_AN++){
    wanA = WhatSpecies_e[iside][GA_AN];
    tnoA = Spe_Total_CNO_e[iside][wanA];
    Anum = MP[GA_AN];
    for (i=0; i<tnoA;i++) {
      for (j=1; j<=NUM; j++) {

	MN = S00_ref(Anum+i,j); 
          
	if ( (fabs(S00_e[iside][MN].r)+fabs(S00_e[iside][MN].i)) < epscutoff ) {
	  S00_e[iside][MN].r = 0.0;
	  S00_e[iside][MN].i = 0.0;
	}
          
	if ( (fabs(S01_e[iside][MN].r)+fabs(S01_e[iside][MN].i)) < epscutoff ) {
	  S01_e[iside][MN].r = 0.0;
	  S01_e[iside][MN].i = 0.0;
	}
          
	for ( spin=0; spin<= SpinP_switch_e[iside]; spin++) {
	  if ( (fabs(H00_e[iside][spin][MN].r)+fabs(H00_e[iside][spin][MN].i)) < epscutoff ) {
	    H00_e[iside][spin][MN].r = 0.0;
	    H00_e[iside][spin][MN].i = 0.0;
	  }
	  if ( (fabs(H01_e[iside][spin][MN].r)+fabs(H01_e[iside][spin][MN].i)) < epscutoff ) {
	    H01_e[iside][spin][MN].r = 0.0;
	    H01_e[iside][spin][MN].i = 0.0;
	  }
	}
          
      } /* j */
    } /* i */ 
  } /* GA_AN */

  /* free arrays */

  free(MP);
} 




void MTRAN_Set_CentOverlap( 
			   MPI_Comm comm1,
			   int job, 
			   int SpinP_switch, 
			   double k2,
			   double k3,
			   int NUM_c,
			   int NUM_e[2],
			   double *****H, 
			   double *****OLP,
			   int atomnum,
			   int atomnum_e[2],
			   int *WhatSpecies,
			   int *WhatSpecies_e[2],
			   int *Spe_Total_CNO,
			   int *Spe_Total_CNO_e[2],
			   int *FNAN,
			   int **natn,
			   int **ncn, 
			   int **atv_ijk,
			   int *TRAN_region,
			   int *TRAN_Original_Id 
			   )
{
  int *MP, *MP_e[2];
  int i;

  /* setup MP */
  MP = (int*)malloc(sizeof(int)*(NUM_c+1));
  MTRAN_Set_MP( 1,  atomnum, WhatSpecies, Spe_Total_CNO, &NUM_c, MP);

  MP_e[0] = (int*)malloc(sizeof(int)*(NUM_e[0]+1));
  MTRAN_Set_MP( 1,  atomnum_e[0], WhatSpecies_e[0], Spe_Total_CNO_e[0], &i, MP_e[0]);

  MP_e[1] = (int*)malloc(sizeof(int)*(NUM_e[1]+1));
  MTRAN_Set_MP( 1,  atomnum_e[1], WhatSpecies_e[1], Spe_Total_CNO_e[1], &i, MP_e[1]);

  if ((job&1)==1) {

    int GA_AN, wanA, tnoA, Anum;
    int LB_AN, GB_AN, wanB, tnoB, l1,l2,l3, Bnum;
    int i,j,k;
    int Rn;
    double kRn,si,co;

    for (i=0;i<NUM_c*NUM_c;i++) {
      SCC[i].r = 0.0;
      SCC[i].i = 0.0;
      for (k=0;k<=SpinP_switch;k++) {
	HCC[k][i].r = 0.0;
	HCC[k][i].i = 0.0;
      }
    }

    /* make Overlap ,  HCC, SCC */
    /*parallel global GA_AN 1:atomnum */

    for (GA_AN=1; GA_AN<=atomnum; GA_AN++){

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

	if (l1!=0 ) continue; /* l1 is the direction to the electrode */

	/*
	 * if ( TRAN_region[GA_AN]==12  || TRAN_region[GB_AN]==12 )  continue;
	 * if ( TRAN_region[GA_AN]==13  || TRAN_region[GB_AN]==13 )  continue;
	 */
	/*      if ( TRAN_region[GA_AN]<10 && TRAN_region[GB_AN]<10 ) { */

        kRn = k2*(double)l2 + k3*(double)l3;
        si = sin(2.0*PI*kRn);
        co = cos(2.0*PI*kRn);

	for (i=0; i<tnoA; i++){
	  for (j=0; j<tnoB; j++){

	    SCC[SCC_ref(Anum+i,Bnum+j)].r += co*OLP[0][GA_AN][LB_AN][i][j];
	    SCC[SCC_ref(Anum+i,Bnum+j)].i += si*OLP[0][GA_AN][LB_AN][i][j];

	    for (k=0;k<=SpinP_switch; k++) {
	      HCC[k][SCC_ref(Anum+i,Bnum+j)].r += co*H[k][GA_AN][LB_AN][i][j];
	      HCC[k][SCC_ref(Anum+i,Bnum+j)].i += si*H[k][GA_AN][LB_AN][i][j];
	    }
	  }
	}

      } /* LB_AN */
    }   /* GA_AN */
  }     /* job&1 */

  if ( (job&2) == 2 ) {

    {
      int GA_AN, wanA, tnoA, Anum;
      int GA_AN_e, Anum_e; 
      int GB_AN, wanB, tnoB, Bnum;
      int GB_AN_e, Bnum_e; 
      int i,j,k;
      int iside;

      /* overwrite CL1 region */

      iside=0;

      /*parallel global GA_AN 1:atomnum */

      for (GA_AN=1; GA_AN<=atomnum; GA_AN++){

	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

	GA_AN_e = TRAN_Original_Id[GA_AN];
	Anum_e = MP_e[iside][GA_AN_e];

	for (GB_AN=1; GB_AN<=atomnum; GB_AN++){

	  if ( TRAN_region[GA_AN]==12  && TRAN_region[GB_AN]==12 )  {

	    wanB = WhatSpecies[GB_AN];
	    tnoB = Spe_Total_CNO[wanB];
	    Bnum = MP[GB_AN];

	    GB_AN_e = TRAN_Original_Id[GB_AN];
	    Bnum_e = MP_e[iside][GB_AN_e];

	    for (i=0; i<tnoA; i++){
	      for (j=0; j<tnoB; j++){
		SCC[SCC_ref(Anum+i,Bnum+j)] = S00_e[iside][S00l_ref(Anum_e+i,Bnum_e+j)];
		for (k=0; k<=SpinP_switch; k++) {
		  HCC[k][SCC_ref(Anum+i,Bnum+j)] = H00_e[iside][k][S00l_ref(Anum_e+i,Bnum_e+j)];
		}
	      }
	    }
	  }
	}
      }
    } 

    {
      int GA_AN, wanA, tnoA, Anum;
      int GA_AN_e, Anum_e;
      int GB_AN, wanB, tnoB, Bnum;
      int GB_AN_e, Bnum_e;
      int i,j,k;
      int iside;

      /* overwrite CR1 region */

      iside=1;

      /*parallel global GA_AN  1:atomnum */

      for (GA_AN=1; GA_AN<=atomnum; GA_AN++){

	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];

	GA_AN_e = TRAN_Original_Id[GA_AN];
	Anum_e = MP_e[iside][GA_AN_e]; /* = Anum */

	for (GB_AN=1; GB_AN<=atomnum; GB_AN++){

	  if ( TRAN_region[GA_AN]==13  && TRAN_region[GB_AN]==13)  {

	    wanB = WhatSpecies[GB_AN];
	    tnoB = Spe_Total_CNO[wanB];
	    Bnum = MP[GB_AN];

	    GB_AN_e = TRAN_Original_Id[GB_AN];
	    Bnum_e = MP_e[iside][GB_AN_e]; /* = Bnum */

	    for (i=0; i<tnoA; i++){
	      for (j=0; j<tnoB; j++){
		SCC[SCC_ref(Anum+i,Bnum+j)] = S00_e[iside][S00r_ref(Anum_e+i,Bnum_e+j)];
		for (k=0;k<=SpinP_switch; k++) {
		  HCC[k][SCC_ref(Anum+i,Bnum+j)] = H00_e[iside][k][S00r_ref(Anum_e+i,Bnum_e+j)];
		}
	      }
	    }
	  }
	}
      }
    }

    {
      int iside;
      int GA_AN, wanA, tnoA, Anum, GA_AN_e, Anum_e;
      int GB_AN_e, wanB_e, tnoB_e, Bnum_e;
      int i,j,k;

      /* make Overlap ,  HCL, SCL from OLP_e, and H_e*/
      iside=0;

      /*parallel global GA_AN  1:atomnum */

      for (GA_AN=1; GA_AN<=atomnum; GA_AN++){

	if (TRAN_region[GA_AN]%10!=2) continue;

	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	Anum = MP[GA_AN];  /* GA_AN is in C */

	GA_AN_e =  TRAN_Original_Id[GA_AN];
	Anum_e = MP_e[iside][GA_AN_e]; 

	for (GB_AN_e=1; GB_AN_e<=atomnum_e[iside]; GB_AN_e++) {

	  wanB_e = WhatSpecies_e[iside][GB_AN_e];
	  tnoB_e = Spe_Total_CNO_e[iside][wanB_e];
          Bnum_e = MP_e[iside][GB_AN_e];

          for (i=0; i<tnoA; i++){
            for (j=0; j<tnoB_e; j++){

              SCL[SCL_ref(Anum+i,Bnum_e+j)] = S01_e[iside][ S00l_ref(Anum_e+i, Bnum_e+j)];

              for (k=0; k<=SpinP_switch; k++) {
                HCL[k][SCL_ref(Anum+i,Bnum_e+j)] = H01_e[iside][k][ S00l_ref(Anum_e+i, Bnum_e+j)];
              }
            }
          }
	}
      }
    }

    {
      int iside;
      int GA_AN, wanA, tnoA, Anum, GA_AN_e, Anum_e;
      int GB_AN_e, wanB_e, tnoB_e, Bnum_e;
      int i,j,k;
      /* make Overlap ,  HCR, SCR from OLP_e, and H_e*/
      iside=1;
      for (GA_AN=1; GA_AN<=atomnum; GA_AN++){

        if (TRAN_region[GA_AN]%10!=3) continue;

        wanA = WhatSpecies[GA_AN];
        tnoA = Spe_Total_CNO[wanA];
        Anum = MP[GA_AN];  /* GA_AN is in C */

        GA_AN_e =  TRAN_Original_Id[GA_AN];
        Anum_e = MP_e[iside][GA_AN_e];

        for (GB_AN_e=1; GB_AN_e<=atomnum_e[iside];GB_AN_e++) {
          wanB_e = WhatSpecies_e[iside][GB_AN_e];
          tnoB_e = Spe_Total_CNO_e[iside][wanB_e];
          Bnum_e = MP_e[iside][GB_AN_e];
          for (i=0; i<tnoA; i++){
            for (j=0; j<tnoB_e; j++){

              SCR[SCR_ref(Anum+i,Bnum_e+j)] = S01_e[iside][ S00r_ref(Anum_e+i, Bnum_e+j)];

              for (k=0;k<=SpinP_switch; k++) {
                HCR[k][SCR_ref(Anum+i,Bnum_e+j)] = H01_e[iside][k][ S00r_ref(Anum_e+i, Bnum_e+j)];
              }
            }
          }
        }
      }
    }

  } /* job&2 */


  /* post-process */
  free(MP);
  free(MP_e[1]);
  free(MP_e[0]);
}



void MTRAN_Allocate_HS(int NUM_c,
                       int NUM_e[2],
                       int SpinP_switch)
{
  int i,side,spin;

  for (side=0;side<=1;side++) {

    S00_e[side] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[side]*NUM_e[side] );
    S01_e[side] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[side]*NUM_e[side] );
    for (i=0; i<(NUM_e[side]*NUM_e[side]); i++) {
      S00_e[side][i].r = 0.0;
      S00_e[side][i].i = 0.0;
      S01_e[side][i].r = 0.0;
      S01_e[side][i].i = 0.0;
    }

    H00_e[side] = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1) );
    H01_e[side] = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1) );
    for (spin=0; spin<=SpinP_switch; spin++) {
      H00_e[side][spin] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[side]*NUM_e[side] );
      H01_e[side][spin] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[side]*NUM_e[side] );
      for (i=0; i<(NUM_e[side]*NUM_e[side]); i++) {
        H00_e[side][spin][i].r = 0.0;
        H00_e[side][spin][i].i = 0.0;
        H01_e[side][spin][i].r = 0.0;
        H01_e[side][spin][i].i = 0.0;
      }
    }

    SCC = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_c ); 
    SCL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_e[0] ); 
    SCR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_e[1] ); 

    for (i=0; i<NUM_c*NUM_c; i++) {
      SCC[i].r = 0.0;
      SCC[i].i = 0.0;
    }
    for (i=0; i<NUM_c*NUM_e[0]; i++) {
      SCL[i].r = 0.0;
      SCL[i].i = 0.0;
    }
    for (i=0; i<NUM_c*NUM_e[1]; i++) {
      SCR[i].r = 0.0;
      SCR[i].i = 0.0;
    }

    HCC = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1)); 
    for (spin=0; spin<=SpinP_switch; spin++) {
      HCC[spin] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_c); 
      for (i=0; i<NUM_c*NUM_c; i++) {
        HCC[spin][i].r = 0.0;
        HCC[spin][i].i = 0.0;
      }
    }

    HCL = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1)); 
    for (spin=0; spin<=SpinP_switch; spin++) {
      HCL[spin] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_e[0]); 
      for (i=0; i<NUM_c*NUM_e[0]; i++) {
        HCL[spin][i].r = 0.0;
        HCL[spin][i].i = 0.0;
      }
    }

    HCR = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1)); 
    for (spin=0; spin<=SpinP_switch; spin++) {
      HCR[spin] = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c*NUM_e[1]); 
      for (i=0; i<NUM_c*NUM_e[1]; i++) {
        HCR[spin][i].r = 0.0;
        HCR[spin][i].i = 0.0;
      }
    }
  } 
}










int main(int argc, char **argv)
{
  int i,j,i2,i3,ii2,ii3,k,iw;
  int Gc_AN,h_AN,Gh_AN;
  int iside,tno0,tno1,Cwan,Hwan;
  int myid, numprocs;
  int **tran_flag;
  double k2,k3;
  MPI_Comm comm1;
  char fnameout[100];

  if (argc<2) {
     printf("usage:  %s   inputname \n",argv[0]);
     exit(0);
  }

  MPI_Init(&argc, &argv);
  comm1 = MPI_COMM_WORLD;   
  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  /**********************************************
                   read system
  **********************************************/

  MTRAN_Input_Sys(argc,argv[1],
            filepath,
            filename);

  /**********************************************
                 read tranb file
  **********************************************/

  MTRAN_Read_Tran_HS( filepath, filename, "tranb" );

  /**********************************************
                 read input file
  **********************************************/

  MTRAN_Input(argc,
              argv[1],
              ChemP_e, 
              /* output */
              &TRAN_TKspace_grid2,
              &TRAN_TKspace_grid3,
	      &SpinP_switch2,
              &E_Temp, 
	      &tran_surfgreen_iteration_max,
	      &tran_surfgreen_eps,
	      &tran_bias_apply,
	      tran_biasvoltage_e,
	      &tran_transmission_on,
	      tran_transmission_energyrange,
	      &tran_transmission_energydiv,
	      &tran_transmission_iv_on,
	      tran_transmission_iv_energyrange,
	      &tran_transmission_iv_energydiv,
	      &tran_transmission,
	      &tran_transmission_iv) ;

  if (SpinP_switch!=SpinP_switch2) {
     printf("SpinP_switch conflicts\n");
     printf("SpinP_switch=%d  SpinP_switch2=%d\n", SpinP_switch,SpinP_switch2);
     exit(0); 
  }

  MTRAN_Allocate_HS(NUM_c,NUM_e,SpinP_switch);

  /**********************************************
         apply the bias voltage if required 
  **********************************************/
  
  printf("ChemP_e %15.12f %15.12f\n",ChemP_e[0],ChemP_e[1]);
  /*  printf("SCF_tran_bias_apply %2d\n",SCF_tran_bias_apply); */

  if (tran_bias_apply){

    if (SCF_tran_bias_apply){

      if (myid==Host_ID){
        printf("\n Not supported to apply an additional bias voltage\n");
        printf(" when a finite bias was applied in the SCF calc.\n\n");
      }

      MPI_Finalize();
      exit(0);
    }

    ChemP_e[0] += tran_biasvoltage_e[0];
    ChemP_e[1] += tran_biasvoltage_e[1];

    for (iside=0; iside<=1; iside++) {
      for (k=0; k<(SpinP_switch+1); k++) {
        for (Gc_AN=1; Gc_AN<=atomnum_e[iside]; Gc_AN++){

  	  Cwan = WhatSpecies_e[iside][Gc_AN];
	  tno0 = Spe_Total_CNO_e[iside][Cwan];

  	  for (h_AN=0; h_AN<=FNAN_e[iside][Gc_AN]; h_AN++){

	    Gh_AN = natn_e[iside][Gc_AN][h_AN];
 	    Hwan = WhatSpecies_e[iside][Gh_AN];
	    tno1 = Spe_Total_CNO_e[iside][Hwan];

  	    for (i=0; i<tno0; i++){
   	      for (j=0; j<tno1; j++){

 	        H_e[iside][k][Gc_AN][h_AN][i][j] +=
                       tran_biasvoltage_e[iside]*OLP_e[iside][0][Gc_AN][h_AN][i][j];

	      } /* j     */
	    }   /* i     */
	  }     /* h_AN  */
	}       /* Gc_AN */
      }         /* k     */
    }           /* iside */
  } /* if (tran_bias_apply) */


  /**********************************************
              calculate transmission
  **********************************************/

  if (tran_transmission_on) {

    /* allocation of arrays */

    tran_flag = (int**)malloc(sizeof(int*)*TRAN_TKspace_grid2); 
    for (i2=0; i2<TRAN_TKspace_grid2; i2++){
      tran_flag[i2] = (int*)malloc(sizeof(int)*TRAN_TKspace_grid3); 
      for (i3=0; i3<TRAN_TKspace_grid3; i3++){
        tran_flag[i2][i3] = 0;
      }
    }

    current = (double***)malloc(sizeof(double**)*TRAN_TKspace_grid2); 
    for (i2=0; i2<TRAN_TKspace_grid2; i2++){
      current[i2] = (double**)malloc(sizeof(double*)*TRAN_TKspace_grid3); 
      for (i3=0; i3<TRAN_TKspace_grid3; i3++){
        current[i2][i3] = (double*)malloc(sizeof(double)*3);
      }
    }

    /* loop for k2 and k3 */

    if (myid==Host_ID) printf("\n  calculating...\n\n"); 

    for (i2=0; i2<TRAN_TKspace_grid2; i2++){

      k2 = -0.5 + (2.0*(double)i2+1.0)/(2.0*(double)TRAN_TKspace_grid2) + Shift_K_Point ;

      for (i3=0; i3<TRAN_TKspace_grid3; i3++){

        k3 = -0.5 + (2.0*(double)i3+1.0)/(2.0*(double)TRAN_TKspace_grid3) - Shift_K_Point;

        if (tran_flag[i2][i3]==0){

          if (myid==Host_ID) printf("  i2=%2d i3=%2d  k2=%8.4f k3=%8.4f\n",i2,i3,k2,k3); 

	  /* set Hamiltonian and overlap matrices of left and right leads */

	  MTRAN_Set_SurfOverlap( "left", k2,k3,SpinP_switch,atomnum_e,
				 OLP_e,H_e,SpeciesNum_e,WhatSpecies_e,Spe_Total_CNO_e,
				 FNAN_e,natn_e,ncn_e,atv_ijk_e,S00_e,S01_e,H00_e,H01_e );
  
	  MTRAN_Set_SurfOverlap( "right", k2,k3,SpinP_switch,atomnum_e,
				 OLP_e,H_e,SpeciesNum_e,WhatSpecies_e,Spe_Total_CNO_e,
				 FNAN_e,natn_e,ncn_e,atv_ijk_e,S00_e,S01_e,H00_e,H01_e );

	  /* set CC, CL and CR */

	  MTRAN_Set_CentOverlap(
				comm1,
				3,
				SpinP_switch, 
				k2,
				k3,
				NUM_c,
				NUM_e,
				H,
				OLP,
				atomnum,
				atomnum_e,
				WhatSpecies,
				WhatSpecies_e,
				Spe_Total_CNO,
				Spe_Total_CNO_e,
				FNAN,
				natn,
				ncn,
				atv_ijk,
				TRAN_region,
				TRAN_Original_Id
				);

	  /* calculate transmission */

	  MTRAN_Transmission(comm1,
			     SpinP_switch, 
                             ChemP_e,
			     NUM_c,
			     NUM_e,
			     H00_e,
			     S00_e,
			     H01_e,
			     S01_e,
			     HCC,
			     HCL,
			     HCR,
			     SCC,
			     SCL,
			     SCR,
			     tran_surfgreen_iteration_max,
			     tran_surfgreen_eps,
			     tran_transmission_energyrange,
			     tran_transmission_energydiv,

			     /* output */
			     tran_transmission[i2][i3]
			     );

	  /* calculate current */

	  MTRAN_Current(comm1,
			SpinP_switch,
			ChemP_e,
			E_Temp,
			tran_transmission_energyrange,
			tran_transmission_energydiv, 
			tran_transmission[i2][i3],
			current[i2][i3] );

	  /* taking account of the inversion symmetry */

          tran_flag[i2][i3] = 1;

          ii2 = TRAN_TKspace_grid2-1-i2;
          ii3 = TRAN_TKspace_grid3-1-i3;
 
          tran_flag[ii2][ii3] = 1;

	  for (k=0; k<=1; k++) {
	    for (iw=0; iw<tran_transmission_energydiv; iw++) {
	      tran_transmission[ii2][ii3][k][iw].r = tran_transmission[i2][i3][k][iw].r;
	      tran_transmission[ii2][ii3][k][iw].i = tran_transmission[i2][i3][k][iw].i; 
	    }
	  }

          current[ii2][ii3][0] = current[i2][i3][0];
          current[ii2][ii3][1] = current[i2][i3][1];
	}
      }
    }
  }

  /**********************************************
                output transmission
  **********************************************/

  if (tran_transmission_on) {

    if (myid==Host_ID) printf("\nTransmission:  files\n\n");
 
    for (i2=0; i2<TRAN_TKspace_grid2; i2++){
       k2 = -0.5 + (2.0*(double)i2+1.0)/(2.0*(double)TRAN_TKspace_grid2) + Shift_K_Point;

      for (i3=0; i3<TRAN_TKspace_grid3; i3++){
        k3 = -0.5 + (2.0*(double)i3+1.0)/(2.0*(double)TRAN_TKspace_grid3) - Shift_K_Point;

        sprintf(fnameout,"%s%s.tran%i_%i",filepath,filename,i2,i3);

	MTRAN_Output_Transmission(
				  comm1,
				  fnameout,
                                  k2,
                                  k3,
				  SpinP_switch, 
				  tran_transmission_energydiv,
				  tran_transmission_energyrange,
				  tran_transmission[i2][i3]
				  );

      }
    } 
  }

  /**********************************************
                   output current
  **********************************************/

  if (tran_transmission_on) {

    if (myid==Host_ID) printf("\nCurrent:  file\n\n");

    sprintf(fnameout,"%s%s.current",filepath,filename);

    MTRAN_Output_Current(
			 comm1,
			 fnameout,
                         TRAN_TKspace_grid2,
                         TRAN_TKspace_grid3,
			 SpinP_switch, 
			 current
			 );

    if (myid==Host_ID) printf("\n");
  }

       
  MPI_Finalize();
  exit(0);
}
