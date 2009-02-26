/**********************************************************************
  SetParaDFT.c:

     SetPara_DFT.c is a subroutine to set several parameters which
     are required in DFT calculations.

  Log of SetPara_DFT.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include "openmx_common.h"
#include "Inputtools.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

static void Set_Atom_Weight();
static void Set_Atom_Symbol();
static void Read_PAO(int spe, char *file);
static void Read_VPS(int spe, char *file);
static void ReadPara_DFT();
static void Set_BasisPara();
static void Set_Comp2Real();
static void Check_InitDensity();
static double V_Hart_atom(int Gensi, double R);
static double Int_phi0_phi1( double *phi0, double *phi1,
                             double *MXV, double *MRV, int Grid_Num );
static void Inverse(int n, double **a, double **ia);
static void Generate_PAO_Mixed_Basis(int spe);
static void Generate_VPS_Mixed_Basis(int spe);

static void output_structures();

void SetPara_DFT()
{
  int i,j,k,l,wanA,spe;    
  int nf,fg;
  int numprocs,myid;
  double MaxCutoffR;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID){
    printf("\n");
    printf("*******************************************************\n");   fflush(stdout);
    printf("                     PAO and VPS                       \n" );  fflush(stdout);
    printf("*******************************************************\n\n"); fflush(stdout); 
  }

  /****************************************************
            Gauss-Legendre quadrature grid
  ****************************************************/

  /* coarse mesh */  

  fg = 1;
  nf = GL_Mesh;
  Gauss_Legendre(GL_Mesh,GL_Abscissae,GL_Weight,&nf,&fg);

  /* fine mesh */  

  fg = 1;
  nf = FineGL_Mesh;
  Gauss_Legendre(FineGL_Mesh,FineGL_Abscissae,FineGL_Weight,&nf,&fg);

  /****************************************************
           Read parameters for DFT calculations
  ****************************************************/

  ReadPara_DFT();

  Check_InitDensity();

  Set_Atom_Weight();
  Set_Atom_Symbol();
  Set_Comp2Real();
  Set_BasisPara();

  /* output structures */

  if (myid==Host_ID) output_structures();

  if (myid==Host_ID){
    printf("\n");
    printf("*******************************************************\n");   fflush(stdout);
    printf("     Fourier transform of PAO and projectors of VNL    \n" );  fflush(stdout);
    printf("*******************************************************\n\n"); fflush(stdout); 
  }

  FT_PAO();
  FT_NLP();
  if (ProExpn_VNA==1){
    FT_ProExpn_VNA();
    FT_VNA();
    FT_ProductPAO();
  }

  /***************************************************
   allocation of arrays:

     dcomplex HOMOs_Coef[YOUSO33][2][YOUSO31]
                        [YOUSO1][YOUSO7];
     dcomplex LUMOs_Coef[YOUSO33][2][YOUSO32]
                        [YOUSO1][YOUSO7];

     int Bulk_Num_HOMOs[YOUSO33];
     int Bulk_Num_LUMOs[YOUSO33];
     int Bulk_HOMO[YOUSO33][2];
  ***************************************************/

  HOMOs_Coef = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[33]);
  for (i=0; i<List_YOUSO[33]; i++){
    HOMOs_Coef[i] = (dcomplex****)malloc(sizeof(dcomplex***)*2);
    for (j=0; j<2; j++){
      HOMOs_Coef[i][j] = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[31]);
      for (k=0; k<List_YOUSO[31]; k++){
        HOMOs_Coef[i][j][k] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[1]);
        for (l=0; l<List_YOUSO[1]; l++){
          HOMOs_Coef[i][j][k][l] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]);
	}
      }
    }
  }

  LUMOs_Coef = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[33]);
  for (i=0; i<List_YOUSO[33]; i++){
    LUMOs_Coef[i] = (dcomplex****)malloc(sizeof(dcomplex***)*2);
    for (j=0; j<2; j++){
      LUMOs_Coef[i][j] = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[32]);
      for (k=0; k<List_YOUSO[32]; k++){
        LUMOs_Coef[i][j][k] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[1]);
        for (l=0; l<List_YOUSO[1]; l++){
          LUMOs_Coef[i][j][k][l] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]);
	}
      }
    }
  }

  Bulk_Num_HOMOs = (int*)malloc(sizeof(int)*List_YOUSO[33]);
  Bulk_Num_LUMOs = (int*)malloc(sizeof(int)*List_YOUSO[33]);

  Bulk_HOMO = (int**)malloc(sizeof(int*)*List_YOUSO[33]);
  for (i=0; i<List_YOUSO[33]; i++){
    Bulk_HOMO[i] = (int*)malloc(sizeof(int)*2);
  }

  /***************************************************
   allocation of arrays:

   double **S;
   double *EV_S;
   double *IEV_S;

   In case of a cluster calculation, a full overlap
   matrix and a vector for the inverse of eigenvalue
   are allocated as global arrays.
  ***************************************************/

  if (Solver==2){
    j = 0;
    for (i=1; i<=atomnum; i++){
      wanA  = WhatSpecies[i];
      j  = j + Spe_Total_CNO[wanA];
    }
    Size_Total_Matrix = j;

    S = (double**)malloc(sizeof(double*)*(Size_Total_Matrix+2));
    for (i=0; i<(Size_Total_Matrix+2); i++){
      S[i]  = (double*)malloc(sizeof(double)*(Size_Total_Matrix+2));
    }

    EV_S = (double*)malloc(sizeof(double)*(Size_Total_Matrix+2));
    IEV_S = (double*)malloc(sizeof(double)*(Size_Total_Matrix+2));

    alloc_first[9] = 0;
  }

  /***************************************************
   Find the maximum value of Spe_Atom_Cut1[spe]
   and compare it to length of the cell vectors.
  ***************************************************/

  MaxCutoffR = 0.0;
  for (spe=0; spe<SpeciesNum; spe++){
    if (MaxCutoffR < Spe_Atom_Cut1[spe]) MaxCutoffR = Spe_Atom_Cut1[spe];
  }

  CellNN_flag = 1;
  MaxCutoffR = MaxCutoffR*MaxCutoffR;

  if ( (tv[1][1]*tv[1][1] + tv[1][2]*tv[1][2] + tv[1][3]*tv[1][3])<MaxCutoffR ) CellNN_flag = 0;
  if ( (tv[2][1]*tv[2][1] + tv[2][2]*tv[2][2] + tv[2][3]*tv[2][3])<MaxCutoffR ) CellNN_flag = 0;
  if ( (tv[3][1]*tv[3][1] + tv[3][2]*tv[3][2] + tv[3][3]*tv[3][3])<MaxCutoffR ) CellNN_flag = 0;

  if (atomnum<2000) CellNN_flag = 0;

}






void ReadPara_DFT()
{

  FILE *fp;
  int spe,po;
  int tmp_remake_headfile;

  char DirPAO[YOUSO10];
  char DirVPS[YOUSO10];
  char ExtPAO[YOUSO10] = ".pao";
  char ExtVPS[YOUSO10] = ".vps";
  char FN_PAO[YOUSO10];
  char FN_VPS[YOUSO10];
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   set DirPAO and DirVPS 
  ****************************************************/

  sprintf(DirPAO,"%s/PAO/",DFT_DATA_PATH);
  sprintf(DirVPS,"%s/VPS/",DFT_DATA_PATH);

  /****************************************************
   Read the data of pseudo atomic orbitals and density
  ****************************************************/

  /*************************
   find 

   List_YOUSO[21]
   List_YOUSO[24]
   List_YOUSO[25]
  *************************/

  remake_headfile = 1; 

  for (spe=0; spe<real_SpeciesNum; spe++){

    fnjoint2(DirPAO,SpeBasisName[spe],ExtPAO,FN_PAO);    

    if ((fp = fopen(FN_PAO,"r")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      Read_PAO(spe,FN_PAO);
      fclose(fp);
    }
    else{
      printf("Could not find %s\n",FN_PAO);
      exit(1);
    }
  }

  remake_headfile = 0;

  Allocate_Arrays(6);

  /*************************
            read
  *************************/

  for (spe=0; spe<real_SpeciesNum; spe++){

    fnjoint2(DirPAO,SpeBasisName[spe],ExtPAO,FN_PAO);    

    if ((fp = fopen(FN_PAO,"r")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      Read_PAO(spe,FN_PAO);
      fclose(fp);

      if (myid==Host_ID){  
        printf("<SetPara_DFT>  PAOs of species %s were normally found.\n",SpeName[spe]);
      }

    }
    else{
      printf("Could not find %s\n",FN_PAO);
      exit(1);
    }
  }

  /********************************
    generate mixed basis functions
  ********************************/

  if (Mixed_Basis_flag==1){
    Generate_PAO_Mixed_Basis(SpeciesNum-1);
  }

  /****************************************************
      Read the data of pseudopotentials and pcc
  ****************************************************/

  /*************************
   find 

   List_YOUSO[19]
   List_YOUSO[20]
   List_YOUSO[22]
   List_YOUSO[30]
   List_YOUSO[37]
  *************************/

  remake_headfile = 1; 

  for (spe=0; spe<real_SpeciesNum; spe++){

    fnjoint2(DirVPS,SpeVPS[spe],ExtVPS,FN_VPS);    
    if ((fp = fopen(FN_VPS,"r")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      Read_VPS(spe,FN_VPS);
      fclose(fp);
    }
    else{
      printf("Could not find %s\n",FN_VPS);
      exit(1);
    }
  }
  remake_headfile = 0;

  Allocate_Arrays(7);

  /*************************
            read 
  *************************/

  for (spe=0; spe<real_SpeciesNum; spe++){

    fnjoint2(DirVPS,SpeVPS[spe],ExtVPS,FN_VPS);

    if ((fp = fopen(FN_VPS,"r")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      Read_VPS(spe,FN_VPS);
      fclose(fp);

      if (myid==Host_ID){  
        printf("<SetPara_DFT>  VPSs of species %s were normally found.\n",
                 SpeName[spe]); 
      }
      if (VPS_j_dependency[spe]==0 && myid==Host_ID) {
        printf("               %s.vps is l-dependent.\n",SpeVPS[spe]);
      }
      else if (SO_switch==1 && VPS_j_dependency[spe]==1 && myid==Host_ID) {
        printf("               %s.vps is j-dependent.\n",SpeVPS[spe]);
      }
      else if (SO_switch==0 && VPS_j_dependency[spe]==1) {

        if (myid==Host_ID){ 
          printf("               %s.vps is j-dependent.\n",SpeVPS[spe]);
          printf("               In case of scf.SpinOrbit.Coupling=off,\n");
          printf("               j-dependent pseudo potentials are averaged by j-degeneracy,\n");
          printf("               which corresponds to a scalar relativistic treatment.\n");
	}

        /* switch off VPS_j_dependency */
        VPS_j_dependency[spe] = 0;
      }

    }
    else{
      printf("Could not find %s\n",FN_VPS);
      exit(1);
    }
  }

  if (2<=level_stdout){
    printf("<ReadPara_DFT>  YOUSO19=%2d YOUSO20=%2d\n",
           List_YOUSO[19],List_YOUSO[20]);
  }

  /*************************
      check j-dependency 
  *************************/

  if (SO_switch==1){
    po = 0;
    for (spe=0; spe<real_SpeciesNum; spe++){
      if (VPS_j_dependency[spe]==1)  po = 1;         
    }
    if (po==0){
      printf("               All VPSs are l-dependent.\n");
      printf("               scf.SpinOrbit.Coupling is changed to OFF.\n");
      SO_switch = 0;
    }
  }

  /********************************
    generate mixed basis functions
  ********************************/

  if (Mixed_Basis_flag==1){
    Generate_VPS_Mixed_Basis(SpeciesNum-1);
  }
}

void Read_PAO(int spe, char *file)
{
  FILE *fp;
  char file1[YOUSO10],file2[YOUSO10],file3[YOUSO10];
  int i,j,L;
  double dum;

  /****************************************************
                      open the file
  ****************************************************/

  input_open(file);

  /****************************************************
                       read data
  ****************************************************/

  input_double("AtomSpecies",&dum,0.0);
  Spe_WhatAtom[spe] = (int)dum; 
  /* dv_EH0 temporaliry is used. */
  dv_EH0[spe] = dum;
  
  input_int("grid.num.output",&Spe_Num_Mesh_PAO[spe],0);

  if (2<=level_stdout){
    printf("<Read_PAO>  spe=%2d  AtomNum=%2d\n",spe,Spe_WhatAtom[spe]);
    printf("<Read_PAO>  spe=%2d  Spe_Num_Mesh_PAO=%2d\n",spe,Spe_Num_Mesh_PAO[spe]);
  }

  if (List_YOUSO[21]<=(Spe_Num_Mesh_PAO[spe]+2))
    List_YOUSO[21] = Spe_Num_Mesh_PAO[spe] + 2;

  /* density */

  if (remake_headfile==0){
    if (fp=input_find("<valence.charge.density")) {
      for (i=0; i<Spe_Num_Mesh_PAO[spe]; i++){
	for (j=0; j<=2; j++){
	  if (fscanf(fp,"%lf",&dum)==EOF){
	    printf("format error of valence.charge.density in %s\n",file);
	    exit(1);
	  }
	  else{
	    if      (j==0) Spe_PAO_XV[spe][i] = dum;
	    else if (j==1) Spe_PAO_RV[spe][i] = dum;
	    else if (j==2) Spe_Atomic_Den[spe][i] = dum;
	  }
	}
      }    
      if (!input_last("valence.charge.density>")) {
	/* format error */
	printf("Format error in valence.charge.density\n");
	exit(1);
      }
    }
  }

  /* atomic orbitals */

  input_double("radial.cutoff.pao",&Spe_Atom_Cut1[spe],(double)5.0);
  input_int("PAO.Lmax",&Spe_PAO_LMAX[spe],0);
  input_int("PAO.Mul",&Spe_PAO_Mul[spe],0);

  if (List_YOUSO[25]<=(Spe_PAO_LMAX[spe]+1)) List_YOUSO[25] = Spe_PAO_LMAX[spe] + 1;
  if (List_YOUSO[24]<=Spe_PAO_Mul[spe])      List_YOUSO[24] = Spe_PAO_Mul[spe];

  if ((Spe_PAO_LMAX[spe]<Spe_MaxL_Basis[spe]) && remake_headfile==0){
    printf("Not enough data for PAO (%s)\n",file);
    exit(1);
  }
  for (L=0; L<=Spe_MaxL_Basis[spe]; L++){
    if ((Spe_PAO_Mul[spe]<Spe_Num_Basis[spe][L]) && remake_headfile==0){
      printf("Not enough data for PAO (%s)\n",file);
      exit(1);
    }
  }
 
  if (remake_headfile==0){
    for (L=0; L<=Spe_PAO_LMAX[spe]; L++){
      sprintf(file1,"<pseudo.atomic.orbitals.L=%i",L);
      sprintf(file2,"pseudo.atomic.orbitals.L=%i>",L);
      if (fp=input_find(file1)) {
	for (i=0; i<Spe_Num_Mesh_PAO[spe]; i++){
	  for (j=0; j<=(Spe_PAO_Mul[spe]+1); j++){
	    if (fscanf(fp,"%lf",&dum)==EOF){
	      printf("File error in pseudo.atomic.orbitals.L=%i\n",L);
	    }
	    else{
	      if (j==0)
		Spe_PAO_XV[spe][i] = dum;
	      else if (j==1)
		Spe_PAO_RV[spe][i] = dum;
	      else
		Spe_PAO_RWF[spe][L][j-2][i] = dum;
	    }
	  }
	}

      }
      if (!input_last(file2)) {
	/* format error */
	printf("Format error in pseudo.atomic.orbitals.L=%i\n",L);
	exit(1);
      }
    } 
  }

  /****************************************************
                    close the file
  ****************************************************/

  input_close();


}



void Read_VPS(int spe, char *file)
{
  FILE *fp;
  int m,n,k,t,i,j,L,ii,jj,NVPS,LVPS,tmp0,maxL,po;
  int number_vps,local_part_vps,charge_pcc_calc;
  double te,ve,VPS_Rcut,dum,dum0,dum1,r,dumping;
  double VNA_width,rmin,rmax,Sr,Dr,sum,dx;
  double *Vcore,**TmpVNL,**TmpMat;
  double **phi,*phi2,*pe;
  int numprocs,myid;
  int i_vec[20];
  char *s_vec[20];
 
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
  allocation of arrays:

  double Vcore[List_YOUSO[22]];
  double TmpVNL[2][List_YOUSO[19]];
  double TmpMat[List_YOUSO[37]][List_YOUSO[37]];
  ****************************************************/

  if (remake_headfile==0){
    Vcore = (double*)malloc(sizeof(double)*List_YOUSO[22]); 

    TmpVNL = (double**)malloc(sizeof(double*)*2);
    for (i=0; i<2; i++){
      TmpVNL[i] = (double*)malloc(sizeof(double)*List_YOUSO[19]);
    }

    TmpMat = (double**)malloc(sizeof(double*)*List_YOUSO[37]);
    for (i=0; i<List_YOUSO[37]; i++){
      TmpMat[i] = (double*)malloc(sizeof(double)*List_YOUSO[37]);
    }
  }

  /****************************************************
                      open the file
  ****************************************************/

  input_open(file);

  /****************************************************
                       read data
  ****************************************************/

  input_double("AtomSpecies",&dum0,0.0);
  tmp0 = (int)dum0;

  if (dv_EH0[spe]!=dum0){
    printf("Not the same atom in PAO and VPS for spe=%2d\n",spe);
    MPI_Finalize();
    exit(1);
  }

  input_double("total.electron",&te,(double)dum0);
  input_double("valence.electron",&ve,(double)0.0);
  Spe_Core_Charge[spe] = ve + dum0 - te;
  input_int("grid.num.output",&Spe_Num_Mesh_VPS[spe],0);

  if (List_YOUSO[22]<=(Spe_Num_Mesh_VPS[spe]+2)){
    List_YOUSO[22] = Spe_Num_Mesh_VPS[spe] + 2;
  }

  /****************************************************
        re-normalization of atomic charge density 
  ****************************************************/

  if (remake_headfile==0){

    sum = 0.0;
    dx = Spe_PAO_XV[spe][1] - Spe_PAO_XV[spe][0];
    for (i=0; i<Spe_Num_Mesh_PAO[spe]; i++){
      sum += Spe_Atomic_Den[spe][i]*exp(3.0*Spe_PAO_XV[spe][i]);
    }
    sum *= 4.0*PI*dx; 

    for (i=0; i<Spe_Num_Mesh_PAO[spe]; i++){
      Spe_Atomic_Den[spe][i] *= (Spe_Core_Charge[spe]/sum);
    }
  }

  /**************************************************
      check j-dependency 
  **************************************************/

  input_logical("j.dependent.pseudo.potentials",&VPS_j_dependency[spe],0);

  /**************************************************
       a new format from adpack1.2 and abred1.3
  **************************************************/
 
  if (fp=input_find("<project.energies")){

    if (2<=level_stdout){
      printf("<Read_VPS>  VPS of %s was a format of ADPACK1.7\n",SpeVPS[spe]);
    }

    fscanf(fp,"%d",&Spe_Num_RVPS[spe]);
    t = 0;
    maxL = 0;

    for (i=0; i<Spe_Num_RVPS[spe]; i++){

      /* l-dependent */
      if (VPS_j_dependency[spe]==0){
	fscanf(fp,"%d %lf",&LVPS,&dum0);
	if (remake_headfile==0){
	  Spe_VNLE[0][spe][i] = dum0;
	  Spe_VPS_List[spe][i+1] = LVPS;
	}
      }

      /* j-dependent */
      else if (VPS_j_dependency[spe]==1){
	fscanf(fp,"%d %lf %lf",&LVPS,&dum0,&dum1);
	if (remake_headfile==0){

	  if (SO_switch==1){
	    Spe_VNLE[0][spe][i] = dum0;
	    Spe_VNLE[1][spe][i] = dum1;
	    Spe_VPS_List[spe][i+1] = LVPS;
	  }

	  /* j-averaging for a scalar relativistic treatment */
	  else if (SO_switch==0){
	    Spe_VNLE[0][spe][i] = ( (double)(LVPS+1)*dum0+(double)LVPS*dum1 )/(double)(2*LVPS+1);
	    Spe_VPS_List[spe][i+1] = LVPS;
	  }

	}
      }

      /* count the number of projectors */
      t = t + 2*LVPS + 1;
      if (maxL<LVPS) maxL = LVPS;
      if (2<=level_stdout && remake_headfile==0){
	printf("<Read_VPS>  i=%2d  Spe_VPS_List=%2d\n",i,Spe_VPS_List[spe][i+1]);
      }
    }
    Spe_Total_VPS_Pro[spe] = t; 

    if (List_YOUSO[30]<=maxL){
      List_YOUSO[30] = maxL + 2;
    }
    if (List_YOUSO[30]==0){
      List_YOUSO[30] = 2;
    }

    if (List_YOUSO[19]<=Spe_Num_RVPS[spe]){
      List_YOUSO[19] = Spe_Num_RVPS[spe] + 2; 
    }
    if (List_YOUSO[19]==0){
      List_YOUSO[19] = 3;
    }

    if (List_YOUSO[20]<=Spe_Total_VPS_Pro[spe]){
      List_YOUSO[20] = Spe_Total_VPS_Pro[spe] + 2;
    }
    if (List_YOUSO[20]==0){
      List_YOUSO[20] = 3;
    }

    if (2<=level_stdout){
      printf("<Read_VPS>  spe=%2d  Spe_Total_VPS_Pro=%2d\n",spe,Spe_Total_VPS_Pro[spe]);
      printf("<Read_VPS>  spe=%2d  Spe_Num_RVPS=%2d\n",spe,Spe_Num_RVPS[spe]);
    }

    if (!input_last("project.energies>")){
      /* format error */
      printf("Format error for project.energies in %s\n",file);
      MPI_Finalize();
      exit(1);
    }

    /* Pseudopotentials */

    if (remake_headfile==0){
      if (fp=input_find("<Pseudo.Potentials")){
        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){

	  for (j=0; j<=( (VPS_j_dependency[spe]+1)*Spe_Num_RVPS[spe]+2 ); j++){

	    if (fscanf(fp,"%lf",&dum)==EOF){
	      printf("format error of Pseudo_Potentials in %s\n",file);
              MPI_Finalize();
	      exit(1);
	    }
            else{
	      if (j==0)
                Spe_VPS_XV[spe][i] = dum;
	      else if (j==1)
                Spe_VPS_RV[spe][i] = dum;
	      else if (j==2)
                Vcore[i] = dum;
              else{

                /* l-dependent */
                if      (VPS_j_dependency[spe]==0){
                  Spe_VNL[0][spe][j-3][i] = dum;
		}
                /* j-dependent */
                else if (VPS_j_dependency[spe]==1){

                  if (SO_switch==1){
                    if (j%2==1)  Spe_VNL[0][spe][ (j-3)/2 ][i] = dum;
                    else         Spe_VNL[1][spe][ (j-3)/2 ][i] = dum;
		  }

                  /* j-averaging for a scalar relativistic treatment */
                  else if (SO_switch==0){
                    if (j%2==1)  TmpVNL[0][ (j-3)/2 ] = dum;
                    else         TmpVNL[1][ (j-3)/2 ] = dum;
		  }
		}
	      }
	    }
	  } /* j */

          /* j-averaging for a scalar relativistic treatment */
          if ( VPS_j_dependency[spe]==1 && SO_switch==0 ){
  	    for (j=0; j<Spe_Num_RVPS[spe]; j++){

              LVPS = Spe_VPS_List[spe][j+1];
              Spe_VNL[0][spe][j][i] = ( (double)(LVPS+1)*TmpVNL[0][j]
                                       +(double)LVPS*TmpVNL[1][j] )/( (double)(2*LVPS+1) );
	    }

	  }

	} /* i */
        if (!input_last("Pseudo.Potentials>")) {
	  /* format error */
	  printf("Format error for Pseudo.Potentials\n");
          MPI_Finalize();
	  exit(1);
        }
      }
    }

  }

  /**************************************************
     an old format before adpack1.2 and abred1.3
  **************************************************/

  else{

    if (2<=level_stdout){
      printf("<Read_VPS>  VPS of %s was a format of ADPACK1.1\n",SpeVPS[spe]);
    }

    input_int("number.vps",&number_vps,0);
    input_int("local.part.vps",&local_part_vps,0);

    if (2<=level_stdout){
      printf("<Read_VPS>  spe=%2d  AtomNum=%2d\n",spe,Spe_WhatAtom[spe]);
      printf("<Read_VPS>  spe=%2d  Spe_Core_Charge=%5.2f\n",spe,Spe_Core_Charge[spe]);
      printf("<Read_VPS>  spe=%2d  Spe_Num_Mesh_VPS=%2d\n",spe,Spe_Num_Mesh_VPS[spe]);
      printf("<Read_VPS>  spe=%2d  number.vps=%2d\n",spe,number_vps);
      printf("<Read_VPS>  spe=%2d  local_part_vps=%2d\n",spe,local_part_vps);
    }

    if (fp=input_find("<pseudo.NandL")) {
      k = 0;
      t = 0;
      maxL = 0;
      for (i=0; i<number_vps; i++){
        fscanf(fp,"%d %d %d %lf",&j,&NVPS,&LVPS,&VPS_Rcut);
        if (local_part_vps!=i){
	  k++;
          t = t + 2*LVPS + 1;
          if (maxL<LVPS) maxL = LVPS;

          if (remake_headfile==0){
	    Spe_VPS_List[spe][k] = LVPS;
	  }

          if (2<=level_stdout && remake_headfile==0){
            printf("<Read_VPS>  k=%2d  Spe_VPS_List=%2d\n",k,Spe_VPS_List[spe][k]);
	  }
        }
      }
      Spe_Total_VPS_Pro[spe] = t; 
      Spe_Num_RVPS[spe] = k;

      if (List_YOUSO[30]<=maxL){
        List_YOUSO[30] = maxL + 2;
      }
      if (List_YOUSO[30]==0){
        List_YOUSO[30] = 2;
      }

      if (List_YOUSO[19]<=Spe_Num_RVPS[spe]){
        List_YOUSO[19] = Spe_Num_RVPS[spe] + 2;
      }
      if (List_YOUSO[19]==0){
        List_YOUSO[19] = 3;
      }

      if (List_YOUSO[20]<=Spe_Total_VPS_Pro[spe]){
        List_YOUSO[20] = Spe_Total_VPS_Pro[spe] + 2;
      }
      if (List_YOUSO[20]==0){
        List_YOUSO[20] = 3;
      }

      if (List_YOUSO[37]==0){
        List_YOUSO[37] = 2;
      }

      if (2<=level_stdout){
        printf("<Read_VPS>  spe=%2d  Spe_Total_VPS_Pro=%2d\n",spe,Spe_Total_VPS_Pro[spe]);
        printf("<Read_VPS>  spe=%2d  Spe_Num_RVPS=%2d\n",spe,Spe_Num_RVPS[spe]);
      }

      if (!input_last("pseudo.NandL>")) {
        /* format error */
        printf("Format error for pseudo.NandL in %s\n",file);
        exit(1);
      }

    }

    /* projection.energies */

    if (remake_headfile==0){
      if (0<Spe_Num_RVPS[spe]){
        if (fp=input_find("<projection.energies")) {
	  for (i=0; i<Spe_Num_RVPS[spe]; i++){

	    fscanf(fp,"%d %lf",&j,&Spe_VNLE[0][spe][i]);
            if (2<=level_stdout){
              printf("<Read_VPS>  i=%2d Spe_VNLE=%15.12f\n",j,Spe_VNLE[0][spe][i]); 
	    }
	  }    
          if (!input_last("projection.energies>")) {
	    /* format error */
	    printf("Format error for projection.energies in %s\n",file);
	    exit(1);
	  }
        }
      }
    }

    /* Pseudo_Potentials */

    if (remake_headfile==0){
      if (fp=input_find("<Pseudo.Potentials")){
        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
  	  for (j=0; j<=(Spe_Num_RVPS[spe]+2); j++){
	    if (fscanf(fp,"%lf",&dum)==EOF){
	      printf("format error of Pseudo_Potentials in %s\n",file);
	      exit(1);
	    }
            else{
	      if (j==0)
                Spe_VPS_XV[spe][i] = dum;
	      else if (j==1)
                Spe_VPS_RV[spe][i] = dum;
	      else if (j==2)
                Vcore[i] = dum;
              else
                Spe_VNL[0][spe][j-3][i] = dum;
	    }
	  }
        }
        if (!input_last("Pseudo.Potentials>")) {
	  /* format error */
	  printf("Format error for Pseudo.Potentials\n");
	  exit(1);
        }
      }
    }		
  }	

  /* PCC */

  if (remake_headfile==0){
    input_logical("charge.pcc.calc",&charge_pcc_calc,0);

    /* initialize Spe_Atomic_PCC */
    for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++) {
      Spe_Atomic_PCC[spe][i]=0.0;
    }

    if (charge_pcc_calc==1){
      if ((fp=input_find("<density.PCC")) != NULL) {
	for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
	  for (j=0; j<=2; j++){
	    if (fscanf(fp,"%lf",&dum)==EOF){
	      printf("format error of charge.pcc.calc in %s\n",file);
	      exit(1);
	    }
	    else{
	      if (j==2){
		Spe_Atomic_PCC[spe][i] = dum;
	      }
	    }
	  }
	}    
	if (!input_last("density.PCC>")) {
	  /* format error */
	  printf("Format error for density.PCC\n");
	  exit(1);
	}
      }
      else{
	printf("There is no data for PCC in %s\n",file);
	exit(1);
      }
    }
  }

  /****************************************************
                    close the file
  ****************************************************/
  
  input_close();
  
  /****************************************************
                calculate VH_Atom and Vna
  ****************************************************/

  if (remake_headfile==0){

    /* calculation of Spe_VH_Atom */

    for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
      r = Spe_VPS_RV[spe][i];
      Spe_VH_Atom[spe][i] = V_Hart_atom(spe,r);

      /*
     printf("A1 %2d %15.12f %15.12f %15.12f %15.12f\n",
              spe,r,Vcore[i],Spe_VH_Atom[spe][i],Spe_Vna[spe][i]);
      */

    }


    /* the nearest point to Spe_Atom_Cut1[spe] */

    dum = 1000.0;

    for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){

      r = Spe_VPS_RV[spe][i];
      dum1 = fabs(r-Spe_Atom_Cut1[spe]);

      if (dum1<dum){
        dum = dum1;
        ii = i;
      }
    }

    /* correct the asymptotic behaviour of Spe_VH_Atom */

    if (1.0e-15<Spe_Core_Charge[spe]){

      dum = -Vcore[ii]/Spe_VH_Atom[spe][ii];
      for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
        Spe_VH_Atom[spe][i] *= dum;
      }
    }

    /* calculation of Spe_Vna */


    for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
      r = Spe_VPS_RV[spe][i];
      dumping = 1.0/(1.0+exp( 20.0*(r-Spe_Atom_Cut1[spe])) );
      Spe_Vna[spe][i] = dumping*(Vcore[i] + Spe_VH_Atom[spe][i]);

      /*
      printf("%2d %15.12f %15.12f %15.12f %15.12f\n",
              spe,r,Vcore[i],Spe_VH_Atom[spe][i],Spe_Vna[spe][i]);
      */

    }
  }

  /****************************************************
              projector expansion of VNA
  ****************************************************/

  if (remake_headfile==0 && ProExpn_VNA==1){

    /* allocaltion of arrays */
    phi = (double**)malloc(sizeof(double*)*List_YOUSO[34]);
    for (k=0; k<List_YOUSO[34]; k++){
      phi[k] = (double*)malloc(sizeof(double)*Spe_Num_Mesh_VPS[spe]);
    }

    for (k=0; k<List_YOUSO[34]; k++){
      phi2 = (double*)malloc(sizeof(double)*Spe_Num_Mesh_VPS[spe]);
    }

    pe = (double*)malloc(sizeof(double)*List_YOUSO[34]);

    for (L=0; L<=List_YOUSO[35]; L++){ 

      /* set initial wave functions */
      for (m=0; m<List_YOUSO[34]; m++){

        if (L<=Spe_PAO_LMAX[spe] && m<Spe_PAO_Mul[spe]){
          for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
            r = Spe_VPS_RV[spe][i];
            phi[m][i] = RadialF(spe,L,m,r);
	  }
	}

        else if (L<=Spe_PAO_LMAX[spe]){
          for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
            r = Spe_VPS_RV[spe][i];
  	    phi[m][i] = pow(0.1*Spe_Vna[spe][i],(double)m)*phi[0][i];
	  }
        }

        else if (m<Spe_PAO_Mul[spe]){
          for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
            r = Spe_VPS_RV[spe][i];
  	    phi[m][i] = RadialF(spe,Spe_PAO_LMAX[spe],m,r)
                        *pow(r,(double)(L-Spe_PAO_LMAX[spe]));
	  }
        }

        else{
          for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
            r = Spe_VPS_RV[spe][i];
  	    phi[m][i] = RadialF(spe,Spe_PAO_LMAX[spe],Spe_PAO_Mul[spe]-1,r)*
                        pow(0.1*Spe_Vna[spe][i],(double)(m-Spe_PAO_Mul[spe]+1));
	  }
        }

      }

      /* Normalization */
      for (m=0; m<List_YOUSO[34]; m++){
        dum0 = Int_phi0_phi1(phi[m], phi[m], Spe_VPS_XV[spe],
                             Spe_VPS_RV[spe], Spe_Num_Mesh_VPS[spe]);

        if (1.0e-17<dum0)  dum0 = 1.0/sqrt(dum0);
        else               dum0 = 0.0;  

        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
          phi[m][i] = dum0*phi[m][i];
        }
      }

      /* Gramm-Schmidt orthogonalization with a norm defined by <f|v|g> */
      for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
        Projector_VNA[spe][L][0][i] = phi[0][i];
        r = Spe_VPS_RV[spe][i];
        phi2[i] = phi[0][i]*Spe_Vna[spe][i];
      }

      dum0 = Int_phi0_phi1(Projector_VNA[spe][L][0], phi2, Spe_VPS_XV[spe], 
                                  Spe_VPS_RV[spe], Spe_Num_Mesh_VPS[spe]); 

      /* empty atom */ 
      if (Spe_WhatAtom[spe]==0){
        pe[0] = 0.0;
      }

      else if (fabs(dum0)<1.0e-15){
        pe[0] = 0.0;
      }

      /* normal atom */ 
      else{   
        pe[0] = 1.0/dum0;
      }

      for (m=1; m<List_YOUSO[34]; m++){
 
        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++) Projector_VNA[spe][L][m][i] = 0.0;

        for (n=0; n<m; n++){

          for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++) phi2[i] = phi[m][i]*Spe_Vna[spe][i];

	  dum0 = Int_phi0_phi1(Projector_VNA[spe][L][n], phi2, Spe_VPS_XV[spe], 
                               Spe_VPS_RV[spe], Spe_Num_Mesh_VPS[spe]);
    
  	  for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
	    Projector_VNA[spe][L][m][i] += Projector_VNA[spe][L][n][i]*pe[n]*dum0;
	  }
        }

        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
	  Projector_VNA[spe][L][m][i] = phi[m][i] - Projector_VNA[spe][L][m][i];
        }
      
        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++) phi2[i] = Projector_VNA[spe][L][m][i]*Spe_Vna[spe][i];

        dum0 = Int_phi0_phi1(Projector_VNA[spe][L][m], phi2, Spe_VPS_XV[spe], 
                                  Spe_VPS_RV[spe], Spe_Num_Mesh_VPS[spe]);

        /* empty atom */ 
        if (Spe_WhatAtom[spe]==0){
          pe[m] = 0.0;
        }

        else if (fabs(dum0)<1.0e-15){
          pe[m] = 0.0;
        }

        /* normal atom */ 
        else{
          pe[m] = 1.0/dum0;
	}
      }

      /* Renormalization */
        
      for (m=0; m<List_YOUSO[34]; m++){

        dum0 = Int_phi0_phi1(Projector_VNA[spe][L][m], Projector_VNA[spe][L][m],
                           Spe_VPS_XV[spe], Spe_VPS_RV[spe], Spe_Num_Mesh_VPS[spe]);

        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
	  Projector_VNA[spe][L][m][i] = dum0*Projector_VNA[spe][L][m][i];
        }

        if (Spe_WhatAtom[spe]==0){
          VNA_proj_ene[spe][L][m] = 0.0;
        }
        else if (fabs(dum0)<1.0e-15) {
          VNA_proj_ene[spe][L][m] = 0.0;
        }
        else{
          VNA_proj_ene[spe][L][m] = pe[m]/(dum0*dum0);
        }

	/*
        printf("spe=%2d L=%2d m=%2d VNA_proj_ene=%15.12f\n",spe,L,m,VNA_proj_ene[spe][L][m]);
	*/
      }

      /* Calc v*VNL_W2 */

      for (m=0; m<List_YOUSO[34]; m++){
        for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
          Projector_VNA[spe][L][m][i] = Spe_Vna[spe][i]*Projector_VNA[spe][L][m][i];
        }
      }
    } /* for (L=0;...) */

    /* freeing of arrays */
    for (k=0; k<List_YOUSO[34]; k++){
      free(phi[k]);
    }  
    free(phi);

    free(phi2);
    free(pe);
  }  

  /****************************************************
  freeing of arrays:

  double Vcore[List_YOUSO[22]];
  double TmpVNL[2][List_YOUSO[19]];
  double TmpMat[List_YOUSO[37]][List_YOUSO[37]];
  ****************************************************/

  if (remake_headfile==0){
    free(Vcore);

    for (i=0; i<2; i++){
      free(TmpVNL[i]);
    }
    free(TmpVNL);

    for (i=0; i<List_YOUSO[37]; i++){
      free(TmpMat[i]);
    }
    free(TmpMat);
  }

}



double Int_phi0_phi1( double *phi0, double *phi1, double *MXV, double *MRV, int Grid_Num )
{

  /*
  int i;
  double x,dx,xmin,xmax,sum,rp;

  xmin = MXV[0];
  xmax = MXV[Grid_Num-1];
  dx = (xmax - xmin)/(double)OneD_Grid;

  sum = 0.0;
  for (i=0; i<=OneD_Grid; i++){
    x = xmin + (double)i*dx;
    rp = exp(x);
    if (i==0 || i==OneD_Grid)
      sum += 0.5*rp*rp*rp*PhiF(rp,phi0,MRV,Grid_Num)*PhiF(rp,phi1,MRV,Grid_Num);
    else
      sum += rp*rp*rp*PhiF(rp,phi0,MRV,Grid_Num)*PhiF(rp,phi1,MRV,Grid_Num);
  }
  sum = sum*dx;
  */


  int i;
  double rmin,rmax,Sr,Dr,r,sum;

  rmin = MRV[0];
  rmax = MRV[Grid_Num-1];
  Sr = rmax + rmin;
  Dr = rmax - rmin;
  sum = 0.0;
  for (i=0; i<GL_Mesh; i++){
    r = 0.50*(Dr*GL_Abscissae[i] + Sr);
    sum += r*r*GL_Weight[i]*PhiF(r,phi0,MRV,Grid_Num)*PhiF(r,phi1,MRV,Grid_Num);
  }
  sum = 0.50*Dr*sum; 


  return sum;
}





double V_Hart_atom(int spe, double R)
{ 
  int i;
  double xmin,xmax,Sx,Dx,x,rp,dx;
  double Inside,Outside,result;
  
  /****************************************************
              contribution from the inside
  ****************************************************/

  /****************
    simple mesh
  ****************/

  /*
  xmin = Spe_PAO_XV[spe][0];
  xmax = log(R);
  dx = (xmax - xmin)/(double)OneD_Grid;

  Inside = 0.0;
  for (i=0; i<=OneD_Grid; i++){
    x = xmin + (double)i*dx;
    rp = exp(x);
    if (i==0 || i==OneD_Grid)
      Inside = Inside + 0.5*AtomicDenF(spe,rp)*rp*rp*rp;
    else
      Inside = Inside + AtomicDenF(spe,rp)*rp*rp*rp;
  }
  Inside = 4.0*PI*Inside*dx/R; 
  */

  /******************
     Gauss-Legendre
  ******************/

  xmin = Spe_PAO_XV[spe][0];
  xmax = log(R);
  Sx = xmax + xmin;
  Dx = xmax - xmin;
  
  Inside = 0.0;
  for (i=0; i<=(FineGL_Mesh-1); i++){
    x = 0.50*(Dx*FineGL_Abscissae[i] + Sx);
    rp = exp(x);
    Inside = Inside + AtomicDenF(spe,rp)*FineGL_Weight[i]*rp*rp*rp;
  } 
  Inside = 0.50*Dx*Inside;
  Inside = 4.0*PI*Inside/R;

  /****************************************************
              Contribution from the outside
  ****************************************************/

  /****************
    simple mesh
  ****************/

  /*
  xmin = log(R);
  xmax = Spe_PAO_XV[spe][Spe_Num_Mesh_PAO[spe]-1];
  dx = (xmax - xmin)/(double)OneD_Grid;

  Outside = 0.0;
  for (i=0; i<=OneD_Grid; i++){
    x = xmin + (double)i*dx;
    rp = exp(x);
    if (i==0 || i==OneD_Grid)
      Outside = Outside + 0.5*AtomicDenF(spe,rp)*rp*rp;
    else
      Outside = Outside + AtomicDenF(spe,rp)*rp*rp;
  }  
  Outside = 4.0*PI*Outside*dx;
  */

  /******************
     Gauss-Legendre
  ******************/

  xmin = log(R);
  xmax = Spe_PAO_XV[spe][Spe_Num_Mesh_PAO[spe]-1];
  Sx = xmax + xmin;
  Dx = xmax - xmin;
  
  Outside = 0.0;
  for (i=0; i<=(FineGL_Mesh-1); i++){
    x = 0.50*(Dx*FineGL_Abscissae[i] + Sx);
    rp = exp(x);
    Outside = Outside + AtomicDenF(spe,rp)*FineGL_Weight[i]*rp*rp;
  } 
  Outside = 0.50*Dx*Outside;
  Outside = 4.0*PI*Outside;

  result = Inside + Outside;
  return result;
} 


void Set_BasisPara()
{

  int spe,num,L,L0,i,j;
  int Mul0,M0,al,p;
  int ***tmp_index;

  /****************************************************
  allocation of arrays:

  int tmp_index[List_YOUSO[25]+1]
                      [List_YOUSO[24]]
                      [2*(List_YOUSO[25]+1)+1];
  ****************************************************/

  tmp_index = (int***)malloc(sizeof(int**)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    tmp_index[i] = (int**)malloc(sizeof(int*)*List_YOUSO[24]);
    for (j=0; j<List_YOUSO[24]; j++){
      tmp_index[i][j] = (int*)malloc(sizeof(int)*(2*(List_YOUSO[25]+1)+1));
    }
  }

  /* allocation of grobal variables */

  Spe_Specified_Num = (int**)malloc(sizeof(int*)*List_YOUSO[18]); 
  Spe_Trans_Orbital = (int***)malloc(sizeof(int**)*List_YOUSO[18]); 

  for (spe=0; spe<SpeciesNum; spe++){

    /****************************************************
                 Correction of Num_CBasis
                    when Cnt_switch==0
    ****************************************************/

    if (Cnt_switch==0){
      for (L=0; L<=Spe_MaxL_Basis[spe]; L++){
        Spe_Num_CBasis[spe][L] = Spe_Num_Basis[spe][L];   
      }
    }

    /****************************************************
             Total number of bases of the element
                  before the contraction
    ****************************************************/

    num = 0;
    for (L=0; L<=Spe_MaxL_Basis[spe]; L++){
      num += Spe_Num_Basis[spe][L]*(2*L + 1);
    }
    Spe_Total_NO[spe] = num;
    if (List_YOUSO[7]<=num) List_YOUSO[7] = num;

    /****************************************************
             Total number of bases of the element
                  after the contraction
    ****************************************************/

    num = 0;
    for (L=0; L<=Spe_MaxL_Basis[spe]; L++){
      num += Spe_Num_CBasis[spe][L]*(2*L + 1);
    }
    Spe_Total_CNO[spe] = num;

    /****************************************************
        allocation of arrays: 

        int Spe_Specified_Num[List_YOUSO[18]]
                             [Spe_Total_NO[spe]];
        int Spe_Trans_Orbital[List_YOUSO[18]]
                             [Spe_Total_NO[spe]]
                             [List_YOUSO[24]];
    ****************************************************/

    Spe_Specified_Num[spe] = (int*)malloc(sizeof(int)*Spe_Total_NO[spe]); 
    Spe_Trans_Orbital[spe] = (int**)malloc(sizeof(int*)*Spe_Total_NO[spe]); 
    for (i=0; i<Spe_Total_NO[spe]; i++){
      Spe_Trans_Orbital[spe][i] = (int*)malloc(sizeof(int)*List_YOUSO[24]); 
    }

    /****************************************************
      Transformation indices from contracted obrbitals
                 to uncontracted orbitals
    ****************************************************/

    if (remake_headfile==0){

      al = -1;
      for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
        for (Mul0=0; Mul0<Spe_Num_Basis[spe][L0]; Mul0++){
          for (M0=0; M0<=2*L0; M0++){
            al++;
            tmp_index[L0][Mul0][M0] = al;             
          }
        }
      } 

      al = -1;
      for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
        for (Mul0=0; Mul0<Spe_Num_CBasis[spe][L0]; Mul0++){
          for (M0=0; M0<=2*L0; M0++){
            al++;
            Spe_Specified_Num[spe][al] = Spe_Num_Basis[spe][L0];
            for (p=0; p<Spe_Specified_Num[spe][al]; p++){
              Spe_Trans_Orbital[spe][al][p] = tmp_index[L0][p][M0];
            }
          }
        }
      }

    }

    /*
    for (al=0; al<=7; al++){
      printf("al=%2d  %2d\n",al,Specified_Num[wan][al]);
    }
    for (al=0; al<=7; al++){
      for (p=0; p<Specified_Num[wan][al]; p++){
        printf("al=%2d p=%2d  %2d\n",
                al,p,Trans_Orbital[wan][al][p]);
      }
    }
    */

  }

  /****************************************************
     Set of initial coefficients of the contraction
  ****************************************************/

  /*
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatAtom[ct_AN];

    for (al=0; al<Total_CNO[wan]; al++){
      for (p=0; p<Specified_Num[wan][al]; p++){
        CntCoes[ct_AN][al][p] = 0.0;
      }
    }

    al = -1;
    for (L0=0; L0<=MaxL_Basis[wan]; L0++){
      for (Mul0=0; Mul0<Num_CBasis[wan][L0]; Mul0++){
        for (M0=0; M0<=2*L0; M0++){
          al++;
          CntCoes[ct_AN][al][Mul0] = 1.0;
	}
      }
    }
  }
  */

  /****************************************************
  freeing of arrays:

  int tmp_index[List_YOUSO[25]+1]
                      [List_YOUSO[24]]
                      [2*(List_YOUSO[25]+1)+1];
  ****************************************************/

  for (i=0; i<(List_YOUSO[25]+1); i++){
    for (j=0; j<List_YOUSO[24]; j++){
      free(tmp_index[i][j]);
    }
    free(tmp_index[i]);
  }
  free(tmp_index);

}


void Set_Comp2Real()
{
  int i,j,L;

  /* s */

  Comp2Real[0][0][0] = Complex( 1.0, 0.0);

  /* p */

  Comp2Real[1][0][0] = Complex(-1.0/sqrt(2.0), 0.0);
  Comp2Real[1][0][1] = Complex( 0.0,           0.0);
  Comp2Real[1][0][2] = Complex( 1.0/sqrt(2.0), 0.0);

  Comp2Real[1][1][0] = Complex( 0.0, -1.0/sqrt(2.0));
  Comp2Real[1][1][1] = Complex( 0.0,           0.0);
  Comp2Real[1][1][2] = Complex( 0.0, -1.0/sqrt(2.0));

  Comp2Real[1][2][0] = Complex( 0.0, 0.0);
  Comp2Real[1][2][1] = Complex(-1.0, 0.0);
  Comp2Real[1][2][2] = Complex( 0.0, 0.0);

  /* d */

  Comp2Real[2][0][0] = Complex( 0.0,           0.0);
  Comp2Real[2][0][1] = Complex( 0.0,           0.0);
  Comp2Real[2][0][2] = Complex( 1.0,           0.0);
  Comp2Real[2][0][3] = Complex( 0.0,           0.0);
  Comp2Real[2][0][4] = Complex( 0.0,           0.0);

  Comp2Real[2][1][0] = Complex( 1.0/sqrt(2.0), 0.0);
  Comp2Real[2][1][1] = Complex( 0.0,           0.0);
  Comp2Real[2][1][2] = Complex( 0.0,           0.0);
  Comp2Real[2][1][3] = Complex( 0.0,           0.0);
  Comp2Real[2][1][4] = Complex( 1.0/sqrt(2.0), 0.0);

  Comp2Real[2][2][0] = Complex( 0.0, 1.0/sqrt(2.0));
  Comp2Real[2][2][1] = Complex( 0.0,           0.0);
  Comp2Real[2][2][2] = Complex( 0.0,           0.0);
  Comp2Real[2][2][3] = Complex( 0.0,           0.0);
  Comp2Real[2][2][4] = Complex( 0.0,-1.0/sqrt(2.0));

  Comp2Real[2][3][0] = Complex( 0.0,           0.0);
  Comp2Real[2][3][1] = Complex( 1.0/sqrt(2.0), 0.0);
  Comp2Real[2][3][2] = Complex( 0.0,           0.0);
  Comp2Real[2][3][3] = Complex(-1.0/sqrt(2.0), 0.0);
  Comp2Real[2][3][4] = Complex( 0.0,           0.0);

  Comp2Real[2][4][0] = Complex( 0.0,           0.0);
  Comp2Real[2][4][1] = Complex( 0.0, 1.0/sqrt(2.0));
  Comp2Real[2][4][2] = Complex( 0.0,           0.0);
  Comp2Real[2][4][3] = Complex( 0.0, 1.0/sqrt(2.0));
  Comp2Real[2][4][4] = Complex( 0.0,           0.0);

  /* f */

  Comp2Real[3][0][0] = Complex( 0.0,           0.0);
  Comp2Real[3][0][1] = Complex( 0.0,           0.0);
  Comp2Real[3][0][2] = Complex( 0.0,           0.0);
  Comp2Real[3][0][3] = Complex(-1.0,           0.0);
  Comp2Real[3][0][4] = Complex( 0.0,           0.0);
  Comp2Real[3][0][5] = Complex( 0.0,           0.0);
  Comp2Real[3][0][6] = Complex( 0.0,           0.0);

  Comp2Real[3][1][0] = Complex( 0.0,           0.0);
  Comp2Real[3][1][1] = Complex( 0.0,           0.0);
  Comp2Real[3][1][2] = Complex(-1.0/sqrt(2.0), 0.0);
  Comp2Real[3][1][3] = Complex( 0.0,           0.0);
  Comp2Real[3][1][4] = Complex( 1.0/sqrt(2.0), 0.0);
  Comp2Real[3][1][5] = Complex( 0.0,           0.0);
  Comp2Real[3][1][6] = Complex( 0.0,           0.0);

  Comp2Real[3][2][0] = Complex( 0.0,           0.0);
  Comp2Real[3][2][1] = Complex( 0.0,           0.0);
  Comp2Real[3][2][2] = Complex( 0.0,-1.0/sqrt(2.0));
  Comp2Real[3][2][3] = Complex( 0.0,           0.0);
  Comp2Real[3][2][4] = Complex( 0.0,-1.0/sqrt(2.0));
  Comp2Real[3][2][5] = Complex( 0.0,           0.0);
  Comp2Real[3][2][6] = Complex( 0.0,           0.0);

  Comp2Real[3][3][0] = Complex( 0.0,           0.0);
  Comp2Real[3][3][1] = Complex(-1.0/sqrt(2.0), 0.0);
  Comp2Real[3][3][2] = Complex( 0.0,           0.0);
  Comp2Real[3][3][3] = Complex( 0.0,           0.0);
  Comp2Real[3][3][4] = Complex( 0.0,           0.0);
  Comp2Real[3][3][5] = Complex(-1.0/sqrt(2.0), 0.0);
  Comp2Real[3][3][6] = Complex( 0.0,           0.0);

  Comp2Real[3][4][0] = Complex( 0.0,           0.0);
  Comp2Real[3][4][1] = Complex( 0.0,-1.0/sqrt(2.0));
  Comp2Real[3][4][2] = Complex( 0.0,           0.0);
  Comp2Real[3][4][3] = Complex( 0.0,           0.0);
  Comp2Real[3][4][4] = Complex( 0.0,           0.0);
  Comp2Real[3][4][5] = Complex( 0.0, 1.0/sqrt(2.0));
  Comp2Real[3][4][6] = Complex( 0.0,           0.0);

  Comp2Real[3][5][0] = Complex(-1.0/sqrt(2.0), 0.0);
  Comp2Real[3][5][1] = Complex( 0.0,           0.0);
  Comp2Real[3][5][2] = Complex( 0.0,           0.0);
  Comp2Real[3][5][3] = Complex( 0.0,           0.0);
  Comp2Real[3][5][4] = Complex( 0.0,           0.0);
  Comp2Real[3][5][5] = Complex( 0.0,           0.0);
  Comp2Real[3][5][6] = Complex( 1.0/sqrt(2.0), 0.0);

  Comp2Real[3][6][0] = Complex( 0.0,-1.0/sqrt(2.0));
  Comp2Real[3][6][1] = Complex( 0.0,           0.0);
  Comp2Real[3][6][2] = Complex( 0.0,           0.0);
  Comp2Real[3][6][3] = Complex( 0.0,           0.0);
  Comp2Real[3][6][4] = Complex( 0.0,           0.0);
  Comp2Real[3][6][5] = Complex( 0.0,           0.0);
  Comp2Real[3][6][6] = Complex( 0.0,-1.0/sqrt(2.0));

  /* g(4), h(5), i(6), j(7), k(8), l(9), and m(10) */

  for (L=4; L<=YOUSO36; L++){

    for (i=0; i<(2*L+1); i++){
      for (j=0; j<(2*L+1); j++){
	Comp2Real[L][i][j]  = Complex(0.0, 0.0);
      }
    }

    Comp2Real[L][0][L] = Complex(-1.0, 0.0);

    j = -1;
    for (i=1; i<(L*2+1); i=i+4){
      j++; 
      Comp2Real[L][i][L-(2*j+1)]  = Complex(-1.0/sqrt(2.0), 0.0);
      Comp2Real[L][i][L+(2*j+1)]  = Complex( 1.0/sqrt(2.0), 0.0);
    }

    j = 0;
    for (i=3; i<(L*2+1); i=i+4){
      j++;
      Comp2Real[L][i][L-2*j]  = Complex(-1.0/sqrt(2.0), 0.0);
      Comp2Real[L][i][L+2*j]  = Complex(-1.0/sqrt(2.0), 0.0);
    }

    j = -1;
    for (i=2; i<(L*2+1); i=i+4){
      j++;
      Comp2Real[L][i][L-(2*j+1)]  = Complex(0.0, -1.0/sqrt(2.0));
      Comp2Real[L][i][L+(2*j+1)]  = Complex(0.0, -1.0/sqrt(2.0));
    }

    j = 0;
    for (i=4; i<(L*2+1); i=i+4){
      j++;
      Comp2Real[L][i][L-2*j]  = Complex(0.0, -1.0/sqrt(2.0));
      Comp2Real[L][i][L+2*j]  = Complex(0.0,  1.0/sqrt(2.0));
    }
  }


}

void Set_Atom_Weight()
{
  /****************************************************
               Atomic weight (Hydrogen=1.0)
  ****************************************************/

  Atom_Weight[  0] =   1.00;
  Atom_Weight[  1] =   1.00;
  Atom_Weight[  2] =   4.00;
  Atom_Weight[  3] =   7.00;
  Atom_Weight[  4] =   9.00;
  Atom_Weight[  5] =  11.00;
  Atom_Weight[  6] =  12.00;
  Atom_Weight[  7] =  14.00;
  Atom_Weight[  8] =  16.00;
  Atom_Weight[  9] =  19.00;
  Atom_Weight[ 10] =  20.00;
  Atom_Weight[ 11] =  23.00;
  Atom_Weight[ 12] =  24.00;
  Atom_Weight[ 13] =  27.00;
  Atom_Weight[ 14] =  28.00;
  Atom_Weight[ 15] =  31.00;
  Atom_Weight[ 16] =  32.00;
  Atom_Weight[ 17] =  35.00;
  Atom_Weight[ 18] =  40.00;
  Atom_Weight[ 19] =  39.00;
  Atom_Weight[ 20] =  40.00;
  Atom_Weight[ 21] =  45.00;
  Atom_Weight[ 22] =  48.00;
  Atom_Weight[ 23] =  51.00;
  Atom_Weight[ 24] =  52.00;
  Atom_Weight[ 25] =  55.00;
  Atom_Weight[ 26] =  56.00;
  Atom_Weight[ 27] =  59.00;
  Atom_Weight[ 28] =  59.00;
  Atom_Weight[ 29] =  64.00;
  Atom_Weight[ 30] =  65.00;
  Atom_Weight[ 31] =  70.00;
  Atom_Weight[ 32] =  73.00;
  Atom_Weight[ 33] =  75.00;
  Atom_Weight[ 34] =  79.00;
  Atom_Weight[ 35] =  80.00;
  Atom_Weight[ 36] =  84.00;
  Atom_Weight[ 37] =  85.00;
  Atom_Weight[ 38] =  88.00;
  Atom_Weight[ 39] =  89.00;
  Atom_Weight[ 40] =  91.00;
  Atom_Weight[ 41] =  93.00;
  Atom_Weight[ 42] =  96.00;
  Atom_Weight[ 43] =  98.00;
  Atom_Weight[ 44] = 101.00;
  Atom_Weight[ 45] = 103.00;
  Atom_Weight[ 46] = 106.00;
  Atom_Weight[ 47] = 108.00;
  Atom_Weight[ 48] = 112.00;
  Atom_Weight[ 49] = 115.00;
  Atom_Weight[ 50] = 119.00;
  Atom_Weight[ 51] = 122.00;
  Atom_Weight[ 52] = 128.00;
  Atom_Weight[ 53] = 127.00;
  Atom_Weight[ 54] = 131.00;
  Atom_Weight[ 55] = 133.00;
  Atom_Weight[ 56] = 137.00;
  Atom_Weight[ 57] = 139.00;
  Atom_Weight[ 58] = 140.00;
  Atom_Weight[ 59] = 141.00;
  Atom_Weight[ 60] = 144.00;
  Atom_Weight[ 61] = 145.00;
  Atom_Weight[ 62] = 150.00;
  Atom_Weight[ 63] = 152.00;
  Atom_Weight[ 64] = 157.00;
  Atom_Weight[ 65] = 159.00;
  Atom_Weight[ 66] = 163.00;
  Atom_Weight[ 67] = 165.00;
  Atom_Weight[ 68] = 167.00;
  Atom_Weight[ 69] = 169.00;
  Atom_Weight[ 70] = 173.00;
  Atom_Weight[ 71] = 175.00;
  Atom_Weight[ 72] = 178.00;
  Atom_Weight[ 73] = 181.00;
  Atom_Weight[ 74] = 184.00;
  Atom_Weight[ 75] = 186.00;
  Atom_Weight[ 76] = 190.00;
  Atom_Weight[ 77] = 192.00;
  Atom_Weight[ 78] = 195.00;
  Atom_Weight[ 79] = 197.00;
  Atom_Weight[ 80] = 201.00;
  Atom_Weight[ 81] = 204.00;
  Atom_Weight[ 82] = 207.00;
  Atom_Weight[ 83] = 209.00;
  Atom_Weight[ 84] = 209.00;
  Atom_Weight[ 85] = 210.00;
  Atom_Weight[ 86] = 222.00;
  Atom_Weight[ 87] = 223.00;
  Atom_Weight[ 88] = 226.00;
  Atom_Weight[ 89] = 227.00;
  Atom_Weight[ 90] = 232.00;
  Atom_Weight[ 91] = 231.00;
  Atom_Weight[ 92] = 238.00;
  Atom_Weight[ 93] = 237.00;
  Atom_Weight[ 94] = 244.00;
  Atom_Weight[ 95] = 243.00;
  Atom_Weight[ 96] = 247.00;
  Atom_Weight[ 97] = 247.00;
  Atom_Weight[ 98] = 251.00;
  Atom_Weight[ 99] = 252.00;
  Atom_Weight[100] = 257.00;
  Atom_Weight[101] = 258.00;
  Atom_Weight[102] = 259.00;
  Atom_Weight[103] = 260.00;
}

void Set_Atom_Symbol()
{
  strcpy(Atom_Symbol[  0], "E");
  strcpy(Atom_Symbol[  1], "H");
  strcpy(Atom_Symbol[  2], "He");
  strcpy(Atom_Symbol[  3], "Li");
  strcpy(Atom_Symbol[  4], "Be");
  strcpy(Atom_Symbol[  5], "B");
  strcpy(Atom_Symbol[  6], "C");
  strcpy(Atom_Symbol[  7], "N");
  strcpy(Atom_Symbol[  8], "O");
  strcpy(Atom_Symbol[  9], "F");
  strcpy(Atom_Symbol[ 10], "Ne");
  strcpy(Atom_Symbol[ 11], "Na");
  strcpy(Atom_Symbol[ 12], "Mg");
  strcpy(Atom_Symbol[ 13], "Al");
  strcpy(Atom_Symbol[ 14], "Si");
  strcpy(Atom_Symbol[ 15], "P");
  strcpy(Atom_Symbol[ 16], "S");
  strcpy(Atom_Symbol[ 17], "Cl");
  strcpy(Atom_Symbol[ 18], "Ar");
  strcpy(Atom_Symbol[ 19], "K");
  strcpy(Atom_Symbol[ 20], "Ca");
  strcpy(Atom_Symbol[ 21], "Sc");
  strcpy(Atom_Symbol[ 22], "Ti");
  strcpy(Atom_Symbol[ 23], "V");
  strcpy(Atom_Symbol[ 24], "Cr");
  strcpy(Atom_Symbol[ 25], "Mn");
  strcpy(Atom_Symbol[ 26], "Fe");
  strcpy(Atom_Symbol[ 27], "Co");
  strcpy(Atom_Symbol[ 28], "Ni");
  strcpy(Atom_Symbol[ 29], "Cu");
  strcpy(Atom_Symbol[ 30], "Zn");
  strcpy(Atom_Symbol[ 31], "Ga");
  strcpy(Atom_Symbol[ 32], "Ge");
  strcpy(Atom_Symbol[ 33], "As");
  strcpy(Atom_Symbol[ 34], "Se");
  strcpy(Atom_Symbol[ 35], "Br");
  strcpy(Atom_Symbol[ 36], "Kr");
  strcpy(Atom_Symbol[ 37], "Rb");
  strcpy(Atom_Symbol[ 38], "Sr");
  strcpy(Atom_Symbol[ 39], "Y");
  strcpy(Atom_Symbol[ 40], "Zr");
  strcpy(Atom_Symbol[ 41], "Nb");
  strcpy(Atom_Symbol[ 42], "Mo");
  strcpy(Atom_Symbol[ 43], "Tc");
  strcpy(Atom_Symbol[ 44], "Ru");
  strcpy(Atom_Symbol[ 45], "Rh");
  strcpy(Atom_Symbol[ 46], "Pd");
  strcpy(Atom_Symbol[ 47], "Ag");
  strcpy(Atom_Symbol[ 48], "Cd");
  strcpy(Atom_Symbol[ 49], "In");
  strcpy(Atom_Symbol[ 50], "Sn");
  strcpy(Atom_Symbol[ 51], "Sb");
  strcpy(Atom_Symbol[ 52], "Te");
  strcpy(Atom_Symbol[ 53], "I");
  strcpy(Atom_Symbol[ 54], "Xe");
  strcpy(Atom_Symbol[ 55], "Cs");
  strcpy(Atom_Symbol[ 56], "Ba");
  strcpy(Atom_Symbol[ 57], "La");
  strcpy(Atom_Symbol[ 58], "Ce");
  strcpy(Atom_Symbol[ 59], "Pr");
  strcpy(Atom_Symbol[ 60], "Nd");
  strcpy(Atom_Symbol[ 61], "Pm");
  strcpy(Atom_Symbol[ 62], "Sm");
  strcpy(Atom_Symbol[ 63], "Eu");
  strcpy(Atom_Symbol[ 64], "Gd");
  strcpy(Atom_Symbol[ 65], "Tb");
  strcpy(Atom_Symbol[ 66], "Dy");
  strcpy(Atom_Symbol[ 67], "Ho");
  strcpy(Atom_Symbol[ 68], "Er");
  strcpy(Atom_Symbol[ 69], "Tm");
  strcpy(Atom_Symbol[ 70], "Yb");
  strcpy(Atom_Symbol[ 71], "Lu");
  strcpy(Atom_Symbol[ 72], "Hf");
  strcpy(Atom_Symbol[ 73], "Ta");
  strcpy(Atom_Symbol[ 74], "W");
  strcpy(Atom_Symbol[ 75], "Re");
  strcpy(Atom_Symbol[ 76], "Os");
  strcpy(Atom_Symbol[ 77], "Ir");
  strcpy(Atom_Symbol[ 78], "Pt");
  strcpy(Atom_Symbol[ 79], "Au");
  strcpy(Atom_Symbol[ 80], "Hg");
  strcpy(Atom_Symbol[ 81], "Tl");
  strcpy(Atom_Symbol[ 82], "Pb");
  strcpy(Atom_Symbol[ 83], "Bi");
  strcpy(Atom_Symbol[ 84], "Po");
  strcpy(Atom_Symbol[ 85], "At");
  strcpy(Atom_Symbol[ 86], "Rn");
  strcpy(Atom_Symbol[ 87], "Fr");
  strcpy(Atom_Symbol[ 88], "Ra");
  strcpy(Atom_Symbol[ 89], "Ac");
  strcpy(Atom_Symbol[ 90], "Th");
  strcpy(Atom_Symbol[ 91], "Pa");
  strcpy(Atom_Symbol[ 92], "U");
  strcpy(Atom_Symbol[ 93], "Np");
  strcpy(Atom_Symbol[ 94], "Pu");
  strcpy(Atom_Symbol[ 95], "Am");
  strcpy(Atom_Symbol[ 96], "Cm");
  strcpy(Atom_Symbol[ 97], "Bk");
  strcpy(Atom_Symbol[ 98], "Cf");
  strcpy(Atom_Symbol[ 99], "Es");
  strcpy(Atom_Symbol[100], "Fm");
  strcpy(Atom_Symbol[101], "Md");
  strcpy(Atom_Symbol[102], "No");
  strcpy(Atom_Symbol[103], "Lr");
}

void Check_InitDensity()
{

  int ct_AN,cwan,po,wan;
  double charge;

  po = 0;  
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    cwan = WhatSpecies[ct_AN];
    charge = Spe_Core_Charge[cwan] - (InitN_USpin[ct_AN] + InitN_DSpin[ct_AN]);

    if (10e-14<fabs(charge)){
      printf("Invalid values for the initial densities of atom %i\n",ct_AN);
      po++;
    }
  }

  if (po!=0){
    MPI_Finalize();
    exit(1);
  }

  Total_Num_Electrons = 0.0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    Total_Num_Electrons += Spe_Core_Charge[wan];
  }
  Total_Num_Electrons -= system_charge;

}


void Inverse(int n, double **a, double **ia)
{
  /****************************************************
                  LU decomposition
                      0 to n
   NOTE:
     This routine does consider the reduction of rank
  ****************************************************/

  int i,j,k,n3;
  double w;
  double *x,*y;
  double **da;

  /***************************************************
    allocation of arrays: 

    x[n+3]
    y[n+3]
    da[n+3][n+3]
  ***************************************************/

  n3 = n + 3;
  x = (double*)malloc(sizeof(double)*n3);
  y = (double*)malloc(sizeof(double)*n3);

  da = (double**)malloc(sizeof(double*)*n3);
  for (i=0; i<n3; i++){
    da[i] = (double*)malloc(sizeof(double)*n3);
  }

  /* start calc. */

  if (n==-1){
    for (i=0; i<n3; i++){
      for (j=0; j<n3; j++){
	a[i][j] = 0.0;
      }
    }
  }
  else{
    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j] = a[i][j];
      }
    }

    /****************************************************
                     LU factorization
    ****************************************************/

    for (k=0; k<=n-1; k++){
      w = 1.0/a[k][k];
      for (i=k+1; i<=n; i++){
	a[i][k] = w*a[i][k];
	for (j=k+1; j<=n; j++){
	  a[i][j] = a[i][j] - a[i][k]*a[k][j];
	}
      }
    }
    for (k=0; k<=n; k++){

      /****************************************************
                             Ly = b
      ****************************************************/

      for (i=0; i<=n; i++){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=0; j<=i-1; j++){
	  y[i] = y[i] - a[i][j]*y[j];
	}
      }

      /****************************************************
                             Ux = y 
      ****************************************************/

      for (i=n; 0<=i; i--){
	x[i] = y[i];
	for (j=n; (i+1)<=j; j--){
	  x[i] = x[i] - a[i][j]*x[j];
	}
	x[i] = x[i]/a[i][i];
	ia[i][k] = x[i];
      }
    }

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	a[i][j] = da[i][j];
      }
    }
  }

  /***************************************************
    freeing of arrays: 

     x[n+3]
     y[n+3]
     da[n+3][n+3]
  ***************************************************/

  free(x);
  free(y);

  for (i=0; i<n3; i++){
    free(da[i]);
  }
  free(da);
}








static void Generate_PAO_Mixed_Basis(int spe)
{
  int i,j,L;
  double dx,xmin,xmax,rs,r; 

  Spe_WhatAtom[spe] = 0;
  Spe_Num_Mesh_PAO[spe] = Spe_Num_Mesh_PAO[spe-1];

  rcut_FEB *= 1.7; 

  xmin = -7.0;
  xmax = log(1.2*rcut_FEB);
  dx = (xmax - xmin)/Spe_Num_Mesh_PAO[spe];

  for (i=0; i<Spe_Num_Mesh_PAO[spe]; i++){
    Spe_PAO_XV[spe][i] = xmin + (double)i*dx;
    Spe_PAO_RV[spe][i] = exp(Spe_PAO_XV[spe][i]);
    Spe_Atomic_Den[spe][i] = 0.0;
  }

  Spe_Atom_Cut1[spe] = rcut_FEB;
  Spe_PAO_LMAX[spe] = Spe_MaxL_Basis[spe];
  Spe_PAO_Mul[spe] = 1;

  for (L=0; L<=Spe_PAO_LMAX[spe]; L++){
    for (j=0; j<Spe_PAO_Mul[spe]; j++){
      for (i=0; i<Spe_Num_Mesh_PAO[spe]; i++){

        r = Spe_PAO_RV[spe][i];
        rs = r/rcut_FEB;

        if (rs<1.0)
          Spe_PAO_RWF[spe][L][j][i] = (1.0 - 3.0*rs*rs + 2.0*rs*rs*rs)*1.0e+3; 
        else 
          Spe_PAO_RWF[spe][L][j][i] = 0.0;
      }
    }
  }
}







static void Generate_VPS_Mixed_Basis(int spe)
{
  int i,j,L,m;
  double dx,xmin,xmax,rs,r; 

  Spe_Core_Charge[spe] = 0.0;
  Spe_Num_Mesh_VPS[spe] = Spe_Num_Mesh_VPS[spe-1];
  VPS_j_dependency[spe] = 0;
  Spe_Num_RVPS[spe] = 0;
  Spe_Total_VPS_Pro[spe] = 0; 

  xmin = -7.0;
  xmax = log(1.2*rcut_FEB);
  dx = (xmax - xmin)/Spe_Num_Mesh_VPS[spe];

  for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
    Spe_VPS_XV[spe][i] = xmin + (double)i*dx;
    Spe_VPS_RV[spe][i] = exp(Spe_VPS_XV[spe][i]);
    Spe_VH_Atom[spe][i] = 0.0;
    Spe_Vna[spe][i] = 0.0;
  }

  if (ProExpn_VNA==1){
    for (L=0; L<=List_YOUSO[35]; L++){ 
      for (m=0; m<List_YOUSO[34]; m++){
	VNA_proj_ene[spe][L][m] = 0.0;
	for (i=0; i<Spe_Num_Mesh_VPS[spe]; i++){
	  Projector_VNA[spe][L][m][i] = 0.0;
	}
      }
    }
  }

}





void output_structures()
{
  int i,j,Gc_AN,k,itmp;
  double tmp[4],Cxyz[4];
  double lena,lenb,lenc,t1;
  double alpha,beta,gamma;
  double CellV,Cell_Volume;
  char fname1[YOUSO10];
  FILE *fp;

  /********************************************
              making of a CIF file
  *********************************************/
  
  sprintf(fname1,"%s%s.cif",filepath,filename);

  Cross_Product(tv[2],tv[3],tmp);
  CellV = Dot_Product(tv[1],tmp); 
  Cell_Volume = fabs(CellV);

  lena = sqrt(fabs(Dot_Product(tv[1],tv[1]))); 
  lenb = sqrt(fabs(Dot_Product(tv[2],tv[2]))); 
  lenc = sqrt(fabs(Dot_Product(tv[3],tv[3]))); 

  Cross_Product(tv[2],tv[3],tmp);
  rtv[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv[1][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[3],tv[1],tmp);
  rtv[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv[2][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[1],tv[2],tmp);
  rtv[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv[3][3] = 2.0*PI*tmp[3]/CellV;

  /* alpha: angle between b and c in Deg. */

  t1 = Dot_Product(tv[2],tv[3]);
  if (fabs(t1)<1.0e-14) 
    alpha = 90.0;
  else 
    alpha = acos(t1/(lenb*lenc))/PI*180.0;

  /* beta: angle between c and a in Deg. */

  t1 = Dot_Product(tv[3],tv[1]);
  if (fabs(t1)<1.0e-14) 
    beta = 90.0;
  else 
    beta = acos(t1/(lenc*lena))/PI*180.0;

  /* gamma: angle between a and b in Deg. */

  t1 = Dot_Product(tv[1],tv[2]);
  if (fabs(t1)<1.0e-14) 
    gamma = 90.0;
  else 
    gamma = acos(t1/(lena*lenb))/PI*180.0;

  if ((fp = fopen(fname1,"w")) != NULL){

    /* write cell infomations */

    fprintf(fp,"data_%s\n",filename);
    fprintf(fp,"_audit_creation_date              2007-10-11\n");
    fprintf(fp,"_audit_creation_method            'Materials Studio'\n");

    fprintf(fp,"_symmetry_space_group_name_H-M    'P1'\n");
    fprintf(fp,"_symmetry_Int_Tables_number       1\n");
    fprintf(fp,"_symmetry_cell_setting            triclinic\n");

    fprintf(fp,"loop_\n");
    fprintf(fp,"_symmetry_equiv_pos_as_xyz\n");
    fprintf(fp,"  x,y,z\n");

    fprintf(fp,"_cell_length_a%26.4f\n",lena*BohrR);
    fprintf(fp,"_cell_length_b%26.4f\n",lenb*BohrR);
    fprintf(fp,"_cell_length_c%26.4f\n",lenc*BohrR);

    fprintf(fp,"_cell_angle_alpha%24.4f\n",alpha);
    fprintf(fp,"_cell_angle_beta %24.4f\n",beta);
    fprintf(fp,"_cell_angle_gamma%24.4f\n",gamma);

    fprintf(fp,"loop_\n");
    fprintf(fp,"_atom_site_label\n");
    fprintf(fp,"_atom_site_type_symbol\n");
    fprintf(fp,"_atom_site_fract_x\n");
    fprintf(fp,"_atom_site_fract_y\n");
    fprintf(fp,"_atom_site_fract_z\n");
    fprintf(fp,"_atom_site_Uiso_or_equiv\n");
    fprintf(fp,"_atom_site_adp_type\n");
    fprintf(fp,"_atom_site_occupancy\n");

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      /* The zero is taken as the origin of the unit cell. */

      Cxyz[1] = Gxyz[Gc_AN][1];
      Cxyz[2] = Gxyz[Gc_AN][2];
      Cxyz[3] = Gxyz[Gc_AN][3];

      Cell_Gxyz[Gc_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
      Cell_Gxyz[Gc_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
      Cell_Gxyz[Gc_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

      /* The fractional coordinates are kept within 0 to 1. */

      for (i=1; i<=3; i++){

	itmp = (int)Cell_Gxyz[Gc_AN][i]; 

	if (1.0<Cell_Gxyz[Gc_AN][i]){
	  Cell_Gxyz[Gc_AN][i] = fabs(Cell_Gxyz[Gc_AN][i] - (double)itmp);
	}
	else if (Cell_Gxyz[Gc_AN][i]<-1.0e-13){
	  Cell_Gxyz[Gc_AN][i] = fabs(Cell_Gxyz[Gc_AN][i] + (double)(abs(itmp)+1));
	}
      }

      k = WhatSpecies[Gc_AN];
      j = Spe_WhatAtom[k];

      fprintf(fp,"A%-6d%-3s%10.5f%10.5f%10.5f%10.5f  Uiso   1.00\n",
              Gc_AN, 
	      Atom_Symbol[j],
	      Cell_Gxyz[Gc_AN][1],
	      Cell_Gxyz[Gc_AN][2],
	      Cell_Gxyz[Gc_AN][3],
              0.0);
    }

    fclose(fp);
  }

  /********************************************
              making of a XYZ file
  *********************************************/

  sprintf(fname1,"%s%s.xyz",filepath,filename);

  if ((fp = fopen(fname1,"w")) != NULL){

    fprintf(fp,"%d\n\n",atomnum);

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      k = WhatSpecies[Gc_AN];
      j = Spe_WhatAtom[k];

      fprintf(fp,"%4s  %18.14f %18.14f %18.14f\n",
	      Atom_Symbol[j],
	      Gxyz[Gc_AN][1]*BohrR,
	      Gxyz[Gc_AN][2]*BohrR,
	      Gxyz[Gc_AN][3]*BohrR);
    }

    fclose(fp);
  }

}
