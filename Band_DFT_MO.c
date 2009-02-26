/**********************************************************************
  Band_DFT_MO.c:

     Band_DFT_MO.c is a subroutine to calculate wave functions
     at given k-points for the file output.

  Log of Band_DFT_MO.c:

     15/May/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

static void Band_DFT_MO_Col(
                      int nkpoint, double **kpoint,
                      int SpinP_switch, 
                      double *****nh,
                      double *****ImNL,
                      double ****CntOLP);

static void Band_DFT_MO_NonCol(
                      int nkpoint, double **kpoint,
                      int SpinP_switch, 
                      double *****nh,
                      double *****ImNL,
                      double ****CntOLP);


void Band_DFT_MO( int nkpoint, double **kpoint,
                  int SpinP_switch, 
                  double *****nh,
                  double *****ImNL,
                  double ****CntOLP)
{
  if (SpinP_switch==0 || SpinP_switch==1){
    Band_DFT_MO_Col( nkpoint, kpoint, SpinP_switch, nh, ImNL, CntOLP);
  }
  else if (SpinP_switch==3){
    Band_DFT_MO_NonCol( nkpoint, kpoint, SpinP_switch, nh, ImNL, CntOLP);
  }
}



static void Band_DFT_MO_Col(
                      int nkpoint, double **kpoint,
                      int SpinP_switch, 
                      double *****nh,
                      double *****ImNL,
                      double ****CntOLP)
{
  int i,j,k,l,n,wan;
  int *MP;
  int i1,j1,po,spin,n1;
  int num2,RnB,l1,l2,l3,kloop;
  int ct_AN,h_AN,wanA,tnoA,wanB,tnoB;
  int GA_AN,Anum,nhomos,nlumos;
  int ii,ij,ik;
  int num0,num1,mul,m,wan1,Gc_AN;
  int LB_AN,GB_AN,Bnum;
  double time0;
  double snum_i,snum_j,snum_k,k1,k2,k3,sum,sumi,Num_State,FermiF;
  double x,Dnum,Dnum2,AcP,ChemP_MAX,ChemP_MIN,EV_cut0;
  double **ko,*M1,***EIGEN;
  double *koS;
  dcomplex ***H,**S,***C;
  dcomplex Ctmp1,Ctmp2;
  double ****Dummy_ImNL;
  double ***Ctmp;
  double u2,v2,uv,vu;
  double dum,sumE,kRn,si,co;
  double Resum,ResumE,Redum,Redum2,Imdum;
  double TStime,TEtime,SiloopTime,EiloopTime;
  double FermiEps = 1.0e-14;
  double x_cut = 30.0;
  double OLP_eigen_cut=Threshold_OLP_Eigen;
  char *Name_Angular[Supported_MaxL+1][2*(Supported_MaxL+1)+1];
  char *Name_Multiple[20];
  char file_EV[YOUSO10];
  FILE *fp_EV;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID) printf("\nBand_DFT_MO start\n");fflush(stdout);
  
  dtime(&TStime);

  /****************************************************
   Allocation

   int       MP[List_YOUSO[1]]
   double    ko[List_YOUSO[23]][n+1]
   double    koS[n+1];
   double    EIGEN[List_YOUSO[33]][List_YOUSO[23]][n+1]
   dcomplex  H[List_YOUSO[23]][n+1][n+1]
   dcomplex  S[n+1][n+1]
   double    M1[n+1]
   dcomplex  C[List_YOUSO[23]][n+1][n+1]
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  
  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n  = n + Spe_Total_CNO[wanA];
  }

  ko = (double**)malloc(sizeof(double*)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    ko[i] = (double*)malloc(sizeof(double)*(n+1));
  }

  koS = (double*)malloc(sizeof(double)*(n+1));

  EIGEN = (double***)malloc(sizeof(double**)*List_YOUSO[33]);
  for (i=0; i<List_YOUSO[33]; i++){
    EIGEN[i] = (double**)malloc(sizeof(double*)*List_YOUSO[23]);
    for (j=0; j<List_YOUSO[23]; j++){
      EIGEN[i][j] = (double*)malloc(sizeof(double)*(n+1));
    }
  }

  H = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    H[i] = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (j=0; j<n+1; j++){
      H[i][j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }
  }

  S = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (i=0; i<n+1; i++){
    S[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }


  M1 = (double*)malloc(sizeof(double)*(n+1));

  C = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    C[i] = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (j=0; j<n+1; j++){
      C[i][j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }
  }

  Ctmp = (double***)malloc(sizeof(double**)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    Ctmp[i] = (double**)malloc(sizeof(double*)*(n+1));
    for (j=0; j<n+1; j++){
      Ctmp[i][j] = (double*)malloc(sizeof(double)*(n+1));
    }
  }

  /* no spin-orbit coupling */
  if (SO_switch==0){
    Dummy_ImNL = (double****)malloc(sizeof(double***)*1);
    Dummy_ImNL[0] = (double***)malloc(sizeof(double**)*1);
    Dummy_ImNL[0][0] = (double**)malloc(sizeof(double*)*1);
    Dummy_ImNL[0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  dtime(&SiloopTime);

  /*****************************************************
         Solve eigenvalue problem at each k-point
  *****************************************************/

  for (kloop=0; kloop<nkpoint; kloop++){

    if (myid==Host_ID) printf("kpoint=%i /%i\n",kloop+1,nkpoint);

    k1 = kpoint[kloop][1];
    k2 = kpoint[kloop][2];
    k3 = kpoint[kloop][3];

    Overlap_Band(Host_ID,CntOLP,S,MP,k1,k2,k3);

    if (myid==Host_ID){

      n = S[0][0].r;
      EigenBand_lapack(S,ko[0],n,1);

      if (3<=level_stdout){
        printf("  kloop %i, k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       kloop,k1,k2,k3);
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

    for (spin=0; spin<=SpinP_switch; spin++){

      if (SO_switch==0)
        Hamiltonian_Band(Host_ID, nh[spin], H[spin],MP,k1,k2,k3);
      else if (SO_switch==1)
        Hamiltonian_Band(Host_ID, nh[spin], H[spin],MP,k1,k2,k3);

      if (myid==Host_ID){

        /* transpose S */

	for (i1=1; i1<=n; i1++){
	  for (j1=i1+1; j1<=n; j1++){
	    Ctmp1 = S[i1][j1];
	    Ctmp2 = S[j1][i1];
	    S[i1][j1] = Ctmp2;
	    S[j1][i1] = Ctmp1;
	  }
	}

        /****************************************************
                        M1 * U^t * H * U * M1
        ****************************************************/

        /* H * U * M1 */
 
        for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n; j1++){

	    sum  = 0.0;
            sumi = 0.0;

	    for (l=1; l<=n; l++){
	      sum  += H[spin][i1][l].r*S[j1][l].r - H[spin][i1][l].i*S[j1][l].i;
  	      sumi += H[spin][i1][l].r*S[j1][l].i + H[spin][i1][l].i*S[j1][l].r;
	    }

  	    C[spin][j1][i1].r = sum;
	    C[spin][j1][i1].i = sumi;
	  }
        }     

        /* M1 * U^+ H * U * M1 */

        for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n; j1++){

	    sum  = 0.0;
            sumi = 0.0;

	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*C[spin][j1][l].r + S[i1][l].i*C[spin][j1][l].i;
	      sumi +=  S[i1][l].r*C[spin][j1][l].i - S[i1][l].i*C[spin][j1][l].r;
	    }

            H[spin][i1][j1].r = sum;
            H[spin][i1][j1].i = sumi;
	  }
        }     

        /* H to C */

        for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n; j1++){
            C[spin][i1][j1] = H[spin][i1][j1];
	  }
        }

	/* penalty for ill-conditioning states */

	EV_cut0 = Threshold_OLP_Eigen;

	for (i1=1; i1<=n; i1++){

	  if (koS[i1]<EV_cut0){
	    C[spin][i1][i1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  }
 
	  /* cutoff the interaction between the ill-conditioned state */
 
	  if (1.0e+3<C[spin][i1][i1].r){
	    for (j1=1; j1<=n; j1++){
	      C[spin][i1][j1] = Complex(0.0,0.0);
	      C[spin][j1][i1] = Complex(0.0,0.0);
	    }
	    C[spin][i1][i1].r = 1.0e+4;
	  }
	}

        n1 = n;
        EigenBand_lapack(C[spin],ko[spin],n1,1);

        for (i1=1; i1<=n1; i1++){
          EIGEN[kloop][spin][i1] = ko[spin][i1];
        }

        /****************************************************
            transformation to the original eigen vectors.
                  NOTE JRCAT-244p and JAIST-2122p 
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
	    Ctmp1 = C[spin][i1][j1];
	    Ctmp2 = C[spin][j1][i1];
	    C[spin][i1][j1] = Ctmp2;
	    C[spin][j1][i1] = Ctmp1;
	  }
	}

        /* shift */
	for (j1=1; j1<=n; j1++){
	  for (l=n; 1<=l; l--){
   	    C[spin][j1][l] = C[spin][j1][l];
	  }
	}

        for (i1=1; i1<=n; i1++){
          for (j1=1; j1<=n; j1++){
            H[spin][i1][j1].r = 0.0;
            H[spin][i1][j1].i = 0.0;
          }
        }

        for (i1=1; i1<=n; i1++){
  	  for (j1=1; j1<=n1; j1++){
	    sum  = 0.0;
            sumi = 0.0;
	    for (l=1; l<=n; l++){
              sum  +=  S[i1][l].r*C[spin][j1][l].r - S[i1][l].i*C[spin][j1][l].i;
              sumi +=  S[i1][l].r*C[spin][j1][l].i + S[i1][l].i*C[spin][j1][l].r;
  	    } 
	    H[spin][i1][j1].r = sum;
	    H[spin][i1][j1].i = sumi;
	  }
        }

        /* find HOMO from eigenvalues */

        for (i1=1; i1<=n1; i1++){
          x = (ko[spin][i1] - ChemP)*Beta;
          if (x<=-x_cut) x = -x_cut;
          if (x_cut<=x)  x = x_cut;
          if (SpinP_switch==0) FermiF = 2.0/(1.0 + exp(x));
          else                 FermiF = 1.0/(1.0 + exp(x));
          if      (SpinP_switch==0 && 1.0<FermiF) Bulk_HOMO[kloop][spin] = i1;
          else if (SpinP_switch==1 && 0.5<FermiF) Bulk_HOMO[kloop][spin] = i1;
        }      

        if (SpinP_switch==0 && 2<=level_stdout){
          printf("k1=%7.3f k2=%7.3f k3=%7.3f  HOMO = %2d\n",
                  k1,k2,k3,Bulk_HOMO[kloop][0]);
        }
        else if (SpinP_switch==1 && 2<=level_stdout){
          printf("k1=%7.3f k2=%7.3f k3=%7.3f  HOMO for up-spin   = %2d\n",
                  k1,k2,k3,Bulk_HOMO[kloop][0]);
          printf("k1=%7.3f k2=%7.3f k3=%7.3f  HOMO for down-spin = %2d\n",
                  k1,k2,k3,Bulk_HOMO[kloop][1]);
        }

      } /* if (myid==Host_ID) */
    }   /* spin */

    /****************************************************
     MPI:

     n1
     Bulk_HOMO
     H
    ****************************************************/

    MPI_Bcast(&n1, 1, MPI_INT, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Bulk_HOMO[kloop][0], 2, MPI_INT, Host_ID, mpi_comm_level1);

    /* H[][][].r to C */
    for (spin=0; spin<=SpinP_switch; spin++){
      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          Ctmp[spin][i1][j1] = H[spin][i1][j1].r; 
	}
      }
    }     

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i1=1; i1<=n; i1++){
         MPI_Bcast(&Ctmp[spin][i1][0], n+1, MPI_DOUBLE, Host_ID, mpi_comm_level1);
      }
    }  

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          C[spin][i1][j1].r = Ctmp[spin][i1][j1]; 
	}
      }
    }     

    /* H[][][].i to C.i */
    for (spin=0; spin<=SpinP_switch; spin++){
      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          Ctmp[spin][i1][j1] = H[spin][i1][j1].i; 
	}
      }
    }     

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i1=1; i1<=n; i1++){
         MPI_Bcast(&Ctmp[spin][i1][0], n+1, MPI_DOUBLE, Host_ID, mpi_comm_level1);
      }
    }  

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          C[spin][i1][j1].i = Ctmp[spin][i1][j1]; 
	}
      }
    }     

    /****************************************************
        LCAO coefficients are stored for calculating
                 values of MOs on grids
    ****************************************************/

    nhomos = num_HOMOs;
    nlumos = num_LUMOs;

    if (SpinP_switch==0){
      if ( (Bulk_HOMO[kloop][0]-nhomos+1)<1 ) nhomos = Bulk_HOMO[kloop][0];
      if ( (Bulk_HOMO[kloop][0]+nlumos)>n1 )  nlumos = n1 - Bulk_HOMO[kloop][0];
    }
    else if (SpinP_switch==1){
      if ( (Bulk_HOMO[kloop][0]-nhomos+1)<1 ) nhomos = Bulk_HOMO[kloop][0];
      if ( (Bulk_HOMO[kloop][1]-nhomos+1)<1 ) nhomos = Bulk_HOMO[kloop][1];
      if ( (Bulk_HOMO[kloop][0]+nlumos)>n1 )  nlumos = n1 - Bulk_HOMO[kloop][0];
      if ( (Bulk_HOMO[kloop][1]+nlumos)>n1 )  nlumos = n1 - Bulk_HOMO[kloop][1];
    }

    /* HOMOs */
    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<nhomos; j++){
        j1 = Bulk_HOMO[kloop][spin] - j;
        for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
          wanA = WhatSpecies[GA_AN];
          tnoA = Spe_Total_CNO[wanA];
          Anum = MP[GA_AN];
          for (i=0; i<tnoA; i++){
            HOMOs_Coef[kloop][spin][j][GA_AN][i].r = C[spin][Anum+i][j1].r;
            HOMOs_Coef[kloop][spin][j][GA_AN][i].i = C[spin][Anum+i][j1].i;
          }
        }
      }        
    }

    /* LUMOs */
    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<nlumos; j++){
        j1 = Bulk_HOMO[kloop][spin] + 1 + j;
        for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
          wanA = WhatSpecies[GA_AN];
          tnoA = Spe_Total_CNO[wanA];
          Anum = MP[GA_AN];
          for (i=0; i<tnoA; i++){
            LUMOs_Coef[kloop][spin][j][GA_AN][i].r = C[spin][Anum+i][j1].r;
            LUMOs_Coef[kloop][spin][j][GA_AN][i].i = C[spin][Anum+i][j1].i;
          }
        }
      }
    }

    Bulk_Num_HOMOs[kloop] = nhomos;
    Bulk_Num_LUMOs[kloop] = nlumos;

    /****************************************************
                          Output
    ****************************************************/

    if (myid==Host_ID){

      strcpy(file_EV,".EV");
      fnjoint(filepath,filename,file_EV);

      if ((fp_EV = fopen(file_EV,"a")) != NULL){

#ifdef xt3
        setvbuf(fp_EV,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        if (kloop==0){

  	  fprintf(fp_EV,"\n");
	  fprintf(fp_EV,"***********************************************************\n");
	  fprintf(fp_EV,"***********************************************************\n");
	  fprintf(fp_EV,"        Eigenvalues (Hartree) and LCAO coefficients        \n");
	  fprintf(fp_EV,"        at the k-points specified in the input file.       \n");
	  fprintf(fp_EV,"***********************************************************\n");
	  fprintf(fp_EV,"***********************************************************\n");
	}

	k1 = kpoint[kloop][1];
	k2 = kpoint[kloop][2];
	k3 = kpoint[kloop][3];

	fprintf(fp_EV,"\n\n");
	fprintf(fp_EV,"   # of k-point = %i\n",kloop);
	fprintf(fp_EV,"   k1=%10.5f k2=%10.5f k3=%10.5f\n\n",k1,k2,k3);
	fprintf(fp_EV,"   Chemical Potential (Hartree) = %18.14f\n",ChemP);

	if (SpinP_switch==0){
	  fprintf(fp_EV,"   HOMO = %i\n\n",Bulk_HOMO[kloop][0]);
	}
	else if (SpinP_switch==1){
	  fprintf(fp_EV,"   HOMO for up-spin   = %i\n",  Bulk_HOMO[kloop][0]);
	  fprintf(fp_EV,"   HOMO for down-spin = %i\n\n",Bulk_HOMO[kloop][1]);
	}

	fprintf(fp_EV,"   Real (Re) and imaginary (Im) parts of LCAO coefficients\n\n");

	num0 = 4;
	num1 = n/num0 + 1*(n%num0!=0);

	for (spin=0; spin<=SpinP_switch; spin++){
	  for (i=1; i<=num1; i++){

	    fprintf(fp_EV,"\n");

	    for (i1=-2; i1<=0; i1++){

	      fprintf(fp_EV,"                     ");

	      for (j=1; j<=num0; j++){

		j1 = num0*(i-1) + j;

		if (j1<=n){ 

		  if (i1==-2){
		    if (spin==0){
		      fprintf(fp_EV,"  %4d (U)",j1);
		      fprintf(fp_EV,"          ");
		    }
		    else if (spin==1){
		      fprintf(fp_EV,"  %4d (D)",j1);
		      fprintf(fp_EV,"          ");
		    }
		  }

		  else if (i1==-1){
		    fprintf(fp_EV,"  %8.4f",EIGEN[kloop][spin][j1]);
		    fprintf(fp_EV,"          ");
		  }

		  else if (i1==0){
		    fprintf(fp_EV,"     Re   ");
		    fprintf(fp_EV,"     Im   ");
		  }
		}
	      }
	      fprintf(fp_EV,"\n");
	      if (i1==-1)  fprintf(fp_EV,"\n");
	      if (i1==0)   fprintf(fp_EV,"\n");
	    }

	    Name_Angular[0][0] = "s          ";
	    Name_Angular[1][0] = "px         ";
	    Name_Angular[1][1] = "py         ";
	    Name_Angular[1][2] = "pz         ";
	    Name_Angular[2][0] = "d3z^2-r^2  ";
	    Name_Angular[2][1] = "dx^2-y^2   ";
	    Name_Angular[2][2] = "dxy        ";
	    Name_Angular[2][3] = "dxz        ";
	    Name_Angular[2][4] = "dyz        ";
	    Name_Angular[3][0] = "f5z^2-3r^2 ";
	    Name_Angular[3][1] = "f5xz^2-xr^2";
	    Name_Angular[3][2] = "f5yz^2-yr^2";
	    Name_Angular[3][3] = "fzx^2-zy^2 ";
	    Name_Angular[3][4] = "fxyz       ";
	    Name_Angular[3][5] = "fx^3-3*xy^2";
	    Name_Angular[3][6] = "f3yx^2-y^3 ";
	    Name_Angular[4][0] = "g1         ";
	    Name_Angular[4][1] = "g2         ";
	    Name_Angular[4][2] = "g3         ";
	    Name_Angular[4][3] = "g4         ";
	    Name_Angular[4][4] = "g5         ";
	    Name_Angular[4][5] = "g6         ";
	    Name_Angular[4][6] = "g7         ";
	    Name_Angular[4][7] = "g8         ";
	    Name_Angular[4][8] = "g9         ";

	    Name_Multiple[0] = "0";
	    Name_Multiple[1] = "1";
	    Name_Multiple[2] = "2";
	    Name_Multiple[3] = "3";
	    Name_Multiple[4] = "4";
	    Name_Multiple[5] = "5";

	    i1 = 1; 

	    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

	      wan1 = WhatSpecies[Gc_AN];
            
	      for (l=0; l<=Supported_MaxL; l++){
		for (mul=0; mul<Spe_Num_CBasis[wan1][l]; mul++){
		  for (m=0; m<(2*l+1); m++){

		    if (l==0 && mul==0 && m==0)
		      fprintf(fp_EV,"%4d %3s %s %s", 
			      Gc_AN,SpeName[wan1],Name_Multiple[mul],Name_Angular[l][m]);
		    else
		      fprintf(fp_EV,"         %s %s", 
			      Name_Multiple[mul],Name_Angular[l][m]);

		    for (j=1; j<=num0; j++){

		      j1 = num0*(i-1) + j;

		      if (0<i1 && j1<=n){
			fprintf(fp_EV,"  %8.5f",C[spin][i1][j1].r);
			fprintf(fp_EV,"  %8.5f",C[spin][i1][j1].i);
		      }
		    }
		    fprintf(fp_EV,"\n");

		    i1++;
		  }
		}
	      }
	    }

	  }
	}

        /* close the file */ 
	fclose(fp_EV);
      }
      else{
	printf("Failure of saving the EV file.\n");
	fclose(fp_EV);
      }
    } /* if (myid==Host_ID) */
  }  /* kloop */

  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);

  for (i=0; i<List_YOUSO[23]; i++){
    free(ko[i]);
  }
  free(ko);

  free(koS);

  for (i=0; i<List_YOUSO[33]; i++){
    for (j=0; j<List_YOUSO[23]; j++){
      free(EIGEN[i][j]);
    }
    free(EIGEN[i]);
  }
  free(EIGEN);  


  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<n+1; j++){
      free(H[i][j]);
    }
    free(H[i]);
  }
  free(H);  


  for (i=0; i<n+1; i++){
    free(S[i]);
  }
  free(S);

  free(M1);

  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<n+1; j++){
      free(C[i][j]);
    }
    free(C[i]);
  }
  free(C);

  for (i=0; i<List_YOUSO[23]; i++){
    for (j=0; j<n+1; j++){
      free(Ctmp[i][j]);
    }
    free(Ctmp[i]);
  }
  free(Ctmp);

  /* no spin-orbit coupling */
  if (SO_switch==0){
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  dtime(&TEtime);
}




static void Band_DFT_MO_NonCol(
                      int nkpoint, double **kpoint,
                      int SpinP_switch, 
                      double *****nh,
                      double *****ImNL,
                      double ****CntOLP)
{
  int i,j,k,l,n,wan,m,ii1,jj1,n2;
  int *MP;
  int i1,j1,po,spin,n1;
  int num2,RnB,l1,l2,l3,kloop;
  int ct_AN,h_AN,wanA,tnoA,wanB,tnoB;
  int GA_AN,Anum,nhomos,nlumos;
  int ii,ij,ik;
  int wan1,mul,Gc_AN,num0,num1;
  double time0;
  int LB_AN,GB_AN,Bnum;
  double snum_i,snum_j,snum_k,k1,k2,k3,sum,sumi,Num_State,FermiF;
  double x,Dnum,Dnum2,AcP,ChemP_MAX,ChemP_MIN;
  double *ko,*M1,**EIGEN;
  double *koS;
  double EV_cut0;
  dcomplex **H,**S,**C;
  double **Ctmp;
  double *****Dummy_ImNL;
  double u2,v2,uv,vu;
  double dum,sumE,kRn,si,co;
  double Resum,ResumE,Redum,Redum2,Imdum;
  double TStime,TEtime,SiloopTime,EiloopTime;
  double FermiEps = 1.0e-14;
  double x_cut = 30.0;
  double OLP_eigen_cut=Threshold_OLP_Eigen;

  char *Name_Angular[Supported_MaxL+1][2*(Supported_MaxL+1)+1];
  char *Name_Multiple[20];
  char file_EV[YOUSO10];
  FILE *fp_EV;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID) printf("\nBand_DFT_MO start\n");fflush(stdout);
  
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
   double    EIGEN[List_YOUSO[33]][n2]
   dcomplex  H[n2][n2]
   dcomplex  S[n2][n2]
   double    M1[n2]
   dcomplex  C[n2][n2]
   double    Ctmp[n2][n2]
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  ko = (double*)malloc(sizeof(double)*n2);
  koS = (double*)malloc(sizeof(double)*(n+1));

  EIGEN = (double**)malloc(sizeof(double*)*List_YOUSO[33]);
  for (i=0; i<List_YOUSO[33]; i++){
    EIGEN[i] = (double*)malloc(sizeof(double)*n2);
  }

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

  Ctmp = (double**)malloc(sizeof(double*)*n2);
  for (j=0; j<n2; j++){
    Ctmp[j] = (double*)malloc(sizeof(double)*n2);
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

  dtime(&SiloopTime);

  /*****************************************************
         Solve eigenvalue problem at each k-point
  *****************************************************/

  for (kloop=0; kloop<nkpoint; kloop++){

    printf("kpoint=%i /%i\n",kloop+1,nkpoint);

    k1 = kpoint[kloop][1];
    k2 = kpoint[kloop][2];
    k3 = kpoint[kloop][3];

    Overlap_Band(Host_ID,CntOLP, S, MP, k1, k2, k3);

    if (myid==Host_ID){

      n = S[0][0].r;
      EigenBand_lapack(S,ko,n,1);

      if (2<=level_stdout){
        printf("  kloop %i, k1 k2 k3 %10.6f %10.6f %10.6f\n",
	       kloop,k1,k2,k3);
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
                      M1 * U^t * H * U * M1
      ****************************************************/

      /* H * U * M1 */

      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=n; j1++){

	  for (m=0; m<=1; m++){

	    sum  = 0.0;
	    sumi = 0.0;
	    for (l=1; l<=n; l++){
	      sum  += (H[i1][l+m*n].r*S[l][j1].r
                       - H[i1][l+m*n].i*S[l][j1].i)*M1[j1];
	      sumi += (H[i1][l+m*n].r*S[l][j1].i
                       + H[i1][l+m*n].i*S[l][j1].r)*M1[j1];
	    }

	    jj1 = 2*j1 - 1 + m;

	    C[i1][jj1].r = sum;
	    C[i1][jj1].i = sumi;
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
	    for (l=1; l<=n; l++){
	      sum  +=  M1[i1]*(S[l][i1].r*C[l+m*n][j1].r +
			       S[l][i1].i*C[l+m*n][j1].i );
	      sumi +=  M1[i1]*(S[l][i1].r*C[l+m*n][j1].i -
			       S[l][i1].i*C[l+m*n][j1].r );
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
      EigenBand_lapack(C, ko, n1,1);

      for (i1=1; i1<=2*n; i1++){
        for (j1=1; j1<=n1; j1++){
          H[i1][j1] = C[i1][j1];
        }
      }

      for (i1=1; i1<=n1; i1++){
        EIGEN[kloop][i1] = ko[i1];
      }

      /****************************************************
          Transformation to the original eigen vectors.
                NOTE JRCAT-244p and JAIST-2122p 
      ****************************************************/

      for (i1=1; i1<=2*n; i1++){
        for (j1=1; j1<=2*n; j1++){
          C[i1][j1].r = 0.0;
          C[i1][j1].i = 0.0;
        }
      }

      for (m=0; m<=1; m++){
	for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n1; j1++){

	    sum  = 0.0; 
	    sumi = 0.0;

	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*M1[l]*H[2*(l-1)+1+m][j1].r
		     - S[i1][l].i*M1[l]*H[2*(l-1)+1+m][j1].i;
	      sumi +=  S[i1][l].r*M1[l]*H[2*(l-1)+1+m][j1].i
		     + S[i1][l].i*M1[l]*H[2*(l-1)+1+m][j1].r;
	    } 
	    C[i1+m*n][j1].r = sum;
	    C[i1+m*n][j1].i = sumi;
	  }
	}
      }

      /* find HOMO from eigenvalues */

      for (i1=1; i1<=n1; i1++){
        x = (ko[i1] - ChemP)*Beta;
        if (x<=-x_cut) x = -x_cut;
        if (x_cut<=x)  x = x_cut;
        FermiF = 1.0/(1.0 + exp(x));
        if (0.5<FermiF) Bulk_HOMO[kloop][0] = i1;
      }      

      if (2<=level_stdout){
        printf("k1=%7.3f k2=%7.3f k3=%7.3f  HOMO = %2d\n",
                k1,k2,k3,Bulk_HOMO[kloop][0]);
      }

    } /* if (myid==Host_ID) */

    /****************************************************
     MPI:

     n1
     Bulk_HOMO
     C
    ****************************************************/

    MPI_Bcast(&n1, 1, MPI_INT, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Bulk_HOMO[kloop][0], 2, MPI_INT, Host_ID, mpi_comm_level1);

    /* C[][].r */
    for (i1=1; i1<=2*n; i1++){
      for (j1=1; j1<=2*n; j1++){
        Ctmp[i1][j1] = C[i1][j1].r; 
      }
    }

    for (i1=1; i1<=2*n; i1++){
      MPI_Bcast(&Ctmp[i1][0], n2, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    }

    for (i1=1; i1<=2*n; i1++){
      for (j1=1; j1<=2*n; j1++){
        C[i1][j1].r = Ctmp[i1][j1]; 
      }
    }

    /* C[][].i */
    for (i1=1; i1<=2*n; i1++){
      for (j1=1; j1<=2*n; j1++){
        Ctmp[i1][j1] = C[i1][j1].i; 
      }
    }

    for (i1=1; i1<=2*n; i1++){
      MPI_Bcast(&Ctmp[i1][0], n2, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    }

    for (i1=1; i1<=2*n; i1++){
      for (j1=1; j1<=2*n; j1++){
        C[i1][j1].i = Ctmp[i1][j1]; 
      }
    }

    /****************************************************
        LCAO coefficients are stored for calculating
                 values of MOs on grids
    ****************************************************/

    nhomos = num_HOMOs;
    nlumos = num_LUMOs;

    if ( (Bulk_HOMO[kloop][0]-nhomos+1)<1 ) nhomos = Bulk_HOMO[kloop][0];
    if ( (Bulk_HOMO[kloop][0]+nlumos)>n1 )  nlumos = n1 - Bulk_HOMO[kloop][0];

    /* HOMOs */

    for (j=0; j<nhomos; j++){
      j1 = Bulk_HOMO[kloop][0] - j;
      for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
        wanA = WhatSpecies[GA_AN];
        tnoA = Spe_Total_CNO[wanA];
        Anum = MP[GA_AN];
        for (i=0; i<tnoA; i++){
          HOMOs_Coef[kloop][0][j][GA_AN][i].r = C[Anum+i][j1].r;
          HOMOs_Coef[kloop][0][j][GA_AN][i].i = C[Anum+i][j1].i;
          HOMOs_Coef[kloop][1][j][GA_AN][i].r = C[Anum+i+n][j1].r;
          HOMOs_Coef[kloop][1][j][GA_AN][i].i = C[Anum+i+n][j1].i;
        }
      }
    }        

    /* LUMOs */
    for (j=0; j<nlumos; j++){
      j1 = Bulk_HOMO[kloop][0] + 1 + j;
      for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
        wanA = WhatSpecies[GA_AN];
        tnoA = Spe_Total_CNO[wanA];
        Anum = MP[GA_AN];
        for (i=0; i<tnoA; i++){
          LUMOs_Coef[kloop][0][j][GA_AN][i].r = C[Anum+i][j1].r;
          LUMOs_Coef[kloop][0][j][GA_AN][i].i = C[Anum+i][j1].i;
          LUMOs_Coef[kloop][1][j][GA_AN][i].r = C[Anum+i+n][j1].r;
          LUMOs_Coef[kloop][1][j][GA_AN][i].i = C[Anum+i+n][j1].i;
        }
      }
    }

    Bulk_Num_HOMOs[kloop] = nhomos;
    Bulk_Num_LUMOs[kloop] = nlumos;

  }  /* kloop */

  /****************************************************
                        Output
  ****************************************************/

  if (myid==Host_ID){

    strcpy(file_EV,".EV");
    fnjoint(filepath,filename,file_EV);

    if ((fp_EV = fopen(file_EV,"a")) != NULL){

#ifdef xt3
      setvbuf(fp_EV,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp_EV,"\n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"        Eigenvalues (Hartree) and LCAO coefficients        \n");
      fprintf(fp_EV,"        at the k-points specified in the input file.       \n");
      fprintf(fp_EV,"***********************************************************\n");
      fprintf(fp_EV,"***********************************************************\n");

      for (kloop=0; kloop<nkpoint; kloop++){
        k1 = kpoint[kloop][1];
        k2 = kpoint[kloop][2];
        k3 = kpoint[kloop][3];
        fprintf(fp_EV,"\n\n");
        fprintf(fp_EV,"   # of k-point = %i\n",kloop);
        fprintf(fp_EV,"   k1=%10.5f k2=%10.5f k3=%10.5f\n\n",k1,k2,k3);
        fprintf(fp_EV,"   Chemical Potential (Hartree) = %18.14f\n",ChemP);
        fprintf(fp_EV,"   HOMO = %i\n\n",Bulk_HOMO[kloop][0]);

	fprintf(fp_EV,"   Real (Re) and imaginary (Im) parts of LCAO coefficients\n\n");

	num0 = 2;
	num1 = 2*n/num0 + 1*((2*n)%num0!=0);
  
	for (i=1; i<=num1; i++){

	  fprintf(fp_EV,"\n");

	  for (i1=-2; i1<=0; i1++){

	    fprintf(fp_EV,"                     ");

	    for (j=1; j<=num0; j++){

	      j1 = num0*(i-1) + j;

	      if (j1<=2*n){ 

		if (i1==-2){
 	          fprintf(fp_EV," %4d",j1);
		  fprintf(fp_EV,"                                   ");
		}

		else if (i1==-1){
		  fprintf(fp_EV,"   %8.5f",EIGEN[kloop][j1]);
		  fprintf(fp_EV,"                             ");
		}

		else if (i1==0){
		  fprintf(fp_EV,"     Re(U)");
		  fprintf(fp_EV,"     Im(U)");
		  fprintf(fp_EV,"     Re(D)");
		  fprintf(fp_EV,"     Im(D)");
		}
	      }
	    }
	    fprintf(fp_EV,"\n");
	    if (i1==-1)  fprintf(fp_EV,"\n");
	    if (i1==0)   fprintf(fp_EV,"\n");
	  }

	  Name_Angular[0][0] = "s          ";
	  Name_Angular[1][0] = "px         ";
	  Name_Angular[1][1] = "py         ";
	  Name_Angular[1][2] = "pz         ";
	  Name_Angular[2][0] = "d3z^2-r^2  ";
	  Name_Angular[2][1] = "dx^2-y^2   ";
	  Name_Angular[2][2] = "dxy        ";
	  Name_Angular[2][3] = "dxz        ";
	  Name_Angular[2][4] = "dyz        ";
	  Name_Angular[3][0] = "f5z^2-3r^2 ";
	  Name_Angular[3][1] = "f5xz^2-xr^2";
	  Name_Angular[3][2] = "f5yz^2-yr^2";
	  Name_Angular[3][3] = "fzx^2-zy^2 ";
	  Name_Angular[3][4] = "fxyz       ";
	  Name_Angular[3][5] = "fx^3-3*xy^2";
	  Name_Angular[3][6] = "f3yx^2-y^3 ";
	  Name_Angular[4][0] = "g1         ";
	  Name_Angular[4][1] = "g2         ";
	  Name_Angular[4][2] = "g3         ";
	  Name_Angular[4][3] = "g4         ";
	  Name_Angular[4][4] = "g5         ";
	  Name_Angular[4][5] = "g6         ";
	  Name_Angular[4][6] = "g7         ";
	  Name_Angular[4][7] = "g8         ";
	  Name_Angular[4][8] = "g9         ";

	  Name_Multiple[0] = "0";
	  Name_Multiple[1] = "1";
	  Name_Multiple[2] = "2";
	  Name_Multiple[3] = "3";
	  Name_Multiple[4] = "4";
	  Name_Multiple[5] = "5";

	  i1 = 1; 

	  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

	    wan1 = WhatSpecies[Gc_AN];
            
	    for (l=0; l<=Supported_MaxL; l++){
	      for (mul=0; mul<Spe_Num_CBasis[wan1][l]; mul++){
		for (m=0; m<(2*l+1); m++){

		  if (l==0 && mul==0 && m==0)
		    fprintf(fp_EV,"%4d %3s %s %s", 
			    Gc_AN,SpeName[wan1],Name_Multiple[mul],Name_Angular[l][m]);
		  else
		    fprintf(fp_EV,"         %s %s", 
			    Name_Multiple[mul],Name_Angular[l][m]);

		  for (j=1; j<=num0; j++){

		    j1 = num0*(i-1) + j;

		    if (0<i1 && j1<=2*n){
		      fprintf(fp_EV,"  %8.5f",C[i1][j1].r);
		      fprintf(fp_EV,"  %8.5f",C[i1][j1].i);
		      fprintf(fp_EV,"  %8.5f",C[i1+n][j1].r);
		      fprintf(fp_EV,"  %8.5f",C[i1+n][j1].i);
		    }
		  }

		  fprintf(fp_EV,"\n");
		  if (i1==-1)  fprintf(fp_EV,"\n");
		  if (i1==0)   fprintf(fp_EV,"\n");

		  i1++;
		}
	      }
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
  } /* if (myid==Host_ID) */

  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);
  free(ko);
  free(koS);

  for (i=0; i<List_YOUSO[33]; i++){
    free(EIGEN[i]);
  }
  free(EIGEN);

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

  for (j=0; j<n2; j++){
    free(Ctmp[j]);
  }
  free(Ctmp);

  /* non-spin-orbit coupling and non-LDA+U */  
  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0
      && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
    free(Dummy_ImNL[0][0][0][0]);
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  dtime(&TEtime);
}




