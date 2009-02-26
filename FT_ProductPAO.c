/**********************************************************************
  FT_ProductPAO.c:

     FT_ProductPAO.c is a subroutine to Fourier transform the
     product of two pseudo atomic orbitals.

  Log of FT_ProductPAO.c:

     18/May/2004  Released by T.Ozaki

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

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif



void FT_ProductPAO()
{
  int numprocs,myid,ID,tag=999;
  int count,NumSpe;
  int L,i,j,kj,l,Lmax;
  int Lspe,spe,GL,GL1,Mul1,GL2,Mul2;
  double Sr,Dr,Sk,Dk,kmin,kmax;
  double norm_k,h,dum0;
  double rmin,rmax,r,sum;
  double sj,sy,sjp,syp;
  double RGL[GL_Mesh+2];
  double Tmp_RGL[GL_Mesh+2];
  double SumTmp;
  double tmp0,tmp1;
  double ***GL_PAO;
  double **SphB;
  double *tmp_SphB,*tmp_SphBp;
  double TStime, TEtime;
  /* for MPI */
  MPI_Status stat;
  MPI_Request request;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  dtime(&TStime);

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<FT_ProductPAO>   Fourier transform of product of PAOs\n");

  for (Lspe=0; Lspe<MSpeciesNum; Lspe++){

    spe = Species_Top[myid] + Lspe;

    /* initalize */

    rmin = Spe_VPS_RV[spe][0];
    rmax = Spe_Atom_Cut1[spe] + 0.5;
    Sr = rmax + rmin;
    Dr = rmax - rmin;
    for (i=0; i<GL_Mesh; i++){
      RGL[i] = 0.50*(Dr*GL_Abscissae[i] + Sr);
      Tmp_RGL[i] = RGL[i]*RGL[i]*GL_Weight[i];
    }

    kmin = Radial_kmin;
    kmax = PAO_Nkmax;
    Sk = kmax + kmin;
    Dk = kmax - kmin;

#pragma omp parallel shared(Spe_ProductRF_Bessel,RGL,Tmp_RGL,Dr,Dk,Sk,GL_Abscissae,List_YOUSO,Spe_MaxL_Basis,spe,Spe_Num_Basis)  private(GL_PAO,i,j,GL1,Mul1,GL2,Lmax,Mul2,l,SumTmp,sj,r,SphB,GL,tmp_SphB,tmp_SphBp,OMPID,Nthrds,Nprocs,kj,norm_k)
    {

    /****************************************************
                \int RL * RL' * jl(k*r) r^2 dr 
    ****************************************************/

    /* allocation of GL_PAO */

    GL_PAO = (double***)malloc(sizeof(double**)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      GL_PAO[i] = (double**)malloc(sizeof(double*)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
        GL_PAO[i][j] = (double*)malloc(sizeof(double)*(GL_Mesh + 2));
      }
    }

    /* calculate GL_PAO */

    for (GL1=0; GL1<=Spe_MaxL_Basis[spe]; GL1++){
      for (Mul1=0; Mul1<Spe_Num_Basis[spe][GL1]; Mul1++){
        for (i=0; i<GL_Mesh; i++){
          r = RGL[i];
          GL_PAO[GL1][Mul1][i] = RadialF(spe,GL1,Mul1,r);
        }
      }
    }

    /* allocate arrays */

    SphB = (double**)malloc(sizeof(double*)*(2*Spe_MaxL_Basis[spe]+3));
    for(GL=0; GL<(2*Spe_MaxL_Basis[spe]+3); GL++){ 
      SphB[GL] = (double*)malloc(sizeof(double)*GL_Mesh);
    }

    tmp_SphB  = (double*)malloc(sizeof(double)*(2*Spe_MaxL_Basis[spe]+3));
    tmp_SphBp = (double*)malloc(sizeof(double)*(2*Spe_MaxL_Basis[spe]+3));

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* kj loop */

    for ( kj=OMPID; kj<GL_Mesh; kj+=Nthrds ){

      norm_k = 0.50*(Dk*GL_Abscissae[kj] + Sk);

      /* calculate SphB */

      for (i=0; i<GL_Mesh; i++){

        r = RGL[i];

	Spherical_Bessel(norm_k*r,2*Spe_MaxL_Basis[spe],tmp_SphB,tmp_SphBp);

	for(GL=0; GL<=2*Spe_MaxL_Basis[spe]; GL++){ 
	  SphB[GL][i]  =  tmp_SphB[GL]; 
	}
      }

      /*  \tilde{R}_{L,L',l}  */

      for (GL1=0; GL1<=Spe_MaxL_Basis[spe]; GL1++){
	for (Mul1=0; Mul1<Spe_Num_Basis[spe][GL1]; Mul1++){
	  for (GL2=0; GL2<=Spe_MaxL_Basis[spe]; GL2++){

	    if (GL1<=GL2){

	      Lmax = 2*GL2;

	      for (Mul2=0; Mul2<Spe_Num_Basis[spe][GL2]; Mul2++){
		for(l=0; l<=Lmax; l++){

                  /* Gauss-Legendre quadrature */

                  SumTmp = 0.0;
                  for (i=0; i<GL_Mesh; i++){
                    r = RGL[i];
                    sj = SphB[l][i];

                    SumTmp += Tmp_RGL[i]*sj*GL_PAO[GL1][Mul1][i]*GL_PAO[GL2][Mul2][i];
                  }

                  Spe_ProductRF_Bessel[spe][GL1][Mul1][GL2][Mul2][l][kj] = 0.5*Dr*SumTmp;

                }
  	      }
	    }
	  }
	}
      }
    } /* kj */

    /* free arrays */

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
        free(GL_PAO[i][j]);
      }
      free(GL_PAO[i]);
    }
    free(GL_PAO);

    /* free SphB */

    for(GL=0; GL<(2*Spe_MaxL_Basis[spe]+3); GL++){ 
      free(SphB[GL]);
    }
    free(SphB);

    free(tmp_SphB);
    free(tmp_SphBp);

#pragma omp flush(Spe_ProductRF_Bessel)

    } /* #pragma omp parallel */
  } /* Lspe */

  /****************************************************
         regenerate radial grids in the k-space
         for the MPI calculation
  ****************************************************/

  for (j=0; j<GL_Mesh; j++){
    kmin = Radial_kmin;
    kmax = PAO_Nkmax;
    Sk = kmax + kmin;
    Dk = kmax - kmin;
    norm_k = 0.50*(Dk*GL_Abscissae[j] + Sk);
    GL_NormK[j] = norm_k;
  }

  /***********************************************************
      sending and receiving of Spe_CrudeVNA_Bessel by MPI
  ***********************************************************/

  MPI_Barrier(mpi_comm_level1);

  for (ID=0; ID<Num_Procs2; ID++) {
    NumSpe = Species_End[ID] - Species_Top[ID] + 1;
    for (Lspe=0; Lspe<NumSpe; Lspe++){
      spe = Species_Top[ID] + Lspe;

      for (GL1=0; GL1<=Spe_MaxL_Basis[spe]; GL1++){
	for (Mul1=0; Mul1<Spe_Num_Basis[spe][GL1]; Mul1++){
	  for (GL2=0; GL2<=Spe_MaxL_Basis[spe]; GL2++){

	    if (GL1<=GL2){

	      Lmax = 2*GL2;

	      for (Mul2=0; Mul2<Spe_Num_Basis[spe][GL2]; Mul2++){
		for(l=0; l<=Lmax; l++){
		  MPI_Bcast(&Spe_ProductRF_Bessel[spe][GL1][Mul1][GL2][Mul2][l][0],
			    GL_Mesh,MPI_DOUBLE,ID,mpi_comm_level1);
                  MPI_Barrier(mpi_comm_level1);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /***********************************************************
                         elapsed time
  ***********************************************************/

  dtime(&TEtime);

  /*
  printf("myid=%2d Elapsed Time (s) = %15.12f\n",myid,TEtime-TStime);
  MPI_Finalize();
  exit(0);
  */

}
