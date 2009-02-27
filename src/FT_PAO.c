/**********************************************************************
  FT_PAO.c:

     FT_PAO.c is a subroutine to Fourier transform pseudo atomic 
     orbitals.

  Log of FT_PAO.c:

     15/Sep/2002  Released by T.Ozaki

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



void FT_PAO()
{
  int numprocs,myid,ID,tag=999;
  int count,NumSpe;
  int i,kj,num_k;
  int Lspe,spe,GL,Mul;
  double dk,norm_k,h;
  double rmin,rmax,r,r2,sum;
  double sy,sjp,syp;
  double Sr,Dr,dum0;
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

  if (myid==Host_ID) printf("<FT_PAO>          Fourier transform of pseudo atomic orbitals\n");

  for (Lspe=0; Lspe<MSpeciesNum; Lspe++){

    spe = Species_Top[myid] + Lspe;

    num_k = Ngrid_NormK;
    dk = PAO_Nkmax/(double)num_k;
    rmin = Spe_PAO_RV[spe][0];
    rmax = Spe_Atom_Cut1[spe] + 0.5;
    h = (rmax - rmin)/(double)OneD_Grid;

    /* kj loop */

#pragma omp parallel shared(num_k,dk,spe,rmin,rmax,h,Spe_PAO_RV,Spe_Atom_Cut1,OneD_Grid,Spe_MaxL_Basis,Spe_Num_Basis,Spe_RF_Bessel)  private(OMPID,Nthrds,Nprocs,kj,norm_k,i,r,r2,tmp_SphB,tmp_SphBp,GL,SphB,Mul,sum)
    {

    /* allocate SphB */

    SphB = (double**)malloc(sizeof(double*)*(Spe_MaxL_Basis[spe]+3));
    for(GL=0; GL<(Spe_MaxL_Basis[spe]+3); GL++){ 
      SphB[GL] = (double*)malloc(sizeof(double)*(OneD_Grid+1));
    }

    tmp_SphB  = (double*)malloc(sizeof(double)*(Spe_MaxL_Basis[spe]+3));
    tmp_SphBp = (double*)malloc(sizeof(double)*(Spe_MaxL_Basis[spe]+3));

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for ( kj=OMPID; kj<num_k; kj+=Nthrds ){

      norm_k = (double)kj*dk;

      /* calculate SphB */
 
      for (i=0; i<=OneD_Grid; i++){

        r = rmin + (double)i*h;

	Spherical_Bessel(norm_k*r,Spe_MaxL_Basis[spe],tmp_SphB,tmp_SphBp);

        r2 = r*r;
	for(GL=0; GL<=Spe_MaxL_Basis[spe]; GL++){ 
	  SphB[GL][i] = tmp_SphB[GL]*r2; 
	}
      }

      for(GL=0; GL<=Spe_MaxL_Basis[spe]; GL++){ 
	SphB[GL][0] *= 0.5;
	SphB[GL][OneD_Grid] *= 0.5;
      }

      /* loof for GL and Mul */

      for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++){
	for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){

	  /****************************************************
                        \int jL(k*r)RL r^2 dr 
	  ****************************************************/

	  /* trapezoidal rule */

          sum = 0.0;

          for (i=0; i<=OneD_Grid; i++){
            r = rmin + (double)i*h;
            sum += RadialF(spe,GL,Mul,r)*SphB[GL][i];
          } 
          sum = sum*h;

          Spe_RF_Bessel[spe][GL][Mul][kj] = sum;

	} /* Mul */
      } /* GL */
    } /* kj */

    /* free SphB */

    for(GL=0; GL<(Spe_MaxL_Basis[spe]+3); GL++){ 
      free(SphB[GL]);
    }
    free(SphB);

    free(tmp_SphB);
    free(tmp_SphBp);

    } /* #pragma omp parallel */

    /*
    for ( kj=0; kj<num_k; kj++ ){
      printf("kj=%4d %25.23f %25.23f\n",kj,Spe_RF_Bessel[spe][0][0][kj],Spe_RF_Bessel[spe][1][0][kj]); 
    }
    */


  } /* Lspe */

  /*  
  MPI_Finalize();
  exit(0);
  */

  /****************************************************
          generate radial grids in the k-space
  ****************************************************/

  dk = PAO_Nkmax/(double)Ngrid_NormK;
  for (i=0; i<Ngrid_NormK; i++){
    NormK[i] = (double)i*dk;
  }

  /***********************************************************
        sending and receiving of Spe_RF_Bessel by MPI
  ***********************************************************/

  for (ID=0; ID<Num_Procs2; ID++){
    NumSpe = Species_End[ID] - Species_Top[ID] + 1;
    for (Lspe=0; Lspe<NumSpe; Lspe++){
      spe = Species_Top[ID] + Lspe;
      for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++){
	for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
	  MPI_Bcast(&Spe_RF_Bessel[spe][GL][Mul][0],
		    List_YOUSO[15],MPI_DOUBLE,ID, mpi_comm_level1);
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



