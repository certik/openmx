#include <stdio.h>
#include <string.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"

void TRAN_Credit(MPI_Comm comm1)
{

   int myid;

   MPI_Comm_rank(comm1,&myid);

   if (myid==Host_ID) {

      printf("\n***********************************************************\n"); 
      printf("***********************************************************\n"); 
      printf(" Welcome to NEGF extension\n");
      printf(" The original fortran code was written by Hisashi Kondo\n");
      printf(" as the Frontier Simulation Software for Industrial Science\n");
      printf(" project, and the C code was written by Hiori Kino as an\n");
      printf(" extension of OpenMX, based on the original fortran code,\n");
      printf(" and a further modification was made by Taisuke Ozaki.\n");
      printf("***********************************************************\n"); 
      printf("***********************************************************\n\n\n\n"); 
   }


}
