#include <stdio.h>
#include <stdlib.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"


void  TRAN_adjust_Ngrid( MPI_Comm comm1, int *Ngrid1,int *Ngrid2, int *Ngrid3)
{
  static int numprocs,myid,ID;
 
  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  if ( Ngrid3_e[0] != Ngrid3_e[1]  ||  Ngrid2_e[0] != Ngrid2_e[1] ) {
    printf("TRAN> internal error, Ngird?_e[0] is different from  Ngrid?_e[1]\n");
    printf("TRAN> Ngrid?[0]=%d %d %d Ngrid?[1]=%d %d %d\n", 
	   Ngrid1_e[0],Ngrid2_e[0], Ngrid3_e[0], 
	   Ngrid1_e[1],Ngrid2_e[1], Ngrid3_e[1] );
    exit(0);

  }

  if (myid==Host_ID){
    printf("<TRAN> adjust Ngrid and Ngrid2 to fit those of electrodes \n");
  }

  *Ngrid2=Ngrid2_e[0];
  *Ngrid3=Ngrid3_e[0];
}
