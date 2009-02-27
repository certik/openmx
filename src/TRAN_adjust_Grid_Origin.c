#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef nompi
#include "mimic_mpi.h" 
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"


void TRAN_adjust_Grid_Origin( MPI_Comm comm1, double Grid_Origin[4])
{
  int po;
  int i;
  static int numprocs,myid,ID;
 
  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  /* test */
  po=0;
  for (i=1;i<=3;i++) {
    if ( 1.0e-14<fabs(Grid_Origin_e[0][i]-Grid_Origin_e[1][i]) ) po++;
  }
  if (po>0) {
    printf("Grid_Origin_e[0] != Grid_Origin_e[1]\n");
    printf("Grid_Origin_e[0]=%lf %lf %lf\n",
            Grid_Origin_e[0][1], Grid_Origin_e[0][2], Grid_Origin_e[0][3]);
    printf("Grid_Origin_e[1]=%lf %lf %lf\n",
            Grid_Origin_e[1][1], Grid_Origin_e[1][2], Grid_Origin_e[1][3]);
    exit(0);
  }

  if (myid==Host_ID){
    printf("Grid_Origin changed\n");
  }

  Grid_Origin[1] = Grid_Origin_e[0][1];  
  Grid_Origin[2] = Grid_Origin_e[0][2];
  Grid_Origin[3] = Grid_Origin_e[0][3];

}
