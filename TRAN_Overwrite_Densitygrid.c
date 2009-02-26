#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_variables.h"
#include "tran_prototypes.h"


/*
 *  [0:TRAN_grid_bound[0]]  and [TRAN_grid_bound[1]:Ngrid3-1]
 *
 *  are from those of electrodes 
 */

void TRAN_Overwrite_Densitygrid(
            MPI_Comm comm1,
            int SpinP_switch,
            int Ngrid1,
            int Ngrid2,
            int Ngrid3,
            int Num_Cells0, 
            int *My_Cell0, 
            int *My_Cell1,
            double **Density_Grid 
)
#define grid_ref(i,j,k)   ( (i)*Ngrid2*Ngrid3+(j)*Ngrid3+(k) )
/*
#define grid_e_ref(i,j,k) ( (i)*Ngrid2*(l3[1]-l3[0]+1)+(j)*(l3[1]-l3[0]+1)+(k)-l3[0] ) 
*/
#define grid_e_ref(i,j,k)  ( ((i)-l1[0]) *Ngrid2*Ngrid3+(j)*Ngrid3+(k) )

{
  int side,l1[2];
  int i,j,k;
  int spin;
  int ie;
  int myid;

  MPI_Comm_rank(comm1,&myid);

  if (myid==Host_ID){
    printf("<TRAN_Overwrite_Densitygrid>\n");
  }

#ifdef DEBUG
  {
    char name[100];
    sprintf(name,"dg0.%d",myid);
    TRAN_Print_Grid_Cell1(name, Num_Cells0,Ngrid2,Ngrid3,My_Cell1, Density_Grid[0]);
  }
#endif

  side = 0;
  l1[0] = 0;
  l1[1] = TRAN_grid_bound[0]; 
  
  for (spin=0; spin<=SpinP_switch; spin++){

    for (i=0; i<Num_Cells0; i++) {

      ie = My_Cell1[i]; 

      if ( l1[0]<=ie && ie<=l1[1] ) {

	for (j=0;j<Ngrid2; j++) {

	  for (k=0; k<Ngrid3;k++) {
	    Density_Grid[spin][ grid_ref(i,j,k) ] = ElectrodeDensity_Grid[side][spin][ grid_e_ref(ie,j,k) ];
	  }

          if (SpinP_switch==0){
	    for (k=0; k<Ngrid3;k++) {
	      Density_Grid[1][ grid_ref(i,j,k) ] = ElectrodeDensity_Grid[side][0][ grid_e_ref(ie,j,k) ];
  	    }
          }
	}
      }
    }

  }  /* spin */

  side = 1;
  l1[0] = TRAN_grid_bound[1];
  l1[1] = Ngrid1-1;
  
  for (spin=0; spin<=SpinP_switch; spin++){

    for (i=0; i<Num_Cells0; i++) {

      ie = My_Cell1[i];

      if ( l1[0]<=ie && ie<=l1[1] ) {
	for (j=0;j<Ngrid2; j++) {

	  for (k=0;k<Ngrid3;k++) {
	    Density_Grid[spin][ grid_ref(i,j,k) ] = ElectrodeDensity_Grid[side][spin][ grid_e_ref(ie,j,k) ];
	  }

          if (SpinP_switch==0){
	    for (k=0;k<Ngrid3;k++) {
	      Density_Grid[1][ grid_ref(i,j,k) ] = ElectrodeDensity_Grid[side][0][ grid_e_ref(ie,j,k) ];
  	    }
	  }

	}
      }
    }
    
  }  /* spin */

#ifdef DEBUG
  {
    char name[100];
    sprintf(name,"dg.%d",myid);
    TRAN_Print_Grid_Cell1(name, Num_Cells0,Ngrid2,Ngrid3,My_Cell1, Density_Grid[0]);
  }
#endif

}
