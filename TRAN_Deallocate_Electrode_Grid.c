#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tran_variables.h"

void TRAN_Deallocate_Electrode_Grid(char *position)
{

  int side,spin;

  if ( strcasecmp(position,"left")==0) {
    side=0;
  } else if ( strcasecmp(position,"right")==0) {
    side=1;
  } 



  for (side=0;side<2;side++) {

  for (spin=0;spin<=SpinP_switch_e[side];spin++) {
    free( ElectrodeDensity_Grid[side][spin] );
  }

    free( ElectrodeDensity_Grid[side] );


    free( ElectrodedVHart_Grid[side] ); 

    free( ElectrodedVHart_Grid_c[side] );

  }

}


