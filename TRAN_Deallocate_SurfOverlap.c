#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tran_variables.h"

void TRAN_Deallocate_SurfOverlap(char *position)
{
   int side;
   double *S00, *S01;
   double  **H00, **H01; 
   int k,SpinP_switch;

  if ( strcasecmp(position,"left")==0) {
    side=0;
  } else if ( strcasecmp(position,"right")==0) {
    side=1;
  } 

  SpinP_switch = SpinP_switch_e[side];

  S00 = S00_e[side];
  S01 = S01_e[side];
  H00 = H00_e[side];
  H01 = H01_e[side];

  free( S00 );
  free( S01 );
  free( H00 );
  free( H01  );
  for (k=0;k<=SpinP_switch; k++) {
    free( H00[k] ); 
    free( H01[k] );
  }

}
