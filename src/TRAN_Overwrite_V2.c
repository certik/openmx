#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"


/*
    Grid(kx,ky,z)
     => copy =>
    ReV2(kx,ky,z)
    ImV2(kx,ky,z)

    z=[0:l3boundary[0]], [l3boundary[1]:Ngrid3-1] 

*/

void  TRAN_Overwrite_V2(
   int Ngrid1,
   int Ngrid2,
   int Ngrid3,
   int l1boundary[2],
   int My_NGrid2_Poisson, 
   int Start2, 
   dcomplex *Grid[2],
   double ***ReV2,  /* output */
   double ***ImV2 ) /* output */
/*#define Grid_e_ref(i,j,k) ( (i)*Ngrid2*(l3[1]-l3[0]+1)+(j)*(l3[1]-l3[0]+1)+(k)-l3[0] )  */
#define Grid_e_ref(i,j,k) ( ( (i)-l1[0] )*Ngrid2*Ngrid3+ (j)*Ngrid3 + (k) )
{
  int i,j,k,l1[2];
  int je;

   int side;

/*left*/
  l1[0]= 0;
  l1[1]= l1boundary[0];
  side=0;
   for (i=l1[0];i<=l1[1];i++){
     for (j=0;j<My_NGrid2_Poisson;j++) {
       je = j+ Start2;
       for (k=0;k<Ngrid3;k++) {
           ReV2[j][i][k] = Grid[side][ Grid_e_ref(i,je,k) ].r;
           ImV2[j][i][k] = Grid[side][ Grid_e_ref(i,je,k) ].i;
       }
     }
   }



/*right */
  l1[0]= l1boundary[1];
  l1[1]= Ngrid1-1;
  side=1;
   for (i=l1[0];i<=l1[1];i++){
     for (j=0;j<My_NGrid2_Poisson;j++) {
       je = j+Start2;
       for (k=0; k<Ngrid3;k++) {
           ReV2[j][i][k] = Grid[side][ Grid_e_ref(i,je,k) ].r;
           ImV2[j][i][k] = Grid[side][ Grid_e_ref(i,je,k) ].i;
       }
     }
   }

}

