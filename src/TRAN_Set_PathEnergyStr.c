#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"


static void TRAN_error_and_exit(
    char *buf
)
{
   printf("%s\n",buf);
   exit(0);
}

/* format:
 *      -1.0 a  0.0 r 0.0  0.1  
 *      auto a  auto r auto auto 
 *     (number|auto)
 */
void TRAN_Set_PathEnergyStr_Square(
   int m,
   char **str,
   double default_relative_ene[4], /* auto values */
   double tran_square_path_ene[4], /* output */
   int    tran_square_path_ene_fix[4]
)
{
   static char *a="auto";
   int i,n,r_a;
   char buf[100];

   printf("TRAN_Set_PathEnergyStr_Square in\n");

   for (i=0;i<2;i++) {

   if (strcasecmp(str[i*2+1],"r")==0) { /* relative */
         r_a=1;
   }
   else if (strcasecmp(str[i*2+1],"a")==0) { /* absolute */
         r_a=0;
   }
   else {
        sprintf(buf,"keyword at %d =<%s> is invalid",2*i+1,str[i*2+1]);
            TRAN_error_and_exit(buf);
   }

   if (strcasecmp(str[i*2],a)==0) {
        tran_square_path_ene[i]=default_relative_ene[i];
        tran_square_path_ene_fix[i]=0;
   }
   else {
       n=sscanf(str[i*2],"%lf",&tran_square_path_ene[i]);
       if (n!=1) {
            sprintf(buf,"keyword at %d =<%s> is invalid",2*i,str[i*2]);
            TRAN_error_and_exit(buf);
       }
       if (r_a==1) {
            tran_square_path_ene_fix[i]=0;
       }
       else {
            tran_square_path_ene_fix[i]=1;
       }
   }

  }

   for (i=2;i<m;i++) {

   if (strcasecmp(str[i+2],a)==0) {
        tran_square_path_ene[i]=default_relative_ene[i];
   }
   else {
       n=sscanf(str[i+2],"%lf",&tran_square_path_ene[i]);
       if (n!=1) {
            sprintf(buf,"keyword at %d =<%s> is invalid",i+2,str[i+2]);
            TRAN_error_and_exit(buf);
       }
   }
   tran_square_path_ene_fix[i]=1;

   }

  printf(" TRAN_Set_PathEnergyStr_Square out \n");

}


void TRAN_Set_PathEnergyStr(
   int tran_integ_pathtype, 
   int m,
   char **str,
   double default_relative_ene[4], /* auto values */
   double tran_square_path_ene[4], /* output */
   int  tran_square_path_ene_fix[4] /* output */
)
{

  if (tran_integ_pathtype==1) {
     TRAN_Set_PathEnergyStr_Square(
     m,
     str,
    default_relative_ene, /* auto values */
    tran_square_path_ene, /* output */
    tran_square_path_ene_fix /* output */
    );
  }
  else if (tran_integ_pathtype==10) {
     TRAN_Set_PathEnergyStr_Square(
     m,
     str,
    default_relative_ene, /* auto values */
    tran_square_path_ene, /* output */
    tran_square_path_ene_fix /* output */
    );
  }
  else {
    printf("tran_integ_pathtype==%d, not supported\n", tran_integ_pathtype);
    exit(0);
  }


}

