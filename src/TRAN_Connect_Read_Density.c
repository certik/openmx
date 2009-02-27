#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tran_prototypes.h"


#define LEN 256


void TRAN_Connect_Read_Density(
    char *filename,
    int SpinP_switch,
    int Ngrid1, 
    int Ngrid2, 
    int Ngrid3,
           /* output */
    double *ChemP,
    double *minE, 
    double *dVHart_Grid,
    double **Vpot_Grid,
    double **Density_Grid
)
{
   FILE *fp;
   char buf[LEN];

   int i_Ngrid1, i_Ngrid2, i_Ngrid3;
   int i_SpinP_switch; 
   double i_ChemP, i_minE;
   int i1,i2,i3;
   int GN;

   if ( (fp=fopen(filename,"r"))==NULL ) {
     printf("TRAN_Connect_Read_Density> Error can not open the file, %s\n",filename);
     exit(10); 
   }

   fgets(buf,LEN,fp);
   sscanf(buf,"%d", &i_SpinP_switch);
     if (SpinP_switch != i_SpinP_switch) {
         printf("Error, SpinP_switch is differnet\n");
         exit(10);
     }
  
   fgets(buf,LEN,fp);
   fgets(buf,LEN,fp);
     sscanf(buf,"%d %d %d",&i_Ngrid1, &i_Ngrid2, &i_Ngrid3);

   if (Ngrid1!=i_Ngrid1 || Ngrid2 !=i_Ngrid2 || Ngrid3 != i_Ngrid3 ) {
      printf("Error, grid is different\n");
      printf("Ngrid=%d %d %d\n",Ngrid1,Ngrid2,Ngrid3);
      printf("i_Ngrid=%d %d %d\n",i_Ngrid1,i_Ngrid2,i_Ngrid3);
      exit(10);
   }

     /* read tv */
   fgets(buf,LEN,fp);
   fgets(buf,LEN,fp);
   fgets(buf,LEN,fp);

   
   fgets(buf,LEN,fp);
   fgets(buf,LEN,fp);
    sscanf(buf,"%lf",&i_ChemP);

   fgets(buf,LEN,fp);
    sscanf(buf,"%lf",&i_minE);

      printf("ChemP=%lf, minE=%lf\n",i_ChemP,i_minE);
        *ChemP = i_ChemP;
        *minE = i_minE;

   fgets(buf,LEN,fp);

  for (i1=0; i1<Ngrid1; i1++){
    for (i2=0; i2<Ngrid2; i2++){
      for (i3=0; i3<Ngrid3; i3++){
        fgets(buf,LEN,fp);
        GN = i1*Ngrid2*Ngrid3 + i2*Ngrid3 + i3;
        sscanf(buf,"%lf \n",&dVHart_Grid[GN]);
      }
    }
  }
          printf("the first and the last of dVHart_Grid=%le %le\n",
               dVHart_Grid[0],dVHart_Grid[Ngrid1*Ngrid2*Ngrid3-1]);

  fgets(buf,LEN,fp);

  for (i1=0; i1<Ngrid1; i1++){
    for (i2=0; i2<Ngrid2; i2++){
      for (i3=0; i3<Ngrid3; i3++){
        GN = i1*Ngrid2*Ngrid3 + i2*Ngrid3 + i3;
        fgets(buf,LEN,fp);
        if ( SpinP_switch == 0 ){
          sscanf(buf,"%lf",&Vpot_Grid[0][GN]);
        }
        else{
          sscanf(buf,"%lf %lf",&Vpot_Grid[0][GN],&Vpot_Grid[1][GN]);
        }
      }
    }
  }
       if ( SpinP_switch==0 ) {
          printf("the first and the last of dVpot_Grid=%le %le\n",
                    Vpot_Grid[0][0],  Vpot_Grid[0][Ngrid1*Ngrid2*Ngrid3-1] );
       }
       else {
          printf("the first and the last of dVpot_Grid=%le %le\n",
                    Vpot_Grid[0][0],  Vpot_Grid[1][Ngrid1*Ngrid2*Ngrid3-1] );
       }


  fgets(buf,LEN,fp);

  for (i1=0; i1<Ngrid1; i1++){
    for (i2=0; i2<Ngrid2; i2++){
      for (i3=0; i3<Ngrid3; i3++){
        GN = i1*Ngrid2*Ngrid3 + i2*Ngrid3 + i3;
        fgets(buf,LEN,fp);
        if ( SpinP_switch == 0 ){
          sscanf(buf,"%lf",&Density_Grid[0][GN]);
        }
        else{
          sscanf(buf,"%lf %lf",&Density_Grid[0][GN],&Density_Grid[1][GN]);
        }
      }
    }
  }
       if ( SpinP_switch==0 ) {
          printf("the first and the last of Density_Grid=%le %le\n",
                    Density_Grid[0][0],  Density_Grid[0][Ngrid1*Ngrid2*Ngrid3-1] );
       }
       else {
          printf("the first and the last of Density_Grid=%le %le\n",
                    Density_Grid[0][0],  Density_Grid[1][Ngrid1*Ngrid2*Ngrid3-1] );
       }


   
  fclose(fp);

}

