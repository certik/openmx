#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HARTREE 27.2113845

static void compare_and_print_error(char *buf,char *str,int val1,int val2)
{

  if ( val1 != val2 ) {
    printf("ERROR %s %d %d\n",str,val1,val2);
    printf("buf=<%s>\n",buf);
    exit(0);
  }

}

#define LEN 256

void TRAN_Connect_Read_Hamiltonian(
    char *filename, 
    int SpinP_switch, 
    int *WhatSpecies,
    int *FNAN, 
    int **natn,
    int **ncn, 
    int *Spe_Total_CNO, 
  /* output */
    double *****H, 
    double ****OLP
)
{

   FILE *fp;
   char buf[LEN];

   int atomnum,Catomnum, Latomnum,Ratomnum;
   int spin;
   int ct_AN, wan1,TNO1,i_fnan, h_AN, Gh_AN, Rn, wan2,TNO2;
   int i,j;
   int i_ct_AN, i_h_AN, i_Gh_AN, i_Rn, i_TNO1, i_TNO2; 
   int i_SpinP_switch;


   printf("TRAN_Connect_Read_Hamiltonian file=%s in\n",filename);

   if  (  (fp=fopen(filename,"r"))==NULL ) {
     printf("can not open a file = %s\n",filename);
     exit(0);
   }


   fgets(buf,LEN,fp);
       sscanf(buf,"%d",&atomnum);
   fgets(buf,LEN,fp);
       sscanf(buf,"%d",&Catomnum);
   fgets(buf,LEN,fp);
       sscanf(buf,"%d",&Latomnum);
   fgets(buf,LEN,fp);
       sscanf(buf,"%d",&Ratomnum);
   fgets(buf,LEN,fp);
       sscanf(buf,"%d",&i_SpinP_switch);
   compare_and_print_error(buf,"i_SpinP_switch,SpinP_switch",i_SpinP_switch,SpinP_switch);


  for (spin=0; spin<=SpinP_switch; spin++){
  /*  fprintf(fp,"          Kohn-Sham Hamiltonian: spin=%i\n",spin); */
    fgets(buf,LEN,fp);
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      /* fprintf(fp,"%i     =FNAN[ct_AN] \n",FNAN[ct_AN]); */
      fgets(buf,LEN,fp);
           sscanf(buf,"%d",&i_fnan);

           compare_and_print_error(buf,"i_fnan,FNAN[ct_AN]",i_fnan,FNAN[ct_AN]);

      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
        Gh_AN = natn[ct_AN][h_AN];
        Rn = ncn[ct_AN][h_AN];
        wan2 = WhatSpecies[Gh_AN];
        TNO2 = Spe_Total_CNO[wan2];
        /*fprintf(fp,"%i     =glbal index \n",ct_AN);
        *fprintf(fp,"%i     =local index \n",h_AN);
        *fprintf(fp,"%i     =grobal \n",Gh_AN);
        *fprintf(fp,"%i     =Rn \n",Rn);
        *fprintf(fp,"%i       %i     =orbitalnum \n",TNO1,TNO2);
        */
        fgets(buf,LEN,fp);
        printf("globalindex: %s",buf);
            sscanf(buf,"%d", &i_ct_AN);
            compare_and_print_error(buf,"i_ct_AN,ct_AN",i_ct_AN,ct_AN);

        fgets(buf,LEN,fp);
        printf("localindex: %s",buf);

            sscanf(buf,"%d", &i_h_AN);
            compare_and_print_error(buf,"i_h_AN,h_AN",i_h_AN,h_AN);

        fgets(buf,LEN,fp);
        printf("global %s",buf);

            sscanf(buf,"%d", &i_Gh_AN);
            compare_and_print_error(buf,"i_Gh_AN,Gh_AN",i_Gh_AN,Gh_AN);

        fgets(buf,LEN,fp);
        printf("Rn: %s",buf);

            sscanf(buf,"%d", &i_Rn);
            compare_and_print_error(buf,"i_Rn,Rn",i_Rn,Rn);

        fgets(buf,LEN,fp);
        printf("TNO1,TNO2: %s",buf);

            sscanf(buf,"%d %d", &i_TNO1,&i_TNO2);
            compare_and_print_error(buf,"i_TNO1,TNO1",i_TNO1,TNO1);
            compare_and_print_error(buf,"i_TNO2,TNO2",i_TNO2,TNO2);


        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO2; j++){
            /* fprintf(fp,"%28.17f ",H[spin][ct_AN][h_AN][i][j]); 
             * fprintf(fp,"\n");
             */
             fgets(buf,LEN,fp);
               sscanf(buf,"%lf",&H[spin][ct_AN][h_AN][i][j]);
             H[spin][ct_AN][h_AN][i][j] *= HARTREE;
          }
        }
      }
    }
  }



/*  fprintf(fp,"          Overlap matrix\n"); */
  fgets(buf,LEN,fp);
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      /*    fprintf(fp,"%i     =FNAN[ct_AN] \n",FNAN[ct_AN]); */
      fgets(buf,LEN,fp);
        sscanf(buf,"%d",&i_fnan);
    for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
      Gh_AN = natn[ct_AN][h_AN];
      Rn = ncn[ct_AN][h_AN];
        wan2 = WhatSpecies[Gh_AN];
        TNO2 = Spe_Total_CNO[wan2];
       /* fprintf(fp,"%i     =glbal index \n",ct_AN);
        * fprintf(fp,"%i     =local index \n",h_AN);
        * fprintf(fp,"%i     =grobal \n",Gh_AN);
        * fprintf(fp,"%i     =Rn \n",Rn);
        * fprintf(fp,"%i       %i     =orbitalnum \n",TNO1,TNO2);
        */
        fgets(buf,LEN,fp);
           sscanf(buf,"%d", &ct_AN);
        fgets(buf,LEN,fp);
           sscanf(buf,"%d", &h_AN);
        fgets(buf,LEN,fp);
           sscanf(buf,"%d", &Gh_AN);
        fgets(buf,LEN,fp);
           sscanf(buf,"%d", &Rn);
        fgets(buf,LEN,fp);
           sscanf(buf,"%d %d", &TNO1, &TNO2);

      for (i=0; i<TNO1; i++){
        for (j=0; j<TNO2; j++){
      /*    fprintf(fp,"%28.17f ",OLP[ct_AN][h_AN][i][j]); 
       *   fprintf(fp,"\n");
       */
         fgets(buf,LEN,fp);
           sscanf(buf,"%lf",&OLP[ct_AN][h_AN][i][j]);
        }
      }
    }
  }


  fclose(fp);

  printf("TRAN_Connect_Read_Hamiltonian end\n");

}
