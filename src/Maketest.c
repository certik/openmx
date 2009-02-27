/**********************************************************************
  Maketest.c:

     Maketest.c is a subroutine to generate *.out files which will be
     used to check whether OpenMX runs normally on many platforms or not.

  Log of Maketest.c:

     25/Oct/2004  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*  end stat section */
#include "openmx_common.h"
#include "Inputtools.h"
 
#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif


void Maketest(char *mode, int argc, char *argv[]) 
{
  FILE *fp0,*fp1,*fp2;
  static int Num_DatFiles,i;

  static char fname0[YOUSO10];
  static char fname1[YOUSO10];
  static char fname_dat[YOUSO10];
  static char fname_dat2[YOUSO10];
  static char operate[800];
  char namemode[YOUSO10];
  char *dir;

  printf("\n*******************************************************\n"); 
  printf("*******************************************************\n"); 
  printf(" Welcome to OpenMX  Ver. %s                            \n",Version_OpenMX); 
  printf(" Copyright (C), 2002-2007, T.Ozaki                     \n"); 
  printf(" OpenMX comes with ABSOLUTELY NO WARRANTY.             \n"); 
  printf(" This is free software, and you are welcome to         \n"); 
  printf(" redistribute it under the constitution of the GNU-GPL.\n");
  printf("*******************************************************\n"); 
  printf("*******************************************************\n\n\n"); 

  if (strcasecmp(mode,"S")==0){  
    dir = "input_example";
    sprintf(namemode,"runtest");
  }
  else if (strcasecmp(mode,"L")==0){  
    dir = "large_example";
    sprintf(namemode,"runtestL");
  }
  else if (strcasecmp(mode,"G")==0){  
    dir = "geoopt_example";
    sprintf(namemode,"runtestG");
  }

  /* print std */

  printf("\n");
  printf(" OpenMX is now in the mode in making of *.out files\n");
  printf(" which will be used in the test mode '%s'.\n",namemode);
  printf("\n");

  sprintf(operate,"ls %s/*.dat > ls_dat000000",dir);
  system(operate);

  sprintf(operate,"wc ls_dat000000 > ls_dat000001");
  system(operate);

  sprintf(fname0,"ls_dat000000");
  sprintf(fname1,"ls_dat000001");

  if ( ((fp0 = fopen(fname0,"r")) != NULL) &&
       ((fp1 = fopen(fname1,"r")) != NULL) )
    {

    fscanf(fp1,"%i",&Num_DatFiles);

    printf(" %2d dat files are found in the directory '%s'.\n\n\n",Num_DatFiles,dir);

    for (i=0; i<Num_DatFiles; i++){
      fscanf(fp0,"%s",fname_dat);

      if      (argc==2)  sprintf(operate,"./openmx %s",fname_dat);
      else if (argc==3)  sprintf(operate,"%s %s",argv[2],fname_dat);

      system(operate);

      input_open(fname_dat);

      input_string("System.Name",fname_dat2,"default");

      sprintf(operate,"cp %s.out %s/",fname_dat2,dir);
      system(operate);
       
      input_close();
    }

    fclose(fp0);
    fclose(fp1);
  }
  else{
    printf("Could not find ls_dat000000 or ls_dat000001\n");
    exit(1);
  }

  sprintf(operate,"rm ls_dat000000");
  system(operate);

  sprintf(operate,"rm ls_dat000001");
  system(operate);

  printf("\n\n\n\n");
  printf("Generated out files are stored in the directory '%s'.\n\n\n",dir);
}






