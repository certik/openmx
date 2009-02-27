/**********************************************************************
  Merge_LogFile.c:

     Merge_LogFile.c is a subrutine to merge several log files, and to
     write the merged data to file.out.

  Log of Merge_LogFile.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*  end stat section */
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

void Merge_LogFile(char *file)
{
  int i,NumFile,Gc_AN,c,n1,k;
  char operate[800];
  char rm_operate[YOUSO10];
  char fname[YOUSO10];
  char fname1[YOUSO10];
  char fname2[YOUSO10];
  char fexn[50][YOUSO10];
  FILE *fp,*fp1,*fp2;
  int numprocs,myid;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){

    NumFile = 15;

    sprintf(fexn[0],"UCell");
    sprintf(fexn[1],"TRN");
    sprintf(fexn[2],"rcn");
    sprintf(fexn[3],"DFTSCF");
    sprintf(fexn[4],"EV");
    sprintf(fexn[5],"OrbOpt");
    sprintf(fexn[6],"SD");
    sprintf(fexn[7],"MC");
    sprintf(fexn[8],"OM");
    sprintf(fexn[9],"VC");
    sprintf(fexn[10],"DM_onsite");
    sprintf(fexn[11],"dpm");
    sprintf(fexn[12],"crd");
    sprintf(fexn[13],"frac");
    sprintf(fexn[14],"CompTime");

    sprintf(fname,"%s%s.out",filepath,filename);
    fp = fopen(fname, "r");   

    if (fp!=NULL){
      fclose(fp); 
      remove(fname);
    }

    /* input file */    

    sprintf(fname1,"%s%s.out",filepath,filename);
    fp1 = fopen(fname1,"a");
    fseek(fp1,0,SEEK_END);

    fp2 = fopen(file,"r");
    if (fp2!=NULL){
      for (c=getc(fp2); c!=EOF; c=getc(fp2))  putc(c,fp1); 
      fclose(fp2); 
    }
    fclose(fp1); 

    /* merge log files */    

    for (i=0; i<NumFile; i++){
      sprintf(fname1,"%s%s.out",filepath,filename);
      fp1 = fopen(fname1,"a");
      fseek(fp1,0,SEEK_END);

      sprintf(fname2,"%s%s.%s",filepath,filename,fexn[i]);
      fp2 = fopen(fname2,"r");

      if (fp2!=NULL){
        for (c=getc(fp2); c!=EOF; c=getc(fp2))  putc(c,fp1); 
	fclose(fp2); 
      }
      fclose(fp1); 
    }   

    for (i=0; i<NumFile; i++){
      sprintf(fname,"%s%s.%s",filepath,filename,fexn[i]);
      fp = fopen(fname, "r");   

      if (fp!=NULL){
	fclose(fp); 
	sprintf(rm_operate,"%s",fname);
	remove(rm_operate);
      }
    }

    /*********************************************
                    delete *.ccs*                
    *********************************************/

    if (Cnt_switch==1){
      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
        sprintf(operate,"%s%s.ccs%i",filepath,filename,Gc_AN);
        remove(operate);
      }
    }

    /*********************************************
                      tar *.rst*
    *********************************************/

    /*
    sprintf(operate,"tar -cf %s%s.RST.tar -C %s %s_rst",
            filepath,filename,filepath,filename);
    system(operate);

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
      sprintf(operate,"%s%s_rst/%s.rst%i",filepath,filename,filename,Gc_AN);
      remove(operate);
    }

    for (k=0; k<=SpinP_switch; k++){
      for (n1=0; n1<Ngrid1; n1++){
        sprintf(operate,"%s%s_rst/%s.crst%i_%i",filepath,filename,filename,k,n1);
        remove(operate);
      }
    }

    sprintf(operate,"%s%s_rst",filepath,filename);
    rmdir(operate); 
    */

  }

}
