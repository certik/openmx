/**********************************************************************
  OutData.c:

     OutData.c is a subrutine to output values of electron densities,
     potetianls, wave functions on the grids in the format of
     Gaussian cube, and atomic cartesian coordinates.

  Log of OutData.c:

     12/May/2003  Released by T.Ozaki
     21/Feb/2006  xsf for non-collinear by F.Ishii

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#define CUBE_EXTENSION ".cube"

static void out_density();
static void out_Veff();
static void out_Vhart();
static void out_Vna();
static void out_Vxc();
static void out_grid();
static void out_atomxyz();
static void out_atomxsf();
static void out_coordinates_bulk();
static void out_Cluster_MO();
static void out_Cluster_NC_MO();
static void out_Bulk_MO();
static void out_OrbOpt(char *inputfile);

static void Print_CubeTitle(FILE *fp);
static void Print_CubeData(FILE *fp, char fext[], double *data, double *data1,char *op);
static void Print_CubeCData_MO(FILE *fp,dcomplex *data,char *op);
static void Print_CubeData_MO(FILE *fp, double *data, double *data1,char *op);
static void Print_VectorData(int flagout, FILE *fp, char fext[],
                             double *data0, double *data1,
                             double *data2, double *data3);

void OutData(char *inputfile)
{
  char operate[YOUSO10];
  int i,c;
  int numprocs,myid;
  char fname1[300];
  char fname2[300];
  FILE *fp1,*fp2;
  char buf[fp_bsize];          /* setvbuf */
  char buf1[fp_bsize];         /* setvbuf */
  char buf2[fp_bsize];         /* setvbuf */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("\noutputting data on grids to files...\n\n");

  out_atomxyz();

  /*
  out_atomxsf();
  */

  if (1<=level_fileout){
    out_density(); 
    out_Veff();
    out_Vhart(); 
  }

  if (2<=level_fileout){
    out_grid();
    if (ProExpn_VNA==0) out_Vna();
    out_Vxc();
  }

  if (specified_system!=0) out_coordinates_bulk();

  if (MO_fileout==1 && specified_system==0){

    /* spin non-collinear */
    if (SpinP_switch==3)
      out_Cluster_NC_MO();
    /* spin collinear */
    else
      out_Cluster_MO();
  }
  else if (MO_fileout==1 && specified_system!=0){
    out_Bulk_MO();
  }

  if (Cnt_switch==1){
    out_OrbOpt(inputfile);    
  } 

  /* divide-conquer, gdc, or Krylov */

  if (Dos_fileout && (Solver==5 || Solver==6 || Solver==8) ) {  

    if (myid==Host_ID){

#ifdef xt3

      sprintf(fname1,"%s%s.Dos.vec",filepath,filename);
      fp1 = fopen(fname1,"a");

      if (fp1!=NULL){
        remove(fname1); 
        fclose(fp1); 
      }
  
      for (i=0; i<numprocs; i++){ 

        sprintf(fname1,"%s%s.Dos.vec",filepath,filename);
        fp1 = fopen(fname1,"a");
        fseek(fp1,0,SEEK_END);
        sprintf(fname2,"%s%s.Dos.vec%i",filepath,filename,i);
        fp2 = fopen(fname2,"r");

        if (fp2!=NULL){
          for (c=getc(fp2); c!=EOF; c=getc(fp2))  putc(c,fp1); 
	  fclose(fp2); 
        }
        fclose(fp1); 
      }

      for (i=0; i<numprocs; i++){
        sprintf(fname2,"%s%s.Dos.vec%i",filepath,filename,i);
        remove(fname2);
      }

#else

      for (i=0; i<numprocs; i++){
        if (i==0) 
          sprintf(operate,"cat %s%s.Dos.vec%i >  tmp1",filepath,filename,i);
        else
          sprintf(operate,"cat %s%s.Dos.vec%i >> tmp1",filepath,filename,i);

        system(operate);
        sprintf(operate,"rm %s%s.Dos.vec%i",filepath,filename,i);
        system(operate);
      }

      sprintf(operate,"mv tmp1 %s%s.Dos.vec",filepath,filename);
      system(operate);

#endif

    }
  }

}


void out_density()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  int MN;
  double x,y,z,vx,vy,vz;
  double phi,theta,sden,oden;
  double xmin,ymin,zmin,xmax,ymax,zmax;
  char fname[YOUSO10];
  char file1[YOUSO10] = ".tden";
  char file2[YOUSO10] = ".den0";
  char file3[YOUSO10] = ".den1";
  char file4[YOUSO10] = ".sden";
  char file5[YOUSO10] = ".ncsden.txt";
  char file6[YOUSO10] = ".nc.txt";
  char file7[YOUSO10] = ".ncoden.txt";
  char file8[YOUSO10] = ".nco.txt";
  char file9[YOUSO10] = ".nc.xsf";
  char file10[YOUSO10] = ".ncsden.xsf";
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID;
  double scaxsf;
  int flagout;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  strcat(file1,CUBE_EXTENSION);
  strcat(file2,CUBE_EXTENSION);
  strcat(file3,CUBE_EXTENSION);
  strcat(file4,CUBE_EXTENSION);

  /****************************************************
                  total electron density
  ****************************************************/

  sprintf(fname,"%s%s%s%i",filepath,filename,file1,myid);

  if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (myid==Host_ID) Print_CubeTitle(fp);

    if (SpinP_switch==0) {
      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[0][MN] = 2.0*Density_Grid[0][MN];
      }
      Print_CubeData(fp,file1,Density_Grid[0],(void*)NULL,(void*)NULL);
    }
    else {
      Print_CubeData(fp,file1,Density_Grid[0],Density_Grid[1],"add");
    }

  }
  else{
    printf("Failure of saving the electron density\n");
  }

  /* spin polization */

  if (SpinP_switch==1 || SpinP_switch==3){

    /****************************************************
                  up-spin electron density
    ****************************************************/

    sprintf(fname,"%s%s%s%i",filepath,filename,file2,myid);

    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (myid==Host_ID) Print_CubeTitle(fp);
      Print_CubeData(fp,file2,Density_Grid[0],NULL,NULL);
    }
    else{
      printf("Failure of saving the electron density\n");
    }

    /****************************************************
                  down-spin electron density
    ****************************************************/

    sprintf(fname,"%s%s%s%i",filepath,filename,file3,myid);

    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (myid==Host_ID) Print_CubeTitle(fp);
      Print_CubeData(fp,file3,Density_Grid[1],NULL,NULL);
    }
    else{
      printf("Failure of saving the electron density\n");
    }

    /****************************************************
                    spin electron density
    ****************************************************/

    sprintf(fname,"%s%s%s%i",filepath,filename,file4,myid);

    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (myid==Host_ID) Print_CubeTitle(fp);
      Print_CubeData(fp,file4,Density_Grid[0],Density_Grid[1],"diff");
    }
    else{
      printf("Failure of saving the electron density\n");
    }
  }

  /****************************************************
        spin electron density with a spin vector
  ****************************************************/

  if (SpinP_switch==3){

     /*for gOpenMol */

    sprintf(fname,"%s%s%s%i",filepath,filename,file5,myid);
    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      flagout=0;

      Print_VectorData(flagout,fp,file5,Density_Grid[0],Density_Grid[1],
                                Density_Grid[2],Density_Grid[3]);
    }
    else{
      printf("Failure of saving the electron spin density with a vector\n");
    }

    /*for XCrysDen */

    /*
    sprintf(fname,"%s%s%s%i",filepath,filename,file10,myid);
    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  
#endif
 
      flagout=1;

      Print_VectorData(flagout,fp,file10,Density_Grid[0],Density_Grid[1],
                                Density_Grid[2],Density_Grid[3]);
    }
    else{
      printf("Failure of saving the electron spin density with a vector\n");
    }
    */

    /* non-collinear spin density by Mulliken population */

    if (myid==Host_ID){

      /* for gOpenMol */

      sprintf(fname,"%s%s%s",filepath,filename,file6);
      if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
        setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        xmin = 10.0e+10;
        ymin = 10.0e+10;
        zmin = 10.0e+10;
        xmax = -10.0e+10;
        ymax = -10.0e+10;
        zmax = -10.0e+10;

        for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
          if (Gxyz[ct_AN][1]<xmin) xmin = Gxyz[ct_AN][1];
          if (Gxyz[ct_AN][2]<ymin) ymin = Gxyz[ct_AN][2];
          if (Gxyz[ct_AN][3]<zmin) zmin = Gxyz[ct_AN][3];
          if (xmax<Gxyz[ct_AN][1]) xmax = Gxyz[ct_AN][1];
          if (ymax<Gxyz[ct_AN][2]) ymax = Gxyz[ct_AN][2];
          if (zmax<Gxyz[ct_AN][3]) zmax = Gxyz[ct_AN][3];
        }

        fprintf(fp,"3 200\n");
        fprintf(fp,"%4d %4d %4d\n",atomnum,0,0);
        fprintf(fp,"%13.3E %13.3E\n",BohrR*xmin,BohrR*xmax);
        fprintf(fp,"%13.3E %13.3E\n",BohrR*ymin,BohrR*ymax);
        fprintf(fp,"%13.3E %13.3E\n",BohrR*zmin,BohrR*zmax);

        for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
          theta = Angle0_Spin[ct_AN];
          phi   = Angle1_Spin[ct_AN];

          sden = InitN_USpin[ct_AN] - InitN_DSpin[ct_AN];

          vx = sden*sin(theta)*cos(phi);
          vy = sden*sin(theta)*sin(phi);
          vz = sden*cos(theta);

          fprintf(fp,"%13.3E %13.3E %13.3E %13.3E %13.3E %13.3E\n",
                  BohrR*Gxyz[ct_AN][1],
                  BohrR*Gxyz[ct_AN][2],
                  BohrR*Gxyz[ct_AN][3],
                  vx,vy,vz);
 	}

        fclose(fp);
      }
      else{
        printf("Failure of saving the mulliken spin vector\n");
      }

      /* for XCrysDen  .xsf*/

      /*
      sprintf(fname,"%s%s%s",filepath,filename,file9);
      if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
        setvbuf(fp,buf,_IOFBF,fp_bsize); 
#endif

        xmin = 10.0e+10;
        ymin = 10.0e+10;
        zmin = 10.0e+10;
        xmax = -10.0e+10;
        ymax = -10.0e+10;
        zmax = -10.0e+10;
        scaxsf=0.02;

        for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
          if (Gxyz[ct_AN][1]<xmin) xmin = Gxyz[ct_AN][1];
          if (Gxyz[ct_AN][2]<ymin) ymin = Gxyz[ct_AN][2];
          if (Gxyz[ct_AN][3]<zmin) zmin = Gxyz[ct_AN][3];
          if (xmax<Gxyz[ct_AN][1]) xmax = Gxyz[ct_AN][1];
          if (ymax<Gxyz[ct_AN][2]) ymax = Gxyz[ct_AN][2];
          if (zmax<Gxyz[ct_AN][3]) zmax = Gxyz[ct_AN][3];
        }

        fprintf(fp,"CRYSTAL\n");
        fprintf(fp,"PRIMVEC\n");
        fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[1][1], BohrR*tv[1][2], BohrR*tv[1][3]);
        fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[2][1], BohrR*tv[2][2], BohrR*tv[2][3]);
        fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[3][1], BohrR*tv[3][2], BohrR*tv[3][3]);
        fprintf(fp,"PRIMCOORD 1\n");
        fprintf(fp,"%4d %4d\n",atomnum, 1);

        for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
            i = WhatSpecies[ct_AN];

          theta = Angle0_Spin[ct_AN];
          phi   = Angle1_Spin[ct_AN];
          sden = InitN_USpin[ct_AN] - InitN_DSpin[ct_AN];

          vx = sden*sin(theta)*cos(phi);
          vy = sden*sin(theta)*sin(phi);
          vz = sden*cos(theta);

          fprintf(fp,"%4d %13.3E %13.3E %13.3E %13.3E %13.3E %13.3E\n",
	           Spe_WhatAtom[i],               
                  BohrR*Gxyz[ct_AN][1],
                  BohrR*Gxyz[ct_AN][2],
                  BohrR*Gxyz[ct_AN][3],
                  scaxsf*vx,scaxsf*vy,scaxsf*vz);
 	}

        fclose(fp);
      }
      else{
        printf("Failure of saving the mulliken spin vector\n");
      }
      */

      /* non-collinear obital magnetic moment by 'on-site' approximation */

      sprintf(fname,"%s%s%s",filepath,filename,file8);
      if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
        setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        xmin =  10.0e+10;
        ymin =  10.0e+10;
        zmin =  10.0e+10;
        xmax = -10.0e+10;
        ymax = -10.0e+10;
        zmax = -10.0e+10;

        for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
          if (Gxyz[ct_AN][1]<xmin) xmin = Gxyz[ct_AN][1];
          if (Gxyz[ct_AN][2]<ymin) ymin = Gxyz[ct_AN][2];
          if (Gxyz[ct_AN][3]<zmin) zmin = Gxyz[ct_AN][3];
          if (xmax<Gxyz[ct_AN][1]) xmax = Gxyz[ct_AN][1];
          if (ymax<Gxyz[ct_AN][2]) ymax = Gxyz[ct_AN][2];
          if (zmax<Gxyz[ct_AN][3]) zmax = Gxyz[ct_AN][3];
        }

        fprintf(fp,"3 200\n");
        fprintf(fp,"%4d %4d %4d\n",atomnum,0,0);
        fprintf(fp,"%13.3E %13.3E\n",BohrR*xmin,BohrR*xmax);
        fprintf(fp,"%13.3E %13.3E\n",BohrR*ymin,BohrR*ymax);
        fprintf(fp,"%13.3E %13.3E\n",BohrR*zmin,BohrR*zmax);

        for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

          theta = Angle0_Orbital[ct_AN];
          phi   = Angle1_Orbital[ct_AN];
          oden  = OrbitalMoment[ct_AN];

          vx = oden*sin(theta)*cos(phi);
          vy = oden*sin(theta)*sin(phi);
          vz = oden*cos(theta);

          fprintf(fp,"%13.3E %13.3E %13.3E %13.3E %13.3E %13.3E\n",
                  BohrR*Gxyz[ct_AN][1],
                  BohrR*Gxyz[ct_AN][2],
                  BohrR*Gxyz[ct_AN][3],
                  vx,vy,vz);
 	}

        fclose(fp);
      }
      else{
        printf("Failure of saving orbital magnetic moments\n");
      }
  
    } /* if (myid==Host_ID) */
  } /* if (SpinP_switch==3) */

}



static void out_Vhart()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char fname[YOUSO10];
  char file1[YOUSO10] = ".vhart";
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  strcat(file1,CUBE_EXTENSION);

  /****************************************************
                     Hartree potential
  ****************************************************/

  sprintf(fname,"%s%s%s%i",filepath,filename,file1,myid);
  if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (myid==Host_ID) Print_CubeTitle(fp);
    Print_CubeData(fp,file1,dVHart_Grid,NULL,NULL);
  }
  else{
    printf("Failure of saving the electron density\n");
  }

}




static void out_Vna()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char fname[YOUSO10];
  char file1[YOUSO10] = ".vna";
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  strcat(file1,CUBE_EXTENSION);

  /****************************************************
                   neutral atom potential
  ****************************************************/

  sprintf(fname,"%s%s%s%i",filepath,filename,file1,myid);

  if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (myid==Host_ID) Print_CubeTitle(fp);
    Print_CubeData(fp,file1,VNA_Grid,NULL,NULL);
  }
  else{
    printf("Failure of saving the electron density\n");
  }
}





static void out_Vxc()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char fname[YOUSO10];
  char file1[YOUSO10] = ".vxc0";
  char file2[YOUSO10] = ".vxc1";
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  strcat(file1,CUBE_EXTENSION);
  strcat(file2,CUBE_EXTENSION);

  /****************************************************
       exchange-correlation potential for up-spin
  ****************************************************/

  sprintf(fname,"%s%s%s%i",filepath,filename,file1,myid);
  if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (myid==Host_ID) Print_CubeTitle(fp);
    Print_CubeData(fp,file1,Vxc_Grid[0],NULL,NULL);
  }
  else{
    printf("Failure of saving the electron density\n");
  }

  /****************************************************
     exchange-correlation potential for down-spin
  ****************************************************/

  if (SpinP_switch==1){

  sprintf(fname,"%s%s%s%i",filepath,filename,file2,myid);

    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (myid==Host_ID) Print_CubeTitle(fp);
      Print_CubeData(fp,file2,Vxc_Grid[1],NULL,NULL);
    }
    else{
      printf("Failure of saving the electron density\n");
    }
  }

}





static void out_Veff()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char fname[YOUSO10];
  char file1[YOUSO10] = ".v0";
  char file2[YOUSO10] = ".v1";
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  strcat(file1,CUBE_EXTENSION);
  strcat(file2,CUBE_EXTENSION);

  /****************************************************
           Kohn-Sham potential for up-spin
  ****************************************************/

  sprintf(fname,"%s%s%s%i",filepath,filename,file1,myid);
  if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (myid==Host_ID) Print_CubeTitle(fp);
    Print_CubeData(fp,file1,Vpot_Grid[0],NULL,NULL);
  }
  else{
    printf("Failure of saving the electron density\n");
  }

  /****************************************************
           Kohn-Sham potential for down-spin
  ****************************************************/

  if (SpinP_switch==1){

    sprintf(fname,"%s%s%s%i",filepath,filename,file2,myid);
    if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (myid==Host_ID) Print_CubeTitle(fp);
      Print_CubeData(fp,file2,Vpot_Grid[1],NULL,NULL);
    }
    else{
      printf("Failure of saving the electron density\n");
    }
  }

}




static void out_grid()
{
  int N;
  char file1[YOUSO10] = ".grid";
  int numprocs,myid,ID;
  double x,y,z;
  double Cxyz[4];
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
     output the real space grids to a file, *.grid
  ****************************************************/

  if (myid==Host_ID){

    fnjoint(filepath,filename,file1);

    if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      for (N=0; N<TNumGrid; N++){
        Get_Grid_XYZ(N,Cxyz);
        x = Cxyz[1];
        y = Cxyz[2];
        z = Cxyz[3];
        fprintf(fp,"%5d  %19.12f %19.12f %19.12f\n", N,BohrR*x,BohrR*y,BohrR*z);
      }
      fclose(fp);
    }
    else{
      printf("Failure of saving grids\n");
    }
  }
}







static void out_atomxyz()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char filexyz[YOUSO10] = ".xyz";
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                cartesian coordinates
  ****************************************************/

  if (myid==Host_ID){

    fnjoint(filepath,filename,filexyz);
    if ((fp = fopen(filexyz,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp,"%i \n\n",atomnum);
      for (k=1; k<=atomnum; k++){
        i = WhatSpecies[k];
        j = Spe_WhatAtom[i];
        fprintf(fp,"%s   %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",
	        Atom_Symbol[j],                
	        Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	        Gxyz[k][24],Gxyz[k][25],Gxyz[k][26]);
      }
      fclose(fp);
    }
    else{
      printf("could not save the xyz file\n");
    }
  }
}



 static void out_atomxsf()
{
  static int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  static char filexsf[YOUSO10] = ".coord.xsf";
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  static int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                cartesian coordinates
  ****************************************************/

  if (myid==Host_ID){

    fnjoint(filepath,filename,filexsf);
    if ((fp = fopen(filexsf,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp,"CRYSTAL\n");
        fprintf(fp,"PRIMVEC\n");
        fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[1][1], BohrR*tv[1][2], BohrR*tv[1][3]);
        fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[2][1], BohrR*tv[2][2], BohrR*tv[2][3]);
        fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[3][1], BohrR*tv[3][2], BohrR*tv[3][3]);
        fprintf(fp,"PRIMCOORD 1\n");
        fprintf(fp,"%4d %d\n",atomnum, 1);

      for (k=1; k<=atomnum; k++){
        i = WhatSpecies[k];
        /*j = Spe_WhatAtom[i];*/
        fprintf(fp,"%4d  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",
	        Spe_WhatAtom[i],               
	        Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	        Gxyz[k][24],Gxyz[k][25],Gxyz[k][26]);
      }
      fclose(fp);
    }
    else{
      printf("failure of saving coord.xsf file\n");
    }
  }
}



void out_coordinates_bulk()
{
  int n,i1,i2,i3,ct_AN,i,j;
  double tx,ty,tz,x,y,z;
  char file1[YOUSO10] = ".bulk.xyz";
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
        atomic coordinates including copied cells
  ****************************************************/

  if (myid==Host_ID){

    n = 1;

    fnjoint(filepath,filename,file1);

    if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp,"%d\n\n",atomnum*(2*n+1)*(2*n+1)*(2*n+1));

      for (i1=-n; i1<=n; i1++){
        for (i2=-n; i2<=n; i2++){
          for (i3=-n; i3<=n; i3++){

            tx = (double)i1*tv[1][1] + (double)i2*tv[2][1] + (double)i3*tv[3][1];
            ty = (double)i1*tv[1][2] + (double)i2*tv[2][2] + (double)i3*tv[3][2];
            tz = (double)i1*tv[1][3] + (double)i2*tv[2][3] + (double)i3*tv[3][3];

            for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
              i = WhatSpecies[ct_AN];
              j = Spe_WhatAtom[i];

              x = BohrR*(Gxyz[ct_AN][1] + tx); 
              y = BohrR*(Gxyz[ct_AN][2] + ty); 
              z = BohrR*(Gxyz[ct_AN][3] + tz); 
              fprintf(fp,"%s %8.5f %8.5f %8.5f\n",Atom_Symbol[j],x,y,z);
            } 
          }
        }
      }

      fclose(fp);
    }
    else{
      printf("Failure of saving atomic coordinates\n");
    }
  }

}
  



void out_Cluster_MO()
{ 
  int Mc_AN,Gc_AN,Cwan,NO0,spin,Nc;
  int orbit,GN,spe,i,i1,i2,i3,so;
  double *MO_Grid_tmp;
  double *MO_Grid;
  char file1[YOUSO10];
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    allocation of arrays:

    double MO_Grid_tmp[TNumGrid];
    double MO_Grid[TNumGrid];
  ****************************************************/

  MO_Grid_tmp = (double*)malloc(sizeof(double)*TNumGrid);
  MO_Grid     = (double*)malloc(sizeof(double)*TNumGrid);

  /*************
      HOMOs
  *************/

  for (so=0; so<=SO_switch; so++){
    for (spin=0; spin<=SpinP_switch; spin++){
      for (orbit=0; orbit<num_HOMOs; orbit++){ 

	/* calc. MO on grids */

	for (GN=0; GN<TNumGrid; GN++) MO_Grid_tmp[GN] = 0.0;

	for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	  Gc_AN = M2G[Mc_AN];    
	  Cwan = WhatSpecies[Gc_AN];
	  NO0 = Spe_Total_CNO[Cwan];

          if (so==0){
  	    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	      GN = GridListAtom[Mc_AN][Nc];
	      for (i=0; i<NO0; i++){
	        MO_Grid_tmp[GN] += HOMOs_Coef[0][spin][orbit][Gc_AN][i].r*
                                                Orbs_Grid[Mc_AN][i][Nc];
	      }
	    }  
	  }

          else if (so==1){
  	    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	      GN = GridListAtom[Mc_AN][Nc];
	      for (i=0; i<NO0; i++){
	        MO_Grid_tmp[GN] += HOMOs_Coef[0][spin][orbit][Gc_AN][i].i*
                                                 Orbs_Grid[Mc_AN][i][Nc];
	      }
	    }  
	  }

	}

	MPI_Reduce(&MO_Grid_tmp[0], &MO_Grid[0], TNumGrid, MPI_DOUBLE,
		   MPI_SUM, Host_ID, mpi_comm_level1);

	/* output HOMO on grids */

	if (myid==Host_ID){ 

          if (so==0)
  	    sprintf(file1,"%s%s.homo%i_%i_r%s",filepath,filename,spin,orbit,CUBE_EXTENSION);
          else if (so==1)
  	    sprintf(file1,"%s%s.homo%i_%i_i%s",filepath,filename,spin,orbit,CUBE_EXTENSION);

	  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
            setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	    Print_CubeTitle(fp);
	    Print_CubeData_MO(fp,MO_Grid,NULL,NULL);
	    fclose(fp);
	  }
	  else{
	    printf("Failure of saving MOs\n");
	  }

	}


	/*
	{

	  double sumr,sumi,kx,ky,kz,x,y,z,tmp,co,si;
          double Cxyz[4];

	  sumr = 0.0;
	  sumi = 0.0;

	  for (GN=0; GN<TNumGrid; GN++){

	    Get_Grid_XYZ(GN,Cxyz);

	    kx = rtv[1][1];
	    ky = rtv[1][2];
	    kz = rtv[1][3];

	    x = Cxyz[1];
	    y = Cxyz[2];
	    z = Cxyz[3];

	    tmp = -(kx*x + ky*y + kz*z);
	    co = cos(tmp);  
	    si = sin(tmp);           

	    sumr += MO_Grid[GN]*MO_Grid[GN]*co;
	    sumi += MO_Grid[GN]*MO_Grid[GN]*si;
	  }

          sumr *= (Cell_Volume/Ngrid1/Ngrid2/Ngrid3);
          sumi *= (Cell_Volume/Ngrid1/Ngrid2/Ngrid3);

          printf("sumr=%15.12f sumi=%15.12f\n",sumr,sumi);

        }
	*/        


      }  /* orbit */ 
    }  /* spin */ 
  }  /* so */

  /*************
      LUMOs
  *************/

  for (so=0; so<=SO_switch; so++){
    for (spin=0; spin<=SpinP_switch; spin++){
      for (orbit=0; orbit<num_LUMOs; orbit++){ 

	/* calc. MO on grids */

	for (GN=0; GN<TNumGrid; GN++) MO_Grid_tmp[GN] = 0.0;

	for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	  Gc_AN = M2G[Mc_AN];    
	  Cwan = WhatSpecies[Gc_AN];
	  NO0 = Spe_Total_CNO[Cwan];

          if (so==0){
	    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	      GN = GridListAtom[Mc_AN][Nc];
	      for (i=0; i<NO0; i++){
		MO_Grid_tmp[GN] += LUMOs_Coef[0][spin][orbit][Gc_AN][i].r*
		                                Orbs_Grid[Mc_AN][i][Nc];
	      }
	    } 
	  }

          else if (so==1){
	    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	      GN = GridListAtom[Mc_AN][Nc];
	      for (i=0; i<NO0; i++){
		MO_Grid_tmp[GN] += LUMOs_Coef[0][spin][orbit][Gc_AN][i].i*
		                                Orbs_Grid[Mc_AN][i][Nc];
	      }
	    } 
	  }
 
	}

	/* output LUMO on grids */

	MPI_Reduce(&MO_Grid_tmp[0], &MO_Grid[0], TNumGrid, MPI_DOUBLE,
		   MPI_SUM, Host_ID, mpi_comm_level1);

	if (myid==Host_ID){ 

          if (so==0)
  	    sprintf(file1,"%s%s.lumo%i_%i_r%s",filepath,filename,spin,orbit,CUBE_EXTENSION);
          else if (so==0)
  	    sprintf(file1,"%s%s.lumo%i_%i_i%s",filepath,filename,spin,orbit,CUBE_EXTENSION);


	  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
            setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	    Print_CubeTitle(fp);
	    Print_CubeData_MO(fp,MO_Grid,NULL,NULL);
	    fclose(fp);
	  }
	  else{
	    printf("Failure of saving MOs\n");
	  }
	}

      }  /* orbit */ 
    }  /* spin */ 
  }  /* so */

  /****************************************************
    freeing of arrays:

    double MO_Grid[TNumGrid];
    double MO_Grid_tmp[TNumGrid];
  ****************************************************/

  free(MO_Grid);
  free(MO_Grid_tmp);
} 


void out_Cluster_NC_MO()
{ 
  int Mc_AN,Gc_AN,Cwan,NO0,spin,Nc;
  int orbit,GN,spe,i,i1,i2,i3;
  dcomplex *MO_Grid;
  double *RMO_Grid_tmp,*IMO_Grid_tmp;
  double *RMO_Grid,*IMO_Grid;
  char file1[YOUSO10];
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    allocation of arrays:

    dcomplex MO_Grid[TNumGrid];
    double RMO_Grid_tmp[TNumGrid];
    double IMO_Grid_tmp[TNumGrid];
    double RMO_Grid[TNumGrid];
    double IMO_Grid[TNumGrid];
  ****************************************************/

  MO_Grid = (dcomplex*)malloc(sizeof(dcomplex)*TNumGrid);
  RMO_Grid_tmp = (double*)malloc(sizeof(double)*TNumGrid);
  IMO_Grid_tmp = (double*)malloc(sizeof(double)*TNumGrid);
  RMO_Grid     = (double*)malloc(sizeof(double)*TNumGrid);
  IMO_Grid     = (double*)malloc(sizeof(double)*TNumGrid);

  /*************
      HOMOs
  *************/

  for (spin=0; spin<=1; spin++){
    for (orbit=0; orbit<num_HOMOs; orbit++){ 

      /* calc. MO on grids */

      for (GN=0; GN<TNumGrid; GN++){
	RMO_Grid_tmp[GN] = 0.0;
	IMO_Grid_tmp[GN] = 0.0;
      }

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];    
	Cwan = WhatSpecies[Gc_AN];
	NO0 = Spe_Total_CNO[Cwan];
	for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	  GN = GridListAtom[Mc_AN][Nc];
	  for (i=0; i<NO0; i++){
 	    RMO_Grid_tmp[GN] += HOMOs_Coef[0][spin][orbit][Gc_AN][i].r*
  	                                     Orbs_Grid[Mc_AN][i][Nc];
	    IMO_Grid_tmp[GN] += HOMOs_Coef[0][spin][orbit][Gc_AN][i].i*
                                             Orbs_Grid[Mc_AN][i][Nc];
	  }
	}  
      }

      MPI_Reduce(&RMO_Grid_tmp[0], &RMO_Grid[0], TNumGrid, MPI_DOUBLE,
		 MPI_SUM, Host_ID, mpi_comm_level1);
      MPI_Reduce(&IMO_Grid_tmp[0], &IMO_Grid[0], TNumGrid, MPI_DOUBLE,
		 MPI_SUM, Host_ID, mpi_comm_level1);

      for (GN=0; GN<TNumGrid; GN++){
	MO_Grid[GN].r = RMO_Grid[GN];
	MO_Grid[GN].i = IMO_Grid[GN];
      }


      if (myid==Host_ID){ 

	/* output the real part of HOMOs on grids */

	sprintf(file1,"%s%s.homo%i_%i_r%s",
		filepath,filename,spin,orbit,CUBE_EXTENSION);

	if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
          setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	  Print_CubeTitle(fp);
	  Print_CubeCData_MO(fp,MO_Grid,"r");
	  fclose(fp);
	}
	else{
	  printf("Failure of saving MOs\n");
	}

	/* output the imaginary part of HOMOs on grids */

	sprintf(file1,"%s%s.homo%i_%i_i%s", 
		filepath,filename,spin,orbit,CUBE_EXTENSION);

	if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
          setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	  Print_CubeTitle(fp);
	  Print_CubeCData_MO(fp,MO_Grid,"i");

	  fclose(fp);
	}
	else{
	  printf("Failure of saving MOs\n");
	}

      }
    } /* orbit */ 
  } /* spin  */

  /*************
      LUMOs
  *************/

  for (spin=0; spin<=1; spin++){
    for (orbit=0; orbit<num_HOMOs; orbit++){ 

      /* calc. MO on grids */

      for (GN=0; GN<TNumGrid; GN++){
	RMO_Grid_tmp[GN] = 0.0;
	IMO_Grid_tmp[GN] = 0.0;
      }

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];    
	Cwan = WhatSpecies[Gc_AN];
	NO0 = Spe_Total_CNO[Cwan];
	for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	  GN = GridListAtom[Mc_AN][Nc];
	  for (i=0; i<NO0; i++){
	    RMO_Grid_tmp[GN] += LUMOs_Coef[0][spin][orbit][Gc_AN][i].r*
                                             Orbs_Grid[Mc_AN][i][Nc];
	    IMO_Grid_tmp[GN] += LUMOs_Coef[0][spin][orbit][Gc_AN][i].i*
                                             Orbs_Grid[Mc_AN][i][Nc];
	  }
	}  
      }

      MPI_Reduce(&RMO_Grid_tmp[0], &RMO_Grid[0], TNumGrid, MPI_DOUBLE,
		 MPI_SUM, Host_ID, mpi_comm_level1);
      MPI_Reduce(&IMO_Grid_tmp[0], &IMO_Grid[0], TNumGrid, MPI_DOUBLE,
		 MPI_SUM, Host_ID, mpi_comm_level1);

      for (GN=0; GN<TNumGrid; GN++){
	MO_Grid[GN].r = RMO_Grid[GN];
	MO_Grid[GN].i = IMO_Grid[GN];
      }


      if (myid==Host_ID){ 

	/* output the real part of LUMOs on grids */

	sprintf(file1,"%s%s.lumo%i_%i_r%s",
		filepath,filename,spin,orbit,CUBE_EXTENSION);

	if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
          setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	  Print_CubeTitle(fp);
	  Print_CubeCData_MO(fp,MO_Grid,"r");
	  fclose(fp);
	}
	else{
	  printf("Failure of saving MOs\n");
	}

	/* output the imaginary part of HOMOs on grids */

	sprintf(file1,"%s%s.lumo%i_%i_i%s", 
		filepath,filename,spin,orbit,CUBE_EXTENSION);

	if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
          setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	  Print_CubeTitle(fp);
	  Print_CubeCData_MO(fp,MO_Grid,"i");

	  fclose(fp);
	}
	else{
	  printf("Failure of saving MOs\n");
	}

      }
    }  /* orbit */ 
  }  /* spin  */

  /****************************************************
    freeing of arrays:

    dcomplex MO_Grid[TNumGrid];
    double RMO_Grid_tmp[TNumGrid];
    double IMO_Grid_tmp[TNumGrid];
    double RMO_Grid[TNumGrid];
    double IMO_Grid[TNumGrid];
  ****************************************************/

  free(MO_Grid);
  free(RMO_Grid_tmp);
  free(IMO_Grid_tmp);
  free(RMO_Grid);
  free(IMO_Grid);
} 




void out_Bulk_MO()
{ 
  int Mc_AN,Gc_AN,Cwan,NO0,spin,Nc;
  int kloop,Mh_AN,h_AN,Gh_AN,Rnh,Hwan;
  int NO1,l1,l2,l3,Nog,RnG,Nh,Rn;
  int orbit,GN,spe,i,j,i1,i2,i3,spinmax;
  double co,si,k1,k2,k3,kRn,ReCoef,ImCoef;
  double *RMO_Grid;
  double *IMO_Grid;
  double *RMO_Grid_tmp;
  double *IMO_Grid_tmp;
  dcomplex *MO_Grid;
  char file1[YOUSO10];
  char buf[fp_bsize];          /* setvbuf */
  FILE *fp;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    allocation of arrays:

    dcomplex MO_Grid[TNumGrid];
    double RMO_Grid[TNumGrid];
    double IMO_Grid[TNumGrid];
    double RMO_Grid_tmp[TNumGrid];
    double IMO_Grid_tmp[TNumGrid];
  ****************************************************/

  MO_Grid = (dcomplex*)malloc(sizeof(dcomplex)*TNumGrid);
  RMO_Grid = (double*)malloc(sizeof(double)*TNumGrid);
  IMO_Grid = (double*)malloc(sizeof(double)*TNumGrid);
  RMO_Grid_tmp = (double*)malloc(sizeof(double)*TNumGrid);
  IMO_Grid_tmp = (double*)malloc(sizeof(double)*TNumGrid);

  if      (SpinP_switch==0) spinmax = 0;
  else if (SpinP_switch==1) spinmax = 1;
  else if (SpinP_switch==3) spinmax = 1;

  /*************
       HOMOs
  *************/

  for (kloop=0; kloop<MO_Nkpoint; kloop++){

    k1 = MO_kpoint[kloop][1];
    k2 = MO_kpoint[kloop][2];
    k3 = MO_kpoint[kloop][3];

    for (spin=0; spin<=spinmax; spin++){
      for (orbit=0; orbit<Bulk_Num_HOMOs[kloop]; orbit++){ 

	/* calc. MO on grids */

	for (GN=0; GN<TNumGrid; GN++){
          RMO_Grid_tmp[GN] = 0.0;
          IMO_Grid_tmp[GN] = 0.0;
	}

        for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
          Gc_AN = M2G[Mc_AN];    
	  Cwan = WhatSpecies[Gc_AN];
	  NO0 = Spe_Total_CNO[Cwan];

          for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
            GN = GridListAtom[Mc_AN][Nc];
            Rn = CellListAtom[Mc_AN][Nc];

            l1 =-atv_ijk[Rn][1];
            l2 =-atv_ijk[Rn][2];
            l3 =-atv_ijk[Rn][3];

            kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;
            si = sin(2.0*PI*kRn);
            co = cos(2.0*PI*kRn);

	    for (i=0; i<NO0; i++){
              ReCoef = co*HOMOs_Coef[kloop][spin][orbit][Gc_AN][i].r 
                      -si*HOMOs_Coef[kloop][spin][orbit][Gc_AN][i].i; 
              ImCoef = co*HOMOs_Coef[kloop][spin][orbit][Gc_AN][i].i
                      +si*HOMOs_Coef[kloop][spin][orbit][Gc_AN][i].r;
              RMO_Grid_tmp[GN] += ReCoef*Orbs_Grid[Mc_AN][i][Nc];
              IMO_Grid_tmp[GN] += ImCoef*Orbs_Grid[Mc_AN][i][Nc];
	    }
          }
	}

        MPI_Reduce(&RMO_Grid_tmp[0], &RMO_Grid[0], TNumGrid, MPI_DOUBLE,
                    MPI_SUM, Host_ID, mpi_comm_level1);
        MPI_Reduce(&IMO_Grid_tmp[0], &IMO_Grid[0], TNumGrid, MPI_DOUBLE,
                    MPI_SUM, Host_ID, mpi_comm_level1);

	for (GN=0; GN<TNumGrid; GN++){
          MO_Grid[GN].r = RMO_Grid[GN];
          MO_Grid[GN].i = IMO_Grid[GN];
	}

        if (myid==Host_ID){ 

  	  /* output the real part of HOMOs on grids */

	  sprintf(file1,"%s%s.homo%i_%i_%i_r%s",
                  filepath,filename,kloop,spin,orbit,CUBE_EXTENSION);
  
  	  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
            setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

            Print_CubeTitle(fp);
            Print_CubeCData_MO(fp,MO_Grid,"r");
	    fclose(fp);
	  }
	  else{
	    printf("Failure of saving MOs\n");
	  }

  	  /* output the imaginary part of HOMOs on grids */

	  sprintf(file1,"%s%s.homo%i_%i_%i_i%s", 
                  filepath,filename,kloop,spin,orbit,CUBE_EXTENSION);

	  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
            setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

            Print_CubeTitle(fp);
            Print_CubeCData_MO(fp,MO_Grid,"i");

	    fclose(fp);
	  }
  	  else{
	    printf("Failure of saving MOs\n");
  	  }
	}

      } /* orbit */
    } /* spin */ 
  } /* kloop */

  /*************
       LUMOs
  *************/

  for (kloop=0; kloop<MO_Nkpoint; kloop++){

    k1 = MO_kpoint[kloop][1];
    k2 = MO_kpoint[kloop][2];
    k3 = MO_kpoint[kloop][3];

    for (spin=0; spin<=spinmax; spin++){
      for (orbit=0; orbit<Bulk_Num_LUMOs[kloop]; orbit++){

	/* calc. MO on grids */

	for (GN=0; GN<TNumGrid; GN++){
          RMO_Grid_tmp[GN] = 0.0;
          IMO_Grid_tmp[GN] = 0.0;
	}

        for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
          Gc_AN = M2G[Mc_AN];    
	  Cwan = WhatSpecies[Gc_AN];
	  NO0 = Spe_Total_CNO[Cwan];

          for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
            GN = GridListAtom[Mc_AN][Nc];
            Rn = CellListAtom[Mc_AN][Nc];

            l1 =-atv_ijk[Rn][1];
            l2 =-atv_ijk[Rn][2];
            l3 =-atv_ijk[Rn][3];

            kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;
            si = sin(2.0*PI*kRn);
            co = cos(2.0*PI*kRn);

	    for (i=0; i<NO0; i++){
              ReCoef = co*LUMOs_Coef[kloop][spin][orbit][Gc_AN][i].r 
                      -si*LUMOs_Coef[kloop][spin][orbit][Gc_AN][i].i; 
              ImCoef = co*LUMOs_Coef[kloop][spin][orbit][Gc_AN][i].i
                      +si*LUMOs_Coef[kloop][spin][orbit][Gc_AN][i].r;
              RMO_Grid_tmp[GN] += ReCoef*Orbs_Grid[Mc_AN][i][Nc];
              IMO_Grid_tmp[GN] += ImCoef*Orbs_Grid[Mc_AN][i][Nc];
	    }
          }
	}

        MPI_Reduce(&RMO_Grid_tmp[0], &RMO_Grid[0], TNumGrid, MPI_DOUBLE,
                    MPI_SUM, Host_ID, mpi_comm_level1);
        MPI_Reduce(&IMO_Grid_tmp[0], &IMO_Grid[0], TNumGrid, MPI_DOUBLE,
                    MPI_SUM, Host_ID, mpi_comm_level1);

	for (GN=0; GN<TNumGrid; GN++){
          MO_Grid[GN].r = RMO_Grid[GN];
          MO_Grid[GN].i = IMO_Grid[GN];
	}

	/* output the real part of LUMOs on grids */

        if (myid==Host_ID){ 

  	  sprintf(file1,"%s%s.lumo%i_%i_%i_r%s",
                  filepath,filename,kloop,spin,orbit,CUBE_EXTENSION);

  	  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
            setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

            Print_CubeTitle(fp);
            Print_CubeCData_MO(fp,MO_Grid,"r");
	    fclose(fp);
	  }
	  else{
	    printf("Failure of saving MOs\n");
	    fclose(fp);
	  }

	  /* output the imaginary part of LUMOs on grids */

	  sprintf(file1,"%s%s.lumo%i_%i_%i_i%s",
                  filepath,filename,kloop,spin,orbit,CUBE_EXTENSION);

  	  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
            setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

            Print_CubeTitle(fp);
            Print_CubeCData_MO(fp,MO_Grid,"i");
	    fclose(fp);
	  }
	  else{
	    printf("Failure of saving MOs\n");
 	  }
	}

      } /* orbit */
    }   /* spin  */
  }     /* kloop */


  /****************************************************
    allocation of arrays:

    dcomplex MO_Grid[TNumGrid];
    double RMO_Grid[TNumGrid];
    double IMO_Grid[TNumGrid];
    double RMO_Grid_tmp[TNumGrid];
    double IMO_Grid_tmp[TNumGrid];
  ****************************************************/

  free(MO_Grid);
  free(RMO_Grid);
  free(IMO_Grid);
  free(RMO_Grid_tmp);
  free(IMO_Grid_tmp);
} 




static void Print_CubeTitle(FILE *fp)
{

   int ct_AN;
   int spe; 

   fprintf(fp," SYS1\n SYS1\n");
   fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
	   atomnum,Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
   fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
	   Ngrid1,gtv[1][1],gtv[1][2],gtv[1][3]);
   fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
	   Ngrid2,gtv[2][1],gtv[2][2],gtv[2][3]);
   fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
	   Ngrid3,gtv[3][1],gtv[3][2],gtv[3][3]);

   for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
     spe = WhatSpecies[ct_AN];
     fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf%12.6lf\n",
	     Spe_WhatAtom[spe],
	     Spe_Core_Charge[spe]-InitN_USpin[ct_AN]-InitN_DSpin[ct_AN],
	     Gxyz[ct_AN][1],Gxyz[ct_AN][2],Gxyz[ct_AN][3]);
   }

}




static void Print_CubeData(FILE *fp, char fext[], double *data, double *data1,char *op)
{
  int i,j,k,i1,i2,i3,c;
  int GN,mul,n1,n2,n3,nn1,nn0;
  int cmd,MN,MN0,MN1,MN2,MN3;
  double ****V;
  double *tmp_array0;
  double *tmp_array1;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  char operate[300];
  char fname1[300];
  char fname2[300];
  FILE *fp1,*fp2;
  char buf[fp_bsize];          /* setvbuf */

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (op==NULL)                    { cmd=0; mul=1; }
  else if (strcmp(op,"add")==0)    { cmd=1; mul=2; }
  else if (strcmp(op,"diff")==0)   { cmd=2; mul=2; }
  else {
    printf("Print_CubeData: op=%s not supported\n",op);
    return;
  }

  /****************************************************
   allocation of arrays:

   double V[mul][My_NGrid1_Poisson][Ngrid2][Ngrid3];
  ****************************************************/

  V = (double****)malloc(sizeof(double***)*mul); 
  for (k=0; k<mul; k++){
    V[k] = (double***)malloc(sizeof(double**)*My_NGrid1_Poisson); 
    for (i=0; i<My_NGrid1_Poisson; i++){
      V[k][i] = (double**)malloc(sizeof(double*)*Ngrid2); 
      for (j=0; j<Ngrid2; j++){
        V[k][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
      }
    }
  }

  /****************************************************
                    set V0 and V1
  ****************************************************/

  /* initialize */
  for (k=0; k<mul; k++){
    for (n1=0; n1<My_NGrid1_Poisson; n1++){
      for (n2=0; n2<Ngrid2; n2++){
        for (n3=0; n3<Ngrid3; n3++){
          V[k][n1][n2][n3] = 0.0;
        }
      }
    }
  }

  /* use their densities using MPI */

  for (k=0; k<mul; k++){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /* Isend */
      if (Num_Snd_Grid1[IDS]!=0){

        tmp_array0 = (double*)malloc(sizeof(double)*Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3); 
  
        for (i=0; i<Num_Snd_Grid1[IDS]; i++){ 
	  n1 = Snd_Grid1[IDS][i];
          nn1 = My_Cell0[n1];
          MN1 = nn1*Ngrid2*Ngrid3;
          MN0 = i*Ngrid2*Ngrid3;
          for (n2=0; n2<Ngrid2; n2++){
            MN2 = n2*Ngrid3;

            if (k==0){
              for (n3=0; n3<Ngrid3; n3++){
                MN = MN1 + MN2 + n3;
                tmp_array0[MN0+MN2+n3] = data[MN];
	      }
	    }
            else if (k==1){
              for (n3=0; n3<Ngrid3; n3++){
                MN = MN1 + MN2 + n3;
                tmp_array0[MN0+MN2+n3] = data1[MN];
	      }
	    }

	  }
        }

        MPI_Isend(&tmp_array0[0], Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3,
                  MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /* Recv */
      if (Num_Rcv_Grid1[IDR]!=0){

        tmp_array1 = (double*)malloc(sizeof(double)*Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3); 

        MPI_Recv(&tmp_array1[0], Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3,
                MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        for (i=0; i<Num_Rcv_Grid1[IDR]; i++){ 
	  n1 = Rcv_Grid1[IDR][i];
          nn1 = My_Cell0[n1];
          nn0 = n1 - Start_Grid1[myid];
          MN0 = i*Ngrid2*Ngrid3;
          for (n2=0; n2<Ngrid2; n2++){
            MN2 = n2*Ngrid3;
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN0 + MN2 + n3;
              V[k][nn0][n2][n3] = tmp_array1[MN];
	    }
	  }
        }

        free(tmp_array1);
      }

      if (Num_Snd_Grid1[IDS]!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array0);
      }

    }

    /* use own densities */
    for (n1=Start_Grid1[myid]; n1<=End_Grid1[myid]; n1++){
      nn1 = My_Cell0[n1];
      nn0 = n1 - Start_Grid1[myid]; 
      if (nn1!=-1){
        MN1 = nn1*Ngrid2*Ngrid3;
        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;

          if (k==0){
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              V[k][nn0][n2][n3] = data[MN];
            }    
	  }
          else if (k==1){
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              V[k][nn0][n2][n3] = data1[MN];
            }    
	  }

        }    
      }
    }

  } /* mul */


  /****************************************************
                    output data 
  ****************************************************/

  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    for (n2=0; n2<Ngrid2; n2++){
      for (n3=0; n3<Ngrid3; n3++){

        switch (cmd) {
        case 0:
          fprintf(fp,"%13.3E",V[0][n1][n2][n3]);
          break;
        case 1:
          fprintf(fp,"%13.3E",V[0][n1][n2][n3]+V[1][n1][n2][n3]);
          break;
        case 2:
          fprintf(fp,"%13.3E",V[0][n1][n2][n3]-V[1][n1][n2][n3]);
          break;
        }

        if ((n3+1)%6==0) { fprintf(fp,"\n"); }
      }
      /* avoid double \n\n when Ngrid3%6 == 0  */
      if (Ngrid3%6!=0) fprintf(fp,"\n");
    }
  }

  /****************************************************
   freeing of arrays:

   double V[mul][My_NGrid1_Poisson][Ngrid2][Ngrid3];
  ****************************************************/

  for (k=0; k<mul; k++){
    for (i=0; i<My_NGrid1_Poisson; i++){
      for (j=0; j<Ngrid2; j++){
        free(V[k][i][j]);
      }
      free(V[k][i]);
    }
    free(V[k]);
  }
  free(V);

  /****************************************************
  fclose(fp);
  ****************************************************/

  fclose(fp);

  /****************************************************
                       merge files
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){

#ifdef xt3

    sprintf(fname1,"%s%s%s",filepath,filename,fext); 
    fp1 = fopen(fname1,"a");

    /*
    setvbuf(fp1,buf,_IOFBF,fp_bsize);  
    */

    if (fp1!=NULL){
      remove(fname1); 
      fclose(fp1); 
    }
 
    for (ID=0; ID<numprocs; ID++){
      sprintf(fname1,"%s%s%s",filepath,filename,fext);
      fp1 = fopen(fname1,"a");

      /*
      setvbuf(fp1,buf,_IOFBF,fp_bsize);  
      */

      fseek(fp1,0,SEEK_END);

      sprintf(fname2,"%s%s%s%i",filepath,filename,fext,ID);
      fp2 = fopen(fname2,"r");

      /*
      setvbuf(fp2,buf,_IOFBF,fp_bsize);  
      */

      if (fp2!=NULL){
        for (c=getc(fp2); c!=EOF; c=getc(fp2))  putc(c,fp1); 
	fclose(fp2); 
      }
      fclose(fp1); 
    }  

    for (ID=0; ID<numprocs; ID++){
      sprintf(operate,"%s%s%s%i",filepath,filename,fext,ID);
      remove(operate);
    }

#else

    sprintf(operate,"cat %s%s%s0 > %s%s%s",
            filepath,filename,fext,filepath,filename,fext);
    system(operate);

    for (ID=1; ID<numprocs; ID++){
      sprintf(operate,"cat %s%s%s%i >> %s%s%s",
              filepath,filename,fext,ID,
              filepath,filename,fext);
      system(operate);
    }

    for (ID=0; ID<numprocs; ID++){
      sprintf(operate,"rm %s%s%s%i",filepath,filename,fext,ID);
      system(operate);
    }

#endif

  }

}


static void Print_VectorData(int flagout, FILE *fp, char fext[],
                             double *data0, double *data1,
                             double *data2, double *data3)
{
  int i,j,k,i1,i2,i3,c,GridNum;
  int GN,mul,n1,n2,n3,nn1,nn0;
  int cmd,MN,MN0,MN1,MN2,MN3;
  double ****V;
  double *tmp_array0;
  double *tmp_array1;
  double x,y,z,vx,vy,vz;
  double xmin,ymin,zmin;
  double xmax,ymax,zmax;
  double sden,Cxyz[4],theta,phi;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  char operate[300];
  char fname1[300];
  char fname2[300];
  FILE *fp1,*fp2;
  char *buf;          /* setvbuf */
  double scaxsf2;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  mul = 4;

  /****************************************************
   allocation of arrays:

   double V[mul][My_NGrid1_Poisson][Ngrid2][Ngrid3];
  ****************************************************/

  /* allocate buf */
  buf = malloc(fp_bsize); /* setvbuf */

  V = (double****)malloc(sizeof(double***)*mul); 
  for (k=0; k<mul; k++){
    V[k] = (double***)malloc(sizeof(double**)*My_NGrid1_Poisson); 
    for (i=0; i<My_NGrid1_Poisson; i++){
      V[k][i] = (double**)malloc(sizeof(double*)*Ngrid2); 
      for (j=0; j<Ngrid2; j++){
        V[k][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
      }
    }
  }

  /****************************************************
                    set V0 and V1
  ****************************************************/

  /* initialize */
  for (k=0; k<mul; k++){
    for (n1=0; n1<My_NGrid1_Poisson; n1++){
      for (n2=0; n2<Ngrid2; n2++){
        for (n3=0; n3<Ngrid3; n3++){
          V[k][n1][n2][n3] = 0.0;
        }
      }
    }
  }

  /* use their densities using MPI */ 

  for (k=0; k<mul; k++){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /* Isend */
      if (Num_Snd_Grid1[IDS]!=0){

        tmp_array0 = (double*)malloc(sizeof(double)*Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3); 
  
        for (i=0; i<Num_Snd_Grid1[IDS]; i++){ 
	  n1 = Snd_Grid1[IDS][i];
          nn1 = My_Cell0[n1];
          MN1 = nn1*Ngrid2*Ngrid3;
          MN0 = i*Ngrid2*Ngrid3;
          for (n2=0; n2<Ngrid2; n2++){
            MN2 = n2*Ngrid3;

            if (k==0){
              for (n3=0; n3<Ngrid3; n3++){
                MN = MN1 + MN2 + n3;
                tmp_array0[MN0+MN2+n3] = data0[MN];
	      }
	    }
            else if (k==1){
              for (n3=0; n3<Ngrid3; n3++){
                MN = MN1 + MN2 + n3;
                tmp_array0[MN0+MN2+n3] = data1[MN];
	      }
	    }
            else if (k==2){
              for (n3=0; n3<Ngrid3; n3++){
                MN = MN1 + MN2 + n3;
                tmp_array0[MN0+MN2+n3] = data2[MN];
	      }
	    }
            else if (k==3){
              for (n3=0; n3<Ngrid3; n3++){
                MN = MN1 + MN2 + n3;
                tmp_array0[MN0+MN2+n3] = data3[MN];
	      }
	    }

	  }
        }

        MPI_Isend(&tmp_array0[0], Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3,
                  MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /* Recv */
      if (Num_Rcv_Grid1[IDR]!=0){

        tmp_array1 = (double*)malloc(sizeof(double)*Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3); 

        MPI_Recv(&tmp_array1[0], Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3,
                MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        for (i=0; i<Num_Rcv_Grid1[IDR]; i++){ 
	  n1 = Rcv_Grid1[IDR][i];
          nn1 = My_Cell0[n1];
          nn0 = n1 - Start_Grid1[myid];
          MN0 = i*Ngrid2*Ngrid3;
          for (n2=0; n2<Ngrid2; n2++){
            MN2 = n2*Ngrid3;
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN0 + MN2 + n3;
              V[k][nn0][n2][n3] = tmp_array1[MN];
	    }
	  }
        }

        free(tmp_array1);
      }

      if (Num_Snd_Grid1[IDS]!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array0);
      }

    }

    /* use own densities */
    for (n1=Start_Grid1[myid]; n1<=End_Grid1[myid]; n1++){
      nn1 = My_Cell0[n1];
      nn0 = n1 - Start_Grid1[myid]; 
      if (nn1!=-1){
        MN1 = nn1*Ngrid2*Ngrid3;
        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;

          if (k==0){
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              V[k][nn0][n2][n3] = data0[MN];
            }    
	  }
          else if (k==1){
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              V[k][nn0][n2][n3] = data1[MN];
            }    
	  }
          else if (k==2){
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              V[k][nn0][n2][n3] = data2[MN];
            }    
	  }
          else if (k==3){
            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              V[k][nn0][n2][n3] = data3[MN];
            }    
	  }
        }    
      }
    }

  } /* mul */


  /****************************************************
                    output data 
  ****************************************************/

  /* for gOpenMol */ 

  if (flagout == 0) {

    if (myid==Host_ID){
      fprintf(fp,"3 200\n");
      fprintf(fp,"%4d %4d %4d\n",Ngrid1,Ngrid2,Ngrid3);

      Get_Grid_XYZ(0,Cxyz);
      xmin = Cxyz[1];
      ymin = Cxyz[2];
      zmin = Cxyz[3];
      Get_Grid_XYZ(TNumGrid-1,Cxyz);
      xmax = Cxyz[1];
      ymax = Cxyz[2];
      zmax = Cxyz[3];

      fprintf(fp,"%13.3E %13.3E\n",BohrR*xmin,BohrR*xmax);
      fprintf(fp,"%13.3E %13.3E\n",BohrR*ymin,BohrR*ymax);
      fprintf(fp,"%13.3E %13.3E\n",BohrR*zmin,BohrR*zmax);
    }

    for (n1=0; n1<My_NGrid1_Poisson; n1++){
      nn1 = Start_Grid1[myid] + n1;
      for (n2=0; n2<Ngrid2; n2++){
	for (n3=0; n3<Ngrid3; n3++){

	  GN = nn1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3;    

	  Get_Grid_XYZ(GN,Cxyz);
	  x = Cxyz[1];
	  y = Cxyz[2];
	  z = Cxyz[3];

	  sden  = V[0][n1][n2][n3] - V[1][n1][n2][n3];
	  theta = V[2][n1][n2][n3];
	  phi   = V[3][n1][n2][n3];

	  vx = sden*sin(theta)*cos(phi);
	  vy = sden*sin(theta)*sin(phi);
	  vz = sden*cos(theta);

	  fprintf(fp,"%13.3E %13.3E %13.3E %13.3E %13.3E %13.3E\n",
		  BohrR*x,BohrR*y,BohrR*z,vx,vy,vz);
	}
      }
    }

  } /* gOpenMol */


  /* for XCrysDen */ 

  if (flagout == 1) {

    scaxsf2=0.02;

    if (myid==Host_ID){
      fprintf(fp,"CRYSTAL\n");
      fprintf(fp,"PRIMVEC\n");
      fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[1][1], BohrR*tv[1][2], BohrR*tv[1][3]);
      fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[2][1], BohrR*tv[2][2], BohrR*tv[2][3]);
      fprintf(fp," %10.6f %10.6f %10.6f\n",BohrR*tv[3][1], BohrR*tv[3][2], BohrR*tv[3][3]);
      fprintf(fp,"PRIMCOORD 1\n");

      /* fprintf(fp,"%4d %4d\n",atomnum, 1); */

      GridNum = (My_NGrid1_Poisson+1)*(Ngrid2+1)*(Ngrid3+1)/64;

      fprintf(fp,"%4d %4d\n",GridNum, 1);

      for (k=1; k<=atomnum; k++){
        i = WhatSpecies[k];
        /*j = Spe_WhatAtom[i];*/
        fprintf(fp,"%4d  %8.5f  %8.5f  %8.5f\n",
	        Spe_WhatAtom[i],               
	        Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR
	        );
      }
    }


    for (n1=0; n1<My_NGrid1_Poisson; n1=n1+4){
      nn1 = Start_Grid1[myid] + n1;
      for (n2=0; n2<Ngrid2; n2=n2+4){
	for (n3=0; n3<Ngrid3; n3=n3+4){

	  /* for (n1=0; n1<My_NGrid1_Poisson; n1++){
	     nn1 = Start_Grid1[myid] + n1;
	     for (n2=0; n2<Ngrid2; n2++){
	     for (n3=0; n3<Ngrid3; n3++){ */

	  GN = nn1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3; 

	  Get_Grid_XYZ(GN,Cxyz);
	  x = Cxyz[1];
	  y = Cxyz[2];
	  z = Cxyz[3];

	  sden  = V[0][n1][n2][n3] - V[1][n1][n2][n3];
	  theta = V[2][n1][n2][n3];
	  phi   = V[3][n1][n2][n3];

	  vx = sden*sin(theta)*cos(phi);
	  vy = sden*sin(theta)*sin(phi);
	  vz = sden*cos(theta);

	  fprintf(fp,"X %13.3E %13.3E %13.3E %13.3E %13.3E %13.3E\n",
		  BohrR*x,BohrR*y,BohrR*z,scaxsf2*vx,scaxsf2*vy,scaxsf2*vz);
	}
      }
    }
  } /* XCrysDen */

  /****************************************************
   freeing of arrays:

   double V[mul][My_NGrid1_Poisson][Ngrid2][Ngrid3];
  ****************************************************/

  for (k=0; k<mul; k++){
    for (i=0; i<My_NGrid1_Poisson; i++){
      for (j=0; j<Ngrid2; j++){
        free(V[k][i][j]);
      }
      free(V[k][i]);
    }
    free(V[k]);
  }
  free(V);

  /****************************************************
  fclose(fp);
  ****************************************************/

  fclose(fp);

  /****************************************************
                   merge files
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){

#ifdef xt3

    sprintf(fname1,"%s%s%s",filepath,filename,fext);
    fp1 = fopen(fname1,"a");
    setvbuf(fp1,buf,_IOFBF,fp_bsize);  /* setvbuf */

    if (fp1!=NULL){
      remove(fname1); 
      fclose(fp1); 
    }

    for (ID=0; ID<numprocs; ID++){
      sprintf(fname1,"%s%s%s",filepath,filename,fext);
      fp1 = fopen(fname1,"a");
      setvbuf(fp1,buf,_IOFBF,fp_bsize);  /* setvbuf */
      fseek(fp1,0,SEEK_END);

      sprintf(fname2,"%s%s%s%i",filepath,filename,fext,ID);
      fp2 = fopen(fname2,"r");
      setvbuf(fp2,buf,_IOFBF,fp_bsize);  /* setvbuf */

      if (fp2!=NULL){
        for (c=getc(fp2); c!=EOF; c=getc(fp2))  putc(c,fp1); 
	fclose(fp2); 
      }
      fclose(fp1); 
    }  

    for (ID=0; ID<numprocs; ID++){
      sprintf(operate,"%s%s%s%i",filepath,filename,fext,ID);
      remove(operate);
    }

#else

    sprintf(operate,"cat %s%s%s0 > %s%s%s",
            filepath,filename,fext,filepath,filename,fext);
    system(operate);

    for (ID=1; ID<numprocs; ID++){
      sprintf(operate,"cat %s%s%s%i >> %s%s%s",
              filepath,filename,fext,ID,
              filepath,filename,fext);
      system(operate);
    }

    for (ID=0; ID<numprocs; ID++){
      sprintf(operate,"rm %s%s%s%i",filepath,filename,fext,ID);
      system(operate);
    }

#endif

  }

  /* free buf */
  free(buf);
}




void out_OrbOpt(char *inputfile)
{
  int i,j,natom,po,al,p,wan,L0,Mul0,M0;
  int num,Mc_AN,Gc_AN;
  double sum;
  double ***tmp_coes;
  char file1[YOUSO10] = ".oopt";
  char DirPAO[YOUSO10];
  char ExtPAO[YOUSO10] = ".pao";
  char FN_PAO[YOUSO10];
  char EndLine[YOUSO10] = "<pseudo.atomic.orbitals.L=0";
  char fname2[YOUSO10];
  char command0[YOUSO10];
  double *Tmp_Vec;
  char *sp0;
  FILE *fp,*fp2;
  char *buf;          /* setvbuf */
  int numprocs,myid,ID,tag=999;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* set DirPAO */

  sprintf(DirPAO,"%s/PAO/",DFT_DATA_PATH);

  /* allocate buf */
  buf = malloc(fp_bsize); /* setvbuf */

  /* allocation of Tmp_Vec */
  Tmp_Vec = (double*)malloc(sizeof(double)*List_YOUSO[7]*List_YOUSO[24]);

  fnjoint(filepath,filename,file1);

  if ((fp = fopen(file1,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      wan = WhatSpecies[Gc_AN];

      ID = G2ID[Gc_AN];

      /**********************************************
                      set Tmp_Vec
      ***********************************************/

      if (myid==ID){

        Mc_AN = F_G2M[Gc_AN];

        al = -1;
        num = 0;
        for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	  for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
	    for (M0=0; M0<=2*L0; M0++){
	      al++;
	      for (p=0; p<Spe_Specified_Num[wan][al]; p++){
	        Tmp_Vec[num] = CntCoes[Mc_AN][al][p];
                num++;
	      }
	    }
	  }
        }

        if (myid!=Host_ID){
          MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
          MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
        }
       
      }

      else if (ID!=myid && myid==Host_ID){
        MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
        MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      }

      /**********************************************
                         write 
      ***********************************************/
      
      if (myid==Host_ID){

        fprintf(fp,"\nAtom=%2d\n",Gc_AN);
        fprintf(fp,"Basis specification  %s\n",SpeBasis[wan]);

        fprintf(fp,"Contraction coefficients  p=");
        for (i=0; i<List_YOUSO[24]; i++){
	  fprintf(fp,"       %i  ",i);
        }
        fprintf(fp,"\n");

        al = -1;
        num = 0;

        for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
  	  for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
	    for (M0=0; M0<=2*L0; M0++){
	      al++;

	      fprintf(fp,"Atom=%3d  L=%2d  Mul=%2d  M=%2d  ",Gc_AN,L0,Mul0,M0);
	      for (p=0; p<Spe_Specified_Num[wan][al]; p++){
	        fprintf(fp,"%9.5f ",Tmp_Vec[num]);
                num++;
	      }
	      for (p=Spe_Specified_Num[wan][al]; p<List_YOUSO[24]; p++){
	        fprintf(fp,"  0.00000 ");
	      }
	      fprintf(fp,"\n");
	    }
	  }
        }
      } /* if (myid==Host_ID) */
    }   /* Gc_AN */
    fclose(fp);
  }
  else{
    printf("Failure of saving contraction coefficients of basis orbitals\n");
  }

  /****************************************************
       outputting of contracted orbitals as *.pao   
  ****************************************************/

  if (CntOrb_fileout==1 && Cnt_switch==1 && RCnt_switch==1){

    /****************************************************
       allocation of array:

       tmp_coes[List_YOUSO[25]+1]
               [List_YOUSO[24]]
               [List_YOUSO[24]]
    ****************************************************/

    tmp_coes = (double***)malloc(sizeof(double**)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      tmp_coes[i] = (double**)malloc(sizeof(double*)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
        tmp_coes[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
      }
    }

    for (natom=0; natom<Num_CntOrb_Atoms; natom++){

      Gc_AN = CntOrb_Atoms[natom];
      ID = G2ID[Gc_AN];
      wan = WhatSpecies[Gc_AN];

      if (myid==ID){

        Mc_AN = F_G2M[Gc_AN];

        al = -1;
        num = 0;
        for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	  for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
	    for (M0=0; M0<=2*L0; M0++){
	      al++;
	      for (p=0; p<Spe_Specified_Num[wan][al]; p++){
	        Tmp_Vec[num] = CntCoes[Mc_AN][al][p];
                num++;
	      }
	    }
	  }
        }

        if (myid!=Host_ID){
          MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
          MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
        }
       
      }

      else if (ID!=myid && myid==Host_ID){
        MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
        MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      }

      /****************************************************
          setting of initial vectors using contraction
             coefficients and unit vectors (1.0)
      ****************************************************/

      if (myid==Host_ID){

        al = -1;
        num = 0;

        for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
          for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
            for (M0=0; M0<=2*L0; M0++){
              al++;
      
              if (M0==0){
                for (p=0; p<Spe_Specified_Num[wan][al]; p++){
	          tmp_coes[L0][Mul0][p] = Tmp_Vec[num];
                  num++; 
	        }
                for (p=Spe_Specified_Num[wan][al]; p<Spe_PAO_Mul[wan]; p++){
	          tmp_coes[L0][Mul0][p] = 0.0;
	        }
              }
              else{
                for (p=0; p<Spe_Specified_Num[wan][al]; p++){
                  num++; 
	        }
              } 
	    }
	  }

          for (Mul0=Spe_Num_CBasis[wan][L0]; Mul0<Spe_PAO_Mul[wan]; Mul0++){
            for (p=0; p<Spe_PAO_Mul[wan]; p++) tmp_coes[L0][Mul0][p] = 0.0;
            tmp_coes[L0][Mul0][Mul0] = 1.0;
	  }
        }

        /****************************************************
                            make *.pao
        ****************************************************/

        sprintf(fname2,"%s%s_%i%s",filepath,SpeName[wan],Gc_AN,ExtPAO);
        if ((fp2 = fopen(fname2,"w")) != NULL){

#ifdef xt3
	  setvbuf(fp2,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

  	  fnjoint2(DirPAO,SpeBasisName[wan],ExtPAO,FN_PAO); 

          /****************************************************
                      Log of the orbital optimization
                                  and 
                        contraction coefficients
          ****************************************************/

          fprintf(fp2,"*****************************************************\n");
          fprintf(fp2,"*****************************************************\n");
          fprintf(fp2," The numerical atomic orbitals were generated\n");
          fprintf(fp2," by the variational optimization.\n");
          fprintf(fp2," The original file of PAOs was %s.\n",FN_PAO);
          fprintf(fp2," Basis specification was %s.\n",SpeBasis[wan]);
          fprintf(fp2," The input file was %s.\n",inputfile);
          fprintf(fp2,"*****************************************************\n");
          fprintf(fp2,"*****************************************************\n\n");

  	  fprintf(fp2,"Contraction coefficients  p=");
	  for (i=0; i<List_YOUSO[24]; i++){
	    fprintf(fp2,"       %i  ",i);
	  }
  	  fprintf(fp2,"\n");

	  al = -1;
          num = 0;

  	  for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	    for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
	      for (M0=0; M0<=2*L0; M0++){
	        al++;

 	        fprintf(fp2,"Atom=%3d  L=%2d  Mul=%2d  M=%2d  ",Gc_AN,L0,Mul0,M0);
	        for (p=0; p<Spe_Specified_Num[wan][al]; p++){
		  fprintf(fp2,"%9.5f ",Tmp_Vec[num]);
                  num++;
	        }
	        for (p=Spe_Specified_Num[wan][al]; p<List_YOUSO[24]; p++){
		  fprintf(fp2,"  0.00000 ");
	        }
	        fprintf(fp2,"\n");
	      }
	    }
	  }

          fprintf(fp2,"\n*****************************************************\n");
          fprintf(fp2,"*****************************************************\n\n");

          /****************************************************
                        from the original file
          ****************************************************/

	  if ((fp = fopen(FN_PAO,"r")) != NULL){

	    po = 0;       
	    do{
	      if (fgets(command0,YOUSO10,fp)!=NULL){
	        command0[strlen(command0)-1] = '\0';
                sp0 = strstr(command0,EndLine);
	        if (sp0!=NULL){
		  po = 1;
	        }
	        else {
		  fprintf(fp2,"%s\n",command0); 
	        }
	      }
	      else{ 
	        po = 1;  
	      }
	    }while(po==0);
                 
	    fclose(fp);
	  }
  	  else{
	    printf("Could not find %s\n",FN_PAO);
	    exit(1);
	  }

          /****************************************************
                          contracted orbitals
          ****************************************************/

          for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
            fprintf(fp2,"<pseudo.atomic.orbitals.L=%d\n",L0);
            for (i=0; i<Spe_Num_Mesh_PAO[wan]; i++){
              fprintf(fp2,"%18.15f %18.15f ",Spe_PAO_XV[wan][i],Spe_PAO_RV[wan][i]); 
 	      for (Mul0=0; Mul0<Spe_PAO_Mul[wan]; Mul0++){
                sum = 0.0;
                for (p=0; p<Spe_PAO_Mul[wan]; p++){
                  sum += tmp_coes[L0][Mul0][p]*Spe_PAO_RWF[wan][L0][p][i];
	        }
                fprintf(fp2,"%18.15f ",sum);
              }
              fprintf(fp2,"\n");
            }
            fprintf(fp2,"pseudo.atomic.orbitals.L=%d>\n",L0);
	  }

          for (L0=(Spe_MaxL_Basis[wan]+1); L0<=Spe_PAO_LMAX[wan]; L0++){
            fprintf(fp2,"<pseudo.atomic.orbitals.L=%d\n",L0);
            for (i=0; i<Spe_Num_Mesh_PAO[wan]; i++){
              fprintf(fp2,"%18.15f %18.15f ",Spe_PAO_XV[wan][i],Spe_PAO_RV[wan][i]); 
 	      for (Mul0=0; Mul0<Spe_PAO_Mul[wan]; Mul0++){
                fprintf(fp2,"%18.15f ",Spe_PAO_RWF[wan][L0][Mul0][i]);
	      }
              fprintf(fp2,"\n");
	    }
            fprintf(fp2,"pseudo.atomic.orbitals.L=%d>\n",L0);
	  }

          fclose(fp2);

        }
        else{
          printf("Could not open %s\n",fname2);
          exit(0);
        }

      } /* if (myid==Host_ID) */
    }   /* for (natom=0; natom<Num_CntOrb_Atoms; natom++) */

    /****************************************************
       freeing of arrays:
    ****************************************************/

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
        free(tmp_coes[i][j]);
      }
      free(tmp_coes[i]);
    }
    free(tmp_coes);

  }

  free(Tmp_Vec);
  free(buf);
}




static void Print_CubeData_MO(FILE *fp, double *data, double *data1,char *op)
{
   int i1,i2,i3;
   int GN;
   int cmd;

   if (op==NULL) { cmd=0; }
   else if (strcmp(op,"add")==0)  { cmd=1; }
   else if (strcmp(op,"diff")==0) { cmd=2; }
   else {
       printf("Print_CubeData: op=%s not supported\n",op);
       return;
   }

    for (i1=0; i1<Ngrid1; i1++){
      for (i2=0; i2<Ngrid2; i2++){
        for (i3=0; i3<Ngrid3; i3++){
          GN = i1*Ngrid2*Ngrid3 + i2*Ngrid3 + i3;
          switch (cmd) {
          case 0:
            fprintf(fp,"%13.3E",data[GN]);
            break;
          case 1:
            fprintf(fp,"%13.3E",data[GN]+data1[GN]);
            break;
          case 2:
            fprintf(fp,"%13.3E",data[GN]-data1[GN]);
            break;
          }
          if ((i3+1)%6==0) { fprintf(fp,"\n"); }
        }
        /* avoid double \n\n when Ngrid3%6 == 0  */
        if (Ngrid3%6!=0) fprintf(fp,"\n");
      }
    }
}



static void Print_CubeCData_MO(FILE *fp, dcomplex *data,char *op)
{
   int i1,i2,i3;
   int GN;
   int cmd;

   if (strcmp(op,"r")==0) { cmd=1; }
   else if (strcmp(op,"i")==0) {cmd=2; }
   else {
       printf("Print_CubeCData: op=%s not supported\n",op);
       return;
   }

    for (i1=0; i1<Ngrid1; i1++){
      for (i2=0; i2<Ngrid2; i2++){
        for (i3=0; i3<Ngrid3; i3++){
          GN = i1*Ngrid2*Ngrid3 + i2*Ngrid3 + i3;
          switch (cmd) {
          case 1:
            fprintf(fp,"%13.3E",data[GN].r);
            break;
          case 2:
            fprintf(fp,"%13.3E",data[GN].i);
            break;
          }
          if ((i3+1)%6==0) { fprintf(fp,"\n"); }
        }
        /* avoid double \n\n when Ngrid3%6 == 0  */
        if (Ngrid3%6!=0) fprintf(fp,"\n");
      }
    }
}
