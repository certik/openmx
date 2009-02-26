/*****************************************************************************

  Ver. 3.4 (June/20/2008)

  OpenMX (Open source package for Material eXplorer) is a program package
  for linear scaling density functional calculations of large-scale materials.
  Almost time-consuming parts can be performed in O(N) operations where N is
  the number of atoms (or basis orbitals). Thus, the program might be useful
  for studies of large-scale materials.
  The distribution of this program package follows the practice of
  the GNU General Public Licence (GPL).

  OpenMX is based on 

   *  Local density and generalized gradient approximation (LDA, LSDA, GGA)
      to the exchange-corellation term
   *  Norm-conserving pseudo potentials
   *  Variationally optimized pseudo atomic basis orbitals
   *  Solution of Poisson's equation using FFT
   *  Evaluation of two-center integrals using Fourier transformation
   *  Evaluation of three-center integrals using fixed real space grids
   *  Simple mixing, direct inversion in the interative subspace (DIIS),
      and Guaranteed-reduction Pulay's methods for SCF calculations.
   *  Solution of the eigenvalue problem using O(N) methods
   *  ...

  See also our website (http://www.openmx-square.org/)
  for recent developments.


    **************************************************************
     Copyright

     Taisuke Ozaki

     Present (June/20/2008) official address

       Research Center for Integrated Sciences (RCIS), 
       Japan Advanced Institute of Science and Technology (JAIST)
       Asahidai 1-1, Nomi, Ishikawa 923-1292, Japan

       e-mail: t-ozaki@jaist.go.jp
    **************************************************************
 
*****************************************************************************/

/**********************************************************************
  openmx.c:

     openmx.c is the main routine of OpenMX.

  Log of openmx.c:

     5/Oct/2003  Released by T.Ozaki

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

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif
  
#ifdef TRAN
#include "tran_prototypes.h"
#include "tran_variables.h"
#endif


int main(int argc, char *argv[]) 
{ 
  static int numprocs,myid;
  static int MD_iter,i,j,po,ip;
  static char fileE[YOUSO10] = ".ene"; 
  static char fileDRC[YOUSO10] = ".md";
  static char fileMemory[YOUSO10]; 
  static char fileRestart[YOUSO10];
  static char operate[200];
  double TStime,TEtime;
  static struct stat statbuf; 

  /* for idle CPUs */
  int tag;
  int complete;
  MPI_Request request;
  MPI_Status  status;

  /* MPI initialize */

  mpi_comm_level1 = MPI_COMM_WORLD; 
  MPI_Init(&argc,&argv);
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  Num_Procs = numprocs;

  /* for measuring elapsed time */

  dtime(&TStime);

  /* check argv */

  if (argc==1){
    printf("\nCould not find an input file.\n\n");
    MPI_Finalize(); 
    exit(0);
  } 

  /* initialize Runtest_flag */

  Runtest_flag = 0;

  /****************************************************
    ./openmx -nt # 

    specifies the number of threads in parallelization
    by OpenMP
  ****************************************************/
  
  openmp_threads_num = 1; /* default */

  po = 0;
  if (myid==Host_ID){
    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-nt")==0 ){
        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, mpi_comm_level1);

  if ( (argc-1)<(ip+1) ){
    if (myid==Host_ID){
      printf("cannot find the number of threads\n");
    }
    MPI_Finalize();
    exit(0);
  }

  if ( po==1 ){
    openmp_threads_num = atoi(argv[ip+1]);

    if (openmp_threads_num<=0){ 
      if (myid==Host_ID){
        printf("check the number of threads\n");
      }
      MPI_Finalize();
      exit(0);
    }

  }

  omp_set_num_threads(openmp_threads_num);  

  if (myid==Host_ID){
    printf("\nThe number of threads in each node for OpenMP parallelization is %d.\n\n",openmp_threads_num);
  }

  /****************************************************
    ./openmx -show 

    making of *.out files in order to check whether 
    OpenMX normally runs on many platforms or not.
  ****************************************************/

  if (strcmp(argv[1],"-show")==0){
    Show_DFT_DATA(argv);
    exit(1);
  }

  /****************************************************
    ./openmx -maketest

    making of *.out files in order to check whether 
    OpenMX normally runs on many platforms or not.
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketest")==0){
    Maketest("S",argc,argv);
    exit(1);
  }

  /****************************************************
   ./openmx -runtest

   check whether OpenMX normally runs on many platforms
   or not by comparing the stored *.out and generated
   *.out on your machine.
  ****************************************************/

  else if ( strcmp(argv[1],"-runtest")==0){
    Runtest("S",argc,argv);
  }

  /****************************************************
   ./openmx -maketestL

    making of *.out files in order to check whether 
    OpenMX normally runs for relatively large systems
    on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestL")==0){
    Maketest("L",argc,argv);
    exit(1);
  }

  /****************************************************
   ./openmx -runtestL

   check whether OpenMX normally runs for relatively
   large systems on many platforms or not by comparing
   the stored *.out and generated *.out on your machine.
  ****************************************************/

  else if (strcmp(argv[1],"-runtestL")==0){
    Runtest("L",argc,argv);
  }

  /*******************************************************
   check memory leak by monitoring actual used memory size
  *******************************************************/

  else if ( (argc==2 || argc==3) && strcmp(argv[1],"-mltest")==0){
    Memory_Leak_test(argc,argv);
    exit(1);
  }

  /****************************************************
   ./openmx -maketestG

    making of *.out files in order to check whether 
    OpenMX normally runs for geometry optimization
    on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestG")==0){
    Maketest("G",argc,argv);
    exit(1);
  }

  /****************************************************
   ./openmx -runtestG

   check whether OpenMX normally runs for geometry 
   optimization on many platforms or not by comparing
   the stored *.out and generated *.out on your machine.
  ****************************************************/

  else if (strcmp(argv[1],"-runtestG")==0){
    Runtest("G",argc,argv);
  }

  /*******************************************************
   check consistency between analytic and numerical forces
  *******************************************************/

  else if ( (argc==3 || argc==4) && strcmp(argv[1],"-forcetest")==0){

    if      (strcmp(argv[2],"0")==0) force_flag = 0; 
    else if (strcmp(argv[2],"1")==0) force_flag = 1; 
    else if (strcmp(argv[2],"2")==0) force_flag = 2; 
    else if (strcmp(argv[2],"3")==0) force_flag = 3; 
    else if (strcmp(argv[2],"4")==0) force_flag = 4; 
    else if (strcmp(argv[2],"5")==0) force_flag = 5; 
    else if (strcmp(argv[2],"6")==0) force_flag = 6; 
    else if (strcmp(argv[2],"7")==0) force_flag = 7;
    else if (strcmp(argv[2],"8")==0) force_flag = 8;
    else {
      printf("unsupported flag for -forcetest\n");
      exit(1);
    }

    Force_test(argc,argv);
    exit(1);
  }

  /* allocation of CompTime */
  CompTime = (double**)malloc(sizeof(double*)*numprocs); 
  for (i=0; i<numprocs; i++){
    CompTime[i] = (double*)malloc(sizeof(double)*20); 
    for (j=0; j<20; j++) CompTime[i][j] = 0.0;
  }

  if (myid==Host_ID){  
    printf("\n*******************************************************\n"); 
    printf("*******************************************************\n"); 
    printf(" Welcome to OpenMX   Ver. %s                           \n",Version_OpenMX); 
    printf(" Copyright (C), 2002-2008, T.Ozaki                     \n"); 
    printf(" OpenMX comes with ABSOLUTELY NO WARRANTY.             \n"); 
    printf(" This is free software, and you are welcome to         \n"); 
    printf(" redistribute it under the constitution of the GNU-GPL.\n");
    printf("*******************************************************\n"); 
    printf("*******************************************************\n\n"); 
  } 

  Init_List_YOUSO();
  remake_headfile = 0;
  ScaleSize = 1.2; 

  /****************************************************
                   Read the input file
  ****************************************************/

  /* setup CPU group */
  setup_CPU_group(argv[1]);
  if (myid>=atomnum)  goto LAST_PROC;  /*  to cut off CPUs */

  init_alloc_first();
  CompTime[myid][1] = readfile(argv);

  MPI_Barrier(mpi_comm_level1);

  /* initialize PrintMemory routine */

  sprintf(fileMemory,"%s%s.memory%i",filepath,filename,myid);
  PrintMemory(fileMemory,0,"init"); 
  PrintMemory_Fix();
 
  /* initialize */
  
  init();
  fnjoint(filepath,filename,fileE);
  fnjoint(filepath,filename,fileDRC);

  /* check "-mltest2" mode */

  po = 0;
  if (myid==Host_ID){
    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-mltest2")==0 ){
        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, mpi_comm_level1);

  if ( po==1 ) ML_flag = 1;
  else         ML_flag = 0;

  /* check "-forcetest2" mode */

  po = 0;
  if (myid==Host_ID){
    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-forcetest2")==0 ){
        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, mpi_comm_level1);

  if ( po==1 ){
    force_flag = atoi(argv[ip+1]);
    ForceConsistency_flag = 1;
  }

  /* check force consistency */

  if (ForceConsistency_flag==1){
    Check_Force(argv);
    OutData(argv[1]);
    Merge_LogFile(argv[1]);
    Free_Arrays(0);
    MPI_Finalize();
    return 1;
  }

  /****************************************************
      SCF-DFT calculations and MD and geometrical
      optimization.
  ****************************************************/

  MD_iter = 1;

  do {


    CompTime[myid][2] += truncation(MD_iter,Solver==6,1);
    if (ML_flag==1 && myid==Host_ID) Get_VSZ(MD_iter);  

#ifdef TRAN
    if (Solver==4) {
      TRAN_Calc_GridBound( mpi_comm_level1, atomnum, WhatSpecies, Spe_Atom_Cut1,
                           Ngrid1, Grid_Origin, Gxyz, tv, gtv, Right_tv );

      /* output: TRAN_region[], TRAN_grid_bound */
    }
#endif

    CompTime[myid][3] += DFT(MD_iter,(MD_iter-1)%orbitalOpt_per_MDIter+1);

    if (myid==Host_ID) iterout(MD_iter,MD_TimeStep*(MD_iter-1),fileE,fileDRC);

    /* MD or geometry optimization */

    if (ML_flag==0) CompTime[myid][4] += MD_pac(MD_iter,argv[1]);

    MD_iter++;

  } while(MD_Opt_OK==0 && MD_iter<=MD_IterNumber);


#ifdef TRAN

   if ( TRAN_output_hks ) {
      /* left is dummy */
      TRAN_RestartFile(mpi_comm_level1, "write","left",filepath,TRAN_hksoutfilename);
   }

#endif

  /****************************************************
               calculate Voronoi charge
  ****************************************************/
 
  if (Voronoi_Charge_flag==1) Voronoi_Charge();

  /****************************************************
        calculate Voronoi orbital magnetic moment
  ****************************************************/
 
  if (Voronoi_OrbM_flag==1) Voronoi_Orbital_Moment();
 
  /****************************************************
                  Making of output files
  ****************************************************/

  OutData(argv[1]);

  /****************************************************
  making of a file *.frac for the fractional coordinates
  ****************************************************/

  Make_FracCoord(argv[1]);

  /****************************************************
    write connectivity, Hamiltonian, overlap, density
    matrices, and etc. to a file, filename.scfout 
  ****************************************************/

  if (HS_fileout==1) SCF2File("write",argv[1]);

  /* elapsed time */

  dtime(&TEtime);
  CompTime[myid][0] = TEtime - TStime;
  Output_CompTime();
  for (i=0; i<numprocs; i++){
    free(CompTime[i]);
  }
  free(CompTime);

  /* merge log files */
  Merge_LogFile(argv[1]);

  /* free arrays */

  Free_Arrays(0);

LAST_PROC: 

  if (myid<atomnum) PrintMemory("total",0,"sum");

  printf("\nThe calculation was normally finished. (proc=%3d)\n",myid);

  MPI_Finalize();

  return 0;
}
