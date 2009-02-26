/**********************************************************************
  readfile.c:

     readfile.c is a subrutine to read a input file or restart file.

  Log of readfile.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef TRAN
#include "tran_prototypes.h"
#endif



double readfile(char *argv[])
{ 
  double time0;
  double TStime,TEtime;
  FILE *fp;
  int numprocs,myid; 
  char fileMemory[YOUSO10]; 
  char buf[fp_bsize];          /* setvbuf */

  dtime(&TStime);

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
           Read a input file or restart file.
  ****************************************************/

  if ((fp = fopen(argv[1],"r")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    Input_std(argv[1]);

    fclose(fp);
  }
  else{
    printf("Failure of reading the input file.\n");
    exit(0);
  }


#if 0
  /* separate CPUs */
  if ( atomnum < numprocs )  {
    int *new_ranks; 
    int i;
    MPI_Group  new_group,old_group; 

    if (myid==Host_ID) {
       printf("****************************************\n");
       printf("Cut off CPUs, New group contains %d CPUs\n",atomnum); 
       printf("****************************************\n");

    }

    new_ranks=(int*)malloc(sizeof(int)*atomnum); 
    for (i=0;i<atomnum; i++) {
     new_ranks[i]=i; /* a new group is made of original rank=0:atomnum-1 */
    }

    MPI_Comm_group(MPI_COMM_WORLD, &old_group);
    /* define a new group */
    MPI_Group_incl(old_group,atomnum,new_ranks,&new_group);
    MPI_Comm_create(MPI_COMM_WORLD,new_group,&mpi_comm_level1);

    free(new_ranks); /* never forget cleaning! */

  }
  /* default mpi_comm_level1 is already set in main() */
#endif

  if (myid<atomnum) {

    Allocate_Arrays(2);
    Set_Allocate_Atom2CPU(0,1,0); /* for species */
    SetPara_DFT();

#ifdef TRAN 

    if ( TRAN_Check_Region_Lead(atomnum, WhatSpecies, Spe_Atom_Cut1, Gxyz, tv)==0 ) {
      printf("\n\nERROR: PAOs of lead atoms can overlap only to the next nearest region.\n\n");
      MPI_Finalize();
      exit(1);
    }
    if ( Solver==4 ) {
      if ( TRAN_Check_Region(atomnum, WhatSpecies, Spe_Atom_Cut1, Gxyz)==0 ) {
	printf("\n\nERROR: PAOs of atoms of L|C|R can overlap only to the next nearest region.\n\n");
	MPI_Finalize();
	exit(1);
      }
    }
#endif

     
    Set_Allocate_Atom2CPU(0,0,0); /* a simple division for atoms (get Matomnum) */

    if (Solver!=6) { /*  except for GDC */

    /*****************************
      final input
      0: atomnum, 
      1: the neighbor 
      2: elapsed time 
    *****************************/    

      Set_Allocate_Atom2CPU(1,0,0); 
    }
  } 

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}




