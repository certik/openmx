#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"

void TRAN_Output_HKS_Write_Grid(
				MPI_Comm comm1,
				int My_NGrid1_Poisson,
				int Ngrid2,
				int Ngrid3,
				int  *Num_Snd_Grid1,
				int **Snd_Grid1,
				int *Num_Rcv_Grid1,
				int **Rcv_Grid1, 
				int *My_Cell0,
				int *Start_Grid1,
				int *End_Grid1, 
				double *data,
				double *data1,
				FILE *fp
				)
{
  double ****V;
  int mul;
  int k,i,j,ID,n1,n2,n3,nn,nn0;
  int nn1,MN1,MN0,MN2,MN;
  double *tmp_array0, *tmp_array1;
  int tag=99;

  int myid,numprocs;
  MPI_Request request;
  MPI_Status  stat;

  MPI_Comm_rank(comm1,&myid);
  MPI_Comm_size(comm1,&numprocs);

  k=0;
  mul=1;

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

  k=0;

  for (ID=0; ID<numprocs; ID++){

    /* Isend */
    if (Num_Snd_Grid1[ID]!=0){

      tmp_array0 = (double*)malloc(sizeof(double)*Num_Snd_Grid1[ID]*Ngrid2*Ngrid3); 
  
      for (i=0; i<Num_Snd_Grid1[ID]; i++){ 
	n1 = Snd_Grid1[ID][i];
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

      MPI_Isend(&tmp_array0[0], Num_Snd_Grid1[ID]*Ngrid2*Ngrid3,
		MPI_DOUBLE, ID, tag, comm1, &request);
    }

    /* Recv */
    if (Num_Rcv_Grid1[ID]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_Rcv_Grid1[ID]*Ngrid2*Ngrid3); 

      MPI_Recv(&tmp_array1[0], Num_Rcv_Grid1[ID]*Ngrid2*Ngrid3,
	       MPI_DOUBLE, ID, tag, comm1, &stat);

      for (i=0; i<Num_Rcv_Grid1[ID]; i++){ 
	n1 = Rcv_Grid1[ID][i];
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

    if (Num_Snd_Grid1[ID]!=0){
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

  { double *v;
  int m;
  v= (double*)malloc(sizeof(double)*Ngrid2*Ngrid3);
  for (ID=0;ID<numprocs;ID++) {
    for (n1=Start_Grid1[ID]; n1<=End_Grid1[ID]; n1++){
      nn0= n1- Start_Grid1[ID];
      if (myid==ID) {
        m=0;
        for (n2=0; n2<Ngrid2; n2++){
	  for (n3=0; n3<Ngrid3; n3++){
	    v[m++]= V[k][nn0][n2][n3];
	  }
        }
        if (myid!=Host_ID)
        MPI_Send(v, Ngrid2*Ngrid3, MPI_DOUBLE, Host_ID, tag,comm1);
      }
      if (myid==Host_ID) {
        if (myid!=ID) 
	MPI_Recv(v,Ngrid2*Ngrid3, MPI_DOUBLE, ID, tag, comm1,&stat);
	fwrite(v, sizeof(double), Ngrid2*Ngrid3, fp); 
      }
    }       
  }      

  if (myid==Host_ID) {
    printf("The last of data=%le\n",v[Ngrid2*Ngrid3-1]);
  }
  free(v);

  }


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



}
