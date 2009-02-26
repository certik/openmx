/**********************************************************************
  TRAN_DFT_Dosout.c:

  TRAN_DFT_Dosout.c is a subroutine to calculate density of states 
  of a central region with left and right infinite leads based on
  a non-equilibrium Green's function method. 

  Log of TRAN_DFT_Dosout.c:

     21/Jan/2006  Released by T.Ozaki

***********************************************************************/

#define MEASURE_TIME  0

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"

 
void dtime(double *);

int Get_OneD_HS_Col(int set_flag, double ****RH, double *H1, int *MP, 
                    int *order_GA, int *My_NZeros, int *is1, int *is2);


static void TRAN_Add_MAT(
    int mode,
    int NUM_c,
    dcomplex w_weight,
    dcomplex *v,
    dcomplex *out
)
{
  int i;
  double dum; 

  switch (mode) {
  case 0:
    for (i=0;i<NUM_c*NUM_c;i++) {
      out[i].r =  0.0;
      out[i].i =  0.0;
    }
    break;
  case 1:
    for (i=0;i<NUM_c*NUM_c;i++) {
      out[i].r += ( v[i].r*w_weight.r - v[i].i*w_weight.i );
      out[i].i += ( v[i].r*w_weight.i + v[i].i*w_weight.r );
    }
    break;
  case 2:
    for (i=0;i<NUM_c*NUM_c;i++) {
      /* real(G^less*weight) */
      dum = v[i].r*w_weight.r - v[i].i*w_weight.i;
      out[i].r += dum;
    }
    break;

  }
}






/*
 *   calculate CDM from nh, CntOLP, ... 
 *
 *       an alternative routine for Cluster_DFT or Band_DFT 
 *
 */


static void TRAN_DFT_Kdependent(
			  /* input */
			  MPI_Comm comm1,
			  int level_stdout,
			  int iter,
			  int SpinP_switch,
                          double k2,
                          double k3,
                          int k_op,
                          int *order_GA,
                          double **H1,
                          double *S1,
			  double *****nh,  /* H */
			  double *****ImNL, /* not used, s-o coupling */
			  double ****CntOLP, 
			  int atomnum,
			  int Matomnum,
			  int *WhatSpecies,
			  int *Spe_Total_CNO,
			  int *FNAN,
			  int **natn, 
			  int **ncn,
			  int *M2G, 
			  int *G2ID, 
			  int **atv_ijk,
			  int *List_YOUSO,
			  /* output */
			  float ****Dos,   /* output, DOS */
			  double *****EDM,  /* not used */
			  double Eele0[2], double Eele1[2]) /* not used */
#define GC_ref(i,j) GC[ NUM_c*((j)-1) + (i)-1 ] 


{
  int i,j,k,iside; 
  int *MP;
  int  iw,iw_method;
  dcomplex w, w_weight;
  dcomplex *GCR,*GCA;
  dcomplex *GRL,*GRR,*SigmaL, *SigmaR; 
  dcomplex *v1;
  double dum,sum,tmpr,tmpi;
  double TStime,TEtime;
  int MA_AN, GA_AN, wanA, tnoA, Anum;
  int LB_AN, GB_AN, wanB, tnoB, Bnum; 
  int l1,l2,l3,Rn;

  static int numprocs,myid,ID;
  int **iwIdx, Miwmax, Miw,iw0; 
  double *r_energy,de;
  
  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);
  
  /* parallel setup */
  iwIdx=(int**)malloc(sizeof(int*)*numprocs);
  Miwmax = (tran_dos_energydiv)/numprocs+1;
  for (i=0;i<numprocs;i++) {
    iwIdx[i]=(int*)malloc(sizeof(int)*Miwmax);
  }

  TRAN_Distribute_Node_Idx(0, tran_dos_energydiv-1, numprocs, Miwmax,
                           iwIdx); /* output */

  /* set up energies where DOS is calculated */
  r_energy = (double*)malloc(sizeof(double)*tran_dos_energydiv);

  de = (tran_dos_energyrange[1]-tran_dos_energyrange[0])/(double)tran_dos_energydiv;
  for (i=0; i<tran_dos_energydiv; i++) {
    r_energy[i] = tran_dos_energyrange[0] + de*(double)i;
  }

  /* setup MP */
  TRAN_Set_MP(0, atomnum, WhatSpecies, Spe_Total_CNO, &NUM_c, MP);
  MP = (int*)malloc(sizeof(int)*(NUM_c+1));
  TRAN_Set_MP(1, atomnum, WhatSpecies, Spe_Total_CNO, &NUM_c, MP);
  
  /* initialize */
  TRAN_Set_Value_double(SCC,NUM_c*NUM_c,    0.0,0.0);
  TRAN_Set_Value_double(SCL,NUM_c*NUM_e[0], 0.0,0.0);
  TRAN_Set_Value_double(SCR,NUM_c*NUM_e[1], 0.0,0.0);
  for (k=0; k<=SpinP_switch; k++) {
    TRAN_Set_Value_double(HCC[k],NUM_c*NUM_c,    0.0,0.0);
    TRAN_Set_Value_double(HCL[k],NUM_c*NUM_e[0], 0.0,0.0);
    TRAN_Set_Value_double(HCR[k],NUM_c*NUM_e[1], 0.0,0.0);
  }

  /* set Hamiltonian and overlap matrices of left and right leads */

  TRAN_Set_SurfOverlap(comm1,"left", k2, k3);
  TRAN_Set_SurfOverlap(comm1,"right",k2, k3);

  /* set CC, CL and CR */

  TRAN_Set_CentOverlap(   comm1,
                          3,
                          SpinP_switch, 
                          k2,
                          k3,
                          order_GA,
                          H1,
                          S1,
                          nh, /* input */
                          CntOLP, /* input */
                          atomnum,
			  Matomnum,
			  M2G,
			  G2ID,
                          WhatSpecies,
                          Spe_Total_CNO,
                          FNAN,
                          natn,
                          ncn,
                          atv_ijk);

  GCR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  GCA = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  GRL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[0]* NUM_e[0]);
  GRR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[1]* NUM_e[1]);
  SigmaL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  v1 = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  
  /*parallel global iw 0: tran_dos_energydiv-1 */
  /*parallel local  Miw 0:Miwmax-1 */
  /*parllel variable iw=iwIdx[myid][Miw] */

  for (Miw=0; Miw<Miwmax; Miw++) {

    iw = iwIdx[myid][Miw];

    for (k=0; k<=SpinP_switch; k++) {

      if (iw>=0) {

        /* w = w.r + i w.i */
        
        w.r = r_energy[iw];
        w.i = tran_dos_energyrange[2];
        
        iside=0;

        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max,
                                   tran_surfgreen_eps, GRL);
        

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL);

        iside=1;

        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				   S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max,
                                   tran_surfgreen_eps, GRR);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR);

        TRAN_Calc_CentGreen(w, NUM_c, SigmaL,SigmaR, HCC[k], SCC, GCR);

        /* w = w.r - i w.i */

        if (TRAN_dos_Kspace_grid2!=1 || TRAN_dos_Kspace_grid3!=1){

          w.r = r_energy[iw];
          w.i =-tran_dos_energyrange[2];

	  iside=0;

	  TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				     S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max,
				     tran_surfgreen_eps, GRL);

	  TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL);

	  iside=1;

	  TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
				     S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max,
				     tran_surfgreen_eps, GRR);

	  TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR);

	  TRAN_Calc_CentGreen(w, NUM_c, SigmaL,SigmaR, HCC[k], SCC, GCA);
	}

        /***********************************************
                  calculate density of states
        ***********************************************/

	for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {

#define v_idx(i,j)   ( ((j)-1)*NUM_c + (i)-1 ) 

	  wanA = WhatSpecies[GA_AN];
	  tnoA = Spe_Total_CNO[wanA];
	  Anum = MP[GA_AN];

          for (i=0; i<tnoA; i++) {

            sum = 0.0;

	    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	      GB_AN = natn[GA_AN][LB_AN];
 	      Rn = ncn[GA_AN][LB_AN];
	      wanB = WhatSpecies[GB_AN];
	      tnoB = Spe_Total_CNO[wanB];
	      Bnum = MP[GB_AN];

              l1 = atv_ijk[Rn][1];
              l2 = atv_ijk[Rn][2];
              l3 = atv_ijk[Rn][3];

              if (TRAN_dos_Kspace_grid2==1 && TRAN_dos_Kspace_grid3==1){
                if (l1==0 && l2==0 && l3==0){
	          for (j=0; j<tnoB; j++) {
		    sum += -GCR[ v_idx( Anum+i, Bnum+j) ].i * SCC[ v_idx( Anum+i, Bnum+j) ].r;
	          }
	        }
	      }
               
              else {
                if (l1==0 && l2==0 && l3==0){
	          for (j=0; j<tnoB; j++){
                    
                    tmpr = 0.5*(GCR[ v_idx( Anum+i, Bnum+j) ].i - GCA[ v_idx( Anum+i, Bnum+j) ].i); 
                    tmpi =-0.5*(GCR[ v_idx( Anum+i, Bnum+j) ].r - GCA[ v_idx( Anum+i, Bnum+j) ].r); 

		    sum += -(tmpr*SCC[ v_idx( Anum+i, Bnum+j) ].r
                            -tmpi*SCC[ v_idx( Anum+i, Bnum+j) ].i);
	          }
	        }
              } 
	    }

            Dos[iw][k][GA_AN][i] += (float)k_op*sum/PI;

          } 
	} 
      } /* iw>=0 */
    }   /* for k */
  }     /* iw    */

  /* free arrays */

  free(v1);
  free(SigmaR);
  free(SigmaL);
  free(GRR);
  free(GRL);
  free(GCR);
  free(GCA);

  free(MP);

  for (i=0;i<numprocs;i++) {
    free(iwIdx[i]);
  }
  free(iwIdx);

  free(r_energy);
}










double TRAN_DFT_Dosout(
		/* input */
                MPI_Comm comm1,
                int level_stdout,
		int iter, 
		int SpinP_switch,
		double *****nh,   /* H */
		double *****ImNL, /* not used, s-o coupling */
		double ****CntOLP, 
		int atomnum,
		int Matomnum,
		int *WhatSpecies,
		int *Spe_Total_CNO,
		int *FNAN,
		int **natn, 
		int **ncn,
		int *M2G, 
		int *G2ID, 
		int **atv_ijk,
		int *List_YOUSO,
                int **Spe_Num_CBasis,
                int SpeciesNum,
                char *filename,
                char *filepath,
		/* output */
		double *****CDM,  /* not used */
		double *****EDM,  /* not used */
		double Eele0[2], double Eele1[2]) /* not used */
{
  int numprocs,myid,ID;
  int i2,i3,k_op,ik;
  int i,j,spin,MA_AN,GA_AN,wanA,tnoA;
  int LB_AN,GB_AN,wanB,tnoB;
  int **op_flag;
  double k2,k3,tmp;
  double TStime,TEtime;
  float ****Dos;

  int *MP; 
  int *order_GA;
  int *My_NZeros;
  int *SP_NZeros;
  int *SP_Atoms;
  int size_H1;
  double **H1,*S1;

  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  dtime(&TStime);

  if (myid==Host_ID){
    printf("<TRAN_DFT_Dosout>\n"); fflush(stdout);
  }

  /* allocate Dos */

  Dos = (float****)malloc(sizeof(float***)*tran_dos_energydiv);
  for (ik=0; ik<tran_dos_energydiv; ik++) {
    Dos[ik] = (float***)malloc(sizeof(float**)*(SpinP_switch+1) );
    for (spin=0; spin<=SpinP_switch; spin++) {
      Dos[ik][spin] = (float**)malloc(sizeof(float*)*(atomnum+1) );
      for (GA_AN=0; GA_AN<=atomnum; GA_AN++) {
        Dos[ik][spin][GA_AN] = (float*)malloc(sizeof(float)*List_YOUSO[7] );
        for (i=0; i<List_YOUSO[7]; i++)  Dos[ik][spin][GA_AN][i] = 0.0;
      }
    }
  }

  /*************************************************************
   one-dimensitonalize H and S and store them in a compact form  
  *************************************************************/

  MP = (int*)malloc(sizeof(int)*(atomnum+1));
  order_GA = (int*)malloc(sizeof(int)*(atomnum+1));

  My_NZeros = (int*)malloc(sizeof(int)*numprocs);
  SP_NZeros = (int*)malloc(sizeof(int)*numprocs);
  SP_Atoms = (int*)malloc(sizeof(int)*numprocs);

  size_H1 = Get_OneD_HS_Col(0, nh[0], S1, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

  H1 = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (spin=0; spin<(SpinP_switch+1); spin++){
    H1[spin] = (double*)malloc(sizeof(double)*size_H1);
  }

  S1 = (double*)malloc(sizeof(double)*size_H1);

  for (spin=0; spin<(SpinP_switch+1); spin++){
    size_H1 = Get_OneD_HS_Col(1, nh[spin], H1[spin], MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);
  }

  size_H1 = Get_OneD_HS_Col(1, CntOLP, S1, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

  /* loop for k2 and k3 */

  op_flag = (int**)malloc(sizeof(int*)*TRAN_dos_Kspace_grid2); 
  for (i2=0; i2<TRAN_dos_Kspace_grid2; i2++){
    op_flag[i2] = (int*)malloc(sizeof(int)*TRAN_dos_Kspace_grid3); 
    for (i3=0; i3<TRAN_dos_Kspace_grid3; i3++){
      op_flag[i2][i3] = 0;
    }
  }

  for (i2=0; i2<TRAN_dos_Kspace_grid2; i2++){

    k2 = -0.5 + (2.0*(double)i2+1.0)/(2.0*(double)TRAN_dos_Kspace_grid2) + Shift_K_Point;

    for (i3=0; i3<TRAN_dos_Kspace_grid3; i3++){

      k3 = -0.5 + (2.0*(double)i3+1.0)/(2.0*(double)TRAN_dos_Kspace_grid3) - Shift_K_Point;

      if (op_flag[i2][i3]==0){  

        if ( (TRAN_dos_Kspace_grid2-1-i2)==i2 && (TRAN_dos_Kspace_grid3-1-i3)==i3 ){
          k_op = 1;
        }
        else {
          k_op = 2;
        }

        TRAN_DFT_Kdependent(comm1, level_stdout, iter, SpinP_switch, k2, k3, k_op,
                            order_GA, H1, S1,
                            nh, ImNL, CntOLP,
	  	  	    atomnum, Matomnum, WhatSpecies, Spe_Total_CNO, FNAN,
			    natn, ncn, M2G, G2ID, atv_ijk, List_YOUSO, Dos, EDM, Eele0, Eele1);

        op_flag[i2][i3] = 1;
        op_flag[TRAN_dos_Kspace_grid2-1-i2][TRAN_dos_Kspace_grid3-1-i3] = 1;
      }

    }
  }

  /* free arrays */

  for (i2=0; i2<TRAN_dos_Kspace_grid2; i2++){
    free(op_flag[i2]);
  }
  free(op_flag);

  free(MP);
  free(order_GA);

  free(My_NZeros);
  free(SP_NZeros);
  free(SP_Atoms);

  for (spin=0; spin<(SpinP_switch+1); spin++){
    free(H1[spin]);
  }
  free(H1);
 
  free(S1);

  /*******************************************************
                 distribution of Dos by MPI 
  *******************************************************/

  {
    int *my_ik2ID,*ik2ID,**iwIdx;
    int Miw,iw,Miwmax;
    double *r_energy,de;
    float *Dos0;

    MPI_Barrier(comm1);

    /* allocation of arrays */

    Miwmax = (tran_dos_energydiv)/numprocs+1;

    iwIdx=(int**)malloc(sizeof(int*)*numprocs);
    for (i=0;i<numprocs;i++) {
      iwIdx[i]=(int*)malloc(sizeof(int)*Miwmax);
    }

    r_energy = (double*)malloc(sizeof(double)*tran_dos_energydiv);

    my_ik2ID = (int*)malloc(sizeof(int)*tran_dos_energydiv);
    ik2ID = (int*)malloc(sizeof(int)*tran_dos_energydiv);

    Dos0 = (float*)malloc(sizeof(float)*(SpinP_switch+1)*(atomnum+1)*List_YOUSO[7]);     

    TRAN_Distribute_Node_Idx(0, tran_dos_energydiv-1, numprocs, Miwmax,
			     iwIdx); /* output */

    /* set up energies where DOS is calculated */

    de = (tran_dos_energyrange[1]-tran_dos_energyrange[0])/(double)tran_dos_energydiv;
    for (i=0; i<tran_dos_energydiv; i++) {
      r_energy[i] = tran_dos_energyrange[0] + de*(double)i;
    }

    for (ik=0; ik<tran_dos_energydiv; ik++) {
      my_ik2ID[ik] = 0;
      ik2ID[ik] = 0;
    }

    for (Miw=0; Miw<Miwmax; Miw++) {
      iw = iwIdx[myid][Miw];
      if (iw>=0) {
	my_ik2ID[iw] = myid; 
      }
    }

    MPI_Barrier(comm1);
    MPI_Allreduce(&my_ik2ID[0], &ik2ID[0], tran_dos_energydiv, MPI_INT, MPI_SUM, comm1);

    for (ik=0; ik<tran_dos_energydiv; ik++) {
      ID = ik2ID[ik];

      if (myid==ID){
        j = 0;
	for (spin=0; spin<=SpinP_switch; spin++) {
	  for (GA_AN=0; GA_AN<=atomnum; GA_AN++) {
	    for (i=0; i<List_YOUSO[7]; i++){
              Dos0[j] = Dos[ik][spin][GA_AN][i];
              j++; 
	    }
	  }
	}
      }

      MPI_Bcast(Dos0, (SpinP_switch+1)*(atomnum+1)*List_YOUSO[7], MPI_FLOAT, ID, comm1);

      if (myid!=ID){
        j = 0;
	for (spin=0; spin<=SpinP_switch; spin++) {
	  for (GA_AN=0; GA_AN<=atomnum; GA_AN++) {
	    for (i=0; i<List_YOUSO[7]; i++){
              Dos[ik][spin][GA_AN][i] = Dos0[j];
              j++; 
	    }
	  }
	}
      }
    }

    MPI_Barrier(comm1);

    /* freeing of arrays */

    for (i=0; i<numprocs; i++) {
      free(iwIdx[i]);
    }
    free(iwIdx);

    free(r_energy);

    free(my_ik2ID);
    free(ik2ID);
    free(Dos0);
  }

  /**********************************************************
   divide Dos by TRAN_dos_Kspace_grid2*TRAN_dos_Kspace_grid3
  **********************************************************/

  tmp = 1.0/(double)(TRAN_dos_Kspace_grid2*TRAN_dos_Kspace_grid3);

  for (ik=0; ik<tran_dos_energydiv; ik++) {
    for (spin=0; spin<=SpinP_switch; spin++) {
      for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {
        wanA = WhatSpecies[GA_AN];
        tnoA = Spe_Total_CNO[wanA];
        for (i=0; i<tnoA; i++){
          Dos[ik][spin][GA_AN][i] *= tmp;
	}
      }
    }
  }

  /**********************************************************
                     save Dos to a file
  **********************************************************/

  if (myid==Host_ID){

    FILE *fp_eig, *fp_ev;
    char file_eig[YOUSO10],file_ev[YOUSO10];
    int l,MaxL;
    double *r_energy,de;
    int i_vec[10];

    r_energy = (double*)malloc(sizeof(double)*tran_dos_energydiv);

    de = (tran_dos_energyrange[1]-tran_dos_energyrange[0])/(double)tran_dos_energydiv;
    for (i=0; i<tran_dos_energydiv; i++) {
      r_energy[i] = tran_dos_energyrange[0] + de*(double)i;
    }

    /* write *.Dos.val */

    sprintf(file_eig,"%s%s.Dos.val",filepath,filename);

    if ( (fp_eig=fopen(file_eig,"w"))==NULL ) {
      printf("can not open a file %s\n",file_eig);
    }
    else {

      printf("  write eigenvalues\n");

      fprintf(fp_eig,"mode        6\n");
      fprintf(fp_eig,"NonCol      0\n");
      fprintf(fp_eig,"N           %d\n",NUM_c);
      fprintf(fp_eig,"Nspin       %d\n",SpinP_switch);
      fprintf(fp_eig,"Erange      %lf %lf\n",tran_dos_energyrange[0],tran_dos_energyrange[1]);
      fprintf(fp_eig,"Kgrid       %d %d %d\n",1,1,1);
      fprintf(fp_eig,"atomnum     %d\n",atomnum);
      fprintf(fp_eig,"<WhatSpecies\n");
      for (i=1; i<=atomnum; i++) {
        fprintf(fp_eig,"%d ",WhatSpecies[i]);
      }
      fprintf(fp_eig,"\nWhatSpecies>\n");
      fprintf(fp_eig,"SpeciesNum     %d\n",SpeciesNum);
      fprintf(fp_eig,"<Spe_Total_CNO\n");
      for (i=0;i<SpeciesNum;i++) {
        fprintf(fp_eig,"%d ",Spe_Total_CNO[i]);
      }
      fprintf(fp_eig,"\nSpe_Total_CNO>\n");
      MaxL=4;
      fprintf(fp_eig,"MaxL           %d\n",4);
      fprintf(fp_eig,"<Spe_Num_CBasis\n");
      for (i=0;i<SpeciesNum;i++) {
        for (l=0;l<=MaxL;l++) {
	  fprintf(fp_eig,"%d ",Spe_Num_CBasis[i][l]);
        }
        fprintf(fp_eig,"\n");
      }
      fprintf(fp_eig,"Spe_Num_CBasis>\n");
      fprintf(fp_eig,"ChemP       %lf\n",ChemP_e[0]);

      fprintf(fp_eig,"irange      %d %d\n",0,tran_dos_energydiv-1);
      fprintf(fp_eig,"<Eigenvalues\n");
      for (spin=0; spin<=SpinP_switch; spin++) {
        fprintf(fp_eig,"%d %d %d ",0,0,0);
        for (ik=0; ik<tran_dos_energydiv; ik++) {
          fprintf(fp_eig,"%lf ",r_energy[ik]);
	}
        fprintf(fp_eig,"\n");
      }  
      fprintf(fp_eig,"Eigenvalues>\n");

      fclose(fp_eig);
    }

    /* write *.Dos.vec */

    printf("  write eigenvectors\n");

    sprintf(file_ev,"%s%s.Dos.vec",filepath,filename);

    if ( (fp_ev=fopen(file_ev,"w"))==NULL ) {
      printf("can not open a file %s\n",file_ev);
    }
    else {

      for (spin=0; spin<=SpinP_switch; spin++) {
        for (ik=0; ik<tran_dos_energydiv; ik++) {

          i_vec[0]=i_vec[1]=i_vec[2]=0;
          if (myid==Host_ID) fwrite(i_vec,sizeof(int),3,fp_ev);

          for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {
	    wanA = WhatSpecies[GA_AN];
	    tnoA = Spe_Total_CNO[wanA];

            fwrite(Dos[ik][spin][GA_AN],sizeof(float),tnoA,fp_ev);
	  }
	}
      }

      fclose(fp_ev);
    }

    /* free arrays */

    free(r_energy);
  }

  /* free Dos */

  for (ik=0; ik<tran_dos_energydiv; ik++) {
    for (spin=0; spin<=SpinP_switch; spin++) {
      for (GA_AN=1; GA_AN<=atomnum; GA_AN++) {
        free(Dos[ik][spin][GA_AN]);
      }
      free(Dos[ik][spin]);
    }
    free(Dos[ik]);
  }
  free(Dos);

  /* for elapsed time */
  dtime(&TEtime);

  /*
  if (myid==Host_ID){
    printf("TRAN_DFT_Dosout time=%12.7f\n",TEtime - TStime);
  }
  */

  return TEtime - TStime;
}
