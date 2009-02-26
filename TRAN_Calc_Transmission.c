#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"
#include "tran_variables.h"

/*
 * calculate transmisson 
 *
 * implicit input:  HCC, HCL, HCR, SCC, SCL, SCR for the C region 
 *    
 * if  (tran_transmission_on) 
 *  at w= 
 * for (iw=0;iw<tran_transmission_energydiv ; iw++) 
 *   w.r = tran_transmission_energyrange[0]+
 *       (tran_transmission_energyrange[1]-tran_transmission_energyrange[0])*
 *       (double)iw/(tran_transmission_energydiv-1);
 *   w.i = tran_transmission_energyrange[2];
 *
 * if (tran_transmission_iv_on)
 *   w.r = tran_transmission_iv_energyrange[0]+
 *     (tran_transmission_iv_energyrange[1]-tran_transmission_iv_energyrange[0])*
 *         (iw+0.5)/(tran_transmission_energydiv-1) 
 *   w.i = tran_transmission_iv_energyrange[2]
 *
 * output: tran_transmission[k][iw]   (dcomplex)
 * 
 */

double TRAN_Calc_Transmission(
		/* input */
		int iter, 
		int SpinP_switch,
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
		int **atv_ijk,
		int *List_YOUSO
)
#define GC_ref(i,j) GC[ NUM_c*((j)-1) + (i)-1 ] 


{
  int i,j,k,iside; 
  int *MP;
  int  iw;
  dcomplex w;
  dcomplex *GC,*GRL,*GRR,*SigmaL, *SigmaR; 
  dcomplex *v1,*v2;
  dcomplex value; 
  double dum; 

  int GA_AN, wanA, tnoA, Anum, LB_AN, GB_AN, wanB, tnoB, Bnum; 



  printf("<TRAN_Calc_Transmission> in\n");

  /* setup MP */
  TRAN_set_MP(0, atomnum, WhatSpecies, Spe_Total_CNO,
	      &NUM_c, MP);
  MP = (int*)malloc(sizeof(int)*(NUM_c+1));
  TRAN_set_MP(1, atomnum, WhatSpecies, Spe_Total_CNO,
	      &NUM_c, MP);



  /* allocate */
  SCC = (double*)malloc(sizeof(double)*NUM_c*NUM_c);
  SCL = (double*)malloc(sizeof(double)*NUM_c*NUM_e[0]);
  SCR = (double*)malloc(sizeof(double)*NUM_c*NUM_e[1]);
  HCC = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  HCL = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  HCR = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (k=0;k<=SpinP_switch; k++) {
    HCC[k] = (double*)malloc(sizeof(double)*NUM_c*NUM_c);
    HCL[k] = (double*)malloc(sizeof(double)*NUM_c*NUM_e[0]);
    HCR[k] = (double*)malloc(sizeof(double)*NUM_c*NUM_e[1]);
  }

  v1 = (dcomplex*) malloc(sizeof(dcomplex)*NUM_c*NUM_c);
  v2 = (dcomplex*) malloc(sizeof(dcomplex)*NUM_c*NUM_c);


  /* initialize */
  TRAN_Set_value_double(SCC,NUM_c*NUM_c,0.0);
  TRAN_Set_value_double(SCL,NUM_c*NUM_e[0],0.0);
  TRAN_Set_value_double(SCR,NUM_c*NUM_e[1],0.0);
  for (k=0;k<=SpinP_switch; k++) {
    TRAN_Set_value_double(HCC[k],NUM_c*NUM_c,0.0);
    TRAN_Set_value_double(HCL[k],NUM_c*NUM_e[0],0.0);
    TRAN_Set_value_double(HCR[k],NUM_c*NUM_e[1],0.0);
  }

#if 0
  TRAN_Connect_Read_Hamiltonian("pol_c_kino.ham_NEGF",
     SpinP_switch,WhatSpecies,FNAN, natn,ncn, Spe_Total_CNO,
       nh,
      CntOLP);
#endif


  /* set CC, CL and CR */
  TRAN_Set_CentOverlap(   mpi_comm_level1, 
                          3,
                          SpinP_switch, 
                          nh, /* input */
                          CntOLP, /* input */
                          atomnum,
                          WhatSpecies,
                          Spe_Total_CNO,
                          FNAN,
                          natn,
                          ncn, 
                          atv_ijk);

#if 0

  TRAN_Read_double("scc",NUM_c, NUM_c,  SCC);
  TRAN_Read_double("scl",NUM_c, NUM_e[0],  SCL);
  TRAN_Read_double("scr",NUM_c, NUM_e[1],  SCR);

  for (k=0;k<=SpinP_switch; k++) {
    TRAN_Read_double("hcc",NUM_c, NUM_c,  HCC[k]);
    TRAN_Read_double("hcl",NUM_c, NUM_e[0],  HCL[k]);
    TRAN_Read_double("hcr",NUM_c, NUM_e[1],  HCR[k]);
  }
#endif

#if 0
  TRAN_Print2_double("scc",NUM_c, NUM_c,  SCC);
  TRAN_Print2_double("scl",NUM_c, NUM_e[0],  SCL);
  TRAN_Print2_double("scr",NUM_c, NUM_e[1],  SCR);

  for (k=0;k<=SpinP_switch; k++) {
    printf("spin=%d of %d\n",k,SpinP_switch);
    TRAN_Print2_double("hcc",NUM_c, NUM_c,  HCC[k]);
    TRAN_Print2_double("hcl",NUM_c, NUM_e[0],  HCL[k]);
    TRAN_Print2_double("hcr",NUM_c, NUM_e[1],  HCR[k]);
  }
#endif



  /* allocate */
  GC = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  GRL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[0]* NUM_e[0]);
  GRR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_e[1]* NUM_e[1]);
  SigmaL = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);
  SigmaR = (dcomplex*)malloc(sizeof(dcomplex)*NUM_c* NUM_c);

  

  if ( tran_transmission_on ) {

    /* initialize */
    for (k=0; k<= SpinP_switch; k++) {
      for (iw=0;iw<tran_transmission_energydiv ; iw++) {

        tran_transmission[k][iw].r = 0.0;
        tran_transmission[k][iw].i = 0.0;

      }
    }

  printf("transmission energy range = [ %le : %le ], imag=%le\n",
          tran_transmission_energyrange[0], 
          tran_transmission_energyrange[1] ,
          tran_transmission_energyrange[2] );


  for (iw=0;iw<tran_transmission_energydiv ; iw++) {

    w.r = tran_transmission_energyrange[0]+
        (tran_transmission_energyrange[1]-tran_transmission_energyrange[0])*
        (double)iw/(tran_transmission_energydiv-1);
    w.i = tran_transmission_energyrange[2];


    printf("iw=%d of %d  w= %le %le \n" ,iw, tran_transmission_energydiv,  w.r,w.i);

    for (k=0; k<= SpinP_switch; k++) {


        iside=0;
        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
             S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max, tran_surfgreen_eps, GRL);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL);


        iside=1;
        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
             S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max, tran_surfgreen_eps, GRR);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR);


        TRAN_Calc_CentGreen(w, NUM_c, SigmaL,SigmaR, HCC[k], SCC, GC);


        TRAN_Calc_OneTransmission(NUM_c, SigmaL, SigmaR, GC,  v1,v2 , &value);

        tran_transmission[k][iw].r = value.r; 
        tran_transmission[k][iw].i = value.i;


    } /* for k */

  } /* iw */


  }





  if ( tran_transmission_iv_on ) {

    /* initialize */
    for (k=0; k<= SpinP_switch; k++) {
      for (iw=0;iw<tran_transmission_iv_energydiv ; iw++) {

        tran_transmission_iv[k][iw].r = 0.0;
        tran_transmission_iv[k][iw].i = 0.0;

      }
    }

  printf("transmission (IV) energy range = [ %le : %le ], imag=%le\n",
          tran_transmission_iv_energyrange[0], 
          tran_transmission_iv_energyrange[1] ,
          tran_transmission_iv_energyrange[2] );


  for (iw=0;iw<tran_transmission_iv_energydiv ; iw++) {

    w.r = tran_transmission_iv_energyrange[0]+
        (tran_transmission_iv_energyrange[1]-tran_transmission_iv_energyrange[0])*
        ((double)iw+0.5)/(tran_transmission_iv_energydiv-1);
    w.i = tran_transmission_iv_energyrange[2];


    printf("iw=%d of %d  w= %le %le \n" ,iw, tran_transmission_iv_energydiv,  w.r,w.i);

    for (k=0; k<= SpinP_switch; k++) {


        iside=0;
        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
             S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max, tran_surfgreen_eps, GRL);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRL, NUM_c, HCL[k], SCL, SigmaL);


        iside=1;
        TRAN_Calc_SurfGreen_direct(w,NUM_e[iside], H00_e[iside][k],H01_e[iside][k],
             S00_e[iside], S01_e[iside], tran_surfgreen_iteration_max, tran_surfgreen_eps, GRR);

        TRAN_Calc_SelfEnergy(w, NUM_e[iside], GRR, NUM_c, HCR[k], SCR, SigmaR);


        TRAN_Calc_CentGreen(w, NUM_c, SigmaL,SigmaR, HCC[k], SCC, GC);


        TRAN_Calc_OneTransmission(NUM_c, SigmaL, SigmaR, GC,  v1,v2 , &value);

        tran_transmission_iv[k][iw].r = value.r; 
        tran_transmission_iv[k][iw].i = value.i;


    } /* for k */

  } /* iw */


  }


  free(SigmaR);
  free(SigmaL);
  free(GRR);
  free(GRL);
  free(GC);
  free(v2);
  free(v1);
  for (k=SpinP_switch;k>=0; k--) {
    free(HCR[k] );
    free(HCL[k] );
    free(HCC[k] );
  }
  free(HCR);
  free(HCL);
  free(HCC);
  free(SCR);
  free(SCL);
  free(SCC);
  free(MP);



}
