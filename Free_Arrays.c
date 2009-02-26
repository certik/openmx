#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

void array0();
void array1();
void array2();
void array3();


void Free_Arrays(int wherefrom)
{

  if      (wherefrom==0) array0();
  else if (wherefrom==1) array1();

}

void array0()
{
  int i,j,k,ii,l,L,n,ct_AN,h_AN,wan,al,tno,Cwan;
  int tno0,tno1,Mc_AN,Gc_AN,Gh_AN,Hwan,m,so,s1,s2;
  int q_AN,Gq_AN,Qwan,Lmax,spe,ns,nc,spin,fan;
  int num,n2,wanA,Gi,Mc_AN_GDC,MAnum;
  int Anum,p,vsize,NUM;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* call from openmx.c */

  /* allocate in truncation.c */

  if (alloc_first[0]==0){

    FNAN[0] = 0; 
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        free(GListTAtoms2[Mc_AN][h_AN]);
        free(GListTAtoms1[Mc_AN][h_AN]);

        if (Allocate_TAtoms0==1){
          free(GListTCells0[Mc_AN][h_AN]);
          free(GListTAtoms3[Mc_AN][h_AN]);
          free(GListTAtoms0[Mc_AN][h_AN]);
	}
      }
      free(GListTAtoms2[Mc_AN]);
      free(GListTAtoms1[Mc_AN]);

      if (Allocate_TAtoms0==1){
        free(GListTCells0[Mc_AN]);
        free(GListTAtoms3[Mc_AN]);
        free(GListTAtoms0[Mc_AN]);
      }
    }
    free(GListTAtoms2);
    free(GListTAtoms1);

    if (Allocate_TAtoms0==1){
      free(GListTCells0);
      free(GListTAtoms3);
      free(GListTAtoms0);
    }
  }

  if (alloc_first[1]==0){
  }

  /* Allocation in UCell_Box() of truncation.c */

  if (alloc_first[2]==0){

    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(MGridListAtom[Mc_AN]);
    }
    free(MGridListAtom);

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      free(GridListAtom[Mc_AN]);
    }
    free(GridListAtom);

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      free(CellListAtom[Mc_AN]);
    }
    free(CellListAtom);
  }

  /* Allocation in truncation.c */
  
  if (alloc_first[3]==0){

    if (SpinP_switch==3){ /* spin non-collinear */

      for (k=0; k<=3; k++){
        free(Density_Grid[k]);
      }
      free(Density_Grid);

    }
    else{ 
      for (k=0; k<=1; k++){
        free(Density_Grid[k]);
      }
      free(Density_Grid);
    }

    free(ADensity_Grid);
    free(PCCDensity_Grid);

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vxc_Grid[k]);
      }
      free(Vxc_Grid);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vxc_Grid[k]);
      }
      free(Vxc_Grid);
    }

    free(VNA_Grid);
    free(dVHart_Grid);

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vpot_Grid[k]);
      }
      free(Vpot_Grid);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vpot_Grid[k]);
      }
      free(Vpot_Grid);
    }

    /* external electric field */
    free(VEF_Grid);

    /* Orbs_Grid */
    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      if (Mc_AN==0){
        tno = 1;
        Gc_AN = 0;
      }
      else{
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno = Spe_Total_NO[Cwan];
      }

      for (i=0; i<tno; i++){
        free(Orbs_Grid[Mc_AN][i]); 
      }
      free(Orbs_Grid[Mc_AN]); 
    }
    free(Orbs_Grid); 

    /* COrbs_Grid */
    if (Cnt_switch!=0){
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
        if (Mc_AN==0){
          tno = 1;
          Gc_AN = 0;
        }
        else{
          Gc_AN = F_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          tno = Spe_Total_CNO[Cwan];
        }

        for (i=0; i<tno; i++){
          free(COrbs_Grid[Mc_AN][i]); 
        }
        free(COrbs_Grid[Mc_AN]); 
      }
      free(COrbs_Grid);
    }
  }

  /****************************************************
    Allocation in truncation.c

    freeing of arrays:

      H0
      CntH0
      HNL
      CntHNL
      OLP
      CntOLP
      H
      CntH
      DS_NL
      CntDS_NL
      DM
      DM_onsite
      NC_OcpN
      NC_v_eff
      ResidualDM
      EDM
      PDM
      IOLP  
      iHNL
      iCntHNL
      H_Hub
  ****************************************************/

  if (alloc_first[4]==0){

    /* H0 */

    for (k=0; k<4; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(H0[k][Mc_AN][h_AN][i]);
	  }

	  free(H0[k][Mc_AN][h_AN]);
	}
	free(H0[k][Mc_AN]);
      }
      free(H0[k]);
    }
    free(H0);

    /* CntH0 */

    if (Cnt_switch==1){
      for (k=0; k<4; k++){
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(CntH0[k][Mc_AN][h_AN][i]);
	    }
	    free(CntH0[k][Mc_AN][h_AN]);
	  }
	  free(CntH0[k][Mc_AN]);
	}
	free(CntH0[k]);
      }
      free(CntH0);
    }

    /* HNL */

    for (k=0; k<List_YOUSO[5]; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
          FNAN[0] = 0;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  for (i=0; i<tno0; i++){
	    free(HNL[k][Mc_AN][h_AN][i]);
	  }
	  free(HNL[k][Mc_AN][h_AN]);
	}
	free(HNL[k][Mc_AN]);
      }
      free(HNL[k]);
    }
    free(HNL);

    /* iHNL */

    if ( SpinP_switch==3 ){

      for (k=0; k<List_YOUSO[5]; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(iHNL[k][Mc_AN][h_AN][i]);
	    }
	    free(iHNL[k][Mc_AN][h_AN]);
	  }
	  free(iHNL[k][Mc_AN]);
	}
	free(iHNL[k]);
      }
      free(iHNL);
    }

    /* iCntHNL */

    if (SO_switch==1 && Cnt_switch==1){

      for (k=0; k<List_YOUSO[5]; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(iCntHNL[k][Mc_AN][h_AN][i]);
	    }
	    free(iCntHNL[k][Mc_AN][h_AN]);
	  }
	  free(iCntHNL[k][Mc_AN]);
	}
	free(iCntHNL[k]);
      }
      free(iCntHNL);
    }

    /* H_Hub  --- added by MJ */  

    if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

      for (k=0; k<=SpinP_switch; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(H_Hub[k][Mc_AN][h_AN][i]);
	    }
            free(H_Hub[k][Mc_AN][h_AN]);
	  }/* h_AN */
          free(H_Hub[k][Mc_AN]);
	}/* Mc_AN */   
        free(H_Hub[k]);
      }/* k */  
      free(H_Hub);
    }

    /* H_Zeeman_NCO */  

    if (Zeeman_NCO_switch==1){

      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (i=0; i<tno0; i++){
	  free(H_Zeeman_NCO[Mc_AN][i]);
	}
        free(H_Zeeman_NCO[Mc_AN]);
      }
      free(H_Zeeman_NCO);
    }

    /* iHNL0 */

    if (SpinP_switch==3){

      for (k=0; k<List_YOUSO[5]; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(iHNL0[k][Mc_AN][h_AN][i]);
	    }
	    free(iHNL0[k][Mc_AN][h_AN]);
	  }
	  free(iHNL0[k][Mc_AN]);
	}
        free(iHNL0[k]);
      }
      free(iHNL0);
    }

    /* OLP_L */  

    for (k=0; k<3; k++){

      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = F_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(OLP_L[k][Mc_AN][h_AN][i]);
	  }
          free(OLP_L[k][Mc_AN][h_AN]);
	}
        free(OLP_L[k][Mc_AN]);
      }
      free(OLP_L[k]);
    }
    free(OLP_L);

    /* OLP */

    for (k=0; k<4; k++){
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  }

	  for (i=0; i<tno0; i++){
	    free(OLP[k][Mc_AN][h_AN][i]);
	  }

	  free(OLP[k][Mc_AN][h_AN]);
	}
	free(OLP[k][Mc_AN]);
      }
      free(OLP[k]);
    }
    free(OLP);

    /* CntOLP */

    if (Cnt_switch==1){

      for (k=0; k<4; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(CntOLP[k][Mc_AN][h_AN][i]);
	    }
	    free(CntOLP[k][Mc_AN][h_AN]);
	  }
	  free(CntOLP[k][Mc_AN]);
	}
	free(CntOLP[k]);
      }
      free(CntOLP);
    }

    /* H */

    for (k=0; k<=SpinP_switch; k++){
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(H[k][Mc_AN][h_AN][i]);
	  }

	  free(H[k][Mc_AN][h_AN]);
	}
	free(H[k][Mc_AN]);
      }
      free(H[k]);
    }
    free(H);

    /* CntH */

    if (Cnt_switch==1){

      for (k=0; k<=SpinP_switch; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_CNO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(CntH[k][Mc_AN][h_AN][i]);
	    }

	    free(CntH[k][Mc_AN][h_AN]);
	  }
	  free(CntH[k][Mc_AN]);
	}
	free(CntH[k]);
      }
      free(CntH);
    }

    /* DS_NL */  

    for (so=0; so<(SO_switch+1); so++){
      for (k=0; k<4; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_VPS_Pro[Hwan] + 2;
	    } 

	    for (i=0; i<tno0; i++){
	      free(DS_NL[so][k][Mc_AN][h_AN][i]);
	    }
	    free(DS_NL[so][k][Mc_AN][h_AN]);
	  }
	  free(DS_NL[so][k][Mc_AN]);
	}
	free(DS_NL[so][k]);
      }
      free(DS_NL[so]);
    }
    free(DS_NL);

    /* CntDS_NL */  

    if (Cnt_switch==1){

      for (so=0; so<(SO_switch+1); so++){
	for (k=0; k<4; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	      FNAN[0] = 0;
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_CNO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno1 = Spe_Total_VPS_Pro[Hwan] + 2;
	      } 

	      for (i=0; i<tno0; i++){
		free(CntDS_NL[so][k][Mc_AN][h_AN][i]);
	      }
	      free(CntDS_NL[so][k][Mc_AN][h_AN]);
	    }
	    free(CntDS_NL[so][k][Mc_AN]);
	  }
	  free(CntDS_NL[so][k]);
	}
	free(CntDS_NL[so]);
      }
      free(CntDS_NL);
    }

    /* DM */  

    for (m=0; m<List_YOUSO[16]; m++){
      for (k=0; k<=SpinP_switch; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
            Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
            Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(DM[m][k][Mc_AN][h_AN][i]);
	    }
            free(DM[m][k][Mc_AN][h_AN]);
	  }
          free(DM[m][k][Mc_AN]);
	}
        free(DM[m][k]);
      }
      free(DM[m]);
    }
    free(DM);

    /* DM_onsite added by MJ */

    if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

      for (m=0; m<2; m++){
	for (k=0; k<=SpinP_switch; k++){

	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
	  
	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }  

	    h_AN = 0;
	      
	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(DM_onsite[m][k][Mc_AN][i]);
	    }

            free(DM_onsite[m][k][Mc_AN]);
	  }
          free(DM_onsite[m][k]);
	}
        free(DM_onsite[m]);
      }
      free(DM_onsite);
    }

    /* v_eff added by MJ */  

    if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){
   
      for (k=0; k<=SpinP_switch; k++){

	FNAN[0] = 0;

        for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }  

	  for (i=0; i<tno0; i++){
	    free(v_eff[k][Mc_AN][i]);
	  }
	  free(v_eff[k][Mc_AN]);

	} /* for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++) */
	free(v_eff[k]);

      } /* for (k=0; k<=SpinP_switch; k++) */
      free(v_eff);

    }

    /*  NC_OcpN */

    if ( (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1)
          && SpinP_switch==3 ){

      for (m=0; m<2; m++){
	for (s1=0; s1<2; s1++){
	  for (s2=0; s2<2; s2++){
	    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	      if (Mc_AN==0){
		Gc_AN = 0;
		tno0 = 1;
	      }
	      else{
		Gc_AN = M2G[Mc_AN];
		Cwan = WhatSpecies[Gc_AN];
		tno0 = Spe_Total_NO[Cwan];  
	      }  

	      for (i=0; i<tno0; i++){
                free(NC_OcpN[m][s1][s2][Mc_AN][i]);
	      }
              free(NC_OcpN[m][s1][s2][Mc_AN]);
	    }
            free(NC_OcpN[m][s1][s2]);
	  }
          free(NC_OcpN[m][s1]);
	}
        free(NC_OcpN[m]);
      }
      free(NC_OcpN);
    }

    /*  NC_v_eff */

    if ( (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) 
         && SpinP_switch==3 ){

      for (s1=0; s1<2; s1++){
	for (s2=0; s2<2; s2++){
	  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }  

	    for (i=0; i<tno0; i++){
	      free(NC_v_eff[s1][s2][Mc_AN][i]);
	    }
            free(NC_v_eff[s1][s2][Mc_AN]);
	  }
          free(NC_v_eff[s1][s2]);
	}
        free(NC_v_eff[s1]);
      }           
      free(NC_v_eff);
    }

    /* ResidualDM */  

    if ( Mixing_switch==0 || Mixing_switch==1 || Mixing_switch==2 ){

      for (m=0; m<List_YOUSO[16]; m++){
	for (k=0; k<=SpinP_switch; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno1 = Spe_Total_NO[Hwan];
	      } 

	      for (i=0; i<tno0; i++){
		free(ResidualDM[m][k][Mc_AN][h_AN][i]);
	      }
	      free(ResidualDM[m][k][Mc_AN][h_AN]);
	    }
	    free(ResidualDM[m][k][Mc_AN]);
	  }
	  free(ResidualDM[m][k]);
	}
	free(ResidualDM[m]);
      }
      free(ResidualDM);

    }

    /* iResidualDM */  

    if ( (Mixing_switch==0 || Mixing_switch==1 || Mixing_switch==2)
	 && SpinP_switch==3 && ( SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1
        || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) ){

      for (m=0; m<List_YOUSO[16]; m++){
        for (k=0; k<2; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno1 = Spe_Total_NO[Hwan];
	      } 

	      for (i=0; i<tno0; i++){
		free(iResidualDM[m][k][Mc_AN][h_AN][i]);
	      }
 	      free(iResidualDM[m][k][Mc_AN][h_AN]);
	    }
            free(iResidualDM[m][k][Mc_AN]);
	  }
          free(iResidualDM[m][k]);
	}
        free(iResidualDM[m]);
      }
      free(iResidualDM);
    }
    else{
      for (m=0; m<List_YOUSO[16]; m++){
        free(iResidualDM[m][0][0][0][0]);
        free(iResidualDM[m][0][0][0]);
        free(iResidualDM[m][0][0]);
        free(iResidualDM[m][0]);
        free(iResidualDM[m]);
      }   
      free(iResidualDM);
    }

    /* EDM */

    for (k=0; k<=SpinP_switch; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(EDM[k][Mc_AN][h_AN][i]);
	  }

	  free(EDM[k][Mc_AN][h_AN]);
	}
	free(EDM[k][Mc_AN]);
      }
      free(EDM[k]);
    }
    free(EDM);

    /* PDM */

    for (k=0; k<=SpinP_switch; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(PDM[k][Mc_AN][h_AN][i]);
	  }

	  free(PDM[k][Mc_AN][h_AN]);
	}
	free(PDM[k][Mc_AN]);
      }
      free(PDM[k]);
    }
    free(PDM);

    /* iDM */

    for (m=0; m<List_YOUSO[16]; m++){
      for (k=0; k<2; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(iDM[m][k][Mc_AN][h_AN][i]);
	    }
	    free(iDM[m][k][Mc_AN][h_AN]);
	  }
	  free(iDM[m][k][Mc_AN]);
	}
	free(iDM[m][k]);
      }
      free(iDM[m]);
    }
    free(iDM);

    /* IOLP */  

    if (Solver==1 || Solver==7){
      FNAN[0] = 0;
      SNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(IOLP[Mc_AN][h_AN][i]);
	  }
	  free(IOLP[Mc_AN][h_AN]);
	}
	free(IOLP[Mc_AN]);
      }
      free(IOLP);
    }

    /* S12 */  

    if (Solver==1 || Solver==5){ /* for recursion, DC and EC */

      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0) n2 = 1;
	else{
	  Gc_AN = M2G[Mc_AN];
	  wan = WhatSpecies[Gc_AN];

	  num = 1;
	  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	    Gi = natn[Gc_AN][i];
	    wanA = WhatSpecies[Gi];
	    num += Spe_Total_CNO[wanA];
	  }
	  n2 = num + 2;
	}

	for (i=0; i<n2; i++){
	  free(S12[Mc_AN][i]);
	}
        free(S12[Mc_AN]);
      }
      free(S12);
    }

    else if (Solver==6){ /* for GDC */

      for (Mc_AN_GDC=0; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

	if (Mc_AN_GDC==0) n2 = 1;
	else{

	  Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
	  Gc_AN = M2G[Mc_AN];
	  wan = WhatSpecies[Gc_AN];

	  num = 1;
	  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	    Gi = natn[Gc_AN][i];
	    wanA = WhatSpecies[Gi];
	    num += Spe_Total_CNO[wanA];
	  }
	  n2 = num + 2;
	}

	for (i=0; i<n2; i++){
	  free(S12[Mc_AN_GDC][i]);
	}
        free(S12[Mc_AN_GDC]);
      }
      free(S12);
    }

    if (Solver==1) { /* for recursion */

      /* Left_U0 */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno1 = 1;
	    vsize = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    wanA = WhatSpecies[Gc_AN];
	    tno1 = Spe_Total_CNO[wanA];

	    Anum = 1;
	    for (p=0; p<=(FNAN[Gc_AN]+SNAN[Gc_AN]); p++){
	      Gi = natn[Gc_AN][p];
	      wanA = WhatSpecies[Gi];
	      Anum += Spe_Total_CNO[wanA];
	    }

	    NUM = Anum - 1;
	    vsize = NUM + 2;
	  }

	  for (k=0; k<List_YOUSO[3]; k++){
	    for (l=0; l<tno1; l++){
	      free(Left_U0[i][j][k][l]);
	    }
            free(Left_U0[i][j][k]);
	  }
          free(Left_U0[i][j]);
	}
        free(Left_U0[i]);
      }
      free(Left_U0);

      /* Right_U0 */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno1 = 1;
	    vsize = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    wanA = WhatSpecies[Gc_AN];
	    tno1 = Spe_Total_CNO[wanA];

	    Anum = 1;
	    for (p=0; p<=(FNAN[Gc_AN]+SNAN[Gc_AN]); p++){
	      Gi = natn[Gc_AN][p];
	      wanA = WhatSpecies[Gi];
	      Anum += Spe_Total_CNO[wanA];
	    }

	    NUM = Anum - 1;
	    vsize = NUM + 2;
	  }

	  for (k=0; k<List_YOUSO[3]; k++){
	    for (l=0; l<tno1; l++){
	      free(Right_U0[i][j][k][l]);
	    }
            free(Right_U0[i][j][k]);
	  }
          free(Right_U0[i][j]);
	}
        free(Right_U0[i]);
      }
      free(Right_U0);
    }

    if (Cnt_switch==1){
      for (i=0; i<=(Matomnum+MatomnumF); i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(CntCoes[i][j]);
	}
	free(CntCoes[i]);
      }
      free(CntCoes);
    }

    if (ProExpn_VNA==1){

      /* HVNA */  

      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(HVNA[Mc_AN][h_AN][i]);
	  }
          free(HVNA[Mc_AN][h_AN]);
	}
        free(HVNA[Mc_AN]);
      }
      free(HVNA);

      /* DS_VNA */  

      for (k=0; k<4; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<(Matomnum+2); Mc_AN++){
          
	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    fan = FNAN[Gc_AN];
	  }
  	  else if ( (Matomnum+1)<=Mc_AN ){
	    fan = List_YOUSO[8];
	    tno0 = List_YOUSO[7];
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	    fan = FNAN[Gc_AN];
	  }
          
	  for (h_AN=0; h_AN<(fan+1); h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      tno1 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	    } 

	    for (i=0; i<tno0; i++){
	      free(DS_VNA[k][Mc_AN][h_AN][i]);
	    }
            free(DS_VNA[k][Mc_AN][h_AN]);
	  }
          free(DS_VNA[k][Mc_AN]);
	}
        free(DS_VNA[k]);
      }
      free(DS_VNA);

      /* CntDS_VNA */  

      if (Cnt_switch==1){

	for (k=0; k<4; k++){

	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<(Matomnum+2); Mc_AN++){
          
	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	      fan = FNAN[Gc_AN];
	    }
	    else if ( Mc_AN==(Matomnum+1) ){
	      fan = List_YOUSO[8];
	      tno0 = List_YOUSO[7];
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_CNO[Cwan]; 
              fan = FNAN[Gc_AN];
	    }

	    for (h_AN=0; h_AN<(fan+1); h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		tno1 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	      } 

	      for (i=0; i<tno0; i++){
		free(CntDS_VNA[k][Mc_AN][h_AN][i]);
	      }
 	      free(CntDS_VNA[k][Mc_AN][h_AN]);
	    }
            free(CntDS_VNA[k][Mc_AN]);
	  }
          free(CntDS_VNA[k]);
	}
        free(CntDS_VNA);
      }

      /* HVNA2 */

      for (k=0; k<4; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(HVNA2[k][Mc_AN][h_AN][i]);
	    }
	    free(HVNA2[k][Mc_AN][h_AN]);
	  }
	  free(HVNA2[k][Mc_AN]);
	}
	free(HVNA2[k]);
      }
      free(HVNA2);

      /* CntHVNA2 */

      if (Cnt_switch==1){

	for (k=0; k<4; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_CNO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	      for (i=0; i<tno0; i++){
		free(CntHVNA2[k][Mc_AN][h_AN][i]);
	      }
	      free(CntHVNA2[k][Mc_AN][h_AN]);
	    }
	    free(CntHVNA2[k][Mc_AN]);
	  }
	  free(CntHVNA2[k]);
	}
	free(CntHVNA2);
      }

    }

    if (Solver==8) { /* embedding cluster method */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

          for (k=0; k<rlmax_EC[j]; k++){
	    for (l=0; l<EKC_core_size[j]; l++){
	      free(Krylov_U[i][j][k][l]);
	    }
            free(Krylov_U[i][j][k]);
	  }
          free(Krylov_U[i][j]);
	}
        free(Krylov_U[i]);
      }
      free(Krylov_U);

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (k=0; k<(rlmax_EC[j]*EKC_core_size[j]+1); k++){
	    free(EC_matrix[i][j][k]);
	  }
          free(EC_matrix[i][j]);
	}
        free(EC_matrix[i]);
      }
      free(EC_matrix);

      free(rlmax_EC);
      free(rlmax_EC2);
      free(EKC_core_size);
      free(scale_rc_EKC);
    }

    /* NEGF */

    if (Solver==4){
      for (spin=0; spin<(SpinP_switch+1); spin++){
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
	  free(TRAN_DecMulP[spin][Mc_AN]);
	}
	free(TRAN_DecMulP[spin]);
      }
      free(TRAN_DecMulP);
    }

  } /*  if (alloc_first[4]==0){ */

  if (alloc_first[5]==0){

    /* NumOLG */
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(NumOLG[Mc_AN]);
    }
    free(NumOLG);

  } /*  if (alloc_first[5]==0){ */

  if (alloc_first[6]==0){

    FNAN[0] = 0;
    SNAN[0] = 0;

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
        free(RMI1[Mc_AN][i]);
      }      
      free(RMI1[Mc_AN]);
    }
    free(RMI1);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
        free(RMI2[Mc_AN][i]);
      }      
      free(RMI2[Mc_AN]);
    }
    free(RMI2);

  } /* if (alloc_first[6]==0){ */

  if (alloc_first[9]==0 && Solver==2){
    for (i=0; i<(Size_Total_Matrix+2); i++){
      free(S[i]);
    }
    free(S);
    free(EV_S);
    free(IEV_S);
  }

  if (alloc_first[10]==0){
    free(M2G);
  }

  if (alloc_first[11]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Snd_MAN[ID]); 
    }
    free(Snd_MAN);

    for (ID=0; ID<numprocs; ID++){
      free(Snd_GAN[ID]); 
    }
    free(Snd_GAN); 
  }

  if (alloc_first[12]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Rcv_GAN[ID]); 
    }
    free(Rcv_GAN); 
  }

  if (alloc_first[13]==0){
    free(F_M2G);
    free(S_M2G);
  }

  if (alloc_first[14]==0){
    free(My_Cell1);
  }

  if (alloc_first[15]==0){
    free(My_Cell0);
    free(Cell_ID0);
    free(edge_block);
  }

  if (alloc_first[16]==0){
    free(Num_Rcv_Grid1);
    free(Num_Snd_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(Rcv_Grid1[ID]);
    }
    free(Rcv_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(Snd_Grid1[ID]);
    }
    free(Snd_Grid1);
  }

  if (alloc_first[17]==0){
    free(Num_IRcv_Grid1);
    free(Num_ISnd_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(IRcv_Grid1[ID]);
    }
    free(IRcv_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(ISnd_Grid1[ID]);
    }
    free(ISnd_Grid1);
  }

  if (alloc_first[18]==0){
    free(Num_Rcv_FNAN2_Grid);
    free(Num_Snd_FNAN2_Grid);
    free(Rcv_FNAN2_MN);
    free(Rcv_FNAN2_GA);
    free(TopMAN2_Grid);

    for (ID=0; ID<numprocs; ID++){
      free(Snd_FNAN2_At[ID]);
    }  
    free(Snd_FNAN2_At);

    for (ID=0; ID<numprocs; ID++){
      free(Snd_FNAN2_Nc[ID]);
    }  
    free(Snd_FNAN2_Nc);
  }

  if (alloc_first[22]==0){

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_GAtom[ID]);
    }
    free(Pro_Snd_GAtom);

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_MAtom[ID]);
    }
    free(Pro_Snd_MAtom);

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_LAtom[ID]);
    }
    free(Pro_Snd_LAtom);

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_LAtom2[ID]);
    }
    free(Pro_Snd_LAtom2);
  }

  /* allocation in Set_BasisPara() of SetPara_DFT.c*/

  for (i=0; i<List_YOUSO[18]; i++){
    for (j=0; j<Spe_Total_NO[i]; j++){
      free(Spe_Trans_Orbital[i][j]);
    }
    free(Spe_Trans_Orbital[i]);
    free(Spe_Specified_Num[i]);
  }
  free(Spe_Trans_Orbital);
  free(Spe_Specified_Num);

  /* Spe_ProductRF_Bessel by Allocate_Arrays(7) in SetPara_DFT.c */

  if (ProExpn_VNA==1){

    for (i=0; i<List_YOUSO[18]; i++){
      for (j=0; j<(Spe_MaxL_Basis[i]+1); j++){
	for (k=0; k<Spe_Num_Basis[i][j]; k++){
	  for (l=0; l<(Spe_MaxL_Basis[i]+1); l++){

	    if (j<=l){
	      Lmax = 2*l;
	    }
	    else{
	      Lmax = 1; 
	    }

	    for (m=0; m<Spe_Num_Basis[i][l]; m++){
	      for (n=0; n<=Lmax; n++){
		free(Spe_ProductRF_Bessel[i][j][k][l][m][n]);
	      }
              free(Spe_ProductRF_Bessel[i][j][k][l][m]);
	    }
            free(Spe_ProductRF_Bessel[i][j][k][l]);
	  }
          free(Spe_ProductRF_Bessel[i][j][k]);
	}
        free(Spe_ProductRF_Bessel[i][j]);
      }
      free(Spe_ProductRF_Bessel[i]);
    }
    free(Spe_ProductRF_Bessel);

  }

  if (alloc_first[23]==0){
    free(NE_T_k_op);
    free(NE_KGrids3);
    free(NE_KGrids2);
    free(NE_KGrids1);
  }

  /* Allocation_Arrays(0) */

  /* arrays for LDA+U added by MJ */
  if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){
    for (i=0; i<SpeciesNum; i++){
      for (l=0; l<(Spe_MaxL_Basis[i]+1); l++){
	free(Hub_U_Basis[i][l]);
      }
      free(Hub_U_Basis[i]);
    }
    free(Hub_U_Basis);

    free(OrbPol_flag);
  }

  for (i=0; i<SpeciesNum; i++){
    free(SpeName[i]);
  }  
  free(SpeName);

  for (i=0; i<SpeciesNum; i++){
    free(SpeBasis[i]);
  }  
  free(SpeBasis);

  for (i=0; i<SpeciesNum; i++){
    free(SpeBasisName[i]);
  }  
  free(SpeBasisName);

  for (i=0; i<SpeciesNum; i++){
    free(SpeVPS[i]);
  }  
  free(SpeVPS);

  free(Spe_MaxL_Basis);

  for (i=0; i<SpeciesNum; i++){
    free(Spe_Num_Basis[i]);
  }  
  free(Spe_Num_Basis);

  for (i=0; i<SpeciesNum; i++){
    free(Spe_Num_CBasis[i]);
  }  
  free(Spe_Num_CBasis);

  free(Spe_Spe2Ban);
  free(Species_Top);
  free(Species_End);
  free(F_Snd_Num);
  free(S_Snd_Num);
  free(F_Rcv_Num);
  free(S_Rcv_Num);
  free(F_Snd_Num_WK);
  free(F_Rcv_Num_WK);
  free(F_TopMAN);
  free(S_TopMAN);
  free(Snd_DS_NL_Size);
  free(Rcv_DS_NL_Size);
  free(Snd_HFS_Size);
  free(Rcv_HFS_Size);
  free(Start_Grid1);
  free(End_Grid1);
  free(Start_Grid2);
  free(End_Grid2);
  free(VPS_j_dependency);

  for (i=0; i<SpeciesNum; i++){
    free(EH0_scaling[i]);
  }
  free(EH0_scaling);

  /* Allocation_Arrays(1) */

  for (i=0; i<(atomnum+1); i++){
    free(Gxyz[i]);
  }
  free(Gxyz);

  num = M_GDIIS_HISTORY + 1;

  for(i=0; i<num; i++) {
    for(j=0; j<(atomnum+1); j++) {
      free(GxyzHistoryIn[i][j]);
    }
    free(GxyzHistoryIn[i]);
  }
  free(GxyzHistoryIn);

  for(i=0; i<num; i++) {
    for(j=0; j<(atomnum+1); j++) {
      free(GxyzHistoryR[i][j]);
    }
    free(GxyzHistoryR[i]);
  }
  free(GxyzHistoryR);

  num = Extrapolated_Charge_History;

  for(i=0; i<num; i++) {
    free(His_Gxyz[i]);
  }
  free(His_Gxyz);

  for(i=0; i<=atomnum; i++){
    free(atom_Fixed_XYZ[i]);
  }
  free(atom_Fixed_XYZ);

  for (i=0; i<(atomnum+1); i++){
    free(Cell_Gxyz[i]);
  }      
  free(Cell_Gxyz);

  free(InitN_USpin);
  free(InitN_DSpin);
  free(WhatSpecies);
  free(GridN_Atom);
  free(RNUM);
  free(RNUM2);
  free(G2ID);
  free(F_G2M);
  free(S_G2M);
  free(time_per_atom);

   /* EF */

  if (MD_switch==4){
    for (i=0; i<(3*atomnum+2); i++){
      free(Hessian[i]);
    }
    free(Hessian);
  }

  /* BFGS */

  if (MD_switch==5){

    for (i=0; i<(3*atomnum+2); i++){
      free(InvHessian[i]);
    }
    free(InvHessian);
  }

  /* RF by hmweng */

  if (MD_switch==6){

    for (i=0; i<(3*atomnum+2); i++){
      free(Hessian[i]);
    }
    free(Hessian);
  }

  /* spin non-collinear */
  if (SpinP_switch==3){
    free(Angle0_Spin);
    free(Angle1_Spin);
    free(InitAngle0_Spin);
    free(InitAngle1_Spin);
    free(Angle0_Orbital);
    free(Angle1_Orbital);
    free(OrbitalMoment);
    free(Constraint_SpinAngle);

    free(InitAngle0_Orbital);
    free(InitAngle1_Orbital);
    for(i=0; i<(atomnum+1); i++){
      free(Orbital_Moment_XYZ[i]);
    }
    free(Orbital_Moment_XYZ);
    free(Constraint_OrbitalAngle);
  }
   /* Allocation_Arrays(2) */

  free(NormK);
  free(Spe_Atom_Cut1);
  free(Spe_Core_Charge);
  free(TGN_EH0);
  free(dv_EH0);
  free(Spe_Num_Mesh_VPS);
  free(Spe_Num_Mesh_PAO);
  free(Spe_Total_VPS_Pro);
  free(Spe_Num_RVPS);
  free(Spe_PAO_LMAX);
  free(Spe_PAO_Mul);
  free(Spe_WhatAtom);
  free(Spe_Total_NO);
  free(Spe_Total_CNO);
  free(FNAN);
  free(SNAN);
  free(SNAN_GDC);
  free(True_SNAN);
  free(zp);
  free(Ep);
  free(Rp);

  /* Allocation_Arrays(3) */

  if (alloc_first[8]==0){

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      free(natn[ct_AN]);
    }
    free(natn);

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      free(ncn[ct_AN]);
    }
    free(ncn);

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      free(Dis[ct_AN]);
    }
    free(Dis);
  }

  /* Allocation_Arrays(4) */

  if (remake_headfile==0){

    for (i=0; i<SpeciesNum; i++){
      free(GridX_EH0[i]);
    }
    free(GridX_EH0);

    for (i=0; i<SpeciesNum; i++){
      free(GridY_EH0[i]);
    }
    free(GridY_EH0);

    for (i=0; i<SpeciesNum; i++){
      free(GridZ_EH0[i]);
    }
    free(GridZ_EH0);

    for (i=0; i<SpeciesNum; i++){
      free(Arho_EH0[i]);
    }
    free(Arho_EH0);

    for (i=0; i<SpeciesNum; i++){
      free(Wt_EH0[i]);
    }
    free(Wt_EH0);

  }

  /* Allocation_Arrays(5) */

  for (i=0; i<(MO_Nkpoint+1); i++){
    free(MO_kpoint[i]);
  }
  free(MO_kpoint);

  /* Set_Periodic() in truncation.c */

  n = 2*CpyCell + 4;
  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      free(ratv[i][j]);
    }
    free(ratv[i]);
  }
  free(ratv);

  for (i=0; i<(TCpyCell+1); i++){
    free(atv[i]);
  }
  free(atv);

  for (i=0; i<(TCpyCell+1); i++){
    free(atv_ijk[i]);
  }
  free(atv_ijk);

  /* Allocate_Arrays(6) in SetPara_DFT.c */

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_PAO_XV[i]);
  }
  free(Spe_PAO_XV);
  
  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_PAO_RV[i]);
  }
  free(Spe_PAO_RV);

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_Atomic_Den[i]);
  }
  free(Spe_Atomic_Den);

  for (i=0; i<List_YOUSO[18]; i++){
    for (j=0; j<=List_YOUSO[25]; j++){
      for (k=0; k<List_YOUSO[24]; k++){
        free(Spe_PAO_RWF[i][j][k]);
      }
      free(Spe_PAO_RWF[i][j]);
    }
    free(Spe_PAO_RWF[i]);
  }
  free(Spe_PAO_RWF);

  for (i=0; i<List_YOUSO[18]; i++){
    for (j=0; j<=List_YOUSO[25]; j++){
      for (k=0; k<List_YOUSO[24]; k++){
        free(Spe_RF_Bessel[i][j][k]);
      }
      free(Spe_RF_Bessel[i][j]);
    }
    free(Spe_RF_Bessel[i]);
  }
  free(Spe_RF_Bessel);

  /* Allocate_Arrays(7) in SetPara_DFT.c */

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_VPS_XV[i]);
  }
  free(Spe_VPS_XV);

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_VPS_RV[i]);
  }
  free(Spe_VPS_RV);

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_Vna[i]);
  }
  free(Spe_Vna);

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_VH_Atom[i]);
  }
  free(Spe_VH_Atom);

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_Atomic_PCC[i]);
  }
  free(Spe_Atomic_PCC);

  for (so=0; so<(SO_switch+1); so++){
    for (i=0; i<List_YOUSO[18]; i++){
      for (j=0; j<List_YOUSO[19]; j++){
        free(Spe_VNL[so][i][j]);
      }
      free(Spe_VNL[so][i]);
    }
    free(Spe_VNL[so]);
  }
  free(Spe_VNL);

  for (so=0; so<(SO_switch+1); so++){
    for (i=0; i<List_YOUSO[18]; i++){
      free(Spe_VNLE[so][i]);
    }
    free(Spe_VNLE[so]);
  }
  free(Spe_VNLE);

  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_VPS_List[i]);
  }
  free(Spe_VPS_List);

  for (so=0; so<(SO_switch+1); so++){
    for (i=0; i<List_YOUSO[18]; i++){
      for (j=0; j<(List_YOUSO[19]+2); j++){
        free(Spe_NLRF_Bessel[so][i][j]);
      }
      free(Spe_NLRF_Bessel[so][i]);
    }
    free(Spe_NLRF_Bessel[so]);
  }
  free(Spe_NLRF_Bessel);

  if (ProExpn_VNA==1){

    for (i=0; i<List_YOUSO[18]; i++){
      for (L=0; L<(List_YOUSO[35]+1); L++){
	for (j=0; j<List_YOUSO[34]; j++){
	  free(Projector_VNA[i][L][j]);
	}
        free(Projector_VNA[i][L]);
      }
      free(Projector_VNA[i]);
    }
    free(Projector_VNA);

    for (i=0; i<List_YOUSO[18]; i++){
      for (L=0; L<(List_YOUSO[35]+1); L++){
	free(VNA_proj_ene[i][L]);
      }
      free(VNA_proj_ene[i]);
    }
    free(VNA_proj_ene);

    for (i=0; i<List_YOUSO[18]; i++){
      for (L=0; L<(List_YOUSO[35]+1); L++){
	for (j=0; j<List_YOUSO[34]; j++){
	  free(Spe_VNA_Bessel[i][L][j]);
	}
        free(Spe_VNA_Bessel[i][L]);
      }
      free(Spe_VNA_Bessel[i]);
    }
    free(Spe_VNA_Bessel);

    for (i=0; i<List_YOUSO[18]; i++){
      free(Spe_CrudeVNA_Bessel[i]);
    }
    free(Spe_CrudeVNA_Bessel);
  }

  for (i=0; i<List_YOUSO[33]; i++){
    for (j=0; j<2; j++){
      for (k=0; k<List_YOUSO[31]; k++){
        for (l=0; l<List_YOUSO[1]; l++){
          free(HOMOs_Coef[i][j][k][l]);
	}
        free(HOMOs_Coef[i][j][k]);
      }
      free(HOMOs_Coef[i][j]);
    }
    free(HOMOs_Coef[i]);
  }
  free(HOMOs_Coef);

  for (i=0; i<List_YOUSO[33]; i++){
    for (j=0; j<2; j++){
      for (k=0; k<List_YOUSO[32]; k++){
        for (l=0; l<List_YOUSO[1]; l++){
          free(LUMOs_Coef[i][j][k][l]);
	}
        free(LUMOs_Coef[i][j][k]);
      }
      free(LUMOs_Coef[i][j]);
    }
    free(LUMOs_Coef[i]);
  }
  free(LUMOs_Coef);

  free(Bulk_Num_HOMOs);
  free(Bulk_Num_LUMOs);
  for (i=0; i<List_YOUSO[33]; i++){
    free(Bulk_HOMO[i]);
  }
  free(Bulk_HOMO);

  if (CntOrb_fileout==1){
    free(CntOrb_Atoms);
  }

  /* allocated in Input_std.c */

  if (Band_Nkpath>0) {

    free(Band_N_perpath);

    for (i=0; i<(Band_Nkpath+1); i++){
      for (j=0; j<3; j++){
	free(Band_kpath[i][j]);
      }
      free(Band_kpath[i]);
    }
    free(Band_kpath);

    for (i=0; i<(Band_Nkpath+1); i++){
      for (j=0; j<3; j++){
	free(Band_kname[i][j]);
      }
      free(Band_kname[i]);
    }
    free(Band_kname);
  }

}

void array1()
{


}

