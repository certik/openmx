#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Inputtools.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"

#define MAXBUF  256




void TRAN_Input_std(
  MPI_Comm comm1, 
  int Solver,          /* input */
  int SpinP_switch,  
  char *filepath,
  double kBvalue,
  double Electronic_Temperature
                      /* no output */
)
{
  FILE *fp;
  int i;
  int i_vec[20],i_vec2[20];
  double r_vec[20];
  char *s_vec[20];
  char buf[MAXBUF];
  int po=0;

  /****************************************************
               parameters for TRANSPORT
  ****************************************************/

  input_logical("Tran.Output_HKS",&TRAN_output_hks,0);
  /* printf("Tran.OutputHKS=%d\n",TRAN_output_hks); */
  input_string("Tran.filename.HKS",TRAN_hksoutfilename,"tran.hks");
  /* printf("TRAN_hksoutfilename=%s\n",TRAN_hksoutfilename); */

  if ( Solver!=4 ) { return; }

  /**** show transport credit ****/
  TRAN_Credit(comm1);

  input_string("Tran.filename.hks.l",TRAN_hksfilename[0],"tran.hks.l");
  input_string("Tran.filename.hks.r",TRAN_hksfilename[1],"tran.hks.r");

  TRAN_RestartFile(comm1, "read","left", filepath,TRAN_hksfilename[0]);
  TRAN_RestartFile(comm1, "read","right",filepath,TRAN_hksfilename[1]);

  /* change chemical potentials */

  i=0;
  r_vec[i++]=ChemP_e[0];
  r_vec[i++]=ChemP_e[1];
  input_doublev("Tran.ChemP",i,ChemP_e,r_vec);

  /* check the conflict of SpinP_switch */

  if ( (SpinP_switch!=SpinP_switch_e[0]) || (SpinP_switch!=SpinP_switch_e[1]) ){
    printf ("scf.SpinPolarization conflicts between leads or lead and center.\n");
    MPI_Finalize();
    exit(0);
  }

  input_double("Tran.temperature",&tran_temperature,300.0); /* unit Kelvin */
  tran_temperature *= kBvalue; 

  input_int(   "Tran.Surfgreen.iterationmax", &tran_surfgreen_iteration_max, 300);
  input_double("Tran.Surfgreen.convergeeps", &tran_surfgreen_eps, 1.0e-12);

  /****  k-points parallel to the layer, which are used for the SCF calc. ****/
  
  i_vec2[0]=1;
  i_vec2[1]=1;
  input_intv("Tran.scf.Kgrid",2,i_vec,i_vec2);
  TRAN_Kspace_grid2 = i_vec[0];
  TRAN_Kspace_grid3 = i_vec[1];

  if (TRAN_Kspace_grid2<=0){
    printf("Tran.scf.Kgrid should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 
  if (TRAN_Kspace_grid3<=0){
    printf("Tran.scf.Kgrid should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 

  /****  k-points parallel to the layer, which are used for the transmission calc. ****/
  
  i_vec2[0]=1;
  i_vec2[1]=1;
  input_intv("Tran.tran.Kgrid",2,i_vec,i_vec2);
  TRAN_TKspace_grid2 = i_vec[0];
  TRAN_TKspace_grid3 = i_vec[1];

  if (TRAN_TKspace_grid2<=0){
    printf("Tran.tran.Kgrid should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 
  if (TRAN_TKspace_grid3<=0){
    printf("Tran.tran.Kgrid should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 

  /**** bias ****/

  input_logical("Tran.bias.apply",&tran_bias_apply,0);
  if ( tran_bias_apply ) {
    i=0;
    r_vec[i++]=0.0;
    r_vec[i++]=0.0;
    input_doublev("Tran.bias.voltage",i,tran_biasvoltage_e,r_vec);
  }
  else {
    tran_biasvoltage_e[0]=0.0;
    tran_biasvoltage_e[1]=0.0;
  }

  if (tran_bias_apply) {
    int side;
    side=0;
    TRAN_Apply_Bias2e(comm1,  tran_biasvoltage_e[side], 
		      SpinP_switch_e[side], atomnum_e[side],
		      WhatSpecies_e[side], Spe_Total_CNO_e[side], FNAN_e[side], natn_e[side],
		      Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side], OLP_e[side][0],
		      &ChemP_e[side],H_e[side], dVHart_Grid_e[side] ); /* output */
    side=1;
    TRAN_Apply_Bias2e(comm1,  tran_biasvoltage_e[side], 
		      SpinP_switch_e[side], atomnum_e[side],
		      WhatSpecies_e[side], Spe_Total_CNO_e[side], FNAN_e[side], natn_e[side],
		      Ngrid1_e[side], Ngrid2_e[side], Ngrid3_e[side], OLP_e[side][0],
		      &ChemP_e[side], H_e[side], dVHart_Grid_e[side] ); /* output */
  }

  /******************************************************
            parameters for the DOS calculation         
  ******************************************************/
  
  i=0;
  r_vec[i++] = -2.0;
  r_vec[i++] =  2.0;
  r_vec[i++] = 1.0e-6;
  input_doublev("Tran.Dos.energyrange",i, tran_dos_energyrange, r_vec);
  input_int("Tran.Dos.energydiv",&tran_dos_energydiv,100);
  
  i_vec2[0]=1;
  i_vec2[1]=1;
  input_intv("Tran.Dos.Kgrid",2,i_vec,i_vec2);
  TRAN_dos_Kspace_grid2 = i_vec[0];
  TRAN_dos_Kspace_grid3 = i_vec[1];

  /******************************************************
                path type for the integration
  ******************************************************/

  s_vec[0]="square"; s_vec[1]="thermalarc"; s_vec[2]="line"; s_vec[3]="GaussHG"; 
  i_vec[0]=1;        i_vec[1]=10;           i_vec[2]= 2;     i_vec[3]= 3;
  input_string2int("Tran.Integ.Pathtype", &tran_integ_pathtype, 4, s_vec,i_vec);

#if 1
  if  ( tran_integ_pathtype==1 ) {  /* square */
    /*   double default_absolute_ene[4]; */
    double default_relative_ene[4];
    char **str,**defaultstr;
    int i,n;

    n=6;
    str =(char**)malloc(sizeof(char*)*n);
    defaultstr=(char**)malloc(sizeof(char*)*n);
    for (i=0;i<6;i++) {
      str[i]=(char*)malloc(sizeof(char)*50);
    }
    defaultstr[0]="auto";defaultstr[1]="r";
    defaultstr[2]="auto";defaultstr[3]="r";
    defaultstr[4]="auto";
    defaultstr[5]="auto";
    input_stringv("Tran.squarePath.energies", 6,str,defaultstr);
    /* format:    -1.0 a auto a 0.0 0.1 , auto|numeric */
    /* default_absolute_ene[0]=0.0;
     * default_absolute_ene[1]=(ChemP_e[0]>=ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
     * default_absolute_ene[2]=1.0e-6;
     * default_absolute_ene[3]=0.1; 
     */
    default_relative_ene[0]=-0.15; 
    default_relative_ene[1]=0.0;
    default_relative_ene[2]=1.0e-6;
    default_relative_ene[3]=0.1;
    TRAN_Set_PathEnergyStr(tran_integ_pathtype, 4, str, 
			   default_relative_ene, 
			   tran_square_path_ene, tran_square_path_ene_fix );
    if ( tran_square_path_ene_fix[1]==0 ) {
      tran_square_path_ene[1]+= (ChemP_e[0]>=ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
      tran_square_path_ene_fix[1]=1;
    }
    /* now tran_square_path_ene_fix[0] maybe 0, then one must add minE */
         

    for (i=0;i<n;i++) {
      free(str[i]);
    }
    free(str);
    free(defaultstr);

#endif
       
    for (i=0;i<3;i++) i_vec[i]=0;
    input_intv("Tran.squarePath.div",3,tran_square_path_div,i_vec);

    input_double("Tran.squarePath.bias.ExpandEnergy",&tran_square_path_bias_expandenergy,0.01);
    input_int("Tran.squarePath.bias.div", & tran_square_path_bias_div,10);

    if (tran_bias_apply) {
      tran_square_path_ene[1] -= fabs(ChemP_e[0]-ChemP_e[1]); 
      /* now tran_square_path_ene[1] = smaller one between ChemP_e[0] and ChemP_e[1] */
      tran_square_path_ene[1] -= tran_square_path_bias_expandenergy; 
      printf("bias=ON, tran_square_path_ene[1]=%le\n",tran_square_path_ene[1]);
    }

    {  
      printf("------------------------------------------------------\n");
      printf("tran_square_path_ene_fix tran_square_path_ene\n");
      for (i=0;i<4;i++) {
	printf("%d %d %le\n", i, tran_square_path_ene_fix[i], tran_square_path_ene[i]);
      }
      printf("------------------------------------------------------\n");

    }

    /* ---------------------------------------------------------------------------------- 
     *  (0,3) ->  (1,3)
     *    A         |
     *    |         V
     *  (0,2)     (1,2)  ->   ( 1+|ChemP_e[0]-ChemP_e[1]|+ tran_square_path_bias_ene,2 ) 
     * ----------------------------------------------------------------------------------
     *  note that 
     *     0 does not contain the effect of 'the minimum eigen energy' 
     *                               in the case of tran_integ_relativeenergy=1
     * ----------------------------------------------------------------------------------
     */


  } /* tran_integ_pathtype == 1 */

  else if ( tran_integ_pathtype==2 ) {  /* line */
    input_int("Tran.linepath.div",&tran_line_path_div,0);
    /* allocate */
    tran_line_path_string = (char**)malloc(sizeof(char*)*(tran_line_path_div*2+(tran_line_path_div-1)*2));
    for (i=0;i< tran_line_path_div*2+(tran_line_path_div-1)*2; i++) {
      tran_line_path_string[i]=(char*)malloc(sizeof(char)*YOUSO10);
    }

    if (fp=input_find("<Tran.linepath") ) {

      for (i=0; i<tran_line_path_div; i++){
	fgets(buf,MAXBUF,fp);
	printf("%d->%s",i,buf);
	if (i==tran_line_path_div-1) {
	  sscanf(buf,"%s %s", 
		 tran_line_path_string[4*i],tran_line_path_string[4*i+1]);
	} 
	else {
	  sscanf(buf,"%s %s %s %s", 
		 tran_line_path_string[4*i],tran_line_path_string[4*i+1],
		 tran_line_path_string[4*i+2],tran_line_path_string[4*i+3]);
	}
      }

      ungetc('\n',fp);
      if ( ! input_last("Tran.linepath>") ) {
	/* format error */
	printf("Format error for Tran.linepath\n");
	po++;
      }

    }
  } /* tran_integ_pathtype==2 */

  else if (tran_integ_pathtype==10) { /* thermal arc path */

/*  
*       ----------
*     /            \ (1-gamma,3)  (1,3)
*    /              -------------+--------(1+5,3) ,   5 = tran_temperature*15 
*    |                           x
*    |                           x
*    |                           x
*    ---+------------------------------>
*    (0,2)
*
*   n = int_C A(x) nF(x) dx 
*
*  (0,2) -> (1-gamma,3)   divided by n0
*  (1-gamma,3)->(1,3)     divided by n1
*  (1,3) -> (1+5,3)       divided by n2 
*
*   if (applied_bias)     additional points
*/

    double default_relative_ene[6];
    char **str,**defaultstr;
    int i,n;

    n=8;
    str =(char**)malloc(sizeof(char*)*n);
    defaultstr=(char**)malloc(sizeof(char*)*n);
    for (i=0;i<7;i++) {
      str[i]=(char*)malloc(sizeof(char)*50);
    }
    defaultstr[0]="auto";defaultstr[1]="r";
    defaultstr[2]="auto";defaultstr[3]="r";
    defaultstr[4]="auto";
    defaultstr[5]="auto";
    defaultstr[6]="auto";

    input_stringv("Tran.thermalarcPath.energies", n,str,defaultstr);
    /* format:    -1.0 a auto a 0.0 0.1 0.1 , auto|numeric */

    default_relative_ene[0]=-0.15; 
    default_relative_ene[1]=tran_temperature*10.0; /* extension to higher energy */
    default_relative_ene[2]=1.0e-6;
    default_relative_ene[3]=0.1;
    default_relative_ene[4]=tran_temperature*10.0+.05; /* gamma */

    TRAN_Set_PathEnergyStr(tran_integ_pathtype,5, str, 
			   default_relative_ene, 
			   tran_thermalarc_path_ene, tran_thermalarc_path_ene_fix );

    if ( tran_thermalarc_path_ene_fix[1]==0 ) {
      tran_thermalarc_path_ene[1]+= (ChemP_e[0]>=ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
      tran_thermalarc_path_ene_fix[1]=1;
    }
    /* now tran_thermal_path_ene_fix[0] maybe 0, then one must add minE */

    for (i=0;i<n;i++) {
      free(str[i]);
    }
    free(str);
    free(defaultstr);

    for (i=0;i<3;i++) i_vec[i]=0;
    input_intv("Tran.thermalarcPath.div",3,tran_thermalarc_path_div,i_vec);

    input_double("Tran.thermalarcPath.bias.ExpandEnergy",&tran_thermalarc_path_bias_expandenergy,0.01);
    input_int("Tran.thermalarcPath.bias.div", & tran_thermalarc_path_bias_div,10);

    printf("thermalarcPath\n");
    for (i=0;i<5;i++) {
      printf("%d %lf\n", tran_thermalarc_path_ene_fix[i], tran_thermalarc_path_ene[i]);
    }
      

    if (tran_bias_apply) {
      tran_thermalarc_path_ene[1] = (ChemP_e[0]<ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
      /* now ene[1] = smaller one between ChemP_e[0] and ChemP_e[1] */
      tran_thermalarc_path_ene[1] -=  tran_thermalarc_path_bias_expandenergy; 

      printf("bias=ON, tran_thermalarc_path_ene[1]=%le\n",
	     tran_thermalarc_path_ene[1]);
    }


  } /* tran_integ_pathtype==10 */

  /* Taisuke Ozaki Copyright (C) */
  /* GaussHG, Gauss's hypergeometric */

  else if ( tran_integ_pathtype==3 ) {  
    input_int("Tran.Num.Poles", &tran_num_poles,150);
  }

  TRAN_Set_IntegPath( kBvalue, Electronic_Temperature );
 
  input_logical("Tran.transmission.on",&tran_transmission_on,1);
  if (tran_transmission_on) {
    /* calculate transmission at (tran_transmission_energyrange[0]+i tran_transmission_energyrange[2],
     *                 tran_transmission_energyrange[1]+i tran_transmission_energyrange[2] )
     */
    i=0;
    r_vec[i++] = -2.0;
    r_vec[i++] =  2.0;
    r_vec[i++] = 1.0e-6;
    input_doublev("Tran.transmission.energyrange",i, tran_transmission_energyrange, r_vec);
    input_int("Tran.transmission.energydiv",&tran_transmission_energydiv,100);

    /*
     *       tran_transmission = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1));
     *       for (i=0;i<= SpinP_switch; i++) {
     *	 tran_transmission[i]=(dcomplex*)malloc(sizeof(dcomplex)*tran_transmission_energydiv);
     *       }
     */
  }
  else {
    tran_transmission_energydiv=0;
    tran_transmission=NULL;
  }

  input_logical("Tran.transmission.IV.on",&tran_transmission_iv_on,
		(tran_bias_apply)? 1:0 );

  if (tran_transmission_iv_on) {
       
    i=0;
    r_vec[i++] = (ChemP_e[0]<ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
    r_vec[i++] = (ChemP_e[0]>ChemP_e[1])? ChemP_e[0]:ChemP_e[1];
    r_vec[i++] = 1.0e-6;
    input_doublev("Tran.transmission.IV.energyrange",i, tran_transmission_iv_energyrange, r_vec);
    input_int("Tran.transmission.IV.energydiv",&tran_transmission_iv_energydiv,15);
       
    /*
     *       tran_transmission_iv = (dcomplex**)malloc(sizeof(dcomplex*)*(SpinP_switch+1));
     *       for (i=0;i<= SpinP_switch; i++) {
     *	 tran_transmission_iv[i]=(dcomplex*)malloc(sizeof(dcomplex)*tran_transmission_iv_energydiv);
     *       }
     */

  }
  else {
    tran_transmission_iv_energydiv=0;
    tran_transmission_iv=NULL;
  }


#if 0
  { int  side;
  side=0;
  TRAN_Connect_Read_Hamiltonian("au3.ham_NEGF",
				SpinP_switch_e[side], WhatSpecies_e[side],FNAN_e[side],
				natn_e[side], ncn_e[side], Spe_Total_CNO_e[side], 
				H_e[side], OLP_e[side][0] );
  side=1;
  TRAN_Connect_Read_Hamiltonian("au3.ham_NEGF",
				SpinP_switch_e[side], WhatSpecies_e[side],FNAN_e[side],
				natn_e[side], ncn_e[side], Spe_Total_CNO_e[side], 
				H_e[side], OLP_e[side][0] );
  }
#endif
     
#if 0
  { int side,k;
  side=0;
  TRAN_FPrint2_double("zs00l",NUM_e[side],NUM_e[side],S00_e[side]);
  TRAN_FPrint2_double("zs01l",NUM_e[side],NUM_e[side],S01_e[side]);
  k=0;
  TRAN_FPrint2_double("zh00l",NUM_e[side],NUM_e[side],H00_e[side][k]);
  TRAN_FPrint2_double("zh01l",NUM_e[side],NUM_e[side],H01_e[side][k]);
  side=1;
  TRAN_FPrint2_double("zs00r",NUM_e[side],NUM_e[side],S00_e[side]);
  TRAN_FPrint2_double("zs01r",NUM_e[side],NUM_e[side],S01_e[side]);
  k=0;
  TRAN_FPrint2_double("zh00r",NUM_e[side],NUM_e[side],H00_e[side][k]);
  TRAN_FPrint2_double("zh01r",NUM_e[side],NUM_e[side],H01_e[side][k]);
  printf("stop after writing h00\n");
  exit(0);
  }
#endif

#if 0
  {int side,k;
  side=0;
  TRAN_Read_double("s00l",NUM_e[side],NUM_e[side],S00_e[side]);
  TRAN_Read_double("s01l",NUM_e[side],NUM_e[side],S01_e[side]);
  for (k=0;k<=SpinP_switch_e[side];k++) {
    TRAN_Read_double("h00l",NUM_e[side],NUM_e[side],H00_e[side][k]);
    TRAN_Read_double("h01l",NUM_e[side],NUM_e[side],H01_e[side][k]);
  }
  side=1;
  TRAN_Read_double("s00r",NUM_e[side],NUM_e[side],S00_e[side]);
  TRAN_Read_double("s01r",NUM_e[side],NUM_e[side],S01_e[side]);
  for (k=0;k<=SpinP_switch_e[side];k++) {
    TRAN_Read_double("h00r",NUM_e[side],NUM_e[side],H00_e[side][k]);
    TRAN_Read_double("h01r",NUM_e[side],NUM_e[side],H01_e[side][k]);
  }
  }
#endif

}



