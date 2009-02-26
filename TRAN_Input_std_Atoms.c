#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Inputtools.h"

#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif


#include "tran_prototypes.h"
#include "tran_variables.h"

void TRAN_Input_std_Atoms( int Solver )
{
  int po=0;

  int Catomnum,Latomnum,Ratomnum;
  FILE *fp;

  char *s_vec[20];
  int i_vec[20];
  double r_vec[20];

  int i,j; 
  char Species[YOUSO10];
  double Length_C, Length_L, Length_R;
  double angleCL, angleCR;
  double Lsign, Rsign; 

 if (Solver!=4) return; 

    /* center */
    input_int("Atoms.Number",&Catomnum,0);
    if (Catomnum<=0){
      printf("Atoms.Number may be wrong.\n");
      po++;
    }

    /* left */
    input_int("LeftLeadAtoms.Number",&Latomnum,0);
    if (Latomnum<=0){
      printf("LeftLeadAtoms.Number may be wrong.\n");
      po++;
    }

    /* right */
    input_int("RightLeadAtoms.Number",&Ratomnum,0);
    if (Ratomnum<=0){
      printf("RightLeadAtoms.Number may be wrong.\n");
      po++;
    }
    
    atomnum = Catomnum + Latomnum + Ratomnum;
    List_YOUSO[1] = atomnum + 1;

    /* memory allocation */

    Allocate_Arrays(1);

    /* memory allocation for TRAN_* */
    TRAN_Allocate_Atoms(atomnum);


    s_vec[0]="Ang";  s_vec[1]="AU";
    i_vec[0]= 0;     i_vec[1]=1;
    input_string2int("Atoms.SpeciesAndCoordinates.Unit",
                     &coordinates_unit,2,s_vec,i_vec);

#ifdef TRAN
    /* left */
    if (fp=input_find("<LeftLeadAtoms.SpeciesAndCoordinates") ) {
      for (i=1; i<=Latomnum; i++){
        fscanf(fp,"%i %s %lf %lf %lf %lf %lf",&j,Species,
               &Gxyz[i][1],&Gxyz[i][2],&Gxyz[i][3],
               &InitN_USpin[i],&InitN_DSpin[i]);
        WhatSpecies[i] = Species2int(Species);
        TRAN_region[i]=2;
        TRAN_Original_Id[i]=j;

        if (i!=j){
          printf("Error of sequential number %i in <LeftLeadAtoms.SpeciesAndCoordinates\n",j);
          po++;
        }

        if (2<=level_stdout){
          printf("<Input_std> L_AN=%2d T_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
                 i,i,
                  WhatSpecies[i],
                  InitN_USpin[i],
                  InitN_DSpin[i]);
        }
      }

      if (!input_last("LeftLeadAtoms.SpeciesAndCoordinates>")) {
        /* format error */
        printf("Format error for LeftLeadAtoms.SpeciesAndCoordinates\n");
        po++;
      }
    }
#endif

    /* center */
    if (fp=input_find("<Atoms.SpeciesAndCoordinates") ) {
      for (i=1; i<=Catomnum; i++){
        fscanf(fp,"%i %s %lf %lf %lf %lf %lf",&j,Species,
#ifdef TRAN
               &Gxyz[Latomnum+i][1],&Gxyz[Latomnum+i][2],&Gxyz[Latomnum+i][3],&InitN_USpin[Latomnum+i],&InitN_DSpin[Latomnum+i]);
        WhatSpecies[Latomnum+i] = Species2int(Species);
        TRAN_region[Latomnum+i]= 1;
        TRAN_Original_Id[Latomnum+i]= j;
#else
	       &Gxyz[i][1],&Gxyz[i][2],&Gxyz[i][3],&InitN_USpin[i],&InitN_DSpin[i]);
        WhatSpecies[i] = Species2int(Species);
#endif

        if (i!=j){
          printf("Error of sequential number %i in <Atoms.SpeciesAndCoordinates\n",j);
          po++;
        }

        if (2<=level_stdout){
          printf("<Input_std>  ct_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
#ifdef TRAN
                  Latomnum+i,WhatSpecies[Latomnum+i],InitN_USpin[Latomnum+i],InitN_DSpin[Latomnum+i]);

#else
                  i,WhatSpecies[i],InitN_USpin[i],InitN_DSpin[i]);
#endif
        }
      }

      if (!input_last("Atoms.SpeciesAndCoordinates>")) {
        /* format error */
        printf("Format error for Atoms.SpeciesAndCoordinates\n");
        po++;
      }
    }

#ifndef TRAN
    /* left */
    if (fp=input_find("<LeftLeadAtoms.SpeciesAndCoordinates") ) {
      for (i=1; i<=Latomnum; i++){
        fscanf(fp,"%i %s %lf %lf %lf %lf %lf",&j,Species,
	       &Gxyz[Catomnum+i][1],&Gxyz[Catomnum+i][2],&Gxyz[Catomnum+i][3],
               &InitN_USpin[Catomnum+i],&InitN_DSpin[Catomnum+i]);
        WhatSpecies[Catomnum+i] = Species2int(Species);

        if (i!=j){
          printf("Error of sequential number %i in <LeftLeadAtoms.SpeciesAndCoordinates\n",j);
          po++;
        }

        if (2<=level_stdout){
          printf("<Input_std> L_AN=%2d T_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
		 i,Catomnum+i,
                  WhatSpecies[Catomnum+i],
                  InitN_USpin[Catomnum+i],
                  InitN_DSpin[Catomnum+i]);
        }
      }

      if (!input_last("LeftLeadAtoms.SpeciesAndCoordinates>")) {
        /* format error */
        printf("Format error for LeftLeadAtoms.SpeciesAndCoordinates\n");
        po++;
      }
    }
#endif


    /* right */
    if (fp=input_find("<RightLeadAtoms.SpeciesAndCoordinates") ) {
      for (i=1; i<=Ratomnum; i++){
        fscanf(fp,"%i %s %lf %lf %lf %lf %lf",&j,Species,
	       &Gxyz[Catomnum+Latomnum+i][1],
               &Gxyz[Catomnum+Latomnum+i][2],
               &Gxyz[Catomnum+Latomnum+i][3],
               &InitN_USpin[Catomnum+Latomnum+i],
               &InitN_DSpin[Catomnum+Latomnum+i]);
        WhatSpecies[Catomnum+Latomnum+i] = Species2int(Species);
#ifdef TRAN
        TRAN_region[Catomnum+Latomnum+i]= 3;
        TRAN_Original_Id[Catomnum+Latomnum+i]= j;
#endif

        if (i!=j){
          printf("Error of sequential number %i in <RightLeadAtoms.SpeciesAndCoordinates\n",j);
          po++;
        }

        if (2<=level_stdout){
          printf("<Input_std> R_AN=%2d T_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
		 i,Catomnum+Latomnum+i,
                  WhatSpecies[Catomnum+Latomnum+i],
                  InitN_USpin[Catomnum+Latomnum+i],
                  InitN_DSpin[Catomnum+Latomnum+i]);
        }
      }

      if (!input_last("RightLeadAtoms.SpeciesAndCoordinates>")) {
        /* format error */
        printf("Format error for RightLeadAtoms.SpeciesAndCoordinates\n");
        po++;
      }
    }

    if (coordinates_unit==0){
      for (i=1; i<=atomnum; i++){
        Gxyz[i][1] = Gxyz[i][1]/BohrR;
        Gxyz[i][2] = Gxyz[i][2]/BohrR;
        Gxyz[i][3] = Gxyz[i][3]/BohrR;
      }
    }

    /****************************************************
                          Unit cell
    ****************************************************/
    
    s_vec[0]="Ang"; s_vec[1]="AU";
    i_vec[1]=0;  i_vec[1]=1;
    input_string2int("Atoms.UnitVectors.Unit",&unitvector_unit,2,s_vec,i_vec);

    /* center */
    if (fp=input_find("<Atoms.Unitvectors")) {
      for (i=1; i<=3; i++){
        fscanf(fp,"%lf %lf %lf",&tv[i][1],&tv[i][2],&tv[i][3]);
      }
      if ( ! input_last("Atoms.Unitvectors>") ) {
        /* format error */
        printf("Format error for Atoms.Unitvectors\n");
        po++;
      }
    }

    /* left */
    if (fp=input_find("<LeftLeadAtoms.Unitvectors")) {
      for (i=1; i<=3; i++){
        fscanf(fp,"%lf %lf %lf",&Left_tv[i][1],&Left_tv[i][2],&Left_tv[i][3]);
      }
      if ( ! input_last("LeftLeadAtoms.Unitvectors>") ) {
        /* format error */
        printf("Format error for LeftLeadAtoms.Unitvectors\n");
        po++;
      }
    }

    /* right */
    if (fp=input_find("<RightLeadAtoms.Unitvectors")) {
      for (i=1; i<=3; i++){
        fscanf(fp,"%lf %lf %lf",&Right_tv[i][1],&Right_tv[i][2],&Right_tv[i][3]);
      }
      if ( ! input_last("RightLeadAtoms.Unitvectors>") ) {
        /* format error */
        printf("Format error for RightLeadAtoms.Unitvectors\n");
        po++;
      }
    }

    if (unitvector_unit==0){
      for (i=1; i<=3; i++){
        tv[i][1] = tv[i][1]/BohrR;
        tv[i][2] = tv[i][2]/BohrR;
        tv[i][3] = tv[i][3]/BohrR;
      }

      for (i=1; i<=3; i++){
        Left_tv[i][1] = Left_tv[i][1]/BohrR;
        Left_tv[i][2] = Left_tv[i][2]/BohrR;
        Left_tv[i][3] = Left_tv[i][3]/BohrR;
      }

      for (i=1; i<=3; i++){
        Right_tv[i][1] = Right_tv[i][1]/BohrR;
        Right_tv[i][2] = Right_tv[i][2]/BohrR;
        Right_tv[i][3] = Right_tv[i][3]/BohrR;
      }
    }

    /* large cell = left lead + central region + right lead */    
    { int idim=1;

    Length_C = sqrt(Dot_Product(tv[idim],tv[idim]));
    Length_L = sqrt(Dot_Product(Left_tv[idim],Left_tv[idim]));
    Length_R = sqrt(Dot_Product(Right_tv[idim],Right_tv[idim]));

    angleCL = Dot_Product(tv[idim], Left_tv[idim])/Length_C/Length_L;
    angleCR = Dot_Product(tv[idim],Right_tv[idim])/Length_C/Length_R;

    if (0.0<angleCL) Lsign =  1.0;
    else             Lsign = -1.0;
    if (0.0<angleCR) Rsign =  1.0;
    else             Rsign = -1.0;

    tv[idim][1] = tv[idim][1] + Lsign*Left_tv[idim][1] + Rsign*Right_tv[idim][1];
    tv[idim][2] = tv[idim][2] + Lsign*Left_tv[idim][2] + Rsign*Right_tv[idim][2];
    tv[idim][3] = tv[idim][3] + Lsign*Left_tv[idim][3] + Rsign*Right_tv[idim][3];
    }

}
