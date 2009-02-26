#include "f77func.h"
#ifdef blaswrap
#define zgemm_ f2c_zgemm
#define zcopy_ f2c_zcopy
#endif

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

#define z_mul_inline(x,y,z) { z.r= x.r*y.r-x.i*y.i;  z.i = x.r * y.i + x.i*y.r; }
#define z_exp_inline(x,z) {double _t; _t=exp(x.r);  z.r = _t*cos(x.i); z.i = _t*sin(x.i); }

#ifndef __sqr_definition___
#define sqr(x)   ( (x)*(x) )
#define __sqr_definition___
#endif


#ifndef Host_ID
#define Host_ID 0
#endif


#ifndef YOUSO10
#define YOUSO10 100
#endif

#ifndef Shift_K_Point
#define Shift_K_Point    1.0e-6      /* disturbance for stabilization of eigenvalue routine */
#endif


int Lapack_LU_Zinverse(int , dcomplex *);



/*** TRAN_PROTOTYPES ***/

/* TRAN_Allocate.c  */
void TRAN_Allocate_Atoms(
    int atomnum    
);

/* TRAN_Allocate.c  */
void TRAN_Deallocate_Atoms( void );

/* TRAN_Allocate.c  */
void TRAN_Allocate_Cregion( 
     MPI_Comm mpi_comm_level1,  
     int  SpinP_switch, 
     int atomnum,
     int *WhatSpecies,
     int *Spe_Total_CNO 
     ) ;

/* TRAN_Allocate.c  */
void TRAN_Deallocate_Cregion(int SpinP_switch);

/* in TRAN_Allocate.c  */
void TRAN_Allocate_Lead_Region( MPI_Comm mpi_comm_level1 ); 

/* in TRAN_Allocate.c  */
void TRAN_Deallocate_Lead_Region();


/* TRAN_Apply_Bias2e.c  */
void  TRAN_Apply_Bias2e(
       MPI_Comm comm1,
       double voltage,  
       int SpinP_switch,
       int atomnum,
       int *WhatSpecies,
       int *Spe_Total_CNO,
       int *FNAN,
       int **natn,
       int Ngrid1,
       int Ngrid2,
       int Ngrid3,
        double ****OLP,
        double *ChemP, 
        double *****H, 
        double *dVHart_Grid
);

/* TRAN_Calc_CentGreen.c  */
void tran_calc_centgreen(
   dcomplex *w, int *nc, dcomplex *sigmaL, dcomplex *sigmaR,
         dcomplex *HCC, dcomplex *SCC, dcomplex *GC);

/* TRAN_Calc_CentGreen.c  */
void TRAN_Calc_CentGreen(
                      dcomplex w,
                      int nc, 
                      dcomplex *sigmaL,
                      dcomplex *sigmaR, 
                      dcomplex *HCC,
                      dcomplex *SCC,
                      dcomplex *GC  
                      );

/* TRAN_Calc_CentGreenLesser.c  */
void TRAN_Calc_CentGreenLesser(
                      dcomplex w,
                      double ChemP_e[2],
                      int nc, 
                      dcomplex *SigmaL,
                      dcomplex *SigmaR, 
                      dcomplex *GC,  
                      dcomplex *v1,  
                      dcomplex *Gless  
                      );


/* Taisuke Ozaki Copyright (C) */
/* TRAN_Calc_GC_LorR.c  */
void TRAN_Calc_GC_LorR(
		       int iw_method,       /*  input */
                       dcomplex w,          /*  input */
                       double ChemP_e[2],   /*  input */
                       int nc,              /*  input */
                       int ne[2],           /*  input */
                       dcomplex *SigmaLorR, /*  input */
                       dcomplex *GC,        /*  input */
                       dcomplex *HCC,         /*  input */ 
                       dcomplex *SCC,         /*  input */ 
                       dcomplex *v1,        /* work, nc*nc */
                       dcomplex *GCLorR     /*  output */
		       );
 
/* TRAN_Calc_GridBound.c  */
void TRAN_Calc_GridBound(
  MPI_Comm mpi_comm_level1,
  int atomnum,
  int *WhatSpecies,
  double *Spe_Atom_Cut1,
  int Ngrid1, 
  double *Grid_Origin, 
  double **Gxyz,
  double tv[4][4],
  double gtv[4][4],
  double Right_tv[4][4]
);

/* TRAN_Calc_OneTransmission.c  */
void TRAN_Calc_OneTransmission(
   int nc, 
   dcomplex *SigmaL_R,   /* at w, changed when exit */
   dcomplex *SigmaL_A,   /* at w, changed when exit */
   dcomplex *SigmaR_R,   /* at w, changed when exit */
   dcomplex *SigmaR_A,   /* at w, changed when exit */
   dcomplex *GC_R,       /* at w, changed when exit */
   dcomplex *GC_A,       /* at w, changed when exit */
   dcomplex *v1,         /* work */
   dcomplex *v2,         /* work */
   dcomplex *value       /* output, transmission */
   );


/* TRAN_Calc_OneTransmission2.c  */
void TRAN_Calc_OneTransmission2(
				dcomplex w,
				int nc, 
				dcomplex *SigmaL,   /* at w, changed when exit */
				dcomplex *SigmaR,   /* at w, changed when exit */
				double *HCC,
				double *SCC,
				dcomplex *GR,       /* at w, changed when exit */
				dcomplex *v1,       /* work */
				dcomplex *v2,       /* work */
				dcomplex *value     /* output, transmission */
				);


/* TRAN_Calc_SelfEnergy.c  */
void tran_calc_selfenergy(
            dcomplex *w,
            int *ne,
            dcomplex *gr,  
            int *nc,
            dcomplex *hce,    
            dcomplex *sce,    
            dcomplex *sigma   
          );

/* TRAN_Calc_SelfEnergy.c  */
void TRAN_Calc_SelfEnergy( 
            dcomplex w,  
            int ne,     
            dcomplex *gr,  
            int nc,
            dcomplex *hce,    
            dcomplex *sce,    
            dcomplex *sigma   
          );

/* TRAN_Calc_SurfGreen.c  */
void tran_calc_surfgreen_direct(
                                dcomplex *w,
                                int *n,
				dcomplex *h00, 
				dcomplex *h01,
				dcomplex *s00,
				dcomplex *s01,
                                int *iteration_max,
                                double *eps,
                                dcomplex *gr  
                                );

/* TRAN_Calc_SurfGreen.c  */
void TRAN_Calc_SurfGreen_direct(
				dcomplex w,
				int n, 
				dcomplex *h00, 
				dcomplex *h01,
				dcomplex *s00,
				dcomplex *s01,
                                int iteration_max,
                                double eps,
	 			dcomplex *gr  
				);

/* TRAN_Calc_Transmission.c  */
double TRAN_Calc_Transmission(
		int iter, 
		int SpinP_switch,
		double *****nh,   
		double *****ImNL,  
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
);

/* TRAN_Calc_VHartree_G0.c  */
void TRAN_Calc_VHartree_G0(dcomplex vl, dcomplex vr,    
                      double vl_c_r,  double vl_c_i,    
                      double vr_c_r,  double vr_c_i,     
                      int ix0, int ixlp1  , double gtv1[4],
                      double *pot_r, double *pot_i   
);

/* TRAN_Calc_VHartree_Gnon0.c  */
void TRAN_Calc_VHartree_Gnon0(
			      int k2, int k3,   
			      int ix0, int ixlp1, double gtv[4][4],
			      double   v1_r, double v1_i, double v2_r, double v2_i, 
			      dcomplex v1_c, dcomplex v2_c,
			      double *ReV2, double *ImV2)  ;

/* TRAN_Connect_Read_Density.c  */
void TRAN_Connect_Read_Density(
    char *filename,
    int SpinP_switch,
    int Ngrid1, 
    int Ngrid2, 
    int Ngrid3,
    double *ChemP,
    double *minE, 
    double *dVHart_Grid,
    double **Vpot_Grid,
    double **Density_Grid
);

/* TRAN_Connect_Read_Hamiltonian.c  */
/* *static
 * void compare_and_print_error(char *buf,char *str,int val1,int val2);*/
/* TRAN_Connect_Read_Hamiltonian.c  */
void TRAN_Connect_Read_Hamiltonian(
    char *filename, 
    int SpinP_switch, 
    int *WhatSpecies,
    int *FNAN, 
    int **natn,
    int **ncn, 
    int *Spe_Total_CNO, 
    double *****H, 
    double ****OLP
);

/* TRAN_Credit.c  */
void TRAN_Credit(MPI_Comm comm1);

/* TRAN_DFT.c  */
/* *static
 * void TRAN_Set_CDM(
 *    MPI_Comm comm1, 
 *    int spin,
 *    int Matomnum,
 *    int *M2G,
 *    int *WhatSpecies,
 *    int *Spe_Total_CNO,
 *    int *MP,
 *    int *FNAN,
 *    int**natn,
 *    dcomplex *v,
 *    int NUM_c,
 *    dcomplex w_weight, 
 *    int mode,   
 *    double *****CDM,
 *    double *****EDM,
 *    double ***TRAN_DecMulP, 
 *    double Eele0[2], double Eele1[2], 
 *    double ChemP_e0[2]
 *);*/
/* TRAN_DFT.c  */
double TRAN_DFT(
                MPI_Comm comm1,
                int level_stdout,
		int iter, 
		int SpinP_switch,
		double *****nh,   
		double *****ImNL,  
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
                int *F_G2M, 
		int **atv_ijk,
		int *List_YOUSO,
		double *****CDM,   
		double *****EDM,   
                double ***TRAN_DecMulP,
		double Eele0[2], double Eele1[2], 
                double ChemP_e0[2]);

/* TRAN_DFT_Dosout.c  */
double TRAN_DFT_Dosout(
                MPI_Comm comm1,
                int level_stdout,
		int iter, 
		int SpinP_switch,
		double *****nh,   
		double *****ImNL,  
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
		double *****CDM,   
		double *****EDM,   
		double Eele0[2], double Eele1[2])  ;

/* TRAN_Distribute_Node.c  */
void TRAN_Distribute_Node(
   int Start,int End,
   int numprocs,
   int *IDStart,
   int *IDEnd
 );

/* TRAN_Distribute_Node.c  */
void TRAN_Distribute_Node_Idx(
   int Start,int End,
   int numprocs,
   int eachiwmax,
   int **Idxlist
 );

/* TRAN_FFT_interpolation3c.c  */
void zfft1(dcomplex *data,   
   int n_in, 
   dcomplex *dataout,   
   int isign0);

/* TRAN_FFT_interpolation3c.c  */
void      zfourier1a(dcomplex *in,   
                     int n_in,       
		     dcomplex *out,   
                     int n_out,      
		     int isign,
		     double alpha,double beta,double gamma) ;

/* TRAN_FFT_interpolation3c.c  */
void TRAN_FFT_Dinterpolation3D(
  MPI_Comm comm1,
  int n_in[3],
  double *din, 
  int n_out[4], 
  double *dout,  
  double alpha[3],double beta[3],double gamma[3]  );

/* TRAN_Input_std.c  */
void TRAN_Input_std(
  MPI_Comm comm1, 
  int Solver,           
  int SpinP_switch,  
  char *filepath,
  double kBvalue,
  double Electronic_Temperature
);

/* TRAN_Input_std_Atoms.c  */
void TRAN_Input_std_Atoms(
   int Solver
);

/* TRAN_Interp_ElectrodeDensity_Grid.c  */
/* *static
 * void    Convert_Cell_tv(double tv[4][4], double pos[4],
 *    double xyz[4]);*/
/* TRAN_Interp_ElectrodeDensity_Grid.c  */
void TRAN_Interp_ElectrodeDensity_Grid(
   MPI_Comm   comm1,
   double offset[4],   
   double *Density_Grid_e ,   
   int Ngrid1_e, int Ngrid2_e, int Ngrid3_e,  
   double tv_e[4][4],   
   double gtv_e[4][4],   
   double *Density_Grid,  
   int Ngrid1, int Ngrid2, int Ngrid3,
   int l1boundary[2],  
   double tv[4][4],  
   double gtv[4][4]  
);

/* TRAN_Output_HKS.c  */
int TRAN_Output_HKS(char *fileHKS);

/* TRAN_Output_HKS_Write_Grid.c  */
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
				);

/* TRAN_Output_Trans_HS.c  */
void  TRAN_Output_Trans_HS(
        MPI_Comm comm1,
        int SpinP_switch, 
        double ChemP ,
        double *****H,
        double *****OLP,
        int atomnum,
        int SpeciesNum,
        int *WhatSpecies,
        int *Spe_Total_CNO,
        int *FNAN,
        int **natn,
        int **ncn,
        int *G2ID,
        int **atv_ijk,
        int Max_FSNAN,
        double ScaleSize,
        int *F_G2M,      
        int TCpyCell,
        int *List_YOUSO,
        char *filepath,
        char *filename,
        char *fname  );

/* TRAN_Output_Transmission.c  */ 
void TRAN_Output_Transmission(int SpinP_switch);

/* TRAN_Overwrite_Densitygrid.c  */
void TRAN_Overwrite_Densitygrid(
            MPI_Comm comm1,
            int SpinP_switch,
            int Ngrid1,
            int Ngrid2,
            int Ngrid3,
            int Num_Cells0, 
            int *My_Cell0, 
            int *My_Cell1,
            double **Density_Grid 
);

/* TRAN_Overwrite_V2.c  */
void  TRAN_Overwrite_V2(
   int Ngrid1,
   int Ngrid2,
   int Ngrid3,
   int l1boundary[2],
   int My_NGrid2_Poisson, 
   int Start2, 
   dcomplex *Grid[2],
   double ***ReV2,   
   double ***ImV2 )  ;

/* TRAN_Poisson.c  */
double TRAN_Poisson(
                    int fft_charge_flag,
		    double ***ReV1, double ***ImV1,
		    double ***ReV2, double ***ImV2);

/* TRAN_Print.c  */
void TRAN_Print2_set_eps(double e1);

/* TRAN_Print.c  */
void TRAN_Print2_set_max(int m1);

/* TRAN_Print.c  */
void TRAN_Print2_dcomplex(char *name, int n1,int n2,dcomplex *gc);

/* TRAN_Print.c  */
void TRAN_Print2_double(char *name, int n1,int n2,double *gc);

/* TRAN_Print.c  */
void TRAN_Print2_dx_dcomplex(char *name, int n1,int dx1, int n2,int dx2, dcomplex *gc);

/* TRAN_Print.c  */
void TRAN_Print2_dx_double(char *name, int n1,int dx1,int n2,int dx2,double *gc);

/* TRAN_Print.c  */
void TRAN_FPrint2_double(char *name, int n1,int n2,double *gc);

/* TRAN_Print.c  */
void TRAN_FPrint2_dcomplex(char *name, int n1,int n2,dcomplex *gc);

/* TRAN_Print.c  */
void TRAN_FPrint2_binary_double(FILE *fp, int n1,int n2,double *gc);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid_Cell1(
  char *filename,
  int n1, int n2, int n3, 
   int *My_Cell1,
   double *realgrid
);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid_Cell0(
   char *filename,
   double origin[4],  
   double gtv[4][4],  
   int Ngrid1, int Ngrid2, int Ngrid3s, int Ngrid3e,   
   double  R[4],   
   int *Cell0,
   double *grid_value
);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid(
   char *filename,
   double origin[4],  
   double gtv[4][4],  
   int Ngrid1, int Ngrid2, int Ngrid3s, int Ngrid3e,   
   double  R[4],   
   double *grid_value 
);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid_z(
   char *filename,
   double origin[4],  
   double gtv[4][4],  
   int Ngrid1, int Ngrid2, int Ngrid3s, int Ngrid3e,   
   double  R[4],   
   double *grid_value
);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid_c(
   char *filenamer,
   char *filenamei,
   double origin[4],  
   double gtv[4][4],  
   int Ngrid1, int Ngrid2, int Ngrid3s, int Ngrid3e,   
   double  R[4],   
   dcomplex *grid_value
);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid_v(
   char *filename,
   double origin[4],  
   double gtv[4][4],  
   int Ngrid1, int Ngrid2, int Ngrid3s, int Ngrid3e,   
   double  R[4],   
   double ***grid_value
);

/* TRAN_Print_Grid.c  */
void TRAN_Print_Grid_Startv(
   char *filename,
   int Ngrid1, int Ngrid2,  int Ngrid3,   
   int Start,
   double ***grid_value
);

/* TRAN_Read.c  */
void TRAN_Read_double(char *str, int n1,int n2, double *a);

/* TRAN_Read.c  */
void TRAN_Read_dcomplex(char *str, int n1,int n2, dcomplex *a);

/* TRAN_Read.c  */
void TRAN_FRead2_binary_double(FILE *fp, int n1, int n2, double *gc);

/* TRAN_RestartFile.c  */
int TRAN_Input_HKS( MPI_Comm comm1, char *fileHKS);

/* TRAN_RestartFile.c  */
int TRAN_RestartFile(MPI_Comm comm1, char *mode, char *position,char *filepath, char *filename);

/* TRAN_Set_CentOverlap.c  */
void TRAN_Set_CentOverlap( 
			   MPI_Comm comm1,
			   int job, 
			   int SpinP_switch, 
                           double k2,
                           double k3,
                           int *order_GA,
			   double **H1,
			   double *S1,
			   double *****H,  
			   double ****OLP,  
			   int atomnum,
			   int Matomnum,
 			   int *M2G, 
			   int *G2ID, 
			   int *WhatSpecies,
			   int *Spe_Total_CNO,
			   int *FNAN,
			   int **natn,
			   int **ncn, 
			   int **atv_ijk
			   );

/* TRAN_Set_Electrode_Grid.c  */
void TRAN_Set_Electrode_Grid(
                MPI_Comm comm1,
                double *Grid_Origin,  
                double tv[4][4],       
                double Left_tv[4][4],  
                double Right_tv[4][4],  
                double gtv[4][4],       
                int Ngrid1,
                int Ngrid2,
                int Ngrid3,      
                double *Density_Grid   
);

/* TRAN_Set_Electrode_Grid.c  */
/* *static
 * void  TRAN_FFT_Electrode_Grid(MPI_Comm comm1, int isign);*/
/* TRAN_Set_IntegPath.c  */
void TRAN_Set_IntegPath_Square(void);

/* TRAN_Set_IntegPath.c  */
void      TRAN_Set_IntegPath_ThermalArc(void);

/* TRAN_Set_IntegPath.c  */
void TRAN_Set_IntegPath( double kBvalue, double Electronic_Temperature );

/* TRAN_Set_MP.c  */
void TRAN_Set_MP(
        int job, 
        int atomnum, int *WhatSpecies, int *Spe_Total_CNO, 
        int *NUM,   
        int *MP     
);

/* TRAN_Set_PathEnergyStr.c  */
/* *static
 * void TRAN_error_and_exit(
 *    char *buf
 *);*/
/* TRAN_Set_PathEnergyStr.c  */
void TRAN_Set_PathEnergyStr_Square(
   int m,
   char **str,
   double default_relative_ene[4],  
   double tran_square_path_ene[4],  
   int    tran_square_path_ene_fix[4]
);

/* TRAN_Set_PathEnergyStr.c  */
void TRAN_Set_PathEnergyStr(
   int tran_integ_pathtype, 
   int m,
   char **str,
   double default_relative_ene[4],  
   double tran_square_path_ene[4],  
   int  tran_square_path_ene_fix[4]  
);

/* TRAN_Set_SurfOverlap.c  */
void TRAN_Set_SurfOverlap( MPI_Comm comm1, char *position, double k2, double k3 );

/* TRAN_Set_Value.c  */
void TRAN_Set_Value_double(dcomplex *A, int n, double a, double b);

/* TRAN_adjust_Grid_Origin.c  */
void TRAN_adjust_Grid_Origin( MPI_Comm comm1, double Grid_Origin[4]);

/* TRAN_adjust_Ngrid.c  */
void  TRAN_adjust_Ngrid( MPI_Comm comm1, int *Ngrid1,int *Ngrid2, int *Ngrid3);


/* TRAN_Check_Region_Lead.c */
int TRAN_Check_Region_Lead(
		  int atomnum,
		    int *WhatSpecies, 
		      double *Spe_Atom_Cut1,
		        double **Gxyz,
			  double tv[4][4]
		);

/* TRAN_Check_Region.c */
int TRAN_Check_Region(
                      int atomnum,
                      int *WhatSpecies, 
                      double *Spe_Atom_Cut1,
                      double **Gxyz
                      );
