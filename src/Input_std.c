#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "openmx_common.h"
#include "Inputtools.h"

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
#endif

#define MAXBUF 1024


void SpeciesString2int(int p);
void kpath_changeunit(double tv[4][4],double tv0[4][4],int Band_Nkpath,
                      double ***Band_kpath);
void kpoint_changeunit(double tv[4][4],double tv0[4][4],int MO_Nkpoint,
                       double **MO_kpoint);
void Set_Cluster_UnitCell(double tv[4][4],int unitflag);

int OrbPol2int(char OrbPol[YOUSO10]);
char *ToCapital(char *s);
int divisible_cheker(int N);
void Setup_Mixed_Basis(char *file, int myid);
static void Set_In_First_Cell();

void Input_std(char *file)
{
  FILE *fp,*fp2;
  int i,j,k,k1,k2,k3,itmp;
  int l,mul; /* added by MJ */
  int po=0;  /* error count */
  double r_vec[20],r_vec2[20];
  int i_vec[20],i_vec2[20];
  char *s_vec[20],Species[YOUSO10];
  char OrbPol[YOUSO10];
  char c; 
  double ecutoff1dfft;
  double sn1,sn2,sn3;
  double mx,my,mz,tmp;
  double S_coordinate[3];
  double Length_C,Length_L,Length_R;
  double angleCL,angleCR,Lsign,Rsign;
  double length,x,y,z;
  double xc,yc,zc,xm,ym,zm;
  int orbitalopt;
  char buf[MAXBUF];
  int numprocs,myid;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID){  
    printf("*******************************************************\n"); 
    printf("       read the input file and initializing            \n");
    printf("*******************************************************\n\n"); 
  }

  /****************************************************
                       open a file
  ****************************************************/

  input_open(file);

  input_string("System.CurrrentDirectory",filepath,"./");
  input_string("System.Name",filename,"default");
  input_string("DATA.PATH",DFT_DATA_PATH,"../DFT_DATA");
  input_int("level.of.stdout", &level_stdout,1);
  input_int("level.of.fileout",&level_fileout,1);

  if (level_stdout<1 || 3<level_stdout){
    printf("Invalid value of level.of.stdout\n");
    po++;
  }

  if (level_fileout<0 || 3<level_fileout){
    printf("Invalid value of level.of.fileout\n");
    po++;
  }

  /****************************************************
   set the number of threads in the OpenMP environment
  ****************************************************/

  /*
  input_logical("OpenMP.Threads.eq.Procs",&openmp_threads_eq_procs,1);   
  input_int("OpenMP.Threads.Number",&openmp_threads_num,1);              

  if (openmp_threads_eq_procs){
    i = omp_get_num_procs();
    omp_set_num_threads(i);
  }
  else {
    omp_set_num_threads(openmp_threads_num);  
  }
  */

  /****************************************************
               projector expansion of VNA
  ****************************************************/

  input_logical("scf.ProExpn.VNA",&ProExpn_VNA,1); /* default=on */
  input_int("scf.BufferL.VNA", &BufferL_ProVNA,5);
  input_int("scf.RadialF.VNA", &List_YOUSO[34],7);
  
  /****************************************************
                 the use of mixed basis
  ****************************************************/

  input_logical("scf.Mixed.Basis",&Mixed_Basis_flag,0); /* default=off */

  /****************************************************
                     cutoff energy 
  ****************************************************/

  /* for cutoff energy */
  
  input_double("scf.energycutoff",&Grid_Ecut,(double)150.0);
  input_logical("scf.MPI.tuned.grids",&MPI_tunedgrid_flag,0);

  /* for fixed Ngrids  */

  i_vec2[0]=0;
  i_vec2[1]=0;
  i_vec2[2]=0;
  input_intv("scf.Ngrid",3,i_vec,i_vec2);
  Ngrid1 = i_vec[0];
  Ngrid2 = i_vec[1];
  Ngrid3 = i_vec[2];

  if (Ngrid1==0 && Ngrid2==0 && Ngrid3==0)
    Ngrid_fixed_flag = 0;
  else 
    Ngrid_fixed_flag = 1;

  if (Ngrid_fixed_flag==1){
    i = divisible_cheker(Ngrid1);
    j = divisible_cheker(Ngrid2);
    k = divisible_cheker(Ngrid3);
   
    if ( (i*j*k)==0 ) {

      printf("scf.Ngrid must be divisible by \n");

      printf("    ");
      for (i=0; i<NfundamentalNum; i++){
        printf("%3d ",fundamentalNum[i]);
      } 
      printf("\n");

      MPI_Finalize(); 
      exit(0);
    }

    /* for mixed basis set */

    if (Mixed_Basis_flag==1){
     
      if ( ((Ngrid1%NG_Mixed_Basis)!=0)
        || ((Ngrid2%NG_Mixed_Basis)!=0)
	|| ((Ngrid3%NG_Mixed_Basis)!=0)){

        printf("scf.Ngrid must be divisible by %2d for the mixed basis\n",NG_Mixed_Basis);

        MPI_Finalize(); 
        exit(0);
      }
    }
  }

  /****************************************************
   adjustment of cutoff energy for mixed basis scheme
  ****************************************************/

  if (Mixed_Basis_flag==1){
    Setup_Mixed_Basis(file,myid);
  }

  /****************************************************
              definition of atomic species
  ****************************************************/

  input_int("Species.Number",&SpeciesNum,0);
  real_SpeciesNum = SpeciesNum;

  if (Mixed_Basis_flag==1) SpeciesNum++;

  if (SpeciesNum<=0){
    printf("Species.Number may be wrong.\n");
    po++;
  }
  List_YOUSO[18] = SpeciesNum;

  /* memory allocation */
  Allocate_Arrays(0);

  /*************************************************************
     for LDA+U
     Hub_U_switch should be called before Allocate_Arrays(1);
  *************************************************************/ 

  input_logical("scf.Hubbard.U",&Hub_U_switch, 0);     /* --- MJ */

  /* default Hub_U_occupation = 2; */

  s_vec[0]="DUAL";            i_vec[0]=2;
  s_vec[1]="ONSITE";          i_vec[1]=0;
  s_vec[2]="FULL" ;           i_vec[2]=1;

  input_string2int("scf.Hubbard.Occupation",&Hub_U_occupation, 3, s_vec,i_vec);

  /*************************************************************
                           read species
  *************************************************************/ 

  if (fp=input_find("<Definition.of.Atomic.Species")) {

    for (i=0; i<SpeciesNum; i++){

      if (Mixed_Basis_flag==1 && i==(SpeciesNum-1)){

        sprintf(SpeName[i], "FEB");
        sprintf(SpeBasis[i],"FEB-s1");
        sprintf(SpeVPS[i],  "FEB");
        SpeciesString2int(i);
      }
      else{
        fscanf(fp,"%s %s %s",SpeName[i],SpeBasis[i],SpeVPS[i]);
        SpeciesString2int(i);
      }
    }

    if (! input_last("Definition.of.Atomic.Species>")) {
      /* format error */

      po++;

      if (myid==Host_ID){
        printf("Format error for Definition.of.Atomic.Species\n");
      }
      MPI_Finalize();
      exit(0);
    }
  }
  if (2<=level_stdout){
    for (i=0; i<SpeciesNum; i++){
      printf("<Input_std>  %i Name  %s\n",i,SpeName[i]);
      printf("<Input_std>  %i Basis %s\n",i,SpeBasis[i]);
      printf("<Input_std>  %i VPS   %s\n",i,SpeVPS[i]);
    }
  }

  List_YOUSO[35] = 0;
  for (i=0; i<SpeciesNum; i++){
    if (List_YOUSO[35]<Spe_MaxL_Basis[i]) List_YOUSO[35] = Spe_MaxL_Basis[i];
  }
  List_YOUSO[35] = List_YOUSO[35] + BufferL_ProVNA;

  /****************************************************
       Molecular dynamics or geometry optimization
  ****************************************************/

  i=0;
  s_vec[i]="NOMD";                    i_vec[i]=0;  i++;
  s_vec[i]="NVE" ;                    i_vec[i]=1;  i++;
  s_vec[i]="NVT_VS";                  i_vec[i]=2;  i++; /* modified by mari */
  s_vec[i]="OPT";                     i_vec[i]=3;  i++;
  s_vec[i]="EF";                      i_vec[i]=4;  i++; 
  s_vec[i]="BFGS";                    i_vec[i]=5;  i++; 
  s_vec[i]="RF";                      i_vec[i]=6;  i++; /* RF method by hmweng */
  s_vec[i]="DIIS";                    i_vec[i]=7;  i++;
  s_vec[i]="Constraint_DIIS";         i_vec[i]=8;  i++; /* not used */
  s_vec[i]="NVT_NH";                  i_vec[i]=9;  i++; 
  s_vec[i]="Opt_LBFGS";               i_vec[i]=10; i++; 

  j = input_string2int("MD.Type",&MD_switch, i, s_vec,i_vec);
  if (j==-1){
    MPI_Finalize();
    exit(0);
  }

  input_int("MD.maxIter",&MD_IterNumber,1);
  if (MD_IterNumber<1){
    printf("MD_IterNumber=%i should be over 0.\n",MD_IterNumber);
    po++;
  }

  input_double("MD.TimeStep",&MD_TimeStep,(double)0.5);
  if (MD_TimeStep<0.0){
    printf("MD.TimeStep=%lf should be over 0.\n",MD_TimeStep);
    po++;
  }

  input_double("MD.Opt.criterion",&MD_Opt_criterion,(double)0.0003);
  input_int("MD.Opt.DIIS.History",&M_GDIIS_HISTORY,3);
  input_int("MD.Opt.StartDIIS",&OptStartDIIS,5);
  input_int("MD.Opt.EveryDIIS",&OptEveryDIIS,200);

  /*
  input_double("MD.Initial.MaxStep",&SD_scaling_user,(double)0.02); 
  */

  /* Ang -> a.u. */
  /*
  SD_scaling_user /= BohrR;
  */

  /*
  input_double("MD.Opt.DIIS.Mixing",&Gdiis_Mixing,(double)0.1);
  */

  if (19<M_GDIIS_HISTORY){
    printf("MD.Opt.DIIS.History should be lower than 19.");
    MPI_Finalize();
    exit(0);
  }

  if (MD_switch==2 || MD_switch==9){ 

    if (fp=input_find("<MD.TempControl")) {

      fscanf(fp,"%i",&TempNum);         

      /* added by mari */
      if (MD_switch==2) { /* NVT_VS */

	NumScale[0] = 0;

	for (i=1; i<=TempNum; i++){  
	  fscanf(fp,"%d %d %lf %lf",&NumScale[i],&IntScale[i],
                                    &TempScale[i],&RatScale[i]);
	  TempPara[i][1] = NumScale[i];
	  TempPara[i][2] = TempScale[i];
	}

        TempPara[0][1] = 0;
	TempPara[0][2] = TempPara[1][2];
      }

      /* added by mari */
      else if (MD_switch==9) { /* NVT_NH */
	for (i=1; i<=TempNum; i++){  
	  fscanf(fp,"%lf %lf",&TempPara[i][1],&TempPara[i][2]);
	}  
      }

      if ( ! input_last("MD.TempControl>") ) {
	/* format error */
	printf("Format error for MD.TempControl\n");
	po++;
      }

    }
  }

  if (fp=input_find("<MD.CellPressureControl")) {
    fscanf(fp,"%i",&PreNum);  
    for (i=1; i<=TempNum; i++){  
      fscanf(fp,"%lf %lf",&PrePara[i][1],&PrePara[i][2]);
    }  
    if ( ! input_last("MD.CellPressureControl>") ) {
      /* format error */
      printf("Format error for MD.CellPressureControl\n");
      po++;
    }
  }

  input_double("NH.Mass.HeatBath",&TempQ,(double)20.0);

  /****************************************************
             solver of the eigenvalue problem
  ****************************************************/

  Solver_DIIS_flag = 0;

  s_vec[0]="Recursion";     s_vec[1]="Cluster"; s_vec[2]="Band";
  s_vec[3]="NEGF";          s_vec[4]="DC";      s_vec[5]="GDC";
  s_vec[6]="Cluster-DIIS";  s_vec[7]="Krylov";
  
  i_vec[0]=1;  i_vec[1]=2;  i_vec[2]=3;
  i_vec[3]=4;  i_vec[4]=5;  i_vec[5]=6;
  i_vec[6]=7;  i_vec[7]=8;

  input_string2int("scf.EigenvalueSolver", &Solver, 8, s_vec,i_vec);

  if (Solver==7){
    Solver = 2;
    Solver_DIIS_flag = 1;
  }

  if (Mixed_Basis_flag==1 && Solver==4){

    if (myid==Host_ID){
      printf("The mixed basis is not supported for NEGF\n");
    }

    MPI_Finalize();
    exit(0);
  }

  /* default=dstevx */
  s_vec[0]="dstevx"; s_vec[1]="dstegr"; s_vec[2]="dstedc"; s_vec[3]="dsteqr"; 
  i_vec[0]=2;        i_vec[1]=0;        i_vec[2]=1;        i_vec[3]=3;        
  input_string2int("scf.lapack.dste", &dste_flag, 4, s_vec,i_vec);

  if (Solver==1){
    if (myid==Host_ID){
      printf("Recursion method is not supported in this version.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /****************************************************
      for generation of Monkhorst-Pack k-points
  ****************************************************/

  input_double("scf.MP.criterion",&Criterion_MP_Special_Kpt,(double)1.0e-5);

  s_vec[0]="REGULAR"; s_vec[1]="MP";
  i_vec[0]=1        ; i_vec[1]=2   ; 
  input_string2int("scf.Generation.Kpoint", &way_of_kpoint, 2, s_vec,i_vec);

  if (Solver==4 && way_of_kpoint==2){
    if (myid==Host_ID){
      printf("The Monkhorst-Pack is not supported for NEGF.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /****************************************************
   flags for saving and reading Fourier transformed
   quantities generated in FT_*.c
  ****************************************************/

  input_logical("FT.files.save",&FT_files_save,0);
  input_logical("FT.files.read",&FT_files_read,0);

  /****************************************************
                SCF or electronic system
  ****************************************************/

  s_vec[0]="LDA"; s_vec[1]="LSDA-CA"; s_vec[2]="LSDA-PW"; s_vec[3]="GGA-PBE";
  i_vec[0]=1    ; i_vec[1]=2     ; i_vec[2]= 3   ; i_vec[3]= 4   ;
  input_string2int("scf.XcType", &XC_switch, 4, s_vec,i_vec);

  s_vec[0]="Off"; s_vec[1]="On"; s_vec[2]="NC";
  i_vec[0]=0    ; i_vec[1]=1   ; i_vec[2]=3;
  input_string2int("scf.SpinPolarization", &SpinP_switch, 3, s_vec,i_vec);
  if      (SpinP_switch==0) List_YOUSO[23] = 1;  
  else if (SpinP_switch==1) List_YOUSO[23] = 2;
  else if (SpinP_switch==3) List_YOUSO[23] = 4;

  if (XC_switch==3 && Solver==8){
    if (myid==Host_ID){
      printf("Krylov subspace method is not supported for non-collinear calculations.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if (XC_switch==1 && 1<=SpinP_switch){
    if (myid==Host_ID){
      printf("SpinP_switch should be OFF for this exchange functional.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /* scf.Constraint.NC.Spin */

  input_logical("scf.Constraint.NC.Spin",&Constraint_NCS_switch,0);

  if (SpinP_switch!=3 && Constraint_NCS_switch==1){
    if (myid==Host_ID){
      printf("The constraint scheme is not supported for a collinear DFT calculation.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if (Constraint_NCS_switch==1 && Hub_U_occupation!=2){
    if (myid==Host_ID){
      printf("The constraint scheme is supported in case of scf.Hubbard.Occupation=dual.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  input_double("scf.Constraint.NC.Spin.V",&Constraint_NCS_V,(double)0.0);
  /* eV to Hartree */
  Constraint_NCS_V = Constraint_NCS_V/eV2Hartree;

  /* scf.SpinOrbit.Coupling */

  input_logical("scf.SpinOrbit.Coupling",&SO_switch,0);

  if (SpinP_switch!=3 && SO_switch==1){
    if (myid==Host_ID){
      printf("Spin-orbit coupling is not supported for collinear DFT calculations.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  input_logical("scf.partialCoreCorrection",&PCC_switch,1);

  if (SpinP_switch==0 && SO_switch==1){
    if (myid==Host_ID){
      printf("scf.SpinOrbit.Coupling must be OFF when scf.SpinPolarization=OFF\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /* scf.NC.Zeeman.Spin */

  input_logical("scf.NC.Zeeman.Spin",&Zeeman_NCS_switch,0);

  if (SpinP_switch!=3 && Zeeman_NCS_switch==1){
    if (myid==Host_ID){
      printf("The Zeeman term is not supported for a collinear DFT calculation.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if (Zeeman_NCS_switch==1 && Hub_U_occupation!=2){
    if (myid==Host_ID){
      printf("The Zeeman term for spin is supported in case of scf.Hubbard.Occupation=dual.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if (Constraint_NCS_switch==1 && Zeeman_NCS_switch==1){
    if (myid==Host_ID){
      printf("For spin magnetic moment, the constraint scheme and the Zeeman term\n");
      printf("are mutually exclusive.  Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  input_double("scf.NC.Mag.Field.Spin",&Mag_Field_Spin,(double)0.0);

  /**************************************
     Tesla to a.u.
     1 Tesla = 1/(2.35051742*10^5) a.u. 
  ***************************************/

  Mag_Field_Spin = Mag_Field_Spin/(2.35051742*100000.0);

  /* scf.NC.Zeeman.Orbital */

  input_logical("scf.NC.Zeeman.Orbital",&Zeeman_NCO_switch,0);

  if (SpinP_switch!=3 && Zeeman_NCO_switch==1){
    if (myid==Host_ID){
      printf("The Zeeman term is not supported for a collinear DFT calculation.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if (Zeeman_NCO_switch==1 && Hub_U_occupation!=2){
    if (myid==Host_ID){
      printf("The Zeeman term for orbital is supported in case of scf.Hubbard.Occupation=dual.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if (Zeeman_NCO_switch==1 && SO_switch==0){
    if (myid==Host_ID){
      printf("The Zeeman term for orbital moment is not supported without the SO term.\n");
      printf("Check your input file.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  input_double("scf.NC.Mag.Field.Orbital",&Mag_Field_Orbital,(double)0.0);

  /**************************************
     Tesla to a.u.
     1 Tesla = 1/(2.35051742*10^5) a.u. 
  ***************************************/

  Mag_Field_Orbital = Mag_Field_Orbital/(2.35051742*100000.0);

  if      (SpinP_switch==0)                 List_YOUSO[5] = 1;
  else if (SpinP_switch==1)                 List_YOUSO[5] = 2;
  else if (SpinP_switch==3 && SO_switch==0) List_YOUSO[5] = 3;
  else if (SpinP_switch==3 && SO_switch==1) List_YOUSO[5] = 3;

  i_vec2[0]=4;
  i_vec2[1]=4;
  i_vec2[2]=4;
  input_intv("scf.Kgrid",3,i_vec,i_vec2);
  Kspace_grid1 = i_vec[0];
  Kspace_grid2 = i_vec[1];
  Kspace_grid3 = i_vec[2];

  if (Kspace_grid1<=0){
    printf("Kspace_grid1 should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 
  if (Kspace_grid2<=0){
    printf("Kspace_grid2 should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 
  if (Kspace_grid3<=0){
    printf("Kspace_grid3 should be over 1\n");
    MPI_Finalize();
    exit(1);
  } 

  if (Solver!=3 && Solver!=4){
    List_YOUSO[27] = 1;
    List_YOUSO[28] = 1;
    List_YOUSO[29] = 1;
  }
  else{
    List_YOUSO[27] = Kspace_grid1;
    List_YOUSO[28] = Kspace_grid2;
    List_YOUSO[29] = Kspace_grid3;
  }

  /* set PeriodicGamma_flag in 1 in the band calc. with only the gamma point */
  PeriodicGamma_flag = 0;
  if (Solver==3 && Kspace_grid1==1 && Kspace_grid2==1 && Kspace_grid3==1){
    printf("When only the gamma point is considered, the eigenvalue solver is changed to 'Cluster' with the periodic boundary condition.\n");
    PeriodicGamma_flag = 1;
    Solver = 2;
  }

  input_double("scf.ElectronicTemperature",&E_Temp,(double)300.0);
  E_Temp = E_Temp/eV2Hartree;
  Original_E_Temp = E_Temp;

  s_vec[0]="Simple"; s_vec[1]="RMM-DIIS"; s_vec[2]="GR-Pulay";
  s_vec[3]="Kerker"; s_vec[4]="RMM-DIISK";
  i_vec[0]=0; i_vec[1]=1; i_vec[2]=2; i_vec[3]=3; i_vec[4]=4;

  input_string2int("scf.Mixing.Type",&Mixing_switch,5,s_vec,i_vec);
  if (Mixing_switch==3 || Mixing_switch==4) Kmixing_flag = 1;
  else Kmixing_flag = 0;

  input_double("scf.Init.Mixing.Weight",&Mixing_weight,(double)0.3);
  input_double("scf.Min.Mixing.Weight",&Min_Mixing_weight,(double)0.001);
  input_double("scf.Max.Mixing.Weight",&Max_Mixing_weight,(double)0.4);
  input_double("scf.Kerker.factor",&Kerker_factor,(double)1.0);
  input_int("scf.Mixing.History",&Num_Mixing_pDM,5);
  input_int("scf.Mixing.StartPulay",&Pulay_SCF,6);
  input_int("scf.Mixing.EveryPulay",&EveryPulay_SCF,5);
  input_int("scf.ExtCharge.History",&Extrapolated_Charge_History,3);

  /* increase electric temperature in case of SCF oscillation, default=off */
  input_logical("scf.Mixing.Control.Temp", &SCF_Control_Temp, 0); 

  if (Mixing_switch==0){
    List_YOUSO[16] = 3;
    List_YOUSO[38] = 1;
  }
  else if (Mixing_switch==1 || Mixing_switch==2) {
    List_YOUSO[16] = Num_Mixing_pDM + 2;
    List_YOUSO[38] = 1;
  }
  else if (Mixing_switch==3) {
    List_YOUSO[16] = 3;
    List_YOUSO[38] = 3;
  }
  else if (Mixing_switch==4) {
    List_YOUSO[16] = 3;
    List_YOUSO[38] = Num_Mixing_pDM + 2;
  }

  input_double("scf.criterion",&SCF_Criterion,(double)1.0e-6);
  if (SCF_Criterion<0.0){
    printf("SCF_Criterion=%10.9f should be larger than 0.\n",SCF_Criterion);
    po++;
  }

  input_double("scf.system.charge",&system_charge,(double)0.0);

  /* scf.fixed.grid */

  r_vec[0]=1.0e+9; r_vec[1]=1.0e+9; r_vec[2]=1.0e+9;
  input_doublev("scf.fixed.grid",3,scf_fixed_origin,r_vec);

#ifdef TRAN
  TRAN_Input_std(mpi_comm_level1, Solver,SpinP_switch,filepath,kB,E_Temp); 
#endif

  /****************************************************
                         atoms
  ****************************************************/

  /* except for NEGF */

  if (Solver!=4){  

    /* atom */

    input_int("Atoms.Number",&atomnum,0);

    if (Mixed_Basis_flag==1) atomnum += Ngrid1_FE*Ngrid2_FE*Ngrid3_FE;

    if (atomnum<=0){
      printf("111 Atoms.Number may be wrong.\n");

      printf("Atoms.Number may be wrong.\n");
      po++;
    }
    List_YOUSO[1] = atomnum + 1;

    /* memory allocation */
    Allocate_Arrays(1);

    /* initialize */

    s_vec[0]="Ang";  s_vec[1]="AU";   s_vec[2]="FRAC";
    i_vec[0]= 0;     i_vec[1]= 1;     i_vec[2]= 2;
    input_string2int("Atoms.SpeciesAndCoordinates.Unit",
                     &coordinates_unit,3,s_vec,i_vec);

    if (fp=input_find("<Atoms.SpeciesAndCoordinates") ) {

      for (i=1; i<=atomnum; i++){
        fgets(buf,MAXBUF,fp);

	/* spin non-collinear */ 
	if (SpinP_switch==3){

	  /*******************************************************
               (1) spin non-collinear
	  *******************************************************/

	  sscanf(buf,"%i %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %s",
		 &j, Species,
		 &Gxyz[i][1],&Gxyz[i][2],&Gxyz[i][3],
		 &InitN_USpin[i],&InitN_DSpin[i],
		 &Angle0_Spin[i], &Angle1_Spin[i], 
		 &Angle0_Orbital[i], &Angle1_Orbital[i],
		 &Constraint_SpinAngle[i],
		 OrbPol );

	  /* atomMove is initialzed as 1 */

	  if (fabs(Angle0_Spin[i])<1.0e-10){
	    Angle0_Spin[i] = Angle0_Spin[i] + rnd(1.0e-5);
	  }

	  Angle0_Spin[i] = PI*Angle0_Spin[i]/180.0;
	  Angle1_Spin[i] = PI*Angle1_Spin[i]/180.0;
	  InitAngle0_Spin[i] = Angle0_Spin[i];
	  InitAngle1_Spin[i] = Angle1_Spin[i];

	  if (fabs(Angle0_Orbital[i])<1.0e-10){
	    Angle0_Orbital[i] = Angle0_Orbital[i] + rnd(1.0e-5);
	  }

	  Angle0_Orbital[i] = PI*Angle0_Orbital[i]/180.0;
	  Angle1_Orbital[i] = PI*Angle1_Orbital[i]/180.0;
	  InitAngle0_Orbital[i] = Angle0_Orbital[i];
	  InitAngle1_Orbital[i] = Angle1_Orbital[i];


          /*************************************************************************
           check whether the Euler angle measured from the direction (1,0) is used,
           if not, change the Euler angle 
	  *************************************************************************/

          if ( (InitN_USpin[i]-InitN_DSpin[i])<0.0 ){

	    mx = -sin(InitAngle0_Spin[i])*cos(InitAngle1_Spin[i]);
	    my = -sin(InitAngle0_Spin[i])*sin(InitAngle1_Spin[i]);
            mz = -cos(InitAngle0_Spin[i]);

            xyz2spherical(mx,my,mz,0.0,0.0,0.0,S_coordinate);

	    Angle0_Spin[i] = S_coordinate[1];
  	    Angle1_Spin[i] = S_coordinate[2];
	    InitAngle0_Spin[i] = Angle0_Spin[i];
  	    InitAngle1_Spin[i] = Angle1_Spin[i];

            tmp = InitN_USpin[i];
            InitN_USpin[i] = InitN_DSpin[i];
            InitN_DSpin[i] = tmp; 

	    mx = -sin(InitAngle0_Orbital[i])*cos(InitAngle1_Orbital[i]);
	    my = -sin(InitAngle0_Orbital[i])*sin(InitAngle1_Orbital[i]);
            mz = -cos(InitAngle0_Orbital[i]);

            xyz2spherical(mx,my,mz,0.0,0.0,0.0,S_coordinate);
             
            Angle0_Orbital[i] = S_coordinate[1];
            Angle1_Orbital[i] = S_coordinate[2]; 
  	    InitAngle0_Orbital[i] = Angle0_Orbital[i];
	    InitAngle1_Orbital[i] = Angle1_Orbital[i];
          }
	}

	/**************************************************
                  (2) spin collinear
	**************************************************/

	else{ 

	  sscanf(buf,"%i %s %lf %lf %lf %lf %lf %s",
		 &j, Species,
		 &Gxyz[i][1],&Gxyz[i][2],&Gxyz[i][3],
		 &InitN_USpin[i],&InitN_DSpin[i], OrbPol );
	}

        WhatSpecies[i] = Species2int(Species);

        if (Hub_U_switch==1) OrbPol_flag[i] = OrbPol2int(OrbPol);
 
        if (i!=j){
          printf("Format error of the sequential number %i in <Atoms.SpeciesAndCoordinates\n",j);
          po++;
        }

        if (2<=level_stdout){
          printf("<Input_std>  ct_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
                  i,WhatSpecies[i],InitN_USpin[i],InitN_DSpin[i]);
        }
      }

      ungetc('\n',fp);
      if (!input_last("Atoms.SpeciesAndCoordinates>")) {
        /* format error */
        printf("Format error for Atoms.SpeciesAndCoordinates\n");
        po++;
      }

    }

    /****************************************************
       the unit of atomic coordinates is transformed 
    ****************************************************/

    /*  Ang to AU */ 
    if (coordinates_unit==0){
      for (i=1; i<=atomnum; i++){
	Gxyz[i][1] = Gxyz[i][1]/BohrR;
	Gxyz[i][2] = Gxyz[i][2]/BohrR;
	Gxyz[i][3] = Gxyz[i][3]/BohrR;
      }
    }

    /****************************************************
                          unit cell
    ****************************************************/

    s_vec[0]="Ang"; s_vec[1]="AU";
    i_vec[0]=0;  i_vec[1]=1;
    input_string2int("Atoms.UnitVectors.Unit",&unitvector_unit,2,s_vec,i_vec);

    if (fp=input_find("<Atoms.Unitvectors")) {
      for (i=1; i<=3; i++){
        fscanf(fp,"%lf %lf %lf",&tv[i][1],&tv[i][2],&tv[i][3]);
      }
      if ( ! input_last("Atoms.Unitvectors>") ) {
        /* format error */
        printf("Format error for Atoms.Unitvectors\n");
        po++;
      }

      /* Ang to AU */
      if (unitvector_unit==0){
        for (i=1; i<=3; i++){
          tv[i][1] = tv[i][1]/BohrR;
          tv[i][2] = tv[i][2]/BohrR;
          tv[i][3] = tv[i][3]/BohrR;
        }
      }
    }

    else {
      /* automatically set up the unit cell */

      if (Solver==3){
	if (myid==Host_ID){
	  printf("You have to give the unit cell for the band calc.\n");
	}
	MPI_Finalize();
	exit(1);
      }

      Set_Cluster_UnitCell(tv,unitvector_unit);
      Determine_Cell_from_ECutoff(tv,Grid_Ecut+0.001);
    }

    /*  FRAC to AU */ 
    if (coordinates_unit==2){

      /* The fractional coordinates should be kept within 0 to 1. */

      for (i=1; i<=atomnum; i++){
        for (k=1; k<=3; k++){

          itmp = (int)Gxyz[i][k]; 

          if (1.0<Gxyz[i][k]){
            Gxyz[i][k] = Gxyz[i][k] - (double)itmp;

            if (myid==Host_ID){
              if (k==1) printf("The fractional coordinate of a-axis for atom %d was translated within the range (0 to 1).\n",i);
              if (k==2) printf("The fractional coordinate of b-axis for atom %d was translated within the range (0 to 1).\n",i);
              if (k==3) printf("The fractional coordinate of c-axis for atom %d was translated within the range (0 to 1).\n",i);
            }
	  }
          else if (Gxyz[i][k]<0.0){
            Gxyz[i][k] = Gxyz[i][k] + (double)(abs(itmp)+1);

            if (myid==Host_ID){
              if (k==1) printf("The fractional coordinate of a-axis for atom %d was translated within the range (0 to 1).\n",i);
              if (k==2) printf("The fractional coordinate of b-axis for atom %d was translated within the range (0 to 1).\n",i);
              if (k==3) printf("The fractional coordinate of c-axis for atom %d was translated within the range (0 to 1).\n",i);
            }
	  }

	}
      }

      /* calculation of xyz-coordinate in A.U. The grid origin is zero. */

      for (i=1; i<=atomnum; i++){
	x = Gxyz[i][1]*tv[1][1] + Gxyz[i][2]*tv[2][1] + Gxyz[i][3]*tv[3][1];
	y = Gxyz[i][1]*tv[1][2] + Gxyz[i][2]*tv[2][2] + Gxyz[i][3]*tv[3][2];
	z = Gxyz[i][1]*tv[1][3] + Gxyz[i][2]*tv[2][3] + Gxyz[i][3]*tv[3][3];
	Gxyz[i][1] = x;
	Gxyz[i][2] = y;
	Gxyz[i][3] = z;
      }
    }

  } /* if (Solver!=4){  */

  /*******************************
                NEGF 
  *******************************/

  else{

#ifdef TRAN
    TRAN_Input_std_Atoms(Solver);
#endif

  } /* else{ */ 

  /****************************************************
   set fixed atomic position in geometry optimization
   and MD:  

      1: fixed 
      0: relaxed
  ****************************************************/

  if (fp=input_find("<MD.Fixed.XYZ")) {

    for (i=1; i<=atomnum; i++){  
      fscanf(fp,"%d %d %d %d",
             &j,&atom_Fixed_XYZ[i][1],&atom_Fixed_XYZ[i][2],&atom_Fixed_XYZ[i][3]);
    }  

    if ( ! input_last("MD.Fixed.XYZ>") ) {
      /* format error */
      printf("Format error for MD.Fixed.XYZ\n");
      po++;
    }
  }

  /****************************************************
             set initial velocities for MD
  ****************************************************/

  MD_Init_Velocity = 0;

  if (fp=input_find("<MD.Init.Velocity")) {

    MD_Init_Velocity = 1;

    for (i=1; i<=atomnum; i++){  

      fscanf(fp,"%d %lf %lf %lf",&j,&Gxyz[i][24],&Gxyz[i][25],&Gxyz[i][26]);

      /***********************************************
          5.291772083*10^{-11} m / 2.418884*10^{-17} s 
          = 2.1876917*10^6 m/s                                                                               
          = 1 a.u. for velocity 

          1 m/s = 0.4571028 * 10^{-6} a.u.
      ***********************************************/
         
      for (j=1; j<=3; j++){        
	if (atom_Fixed_XYZ[i][j]==0){
	  Gxyz[i][23+j] = Gxyz[i][23+j]*0.4571028*0.000001;
	  Gxyz[i][26+j] = Gxyz[i][23+j];
	}

	else{
	  Gxyz[i][23+j] = 0.0;
	  Gxyz[i][26+j] = 0.0;
	}
      }

  printf("ABC0 %15.12f %15.12f %15.12f\n",Gxyz[1][24],Gxyz[1][25],Gxyz[1][26]);


    }  

    if ( ! input_last("MD.Init.Velocity>") ) {
      /* format error */
      printf("Format error for MD.Init.Velocity\n");
      po++;
    }
  }

  /****************************************************
         Starting point of  LDA+U    --- by MJ
  ****************************************************/

  /* for LDA+U */ 

  if (Hub_U_switch == 1){                              /* --- MJ */

    if (fp=input_find("<Hubbard.U.values")) {    

      /* initialize the U-values */ 
      for (i=0; i<SpeciesNum; i++){
	for (l=0; l<=Spe_MaxL_Basis[i]; l++){
	  for (mul=0; mul<Spe_Num_Basis[i][l]; mul++){
	    Hub_U_Basis[i][l][mul]=0.0 ;
	  }
	}
      }
	
      /* read the U-values from the '.dat' file  */    /* --- MJ */
      for (i=0; i<SpeciesNum; i++){
	fscanf(fp,"%s",Species);
        j = Species2int(Species);

	for (l=0; l<=Spe_MaxL_Basis[j]; l++){
	  for (mul=0; mul<Spe_Num_Basis[j][l]; mul++){
	    fscanf(fp,"%s %lf", buf, &Hub_U_Basis[j][l][mul]) ;
	  }
	}
      }

      if (! input_last("Hubbard.U.values>") ) {
	/* format error */
	printf("Format error for Hubbard.U.values\n");
	po++;
      }

    }   /*  if (fp=input_find("<Hubbard.U.values"))  */

    for (i=0; i<SpeciesNum; i++){
      for (l=0; l<=Spe_MaxL_Basis[i]; l++){
	for (mul=0; mul<Spe_Num_Basis[i][l]; mul++){
          Hub_U_Basis[i][l][mul] = Hub_U_Basis[i][l][mul]/eV2Hartree;
	}
      }
    }

  }   /*  if (Hub_U_switch == 1)  */ 

  /****************************************************
           Ending point of  LDA+U    --- by MJ
  ****************************************************/

  /****************************************************
           the maximum number of SCF iterations
  ****************************************************/

  input_int("scf.maxIter",&DFTSCF_loop,40);
  if (DFTSCF_loop<0) {
    printf("DFTSCF_loop should be over 0.\n");
    po++;
  }

  /****************************************************
                 Net charge of the system
  ****************************************************/

  input_double("Atoms.NetCharge",&Given_Total_Charge,(double)0.0);

  /*****************************************************
    if DM file exists, read data from it. default = off
  *****************************************************/

  input_logical("scf.restart",&Scf_RestartFromFile, 0); 

  /****************************************************
                    Band dispersion
  ****************************************************/

  input_logical("Band.dispersion",&Band_disp_switch, 0);

  Band_kpath=NULL;
  Band_kname=NULL;
  Band_N_perpath=NULL;
  input_int("Band.Nkpath",&Band_Nkpath,0);
  if (2<=level_stdout) printf("Band.Nkpath=%d\n",Band_Nkpath);

  if (Band_Nkpath>0) {

       Band_kPathUnit=0;
       if (fp=input_find("<Band.kpath.UnitCell") ) {   
          Band_kPathUnit=1;
          for (i=1;i<=3;i++) {
            fscanf(fp,"%lf %lf %lf",&Band_UnitCell[i][1],
                      &Band_UnitCell[i][2],&Band_UnitCell[i][3]);
          }
          if ( ! input_last("Band.kpath.UnitCell>") ) {
          /* format error */
           printf("Format error near Band.kpath.UnitCell>\n");
           po++;
          }
          for (i=1;i<=3;i++) {
            printf("%lf %lf %lf\n",Band_UnitCell[i][1],
                      Band_UnitCell[i][2],Band_UnitCell[i][3]);
          }
          

          if (unitvector_unit==0) {
            for (i=1;i<=3;i++)
              for (j=1;j<=3;j++)
               Band_UnitCell[i][j]=Band_UnitCell[i][j]/BohrR;
          }
       }
       else {
          for (i=1;i<=3;i++) 
            for (j=1;j<=3;j++) 
               Band_UnitCell[i][j]=tv[i][j];
       }

    /* allocate */

    Band_N_perpath=(int*)malloc(sizeof(int)*(Band_Nkpath+1));
    for (i=0; i<(Band_Nkpath+1); i++) Band_N_perpath[i] = 0;

    Band_kpath = (double***)malloc(sizeof(double**)*(Band_Nkpath+1));
    for (i=0; i<(Band_Nkpath+1); i++){
      Band_kpath[i] = (double**)malloc(sizeof(double*)*3);
      for (j=0; j<3; j++){
        Band_kpath[i][j] = (double*)malloc(sizeof(double)*4);
        for (k=0; k<4; k++) Band_kpath[i][j][k] = 0.0;
      }
    }

    Band_kname = (char***)malloc(sizeof(char**)*(Band_Nkpath+1));
    for (i=0; i<(Band_Nkpath+1); i++){
      Band_kname[i] = (char**)malloc(sizeof(char*)*3);
      for (j=0; j<3; j++){
        Band_kname[i][j] = (char*)malloc(sizeof(char)*YOUSO10);
      }
    }

    /* end of allocation */

    if (myid==Host_ID && 2<=level_stdout) printf("kpath\n");

    if (fp=input_find("<Band.kpath") ) {
      for (i=1; i<=Band_Nkpath; i++) {
        fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %s %s",
         &Band_N_perpath[i]  , 
         &Band_kpath[i][1][1], &Band_kpath[i][1][2],&Band_kpath[i][1][3],
         &Band_kpath[i][2][1], &Band_kpath[i][2][2],&Band_kpath[i][2][3],
         Band_kname[i][1],Band_kname[i][2]);

	if (myid==Host_ID && 2<=level_stdout){
          printf("%d (%lf %lf %lf) (%lf %lf %lf) %s %s\n",
          Band_N_perpath[i]  ,
          Band_kpath[i][1][1], Band_kpath[i][1][2],Band_kpath[i][1][3],
          Band_kpath[i][2][1], Band_kpath[i][2][2],Band_kpath[i][2][3],
          Band_kname[i][1],Band_kname[i][2]);
	}

      }
      if ( ! input_last("Band.kpath>") ) {
         /* format error */
         printf("Format error near Band.kpath>\n");
         po++;
      }
    }
    else {
      /* format error */
            printf("<Band.kpath is necessary.\n");
            po++;
    }
    if (Band_kPathUnit){
      kpath_changeunit( tv, Band_UnitCell, Band_Nkpath, Band_kpath );
    }
 
  }

  /****************************************************
                   One dimentional FFT
  ****************************************************/

  input_int("1DFFT.NumGridK",&Ngrid_NormK,900);

  if (Ngrid_NormK<GL_Mesh)
    List_YOUSO[15] = GL_Mesh + 10;
  else 
    List_YOUSO[15] = Ngrid_NormK + 2;

  input_double("1DFFT.EnergyCutoff",&ecutoff1dfft,(double)3600.0);
  PAO_Nkmax=sqrt(ecutoff1dfft);

  input_int("1DFFT.NumGridR",&OneD_Grid,900);

  /****************************************************
                   Orbital optimization
  ****************************************************/

  input_double("orbitalOpt.initPrefactor",&Cnt_scaling,(double)0.1);

  s_vec[0]="OFF";        s_vec[1]="Unrestricted"; 
  s_vec[2]="Restricted"; s_vec[3]="Species";
  i_vec[0]=0; i_vec[1]=1; i_vec[2]=2; i_vec[3]=3;
  input_string2int("orbitalOpt.Method",&orbitalopt,4,s_vec,i_vec);
  switch (orbitalopt) {
    case 0: { Cnt_switch=0; }                               break;
    case 1: { Cnt_switch=1; RCnt_switch=0; }                break;
    case 2: { Cnt_switch=1; RCnt_switch=1; }                break;
    case 3: { Cnt_switch=1; RCnt_switch=1; ACnt_switch=1; } break;
  }

  /************************************************************
    it it not allowed to make the simultaneuos specification
    of spin non-collinear and orbital optimization
  ************************************************************/

  if (SpinP_switch==3 && Cnt_switch==1){
    if (myid==Host_ID){
      printf("unsupported:\n");
      printf("simultaneuos specification of spin non-collinear\n");
      printf("and orbital optimization\n");
    }

    MPI_Finalize();
    exit(1);
  }

  /************************************************************
    it it not allowed to make the simultaneuos specification
    of NEGF and orbital optimization
  ************************************************************/

  if (Solver==4 && Cnt_switch==1){
    if (myid==Host_ID){
      printf("unsupported:\n");
      printf("simultaneuos specification of NEGF and orbital optimization\n");
    }

    MPI_Finalize();
    exit(1);
  }

  s_vec[0]="Symmetrical"; s_vec[1]="Free";
  i_vec[0]=1; i_vec[1]=0;
  input_string2int("orbitalOpt.InitCoes",&SICnt_switch,2,s_vec,i_vec);
  input_int("orbitalOpt.scf.maxIter",&orbitalOpt_SCF,12);
  input_int("orbitalOpt.MD.maxIter",&orbitalOpt_MD,5);
  input_int("orbitalOpt.per.MDIter",&orbitalOpt_per_MDIter,1000000);
  input_double("orbitalOpt.criterion",&orbitalOpt_criterion,(double)1.0e-4);

  /****************************************************
                  order-N method for SCF
  ****************************************************/

  input_double("orderN.HoppingRanges",&BCR,(double)5.0);
  BCR=BCR/BohrR;

  input_int("orderN.NumHoppings",&NOHS_L,2);
  if (Solver==2 || Solver==3 || Solver==4 || Solver==7){
    NOHS_L = 1;
    BCR = 1.0;
  }

  NOHS_C = NOHS_L;

  /* start Krylov */

  /* input_int("orderN.Kgrid",&orderN_Kgrid,5); */

  input_int("orderN.KrylovH.order",&KrylovH_order,400);
  input_int("orderN.KrylovS.order",&KrylovS_order,4*KrylovH_order);
  input_logical("orderN.Recalc.Buffer",&recalc_EM,0);
  input_logical("orderN.Exact.Inverse.S",&EKC_Exact_invS_flag,1);
  input_logical("orderN.Expand.Core",&EKC_expand_core_flag,1);
  input_logical("orderN.Inverse.S",&EKC_invS_flag,0); /* ?? */

  if (EKC_Exact_invS_flag==0){
    EKC_invS_flag = 1;
  }
  else {
    EKC_invS_flag = 0;
  }

  /* end Krylov */

  input_int("orderN.RecursiveLevels",&rlmax,10);
  List_YOUSO[3] = rlmax + 3;

  s_vec[0]="NO"; s_vec[1]="SQRT";
  i_vec[0]=1   ; i_vec[1]=2;
  input_string2int("orderN.TerminatorType",&T_switch,2, s_vec,i_vec);

  s_vec[0]="Recursion"; s_vec[1]="Divide";
  i_vec[0]=1; i_vec[1]=2;
  input_string2int("orderN.InverseSolver",&IS_switch,2,s_vec,i_vec);

  input_int("orderN.InvRecursiveLevels",&rlmax_IS,30);
  List_YOUSO[9] = rlmax_IS + 3;

  input_double("orderN.ChargeDeviation",&CN_Error,(double)1.0e-7);
  input_double("orderN.InitChemPot",&ChemP,(double)-1.0);
  ChemP = ChemP/eV2Hartree;

  input_int("orderN.AvNumTerminater",&Av_num,1);
  input_int("orderN.NumPoles",&POLES,300);

  /****************************************************
    control patameters for outputting wave functions
  ****************************************************/

  input_logical("MO.fileout",&MO_fileout,0);
  input_int("num.HOMOs",&num_HOMOs,1);
  input_int("num.LUMOs",&num_LUMOs,1);

  if (MO_fileout==0){
    num_HOMOs = 1;
    num_LUMOs = 1;
  }

  if ((Solver!=2 && Solver!=3 && Solver!=7) && MO_fileout==1){

    s_vec[0]="Recursion"; s_vec[1]="Cluster"; s_vec[2]="Band";
    s_vec[3]="NEGF";      s_vec[4]="DC";      s_vec[5]="GDC";
    s_vec[6]="Cluster2";  s_vec[7]="EC";

    printf("MO.fileout=ON is not supported in case of scf.EigenvalueSolver=%s\n",
           s_vec[Solver-1]);  
    printf("MO.fileout is changed to OFF\n");
    MO_fileout = 0;
  }

  if ( (Solver==2 || Solver==3 || Solver==7) && MO_fileout==1 ){
    List_YOUSO[31] = num_HOMOs;
    List_YOUSO[32] = num_LUMOs;
  }
  else{
    List_YOUSO[31] = 1;
    List_YOUSO[32] = 1;
  }

  /* for bulk */  

  input_int("MO.Nkpoint",&MO_Nkpoint,1);
  if (MO_Nkpoint<0){
    printf("MO_Nkpoint must be positive.\n");
    po++;
  }

  if ( (Solver==2 || Solver==3 || Solver==7) && MO_fileout==1 ){
    List_YOUSO[33] = MO_Nkpoint;
  }
  else{
    List_YOUSO[33] = 1;
  }

  /* memory allocation */
  Allocate_Arrays(5);

  if (fp=input_find("<MO.kpoint")) {
    for (i=0; i<MO_Nkpoint; i++){
      fscanf(fp,"%lf %lf %lf",&MO_kpoint[i][1],&MO_kpoint[i][2],&MO_kpoint[i][3]);

      if (2<=level_stdout){
        printf("<Input_std>  MO_kpoint %2d  %10.6f  %10.6f  %10.6f\n",
	       i,MO_kpoint[i][1],MO_kpoint[i][2],MO_kpoint[i][3]);
      }
    }
    if (!input_last("MO.kpoint>")) {
      /* format error */
      printf("Format error for MO.kpoint\n");
      po++;
    }

    if (Band_kPathUnit){
      kpoint_changeunit(tv,Band_UnitCell,MO_Nkpoint,MO_kpoint);
    }
  }

  /****************************************************
           for output of contracted orbitals    
  ****************************************************/

  input_logical("CntOrb.fileout",&CntOrb_fileout,0);
  if ((!(Cnt_switch==1 && RCnt_switch==1)) && CntOrb_fileout==1){
    printf("CntOrb.fileout=on is valid in case of orbitalOpt.Method=Restricted or Speciese\n");
    po++;    
  }

  if (CntOrb_fileout==1){
    input_int("Num.CntOrb.Atoms",&Num_CntOrb_Atoms,1);
    CntOrb_Atoms = (int*)malloc(sizeof(int)*Num_CntOrb_Atoms);      

    if (fp=input_find("<Atoms.Cont.Orbitals")) {
      for (i=0; i<Num_CntOrb_Atoms; i++){
        fscanf(fp,"%i",&CntOrb_Atoms[i]);
        if (CntOrb_Atoms[i]<1 || atomnum<CntOrb_Atoms[i]){
          printf("Invalid atom in <Atoms.Cont.Orbitals\n" ); 
          po++;
        }
      }

      if (!input_last("Atoms.Cont.Orbitals>")) {
        /* format error */
        printf("Format error for Atoms.Cont.Orbitals\n");
        po++;
      }
    }
  }

  /****************************************************
                 external electric field
  ****************************************************/

  r_vec[0]=0.0; r_vec[1]=0.0; r_vec[2]=0.0;
  input_doublev("scf.Electric.Field",3,E_Field,r_vec);

  if ( 1.0e-50<fabs(E_Field[0]) || 1.0e-50<fabs(E_Field[1]) || 1.0e-50<fabs(E_Field[2]) ){

    E_Field_switch = 1;    

    /*******************************************
                 unit transformation
         
        V/m = J/C/m
            = 1/(4.35975*10^{-18}) Hatree
	      *(1/(1.602177*10^{-19}) e )^{-1}
              *(1/(0.5291772*10^{-10}) a0 )^{-1}
            = 0.1944688 * 10^{-11} Hartree/e/a0

       input unit:  GV/m = 10^9 V/m
       used unit:   Hartree/e/a0

       GV/m = 0.1944688 * 10^{-2} Hartree/e/a0 
    *******************************************/

    length = sqrt( Dot_Product(tv[1], tv[1]) );
    x = E_Field[0]*tv[1][1]/length;
    y = E_Field[0]*tv[1][2]/length;
    z = E_Field[0]*tv[1][3]/length;

    length = sqrt( Dot_Product(tv[2], tv[2]) );
    x += E_Field[1]*tv[2][1]/length;
    y += E_Field[1]*tv[2][2]/length;
    z += E_Field[1]*tv[2][3]/length;

    length = sqrt( Dot_Product(tv[3], tv[3]) );
    x += E_Field[2]*tv[3][1]/length;
    y += E_Field[2]*tv[3][2]/length;
    z += E_Field[2]*tv[3][3]/length;

    length = sqrt( x*x + y*y + z*z );
    x = x/length;
    y = y/length;
    z = z/length;

    if (myid==Host_ID){  
      printf("\n");
      printf("<Applied External Electric Field>\n");
      printf("  direction (x, y, z)  = %10.5f %10.5f %10.5f\n",x,y,z);
      printf("  magnitude (10^9 V/m) = %10.5f\n\n",length);
    }

    E_Field[0] = 0.1944688*0.01*E_Field[0];
    E_Field[1] = 0.1944688*0.01*E_Field[1];
    E_Field[2] = 0.1944688*0.01*E_Field[2];
  }
  else {
    E_Field_switch = 0;    
  }

  /****************************************************
                      DOS, PDOS
  ****************************************************/

  input_logical("OpticalConductivity.fileout",&Opticalconductivity_fileout,0);   

  input_logical("Dos.fileout",&Dos_fileout,0);
  input_logical("DosGauss.fileout",&DosGauss_fileout,0);
  input_int("DosGauss.Num.Mesh",&DosGauss_Num_Mesh,200);
  input_double("DosGauss.Width",&DosGauss_Width,0.2);

  /* change the unit from eV to Hartree */
  DosGauss_Width = DosGauss_Width/eV2Hartree;

  if (Dos_fileout && DosGauss_fileout){

    if (myid==Host_ID){
      printf("Dos.fileout and DosGauss.fileout are mutually exclusive.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  if ( DosGauss_fileout && (Solver!=3 && PeriodicGamma_flag!=1) ){  /* band */

    if (myid==Host_ID){
      printf("DosGauss.fileout is supported for only band calculation.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  r_vec[0]=-20.0; r_vec[1]=20.0;
  input_doublev("Dos.Erange",2,Dos_Erange,r_vec);
  /* change the unit from eV to Hartree */
  Dos_Erange[0]= Dos_Erange[0]/eV2Hartree;
  Dos_Erange[1]= Dos_Erange[1]/eV2Hartree;

  i_vec[0]=Kspace_grid1; i_vec[1]=Kspace_grid2, i_vec[2]=Kspace_grid3;
  input_intv("Dos.Kgrid",3,Dos_Kgrid,i_vec);

  /****************************************************
   write a binary file, filename.scfout, which includes
   Hamiltonian, overlap, and density matrices.
  ****************************************************/
 
  input_logical("HS.fileout",&HS_fileout,0);

  /****************************************************
                   Voronoi charge
  ****************************************************/

  input_logical("Voronoi.charge",&Voronoi_Charge_flag,0);

  /****************************************************
                 Voronoi orbital moment
  ****************************************************/

  input_logical("Voronoi.orbital.moment",&Voronoi_OrbM_flag,0);

  /****************************************************
                       input_close
  ****************************************************/

  input_close();

  if (po>0 || input_errorCount()>0) {
    printf("errors in the inputfile\n");
    MPI_Finalize();
    exit(1);
  } 

  /****************************************************
           adjustment of atomic position
  ****************************************************/

  if (Solver!=4) Set_In_First_Cell();

  /****************************************************
                 Gxyz -> His_Gxyz[0]
  ****************************************************/

  k = 0;
  for (i=1; i<=atomnum; i++){
    His_Gxyz[0][k] = Gxyz[i][1]; k++;       
    His_Gxyz[0][k] = Gxyz[i][2]; k++;       
    His_Gxyz[0][k] = Gxyz[i][3]; k++;       
  }

  /****************************************************
          generate Monkhorst-Pack k-points 
  ****************************************************/

  if (Solver==3 && way_of_kpoint==2){

    Generating_MP_Special_Kpt(/* input */
			      atomnum, SpeciesNum,
			      tv, Gxyz,
			      InitN_USpin, InitN_DSpin,
			      Criterion_MP_Special_Kpt,
			      SpinP_switch,
			      WhatSpecies, 
			      Kspace_grid1, Kspace_grid2, Kspace_grid3
			      /* implicit output 
			      num_non_eq_kpt,
			      NE_KGrids1, NE_KGrids2, NE_KGrids3,
			      NE_T_k_op */ );
  }

  /****************************************************
                   print out to std
  ****************************************************/

  if (myid==Host_ID){  
    printf("\n\n<Input_std>  Your input file was normally read.\n");
    printf("<Input_std>  The system includes %i species and %i atoms.\n",
            real_SpeciesNum,atomnum);
  }

}




static void Set_In_First_Cell()
{
  int i,Gc_AN;
  int itmp;
  double Cxyz[4];
  double tmp[4];
  double xc,yc,zc;
  double CellV;

  /* calculate the reciprocal vectors */

  Cross_Product(tv[2],tv[3],tmp);
  CellV = Dot_Product(tv[1],tmp); 
  Cell_Volume = fabs(CellV);

  Cross_Product(tv[2],tv[3],tmp);
  rtv[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv[1][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv[3],tv[1],tmp);
  rtv[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv[2][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[1],tv[2],tmp);
  rtv[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv[3][3] = 2.0*PI*tmp[3]/CellV;

  /* find the center of coordinates */

  xc = 0.0;
  yc = 0.0;
  zc = 0.0;

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    xc += Gxyz[Gc_AN][1];
    yc += Gxyz[Gc_AN][2];
    zc += Gxyz[Gc_AN][3];
  }

  xc = xc/(double)atomnum;
  yc = yc/(double)atomnum;
  zc = zc/(double)atomnum;

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

    Cxyz[1] = Gxyz[Gc_AN][1] - xc;
    Cxyz[2] = Gxyz[Gc_AN][2] - yc;
    Cxyz[3] = Gxyz[Gc_AN][3] - zc;

    Cell_Gxyz[Gc_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
    Cell_Gxyz[Gc_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
    Cell_Gxyz[Gc_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

    for (i=1; i<=3; i++){
      if (1.0<fabs(Cell_Gxyz[Gc_AN][i])){
        if (0.0<=Cell_Gxyz[Gc_AN][i]){
          itmp = (int)Cell_Gxyz[Gc_AN][i]; 
          Cell_Gxyz[Gc_AN][i] = Cell_Gxyz[Gc_AN][i] - (double)itmp;
	}
        else{
          itmp = abs((int)Cell_Gxyz[Gc_AN][i]) + 1; 
          Cell_Gxyz[Gc_AN][i] = Cell_Gxyz[Gc_AN][i] + (double)itmp;
        }
      }
    }

    Gxyz[Gc_AN][1] =  Cell_Gxyz[Gc_AN][1]*tv[1][1]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][1]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][1] + xc;

    Gxyz[Gc_AN][2] =  Cell_Gxyz[Gc_AN][1]*tv[1][2]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][2]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][2] + yc;

    Gxyz[Gc_AN][3] =  Cell_Gxyz[Gc_AN][1]*tv[1][3]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][3]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][3] + zc;
  }

}





void SpeciesString2int(int p)
{
  int i,l,n,po;
  char c,cstr[YOUSO10*3];
  
  /* Get basis name */

  sprintf(SpeBasisName[p],"");

  i = 0;
  po = 0;
  while( ((c=SpeBasis[p][i])!='\0' || po==0) && i<YOUSO10 ){
    if (c=='-'){
      po = 1;
      SpeBasisName[p][i] = '\0';  
    }
    if (po==0) SpeBasisName[p][i] = SpeBasis[p][i]; 
    i++;
  }


  if (2<=level_stdout){
    printf("<Input_std>  SpeBasisName=%s\n",SpeBasisName[p]);
  }

  /* Get basis type */

  for (l=0; l<5; l++){
    Spe_Num_Basis[p][l] = 0;
  }

  i = 0;
  po = 0;

  while((c=SpeBasis[p][i])!='\0'){
    if (po==1){
      if      (c=='s'){ l=0; n=0; }
      else if (c=='p'){ l=1; n=0; }
      else if (c=='d'){ l=2; n=0; }
      else if (c=='f'){ l=3; n=0; }
      else if (c=='g'){ l=4; n=0; }
      else{

        if (n==0){
          cstr[0] = c;
          cstr[1] = '\0';
          Spe_Num_Basis[p][l]  = atoi(cstr);
          Spe_Num_CBasis[p][l] = atoi(cstr);
          n++;
        }
        else if (n==1){
          cstr[0] = c;
          cstr[1] = '\0';
          Spe_Num_CBasis[p][l] = atoi(cstr);   
          if (Spe_Num_Basis[p][l]<Spe_Num_CBasis[p][l]){
            printf("# of contracted orbitals are larger than # of primitive oribitals\n");
            MPI_Finalize();
            exit(1); 
          } 

          n++;
        }
        else {
          printf("Format error in Definition of Atomic Species\n");
          MPI_Finalize();
          exit(1);
	}
      } 
    }  

    if (SpeBasis[p][i]=='-') po = 1;
    i++;
  }

  for (l=0; l<5; l++){
    if (Spe_Num_Basis[p][l]!=0) Spe_MaxL_Basis[p] = l;

    if (2<=level_stdout){
      printf("<Input_std>  p=%2d l=%2d %2d %2d\n",
              p,l,Spe_Num_Basis[p][l],Spe_Num_CBasis[p][l]);
    }
  }
}



int Species2int(char Species[YOUSO10])
{
  int i,po;

  i = 0;
  po = 0; 
  while (i<SpeciesNum && po==0){
    if (SEQ(Species,SpeName[i])==1){
      po = 1;
    }
    if (po==0) i++;
  };

  if (po==1) return i;
  else {
    printf("Found an undefined species name %s\n",Species);
    printf("in Atoms.SpeciesAndCoordinates or Hubbard.U.values\n");

    printf("Please check your input file\n");
    MPI_Finalize();
    exit(1);
  }
}



int OrbPol2int(char OrbPol[YOUSO10])
{
  int i,po;
  char opns[3][YOUSO10]={"OFF","ON","EX"};

  i = 0;
  po = 0; 

  ToCapital(OrbPol);  

  while (i<3 && po==0){
    if (SEQ(OrbPol,opns[i])==1){
      po = 1;
    }
    if (po==0) i++;
  };

  if (po==1) return i;

  else {
    printf("Invalid flag for LDA+U (Atoms.SpeciesAndCoordinates)  %s\n",OrbPol);
    printf("Please check your input file\n");
    MPI_Finalize();
    exit(1);
  }
}



char *ToCapital(char *s)
{
  char *p;
  for (p=s; *p; p++) *p = toupper(*p);
  return (s);  
}




void kpath_changeunit( double tv[4][4], double tv0[4][4], int Band_Nkpath,
                       double ***Band_kpath )
{
  /***********************************************************************
    k1 rtv0[0] + k2 rtv0[1] + k3 rtv0[2] = l rtv[0] + m rtv[1] + n rtv[2] 
      rtv = reciptical vector of tv
      rtv0 = reciptical vector of tv0
    e.g.   l is given by 
     tv[0] ( k1 rtv0[0] + k2 rtv0[1] + k3 rtv0[2]) = l tv[0] rtv[0] 
  ************************************************************************/

  double tmp[4], CellV;
  double rtv[4][4],rtv0[4][4];
  int i,j;
  double r;
  
  Cross_Product(tv[2],tv[3],tmp);
  CellV = Dot_Product(tv[1],tmp); 
  
  Cross_Product(tv[2],tv[3],tmp);
  rtv[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv[1][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[3],tv[1],tmp);
  rtv[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv[2][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[1],tv[2],tmp);
  rtv[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv[3][3] = 2.0*PI*tmp[3]/CellV; 

  Cross_Product(tv0[2],tv0[3],tmp);
  CellV = Dot_Product(tv0[1],tmp);

  Cross_Product(tv0[2],tv0[3],tmp);
  rtv0[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv0[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv0[1][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv0[3],tv0[1],tmp);
  rtv0[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv0[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv0[2][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv0[1],tv0[2],tmp);
  rtv0[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv0[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv0[3][3] = 2.0*PI*tmp[3]/CellV;

  printf("kpath (converted)\n");
  for (i=1;i<=Band_Nkpath;i++) {
    for (j=1;j<=3;j++) tmp[j]=Band_kpath[i][1][j];
    for (j=1;j<=3;j++) {
      r =    tmp[1]* Dot_Product(tv[j],rtv0[1]) 
           + tmp[2]* Dot_Product(tv[j],rtv0[2])
	   + tmp[3]* Dot_Product(tv[j],rtv0[3]);
      Band_kpath[i][1][j] = r/ Dot_Product(tv[j],rtv[j]);
    }
    for (j=1;j<=3;j++) tmp[j]=Band_kpath[i][2][j];
    for (j=1;j<=3;j++) {
      r =    tmp[1]* Dot_Product(tv[j],rtv0[1])
             + tmp[2]* Dot_Product(tv[j],rtv0[2])
             + tmp[3]* Dot_Product(tv[j],rtv0[3]);
      Band_kpath[i][2][j] = r/ Dot_Product(tv[j],rtv[j]);
    }

    printf("(%lf %lf %lf) (%lf %lf %lf)\n",
           Band_kpath[i][1][1],Band_kpath[i][1][2],Band_kpath[i][1][3],
	   Band_kpath[i][2][1],Band_kpath[i][2][2],Band_kpath[i][2][3]);
  }   

}


void kpoint_changeunit(double tv[4][4],double tv0[4][4],int MO_Nkpoint,
                       double **MO_kpoint)
{
  /***********************************************************************
    k1 rtv0[0] + k2 rtv0[1] + k3 rtv0[2] = l rtv[0] + m rtv[1] + n rtv[2] 
      rtv = reciptical vector of tv
      rtv0 = reciptical vector of tv0
    e.g.   l is given by 
     tv[0] ( k1 rtv0[0] + k2 rtv0[1] + k3 rtv0[2]) = l tv[0] rtv[0] 
  ************************************************************************/

  double tmp[4], CellV;
  double rtv[4][4],rtv0[4][4];
  int i,j;
  double r;

  Cross_Product(tv[2],tv[3],tmp);
  CellV = Dot_Product(tv[1],tmp); 
  
  Cross_Product(tv[2],tv[3],tmp);
  rtv[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv[1][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[3],tv[1],tmp);
  rtv[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv[2][3] = 2.0*PI*tmp[3]/CellV;
  
  Cross_Product(tv[1],tv[2],tmp);
  rtv[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv[3][3] = 2.0*PI*tmp[3]/CellV; 

  Cross_Product(tv0[2],tv0[3],tmp);
  CellV = Dot_Product(tv0[1],tmp);

  Cross_Product(tv0[2],tv0[3],tmp);
  rtv0[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv0[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv0[1][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv0[3],tv0[1],tmp);
  rtv0[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv0[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv0[2][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv0[1],tv0[2],tmp);
  rtv0[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv0[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv0[3][3] = 2.0*PI*tmp[3]/CellV;

  printf("kpoint at which wave functions are calculated (converted)\n");

  for (i=0;i<MO_Nkpoint;i++){
    for (j=1;j<=3;j++) tmp[j]= MO_kpoint[i][j];
    for (j=1;j<=3;j++) {
      r =    tmp[1]* Dot_Product(tv[j],rtv0[1]) 
           + tmp[2]* Dot_Product(tv[j],rtv0[2])
           + tmp[3]* Dot_Product(tv[j],rtv0[3]);
      MO_kpoint[i][j] = r/ Dot_Product(tv[j],rtv[j]);
    }
    printf("%lf %lf %lf\n",MO_kpoint[i][1],MO_kpoint[i][2],MO_kpoint[i][3]);
  }
}




/* *** calculate an unit cell of a cluster,    ***
 * assuming that unit of Gxyz and Rc is A.U. 
 * cell size = (max[ xyz-Rc ] - min[ xyz+Rc ])*1.01
*/

void Set_Cluster_UnitCell(double tv[4][4], int unitflag)
{
/* 
 * Species: int SpeciesNum, char SpeName[], char SpeBasis[]
 *
 * Coordinate:   int WhatSpecies[]; double Gxyz[][]
 *
 * unitflag=0 (Ang.)  unitflag=1 (a.u.),  used only to print them
 *
 * tv[][] is always in a.u.
*/
  FILE *fp;
  int i,id,spe,myid; 
  double *sperc;
  double min[4],max[4],gmin[4],gmax[4]; 
  char FN_PAO[YOUSO10];
  char ExtPAO[YOUSO10] = ".pao";
  char DirPAO[YOUSO10];

  double margin=1.10;
  double unit;
  char *unitstr[2]={"Ang.","a.u."};

  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* set DirPAO */

  sprintf(DirPAO,"%s/PAO/",DFT_DATA_PATH);

  unit=1.0;
  if (unitflag==0) unit=BohrR;

  sperc=(double*)malloc(sizeof(double)*SpeciesNum);

  for (spe=0; spe<SpeciesNum; spe++){
    fnjoint2(DirPAO,SpeBasisName[spe],ExtPAO,FN_PAO);    
    if ((fp = fopen(FN_PAO,"r")) != NULL){
      input_open(FN_PAO);
      input_double("radial.cutoff.pao",&sperc[spe],(double)0.0);
      input_close();
      fclose(fp);
    }
    else {

      if (myid==Host_ID){ 
        printf("Set_Cluster_UnitCell: can not open %s\n", FN_PAO); 
      }

      MPI_Finalize();     
      exit(10); 
    }
  }

  if (level_stdout>=2) {
    for (spe=0;spe<SpeciesNum;spe++) {
    printf("<Set_Cluster_UnitCell> %d %s rc=%lf\n",spe,SpeName[spe],sperc[spe]);
    }
  }


  if (level_stdout>=2) {
   printf("<Set_Cluster_UnitCell> x y z   Rc\n");
  }

  for (i=1;i<=3;i++) {
     gmax[i]=gmin[i]=Gxyz[1][i];
  }

  for (i=1;i<=atomnum;i++) {
    id=WhatSpecies[i];   
    /*  printf("%d %d %lf\n",i,id,sperc[id]); */
    if (level_stdout>=2) {
      printf("<Set_Cluster_UnitCell> %lf %lf %lf %lf\n",
              Gxyz[i][1],Gxyz[i][2],Gxyz[i][3],sperc[id]);
    }
    min[1]=Gxyz[i][1]-sperc[id];
    min[2]=Gxyz[i][2]-sperc[id];
    min[3]=Gxyz[i][3]-sperc[id];
    max[1]=Gxyz[i][1]+sperc[id];
    max[2]=Gxyz[i][2]+sperc[id];
    max[3]=Gxyz[i][3]+sperc[id];
    if (min[1]<gmin[1]) gmin[1]=min[1];
    if (min[2]<gmin[2]) gmin[2]=min[2];
    if (min[3]<gmin[3]) gmin[3]=min[3];
    if (max[1]>gmax[1]) gmax[1]=max[1];
    if (max[2]>gmax[2]) gmax[2]=max[2];
    if (max[3]>gmax[3]) gmax[3]=max[3];
  }

  /* initialize */
  for (id=1;id<=3;id++){
    for(i=1;i<=3;i++) {
      tv[id][i]=0.0;
    }
  }

  tv[1][1]=(gmax[1]-gmin[1])*margin;
  tv[2][2]=(gmax[2]-gmin[2])*margin;
  tv[3][3]=(gmax[3]-gmin[3])*margin;

  if (level_stdout>=2) {
    printf("<Set_Cluster_UnitCell> to determine the unit cell, min and max includes effects of Rc\n");
    for (i=1;i<=3;i++) 
    printf("<Set_Cluster_UnitCell> axis=%d min,max=%lf %lf\n", i, gmin[i]*unit,gmax[i]*unit);
  } 

  if (myid==Host_ID){

    printf("<Set_Cluster_UnitCell> automatically determied UnitCell(%s)\n<Set_Cluster_UnitCell> from atomic positions and Rc of PAOs (margin= %5.2lf%%)\n",unitstr[unitflag],(margin-1.0)*100.0);
    for(i=1;i<=3;i++) {
      printf("<Set_Cluster_UnitCell> %lf %lf %lf\n",
               tv[i][1]*unit,tv[i][2]*unit, tv[i][3]*unit);
    }
    printf("\n");
  }

  free(sperc); 
}



int divisible_cheker(int N)
{
  /************************
   return 0; non-divisible 
   return 1; divisible 
  ************************/

  int i,po;

  if (N!=1){
    po = 1; 
    for (i=0; i<NfundamentalNum; i++){
      if ( N!=1 && (N % fundamentalNum[i])==0 ){
	po = 0;
	N = N/fundamentalNum[i];
      }
    }
  }
  else{
    po = 0;
  }

  if (po==0 && N!=1){
    divisible_cheker(N);
  }

  if (po==0) return 1;
  else       return 0;
}




void Setup_Mixed_Basis(char *file, int myid)
{
  int i,po;
  FILE *fp;
  char *s_vec[20];
  int i_vec[20];
  double LngA,LngB,LngC;
  double A2,B2,C2;
  double GridV,GVolume;
  double tmp[4];

  po = 0;

  s_vec[0]="Ang"; s_vec[1]="AU";
  i_vec[0]=0;  i_vec[1]=1;
  input_string2int("Atoms.UnitVectors.Unit",&unitvector_unit,2,s_vec,i_vec);

  if (fp=input_find("<Atoms.Unitvectors")) {
    for (i=1; i<=3; i++){
      fscanf(fp,"%lf %lf %lf",&tv[i][1],&tv[i][2],&tv[i][3]);
    }
    if ( ! input_last("Atoms.Unitvectors>") ) {
      /* format error */
      printf("Format error for Atoms.Unitvectors\n");
      po++;
    }

    /* Ang to AU */
    if (unitvector_unit==0){
      for (i=1; i<=3; i++){
	tv[i][1] = tv[i][1]/BohrR;
	tv[i][2] = tv[i][2]/BohrR;
	tv[i][3] = tv[i][3]/BohrR;
      }
    }
  }

  else {
    /* for automatically set up the unit cell */

    if (myid==Host_ID){
      printf("You have to give the unit cell for the mixed basis scheme.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /***************************
     set up:
     Finite_Elements_Ecut
     Ngrid1_FE
     Ngrid2_FE
     Ngrid3_FE
  ***************************/

  Finite_Elements_Ecut = Grid_Ecut/(double)(NG_Mixed_Basis*NG_Mixed_Basis);

  printf("Finite_Elements_Ecut=%15.12f\n",Finite_Elements_Ecut);  

  Ngrid1_FE = 4;
  Ngrid2_FE = 4;
  Ngrid3_FE = 4;

  gtv_FE[1][1] = tv[1][1]/((double)Ngrid1_FE);
  gtv_FE[1][2] = tv[1][2]/((double)Ngrid1_FE);
  gtv_FE[1][3] = tv[1][3]/((double)Ngrid1_FE);

  gtv_FE[2][1] = tv[2][1]/((double)Ngrid2_FE);
  gtv_FE[2][2] = tv[2][2]/((double)Ngrid2_FE);
  gtv_FE[2][3] = tv[2][3]/((double)Ngrid2_FE);
    
  gtv_FE[3][1] = tv[3][1]/((double)Ngrid3_FE);
  gtv_FE[3][2] = tv[3][2]/((double)Ngrid3_FE);
  gtv_FE[3][3] = tv[3][3]/((double)Ngrid3_FE);
    
  Cross_Product(gtv_FE[2],gtv_FE[3],tmp);

  GridV = Dot_Product(gtv_FE[1],tmp); 
  GVolume = fabs( GridV );

  Cross_Product(gtv_FE[2],gtv_FE[3],tmp);
  rgtv_FE[1][1] = 2.0*PI*tmp[1]/GridV;
  rgtv_FE[1][2] = 2.0*PI*tmp[2]/GridV;
  rgtv_FE[1][3] = 2.0*PI*tmp[3]/GridV;

  Cross_Product(gtv_FE[3],gtv_FE[1],tmp);
  rgtv_FE[2][1] = 2.0*PI*tmp[1]/GridV;
  rgtv_FE[2][2] = 2.0*PI*tmp[2]/GridV;
  rgtv_FE[2][3] = 2.0*PI*tmp[3]/GridV;
    
  Cross_Product(gtv_FE[1],gtv_FE[2],tmp);
  rgtv_FE[3][1] = 2.0*PI*tmp[1]/GridV;
  rgtv_FE[3][2] = 2.0*PI*tmp[2]/GridV;
  rgtv_FE[3][3] = 2.0*PI*tmp[3]/GridV;

  A2 = rgtv_FE[1][1]*rgtv_FE[1][1] + rgtv_FE[1][2]*rgtv_FE[1][2] + rgtv_FE[1][3]*rgtv_FE[1][3];
  B2 = rgtv_FE[2][1]*rgtv_FE[2][1] + rgtv_FE[2][2]*rgtv_FE[2][2] + rgtv_FE[2][3]*rgtv_FE[2][3];
  C2 = rgtv_FE[3][1]*rgtv_FE[3][1] + rgtv_FE[3][2]*rgtv_FE[3][2] + rgtv_FE[3][3]*rgtv_FE[3][3];

  A2 = A2/4.0;
  B2 = B2/4.0;
  C2 = C2/4.0;

  Find_ApproxFactN(tv,&Finite_Elements_Ecut,&Ngrid1_FE,&Ngrid2_FE,&Ngrid3_FE,&A2,&B2,&C2);


  /*
  Ngrid1_FE = 2;
  Ngrid2_FE = 2;
  Ngrid3_FE = 2;
  */



  printf("Ngrid1_FE=%2d A2=%15.12f\n",Ngrid1_FE,A2);  
  printf("Ngrid2_FE=%2d B2=%15.12f\n",Ngrid2_FE,B2);  
  printf("Ngrid3_FE=%2d C2=%15.12f\n",Ngrid3_FE,C2);

  gtv_FE[1][1] = tv[1][1]/(double)Ngrid1_FE;
  gtv_FE[1][2] = tv[1][2]/(double)Ngrid1_FE;
  gtv_FE[1][3] = tv[1][3]/(double)Ngrid1_FE;

  gtv_FE[2][1] = tv[2][1]/(double)Ngrid2_FE;
  gtv_FE[2][2] = tv[2][2]/(double)Ngrid2_FE;
  gtv_FE[2][3] = tv[2][3]/(double)Ngrid2_FE;
    
  gtv_FE[3][1] = tv[3][1]/(double)Ngrid3_FE;
  gtv_FE[3][2] = tv[3][2]/(double)Ngrid3_FE;
  gtv_FE[3][3] = tv[3][3]/(double)Ngrid3_FE;

  LngA = sqrt( Dot_Product(gtv_FE[1],gtv_FE[1]) );
  LngB = sqrt( Dot_Product(gtv_FE[2],gtv_FE[2]) );
  LngC = sqrt( Dot_Product(gtv_FE[3],gtv_FE[3]) );

  rcut_FEB = LngA;
  if (LngB<rcut_FEB) rcut_FEB = LngB;
  if (LngC<rcut_FEB) rcut_FEB = LngC;
  rcut_FEB -= 1.0e-10;

  printf("rcut_FEB=%15.12f\n",rcut_FEB);

  /* set up Grid_Ecut */

  Ngrid1 = Ngrid1_FE*NG_Mixed_Basis;
  Ngrid2 = Ngrid2_FE*NG_Mixed_Basis;
  Ngrid3 = Ngrid3_FE*NG_Mixed_Basis;

  gtv[1][1] = tv[1][1]/((double)Ngrid1);
  gtv[1][2] = tv[1][2]/((double)Ngrid1);
  gtv[1][3] = tv[1][3]/((double)Ngrid1);

  gtv[2][1] = tv[2][1]/((double)Ngrid2);
  gtv[2][2] = tv[2][2]/((double)Ngrid2);
  gtv[2][3] = tv[2][3]/((double)Ngrid2);
    
  gtv[3][1] = tv[3][1]/((double)Ngrid3);
  gtv[3][2] = tv[3][2]/((double)Ngrid3);
  gtv[3][3] = tv[3][3]/((double)Ngrid3);
    
  Cross_Product(gtv[2],gtv[3],tmp);

  GridV = Dot_Product(gtv[1],tmp); 
  GVolume = fabs( GridV );

  Cross_Product(gtv[2],gtv[3],tmp);
  rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

  Cross_Product(gtv[3],gtv[1],tmp);
  rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[2][3] = 2.0*PI*tmp[3]/GridV;
    
  Cross_Product(gtv[1],gtv[2],tmp);
  rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

  A2 = rgtv[1][1]*rgtv[1][1] + rgtv[1][2]*rgtv[1][2] + rgtv[1][3]*rgtv[1][3];
  B2 = rgtv[2][1]*rgtv[2][1] + rgtv[2][2]*rgtv[2][2] + rgtv[2][3]*rgtv[2][3];
  C2 = rgtv[3][1]*rgtv[3][1] + rgtv[3][2]*rgtv[3][2] + rgtv[3][3]*rgtv[3][3];

  A2 = A2/4.0;
  B2 = B2/4.0;
  C2 = C2/4.0;

  Grid_Ecut = (A2 + B2 + C2 - 1.0)/3.0;
}






