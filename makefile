###################################################################
#                                                                 #
#  Please set a proper CC and LIB for the compilation.            #
#  Examples of CC and LIB on several platforms are shown below.   #
#                                                                 #
###################################################################

#
#  Cygwin, IBM ThinkPAD X40 (Pentium M 1.0GHz)
#
# CC      = gcc -Dnompi -Dblaswrap -O3 -I/home/ozaki/include
# STACK   = -Wl,--heap,9000000,--stack,9000000  
# LIB     =-L/home/ozaki/lib -llapack -lblas -lg2c -lI77 -lfftw3 -static
#

#
# Cray-XT3 (Opteron 2.4GHz)
#
# CC      = mpicc -Dxt3 -tp k8-64 -O3 -mcmodel=medium -I/work/t-ozaki/include
# LIB     =-L/opt/gcc-catamount/3.3/lib64 -lg2c -L/work/t-ozaki/lib -lacml -L/ufs/home/fs3024/t-ozaki/XT3/lib -lfftw3
#

# 
# Altix 4700
#
# CC = icc -openmp -O2 -I/opt/intel/mkl/10.0.1.014/include/fftw
# LIB= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lfftw3xc_intel -lmpi
#

# 
# P32 (Opteron 2.0GHz)
#
# CC       = mpicc -tp k8-64 -fast -mcmodel=medium -I/home/t-ozaki/include
# LIB      = -L/home/t-ozaki/lib -lfftw3 /opt/acml3.1.0/gnu64/lib/libacml.a -L/usr/lib64/ -lg2c
#

#
# mx73 (Xeon)
#
# mpicc = /usr/local/mpich-1.2.5-intel/bin/mpicc
# CC      = $(mpicc) -O3 -DTRAN -I/usr/local/include -I/home/ozaki/include
# LIB     = -L/home/ozaki/lib -latlas_p4 -lfftw3 -static
#

#
# abacus2 (Opteron 2.6 GHz)
#
# CC=/usr/local/mpich-1.2.7p1/bin/mpicc -tp amd64e -O3 -mcmodel=medium -I/usr/local/fftw3/include
# LIB= -L/usr/local/fftw3/lib -lfftw3 /usr/local/acml/gnu64/lib/libacml.a /usr/lib64/libg2c.a
#

#
# lauda (Pentium 4)
#
# CC      = mpicc -O3 -DTRAN -Dblaswrap -I/home/ozaki/include
# LIB     = -L/home/ozaki/lib -lfftw3 -llapack -lblas -lg2c -lI77 -static
#

#
# manfath1 (Pentium D)
#
# CC      = mpicc -O1 -I/home/ozaki/include -I/usr/include
# LIB     = /usr/lib64/liblapack.a /usr/lib64/libblas.a /usr/lib/gcc-lib/x86_64-redhat-linux/3.2.3/libg2c.a -lfftw3 -static
#

#
# pcc00 (AMD Opteron DP Model 250(2.4GHz) x 32)
#
# CC       = /opt/mpich/mpich-pgi-1.2.7/bin/mpicc -tp k8-64 -fast -mcmodel=medium -I/work/t-ozaki/include
# LIB      = -L/work/t-ozaki/lib -lfftw3 /work/t-ozaki/acml/gnu64/lib/libacml.a /opt/gcc33/lib64/libg2c.a
#
# for parallel running 
#  /opt/mpich/mpich-pgi-1.2.7/bin/mpirun openmx Methane.dat
#

# 
# Sun Fire V890 (ltraSPARC-IV(1.35GHz)x8CPU, 64G memory)
#
# CC      = gcc -O3 -Dnompi -Dblaswrap -I/work/t-ozaki/include
# LIB     = -L/work/t-ozaki/lib -lfftw3 -llapack -lblas -lg2c -lI77
#

#
# Hitachi SR11K (IBM Power4+) in IMS
#
# Compilation on SR11k is tricky.
# After compiling sources using mpcc, then link using mpif90
# CC	= mpif90 -64 -Os -noparallel -I/home/o000/nak1/SR11K_cc64/include
# CC	= mpcc -O3 -q64 -Df77 -bnoquiet -qstrict -I/home/o000/nak1/SR11K_cc64/include
# LIB     = /usr/local/lib/liblapack.a /usr/local/lib/libblas.a /usr/local/lib/libg2c.a /home/o000/nak1/SR11K_cc64/lib/libfftw3.a
#
# FFTW should be also compiled by the same compiler as in compilation of openmx in 64bit mode.
# compilation of FFTW Ver.3.0.1 in 64 bit mode can be done as follows:
#
#  (1) In a file "configure", change as follows:
#      test -z "$AR_FLAGS" && AR_FLAGS="-X64 cru"
#  (2) ./configure CC="cc" CFLAGS="-O3 -qstrict -q64" AR_FLAGS="-X64 cru" --prefix=/home/o000/nak1/SR11K_cc64/
#  (3) In a file "libbench2/Makefile", change as follows:
#      libbench2_a_AR = $(AR) -X64 cru
#  (4) make install
#

#
# HA8000 (Xeon in IMS)
#
# CC      = mpicc -nocheckpoint -compiler intel7 -Dblaswrap -O3 -I/home/o000/nak1/HA_gcc/include  
# LIB     = -L/home/o000/nak1/HA_gcc/lib -lfftw3 -llapack -lblas -lg2c -lI77 -static
#


CC      = gcc -O3 -Dnompi -Dnoomp
LIB     = -lfftw3 -llapack -lblas -lgfortran


CFLAGS  = -g 

OBJS    = openmx.o openmx_common.o Input_std.o Inputtools.o \
          init.o LU_inverse.o ReLU_inverse.o \
          truncation.o Find_ApproxFactN.o readfile.o FT_PAO.o FT_NLP.o \
          FT_ProExpn_VNA.o FT_VNA.o FT_ProductPAO.o \
          Hamiltonian_Cluster.o Overlap_Cluster.o Hamiltonian_Band.o \
          Overlap_Band.o Hamiltonian_Cluster_NC.o Hamiltonian_Band_NC.o \
          Hamiltonian_Cluster_SO.o Get_OneD_HS_Col.o SetPara_DFT.o \
          XC_Ceperly_Alder.o XC_CA_LSDA.o XC_PW91C.o XC_PBE.o XC_EX.o \
          DFT.o Mixing_DM.o Force.o Poisson.o \
          Cluster_DFT.o Cluster_DFT_Dosout.o \
          Band_DFT_Col.o  Band_DFT_NonCol.o \
          Band_DFT_kpath.o Band_DFT_MO.o Band_DFT_Dosout.o \
          Set_Density_Grid.o Set_Orbitals_Grid.o Set_Aden_Grid.o \
          Gauss_Legendre.o zero_cfrac.o xyz2spherical.o AngularF.o \
          RadialF.o Dr_RadialF.o PhiF.o  VNAF.o Dr_VNAF.o VH_AtomF.o \
          Dr_VH_AtomF.o Dr_AtomicDenF.o RF_BesselF.o QuickSort.o \
          Nonlocal_RadialF.o AtomicDenF.o Mulliken_Charge.o \
          Occupation_Number_LDA_U.o Eff_Hub_Pot.o \
          EulerAngle_Spin.o Smoothing_Func.o Orbital_Moment.o \
          AtomicPCCF.o Dr_AtomicPCCF.o Pot_NeutralAtom.o \
          Simple_Mixing_DM.o DIIS_Mixing_DM.o GR_Pulay_DM.o \
          Kerker_Mixing_Rhok.o DIIS_Mixing_Rhok.o \
          Total_Energy.o Contract_Hamiltonian.o Contract_iHNL.o \
          Cont_Matrix0.o Cont_Matrix1.o Cont_Matrix2.o Cont_Matrix3.o \
          Opt_Contraction.o Initial_CntCoes.o Set_XC_Grid.o \
          Get_Orbitals.o Get_dOrbitals.o Get_Cnt_Orbitals.o \
          Get_Cnt_dOrbitals.o Gaunt.o Find_CGrids.o MD_pac.o \
          RestartFileDFT.o Output_CompTime.o Merge_LogFile.o Make_FracCoord.o \
          Make_InputFile_with_FinalCoord.o \
          Divide_Conquer.o GDivide_Conquer.o Krylov.o Krylov_Dosout.o \
          Divide_Conquer_Dosout.o GDivide_Conquer_Dosout.o \
          Eigen_lapack.o EigenBand_lapack.o \
          Eigen_PReHH.o BroadCast_ReMatrix.o \
          Eigen_PHH.o BroadCast_ComplexMatrix.o \
          lapack_dstedc1.o lapack_dstedc2.o lapack_dstegr1.o lapack_dstegr2.o \
          lapack_dstevx1.o lapack_dstevx2.o lapack_dsteqr1.o \
          Nonlocal_Basis.o Set_OLP_Kin.o Set_Nonlocal.o Set_ProExpn_VNA.o \
          Set_Hamiltonian.o Set_Vpot.o \
          Voronoi_Charge.o Voronoi_Orbital_Moment.o Fuzzy_Weight.o \
          dampingF.o deri_dampingF.o Spherical_Bessel.o \
          iterout.o iterout_md.o Allocate_Arrays.o Free_Arrays.o \
          Init_List_YOUSO.o outputfile1.o \
          malloc_multidimarray.o PrintMemory.o PrintMemory_Fix.o \
          dtime.o OutData.o init_alloc_first.o File_CntCoes.o \
          SCF2File.o mimic_mpi.o mimic_omp.o Make_Comm_Worlds.o \
          Set_Allocate_Atom2CPU.o Cutoff.o Generating_MP_Special_Kpt.o \
          setup_CPU_group.o Maketest.o Runtest.o Memory_Leak_test.o Force_test.o Show_DFT_DATA.o \
          TRAN_Allocate.o TRAN_DFT.o TRAN_DFT_Dosout.o TRAN_Apply_Bias2e.o \
          TRAN_FFT_interpolation3c.o TRAN_RestartFile.o \
          TRAN_Calc_CentGreen.o TRAN_Input_std.o TRAN_Set_CentOverlap.o \
          TRAN_Calc_CentGreenLesser.o TRAN_Calc_GC_LorR.o TRAN_Input_std_Atoms.o TRAN_Set_Electrode_Grid.o \
          TRAN_Calc_GridBound.o TRAN_Interp_ElectrodeDensity_Grid.o TRAN_Set_IntegPath.o \
          TRAN_Output_HKS.o TRAN_Set_MP.o TRAN_Calc_SelfEnergy.o \
          TRAN_Output_Trans_HS.o TRAN_Set_PathEnergyStr.o \
          TRAN_Calc_SurfGreen.o TRAN_Output_Transmission.o TRAN_Set_SurfOverlap.o \
          TRAN_Overwrite_Densitygrid.o TRAN_Set_Value.o TRAN_Calc_VHartree_G0.o \
          TRAN_Overwrite_V2.o TRAN_adjust_Grid_Origin.o TRAN_Calc_VHartree_Gnon0.o \
          TRAN_Poisson.o TRAN_adjust_Ngrid.o TRAN_Print.o \
          TRAN_Print_Grid.o Lapack_LU_inverse.o TRAN_Distribute_Node.o \
          TRAN_Output_HKS_Write_Grid.o TRAN_Credit.o \
          TRAN_Check_Region_Lead.o TRAN_Check_Region.o	\
          RecursionS_H.o 

#         RecursionS_B.o RecursionS_C.o RecursionS_E.o \
#         RecursionS_F.o RecursionS_G.o RecursionS_H.o RecursionS_I.o \
#         IS_Lanczos.o IS_Taylor.o IS_Hotelling.o IS_LU.o \

# PROG    = openmx.exe
# PROG    = openmx

PROG    = openmx
DESTDIR = ../work

All:	$(PROG)

openmx:	$(OBJS)
	$(CC) $(OBJS) $(STACK) $(LIB) -lm -o openmx

openmx.o: openmx.c openmx_common.h
	$(CC) -c openmx.c
openmx_common.o: openmx_common.c openmx_common.h
	$(CC) -c openmx_common.c
Input_std.o: Input_std.c openmx_common.h Inputtools.h
	$(CC) -c Input_std.c
Inputtools.o: Inputtools.c
	$(CC) -c Inputtools.c

init.o: init.c openmx_common.h
	$(CC) -c init.c
LU_inverse.o: LU_inverse.c openmx_common.h
	$(CC) -c LU_inverse.c
ReLU_inverse.o: ReLU_inverse.c openmx_common.h
	$(CC) -c ReLU_inverse.c
truncation.o: truncation.c openmx_common.h
	$(CC) -c truncation.c
Find_ApproxFactN.o: Find_ApproxFactN.c openmx_common.h
	$(CC) -c Find_ApproxFactN.c
Find_CGrids.o: Find_CGrids.c openmx_common.h
	$(CC) -c Find_CGrids.c
readfile.o: readfile.c openmx_common.h
	$(CC) -c readfile.c
#
#
#
Hamiltonian_Cluster.o: Hamiltonian_Cluster.c openmx_common.h
	$(CC) -c Hamiltonian_Cluster.c
Overlap_Cluster.o: Overlap_Cluster.c openmx_common.h
	$(CC) -c Overlap_Cluster.c
Hamiltonian_Band.o: Hamiltonian_Band.c openmx_common.h
	$(CC) -c Hamiltonian_Band.c
Overlap_Band.o: Overlap_Band.c openmx_common.h
	$(CC) -c Overlap_Band.c
Hamiltonian_Cluster_NC.o: Hamiltonian_Cluster_NC.c openmx_common.h
	$(CC) -c Hamiltonian_Cluster_NC.c
Hamiltonian_Cluster_SO.o: Hamiltonian_Cluster_SO.c openmx_common.h
	$(CC) -c Hamiltonian_Cluster_SO.c
Hamiltonian_Band_NC.o: Hamiltonian_Band_NC.c openmx_common.h
	$(CC) -c Hamiltonian_Band_NC.c
Get_OneD_HS_Col.o: Get_OneD_HS_Col.c openmx_common.h
	$(CC) -c Get_OneD_HS_Col.c
#
#
#  
SetPara_DFT.o: SetPara_DFT.c openmx_common.h
	$(CC) -c SetPara_DFT.c
XC_Ceperly_Alder.o: XC_Ceperly_Alder.c openmx_common.h
	$(CC) -c XC_Ceperly_Alder.c
XC_CA_LSDA.o: XC_CA_LSDA.c openmx_common.h
	$(CC) -c XC_CA_LSDA.c
XC_PW91C.o: XC_PW91C.c openmx_common.h
	$(CC) -c XC_PW91C.c
XC_PBE.o: XC_PBE.c openmx_common.h
	$(CC) -c XC_PBE.c
XC_EX.o: XC_EX.c openmx_common.h
	$(CC) -c XC_EX.c
#
# SCF
#
DFT.o: DFT.c openmx_common.h
	$(CC) -c DFT.c
Cluster_DFT.o: Cluster_DFT.c openmx_common.h
	$(CC) -c Cluster_DFT.c
Cluster_DFT_Dosout.o: Cluster_DFT_Dosout.c openmx_common.h
	$(CC) -c Cluster_DFT_Dosout.c
Band_DFT_Col.o: Band_DFT_Col.c openmx_common.h
	$(CC) -c Band_DFT_Col.c
Band_DFT_NonCol.o: Band_DFT_NonCol.c openmx_common.h
	$(CC) -c Band_DFT_NonCol.c
Band_DFT_kpath.o: Band_DFT_kpath.c openmx_common.h
	$(CC) -c Band_DFT_kpath.c
Band_DFT_MO.o: Band_DFT_MO.c openmx_common.h
	$(CC) -c Band_DFT_MO.c
Band_DFT_Dosout.o: Band_DFT_Dosout.c openmx_common.h
	$(CC) -c Band_DFT_Dosout.c
Mixing_DM.o: Mixing_DM.c openmx_common.h
	$(CC) -c Mixing_DM.c
Force.o: Force.c openmx_common.h
	$(CC) -c Force.c
Poisson.o: Poisson.c openmx_common.h
	$(CC) -c Poisson.c
Mulliken_Charge.o: Mulliken_Charge.c openmx_common.h
	$(CC) -c Mulliken_Charge.c
Occupation_Number_LDA_U.o: Occupation_Number_LDA_U.c openmx_common.h
	$(CC) -c Occupation_Number_LDA_U.c
Eff_Hub_Pot.o: Eff_Hub_Pot.c openmx_common.h
	$(CC) -c Eff_Hub_Pot.c
EulerAngle_Spin.o: EulerAngle_Spin.c openmx_common.h
	$(CC) -c EulerAngle_Spin.c
Orbital_Moment.o: Orbital_Moment.c openmx_common.h
	$(CC) -c Orbital_Moment.c
Smoothing_Func.o: Smoothing_Func.c openmx_common.h
	$(CC) -c Smoothing_Func.c
Gauss_Legendre.o: Gauss_Legendre.c openmx_common.h
	$(CC) -c Gauss_Legendre.c
zero_cfrac.o: zero_cfrac.c openmx_common.h
	$(CC) -c zero_cfrac.c
xyz2spherical.o: xyz2spherical.c openmx_common.h
	$(CC) -c xyz2spherical.c
AngularF.o: AngularF.c openmx_common.h
	$(CC) -c AngularF.c
RadialF.o: RadialF.c openmx_common.h
	$(CC) -c RadialF.c
Dr_RadialF.o: Dr_RadialF.c openmx_common.h
	$(CC) -c Dr_RadialF.c
PhiF.o: PhiF.c openmx_common.h
	$(CC) -c PhiF.c
VNAF.o: VNAF.c openmx_common.h
	$(CC) -c VNAF.c
Dr_VNAF.o: Dr_VNAF.c openmx_common.h
	$(CC) -c Dr_VNAF.c
VH_AtomF.o: VH_AtomF.c openmx_common.h
	$(CC) -c VH_AtomF.c
Dr_VH_AtomF.o: Dr_VH_AtomF.c openmx_common.h
	$(CC) -c Dr_VH_AtomF.c

RF_BesselF.o: RF_BesselF.c openmx_common.h
	$(CC) -c RF_BesselF.c
Nonlocal_RadialF.o: Nonlocal_RadialF.c openmx_common.h
	$(CC) -c Nonlocal_RadialF.c

Set_Orbitals_Grid.o: Set_Orbitals_Grid.c openmx_common.h
	$(CC) -c Set_Orbitals_Grid.c
Set_Density_Grid.o: Set_Density_Grid.c openmx_common.h
	$(CC) -c Set_Density_Grid.c
Set_Aden_Grid.o: Set_Aden_Grid.c openmx_common.h
	$(CC) -c Set_Aden_Grid.c

AtomicDenF.o: AtomicDenF.c openmx_common.h
	$(CC) -c AtomicDenF.c
Dr_AtomicDenF.o: Dr_AtomicDenF.c openmx_common.h
	$(CC) -c Dr_AtomicDenF.c
AtomicPCCF.o: AtomicPCCF.c openmx_common.h
	$(CC) -c AtomicPCCF.c
Dr_AtomicPCCF.o: Dr_AtomicPCCF.c openmx_common.h
	$(CC) -c Dr_AtomicPCCF.c
Pot_NeutralAtom.o: Pot_NeutralAtom.c openmx_common.h
	$(CC) -c Pot_NeutralAtom.c
Simple_Mixing_DM.o: Simple_Mixing_DM.c openmx_common.h
	$(CC) -c Simple_Mixing_DM.c
DIIS_Mixing_DM.o: DIIS_Mixing_DM.c openmx_common.h
	$(CC) -c DIIS_Mixing_DM.c
GR_Pulay_DM.o: GR_Pulay_DM.c openmx_common.h
	$(CC) -c GR_Pulay_DM.c
Kerker_Mixing_Rhok.o: Kerker_Mixing_Rhok.c openmx_common.h
	$(CC) -c Kerker_Mixing_Rhok.c
DIIS_Mixing_Rhok.o: DIIS_Mixing_Rhok.c openmx_common.h
	$(CC) -c DIIS_Mixing_Rhok.c
Total_Energy.o: Total_Energy.c openmx_common.h
	$(CC) -c Total_Energy.c
Contract_Hamiltonian.o: Contract_Hamiltonian.c openmx_common.h
	$(CC) -c Contract_Hamiltonian.c
Contract_iHNL.o: Contract_iHNL.c openmx_common.h
	$(CC) -c Contract_iHNL.c
Cont_Matrix0.o: Cont_Matrix0.c openmx_common.h
	$(CC) -c Cont_Matrix0.c
Cont_Matrix1.o: Cont_Matrix1.c openmx_common.h
	$(CC) -c Cont_Matrix1.c
Cont_Matrix2.o: Cont_Matrix2.c openmx_common.h
	$(CC) -c Cont_Matrix2.c
Cont_Matrix3.o: Cont_Matrix3.c openmx_common.h
	$(CC) -c Cont_Matrix3.c
Opt_Contraction.o: Opt_Contraction.c openmx_common.h
	$(CC) -c Opt_Contraction.c
Initial_CntCoes.o: Initial_CntCoes.c openmx_common.h
	$(CC) -c Initial_CntCoes.c


Set_XC_Grid.o: Set_XC_Grid.c openmx_common.h
	$(CC) -c Set_XC_Grid.c
Get_Orbitals.o: Get_Orbitals.c openmx_common.h
	$(CC) -c Get_Orbitals.c
Get_dOrbitals.o: Get_dOrbitals.c openmx_common.h
	$(CC) -c Get_dOrbitals.c
Get_Cnt_Orbitals.o: Get_Cnt_Orbitals.c openmx_common.h
	$(CC) -c Get_Cnt_Orbitals.c
Get_Cnt_dOrbitals.o: Get_Cnt_dOrbitals.c openmx_common.h
	$(CC) -c Get_Cnt_dOrbitals.c
Gaunt.o: Gaunt.c openmx_common.h
	$(CC) -c Gaunt.c
RestartFileDFT.o: RestartFileDFT.c openmx_common.h
	$(CC) -c RestartFileDFT.c
Output_CompTime.o: Output_CompTime.c openmx_common.h
	$(CC) -c Output_CompTime.c
Merge_LogFile.o: Merge_LogFile.c openmx_common.h
	$(CC) -c Merge_LogFile.c
Make_FracCoord.o: Make_FracCoord.c openmx_common.h
	$(CC) -c Make_FracCoord.c
Make_InputFile_with_FinalCoord.o: Make_InputFile_with_FinalCoord.c openmx_common.h
	$(CC) -c Make_InputFile_with_FinalCoord.c
#
#
#
QuickSort.o: QuickSort.c openmx_common.h
	$(CC) -c QuickSort.c
Eigen_lapack.o: Eigen_lapack.c openmx_common.h lapack_prototypes.h
	$(CC) -c Eigen_lapack.c
EigenBand_lapack.o: EigenBand_lapack.c openmx_common.h lapack_prototypes.h
	$(CC) -c EigenBand_lapack.c
Eigen_PReHH.o: Eigen_PReHH.c openmx_common.h
	$(CC) -c Eigen_PReHH.c
Eigen_PHH.o: Eigen_PHH.c openmx_common.h
	$(CC) -c Eigen_PHH.c
BroadCast_ReMatrix.o: BroadCast_ReMatrix.c openmx_common.h
	$(CC) -c BroadCast_ReMatrix.c
BroadCast_ComplexMatrix.o: BroadCast_ComplexMatrix.c openmx_common.h
	$(CC) -c BroadCast_ComplexMatrix.c
lapack_dstedc1.o: lapack_dstedc1.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dstedc1.c
lapack_dstedc2.o: lapack_dstedc2.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dstedc2.c
lapack_dstegr1.o: lapack_dstegr1.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dstegr1.c
lapack_dstegr2.o: lapack_dstegr2.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dstegr2.c
lapack_dstevx1.o: lapack_dstevx1.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dstevx1.c
lapack_dstevx2.o: lapack_dstevx2.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dstevx2.c
lapack_dsteqr1.o: lapack_dsteqr1.c openmx_common.h lapack_prototypes.h
	$(CC) -c lapack_dsteqr1.c
Nonlocal_Basis.o: Nonlocal_Basis.c openmx_common.h
	$(CC) -c Nonlocal_Basis.c
Set_OLP_Kin.o: Set_OLP_Kin.c openmx_common.h 
	$(CC) -c Set_OLP_Kin.c
Set_Nonlocal.o: Set_Nonlocal.c openmx_common.h
	$(CC) -c Set_Nonlocal.c
Set_ProExpn_VNA.o: Set_ProExpn_VNA.c openmx_common.h
	$(CC) -c Set_ProExpn_VNA.c
Set_Hamiltonian.o: Set_Hamiltonian.c openmx_common.h
	$(CC) -c Set_Hamiltonian.c
Set_Vpot.o: Set_Vpot.c openmx_common.h
	$(CC) -c Set_Vpot.c
#
#
#
FT_PAO.o: FT_PAO.c openmx_common.h
	$(CC) -c FT_PAO.c
FT_NLP.o: FT_NLP.c openmx_common.h
	$(CC) -c FT_NLP.c
FT_ProExpn_VNA.o: FT_ProExpn_VNA.c openmx_common.h
	$(CC) -c FT_ProExpn_VNA.c
FT_VNA.o: FT_VNA.c openmx_common.h
	$(CC) -c FT_VNA.c
FT_ProductPAO.o: FT_ProductPAO.c openmx_common.h
	$(CC) -c FT_ProductPAO.c
#
#
#
Divide_Conquer.o: Divide_Conquer.c openmx_common.h
	$(CC) -c Divide_Conquer.c
GDivide_Conquer.o: GDivide_Conquer.c openmx_common.h
	$(CC) -c GDivide_Conquer.c
Divide_Conquer_Dosout.o: Divide_Conquer_Dosout.c openmx_common.h
	$(CC) -c Divide_Conquer_Dosout.c
GDivide_Conquer_Dosout.o: GDivide_Conquer_Dosout.c openmx_common.h
	$(CC) -c GDivide_Conquer_Dosout.c
Krylov.o: Krylov.c openmx_common.h
	$(CC) -c Krylov.c
Krylov_Dosout.o: Krylov_Dosout.c openmx_common.h
	$(CC) -c Krylov_Dosout.c
#
#
#
RecursionS_H.o: RecursionS_H.c openmx_common.h lapack_prototypes.h
	$(CC) -c RecursionS_H.c
IS_Lanczos.o: IS_Lanczos.c openmx_common.h
	$(CC) -c IS_Lanczos.c
IS_Taylor.o: IS_Taylor.c openmx_common.h
	$(CC) -c IS_Taylor.c
IS_Hotelling.o: IS_Hotelling.c openmx_common.h
	$(CC) -c IS_Hotelling.c
IS_LU.o: IS_LU.c openmx_common.h lapack_prototypes.h
	$(CC) -c IS_LU.c
#
#
#
MD_pac.o: MD_pac.c openmx_common.h lapack_prototypes.h
	$(CC) -c MD_pac.c
#
#
#
iterout.o: iterout.c openmx_common.h
	$(CC) -c iterout.c
iterout_md.o: iterout_md.c openmx_common.h
	$(CC) -c iterout_md.c
Allocate_Arrays.o: Allocate_Arrays.c openmx_common.h
	$(CC) -c Allocate_Arrays.c
Free_Arrays.o: Free_Arrays.c openmx_common.h
	$(CC) -c Free_Arrays.c
Init_List_YOUSO.o: Init_List_YOUSO.c openmx_common.h
	$(CC) -c Init_List_YOUSO.c
outputfile1.o: outputfile1.c openmx_common.h
	$(CC) -c outputfile1.c
malloc_multidimarray.o: malloc_multidimarray.c
	$(CC) -c malloc_multidimarray.c 
PrintMemory.o: PrintMemory.c
	$(CC) -c PrintMemory.c 
PrintMemory_Fix.o: PrintMemory_Fix.c openmx_common.h
	$(CC) -c PrintMemory_Fix.c 
dtime.o: dtime.c
	$(CC) -c dtime.c 
OutData.o: OutData.c openmx_common.h
	$(CC) -c OutData.c
init_alloc_first.o: init_alloc_first.c openmx_common.h
	$(CC) -c init_alloc_first.c
File_CntCoes.o: File_CntCoes.c openmx_common.h
	$(CC) -c File_CntCoes.c
SCF2File.o: SCF2File.c openmx_common.h
	$(CC) -c SCF2File.c
Cutoff.o: Cutoff.c openmx_common.h
	$(CC) -c Cutoff.c
setup_CPU_group.o: setup_CPU_group.c openmx_common.h
	$(CC) -c setup_CPU_group.c
Voronoi_Charge.o: Voronoi_Charge.c openmx_common.h
	$(CC) -c Voronoi_Charge.c
Voronoi_Orbital_Moment.o: Voronoi_Orbital_Moment.c openmx_common.h
	$(CC) -c Voronoi_Orbital_Moment.c
Fuzzy_Weight.o: Fuzzy_Weight.c openmx_common.h
	$(CC) -c Fuzzy_Weight.c
dampingF.o: dampingF.c openmx_common.h
	$(CC) -c dampingF.c
deri_dampingF.o: deri_dampingF.c openmx_common.h
	$(CC) -c deri_dampingF.c
Spherical_Bessel.o: Spherical_Bessel.c openmx_common.h
	$(CC) -c Spherical_Bessel.c
Generating_MP_Special_Kpt.o: Generating_MP_Special_Kpt.c openmx_common.h
	$(CC) -c Generating_MP_Special_Kpt.c

#
#
#
mimic_mpi.o: mimic_mpi.c mimic_mpi.h
	$(CC) -c mimic_mpi.c
mimic_omp.o: mimic_omp.c mimic_omp.h
	$(CC) -c mimic_omp.c
Make_Comm_Worlds.o: Make_Comm_Worlds.c
	$(CC) -c Make_Comm_Worlds.c
Set_Allocate_Atom2CPU.o: Set_Allocate_Atom2CPU.c openmx_common.h
	$(CC) -c Set_Allocate_Atom2CPU.c

#
#
# Maketest, Runtest, Memory_Leak_test, Force_test, Show_DFT_DATA
#
#

Maketest.o: Maketest.c openmx_common.h Inputtools.h
	$(CC) -c Maketest.c
Runtest.o: Runtest.c openmx_common.h Inputtools.h
	$(CC) -c Runtest.c
Memory_Leak_test.o: Memory_Leak_test.c openmx_common.h Inputtools.h
	$(CC) -c Memory_Leak_test.c
Force_test.o: Force_test.c openmx_common.h Inputtools.h
	$(CC) -c Force_test.c
Show_DFT_DATA.o: Show_DFT_DATA.c openmx_common.h Inputtools.h
	$(CC) -c Show_DFT_DATA.c


#
#
# install
#
#

install:	$(PROG)
	strip $(PROG)
	cp $(PROG) $(DESTDIR)/usr/bin/$(PROG)

#
#
# clean executable and object files 
#
#

clean:
	rm -f $(PROG) $(OBJS)

#
#
# programs for generating DOS from files *.Dos.val and *.Dos.vec
#
#

DosMain: DosMain.o Inputtools.o malloc_multidimarray.o Tetrahedron_Blochl.o 
	$(CC) -o $@ DosMain.o Inputtools.o malloc_multidimarray.o Tetrahedron_Blochl.o -lm 
	cp DosMain $(DESTDIR)/DosMain

DosMain.o :DosMain.c openmx_common.h
	$(CC) -o $@ -c DosMain.c
Tetrahedron_Blochl.o : Tetrahedron_Blochl.c
	$(CC) -o $@ -c Tetrahedron_Blochl.c 

#
#
#  exchange interaction coupling constant J between two atoms
#
#

jx: jx.o read_scfout.o Eigen_lapack.o 
	$(CC) jx.o read_scfout.o $(LIB) -lm -o jx
	cp jx $(DESTDIR)/jx

jx.o: jx.c read_scfout.h 
	$(CC) -c jx.c

#
#
# analysis_example
#
#

analysis_example: analysis_example.o read_scfout.o
	$(CC) analysis_example.o read_scfout.o $(LIB)  -lm -o analysis_example
	cp analysis_example $(DESTDIR)/analysis_example

analysis_example.o: analysis_example.c read_scfout.h 
	$(CC) -c analysis_example.c

read_scfout.o: read_scfout.c read_scfout.h 
	$(CC) -c read_scfout.c

#
#
# program for generating EPS from files *.out and *.vhart
#
#

OBJS_ESP  = esp.o Inputtools.o mimic_mpi.o 
esp:	$(OBJS_ESP)
	$(CC) $(OBJS_ESP) $(LIB) -lm -o $@
	cp esp $(DESTDIR)/esp
esp.o : esp.c Inputtools.h
	$(CC) -o $@ -c esp.c

#
#
# check_lead
#
#

check_lead: check_lead.o Inputtools.o
	$(CC) check_lead.o Inputtools.o -lm -o check_lead
	cp check_lead $(DESTDIR)/check_lead

check_lead.o: check_lead.c Inputtools.h 
	$(CC) -c check_lead.c

#
#
#  optical conductivity 
#
#

OpticalConductivityMain: OpticalConductivityMain.o \
              Inputtools.o  malloc_multidimarray.o
	$(CC) -o $@   OpticalConductivityMain.o  Inputtools.o  malloc_multidimarray.o -lm 
	cp OpticalConductivityMain $(DESTDIR)/OpticalConductivityMain

#
#
#  electric polarization using Berry's phase
#
#

OBJS_polB = polB.o read_scfout.o mimic_mpi.o
polB:	$(OBJS_polB)
	$(CC) $(OBJS_polB) $(LIB) -lm -o polB
	cp polB $(DESTDIR)/polB

polB.o: polB.c read_scfout.h 
	$(CC) -c polB.c

#
#
# test_mpi
#
#

test_mpi: test_mpi.o
	$(CC) test_mpi.o $(LIB) -lm -o test_mpi
	cp test_mpi $(DESTDIR)/test_mpi

test_mpi.o: test_mpi.c
	$(CC) -c test_mpi.c


#
#
# transport calculation using non-equilibrium Green's function method
#
#

TranMain_OBJ = \
   TranMain.o TRAN_Read.o Inputtools.o TRAN_Set_Value.o \
   TRAN_Calc_SurfGreen.o TRAN_Calc_SelfEnergy.o TRAN_Calc_CentGreen.o \
   TRAN_Calc_OneTransmission.o Lapack_LU_inverse.o \
   TRAN_Print.o TRAN_Distribute_Node.o mimic_mpi.o

TranMain: $(TranMain_OBJ)
	$(CC) $(TranMain_OBJ) $(LIB) -lm -o $@ 




MAIN_TRAN_Display_Gridvalue: MAIN_TRAN_Display_Gridvalue.o TRAN_Read.o TRAN_Print.o
	$(CC) -o $@  MAIN_TRAN_Display_Gridvalue.o TRAN_Read.o TRAN_Print.o -lm $(LIB)  


TRAN_Allocate.o: tran_variables.h 
TRAN_Calc_GridBound.o: tran_variables.h 
TRAN_Calc_Transmission.o: tran_variables.h 
TRAN_DFT.o: tran_variables.h 
TRAN_DFT_Dosout.o: tran_variables.h
TRAN_Input_std.o: tran_variables.h 
TRAN_Input_std_Atoms.o: tran_variables.h 
TRAN_Output_HKS.o: tran_variables.h 
TRAN_Output_Trans_HS.o: tran_variables.h 
TRAN_Output_Transmission.o: tran_variables.h 
TRAN_Overwrite_Densitygrid.o: tran_variables.h 
TRAN_Poisson.o: tran_variables.h 
TRAN_RestartFile.o: tran_variables.h 
TRAN_Set_CentOverlap.o: tran_variables.h 
TRAN_Set_Electrode_Grid.o: tran_variables.h 
TRAN_Set_IntegPath.o: tran_variables.h lapack_prototypes.h  
TRAN_Set_SurfOverlap.o: tran_variables.h 
TRAN_adjust_Grid_Origin.o: tran_variables.h 
TRAN_adjust_Ngrid.o: tran_variables.h 
openmx.o: tran_variables.h 
DFT.o: tran_prototypes.h 
Input_std.o: tran_prototypes.h 
Lapack_LU_inverse.o: tran_prototypes.h 
MAIN_TRAN_Calc_Transmission.o: tran_prototypes.h 
TranMain.o: tran_prototypes.h 
Poisson.o: tran_prototypes.h 
TRAN_Allocate.o: tran_prototypes.h 
TRAN_Apply_Bias2e.o: tran_prototypes.h 
TRAN_Calc_CentGreen.o: tran_prototypes.h 
TRAN_Calc_CentGreenLesser.o: tran_prototypes.h 
TRAN_Calc_GC_LorR.o: tran_prototypes.h 
TRAN_Calc_OneTransmission.o: tran_prototypes.h 
TRAN_Calc_SelfEnergy.o: tran_prototypes.h 
TRAN_Calc_SurfGreen.o: tran_prototypes.h 
TRAN_Calc_Transmission.o: tran_prototypes.h 
TRAN_Calc_VHartree_G0.o: tran_prototypes.h 
TRAN_Calc_VHartree_Gnon0.o: tran_prototypes.h 
TRAN_Connect_Read_Density.o: tran_prototypes.h 
TRAN_Credit.o: tran_prototypes.h 
TRAN_DFT.o: tran_prototypes.h
TRAN_DFT_Dosout.o: tran_prototypes.h
TRAN_FFT_interpolation3c.o: tran_prototypes.h 
TRAN_Input_std.o: tran_prototypes.h 
TRAN_Input_std_Atoms.o: tran_prototypes.h 
TRAN_Interp_ElectrodeDensity_Grid.o: tran_prototypes.h 
TRAN_Output_HKS.o: tran_prototypes.h  
TRAN_Output_HKS_Write_Grid.o: tran_prototypes.h 
TRAN_Output_Trans_HS.o: tran_prototypes.h 
TRAN_Output_Transmission.o: tran_prototypes.h 
TRAN_Overwrite_Densitygrid.o: tran_prototypes.h 
TRAN_Overwrite_V2.o: tran_prototypes.h 
TRAN_Poisson.o: tran_prototypes.h 
TRAN_Print.o: tran_prototypes.h 
TRAN_Print_Grid.o: tran_prototypes.h 
TRAN_Read.o: tran_prototypes.h 
TRAN_RestartFile.o: tran_prototypes.h  
TRAN_Set_CentOverlap.o: tran_prototypes.h 
TRAN_Set_Electrode_Grid.o: tran_prototypes.h  
TRAN_Set_IntegPath.o: tran_prototypes.h lapack_prototypes.h 
TRAN_Set_PathEnergyStr.o: tran_prototypes.h 
TRAN_Set_SurfOverlap.o: tran_prototypes.h 
TRAN_Set_Value.o: tran_prototypes.h 
TRAN_Check_Region_Lead.o: tran_prototypes.h 
TRAN_Check_Region.o: tran_prototypes.h 
openmx.o: tran_prototypes.h 
truncation.o: tran_prototypes.h 

