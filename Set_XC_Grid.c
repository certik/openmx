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

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif




void Set_XC_Grid(int XC_P_switch, int XC_switch)
{
  /****************************************************
        XC_P_switch:
            0  \epsilon_XC (XC energy density)  
            1  \mu_XC      (XC potential)  
            2  \epsilon_XC - \mu_XC
  ****************************************************/

  int MN,MN1,MN2,i,j,k,ri,ri1,ri2;
  int i1,i2,j1,j2,k1,k2,n,nmax;
  double den_min=1.0e-14; 
  double Ec_unif[1],Vc_unif[2],Exc[2],Vxc[2];
  double Ex_unif[1],Vx_unif[2],tot_den;
  double ED[2],GDENS[3][2];
  double DEXDD[2],DECDD[2];
  double DEXDGD[3][2],DECDGD[3][2];
  double ***dEXC_dGD,***dDen_Grid;
  double up_x_a,up_x_b,up_x_c;
  double up_y_a,up_y_b,up_y_c;
  double up_z_a,up_z_b,up_z_c;
  double dn_x_a,dn_x_b,dn_x_c;
  double dn_y_a,dn_y_b,dn_y_c;
  double dn_z_a,dn_z_b,dn_z_c;
  double up_a,up_b,up_c;
  double dn_a,dn_b,dn_c;

  double tmp0,tmp1;
  double cot,sit,sip,cop,phi,theta;
  double detA,igtv[4][4];
  int numprocs,myid;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /****************************************************
   when GGA, allocation

   double dEXC_dGD[2][3][My_NumGrid1]
   double dDen_Grid[2][3][My_NumGrid1]
  ****************************************************/

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (XC_switch==4){

    dDen_Grid = (double***)malloc(sizeof(double**)*2); 
    for (k=0; k<=1; k++){
      dDen_Grid[k] = (double**)malloc(sizeof(double*)*3); 
      for (i=0; i<3; i++){
        dDen_Grid[k][i] = (double*)malloc(sizeof(double)*My_NumGrid1); 
        for (j=0; j<My_NumGrid1; j++) dDen_Grid[k][i][j] = 0.0;
      }
    }

    if (XC_P_switch!=0){
      dEXC_dGD = (double***)malloc(sizeof(double**)*2); 
      for (k=0; k<=1; k++){
        dEXC_dGD[k] = (double**)malloc(sizeof(double*)*3); 
        for (i=0; i<3; i++){
          dEXC_dGD[k][i] = (double*)malloc(sizeof(double)*My_NumGrid1); 
          for (j=0; j<My_NumGrid1; j++) dEXC_dGD[k][i][j] = 0.0;
        }
      }
    }

    /* PrintMemory */
    PrintMemory("Set_XC_Grid: dDen_Grid", sizeof(double)*6*My_NumGrid1, NULL);
    PrintMemory("Set_XC_Grid: dEXC_dGD",  sizeof(double)*6*My_NumGrid1, NULL);

    /****************************************************
     calculate dDen_Grid
    ****************************************************/
 
    detA =   gtv[1][1]*gtv[2][2]*gtv[3][3]
           + gtv[1][2]*gtv[2][3]*gtv[3][1]
           + gtv[1][3]*gtv[2][1]*gtv[3][2]
           - gtv[1][3]*gtv[2][2]*gtv[3][1]
           - gtv[1][2]*gtv[2][1]*gtv[3][3]
           - gtv[1][1]*gtv[2][3]*gtv[3][2];     

    igtv[1][1] =  (gtv[2][2]*gtv[3][3] - gtv[2][3]*gtv[3][2])/detA;
    igtv[2][1] = -(gtv[2][1]*gtv[3][3] - gtv[2][3]*gtv[3][1])/detA;
    igtv[3][1] =  (gtv[2][1]*gtv[3][2] - gtv[2][2]*gtv[3][1])/detA; 

    igtv[1][2] = -(gtv[1][2]*gtv[3][3] - gtv[1][3]*gtv[3][2])/detA;
    igtv[2][2] =  (gtv[1][1]*gtv[3][3] - gtv[1][3]*gtv[3][1])/detA;
    igtv[3][2] = -(gtv[1][1]*gtv[3][2] - gtv[1][2]*gtv[3][1])/detA; 

    igtv[1][3] =  (gtv[1][2]*gtv[2][3] - gtv[1][3]*gtv[2][2])/detA;
    igtv[2][3] = -(gtv[1][1]*gtv[2][3] - gtv[1][3]*gtv[2][1])/detA;
    igtv[3][3] =  (gtv[1][1]*gtv[2][2] - gtv[1][2]*gtv[2][1])/detA; 

#pragma omp parallel shared(igtv,dDen_Grid,PCCDensity_Grid,PCC_switch,Density_Grid,den_min,My_Cell0,My_Cell1,Ngrid3,Ngrid2,Num_Cells0) private(OMPID,Nthrds,Nprocs,nmax,n,i,j,k,ri,ri1,ri2,i1,i2,j1,j2,k1,k2,MN,MN1,MN2,up_a,dn_a,up_b,dn_b,up_c,dn_c)
    {

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();
      nmax = Num_Cells0*Ngrid2*Ngrid3; 

      for (n=OMPID*nmax/Nthrds; n<(OMPID+1)*nmax/Nthrds; n++){

	i = n/(Ngrid2*Ngrid3);
	j = (n-i*Ngrid2*Ngrid3)/Ngrid3;
	k = n - i*Ngrid2*Ngrid3 - j*Ngrid3; 
	ri = My_Cell1[i];

	/* find ri1, ri2, i1, and i2 */

	if (ri==0){
	  ri1 = Ngrid1 - 1;
	  ri2 = 1;        
	  i1 = My_Cell0[ri1];
	  i2 = My_Cell0[ri2];
	}
	else if (ri==(Ngrid1-1)){
	  ri1 = Ngrid1 - 2;
	  ri2 = 0;
	  i1 = My_Cell0[ri1];
	  i2 = My_Cell0[ri2];
	}      
	else{
	  ri1 = ri - 1;
	  ri2 = ri + 1;
	  i1 = My_Cell0[ri1];
	  i2 = My_Cell0[ri2];
	}

	/* because we have +-1 buffer cells. */

	if (i1!=-1 && i2!=-1){

	  /* find j1 and j2 */

	  if (j==0){
	    j1 = Ngrid2 - 1;
	    j2 = 1;
	  }
	  else if (j==(Ngrid2-1)){
	    j1 = Ngrid2 - 2;
	    j2 = 0;
	  }
	  else{
	    j1 = j - 1;
	    j2 = j + 1;
	  }

	  /* find k1 and k2 */

	  if (k==0){
	    k1 = Ngrid3 - 1;
	    k2 = 1;
	  }
	  else if (k==(Ngrid3-1)){
	    k1 = Ngrid3 - 2;
	    k2 = 0;
	  }
	  else{
	    k1 = k - 1;
	    k2 = k + 1;
	  }  

	  /* set MN */

	  MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 

	  /* set dDen_Grid */

	  if ( den_min<(Density_Grid[0][MN]+Density_Grid[1][MN]) ){

	    /* a-axis */

	    MN1 = i1*Ngrid2*Ngrid3 + j*Ngrid3 + k;
	    MN2 = i2*Ngrid2*Ngrid3 + j*Ngrid3 + k;

	    if (PCC_switch==0) {
	      up_a = Density_Grid[0][MN2] - Density_Grid[0][MN1];
	      dn_a = Density_Grid[1][MN2] - Density_Grid[1][MN1];
	    }
	    else if (PCC_switch==1) {
	      up_a = Density_Grid[0][MN2] + PCCDensity_Grid[MN2]
	           - Density_Grid[0][MN1] - PCCDensity_Grid[MN1];
	      dn_a = Density_Grid[1][MN2] + PCCDensity_Grid[MN2]
	           - Density_Grid[1][MN1] - PCCDensity_Grid[MN1];
	    }

	    /* b-axis */

	    MN1 = i*Ngrid2*Ngrid3 + j1*Ngrid3 + k; 
	    MN2 = i*Ngrid2*Ngrid3 + j2*Ngrid3 + k; 

	    if (PCC_switch==0) {
	      up_b = Density_Grid[0][MN2] - Density_Grid[0][MN1];
	      dn_b = Density_Grid[1][MN2] - Density_Grid[1][MN1];
	    }
	    else if (PCC_switch==1) {
	      up_b = Density_Grid[0][MN2] + PCCDensity_Grid[MN2]
	           - Density_Grid[0][MN1] - PCCDensity_Grid[MN1];
	      dn_b = Density_Grid[1][MN2] + PCCDensity_Grid[MN2]
	           - Density_Grid[1][MN1] - PCCDensity_Grid[MN1];
	    }

	    /* c-axis */

	    MN1 = i*Ngrid2*Ngrid3 + j*Ngrid3 + k1; 
	    MN2 = i*Ngrid2*Ngrid3 + j*Ngrid3 + k2; 

	    if (PCC_switch==0) {
	      up_c = Density_Grid[0][MN2] - Density_Grid[0][MN1];
	      dn_c = Density_Grid[1][MN2] - Density_Grid[1][MN1];
	    }
	    else if (PCC_switch==1) {
	      up_c = Density_Grid[0][MN2] + PCCDensity_Grid[MN2]
	           - Density_Grid[0][MN1] - PCCDensity_Grid[MN1];
	      dn_c = Density_Grid[1][MN2] + PCCDensity_Grid[MN2]
	           - Density_Grid[1][MN1] - PCCDensity_Grid[MN1];
	    }

	    /* up */
	    dDen_Grid[0][0][MN] = 0.5*(igtv[1][1]*up_a + igtv[1][2]*up_b + igtv[1][3]*up_c);
	    dDen_Grid[0][1][MN] = 0.5*(igtv[2][1]*up_a + igtv[2][2]*up_b + igtv[2][3]*up_c);
	    dDen_Grid[0][2][MN] = 0.5*(igtv[3][1]*up_a + igtv[3][2]*up_b + igtv[3][3]*up_c);

	    /* down */
	    dDen_Grid[1][0][MN] = 0.5*(igtv[1][1]*dn_a + igtv[1][2]*dn_b + igtv[1][3]*dn_c);
	    dDen_Grid[1][1][MN] = 0.5*(igtv[2][1]*dn_a + igtv[2][2]*dn_b + igtv[2][3]*dn_c);
	    dDen_Grid[1][2][MN] = 0.5*(igtv[3][1]*dn_a + igtv[3][2]*dn_b + igtv[3][3]*dn_c);
	  }

	  else{
	    dDen_Grid[0][0][MN] = 0.0;
	    dDen_Grid[0][1][MN] = 0.0;
	    dDen_Grid[0][2][MN] = 0.0;
	    dDen_Grid[1][0][MN] = 0.0;
	    dDen_Grid[1][1][MN] = 0.0;
	    dDen_Grid[1][2][MN] = 0.0;
	  }

	} /* if (i1!=-1 && i2!=-1) */
      } /* n */

#pragma omp flush(dDen_Grid)

    } /* #pragma omp parallel */
  } /* if (XC_switch==4) */ 

  /****************************************************
   loop MN
  ****************************************************/

#pragma omp parallel shared(dDen_Grid,dEXC_dGD,den_min,Vxc_Grid,My_NumGrid1,XC_P_switch,XC_switch,Density_Grid,PCC_switch,PCCDensity_Grid) private(OMPID,Nthrds,Nprocs,MN,tot_den,tmp0,ED,Exc,Ec_unif,Vc_unif,Vxc,Ex_unif,Vx_unif,GDENS,DEXDD,DECDD,DEXDGD,DECDGD)
  {

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (MN=OMPID*My_NumGrid1/Nthrds; MN<(OMPID+1)*My_NumGrid1/Nthrds; MN++){

      switch(XC_switch){
        
	/******************************************************************
         LDA (Ceperly-Alder)

         constructed by Ceperly and Alder,
         ref.
         D. M. Ceperley, Phys. Rev. B18, 3126 (1978)
         D. M. Ceperley and B. J. Alder, Phys. Rev. Lett., 45, 566 (1980) 

         and parametrized by Perdew and Zunger.
         ref.
         J. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)
	******************************************************************/
        
      case 1:
        
	tot_den = Density_Grid[0][MN] + Density_Grid[1][MN];

	/* partial core correction */
	if (PCC_switch==1) {
	  tot_den += PCCDensity_Grid[MN]*2.0;
	}

	tmp0 = XC_Ceperly_Alder(tot_den,XC_P_switch);
	Vxc_Grid[0][MN] = tmp0;
	Vxc_Grid[1][MN] = tmp0;
        
	break;

	/******************************************************************
         LSDA-CA (Ceperly-Alder)

         constructed by Ceperly and Alder,
         ref.
         D. M. Ceperley, Phys. Rev. B18, 3126 (1978)
         D. M. Ceperley and B. J. Alder, Phys. Rev. Lett., 45, 566 (1980) 

         and parametrized by Perdew and Zunger.
         ref.
         J. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)
	******************************************************************/

      case 2:

	ED[0] = Density_Grid[0][MN];
	ED[1] = Density_Grid[1][MN];

	/* partial core correction */
	if (PCC_switch==1) {
	  ED[0] += PCCDensity_Grid[MN];
	  ED[1] += PCCDensity_Grid[MN];
	}

	XC_CA_LSDA(ED[0], ED[1], Exc, XC_P_switch);
	Vxc_Grid[0][MN] = Exc[0];
	Vxc_Grid[1][MN] = Exc[1];

	break;

	/******************************************************************
         LSDA-PW (PW91)
         used as Grad\rho = 0 in their GGA formalism

         ref.
         J.P.Perdew and Yue Wang, Phys. Rev. B45, 13244 (1992) 
	******************************************************************/

      case 3:

	ED[0] = Density_Grid[0][MN];
	ED[1] = Density_Grid[1][MN];

	/* partial core correction */
	if (PCC_switch==1) {
	  ED[0] += PCCDensity_Grid[MN];
	  ED[1] += PCCDensity_Grid[MN];
	}

	if ((ED[0]+ED[1])<den_min){
	  Vxc_Grid[0][MN] = 0.0;
	  Vxc_Grid[1][MN] = 0.0;
	}
	else{
      
	  if (XC_P_switch==0){

	    XC_PW91C(ED,Ec_unif,Vc_unif);

	    Vxc[0] = Vc_unif[0];
	    Vxc[1] = Vc_unif[1];
	    Exc[0] = Ec_unif[0];

	    XC_EX(1,2.0*ED[0],ED,Ex_unif,Vx_unif);
	    Vxc[0] = Vxc[0] + Vx_unif[0];
	    Exc[1] = 2.0*ED[0]*Ex_unif[0];

	    XC_EX(1,2.0*ED[1],ED,Ex_unif,Vx_unif);
	    Vxc[1] += Vx_unif[0];
	    Exc[1] += 2.0*ED[1]*Ex_unif[0];

	    Exc[1] = 0.5*Exc[1]/(ED[0]+ED[1]);

	    Vxc_Grid[0][MN] = Exc[0] + Exc[1];
	    Vxc_Grid[1][MN] = Exc[0] + Exc[1];
	  }

	  else if (XC_P_switch==1){
	    XC_PW91C(ED,Ec_unif,Vc_unif);
	    Vxc_Grid[0][MN] = Vc_unif[0];
	    Vxc_Grid[1][MN] = Vc_unif[1];

	    XC_EX(1,2.0*ED[0],ED,Ex_unif,Vx_unif);
	    Vxc_Grid[0][MN] = Vxc_Grid[0][MN] + Vx_unif[0];

	    XC_EX(1,2.0*ED[1],ED,Ex_unif,Vx_unif);
	    Vxc_Grid[1][MN] = Vxc_Grid[1][MN] + Vx_unif[0];
	  }

	  else if (XC_P_switch==2){

	    XC_PW91C(ED,Ec_unif,Vc_unif);

	    Vxc[0] = Vc_unif[0];
	    Vxc[1] = Vc_unif[1];
	    Exc[0] = Ec_unif[0];

	    XC_EX(1,2.0*ED[0],ED,Ex_unif,Vx_unif);
	    Vxc[0]  = Vxc[0] + Vx_unif[0];
	    Exc[1]  = 2.0*ED[0]*Ex_unif[0];

	    XC_EX(1,2.0*ED[1],ED,Ex_unif,Vx_unif);
	    Vxc[1] += Vx_unif[0];
	    Exc[1] += 2.0*ED[1]*Ex_unif[0];

	    Exc[1] = 0.5*Exc[1]/(ED[0]+ED[1]);

	    Vxc_Grid[0][MN] = Exc[0] + Exc[1] - Vxc[0];
	    Vxc_Grid[1][MN] = Exc[0] + Exc[1] - Vxc[1];
	  }
	}

	break;

	/******************************************************************
         GGA-PBE
         ref.
         J. P. Perdew, K. Burke, and M. Ernzerhof,
         Phys. Rev. Lett. 77, 3865 (1996).
	******************************************************************/

      case 4:

	/****************************************************
         ED[0]       density of up spin:     n_up   
         ED[1]       density of down spin:   n_down

         GDENS[0][0] derivative (x) of density of up spin
         GDENS[1][0] derivative (y) of density of up spin
         GDENS[2][0] derivative (z) of density of up spin
         GDENS[0][1] derivative (x) of density of down spin
         GDENS[1][1] derivative (y) of density of down spin
         GDENS[2][1] derivative (z) of density of down spin

         DEXDD[0]    d(fx)/d(n_up) 
         DEXDD[1]    d(fx)/d(n_down) 
         DECDD[0]    d(fc)/d(n_up) 
         DECDD[1]    d(fc)/d(n_down) 

         n'_up_x   = d(n_up)/d(x)
         n'_up_y   = d(n_up)/d(y)
         n'_up_z   = d(n_up)/d(z)
         n'_down_x = d(n_down)/d(x)
         n'_down_y = d(n_down)/d(y)
         n'_down_z = d(n_down)/d(z)
       
         DEXDGD[0][0] d(fx)/d(n'_up_x) 
         DEXDGD[1][0] d(fx)/d(n'_up_y) 
         DEXDGD[2][0] d(fx)/d(n'_up_z) 
         DEXDGD[0][1] d(fx)/d(n'_down_x) 
         DEXDGD[1][1] d(fx)/d(n'_down_y) 
         DEXDGD[2][1] d(fx)/d(n'_down_z) 

         DECDGD[0][0] d(fc)/d(n'_up_x) 
         DECDGD[1][0] d(fc)/d(n'_up_y) 
         DECDGD[2][0] d(fc)/d(n'_up_z) 
         DECDGD[0][1] d(fc)/d(n'_down_x) 
         DECDGD[1][1] d(fc)/d(n'_down_y) 
         DECDGD[2][1] d(fc)/d(n'_down_z) 
	****************************************************/

	ED[0] = Density_Grid[0][MN];
	ED[1] = Density_Grid[1][MN];

	if ((ED[0]+ED[1])<den_min){
	  Vxc_Grid[0][MN] = 0.0;
	  Vxc_Grid[1][MN] = 0.0;

	  /* later add its derivatives */
	  if (XC_P_switch!=0){
	    dEXC_dGD[0][0][MN] = 0.0;
	    dEXC_dGD[0][1][MN] = 0.0;
	    dEXC_dGD[0][2][MN] = 0.0;

	    dEXC_dGD[1][0][MN] = 0.0;
	    dEXC_dGD[1][1][MN] = 0.0;
	    dEXC_dGD[1][2][MN] = 0.0;
	  }
	}
     
	else{

	  GDENS[0][0] = dDen_Grid[0][0][MN];
	  GDENS[1][0] = dDen_Grid[0][1][MN];
	  GDENS[2][0] = dDen_Grid[0][2][MN];
	  GDENS[0][1] = dDen_Grid[1][0][MN];
	  GDENS[1][1] = dDen_Grid[1][1][MN];
	  GDENS[2][1] = dDen_Grid[1][2][MN];

	  if (PCC_switch==1) {
	    ED[0] += PCCDensity_Grid[MN];
	    ED[1] += PCCDensity_Grid[MN];
	  }

	  XC_PBE(ED, GDENS, Exc, DEXDD, DECDD, DEXDGD, DECDGD);

	  /* XC energy density */
	  if      (XC_P_switch==0){
	    Vxc_Grid[0][MN] = Exc[0] + Exc[1];
	    Vxc_Grid[1][MN] = Exc[0] + Exc[1];
	  }

	  /* XC potential */
	  else if (XC_P_switch==1){
	    Vxc_Grid[0][MN] = DEXDD[0] + DECDD[0];
	    Vxc_Grid[1][MN] = DEXDD[1] + DECDD[1];
	  }

	  /* XC energy density - XC potential */
	  else if (XC_P_switch==2){
	    Vxc_Grid[0][MN] = Exc[0] + Exc[1] - DEXDD[0] - DECDD[0];
	    Vxc_Grid[1][MN] = Exc[0] + Exc[1] - DEXDD[1] - DECDD[1];
	  }

	  /* later add its derivatives */
	  if (XC_P_switch!=0){
	    dEXC_dGD[0][0][MN] = DEXDGD[0][0] + DECDGD[0][0];
	    dEXC_dGD[0][1][MN] = DEXDGD[1][0] + DECDGD[1][0];
	    dEXC_dGD[0][2][MN] = DEXDGD[2][0] + DECDGD[2][0];

	    dEXC_dGD[1][0][MN] = DEXDGD[0][1] + DECDGD[0][1];
	    dEXC_dGD[1][1][MN] = DEXDGD[1][1] + DECDGD[1][1];
	    dEXC_dGD[1][2][MN] = DEXDGD[2][1] + DECDGD[2][1];
	  }
	}

	break;

      } /* switch(XC_switch) */
    }   /* MN */

#pragma omp flush(dEXC_dGD)

  } /* #pragma omp parallel */

  /****************************************************
        calculate the second part of XC potential
               when GGA and XC_P_switch!=0
  ****************************************************/

  if (XC_switch==4 && XC_P_switch!=0){

#pragma omp parallel shared(XC_P_switch,Vxc_Grid,igtv,dEXC_dGD,Density_Grid,den_min,My_Cell0,My_Cell1,Num_Cells0,Ngrid2,Ngrid3) private(OMPID,Nthrds,Nprocs,nmax,n,i,j,k,ri,ri1,ri2,i1,i2,j1,j2,k1,k2,MN,MN1,MN2,up_x_a,up_y_a,up_z_a,dn_x_a,dn_y_a,dn_z_a,up_x_b,up_y_b,up_z_b,dn_x_b,dn_y_b,dn_z_b,up_x_c,up_y_c,up_z_c,dn_x_c,dn_y_c,dn_z_c,tmp0,tmp1)
    {

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();
      nmax = Num_Cells0*Ngrid2*Ngrid3; 

      for (n=OMPID*nmax/Nthrds; n<(OMPID+1)*nmax/Nthrds; n++){

	i = n/(Ngrid2*Ngrid3);
	j = (n-i*Ngrid2*Ngrid3)/Ngrid3;
	k = n - i*Ngrid2*Ngrid3 - j*Ngrid3; 
	ri = My_Cell1[i];

	/* find ri1, ri2, i1, and i2 */

	if (ri==0){
	  ri1 = Ngrid1 - 1;
	  ri2 = 1;        
	  i1 = My_Cell0[ri1];
	  i2 = My_Cell0[ri2];
	}
	else if (ri==(Ngrid1-1)){
	  ri1 = Ngrid1 - 2;
	  ri2 = 0;
	  i1 = My_Cell0[ri1];
	  i2 = My_Cell0[ri2];
	}      
	else{
	  ri1 = ri - 1;
	  ri2 = ri + 1;
	  i1 = My_Cell0[ri1];
	  i2 = My_Cell0[ri2];
	}

	if (i1!=-1 && i2!=-1){

	  /* find j1 and j2 */

	  if (j==0){
	    j1 = Ngrid2 - 1;
	    j2 = 1;
	  }
	  else if (j==(Ngrid2-1)){
	    j1 = Ngrid2 - 2;
	    j2 = 0;
	  }
	  else{
	    j1 = j - 1;
	    j2 = j + 1;
	  }

	  /* find k1 and k2 */

	  if (k==0){
	    k1 = Ngrid3 - 1;
	    k2 = 1;
	  }
	  else if (k==(Ngrid3-1)){
	    k1 = Ngrid3 - 2;
	    k2 = 0;
	  }
	  else{
	    k1 = k - 1;
	    k2 = k + 1;
	  }  

	  /* set MN */

	  MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 

	  /* set Vxc_Grid */

	  if ( den_min<(Density_Grid[0][MN]+Density_Grid[1][MN]) ){

	    /* a-axis */

	    MN1 = i1*Ngrid2*Ngrid3 + j*Ngrid3 + k;
	    MN2 = i2*Ngrid2*Ngrid3 + j*Ngrid3 + k;

	    up_x_a = dEXC_dGD[0][0][MN2] - dEXC_dGD[0][0][MN1];
	    up_y_a = dEXC_dGD[0][1][MN2] - dEXC_dGD[0][1][MN1];
	    up_z_a = dEXC_dGD[0][2][MN2] - dEXC_dGD[0][2][MN1];

	    dn_x_a = dEXC_dGD[1][0][MN2] - dEXC_dGD[1][0][MN1];
	    dn_y_a = dEXC_dGD[1][1][MN2] - dEXC_dGD[1][1][MN1];
	    dn_z_a = dEXC_dGD[1][2][MN2] - dEXC_dGD[1][2][MN1];

	    /* b-axis */

	    MN1 = i*Ngrid2*Ngrid3 + j1*Ngrid3 + k; 
	    MN2 = i*Ngrid2*Ngrid3 + j2*Ngrid3 + k; 

	    up_x_b = dEXC_dGD[0][0][MN2] - dEXC_dGD[0][0][MN1];
	    up_y_b = dEXC_dGD[0][1][MN2] - dEXC_dGD[0][1][MN1];
	    up_z_b = dEXC_dGD[0][2][MN2] - dEXC_dGD[0][2][MN1];

	    dn_x_b = dEXC_dGD[1][0][MN2] - dEXC_dGD[1][0][MN1];
	    dn_y_b = dEXC_dGD[1][1][MN2] - dEXC_dGD[1][1][MN1];
	    dn_z_b = dEXC_dGD[1][2][MN2] - dEXC_dGD[1][2][MN1];

	    /* c-axis */

	    MN1 = i*Ngrid2*Ngrid3 + j*Ngrid3 + k1; 
	    MN2 = i*Ngrid2*Ngrid3 + j*Ngrid3 + k2; 

	    up_x_c = dEXC_dGD[0][0][MN2] - dEXC_dGD[0][0][MN1];
	    up_y_c = dEXC_dGD[0][1][MN2] - dEXC_dGD[0][1][MN1];
	    up_z_c = dEXC_dGD[0][2][MN2] - dEXC_dGD[0][2][MN1];

	    dn_x_c = dEXC_dGD[1][0][MN2] - dEXC_dGD[1][0][MN1];
	    dn_y_c = dEXC_dGD[1][1][MN2] - dEXC_dGD[1][1][MN1];
	    dn_z_c = dEXC_dGD[1][2][MN2] - dEXC_dGD[1][2][MN1];

	    /* up */

	    tmp0 = igtv[1][1]*up_x_a + igtv[1][2]*up_x_b + igtv[1][3]*up_x_c
	      + igtv[2][1]*up_y_a + igtv[2][2]*up_y_b + igtv[2][3]*up_y_c
	      + igtv[3][1]*up_z_a + igtv[3][2]*up_z_b + igtv[3][3]*up_z_c;
	    tmp0 = 0.5*tmp0;

	    /* down */

	    tmp1 = igtv[1][1]*dn_x_a + igtv[1][2]*dn_x_b + igtv[1][3]*dn_x_c
	      + igtv[2][1]*dn_y_a + igtv[2][2]*dn_y_b + igtv[2][3]*dn_y_c
	      + igtv[3][1]*dn_z_a + igtv[3][2]*dn_z_b + igtv[3][3]*dn_z_c;
	    tmp1 = 0.5*tmp1;

	    /* XC potential */

	    if (XC_P_switch==1){
	      Vxc_Grid[0][MN] -= tmp0; 
	      Vxc_Grid[1][MN] -= tmp1;
	    }

	    /* XC energy density - XC potential */

	    else if (XC_P_switch==2){
	      Vxc_Grid[0][MN] += tmp0; 
	      Vxc_Grid[1][MN] += tmp1;
	    }

	  }
	}
      }

#pragma omp flush(Vxc_Grid)

    } /* #pragma omp parallel */
  } /* if (XC_switch==4 && XC_P_switch!=0) */

  /****************************************************
            In case of non-collinear spin DFT 
  ****************************************************/

  if (SpinP_switch==3 && XC_P_switch!=2){

#pragma omp parallel shared(Density_Grid,Vxc_Grid,My_NumGrid1) private(OMPID,Nthrds,Nprocs,MN,tmp0,tmp1,theta,phi,sit,cot,sip,cop)
    {

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (MN=OMPID*My_NumGrid1/Nthrds; MN<(OMPID+1)*My_NumGrid1/Nthrds; MN++){

	tmp0 = 0.5*(Vxc_Grid[0][MN] + Vxc_Grid[1][MN]);
	tmp1 = 0.5*(Vxc_Grid[0][MN] - Vxc_Grid[1][MN]);
	theta = Density_Grid[2][MN];
	phi   = Density_Grid[3][MN];
	sit = sin(theta);
	cot = cos(theta);
	sip = sin(phi);
	cop = cos(phi);

	Vxc_Grid[0][MN] =  tmp0 + cot*tmp1;  /* Re Vxc11 */
	Vxc_Grid[1][MN] =  tmp0 - cot*tmp1;  /* Re Vxc22 */
	Vxc_Grid[2][MN] =  tmp1*sit*cop;     /* Re Vxc12 */
	Vxc_Grid[3][MN] = -tmp1*sit*sip;     /* Im Vxc12 */ 
      }

#pragma omp flush(Vxc_Grid)

    } /* #pragma omp parallel */ 
  }

  /*
  {
    int hN1,hN2,hN3,i;
    double Re11,Re22,Re12,Im12;

    hN1 = Ngrid1/2;
    hN2 = Ngrid2/2;
    hN3 = Ngrid3/2;

    for (i=0; i<Num_Cells0; i++){

    MN = i*Ngrid2*Ngrid3 + hN2*Ngrid3 + hN3;
 

    Re11 = Vxc_Grid[0][MN];
    Re22 = Vxc_Grid[1][MN];
    Re12 = Vxc_Grid[2][MN];
    Im12 = Vxc_Grid[3][MN];

    printf("MN=%4d %15.12f %15.12f %15.12f %15.12f\n",
           MN,Re11,Re22,Re12,Im12);
    }
  }


  MPI_Finalize();
  exit(0);
  */


  /****************************************************
   In case of GGA,
   free arrays
   double dEXC_dGD[2][3][My_NumGrid1]
   double dDen_Grid[2][3][My_NumGrid1]
  ****************************************************/

  if (XC_switch==4){

    for (k=0; k<=1; k++){
      for (i=0; i<3; i++){
        free(dDen_Grid[k][i]);
      }
      free(dDen_Grid[k]);
    }
    free(dDen_Grid);

    if (XC_P_switch!=0){
      for (k=0; k<=1; k++){
        for (i=0; i<3; i++){
          free(dEXC_dGD[k][i]);
        }
        free(dEXC_dGD[k]);
      }
      free(dEXC_dGD);
    }
  }
}
