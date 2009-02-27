#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*typedef struct {double r,i; } dcomplex; */


#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"


#define print_stdout 0



/********************************************************************************
 fourier interpolation 
 ----------------------
   input H( (n+gamma)/Nin ) => output H( (alpha+beta m)/Nout )
      n=0,Nin-1                           m can be beyound [0:1] 

   Note that 
      input H(1) = output H(1). 
   Please take account of the change of the unit vector 




   forward

   H( (n+gamma)/Nin ) , n=0,Nin-1, assuming periodic boundary condition

   h(k) = sum_{n=0,Nin-1} H( n+gamma) exp(-2pi i k (n+gamma)/Nin )
        = exp(-2pi i k gamma/Nin ) sum_{n=0,Nin-1} H(n+gamma) exp(-2pi i k n/Nin)
        = exp(-2 pi i k gamma/Nin) h^0(k)                  ---(1)
  where 
  h^0(k) = sum_{n=0,Nin-1} H(n+gamma) exp(-2pi i k n/Nin)            ---(2)

  h^0(k) is calculated using FFT. 
  h^0(k), n=0,Nin-1

  backward

  H( (alpha+beta n)Nout ), (alpha+beta n)/Nout can be beyond the range of [0:1]

  H'(alpha+beta n ) = sum_{k=0,Nin-1} h(k) exp(2 pi i k (alpha+ beta n )/Nout )
                   = sum_{k=0,Nin-1}  exp(-2 pi i k gamma/Nin) h^0(k) exp(2 pi i k (alpha+ beta n )/Nout )
                   = sum_{k=0,Nin-1}  h^0(k) exp( 2 pi i k ( (alpha+beta n)/Nout - gamma/Nin ) )  ---(3)

  Note, the transformation bove can not give all of (alpha+beta n)Nout correctly.
   In order to interpolate data, use the formula below  
     H(alpha+beta n ) = sum_{k=0,Nin/2}  h^0(k) exp( 2 pi i k ( (alpha+beta n)/Nout - gamma/Nin ) )
                      + sum_{k=-Nin/2,0}  h^0(k) exp( -2 pi i k ( (alpha+beta n)/Nout - gamma/Nin ) )
                                               ---(4)
     h^0(k-Nin) = h^0(k)


*********************************************************************************/



/******************************
  transformation (2) 
********************************/

#if 1


#ifdef fftw2 
#include <fftw.h>
#else
#include <fftw3.h> 
#endif


void zfft1(dcomplex *data,  /* size n_in */
	   int n_in, 
	   dcomplex *dataout,  /* output */
	   int isign0)
{
  fftw_complex *in, *out;
  fftw_plan p;

  int i;
  double scale;
  int isign;

  if (isign0>0) { isign=-1; }
  else { isign=1; }

#ifdef fftw2
  in  = (fftw_complex*)malloc(sizeof(fftw_complex)*n_in); 
  out = (fftw_complex*)malloc(sizeof(fftw_complex)*n_in); 
#else
  in  = fftw_malloc(sizeof(fftw_complex)*n_in);
  out = fftw_malloc(sizeof(fftw_complex)*n_in);
#endif

  for (i=0;i<n_in;i++) {

#ifdef fftw2
    c_re(in[i]) = data[i].r;
    c_im(in[i]) = data[i].i;
#else
    in[i][0]=data[i].r;
    in[i][1]=data[i].i;
#endif

  }

#ifdef fftw2
  p = fftw_create_plan(n_in, isign, FFTW_ESTIMATE);
  fftw_one(p,in,out); 
#else
  p = fftw_plan_dft_1d(n_in,in,out,isign,FFTW_ESTIMATE);
  fftw_execute(p);
#endif

  if (isign==-1) {
    scale = 1.0/( (double)n_in);

    for (i=0;i<n_in;i++) {

#ifdef fftw2
      dataout[i].r=c_re(out[i])*scale;
      dataout[i].i=c_im(out[i])*scale;
#else
      dataout[i].r=out[i][0]*scale;
      dataout[i].i=out[i][1]*scale;
#endif

    }
  }
  else {
    for (i=0;i<n_in;i++) {

#ifdef fftw2
      dataout[i].r=c_re(out[i]);
      dataout[i].i=c_im(out[i]);
#else
      dataout[i].r=out[i][0];
      dataout[i].i=out[i][1];
#endif

    }
  }


  fftw_destroy_plan(p);  

#ifdef fftw2
  free(out);
  free(in);
#else
  fftw_free(out);
  fftw_free(in);
#endif

}


#else
/* ESSL */
/* from essl manual */
/*GLOBAL STATEMENTS FOR ESSL ERROR HANDLING*/
#define __ESVERR
#include <essl.h>
extern int enotrm();
/*DECLARE THE VARIABLES*/

void zfft1(dcomplex *data, int n_in, dcomplex *dataout, int isign)
{
  int naux1, naux2;
  double *aux1, *aux2;
  int init;
  double scale;
  int inc1x, inc2x, inc1y, inc2y;
  int iret;
  int i;
  
  inc1x=1;
  inc2x=1;
  inc1y=1;
  inc2y=1;

  if ( n_in <= 2048 ) { naux1=20000; }
  else { naux1=20000+ 2.28*n_in; }
  naux2=naux1;
  aux1=(double*)malloc(sizeof(double)*naux1);
  aux2=(double*)malloc(sizeof(double)*naux2);

  init=1;
  if (isign==1) {
    scale = 1.0/( (double)n_in);
  }
  else {
    scale = 1.0;
  }
  iret=dcft(init,data,inc1x,inc2x,dataout,inc1y,inc2y, 
	    &n_in, 
	    1,
	    isign, 
	    scale, aux1,&naux1,aux2, &naux2);
  init=0;
  iret=dcft(init,data,inc1x,inc2x,dataout,inc1y,inc2y, &n_in,1, isign,
	    scale, aux1,&naux1,aux2, &naux2);

  free(aux2);
  free(aux1);

}
#endif


/*************************************************
    transformation (4) 
*************************************************/

/*
  isign=-1
  ->       z * exp(i gx) = (z.r, z.i ) * (cos(gx)+isin(gx))
  (z.r*cos(gx)-z.i*sin(gx), z.r*sin(gx)+z.i*cos(gx)  )
  isign=1
  ->   z * exp(i gx) = (z.r, z.i ) * (cos(gx)-isin(gx))
  (z.r*cos(gx)+z.i*sin(gx), -z.r*sin(gx)+z.i*cos(gx)  )
*/

void      zfourier1a(dcomplex *in,  /* input data */
                     int n_in,      /* size of in */
		     dcomplex *out,  /* output */
                     int n_out,     /* size of out */
		     int isign,
		     double alpha,double beta,double gamma) 
{
  
  double PI;
  double asign,bsign;
  int iout,k;
  double gx,kk,c,s;
  dcomplex z;

  PI = atan(1.0)*4.0;
  
  if (isign<0) {
    asign=-1.0;
    bsign=1.0;
  }
  else  {
    asign=1.0;
    bsign=-1.0;
  }
#if 0
  for (iout=0; iout<n_out;iout++) {
    z.r=0; z.i=0;
    gx= 2.0*PI*( (alpha+beta*(double)iout)/(double)n_out-gamma/(double)n_in ); 
    for (k=0;k<n_in;k++) {
      kk=(double)k*gx;
      c=cos(kk);
      s=sin(kk);
      z.r += in[k].r*c+in[k].i*s*asign;
      z.i += in[k].r*s*bsign+in[k].i*c;
    }
    out[iout]= z;
  }
#else
  if (n_in%2==0) { /* even */
      
    for (iout=0; iout<n_out;iout++) {
      z.r=0; z.i=0;
      gx= 2.0*PI*( (alpha+beta*(double)iout)/(double)n_out-gamma/(double)n_in ); 
      /* exp(ikx) */
      for (k=0;k<=n_in/2-1;k++) {
	kk=(double)k*gx;
	c=cos(kk);
	s=sin(kk);
	z.r += in[k].r*c+in[k].i*s*asign;
	z.i += in[k].r*s*bsign+in[k].i*c;
      }
      /* sum up as exp(-ikx) */
      for (k=n_in/2+1;k<n_in;k++) {
	kk=(double)(k-n_in)*gx;
	c=cos(kk);
	s=sin(kk);
	z.r += in[k].r*c+in[k].i*s*asign;
	z.i += in[k].r*s*bsign+in[k].i*c;
      }
      /* boundary */
      k=n_in/2;
      kk=(double)k*gx;
      c=cos(kk);
      s=sin(kk);
      z.r += (in[k].r*c+in[k].i*s*asign)*0.5;
      z.i += (in[k].r*s*bsign+in[k].i*c)*0.5;

      kk=(double)(k-n_in)*gx;
      c=cos(kk);
      s=sin(kk);
      z.r += (in[k].r*c+in[k].i*s*asign)*0.5;
      z.i += (in[k].r*s*bsign+in[k].i*c)*0.5;

      out[iout]= z;
    }


  } else { /* odd */

    for (iout=0; iout<n_out;iout++) {
      z.r=0; z.i=0;
      gx= 2.0*PI*( (alpha+beta*(double)iout)/(double)n_out-gamma/(double)n_in ); 
      /* exp(ikx) */
      for (k=0;k<=n_in/2;k++) {
	kk=(double)k*gx;
	c=cos(kk);
	s=sin(kk);
	z.r += in[k].r*c+in[k].i*s*asign;
	z.i += in[k].r*s*bsign+in[k].i*c;
      }
      /* sum up as exp(-ikx) */
      for (k=n_in/2+1;k<n_in;k++) {
	kk=(double)(k-n_in)*gx;
	c=cos(kk);
	s=sin(kk);
	z.r += in[k].r*c+in[k].i*s*asign;
	z.i += in[k].r*s*bsign+in[k].i*c;
      }
      out[iout]= z;
    }


  }

#endif
}



/*  din  size n_in[0]*n_in[1]*n_in[2], C order */
/*  dout size (n_out[3]+1)*n_out[1]*n_out[2], C order */
/*         n_out[0], n_out[1], n_out[2] is dimension for fourier transformation, 
	   you can find the explanation at top of this file .
	   But only 0:n_out[3] is stored */
/*   dout [0:n_out[3], 0:n_out[1]-1, 0:n_out[2]-1] */

void TRAN_FFT_Dinterpolation3D(
			       MPI_Comm comm1,
			       int n_in[3],
			       double *din, 
			       int n_out[4], 
			       double *dout, /* output */
			       double alpha[3],double beta[3],double gamma[3]  )
#define din_ref(i,j,k)    n_in[1]*n_in[2]*(i)+n_in[2]*(j)+k
#define dout_ref(i,j,k)  n_out[1]*n_out[2]*(i)+n_out[2]*(j)+k
{

  int n_max;

  int i;
  int i0,i1,i2;

  dcomplex *in;
  dcomplex *out;

  dcomplex *fft;
  int isign;
  int n_v[3];
  dcomplex ***zv;

  /*  FILE *fp; */

  static int numprocs,myid;
  int *idx0Start, *idx0End; 
  int *idx1Start, *idx1End;
  int *idx00Start, *idx00End;

  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  idx0Start=(int*)malloc(sizeof(int)*numprocs);
  idx0End=(int*)malloc(sizeof(int)*numprocs);
  idx1Start=(int*)malloc(sizeof(int)*numprocs);
  idx1End=(int*)malloc(sizeof(int)*numprocs);
  idx00Start=(int*)malloc(sizeof(int)*numprocs);
  idx00End=(int*)malloc(sizeof(int)*numprocs);


  TRAN_Distribute_Node(0,n_in[0]-1,numprocs,idx0Start, idx0End);
  TRAN_Distribute_Node(0,n_in[1]-1,numprocs,idx1Start, idx1End);
  TRAN_Distribute_Node(0,n_out[3],numprocs,idx00Start, idx00End);


  for (i=0;i<3;i++) {
    n_v[i] = ( n_in[i] > n_out[i] )? n_in[i]:n_out[i];
  }
  zv = (dcomplex***)malloc(sizeof(dcomplex**)*n_v[0]);
  for (i0=0;i0<n_v[0];i0++) {
    zv[i0]=(dcomplex**)malloc(sizeof(dcomplex*)*n_v[1]);
    for (i1=0;i1<n_v[1];i1++) {
      zv[i0][i1]=(dcomplex*)malloc(sizeof(dcomplex)*n_v[2]);
    }
  }


  n_max=0;
  for (i=0;i<3;i++) {
    if (n_max < n_in[i]) n_max = n_in[i];
  }
  for (i=0;i<3;i++) {
    if (n_max < n_out[i]) n_max = n_out[i];
  }

  in  = (dcomplex*)malloc(sizeof(dcomplex)*n_max);
  out = (dcomplex*)malloc(sizeof(dcomplex)*n_max);



  /* forward FFT */
  isign=1;

  /* dim(x,y,z) */

  /*MPI parallel   global: i0  0:n_in[0]-1 */
  /*MPI variable: ID=myid; */
  /*MPI local: i0:idx0Start[ID]:idx0End[ID] */
  if ( idx0Start[myid]>=0) {
    for (i0=idx0Start[myid];i0<=idx0End[myid];i0++) {
      for (i1=0;i1<n_in[1];i1++) {
	for (i2=0;i2<n_in[2];i2++) {
	  in[i2].r = din[din_ref(i0,i1,i2)];
	  in[i2].i = 0.0;
	}
	zfft1( in,n_in[2],out,isign);
	for (i2=0;i2<n_in[2];i2++) {
	  zv[i0][i1][i2] = out[i2];
	}
      }
    }
  }

  /* zv(Lx,y,kz) */


  /*MPI parallel   global: i0  0:n_in[0]-1 */
  /*MPI variable: ID=myid; */
  /*MPI local: i0:idx0Start[ID]:idx0End[ID] */
  if ( idx0Start[myid]>=0) {
    for (i0=idx0Start[myid];i0<=idx0End[myid];i0++) {
      for (i2=0;i2<n_in[2];i2++) {
	for (i1=0;i1<n_in[1];i1++) {
	  in[i1]= zv[i0][i1][i2];
	}
	zfft1( in,n_in[1],out,isign);
	for (i1=0;i1<n_in[1];i1++) {
	  zv[i0][i1][i2] = out[i1];
	}
      }
    }
  }

  /* zv(Lx,ky,kz) */
  { int ID;
  for (ID=0;ID<numprocs;ID++) {
    if ( idx0Start[ID]>=0 ) {
      for (i0=idx0Start[ID];  i0<=idx0End[ID]; i0++) {
	for (i1=0;i1<n_in[1];i1++) {
	  /*MPI Bcast: zv[i0][i1]  size: n_in[2]  */
	  MPI_Bcast(zv[i0][i1], n_in[2]*2, MPI_DOUBLE, ID,comm1);
	}
      }
    }
  }
  }

  /* zv(x,ky,kz) */

  if ( idx1Start[myid]>=0 ) {
    for (i2=0;i2<n_in[2];i2++) {
      /*MPI parallel */
      /*MPI global: i1   0:n_in[1]-1 */
      /*MPI lobal: i1   idx1Start[ID]:idx1End[ID] */
      for (i1=idx1Start[myid];i1<=idx1End[myid];i1++) {
	for (i0=0;i0<n_in[0];i0++) {
	  in[i0]= zv[i0][i1][i2];
	}
	zfft1( in,n_in[0],out,isign);
	for (i0=0;i0<n_in[0];i0++) {
	  zv[i0][i1][i2] = out[i0];
	}
      }
    }
  }

  /* zv(kx,Lky,kz) */

#define FOURIER 
  /* inverse fourier transformation  */
  isign=-1;

  /* zv(kx,Lky,kz) */

  if ( idx1Start[myid]>=0 ) {
    for (i0=0;i0<n_in[0];i0++) {
      /*MPI parallel */
      /*MPI global: i1   0:n_in[1]-1 */
      /*MPI lobal: i1   idx1Start[ID]:idx1End[ID] */
      for (i1=idx1Start[myid];i1<=idx1End[myid];i1++) {
	for (i2=0;i2<n_in[2];i2++) {
	  in[i2]= zv[i0][i1][i2];
	}
	zfourier1a( in,n_in[2], out,n_out[2], isign, alpha[2],beta[2],gamma[2]);
	for (i2=0;i2<n_out[2];i2++) {
	  zv[i0][i1][i2] = out[i2];

	}
      }
    }
  }
  /* zv(kx,Lky,Z)*/

  if ( idx1Start[myid]>=0 ) {
    for (i2=0;i2<n_out[2];i2++) {
      /*MPI parallel */
      /*MPI global: i1   0:n_in[1]-1 */
      /*MPI lobal: i1   idx1Start[ID]:idx1End[ID] */
      for (i1=idx1Start[myid];i1<=idx1End[myid];i1++) {
	for (i0=0;i0<n_in[0];i0++) {
	  in[i0]= zv[i0][i1][i2];
	}
	zfourier1a( in,n_in[0],out,n_out[0],isign,alpha[0],beta[0],gamma[0]);
	for (i0=0;i0<=n_out[3];i0++) {
	  zv[i0][i1][i2] = out[i0];
	}
      }
    }
  }

  /* zv(X,Lky,Z) */
  { int ID;
  for (ID=0;ID<numprocs;ID++) {
    for (i0=0;  i0<=n_out[3]; i0++) {
      /*MPI parallel */
      /*MPI global: i1   0:n_in[1]-1 */
      /*MPI lobal: i1   idx1Start[ID]:idx1End[ID] */
      if ( idx1Start[ID]>=0 ) {
	for (i1=idx1Start[ID];i1<=idx1End[ID];i1++) {
	  /*MPI Bcast: zv[i0][i1]  size: n_in[2]  */
	  MPI_Bcast(zv[i0][i1], n_in[2]*2, MPI_DOUBLE, ID,comm1);
	}
      }
    }
  }
  }

  /* zv(X,ky,Z) */

  if ( idx00Start[myid] >=0 ) {
    for (i0=idx00Start[myid]; i0<=idx00End[myid]; i0++) {
      for (i2=0; i2<n_out[2]; i2++) {

	for (i1=0; i1<n_in[1]; i1++) {
	  in[i1]= zv[i0][i1][i2];
	}

	zfourier1a( in,n_in[1],out,n_out[1],isign,alpha[1],beta[1],gamma[1]);

	for (i1=0; i1<n_out[1]; i1++) {

	  dout[ dout_ref(i0,i1,i2) ] = out[i1].r;

	  /*
          printf("myid=%2d dout_ref(i0,i1,i2)=%3d\n",myid,dout_ref(i0,i1,i2));
	  */
          
	}
      }
    }
  }

  /* dout(LX,Y,Z) */
  { int ID;
  for (ID=0;ID<numprocs;ID++) {
    if ( idx00Start[ID]>=0 )  {
      for (i0=idx00Start[ID];  i0<=idx00End[ID]; i0++) {
        /*MPI Bcast: dout_ref(i0,0,0)  size: n_out[1]*n_out[2]  */
        MPI_Bcast( &dout[ dout_ref(i0,0,0) ], n_out[1]*n_out[2], MPI_DOUBLE, ID, comm1);
      }
    }
  }
  }

  /* dout(X,Y,Z) */

  if (myid==Host_ID && print_stdout){

    printf("the last of dout= %d %d %d, %d %le\n",
	   n_out[0], n_out[1], n_out[2],
	   dout_ref(n_out[3],n_out[1]-1,n_out[2]-1), 
	   dout[dout_ref(n_out[3],n_out[1]-1,n_out[2]-1)]);
  }

  free(idx0Start);
  free(idx0End);
  free(idx1Start);
  free(idx1End);
  free(idx00Start);
  free(idx00End);



  free(out);
  free(in);

  for (i0=0;i0<n_v[0];i0++) {
    for (i1=0;i1<n_v[1];i1++) {
      free(zv[i0][i1]);
    }
    free(zv[i0]);
  }
  free(zv);



}



#if 0

#define n0  2
#define n1  8
#define n2  16
#define nout0 2
#define nout1 8
#define nout2 16


#define in_ref(i,j,k) in[n_in[1]*n_in[2]*(i)+n_in[2]*(j)+k]
#define out_ref(i,j,k) out[n_out[1]*n_out[2]*(i)+n_out[2]*(j)+k]

main()
{
  double alpha[3],beta[3],gamma[3];
  double x0,x1,x2;

  int n_in[3], n_out[4];
  double *in=NULL;
  double *out=NULL; 

  int i0,i1,i2,i2b;

  FILE *fin, *fout;

  n_in[0]= n0;
  n_in[1]= n1;
  n_in[2]= n2;
  n_out[0]=nout0;
  n_out[1]=nout1;
  n_out[2]=nout2;
  n_out[3]=nout2-3;

  in=(double*)malloc(sizeof(double)*n_in[0]*n_in[1]*n_in[2]);
  out=(double*)malloc(sizeof(double)*n_out[0]*n_out[1]*n_out[2]);
   

  alpha[0]=0.5;
  alpha[1]=0.5;
  alpha[2]=0.5*0.8;
  beta[0]=1.0;
  beta[1]=1.0;
  beta[2]=1.0*0.8;
  gamma[0]=0.5;
  gamma[1]=0.5;
  gamma[2]=0.5;

  for (i0=0;i0<n_in[0];i0++) {
    x0= gamma[0]+i0-n_in[0]/2;
    for (i1=0;i1<n_in[1];i1++) {
      x1 =gamma[1]+i1-n_in[1]/2;
      for (i2=0;i2<n_in[2];i2++) {
        x2= gamma[2]+i2-n_in[2]/2;
	in_ref(i0,i1,i2)= x0*x0+x1*x1+x2*x2;
      }
    }
  }

  fin= fopen("in","w");
  for (i0=0;i0<n_in[0];i0++) {
    x0=(gamma[0]+i0)/n_in[0]; 
    for (i1=0;i1<n_in[1];i1++) {
      x1=(gamma[1]+i1)/n_in[1];
#if 0
      for (i2=0;i2<n_in[2];i2++) {
	x2= (gamma[2]+i2)/n_in[2];
#else
	for (i2b=0;i2b<=n_in[2];i2b++) {
	  x2= (gamma[2]+i2b)/n_in[2];
	  i2 = i2b % n_in[2];
#endif
	  /*  printf("%lf %lf ",in[i0][i1][i2].r, in[i0][i1][i2].i); */
	  fprintf(fin,"%lf %lf %lf %lf\n",x0,x1,x2,in_ref(i0,i1,i2));
	}
	fprintf(fin,"\n");
      }
    }
    fprintf(fin,"\n");
    fclose(fin);


    TRAN_FFT_Dinterpolation3D( n_in, in, n_out, out, alpha,beta,gamma);


    fout= fopen("out","w");
    for (i0=0;i0<n_out[0];i0++) {
      x0=(alpha[0]+beta[0]*i0)/n_out[0];
      for (i1=0;i1<n_out[1];i1++) {
	x1=(alpha[1]+beta[1]*i1)/n_out[1];
	for (i2=0;i2<n_out[2];i2++) {
	  x2=(alpha[2]+beta[2]*i2)/n_out[2];
	  /*printf("%lf %lf ",out[i0][i1][i2].r, out[i0][i1][i2].i);*/
	  fprintf(fout,"%lf %lf %lf %lf\n",x0,x1,x2,out_ref(i0,i1,i2));

	}
	fprintf(fout,"\n");
      }
    }
    fclose(fout);

    free(out);
    free(in);

  }

#endif
