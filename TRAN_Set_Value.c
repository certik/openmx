#include <stdio.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"

/*
 * used to initialize A
 */
void TRAN_Set_Value_double(dcomplex *A, int n, double a, double b)
{
  int i;
  for(i=0; i<n; i++) {
    A[i].r = a;
    A[i].i = b;
  }
}
