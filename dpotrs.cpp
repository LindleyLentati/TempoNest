#include <math.h>
#include <algorithm>
#include "dpotrs.h"

void dpotrs(double **A, double *B, int msize)
{

	int info;
	double *a;
	char UPLO='L';
	int NRHS=1;
	a = dpotrs_ctof(A, msize, msize); 
	
	dpotrs_(&UPLO, &msize, &NRHS, a, &msize, B, &msize, &info);
	dpotrs_ftoc(a, A, msize, msize);
  
  delete a;

}

void dpotrsInfo(double **A, double *B, int msize, int &info)
{


	double *a;
	char UPLO='L';
	int NRHS=1;
	a = dpotrs_ctof(A, msize, msize); 
	
	dpotrs_(&UPLO, &msize, &NRHS, a, &msize, B, &msize, &info);
	dpotrs_ftoc(a, A, msize, msize);
  
  delete a;

}


double* dpotrs_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dpotrs_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
