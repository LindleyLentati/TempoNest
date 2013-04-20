#include <math.h>
#include <algorithm>
#include "dpotri.h"

void dpotri(double **A, int msize)
{

	int info;
	double *a;
	char UPLO='L';

	a = dpotri_ctof(A, msize, msize); 
	
	dpotri_(&UPLO, &msize, a, &msize, &info);
	dpotri_ftoc(a, A, msize, msize);

	for(int i=0;i<msize;i++){
		for(int j=0;j<i;j++){
			A[j][i]=A[i][j];
		}
	}

  
  delete a;

}


double* dpotri_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dpotri_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}