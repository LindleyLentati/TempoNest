#include <math.h>
#include <algorithm>
#include "dpotrf.h"
#include <stdio.h>
void dpotrf(double **A, int msize, double &det)
{

	int info;
	double *a;
	char UPLO='L';

	a = dpotrf_ctof(A, msize, msize); 
	
	dpotrf_(&UPLO, &msize, a, &msize, &info);
	dpotrf_ftoc(a, A, msize, msize);

	det=0;
	for(int i=0;i<msize;i++){
		det+=log(A[i][i]);
	}

	det=det*2;

  	//printf("info: %i \n", info);
  delete a;

}

void dpotrfU(double **A, int msize, double &det)
{

	int info;
	double *a;
	char UPLO='U';

	a = dpotrf_ctof(A, msize, msize); 
	
	dpotrf_(&UPLO, &msize, a, &msize, &info);
	dpotrf_ftoc(a, A, msize, msize);

	det=0;
	for(int i=0;i<msize;i++){
		det+=log(A[i][i]);
	}

	det=det*2;

  
  delete a;

}


void dpotrfInfo(double **A, int msize, double &det, int &info)
{

        double *a;
        char UPLO='L';

        a = dpotrf_ctof(A, msize, msize);

        dpotrf_(&UPLO, &msize, a, &msize, &info);
        dpotrf_ftoc(a, A, msize, msize);

        det=0;
        for(int i=0;i<msize;i++){
		//printf("Det %i %g \n", i, A[i][i]);
                det+=log(A[i][i]);
        }
		//printf("in func info %i \n", info);
        det=det*2;


  delete a;

}



double* dpotrf_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dpotrf_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
