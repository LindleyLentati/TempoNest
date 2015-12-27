#include <math.h>
#include <stdio.h>
#include <algorithm>
#include "dgesvd.h"
//#include "mkl_lapacke.h"
//#include "mkl.h"

//#define min(a,b) ((a)>(b)?(b):(a))


void vector_dgesvd(double *A, int M, int N)
{

	char jobu, jobvt;
	int m=M, n=N ,lda=N, ldu=M, ldvt=N;
	double *a, *s, *u, *vt;
	int lwork, info;
	double wkopt;
	double *work;
	int *iwork = new int[8*N];

	vt = new double [N*N];
	s = new double[N];


	jobu = 'O'; /* Specifies options for computing U.
		 A: all M columns of U are returned in array U;
		 S: the first min(m,n) columns of U (the left
		    singular vectors) are returned in the array U;
		 O: the first min(m,n) columns of U (the left
		    singular vectors) are overwritten on the array A;
		 N: no columns of U (no left singular vectors) are
		    computed. */
	lwork = -1;
	dgesdd_( "O", &m, &n, A, &ldu, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info );
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Compute SVD */
	dgesdd_( "O", &m, &n, A, &ldu, s, u, &ldu, vt, &ldvt, work,  &lwork, iwork, &info );


	delete vt;
	delete s;
	delete iwork;
	delete work;
}
void dgesvd(double **A, int M, int N, double *S, double **U, double **VT)
{

	char jobu, jobvt;
	int m=M, n=N ,lda=N, ldu=M, ldvt=N;
	double *a, *s, *u, *vt;

	int lwork, info;
	double wkopt;
	double *work;
	int *iwork = new int[8*N];



  a = new double [M*N];
  u = new double [M*M];
  vt = new double [N*N];

  
  a = dgesvd_ctof(A, lda, n); /* Convert the matrix A from double pointer
			  C form to single pointer Fortran form. */

  
	lwork = -1;
	dgesdd_( "O", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info );
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Compute SVD */
	dgesdd_( "O", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work,  &lwork, iwork, &info );


	dgesvd_ftoc(u, U, ldu, ldu);
	dgesvd_ftoc(vt, VT, ldvt, n);

  delete a;
  delete u;
  delete vt;
  delete iwork;
  delete work;
}


double* dgesvd_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dgesvd_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
