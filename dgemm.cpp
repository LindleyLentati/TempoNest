#include <math.h>
#include <algorithm>
#include "dgemm.h"

void dgemm(double **A, double **B,double **C,int rowa, int cola, int rowb, int colb, char AT, char BT)
{

  int M,N,K;
  double *a, *b, *c;

	double alpha=1;
	double beta=0;


	if(AT == 'N'){
		M=rowa;
		K=cola;
	}
	else if(AT == 'T'){
		M=cola;
		K=rowa;
	}

	if(BT == 'N'){
		N=colb;
	}
	else if(BT == 'T'){
		N=rowa;
	}	
/*
(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)*/


	a = dgemm_ctof(A, rowa, cola); 
	b = dgemm_ctof(B, rowb, colb);
	c = new double[M*N];
	
	dgemm_(&AT, &BT, &M, &N, &K, &alpha, a, &rowa,b, &rowb, &beta, c, &M);
	dgemm_ftoc(c, C, M, N);

  
  delete a;
  delete b;
  delete c;

}


double* dgemm_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dgemm_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
