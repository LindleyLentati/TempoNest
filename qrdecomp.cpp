#include <math.h>
#include <algorithm>
#include "qrdecomp.h"
#include <stdio.h>

    int geqrf(int m, int n, double* A, int lda, double *tau) {
        int info=0;
        int lwork=-1;
        double iwork;
        dgeqrf_(&m, &n, A, &lda, tau, &iwork, &lwork, &info);
        lwork = (int)iwork;
        double* work = new double[lwork];
        dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
        delete[] work;
        return info;
    }

    int ormqr(char side, char trans, int m, int n, int k, double *A, int lda, double *tau, double* C, int ldc) {

        int info=0;
        int lwork=-1;
        double iwork;
        dormqr_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, &iwork, &lwork, &info);
        lwork = (int)iwork;
        double* work = new double[lwork];
        dormqr_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, &info);
        delete[] work;
        return info;
    }

int trtrs(char uplo, char trans, char diag, int n, int nrhs, double* A, int lda, double* B, int ldb) {

        int info = 0;
        dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, A, &lda, B, &ldb, &info);
        return info;
    }

int trrfs(char uplo, char trans, char diag, int n, int nrhs, double *A, int lda, double *B, int ldb, double *X, int ldx){


	int info = 0;
	double *work = new double[3*n];
	int *work2 = new int[n];
	double *Ferr = new double[n];
	double *Berr = new double[n];
	printf("Here \n");
	dtrrfs_(&uplo, &trans, &diag, &n, &nrhs, A, &lda, B, &ldb, X, &ldx, Ferr, Berr, work, work2, &info);


	for(int i =0; i < n; i ++){
		printf("Ferr %i %g %g \n", i, Ferr[i], Ferr[i]*X[i]);
	}

	delete[] work;
	delete[] work2;
	delete[] Ferr;
	delete[] Berr;

	return info;
	}



void vector_TNqrsolve(double *aIn, double *bIn, double *xOut, int msize, double &det, int &info){

	double *a = new double[msize*msize];
	for(int i =0 ; i < msize*msize; i++){	
		a[i] = aIn[i];
	}

	double *Atemp = new double [msize*msize];

	double *b = new double[msize];
	for(int i =0 ; i < msize; i++){
		b[i]=bIn[i];
	}

	double* tau = new double[msize];

	int tempinfo=0;
	tempinfo=geqrf(msize, msize, a, msize, tau);

	if(tempinfo != 0)info=tempinfo;
	for(int i =0 ; i < msize*msize; i++){
                Atemp[i] = a[i];
        }

	det=0;
	for(int i =0; i < msize; i++){
		det += log(fabs(Atemp[i + i*msize]));
	}

	tempinfo=0;
	tempinfo=ormqr('L', 'T', msize, 1, msize, a, msize, tau, b, msize);
	if(tempinfo != 0)info=tempinfo;
	tempinfo=0;
	tempinfo=trtrs('U', 'N', 'N', msize, 1, a, msize, b, msize);
	if(tempinfo != 0)info=tempinfo;
	for(int i =0; i < msize; i++){
		xOut[i] = b[i];
	}
	//tempinfo = trrfs('U', 'N', 'N', msize, 1, a, msize, bIn, msize, xOut, msize);

	delete[] a;
	delete[] tau;
	delete[] b;
	delete[] Atemp;
}

void TNqrsolve(double **aIn, double *bIn, double *xOut, int msize, double &det, int &info){

	
	double *a;
	a = qr_ctof(aIn, msize, msize);

	double **Atemp = new double *[msize];
	for(int i = 0 ; i <msize; i++){
		Atemp[i] = new double[msize];
	}

	double *b = new double[msize];
	for(int i =0 ; i < msize; i++){
		b[i]=bIn[i];
	}

	double* tau = new double[msize];

	int tempinfo=0;
	tempinfo=geqrf(msize, msize, a, msize, tau);

	if(tempinfo != 0)info=tempinfo;
	qr_ftoc(a, Atemp, msize, msize);
	det=0;
	for(int i =0; i < msize; i++){
		det += log(fabs(Atemp[i][i]));
	}

	tempinfo=0;
	tempinfo=ormqr('L', 'T', msize, 1, msize, a, msize, tau, b, msize);
	if(tempinfo != 0)info=tempinfo;
	tempinfo=0;
	tempinfo=trtrs('U', 'N', 'N', msize, 1, a, msize, b, msize);
	if(tempinfo != 0)info=tempinfo;
	for(int i =0; i < msize; i++){
		xOut[i] = b[i];
	}
	//tempinfo = trrfs('U', 'N', 'N', msize, 1, a, msize, bIn, msize, xOut, msize);

	delete[] a;
	delete[] tau;
	delete[] b;
	for(int i = 0 ; i <msize; i++){	
		delete[] Atemp[i];
	}
	delete[] Atemp;
}
/*

void qrsolveInfo(double **A, double *B, int msize, int &info)
{


	double *a;
	char UPLO='U';
	char side='L';
	char trans='T';
	char notrans='N';
	char diag='N';
	int NRHS=1;
	int LWORK = 2*msize + (msize+1)*msize;
	a = qr_ctof(A, msize, msize); 


	double *tau = new double[msize];
	int *jpvt = new int[msize];
	double *work = new double[LWORK];
	double *C = new double[msize];
	

	for(int i =0; i < msize; i++){
		tau[i]=0;
		jpvt[i]=0;
	}

	info=0;
	dgeqp3_(&msize, &msize, a, &msize, jpvt, tau, work, &LWORK, &info);
	info=0;
	dormqr_(&side, &trans, &msize, &NRHS, &msize, a, &msize, tau, B, &msize, work, &LWORK, &info);
	info=0;
	dtrtrs_(&UPLO, &notrans, &diag, &msize, &NRHS, a, &msize, B, &msize, &info);
	
	qr_ftoc(a, A, msize, msize);


 	delete[] tau;
	delete[] jpvt;  
	delete[] work;
	delete[] C;
	delete a;

}
*/

double* qr_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void qr_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
