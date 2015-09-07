#include <math.h>
#include <stdio.h>
#include <algorithm>
#include "dgesvd.h"
#include "mkl_lapacke.h"
#include "mkl.h"

#define min(a,b) ((a)>(b)?(b):(a))


void vector_dgesvd(double *A, int M, int N)
{

  char jobu, jobvt;
  MKL_INT m=M, n=N ,lda=N, ldu=M, ldvt=N;
  double *a, *s, *u, *vt;
//  a = new double [M*N];
//  u = new double [M*M];
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

  jobvt = 'N'; /* Specifies options for computing VT.
		  A: all N rows of V**T are returned in the array
		     VT;
		  S: the first min(m,n) rows of V**T (the right
		     singular vectors) are returned in the array VT;
		  O: the first min(m,n) rows of V**T (the right
		     singular vectors) are overwritten on the array A;
		  N: no rows of V**T (no right singular vectors) are
		     computed. */

//	dgesvd_2Dto1D(A, a, lda, n);

	MKL_INT info;
	//info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, jobu, jobvt, m, n, a, lda,
        //                S, u, ldu, vt, ldvt, superb );
	info = LAPACKE_dgesdd( LAPACK_ROW_MAJOR, 'O', m, n, A, lda, s, u, ldu, vt, ldvt);

//	dgesvd_1Dto2D(u, U, ldu, ldu);
//	dgesvd_1Dto2D(vt, VT, ldvt, n);
  
//  delete a;
//delete u;
  delete vt;
  delete s;
  //delete work;
}
void dgesvd(double **A, int M, int N, double *S, double **U, double **VT)
{

  char jobu, jobvt;
  MKL_INT m=M, n=N ,lda=N, ldu=M, ldvt=N;
  double *a, *s, *u, *vt;

  a = new double [M*N];
  u = new double [M*M];
  vt = new double [N*N];

  

  jobu = 'A'; /* Specifies options for computing U.
		 A: all M columns of U are returned in array U;
		 S: the first min(m,n) columns of U (the left
		    singular vectors) are returned in the array U;
		 O: the first min(m,n) columns of U (the left
		    singular vectors) are overwritten on the array A;
		 N: no columns of U (no left singular vectors) are
		    computed. */

  jobvt = 'A'; /* Specifies options for computing VT.
		  A: all N rows of V**T are returned in the array
		     VT;
		  S: the first min(m,n) rows of V**T (the right
		     singular vectors) are returned in the array VT;
		  O: the first min(m,n) rows of V**T (the right
		     singular vectors) are overwritten on the array A;
		  N: no rows of V**T (no right singular vectors) are
		     computed. */

	dgesvd_2Dto1D(A, a, lda, n);

	MKL_INT info;
	//info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, jobu, jobvt, m, n, a, lda,
        //                S, u, ldu, vt, ldvt, superb );
	info = LAPACKE_dgesdd( LAPACK_ROW_MAJOR, 'S', m, n, a, lda,
                        S, u, ldu, vt, ldvt);

	dgesvd_1Dto2D(u, U, ldu, ldu);
	dgesvd_1Dto2D(vt, VT, ldvt, n);
  
  delete a;
  delete u;
  delete vt;
  //delete work;
}
// 		LWMAX=100000
// 		allocate(A(n,nf),U(n,n),Sig(nf),V(nf,nf),IWORK(8*nf),WORK(LWMAX))
// 		A=0
// 		A=onedesign
// 		U=0
// 		Sig=0
// 		V=0
// 		Work=0
// 		IWork=0
// 		LWORK=-1
// 		Info=0
// 		call DGESVD('A','A',n,nf,A,n,Sig,U,n,V,nf,Work,LWork,Info)

void dgesvd_2Dto1D(double **in, double *out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i + j*rows] = in[i][j];
}


void dgesvd_1Dto2D(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}

