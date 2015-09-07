#include <math.h>
#include <algorithm>
//#include "dgesvd.h"
#include "mkl_lapacke.h"
#include "mkl.h"

#define min(a,b) ((a)>(b)?(b):(a))

void dgesvd(double **A, int m, int n, double *S, double **U, double **VT)
{
  char jobu, jobvt;
  int lda, ldu, ldvt, lwork;
  double *a, *u, *vt, *work;

  int minmn, maxmn;

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

  lda = m; // The leading dimension of the matrix a.
#if 0
  a = dgesvd_ctof(A, lda, n); /* Convert the matrix A from double pointer
			  C form to single pointer Fortran form. */
#endif

  ldu = m;


	maxmn = m;
	minmn = n;

	ldu = m; // Left singular vector matrix
	u = new double[ldu*ldu];
	
	ldvt = n; // Right singular vector matrix
	vt = new double[ldvt*n];

	double superb[min(m,n)-1];
/*
	int LMAX=100000;
	
	work = new double[LMAX];
	lwork = -1; // Set up the work array, larger than needed.

// 	printf("parm 11 %i %i\n",ldu,ldvt);
	dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, S, u,&ldu, vt, &ldvt, work, &lwork, &info);
	
	lwork = std::min(LMAX,int(work[0]));
	
	dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, S, u,&ldu, vt, &ldvt, work, &lwork, &info);
// 	printf("parm 11 out %i %i\n",ldu,ldvt);
*/
	MKL_INT info;
	//info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, jobu, jobvt, m, n, a, lda,
        //                S, u, ldu, vt, ldvt, superb );
	info = LAPACKE_dgesdd( LAPACK_ROW_MAJOR, 'S', m, n, a, lda,
                        S, u, ldu, vt, ldvt);

	//dgesvd_ftoc(u, U, ldu, ldu);
	//dgesvd_ftoc(vt, VT, ldvt, n);
  
  delete a;
  delete u;
  delete vt;
  delete work;
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

