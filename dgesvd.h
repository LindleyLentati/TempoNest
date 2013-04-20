

#include <math.h>
#include <algorithm>

void dgesvd(double **A, int m, int n, double *S, double **U, double **VT);

double *dgesvd_ctof(double **in, int rows, int cols);
void dgesvd_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
			double *a, int *lda, double *s, double *u,
			int *ldu, double *vt, int *ldvt, double *work,
			int *lwork, int *info);


