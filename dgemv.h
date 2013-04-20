
#include <math.h>
#include <algorithm>

void dgemv(double **A, double *vecin,double *vecout,int rowa, int cola, char AT);

double *dgemv_ctof(double **in, int rows, int cols);
void dgemv_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dgemv_(char *jobu, int *m, int *n,
			double *alpha, double *a, int *lda,
			double *x, int *incx, double *beta, double *y, int *incy);


