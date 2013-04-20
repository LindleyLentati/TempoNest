#include <math.h>
#include <algorithm>

void dgemm(double **A, double **B,double **C,int rowa, int cola, int rowb, int colb, char AT, char BT);

double *dgemm_ctof(double **in, int rows, int cols);
void dgemm_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dgemm_(char *jobu, char *jobvt, int *m, int *n,
			int *k, double *alpha, double *a, int *lda,
			double *b, int *ldb, double *beta, double *c,
			int *ldc);
