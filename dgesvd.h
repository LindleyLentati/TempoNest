#include <math.h>
#include <algorithm>

void dgesvd(double **A, int m, int n, double *S, double **U, double **VT);
void vector_dgesvd(double *A, int M, int N);

void dgesvd_2Dto1D(double **in, double *out, int rows, int cols);
void dgesvd_1Dto2D(double *in, double **out, int rows, int cols);

extern "C" void dgesdd_(char *JOBZ, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *IWORK, int *INFO);
//
double* dgesvd_ctof(double **in, int rows, int cols);
void dgesvd_ftoc(double *in, double **out, int rows, int cols);
//extern "C" void dgesvd_( char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *S, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);





