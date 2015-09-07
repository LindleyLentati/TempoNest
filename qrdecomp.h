

#include <math.h>
#include <algorithm>

void TNqrsolve(double **aIn, double *bIn, double *xOut, int msize, double &det, int &info);
void vector_TNqrsolve(double *aIn, double *bIn, double *xOut, int msize, double &det, int &info);
void qrsolveInfo(double **A, double *B, int msize, int &info);
double *qr_ctof(double **in, int rows, int cols);
void qr_ftoc(double *in, double **out, int rows, int cols);


extern "C" void dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
extern "C" void dormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);
extern "C" void dtrtrs_(char *uplo, char *trans, char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
extern "C" void dgeqrf_(int* M, int* N, double* A, int* LDA, double* TAU, double* WORK, int* LWORK, int* INFO );
extern "C" void dtrrfs_(char *uplo, char *trans, char *diag, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, double *X,  int *ldx, double *Ferr, double *Berr, double *work, int *work2, int *info);
