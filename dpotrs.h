

#include <math.h>
#include <algorithm>

void dpotrs(double **A, double *B, int msize);
void dpotrsInfo(double **A, double *B, int msize, int &info);
double *dpotrs_ctof(double **in, int rows, int cols);
void dpotrs_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dpotrs_(char *UPLO, int *msize, int *rhs, double *a, int *lda, double *b, int *ldb, int *info);

