

#include <math.h>
#include <algorithm>

void dpotrf(double **A, int msize, double &det);
void dpotrfU(double **A, int msize, double &det);
void dpotrfInfo(double **A, int msize, double &det, int &info);
double *dpotrf_ctof(double **in, int rows, int cols);
void dpotrf_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dpotrf_(char *UPLO, int *msize, double *a, int *lda, int *info);

