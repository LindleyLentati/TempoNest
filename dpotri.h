

#include <math.h>
#include <algorithm>

void dpotri(double **A, int msize);

double *dpotri_ctof(double **in, int rows, int cols);
void dpotri_ftoc(double *in, double **out, int rows, int cols);

extern "C" void dpotri_(char *UPLO, int *msize, double *a, int *lda, int *info);


