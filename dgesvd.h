#include <math.h>
#include <algorithm>

void dgesvd(double **A, int m, int n, double *S, double **U, double **VT);
void vector_dgesvd(double *A, int M, int N);

void dgesvd_2Dto1D(double **in, double *out, int rows, int cols);
void dgesvd_1Dto2D(double *in, double **out, int rows, int cols);




