#include "config.h"
extern char useT2accel;


#ifdef HAVE_LAPACK
#define ACCEL_UINV
void accel_uinv(double* _m, int n);
#endif

#ifdef HAVE_BLAS
#define ACCEL_MULTMATRIX
void accel_multMatrix(double* m1,double* m2, int ndata,int ndata2,int npol, double* out);
#endif
