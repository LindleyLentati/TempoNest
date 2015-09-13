#include "config.h"
extern char useT2accel;


#ifdef HAVE_LAPACK
#define ACCEL_UINV
#define ACCEL_LSQ
int accel_uinv(double* _m, int n);
double accel_lsq_qr(double** dm, double* data, double* oparm, int ndata, int nparam, double** Ocvm);
#endif

#ifdef HAVE_BLAS
#define ACCEL_MULTMATRIX
void accel_multMatrixVec(double* m1,double* v, int ndata,int npol, double* out);
void accel_multMatrix(double* m1,double* m2, int ndata,int ndata2,int npol, double* out);
#endif
