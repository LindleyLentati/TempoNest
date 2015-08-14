#ifdef __cplusplus
extern "C" {
#endif

void cholesky_formUinvL(double **uinv, __float128 **m,int np, double &det, int &info);
void TNcholDecompositionL(__float128 **a, int n, double &det, int &info);
void T2cholDecompositionL(__float128 **a, int n,  __float128 *p, int &info);
#ifdef __cplusplus
}
#endif


