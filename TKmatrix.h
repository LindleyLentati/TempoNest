//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

#ifndef __TKmatrix_h
#define __TKmatrix_h


/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/
#ifdef __cplusplus
extern "C" {
#endif
void TKmultMatrix_sq( double **idcm, double **u,int ndata,int npol,double **uout);
void TKmultMatrixVec_sq( double **idcm, double *b,int ndata,double *bout);
void TKmultMatrix( double **idcm, double **u,int ndata,int ndata2,int npol,double **uout);
void TKmultMatrixVec( double **idcm, double *b,int ndata,int ndata2,double *bout);

double** malloc_uinv(int n);
double **malloc_blas(int n,int m);
void free_blas(double** matrix);
void free_uinv(double** uinv);
int get_blas_rows(double** uinv);
int get_blas_cols(double** uinv);
#ifdef __cplusplus
}
#endif


#ifdef __Tempo2_h
void free_2dLL(longdouble** m);
longdouble** malloc_2dLL(int rows,int cols);
#endif
#endif
