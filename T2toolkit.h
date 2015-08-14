//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

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



/* ************************************************************** */
/* Routines to convert doubles to floats                          */
/* ************************************************************** */

void TKconvertFloat1(double *x,float *ox,int n);
void TKconvertFloat2(double *x,double *y,float *ox,float *oy,int n);


/* ************************************************************** */
/* Basic statistics                                               */
/* ************************************************************** */

float  TKfindMin_f(float *x,int n);
float  TKfindMedian_f(float *val,int count);
double TKfindMedian_d(double *val,int count);
float  TKfindRMS_f(float *x,int n);
float  TKfindRMS_d(double *x,int n);
float TKfindRMSweight_d(double *x,double *e,int n);
float  TKfindMax_f(float *x,int n);
float  TKmean_f(float *x,int n);
double TKmean_d(double *x,int n);
double TKvariance_d(double *x,int n);
double TKrange_d(double *x,int n);
float TKrange_f(float *x,int n);
double TKfindMin_d(double *x,int n);
double TKfindMax_d(double *x,int n);

double TKmean_d(double *x,int n);
double TKfindMin_d(double *x,int n);
double TKsign_d(double a,double b);
double TKretMax_d(double a,double b);
double TKretMin_d(double a,double b);
float TKretMax_f(float a,float b);
float TKretMin_f(float a,float b);
int TKretMin_i(int a,int b);

/* ************************************************************** */
/* Sorting                                                        */
/* ************************************************************** */

void TKsort_f(float *val,int nobs);
void TKsort_d(double *val,int nobs);
void TKsort_2f(float *val,float *val2,int nobs);
void TKsort_3d(double *val,double *val2,double *val3,int nobs);

/* ************************************************************** */
/* Routines that modify a data array                              */
/* ************************************************************** */

void TKzeromean_d(int n,double *y);

/* ************************************************************** */
/* Random number routines                                         */
/* ************************************************************** */

double TKranDev(long *seed);
double TKgaussDev(long *seed);
long TKsetSeed();
void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double genrand_real1(void);

