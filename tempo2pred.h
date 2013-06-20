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

//-*-C++-*-
#ifndef TEMPO2PRED_H
#define TEMPO2PRED_H

#include <stdio.h>

#ifdef __cplusplus
#include <string>
extern "C" {
#endif

  /*******************************************************************/
  /* TEMPO2 Predictive mode library interface for C/C++              */
  /* Initial version: Russell Edwards 14 Feb 2006.                   */
  /*******************************************************************/
  
  /************************* C TYPES *********************************/
  /* Note: Most of these types are for internal use only. The        */
  /*       external caller need only know about T2Predictor and      */
  /*       its associated functions declared further down            */
  /*******************************************************************/

  /* Basic Chebyshev polynomial */
  typedef struct 
  {
    int nx, ny;
    long double *coeff;
  } Cheby2D;
  
  /* Chebyshev predictive phase model */
  typedef struct
  {
    char psrname[64];
    char sitename[64];
    long double mjd_start, mjd_end;
    long double freq_start, freq_end;
    long double dispersion_constant; // phase = polyco + d_c/f^2
    Cheby2D cheby;
    Cheby2D frequency_cheby;
  } ChebyModel;

  /* Set of Chebyshev predictive phase models */
  typedef struct
  {
    ChebyModel *segments;
    int nsegments;
  } ChebyModelSet;

  /* Flag set by ChebyModelSet_GetNearest when mjd is out of range */
  extern int ChebyModelSet_OutOfRange;

  /* TEMPO1-style Taylor series polynomial predictive phase model */
  typedef struct
  {
    char psrname[64];
    char date_string[10];
    char utc_string[13];
    long double mjd_mid;
    double dm;
    double doppler;
    double log10rms;
    long double reference_phase;
    long double frequency_psr_0;
    char sitename[5];
    int span;
    int ncoeff;
    double frequency_obs;
    double binary_phase;
    double binary_frequency;
    long double coeff[32];
  } T1Polyco;
  
  /* Set of TEMPO1-style polycos */
  typedef struct
  {
    T1Polyco *segments;
    int nsegments;
  } T1PolycoSet;

  /* Types of models available */
  typedef enum  {None, Cheby, T1} T2PredictorKind;

  /* Generic set-of-predictive models type */
  typedef struct
  {
    T2PredictorKind kind;
    union
    {
      ChebyModelSet cheby;
      T1PolycoSet t1;
    } modelset;
  } T2Predictor;

  /************************* C FUNCTIONS ******************************/

  /* Initialization, Copying, Building, etc */
  void T2Predictor_Init(T2Predictor *t2p);
  void T2Predictor_Copy(T2Predictor *into_t2p, const T2Predictor *from_t2p);
  int  T2Predictor_Insert(T2Predictor *into_t2p, const T2Predictor *from_t2p);
  void T2Predictor_Destroy(T2Predictor *t2p);

  /* I/O etc */
  int T2Predictor_Read(T2Predictor *t2p, char *fname);
  int T2Predictor_FRead(T2Predictor *t2p, FILE *f);
  void T2Predictor_Write(const T2Predictor *t2p, char *fname);
  void T2Predictor_FWrite(const T2Predictor *t2p, FILE *f);

  /* Information */
  char * T2Predictor_GetPSRName(T2Predictor *t2p);
  char * T2Predictor_GetSiteName(T2Predictor *t2p);
  long double T2Predictor_GetStartMJD(T2Predictor *t2p);
  long double T2Predictor_GetEndMJD(T2Predictor *t2p);
  long double T2Predictor_GetStartFreq(T2Predictor *t2p); // MHz
  long double T2Predictor_GetEndFreq(T2Predictor *t2p);  // MHz
  T2PredictorKind T2Predictor_Kind(T2Predictor *t2p);

  /* Prediction */
  long double T2Predictor_GetPhase(const T2Predictor *t2p, long double mjd,
				   long double freq); // freq in MHz
  long double T2Predictor_GetFrequency(const T2Predictor *t2p, long double mjd,
				       long double freq); // freq in MHz
  /* This function makes a piecewise approximation. At the moment it just
     interpolates between phase evaluations. In future it will minimise
     the mean offset. Returns 0 on success */
  int T2Predictor_GetPlan(char *filename,
			   long double mjd_start,
			   long double mjd_end,
			   long double step, // seconds
			   long double freq, // MHz
			   // remaining arguments are returned
			   long double *phase0,
			   int *nsegments,
			   long double *pulse_frequencies);
   int T2Predictor_GetPlan_Ext(char *filename,
			   long double mjd_start,
			   long double mjd_end,
			   long double step, // seconds
			   long double freq, // MHz
			   // remaining arguments are returned
			       char *psrname, char *sitename,
			   long double *phase0,
			   int *nsegments,
			   long double *pulse_frequencies);
 

#ifdef __cplusplus
}
#endif

#endif
