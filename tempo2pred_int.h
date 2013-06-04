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

/* Tempo2 prediction mode -- internal routines, i.e. those that don't
 * need to be available (i.e. declared) to callers from external programs;
   Note, some of these will need *linkage* into the tempo2pred library. */

#include "tempo2.h"
#include "tempo2pred.h"

void ChebyModel_Construct(ChebyModel *cm, const pulsar *psr);
void ChebyModel_Test(ChebyModel *cm, const pulsar *psr, int nmjd, int nfreq,
		     long double *residualRMS, long double *residualMAV);
void ChebyModelSet_Construct(ChebyModelSet *cms, const pulsar *psr,
			     const char *sitename,
			     long double mjd_start, long double mjd_end,
			     long double segment_length, long double overlap,
			     long double freq_start,long double freq_end,
			     int nmjdcoeff, int nfreqcoeff);
void ChebyModelSet_Test(ChebyModelSet *cms, const pulsar *psr,
			int nmjd, int nfreq,
			long double *residualRMS, long double *residualMAV);

#ifdef __cplusplus
extern "C" {
#endif
void Cheby2D_Construct(Cheby2D *cheby,
		  void
		  (*func)(long double *x, long double *y, 
			  int nx, int ny, long double *z, void *info),
		  void *info);
void Cheby2D_Construct_x_Derivative(Cheby2D *dcheby, const Cheby2D *cheby);
void Cheby2D_Test(Cheby2D *cheby, int nx_test, int ny_test,
	     void
	     (*func)(long double *x, long double *y, 
		     int nx, int ny, long double *z, void *info),
	     void *info,
		  long double *residualRMS, long double *residualMAV);

void ChebyModel_Init(ChebyModel *cmodel, int nmjdcoeff, int nfreqcoeff);
void ChebyModel_Copy(ChebyModel *cm, ChebyModel *from);
void ChebyModel_Destroy(ChebyModel *cm);
long double ChebyModel_GetPhase(const ChebyModel *cm, long double mjd, long double freq);
long double ChebyModel_GetFrequency(const ChebyModel *cm, long double mjd, long double freq);
void ChebyModel_Write(const ChebyModel *cm, FILE *f);
int ChebyModel_Read(ChebyModel *cm, FILE *f);
ChebyModel *ChebyModelSet_GetNearest(const ChebyModelSet *cms, long double mjd);
long double ChebyModelSet_GetPhase(const ChebyModelSet *cms, long double mjd, long double freq);
long double ChebyModelSet_GetFrequency(const ChebyModelSet *cms, long double mjd, long double freq);
void ChebyModelSet_Write(const ChebyModelSet *cms, FILE *f);
int ChebyModelSet_Read(ChebyModelSet *cms, FILE *f);

void ChebyModelSet_Init(ChebyModelSet *cms);
int ChebyModelSet_Insert(ChebyModelSet *cms, const ChebyModelSet *from);
void ChebyModelSet_Destroy(ChebyModelSet *cms);



long double T1Polyco_GetPhase(const T1Polyco *t1p, long double mjd, long double freq);
long double T1Polyco_GetFrequency(const T1Polyco *t1p, long double mjd, long double freq);
void T1Polyco_Write(const T1Polyco *t1p, FILE *f);
int T1Polyco_Read(T1Polyco *t1p, FILE *f);
T1Polyco *T1PolycoSet_GetNearest(long double mjd);
long double T1PolycoSet_GetPhase(const T1PolycoSet *t1ps, long double mjd, long double freq);
long double T1PolycoSet_GetFrequency(const T1PolycoSet *t1ps, long double mjd, long double freq);
void T1PolycoSet_Write(const T1PolycoSet *t1ps, FILE *f);
int T1PolycoSet_Read(T1PolycoSet *t1ps, FILE *f);
void T1PolycoSet_Destroy(T1PolycoSet *t1ps);

#ifdef __cplusplus
}
#endif
