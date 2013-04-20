//  Copyright (C) 2013 Lindley Lentati

/*
*    This file is part of TempoNest 
* 
*    TempoNest is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TempoNest  is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TempoNest.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TempoNest and as a byproduct both Tempo2 and MultiNest
*    then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for TempoNest,
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model and MultiNest Papers here.
*/


#include <vector>

typedef struct {
	pulsar *pulse;
	long double **LDpriors;
	double **Dpriors;
	double **DMatrix;
	double **FMatrix;
	double **GMatrix;
	int numberpulsars;
	int numFitJumps;
	int numFitTiming;
	int numFitEFAC;
	int numFitEQUAD;
	int numFitRedCoeff;
	int numdims;
	int incRED;
	int Gsize;
	int Dsize;
	int **TempoFitNums;
	int *TempoJumpNums;
	int *sysFlags;
	int TimeMargin;
	int JumpMargin;
	
} MNStruct;



double iter_factorial(unsigned int n);
void store_factorial();
void fastephemeris_routines(pulsar *psr,int npsr);
void fastformBatsAll(pulsar *psr,int npsr);
void TNtextOutput(pulsar *psr, int npsr, int newpar, long double *Tempo2Fit, void *context, int incRED, int ndims, std::vector<double> paramlist, double Evidence, int MarginTime, int MarginJumps, int doLinear);

void NelderMeadOptimum(int nParameters, double *pdParameters, void *context);

void doSim(int argc, char **commandLine, pulsar *psr, char timFile[][MAX_FILELEN], char parFile[][MAX_FILELEN]);
void TNSimRedfromTim(int argc, char **commandLine, pulsar *psr, char timFile[][MAX_FILELEN], char parFile[][MAX_FILELEN], double EFAC, double EQUAD, int doRed, double redlogamp, double redslope, int updateEFAC, int updateEQUAD, int doDM, double DMlogamp, double DMslope,long idum);


//Linear timing model likelihood functions and utilities
void getLinearPriors(pulsar *psr,  double **Dmatrix, long double **LDpriors, double **Dpriors, int numtofit, int fitsig);
void convertFromLinear(pulsar *psr, std::string longname, int ndim, void *context);
void TNupdateParameters(pulsar *psr,int p,double *val,double *error, double *outval);
void TNIupdateParameters(pulsar *psr,int p,double *val,double *error, double *outval);

void WhiteLinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void WhiteMarginLinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedLinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedMarginLinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedLinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedMarginLinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);


//non linear timing model likelihood functions
void WhiteLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void WhiteMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

void LRedDMMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//GPU linear timing model likelihood functions
void WhiteMarginGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedMarginGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedMarginGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

//GPU non linear timing model likelihood functions
void WhiteMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void vHRedMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

void LRedDMMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);


//Functions to calculate the design matrices or 'G' marginalisation matrices
void makeGDesign(pulsar *pulse, int &Gsize, int numtofit, double** staticGMatrix, double **oneDesign);
void getDMatrix(pulsar *pulse, int TimeToFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM);
void getMarginDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM, int linearFit);



void readsummary(pulsar *psr, std::string longname, int ndim, void *context, long double *Tempo2Fit, int incRED, int ndims, int MarginTime, int MarginJumps, int doLinear);

void setupMNparams(int &IS, int &modal, int &ceff, int &nlive, double &efr);
void setupparams(char *Type,
		int &numTempo2its, 
		int &doLinearFit, 
		int &doMax,
		int &incEFAC,
		int &incEQUAD,
		int &incRED,
		int &doTimeMargin,
		int &doJumpMargin,
		int &FitSig,
		int &diffPriors,
		double *EFACPrior,
		double *EQUADPrior,
		double *AlphaPrior,
		double *AmpPrior,
		int &numCoeff,
		double *CoeffPrior);


void setTNPriors(double **Dpriors, long double **TempoPriors);
