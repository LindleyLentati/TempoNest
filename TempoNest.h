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
	double Tspan;
	int TimetoMargin;
	int totRedShapeCoeff;
	int totalsize;
	
	long double **LDpriors;
	double **Dpriors;
	double **DMatrix;
	double **FMatrix;
	double **GMatrix;
	double **staticGMatrix;  //staticGMatrix and staicTimedet for speeding up when no fitted white noise
	double staticTimeDet;
	double **UMatrix;     //UMatrix and SVec are for speeding up when only have 1EFAC and 1EQUAD
	double *SVec; 
	double *maxLikeRes;
	int numberpulsars;
	int doLinear;
	int numFitJumps;
	int numFitTiming;
	int numFitEFAC;
	int EPolTerms;
	int numFitEQUAD;
	int incDMEQUAD;
	int *includeEQsys;
	int numFitRedCoeff;
	int numFitDMCoeff;
	int totCoeff;
	int numFitRedPL;
	int numFitDMPL;
	double *sampleFreq;
	int numdims;
	int incRED;
	int incDM;
	int incFloatDM;
	int incFloatRed;
	int yearlyDM;
	int incsinusoid;
	int FloatDMstart;
    	int FloatRedstart;
	int Gsize;
	int Dsize;
	int **TempoFitNums;
	int *TempoJumpNums;
	int *sysFlags;
	double *TobsInfo;
	int systemcount;
	int TimeMargin;
	int JumpMargin;
	std::string *name;
	int incStep;
	char *whiteflag;
	int whitemodel;
	int varyRedCoeff;
	int varyDMCoeff;
	int incGWB;
	int incDMEvent;
	int incDMShapeEvent;
	int numDMShapeCoeff; 
	int incRedShapeEvent;
	int MarginRedShapeCoeff;
	int numRedShapeCoeff;
	int incDMScatterShapeEvent;
	int numDMScatterShapeCoeff;
	int RedPriorType;
	int DMPriorType;
	int EQUADPriorType;
	int EFACPriorType;
	int useOriginalErrors;
	int incShannonJitter;
	int incNGJitter;
	int numNGJitterEpochs;
	double **NGJitterMatrix;
	int *NGJitterSysFlags;
	int incGlitch;
	int incGlitchTerms;
	int incBreakingIndex;
	int FitLowFreqCutoff;
	int uselongdouble;
	int incGroupNoise;
	int numFitGroupNoiseCoeff;
	int **FitForGroup;
	int numGroupstoFit;
	double **GroupStartTimes;
	int incBandNoise;
	int numFitBandNoiseCoeff;
	int **FitForBand;
	int printResiduals;
	int *GroupNoiseFlags;
	int FitSolarWind;
	int FitWhiteSolarWind;
	int storeFMatrices;
	double *StoredTMatrix;
	double *StoredDMVals;
	int *TimingGradientSigns;
	/*GPTA stuff*/

	int incWideBandNoise;
	int *numshapecoeff;
	int totshapecoeff;
	int *numshapestoccoeff;
	int totalshapestoccoeff;
	int TOAnumber;
	int FixProfile;
	int InterpolateProfile;
	int NumToInterpolate;
	double InterpolatedTime;
	//double ***InterpolatedShapelets;
	double **InterpolatedShapeletsVec;
	double **InterpolatedJitterProfileVec;
	double **InterpolatedMeanProfile;
	double **InterpolatedJitterProfile;
	double **InterpolatedWidthProfile;
	int *numChanPerInt;
	int numProfileEpochs;
	int TotalProfiles;
	double *MeanProfileShape;
	double **MeanProfileEvo;
	double MeanProfileBeta;
	long double **ProfileInfo;
	long double ***ProfileData;
	long double ReferencePeriod;
	double *Factorials;
	double *Binomial;
	double MaxShapeAmp;	
	double offPulseLevel;
	int incHighFreqStoc;
	int incProfileEvo;
	int FitEvoExponent;
	double EvoRefFreq;
	int *numEvoCoeff;
	int totalEvoCoeff;
	int *numEvoFitCoeff;
	int totalEvoFitCoeff;
	int incProfileFit;
	int *numProfileFitCoeff;
	int totalProfileFitCoeff;
	int totalCoeffForMult;
	int FitLinearProfileWidth;
	int numProfComponents;
	int incWidthJitter;
	int JitterProfComp;
	int incProfileEnergyEvo;
	int ProfileBaselineTerms;
	int incProfileNoise;
	int ProfileNoiseCoeff;
	int SubIntToFit;
	int ChannelToFit;
	double **MLProfileNoise;
	int NProfileEvoPoly;
	double *ProfCompSeps;
	int diagonalGHS;
	double PhasePrior;
	int WriteNewML;

	/*Template Stuff*/

	int numTempFreqs;
	double *TemplateFreqs;
	double *TemplateChans;

	/*Grade Stuff: Need to track Grades and store previous likelihood things for hierarchial evaluation*/

	int sampler;
	int doGrades;
	int *PolyChordGrades;
	double *PriorsArray;
	int PreviousInfo;
	double PreviousJointDet;
	double PreviousFreqDet;
	double PreviousUniformPrior;
        double *LastParams;
	double **PreviousTNT;
	double **PreviousNT;
	double *PreviousNoise;
	int *hypercube_indices;

	int debug;

} MNStruct;




void assigncontext(void *context);
void assignGPUcontext(void *context);
void assignGHScontext(void *context);

void callGHS();

double iter_factorial(unsigned int n);
void store_factorial();
void fastephemeris_routines(pulsar *psr,int npsr);
void fastSubIntephemeris_routines(pulsar *psr,int npsr);
void fastformBatsAll(pulsar *psr,int npsr);
void fastformSubIntBatsAll(pulsar *psr,int npsr);
void outputProfile(int ndim);

void TNtextOutput(pulsar *psr, int npsr, int newpar, long double *Tempo2Fit, void *context, int incRED, int ndims, std::vector<double> paramlist, double Evidence, int MarginTime, int MarginJumps, int doLinear, std::string longname, double **paramarray);
void getmaxlikeDM(pulsar *pulse,std::string longname, int ndim, void *context, double **paramsarray);


double  AllTOALike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double  AllTOAJitterLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double  AllTOASim(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double AllTOAStocProfLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double AllTOAMaxLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double  AllTOAWriteMaxLike(std::string longname, int &ndim);
double AllTOAMarginStocProfLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double TemplateProfLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
void  WriteMaxTemplateProf(std::string longname, int &ndim);
void  WriteSubIntStocProfLike(std::string longname, int &ndim);
void PreComputeShapelets(double **StoredShapelet, double **StoredJitter, double **InterpolatedMeanProfile, double **InterpolatedJitterProfile, double **InterpolatedWidthProfile, long double finalInterpTime, int numtointerpolate, double MeanBeta, double &MaxShapeAmp);
void Tscrunch(void *globalcontext);
void getNumTempFreqs(int &NumFreqs, void *context);

double ProfileDomainLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
void ProfileDomainLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void  WriteProfileDomainLike(std::string longname, int &ndim);


double  Template2DProfLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
void Template2DProfLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void WriteTemplate2DProfLike(std::string longname, int &ndim);
void  WriteDMSignal(std::string longname, int &ndim);
void TemplateProfLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew, void *context);
double SubIntStocProfLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
void SubIntStocProfLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//non linear timing model likelihood functions
//double  WhiteLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
double NewLRedMarginLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);
//void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//double LRedNumericalLogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context);

void LRedLikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void LRedGPULikeMNWrap(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//non-Gaussian noise likelihoods
void subtractMLsolution(void *context);
void processPDF(void *context, std::string longname);
void TempoNestNGConvolvedLikeFunc(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void TempoNestNGLikeFunc(double *Cube, int &ndim, int &npars, double &lnew, void *context);


//GPU non linear timing model likelihood functions
//void WhiteMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//void LRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
double NewLRedMarginGPULogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *globalcontext);
void NelderMeadOptimum(int nParameters);



//Functions to calculate the design matrices or 'G' marginalisation matrices
//void makeGDesign(pulsar *pulse, int &Gsize, int numtofit, double** staticGMatrix, double **oneDesign);
//void getDMatrix(pulsar *pulse, int TimeToFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM);
//void getMarginDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM, int linearFit);
void getCustomDMatrix(pulsar *pulse, int *MarginList, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int incDM, int TimetoFit, int JumptoFit);
//void makeStaticGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** staticGMatrix, double &tdet);
//void makeStaticDiagGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** UMatrix, double *SVec);
void getCustomDMatrixLike(void *context, double **TNDM);
void getCustomDVectorLike(void *context, double *TNDM, int nobs, int TimeToMargin, int TotalSize);
void getNGJitterMatrix(pulsar *pulse, double **JitterMatrix, int &NumEpochs);
void getNGJitterMatrixEpochs(pulsar *pulse, int &NumEpochs);
void StoreTMatrix(double *TMatrix, void *context);
void getArraySizeInfo(void *context);
void getPhysDVector(void *context, double **TNDM, int Nobs, int *TimingGradientSigns);
void UpdatePhysDVector(void *context, double **TNDM, int Nobs);

void readsummary(pulsar *psr, std::string longname, int ndim, void *context, long double *Tempo2Fit, int incRED, int ndims, int MarginTime, int MarginJumps, int doLinear);

void setupMNparams(int &sampler, int &IS, int &modal, int &ceff, int &nlive, double &efr, int &sample, int &updInt, int &nClsPar, int &Nchords);
void setupparams(int &useGPUS,
		char *root,
		int &numTempo2its,
		int &doLinearFit, 
		int &doMax,
		int &incEFAC,
		int &EPolyTerms,
		int &incEQUAD,
		int &incRED,
		int &incDM,
		int &doTimeMargin,
		int &doJumpMargin,
		double &FitSig,
		int &customPriors,
		double *EFACPrior,
		double *EPolyPriors,
		double *EQUADPrior,
		double *AlphaPrior,
		double *AmpPrior,
		double *DMAlphaPrior,
		double *DMAmpPrior,
		double &numRedCoeff,
		double &numDMCoeff,
		int &numRedPL,
		int &numDMPL,
		double *RedCoeffPrior,
		double *DMCoeffPrior,
		int &FloatingDM,
		double *DMFreqPrior,
		int &yearlyDM,
		int &incsinusoid,
		int &FloatingRed,
		double *RedFreqPrior,
		double &FourierSig,
		int &incStep,
		double *StepAmpPrior,
		char *whiteflag,
		int &whitemodel,
		int &varyRedCoeff,
		int &varyDMCoeff,
		int &incGWB,
		double *GWBAmpPrior,
		int &RedPriorType,
		int &DMPriorType,
		int &EQUADPriorType,
		int &EFACPriorType,
		int &useOriginalErrors,
		int &incShannonJitter,
		int &incDMEvent,
		double *DMEventStartPrior,
		double *DMEventLengthPrior,
		int &incDMShapeEvent,
		int &numDMShapeCoeff,
		double *DMShapeCoeffPrior,
		int &incRedShapeEvent,
		int &numRedShapeCoeff,
		int &MarginRedShapeCoeff,
		double *RedShapeCoeffPrior,
		int &incDMScatterShapeEvent,
		int &numDMScatterShapeCoeff,
		double *DMScatterShapeCoeffPrior,
		int &incBandNoise,
		int &numBandNoiseCoeff,
		double *BandNoiseAmpPrior,
		double *BandNoiseAlphaPrior,
		int &incNGJitter,
		int &incGlitch,
                int &incGlitchTerms,
                double &GlitchFitSig,
		int &incBreakingIndex,
		int &FitLowFreqCutoff,
		int &uselongdouble,
		int &incGroupNoise,
		int &numGroupCoeff,
		double *GroupNoiseAmpPrior,
		double *GroupNoiseAlphaPrior,
		int &FitSolarWind,
		int &FitWhiteSolarWind,
		double *SolarWindPrior,
		double *WhiteSolarWindPrior,
		int &GPTA,
		char *GroupNoiseFlag,
		int &FixProfile,
		int &FitTemplate,
		int &InterpolateProfile,
		double &InterpolatedTime,
		int &StoreTMatrix, 
		int &incHighFreqStoc,
		double *HighFreqStocPrior,
		int &incProfileEvo,
		double &EvoRefFreq, 
		double *ProfileEvoPrior,
		int &FitEvoExponent,
		int &incWideBandNoise,
		int &incProfileFit,
		double *ProfileFitPrior,
		int &FitLinearProfileWidth,
		double *LinearProfileWidthPrior,
		int &incDMEQUAD,
		double *DMEQUADPrior,
		double &offPulseLevel,
		char *ProfFile,
		int &numProfComponents,
		int &incWidthJitter,
		double *WidthJitterPrior,
		int &JitterProfCompint,
		int &incProfileEnergyEvo,
		double *ProfileEnergyEvoPrior,
		int &debug,
		int &ProfileBaselineTerms,
		int &incProfileNoise,
		int &ProfileNoiseCoeff,
		double *ProfileNoiseAmpPrior,
		double *ProfileNoiseSpecPrior,
		int &SubIntToFit,
		int &ChannelToFit,
		int &NProfileEvoPoly);

void setTNPriors(double **Dpriors, long double **TempoPriors, int TPsize, int DPsize);
void setFrequencies(double *SampleFreq, int numRedfreqs, int numDMfreqs, int numRedLogFreqs, int numDMLogFreqs, double RedLowFreq, double DMLowFreq, double RedMidFreq, double DMMidFreq);
void GetGroupsToFit(int incGroupNoise, int **FitForGroup, int incBandNoise, int **FitForBand);
void setShapePriors(double **ShapePriors, double *BetaPrior, int numcoeff);
void GetProfileFitInfo(int numProfComponents, int *numGPTAshapecoeff, int *numProfileFitCoeff, int *numEvoCoeff, int *numFitEvoCoeff, 	int *numGPTAstocshapecoeff, double *ProfCompSeps);
