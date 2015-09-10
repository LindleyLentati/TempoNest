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


#include <string.h>
#include <stdio.h>
#include <math.h>
#include "configfile.h"

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
		int &numGPTAshapecoeff,
                int &numGPTAstocshapecoeff,
		char *GroupNoiseFlag,
		int &FixProfile,
		int &FitTemplate,
		int &InterpolateProfile,
		double &InterpolatedTime,
		int &StoreFMatrices,
		int &incHighFreqStoc,
		double *HighFreqStocPrior){

	//General parameters:
	//Use GPUs 0=No, 1=Yes
		
	useGPUS=0;
	uselongdouble=0;
	GPTA = 0;
	numGPTAshapecoeff=0;
	numGPTAstocshapecoeff=0;
	FixProfile = 0;	
	FitTemplate = 0;
	InterpolateProfile = 0;
	InterpolatedTime = 1; //in nanoseconds
	StoreFMatrices = 0; // Recompute FMatrices when computing new bats - default is dont just precompute and use those


    //Root of the results files,relative to the directory in which TempoNest is run. This will be followed by the pulsar name, and then the individual output file extensions.
    strcpy( root, "results/Example1-");

    //numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.  Should only be set to 0 if all the priors are set in setTNPriors
    numTempo2its=1;

    //doLinearFit:  Switches between the full non linear timing model (doLinearFit=0) and the linear approximation for the timing model based on the initial Tempo2 Fit (= 1).

    doLinearFit=0;

    //doMax: Find maximum likelihood values for non linear timing model for the stochastic model chosen.
    //Will find maximum in the full un marginalised problem in order to find the best values to then marginalise over.
    //Start point of non linear search will be performed at this maximum if chosen unless set otherwise in custom priors.
    //Central point of prior for non linear search will be set at this value unless set otherwise in custom priors.
    doMax=0;

    //ModelChoices

    //White noise

    useOriginalErrors = 0;  // Use tempo2 errors before modification by TNEF/TNEQ
    incEFAC=0; //include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each observing system
    EPolyTerms = 1; //Number of terms to include in EFAC polynomial (A*TOAerr + B*TOAERR^2 etc)
    incEQUAD=0; //include EQUAD: 0 = no, 1 = yes
    incNGJitter = 0; //Include NG Jitter for systems flagged in par file

	FitSolarWind = 0; // Basically for for ne_sw
	FitWhiteSolarWind = 0; //Fit for a white component proportional to tdis2 with ne_sw=1


    strcpy( whiteflag, "-sys");
    whitemodel=0;
    
    incShannonJitter=0;

    incRED=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013)
    incDM=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013)

	FitLowFreqCutoff = 0; //Include f_low as a free parameter

    doTimeMargin=0 ; //0=No Analytical Marginalisation over Timing Model. 1=Marginalise over QSD. 2=Marginalise over all Model params excluding jumps.
    doJumpMargin=0; //0=No Analytical Marginalisation over Jumps. 1=Marginalise over Jumps.
    incBreakingIndex=0; //Use breaking index to constrain F1/F2 etc



    //Priors

    //Which priors to use: customPriors=0 uses the Priors from tempo2 fit, along with values set in this function, =1:set priors for specific parameters in setTNPriors
    customPriors=0; 


	//FitSig sets the priors for all timing model and jump parameters for both non linear and linear timing models.
	//For the non linear fit, Fitsig multiples the error returned by Tempo2, and sets the prior to be the best fit value returned by tempo2 +/- the scaled error.
	// For the linear fit, multiplies the ratio of the rms of the designmatrix vector for each timing model parameter, and the rms of the residuals returned by Tempo2.
	FitSig=5;
	
	//Remaining priors for the stochastic parameters.  

	RedPriorType = 0; // 0 = Log, 1 = Uniform
	DMPriorType = 0;   // 0 = Log, 1 = Uniform
	EQUADPriorType = 0;   // 0 = Log, 1 = Uniform
	EFACPriorType = 0;   // 0 = Log, 1 = Uniform



	EFACPrior[0]=0.1;
	EFACPrior[1]=10;

	EPolyPriors[0]=-20;
	EPolyPriors[1]=20;   //Prior on terms in EPoly > linear
	
	
	EQUADPrior[0]=-10;
	EQUADPrior[1]=-5;


	SolarWindPrior[0] = 0;
	SolarWindPrior[1] = 10;
	WhiteSolarWindPrior[0] = -3;
	WhiteSolarWindPrior[1] = 3;
	
	numRedPL=0;
	numDMPL=0;

	numRedCoeff=10;
	numDMCoeff=10;

	varyRedCoeff=0;
	varyDMCoeff=0;
	
	AlphaPrior[0]=1.1;
	AlphaPrior[1]=6.1;
	
	AmpPrior[0]=-20;
	AmpPrior[1]=-10;
	
	DMAlphaPrior[0]=1.1;
	DMAlphaPrior[1]=6.1;
	
	DMAmpPrior[0]=-18;
	DMAmpPrior[1]=-8;

	
	RedCoeffPrior[0]=-10;
	RedCoeffPrior[1]=0;
	
	DMCoeffPrior[0]=-10;
	DMCoeffPrior[1]=0;
	
	FourierSig = 5;


	
	FloatingDM = 0;
	DMFreqPrior[0]=1;
	DMFreqPrior[1]=100;
	
	FloatingRed = 0;
	RedFreqPrior[0]=1;
	RedFreqPrior[1]=100;

	yearlyDM=0;
	incsinusoid=0;
	
	incStep = 0;
	StepAmpPrior[0] = -1;
	StepAmpPrior[1] = 1;

	incGWB=0;
	GWBAmpPrior[0] = -20;
	GWBAmpPrior[1] = -10;


	incDMEvent = 0;

	DMEventStartPrior[0] = 50000;
	DMEventStartPrior[1] = 51000;

	DMEventLengthPrior[0] = 14;
	DMEventLengthPrior[1] = 365.25;


	incDMShapeEvent = 0;
	numDMShapeCoeff = 0;
	DMShapeCoeffPrior[0] = -0.01;
	DMShapeCoeffPrior[1] = 0.01;

	incRedShapeEvent = 0;
	numRedShapeCoeff = 0;
	MarginRedShapeCoeff = 0;
	RedShapeCoeffPrior[0] = -0.01;
	RedShapeCoeffPrior[1] = 0.01;

	incDMScatterShapeEvent = 0;
	numDMScatterShapeCoeff = 0;
	DMScatterShapeCoeffPrior[0] = -0.01;
	DMScatterShapeCoeffPrior[1] = 0.01;


	incBandNoise = 0;
	numBandNoiseCoeff = 10;
	BandNoiseAmpPrior[0] = -20;
	BandNoiseAmpPrior[1] = -8;
	BandNoiseAlphaPrior[0] = 0;
	BandNoiseAlphaPrior[1] = 7;


	strcpy( GroupNoiseFlag, "-sys");
	incGroupNoise = 0;
	numGroupCoeff = 10;
	GroupNoiseAmpPrior[0] = -20;
	GroupNoiseAmpPrior[1] = -8;
	GroupNoiseAlphaPrior[0] = 0;
	GroupNoiseAlphaPrior[1] = 7;



	//GPTA Params
	
	incHighFreqStoc = 0;
	HighFreqStocPrior[0] = -10;
	HighFreqStocPrior[1] = 1;



    // Use a configfile, if we can, to overwrite the defaults set in this file.
    try {
        string strBuf;
        strBuf = string("defaultparameters.conf");
        ConfigFile parameters(strBuf);

        /* We can check whether a value is not set in the file by doing
         * if(! parameters.readInto(variable, "name", default)) {
         *   printf("WARNING");
         * }
         *
         * At the moment I was too lazy to print warning messages, and the
         * default value from this file is used in that case.
         *
         * Note: the timing model parameters are not done implemented yet
         */
	parameters.readInto(useGPUS, "useGPUS", useGPUS);
	parameters.readInto(uselongdouble, "uselongdouble", uselongdouble);
	parameters.readInto(GPTA, "GPTA", GPTA);
	parameters.readInto(numGPTAshapecoeff, "numGPTAshapecoeff", numGPTAshapecoeff);
	parameters.readInto(numGPTAstocshapecoeff, "numGPTAstocshapecoeff", numGPTAstocshapecoeff);
	parameters.readInto(FixProfile, "FixProfile", FixProfile);
	parameters.readInto(FitTemplate, "FitTemplate", FitTemplate);
	parameters.readInto(InterpolateProfile, "InterpolateProfile", InterpolateProfile);
	parameters.readInto(InterpolatedTime, "InterpolatedTime", InterpolatedTime);
        parameters.readInto(StoreFMatrices, "StoreFMatrices",StoreFMatrices );


        parameters.readInto(strBuf, "root", string("results/Example1"));
        strcpy(root, strBuf.data());
        parameters.readInto(numTempo2its, "numTempo2its", numTempo2its);
        parameters.readInto(doLinearFit, "doLinearFit", doLinearFit);
        parameters.readInto(doMax, "doMax", doMax);
        parameters.readInto(incEFAC, "incEFAC", incEFAC);
        parameters.readInto(EPolyTerms, "EPolyTerms", EPolyTerms);
        parameters.readInto(incEQUAD, "incEQUAD", incEQUAD);


	parameters.readInto(FitSolarWind, "FitSolarWind", FitSolarWind);
        parameters.readInto(FitWhiteSolarWind, "FitWhiteSolarWind", FitWhiteSolarWind);
        parameters.readInto(SolarWindPrior[0], "SolarWindPrior[0]", SolarWindPrior[0]);
        parameters.readInto(SolarWindPrior[1], "SolarWindPrior[1]", SolarWindPrior[1]);
        parameters.readInto(WhiteSolarWindPrior[0], "WhiteSolarWindPrior[0]", WhiteSolarWindPrior[0]);
        parameters.readInto(WhiteSolarWindPrior[1], "WhiteSolarWindPrior[1]", WhiteSolarWindPrior[1]);




        parameters.readInto(incRED, "incRED", incRED);
		parameters.readInto(incDM, "incDM", incDM);
        parameters.readInto(doTimeMargin, "doTimeMargin", doTimeMargin);
        parameters.readInto(doJumpMargin, "doJumpMargin", doJumpMargin);
        parameters.readInto(customPriors, "customPriors", customPriors);
        parameters.readInto(FitSig, "FitSig", FitSig);
        parameters.readInto(EFACPrior[0], "EFACPrior[0]", EFACPrior[0]);
        parameters.readInto(EFACPrior[1], "EFACPrior[1]", EFACPrior[1]);
	parameters.readInto(EPolyPriors[0], "EPolyPriors[0]", EPolyPriors[0]);
	parameters.readInto(EPolyPriors[1], "EPolyPriors[1]", EPolyPriors[1]);
        parameters.readInto(EQUADPrior[0], "EQUADPrior[0]", EQUADPrior[0]);
        parameters.readInto(EQUADPrior[1], "EQUADPrior[1]", EQUADPrior[1]);
        parameters.readInto(AlphaPrior[0], "AlphaPrior[0]", AlphaPrior[0]);
        parameters.readInto(AlphaPrior[1], "AlphaPrior[1]", AlphaPrior[1]);
        parameters.readInto(AmpPrior[0], "AmpPrior[0]", AmpPrior[0]);
        parameters.readInto(AmpPrior[1], "AmpPrior[1]", AmpPrior[1]);
        parameters.readInto(numRedCoeff, "numRedCoeff", numRedCoeff);
		parameters.readInto(numDMCoeff, "numDMCoeff", numDMCoeff);
        parameters.readInto(numRedPL, "numRedPL", numRedPL);
		parameters.readInto(numDMPL, "numDMPL", numDMPL);
        parameters.readInto(RedCoeffPrior[0], "RedCoeffPrior[0]", RedCoeffPrior[0]);
        parameters.readInto(RedCoeffPrior[1], "RedCoeffPrior[1]", RedCoeffPrior[1]);
        parameters.readInto(DMCoeffPrior[0], "DMCoeffPrior[0]", DMCoeffPrior[0]);
        parameters.readInto(DMCoeffPrior[1], "DMCoeffPrior[1]", DMCoeffPrior[1]);
        parameters.readInto(DMAlphaPrior[0], "DMAlphaPrior[0]", DMAlphaPrior[0]);
        parameters.readInto(DMAlphaPrior[1], "DMAlphaPrior[1]", DMAlphaPrior[1]);
        parameters.readInto(DMAmpPrior[0], "DMAmpPrior[0]", DMAmpPrior[0]);
        parameters.readInto(DMAmpPrior[1], "DMAmpPrior[1]", DMAmpPrior[1]);
        
        parameters.readInto(FloatingDM, "FloatingDM", FloatingDM);
	parameters.readInto(yearlyDM, "yearlyDM", yearlyDM);
	parameters.readInto(incsinusoid, "incsinusoid", incsinusoid);
        parameters.readInto(DMFreqPrior[0], "DMFreqPrior[0]", DMFreqPrior[0]);
        parameters.readInto(DMFreqPrior[1], "DMFreqPrior[1]", DMFreqPrior[1]);
        parameters.readInto(FloatingRed, "FloatingRed", FloatingRed);
        parameters.readInto(RedFreqPrior[0], "RedFreqPrior[0]", RedFreqPrior[0]);
        parameters.readInto(RedFreqPrior[1], "RedFreqPrior[1]", RedFreqPrior[1]);

        parameters.readInto(incStep, "incStep", incStep);
        parameters.readInto(StepAmpPrior[0], "StepAmpPrior[0]", StepAmpPrior[0]);
        parameters.readInto(StepAmpPrior[1], "StepAmpPrior[1]", StepAmpPrior[1]);

	parameters.readInto(strBuf, "whiteflag", string("-sys"));
        strcpy(whiteflag, strBuf.data());
	parameters.readInto(whitemodel, "whitemodel", whitemodel);

        parameters.readInto(FourierSig, "FourierSig", FourierSig);

        parameters.readInto(varyRedCoeff, "varyRedCoeff", varyRedCoeff);
        parameters.readInto(varyDMCoeff, "varyDMCoeff", varyDMCoeff);
	parameters.readInto(incGWB, "incGWB", incGWB);
	parameters.readInto(GWBAmpPrior[0], "GWBAmpPrior[0]", GWBAmpPrior[0]);
	parameters.readInto(GWBAmpPrior[1], "GWBAmpPrior[1]", GWBAmpPrior[1]);
	parameters.readInto(RedPriorType, "RedPriorType", RedPriorType);
	parameters.readInto(DMPriorType, "DMPriorType", DMPriorType);
	parameters.readInto(EFACPriorType, "EFACPriorType", EFACPriorType);
	parameters.readInto(EQUADPriorType, "EQUADPriorType", EQUADPriorType);
	parameters.readInto(useOriginalErrors, "useOriginalErrors", useOriginalErrors);
	parameters.readInto(incShannonJitter, "incShannonJitter", incShannonJitter);
	parameters.readInto(incNGJitter, "incNGJitter", incNGJitter);
	parameters.readInto(incGlitch, "incGlitch", incGlitch);
        parameters.readInto(incGlitchTerms, "incGlitchTerms", incGlitchTerms);        
	parameters.readInto(GlitchFitSig, "GlitchFitSig", GlitchFitSig);	
	parameters.readInto(incBreakingIndex, "incBreakingIndex", incBreakingIndex);
	parameters.readInto(FitLowFreqCutoff, "FitLowFreqCutoff", FitLowFreqCutoff);


	parameters.readInto(incDMEvent, "incDMEvent", incDMEvent);
	parameters.readInto(DMEventStartPrior[0], "DMEventStartPrior[0]", DMEventStartPrior[0]);
	parameters.readInto(DMEventStartPrior[1], "DMEventStartPrior[1]", DMEventStartPrior[1]);
	parameters.readInto(DMEventLengthPrior[0], "DMEventLengthPrior[0]", DMEventLengthPrior[0]);
	parameters.readInto(DMEventLengthPrior[1], "DMEventLengthPrior[1]", DMEventLengthPrior[1]);


	parameters.readInto(incDMShapeEvent, "incDMShapeEvent", incDMShapeEvent);
	parameters.readInto(numDMShapeCoeff, "numDMShapeCoeff", numDMShapeCoeff);
	parameters.readInto(DMShapeCoeffPrior[0], "DMShapeCoeffPrior[0]", DMShapeCoeffPrior[0]);
	parameters.readInto(DMShapeCoeffPrior[1], "DMShapeCoeffPrior[1]", DMShapeCoeffPrior[1]);

	parameters.readInto(incRedShapeEvent, "incRedShapeEvent", incRedShapeEvent);
	parameters.readInto(numRedShapeCoeff, "numRedShapeCoeff", numRedShapeCoeff);
	parameters.readInto(MarginRedShapeCoeff, "MarginRedShapeCoeff", MarginRedShapeCoeff);
	parameters.readInto(RedShapeCoeffPrior[0], "RedShapeCoeffPrior[0]", RedShapeCoeffPrior[0]);
	parameters.readInto(RedShapeCoeffPrior[1], "RedShapeCoeffPrior[1]", RedShapeCoeffPrior[1]);

	parameters.readInto(incDMScatterShapeEvent, "incDMScatterShapeEvent", incDMScatterShapeEvent);
	parameters.readInto(numDMScatterShapeCoeff, "numDMScatterShapeCoeff", numDMScatterShapeCoeff);
	parameters.readInto(DMScatterShapeCoeffPrior[0], "DMScatterShapeCoeffPrior[0]", DMScatterShapeCoeffPrior[0]);
	parameters.readInto(DMScatterShapeCoeffPrior[1], "DMScatterShapeCoeffPrior[1]", DMScatterShapeCoeffPrior[1]);


	parameters.readInto(incBandNoise, "incBandNoise", incBandNoise);
	parameters.readInto(numBandNoiseCoeff, "numBandNoiseCoeff", numBandNoiseCoeff);
	parameters.readInto(BandNoiseAmpPrior[0], "BandNoiseAmpPrior[0]", BandNoiseAmpPrior[0]);
	parameters.readInto(BandNoiseAmpPrior[1], "BandNoiseAmpPrior[1]", BandNoiseAmpPrior[1]);
	parameters.readInto(BandNoiseAlphaPrior[0], "BandNoiseAlphaPrior[0]", BandNoiseAlphaPrior[0]);
	parameters.readInto(BandNoiseAlphaPrior[1], "BandNoiseAlphaPrior[1]", BandNoiseAlphaPrior[1]);


	parameters.readInto(strBuf, "GroupNoiseFlag", string("-sys"));
        strcpy(GroupNoiseFlag, strBuf.data());
	parameters.readInto(incGroupNoise, "incGroupNoise", incGroupNoise);
	parameters.readInto(numGroupCoeff, "numGroupCoeff", numGroupCoeff);
	parameters.readInto(GroupNoiseAmpPrior[0], "GroupNoiseAmpPrior[0]", GroupNoiseAmpPrior[0]);
	parameters.readInto(GroupNoiseAmpPrior[1], "GroupNoiseAmpPrior[1]", GroupNoiseAmpPrior[1]);
	parameters.readInto(GroupNoiseAlphaPrior[0], "GroupNoiseAlphaPrior[0]", GroupNoiseAlphaPrior[0]);
	parameters.readInto(GroupNoiseAlphaPrior[1], "GroupNoiseAlphaPrior[1]", GroupNoiseAlphaPrior[1]);

        parameters.readInto(incHighFreqStoc, "incHighFreqStoc",incHighFreqStoc);
        parameters.readInto(HighFreqStocPrior[0], "HighFreqStocPrior[0]", HighFreqStocPrior[0]);
        parameters.readInto(HighFreqStocPrior[1], "HighFreqStocPrior[1]", HighFreqStocPrior[1]);

	
	
    } catch(ConfigFile::file_not_found oError) {
        printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
    } // try

}

void setTNPriors(double **Dpriors, long double **TempoPriors, int TPsize, int DPsize){

//This function overwrites the default values for the priors sent to multinest, and the long double priors used by tempo2, you need to be aware of what dimension is what if you use this function.

//THe order of the parameters is always the same:
//Timing Model parameters (linear or non linear)
//Jumps
//EFAC(s)
//EQUAD
//Red Noise Parameters (Amplitude then Alpha for incRed=1, coefficients 1..n for incRed=2)

	for(int i =0;i<TPsize; i++){	
	//	printf("TP %i \n", i);

	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
  		n=sprintf (buffer, "TempoPriors[%i][0]", i);
		parameters.readInto(TempoPriors[i][0], buffer, TempoPriors[i][0]);
		n=sprintf (buffer, "TempoPriors[%i][1]", i);
                parameters.readInto(TempoPriors[i][1], buffer, TempoPriors[i][1]);
		n=sprintf (buffer, "TempoPriors[%i][2]", i);
		parameters.readInto(TempoPriors[i][2], buffer, TempoPriors[i][2]);
		n=sprintf (buffer, "TempoPriors[%i][3]", i);
                parameters.readInto(TempoPriors[i][3], buffer, TempoPriors[i][3]);


	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}



	for(int i =0;i<DPsize; i++){	


	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
                char buffer [50];
                int n;
                n=sprintf (buffer, "Dpriors[%i][0]", i);
				parameters.readInto(Dpriors[i][0], buffer, Dpriors[i][0]);
				n=sprintf (buffer, "Dpriors[%i][1]", i);
				parameters.readInto(Dpriors[i][1], buffer, Dpriors[i][1]);

	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}





}


void setFrequencies(double *SampleFreq, int numRedfreqs, int numDMfreqs, int numRedLogFreqs, int numDMLogFreqs, double RedLowFreq, double DMLowFreq, double RedMidFreq, double DMMidFreq){

//This function sets or overwrites the default values for the sampled frequencies sent to multinest


	int startpoint=0;
	double RedLogDiff = log10(RedMidFreq) - log10(RedLowFreq);
	for(int i =0; i < numRedLogFreqs; i++){
		SampleFreq[startpoint]=pow(10.0, log10(RedLowFreq) + i*RedLogDiff/numRedLogFreqs);
		startpoint++;
		printf("%i %g %g \n", i, log10(RedLowFreq) - i*log10(RedLowFreq)/numRedLogFreqs, SampleFreq[startpoint-1]);
		
	}
	
	for(int i =0;i < numRedfreqs-numRedLogFreqs; i++){	
		SampleFreq[startpoint]=i+RedMidFreq;
		startpoint++;
		printf("making freqs %i %g\n", startpoint-1, SampleFreq[startpoint-1]);
	
	}
        for(int i =0;i < numDMfreqs; i++){
                SampleFreq[startpoint]=i+1;
		startpoint++;
		//printf("making freqs %i %g", startpoint+i, SampleFreq[startpoint+i]);
        }

	startpoint=0;
	for(int i =0;i<numRedfreqs; i++){	


	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
  		n=sprintf (buffer, "SampleFreq[%i]", i);
		parameters.readInto(SampleFreq[i], buffer, SampleFreq[i]);

	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}


	startpoint=startpoint+numRedfreqs;
	for(int i =0;i<numDMfreqs; i++){	


	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
  		n=sprintf (buffer, "SampleFreq[%i]", startpoint+i);
		parameters.readInto(SampleFreq[startpoint+i], buffer, SampleFreq[startpoint+i]);

	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}
}



void GetGroupsToFit(int incGroupNoise, int **FitForGroup, int incBandNoise, int **FitForBand){

//This function reads in the groups that will be fit as Group Noise terms



	for(int i =0;i<incGroupNoise; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			/* We can check whether a value is not set in the file by doing
			* if(! parameters.readInto(variable, "name", default)) {
			*   printf("WARNING");
			* }
			*
			* At the moment I was too lazy to print warning messages, and the
			* default value from this file is used in that case.
			*
			* Note: the timing model parameters are not done implemented yet
			*/
			char buffer [50];
			int n;
			n=sprintf (buffer, "FitForGroup[%i][0]", i);
			parameters.readInto(FitForGroup[i][0], buffer, FitForGroup[i][0]);
                        n=sprintf (buffer, "FitForGroup[%i][1]", i);
                        parameters.readInto(FitForGroup[i][1], buffer, FitForGroup[i][1]);
                        n=sprintf (buffer, "FitForGroup[%i][2]", i);
                        parameters.readInto(FitForGroup[i][2], buffer, FitForGroup[i][2]);
                        n=sprintf (buffer, "FitForGroup[%i][3]", i);
                        parameters.readInto(FitForGroup[i][3], buffer, FitForGroup[i][3]);
                        n=sprintf (buffer, "FitForGroup[%i][4]", i);
                        parameters.readInto(FitForGroup[i][4], buffer, FitForGroup[i][4]);
                        n=sprintf (buffer, "FitForGroup[%i][5]", i);
                        parameters.readInto(FitForGroup[i][5], buffer, FitForGroup[i][5]);


		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}


	for(int i =0;i<incBandNoise; i++){	


		// Use a configfile, if we can, to overwrite the defaults set in this file.
		try {
			string strBuf;
			strBuf = string("defaultparameters.conf");
			ConfigFile parameters(strBuf);

			/* We can check whether a value is not set in the file by doing
			* if(! parameters.readInto(variable, "name", default)) {
			*   printf("WARNING");
			* }
			*
			* At the moment I was too lazy to print warning messages, and the
			* default value from this file is used in that case.
			*
			* Note: the timing model parameters are not done implemented yet
			*/
			char buffer [50];
			int n;
			n=sprintf (buffer, "FitForBand[%i][0]", i);
			parameters.readInto(FitForBand[i][0], buffer, FitForBand[i][0]);
                        n=sprintf (buffer, "FitForBand[%i][1]", i);
                        parameters.readInto(FitForBand[i][1], buffer, FitForBand[i][1]);
                        n=sprintf (buffer, "FitForBand[%i][2]", i);
                        parameters.readInto(FitForBand[i][2], buffer, FitForBand[i][2]);
                        n=sprintf (buffer, "FitForBand[%i][3]", i);
                        parameters.readInto(FitForBand[i][3], buffer, FitForBand[i][3]);

		} 
		catch(ConfigFile::file_not_found oError) {
			printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}
}




void setShapePriors(double **ShapePriors, double *BetaPrior, int numcoeff){


	for(int i =0;i<numcoeff; i++){	
	//	printf("TP %i \n", i);

	    // Use a configfile, if we can, to overwrite the defaults set in this file.
	    try {
		string strBuf;
		strBuf = string("defaultparameters.conf");
		ConfigFile parameters(strBuf);

		/* We can check whether a value is not set in the file by doing
		 * if(! parameters.readInto(variable, "name", default)) {
		 *   printf("WARNING");
		 * }
		 *
		 * At the moment I was too lazy to print warning messages, and the
		 * default value from this file is used in that case.
		 *
		 * Note: the timing model parameters are not done implemented yet
		 */
		char buffer [50];
		int n;
		
  		n=sprintf (buffer, "ShapePriors[%i][0]", i);
		parameters.readInto(ShapePriors[i][0], buffer, ShapePriors[i][0]);
		n=sprintf (buffer, "ShapePriors[%i][1]", i);
                parameters.readInto(ShapePriors[i][1], buffer, ShapePriors[i][1]);


	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

	}



    // Use a configfile, if we can, to overwrite the defaults set in this file.
    try {
	string strBuf;
	strBuf = string("defaultparameters.conf");
	ConfigFile parameters(strBuf);

	/* We can check whether a value is not set in the file by doing
	 * if(! parameters.readInto(variable, "name", default)) {
	 *   printf("WARNING");
	 * }
	 *
	 * At the moment I was too lazy to print warning messages, and the
	 * default value from this file is used in that case.
	 *
	 * Note: the timing model parameters are not done implemented yet
	 */
        parameters.readInto(BetaPrior[0], "BetaPrior[0]", BetaPrior[0]);
        parameters.readInto(BetaPrior[1], "BetaPrior[1]", BetaPrior[1]);


    } catch(ConfigFile::file_not_found oError) {
	printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
    } // try


}
