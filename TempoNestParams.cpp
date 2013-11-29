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
#include "configfile.h"

void setupparams(char *root,
		int &numTempo2its,
		int &doLinearFit, 
		int &doMax,
		int &incEFAC,
		int &incEQUAD,
		int &incRED,
		int &incDM,
		int &doTimeMargin,
		int &doJumpMargin,
		double &FitSig,
		int &customPriors,
		double *EFACPrior,
		double *EQUADPrior,
		double *AlphaPrior,
		double *AmpPrior,
		double *DMAlphaPrior,
		double *DMAmpPrior,
		int &numRedCoeff,
		int &numDMCoeff,
		int &numRedPL,
		int &numDMPL,
		double *RedCoeffPrior,
		double *DMCoeffPrior,
		int &FloatingDM,
		double *DMFreqPrior,
		int &FloatingRed,
		double *RedFreqPrior,
		double &FourierSig,
		int &incStep,
		double *StepAmpPrior,
		char *whiteflag,
		int &whitemodel){

    //General parameters:
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


    incEFAC=0; //include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each observing system
    incEQUAD=0; //include EQUAD: 0 = no, 1 = yes
    strcpy( whiteflag, "-sys");
    whitemodel=0;

    incRED=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013)
    incDM=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independent (L2013)

    doTimeMargin=0 ; //0=No Analytical Marginalisation over Timing Model. 1=Marginalise over QSD. 2=Marginalise over all Model params excluding jumps.
    doJumpMargin=0; //0=No Analytical Marginalisation over Jumps. 1=Marginalise over Jumps.



    //Priors

    //Which priors to use: customPriors=0 uses the Priors from tempo2 fit, along with values set in this function, =1:set priors for specific parameters in setTNPriors
    customPriors=0; 


	//FitSig sets the priors for all timing model and jump parameters for both non linear and linear timing models.
	//For the non linear fit, Fitsig multiples the error returned by Tempo2, and sets the prior to be the best fit value returned by tempo2 +/- the scaled error.
	// For the linear fit, multiplies the ratio of the rms of the designmatrix vector for each timing model parameter, and the rms of the residuals returned by Tempo2.
	FitSig=5;
	
	//Remaining priors for the stochastic parameters.  
	EFACPrior[0]=0.1;
	EFACPrior[1]=10;
	
	
	EQUADPrior[0]=-10;
	EQUADPrior[1]=-5;
	
	numRedPL=0;
	numDMPL=0;

	numRedCoeff=10;
	numDMCoeff=10;

// 	varyRedCoeff=0;
// 	varyDMCoeff=0;
	
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
	
	incStep = 0;
	StepAmpPrior[0] = -1;
	StepAmpPrior[1] = 1;


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

        parameters.readInto(strBuf, "root", string("results/Example1"));
        strcpy(root, strBuf.data());
        parameters.readInto(numTempo2its, "numTempo2its", numTempo2its);
        parameters.readInto(doLinearFit, "doLinearFit", doLinearFit);
        parameters.readInto(doMax, "doMax", doMax);
        parameters.readInto(incEFAC, "incEFAC", incEFAC);
        parameters.readInto(incEQUAD, "incEQUAD", incEQUAD);
        parameters.readInto(incRED, "incRED", incRED);
		parameters.readInto(incDM, "incDM", incDM);
        parameters.readInto(doTimeMargin, "doTimeMargin", doTimeMargin);
        parameters.readInto(doJumpMargin, "doJumpMargin", doJumpMargin);
        parameters.readInto(customPriors, "customPriors", customPriors);
        parameters.readInto(FitSig, "FitSig", FitSig);
        parameters.readInto(EFACPrior[0], "EFACPrior[0]", EFACPrior[0]);
        parameters.readInto(EFACPrior[1], "EFACPrior[1]", EFACPrior[1]);
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

//         parameters.readInto(varyRedCoeff, "varyRedCoeff", varyRedCoeff);
//         parameters.readInto(varyDMCoeff, "varyDMCoeff", varyDMCoeff);


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


void setFrequencies(double *SampleFreq, int numRedfreqs, int numDMfreqs){

//This function sets or overwrites the default values for the sampled frequencies sent to multinest


	int startpoint=0;
	for(int i =0;i < numRedfreqs; i++){	
		SampleFreq[startpoint+i]=i+1;
		//printf("making freqs %i %g", startpoint+i, SampleFreq[startpoint+i]);
	
	}
	startpoint=startpoint+numRedfreqs;
        for(int i =0;i < numDMfreqs; i++){
                SampleFreq[startpoint+i]=i+1;
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
