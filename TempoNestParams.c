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
		int &numCoeff,
		double *CoeffPrior){

//General parameters:
//Root of the results files,relative to the directory in which TempoNest is run. This will be followed by the pulsar name, and then the individual output file extensions.
strcpy( root, "results/Example1-NL-");

//numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.  Should only be set to 0 if all the priors are set in setTNPriors
numTempo2its=10;

//doLinearFit:  Switches between the full non linear timing model (doLinearFit=0) and the linear approximation for the timing model based on the initial Tempo2 Fit (= 1).

doLinearFit=0;

//doMax: Find maximum likelihood values for non linear timing model for the stochastic model chosen.
//Will find maximum in the full un marginalised problem in order to find the best values to then marginalise over.
//Start point of non linear search will be performed at this maximum if chosen unless set otherwise in custom priors.
//Central point of prior for non linear search will be set at this value unless set otherwise in custom priors.
doMax=0;

//Model Choice

incEFAC=0; //include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each observing system
incEQUAD=0; //include EQUAD: 0 = no, 1 = yes
incRED=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independant (L2013)
incDM=0; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independant (L2013)

doTimeMargin=0; //0=No Analytical Marginalisation over Timing Model. 1=Marginalise over QSD. 2=Marginalise over all Model params excluding jumps.
doJumpMargin=0; //0=No Analytical Marginalisation over Jumps. 1=Marginalise over Jumps.



//Priors

//Which priors to use: customPriors=0 uses the Priors from tempo2 fit, along with values set in this function, =1:set priors for specific parameters in setTNPriors
customPriors=0; 


//FitSig sets the priors for all timing model and jump parameters for both non linear and linear timing models.
//For the non linear fit, Fitsig multiples the error returned by Tempo2, and sets the prior to be the best fit value returned by tempo2 +/- the scaled error.
// For the linear fit, multiplies the ratio of the rms of the designmatrix vector for each timing model parameter, and the rms of the residuals returned by Tempo2.
FitSig=5;

//Remaining priors for the stochastic parameters.  
EFACPrior[0]=0.0;
EFACPrior[1]=5;


EQUADPrior[0]=-10;
EQUADPrior[1]=-5;



AlphaPrior[0]=1.1;
AlphaPrior[1]=6.1;


AmpPrior[0]=-20;
AmpPrior[1]=-10;

numCoeff=0;

CoeffPrior[0]=-10;
CoeffPrior[1]=0;

DMAlphaPrior[0]=1.1;
DMAlphaPrior[1]=6.1;


DMAmpPrior[0]=-18;
DMAmpPrior[1]=-8;


}

void setTNPriors(double **Dpriors, long double **TempoPriors){

//This function overwrites the default values for the priors sent to multinest, and the long double priors used by tempo2, you need to be aware of what dimension is what if you use this function.

//THe order of the parameters is always the same:
//Timing Model parameters (linear or non linear)
//Jumps
//EFAC(s)
//EQUAD
//Red Noise Parameters (Amplitude then Alpha for incRed=1, coefficients 1..n for incRed=2)



}
