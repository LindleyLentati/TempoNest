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
		int &doTimeMargin,
		int &doJumpMargin,
		int &FitSig,
		int &customPriors,
		double *EFACPrior,
		double *EQUADPrior,
		double *AlphaPrior,
		double *AmpPrior,
		int &numCoeff,
		double *CoeffPrior){

//General parameters:
//Root of the results files,relative to the directory in which TempoNest is run. This will be followed by the pulsar name, and then the individual output file extensions.
strcpy( root, "results/vHRedLinear-");

//numTempo2its - sets the number of iterations Tempo2 should do before setting the priors.  Should only be set to 0 if all the priors are set in setTNPriors
numTempo2its=10;

//doLinearFit:  Switches between the full non linear timing model (doLinearFit=0) and the linear approximation for the timing model based on the initial Tempo2 Fit (= 1).

doLinearFit=1;

//doMax: Find maximum likelihood values for non linear timing model for the stochastic model chosen.
//Will find maximum in the full un marginalised problem in order to find the best values to then marginalise over.
//Start point of non linear search will be performed at this maximum if chosen unless set otherwise in custom priors.
//Central point of prior for non linear search will be set at this value unless set otherwise in custom priors.
doMax=0;

//Model Choice

incEFAC=0; //include EFAC: 0 = none, 1 = one for all residuals, 2 = one for each observing system
incEQUAD=0; //include EQUAD: 0 = no, 1 = yes
incRED=1; //include Red Noise model: 0 = no, 1 = power law model (vHL2013), 2 = model independant (L2013)

doTimeMargin=1; //0=No Analytical Marginalisation over Timing Model. 1=Marginalise over QSD. 2=Marginalise over all Model params excluding jumps.
doJumpMargin=1; //0=No Analytical Marginalisation over Jumps. 1=Marginalise over Jumps.



//Priors

//Which priors to use: customPriors=0 uses the Priors from tempo2 fit, along with values set in this function, =1:set priors for specific parameters in setTNPriors
customPriors=1; 


//FitSig sets the priors for all timing model and jump parameters for both non linear and linear timing models.
//For the non linear fit, Fitsig multiples the error returned by Tempo2, and sets the prior to be the best fit value returned by tempo2 +/- the scaled error.
// For the linear fit, multiplies the ratio of the rms of the designmatrix vector for each timing model parameter, and the rms of the residuals returned by Tempo2.
FitSig=40;

//Remaining priors for the stochastic parameters.  
EFACPrior[0]=0.1;
EFACPrior[1]=10;


EQUADPrior[0]=-10;
EQUADPrior[1]=-0;


AlphaPrior[0]=1.1;
AlphaPrior[1]=6.1;


AmpPrior[0]=-12.8;
AmpPrior[1]=-14.0;

numCoeff=3;

CoeffPrior[0]=-20;
CoeffPrior[1]=0;

}

void setTNPriors(double **Dpriors, long double **TempoPriors){

//This function overwrites the default values for the priors sent to multinest, and the long double priors used by tempo2, you need to be aware of what dimension is what if you use this function.

//THe order of the parameters is always the same:
//Timing Model parameters (linear or non linear)
//Jumps
//EFAC(s)
//EQUAD
//Red Noise Parameters (Amplitude then Alpha for incRed=1, coefficients 1..n for incRed=2)

TempoPriors[0][0]=4.51091490228764;
//TempoPriors[0][1]=50.0;
TempoPriors[1][0]=0.136026598480345;
//TempoPriors[1][1]=50.0;
// TempoPriors[2][0]=0.0;
// TempoPriors[2][1]=50.0;
// TempoPriors[3][0]=0.0;
// TempoPriors[3][1]=50.0;
TempoPriors[4][0]=15.993628344346177451;
//TempoPriors[4][1]=0.00042845996245686584;
TempoPriors[5][0]=4.9161262510066621634;
//TempoPriors[5][1]=0.00337318189604535103;
TempoPriors[6][0]=-3.9208688084104856911;
//TempoPriors[6][1]=0.00723446217934695449;
TempoPriors[7][0]=0.93590451813572095969;
//TempoPriors[7][1]=0.04157845195654089748;
TempoPriors[8][0]=67.825129683919098071;
//TempoPriors[8][1]=0.00001317880165501870;
TempoPriors[9][0]=54303.635387746974406;
//TempoPriors[9][1]=0.00065029991091978233;
TempoPriors[10][0]=32.342422339041552736;
//TempoPriors[10][1]=0.00000030541335920403;
TempoPriors[11][0]=176.2041567162837821;
//TempoPriors[11][1]=0.00000000134802713855;
TempoPriors[12][0]=7.4940259711370419016e-05;
//TempoPriors[12][1]=0.00000000134802713855;
TempoPriors[13][0]=3.0016648691769097383e-13;
//TempoPriors[13][1]=6.6088722303607495049e-13;
TempoPriors[14][0]=-3.6932521274153480851e-05;
//TempoPriors[14][1]=0.00037669664480641853;
TempoPriors[15][0]=0.31122975057229583705;
//TempoPriors[15][1]=0.02917726591964311250;
TempoPriors[16][0]=93.905813630806996208;
//TempoPriors[16][1]=3.46689182108917082203;
TempoPriors[17][0]=71.139153996866426141;
//TempoPriors[17][1]=1.41137258470735016402;

 Dpriors[0][0]=0;
 Dpriors[0][1]=0;
 Dpriors[1][0]=-10;
 Dpriors[1][1]=10;
 Dpriors[2][0]=-15;
 Dpriors[2][1]=10;
 Dpriors[3][0]=0;
 Dpriors[3][1]=0;
 Dpriors[4][0]=0;
 Dpriors[4][1]=0;
 Dpriors[5][0]=-8;
 Dpriors[5][1]=8;
 Dpriors[6][0]=-20;
 Dpriors[6][1]=20;
 Dpriors[7][0]=-35;
 Dpriors[7][1]=15;
 Dpriors[8][0]=-10;
 Dpriors[8][1]=10;
 Dpriors[9][0]=-8;
 Dpriors[9][1]=8;
 Dpriors[10][0]=-10;
 Dpriors[10][1]=5;
 Dpriors[11][0]=-8;
 Dpriors[11][1]=10;
 Dpriors[12][0]=-10;
 Dpriors[12][1]=5;
 Dpriors[13][0]=-8;
 Dpriors[13][1]=8;
 Dpriors[14][0]=-8;
 Dpriors[14][1]=8;
 Dpriors[15][0]=-8;
 Dpriors[15][1]=8;
 Dpriors[16][0]=-15;
 Dpriors[16][1]=10;
 Dpriors[17][0]=-5;
 Dpriors[17][1]=5;
 Dpriors[18][0]=-8;
 Dpriors[18][1]=8;
 Dpriors[19][0]=0;
 Dpriors[19][1]=0;
}
