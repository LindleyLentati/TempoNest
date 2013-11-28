#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
// Copyright (C) 2013 Lindley Lentati

/*
* This file is part of TempoNest
*
* TempoNest is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* TempoNest is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with TempoNest. If not, see <http://www.gnu.org/licenses/>.
*/

/*
* If you use TempoNest and as a byproduct both Tempo2 and MultiNest
* then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for TempoNest,
* Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2,
* pp. 655-672 (bibtex: 2006MNRAS.369..655H)
* or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
* pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
* timing model and MultiNest Papers here.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "tempo2.h"
#include "tempo2Util.h"
#include "tempo2pred.h"
#include "tempo2pred_int.h"
#include "T2accel.h"
#include <dlfcn.h>
// #include "T2toolkit.h"
#include <vector>
#include <algorithm>
#include <float.h>
#include "multinest.h"
#include "TempoNest.h"

#include <string.h>
#include <sstream>
#include <iterator>
#include <cstring>
#include <iostream>
#include <fstream>

#include <gsl/gsl_sf_gamma.h>
#include "dpotri.h"
#include "dpotrf.h"

#ifdef HAVE_CULA
#include <cula_lapack_device.h>
#endif /* HAVE_CULA */

void ephemeris_routines(pulsar *psr,int npsr);
void clock_corrections(pulsar *psr,int npsr);
void extra_delays(pulsar *psr,int npsr);

#ifdef HAVE_CULA
extern "C" void copy_gmat_(double *G, int N);
extern "C" void copy_floatgmat_(float *G, int N);
extern "C" void copy_staticgmat_(double *G, int M, int N);
extern "C" void copy_staticumat_(double *G, int M, int N);
#endif /* HAVE_CULA */


void fastephemeris_routines(pulsar *psr,int npsr)
{
vectorPulsar(psr,npsr); /* 1. Form a vector pointing at the pulsar */
//readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
//tt2tb(psr,npsr); /* Observatory/time-dependent part of TT-TB */
//readEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */

}

void fastformBatsAll(pulsar *psr,int npsr)
{
    //clock_corrections(psr,npsr); /* Clock corrections ... */
    fastephemeris_routines(psr,npsr); /* Ephemeris routines ... */
  	extra_delays(psr,npsr); /* Other time delays ... */
	formBats(psr,npsr); /* Form Barycentric arrival times */
	secularMotion(psr,npsr);

}

MNStruct* init_struct(pulsar *pulseval,	 long double **LDpriorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval, int systemcountval, int numFitEFACval, int numFitEQUADval, int numFitRedCoeffval, int numFitDMCoeffval,int numFitRedPLval, int numFitDMPLval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval, int numdimsval, int incREDval, int incDMval, int incFloatDMval, int incFloatRedval, int DMFloatstartval, int RedFloatstartval, int TimeMarginVal, int JumpMarginVal, int doLinearVal, double *SampleFreqsVal, int incStepVal, char *whiteflagval, int whitemodelval)
{
    MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));

	MNS->pulse=pulseval;
	MNS->LDpriors=LDpriorsval;
	MNS->numberpulsars=numberpulsarsval;
	MNS->numFitJumps=numFitJumpsval;
	MNS->numFitTiming=numFitTimingval;
	MNS->systemcount=systemcountval;
	MNS->numFitEFAC=numFitEFACval;
	MNS->numFitEQUAD=numFitEQUADval;
	MNS->numFitRedCoeff=numFitRedCoeffval;
	MNS->numFitDMCoeff=numFitDMCoeffval;
	MNS->numFitRedPL=numFitRedPLval;
	MNS->numFitDMPL=numFitDMPLval;
	MNS->TempoFitNums=TempoFitNumsval;
	MNS->TempoJumpNums=TempoJumpNumsval;
	MNS->sysFlags=sysFlagsval;
	MNS->numdims=numdimsval;
	MNS->incRED=incREDval;
	MNS->incDM=incDMval;
	MNS->incFloatRed=incFloatRedval;
	MNS->incFloatDM=incFloatDMval;
	MNS->FloatRedstart=RedFloatstartval;
        MNS->FloatDMstart=DMFloatstartval;
	MNS->TimeMargin=TimeMarginVal;
	MNS->JumpMargin=JumpMarginVal;
	MNS->doLinear=doLinearVal;
	MNS->sampleFreq=SampleFreqsVal;
	MNS->incStep=incStepVal;
	MNS->whiteflag=whiteflagval;
	MNS->whitemodel=whitemodelval;

	return MNS;
}

void printPriors(pulsar *psr, long double **TempoPriors, double **Dpriors, int incEFAC, int incEQUAD, int incRED, int incDM, int numRedCoeff, int numDMCoeff, int numFloatRed, int numFloatDM, int fitDMModel, std::string longname, int incStep){


	std::ofstream getdistparamnames;
	std::string gdpnfname = longname+".paramnames";
	getdistparamnames.open(gdpnfname.c_str());
	
	int getdistlabel=1;

	if(TempoPriors[0][2] == 0){
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "Phase\n";
		getdistlabel++;
	}
	
	printf("\nPriors:\n");
	int paramsfitted=0;
	
	printf("Prior on Phase : %.25Lg -> %.25Lg\n",TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
	paramsfitted++;
	
	for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<psr[0].param[p].aSize;k++){
				if(psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel){
				
					printf("Prior on %s : %.25Lg -> %.25Lg\n",psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
					
					
					if(TempoPriors[paramsfitted][2] == 0){
						getdistparamnames << getdistlabel;
						getdistparamnames << " ";
						getdistparamnames <<  psr[0].param[p].shortlabel[k];
						getdistparamnames << "\n";
						getdistlabel++;
					}

					paramsfitted++;

				}
			}
	}
		
	int jumpsfitted=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1){
		
			printf("Prior on Jump %i : %.25Lg -> %.25Lg\n",jumpsfitted+1, TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
			
			if(TempoPriors[paramsfitted][2] == 0){
				getdistparamnames << getdistlabel;
				getdistparamnames << " ";
				getdistparamnames <<  "Jump";
				getdistparamnames << jumpsfitted+1;
				getdistparamnames << "\n";
			}

			paramsfitted++;
			jumpsfitted++;
		}
	} 


	if(incStep>0){
		for(int i =0; i < incStep; i++){
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "StepAmp";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "StepTime";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

			printf("Prior for Step Amp: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior for Step Time: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
		}
	}
	

	
		
	if(incEFAC>0){
		int EFACnum=1;
		for(int i =0;i<incEFAC;i++){
			printf("Prior on EFAC %i : %.5g -> %.5g\n",EFACnum, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "EFAC";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			getdistlabel++;
	
			paramsfitted++;
			EFACnum++;
		}
	}

	if(incEQUAD>0){
		int EQUADnum=1;	
		for(int i =0;i<incEQUAD;i++){
		
			printf("Prior on EQUAD %i: %.5g -> %.5g\n",EQUADnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "EQUAD";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			getdistlabel++;
				
			paramsfitted++;
			EQUADnum++;
		}
	}
			
	if(incRED==1 || incRED==3){
	
		printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;		
		
		printf("Prior on Red Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;
		
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "RedAmp";
		getdistparamnames << "\n";
		getdistlabel++;	
		
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "RedSlope";
		getdistparamnames << "\n";
		getdistlabel++;
			
	}
	else if(incRED==2){
		int Coeffnum=1;
		for(int i =0;i<numRedCoeff;i++){
			printf("Prior on Red Noise Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
				
			paramsfitted++;	
		}
	}
	else if(incRED==4){
		 for(int i =0;i < 2*numRedCoeff;i++){
			printf("Prior on fourier coefficients: %i %g %g \n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "FourierC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			
			
			paramsfitted++;	
		}

		for(int i =0; i< numRedCoeff; i++){
		 	printf("Prior on Red Noise coefficients: %i %g %g \n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		 	
		 	getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			
			
			paramsfitted++;	
	 	}
	}
	else if(incRED==5){

		for(int i =0;i < 2*numRedCoeff;i++){
				printf("Prior on fourier coefficients: %i %g %g \n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				
				getdistparamnames << getdistlabel;
				getdistparamnames << " ";
				getdistparamnames <<  "FourierC";
				getdistparamnames <<  i+1;
				getdistparamnames << "\n";
				getdistlabel++;
			
			
				paramsfitted++;	
		}
		
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedAmp";
			getdistparamnames << "\n";
			getdistlabel++;	
		
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedSlope";
			getdistparamnames << "\n";
			getdistlabel++;

			printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;		
			printf("Prior on Red Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;
	}	

        if(numFloatRed>0){
        	for(int i =0; i < numFloatRed; i++){
		 	getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RFF ";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

		 	getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RFA ";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

			printf("Prior on Floating Red Noise Frequency %i : %.5g -> %.5g\n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;		
			printf("Prior on Floating Red Noise Log Amplitude %i : %.5g -> %.5g\n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;
		    }
        }
			
	if(incDM==1 || incDM==3){
	
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "DMAmp";
		getdistparamnames << "\n";
		getdistlabel++;	
		
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "DMSlope";
		getdistparamnames << "\n";
		getdistlabel++;
		
		
		printf("Prior on DM Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;		
		printf("Prior on DM Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;
	}
	else if(incDM==2){
		int Coeffnum=1;
		for(int i =0;i<numDMCoeff;i++){
			printf("Prior on DM Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			
			
			paramsfitted++;	
		}
	}

        if(numFloatDM>0){
        	for(int i =0; i < numFloatDM; i++){
		 	getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMFF ";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

		 	getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMFA ";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

			printf("Prior on Floating DM Frequency %i : %.5g -> %.5g\n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;		
			printf("Prior on Floating DM Log Amplitude %i : %.5g -> %.5g\n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;
		    }
        }

			
	if(fitDMModel==1){
		for(int i=0;i<psr[0].dmoffsDMnum;i++){
		
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMModel";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			
			
			printf("Prior for DMModel: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
		}
	}


	
	getdistparamnames.close(); 
			
}


void update_MNPriors(MNStruct* MNS, double ** DPriorsval, long double **priorsval, int linearPriors){


	int paramsfitted=1;

	
	if(linearPriors != 2){
	//printf("lin p is 0\n");
	for (int p=0;p<MAX_PARAMS;p++) {
        	for (int k=0;k<MNS->pulse->param[p].aSize;k++){
                	if(MNS->pulse->param[p].fitFlag[k] == 1  && p != param_dmmodel){
				if(p == param_ecc || p ==param_px || p == param_m2 || p==param_dm && k==0){
					long double minprior=priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1];
					if(minprior < 0 && priorsval[paramsfitted][2]==0){
		//				printf("%.10Lg %.10Lg %g %g \n",priorsval[paramsfitted][0],priorsval[paramsfitted][1],DPriorsval[paramsfitted+linearPriors][0],DPriorsval[paramsfitted+linearPriors][1]);
						long double newprior=-priorsval[paramsfitted][0]/priorsval[paramsfitted][1];
						DPriorsval[paramsfitted+linearPriors][0]=(double)newprior;
	                        		printf("Prior on %s updated to be physical (was <0) : %.25Lg -> %.25Lg\n", MNS->pulse->param[p].shortlabel[k], priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1],priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][1]*priorsval[paramsfitted][1]);
	          //              		printf("%.10Lg %.10Lg %g %g \n",priorsval[paramsfitted][0],priorsval[paramsfitted][1],DPriorsval[paramsfitted+linearPriors][0],DPriorsval[paramsfitted+linearPriors][1]);
					}
				}
                                paramsfitted++;
			 }
                }
       }

	}
	

	MNS->Dpriors=DPriorsval;

	MNS->LDpriors=priorsval;

}




	void update_MNGdata(MNStruct* MNS, int Gsizeval, double **GMatrixval)
	{

		MNS->Gsize=Gsizeval;
		MNS->GMatrix=GMatrixval;


	}

	void update_MNDdata(MNStruct* MNS, int Dsizeval, double **DMatrixval)
	{

  
		MNS->Dsize=Dsizeval;
		MNS->DMatrix=DMatrixval;

	}

	void update_MNGDdata(MNStruct* MNS, int Dsizeval, double **DMatrixval,int Gsizeval, double **GMatrixval)
	{


		MNS->Dsize=Dsizeval;
		MNS->DMatrix=DMatrixval;
		MNS->Gsize=Gsizeval;
		MNS->GMatrix=GMatrixval;



	}

	void update_MNstaticG(MNStruct* MNS, double **staticGMatrixVal, double staticdetval)
	{

		MNS->staticGMatrix=staticGMatrixVal;
		MNS->staticTimeDet=staticdetval;
	}
	
	void update_MNstaticDiagG(MNStruct* MNS, double **UMatrixVal, double *SVecval)
	{

		MNS->UMatrix=UMatrixVal;
		MNS->SVec=SVecval;
	}



/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

/* The main function of a plugin called from Tempo2 is 'graphicalInterface'
*/
extern "C" int graphicalInterface(int argc, char **argv,
    pulsar *psr, int *pnpsr) {
	int iteration; 
	int listparms;
	int outRes=0;
	int writeModel=0;
	char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
	char outputSO[MAX_FILELEN];
	char str[MAX_FILELEN];
	char newparname[MAX_FILELEN];
	longdouble coeff[MAX_COEFF]; /* For polynomial coefficients in polyco */
	int npsr=*pnpsr;      /* The number of pulsars */
	int noWarnings=1;
	double globalParameter=0.0;
	int  displayParams,p;
	int nGlobal,i,flagPolyco=0,it,k;
	char polyco_args[128];
	char polyco_file[128];
	int newpar=0;
	int onlypre=0;
	//  char tempo2MachineType[MAX_FILELEN]="";
	FILE *alias;
	char **commandLine;
	clock_t startClock,endClock;
	const char *CVS_verNum = "$Revision: 1.28 $";
	int numFitJumps;
	int numToMargin=0;


  printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome to redistribute it\n");
  printf("under conditions of GPL license.\n\n");

  //
  startClock = clock();


  commandLine = (char **)malloc(1000*sizeof(char *));

  for (i=0;i<1000;i++)
    commandLine[i] = (char *)malloc(sizeof(char)*1000);

  /* Parse input line for machine type */
  for (i=0;i<argc;i++)    
    {
      strcpy(commandLine[i],argv[i]);
    }

  strcpy(outputSO, "");
  npsr = 0;   /* Initialise the number of pulsars */
  displayParams=0;
  nGlobal=0;
  /* Obtain command line arguments */
  logdbg("Running getInputs %d",psr[0].nits);
  logdbg("Completed getInputs");
  getInputs(psr,argc, commandLine, timFile,parFile,&listparms,&npsr,&nGlobal,&outRes,&writeModel,
            outputSO,&flagPolyco,polyco_args,polyco_file,&newpar,&onlypre,dcmFile,covarFuncFile,newparname);

  int doMaxLike=0;

  logdbg("Reading par file");
  readParfile(psr,parFile,timFile,npsr); /* Read .par file to define the pulsar's initial parameters */  
  logdbg("Finished reading par file %d",psr[0].nits);
  if (flagPolyco==0)
    {
      logdbg("Running readTimfile");
      readTimfile(psr,timFile,npsr); /* Read .tim file to define the site-arrival-times */
      logdbg("Completed readTimfile %d",psr[0].param[param_ecc].paramSet[1]);
    }

  logdbg("Running preProcess %d",psr[0].nits);
  preProcess(psr,npsr,argc,commandLine);
  logdbg("Completed preProcess %d",psr[0].nits);
  
   if (debugFlag==1)
    {
      logdbg("Number of iterations = %d",psr[0].nits);
      logdbg("Maximum number of parameters = %d",MAX_PARAMS);
      logdbg("Number of pulsars = %d",npsr);
    }

	
	char root[100]; 
	int numTempo2its;
	int doLinearFit;
	int doMax;
	int incEFAC;
	int incEQUAD;
	int incRED;
	int incDM;
	int incFloatRed=0;
	int incFloatDM=0;
	int doTimeMargin;
	int doJumpMargin;
	double FitSig;
	int customPriors;
	int Reddims=0;
	int DMdims=0;
	int DMModeldims=0;
	int fitDMModel=0;
	double *EFACPrior;
	double *EQUADPrior;
	double *AlphaPrior;
	double *AmpPrior;
	double *DMAlphaPrior;
	double *DMAmpPrior;
	double *DMFreqPrior;
	double *RedFreqPrior;
	int numRedCoeff;
	int numDMCoeff;
	int numRedPL;
	int numDMPL;
	double *RedCoeffPrior;
	double *DMCoeffPrior;
	double FourierSig;
	double *SampleFreq;
	int numEFAC=0;
	int numEQUAD=0;
	int numStep=0;
	double *StepAmpPrior;
	double *StepTimePrior;
	char wflag[100];
	int whitemodel;

	char *Type = new char[100];
	char *WhiteName = new char[100];
	EFACPrior=new double[2];
	EQUADPrior=new double[2];
	AlphaPrior=new double[2];
	AmpPrior=new double[2];
	DMAlphaPrior=new double[2];
	DMAmpPrior=new double[2];
	RedCoeffPrior=new double[2];
	DMCoeffPrior=new double[2];
	RedFreqPrior=new double[2];
	DMFreqPrior=new double[2];
	StepAmpPrior=new double[2];
	StepTimePrior=new double[2];

	setupparams(Type, numTempo2its, doLinearFit, doMax, incEFAC, incEQUAD, incRED, incDM, doTimeMargin, doJumpMargin, FitSig, customPriors, EFACPrior, EQUADPrior, AlphaPrior, AmpPrior, DMAlphaPrior, DMAmpPrior, numRedCoeff, numDMCoeff, numRedPL, numDMPL, RedCoeffPrior, DMCoeffPrior, incFloatDM, DMFreqPrior, incFloatRed, RedFreqPrior, FourierSig, numStep, StepAmpPrior, WhiteName,whitemodel); 
	
	
	if(incRED < 2)numRedCoeff=0;
	if(incDM < 2)numDMCoeff=0;	
	SampleFreq=new double[numRedCoeff+numDMCoeff];
	setFrequencies(SampleFreq,numRedCoeff, numDMCoeff);


	
  formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
  logdbg("calling formResiduals");
  formResiduals(psr,npsr,1);       /* Form residuals */

  for (it=0;it<numTempo2its;it++) /* Why pulsar 0 should select the iterations? */
    {
      if (it>0) /* Copy post-fit values to pre-fit values */
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (p=0;p<npsr;p++)
		{
		  for (k=0;k<psr[p].param[i].aSize;k++)
		    {
		      psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
		      psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
		    }
		}
	    }
	}
      //      long seed = TKsetSeed();
      for (iteration=0;iteration<2;iteration++) /* Do pre- and post- fit analysis */
	{
	  logdbg("iteration %d",iteration);
	  logdbg("calling formBatsAll");
	  //	  printf("Calling formBats\n");
	  formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
	  logdbg("calling formResiduals");
	  formResiduals(psr,npsr,1);       /* Form residuals */


	  if (listparms==1 && iteration==0)displayParameters(13,timFile,parFile,psr,npsr); /* List out all the parameters */  
	  if (iteration==0)          /* Only fit to pre-fit residuals */
	    {
	      logdbg("calling doFit");

	      if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
		doFit(psr,npsr,writeModel); /* Fit to the residuals to obtain updated parameters */
	      else
		doFitDCM(psr,dcmFile,covarFuncFile,npsr,writeModel);
	      printf("Complete return\n");
	      /* doFitGlobal(psr,npsr,&globalParameter,nGlobal,writeModel);*/ /* Fit to the residuals to obtain updated parameters  */
	      logdbg("completed doFit");
	    }
	  if (iteration==1 || onlypre==1)
	    {
	      if (strlen(outputSO)==0){
	      //printf("CHI SQ IS: %g \n",psr->fitChisq);
		textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,newparname); /* Output results to the screen */
		//printf("CHI SQ IS: %g %g %g \n",psr->fitChisq,psr[0].offset,psr[0].offset_e);
		}
	      else  /* Use a plug in for the output */
		{
		  char *(*entry)(int,char **,pulsar *,int);
		  void * module;
		  for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
			  sprintf(str,"%s/%s_%s_plug.t2",tempo2_plug_path[iplug],
					  outputSO,tempo2MachineType);
			  printf("Looking for %s\n",str);
			  module = dlopen(str, RTLD_NOW); 
			  if(module==NULL){	  
				  printf("dlerror() = %s\n",dlerror());
			  } else break;
		  }
		  if(!module)  {
		    fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
		    return -1;
		  }
		  /*
		   * Check that the plugin is compiled against the same version of tempo2.h
		   */
		  char ** pv  = (char**)dlsym(module, "plugVersionCheck");
		  if(pv!=NULL){
			  // there is a version check for this plugin
			  if(strcmp(TEMPO2_h_VER,*pv)){
				  fprintf(stderr, "[error]: Plugin version mismatch\n");
				  fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
				  fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
				  dlclose(module);
				  return -1;
			  }
		  }


		  entry = (char*(*)(int,char **,pulsar *,int))dlsym(module, "tempoOutput");
		  if( entry == NULL ) {
		    dlclose(module);
		    fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
		    return -1;
		  }
		  entry(argc,argv,psr,npsr);
		}
	    }
	  psr[0].noWarnings=2;
	  if (onlypre==1) iteration=2;
	  /*	  textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,"new.par");*/ /* Output results to the screen */
	  /*	  printf("Next iteration\n");*/
	}
    }

	
	
	if(incRED==0)Reddims=0;
	if(incRED==1)Reddims=2;
	if(incRED==2)Reddims=numRedCoeff;
	if(incRED==3)Reddims=2*numRedPL;
	if(incRED==4)Reddims=3*numRedCoeff;
	if(incRED==5)Reddims=2*numRedCoeff+2;
	if(incDM==0)DMdims=0;
	if(incDM==1)DMdims=2;
	if(incDM==2)DMdims=numDMCoeff;
	if(incDM==3)DMdims=2*numDMPL;
	if(incDM==4)DMdims=3*numDMCoeff;
	if(incDM==5)DMdims=2*numDMCoeff+2;
	if(incFloatDM>0)DMdims+=2*incFloatDM;
    if(incFloatRed>0)Reddims+=2*incFloatRed;
    
	//printf("DMModel flag: %i \n",psr[0].param[param_dmmodel].fitFlag[0]);
	if(psr[0].param[param_dmmodel].fitFlag[0]==1){
	//	printf("Fitting for DMModel using %i structure functions \n",psr[0].dmoffsDMnum);
		DMModeldims=psr[0].dmoffsDMnum; 
		fitDMModel=1;
	}
	
	if(fitDMModel==1 && ( doTimeMargin != 0 || doJumpMargin != 0 || doLinearFit != 0)){
			printf("DMModel currently only supports doTimeMargin= 0, doJumpMargin = 0, doLinear = 0");
// 			return 0;
	}

	if((incRED == 1 && incDM != 1 && incDM != 0 )|| (incRED != 1 && incRED != 0 && incDM == 1)){
		    printf("Different methods for DM and red noise not currently supported, please use the same option for both");
		    return 0;
	}
        
	std::string pulsarname=psr[0].name;
	std::string longname=Type+pulsarname+"-";

	if(longname.size() >= 100){printf("Root Name is too long, needs to be less than 100 characters, currently %i .\n",longname.size());return 0;}

	for(int r=0;r<=longname.size();r++){root[r]=longname[r];}

	std::string wname=WhiteName;

	if(wname.size() >= 100){printf("white noise flag is too long, needs to be less than 100 characters, currently %i .\n",wname.size());return 0;}

	for(int r=0;r<=wname.size();r++){wflag[r]=wname[r];}



    printf("Graphical Interface: TempoNest\n");
    printf("Author:              L. Lentati\n");
    printf("Version:             1.0\n");
    printf("----------------------------------------------------------------\n");
    printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
    printf("This is free software, and you are welcome to redistribute it\n");
    printf("under conditions of GPL license.\n\n");
    printf("----------------------------------------------------------------\n");

	
	printf("\n\n\n\n*****************************************************\n");
	printf("Starting TempoNest\n");
	printf("*****************************************************\n\n\n\n");
	printf("Details of the fit:\n");
#ifdef HAVE_CULA
	printf("Using GPUs\n");
#endif
	printf("file root set to %s \n",root);
	
	int systemcount=0;
	int *numFlags=new int[psr[0].nobs];
	for(int o=0;o<psr[0].nobs;o++){
		numFlags[o]=0;
	}
//	if(incEFAC == 0 || incEFAC == 1){
//		if(incEFAC == 0){printf("Not Including EFAC\n");systemcount=0;}
//		if(incEFAC == 1){printf("Including One EFAC for all observations\n");systemcount=1;}
		
//		std::string whiteflagname=pulsarname+"-whiteflags.txt";
//		std::ofstream whiteflagoutput;
//		whiteflagoutput.open(whiteflagname.c_str());
//		whiteflagoutput << 1;
//		whiteflagoutput <<  "\n";
		
	
		
//		for(int o=0;o<psr[0].nobs;o++){
//			numFlags[o]=0;
//			whiteflagoutput << numFlags[o];
//			whiteflagoutput <<  "\n";
//		}
		
//		whiteflagoutput.close();
		
//	}
//	else{
		if(incEFAC > 0 || incEQUAD > 0){printf("using white noise model %i\n", whitemodel);}
// 		if(incEFAC == 2){printf("Including One EFAC for each %s\n", wflag);}
// 		if(incEQUAD == 2){printf("Including One EQUAD for each %s\n", wflag);}
		std::vector<std::string>systemnames;
		for(int o=0;o<psr[0].nobs;o++){
			int found=0;
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
				if(strcasecmp(psr[0].obsn[o].flagID[f],wflag)==0){
					
					if(std::find(systemnames.begin(), systemnames.end(), psr[0].obsn[o].flagVal[f]) != systemnames.end()) {
	// 				/* systemnames contains x */
					} else {
	// 	
	// 
		// 				/* systemnames does not contain x */
						if(incEFAC==2 || incEQUAD==2){printf("Found %s %s \n",wflag, psr[0].obsn[o].flagVal[f]);}
						systemnames.push_back(psr[0].obsn[o].flagVal[f]);
						systemcount++;
					}
					found=1;
				}

			}
			if(found==0 && (incEFAC==2 || incEQUAD==2)){printf("Observation %i is missing the %s flag, please check before continuing\n",wflag,o);return 0;}
		}
		if(incEFAC==2 || incEQUAD==2){printf("total number of systems: %i \n",systemcount);}
		if(systemcount==0 && (incEFAC < 2 && incEQUAD < 2)){systemcount=1;}
		
// 		std::string whiteflagname=pulsarname+"-whiteflags.txt";
// 		std::ofstream whiteflagoutput;
// 		whiteflagoutput.open(whiteflagname.c_str());
// 		whiteflagoutput << systemcount;
// 		whiteflagoutput <<  "\n";
		
		for(int o=0;o<psr[0].nobs;o++){
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
			
				if(strcasecmp(psr[0].obsn[o].flagID[f],wflag)==0){
					for (int l=0;l<systemcount;l++){
						if(psr[0].obsn[o].flagVal[f] == systemnames[l]){
							numFlags[o]=l;
// 							whiteflagoutput << numFlags[o];
// 							whiteflagoutput <<  "\n";
						}
					}
				}
	
			}
		}
// 		whiteflagoutput.close();
//	}


    if(incEFAC == 0){printf("Not Including EFAC\n");numEFAC=0;}
    if(incEFAC == 1){printf("Including One EFAC for all observations\n");numEFAC=1;}
    if(incEFAC == 2){printf("Including One EFAC for each %s\n", wflag);numEFAC=systemcount;}	
    if(incEQUAD == 0){printf("Not Including EQUAD\n");numEQUAD=0;}
    if(incEQUAD == 1){printf("Including One EQUAD for all observations\n");numEQUAD=1;}
    if(incEQUAD == 2){printf("Including One EQUAD for each %s\n", wflag);numEQUAD=systemcount;}
    if(incRED == 0){printf("Not Including Red Noise\n");}
    if(incRED == 1){printf("Including Red Noise : Power Law Model\n");}
    if(incRED == 2){printf("Including Red Noise : Model Independant - Fitting %i Coefficients\n", numRedCoeff);}
    if(incRED ==3){printf("Including Red Noise: Power Law Model to %i Coefficients \n", numRedCoeff);}
    if(incRED ==4){printf("Including Red Noise Numerically: Model Independant - Fitting %i Coefficients\n \n", numRedCoeff);}
    if(incRED ==5){printf("Including Red Noise Numerically: Power Law Model to %i Coefficients \n", numRedCoeff);}
    if(incDM == 0){printf("Not Including DM\n");}
    if(incDM == 1){printf("Including DM : Power Law Model\n");}
    if(incDM == 2){printf("Including DM : Model Independant - Fitting %i Coefficients\n", numDMCoeff);}
    if(incDM ==3){printf("Including DM: Power Law Model to %i Coefficients \n", numDMCoeff);}
	if(incFloatDM==1){printf("Including Floating DM power spectrum coefficient\n");}
    if(incFloatRed==1){printf("Including Floating Red noise power spectrum coefficient\n");}
	
	int fitcount=0;
	printf("fitting for: Arbitrary Phase \n");
	fitcount++;
	for (int p=0;p<MAX_PARAMS;p++) {
	      for (int k=0;k<psr[0].param[p].aSize;k++){
			if(psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel){
				printf("fitting for: %s \n",psr[0].param[p].shortlabel[k]);
				fitcount++;
	    		}
			if(psr[0].param[p].fitFlag[k] == 1 && p == param_dmmodel){
                    printf("fitting for: %s \n",psr[0].param[p].shortlabel[k]);
            }

		}
	}
// 	printf("fitting for: %i \n",fitcount);

// 	printf("num jumps: %i \n",psr[0].nJumps);

	numFitJumps=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1)numFitJumps++;
// 		  printf("%i %i %g %g \n",i,psr[0].fitJump[i], psr[0].jumpVal[i],psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree));
	}  
	printf("Found %i jumps to fit \n",numFitJumps);
	printf("total Timing Model params to fit:  %i \n",numFitJumps+fitcount);

	long double **TempoPriors;
	long double *Tempo2Fit = new long double[numFitJumps+fitcount];
	int **TempoFitNums;
	int *TempoJumpNums;
	TempoPriors=new long double*[numFitJumps+fitcount];
	for(int i=0;i<numFitJumps+fitcount;i++){
		TempoPriors[i]=new long double[3];
		for(int j=0; j< 3; j++){
				TempoPriors[i][j]=0;
		}
		
	}
	TempoFitNums=new int*[fitcount];for(int i=0;i<fitcount;i++){TempoFitNums[i]=new int[2];}
	TempoJumpNums=new int[numFitJumps];
	//printf("allocated\n");
	int paramsfitted=0;

	TempoPriors[paramsfitted][0]=0;
	TempoPriors[paramsfitted][1]=psr[0].rmsPre*pow(10.0,-6);//offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree);
	TempoFitNums[paramsfitted][0]=0;
	TempoFitNums[paramsfitted][1]=0;
	if(doTimeMargin != 0 || doJumpMargin != 0)TempoPriors[paramsfitted][2]=1;
	paramsfitted++;
	//printf("p plus\n");
	for (int p=0;p<MAX_PARAMS;p++) {
	      for (int k=0;k<psr[0].param[p].aSize;k++){
			if(psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel){
				
				TempoPriors[paramsfitted][0]=psr[0].param[p].prefit[k];
				TempoPriors[paramsfitted][1]=psr[0].param[p].err[k]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
				Tempo2Fit[paramsfitted]=psr[0].param[p].val[k];

				TempoFitNums[paramsfitted][0]=p;
				TempoFitNums[paramsfitted][1]=k;	

				if(strcasecmp(psr[0].param[p].shortlabel[0],"F0")== 0){
					if(doTimeMargin != 0)TempoPriors[paramsfitted][2]=1;
				}
				if(doTimeMargin == 2)TempoPriors[paramsfitted][2]=1;				 


				paramsfitted++;


	    		}
		}
	}
	//printf("and time\n");
	int jumpsfitted=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1){

			//printf("gonna read ump %i %s \n",i,psr[0].jumpStr[i]);	
			char str1[100],str2[100],str3[100],str4[100],str5[100];
			int nread=sscanf(psr[0].jumpStr[i],"%s %s %s %s %s",str1,str2,str3,str4,str5);
			double prejump=atof(str3);
			//printf("Pre jump %g \n",prejump);
			
			TempoPriors[paramsfitted][0]=prejump;
			TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
			Tempo2Fit[paramsfitted]=psr[0].jumpVal[i];
			TempoJumpNums[jumpsfitted]=i;
			
			if(doJumpMargin == 1)TempoPriors[paramsfitted][2]=1;
			
			paramsfitted++;
			jumpsfitted++;

		}
	}  

	// set the MultiNest sampling parameters
	
// 	return 0;
	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 500;				// number of live points
	double efr = 0.1;				// set the required efficiency

	setupMNparams(IS, mmodal, ceff, nlive, efr);


	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = numFitJumps+fitcount+numEFAC+numEQUAD+Reddims+DMdims+ numStep*2 + DMModeldims; // dimensionality (no. of free parameters)
	int nPar = ndims;					// total no. of parameters including free & derived parameters
	int nClsPar = 2;				// no. of parameters to do mode separation on
	int updInt = 500;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI set it to F if you want your main program to handle MPI initialization
	double logZero = -DBL_MAX;			// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass
	//printf("Here \n");
	
	int FloatRedstart=numFitJumps+fitcount+numEFAC+numEQUAD+Reddims-incFloatRed*2;
	int FloatDMstart=numFitJumps+fitcount+numEFAC+numEQUAD+Reddims+DMdims-incFloatDM*2;

	MNStruct *MNS = init_struct(psr,TempoPriors,npsr,numFitJumps,fitcount,systemcount,numEFAC,numEQUAD, numRedCoeff, numDMCoeff, numRedPL, numDMPL, TempoFitNums,TempoJumpNums,numFlags, ndims, incRED,incDM, incFloatDM,incFloatRed, FloatDMstart, FloatRedstart, doTimeMargin,doJumpMargin, doLinearFit, SampleFreq, numStep, wflag, whitemodel);
	
	
	context=MNS;

#ifdef HAVE_CULA

	culaStatus status;
    	status = culaInitialize();
    	store_factorial();
    
#endif /* HAVE_CULA */
 // 	printf("Di\n"); 
    double **Dpriors;
    Dpriors = new double*[ndims]; for(int i = 0; i < ndims; i++){Dpriors[i]=new double[2];};

  	
     long double *TNMaxParameters = new long double[ndims];
    //If using custompriors for errors incase T2 doesnt converge, get those values before doing anything else
    if(customPriors == 1){
		setTNPriors(Dpriors, TempoPriors, ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
		
		int pcount=1;
		for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			psr[0].param[((MNStruct *)context)->TempoFitNums[j][0]].val[((MNStruct *)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
		//	printf("TP: %i %.10Lg %.10Lg \n",pcount,TempoPriors[pcount][0],TempoPriors[pcount][1]);
			pcount++;
		}

		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
			psr[0].jumpVal[((MNStruct *)context)->TempoJumpNums[j]] =  TempoPriors[pcount][0];
	//		printf("TP: %i %.10Lg %.10Lg \n",pcount,TempoPriors[pcount][0],TempoPriors[pcount][1]);
			pcount++;
		}
	}

	if(doMaxLike==1){
#ifdef HAVE_CULA
		GPUFindMLHypervisor(ndims, context,longname);
#else
		FindMLHypervisor(ndims, context,longname);
#endif
		return 0;
	}
	
	//If wanting to find the max do so now
	if(doMax==1){
		
		NelderMeadOptimum(ndims, TNMaxParameters, context);
		for(int i =0; i< numFitJumps+fitcount; i++){
			TempoPriors[i][0]=TNMaxParameters[i];
		}
		
		int pcount=1;
		for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			psr[0].param[((MNStruct *)context)->TempoFitNums[j][0]].val[((MNStruct *)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
			pcount++;
		}

		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
			psr[0].jumpVal[((MNStruct *)context)->TempoJumpNums[j]] =  TempoPriors[pcount][0];
			pcount++;
		}
	}
	
	//reupdate any of the priors from custom priors that were overwritten by findmax
    if(customPriors == 1){
		setTNPriors(Dpriors, TempoPriors,  ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
		
		int pcount=1;
		for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			psr[0].param[((MNStruct *)context)->TempoFitNums[j][0]].val[((MNStruct *)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
//			printf("TP2: %i %.10Lg %.10Lg \n",pcount,TempoPriors[pcount][0],TempoPriors[pcount][1]);
			pcount++;
		}

		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
			psr[0].jumpVal[((MNStruct *)context)->TempoJumpNums[j]] =  TempoPriors[pcount][0];
			pcount++;
		}
	}

	if(doLinearFit != 1){

		//Combine all the priors into one aray: Dpriors
		int pcount=0;

		for(int i =0; i< numFitJumps+fitcount; i++){
			Dpriors[pcount][0]=-FitSig;
			Dpriors[pcount][1]=FitSig;
			pcount++;
		}


		if(numStep>0){
			for(int i =0; i < numStep; i++){
				Dpriors[pcount][0]=StepAmpPrior[0]*psr[0].rmsPre*pow(10.0,-6);
				Dpriors[pcount][1]=StepAmpPrior[1]*psr[0].rmsPre*pow(10.0,-6);
				pcount++;
				Dpriors[pcount][0]=psr[0].obsn[1].bat;
				Dpriors[pcount][1]=psr[0].obsn[psr[0].nobs-1].bat;
				pcount++;
			}
		}			



		for(int i =0; i< numEFAC; i++){
			Dpriors[pcount][0]=EFACPrior[0];
			Dpriors[pcount][1]=EFACPrior[1];
			pcount++;
		}
		for(int i =0; i< numEQUAD; i++){
			Dpriors[pcount][0]=EQUADPrior[0];
			Dpriors[pcount][1]=EQUADPrior[1];
			pcount++;
		}	
		if(incRED==1 || incRED==3){
			for(int i =0; i< numRedPL; i++){
				Dpriors[pcount][0]=AmpPrior[0];
				Dpriors[pcount][1]=AmpPrior[1];
				pcount++;
				Dpriors[pcount][0]=AlphaPrior[0];
				Dpriors[pcount][1]=AlphaPrior[1];
				pcount++;
			}
		}	
		else if(incRED==2){	
			for(int i =0; i< numRedCoeff; i++){
				Dpriors[pcount][0]=RedCoeffPrior[0];
				Dpriors[pcount][1]=RedCoeffPrior[1];
				pcount++;
			}
		}	
		else if(incRED==4){

            for(int i =0;i < 2*numRedCoeff;i++){
                    Dpriors[pcount][0]=-FourierSig*psr[0].rmsPre*pow(10.0,-6);
                    Dpriors[pcount][1]=FourierSig*psr[0].rmsPre*pow(10.0,-6);
 
                    pcount++;
            }
    
    		for(int i =0; i< numRedCoeff; i++){
                    Dpriors[pcount][0]=RedCoeffPrior[0];
                    Dpriors[pcount][1]=RedCoeffPrior[1];
                    pcount++;
            }
		}
        else if(incRED==5){

			for(int i =0;i < 2*numRedCoeff;i++){
				Dpriors[pcount][0]=-FourierSig*psr[p].rmsPre*pow(10.0,-6);
	            Dpriors[pcount][1]=FourierSig*psr[p].rmsPre*pow(10.0,-6);

		        pcount++;
			}

            Dpriors[pcount][0]=AmpPrior[0];
            Dpriors[pcount][1]=AmpPrior[1];
            pcount++;
            Dpriors[pcount][0]=AlphaPrior[0];
            Dpriors[pcount][1]=AlphaPrior[1];
            pcount++;
        }	
        
        if(incFloatRed>0){
        	for(int i =0; i < incFloatRed; i++){
		        Dpriors[pcount][0]=RedFreqPrior[0];
		        Dpriors[pcount][1]=RedFreqPrior[1];
		        pcount++;
		        Dpriors[pcount][0]=RedCoeffPrior[0];
		        Dpriors[pcount][1]=RedCoeffPrior[1];
		        pcount++;
		    }
		        
        }
        
        
        		
		if(incDM==1 || incDM==3){
			for(int i =0; i< numDMPL; i++){
				Dpriors[pcount][0]=DMAmpPrior[0];
				Dpriors[pcount][1]=DMAmpPrior[1];
				pcount++;
				Dpriors[pcount][0]=DMAlphaPrior[0];
				Dpriors[pcount][1]=DMAlphaPrior[1];
				pcount++;
			}
		}	
        else if(incDM==2){
                for(int i =0; i< numDMCoeff; i++){
                        Dpriors[pcount][0]=DMCoeffPrior[0];
                        Dpriors[pcount][1]=DMCoeffPrior[1];
                        pcount++;
                }
        }
        if(incFloatDM>0){
        	for(int i =0; i < incFloatDM; i++){
		        Dpriors[pcount][0]=DMFreqPrior[0];
		        Dpriors[pcount][1]=DMFreqPrior[1];
		        pcount++;
		        Dpriors[pcount][0]=DMCoeffPrior[0];
		        Dpriors[pcount][1]=DMCoeffPrior[1];
		        pcount++;
		    }
        }

	if(fitDMModel==1){
		for(int i=0;i<psr[0].dmoffsDMnum;i++){
			Dpriors[pcount][0]=psr[0].dmoffsDM[i] - FitSig*psr[0].dmoffsDM_error[i];
			Dpriors[pcount][1]=psr[0].dmoffsDM[i] + FitSig*psr[0].dmoffsDM_error[i];
			printf("DMModel: %i %g %g \n",i,psr[0].dmoffsDM[i],psr[0].dmoffsDM_error[i]);
			pcount++;
		}
	}

		//printf("set up priors, pcount: %i \n",pcount);
		
		int MarginDMQuad=0;	
		numToMargin=0;
		//if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5 && MarginDMQuad == 1)numToMargin +=2;
		for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
			if(TempoPriors[i][2]==1){
				numToMargin++;
			}
		}

		if(numToMargin>0 ){
	
			int *FitList=new int[ndims];
			for(int i=0;i<ndims;i++){
				FitList[i]=0;
			}
			numToMargin=0;
			//if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5 &&  MarginDMQuad == 1)numToMargin +=2;
			for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
				if(TempoPriors[i][2]==1){
					FitList[i]=1;
					printf("marginalising over param: %i \n",i);
					numToMargin++;
				}
			}
			int Gsize=psr[0].nobs-numToMargin;
			double **TNDM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNDM[i]=new double[numToMargin];
			}
	
			double **TNGM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNGM[i]=new double[psr[0].nobs];

			}

			getCustomDMatrix(psr, FitList, TNDM, TempoFitNums, TempoJumpNums, Dpriors, incDM, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps);
			//Get DMatrix for marginalisation, using T2 values if not custom or max, max if not custom, or else custom position.
			//getMarginDMatrix(psr, fitcount, numFitJumps, numToMargin, TempoFitNums, TempoJumpNums, Dpriors, doJumpMargin, doTimeMargin, TNDM, 0);

			makeGDesign(psr, Gsize, numToMargin, TNGM, TNDM);

			update_MNPriors(MNS,Dpriors, TempoPriors,0);
			update_MNGdata(MNS, Gsize,TNGM);

			context=MNS;

			
			//Finally after doing everything reget custom priors in case overwritten by previous steps.
			if(customPriors == 1){
				printf("Set to use custom priors, updating from setPriors function \n");
				setTNPriors(Dpriors, TempoPriors,  ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
				paramsfitted=0;
			}
	//		printf("Step size: %.20Lg \n", TempoPriors[3][1]);
			update_MNPriors(MNS,Dpriors, TempoPriors,0);
			context=MNS;
			
			
#ifdef HAVE_CULA
	double *GMatrixVec=new double[psr[0].nobs*Gsize];


	for(int g=0;g<Gsize; g++){
		for(int o=0;o<psr[0].nobs; o++){

			GMatrixVec[g*psr[0].nobs + o]=TNGM[o][g];
		}
	}

	copy_gmat_(GMatrixVec, psr[0].nobs*Gsize);


#endif /* HAVE_CULA */

			printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff, numDMCoeff, incFloatRed,incFloatDM, fitDMModel, longname, numStep);

			printf("\n\n");
			ndims=ndims-numToMargin;
			//if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5 && MarginDMQuad == 1)ndims +=2;
			nPar=ndims;

			if(numEFAC==0 && numEQUAD==0){
				printf("Not Fitting for white noise: Pre computing Matrices.\n");
				double **staticG=new double*[psr[0].nobs];
				for(int i=0;i<psr[0].nobs;i++){
					staticG[i]=new double[psr[0].nobs];
				}
				double statictdet=0;
				makeStaticGMatrix(psr, Gsize, TNGM, staticG, statictdet);
				update_MNstaticG(MNS, staticG, statictdet);
				printf("Matrices generated, static matrix determinant: %g \n",statictdet);
				context=MNS;
				
				
#ifdef HAVE_CULA
				double *staticGMatrixVec=new double[psr[0].nobs*psr[0].nobs];


				for(int g=0;g<psr[0].nobs; g++){
					for(int o=0;o<psr[0].nobs; o++){

						staticGMatrixVec[g*psr[0].nobs + o]=staticG[o][g];
					}
				}

				copy_staticgmat_(staticGMatrixVec, psr[0].nobs, psr[0].nobs);


#endif /* HAVE_CULA */


			}
			
			
			if(numEFAC==1 || numEQUAD==1 && numEFAC<2 && numEQUAD < 2){
			
				printf("Fitting for only one EFAC or EQUAD: Pre computing Matrices.\n");
				double **UM=new double*[Gsize];
				for(int i=0;i<Gsize;i++){
					UM[i]=new double[psr[0].nobs];
				}
				double *SM=new double[Gsize];
				makeStaticDiagGMatrix(psr, Gsize, TNGM, UM, SM);
				update_MNstaticDiagG(MNS, UM, SM);
				printf("Matrices generated\n");
				context=MNS;
				
#ifdef HAVE_CULA
				double *staticUMatrixVec=new double[Gsize*psr[0].nobs];


				for(int g=0;g<psr[0].nobs; g++){
					for(int o=0;o<Gsize; o++){

						staticUMatrixVec[g*Gsize + o]=UM[o][g];
					}
				}

				copy_staticumat_(staticUMatrixVec, Gsize, psr[0].nobs);


#endif /* HAVE_CULA */


			}			
			
			
			
		if(incRED==0 && incDM == 0){
#ifdef HAVE_CULA
			nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginGPULogLike, dumper, context);
#else
			nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
		else if(incRED==1 || incDM==1){
#ifdef HAVE_CULA
			nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginGPULogLike, dumper, context);
#else
			nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
		else if(incRED==2 || incDM==2 || incRED==3 || incDM==3 || incRED==4 || incRED==5){
#ifdef HAVE_CULA
			nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginGPULogLike, dumper, context);
#else
			nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
	}
	else if(doJumpMargin == 0 && doTimeMargin == 0 ){

		//If not marginalising, only need to reget custom priors to overwrite anything done previously
		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
		setTNPriors(Dpriors, TempoPriors,  ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
			}

		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
		
		
		printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff, numDMCoeff, incFloatRed,incFloatDM,fitDMModel, longname, numStep);

		if(incRED==0 && incDM == 0){
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteLogLike, dumper, context);
		}
			else if(incRED==1 || incDM==1 ){
#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedGPULogLike, dumper, context);
#else
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
			else if(incRED==2 || incDM==2 || incRED==3 || incDM==3 || incRED==4 || incRED==5){
#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedGPULogLike, dumper, context);
#else
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
		}
	}
	else if(doLinearFit == 1){

		//Combine all the priors into one aray: Dpriors
		int pcount=0;
		for(int i =0; i< numFitJumps+fitcount; i++){
			Dpriors[pcount][0]=-FitSig;
			Dpriors[pcount][1]=FitSig;
			pcount++;
		}

		if(numStep>0){
			for(int i =0; i < numStep; i++){
				Dpriors[pcount][0]=StepAmpPrior[0];
				Dpriors[pcount][1]=StepAmpPrior[1];
				pcount++;
				Dpriors[pcount][0]=psr[0].obsn[1].bat;
				Dpriors[pcount][1]=psr[0].obsn[psr[0].nobs-1].bat;
				pcount++;
			}
		}			



		for(int i =0; i< numEFAC; i++){
			Dpriors[pcount][0]=EFACPrior[0];
			Dpriors[pcount][1]=EFACPrior[1];
			pcount++;
		}
		for(int i =0; i< numEQUAD; i++){
			Dpriors[pcount][0]=EQUADPrior[0];
			Dpriors[pcount][1]=EQUADPrior[1];
			pcount++;
		}	
		if(incRED==1 || incRED==3){
			for(int i =0; i< numRedPL; i++){
				Dpriors[pcount][0]=AmpPrior[0];
				Dpriors[pcount][1]=AmpPrior[1];
				pcount++;
				Dpriors[pcount][0]=AlphaPrior[0];
				Dpriors[pcount][1]=AlphaPrior[1];
				pcount++;
			}
		}	
		else if(incRED==2){	
			for(int i =0; i< numRedCoeff; i++){
				Dpriors[pcount][0]=RedCoeffPrior[0];
				Dpriors[pcount][1]=RedCoeffPrior[1];
				pcount++;
			}
		}	
        else if(incRED==4){

            for(int i =0;i < 2*numRedCoeff;i++){
                    Dpriors[pcount][0]=-FourierSig*psr[p].rmsPre*pow(10.0,-6);
                    Dpriors[pcount][1]=FourierSig*psr[p].rmsPre*pow(10.0,-6);
 
                    pcount++;
            }
    
    		for(int i =0; i< numRedCoeff; i++){
                    Dpriors[pcount][0]=RedCoeffPrior[0];
                    Dpriors[pcount][1]=RedCoeffPrior[1];
                    pcount++;
            }
		}
        else if(incRED==5){

			for(int i =0;i < 2*numRedCoeff;i++){
				Dpriors[pcount][0]=-FourierSig*psr[p].rmsPre*pow(10.0,-6);
	            Dpriors[pcount][1]=FourierSig*psr[p].rmsPre*pow(10.0,-6);

		        pcount++;
			}

            Dpriors[pcount][0]=AmpPrior[0];
            Dpriors[pcount][1]=AmpPrior[1];
            pcount++;
            Dpriors[pcount][0]=AlphaPrior[0];
            Dpriors[pcount][1]=AlphaPrior[1];
            pcount++;
        }	
        if(incFloatRed>0){
        	for(int i =0; i < incFloatRed; i++){
		        Dpriors[pcount][0]=RedFreqPrior[0];
		        Dpriors[pcount][1]=RedFreqPrior[1];
		        pcount++;
		        Dpriors[pcount][0]=RedCoeffPrior[0];
		        Dpriors[pcount][1]=RedCoeffPrior[1];
		        pcount++;
		    }
		        
        }
        
        		
		if(incDM==1 || incDM==3){
			for(int i =0; i< numDMPL; i++){
				Dpriors[pcount][0]=DMAmpPrior[0];
				Dpriors[pcount][1]=DMAmpPrior[1];
				pcount++;
				Dpriors[pcount][0]=DMAlphaPrior[0];
				Dpriors[pcount][1]=DMAlphaPrior[1];
				pcount++;
			}
		}	
		else if(incDM==2){	
			for(int i =0; i< numDMCoeff; i++){
				Dpriors[pcount][0]=DMCoeffPrior[0];
				Dpriors[pcount][1]=DMCoeffPrior[1];
				pcount++;
			}
		}
		if(incFloatDM>0){
			for(int i =0; i < incFloatDM; i++){
				Dpriors[pcount][0]=DMFreqPrior[0];
				Dpriors[pcount][1]=DMFreqPrior[1];
				pcount++;
				Dpriors[pcount][0]=DMCoeffPrior[0];
				Dpriors[pcount][1]=DMCoeffPrior[1];
				pcount++;
			}
		}


		int linearNum=numFitJumps+fitcount;
		double **TNDM=new double*[psr[0].nobs];
		for(int i=0;i<psr[0].nobs;i++){
			TNDM[i]=new double[linearNum];
		}
		
					
		//Update Pulsar position to reflect values in TempoPriors which will be either the Tempo2 fit, the max, or the value set in 			custom priors
		
		paramsfitted=1;
		for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<psr[0].param[p].aSize;k++){
				if(psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel){
					((MNStruct *)context)->pulse->param[p].val[k]=TempoPriors[paramsfitted][0];
					//printf("here: %i %.10Lg %.10Lg\n",paramsfitted,TempoPriors[paramsfitted][0],TempoPriors[paramsfitted][1]);
					paramsfitted++;
	
				}
			}
		}
	

		for(int i=0;i<=psr[0].nJumps;i++){
			if(psr[0].fitJump[i] == 1){
				((MNStruct *)context)->pulse->jumpVal[i]=TempoPriors[paramsfitted][0];
				paramsfitted++;
			}
		} 
		
			
	  formBatsAll(((MNStruct *)context)->pulse,npsr);                /* Form Barycentric arrival times */
	  logdbg("calling formResiduals");
	  formResiduals(((MNStruct *)context)->pulse, npsr,1);       /* Form residuals */
		
		getDMatrix(((MNStruct *)context)->pulse, fitcount, numFitJumps, linearNum, TempoFitNums, TempoJumpNums, Dpriors, 1, 2, TNDM);
		
		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors, ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
			paramsfitted=0;
		}
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		getLinearPriors(((MNStruct *)context)->pulse, TNDM, TempoPriors, Dpriors, numFitJumps+fitcount, FitSig);
		//printf("set up priors, pcount: %i \n",pcount);
		//
                numToMargin=0;
             //   if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5)numToMargin +=2;
                for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
                        if(TempoPriors[i][2]==1){
                                numToMargin++;
                        }
                }


		if(numToMargin>0){

                        int *FitList=new int[ndims];
                        for(int i=0;i<ndims;i++){
                                FitList[i]=0;
                        }
                        numToMargin=0;
		//	if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5)numToMargin +=2;
                        for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
                                if(TempoPriors[i][2]==1){
                                        FitList[i]=1;
                                        printf("marginalising over param: %i \n",i);
                                        numToMargin++;
                                }
                        }

                        int Gsize=psr[0].nobs-numToMargin;
                        double **TNGDM=new double*[psr[0].nobs];
                        for(int i=0;i<psr[0].nobs;i++){
                                TNGDM[i]=new double[numToMargin];
                        }

                        double **TNGM=new double*[psr[0].nobs];
                        for(int i=0;i<psr[0].nobs;i++){
                                TNGM[i]=new double[psr[0].nobs];

                        }

                        getCustomDMatrix(psr, FitList, TNGDM, TempoFitNums, TempoJumpNums, Dpriors, incDM, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps);
                        //Get DMatrix for marginalisation, using T2 values if not custom or max, max if not custom, or else custom position.
                        //getMarginDMatrix(psr, fitcount, numFitJumps, numToMargin, TempoFitNums, TempoJumpNums, Dpriors, doJumpMargin, doTimeMargin, TNDM, 0);
                        makeGDesign(psr, Gsize, numToMargin, TNGM, TNGDM);

                        update_MNPriors(MNS,Dpriors, TempoPriors,2);
			update_MNGDdata(MNS, numFitJumps+fitcount, TNDM, Gsize,TNGM);
	
			context=MNS;
		

#ifdef HAVE_CULA
			double *GMatrixVec=new double[psr[0].nobs*Gsize];


			for(int g=0;g<Gsize; g++){
				for(int o=0;o<psr[0].nobs; o++){

					GMatrixVec[g*psr[0].nobs + o]=TNGM[o][g];


				}
				}

			copy_gmat_(GMatrixVec, psr[0].nobs*Gsize);

#endif

			ndims=ndims-numToMargin;
          //  if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5)ndims +=2;
            nPar=ndims;
                        
            printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff, numDMCoeff,incFloatRed,incFloatDM, fitDMModel, longname, numStep);
            
              
			if(numEFAC==0 && numEQUAD==0){
				printf("Not Fitting for white noise: Pre computing Matrices.\n");
				double **staticG=new double*[psr[0].nobs];
				for(int i=0;i<psr[0].nobs;i++){
					staticG[i]=new double[psr[0].nobs];
				}
				double statictdet=0;
				makeStaticGMatrix(psr, Gsize, TNGM, staticG, statictdet);
				update_MNstaticG(MNS, staticG, statictdet);
				printf("Matrices generated, static matrix determinant: %g \n",statictdet);
				context=MNS;
				
				
#ifdef HAVE_CULA
				double *staticGMatrixVec=new double[psr[0].nobs*psr[0].nobs];


				for(int g=0;g<psr[0].nobs; g++){
					for(int o=0;o<psr[0].nobs; o++){

						staticGMatrixVec[g*psr[0].nobs + o]=staticG[o][g];
					}
				}

				copy_staticgmat_(staticGMatrixVec, psr[0].nobs, psr[0].nobs);


#endif /* HAVE_CULA */


			}
			
			
			if(numEFAC==1 || numEQUAD==1 && numEFAC < 2 && numEQUAD < 2){
			
				printf("Fitting for only one EFAC or EQUAD: Pre computing Matrices.\n");
				double **UM=new double*[Gsize];
				for(int i=0;i<Gsize;i++){
					UM[i]=new double[psr[0].nobs];
				}
				double *SM=new double[Gsize];
				makeStaticDiagGMatrix(psr, Gsize, TNGM, UM, SM);
				update_MNstaticDiagG(MNS, UM, SM);
				printf("Matrices generated");
				context=MNS;
				
#ifdef HAVE_CULA
				double *staticUMatrixVec=new double[Gsize*psr[0].nobs];


				for(int g=0;g<psr[0].nobs; g++){
					for(int o=0;o<Gsize; o++){

						staticUMatrixVec[g*Gsize + o]=UM[o][g];
					}
				}

				copy_staticumat_(staticUMatrixVec, Gsize, psr[0].nobs);


#endif /* HAVE_CULA */


			}
			
			
//			printf("incRED: %i \n",incRED);
			if(incRED==0 && incDM ==0){

#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginGPULogLike, dumper, context);
#else
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
			else if(incRED==1 || incDM==1 ){
#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginGPULogLike, dumper, context);
#else
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
			else if(incRED==2 || incDM==2 || incRED==3 || incDM==3 || incRED==4 || incRED==5){
#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginGPULogLike, dumper, context);
#else
 				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
			}

		}
		else if(doJumpMargin == 0 && doTimeMargin == 0 ){


			update_MNPriors(MNS,Dpriors, TempoPriors,2);

			update_MNDdata(MNS, numFitJumps+fitcount,TNDM);

			context=MNS;
	
			printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff, numDMCoeff,incFloatRed,incFloatDM, fitDMModel, longname, numStep);
			
			if(incRED==0 && incDM==0){
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteLogLike, dumper, context);
		}
			else if(incRED==1 || incDM ==1 ){
#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedGPULogLike, dumper, context);
#else
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
			else if(incRED==2 || incDM==2 || incRED==3 || incDM==3 || incRED==4 || incRED==5){
#ifdef HAVE_CULA
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedGPULogLike, dumper, context);
#else
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
		}

		convertFromLinear(psr, longname, ndims, context);
	}


	readsummary(psr,longname, ndims,context,  Tempo2Fit,incRED, ndims, doTimeMargin, doJumpMargin,doLinearFit);

// 	printf("num its %i \n",psr[0].nits);
	endClock = clock();
  	printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
// 	 exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
} 





// redwards function to force linkage with library functions used by
// plugins
void
thwart_annoying_dynamic_library_stuff(int never_call_me, float or_sink)
{
  ChebyModel *cm;
  T2Predictor *t2p;
  ChebyModel_Init(cm, 0, 0);
  T2Predictor_GetPhase(t2p, 0, 0);
}
