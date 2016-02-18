#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2013 Lindley Lentati

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
#include <sys/time.h>
#include "tempo2.h"
//#include "tempo2Util.h"
#include "tempo2pred.h"
#include "tempo2pred_int.h"
//#include "T2accel.h"
#include <dlfcn.h>
// #include "T2toolkit.h"
#include <vector>
#include <algorithm>
#include <float.h>
#include "multinest.h"
#include "polychord.h"
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
#include "dgesvd.h"
#ifdef HAVE_CULA
#include <cula_lapack_device.h>
#endif /* HAVE_CULA */

#ifdef HAVE_PSRCHIVE
// psrchive stuff
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/FaradayRotation.h"
#include "Pulsar/PolnProfileStats.h"
#include "T2Observatory.h"
using namespace Pulsar;
#endif


void ephemeris_routines(pulsar *psr,int npsr);
void clock_corrections(pulsar *psr,int npsr);
void extra_delays(pulsar *psr,int npsr);

#ifdef HAVE_CULA
//extern "C" void copy_gmat_(double *G, int N);
//extern "C" void copy_floatgmat_(float *G, int N);
//extern "C" void copy_staticgmat_(double *G, int M, int N);
//extern "C" void copy_staticumat_(double *G, int M, int N);
//extern "C" void copy_staticdmat_(double **TNDM, double *D, int M, int N);
//extern "C" void copy_staticECorrmat_(double *EVec, int E, int N);
extern "C" void copy_staticTmat_(double *T, int totalsize, int Nobs);
#endif /* HAVE_CULA */


void fastephemeris_routines(pulsar *psr,int npsr)
{
vectorPulsar(psr,npsr); /* 1. Form a vector pointing at the pulsar */
//readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
//get_obsCoord(psr,npsr);   /* 3. Get Coordinate of observatory relative to Earth's centre */
//tt2tb(psr,npsr); /* Observatory/time-dependent part of TT-TB */
readEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */
//readShortEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */
}


void fastSubIntephemeris_routines(pulsar *psr,int npsr)
{
//vectorPulsar(psr,npsr); /* 1. Form a vector pointing at the pulsar */
//readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
get_obsCoord(psr,npsr);   /* 3. Get Coordinate of observatory relative to Earth's centre */
//tt2tb(psr,npsr); /* Observatory/time-dependent part of TT-TB */
//readEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */

}

void fastformBatsAll(pulsar *psr,int npsr)
{
	//clock_corrections(psr,npsr); /* Clock corrections ... */

	int dotime = 0;
	struct timeval tval_before, tval_after, tval_resultone, tval_resulttwo;
	if(dotime == 1){
		gettimeofday(&tval_before, NULL);
	}

	fastephemeris_routines(psr,npsr); /* Ephemeris routines ... */

	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to end of ephem: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

	}


	extra_delays(psr,npsr); /* Other time delays ... */

	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to end of extra: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

	}

	formBats(psr,npsr); /* Form Barycentric arrival times */

	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to end of formbats: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

	}

	secularMotion(psr,npsr);
	//updateBatsAll(psr, npsr);
	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to end of secular: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

	}
}


void fastformSubIntBatsAll(pulsar *psr,int npsr)
{
        //clock_corrections(psr,npsr); /* Clock corrections ... */
        fastSubIntephemeris_routines(psr,npsr); /* Ephemeris routines ... */
  	extra_delays(psr,npsr); /* Other time delays ... */
	formBats(psr,npsr); /* Form Barycentric arrival times */
	secularMotion(psr,npsr);

}


void SaveNewArchive(Archive *archive, std::string NewArchiveName){

	

	archive->unload(NewArchiveName);


}

void readDMVals(double *StoredDMVals){


	std::ifstream MeanProfileFile;
	std::string MeanProfileFilename = "/home/ltl21/scratch/Pulsars/ProfileData/FrequencyData/Sims/WideBandPaper/Sim4/RedSignal.dat";
	MeanProfileFile.open(MeanProfileFilename.c_str());

	for(int i = 0; i < 100; i++){
		std::string line;
		getline(MeanProfileFile,line);
		std::istringstream myStream( line );
		std::istream_iterator< double > begin(myStream),eof;
		std::vector<double> ProfileVal(begin,eof);
		StoredDMVals[i] = ProfileVal[0];
	}
	MeanProfileFile.close();

}

void SetNewAmps(Archive *archive){



	std::ifstream MeanProfileFile;
	std::string MeanProfileFilename = "/scratch/ltl21/Pulsars/ProfileData/FrequencyData/J1713+0747/results/OneProfStoc.dat";
	MeanProfileFile.open(MeanProfileFilename.c_str());

	float *readvals = new float[1024];
	for(int i = 0; i < 1024; i++){
		std::string line;
		getline(MeanProfileFile,line);
		std::istringstream myStream( line );
		std::istream_iterator< double > begin(myStream),eof;
		std::vector<double> ProfileVal(begin,eof);
		readvals[i] = ProfileVal[0];
		//printf("Read File: coeff %i = %g \n", i, readvals[i]);
		
	}
	MeanProfileFile.close();

	for(int i=0; i < archive->get_nsubint(); i++){

		Integration *subint = archive->get_Integration(i);

		int nbins = subint->get_nbin();
		int nchans = subint->get_nchan();
		int npols = subint->get_npol();

		for(int j = 0; j < nchans; j++){

			int ipol = 0; //only want total Intensity
			Profile *prof = archive->get_Profile(i, ipol, j);
			float *pvalues = prof->get_amps();

			double oneflux = 0;	
			for(int b =0; b < nbins; b++){
				oneflux = oneflux + pvalues[b];
			}  
			//float *newpvalues = new float[nbins]();
			if(oneflux != 0){prof->set_amps(readvals);}
		}
	}
}


MNStruct* init_struct(pulsar *pulseval,	 long double **LDpriorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval, int systemcountval, int numFitEFACval, int numFitEQUADval, int numFitRedCoeffval, int numFitDMCoeffval,int numFitRedPLval, int numFitDMPLval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval, int numdimsval, int incREDval, int incDMval, int incFloatDMval, int incFloatRedval, int DMFloatstartval, int RedFloatstartval, int TimeMarginVal, int JumpMarginVal, int doLinearVal, double *SampleFreqsVal, int incStepVal, char *whiteflagval, int whitemodelval, int varyRedCoeffval, int varyDMCoeffval, int yearlyDMval, int incsinusoidval, int EPolTermsval, int incGWBval,int RedPriorType,int DMPriorType,int EQUADPriorType,int EFACPriorType,int useOriginalErrors, int incShannonJitter, int incDMEvent, int incDMShapeEvent, int numDMShapeCoeff, int incBandNoise, int numFitBandNoiseCoeff, int incRedShapeEvent, int numRedShapeCoeff, int MarginRedShapeCoeff, int incDMScatterShapeEvent, int numDMScatterShapeCoeff, int incNGJitter, int incGlitch, int incGlitchTerms, int incBreakingIndex, int FitLowFreqCutoff, int uselongdouble, int incGroupNoise, int numFitGroupNoiseCoeff, int **FitForGroup, int numGroupstoFit, int *GroupNoiseFlags, int FitSolarWind, int FitWhiteSolarWind, int interpolateProfile, double InterpolatedTime, int sampler, int *GPTAnumstoccoeff, int totalshapestoccoeff, int StoreFMatrices, int incHighFreqStoc, int numNGJitter, int **FitForBand, int incProfileEvo, double EvoRefFreq, int *numEvoFitCoeff, int incWideBandNoise, int incProfileFit, int *numProfileFitCoeff, int incDMEQUAD, int FitLinearProfileWidth, double offPulseLevel, double **GroupStartTimes, int FitEvoExponent, int numProfComponents, int totalEvoFitCoeff, int totalProfileFitCoeff, int *numEvoCoeff, int totalEvoCoeff, int incWidthJitter, int JitterProfComp, int incProfileEnergyEvo, int debug, int ProfileBaselineTerms, int incProfileNoise, int ProfileNoiseCoeff, int SubIntToFit, int ChannelToFit, int NProfileEvoPoly, double *ProfCompSeps, int usecosiprior)
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
	MNS->yearlyDM=yearlyDMval;
	MNS->incsinusoid=incsinusoidval;
	MNS->FloatRedstart=RedFloatstartval;
   	MNS->FloatDMstart=DMFloatstartval;
	MNS->TimeMargin=TimeMarginVal;
	MNS->JumpMargin=JumpMarginVal;
	MNS->doLinear=doLinearVal;
	MNS->sampleFreq=SampleFreqsVal;
	MNS->incStep=incStepVal;
	MNS->whiteflag=whiteflagval;
	MNS->whitemodel=whitemodelval;
	MNS->varyRedCoeff=varyRedCoeffval;
	MNS->varyDMCoeff=varyDMCoeffval;
	MNS->EPolTerms=EPolTermsval;
	MNS->incGWB=incGWBval;
	MNS->incDMEvent=incDMEvent;
	MNS->RedPriorType=RedPriorType;
	MNS->DMPriorType=DMPriorType;
	MNS->EQUADPriorType=EQUADPriorType;
	MNS->EFACPriorType=EFACPriorType;
	MNS->useOriginalErrors=useOriginalErrors;
	MNS->incShannonJitter=incShannonJitter;
	MNS->incDMShapeEvent=incDMShapeEvent;
	MNS->numDMShapeCoeff=numDMShapeCoeff;
	MNS->incRedShapeEvent=incRedShapeEvent;
	MNS->numRedShapeCoeff=numRedShapeCoeff;
	MNS->MarginRedShapeCoeff=MarginRedShapeCoeff;
	MNS->incDMScatterShapeEvent=incDMScatterShapeEvent;
	MNS->numDMScatterShapeCoeff=numDMScatterShapeCoeff;
	MNS->incBandNoise=incBandNoise;
	MNS->numFitBandNoiseCoeff=numFitBandNoiseCoeff;
	MNS->incNGJitter=incNGJitter;
	MNS->incGlitch=incGlitch;
	MNS->incGlitchTerms=incGlitchTerms;
	MNS->incBreakingIndex=incBreakingIndex;
	MNS->FitLowFreqCutoff=FitLowFreqCutoff;
	MNS->uselongdouble=uselongdouble;
	MNS->incGroupNoise=incGroupNoise;
	MNS->numFitGroupNoiseCoeff=numFitGroupNoiseCoeff;	
	MNS->FitForGroup=FitForGroup;
	MNS->numGroupstoFit=numGroupstoFit;
	MNS->GroupNoiseFlags=GroupNoiseFlags;
	MNS->FitSolarWind = FitSolarWind;
	MNS->FitWhiteSolarWind = FitWhiteSolarWind;
	MNS->InterpolateProfile = interpolateProfile;
	MNS->InterpolatedTime = InterpolatedTime*pow(10.0, -9);
	MNS->sampler = sampler;
	MNS->numshapestoccoeff = GPTAnumstoccoeff;
	MNS->storeFMatrices = StoreFMatrices;
	MNS->incHighFreqStoc = incHighFreqStoc;
	MNS->incNGJitter=numNGJitter;
	MNS->FitForBand = FitForBand;
	MNS->incProfileEvo = incProfileEvo;
	MNS->EvoRefFreq = EvoRefFreq;
	MNS->numEvoFitCoeff = numEvoFitCoeff;
	MNS->incWideBandNoise = incWideBandNoise;
	MNS->incProfileFit = incProfileFit;
	MNS->numProfileFitCoeff = numProfileFitCoeff;
	MNS->incDMEQUAD = incDMEQUAD;
	MNS->FitLinearProfileWidth = FitLinearProfileWidth;
	MNS->offPulseLevel = offPulseLevel;
	MNS->GroupStartTimes=GroupStartTimes;
	MNS->FitEvoExponent = FitEvoExponent;
	MNS->numProfComponents = numProfComponents;
	MNS->totalEvoFitCoeff = totalEvoFitCoeff;
	MNS->totalProfileFitCoeff = totalProfileFitCoeff;
	MNS->totalEvoCoeff = totalEvoCoeff;
	MNS->numEvoCoeff = numEvoCoeff;
	MNS->incWidthJitter = incWidthJitter;
	MNS->JitterProfComp = JitterProfComp;
	MNS->totalshapestoccoeff = totalshapestoccoeff;
	MNS->incProfileEnergyEvo = incProfileEnergyEvo;
	MNS->debug = debug;
	MNS->ProfileBaselineTerms = ProfileBaselineTerms;
	MNS->incProfileNoise = incProfileNoise;
	MNS->ProfileNoiseCoeff = ProfileNoiseCoeff;
	MNS->SubIntToFit = SubIntToFit;
	MNS->ChannelToFit = ChannelToFit;
	MNS->NProfileEvoPoly = NProfileEvoPoly;
	MNS->ProfCompSeps=ProfCompSeps;
	MNS->usecosiprior=usecosiprior;


        MNS->Tspan = 0;
        MNS->TimetoMargin = 0;
        MNS->totCoeff = 0;
        MNS->totRedShapeCoeff = 0;
        MNS->totalsize = 0;
	



	return MNS;
}

void printPriors(pulsar *psr, long double **TempoPriors, double **Dpriors, int incEFAC, int incEQUAD, int incRED, int incDM, int numRedCoeff, int numDMCoeff, int numFloatRed, int numFloatDM, int fitDMModel, std::string longname, int incStep, int varyRedCoeff, int varyDMCoeff, int yearlyDM, int incsinusoid, int numEPolTerms, int incShannonJitter, void *context){


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
	if(TempoPriors[0][2] == 0){	
	printf("Prior on Phase : %.25Lg -> %.25Lg\n",TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
	}
	paramsfitted++;
	
	for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<psr[0].param[p].aSize;k++){
				if(psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel){
					if(TempoPriors[paramsfitted][2] == 0){			
					printf("Prior on %s : %.25Lg -> %.25Lg\n",psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
					}
					
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
	
			if(TempoPriors[paramsfitted][2] == 0){
	
				printf("Prior on Jump %i : %.25Lg -> %.25Lg\n",jumpsfitted+1, TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
			
				getdistparamnames << getdistlabel;
				getdistparamnames << " ";
				getdistparamnames <<  "Jump";
				getdistparamnames << jumpsfitted+1;
				getdistparamnames << "\n";
				getdistlabel++;

			}

			paramsfitted++;
			jumpsfitted++;
		}
	} 

	if(incStep>0){
		for(int i =0; i < incStep; i++){
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "StepG";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "StepT";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "StepA";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

			printf("Prior for Step Group: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior for Step Time: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior for Step Amp: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	

		}
	}
	

	if(((MNStruct *)context)->incGlitch>0){
		for(int i =0; i < ((MNStruct *)context)->incGlitch; i++){
			                                getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "GlitchMJD";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;

			printf("Prior for Glitch Epoch: %i %g %g \n",i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;

			for(int j=0; j < ((MNStruct *)context)->incGlitchTerms; j++){
				getdistparamnames << getdistlabel;
				getdistparamnames << " ";
				getdistparamnames <<  "GlitchF"<<j<<" ";
				getdistparamnames <<  i+1;
				getdistparamnames << "\n";
				getdistlabel++;

				printf("Prior for Glitch F%i: %i %g %g \n",j, i,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
                }
        }


		
	if(incEFAC>0){
		for(int n=1;n<=numEPolTerms; n++){
			int EFACnum=1;
			for(int i =0;i<incEFAC;i++){
				printf("Prior on EPol %i EFAC %i : %.5g -> %.5g\n",n, EFACnum, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				
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
	
	if(incShannonJitter>0){
		int SQUADnum=1;	
		for(int i =0;i<incShannonJitter;i++){
			printf("Prior on SQUAD %i: %.5g -> %.5g\n",SQUADnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "SQUAD";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			getdistlabel++;
				
			paramsfitted++;	
			SQUADnum++;
		}
	}

	if(((MNStruct *)context)->incHighFreqStoc>0){
		int HF=1;	
		for(int i =0;i<((MNStruct *)context)->incHighFreqStoc;i++){
			printf("Prior on High Frequency Profile Stochasticity %i: %.5g -> %.5g\n",HF,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "HF";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			getdistlabel++;
				
			paramsfitted++;	
			HF++;
		}
	}

	if(((MNStruct *)context)->FitSolarWind==1){
		printf("Prior on Solar Wind %.5g -> %.5g\n", Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "SW";
		getdistparamnames << "\n";
		getdistlabel++;

		paramsfitted++;
	}	


        if(((MNStruct *)context)->FitWhiteSolarWind==1){
                printf("Prior on White Solar Wind %.5g -> %.5g\n", Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                getdistparamnames << getdistlabel;
                getdistparamnames << " ";
                getdistparamnames <<  "WSW";
                getdistparamnames << "\n";
                getdistlabel++;

                paramsfitted++;
        }


	if(((MNStruct *)context)->incNGJitter>0){
		int ECORRnum=1;	
		for(int i =0;i<((MNStruct *)context)->incNGJitter;i++){
			printf("Prior on ECORR %i: %.5g -> %.5g\n",ECORRnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "ECORR";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			getdistlabel++;
				
			paramsfitted++;	
			ECORRnum++;
		}
	}

	if(((MNStruct *)context)->incDMEQUAD == 1){

		printf("Prior on DMEQUAD: %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "DMEQ";
		getdistparamnames << "\n";
		getdistlabel++;
			
		paramsfitted++;	

	}

	if(((MNStruct *)context)->incWidthJitter == 1){

		printf("Prior on WidthJitter: %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "WJ";
		getdistparamnames << "\n";
		getdistlabel++;
			
		paramsfitted++;	

	}


	if(((MNStruct *)context)->incDMEvent != 0){
                for(int i =0; i < ((MNStruct *)context)->incDMEvent; i++){
		
			printf("Prior on DM Event %i Start Point %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;		
			printf("Prior on DM Event %i Length %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Event %i Log Amp %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Event %i Spectral Index %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Event %i Offset Term %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Event %i Linear Term %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Event %i Quadratic Term %g -> %g \n", i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	

			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "Start"<< "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "Length"<< "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "LogAmp"<< "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "Index"<< "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "Offset"<< "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "Linear"<< "\n";
			getdistlabel++;
			getdistparamnames << getdistlabel <<  " DM Event"<<i+1<< "Quad"<< "\n";
			getdistlabel++;

		}
	}


	if(((MNStruct *)context)->FitLowFreqCutoff > 0){
		printf("Prior on Red LF Cutoff %g -> %g \n", Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;
		getdistparamnames << getdistlabel <<  "LFC\n";
		getdistlabel++;
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

		if(varyRedCoeff==1){
	
			printf("Varying RedCoeffs : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedPLCoeff";
			getdistparamnames << "\n";
			getdistlabel++;
	
			paramsfitted++;	
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedPLCoeffAmp";
			getdistparamnames << "\n";
			getdistlabel++;	
		}
			
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
                printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                paramsfitted++;

                printf("Prior on Red Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                paramsfitted++;
        
                printf("Prior on Red Noise Corner Freq : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
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

                getdistparamnames << getdistlabel;
                getdistparamnames << " ";
                getdistparamnames <<  "RedCorner";
                getdistparamnames << "\n";
                getdistlabel++; 

                if(varyRedCoeff==1){

                        printf("Varying RedCoeffs : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                        paramsfitted++;
                        getdistparamnames << getdistlabel;
                        getdistparamnames << " ";
                        getdistparamnames <<  "RedPLCoeff";
                        getdistparamnames << "\n";
                        getdistlabel++;

                        paramsfitted++;
                        getdistparamnames << getdistlabel;
                        getdistparamnames << " ";
                        getdistparamnames <<  "RedPLCoeffAmp";
                        getdistparamnames << "\n";
                        getdistlabel++;
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

	if(((MNStruct *)context)->incGWB == 1){
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "GWBAmp ";
		getdistparamnames << "\n";
		getdistlabel++;


		printf("Prior on Log GWB Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;
	}


	if(((MNStruct *)context)->incRedShapeEvent >0){
		for(int i =0; i < ((MNStruct *)context)->incRedShapeEvent; i++){
			getdistparamnames << getdistlabel << " RedShape"<<i+1<<"Pos\n";
			getdistlabel++;	
			getdistparamnames << getdistlabel << " RedShape"<<i+1<<"Width\n";
			getdistlabel++;	

			printf("Prior on Red Shape Event %i Start : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on Red Shape Event %i Width : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	

			if(((MNStruct *)context)->MarginRedShapeCoeff == 0){

				for(int c =0; c< ((MNStruct *)context)->numRedShapeCoeff; c++){
					getdistparamnames << getdistlabel << " RedShape"<<i+1<<"C"<<c+1<<"\n";
					getdistlabel++;	
					printf("Prior on Red Shape Event %i Coeff %i : %.5g -> %.5g\n",i+1,c+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;	
				}
			}
		}
	}

        if(incsinusoid==1){
	 	getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "sineAmp ";
		getdistparamnames << "\n";
		getdistlabel++;
	 	getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "sinePhase ";
		getdistparamnames << "\n";
		getdistlabel++;
	 	getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "sineFreq ";
		getdistparamnames << "\n";
		getdistlabel++;


		printf("Prior on Sine Log Amplitude : %.5g -> %.5g\n", Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;		
		printf("Prior on Sine phase : %.5g -> %.5g\n", Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;
		printf("Prior on Log Sine Freq : %.5g -> %.5g\n", Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
		paramsfitted++;
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

		if(varyDMCoeff==1){
	
			printf("Varying DMCoeffs : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMPLCoeff";
			getdistparamnames << "\n";
			getdistlabel++;
	
			paramsfitted++;	
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMPLCoeffAmp";
			getdistparamnames << "\n";
			getdistlabel++;	
		}
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

        else if(incDM==5){

                for(int i =0;i < 2*numDMCoeff;i++){
                                printf("Prior on DM fourier coefficients: %i %g %g \n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);

                                getdistparamnames << getdistlabel;
                                getdistparamnames << " ";
                                getdistparamnames <<  "DMFourierC";
                                getdistparamnames <<  i+1;
                                getdistparamnames << "\n";
                                getdistlabel++;


                                paramsfitted++;
                }

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

                        printf("Prior on DM Variations Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                        paramsfitted++;
                        printf("Prior on DM Variations Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                        paramsfitted++;
        }

        else if(incDM==6){

                for(int i =0;i < 2*numDMCoeff;i++){
                                printf("Prior on DM fourier coefficients: %i %g %g \n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);

                                getdistparamnames << getdistlabel;
                                getdistparamnames << " ";
                                getdistparamnames <<  "DMFourierC";
                                getdistparamnames <<  i+1;
                                getdistparamnames << "\n";
                                getdistlabel++;


                                paramsfitted++;
                }

                for(int i =0;i < numDMCoeff;i++){
                                printf("Prior on DM PS coefficients: %i %g %g \n",i, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);

                                getdistparamnames << getdistlabel;
                                getdistparamnames << " ";
                                getdistparamnames <<  "DMFourierPC";
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

	if(((MNStruct *)context)->incDMShapeEvent >0){
		for(int i =0; i < ((MNStruct *)context)->incDMShapeEvent; i++){
			getdistparamnames << getdistlabel << " DMShape"<<i+1<<"Pos\n";
			getdistlabel++;	
			getdistparamnames << getdistlabel << " DMShape"<<i+1<<"Width\n";
			getdistlabel++;	

			printf("Prior on DM Shape Event %i Start : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Shape Event %i Width : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	


			for(int c =0; c< ((MNStruct *)context)->numDMShapeCoeff; c++){
				getdistparamnames << getdistlabel << " DMShape"<<i+1<<"C"<<c+1<<"\n";
				getdistlabel++;	
				printf("Prior on DM Shape Event %i Coeff %i : %.5g -> %.5g\n",i+1,c+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;	
			}
		}
	}

	if(((MNStruct *)context)->incDMScatterShapeEvent >0){
		for(int i =0; i < ((MNStruct *)context)->incDMScatterShapeEvent; i++){
			getdistparamnames << getdistlabel << " DMScatterShape"<<i+1<<"Pos\n";
			getdistlabel++;	
			getdistparamnames << getdistlabel << " DMScatterShape"<<i+1<<"Width\n";
			getdistlabel++;	
			getdistparamnames << getdistlabel << " DMScatterShape"<<i+1<<"Freq\n";
			getdistlabel++;	

			printf("Prior on DM Scatter Shape Event %i Start : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Scatter Shape Event %i Width : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	
			printf("Prior on DM Scatter Shape Event %i Frequency Dependence : %.5g -> %.5g\n",i+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;	

			for(int c =0; c< ((MNStruct *)context)->numDMScatterShapeCoeff; c++){
				getdistparamnames << getdistlabel << " DMScatterShape"<<i+1<<"C"<<c+1<<"\n";
				getdistlabel++;	
				printf("Prior on DM Scatter Shape Event %i Coeff %i : %.5g -> %.5g\n",i+1,c+1,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;	
			}
		}
	}

	if(((MNStruct *)context)->incBandNoise > 0 ){

		for(int b = 0; b < ((MNStruct *)context)->incBandNoise; b++){

			getdistparamnames << getdistlabel << " BNAmp\n";
			getdistlabel++;	
	
			getdistparamnames << getdistlabel << " BNSpec\n";
			getdistlabel++;
	

			
			printf("Prior on %i - %i Noise Log Amplitude : %.5g -> %.5g\n", ((MNStruct *)context)->FitForBand[b][0], ((MNStruct *)context)->FitForBand[b][1], Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;		
			printf("Prior on %i - %i Noise Slope : %.5g -> %.5g\n", ((MNStruct *)context)->FitForBand[b][0], ((MNStruct *)context)->FitForBand[b][1], Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;
		}

	}




	if(((MNStruct *)context)->incGroupNoise > 0 ){

		for(int i = 0 ; i <((MNStruct *)context)->incGroupNoise; i++){

			if(((MNStruct *)context)->FitForGroup[i][0] == -1){
				getdistparamnames << getdistlabel << " Group\n";
				getdistlabel++;	
			}
                        if(((MNStruct *)context)->FitForGroup[i][1] == 1){
                                getdistparamnames << getdistlabel << " GroupStart\n";
                                getdistlabel++; 
                                getdistparamnames << getdistlabel << " GroupLength\n";
                                getdistlabel++;

                        }

			getdistparamnames << getdistlabel << " GroupAmp\n";
			getdistlabel++;	
	
			getdistparamnames << getdistlabel << " GroupSpec\n";
			getdistlabel++;
	

			if(((MNStruct *)context)->FitForGroup[i][0] == -1){
				printf("Prior on Group Noise Group : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;	
			}
                        if(((MNStruct *)context)->FitForGroup[i][1] == 1){
                                printf("Prior on Group Noise Start : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                                paramsfitted++;
                                printf("Prior on Group Noise Length : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
                                paramsfitted++;

                        }


			printf("Prior on Group Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
			paramsfitted++;		
			printf("Prior on Group Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
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

	for(int c = 0; c < ((MNStruct *)context)->numProfComponents; c++){
		for(int i =0; i<((MNStruct *)context)->numshapestoccoeff[c]; i++){
			printf("Prior on Profile Shape Stochasticitiy Coeff %i for comp %i Log Amplitude : %.5g -> %.5g\n",i, c, -10.0, 1.0);

			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "PS";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;


			paramsfitted++;
		}
	}

	for(int c = 0; c < ((MNStruct *)context)->numProfComponents; c++){
		for(int i =0; i<((MNStruct *)context)->numProfileFitCoeff[c]; i++){
			printf("Prior on Profile Fit Coeff %i for comp %i: %.5g -> %.5g\n",i,-0.1, 0.1, c);

			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "PF";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			getdistlabel++;


			paramsfitted++;
		}
	}

	if(((MNStruct *)context)->FitLinearProfileWidth == 1){

			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "PW";
			getdistparamnames << "\n";
			getdistlabel++;


			paramsfitted++;
	}



	if(((MNStruct *)context)->FitEvoExponent == 1){

			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "EvoExp";
			getdistparamnames << "\n";
			getdistlabel++;


			paramsfitted++;
	}

	for(int c = 0; c < ((MNStruct *)context)->numProfComponents; c++){
		for(int i =0; i<((MNStruct *)context)->numEvoFitCoeff[c]; i++){
			for(int j =0; j<((MNStruct *)context)->NProfileEvoPoly; j++){
				printf("Prior on Profile Evolution Coeff %i %i : %.5g -> %.5g\n",i,j, -0.1, 0.1);

				getdistparamnames << getdistlabel;
				getdistparamnames << " ";
				getdistparamnames <<  "PE";
				getdistparamnames <<  i+1;
				getdistparamnames << "\n";
				getdistlabel++;


				paramsfitted++;
			}
		}
	}	

	if(((MNStruct *)context)->incProfileEvo == 2){
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "EPW";
			getdistparamnames << "\n";
			getdistlabel++;


			paramsfitted++;
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
				if(p==param_stig){
					long double minprior=priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1];
					long double maxprior=priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][1]*priorsval[paramsfitted][1];

					if(minprior < 0 && priorsval[paramsfitted][2]==0 ){
						long double newprior=-priorsval[paramsfitted][0]/priorsval[paramsfitted][1];
						DPriorsval[paramsfitted+linearPriors][0]=(double)newprior;
                                                printf("Prior on %s updated to be physical (was <0) : %.25Lg -> %.25Lg\n", MNS->pulse->param[p].shortlabel[k], priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1],priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][1]*priorsval[paramsfitted][1]);
					}
					if(maxprior > 1 && priorsval[paramsfitted][2]==0 ){
                                                long double newprior=(1.0-priorsval[paramsfitted][0])/priorsval[paramsfitted][1];
                                                DPriorsval[paramsfitted+linearPriors][1]=(double)newprior;
                                                printf("Prior on %s updated to be physical (was >1) : %.25Lg -> %.25Lg\n", MNS->pulse->param[p].shortlabel[k], priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1],priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][1]*priorsval[paramsfitted][1]);
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
//printf("enter dump %i %i \n", nSamples, nPar);	
//	int i, j;
//	
//	double **postdist=new double*[nSamples];
//	for(int m1=0; m1<nSamples; m1++){
//		postdist[m1]=new double[nPar+2];
//	}
//	for( i = 0; i < nPar + 2; i++ ){
//		printf("%i\n",i);
//		for( j = 0; j < nSamples; j++ )
//			postdist[j][i] = posterior[0][i * nSamples + j];
//	}
	
//printf("mid dump \n");	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
//	double pLivePts[nlive][nPar + 1];
//	for( i = 0; i < nPar + 1; i++ )
//		for( j = 0; j < nlive; j++ )
//			pLivePts[j][i] = physLive[0][i * nlive + j];
//
////	printf("exit dump\n");
//	 for(int m1=0; m1<nSamples; m1++){
  //              free(postdist[m1]);
    //    }
//	free(postdist);

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
	time_t rawstarttime, rawstoptime;
	struct tm *rawstarttimeinfo;
	struct tm *rawstoptimeinfo;
	const char *CVS_verNum = "$Revision: 1.28 $";
	int numFitJumps;
	int numToMargin=0;
	int myrank = 0;
	int nprocs = 0;



  printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome to redistribute it\n");
  printf("under conditions of GPL license.\n\n");

  //
  startClock = clock();

	time(&rawstarttime);
	rawstarttimeinfo = localtime(&rawstarttime);
	time_t starttime = mktime(rawstarttimeinfo);
	printf ( "The start date/time is: %s", asctime(rawstarttimeinfo));
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

	int useGPUS;
	int uselongdouble =0;
	int GPTA=0;
	int FixProfile = 0;
	int FitTemplate = 0;
	int interpolateProfile = 0;
	double InterpolatedTime = 1;

	int StoreFMatrices = 1;
	char root[100]; 
	int numTempo2its;
	int doLinearFit;
	int doMax;
	int incEFAC;
	int numEPolTerms;
	int incEQUAD;

	int swdims=0;
	int FitSolarWind=0;
	int FitWhiteSolarWind=0;
	double *SolarWindPrior = new double[2];
	double *WhiteSolarWindPrior = new double[2];

	int incRED;
	int FitLowFreqCutoff=0;
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
	double *EPolPrior;
	double *EQUADPrior;
	double *AlphaPrior;
	double *AmpPrior;
	double *DMAlphaPrior;
	double *DMAmpPrior;
	double *DMFreqPrior;
	double *RedFreqPrior;
	double numRedCoeff;
	double numDMCoeff;
	int numRedPL;
	int numDMPL;
	double *RedCoeffPrior;
	double *DMCoeffPrior;
	double FourierSig;
	double *SampleFreq;
	int numEFAC=0;
	int numEQUAD=0;
	int numSQUAD=0;
	int incBreakingIndex=0;
	int numStep=0;
	
	int incGlitch=0;
	int incGlitchTerms=0;
	double GlitchFitSig=0;
	int Glitchdims=0;

	double *StepAmpPrior;
	double *StepTimePrior;
	char wflag[100];
	int whitemodel;
	int varyRedCoeff;
	int varyDMCoeff;
	int yearlyDM;
	int incsinusoid;
	int RedPriorType;
	int DMPriorType;
	int EQUADPriorType;
	int EFACPriorType;
	int useOriginalErrors;
	int incShannonJitter;
	int incNGJitter;
	int numNGJitter=0;

	
	char *Type = new char[100];
	char *WhiteName = new char[100];
	EFACPrior=new double[2];
	EPolPrior=new double[2];
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

	int incGWB;
	double *GWBAmpPrior;
	GWBAmpPrior=new double[2];

	int incDMEvent;
	double *DMEventStartPrior;
	double *DMEventLengthPrior;
	DMEventStartPrior=new double[2];
	DMEventLengthPrior=new double[2];

	int incDMShapeEvent;
	int numDMShapeCoeff;
	double *DMShapeCoeffPrior;
	DMShapeCoeffPrior=new double[2];

	int incRedShapeEvent;
	int numRedShapeCoeff;
	int MarginRedShapeCoeff;
	double *RedShapeCoeffPrior;
	RedShapeCoeffPrior=new double[2];

	int incDMScatterShapeEvent;
	int numDMScatterShapeCoeff;
	double *DMScatterShapeCoeffPrior;
	DMScatterShapeCoeffPrior=new double[2];


	int incBandNoise;
	int numBandNoiseCoeff = 0;
	double *BandNoiseAmpPrior;
	double *BandNoiseAlphaPrior;
	BandNoiseAmpPrior=new double[2];
	BandNoiseAlphaPrior=new double[2];

	char *GroupNoiseName = new char[100];
	char GroupNoiseSysFlag[100];
	int **FitForGroup=0;
	int incGroupNoise=0;
	int numGroupstoFit=0;
	int numGroupTimestoFit=0;
	int numGroupCoeff = 0;
	double *GroupNoiseAmpPrior;
	double *GroupNoiseAlphaPrior;	
	GroupNoiseAmpPrior=new double[2];
	GroupNoiseAlphaPrior=new double[2];

	int incHighFreqStoc = 0;
	double *HighFreqStocPrior = new double[2];

	int incProfileEvo = 0;
	int NProfileEvoPoly = 0;
	int FitEvoExponent = 0;
	double EvoRefFreq = 0;
	double *ProfileEvoPrior = new double[2];

	int incProfileFit = 0;
	double *ProfileFitPrior = new double[2];
	
	int incWideBandNoise = 0;

	int incDMEQUAD = 0;
	double *DMEQUADPrior = new double[2];

	int FitLinearProfileWidth = 0;
	double *LinearProfileWidthPrior = new double[2];

	double offPulseLevel=0;

	char *ProfFile = new char[100];
	char ProfFileName[100]; 

	int numProfComponents = 0;

	int incWidthJitter = 0;
	double *WidthJitterPrior = new double[2];	

	int JitterProfComp = 0;

	int incProfileEnergyEvo = 0;
	double *ProfileEnergyEvoPrior = new double[2];

	int debug = 0;
	int ProfileBaselineTerms = 0;

	int  pnoisedims = 0;
	int incProfileNoise = 0;
	int ProfileNoiseCoeff = 0;
	double *ProfileNoiseAmpPrior = new double[2];	
	double *ProfileNoiseSpecPrior = new double[2];


	int SubIntToFit = 0;
	int ChannelToFit = 0;

	int usecosiprior = 0;

	setupparams(useGPUS, Type, numTempo2its, doLinearFit, doMax, incEFAC, numEPolTerms, incEQUAD, incRED, incDM, doTimeMargin, doJumpMargin, FitSig, customPriors, EFACPrior, EPolPrior, EQUADPrior, AlphaPrior, AmpPrior, DMAlphaPrior, DMAmpPrior, numRedCoeff, numDMCoeff, numRedPL, numDMPL, RedCoeffPrior, DMCoeffPrior, incFloatDM, DMFreqPrior, yearlyDM, incsinusoid, incFloatRed, RedFreqPrior, FourierSig, numStep, StepAmpPrior, WhiteName,whitemodel, varyRedCoeff, varyDMCoeff, incGWB, GWBAmpPrior, RedPriorType, DMPriorType, EQUADPriorType, EFACPriorType, useOriginalErrors, incShannonJitter, incDMEvent, DMEventStartPrior, DMEventLengthPrior,incDMShapeEvent, numDMShapeCoeff, DMShapeCoeffPrior, incRedShapeEvent, numRedShapeCoeff, MarginRedShapeCoeff, RedShapeCoeffPrior, incDMScatterShapeEvent, numDMScatterShapeCoeff, DMScatterShapeCoeffPrior,incBandNoise, numBandNoiseCoeff, BandNoiseAmpPrior, BandNoiseAlphaPrior, incNGJitter, incGlitch, incGlitchTerms, GlitchFitSig, incBreakingIndex, FitLowFreqCutoff, uselongdouble, incGroupNoise, numGroupCoeff, GroupNoiseAmpPrior, GroupNoiseAlphaPrior, FitSolarWind, FitWhiteSolarWind, SolarWindPrior, WhiteSolarWindPrior,  GPTA,  GroupNoiseName, FixProfile, FitTemplate, interpolateProfile, InterpolatedTime, StoreFMatrices, incHighFreqStoc, HighFreqStocPrior, incProfileEvo, EvoRefFreq, ProfileEvoPrior, FitEvoExponent, incWideBandNoise, incProfileFit, ProfileFitPrior, FitLinearProfileWidth, LinearProfileWidthPrior, incDMEQUAD, DMEQUADPrior, offPulseLevel,ProfFile,numProfComponents,incWidthJitter,WidthJitterPrior,JitterProfComp,incProfileEnergyEvo,ProfileEnergyEvoPrior, debug, ProfileBaselineTerms, incProfileNoise, ProfileNoiseCoeff, ProfileNoiseAmpPrior, ProfileNoiseSpecPrior, SubIntToFit, ChannelToFit, NProfileEvoPoly, usecosiprior); 



	FitForGroup = new int*[incGroupNoise];
	for(int i =0; i < incGroupNoise; i++){
		FitForGroup[i] = new int[6];
		for(int j=0; j < 6; j++){
			FitForGroup[i][j] = 0;
		}
	}


	int **FitForBand = new int*[incBandNoise];
	for(int i =0; i < incBandNoise; i++){
		FitForBand[i] = new int[4];
		for(int j=0; j < 4; j++){
			FitForBand[i][j] = 0;
		}
	}

	GetGroupsToFit(incGroupNoise, FitForGroup, incBandNoise, FitForBand);
	for(int i =0; i < incGroupNoise; i++){
//		printf("Check Groups %i %i \n", i, FitForGroup[i]);
		if(FitForGroup[i][0] ==-1)numGroupstoFit++;
		if(FitForGroup[i][1] == 1)numGroupTimestoFit++;
	}
	printf("Fitting for %i group noise terms, and for %i groups, %i Times \n", incGroupNoise, numGroupstoFit, numGroupTimestoFit);
	printf("Fitting for %i band noise terms \n", incBandNoise);

	int *GPTAnumshapecoeff;
	int *numEvoCoeff;
	int *numEvoFitCoeff;
	int *numProfileFitCoeff;
	int *GPTAnumstocshapecoeff;
	double *ProfCompSeps;

	int totalEvoCoeff = 0;
	int totalEvoFitCoeff = 0;
	int totalProfileFitCoeff = 0;
	int totalshapestoccoeff = 0;
	if(GPTA ==1){
		GPTAnumshapecoeff= new int[numProfComponents];
		numEvoCoeff = new int[numProfComponents];
		numEvoFitCoeff = new int[numProfComponents];
		numProfileFitCoeff = new int[numProfComponents];
		GPTAnumstocshapecoeff  = new int[numProfComponents];
		ProfCompSeps   = new double[numProfComponents];
		for(int i = 0; i < numProfComponents; i++){
			GPTAnumshapecoeff[i] = 0;
			numEvoCoeff[i] = 0;
        	        numEvoFitCoeff[i] = 0;
	                numProfileFitCoeff[i] = 0;
			GPTAnumstocshapecoeff[i] = 0;
			ProfCompSeps[i] = 0;
		}


		GetProfileFitInfo(numProfComponents, GPTAnumshapecoeff, numProfileFitCoeff, numEvoCoeff, numEvoFitCoeff, GPTAnumstocshapecoeff, ProfCompSeps);
		for(int i = 0; i < numProfComponents; i++){
			totalEvoFitCoeff +=  numEvoFitCoeff[i];
			totalProfileFitCoeff += numProfileFitCoeff[i];
			totalEvoCoeff += numEvoCoeff[i];
			totalshapestoccoeff += GPTAnumstocshapecoeff[i];

			if(i==0 && numProfileFitCoeff[0] > 0){
				numProfileFitCoeff[0]--;
				totalProfileFitCoeff--;
			}

			if(i==0 && numEvoCoeff[0] > 0){
				numEvoCoeff[0]--;
				totalEvoCoeff--;

				if(numEvoFitCoeff[0] > 0){
					numEvoFitCoeff[0]--;
					totalEvoFitCoeff--;
				}
			}

			if(i==0 && GPTAnumstocshapecoeff[0] > 0){
				GPTAnumstocshapecoeff[0]--;
				totalshapestoccoeff--;
			}


			printf("including %i shapelets in component %i Fit: %i Evo: %i Stoc: %i\n", GPTAnumshapecoeff[i], i, numProfileFitCoeff[i], numEvoFitCoeff[i], GPTAnumstocshapecoeff[i]);
		}


		if(incProfileEvo != 1){
			NProfileEvoPoly = 0;
			for(int i = 0; i < numProfComponents; i++){
				numEvoFitCoeff[i]=0;
			}
		}
		if(incProfileFit == 0){
			for(int i = 0; i < numProfComponents; i++){
				numProfileFitCoeff[i]=0;
			}	
		}
	}



  formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
  logdbg("calling formResiduals");
  formResiduals(psr,npsr,1);       /* Form residuals */


/*Work out data time span to get maximum number of coefficients*/
	
        double start,end;
        int go=0;
        for (int i=0;i<psr[0].nobs;i++)
          {
            if (psr[0].obsn[i].deleted==0)
              {
                if (go==0)
                  {
                    go = 1;
                    start = (double)psr[0].obsn[i].bat;
                    end  = start;
                  }
                else
                  {
                    if (start > (double)psr[0].obsn[i].bat)
                      start = (double)psr[0].obsn[i].bat;
                    if (end < (double)psr[0].obsn[i].bat)
                      end = (double)psr[0].obsn[i].bat;

                  }
              }
          }





        double maxtspan=1*(end-start);

	if(maxtspan < 1){
		maxtspan = maxtspan*24*60;
		printf("Assume less than a day, Tspan is now in minutes\n");
	}

        int mindays =int(floor(1+2*maxtspan/psr[0].nobs));
        int mincoeff = int(floor(1+maxtspan/mindays));

        int Reddaysincoeffs=int(floor(maxtspan/numRedCoeff));
	int DMdaysincoeffs=int(floor(maxtspan/numDMCoeff));
	int BandDMdaysincoeffs=int(floor(maxtspan/numBandNoiseCoeff));
	int Groupdaysincoeffs=int(floor(maxtspan/numGroupCoeff));

        if(numRedCoeff < mindays){
                numRedCoeff=mincoeff;
		if(FitLowFreqCutoff  > 0){
			numRedCoeff=mincoeff+10;//int(Reddaysincoeffs);
		}
        }
        else{
                numRedCoeff=int(Reddaysincoeffs);//Reddaysincoeffs;
        }

        if(numDMCoeff < mindays){
                numDMCoeff=mincoeff;
        }
        else{
                numDMCoeff=int(DMdaysincoeffs);//DMdaysincoeffs;
        }

        if(numBandNoiseCoeff < mindays){
                numBandNoiseCoeff=mincoeff;
        }
        else{
                numBandNoiseCoeff=int(BandDMdaysincoeffs);//DMdaysincoeffs;
        }

        if(numGroupCoeff < mindays){
                numGroupCoeff=mincoeff;
        }
        else{
                numGroupCoeff=int(Groupdaysincoeffs);//DMdaysincoeffs;
        }


	if(incRED == 0 && incGWB == 0)numRedCoeff=0;
	if(incDM == 0)numDMCoeff=0;	
	SampleFreq=new double[int(numRedCoeff+numDMCoeff)];
	setFrequencies(SampleFreq,numRedCoeff, numDMCoeff, 0, 0, 1, 1, 1, 1);
	





  printf("Num T2 its %i \n", numTempo2its);
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


	
	if(incProfileNoise == 1)pnoisedims = 2;
	if(incProfileNoise == 2)pnoisedims = ProfileNoiseCoeff;
	if(incProfileNoise == 3)pnoisedims = ProfileNoiseCoeff;

	if(FitSolarWind == 1)swdims++;
	if(FitWhiteSolarWind == 1)swdims++;

	if(incRED==0)Reddims=0;
	if(incRED==1)Reddims=2;
	if(incRED==2)Reddims=numRedCoeff;
	if(incRED==3)Reddims=2*numRedPL;
	if(incRED==4)Reddims=3*numRedPL;
	if(incRED==5)Reddims=2*numRedCoeff+2;
	if(incGWB==1)Reddims+=1;
        if(incFloatRed>0)Reddims+=2*incFloatRed;
	if(incRedShapeEvent>0 && MarginRedShapeCoeff == 0){ Reddims+=incRedShapeEvent*(2+numRedShapeCoeff);}
	if(incRedShapeEvent>0 && MarginRedShapeCoeff == 1){ Reddims+=incRedShapeEvent*(2);}
	if(FitLowFreqCutoff>0)Reddims+=1;

	if(incDM==0)DMdims=0;
	if(incDM==1)DMdims=2;
	if(incDM==2)DMdims=numDMCoeff;
	if(incDM==3)DMdims=2*numDMPL;
	if(incDM==4)DMdims=3*numDMCoeff;
	if(incDM==5)DMdims=2*numDMCoeff+2;
	if(incDM==6)DMdims=2*numDMCoeff+numDMCoeff;

	if(incBandNoise > 0) DMdims+= 2*incBandNoise;
	if(incGroupNoise > 0) Reddims += 2*incGroupNoise + numGroupstoFit + 2*numGroupTimestoFit;
	

	if(incDMShapeEvent>0)DMdims+=incDMShapeEvent*(2+numDMShapeCoeff);
	if(incDMScatterShapeEvent>0)DMdims+=incDMScatterShapeEvent*(3+numDMScatterShapeCoeff);	
	if(incFloatDM>0)DMdims+=2*incFloatDM;

	if(incGlitch > 0)Glitchdims=(incGlitchTerms+1)*incGlitch;

    
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

	std::string GNname=GroupNoiseName;

	if(GNname.size() >= 100){printf("Group Noise flag is too long, needs to be less than 100 characters, currently %i .\n",GNname.size());return 0;}

	for(int r=0;r<=GNname.size();r++){GroupNoiseSysFlag[r]=GNname[r];}

	std::string ProfileFileName = ProfFile;
	if(ProfileFileName.size() >= 100){printf("Profile File Name is too long, needs to be less than 100 characters, currently %i .\n",ProfileFileName.size());return 0;}
	for(int r=0;r<=ProfileFileName.size();r++){ProfFileName[r]=ProfileFileName[r];}

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
	if(useGPUS == 1){
		printf("Using GPUs\n");
	}
#endif

#ifdef HAVE_MLAPACK
	printf("Have access to mlapack\n");
	if(useGPUS == 1){
		printf("Using GPUs\n");
	}
#endif


	printf("file root set to %s \n",root);
	
	int systemcount=0;
	int *numFlags=new int[psr[0].nobs];
	for(int o=0;o<psr[0].nobs;o++){
		numFlags[o]=0;
	}
	
	
	for(int o=0;o<psr[0].nobs;o++){
		psr[0].obsn[o].snr = 1;
		psr[0].obsn[o].tobs = 1;
		for (int f=0;f< psr[0].obsn[o].nFlags;f++){
                	if(strcasecmp( psr[0].obsn[o].flagID[f],"-snr")==0){
				psr[0].obsn[o].snr =  atof(psr[0].obsn[o].flagVal[f]);
			}
			if(strcasecmp( psr[0].obsn[o].flagID[f],"-tobs")==0){
				psr[0].obsn[o].tobs =  atof(psr[0].obsn[o].flagVal[f]);
			}
		}
			
	}

//	if(incEFAC > 0 || incEQUAD > 0 || incShannonJitter > 0){printf("using white noise model %i\n", whitemodel);}

	std::vector<std::string>systemnames;
	for(int o=0;o<psr[0].nobs;o++){
		int found=0;
		for (int f=0;f<psr[0].obsn[o].nFlags;f++){
			if(strcasecmp(psr[0].obsn[o].flagID[f],wflag)==0){
				
				if(std::find(systemnames.begin(), systemnames.end(), psr[0].obsn[o].flagVal[f]) != systemnames.end()) {
 				/* systemnames contains x */
				} else {
 	
 
	 				/* systemnames does not contain x */
					if(incEFAC==2 || incEQUAD==2 || incShannonJitter==2){printf("Found %s %s \n",wflag, psr[0].obsn[o].flagVal[f]);}
					systemnames.push_back(psr[0].obsn[o].flagVal[f]);
					systemcount++;
				}
				found=1;
			}
			

		}
		if(found==0 && (incEFAC==2 || incEQUAD==2 || incShannonJitter==2)){
			printf("Observation %i is missing the %s flag, please check before continuing\n",o,wflag);return 0;
		}
	}

	if(incEFAC==2 || incEQUAD==2 ||incShannonJitter==2 ){printf("total number of systems: %i \n",systemcount);}
	if(systemcount==0 && (incEFAC < 2 && incEQUAD < 2 && incShannonJitter < 2)){systemcount=1;}

		
	for(int o=0;o<psr[0].nobs;o++){
		for (int f=0;f<psr[0].obsn[o].nFlags;f++){
		
			if(strcasecmp(psr[0].obsn[o].flagID[f],wflag)==0){
				for (int l=0;l<systemcount;l++){
					if(psr[0].obsn[o].flagVal[f] == systemnames[l]){
						numFlags[o]=l;
					}
				}
			}

		}
	}

	int *includeEQsys=new int[systemcount];
	for (int l=0;l<systemcount;l++){
		includeEQsys[l] = 1;
	}


	if(incEQUAD==3){
		for(int i =0; i < psr[0].nTNEQ; i++){
			printf("Including for EQUAD==3 %i %s \n", i, psr->TNEQFlagVal[i]);
			for (int l=0;l<systemcount;l++){
				if(psr->TNEQFlagVal[i] == systemnames[l]){
					printf("This corresponds to system %i \n", l);
					includeEQsys[l] = 0;
				}
			}
		}

		for (int l=0;l<systemcount;l++){
			int changed = 0;
			if(includeEQsys[l] == 0){ includeEQsys[l] = 1; changed=1;}			
			if(includeEQsys[l] == 1 && changed == 0){ includeEQsys[l] = 0;}
		}

		for (int l=0;l<systemcount;l++){
			printf("Including? %i %i \n", l, includeEQsys[l]);
		}
		
	}



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////Get Group Noise System Positions/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int GroupNoiseSysCount=0;
	int *GroupNoiseSys=new int[psr[0].nobs];
	double **GroupStartTimes;
	for(int o=0;o<psr[0].nobs;o++){
		GroupNoiseSys[o]=0;
	}
	
	if(incGroupNoise > 0 || numStep > 0){
		std::vector<std::string>GroupNoiseSystemNames;
		for(int o=0;o<psr[0].nobs;o++){
			int found=0;
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
				if(strcasecmp(psr[0].obsn[o].flagID[f],GroupNoiseSysFlag)==0){
				
					if(std::find(GroupNoiseSystemNames.begin(), GroupNoiseSystemNames.end(), psr[0].obsn[o].flagVal[f]) != GroupNoiseSystemNames.end()) {
	 				/* systemnames contains x */
					} else {
	 	
	 
		 				/* systemnames does not contain x */
						printf("Found %s %s \n",GroupNoiseSysFlag, psr[0].obsn[o].flagVal[f]);
						GroupNoiseSystemNames.push_back(psr[0].obsn[o].flagVal[f]);
						GroupNoiseSysCount++;
					}
					found=1;
				}
			

			}

			if(found==0 ){
				printf("Observation %i is missing the %s flag, please check before continuing\n",o,GroupNoiseSysFlag);return 0;
			}
		}

		printf("total number of Group Noise Systems: %i \n",GroupNoiseSysCount);

		GroupStartTimes = new double*[GroupNoiseSysCount];
		for(int o = 0; o < GroupNoiseSysCount; o++){
			GroupStartTimes[o] = new double[2];

			GroupStartTimes[o][0] = 100000;
			GroupStartTimes[o][1] = 0;
		}
		
		for(int o=0;o<psr[0].nobs;o++){
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
		
				if(strcasecmp(psr[0].obsn[o].flagID[f],GroupNoiseSysFlag)==0){
					for (int l=0;l<GroupNoiseSysCount;l++){
						if(psr[0].obsn[o].flagVal[f] == GroupNoiseSystemNames[l]){
							GroupNoiseSys[o]=l;
							if((double)psr[0].obsn[o].sat < GroupStartTimes[l][0]){GroupStartTimes[l][0] = (double)psr[0].obsn[o].sat;}
							if((double)psr[0].obsn[o].sat > GroupStartTimes[l][1]){GroupStartTimes[l][1] = (double)psr[0].obsn[o].sat;}
						}
					}
				}

			}
		}

		for (int l=0;l<GroupNoiseSysCount;l++){
			printf("Group Times: %i %g %g \n", l, GroupStartTimes[l][0], GroupStartTimes[l][1]);
		}

	}



//	if(incEFAC == 0){printf("Not Including EFAC\n");numEFAC=0;}
	if(incEFAC == 1){printf("Including One EFAC for all observations\n");numEFAC=1;}
	if(incEFAC == 2){printf("Including One EFAC for each %s\n", wflag);numEFAC=systemcount;}	
//	if(incEQUAD == 0){printf("Not Including EQUAD\n");numEQUAD=0;}
	if(incEQUAD == 1){printf("Including One EQUAD for all observations\n");numEQUAD=1;}
	if(incEQUAD == 2){printf("Including One EQUAD for each %s\n", wflag);numEQUAD=systemcount;}
	if(incEQUAD == 3){printf("Including One EQUAD for the %i included TNEQ in par file\n",psr[0].nTNEQ);numEQUAD=psr[0].nTNEQ;}
	if(incHighFreqStoc > 0){printf("Including High frequency profile stochasticity \n");}
	if(incShannonJitter==1){printf("Including One SQUAD for all observations \n");numSQUAD=1;}
	if(incShannonJitter == 2){printf("Including One SQUAD for each %s\n", wflag);numSQUAD=systemcount;}
	if(incNGJitter == 1){printf("Including ECORR for %i systems\n", psr->nTNECORR);numNGJitter=psr->nTNECORR;}
	if(incDMEQUAD == 1){printf("Including DM EQUAD\n");}
	if(FitSolarWind == 1){printf("Including Deterministic Solar Wind\n");}
	if(FitWhiteSolarWind == 1){printf("Including Stochastic Solar Wind \n");}
//	if(incRED == 0){printf("Not Including Red Noise\n");}
	if(incRED == 1){printf("Including Red Noise : Power Law Model\n");}
	if(incRED == 2){printf("Including Red Noise : Model Independant - Fitting %i Coefficients\n", int(numRedCoeff));}
	if(incRED ==3){printf("Including Red Noise: Power Law Model to %i Coefficients \n", int(numRedCoeff));}
	if(incRED ==4){printf("Including Red Noise: Power Law Model with Corner Freq to %i Coefficients\n \n", int(numRedCoeff));}
	if(incRED ==5){printf("Including Red Noise Numerically: Power Law Model to %i Coefficients \n", int(numRedCoeff));}
//	if(incDM == 0){printf("Not Including DM\n");}
	if(incDM == 1){printf("Including DM : Power Law Model\n");}
	if(incDM == 2){printf("Including DM : Model Independant - Fitting %i Coefficients\n", int(numDMCoeff));}
	if(incDM ==3){printf("Including DM: %i Component Power Law Model to %i Coefficients \n", numDMPL, int(numDMCoeff));}
        if(incDM ==5){printf("Including DM Variations Numerically: Power Law Model to %i Coefficients \n", int(numDMCoeff));}
	if(incFloatDM==1){printf("Including Floating DM power spectrum coefficient\n");}
	if(incFloatRed==1){printf("Including Floating Red noise power spectrum coefficient\n");}
	if(yearlyDM==1){printf("Including yearly DM variations\n");}
	if(incsinusoid==1){printf("Including additional sinusoid\n");}	
	if(incDMEvent != 0){printf("Including %i DM Events\n",incDMEvent);}
	if(incDMShapeEvent != 0){printf("Including %i DM Shape Events\n",incDMShapeEvent);}
	if(incDMScatterShapeEvent != 0){printf("Including %i DM Scatter Shape Events\n",incDMShapeEvent);}
	if(incRedShapeEvent != 0){printf("Including %i Red Shape Events\n",incRedShapeEvent);}
	if(incBandNoise > 0){printf("Including Band Noise for %i bands\n", incBandNoise);}
	if(incGroupNoise > 0){printf("Including Group Noise: Power Law Model to %i Coefficients \n", int(numGroupCoeff));}
	if(incProfileFit > 0){printf("Including profile fit\n");}
	if(incProfileEvo > 0){printf("Including profile evolution with %i polynoimals\n", NProfileEvoPoly);}
	if(FitLinearProfileWidth == 1){printf("Including Linear Profile Width\n");}
	if(incWidthJitter == 1){printf("Including Width Jitter \n");}

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
		TempoPriors[i]=new long double[4];
		for(int j=0; j< 4; j++){
				TempoPriors[i][j]=0;
		}
		
	}
	TempoFitNums=new int*[fitcount];for(int i=0;i<fitcount;i++){TempoFitNums[i]=new int[2];}
	TempoJumpNums=new int[numFitJumps];
	//printf("allocated\n");
	int paramsfitted=0;
//	printf("Offset and error: %g %g %g\n", psr[0].offset, psr[0].offset_e, psr[0].offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree));
	TempoPriors[paramsfitted][0]=0;//psr[0].offset;
	TempoPriors[paramsfitted][1]=psr[0].offset_e;//offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree);

	if(GPTA == 1){
		TempoPriors[paramsfitted][0]=0;
        	TempoPriors[paramsfitted][1]=1;//offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree);
	}


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

//				printf("TempoPriors %i %Lg %Lg\n", paramsfitted-1, TempoPriors[paramsfitted-1][0], TempoPriors[paramsfitted-1][1]);
	    		}
		}
	}
	//printf("and time\n");
	int jumpsfitted=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1){

//			printf("gonna read jump %i %s \n",i,psr[0].jumpStr[i]);	
			char str1[100],str2[100],str3[100],str4[100],str5[100];
			int nread=sscanf(psr[0].jumpStr[i],"%s %s %s %s %s",str1,str2,str3,str4,str5);
			double prejump=atof(str3);
			//printf("Pre jump %i %g %g\n",i,prejump,psr[0].jumpVal[i]);
			//if(GPTA==0){	
			if(numTempo2its > 1){
				TempoPriors[paramsfitted][0]=psr[0].jumpVal[i];
				TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
				Tempo2Fit[paramsfitted]=psr[0].jumpVal[i];
			}
			else{
				TempoPriors[paramsfitted][0]=prejump;
				TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
				Tempo2Fit[paramsfitted]=prejump;
			}
			//}
			//else{

                        //        TempoPriors[paramsfitted][0]=prejump;
                        ///        TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
                       //         Tempo2Fit[paramsfitted]=prejump;

			//}
			TempoJumpNums[jumpsfitted]=i;
			
			if(doJumpMargin == 1)TempoPriors[paramsfitted][2]=1;
			
			paramsfitted++;
			jumpsfitted++;

		}
	}  

	// set the MultiNest sampling parameters
	
// 	return 0;
	int sampler = 1;
	int Nchords=1; 	
	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 500;				// number of live points
	double efr = 0.1;				// set the required efficiency
	int sample = 1;
	int nClsPar = 1;				// no. of parameters to do mode separation on
	int updInt = 2000;				// after how many iterations feedback is required & the output files should be updated

	setupMNparams(sampler, IS, mmodal, ceff, nlive, efr, sample, updInt, nClsPar, Nchords);



	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = numFitJumps+fitcount+numEFAC*numEPolTerms +numEQUAD+numNGJitter+Reddims+DMdims+ numStep*3 + DMModeldims + 2*varyDMCoeff + 2*varyRedCoeff + 2*yearlyDM + 3*incsinusoid+numSQUAD + 7*incDMEvent + Glitchdims +swdims + incHighFreqStoc + incDMEQUAD + incWidthJitter + pnoisedims; // dimensionality (no. of free parameters)
	int nPar = ndims;					// total no. of parameters including free & derived parameters
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

	MNStruct *MNS = init_struct(psr,TempoPriors,npsr,numFitJumps,fitcount,systemcount,numEFAC,numEQUAD, int(numRedCoeff), int(numDMCoeff), numRedPL, numDMPL, TempoFitNums,TempoJumpNums,numFlags, ndims, incRED,incDM, incFloatDM,incFloatRed, FloatDMstart, FloatRedstart, doTimeMargin,doJumpMargin, doLinearFit, SampleFreq, numStep, wflag, whitemodel,varyRedCoeff, varyDMCoeff,yearlyDM, incsinusoid, numEPolTerms, incGWB,RedPriorType, DMPriorType, EQUADPriorType,EFACPriorType,useOriginalErrors,numSQUAD, incDMEvent, incDMShapeEvent, numDMShapeCoeff, incBandNoise, numBandNoiseCoeff, incRedShapeEvent, numRedShapeCoeff, MarginRedShapeCoeff, incDMScatterShapeEvent, numDMScatterShapeCoeff, incNGJitter, incGlitch, incGlitchTerms, incBreakingIndex, FitLowFreqCutoff, uselongdouble, incGroupNoise, numGroupCoeff, FitForGroup, numGroupstoFit,GroupNoiseSys, FitSolarWind, FitWhiteSolarWind, interpolateProfile, InterpolatedTime, sampler, GPTAnumstocshapecoeff, totalshapestoccoeff, StoreFMatrices, incHighFreqStoc, numNGJitter, FitForBand, incProfileEvo, EvoRefFreq, numEvoFitCoeff, incWideBandNoise, incProfileFit, numProfileFitCoeff, incDMEQUAD, FitLinearProfileWidth, offPulseLevel,GroupStartTimes, FitEvoExponent,numProfComponents, totalEvoFitCoeff, totalProfileFitCoeff, numEvoCoeff, totalEvoCoeff,incWidthJitter, JitterProfComp, incProfileEnergyEvo, debug, ProfileBaselineTerms, incProfileNoise, ProfileNoiseCoeff, SubIntToFit, ChannelToFit, NProfileEvoPoly, ProfCompSeps,usecosiprior);
	
	MNS->includeEQsys = includeEQsys;	




/*
//	printf("Check NG: %i %i\n", incNGJitter, ndims );
	if(incNGJitter ==1){
		//printf("Calling\n");
		double **NGJitterMatrix;
		int NumNGEpochs;
		getNGJitterMatrixEpochs(psr, NumNGEpochs);

		NGJitterMatrix = new double*[psr->nobs];
		for(int i =0; i < psr->nobs; i++){
			NGJitterMatrix[i] = new double[NumNGEpochs];
			for(int j =0; j < NumNGEpochs;  j++){
				NGJitterMatrix[i][j] = 0;
			}
		}
		int *NGJitterSysFlags = new int[psr->nobs];
		

		for (int i=0;i<psr->nobs;i++){
			for (int j=0;j<psr->obsn[i].nFlags;j++){
				for (int k=0;k<psr->nTNECORR;k++){
					if (strcmp(psr->obsn[i].flagID[j], psr->TNECORRFlagID[k])==0){
						if (strcmp(psr->obsn[i].flagVal[j],psr->TNECORRFlagVal[k])==0){
							NGJitterSysFlags[i] = k;
							//printf("NGFLag? %i %i %s \n", i, k, psr->TNECORRFlagVal[k]);
						}
					}
				}
			}
		}


		getNGJitterMatrix(psr, NGJitterMatrix, NumNGEpochs);

		int *NGJitterEpochFlags = new int[NumNGEpochs];

		printf("Check: %i %g \n", NumNGEpochs, NGJitterMatrix[0][0]);
		for(int i =0; i < NumNGEpochs; i++){
			for(int j=0; j < psr->nobs; j++){
				if(NGJitterMatrix[j][i] != 0) {
				//	printf("%i %i %i\n", j, i, NGJitterSysFlags[j]);
					NGJitterEpochFlags[i]=NGJitterSysFlags[j];
				}
			}
		}
		delete[] NGJitterSysFlags;
		MNS->NGJitterSysFlags=NGJitterEpochFlags;
		MNS->NGJitterMatrix=NGJitterMatrix;
		MNS->numNGJitterEpochs=NumNGEpochs;
		MNS->incNGJitter=numNGJitter;
	}

*/
	
	//return 0;
	context=MNS;



#ifdef HAVE_CULA
	if(useGPUS==1){
		culaStatus status;
		status = culaInitialize();
		store_factorial();
	}
    
#endif /* HAVE_CULA */

	double **Dpriors;
	Dpriors = new double*[ndims]; for(int i = 0; i < ndims; i++){Dpriors[i]=new double[2];};

	int findbrakeparam=0;
	if(incBreakingIndex==1){
		
		for (int p=0;p<MAX_PARAMS;p++) {
		      for (int k=0;k<psr[0].param[p].aSize;k++){
				if(psr[0].param[p].fitFlag[k] == 1 && p != param_dmmodel){
					findbrakeparam++;			
					if(p == param_brake){


						printf("fitting for brake %i \n", p);

						TempoPriors[findbrakeparam][0]=0;
						TempoPriors[findbrakeparam][1]=1;
						TempoPriors[findbrakeparam][2]=0;
						TempoPriors[findbrakeparam][3]=1;
						Dpriors[findbrakeparam][0] = -3;
						Dpriors[findbrakeparam][1] = 3;


						Tempo2Fit[findbrakeparam]=0;

						TempoFitNums[findbrakeparam][0]=p;
						TempoFitNums[findbrakeparam][1]=k;

					}
				}
			}

		}
	}

	//If using custompriors for errors incase T2 doesnt converge, get those values before doing anything else
	if(customPriors == 1){
		setTNPriors(Dpriors, TempoPriors, ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
		update_MNPriors(MNS,Dpriors, TempoPriors,0);

		context=MNS;
		
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




	//Combine all the priors into one aray: Dpriors
	int pcount=0;

	for(int i =0; i< numFitJumps+fitcount; i++){
		Dpriors[pcount][0]=-FitSig;
		Dpriors[pcount][1]=FitSig;
		pcount++;
	}

	if(incBreakingIndex ==1){
		Dpriors[findbrakeparam][0]=-3;
		Dpriors[findbrakeparam][1]=3;
	}


	if(numStep>0){
		for(int i =0; i < numStep; i++){
			Dpriors[pcount][0]=0;
			Dpriors[pcount][1]=GroupNoiseSysCount;
			pcount++;
			Dpriors[pcount][0]=0;
			Dpriors[pcount][1]=1;
			pcount++;
			Dpriors[pcount][0]=-5*pow(10.0, -7);
			Dpriors[pcount][1]= 5*pow(10.0, -7);
			pcount++;

		}
	}			

	for(int i =0; i < incGlitch; i++){
		Dpriors[pcount][0]=start;
		Dpriors[pcount][1]=end;
		pcount++;
		for(int j=0; j < incGlitchTerms; j++){
			Dpriors[pcount][0]=-GlitchFitSig*psr[0].param[param_f].err[j];
			Dpriors[pcount][1]=GlitchFitSig*psr[0].param[param_f].err[j];
			pcount++;
		}
	}


	for(int i =0; i< numEFAC; i++){
		for(int n=1; n<=numEPolTerms; n++){
			if(n==1){
				Dpriors[pcount][0]=EFACPrior[0];
				Dpriors[pcount][1]=EFACPrior[1];
				pcount++;
			}
			else{
				Dpriors[pcount][0]=EPolPrior[0];
            Dpriors[pcount][1]=EPolPrior[1];
            pcount++;
			}
		}
	}
	for(int i =0; i< numEQUAD; i++){
		Dpriors[pcount][0]=EQUADPrior[0];
		Dpriors[pcount][1]=EQUADPrior[1];
		pcount++;
	}	
	
	for(int i =0; i< numSQUAD; i++){
		Dpriors[pcount][0]=EQUADPrior[0];
		Dpriors[pcount][1]=EQUADPrior[1];
		pcount++;
	}


	for(int i =0; i< incHighFreqStoc; i++){
		Dpriors[pcount][0]=HighFreqStocPrior[0];
		Dpriors[pcount][1]=HighFreqStocPrior[1];
		pcount++;
	}

	for(int i =0; i< numNGJitter; i++){
		Dpriors[pcount][0]=EQUADPrior[0];
		Dpriors[pcount][1]=EQUADPrior[1];
		pcount++;
	}

	if(incDMEQUAD == 1){
		Dpriors[pcount][0]=DMEQUADPrior[0];
		Dpriors[pcount][1]=DMEQUADPrior[1];
		pcount++;
	}

	if(incWidthJitter == 1){
		Dpriors[pcount][0]=WidthJitterPrior[0];
		Dpriors[pcount][1]=WidthJitterPrior[1];
		pcount++;
	}


	if(incProfileNoise == 1){
                Dpriors[pcount][0]=ProfileNoiseAmpPrior[0];
                Dpriors[pcount][1]=ProfileNoiseAmpPrior[1];
                pcount++;
                Dpriors[pcount][0]=ProfileNoiseSpecPrior[0];
                Dpriors[pcount][1]=ProfileNoiseSpecPrior[1];
                pcount++;
	}
	
	if(incProfileNoise == 2 || incProfileNoise == 3){
		for(int q = 0; q < ProfileNoiseCoeff; q++){
			Dpriors[pcount][0] = ProfileNoiseAmpPrior[0];
			Dpriors[pcount][1] = ProfileNoiseAmpPrior[1];
			pcount++;
		}
	}	

	if(FitSolarWind == 1){
                Dpriors[pcount][0]=SolarWindPrior[0];
                Dpriors[pcount][1]=SolarWindPrior[1];
                pcount++;
	}

	if(FitWhiteSolarWind == 1){
                Dpriors[pcount][0]=WhiteSolarWindPrior[0];
                Dpriors[pcount][1]=WhiteSolarWindPrior[1];
                pcount++;
	}
	
	if(incDMEvent > 0 ){
		for(int i =0;i < incDMEvent; i++){

                        Dpriors[pcount][0]=DMEventStartPrior[0];
                        Dpriors[pcount][1]=DMEventStartPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=DMEventLengthPrior[0];
                        Dpriors[pcount][1]=DMEventLengthPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=DMAmpPrior[0];
                        Dpriors[pcount][1]=DMAmpPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=DMAlphaPrior[0];
                        Dpriors[pcount][1]=DMAlphaPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=-pow(10.0,-2);
                        Dpriors[pcount][1]=pow(10.0,-2);
                        pcount++;
                        Dpriors[pcount][0]=-pow(10.0,-2);
                        Dpriors[pcount][1]=pow(10.0,-2);
                        pcount++;
                        Dpriors[pcount][0]=-pow(10.0,-2);
                        Dpriors[pcount][1]=pow(10.0,-2);
                        pcount++;
		}
	}

	if(FitLowFreqCutoff  > 0){
		Dpriors[pcount][0]=-1;
		Dpriors[pcount][1]=0;
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
		if(varyRedCoeff==1){
			Dpriors[pcount][0]=0;
			Dpriors[pcount][1]=log10(numRedCoeff);
			pcount++;
			Dpriors[pcount][0]=RedCoeffPrior[0];
			Dpriors[pcount][1]=RedCoeffPrior[1];
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
                for(int i =0; i< numRedPL; i++){

                        Dpriors[pcount][0]=AmpPrior[0];
                        Dpriors[pcount][1]=AmpPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=AlphaPrior[0];
                        Dpriors[pcount][1]=AlphaPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=-2;
                        Dpriors[pcount][1]=1;
                        pcount++;

                }
                if(varyRedCoeff==1){
                        Dpriors[pcount][0]=0;
                        Dpriors[pcount][1]=log10(numRedCoeff);
                        pcount++;
                        Dpriors[pcount][0]=RedCoeffPrior[0];
                        Dpriors[pcount][1]=RedCoeffPrior[1];
                        pcount++;
                }

	}
	else if(incRED==5){

		for(int i =0; i < 2*numRedCoeff;i++){
			//printf("Red C %i %g \n", i, psr[0].TNRedCoeffs[i+100]);
			Dpriors[pcount][0]=-500*pow(10.0,-9);
			Dpriors[pcount][1]=500*pow(10.0,-9);;
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
	if(incGWB==1){
		 Dpriors[pcount][0]=GWBAmpPrior[0];
	         Dpriors[pcount][1]=GWBAmpPrior[1];
	         pcount++;
	}



	if(incRedShapeEvent > 0 ){
		for(int i =0;i < incRedShapeEvent; i++){

	                Dpriors[pcount][0]=DMEventStartPrior[0];
	                Dpriors[pcount][1]=DMEventStartPrior[1];
	                pcount++;
	                Dpriors[pcount][0]=DMEventLengthPrior[0];
	                Dpriors[pcount][1]=DMEventLengthPrior[1];
		        pcount++;
			if(MarginRedShapeCoeff == 0){
				for(int c =0;c < numRedShapeCoeff; c++){
					Dpriors[pcount][0]=RedShapeCoeffPrior[0];
					Dpriors[pcount][1]=RedShapeCoeffPrior[1];
					pcount++;
				}
			}
		}
	}

	if(incsinusoid==1){
		Dpriors[pcount][0]=RedCoeffPrior[0];
		Dpriors[pcount][1]=RedCoeffPrior[1];
		pcount++;
		Dpriors[pcount][0]=0;
		Dpriors[pcount][1]=2*M_PI;
		pcount++;
		Dpriors[pcount][0]=-2;
		Dpriors[pcount][1]=log10(numRedCoeff);
       	 	pcount++;
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
		if(varyDMCoeff==1){
			Dpriors[pcount][0]=0;
			Dpriors[pcount][1]=log10(numDMCoeff);
			pcount++;
			Dpriors[pcount][0]=DMCoeffPrior[0];
			Dpriors[pcount][1]=DMCoeffPrior[1];
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

	else if(incDM==5){

	        for(int i =0; i < 2*numDMCoeff;i++){
	                //Dpriors[pcount][0]= -0.005;
	                //Dpriors[pcount][1]=  0.005;

	                Dpriors[pcount][0]= -10;
	                Dpriors[pcount][1]= 10;
	                pcount++;
	        }

	        Dpriors[pcount][0]=DMAmpPrior[0];
	        Dpriors[pcount][1]=DMAmpPrior[1];
	        pcount++;
	        Dpriors[pcount][0]=DMAlphaPrior[0];
	        Dpriors[pcount][1]=DMAlphaPrior[1];
	        pcount++;
	}

	else if(incDM==6){

	        for(int i =0; i < 2*numDMCoeff;i++){
	                //Dpriors[pcount][0]= -0.001;
	                //Dpriors[pcount][1]=  0.001;

	                Dpriors[pcount][0]= -10;
	                Dpriors[pcount][1]= 10;
	                pcount++;
	        }

	        for(int i =0; i < numDMCoeff;i++){
	                Dpriors[pcount][0]= DMCoeffPrior[0];
	                Dpriors[pcount][1]=  DMCoeffPrior[1];

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




	if(incDMShapeEvent > 0 ){
		for(int i =0;i < incDMShapeEvent; i++){

	                Dpriors[pcount][0]=DMEventStartPrior[0];
	                Dpriors[pcount][1]=DMEventStartPrior[1];
	                pcount++;
	                Dpriors[pcount][0]=DMEventLengthPrior[0];
	                Dpriors[pcount][1]=DMEventLengthPrior[1];
		        pcount++;
			for(int c =0;c < numDMShapeCoeff; c++){
			        Dpriors[pcount][0]=DMShapeCoeffPrior[0];
			        Dpriors[pcount][1]=DMShapeCoeffPrior[1];
			        pcount++;
			}
		}
	}

	if(incDMScatterShapeEvent > 0 ){
		for(int i =0;i < incDMScatterShapeEvent; i++){

		
	                Dpriors[pcount][0]=DMEventStartPrior[0];
	                Dpriors[pcount][1]=DMEventStartPrior[1];
	                pcount++;
	                Dpriors[pcount][0]=DMEventLengthPrior[0];
	                Dpriors[pcount][1]=DMEventLengthPrior[1];
		        pcount++;
	                //Dpriors[pcount][0]=DMScatterFreqPrior[0];
	                //Dpriors[pcount][1]=DMScatterFreqPrior[1];
		        pcount++;
			for(int c =0;c < numDMScatterShapeCoeff; c++){
				printf("P %i %i %i %i\n", i,c, pcount, ndims);
			        Dpriors[pcount][0]=DMScatterShapeCoeffPrior[0];
			        Dpriors[pcount][1]=DMScatterShapeCoeffPrior[1];
			        pcount++;
			}
		}
	}

	if(yearlyDM==1){
		Dpriors[pcount][0]= -0.01;//DMCoeffPrior[0];
		Dpriors[pcount][1]=  0.01;//DMCoeffPrior[1];
		pcount++;
		Dpriors[pcount][0]=-0.01;//0;
		Dpriors[pcount][1]=0.01;//2*M_PI;
		pcount++;
	}

	printf("about to do priors for band noise %i %i\n", pcount, ndims);
	if(incBandNoise > 0){
		for(int i =0; i < incBandNoise; i++){
			printf("priors for band noise %i %i %i\n", i, pcount, ndims);
			Dpriors[pcount][0]=BandNoiseAmpPrior[0];
			Dpriors[pcount][1]=BandNoiseAmpPrior[1];
			pcount++;
			Dpriors[pcount][0]=BandNoiseAlphaPrior[0];
			Dpriors[pcount][1]=BandNoiseAlphaPrior[1];
			pcount++;
		}
	}



	if(incGroupNoise > 0){
		for(int i =0; i < incGroupNoise; i++){

			if(FitForGroup[i][0]==-1){
				Dpriors[pcount][0]=0;
				Dpriors[pcount][1]=GroupNoiseSysCount;
				pcount++;
			}
	                if(FitForGroup[i][1]== 1){
	                        Dpriors[pcount][0]=0;
	                        Dpriors[pcount][1]=1;
	                        pcount++;
	                        Dpriors[pcount][0]=0;
	                        Dpriors[pcount][1]=1;
	                        pcount++;

	                }
		


			Dpriors[pcount][0]=GroupNoiseAmpPrior[0];
			Dpriors[pcount][1]=GroupNoiseAmpPrior[1];
			pcount++;
			Dpriors[pcount][0]=GroupNoiseAlphaPrior[0];
			Dpriors[pcount][1]=GroupNoiseAlphaPrior[1];
			pcount++;	
		}
	}	

/*
	numToMargin=0;
	for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
		if(TempoPriors[i][2]==1){
			numToMargin++;
		}
	}
*/
	if(GPTA == 0){




		int *FitList=new int[ndims];
		for(int i=0;i<ndims;i++){
			FitList[i]=0;
		}
		numToMargin=0;
		for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
			if(TempoPriors[i][2]==1){
				FitList[i]=1;
				printf("marginalising over timing param: %i \n",i);
				numToMargin++;
			}
		}

		getCustomDMatrix(psr, FitList, TempoFitNums, TempoJumpNums, Dpriors, incDM, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps);

		update_MNPriors(MNS,Dpriors, TempoPriors,0);

		context=MNS;

		
		//Finally after doing everything reset custom priors in case overwritten by previous steps.
		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors,  ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
			paramsfitted=0;
		}
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
		

		printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff, numDMCoeff, incFloatRed,incFloatDM, fitDMModel, longname, numStep, varyRedCoeff, varyDMCoeff,yearlyDM, incsinusoid, numEPolTerms,numSQUAD,context);

		printf("\n\n");
		ndims=ndims-numToMargin;
		nPar=ndims;


		
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////get TotalMatrix////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////   



		getArraySizeInfo(context);
		printf("Pre-Computing Matrices\n");
	
	
		for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
	
			((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = ((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);

		}

		formBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);      
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);      

		double *TotalMatrix = new double[((MNStruct *)context)->pulse->nobs*((MNStruct *)context)->totalsize]();	

		StoreTMatrix(TotalMatrix, context);
	
#ifdef HAVE_CULA	
		if(useGPUS==1){
			copy_staticTmat_(TotalMatrix, ((MNStruct *)context)->totalsize, ((MNStruct *)context)->pulse->nobs);
			delete[] TotalMatrix;
		}
#endif



		((MNStruct *)context)->StoredTMatrix = TotalMatrix;

//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////call PolyChord/////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

		double *PriorsArray=new double[2*ndims];
		int NDerived = 0;



		int p = 0;
		pcount = 0;
		while(pcount < ndims){
			if(((MNStruct *)context)->Dpriors[p][0] != ((MNStruct *)context)->Dpriors[p][1]){
				PriorsArray[pcount]=((MNStruct *)context)->Dpriors[p][0];
				PriorsArray[pcount+ndims]=((MNStruct *)context)->Dpriors[p][1];
				printf("param %i in, %g %g \n", pcount, PriorsArray[pcount], PriorsArray[pcount+ndims]); 
				pcount++;
			}
	
			p++;
		}

		double *output = new double[5];

		/////////////////////////Sort Grade stuff////////////////////////////   


		int do_grades = 0;
		int *grades = new int[ndims];
		for(int i = 0; i < ndims; i++){
			grades[i] = 1;
		}

		int maxgrade = 1;


		int *hypercube_indices = new int[ndims];
		int *physical_indices = new int[ndims];

		//Set indices, note this is fortran convention -> starts from 1/////

		for(int i = 0; i < ndims; i++){
	        	hypercube_indices[i] = i+1;
		        physical_indices[i] = i+1;
		}


	       int numfast = 0;
	       int numslow = 0;

	       for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
			if(((MNStruct *)context)->LDpriors[i][2]==0){
				numfast++;
			}
	       }

		numslow = ndims - numfast;

		//if(numfast > 0){
		//	do_grades = 1;
		//	maxgrade = 2;
		//}

  	        int *grade_repeats = new int[maxgrade];
	        for(int i = 0; i < maxgrade; i++){
	                grade_repeats[i] = 1;
	        }




		if(do_grades == 1){


			int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
			int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);
			int FitBandNoiseCoeff=2*(((MNStruct *)context)->numFitBandNoiseCoeff);
			int FitGroupNoiseCoeff = 2*((MNStruct *)context)->numFitGroupNoiseCoeff;


			int totCoeff=0;
			if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
			if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;
			if(((MNStruct *)context)->incBandNoise > 0) totCoeff += ((MNStruct *)context)->incBandNoise*FitBandNoiseCoeff;
			if(((MNStruct *)context)->incNGJitter >0)totCoeff+=((MNStruct *)context)->numNGJitterEpochs;
			if(((MNStruct *)context)->incGroupNoise > 0)totCoeff += ((MNStruct *)context)->incGroupNoise*FitGroupNoiseCoeff;


			((MNStruct *)context)->totCoeff = totCoeff;

			int TimetoMargin=0;
			for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
				if(((MNStruct *)context)->LDpriors[i][2]==1)TimetoMargin++;
			}

			int totalsize = totCoeff + TimetoMargin;


			double *LastParams = new double[ndims];
			for(int i = 0; i < ndims; i++){ LastParams[i] = 0;}

			double **PreviousTNT = new double*[totalsize];
			for(int i = 0; i < totalsize; i++){
				PreviousTNT[i] = new double[totalsize];
				for(int j = 0; j < totalsize; j++){
					PreviousTNT[i][j] = 0;
				}
			}

			double **PreviousNT = new double*[((MNStruct *)context)->pulse->nobs];
			for(int i = 0; i < ((MNStruct *)context)->pulse->nobs; i++){
				PreviousNT[i] = new double[totalsize];
				for(int j = 0; j < totalsize; j++){
					PreviousNT[i][j] = 0;
				}
			}

			double *PreviousNoise = new double[((MNStruct *)context)->pulse->nobs];
			for(int i = 0; i < ((MNStruct *)context)->pulse->nobs; i++){
				PreviousNoise[i] = 0;
			}


			((MNStruct *)context)->LastParams = LastParams;
			((MNStruct *)context)->PreviousTNT = PreviousTNT;
			((MNStruct *)context)->PreviousNT = PreviousNT;
			((MNStruct *)context)->PreviousNoise = PreviousNoise;
			((MNStruct *)context)->PreviousInfo = 0;

			pcount = 0;
			for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
	        		if(((MNStruct *)context)->LDpriors[i][2]==0){
					grades[numslow + pcount] = 2;
					pcount++;
				}
			}


			grade_repeats[1] = 10;

			for(int i = 0; i < numslow; i++){
				hypercube_indices[i] = 1 + numfast + i;
			}
			for(int i = 0; i < numfast; i++){
				hypercube_indices[numslow+i] = 1 + i;
			}


		}
		((MNStruct *)context)->PriorsArray = PriorsArray;
		((MNStruct *)context)->hypercube_indices = hypercube_indices;	
		((MNStruct *)context)->doGrades = do_grades;
		((MNStruct *)context)->PolyChordGrades = grades;	

		assigncontext(context);

#ifdef HAVE_CULA
		assignGPUcontext(context);
#endif



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////Call Samplers////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		if(sample==1){

#ifdef HAVE_CULA
			if(useGPUS==1){
				if(sampler == 0){
					 nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedGPULikeMNWrap, dumper, context);
				}
				if(sampler == 1){
					chord::Sample(NewLRedMarginGPULogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
				}

			}
			else{
				if(sampler == 0){
					 nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLikeMNWrap, dumper, context);
				}
				if(sampler == 1){
					chord::Sample(NewLRedMarginLogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
				}
	      		}   
#else
			if(sampler == 0){
				nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLikeMNWrap, dumper, context);
                        }
                        if(sampler == 1){
				chord::Sample(NewLRedMarginLogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
			}
#endif 
		}
	}
	else if(GPTA == 1){

#ifdef HAVE_PSRCHIVE


		//If not marginalising, only need to reget custom priors to overwrite anything done previously
		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors,  ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps,ndims);
		}

		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
	
		printPriors(psr, TempoPriors, Dpriors, numEFAC, numEQUAD, incRED, incDM, numRedCoeff, numDMCoeff, incFloatRed,incFloatDM,fitDMModel, longname, numStep, varyRedCoeff, varyDMCoeff,yearlyDM, incsinusoid, numEPolTerms,numSQUAD,context);


		printf("\n\n");
		

		//if(GPTA == 1){


		long double ***ProfileData = new long double**[((MNStruct *)context)->pulse->nobs];
		long double **ProfileInfo = new long double*[((MNStruct *)context)->pulse->nobs];
		std::vector<int> NumChanPerInt;
		//double weightsum=0;
		//double weightedfreq = 0;
		int t=0;
		//for(int t = 0; t < ((MNStruct *)context)->pulse->nobs; t++){
		while(t < ((MNStruct *)context)->pulse->nobs){
			int FoundProfile = 0;
			std::string ProfileName =  ((MNStruct *)context)->pulse->obsn[t].fname;

			if(debug == 1){printf("loading profile %i %s \n", t, ((MNStruct *)context)->pulse->obsn[t].fname);}

			//std::string fname = "profiles/"+ProfileName+".ASCII"; //"/home/ltl21/scratch/Pulsars/ProfileData/J1909-10cm/profiles/"+ProfileName+".ASCII";

			//std::string archivename = "profiles/"+ProfileName+".new_code"; //"/home/ltl21/scratch/Pulsars/ProfileData/J1909-10cm/profiles/"+ProfileName+".ASCII";
			Reference::To< Archive > archive = Archive::load(ProfileName.c_str());

			///std::string NewArchiveName = "NewArchive.rf.pcm.dzTf32";						
			//archive->unload(NewArchiveName);
			//SaveNewArchive(archive, NewArchiveName);


			if( !archive ){
				printf("Problem opening archive\n");			
				 return 0;	

			}
			int npols = 0;
			npols = archive->get_npol();
			if(npols > 1){printf("PScrunching Archive to get total Intensity\n"); archive->pscrunch();}
			//archive->remove_baseline();

			//SetNewAmps(archive);	
	
			//return 0;

			npols = archive->get_npol();
			int nsub = archive->get_nsubint();
			int nbins = 0;
			int nchans = 0;
			
			double foldingperiod = 0;
			double inttime = 0;
			double centerfreq = 0;
			double chanfreq = 0;
			int intday = 0;
			double fracday = 0;
			int intsec = 0;
			double fracsecs = 0;
			bool isdedispersed = 0;


			if(debug == 1){printf("Opened: %s, nsub = %i, npol = %i\n", ProfileName.c_str(), nsub, npols);}
			//std::string telescope = archive->get_telescope();
			//Tempo2::Observatory *obs = new Tempo2::Observatory(telescope);
			//printf("telescope code: %s \n", telescope.c_str());
			if(nsub > 1){
				//printf("Only one subint per archive supported currently: nsub %i \n", nsub);
				//return 0;
			}

			if(FitTemplate == 1){printf("Opening Template to fit: %i subs \n", nsub);}
			float *pvalues;
			for(int i=0; i < nsub; i++){

				if(FoundProfile == 1){break;}
				Integration *subint = archive->get_Integration(i);

				nbins = subint->get_nbin();
				nchans = subint->get_nchan();
				npols = subint->get_npol();
				foldingperiod = subint->get_folding_period();
				inttime = subint->get_duration();
				centerfreq = subint->get_centre_frequency(); 
				MJD firstbin = subint->get_epoch();
				intday = firstbin.intday();
				fracday = firstbin.fracday();
				intsec = firstbin.get_secs();
				fracsecs = firstbin.get_fracsec();
				isdedispersed = subint->get_dedispersed();




				if(FitTemplate == 1){printf("Template Details: %i subs %i chans %i pols \n", nsub, nchans, npols);}

                                long double ProfileMJD = intday;    //ds[0];
                                long double FirstBinSec = intsec+(long double)fracsecs;   //ds[1];
                                long double FoldingPeriod = foldingperiod;  //ds[2];
                                int Nbins = nbins;  //(int)ds[6];

	
			
				long double SubIntTimeDiff = FirstBinSec-((MNStruct *)context)->pulse->obsn[t].sat_sec*SECDAY;
				long double PeriodDiff = SubIntTimeDiff*((MNStruct *)context)->pulse->param[param_f].val[0];
				double DPeriodDiff = (double)PeriodDiff;			
				int usedchan=0;
				//printf("in profile %i of %i, channels: %i with %i pol %i subs %g\n", t, ((MNStruct *)context)->pulse->nobs, nchans, npols, nsub, DPeriodDiff);

				if(FitTemplate == 1){printf("Template Details: MJD %.15Lg %.15Lg DPdiff: %.10g \n", ProfileMJD, FirstBinSec/SECDAY, DPeriodDiff);}


				if(debug == 1){
					printf("Details:\n");
					printf("%i subs %i chans %i pols %i DD \n", nsub, nchans, npols, isdedispersed);
					printf("TOA: %.15Lg %.15Lg\n", ((MNStruct *)context)->pulse->obsn[t].sat_day, ((MNStruct *)context)->pulse->obsn[t].sat_sec);
					printf("MJD %.15Lg %.15Lg DPdiff: %.10g \n", ProfileMJD, FirstBinSec/SECDAY, DPeriodDiff);	
				}



				if(fabs(DPeriodDiff) < 2.0){
					
					for(int j = 0; j < nchans; j++){
				//		printf("in profile %i of %i, channel: %i of %i with %i pol\n", t, ((MNStruct *)context)->pulse->nobs, j, nchans, npols);

						int ipol = 0; //only want total Intensity
						Profile *prof = subint->get_Profile(ipol, j);
						
						//prof->remove_baseline();

						pvalues = prof->get_amps();
						chanfreq = subint->get_centre_frequency(j);
						double toafreq = ((MNStruct *)context)->pulse->obsn[t].freq;
						double checkfreq = 0;
	//					if(isdedispersed == 0){
							checkfreq = chanfreq;				     
	//					}
	//					else{
	//						checkfreq = centerfreq;
	//					}
	 
						double oneflux = 0;	
						for(int b =0; b < Nbins; b++){
							oneflux = oneflux + pvalues[b];
						}   

						if(FitTemplate == 1){printf("Template Details: %.10g Freq %.10g flux \n", chanfreq, oneflux);}


						if(debug == 1){
							printf("Chan Details:\n");
							printf("%.10g Freq %.10g flux \n", chanfreq, oneflux);
							printf("freq toa %.10g prof %.10g flux %g\n", toafreq, checkfreq, oneflux);	
						}


						//printf("checks: freq toa %.10g prof %.10g flux %g\n", toafreq, checkfreq, oneflux);
						if(oneflux == 0){
			//				printf("Profile %i: subint %i, channel %i is all zeros, not including it in the analysis\n", t, i, j);
						}
						else{
							if(fabs(toafreq - checkfreq) < 0.001){
								if(debug ==1){printf("Including: ToA %i SubInt %i Chan %i details: %.20Lg %.20Lg %.10g\n", t, i, j, ((MNStruct *)context)->pulse->obsn[t].sat_day, ((MNStruct *)context)->pulse->obsn[t].sat_sec, toafreq);}

								ProfileData[t] = new long double*[Nbins];
								for(int b = 0; b < Nbins; b++){
									ProfileData[t][b] = new long double[2];
								}


								ProfileInfo[t] = new long double[7];

								ProfileInfo[t][0] = FoldingPeriod;
								ProfileInfo[t][1] = (long double)nbins;     //ds[6];
								
								double minval=pow(10.0, 100);
								for(int b =0; b < Nbins; b++){
									if(pvalues[b] < minval){ minval=pvalues[b];}
								}
								for(int b =0; b < Nbins; b++){
									ProfileData[t][b][1] = pvalues[b];    //prof[1];
									//printf("%i %i %g \n", t, b, pvalues[b]);
								}   


								double Tobs = inttime;

								double PNoiseVal = 0;

								ProfileInfo[t][2] = (long double) Tobs;
								ProfileInfo[t][3] = (long double) PNoiseVal;
								ProfileInfo[t][4] = ((MNStruct *)context)->pulse->obsn[t].sat_day;
								ProfileInfo[t][5] = ((MNStruct *)context)->pulse->obsn[t].sat_sec;
								ProfileInfo[t][6] = ProfileMJD;


								long double pulsesamplerate = FoldingPeriod/Nbins/SECDAY;
								long double pulsesamplerateSec = FoldingPeriod/Nbins;
								long double *ProfileSATS = new long double[Nbins];

								if(FirstBinSec/SECDAY > 1){ printf("GREATER THAN A DAY! \n");}

								for(int b =0; b < Nbins; b++){
									ProfileData[t][b][0] = FirstBinSec/SECDAY + pulsesamplerate*b + pulsesamplerate*0.5;
								}

								t=t+1;
								usedchan++;
								if(t ==  ((MNStruct *)context)->pulse->nobs){printf("got them all \n"); FoundProfile = 1; break;}
							}
						}
					}
					if(usedchan > 0){
						//printf("num used chan is %i \n", usedchan);
						NumChanPerInt.push_back(usedchan);
					}
				}
				
				
			}
		}

		//printf("weighted average freq: %g \n", weightedfreq/weightsum);
		printf("Loaded  %i  Profiles, %i Epochs\n", t, NumChanPerInt.size());
		//return 0;

		int *numchanperint = new int[NumChanPerInt.size()];
		for(int p = 0; p < NumChanPerInt.size(); p++){
			numchanperint[p] = NumChanPerInt[p];
		}
	
		for(int p = 0; p < ((MNStruct *)context)->pulse->nobs; p++){
			((MNStruct *)context)->pulse->obsn[p].origsat = ((MNStruct *)context)->pulse->obsn[p].sat;
		}

		((MNStruct *)context)->ProfileData = ProfileData;
		((MNStruct *)context)->ProfileInfo = ProfileInfo;
		((MNStruct *)context)->ReferencePeriod = ProfileInfo[0][0];
		((MNStruct *)context)->numChanPerInt = numchanperint;
		((MNStruct *)context)->numProfileEpochs = NumChanPerInt.size();

		if(SubIntToFit >= ((MNStruct *)context)->numProfileEpochs){ printf("Not that many subints %i %i \n", SubIntToFit, ((MNStruct *)context)->numProfileEpochs); return 0;}

		int totshapecoeff = 0;	
		for(int i = 0; i < numProfComponents; i++){
			totshapecoeff += GPTAnumshapecoeff[i];	
		}
		//int nprofcoeff = GPTAnumshapecoeff;
		//int numshapestoccoeff = totalshapestoccoeff;

		int shapedims = totshapecoeff + totalshapestoccoeff;

		int maxshapecoeff = 0;
		if(totshapecoeff+1>=totalshapestoccoeff+1){
			maxshapecoeff=totshapecoeff+1;
		}
		if(totalshapestoccoeff+1 > totshapecoeff+1){
			maxshapecoeff=totalshapestoccoeff+1;
		}


		//int nprofcoeff = 1;
		((MNStruct *)context)->numshapecoeff = GPTAnumshapecoeff;
		((MNStruct *)context)->totshapecoeff = totshapecoeff;
		//((MNStruct *)context)->numshapestoccoeff = numshapestoccoeff;

		double *Factorials=new double[maxshapecoeff];
		double *Binomials=new double[maxshapecoeff];
		for(int i = 0; i < maxshapecoeff; i++){
			Factorials[i] = iter_factorial(i);
//			printf("%i %g \n", i, Factorials[i]);
		}
		for(int i =0; i < maxshapecoeff; i=i+2){
			double bsum=1;
			for(int j = 1; j <= i/2; j++){
				bsum *= double(i + 1 - j)/j;
//				printf("%i %i %g \n", i, j, bsum);
			}
			Binomials[i] = bsum;
		}
		((MNStruct *)context)->Factorials = Factorials;
		((MNStruct *)context)->Binomial = Binomials;

		int profevodims = 0;
		if(incProfileEvo == 0){
			profevodims = 0; 
			totalEvoFitCoeff = 0;
			for(int i =0; i < numProfComponents; i++){
				numEvoFitCoeff[i] = 0;
			}
		}
		if(incProfileEvo == 1){profevodims = totalEvoFitCoeff*NProfileEvoPoly;}
		if(incProfileEvo == 2){profevodims = 1*NProfileEvoPoly;}
		if(FitEvoExponent == 1 && incProfileEvo > 0){profevodims++;}


		int profEnergyEvoDims=0;
		if(incProfileEnergyEvo == 2){profEnergyEvoDims = 1;}

		int proffitdims = 0;
		if(incProfileFit == 0){
			proffitdims = 0; 
			totalProfileFitCoeff = 0;
			for(int i =0; i < numProfComponents; i++){
                                numProfileFitCoeff[i] = 0;
                        }
		}
		if(incProfileFit == 1){proffitdims = totalProfileFitCoeff;}


		int olddims=numFitJumps+fitcount+numEFAC+numEQUAD+numSQUAD+Reddims+DMdims+incHighFreqStoc+incDMEQUAD+incWidthJitter + pnoisedims + 2*yearlyDM;
		ndims = totshapecoeff+1+numFitJumps+fitcount+numEFAC+numEQUAD+numSQUAD+totalshapestoccoeff+Reddims+DMdims+incHighFreqStoc+profevodims+proffitdims+incDMEQUAD+FitLinearProfileWidth+incWidthJitter+profEnergyEvoDims + pnoisedims  + 2*yearlyDM;
		nPar = 0;


		if(FixProfile == 1){ndims -= (totshapecoeff+1);}
		if(FitTemplate > 0){ndims = 2; shapedims=0;}
		if(FitTemplate > 0 && numProfComponents > 1){ndims+=numProfComponents-1;}


		printf("Starting GPTA %i %i %i\n", ndims, totshapecoeff, totalshapestoccoeff);
		double *PriorsArray=new double[2*ndims];
		int NDerived = 0;
		//assigncontext(context);


		pcount = 0;
		for(int p = 0; p < olddims; p++){

		        printf("Prior %i %.10g %.10g \n", p,Dpriors[p][0], Dpriors[p][1]);

	                PriorsArray[p]=((MNStruct *)context)->Dpriors[p][0];
	                PriorsArray[p+ndims]=((MNStruct *)context)->Dpriors[p][1];
			pcount++;

		}


		double *output = new double[5];

		/////////////////////////Sort Grade stuff////////////////////////////   


		int do_grades = 0;
		int *grades = new int[ndims];
		for(int i = 0; i < ndims; i++){
			grades[i] = 1;
		}

		int maxgrade = 1;
		int *grade_repeats = new int[maxgrade];
		for(int i = 0; i < maxgrade; i++){
			grade_repeats[i] = 1;
		}


		int *hypercube_indices = new int[ndims];
		int *physical_indices = new int[ndims];

		//Set indices, note this is fortran convention -> starts from 1/////

		for(int i = 0; i < ndims; i++){
			hypercube_indices[i] = i+1;
			physical_indices[i] = i+1;
		}

	

		double **shapecoeffprior = new double*[shapedims];
		for(int i=0; i < shapedims; i++){
			shapecoeffprior[i] = new double[2];
			for(int j=0; j < 2; j++){
				shapecoeffprior[i][j] = 0.0;
			}
		}

		double betaminp = 0;
		double betamaxp = 0;
		double *BetaPrior = new double[2];
		double MeanBeta=0;

		int ReadProfFile = 1;
		if(FitTemplate > 0){ReadProfFile = 0;}


		double *MeanProfileShape;
		double **MeanProfileEvo;
		double *MeanProfileStoc;

		if(ReadProfFile == 1){
			std::ifstream MeanProfileFile;
			std::string MeanProfileFilename = ProfileFileName;
			MeanProfileFile.open(MeanProfileFilename.c_str());

			std::string line;
		        getline(MeanProfileFile,line);
		        std::istringstream myStream( line );
		        std::istream_iterator< double > begin(myStream),eof;
		        std::vector<double> ProfileVal(begin,eof);

			int numMeanCoeff = ProfileVal[0];
			int numMeanEvoCoeff = ProfileVal[1];
			//printf("Read File: num coeffs = %i %i\n", numMeanCoeff,numEvoCoeff);

			MeanProfileShape  = new double[numMeanCoeff];
			MeanProfileEvo    = new double*[NProfileEvoPoly];
			for(int i = 0; i < NProfileEvoPoly; i++){	
				MeanProfileEvo[i] = new double[numMeanEvoCoeff];
			}
			MeanProfileStoc = new double[totalshapestoccoeff];
			for(int i = 0; i < numMeanCoeff; i++){
				std::string line;
				getline(MeanProfileFile,line);
				std::istringstream myStream( line );
				std::istream_iterator< double > begin(myStream),eof;
				std::vector<double> ProfileVal(begin,eof);
				MeanProfileShape[i] = ProfileVal[0];
				//printf("Read File: coeff %i = %g \n", i, MeanProfileShape[i]);
				
			}

			for(int i = 0; i < 1; i++){
				getline(MeanProfileFile,line);
				std::istringstream myStream( line );
				std::istream_iterator< double > begin(myStream),eof;
				std::vector<double> ProfileVal(begin,eof);
				MeanBeta = ProfileVal[0];
				//printf("Read File: Beta = %g \n", MeanBeta);
				
			}
			for(int p = 0; p < NProfileEvoPoly; p++){		
				for(int i = 0; i < numMeanEvoCoeff; i++){
					std::string line;
					getline(MeanProfileFile,line);
					std::istringstream myStream( line );
					std::istream_iterator< double > begin(myStream),eof;
					std::vector<double> ProfileVal(begin,eof);
					MeanProfileEvo[p][i] = ProfileVal[0];
					printf("Read File: Evo coeff %i = %g \n", i, MeanProfileEvo[p][i]);

				}
			}
			for(int i = 0; i < totalshapestoccoeff; i++){
				std::string line;
				getline(MeanProfileFile,line);
				std::istringstream myStream( line );
				std::istream_iterator< double > begin(myStream),eof;
				std::vector<double> ProfileVal(begin,eof);
				MeanProfileStoc[i] = ProfileVal[0];
				printf("Read File: Stoc coeff %i = %g \n", i, MeanProfileStoc[i]);
				
			}
			//return 0;
							
		}

		setShapePriors(shapecoeffprior, BetaPrior, shapedims);




		for(int p = 0; p < totshapecoeff; p++){

			if(FixProfile == 0 && FitTemplate == 0){
				printf("Shape Prior %i %.10g %.10g \n", pcount,shapecoeffprior[p][0], shapecoeffprior[p][1]);
				PriorsArray[pcount]=shapecoeffprior[p][0];
				PriorsArray[pcount+ndims]=shapecoeffprior[p][1];
				pcount++;
			}
			//else if(FixProfile == 1 && FitTemplate == 0){
			//	MeanProfileShape[p] = shapecoeffprior[p][0];
			//}

		}
		if(FixProfile == 0 || FitTemplate > 0){
			printf("Beta Prior: %i %i %g %g\n", totshapecoeff, pcount, BetaPrior[0], BetaPrior[1]);
			PriorsArray[pcount] = BetaPrior[0];
			PriorsArray[pcount+ndims] = BetaPrior[1];	
			pcount++;
		}
		//else{
			//MeanBeta = BetaPrior[0];
		//}

		if(FitTemplate > 0 && numProfComponents > 1){
			
			printf("Additional component \n");
			PriorsArray[pcount] =0.3;
			PriorsArray[pcount+ndims] = 0.7;	
			pcount++;

			printf("Additional component \n");
			PriorsArray[pcount] = -0.25;
			PriorsArray[pcount+ndims] = -0.15;	
			pcount++;
		}
			
		for(int p = 0; p < totalshapestoccoeff; p++){
			printf("p %i %i %i\n", p, pcount, ndims);
			printf("Shape Power Prior %i %i %i %.10g %.10g \n", p,pcount, totalshapestoccoeff, shapecoeffprior[totshapecoeff+p][0], shapecoeffprior[totshapecoeff+p][1]);

                        PriorsArray[pcount]=-2;//shapecoeffprior[totshapecoeff+p][0];
                        PriorsArray[pcount+ndims]=2;//shapecoeffprior[totshapecoeff+p][1];

			//if(p==0){
	                    //    PriorsArray[pcount] = -10;
                       		//PriorsArray[pcount+ndims] = -10;	
			//}
                        pcount++;

                }
 
		int psum=0;
		for(int p2 = 0; p2 < numProfComponents; p2++){
			for(int p3 = 0; p3 < numProfileFitCoeff[p2]; p3++){
				printf("fit %i %i\n", p2, p3);
				printf("fit Prior %i %i %.10g %.10g %.10g %.10g\n", pcount, proffitdims, ProfileFitPrior[0], ProfileFitPrior[1], shapecoeffprior[psum+p3][0], shapecoeffprior[psum+p3][1]);

		                PriorsArray[pcount] = ProfileFitPrior[0];
		                PriorsArray[pcount+ndims] = ProfileFitPrior[1];

				if(shapecoeffprior[psum+p3][0] != 0 || shapecoeffprior[psum+p3][1] != 0){
			                PriorsArray[pcount] = shapecoeffprior[psum+p3][0];
		     		        PriorsArray[pcount+ndims] = shapecoeffprior[psum+p3][1];	
					printf("change fit Prior %i %i %.10g %.10g %.10g %.10g\n", pcount, proffitdims, ProfileFitPrior[0], ProfileFitPrior[1], shapecoeffprior[psum+p3][0], shapecoeffprior[psum+p3][1]);
				}		

		                pcount++;
			}
			psum += GPTAnumshapecoeff[p2];

                }

		if(FitLinearProfileWidth == 1){
			printf("fit Linear Width Prior %.10g %.10g \n", LinearProfileWidthPrior[0], LinearProfileWidthPrior[1]);

                        PriorsArray[pcount] = LinearProfileWidthPrior[0];
                        PriorsArray[pcount+ndims] = LinearProfileWidthPrior[1];
			pcount++;
		}

		int evopcount = 0;
		if(FitEvoExponent == 1){
			printf("Fit Evo Exponent\n");
                        PriorsArray[pcount] = 0;
     		        PriorsArray[pcount+ndims] = 4;
			pcount++;
			evopcount++;
		}

		for(int p1 = 0; p1 < NProfileEvoPoly; p1++){
			for(int p2 = 0; p2 < numProfComponents; p2++){
				for(int p3 = 0; p3 < numEvoFitCoeff[p2]; p3++){

					printf("evo %i %i %i\n", p1, p2, p3);
					printf("evo Prior %i %i %.10g %.10g \n", pcount, profevodims, ProfileEvoPrior[0], ProfileEvoPrior[1]);

					PriorsArray[pcount] = ProfileEvoPrior[0];
					PriorsArray[pcount+ndims] = ProfileEvoPrior[1];

				        pcount++;
				}
			}

                }

		if(incProfileEnergyEvo == 2){
			printf("fit energy Width Prior %.10g %.10g \n", ProfileEnergyEvoPrior[0], ProfileEnergyEvoPrior[1]);

                        PriorsArray[pcount] = ProfileEnergyEvoPrior[0];
                        PriorsArray[pcount+ndims] = ProfileEnergyEvoPrior[1];
			pcount++;
		}


		printf("set priors \n");
		((MNStruct *)context)->PriorsArray = PriorsArray;
		((MNStruct *)context)->FixProfile = FixProfile;
		((MNStruct *)context)->MeanProfileShape = MeanProfileShape;	
		((MNStruct *)context)->MeanProfileEvo = MeanProfileEvo;	
		((MNStruct *)context)->MeanProfileBeta = MeanBeta;
		((MNStruct *)context)->MeanProfileStoc = MeanProfileStoc;
		//outputProfile(ndims);
		//return 0;

		assigncontext(context);
		printf("totalProfileFitCoeff %i \n", ((MNStruct *)context)->totalProfileFitCoeff);
		//return 0;
		if(interpolateProfile == 1){

			int Nbin = (int)ProfileInfo[0][1];
			long double timetointerpolate = pow(10.0, -9)*InterpolatedTime;
			int numtointerpolate = int(ceil(((MNStruct *)context)->ReferencePeriod/Nbin/timetointerpolate));
			long double finalInterpTime = ((MNStruct *)context)->ReferencePeriod/Nbin/numtointerpolate;
			((MNStruct *)context)->InterpolatedTime = finalInterpTime;

			int *numShapeToSave = new int[numProfComponents];
			int totalShapeToSave = 0;
			for(int i = 0; i < numProfComponents; i++){
				printf("num to save: %i %i %i %i \n", i, totalshapestoccoeff, numEvoCoeff[i],numProfileFitCoeff[i]); 
				numShapeToSave[i] = totalshapestoccoeff;
				if(numEvoCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numEvoCoeff[i];}
				if(numProfileFitCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numProfileFitCoeff[i];}
				totalShapeToSave += numShapeToSave[i];
			}

			printf("Final interp time %.10g %.10Lg %i\n", InterpolatedTime, finalInterpTime, numtointerpolate);

			//double ***StoredShapelets = new double**[numtointerpolate];
			//for(int i = 0; i < numtointerpolate; i++){
				//StoredShapelets[i] = new double*[Nbin];
				//for(int j = 0; j < Nbin; j++){
					// StoredShapelets[i][j] = new double[numShapeToSave];
					//for(int k = 0; k < numShapeToSave; k++){
						//StoredShapelets[i][j][k] = 0;
					//}//
				//}
			//}

			printf("allocate 1 %i %i %i \n", numtointerpolate, totalShapeToSave, Nbin);
			double **StoredShapeletsVec = new double*[numtointerpolate];
			double **StoredShapeletsJitterVec = new double*[numtointerpolate];
			for(int i = 0; i < numtointerpolate; i++){
				StoredShapeletsJitterVec[i] = new double[Nbin*totalShapeToSave];
				StoredShapeletsVec[i] = new double[Nbin*totalShapeToSave];
				for(int j = 0; j < Nbin*totalShapeToSave; j++){
					StoredShapeletsJitterVec[i][j] = 0;
					StoredShapeletsVec[i][j] = 0;
				}
			}

			printf("allocate 2 %i %i %i \n", numtointerpolate, totalShapeToSave, Nbin);
			double **InterpolatedMeanProfile = new double *[numtointerpolate];
			double **InterpolatedJitterProfile = new double *[numtointerpolate];
			double **InterpolatedWidthProfile = new double *[numtointerpolate];
			for(int i = 0; i < numtointerpolate; i++){
				InterpolatedMeanProfile[i] = new double[Nbin];
				InterpolatedJitterProfile[i] = new double[Nbin];
				InterpolatedWidthProfile[i] = new double[Nbin];
				for(int j = 0; j < Nbin; j++){
					InterpolatedMeanProfile[i][j] = 0;
					InterpolatedJitterProfile[i][j] = 0;
					InterpolatedWidthProfile[i][j] = 0;
				}	
			}

			double MaxShapeAmp = 0;
			PreComputeShapelets(StoredShapeletsVec, StoredShapeletsJitterVec, InterpolatedMeanProfile, InterpolatedJitterProfile, InterpolatedWidthProfile, finalInterpTime, numtointerpolate, MeanBeta, MaxShapeAmp);

			//for(int i = 0; i < numtointerpolate; i++){
			//	for(int j = 0; j < Nbin; j++){
				//	for(int k = 0; k < numShapeToSave; k++){	
					//	StoredShapeletsVec[i][j+k*Nbin] = StoredShapelets[i][j][k];	
					//}
				//}
			//}

			((MNStruct *)context)->InterpolatedShapeletsVec = StoredShapeletsVec;
			((MNStruct *)context)->InterpolatedJitterProfileVec = StoredShapeletsJitterVec;
			//((MNStruct *)context)->InterpolatedShapelets = StoredShapelets;
			((MNStruct *)context)->InterpolatedMeanProfile = InterpolatedMeanProfile;
			((MNStruct *)context)->InterpolatedJitterProfile = InterpolatedJitterProfile;
			((MNStruct *)context)->InterpolatedWidthProfile = InterpolatedWidthProfile;
			((MNStruct *)context)->NumToInterpolate = numtointerpolate;
			((MNStruct *)context)->MaxShapeAmp = MaxShapeAmp;
			//return 0;
		}

		if(FitTemplate > 0){

			 int Nbin = (int)ProfileInfo[0][1];
			long double timetointerpolate = pow(10.0, -9)*InterpolatedTime;
			int numtointerpolate = int(ceil(((MNStruct *)context)->ReferencePeriod/Nbin/timetointerpolate));
			long double finalInterpTime = ((MNStruct *)context)->ReferencePeriod/Nbin/numtointerpolate;
			((MNStruct *)context)->InterpolatedTime = finalInterpTime;
			((MNStruct *)context)->NumToInterpolate = numtointerpolate;


			int numTempFreqs = 0;
			getNumTempFreqs(numTempFreqs, context);

			double *TemplateFreqs = new double[numTempFreqs];
			double *TemplateChans = new double[numTempFreqs*Nbin];

			((MNStruct *)context)->TemplateFreqs = TemplateFreqs;
			((MNStruct *)context)->TemplateChans = TemplateChans;
			((MNStruct *)context)->numTempFreqs = numTempFreqs;

			
			getNumTempFreqs(numTempFreqs, context);
			Tscrunch(context);
			//return 0;
		}

		if(incProfileNoise == 3 || incProfileNoise == 4){

			std::ifstream ProfileNoiseFile;
			std::string ProfileNoiseFilename = "Powers.dat";
			ProfileNoiseFile.open(ProfileNoiseFilename.c_str());

			double **MLProfileNoise  = new double*[((MNStruct *)context)->numProfileEpochs];
			for(int i = 0; i < ((MNStruct *)context)->numProfileEpochs; i++){
				MLProfileNoise[i] = new double[ProfileNoiseCoeff+1];
			}

			for(int i = 0; i <  ((MNStruct *)context)->numProfileEpochs; i++){
				std::string line;
				getline(ProfileNoiseFile,line);
				std::istringstream myStream( line );
				std::istream_iterator< double > begin(myStream),eof;
				std::vector<double> ProfileNoiseVal(begin,eof);
				for(int j = 0; j < ProfileNoiseCoeff + 1; j++){
					MLProfileNoise[i][j] = ProfileNoiseVal[j];
					if(j==0)MLProfileNoise[i][j] = pow(10.0, MLProfileNoise[i][j]);	
					printf("Read NoiseFile: %i coeff %i = %g \n", i, j, MLProfileNoise[i][j]);
				}
			}
			((MNStruct *)context)->MLProfileNoise = MLProfileNoise;
		}
		else{
                        double **MLProfileNoise  = new double*[((MNStruct *)context)->numProfileEpochs];
                        for(int i = 0; i < ((MNStruct *)context)->numProfileEpochs; i++){
                                MLProfileNoise[i] = new double[1];
				MLProfileNoise[i][0] = 1;
                        }
			((MNStruct *)context)->MLProfileNoise = MLProfileNoise;

		}

		 ((MNStruct *)context)->Tspan = maxtspan;
		((MNStruct *)context)->numdims = ndims;

	
//		double *StoredDMVals = new double[100];
//		readDMVals(StoredDMVals);
//		((MNStruct *)context)->StoredDMVals = StoredDMVals;

		assigncontext(context);

		if(doMax == 1){
			NelderMeadOptimum(ndims);
			return 0;
		}

		if(sample==1 && sampler ==2){
			assignGHScontext(context);
			callGHS();
			return 0;
		}

		if(sample==1){
			if(FitTemplate > 0){
				if(sampler == 0){
					nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, Template2DProfLikeMNWrap, dumper, context);
					WriteTemplate2DProfLike(longname, ndims);
					return 0;

				}
				if(sampler == 1){
					chord::Sample(TemplateProfLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
				}
			}
			else{

				if(sampler == 0){
					//nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, SubIntStocProfLikeMNWrap, dumper, context);
					nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, ndims, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, ProfileDomainLikeMNWrap, dumper, context);
				}
				if(sampler == 1){
					chord::Sample(ProfileDomainLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);

					//chord::Sample(SubIntStocProfLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
				}
			}
		}

		if(FitTemplate > 0){
			WriteMaxTemplateProf(longname, ndims);
		}
		else{
			WriteProfileDomainLike(longname, ndims);
			//AllTOAWriteMaxLike(longname, ndims);
			//WriteSubIntStocProfLike(longname, ndims);
		}


#else
		printf("No installation of PSRCHIVE detected, cannot use GPTA likelihoods.  Please install PSRCHIVE, and set $PSRCHIVE to be the location of the installation directory\n");

#endif
//		return 0;
		

		//}

/*

	        printf("Not marginalising over anything %i \n", ndims);
	        double *PriorsArray=new double[2*ndims];
	        int NDerived = 0;
	        assigncontext(context);

	        int p = 0;
	        pcount = 0;
	        while(p < ndims){
	                printf("Prior %i %.10g %.10g %.10g \n", p,Dpriors[p][0], Dpriors[p][1], psr[0].rmsPre*pow(10.0,-6));
	                if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
	                        PriorsArray[pcount]=((MNStruct *)context)->Dpriors[p][0];
	                        PriorsArray[pcount+ndims]=((MNStruct *)context)->Dpriors[p][1];
	                        printf("param %i in, dim %i \n", p, pcount);
	                        pcount++;
	                }
	                 printf("param %i out\n", p);
	                p++;
	        }


		double *output = new double[5];

		/////////////////////////Sort Grade stuff////////////////////////////   


		int do_grades = 0;
		int *grades = new int[ndims];
		for(int i = 0; i < ndims; i++){
			grades[i] = 1;
		}

		int maxgrade = 1;
		int *grade_repeats = new int[maxgrade];
		for(int i = 0; i < maxgrade; i++){
			grade_repeats[i] = 1;
		}


		int *hypercube_indices = new int[ndims];
		int *physical_indices = new int[ndims];

		//Set indices, note this is fortran convention -> starts from 1/////

		for(int i = 0; i < ndims; i++){
			hypercube_indices[i] = i+1;
			physical_indices[i] = i+1;
		}


		if(incRED==0 && incDM == 0){


           	     chord::Sample(WhiteLogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);

		}
		else if(incRED==1 || incDM==1 ){

			printf("Time domain no longer supported, use incRed==3 for power law\n");

		}
		else if(incRED==2 || incDM==2 || incRED==3 || incDM==3 || incRED==4){
#ifdef HAVE_CULA
			if(useGPUS==1){

				chord::Sample(WhiteLogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
			}
			else{
				chord::Sample(WhiteLogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
			}
#else
		//	nested::run(IS,mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
#endif
		}
		else if(incRED==5 || incDM==5){

//				chord::Sample(LRedNumericalLogLike, ndims, NDerived, nlive, Nchords,  PriorsArray, root, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
		}


	*/
	
	}

	

	if(myrank == 0){
		readsummary(psr,longname, ndims,context,  Tempo2Fit,incRED, ndims, doTimeMargin, doJumpMargin,doLinearFit);



		time(&rawstoptime);
		rawstoptimeinfo = localtime(&rawstoptime);
		time_t stoptime = mktime(rawstoptimeinfo);
		double seconds = difftime(stoptime, starttime);
	        printf ( "The stop date/time was: %s", asctime(rawstoptimeinfo) );

		printf("Total Wall clock run time: %g \n", seconds);
	}


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
