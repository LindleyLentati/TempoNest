#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
 
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

#include <lbfgs.h>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "dpotrs.h"
#include "tempo2.h"
#include "TempoNest.h"
#include "dgesvd.h"
#include "qrdecomp.h"
#include "cholesky.h"
#include "T2toolkit.h"
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <cstring>
#include <limits>
#include <iostream>

#ifdef HAVE_MLAPACK
#include <mpack/mblas_qd.h>
#include <mpack/mlapack_qd.h>
#include <mpack/mblas_dd.h>
#include <mpack/mlapack_dd.h>
#endif

#include <guided_hmc.h>
# include <omp.h>
#include <fftw3.h>

FILE* test_gauss_outfile;
FILE* test_gauss_goutfile;

void write_gauss_ghs_extract(int* ndim,double* x,double* val,double* g);
void write_gauss_ghs_extract_with_logpostval(int* ndim,double* x,double* logpostval,double* g);
void write_gauss_ghs_extract_with_logpostval_and_grad(int* ndim,double* x,double* logpostval,double* g);
void nd_uncorr_gauss_neg_log_post(int* ndim,double* x,double* v,double* g);


void NewGHSProfileDomainLike(int* ndim,double* x,double* v,double* g);

void GetMaxAmps(int ndim, double *MaxAmps);
void GetMaxSteps(int ndim, double *MaxAmps, double *StepSize, double **SSM);
void GHSWriteProf(std::string ename);
//void WriteGHSProfileDomainLike(std::string longname, int ndim, double* Cube);

extern "C" void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);

void PosDefEVD(double *a, double *ev, int asize);
void rotate2Physical(double *xPrin, double *xPhy);
void rotate2Principal(double *gradPrin, double *gradPhy);

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    );
static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    );




using namespace std;

double GlobalMaxLike;
double *GlobalMaxLikeVec;
double **GlobalDMatrix;
double *GlobalDMatrixVec;
void *GHSglobalcontext;
double *GlobalStartPoint;
double **GlobalHessian;
double *GlobalParmSet;
double *GlobalStepSize;



void readextract(std::string ename);
void assignGHScontext(void *context){
        GHSglobalcontext=context;
}

void callFourierDomainGHS(int NBurn, int NSamp, int GHSresume){


	int id;



	double* StartPoint;
	double* StepSize;
	double **SSM;
	void (*nlp)(int*,double*,double*,double*);
	void (*wrt_ext)(int*,double*,double*,double*);
	double ScaleFactor;

	char* fl_pfx =    ((MNStruct *)GHSglobalcontext)->rootName; 
	int seed;
	int fb_int;
	int max_stp;
	int resume;
	char ext_file_name[128];
	char gext_file_name[128];
	char prof_file_name[128];
	int nburn=NBurn;
	int nsamp=NSamp;


	((MNStruct *)GHSglobalcontext)->diagonalGHS = 0; 
	((MNStruct *)GHSglobalcontext)->WriteNewML = 1;
	

	int readMLFile = 0;

	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = ((MNStruct *)GHSglobalcontext)->pulse->nobs;
	((MNStruct *)GHSglobalcontext)->TotalProfiles = TotalProfs;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////////////////////////////Get dimensionality/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int perProfDims = 2;
	int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*perProfDims;

	int perEpochDims = 0;
	int epochpriordims = 0;
	int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*perEpochDims;


	

	int globaldims = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + (1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly)*(((MNStruct *)GHSglobalcontext)->totshapecoeff-1)+((MNStruct *)GHSglobalcontext)->pulse->nSx;


	((MNStruct *)GHSglobalcontext)->GHSperProfDims = perProfDims;
	((MNStruct *)GHSglobalcontext)->GHSperEpochDims = perEpochDims;
	((MNStruct *)GHSglobalcontext)->GHSepochpriordims = epochpriordims;
	((MNStruct *)GHSglobalcontext)->GHSglobaldims = globaldims;


	int TimingDims = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;

	const int ndim = profdims+epochdims+epochpriordims+globaldims;

	double **PhysDMatrix = new double*[TotalProfs];
	for(int i = 0; i < TotalProfs; i++){ PhysDMatrix[i] = new double[TimingDims];}
	double *PhysDMatrixVec = new double[TotalProfs*TimingDims];

	int *TimingGradientSigns = new int[TimingDims]();

	getPhysDVector(GHSglobalcontext, PhysDMatrix, TotalProfs,TimingGradientSigns);
	GlobalDMatrix = PhysDMatrix;
	

	for(int i = 0; i < TotalProfs; i++){
		for(int j = 0; j < TimingDims; j++){
			PhysDMatrixVec[i+j*TotalProfs] = PhysDMatrix[i][j];
		}
	}
	GlobalDMatrixVec = PhysDMatrixVec;
	
	printf("dimensionality is: %i %i %i %i %i %i \n", profdims,epochdims,epochpriordims,globaldims, ndim, TotalProfs);	



	StartPoint = new double[ndim]();
	
	StepSize = new double[ndim]();
	GlobalMaxLikeVec = new double[ndim]();


	int HessianSize = TotalProfs;
	if(epochdims > 0)HessianSize+=((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	if(epochpriordims > 0)HessianSize+=1;
	if(globaldims > 0)HessianSize+=1;
	
	SSM = new double*[HessianSize];

	int HessCount = 0;
	for(int i = 0; i < TotalProfs; i++){ SSM[HessCount] = new double[perProfDims*perProfDims](); HessCount++;} 
	if(epochpriordims > 0){ SSM[HessCount] = new double[(epochpriordims+epochdims+globaldims)*(epochpriordims+epochdims+globaldims)](); HessCount++;}
	if(epochdims > 0){ for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfileEpochs; i++){ SSM[HessCount] = new double[perEpochDims*perEpochDims](); HessCount++;} }
	

	if(globaldims > 0){ SSM[HessCount] = new double[globaldims*globaldims](); HessCount++; }

	printf("Getting Max and StepSize %i \n", ndim);
	GetMaxAmps(ndim, StartPoint);
	GetMaxSteps(ndim, StartPoint, StepSize, SSM);

	GlobalMaxLike = -1*pow(10, 100);
	
	GlobalStartPoint = StartPoint;
	GlobalHessian = SSM;

	for(int i = 0; i < ndim; i++){
		GlobalMaxLikeVec[i] = StartPoint[i];
	}


	double *PrincipalStartPoint = new double[ndim]();

	rotate2Principal(PrincipalStartPoint, StartPoint);


	for(int i = 0; i < ndim; i++){
		StepSize[i] = 1.0/sqrt(StepSize[i]);
	}

	for(int i = 0; i < ndim; i++){
		if(i>=0){printf("Final Step Size %i %g %g %g \n", i, StartPoint[i], StepSize[i], PrincipalStartPoint[i]);}
	}
	GlobalStepSize = new double[ndim];
	GlobalStepSize = StepSize;

	
	
	nlp=&NewGHSProfileDomainLike;
	wrt_ext=&write_gauss_ghs_extract_with_logpostval_and_grad; 
	

/*  dimensionality scaling factor (a value between 0 and 1) ideal choise is the one which gives you ~68% acceptance rate */
	ScaleFactor = 0.4;


/*  feed back to console interval */

	fb_int=100;

/*  maximum number of steps in the leapforg (10 is fine) */

	max_stp=10;

/* resume from previous run? 0(no) or 1(yes) */

	resume=GHSresume;
	printf("GHS Resume is: %i \n", resume);
/* random number generator seed */

	seed=1235;
	
	strcpy(ext_file_name,fl_pfx);
	strcat(ext_file_name,".extract.dat");
	
	strcpy(gext_file_name,fl_pfx);
	strcat(gext_file_name,".gextract.dat");
	

	strcpy(prof_file_name,fl_pfx);


	if(resume==1)
	{
		test_gauss_outfile=fopen(ext_file_name,"a");
		test_gauss_goutfile=fopen(gext_file_name,"a");
	}
	else
	{
		test_gauss_outfile=fopen(ext_file_name,"w");
		test_gauss_goutfile=fopen(gext_file_name,"w");
	}
	run_guided_hmc(ndim, PrincipalStartPoint, ScaleFactor, max_stp, StepSize, fl_pfx, seed, resume, fb_int, nlp, wrt_ext, nburn, nsamp);
		
	if(test_gauss_outfile) fclose(test_gauss_outfile);


	for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
		if(((MNStruct *)GHSglobalcontext)->pulse->fitJump[k] == 0){

			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k] = ((MNStruct *)GHSglobalcontext)->PreJumpVals[k];
		}
	}

	

	readextract(ext_file_name);
	GHSWriteProf(ext_file_name);
	OutputMLFiles(ndim, GlobalMaxLikeVec, GlobalMaxLike, 2*TotalProfs+epochdims+epochpriordims);
	//WriteGHSProfileDomainLike(prof_file_name, ndim, GlobalMaxLikeVec);

	delete[] StartPoint;
	delete[] StepSize;

}

void write_gauss_ghs_extract(int* ndim,double* x,double* val,double* g)
{
	int i;
	if(test_gauss_outfile != NULL)
	{
		for(i=0;i<*ndim-1;++i)
		{
			fprintf(test_gauss_outfile,"%e ",x[i]);
		}
		fprintf(test_gauss_outfile,"%e\n",x[*ndim-1]);
	}
	else
	{
		printf("\nERROR IN WRITING GHS EXTRACT!\n");
	}
}

void write_newMLLike(int ndim,double* x, double logpostval)
{
	FILE* NewMLoutfile;
	char* ml_pfx="NewML.dat";
	NewMLoutfile=fopen(ml_pfx,"w");

	for(int i=0;i<ndim;++i){
		fprintf(NewMLoutfile,"%.10g \n",x[i]);
		GlobalMaxLikeVec[i] = x[i];
	}
	fprintf(NewMLoutfile,"%.10g\n",logpostval); 



	fclose(NewMLoutfile);
}


void write_gauss_ghs_extract_with_logpostval(int* ndim,double* x,double* logpostval,double* g)
{

	double *xPhys = new double[*ndim];

	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Physical(x, xPhys);
	}
	else{
		for(int i = 0; i < *ndim; i++){
			xPhys[i] = x[i];
		}
	}


	int i;
	if(test_gauss_outfile != NULL)
	{
		for(i=0;i<*ndim;++i)
		{
			fprintf(test_gauss_outfile,"%.10e ",xPhys[i]);
		}
		fprintf(test_gauss_outfile,"%.10e\n",*logpostval); 
	}
	else
	{
		printf("\nERROR IN WRITING GHS EXTRACT!\n");
	}

	delete[] xPhys;
}

void write_gauss_ghs_extract_with_logpostval_and_grad(int* ndim,double* x,double* logpostval,double* g)
{

	double *xPhys = new double[*ndim];
	double *gradPhys = new double[*ndim];

	//int start = 0;
	int start = 2*((MNStruct *)GHSglobalcontext)->TotalProfiles;
	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Physical(x, xPhys);
		rotate2Physical(g, gradPhys);
	}
	else{
		for(int i = 0; i < *ndim; i++){
			xPhys[i] = x[i];
			gradPhys[i] = g[i];
		}
	}


	int i;
	if(test_gauss_outfile != NULL && test_gauss_goutfile != NULL)
	{

		//printf("Writing Sample: %i %i \n", start, *ndim);
		for(i=start;i<*ndim;++i)
		{
			//printf("Writing Sample: %i %g \n", i,xPhys[i] );	
			fprintf(test_gauss_outfile,"%.12g ",xPhys[i]);
			fprintf(test_gauss_goutfile,"%.12g ",gradPhys[i]);
		}
		fprintf(test_gauss_outfile,"%.12g \n",*logpostval);
		fprintf(test_gauss_goutfile,"%.12g \n",*logpostval);
	}
	else
	{
		printf("\nERROR IN WRITING GHS EXTRACT!\n");
	}

	delete[] xPhys;
	delete[] gradPhys;
}




void GetMaxAmps(int ndim, double *MaxAmps){





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////////////////////////////Get dimensionality/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        int perProfDims = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int perEpochDims = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;

	int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*perProfDims;
	int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*perEpochDims;




	int dotime = 0;
	int debug = ((MNStruct *)GHSglobalcontext)->debug; 

	for(int i = 0; i < ndim; i++){
		MaxAmps[i] = 0;
	}


	
	if(debug == 1)printf("In Max Amps \n");



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the BATS///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = ((MNStruct *)GHSglobalcontext)->pulse->nobs;


	 long double *ProfileBats=new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

		ProfileBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
	      
		ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;

	 }


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	long double phase = ((MNStruct *)GHSglobalcontext)->LDpriors[0][0];	
	phase = phase*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;

	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

     		 ModelBats[i] -=  phase;
	}	


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Profile Parameters//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int NFBasis = ((MNStruct *)GHSglobalcontext)->NFBasis;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
	}


	double **ProfCoeffs=new double*[1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){ProfCoeffs[i] = new double[totshapecoeff]();}

	for(int i = 0; i < totshapecoeff; i++){
		ProfCoeffs[0][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[i];
	}
	for(int p = 1; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){	
		for(int i = 1; i < totshapecoeff; i++){
			ProfCoeffs[p][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i];
		}
	}

	double *ProfileScatter;
	if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){
		ProfileScatter =  new double[((MNStruct *)GHSglobalcontext)->pulse->nSx];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nSx; i++){
			printf("Scattering Params: %i %g \n", i, (double)((MNStruct *)GHSglobalcontext)->pulse->param[param_sx].val[i]);
			ProfileScatter[i] =  pow(10.0, ((MNStruct *)GHSglobalcontext)->pulse->param[param_sx].val[i]); 
		}
	}




/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;


	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;

	
	for(int ep = minep; ep < maxep; ep++){

		double *TotalProfCoeffs = new double[totshapecoeff];


		int t = 0;

		for(int sep = 0; sep < ep; sep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[sep];
		}
		

		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];



		for(int ch = 0; ch < NChanInEpoch; ch++){

			int nTOA = t;



			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];


			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;


			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;

			double *RolledData  = new double[2*NFBasis]();
			double *ProfileVec  = new double[2*NFBasis]();
			double *ProfGradVec = new double[2*NFBasis]();
			double *RolledData2  = new double[2*NFBasis]();



			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get and modify binpos///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////
		    
			long double binpos = ModelBats[nTOA]; 


			long double binposSec = (ProfileBats[t] - ModelBats[nTOA])*SECDAY;
			double DFP = double(FoldingPeriod);
			long double WrapBinPos = -FoldingPeriod/2 + fmod(FoldingPeriod + fmod(binposSec+FoldingPeriod/2, FoldingPeriod), FoldingPeriod);

			long double maxbinvalue = (FoldingPeriod/Nbins);
			double DNewInterpBin = double(fmod(maxbinvalue + fmod(-WrapBinPos, maxbinvalue), maxbinvalue)/((MNStruct *)GHSglobalcontext)->InterpolatedTime);
			int InterpBin = floor(DNewInterpBin+0.5);
			InterpBin = InterpBin%((MNStruct *)GHSglobalcontext)->NumToInterpolate;

			long double NewWholeBinTime = -WrapBinPos-InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime;
			int RollBin = int(round(NewWholeBinTime/maxbinvalue));
			int WrapRollBin = ((Nbins) + -RollBin % (Nbins)) % (Nbins);


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freq - reffreq)/1000.0;


			for(int i =0; i < totshapecoeff; i++){
				TotalProfCoeffs[i]=0;	
			}				

			int cpos = 0;
	
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				for(int p = 0; p < 1 + ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i = 0; i < numcoeff[c]; i++){
						TotalProfCoeffs[i+cpos] += ProfCoeffs[p][i+cpos]*pow(freqdiff, p);						
					}
				}

				cpos += numcoeff[c];
			}




			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[InterpBin], TotalProfCoeffs, ProfileVec,2*NFBasis,totshapecoeff,'N');
			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFJitterBasis[InterpBin], TotalProfCoeffs, ProfGradVec,2*NFBasis,totshapecoeff,'N');         



		///////////////////////////////////////////Marginalise over arbitrary offset and absolute amplitude////////////////////////////////////////////////////////////


			double MNM = 0;
			double dM = 0;


			for(int j = 0; j < NFBasis; j++){

				double freq = (j+1.0)/ProfNbins;
				double theta = 2*M_PI*WrapRollBin*freq;
				double RealShift = cos(theta);
				double ImagShift = sin(-theta);

				double RData = ((MNStruct *)GHSglobalcontext)->FourierProfileData[t][1+j];
				double IData = ((MNStruct *)GHSglobalcontext)->FourierProfileData[t][1+j+ProfNbins/2+1];

				double RealRolledData = RData*RealShift - IData*ImagShift;
				double ImagRolledData = RData*ImagShift + IData*RealShift;

				RolledData[j] = RealRolledData;
				RolledData[j+NFBasis] = ImagRolledData;


			}

			if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){

				int Sindex = ((MNStruct *)GHSglobalcontext)->ScatterIndex[nTOA];

				double ISS = 1.0/(pow( ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB, 4)/pow(10.0, 9.0*4.0));
				double tau = ProfileScatter[Sindex];
				double STime = tau*ISS;

				
				for(int j = 0; j < NFBasis; j++){


					double f = (j+1.0)/((MNStruct *)GHSglobalcontext)->ReferencePeriod;
					double w = 2.0*M_PI*f;
					double RConv = 1.0/(w*w*STime*STime+1); 
					double IConv = -w*STime/(w*w*STime*STime+1);

					double RProf = ProfileVec[j];
					double IProf = ProfileVec[j+NFBasis];

					double ScatterRealProf =  RProf*RConv - IProf*IConv;
					double ScatterImagProf =  RProf*IConv + IProf*RConv;

					 ProfileVec[j] = ScatterRealProf;
					 ProfileVec[j+NFBasis] = ScatterImagProf;

				}
			}


			for(int j = 0; j < NFBasis; j++){


				MNM +=  ProfileVec[j]*ProfileVec[j]  + ProfileVec[j+NFBasis]*ProfileVec[j+NFBasis];
				dM += ProfileVec[j]*RolledData[j] + ProfileVec[j+NFBasis]*RolledData[j+NFBasis];


			}	



			double ProfileAmp = dM/MNM;
			double ProfileSigma = 0;
			double Chisq=0;
			
			for(int j = 0; j < NFBasis; j++){

				double RealRes = ProfileAmp*ProfileVec[j] -  RolledData[j];
				double ImRes = ProfileAmp*ProfileVec[j+NFBasis] - RolledData[j+NFBasis];

				Chisq += RealRes*RealRes + ImRes*ImRes;

			}
			ProfileSigma = sqrt(Chisq/(2*NFBasis));
//			printf("Max: %i %g %g  \n", t, ProfileSigma, ProfileAmp);

			MaxAmps[TotalProfs + t] = ProfileAmp;
			MaxAmps[t] = ProfileSigma;
			

			delete[] ProfileVec;
			delete[] ProfGradVec;
			delete[] RolledData;

			t++;
		}



///////////////////////////////////////////


		delete[] TotalProfCoeffs;

////////////////////////////////////////////
	
	}


	int pcount = 2*TotalProfs;

	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + (1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly)*(((MNStruct *)GHSglobalcontext)->totshapecoeff-1); i++){
		MaxAmps[pcount] = 0;
		pcount++;
	}
	
	if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nSx; i++){
			MaxAmps[pcount] = ((MNStruct *)GHSglobalcontext)->pulse->param[param_sx].val[i];
			pcount++; 
		}
	}



	delete[] ProfileBats;
	delete[] ModelBats;


	for(int j = 0; j < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] ProfCoeffs[j];
	}
	delete[] ProfCoeffs;



}

void GetMaxSteps(int ndim, double *MaxAmps, double *StepSize, double **SSM){


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////////////////////////////Get dimensionality/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        int perProfDims = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int perEpochDims = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;

	int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*perProfDims;
	int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*perEpochDims;

	int LinearParams = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + (1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly)*(((MNStruct *)GHSglobalcontext)->totshapecoeff-1);
	int HessSize = epochpriordims+epochdims+globaldims;



	int dotime = 0;
	int debug = ((MNStruct *)GHSglobalcontext)->debug; 

	for(int i = 0; i < ndim; i++){
		StepSize[i] = 0;
	}


	
	if(debug == 1)printf("In Max Steps \n");



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the BATS///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = ((MNStruct *)GHSglobalcontext)->pulse->nobs;


	 long double *ProfileBats=new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

		ProfileBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
	      
		ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;

	 }


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	long double phase = ((MNStruct *)GHSglobalcontext)->LDpriors[0][0];	
	phase = phase*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;

	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

     		 ModelBats[i] -=  phase;
	}	


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Profile Parameters//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int NFBasis = ((MNStruct *)GHSglobalcontext)->NFBasis;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
	}


	double **ProfCoeffs=new double*[1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){ProfCoeffs[i] = new double[totshapecoeff]();}

	for(int i = 0; i < totshapecoeff; i++){
		ProfCoeffs[0][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[i];
	}
	for(int p = 1; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){	
		for(int i =1; i < totshapecoeff; i++){
			ProfCoeffs[p][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i];
		}
	}


	double *ProfileScatter;
	if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){
		ProfileScatter =  new double[((MNStruct *)GHSglobalcontext)->pulse->nSx];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nSx; i++){
			printf("Scattering Params: %i %g \n", i, (double)((MNStruct *)GHSglobalcontext)->pulse->param[param_sx].val[i]);
			ProfileScatter[i] =  pow(10.0, ((MNStruct *)GHSglobalcontext)->pulse->param[param_sx].val[i]); 
		}
	}



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;


	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;

	
	for(int ep = minep; ep < maxep; ep++){

		double *TotalProfCoeffs = new double[totshapecoeff];


		int t = 0;

		for(int sep = 0; sep < ep; sep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[sep];
		}
		

		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];



		for(int ch = 0; ch < NChanInEpoch; ch++){

			int nTOA = t;

			double ProfileAmp = MaxAmps[TotalProfs + t];
			double ProfileSigma = MaxAmps[t];

			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];


			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;


			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;

			double *RolledData  = new double[2*NFBasis]();
			double *ProfileVec  = new double[2*NFBasis]();
			double *ProfGradVec = new double[2*NFBasis]();

			double *LinearMatrix = new double[LinearParams*2*NFBasis]();
			double *LinearHess = new double[LinearParams*LinearParams]();



			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get and modify binpos///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////
		    
			long double binpos = ModelBats[nTOA]; 


			long double binposSec = (ProfileBats[t] - ModelBats[nTOA])*SECDAY;
			double DFP = double(FoldingPeriod);
			long double WrapBinPos = -FoldingPeriod/2 + fmod(FoldingPeriod + fmod(binposSec+FoldingPeriod/2, FoldingPeriod), FoldingPeriod);

			long double maxbinvalue = (FoldingPeriod/Nbins);
			double DNewInterpBin = double(fmod(maxbinvalue + fmod(-WrapBinPos, maxbinvalue), maxbinvalue)/((MNStruct *)GHSglobalcontext)->InterpolatedTime);
			int InterpBin = floor(DNewInterpBin+0.5);
			InterpBin = InterpBin%((MNStruct *)GHSglobalcontext)->NumToInterpolate;

			long double NewWholeBinTime = -WrapBinPos-InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime;
			int RollBin = int(round(NewWholeBinTime/maxbinvalue));
			int WrapRollBin = ((Nbins) + -RollBin % (Nbins)) % (Nbins);


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freq - reffreq)/1000.0;


			for(int i =0; i < totshapecoeff; i++){
				TotalProfCoeffs[i]=0;	
			}				

			int cpos = 0;
	
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				for(int p = 0; p < 1 + ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i = 0; i < numcoeff[c]; i++){
						TotalProfCoeffs[i+cpos] += ProfCoeffs[p][i+cpos]*pow(freqdiff, p);						
					}
				}

				cpos += numcoeff[c];
			}



			double *OneInterpMatrix = new double[totshapecoeff*2*NFBasis];
			for(int k = 0; k < totshapecoeff; k++){
				for(int j = 0; j < 2*NFBasis; j++){
					OneInterpMatrix[j+k*2*NFBasis] = ((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[InterpBin][j+k*2*NFBasis];
				}
			}


			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[InterpBin], TotalProfCoeffs, ProfileVec,2*NFBasis,totshapecoeff,'N');
			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFJitterBasis[InterpBin], TotalProfCoeffs, ProfGradVec,2*NFBasis,totshapecoeff,'N');         



		///////////////////////////////////////////Marginalise over arbitrary offset and absolute amplitude////////////////////////////////////////////////////////////


			double AmpStep = 0;
			double SigmaStep = 0;

			for(int j = 0; j < NFBasis; j++){

				double freq = (j+1.0)/ProfNbins;
				double theta = 2*M_PI*WrapRollBin*freq;
				double RealShift = cos(theta);
				double ImagShift = sin(-theta);

				double RData = ((MNStruct *)GHSglobalcontext)->FourierProfileData[t][1+j];
				double IData = ((MNStruct *)GHSglobalcontext)->FourierProfileData[t][1+j+ProfNbins/2+1];

				double RealRolledData = RData*RealShift - IData*ImagShift;
				double ImagRolledData = RData*ImagShift + IData*RealShift;

				RolledData[j] = RealRolledData;
				RolledData[j+NFBasis] = ImagRolledData;


			}
			double *OrigProf;
			if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){

				OrigProf = new double[2*NFBasis];
				int Sindex = ((MNStruct *)GHSglobalcontext)->ScatterIndex[nTOA];

				double ISS = 1.0/(pow( ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB, 4)/pow(10.0, 9.0*4.0));
				double tau = ProfileScatter[Sindex];
				double STime = tau*ISS;

				
				for(int j = 0; j < NFBasis; j++){


					OrigProf[j] = ProfileVec[j];
					OrigProf[j+NFBasis] = ProfileVec[j+NFBasis];
					double f = (j+1.0)/((MNStruct *)GHSglobalcontext)->ReferencePeriod;
					double w = 2.0*M_PI*f;
					double RConv = 1.0/(w*w*STime*STime+1); 
					double IConv = -w*STime/(w*w*STime*STime+1);

					double RProf = ProfileVec[j];
					double IProf = ProfileVec[j+NFBasis];

					double ScatterRealProf =  RProf*RConv - IProf*IConv;
					double ScatterImagProf =  RProf*IConv + IProf*RConv;

					ProfileVec[j] = ScatterRealProf;
					ProfileVec[j+NFBasis] = ScatterImagProf;


					RProf = ProfGradVec[j];
					IProf = ProfGradVec[j+NFBasis];

					ScatterRealProf =  RProf*RConv - IProf*IConv;
					ScatterImagProf =  RProf*IConv + IProf*RConv;

					ProfGradVec[j] = ScatterRealProf;
					ProfGradVec[j+NFBasis] = ScatterImagProf;

				}

				for(int k = 0; k < totshapecoeff; k++){
					for(int j = 0; j < NFBasis; j++){


						double f = (j+1.0)/((MNStruct *)GHSglobalcontext)->ReferencePeriod;
						double w = 2.0*M_PI*f;
						double RConv = 1.0/(w*w*STime*STime+1); 
						double IConv = -w*STime/(w*w*STime*STime+1);

						double RProf = OneInterpMatrix[j+k*2*NFBasis];
						double IProf = OneInterpMatrix[(j+NFBasis) + k*2*NFBasis];

						double ScatterRealProf =  RProf*RConv - IProf*IConv;
						double ScatterImagProf =  RProf*IConv + IProf*RConv;

						 OneInterpMatrix[j+k*2*NFBasis] = ScatterRealProf;
						 OneInterpMatrix[(j+NFBasis) + k*2*NFBasis] = ScatterImagProf;

					}
				}
			}



			for(int j = 0; j < NFBasis; j++){

				double RealRes = ProfileAmp*ProfileVec[j] -  RolledData[j];
				double ImRes = ProfileAmp*ProfileVec[j+NFBasis] - RolledData[j+NFBasis];


				AmpStep +=  (ProfileVec[j]*ProfileVec[j]  + ProfileVec[j+NFBasis]*ProfileVec[j+NFBasis])/(ProfileSigma*ProfileSigma);
				SigmaStep += 3*(RealRes*RealRes + ImRes*ImRes)/(ProfileSigma*ProfileSigma*ProfileSigma*ProfileSigma) - 1.0/(ProfileSigma*ProfileSigma);

			}


			StepSize[TotalProfs + nTOA] = AmpStep;
			StepSize[nTOA] = SigmaStep;
			


			printf("Step Size: %i %g %g %g %g\n", t, AmpStep, SigmaStep,  1.0/sqrt(AmpStep), 1.0/sqrt(SigmaStep));


			double Scale = ProfileAmp/ProfileSigma;

			int Lincount = 0;
			for(int j = 0; j < 2*NFBasis; j++){
				LinearMatrix[j + Lincount*2*NFBasis] = Scale*ProfGradVec[j]*ReferencePeriod;
			}
			Lincount += 1;

			//Hessian for Timing Model

			for(int c = 1; c < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; c++){
				double gradsquare = 0;
				for(int j = 0; j < 2*NFBasis; j++){
					//if(t==1){printf("Timing: %i %i %g %g %g \n", InterpBin, j,  ProfGradVec[j],Scale,GlobalDMatrix[t][c]);}
					gradsquare += ProfGradVec[j]*ProfGradVec[j];
					LinearMatrix[j + Lincount*2*NFBasis] = ProfGradVec[j]*Scale*GlobalDMatrix[t][c];
				}
				Lincount += 1;
				//printf("grad square: %i %g %g \n", t, gradsquare, gradsquare*Scale*Scale);
			}

			//Hessian for Shapelet parameters



			for(int p = 0; p < 1 + ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
				double fval = pow(freqdiff, p);	
				for(int k = 1; k < totshapecoeff; k++){
					for(int j = 0; j < 2*NFBasis; j++){
						LinearMatrix[j + Lincount*2*NFBasis] = fval*OneInterpMatrix[j+k*2*NFBasis]*Scale;
					}
					Lincount += 1;				
				}
			}
			delete[] OneInterpMatrix;


			vector_dgemm(LinearMatrix, LinearMatrix , LinearHess, 2*NFBasis, LinearParams, 2*NFBasis, LinearParams, 'T', 'N');

			int HessOffset = 0;

			for(int j = 0; j < LinearParams; j++){
				for(int k = 0; k < LinearParams; k++){
					SSM[TotalProfs][(j+HessOffset) + (k+HessOffset)*HessSize] += LinearHess[j + k*LinearParams];
					//if(t<2){printf("Linear Hess: %i %i %i %g \n", t, j, k, LinearHess[j + k*LinearParams]);}
				}
			}


			int pcount = LinearParams+HessOffset;
			
			if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){

				int Sindex = ((MNStruct *)GHSglobalcontext)->ScatterIndex[nTOA];

				double ISS = 1.0/(pow( ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB, 4)/pow(10.0, 9.0*4.0));
				double tau = ProfileScatter[Sindex];
				double STime = tau*ISS;
				double pnoise = ProfileSigma;


				double FullRealHess = 0;
				double FullImagHess = 0;

				double *ScatterGrad = new double[2*NFBasis]();
				double *LinearScatterCross = new double[LinearParams]();
				
				for(int j = 0; j < NFBasis; j++){


					double f = (j+1.0)/((MNStruct *)GHSglobalcontext)->ReferencePeriod;
					double w = 2.0*M_PI*f;
					double RConv = 1.0/(w*w*STime*STime+1); 
					double IConv = -w*STime/(w*w*STime*STime+1);

					double RProf = ProfileAmp*OrigProf[j];
					double IProf = ProfileAmp*OrigProf[j+NFBasis];

					double HessDenom = 1.0/pow(1.0 + tau*tau*w*w*ISS*ISS, 3);
					double GradDenom = 1.0/pow(1.0 + tau*tau*w*w*ISS*ISS, 2);

					double RealFunc =  RolledData[j] - RProf*RConv + IProf*IConv;

					

					double RealGrad = 2*tau*tau*ISS*ISS*w*w*log(10.0)*GradDenom*RProf + tau*ISS*w*(tau*tau*ISS*ISS*w*w - 1)*log(10.0)*GradDenom*IProf;
					double RealHess = -(4*tau*tau*ISS*ISS*w*w*(tau*tau*ISS*ISS*w*w - 1)*log(10.0)*log(10.0))*HessDenom*RProf - tau*ISS*w*(1+tau*tau*ISS*ISS*w*w*(tau*tau*ISS*ISS*w*w - 6))*log(10.0)*log(10.0)*HessDenom*IProf;

					FullRealHess += 1*(RealHess*RealFunc + RealGrad*RealGrad)*(1.0/pnoise/pnoise);

					//printf("RHEss: %i %g %g %g %g \n", j, RealFunc, RealGrad, RealHess, FullRealHess);
					

					double ImagFunc = RolledData[j+NFBasis] - RProf*IConv - IProf*RConv;
					double ImagGrad = 2*tau*tau*ISS*ISS*w*w*log(10.0)*GradDenom*IProf - tau*ISS*w*(tau*tau*ISS*ISS*w*w - 1)*log(10.0)*GradDenom*RProf;
					double ImagHess = -(4*tau*tau*ISS*ISS*w*w*(tau*tau*ISS*ISS*w*w - 1)*log(10.0)*log(10.0))*HessDenom*IProf + tau*ISS*w*(1+tau*tau*ISS*ISS*w*w*(tau*tau*ISS*ISS*w*w - 6))*log(10.0)*log(10.0)*HessDenom*RProf;

					//printf("IFunc: %i %g %g \n", j, RolledData[j+NFBasis] , RProf*IConv + IProf*RConv);

					FullImagHess += 1*(ImagHess*ImagFunc + ImagGrad*ImagGrad)*(1.0/pnoise/pnoise);
					//printf("RHEss: %i %g %g %g %g \n", j, ImagFunc, ImagGrad, ImagHess, FullImagHess);

					ScatterGrad[j] = RealGrad/pnoise;
					ScatterGrad[j+NFBasis] = ImagGrad/pnoise;
				}	

		
				SSM[TotalProfs][pcount+Sindex + (pcount+Sindex)*HessSize] += FullImagHess + FullRealHess;
	

				vector_dgemv(LinearMatrix, ScatterGrad, LinearScatterCross,2*NFBasis,LinearParams,'T');
				
				for(int j = 0; j < LinearParams; j++){
					//printf("Cross: %i %g \n", j, LinearScatterCross[j]);

					SSM[TotalProfs][(HessOffset+j) + (pcount+Sindex)*HessSize] += -LinearScatterCross[j];
					SSM[TotalProfs][(pcount+Sindex) + (HessOffset+j)*HessSize] += -LinearScatterCross[j];
				}

				//Hessian[pcount+c,pcount+c] += np.sum(profhess)
				//Hessian[:LinearSize, pcount+c] += LinearScatterCross
				//Hessian[pcount+c, :LinearSize] += LinearScatterCross

				delete[] OrigProf;
				delete[] ScatterGrad;
				delete[] LinearScatterCross;
			}
			pcount++;

			delete[] ProfileVec;
			delete[] ProfGradVec;
			delete[] RolledData;
			delete[] LinearMatrix;
			delete[] LinearHess;

			t++;
		}



///////////////////////////////////////////


		delete[] TotalProfCoeffs;

////////////////////////////////////////////
	
	}


	int pcount = 2*TotalProfs;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Get Epoch Prior EVDs///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	int HessCount = TotalProfs;
	int dimCount = TotalProfs*2;


	int PhaseDim = epochpriordims+epochdims;


	double expansionfactor = ((MNStruct *)GHSglobalcontext)->phasePriorExpansion;
	((MNStruct *)GHSglobalcontext)->PhasePrior = (1.0/sqrt(SSM[HessCount][PhaseDim + PhaseDim*(epochpriordims+epochdims+globaldims)]))*expansionfactor;
	printf("using prior on phase: %i %.15g %.15g \n", PhaseDim + PhaseDim*(epochpriordims+epochdims+globaldims), ((MNStruct *)GHSglobalcontext)->PhasePrior, ((MNStruct *)GHSglobalcontext)->phasePriorExpansion);
	SSM[HessCount][PhaseDim + PhaseDim*(epochpriordims+epochdims+globaldims)] += 1.0/(((MNStruct *)GHSglobalcontext)->PhasePrior*((MNStruct *)GHSglobalcontext)->PhasePrior);


	double *PEVals = new double[epochpriordims+epochdims+globaldims]();
	
	//for(int c = 0; c < epochpriordims+epochdims+globaldims; c++){
	//	for(int d = 0; d < epochpriordims+epochdims+globaldims; d++){
	//		printf("A[%i][%i] =  %.15g \n", c, d, SSM[HessCount][c + d*(epochpriordims+epochdims+globaldims)]);
	//	}
	//}

	PosDefEVD(SSM[HessCount], PEVals, epochpriordims+epochdims+globaldims);

//	for(int c = 0; c < epochpriordims+epochdims+globaldims; c++){
//		for(int d = 0; d < epochpriordims+epochdims+globaldims; d++){
//			printf("after: %i %i %g \n", c, d, SSM[HessCount][c + d*(epochpriordims+epochdims+globaldims)]);
//		}
//	}


	for(int j = 0; j < epochpriordims+epochdims+globaldims; j++){
		StepSize[dimCount + j] = PEVals[j];
		printf("Step: %i %g %g \n", j,StepSize[dimCount + j], 1.0/sqrt(StepSize[dimCount + j]) );

	}
	delete[] PEVals;
	


	HessCount += 1+((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	dimCount += epochpriordims+epochdims+globaldims;



	



	delete[] ProfileBats;
	delete[] ModelBats;


	for(int j = 0; j < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] ProfCoeffs[j];
	}
	delete[] ProfCoeffs;
}


void PosDefEVD(double *a, double *ev, int asize){

	//printf("in EVD with %i \n", asize);
        int iwkopt;
        double wkopt;

	double VL = 0.0;
	double VU = 0.0;
	int IL = 0;
	int IU = 0;
	int m = 0;

	double abstol = -1.0;

	double *Z = new double[asize*asize]();
	int *ISUPPZ = new int[2*asize]();
	int info = 0;

	

        /* Query and allocate the optimal workspace */
        int LWORK = -1;
        int LIWORK = -1;

	//printf("About to  EVD with %i \n", asize);
	dsyevr_("V", "A", "L", &asize, a, &asize, &VL, &VU, &IL, &IU, &abstol,  &m, ev, Z, &asize, ISUPPZ, &wkopt, &LWORK, &iwkopt, &LIWORK, &info);
	//printf("Just did EVD with %i \n", asize);
	//printf("got work size %i %i %i \n", (int)wkopt, (int)iwkopt, asize);

        LWORK = (int)wkopt;
        LIWORK = iwkopt;


        int* iwork = new int[LIWORK]();
        double* work = new double[LWORK]();

        /* Solve eigenproblem */
	dsyevr_("V", "A", "L", &asize, a, &asize, &VL, &VU, &IL, &IU, &abstol,  &m, ev, Z, &asize, ISUPPZ, work, &LWORK, iwork, &LIWORK, &info);

        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }

	int negeigen = 0;
	for(int i = 0; i < asize; i++){
		if( ev[i] < 0 || isnan(ev[i])){negeigen = 1; printf("this eigen is negative %i %g \n", i, ev[i]);}
	}

	if(negeigen == 0){
		for(int i = 0; i < asize; i++){
			for(int j = 0; j < asize; j++){
				a[i + j*asize] = Z[i+j*asize];
				//printf("EVev: %i %i %g \n", i,j,a[i + j*asize]);
			}
		}
	}
	else{
		for(int i = 0; i < asize; i++){
			for(int j = 0; j < asize; j++){
				
				if(i==j){
					
					if( ev[i] < 0 || isnan(ev[i])){
						printf("Bad Eigen: %i %i %g \n", i,j, ev[i]);
						ev[i] = 1.0/(0.3*0.3);
					}
					else{
						ev[i] = a[i + j*asize];
						if(ev[i] > 1)ev[i] = 1;
					}

					a[i + j*asize] = 1; 
				}
				else{
					a[i + j*asize] = 0;
				}
				
			}
		}
	}

	for(int i = 0; i < asize; i++){
//		printf("EVal: %i %g \n", i, ev[i]);
	}


	delete[] iwork;
	delete[] work;
	delete[] ISUPPZ;
	delete[] Z;

	//sleep(5);

}


void rotate2Principal(double *gradPrin, double *gradPhy){

        int DimsPerProf = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int DimsPerEpoch = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;

        int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*DimsPerProf;
        int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*DimsPerEpoch;






	int DimCount = 0;
	int HessCount = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////First the per profile parameters/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	if(DimsPerProf > 0){
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->TotalProfiles; i++){
			for(int j = 0; j < DimsPerProf; j++){
//				printf("RotateToPrin: %i %i %g \n", i,j,gradPhy[i*DimsPerProf + j]);
				gradPrin[i*DimsPerProf + j] = gradPhy[i*DimsPerProf + j];
			}
		}



		DimCount += ((MNStruct *)GHSglobalcontext)->TotalProfiles*DimsPerProf;
		HessCount += ((MNStruct *)GHSglobalcontext)->TotalProfiles;
	}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Then Prior Parameters////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(globaldims+epochdims+epochpriordims > 0){

		int thissize = globaldims+epochpriordims+epochdims;
		double *oneGrad = new double[thissize];
		double *RGrad = new double[thissize];

		for(int j = 0; j < thissize; j++){
			oneGrad[j] = gradPhy[DimCount + j];
		}
		vector_dgemv(GlobalHessian[HessCount], oneGrad,RGrad, thissize, thissize,'T');
		for(int j = 0; j < thissize; j++){
			gradPrin[DimCount + j] = RGrad[j];
		}
	
		delete[] oneGrad;	
		delete[] RGrad;

		DimCount += thissize;
		HessCount += 1+ ((MNStruct *)GHSglobalcontext)->numProfileEpochs;

	}


}



void rotate2Physical(double *xPrin, double *xPhy){

        int DimsPerProf = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int DimsPerEpoch = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;

        int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*DimsPerProf;
        int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*DimsPerEpoch;



	int DimCount = 0;
	int HessCount = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////First the per profile parameters/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	if(DimsPerProf > 0){

		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->TotalProfiles; i++){
			for(int j = 0; j < DimsPerProf; j++){
		//		printf("RotateToPhys: %i %i %g \n", i,j,xPrin[i*DimsPerProf + j]);
				xPhy[i*DimsPerProf + j] =  xPrin[i*DimsPerProf + j];
			}
		}

		DimCount += ((MNStruct *)GHSglobalcontext)->TotalProfiles*DimsPerProf;
		HessCount += ((MNStruct *)GHSglobalcontext)->TotalProfiles;
	
	}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Then Prior Parameters////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(globaldims+epochdims+epochpriordims > 0){

		int thissize = globaldims+epochpriordims+epochdims;

		double *oneGrad = new double[thissize];
		double *RGrad = new double[thissize];

		for(int j = 0; j < thissize; j++){
			oneGrad[j] = xPrin[DimCount + j];
		}
		vector_dgemv(GlobalHessian[HessCount], oneGrad,RGrad, thissize, thissize,'N');
		for(int j = 0; j < thissize; j++){
			xPhy[DimCount + j] = RGrad[j];
		}
	
		delete[] oneGrad;	
		delete[] RGrad;

		DimCount += thissize;
		HessCount += 1+((MNStruct *)GHSglobalcontext)->numProfileEpochs;

	}

}



void OutputMLFiles(int nParameters, double* pdParameterEstimates, double MLike, int startDim){

        std::ofstream designfile;
	std::string rootname = ((MNStruct *)GHSglobalcontext)->rootName;
        std::string T2Sname = rootname+"T2scaling.txt";

        designfile.open(T2Sname.c_str());
        int numtofit=((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps;
        for(int i =1;i<((MNStruct *)GHSglobalcontext)->numFitTiming; i++){
                designfile << ((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[i][0]].label[((MNStruct *)GHSglobalcontext)->TempoFitNums[i][1]];
                designfile << " ";
                std::stringstream ss;
                ss.precision(std::numeric_limits<long double>::digits);//override the default

                ss << ((MNStruct *)GHSglobalcontext)->LDpriors[i][0];
                ss << " ";
                ss << ((MNStruct *)GHSglobalcontext)->LDpriors[i][1];
                designfile << ss.str();
                designfile << "\n";
        }

        designfile.close();

	int incRED = 0;
	
	long double Tempo2Fit[((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps];
	std::string longname="GHSMLVals";

	std::vector<double> paramlist(2*nParameters);

	double **paramarray = new double*[nParameters];
	for(int p =0;p < nParameters; p++){
		paramarray[p]=new double[4]();
	}


	int pcount = startDim;


	int fitcount = 0;
	Tempo2Fit[fitcount] = 0;

	long double val = pdParameterEstimates[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[0][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[0][0]);
	printf("Phase %g %.20Lg \n", pdParameterEstimates[pcount], val);
	pcount++;
	fitcount++;

	for(int i =1; i< ((MNStruct *)GHSglobalcontext)->numFitTiming; i++){
		long double val = pdParameterEstimates[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[i][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[i][0]);

		((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[i][0]].val[((MNStruct *)GHSglobalcontext)->TempoFitNums[i][1]] = val;
		printf("%s   %g    %.20Lg \n", ((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[i][0]].shortlabel[0], pdParameterEstimates[pcount], val);

		
		
		Tempo2Fit[fitcount] = ((MNStruct *)GHSglobalcontext)->LDpriors[i][0];

		pcount++;
		fitcount++;
	}

	for(int i =0; i<((MNStruct *)GHSglobalcontext)->numFitJumps; i++){
		long double val = pdParameterEstimates[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[fitcount+i][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[fitcount+i][0]);
		((MNStruct *)GHSglobalcontext)->pulse->jumpVal[((MNStruct *)GHSglobalcontext)->TempoJumpNums[i]] = val;
		((MNStruct *)GHSglobalcontext)->pulse->jumpValErr[((MNStruct *)GHSglobalcontext)->TempoJumpNums[i]] = 0;

		Tempo2Fit[fitcount] = ((MNStruct *)GHSglobalcontext)->LDpriors[fitcount+i][0];
		printf("Jump %i %g %.20Lg \n", i, pdParameterEstimates[pcount], val);
		pcount++;
	}

	std::ofstream profilefile;
	std::string dname = "GHSNewProfInfo.dat";

	profilefile.open(dname.c_str());
	profilefile << ((MNStruct *)GHSglobalcontext)->totshapecoeff << " " << ((MNStruct *)GHSglobalcontext)->totalEvoCoeff <<"\n";


	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff;
	profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileShape[0] <<"\n";
	for(int i = 1; i < totshapecoeff; i++){
		double newPVal = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[i] +  pdParameterEstimates[pcount];
		//printf("B: %g \n", ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[c]);
		profilefile << std::setprecision(10) << newPVal <<"\n";
		pcount++;
	}

	for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
		printf("B: %g \n", ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[c]);
		profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[c] <<"\n";
	}
	

	for(int p = 1; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){
		profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][0] <<"\n";	
		for(int i =1; i < totshapecoeff; i++){
			double newPVal = ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i]  + pdParameterEstimates[pcount];
			profilefile << std::setprecision(10) << newPVal <<"\n";
			pcount++;
		}
	}

	profilefile.close();

	double *ProfileScatter;
	if(((MNStruct *)GHSglobalcontext)->incProfileScatter == 1){
		ProfileScatter =  new double[((MNStruct *)GHSglobalcontext)->pulse->nSx];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nSx; i++){
			ProfileScatter[i] = pdParameterEstimates[pcount]; 
			((MNStruct *)GHSglobalcontext)->pulse->param[param_sx].val[i] = ProfileScatter[i];
			pcount++;
		}
	}

	
	printf("   Max Like %.10g \n", MLike);
	TNtextOutput(((MNStruct *)GHSglobalcontext)->pulse, 1, 0, Tempo2Fit,  GHSglobalcontext, 0, nParameters, paramlist, 0.0, 0, 0, 0, longname, paramarray);


}


void readextract(std::string ename){

        int number_of_lines = 0;

        std::ifstream checkfile;
        std::string checkname;
        checkfile.open(ename.c_str());
        std::string line;
        while (getline(checkfile, line))
                ++number_of_lines;

        checkfile.close();

	std::ifstream summaryfile;


	summaryfile.open(ename.c_str());



//	printf("Getting ML \n");
	double maxlike = -1.0*pow(10.0,10);
	for(int i=0;i<number_of_lines;i++){

		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> paramlist(begin,eof);

		int ndim=paramlist.size()-1;

		double like = paramlist[ndim];
		int start = 2*((MNStruct *)GHSglobalcontext)->TotalProfiles;
		if(i > number_of_lines-100 && like > maxlike){
			maxlike = like;
			 for(int i = 0; i < ndim; i++){
				GlobalMaxLikeVec[start+i] = paramlist[i];

                	 }
		}

	}
	summaryfile.close();

}





void NewGHSProfileDomainLike(int* ndim, double* PrinCube, double* likelihood, double* PrinGrad){


	int debug = ((MNStruct *)GHSglobalcontext)->debug; 
	int threads = omp_get_max_threads();


	
	if(debug == 1)printf("In like \n");

	double *grad = new double[*ndim]();
	double *Cube = new double[*ndim];


	rotate2Physical(PrinCube, Cube);

//	for(int i = 0; i < *ndim-1; i++){
//		Cube[i]=GlobalMaxLikeVec[i];
//	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Set Constants///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	double LogTen = log(10.0);
	double TwoPi = 2*M_PI;

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get dimensionality//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


        int perProfDims = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int DimsPerEpoch = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;


	int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*perProfDims;	
	int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*DimsPerEpoch;



	
	int globalparams = globaldims;

	double **threadgrads = new double*[globalparams];
	for(int i = 0; i < globalparams; i++){ threadgrads[i] = new double[threads]();}


	double **PriorThreadGrads = new double*[epochpriordims];
	for(int i = 0; i < epochpriordims; i++){ PriorThreadGrads[i] = new double[threads]();}


	double priorterm = 0;

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the BATS///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = ((MNStruct *)GHSglobalcontext)->pulse->nobs;


	 long double *ProfileBats=new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

		ProfileBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
	      
		ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;

	 }






/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int pcount = 2*TotalProfs+epochdims+epochpriordims;
	int phaseDim = 2*TotalProfs+epochdims+epochpriordims;

	long double phase = Cube[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[0][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[0][0]);	
	phase = phase*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;
	pcount++;

	double PhasePrior = ((MNStruct *)GHSglobalcontext)->PhasePrior;
	double phasePriorTerm = 0.5*Cube[phaseDim]*Cube[phaseDim]/(PhasePrior*PhasePrior);
	double phasePriorGrad = Cube[phaseDim]/(PhasePrior*PhasePrior);


	double *TimingParams = new double[((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps]();
	double *TimingSignal = new double[((MNStruct *)GHSglobalcontext)->pulse->nobs];

	int fitcount = 1;
	for(int i = phaseDim+1; i < phaseDim+((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps; i++){
		TimingParams[fitcount] = Cube[i];
		fitcount++;
		pcount++;
	}

	vector_dgemv(GlobalDMatrixVec, TimingParams,TimingSignal,((MNStruct *)GHSglobalcontext)->pulse->nobs,((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps,'N');

	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){
		//if(i==0){printf("MB: %.15Lg %.15Lg %.15Lg\n", ModelBats[i] , ModelBats[i] - phase, phase);}
     		 ModelBats[i] -=  (long double) TimingSignal[i]/SECDAY + phase;
	}	

	delete[] TimingSignal;
	delete[] TimingParams;	



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Profile Parameters//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int NFBasis = ((MNStruct *)GHSglobalcontext)->NFBasis;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
	}


	double **ProfCoeffs=new double*[1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){ProfCoeffs[i] = new double[totshapecoeff]();}

	ProfCoeffs[0][0] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[0];
	for(int i = 1; i < totshapecoeff; i++){
		ProfCoeffs[0][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[i] +  Cube[pcount];
		pcount++;
	}

	for(int p = 1; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){	
		for(int i =1; i < totshapecoeff; i++){
			ProfCoeffs[p][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i]  + Cube[pcount];
			//printf("Eco: %i %g %g \n", i, ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i]  , Cube[pcount]);
			pcount++;
		}
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Scattering//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	double *ProfileScatter;
	if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){

		int sthreaddim = ((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps+(((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1)*(totshapecoeff-1);
		ProfileScatter =  new double[((MNStruct *)GHSglobalcontext)->pulse->nSx];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nSx; i++){
			ProfileScatter[i] =  pow(10.0, Cube[pcount]); 
			
			if(Cube[pcount] < -10){
                        	priorterm -= 2*log(10.0)*(Cube[pcount]+10);

				for(int i = 0; i < threads; i++){
                                        threadgrads[sthreaddim + i][i] -= log(10.0)/threads;
                                }
			}
			pcount++;
		}
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;


	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;

	
	double *eplikes = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();

	#pragma omp parallel for 
	for(int ep = minep; ep < maxep; ep++){

		double *TotalProfCoeffs = new double[totshapecoeff];
		int thisthread = omp_get_thread_num();


		int t = 0;

		for(int sep = 0; sep < ep; sep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[sep];
		}
		

		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];


		double EpochChisq = 0;	
		double EpochdetN = 0;
		double EpochLike = 0;


		int EDimcount = 2*TotalProfs + epochpriordims + ep*DimsPerEpoch;
		int OriginalEDimCount = EDimcount;



		for(int ch = 0; ch < NChanInEpoch; ch++){

			int nTOA = t;
			

			double ProfileAmp = Cube[TotalProfs + nTOA];
			double ProfileSigma = Cube[nTOA];

			//printf("Amp and Noise: %i %g %g \n", nTOA, ProfileAmp, ProfileSigma);

			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];


			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;


			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;

			double *RolledData  = new double[2*NFBasis]();
			double *ProfileVec  = new double[2*NFBasis]();
			double *ProfGradVec = new double[2*NFBasis]();




			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get and modify binpos///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////
		    
			long double binpos = ModelBats[nTOA]; 


			long double binposSec = (ProfileBats[nTOA] - ModelBats[nTOA])*SECDAY;
			double DFP = double(FoldingPeriod);
			long double WrapBinPos = -FoldingPeriod/2 + fmod(FoldingPeriod + fmod(binposSec+FoldingPeriod/2, FoldingPeriod), FoldingPeriod);

			long double maxbinvalue = (FoldingPeriod/Nbins);
			double DNewInterpBin = double(fmod(maxbinvalue + fmod(-WrapBinPos, maxbinvalue), maxbinvalue)/((MNStruct *)GHSglobalcontext)->InterpolatedTime);
			int InterpBin = floor(DNewInterpBin+0.5);
			InterpBin = InterpBin%((MNStruct *)GHSglobalcontext)->NumToInterpolate;



			long double NewWholeBinTime = -WrapBinPos-InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime;
			int RollBin = int(round(NewWholeBinTime/maxbinvalue));
			int WrapRollBin = ((Nbins) + -RollBin % (Nbins)) % (Nbins);


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].freq - reffreq)/1000.0;


			for(int i =0; i < totshapecoeff; i++){
				TotalProfCoeffs[i]=0;	
			}				

			int cpos = 0;
	
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				for(int p = 0; p < 1 + ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i = 0; i < numcoeff[c]; i++){
						TotalProfCoeffs[i+cpos] += ProfCoeffs[p][i+cpos]*pow(freqdiff, p);						
					}
				}

				cpos += numcoeff[c];
			}




			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[InterpBin], TotalProfCoeffs, ProfileVec,2*NFBasis,totshapecoeff,'N');
			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFJitterBasis[InterpBin], TotalProfCoeffs, ProfGradVec,2*NFBasis,totshapecoeff,'N');         



		///////////////////////////////////////////Rotate Data Vector////////////////////////////////////////////////////////////





			for(int j = 0; j < NFBasis; j++){

				double freq = (j+1.0)/ProfNbins;
				double theta = 2*M_PI*WrapRollBin*freq;
				double RealShift = cos(theta);
				double ImagShift = sin(-theta);

				double RData = ((MNStruct *)GHSglobalcontext)->FourierProfileData[t][1+j];
				double IData = ((MNStruct *)GHSglobalcontext)->FourierProfileData[t][1+j+ProfNbins/2+1];

				double RealRolledData = RData*RealShift - IData*ImagShift;
				double ImagRolledData = RData*ImagShift + IData*RealShift;

				RolledData[j] = RealRolledData;
				RolledData[j+NFBasis] = ImagRolledData;

			}


			double *OneInterpMatrix = new double[totshapecoeff*2*NFBasis]();
			std::memcpy(OneInterpMatrix, ((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[InterpBin], totshapecoeff*2*NFBasis*sizeof(double));

			double ScatterGrad = 0;
			if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){



	
				int Sindex = ((MNStruct *)GHSglobalcontext)->ScatterIndex[nTOA];

				double ISS = 1.0/(pow( ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB, 4)/pow(10.0, 9.0*4.0));
				double tau = ProfileScatter[Sindex];
				double STime = tau*ISS;


				for(int j = 0; j < NFBasis; j++){


					//OrigProf[j] = ProfileVec[j];
					//OrigProf[j+NFBasis] = ProfileVec[j+NFBasis];
					double f = (j+1.0)/((MNStruct *)GHSglobalcontext)->ReferencePeriod;
					double w = TwoPi*f;
					double RConv = 1.0/(w*w*STime*STime+1); 
					double IConv = -w*STime/(w*w*STime*STime+1);

					double RProf = ProfileVec[j];
					double IProf = ProfileVec[j+NFBasis];


//					double ScatterRealProf =  RProf*RConv - IProf*IConv;
//					double ScatterImagProf =  RProf*IConv + IProf*RConv;

					ProfileVec[j] = RProf*RConv - IProf*IConv; //ScatterRealProf;
					ProfileVec[j+NFBasis] = RProf*IConv + IProf*RConv; //ScatterImagProf;


					//Get the Gradient while you are here

					double ROrigProf = ProfileAmp*RProf;
					double IOrigProf = ProfileAmp*IProf;

					double GradDenom = 1.0/pow(1.0 + STime*STime*w*w, 2);

					double RealFunc =  RolledData[j] - ROrigProf*RConv + IOrigProf*IConv;
					double ImagFunc = RolledData[j+NFBasis] - ROrigProf*IConv - IOrigProf*RConv;

					

					double RealGrad = 2*STime*STime*w*w*LogTen*GradDenom*ROrigProf + STime*w*(STime*STime*w*w - 1)*LogTen*GradDenom*IOrigProf;
					double ImagGrad = 2*STime*STime*w*w*LogTen*GradDenom*IOrigProf - STime*w*(STime*STime*w*w - 1)*LogTen*GradDenom*ROrigProf;

					ScatterGrad += (RealGrad*RealFunc + ImagGrad*ImagFunc)/ProfileSigma/ProfileSigma;


					RProf = ProfGradVec[j];
					IProf = ProfGradVec[j+NFBasis];

//					ScatterRealProf =  RProf*RConv - IProf*IConv;
//					ScatterImagProf =  RProf*IConv + IProf*RConv;

					ProfGradVec[j] = RProf*RConv - IProf*IConv; //ScatterRealProf;
					ProfGradVec[j+NFBasis] = RProf*IConv + IProf*RConv; //ScatterImagProf;

					for(int k = 0; k < totshapecoeff; k++){


						RProf = OneInterpMatrix[j+k*2*NFBasis];
                                                IProf = OneInterpMatrix[(j+NFBasis) + k*2*NFBasis];

  //                                              ScatterRealProf =  RProf*RConv - IProf*IConv;
//                                                ScatterImagProf =  RProf*IConv + IProf*RConv;

                                                OneInterpMatrix[j+k*2*NFBasis] = RProf*RConv - IProf*IConv; //ScatterRealProf;
                                                OneInterpMatrix[(j+NFBasis) + k*2*NFBasis] = RProf*IConv + IProf*RConv; //ScatterImagProf;
					}

				}

				/*for(int k = 0; k < totshapecoeff; k++){
					for(int j = 0; j < NFBasis; j++){


						double f = (j+1.0)/((MNStruct *)GHSglobalcontext)->ReferencePeriod;
						double w = 2.0*M_PI*f;
						double RConv = 1.0/(w*w*STime*STime+1); 
						double IConv = -w*STime/(w*w*STime*STime+1);

						double RProf = OneInterpMatrix[j+k*2*NFBasis];
						double IProf = OneInterpMatrix[(j+NFBasis) + k*2*NFBasis];

						double ScatterRealProf =  RProf*RConv - IProf*IConv;
						double ScatterImagProf =  RProf*IConv + IProf*RConv;

						OneInterpMatrix[j+k*2*NFBasis] = ScatterRealProf;
						OneInterpMatrix[(j+NFBasis) + k*2*NFBasis] = ScatterImagProf;

					}
				}*/


			}




		///////////////////////////////////////////Calculate Likelihood/////////////////////////////////////////////////

			      
			double Chisq = 0;
			double detN = 0;
			double profilelike=0;
	

			double noise = ProfileSigma*ProfileSigma;
			double OffPulsedetN = log(noise);
			detN = 2*NFBasis*OffPulsedetN;
			noise = 1.0/noise;

			double *NResVec = new double[2*NFBasis]();
		

			double BaseGrad = 0;
			double AmpGrad = 0;

			for(int j = 0; j < 2*NFBasis; j++){
				//if(t==0){printf("Like %i %i %g %g %g %g %g %g\n", nTOA, j, RolledData[j], ProfileAmp*ProfileVec[j], Chisq, ProfileAmp, ProfileVec[j], ProfileSigma );}
				double pres = RolledData[j] - ProfileAmp*ProfileVec[j];
				NResVec[j] = pres*noise;

				Chisq += pres*NResVec[j];
				AmpGrad -= ProfileVec[j]*NResVec[j];
			}
			//sleep(1);

			double NoiseGrad = (-Chisq+2*NFBasis)/ProfileSigma;

			grad[TotalProfs+nTOA] = AmpGrad;
			grad[nTOA] = NoiseGrad;

			//if(t==0){printf("AmpGrads: %g %g \n", AmpGrad, NoiseGrad);}


			int writeprofiles = 0;

			if(writeprofiles == 1){
				std::ofstream profilefile;
				char toanum[15];
				sprintf(toanum, "%d", nTOA);
				std::string ProfileName =  toanum;
				char* longname =    ((MNStruct *)GHSglobalcontext)->rootName;
				std::string dname = longname+ProfileName+"-Profile.txt";

				double *shapevec = new double[Nbins];

				fftw_plan plan;
				fftw_complex *FFTData;
				FFTData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nbins/2+1));

				plan = fftw_plan_dft_c2r_1d(Nbins, FFTData, shapevec, FFTW_ESTIMATE);	

				for(int b = 0; b < Nbins/2+1; b++){
					FFTData[b][0] = 0;
					FFTData[b][1] = 0;
				}
			       for(int b = 1; b < NFBasis+1; b++){
					FFTData[b][0] = ProfileVec[b-1];
					FFTData[b][1] = ProfileVec[b+NFBasis-1];
				}	

				fftw_execute(plan);

				profilefile.open(dname.c_str());
				double MLChisq = 0;
				profilefile << "#" << std::setprecision(20) << ((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].bat << " " << ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0] << " " << ((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].freq <<"\n";

				for(int b = 0; b < Nbins; b++){

					profilefile << b << " " << ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][b][1]/ProfileAmp  <<  " " << shapevec[b] << "\n";

				}

				delete[] shapevec;
				fftw_free(FFTData);
				fftw_destroy_plan ( plan );

				profilefile.close();
			}

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Get EpochGrads//////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int B2C = 0;
			double PhaseGrad = 0;
			for(int j = 0; j < 2*NFBasis; j++){
//				if(t==0){printf("NRV %i %g %g \n",j,NResVec[j], ProfGradVec[j]);}
				PhaseGrad += -1*NResVec[j]*ProfGradVec[j]*ProfileAmp;		
			}
			
	
			threadgrads[B2C][thisthread] += double(((MNStruct *)GHSglobalcontext)->ReferencePeriod)*PhaseGrad;
			B2C++;


			for(int j = 1; j < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; j++){	
				threadgrads[B2C][thisthread] += GlobalDMatrix[nTOA][j]*PhaseGrad;
				B2C++;

			}


			double *GradResVec = new double[totshapecoeff];
	
			vector_dgemv(OneInterpMatrix, NResVec, GradResVec, 2*NFBasis, totshapecoeff, 'T');
			

			for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){

				double fval = -1*ProfileAmp*pow(freqdiff, p);

				for(int c = 1; c < totshapecoeff; c++){
			//		if(t==0){printf("GRV: %i %i %g \n", p,c,GradResVec[c]);}
					threadgrads[B2C + p*(totshapecoeff-1) + c-1][thisthread] += fval*GradResVec[c];
					
				}
			}
			B2C += (((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1)*(totshapecoeff-1);		

			delete[] GradResVec;

			if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){
				threadgrads[B2C + ((MNStruct *)GHSglobalcontext)->ScatterIndex[nTOA]][thisthread] +=  ScatterGrad;
				B2C += 1;
			}
		

			profilelike = -0.5*(detN + Chisq);



			EpochChisq+=Chisq;
			EpochdetN+=detN;
			EpochLike+=profilelike;

			delete[] ProfileVec;
			delete[] ProfGradVec;
			delete[] RolledData;
			delete[] NResVec;


			delete[] OneInterpMatrix;
			


			t++;
		}



///////////////////////////////////////////


		if(((MNStruct *)GHSglobalcontext)->incWideBandNoise == 0){
			eplikes[ep] += EpochLike;
			lnew += EpochLike;
		}
		delete[] TotalProfCoeffs;

////////////////////////////////////////////
	
	}


	





	double omplnew = 0;
	for(int ep = minep; ep < maxep; ep++){
		omplnew += eplikes[ep];
	}
	delete[] eplikes;

	for(int i = 0; i < epochpriordims; i++){
		double gsum = 0;
                for(int j = 0; j < threads; j++){
                        gsum += PriorThreadGrads[i][j];
                }
		 grad[2*TotalProfs + i] = gsum;
	}

	for(int i = 0; i < globalparams; i++){
		double gsum = 0;
		for(int j = 0; j < threads; j++){
			gsum += threadgrads[i][j];
		}
	 	// if(i==globalparams-1){printf("global grad: %i %g %g %g \n", i, Cube[phaseDim+i], gsum, -1*omplnew + phasePriorTerm);}
		 grad[2*TotalProfs + epochdims + epochpriordims + i] = gsum;
		
	}



	grad[2*TotalProfs + epochdims + epochpriordims]+=phasePriorGrad;

//	for(int i = 0; i < *ndim; i++){
//		printf("par grad: %i %g %g %g \n", i, Cube[i], grad[i], -1*omplnew + phasePriorTerm);
//	}


	

	if(((MNStruct *)GHSglobalcontext)->incProfileScatter > 0){
		delete[] ProfileScatter;
	}


	for(int j = 0; j < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] ProfCoeffs[j];
	}
	delete[] ProfCoeffs;
	delete[] numcoeff;

	delete[] ProfileBats;
	delete[] ModelBats;


	for(int i = 0; i < globalparams; i++){
		delete[] threadgrads[i];
	}
	delete[] threadgrads;


	for(int i = 0; i < epochpriordims; i++){
		delete[] PriorThreadGrads[i];
	}
	delete[] PriorThreadGrads;








	double finallikelihood = -1*omplnew + phasePriorTerm + 0.5*priorterm;
	*likelihood=finallikelihood;

	if(-1*finallikelihood > GlobalMaxLike && ((MNStruct *)GHSglobalcontext)->WriteNewML == 1){
		GlobalMaxLike = -1*finallikelihood;
		printf("Like: %.16g \n", -1*finallikelihood);
		write_newMLLike(*ndim, Cube, GlobalMaxLike);
		if(debug == 1)printf("Written new ML file\n");
	
	}

//	printf("Like: %.16g \n", -1*finallikelihood);

	

	rotate2Principal(PrinGrad, grad);


	delete[] Cube;
	delete[] grad; 


	
}



void GHSWriteProf(std::string ename){


        int number_of_lines = 0;

        std::ifstream checkfile;
        std::string checkname;
        checkfile.open(ename.c_str());
        std::string line;
        while (getline(checkfile, line))
                ++number_of_lines;

        checkfile.close();

	std::ifstream summaryfile;


	summaryfile.open(ename.c_str());


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get dimensionality//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


        int perProfDims = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int DimsPerEpoch = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;


	int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*perProfDims;	
	int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*DimsPerEpoch;
	
	int globalparams = globaldims;

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	double chanfreq[3];
	chanfreq[0] = 100;
	chanfreq[1] = 145;
	chanfreq[2] = 190;


	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int Nbins = GlobalNBins;


	double **AverageProfile = new double*[3];
	double **AverageErrs = new double*[3];

	for(int i = 0; i < 3; i++){
		AverageProfile[i]= new double[Nbins]();
		AverageErrs[i] = new double[Nbins]();
	}

	for(int i=0;i<number_of_lines;i++){

		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> Cube(begin,eof);

		int ndim=Cube.size()-1;

	/////////////////////////////////////////////////////////////////////////////////////////////  
	/////////////////////////Get Profile Parameters//////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


		int NFBasis = ((MNStruct *)GHSglobalcontext)->NFBasis;
		int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

		int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
			numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
		}


		double **ProfCoeffs=new double*[1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
		for(int i = 0; i < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){ProfCoeffs[i] = new double[totshapecoeff]();}

		ProfCoeffs[0][0] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[0];
		int pcount = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;
		for(int i = 1; i < totshapecoeff; i++){
			ProfCoeffs[0][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[i] +  Cube[pcount];
			pcount++;
		}

		for(int p = 1; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){	
			for(int i =1; i < totshapecoeff; i++){
				ProfCoeffs[p][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i]  + Cube[pcount];
				pcount++;
			}
		}


		double *ProfileVec  = new double[2*NFBasis]();
		double *shapevec = new double[Nbins];

		fftw_plan plan;
		fftw_complex *FFTData;
		FFTData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nbins/2+1));

		plan = fftw_plan_dft_c2r_1d(Nbins, FFTData, shapevec, FFTW_ESTIMATE);

			
			
			




		for(int f = 0; f < 3; f++){


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes//////////////////////////


			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (chanfreq[f] - reffreq)/1000.0;

			double *TotalProfCoeffs = new double[totshapecoeff];
			for(int i =0; i < totshapecoeff; i++){
				TotalProfCoeffs[i]=0;	
			}				

			int cpos = 0;

			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				for(int p = 0; p < 1 + ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i = 0; i < numcoeff[c]; i++){
						TotalProfCoeffs[i+cpos] += ProfCoeffs[p][i+cpos]*pow(freqdiff, p);						
					}
				}

				cpos += numcoeff[c];
			}




			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[0], TotalProfCoeffs, ProfileVec,2*NFBasis,totshapecoeff,'N');

			for(int j = 0; j < NFBasis; j++){

				double freq = (j+1.0)/Nbins;
				double theta = 2*M_PI*(Nbins/2)*freq;
				double RealShift = cos(theta);
				double ImagShift = sin(-theta);

				double RData = ProfileVec[j];
				double IData = ProfileVec[j+NFBasis];

				double RealRolledData = RData*RealShift - IData*ImagShift;
				double ImagRolledData = RData*ImagShift + IData*RealShift;

				ProfileVec[j] = RealRolledData;
				ProfileVec[j+NFBasis] = ImagRolledData;

			}


			for(int b = 0; b < Nbins/2+1; b++){
				FFTData[b][0] = 0;
				FFTData[b][1] = 0;
			}
		       for(int b = 1; b < NFBasis+1; b++){
				FFTData[b][0] = ProfileVec[b-1];
				FFTData[b][1] = ProfileVec[b+NFBasis-1];
			}
		       



			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Fill Arrays with interpolated state//////////////////////////////////////////////////////



			fftw_execute(plan);

			for(int j = 0; j < Nbins; j++){
				AverageProfile[f][j] += shapevec[j];
			}

			delete[] TotalProfCoeffs;
		}


		for(int j = 0; j < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
		    delete[] ProfCoeffs[j];
		}
		delete[] ProfCoeffs;
		delete[] numcoeff;
		delete[] ProfileVec;
		delete[] shapevec;
		fftw_free(FFTData);
		fftw_destroy_plan ( plan );
	
		
	}

	for(int f = 0; f < 3; f++){
		for(int j = 0; j < Nbins; j++){
			AverageProfile[f][j] /= number_of_lines;

	//		printf("Mean: %i %i %g \n", f, j, AverageProfile[f][j]);
		}
	}

	summaryfile.close();
	summaryfile.open(ename.c_str());

	for(int i=0;i<number_of_lines;i++){

		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> Cube(begin,eof);

		int ndim=Cube.size()-1;

	/////////////////////////////////////////////////////////////////////////////////////////////  
	/////////////////////////Get Profile Parameters//////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


		int NFBasis = ((MNStruct *)GHSglobalcontext)->NFBasis;
		int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

		int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
			numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
		}


		double **ProfCoeffs=new double*[1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
		for(int i = 0; i < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){ProfCoeffs[i] = new double[totshapecoeff]();}

		ProfCoeffs[0][0] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[0];
		int pcount = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;
		for(int i = 1; i < totshapecoeff; i++){
			ProfCoeffs[0][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileShape[i] +  Cube[pcount];
			pcount++;
		}

		for(int p = 1; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){	
			for(int i =1; i < totshapecoeff; i++){
				ProfCoeffs[p][i] = ((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p-1][i]  + Cube[pcount];
				pcount++;
			}
		}


		double *ProfileVec  = new double[2*NFBasis]();
		double *shapevec = new double[Nbins];

		fftw_plan plan;
		fftw_complex *FFTData;
		FFTData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nbins/2+1));

		plan = fftw_plan_dft_c2r_1d(Nbins, FFTData, shapevec, FFTW_ESTIMATE);


		for(int f = 0; f < 3; f++){



			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes//////////////////////////


			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (chanfreq[f] - reffreq)/1000.0;

			double *TotalProfCoeffs = new double[totshapecoeff];
			for(int i =0; i < totshapecoeff; i++){
				TotalProfCoeffs[i]=0;	
			}				

			int cpos = 0;

			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				for(int p = 0; p < 1 + ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i = 0; i < numcoeff[c]; i++){
						TotalProfCoeffs[i+cpos] += ProfCoeffs[p][i+cpos]*pow(freqdiff, p);						
					}
				}

				cpos += numcoeff[c];
			}




			vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedFBasis[0], TotalProfCoeffs, ProfileVec,2*NFBasis,totshapecoeff,'N');

			for(int j = 0; j < NFBasis; j++){

				double freq = (j+1.0)/Nbins;
				double theta = 2*M_PI*(Nbins/2)*freq;
				double RealShift = cos(theta);
				double ImagShift = sin(-theta);

				double RData = ProfileVec[j];
				double IData = ProfileVec[j+NFBasis];

				double RealRolledData = RData*RealShift - IData*ImagShift;
				double ImagRolledData = RData*ImagShift + IData*RealShift;

				ProfileVec[j] = RealRolledData;
				ProfileVec[j+NFBasis] = ImagRolledData;

			}
		       

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Fill Arrays with interpolated state//////////////////////////////////////////////////////

			for(int b = 0; b < Nbins/2+1; b++){
				FFTData[b][0] = 0;
				FFTData[b][1] = 0;
			}
		       for(int b = 1; b < NFBasis+1; b++){
				FFTData[b][0] = ProfileVec[b-1];
				FFTData[b][1] = ProfileVec[b+NFBasis-1];
			}

			fftw_execute(plan);

			for(int j = 0; j < Nbins; j++){
				AverageErrs[f][j] += pow(shapevec[j]-AverageProfile[f][j],2);
			}

			delete[] TotalProfCoeffs;

		}

		for(int j = 0; j < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
		    delete[] ProfCoeffs[j];
		}
		delete[] ProfCoeffs;
		delete[] numcoeff;
		delete[] ProfileVec;
		delete[] shapevec;
		fftw_free(FFTData);
		fftw_destroy_plan ( plan );
	
		
	}

	for(int f = 0; f < 3; f++){
		for(int j = 0; j < Nbins; j++){
			AverageErrs[f][j] = sqrt(AverageErrs[f][j]/number_of_lines);

			///printf("Mean: %i %i %g %g\n", f, j, AverageProfile[f][j], AverageErrs[f][j]);
		}
	}


	for(int f = 0; f < 3; f++){
		double max = 0;
		for(int j = 0; j < Nbins; j++){
			if(AverageProfile[f][j] > max){max = AverageProfile[f][j];}
		}
		for(int j = 0; j < Nbins; j++){
			AverageErrs[f][j] /= max;
			AverageProfile[f][j] /= max;
		}

	}

	for(int j = 0; j < Nbins; j++){

		printf("Mean: %i %g %g %g %g %g %g \n", j, AverageProfile[0][j], AverageErrs[0][j], AverageProfile[1][j], AverageErrs[1][j], AverageProfile[2][j], AverageErrs[2][j]);
	}
	
}
