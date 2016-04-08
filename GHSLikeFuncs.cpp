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


FILE* test_gauss_outfile;
FILE* test_gauss_goutfile;

void write_gauss_ghs_extract(int* ndim,double* x,double* val,double* g);
void write_gauss_ghs_extract_with_logpostval(int* ndim,double* x,double* logpostval,double* g);
void write_gauss_ghs_extract_with_logpostval_and_grad(int* ndim,double* x,double* logpostval,double* g);
void nd_uncorr_gauss_neg_log_post(int* ndim,double* x,double* v,double* g);

void GHSProfileDomainLike(int* ndim,double* x,double* v,double* g);
void GHSProfileDomainLike2(int* ndim,double* x,double &v,double* g);
void GetMaxAmps(int ndim, double *MaxAmps);
void GetMaxSteps(int ndim, double *MaxAmps, double *StepSize, double **SSM);

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

void assignGHScontext(void *context){
        GHSglobalcontext=context;
}

int GHSWrap(int kX, int const kLowerBound, int const kUpperBound)
{
    int range_size = kUpperBound - kLowerBound + 1;

    if (kX < kLowerBound)
        kX += range_size * ((kLowerBound - kX) / range_size + 1);

    return kLowerBound + (kX - kLowerBound) % range_size;
}

void callGHS(int NBurn, int NSamp, int GHSresume){


	int id;



	double* StartPoint;
	double* StepSize;
	double **SSM;
	void (*nlp)(int*,double*,double*,double*);
	void (*wrt_ext)(int*,double*,double*,double*);
	double ScaleFactor;
//	std::string temproot = ((MNStruct *)GHSglobalcontext)->rootName;
//	char *chartroot = new char[temproot.length() + 1];
//	std::strcpy(chartroot, temproot.c_str());
	char* fl_pfx =    ((MNStruct *)GHSglobalcontext)->rootName; //chartroot; //"GHSResults/102050cm-FitProfile-EPS-";
	int seed;
	int fb_int;
	int max_stp;
	int resume;
	char ext_file_name[128];
	char gext_file_name[128];
	int nburn=NBurn;
	int nsamp=NSamp;
	/*
	Note that if nsamp>0, the sampler will be stopped forcefully after nsamp samples
	have been taken. This does not mean that the algorithm will have converged.
	Check the Hanson values for convergence.
	*/
	
	printf("\n=================================================\n");
	printf("|    Example program in C: 1                    |\n");
	printf("|    Uncorrelated N-D Gaussian                  |"); 
	printf("\n=================================================\n");


	((MNStruct *)GHSglobalcontext)->ProfEvoTimeRef = 56100;

	((MNStruct *)GHSglobalcontext)->FitProfParams = 1;
	((MNStruct *)GHSglobalcontext)->diagonalGHS = 0; 

	if(((MNStruct *)GHSglobalcontext)->doMax == 2){
		((MNStruct *)GHSglobalcontext)->diagonalGHS = 1;
	}
	((MNStruct *)GHSglobalcontext)->WriteNewML = 1;
	int readMLFile = 0;

	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = 0;

	for(int ep = 0; ep < NEpochs; ep++){ TotalProfs +=  ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep]; } 

	((MNStruct *)GHSglobalcontext)->TotalProfiles = TotalProfs;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////////////////////////////Get dimensionality/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int perProfDims = 3;
	int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*perProfDims;

	int perEpochDims = ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff;
	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){ perEpochDims++; }

	int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*perEpochDims;

	int epochpriordims = ((MNStruct *)GHSglobalcontext)->numFitEQUAD + ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff;

	int globaldims = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly);

	if(((MNStruct *)GHSglobalcontext)->incDM==5){
		epochpriordims += 2;
		globaldims += 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;
	}


        if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1){
                globaldims++;
        }
        if(((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
                globaldims+=((MNStruct *)GHSglobalcontext)->incWidthEvoTime;
        }



	((MNStruct *)GHSglobalcontext)->GHSperProfDims = perProfDims;
	((MNStruct *)GHSglobalcontext)->GHSperEpochDims = perEpochDims;
	((MNStruct *)GHSglobalcontext)->GHSepochpriordims = epochpriordims;
	((MNStruct *)GHSglobalcontext)->GHSglobaldims = globaldims;

	int ProfileDims = profdims + ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly);
	int TimingDims = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;

	const int ndim = profdims+epochdims+epochpriordims+globaldims;

	double **PhysDMatrix = new double*[TotalProfs];
	for(int i = 0; i < TotalProfs; i++){ PhysDMatrix[i] = new double[TimingDims];}
	double *PhysDMatrixVec = new double[TotalProfs*TimingDims];

	int *TimingGradientSigns = new int[TimingDims]();

	getPhysDVector(GHSglobalcontext, PhysDMatrix, TotalProfs,TimingGradientSigns);
	((MNStruct *)GHSglobalcontext)->TimingGradientSigns=TimingGradientSigns;
//	UpdatePhysDVector(GHSglobalcontext, PhysDMatrix, TotalProfs);


	GlobalDMatrix = PhysDMatrix;
	

	for(int i = 0; i < TotalProfs; i++){
		for(int j = 0; j < TimingDims; j++){
			PhysDMatrixVec[i+j*TotalProfs] = PhysDMatrix[i][j];
		}
	}
	GlobalDMatrixVec = PhysDMatrixVec;
	
	printf("dimensionality is: %i %i %i %i %i %i \n", profdims,epochdims,epochpriordims,globaldims, ndim, TotalProfs);	
	
/*	
  ! start point and step sizes for each parameter
  ! We need to start from a point as close to the 
  ! peak as possible. Here it is (0,0,...,0).
  ! Step-sizes should be approximately the width of
  ! the posterior distribution of the parameter.
  ! Here it is (1,1,...,1).
*/


	//int totalCoeffForMult =  ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff;
	


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


	GlobalParmSet = new double[epochpriordims+epochdims+globaldims]();
	double GPSval = 100;
/**	GlobalParmSet[0]  = GPSval/8117; 2785;*4589.0;
	GlobalParmSet[1]  = GPSval/342
	GlobalParmSet[6]  = GPSval;
	GlobalParmSet[8]  = GPSval;
	GlobalParmSet[9]  = GPSval;
	GlobalParmSet[10]  = GPSval;
	GlobalParmSet[11]  = GPSval;
	GlobalParmSet[12]  = GPSval;
	GlobalParmSet[13]  = GPSval;
	GlobalParmSet[14]  = GPSval;
	GlobalParmSet[15]  = GPSval;
	GlobalParmSet[16]  = GPSval;
	GlobalParmSet[17]  = GPSval;
	GlobalParmSet[18]  = GPSval;
	GlobalParmSet[19]  = GPSval;
	GlobalParmSet[20]  = GPSval;
	GlobalParmSet[21]  = GPSval;
	GlobalParmSet[22]  = GPSval;
	GlobalParmSet[23]  = GPSval;
	GlobalParmSet[24]  = GPSval;

/*
	GlobalParmSet[1]  = 5;
	GlobalParmSet[6]  = 5;
	GlobalParmSet[8]  = 6;
	GlobalParmSet[9]  = 26;
	GlobalParmSet[10]  = 7;
	GlobalParmSet[11]  = 5;
	GlobalParmSet[12]  = 42;
	GlobalParmSet[13]  = 18;
	GlobalParmSet[14]  = 7;
	GlobalParmSet[15]  = 10;
	GlobalParmSet[16]  = 7;
	GlobalParmSet[17]  = 47;
	GlobalParmSet[18]  = 8;
	GlobalParmSet[19]  = 5;
	GlobalParmSet[20]  = 7;
	GlobalParmSet[21]  = 32;
	GlobalParmSet[22]  = 8;
	GlobalParmSet[23]  = 7;
	GlobalParmSet[24]  = 5;
*/

	printf("Getting Max and StepSize %i \n", ndim);
	GetMaxAmps(ndim, StartPoint);

	

	if(readMLFile == 1){
		FILE* MLoutfile;
		char* ml_pfx="NewML.dat";
		MLoutfile=fopen(ml_pfx,"r");

		double* FileStartPoint = new double[ndim]();
		
		for(int i=0;i<ndim;++i){
			fscanf(MLoutfile, "%lf",&FileStartPoint[i]);
			//printf("New start point from file: %i %.15g %.15g\n", i, StartPoint[i], FileStartPoint[i]);
			 StartPoint[i]=FileStartPoint[i];
		}
		fclose(MLoutfile);
		delete[] FileStartPoint;
	}

	GetMaxSteps(ndim, StartPoint, StepSize, SSM);
	GlobalMaxLike = -1*pow(10, 100);
	
	GlobalStartPoint = StartPoint;
	GlobalHessian = SSM;



	if(((MNStruct *)GHSglobalcontext)->FitProfParams == 0){
	//	ndim -= 3*((MNStruct *)GHSglobalcontext)->TotalProfiles; 
		profdims = 0;
	}
	double *PrincipalStartPoint = new double[ndim]();
	
	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Principal(PrincipalStartPoint, StartPoint);
	}
	else{
		for(int i = 0; i < ndim; i++){
			PrincipalStartPoint[i] = StartPoint[i];
		}
	}
//	UpdatePhysDVector(GHSglobalcontext, GlobalDMatrix, TotalProfs);

//	getPhysDVector(GHSglobalcontext, PhysDMatrix, TotalProfs,TimingGradientSigns);
//	GlobalDMatrix = PhysDMatrix;

//	for(int i = 0; i < TotalProfs; i++){
//		for(int j = 0; j < TimingDims; j++){
//			PhysDMatrixVec[i+j*TotalProfs] = PhysDMatrix[i][j];
//		}
//	}
//	GlobalDMatrixVec = PhysDMatrixVec;


/*
      ! find its eigen values (PCA)
      ! the feature matrix is a matrix of eigen values
      ! allows us to sample in principal coordinates
      gauss_featureMatrix=gauss_inverCovMatrix
      call SymPosDefEigenDecomp(gauss_featureMatrix,gauss_eigenValues,gauss_numDims)

      ! assign step size as the inverse square root of eigen values
      ! inverse square root of the eigen values gives the width of the
      ! distribution in the principal coordinates
*/
	for(int i = 0; i < ndim; i++){
		StepSize[i] = 1.0/sqrt(StepSize[i]);
		if(i>=3*TotalProfs){printf("Final Step Size %i %g %g %G \n", i, StartPoint[i], StepSize[i], PrincipalStartPoint[i]);}
	}
	GlobalStepSize = new double[ndim];
	GlobalStepSize = StepSize;
	
	//return;
/*
      ! rotate the start point to the principal coordinates
      ! we need to start from the peak in the principal coordinates
      call rotate2Principal(gauss_numDims,gauss_startPoint)	
*/
	
	
	nlp=&GHSProfileDomainLike;
	/*wrt_ext=&write_gauss_ghs_extract;*/
	wrt_ext=&write_gauss_ghs_extract_with_logpostval_and_grad; 
	

/*  dimensionality scaling factor (a value between 0 and 1) ideal choise is the one which gives you ~68% acceptance rate */
	ScaleFactor = 0.4;


/*  feed back to console interval */

	fb_int=10;

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
	


	if(((MNStruct *)GHSglobalcontext)->doMax == 2){

		double GradML;
		lbfgs_parameter_t param;
		lbfgs_parameter_init(&param);
		double *start = new double[ndim]();
		printf("dims: %i \n", ndim);

		for(int i = 0; i < profdims + epochpriordims+epochdims; i++){
			PrincipalStartPoint[i] = 0;
		}
		GlobalStartPoint[profdims + epochpriordims] = 0;
		int maxruns=3;
		int currentrun = 0;
		int ret = -1;
	
		while(currentrun < maxruns && ret != 0){




			if(currentrun !=0){
				for(int i = 0; i < ndim; i++){
					if(i>=profdims){printf("New Start: %i %g %g %g \n", i, PrincipalStartPoint[i], GlobalMaxLikeVec[i], GlobalStartPoint[i]);}
					PrincipalStartPoint[i] = (GlobalMaxLikeVec[i] - GlobalStartPoint[i])/GlobalStepSize[i];
				}
			
			}



			ret = lbfgs(ndim, PrincipalStartPoint, &GradML, evaluate, progress, NULL, &param);
			currentrun++;

		}


		/* Report the result. */
		printf("L-BFGS optimization terminated with status code = %d, ML %g \n", ret, GradML);

		for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			if(((MNStruct *)GHSglobalcontext)->pulse->fitJump[k] == 0){

				((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k] = ((MNStruct *)GHSglobalcontext)->PreJumpVals[k];
			}
		}

		OutputMLFiles(ndim, GlobalMaxLikeVec, GlobalMaxLike, profdims+epochdims+epochpriordims);

		delete[] StartPoint;
		delete[] StepSize;

    		return ;

	}
	


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

	

	
	OutputMLFiles(ndim, GlobalMaxLikeVec, GlobalMaxLike, profdims+epochdims+epochpriordims);

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
	int start = 3*((MNStruct *)GHSglobalcontext)->TotalProfiles;
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

/*

 A subroutine for calculating n-d gaussian posteior
 The sampler requires -log(posterior) = v
 ndim -> number of dimensions
 x    -> input vector at which the posterior is be calculated
 g    -> gradient of the -log(posterior) = dPsi/dx
*/

void nd_uncorr_gauss_neg_log_post(int* ndim,double* x,double* v,double* g)
{
	int i;
	(*v)=0;
	for(i=0;i<*ndim;++i)
	{
		(*v)+=x[i]*x[i];
		g[i]=x[i];
	}
	(*v)*=0.5;
}




void GHSProfileDomainLike(int* ndim, double* PrinCube, double* likelihood, double* PrinGrad){

	int dotime = 0;
	int debug = ((MNStruct *)GHSglobalcontext)->debug; 
	int testdim = 6;
	double wtime = 0;

	int testuniformpriors = 0;


	int threads = omp_get_max_threads();

   	

	struct timeval tval_before, tval_after, tval_resultone, tval_resulttwo;

	if(dotime == 1){
		gettimeofday(&tval_before, NULL);
		wtime = omp_get_wtime ();
	}
	
	if(debug == 1)printf("In like \n");

	double *grad = new double[*ndim];
	for(int i = 0; i < *ndim; i++){
		//if(i!=testdim){Cube[i]=GlobalStartPoint[i];}
		//if(i==testdim){Cube[i]=100;}
		grad[i] = 0;
	}
	
	double *Cube = new double[*ndim];

	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Physical(PrinCube, Cube);
	}
	else{
		for(int i = 0; i < *ndim; i++){
			Cube[i] = PrinCube[i];
		}
	}



	if(dotime == 3){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to First Transform: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);
	}	


	//Cube[((MNStruct *)GHSglobalcontext)->TotalProfiles*3] = -6.45022;
	//Cube[((MNStruct *)GHSglobalcontext)->TotalProfiles*3+1] = -1.1342;

	//for(int i = ((MNStruct *)GHSglobalcontext)->TotalProfiles*3; i < *ndim; i++){
	//	printf("change: %i %g %g \n", i, PrinCube[i], Cube[i]);
	//}
	//sleep(5);



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

	for(int i = 0; i < profdims; i++){
//		Cube[i]=GlobalStartPoint[i];
	}


	for(int i =  profdims+ epochpriordims; i < *ndim; i++){
//		Cube[i]=GlobalStartPoint[i];
	}

	double **threadgrads = new double*[globalparams];
	for(int i = 0; i < globalparams; i++){ threadgrads[i] = new double[threads]();}


	double **PriorThreadGrads = new double*[epochpriordims];
	for(int i = 0; i < epochpriordims; i++){ PriorThreadGrads[i] = new double[threads]();}



	double uniformpriorterm = 0;
	double priorterm = 0;
	double priorDet = 0;
	double JacobianTerm = 0;


	//if(((MNStruct *)GHSglobalcontext)->numFitEQUAD == 1){Cube[profdims] = -7;}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	long double LDparams[((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps];
	int pcount = profdims+epochdims+epochpriordims;
	int phaseDim = profdims+epochdims+epochpriordims;
	int fitcount = 0;

	for(int p=0;p< ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; p++){
		//printf("Timing: %i %g %g \n", p, Cube[pcount], double(Cube[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[p][0])));
		LDparams[p]=Cube[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[p][0]);
		pcount++;

	}	
	
	long double phase = LDparams[0]*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;
	fitcount++;

	double PhasePrior = ((MNStruct *)GHSglobalcontext)->PhasePrior;
	double phasePriorTerm = 0.5*Cube[phaseDim]*Cube[phaseDim]/(PhasePrior*PhasePrior);
	double phasePriorGrad = Cube[phaseDim]/(PhasePrior*PhasePrior);


	if(((MNStruct *)GHSglobalcontext)->doLinear==0){

		for(int p=1;p<((MNStruct *)GHSglobalcontext)->numFitTiming;p++){
			((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][0]].val[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][1]] = LDparams[fitcount];
			fitcount++;
		}
		for(int p=0;p<((MNStruct *)GHSglobalcontext)->numFitJumps;p++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[((MNStruct *)GHSglobalcontext)->TempoJumpNums[p]]= LDparams[fitcount];
			fitcount++;
		}

		long double *JumpVec = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
		   JumpVec[p] = 0;
		   for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			for (int l=0;l<((MNStruct *)GHSglobalcontext)->pulse->obsn[p].obsNjump;l++){
				if(((MNStruct *)GHSglobalcontext)->pulse->obsn[p].jump[l]==k){
					if(((MNStruct *)GHSglobalcontext)->pulse->fitJump[k] == 0){

						JumpVec[p] += (long double)((MNStruct *)GHSglobalcontext)->PreJumpVals[k]/SECDAY;

					}
					else{

						JumpVec[p] += (long double) ((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]/SECDAY;
					}
				}
			 }
		   }


		}



		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
			((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat = ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase + JumpVec[p];
			//printf("Mod SAT: %.20Lg %.20Lg %.20Lg\n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase,  JumpVec[p], ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat);
		}
		for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]= 0;
		}


	
		fastformBatsAll(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars);       
		formResiduals(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars,0);      

		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
	//		printf("JumpVec: %i %g %g %g %g \n", p, (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat, (double)JumpVec[p], (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[p].bat, (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[p].residual);
		}
		delete[] JumpVec;
		//UpdatePhysDVector(GHSglobalcontext, GlobalDMatrix, ((MNStruct *)GHSglobalcontext)->pulse->nobs);

	}

	
	if(debug == 1)printf("Phase: %g \n", (double)phase);
	if(debug == 1)printf("Formed Residuals \n");



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Epoch Signals///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	int EPriorcount = profdims;
	double *EQUAD;
	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){ EQUAD = new double[((MNStruct *)GHSglobalcontext)->numFitEQUAD]; }

	//double *EQUADSignal = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();
	
//	double EQUADGrad = 0;
//	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD == 1){
	for(int eq = 0; eq < ((MNStruct *)GHSglobalcontext)->numFitEQUAD; eq++){

		if(debug == 1){printf("EQUAD %i : %i %g \n",eq, EPriorcount,Cube[EPriorcount]);}
		if(testuniformpriors == 0){
			if(Cube[EPriorcount] < -10){ 
				priorterm -= 2*log(10.0)*(Cube[EPriorcount]+10); 
//				EQUADGrad -= log(10.0);

				for(int i = 0; i < threads; i++){
					PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
				}	
			}
		}
		else{
			priorterm -= 2*log(10.0)*Cube[EPriorcount]; 
//			EQUADGrad -= log(10.0);

			for(int i = 0; i < threads; i++){
				PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
			}
		}

		EQUAD[eq]=pow(10.0,Cube[EPriorcount]);
		EPriorcount++;

		if(((MNStruct *)GHSglobalcontext)->EQUADPriorType ==1) {uniformpriorterm += log(EQUAD[eq]);}


		
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////DM Variations///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	double *SignalVec;
	double *PowerSpec;
	double *FMatrix;

	double FreqLike = 0;
	int startpos = 0;
	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;

		FMatrix = new double[FitDMCoeff*((MNStruct *)GHSglobalcontext)->numProfileEpochs]();
		PowerSpec = new double[FitDMCoeff]();
		SignalVec = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();

		
		double *SignalCoeff = new double[FitDMCoeff];
		for(int i = 0; i < FitDMCoeff; i++){
			SignalCoeff[startpos + i] = Cube[pcount];
			pcount++;
		}
			
		

		double Tspan = ((MNStruct *)GHSglobalcontext)->Tspan;
		double f1yr = 1.0/3.16e7;


		if(((MNStruct *)GHSglobalcontext)->incDM==5){

			double DMamp=Cube[EPriorcount];



			
			DMamp=pow(10.0, DMamp); 
			if(((MNStruct *)GHSglobalcontext)->DMPriorType == 1) { 
				priorterm -= 2*log(DMamp); 
				for(int i = 0; i < threads; i++){
					PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
				}
			}

			EPriorcount++;
			double DMindex=Cube[EPriorcount];
			EPriorcount++;

			for (int i=0; i< FitDMCoeff/2; i++){

				FreqLike += SignalCoeff[i]*SignalCoeff[i] + SignalCoeff[i+FitDMCoeff/2]*SignalCoeff[i+FitDMCoeff/2];
				
				double freq = ((double)(i+1.0))/Tspan;
				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freq*365.25,(-DMindex))/(Tspan*24*60*60);
		
				PowerSpec[startpos + i] = sqrt(rho);
				PowerSpec[startpos + i + FitDMCoeff/2] = sqrt(rho);

				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  
				
			}
		}


		if(((MNStruct *)GHSglobalcontext)->incDM==6){




			for (int i=0; i< FitDMCoeff/2; i++){
			

				double DMAmp = pow(10.0, Cube[EPriorcount]);	
				double freq = ((double)(i+1.0))/Tspan;
				
				if(((MNStruct *)GHSglobalcontext)->DMPriorType == 1) { 
					priorterm -= 2*log(DMAmp); 
					for(int i = 0; i < threads; i++){
						PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
					}
				}

				EPriorcount++;

				FreqLike += SignalCoeff[i]*SignalCoeff[i] + SignalCoeff[i+FitDMCoeff/2]*SignalCoeff[i+FitDMCoeff/2];

				double rho = (DMAmp*DMAmp);
				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  




			}
		}

		
		for(int i=0;i< FitDMCoeff/2;i++){
			int DMt = 0;
			for(int k=0;k<((MNStruct *)GHSglobalcontext)->numProfileEpochs;k++){
				double time=(double)((MNStruct *)GHSglobalcontext)->pulse->obsn[DMt].bat;
	
				FMatrix[k + (i+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]=cos(2*M_PI*(double(i+1)/Tspan)*time);
				FMatrix[k + (i+FitDMCoeff/2+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs] = sin(2*M_PI*(double(i+1)/Tspan)*time);
				DMt += ((MNStruct *)GHSglobalcontext)->numChanPerInt[k];

			}
		}

		
		vector_dgemv(FMatrix,SignalCoeff,SignalVec,((MNStruct *)GHSglobalcontext)->numProfileEpochs, FitDMCoeff,'N');
		startpos=FitDMCoeff;
		delete[] SignalCoeff;	
		

    	}


	

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	if(dotime == 3){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to ProfileParams: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);
	}



	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = 0;
	for(int ep = 0; ep < NEpochs; ep++){ TotalProfs +=  ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep]; }



	 long double **ProfileBats=new long double*[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

	        int ProfNbin  = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[i][1];

		ProfileBats[i] = new long double[2];

		ProfileBats[i][0] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
		ProfileBats[i][1] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][ProfNbin-1][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;

	      
	      ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - phase - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;
	//	printf("MBat %i %.15Lg \n", i, ModelBats[i]);

		//printf("res: %i %g \n", i, (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual);
	 }


	if(((MNStruct *)GHSglobalcontext)->doLinear>0){

		//printf("in do linear %i %i\n", phaseDim+1, phaseDim+((MNStruct *)GHSglobalcontext)->numFitTiming);

		double *TimingParams = new double[((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps]();
		double *TimingSignal = new double[((MNStruct *)GHSglobalcontext)->pulse->nobs];

		int fitcount = 1;
		for(int i = phaseDim+1; i < phaseDim+((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps; i++){
			//printf("in do linear %i %i \n", i, fitcount);
			///printf("in do linear %i %i %g \n", i, fitcount, Cube[i]);
			TimingParams[fitcount] = Cube[i];
			fitcount++;
		}

		vector_dgemv(GlobalDMatrixVec, TimingParams,TimingSignal,((MNStruct *)GHSglobalcontext)->pulse->nobs,((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps,'N');

		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){
//			printf("LinTime: %i %g %g \n", i, (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[i].sat, TimingSignal[i] );
	     		 ModelBats[i] -=  (long double) TimingSignal[i]/SECDAY;

		}	

		delete[] TimingSignal;
		delete[] TimingParams;	

	}
	

	if(dotime == 3){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to sorting bats: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);
	}


	int maxshapecoeff = 0;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
		if(debug == 1){printf("num coeff in comp %i: %i \n", i, numcoeff[i]);}
	}

	
        int *numProfileStocCoeff = ((MNStruct *)GHSglobalcontext)->numshapestoccoeff;
        int totalshapestoccoeff = ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff;


	int *numEvoCoeff = ((MNStruct *)GHSglobalcontext)->numEvoCoeff;
	int totalEvoCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoCoeff;

	int *numEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->numEvoFitCoeff;
	int totalEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoFitCoeff;

	int *numProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->numProfileFitCoeff;
	int totalProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff;



	int totalCoeffForMult = 0;
	int *NumCoeffForMult = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		NumCoeffForMult[i] = numProfileStocCoeff[i];
		if(numEvoCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numEvoCoeff[i];}
		if(numProfileFitCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numProfileFitCoeff[i];}
		totalCoeffForMult += NumCoeffForMult[i];
		if(debug == 1){printf("num coeff for mult from comp %i: %i \n", i, NumCoeffForMult[i]);}
	}

        int *numShapeToSave = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                numShapeToSave[i] = numProfileStocCoeff[i];
                if(numEvoCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numEvoCoeff[i];}
                if(numProfileFitCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numProfileFitCoeff[i];}
                if(debug == 1){printf("saved %i %i \n", i, numShapeToSave[i]);}
        }
        int totShapeToSave = 0;
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                totShapeToSave += numShapeToSave[i];
        }

	double shapecoeff[totshapecoeff];
	double StocProfCoeffs[totalshapestoccoeff];
	double **EvoCoeffs=new double*[((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){EvoCoeffs[i] = new double[totalEvoCoeff];}

	double ProfileFitCoeffs[totalProfileFitCoeff];
	

	for(int i =0; i < totshapecoeff; i++){
		shapecoeff[i]=((MNStruct *)GHSglobalcontext)->MeanProfileShape[i];
		//printf("loaded shape coeff %i %g \n", i, shapecoeff[i]);
	}
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		for(int i =0; i < totalEvoCoeff; i++){
			EvoCoeffs[p][i]=((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p][i];
			//printf("loaded evo coeff %i %g \n", i, EvoCoeffs[i]);
		}
	}
	if(debug == 1){printf("Filled %i Coeff, %i EvoCoeff \n", totshapecoeff, totalEvoCoeff);}
	double *betas = new double[((MNStruct *)GHSglobalcontext)->numProfComponents]();
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		betas[i] = ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[i]*((MNStruct *)GHSglobalcontext)->ReferencePeriod;
	}


	//double *StocProfPriorGrad = new double[totalshapestoccoeff]();
	for(int i =0; i < totalshapestoccoeff; i++){
		//Cube[EPriorcount]=GlobalStartPoint[EPriorcount];
		StocProfCoeffs[i]= pow(10.0, Cube[EPriorcount]);

		if(debug == 1){printf("Stoc: %i %g \n",EPriorcount,Cube[EPriorcount]);}



		if(testuniformpriors == 0){
			if(Cube[EPriorcount] < -5){ 
				priorterm -= 2*log(10.0)*(Cube[EPriorcount]+5); 

				for(int i = 0; i < threads; i++){
					PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;
				}

			}
		}
		else{

			priorterm -= 2*log(10.0)*Cube[EPriorcount]; 

			for(int i = 0; i < threads; i++){
				PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
			}
		}

		
		EPriorcount++;
	}

	for(int i =0; i < totalProfileFitCoeff; i++){
		
		ProfileFitCoeffs[i]= Cube[pcount];
		pcount++;
	}


	int  EvoFreqExponent = 1;
	if(((MNStruct *)GHSglobalcontext)->FitEvoExponent == 1){
		EvoFreqExponent =  floor(Cube[pcount]);
		pcount++;
	}
	
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		int cpos = 0;
		for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
			for(int i =0; i < numEvoFitCoeff[c]; i++){
				//printf("Evo: %i %g %g \n", i, EvoCoeffs[p][i+cpos] ,Cube[pcount]);
				EvoCoeffs[p][i+cpos] += Cube[pcount];
				pcount++;
			}
			cpos += numEvoCoeff[c];
		
		}
	}

	double LinearProfileWidth=0;
	if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1){ 
		LinearProfileWidth = Cube[pcount];
		pcount++;
	}
        double *LinearWidthEvoTime;
        if(((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
		LinearWidthEvoTime = new double[((MNStruct *)GHSglobalcontext)->incWidthEvoTime];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->incWidthEvoTime; i++){
			LinearWidthEvoTime[i] = Cube[pcount];
			pcount++;
		}
        }
	double EvoProfileWidth=0;
	if(((MNStruct *)GHSglobalcontext)->incProfileEvo == 2){
		EvoProfileWidth = Cube[pcount];
		pcount++;
	}

	double EvoEnergyProfileWidth=0;
	if(((MNStruct *)GHSglobalcontext)->incProfileEnergyEvo == 2){
		EvoEnergyProfileWidth = Cube[pcount];
		pcount++;
	}

	if(totshapecoeff+1>=totalshapestoccoeff+1){
		maxshapecoeff=totshapecoeff+1;
	}
	if(totalshapestoccoeff+1 > totshapecoeff+1){
		maxshapecoeff=totalshapestoccoeff+1;
	}


	double modelflux=0;
	int fluxpos = 0;
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		for(int j =0; j < numcoeff[i]; j=j+2){
			modelflux+=sqrt(sqrt(M_PI))*sqrt(betas[i])*pow(2.0, 0.5*(1.0-j))*sqrt(((MNStruct *)GHSglobalcontext)->Binomial[j])*shapecoeff[fluxpos+j];
		}
		fluxpos+= numcoeff[i];
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;





	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;


	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to Start of main loop: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

		//wtime = omp_get_wtime ( ) - wtime;
		//cout << "  Elapsed wall clock time = " << wtime << "\n";

		wtime = omp_get_wtime ( ); 
	}	

	
	double *eplikes = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();


	#pragma omp parallel for 
	for(int ep = minep; ep < maxep; ep++){

		double *ProfileModCoeffs = new double[totalCoeffForMult];

		int thisthread = omp_get_thread_num();
	
		int t = 0;

		for(int sep = 0; sep < ep; sep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[sep];
		}
		


		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];
		int NEpochBins = NChanInEpoch*GlobalNBins;


		double EpochChisq = 0;	
		double EpochdetN = 0;
		double EpochLike = 0;

		double EQUADSignal = 0;

		int EDimcount = profdims + epochpriordims + ep*DimsPerEpoch;
		int OriginalEDimCount = EDimcount;
		if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){


			int EQIndex = 0;
                        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[t];}


			if(debug == 1){printf("EQUAD Signal: %i %g \n",EDimcount,Cube[EDimcount]);}


			EQUADSignal = Cube[EDimcount]*EQUAD[EQIndex];
			eplikes[ep] += -0.5*Cube[EDimcount]*Cube[EDimcount];

			EDimcount++;
		}

		double *OneEpochStochAmps = new double[totalshapestoccoeff];
		for(int i = 0; i < totalshapestoccoeff; i++){

			if(debug == 1){printf("Stoc Signal: %i %g \n",EDimcount,Cube[EDimcount]);}

			OneEpochStochAmps[i] = Cube[EDimcount]*StocProfCoeffs[i];
			eplikes[ep] += -0.5*Cube[EDimcount]*Cube[EDimcount];
                        EDimcount++;	
		}
		
		



		for(int ch = 0; ch < NChanInEpoch; ch++){

			EDimcount = OriginalEDimCount;

			if(dotime == 2){	
				gettimeofday(&tval_before, NULL);
			}
			if(debug == 1){
				printf("In toa %i \n", t);
				printf("sat: %.15Lg \n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].sat);
			}


			int nTOA = t;


			double ProfileBaseline = Cube[nTOA*perProfDims + 0];
			double ProfileAmp = Cube[nTOA*perProfDims + 1];
			double ProfileSigma = Cube[nTOA*perProfDims + 2];

			double detN  = 0;
			double Chisq  = 0;
  

			double profilelike=0;

			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];
			long double FoldingPeriodDays = FoldingPeriod/SECDAY;

			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;

			double Tobs = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][2];
			double noiseval = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][3];
			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;


			double *shapevec  = new double[Nbins];
			double *ProfileModVec = new double[Nbins]();
			double *ProfileJitterModVec = new double[Nbins]();
			double *ProfileWidthModVec = new double[Nbins]();


			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get and modify binpos///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////
		    
			long double binpos = ModelBats[nTOA]; 


			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
				binpos += (long double) EQUADSignal/SECDAY;
			}

			if(((MNStruct *)GHSglobalcontext)->incDM > 4 || ((MNStruct *)GHSglobalcontext)->yearlyDM == 1){
	                	double DMKappa = 2.410*pow(10.0,-16);
        		        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				long double DMshift = (long double)(SignalVec[ep]*DMScale);
				
		 		binpos+=DMshift/SECDAY;

			}


			if(binpos < ProfileBats[nTOA][0])binpos+=FoldingPeriodDays;
			if(binpos > ProfileBats[nTOA][1])binpos-=FoldingPeriodDays;


			if(binpos > ProfileBats[nTOA][1]){printf("OverBoard! %.10Lg %.10Lg %.10Lg\n", binpos, ProfileBats[nTOA][1], (binpos-ProfileBats[nTOA][1])/FoldingPeriodDays);}

			long double minpos = binpos - FoldingPeriodDays/2;
			if(minpos < ProfileBats[nTOA][0])minpos=ProfileBats[nTOA][0];
			long double maxpos = binpos + FoldingPeriodDays/2;
			if(maxpos> ProfileBats[nTOA][1])maxpos =ProfileBats[nTOA][1];

			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get Interpolation Bin///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////

			int InterpBin = 0;
			double FirstInterpTimeBin = 0;
			int  NumWholeBinInterpOffset = 0;

			if(((MNStruct *)GHSglobalcontext)->InterpolateProfile == 1){

		
				long double timediff = 0;
				long double bintime = ProfileBats[t][0];


				if(bintime  >= minpos && bintime <= maxpos){
				    timediff = bintime - binpos;
				}
				else if(bintime < minpos){
				    timediff = FoldingPeriodDays+bintime - binpos;
				}
				else if(bintime > maxpos){
				    timediff = bintime - FoldingPeriodDays - binpos;
				}

				timediff=timediff*SECDAY;

				double OneBin = FoldingPeriod/Nbins;
				int NumBinsInTimeDiff = floor(timediff/OneBin + 0.5);
				double WholeBinsInTimeDiff = NumBinsInTimeDiff*FoldingPeriod/Nbins;
				double OneBinTimeDiff = -1*((double)timediff - WholeBinsInTimeDiff);

				double PWrappedTimeDiff = (OneBinTimeDiff - floor(OneBinTimeDiff/OneBin)*OneBin);

				if(debug == 1)printf("Making InterpBin: %g %g %i %g %g %g\n", (double)timediff, OneBin, NumBinsInTimeDiff, WholeBinsInTimeDiff, OneBinTimeDiff, PWrappedTimeDiff);

				InterpBin = floor(PWrappedTimeDiff/((MNStruct *)GHSglobalcontext)->InterpolatedTime+0.5);
				if(InterpBin >= ((MNStruct *)GHSglobalcontext)->NumToInterpolate)InterpBin -= ((MNStruct *)GHSglobalcontext)->NumToInterpolate;

				FirstInterpTimeBin = -1*(InterpBin-1)*((MNStruct *)GHSglobalcontext)->InterpolatedTime;

				if(debug == 1)printf("Interp Time Diffs: %g %g %g %g \n", ((MNStruct *)GHSglobalcontext)->InterpolatedTime, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime, PWrappedTimeDiff, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime-PWrappedTimeDiff);

				double FirstBinOffset = timediff-FirstInterpTimeBin;
				double dNumWholeBinOffset = FirstBinOffset/(FoldingPeriod/Nbins);
				int  NumWholeBinOffset = 0;

				NumWholeBinInterpOffset = floor(dNumWholeBinOffset+0.5);
	
				if(debug == 1)printf("Interp bin %i is: %i , First Bin is %g, Offset is %i \n", t, InterpBin, FirstInterpTimeBin, NumWholeBinInterpOffset);


			}



			//if(nTOA==1)printf("Jump Check: %i %i %g %.20Lg %g %g\n", (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].residual, ModelBats[nTOA], InterpBin, NumWholeBinInterpOffset, Cube[profdims+1], GlobalDMatrix[nTOA][1]*Cube[profdims+1]);


	
			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed up to start of Interp: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freq - reffreq)/1000.0;
			double freqscale = pow(freqdiff, EvoFreqExponent);


			double snr = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].snr;
			double tobs = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].tobs; 


			snr = snr*3600/tobs;

			double refSN = 1000;
			double SNdiff =  snr/refSN;
			double SNscale = snr-refSN;




			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(int i =0; i < totalCoeffForMult; i++){
				ProfileModCoeffs[i]=0;	
			}				

			int cpos = 0;
			int epos = 0;
			int fpos = 0;
			int spos = 0;
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				
				for(int i =0; i < numProfileFitCoeff[c]; i++){
					ProfileModCoeffs[i+cpos] += ProfileFitCoeffs[i+fpos];

				}
				for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i =0; i < numEvoCoeff[c]; i++){
						ProfileModCoeffs[i+cpos] += EvoCoeffs[p][i+epos]*pow(freqscale, p+1);						
					}
				}

				for(int i = 0; i < numProfileStocCoeff[c]; i++){
					ProfileModCoeffs[i+cpos] += OneEpochStochAmps[i+spos]*modelflux;
					//printf("adding stoc: %i %i %i %g \n", ep, ch, i,OneEpochStochAmps[i+spos] );
				}

				cpos += NumCoeffForMult[c];
				epos += numEvoCoeff[c];
				fpos += numProfileFitCoeff[c];	
				spos += numProfileStocCoeff[c];
			}


			

			if(totalCoeffForMult > 0){
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ProfileModCoeffs,ProfileModVec,Nbins,totalCoeffForMult,'N');
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfileVec[InterpBin], ProfileModCoeffs,ProfileJitterModVec,Nbins,totalCoeffForMult,'N');         
				if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1 || ((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
                                        vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfileVec[InterpBin], ProfileModCoeffs,ProfileWidthModVec,Nbins,totalCoeffForMult,'N');
                                }


			}


	

			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Modded Profile: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}





			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Fill Arrays with interpolated state//////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			
			int ZeroWrap = GHSWrap(0 + NumWholeBinInterpOffset, 0, Nbins-1);


			//double corrAmpGrad = 0;
	
			for(int j =0; j < Nbins-ZeroWrap; j++){

				double NewIndex = (j + NumWholeBinInterpOffset);
				int Nj =  ZeroWrap+j;

				

				double widthTerm = (ProfileWidthModVec[Nj] + ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj])*(LinearProfileWidth); 

				for(int ew = 0; ew < ((MNStruct *)GHSglobalcontext)->incWidthEvoTime; ew++){
					widthTerm += (ProfileWidthModVec[Nj] + ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj])*LinearWidthEvoTime[ew]*pow( (((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].bat - ((MNStruct *)GHSglobalcontext)->ProfEvoTimeRef)/365.25, ew+1);
				}



				double evoWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoProfileWidth*freqscale;
				double SNWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoEnergyProfileWidth*SNscale;

				shapevec[j] = ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj] + widthTerm + ProfileModVec[Nj] + evoWidthTerm + SNWidthTerm;
//				if(t==0 || t == 10 || t == 100 || t == 300){printf("Shape: %i %i %g\n", t, j, shapevec[j]);}

			}
	
			for(int j = Nbins-ZeroWrap; j < Nbins; j++){

				double NewIndex = (j + NumWholeBinInterpOffset);
				int Nj =  j-Nbins+ZeroWrap;


				double widthTerm = (ProfileWidthModVec[Nj] + ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj])*(LinearProfileWidth); 

				for(int ew = 0; ew < ((MNStruct *)GHSglobalcontext)->incWidthEvoTime; ew++){
					widthTerm += (ProfileWidthModVec[Nj] + ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj])*LinearWidthEvoTime[ew]*pow( (((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].bat - ((MNStruct *)GHSglobalcontext)->ProfEvoTimeRef)/365.25, ew+1);
				}

 
				double evoWidthTerm = 0;// ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoProfileWidth*freqscale;
				double SNWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoEnergyProfileWidth*SNscale;

				shapevec[j] = ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj] + widthTerm + ProfileModVec[Nj] + evoWidthTerm + SNWidthTerm;
//				if(t==0 || t == 10 || t == 100 || t == 300){printf("Shape: %i %i %g\n", t, j, shapevec[j]);}

			}


			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Filled Arrays: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}


			


		///////////////////////////////////////////Marginalise over arbitrary offset and absolute amplitude////////////////////////////////////////////////////////////

			      
			Chisq = 0;


			double noise = ProfileSigma*ProfileSigma;
			double OffPulsedetN = log(noise);
			detN = ProfNbins*OffPulsedetN;
			noise=1.0/noise;

			double *NResVec = new double[Nbins]();
		

			ZeroWrap = GHSWrap(0 + (Nbins - NumWholeBinInterpOffset), 0, Nbins-1);

			int ZeroOffset = (Nbins-GHSWrap((Nbins - NumWholeBinInterpOffset), 0, Nbins-1))%BinRatio;


			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Allocated Arrays: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}
			

			int startpoint = 0;
			int stoppoint = Nbins - ZeroWrap;

			//for(int i = 0; i < ProfNbins; i++){
			//	if(ProfNbins < Nbins){printf("Prof bins: %i %i %.10g\n", t, i, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1]);}
			//}

			for(int i = startpoint+ZeroOffset; i < stoppoint; i+=BinRatio){

				
			//	printf("Comp: %i %i %g %g \n", i, ZeroWrap+i, shapevec[ZeroWrap+i],((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][i]);	
				int Nj = ZeroWrap+i;


				double gres = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1] - ProfileAmp*shapevec[Nj] - ProfileBaseline;
				NResVec[i] = gres*noise;
				Chisq +=  gres*NResVec[i];


				grad[nTOA*perProfDims + 0] +=   -NResVec[i]; //Baseline
				grad[nTOA*perProfDims + 1] +=   -shapevec[Nj]*NResVec[i]; //Amp
				grad[nTOA*perProfDims + 2] +=   -gres*NResVec[i]/ProfileSigma + 1.0/ProfileSigma;    //Sigma

				//if(ProfNbins  == 256)printf("Make NRVec %i %i %i %i %i %i %.10g %.10g %.10g \n", startpoint+ZeroOffset, ZeroOffset, stoppoint, t, i, Nj, double(Nj)/BinRatio, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1], ProfileAmp*shapevec[Nj] + ProfileBaseline);


			}

			startpoint = Nbins-ZeroWrap; // + BinRatio - (Nbins-ZeroWrap-1)%BinRatio - 1; 
			stoppoint = Nbins;
	
			for(int i = startpoint; i < stoppoint; i+=BinRatio){


				
			
				int Nj = ZeroWrap+i-Nbins;

				double gres = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1] - ProfileAmp*shapevec[Nj] - ProfileBaseline;
	
				NResVec[i] = gres*noise;
				Chisq +=  gres*NResVec[i];


				grad[nTOA*perProfDims + 0] +=   -NResVec[i]; //Baseline
				grad[nTOA*perProfDims + 1] +=   -shapevec[Nj]*NResVec[i]; //Amp
				grad[nTOA*perProfDims + 2] +=   -gres*NResVec[i]/ProfileSigma + 1.0/ProfileSigma;    //Sigma

				//if(ProfNbins == 256)printf("Make NRVec %i %i %i %i %i %i %.10g %.10g %.10g \n", startpoint, ZeroOffset, stoppoint, t, i, Nj, double(Nj)/BinRatio, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1], ProfileAmp*shapevec[Nj] + ProfileBaseline);

			}



			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Get EpochGrads///////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int B2C = 0;
			double PhaseGrad = 0;
			for(int j = 0; j < Nbins; j++){
				ProfileJitterModVec[j] += ((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][j];
//				PhaseGrad += -1*NResVec[j]*(((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][j] + ProfileJitterModVec[j])*ProfileAmp;
				PhaseGrad += -1*NResVec[j]*ProfileJitterModVec[j]*ProfileAmp;
				
			}

			for(int j = 0; j < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; j++){	
				double phasefac = 1;
				if(j==0){phasefac *= double(((MNStruct *)GHSglobalcontext)->ReferencePeriod);}		
				threadgrads[B2C][thisthread] += GlobalDMatrix[nTOA][j]*phasefac*PhaseGrad;
				B2C++;
//				if(j==1 &&  GlobalDMatrix[nTOA][j] != 0 )printf("JG: %i %i %g %g %g \n", nTOA, j, PhaseGrad,  GlobalDMatrix[nTOA][j], threadgrads[B2C-1][thisthread]);
			}


			int PriorGradLoopCount = 0;
			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

				int EQIndex = 0;
	                        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[nTOA];}


				grad[EDimcount] += 1*EQUADSignal/EQUAD[EQIndex]/NChanInEpoch - PhaseGrad*EQUAD[EQIndex];
                                PriorThreadGrads[PriorGradLoopCount+EQIndex][thisthread] += (-PhaseGrad*EQUADSignal/EQUAD[EQIndex])*log(10.0)*EQUAD[EQIndex];

				EDimcount++;
				PriorGradLoopCount+=((MNStruct *)GHSglobalcontext)->numFitEQUAD;
			}
			

			if(((MNStruct *)GHSglobalcontext)->incDM==5){

				int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;



				double DMKappa = 2.410*pow(10.0,-16);
        		        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				double DMshift = (SignalVec[ep]*DMScale);

				
				//Gradient for Amplitude of DM
				PriorThreadGrads[PriorGradLoopCount][thisthread] += -PhaseGrad*DMshift*log(10.0);
				PriorGradLoopCount++;

				
				int DMAmpDim = profdims + epochdims + epochpriordims + ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;

				double DMSpecGradTerm = 0;

				for (int i=0; i< FitDMCoeff/2; i++){

					double KPF1 = DMScale*FMatrix[ep + i*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[i];
					double KPF2 = DMScale*FMatrix[ep + (i+FitDMCoeff/2)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[i+FitDMCoeff/2];

					threadgrads[B2C+i][thisthread] +=  Cube[DMAmpDim+i]/TotalProfs - PhaseGrad*KPF1;
					threadgrads[B2C+i+FitDMCoeff/2][thisthread] +=  Cube[DMAmpDim+i+FitDMCoeff/2]/TotalProfs - PhaseGrad*KPF2;

					double freq = ((double)(i+1.0))/((MNStruct *)GHSglobalcontext)->Tspan;

					DMSpecGradTerm += KPF1*(-0.5*log(freq*365.25))*Cube[DMAmpDim+i];
					DMSpecGradTerm += KPF2*(-0.5*log(freq*365.25))*Cube[DMAmpDim+i+FitDMCoeff/2];

				
				}

				//Gradient for Spectral Index of DM
				PriorThreadGrads[PriorGradLoopCount][thisthread] += -PhaseGrad*DMSpecGradTerm;


				B2C += FitDMCoeff;
				PriorGradLoopCount++;
			}


			if(totalCoeffForMult > 0){

				double *GradResVec = new double[totalCoeffForMult];

				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], NResVec,GradResVec,Nbins,totalCoeffForMult,'T');

				for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){

					double fval = -1*ProfileAmp*pow(freqscale, p);

					int cpos = 0;
					int fpos = 0;
					for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
						for(int i = 0; i < numProfileFitCoeff[c]; i++){
			
							threadgrads[B2C + p*totalProfileFitCoeff + fpos + i][thisthread] += fval*GradResVec[i+cpos];
						}
						fpos += numProfileFitCoeff[c];	
						cpos += NumCoeffForMult[c];
					}
				}
				B2C += (((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1)*totalProfileFitCoeff;


				int cpos = 0;
				int spos = 0;
				for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
					for(int i = 0; i < numProfileStocCoeff[c]; i++){
						grad[EDimcount] += OneEpochStochAmps[i + spos]/StocProfCoeffs[i + spos]/NChanInEpoch - modelflux*ProfileAmp*GradResVec[i+cpos]*StocProfCoeffs[i + spos];
						EDimcount++;

						//StocProfPriorGrad[i+spos] += -ProfileAmp*OneEpochStochAmps[i + spos]*GradResVec[i+cpos]*log(10.0);
						PriorThreadGrads[PriorGradLoopCount][thisthread] +=  -modelflux*ProfileAmp*OneEpochStochAmps[i + spos]*GradResVec[i+cpos]*log(10.0);
						PriorGradLoopCount++;
					}
					spos += numProfileStocCoeff[c];
					cpos += NumCoeffForMult[c];
				}

				delete[] GradResVec;
			}


			if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1){
				double LinearPGrad = 0;
				for(int i = 0; i < Nbins; i++){	
					LinearPGrad += -1*ProfileAmp*NResVec[i]*(((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][i]+ProfileWidthModVec[i]);
				}
				threadgrads[B2C][thisthread] +=	LinearPGrad;
				B2C++;
			}

			if(((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
				for(int ew = 0; ew < ((MNStruct *)GHSglobalcontext)->incWidthEvoTime; ew++){
					double LinearPGrad = 0;
					for(int i = 0; i < Nbins; i++){	
						LinearPGrad += -1*ProfileAmp*NResVec[i]*(((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][i]+ProfileWidthModVec[i])*pow((((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].bat - ((MNStruct *)GHSglobalcontext)->ProfEvoTimeRef)/365.25, ew+1);
					}
					threadgrads[B2C][thisthread] +=	LinearPGrad;
					B2C++;
				}
			}




			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Up to LA: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			//int id;

			  //  id = omp_get_thread_num ( );
			  //  cout << "  This is process " << id << " " << t << " " << detN << " " << Chisq << "\n";
			  			

			profilelike = -0.5*(detN + Chisq);
			//printf("Like: %i %.15g %.15g \n", nTOA, detN, Chisq);
			if(debug == 1)printf("Like: %i %.15g %.15g %.15g \n", nTOA, lnew, detN, Chisq);

			EpochChisq+=Chisq;
			EpochdetN+=detN;
			EpochLike+=profilelike;

			delete[] shapevec;
			delete[] ProfileModVec;
			delete[] ProfileJitterModVec;
			delete[] ProfileWidthModVec;
			
			delete[] NResVec;

			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  End of Epoch: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			t++;
		}



///////////////////////////////////////////


		if(((MNStruct *)GHSglobalcontext)->incWideBandNoise == 0){
			eplikes[ep] += EpochLike;
			lnew += EpochLike;
		}
		delete[] ProfileModCoeffs;
		delete[] OneEpochStochAmps;

////////////////////////////////////////////
	
	}




	double omplnew = 0;
	for(int ep = minep; ep < maxep; ep++){
		omplnew += eplikes[ep];
	}

	for(int i = 0; i < epochpriordims; i++){
		double gsum = 0;
                for(int j = 0; j < threads; j++){
                        gsum += PriorThreadGrads[i][j];
                }
		 grad[profdims + i] = gsum;
	}

	for(int i = 0; i < globalparams; i++){
		double gsum = 0;
		for(int j = 0; j < threads; j++){
			gsum += threadgrads[i][j];
		}
		 grad[profdims + epochdims + epochpriordims + i] = gsum;
	}
	grad[profdims + epochdims + epochpriordims]+=phasePriorGrad;
	///printf("phase prior grad %g \n", phasePriorGrad);
	if(dotime == 1){
		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed,  End of Loop: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

		wtime = omp_get_wtime ( ) - wtime;
		cout << "  Elapsed wall clock time = " << wtime << "\n";

		wtime = omp_get_wtime ( ) ;
	}

	
	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->pulse->nobs; j++){
	    delete[] ProfileBats[j];
	}
	delete[] ProfileBats;
	delete[] ModelBats;
	delete[] numcoeff;
	delete[] NumCoeffForMult;
	delete[] numShapeToSave;

	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] EvoCoeffs[j];
	}
	delete[] EvoCoeffs;
	delete[] eplikes;


	for(int i = 0; i < globalparams; i++){
		delete[] threadgrads[i];
	}
	delete[] threadgrads;


	for(int i = 0; i < epochpriordims; i++){
		delete[] PriorThreadGrads[i];
	}
	delete[] PriorThreadGrads;
	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){ delete[] EQUAD; }

	if(((MNStruct *)GHSglobalcontext)->incWidthEvoTime){
		delete[] LinearWidthEvoTime;
	}
	

	double finallikelihood = -1*omplnew + 0.5*FreqLike + phasePriorTerm + 0.5*priorterm;
	*likelihood=finallikelihood;
        //printf("in like %g %g \n", finallikelihood, liketest[0]);	
	if(-1*finallikelihood > GlobalMaxLike && ((MNStruct *)GHSglobalcontext)->WriteNewML == 1){
		GlobalMaxLike = -1*finallikelihood;
		printf("Like: %.16g \n", -1*finallikelihood);
		write_newMLLike(*ndim, Cube, GlobalMaxLike);
		if(debug == 1)printf("Written new ML file\n");
	
	}
	//sleep(5);

	if(debug == 0){
		for(int i =profdims; i < profdims+1; i++){
		//	printf("Parameters and Grads: %i %g %g %.10g \n", i, Cube[i], grad[i], omplnew-phasePriorTerm-0.5*priorterm);
		}
	}

	for(int i = 0; i < profdims; i++){
	//	grad[i]=0;
	}


	for(int i =  profdims + epochpriordims ; i < *ndim; i++){
	//grad[i]=0;
	}

	if(debug == 1)printf("End Like: %.10g \n", omplnew-phasePriorTerm-0.5*priorterm );
	
	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Principal(PrinGrad, grad);
	}
	else{
		for(int i = 0; i < *ndim; i++){
			PrinGrad[i] = grad[i];
		}
	}
	



	if(dotime == 1){
		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed,  Up to End: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

		wtime = omp_get_wtime ( ) - wtime;
		cout << "  Elapsed wall clock time = " << wtime << "\n";


	}




	delete[] Cube;
	delete[] grad;
	delete[] betas;

	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		delete[] FMatrix;
		delete[] PowerSpec;
		delete[] SignalVec;
	}
	
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
		//StepSize[i] = 0;
		MaxAmps[i] = 0;
	}


	struct timeval tval_before, tval_after, tval_resultone, tval_resulttwo;

	if(dotime == 1){
		gettimeofday(&tval_before, NULL);
	}
	
	if(debug == 1)printf("In like \n");

	double *EpochPriorWeights = new double[epochpriordims]();



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	
	long double LDparams[((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps];	
	int pcount=0;
	int fitcount = 0;

	for(int p=0;p< ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; p++){
		//printf("Timing: %i %g %g \n", p, Cube[pcount], double(Cube[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[p][0])));
		LDparams[p] = ((MNStruct *)GHSglobalcontext)->LDpriors[p][0];
		pcount++;

	}	
	
	long double phase = LDparams[0]*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;
	fitcount++;

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form SATS///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	if(((MNStruct *)GHSglobalcontext)->doLinear==0){

		for(int p=1;p<((MNStruct *)GHSglobalcontext)->numFitTiming;p++){
			((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][0]].val[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][1]] = LDparams[fitcount];
			fitcount++;
		}
		for(int p=0;p<((MNStruct *)GHSglobalcontext)->numFitJumps;p++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[((MNStruct *)GHSglobalcontext)->TempoJumpNums[p]]= LDparams[fitcount];
			fitcount++;
		}

		long double *JumpVec = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
		   JumpVec[p] = 0;
		   for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			for (int l=0;l<((MNStruct *)GHSglobalcontext)->pulse->obsn[p].obsNjump;l++){
				if(((MNStruct *)GHSglobalcontext)->pulse->obsn[p].jump[l]==k){
					if(((MNStruct *)GHSglobalcontext)->pulse->fitJump[k] == 0){

						JumpVec[p] += (long double)((MNStruct *)GHSglobalcontext)->PreJumpVals[k]/SECDAY;

					}
					else{

						JumpVec[p] += (long double) ((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]/SECDAY;
					}
				}
			 }
		   }

		}



		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
			((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat = ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase + JumpVec[p];
			//printf("Mod SAT: %.20Lg %.20Lg %.20Lg\n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase,  JumpVec[p], ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat);
		}
		for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]= 0;
		}

		delete[] JumpVec;

	
		fastformBatsAll(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars,0);       /* Form residuals */


		//UpdatePhysDVector(GHSglobalcontext, GlobalDMatrix, ((MNStruct *)GHSglobalcontext)->pulse->nobs);

	}



	
	if(debug == 1)printf("Phase: %g \n", (double)phase);
	if(debug == 1)printf("Formed Residuals \n");




/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////DM Variations///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	double *SignalVec;
	double *PowerSpec;
	double *FMatrix;
	double *SignalCoeff;

	double FreqLike = 0;
	int startpos = 0;
	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;

		FMatrix = new double[FitDMCoeff*((MNStruct *)GHSglobalcontext)->numProfileEpochs]();
		PowerSpec = new double[FitDMCoeff]();
		SignalVec = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();

		
		SignalCoeff = new double[FitDMCoeff];
		SignalCoeff[startpos+0] =  0.676197;
		SignalCoeff[startpos+1] =  -0.38914;
		SignalCoeff[startpos+2] =  -1.02693;
		SignalCoeff[startpos+3] =  0.0234188;
		SignalCoeff[startpos+4] =  0.0522272;
		SignalCoeff[startpos+5] =  -0.987655;
		SignalCoeff[startpos+6] =  -0.67628;
		SignalCoeff[startpos+7] =  -0.561906;
		SignalCoeff[startpos+8] =  -0.0522284;
		SignalCoeff[startpos+9] =  0.194622;
		SignalCoeff[startpos+10] =  0.0613292;
		SignalCoeff[startpos+11] =  -0.651747;
		SignalCoeff[startpos+12] =  0.587244;
		SignalCoeff[startpos+13] =  0.142931;
		SignalCoeff[startpos+14] =  -0.136876;
		SignalCoeff[startpos+15] =  0.584652;
		SignalCoeff[startpos+16] =  0.747795;
		SignalCoeff[startpos+17] =  2.37543;
		SignalCoeff[startpos+18] =  -0.521506;
		SignalCoeff[startpos+19] =  0.724316;
		SignalCoeff[startpos+20] =  -3.05028;
		SignalCoeff[startpos+21] =  0.278179;
		SignalCoeff[startpos+22] =  -1.07729;
		SignalCoeff[startpos+23] =  -1.04139;
		SignalCoeff[startpos+24] =  0.157578;
		SignalCoeff[startpos+25] =  -0.0178301;
		SignalCoeff[startpos+26] =  -0.0457155;
		SignalCoeff[startpos+27] =  0.249242;
		SignalCoeff[startpos+28] =  0.512135;
		SignalCoeff[startpos+29] =  -0.276627;



		double Tspan = ((MNStruct *)GHSglobalcontext)->Tspan;
		double f1yr = 1.0/3.16e7;


		if(((MNStruct *)GHSglobalcontext)->incDM==5){

			double DMamp = ((MNStruct *)GHSglobalcontext)->pulse->TNDMAmp; 



			
			DMamp=pow(10.0, DMamp); 
			double DMindex = ((MNStruct *)GHSglobalcontext)->pulse->TNDMGam;

			for (int i=0; i< FitDMCoeff/2; i++){

				FreqLike += SignalCoeff[i]*SignalCoeff[i] + SignalCoeff[i+FitDMCoeff/2]*SignalCoeff[i+FitDMCoeff/2];
				
				double freq = ((double)(i+1.0))/Tspan;
				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freq*365.25,(-DMindex))/(Tspan*24*60*60);
		
				PowerSpec[startpos + i] = sqrt(rho);
				PowerSpec[startpos + i + FitDMCoeff/2] = sqrt(rho);

				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  
				
			}
		}


		if(((MNStruct *)GHSglobalcontext)->incDM==6){




			for (int i=0; i< FitDMCoeff/2; i++){
			

//	//			double DMAmp = pow(10.0, Cube[EPriorcount]);	
//				double freq = ((double)(i+1.0))/Tspan;
				


//				FreqLike += SignalCoeff[i]*SignalCoeff[i] + SignalCoeff[i+FitDMCoeff/2]*SignalCoeff[i+FitDMCoeff/2];

//				double rho = (DMAmp*DMAmp);
//				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
//				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  




			}
		}

		
		for(int i=0;i< FitDMCoeff/2;i++){
			int DMt = 0;
			for(int k=0;k<((MNStruct *)GHSglobalcontext)->numProfileEpochs;k++){
				double time=(double)((MNStruct *)GHSglobalcontext)->pulse->obsn[DMt].bat;
	
				FMatrix[k + (i+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]=cos(2*M_PI*(double(i+1)/Tspan)*time);
				FMatrix[k + (i+FitDMCoeff/2+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs] = sin(2*M_PI*(double(i+1)/Tspan)*time);
				DMt += ((MNStruct *)GHSglobalcontext)->numChanPerInt[k];

			}
		}

		
		vector_dgemv(FMatrix,SignalCoeff,SignalVec,((MNStruct *)GHSglobalcontext)->numProfileEpochs, FitDMCoeff,'N');


		if(((MNStruct *)GHSglobalcontext)->incDM==5){

			double DMamp = ((MNStruct *)GHSglobalcontext)->pulse->TNDMAmp; 



			
			DMamp=pow(10.0, DMamp); 
			double DMindex = ((MNStruct *)GHSglobalcontext)->pulse->TNDMGam;

			for (int i=0; i< FitDMCoeff/2; i++){

				
				double freq = ((double)(i+1.0))/Tspan;
				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freq*365.25,(-DMindex))/(Tspan*24*60*60);
		

				SignalCoeff[i] = SignalCoeff[i]/sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]/sqrt(rho);  
				
			}
		}
		startpos=FitDMCoeff;
		

    	}


	

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int TotalProfs = ((MNStruct *)GHSglobalcontext)->TotalProfiles;




	 long double **ProfileBats=new long double*[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

	        int ProfNbin  = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[i][1];

		ProfileBats[i] = new long double[2];

		ProfileBats[i][0] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
		ProfileBats[i][1] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][ProfNbin-1][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;

	      
	      ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - phase - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;

	 }


	int maxshapecoeff = 0;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
		if(debug == 1){printf("num coeff in comp %i: %i \n", i, numcoeff[i]);}
	}

	
        int *numProfileStocCoeff = ((MNStruct *)GHSglobalcontext)->numshapestoccoeff;
        int totalshapestoccoeff = ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff;


	int *numEvoCoeff = ((MNStruct *)GHSglobalcontext)->numEvoCoeff;
	int totalEvoCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoCoeff;

	int *numEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->numEvoFitCoeff;
	int totalEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoFitCoeff;

	int *numProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->numProfileFitCoeff;
	int totalProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff;



	int totalCoeffForMult = 0;
	int *NumCoeffForMult = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		NumCoeffForMult[i] = 0;
		if(numEvoCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numEvoCoeff[i];}
		if(numProfileFitCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numProfileFitCoeff[i];}
		totalCoeffForMult += NumCoeffForMult[i];
		if(debug == 1){printf("num coeff for mult from comp %i: %i \n", i, NumCoeffForMult[i]);}
	}

	((MNStruct *)GHSglobalcontext)->totalCoeffForMult = totalCoeffForMult;

        int *numShapeToSave = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                numShapeToSave[i] = numProfileStocCoeff[i];
                if(numEvoCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numEvoCoeff[i];}
                if(numProfileFitCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numProfileFitCoeff[i];}
                if(debug == 1){printf("saved %i %i \n", i, numShapeToSave[i]);}
        }
        int totShapeToSave = 0;
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                totShapeToSave += numShapeToSave[i];
        }

	double shapecoeff[totshapecoeff];
	double StocProfCoeffs[totalshapestoccoeff];
	double **EvoCoeffs=new double*[((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){EvoCoeffs[i] = new double[totalEvoCoeff];}

	double ProfileFitCoeffs[totalProfileFitCoeff];
	double ProfileModCoeffs[totalCoeffForMult];
	
	for(int i = 0; i < totalProfileFitCoeff; i++){ProfileFitCoeffs[i] = 0;}
	for(int i = 0; i < totalCoeffForMult; i++){ProfileModCoeffs[i] = 0;}

	for(int i =0; i < totshapecoeff; i++){
		shapecoeff[i]=((MNStruct *)GHSglobalcontext)->MeanProfileShape[i];
	}
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		for(int i =0; i < totalEvoCoeff; i++){
			EvoCoeffs[p][i]=((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p][i];
			//printf("loaded evo coeff %i %g \n", i, EvoCoeffs[i]);
		}
	}
	if(debug == 1){printf("Filled %i Coeff, %i EvoCoeff \n", totshapecoeff, totalEvoCoeff);}
	double *betas = new double[((MNStruct *)GHSglobalcontext)->numProfComponents]();
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		betas[i] = ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[i]*((MNStruct *)GHSglobalcontext)->ReferencePeriod;
	}

	int cpos = 0;
	int epos = 0;
	int fpos = 0;

	int  EvoFreqExponent = 1;	
	if(totshapecoeff+1>=totalshapestoccoeff+1){
		maxshapecoeff=totshapecoeff+1;
	}
	if(totalshapestoccoeff+1 > totshapecoeff+1){
		maxshapecoeff=totalshapestoccoeff+1;
	}




	double modelflux=0;
	cpos = 0;
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		for(int j =0; j < numcoeff[i]; j=j+2){
			modelflux+=sqrt(sqrt(M_PI))*sqrt(betas[i])*pow(2.0, 0.5*(1.0-j))*sqrt(((MNStruct *)GHSglobalcontext)->Binomial[j])*shapecoeff[cpos+j];
		}
		cpos+= numcoeff[i];
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;


	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to Start of main loop: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		
	}	


	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int t = 0;
	if(((MNStruct *)GHSglobalcontext)->SubIntToFit != -1){
		minep = ((MNStruct *)GHSglobalcontext)->SubIntToFit;
		maxep = minep+1;
		for(int ep = 0; ep < minep; ep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];
		}
	}


	

	double EQStart = 0;
	double *EpochShapeStart = new double[totalshapestoccoeff]();
	for(int ep = minep; ep < maxep; ep++){
	



		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];
		int NEpochBins = NChanInEpoch*GlobalNBins;


		double EpochChisq = 0;	
		double EpochdetN = 0;
		double EpochLike = 0;

		double EpochPhaseHess = 0;
		double EpochPhaseGrad = 0;

		double *EpochShapeGrad = new double[totalshapestoccoeff]();
		double *EpochShapeHess = new double[totalshapestoccoeff]();

		int EpochMsize = 0;

		EpochMsize += totalshapestoccoeff;

		if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
			EpochMsize++;
		}
		if(((MNStruct *)GHSglobalcontext)->incDMEQUAD > 0){
			EpochMsize++;
		}
		if(((MNStruct *)GHSglobalcontext)->incWidthJitter > 0){
			EpochMsize++;
		}

		double *EpochMNMHess = new double[EpochMsize*EpochMsize]();
		double *EpochMMatrix = new double[EpochMsize*GlobalNBins]();
		double *EpochdNM = new double[EpochMsize]();
		double *EpochTempdNM = new double[EpochMsize]();
		

		for(int ch = 0; ch < NChanInEpoch; ch++){

			if(dotime == 1){	
				gettimeofday(&tval_before, NULL);
			}
			if(debug == 1){
				printf("In toa %i \n", t);
				printf("sat: %.15Lg \n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].sat);
			}


			int nTOA = t;




			double detN  = 0;
			double Chisq  = 0;
  

			double profilelike=0;

			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];
			long double FoldingPeriodDays = FoldingPeriod/SECDAY;

			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;

			double Tobs = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][2];
			double noiseval = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][3];
			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;

			double *shapevec  = new double[Nbins];
			double *ProfileModVec = new double[Nbins]();
		   	double *ProfileJitterModVec = new double[Nbins](); 
			long double binpos = ModelBats[nTOA];


 
			if(((MNStruct *)GHSglobalcontext)->incDM > 4 || ((MNStruct *)GHSglobalcontext)->yearlyDM == 1){
	                	double DMKappa = 2.410*pow(10.0,-16);
        		        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				long double DMshift = (long double)(SignalVec[ep]*DMScale);
				
		 		binpos+=DMshift/SECDAY;

			}

			if(binpos < ProfileBats[nTOA][0])binpos+=FoldingPeriodDays;
			if(binpos > ProfileBats[nTOA][1])binpos-=FoldingPeriodDays;


			if(binpos > ProfileBats[nTOA][1]){printf("OverBoard! %.10Lg %.10Lg %.10Lg\n", binpos, ProfileBats[nTOA][1], (binpos-ProfileBats[nTOA][1])/FoldingPeriodDays);}

			long double minpos = binpos - FoldingPeriodDays/2;
			if(minpos < ProfileBats[nTOA][0])minpos=ProfileBats[nTOA][0];
			long double maxpos = binpos + FoldingPeriodDays/2;
			if(maxpos> ProfileBats[nTOA][1])maxpos =ProfileBats[nTOA][1];

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Interpolation Bin///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

			int InterpBin = 0;
			double FirstInterpTimeBin = 0;
			int  NumWholeBinInterpOffset = 0;

			if(((MNStruct *)GHSglobalcontext)->InterpolateProfile == 1){

		
				long double timediff = 0;
				long double bintime = ProfileBats[t][0];


				if(bintime  >= minpos && bintime <= maxpos){
				    timediff = bintime - binpos;
				}
				else if(bintime < minpos){
				    timediff = FoldingPeriodDays+bintime - binpos;
				}
				else if(bintime > maxpos){
				    timediff = bintime - FoldingPeriodDays - binpos;
				}

				timediff=timediff*SECDAY;

				double OneBin = FoldingPeriod/Nbins;
				int NumBinsInTimeDiff = floor(timediff/OneBin + 0.5);
				double WholeBinsInTimeDiff = NumBinsInTimeDiff*FoldingPeriod/Nbins;
				double OneBinTimeDiff = -1*((double)timediff - WholeBinsInTimeDiff);

				double PWrappedTimeDiff = (OneBinTimeDiff - floor(OneBinTimeDiff/OneBin)*OneBin);

				if(debug == 1)printf("Making InterpBin: %g %g %i %g %g %g\n", (double)timediff, OneBin, NumBinsInTimeDiff, WholeBinsInTimeDiff, OneBinTimeDiff, PWrappedTimeDiff);

				InterpBin = floor(PWrappedTimeDiff/((MNStruct *)GHSglobalcontext)->InterpolatedTime+0.5);
				if(InterpBin >= ((MNStruct *)GHSglobalcontext)->NumToInterpolate)InterpBin -= ((MNStruct *)GHSglobalcontext)->NumToInterpolate;

				FirstInterpTimeBin = -1*(InterpBin-1)*((MNStruct *)GHSglobalcontext)->InterpolatedTime;

				if(debug == 1)printf("Interp Time Diffs: %g %g %g %g \n", ((MNStruct *)GHSglobalcontext)->InterpolatedTime, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime, PWrappedTimeDiff, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime-PWrappedTimeDiff);

				double FirstBinOffset = timediff-FirstInterpTimeBin;
				double dNumWholeBinOffset = FirstBinOffset/(FoldingPeriod/Nbins);
				int  NumWholeBinOffset = 0;

				NumWholeBinInterpOffset = floor(dNumWholeBinOffset+0.5);
	
				if(debug == 1)printf("Interp bin is: %i , First Bin is %g, Offset is %i \n", InterpBin, FirstInterpTimeBin, NumWholeBinInterpOffset);


			}

	
			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed up to start of Interp: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}









			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freq - reffreq)/1000.0;
			double freqscale = pow(freqdiff, EvoFreqExponent);


			double snr = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].snr;
			double tobs = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].tobs; 


			snr = snr*3600/tobs;

			double refSN = 1000;
			double SNdiff =  snr/refSN;
			double SNscale = snr-refSN;




			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(int i =0; i < totalCoeffForMult; i++){
				ProfileModCoeffs[i]=0;	
			}				

			cpos = 0;
			epos = 0;
			fpos = 0;

	
			
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				
				for(int i =0; i < numProfileFitCoeff[c]; i++){
					//printf("PFD: %i %g \n",ProfileFitCoeffs[i+fpos]);
					ProfileModCoeffs[i+cpos] += ProfileFitCoeffs[i+fpos];

				}
				for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i =0; i < numEvoCoeff[c]; i++){
						ProfileModCoeffs[i+cpos] += EvoCoeffs[p][i+epos]*pow(freqscale, p+1);						
					}
				}
				cpos += NumCoeffForMult[c];
				epos += numEvoCoeff[c];
				fpos += numProfileFitCoeff[c];	
			}


			
			double *OneProfileParamHessian = new double[totalCoeffForMult*totalCoeffForMult];
			if(totalCoeffForMult > 0){
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ProfileModCoeffs,ProfileModVec,Nbins,totalCoeffForMult,'N');
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfileVec[InterpBin], ProfileModCoeffs,ProfileJitterModVec,Nbins,totalCoeffForMult,'N');
				vector_dgemm(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin] , OneProfileParamHessian, Nbins, totalCoeffForMult,Nbins, totalCoeffForMult, 'T', 'N');

			}


			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Modded Profile: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Fill Arrays with interpolated state//////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


				
			for(int j =0; j < Nbins; j++){

				double NewIndex = (j + NumWholeBinInterpOffset);
				int Nj =  GHSWrap(j + NumWholeBinInterpOffset, 0, Nbins-1);//(int)(NewIndex - floor(NewIndex/Nbins)*Nbins);


				//printf("Shape: %i %i %g %g \n", t, j, ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj], ProfileModVec[Nj] );
				shapevec[j] = ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj] + ProfileModVec[Nj];
			}


			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Filled Arrays: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////Get ML Amps and Baselines and Noise/////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			double maxamp = 0;
			for(int i =0; i < Nbins; i++){
				if(shapevec[i] > maxamp){maxamp = shapevec[i];}
			}

			int AmpZeroWrap = Nbins-GHSWrap((Nbins - NumWholeBinInterpOffset), 0, Nbins-1);
			int AmpZeroOffset = 0;//AmpZeroWrap%BinRatio;



			double *datavec = new double[ProfNbins]();
			double *ModelVec = new double[2*ProfNbins]();
			double *MM = new double[2*2]();
			double *Md = new double[2]();
			double *TempMd = new double[2]();

			for(int i =0; i < ProfNbins; i++){
				datavec[i] = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1];
				ModelVec[i + 0*ProfNbins] = 1;
				ModelVec[i + 1*ProfNbins] = shapevec[AmpZeroOffset+i*BinRatio];
			}

			vector_dgemm(ModelVec, ModelVec, MM, ProfNbins, 2,ProfNbins, 2, 'T', 'N');
			vector_dgemv(ModelVec, datavec,Md,ProfNbins,2,'T');

			for(int i =0; i < 2; i++){
				TempMd[i] = Md[i];
			}

			int info=0;
			double Margindet = 0;
			vector_dpotrfInfo(MM, 2, Margindet, info);
			vector_dpotrsInfo(MM, TempMd, 2, info);

			double MLSigma = 0;
			for(int i =0; i < ProfNbins; i++){
				//double res = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1] - shapevec[AmpZeroOffset+i*BinRatio]*MLAmp;
				//MLSigma += (res - noisemean)*(res - noisemean);

				double res2 = datavec[i] - TempMd[0] - TempMd[1]*shapevec[AmpZeroOffset+i*BinRatio];
				MLSigma += res2*res2;
				//if(t==15657){printf("MaxAmp 15657: %i %g %g \n", i, datavec[i], TempMd[0]+ TempMd[1]*shapevec[AmpZeroOffset+i*BinRatio]);}
			}

			MLSigma = sqrt(MLSigma/ProfNbins);
			double noisemean = TempMd[0];
			double MLAmp = TempMd[1];

			MaxAmps[nTOA*3 + 0] = noisemean;
			MaxAmps[nTOA*3 + 1] = MLAmp;
			MaxAmps[nTOA*3 + 2] = MLSigma;

			//printf("Max: %i %g %g %g %g %g %g\n", t, MaxAmps[nTOA*3 + 0], MaxAmps[nTOA*3 + 1], MaxAmps[nTOA*3 + 2]);


			delete[] datavec;
			delete[] ModelVec;
			delete[] MM;
			delete[] Md;
			delete[] TempMd;

			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Up to LA: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}









/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Phase Term ////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		int HessCount = TotalProfs;


		double PhaseScale = MLAmp/MLSigma;
		double OnePhaseHess = 0;
		double OnePhaseGrad = 0;

		int ZeroWrap = Nbins-GHSWrap((Nbins - NumWholeBinInterpOffset), 0, Nbins-1);
		int ZeroOffset = ZeroWrap%BinRatio;

		double *NGRes = new double[Nbins]();
		for(int i =ZeroOffset; i < Nbins; i+=BinRatio){
			ProfileJitterModVec[i] += ((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i];
//			OnePhaseHess += PhaseScale*PhaseScale*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i]*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i];
			OnePhaseHess += PhaseScale*PhaseScale*ProfileJitterModVec[i]*ProfileJitterModVec[i];



			int RevWrap = GHSWrap(i + (Nbins - NumWholeBinInterpOffset), 0, Nbins-1);

			//printf("PhaseHess: %i %i %i %i %i %g \n", ZeroWrap, ZeroOffset, t, i, RevWrap, double(RevWrap)/BinRatio);

			double gres = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][RevWrap/BinRatio][1] - MLAmp*shapevec[RevWrap] - noisemean;
			NGRes[i] = gres/MLSigma/MLSigma;
			OnePhaseGrad += ProfileJitterModVec[i]*PhaseScale*gres/MLSigma;

		}

		EpochPhaseHess += OnePhaseHess;
		EpochPhaseGrad += OnePhaseGrad;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////Make EpochMatrix/////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int EpochDimCounter = 0;
			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
				

				for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
					EpochMMatrix[i + EpochDimCounter*Nbins] = ProfileJitterModVec[i]*MLAmp;
				}
				EpochDimCounter++;
			}


			int matrixcpos = 0;
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){

				for(int i = 0; i < numProfileStocCoeff[c]; i++){
					for(int j = ZeroOffset; j < Nbins; j+=BinRatio){
						EpochMMatrix[j + EpochDimCounter*Nbins] = ((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin][j+(i+matrixcpos)*Nbins]*MLAmp*modelflux;
					}
					EpochDimCounter++;

				}

				matrixcpos += NumCoeffForMult[c];
			}

			if(EpochDimCounter > 0){
				double *TempEpochMNMHess = new double[EpochMsize*EpochMsize]();
				vector_dgemm(EpochMMatrix, EpochMMatrix , TempEpochMNMHess, Nbins, EpochMsize,Nbins, EpochMsize, 'T', 'N');
				for(int i = 0; i < EpochMsize; i++){
					for(int j =0; j < EpochMsize; j++){
						EpochMNMHess[i+j*EpochMsize] += TempEpochMNMHess[i+j*EpochMsize]/MLSigma/MLSigma;
					}
				}

				delete[] TempEpochMNMHess;


				vector_dgemv(EpochMMatrix, NGRes, EpochTempdNM,Nbins,EpochMsize,'T');

				for(int i = 0; i < EpochMsize; i++){
					EpochdNM[i] += EpochTempdNM[i];
					//printf("Epoch dNM %i %i %g \n", t, i, EpochdNM[i]);
				}

			}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Per Epoch Parameters///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(epochpriordims > 0){ 



			HessCount += 1;
		}




		if(epochdims > 0){


			if(totalshapestoccoeff > 0) {
				double *BNRes = new double[totalCoeffForMult];
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], NGRes, BNRes,Nbins,totalCoeffForMult,'T');

				int cpos = 0;
				int spos = 0;
				for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
					for(int i = 0; i < numProfileStocCoeff[c]; i++){
	
						EpochShapeGrad[i+spos] += BNRes[i+cpos]*MLAmp*modelflux;
						EpochShapeHess[i+spos] += OneProfileParamHessian[i+cpos + (i+cpos)*totalCoeffForMult]*MLAmp*MLAmp*modelflux*modelflux/(MLSigma*MLSigma);
					}
					spos += numProfileStocCoeff[c];
					cpos += NumCoeffForMult[c];
				}

				delete[] BNRes;
			}

		} 

		delete[] NGRes;
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Remaining Global Parameters////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		int B2S = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + totalCoeffForMult*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly);
		int B2C = 0;


/*
		double *PhaseShapeCross = new double[totalCoeffForMult];

		vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin],PhaseShapeCross,Nbins,totalCoeffForMult,'T');

		for(int p = 0; p < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){
			for(int i =0; i < totalCoeffForMult; i++){
			
				PhaseScale = (pow(freqscale, p)*MLAmp)*(((MNStruct *)GHSglobalcontext)->ReferencePeriod*MLAmp)/(MLSigma*MLSigma);
				double val = PhaseShapeCross[i]*PhaseScale;

			}
		}

		delete[] PhaseShapeCross;
*/
		//printf("Phase Hessian: %g %g\n", OnePhaseHess, SSM[TotalProfs][0]);

		B2C++;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// One Profile Shape parameters///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
			for(int p = 0; p < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){
				for(int q = 0; q < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; q++){
					for(int i =0; i < totalCoeffForMult; i++){
						for(int j =0; j < totalCoeffForMult; j++){
							//if(p >= q){printf("val: %i %i %i %g %g \n", t, p, q, freqscale, pow(freqscale, p)*pow(freqscale, q));}
							double val = pow(freqscale, p)*pow(freqscale, q)*OneProfileParamHessian[j+i*totalCoeffForMult]*MLAmp*MLAmp/(MLSigma*MLSigma);
						//printf("Shape part: %i %i %g %g \n", p*totalCoeffForMult + B2C + i, q*totalCoeffForMult + B2C + j, SSM[TotalProfs][(p*totalCoeffForMult + B2C + i) + (q*totalCoeffForMult + B2C + j)*B2S], val);	
					//	SSM[HessCount][(p*totalCoeffForMult + B2C + i) + (q*totalCoeffForMult + B2C + j)*B2S] += val;
							
						}
					}
				}
			}
*/
			profilelike = -0.5*(detN + Chisq);

			if(debug == 1)printf("Like: %i %.15g %.15g %.15g \n", nTOA, lnew, detN, Chisq);

			EpochChisq+=Chisq;
			EpochdetN+=detN;
			EpochLike+=profilelike;

			delete[] shapevec;
			delete[] ProfileModVec;
			delete[] ProfileJitterModVec;
			delete[] OneProfileParamHessian;

			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  End of Epoch: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
			}

			t++;
		}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Finish Per Epoch Parameters////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		


		if(epochdims > 0){ 


			double *TempMNM = new double[EpochMsize*EpochMsize]();

			for(int i =0; i < EpochMsize; i++){
				EpochTempdNM[i] = EpochdNM[i];
			}

			int info=0;
			double Margindet = 0;

			//printf("Hess1: %g %g \n", EpochMNMHess[0],1.0/(pow(10.0, 2*0.0095)));

			int epochdimcount = 0;
			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

				int EQIndex = 0;
	                        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[t-1];}

//				printf("Epoch EQ: %i %i %i %s %g \n", ((MNStruct *)GHSglobalcontext)->numFitEQUAD, ep, t-1,((MNStruct *)GHSglobalcontext)->pulse->TNEQFlagVal[EQIndex], ((MNStruct *)GHSglobalcontext)->pulse->TNEQVal[EQIndex]);
				double EQ = ((MNStruct *)GHSglobalcontext)->pulse->TNEQVal[EQIndex];
				EpochMNMHess[epochdimcount+epochdimcount*EpochMsize] += 1.0/(pow(10.0,2*EQ));
				epochdimcount++;
			}

			for(int i = 0; i < totalshapestoccoeff; i++){
				//printf("ML Shape Stoc: %i %.10g \n", i,  ((MNStruct *)GHSglobalcontext)->MeanProfileStoc[i]);
				EpochMNMHess[epochdimcount+epochdimcount*EpochMsize] += 1.0/(pow(10.0, 2*((MNStruct *)GHSglobalcontext)->MeanProfileStoc[i]));
				epochdimcount++;
			}

			for(int i =0; i < EpochMsize*EpochMsize; i++){
				TempMNM[i] = EpochMNMHess[i];
			}

			//printf("Hess2: %g\n", EpochMNMHess[0] );	

			vector_dpotrfInfo(TempMNM, EpochMsize, Margindet, info);
			vector_dpotri(TempMNM, EpochMsize);


			
			vector_dpotrfInfo(EpochMNMHess, EpochMsize, Margindet, info);
			vector_dpotrsInfo(EpochMNMHess, EpochTempdNM, EpochMsize, info);



			

			int PerEpochDimCount = 0;
			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

				double MLEpochPShift =  EpochTempdNM[PerEpochDimCount];//EpochPhaseGrad/EpochPhaseHess;
				double MLErr = sqrt(TempMNM[PerEpochDimCount+PerEpochDimCount*EpochMsize]);


				MaxAmps[profdims + epochpriordims + ep*perEpochDims+PerEpochDimCount] 	= -1*MLEpochPShift;
//				printf("ML Phase shift for epoch %i : %g %g \n", ep, MLEpochPShift, MLErr);



				double sqErr = sqrt(2.0)*MLErr;
				double sqweight = 1.0/sqErr/sqErr;

				EpochPriorWeights[PerEpochDimCount] += sqweight;

				EQStart += MLEpochPShift*MLEpochPShift*sqweight;

				PerEpochDimCount++;
				
			}
			for(int i = 0; i < totalshapestoccoeff; i++){
				double MLEpochShapeShift = EpochTempdNM[PerEpochDimCount];//EpochShapeGrad[i]/EpochShapeHess[i];
				double MLErr = sqrt(TempMNM[PerEpochDimCount+PerEpochDimCount*EpochMsize]);


				MaxAmps[profdims + epochpriordims + ep*perEpochDims+PerEpochDimCount] 	= MLEpochShapeShift;
				//printf("Per Epoch Max Shape change: %i %i %i %g %g \n", ep, i, profdims + epochpriordims + ep*perEpochDims+PerEpochDimCount, MLEpochShapeShift, MLErr);

	
				double sqErr = sqrt(2.0)*MLErr;
				double sqweight = 1.0/sqErr/sqErr;

				EpochPriorWeights[PerEpochDimCount] += sqweight;

				EpochShapeStart[i] += MLEpochShapeShift*MLEpochShapeShift*sqweight;

				PerEpochDimCount++;
			}

			delete[] TempMNM;

		} 


		



		delete[] EpochShapeGrad;
		delete[] EpochShapeHess;

		delete[] EpochMNMHess;
		delete[] EpochdNM;
		delete[] EpochTempdNM;
		delete[] EpochMMatrix;
		


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// End of Epoch///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
	}




	
	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->pulse->nobs; j++){
	    delete[] ProfileBats[j];
	}
	delete[] ProfileBats;
	delete[] ModelBats;
	delete[] numcoeff;
	delete[] NumCoeffForMult;
	delete[] numShapeToSave;

	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] EvoCoeffs[j];
	}
	delete[] EvoCoeffs;


	if(debug == 1)printf("End Like: %.10g \n", lnew);



	int dimCount = 0;
	dimCount += profdims;




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Get Epoch Prior EVDs///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	if(epochpriordims > 0){ 

		int PerEpochDimCount = 0;
//		if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
		for(int eq = 0; eq < ((MNStruct *)GHSglobalcontext)->numFitEQUAD; eq++){
			double Average = sqrt(EQStart/EpochPriorWeights[PerEpochDimCount]);

			MaxAmps[dimCount + PerEpochDimCount] = ((MNStruct *)GHSglobalcontext)->pulse->TNEQVal[eq]; 
			PerEpochDimCount++;
			
		}

		if(((MNStruct *)GHSglobalcontext)->incDM == 5){

			MaxAmps[dimCount + PerEpochDimCount] = ((MNStruct *)GHSglobalcontext)->pulse->TNDMAmp;
			PerEpochDimCount++;
			MaxAmps[dimCount + PerEpochDimCount] = ((MNStruct *)GHSglobalcontext)->pulse->TNDMGam;
			PerEpochDimCount++;
		}		
		

		for(int i = 0; i < totalshapestoccoeff; i++){

			double Average = sqrt(EpochShapeStart[i]/EpochPriorWeights[PerEpochDimCount]);
			//printf("ML Shape Stoc: %i %.10g \n", i,  ((MNStruct *)GHSglobalcontext)->MeanProfileStoc[i]);
			MaxAmps[dimCount + PerEpochDimCount] =  ((MNStruct *)GHSglobalcontext)->MeanProfileStoc[i];//log10(Average);
			PerEpochDimCount++;
		}

	
		dimCount += epochpriordims;
	}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Get Epoch EVDs/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	if(epochdims > 0){ 

		int epT = 0;
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfileEpochs; i++){


			int EDimPCount = 0;
			int EDimAmpCount = 0;
			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

				int EQIndex = 0;
				if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[epT];}

				
				MaxAmps[dimCount + i*perEpochDims + EDimAmpCount] = MaxAmps[dimCount + i*perEpochDims + EDimAmpCount]/pow(10.0, MaxAmps[profdims+EDimPCount+EQIndex]);
				printf("Signal %i %s %g %i %g %g %g \n", i, ((MNStruct *)GHSglobalcontext)->pulse->obsn[epT].fname, (double) ((MNStruct *)GHSglobalcontext)->pulse->obsn[epT].sat, EDimAmpCount,  MaxAmps[profdims+EDimPCount+EQIndex], MaxAmps[dimCount + i*perEpochDims + EDimAmpCount], MaxAmps[dimCount + i*perEpochDims + EDimAmpCount]/pow(10.0, MaxAmps[profdims+EDimPCount+EQIndex]));


				EDimPCount += ((MNStruct *)GHSglobalcontext)->numFitEQUAD;
				EDimAmpCount += 1;
				

			}


			for(int j = 0; j < totalshapestoccoeff; j++){

				MaxAmps[dimCount + i*perEpochDims + EDimAmpCount] = MaxAmps[dimCount + i*perEpochDims + EDimAmpCount]/pow(10.0, MaxAmps[profdims+EDimPCount]);

				EDimPCount++;
				EDimAmpCount++;

			}

/*
			for(int j = 0; j < perEpochDims; j++){
				printf("Signal %i %i %g %g %g \n", i, j,  MaxAmps[profdims+j], MaxAmps[dimCount + i*perEpochDims + j], MaxAmps[dimCount + i*perEpochDims + j]/pow(10.0, MaxAmps[profdims+j]));
				MaxAmps[dimCount + i*perEpochDims + j] = MaxAmps[dimCount + i*perEpochDims + j]/pow(10.0, MaxAmps[profdims+j]);
			}

*/
			epT += ((MNStruct *)GHSglobalcontext)->numChanPerInt[i]; 
		}
		
		dimCount += epochdims;
	}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Finish Shape parameters////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int vsize = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + totalCoeffForMult*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly);
	if(vsize > 0){

		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; i++){
			MaxAmps[dimCount + i] = 0;
		}
		dimCount +=  ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;
		for(int i = 0; i < 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff; i++){
			MaxAmps[dimCount + i] = SignalCoeff[i];
		}
		dimCount += 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;

		 for(int i = 0; i < totalProfileFitCoeff*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly); i++){
			MaxAmps[dimCount + i] = 0;
                }
		

		dimCount += totalProfileFitCoeff*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly);

		if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1){
			MaxAmps[dimCount] = 0;
			dimCount++;
		}
		if(((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
			for(int ew = 0; ew < ((MNStruct *)GHSglobalcontext)->incWidthEvoTime; ew++){
				MaxAmps[dimCount] = 0;
				dimCount++;
			}
		}

	}


	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		delete[] FMatrix;
		delete[] PowerSpec;
		delete[] SignalVec;
		delete[] SignalCoeff;
	}

}

void GetMaxSteps(int ndim, double *MaxAmps, double *StepSize, double **SSM){


	int crossScale = 0.9;//0.99;

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
		StepSize[i] = 0;
	}


	struct timeval tval_before, tval_after, tval_resultone, tval_resulttwo;

	if(dotime == 1){
		gettimeofday(&tval_before, NULL);
	}
	
	if(debug == 1)printf("In like \n");



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	long double LDparams[((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps];
	int pcount = profdims+epochdims+epochpriordims;
	int fitcount = 0;
	//long double phase = -0.000180659229620516726*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;
	

	for(int p=0;p< ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; p++){
		//printf("Timing: %i %g\n", p, MaxAmps[pcount+p]);
		LDparams[p]=MaxAmps[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[p][0]);
		pcount++;

	}	
	
	long double phase = LDparams[0]*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;
	fitcount++;
	
	
//	int pcount=0;
//	long double phase = ((MNStruct *)GHSglobalcontext)->LDpriors[0][0]*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY; 
	//printf("Phase: %g %g\n", (double)phase, (double)((MNStruct *)GHSglobalcontext)->LDpriors[0][0]);


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form SATS///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	if(((MNStruct *)GHSglobalcontext)->doLinear==0){

		for(int p=1;p<((MNStruct *)GHSglobalcontext)->numFitTiming;p++){
			((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][0]].val[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][1]] = LDparams[fitcount];
			fitcount++;
		}
		for(int p=0;p<((MNStruct *)GHSglobalcontext)->numFitJumps;p++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[((MNStruct *)GHSglobalcontext)->TempoJumpNums[p]]= LDparams[fitcount];
			fitcount++;
		}

		long double *JumpVec = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
		   JumpVec[p] = 0;
		   for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			for (int l=0;l<((MNStruct *)GHSglobalcontext)->pulse->obsn[p].obsNjump;l++){
				if(((MNStruct *)GHSglobalcontext)->pulse->obsn[p].jump[l]==k){
					if(((MNStruct *)GHSglobalcontext)->pulse->fitJump[k] == 0){

						JumpVec[p] += (long double)((MNStruct *)GHSglobalcontext)->PreJumpVals[k]/SECDAY;

					}
					else{

						JumpVec[p] += (long double) ((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]/SECDAY;
					}
				}
			 }
		   }

		}



		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
			((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat = ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase + JumpVec[p];
			//printf("Mod SAT: %.20Lg %.20Lg %.20Lg\n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase,  JumpVec[p], ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat);
		}
		for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]= 0;
		}

		delete[] JumpVec;

	
		fastformBatsAll(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars,0);       /* Form residuals */


		//UpdatePhysDVector(GHSglobalcontext, GlobalDMatrix, ((MNStruct *)GHSglobalcontext)->pulse->nobs);

	}

	
	if(debug == 1)printf("Phase: %g \n", (double)phase);
	if(debug == 1)printf("Formed Residuals \n");



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Epoch Signals///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	int EPriorcount = profdims;
	double *EQUAD;
	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){ EQUAD = new double[((MNStruct *)GHSglobalcontext)->numFitEQUAD]; }

	for(int eq = 0; eq < ((MNStruct *)GHSglobalcontext)->numFitEQUAD; eq++){


		if(debug == 1){printf("EQUAD: %i %g \n",EPriorcount,MaxAmps[EPriorcount]);}

		EQUAD[eq]=pow(10.0,MaxAmps[EPriorcount]);
		EPriorcount++;

		
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////DM Variations///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	double *SignalVec;
	double *PowerSpec;
	double *FMatrix;

	int startpos = 0;
	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;

		FMatrix = new double[FitDMCoeff*((MNStruct *)GHSglobalcontext)->numProfileEpochs]();
		PowerSpec = new double[FitDMCoeff]();
		SignalVec = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();

		
		double *SignalCoeff = new double[FitDMCoeff];
		for(int i = 0; i < FitDMCoeff; i++){
			SignalCoeff[startpos + i] = MaxAmps[pcount];
			pcount++;
		}
			
		

		double Tspan = ((MNStruct *)GHSglobalcontext)->Tspan;
		double f1yr = 1.0/3.16e7;


		if(((MNStruct *)GHSglobalcontext)->incDM==5){

			double DMamp=MaxAmps[EPriorcount];


			DMamp=pow(10.0, DMamp); 
			EPriorcount++;
			double DMindex=MaxAmps[EPriorcount];
			EPriorcount++;

			for (int i=0; i< FitDMCoeff/2; i++){

				double freq = ((double)(i+1.0))/Tspan;
				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freq*365.25,(-DMindex))/(Tspan*24*60*60);
		
				PowerSpec[startpos + i] = sqrt(rho);
				PowerSpec[startpos + i + FitDMCoeff/2] = sqrt(rho);

				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  
				
			}
		}


		if(((MNStruct *)GHSglobalcontext)->incDM==6){




			for (int i=0; i< FitDMCoeff/2; i++){
			

				double DMAmp = pow(10.0, MaxAmps[EPriorcount]);	
				double freq = ((double)(i+1.0))/Tspan;


				EPriorcount++;


				double rho = (DMAmp*DMAmp);
				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  




			}
		}

		
		for(int i=0;i< FitDMCoeff/2;i++){
			int DMt = 0;
			for(int k=0;k<((MNStruct *)GHSglobalcontext)->numProfileEpochs;k++){
				double time=(double)((MNStruct *)GHSglobalcontext)->pulse->obsn[DMt].bat;
	
				FMatrix[k + (i+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]=cos(2*M_PI*(double(i+1)/Tspan)*time);
				FMatrix[k + (i+FitDMCoeff/2+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs] = sin(2*M_PI*(double(i+1)/Tspan)*time);
				DMt += ((MNStruct *)GHSglobalcontext)->numChanPerInt[k];

			}
		}

		
		vector_dgemv(FMatrix,SignalCoeff,SignalVec,((MNStruct *)GHSglobalcontext)->numProfileEpochs, FitDMCoeff,'N');
		startpos=FitDMCoeff;
		delete[] SignalCoeff;	
		

    	}





/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	int TotalProfs = ((MNStruct *)GHSglobalcontext)->TotalProfiles;




	 long double **ProfileBats=new long double*[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

	        int ProfNbin  = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[i][1];

		ProfileBats[i] = new long double[2];

		ProfileBats[i][0] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
		ProfileBats[i][1] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][ProfNbin-1][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;

	      
	      ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - phase - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;
	 }


	int maxshapecoeff = 0;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
		if(debug == 1){printf("num coeff in comp %i: %i \n", i, numcoeff[i]);}
	}

	
        int *numProfileStocCoeff = ((MNStruct *)GHSglobalcontext)->numshapestoccoeff;
        int totalshapestoccoeff = ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff;


	int *numEvoCoeff = ((MNStruct *)GHSglobalcontext)->numEvoCoeff;
	int totalEvoCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoCoeff;

	int *numEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->numEvoFitCoeff;
	int totalEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoFitCoeff;

	int *numProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->numProfileFitCoeff;
	int totalProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff;



	int totalCoeffForMult = 0;
	int *NumCoeffForMult = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		NumCoeffForMult[i] = 0;
		if(numEvoCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numEvoCoeff[i];}
		if(numProfileFitCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numProfileFitCoeff[i];}
		totalCoeffForMult += NumCoeffForMult[i];
		if(debug == 1){printf("num coeff for mult from comp %i: %i \n", i, NumCoeffForMult[i]);}
	}

	((MNStruct *)GHSglobalcontext)->totalCoeffForMult = totalCoeffForMult;

        int *numShapeToSave = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                numShapeToSave[i] = numProfileStocCoeff[i];
                if(numEvoCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numEvoCoeff[i];}
                if(numProfileFitCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numProfileFitCoeff[i];}
                if(debug == 1){printf("saved %i %i \n", i, numShapeToSave[i]);}
        }
        int totShapeToSave = 0;
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                totShapeToSave += numShapeToSave[i];
        }

	double shapecoeff[totshapecoeff];
	double StocProfCoeffs[totalshapestoccoeff];
	double **EvoCoeffs=new double*[((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){EvoCoeffs[i] = new double[totalEvoCoeff];}

	double ProfileFitCoeffs[totalProfileFitCoeff];
	double ProfileModCoeffs[totalCoeffForMult];
	
	for(int i = 0; i < totalProfileFitCoeff; i++){ProfileFitCoeffs[i] = 0;}
	for(int i = 0; i < totalCoeffForMult; i++){ProfileModCoeffs[i] = 0;}

	for(int i =0; i < totshapecoeff; i++){
		shapecoeff[i]=((MNStruct *)GHSglobalcontext)->MeanProfileShape[i];
	}
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		for(int i =0; i < totalEvoCoeff; i++){
			EvoCoeffs[p][i]=((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p][i];
			//printf("loaded evo coeff %i %g \n", i, EvoCoeffs[i]);
		}
	}
	if(debug == 1){printf("Filled %i Coeff, %i EvoCoeff \n", totshapecoeff, totalEvoCoeff);}
	double *betas = new double[((MNStruct *)GHSglobalcontext)->numProfComponents]();
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		betas[i] = ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[i]*((MNStruct *)GHSglobalcontext)->ReferencePeriod;
	}


	for(int i =0; i < totalshapestoccoeff; i++){
		//Cube[EPriorcount]=GlobalStartPoint[EPriorcount];
		StocProfCoeffs[i]= pow(10.0, MaxAmps[EPriorcount]);

		if(debug == 1){printf("Stoc: %i %g \n",EPriorcount,MaxAmps[EPriorcount]);}
		
		EPriorcount++;
	}




	int cpos = 0;
	int epos = 0;
	int fpos = 0;

	for(int i =0; i < totalProfileFitCoeff; i++){
		
		ProfileFitCoeffs[i]= MaxAmps[pcount];
		pcount++;
	}

	
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		int cpos = 0;
		for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
			for(int i =0; i < numEvoFitCoeff[c]; i++){
				EvoCoeffs[p][i+cpos] += MaxAmps[pcount];
				pcount++;
			}
			cpos += numEvoCoeff[c];
		
		}
	}


	int  EvoFreqExponent = 1;	
	if(totshapecoeff+1>=totalshapestoccoeff+1){
		maxshapecoeff=totshapecoeff+1;
	}
	if(totalshapestoccoeff+1 > totshapecoeff+1){
		maxshapecoeff=totalshapestoccoeff+1;
	}


	double modelflux=0;
	cpos = 0;
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		for(int j =0; j < numcoeff[i]; j=j+2){
			modelflux+=sqrt(sqrt(M_PI))*sqrt(betas[i])*pow(2.0, 0.5*(1.0-j))*sqrt(((MNStruct *)GHSglobalcontext)->Binomial[j])*shapecoeff[cpos+j];
		}
		cpos+= numcoeff[i];
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;


	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to Start of main loop: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		
	}	


	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int t = 0;
	if(((MNStruct *)GHSglobalcontext)->SubIntToFit != -1){
		minep = ((MNStruct *)GHSglobalcontext)->SubIntToFit;
		maxep = minep+1;
		for(int ep = 0; ep < minep; ep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];
		}
	}
	

	double EQStart = 0;
	for(int ep = minep; ep < maxep; ep++){
	



		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];
		int NEpochBins = NChanInEpoch*GlobalNBins;


		double EpochChisq = 0;	
		double EpochdetN = 0;
		double EpochLike = 0;

		double EpochPhaseHess = 0;
		double EpochPhaseGrad = 0;

		int EDimcount = profdims + epochpriordims + ep*perEpochDims;
		int OriginalEDimCount = EDimcount;
		double EQUADSignal = 0;
			

	
		if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

			int EQIndex = 0;
                        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[t];}


			if(debug == 1){printf("EQUAD Signal: %i %i %g \n",ep, EDimcount,MaxAmps[EDimcount]);}
			EQUADSignal = MaxAmps[EDimcount]*EQUAD[EQIndex];
			EDimcount++;
		}

		double *OneEpochStochAmps = new double[totalshapestoccoeff];
		for(int i = 0; i < totalshapestoccoeff; i++){

			if(debug == 1){printf("Stoc Signal: %i %g \n",EDimcount,MaxAmps[EDimcount]);}

			OneEpochStochAmps[i] = MaxAmps[EDimcount]*StocProfCoeffs[i];
                        EDimcount++;	
		}
		
		



		for(int ch = 0; ch < NChanInEpoch; ch++){

			EDimcount = OriginalEDimCount;


			if(dotime == 1){	
				gettimeofday(&tval_before, NULL);
			}
			if(debug == 1){
				printf("In toa %i \n", t);
				printf("sat: %.15Lg \n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].sat);
			}


			int nTOA = t;


			

			double detN  = 0;
			double Chisq  = 0;
  

			double profilelike=0;

			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];
			long double FoldingPeriodDays = FoldingPeriod/SECDAY;

			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;

			double Tobs = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][2];
			double noiseval = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][3];
			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;


	

			double *shapevec  = new double[Nbins];
			double *ProfileModVec = new double[Nbins]();
			double *ProfileJitterModVec = new double[Nbins]();		   
			double *ProfileWidthModVec = new double[Nbins](); 
			long double binpos = ModelBats[nTOA]; 


			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
				binpos += (long double) EQUADSignal/SECDAY;
			}

			if(((MNStruct *)GHSglobalcontext)->incDM > 4 || ((MNStruct *)GHSglobalcontext)->yearlyDM == 1){
	                	double DMKappa = 2.410*pow(10.0,-16);
        		        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				long double DMshift = (long double)(SignalVec[ep]*DMScale);
				
		 		binpos+=DMshift/SECDAY;

			}

			if(binpos < ProfileBats[nTOA][0])binpos+=FoldingPeriodDays;
			if(binpos > ProfileBats[nTOA][1])binpos-=FoldingPeriodDays;


			if(binpos > ProfileBats[nTOA][1]){printf("OverBoard! %.10Lg %.10Lg %.10Lg %.10Lg\n", binpos, ModelBats[nTOA], ProfileBats[nTOA][1], (binpos-ProfileBats[nTOA][1])/FoldingPeriodDays);}

			long double minpos = binpos - FoldingPeriodDays/2;
			if(minpos < ProfileBats[nTOA][0])minpos=ProfileBats[nTOA][0];
			long double maxpos = binpos + FoldingPeriodDays/2;
			if(maxpos> ProfileBats[nTOA][1])maxpos =ProfileBats[nTOA][1];

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Interpolation Bin///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

			int InterpBin = 0;
			double FirstInterpTimeBin = 0;
			int  NumWholeBinInterpOffset = 0;

			if(((MNStruct *)GHSglobalcontext)->InterpolateProfile == 1){

		
				long double timediff = 0;
				long double bintime = ProfileBats[t][0];


				if(bintime  >= minpos && bintime <= maxpos){
				    timediff = bintime - binpos;
				}
				else if(bintime < minpos){
				    timediff = FoldingPeriodDays+bintime - binpos;
				}
				else if(bintime > maxpos){
				    timediff = bintime - FoldingPeriodDays - binpos;
				}

				timediff=timediff*SECDAY;

				double OneBin = FoldingPeriod/Nbins;
				int NumBinsInTimeDiff = floor(timediff/OneBin + 0.5);
				double WholeBinsInTimeDiff = NumBinsInTimeDiff*FoldingPeriod/Nbins;
				double OneBinTimeDiff = -1*((double)timediff - WholeBinsInTimeDiff);

				double PWrappedTimeDiff = (OneBinTimeDiff - floor(OneBinTimeDiff/OneBin)*OneBin);

				if(debug == 1)printf("Making InterpBin: %g %g %i %g %g %g\n", (double)timediff, OneBin, NumBinsInTimeDiff, WholeBinsInTimeDiff, OneBinTimeDiff, PWrappedTimeDiff);

				InterpBin = floor(PWrappedTimeDiff/((MNStruct *)GHSglobalcontext)->InterpolatedTime+0.5);
				if(InterpBin >= ((MNStruct *)GHSglobalcontext)->NumToInterpolate)InterpBin -= ((MNStruct *)GHSglobalcontext)->NumToInterpolate;

				FirstInterpTimeBin = -1*(InterpBin-1)*((MNStruct *)GHSglobalcontext)->InterpolatedTime;

				if(debug == 1)printf("Interp Time Diffs: %g %g %g %g \n", ((MNStruct *)GHSglobalcontext)->InterpolatedTime, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime, PWrappedTimeDiff, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime-PWrappedTimeDiff);

				double FirstBinOffset = timediff-FirstInterpTimeBin;
				double dNumWholeBinOffset = FirstBinOffset/(FoldingPeriod/Nbins);
				int  NumWholeBinOffset = 0;

				NumWholeBinInterpOffset = floor(dNumWholeBinOffset+0.5);
	
				if(debug == 1)printf("Interp bin is: %i , First Bin is %g, Offset is %i \n", InterpBin, FirstInterpTimeBin, NumWholeBinInterpOffset);


			}

	
			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed up to start of Interp: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freq - reffreq)/1000.0;
			double freqscale = pow(freqdiff, EvoFreqExponent);


			double snr = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].snr;
			double tobs = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].tobs; 


			snr = snr*3600/tobs;

			double refSN = 1000;
			double SNdiff =  snr/refSN;
			double SNscale = snr-refSN;


			double noisemean = MaxAmps[nTOA*perProfDims + 0];
			double MLAmp = MaxAmps[nTOA*perProfDims + 1];
			double MLSigma = MaxAmps[nTOA*perProfDims + 2];

			
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(int i =0; i < totalCoeffForMult; i++){
				ProfileModCoeffs[i]=0;	
			}				

			cpos = 0;
			epos = 0;
			fpos = 0;
			int spos = 0;
	
			
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				
				for(int i =0; i < numProfileFitCoeff[c]; i++){
					//printf("Making steps: PFD: %i %g \n",i, ProfileFitCoeffs[i+fpos]);
					ProfileModCoeffs[i+cpos] += ProfileFitCoeffs[i+fpos];

				}
				for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i =0; i < numEvoCoeff[c]; i++){
						ProfileModCoeffs[i+cpos] += EvoCoeffs[p][i+epos]*pow(freqscale, p+1);						
					}
				}


				for(int i = 0; i < numProfileStocCoeff[c]; i++){
					ProfileModCoeffs[i+cpos] += OneEpochStochAmps[i+spos]*modelflux;
					//printf("adding stoc: %i %i %i %g \n", ep, ch, i,OneEpochStochAmps[i+spos] );
				}
				

				cpos += NumCoeffForMult[c];
				epos += numEvoCoeff[c];
				fpos += numProfileFitCoeff[c];	
				spos += numProfileStocCoeff[c];
			}


			
			double *OneProfileParamHessian = new double[totalCoeffForMult*totalCoeffForMult];
			if(totalCoeffForMult > 0){
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ProfileModCoeffs,ProfileModVec,Nbins,totalCoeffForMult,'N');
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfileVec[InterpBin], ProfileModCoeffs,ProfileJitterModVec,Nbins,totalCoeffForMult,'N');
				vector_dgemm(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin] , OneProfileParamHessian, Nbins, totalCoeffForMult,Nbins, totalCoeffForMult, 'T', 'N');

				if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1 || ((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
					vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfileVec[InterpBin], ProfileModCoeffs,ProfileWidthModVec,Nbins,totalCoeffForMult,'N');
				}



			}


			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Modded Profile: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Fill Arrays with interpolated state//////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


				
			for(int j =0; j < Nbins; j++){

				double NewIndex = (j + NumWholeBinInterpOffset);
				int Nj =  GHSWrap(j + NumWholeBinInterpOffset, 0, Nbins-1);//(int)(NewIndex - floor(NewIndex/Nbins)*Nbins);


				//printf("Shape: %i %i %g %g \n", t, j, ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj], ProfileModVec[Nj] );
				shapevec[j] = ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj] + ProfileModVec[Nj];
			}


			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Filled Arrays: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////Get ML Amps and Baselines and Noise/////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




			double BaselineStep = 0;
			double AmpStep = 0;
			double SigmaStep = 0;

			int AmpZeroWrap = Nbins-GHSWrap((Nbins - NumWholeBinInterpOffset), 0, Nbins-1);
			int AmpZeroOffset = 0;//AmpZeroWrap%BinRatio;

			//printf("MLSigma %i %i %g\n", nTOA, nTOA*3, MLSigma);
			
			for(int i =0; i < ProfNbins; i++){

				double res = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1] - shapevec[AmpZeroOffset+i*BinRatio]*MLAmp - noisemean;
				//if(t==15657){printf("%i %i %g %g %g %g %g \n", t, i, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1], shapevec[AmpZeroOffset+i*BinRatio]*MLAmp + noisemean, res, 3*res*res/(MLSigma*MLSigma*MLSigma*MLSigma), 1.0/(MLSigma*MLSigma));}


				BaselineStep += 1.0/(MLSigma*MLSigma);	
				AmpStep += shapevec[AmpZeroOffset+i*BinRatio]*shapevec[AmpZeroOffset+i*BinRatio]/(MLSigma*MLSigma);	
				SigmaStep += 3*res*res/(MLSigma*MLSigma*MLSigma*MLSigma) - 1.0/(MLSigma*MLSigma);
				
				//Basline-Baseline Hessian
				SSM[nTOA][0 + 0*perProfDims] +=  1.0/(MLSigma*MLSigma);

				//Basline-Amp Hessian
				SSM[nTOA][0 + 1*perProfDims] +=  shapevec[AmpZeroOffset+i*BinRatio]/(MLSigma*MLSigma);	
				SSM[nTOA][1 + 0*perProfDims] +=  shapevec[AmpZeroOffset+i*BinRatio]/(MLSigma*MLSigma);	

				///Baseline-Sigma Hessian
				SSM[nTOA][2 + 0*perProfDims] += 2*res/(MLSigma*MLSigma*MLSigma);
				SSM[nTOA][0 + 2*perProfDims] += 2*res/(MLSigma*MLSigma*MLSigma);

				//Amp-Amp Hessian
				SSM[nTOA][1 + 1*perProfDims] += shapevec[AmpZeroOffset+i*BinRatio]*shapevec[AmpZeroOffset+i*BinRatio]/(MLSigma*MLSigma);

	
				//Amp-Sigma Hessian
				SSM[nTOA][2 + 1*perProfDims] += 2*shapevec[AmpZeroOffset+i*BinRatio]*res/(MLSigma*MLSigma*MLSigma);
				SSM[nTOA][1 + 2*perProfDims] += 2*shapevec[AmpZeroOffset+i*BinRatio]*res/(MLSigma*MLSigma*MLSigma);

				//Sigma-Sigma Hessian
				SSM[nTOA][2 + 2*perProfDims] += 3*res*res/(MLSigma*MLSigma*MLSigma*MLSigma) - 1.0/(MLSigma*MLSigma);
			
			}

		
			StepSize[nTOA*perProfDims + 0] = BaselineStep;
			StepSize[nTOA*perProfDims + 1] = AmpStep;
			StepSize[nTOA*perProfDims + 2] = SigmaStep;

//			printf("Max: %i %g %g %g \n", t, MaxAmps[nTOA*3 + 0], MaxAmps[nTOA*3 + 1], MaxAmps[nTOA*3 + 2] );
//			printf("Step: %i %g %g %g \n", t, BaselineStep, AmpStep, SigmaStep);


			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Up to LA: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Form Epoch Matrix//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			for(int i = 0; i < Nbins; i++){
				ProfileJitterModVec[i] += ((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i];
			}

			if(debug == 1){printf("Forming Epoch Matrix in MaxSteps\n");}
			int ZeroWrap = Nbins-GHSWrap((Nbins - NumWholeBinInterpOffset), 0, Nbins-1);
			int ZeroOffset = ZeroWrap%BinRatio;

			int HessCount = TotalProfs;
			int EpochMsize = globaldims + perEpochDims;

			double *EpochMNMHess = new double[EpochMsize*EpochMsize]();
			double *EpochMMatrix = new double[EpochMsize*Nbins]();
//			printf("allocated %i %i %i \n", EpochMsize, globaldims, perEpochDims);
			if(globaldims+epochdims+epochpriordims > 0){ 
//				printf("EQUAD1: %i %i\n", ((MNStruct *)GHSglobalcontext)->numFitEQUAD, t);

				double PhaseScale = -1*MLAmp/MLSigma;
				double ShapeScale = MLAmp*modelflux/MLSigma;

				int EpochPriorCount = 0;
				int EpochDimCount = epochpriordims+ep*perEpochDims;

				int EpochCount = 0;
				if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

//					printf("EQUAD: %i %i\n", ((MNStruct *)GHSglobalcontext)->numFitEQUAD, t);
					int EQIndex = 0;
					if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[t];}

					double MLEQUAD = pow(10.0, MaxAmps[profdims+EpochPriorCount+EQIndex]);
//					printf("ML EQUAD: %g %g \n", EQIndex, MaxAmps[profdims+EpochPriorCount+EQIndex], MLEQUAD);
					for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
			//			EpochMMatrix[i + EpochCount*Nbins] = MLEQUAD*PhaseScale*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i];
						EpochMMatrix[i + EpochCount*Nbins] = MLEQUAD*PhaseScale*ProfileJitterModVec[i];
					}

					EpochCount++;
					EpochPriorCount+=((MNStruct *)GHSglobalcontext)->numFitEQUAD;
					EpochDimCount++;				
				}
				//sleep(10);	

				if(totalshapestoccoeff > 0) {

					int cpos = 0;
					int spos = 0;

					double ShapeScale = MLAmp*modelflux/MLSigma;

					for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
						for(int i = 0; i < numProfileStocCoeff[c]; i++){

							double MLShapeP = pow(10.0, MaxAmps[profdims+EpochPriorCount]);
							for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
								EpochMMatrix[i + EpochCount*Nbins] = MLShapeP*ShapeScale*((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin][i];
							}

							EpochCount++;
							EpochPriorCount++;
							EpochDimCount++;
						}
						cpos += NumCoeffForMult[c];
						spos += numProfileStocCoeff[c];	
					}
				}

				for(int m = 0; m < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; m++){
					double phasefac = -1;
					if(m == 0){phasefac *= ((MNStruct *)GHSglobalcontext)->ReferencePeriod;}

					for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
						EpochMMatrix[i + EpochCount*Nbins] = GlobalDMatrix[nTOA][m]*ProfileJitterModVec[i]*PhaseScale*phasefac;

					}
					EpochCount++;
				}

				if(((MNStruct *)GHSglobalcontext)->incDM > 4){


					int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;
				        double DMKappa = 2.410*pow(10.0,-16);
        		  	        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));

					for(int m = 0; m < FitDMCoeff; m++){
						for(int i = ZeroOffset; i < Nbins; i+=BinRatio){

				//			EpochMMatrix[i + EpochCount*Nbins] = DMScale*PowerSpec[m]*FMatrix[ep + m*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i]*PhaseScale;
							EpochMMatrix[i + EpochCount*Nbins] = DMScale*PowerSpec[m]*FMatrix[ep + m*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*ProfileJitterModVec[i]*PhaseScale;




						}
						EpochCount++;
					}

				}


				for(int p = 0; p < 1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){

					double fval = pow(freqscale, p);
					int cpos = 0;
					int fpos = 0;

					for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
						for(int k = 0; k < numProfileFitCoeff[c]; k++){

							for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
								EpochMMatrix[i + EpochCount*Nbins] =  fval*((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin][i+(k+cpos)*Nbins]*MLAmp/MLSigma;
								//if(p==0 && c==0 && ProfNbins < Nbins){printf("Shape bit of matrix: %i %i %i %i %i %i %g \n", ep, ch, ZeroOffset, p, c, i,pow(freqscale, p)*((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin][i+c*Nbins]);}
							}
							EpochCount++;	
						}
						cpos += NumCoeffForMult[c];
						fpos +=	numProfileFitCoeff[c];	
						
					}	
					
				}

				if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1){
					for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
						EpochMMatrix[i + EpochCount*Nbins] = (ProfileWidthModVec[i] + ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][i])*MLAmp/MLSigma;
					}
					EpochCount++;
				}

				if(((MNStruct *)GHSglobalcontext)->incWidthEvoTime > 0){
					for(int ew = 0; ew < ((MNStruct *)GHSglobalcontext)->incWidthEvoTime; ew++){
						for(int i = ZeroOffset; i < Nbins; i+=BinRatio){
							EpochMMatrix[i + EpochCount*Nbins] = (ProfileWidthModVec[i] + ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][i])*pow(((((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].bat - ((MNStruct *)GHSglobalcontext)->ProfEvoTimeRef)/365.25), ew+1)*MLAmp/MLSigma;
						}
						EpochCount++;
					}
				}


				vector_dgemm(EpochMMatrix, EpochMMatrix , EpochMNMHess, Nbins, EpochMsize,Nbins, EpochMsize, 'T', 'N');

				EpochCount = 0;
				int GPSCount = 0;
				if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
					if(GlobalParmSet[GPSCount+((MNStruct *)GHSglobalcontext)->sysFlags[t]] > 0){
//						printf("Cut: %i %i \n", t, ((MNStruct *)GHSglobalcontext)->sysFlags[t]);
						for(int q = 0; q < EpochMsize; q++){
                                                        EpochMNMHess[EpochCount + q*EpochMsize] = 0;
                                                        EpochMNMHess[q + EpochCount*EpochMsize] = 0;
                                                }
                                        }
					else{
//						printf("No Cut: %i %i \n", t, ((MNStruct *)GHSglobalcontext)->sysFlags[t]);	
					}
					EpochCount++;
					GPSCount+=((MNStruct *)GHSglobalcontext)->numFitEQUAD;
				}


				for(int c = 0; c < totalshapestoccoeff; c++){
					if(GlobalParmSet[GPSCount] > 0){
						for(int q = 0; q < EpochMsize; q++){
                                                        EpochMNMHess[EpochCount + q*EpochMsize] = 0;
                                                        EpochMNMHess[q + EpochCount*EpochMsize] = 0;
                                                }
                                        }
                                        EpochCount++;
					GPSCount++;
				}
				
				EpochCount += ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;
				for(int m = 0; m < 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff; m++){
					if(GlobalParmSet[GPSCount] > 0){
                                                for(int q = 0; q < EpochMsize; q++){
                                                        EpochMNMHess[EpochCount + q*EpochMsize] = 0;
                                                        EpochMNMHess[q + EpochCount*EpochMsize] = 0;
                                                }
                                        }
                                        EpochCount++;
                                        GPSCount++;
                                }



/*


				for(int p = 0; p < perEpochDims; p++){
					if(GlobalParmSet[GPSCount] > 0.01){
						for(int q = 0; q < EpochMsize; q++){
							EpochMNMHess[p + q*EpochMsize] = 0;
							EpochMNMHess[q + p*EpochMsize] = 0;
						}
					}
				}
*/

				EpochCount=0;
				if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
					EpochMNMHess[EpochCount + EpochCount*EpochMsize] +=  1.0/NChanInEpoch;
					EpochCount++;	
				}
		
				for(int c = 0; c < totalshapestoccoeff; c++){
					EpochMNMHess[EpochCount + EpochCount*EpochMsize] +=   1.0/NChanInEpoch;
					EpochCount++;
				}

				EpochCount += ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;

				for(int m = 0; m < 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff; m++){
					//printf("DM: %i %g %g\n", m, EpochMNMHess[EpochCount + EpochCount*EpochMsize], 1.0/NChanInEpoch);
					EpochMNMHess[EpochCount + EpochCount*EpochMsize] +=   1.0/NChanInEpoch;
					
					EpochCount++;
				}

				EpochDimCount = epochpriordims+ep*perEpochDims;
				int globalDimCount = epochpriordims+epochdims;
				int thisHessSize = globaldims+epochdims+epochpriordims;

				for(int i = 0; i < perEpochDims; i++){
					for(int j = 0; j < perEpochDims; j++){
						SSM[HessCount][EpochDimCount+i + (EpochDimCount+j)*thisHessSize] += EpochMNMHess[i + j*EpochMsize];
					}
				}

				for(int i = 0; i < globaldims; i++){
					double crossfac = 1;
					for(int j = 0; j < globaldims; j++){
						double MVal = EpochMNMHess[perEpochDims+i + (perEpochDims+j)*EpochMsize];
						if(i != j){MVal *= crossfac;}
						SSM[HessCount][globalDimCount+i + (globalDimCount+j)*thisHessSize] += MVal;
					}
				}

				for(int i = 0; i < perEpochDims; i++){
					for(int j = 0; j < globaldims; j++){
						SSM[HessCount][EpochDimCount+i + (globalDimCount+j)*thisHessSize] += EpochMNMHess[i + (perEpochDims+j)*EpochMsize];
						SSM[HessCount][(globalDimCount+j) + (EpochDimCount+i)*thisHessSize] += EpochMNMHess[i + (perEpochDims+j)*EpochMsize];
						//if(j<6)printf("Timing Cross: %i %i %i %i %g \n", ep, ch, i, j,  EpochMNMHess[i + (perEpochDims+j)*EpochMsize]);
					}
				}
				

			}

			delete[] EpochMMatrix;

			if(debug == 1){printf("Made EPoch Matrix \n");}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Phase Term ////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		


		double PhaseScale = MLAmp/MLSigma;
		double OnePhaseHess = 0;
		double OnePhaseGrad = 0;

		double *NGRes = new double[Nbins]();
		for(int i =ZeroOffset; i < Nbins; i+=BinRatio){
		//	OnePhaseHess += PhaseScale*PhaseScale*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i]*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i];
			OnePhaseHess += PhaseScale*PhaseScale*ProfileJitterModVec[i]*ProfileJitterModVec[i];


			int RevWrap = GHSWrap(i + (Nbins - NumWholeBinInterpOffset), 0, Nbins-1);
			double gres = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][RevWrap/BinRatio][1] - MLAmp*shapevec[RevWrap] - noisemean;
//			OnePhaseGrad += -1*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][i]*PhaseScale*gres/MLSigma;
			OnePhaseGrad += -1*ProfileJitterModVec[i]*PhaseScale*gres/MLSigma;



			NGRes[i] = gres/MLSigma/MLSigma;

		}
		if(debug == 1){printf("Made Phase Term\n");}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Epoch Terms////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double EQUADGrad = 0;
		int EpochCount = 0;
		if(epochpriordims > 0){ 

			
			int EpochPriorCount = 0;
			int EpochDimCount = epochpriordims+ep*perEpochDims;
			
			double MLEQUAD = 0;
			double MLEQSignal = 0;
			int EQUADDim = 0;
			int EQSignalDim = 0;

			int EQIndex = 0;
                        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[t];}


			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

				EQUADDim = EpochPriorCount + EQIndex;
				EQSignalDim = EpochDimCount;
				MLEQUAD = pow(10.0, MaxAmps[profdims+EpochPriorCount+EQIndex]);
				MLEQSignal = MaxAmps[profdims+EpochDimCount];
				

				double PriorTerm = log(10.0)*log(10.0)*(MLEQSignal*MLEQUAD*OnePhaseHess*MLEQSignal*MLEQUAD);// - MLEQSignal*MLEQUAD*OnePhaseGrad);
				if(GlobalParmSet[EpochCount+EQIndex] > 0){PriorTerm = GlobalParmSet[EpochCount+EQIndex];}

				SSM[HessCount][EQUADDim + EQUADDim*(epochdims+epochpriordims+globaldims)] +=  PriorTerm;


				
				double PriorCrossTerm = log(10.0)*(MLEQUAD*OnePhaseHess*MLEQUAD*MLEQSignal - MLEQUAD*OnePhaseGrad);
				
				if(GlobalParmSet[EpochCount+EQIndex] > 0){PriorCrossTerm = 0;}
				


				SSM[HessCount][EQUADDim + EpochDimCount*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
				SSM[HessCount][EpochDimCount + EQUADDim*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;





				//////////////////////////////////////////////Add in Cross Terms for Timing Model///////////////////////////////////////////////////


				
				for(int m = 1; m < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; m++){
					double phasefac = -1;
					if(m == 0){phasefac *= ((MNStruct *)GHSglobalcontext)->ReferencePeriod;}


					double TimingCrossTerm = 0;//log(10.0)*MLEQUAD*MLEQSignal*OnePhaseHess*GlobalDMatrix[nTOA][m]*phasefac;

					if(GlobalParmSet[EpochCount+EQIndex] > 0){
						TimingCrossTerm = 0;
					}


					int TimingDim = epochpriordims+epochdims+m;
					SSM[HessCount][EQUADDim + TimingDim*(epochpriordims+epochdims+globaldims)] += 0.5*TimingCrossTerm;
        	                        SSM[HessCount][TimingDim + EQUADDim*(epochpriordims+epochdims+globaldims)] += 0.5*TimingCrossTerm;
	


				}

				EpochPhaseHess += OnePhaseHess;
				EpochPhaseGrad += EQUADGrad;


				EpochCount+=((MNStruct *)GHSglobalcontext)->numFitEQUAD;
				EpochPriorCount+=((MNStruct *)GHSglobalcontext)->numFitEQUAD;
				EpochDimCount++;
			}

			if(totalshapestoccoeff > 0) {

				double *BNRes = new double[totalCoeffForMult];
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], NGRes, BNRes,Nbins,totalCoeffForMult,'T');

				double *PhaseShapeCross = new double[totalCoeffForMult];

//				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin],PhaseShapeCross,Nbins,totalCoeffForMult,'T');
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ProfileJitterModVec ,PhaseShapeCross,Nbins,totalCoeffForMult,'T');



				int cpos = 0;
				int spos = 0;

				for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
					for(int i = 0; i <numProfileStocCoeff[c]; i++){

						int crossScale = 0.5;

						double OneShapeGrad = BNRes[i+cpos]*MLAmp*modelflux;
						double OneShapeHess = OneProfileParamHessian[i+cpos + (i+cpos)*totalCoeffForMult]*MLAmp*MLAmp*modelflux*modelflux/(MLSigma*MLSigma);
					
						double MLShapeP = pow(10.0, MaxAmps[profdims+EpochPriorCount]);
						double MLShapeSignal = MaxAmps[profdims+EpochDimCount];

						double PriorTerm = log(10.0)*log(10.0)*(MLShapeSignal*MLShapeP*OneShapeHess*MLShapeSignal*MLShapeP);// - MLShapeSignal*MLShapeP*OneShapeGrad);

						if(t==0){printf("New p: %i %g \n", i,  GlobalParmSet[EpochCount]);}
						if(GlobalParmSet[EpochCount] > 0){if(t==0){printf("New p: %i %g \n", i,  GlobalParmSet[EpochCount]);}PriorTerm = GlobalParmSet[EpochCount]/TotalProfs;}

					
						SSM[HessCount][EpochPriorCount + EpochPriorCount*(epochdims+epochpriordims+globaldims)] += PriorTerm;

						double PriorCrossTerm = log(10.0)*(MLShapeP*OneShapeHess*MLShapeP*MLShapeSignal - MLShapeP*OneShapeGrad);

						if(GlobalParmSet[EpochCount] > 0){PriorCrossTerm = 0;}

				
						SSM[HessCount][EpochPriorCount + EpochDimCount*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
						SSM[HessCount][EpochDimCount + EpochPriorCount*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;

						if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0 && GlobalParmSet[EpochCount] == 0 && GlobalParmSet[EQIndex] == 0){

							//Need -ve for phase component as with gradient as sign is backwards
							double CrossFac = (-MLAmp)*MLAmp*modelflux/(MLSigma*MLSigma);

							double EQAmpShapePCrossTerm = log(10.0)*CrossFac*MLEQUAD*PhaseShapeCross[i+cpos]*MLShapeP*MLShapeSignal;



							SSM[HessCount][EQSignalDim + EpochPriorCount*(epochpriordims+epochdims+globaldims)] += crossScale*EQAmpShapePCrossTerm;
							SSM[HessCount][EpochPriorCount + EQSignalDim*(epochpriordims+epochdims+globaldims)] += crossScale*EQAmpShapePCrossTerm;

							double ShapeAmpEQPCrossTerm = log(10.0)*CrossFac*MLEQSignal*MLEQUAD*PhaseShapeCross[i+cpos]*MLShapeP;

							SSM[HessCount][EpochDimCount + EQUADDim*(epochpriordims+epochdims+globaldims)] += crossScale*ShapeAmpEQPCrossTerm;
							SSM[HessCount][EQUADDim + EpochDimCount*(epochpriordims+epochdims+globaldims)] += crossScale*ShapeAmpEQPCrossTerm;
//

						//	double ShapePEQPCrossTerm = log(10.0)*log(10.0)*CrossFac*MLEQSignal*MLEQUAD*PhaseShapeCross[i+cpos]*MLShapeP*MLShapeSignal;

						//	SSM[HessCount][EpochPriorCount + EQUADDim*(epochpriordims+epochdims)] = ShapePEQPCrossTerm;
						//	SSM[HessCount][EQUADDim + EpochPriorCount*(epochpriordims+epochdims)] = ShapePEQPCrossTerm;

	


						}
						


						EpochPriorCount++;
						EpochDimCount++;
						EpochCount++;

					}

					cpos += NumCoeffForMult[c];
					spos += numProfileStocCoeff[c];
				}

				delete[] PhaseShapeCross;
				delete[] BNRes;
			}

			if(debug == 1){printf("Done EPoch Terms\n");}


			if(((MNStruct *)GHSglobalcontext)->incDM > 4){

				int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;
			        double DMKappa = 2.410*pow(10.0,-16);
		  	        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				double DMshift = SignalVec[ep]*DMScale;


				

				

				/////////////Hessian for the Amplitude////////////////////


				double PriorTerm = log(10.0)*log(10.0)*(DMshift*OnePhaseHess*DMshift - DMshift*OnePhaseGrad);
				if(GlobalParmSet[EpochCount] > 0){PriorTerm = GlobalParmSet[EpochCount]/TotalProfs;}

				SSM[HessCount][EpochPriorCount + EpochPriorCount*(epochdims+epochpriordims+globaldims)] +=  PriorTerm;

				int DMAmpDim = profdims + epochdims + epochpriordims + ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;
				int DMAmpHessDim = epochdims + epochpriordims + ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;


				for(int m = 0; m < FitDMCoeff; m++){

					double MLDMSignal = MaxAmps[DMAmpDim+m];

					double PF = FMatrix[ep + m*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[m]*DMScale;
				
					double GradTerm = -PF*OnePhaseGrad;
					double HessTerm = PF*OnePhaseHess*PF*MLDMSignal;

					double PriorCrossTerm = log(10.0)*(HessTerm + GradTerm);
				
					if(GlobalParmSet[EpochCount] > 0){PriorCrossTerm = 0;}
				


					SSM[HessCount][EpochPriorCount + (DMAmpHessDim+m)*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
					SSM[HessCount][(DMAmpHessDim+m) + EpochPriorCount*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
				}

				EpochCount++;
				EpochPriorCount++;

				/////////////Hessian for the Spectral Index////////////////////


				double DMSpecHessTerm = 0;
				double DMSpecAmpHessTerm = 0;
				for(int m = 0; m < FitDMCoeff/2; m++){

					double freq = ((double)(m+1.0))/((MNStruct *)GHSglobalcontext)->Tspan;

					double MLDMSignal1 = MaxAmps[DMAmpDim+m];
					double KPF1 = FMatrix[ep + m*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[m]*DMScale;
					DMSpecHessTerm += (-0.5*log(freq*365.25))*KPF1*MLDMSignal1*OnePhaseHess*KPF1*MLDMSignal1*(-0.5*log(freq*365.25)) - KPF1*MLDMSignal1*(0.25*log(freq*365.25)*log(freq*365.25))*OnePhaseGrad;


					DMSpecAmpHessTerm += (-0.5*log(freq*365.25))*KPF1*MLDMSignal1*OnePhaseHess*KPF1*MLDMSignal1*(log(10.0)) - KPF1*MLDMSignal1*(-0.5*log(freq*365.25)*log(10.0))*OnePhaseGrad;

					double MLDMSignal2 = MaxAmps[DMAmpDim+m+FitDMCoeff/2];
	 			        double KPF2 = FMatrix[ep + (m+FitDMCoeff/2)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[m+FitDMCoeff/2]*DMScale;
					DMSpecHessTerm += (-0.5*log(freq*365.25))*KPF2*MLDMSignal2*OnePhaseHess*KPF2*MLDMSignal2*(-0.5*log(freq*365.25)) - KPF2*MLDMSignal2*(0.25*log(freq*365.25)*log(freq*365.25))*OnePhaseGrad;

					DMSpecAmpHessTerm += (-0.5*log(freq*365.25))*KPF2*MLDMSignal2*OnePhaseHess*KPF2*MLDMSignal2*(log(10.0)) - KPF2*MLDMSignal2*(-0.5*log(freq*365.25)*log(10.0))*OnePhaseGrad;
				
					double GradTerm = -1*(-0.5*log(freq*365.25))*KPF1*OnePhaseGrad;
					double HessTerm =  (-0.5*log(freq*365.25))*KPF1*OnePhaseHess*KPF1*MLDMSignal1;

					double PriorCrossTerm = (HessTerm + GradTerm);
				
					if(GlobalParmSet[EpochCount] > 0){PriorCrossTerm = 0;}
	
					SSM[HessCount][EpochPriorCount + (DMAmpHessDim+m)*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
					SSM[HessCount][(DMAmpHessDim+m) + EpochPriorCount*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;


					GradTerm = -1*(-0.5*log(freq*365.25))*KPF2*OnePhaseGrad;
					HessTerm =  (-0.5*log(freq*365.25))*KPF2*OnePhaseHess*KPF2*MLDMSignal2;

					PriorCrossTerm = (HessTerm + GradTerm);
				
					if(GlobalParmSet[EpochCount] > 0){PriorCrossTerm = 0;}
	
					SSM[HessCount][EpochPriorCount + (DMAmpHessDim+m+FitDMCoeff/2)*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
					SSM[HessCount][(DMAmpHessDim+m+FitDMCoeff/2) + EpochPriorCount*(epochpriordims+epochdims+globaldims)] += crossScale*PriorCrossTerm;
				}

				PriorTerm = DMSpecHessTerm;

				if(GlobalParmSet[EpochCount] > 0){PriorTerm = GlobalParmSet[EpochCount]/TotalProfs;}
				SSM[HessCount][EpochPriorCount + EpochPriorCount*(epochdims+epochpriordims+globaldims)] +=  PriorTerm;

				if(GlobalParmSet[EpochCount] > 0){DMSpecAmpHessTerm = 0;}

				SSM[HessCount][EpochPriorCount-1 + EpochPriorCount*(epochdims+epochpriordims+globaldims)] +=  DMSpecAmpHessTerm;
				SSM[HessCount][EpochPriorCount + (EpochPriorCount-1)*(epochdims+epochpriordims+globaldims)] +=  DMSpecAmpHessTerm;

				EpochCount++;
				EpochPriorCount++;
				



			}

			HessCount += ((MNStruct *)GHSglobalcontext)->numProfileEpochs+1;

		}


			profilelike = -0.5*(detN + Chisq);

			if(debug == 1)printf("Like: %i %.15g %.15g %.15g \n", nTOA, lnew, detN, Chisq);

			EpochChisq+=Chisq;
			EpochdetN+=detN;
			EpochLike+=profilelike;



			
			delete[] EpochMNMHess;


			
			delete[] OneProfileParamHessian;
			delete[] shapevec;
			delete[] ProfileModVec;
			delete[] ProfileJitterModVec;
			delete[] ProfileWidthModVec;
			delete[] NGRes;

			if(dotime == 1){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  End of Epoch: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
			}

			

			t++;
		}



///////////////////////////////////////////

		if(dotime == 1){
			gettimeofday(&tval_before, NULL);
		}

		if(((MNStruct *)GHSglobalcontext)->incWideBandNoise == 0){

			lnew += EpochLike;
		}

		delete[] OneEpochStochAmps;


////////////////////////////////////////////
	
	}

	
	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->pulse->nobs; j++){
	    delete[] ProfileBats[j];
	}
	delete[] ProfileBats;
	delete[] ModelBats;
	delete[] numcoeff;
	delete[] NumCoeffForMult;
	delete[] numShapeToSave;

	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] EvoCoeffs[j];
	}
	delete[] EvoCoeffs;
	delete[] betas;
	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){ delete[] EQUAD; }

	if(debug == 1)printf("End Like: %.10g \n", lnew);



	int HessCount = 0;
	int dimCount = 0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Get Profile EVDs///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		double *PEVals = new double[perProfDims]();
		for(int i = 0; i < TotalProfs; i++){
			PosDefEVD(SSM[i], PEVals, perProfDims);
			for(int j = 0; j < perProfDims; j++){
				StepSize[i*perProfDims + j] = PEVals[j];
			}
		}
		delete[] PEVals;
	}
	HessCount += TotalProfs;
	dimCount += profdims;




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Get Epoch Prior EVDs///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	if(globaldims+epochdims+epochpriordims > 0){ 


		double expansionfactor = ((MNStruct *)GHSglobalcontext)->phasePriorExpansion;
		((MNStruct *)GHSglobalcontext)->PhasePrior = (1.0/sqrt(SSM[HessCount][epochpriordims+epochdims + (epochpriordims+epochdims)*(epochpriordims+epochdims+globaldims)]))*expansionfactor;
		printf("using prior on phase: %g %g \n", ((MNStruct *)GHSglobalcontext)->PhasePrior, ((MNStruct *)GHSglobalcontext)->phasePriorExpansion);
		SSM[HessCount][epochpriordims+epochdims + (epochpriordims+epochdims)*(epochpriordims+epochdims+globaldims)] += 1.0/(((MNStruct *)GHSglobalcontext)->PhasePrior*((MNStruct *)GHSglobalcontext)->PhasePrior);



		//for(int j = 0; j < epochpriordims; j++){
		//	SSM[HessCount][j + j*epochpriordims] = SSM[HessCount][j + j*epochpriordims]/sqrt(1.0*((MNStruct *)GHSglobalcontext)->numProfileEpochs);
		//}


		for(int j = epochdims+epochpriordims; j < globaldims+epochdims+epochpriordims; j++){
			for(int k = 0; k < epochpriordims+globaldims+epochdims; k++){
//				printf("A[%i][%i] = %g \n", j,k, SSM[HessCount][j + k*(epochpriordims+epochdims+globaldims)]);	
			}
		}

		printf("Prior Step: %g \n", SSM[HessCount][epochpriordims+epochdims + (epochpriordims+epochdims)*(epochpriordims+epochdims+globaldims)]);
		if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
			double *PEVals = new double[epochpriordims+epochdims+globaldims]();
			PosDefEVD(SSM[HessCount], PEVals, epochpriordims+epochdims+globaldims);
			for(int j = 0; j < epochpriordims+epochdims+globaldims; j++){
//				printf("Prior Step2: %g \n", PEVals[j]);
				StepSize[dimCount + j] = PEVals[j];

			}
			delete[] PEVals;
		}
		else{

			for(int j = 0; j < epochpriordims+epochdims+globaldims; j++){
				StepSize[dimCount + j] = SSM[HessCount][j + j*(epochpriordims+epochdims+globaldims)];
				printf("prior step %i %g \n", j, SSM[HessCount][j + j*(epochpriordims+epochdims+globaldims)]);
			}

		}
		HessCount += 1+((MNStruct *)GHSglobalcontext)->numProfileEpochs;
		dimCount += epochpriordims+epochdims+globaldims;
	}



//Dont use below
/*

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Get Epoch EVDs/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	if(epochdims > 0){ 

		if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
			double *PEVals = new double[perEpochDims]();
			for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfileEpochs; i++){
				PosDefEVD(SSM[HessCount], PEVals, perEpochDims);
				for(int j = 0; j < perEpochDims; j++){
					//printf("Step comp %i %i %g %g \n", i, j, StepSize[i*3 + j], PEVals[j]);
					StepSize[dimCount + i*perEpochDims + j] = PEVals[j];
					//MaxAmps[dimCount + i*perEpochDims + j] = 0;
				}
				HessCount++;
			}
			delete[] PEVals;
		}
		else{
			for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfileEpochs; i++){
				for(int j = 0; j < perEpochDims; j++){
					StepSize[dimCount + i*perEpochDims + j] = SSM[HessCount][j + j*perEpochDims];//1;
					printf("Actual Epoch Dim Hessian: %i %i %g %g\n", i, j, SSM[HessCount][j + j*perEpochDims], 1.0/sqrt(SSM[HessCount][j + j*perEpochDims]) );
					//MaxAmps[dimCount + i*perEpochDims + j] = 0;
				}
				HessCount++;
			}
		}


		dimCount += epochdims;
	}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// Finish Shape parameters////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int vsize = ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps + totalCoeffForMult*(1+((MNStruct *)GHSglobalcontext)->NProfileEvoPoly);
	if(vsize > 0){

		//SSM[TotalProfs][0] += 1.0/(4.2e-05)/(4.2e-05);
		double expansionfactor = 1;
		((MNStruct *)GHSglobalcontext)->PhasePrior = (1.0/sqrt(SSM[HessCount][0]))*expansionfactor;
		printf("using prior on phase: %g \n", ((MNStruct *)GHSglobalcontext)->PhasePrior);
		SSM[HessCount][0] += 1.0/(((MNStruct *)GHSglobalcontext)->PhasePrior*((MNStruct *)GHSglobalcontext)->PhasePrior);

		for(int i = 0; i < vsize; i++){
			for(int j = i; j < vsize; j++){
				//printf("A[%i][%i] = %g\n", i, j, SSM[HessCount][i + j*vsize]);
			}

		}

		if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){

			double *PPHEVals = new double[vsize]();
		//	printf("calling EVD with %i \n", vsize);
			PosDefEVD(SSM[HessCount], PPHEVals, vsize);



			for(int i = 0; i < vsize; i++){
				//MaxAmps[dimCount + i] = 0;	
				StepSize[dimCount + i] = PPHEVals[i];	
				//printf("Shape steps: %i %g %g \n", t*3+i, MaxAmps[t*3 + i], StepSize[t*3 + i]);

			}

			delete[] PPHEVals;
		}
		else{

			for(int i = 0; i < vsize; i++){
				//MaxAmps[dimCount + i] = 0;
				StepSize[dimCount + i] = SSM[HessCount][i + i*vsize];
				//if(i==1){StepSize[dimCount + i] = 30;}
				
			}
		}


		HessCount += 1;
		dimCount += vsize;

	}

*/
//Dont use above

	printf("Finished with steps\n");

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

		double *oneGrad = new double[DimsPerProf];
		double *RGrad = new double[DimsPerProf];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->TotalProfiles; i++){
			for(int j = 0; j < DimsPerProf; j++){
				oneGrad[j] = gradPhy[i*DimsPerProf + j];
			}
			vector_dgemv(GlobalHessian[i], oneGrad,RGrad, DimsPerProf, DimsPerProf,'T');
			for(int j = 0; j < DimsPerProf; j++){
				gradPrin[i*DimsPerProf + j] = RGrad[j];
			}
		}


		delete[] oneGrad;	
		delete[] RGrad;

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


/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Then Epoch Parameters////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(DimsPerEpoch > 0){

		double *oneGrad = new double[DimsPerEpoch];
		double *RGrad = new double[DimsPerEpoch];
	
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfileEpochs; i++){
			for(int j = 0; j < DimsPerEpoch; j++){
				oneGrad[j] = gradPhy[DimCount + j];
			}
			vector_dgemv(GlobalHessian[HessCount], oneGrad,RGrad, DimsPerEpoch, DimsPerEpoch,'T');
			for(int j = 0; j < DimsPerEpoch; j++){
				gradPrin[DimCount + j] = RGrad[j];
			}

			HessCount++;
			DimCount += DimsPerEpoch;
		}
	
		delete[] oneGrad;	
		delete[] RGrad;


	}
*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Then Global Parameters////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

	if(globaldims > 0){

		double *oneGrad = new double[globaldims];
		double *RGrad = new double[globaldims];

		for(int j = 0; j < globaldims; j++){
			oneGrad[j] = gradPhy[DimCount + j];
		}

		vector_dgemv(GlobalHessian[HessCount], oneGrad, RGrad,  globaldims,  globaldims,'T');

		for(int j = 0; j < globaldims; j++){
			gradPrin[DimCount + j] = RGrad[j];
		}

		delete[] oneGrad;	
		delete[] RGrad;

	}

*/	
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

		double *oneGrad = new double[DimsPerProf];
		double *RGrad = new double[DimsPerProf];
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->TotalProfiles; i++){
			for(int j = 0; j < DimsPerProf; j++){
				oneGrad[j] = xPrin[i*DimsPerProf + j];
			}
			vector_dgemv(GlobalHessian[i], oneGrad,RGrad, DimsPerProf, DimsPerProf,'N');
			for(int j = 0; j < DimsPerProf; j++){
				xPhy[i*DimsPerProf + j] = RGrad[j];
			}
		}


		delete[] oneGrad;	
		delete[] RGrad;

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
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Then Epoch Parameters////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(DimsPerEpoch > 0){

		double *oneGrad = new double[DimsPerEpoch];
		double *RGrad = new double[DimsPerEpoch];
	
		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfileEpochs; i++){
			for(int j = 0; j < DimsPerEpoch; j++){
				oneGrad[j] = xPrin[DimCount + j];
			}
			vector_dgemv(GlobalHessian[HessCount], oneGrad,RGrad, DimsPerEpoch, DimsPerEpoch,'T');
			for(int j = 0; j < DimsPerEpoch; j++){
				xPhy[DimCount + j] = RGrad[j];
			}

			HessCount++;
			DimCount += DimsPerEpoch;
		}
	
		delete[] oneGrad;	
		delete[] RGrad;



	}


*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Then Global Parameters////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
	if(globaldims > 0){

		double *oneGrad = new double[globaldims];
		double *RGrad = new double[globaldims];

		for(int j = 0; j < globaldims; j++){
			oneGrad[j] = xPrin[DimCount + j];
		}

		vector_dgemv(GlobalHessian[HessCount], oneGrad,RGrad,  globaldims, globaldims,'N');

		for(int j = 0; j < globaldims; j++){
			xPhy[DimCount + j] = RGrad[j];
		}

		delete[] oneGrad;	
		delete[] RGrad;

	}

*/

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

	int fitpsum=0;
	int nofitpsum=0;
	for(int p2 = 0; p2 < ((MNStruct *)GHSglobalcontext)->numProfComponents; p2++){
		int skipone=0;
		if(p2==0 && nofitpsum == 0){
			printf("P0: %i %g %g\n", nofitpsum, ((MNStruct *)GHSglobalcontext)->MeanProfileShape[nofitpsum], 0.0);
			profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileShape[nofitpsum] <<"\n";
			skipone=1;
			nofitpsum++;
		}
		for(int p3 = 0; p3 < ((MNStruct *)GHSglobalcontext)->numProfileFitCoeff[p2]; p3++){
			printf("PF: %i %g %g\n", nofitpsum+p3, ((MNStruct *)GHSglobalcontext)->MeanProfileShape[nofitpsum+p3]+pdParameterEstimates[pcount], pdParameterEstimates[pcount]);
			profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileShape[nofitpsum+p3]+pdParameterEstimates[pcount] <<"\n";
			pcount++;
			
		}
		for(int p3 = ((MNStruct *)GHSglobalcontext)->numProfileFitCoeff[p2]; p3 < ((MNStruct *)GHSglobalcontext)->numshapecoeff[p2]-skipone; p3++){
			printf("PNF: %i %g\n", nofitpsum+p3, ((MNStruct *)GHSglobalcontext)->MeanProfileShape[nofitpsum+p3]);
			profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileShape[nofitpsum+p3] <<"\n";
		}
		fitpsum+=((MNStruct *)GHSglobalcontext)->numProfileFitCoeff[p2];
		nofitpsum+=((MNStruct *)GHSglobalcontext)->numshapecoeff[p2]-skipone;
	}

	for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
		printf("B: %g \n", ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[c]);
		profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[c] <<"\n";
	}
	
	double **EvoCoeffs=new double*[((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){EvoCoeffs[i] = new double[((MNStruct *)GHSglobalcontext)->totalEvoCoeff];}

	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		for(int i =0; i < ((MNStruct *)GHSglobalcontext)->totalEvoCoeff; i++){
			EvoCoeffs[p][i]=((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p][i];
		}
	}


	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		int cpos = 0;
		for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
			for(int i =0; i < ((MNStruct *)GHSglobalcontext)->numEvoFitCoeff[c]; i++){
				EvoCoeffs[p][i+cpos] += pdParameterEstimates[pcount];
				printf("Evo check: %i %i %i %i %g \n",p, c, i, pcount,  pdParameterEstimates[pcount]); 
				pcount++;
			}
			cpos += ((MNStruct *)GHSglobalcontext)->numEvoCoeff[c];
		
		}
	}
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		for(int i =0; i < ((MNStruct *)GHSglobalcontext)->totalEvoCoeff; i++){
			printf("E: %g \n", EvoCoeffs[p][i]);
			profilefile << std::setprecision(10) << EvoCoeffs[p][i] <<"\n";
		}
	}
	for(int i =0; i < ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff; i++){
			printf("S: %g \n", ((MNStruct *)GHSglobalcontext)->MeanProfileStoc[i]);
			profilefile << std::setprecision(10) << ((MNStruct *)GHSglobalcontext)->MeanProfileStoc[i] <<"\n";	
	}
	profilefile.close();
	printf("   Max Like %.10g \n", MLike);
	TNtextOutput(((MNStruct *)GHSglobalcontext)->pulse, 1, 0, Tempo2Fit,  GHSglobalcontext, 0, nParameters, paramlist, 0.0, 0, 0, 0, longname, paramarray);


}

static lbfgsfloatval_t evaluate( void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){


	int ndim2 = n;
        int *ndim = &ndim2;
	double 	like = 0;
	double mylike = 0;
	
	double *Myx = new double[n];
        double *Myg = new double[n];

         for (int i = 0;i < n;i++) {
		if(i >= ((MNStruct *)GHSglobalcontext)->TotalProfiles*3){
			//printf("stepsize? %i %g \n", i, GlobalStepSize[i]);
		}
                Myx[i] = GlobalStartPoint[i] + x[i]*GlobalStepSize[i];
                Myg[i] = 0;
        }
	GHSProfileDomainLike2(ndim, Myx, mylike, Myg);



	printf("in function %g %g %.16g\n", x[0], Myx[0], mylike);
        
/*        lbfgsfloatval_t fx = 0.0;

    for (int i = 0;i < n; i += 2) {
        double t1 = 1.0 - Myx[i];
        double t2 = 10.0 * (Myx[i+1] - Myx[i] * Myx[i]);
        Myg[i+1] = 20.0 * t2;
        Myg[i] = -2.0 * (Myx[i] * Myg[i+1] + t1);
        like += t1 * t1 + t2 * t2;
    
	}*/

         for (int i = 0;i < n;i++) {
//		printf("grad %i %g %g %g\n", i, Myg[i], GlobalStepSize[i],  Myg[i]*GlobalStepSize[i]);
		g[i] = Myg[i]*GlobalStepSize[i];

		if(isnan(Myg[i])){printf("Grad %i was a nan \n", i); g[i] = 0;}
		if(isnan(GlobalStepSize[i])){printf("Step %i was a nan \n", i); g[i] = 0;}
        }

        delete[] Myg;
        delete[] Myx;

    return mylike;
}

/*
static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    double  *fx;
	int testn = n;
	int *testdim = &testn;
	int *ndim = testdim;
	double *Myx = new double[n];
	double *Myg = new double[n];

	 for (int i = 0;i < n;i++) {
		Myx[i] = 0;
		Myg[i] = 0;
	}

    lbfgsfloatval_t thisfx = 0.0;
	
    for (int i = 0;i < *ndim;i += 2) {
       // lbfgsfloatval_t t1 = 1.0 - x[i];
        //lbfgsfloatval_t t2 = 10.0 * (x[i+1] - x[i] * x[i]);
        //g[i+1] = 20.0 * t2;
        //g[i] = -2.0 * (x[i] * g[i+1] + t1);
        //thisfx += t1 * t1 + t2 * t2;
    }

//	GHSProfileDomainLike(ndim, Myx, fx, Myg);

    for (int i = 0;i < n;i ++) {
        //g[i] = 0;//Myg[i];
    }
	printf("ML: %g \n", *fx);
	delete[] Myg;
	delete[] Myx;

	printf("done \n");
    return thisfx;
}
*/
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
    )
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f \n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}


void GHSProfileDomainLike2(int* ndim, double* PrinCube, double &likelihood, double* PrinGrad){

	int dotime = 0;
	int debug = ((MNStruct *)GHSglobalcontext)->debug; 
	int testdim = 6;
	double wtime = 0;

	int testuniformpriors = 0;


	int threads = omp_get_max_threads();

   	

	struct timeval tval_before, tval_after, tval_resultone, tval_resulttwo;

	if(dotime == 1){
		gettimeofday(&tval_before, NULL);
		wtime = omp_get_wtime ();
	}
	
	if(debug == 1)printf("In like \n");

	double *grad = new double[*ndim];
	for(int i = 0; i < *ndim; i++){
		//if(i!=testdim){Cube[i]=GlobalStartPoint[i];}
		//if(i==testdim){Cube[i]=100;}
		grad[i] = 0;
	}
	
	double *Cube = new double[*ndim];

	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Physical(PrinCube, Cube);
	}
	else{
		for(int i = 0; i < *ndim; i++){
			Cube[i] = PrinCube[i];
		}
	}



	if(dotime == 3){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to First Transform: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);
	}	


	//Cube[((MNStruct *)GHSglobalcontext)->TotalProfiles*3] = -6.45022;
	//Cube[((MNStruct *)GHSglobalcontext)->TotalProfiles*3+1] = -1.1342;

	//for(int i = ((MNStruct *)GHSglobalcontext)->TotalProfiles*3; i < *ndim; i++){
	//	printf("change: %i %g %g \n", i, PrinCube[i], Cube[i]);
	//}
	//sleep(5);



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get dimensionality//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

        int DimsPerProf = ((MNStruct *)GHSglobalcontext)->GHSperProfDims;
        int DimsPerEpoch = ((MNStruct *)GHSglobalcontext)->GHSperEpochDims;
        int epochpriordims = ((MNStruct *)GHSglobalcontext)->GHSepochpriordims;
        int globaldims = ((MNStruct *)GHSglobalcontext)->GHSglobaldims;

        int profdims = ((MNStruct *)GHSglobalcontext)->TotalProfiles*DimsPerProf;
        int epochdims = ((MNStruct *)GHSglobalcontext)->numProfileEpochs*DimsPerEpoch;





	
	int globalparams = globaldims;

	for(int i = 0; i < profdims+epochpriordims; i++){
		Cube[i]=GlobalStartPoint[i];
	}


	for(int i =  profdims + epochdims + epochpriordims + ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; i < *ndim; i++){
		//Cube[i]=GlobalStartPoint[i];
	}

	double **threadgrads = new double*[globalparams];
	for(int i = 0; i < globalparams; i++){ threadgrads[i] = new double[threads]();}


	double **PriorThreadGrads = new double*[epochpriordims];
	for(int i = 0; i < epochpriordims; i++){ PriorThreadGrads[i] = new double[threads]();}



	double uniformpriorterm = 0;
	double priorterm = 0;
	double priorDet = 0;
	double JacobianTerm = 0;


	//if(((MNStruct *)GHSglobalcontext)->numFitEQUAD == 1){Cube[profdims] = -7;}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	long double LDparams[((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps];
	int pcount = profdims+epochdims+epochpriordims;
	int phaseDim = profdims+epochdims+epochpriordims;
	int fitcount = 0;

	for(int p=0;p< ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; p++){
		//printf("Timing: %i %g %g \n", p, Cube[pcount], double(Cube[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[p][0])));
		LDparams[p]=Cube[pcount]*(((MNStruct *)GHSglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GHSglobalcontext)->LDpriors[p][0]);
		pcount++;

	}	
	
	long double phase = LDparams[0]*((MNStruct *)GHSglobalcontext)->ReferencePeriod/SECDAY;
	fitcount++;

	double PhasePrior = ((MNStruct *)GHSglobalcontext)->PhasePrior;
	double phasePriorTerm = 0.5*Cube[phaseDim]*Cube[phaseDim]/(PhasePrior*PhasePrior);
	double phasePriorGrad = Cube[phaseDim]/(PhasePrior*PhasePrior);


	if(((MNStruct *)GHSglobalcontext)->doLinear==0){

		for(int p=1;p<((MNStruct *)GHSglobalcontext)->numFitTiming;p++){
			((MNStruct *)GHSglobalcontext)->pulse->param[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][0]].val[((MNStruct *)GHSglobalcontext)->TempoFitNums[p][1]] = LDparams[fitcount];
			fitcount++;
		}
		for(int p=0;p<((MNStruct *)GHSglobalcontext)->numFitJumps;p++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[((MNStruct *)GHSglobalcontext)->TempoJumpNums[p]]= LDparams[fitcount];
			fitcount++;
		}

		long double *JumpVec = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
		   JumpVec[p] = 0;
		   for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			for (int l=0;l<((MNStruct *)GHSglobalcontext)->pulse->obsn[p].obsNjump;l++){
				if(((MNStruct *)GHSglobalcontext)->pulse->obsn[p].jump[l]==k){
					if(((MNStruct *)GHSglobalcontext)->pulse->fitJump[k] == 0){

						JumpVec[p] += (long double)((MNStruct *)GHSglobalcontext)->PreJumpVals[k]/SECDAY;

					}
					else{

						JumpVec[p] += (long double) ((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]/SECDAY;
					}
				}
			 }
		   }

		}



		for(int p=0;p< ((MNStruct *)GHSglobalcontext)->pulse->nobs; p++){
			((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat = ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase + JumpVec[p];
			//printf("Mod SAT: %.20Lg %.20Lg %.20Lg\n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].origsat-phase,  JumpVec[p], ((MNStruct *)GHSglobalcontext)->pulse->obsn[p].sat);
		}
		for (int k=1;k<=((MNStruct *)GHSglobalcontext)->pulse->nJumps;k++){
			((MNStruct *)GHSglobalcontext)->pulse->jumpVal[k]= 0;
		}

		delete[] JumpVec;

	
		fastformBatsAll(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars);       
		formResiduals(((MNStruct *)GHSglobalcontext)->pulse,((MNStruct *)GHSglobalcontext)->numberpulsars,0);      


		//UpdatePhysDVector(GHSglobalcontext, GlobalDMatrix, ((MNStruct *)GHSglobalcontext)->pulse->nobs);

	}

	
	if(debug == 1)printf("Phase: %g \n", (double)phase);
	if(debug == 1)printf("Formed Residuals \n");



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Epoch Signals///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	int EPriorcount = profdims;
	double *EQUAD;

	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){EQUAD = new double[((MNStruct *)GHSglobalcontext)->numFitEQUAD]; }
	//double *EQUADSignal = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();
	
//	if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
	for(int eq = 0; eq < ((MNStruct *)GHSglobalcontext)->numFitEQUAD; eq++){

		if(debug == 1){printf("EQUAD: %i %g \n",EPriorcount,Cube[EPriorcount]);}
		if(testuniformpriors == 0){
			if(Cube[EPriorcount] < -10){ 
				priorterm -= 2*log(10.0)*(Cube[EPriorcount]+10); 

				for(int i = 0; i < threads; i++){
					PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
				}	
			}
		}
		else{
			priorterm -= 2*log(10.0)*Cube[EPriorcount]; 

			for(int i = 0; i < threads; i++){
				PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
			}
		}

		EQUAD[eq]=pow(10.0,Cube[EPriorcount]);
		EPriorcount++;

		if(((MNStruct *)GHSglobalcontext)->EQUADPriorType ==1) {uniformpriorterm += log(EQUAD[eq]);}


		
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////DM Variations///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	double *SignalVec;
	double *PowerSpec;
	double *FMatrix;

	double FreqLike = 0;
	int startpos = 0;
	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;

		FMatrix = new double[FitDMCoeff*((MNStruct *)GHSglobalcontext)->numProfileEpochs]();
		PowerSpec = new double[FitDMCoeff]();
		SignalVec = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();

		
		double *SignalCoeff = new double[FitDMCoeff];
		for(int i = 0; i < FitDMCoeff; i++){
			SignalCoeff[startpos + i] = Cube[pcount];
			pcount++;
		}
			
		

		double Tspan = ((MNStruct *)GHSglobalcontext)->Tspan;
		double f1yr = 1.0/3.16e7;


		if(((MNStruct *)GHSglobalcontext)->incDM==5){

			double DMamp=Cube[EPriorcount];



			
			DMamp=pow(10.0, DMamp); 
			if(((MNStruct *)GHSglobalcontext)->DMPriorType == 1) { 
				priorterm -= 2*log(DMamp); 
				for(int i = 0; i < threads; i++){
					PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
				}
			}

			EPriorcount++;
			double DMindex=Cube[EPriorcount];
			EPriorcount++;

			for (int i=0; i< FitDMCoeff/2; i++){

				FreqLike += SignalCoeff[i]*SignalCoeff[i] + SignalCoeff[i+FitDMCoeff/2]*SignalCoeff[i+FitDMCoeff/2];
				
				double freq = ((double)(i+1.0))/Tspan;
				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freq*365.25,(-DMindex))/(Tspan*24*60*60);
		
				PowerSpec[startpos + i] = sqrt(rho);
				PowerSpec[startpos + i + FitDMCoeff/2] = sqrt(rho);

				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  
				
			}
		}


		if(((MNStruct *)GHSglobalcontext)->incDM==6){




			for (int i=0; i< FitDMCoeff/2; i++){
			

				double DMAmp = pow(10.0, Cube[EPriorcount]);	
				double freq = ((double)(i+1.0))/Tspan;
				
				if(((MNStruct *)GHSglobalcontext)->DMPriorType == 1) { 
					priorterm -= 2*log(DMAmp); 
					for(int i = 0; i < threads; i++){
						PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
					}
				}

				EPriorcount++;

				FreqLike += SignalCoeff[i]*SignalCoeff[i] + SignalCoeff[i+FitDMCoeff/2]*SignalCoeff[i+FitDMCoeff/2];

				double rho = (DMAmp*DMAmp);
				SignalCoeff[i] = SignalCoeff[i]*sqrt(rho);
				SignalCoeff[i+FitDMCoeff/2] = SignalCoeff[i+FitDMCoeff/2]*sqrt(rho);  




			}
		}

		
		for(int i=0;i< FitDMCoeff/2;i++){
			int DMt = 0;
			for(int k=0;k<((MNStruct *)GHSglobalcontext)->numProfileEpochs;k++){
				double time=(double)((MNStruct *)GHSglobalcontext)->pulse->obsn[DMt].bat;
	
				FMatrix[k + (i+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]=cos(2*M_PI*(double(i+1)/Tspan)*time);
				FMatrix[k + (i+FitDMCoeff/2+startpos)*((MNStruct *)GHSglobalcontext)->numProfileEpochs] = sin(2*M_PI*(double(i+1)/Tspan)*time);
				DMt += ((MNStruct *)GHSglobalcontext)->numChanPerInt[k];

			}
		}

		
		vector_dgemv(FMatrix,SignalCoeff,SignalVec,((MNStruct *)GHSglobalcontext)->numProfileEpochs, FitDMCoeff,'N');
		startpos=FitDMCoeff;
		delete[] SignalCoeff;	
		

    	}


	

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	if(dotime == 3){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to ProfileParams: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);
	}



	int NEpochs = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;
	int TotalProfs = 0;
	for(int ep = 0; ep < NEpochs; ep++){ TotalProfs +=  ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep]; }



	 long double **ProfileBats=new long double*[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)GHSglobalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

	        int ProfNbin  = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[i][1];

		ProfileBats[i] = new long double[2];

		ProfileBats[i][0] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][0][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;
		ProfileBats[i][1] = ((MNStruct *)GHSglobalcontext)->ProfileData[i][ProfNbin-1][0] + ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr;

	      
	      ModelBats[i] = ((MNStruct *)GHSglobalcontext)->ProfileInfo[i][5]+((MNStruct *)GHSglobalcontext)->pulse->obsn[i].batCorr - phase - ((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual/SECDAY;

		//printf("res: %i %g \n", i, (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[i].residual);
	 }


	if(((MNStruct *)GHSglobalcontext)->doLinear>0){

		//printf("in do linear %i %i\n", phaseDim+1, phaseDim+((MNStruct *)GHSglobalcontext)->numFitTiming);

		double *TimingParams = new double[((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps]();
		double *TimingSignal = new double[((MNStruct *)GHSglobalcontext)->pulse->nobs];

		int fitcount = 1;
		for(int i = phaseDim+1; i < phaseDim+((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps; i++){
			//printf("in do linear %i %i \n", i, fitcount);
			///printf("in do linear %i %i %g \n", i, fitcount, Cube[i]);
			TimingParams[fitcount] = Cube[i];
			fitcount++;
		}

		vector_dgemv(GlobalDMatrixVec, TimingParams,TimingSignal,((MNStruct *)GHSglobalcontext)->pulse->nobs,((MNStruct *)GHSglobalcontext)->numFitTiming+((MNStruct *)GHSglobalcontext)->numFitJumps,'N');

		for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->pulse->nobs; i++){

	     		 ModelBats[i] -=  (long double) TimingSignal[i]/SECDAY;

		}	

		delete[] TimingSignal;
		delete[] TimingParams;	

	}
	

	if(dotime == 3){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to sorting bats: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);
	}


	int maxshapecoeff = 0;
	int totshapecoeff = ((MNStruct *)GHSglobalcontext)->totshapecoeff; 

	int *numcoeff= new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		numcoeff[i] =  ((MNStruct *)GHSglobalcontext)->numshapecoeff[i];
		if(debug == 1){printf("num coeff in comp %i: %i \n", i, numcoeff[i]);}
	}

	
        int *numProfileStocCoeff = ((MNStruct *)GHSglobalcontext)->numshapestoccoeff;
        int totalshapestoccoeff = ((MNStruct *)GHSglobalcontext)->totalshapestoccoeff;


	int *numEvoCoeff = ((MNStruct *)GHSglobalcontext)->numEvoCoeff;
	int totalEvoCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoCoeff;

	int *numEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->numEvoFitCoeff;
	int totalEvoFitCoeff = ((MNStruct *)GHSglobalcontext)->totalEvoFitCoeff;

	int *numProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->numProfileFitCoeff;
	int totalProfileFitCoeff = ((MNStruct *)GHSglobalcontext)->totalProfileFitCoeff;



	int totalCoeffForMult = 0;
	int *NumCoeffForMult = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		NumCoeffForMult[i] = numProfileStocCoeff[i];
		if(numEvoCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numEvoCoeff[i];}
		if(numProfileFitCoeff[i] > NumCoeffForMult[i]){NumCoeffForMult[i]=numProfileFitCoeff[i];}
		totalCoeffForMult += NumCoeffForMult[i];
		if(debug == 1){printf("num coeff for mult from comp %i: %i \n", i, NumCoeffForMult[i]);}
	}

        int *numShapeToSave = new int[((MNStruct *)GHSglobalcontext)->numProfComponents];
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                numShapeToSave[i] = numProfileStocCoeff[i];
                if(numEvoCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numEvoCoeff[i];}
                if(numProfileFitCoeff[i] >  numShapeToSave[i]){numShapeToSave[i] = numProfileFitCoeff[i];}
                if(debug == 1){printf("saved %i %i \n", i, numShapeToSave[i]);}
        }
        int totShapeToSave = 0;
        for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
                totShapeToSave += numShapeToSave[i];
        }

	double shapecoeff[totshapecoeff];
	double StocProfCoeffs[totalshapestoccoeff];
	double **EvoCoeffs=new double*[((MNStruct *)GHSglobalcontext)->NProfileEvoPoly]; 
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; i++){EvoCoeffs[i] = new double[totalEvoCoeff];}

	double ProfileFitCoeffs[totalProfileFitCoeff];
	

	for(int i =0; i < totshapecoeff; i++){
		shapecoeff[i]=((MNStruct *)GHSglobalcontext)->MeanProfileShape[i];
		//printf("loaded shape coeff %i %g \n", i, shapecoeff[i]);
	}
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		for(int i =0; i < totalEvoCoeff; i++){
			EvoCoeffs[p][i]=((MNStruct *)GHSglobalcontext)->MeanProfileEvo[p][i];
			//printf("loaded evo coeff %i %g \n", i, EvoCoeffs[i]);
		}
	}
	if(debug == 1){printf("Filled %i Coeff, %i EvoCoeff \n", totshapecoeff, totalEvoCoeff);}
	double *betas = new double[((MNStruct *)GHSglobalcontext)->numProfComponents]();
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		betas[i] = ((MNStruct *)GHSglobalcontext)->MeanProfileBeta[i]*((MNStruct *)GHSglobalcontext)->ReferencePeriod;
	}


	//double *StocProfPriorGrad = new double[totalshapestoccoeff]();
	for(int i =0; i < totalshapestoccoeff; i++){
		//Cube[EPriorcount]=GlobalStartPoint[EPriorcount];
		StocProfCoeffs[i]= pow(10.0, Cube[EPriorcount]);

		if(debug == 1){printf("Stoc: %i %g \n",EPriorcount,Cube[EPriorcount]);}



		if(testuniformpriors == 0){
			if(Cube[EPriorcount] < -5){ 
				priorterm -= 2*log(10.0)*(Cube[EPriorcount]+5); 

				for(int i = 0; i < threads; i++){
					PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;
				}

			}
		}
		else{

			priorterm -= 2*log(10.0)*Cube[EPriorcount]; 

			for(int i = 0; i < threads; i++){
				PriorThreadGrads[EPriorcount-profdims][i] -= log(10.0)/threads;	
			}
		}

		
		EPriorcount++;
	}

	for(int i =0; i < totalProfileFitCoeff; i++){
		
		ProfileFitCoeffs[i]= Cube[pcount];
		pcount++;
	}
	double LinearProfileWidth=0;
	if(((MNStruct *)GHSglobalcontext)->FitLinearProfileWidth == 1){ 
		LinearProfileWidth = Cube[pcount];
		pcount++;
	}


	int  EvoFreqExponent = 1;
	if(((MNStruct *)GHSglobalcontext)->FitEvoExponent == 1){
		EvoFreqExponent =  floor(Cube[pcount]);
		pcount++;
	}
	
	for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
		int cpos = 0;
		for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
			for(int i =0; i < numEvoFitCoeff[c]; i++){
				//printf("Evo: %i %g %g \n", i, EvoCoeffs[p][i+cpos] ,Cube[pcount]);
				EvoCoeffs[p][i+cpos] += Cube[pcount];
				pcount++;
			}
			cpos += numEvoCoeff[c];
		
		}
	}

	double EvoProfileWidth=0;
	if(((MNStruct *)GHSglobalcontext)->incProfileEvo == 2){
		EvoProfileWidth = Cube[pcount];
		pcount++;
	}

	double EvoEnergyProfileWidth=0;
	if(((MNStruct *)GHSglobalcontext)->incProfileEnergyEvo == 2){
		EvoEnergyProfileWidth = Cube[pcount];
		pcount++;
	}

	if(totshapecoeff+1>=totalshapestoccoeff+1){
		maxshapecoeff=totshapecoeff+1;
	}
	if(totalshapestoccoeff+1 > totshapecoeff+1){
		maxshapecoeff=totalshapestoccoeff+1;
	}


	double modelflux=0;
	int fluxpos = 0;
	for(int i = 0; i < ((MNStruct *)GHSglobalcontext)->numProfComponents; i++){
		for(int j =0; j < numcoeff[i]; j=j+2){
			modelflux+=sqrt(sqrt(M_PI))*sqrt(betas[i])*pow(2.0, 0.5*(1.0-j))*sqrt(((MNStruct *)GHSglobalcontext)->Binomial[j])*shapecoeff[fluxpos+j];
		}
		fluxpos+= numcoeff[i];
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double lnew = 0;
	int GlobalNBins = (int)((MNStruct *)GHSglobalcontext)->LargestNBins;
	int ProfileBaselineTerms = ((MNStruct *)GHSglobalcontext)->ProfileBaselineTerms;





	int minep = 0;
	int maxep = ((MNStruct *)GHSglobalcontext)->numProfileEpochs;


	if(dotime == 1){

		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed Up to Start of main loop: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

		//wtime = omp_get_wtime ( ) - wtime;
		//cout << "  Elapsed wall clock time = " << wtime << "\n";

		wtime = omp_get_wtime ( ); 
	}	

	
	double *eplikes = new double[((MNStruct *)GHSglobalcontext)->numProfileEpochs]();


	#pragma omp parallel for 
	for(int ep = minep; ep < maxep; ep++){

		double *ProfileModCoeffs = new double[totalCoeffForMult];

		int thisthread = omp_get_thread_num();
	
		int t = 0;

		for(int sep = 0; sep < ep; sep++){	
			t += ((MNStruct *)GHSglobalcontext)->numChanPerInt[sep];
		}
		


		int NChanInEpoch = ((MNStruct *)GHSglobalcontext)->numChanPerInt[ep];
		int NEpochBins = NChanInEpoch*GlobalNBins;


		double EpochChisq = 0;	
		double EpochdetN = 0;
		double EpochLike = 0;

		double EQUADSignal = 0;

		int EDimcount = profdims + epochpriordims + ep*DimsPerEpoch;
		int OriginalEDimCount = EDimcount;
		if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){

			if(debug == 1){printf("EQUAD Signal: %i %g \n",EDimcount,Cube[EDimcount]);}

			int EQIndex = 0;
                        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[t];}


			EQUADSignal = Cube[EDimcount]*EQUAD[EQIndex];
			eplikes[ep] += -0.5*Cube[EDimcount]*Cube[EDimcount];

			EDimcount++;
		}

		double *OneEpochStochAmps = new double[totalshapestoccoeff];
		for(int i = 0; i < totalshapestoccoeff; i++){

			if(debug == 1){printf("Stoc Signal: %i %g \n",EDimcount,Cube[EDimcount]);}

			OneEpochStochAmps[i] = Cube[EDimcount]*StocProfCoeffs[i];
			eplikes[ep] += -0.5*Cube[EDimcount]*Cube[EDimcount];
                        EDimcount++;	
		}
		
		



		for(int ch = 0; ch < NChanInEpoch; ch++){

			EDimcount = OriginalEDimCount;

			if(dotime == 2){	
				gettimeofday(&tval_before, NULL);
			}
			if(debug == 1){
				printf("In toa %i \n", t);
				printf("sat: %.15Lg \n", ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].sat);
			}


			int nTOA = t;


			double ProfileBaseline = Cube[nTOA*3 + 0];
			double ProfileAmp = Cube[nTOA*3 + 1];
			double ProfileSigma = Cube[nTOA*3 + 2];

			double detN  = 0;
			double Chisq  = 0;
  

			double profilelike=0;

			long double FoldingPeriod = ((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][0];
			long double FoldingPeriodDays = FoldingPeriod/SECDAY;

			int Nbins = GlobalNBins;
			int ProfNbins = (int)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][1];
			int BinRatio = Nbins/ProfNbins;

			double Tobs = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][2];
			double noiseval = (double)((MNStruct *)GHSglobalcontext)->ProfileInfo[nTOA][3];
			long double ReferencePeriod = ((MNStruct *)GHSglobalcontext)->ReferencePeriod;


			double *shapevec  = new double[Nbins];
			double *ProfileModVec = new double[Nbins]();
			double *ProfileJitterModVec = new double[Nbins]();

	

			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get and modify binpos///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////
		    
			long double binpos = ModelBats[nTOA]; 


			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){
				binpos += (long double) EQUADSignal/SECDAY;
			}

			if(((MNStruct *)GHSglobalcontext)->incDM > 4 || ((MNStruct *)GHSglobalcontext)->yearlyDM == 1){
	                	double DMKappa = 2.410*pow(10.0,-16);
        		        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				long double DMshift = (long double)(SignalVec[ep]*DMScale);
				
		 		binpos+=DMshift/SECDAY;

			}


			if(binpos < ProfileBats[nTOA][0])binpos+=FoldingPeriodDays;
			if(binpos > ProfileBats[nTOA][1])binpos-=FoldingPeriodDays;


			if(binpos > ProfileBats[nTOA][1]){printf("OverBoard! %.10Lg %.10Lg %.10Lg\n", binpos, ProfileBats[nTOA][1], (binpos-ProfileBats[nTOA][1])/FoldingPeriodDays);}

			long double minpos = binpos - FoldingPeriodDays/2;
			if(minpos < ProfileBats[nTOA][0])minpos=ProfileBats[nTOA][0];
			long double maxpos = binpos + FoldingPeriodDays/2;
			if(maxpos> ProfileBats[nTOA][1])maxpos =ProfileBats[nTOA][1];

			/////////////////////////////////////////////////////////////////////////////////////////////  
			/////////////////////////Get Interpolation Bin///////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////

			int InterpBin = 0;
			double FirstInterpTimeBin = 0;
			int  NumWholeBinInterpOffset = 0;

			if(((MNStruct *)GHSglobalcontext)->InterpolateProfile == 1){

		
				long double timediff = 0;
				long double bintime = ProfileBats[t][0];


				if(bintime  >= minpos && bintime <= maxpos){
				    timediff = bintime - binpos;
				}
				else if(bintime < minpos){
				    timediff = FoldingPeriodDays+bintime - binpos;
				}
				else if(bintime > maxpos){
				    timediff = bintime - FoldingPeriodDays - binpos;
				}

				timediff=timediff*SECDAY;

				double OneBin = FoldingPeriod/Nbins;
				int NumBinsInTimeDiff = floor(timediff/OneBin + 0.5);
				double WholeBinsInTimeDiff = NumBinsInTimeDiff*FoldingPeriod/Nbins;
				double OneBinTimeDiff = -1*((double)timediff - WholeBinsInTimeDiff);

				double PWrappedTimeDiff = (OneBinTimeDiff - floor(OneBinTimeDiff/OneBin)*OneBin);

				if(debug == 1)printf("Making InterpBin: %g %g %i %g %g %g\n", (double)timediff, OneBin, NumBinsInTimeDiff, WholeBinsInTimeDiff, OneBinTimeDiff, PWrappedTimeDiff);

				InterpBin = floor(PWrappedTimeDiff/((MNStruct *)GHSglobalcontext)->InterpolatedTime+0.5);
				if(InterpBin >= ((MNStruct *)GHSglobalcontext)->NumToInterpolate)InterpBin -= ((MNStruct *)GHSglobalcontext)->NumToInterpolate;

				FirstInterpTimeBin = -1*(InterpBin-1)*((MNStruct *)GHSglobalcontext)->InterpolatedTime;

				if(debug == 1)printf("Interp Time Diffs: %g %g %g %g \n", ((MNStruct *)GHSglobalcontext)->InterpolatedTime, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime, PWrappedTimeDiff, InterpBin*((MNStruct *)GHSglobalcontext)->InterpolatedTime-PWrappedTimeDiff);

				double FirstBinOffset = timediff-FirstInterpTimeBin;
				double dNumWholeBinOffset = FirstBinOffset/(FoldingPeriod/Nbins);
				int  NumWholeBinOffset = 0;

				NumWholeBinInterpOffset = floor(dNumWholeBinOffset+0.5);
	
				if(debug == 1)printf("Interp bin is: %i , First Bin is %g, Offset is %i \n", InterpBin, FirstInterpTimeBin, NumWholeBinInterpOffset);


			}



			//if(nTOA==1)printf("Jump Check: %i %i %g %.20Lg %g %g\n", (double)((MNStruct *)GHSglobalcontext)->pulse->obsn[nTOA].residual, ModelBats[nTOA], InterpBin, NumWholeBinInterpOffset, Cube[profdims+1], GlobalDMatrix[nTOA][1]*Cube[profdims+1]);


	
			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed up to start of Interp: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			double reffreq = ((MNStruct *)GHSglobalcontext)->EvoRefFreq;
			double freqdiff =  (((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freq - reffreq)/1000.0;
			double freqscale = pow(freqdiff, EvoFreqExponent);


			double snr = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].snr;
			double tobs = ((MNStruct *)GHSglobalcontext)->pulse->obsn[t].tobs; 


			snr = snr*3600/tobs;

			double refSN = 1000;
			double SNdiff =  snr/refSN;
			double SNscale = snr-refSN;




			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Add in any Profile Changes///////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(int i =0; i < totalCoeffForMult; i++){
				ProfileModCoeffs[i]=0;	
			}				

			int cpos = 0;
			int epos = 0;
			int fpos = 0;
			int spos = 0;
			for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
				
				for(int i =0; i < numProfileFitCoeff[c]; i++){
					ProfileModCoeffs[i+cpos] += ProfileFitCoeffs[i+fpos];

				}
				for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; p++){	
					for(int i =0; i < numEvoCoeff[c]; i++){
						ProfileModCoeffs[i+cpos] += EvoCoeffs[p][i+epos]*pow(freqscale, p+1);						
					}
				}

				for(int i = 0; i < numProfileStocCoeff[c]; i++){
					ProfileModCoeffs[i+cpos] += OneEpochStochAmps[i+spos]*modelflux;
					//printf("adding stoc: %i %i %i %g \n", ep, ch, i,OneEpochStochAmps[i+spos] );
				}

				cpos += NumCoeffForMult[c];
				epos += numEvoCoeff[c];
				fpos += numProfileFitCoeff[c];	
				spos += numProfileStocCoeff[c];
			}


			

			if(totalCoeffForMult > 0){
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], ProfileModCoeffs,ProfileModVec,Nbins,totalCoeffForMult,'N');
				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfileVec[InterpBin], ProfileModCoeffs,ProfileJitterModVec,Nbins,totalCoeffForMult,'N');  

			}


	

			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Modded Profile: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}





			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Fill Arrays with interpolated state//////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			
			int ZeroWrap = GHSWrap(0 + NumWholeBinInterpOffset, 0, Nbins-1);


			//double corrAmpGrad = 0;
	
			for(int j =0; j < Nbins-ZeroWrap; j++){

				double NewIndex = (j + NumWholeBinInterpOffset);
				int Nj =  ZeroWrap+j;

				

				double widthTerm = ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*LinearProfileWidth; 
				double evoWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoProfileWidth*freqscale;
				double SNWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoEnergyProfileWidth*SNscale;

				shapevec[j] = ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj] + widthTerm + ProfileModVec[Nj] + evoWidthTerm + SNWidthTerm;

			}
	
			for(int j = Nbins-ZeroWrap; j < Nbins; j++){

				double NewIndex = (j + NumWholeBinInterpOffset);
				int Nj =  j-Nbins+ZeroWrap;


				double widthTerm = ((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*LinearProfileWidth; 
				double evoWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoProfileWidth*freqscale;
				double SNWidthTerm = 0;//((MNStruct *)GHSglobalcontext)->InterpolatedWidthProfile[InterpBin][Nj]*EvoEnergyProfileWidth*SNscale;

				shapevec[j] = ((MNStruct *)GHSglobalcontext)->InterpolatedMeanProfile[InterpBin][Nj] + widthTerm + ProfileModVec[Nj] + evoWidthTerm + SNWidthTerm;

			}


			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Filled Arrays: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}


			
			double *datavec = new double[ProfNbins]();
			double *ModelVec = new double[2*ProfNbins]();
			double *MM = new double[2*2]();
			double *Md = new double[2]();
			double *TempMd = new double[2]();

			for(int i =0; i < ProfNbins; i++){
				datavec[i] = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1];
				ModelVec[i + 0*ProfNbins] = 1;
				ModelVec[i + 1*ProfNbins] = shapevec[i*BinRatio];
			}

			vector_dgemm(ModelVec, ModelVec, MM, ProfNbins, 2,ProfNbins, 2, 'T', 'N');
			vector_dgemv(ModelVec, datavec,Md,ProfNbins,2,'T');

			for(int i =0; i < 2; i++){
				TempMd[i] = Md[i];
			}

			int info=0;
			double Margindet = 0;
			vector_dpotrfInfo(MM, 2, Margindet, info);
			vector_dpotrsInfo(MM, TempMd, 2, info);

			double MLSigma = 0;
			for(int i =0; i < ProfNbins; i++){
				double res2 = datavec[i] - TempMd[0] - TempMd[1]*shapevec[i*BinRatio];
				MLSigma += res2*res2;
			}

			MLSigma = sqrt(MLSigma/ProfNbins);
			double noisemean = TempMd[0];
			double MLAmp = TempMd[1];

                        Cube[nTOA*3 + 0]=noisemean;
                        Cube[nTOA*3 + 1]=MLAmp;
                        Cube[nTOA*3 + 2]=MLSigma;


			ProfileBaseline = noisemean;
			ProfileAmp = MLAmp;
			ProfileSigma = MLSigma;

			delete[] datavec;
			delete[] ModelVec;
			delete[] MM;
			delete[] Md;
			delete[] TempMd;


		///////////////////////////////////////////Marginalise over arbitrary offset and absolute amplitude////////////////////////////////////////////////////////////

			      
			Chisq = 0;


			double noise = ProfileSigma*ProfileSigma;
			double OffPulsedetN = log(noise);
			detN = ProfNbins*OffPulsedetN;
			noise=1.0/noise;

			double *NResVec = new double[Nbins]();
		

			ZeroWrap = GHSWrap(0 + (Nbins - NumWholeBinInterpOffset), 0, Nbins-1);

			int ZeroOffset = (Nbins-GHSWrap((Nbins - NumWholeBinInterpOffset), 0, Nbins-1))%BinRatio;


			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Allocated Arrays: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}
			

			int startpoint = 0;
			int stoppoint = Nbins - ZeroWrap;

			//for(int i = 0; i < ProfNbins; i++){
			//	if(ProfNbins < Nbins){printf("Prof bins: %i %i %.10g\n", t, i, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][i][1]);}
			//}

			for(int i = startpoint+ZeroOffset; i < stoppoint; i+=BinRatio){

				
				
				int Nj = ZeroWrap+i;


				double gres = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1] - ProfileAmp*shapevec[Nj] - ProfileBaseline;
				NResVec[i] = gres*noise;
				Chisq +=  gres*NResVec[i];


				grad[nTOA*3 + 0] +=   -NResVec[i]; //Baseline
				grad[nTOA*3 + 1] +=   -shapevec[Nj]*NResVec[i]; //Amp
				grad[nTOA*3 + 2] +=   -gres*NResVec[i]/ProfileSigma + 1.0/ProfileSigma;    //Sigma

				//if(ProfNbins  == 256)printf("Make NRVec %i %i %i %i %i %i %.10g %.10g %.10g \n", startpoint+ZeroOffset, ZeroOffset, stoppoint, t, i, Nj, double(Nj)/BinRatio, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1], ProfileAmp*shapevec[Nj] + ProfileBaseline);


			}

			startpoint = Nbins-ZeroWrap; // + BinRatio - (Nbins-ZeroWrap-1)%BinRatio - 1; 
			stoppoint = Nbins;
	
			for(int i = startpoint; i < stoppoint; i+=BinRatio){


				
			
				int Nj = ZeroWrap+i-Nbins;

				double gres = ((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1] - ProfileAmp*shapevec[Nj] - ProfileBaseline;
	
				NResVec[i] = gres*noise;
				Chisq +=  gres*NResVec[i];


				grad[nTOA*3 + 0] +=   -NResVec[i]; //Baseline
				grad[nTOA*3 + 1] +=   -shapevec[Nj]*NResVec[i]; //Amp
				grad[nTOA*3 + 2] +=   -gres*NResVec[i]/ProfileSigma + 1.0/ProfileSigma;    //Sigma

				//if(ProfNbins == 256)printf("Make NRVec %i %i %i %i %i %i %.10g %.10g %.10g \n", startpoint, ZeroOffset, stoppoint, t, i, Nj, double(Nj)/BinRatio, (double)((MNStruct *)GHSglobalcontext)->ProfileData[nTOA][Nj/BinRatio][1], ProfileAmp*shapevec[Nj] + ProfileBaseline);

			}



			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Get EpochGrads///////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int B2C = 0;
			double PhaseGrad = 0;
			for(int j = 0; j < Nbins; j++){
				ProfileJitterModVec[j] += ((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][j];
//				PhaseGrad += -1*NResVec[j]*((MNStruct *)GHSglobalcontext)->InterpolatedJitterProfile[InterpBin][j]*ProfileAmp;
				PhaseGrad += -1*NResVec[j]*ProfileJitterModVec[j]*ProfileAmp;

			}

			for(int j = 0; j < ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps; j++){	
				double phasefac = 1;
				if(j==0){phasefac *= double(((MNStruct *)GHSglobalcontext)->ReferencePeriod);}		
				threadgrads[B2C][thisthread] += GlobalDMatrix[nTOA][j]*phasefac*PhaseGrad;
				B2C++;
//				if(j==1 &&  GlobalDMatrix[nTOA][j] != 0 )printf("JG: %i %i %g %g %g \n", nTOA, j, PhaseGrad,  GlobalDMatrix[nTOA][j], threadgrads[B2C-1][thisthread]);
			}


			int PriorGradLoopCount = 0;
			if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){


                                int EQIndex = 0;
                                if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 1){EQIndex = ((MNStruct *)GHSglobalcontext)->sysFlags[nTOA];}


                                grad[EDimcount] += 1*EQUADSignal/EQUAD[EQIndex]/NChanInEpoch - PhaseGrad*EQUAD[EQIndex];
                                PriorThreadGrads[PriorGradLoopCount+EQIndex][thisthread] += (-PhaseGrad*EQUADSignal/EQUAD[EQIndex])*log(10.0)*EQUAD[EQIndex];

                                EDimcount++;
                                PriorGradLoopCount+=((MNStruct *)GHSglobalcontext)->numFitEQUAD;

			}
			

			if(((MNStruct *)GHSglobalcontext)->incDM==5){

				int FitDMCoeff = 2*((MNStruct *)GHSglobalcontext)->numFitDMCoeff;



				double DMKappa = 2.410*pow(10.0,-16);
        		        double DMScale = 1.0/(DMKappa*pow((double)((MNStruct *)GHSglobalcontext)->pulse->obsn[t].freqSSB,2));
				double DMshift = (SignalVec[ep]*DMScale);

				
				//Gradient for Amplitude of DM
				PriorThreadGrads[PriorGradLoopCount][thisthread] += -PhaseGrad*DMshift*log(10.0);
				PriorGradLoopCount++;

				
				int DMAmpDim = profdims + epochdims + epochpriordims + ((MNStruct *)GHSglobalcontext)->numFitTiming + ((MNStruct *)GHSglobalcontext)->numFitJumps;

				double DMSpecGradTerm = 0;

				for (int i=0; i< FitDMCoeff/2; i++){

					double KPF1 = DMScale*FMatrix[ep + i*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[i];
					double KPF2 = DMScale*FMatrix[ep + (i+FitDMCoeff/2)*((MNStruct *)GHSglobalcontext)->numProfileEpochs]*PowerSpec[i+FitDMCoeff/2];

					threadgrads[B2C+i][thisthread] +=  Cube[DMAmpDim+i]/TotalProfs - PhaseGrad*KPF1;
					threadgrads[B2C+i+FitDMCoeff/2][thisthread] +=  Cube[DMAmpDim+i+FitDMCoeff/2]/TotalProfs - PhaseGrad*KPF2;

					double freq = ((double)(i+1.0))/((MNStruct *)GHSglobalcontext)->Tspan;

					DMSpecGradTerm += KPF1*(-0.5*log(freq*365.25))*Cube[DMAmpDim+i];
					DMSpecGradTerm += KPF2*(-0.5*log(freq*365.25))*Cube[DMAmpDim+i+FitDMCoeff/2];

				
				}

				//Gradient for Spectral Index of DM
				PriorThreadGrads[PriorGradLoopCount][thisthread] += -PhaseGrad*DMSpecGradTerm;


				B2C += FitDMCoeff;
				PriorGradLoopCount++;
			}


			if(totalCoeffForMult > 0){

				double *GradResVec = new double[totalCoeffForMult];

				vector_dgemv(((MNStruct *)GHSglobalcontext)->InterpolatedShapeletsVec[InterpBin], NResVec,GradResVec,Nbins,totalCoeffForMult,'T');

				for(int p = 0; p < ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly+1; p++){

					double fval = -1*ProfileAmp*pow(freqscale, p);

					int cpos = 0;
					int fpos = 0;
					for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
						for(int i = 0; i < numProfileFitCoeff[c]; i++){
			
							threadgrads[B2C + p*totalProfileFitCoeff + fpos + i][thisthread] += fval*GradResVec[i+cpos];
						}
						fpos += numProfileFitCoeff[c];	
						cpos += NumCoeffForMult[c];
					}
				}

				int cpos = 0;
				int spos = 0;
				for(int c = 0; c < ((MNStruct *)GHSglobalcontext)->numProfComponents; c++){
					for(int i = 0; i < numProfileStocCoeff[c]; i++){
						grad[EDimcount] += OneEpochStochAmps[i + spos]/StocProfCoeffs[i + spos]/NChanInEpoch - modelflux*ProfileAmp*GradResVec[i+cpos]*StocProfCoeffs[i + spos];
						EDimcount++;

						//StocProfPriorGrad[i+spos] += -ProfileAmp*OneEpochStochAmps[i + spos]*GradResVec[i+cpos]*log(10.0);
						PriorThreadGrads[PriorGradLoopCount][thisthread] +=  -modelflux*ProfileAmp*OneEpochStochAmps[i + spos]*GradResVec[i+cpos]*log(10.0);
						PriorGradLoopCount++;
					}
					spos += numProfileStocCoeff[c];
					cpos += NumCoeffForMult[c];
				}

				delete[] GradResVec;
			}



			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  Up to LA: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			//int id;

			  //  id = omp_get_thread_num ( );
			  //  cout << "  This is process " << id << " " << t << " " << detN << " " << Chisq << "\n";
			  			

			profilelike = -0.5*(detN + Chisq);
			//printf("Like: %i %.15g %.15g \n", nTOA, detN, Chisq);
			if(debug == 1)printf("Like: %i %.15g %.15g %.15g \n", nTOA, lnew, detN, Chisq);

			EpochChisq+=Chisq;
			EpochdetN+=detN;
			EpochLike+=profilelike;

			delete[] shapevec;
			delete[] ProfileModVec;
			delete[] ProfileJitterModVec;
			
			delete[] NResVec;

			if(dotime == 2){
				gettimeofday(&tval_after, NULL);
				timersub(&tval_after, &tval_before, &tval_resultone);
				printf("Time elapsed,  End of Epoch: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
				gettimeofday(&tval_before, NULL);
			}

			t++;
		}



///////////////////////////////////////////


		if(((MNStruct *)GHSglobalcontext)->incWideBandNoise == 0){
			eplikes[ep] += EpochLike;
			lnew += EpochLike;
		}
		delete[] ProfileModCoeffs;
		delete[] OneEpochStochAmps;

////////////////////////////////////////////
	
	}




	double omplnew = 0;
	for(int ep = minep; ep < maxep; ep++){
		omplnew += eplikes[ep];
	}

	for(int i = 0; i < epochpriordims; i++){
		double gsum = 0;
                for(int j = 0; j < threads; j++){
                        gsum += PriorThreadGrads[i][j];
                }
		 grad[profdims + i] = gsum;
	}

	for(int i = 0; i < globalparams; i++){
		double gsum = 0;
		for(int j = 0; j < threads; j++){
			gsum += threadgrads[i][j];
		}
		 grad[profdims + epochdims + epochpriordims + i] = gsum;
	}
	grad[profdims + epochdims + epochpriordims]+=phasePriorGrad;
	///printf("phase prior grad %g \n", phasePriorGrad);
	if(dotime == 1){
		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed,  End of Loop: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

		wtime = omp_get_wtime ( ) - wtime;
		cout << "  Elapsed wall clock time = " << wtime << "\n";

		wtime = omp_get_wtime ( ) ;
	}

	
	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->pulse->nobs; j++){
	    delete[] ProfileBats[j];
	}
	delete[] ProfileBats;
	delete[] ModelBats;
	delete[] numcoeff;
	delete[] NumCoeffForMult;
	delete[] numShapeToSave;

	for(int j =0; j< ((MNStruct *)GHSglobalcontext)->NProfileEvoPoly; j++){
	    delete[] EvoCoeffs[j];
	}
	delete[] EvoCoeffs;
	delete[] eplikes;


	for(int i = 0; i < globalparams; i++){
		delete[] threadgrads[i];
	}
	delete[] threadgrads;


	for(int i = 0; i < epochpriordims; i++){
		delete[] PriorThreadGrads[i];
	}
	delete[] PriorThreadGrads;
	


	

	double finallikelihood = -1*omplnew + 0.5*FreqLike + phasePriorTerm + 0.5*priorterm;
	likelihood=finallikelihood;
        printf("in like %g  \n", finallikelihood);	
	if(omplnew-phasePriorTerm-0.5*priorterm > GlobalMaxLike && ((MNStruct *)GHSglobalcontext)->WriteNewML == 1){
		GlobalMaxLike = omplnew-phasePriorTerm-0.5*priorterm;
		printf("Like: %.16g \n", GlobalMaxLike);
		write_newMLLike(*ndim, Cube, GlobalMaxLike);
		if(debug == 1)printf("Written new ML file\n");
	
	}
	//sleep(5);

	if(debug == 0){
		for(int i =profdims; i < profdims+epochpriordims+1; i++){
			printf("Parameters and Grads: %i %g %g %.10g\n", i, Cube[i], grad[i], omplnew-phasePriorTerm-0.5*priorterm );
		}
	}

	for(int i = 0; i < profdims+epochpriordims; i++){
		grad[i]=0;
	}


	for(int i =  profdims + epochpriordims; i < *ndim; i++){
		//printf("grad: %i %g\n", i, grad[i]);
//		grad[i]=0;
	}

	if(debug == 1)printf("End Like: %.10g \n", omplnew-phasePriorTerm-0.5*priorterm );
	
	if(((MNStruct *)GHSglobalcontext)->diagonalGHS == 0){
		rotate2Principal(PrinGrad, grad);
	}
	else{
		for(int i = 0; i < *ndim; i++){
			PrinGrad[i] = grad[i];
		}
	}
	



	if(dotime == 1){
		gettimeofday(&tval_after, NULL);
		timersub(&tval_after, &tval_before, &tval_resultone);
		printf("Time elapsed,  Up to End: %ld.%06ld\n", (long int)tval_resultone.tv_sec, (long int)tval_resultone.tv_usec);
		gettimeofday(&tval_before, NULL);

		wtime = omp_get_wtime ( ) - wtime;
		cout << "  Elapsed wall clock time = " << wtime << "\n";


	}




	delete[] Cube;
	delete[] grad;
	delete[] betas;
        if(((MNStruct *)GHSglobalcontext)->numFitEQUAD > 0){ delete[] EQUAD; }

	if(((MNStruct *)GHSglobalcontext)->incDM > 4){

		delete[] FMatrix;
		delete[] PowerSpec;
		delete[] SignalVec;
	}
	
}

