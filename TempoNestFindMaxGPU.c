#include <stdio.h>
#include <vector>
#include "/usr/include/gsl/gsl_sf_gamma.h"
#include <gsl/gsl_multimin.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "tempo2.h"
#include "TempoNest.h"

double *TNMaxFactorialList = new double[21];


extern "C" void WhiteMarginGPUWrapper_(double *Noise, double *Res, double *likeInfo, int N, int G);

extern "C" void vHRedGPUWrapper_(double *SpecInfo, double *BatVec,  double *DMVec, double *Res, double *NoiseVec, double *likeInfo, int N);
extern "C" void vHRedMarginGPUWrapper_(double *Res, double *BatVec, double *DMVec, double *NoiseVec, double *SpecInfo, double *likeInfo, double *FactorialList, int N, int G);

extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, int N, int F);
extern "C" void LRedMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int F, int G);

extern "C" void vHRedDMMarginGPUWrapper_(double *Res, double *BatVec, double *NoiseVec, double *DMVec, double *SpecInfo, double *likeInfo, double *FactorialList, int N, int G);
extern "C" void vHRedDMGPUWrapper_(double *SpecInfo, double *BatVec, double *Res, double *NoiseVec, double *DMVec, double *likeInfo, int N);
extern "C" void LRedDMMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int F, int G);



void WhiteLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;
	
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


	
	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}

	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	}


	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		pcount++;
	}







	double Chisq=0;
	double noiseval=0;
	double detN=0;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		noiseval=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		Chisq += pow(Resvec[o],2)/noiseval;
		detN += log(noiseval);
	}

	if(isnan(detN) || isinf(detN) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->pulse->nobs*log(2*M_PI) + detN + Chisq);	
	}

	delete[] EFAC;
	delete[] Resvec;


}



void WhiteMarginGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}

	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	}
	

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
 		
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		pcount++;

	}


	double Chisq=0;
	double *Noise=new double[((MNStruct *)context)->pulse->nobs];

	double *likeInfo=new double[2];



	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;

	}
	
	WhiteMarginGPUWrapper_(Noise, Resvec, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

	double det=likeInfo[0];
	Chisq=likeInfo[1];

	if(isnan(det) || isinf(det) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
	
	}
	else{
		lnew = -0.5*(det+Chisq);	


	}

	delete[] EFAC;
	delete[] Noise;
	delete[] Resvec;
	delete[] likeInfo;
	
	//printf("Like, %g %g \n", det, Chisq);

	

}




void vHRedMarginGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD, redamp, redalpha, DMamp, DMalpha;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}
	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	}


	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
// 		
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;

	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		pcount++;


	}

	double *SpecParm=new double[6];
	for(int o=0;o<6; o++){
		SpecParm[o]=0;
	}
	
	if(((MNStruct *)context)->incRED==1){
		redamp=Cube[pcount];
		pcount++;
		redalpha=Cube[pcount];
		pcount++;
		SpecParm[0]=redamp;
		SpecParm[1]=redalpha;
	}
	
  	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);
	if(((MNStruct *)context)->incDM==1){
		DMamp=Cube[pcount];
		pcount++;
		DMalpha=Cube[pcount];
		pcount++;
		SpecParm[3]=DMamp;
		SpecParm[4]=DMalpha;
		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
		}
	}

	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noisevec=new double[((MNStruct *)context)->pulse->nobs];
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		Noisevec[o]=(double)(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD);
	}
	 
	



	
	double *likeInfo=new double[2];
	//printf("Entering %g %g \n", SpecParm[0], SpecParm[1]);
	vHRedMarginGPUWrapper_(Resvec, BATvec, DMVec, Noisevec, SpecParm, likeInfo, TNMaxFactorialList, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];

	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(covdet+Chisq);	


	}
// 	endClock = clock();
// //   	printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(double)CLOCKS_PER_SEC);
	

	delete[] EFAC;
	delete[] BATvec;
	delete[] Noisevec;
	delete[] SpecParm;
	delete[] Resvec;
	delete[] DMVec;

}

void vHRedGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD, redamp, redalpha, DMamp, DMalpha;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}

	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		pcount++;
	}


	double *SpecParm=new double[6];
	for(int o=0;o<6; o++){
		SpecParm[o]=0;
	}
	
	if(((MNStruct *)context)->incRED==1){
		redamp=Cube[pcount];
		pcount++;
		redalpha=Cube[pcount];
		pcount++;
		SpecParm[0]=redamp;
		SpecParm[1]=redalpha;
	}
	
  	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);
	if(((MNStruct *)context)->incDM==1){
		DMamp=Cube[pcount];
		pcount++;
		DMalpha=Cube[pcount];
		pcount++;
		SpecParm[3]=DMamp;
		SpecParm[4]=DMalpha;
		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
		}
	}
	

	double *BatVec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		BatVec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
	}


	double *likeInfo=new double[2];


	vHRedGPUWrapper_(SpecParm, BatVec,DMVec, Resvec, Noise, likeInfo,((MNStruct *)context)->pulse->nobs);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];


	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(covdet+Chisq);	
	}


	delete[] EFAC;
	delete[] Resvec;
	delete[] BatVec;
	delete[] Noise;
	delete[] SpecParm;
	delete[] DMVec;

		//printf("GPU Like: %g %g %g \n",lnew,Chisq,covdet);
}



void LRedGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}

	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	}


	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
// 		
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
// 		printf("EQUAD: %g \n",EQUAD);
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		pcount++;
// 		printf("EQUAD: %g %g %g %i \n",EQUAD,EQUADPrior[0],EQUADPrior[1],((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC);
	}

  	


	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[FitCoeff];

	double tdet=0;
	double timelike=0;



	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow((((MNStruct *)context)->pulse->obsn[o].residual),2)*WorkNoise[o];

			BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;


	}

	double *NFd = new double[FitCoeff];


	double **FNF=new double*[FitCoeff];
	for(int i=0;i<FitCoeff;i++){
		FNF[i]=new double[FitCoeff];
	}

	double start,end;
	int go=0;
	for (int i=0;i<((MNStruct *)context)->pulse->nobs;i++)
	  {
	    if (((MNStruct *)context)->pulse->obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    if (end < (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      end = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		  }
	      }
	  }

	double maxtspan=end-start;


	double freqdet=0;
	for (int i=0; i<FitCoeff/2; i++){
		int pnum=pcount;
		double pc=Cube[pcount];
		
		powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
		powercoeff[i+FitCoeff/2]=powercoeff[i];
		freqdet=freqdet+2*log(powercoeff[i]);
		pcount++;
	}


	int coeffsize=FitCoeff/2;
	double *freqs = new double[FitCoeff/2];
	for(int i=0;i<FitCoeff/2;i++){
		freqs[i]=double(i+1)/maxtspan;
	}	
	

	LRedGPUWrapper_(freqs, Resvec, BATvec, WorkNoise, FNF, NFd, ((MNStruct *)context)->pulse->nobs, FitCoeff);


	double **PPFM=new double*[FitCoeff];
	for(int i=0;i<FitCoeff;i++){
		PPFM[i]=new double[FitCoeff];
		for(int j=0;j<FitCoeff;j++){
			PPFM[i][j]=0;

		}
	}


	for(int c1=0; c1<FitCoeff; c1++){
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}



	for(int j=0;j<FitCoeff;j++){
		for(int k=0;k<FitCoeff;k++){

			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}

 	
	double jointdet=0;
	dpotrf(PPFM, FitCoeff, jointdet);
    dpotri(PPFM,FitCoeff);

	double freqlike=0;
	for(int i=0;i<FitCoeff;i++){
		for(int j=0;j<FitCoeff;j++){
			freqlike=freqlike+NFd[i]*PPFM[i][j]*NFd[j];
		}
	}
	
	lnew=-0.5*(jointdet+tdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);

		
	}
 	//printf("Like: %g %g %g %g %g %g\n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

	delete[] EFAC;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] BATvec;
	delete[] NFd;
	delete[] freqs;

	for (int j = 0; j < FitCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;



	for (int j = 0; j < FitCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;



}


void LRedMarginGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}

	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
// 		
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
// 		printf("EQUAD: %g \n",EQUAD);
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		pcount++;
// 		printf("EQUAD: %g %g %g %i \n",EQUAD,EQUADPrior[0],EQUADPrior[1],((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC);
	}



	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	double *powercoeff=new double[FitCoeff];



	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	}

	double *NFd = new double[FitCoeff];
	double **FNF=new double*[FitCoeff];
	for(int i=0;i<FitCoeff;i++){
		FNF[i]=new double[FitCoeff];
	}	
	
	double start,end;
	int go=0;
	for (int i=0;i<((MNStruct *)context)->pulse->nobs;i++)
	  {
	    if (((MNStruct *)context)->pulse->obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      start = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		    if (end < (double)((MNStruct *)context)->pulse->obsn[i].bat)
		      end = (double)((MNStruct *)context)->pulse->obsn[i].bat;
		  }
	      }
	  }

	double maxtspan=end-start;

	double freqdet=0;
	for (int i=0; i<FitCoeff/2; i++){
		int pnum=pcount;
		double pc=Cube[pcount];
		
		powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
		powercoeff[i+FitCoeff/2]=powercoeff[i];
		freqdet=freqdet+2*log(powercoeff[i]);
		pcount++;
	}


	int coeffsize=FitCoeff/2;
	double *freqs = new double[(FitCoeff/2)];
	for(int i=0;i<FitCoeff/2;i++){
		freqs[i]=double(i+1)/maxtspan;
	}	
	
	double *likeVals=new double[2];
	LRedMarginGPUWrapper_(freqs, Resvec, BATvec, Noise, FNF, NFd, likeVals, ((MNStruct *)context)->pulse->nobs, FitCoeff,((MNStruct *)context)->Gsize);
	
	double tdet=likeVals[0];
	double timelike=likeVals[1];
	
	



	double **PPFM=new double*[FitCoeff];
	for(int i=0;i<FitCoeff;i++){
		PPFM[i]=new double[FitCoeff];
		for(int j=0;j<FitCoeff;j++){
			PPFM[i][j]=0;
		}
	}


	for(int c1=0; c1<FitCoeff; c1++){
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}



	for(int j=0;j<FitCoeff;j++){
		for(int k=0;k<FitCoeff;k++){
			
			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}

 	
	double jointdet=0;
	dpotrf(PPFM, FitCoeff, jointdet);
    dpotri(PPFM,FitCoeff);

	double freqlike=0;
	for(int i=0;i<FitCoeff;i++){
		for(int j=0;j<FitCoeff;j++){
// 			printf("%i %i %g %g\n",i,j,NFd[i],PPFM[i][j]);
			freqlike=freqlike+NFd[i]*PPFM[i][j]*NFd[j];
		}
	}
	
	lnew=-0.5*(tdet+jointdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
	
	}


	delete[] EFAC;
	delete[] powercoeff;
	delete[] NFd;

	for (int j = 0; j < FitCoeff; j++){
		delete[]PPFM[j];
	}
	delete[]PPFM;


	for (int j = 0; j < FitCoeff; j++){
		delete[]FNF[j];
	}
	delete[]FNF;

	delete[] Noise;
	delete[] Resvec;
	
	
	//printf("GPULike: %g %g %g %g %g %g\n",lnew,jointdet,tdet,freqdet,timelike,freqlike);


}

/* The negative likelihood function in GSL understandable language
 *
 * Input
 * pvP:	 	Parameters of the likelihood function
 * params:	Extra information. Here a pointer to the model object
 *
 * Output
 *
 * Return:	Log-Likelihood value
 * */
double negloglikelihood(const gsl_vector *pvP, void *context) {

	int ndims=((MNStruct *)context)->numdims;
	double pdParameters[ndims];
	int nParameters;
	double lval;

	
	
	// Obtain the pointer to the model
	
	for(int i=0; i<ndims; i++) {
		pdParameters[i] = gsl_vector_get(pvP, i);
  	} // for i

	if(((MNStruct *)context)->JumpMargin != 0 || ((MNStruct *)context)->TimeMargin != 0 ){
		//printf("Margin \n");
		if(((MNStruct *)context)->incRED==0){
			WhiteLogLikeMax(pdParameters, ndims, ndims, lval, context);
		}
		else if(((MNStruct *)context)->incRED==1){
			vHRedGPULogLikeMax(pdParameters, ndims, ndims, lval, context);
		}
		else if(((MNStruct *)context)->incRED==2){
			LRedGPULogLikeMax(pdParameters, ndims, ndims, lval, context);
		}
	}

	else if(((MNStruct *)context)->JumpMargin == 0 && ((MNStruct *)context)->TimeMargin == 0 ){
		//printf("No Margin \n");
		if(((MNStruct *)context)->incRED==0){
			WhiteLogLikeMax(pdParameters, ndims, ndims, lval, context);
		}
		else if(((MNStruct *)context)->incRED==1){
			vHRedGPULogLikeMax(pdParameters, ndims, ndims, lval, context);
		}
		else if(((MNStruct *)context)->incRED==2){
			LRedGPULogLikeMax(pdParameters, ndims, ndims, lval, context);
		}
	}
	
	//printf("Like: %g \n",lval);
	// The loglikelihood function is virtual, so the correct one is called
	return -lval;
} // loglikelihood

/* This function finds the maximum likelihood parameters for the stochastic
 * signal (here power-law power spectral density with a white high-frequency
 * tail)
 *
 * Input
 * nParameters:		The number of stochastic signal parameters
 * 
 * Output:
 * pdParameters:	The maximum likelihood value of the parameters
 * */
void NelderMeadOptimum(int nParameters, long double *LdParameters, void *context) {

	printf("\n\n Performing Minimisation over current model parameters \n\n");

  const gsl_multimin_fminimizer_type *pmtMinimiserType = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *pmMinimiser = NULL;
  gsl_vector *vStepSize, *vStart;
  gsl_multimin_function MinexFunc;

  double *pdParameterEstimates;
  double dSize;
  int nIteration=0, nStatus;

  pdParameterEstimates = new double[nParameters];
  vStepSize = gsl_vector_alloc(nParameters);
  vStart = gsl_vector_alloc(nParameters);

	int pcount=0;
  // Obtain the starting point as rough parameter estimates
  for(int i=0; i<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++) {
    gsl_vector_set(vStepSize, pcount, 10.0);
    gsl_vector_set(vStart, pcount, 0.0);
	pcount++;
  } // for i

  for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){
		//printf("Setting EFAC: %i\n",i);
	    gsl_vector_set(vStepSize, pcount, 1.0);
   	    gsl_vector_set(vStart, pcount, 2.0);
		pcount++;
	}
  for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
		//printf("Setting EQUD: %i\n",i);
	    gsl_vector_set(vStepSize, pcount, 0.5);
   	    gsl_vector_set(vStart, pcount, -5.0);
		pcount++;
	}
	if(((MNStruct *)context)->incRED == 1){
		//printf("Setting Red\n");
		gsl_vector_set(vStart, pcount, -12.0);
		gsl_vector_set(vStepSize, pcount, 0.5);
		pcount++;
		gsl_vector_set(vStart, pcount, 3.1);
		gsl_vector_set(vStepSize, pcount, 0.25);
		pcount++;
	}

	if(((MNStruct *)context)->incRED == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
			gsl_vector_set(vStepSize, pcount, 0.5);
			gsl_vector_set(vStart, pcount, -3.0);
			pcount++;
		}
	}	

  // Initial GSL vector of the vertex sizes, and the starting point



  // Initialise the iteration procedure
  MinexFunc.f = &negloglikelihood;
  MinexFunc.n = nParameters;
  MinexFunc.params = context;
  pmMinimiser = gsl_multimin_fminimizer_alloc(pmtMinimiserType, nParameters);
  gsl_multimin_fminimizer_set(pmMinimiser, &MinexFunc, vStart, vStepSize);

  // Iterate to the maximum likelihood
  do {
    nIteration++;
    nStatus = gsl_multimin_fminimizer_iterate(pmMinimiser);

    if(nStatus) break;

    // Check whether we are close enough to the minimum (1e-3 error)
    dSize = gsl_multimin_fminimizer_size(pmMinimiser);
    nStatus = gsl_multimin_test_size(dSize, 1e-3);

      for(int i=0; i<nParameters; i++) {
	pdParameterEstimates[i] = gsl_vector_get(pmMinimiser->x, i);
	//printf("%i %g \n", i,pdParameters[i]);
      } // for i

    if(nStatus == GSL_SUCCESS) {
      fprintf(stderr, "Converged to maximum likelihood with downhill simplex                               \n");
    } // if nStatus

    // Print iteration values
	if(nIteration % 100 ==0){
    		printf("Step[%i]: Convergence: %g Current Minimum: %g \n", nIteration,dSize,gsl_multimin_fminimizer_minimum(pmMinimiser));
	}
  } while (nStatus == GSL_CONTINUE);


	printf("\n");

	pcount=0;
	for(int j=0;j<((MNStruct *)context)->numFitTiming;j++){
	
		long double value=pdParameterEstimates[j]*(((MNStruct *)context)->LDpriors[j][1])+(((MNStruct *)context)->LDpriors[j][0]);
		LdParameters[pcount]=value;
		printf("   Max %s : %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], value);
		pcount++;
	}

	for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){

		long double value=pdParameterEstimates[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
		LdParameters[pcount]=value;
		printf("   Max Jump %i :  %.10Lg \n", j,value);
		pcount++;
	}
	
	  for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){
		printf("   Max EFAC %i :  %.10g \n", i,pdParameterEstimates[pcount]);
		pcount++;
	}
  for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
		printf("   Max EQUAD %i :  %.10g \n", i,pdParameterEstimates[pcount]);
		pcount++;
	}
	if(((MNStruct *)context)->incRED == 1){
		printf("   Max Red Amp :  %.10g \n", pdParameterEstimates[pcount]);
		pcount++;
		printf("   Max Red Index :  %.10g \n",pdParameterEstimates[pcount]);
		pcount++;
	}

	if(((MNStruct *)context)->incRED == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
			printf("   Max Red Coeff %i :  %.10g \n", i,pdParameterEstimates[pcount]);
			pcount++;
		}
	}	



  gsl_vector_free(vStart);
  gsl_vector_free(vStepSize);
  gsl_multimin_fminimizer_free(pmMinimiser);
  delete[] pdParameterEstimates;
} // NelderMeadOptimum

