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
extern "C" void vHRedGPUWrapper_(double *SpecInfo, double *BatVec, double *Res, double *NoiseVec, double *likeInfo, int N);
extern "C" void vHRedMarginGPUWrapper_(double *CovMatrix, double *Res, double *likeInfo, int N, int G);
extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, int N, int F);
extern "C" void LRedMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int F, int G);

extern "C" void vHRedGPUWrapper2_(double *Res, double *BATvec, double *Noisevec, double *SpecParm, double *likeInfo, int N);
extern "C" void vHRedMarginGPUWrapper2_(double *Res, double *BATvec, double *Noisevec, double *SpecParm, double *likeInfo,double *FactorialList, int N, int G);



void WhiteLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD;
	int pcount=0;
// 	printf("here1\n");

// 	printf("here1.5\n");
	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
		//printf("CubeTM: %i %g \n", p, Cube[p]); 
	}
// 	printf("here2\n");

	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
// 	printf("here3\n");
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
			//printf("CubeEFAC: %i %g \n", pcount, Cube[pcount]); 
			pcount++;
			
		}
	}				
// 	printf("here4\n");
	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
// 		printf("E: %g %g\n",EFAC[0],EQUAD);
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		//printf("CubeEQUAD: %i %g \n", pcount, Cube[pcount]);
		pcount++;
		 

	}
// 	printf("here\n",EFAC[0],EQUAD);
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */




	double Chisq=0;
	double noiseval=0;
	double detN=0;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		noiseval=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		Chisq += pow((((MNStruct *)context)->pulse->obsn[o].residual),2)/noiseval;
		detN += log(noiseval);
// 		printf("detn: %g %g \n",noiseval,detN);
	}

	if(isnan(detN) || isinf(detN) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->pulse->nobs*log(2*M_PI) + detN + Chisq);	
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);

	}

	delete[] EFAC;


}



void WhiteMarginGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD;
	int pcount=0;


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
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

	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */




	double Chisq=0;
	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	double *likeInfo=new double[2];



	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		Res[o]=((MNStruct *)context)->pulse->obsn[o].residual;
	}
	
	WhiteMarginGPUWrapper_(Noise, Res, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

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
	delete[] Res;
	delete[] likeInfo;
	
	printf("Like, %g %g \n", det, Chisq);


	

}




void vHRedMarginGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD, redamp, redalpha;
	int pcount=0;



	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
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
//		printf("E: %g %g \n",EQUAD,EFAC[0]);

	}


	redamp=Cube[pcount];
	pcount++;
	redalpha=Cube[pcount];
	pcount++;
  	

// 	startClock = clock();
	
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noisevec=new double[((MNStruct *)context)->pulse->nobs];
	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Res[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual;
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		Noisevec[o]=(double)(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD);
	}
	 
	double *SpecParm=new double[3];

	SpecParm[0]=redamp;
	SpecParm[1]=redalpha;

	for(int k=0; k <=20; k++){
		TNMaxFactorialList[k]=iter_factorial(k);
	}
	
	double *likeInfo=new double[2];
	//printf("Entering %g %g \n", SpecParm[0], SpecParm[1]);
	vHRedMarginGPUWrapper2_(Res, BATvec, Noisevec, SpecParm, likeInfo, TNMaxFactorialList, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];

	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
 		//printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
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
	delete[] Res;
	//printf("Like: %g %g %g \n",lnew,Chisq,covdet);

}

void vHRedGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD, redamp, redalpha;
	int pcount=0;


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
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


	redamp=Cube[pcount];
	pcount++;
	redalpha=Cube[pcount];
	pcount++;

	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	
	double *BatVec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Res[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual;
		BatVec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
	}


	double *likeInfo=new double[2];
	double *SpecInfo=new double[2];
	
	SpecInfo[0]=redamp;
	SpecInfo[1]=redalpha;
	vHRedGPUWrapper_(SpecInfo, BatVec,Res, Noise, likeInfo,((MNStruct *)context)->pulse->nobs);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];


	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(covdet+Chisq);	
	}


	delete[] EFAC;
	delete[] Res;
	delete[] BatVec;
	delete[] Noise;
	delete[] SpecInfo;

		//printf("GPU Like: %g %g %g \n",lnew,Chisq,covdet);
}



void LRedGPULogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD;
	int pcount=0;



	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
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

  	

// 	startClock = clock();
	
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */

	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[FitCoeff];

	double freqdet=0;
	double tdet=0;
	double timelike=0;
	for (int i=0; i<FitCoeff/2; i++){
		int pnum=pcount;
		double pc=Cube[pcount];
		
		powercoeff[i]=pow(10.0,pc)/(365.25*24*60*60)/4;
		powercoeff[i+FitCoeff/2]=powercoeff[i];
		freqdet=freqdet+2*log(powercoeff[i]);
		pcount++;
		//prior=prior+pc;
	}

	double *resvec=new double[((MNStruct *)context)->pulse->nobs];
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow((((MNStruct *)context)->pulse->obsn[o].residual),2)*WorkNoise[o];
			resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual;
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


	int coeffsize=FitCoeff/2;
	double *freqs = new double[FitCoeff/2];
	for(int i=0;i<FitCoeff/2;i++){
		freqs[i]=double(i+1)/maxtspan;
	}	
	

	LRedGPUWrapper_(freqs, resvec, BATvec, WorkNoise, FNF, NFd, ((MNStruct *)context)->pulse->nobs, FitCoeff);


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
	delete[] resvec;
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
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD;
	int pcount=0;



	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
	}


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
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

  	

// 	startClock = clock();
	
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */

	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	double *powercoeff=new double[FitCoeff];

	double freqdet=0;
	for (int i=0; i<FitCoeff/2; i++){
		int pnum=pcount;
		double pc=Cube[pcount];
		
		powercoeff[i]=pow(10.0,pc)/(365.25*24*60*60)/4;
		powercoeff[i+FitCoeff/2]=powercoeff[i];
		freqdet=freqdet+2*log(powercoeff[i]);
		pcount++;
	}

	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];


	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		Res[o]=((MNStruct *)context)->pulse->obsn[o].residual;
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


	int coeffsize=FitCoeff/2;
	double *freqs = new double[(FitCoeff/2)];
	for(int i=0;i<FitCoeff/2;i++){
		freqs[i]=double(i+1)/maxtspan;
	}	
	
	double *likeVals=new double[2];
	LRedMarginGPUWrapper_(freqs, Res, BATvec, Noise, FNF, NFd, likeVals, ((MNStruct *)context)->pulse->nobs, FitCoeff,((MNStruct *)context)->Gsize);
	
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
	delete[] Res;
	
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
void NelderMeadOptimum(int nParameters, double *pdParameters, void *context) {

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
	pdParameters[i] = gsl_vector_get(pmMinimiser->x, i);
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
		

	for(int i=0; i<((MNStruct *)context)->numdims; i++) {
		printf("%i %g \n", i,pdParameters[i]);
		pcount++;
	} 
	printf("\n");

	pcount=0;
	for(int j=0;j<((MNStruct *)context)->numFitTiming;j++){
	
		long double value=pdParameters[j]*(((MNStruct *)context)->LDpriors[j][1])+(((MNStruct *)context)->LDpriors[j][0]);
		pdParameters[pcount]=value;
		printf("Max Param: %i %.10Lg \n", pcount,value);
		pcount++;
	}

	for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){

		long double value=pdParameters[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
		pdParameters[pcount]=value;
		printf("Param:  %.10Lg \n", pcount,value);
		pcount++;
	}


  gsl_vector_free(vStart);
  gsl_vector_free(vStepSize);
  gsl_multimin_fminimizer_free(pmMinimiser);
  delete[] pdParameterEstimates;
} // NelderMeadOptimum


