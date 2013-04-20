#include <stdio.h>
#include <vector>
#include "/usr/include/gsl/gsl_sf_gamma.h"
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "tempo2.h"
#include "TempoNest.h"

extern "C" void WhiteMarginGPUWrapper_(double *Noise, double *Res, double *likeInfo, int N, int G);
extern "C" void vHRedGPUWrapper_(double *CovMatrix, double *Res, double *likeInfo, int N);
extern "C" void vHRedMarginGPUWrapper_(double *CovMatrix, double *Res, double *likeInfo, int N, int G);
extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, int N, int F);
extern "C" void LRedMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int F, int G);

void WhiteMarginGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	clock_t startClock,endClock;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps+1;
	double Fitparams[numfit];
	double *EFAC;
	double EQUAD;
	int pcount=0;
	

	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

	for(int p=0;p < numfit; p++){
	
		Fitparams[p]=Cube[p];
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



	double Chisq=0;
	double *likeInfo=new double[2];

	double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Diffvec=new double[((MNStruct *)context)->pulse->nobs];
	
//	startClock = clock();
	dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Diffvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
	}
	//endClock = clock();

	double *Noise=new double[((MNStruct *)context)->pulse->nobs];

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
	}
//	printf("MVmult: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
	
//	startClock = clock();
	WhiteMarginGPUWrapper_(Noise, Diffvec, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);
//	endClock = clock();
	//printf("GPU: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
	
	double det=likeInfo[0];
	Chisq=likeInfo[1];

	if(isnan(det) || isinf(det) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
	
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->pulse->nobs*log(2*M_PI) + det+Chisq);	


	}

	delete[] EFAC;
	delete[] Noise;
	delete[] Fitvec;
	delete[] Diffvec;


	//printf("%g %g %g \n",lnew,det,Chisq);

}


void vHRedGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps+1;
	double Fitparams[numfit];
	double *EFAC;
	double EQUAD, redamp, redalpha;
	int pcount=0;


	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

	for(int p=0;p < numfit; p++){
		Fitparams[p]=Cube[p];
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

	double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Diffvec=new double[((MNStruct *)context)->pulse->nobs];
	dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Diffvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
	}

	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=redalpha;
	double gwamp=pow(10.0,redamp);
	double gwampsquared=gwamp*gwamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-modelalpha)))/(pow(flo,(modelalpha-1)));

	double timdiff=0;

	double covconst=gsl_sf_gamma(1-modelalpha)*sin(0.5*M_PI*modelalpha);
// 	printf("constants: %g %g \n",gwampsquared,covconst);


	
	double *CovMatrix = new double[((MNStruct *)context)->pulse->nobs*((MNStruct *)context)->pulse->nobs];

	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){

		for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
			timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
			double tau=2.0*M_PI*fabs(timdiff);
			double covsum=0;

			for(int k=0; k <=10; k++){
				covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-modelalpha));

			}

			CovMatrix[o1*((MNStruct *)context)->pulse->nobs + o2]=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);
// 			printf("%i %i %g %g %g\n",o1,o2,CovMatrix[o1][o2],fabs(timdiff),covsum);

			if(o1==o2){
				CovMatrix[o1*((MNStruct *)context)->pulse->nobs + o2] += pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD;
			}

		}
	}
	


	double *likeInfo=new double[2];
	vHRedGPUWrapper_(CovMatrix, Diffvec, likeInfo,((MNStruct *)context)->pulse->nobs);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];


	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(covdet+Chisq);	
	}


	delete[] EFAC;
	delete[] Diffvec;
	delete[] Fitvec;
	delete[] likeInfo;
	delete[] CovMatrix;

		//printf("GPU Like: %g %g %g \n",lnew,Chisq,covdet);
}



void vHRedMarginGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps+1;
	double Fitparams[numfit];
	double *EFAC;
	double EQUAD, redamp, redalpha;
	int pcount=0;


	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

	for(int p=0;p < numfit; p++){
		Fitparams[p]=Cube[p];
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
	
	
	double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Diffvec=new double[((MNStruct *)context)->pulse->nobs];
	dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Diffvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
	}

	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=redalpha;
	double gwamp=pow(10.0,redamp);
	double gwampsquared=gwamp*gwamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-modelalpha)))/(pow(flo,(modelalpha-1)));

	double timdiff=0;

	double covconst=gsl_sf_gamma(1-modelalpha)*sin(0.5*M_PI*modelalpha);
// 	printf("constants: %g %g \n",gwampsquared,covconst);


	
	double *CovMatrix = new double[((MNStruct *)context)->pulse->nobs*((MNStruct *)context)->pulse->nobs];

	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){

		for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
			timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
			double tau=2.0*M_PI*fabs(timdiff);
			double covsum=0;

			for(int k=0; k <=10; k++){
				covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-modelalpha));

			}

			CovMatrix[o1*((MNStruct *)context)->pulse->nobs + o2]=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);
// 			printf("%i %i %g %g %g\n",o1,o2,CovMatrix[o1][o2],fabs(timdiff),covsum);

			if(o1==o2){
				CovMatrix[o1*((MNStruct *)context)->pulse->nobs + o2] += pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD;
			}

		}
	}
	
	double *likeInfo=new double[2];
	vHRedMarginGPUWrapper_(CovMatrix, Diffvec, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];

	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}
	else{
		lnew = -0.5*(covdet+Chisq);	


	}
// 	endClock = clock();
// //   	printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
	

	delete[] EFAC;
	delete[] CovMatrix;
	delete[] Diffvec;
	delete[] Fitvec;
	delete[] likeInfo;

}


void LRedGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps+1;
	double Fitparams[numfit];
	double *EFAC;
	double EQUAD;
	int pcount=0;
	
// 	printf("here1\n");
	for(int p=0;p<ndim;p++){
// 		printf("param %i %g %g\n",p,((MNStruct *)context)->Dpriors[p][0],((MNStruct *)context)->Dpriors[p][1]);
		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}
// 	printf("here1.5\n");
	for(int p=0;p < numfit; p++){
		Fitparams[p]=Cube[p];
		pcount++;
// 		printf("param: %i %g \n",p,Fitparams[p]);
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

  	

	double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Diffvec=new double[((MNStruct *)context)->pulse->nobs];
	dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Diffvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
	}

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


	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Diffvec[o],2)*WorkNoise[o];
	
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
	

	LRedGPUWrapper_(freqs, Diffvec, BATvec, WorkNoise, FNF, NFd, ((MNStruct *)context)->pulse->nobs, FitCoeff);


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
	delete[] Diffvec;
	delete[] BATvec;
	delete[] NFd;
	delete[] freqs;
	delete[] Fitvec;

	for (int j = 0; j < FitCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;



	for (int j = 0; j < FitCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;


}


void LRedMarginGPULinearLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps+1;
	double Fitparams[numfit];
	double *EFAC;
	double EQUAD;
	int pcount=0;
	
// 	printf("here1\n");
	for(int p=0;p<ndim;p++){
// 		printf("param %i %g %g\n",p,((MNStruct *)context)->Dpriors[p][0],((MNStruct *)context)->Dpriors[p][1]);
		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}
// 	printf("here1.5\n");
	for(int p=0;p < numfit; p++){
		Fitparams[p]=Cube[p];
		pcount++;
// 		printf("param: %i %g \n",p,Fitparams[p]);
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

  	

	double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Diffvec=new double[((MNStruct *)context)->pulse->nobs];
	dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Diffvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
	}

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


	int coeffsize=FitCoeff/2;
	double *freqs = new double[(FitCoeff/2)];
	for(int i=0;i<FitCoeff/2;i++){
		freqs[i]=double(i+1)/maxtspan;
	}	
	
	double *likeVals=new double[2];
	LRedMarginGPUWrapper_(freqs, Diffvec, BATvec, Noise, FNF, NFd, likeVals, ((MNStruct *)context)->pulse->nobs, FitCoeff,((MNStruct *)context)->Gsize);
	
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
	delete[] Diffvec;
	delete[] Fitvec;


}
