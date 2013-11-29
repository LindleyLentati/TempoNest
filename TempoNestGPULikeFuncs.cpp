#include <stdio.h>
#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "dpotrs.h"
#include "tempo2.h"
#include "TempoNest.h"

double *TNFactorialList = new double[21];



extern "C" void WhiteMarginGPUWrapper_(double *Noise, double *Res, double *likeInfo, int N, int G, int incEFAC, int incEQUAD);

extern "C" void vHRedGPUWrapper_(double *SpecInfo, double *BatVec,  double *DMVec, double *Res, double *NoiseVec, double *likeInfo, int N);
extern "C" void vHRedMarginGPUWrapper_(double *Res, double *BatVec, double *DMVec, double *NoiseVec, double *SpecInfo, double *likeInfo, double *FactorialList, int N, int G);

extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise, double **FNF, double *NFd, int N, int RF, int DMF, int F, int incRED, int incDM);
extern "C" void LRedMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int RF, int DMF, int F, int G, int incRED, int incDM, int incEFAC, int incEQUAD);

extern "C" void vHRedDMMarginGPUWrapper_(double *Res, double *BatVec, double *NoiseVec, double *DMVec, double *SpecInfo, double *likeInfo, double *FactorialList, int N, int G);
extern "C" void vHRedDMGPUWrapper_(double *SpecInfo, double *BatVec, double *Res, double *NoiseVec, double *DMVec, double *likeInfo, int N);
extern "C" void LRedDMMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int F, int G);


void store_factorial(){

	for(int k=0; k <=20; k++){
		TNFactorialList[k]=iter_factorial(k);
	}
}

void WhiteMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;

	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;
	for(int p=0;p<totdims;p++){
		//printf("DP: %i %g %g \n", p, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1]);
		if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
		Cube[pcount]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[pcount]+((MNStruct *)context)->Dpriors[p][0];
		//printf("Cube: %i %g %g %g \n", pcount, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1], Cube[pcount]);
		pcount++;
		}
	}
	pcount=0;
	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){

				LDparams[p]=Cube[fitcount]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				LDparams[p]=((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
			}
			//printf("LD: %i %Lg \n", p, LDparams[p]);

		}
		pcount=0;
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
			//printf("Res %i %g\n", o, Resvec[o]);
		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
		fitcount=0;

		
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=Cube[fitcount];
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=0;
			}
		}
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time = ((MNStruct *)context)->pulse->obsn[o1].bat;
				if(time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
        }

	double Chisq=0;
	double det = 0;
	double *Noise;

	double *likeInfo=new double[2];

	if(((MNStruct *)context)->numFitEFAC >1 || ((MNStruct *)context)->numFitEQUAD >1){
		Noise=new double[((MNStruct *)context)->pulse->nobs];
		if(((MNStruct *)context)->whitemodel == 0){
			for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
				Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
	
			}
		}
		else if(((MNStruct *)context)->whitemodel == 1){
			for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
				Noise[o]=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
	
			}
		}
	
		WhiteMarginGPUWrapper_(Noise, Resvec, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->numFitEFAC,((MNStruct *)context)->numFitEQUAD);

		det=likeInfo[0];
		Chisq=likeInfo[1];
	}	
	else if(((MNStruct *)context)->numFitEFAC == 1 || ((MNStruct *)context)->numFitEQUAD==1 && ((MNStruct *)context)->numFitEFAC < 2 && ((MNStruct *)context)->numFitEQUAD < 2){
		Noise=new double[((MNStruct *)context)->Gsize];
		if(((MNStruct *)context)->whitemodel == 0){
			for(int o=0;o<((MNStruct *)context)->Gsize; o++){
				Noise[o]=pow(EFAC[0],2)*((MNStruct *)context)->SVec[o] + EQUAD[0];
				det += log(Noise[o]);
				Noise[o] = 1.0/Noise[o];
	
			}
		}
		else if(((MNStruct *)context)->whitemodel == 1){
			for(int o=0;o<((MNStruct *)context)->Gsize; o++){
				Noise[o]=pow(EFAC[0],2)*(((MNStruct *)context)->SVec[o] + EQUAD[0]);
				det += log(Noise[o]);
				Noise[o] = 1.0/Noise[o];
	
			}
		}
		
		WhiteMarginGPUWrapper_(Noise, Resvec, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->numFitEFAC,((MNStruct *)context)->numFitEQUAD);

		Chisq=likeInfo[1];
	}	
	else if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD == 0){
		Noise=new double[((MNStruct *)context)->pulse->nobs];
		det=((MNStruct *)context)->staticTimeDet;
		
		WhiteMarginGPUWrapper_(Noise, Resvec, likeInfo, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->numFitEFAC,((MNStruct *)context)->numFitEQUAD);

		Chisq=likeInfo[1];
	}

	if(isnan(det) || isinf(det) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
	
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->Gsize*log(2.0*M_PI) + det + Chisq);	


	}

	delete[] EFAC;
	delete[] EQUAD;
	delete[] Noise;
	delete[] Resvec;
	delete[] likeInfo;
	
	//printf("Like, %g %g \n", det, Chisq);


	

}




void vHRedMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
//	printf("start like %i\n", ndim);
	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	double redamp, redalpha, DMamp, DMalpha;
	int pcount=0;
	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	totdims+=2*((MNStruct *)context)->incStep;
	if(((MNStruct *)context)->incRED==1)totdims+=2;
	if(((MNStruct *)context)->incDM==1)totdims+=2;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;
	for(int p=0;p<totdims;p++){
		//printf("DP: %i %g %g \n", p, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1]);
		if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
		Cube[pcount]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[pcount]+((MNStruct *)context)->Dpriors[p][0];
		//printf("Cube: %i %g %g %g \n", pcount, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1], Cube[pcount]);
		pcount++;
		}
	}
	pcount=0;
	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){

				LDparams[p]=Cube[fitcount]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				LDparams[p]=((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
			}
			//printf("LD: %i %Lg \n", p, LDparams[p]);

		}
		pcount=0;
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
			//printf("Res %i %g\n", o, Resvec[o]);
		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
		fitcount=0;

		
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=Cube[fitcount];
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=0;
			}
		}
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time = (double)((MNStruct *)context)->pulse->obsn[o1].bat;
				if( time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
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
//	printf("red: %g %g %g %g\n", redamp, redalpha, DMamp, DMalpha);
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noisevec=new double[((MNStruct *)context)->pulse->nobs];
	
	if(((MNStruct *)context)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	
			if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
				
				Noisevec[o]=(double) (pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0]);
			}
			else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
				Noisevec[o]=(double) (pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
			}
	
		}
	}
	else if(((MNStruct *)context)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	
			if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
				
				Noisevec[o]=(double) EFAC[0]*EFAC[0]*((pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[0]));
			}
			else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
				Noisevec[o]=(double) EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*((pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
			}
	
		}
	}
	 
	



	
	double *likeInfo=new double[2];
	//printf("Entering %g %g \n", SpecParm[0], SpecParm[1]);
	vHRedMarginGPUWrapper_(Resvec, BATvec, DMVec, Noisevec, SpecParm, likeInfo, TNFactorialList, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];

	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->Gsize*log(2.0*M_PI)+covdet+Chisq);	


	}
// 	endClock = clock();
// //   	printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(double)CLOCKS_PER_SEC);
	
//	printf("Like %g \n",lnew);
	delete[] EFAC;
	delete[] EQUAD;
	delete[] BATvec;
	delete[] Noisevec;
	delete[] SpecParm;
	delete[] Resvec;
	delete[] DMVec;


}

void vHRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	double redamp, redalpha, DMamp, DMalpha;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	double NLphase=0;
	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
		}

		NLphase=(double)LDparams[0];
		pcount++;
		for(int p=1;p<((MNStruct *)context)->numFitTiming;p++){
			((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = LDparams[pcount];	
			pcount++;
		}
		for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
			((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]= LDparams[pcount];
			pcount++;
		}
		
	//	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	//	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
		
	//	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
	//		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
	//	}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
	
		for(int p=0;p < numfit; p++){
			Fitparams[p]=Cube[p];
			pcount++;
		}
		
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}

	double **Steps;
	if(((MNStruct *)context)->incStep > 0){
		
		Steps=new double*[((MNStruct *)context)->incStep];
		
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			Steps[i]=new double[2];
			Steps[i][0] = Cube[pcount];
			pcount++;
			Steps[i][1] = Cube[pcount];
			pcount++;
		}
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
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
	
	//DO DMModel stuff if fitting for that
	if(((MNStruct *)context)->doLinear==0){
		if(((MNStruct *)context)->pulse->param[param_dmmodel].fitFlag[0] == 1){
			for(int i =0; i < ((MNStruct *)context)->pulse[0].dmoffsDMnum; i++){
				((MNStruct *)context)->pulse[0].dmoffsDM[i]=Cube[pcount];
				pcount++;
			}
		}

		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			  Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+NLphase;
		}
	}

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Steps[i][0];
			double StepTime = Steps[i][1];

			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time=(double)((MNStruct *)context)->pulse->obsn[o1].bat ;
				if( time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}

	double *BatVec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noise=new double[((MNStruct *)context)->pulse->nobs];

	if(((MNStruct *)context)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			BatVec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
			if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
				Noise[o]=(double) (pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0]);
			}
			else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
				Noise[o]=(double) (pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
			}
	
	
		}
	}
	if(((MNStruct *)context)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			BatVec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
			if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
				Noise[o]=(double) EFAC[0]*EFAC[0]*((pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[0]));
			}
			else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
				Noise[o]=(double) EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*((pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
			}
	
	
		}
	}


	double *likeInfo=new double[2];


	vHRedGPUWrapper_(SpecParm, BatVec,DMVec, Resvec, Noise, likeInfo,((MNStruct *)context)->pulse->nobs);

	double covdet=likeInfo[0];
	double Chisq=likeInfo[1];


	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(((double)((MNStruct *)context)->pulse->nobs)*log(2.0*M_PI) + covdet+Chisq);	
	}


	delete[] EFAC;
	delete[] EQUAD;
	delete[] Resvec;
	delete[] BatVec;
	delete[] Noise;
	delete[] SpecParm;
	delete[] DMVec;
	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			delete[] Steps[i];
		}
		delete[] Steps;
	}

	//printf("GPU Like: %g %g %g \n",lnew,Chisq,covdet);
}


void LRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	//printf("hereNM");
	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];

	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}


	if(((MNStruct *)context)->doLinear==0){
	
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
		
		if(((MNStruct *)context)->pulse->param[param_dmmodel].fitFlag[0] == 1){
			int DMnum=((MNStruct *)context)->pulse[0].dmoffsDMnum;
			for(int i =0; i < DMnum; i++){
				((MNStruct *)context)->pulse[0].dmoffsDM[i]=Cube[ndim-DMnum+i];
			}
		}
		
		
		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
		
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
	
		for(int p=0;p < numfit; p++){
			Fitparams[p]=Cube[p];
			pcount++;
		}
		
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time = (double)((MNStruct *)context)->pulse->obsn[o1].bat;
				if(time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
        }

  	
        int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);

        int totCoeff=0;
        if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
        if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

        double *powercoeff=new double[totCoeff];
        for(int o=0;o<totCoeff; o++){
                powercoeff[o]=0;
        }

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];

	double tdet=0;
	double timelike=0;



	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];

	if(((MNStruct *)context)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
				WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
				
				tdet=tdet+log(WorkNoise[o]);
				WorkNoise[o]=1.0/WorkNoise[o];
				timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	
				BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	
	
		}
	}
	else if(((MNStruct *)context)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
				WorkNoise[o]=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
				
				tdet=tdet+log(WorkNoise[o]);
				WorkNoise[o]=1.0/WorkNoise[o];
				timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	
				BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
	
	
		}
	}

	double *NFd = new double[totCoeff];
        double **FNF=new double*[totCoeff];
        for(int i=0;i<totCoeff;i++){
                FNF[i]=new double[totCoeff];
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


	double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
	

        if(((MNStruct *)context)->incRED==2){
        
        	if(((MNStruct *)context)->incFloatRed == 0){
                for (int i=0; i<FitRedCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
            }
            else if(((MNStruct *)context)->incFloatRed >0){

                for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
                
                		int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
                
                 for (int i=FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed; i<FitRedCoeff/2; i++){
                		//printf("Freq: %g \n", Cube[pcount]);
                		
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
                }
                
            }
                
                
                startpos=FitRedCoeff;

        }
        else if(((MNStruct *)context)->incRED==3){

			freqdet=0;
			
			for(int pl = 0; pl < ((MNStruct *)context)->numFitRedPL; pl ++){
			
                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;

                
   			   redamp=pow(10.0, redamp);

				for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){

		                freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
		                freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
		                
		                double PLcomp=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);

		                powercoeff[i]+= PLcomp;
		                powercoeff[i+FitRedCoeff/2]+= PLcomp;
				}
			}
				
				
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
				freqdet=freqdet+2*log(powercoeff[i]);
		//		printf("%i %g %g \n",i,powercoeff[i], freqdet);
			}

                 for (int i=FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed; i<FitRedCoeff/2; i++){
                                
                    //    Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                startpos=FitRedCoeff;

        }
// 		printf("DM\n");
        double nlist[((MNStruct *)context)->incFloatDM][2];
	if(((MNStruct *)context)->incDM==2){
        
        	if(((MNStruct *)context)->incFloatDM == 0){

                for (int i=0; i<FitDMCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }
            }
            else if(((MNStruct *)context)->incFloatDM >0){

                for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }
                
                for (int i=FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM; i<FitDMCoeff/2; i++){

                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        
						
                        
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
                }




                
            }            

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


        }
        else if(((MNStruct *)context)->incDM==3){
        
        
				for(int pl = 0; pl < ((MNStruct *)context)->numFitDMPL; pl ++){
			
		            double DMamp=Cube[pcount];
		            pcount++;
		            double DMindex=Cube[pcount];
		            pcount++;

		            DMamp=pow(10.0, DMamp);

					for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){

				            freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
		                    freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
		                    
		                     double PLcomp=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);
				            
				            powercoeff[startpos+i]+=PLcomp;
		                    powercoeff[startpos+i+FitDMCoeff/2]+=PLcomp;
					}
				}
			
				for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
					freqdet=freqdet+2*log(powercoeff[startpos+i]);
//					printf("%i %g %g \n", i, powercoeff[startpos+i], freqdet);
				}
			



                 for (int i= FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i<FitDMCoeff/2; i++){

						//Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

        }
	LRedGPUWrapper_(freqs, Resvec, BATvec, DMVec, WorkNoise, FNF, NFd, ((MNStruct *)context)->pulse->nobs, FitRedCoeff, FitDMCoeff, totCoeff,((MNStruct *)context)->incRED, ((MNStruct *)context)->incDM);



        double **PPFM=new double*[totCoeff];
        for(int i=0;i<totCoeff;i++){
                PPFM[i]=new double[totCoeff];
                for(int j=0;j<totCoeff;j++){
                        PPFM[i][j]=0;
                }
        }

        for(int c1=0; c1<totCoeff; c1++){
                PPFM[c1][c1]=1.0/powercoeff[c1];
        }



        for(int j=0;j<totCoeff;j++){
                for(int k=0;k<totCoeff;k++){
                        PPFM[j][k]=PPFM[j][k]+FNF[j][k];
                }
        }
        
        double jointdet=0;
        double freqlike=0;
       double *WorkCoeff = new double[totCoeff];
       for(int o1=0;o1<totCoeff; o1++){
                WorkCoeff[o1]=NFd[o1];
        }



        int info=0;
        dpotrfInfo(PPFM, totCoeff, jointdet, info);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
        }
        lnew=-0.5*(((double)((MNStruct *)context)->pulse->nobs)*log(2.0*M_PI) + tdet+jointdet+freqdet+timelike-freqlike);


	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);

		
	}
 	//printf("Like: %g %g %g %g %g %g\n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
	 //printf("CPULIKE: %g %g %g %g %g %g \n", lnew, jointdet,tdet,freqdet,timelike,freqlike);

	delete[] EFAC;
	delete[] EQUAD;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] BATvec;
	delete[] NFd;
	delete[] freqs;
	delete[] DMVec;
	delete[] WorkCoeff;

	for (int j = 0; j < totCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;



	for (int j = 0; j < totCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;


}




void LRedMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

//	printf("hereM");
	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;
        int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;
        if(((MNStruct *)context)->incRED==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
        if(((MNStruct *)context)->incDM==2)totdims+=((MNStruct *)context)->numFitDMCoeff;
        if(((MNStruct *)context)->incRED==3)totdims+=2*((MNStruct *)context)->numFitRedPL;
        if(((MNStruct *)context)->incDM==3)totdims+=2*((MNStruct *)context)->numFitDMPL;
	if(((MNStruct *)context)->incFloatDM != 0)totdims+=2*((MNStruct *)context)->incFloatDM;
	if(((MNStruct *)context)->incFloatRed != 0)totdims+=2*((MNStruct *)context)->incFloatRed;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;
	int reject=0;
	for(int p=0;p<totdims;p++){
		//printf("DP: %i %g %g \n", p, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1]);
		if(((MNStruct *)context)->Dpriors[p][0] != ((MNStruct *)context)->Dpriors[p][1]){
		//	if(((MNStruct *)context)->incFloatRed > 1 && p > ((MNStruct *)context)->FloatRedstart && (p-((MNStruct *)context)->FloatRedstart)%2 ==0 && p < ((MNStruct *)context)->FloatRedstart+2*((MNStruct *)context)->incFloatRed) {
                        //printf("%i %i %i %i %i\n", p,((MNStruct *)context)->incFloatRed,((MNStruct *)context)->FloatRedstart,(p-((MNStruct *)context)->FloatRedstart)%2,((MNStruct *)context)->FloatRedstart+2*((MNStruct *)context)->incFloatRed);

		//		((MNStruct *)context)->Dpriors[p][0] = Cube[pcount-2]+0.1;
		//		if(((MNStruct *)context)->Dpriors[p][0] > ((MNStruct *)context)->Dpriors[p][1]){reject=1;}
		//	}
		Cube[pcount]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[pcount]+((MNStruct *)context)->Dpriors[p][0];
		//printf("DP: %i %i %g %g %g \n", p, pcount, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1], Cube[pcount]);
		pcount++;
		}
	}
	pcount=0;
	if(((MNStruct *)context)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){

				LDparams[p]=Cube[fitcount]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				LDparams[p]=((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
			}
			//printf("LD: %i %Lg \n", p, LDparams[p]);

		}
		pcount=0;
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
			//printf("Res %i %g\n", o, Resvec[o]);
		}
	
	}
	else if(((MNStruct *)context)->doLinear==1){
		fitcount=0;

		
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=Cube[fitcount];
				fitcount++;
			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				Fitparams[p]=0;
			}
// 			printf("%i %g \n",p,Fitparams[p]);
		}
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
// 			printf("%i %g %g \n", o, (double)((MNStruct *)context)->pulse->obsn[o].residual,Fitvec[o]);
// 			for(int p=0;p<8; p++){
// 				printf("DM %i %g \n",p,((MNStruct *)context)->DMatrix[o][p]);
// 			}
			
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;

	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				if(((MNStruct *)context)->pulse->obsn[o1].bat > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}	


	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->systemcount];
		for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
			//printf("EQ %i %g \n", o,  EQUAD[o]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
        }


	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);
	
    if(((MNStruct *)context)->incFloatDM != 0)FitDMCoeff+=2*((MNStruct *)context)->incFloatDM;
    if(((MNStruct *)context)->incFloatRed != 0)FitRedCoeff+=2*((MNStruct *)context)->incFloatRed;

    int totCoeff=0;
    if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
    if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;
   

	double *powercoeff=new double[totCoeff];
	for(int o=0;o<totCoeff; o++){
		powercoeff[o]=0;
	}

	double tdet=0;
	double timelike=0;
	double temptimedet=0;

	double *Noise;
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	//printf("white noise %i %i \n",((MNStruct *)context)->numFitEFAC,((MNStruct *)context)->numFitEQUAD);
	if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD==0){
		//printf("loop2\n");
		Noise=new double[((MNStruct *)context)->pulse->nobs];
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2);
			BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
			
		}
		
	}
	else if(((MNStruct *)context)->numFitEFAC == 1 || ((MNStruct *)context)->numFitEQUAD==1 && ((MNStruct *)context)->numFitEFAC <2 && ((MNStruct *)context)->numFitEQUAD<2){

		Noise=new double[((MNStruct *)context)->Gsize];
		if(((MNStruct *)context)->whitemodel == 0){
			for(int o=0;o<((MNStruct *)context)->Gsize; o++){
				//	  printf("%i %g %g\n",o, EFAC[0], EQUAD);
				Noise[o]=pow(EFAC[0],2)*((MNStruct *)context)->SVec[o] + EQUAD[0];
				temptimedet += log(Noise[o]);
				Noise[o] = 1.0/Noise[o];
			
			}
		}
		if(((MNStruct *)context)->whitemodel == 1){
			for(int o=0;o<((MNStruct *)context)->Gsize; o++){

				Noise[o]=pow(EFAC[0],2)*(((MNStruct *)context)->SVec[o] + EQUAD[0]);
				temptimedet += log(Noise[o]);
				Noise[o] = 1.0/Noise[o];
			
			}
		}

		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		}
		
	}
	
	else if(((MNStruct *)context)->numFitEFAC > 1 || ((MNStruct *)context)->numFitEQUAD > 1){
		//printf("loop2\n");
		Noise=new double[((MNStruct *)context)->pulse->nobs];

		if(((MNStruct *)context)->whitemodel == 0){
			for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
				Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
				BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
				
			}
		}
		if(((MNStruct *)context)->whitemodel == 1){
			for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
				Noise[o]=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
				BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
				
			}
		}

	}

// 	printf("Noise stuff\n");
	double *NFd = new double[totCoeff];
	double **FNF=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		FNF[i]=new double[totCoeff];
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


		double maxtspan=1*(end-start);

        double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
		int tooclose=0;
		
        if(((MNStruct *)context)->incRED==2){
        
        	if(((MNStruct *)context)->incFloatRed == 0){
                for (int i=0; i<FitRedCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
            }
            else if(((MNStruct *)context)->incFloatRed >0){

                for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
                
                		int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
                
                 for (int i=FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed; i<FitRedCoeff/2; i++){
                		//printf("Freq: %g \n", Cube[pcount]);
                		
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
                }
                
            }
                
                
                startpos=FitRedCoeff;

        }
        else if(((MNStruct *)context)->incRED==3){

		freqdet=0;
		
		for(int pl = 0; pl < ((MNStruct *)context)->numFitRedPL; pl ++){
			
			double redamp=Cube[pcount];
			pcount++;
			double redindex=Cube[pcount];
			pcount++;

   			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
    			
			
			redamp=pow(10.0, redamp);

			double Agw=redamp;
	
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){

		

// 				printf("max tspan: %g \n", maxtspan*(24*60*60));
	
				freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/(maxtspan);
				freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
				
				double PLcomp=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);
				powercoeff[i]+= PLcomp;
				powercoeff[i+FitRedCoeff/2]+= PLcomp;


// 				double rho = (Agw*Agw/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-redindex))/(maxtspan*24*60*60);
// 				double rho = (Agw*Agw)*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-redindex))/(maxtspan*24*60*60);
// 				powercoeff[i]+= rho;
// 				powercoeff[i+FitRedCoeff/2]+= rho;
			}
		}
				
				
		for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
			freqdet=freqdet+2*log(powercoeff[i]);
		}

                 for (int i=FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed; i<FitRedCoeff/2; i++){
                                
                    //    Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                startpos=FitRedCoeff;

        }
// 		printf("DM\n");
        double nlist[((MNStruct *)context)->incFloatDM][2];
	if(((MNStruct *)context)->incDM==2){
        
        	if(((MNStruct *)context)->incFloatDM == 0){

			for (int i=0; i<FitDMCoeff/2; i++){
				int pnum=pcount;
				double pc=Cube[pcount];
				freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed+i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
	
				powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
				powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
				pcount++;
			}
           	 }
           	else if(((MNStruct *)context)->incFloatDM >0){

			for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
				int pnum=pcount;
				double pc=Cube[pcount];
				freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
	
				powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
				powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
				pcount++;
			}
                
			for (int i=FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM; i<FitDMCoeff/2; i++){
	
				freqs[startpos+i]=Cube[pcount]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
	
				pcount++;
				
				int pnum=pcount;
				double pc=Cube[pcount];
				pcount++;
				
							
				
				
	
				powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
				powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
			
			}




                
		}            

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


        }
        else if(((MNStruct *)context)->incDM==3){
        
        
		for(int pl = 0; pl < ((MNStruct *)context)->numFitDMPL; pl ++){
			
			double DMamp=Cube[pcount];
			pcount++;
			double DMindex=Cube[pcount];
			pcount++;

			double f1yr = 1.0/3.16e7;
			DMamp=pow(10.0, DMamp);

			for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
	
				freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
				
				double PLcomp=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);
				powercoeff[startpos+i]+= PLcomp;
				powercoeff[startpos+i+FitDMCoeff/2]+= PLcomp;


// 				double rho = (DMamp*DMamp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMindex))/(maxtspan*24*60*60);	
// 				powercoeff[startpos+i]+=rho;
// 				powercoeff[startpos+i+FitDMCoeff/2]+=rho;
			}
		}
		
		for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i++){
			freqdet=freqdet+2*log(powercoeff[startpos+i]);
		}

                 for (int i= FitDMCoeff/2 - ((MNStruct *)context)->incFloatDM ; i<FitDMCoeff/2; i++){

						//Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

        }




// 	printf("Entering gpu \n");
	double *likeVals=new double[2];
	LRedMarginGPUWrapper_(freqs, Resvec, BATvec, DMVec, Noise, FNF, NFd, likeVals, ((MNStruct *)context)->pulse->nobs, FitRedCoeff,FitDMCoeff, totCoeff,((MNStruct *)context)->Gsize, ((MNStruct *)context)->incRED,((MNStruct *)context)->incDM,((MNStruct *)context)->numFitEFAC,((MNStruct *)context)->numFitEQUAD);
	
	tdet=likeVals[0];
	timelike=likeVals[1];
	
	if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD == 0){
			tdet=((MNStruct *)context)->staticTimeDet;
	}
	else if(((MNStruct *)context)->numFitEFAC == 1 || ((MNStruct *)context)->numFitEQUAD==1 && ((MNStruct *)context)->numFitEFAC <2 && ((MNStruct *)context)->numFitEQUAD<2){
			tdet=temptimedet;
	}

	double **PPFM=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		PPFM[i]=new double[totCoeff];
		for(int j=0;j<totCoeff;j++){
			PPFM[i][j]=0;
		}
	}

	for(int c1=0; c1<totCoeff; c1++){
		PPFM[c1][c1]=1.0/powercoeff[c1];
// 		printf("PPRM: %i %g\n", c1, PPFM[c1][c1]);
	}



	for(int j=0;j<totCoeff;j++){
		for(int k=0;k<totCoeff;k++){
			//printf("FNF: %i %i %g  %g\n", j,k,PPFM[j][k],FNF[j][k]);		
			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}
	
	double jointdet=0;
	double freqlike=0;
       double *WorkCoeff = new double[totCoeff];
       for(int o1=0;o1<totCoeff; o1++){
                WorkCoeff[o1]=NFd[o1];
		//printf("NFd: %i %g \n", o1, WorkCoeff[o1]);
        }



	int info=0;
	dpotrfInfo(PPFM, totCoeff, jointdet, info);
	dpotrs(PPFM, WorkCoeff, totCoeff);
	for(int j=0;j<totCoeff;j++){
		freqlike += NFd[j]*WorkCoeff[j];
	//	printf("Freqlike : %i %g %g %g \n", j, freqlike, NFd[j], WorkCoeff[j]);
	}

	lnew=-0.5*(((double)((MNStruct *)context)->Gsize)*log(2.0*M_PI) + tdet+jointdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew) || reject == 1){

		lnew=-pow(10.0,200);
	
	}

	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] powercoeff;
	delete[] NFd;
	delete[] DMVec;
	for (int j = 0; j < totCoeff; j++){
		delete[]PPFM[j];
	}
	delete[]PPFM;


	for (int j = 0; j < totCoeff; j++){
		delete[]FNF[j];
	}
	delete[]FNF;

	delete[] Noise;
	delete[] Resvec;
	delete[] likeVals;
	delete[] freqs;
	delete[] BATvec;
	
 	//printf("GPULike: %g %g %g %g %g %g\n",lnew,timelike, tdet, freqlike, jointdet, freqdet);
//
	// printf("CPULIKE: %g %g %g %g %g %g \n", lnew, jointdet,tdet,freqdet,timelike,freqlike);

}





void vHRedDMMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD, redamp, redalpha, dmamp, dmalpha;
	int pcount=0;

	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];

	}

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
	
	
	dmamp=Cube[pcount];
	pcount++;
	dmalpha=Cube[pcount];
	pcount++;
  	

// 	startClock = clock();
	
	
	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noisevec=new double[((MNStruct *)context)->pulse->nobs];
	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	
	double DMKappa = 2.410*pow(10.0,-16);


	
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
	
		Res[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual;
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		Noisevec[o]=(double)(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD);
		DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
				//printf("%i %g %g %g %g \n",o,BATvec[o],Res[o],Noisevec[o],DMVec[o]);
	}
	 
	double *SpecParm=new double[6];

	SpecParm[0]=redamp;
	SpecParm[1]=redalpha;
	SpecParm[2]=dmamp;
	SpecParm[3]=dmalpha;
	
	double *likeInfo=new double[2];
	//printf("Entering %g %g \n", SpecParm[0], SpecParm[1]);
	vHRedDMMarginGPUWrapper_(Res, BATvec, Noisevec, DMVec, SpecParm, likeInfo, TNFactorialList, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize);

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
	delete[] DMVec;
	//printf("Like: %g %g %g \n",lnew,Chisq,covdet);

}


void vHRedDMGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD, redamp, redalpha, dmamp, dmalpha;
	int pcount=0;

	for(int p=0;p<ndim;p++){

			Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	
	}

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
	
	dmamp=Cube[pcount];
	pcount++;
	dmalpha=Cube[pcount];
	pcount++;

	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);                /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
	
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	
	double *BatVec=new double[((MNStruct *)context)->pulse->nobs];
	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	
	double DMKappa = 2.410*pow(10.0,-16);
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Res[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual;
		BatVec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));

	}


	double *likeInfo=new double[2];
	double *SpecParm=new double[6];

	SpecParm[0]=redamp;
	SpecParm[1]=redalpha;
	SpecParm[2]=dmamp;
	SpecParm[3]=dmalpha;
	
	vHRedDMGPUWrapper_(SpecParm, BatVec,Res, Noise, DMVec, likeInfo,((MNStruct *)context)->pulse->nobs);

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
	delete[] SpecParm;
	delete[] DMVec;

		//printf("GPU Like: %g %g %g \n",lnew,Chisq,covdet);
}



void LRedDMMarginGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

//	printf("hereM");
	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD;
	int pcount=0;

	for(int p=0;p<ndim;p++){
		Cube[p]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[p]+((MNStruct *)context)->Dpriors[p][0];
	}

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

	int RedCoeff=(((MNStruct *)context)->numFitRedCoeff);
	int DMCoeff=(((MNStruct *)context)->numFitRedCoeff);
	int totCoeff=RedCoeff+DMCoeff;
	double *powercoeff=new double[totCoeff];



	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *Res=new double[((MNStruct *)context)->pulse->nobs];
	double *BATvec=new double[((MNStruct *)context)->pulse->nobs];
	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	
	double DMKappa = 2.410*pow(10.0,-16);

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
		Res[o]=((MNStruct *)context)->pulse->obsn[o].residual;
		BATvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].bat;
		DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
	}

	double *NFd = new double[totCoeff];
	double **FNF=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		FNF[i]=new double[totCoeff];
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
	for (int i=0; i<totCoeff/2; i++){
		int pnum=pcount;
		double pc=Cube[pcount];
		
		powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
		powercoeff[i+totCoeff/2]=powercoeff[i];
		freqdet=freqdet+2*log(powercoeff[i]);
		pcount++;
	}


	int coeffsize=totCoeff/4;
	double *freqs = new double[(totCoeff/4)];
	for(int i=0;i<totCoeff/4;i++){
		freqs[i]=double(i+1)/maxtspan;
	}	
	
	double *likeVals=new double[2];
	LRedDMMarginGPUWrapper_(freqs, Res, BATvec, DMVec,  Noise, FNF, NFd, likeVals, ((MNStruct *)context)->pulse->nobs, totCoeff,((MNStruct *)context)->Gsize);
	
	double tdet=likeVals[0];
	double timelike=likeVals[1];
	
	



	double **PPFM=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		PPFM[i]=new double[totCoeff];
		for(int j=0;j<totCoeff;j++){
			PPFM[i][j]=0;
		}
	}


	for(int c1=0; c1<totCoeff; c1++){
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}



	for(int j=0;j<totCoeff;j++){
		for(int k=0;k<totCoeff;k++){
			
			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}

 	
	double jointdet=0;
	dpotrf(PPFM, totCoeff, jointdet);
    dpotri(PPFM,totCoeff);

	double freqlike=0;
	for(int i=0;i<totCoeff;i++){
		for(int j=0;j<totCoeff;j++){
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

	for (int j = 0; j < totCoeff; j++){
		delete[]PPFM[j];
	}
	delete[]PPFM;


	for (int j = 0; j < totCoeff; j++){
		delete[]FNF[j];
	}
	delete[]FNF;

	delete[] Noise;
	delete[] Res;
	
	printf("GPULike: %g %g %g %g %g %g\n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
//

}


