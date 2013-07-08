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


#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "dpotrs.h"
#include "tempo2.h"
#include "TempoNest.h"


void vHRedLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD, redamp, redalpha, DMamp, DMalpha;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual+Fitvec[o];
		}
		
		delete[] Fitvec;
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

	double redampsquared=0;
	double redcovconst=0;
	
	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;
	
	if(((MNStruct *)context)->incRED==1){
		redamp=Cube[pcount];
		pcount++;
		redalpha=Cube[pcount];
		pcount++;
		
		redamp=pow(10.0,redamp);
		redampsquared=redamp*redamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-redalpha)))/(pow(flo,(redalpha-1)));
		redcovconst=gsl_sf_gamma(1-redalpha)*sin(0.5*M_PI*redalpha);
	}
	
	double DMampsquared=0;
	double DMcovconst=0;

  	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);

	if(((MNStruct *)context)->incDM==1){
		DMamp=Cube[pcount];
		pcount++;
		DMalpha=Cube[pcount];
		pcount++;
		
		DMamp=pow(10.0,DMamp);
		DMampsquared=DMamp*DMamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-DMalpha)))/(pow(flo,(DMalpha-1)));
		DMcovconst=gsl_sf_gamma(1-DMalpha)*sin(0.5*M_PI*DMalpha);

		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
		}
	}


	double timdiff=0;

	
	double **CovMatrix = new double*[((MNStruct *)context)->pulse->nobs]; 
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs;o1++){
		 CovMatrix[o1]=new double[((MNStruct *)context)->pulse->nobs];
		 for(int o2=0;o2<((MNStruct *)context)->pulse->nobs;o2++){
		 	CovMatrix[o1][o2]=0;
		 }
	}
	
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
		CovMatrix[o1][o1] += pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD;
	}

		

	if(((MNStruct *)context)->incRED==1){
		for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
			for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
		
				timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
				double tau=2.0*M_PI*fabs(timdiff);
			
			
				double redcovsum=0;

				for(int k=0; k <=4; k++){
					redcovsum=redcovsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-redalpha));

				}

				CovMatrix[o1][o2] += redampsquared*(redcovconst*pow((flo*tau),(redalpha-1)) - redcovsum);
			}
		}
	}
	
	if(((MNStruct *)context)->incDM==1){
		for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
			for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
		
				timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
				double tau=2.0*M_PI*fabs(timdiff);
				double DMcovsum=0;

				for(int k=0; k <=4; k++){
					DMcovsum=DMcovsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-DMalpha));

				}

				CovMatrix[o1][o2] += DMampsquared*(DMcovconst*pow((flo*tau),(DMalpha-1)) - DMcovsum)*DMVec[o1]*DMVec[o2];;
			}
		}
	}	

	double covdet=0;
	double *WorkRes = new double[((MNStruct *)context)->pulse->nobs];
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
		WorkRes[o1]=Resvec[o1];
	}
	dpotrf(CovMatrix, ((MNStruct *)context)->pulse->nobs, covdet);
    dpotrs(CovMatrix, WorkRes, ((MNStruct *)context)->pulse->nobs);


	double Chisq=0;


	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
		Chisq += Resvec[o1]*WorkRes[o1];
	}

	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->pulse->nobs*log(2*M_PI) + covdet + Chisq);	

	}
	
	//printf("Like: %g %g %g\n", lnew, covdet, Chisq);
	delete[] EFAC;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs;o++){delete[] CovMatrix[o];}
	delete[] CovMatrix;
	delete[] WorkRes;
	delete[] Resvec;


}

void vHRedMarginLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD, redamp, redalpha, DMamp, DMalpha;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual+Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;	
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


	double redampsquared=0;
	double redcovconst=0;
	
	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;
	
	if(((MNStruct *)context)->incRED==1){
		redamp=Cube[pcount];
		pcount++;
		redalpha=Cube[pcount];
		pcount++;
		redamp=pow(10.0,redamp);
		redampsquared=redamp*redamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-redalpha)))/(pow(flo,(redalpha-1)));
		redcovconst=gsl_sf_gamma(1-redalpha)*sin(0.5*M_PI*redalpha);
	}
	
	double DMampsquared=0;
	double DMcovconst=0;

  	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);

	if(((MNStruct *)context)->incDM==1){
		DMamp=Cube[pcount];
		pcount++;
		DMalpha=Cube[pcount];
		pcount++;

		
		DMamp=pow(10.0,DMamp);
		DMampsquared=DMamp*DMamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-DMalpha)))/(pow(flo,(DMalpha-1)));
		DMcovconst=gsl_sf_gamma(1-DMalpha)*sin(0.5*M_PI*DMalpha);

		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
		}
	}


	double timdiff=0;

	
	double **CovMatrix = new double*[((MNStruct *)context)->pulse->nobs]; 
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs;o1++){
		 CovMatrix[o1]=new double[((MNStruct *)context)->pulse->nobs];
		 for(int o2=0;o2<((MNStruct *)context)->pulse->nobs;o2++){
		 	CovMatrix[o1][o2]=0;
		 }
	}
	
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
		CovMatrix[o1][o1] += pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD;
	}

		

	if(((MNStruct *)context)->incRED==1){
		for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
			for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
		
				timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
				double tau=2.0*M_PI*fabs(timdiff);
			
			
				double redcovsum=0;

				for(int k=0; k <=4; k++){
					redcovsum=redcovsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-redalpha));

				}

				CovMatrix[o1][o2] += redampsquared*(redcovconst*pow((flo*tau),(redalpha-1)) - redcovsum);
			}
		}
	}
	
	if(((MNStruct *)context)->incDM==1){
		for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
			for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
		
				timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
				double tau=2.0*M_PI*fabs(timdiff);
				double DMcovsum=0;

				for(int k=0; k <=4; k++){
					DMcovsum=DMcovsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-DMalpha));

				}

				CovMatrix[o1][o2] += DMampsquared*(DMcovconst*pow((flo*tau),(DMalpha-1)) - DMcovsum)*DMVec[o1]*DMVec[o2];;
			}
		}
	}	

	double **CG = new double*[((MNStruct *)context)->pulse->nobs]; for(int o1=0;o1<((MNStruct *)context)->pulse->nobs;o1++)CG[o1]=new double[((MNStruct *)context)->Gsize];

	double **GCG= new double*[((MNStruct *)context)->Gsize]; for(int o1=0;o1<((MNStruct *)context)->Gsize;o1++)GCG[o1]=new double[((MNStruct *)context)->Gsize];
	
	double *GRes=new double[((MNStruct *)context)->Gsize];


	dgemm(CovMatrix,((MNStruct *)context)->GMatrix, CG, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'N', 'N');

	dgemm(((MNStruct *)context)->GMatrix,CG, GCG, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'T', 'N');

	dgemv(((MNStruct *)context)->GMatrix,Resvec,GRes,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'T');

	double covdet=0;
	double *WorkGRes = new double[((MNStruct *)context)->Gsize];
	for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
		WorkGRes[o1]=GRes[o1];
	}
	dpotrf(GCG, ((MNStruct *)context)->Gsize, covdet);
    dpotrs(GCG, WorkGRes, ((MNStruct *)context)->Gsize);



	double Chisq=0;


	for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
		Chisq += GRes[o1]*WorkGRes[o1];
	}

	if(isnan(covdet) || isinf(covdet) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->Gsize*log(2*M_PI) + covdet + Chisq);	

	}


	delete[] EFAC;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs;o++){delete[] CovMatrix[o];}
	delete[] CovMatrix;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs;o++){delete[] CG[o];}
	delete[] CG;
	for(int o=0;o<((MNStruct *)context)->Gsize;o++){delete[] GCG[o];}
	delete[] GCG;


	delete[] Resvec;
	delete[] WorkGRes;
	delete[] GRes;
	
}


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
			//printf("FitP: %i %g \n",p,Cube[p]);
			pcount++;
		}
		
		double *Fitvec=new double[((MNStruct *)context)->pulse->nobs];

		dgemv(((MNStruct *)context)->DMatrix,Fitparams,Fitvec,((MNStruct *)context)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual+Fitvec[o];
			//printf("FitVec: %i %g \n",o,Fitvec[o]);
		}
		
		delete[] Fitvec;
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


void WhiteMarginLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;


	for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
		LDparams[p]=Cube[p]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual+Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;
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

	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *GRes=new double[((MNStruct *)context)->Gsize];

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
	}

	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[((MNStruct *)context)->Gsize];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<((MNStruct *)context)->Gsize; j++){
			NG[i][j]=((MNStruct *)context)->GMatrix[i][j]*Noise[i];
		}
	}

	double** GG = new double*[((MNStruct *)context)->Gsize]; for (int k=0; k<((MNStruct *)context)->Gsize; k++) GG[k] = new double[((MNStruct *)context)->Gsize];

	dgemm(((MNStruct *)context)->GMatrix, NG,GG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'T','N');


	dgemv(((MNStruct *)context)->GMatrix,Resvec,GRes,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'T');



	double det=0;
	double *WorkGRes = new double[((MNStruct *)context)->Gsize];
	for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
		WorkGRes[o1]=GRes[o1];
	}
	dpotrf(GG, ((MNStruct *)context)->Gsize, det);
    dpotrs(GG, WorkGRes, ((MNStruct *)context)->Gsize);



	double Chisq=0;
	for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
		Chisq += GRes[o1]*WorkGRes[o1];
	}



	if(isnan(det) || isinf(det) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
		
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->Gsize*log(2.0*M_PI) + det + Chisq);	

	}

	delete[] EFAC;
	delete[] Noise;
	delete[] Resvec;
	delete[] GRes;
	delete[] WorkGRes;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] NG[j];
	}
	delete[] NG;

	for (int j = 0; j < ((MNStruct *)context)->Gsize; j++){
		delete[] GG[j];
	}
	delete[] GG;

	
	//printf("Chisq: %g %g \n",det, Chisq);


}


void LRedLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;
	int bad=0;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];


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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual+Fitvec[o];
		}
		
		delete[] Fitvec;
	}
		

	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double[1];
		EFAC[0]=1;
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double[1];
		EFAC[0]=Cube[pcount];
		if(EFAC[0]<0)bad=1;
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)context)->numFitEFAC];
		for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			EFAC[p]=Cube[pcount];
			if(EFAC[p]<0.1)bad=1;
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



	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
        int totCoeff=0;
        if(((MNStruct *)context)->incRED != 0)totCoeff+=FitCoeff;
        if(((MNStruct *)context)->incDM != 0)totCoeff+=FitCoeff;

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[totCoeff];

	double tdet=0;
	double timelike=0;


	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
			//printf("Tlike: %i %g %g %g %g %g\n",o,timelike,(double)((MNStruct *)context)->pulse->obsn[o].residual,WorkNoise[o],tdet,EFAC[((MNStruct *)context)->sysFlags[o]]);
	}

	double *NFd = new double[totCoeff];
	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
	}

	double **NF=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		NF[i]=new double[totCoeff];
	}

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
// 	printf("Total time span = %.6f days = %.6f years\n",end-start,(end-start)/365.25);
	double maxtspan=end-start;

        double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
        if(((MNStruct *)context)->incRED==2){
                for (int i=0; i<FitCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[i]=double(i+1)/maxtspan;
                        freqs[i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }


                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }



                startpos=FitCoeff;

        }
   else if(((MNStruct *)context)->incRED==3){

                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;
		//printf("red: %g %g \n", redamp, redindex);
                 redamp=pow(10.0, redamp);

                freqdet=0;
                 for (int i=0; i<FitCoeff/2; i++){

                        freqs[i]=double(i+1)/maxtspan;
                        freqs[i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[i]=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);


                 }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }


                startpos=FitCoeff;

        }


       if(((MNStruct *)context)->incDM==2){

                for (int i=0; i<FitCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=double(i+1)/maxtspan;
                        freqs[startpos+i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }



        }
        else if(((MNStruct *)context)->incDM==3){
                double DMamp=Cube[pcount];
                pcount++;
                double DMindex=Cube[pcount];
                pcount++;

                DMamp=pow(10.0, DMamp);

                 for (int i=0; i<FitCoeff/2; i++){
                        freqs[startpos + i]=double(i+1)/maxtspan;
                        freqs[startpos + i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[startpos+i]=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);


                 }
                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }



        }





// 	makeFourier(((MNStruct *)context)->pulse, FitCoeff, FMatrix);

	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<totCoeff;j++){
// 			printf("%i %i %g %g \n",i,j,WorkNoise[i],FMatrix[i][j]);
			NF[i][j]=WorkNoise[i]*FMatrix[i][j];
		}
	}
	dgemv(NF,Resvec,NFd,((MNStruct *)context)->pulse->nobs,totCoeff,'T');
	dgemm(FMatrix, NF , FNF, ((MNStruct *)context)->pulse->nobs, totCoeff, ((MNStruct *)context)->pulse->nobs, totCoeff, 'T', 'N');


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



        dpotrf(PPFM, totCoeff, jointdet);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
        }
	
	lnew=-0.5*(jointdet+tdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
		
	}
	if(bad==1)lnew=-pow(10.0,200);
//	printf("Like: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
	delete[] DMVec;
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] NFd;

	for (int j = 0; j < FitCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] NF[j];
	}
	delete[] NF;

	for (int j = 0; j < FitCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] FMatrix[j];
	}
	delete[] FMatrix;

}

void LRedMarginLogLikeMax(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{


	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;


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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual+Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;	
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


	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
        int totCoeff=0;
        if(((MNStruct *)context)->incRED != 0)totCoeff+=FitCoeff;
        if(((MNStruct *)context)->incDM != 0)totCoeff+=FitCoeff;

	double *powercoeff=new double[totCoeff];


	double *Noise=new double[((MNStruct *)context)->pulse->nobs];
	double *GRes=new double[((MNStruct *)context)->Gsize];

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
	}

	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[((MNStruct *)context)->Gsize];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<((MNStruct *)context)->Gsize; j++){
			//printf("NS: %i %i %g %g \n", i,j,((MNStruct *)context)->GMatrix[i][j], Noise[i]);
			NG[i][j]=((MNStruct *)context)->GMatrix[i][j]*Noise[i];

		}
	}

	double** GG = new double*[((MNStruct *)context)->Gsize]; for (int k=0; k<((MNStruct *)context)->Gsize; k++) GG[k] = new double[((MNStruct *)context)->Gsize];

	dgemm(((MNStruct *)context)->GMatrix, NG,GG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'T','N');

	
	double tdet=0;
	dpotrf(GG, ((MNStruct *)context)->Gsize, tdet);
	dpotri(GG,((MNStruct *)context)->Gsize);
	


	dgemm(((MNStruct *)context)->GMatrix, GG,NG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, 'N','N');

	double **GNG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) GNG[k] = new double[((MNStruct *)context)->pulse->nobs];	

	dgemm(NG, ((MNStruct *)context)->GMatrix, GNG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'N','T');
	
	double *dG=new double[((MNStruct *)context)->pulse->nobs];
	dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
	double timelike=0;
	for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
		timelike+=Resvec[o1]*dG[o1];
	//	printf("timlike: %i %g %g \n", o1, Resvec[o1], dG[o1]);
	}


	double *NFd = new double[totCoeff];
	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
	}

	double **NF=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		NF[i]=new double[totCoeff];
	}

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
                for (int i=0; i<FitCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[i]=double(i+1)/maxtspan;
                        freqs[i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }


                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }



                startpos=FitCoeff;

        }
   else if(((MNStruct *)context)->incRED==3){

                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;
	//	printf("red: %g %g \n", redamp, redindex);
                 redamp=pow(10.0, redamp);

                freqdet=0;
                 for (int i=0; i<FitCoeff/2; i++){

                        freqs[i]=double(i+1)/maxtspan;
                        freqs[i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[i]=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);


                 }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
		//		if(k==0)printf("FM %i %g %g \n", i, freqs[i], time);

                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }


                startpos=FitCoeff;

        }


       if(((MNStruct *)context)->incDM==2){

                for (int i=0; i<FitCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=double(i+1)/maxtspan;
                        freqs[startpos+i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }



        }
        else if(((MNStruct *)context)->incDM==3){
                double DMamp=Cube[pcount];
                pcount++;
                double DMindex=Cube[pcount];
                pcount++;
	//	printf("DM: %g %g \n", DMamp, DMindex);
                DMamp=pow(10.0, DMamp);

                 for (int i=0; i<FitCoeff/2; i++){
                        freqs[startpos + i]=double(i+1)/maxtspan;
                        freqs[startpos + i+FitCoeff/2]=double(i+1)/maxtspan;

                        powercoeff[startpos+i]=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);


                 }
                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitCoeff/2]=sin(2*M_PI*freqs[i]*time)*DMVec[k];
                        }
                }



        }
//	printf("made \n");
	dgemm(GNG, FMatrix , NF, ((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->pulse->nobs, totCoeff, 'N', 'N');

	dgemm(FMatrix, NF , FNF, ((MNStruct *)context)->pulse->nobs, totCoeff, ((MNStruct *)context)->pulse->nobs, totCoeff, 'T', 'N');

	dgemv(NF,Resvec,NFd,((MNStruct *)context)->pulse->nobs,totCoeff,'T');

	double **PPFM=new double*[totCoeff];
	for(int i=0;i<totCoeff;i++){
		PPFM[i]=new double[totCoeff];
		for(int j=0;j<totCoeff;j++){
			PPFM[i][j]=0;
		}
	}


	for(int c1=0; c1<totCoeff; c1++){
//		printf("PC %i %g \n",c1,powercoeff[c1]);
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}



	for(int j=0;j<totCoeff;j++){
		for(int k=0;k<totCoeff;k++){
		        //printf("CPU: %i %i %g %g \n", j,k,PPFM[j][k], FNF[j][k]);	
			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}

        double jointdet=0;
        double freqlike=0;
       double *WorkCoeff = new double[totCoeff];
       for(int o1=0;o1<totCoeff; o1++){
                WorkCoeff[o1]=NFd[o1];
	//	printf("WOrk: %i %g \n", o1, WorkCoeff[o1]);
        }



        dpotrf(PPFM, totCoeff, jointdet);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
	//	printf("freqlike: %i %g %g %g\n",j, freqlike, NFd[j], WorkCoeff[j]);
        }
	
	lnew=-0.5*(tdet+jointdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}

	delete[] DMVec;
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] powercoeff;
	delete[] NFd;
	delete[] dG;

	for (int j = 0; j < FitCoeff; j++){
		delete[]PPFM[j];
	}
	delete[]PPFM;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]NF[j];
	}
	delete[]NF;

	for (int j = 0; j < FitCoeff; j++){
		delete[]FNF[j];
	}
	delete[]FNF;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]FMatrix[j];
	}
	delete[]FMatrix;

	delete[] Noise;
	delete[] Resvec;
	delete[] GRes;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] NG[j];
	}
	delete[] NG;

	for (int j = 0; j < ((MNStruct *)context)->Gsize; j++){
		delete[]GG[j];
	}
	delete[] GG;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] GNG[j];
	}
	delete[] GNG;
	
	//printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}
	
	//printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
// 	printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

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


	if(((MNStruct *)context)->incRED==0 && ((MNStruct *)context)->incDM==0){
		WhiteLogLikeMax(pdParameters, ndims, ndims, lval, context);
	}
	else if(((MNStruct *)context)->incRED==1 || ((MNStruct *)context)->incDM==1){
		vHRedLogLikeMax(pdParameters, ndims, ndims, lval, context);
	}
	else if(((MNStruct *)context)->incRED==2 || ((MNStruct *)context)->incRED==3 || ((MNStruct *)context)->incDM==2 || ((MNStruct *)context)->incDM==3){
		LRedLogLikeMax(pdParameters, ndims, ndims, lval, context);
	}
	

	
	//printf("Like: %g \n",lval);
	// The loglikelihood function is virtual, so the correct one is called
	return -lval;
} // loglikelihood


double neglogMarginlikelihood(const gsl_vector *pvP, void *context) {

	int ndims=((MNStruct *)context)->numdims;
	double pdParameters[ndims];
	int nParameters;
	double lval;

	
	
	// Obtain the pointer to the model
	
	for(int i=0; i<ndims; i++) {
		pdParameters[i] = gsl_vector_get(pvP, i);
  	} // for i


	if(((MNStruct *)context)->incRED==0 && ((MNStruct *)context)->incDM==0){
		WhiteMarginLogLikeMax(pdParameters, ndims, ndims, lval, context);
	}
	else if(((MNStruct *)context)->incRED==1 || ((MNStruct *)context)->incDM==1){
		vHRedMarginLogLikeMax(pdParameters, ndims, ndims, lval, context);
	}
	else if(((MNStruct *)context)->incRED==2 || ((MNStruct *)context)->incRED==3 || ((MNStruct *)context)->incDM==2 || ((MNStruct *)context)->incDM==3){
		LRedMarginLogLikeMax(pdParameters, ndims, ndims, lval, context);
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
	    gsl_vector_set(vStepSize, pcount, 0.2);
   	    gsl_vector_set(vStart, pcount, 1.0);
		pcount++;
	}
  for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
		//printf("Setting EQUD: %i\n",i);
	    gsl_vector_set(vStepSize, pcount, 0.5);
   	    gsl_vector_set(vStart, pcount, -5.0);
		pcount++;
	}
	if(((MNStruct *)context)->incRED == 1 ){
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
	
	if(((MNStruct *)context)->incRED == 3 ){
		//printf("Setting Red\n");
		gsl_vector_set(vStart, pcount, -3.0);
		gsl_vector_set(vStepSize, pcount, 0.5);
		pcount++;
		gsl_vector_set(vStart, pcount, 3.1);
		gsl_vector_set(vStepSize, pcount, 0.25);
		pcount++;
	}
	
	if(((MNStruct *)context)->incDM == 1 ){
		//printf("Setting Red\n");
		gsl_vector_set(vStart, pcount, -10.0);
		gsl_vector_set(vStepSize, pcount, 0.5);
		pcount++;
		gsl_vector_set(vStart, pcount, 3.1);
		gsl_vector_set(vStepSize, pcount, 0.25);
		pcount++;
	}

	if(((MNStruct *)context)->incDM == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
			gsl_vector_set(vStepSize, pcount, 0.5);
			gsl_vector_set(vStart, pcount, -1.0);
			pcount++;
		}
	}	
	
	if(((MNStruct *)context)->incDM == 3 ){
		//printf("Setting Red\n");
		gsl_vector_set(vStart, pcount, -0.0);
		gsl_vector_set(vStepSize, pcount, 0.5);
		pcount++;
		gsl_vector_set(vStart, pcount, 3.1);
		gsl_vector_set(vStepSize, pcount, 0.25);
		pcount++;
	}
	
	if(((MNStruct *)context)->pulse->param[param_dmmodel].fitFlag[0] == 1){
			int DMnum=((MNStruct *)context)->pulse[0].dmoffsDMnum;
			for(int i =0; i < DMnum; i++){
				gsl_vector_set(vStart, pcount, ((MNStruct *)context)->pulse->dmoffsDM[i]);
				gsl_vector_set(vStepSize, pcount, 5*((MNStruct *)context)->pulse->dmoffsDM_error[i]);
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
	if(((MNStruct *)context)->incRED == 1 || ((MNStruct *)context)->incRED == 3){
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
	if(((MNStruct *)context)->incDM == 1 || ((MNStruct *)context)->incDM == 3){
		printf("   Max DM Amp :  %.10g \n", pdParameterEstimates[pcount]);
		pcount++;
		printf("   Max DM Index :  %.10g \n",pdParameterEstimates[pcount]);
		pcount++;
	}

	if(((MNStruct *)context)->incDM == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
			printf("   Max DM Coeff %i :  %.10g \n", i,pdParameterEstimates[pcount]);
			pcount++;
		}
	}	


  gsl_vector_free(vStart);
  gsl_vector_free(vStepSize);
  gsl_multimin_fminimizer_free(pmMinimiser);
  delete[] pdParameterEstimates;
} // NelderMeadOptimum



void getnewMarginDMatrix(pulsar *pulse, int *MarginList, double **TNDM, int **TempoFitNums, int *TempoJumpNums,int TimetoFit, int JumptoFit){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix


		int pcount=1;
		int numToMargin=1;
		
		for (int p=1;p<TimetoFit;p++) {
			if(MarginList[pcount]!=1){
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
			}
			else if(MarginList[pcount]==1){
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < JumptoFit; i++){
			if(MarginList[pcount]!=1){
				pulse[0].fitJump[TempoJumpNums[i]]=0;
			}
			else if(MarginList[pcount]==1){
				pulse[0].fitJump[TempoJumpNums[i]]=1;
				numToMargin++;
			}
		}
	
	
		for(int i=0; i < pulse->nobs; i++) {
			FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin, pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
				TNDM[i][j]=pdParamDeriv[j];
			} 
		} 
	
	
		//Now set fit flags back to how they were
	
		for (int p=1;p<TimetoFit;p++) {
			pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
		}
	
		for(int i=0; i < JumptoFit; i++){
			pulse[0].fitJump[TempoJumpNums[i]]=1;
		}
}


int MLlongturn_hms(long double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  long double sec;

  hh = (int)(turn*24.);
  mm = (int)((turn*24.-hh)*60.);
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (int)((sec*10000. +0.5)/10000);
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

  sprintf(hms," %02d:%02d:%.12Lg",hh,mm,sec);
  return 0;
}

int MLlongturn_dms(long double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  long double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = (int)(trn*360.);
  mm = (int)((trn*360.-dd)*60.);
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (int)((sec*1000. +0.5)/1000);
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
  sprintf(dms,"%c%02d:%02d:%.12Lg",sign,dd,mm,sec);
  return 0;
}

void writeMLPar(std::string longname, pulsar *psr){

	double rms_pre=0.0,rms_post=0.0;
	double mean_pre=0.0,mean_post=0.0,chisqr;
	int i,p,count,k;
	char *fname;
	int npsr=1;

	for(p=0;p<npsr;p++){
		
		std::string parname=longname+"ML.par";

		fname=(char*)parname.c_str();
		//fname="newpar.par";
		printf("writing par file %s. \n",fname);
		FILE *fout2;
		char fname2[1000];
		char str1[100],str2[100],str3[100],str4[100],str5[100];
		int nread;
	
			char hmsstr[100];
			MLlongturn_hms(psr[p].param[param_raj].val[0]/(2*M_PI), hmsstr);
			
	
			strcpy(psr[p].rajStrPost,hmsstr);
				
			MLlongturn_dms(psr[p].param[param_decj].val[0]/(2*M_PI), hmsstr);
			strcpy(psr[p].decjStrPost,hmsstr);
	
		if (strlen(fname)==0)
		{
		printf("Enter filename for new parameter file ");
		scanf("%s",fname2);
		}
		else
		{
		char fname3[1000];
		if (npsr > 1)
			sprintf(fname2,"%s_%d",fname,p+1);
		else
			strcpy(fname2,fname);
		}
		if (!(fout2 = fopen(fname2,"w")))
		{
		printf("Sorry, unable to write to file %s\n",fname2);
		}
		else
		{
		fprintf(fout2,"%-15.15s%s\n","PSRJ",psr[p].name);
		for (i=0;i<MAX_PARAMS;i++)
			{
			for (k=0;k<psr[p].param[i].aSize;k++)
			{
			if (psr[p].param[i].paramSet[k]==1 && i!=param_wave_om 
			&& i!=param_waveepoch && i!=param_ifunc && i!=param_dmmodel &&
			(psr[p].tempo1==0 || (i!=param_dmepoch)))
			{
			if (strcmp(psr[p].param[i].shortlabel[k],"PB")==0 || strcmp(psr[p].param[i].shortlabel[k],"FB0")==0)
				fprintf(fout2,"%-15.15s%s\n","BINARY",psr[p].binaryModel);
			
			if (i == param_raj && psr[p].eclCoord==1)
				fprintf(fout2,"%-15.15s","ELONG");
			else if (i == param_decj && psr[p].eclCoord==1)
				fprintf(fout2,"%-15.15s","ELAT");
			else if (i == param_pmra && psr[p].eclCoord==1)
				fprintf(fout2,"%-15.15s","PMELONG");
			else if (i == param_pmdec && psr[p].eclCoord==1)
				fprintf(fout2,"%-15.15s","PMELAT");		      
			else
				{
				if (i==param_track && psr[p].param[i].val[k]==0)
				{// Do nothing
				}
				else fprintf(fout2,"%-15.15s",psr[p].param[i].shortlabel[k]);
				}
			if (psr[p].eclCoord==0 && i==param_raj)
				fprintf(fout2,"%-25.25s",psr[p].rajStrPost);
			else if (psr[p].eclCoord==1 && i==param_raj)
				fprintf(fout2,"%-25.25Lf",psr[p].param[i].val[k]*180.0/M_PI);
			else if (psr[p].eclCoord==0 && i==param_decj)
				fprintf(fout2,"%-25.25s",psr[p].decjStrPost);
			else if (psr[p].eclCoord==1 && i==param_decj)
				fprintf(fout2,"%-25.25Lf",psr[p].param[i].val[k]*180.0/M_PI);
			else if (i == param_sini && psr[p].param[i].nLinkTo>0){
				fprintf(fout2," KIN\n"); 
				fprintf(fout2,"#SINI\t\t%-25.20Lg",psr[p].param[i].val[k]);
			}
			else if (i==param_tres)
				fprintf(fout2,"%-10.3Lf",psr[p].param[i].val[k]);
			else if (i==param_tzrfrq)
				{
				fprintf(fout2,"%-25.20Lg\n",psr[p].param[i].val[k]);
				if (strcmp(psr[p].tzrsite,"NULL")!=0) 
				fprintf(fout2,"%-15.15s%s","TZRSITE",psr[p].tzrsite);
				}
			else 
				{
				if (i==param_track && psr[p].param[i].val[k]==0)
				{
				// Do nothing
				}
				else
				fprintf(fout2,"%-25.20Lg",psr[p].param[i].val[k]);
				}
			if (psr[p].param[i].fitFlag[k]==1)
				{
				fprintf(fout2," 1 ");
				if (i==param_raj)
				{
				if (psr[p].param[i].err[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].err[k]*12.0*60.0*60.0/M_PI);
				else if (psr[p].param[i].err[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].err[k]*12.0*60.0*60.0/M_PI);
				}
				else if (i==param_decj)
				{
				if (psr[p].param[i].err[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].err[k]*180.0*60.0*60.0/M_PI);
				else if (psr[p].param[i].err[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].err[k]*180.0*60.0*60.0/M_PI);
				}
				else
				{
				if (psr[p].param[i].err[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].err[k]);
				else if (psr[p].param[i].err[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].err[k]);
				}
				}
			else
				{
				fprintf(fout2,"   ");
				if (i==param_raj)
				{
				if (psr[p].param[i].prefitErr[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].prefitErr[k]*12.0*60.0*60.0/M_PI);
				else if (psr[p].param[i].prefitErr[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].prefitErr[k]*12.0*60.0*60.0/M_PI);
				}
				else if (i==param_decj)
				{
				if (psr[p].param[i].prefitErr[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].prefitErr[k]*180.0*60.0*60.0/M_PI);
				else if (psr[p].param[i].prefitErr[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].prefitErr[k]*180.0*60.0*60.0/M_PI);
				}
				else
				{
				if (psr[p].param[i].prefitErr[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].prefitErr[k]);
				else if (psr[p].param[i].prefitErr[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].prefitErr[k]);
				}
				}
			fprintf(fout2,"\n");
			}
			}
		}
		if (psr[p].tempo1 == 1)
		{
		if (!strcmp(psr[p].clock, "TT(UTC(NIST))")) fprintf(fout2,"%-15.15s%s\n","CLK","UTC(NIST)");
		else if (!strcmp(psr[p].clock, "TT(TAI))")) fprintf(fout2,"%-15.15s%s\n","CLK","UTC(BIPM)");
		else if (!strcmp(psr[p].clock, "TT(UTC(PTB))")) fprintf(fout2,"%-15.15s%s\n","CLK","PTB");
		else if (!strcmp(psr[p].clock, "TT(TA(NIST))")) fprintf(fout2,"%-15.15s%s\n","CLK","AT1");
		else fprintf(fout2,"%-15.15s%s\n","CLK",psr[p].clock);
		}
		else	   
		fprintf(fout2,"%-15.15s%s\n","CLK",psr[p].clock);
	
		if (psr[p].fitMode==1) fprintf(fout2,"MODE 1\n");
		if (psr[p].units != SI_UNITS)
		fprintf(fout2, "%-15.15s%s\n", "UNITS", "TDB");
		if (psr[p].timeEphemeris != IF99_TIMEEPH)
		fprintf(fout2, "%-15.15s%s\n", "TIMEEPH", "FB90");
		if (!psr[p].dilateFreq)
		fprintf(fout2, "%-15.15s%s\n", "DILATEFREQ", "N");
		if (!psr[p].planetShapiro)
		fprintf(fout2, "%-15.15s%s\n", "PLANET_SHAPIRO", "N");
		if (psr[p].t2cMethod != T2C_IAU2000B)
		fprintf(fout2, "%-15.15s%s\n", "T2CMETHOD", "TEMPO");
		if (!psr[p].correctTroposphere)
		fprintf(fout2, "%-21.21s%s\n", "CORRECT_TROPOSPHERE", "N");
	
		fprintf(fout2,"%-15.15s%s\n","EPHEM",psr[p].ephemeris);
		fprintf(fout2,"%-15.15s%s\n","NITS","1");
		fprintf(fout2,"%-15.15s%d\n","NTOA",psr[p].nFit);
		fprintf(fout2,"%-15.15s%.4f %d\n","CHI2R",chisqr/(double)psr[p].fitNfree,psr[p].fitNfree);
		/*	  if (psr[p].tempo1 == 1)
		fprintf(fout2,"EPHVER         2\n");
		else
		fprintf(fout2,"EPHVER         5\n"); */
	
		/* Add jumps */
		for (i=1;i<=psr[p].nJumps;i++)
		{
		nread = sscanf(psr[p].jumpStr[i],"%s %s %s %s %s",str1,str2,str3,str4,str5);
		if (strcasecmp(str1,"FREQ")==0 || strcasecmp(str1,"MJD")==0)
			fprintf(fout2,"JUMP %s %s %s %.14g %d\n",str1,str2,str3,psr[p].jumpVal[i],psr[p].fitJump[i]);
		else if (strcasecmp(str1,"NAME")==0 || strcasecmp(str1,"TEL")==0 || str1[0]=='-')
			fprintf(fout2,"JUMP %s %s %.14g %d\n",str1,str2,psr[p].jumpVal[i],psr[p].fitJump[i]);
		}	  
	
		/* Add whitening flags */
		if (psr[p].param[param_wave_om].paramSet[0]==1)
		{
		fprintf(fout2,"WAVEEPOCH %.14Lg\n",psr[p].param[param_waveepoch].val[0]);
		fprintf(fout2,"WAVE_OM %.14Lg 0\n",psr[p].param[param_wave_om].val[0]);
		if (psr[p].waveScale!=0) fprintf(fout2,"WAVE_SCALE %g\n",psr[p].waveScale);
		for (i=0;i<psr[p].nWhite;i++)
			fprintf(fout2,"WAVE%d %.14g %.14g\n",i+1,psr[p].wave_sine[i],psr[p].wave_cos[i]);
		}
	
		if (psr[p].param[param_ifunc].paramSet[0]==1)
		{
		fprintf(fout2,"SIFUNC %d %d\n",(int)psr[p].param[param_ifunc].val[0],
			(int)psr[p].param[param_ifunc].fitFlag[0]);
		for (i=0;i<psr[p].ifuncN;i++)
			fprintf(fout2,"IFUNC%d %.15f %.10f %.10f\n",i+1,psr[p].ifuncT[i],
				psr[p].ifuncV[i],psr[p].ifuncE[i]);
		}
		/* Add phase jumps */
		for (i=0;i<psr[p].nPhaseJump;i++)
		{
		if (psr[p].phaseJumpDir[i]!=0)
			fprintf(fout2,"PHASE %+d %.14g\n",psr[p].phaseJumpDir[i],(double)(psr[p].obsn[psr[p].phaseJumpID[i]].sat+1.0/SECDAY));
		}
		// Add DM value parameters
		if (psr[p].param[param_dmmodel].paramSet[0]==1)
		{
	
		if (psr->param[param_dmmodel].linkTo[0] == param_dm){
			fprintf(fout2,"DMMODEL DM %d\n",(int)psr[p].param[param_dmmodel].fitFlag[0]);
		} else {
			fprintf(fout2,"DMMODEL %.14Lg %d\n",psr[p].param[param_dmmodel].val[0],(int)psr[p].param[param_dmmodel].fitFlag[0]);
		}
			bool useDMOFF=psr[p].dmoffsDMnum==psr[p].dmoffsCMnum;
			if(useDMOFF)for (i=0;i<psr[p].dmoffsDMnum;i++){
				if(psr[p].dmoffsDM_mjd[i]!=psr[p].dmoffsCM_mjd[i])useDMOFF=false;
				break;
			}
			if(useDMOFF){
			for (i=0;i<psr[p].dmoffsDMnum;i++)
				fprintf(fout2,"DMOFF\t %.15g %.15g %.15g\n",psr[p].dmoffsDM_mjd[i],psr[p].dmoffsDM[i],psr[p].dmoffsDM_error[i]);
	
			}else{
				for (i=0;i<psr[p].dmoffsDMnum;i++)
					fprintf(fout2,"_DM\t %.15g %.15g %.15g\n",psr[p].dmoffsDM_mjd[i],psr[p].dmoffsDM[i],psr[p].dmoffsDM_error[i]);
				for (i=0;i<psr[p].dmoffsCMnum;i++)
					fprintf(fout2,"_CM\t %.15g %.15g %.15g\n",psr[p].dmoffsCM_mjd[i],psr[p].dmoffsCM[i],psr[p].dmoffsCM_error[i]);
			}
			}
	
		// add constraints
		for (int i = 0; i < psr[p].nconstraints; i++){
			if (psr[p].constraints[i]==constraint_dmmodel_mean){
				fprintf(fout2,"CONSTRAIN DMMODEL\n");
			}
		}
		fclose(fout2);	 
		}
	}
}
	

void FindMLMasterFunc(int nParameters, double *startvec, double *endvec, double *stepvec, int *FitList, void *context, FILE *History) {

	printf("\n\n Entered Master ML Function \n\n");
	fprintf(History, "\n\n Entered Master ML Function \n\n");

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



	printf("Processing FitList...\n");
	fprintf(History, "Processing FitList...\n");

	int HaveMargin=0;
	for(int i =0; i < nParameters; i++){
		if(FitList[i]==0){
			printf("Fitting parameter %i \n",i);
			fprintf(History, "Fitting parameter %i \n",i);
		}
		if(FitList[i]==1){
			printf("Marginalising over parameter %i \n",i);
			fprintf(History, "Marginalising over parameter %i \n",i);
			HaveMargin++;
		}
		if(FitList[i]==2){
			printf("Fixing parameter %i \n",i);
			fprintf(History, "Fixing parameter %i \n",i);
		}
	}

	if(HaveMargin > 0){
		
		printf("Making new G Matrix for marginalisation\n");
		fprintf(History, "Making new G Matrix for marginalisation\n");
		if(FitList[0] != 1){
			FitList[0] = 1;
			HaveMargin++;
		}
		int Gsize=((MNStruct *)context)->pulse->nobs-HaveMargin;
		double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			TNDM[i]=new double[HaveMargin];
		}
	
		double **TNGM=new double*[((MNStruct *)context)->pulse->nobs];
		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			TNGM[i]=new double[((MNStruct *)context)->pulse->nobs];
		}
	
		
		//Get DMatrix for marginalisation
		getnewMarginDMatrix(((MNStruct *)context)->pulse, FitList, TNDM, ((MNStruct *)context)->TempoFitNums, ((MNStruct *)context)->TempoJumpNums,((MNStruct *)context)->numFitTiming,((MNStruct *)context)->numFitJumps);
		
		makeGDesign(((MNStruct *)context)->pulse, Gsize, HaveMargin, TNGM, TNDM);
		
		//Update context
		((MNStruct *)context)->Gsize=Gsize;
		((MNStruct *)context)->GMatrix=TNGM;
	}


	int pcount=0;
	// Obtain the starting point as rough parameter estimates
	for(int i=0; i<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++) {
		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		
		pcount++;
	} // for i

	for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){
		//printf("Setting EFAC: %i\n",i);
		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		
		pcount++;
	}

  	for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
		//printf("Setting EQUD: %i\n",i);
		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		
		pcount++;
	}
	if(((MNStruct *)context)->incRED == 1){
		//printf("Setting Red\n");
		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		pcount++;

		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		pcount++;
	}

	if(((MNStruct *)context)->incRED == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
			gsl_vector_set(vStart, pcount, startvec[pcount]);
			if(FitList[pcount] == 0){
				gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
			}
			else if(FitList[pcount] != 0){
				gsl_vector_set(vStepSize, pcount, 0.0);
			}
			
			pcount++;
		}
	}	
	if(((MNStruct *)context)->incDM  == 1){
		//printf("Setting Red\n");
		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		pcount++;

		gsl_vector_set(vStart, pcount, startvec[pcount]);
		if(FitList[pcount] == 0){
			gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
		}
		else if(FitList[pcount] != 0){
			gsl_vector_set(vStepSize, pcount, 0.0);
		}
		pcount++;
	}

	if(((MNStruct *)context)->incDM == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
			gsl_vector_set(vStart, pcount, startvec[pcount]);
			if(FitList[pcount] == 0){
				gsl_vector_set(vStepSize, pcount, stepvec[pcount]);
			}
			else if(FitList[pcount] != 0){
				gsl_vector_set(vStepSize, pcount, 0.0);
			}
			
			pcount++;
		}
	}	
  // Initial GSL vector of the vertex sizes, and the starting point



  // Initialise the iteration procedure
	if(HaveMargin == 0){
  		MinexFunc.f = &negloglikelihood;
	}
	if(HaveMargin > 0){
		MinexFunc.f = &neglogMarginlikelihood;
	}

  MinexFunc.n = nParameters;
  MinexFunc.params = context;
  pmMinimiser = gsl_multimin_fminimizer_alloc(pmtMinimiserType, nParameters);
  gsl_multimin_fminimizer_set(pmMinimiser, &MinexFunc, vStart, vStepSize);

	printf("\n\n Performing Minimisation over current model parameters \n\n");
	fprintf(History, "\n\n Performing Minimisation over current model parameters \n\n");
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
      fprintf(stderr, "Converged to maximum likelihood with downhill simplex: Convergence: %g Current Minimum: %g \n", nIteration,dSize,gsl_multimin_fminimizer_minimum(pmMinimiser));
     fprintf(History, "Converged to maximum likelihood with downhill simplex: Convergence: %g Current Minimum: %g \n", nIteration,dSize,gsl_multimin_fminimizer_minimum(pmMinimiser));
    } // if nStatus

    // Print iteration values
	if(nIteration == 1){
    		printf("First Step Likelihood: %g \n", gsl_multimin_fminimizer_minimum(pmMinimiser));
		fprintf(History, "First Step Likelihood: %g \n", gsl_multimin_fminimizer_minimum(pmMinimiser));
	}	
	if(nIteration % 100 ==0){
    		printf("Step[%i]: Convergence: %g Current Minimum: %g \n", nIteration,dSize,gsl_multimin_fminimizer_minimum(pmMinimiser));
		fprintf(History, "Step[%i]: Convergence: %g Current Minimum: %g \n", nIteration,dSize,gsl_multimin_fminimizer_minimum(pmMinimiser));
	}
  } while (nStatus == GSL_CONTINUE);


	printf("\n");

	
	for(int i =0; i < nParameters; i++){
		endvec[i] = pdParameterEstimates[i];
	}


	printf("\n");

	pcount=0;
	
	if(FitList[0]==0){
		printf("ML Phase: %.20Lg \n", endvec[0]*(((MNStruct *)context)->LDpriors[0][1])+(((MNStruct *)context)->LDpriors[0][0]));
		fprintf(History, "ML Phase: %.20Lg \n", endvec[0]*(((MNStruct *)context)->LDpriors[0][1])+(((MNStruct *)context)->LDpriors[0][0]));
	}
	if(FitList[0]==1){
		printf("Marginalised Phase at: %.20Lg \n", startvec[0]*(((MNStruct *)context)->LDpriors[0][1])+(((MNStruct *)context)->LDpriors[0][0]));
		fprintf(History, "Marginalised Phase at: %.20Lg \n", startvec[0]*(((MNStruct *)context)->LDpriors[0][1])+(((MNStruct *)context)->LDpriors[0][0]));
	}
	if(FitList[0]==2){
		printf("Fixed Phase at: %.20Lg \n", startvec[0]*(((MNStruct *)context)->LDpriors[0][1])+(((MNStruct *)context)->LDpriors[0][0]));
		fprintf(History, "Fixed Phase at: %.20Lg \n", startvec[0]*(((MNStruct *)context)->LDpriors[0][1])+(((MNStruct *)context)->LDpriors[0][0]));
	}		
	pcount++;

	for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){

		if(FitList[pcount]==0){
			printf("ML %s: %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], endvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
			fprintf(History, "ML %s: %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], endvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
		}
		if(FitList[pcount]==1){
			printf("Marginalised %s at: %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
			fprintf(History, "Marginalised %s at: %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
		}
		if(FitList[pcount]==2){
			printf("Fixed %s at: %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
			fprintf(History, "Fixed %s at: %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
		}

		pcount++;
	}

	for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){

		if(FitList[pcount]==0){
			printf("ML Jump %i: %.20Lg \n", j, endvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
			fprintf(History, "ML Jump %i: %.20Lg \n", j, endvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
		}
		if(FitList[pcount]==1){
			printf("Marginalised Jump %i at: %.20Lg \n", j, startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
			fprintf(History, "Marginalised Jump %i at: %.20Lg \n", j, startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
		}
		if(FitList[pcount]==2){
			printf("Fixed Jump %i at: %.20Lg \n", j, startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
			fprintf(History, "Fixed Jump %i at: %.20Lg \n", j, startvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]));
		}

		pcount++;

	}
	
	for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){

		if(FitList[pcount]==0){
			printf("ML EFAC %i: %.10g \n", i, endvec[pcount]);
			fprintf(History, "ML EFAC %i: %.10g \n", i, endvec[pcount]);
		}
		if(FitList[pcount]==1){
			printf("Marginalised EFAC %i, shouldn't have done this\n");
			fprintf(History, "Marginalised EFAC %i, shouldn't have done this\n");
		}
		if(FitList[pcount]==2){
			printf("Fixed EFAC %i at: %.10g \n", i, startvec[pcount]);
			fprintf(History, "Fixed EFAC %i at: %.10g \n", i, startvec[pcount]);
		}

		pcount++;
	}
 	 for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
		
		if(FitList[pcount]==0){
			printf("ML EQUAD %i: %.10g \n", i, endvec[pcount]);
			fprintf(History, "ML EQUAD %i: %.10g \n", i, endvec[pcount]);
		}
		if(FitList[pcount]==1){
			printf("Marginalised EQUAD %i, shouldn't have done this\n");
			fprintf(History, "Marginalised EQUAD %i, shouldn't have done this\n");
		}
		if(FitList[pcount]==2){
			printf("Fixed EQUAD %i at: %.10g \n", i, startvec[pcount]);
			fprintf(History, "Fixed EQUAD %i at: %.10g \n", i, startvec[pcount]);
		}

		pcount++;
	}
	if(((MNStruct *)context)->incRED == 1){

		if(FitList[pcount]==0){
			printf("ML Red Amp: %.10g \n", endvec[pcount]);
			fprintf(History, "ML Red Amp: %.10g \n", endvec[pcount]);
		}
		if(FitList[pcount]==1){
			printf("Marginalised Red Amp, shouldn't have done this\n");
			fprintf(History, "Marginalised Red Amp, shouldn't have done this\n");
		}
		if(FitList[pcount]==2){
			printf("Fixed Red Amp at: %.10g \n", startvec[pcount]);
			fprintf(History, "Fixed Red Amp at: %.10g \n", startvec[pcount]);
		}
		pcount++;


		if(FitList[pcount]==0){
			printf("ML Red Index: %.10g \n", endvec[pcount]);
			fprintf(History, "ML Red Index: %.10g \n", endvec[pcount]);
		}
		if(FitList[pcount]==1){
			printf("Marginalised Red Index, shouldn't have done this\n");
			fprintf(History, "Marginalised Red Index, shouldn't have done this\n");
		}
		if(FitList[pcount]==2){
			printf("Fixed Red Index at: %.10g \n", startvec[pcount]);
			fprintf(History, "Fixed Red Index at: %.10g \n", startvec[pcount]);
		}
		pcount++;
	}

	if(((MNStruct *)context)->incRED == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){

			if(FitList[pcount]==0){
				printf("ML Red Coeff %i: %.10g \n", i, endvec[pcount]);
				fprintf(History, "ML Red Coeff %i: %.10g \n", i, endvec[pcount]);
			}
			if(FitList[pcount]==1){
				printf("Marginalised Red Coeff %i, shouldn't have done this\n");
				fprintf(History, "Marginalised Red Coeff %i, shouldn't have done this\n");
			}
			if(FitList[pcount]==2){
				printf("Fixed Red Coeff %i at: %.10g \n", i, startvec[pcount]);
				fprintf(History, "Fixed Red Coeff %i at: %.10g \n", i, startvec[pcount]);
			}
			pcount++;
		}
	}	

	if(((MNStruct *)context)->incDM == 1){

		if(FitList[pcount]==0){
			printf("ML DM Amp: %.10g \n", endvec[pcount]);
			fprintf(History, "ML DM Amp: %.10g \n", endvec[pcount]);
		}
		if(FitList[pcount]==1){
			printf("Marginalised DM Amp, shouldn't have done this\n");
			fprintf(History, "Marginalised DM Amp, shouldn't have done this\n");
		}
		if(FitList[pcount]==2){
			printf("Fixed DM Amp at: %.10g \n", startvec[pcount]);
			fprintf(History, "Fixed DM Amp at: %.10g \n", startvec[pcount]);
		}
		pcount++;


		if(FitList[pcount]==0){
			printf("ML DM Index: %.10g \n", endvec[pcount]);	
			fprintf(History, "ML DM Index: %.10g \n", endvec[pcount]);
		}
		if(FitList[pcount]==1){
			printf("Marginalised DM Index, shouldn't have done this\n");
			fprintf(History, "Marginalised DM Index, shouldn't have done this\n");
		}
		if(FitList[pcount]==2){
			printf("Fixed DM Index at: %.10g \n", startvec[pcount]);
			fprintf(History, "Fixed DM Index at: %.10g \n", startvec[pcount]);
		}
		pcount++;
	}

	if(((MNStruct *)context)->incDM == 2){
		for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){

			if(FitList[pcount]==0){
				printf("ML DM Coeff %i: %.10g \n", i, endvec[pcount]);
				fprintf(History, "ML DM Coeff %i: %.10g \n", i, endvec[pcount]);
			}
			if(FitList[pcount]==1){
				printf("Marginalised DM Coeff %i, shouldn't have done this\n");
				fprintf(History, "Marginalised DM Coeff %i, shouldn't have done this\n");
			}
			if(FitList[pcount]==2){
				printf("Fixed DM Coeff %i at: %.10g \n", i, startvec[pcount]);
				fprintf(History, "Fixed DM Coeff %i at: %.10g \n", i, startvec[pcount]);
			}
			pcount++;
		}
	}


  gsl_vector_free(vStart);
  gsl_vector_free(vStepSize);
  gsl_multimin_fminimizer_free(pmMinimiser);
  delete[] pdParameterEstimates;
} 


void FindMLHypervisor(int nParameters, void *context, std::string longname){

	int *FitList=new int[nParameters];
	double *startvec = new double[nParameters];
	double *stepvec = new double[nParameters];
	double *endvec = new double[nParameters];

	FILE *pFile;
	std::string MLHistoryname = longname+"MLHistory.txt";
	pFile = fopen(MLHistoryname.c_str(), "w+");

     
     
	int pcount=0;
	printf("\n\n Entered ML Hypervisor Function\n\n");
	fprintf(pFile,"\n\n Entered ML Hypervisor Function\n\n");

	for(int it=0;it<2;it++){

		printf("\n Beginning Iteration %i of %i\n",it+1,2);
		fprintf(pFile, "\n Beginning Iteration %i of %i\n",it+1,2);


		if(((MNStruct *)context)->numFitEFAC > 0 || ((MNStruct *)context)->numFitEQUAD > 0 || ((MNStruct *)context)->incRED > 0 || ((MNStruct *)context)->incDM >0){
	
			printf("\n Stochastic Parameters Detected \n");
			fprintf(pFile,"\n Stochastic Parameters Detected \n");	
			
			printf("\n\n Marginalising over Current Timing Model, Finding ML Stochastic Parameters \n\n");
			fprintf(pFile,"\n\n Marginalising over Current Timing Model, Finding ML Stochastic Parameters \n\n");
		
			pcount=0;
			for(int i=0; i<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++) {
				startvec[pcount]=0;
				stepvec[pcount]=10.0;
				FitList[pcount]=1;
				pcount++;
			} 
			
			for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){
		
				startvec[pcount]=2.0;
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
			}
			for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
		
				startvec[pcount]=-5;
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
			}
			if(((MNStruct *)context)->incRED == 1){
				//printf("Setting Red\n");
				startvec[pcount]=-12.0;
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
				startvec[pcount]=3.1;
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
			}
		
			if(((MNStruct *)context)->incRED == 2){
				for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
					startvec[pcount]=-3.0;
					stepvec[pcount]=0.5;
					FitList[pcount]=0;
					pcount++;
				}
			}
		
			if(((MNStruct *)context)->incDM == 1){
				//printf("Setting Red\n");
				startvec[pcount]=-10.0;
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
				startvec[pcount]=3.1;
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
			}
		
			if(((MNStruct *)context)->incDM == 2){
				for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
					startvec[pcount]=0.0;
					stepvec[pcount]=0.5;
					FitList[pcount]=0;
					pcount++;
				}
			}	
		
			FindMLMasterFunc(nParameters, startvec, endvec,  stepvec,  FitList,  context, pFile);
	
		}
		else if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD == 0 && ((MNStruct *)context)->incRED == 0 && ((MNStruct *)context)->incDM ==0){
	
			for(int i =0; i < nParameters; i++){
				endvec[i]=startvec[i];
			}
		}
	
		printf("\n\n Fixing Any Stochastic Parameters, Marginalising over Jumps, Finding ML Timing model \n\n");
		fprintf(pFile, "\n\n Fixing Any Stochastic Parameters, Marginalising over Jumps, Finding ML Timing model \n\n");
		pcount=0;
		for(int i=0; i<((MNStruct *)context)->numFitTiming; i++) {
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=10.0;
			FitList[pcount]=0;
			pcount++;
		} 
	
		for(int i=0; i<((MNStruct *)context)->numFitJumps; i++) {
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=10.0;
			FitList[pcount]=1;
			pcount++;
		} 
		
		for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){
	
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=2;
			pcount++;
		}
		for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
	
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=2;
			pcount++;
		}
		if(((MNStruct *)context)->incRED == 1){
			//printf("Setting Red\n");
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=2;
			pcount++;
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=2;
			pcount++;
		}
	
		if(((MNStruct *)context)->incRED == 2){
			for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
				startvec[pcount]=endvec[pcount];
				stepvec[pcount]=0.5;
				FitList[pcount]=2;
				pcount++;
			}
		}
	
		if(((MNStruct *)context)->incDM == 1){
			//printf("Setting Red\n");
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=2;
			pcount++;
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=2;
			pcount++;
		}
	
		if(((MNStruct *)context)->incDM == 2){
			for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
				startvec[pcount]=endvec[pcount];
				stepvec[pcount]=0.5;
				FitList[pcount]=2;
				pcount++;
			}
		}
	
		FindMLMasterFunc(nParameters, startvec, endvec,  stepvec,  FitList,  context, pFile);
	
	
		printf("\n\n Marginalising over Jumps, Finding ML Timing model and Stochastic Parameters \n\n");
		fprintf(pFile, "\n\n Marginalising over Jumps, Finding ML Timing model and Stochastic Parameters \n\n");
		pcount=0;
		for(int i=0; i<((MNStruct *)context)->numFitTiming; i++) {
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=10.0;
			FitList[pcount]=0;
			pcount++;
		} 
	
		for(int i=0; i<((MNStruct *)context)->numFitJumps; i++) {
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=10.0;
			FitList[pcount]=1;
			pcount++;
		} 
		
		for(int i =0; i< ((MNStruct *)context)->numFitEFAC; i++){
	
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=0;
			pcount++;
		}
		for(int i =0; i< ((MNStruct *)context)->numFitEQUAD; i++){
	
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=0;
			pcount++;
		}
		if(((MNStruct *)context)->incRED == 1){
			//printf("Setting Red\n");
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=0;
			pcount++;
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=0;
			pcount++;
		}
	
		if(((MNStruct *)context)->incRED == 2){
			for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
				startvec[pcount]=endvec[pcount];
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
			}
		}
	
		if(((MNStruct *)context)->incDM == 1){
			//printf("Setting Red\n");
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=0;
			pcount++;
			startvec[pcount]=endvec[pcount];
			stepvec[pcount]=0.5;
			FitList[pcount]=0;
			pcount++;
		}
	
		if(((MNStruct *)context)->incDM == 2){
			for(int i =0; i< ((MNStruct *)context)->numFitRedCoeff; i++){
				startvec[pcount]=endvec[pcount];
				stepvec[pcount]=0.5;
				FitList[pcount]=0;
				pcount++;
			}
		}
	
		FindMLMasterFunc(nParameters, startvec, endvec,  stepvec,  FitList,  context, pFile);

		for(int i =0; i < nParameters; i++){
			startvec[i]=endvec[i];
		}
	}


	//Update Pulsar to new ML values

	pcount=1;
	printf("\nFinal ML Timing model: \n");	
	fprintf(pFile, "\nFinal ML Timing model: \n");	
	for(int i=1; i<((MNStruct *)context)->numFitTiming; i++) {
		long double value = endvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
		long double error = 0;
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].val[((MNStruct *)context)->TempoFitNums[pcount][1]] = value;
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].err[((MNStruct *)context)->TempoFitNums[pcount][1]] = error;
		printf("Max %s : %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], value);
		fprintf(pFile, "Max %s : %.20Lg \n", ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]], value);
		
		pcount++;
	} 

	for(int i=0; i<((MNStruct *)context)->numFitJumps; i++) {
		long double value = endvec[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
		long double error = 0;
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[i]] = value;
		((MNStruct *)context)->pulse->jumpValErr[((MNStruct *)context)->TempoJumpNums[i]] = error;
		printf("Max Jump %i :  %.10Lg \n", i,value);
		fprintf(pFile, "Max Jump %i :  %.10Lg \n", i,value);
		pcount++;
	} 

	printf("\n\n Finished ML analysis \n\n");
	fprintf( pFile, "\n\n Finished ML analysis \n\n");

	writeMLPar(longname, ((MNStruct *)context)->pulse);

	delete[] startvec;
	delete[] endvec;
	delete[] FitList;
	fclose(pFile);	

}
