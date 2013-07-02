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
#include <gsl/gsl_sf_gamma.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "dpotrs.h"
#include "tempo2.h"
#include "TempoNest.h"

void vHRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD, redamp, redalpha, DMamp, DMalpha;
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

void vHRedMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD, redamp, redalpha, DMamp, DMalpha;
	int pcount=0;
	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
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


void WhiteLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
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
	
	//printf("White: %g %g \n", detN, Chisq);

	if(isnan(detN) || isinf(detN) || isnan(Chisq) || isinf(Chisq)){

		lnew=-pow(10.0,200);
	}
	else{
		lnew = -0.5*(((MNStruct *)context)->pulse->nobs*log(2*M_PI) + detN + Chisq);	
	}

	delete[] EFAC;
	delete[] Resvec;


}


void WhiteMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;

	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;

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


void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	double phase=0;
	double *EFAC;
	double EQUAD;
	int pcount=0;
	int bad=0;
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

		phase=(double)LDparams[0];
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
		//printf("Phase: %g \n", phase);
	
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



	int FitCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
        int totCoeff=0;
        if(((MNStruct *)context)->incRED != 0)totCoeff+=FitCoeff;
        if(((MNStruct *)context)->incDM != 0)totCoeff+=FitCoeff;

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[totCoeff];

	double tdet=0;
	double timelike=0;
	double timelike2=0;

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
			timelike2=timelike2+pow((double)((MNStruct *)context)->pulse->obsn[o].residual,2)*WorkNoise[o];
			//printf("CPUTlike: %i %g %g %g\n",o,timelike,timelike2,phase);
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


	    int info=0;
        dpotrfInfo(PPFM, totCoeff, jointdet,info);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
        }
	
	lnew=-0.5*(jointdet+tdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
		
	}
	
	//printf("CPULike: %g %g %g %g %g %g \n", lnew, jointdet, tdet, freqdet, timelike, freqlike);
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

void LRedMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double EQUAD;
	int pcount=0;
	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	if(((MNStruct *)context)->incRED==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
	if(((MNStruct *)context)->incDM==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
        if(((MNStruct *)context)->incRED==3)totdims+=2;
        if(((MNStruct *)context)->incDM==3)totdims+=2;
	int numfit=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	int fitcount=0;
	for(int p=0;p<totdims;p++){
	//	printf("DP: %i %g %g \n", p, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1]);
		if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){
		Cube[pcount]=(((MNStruct *)context)->Dpriors[p][1]-((MNStruct *)context)->Dpriors[p][0])*Cube[pcount]+((MNStruct *)context)->Dpriors[p][0];
	//	printf("Cube: %i %g %g %g \n", pcount, ((MNStruct *)context)->Dpriors[p][0], ((MNStruct *)context)->Dpriors[p][1], Cube[pcount]);
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
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

}

