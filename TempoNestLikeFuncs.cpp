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

	double *EFAC, *EQUAD;
	double redamp, redalpha, DMamp, DMalpha;
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
// 			printf("step: %g %g \n", StepAmp, StepTime);
			for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
				double time=((MNStruct *)context)->pulse->obsn[o1].bat;
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
	double noiseval=0;	
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
                if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
			if(((MNStruct *)context)->whitemodel == 0){
				noiseval=pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
				noiseval= EFAC[0]*EFAC[0]*(pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6)),2) + EQUAD[0]);
			}
                }
                else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
			if(((MNStruct *)context)->whitemodel == 0){
                       	 	noiseval=pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD[((MNStruct *)context)->sysFlags[o1]];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
                       	 	noiseval=EFAC[((MNStruct *)context)->sysFlags[o1]]*EFAC[((MNStruct *)context)->sysFlags[o1]]*(pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o1]]);
			}
                }

		CovMatrix[o1][o1] += noiseval;
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


	double noiseval=0;	
	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
                if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
			if(((MNStruct *)context)->whitemodel == 0){
				noiseval=pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
				noiseval= EFAC[0]*EFAC[0]*(pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6)),2) + EQUAD[0]);
			}
                }
                else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
			if(((MNStruct *)context)->whitemodel == 0){
                       	 	noiseval=pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD[((MNStruct *)context)->sysFlags[o1]];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
                       	 	noiseval=EFAC[((MNStruct *)context)->sysFlags[o1]]*EFAC[((MNStruct *)context)->sysFlags[o1]]*(pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o1]]);
			}
                }

		CovMatrix[o1][o1] += noiseval;
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
	double *EQUAD;
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
		
//		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
//		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
//		
//		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
//			Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;
//		}
	
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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
			//printf("FitVec: %i %g \n",o,Fitvec[o]);
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





	double Chisq=0;
	double noiseval=0;
	double detN=0;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

                if(((MNStruct *)context)->numFitEQUAD < 2 && ((MNStruct *)context)->numFitEFAC < 2){
			if(((MNStruct *)context)->whitemodel == 0){
				noiseval=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
				noiseval= EFAC[0]*EFAC[0]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[0]);
			}
                }
                else if(((MNStruct *)context)->numFitEQUAD > 1 || ((MNStruct *)context)->numFitEFAC > 1){
			if(((MNStruct *)context)->whitemodel == 0){
                       	 	noiseval=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
			}
			else if(((MNStruct *)context)->whitemodel == 1){
                       	 	noiseval=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
			}
                }
	
	

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
	if(((MNStruct *)context)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)context)->incStep; i++){
			delete[] Steps[i];
		}
		delete[] Steps;
	}


}


void WhiteMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
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
// 	printf("pcount: %i \n",pcount);
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




	double *Noise;
	double *GRes=new double[((MNStruct *)context)->Gsize];
	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[((MNStruct *)context)->Gsize];
	double** GG = new double*[((MNStruct *)context)->Gsize]; for (int k=0; k<((MNStruct *)context)->Gsize; k++) GG[k] = new double[((MNStruct *)context)->Gsize];
	double *WorkGRes = new double[((MNStruct *)context)->Gsize];
		
	double Chisq=0;
	double det=0;
	
	if(((MNStruct *)context)->numFitEFAC >1 || ((MNStruct *)context)->numFitEQUAD > 1){
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

		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			for(int j=0;j<((MNStruct *)context)->Gsize; j++){
				NG[i][j]=((MNStruct *)context)->GMatrix[i][j]*Noise[i];
			}
		}



		dgemm(((MNStruct *)context)->GMatrix, NG,GG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'T','N');
	

		dgemv(((MNStruct *)context)->GMatrix,Resvec,GRes,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'T');



	
	
		for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
			WorkGRes[o1]=GRes[o1];
		}
		dpotrf(GG, ((MNStruct *)context)->Gsize, det);
		dpotrs(GG, WorkGRes, ((MNStruct *)context)->Gsize);

		

	
		for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
			Chisq += GRes[o1]*WorkGRes[o1];
		}
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

		double *UGRes=new double[((MNStruct *)context)->Gsize];
		dgemv(((MNStruct *)context)->UMatrix,Resvec,UGRes,((MNStruct *)context)->Gsize,((MNStruct *)context)->pulse->nobs,'N');
		

		for(int o1=0;o1<((MNStruct *)context)->Gsize; o1++){
			Chisq += UGRes[o1]*UGRes[o1]*Noise[o1];

		}
		delete[] UGRes;

	
	}
	
	else if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD == 0){
		Noise=new double[((MNStruct *)context)->Gsize];
		det=((MNStruct *)context)->staticTimeDet;
		double *GVec=new double[((MNStruct *)context)->pulse->nobs];
		dgemv(((MNStruct *)context)->staticGMatrix,Resvec,GVec,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
		for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){
			Chisq += GVec[o1]*Resvec[o1];
		}
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

	


}


void LRedLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;
	double phase=0;
	double *EFAC;
	double *EQUAD;
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
// 		printf("Phase: %g \n", phase);
	
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
	if(((MNStruct *)context)->incFloatDM != 0)FitDMCoeff+=2*((MNStruct *)context)->incFloatDM;
	if(((MNStruct *)context)->incFloatRed != 0)FitRedCoeff+=2*((MNStruct *)context)->incFloatRed;
    	int totCoeff=0;
    	if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
    	if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[totCoeff];

	double tdet=0;
	double timelike=0;
	double timelike2=0;

	if(((MNStruct *)context)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			
			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
			timelike2=timelike2+pow((double)((MNStruct *)context)->pulse->obsn[o].residual,2)*WorkNoise[o];

		}
	}
	else if(((MNStruct *)context)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			
			WorkNoise[o]=EFAC[((MNStruct *)context)->sysFlags[o]]*EFAC[((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]);
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
			timelike2=timelike2+pow((double)((MNStruct *)context)->pulse->obsn[o].residual,2)*WorkNoise[o];

		}
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
                for (int i=0; i<FitRedCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[i+FitRedCoeff/2]=freqs[i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }


                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }



                startpos=FitRedCoeff;

        }
   else if(((MNStruct *)context)->incRED==3){

                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;
// 		printf("red: %g %g \n", redamp, redindex);
                 redamp=pow(10.0, redamp);

                freqdet=0;
                 for (int i=0; i<FitRedCoeff/2; i++){

                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);


                 }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }


                startpos=FitRedCoeff;

        }


       if(((MNStruct *)context)->incDM==2){

                for (int i=0; i<FitDMCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }



        }
        else if(((MNStruct *)context)->incDM==3){
                double DMamp=Cube[pcount];
                pcount++;
                double DMindex=Cube[pcount];
                pcount++;

                DMamp=pow(10.0, DMamp);

                 for (int i=0; i<FitDMCoeff/2; i++){
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);


                 }
                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
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
	delete[] EQUAD;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] NFd;
	delete[] freqs;

	for (int j = 0; j < totCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] NF[j];
	}
	delete[] NF;

	for (int j = 0; j < totCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[] FMatrix[j];
	}
	delete[] FMatrix;

}


void LRedMarginLogLike2(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;
	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;
	if(((MNStruct *)context)->incRED==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
	if(((MNStruct *)context)->incDM==2)totdims+=((MNStruct *)context)->numFitDMCoeff;
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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}
	pcount=fitcount;	

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

	EFAC[0]=0.113838908428069896E+01;
	EQUAD[0]=pow(10.0,2*-0.738979250214491135E+01);
	

	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);
        int totCoeff=0;
        if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
        if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

	double *powercoeff=new double[totCoeff];


	double *Noise;
	double *GRes=new double[((MNStruct *)context)->Gsize];
	
	double tdet=0;
	double timelike=0;
	double *dG=new double[((MNStruct *)context)->pulse->nobs];
	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[((MNStruct *)context)->Gsize];
	double** GG = new double*[((MNStruct *)context)->Gsize]; for (int k=0; k<((MNStruct *)context)->Gsize; k++) GG[k] = new double[((MNStruct *)context)->Gsize];
	double **GNG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) GNG[k] = new double[((MNStruct *)context)->pulse->nobs];	
	
	if(((MNStruct *)context)->numFitEFAC >1 || ((MNStruct *)context)->numFitEQUAD > 1){
		Noise=new double[((MNStruct *)context)->pulse->nobs];
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			Noise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD[((MNStruct *)context)->sysFlags[o]];
		}

		
		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			for(int j=0;j<((MNStruct *)context)->Gsize; j++){

				NG[i][j]=((MNStruct *)context)->GMatrix[i][j]*Noise[i];

			}
		}

		

		dgemm(((MNStruct *)context)->GMatrix, NG,GG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'T','N');

	
		tdet=0;
		dpotrf(GG, ((MNStruct *)context)->Gsize, tdet);
		dpotri(GG,((MNStruct *)context)->Gsize);
		//printf("det: %g\n", tdet);


		dgemm(((MNStruct *)context)->GMatrix, GG,NG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, 'N','N');

		

		dgemm(NG, ((MNStruct *)context)->GMatrix, GNG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'N','T');
		
		
		dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
		timelike=0;
		for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
			timelike+=Resvec[o1]*dG[o1];
		//	printf("timlike: %i %g %g \n", o1, Resvec[o1], dG[o1]);
		}
		
	}

        else if(((MNStruct *)context)->numFitEFAC == 1 || ((MNStruct *)context)->numFitEQUAD==1 && ((MNStruct *)context)->numFitEFAC < 2 && ((MNStruct *)context)->numFitEQUAD < 2){
		Noise=new double[((MNStruct *)context)->Gsize];
                for(int o=0;o<((MNStruct *)context)->Gsize; o++){
                        Noise[o]=pow(EFAC[0],2)*((MNStruct *)context)->SVec[o] + EQUAD[0];
                      //  printf("tdet: %i %g %g %g \n",o, EFAC[0],   EQUAD[0], ((MNStruct *)context)->SVec[o] );
                        tdet += log(Noise[o]);
                        Noise[o] = 1.0/Noise[o];
                }
                
        		for(int i=0;i<((MNStruct *)context)->Gsize;i++){
					for(int j=0;j<((MNStruct *)context)->pulse->nobs; j++){

						NG[j][i]= ((MNStruct *)context)->UMatrix[i][j]*Noise[i];

					}
				}
		
				dgemm(NG, ((MNStruct *)context)->UMatrix, GNG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->pulse->nobs, 'N','N');
				
				

               dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
				timelike=0;
				for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
					timelike+=Resvec[o1]*dG[o1];
				}

        }

	else if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD == 0){
		Noise=new double[((MNStruct *)context)->pulse->nobs];
		tdet=((MNStruct *)context)->staticTimeDet;
		for (int j=0; j<((MNStruct *)context)->pulse->nobs; j++) {
			for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) {
				GNG[j][k] = ((MNStruct *)context)->staticGMatrix[j][k];
			}
		}
		
		dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
		timelike=0;
		for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
			timelike+=Resvec[o1]*dG[o1];
		}
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

	double maxtspan=1*(end-start);



	



       double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
        if(((MNStruct *)context)->incRED==2){
                for (int i=0; i<FitRedCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[i+FitRedCoeff/2]=freqs[i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }


                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }



                startpos=FitRedCoeff;

        }
   else if(((MNStruct *)context)->incRED==3){

                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;
	//	printf("red: %g %g \n", redamp, redindex);
	
				redamp = -0.380723655956649498E+01;
				redindex = 0.406554824919563018E+00;
	
                 redamp=pow(10.0, redamp);

                freqdet=0;
                 for (int i=0; i<FitRedCoeff/2; i++){

                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        
                        printf("Red: %g %g \n",log10(freqs[startpos+i]), log10(redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)));
                        
                        freqdet=freqdet+2*log(powercoeff[i]);


                 }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
		//		if(k==0)printf("FM %i %g %g \n", i, freqs[i], time);

                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }


                startpos=FitRedCoeff;

        }



 // Cube[pcount]=   -0.990748409075186132E-01;
  //Cube[pcount+1]=    -0.112404079564142467E+01;
  //Cube[pcount+2]=   -0.498555523093632935E+00;
  //Cube[pcount+3]=   -0.929347396735936471E+00;
  //Cube[pcount+4]=   -0.146576763218658490E+01;
  //Cube[pcount+5]=   -0.111278196533701568E+01;
  //Cube[pcount+6]=   -0.175820832419789097E+01;
  //Cube[pcount+7]=   -0.122412406111443017E+01;
  //Cube[pcount+8]=  -0.114232473771063381E+01;
  //Cube[pcount+9]=   -0.425980556384346443E+01;
  //Cube[pcount+10]=   -0.371279932392764245E+00;



       if(((MNStruct *)context)->incDM==2){

                for (int i=0; i<FitDMCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
			printf("DM PC: %i %g \n", pcount, pc);
                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }



        }
        else if(((MNStruct *)context)->incDM==3){
                double DMamp=Cube[pcount];
                pcount++;
                double DMindex=Cube[pcount];
                pcount++;
                
                
                DMamp= -0.532817247563341700E+00;
                DMindex=  0.134328061336200943E+01;
                
	//	printf("DM: %g %g \n", DMamp, DMindex);
                DMamp=pow(10.0, DMamp);

                 for (int i=0; i<FitDMCoeff/2; i++){
  						freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
                        powercoeff[startpos+i]=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        printf("DM: %g %g \n",log10(freqs[startpos+i]), log10(DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)));
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);


                 }
                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }


                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
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

	double *signal=new double[((MNStruct *)context)->pulse->nobs];
	double *Gsignal=new double[((MNStruct *)context)->Gsize];
	double *GGsignal=new double[((MNStruct *)context)->pulse->nobs];
	double *GGdata=new double[((MNStruct *)context)->pulse->nobs];
	double *maxvec=new double[totCoeff];


        dpotrf(PPFM, totCoeff, jointdet);
        
        
    	dpotri(PPFM,totCoeff);

	printf("start1\n");
 	dgemv(PPFM,NFd,maxvec,totCoeff,totCoeff,'T');
 	dgemv(FMatrix,maxvec,signal,((MNStruct *)context)->pulse->nobs,totCoeff,'N');
 	printf("start2\n");
 	dgemv(((MNStruct *)context)->GMatrix,signal,Gsignal,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'T');
 	dgemv(((MNStruct *)context)->GMatrix,Gsignal,GGsignal,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'N');
 	printf("start3\n");
 	dgemv(((MNStruct *)context)->GMatrix,Resvec,GRes,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'T');
 	dgemv(((MNStruct *)context)->GMatrix,GRes,GGdata,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->Gsize,'N');

 	for(int i=0; i< ((MNStruct *)context)->pulse->nobs; i++){
 		Noise[i]=sqrt(pow(((((MNStruct *)context)->pulse->obsn[i].toaErr)*pow(10.0,-6))*EFAC[0],2) + EQUAD[0]);
 		//printf("err %i %g %g %g \n",i,EFAC[((MNStruct *)context)->sysFlags[i]],EQUAD,Noise[i]);
 	}
 
 	
 	for(int i=0; i< ((MNStruct *)context)->pulse->nobs; i++){
 		printf("%g %g %g %g %g\n",(double)((MNStruct *)context)->pulse->obsn[i].bat, GGdata[i], GGsignal[i], Noise[i],(double)((MNStruct *)context)->pulse->obsn[i].freqSSB );
 	}
 	printf("finish\n");


        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
	//	printf("freqlike: %i %g %g %g\n",j, freqlike, NFd[j], WorkCoeff[j]);
        }
	
	lnew=-0.5*(((double)((MNStruct *)context)->Gsize)*log(2.0*M_PI) + tdet+jointdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}

	delete[] DMVec;
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] powercoeff;
	delete[] NFd;
	delete[] dG;
	delete[] freqs;

	for (int j = 0; j < totCoeff; j++){
		delete[]PPFM[j];
	}
	delete[]PPFM;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]NF[j];
	}
	delete[]NF;

	for (int j = 0; j < totCoeff; j++){
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
	
//	printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

}



void LRedMarginLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{

	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;
	int totdims=((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps + ((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD;
	totdims +=2*((MNStruct *)context)->incStep;
	if(((MNStruct *)context)->incRED==2)totdims+=((MNStruct *)context)->numFitRedCoeff;
	if(((MNStruct *)context)->incDM==2)totdims+=((MNStruct *)context)->numFitDMCoeff;
	if(((MNStruct *)context)->incFloatDM != 0)totdims+=2*((MNStruct *)context)->incFloatDM;
	if(((MNStruct *)context)->incFloatRed != 0)totdims+=2*((MNStruct *)context)->incFloatRed;
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

	double *Noise;
	double *GRes=new double[((MNStruct *)context)->Gsize];
	
	double tdet=0;
	double timelike=0;
	double *dG=new double[((MNStruct *)context)->pulse->nobs];
	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[((MNStruct *)context)->Gsize];
	double** GG = new double*[((MNStruct *)context)->Gsize]; for (int k=0; k<((MNStruct *)context)->Gsize; k++) GG[k] = new double[((MNStruct *)context)->Gsize];
	double **GNG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) GNG[k] = new double[((MNStruct *)context)->pulse->nobs];	
	
	if(((MNStruct *)context)->numFitEFAC >1 || ((MNStruct *)context)->numFitEQUAD > 1){
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


		
		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			for(int j=0;j<((MNStruct *)context)->Gsize; j++){

				NG[i][j]=((MNStruct *)context)->GMatrix[i][j]*Noise[i];

			}
		}

		

		dgemm(((MNStruct *)context)->GMatrix, NG,GG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'T','N');

	
		tdet=0;
		dpotrf(GG, ((MNStruct *)context)->Gsize, tdet);
		dpotri(GG,((MNStruct *)context)->Gsize);
		//printf("det: %g\n", tdet);


		dgemm(((MNStruct *)context)->GMatrix, GG,NG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, 'N','N');

		

		dgemm(NG, ((MNStruct *)context)->GMatrix, GNG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, 'N','T');
		
		
		dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
		timelike=0;
		for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
			timelike+=Resvec[o1]*dG[o1];
		//	printf("timlike: %i %g %g \n", o1, Resvec[o1], dG[o1]);
		}
		
	}

        else if(((MNStruct *)context)->numFitEFAC == 1 || ((MNStruct *)context)->numFitEQUAD==1 && ((MNStruct *)context)->numFitEFAC < 2 && ((MNStruct *)context)->numFitEQUAD < 2){
		Noise=new double[((MNStruct *)context)->Gsize];
		if(((MNStruct *)context)->whitemodel == 0){
			for(int o=0;o<((MNStruct *)context)->Gsize; o++){
				Noise[o]=pow(EFAC[0],2)*((MNStruct *)context)->SVec[o] + EQUAD[0];
				tdet += log(Noise[o]);
				Noise[o] = 1.0/Noise[o];
			}
		}
		if(((MNStruct *)context)->whitemodel == 1){
			for(int o=0;o<((MNStruct *)context)->Gsize; o++){
				Noise[o]=pow(EFAC[0],2)*(((MNStruct *)context)->SVec[o] + EQUAD[0]);
				tdet += log(Noise[o]);
				Noise[o] = 1.0/Noise[o];
			}
		}
                
        		for(int i=0;i<((MNStruct *)context)->Gsize;i++){
					for(int j=0;j<((MNStruct *)context)->pulse->nobs; j++){

						NG[j][i]= ((MNStruct *)context)->UMatrix[i][j]*Noise[i];

					}
				}
		
				dgemm(NG, ((MNStruct *)context)->UMatrix, GNG,((MNStruct *)context)->pulse->nobs, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->Gsize, ((MNStruct *)context)->pulse->nobs, 'N','N');
				
				

               dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
				timelike=0;
				for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
					timelike+=Resvec[o1]*dG[o1];
				}

        }

	else if(((MNStruct *)context)->numFitEFAC == 0 && ((MNStruct *)context)->numFitEQUAD == 0){
		Noise=new double[((MNStruct *)context)->pulse->nobs];
		tdet=((MNStruct *)context)->staticTimeDet;
		for (int j=0; j<((MNStruct *)context)->pulse->nobs; j++) {
			for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) {
				GNG[j][k] = ((MNStruct *)context)->staticGMatrix[j][k];
			}
		}
		
		dgemv(GNG,Resvec,dG,((MNStruct *)context)->pulse->nobs,((MNStruct *)context)->pulse->nobs,'T');
	
		timelike=0;
		for(int o1=0; o1<((MNStruct *)context)->pulse->nobs; o1++){
			timelike+=Resvec[o1]*dG[o1];
		}
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

	double maxtspan=1*(end-start);


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
//                 for (int i=0; i<FitRedCoeff/2; i++){
//                         int pnum=pcount;
//                         double pc=Cube[pcount];
//                         freqs[i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
//                         freqs[i+FitRedCoeff/2]=freqs[i];
// 
//                         powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
//                         powercoeff[i+FitRedCoeff/2]=powercoeff[i];
//                         freqdet=freqdet+2*log(powercoeff[i]);
//                         pcount++;
//                 }


                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }



                startpos=FitRedCoeff;

        }
   else if(((MNStruct *)context)->incRED==3){
/*
                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;
	//	printf("red: %g %g \n", redamp, redindex);
                 redamp=pow(10.0, redamp);

                freqdet=0;
                 for (int i=0; i<FitRedCoeff/2; i++){

                        freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);


                 }*/

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

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i]=cos(2*M_PI*freqs[i]*time);
		//		if(k==0)printf("FM %i %g %g \n", i, freqs[i], time);

                        }
                }

                for(int i=0;i<FitRedCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][i+FitRedCoeff/2]=sin(2*M_PI*freqs[i]*time);
                        }
                }


                startpos=FitRedCoeff;

        }


       if(((MNStruct *)context)->incDM==2){

//                 for (int i=0; i<FitDMCoeff/2; i++){
//                         int pnum=pcount;
//                         double pc=Cube[pcount];
//                         freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2+i]/maxtspan;
//                         freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
// 
//                         powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
//                         powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
//                         freqdet=freqdet+2*log(powercoeff[startpos+i]);
//                         pcount++;
//                 }
// 
//                 for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
//                         DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
//                 }

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

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }



        }
        else if(((MNStruct *)context)->incDM==3){
//                 double DMamp=Cube[pcount];
//                 pcount++;
//                 double DMindex=Cube[pcount];
//                 pcount++;
// 	//	printf("DM: %g %g \n", DMamp, DMindex);
//                 DMamp=pow(10.0, DMamp);
// 
//                  for (int i=0; i<FitDMCoeff/2; i++){
//   						freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2+i]/maxtspan;
//                         freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
//                         powercoeff[startpos+i]=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);///(365.25*24*60*60)/4;
//                         powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
//                         freqdet=freqdet+2*log(powercoeff[startpos+i]);
// 
// 
//                  }
//                 for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
//                         DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
//                 }

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

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
                        }
                }

                for(int i=0;i<FitDMCoeff/2;i++){
                        for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                                double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
                                FMatrix[k][startpos+i+FitDMCoeff/2]=sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
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
	
	lnew=-0.5*(((double)((MNStruct *)context)->Gsize)*log(2.0*M_PI) + tdet+jointdet+freqdet+timelike-freqlike);

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}

	delete[] DMVec;
	delete[] WorkCoeff;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] powercoeff;
	delete[] NFd;
	delete[] dG;
	delete[] freqs;

	for (int j = 0; j < totCoeff; j++){
		delete[]PPFM[j];
	}
	delete[]PPFM;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]NF[j];
	}
	delete[]NF;

	for (int j = 0; j < totCoeff; j++){
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
	
//	printf("CPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

}


void LRedNumLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
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
		
		if(((MNStruct *)context)->pulse->param[param_dmmodel].fitFlag[0] == 1){
			int DMnum=((MNStruct *)context)->pulse[0].dmoffsDMnum;
			for(int i =0; i < DMnum; i++){
				((MNStruct *)context)->pulse[0].dmoffsDM[i]=Cube[ndim-DMnum+i];
				pcount++;
			}
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
			Resvec[o]=((MNStruct *)context)->pulse->obsn[o].residual-Fitvec[o];
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
    int MarginCoeff=0;
    if(((MNStruct *)context)->incRED == 2 || ((MNStruct *)context)->incRED == 3)MarginCoeff+=FitCoeff;
    if(((MNStruct *)context)->incDM == 2 || ((MNStruct *)context)->incDM == 3)MarginCoeff+=FitCoeff;

	double *WorkNoise=new double[((MNStruct *)context)->pulse->nobs];
	double *powercoeff=new double[MarginCoeff];

	double tdet=0;
	double timelike=0;

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

			WorkNoise[o]=pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o]],2) + EQUAD;
			
			tdet=tdet+log(WorkNoise[o]);
			WorkNoise[o]=1.0/WorkNoise[o];
			timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	}

	double *NFd = new double[MarginCoeff];
	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[MarginCoeff];
	}

	double **NF=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		NF[i]=new double[MarginCoeff];
	}

	double **FNF=new double*[MarginCoeff];
	for(int i=0;i<MarginCoeff;i++){
		FNF[i]=new double[MarginCoeff];
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

    double *freqs = new double[MarginCoeff];

    double *DMVec=new double[((MNStruct *)context)->pulse->nobs];
    double DMKappa = 2.410*pow(10.0,-16);
    int startpos=0;
    double freqdet=0;
    double numfreqlike=0;

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

                 redamp=pow(10.0, redamp);


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
        else if(((MNStruct *)context)->incRED==4){


	        double **TempFMatrix=new double *[((MNStruct *)context)->pulse->nobs];
        	for(int i =0; i<((MNStruct *)context)->pulse->nobs; i++){
                	TempFMatrix[i]=new double[FitCoeff];
        	}

        	for(int i=0;i<FitCoeff/2;i++){
                	for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
							double freq=double(i+1)/maxtspan;
                        	TempFMatrix[k][i]=cos(2*M_PI*freq*((double)((MNStruct *)context)->pulse->obsn[k].bat));
                	}
        	}

       		for(int i=0;i<FitCoeff/2;i++){
                	for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
							double freq=double(i+1)/maxtspan;
                        	TempFMatrix[k][i+FitCoeff/2]=sin(2*M_PI*freq*((double)((MNStruct *)context)->pulse->obsn[k].bat));
                	}
        	}

			double *redcoeff=new double[FitCoeff];
			for(int i=0;i < FitCoeff; i++){
				redcoeff[i]=Cube[pcount];
				pcount++;
			}
			double *redsignal=new double[((MNStruct *)context)->pulse->nobs];
			dgemv(TempFMatrix,redcoeff,redsignal,((MNStruct *)context)->pulse->nobs,FitCoeff,'N');
			for(int i =0; i<((MNStruct *)context)->pulse->nobs; i++){
				Resvec[i] = Resvec[i]-redsignal[i];
			//	printf("signal: %i %g \n", i, redsignal[i]);
			}



           
			for (int i=0; i<FitCoeff/2; i++){

				double freq=double(i+1)/maxtspan;

				double Temppowercoeff=pow(10.0, Cube[pcount])/(maxtspan*24*60*60);
				freqdet=freqdet+2*log(Temppowercoeff);
				numfreqlike+=redcoeff[i]*redcoeff[i]/Temppowercoeff + redcoeff[i+FitCoeff/2]*redcoeff[i+FitCoeff/2]/Temppowercoeff;
				pcount++;
			}
			delete[] redcoeff;
			delete[] redsignal;
		    for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
			        delete[]TempFMatrix[j];
	  		}
			delete[]TempFMatrix;


        }
        
        else if(((MNStruct *)context)->incRED==5){
        
	        double **TempFMatrix=new double *[((MNStruct *)context)->pulse->nobs];
        	for(int i =0; i<((MNStruct *)context)->pulse->nobs; i++){
                	TempFMatrix[i]=new double[FitCoeff];
        	}

        	for(int i=0;i<FitCoeff/2;i++){
                	for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
							double freq=double(i+1)/maxtspan;
                        	TempFMatrix[k][i]=cos(2*M_PI*freq*((double)((MNStruct *)context)->pulse->obsn[k].bat));
                	}
        	}

       		for(int i=0;i<FitCoeff/2;i++){
                	for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
							double freq=double(i+1)/maxtspan;
                        	TempFMatrix[k][i+FitCoeff/2]=sin(2*M_PI*freq*((double)((MNStruct *)context)->pulse->obsn[k].bat));
                	}
        	}

			double *redcoeff=new double[FitCoeff];
			for(int i=0;i < FitCoeff; i++){
				redcoeff[i]=Cube[pcount];
				pcount++;
			}
			double *redsignal=new double[((MNStruct *)context)->pulse->nobs];
			dgemv(TempFMatrix,redcoeff,redsignal,((MNStruct *)context)->pulse->nobs,FitCoeff,'N');
			for(int i =0; i<((MNStruct *)context)->pulse->nobs; i++){
				Resvec[i] = Resvec[i]-redsignal[i];
			}

            double redamp=Cube[pcount];
            pcount++;
            double redindex=Cube[pcount];
            pcount++;

			redamp=pow(10.0, redamp);


            for (int i=0; i<FitCoeff/2; i++){

                double freq=double(i+1)/maxtspan;

                double Temppowercoeff=redamp*redamp*pow((freq*365.25),-1.0*redindex)/(maxtspan*24*60*60);
                freqdet=freqdet+2*log(Temppowercoeff);
				numfreqlike+=redcoeff[i]*redcoeff[i]/Temppowercoeff + redcoeff[i+FitCoeff/2]*redcoeff[i+FitCoeff/2]/Temppowercoeff;

             }
			delete[] redcoeff;
			delete[] redsignal;
	        for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
        	        delete[]TempFMatrix[j];
      		}
        	delete[]TempFMatrix;

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
		for(int j=0;j<MarginCoeff;j++){
// 			printf("%i %i %g %g \n",i,j,WorkNoise[i],FMatrix[i][j]);
			NF[i][j]=WorkNoise[i]*FMatrix[i][j];
		}
	}
	dgemv(NF,Resvec,NFd,((MNStruct *)context)->pulse->nobs,MarginCoeff,'T');
	dgemm(FMatrix, NF , FNF, ((MNStruct *)context)->pulse->nobs, MarginCoeff, ((MNStruct *)context)->pulse->nobs, MarginCoeff, 'T', 'N');


	double **PPFM=new double*[MarginCoeff];
	for(int i=0;i<MarginCoeff;i++){
		PPFM[i]=new double[MarginCoeff];
		for(int j=0;j<MarginCoeff;j++){
			PPFM[i][j]=0;
		}
	}


	for(int c1=0; c1<MarginCoeff; c1++){
		PPFM[c1][c1]=1.0/powercoeff[c1];
	}



	for(int j=0;j<MarginCoeff;j++){
		for(int k=0;k<MarginCoeff;k++){
			PPFM[j][k]=PPFM[j][k]+FNF[j][k];
		}
	}

       double jointdet=0;
       double Marginfreqlike=0;
       double *WorkCoeff = new double[MarginCoeff];
       for(int o1=0;o1<MarginCoeff; o1++){
                WorkCoeff[o1]=NFd[o1];
        }


	    int info=0;
        dpotrfInfo(PPFM, MarginCoeff, jointdet,info);
        dpotrs(PPFM, WorkCoeff, MarginCoeff);
        for(int j=0;j<MarginCoeff;j++){
                Marginfreqlike += NFd[j]*WorkCoeff[j];
        }
	
	lnew=-0.5*(jointdet+tdet+freqdet+timelike-Marginfreqlike+numfreqlike);

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


