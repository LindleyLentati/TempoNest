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


#include "tempo2.h"
#include "TempoNest.h"
#include "dgemm.h"
#include "dgesvd.h"
#include "dpotrf.h"
#include "dpotri.h"
#include "dgemv.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_sf_gamma.h>


double iter_factorial(unsigned int n)
{
    double ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

void readtxtoutput(std::string longname, int ndim, double **paramarray){

        int number_of_lines = 0;
	double weightsum=0;

        std::ifstream checkfile;
        std::string checkname = longname+".txt";
        checkfile.open(checkname.c_str());
        std::string line;
        while (getline(checkfile, line))
                ++number_of_lines;

//        printf("number of lines %i \n",number_of_lines);
        checkfile.close();

	std::ifstream summaryfile;
	std::string fname = longname+".txt";
	summaryfile.open(fname.c_str());

	for(int i = 0; i < ndim; i++){
		paramarray[i][0] = 0;
		paramarray[i][1] = 0;
	}

	printf("Processing Posteriors \n");
//	printf("Getting Means \n");
	double maxlike = -1.0*pow(10.0,10);
	double MAP = 0;
	for(int i=0;i<number_of_lines;i++){

		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> paramlist(begin,eof);

		weightsum += paramlist[0];
		double like = paramlist[1];

		if(like > maxlike){
			maxlike = like;
			 for(int i = 0; i < ndim; i++){
                        	 paramarray[i][2] = paramlist[i+2];
                	 }
		}

		
                if(paramlist[0] > MAP){
                         MAP = paramlist[0];
                         for(int i = 0; i < ndim; i++){
                                 paramarray[i][3] = paramlist[i+2];
                         }
                }


		for(int i = 0; i < ndim; i++){
	                paramarray[i][0] += paramlist[i+2]*paramlist[0];
        	}
	}

        for(int i = 0; i < ndim; i++){
		paramarray[i][0] = paramarray[i][0]/weightsum;
        }

	summaryfile.close();

	summaryfile.open(fname.c_str());

//	printf("Getting Errors \n");

        for(int i=0;i<number_of_lines;i++){

                std::string line;
                getline(summaryfile,line);
                std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> paramlist(begin,eof);

                for(int i = 0; i < ndim; i++){
                        paramarray[i][1] += paramlist[0]*(paramlist[i+2] - paramarray[i][0])*(paramlist[i+2] - paramarray[i][0]);
                }
        }


	for(int i = 0; i < ndim; i++){
		paramarray[i][1] = sqrt(paramarray[i][1]/weightsum);
	}

	summaryfile.close();

}




void readphyslive(std::string longname, int ndim, double **paramarray, int sampler){

        int number_of_lines = 0;

        std::ifstream checkfile;
        std::string checkname;
	if(sampler == 0){
                checkname = longname+"phys_live.points";
        }
        if(sampler == 1){
                checkname = longname+"_phys_live.txt";
        }
        checkfile.open(checkname.c_str());
        std::string line;
        while (getline(checkfile, line))
                ++number_of_lines;

        checkfile.close();

	std::ifstream summaryfile;
	std::string fname;
        if(sampler == 0){
                fname = longname+"phys_live.points";
        }
        if(sampler == 1){
                fname = longname+"_phys_live.txt";
        }


	summaryfile.open(fname.c_str());

	double *TempCube = new double[ndim];
	for(int i = 0; i < ndim; i++){
		paramarray[i][2] = 0;
	}


//	printf("Getting ML \n");
	double maxlike = -1.0*pow(10.0,10);
	for(int i=0;i<number_of_lines;i++){

		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> paramlist(begin,eof);

		double like = paramlist[ndim];

		if(like > maxlike){
			maxlike = like;
			 for(int i = 0; i < ndim; i++){
                        	 paramarray[i][2] = paramlist[i];
				TempCube[i] = paramlist[i];

                	 }
		}

	}
	summaryfile.close();

	void *tempcontext;
	double *DerivedParams=new double[1];
	int np = 1;
	double like = NewLRedMarginLogLike(ndim, TempCube, np, DerivedParams, tempcontext);

}

void readsummary(pulsar *psr, std::string longname, int ndim, void *context, long double *Tempo2Fit, int incRED, int ndims, int doTimeMargin, int doJumpMargin, int doLinear){

	int pcount=0;
	int fitcount=0;

			
	std::vector<double> paramlist(2*ndims);

	double **paramarray = new double*[ndims];
	for(int p =0;p < ndims; p++){
		paramarray[p]=new double[4];
	}


	readtxtoutput(longname, ndim, paramarray);
	readphyslive(longname, ndim, paramarray, ((MNStruct *)context)->sampler);

	int numlongparams=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;
	long double *LDP = new long double[numlongparams];
	pcount=0;
	fitcount=0;
	for(int j=0;j<((MNStruct *)context)->numFitTiming;j++){
		if(((MNStruct *)context)->Dpriors[pcount][0] != ((MNStruct *)context)->Dpriors[pcount][1]){
			double val = 0;
			if(((MNStruct *)context)->LDpriors[pcount][3] == 0){
				val = paramarray[fitcount][2];
			}
			if(((MNStruct *)context)->LDpriors[pcount][3] == 1){
				val = pow(10.0,paramarray[fitcount][2]);
			}

			LDP[j]=val*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
			fitcount++;
		}
		else if(((MNStruct *)context)->Dpriors[pcount][0] == ((MNStruct *)context)->Dpriors[pcount][1]){  
			LDP[j]=((MNStruct *)context)->Dpriors[pcount][0]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
		}
	//	printf("LD: %.25Lg %g %g \n",LDP[j], ((MNStruct *)context)->Dpriors[pcount][0],((MNStruct *)context)->Dpriors[pcount][1]); 
		pcount++;
	}

	for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
		if(((MNStruct *)context)->Dpriors[pcount][0] != ((MNStruct *)context)->Dpriors[pcount][1]){ 
			LDP[pcount]=paramarray[fitcount][2]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
			fitcount++;
		}
		else if(((MNStruct *)context)->Dpriors[pcount][0] == ((MNStruct *)context)->Dpriors[pcount][1]){
			LDP[pcount]=((MNStruct *)context)->Dpriors[pcount][0]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
		}
		pcount++;
	}
	

		
	pcount=0;
	int jj=0;
	pcount=1;
	fitcount=0;
	if(((MNStruct *)context)->LDpriors[0][2]==0)fitcount++;
	for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
	
		long double value;
		long double error;
		
		value=LDP[pcount];

		if(((MNStruct *)context)->LDpriors[j][2]==0){
			error=paramarray[fitcount][1]*(((MNStruct *)context)->LDpriors[pcount][1]);
			fitcount++;
		}
		else if(((MNStruct *)context)->LDpriors[j][2]==1){
			error=0;
		}

		
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].val[((MNStruct *)context)->TempoFitNums[pcount][1]] = value;
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].err[((MNStruct *)context)->TempoFitNums[pcount][1]] = error;
		pcount++;
	}
	
	for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
		
		long double value;
		long double error;
		
		value=LDP[pcount];
		

		if(((MNStruct *)context)->LDpriors[pcount][2]==0){
			error=paramarray[fitcount][1]*(((MNStruct *)context)->LDpriors[pcount][1]);
			fitcount++;
		}
		else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
			error=0;
		}
		
		((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[j]] = value;
		((MNStruct *)context)->pulse->jumpValErr[((MNStruct *)context)->TempoJumpNums[j]] = error;
		pcount++;
	}	
	
	

	if(((MNStruct *)context)->incBreakingIndex==1){
		long double F0 = ((MNStruct *)context)->pulse->param[param_f].val[0];
		long double F1 = ((MNStruct *)context)->pulse->param[param_f].val[1];
		long double BIndex = (long double)paramarray[fitcount][2];
		fitcount++;
		long double Kappa = -F1/pow(F0, BIndex);//(long double)pow(10.0, Cube[pcount]);
		fitcount++;
		long double F2 = -Kappa*BIndex*F1*pow(F0, BIndex-1);

		printf("%.15Lg %.15Lg %.15Lg \n", F0, F1, F2);
		((MNStruct *)context)->pulse->param[param_f].val[2] = F2;

	}

	
		

	formBatsAll(((MNStruct *)context)->pulse,1);           // Form Barycentric arrival times 
	//printf("formed bats \n");
	formResiduals(((MNStruct *)context)->pulse,1,1);       //Form residuals 
	//printf("done bats and stuff \n");	
	std::ofstream designfile;
	std::string dname = longname+"T2scaling.txt";
	
	designfile.open(dname.c_str());
	double pdParamDeriv[MAX_PARAMS];
	int numtofit=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;
	for(int i =1;i<((MNStruct *)context)->numFitTiming; i++){
		designfile << psr->param[((MNStruct *)context)->TempoFitNums[i][0]].label[((MNStruct *)context)->TempoFitNums[i][1]];
		designfile << " ";
		std::stringstream ss;
		ss.precision(std::numeric_limits<long double>::digits);//override the default
		
		ss << ((MNStruct *)context)->LDpriors[i][0];
		ss << " ";
		ss << ((MNStruct *)context)->LDpriors[i][1];
		designfile << ss.str();
		designfile << "\n";		
	} 
	
	designfile.close();
	
	double Evidence=0;
	TNtextOutput(((MNStruct *)context)->pulse, 1, 0, Tempo2Fit,  context,incRED,ndims,paramlist, Evidence, doTimeMargin, doJumpMargin, doLinear, longname, paramarray);

	
	
	
//	printf("finished output \n");

}


void getmaxlikeDM(pulsar *pulse,std::string longname, int ndim, void *context, double **paramsarray){

	formBatsAll(pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(pulse,((MNStruct *)context)->numberpulsars,1);       /* Form residuals */
//	printf("Jump1: %g \n",pulse[0].jumpVal[((MNStruct *)context)->TempoJumpNums[0]]);


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	
	int TimetoFit=((MNStruct *)context)->numFitTiming;
	int JumpstoFit=((MNStruct *)context)->numFitJumps;
	int TimetoMargin=TimetoFit+JumpstoFit;
	int numtofit= TimetoFit+JumpstoFit;
	
	printf("num params, %i %i \n", TimetoFit, JumpstoFit);

	double **TNDM=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}
	


	//getDMatrix(((MNStruct *)context)->pulse, TimetoFit, JumpstoFit, numtofit, ((MNStruct *)context)->TempoFitNums,((MNStruct *)context)->TempoJumpNums, ((MNStruct *)context)->Dpriors,1, 2, TNDM);
	getCustomDMatrixLike(context, TNDM);
	
	
	std::ofstream designfile;
	std::string dname = longname+"DesignMatrix.txt";

	designfile.open(dname.c_str());

	designfile << ((MNStruct *)context)->pulse->nobs;
	designfile << " ";
	designfile << numtofit;
	designfile << "\n";
	designfile << ((MNStruct *)context)->pulse->param[param_raj].val[0];
	designfile << " ";
	designfile << ((MNStruct *)context)->pulse->param[param_decj].val[0];
	designfile << "\n";		
	for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {

		for(int j=0; j<numtofit; j++) {
			std::stringstream ss;
			ss.precision(std::numeric_limits<double>::digits10);//override the default
			ss << TNDM[i][j];
			designfile << ss.str();
			designfile << " ";
		//	if(j==6)printf("%i %g \n", i, TNDM[i][j]);
		} 
		designfile << "\n";

	} 

	designfile.close();


	double* S = new double[TimetoMargin];
	double** U = new double*[((MNStruct *)context)->pulse->nobs];
	for(int k=0; k < ((MNStruct *)context)->pulse->nobs; k++){
		U[k] = new double[((MNStruct *)context)->pulse->nobs];
	}
	double** VT = new double*[TimetoMargin]; 
	for (int k=0; k<TimetoMargin; k++) VT[k] = new double[TimetoMargin];

	dgesvd(TNDM,((MNStruct *)context)->pulse->nobs, TimetoMargin, S, U, VT);

	double **V=new double*[numtofit];

	for(int i=0;i<numtofit;i++){
		V[i]=new double[numtofit];
	}


	for(int j=0;j < numtofit;j++){
		for(int k=0;k < numtofit;k++){
			V[j][k]=VT[k][j];
		}
	}
		
	

	for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
		for(int k=0;k < TimetoMargin;k++){
				TNDM[j][k]=U[j][k];
		}
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get White Noise vector//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  

	double **EFAC;
	double *EQUAD;

	int pcount=0;
	if(((MNStruct *)context)->numFitEFAC == 0){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)context)->systemcount; o++){
					EFAC[n-1][o]=1;
				}
			}
			else{
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
                   EFAC[n-1][o]=0;
                }
			}
		}
	}
	else if(((MNStruct *)context)->numFitEFAC == 1){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)context)->systemcount; o++){
					
					EFAC[n-1][o]=pow(10.0,paramsarray[pcount][2]);
				}
				pcount++;
			}
			else{
                                for(int o=0;o<((MNStruct *)context)->systemcount; o++){

                                        EFAC[n-1][o]=pow(10.0,paramsarray[pcount][2]);
                                }
                                pcount++;
                        }
		}
		
	}
	else if(((MNStruct *)context)->numFitEFAC > 1){
		EFAC=new double*[((MNStruct *)context)->EPolTerms];
		for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)context)->systemcount];
			if(n==1){
				for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
					EFAC[n-1][p]=pow(10.0,paramsarray[pcount][2]);
					pcount++;
				}
			}
			else{
                                for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
                                        EFAC[n-1][p]=pow(10.0,paramsarray[pcount][2]);
                                        pcount++;
                                }
                        }
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
			EQUAD[o]=pow(10.0,2*paramsarray[pcount][2]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
        EQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            EQUAD[o]=pow(10.0,2*paramsarray[pcount][2]);
	    pcount++;
        }
    }
    

    	double *SQUAD;
	if(((MNStruct *)context)->incShannonJitter == 0){
		SQUAD=new double[((MNStruct *)context)->systemcount];
		for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			SQUAD[o]=0;
		}
	}
	else if(((MNStruct *)context)->incShannonJitter == 1){
		SQUAD=new double[((MNStruct *)context)->systemcount];
                for(int o=0;o<((MNStruct *)context)->systemcount; o++){
			SQUAD[o]=pow(10.0,2*paramsarray[pcount][2]);
		}
		pcount++;
	}
	else if(((MNStruct *)context)->incShannonJitter > 1){
        SQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            SQUAD[o]=pow(10.0,2*paramsarray[pcount][2]);
	    pcount++;
        }
    }
    

	double *Noise;	
	Noise=new double[((MNStruct *)context)->pulse->nobs];

	if(((MNStruct *)context)->whitemodel == 0){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
			double EFACterm=0;
			double noiseval=0;
			double ShannonJitterTerm=0;
			
			
			if(((MNStruct *)context)->useOriginalErrors==0){
				noiseval=((MNStruct *)context)->pulse->obsn[o].toaErr;
			}
			else if(((MNStruct *)context)->useOriginalErrors==1){
				noiseval=((MNStruct *)context)->pulse->obsn[o].origErr;
			}


			for(int n=1; n <=((MNStruct *)context)->EPolTerms; n++){
				EFACterm=EFACterm + pow((noiseval*pow(10.0,-6))/pow(pow(10.0,-7),n-1),n)*EFAC[n-1][((MNStruct *)context)->sysFlags[o]];
			}	
			if(((MNStruct *)context)->incShannonJitter >0){		
		 	ShannonJitterTerm=SQUAD[((MNStruct *)context)->sysFlags[o]]*((MNStruct *)context)->TobsInfo[o]/1000.0;
			}
			Noise[o]= 1.0/(pow(EFACterm,2) + EQUAD[((MNStruct *)context)->sysFlags[o]]+ShannonJitterTerm);

		}
		
	}
	else if(((MNStruct *)context)->whitemodel == 1){
	
		for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

			Noise[o]=1.0/(EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]]*(pow(((((MNStruct *)context)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)context)->sysFlags[o]]));
		}
		
	}
	

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get F Matrix////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  

	int RedAmpParam=0;
	int RedSpecParam=0;
	int DMAmpParam=0;
	int DMSpecParam=0;

	double RedAmp=0;
	double RedIndex=0;

	double DMAmp=0;
	double DMIndex=0;

	if(((MNStruct *)context)->incRED==3){
		RedAmpParam=((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD + ((MNStruct *)context)->incShannonJitter;
		RedSpecParam=RedAmpParam+1;

		RedAmp=paramsarray[RedAmpParam][2];
		RedIndex=paramsarray[RedSpecParam][2];

		if(((MNStruct *)context)->incDM==3){
			DMAmpParam=RedSpecParam+1;
			DMSpecParam=DMAmpParam+1;

			DMAmp=paramsarray[DMAmpParam][2];
			DMIndex=paramsarray[DMSpecParam][2];
		}
	}
	else{
		if(((MNStruct *)context)->incDM==3){
			DMAmpParam=((MNStruct *)context)->numFitEFAC + ((MNStruct *)context)->numFitEQUAD + ((MNStruct *)context)->incShannonJitter;
			DMSpecParam=DMAmpParam+1;

			DMAmp=paramsarray[DMAmpParam][2];
			DMIndex=paramsarray[DMSpecParam][2];
		}
	}





	printf("params %i %i %i %i\n", RedAmpParam, RedSpecParam, DMAmpParam, DMSpecParam);
	printf("params %g %g %g %g\n", RedAmp, RedIndex, DMAmp, DMIndex);

	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);

   	int totCoeff=0;
    	if(((MNStruct *)context)->incRED != 0)totCoeff+=FitRedCoeff;
    	if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;

	double **FMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		FMatrix[i]=new double[totCoeff];
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

	double *powercoeff=new double[totCoeff];
	for(int o=0;o<totCoeff; o++){
		powercoeff[o]=0;
	}

	double Tspan = maxtspan;
	double f1yr = 1.0/3.16e7;

	if(((MNStruct *)context)->incRED==3){

		RedAmp=pow(10.0, RedAmp);

		for (int i=0; i<FitRedCoeff/2; i++){

			freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
			freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
		
			double rho = (RedAmp*RedAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-RedIndex))/(maxtspan*24*60*60);
			powercoeff[i]+= rho;
			powercoeff[i+FitRedCoeff/2]+= rho;
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


	if(((MNStruct *)context)->incDM==3){
		DMAmp=pow(10.0, DMAmp);

		for (int i=0; i<FitDMCoeff/2; i++){

			freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed +i]/maxtspan;
			freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
		
			double rho = (DMAmp*DMAmp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMIndex))/(maxtspan*24*60*60);	
			powercoeff[startpos+i]+=rho;
			powercoeff[startpos+i+FitDMCoeff/2]+=rho;
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


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get T Matrix////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	int totalsize=numtofit+totCoeff;
	double **TotalMatrix=new double*[((MNStruct *)context)->pulse->nobs];
	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		TotalMatrix[i]=new double[totalsize];
		for(int j =0;j<totalsize; j++){
			TotalMatrix[i][j]=0;
		}
	}
	

	for(int i =0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j =0;j<numtofit; j++){
			TotalMatrix[i][j]=TNDM[i][j];
		}
		
		for(int j =0;j<totCoeff; j++){
			TotalMatrix[i][j+numtofit]=FMatrix[i][j];
		}
	}


// 	printf("made NMatrix\n");
	
	double **NG = new double*[((MNStruct *)context)->pulse->nobs]; for (int k=0; k<((MNStruct *)context)->pulse->nobs; k++) NG[k] = new double[totalsize];
	double **GNG = new double*[totalsize]; for (int k=0; k<totalsize; k++) GNG[k] = new double[totalsize];



	for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		for(int j=0;j<totalsize; j++){

			NG[i][j]=TotalMatrix[i][j]*Noise[i];

		}
	}


	dgemm(TotalMatrix, NG,GNG,((MNStruct *)context)->pulse->nobs, totalsize,((MNStruct *)context)->pulse->nobs, totalsize, 'T','N');


	for(int j =0;j<totCoeff; j++){
		GNG[numtofit+j][numtofit+j]+=PPFM[j][j];
	}


	double tdet=0;
	dpotrf(GNG, totalsize, tdet);
	dpotri(GNG,totalsize);


	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)pulse->obsn[o].residual;
	}


	double *dG=new double[totalsize];
	dgemv(NG,Resvec,dG,((MNStruct *)context)->pulse->nobs,totalsize,'T');

	double *maxcoeff=new double[totalsize];
	dgemv(GNG,dG,maxcoeff,totalsize,totalsize,'N');

	double *Errorvec=new double[totalsize];

	for(int i =0; i < totalsize; i++){
		Errorvec[i]=pow(GNG[i][i], 0.5);
	}

        double *Scoeff=new double[numtofit];
        double *Serr=new double[numtofit];
        for(int i =0; i < numtofit; i++){
                Scoeff[i]=maxcoeff[i]/S[i];
                Serr[i]=Errorvec[i]/S[i];
        }

        double *TempoCoeff = new double[numtofit];
        double *TempoErr =  new double[numtofit];
        dgemv(V,Scoeff,TempoCoeff,numtofit,numtofit, 'N');
        dgemv(V,Serr,TempoErr,numtofit,numtofit, 'N');

        for(int i=0;i<numtofit; i++){
                double errsum=0;
          for(int j=0;j<numtofit; j++){
                                errsum += pow(V[i][j]*Serr[j],2);
                }
          TempoErr[i]=pow(errsum,0.5);
        }
        //updateParameters(((MNStruct *)context)->pulse,0,TempoCoeff,TempoErr);
	 for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
                        //printf("%i %.25Lg \n", p, ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]]);
                }
                for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
                        //printf("%i %g \n", p,((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[p]]);
                }




 	char name[1000];

        std::ofstream Finalfile;
        std::string Finalfilename = longname+"Final.dat";
        Finalfile.open(Finalfilename.c_str());


        for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
                double dsum=0;
                for(int j=0;j<totalsize; j++){
                        dsum=dsum+TotalMatrix[i][j]*maxcoeff[j];
                }
		Finalfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].bat << " ";
        	Finalfile << std::scientific << Resvec[i] << " " << dsum <<  " " << pow(1.0/Noise[i],0.5) << " ";
		Finalfile << std::fixed  << std::setprecision(1)  <<(double)((MNStruct *)context)->pulse->obsn[i].freqSSB << " ";
		Finalfile << ((MNStruct *)context)->sysFlags[i]  << "\n";


        }
	Finalfile.close();

        std::ofstream Resfile;
        std::string Resfilename = longname+"Res.dat";
        Resfile.open(Resfilename.c_str());
        Resfile << ((MNStruct *)context)->pulse->name << "\n";
        Resfile << "obsid    Freq (MHz)   SAT    Res (s)    Error (s)\n";

        for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
                double dsum=0;
                for(int j=0;j<numtofit; j++){
                        dsum=dsum+TotalMatrix[i][j]*maxcoeff[j];
                }
                sscanf(((MNStruct *)context)->pulse->obsn[i].fname,"%s",name);
                Resfile <<  name << " ";
                Resfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " ";
                Resfile << std::scientific << Resvec[i]-dsum <<  " " << pow(1.0/Noise[i],0.5) << "\n";
        }

        Resfile.close();



	if(((MNStruct *)context)->incDM==3){
		std::ofstream DMfile;
		std::string DMfilename = longname+"DM.dat";
		DMfile.open(DMfilename.c_str());
		DMfile << ((MNStruct *)context)->pulse->name << "\n";
		DMfile << "obsid    Freq (MHz)   SAT    DMRes (s)    Error (s)  DMVal    Error\n";


		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			double dsum=0;
			double scaledsum=0;
			double errsum=0;
			double scaleerrsum=0;
			for(int j=0;j<FitDMCoeff; j++){
				dsum=dsum+TotalMatrix[i][j+numtofit+FitRedCoeff]*maxcoeff[j+numtofit+FitRedCoeff];
				scaledsum=scaledsum+TotalMatrix[i][j+numtofit+FitRedCoeff]*maxcoeff[j+numtofit+FitRedCoeff]/DMVec[i];
				errsum=errsum+pow(Errorvec[j+numtofit+FitRedCoeff]*TotalMatrix[i][j+numtofit+FitRedCoeff],2);
				scaleerrsum=scaleerrsum+pow(Errorvec[j+numtofit+FitRedCoeff]*TotalMatrix[i][j+numtofit+FitRedCoeff]/DMVec[i],2);
			}
		
			double freq=(double)((MNStruct *)context)->pulse->obsn[i].freqSSB;
			long double yrs = (((MNStruct *)context)->pulse->obsn[i].sat - ((MNStruct *)context)->pulse->param[param_dmepoch].val[0])/365.25;
	 	        long double arg = 1.0;
			double dmDot=0;
			double dmDotErr=0;
		        for (int d=1;d<9;d++){	
		  		arg *= yrs;
		  		if (((MNStruct *)context)->pulse->param[param_dm].paramSet[d]==1){
		    			dmDot+=(double)(((MNStruct *)context)->pulse->param[param_dm].val[d]*arg);
					dmDotErr+=pow((double)(((MNStruct *)context)->pulse->param[param_dm].err[d]*arg),2);
				}
			}
			dsum=dsum+dmDot*DMVec[i];
			scaledsum=scaledsum+dmDot;
			errsum=errsum+dmDotErr*pow(DMVec[i],2);
			scaleerrsum=scaleerrsum+dmDotErr;
			double restsum=0;
			for(int j=0;j<numtofit+FitRedCoeff; j++){
				restsum=restsum+TotalMatrix[i][j]*maxcoeff[j];
			}
			sscanf(((MNStruct *)context)->pulse->obsn[i].fname,"%s",name);
			DMfile <<  name << " ";
			DMfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " "; 
			DMfile << std::scientific << dsum <<  " " << pow(errsum,0.5) << " " << scaledsum << " " << pow(scaleerrsum,0.5) << "\n";
		}

		DMfile.close();

	}

	if(((MNStruct *)context)->incRED==3){

		std::ofstream Redfile;
		std::string Redfilename = longname+"Red.dat";
		Redfile.open(Redfilename.c_str());
		Redfile << ((MNStruct *)context)->pulse->name << "\n";
		Redfile << "obsid    Freq (MHz)   SAT    RedRes (s)    Error (s) \n";

		for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
		        double dsum=0;
		        double errsum=0;
		        for(int j=0;j<FitRedCoeff; j++){
		                dsum=dsum+TotalMatrix[i][j+numtofit]*maxcoeff[j+numtofit];
		                errsum=errsum+pow(Errorvec[j+numtofit]*TotalMatrix[i][j+numtofit],2);
		        }
		        sscanf(((MNStruct *)context)->pulse->obsn[i].fname,"%s",name);
		        Redfile <<  name << " ";
		        Redfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " ";
		        Redfile << std::scientific << dsum <<  " " << pow(errsum,0.5) << "\n";
		}

		Redfile.close();
	}


	delete[]S;	

	for (int j = 0; j < TimetoMargin; j++){
		delete[]VT[j];
	}
	
	delete[]VT;

	for (int j = 0; j < ((MNStruct *)context)->pulse->nobs; j++){
		delete[]U[j];
	}
	delete[]U;

}



void getNGJitterMatrixEpochs(pulsar *pulse, int &NumEpochs){


	// count ECORR values
	if(pulse->nTNECORR > 0){
		for(int i=0; i<pulse->nTNECORR; i++){
			printf("\nIncluding ECORR value for backend %s: %g mus", \
					pulse->TNECORRFlagVal[i], pulse->TNECORRVal[i]);

		}
	}


	// find number of epochs (default dt= 10 s)
	int *Processed = new int[pulse->nobs];

	// initialize processed flags
	for (int i=0;i<pulse->nobs;i++){
		Processed[i] = 1;
	}

	// make sure we only process the epochs with the chosen flags
	for (int i=0;i<pulse->nobs;i++){
		for (int j=0;j<pulse->obsn[i].nFlags;j++){
			for (int k=0;k<pulse->nTNECORR;k++){
				if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
					if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
						Processed[i] = 0;
					}
				}
			}
		}
	}


	double dt = 10.0 / SECDAY;
	double satmin;
	double satmax;
	int nepoch = 0;
	int in = 0;
	int allProcessed = 0;
	while (!allProcessed){
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				satmin = (double)pulse->obsn[i].bat - dt;
				satmax = (double)pulse->obsn[i].bat + dt;
				break;
			}
		}
		for (int i=0;i<pulse->nobs;i++){
			for (int j=0;j<pulse->obsn[i].nFlags;j++){
				for (int k=0;k<pulse->nTNECORR;k++){
					if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
						if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
							if ((double)pulse->obsn[i].bat > satmin && \
									(double)pulse->obsn[i].bat < satmax){
								Processed[i] = 1;
								in++;
							}
						}
					}
				}
			}
		}
		if (in != 0){
			nepoch++;
			in = 0;
		}
		allProcessed = 1;
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				allProcessed = 0;
				break;
			}
		}
	}


	if (nepoch > 0){
		printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
	}

	NumEpochs=nepoch;


}

void getNGJitterMatrix(pulsar *pulse, double **JitterMatrix, int &NumEpochs){


	// count ECORR values
	if(pulse->nTNECORR > 0){
		for(int i=0; i<pulse->nTNECORR; i++){
			printf("\nIncluding ECORR value for backend %s: %g mus", \
					pulse->TNECORRFlagVal[i], pulse->TNECORRVal[i]);

		}
	}


	// find number of epochs (default dt= 10 s)
	int *Processed = new int[pulse->nobs];

	// initialize processed flags
	for (int i=0;i<pulse->nobs;i++){
		Processed[i] = 1;
	}

	// make sure we only process the epochs with the chosen flags
	for (int i=0;i<pulse->nobs;i++){
		for (int j=0;j<pulse->obsn[i].nFlags;j++){
			for (int k=0;k<pulse->nTNECORR;k++){
				if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
					if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
						Processed[i] = 0;
					}
				}
			}
		}
	}


	double dt = 10.0 / SECDAY;
	double satmin;
	double satmax;
	int nepoch = 0;
	int in = 0;
	int allProcessed = 0;
	while (!allProcessed){
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				satmin = (double)pulse->obsn[i].bat - dt;
				satmax = (double)pulse->obsn[i].bat + dt;
				break;
			}
		}
		for (int i=0;i<pulse->nobs;i++){
			for (int j=0;j<pulse->obsn[i].nFlags;j++){
				for (int k=0;k<pulse->nTNECORR;k++){
					if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
						if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
							if ((double)pulse->obsn[i].bat > satmin && \
									(double)pulse->obsn[i].bat < satmax){
								Processed[i] = 1;
								in++;
							}
						}
					}
				}
			}
		}
		if (in != 0){
			nepoch++;
			in = 0;
		}
		allProcessed = 1;
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				allProcessed = 0;
				break;
			}
		}
	}


	if (nepoch > 0){
		printf("\n\nUsing %d epochs for PSR %s\n\n", nepoch, pulse->name);
	}

	NumEpochs=nepoch;



	// initialize processed flags
	for (int i=0;i<pulse->nobs;i++){
		Processed[i] = 1;
	}

	// make sure we only process the epochs with the chosen flags
	for (int i=0;i<pulse->nobs;i++){
		for (int j=0;j<pulse->obsn[i].nFlags;j++){
			for (int k=0;k<pulse->nTNECORR;k++){
				if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
					if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
						Processed[i] = 0;
					}
				}
			}
		}
	}


	nepoch = 0;
	in = 0;
	allProcessed = 0;
	while (!allProcessed){
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				satmin = (double)pulse->obsn[i].bat - dt;
				satmax = (double)pulse->obsn[i].bat + dt;
				break;
			}
		}
		for (int i=0;i<pulse->nobs;i++){
			for (int j=0;j<pulse->obsn[i].nFlags;j++){
				for (int k=0;k<pulse->nTNECORR;k++){
					if (strcmp(pulse->obsn[i].flagID[j], pulse->TNECORRFlagID[k])==0){
						if (strcmp(pulse->obsn[i].flagVal[j],pulse->TNECORRFlagVal[k])==0){
							if ((double)pulse->obsn[i].bat > satmin && \
									(double)pulse->obsn[i].bat < satmax){
								Processed[i] = 1;
								JitterMatrix[i][nepoch] = 1;
								in++;
							}
							else{
								JitterMatrix[i][nepoch] = 0;
							}
						}
					}
				}
			}
		}
		if (in != 0){
			nepoch++;
			in = 0;
		}
		allProcessed = 1;
		for (int i=0;i<pulse->nobs;i++){
			if (Processed[i]==0){
				allProcessed = 0;
				break;
			}
		}
	}




}


void getCustomDMatrix(pulsar *pulse, int *MarginList, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int incDM, int TimetoFit, int JumptoFit){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix


		int pcount=1;
		int numToMargin=1;
		
		Dpriors[0][0]=0;
		Dpriors[0][1]=0;
		
		for (int p=1;p<TimetoFit;p++) {
			if(MarginList[pcount]!=1){
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
			}
			else if(MarginList[pcount]==1){
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
				Dpriors[pcount][0]=0;
				Dpriors[pcount][1]=0;
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
				Dpriors[pcount][0]=0;
				Dpriors[pcount][1]=0;
				numToMargin++;
			}
			pcount++;
		}
	
	
//		for(int i=0; i < pulse->nobs; i++) {
//			FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin, pulse, i,0);
//			for(int j=0; j<numToMargin; j++) {
//				TNDM[i][j]=pdParamDeriv[j];
//					//printf("Dmatrix: %i %i %22.20e \n", i,j,pdParamDeriv[j]);
//			} 
//		} 

		//Now set fit flags back to how they were
	
		for (int p=1;p<TimetoFit;p++) {
			pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 1;
		}
	
		for(int i=0; i < JumptoFit; i++){
			pulse[0].fitJump[TempoJumpNums[i]]=1;
		}
}


void getCustomDVectorLike(void *context, double *TNDM){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix

		int pcount=1;
		int numToMargin=1;
		
		for (int p=1;p<((MNStruct *)context)->numFitTiming;p++) {
			if(((MNStruct *)context)->LDpriors[pcount][2]!=1){
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 0;
			}
			else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 1;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < ((MNStruct *)context)->numFitJumps; i++){
			if(((MNStruct *)context)->LDpriors[pcount][2]!=1){
				((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=0;
			}
			else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
				((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=1;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
			FITfuncs(((MNStruct *)context)->pulse->obsn[i].bat - ((MNStruct *)context)->pulse->param[param_pepoch].val[0], pdParamDeriv, numToMargin, ((MNStruct *)context)->pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
				//printf("CDML: %i %i %i %g\n", i,j,numToMargin,pdParamDeriv[j]);
				TNDM[i + j*((MNStruct *)context)->pulse->nobs]=pdParamDeriv[j];
			} 
		} 
		//Now set fit flags back to how they were
	

		for (int p=1;p<((MNStruct *)context)->numFitTiming;p++) {
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 1;

		}
	
		for(int i=0; i < ((MNStruct *)context)->numFitJumps; i++){
			((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=1;
		}
	
}
void getCustomDMatrixLike(void *context, double **TNDM){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix

		int pcount=1;
		int numToMargin=1;
		
		for (int p=1;p<((MNStruct *)context)->numFitTiming;p++) {
			if(((MNStruct *)context)->LDpriors[pcount][2]!=1){
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 0;
			}
			else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 1;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < ((MNStruct *)context)->numFitJumps; i++){
			if(((MNStruct *)context)->LDpriors[pcount][2]!=1){
				((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=0;
			}
			else if(((MNStruct *)context)->LDpriors[pcount][2]==1){
				((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=1;
				numToMargin++;
			}
			pcount++;
		}
	
		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
			FITfuncs(((MNStruct *)context)->pulse->obsn[i].bat - ((MNStruct *)context)->pulse->param[param_pepoch].val[0], pdParamDeriv, numToMargin, ((MNStruct *)context)->pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
//				printf("CDML: %i %i %i %g\n", i,j,numToMargin,pdParamDeriv[j]);
				TNDM[i][j]=pdParamDeriv[j];
			} 
		} 
		//Now set fit flags back to how they were
	

		for (int p=1;p<((MNStruct *)context)->numFitTiming;p++) {
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].fitFlag[((MNStruct *)context)->TempoFitNums[p][1]] = 1;

		}
	
		for(int i=0; i < ((MNStruct *)context)->numFitJumps; i++){
			((MNStruct *)context)->pulse->fitJump[((MNStruct *)context)->TempoJumpNums[i]]=1;
		}
	
}




void getDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix
	//If we are marginalising, change the prior so you dont sample it
	if(JumptoFit>0){
		if(doJumpMargin==0){
			for(int i=0; i < JumptoFit; i++){
				pulse[0].fitJump[TempoJumpNums[i]]=0;
			}
		} 
		else if(doJumpMargin==1){
			for(int i=0; i < JumptoFit; i++){
// 				Dpriors[TimetoFit+i][0]=0;
// 				Dpriors[TimetoFit+i][1]=0;
				
			}
		}
	}


	if(doTimeMargin==0){
		for (int p=1;p<TimetoFit;p++) {
			pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
		}
	}
	//If only marginalising over QSD (i.e. F0 and F1) then set all other flags to 0
	else if(doTimeMargin==1){	
		int pcount=0;
		for (int p=0;p<MAX_PARAMS;p++) {
			if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")!=0){
				for (int k=0;k<pulse[0].param[p].aSize;k++){
					if(pulse[0].param[p].fitFlag[k] == 1){
						//printf("in getD %s \n",pulse[0].param[p].shortlabel[k]);
						pulse[0].param[p].fitFlag[k]=0;
						pcount++;
					}
				}
			}
			else if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")==0){
				for (int k=0;k<pulse[0].param[p].aSize;k++){
					if(pulse[0].param[p].fitFlag[k] == 1){
						//printf("in getD %s setting prior %i to zero\n",pulse[0].param[p].shortlabel[k], pcount);
// 						Dpriors[pcount][0]=0;
// 						Dpriors[pcount][1]=0;
						pcount++;
					}
				}
			}
					
		}
	}
	else if(doTimeMargin==2){	
		for(int i=0; i < TimetoFit; i++){
// 			Dpriors[i][0]=0;
// 			Dpriors[i][1]=0;
		}
	}


	for(int i=0; i < pulse->nobs; i++) {
		FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin-1, pulse, i,0);
		for(int j=0; j<numToMargin; j++) {
			TNDM[i][j]=pdParamDeriv[j];
 			printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);

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


void getMarginDMatrix(pulsar *pulse, int TimetoFit, int JumptoFit, int numToMargin, int **TempoFitNums, int *TempoJumpNums, double **Dpriors, int doJumpMargin, int doTimeMargin, double **TNDM, int linearFit){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//Unset all fit flags for parameters we arn't marginalising over so they arn't in the design Matrix
	//If we are marginalising, change the prior so you dont sample it


		if(JumptoFit>0){
		
				Dpriors[0][0]=0;  //Prior on Phase is zero
				Dpriors[0][1]=0;
				
				
			if(doJumpMargin==0){
				for(int i=0; i < JumptoFit; i++){
					pulse[0].fitJump[TempoJumpNums[i]]=0;
				}
			} 
			else if(doJumpMargin==1){
				for(int i=0; i < JumptoFit; i++){
					Dpriors[TimetoFit+i][0]=0;
					Dpriors[TimetoFit+i][1]=0;
					
				}
			}
		}
	
	
		if(doTimeMargin==0){
			for (int p=1;p<TimetoFit;p++) {
				pulse[0].param[TempoFitNums[p][0]].fitFlag[TempoFitNums[p][1]] = 0;
			}
		}
		//If only marginalising over QSD (i.e. F0 and F1) then set all other flags to 0
		else if(doTimeMargin==1){	
			int pcount=0;
			Dpriors[pcount][0]=0;
			Dpriors[pcount][1]=0;
			pcount++;
			for (int p=0;p<MAX_PARAMS;p++) {
				if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")!=0){
					for (int k=0;k<pulse[0].param[p].aSize;k++){
						if(pulse[0].param[p].fitFlag[k] == 1){
							//printf("in getDNM %s \n",pulse[0].param[p].shortlabel[k]);
							pulse[0].param[p].fitFlag[k]=0;
							pcount++;
						}
					}
				}
				else if(strcasecmp(pulse[0].param[p].shortlabel[0],"F0")==0){
					for (int k=0;k<pulse[0].param[p].aSize;k++){
						if(pulse[0].param[p].fitFlag[k] == 1){
							//printf("in getDNM %s setting prior %i to zero\n",pulse[0].param[p].shortlabel[k], pcount);
							Dpriors[pcount][0]=0;
							Dpriors[pcount][1]=0;
							pcount++;
						}
					}
				}
						
			}
		}
		else if(doTimeMargin==2){	
			for(int i=0; i < TimetoFit; i++){
				Dpriors[i][0]=0;
				Dpriors[i][1]=0;
			}
		}
	
	
		for(int i=0; i < pulse->nobs; i++) {
			FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin, pulse, i,0);
			for(int j=0; j<numToMargin; j++) {
				TNDM[i][j]=pdParamDeriv[j];
 			//	printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);
	
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



void makeGDesign(pulsar *pulse, int &Gsize, int numtofit, double** staticGMatrix, double **oneDesign){

	
	int numpulsars=1;
	double*** Umatrices=new double**[numpulsars];
	for(int i=0;i<numpulsars;i++){
		Umatrices[i]=new double*[pulse->nobs];
		for(int j=0;j<pulse->nobs;j++){
			Umatrices[i][j]=new double[pulse->nobs];
		}
	}

	for(int i=0;i<numpulsars;i++){

		double* S = new double[numtofit];
		double** U = new double*[pulse->nobs]; for (int k=0; k<pulse->nobs; k++) U[k] = new double[pulse->nobs];
		double** VT = new double*[numtofit]; for (int k=0; k<numtofit; k++) VT[k] = new double[numtofit];

		dgesvd(oneDesign,pulse->nobs, numtofit, S, U, VT);

		for(int j=0;j<pulse->nobs;j++){
// 			printf("h %i %i\n",i,j);
			for(int k=0;k<pulse->nobs;k++){
 				//printf("U %i %i %i %i %g \n",pulse->nobs, numtofit, j,k,U[j][k]);
				Umatrices[i][j][k]=U[j][k];
				
			}
		}

		delete[]S;	
		for (int j = 0; j < pulse->nobs; j++){
			delete[]U[j];
		}
		for (int j = 0; j < numtofit; j++){
			delete[]VT[j];
		}
		delete[]U;
		delete[]VT;
	
		
	}

// 	printf("Done SVD ALL");
	int Dsize=0;
	for(int i=0;i<numpulsars;i++){
		Dsize += numtofit;
	}




	int Osum=0;
	int Gsum=0;
	for(int i=0;i<numpulsars;i++){
		for(int j=0;j<pulse->nobs;j++){
			int nfsum=0;
			for(int k=0;k<pulse->nobs;k++){
				if(k>=numtofit){
					staticGMatrix[Osum+j][Osum-Gsum+nfsum]=Umatrices[i][j][k];
 					//printf("GBack %i %i %g \n",Osum+j,Osum-Gsum+nfsum,staticGMatrix[Osum+j][Osum-Gsum+nfsum]);
					nfsum++;

				}
			}
		}
		Gsum+=numtofit;
		Osum+=pulse->nobs;
	}


	Gsize=Osum-Gsum;
/*	for(int i=0;i<pulse->nobs;i++){
		for(int j=0;j<Gsize;j++){
			printf("%i %i %g \n ",i,j,staticGMatrix[i][j]);
		}
	}
	printf("%i %i \n",pulse->nobs,Gsize)*/;

// 	printf("Done SVD ALL");
}

void makeStaticGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** staticGMatrix, double &tdet){


	double *Noise=new double[pulse->nobs];
	for(int o=0;o < pulse->nobs; o++){
		Noise[o]=pow(pulse->obsn[o].toaErr*pow(10.0,-6),2);
	}

	double **NG = new double*[pulse->nobs]; for (int k=0; k< pulse->nobs; k++) NG[k] = new double[Gsize];
	for(int i=0;i< pulse->nobs;i++){
		for(int j=0;j< Gsize; j++){
			NG[i][j] = GMatrix[i][j]*Noise[i];
		}
	}

	double **GG = new double*[Gsize]; for (int k=0; k< Gsize; k++) GG[k] = new double[Gsize];

	dgemm(GMatrix, NG,GG,pulse->nobs, Gsize,pulse->nobs, Gsize, 'T','N');
	
	tdet=0;
	dpotrf(GG, Gsize, tdet);
	dpotri(GG,Gsize);
	


	dgemm(GMatrix, GG,NG, pulse->nobs, Gsize, Gsize, Gsize, 'N','N');

	dgemm(NG, GMatrix, staticGMatrix, pulse->nobs, Gsize, pulse->nobs, Gsize, 'N','T');
	
	delete[] Noise;
	
	for (int j = 0; j < pulse->nobs; j++){
		delete[] NG[j];
	}
	delete[] NG;

	for (int j = 0; j < Gsize; j++){
		delete[]GG[j];
	}
	delete[] GG;	
	
	
}

void makeStaticDiagGMatrix(pulsar *pulse, int Gsize, double **GMatrix, double** GNMatrix, double *SVec){


		




	double *Noise=new double[pulse->nobs];
	for(int o=0;o < pulse->nobs; o++){
		Noise[o]=pow(pulse->obsn[o].toaErr*pow(10.0,-6),2);
	}

	double **NG = new double*[pulse->nobs]; for (int k=0; k< pulse->nobs; k++) NG[k] = new double[Gsize];
	for(int i=0;i< pulse->nobs;i++){
		for(int j=0;j< Gsize; j++){
			NG[i][j] = GMatrix[i][j]*Noise[i];

		}
	}

	double **GG = new double*[Gsize]; for (int k=0; k< Gsize; k++) GG[k] = new double[Gsize];

	dgemm(GMatrix, NG,GG,pulse->nobs, Gsize,pulse->nobs, Gsize, 'T','N');
	

	double** U = new double*[Gsize]; for (int k=0; k<Gsize; k++) U[k] = new double[Gsize];
	double** VT = new double*[Gsize]; for (int k=0; k<Gsize; k++) VT[k] = new double[Gsize];

	dgesvd(GG,Gsize, Gsize, SVec, U, VT);
	
	double **GT = new double*[Gsize]; for (int k=0; k< Gsize; k++) GT[k] = new double[pulse->nobs];
	
	for(int i=0;i< Gsize;i++){
			for(int j=0;j< pulse->nobs; j++){
				GT[i][j]=GMatrix[j][i];
				
			}
		}
		
		
		
	dgemm(VT, GT,GNMatrix,Gsize,Gsize, Gsize, pulse->nobs, 'N','N');
	

	delete[] Noise;
	
	for (int j = 0; j < pulse->nobs; j++){
		delete[] NG[j];
		
	}
	delete[] NG;

	for (int j = 0; j < Gsize; j++){
		delete[]GG[j];
		delete[]VT[j];
		delete[] U[j];
		delete[] GT[j];

	}
	delete[] GG;	
	delete[]VT;
	delete[] U;
	delete[] GT;
	
}




