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
#include "/usr/include/gsl/gsl_sf_gamma.h"

void writeCov(std::vector<double> Cube, int &ndim,void *context, std::string longname, int outputoption);

double iter_factorial(unsigned int n)
{
    double ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}


	

void readsummary(pulsar *psr, std::string longname, int ndim, void *context, long double *Tempo2Fit, int incRED, int ndims, int doTimeMargin, int doJumpMargin, int doLinear){

	int number_of_lines = 0;
	int outputoption=0;
	char *outname;
	int pcount=0;

	std::ifstream checkfile;
	std::string checkname = longname+"summary.txt";
	checkfile.open(checkname.c_str());
	std::string line;
	while (getline(checkfile, line))
		++number_of_lines;

	printf("number of lines %i \n",number_of_lines);
	checkfile.close();
	
	if(number_of_lines == 1){
		std::ifstream summaryfile;
		std::string fname = longname+"summary.txt";
		summaryfile.open(fname.c_str());
	
		for(int i=0;i<number_of_lines;i++){
			
			std::string line;
			getline(summaryfile,line);
			//printf("line %s \n", line);
			std::istringstream myStream( line );
			std::istream_iterator< double > begin(myStream),eof;
			std::vector<double> paramlist(begin,eof);

	
			long double *LDP = new long double[ndim];

			pcount=0;
			for(int j=0;j<((MNStruct *)context)->numFitTiming;j++){
						
					LDP[j]=paramlist[j+outputoption*ndim]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
					pcount++;
			}
			for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
					LDP[pcount]=paramlist[pcount+outputoption*ndim]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
					pcount++;
			}
			
	
				
			pcount=0;
			int jj=0;
			if(doLinear==1){
			
				double *linearParams=new double[ndim];
				double *nonLinearParams=new double[ndim];
				double *errorvec=new double[ndim];
				double *nonLinearerror=new double[ndim];
				for(int j=0;j<ndim;j++){
					linearParams[j]=paramlist[j+outputoption*ndim];
					errorvec[j]=paramlist[j+ndim];
					printf("%i %g %g \n",j,linearParams[j],errorvec[j]);
					nonLinearerror[j]=0;
					nonLinearParams[j]=0;
				}
				TNupdateParameters(psr,0,linearParams,errorvec, nonLinearParams);
				TNupdateParameters(psr,0,errorvec,errorvec, nonLinearerror);
                for(int j=0;j<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j++){
                        if(((MNStruct *)context)->LDpriors[j][1] != 0){
						 printf("before %i %g %g %Lg\n",j,nonLinearParams[j],nonLinearerror[j],((MNStruct *)context)->LDpriors[j][1]);
						nonLinearParams[j]=nonLinearParams[j]/((MNStruct *)context)->LDpriors[j][1];
						nonLinearerror[j]=nonLinearerror[j]/((MNStruct *)context)->LDpriors[j][1];
						printf("after %i %g %g %Lg\n",j,nonLinearParams[j],nonLinearerror[j],((MNStruct *)context)->LDpriors[j][1]);
					}
					else{
						nonLinearParams[j]=0;
						nonLinearerror[j]=0;
					}
                }

				pcount=1;
                for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){

                        long double value=nonLinearParams[j]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
                        long double error=nonLinearerror[j]*(((MNStruct *)context)->LDpriors[pcount][1]);
                        ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].val[((MNStruct *)context)->TempoFitNums[pcount][1]] = value;
                        ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].err[((MNStruct *)context)->TempoFitNums[pcount][1]] = error;
                        printf("%i %Lg %Lg %Lg %Lg \n",j,(((MNStruct *)context)->LDpriors[pcount][0]),(((MNStruct *)context)->LDpriors[pcount][1]),value,error);
                        printf("%i %g \n", j, paramlist[j]);
                        pcount++;
                }

                       for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){

                                long double value=nonLinearParams[pcount]*(((MNStruct *)context)->LDpriors[pcount][1])+(((MNStruct *)context)->LDpriors[pcount][0]);
                                long double error=nonLinearerror[pcount]*(((MNStruct *)context)->LDpriors[pcount][1]);
                                ((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[j]] = value;
                                ((MNStruct *)context)->pulse->jumpValErr[((MNStruct *)context)->TempoJumpNums[j]] = error;
                                pcount++;
                        }
			}
			else{
			pcount=1;
			for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			
				long double value;
				long double error;
				
				value=LDP[pcount];
	
				if(doTimeMargin==0){

				error=paramlist[j+ndim]*(((MNStruct *)context)->LDpriors[pcount][1]);
				}
				else if(doTimeMargin==1){
					if(strcasecmp(((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]],"F0")==0 || strcasecmp(((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].shortlabel[((MNStruct *)context)->TempoFitNums[pcount][1]],"F1")==0){
						error=0;
					}
					else{


						error=paramlist[j+ndim]*(((MNStruct *)context)->LDpriors[pcount][1]);
					}
				}
				else if(doTimeMargin==2){
					
					error=0;
				}
				
				
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].val[((MNStruct *)context)->TempoFitNums[pcount][1]] = value;
				((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[pcount][0]].err[((MNStruct *)context)->TempoFitNums[pcount][1]] = error;
				printf("%i %Lg %Lg %Lg %Lg \n",j,(((MNStruct *)context)->LDpriors[pcount][0]),(((MNStruct *)context)->LDpriors[pcount][1]),value,error);
				printf("%i %.25Lg %.25Lg \n", j, value, error);
				pcount++;
			}
			
			for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
				
				long double value;
				long double error;
				
				value=LDP[pcount];
				
				if(doJumpMargin==0){

					error=paramlist[pcount+ndim]*(((MNStruct *)context)->LDpriors[pcount][1]);
				}
				else{
				
				error=0;
				}
				
				((MNStruct *)context)->pulse->jumpVal[((MNStruct *)context)->TempoJumpNums[j]] = value;
				((MNStruct *)context)->pulse->jumpValErr[((MNStruct *)context)->TempoJumpNums[j]] = error;
				printf("%i %Lg %Lg %Lg %Lg \n",pcount,(((MNStruct *)context)->LDpriors[pcount][0]),(((MNStruct *)context)->LDpriors[pcount][1]),value,error);
				pcount++;
			}	
			}
			
				
		 
		  	
			
	
	
			formBatsAll(((MNStruct *)context)->pulse,1);           // Form Barycentric arrival times 
			formResiduals(((MNStruct *)context)->pulse,1,1);       //Form residuals 
			
			std::ofstream designfile;
			std::string dname = longname+"designMatrix.txt";
			std::ofstream resfile;
			std::string rname = longname+"res.dat";
			printf("Writing Timing Model Design Matrix and Residuals\n");

			designfile.open(dname.c_str());
			resfile.open(rname.c_str());
			double pdParamDeriv[MAX_PARAMS];
			int numtofit=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;
			designfile << psr->nobs;
			designfile << " ";
			designfile << numtofit;
			designfile << "\n";
			designfile << psr->param[param_raj].val[0];
			designfile << " ";
			designfile << psr->param[param_decj].val[0];
			designfile << "\n";		
			for(int i=0; i < psr->nobs; i++) {

				std::stringstream rs;
				rs.precision(std::numeric_limits<double>::digits10);//override the default

				rs << psr->obsn[i].bat;
				rs << " ";
				rs << psr->obsn[i].residual;
				rs << " ";
				rs << psr->obsn[i].toaErr*pow(10.0,-6);
				resfile << rs.str();
				resfile << "\n";		

				FITfuncs(psr[0].obsn[i].bat - psr[0].param[param_pepoch].val[0], pdParamDeriv, numtofit, psr, i,0);
				for(int j=0; j<numtofit; j++) {
					std::stringstream ss;
					ss.precision(std::numeric_limits<double>::digits10);//override the default
					ss << pdParamDeriv[j];
	
					designfile << ss.str();
					designfile << " ";
	 				//printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);
	
				} 
				designfile << "\n";
	
			} 

			designfile.close();
			resfile.close();
		
			writeCov(paramlist, ndim,context, longname,outputoption);

			double Evidence=paramlist[4*ndim];
			TNtextOutput(((MNStruct *)context)->pulse, 1, 0, Tempo2Fit,  context,incRED,ndims,paramlist, Evidence, doTimeMargin, doJumpMargin, doLinear, longname);
			
			
			printf("finished output \n");
		}
		summaryfile.close();

	}
	else{
		printf("More than one mode has been detected, you'll have to handle that on your own for now.");
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
		FITfuncs(pulse[0].obsn[i].bat - pulse[0].param[param_pepoch].val[0], pdParamDeriv, numToMargin, pulse, i,0);
		for(int j=0; j<numToMargin; j++) {
			TNDM[i][j]=pdParamDeriv[j];
 			//printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);

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
// 				printf("%i %i %22.20e \n", i,j,pdParamDeriv[j]);
	
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
// 				printf("%i %i %i %g \n",i,j,k,U[j][k]);
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
// 					printf("GBack %i %i %g \n",Osum+j,Osum-Gsum+nfsum,GBack[Osum+j][Osum-Gsum+nfsum]);
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


void convertFromLinear(pulsar *psr, std::string longname, int ndim, void *context){

	int number_of_lines = 0;
	char *outname;

	std::ifstream checkfile;
	std::string checkname = longname+".txt";
	checkfile.open(checkname.c_str());
	std::string line;
	while (getline(checkfile, line))
		++number_of_lines;

	printf("CFL number of lines %i \n",number_of_lines);
	checkfile.close();
	
	std::ifstream summaryfile;
	std::string fname = longname+".txt";
	summaryfile.open(fname.c_str());

	std::ofstream newtxtfile;
	std::string newtxtname = longname+"NL.txt";
	newtxtfile.open(newtxtname.c_str());

	double *linearParams = new double[((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps+1];
	double *nonLinearParams = new double[((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps+1];
	double *errorvec = new double[((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps+1];
	for(int i=0;i<number_of_lines;i++){
		
		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
		std::istream_iterator< double > begin(myStream),eof;
		std::vector<double> paramlist(begin,eof);

		
		for(int j=0;j<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j++){
			linearParams[j]=paramlist[2+j];
			nonLinearParams[j]=0;
		}
		TNupdateParameters(psr,0,linearParams,errorvec, nonLinearParams);
		
		std::stringstream txtstream;
		txtstream.precision(std::numeric_limits<double>::digits10);//override the default


		txtstream << paramlist[0];
		txtstream << "\t";
		txtstream << paramlist[1];
		txtstream << "\t";
		for(int j=0;j<((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j++){
			if(j==0){
				txtstream << nonLinearParams[j];
				txtstream << "\t";
			}
			else{
				if(((MNStruct *)context)->LDpriors[j][1] != 0){
					txtstream << nonLinearParams[j]/((MNStruct *)context)->LDpriors[j][1];
					txtstream << "\t";
				}
				else{
					txtstream << 0;
                    txtstream << "\t";
				}
// 				printf("%i %i %g %g %g \n",i,j,linearParams[j],nonLinearParams[j],(double)((MNStruct *)context)->LDpriors[j-1][1]);
			}
		}
		for(int j=((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps;j<ndim;j++){

				txtstream << paramlist[j+2];
				txtstream << "\t";


		}
		newtxtfile << txtstream.str();
		newtxtfile << "\n";

	}
	summaryfile.close();
	newtxtfile.close();



}


void writeCov(std::vector<double> Cube, int &ndim,void *context, std::string longname, int outputoption)
{

	clock_t startClock,endClock;
	long double LDparams[ndim];
	double *EFAC;
	double EQUAD, redamp, redalpha;
	
			  for(int i =0; i< 4*ndim; i++){
		  	printf("bparm: %i %g \n ", i, Cube[i]);
			}	
	
	
	int pcount=outputoption*ndim;


	for(int p=0;p<((MNStruct *)context)->numFitTiming;p++){
		pcount++;
	}
	for(int p=0;p<((MNStruct *)context)->numFitJumps;p++){
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
			printf("EF: %i %g \n",p,EFAC[p]);
			pcount++;
			
		}
	}				

	if(((MNStruct *)context)->numFitEQUAD == 0){
		EQUAD=0;
// 		printf("EQUAD: %g \n",EQUAD);
	}
	else{
		
		EQUAD=pow(10.0,2*Cube[pcount]);
		printf("EQ: %g \n",Cube[pcount]);
		pcount++;
		

	}


	redamp=Cube[pcount];
	pcount++;
	redalpha=Cube[pcount];
	pcount++;

	printf("red: %g %g \n",redamp, redalpha);
	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=redalpha;
	double gwamp=pow(10.0,redamp);
	double gwampsquared=gwamp*gwamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-modelalpha)))/(pow(flo,(modelalpha-1)));

	double timdiff=0;

	double covconst=gsl_sf_gamma(1-modelalpha)*sin(0.5*M_PI*modelalpha);

	
	double **CovMatrix = new double*[((MNStruct *)context)->pulse->nobs]; for(int o1=0;o1<((MNStruct *)context)->pulse->nobs;o1++)CovMatrix[o1]=new double[((MNStruct *)context)->pulse->nobs];

	for(int o1=0;o1<((MNStruct *)context)->pulse->nobs; o1++){

		for(int o2=0;o2<((MNStruct *)context)->pulse->nobs; o2++){
			timdiff=((MNStruct *)context)->pulse->obsn[o1].bat-((MNStruct *)context)->pulse->obsn[o2].bat;	
			double tau=2.0*M_PI*fabs(timdiff);
			double covsum=0;

			for(int k=0; k <=4; k++){
				covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-modelalpha));

			}

			CovMatrix[o1][o2]=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);
// 			printf("%i %i %g %g %g\n",o1,o2,CovMatrix[o1][o2],fabs(timdiff),covsum);

			if(o1==o2){
				CovMatrix[o1][o2] += pow(((((MNStruct *)context)->pulse->obsn[o1].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)context)->sysFlags[o1]],2) + EQUAD;
			}

		}
	}
	
	
	FILE *pFile;
	std::string covname = longname+"covMatrix.txt";
	pFile = fopen(covname.c_str(), "w+");
	
	for(int i=0; i<((MNStruct *)context)->pulse->nobs; i++) {

    	for(int j=0; j<((MNStruct *)context)->pulse->nobs; j++) {
     		 fprintf(pFile, "%22.20e", CovMatrix[i][j]);
     		 fprintf(pFile, " ");
     	}
        fprintf(pFile, "\n");
     }
     	
     fclose(pFile);	
     

	std::string maxstoc = longname+"maxStoc.txt";
	pFile = fopen(maxstoc.c_str(), "w+");
	
	for(int p=0;p< ((MNStruct *)context)->numFitEFAC; p++){
			fprintf(pFile, "Max EFAC %i %22.20e \n",p, EFAC[p]);
		}

	fprintf(pFile, "Max EQUAD 0 %22.20e \n",0.50*log10(EQUAD));
	fprintf(pFile, "Max Red Amp %22.20e \n",redamp);
	fprintf(pFile, "Max Red Index %22.20e \n",redalpha);

     fclose(pFile);	
}	 
     		 
