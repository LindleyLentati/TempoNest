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

#ifdef HAVE_PSRCHIVE
// psrchive stuff
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/FaradayRotation.h"
#include "Pulsar/PolnProfileStats.h"
#include "T2Observatory.h"
using namespace Pulsar;
#endif

int UtWrap(int kX, int const kLowerBound, int const kUpperBound)
{
    int range_size = kUpperBound - kLowerBound + 1;

    if (kX < kLowerBound)
        kX += range_size * ((kLowerBound - kX) / range_size + 1);

    return kLowerBound + (kX - kLowerBound) % range_size;
}

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

	printf("Getting Errors \n");

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

	//void *tempcontext;
	//double *DerivedParams=new double[1];
	//int np = 1;
	//double like = NewLRedMarginLogLike(ndim, TempCube, np, DerivedParams, tempcontext);

}

void read_ddgr(std::string longname, int ndim, int sampler,  void *context){

        int number_of_lines = 0;

        std::ifstream checkfile;
        std::string checkname;
	if(sampler == 0){
                checkname = longname+"post_equal_weights.dat";
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
                fname = longname+"post_equal_weights.dat";
        }
        if(sampler == 1){
                fname = longname+"_phys_live.txt";
        }


	summaryfile.open(fname.c_str());

	double *TempCube = new double[ndim];


//	printf("Getting ML \n");
	double maxlike = -1.0*pow(10.0,10);
	for(int i=0;i<number_of_lines;i++){

		std::string line;
		getline(summaryfile,line);
		std::istringstream myStream( line );
                std::istream_iterator< double > begin(myStream),eof;
                std::vector<double> paramlist(begin,eof);

		double like = paramlist[ndim];

		 for(int i = 0; i < ndim; i++){
			TempCube[i] = paramlist[i];

		 }
		long double LDparams[((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps];
		int fitcount = 0;
		for(int p=0;p< ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps; p++){
			if(((MNStruct *)context)->Dpriors[p][1] != ((MNStruct *)context)->Dpriors[p][0]){

				double val = 0;
				if((((MNStruct *)context)->LDpriors[p][3]) == 0){
					val = TempCube[fitcount];
				}
				if((((MNStruct *)context)->LDpriors[p][3]) == 1){
					val = pow(10.0,TempCube[fitcount]);
				}
				if((((MNStruct *)context)->LDpriors[p][3]) == 2){
					val = pow(10.0,TempCube[fitcount]);
				}
				LDparams[p]=val*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
				
				if(((MNStruct *)context)->TempoFitNums[p][0] == param_sini && ((MNStruct *)context)->usecosiprior == 1){
						val = TempCube[fitcount];
						LDparams[p] = sqrt(1.0 - val*val);
				}

				fitcount++;

			}
			else if(((MNStruct *)context)->Dpriors[p][1] == ((MNStruct *)context)->Dpriors[p][0]){
				LDparams[p]=((MNStruct *)context)->Dpriors[p][0]*(((MNStruct *)context)->LDpriors[p][1]) + (((MNStruct *)context)->LDpriors[p][0]);
			}


		}
		int pcount=0;
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
		DDGRmodel(((MNStruct *)context)->pulse,  0, 0, -2);

		printf("PBDot %.20Lg  OMDot %.20Lg Gamma %.20Lg \n", ((MNStruct *)context)->pulse->param[param_pbdot].val[0], ((MNStruct *)context)->pulse->param[param_omdot].val[0], ((MNStruct *)context)->pulse->param[param_gamma].val[0]);


	}
	summaryfile.close();

	//void *tempcontext;
	//double *DerivedParams=new double[1];
	//int np = 1;
	//double like = NewLRedMarginLogLike(ndim, TempCube, np, DerivedParams, tempcontext);

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
//	read_ddgr( longname, ndim, ((MNStruct *)context)->sampler,  context);
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

			if(((MNStruct *)context)->TempoFitNums[pcount][0] == param_sini && ((MNStruct *)context)->usecosiprior == 1){
				val = paramarray[fitcount][2];
				LDP[j]  = sqrt(1.0 - val*val);
			}


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
	printf("text output\n");	
	double Evidence=0;
	TNtextOutput(((MNStruct *)context)->pulse, 1, 0, Tempo2Fit,  context,incRED,ndims,paramlist, Evidence, doTimeMargin, doJumpMargin, doLinear, longname, paramarray);

	
	
	
	printf("finished output \n");

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

void getPhysDVector(void *context, double **TNDM, int Nobs, int *TimingGradientSigns){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//In order to get design matrix in physical units, need to use update parameters function.  Current fit values/errors are backed up in TempoPriors. Set value/error to 0,1 then call update params.
	//The new error will be the scaling required for the design matrix, then set values/errors back to old values.

	//long double newF0 = 3.5*(((MNStruct *)context)->LDpriors[1][1]) + (((MNStruct *)context)->LDpriors[1][0]);
	//((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[1][0]].val[((MNStruct *)context)->TempoFitNums[1][1]] = newF0;

	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,0);       /* Form residuals */



	int numToFit = ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	double *newVals = new double[numToFit]();
	double *newErrs = new double[numToFit]();
	double *oldErrs = new double[numToFit]();

	int fitcount=0;
	oldErrs[fitcount] = 0;
	fitcount++;
	for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<((MNStruct *)context)->pulse->param[p].aSize;k++){
			if(((MNStruct *)context)->pulse->param[p].fitFlag[k] == 1 && p != param_dmmodel){
				//printf("Old Error for: %s %g\n",((MNStruct *)context)->pulse->param[p].shortlabel[k], (double)((MNStruct *)context)->pulse->param[p].err[k]);
				oldErrs[fitcount] = (double)((MNStruct *)context)->pulse->param[p].err[k];
				newErrs[fitcount] = 1;
				fitcount++;
	    		}
		}
	}

	updateParameters(((MNStruct *)context)->pulse, 0,newVals,newErrs);



	for (int p=0;p<numToFit;p++) {
		TimingGradientSigns[p] = 1;
	}

	fitcount=0;
	newVals[fitcount] = 1;
	fitcount++;
	for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<((MNStruct *)context)->pulse->param[p].aSize;k++){
			if(((MNStruct *)context)->pulse->param[p].fitFlag[k] == 1 && p != param_dmmodel){
				//printf("Scaling factor for: %s %g\n",((MNStruct *)context)->pulse->param[p].shortlabel[k], (double)((MNStruct *)context)->pulse->param[p].err[k]);
				newVals[fitcount] =  (double)((MNStruct *)context)->pulse->param[p].err[k];
				//printf("Resetting to old error: %g\n",oldErrs[fitcount]);
				((MNStruct *)context)->pulse->param[p].err[k] = oldErrs[fitcount];
				fitcount++;
	    		}
		}
	}

	for(int i = 0; i <  ((MNStruct *)context)->numFitJumps; i++){
		newVals[((MNStruct *)context)->numFitTiming+i] = 1;
	}
	
	int DMparamnum=3;
	double averageDMDelay = 0;
	
	for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
		FITfuncs(((MNStruct *)context)->pulse->obsn[i].bat - ((MNStruct *)context)->pulse->param[param_pepoch].val[0], pdParamDeriv, numToFit, ((MNStruct *)context)->pulse, i,0);
		for(int j=0; j<numToFit; j++) {
			double PriorFac = (double)((MNStruct *)context)->LDpriors[j][1];
	//		printf("CDML: %i %i %i %g %g %g\n", i,j,numToFit,pdParamDeriv[j], pdParamDeriv[j]/newVals[j], (pdParamDeriv[j]/newVals[j])*PriorFac);
			TNDM[i][j]= pdParamDeriv[j];//(pdParamDeriv[j]/newVals[j])*PriorFac;
//			if(j< 20){printf("%i %i %g \n", i, j, TNDM[i][j]);}
			if(j==DMparamnum){averageDMDelay += TNDM[i][j];}
		} 
	} 

	averageDMDelay/=((MNStruct *)context)->pulse->nobs;
	printf("AverageDMDelay is : %g\n", averageDMDelay);
//	if(((MNStruct *)context)->doLinear==1){
//		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
//			TNDM[i][DMparamnum] -= averageDMDelay;
//		}
//	}

	if(((MNStruct *)context)->doLinear==0){
		for(int p= 1; p < ((MNStruct *)context)->numFitTiming; p++){
			double zeroval=	TNDM[0][p];
			for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
				TNDM[i][p] = TNDM[i][p]-zeroval;
			}
		}
	}

	



	double *oldRes = new double[((MNStruct *)context)->pulse->nobs];
	for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {	
		oldRes[i] = (double)((MNStruct *)context)->pulse->obsn[i].residual;
	}

	for(int p= 1; p < ((MNStruct *)context)->numFitTiming; p++){

		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,0);       /* Form residuals */


		double nobatsSTD1 = 0;
		double nobatsSTD2 = 0;

		long double newVal = ((MNStruct *)context)->LDpriors[p][1] + ((MNStruct *)context)->LDpriors[p][0];
		long double oldVal = ((MNStruct *)context)->LDpriors[p][0];
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = newVal;


//		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,0);       /* Form residuals */
		double sign = 0;
		double DMSub = 0;
		if(strcasecmp(((MNStruct *)context)->pulse[0].param[p].shortlabel[0],"DM")!=0 || strcasecmp(((MNStruct *)context)->pulse[0].param[p].shortlabel[0],"DMX_001")!=0){
			//printf("This is a DM parameter %i\n",p);
			//DMSub = averageDMDelay;
		}
		double zeroval= TNDM[0][p];

		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {	
			double PriorFac = (double)((MNStruct *)context)->LDpriors[p][1];
			TNDM[i][p] = ((TNDM[i][p]-zeroval)/newVals[p])*PriorFac;
			//if(p < 20){printf("Vecs: %i %i %.16g %.16g %.16g %.16g\n", p, i, (double)((MNStruct *)context)->pulse->obsn[i].bat, (double)((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i],TNDM[i][p], (TNDM[i][p]-(TNDM[0][p])));}
			sign += double(((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i])*(TNDM[i][p]-(TNDM[0][p]));

			nobatsSTD1 += (((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i])*(((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i]);
			nobatsSTD2 += TNDM[i][p]*TNDM[i][p];
			//TNDM[i][p] = (double)((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i];
		}

                double withbatsSTD1 = 0;
                double withbatsSTD2 = 0;

		fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars,0);       /* Form residuals */

		double average = 0;
		double sign2 = 0;
		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {	
			if(i<5){printf("Vecs: %i %i %.16g %.16g %.16g\n", p, i, (double)((MNStruct *)context)->pulse->obsn[i].bat, (double)((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i],TNDM[i][p]);}
			sign2 += double(((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i])*(TNDM[i][p]-(TNDM[0][p]));
			withbatsSTD1 += (((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i])*(((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i]);
			withbatsSTD2 += TNDM[i][p]*TNDM[i][p];
			average += TNDM[i][p];
			//TNDM[i][p] = (double)((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i];
		}

		printf("BATS Ratio: %i %s %g %g \n", p, ((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].shortlabel[0], nobatsSTD1/withbatsSTD1, withbatsSTD2/nobatsSTD2);


		printf("signs? %i %g %g %g\n", p, sign, sign2, average/((MNStruct *)context)->pulse->nobs);

//		printf("Sign? %i %g %g\n", p, sign, sign/fabs(sign));
		TimingGradientSigns[p] = int(sign2/fabs(sign2));
		((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[p][0]].val[((MNStruct *)context)->TempoFitNums[p][1]] = oldVal;

		if(fabs(1-(nobatsSTD1/withbatsSTD1)) > 0.01 && ((MNStruct *)context)->doLinear==0){

			printf("param %i changes bats alot, using residuals to form gradient \n", p);
			for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
				TNDM[i][p] = (double)((MNStruct *)context)->pulse->obsn[i].residual-oldRes[i];
			}
		}
		else{
			printf("param %i only changes residuals, using design matrix for gradient \n",p);
			if(sign2 < 0){
				for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {	
					TNDM[i][p] *=-1;
				}
			}
		}	
		
	}

	for(int p = 0; p < ((MNStruct *)context)->numFitJumps; p++){


		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
//			printf("Jumps in TNDM: %i %i %g %g \n", p, i,  (double)((MNStruct *)context)->LDpriors[p+((MNStruct *)context)->numFitTiming][1], TNDM[i][p+((MNStruct *)context)->numFitTiming]);
			TNDM[i][p+((MNStruct *)context)->numFitTiming] = ((MNStruct *)context)->LDpriors[p+((MNStruct *)context)->numFitTiming][1]*TNDM[i][p+((MNStruct *)context)->numFitTiming];
		}
	

	}
	printf("do linear is: %i \n", ((MNStruct *)context)->doLinear);
	if(((MNStruct *)context)->doLinear==2){

		double *TempTNDM = new double[((MNStruct *)context)->pulse->nobs*numToFit]();
		double *NU = new double[((MNStruct *)context)->pulse->nobs*(numToFit)]();
		double *UU = new double[(numToFit)*(numToFit)];

		for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
			for(int k=0;k<numToFit;k++){
			
				if(k==0){TNDM[j][k] = ((MNStruct *)context)->ReferencePeriod;}
				TempTNDM[j + k*((MNStruct *)context)->pulse->nobs] = TNDM[j][k];
			}
		}

	//	vector_dgesvd2(TempTNDM, ((MNStruct *)context)->pulse->nobs, numToFit-1, S, U, VT);
		vector_dgesvd(TempTNDM,((MNStruct *)context)->pulse->nobs, numToFit);

		for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
			for(int k=0;k<numToFit;k++){
				NU[j + k*((MNStruct *)context)->pulse->nobs] = TempTNDM[j + k*((MNStruct *)context)->pulse->nobs]/pow(((MNStruct *)context)->pulse->obsn[j].toaErr*pow(10.0,-6),2);
			}
		}

		double *dNU= new double[numToFit];
		double *ML= new double[numToFit];
		double *residuals = new double[((MNStruct *)context)->pulse->nobs];

		 for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
			residuals[j] = ((MNStruct *)context)->pulse->obsn[j].residual;
		}

		vector_dgemm(NU, TempTNDM, UU, ((MNStruct *)context)->pulse->nobs, numToFit, ((MNStruct *)context)->pulse->nobs, numToFit, 'T', 'N');
		vector_dgemv(NU, residuals, dNU, ((MNStruct *)context)->pulse->nobs, numToFit, 'T');

		double Margindet = 0;
		int info = 0;
	        vector_dpotrfInfo(UU, numToFit, Margindet, info);
                vector_dpotri(UU,numToFit);

		vector_dgemv(UU, dNU, ML, numToFit, numToFit, 'N');
		vector_dgemv(TempTNDM, ML, residuals, ((MNStruct *)context)->pulse->nobs, numToFit, 'N');

		for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
			 ((MNStruct *)context)->pulse->obsn[j].residual -= residuals[j];
		}	

		double *ErrScale = new double[numToFit];

		for(int k=0;k<numToFit;k++){
			printf("ErrScale: %i %g  %g %g \n", k, sqrt(UU[k + k*(numToFit)]), ML[k], ML[k]/sqrt(UU[k + k*(numToFit)]));
			ErrScale[k] = sqrt(UU[k + k*(numToFit)]);
		}

		for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
			for(int k=0;k<numToFit;k++){
				if(k==0)printf("Lin2: %i %i %g %g %g\n", j, k, TNDM[j][k], TempTNDM[j + k*((MNStruct *)context)->pulse->nobs], TempTNDM[j + k*((MNStruct *)context)->pulse->nobs]*ErrScale[k]);
				TNDM[j][k]=TempTNDM[j + k*((MNStruct *)context)->pulse->nobs]*ErrScale[k];
			}
		}

		for(int k=1;k<numToFit;k++){
			double zeroval = TNDM[0][k];
			for(int j=0;j<((MNStruct *)context)->pulse->nobs;j++){
//				TNDM[j][k] -= zeroval;
			}
		}

		delete[]TempTNDM;
		delete[] ErrScale;
	
		
	}





	for(int j=0; j<numToFit; j++) {
		double dsum = 0;
		double gsum = 0;
		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {	
			dsum += TNDM[i][j]*TNDM[i][j];
			gsum+= TNDM[i][j]*(((MNStruct *)context)->pulse->obsn[i].residual)/(((MNStruct *)context)->pulse->obsn[i].toaErr*pow(10.0, -6))/(((MNStruct *)context)->pulse->obsn[i].toaErr*pow(10.0, -6));
		}
		printf("TOA Hess and grad: %i %g %g \n", j, dsum, gsum);
	}

	delete[] newVals;
	delete[] newErrs;
	delete[] oldErrs;
	delete[] oldRes;
}


void UpdatePhysDVector(void *context, double **TNDM, int Nobs){
	
	double pdParamDeriv[MAX_PARAMS], dMultiplication;


	//In order to get design matrix in physical units, need to use update parameters function.  Current fit values/errors are backed up in TempoPriors. Set value/error to 0,1 then call update params.
	//The new error will be the scaling required for the design matrix, then set values/errors back to old values.

	//long double newF0 = 3.5*(((MNStruct *)context)->LDpriors[1][1]) + (((MNStruct *)context)->LDpriors[1][0]);
	//((MNStruct *)context)->pulse->param[((MNStruct *)context)->TempoFitNums[1][0]].val[((MNStruct *)context)->TempoFitNums[1][1]] = newF0;



	int numToFit = ((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;
	double *newVals = new double[numToFit]();
	double *newErrs = new double[numToFit]();

	for(int p= 0; p < numToFit; p++){
		newErrs[p] = 1;
	}

	updateParameters(((MNStruct *)context)->pulse, 0,newVals,newErrs);


	int fitcount=1;

	for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<((MNStruct *)context)->pulse->param[p].aSize;k++){
			if(((MNStruct *)context)->pulse->param[p].fitFlag[k] == 1 && p != param_dmmodel){
				//printf("Scaling factor for: %s %g\n",((MNStruct *)context)->pulse->param[p].shortlabel[k], (double)((MNStruct *)context)->pulse->param[p].err[k]);
				newVals[fitcount] =  (double)((MNStruct *)context)->pulse->param[p].err[k];
				fitcount++;
	    		}
		}
	}

	
	for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
		FITfuncs(((MNStruct *)context)->pulse->obsn[i].bat - ((MNStruct *)context)->pulse->param[param_pepoch].val[0], pdParamDeriv, numToFit, ((MNStruct *)context)->pulse, i,0);
		for(int j=1; j<((MNStruct *)context)->numFitTiming; j++) {
			TNDM[i][j]= pdParamDeriv[j];
		} 
	} 

	for(int p= 1; p < ((MNStruct *)context)->numFitTiming; p++){
		double zeroval=	TNDM[0][p];
		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {
			TNDM[i][p] = TNDM[i][p]-zeroval;
		}
	}




	for(int p= 1; p < ((MNStruct *)context)->numFitTiming; p++){

		int sign = ((MNStruct *)context)->TimingGradientSigns[p];
		double PriorFac = (double)((MNStruct *)context)->LDpriors[p][1];
		PriorFac *= sign/newVals[p];
		//printf("updating %i %i %g %g \n", p, sign, newVals[p], PriorFac);

		for(int i=0; i < ((MNStruct *)context)->pulse->nobs; i++) {	

			TNDM[i][p] = TNDM[i][p]*PriorFac;

		}
	}

	delete[] newVals;
	delete[] newErrs;

}



void getCustomDVectorLike(void *context, double *TNDM, int Nobs, int TimeToMargin, int TotalSize){
	
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

		//printf("compare: %i %i \n", numToMargin, TimeToMargin);
	
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


/*

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
	for(int i=0;i<pulse->nobs;i++){
		for(int j=0;j<Gsize;j++){
			printf("%i %i %g \n ",i,j,staticGMatrix[i][j]);
		}
	}
	printf("%i %i \n",pulse->nobs,Gsize);

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
*/

void StoreTMatrix(double *TotalMatrix, void *context){

	int totalsize = ((MNStruct *)context)->totalsize;

	for(int i = 0; i < ((MNStruct *)context)->pulse->nobs*totalsize; i++){
		TotalMatrix[i] = 0;
	}
		
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	int TimetoMargin=((MNStruct *)context)->TimetoMargin;

	getCustomDVectorLike(context, TotalMatrix, ((MNStruct *)context)->pulse->nobs, TimetoMargin, totalsize);
	vector_dgesvd(TotalMatrix,((MNStruct *)context)->pulse->nobs, TimetoMargin);
		



//////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////Set up Coefficients///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  
	double maxtspan=((MNStruct *)context)->Tspan;


	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);
	int FitBandCoeff=2*(((MNStruct *)context)->numFitBandNoiseCoeff);
	int FitGroupNoiseCoeff = 2*((MNStruct *)context)->numFitGroupNoiseCoeff;

	int totCoeff=((MNStruct *)context)->totCoeff;



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Red Noise///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	double *freqs = new double[totCoeff];
	double *DMVec=new double[((MNStruct *)context)->pulse->nobs];

	double DMKappa = 2.410*std::pow(10.0,-16);
	int startpos=0;

	if(((MNStruct *)context)->incRED > 0 || ((MNStruct *)context)->incGWB == 1){
		for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)context)->incFloatRed ; i++){
			
			freqs[startpos+i]=(double)((MNStruct *)context)->sampleFreq[i]/maxtspan;
			freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

		}
	

		for(int i=0;i<FitRedCoeff/2;i++){
			for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
				double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
				TotalMatrix[k + (i+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs]=cos(2*M_PI*freqs[i]*time);
				TotalMatrix[k + (i+FitRedCoeff/2+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs] = sin(2*M_PI*freqs[i]*time);


			}
		}
			
			    
		startpos += FitRedCoeff;

	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////DM Variations////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////// 


	if(((MNStruct *)context)->incDM > 0){

                for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*std::pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                }

        	for(int i=0;i<FitDMCoeff/2;i++){

			freqs[startpos+i]=((MNStruct *)context)->sampleFreq[startpos/2 - ((MNStruct *)context)->incFloatRed+i]/maxtspan;
			freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];


                	for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
                        	double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;

				TotalMatrix[k + (i+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs]=cos(2*M_PI*freqs[startpos+i]*time)*DMVec[k];
				TotalMatrix[k + (i+FitDMCoeff/2+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs] = sin(2*M_PI*freqs[startpos+i]*time)*DMVec[k];

        	        }
		}

		startpos += FitDMCoeff;
	} 


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Band DM/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


        if(((MNStruct *)context)->incBandNoise > 0){

                if(((MNStruct *)context)->incDM == 0){
                        for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
                                DMVec[o]=1.0/(DMKappa*std::pow((double)((MNStruct *)context)->pulse->obsn[o].freqSSB,2));
                        }
                }

		for(int b = 0; b < ((MNStruct *)context)->incBandNoise; b++){

			double startfreq = ((MNStruct *)context)->FitForBand[b][0];
			double stopfreq = ((MNStruct *)context)->FitForBand[b][1];
			double BandScale = ((MNStruct *)context)->FitForBand[b][2];
			int BandPriorType = ((MNStruct *)context)->FitForBand[b][3];


			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;

			for (int i=0; i<FitBandCoeff/2; i++){

				freqs[startpos+i]=((double)(i+1))/maxtspan;
				freqs[startpos+i+FitBandCoeff/2]=freqs[startpos+i];
				
			}
			

			for(int i=0;i<FitBandCoeff/2;i++){
				for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
					if(((MNStruct *)context)->pulse->obsn[k].freq > startfreq && ((MNStruct *)context)->pulse->obsn[k].freq < stopfreq){
						double time=(double)((MNStruct *)context)->pulse->obsn[k].bat;
						TotalMatrix[k + (i+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs]=cos(2*M_PI*freqs[startpos+i]*time);
                                                TotalMatrix[k + (i+TimetoMargin+startpos+FitBandCoeff/2)*((MNStruct *)context)->pulse->nobs]=sin(2*M_PI*freqs[startpos+i]*time);
					}
					else{	
                                                TotalMatrix[k + (i+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs]=0;
                                                TotalMatrix[k + (i+TimetoMargin+startpos+FitBandCoeff/2)*((MNStruct *)context)->pulse->nobs]=0;

					}


				}
			}


			startpos += FitBandCoeff;
		}

    	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add Group Noise/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

        if(((MNStruct *)context)->incGroupNoise > 0){

		for(int g = 0; g < ((MNStruct *)context)->incGroupNoise; g++){

			startpos=startpos+FitGroupNoiseCoeff;
		}

    }

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add ECORR Coeffs////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	if(((MNStruct *)context)->incNGJitter > 0){

		//printf("Calling\n");
		double **NGJitterMatrix;
		int NumNGEpochs = ((MNStruct *)context)->numNGJitterEpochs;
		//getNGJitterMatrixEpochs(psr, NumNGEpochs);

		NGJitterMatrix = new double*[((MNStruct *)context)->pulse->nobs];
		for(int i =0; i < ((MNStruct *)context)->pulse->nobs; i++){
			NGJitterMatrix[i] = new double[NumNGEpochs];
			for(int j =0; j < NumNGEpochs;  j++){
				NGJitterMatrix[i][j] = 0;
			}
		}
		int *NGJitterSysFlags = new int[((MNStruct *)context)->pulse->nobs];
		

		for (int i=0;i<((MNStruct *)context)->pulse->nobs;i++){
			for (int j=0;j<((MNStruct *)context)->pulse->obsn[i].nFlags;j++){
				for (int k=0;k<((MNStruct *)context)->pulse->nTNECORR;k++){
					if (strcmp(((MNStruct *)context)->pulse->obsn[i].flagID[j], ((MNStruct *)context)->pulse->TNECORRFlagID[k])==0){
						if (strcmp(((MNStruct *)context)->pulse->obsn[i].flagVal[j],((MNStruct *)context)->pulse->TNECORRFlagVal[k])==0){
							NGJitterSysFlags[i] = k;
							//printf("NGFLag? %i %i %s \n", i, k, psr->TNECORRFlagVal[k]);
						}
					}
				}
			}
		}


		getNGJitterMatrix(((MNStruct *)context)->pulse, NGJitterMatrix, NumNGEpochs);

		int *NGJitterEpochFlags = new int[NumNGEpochs];

		for(int i =0; i < NumNGEpochs; i++){
			for(int j=0; j < ((MNStruct *)context)->pulse->nobs; j++){
				if(NGJitterMatrix[j][i] != 0) {
					NGJitterEpochFlags[i]=NGJitterSysFlags[j];
				}
			}
		}

		((MNStruct *)context)->NGJitterSysFlags=NGJitterEpochFlags;

		for(int k=0;k<((MNStruct *)context)->pulse->nobs;k++){
			for(int i=0; i < ((MNStruct *)context)->numNGJitterEpochs; i++){
				TotalMatrix[k + (i+TimetoMargin+startpos)*((MNStruct *)context)->pulse->nobs] = NGJitterMatrix[k][i];///sqrt((((MNStruct *)context)->TobsInfo[k]/3600.0));
			}
		}

		delete[] NGJitterSysFlags;
		for(int j =0; j < ((MNStruct *)context)->pulse->nobs;  j++){
			delete[] NGJitterMatrix[j];
		}
		delete[]NGJitterMatrix;


	}





	delete[] DMVec;
	delete[] freqs;
	
	

	
}


void getArraySizeInfo(void *context){


	int TimetoMargin=0;
	for(int i =0; i < ((MNStruct *)context)->numFitTiming+((MNStruct *)context)->numFitJumps; i++){
		if(((MNStruct *)context)->LDpriors[i][2]==1)TimetoMargin++;
	}


//////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////Set up Coefficients///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  



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


	int FitRedCoeff=2*(((MNStruct *)context)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)context)->numFitDMCoeff);
	int FitBandNoiseCoeff=2*(((MNStruct *)context)->numFitBandNoiseCoeff);
	int FitGroupNoiseCoeff = 2*((MNStruct *)context)->numFitGroupNoiseCoeff;

	int NumNGEpochs = 0;
	if(((MNStruct *)context)->incNGJitter > 0){
		getNGJitterMatrixEpochs(((MNStruct *)context)->pulse, NumNGEpochs);
	}
	((MNStruct *)context)->numNGJitterEpochs = NumNGEpochs;

	int totCoeff=0;
	if(((MNStruct *)context)->incRED != 0 || ((MNStruct *)context)->incGWB == 1)totCoeff+=FitRedCoeff;
	if(((MNStruct *)context)->incDM != 0)totCoeff+=FitDMCoeff;
	if(((MNStruct *)context)->incBandNoise > 0)totCoeff+= ((MNStruct *)context)->incBandNoise*FitBandNoiseCoeff;
	if(((MNStruct *)context)->incNGJitter >0)totCoeff+=((MNStruct *)context)->numNGJitterEpochs;
	if(((MNStruct *)context)->incGroupNoise > 0)totCoeff += ((MNStruct *)context)->incGroupNoise*FitGroupNoiseCoeff;


	int MarginRedShapelets = ((MNStruct *)context)->MarginRedShapeCoeff;
	int totalredshapecoeff = 0;


	if(((MNStruct *)context)->incRedShapeEvent != 0){
		if(MarginRedShapelets == 1){
			totalredshapecoeff = ((MNStruct *)context)->numRedShapeCoeff*((MNStruct *)context)->incRedShapeEvent;
		}
	}

	int totalsize=TimetoMargin+totCoeff+totalredshapecoeff;

        ((MNStruct *)context)->Tspan = maxtspan;
        ((MNStruct *)context)->TimetoMargin = TimetoMargin;
        ((MNStruct *)context)->totCoeff= totCoeff;
        ((MNStruct *)context)->totRedShapeCoeff = totalredshapecoeff;
        ((MNStruct *)context)->totalsize = totalsize;


//	printf("TimetoMargin %i, totCoeff %i, totalredshapecoeff %i, totalsize %i \n", TimetoMargin, totCoeff, totalredshapecoeff, totalsize);

}

void getNumTempFreqs(int &NumFreqs, void *globalcontext, double TemplateChanWidth){


	double chanwidth  = TemplateChanWidth;
	double minfreq=pow(10.0, 100);

	for(int p=0;p< ((MNStruct *)globalcontext)->pulse->nobs; p++){
		double freq = floor(((MNStruct *)globalcontext)->pulse->obsn[p].freq);
		if(freq < minfreq){minfreq=freq;}
	}

	if(NumFreqs == 0){

		printf("getting freqs \n");
		std::vector <double> freqs;
		for(int p=0;p< ((MNStruct *)globalcontext)->pulse->nobs; p++){
			double freq = floor(((MNStruct *)globalcontext)->pulse->obsn[p].freq);
			freq = floor((freq-minfreq)/chanwidth)*chanwidth + minfreq;
			//printf("freq: %i %g %g %g \n", p, chanwidth, freq, ((MNStruct *)globalcontext)->pulse->obsn[p].freq);
			if(std::find(freqs.begin(),freqs.end(), freq) != freqs.end()) {
	
			}
			else{
				freqs.push_back(freq);
			}
		}

		NumFreqs = int(freqs.size());
		printf("have %i unique frequencies \n", NumFreqs);
	}
	else{
		std::vector <double> freqs;
		for(int p=0;p< ((MNStruct *)globalcontext)->pulse->nobs; p++){
			double freq = floor(((MNStruct *)globalcontext)->pulse->obsn[p].freq);
			freq = floor((freq-minfreq)/chanwidth)*chanwidth + minfreq;
	//		printf("freq: %i %g %g \n", p, freq, ((MNStruct *)globalcontext)->pulse->obsn[p].freq);
			if(std::find(freqs.begin(),freqs.end(), freq) != freqs.end()) {
	
			}
			else{
				((MNStruct *)globalcontext)->TemplateFreqs[int(freqs.size())] = freq;
				freqs.push_back(freq);
				
			}
		}

		NumFreqs = int(freqs.size());
		printf("have set %i unique frequencies \n", NumFreqs);

	}

}


void Tscrunch(void *globalcontext, double TemplateChanWidth){


	double chanwidth  = TemplateChanWidth;
	double minfreq=pow(10.0, 100);

	for(int p=0;p< ((MNStruct *)globalcontext)->pulse->nobs; p++){
		double freq = floor(((MNStruct *)globalcontext)->pulse->obsn[p].freq);
		if(freq < minfreq){minfreq=freq;}
	}


        long double LDparams[((MNStruct *)globalcontext)->numFitTiming + ((MNStruct *)globalcontext)->numFitJumps];
	int pcount = 0;
	int fitcount=0;

	int debug = 0;

	int numfreqs = ((MNStruct *)globalcontext)->numTempFreqs;
	double *freqs = ((MNStruct *)globalcontext)->TemplateFreqs;
	double *fweights = new double[numfreqs]();

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	long double phase = 0.0*((MNStruct *)globalcontext)->ReferencePeriod/SECDAY;


	for(int p=0;p< ((MNStruct *)globalcontext)->pulse->nobs; p++){
		((MNStruct *)globalcontext)->pulse->obsn[p].sat = ((MNStruct *)globalcontext)->pulse->obsn[p].origsat-phase;
	}



	fastformBatsAll(((MNStruct *)globalcontext)->pulse,((MNStruct *)globalcontext)->numberpulsars);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)globalcontext)->pulse,((MNStruct *)globalcontext)->numberpulsars,0);       /* Form residuals */
	



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	 long double **ProfileBats=new long double*[((MNStruct *)globalcontext)->pulse->nobs];
	 long double *ModelBats = new long double[((MNStruct *)globalcontext)->pulse->nobs];
	 for(int i = 0; i < ((MNStruct *)globalcontext)->pulse->nobs; i++){

	      int Nbin  = (int)((MNStruct *)globalcontext)->ProfileInfo[i][1];
	      ProfileBats[i] = new long double[Nbin];
	      for(int j = 0; j < Nbin; j++){
		    ProfileBats[i][j] = ((MNStruct *)globalcontext)->ProfileData[i][j][0] + ((MNStruct *)globalcontext)->pulse->obsn[i].batCorr;
		
	      }
	      
	      
	      ModelBats[i] = ((MNStruct *)globalcontext)->ProfileInfo[i][5]+((MNStruct *)globalcontext)->pulse->obsn[i].batCorr - phase - ((MNStruct *)globalcontext)->pulse->obsn[i].residual/SECDAY ;

	 }

	int GlobalNBins = (int)((MNStruct *)globalcontext)->ProfileInfo[0][1];
	int interpFactor = 1;
	int interpNBins = GlobalNBins*interpFactor;

	double *TScrunched = new double[numfreqs*GlobalNBins]();
	double **interpTScrunched = new double*[interpFactor];
	for(int i = 0; i < interpFactor; i++){interpTScrunched[i]= new double[numfreqs*GlobalNBins]();}
	double **interpWeights = new double*[interpFactor];
	for(int i = 0; i < interpFactor; i++){interpWeights[i] = new double[numfreqs*GlobalNBins]();}

	int t = 0;
	int nTOA = 0;

	double WeightSum = 0.0;
	
	for(int ep = 0; ep < ((MNStruct *)globalcontext)->numProfileEpochs; ep++){
	
		int NChanInEpoch = ((MNStruct *)globalcontext)->numChanPerInt[ep];

		for(int ch = 0; ch < NChanInEpoch; ch++){


			nTOA = t;

			long double FoldingPeriod = ((MNStruct *)globalcontext)->ProfileInfo[nTOA][0];
			long double FoldingPeriodDays = FoldingPeriod/SECDAY;
			int Nbins = (int)((MNStruct *)globalcontext)->ProfileInfo[nTOA][1];
			double Tobs = (double)((MNStruct *)globalcontext)->ProfileInfo[nTOA][2];
			double noiseval = (double)((MNStruct *)globalcontext)->ProfileInfo[nTOA][3];
			long double ReferencePeriod = ((MNStruct *)globalcontext)->ReferencePeriod;

			double *Betatimes = new double[Nbins];
	

		

			double noisemean=0;
			double MLSigma = 0;



			long double binpos = ModelBats[nTOA];

			if(binpos < ProfileBats[nTOA][0])binpos+=FoldingPeriodDays;

			long double minpos = binpos - FoldingPeriodDays/2;
			if(minpos < ProfileBats[nTOA][0])minpos=ProfileBats[nTOA][0];
			long double maxpos = binpos + FoldingPeriodDays/2;
			if(maxpos> ProfileBats[nTOA][Nbins-1])maxpos =ProfileBats[nTOA][Nbins-1];






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Work Out bin position for NBins/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			int InterpBin = 0;
			double FirstInterpTimeBin = 0;
			int  NumWholeBinInterpOffset = 0;
		
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

			InterpBin = floor(PWrappedTimeDiff/((MNStruct *)globalcontext)->InterpolatedTime+0.5);
			if(InterpBin >= ((MNStruct *)globalcontext)->NumToInterpolate)InterpBin -= ((MNStruct *)globalcontext)->NumToInterpolate;

			FirstInterpTimeBin = -1*(InterpBin-1)*((MNStruct *)globalcontext)->InterpolatedTime;

			if(debug == 1)printf("Interp Time Diffs: %g %g %g %g \n", ((MNStruct *)globalcontext)->InterpolatedTime, InterpBin*((MNStruct *)globalcontext)->InterpolatedTime, PWrappedTimeDiff, InterpBin*((MNStruct *)globalcontext)->InterpolatedTime-PWrappedTimeDiff);

			double FirstBinOffset = timediff-FirstInterpTimeBin;
			double dNumWholeBinOffset = FirstBinOffset/(FoldingPeriod/Nbins);
			int  NumWholeBinOffset = 0;

			NumWholeBinInterpOffset = floor(dNumWholeBinOffset+0.5);


			if(debug == 1)printf("Interp bin is: %i , First Bin is %g, Offset is %i \n", InterpBin, FirstInterpTimeBin, NumWholeBinInterpOffset);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Work Out bin position for InterpBins////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			InterpBin = 0;
			FirstInterpTimeBin = 0;
			int  InterpFacNumWholeBinInterpOffset = 0;
		
			timediff = 0;
		        bintime = ProfileBats[t][0];


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

			OneBin = FoldingPeriod/interpNBins;
			NumBinsInTimeDiff = floor(timediff/OneBin + 0.5);
			WholeBinsInTimeDiff = NumBinsInTimeDiff*FoldingPeriod/interpNBins;
			OneBinTimeDiff = -1*((double)timediff - WholeBinsInTimeDiff);

			PWrappedTimeDiff = (OneBinTimeDiff - floor(OneBinTimeDiff/OneBin)*OneBin);

			if(debug == 1)printf("Making InterpBin: %g %g %i %g %g %g\n", (double)timediff, OneBin, NumBinsInTimeDiff, WholeBinsInTimeDiff, OneBinTimeDiff, PWrappedTimeDiff);

			InterpBin = floor(PWrappedTimeDiff/((MNStruct *)globalcontext)->InterpolatedTime+0.5);
			if(InterpBin >= ((MNStruct *)globalcontext)->NumToInterpolate)InterpBin -= ((MNStruct *)globalcontext)->NumToInterpolate;

			FirstInterpTimeBin = -1*(InterpBin-1)*((MNStruct *)globalcontext)->InterpolatedTime;

			if(debug == 1)printf("Interp Time Diffs: %g %g %g %g \n", ((MNStruct *)globalcontext)->InterpolatedTime, InterpBin*((MNStruct *)globalcontext)->InterpolatedTime, PWrappedTimeDiff, InterpBin*((MNStruct *)globalcontext)->InterpolatedTime-PWrappedTimeDiff);

			FirstBinOffset = timediff-FirstInterpTimeBin;
			dNumWholeBinOffset = FirstBinOffset/(FoldingPeriod/interpNBins);
			NumWholeBinOffset = 0;

			InterpFacNumWholeBinInterpOffset = floor(dNumWholeBinOffset+0.5);


			if(debug == 1)printf("Interp bin is: %i , First Bin is %g, Offset is %i \n", InterpBin, FirstInterpTimeBin, NumWholeBinInterpOffset);		   


			int interpIndex=int(abs(InterpFacNumWholeBinInterpOffset-interpFactor*NumWholeBinInterpOffset));
			//printf("InterpBins: %i %i %i\n", InterpFacNumWholeBinInterpOffset, NumWholeBinInterpOffset, ((InterpFacNumWholeBinInterpOffset-interpFactor*NumWholeBinInterpOffset)));

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////Get Noise Level//////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//printf("%i %i \n", nTOA, Nbins);
			noisemean=0;
			int MLnoisemean=0;
			int numoffpulse = 0;
			MLSigma = 0;
			int MLSigmaCount = 0;

			double MinSigma = pow(10.0, 100);
			for(int b = 0; b < Nbins-100; b+=100){
				MLSigmaCount = 0;
				for(int j =b; j < b+100; j++){

						noisemean += (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1];
						numoffpulse += 1;
				}

				noisemean = noisemean/(numoffpulse);

				for(int j =b; j < b+100; j++){

					double res = (double)((MNStruct *)globalcontext)->ProfileData[nTOA][j][1] - noisemean;
					MLSigma += res*res; MLSigmaCount += 1;
				}

				MLSigma = sqrt(MLSigma/MLSigmaCount);

				if(MLSigma < MinSigma){
					MinSigma = MLSigma;
					MLnoisemean = noisemean;
				}
			}


			if(debug == 1)printf("noise mean %g num off pulse %i\n", MLnoisemean, numoffpulse);


			double cfreq = floor(((MNStruct *)globalcontext)->pulse->obsn[t].freq);
			cfreq = floor((cfreq-minfreq)/chanwidth)*chanwidth + minfreq;
			int freqpos = 0;

			for(int f = 0; f < numfreqs; f++){
				if(fabs(cfreq -freqs[f]) < 0.01){freqpos=f;}
			}
				
			for(int i =0; i < Nbins; i++){
				int Nj =  UtWrap(i+Nbins/2 + NumWholeBinInterpOffset, 0, Nbins-1);
				int INj =  UtWrap(interpFactor*i+interpNBins/2 + InterpFacNumWholeBinInterpOffset, 0, interpNBins-1);

				//if(t==0){printf("Interped Bins: %i %i \n", Nj, INj);}
				TScrunched[freqpos*GlobalNBins + Nj] += (double)((MNStruct *)globalcontext)->ProfileData[nTOA][i][1]*(1.0/(MinSigma*MinSigma));

				interpTScrunched[interpIndex][freqpos*GlobalNBins + Nj] += (double)(((MNStruct *)globalcontext)->ProfileData[nTOA][i][1]-MLnoisemean)*(1.0/(MinSigma*MinSigma));

				interpWeights[interpIndex][freqpos*GlobalNBins + Nj] += (1.0/(MinSigma*MinSigma));

			}

			fweights[freqpos] += (1.0/(MinSigma*MinSigma));
			t++;
		
		}		

	}

	for(int f = 0; f < numfreqs; f++){
		 for(int i =0; i < GlobalNBins; i++){

                        TScrunched[f*GlobalNBins + i] /= fweights[f];
                }

		double minval=pow(10.0, 100);
		for(int i =0; i < GlobalNBins; i++){
			if(TScrunched[f*GlobalNBins + i] < minval){minval=TScrunched[f*GlobalNBins + i];}
		}
		for(int i =0; i < GlobalNBins; i++){
			((MNStruct *)globalcontext)->TemplateChans[f*GlobalNBins + i] = TScrunched[f*GlobalNBins + i]-minval;
                        //if(f==5){printf("Tscrunched: %i %i %.10g\n", i, f, TScrunched[f*GlobalNBins + i]-minval);}
                }

		for(int ifac = 0; ifac < interpFactor; ifac++){

	 	 	for(int i =0; i < GlobalNBins; i++){
				//printf("index %i %i %i\n", i, f*interpNBins+i, numfreqs*interpNBins);
		                interpTScrunched[ifac][f*GlobalNBins + i] /= interpWeights[ifac][f*GlobalNBins + i];
		        }

			double noisemean=0;
			double MLNoiseMean=0;
			int numoffpulse = 0;
			double MLSigma = 0;
			double MLSigmaCount = 0;

			double MinSigma = pow(10.0, 100);
			for(int b = 0; b < GlobalNBins-100; b+=100){
				MLSigmaCount = 0;
				for(int j =b; j < b+100; j++){

						noisemean += interpTScrunched[ifac][f*GlobalNBins + j];
						numoffpulse += 1;
				}

				noisemean = noisemean/(numoffpulse);

				for(int j =b; j < b+100; j++){

					double res = interpTScrunched[ifac][f*GlobalNBins + j] - noisemean;
					MLSigma += res*res; MLSigmaCount += 1;
				}

				MLSigma = sqrt(MLSigma/MLSigmaCount);

				if(MLSigma < MinSigma){
					MinSigma = MLSigma;
					MLNoiseMean = noisemean;
				}
			}
			
			double dsum=0;
			for(int i =0; i < GlobalNBins; i++){
				interpTScrunched[ifac][f*GlobalNBins + i] -= MLNoiseMean;
				if(interpTScrunched[ifac][f*GlobalNBins + i] > 3*MinSigma){
					dsum+=interpTScrunched[ifac][f*GlobalNBins + i];
				}
		        }
			printf("Noise mean: %i %i %g %g \n", ifac, f, MLNoiseMean, dsum);
			for(int i =0; i < GlobalNBins; i++){
				//((MNStruct *)globalcontext)->TemplateChans[f*GlobalNBins + i] = TScrunched[f*GlobalNBins + i]-minval;
		                if(f==5){printf("ITscrunched: %g %.10g %g\n", (double(i))+double(ifac)/interpFactor, interpTScrunched[ifac][f*GlobalNBins + i]/dsum, MinSigma/dsum);}
		        }
		}

	}

	printf("made TScrunched Data \n");

	delete[] TScrunched;	
	delete[] interpTScrunched;
	delete[] interpWeights;


	printf("done \n");
}



