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
#include "dgesvd.h"
#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <limits>
#include <gsl/gsl_integration.h>


double MarginEQUAD(double d, double cg, double cg1, double *a);
double evalNGPDF(double j, double *a);

void othpl(int n,double x,double *pl){


        double a=2.0;
        double b=0.0;
        double c=1.0;
        double y0=1.0;
        double y1=2.0*x;
        pl[0]=1.0;
        pl[1]=2.0*x;


//	printf("I AM IN OTHPL %i \n", n);
        for(int k=2;k<n;k++){

                double c=2.0*(k-1.0);
//		printf("%i %g\n", k, sqrt(double(k*1.0)));
		y0=y0/sqrt(double(k*1.0));
		y1=y1/sqrt(double(k*1.0));
                double yn=(a*x+b)*y1-c*y0;
		yn=yn;///sqrt(double(k));
                pl[k]=yn;///sqrt(double(k));
                y0=y1;
                y1=yn;

        }

}

void processPDF(void *context, std::string longname){

	printf("in process pdf \n");
	std::ifstream resultsfile;
	std::string resultsfilename = longname+"post_equal_weights.dat";
	resultsfile.open(resultsfilename.c_str());

	FILE *PDFFile;
        std::string PDFFilename = longname+"NGPDFs.txt";

        PDFFile = fopen(PDFFilename.c_str(), "a+");



	double *NGCoeffs=new double[5];


	std::string line;
		
	while(getline(resultsfile,line)){
		std::istringstream myStream( line );
		std::istream_iterator< double > begin(myStream),eof;
		std::vector<double> paramlist(begin,eof);

	
		for(int j=0;j < 5;j++){
			NGCoeffs[j]=paramlist[((MNStruct *)context)->systemcount*2 + j];
		}

		for(int j=0;j < 400;j++){
			double d=-20+((double)j)*0.1;
			double P= evalNGPDF(d, NGCoeffs);
			fprintf(PDFFile, "%g", P);
			fprintf(PDFFile, " ");
		}
		fprintf(PDFFile, "\n");
	
	}
	fclose(PDFFile);

}



void subtractMLsolution(void *context){

	printf("getting residuals");

	formBatsAll(((MNStruct *)context)->pulse,1);       /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,1,1);       /* Form residuals */



        std::ofstream Resfile;
        std::string Resfilename = "MLGaussRes.dat";
        Resfile.open(Resfilename.c_str());


	double *Resvec = new double[((MNStruct *)context)->pulse->nobs];

        for(int i=0;i<((MNStruct *)context)->pulse->nobs;i++){

		//printf("%i %g %g %g %g\n", i, (double)((MNStruct *)context)->pulse->obsn[i].bat, Resvec[i], dsum,(double)((MNStruct *)context)->pulse->obsn[i].toaErr);	
		Resvec[i] = (double)((MNStruct *)context)->pulse->obsn[i].residual;

 		Resfile << std::fixed  << std::setprecision(8)  << (double)((MNStruct *)context)->pulse->obsn[i].freq << " " << (double)((MNStruct *)context)->pulse->obsn[i].sat << " ";
                Resfile << std::scientific << Resvec[i] <<  " " << (double)((MNStruct *)context)->pulse->obsn[i].toaErr*pow(10.0,-6) << "\n";

	}
        Resfile.close();
	((MNStruct *)context)->maxLikeRes=Resvec;

	delete[] Resvec;


}

void TempoNestNGConvolvedLikeFunc(double *Cube, int &ndim, int &npars, double &lnew, void *context){


	int numNGCoeff = 5;


	double **DPriors = new double *[ndim];
	for(int i =0; i < ndim;i++){
		DPriors[i]=new double[2];
	}
	
	for(int i =0; i < ((MNStruct *)context)->systemcount; i++){
		DPriors[i][0] = -3;
		DPriors[i][1] = 1;
	}

	
	for(int i =0; i < ((MNStruct *)context)->systemcount; i++){
		DPriors[((MNStruct *)context)->systemcount+i][0] = -10;
		DPriors[((MNStruct *)context)->systemcount+i][1] = -5;

	}

	for(int i =0; i < numNGCoeff; i++){
 		DPriors[2*((MNStruct *)context)->systemcount+i][0]=-1;
 		DPriors[2*((MNStruct *)context)->systemcount+i][1]=1;
	}
	
	DPriors[2*((MNStruct *)context)->systemcount+numNGCoeff][0]=-5*pow(10.0,-7);
	DPriors[2*((MNStruct *)context)->systemcount+numNGCoeff][1]=5*pow(10.0,-7);


	for(int i =0; i < ndim;i++){
		Cube[i]=(DPriors[i][1]-DPriors[i][0])*Cube[i]+DPriors[i][0];
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get White Noise vectors/////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  

	double **EFAC;
	double *EQUAD;
	int pcount=0;

	double uniformpriorterm=0;

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
					
					EFAC[n-1][o]=pow(10.0,Cube[pcount]);
					if(((MNStruct *)context)->EFACPriorType ==1) {uniformpriorterm += log(EFAC[n-1][o]);}
				}
				pcount++;
			}
			else{
                                for(int o=0;o<((MNStruct *)context)->systemcount; o++){

                                        EFAC[n-1][o]=pow(10.0,Cube[pcount]);
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
					EFAC[n-1][p]=pow(10.0,Cube[pcount]);
					if(((MNStruct *)context)->EFACPriorType ==1) {uniformpriorterm += log(EFAC[n-1][p]);}
					pcount++;
				}
			}
			else{
                                for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
                                        EFAC[n-1][p]=pow(10.0,Cube[pcount]);
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
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
			if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount]));}
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
        EQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            EQUAD[o]=pow(10.0,2*Cube[pcount]);
	    if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount])); }

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
			SQUAD[o]=pow(10.0,2*Cube[pcount]);
			if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount]));}
		}
		pcount++;
	}
	else if(((MNStruct *)context)->incShannonJitter > 1){
        SQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            SQUAD[o]=pow(10.0,2*Cube[pcount]);
	    if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount])); }

			pcount++;
        }
    }



	double *NGCoeff=new double[numNGCoeff];

	double checksum=0;
	for(int i =0; i < numNGCoeff; i++){
		NGCoeff[i]=Cube[2*((MNStruct *)context)->systemcount+i];
		checksum+=NGCoeff[i]*NGCoeff[i];

	}

	int badlike=0;
	if(checksum>1){
		badlike=1;
	}

	double offset=Cube[2*((MNStruct *)context)->systemcount+numNGCoeff];

	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);      /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,1,1);       /* Form residuals */


	
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]= (double)((MNStruct *)context)->pulse->obsn[o].residual-offset;
		//printf("Res: %i %g %g \n", o, Resvec[o],((MNStruct *)context)->pulse->obsn[o].origErr*pow(10.0,-6));

	}





	
	double Chisq=0;
	double detN=0;
	double NGLike=0;

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

		double EF = EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]];
		double EQ = EQUAD[((MNStruct *)context)->sysFlags[o]];
		double ShannonJitterTerm=0;

		double toaerr=0;
		if(((MNStruct *)context)->useOriginalErrors==0){
			toaerr=((MNStruct *)context)->pulse->obsn[o].toaErr*pow(10.0,-6);
		}
		else if(((MNStruct *)context)->useOriginalErrors==1){
			toaerr=((MNStruct *)context)->pulse->obsn[o].origErr*pow(10.0,-6);
		}

		double noiseval = EF*toaerr*toaerr + EQ;
		double chisqterm = Resvec[o]*Resvec[o]/noiseval;

		Chisq += 0.5*chisqterm;
		detN += (25.0/2.0)*log(noiseval);

		double NGTerm = MarginEQUAD(Resvec[o], sqrt(EF*toaerr*toaerr), sqrt(EQ), NGCoeff);
		NGLike+= log(NGTerm);


	}


	if(isnan(detN) || isinf(detN) || isnan(Chisq) || isinf(Chisq) || badlike==1){

		lnew=-pow(10.0,20);
	}
	else{
		lnew = -Chisq - detN + NGLike+uniformpriorterm;
	}




	delete[] Resvec;
	delete[] NGCoeff;

	for(int i =0; i < ndim;i++){
		delete[] DPriors[i];
	}
	for(int n=0; n <((MNStruct *)context)->EPolTerms; n++){
		delete[] EFAC[n];
	}
	delete[] EFAC;
	delete[] EQUAD;
	delete[] SQUAD;
	delete[] DPriors;
}


void TempoNestNGLikeFunc(double *Cube, int &ndim, int &npars, double &lnew, void *context){


	int numNGCoeff = 6;


	double **DPriors = new double *[ndim];
	for(int i =0; i < ndim;i++){
		DPriors[i]=new double[2];
	}
	
	for(int i =0; i < ((MNStruct *)context)->systemcount; i++){
		DPriors[i][0] = -3;
		DPriors[i][1] = 1;
	}

	
	for(int i =0; i < ((MNStruct *)context)->systemcount; i++){
		DPriors[((MNStruct *)context)->systemcount+i][0] = -10;
		DPriors[((MNStruct *)context)->systemcount+i][1] = -5;

	}

	for(int i =0; i < numNGCoeff-1; i++){
 		DPriors[2*((MNStruct *)context)->systemcount+i][0]=-1;
 		DPriors[2*((MNStruct *)context)->systemcount+i][1]=1;
	}
	
	DPriors[2*((MNStruct *)context)->systemcount+numNGCoeff-1][0]=-5*pow(10.0,-7);
	DPriors[2*((MNStruct *)context)->systemcount+numNGCoeff-1][1]=5*pow(10.0,-7);


	for(int i =0; i < ndim;i++){
		Cube[i]=(DPriors[i][1]-DPriors[i][0])*Cube[i]+DPriors[i][0];
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get White Noise vectors/////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  

	double **EFAC;
	double *EQUAD;
	int pcount=0;

	double uniformpriorterm=0;

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
					
					EFAC[n-1][o]=pow(10.0,Cube[pcount]);
					if(((MNStruct *)context)->EFACPriorType ==1) {uniformpriorterm += log(EFAC[n-1][o]);}
				}
				pcount++;
			}
			else{
                                for(int o=0;o<((MNStruct *)context)->systemcount; o++){

                                        EFAC[n-1][o]=pow(10.0,Cube[pcount]);
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
					EFAC[n-1][p]=pow(10.0,Cube[pcount]);
					if(((MNStruct *)context)->EFACPriorType ==1) {uniformpriorterm += log(EFAC[n-1][p]);}
					pcount++;
				}
			}
			else{
                                for(int p=0;p< ((MNStruct *)context)->systemcount; p++){
                                        EFAC[n-1][p]=pow(10.0,Cube[pcount]);
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
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
			if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount]));}
		}
		pcount++;
	}
	else if(((MNStruct *)context)->numFitEQUAD > 1){
        EQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            EQUAD[o]=pow(10.0,2*Cube[pcount]);
	    if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount])); }

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
			SQUAD[o]=pow(10.0,2*Cube[pcount]);
			if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount]));}
		}
		pcount++;
	}
	else if(((MNStruct *)context)->incShannonJitter > 1){
        SQUAD=new double[((MNStruct *)context)->systemcount];
        for(int o=0;o<((MNStruct *)context)->systemcount; o++){
            SQUAD[o]=pow(10.0,2*Cube[pcount]);
	    if(((MNStruct *)context)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount])); }

			pcount++;
        }
    }



	double *NGCoeff=new double[numNGCoeff];

	double checksum=0;
//	NGCoeff[1] = 0;
	for(int i =1; i < numNGCoeff; i++){
		NGCoeff[i]=Cube[2*((MNStruct *)context)->systemcount+i-1];
		checksum+=NGCoeff[i]*NGCoeff[i];

	}
	


	int badlike=0;
	if(checksum>1){
		badlike=1;
	}
	else{	
		NGCoeff[0]=sqrt(1.0-checksum);
	}
	double offset=Cube[2*((MNStruct *)context)->systemcount+numNGCoeff-1];

	fastformBatsAll(((MNStruct *)context)->pulse,((MNStruct *)context)->numberpulsars);      /* Form Barycentric arrival times */
	formResiduals(((MNStruct *)context)->pulse,1,1);       /* Form residuals */


	
	double *Resvec=new double[((MNStruct *)context)->pulse->nobs];

	double meanphase=0;
	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		meanphase += (double)((MNStruct *)context)->pulse->obsn[o].residual;
		//Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual+phase;

	}
	meanphase = meanphase/((MNStruct *)context)->pulse->nobs;

	for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)context)->pulse->obsn[o].residual-meanphase;

	}


	
	double Chisq=0;
	double detN=0;
	double NGLike=0;
	double *HVec=new double[numNGCoeff];



        for(int o=0;o<((MNStruct *)context)->pulse->nobs; o++){

		double EF = EFAC[0][((MNStruct *)context)->sysFlags[o]]*EFAC[0][((MNStruct *)context)->sysFlags[o]];
		double EQ = EQUAD[((MNStruct *)context)->sysFlags[o]];
		double ShannonJitterTerm=0;

		double toaerr=0;
		if(((MNStruct *)context)->useOriginalErrors==0){
			toaerr=((MNStruct *)context)->pulse->obsn[o].toaErr*pow(10.0,-6);
		}
		else if(((MNStruct *)context)->useOriginalErrors==1){
			toaerr=((MNStruct *)context)->pulse->obsn[o].origErr*pow(10.0,-6);
		}

		double noiseval = sqrt(EF*toaerr*toaerr + EQ);


                double dsum=0;

                othpl(numNGCoeff,Resvec[o]/(sqrt(2.0)*noiseval),HVec);
                Chisq += -0.5*Resvec[o]*Resvec[o]/(noiseval*noiseval);
                for(int k =0; k < numNGCoeff; k++){

                        double Cconst=1.0/sqrt(sqrt(noiseval*noiseval)*sqrt(2.0*M_PI)*pow(2.0,k)*iter_factorial(k));
                        dsum += HVec[k]*Cconst*NGCoeff[k];
                }

                NGLike += 2*log(fabs(dsum));


        }

        



	if(isnan(detN) || isinf(detN) || isnan(Chisq) || isinf(Chisq) || badlike==1){

		lnew=-pow(10.0,20);
	}
	else{
		lnew = Chisq + NGLike + uniformpriorterm;
	}




	delete[] Resvec;
	delete[] NGCoeff;

	for(int i =0; i < ndim;i++){
		delete[] DPriors[i];
	}
	for(int n=0; n <((MNStruct *)context)->EPolTerms; n++){
		delete[] EFAC[n];
	}
	delete[] EFAC;
	delete[] EQUAD;
	delete[] SQUAD;
	delete[] DPriors;
	delete[] HVec;
}

double MarginEQUAD(double d, double cg, double cg1, double *a){

//d = residual (data - Me - Fa)
//cg = efac*toaerr
//cg1=equad
//a = non gaussian coefficients

double t1 = a[1] * a[2];
double  t2 = pow(cg1, 0.2e1);
double  t3 = t2 * cg1;
double  t4 = t2 * t2;
double  t5 = t4 * t3;
double  t6 = t4 * t4;
double  t7 = t6 * t6;
double  t9 = t7 * t5 * d;
double  t12 = pow(a[1], 0.2e1);
double  t13 = pow(cg, 0.2e1);
double  t14 = t13 * t13;
double  t15 = t14 * t13;
double  t16 = t14 * t14;
double  t17 = t16 * t16;
double  t18 = t17 * t15;
double  t22 = pow(a[4], 0.2e1);
double  t24 = t7 * t4;
double  t27 = pow(a[0], 0.2e1);
double  t28 = t17 * t13;
double  t30 = t4 * t2;
double  t33 = d * d;
double  t34 = t33 * t33;
double  t35 = t34 * t22;
double  t38 = pow(a[3], 0.2e1);
double  t45 = pow(a[2], 0.2e1);
double  t46 = t17 * t14;
double  t50 = t16 * t15;
double  t52 = t6 * t2;
double  t55 = t16 * t13;
double  t57 = t6 * t30;
double  t60 = t34 * t38;
double  t63 = -0.1080e4 * t9 * t1 - 0.2160e4 * t2 * t18 * t12 - 0.40095e5 * t24 * t14 * t22 - 0.68400e5 * t30 * t28 * t27 + 0.2475e4 * t24 * t35 - 0.2839500000e6 * t7 * t16 * t38 - 0.1116000000e6 * t30 * t28 * t22 - 0.25200e5 * t4 * t46 * t45 - 0.3456000000e6 * t52 * t50 * t27 - 0.4466700000e6 * t57 * t55 * t38 - 0.1800e4 * t24 * t60;
double  t64 = t16 * t14;
double  t66 = t6 * t4;
double  t69 = t33 * t12;
double  t70 = t7 * t30;
double  t82 = t33 * t22;
double  t94 = a[2] * a[0];
double  t95 = t7 * t6;
double  t98 = t34 * t33;
double  t99 = t98 * t38;
double  t100 = t7 * t2;
double  t106 = -0.4863600000e6 * t66 * t64 * t45 + 0.1080e4 * t70 * t69 - 0.30600e5 * t24 * t14 * t27 - 0.20520e5 * t4 * t46 * t12 - 0.1008000000e6 * t30 * t28 * t45 - 0.1350e4 * t70 * t82 - 0.1890000000e6 * t6 * t17 * t27 - 0.6480e4 * t70 * t13 * t45 - 0.4554900000e6 * t57 * t55 * t22 - 0.6235382907e3 * t95 * t94 + 0.780e3 * t100 * t99 - 0.2452500000e6 * t6 * t17 * t45;
double  t111 = a[0] * a[4];
double  t132 = t33 * t38;
double  t141 = -0.4284000000e6 * t57 * t55 * t45 + 0.5692099788e3 * t95 * t111 - 0.4068000000e6 * t52 * t50 * t45 - 0.4336200000e6 * t52 * t50 * t22 - 0.7290e4 * t70 * t13 * t38 - 0.2484000000e6 * t7 * t16 * t27 - 0.3931200000e6 * t57 * t55 * t27 - 0.5058900000e6 * t66 * t64 * t38 + 0.1350e4 * t70 * t132 - 0.2268000000e6 * t6 * t17 * t12 - 0.4384800000e6 * t66 * t64 * t27;
double  t148 = t34 * t34;
double  t167 = sqrt(0.1e1 - 0.1e1 * t27 - 0.1e1 * t12 - 0.1e1 * t45 - 0.1e1 * t38 - 0.1e1 * t22);
double  t168 = a[2] * t167;
double  t177 = t34 * t12;
double  t186 = -0.1264500000e6 * t100 * t15 * t38 - 0.37800e5 * t24 * t14 * t12 + t66 * t148 * t34 * t22 - 0.4082400000e6 * t57 * t55 * t12 - 0.5225850000e6 * t66 * t64 * t22 - 0.4687200000e6 * t66 * t64 * t12 + 0.8818163074e3 * t95 * t168 - 0.3888000000e6 * t52 * t50 * t12 - 0.28800e5 * t4 * t46 * t38 - 0.720e3 * t24 * t177 - 0.88200e5 * t30 * t28 * t12 - 0.2605500000e6 * t6 * t17 * t22;
double  t210 = t148 * t22;
double  t213 = t34 * t45;
double  t216 = t148 * t45;
double  t222 = -0.1296000000e6 * t100 * t15 * t22 - 0.1440e4 * t2 * t18 * t27 - 0.1260000000e6 * t100 * t15 * t45 - 0.1188000000e6 * t100 * t15 * t12 - 0.4196700000e6 * t52 * t50 * t38 - 0.2868750000e6 * t7 * t16 * t22 - 0.38070e5 * t24 * t14 * t38 + 0.315e3 * t7 * t210 + 0.1260e4 * t24 * t213 + 0.30e2 * t7 * t216 - 0.2763000000e6 * t7 * t16 * t45;
double  t226 = a[2] * a[4];
double  t238 = t98 * t45;
double  t241 = a[0] * t167;
double  t247 = t98 * t12;
double  t256 = t98 * t22;
double  t259 = -0.37800e5 * t24 * t14 * t45 - 0.4929503018e3 * t95 * t226 - 0.14760e5 * t4 * t46 * t27 - 0.7560e4 * t70 * t13 * t12 - 0.7290e4 * t70 * t13 * t22 - 0.360e3 * t100 * t238 - 0.1018233765e4 * t95 * t241 - 0.31320e5 * t4 * t46 * t22 + 0.120e3 * t100 * t247 - 0.5040e4 * t70 * t13 * t27 - 0.2542500000e6 * t6 * t17 * t38 - 0.1380e4 * t100 * t256;
double  t273 = t148 * t38;
double  t276 = t148 * t33;
double  t277 = t276 * t22;
double  t280 = t34 * t27;
double  t283 = t33 * t27;
double  t286 = t33 * t45;
double  t289 = t276 * t38;
double  t292 = a[4] * t167;
double  t295 = -0.2592000000e6 * t7 * t16 * t12 - 0.1080000000e6 * t100 * t15 * t27 - 0.4320e4 * t2 * t18 * t22 - 0.3600e4 * t2 * t18 * t38 - 0.120e3 * t7 * t273 - 0.30e2 * t57 * t277 + 0.360e3 * t24 * t280 - 0.720e3 * t70 * t283 - 0.1080e4 * t70 * t286 + 0.6e1 * t57 * t289 - 0.8049844719e3 * t95 * t292;
double  t304 = a[0] * a[1];
double  t305 = t33 * d;
double  t306 = t4 * cg1;
double  t308 = t17 * t306 * t305;
double  t312 = t50 * t5 * t305;
double  t315 = t6 * cg1;
double  t317 = t64 * t315 * t305;
double  t320 = t6 * t3;
double  t322 = t55 * t320 * t305;
double  t325 = t6 * t306;
double  t327 = t16 * t325 * t305;
double  t330 = t6 * t5;
double  t332 = t15 * t330 * t305;
double  t335 = t7 * cg1;
double  t337 = t14 * t335 * t305;
double  t340 = t7 * t3;
double  t342 = t13 * t340 * t305;
double  t346 = t18 * cg1 * d;
double  t350 = t46 * t3 * d;
double  t353 = 0.1095445115e2 * t14 * t52 * t276 * t226 + 0.2190890230e2 * t13 * t66 * t276 * t226 + 0.1829045653e5 * t308 * t304 + 0.5653413836e5 * t312 * t304 + 0.9311505141e5 * t317 * t304 + 0.8147566999e5 * t322 * t304 + 0.2327876285e5 * t327 * t304 - 0.2327876285e5 * t332 * t304 - 0.2660430040e5 * t337 * t304 - 0.1080799704e5 * t342 * t304 + 0.2494153163e4 * t346 * t304 + 0.1995322530e5 * t350 * t304;
double  t358 = t28 * t306 * d;
double  t362 = t17 * t5 * d;
double  t366 = t50 * t315 * d;
double  t370 = t64 * t320 * d;
double  t373 = a[2] * a[3];
double  t374 = t148 * d;
double  t376 = t15 * t315 * t374;
double  t380 = t14 * t320 * t374;
double  t384 = t13 * t325 * t374;
double  t387 = t34 * t305;
double  t389 = t55 * t5 * t387;
double  t393 = t16 * t315 * t387;
double  t397 = t15 * t320 * t387;
double  t401 = t14 * t325 * t387;
double  t404 = 0.6858921198e5 * t358 * t304 + 0.1309430411e6 * t362 * t304 + 0.1496491898e6 * t366 * t304 + 0.1047544328e6 * t370 * t304 + 0.2683281573e2 * t376 * t373 + 0.8049844719e2 * t380 * t373 + 0.8049844719e2 * t384 * t373 + 0.5366563146e3 * t389 * t373 + 0.1717300207e4 * t393 * t373 + 0.1502637681e4 * t397 * t373 - 0.4293250517e3 * t401 * t373;
double  t406 = t13 * t330 * t387;
double  t409 = a[0] * a[3];
double  t420 = a[3] * t167;
double  t422 = t14 * t340 * d;
double  t425 = t34 * d;
double  t427 = t14 * t330 * t425;
double  t431 = t13 * t335 * t425;
double  t435 = t64 * t5 * t425;
double  t439 = t55 * t315 * t425;
double  t442 = a[1] * t167;
double  t444 = t15 * t335 * d;
double  t447 = -0.1180643892e4 * t406 * t373 + 0.9295160031e2 * t389 * t409 + 0.4647580015e3 * t393 * t409 + 0.9295160031e3 * t397 * t409 + 0.9295160031e3 * t401 * t409 + 0.4647580015e3 * t406 * t409 + 0.7098484345e5 * t422 * t420 + 0.2760521690e4 * t427 * t420 + 0.9201738966e3 * t431 * t420 + 0.9201738966e3 * t435 * t420 + 0.2760521690e4 * t439 * t420 - 0.2116359138e6 * t444 * t442;
double  t454 = t17 * t4 * t34;
double  t458 = t55 * t34 * t52;
double  t462 = t15 * t325 * t425;
double  t466 = t50 * t306 * t425;
double  t479 = -0.7936346767e5 * t362 * t442 - 0.2116359138e6 * t366 * t442 + 0.2939387691e3 * t454 * t168 + 0.1646057107e5 * t458 * t168 + 0.4600869483e4 * t462 * t420 + 0.1314534138e3 * t466 * t420 - 0.1314534138e4 * t308 * t420 - 0.1051627310e5 * t312 * t420 - 0.3680695586e5 * t317 * t420 - 0.7361391173e5 * t322 * t420 + 0.1971801207e4 * t358 * t420;
double  t481 = t16 * t320 * t425;
double  t493 = t16 * t330 * d;
double  t503 = t14 * t7 * t34;
double  t507 = t13 * t100 * t34;
double  t511 = t64 * t6 * t34;
double  t514 = 0.4600869483e4 * t481 * t420 - 0.9201738966e5 * t327 * t420 - 0.7361391173e5 * t332 * t420 - 0.3680695586e5 * t337 * t420 - 0.1051627310e5 * t342 * t420 + 0.2484469521e6 * t493 * t420 + 0.1774621086e5 * t362 * t420 + 0.7098484345e5 * t366 * t420 + 0.4938171321e5 * t332 * t442 + 0.8230285536e4 * t503 * t168 + 0.2351510153e4 * t507 * t168 + 0.8230285536e4 * t511 * t168;
double  t526 = t55 * t325 * d;
double  t530 = t28 * t3 * t305;
double  t536 = t28 * t4 * t33;
double  t545 = 0.2116359138e5 * t337 * t442 + 0.5290897844e4 * t342 * t442 + 0.5290897844e4 * t308 * t442 + 0.1656313014e6 * t370 * t420 + 0.2484469521e6 * t526 * t420 + 0.5878775383e3 * t530 * t442 - 0.1763632615e4 * t350 * t442 - 0.1763632615e4 * t536 * t168 + 0.2880e4 * t346 * t1 + 0.18720e5 * t350 * t1 + 0.50400e5 * t358 * t1;
double  t558 = t7 * t306;
double  t560 = t13 * t558 * d;
double  t573 = 0.73800e5 * t362 * t1 + 0.72000e5 * t366 * t1 + 0.70560e5 * t370 * t1 + 0.80640e5 * t526 * t1 + 0.68400e5 * t493 * t1 + 0.28800e5 * t444 * t1 - 0.4320e4 * t560 * t1 + 0.1440e4 * t466 * t1 + 0.7560e4 * t435 * t1 + 0.15120e5 * t439 * t1 + 0.12600e5 * t481 * t1 - 0.7560e4 * t427 * t1;
double  t599 = -0.5040e4 * t431 * t1 + 0.4320e4 * t530 * t1 + 0.24480e5 * t308 * t1 + 0.52920e5 * t312 * t1 + 0.47880e5 * t317 * t1 + 0.2520e4 * t322 * t1 - 0.22680e5 * t327 * t1 - 0.2520e4 * t332 * t1 + 0.16920e5 * t337 * t1 + 0.11880e5 * t342 * t1 + 0.120e3 * t389 * t1 + 0.600e3 * t393 * t1;
double  t607 = t64 * t30 * t98;
double  t614 = a[1] * a[3];
double  t618 = t55 * t6 * t98;
double  t622 = t16 * t52 * t98;
double  t626 = t15 * t66 * t98;
double  t630 = t14 * t57 * t98;
double  t634 = t13 * t7 * t98;
double  t637 = 0.1200e4 * t397 * t1 + 0.1200e4 * t401 * t1 + 0.600e3 * t406 * t1 + 0.4553679831e3 * t607 * t111 + 0.6439875775e4 * t337 * t373 - 0.8049844719e4 * t342 * t373 + 0.8049844719e3 * t607 * t614 + 0.3327269151e4 * t618 * t614 + 0.4561578674e4 * t622 * t614 + 0.1073312629e4 * t626 * t614 - 0.2951609730e4 * t630 * t614 - 0.2683281573e4 * t634 * t614;
double  t645 = t50 * t30 * t34;
double  t653 = t16 * t66 * t34;
double  t657 = t15 * t57 * t34;
double  t670 = 0.1662768775e4 * t454 * t94 + 0.1018445875e5 * t645 * t94 + 0.2473368553e5 * t511 * t94 + 0.2764353089e5 * t458 * t94 + 0.7274613392e4 * t653 * t94 - 0.1600414946e5 * t657 * t94 - 0.1891399482e5 * t503 * t94 - 0.8521689973e4 * t507 * t94 + 0.3219937888e4 * t346 * t373 + 0.1609968944e5 * t350 * t373 + 0.3219937888e5 * t358 * t373;
double  t695 = 0.4024922360e5 * t362 * t373 + 0.5755638974e5 * t366 * t373 + 0.8734081520e5 * t370 * t373 + 0.8170592390e5 * t526 * t373 + 0.3421184006e5 * t493 * t373 + 0.2012461180e4 * t444 * t373 + 0.1207476708e4 * t422 * t373 + 0.3622430124e4 * t560 * t373 + 0.5634891303e4 * t431 * t373 + 0.2078460969e3 * t607 * t94 + 0.1247076581e4 * t618 * t94 + 0.3117691454e4 * t622 * t94;
double  t704 = t16 * t6 * t148;
double  t708 = t15 * t52 * t148;
double  t712 = t14 * t66 * t148;
double  t716 = t13 * t57 * t148;
double  t727 = 0.4156921938e4 * t626 * t94 + 0.3117691454e4 * t630 * t94 + 0.1247076581e4 * t634 * t94 + 0.5366563146e2 * t704 * t614 + 0.2146625258e3 * t708 * t614 + 0.3219937888e3 * t712 * t614 + 0.2146625258e3 * t716 * t614 + 0.5237721642e5 * t526 * t304 + 0.3741229744e5 * t493 * t304 + 0.3741229744e5 * t444 * t304 + 0.2494153163e5 * t422 * t304;
double  t734 = t55 * t66;
double  t737 = t16 * t57;
double  t740 = t15 * t7;
double  t743 = t14 * t100;
double  t746 = t13 * t24;
double  t749 = t64 * t30;
double  t752 = t55 * t6;
double  t755 = t16 * t52;
double  t758 = t15 * t66;
double  t761 = 0.8729536070e4 * t560 * t304 + 0.4156921938e3 * t466 * t304 + 0.2909845357e4 * t435 * t304 + 0.43470e5 * t734 * t82 + 0.9450e4 * t737 * t82 + 0.8100e4 * t740 * t82 + 0.10800e5 * t743 * t82 + 0.1350e4 * t746 * t82 + 0.480e3 * t749 * t238 + 0.2040e4 * t752 * t238 + 0.3000e4 * t755 * t238 + 0.1200e4 * t758 * t238;
double  t764 = t14 * t57;
double  t767 = t13 * t7;
double  t770 = t17 * t4;
double  t773 = t50 * t30;
double  t776 = t64 * t6;
double  t779 = t55 * t52;
double  t782 = t15 * t57;
double  t785 = t14 * t7;
double  t788 = t13 * t100;
double  t791 = t16 * t6;
double  t794 = t15 * t52;
double  t797 = -0.1200e4 * t764 * t238 - 0.1320e4 * t767 * t238 + 0.3600e4 * t770 * t60 + 0.10800e5 * t773 * t60 + 0.6300e4 * t776 * t60 - 0.6300e4 * t779 * t60 + 0.12600e5 * t782 * t60 + 0.6300e4 * t785 * t60 - 0.2700e4 * t788 * t60 + 0.30e2 * t791 * t216 + 0.120e3 * t794 * t216;
double  t798 = t14 * t66;
double  t801 = t13 * t57;
double  t820 = t16 * t66;
double  t825 = 0.180e3 * t798 * t216 + 0.120e3 * t801 * t216 + 0.2160e4 * t770 * t213 + 0.10080e5 * t773 * t213 + 0.16380e5 * t776 * t213 + 0.7560e4 * t779 * t213 + 0.360e3 * t770 * t280 + 0.2880e4 * t773 * t280 + 0.10080e5 * t776 * t280 + 0.20160e5 * t779 * t280 + 0.25200e5 * t820 * t280 + 0.20160e5 * t782 * t280;
double  t831 = t46 * t2;
double  t834 = t28 * t4;
double  t837 = t17 * t30;
double  t840 = t50 * t6;
double  t843 = t64 * t52;
double  t856 = 0.10080e5 * t785 * t280 + 0.2880e4 * t788 * t280 + 0.1440e4 * t831 * t283 + 0.12240e5 * t834 * t283 + 0.45360e5 * t837 * t283 + 0.95040e5 * t840 * t283 + 0.1209600000e6 * t843 * t283 + 0.90720e5 * t734 * t283 + 0.30240e5 * t737 * t283 - 0.8640e4 * t740 * t283 - 0.12960e5 * t743 * t283 - 0.5040e4 * t746 * t283;
double  t857 = t14 * t52;
double  t860 = t13 * t66;
double  t863 = t100 * t98;
double  t866 = t24 * t34;
double  t873 = a[3] * a[4];
double  t880 = t57 * t55;
double  t883 = t7 * t16;
double  t886 = t100 * t15;
double  t889 = 0.6e1 * t857 * t289 + 0.12e2 * t860 * t289 - 0.6976532090e3 * t863 * t614 - 0.1454922678e4 * t866 * t94 + 0.1207476708e4 * t9 * t373 + 0.1247076581e4 * t9 * t304 - 0.1102270384e4 * t9 * t873 - 0.1394274005e4 * t9 * t409 - 0.6071573108e3 * t863 * t111 - 0.2449963000e6 * t880 * t226 - 0.1454203390e6 * t883 * t226 - 0.5866108591e5 * t886 * t226;
double  t893 = t24 * t14;
double  t896 = a[1] * a[4];
double  t897 = t335 * t387;
double  t900 = t330 * t374;
double  t905 = t70 * t13;
double  t908 = t6 * t17;
double  t911 = t52 * t50;
double  t914 = t66 * t64;
double  t917 = t340 * t425;
double  t920 = t7 * t148;
double  t923 = t558 * t305;
double  t926 = -0.1626735996e5 * t893 * t226 - 0.3943602414e3 * t897 * t896 + 0.2190890230e2 * t900 * t896 + 0.9859006035e3 * t9 * t896 - 0.2414953416e4 * t905 * t614 - 0.1794339098e6 * t908 * t226 - 0.2627425108e6 * t911 * t226 - 0.2933054295e6 * t914 * t226 + 0.1971801207e4 * t917 * t896 - 0.2300434742e3 * t920 * t226 - 0.3286335345e4 * t923 * t896;
double  t929 = t70 * t33;
double  t939 = t2 * t18;
double  t942 = t4 * t46;
double  t945 = t30 * t28;
double  t956 = -0.3450652112e4 * t905 * t226 + 0.2464751509e4 * t929 * t226 + 0.1511714259e4 * t863 * t226 + 0.1095445115e2 * t57 * t276 * t226 - 0.3614968880e4 * t866 * t226 - 0.3943602414e4 * t939 * t226 - 0.2760521690e5 * t942 * t226 - 0.8873105432e5 * t945 * t226 - 0.3219937888e4 * t923 * t373 - 0.4293250517e3 * t897 * t373 + 0.2092959627e4 * t917 * t373 + 0.2078460969e3 * t863 * t94;
double  t980 = -0.2414953416e4 * t929 * t614 - 0.1080e4 * t917 * t1 + 0.2520e4 * t923 * t1 + 0.120e3 * t897 * t1 + 0.150e3 * t791 * t273 + 0.330e3 * t794 * t273 + 0.90e2 * t798 * t273 - 0.210e3 * t801 * t273 - 0.6761869564e5 * t914 * t292 - 0.5600285707e5 * t893 * t241 + 0.2414953416e4 * t929 * t292;
double  t1005 = 0.5366563146e2 * t863 * t292 - 0.6761869564e5 * t886 * t292 - 0.4704239994e6 * t880 * t241 - 0.3360171424e6 * t883 * t241 - 0.1680085712e6 * t886 * t241 - 0.1018233765e4 * t939 * t241 - 0.8049844719e3 * t866 * t292 - 0.7244860247e4 * t905 * t292 - 0.1314534138e4 * t923 * t420 + 0.1314534138e3 * t917 * t420 + 0.1058179569e6 * t911 * t168 + 0.1851814246e6 * t914 * t168;
double  t1030 = 0.2222177095e6 * t880 * t168 + 0.1851814246e6 * t883 * t168 + 0.1058179569e6 * t886 * t168 - 0.1763632615e4 * t9 * t442 + 0.8818163074e3 * t942 * t168 + 0.8818163074e4 * t945 * t168 + 0.5878775383e3 * t923 * t442 + 0.3968173383e5 * t908 * t168 - 0.1763632615e4 * t929 * t168 + 0.2939387691e3 * t866 * t168 + 0.2732207898e5 * t945 * t111;
double  t1055 = 0.9619648642e5 * t908 * t111 + 0.1958082327e6 * t911 * t111 + 0.2550060705e6 * t914 * t111 + 0.2231303117e6 * t880 * t111 + 0.1354719750e6 * t883 * t111 + 0.5919783780e5 * t886 * t111 + 0.1935313928e5 * t893 * t111 + 0.4553679831e4 * t905 * t111 + 0.2323790008e4 * t923 * t409 + 0.5366563146e2 * t920 * t614 + 0.2683281573e2 * t900 * t373 - 0.2494153163e4 * t939 * t94;
double  t1081 = -0.2244737847e5 * t942 * t94 - 0.9041305216e5 * t945 * t94 - 0.2151207103e6 * t908 * t94 - 0.3367106770e6 * t911 * t94 - 0.3666405149e6 * t914 * t94 - 0.2880746903e6 * t880 * t94 - 0.1683553385e6 * t883 * t94 - 0.7482459489e5 * t886 * t94 - 0.2494153163e5 * t893 * t94 - 0.5611844617e4 * t905 * t94 - 0.1662768775e4 * t923 * t304 + 0.4156921938e3 * t917 * t304;
double  t1106 = -0.3219937888e4 * t939 * t614 - 0.2575950310e5 * t942 * t614 - 0.9257321427e5 * t945 * t614 - 0.1996361490e6 * t908 * t614 - 0.2930143478e6 * t911 * t614 - 0.3155539130e6 * t914 * t614 - 0.2592050000e6 * t880 * t614 - 0.1609968944e6 * t883 * t614 - 0.7083863353e5 * t886 * t614 + 0.2414953416e4 * t866 * t614 - 0.1931962733e5 * t893 * t614 + 0.1870614872e4 * t929 * t94;
double  t1134 = -0.6300e4 * t820 * t213 - 0.5040e4 * t782 * t213 + 0.3780e4 * t785 * t213 + 0.2880e4 * t831 * t286 + 0.15840e5 * t834 * t286 + 0.34560e5 * t837 * t286 + 0.39240e5 * t840 * t286 + 0.32760e5 * t843 * t286 + 0.37800e5 * t734 * t286 + 0.42840e5 * t737 * t286 + 0.25560e5 * t740 * t286;
double  t1159 = 0.5400e4 * t770 * t35 + 0.7200e4 * t773 * t35 - 0.6300e4 * t776 * t35 + 0.20475e5 * t820 * t35 + 0.6300e4 * t782 * t35 - 0.9450e4 * t785 * t35 - 0.900e3 * t788 * t35 + 0.1200e4 * t749 * t99 + 0.3000e4 * t752 * t99 + 0.780e3 * t755 * t99 - 0.2880e4 * t758 * t99 - 0.1320e4 * t764 * t99;
double  t1183 = 0.1320e4 * t767 * t99 + 0.120e3 * t749 * t247 + 0.720e3 * t752 * t247 + 0.1800e4 * t755 * t247 + 0.2400e4 * t758 * t247 + 0.1800e4 * t764 * t247 + 0.720e3 * t767 * t247 + 0.1080e4 * t770 * t177 + 0.6840e4 * t773 * t177 + 0.17640e5 * t776 * t177 + 0.22680e5 * t779 * t177;
double  t1208 = 0.12600e5 * t820 * t177 - 0.2520e4 * t782 * t177 - 0.7560e4 * t785 * t177 - 0.3960e4 * t788 * t177 + 0.2160e4 * t831 * t69 + 0.15120e5 * t834 * t69 + 0.44280e5 * t837 * t69 + 0.69120e5 * t840 * t69 + 0.60480e5 * t843 * t69 + 0.30240e5 * t734 * t69 + 0.15120e5 * t737 * t69 + 0.17280e5 * t740 * t69;
double  t1233 = 0.15120e5 * t743 * t69 + 0.4680e4 * t788 * t213 + 0.3240e4 * t743 * t286 - 0.3240e4 * t746 * t286 + 0.450e3 * t791 * t210 + 0.180e3 * t794 * t210 - 0.675e3 * t798 * t210 - 0.90e2 * t801 * t210 + 0.2400e4 * t749 * t256 + 0.1800e4 * t752 * t256 - 0.4320e4 * t755 * t256;
double  t1258 = -0.1140e4 * t758 * t256 + 0.4500e4 * t764 * t256 + 0.540e3 * t767 * t256 + 0.6480e4 * t746 * t69 + 0.36e2 * t857 * t277 + 0.6e1 * t860 * t277 + 0.3600e4 * t831 * t132 + 0.14400e5 * t834 * t132 + 0.21600e5 * t837 * t132 + 0.23400e5 * t840 * t132 + 0.40950e5 * t843 * t132 + 0.56700e5 * t734 * t132;
double  t1283 = 0.34650e5 * t737 * t132 + 0.3600e4 * t740 * t132 - 0.1350e4 * t743 * t132 + 0.2700e4 * t746 * t132 + 0.4320e4 * t831 * t82 + 0.10800e5 * t834 * t82 + 0.10800e5 * t837 * t82 + 0.27000e5 * t840 * t82 + 0.56700e5 * t843 * t82 - 0.1080000000e6 * t30 * t28 * t38 + 0.9295160031e2 * t897 * t409;
double  t1294 = t148 * t305;
double  t1310 = 0.4041658076e4 * t923 * t873 + 0.1028785692e4 * t897 * t873 - 0.1224744871e3 * t900 * t873 - 0.3380295845e4 * t917 * t873 - 0.2276839915e4 * t929 * t111 + 0.4898979486e1 * t325 * t1294 * t873 + 0.2276839915e4 * t866 * t111 - 0.1022467603e4 * t917 * t409 + 0.3794733192e2 * t920 * t111 + 0.3415259873e4 * t942 * t111 - 0.2880e4 * t2 * t18 * t45 + 0.720e3 * t95;
double  t1325 = t13 * t24 * t33;
double  t1338 = 0.720e3 * t17 * t16 + 0.1028785692e5 * t435 * t873 - 0.6172714152e4 * t439 * t873 - 0.1337421400e5 * t481 * t873 + 0.8230285536e4 * t462 * t873 - 0.9107359661e4 * t1325 * t111 + 0.8818163074e4 * t530 * t873 + 0.2057571384e5 * t308 * t873 + 0.1469693846e4 * t312 * t873 - 0.1028785692e5 * t317 * t873 + 0.3343553499e5 * t322 * t873;
double  t1363 = 0.4372339191e5 * t327 * t873 - 0.5143928460e4 * t332 * t873 - 0.1690147923e5 * t337 * t873 + 0.2571964230e4 * t342 * t873 + 0.3527265230e4 * t346 * t873 + 0.1234542830e5 * t350 * t873 + 0.1763632615e5 * t358 * t873 + 0.3086357076e5 * t362 * t873 + 0.6834076382e5 * t366 * t873 + 0.8178846251e5 * t370 * t873 + 0.4320899906e5 * t526 * t873 + 0.1432951500e5 * t493 * t873;
double  t1387 = 0.1543178538e5 * t444 * t873 + 0.9920433458e4 * t422 * t873 + 0.1859032006e4 * t530 * t409 + 0.7436128025e4 * t308 * t409 + 0.2323790008e4 * t312 * t409 - 0.3578636612e5 * t317 * t409 - 0.8133265027e5 * t322 * t409 - 0.7482603825e5 * t327 * t409 - 0.2277314208e5 * t332 * t409 + 0.1161895004e5 * t337 * t409 + 0.1068943404e5 * t342 * t409;
double  t1412 = -0.5577096019e4 * t350 * t409 - 0.4182822014e5 * t358 * t409 - 0.1352445784e6 * t362 * t409 - 0.2453922248e6 * t366 * t409 - 0.2732777049e6 * t370 * t409 - 0.1951983606e6 * t526 * t409 - 0.9759918032e5 * t493 * t409 - 0.4461676815e5 * t444 * t409 - 0.2230838407e5 * t422 * t409 - 0.8365644028e4 * t560 * t409 - 0.495e3 * t95 * t22 - 0.360e3 * t95 * t27;
double  t1427 = -0.450e3 * t95 * t45 + 0.8640e4 * t939 + 0.47520e5 * t942 + 0.1584000000e6 * t945 + 0.3564000000e6 * t908 + 0.5702400000e6 * t911 + 0.6652800000e6 * t914 + 0.5702400000e6 * t880 + 0.3564000000e6 * t883 + 0.1584000000e6 * t886 + 0.47520e5 * t893;
double  t1451 = 0.8640e4 * t905 - 0.720e3 * t95 * t38 - 0.6506612022e4 * t462 * t409 - 0.9759918032e4 * t427 * t409 - 0.5205289617e4 * t431 * t409 + 0.3794733192e2 * t704 * t111 + 0.1517893277e3 * t708 * t111 + 0.2276839915e3 * t712 * t111 + 0.1517893277e3 * t716 * t111 + 0.7348469228e3 * t406 * t873 + 0.1234542830e5 * t427 * t873 - 0.2057571384e4 * t431 * t873;
double  t1477 = 0.2629068276e3 * t704 * t226 + 0.5586770087e3 * t708 * t226 + 0.9859006035e2 * t712 * t226 - 0.4272235949e3 * t716 * t226 + 0.1971801207e4 * t607 * t226 + 0.4469416069e4 * t618 * t226 - 0.3286335345e3 * t622 * t226 - 0.6572670690e4 * t626 * t226 - 0.2629068276e4 * t630 * t226 + 0.2629068276e4 * t634 * t226 - 0.9859006035e4 * t342 * t896 - 0.1478850905e5 * t422 * t896;
double  t1497 = t46 * t2 * t33;
double  t1503 = t17 * t30 * t33;
double  t1506 = 0.9859006035e3 * t560 * t896 + 0.5258136552e4 * t454 * t226 + 0.1248807431e5 * t645 * t226 - 0.4600869483e4 * t511 * t226 - 0.2990565164e5 * t458 * t226 - 0.1150217371e5 * t653 * t226 + 0.1840347793e5 * t657 * t226 + 0.9201738966e4 * t503 * t226 - 0.6244037156e4 * t507 * t226 + 0.3943602414e4 * t1497 * t226 + 0.7887204828e4 * t536 * t226 - 0.1774621086e5 * t1503 * t226;
double  t1512 = t50 * t6 * t33;
double  t1516 = t64 * t52 * t33;
double  t1520 = t55 * t66 * t33;
double  t1524 = t16 * t57 * t33;
double  t1528 = t15 * t7 * t33;
double  t1532 = t14 * t100 * t33;
double  t1545 = -0.6309763862e5 * t1512 * t226 - 0.5866108591e5 * t1516 * t226 - 0.2070391267e5 * t1520 * t226 - 0.2415456479e5 * t1524 * t226 - 0.3746422293e5 * t1528 * t226 - 0.1626735996e5 * t1532 * t226 + 0.2957701811e4 * t1325 * t226 + 0.3943602414e3 * t389 * t896 + 0.1183080724e4 * t393 * t896 + 0.7887204828e3 * t397 * t896 - 0.7887204828e3 * t401 * t896;
double  t1570 = -0.1183080724e4 * t406 * t896 + 0.1971801207e4 * t466 * t896 + 0.5521043380e4 * t435 * t896 - 0.1380260845e5 * t481 * t896 - 0.1380260845e5 * t462 * t896 + 0.5521043380e4 * t431 * t896 + 0.2629068276e4 * t530 * t896 + 0.3943602414e4 * t308 * t896 - 0.2168981328e5 * t312 * t896 - 0.6901304225e5 * t317 * t896 - 0.6901304225e5 * t322 * t896 - 0.1380260845e5 * t327 * t896;
double  t1594 = 0.1380260845e5 * t332 * t896 - 0.1971801207e4 * t337 * t896 - 0.7887204828e4 * t350 * t896 - 0.4929503018e5 * t358 * t896 - 0.1301388797e6 * t362 * t896 - 0.1922506177e6 * t366 * t896 - 0.1863352141e6 * t370 * t896 - 0.1449273887e6 * t526 * t896 + 0.9295160031e3 * t466 * t409 + 0.4554628415e4 * t435 * t409 + 0.7807934426e4 * t439 * t409;
double  t1619 = 0.3253306011e4 * t481 * t409 - 0.1035195634e6 * t493 * t896 - 0.5619633440e5 * t444 * t896 + 0.2190890230e2 * t376 * t896 + 0.6572670690e2 * t380 * t896 + 0.6572670690e2 * t384 * t896 + 0.1669682605e4 * t618 * t111 + 0.1517893277e4 * t622 * t111 - 0.1517893277e4 * t626 * t111 - 0.3794733192e4 * t630 * t111 - 0.2580418571e4 * t634 * t111 + 0.1138419958e4 * t454 * t111;
double  t1644 = 0.2276839915e4 * t645 * t111 - 0.7968939704e4 * t511 * t111 - 0.3187575881e5 * t458 * t111 - 0.3984469852e5 * t653 * t111 - 0.1593787941e5 * t657 * t111 + 0.7968939704e4 * t503 * t111 + 0.9107359661e4 * t507 * t111 - 0.6830519746e4 * t536 * t111 - 0.4098311848e5 * t1503 * t111 - 0.9790411636e5 * t1512 * t111 - 0.1115651559e6 * t1516 * t111;
double  t1671 = -0.4781363822e5 * t1520 * t111 + 0.1593787941e5 * t1524 * t111 + 0.1593787941e5 * t1528 * t111 - 0.6830519746e4 * t1532 * t111 + 0.4898979486e1 * t13 * t320 * t1294 * t873 + 0.1469693846e3 * t376 * t873 + 0.1714642820e3 * t380 * t873 - 0.9797958971e2 * t384 * t873 + 0.1469693846e4 * t389 * t873 + 0.2057571384e4 * t393 * t873 - 0.1616663230e4 * t397 * t873 - 0.2498479538e4 * t401 * t873;
double  t1697 = 0.5878775383e4 * t466 * t873 - 0.5634891303e4 * t458 * t614 - 0.2817445652e5 * t653 * t614 - 0.1690467391e5 * t657 * t614 + 0.5634891303e4 * t503 * t614 + 0.8854829191e4 * t507 * t614 + 0.3219937888e4 * t1497 * t614 + 0.1287975155e5 * t536 * t614 + 0.7244860247e4 * t1503 * t614 - 0.4185919254e5 * t1512 * t614 - 0.9015826085e5 * t1516 * t614 - 0.6761869564e5 * t1520 * t614;
double  t1722 = -0.1126978261e5 * t1524 * t614 + 0.3219937888e4 * t1528 * t614 - 0.9659813663e4 * t1532 * t614 - 0.9659813663e4 * t1325 * t614 + 0.2494153163e4 * t1497 * t94 + 0.1496491898e5 * t536 * t94 + 0.3180045283e5 * t1503 * t94 + 0.1496491898e5 * t1512 * t94 - 0.5237721642e5 * t1516 * t94 - 0.1047544328e6 * t1520 * t94 - 0.7856582463e5 * t1524 * t94 - 0.1496491898e5 * t1528 * t94;
double  t1748 = 0.1496491898e5 * t1532 * t94 + 0.3219937888e4 * t466 * t373 + 0.1126978261e5 * t435 * t373 + 0.1014280435e5 * t439 * t373 - 0.5634891303e4 * t481 * t373 - 0.1126978261e5 * t462 * t373 + 0.6439875775e4 * t530 * t373 + 0.2575950310e5 * t308 * t373 + 0.3058940993e5 * t312 * t373 - 0.1126978261e5 * t322 * t373 + 0.2253956521e5 * t327 * t373;
double  t1773 = 0.3380934782e5 * t332 * t373 + 0.2494153163e4 * t530 * t304 - 0.7936346767e5 * t422 * t442 - 0.1763632615e5 * t560 * t442 + 0.4938171321e5 * t317 * t442 + 0.7407256982e5 * t322 * t442 - 0.4444354189e6 * t526 * t442 - 0.3703628491e6 * t493 * t442 - 0.1481451396e6 * t1528 * t168 - 0.6349077413e5 * t1532 * t168 - 0.1587269353e5 * t1325 * t168 - 0.1481451396e6 * t1516 * t168;
double  t1797 = -0.2222177095e6 * t1520 * t168 - 0.2222177095e6 * t1524 * t168 - 0.1587269353e5 * t1503 * t168 - 0.6349077413e5 * t1512 * t168 - 0.1763632615e5 * t358 * t442 + 0.2116359138e5 * t312 * t442 + 0.3219937888e3 * t634 * t292 + 0.5366563146e2 * t607 * t292 + 0.3219937888e3 * t618 * t292 + 0.1931962733e5 * t1512 * t292 + 0.1352373913e6 * t1528 * t292;
double  t1822 = 0.2414953416e4 * t1503 * t292 + 0.8049844719e3 * t622 * t292 + 0.6761869564e5 * t1532 * t292 + 0.1931962733e5 * t1325 * t292 + 0.8049844719e3 * t630 * t292 - 0.1690467391e5 * t458 * t292 - 0.8049844719e3 * t645 * t292 + 0.6761869564e5 * t1516 * t292 + 0.1352373913e6 * t1520 * t292 + 0.1018233765e4 * t1497 * t241 + 0.1221880518e6 * t1528 * t241 + 0.4582051942e5 * t1532 * t241;
double  t1847 = 0.1018233765e5 * t1325 * t241 + 0.2138290906e6 * t1516 * t241 + 0.1073312629e4 * t626 * t292 + 0.2565949088e6 * t1520 * t241 + 0.2138290906e6 * t1524 * t241 + 0.1018233765e5 * t536 * t241 + 0.4582051942e5 * t1503 * t241 + 0.1221880518e6 * t1512 * t241 - 0.5634891303e4 * t507 * t292 + 0.1690467391e6 * t1524 * t292 + 0.2057571384e5 * t653 * t168;
double  t1872 = 0.2351510153e4 * t645 * t168 - 0.2817445652e5 * t653 * t292 - 0.2817445652e5 * t657 * t292 - 0.1690467391e5 * t503 * t292 - 0.5634891303e4 * t511 * t292 - 0.3703628491e6 * t370 * t442 + 0.7407256982e5 * t327 * t442 + 0.1646057107e5 * t657 * t168 + 0.1774621086e5 * t560 * t420 + 0.1656313014e6 * t444 * t420 + 0.9976612652e4 * t1325 * t94 + 0.3968173383e5 * t893 * t168;
double  t1898 = 0.8818163074e4 * t905 * t168 + 0.1971801207e4 * t9 * t420 - 0.4704239994e6 * t914 * t241 - 0.3360171424e6 * t911 * t241 - 0.2897944099e5 * t893 * t292 - 0.1014280435e6 * t880 * t292 - 0.1014280435e6 * t883 * t292 - 0.2897944099e5 * t911 * t292 + 0.1018233765e4 * t929 * t241 - 0.8049844719e3 * t945 * t292 - 0.7244860247e4 * t908 * t292 - 0.1120057141e5 * t942 * t241;
double  t1923 = -0.5600285707e5 * t945 * t241 - 0.1680085712e6 * t908 * t241 - 0.1120057141e5 * t905 * t241 + 0.8729536070e4 * t439 * t304 + 0.1454922678e5 * t481 * t304 + 0.1454922678e5 * t462 * t304 + 0.8729536070e4 * t427 * t304 + 0.2909845357e4 * t431 * t304 + 0.3219937888e4 * t454 * t614 + 0.1368473602e5 * t645 * t614 + 0.1690467391e5 * t511 * t614 - 0.720e3 * t95 * t12;
double  t1929 = t825 + t353 + t1671 + t1923 + t1644 + t479 + t514 + t259 + t1233 + t106 + t447 + t1545 + t1081 + t1451 + t1055 + t599 + t1412 + t1773 + t1283 + t404 + t1106 + t573 + t926 + t1748 + t1594 + t1005 + t1310 + t1363 + t889 + t1847 + t1338 + t1872 + t1159 + t222 + t1477 + t1030 + t295 + t1697 + t761 + t797 + t1822 + t1619 + t1387 + t956 + t856 + t1570 + t1722 + t186 + t1208 + t637 + t727 + t1258 + t1134 + t980 + t545 + t1898 + t1427 + t1183 + t695 + t141 + t1797 + t63 + t1506 + t670;


return t1929;

}


double evalNGPDF(double j, double *a){



double t1 = pow(a[0], 0.2e1);
double t2 = pow(a[1], 0.2e1);
double t3 = pow(a[2], 0.2e1);
double t4 = pow(a[3], 0.2e1);
double t5 = pow(a[4], 0.2e1);
double t7 = sqrt(0.1e1 - t1 - t2 - t3 - t4 - t5);
double t8 = sqrt(0.2e1);
double t9 = sqrt(M_PI);
double t11 = sqrt(t9 * t8);
double t12 = 0.1e1 / t11;
double t14 = j * j;
double t16 = exp(-t14 / 0.4e1);
double t24 = sqrt(0.3e1);
double t27 = t8 * t14 * j;
double t29 = t8 * j;
double t36 = sqrt(0.6e1);
double t38 = t14 * t14;
double t46 = sqrt(0.15e2);
double t58 = sqrt(0.5e1);
double t70 = pow(t16 * t12 * t7 + t16 * (-0.1e1 + t14) * t12 * t8 * a[0] / 0.2e1 + t16 * (0.2e1 * t27 - 0.6e1 * t29) * t12 * t24 * a[1] / 0.12e2 + t16 * (0.12e2 + 0.4e1 * t38 - 0.24e2 * t14) * t12 * t36 * a[2] / 0.48e2 + t16 * (0.4e1 * t8 * t38 * j - 0.40e2 * t27 + 0.60e2 * t29) * t12 * t46 * a[3] / 0.240e3 + t16 * (-0.120e3 + 0.8e1 * t38 * t14 - 0.120e3 * t38 + 0.360e3 * t14) * t12 * t58 * a[4] / 0.480e3, 0.2e1);


return t70;

}
