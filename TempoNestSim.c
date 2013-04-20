#include <stdio.h>
#include <vector>
#include <string.h>
#include "/usr/include/gsl/gsl_sf_gamma.h"
#include <gsl/gsl_multimin.h>
#include "dgemm.h"
#include "dgemv.h"
#include "dpotri.h"
#include "dpotrf.h"
#include "dpotrs.h"
#include "tempo2.h"
#include "T2toolkit.h"
#include "TKfit.h"
#include "TempoNest.h"


void callFit(pulsar *psr,int npsr){
  int iteration;
  double globalParameter = 0.0;

  for (iteration=0;iteration<2;iteration++)
    {
      formBatsAll(psr,npsr);
            
      /* Form residuals */
      formResiduals(psr,npsr,0);
      
      /* Do the fitting */
      if (iteration==0) doFit(psr,npsr,0);
      else textOutput(psr,npsr,globalParameter,0,0,0,"");
    }
}


void doSim(int argc, char **argv, pulsar *psr, char timFile[][MAX_FILELEN], char parFile[][MAX_FILELEN]){

	int doRed=0;
	int updateEFAC;
	int updateEQUAD;
	double redamp=-1000;
	double redalpha=-1000;
	double EFAC=1;
	double logequad=-1000;
	double EQUAD=0;
	int doDM=0;
	double DMamp=-1000;
	double DMalpha=-1000;
	long idum=-1;

	printf("\n\n Entering TempoNest Simulation Mode.\n\n");
	printf("Allows for the creation of simulated data using existing par/tim files\n");
	printf("Includes support for EFAC/EQUAD, red noise and DM additions to existing toa error bars\n");
	printf("For more information use TempoNest -sim -h\n");


	for (int i=0;i<argc;i++){
		if(strcmp(argv[i],"-incred")==0){
			doRed=1;
			printf("Including Red Noise in simulation\n");
		}
		if(strcmp(argv[i],"-redamp")==0){
			 sscanf(argv[i+1],"%lf",&redamp);
			printf("Including Red Noise Amplitude %g\n",redamp);
		}
		if(strcmp(argv[i],"-redindex")==0){
			sscanf(argv[i+1],"%lf",&redalpha);
			printf("Including Red Noise Index %g\n",redalpha);
		}
		if(strcmp(argv[i],"-incDM")==0){
			doDM=1;
			printf("Including DM in simulation\n");
		}
		if(strcmp(argv[i],"-DMamp")==0){
			 sscanf(argv[i+1],"%lf",&DMamp);
			printf("Including DM Amplitude %g\n",DMamp);
		}
		if(strcmp(argv[i],"-DMindex")==0){
			sscanf(argv[i+1],"%lf",&DMalpha);
			printf("Including DM Index %g\n",DMalpha);
		}
		if(strcmp(argv[i],"-efac")==0){
			sscanf(argv[i+1],"%lf",&EFAC);
			printf("Including EFAC %g\n",EFAC);
		}
		if(strcmp(argv[i],"-equad")==0){
			sscanf(argv[i+1],"%lf",&logequad);
			EQUAD=pow(10.0, logequad);
			printf("Including log_10 EQUAD in simulation\n",logequad);
		}
		if(strcmp(argv[i],"-updateefac")==0){
			
			updateEFAC=1;
			printf("Updating TOA errors to reflect included EFAC\n");
		}
		if(strcmp(argv[i],"-updateequad")==0){
			
			updateEQUAD=1;
			printf("Updating TOA errors to reflect included EQUAD\n");
		}
		if(strcmp(argv[i],"-seed")==0){
			sscanf(argv[i+1],"%d",&idum);
			printf("Using seed %d\n",idum);
		}
		if(strcmp(argv[i],"-h")==0){
		printf("==========================================================================================\n");
		printf(" Sim TempoNest - usage instructions.\n");
		printf(" TempoNest -sim [options] -f file.par file.tim.\n");
		printf("\n options:\n");
		printf(" \t -seed X : specify seed X for random numbers, default=system time\n");
		printf(" \t -incred: include power law red noise of the form S(f) = A^2.1yr^3/(12Pi^2) . (f/1yr^-1)^-y\n");
		printf(" \t -redamp xxx: specify log_10 amplitude of red noise\n");
		printf(" \t -redindex xxx: specify spectral index of red noise y.\n");
		printf(" \t -incDM: include power law DM model of the form S(f) = A^2.1yr^3/(12Pi^2) . (f/1yr^-1)^-y\n");
		printf(" \t -DMamp xxx: specify log_10 amplitude of DM\n");
		printf(" \t -DMindex xxx: specify spectral index of DM y.\n");
		printf(" \t -efac xxx: include global EFAC value for toa error bars of value xxx\n");
		printf(" \t -equad xxx: specify value of log_10 amplitude EQUAD value xxx\n");
		printf(" \t -updateefac : update toa errors to reflect inclusion of EFAC\n");
		printf(" \t -updateequad : update toa errors to reflect inclusion of EQUAD\n");
		printf("\n\n Have a nice day!\n\n");
		printf("==========================================================================================\n");
		exit(0);
		}
	}

	if(doRed ==1 && redamp == -1000){
		printf("Set to include red noise, but no amplitude given\n");
		exit(0);
	}
	if(doRed ==1 && redalpha == -1000){
		printf("Set to include red noise, but no index given\n");	
		exit(0);
	}

	if(doDM ==1 && DMamp == -1000){
		printf("Set to include DM, but no amplitude given\n");
		exit(0);
	}
	if(doDM ==1 && DMalpha == -1000){
		printf("Set to include DM, but no index given\n");
		exit(0);
	}

	TNSimRedfromTim(argc, argv, psr, timFile, parFile, EFAC, EQUAD, doRed, redamp, redalpha, updateEFAC, updateEQUAD, doDM, DMamp, DMalpha,idum);

}


void TNSimRedfromTim(int argc, char **commandLine, pulsar *psr, char timFile[][MAX_FILELEN], char parFile[][MAX_FILELEN], double EFAC, double EQUAD, int doRed, double redlogamp, double redslope, int updateEFAC, int updateEQUAD, int doDM, double DMlogamp, double DMslope,long idum)
{

	char str[MAX_FILELEN];

	readParfile(psr,parFile,timFile,1); // Load the parameters
	readTimfile(psr,timFile,1); /* Read .tim file to define the site-arrival-times */
	preProcess(psr,1,argc,commandLine);
	callFit(psr,1);             /* Do all the fitting routines */

	int nit=4;
	strcpy(str,parFile[0]);
	str[strlen(str)-4]='\0';
	strcat(str,".simulate");
	strcpy(timFile[0],str);

	for (int j=0;j<nit;j++)
	{
	  /* Now update the site arrival times depending upon the residuals */
	  
	  for (int i=0;i<psr[0].nobs;i++)  
	    {
	      psr[0].obsn[i].sat -= psr[0].obsn[i].prefitResidual/SECDAY; 
	      //psr->obsn[i].nFlags = 0;
	    } 
	  writeTim(str,psr,"tempo2");
	  //	  initialise(&psr[ii],0);
	  // Reset the jumps
	  psr[0].nJumps = 0;
	  for(int kk=0;kk<MAX_JUMPS;kk++){
	      psr[0].jumpVal[kk] = 0.0;
	      psr[0].jumpValErr[kk] = 0.0;
	  }
	  for(int jj=0;jj<MAX_PARAMS;jj++){
	      psr[0].param[jj].nLinkTo = 0;
	      psr[0].param[jj].nLinkFrom = 0;
	  }
	  readParfile(psr,parFile,timFile,1); /* Load the parameters       */
	  readTimfile(psr,timFile,1); 
	  preProcess(psr,1,argc,commandLine);
	  /* Now run the superTEMPO code again */
	  callFit(psr,1);             /* Do all the fitting routines */
	}

     	 printf("Complete %d iterations %s\n",nit, timFile[0]);

	for (int i=0;i<psr[0].nobs;i++)  
	{
	  psr[0].obsn[i].sat -= psr[0].obsn[i].prefitResidual/SECDAY;  
	 // psr->obsn[i].nFlags = 0;
	}

	formBatsAll(psr,1);
	double *bats=new double[psr->nobs];

	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;
	





	double timdiff=0;

/* ************************Calculate DM Matrix if to be included ************** */

	double DMconst=0;
	double DMampsquared=0;
	double DMalpha=0;
	if(doDM ==1){
		DMalpha=DMslope;
		double DMamp=pow(10.0,DMlogamp);
		DMampsquared=DMamp*DMamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-DMalpha)))/(pow(flo,(DMalpha-1)));
		DMconst=gsl_sf_gamma(1-DMalpha)*sin(0.5*M_PI*DMalpha);
	}

	double **DMMatrix = new double*[psr->nobs]; for(int o1=0;o1<psr->nobs;o1++)DMMatrix[o1]=new double[psr->nobs];

	for(int o1=0;o1<psr->nobs; o1++){
		for(int o2=0; o2< psr->nobs; o2++){

			timdiff=(double)psr->obsn[o1].bat - (double)psr->obsn[o2].bat;	
			double tau=2.0*M_PI*fabs(timdiff);
			double DMsum=0;

			for(int k=0; k <=10; k++){
				DMsum=DMsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-DMalpha));

			}
			if(doDM ==1){
			DMMatrix[o1][o2]=DMampsquared*(DMconst*pow((flo*tau),(DMalpha-1)) - DMsum);
			}
			else if(doDM == 0){
			DMMatrix[o1][o2]=0;
			}
		}
	}

	if(doDM ==1){
		//Get vector of frequencies
		double *DMVec = new double[psr->nobs]; 
		double DMKappa = 2.410*pow(10.0,-16);
		for(int o1=0;o1<psr->nobs; o1++){
			DMVec[o1]=1.0/(DMKappa*pow((double)psr->obsn[o1].freqSSB,2));
			printf("freqs: %i %g %g %g \n",o1,(double)psr->obsn[o1].freq,(double)psr->obsn[o1].freqSSB,DMVec[o1]);
		}

		//multiply matrices
		//First DM*FMatrix
		for(int o1=0;o1<psr->nobs; o1++){
			for(int o2=0; o2< psr->nobs; o2++){
				DMMatrix[o1][o2] = DMMatrix[o1][o2]*DMVec[o2];
			}
		}
		//Then FMatrix*DM
		for(int o1=0;o1<psr->nobs; o1++){
			for(int o2=0; o2< psr->nobs; o2++){
				DMMatrix[o1][o2] = DMVec[o1]*DMMatrix[o1][o2];
			}
		}	
	}

/* ************************Calculate Red Matrix if to be included ************** */

	double Redconst=0;
	double Redampsquared=0;
	double Redalpha=0;
	if(doRed ==1){
		Redalpha=redslope;
		double Redamp=pow(10.0,redlogamp);
		Redampsquared=Redamp*Redamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-Redalpha)))/(pow(flo,(Redalpha-1)));
		Redconst=gsl_sf_gamma(1-Redalpha)*sin(0.5*M_PI*Redalpha);
		
	}

	double **RedMatrix = new double*[psr->nobs]; for(int o1=0;o1<psr->nobs;o1++)RedMatrix[o1]=new double[psr->nobs];

	for(int o1=0;o1<psr->nobs; o1++){
		for(int o2=0; o2< psr->nobs; o2++){

			timdiff=(double)psr->obsn[o1].bat - (double)psr->obsn[o2].bat;	
			double tau=2.0*M_PI*fabs(timdiff);
			double Redsum=0;

			for(int k=0; k <=10; k++){
				Redsum=Redsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(iter_factorial(2*k)*(2*k+1-Redalpha));

			}
			if(doRed ==1){
			RedMatrix[o1][o2]=Redampsquared*(Redconst*pow((flo*tau),(Redalpha-1)) - Redsum);
			}
			else if(doRed == 0){
			RedMatrix[o1][o2]=0;
			}
		}
	}

/* ************************Calculate total Cov Matrix ************** */

	double **CovMatrix = new double*[psr->nobs]; for(int o1=0;o1<psr->nobs;o1++)CovMatrix[o1]=new double[psr->nobs];

	for(int o1=0;o1<psr->nobs; o1++){

		for(int o2=0; o2< psr->nobs; o2++){


			CovMatrix[o1][o2]=DMMatrix[o1][o2]+RedMatrix[o1][o2];


			if(o1==o2){
				double whitenoise = pow(psr->obsn[o1].toaErr*pow(10.0,-6)*EFAC,2) + EQUAD;
				CovMatrix[o1][o2] += whitenoise;
				if(updateEFAC==1 && updateEQUAD==0){
					psr[0].obsn[o1].origErr=psr->obsn[o1].toaErr*EFAC;
				}
				if(updateEFAC==0 && updateEQUAD==1){
					psr[0].obsn[o1].origErr=sqrt(pow(psr->obsn[o1].toaErr*pow(10.0,-6),2) + EQUAD)*pow(10.0,6);
				}
				if(updateEFAC==1 && updateEQUAD==1){
					psr[0].obsn[o1].origErr=sqrt(whitenoise)*pow(10.0,6);
				}
				printf("%i %i %g \n",o1,o2,whitenoise);
			}
			

		}
	}

	for(int o1=0;o1<psr->nobs; o1++){
		for(int o2=0; o2 < o1; o2++){
			CovMatrix[o1][o2]=0;
		}
	}

	double covdet=0;
	dpotrfU(CovMatrix, psr->nobs, covdet);
	long seed = 1;
	if(idum==-1){
		seed=TKsetSeed();
	}
	else{
		seed=idum;
	}
	double *RanVec = new double[psr->nobs];
	double *NoiseVec = new double[psr->nobs];
	for(int o1=0;o1 < psr->nobs; o1++){
		
		RanVec[o1]=TKgaussDev(&seed);
	}

	dgemv(CovMatrix,RanVec,NoiseVec,psr->nobs,psr->nobs,'T');
	double sum=0;
	for(int o1=0;o1 < psr->nobs; o1++){
		sum+=NoiseVec[o1];		
	}

	for(int o1=0;o1 < psr->nobs; o1++){
		NoiseVec[o1]-= sum/psr->nobs;
		bats[o1]=(double)psr[0].obsn[o1].bat;		
	}

	TKremovePoly_d(bats,NoiseVec,psr[0].nobs,2); // remove a quadratic to reduce the chances of phase wraps

	for(int o1=0;o1 < psr->nobs; o1++){
		psr[0].obsn[o1].sat +=	NoiseVec[o1]/(24*60*60);
		printf("NoiseVec: %i %g \n",o1,NoiseVec[o1]);	
	}

	for(int o=0;o<psr->nobs;o++){delete[] CovMatrix[o];}
	delete[] CovMatrix;
	delete[] RanVec;
	delete[] NoiseVec;
	//printf("Like: %g %g %g \n",lnew,Chisq,covdet);

	writeTim(str, psr, "tempo2");

}



