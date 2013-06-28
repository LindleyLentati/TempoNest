#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
// Copyright (C) 2013 Lindley Lentati

/*
* This file is part of TempoNest
*
* TempoNest is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* TempoNest is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with TempoNest. If not, see <http://www.gnu.org/licenses/>.
*/

/*
* If you use TempoNest and as a byproduct both Tempo2 and MultiNest
* then please acknowledge it by citing Lentati L., Alexander P., Hobson M. P. (2013) for TempoNest,
* Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2,
* pp. 655-672 (bibtex: 2006MNRAS.369..655H)
* or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
* pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
* timing model and MultiNest Papers here.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "tempo2.h"
#include "tempo2Util.h"
#include "tempo2pred.h"
#include "tempo2pred_int.h"
#include "T2accel.h"
#include <dlfcn.h>
// #include "T2toolkit.h"
#include <vector>
#include <algorithm>
#include <float.h>
#include "multinest.h"
#include "TempoNest.h"

#include <string.h>
#include <sstream>
#include <iterator>
#include <cstring>
#include <iostream>
#include <fstream>

#include <gsl/gsl_sf_gamma.h>
#include "dpotri.h"
#include "dpotrf.h"

#ifdef HAVE_CULA
#include <cula_lapack_device.h>
#endif /* HAVE_CULA */

void ephemeris_routines(pulsar *psr,int npsr);
void clock_corrections(pulsar *psr,int npsr);
void extra_delays(pulsar *psr,int npsr);

#ifdef HAVE_CULA
extern "C" void copy_gmat_(double *G, int N);
extern "C" void copy_floatgmat_(float *G, int N);
#endif /* HAVE_CULA */


void fastephemeris_routines(pulsar *psr,int npsr)
{
vectorPulsar(psr,npsr); /* 1. Form a vector pointing at the pulsar */
//readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
//tt2tb(psr,npsr); /* Observatory/time-dependent part of TT-TB */
//readEphemeris(psr,npsr,0); /* Re-evaluate ephemeris with correct TB */

}

void fastformBatsAll(pulsar *psr,int npsr)
{
    //clock_corrections(psr,npsr); /* Clock corrections ... */
    fastephemeris_routines(psr,npsr); /* Ephemeris routines ... */
  	extra_delays(psr,npsr); /* Other time delays ... */
	formBats(psr,npsr); /* Form Barycentric arrival times */
	secularMotion(psr,npsr);

}

MNStruct* init_struct(pulsar *pulseval,	 long double **LDpriorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval,int numFitEFACval, int numFitEQUADval, int numFitCoeffval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval, int numdimsval, int incREDval, int incDMval, int TimeMarginVal, int JumpMarginVal, int doLinearVal)
{
    MNStruct* MNS = (MNStruct*)malloc(sizeof(MNStruct));

	MNS->pulse=pulseval;
	MNS->LDpriors=LDpriorsval;
	MNS->numberpulsars=numberpulsarsval;
	MNS->numFitJumps=numFitJumpsval;
	MNS->numFitTiming=numFitTimingval;
	MNS->numFitEFAC=numFitEFACval;
	MNS->numFitEQUAD=numFitEQUADval;
	MNS->numFitRedCoeff=numFitCoeffval;
	MNS->TempoFitNums=TempoFitNumsval;
	MNS->TempoJumpNums=TempoJumpNumsval;
	MNS->sysFlags=sysFlagsval;
	MNS->numdims=numdimsval;
	MNS->incRED=incREDval;
	MNS->incDM=incDMval;
	MNS->TimeMargin=TimeMarginVal;
	MNS->JumpMargin=JumpMarginVal;
	MNS->doLinear=doLinearVal;

	return MNS;
}


void update_MNPriors(MNStruct* MNS, double ** DPriorsval, long double **priorsval, int linearPriors){


	int paramsfitted=1;

	
	if(linearPriors != 2){
	//printf("lin p is 0\n");
	for (int p=0;p<MAX_PARAMS;p++) {
        	for (int k=0;k<MNS->pulse->param[p].aSize;k++){
                	if(MNS->pulse->param[p].fitFlag[k] == 1){
				if(p == param_ecc || p ==param_px || p == param_m2 || p==param_dm){
					long double minprior=priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1];
					if(minprior < 0 && priorsval[paramsfitted][2]==0){
						printf("%.10Lg %.10Lg %g %g \n",priorsval[paramsfitted][0],priorsval[paramsfitted][1],DPriorsval[paramsfitted+linearPriors][0],DPriorsval[paramsfitted+linearPriors][1]);
						long double newprior=-priorsval[paramsfitted][0]/priorsval[paramsfitted][1];
						DPriorsval[paramsfitted+linearPriors][0]=(double)newprior;
	                        		printf("Prior on %s updated to be physical (was <0) : %.25Lg -> %.25Lg\n", MNS->pulse->param[p].shortlabel[k], priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][0]*priorsval[paramsfitted][1],priorsval[paramsfitted][0]+DPriorsval[paramsfitted+linearPriors][1]*priorsval[paramsfitted][1]);
	                        		printf("%.10Lg %.10Lg %g %g \n",priorsval[paramsfitted][0],priorsval[paramsfitted][1],DPriorsval[paramsfitted+linearPriors][0],DPriorsval[paramsfitted+linearPriors][1]);
					}
				}
                                paramsfitted++;
			 }
                }
       }

	}
	

	MNS->Dpriors=DPriorsval;

	MNS->LDpriors=priorsval;

}




	void update_MNGdata(MNStruct* MNS, int Gsizeval, double **GMatrixval)
	{

		MNS->Gsize=Gsizeval;
		MNS->GMatrix=GMatrixval;


	}

	void update_MNDdata(MNStruct* MNS, int Dsizeval, double **DMatrixval)
	{

  
		MNS->Dsize=Dsizeval;
		MNS->DMatrix=DMatrixval;

	}

	void update_MNGDdata(MNStruct* MNS, int Dsizeval, double **DMatrixval,int Gsizeval, double **GMatrixval)
	{


		MNS->Dsize=Dsizeval;
		MNS->DMatrix=DMatrixval;
		MNS->Gsize=Gsizeval;
		MNS->GMatrix=GMatrixval;



	}



//
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context)
{
// convert the 2D Fortran arrays to C arrays


// the posterior distribution
// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns

int i, j;

double postdist[nSamples][nPar + 2];
for( i = 0; i < nPar + 2; i++ )
for( j = 0; j < nSamples; j++ )
postdist[j][i] = posterior[0][i * (nSamples) + j];



// last set of live points
// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column

double pLivePts[nlive][nPar + 1];
for( i = 0; i < nPar + 1; i++ )
for( j = 0; j < nlive; j++ )
pLivePts[j][i] = physLive[0][i * (nlive) + j];
}



/* The main function of a plugin called from Tempo2 is 'graphicalInterface'
*/
int main(int argc, char *argv[])
{
	int iteration; 
	int listparms;
	int outRes=0;
	int writeModel=0;
	char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
	char outputSO[MAX_FILELEN];
	char str[MAX_FILELEN];
	char newparname[MAX_FILELEN];
	longdouble coeff[MAX_COEFF]; /* For polynomial coefficients in polyco */
	int npsr;      /* The number of pulsars */
	int noWarnings=1;
	double globalParameter=0.0;
	int  displayParams,p;
	int nGlobal,i,flagPolyco=0,it,k;
	char polyco_args[128];
	int newpar=0;
	int onlypre=0;
	//  char tempo2MachineType[MAX_FILELEN]="";
	FILE *alias;
	char **commandLine;
	clock_t startClock,endClock;
	const char *CVS_verNum = "$Revision: 1.28 $";
	int numFitJumps;
	int numToMargin=0;


  printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome to redistribute it\n");
  printf("under conditions of GPL license.\n\n");

  //
  startClock = clock();

  pulsar *psr;

  commandLine = (char **)malloc(1000*sizeof(char *));

  for (i=0;i<1000;i++)
    commandLine[i] = (char *)malloc(sizeof(char)*1000);

  /* Parse input line for machine type */
  for (i=0;i<argc;i++)    
    {
      if (strcasecmp(argv[i],"-npsr")==0)
	sscanf(argv[i+1],"%d",&MAX_PSR);
      else if (strcasecmp(argv[i],"-displayVersion")==0)
	displayCVSversion = 1;
      else if (strcasecmp(argv[i],"-nobs")==0)
	sscanf(argv[i+1],"%d",&MAX_OBSN);
      else if (strcasecmp(argv[i],"-debug")==0){
	debugFlag=1;
	tcheck=1;
	writeResiduals=1;
	  }
	  else if (strcasecmp(argv[i],"-tcheck")==0)
	tcheck=1;
	  else if (strcasecmp(argv[i],"-noaccel")==0)
	useT2accel=0;
	  else if (strcasecmp(argv[i],"-writeres")==0)
	writeResiduals=1;
      else if (strcasecmp(argv[i],"-veryfast")==0)
	veryFast=1;

      strcpy(commandLine[i],argv[i]);
    }
  if (displayCVSversion == 1) CVSdisplayVersion("tempo2.C","main()",CVS_verNum);

  if ((psr = (pulsar *)malloc(sizeof(pulsar)*MAX_PSR))==NULL)
    {
      printf("Not enough memory to allocate room for %d pulsars\n",MAX_PSR);
      printf("Please decrease the value of MAX_PSR_VAL in tempo2.h\n"); 
      exit(1); 
    }
  logdbg("Have allocated memory for pulsar");
  psr[0].jboFormat = 0;

  for (i=1;i<argc;i++)
    {
      if (strcmp(commandLine[i],"-machine")==0)
	strcpy(tempo2MachineType,commandLine[++i]);
      else if (strcasecmp(commandLine[i],"-noWarnings")==0)
	noWarnings=2;
      else if (strcasecmp(commandLine[i],"-allInfo")==0)
	noWarnings=0;
      else if (strcmp(commandLine[i],"-jbo")==0)
	psr[0].jboFormat=1;
      else if (strcmp(commandLine[i],"-test")==0) /* Use TEMPO2_TEST environment variable */
	strcpy(TEMPO2_ENVIRON,"TEMPO2_TEST");
      else
	{
	  char oldCommandLine[1000][1000];
	  int  oldArgc = argc;
	  int  oldI = i;

	  for (k=i;k<argc;k++)
	    strcpy(oldCommandLine[k-i],commandLine[k]);

	  sprintf(str,"%s/alias.dat",getenv(TEMPO2_ENVIRON));	  
	  if ((alias = fopen(str,"r")))
	    {
	      while (!feof(alias))
		{
		  char *ret;
		  fgets(str,MAX_FILELEN,alias);
		  ret = strtok(str," ");
		  if (strcmp(ret,commandLine[i])==0)
		    {
		      printf("Using alias: %s\n",ret);
		      /* Have an alias: now split up the remainder of the line */
		      while ((ret = strtok(NULL," "))!=NULL)
			{
			  strcpy(commandLine[i],ret);
			  if (commandLine[i][strlen(commandLine[i])-1]=='\n')
			    commandLine[i][strlen(commandLine[i])-1] = '\0';
			  i++;
			  argc++;
			}
		      argc--;
		      for (k=1;k<oldArgc-oldI;k++)
			strcpy(commandLine[i+k-1],oldCommandLine[k]);
		      i=oldI-1;

		    }
		}
	      fclose(alias);
	    }
	}
    }
  int ii;
  for(ii=1;ii<argc;ii++){
      if (strcasecmp(commandLine[ii],"-reminder")==0){
	  // Writing command line to log file
	  char commandfile[200] = "T2command.input";
	  FILE *fout;
	  time_t rawtime;
	  struct tm * timeinfo;
	  time ( &rawtime );
	  timeinfo = localtime ( &rawtime );
	  char timeval[200];
	  strcpy(timeval,asctime (timeinfo));
	  strcpy(&timeval[(int)strlen(timeval)-1],"");
	  fout = fopen(commandfile,"a");
	  fprintf(fout,"[%s]>> ",timeval);
	  for(i=0;i<argc;i++){
	      fprintf(fout," %s ",commandLine[i]);
	  }
	  fprintf(fout,"\n");
	  fclose(fout);
      }
  }

  if (getenv(TEMPO2_ENVIRON)==NULL)
    {
      printf("Environment variable >%s< not set\n",TEMPO2_ENVIRON);
      exit(1);
    }

  /* get path to look for plugins */
  setPlugPath();


  /* If running from the command line ... */
  logdbg("Running initialise");
  initialise(psr,noWarnings); /* Initialise all */
  logdbg("Completed running initialise %d",psr[0].nits);
  /* Obtain login architecture */
  if (strlen(tempo2MachineType)==0)
    {
#ifdef  TEMPO2_ARCH 
      strcpy(tempo2MachineType, TEMPO2_ARCH);
#else
      if (getenv("LOGIN_ARCH")==NULL)
	{
	  printf("Unable to determine machine type: You must do one of the following:\n"
                 "Re-compile tempo2 with the standard export distrubution, or\n"
                 "Set the LOGIN_ARCH environment variable, or\n"
                 "Use -machine on the command line\n");
	  exit(1);
	}
      strcpy(tempo2MachineType, getenv("LOGIN_ARCH"));
#endif
    }

  if (sizeof(longdouble)!=16 && noWarnings<1)
    {
      printf("Warning: the size of a longdouble is only %d bytes\n",sizeof(longdouble));
      printf(" --- the size of a double is %d bytes\n",sizeof(double));
    }
  strcpy(outputSO,"");
  if (argc==1) /* No command line arguments */
    {
      printf("%s\n",PACKAGE_STRING);
      printf("  Usage:     %s -f XXX.par XXX.tim\n",argv[0]);
      printf("Plugin search paths:\n");
      for (i=0; i < tempo2_plug_path_len; i++){
	      printf(" -- %s/*.t2\n",tempo2_plug_path[i]);
      }
	  printf("\n");
#ifdef HAVE_LAPACK
	  printf("* Using LAPACK acceleration for Cholesky decomposition\n");
#endif
#ifdef HAVE_BLAS
	  printf("* Using BLAS acceleration for matrix mulitplication\n");
#endif
      printf("\nFor more help, use %s -h\n",argv[0]);
      exit(1);
    }
  npsr = 0;   /* Initialise the number of pulsars */
  displayParams=0;
  nGlobal=0;
  /* Obtain command line arguments */
  logdbg("Running getInputs %d",psr[0].nits);
  getInputs(psr,argc, commandLine, timFile,parFile,&listparms,&npsr,&nGlobal,&outRes,&writeModel,
	    outputSO,&flagPolyco,polyco_args,&newpar,&onlypre,dcmFile,covarFuncFile,newparname);
  logdbg("Completed getInputs");

	int doMaxLike=0;
  for (i=1;i<argc;i++)
    {
	if (strcmp(commandLine[i],"-sim")==0 ){ // Doing Sim?
		//printf("Doing sim\n");
		doSim(argc, commandLine, psr, timFile, parFile);
		return 0;
	}

	if (strcmp(commandLine[i],"-maxlike")==0 ){ // Doing Sim?
		//printf("Doing sim\n");
		doMaxLike=1;
	}

      if (strcmp(commandLine[i],"-gr")==0 || strcmp(commandLine[i],"-gr2")==0) 
	/* Running from a graphical interface? */
	{      
	  char *(*entry)(int,char **,pulsar *,int *);
	  void * module;

	  if (strcmp(commandLine[i],"-gr2")==0){
		  sprintf(str,"./%s_%s_plug.t2",commandLine[i+1],tempo2MachineType);
		  printf("Looking for %s\n",str);
		  module = dlopen(str, RTLD_NOW);
	  } else{
		  for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
			  sprintf(str,"%s/%s_%s_plug.t2",tempo2_plug_path[iplug],
					  commandLine[i+1],tempo2MachineType);
			  printf("Looking for %s\n",str);
			  module = dlopen(str, RTLD_NOW); 
			  if(module==NULL){	  
				  printf("dlerror() = %s\n",dlerror());
			  } else break;
		  }
	  }
	  if(!module)  {
	    fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
//	    fprintf(stderr, "dlerror() = %s\n",dlerror());
	    return -1;
	  }

	  /*
	   * Check that the plugin is compiled against the same version of tempo2.h
	   */
	  char ** pv  = (char**)dlsym(module, "plugVersionCheck");
	  if(pv!=NULL){
		  // there is a version check for this plugin
		  if(strcmp(TEMPO2_h_VER,*pv)){
			  fprintf(stderr, "[error]: Plugin version mismatch\n");
			  fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
			  fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
			  dlclose(module);
			  return -1;
		  }
	  }

	  entry = (char*(*)(int,char **,pulsar *,int *))dlsym(module, "graphicalInterface");
	  if( entry == NULL ) {
	    dlclose(module);
	    fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
	    fprintf(stderr, "dlerror() = %s\n",dlerror());
	    return -1;
	  }
	  logdbg("--ENTER GRAPHICAL PLUGIN--");
	  entry(argc,commandLine,psr,&npsr);
	  return 0;
	}
    }
	
  logdbg("Reading par file");
  readParfile(psr,parFile,timFile,npsr); /* Read .par file to define the pulsar's initial parameters */  
  logdbg("Finished reading par file %d",psr[0].nits);
  if (flagPolyco==0)
    {
      logdbg("Running readTimfile");
      readTimfile(psr,timFile,npsr); /* Read .tim file to define the site-arrival-times */
      logdbg("Completed readTimfile %d",psr[0].param[param_ecc].paramSet[1]);
    }

  logdbg("Running preProcess %d",psr[0].nits);
  preProcess(psr,npsr,argc,commandLine);
  logdbg("Completed preProcess %d",psr[0].nits);
  if (flagPolyco> 0)  /* Running tempo2 in polyco mode? */
  {
    if (flagPolyco == 1)
    {
      longdouble mjd1, mjd2,maxha,freq;
      int ncoeff,nspan; 
      double seg_length; 
      char sitename[128];
      if (sscanf(polyco_args, "%Lf %Lf %d %d %Lf %s %Lf", &mjd1, &mjd2, &nspan,
		 &ncoeff, &maxha,sitename,&freq)!=7)
      {
	fprintf(stderr, "Error parsing -polyco arguments! See tempo2 -h.\n");
	printf("Have: %s\n",polyco_args);
	exit(1);
      }
      if (psr[0].tempo1 == 0)
	printf("WARNING: Should probably use -tempo1 option\n");

      if (psr[0].param[param_tzrmjd].paramSet[0]==0)
	{
	  printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)psr[0].param[param_pepoch].val[0]);
	  psr[0].param[param_tzrmjd].paramSet[0]=1;
	  psr[0].param[param_tzrmjd].val[0] = psr[0].param[param_pepoch].val[0];
	}
      if (psr[0].param[param_tzrfrq].paramSet[0]==0)
	{
	  printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)1400.0);
	  psr[0].param[param_tzrfrq].paramSet[0]=1;
	  psr[0].param[param_tzrfrq].val[0] = 1400.0;
	}
      if (strlen(psr[0].tzrsite)<1)
	{
	  printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
	  strcpy(psr[0].tzrsite,sitename);
	}


      polyco(psr,npsr,mjd1, mjd2,nspan,ncoeff,maxha,sitename,freq,coeff,1);
    }
    else if (flagPolyco==2)
    {
      long double seg_length;
      int ntimecoeff, nfreqcoeff;
      char sitename[64];
      long double mjd_start, mjd_end;
      long double freq_start, freq_end;

      if (sscanf(polyco_args, "%s %Lf %Lf %Lf %Lf %d %d %Lf", sitename,
		 &mjd_start, &mjd_end, &freq_start, &freq_end,
		 &ntimecoeff, &nfreqcoeff, &seg_length)!=8)
      {
	fprintf(stderr, "Error parsing -pred arguments! See tempo2 -h.\n");
	printf("Have: %s\n",polyco_args);
	exit(1);
      }
      /* Actually want seg_length in days */
      seg_length /= SECDAY;

      if (ntimecoeff%2)
      {
	ntimecoeff++;
	printf("Adjusting number of coefficients (time axis) to %d (must be even)\n", ntimecoeff);
      }
      if (nfreqcoeff%2)
      {
	nfreqcoeff++;
	printf("Adjusting number of coefficients (frequency axis) to %d (must be even)\n", nfreqcoeff);
      }
      if (psr[0].param[param_tzrmjd].paramSet[0]==0)
	{
	  printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)psr[0].param[param_pepoch].val[0]);
	  psr[0].param[param_tzrmjd].paramSet[0]=1;
	  psr[0].param[param_tzrmjd].val[0] = psr[0].param[param_pepoch].val[0];
	}
      if (psr[0].param[param_tzrfrq].paramSet[0]==0)
	{
	  printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)1400.0);
	  psr[0].param[param_tzrfrq].paramSet[0]=1;
	  psr[0].param[param_tzrfrq].val[0] = 1400.0;
	}
      if (strlen(psr[0].tzrsite)<1)
	{
	  printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
	  strcpy(psr[0].tzrsite,sitename);
	}

      ChebyModelSet cms;
      ChebyModelSet_Construct(&cms, psr, sitename, mjd_start, mjd_end,
			      seg_length, seg_length*0.1, 
			      freq_start, freq_end, ntimecoeff, nfreqcoeff);
      FILE *f = fopen("t2pred.dat", "w");
      if (!f)
      {
	fprintf(stderr, "Could not open t2pred.dat for writing!\n");
	exit(1);
      }
      ChebyModelSet_Write(&cms, f);
      fclose(f);
      long double rms, mav;
      ChebyModelSet_Test(&cms, psr, ntimecoeff*5*cms.nsegments, 
			 nfreqcoeff*5*cms.nsegments, &rms, &mav);
      printf("Predictive model constructed and written to t2pred.dat.\n");
      printf("RMS error = %.3Lg s MAV= %.3Lg s\n", 
	     rms/psr[0].param[param_f].val[0], mav/psr[0].param[param_f].val[0]);
      ChebyModelSet_Destroy(&cms);
    }

    return 0;
  }
  if (debugFlag==1)
    {
      logdbg("Number of iterations = %d",psr[0].nits);
      logdbg("Maximum number of parameters = %d",MAX_PARAMS);
      logdbg("Number of pulsars = %d",npsr);
    }

	
	char root[100]; 
	int numTempo2its;
	int doLinearFit;
	int doMax;
	int incEFAC;
	int incEQUAD;
	int incRED;
	int incDM;
	int doTimeMargin;
	int doJumpMargin;
	double FitSig;
	int customPriors;
	int Reddims=0;
	int DMdims=0;
	double *EFACPrior;
	double *EQUADPrior;
	double *AlphaPrior;
	double *AmpPrior;
	double *DMAlphaPrior;
	double *DMAmpPrior;
	int numCoeff;
	double *CoeffPrior;


	char *Type = new char[100];
	EFACPrior=new double[2];
	EQUADPrior=new double[2];
	AlphaPrior=new double[2];
	AmpPrior=new double[2];
	DMAlphaPrior=new double[2];
	DMAmpPrior=new double[2];
	CoeffPrior=new double[2];



	setupparams(Type, numTempo2its, doLinearFit, doMax, incEFAC, incEQUAD, incRED, incDM, doTimeMargin, doJumpMargin, FitSig, customPriors, EFACPrior, EQUADPrior, AlphaPrior, AmpPrior, DMAlphaPrior, DMAmpPrior, numCoeff, CoeffPrior); 



	
  formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
  logdbg("calling formResiduals");
  formResiduals(psr,npsr,1);       /* Form residuals */

  for (it=0;it<numTempo2its;it++) /* Why pulsar 0 should select the iterations? */
    {
      if (it>0) /* Copy post-fit values to pre-fit values */
	{
	  for (i=0;i<MAX_PARAMS;i++)
	    {
	      for (p=0;p<npsr;p++)
		{
		  for (k=0;k<psr[p].param[i].aSize;k++)
		    {
		      psr[p].param[i].prefit[k] = psr[p].param[i].val[k];
		      psr[p].param[i].prefitErr[k] = psr[p].param[i].err[k];
		    }
		}
	    }
	}
      //      long seed = TKsetSeed();
      for (iteration=0;iteration<2;iteration++) /* Do pre- and post- fit analysis */
	{
	  logdbg("iteration %d",iteration);
	  logdbg("calling formBatsAll");
	  //	  printf("Calling formBats\n");
	  formBatsAll(psr,npsr);                /* Form Barycentric arrival times */
	  logdbg("calling formResiduals");
	  formResiduals(psr,npsr,1);       /* Form residuals */


	  if (listparms==1 && iteration==0)displayParameters(13,timFile,parFile,psr,npsr); /* List out all the parameters */  
	  if (iteration==0)          /* Only fit to pre-fit residuals */
	    {
	      logdbg("calling doFit");

	      if (strcmp(dcmFile,"NULL")==0 && strcmp(covarFuncFile,"NULL")==0)
		doFit(psr,npsr,writeModel); /* Fit to the residuals to obtain updated parameters */
	      else
		doFitDCM(psr,dcmFile,covarFuncFile,npsr,writeModel);
	      printf("Complete return\n");
	      /* doFitGlobal(psr,npsr,&globalParameter,nGlobal,writeModel);*/ /* Fit to the residuals to obtain updated parameters  */
	      logdbg("completed doFit");
	    }
	  if (iteration==1 || onlypre==1)
	    {
	      if (strlen(outputSO)==0){
	      //printf("CHI SQ IS: %g \n",psr->fitChisq);
		textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,newparname); /* Output results to the screen */
		//printf("CHI SQ IS: %g %g %g \n",psr->fitChisq,psr[0].offset,psr[0].offset_e);
		}
	      else  /* Use a plug in for the output */
		{
		  char *(*entry)(int,char **,pulsar *,int);
		  void * module;
		  for (int iplug=0; iplug < tempo2_plug_path_len; iplug++){
			  sprintf(str,"%s/%s_%s_plug.t2",tempo2_plug_path[iplug],
					  outputSO,tempo2MachineType);
			  printf("Looking for %s\n",str);
			  module = dlopen(str, RTLD_NOW); 
			  if(module==NULL){	  
				  printf("dlerror() = %s\n",dlerror());
			  } else break;
		  }
		  if(!module)  {
		    fprintf(stderr, "[error]: dlopen() failed while resolving symbols.\n" );
		    return -1;
		  }
		  /*
		   * Check that the plugin is compiled against the same version of tempo2.h
		   */
		  char ** pv  = (char**)dlsym(module, "plugVersionCheck");
		  if(pv!=NULL){
			  // there is a version check for this plugin
			  if(strcmp(TEMPO2_h_VER,*pv)){
				  fprintf(stderr, "[error]: Plugin version mismatch\n");
				  fprintf(stderr, " '%s' != '%s'\n",TEMPO2_h_VER,*pv);
				  fprintf(stderr, " Please recompile plugin against same tempo2 version!\n");
				  dlclose(module);
				  return -1;
			  }
		  }


		  entry = (char*(*)(int,char **,pulsar *,int))dlsym(module, "tempoOutput");
		  if( entry == NULL ) {
		    dlclose(module);
		    fprintf(stderr, "[error]: dlerror() failed while  retrieving address.\n" ); 
		    return -1;
		  }
		  entry(argc,argv,psr,npsr);
		}
	    }
	  psr[0].noWarnings=2;
	  if (onlypre==1) iteration=2;
	  /*	  textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,"new.par");*/ /* Output results to the screen */
	  /*	  printf("Next iteration\n");*/
	}
    }

	if(incRED==0)Reddims=0;
	if(incRED==1)Reddims=2;
	if(incRED==2)Reddims=numCoeff;
 	if(incRED==3)Reddims=2;	
        if(incRED==4)Reddims=3*numCoeff;
        if(incRED==5)Reddims=2*numCoeff+2;
	if(incDM==0)DMdims=0;
	if(incDM==1)DMdims=2;
	if(incDM==2)DMdims=numCoeff;
	if(incDM==3)DMdims=2;
        if(incDM==4)DMdims=3*numCoeff;
        if(incDM==5)DMdims=2*numCoeff+2;

	if(incRED == 1 && incDM != 1 || incRED != 1 && incDM == 1){
		printf("Different methods for DM and red noise not currently supported, please use the same option for both");
		return 0;
	}
	
	
	std::string pulsarname=psr[0].name;
	std::string longname=Type+pulsarname+"-";

	if(longname.size() >= 100){printf("Root Name is too long, needs to be less than 100 characters, currently %i .\n",longname.size());return 0;}

	for(int r=0;r<=longname.size();r++){root[r]=longname[r];}

	std::ofstream getdistparamnames;
	std::string gdpnfname = longname+".paramnames";
	getdistparamnames.open(gdpnfname.c_str());




	
	printf("\n\n\n\n*****************************************************\n");
	printf("Starting TempoNest\n");
	printf("*****************************************************\n\n\n\n");
	printf("Details of the fit:\n");
	printf("file root set to %s \n",root);

	int systemcount=0;
	int *numFlags=new int[psr[0].nobs];
	if(incEFAC == 0 || incEFAC == 1){
		if(incEFAC == 0){printf("Not Including EFAC\n");systemcount=0;}
		if(incEFAC == 1){printf("Including One EFAC for all observations\n");systemcount=1;}
		
			std::string whiteflagname=pulsarname+"-whiteflags.txt";
			std::ofstream whiteflagoutput;
			whiteflagoutput.open(whiteflagname.c_str());
			whiteflagoutput << 1;
			whiteflagoutput <<  "\n";
		
	
		
		for(int o=0;o<psr[0].nobs;o++){
			numFlags[o]=0;
			whiteflagoutput << numFlags[o];
			whiteflagoutput <<  "\n";
		}
		
		whiteflagoutput.close();
		
	}
	else{
		if(incEFAC == 2){printf("Including One EFAC for observing system\n");}
		std::vector<std::string>systemnames;
		for(int o=0;o<psr[0].nobs;o++){
			int found=0;
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
				
				if(strcasecmp(psr[0].obsn[o].flagID[f],"-sys")==0){
					
					if(std::find(systemnames.begin(), systemnames.end(), psr[0].obsn[o].flagVal[f]) != systemnames.end()) {
	// 				/* systemnames contains x */
					} else {
	// 	
	// 
	// 				/* systemnames does not contain x */
					printf("Found system %s \n",psr[0].obsn[o].flagVal[f]);
					systemnames.push_back(psr[0].obsn[o].flagVal[f]);
					systemcount++;
					}
					found=1;
				}

			}
			if(found==0){printf("Observation %i is missing the -sys flag, please check before continuing or set incEFAC=1\n",o);return 0;}
		}
		printf("total number of systems: %i \n",systemcount);
		
		std::string whiteflagname=pulsarname+"-whiteflags.txt";
		std::ofstream whiteflagoutput;
		whiteflagoutput.open(whiteflagname.c_str());
		whiteflagoutput << systemcount;
		whiteflagoutput <<  "\n";
		
		for(int o=0;o<psr[0].nobs;o++){
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
			
				if(strcasecmp(psr[0].obsn[o].flagID[f],"-sys")==0){
					for (int l=0;l<systemcount;l++){
						if(psr[0].obsn[o].flagVal[f] == systemnames[l]){
							numFlags[o]=l;
							whiteflagoutput << numFlags[o];
							whiteflagoutput <<  "\n";
						}
					}
				}
	
			}
		}
		whiteflagoutput.close();
	}

	if(incEQUAD == 0){printf("Not Including EQUAD\n");}
	if(incEQUAD == 1){printf("Including One EQUAD for all observations\n");}
	if(incRED == 0){printf("Not Including Red Noise\n");}
	if(incRED == 1){printf("Including Red Noise : Power Law Model\n");}
	if(incRED == 2){printf("Including Red Noise : Model Independant - Fitting %i Coefficients\n", numCoeff);}
	if(incRED ==3){printf("Including Red Noise: Power Law Model to %i Coefficients \n", numCoeff);}
	if(incDM == 0){printf("Not Including DM\n");}
	if(incDM == 1){printf("Including DM : Power Law Model\n");}
	if(incDM == 2){printf("Including DM : Model Independant - Fitting %i Coefficients\n", numCoeff);}
 	if(incDM ==3){printf("Including DM: Power Law Model to %i Coefficients \n", numCoeff);}	
	
	int fitcount=0;
	printf("fitting for: Arbitrary Phase \n");
	fitcount++;
	for (int p=0;p<MAX_PARAMS;p++) {
	      for (int k=0;k<psr[0].param[p].aSize;k++){
			if(psr[0].param[p].fitFlag[k] == 1){
				printf("fitting for: %s \n",psr[0].param[p].shortlabel[k]);
				fitcount++;
	    		}
		}
	}
// 	printf("fitting for: %i \n",fitcount);

// 	printf("num jumps: %i \n",psr[0].nJumps);

	numFitJumps=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1)numFitJumps++;
// 		  printf("%i %i %g %g \n",i,psr[0].fitJump[i], psr[0].jumpVal[i],psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree));
	}  
	printf("Found %i jumps to fit \n",numFitJumps);
	printf("total Timing Model params to fit:  %i \n",numFitJumps+fitcount);

	long double **TempoPriors;
	long double *Tempo2Fit = new long double[numFitJumps+fitcount];
	int **TempoFitNums;
	int *TempoJumpNums;
	TempoPriors=new long double*[numFitJumps+fitcount];
	for(int i=0;i<numFitJumps+fitcount;i++){
		TempoPriors[i]=new long double[3];
		for(int j=0; j< 3; j++){
				TempoPriors[i][j]=0;
		}
		
	}
	TempoFitNums=new int*[fitcount];for(int i=0;i<fitcount;i++){TempoFitNums[i]=new int[2];}
	TempoJumpNums=new int[numFitJumps];
	//printf("allocated\n");
	int paramsfitted=0;
	int getdistlabel=1;

	getdistparamnames << getdistlabel;
	getdistparamnames << " ";
	getdistparamnames <<  "Phase\n";
	getdistlabel++;
	TempoPriors[paramsfitted][0]=0;
	TempoPriors[paramsfitted][1]=psr[0].offset_e/sqrt(psr[0].fitChisq/psr[0].fitNfree);
	TempoFitNums[paramsfitted][0]=0;
	TempoFitNums[paramsfitted][1]=0;
	if(doTimeMargin != 0 || doJumpMargin != 0)TempoPriors[paramsfitted][2]=1;
	paramsfitted++;
	//printf("p plus\n");
	for (int p=0;p<MAX_PARAMS;p++) {
	      for (int k=0;k<psr[0].param[p].aSize;k++){
			if(psr[0].param[p].fitFlag[k] == 1){
				getdistparamnames << getdistlabel;
				getdistparamnames << " ";
				getdistparamnames <<  psr[0].param[p].shortlabel[k];
				getdistparamnames << "\n";
				TempoPriors[paramsfitted][0]=psr[0].param[p].prefit[k];
				TempoPriors[paramsfitted][1]=psr[0].param[p].err[k]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
				Tempo2Fit[paramsfitted]=psr[0].param[p].val[k];

				TempoFitNums[paramsfitted][0]=p;
				TempoFitNums[paramsfitted][1]=k;	

				if(strcasecmp(psr[0].param[p].shortlabel[0],"F0")== 0){
					if(doTimeMargin != 0)TempoPriors[paramsfitted][2]=1;
				}
				if(doTimeMargin == 2)TempoPriors[paramsfitted][2]=1;				 


				paramsfitted++;
				getdistlabel++;

	    		}
		}
	}
	//printf("and time\n");
	int jumpsfitted=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1){
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "Jump";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			//printf("gonna read ump %i %s \n",i,psr[0].jumpStr[i]);	
			char str1[100],str2[100],str3[100],str4[100],str5[100];
			int nread=sscanf(psr[0].jumpStr[i],"%s %s %s %s %s",str1,str2,str3,str4,str5);
			double prejump=atof(str3);
			//printf("Pre jump %g \n",prejump);
			
			TempoPriors[paramsfitted][0]=prejump;
			TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
			Tempo2Fit[paramsfitted]=psr[0].jumpVal[i];
			TempoJumpNums[jumpsfitted]=i;
			
			if(doJumpMargin == 1)TempoPriors[paramsfitted][2]=1;
			
			paramsfitted++;
			jumpsfitted++;
			getdistlabel++;
		}
	}  

	//printf("dont this\n");
	for(int i =0;i<systemcount;i++){
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "EFAC";
		getdistparamnames << i+1;
		getdistparamnames << "\n";
		paramsfitted++;
		getdistlabel++;
	}
	for(int i =0;i<incEQUAD;i++){
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "EQUAD";
		getdistparamnames << "\n";
		paramsfitted++;
		getdistlabel++;
	}
	if(incRED==1 || incRED == 3){
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "RedAmp";
		getdistparamnames << "\n";
		paramsfitted++;	
		getdistlabel++;	
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "RedSlope";
		getdistparamnames << "\n";
		paramsfitted++;
		getdistlabel++;
	}
	else if(incRED==2){
		for(int i =0;i<numCoeff;i++){
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "RedC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			paramsfitted++;	
			getdistlabel++;
		}
	}
	if(incDM==1 || incDM == 3){
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "DMAmp";
		getdistparamnames << "\n";
		paramsfitted++;	
		getdistlabel++;	
		getdistparamnames << getdistlabel;
		getdistparamnames << " ";
		getdistparamnames <<  "DMSlope";
		getdistparamnames << "\n";
		paramsfitted++;
		getdistlabel++;
	}
	else if(incDM==2){
		for(int i =0;i<numCoeff;i++){
			getdistparamnames << getdistlabel;
			getdistparamnames << " ";
			getdistparamnames <<  "DMC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			paramsfitted++;	
			getdistlabel++;
		}
	}

	getdistparamnames.close();

	// set the MultiNest sampling parameters
	
// 	return 0;
	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 500;				// number of live points
	double efr = 0.1;				// set the required efficiency

	setupMNparams(IS, mmodal, ceff, nlive, efr);


	double tol = 0.1;				// tol, defines the stopping criteria
	int ndims = numFitJumps+fitcount+systemcount+incEQUAD+Reddims+DMdims;					// dimensionality (no. of free parameters)
	int nPar = ndims;					// total no. of parameters including free & derived parameters
	int nClsPar = 2;				// no. of parameters to do mode separation on
	int updInt = 500;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI set it to F if you want your main program to handle MPI initialization
	double logZero = -DBL_MAX;			// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass
	//printf("Here \n");
	MNStruct *MNS = init_struct(psr,TempoPriors,npsr,numFitJumps,fitcount,systemcount,incEQUAD, numCoeff, TempoFitNums,TempoJumpNums,numFlags, ndims, incRED,incDM, doTimeMargin,doJumpMargin, doLinearFit);
	
	
	context=MNS;

#ifdef HAVE_CULA

	culaStatus status;
    status = culaInitialize();
    store_factorial();
    
#endif /* HAVE_CULA */
   
    double **Dpriors;
    Dpriors = new double*[ndims]; for(int i = 0; i < ndims; i++){Dpriors[i]=new double[2];};

  	
     long double *TNMaxParameters = new long double[ndims];
    //If using custompriors for errors incase T2 doesnt converge, get those values before doing anything else
    if(customPriors == 1){
		setTNPriors(Dpriors, TempoPriors);
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
		
		int pcount=1;
		for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			psr[0].param[((MNStruct *)context)->TempoFitNums[j][0]].val[((MNStruct *)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
			//printf("TP: %i %.10Lg %.10Lg \n",pcount,TempoPriors[pcount][0],TempoPriors[pcount][1]);
			pcount++;
		}

		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
			psr[0].jumpVal[((MNStruct *)context)->TempoJumpNums[j]] =  TempoPriors[pcount][0];
	//		printf("TP: %i %.10Lg %.10Lg \n",pcount,TempoPriors[pcount][0],TempoPriors[pcount][1]);
			pcount++;
		}
	}

	if(doMaxLike==1){
#ifdef HAVE_CULA
		GPUFindMLHypervisor(ndims, context,longname);
#else
		FindMLHypervisor(ndims, context,longname);
#endif
		return 0;
	}
	
	//If wanting to find the max do so now
	if(doMax==1){
		
		NelderMeadOptimum(ndims, TNMaxParameters, context);
		for(int i =0; i< numFitJumps+fitcount; i++){
			TempoPriors[i][0]=TNMaxParameters[i];
		}
		
		int pcount=1;
		for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			psr[0].param[((MNStruct *)context)->TempoFitNums[j][0]].val[((MNStruct *)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
			pcount++;
		}

		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
			psr[0].jumpVal[((MNStruct *)context)->TempoJumpNums[j]] =  TempoPriors[pcount][0];
			pcount++;
		}
	}
	
	//reupdate any of the priors from custom priors that were overwritten by findmax
    if(customPriors == 1){
		setTNPriors(Dpriors, TempoPriors);
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
		
		int pcount=1;
		for(int j=1;j<((MNStruct *)context)->numFitTiming;j++){
			psr[0].param[((MNStruct *)context)->TempoFitNums[j][0]].val[((MNStruct *)context)->TempoFitNums[j][1]] = TempoPriors[pcount][0];
	//		printf("TP2: %i %.10Lg %.10Lg \n",pcount,TempoPriors[pcount][0],TempoPriors[pcount][1]);
			pcount++;
		}

		for(int j=0;j<((MNStruct *)context)->numFitJumps;j++){
			psr[0].jumpVal[((MNStruct *)context)->TempoJumpNums[j]] =  TempoPriors[pcount][0];
			pcount++;
		}
	}

	if(doLinearFit != 1){

		//Combine all the priors into one aray: Dpriors
		int pcount=0;
		printf("Setting up priors %i\n",ndims);
		for(int i =0; i< numFitJumps+fitcount; i++){
			Dpriors[pcount][0]=-FitSig;
			Dpriors[pcount][1]=FitSig;
			pcount++;
		}
		for(int i =0; i< systemcount; i++){
			Dpriors[pcount][0]=EFACPrior[0];
			Dpriors[pcount][1]=EFACPrior[1];
			pcount++;
		}
		for(int i =0; i< incEQUAD; i++){
			Dpriors[pcount][0]=EQUADPrior[0];
			Dpriors[pcount][1]=EQUADPrior[1];
			pcount++;
		}	
		if(incRED==1 || incRED == 3){
			Dpriors[pcount][0]=AmpPrior[0];
			Dpriors[pcount][1]=AmpPrior[1];
			pcount++;
			Dpriors[pcount][0]=AlphaPrior[0];
			Dpriors[pcount][1]=AlphaPrior[1];
			pcount++;
		}	
		else if(incRED==2){	
			for(int i =0; i< numCoeff; i++){
				Dpriors[pcount][0]=CoeffPrior[0];
				Dpriors[pcount][1]=CoeffPrior[1];
				pcount++;
			}
		}
/*                if(incRED==4){

                        for(int i =0;i < 2*numCoeff;i++){
                                Dpriors[pcount][0]=-5*psr[p].rmsPre*pow(10.0,-6);
                                Dpriors[pcount][1]=5*psr[p].rmsPre*pow(10.0,-6);
                                printf("Prior on fourier coefficients: %g\n",psr[p].rmsPre*pow(10.0,-6));
                                pcount++;
                        }
                
		        for(int i =0; i< numCoeff; i++){
                                Dpriors[pcount][0]=CoeffPrior[0];
                                Dpriors[pcount][1]=CoeffPrior[1];
                                pcount++;
                        }
}

                if(incRED==5){

			for(int i =0;i < 2*numCoeff;i++){
				Dpriors[pcount][0]=-5*psr[p].rmsPre*pow(10.0,-6);
	                        Dpriors[pcount][1]=5*psr[p].rmsPre*pow(10.0,-6);
				printf("Prior on fourier coefficients: %g\n",psr[p].rmsPre*pow(10.0,-6));
        	                pcount++;
			}

                        Dpriors[pcount][0]=AmpPrior[0];
                        Dpriors[pcount][1]=AmpPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=AlphaPrior[0];
                        Dpriors[pcount][1]=AlphaPrior[1];
                        pcount++;
                }	 */
		if(incDM==1 || incDM == 3){
			printf("DMP %i\n",pcount);
			Dpriors[pcount][0]=DMAmpPrior[0];
			Dpriors[pcount][1]=DMAmpPrior[1];
			pcount++;
			Dpriors[pcount][0]=DMAlphaPrior[0];
			Dpriors[pcount][1]=DMAlphaPrior[1];
			pcount++;
		}
                else if(incDM==2){
                        for(int i =0; i< numCoeff; i++){
                                Dpriors[pcount][0]=CoeffPrior[0];
                                Dpriors[pcount][1]=CoeffPrior[1];
                                pcount++;
                        }
                }

/*
               if(incDM==4){

                        for(int i =0;i < 2*numCoeff;i++){
                                Dpriors[pcount][0]=-5*psr[p].rmsPre*pow(10.0,-6);
                                Dpriors[pcount][1]=5*psr[p].rmsPre*pow(10.0,-6);
                                printf("Prior on fourier coefficients: %g\n",psr[p].rmsPre*pow(10.0,-6));
                                pcount++;
                        }

                        for(int i =0; i< numCoeff; i++){
                                Dpriors[pcount][0]=CoeffPrior[0];
                                Dpriors[pcount][1]=CoeffPrior[1];
                                pcount++;
                        }
}

                if(incDM=5){

                        for(int i =0;i < 2*numCoeff;i++){
                                Dpriors[pcount][0]=-1000*psr[p].rmsPre*pow(10.0,-6);
                                Dpriors[pcount][1]=1000*psr[p].rmsPre*pow(10.0,-6);
                                printf("Prior on fourier coefficients: %g\n",psr[p].rmsPre*pow(10.0,-6));
                                pcount++;
                        }

                        Dpriors[pcount][0]=DMAmpPrior[0];
                        Dpriors[pcount][1]=DMAmpPrior[1];
                        pcount++;
                        Dpriors[pcount][0]=DMAlphaPrior[0];
                        Dpriors[pcount][1]=DMAlphaPrior[1];
                        pcount++;
                }
*/
		printf("set up priors, pcount: %i \n",pcount);


		
		if(doJumpMargin != 0 || doTimeMargin != 0 ){
	
			int *FitList=new int[ndims];
			for(int i=0;i<ndims;i++){
				FitList[i]=0;
			}
			numToMargin=0;
			if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5)numToMargin +=2;
			for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
				if(TempoPriors[i][2]==1){
					FitList[i]=1;
					printf("marginalising over param: %i \n",i);
					numToMargin++;
				}
	}
	
			int Gsize=psr[0].nobs-numToMargin;
			double **TNDM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNDM[i]=new double[numToMargin];
			}
	
			double **TNGM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNGM[i]=new double[psr[0].nobs];

			}

			getCustomDMatrix(psr, FitList, TNDM, TempoFitNums, TempoJumpNums, Dpriors, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps);
			//Get DMatrix for marginalisation, using T2 values if not custom or max, max if not custom, or else custom position.
			//getMarginDMatrix(psr, fitcount, numFitJumps, numToMargin, TempoFitNums, TempoJumpNums, Dpriors, doJumpMargin, doTimeMargin, TNDM, 0);
			makeGDesign(psr, Gsize, numToMargin, TNGM, TNDM);

			update_MNPriors(MNS,Dpriors, TempoPriors,0);
			update_MNGdata(MNS, Gsize,TNGM);

			context=MNS;

			
			//Finally after doing everything reget custom priors in case overwritten by previous steps.
			if(customPriors == 1){
				printf("Set to use custom priors, updating from setPriors function \n");
				setTNPriors(Dpriors, TempoPriors);
				paramsfitted=0;
			}
	//		printf("Step size: %.20Lg \n", TempoPriors[3][1]);
			update_MNPriors(MNS,Dpriors, TempoPriors,0);
			context=MNS;
			
			
#ifdef HAVE_CULA
	double *GMatrixVec=new double[psr[0].nobs*Gsize];


	for(int g=0;g<Gsize; g++){
		for(int o=0;o<psr[0].nobs; o++){

			GMatrixVec[g*psr[0].nobs + o]=TNGM[o][g];
		}
	}

	copy_gmat_(GMatrixVec, psr[0].nobs*Gsize);


#endif /* HAVE_CULA */

			printf("\nPriors:\n");
			paramsfitted=0;
			printf("Prior on Phase : %.25Lg -> %.25Lg\n",TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
			paramsfitted++;
			for (int p=0;p<MAX_PARAMS;p++) {
			for (int k=0;k<psr[0].param[p].aSize;k++){
					if(psr[0].param[p].fitFlag[k] == 1){
						printf("Prior on %s : %.25Lg -> %.25Lg\n",psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
						paramsfitted++;
		
					}
				}
			}
		
			jumpsfitted=0;
			for(int i=0;i<=psr[0].nJumps;i++){
				if(psr[0].fitJump[i] == 1){
					printf("Prior on Jump %i : %.25Lg -> %.25Lg\n",jumpsfitted+1, TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);

					paramsfitted++;
					jumpsfitted++;
				}
			}  
		
			if(incEFAC>0){
				int EFACnum=1;
				for(int i =0;i<systemcount;i++){
					printf("Prior on EFAC %i : %.5g -> %.5g\n",EFACnum, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;
				}
			}
			
			for(int i =0;i<incEQUAD;i++){
				printf("Prior on EQUAD : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			
			if(incRED==1 || incRED == 3){
				printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;		
				printf("Prior on Red Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			else if(incRED==2){
				int Coeffnum=1;
				for(int i =0;i<numCoeff;i++){
					printf("Prior on Red Noise Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;	
				}
			}
			
			if(incDM==1 || incDM == 3){
				printf("Prior on DM Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;		
				printf("Prior on DM Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			else if(incDM==2){
				int Coeffnum=1;
				for(int i =0;i<numCoeff;i++){
					printf("Prior on DM Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;	
				}
			}
			printf("%i %i\n\n", incRED, incDM);
			ndims=ndims-numToMargin;
			if(incDM == 2 || incDM == 3 || incDM == 4 || incDM == 5)ndims +=2;
			nPar=ndims;
		if(incRED==0 && incDM == 0){
#ifdef HAVE_CULA
			nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginGPULogLike, dumper, context);
#else
			nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
		else if(incRED==1 || incDM==1){
#ifdef HAVE_CULA
			nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginGPULogLike, dumper, context);
#else
			nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
		else if(incRED==2 || incRED == 3 || incDM==2 || incDM == 3){
#ifdef HAVE_CULA
			nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginGPULogLike, dumper, context);
#else
			nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
                else if(incRED==4 || incRED == 5 || incDM==4 || incDM == 5){
#ifdef HAVE_CULA
                        nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginNumGPULogLike, dumper, context);
#else
                        nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginNumLogLike, dumper, context);
#endif /* HAVE_CULA */
                }

	}
	else if(doJumpMargin == 0 && doTimeMargin == 0 ){

		//If not marginalising, only need to reget custom priors to overwrite anything done previously
		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors);
		}

		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		context=MNS;
			printf("\nPriors:\n");
			paramsfitted=0;
			printf("Prior on Phase : %.25Lg -> %.25Lg\n",TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
			paramsfitted++;
			for (int p=0;p<MAX_PARAMS;p++) {
			for (int k=0;k<psr[0].param[p].aSize;k++){
					if(psr[0].param[p].fitFlag[k] == 1){
						printf("Prior on %s : %.25Lg -> %.25Lg\n",psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);
						paramsfitted++;
		
					}
				}
			}
		
			jumpsfitted=0;
			for(int i=0;i<=psr[0].nJumps;i++){
				if(psr[0].fitJump[i] == 1){
					printf("Prior on Jump %i : %.25Lg -> %.25Lg\n",jumpsfitted+1, TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][1]*TempoPriors[paramsfitted][1]);

					paramsfitted++;
					jumpsfitted++;
				}
			}  
		
			if(incEFAC>0){
				int EFACnum=1;
				for(int i =0;i<systemcount;i++){
					printf("Prior on EFAC %i : %.5g -> %.5g\n",EFACnum, Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;
				}
			}
			
			for(int i =0;i<incEQUAD;i++){
				printf("Prior on EQUAD : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			
			if(incRED==1 || incRED==3){
				printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;		
				printf("Prior on Red Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			else if(incRED==2){
				int Coeffnum=1;
				for(int i =0;i<numCoeff;i++){
					printf("Prior on Red Noise Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;	
				}
			}
			
			if(incDM==1 || incDM==3){
				printf("Prior on DM Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;		
				printf("Prior on DM Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			else if(incDM==2){
				int Coeffnum=1;
				for(int i =0;i<numCoeff;i++){
					printf("Prior on DM Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
					paramsfitted++;	
				}
			}
		//	printf("dims are: %i \n\n", ndims);

			if(incRED==0){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteLogLike, dumper, context);
		}
			else if(incRED==1){
#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
			else if(incRED==2 || incRED==3 || incDM==2 || incDM==3){
#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
		}
	}
	else if(doLinearFit == 1){

		//Combine all the priors into one aray: Dpriors
		int pcount=0;
		for(int i =0; i< numFitJumps+fitcount; i++){
			Dpriors[pcount][0]=-FitSig;
			Dpriors[pcount][1]=FitSig;
			pcount++;
		}
		for(int i =0; i< systemcount; i++){
			Dpriors[pcount][0]=EFACPrior[0];
			Dpriors[pcount][1]=EFACPrior[1];
			pcount++;
		}
		for(int i =0; i< incEQUAD; i++){
			Dpriors[pcount][0]=EQUADPrior[0];
			Dpriors[pcount][1]=EQUADPrior[1];
			pcount++;
		}	
		if(incRED==1 || incRED ==3){
			Dpriors[pcount][0]=AmpPrior[0];
			Dpriors[pcount][1]=AmpPrior[1];
			pcount++;
			Dpriors[pcount][0]=AlphaPrior[0];
			Dpriors[pcount][1]=AlphaPrior[1];
			pcount++;
		}	
		else if(incRED==2){	
			for(int i =0; i< numCoeff; i++){
				Dpriors[pcount][0]=CoeffPrior[0];
				Dpriors[pcount][1]=CoeffPrior[1];
				pcount++;
			}
		}	
		
		if(incDM==1 || incDM==3){
			Dpriors[pcount][0]=DMAmpPrior[0];
			Dpriors[pcount][1]=DMAmpPrior[1];
			pcount++;
			Dpriors[pcount][0]=DMAlphaPrior[0];
			Dpriors[pcount][1]=DMAlphaPrior[1];
			pcount++;
		}	
		else if(incDM==2){	
			for(int i =0; i< numCoeff; i++){
				Dpriors[pcount][0]=CoeffPrior[0];
				Dpriors[pcount][1]=CoeffPrior[1];
				pcount++;
			}
		}

		int linearNum=numFitJumps+fitcount;
		double **TNDM=new double*[psr[0].nobs];
		for(int i=0;i<psr[0].nobs;i++){
			TNDM[i]=new double[linearNum];
		}
		
					
		//Update Pulsar position to reflect values in TempoPriors which will be either the Tempo2 fit, the max, or the value set in 			custom priors
		
		paramsfitted=1;
		for (int p=0;p<MAX_PARAMS;p++) {
		for (int k=0;k<psr[0].param[p].aSize;k++){
				if(psr[0].param[p].fitFlag[k] == 1){
					((MNStruct *)context)->pulse->param[p].val[k]=TempoPriors[paramsfitted][0];
					//printf("here: %i %.10Lg %.10Lg\n",paramsfitted,TempoPriors[paramsfitted][0],TempoPriors[paramsfitted][1]);
					paramsfitted++;
	
				}
			}
		}
	

		for(int i=0;i<=psr[0].nJumps;i++){
			if(psr[0].fitJump[i] == 1){
				((MNStruct *)context)->pulse->jumpVal[i]=TempoPriors[paramsfitted][0];
				paramsfitted++;
			}
		} 
		
			
	  formBatsAll(((MNStruct *)context)->pulse,npsr);                /* Form Barycentric arrival times */
	  logdbg("calling formResiduals");
	  formResiduals(((MNStruct *)context)->pulse, npsr,1);       /* Form residuals */
		
		getDMatrix(((MNStruct *)context)->pulse, fitcount, numFitJumps, linearNum, TempoFitNums, TempoJumpNums, Dpriors, 1, 2, TNDM);
		
		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors);
			paramsfitted=0;
		}
		update_MNPriors(MNS,Dpriors, TempoPriors,0);
		getLinearPriors(((MNStruct *)context)->pulse, TNDM, TempoPriors, Dpriors, numFitJumps+fitcount, FitSig);
		//printf("set up priors, pcount: %i \n",pcount);

		if(doJumpMargin != 0 || doTimeMargin != 0 ){

                        int *FitList=new int[ndims];
                        for(int i=0;i<ndims;i++){
                                FitList[i]=0;
                        }
                        numToMargin=0;
                        for(int i=0;i<((MNStruct *)context)->numFitTiming + ((MNStruct *)context)->numFitJumps;i++){
                                if(TempoPriors[i][2]==1){
                                        FitList[i]=1;
                                        printf("marginalising over param: %i \n",i);
                                        numToMargin++;
                                }
                        }

                        int Gsize=psr[0].nobs-numToMargin;
                        double **TNDM=new double*[psr[0].nobs];
                        for(int i=0;i<psr[0].nobs;i++){
                                TNDM[i]=new double[numToMargin];
                        }

                        double **TNGM=new double*[psr[0].nobs];
                        for(int i=0;i<psr[0].nobs;i++){
                                TNGM[i]=new double[psr[0].nobs];

                        }

                        getCustomDMatrix(psr, FitList, TNDM, TempoFitNums, TempoJumpNums, Dpriors, ((MNStruct *)context)->numFitTiming, ((MNStruct *)context)->numFitJumps);
                        //Get DMatrix for marginalisation, using T2 values if not custom or max, max if not custom, or else custom position.
                        //getMarginDMatrix(psr, fitcount, numFitJumps, numToMargin, TempoFitNums, TempoJumpNums, Dpriors, doJumpMargin, doTimeMargin, TNDM, 0);
                        makeGDesign(psr, Gsize, numToMargin, TNGM, TNDM);

                        update_MNPriors(MNS,Dpriors, TempoPriors,2);
			update_MNGDdata(MNS, numFitJumps+fitcount, TNDM, Gsize,TNGM);
	
			context=MNS;
		

#ifdef HAVE_CULA
			double *GMatrixVec=new double[psr[0].nobs*Gsize];


			for(int g=0;g<Gsize; g++){
				for(int o=0;o<psr[0].nobs; o++){

					GMatrixVec[g*psr[0].nobs + o]=TNGM[o][g];


				}
				}

			copy_gmat_(GMatrixVec, psr[0].nobs*Gsize);

#endif


			if(incRED==0){

#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
			else if(incRED==1){
#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
			else if(incRED==2 || incRED==3){
#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginLogLike, dumper, context);
#endif /* HAVE_CULA */
			}

		}
		else if(doJumpMargin == 0 && doTimeMargin == 0 ){


			update_MNPriors(MNS,Dpriors, TempoPriors,2);

			update_MNDdata(MNS, numFitJumps+fitcount,TNDM);

			context=MNS;
	

			if(incRED==0){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteLogLike, dumper, context);
		}
			else if(incRED==1){
#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedLogLike, dumper, context);
#endif /* HAVE_CULA */
		}
			else if(incRED==2 || incRED==3){
#ifdef HAVE_CULA
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedGPULogLike, dumper, context);
#else
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
#endif /* HAVE_CULA */
			}
		}

		convertFromLinear(psr, longname, ndims, context);
	}


	readsummary(psr,longname, ndims,context,  Tempo2Fit,incRED, ndims, doTimeMargin, doJumpMargin,doLinearFit);

// 	printf("num its %i \n",psr[0].nits);
	endClock = clock();
  	printf("Finishing off: time taken = %.2f (s)\n",(endClock-startClock)/(float)CLOCKS_PER_SEC);
 	 exit(EXIT_SUCCESS);
} 





// redwards function to force linkage with library functions used by
// plugins
void
thwart_annoying_dynamic_library_stuff(int never_call_me, float or_sink)
{
  ChebyModel *cm;
  T2Predictor *t2p;
  ChebyModel_Init(cm, 0, 0);
  T2Predictor_GetPhase(t2p, 0, 0);
}
