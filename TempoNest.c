#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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

/* TempoNest.c                                         */
/* --------                                         */
/*                                                  */
/* TEMPO2.C is a clone of TEMPO, but is written     */
/* in C/C++ and contains full descriptions of the   */
/* precision available from each routine.           */
/* Each variable is documented and can be           */
/* tabulated.                                       */
/*                                                  */
/* V1.0 G. Hobbs, R. Edwards                        */

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

#include "/usr/include/gsl/gsl_sf_gamma.h"
#include "dpotri.h"
#include "dpotrf.h"

void ephemeris_routines(pulsar *psr,int npsr);
void clock_corrections(pulsar *psr,int npsr);
void extra_delays(pulsar *psr,int npsr);


void fastephemeris_routines(pulsar *psr,int npsr)
{ 
	vectorPulsar(psr,npsr);   /* 1. Form a vector pointing at the pulsar */
	readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
	tt2tb(psr,npsr);          /* Observatory/time-dependent part of TT-TB */
	readEphemeris(psr,npsr,0);  /* Re-evaluate ephemeris with correct TB */ 

}

void fastformBatsAll(pulsar *psr,int npsr)
{
  	clock_corrections(psr,npsr);          /* Clock corrections  ... */  
 	fastephemeris_routines(psr,npsr);         /* Ephemeris routines ... */
  	extra_delays(psr,npsr);               /* Other time delays  ... */
	formBats(psr,npsr);    /* Form Barycentric arrival times */
	secularMotion(psr,npsr); 

}



MNStruct* insert_data(pulsar *pulseval,	double ** DPriorsval, long double **priorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval,int numFitEFACval, int numFitEQUADval, int numFitCoeffval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval)
{
    	MNStruct* MN = (MNStruct*)malloc(sizeof(MNStruct));

	MN->pulse=pulseval;
	MN->Dpriors=DPriorsval;
	MN->LDpriors=priorsval;
	MN->numberpulsars=numberpulsarsval;
	MN->numFitJumps=numFitJumpsval;
	MN->numFitTiming=numFitTimingval;
	MN->numFitEFAC=numFitEFACval;
	MN->numFitEQUAD=numFitEQUADval;
	MN->numFitRedCoeff=numFitCoeffval;
	MN->TempoFitNums=TempoFitNumsval;
	MN->TempoJumpNums=TempoJumpNumsval;
	MN->sysFlags=sysFlagsval;

	return MN;
}

MNStruct* insert_Gdata(pulsar *pulseval,double ** DPriorsval, long double **priorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval,int numFitEFACval, int numFitEQUADval, int numFitCoeffval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval, int Gsizeval, double **GMatrixval)
{

    	MNStruct* MN = (MNStruct*)malloc(sizeof(MNStruct));

	MN->pulse=pulseval;
	MN->Dpriors=DPriorsval;
	MN->LDpriors=priorsval;
	MN->numberpulsars=numberpulsarsval;
	MN->numFitJumps=numFitJumpsval;
	MN->numFitTiming=numFitTimingval;

	MN->numFitEFAC=numFitEFACval;
	MN->numFitEQUAD=numFitEQUADval;
	MN->numFitRedCoeff=numFitCoeffval;
	MN->TempoFitNums=TempoFitNumsval;
	MN->TempoJumpNums=TempoJumpNumsval;
	MN->sysFlags=sysFlagsval;
	MN->Gsize=Gsizeval;
	MN->GMatrix=GMatrixval;

	return MN;
}

MNStruct* insert_Ddata(pulsar *pulseval,double ** DPriorsval, long double **priorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval,int numFitEFACval, int numFitEQUADval, int numFitCoeffval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval, int Dsizeval, double **DMatrixval)
{

    	MNStruct* MN = (MNStruct*)malloc(sizeof(MNStruct));

	MN->pulse=pulseval;
	MN->Dpriors=DPriorsval;
	MN->LDpriors=priorsval;
	MN->numberpulsars=numberpulsarsval;
	MN->numFitJumps=numFitJumpsval;
	MN->numFitTiming=numFitTimingval;

	MN->numFitEFAC=numFitEFACval;
	MN->numFitEQUAD=numFitEQUADval;
	MN->numFitRedCoeff=numFitCoeffval;
	MN->TempoFitNums=TempoFitNumsval;
	MN->TempoJumpNums=TempoJumpNumsval;
	MN->sysFlags=sysFlagsval;
	MN->Dsize=Dsizeval;
	MN->DMatrix=DMatrixval;

	return MN;
}

MNStruct* insert_GDdata(pulsar *pulseval,double ** DPriorsval, long double **priorsval, int numberpulsarsval,int numFitJumpsval,int numFitTimingval,int numFitEFACval, int numFitEQUADval, int numFitCoeffval, int **TempoFitNumsval,int *TempoJumpNumsval, int *sysFlagsval, int Dsizeval, double **DMatrixval,int Gsizeval, double **GMatrixval)
{

    	MNStruct* MN = (MNStruct*)malloc(sizeof(MNStruct));

	MN->pulse=pulseval;
	MN->Dpriors=DPriorsval;
	MN->LDpriors=priorsval;
	MN->numberpulsars=numberpulsarsval;
	MN->numFitJumps=numFitJumpsval;
	MN->numFitTiming=numFitTimingval;

	MN->numFitEFAC=numFitEFACval;
	MN->numFitEQUAD=numFitEQUADval;
	MN->numFitRedCoeff=numFitCoeffval;
	MN->TempoFitNums=TempoFitNumsval;
	MN->TempoJumpNums=TempoJumpNumsval;
	MN->sysFlags=sysFlagsval;
	MN->Dsize=Dsizeval;
	MN->DMatrix=DMatrixval;
	MN->Gsize=Gsizeval;
	MN->GMatrix=GMatrixval;

	return MN;
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

  for (i=1;i<argc;i++)
    {
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
  for (it=0;it<10;it++) /* Why pulsar 0 should select the iterations? */
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
	      if (strlen(outputSO)==0)
		textOutput(psr,npsr,globalParameter,nGlobal,outRes,newpar,newparname); /* Output results to the screen */
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

	char root[100]; 
	int doLinearFit;
	int incEFAC;
	int incEQUAD;
	int incRED;
	int doTimeMargin;
	int doJumpMargin;
	int FitSig;
	int customPriors;
	int Reddims=0;
	double *EFACPrior;
	double *EQUADPrior;
	double *AlphaPrior;
	double *AmpPrior;
	int numCoeff;
	double *CoeffPrior;


	char *Type = new char[100];
	EFACPrior=new double[2];
	EQUADPrior=new double[2];
	AlphaPrior=new double[2];
	AmpPrior=new double[2];
	CoeffPrior=new double[2];



	setupparams(Type, doLinearFit, incEFAC, incEQUAD, incRED, doTimeMargin, doJumpMargin, FitSig, customPriors, EFACPrior, EQUADPrior, AlphaPrior, AmpPrior, numCoeff, CoeffPrior); 

	if(incRED==0)Reddims=0;
	if(incRED==1)Reddims=2;
	if(incRED==2)Reddims=numCoeff;
	std::string pulsarname=psr[0].name;
	std::string longname=Type+pulsarname+"-";
// 	printf("string size %i \n",longname.size());
	if(longname.size() >= 100){printf("Root Name is too long, needs to be less than 100 characters, currently %i .\n",longname.size());return 0;}
	
// 	root=longname;
// 	printf("%i %i %i %i %i \n",incEFAC,incEQUAD,incRED,FitSig,numCoeff);
// 	printf(" %g %g \n",EQUADPrior[0],EQUADPrior[1]);
// 	printf("long name %s \n",longname);

	for(int r=0;r<=longname.size();r++){root[r]=longname[r];}
// 	printf("long name %s \n",root);
// 	return 0;
	//Open file to use for getdist paramnames
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
		for(int o=0;o<psr[0].nobs;o++){
			numFlags[o]=0;
		}
		
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
		
		for(int o=0;o<psr[0].nobs;o++){
			for (int f=0;f<psr[0].obsn[o].nFlags;f++){
			
				if(strcasecmp(psr[0].obsn[o].flagID[f],"-sys")==0){
					for (int l=0;l<systemcount;l++){
						if(psr[0].obsn[o].flagVal[f] == systemnames[l]){
							numFlags[o]=l;
						}
					}
				}
	
			}
		}
	}
	if(incEQUAD == 0){printf("Not Including EQUAD\n");}
	if(incEQUAD == 1){printf("Including One EQUAD for all observations\n");}
	if(incRED == 0){printf("Not Including Red Noise\n");}
	if(incRED == 1){printf("Including Red Noise : Power Law Model\n");}
	if(incRED == 2){printf("Including Red Noise : Model Independant - Fitting %i Coefficients\n", numCoeff);}
		
	int fitcount=0;
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
	TempoPriors=new long double*[numFitJumps+fitcount];for(int i=0;i<numFitJumps+fitcount;i++){TempoPriors[i]=new long double[2];}
	TempoFitNums=new int*[fitcount];for(int i=0;i<fitcount;i++){TempoFitNums[i]=new int[2];}
	TempoJumpNums=new int[numFitJumps];

	int paramsfitted=0;
	for (int p=0;p<MAX_PARAMS;p++) {
	      for (int k=0;k<psr[0].param[p].aSize;k++){
			if(psr[0].param[p].fitFlag[k] == 1){
				getdistparamnames << paramsfitted+1;
				getdistparamnames << " ";
				getdistparamnames <<  psr[0].param[p].label[k];
				getdistparamnames << "\n";
				TempoPriors[paramsfitted][0]=psr[0].param[p].val[k];
				TempoPriors[paramsfitted][1]=psr[0].param[p].err[k]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
				Tempo2Fit[paramsfitted]=psr[0].param[p].val[k];
// 				if(k == 2 || k == 3){TempoPriors[paramsfitted][1] *= 100;}
				TempoFitNums[paramsfitted][0]=p;
				TempoFitNums[paramsfitted][1]=k;	
				//printf("%i Sigma Prior on %s : %.25Lg +/- %.25Lg\n",FitSig,psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0],FitSig*TempoPriors[paramsfitted][1]);
				paramsfitted++;

	    		}
		}
	}

	int jumpsfitted=0;
	for(int i=0;i<=psr[0].nJumps;i++){
		if(psr[0].fitJump[i] == 1){
			getdistparamnames << paramsfitted+1;
			getdistparamnames << " ";
			getdistparamnames <<  "Jump";
			getdistparamnames << i+1;
			getdistparamnames << "\n";
			TempoPriors[paramsfitted][0]=psr[0].jumpVal[i];
			TempoPriors[paramsfitted][1]=psr[0].jumpValErr[i]/sqrt(psr[0].fitChisq/psr[0].fitNfree);
			TempoJumpNums[jumpsfitted]=i;
			paramsfitted++;
			jumpsfitted++;
		}
	}  


	for(int i =0;i<systemcount;i++){
		getdistparamnames << paramsfitted+1;
		getdistparamnames << " ";
		getdistparamnames <<  "EFAC";
		getdistparamnames << i+1;
		getdistparamnames << "\n";
		paramsfitted++;
	}
	for(int i =0;i<incEQUAD;i++){
		getdistparamnames << paramsfitted+1;
		getdistparamnames << " ";
		getdistparamnames <<  "EQUAD";
		getdistparamnames << "\n";
		paramsfitted++;
	}
	if(incRED==1){
		getdistparamnames << paramsfitted+1;
		getdistparamnames << " ";
		getdistparamnames <<  "Amp";
		getdistparamnames << "\n";
		paramsfitted++;		
		getdistparamnames << paramsfitted+1;
		getdistparamnames << " ";
		getdistparamnames <<  "Slope";
		getdistparamnames << "\n";
		paramsfitted++;
	}
	else if(incRED==2){
		for(int i =0;i<numCoeff;i++){
			getdistparamnames << paramsfitted+1;
			getdistparamnames << " ";
			getdistparamnames <<  "RedC";
			getdistparamnames <<  i+1;
			getdistparamnames << "\n";
			paramsfitted++;	
		}
	}


	getdistparamnames.close();

	// set the MultiNest sampling parameters
	
// 	return 0;
	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 50;				// number of live points
	double efr = 0.5;				// set the required efficiency
	double tol = 0.1;				// tol, defines the stopping criteria
	int ndims = numFitJumps+fitcount+systemcount+incEQUAD+Reddims;					// dimensionality (no. of free parameters)
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





	if(doLinearFit != 1){

		//Combine all the priors into one aray: Dpriors
		int pcount=0;
		double **Dpriors=new double*[ndims]; for(int i = 0; i < ndims; i++){Dpriors[i]=new double[2];}
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
		if(incRED==1){
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
		//printf("set up priors, pcount: %i \n",pcount);

		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors);
			paramsfitted=0;
			for (int p=0;p<MAX_PARAMS;p++) {
				for (int k=0;k<psr[0].param[p].aSize;k++){
					if(psr[0].param[p].fitFlag[k] == 1){
						//printf("New Prior on %s : %.25Lg +/- %i x %.25Lg\n",psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0],FitSig,TempoPriors[paramsfitted][1]);
						paramsfitted++;	
					}
				}
			}
		}

		if(doJumpMargin != 0 || doTimeMargin != 0 ){
	
			int numToMargin=1;
			if(doJumpMargin==1)numToMargin += numFitJumps;
			if(doTimeMargin==1)numToMargin += 2;
			if(doTimeMargin==2)numToMargin += fitcount;
	
			int Gsize=psr[0].nobs-numToMargin;
			double **TNDM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNDM[i]=new double[numToMargin];
			}
	
			double **TNGM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNGM[i]=new double[psr[0].nobs];
	// 			TNGM[i]=new double[psr[0].nobs-numToMargin];
			}
	
	
			getMarginDMatrix(psr, fitcount, numFitJumps, numToMargin, TempoFitNums, TempoJumpNums, Dpriors, doJumpMargin, doTimeMargin, TNDM, 0);
	// 		makeGFromDMatrix(psr, Gsize, numToMargin, TNGM, TNDM);
			makeGDesign(psr, Gsize, numToMargin, TNGM, TNDM);

			MNStruct *MNS = insert_Gdata(psr,Dpriors, TempoPriors,npsr,numFitJumps,fitcount,systemcount,incEQUAD, numCoeff, TempoFitNums,TempoJumpNums,numFlags,Gsize,TNGM);

			context=MNS;


			printf("\nPriors:\n");
			paramsfitted=0;
			for (int p=0;p<MAX_PARAMS;p++) {
			for (int k=0;k<psr[0].param[p].aSize;k++){
					if(psr[0].param[p].fitFlag[k] == 1){
						printf("Prior on %s : %.25Lg -> %.25Lg\n",psr[0].param[p].shortlabel[k], TempoPriors[paramsfitted][0]-Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1]);
						paramsfitted++;
		
					}
				}
			}
		
			jumpsfitted=0;
			for(int i=0;i<=psr[0].nJumps;i++){
				if(psr[0].fitJump[i] == 1){
					printf("Prior on Jump %i : %.25Lg -> %.25Lg\n",jumpsfitted+1, TempoPriors[paramsfitted][0]-Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1],TempoPriors[paramsfitted][0]+Dpriors[paramsfitted][0]*TempoPriors[paramsfitted][1]);

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
			if(incRED==1){
				printf("Prior on Red Noise Log Amplitude : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;		
				printf("Prior on Red Noise Slope : %.5g -> %.5g\n",Dpriors[paramsfitted][0],Dpriors[paramsfitted][1]);
				paramsfitted++;
			}
			else if(incRED==2){
				int Coeffnum=1;
				for(int i =0;i<numCoeff;i++){
					printf("Prior on Red Noise Coefficient %i Log Amplitude : %.5g -> %.5g\n",Coeffnum,Dpriors[paramsfitted][0],Dpriors[paramsfitted][0]);
					paramsfitted++;	
				}
			}
			printf("\n\n");
	
	
			if(incRED==0){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginLogLike, dumper, context);
			}
			else if(incRED==1){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginLogLike, dumper, context);
			}
			else if(incRED==2){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedMarginLogLike, dumper, context);
			}
		}
		else if(doJumpMargin == 0 && doTimeMargin == 0 ){
	
			MNStruct *MNS = insert_data(psr,Dpriors, TempoPriors,npsr,numFitJumps,fitcount,systemcount,incEQUAD, numCoeff, TempoFitNums,TempoJumpNums,numFlags);
	
			context=MNS;
	
			if(incRED==0){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteLogLike, dumper, context);
			}
			else if(incRED==1){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedLogLike, dumper, context);
			}
			else if(incRED==2){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
			}
		}
	}
	else if(doLinearFit == 1){

		ndims++;
		nPar++;
		//Combine all the priors into one aray: Dpriors
		int pcount=0;
		double **Dpriors=new double*[ndims]; for(int i = 0; i < ndims; i++){Dpriors[i]=new double[2];}
		for(int i =0; i< numFitJumps+fitcount+1; i++){
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
		if(incRED==1){
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

		int linearNum=1+numFitJumps+fitcount;
		double **TNDM=new double*[psr[0].nobs];
		for(int i=0;i<psr[0].nobs;i++){
			TNDM[i]=new double[linearNum];
		}
		
		getDMatrix(psr, fitcount, numFitJumps, linearNum, TempoFitNums, TempoJumpNums, Dpriors, 1, 2, TNDM);

		getLinearPriors(psr, TNDM, TempoPriors, Dpriors, numFitJumps+fitcount+1, FitSig);
		printf("set up priors, pcount: %i \n",pcount);

		if(customPriors == 1){
			printf("Set to use custom priors, updating from setPriors function \n");
			setTNPriors(Dpriors, TempoPriors);
			paramsfitted=0;
		}

		if(doJumpMargin != 0 || doTimeMargin != 0 ){


	
			int numToMargin=1;
			if(doJumpMargin==1)numToMargin += numFitJumps;
			if(doTimeMargin==1)numToMargin += 2;
			if(doTimeMargin==2)numToMargin += fitcount;
	
			int Gsize=psr[0].nobs-numToMargin;
			double **DMargin=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				DMargin[i]=new double[numToMargin];
			}
	
			double **TNGM=new double*[psr[0].nobs];
			for(int i=0;i<psr[0].nobs;i++){
				TNGM[i]=new double[psr[0].nobs];
	// 			TNGM[i]=new double[psr[0].nobs-numToMargin];
			}
	
	
			getMarginDMatrix(psr, fitcount, numFitJumps, numToMargin, TempoFitNums, TempoJumpNums, Dpriors, doJumpMargin, doTimeMargin, DMargin, 1);
			makeGDesign(psr, Gsize, numToMargin, TNGM, DMargin);



	
			MNStruct *MNS = insert_GDdata(psr,Dpriors, TempoPriors,npsr,numFitJumps,fitcount,systemcount,incEQUAD, numCoeff, TempoFitNums,TempoJumpNums,numFlags,numFitJumps+fitcount,TNDM, Gsize, TNGM);
			context=MNS;
	
			if(incRED==0){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteMarginLinearLogLike, dumper, context);
			}
			else if(incRED==1){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedMarginLinearLogLike, dumper, context);
			}
			else if(incRED==2){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
			}

		}
		else if(doJumpMargin == 0 && doTimeMargin == 0 ){


			MNStruct *MNS = insert_Ddata(psr,Dpriors, TempoPriors,npsr,numFitJumps,fitcount,systemcount,incEQUAD, numCoeff, TempoFitNums,TempoJumpNums,numFlags,numFitJumps+fitcount,TNDM);
			context=MNS;

			if(incRED==0){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, WhiteLinearLogLike, dumper, context);
			}
			else if(incRED==1){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, vHRedLinearLogLike, dumper, context);
			}
			else if(incRED==2){
				nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LRedLogLike, dumper, context);
			}
		}

		convertFromLinear(psr, longname, ndims, context);
	}


	readsummary(psr,longname, ndims,context,  Tempo2Fit,incRED, ndims, doTimeMargin, doJumpMargin);

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
