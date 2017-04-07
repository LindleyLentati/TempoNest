#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cpgplot.h"

void maxmin(float val, float *max, float *min)
{

  if (val > *max)
    {
      *max =val;
    }
  
  if (val < *min)
    {
      *min = val;
    }

  return;
}

void headspace(float *min, float *max, float fac)
{
  float tmin, tmax;
  tmin = *min;
  tmax = *max;

  *min = tmin - fac*(tmax-tmin);
  *max = tmax + fac*(tmax-tmin);


  return;


}

void waterfall_plot(char *proflist, char *outname, float pmin, float pmax)
{
  
  FILE *listfile, *proffile;

  int iprof,nprof;
  int ibin, nbin;
  int NCHAR_MAX=256;

  listfile=fopen(proflist, "r");
  fscanf(listfile, "%d %d", &nprof, &nbin);

  char fname[160];
  char line[NCHAR_MAX];
  char dumchar;

  float *mjd, *prof, *model, *resid;
  int dbin;
  float dprof, dmodel;
  float dum1, dum2, dum3, dnoise;
  float offset;
  float mmax, mmin;
  float profmax, profmin;
  float resmax, resmin;
  int ep; 

  profmax=-1e20;
  profmin=1e20;
  
  resmax=-1e20;
  resmin=1e20;
  float *noise = (float *) malloc(nprof*sizeof(float));
  mjd= (float*) malloc(nprof*sizeof(float));
  prof = (float*) malloc(nbin*nprof*sizeof(float));
  model = (float*) malloc(nbin*nprof*sizeof(float));
  resid = (float*) malloc(nbin*nprof*sizeof(float));
 
 
  
  for(iprof=0;iprof<nprof;iprof++)
    {
      fscanf(listfile, "%s", fname);
      
      proffile=fopen(fname, "r");
      fgets(line, NCHAR_MAX, proffile);
      // read first line of string which is #
      //sscanf(line,"%c", &dumchar);
      sscanf(line+1,"%f", &mjd[iprof]);
      mjd[iprof]=(mjd[iprof]-53005)/365.25+2004;

      mmax=-1e20;
      mmin=1e20;
      for(ibin=0;ibin<nbin;ibin++)
	{
	  fgets(line, NCHAR_MAX,proffile);
	  // baseline is column 7
	  sscanf(line, "%d %f %f %f", &dbin, &dprof,&dmodel,&dnoise);
	  offset=0;

//	  sscanf(line, "%d %f %f %f %f %f %f", &dbin, &dprof,&dmodel, &dum1, &dum2, &dum3, &offset);
	  //printf("%d %.3f %.3f\n", dbin,dmodel, dnoise);
	  //offset=0;
	  noise[iprof] = dnoise; 
	  prof[dbin+nbin*iprof]=dprof-offset;
	  model[dbin+nbin*iprof]=dmodel-offset;
	  resid[dbin+nbin*iprof]=dprof-dmodel;
	  maxmin(model[dbin+nbin*iprof], &mmax, &mmin);
	 
	}
      fclose(proffile);
      // rescale profile to all have same max
      for(ibin=0;ibin<nbin;ibin++)
	{
//	  printf("noise: %i %i %g \n", iprof, ibin, noise[iprof]);
	  model[ibin+nbin*iprof]/=mmax;
	  prof[ibin+nbin*iprof]/=mmax;
	  resid[ibin+nbin*iprof]/=mmax;
	  if(ibin==0){noise[iprof] /= mmax;}	  
	  //maxmin(model[ibin+nbin*iprof], &profmax, &profmin);
	  maxmin(prof[ibin+nbin*iprof], &profmax, &profmin);
	  maxmin(resid[ibin+nbin*iprof], &resmax, &resmin);
	  
	}


    }
  fclose(listfile);

  
  float mjdmax, mjdmin;
  mjdmin=mjd[0];
  mjdmax=mjd[nprof-1];

////////////////////////////////////Work out Time Averaging////////////////////////////////////////

	double window = 1/365.25;
	int MaxEpochs = floor((mjdmax-mjdmin)/window+0.5)+1;
	printf("max epochs for window %g is %i\n", window, MaxEpochs);

	int *EpochCount = (int *) malloc(MaxEpochs*sizeof(int));
	for(iprof = 0; iprof < MaxEpochs; iprof++){
		EpochCount[iprof] = 0;
	}
 
	for(iprof=0;iprof<nprof;iprof++){
		int mjdbin = floor((mjd[iprof]-mjdmin)/window);
		EpochCount[mjdbin]++;
	}
	int NumEpochs = 0;
	for(ep = 0; ep < MaxEpochs; ep++){
		if(EpochCount[ep] > 0){
			NumEpochs++;
		}
	}
	printf("Using : %i epochs \n", NumEpochs);

	float *EpochAveragedData = (float *) malloc(NumEpochs*nbin*sizeof(float));
	float *EpochAveragedModel = (float *) malloc(NumEpochs*nbin*sizeof(float));
	float *EpochAveragedRes = (float *) malloc(NumEpochs*nbin*sizeof(float));
	float *EpochWeights = (float *) malloc(NumEpochs*sizeof(float));
	float *EpochMJDs = (float *) malloc(NumEpochs*sizeof(float));


	for(iprof=0;iprof<NumEpochs*nbin;iprof++){
		EpochAveragedData[iprof] = 0;
		EpochAveragedModel[iprof] = 0;
		EpochAveragedRes[iprof] = 0;
	}
	for(iprof=0;iprof<NumEpochs; iprof++){
		EpochWeights[iprof] = 0;	
	}

	int realmjdbin = 0;
	int prevbin = 0;
	for(iprof=0;iprof<nprof;iprof++){
		int mjdbin = floor((mjd[iprof]-mjdmin)/window);
		if(mjdbin != prevbin){
			prevbin = mjdbin;
			realmjdbin++;
		}
//		printf("weight: %i %i %i %g \n", iprof, mjdbin, realmjdbin, 1.0/noise[iprof]/noise[iprof]);
		for(ibin=0;ibin<nbin;ibin++){
			//if(model[ibin+nbin*iprof] > 0.9)printf("model: %i %i %g \n", iprof, ibin,  model[ibin+nbin*iprof]);
			EpochAveragedData[ibin+nbin*realmjdbin] += prof[ibin+nbin*iprof]/noise[iprof]/noise[iprof];
			EpochAveragedModel[ibin+nbin*realmjdbin] += model[ibin+nbin*iprof]/noise[iprof]/noise[iprof];
			EpochAveragedRes[ibin+nbin*realmjdbin] += resid[ibin+nbin*iprof]/noise[iprof]/noise[iprof];
			
		}
		EpochWeights[realmjdbin] += 1.0/noise[iprof]/noise[iprof];
		EpochMJDs[realmjdbin] = mjd[iprof];
	}

	for(iprof=0;iprof<NumEpochs; iprof++){
		 for(ibin=0;ibin<nbin;ibin++){

			EpochAveragedData[ibin+nbin*iprof] /= EpochWeights[iprof];
			EpochAveragedModel[ibin+nbin*iprof] /= EpochWeights[iprof];
			EpochAveragedRes[ibin+nbin*iprof] /= EpochWeights[iprof];
			printf("%g %i %g %g %g %g \n",  EpochMJDs[iprof], ibin, EpochAveragedData[ibin+nbin*iprof], EpochAveragedModel[ibin+nbin*iprof] , EpochAveragedRes[ibin+nbin*iprof], 1.0/sqrt(EpochWeights[iprof]));
		}
	}

	
  
  fprintf(stderr, "%.3e %.3e\n", mjdmax, mjdmin);

  //headspace(&mjdmax, &mjdmin, 0.1);
  
  
  float tr[6];
  

  

  cpgopen(outname);
  cpgpap(0.,1.);
   cpgsvp(0.2, 0.8, 0.5, 0.8);
  cpgsch(1.);
  cpgslw(2);

  cpgswin(mjdmin, mjdmax, pmin, pmax);
  

  
  for(iprof=0;iprof<NumEpochs;iprof++)
    {
      
      if(iprof==0)
	{
	  tr[0]=0.25*(3.*mjdmin-3+2*EpochMJDs[iprof]-EpochMJDs[iprof+1]);
	  tr[1]=0.5*(EpochMJDs[iprof+1]-mjdmin+1);  
      //	tr[0]=mjd[iprof-1];
      //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
      //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
      //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      else if(iprof == NumEpochs-1)
	{
	  tr[0]=0.25*(3.*EpochMJDs[iprof-1]+2*EpochMJDs[iprof]-mjdmax);
	  tr[1]=0.5*(mjdmax+1-EpochMJDs[iprof-1]+1);  
	  //	tr[0]=mjd[iprof-1];
	  //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
	  //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
	  //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      else
	{

	  tr[0]=0.25*(3.*EpochMJDs[iprof-1]+2*EpochMJDs[iprof]-EpochMJDs[iprof+1]-0.03);
	  tr[1]=0.5*(EpochMJDs[iprof+1]-EpochMJDs[iprof-1]+0.01);  
	  //	tr[0]=mjd[iprof-1];
	  //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
	  //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
	  //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      // 
      cpggray(&EpochAveragedData[nbin*iprof], 1, nbin, 1,1, 1,nbin,profmax, profmin, &tr[0]);


    }
  cpgbox("BCTS", 0, 0, "BCNTS", 0,0);
 
  cpgsvp(0.2, 0.8, 0.2, 0.5);
  cpgsch(1.);
  cpgslw(2);
  cpgswin(mjdmin, mjdmax, pmin, pmax);

  for(iprof=0;iprof<NumEpochs;iprof++)
    {

        if(iprof==0)
	{
	  tr[0]=0.25*(3.*mjdmin-3+2*EpochMJDs[iprof]-EpochMJDs[iprof+1]);
	  tr[1]=0.5*(EpochMJDs[iprof+1]-mjdmin+1);  
      //	tr[0]=mjd[iprof-1];
      //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
      //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
      //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      else if(iprof == NumEpochs-1)
	{
	  tr[0]=0.25*(3.*EpochMJDs[iprof-1]+2*EpochMJDs[iprof]-mjdmax);
	  tr[1]=0.5*(mjdmax+1-EpochMJDs[iprof-1]);  
	  //	tr[0]=mjd[iprof-1];
	  //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
	  //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
	  //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      else
	{

	  tr[0]=0.25*(3.*EpochMJDs[iprof-1]+2*EpochMJDs[iprof]-EpochMJDs[iprof+1]-0.03);
	  tr[1]=0.5*(EpochMJDs[iprof+1]-EpochMJDs[iprof-1]+0.01);  
	  //	tr[0]=mjd[iprof-1];
	  //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
	  //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
	  //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}

    ;

      // 
	cpggray(&EpochAveragedRes[nbin*iprof], 1, nbin, 1,1, 1,nbin,resmax, resmin, &tr[0]);


    }
  cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
  
  cpgsch(1.5);
  cpgsvp(0.2, 0.8, 0.2, 0.8);
  cpglab("T (yr)", "\\gD\\gF (cy)", "");

  cpgend();

  

  free(mjd);
  free(prof);
  free(model);
  free(resid);


  return;
}

int main(int narg, char **argv)
{

  // input file, outputfile, phase range

	double phasemin=0;
	double phasemax = 0;
	char *inputfile;
	char *outputfile; 

	
	inputfile = argv[0+1];
	outputfile = argv[1+1];
	 sscanf(argv[2+1],"%lf",&phasemin);
	 sscanf(argv[3+1],"%lf",&phasemax);

                
 	printf("making plot with: %s %s %g %g \n", inputfile, outputfile, phasemin, phasemax);
 
  waterfall_plot(inputfile, "lwaterfall.eps/VCPS", phasemin,phasemax);
  

  return 1;
}

