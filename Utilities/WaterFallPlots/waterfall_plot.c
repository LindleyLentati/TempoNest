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
  

  profmax=-1e20;
  profmin=1e20;
  
  resmax=-1e20;
  resmin=1e20;

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
//	  printf("%d %.3f %.3f\n", dbin,dmodel, offset);
	  //offset=0;
	  prof[dbin+nbin*iprof]=dprof-offset;
	  model[dbin+nbin*iprof]=dmodel-offset;
	  resid[dbin+nbin*iprof]=dprof-dmodel;
	  maxmin(model[dbin+nbin*iprof], &mmax, &mmin);
	 
	}
      fclose(proffile);
      // rescale profile to all have same max
      for(ibin=0;ibin<nbin;ibin++)
	{
	  model[ibin+nbin*iprof]/=mmax;
	  prof[ibin+nbin*iprof]/=mmax;
	  resid[ibin+nbin*iprof]/=mmax;
	  
	  //maxmin(model[ibin+nbin*iprof], &profmax, &profmin);
	  maxmin(prof[ibin+nbin*iprof], &profmax, &profmin);
	  maxmin(resid[ibin+nbin*iprof], &resmax, &resmin);
	  
	}


    }
  fclose(listfile);

  
  
  float mjdmax, mjdmin;
  mjdmin=mjd[0];
  mjdmax=mjd[nprof-1];
  
  fprintf(stderr, "%.3e %.3e\n", mjdmax, mjdmin);

  //headspace(&mjdmax, &mjdmin, 0.1);
  
  
  float tr[6];
  

  

  cpgopen(outname);
  cpgpap(0.,1.);
   cpgsvp(0.2, 0.8, 0.5, 0.8);
  cpgsch(1.);
  cpgslw(2);

  cpgswin(mjdmin, mjdmax, pmin, pmax);
  

  
  for(iprof=0;iprof<nprof;iprof++)
    {
      
      if(iprof==0)
	{
	  tr[0]=0.25*(3.*mjdmin-3+2*mjd[iprof]-mjd[iprof+1]);
	  tr[1]=0.5*(mjd[iprof+1]-mjdmin+1);  
      //	tr[0]=mjd[iprof-1];
      //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
      //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
      //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      else if(iprof == nprof-1)
	{
	  tr[0]=0.25*(3.*mjd[iprof-1]+2*mjd[iprof]-mjdmax);
	  tr[1]=0.5*(mjdmax+1-mjd[iprof-1]+1);  
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

	  tr[0]=0.25*(3.*mjd[iprof-1]+2*mjd[iprof]-mjd[iprof+1]-0.03);
	  tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]+0.01);  
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
      cpggray(&prof[nbin*iprof], 1, nbin, 1,1, 1,nbin,profmax, profmin, &tr[0]);


    }
  cpgbox("BCTS", 0, 0, "BCNTS", 0,0);
 
  cpgsvp(0.2, 0.8, 0.2, 0.5);
  cpgsch(1.);
  cpgslw(2);
  cpgswin(mjdmin, mjdmax, pmin, pmax);

  for(iprof=0;iprof<nprof;iprof++)
    {

        if(iprof==0)
	{
	  tr[0]=0.25*(3.*mjdmin-3+2*mjd[iprof]-mjd[iprof+1]);
	  tr[1]=0.5*(mjd[iprof+1]-mjdmin+1);  
      //	tr[0]=mjd[iprof-1];
      //tr[1]=(mjd[iprof+1]-mjd[iprof-1]);
      //tr[0]=0.5*(mjd[iprof-1]+mjd[iprof]);
      //tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]);
	  tr[2]=0;
	  tr[3]=-0.5;
	  tr[4]=0;
	  tr[5]=1./nbin;
	}
      else if(iprof == nprof-1)
	{
	  tr[0]=0.25*(3.*mjd[iprof-1]+2*mjd[iprof]-mjdmax);
	  tr[1]=0.5*(mjdmax+1-mjd[iprof-1]);  
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

	  tr[0]=0.25*(3.*mjd[iprof-1]+2*mjd[iprof]-mjd[iprof+1]-0.03);
	  tr[1]=0.5*(mjd[iprof+1]-mjd[iprof-1]+0.01);  
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
	cpggray(&resid[nbin*iprof], 1, nbin, 1,1, 1,nbin,resmax, resmin, &tr[0]);


    }
  cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
  
  cpgsch(1.5);
  cpgsvp(0.2, 0.8, 0.2, 0.8);
  cpglab("T (yr)", "\\gD\\gF (cy)", "");

  cpgend();

  for(iprof=0;iprof<nprof;iprof++)
    {
      for(ibin=0;ibin<nbin;ibin++)
  	{
  	//  fprintf(stdout, "%d %d %.3e %.3e %.3e\n", iprof, ibin,   prof[ibin+nbin*iprof],  model[ibin+nbin*iprof], resid[ibin+nbin*iprof]);
  	}
   }
  

  free(mjd);
  free(prof);
  free(model);
  free(resid);


  return;
}

int main(int narg, char **argv)
{

  // input file, outputfile, phase range
  
  waterfall_plot("llist.dat", "lwaterfall.eps/VCPS", -0.5,0.5);
  

  return 1;
}

