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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "tempo2.h"
#include "constraints.h"
#include "TKfit.h"
#include "TempoNest.h"
#include <algorithm>

//#define TSUN (4.925490947e-6L) (Should be tempo2.h now).

double TNm2(longdouble mf, longdouble sini, longdouble m1);
void TNprintGlitch(pulsar psr);
double TNdglep(pulsar psr,int gn,double fph);
double TNcalcRMS(pulsar *psr,int p);

/* ******************************************** */
/* textOutput                                   */
/* Author:  G. Hobbs (17 May 2003)              */
/* Purpose: Displays the fitted results on the  */
/*          screen.                             */
/* Inputs:  pulsar structure                    */
/*                                              */
/* Outputs:                                     */
/*                                              */
/* Notes:                                       */
/*                                              */
/* Changes:                                     */
/* ******************************************** */

int longturn_hms(long double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  long double sec;

  hh = (int)(turn*24.);
  mm = (int)((turn*24.-hh)*60.);
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (int)((sec*10000. +0.5)/10000);
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

  sprintf(hms," %02d:%02d:%.12Lg",hh,mm,sec);
  return 0;
}

int longturn_dms(long double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  long double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = (int)(trn*360.);
  mm = (int)((trn*360.-dd)*60.);
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (int)((sec*1000. +0.5)/1000);
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
  sprintf(dms,"%c%02d:%02d:%.12Lg",sign,dd,mm,sec);
  return 0;
}



void TNtextOutput(pulsar *psr, int npsr, int newpar, long double *Tempo2Fit, void *context,int incRED, int ndim, std::vector<double> paramlist, double Evidence,int doTimeMargin, int doJumpMargin, int doLinear, std::string longname)
{
  double rms_pre=0.0,rms_post=0.0;
  double mean_pre=0.0,mean_post=0.0,chisqr;
  int i,p,count,k;
  FILE *fout;
	char *fname;


    logdbg("In textOutput");

  for (p=0;p<npsr;p++)
    {
      rms_pre=0.0;
      rms_post=0.0;
      mean_pre=0.0;
      mean_post=0.0;
      count=0;
      TNcalcRMS(psr,p);
      {
	long double centrePos;
	long double closestV,check;
	int closestI=-1;
	centrePos = (psr[p].param[param_start].val[0]+
		     psr[p].param[param_finish].val[0])/2.0L;
	
	
	for (i=0;i<psr[p].nobs;i++)
	  {
	    if (psr[p].obsn[i].deleted==0)
	      {
		check = fabs(psr[p].obsn[i].sat-centrePos);
		if (closestI==-1 || (check < closestV))
		  {
		    closestV = check;
		    closestI = i;
		  }
	      }
	  }
	if (closestI==-1)
	  {
	    printf("WARNING: Cannot calculate TZRMJD\n");
	    psr[p].param[param_tzrmjd].paramSet[0]=0;
	  }
	else
	  {
	    psr[p].param[param_tzrmjd].val[0] = psr[p].obsn[closestI].sat-psr[p].obsn[closestI].residual/86400.0L;
	    psr[p].param[param_tzrmjd].paramSet[0] = 1;
	    psr[p].param[param_tzrfrq].val[0] = psr[p].obsn[closestI].freq;
	    psr[p].param[param_tzrfrq].paramSet[0] = 1;
	    strcpy(psr[p].tzrsite,psr[p].obsn[closestI].telID);
	  }
      }

      
      /*	  wrms_pre = sqrt((sumsq_pre/(double)count-sum_pre*sum_pre/pow(double(count),2))/(sumwt/(double)count)/(sumwt/double(count)))*1e6; */
      
      printf("\n\n");
      printf("TempoNest Results for PSR %s\n",psr[p].name);
      printf("\n\n");


      if (psr[p].fitMode==1)
	{
	  printf("Evidence = %.10g\t",Evidence);
	  /*      printf("wmean = %g\n",sumwt/(double)psr[p].nobs); */
	  /* chisqr = psr[p].fitChisq*pow(psr[p].param[param_f0].val,2)*sumwt/(double)psr[p].nobs/(double)psr[p].fitNfree; */
	  chisqr = psr[p].fitChisq;

	}
	int pcount=1;
	int fitcount=0;
	 if(((MNStruct *)context)->LDpriors[0][2]==0)fitcount++;

      printf("\n\n");
      printf("PARAMETER          Tempo2-Fit              TempoNest-Fit          Uncertainty   Difference   Fit\n");
      printf("---------------------------------------------------------------------------------------------------\n");
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1)
		{
		  /* PARAMETER (name) */
		  if (i == param_raj && psr[p].eclCoord==1)
		    printf("%-15.15s ","ELONG");
		  else if (i == param_decj && psr[p].eclCoord==1)
		    printf("%-15.15s ","ELAT");
 		  else if (i == param_pmra && psr[p].eclCoord==1)
		    printf("%-15.15s ","PMELONG");
 		  else if (i == param_pmdec && psr[p].eclCoord==1)
		    printf("%-15.15s ","PMELAT");
		  else
		    printf("%-15.15s ",psr[p].param[i].label[k]);
		  
		  /* Pre-fit value */
		  if ((i == param_raj || i == param_decj) && psr[p].eclCoord==1)
		    printf("%-25.15g ",(double)Tempo2Fit[pcount]*180.0/M_PI);
      else if (i == param_sini && psr[p].param[param_sini].nLinkTo != 0)
        printf("%-25.15g ",(double)sin(Tempo2Fit[pcount]/180.0*M_PI));
		  else if (i == param_ephver)
		    {
		      if (psr[p].param[i].prefit[k]==2)
			printf("%-25.25s ","TEMPO1");
		      if (psr[p].param[i].prefit[k]==5)
			printf("%-25.25s ","TEMPO2");
		      else
			printf("%-25.25g ",(double)Tempo2Fit[pcount]);
		    }
		  else
		    printf("%-25.15g ",(double)Tempo2Fit[pcount]);
		  
		  /* Post-fit value */
		  if (i==param_pmra) /* Convert from radian/sec to mas/yr */
		    printf("%-25.15g ",(double)psr[p].param[i].val[k]);
		  /* *180.0/M_PI*60.0*60.0* 1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val)); */
		  else if (i==param_pmdec) /* Convert from radian/sec to mas/yr */
		    printf("%-25.15g ",(double)psr[p].param[i].val[k]);
		  /* *180.0/M_PI*60.0*60.0*1000.0*SECDAY*365.25/24.0/3600.0);*/
		  else if ((i == param_raj || i == param_decj) && psr[p].eclCoord==1)
		    printf("%-25.15g ",(double)psr[p].param[i].val[k]*180.0/M_PI);
		  else if (i == param_ephver)
		    {
		      if (psr[p].param[i].val[k]==2)
			printf("%-25.25s ","TEMPO1");
		      if (psr[p].param[i].val[k]==5)
			printf("%-25.25s ","TEMPO2");
		      else
			printf("%-25.25g ",(double)psr[p].param[i].val[k]);
		    }
		  else printf("%-25.15g ",(double)getParameterValue(&psr[p],i,k));
		  
		  if ((i == param_raj || i == param_decj) && psr[p].eclCoord==1)
		    {
		      printf("%-13.5g ", (double)psr[p].param[i].err[k]*180.0/M_PI);
		      printf("%-13.5g ", 180.0/M_PI*((double)psr[p].param[i].val[k]-(double)Tempo2Fit[pcount]));
		    }
		  else
		    {
		      printf("%-13.5g ", (double)psr[p].param[i].err[k]);
		      printf("%-13.5g ", (double)getParameterValue(&psr[p],i,k)-(double)Tempo2Fit[pcount]);
		    }
		  if (doTimeMargin==0){ 
				printf("Y");
			}
		  else if (doTimeMargin==1) {
			if (i == param_f){	
				printf("M");
			}
			else{
				printf("Y");
			}
		  }
		  else if (doTimeMargin==2) {

				printf("M");
			}


		  printf("\n");

		  if (i==param_tzrfrq)
		    {
		      if (strcmp(psr[p].tzrsite,"NULL")!=0) 
			printf("%-15.15s %-25.25s\n","TZRSITE",psr[p].tzrsite);
		    }
		  if (i==param_raj && psr[p].eclCoord==0)
		    {
	                char hmsstr[100];
        	        longturn_hms(psr[p].param[param_raj].val[0]/(2*M_PI), hmsstr);
                	strcpy(psr[p].rajStrPost,hmsstr);

		      printf("%-15.15s %-25.25s %-25.25s %-13.5g %-13.5g\n","RAJ (hms)",psr[p].rajStrPre,psr[p].rajStrPost,
			     (double)psr[p].param[i].err[k]*12.0*60.0*60.0/M_PI,
			     ((double)psr[p].param[i].val[k]-(double)Tempo2Fit[pcount])*12.0*60.0*60.0/M_PI);
		    }
		  else if (i==param_decj && psr[p].eclCoord==0)
		    {
			char hmsstr[100];
        	        longturn_dms(psr[p].param[param_decj].val[0]/(2*M_PI), hmsstr);
	                strcpy(psr[p].decjStrPost,hmsstr);

		      printf("%-15.15s %-25.25s %-25.25s %-13.5g %-13.5g\n","DECJ (dms)",psr[p].decjStrPre,psr[p].decjStrPost,
			     (double)psr[p].param[i].err[k]*180.0*60.0*60.0/M_PI,
			     ((double)psr[p].param[i].val[k]-(double)Tempo2Fit[pcount])*180.0*60.0*60.0/M_PI);
		    }
		
			if(((MNStruct *)context)->LDpriors[pcount][2]==0)fitcount++;
			pcount++;
		}
	    }
	}      
      printf("---------------------------------------------------------------------------------------------------\n");
      if (psr->rescaleErrChisq == 1 && psr->fitMode==1)
      if (psr[p].nconstraints>0){
	printf("\nCONSTRAINTS:\n");
	for (i=0; i < psr[p].nconstraints; i++){
	  printf("%s\n",get_constraint_name(psr[p].constraints[i]).c_str());
	}
	printf("\n");
      }

      /* JUMPS */
      for (i=1;i<=psr[p].nJumps;i++){
	{
	  printf("Jump %d (%s): %.14g %.14g ",i,psr[p].jumpStr[i],psr[p].jumpVal[i],psr[p].jumpValErr[i]);
	  if (psr[p].fitJump[i]==1) {
		if(doJumpMargin==0){
			printf("Y\n");
			fitcount++;
		}
		else if(doJumpMargin==1){
			printf("M\n");
		}
		pcount++;
		
	  }
	  else printf("N\n");
	}
      }

	if(((MNStruct *)context)->incStep > 0){
		printf("%i Step Functions used:\n",((MNStruct *)context)->incStep);
		for(int i =0; i < ((MNStruct *)context)->incStep; i++){

			printf("Step Amp %i: %g +/- %g\n", i+1, paramlist[fitcount],paramlist[fitcount+ndim]);
			fitcount++;
			printf("Step Time %i: %g +/- %g\n", i+1, paramlist[fitcount],paramlist[fitcount+ndim]);
			fitcount++;
		}	
	}
                std::vector<int>systempos;
                std::vector<int>sysflag;
                std::vector<std::string>systemnames;
		
                if(((MNStruct *)context)->numFitEFAC > 1 || ((MNStruct *)context)->numFitEQUAD > 1){
                        for(int o=0;o<((MNStruct *)context)->pulse[0].nobs;o++){
                                int found=0;
                                for (int f=0;f<((MNStruct *)context)->pulse[0].obsn[o].nFlags;f++){

                                        if(strcasecmp(((MNStruct *)context)->pulse[0].obsn[o].flagID[f],((MNStruct *)context)->whiteflag)==0){

                                                if(std::find(systemnames.begin(), systemnames.end(), ((MNStruct *)context)->pulse[0].obsn[o].flagVal[f]) != systemnames.end()) {
                                                } else {

                                                        systempos.push_back(o);
                                                        sysflag.push_back(f);
                                                        systemnames.push_back(((MNStruct *)context)->pulse[0].obsn[o].flagVal[f]);

                                                }
                                                found=1;
                                        }

                                }
                        }
                }

  	
	if(incRED != 0 || ((MNStruct *)context)->incDM !=0 ||((MNStruct *)context)->numFitEFAC > 0 || ((MNStruct *)context)->numFitEQUAD > 0){
                printf("------------------------------------------------------------------------------\n");
                printf("Stochastic Parameters:\n");
                if(((MNStruct *)context)->numFitEFAC == 1){
                        printf("Global EFAC: %g +/- %g\n", paramlist[fitcount],paramlist[fitcount+ndim]);
                        fitcount++;
                }
                else if(((MNStruct *)context)->numFitEFAC > 1){
                        int system=1;
                        for(int i =0;i<((MNStruct *)context)->numFitEFAC; i++){
                                printf("EFAC for %s %s: %g +/- %g\n",((MNStruct *)context)->whiteflag, ((MNStruct *)context)->pulse[0].obsn[systempos[i]].flagVal[sysflag[i]] ,paramlist[fitcount],paramlist[fitcount+ndim]);
                                fitcount++;
                                system++;
                        }
                }

                if(((MNStruct *)context)->numFitEQUAD ==1){
                        printf("Global EQUAD: %g +/- %g\n", paramlist[fitcount],paramlist[fitcount+ndim]);
                        fitcount++;
                }
                else if(((MNStruct *)context)->numFitEQUAD > 1){
                        int system=1;
                        for(int i =0;i<((MNStruct *)context)->numFitEFAC; i++){
                                printf("EQUAD for %s %s: %g +/- %g\n", ((MNStruct *)context)->whiteflag, ((MNStruct *)context)->pulse[0].obsn[systempos[i]].flagVal[sysflag[i]],paramlist[fitcount],paramlist[fitcount+ndim]);
                                fitcount++;
                                system++;
                        }
                }

		if(incRED ==1 || incRED ==3){
			printf("Power Law Red Noise Model:\n");
			printf("Log Amplitude: %g +/- %g\n",paramlist[fitcount],paramlist[fitcount+ndim]);
			fitcount++;
			printf("Spectral Index: %g +/- %g\n",paramlist[fitcount],paramlist[fitcount+ndim]);	
			fitcount++;
		}
		else if(incRED ==2){	
			printf("Model Independant Red Noise, %i Coefficients used:\n",((MNStruct *)context)->numFitRedCoeff);
			int coeff=1;
			for(int i =0; i < ((MNStruct *)context)->numFitRedCoeff; i++){
				printf("Log Amplitude Coefficient %i: %g +/- %g\n", coeff,paramlist[fitcount],paramlist[fitcount+ndim]);	
				fitcount++;
				coeff++;
			}
		}

		if(((MNStruct *)context)->incFloatRed > 0){
			printf("Floating Red Noise, %i Coefficients used:\n",((MNStruct *)context)->incFloatRed);
			for(int i =0; i < ((MNStruct *)context)->incFloatRed; i++){

				printf("Frequency %i: %g +/- %g\n", i+1, paramlist[fitcount],paramlist[fitcount+ndim]);
				fitcount++;
				printf("Log Amplitude %i: %g +/- %g\n", i+1, paramlist[fitcount],paramlist[fitcount+ndim]);
				fitcount++;
		    }
		}


		if(((MNStruct *)context)->incDM ==1 ||((MNStruct *)context)->incDM ==3 ){
			printf("Power Law DM Model:\n");
			printf("Log Amplitude: %g +/- %g\n",paramlist[fitcount],paramlist[fitcount+ndim]);
			fitcount++;
			printf("Spectral Index: %g +/- %g\n",paramlist[fitcount],paramlist[fitcount+ndim]);	
			fitcount++;
		}
		else if(((MNStruct *)context)->incDM==2){	
			printf("Model Independant DM, %i Coefficients used:\n",((MNStruct *)context)->numFitDMCoeff);
			int coeff=1;
			for(int i =0; i < ((MNStruct *)context)->numFitDMCoeff; i++){
				printf("Log Amplitude Coefficient %i: %g +/- %g\n", coeff,paramlist[fitcount],paramlist[fitcount+ndim]);	
				fitcount++;
				coeff++;
			}
		}

		if(((MNStruct *)context)->incFloatDM > 0){
			printf("Floating DM, %i Coefficients used:\n",((MNStruct *)context)->incFloatDM);
			for(int i =0; i < ((MNStruct *)context)->incFloatDM; i++){

				printf("Frequency %i: %g +/- %g\n", i+1, paramlist[fitcount],paramlist[fitcount+ndim]);
				fitcount++;
				printf("Log Amplitude %i: %g +/- %g\n", i+1, paramlist[fitcount],paramlist[fitcount+ndim]);
				fitcount++;
		    }
		}

	}


      /* Whitening */
      if (psr[p].param[param_wave_om].paramSet[0]==1)
	{
	  double perr,pwr;
	  double xval,yval;
	  int j;
	  FILE *fout;
	  printf("      \t%-15.15s %-15.15s %-15.15s %-15.15s %-15.15s\n","Freq","Period","Cosine amp","Sine amp","Power");
	  printf("      \t%-15.15s %-15.15s %-15.15s %-15.15s %-15.15s\n","(yr^-1)","(yr)","(s)","(s)","(s^2)");
	  printf("------------------------------------------------------------------------------\n");
	  for (i=0;i<psr[p].nWhite;i++)
	    {
	      pwr = pow(psr[p].wave_cos[i],2)+pow(psr[p].wave_sine[i],2);
	      perr = sqrt(pow(2*psr[p].wave_cos[i]*psr[p].wave_cos_err[i],2)+pow(2*psr[p].wave_sine[i]*psr[p].wave_sine_err[i],2));
	      printf("WAVE%d\t%-15.5Lg %-15.5Lg %-+10.5g %-+10.5g %-+10.5g %-+10.5g %-+10.5g %-+10.5g\n",
               i+1, // Wave number (counter starting at 1 - i.e. 'i' starts at 0)
		     (i+1)*psr[p].param[param_wave_om].val[0]/2.0/M_PI*365.25,       // Wave frequency (yr^-1) JORIS
               1.0/((i+1)*psr[p].param[param_wave_om].val[0]/2.0/M_PI*365.25), // Wave period (yrs)
               psr[p].wave_cos[i],       // Wave cosine amplitude
               psr[p].wave_cos_err[i],   // Wave cosine amplitude uncertainty
               psr[p].wave_sine[i],      // Wave sine amplitude
               psr[p].wave_sine_err[i],  // Wave sine amplitude uncertainty
               pwr,  // Wave power
               perr); // Wave power uncertainty
	    }
	  printf("------------------------------------------------------------------------------\n");
	}
      if (psr[p].param[param_gwsingle].paramSet[0]==1)
	{
	  printf("GW single source:\n");
	  printf("Omega: %g\n",(double)psr[p].param[param_gwsingle].val[0]);
	  printf("Aplus = %g (%g) %g (%g)\n",psr[p].gwsrc_aplus_r,psr[p].gwsrc_aplus_r_e,
		 psr[p].gwsrc_aplus_i,psr[p].gwsrc_aplus_i_e);
	  printf("Across = %g (%g) %g (%g)\n",psr[p].gwsrc_across_r,psr[p].gwsrc_across_r_e,
		 psr[p].gwsrc_across_i,psr[p].gwsrc_across_i_e);
	}
      if (psr[p].param[param_quad_om].paramSet[0]==1)
	{
	  int j;
	  for (j=0;j<psr[p].nQuad;j++)
	    {
	      printf("QUAD%d %g %g %g %g %g %g %g %g %g %g\n",j+1,(double)(psr[p].param[param_quad_om].val[0])*(j+1),
		     (double)(psr[p].param[param_quad_om].val[0]*(j+1)/2.0/M_PI),(double)psr[p].quad_aplus_r[j],
		     (double)psr[p].quad_aplus_i[j],(double)psr[p].quad_across_r[j],
		     (double)psr[p].quad_across_i[j],(double)psr[p].quad_aplus_r_e[j],(double)psr[p].quad_aplus_i_e[j],
		     (double)psr[p].quad_across_r_e[j],(double)psr[p].quad_across_i_e[j]);
	    }
	}
	  
      if (psr[p].param[param_dmmodel].paramSet[0]==1)
	{
	  printf("\n");
	  printf("Dispersion measure values:\n");
	  char dmfname[MAX_STRLEN];
	  sprintf(dmfname,"%s.dm",psr[p].name);
	  FILE *dmfile=fopen(dmfname,"w");
	  double sum=0;
	  for (i=0;i<(int)psr[p].dmoffsDMnum;i++){
	    printf("_DM\t % 7.1f % 10.3g % 10.3g % 5.2f\n",psr[p].dmoffsDM_mjd[i],psr[p].dmoffsDM[i],psr[p].dmoffsDM_error[i],psr[p].dmoffsDM_weight[i]);
		if(dmfile)fprintf(dmfile,"% 7.1f % 15.7g % 15.7g % 7.4f\n",psr[p].dmoffsDM_mjd[i],psr[p].dmoffsDM[i],psr[p].dmoffsDM_error[i],psr[p].dmoffsDM_weight[i]);
		sum+=psr[p].dmoffsDM[i];
	  }
	  if(dmfile)fclose(dmfile);
      printf("Mean DMOFF = %lg\n",sum/(double)psr[p].dmoffsDMnum);

	  sprintf(dmfname,"%s.cm",psr[p].name);
	  dmfile=fopen(dmfname,"w");
	  sum=0;
	  double eee;
	  for (i=0;i<(int)psr[p].dmoffsCMnum;i++){
	    printf("_CM\t % 7.1f % 10.3g % 10.3g % 5.2f\n",psr[p].dmoffsCM_mjd[i],psr[p].dmoffsCM[i],psr[p].dmoffsCM_error[i],psr[p].dmoffsCM_weight[i]);
		if(dmfile)fprintf(dmfile,"% 7.1f % 15.7g % 15.7g % 7.4f\n",psr[p].dmoffsCM_mjd[i],psr[p].dmoffsCM[i],psr[p].dmoffsCM_error[i],psr[p].dmoffsCM_weight[i]);
		eee=psr[p].dmoffsCM_error[i]/86400.0/365.25;
		sum+=1.0/(eee*eee);
	  }
	  double tobs=(psr[p].dmoffsCM_mjd[psr[p].dmoffsCMnum-1]-psr[p].dmoffsCM_mjd[0])/365.25;
	  if(sum > 0)
	  printf("CM white PSD estimate: %lg\n",2*tobs/sum);
	  if(dmfile)fclose(dmfile);

	}
      if (psr[p].param[param_clk_offs].paramSet[0]==1)
	{
	  FILE *fout;
	  printf("Clock offsets\n");
	  for (i=0;i<psr[p].clkOffsN;i++)
	    printf("%.2f %.10g %.10g\n",psr[p].clk_offsT[i],psr[p].clk_offsV[i],psr[p].clk_offsE[i]);
	  fout = fopen("clockOffset.dat","w");
	  for (i=0;i<psr[p].clkOffsN;i++)
	    fprintf(fout,"%.2f %.10g %.10g\n",psr[p].clk_offsT[i],psr[p].clk_offsV[i],psr[p].clk_offsE[i]);
	  fclose(fout);
	  
	}
	  if (psr[p].param[param_ifunc].paramSet[0]==1)
	    {
	      printf("Interpolated function\n");
	      for (i=0;i<psr[p].ifuncN;i++)
		printf("%.2f %.10g %.10g %.3f\n",psr[p].ifuncT[i],psr[p].ifuncV[i],psr[p].ifuncE[i],psr[p].ifunc_weights[i]);
	    }
	  if (psr[p].param[param_tel_dx].paramSet[0]==1)
	    {
	      FILE *fout;
	      printf("Telescope x function\n");
	      for (i=0;i<psr[p].nTelDX;i++)
		printf("%.2f %.10g %.10g\n",psr[p].telDX_t[i],psr[p].telDX_v[i],psr[p].telDX_e[i]);
	      fout = fopen("telescopeXYZ.dat","w");
	      for (i=0;i<psr[p].nTelDX;i++)
		fprintf(fout,"%.2f %.10g %.10g %.10g %.10g %.10g %.10g\n",psr[p].telDX_t[i],psr[p].telDX_v[i],psr[p].telDX_e[i],psr[p].telDY_v[i],psr[p].telDY_e[i],psr[p].telDZ_v[i],psr[p].telDZ_e[i]);
	      fclose(fout);
	    }
	  if (psr[p].param[param_tel_dy].paramSet[0]==1)
	    {
	      printf("Telescope y function\n");
	      for (i=0;i<psr[p].nTelDY;i++)
		printf("%.2f %.10g %.10g\n",psr[p].telDX_t[i],psr[p].telDY_v[i],psr[p].telDY_e[i]);
	    }
	  if (psr[p].param[param_tel_dz].paramSet[0]==1)
	    {
	      printf("Telescope z function\n");
	      for (i=0;i<psr[p].nTelDZ;i++)
		printf("%.2f %.10g %.10g\n",psr[p].telDZ_t[i],psr[p].telDZ_v[i],psr[p].telDZ_e[i]);
	    }
	  if (psr[p].param[param_quad_ifunc_p].paramSet[0]==1)
	    {
	      FILE *fout;
	      int fileout=1;
	      if (!(fout = fopen("aplus_t2.dat","w")))
		fileout=0;
	      printf("---------------------------------------\n");
	      printf("Interpolated quadrupolar plus function\n");
	      printf("---------------------------------------\n");
	      for (i=0;i<psr[p].quad_ifuncN_p;i++)
		{
		  printf("%.2f %.10g %.10g\n",psr[p].quad_ifuncT_p[i],psr[p].quad_ifuncV_p[i],psr[p].quad_ifuncE_p[i]);
		  if (fileout==1)
		    fprintf(fout,"%.2f %.10g %.10g\n",psr[p].quad_ifuncT_p[i],psr[p].quad_ifuncV_p[i],psr[p].quad_ifuncE_p[i]);
		}
	      if (fileout==1) fclose(fout);
	    }
	  if (psr[p].param[param_quad_ifunc_c].paramSet[0]==1)
	    {
	      FILE *fout;
	      int fileout=1;
	      if (!(fout = fopen("across_t2.dat","w")))
		fileout=0;

	      printf("---------------------------------------\n");
	      printf("Interpolated quadrupolar cross function\n");
	      printf("---------------------------------------\n");
	      for (i=0;i<psr[p].quad_ifuncN_c;i++)
		{
		  printf("%.2f %.10g %.10g\n",psr[p].quad_ifuncT_c[i],psr[p].quad_ifuncV_c[i],psr[p].quad_ifuncE_c[i]);
		  if (fileout==1)
		    fprintf(fout,"%.2f %.10g %.10g\n",psr[p].quad_ifuncT_c[i],psr[p].quad_ifuncV_c[i],psr[p].quad_ifuncE_c[i]);
		}
	      if (fileout==1) fclose(fout);
	    }

      if (psr[0].param[param_f].paramSet[0]==1 && psr[0].param[param_f].paramSet[1]==1)
	{
	  longdouble p0,p1,age,bs,p0e,p1e;
	  longdouble f0,f1,f0e,f1e;





	  printf("------------------------------------------------------------------------------\n");
	  printf("\nDerived parameters:\n\n");

	  f0 = psr[0].param[param_f].val[0];
	  f1 = psr[0].param[param_f].val[1];
	  f0e = psr[0].param[param_f].err[0];
	  f1e = psr[0].param[param_f].err[1];
	  
	  p0 = (1.0/f0);
	  p0e = 1.0/f0/f0*f0e;
	  p1 = (-1.0/(f0*f0)*f1);
	  p1e = sqrt(pow(p0*p0*f1e,2)+pow(2*p0*p0*p0*f1*f0e,2));
	  age = p0/2.0/p1/60.0/60.0/24.0/365.0/1.0e6;
	  bs = (sqrt(-f1/pow(f0,3))*3.2e19);
	  
	  printf("P0 (s)      = %-25.15g %-13.5g\n",(double)p0,(double)p0e);
	  printf("P1          = %-25.15g %-13.5g\n",(double)p1,(double)p1e);
	  printf("tau_c (Myr) = %.5g\n",(double)age);
	  printf("bs (G)      = %.5g\n\n",(double)bs);
	}
      /* Binary parameters */
      if (psr[p].param[param_pb].paramSet[0]==1 || psr[p].param[param_fb].paramSet[0]==1)
	{
	  double err;
	  longdouble pb,a1,pbe,a1e,si=-2,si_lo=-2,si_hi=-2;

	  pb  = psr[p].param[param_pb].val[0]*SECDAY;
	  pbe = psr[p].param[param_pb].err[0]*SECDAY;
	  a1  = psr[p].param[param_a1].val[0]*SPEED_LIGHT;
	  a1e = psr[p].param[param_a1].err[0]*SPEED_LIGHT;

          printf("Binary model: %s\n",psr[p].binaryModel);

	  /* Display mass function */
	  if (psr[p].param[param_a1].paramSet[0]==1)
	    {
	      int c1=0,c2=0,j,a,count=0;
	      longdouble fn;

	      fn  = 4*M_PI*M_PI/GM*pow(a1,3)/pow(pb,2);	      
	      printf("Mass function                  = %.12f ",(double)fn);

	      for (j=0;j<MAX_PARAMS;j++)
		{
		  for (a=0;a<psr[p].param[j].aSize;a++)
		    {
		      if (j==param_a1 && a==0)
			c1 = count;
		      if (j==param_pb && a==0)
			c2 = count;
		      if (psr[p].param[j].fitFlag[a]==1)
			count++;
		    }
		}
	      /* use covariance matrix to determine error */
	      err = fn * sqrt(9.0*pow(sqrt(psr[p].covar[c1][c1])*SPEED_LIGHT/a1,2) 
			      +  4.0*pow(sqrt(psr[p].covar[c2][c2])/pb,2) 
			      - 12.0*psr[p].covar[c1][c2]*SPEED_LIGHT/(a1*pb));
	      printf(" +- %.12f solar masses \n",err);
	      printf("Minimum, median and maximum companion mass: %.4f < %.4f < %.4f solar masses\n",
		     TNm2(fn,1.0,1.35),TNm2(fn,0.866025403,1.35),TNm2(fn,0.4358898944,1.35));

	      
// M2 and SINI from DDH model (FW10)
        if( psr[p].param[param_h3].paramSet[0] == 1 ){ 
          long double h3 = psr[p].param[param_h3].val[0];
          long double m2, sini;
          if( psr[p].param[param_stig].paramSet[0] == 1 ){
            long double stig = psr[p].param[param_stig].val[0];
            // Freire & Wex, Eq. 20:
            m2 = h3 / 4.925490947e-6 * pow( stig, -3.0 );
            // Freire & Wex, Eq. 22:
            sini = 2.0 * stig / ( 1.0 + pow( stig, 2.0 ) );
            printf( " DDH-model Derived parameters:  M2 = %lg\n", (double)m2 );
            printf( "                              SINI = %lg\n", (double)sini );
          }else if( psr[p].param[param_h4].paramSet[0] == 1 ){
            long double h4 = psr[p].param[param_h4].val[0];
            // Freire & Wex, Eq. 26:
            if( h4 != 0.0 )
              m2 = pow( h3, 4.0 ) / pow( h4, 3.0 ) / 4.925490947e-6;
            // Freire & Wex, Eq. 25:
            sini = 2.0 * h3 * h4 / ( h3 * h3 + h4 * h4 );
            printf( " DDH-model Derived parameters:  M2 = %lg\n", (double)m2 );
            printf( "                              SINI = %lg\n", (double)sini );
          }
        }        

	      // Joris' mass calculations.

	      if(psr[p].param[param_sini].paramSet[0]*
           psr[p].param[param_m2].paramSet[0]*
           psr[p].param[param_a1].paramSet[0]*
           psr[p].param[param_pb].paramSet[0]==1
           && psr[p].param[param_sini].nLinkTo==0){
          long double mp[2];
          double DAY2S = (24.0L*3600.0L);
          mp[0] = -psr[p].param[param_m2].val[0]+
            sqrt(TSUN*pow(psr[p].param[param_pb].val[0]*DAY2S/2.0/M_PI,2.0)*
                 pow(psr[p].param[param_m2].val[0]*
                     psr[p].param[param_sini].val[0]/
                     psr[p].param[param_a1].val[0],3.0));
          
          longdouble Cte = sqrt(TSUN*pow(1/2.0L/M_PI,2.0));
          mp[1] = sqrt(pow(psr[p].param[param_m2].err[0]*
                           (-1.0+1.5*Cte*DAY2S*psr[p].param[param_pb].val[0]*
                            pow(psr[p].param[param_sini].val[0]/
                                psr[p].param[param_a1].val[0],1.5)*
                            sqrt(psr[p].param[param_m2].val[0])),2.0)+
                       pow(psr[p].param[param_pb].err[0]*Cte*DAY2S*
                           pow(psr[p].param[param_m2].val[0]*
                               psr[p].param[param_sini].val[0]/
                               psr[p].param[param_a1].val[0],1.5),2.0)+
                       pow(psr[p].param[param_sini].err[0]*1.5L*
                           sqrt(psr[p].param[param_sini].val[0])*Cte*
                           psr[p].param[param_pb].val[0]*DAY2S*
                           pow(psr[p].param[param_m2].val[0]/
                               psr[p].param[param_a1].val[0],1.5),2.0)+
                       pow(psr[p].param[param_a1].err[0]*1.5/
                           pow(psr[p].param[param_a1].val[0],2.5)*Cte*
                           psr[p].param[param_pb].val[0]*DAY2S*
                           pow(psr[p].param[param_m2].val[0]*
                               psr[p].param[param_sini].val[0],1.5),2.0));
          

		printf("Pulsar Mass (Shapiro Delay): %lg (+/- %lg) Msun.\n",
           (double)mp[0],(double)mp[1]);
	      }
	      if(psr[p].param[param_kin].paramSet[0]*
           psr[p].param[param_m2].paramSet[0]*
           psr[p].param[param_a1].paramSet[0]*
           psr[p].param[param_pb].paramSet[0]==1){
          longdouble mp[2];
          double DAY2S = (24.0L*3600.0L);
          mp[0] = -psr[p].param[param_m2].val[0]+
            sqrt(TSUN*pow(psr[p].param[param_pb].val[0]*DAY2S/2.0/M_PI,2.0)*
                 pow(psr[p].param[param_m2].val[0]*
                     sin(psr[p].param[param_kin].val[0]/180.0*M_PI)/
                     psr[p].param[param_a1].val[0],3.0));
          
          longdouble Cte = sqrt(TSUN*pow(1/2.0L/M_PI,2.0));
          mp[1] = sqrt(pow(psr[p].param[param_m2].err[0]*
                           (-1.0+1.5*Cte*DAY2S*psr[p].param[param_pb].val[0]*
                            pow(sin(psr[p].param[param_kin].val[0]/180.0*M_PI)/
                                psr[p].param[param_a1].val[0],1.5)*
                            sqrt(psr[p].param[param_m2].val[0])),2.0)+
                       pow(psr[p].param[param_pb].err[0]*Cte*DAY2S*
                           pow(psr[p].param[param_m2].val[0]*
                               sin(psr[p].param[param_kin].val[0]/180.0*M_PI)/
                               psr[p].param[param_a1].val[0],1.5),2.0)+
                       pow(psr[p].param[param_kin].err[0]/180.0*M_PI*1.5L*
                           cos(psr[p].param[param_kin].val[0]/180.0*M_PI)*
                           sqrt(sin(psr[p].param[param_kin].val[0]
                                    /180.0*M_PI))*Cte*
                           psr[p].param[param_pb].val[0]*DAY2S*
                           pow(psr[p].param[param_m2].val[0]/
                               psr[p].param[param_a1].val[0],1.5),2.0)+
                       pow(psr[p].param[param_a1].err[0]*1.5/
                           pow(psr[p].param[param_a1].val[0],2.5)*Cte*
                           psr[p].param[param_pb].val[0]*DAY2S*
                           pow(psr[p].param[param_m2].val[0]*
                               sin(psr[p].param[param_kin].val[0]
                                   /180.0*M_PI),1.5),2.0));
          
          
          printf("Pulsar Mass (annual orbital parallax): %lg (+/- %lg) Msun.\n",
                 (double)mp[0],(double)mp[1]);
	      }

	      /* Joris' distance calculations */
	      longdouble transV[4]; 
	      /* transverse velocity. 
		 0: pxdistvalue 
		 1: pxdisterror 
		 2: pbdotval
		 3: pbdoterr */
	      if(psr[p].param[param_px].paramSet[0]==1){
		longdouble pxdist[2]; // 0: value; 1: error
		longdouble pmsqrd[2];
		pxdist[0] = (longdouble)(1.0/psr[p].param[param_px].val[0]*1000.0);
		pxdist[1] = 1/powl(psr[p].param[param_px].val[0],2.0)*
		  psr[p].param[param_px].err[0]*1000.0;
		printf("\nParallax distance is %lg (+/- %lg) pc.\n",(double)pxdist[0],
		       (double)pxdist[1]);
		pmsqrd[0] = powl(psr[p].param[param_pmra].val[0]*MASYR2RADS,2.0)+
		  powl(psr[p].param[param_pmdec].val[0]*MASYR2RADS,2.0);
		pmsqrd[1] = powl(2*psr[p].param[param_pmra].val[0]*
				 psr[p].param[param_pmra].err[0]*powl(MASYR2RADS,2.0),2.0)
		  +powl(2*psr[p].param[param_pmdec].val[0]*psr[p].param[param_pmdec].err[0]
			*powl(MASYR2RADS,2.0),2.0);
		/*transV[0] = sqrtl(pmsqrd[0])*pxdist[0]*PCM/1000.0; // now in km/s
		transV[1] = sqrtl(powl(pmsqrd[1]*pxdist[0]*PCM/(2.0*pmsqrd[0]),2.0)+
				  pmsqrd[0]*powl(pxdist[1]*PCM,2.0))/1000.0;
		printf("\tTransverse velocity based on parallax distance: %lg +/- %lg km/s.\n",
		       (double)transV[0],(double)transV[1]);*/
	      }

	      if(psr[p].param[param_pbdot].paramSet[0]==1){
		longdouble pbdotdist[2];
		longdouble pmsqrd[2];
		pmsqrd[0] = powl(psr[p].param[param_pmra].val[0]*MASYR2RADS,2.0)+
		  powl(psr[p].param[param_pmdec].val[0]*MASYR2RADS,2.0);
		pmsqrd[1] = powl(2*psr[p].param[param_pmra].val[0]*
				 psr[p].param[param_pmra].err[0]*powl(MASYR2RADS,2.0),2.0)+
		  powl(2*psr[p].param[param_pmdec].val[0]*
		       psr[p].param[param_pmdec].err[0]*powl(MASYR2RADS,2.0),2.0);
		pbdotdist[0] = psr[p].param[param_pbdot].val[0]*SPEED_LIGHT/
		  (psr[p].param[param_pb].val[0]*SECDAY)/
		  pmsqrd[0]/PCM;
		pbdotdist[1] = sqrtl(powl(SPEED_LIGHT/(psr[p].param[param_pb].val[0]
						       *SECDAY*pmsqrd[0])
					  *psr[p].param[param_pbdot].err[0],2.0)+
				     powl(psr[p].param[param_pb].err[0]*SECDAY*SPEED_LIGHT*
					  psr[p].param[param_pbdot].val[0]/
					  (pmsqrd[0]*powl(psr[p].param[param_pb].val[0]
							  *SECDAY,2.0)),2.0)+
				     powl(pmsqrd[1]*SPEED_LIGHT*psr[p].param[param_pbdot].val[0]/
					  (psr[p].param[param_pb].val[0]*SECDAY*
					   powl(pmsqrd[0],2.0)),2.0))/PCM;
		printf("Pbdot distance is %lg (+/- %lg) pc.\n", 
		       (double)pbdotdist[0],(double)pbdotdist[1]);
		
		/*transV[2] = sqrtl(pmsqrd[0])*pbdotdist[0]*PCM/1000.0; // now in km/s
		transV[3] = sqrtl(powl(pmsqrd[1]*pbdotdist[0]*PCM/(2.0*pmsqrd[0]),2.0)+
				  pmsqrd[0]*powl(pbdotdist[1]*PCM,2.0))/1000.0;
		printf("\tTransverse velocity based on pbdot distance: %lg +/- %lg km/s.\n\n",
		       (double)transV[2],(double)transV[3]);*/		
	      }
	      if(psr[p].param[param_daop].paramSet[0]==1)
		printf("Used daop distance of %g for Annual-Orbital parallax delays.\n",
		       (double)getParameterValue(&psr[p],param_daop,0));
	    }
	  /* Calculate sin i */
	  if (psr[p].param[param_shapmax].paramSet[0]==1)
	    {
	      longdouble smax;
	      double err;

	      smax = psr[p].param[param_shapmax].val[0];
	      err  = psr[p].param[param_shapmax].err[0];

	      si = 1.0 - exp(-1.0*smax);
	      printf("sini derived from SHAPMAX      = %.14Lf ",si);

	      si_lo = 1.0-exp(-1.0*(smax-err));
	      si_hi = 1.0-exp(-1.0*(smax+err));
	      printf(" (+ %.8f  - %.8f)\n",(double)(si_hi-si),(double)(si-si_lo));

	    }
	  /* mtot derived from m2 and sini */
	  if (psr[p].param[param_m2].paramSet[0]==1 && (si!=-2 || psr[p].param[param_sini].paramSet[0]==1))
	    {
	      longdouble m2,amtot;

	      if (si==-2)
		si = getParameterValue(&psr[p],param_sini,0);
	      
	      m2 = psr[p].param[param_m2].val[0]*SOLAR_MASS;
	      amtot = pb/2.0/M_PI*sqrt(pow(m2*si,3)*BIG_G/pow(a1,3))/SOLAR_MASS;
	      printf("MTOT derived from sin i and m2 = %.14g\n",(double)amtot);
	    }
	  /* inclination angle derived from sini */
	  if ((si!=-2 || psr[p].param[param_sini].paramSet[0]==1)
                && psr[p].param[param_sini].nLinkTo == 0) 
	    {
	      longdouble inc,inc_lo,inc_hi;
	      if (si==-2)
		si = getParameterValue(&psr[p],param_sini,0);
	      if (si_lo==-2)
		si_lo = getParameterValue(&psr[p],param_sini,0)-psr[p].param[param_sini].err[0];
	      if (si_hi==-2)
		si_hi = getParameterValue(&psr[p],param_sini,0)+psr[p].param[param_sini].err[0];
	      
	      inc    = asin(si)*180.0/M_PI;
	      inc_lo = asin(si_lo)*180.0/M_PI;
	      inc_hi = asin(si_hi)*180.0/M_PI;
	      printf("Inclination angle (deg)        = %.14g (+ %.7g - %.7g)\n",
		     (double)inc,(double)(inc_hi-inc),(double)(inc-inc_lo));
	    }
	  // Print out OM and ECC if using the ELL1 model
	  if (psr[p].param[param_eps1].paramSet[0]==1 &&
	      psr[p].param[param_eps2].paramSet[0]==1)
	    {
	      long double om,ecc,t0,pb;
	      long double ecc_err,om_err,t0_err;
	      long double eps1,eps2,tasc,err1,err2,err3;

	      eps1 = psr[p].param[param_eps1].val[0];
	      eps2 = psr[p].param[param_eps2].val[0];
	      tasc = psr[p].param[param_tasc].val[0];
	      pb   = psr[p].param[param_pb].val[0];
	      err1 = psr[p].param[param_eps1].err[0];
	      err2 = psr[p].param[param_eps2].err[0];
	      err3 = psr[p].param[param_tasc].err[0];

	      ecc = sqrt(powl(eps1,2.0)+powl(eps2,2.0));
	      om  = atan2(eps1,eps2)*180.0/M_PI;
	      if (om<0) om+=360.0;

	      t0  = tasc + pb/360*om;
	      om_err  = pow(pow(eps2*err1,2)+pow(eps1*err2,2),0.5)/ecc/ecc*180.0/M_PI;		  
	      /* What should this be if EPS1 and EPS2 = 0 */
	      if (eps1==0.0 && eps2==0.0)
		ecc_err = sqrt(pow(err1,2)+pow(err2,2));
	      else
		ecc_err = pow(pow(eps1*err1,2)+pow(eps2*err2,2),0.5)/ecc;

	      t0_err = sqrt(err3*err3+pow(pb/360.0*om_err,2));

	      printf("\nConversion from ELL1 parameters: \n");
	      printf("------------------------------------\n");
	      printf("ECC = %.15lg +/- %.15g\n",(double)ecc,(double)ecc_err);
	      printf("OM  = %.15g +/- %.15g degrees\n",(double)om,(double)om_err);
	      printf("T0  = %.15Lg +/- %.15g\n",t0,(double)t0_err);
	      printf("------------------------------------\n");
	    }
	}
      if (psr[p].param[param_pmra].paramSet[0]==1 &&
	  psr[p].param[param_pmdec].paramSet[0]==1)
	{
	  double pmtot,epmtot,val1,val2,err1,err2;
	  val1 = psr[p].param[param_pmra].val[0];
	  val2 = psr[p].param[param_pmdec].val[0];
	  err1 = psr[p].param[param_pmra].err[0];
	  err2 = psr[p].param[param_pmdec].err[0];
	  pmtot = sqrt(pow(val1,2)+pow(val2,2));
	  epmtot = sqrt((pow(val1*err1,2)+pow(val2*err2,2))/(val1*val1 + val2*val2));
	  printf("Total proper motion = %.5g +/- %.5g mas/yr\n",pmtot,epmtot);
	}
      // Glitch parameters
      if (psr[p].param[param_glep].paramSet[0]==1)
	TNprintGlitch(psr[p]);
      if (psr[p].param[param_dmassplanet].paramSet[4]==1)
	{
	  long double diff,err;
	  if (strstr(psr[p].JPL_EPHEMERIS,"DE405")!=NULL)
	    {
	      printf("M_Jupiter fit   = %.15Lg +/- %.15Lg (Solar masses)\n",psr[p].param[param_dmassplanet].val[4]+0.00095479193842432214L,psr[p].param[param_dmassplanet].err[4]);
	      printf("M_Jupiter DE405 = %.15Lg (Solar masses)\n",0.00095479193842432214L);
	      printf("M_Jupiter best  = 0.000954791915(11) (Solar masses)\n");
	      diff = psr[p].param[param_dmassplanet].val[4]+0.00095479193842432214L-0.000954791915;
	      err = sqrtl(powl(psr[p].param[param_dmassplanet].err[4],2)+powl(0.000000000011,2));
								    
	      printf("diff            = %.15Lg +/- %.15Lg (Solar masses)\n",diff,err);
	    }
	  else if (strstr(psr[p].JPL_EPHEMERIS,"DE200")!=NULL)
	    {
	      printf("DE200 Jupiter mass not set\n");
	    }
	}

      // Calculate time span
      {
	double start,end;
	int go=0;
	for (i=0;i<psr[p].nobs;i++)
	  {
	    if (psr[p].obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)psr[p].obsn[i].sat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)psr[p].obsn[i].sat)
		      start = (double)psr[p].obsn[i].sat;
		    if (end < (double)psr[p].obsn[i].sat)
		      end = (double)psr[p].obsn[i].sat;
		  }
	      }
	  }
	printf("Total time span = %.3f days = %.3f years\n",end-start,(end-start)/365.25);
      }
    
      if (1==1)  /* Write a new .par file */
	{
		
		std::string parname=longname+".par";
		printf("name size: %i \n",parname.size());
		fname=(char*)parname.c_str();
		//fname="newpar.par";
		printf("writing par file %s. \n",fname);
	  FILE *fout2;
	  char fname2[1000];
	  char str1[100],str2[100],str3[100],str4[100],str5[100];
	  int nread;

		char hmsstr[100];
		longturn_hms(psr[p].param[param_raj].val[0]/(2*M_PI), hmsstr);
		

		strcpy(psr[p].rajStrPost,hmsstr);
			
		longturn_dms(psr[p].param[param_decj].val[0]/(2*M_PI), hmsstr);
		strcpy(psr[p].decjStrPost,hmsstr);

	  if (strlen(fname)==0)
	    {
	      printf("Enter filename for new parameter file ");
	      scanf("%s",fname2);
	    }
	  else
	    {
	      char fname3[1000];
	      if (npsr > 1)
		sprintf(fname2,"%s_%d",fname,p+1);
	      else
		strcpy(fname2,fname);
	    }
	  if (!(fout2 = fopen(fname2,"w")))
	    {
	      printf("Sorry, unable to write to file %s\n",fname2);
	    }
	  else
	    {
	      fprintf(fout2,"%-15.15s%s\n","PSRJ",psr[p].name);
	      for (i=0;i<MAX_PARAMS;i++)
		{
		  for (k=0;k<psr[p].param[i].aSize;k++)
		    {
		  if (psr[p].param[i].paramSet[k]==1 && i!=param_wave_om 
		      && i!=param_waveepoch && i!=param_ifunc && i!=param_dmmodel &&
		      (psr[p].tempo1==0 || (i!=param_dmepoch)))
		    {
		      if (strcmp(psr[p].param[i].shortlabel[k],"PB")==0 || strcmp(psr[p].param[i].shortlabel[k],"FB0")==0)
			fprintf(fout2,"%-15.15s%s\n","BINARY",psr[p].binaryModel);
		      
		      if (i == param_raj && psr[p].eclCoord==1)
			fprintf(fout2,"%-15.15s","ELONG");
		      else if (i == param_decj && psr[p].eclCoord==1)
			fprintf(fout2,"%-15.15s","ELAT");
		      else if (i == param_pmra && psr[p].eclCoord==1)
			fprintf(fout2,"%-15.15s","PMELONG");
		      else if (i == param_pmdec && psr[p].eclCoord==1)
			fprintf(fout2,"%-15.15s","PMELAT");		      
		      else
			{
			  if (i==param_track && psr[p].param[i].val[k]==0)
			    {// Do nothing
			    }
			  else fprintf(fout2,"%-15.15s",psr[p].param[i].shortlabel[k]);
			}
		      if (psr[p].eclCoord==0 && i==param_raj)
			fprintf(fout2,"%-25.25s",psr[p].rajStrPost);
		      else if (psr[p].eclCoord==1 && i==param_raj)
			fprintf(fout2,"%-25.25Lf",psr[p].param[i].val[k]*180.0/M_PI);
		      else if (psr[p].eclCoord==0 && i==param_decj)
			fprintf(fout2,"%-25.25s",psr[p].decjStrPost);
		      else if (psr[p].eclCoord==1 && i==param_decj)
			fprintf(fout2,"%-25.25Lf",psr[p].param[i].val[k]*180.0/M_PI);
		      else if (i == param_sini && psr[p].param[i].nLinkTo>0){
			fprintf(fout2," KIN\n"); 
			fprintf(fout2,"#SINI\t\t%-25.20Lg",psr[p].param[i].val[k]);
		      }
		      else if (i==param_tres)
			fprintf(fout2,"%-10.3Lf",psr[p].param[i].val[k]);
		      else if (i==param_tzrfrq)
			{
			  fprintf(fout2,"%-25.20Lg\n",psr[p].param[i].val[k]);
			  if (strcmp(psr[p].tzrsite,"NULL")!=0) 
			    fprintf(fout2,"%-15.15s%s","TZRSITE",psr[p].tzrsite);
			}
		      else 
			{
			  if (i==param_track && psr[p].param[i].val[k]==0)
			    {
			      // Do nothing
			    }
			  else
			    fprintf(fout2,"%-25.20Lg",psr[p].param[i].val[k]);
			}
		      if (psr[p].param[i].fitFlag[k]==1)
			{
			  fprintf(fout2," 1 ");
			  if (i==param_raj)
			    {
			      if (psr[p].param[i].err[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].err[k]*12.0*60.0*60.0/M_PI);
			      else if (psr[p].param[i].err[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].err[k]*12.0*60.0*60.0/M_PI);
			    }
			  else if (i==param_decj)
			    {
			      if (psr[p].param[i].err[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].err[k]*180.0*60.0*60.0/M_PI);
			      else if (psr[p].param[i].err[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].err[k]*180.0*60.0*60.0/M_PI);
			    }
			  else
			    {
			      if (psr[p].param[i].err[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].err[k]);
			      else if (psr[p].param[i].err[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].err[k]);
			    }
			}
		      else
			{
			  fprintf(fout2,"   ");
			  if (i==param_raj)
			    {
			      if (psr[p].param[i].prefitErr[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].prefitErr[k]*12.0*60.0*60.0/M_PI);
			      else if (psr[p].param[i].prefitErr[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].prefitErr[k]*12.0*60.0*60.0/M_PI);
			    }
			  else if (i==param_decj)
			    {
			      if (psr[p].param[i].prefitErr[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].prefitErr[k]*180.0*60.0*60.0/M_PI);
			      else if (psr[p].param[i].prefitErr[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].prefitErr[k]*180.0*60.0*60.0/M_PI);
			    }
			  else
			    {
			      if (psr[p].param[i].prefitErr[k]>1e-12) fprintf(fout2," %-25.20Lf",psr[p].param[i].prefitErr[k]);
			      else if (psr[p].param[i].prefitErr[k]>0) fprintf(fout2," %-25.20Lg",psr[p].param[i].prefitErr[k]);
			    }
			}
		      fprintf(fout2,"\n");
		    }
		}
	    }
	  if (psr[p].tempo1 == 1)
	    {
	      if (!strcmp(psr[p].clock, "TT(UTC(NIST))")) fprintf(fout2,"%-15.15s%s\n","CLK","UTC(NIST)");
	      else if (!strcmp(psr[p].clock, "TT(TAI))")) fprintf(fout2,"%-15.15s%s\n","CLK","UTC(BIPM)");
	      else if (!strcmp(psr[p].clock, "TT(UTC(PTB))")) fprintf(fout2,"%-15.15s%s\n","CLK","PTB");
	      else if (!strcmp(psr[p].clock, "TT(TA(NIST))")) fprintf(fout2,"%-15.15s%s\n","CLK","AT1");
	      else fprintf(fout2,"%-15.15s%s\n","CLK",psr[p].clock);
	    }
	  else	   
	    fprintf(fout2,"%-15.15s%s\n","CLK",psr[p].clock);

	  if (psr[p].fitMode==1) fprintf(fout2,"MODE 1\n");
	  if (psr[p].units != SI_UNITS)
	    fprintf(fout2, "%-15.15s%s\n", "UNITS", "TDB");
	  if (psr[p].timeEphemeris != IF99_TIMEEPH)
	    fprintf(fout2, "%-15.15s%s\n", "TIMEEPH", "FB90");
	  if (!psr[p].dilateFreq)
	    fprintf(fout2, "%-15.15s%s\n", "DILATEFREQ", "N");
	  if (!psr[p].planetShapiro)
	    fprintf(fout2, "%-15.15s%s\n", "PLANET_SHAPIRO", "N");
	  if (psr[p].t2cMethod != T2C_IAU2000B)
	    fprintf(fout2, "%-15.15s%s\n", "T2CMETHOD", "TEMPO");
	  if (!psr[p].correctTroposphere)
	    fprintf(fout2, "%-21.21s%s\n", "CORRECT_TROPOSPHERE", "N");

	  fprintf(fout2,"%-15.15s%s\n","EPHEM",psr[p].ephemeris);
	  fprintf(fout2,"%-15.15s%s\n","NITS","1");
	  fprintf(fout2,"%-15.15s%d\n","NTOA",psr[p].nFit);
	  fprintf(fout2,"%-15.15s%.4f %d\n","CHI2R",chisqr/(double)psr[p].fitNfree,psr[p].fitNfree);
	  /*	  if (psr[p].tempo1 == 1)
	    fprintf(fout2,"EPHVER         2\n");
	  else
	  fprintf(fout2,"EPHVER         5\n"); */

	  /* Add jumps */
	  for (i=1;i<=psr[p].nJumps;i++)
	    {
	      nread = sscanf(psr[p].jumpStr[i],"%s %s %s %s %s",str1,str2,str3,str4,str5);
	      if (strcasecmp(str1,"FREQ")==0 || strcasecmp(str1,"MJD")==0)
		fprintf(fout2,"JUMP %s %s %s %.14g %d\n",str1,str2,str3,psr[p].jumpVal[i],psr[p].fitJump[i]);
	      else if (strcasecmp(str1,"NAME")==0 || strcasecmp(str1,"TEL")==0 || str1[0]=='-')
		fprintf(fout2,"JUMP %s %s %.14g %d\n",str1,str2,psr[p].jumpVal[i],psr[p].fitJump[i]);
	    }	  

	  /* Add whitening flags */
	  if (psr[p].param[param_wave_om].paramSet[0]==1)
	    {
	      fprintf(fout2,"WAVEEPOCH %.14Lg\n",psr[p].param[param_waveepoch].val[0]);
	      fprintf(fout2,"WAVE_OM %.14Lg 0\n",psr[p].param[param_wave_om].val[0]);
	      if (psr[p].waveScale!=0) fprintf(fout2,"WAVE_SCALE %g\n",psr[p].waveScale);
	      for (i=0;i<psr[p].nWhite;i++)
		fprintf(fout2,"WAVE%d %.14g %.14g\n",i+1,psr[p].wave_sine[i],psr[p].wave_cos[i]);
	    }
 
	  if (psr[p].param[param_ifunc].paramSet[0]==1)
	    {
	      fprintf(fout2,"SIFUNC %d %d\n",(int)psr[p].param[param_ifunc].val[0],
		      (int)psr[p].param[param_ifunc].fitFlag[0]);
	      for (i=0;i<psr[p].ifuncN;i++)
		fprintf(fout2,"IFUNC%d %.15f %.10f %.10f\n",i+1,psr[p].ifuncT[i],
			psr[p].ifuncV[i],psr[p].ifuncE[i]);
	    }
	  /* Add phase jumps */
	  for (i=0;i<psr[p].nPhaseJump;i++)
	    {
	      if (psr[p].phaseJumpDir[i]!=0)
		fprintf(fout2,"PHASE %+d %.14g\n",psr[p].phaseJumpDir[i],(double)(psr[p].obsn[psr[p].phaseJumpID[i]].sat+1.0/SECDAY));
	    }
	  // Add DM value parameters
	  if (psr[p].param[param_dmmodel].paramSet[0]==1)
	    {

	      if (psr->param[param_dmmodel].linkTo[0] == param_dm){
		      fprintf(fout2,"DMMODEL DM %d\n",(int)psr[p].param[param_dmmodel].fitFlag[0]);
	      } else {
		      fprintf(fout2,"DMMODEL %.14Lg %d\n",psr[p].param[param_dmmodel].val[0],(int)psr[p].param[param_dmmodel].fitFlag[0]);
	      }
		  bool useDMOFF=psr[p].dmoffsDMnum==psr[p].dmoffsCMnum;
		  if(useDMOFF)for (i=0;i<psr[p].dmoffsDMnum;i++){
			 if(psr[p].dmoffsDM_mjd[i]!=psr[p].dmoffsCM_mjd[i])useDMOFF=false;
			 break;
		  }
		  if(useDMOFF){
		  for (i=0;i<psr[p].dmoffsDMnum;i++)
			 fprintf(fout2,"DMOFF\t %.15g %.15g %.15g\n",psr[p].dmoffsDM_mjd[i],psr[p].dmoffsDM[i],psr[p].dmoffsDM_error[i]);

		  }else{
			 for (i=0;i<psr[p].dmoffsDMnum;i++)
				fprintf(fout2,"_DM\t %.15g %.15g %.15g\n",psr[p].dmoffsDM_mjd[i],psr[p].dmoffsDM[i],psr[p].dmoffsDM_error[i]);
			 for (i=0;i<psr[p].dmoffsCMnum;i++)
				fprintf(fout2,"_CM\t %.15g %.15g %.15g\n",psr[p].dmoffsCM_mjd[i],psr[p].dmoffsCM[i],psr[p].dmoffsCM_error[i]);
		  }
		}

	  // add constraints
	  for (int i = 0; i < psr[p].nconstraints; i++){
		  if (psr[p].constraints[i]==constraint_dmmodel_mean){
			fprintf(fout2,"CONSTRAIN DMMODEL\n");
		  }
	  }
	  fclose(fout2);	 
	    }
	}
      /* printf("Precision: routine, precision, comment\n");
	 for (i=0;i<psr[p].nStorePrecision;i++)
	 {
	 printf("%s\t%Lg\t%s\n",psr[p].storePrec[i].routine,psr[p].storePrec[i].minPrec,
	 psr[p].storePrec[i].comment);
	 } */
    }
  

}

double TNcalcRMS(pulsar *psr,int p)
{
  double sumwt=0.0,rms_pre=0.0,rms_post=0.0,wgt=0.0,sumsq_post=0.0;
  double sumsq_pre=0.0,sum_pre=0.0,sum_post=0.0,mean_post=0.0,mean_pre=0.0;
  int i,count=0;

  for (i=0;i<psr[p].nobs;i++)
    {
      if (psr[p].obsn[i].deleted==0 &&
	  (psr[p].param[param_start].paramSet[0]!=1 || psr[p].param[param_start].fitFlag[0]!=1 ||	  
	   psr[p].param[param_start].val[0] < psr[p].obsn[i].bat) &&
	  (psr[p].param[param_finish].paramSet[0]!=1 || psr[p].param[param_finish].fitFlag[0]!=1 ||
	   psr[p].param[param_finish].val[0] > psr[p].obsn[i].bat))
	{
	  mean_pre += (double)psr[p].obsn[i].prefitResidual;	  
	  if (psr[p].fitMode==1) 
	    wgt = 1.0 /
	      (1.0e-6*psr[p].obsn[i].toaErr*psr[p].param[param_f].val[0]*
		1.0e-6*psr[p].obsn[i].toaErr*psr[p].param[param_f].val[0]);
	  else wgt=1.0/(1.0e-6*psr[p].param[param_f].val[0]*1.0e-6*psr[p].param[param_f].val[0]);
	  sumsq_pre += (double)(wgt*psr[p].obsn[i].prefitResidual*psr[p].param[param_f].val[0]*psr[p].obsn[i].prefitResidual*psr[p].param[param_f].val[0]);
	  sum_pre   += (double)(wgt*psr[p].obsn[i].prefitResidual*psr[p].param[param_f].val[0]);
	  
	  sumsq_post += (double)(wgt*psr[p].obsn[i].residual*psr[p].param[param_f].val[0]*psr[p].obsn[i].residual*psr[p].param[param_f].val[0]);
	  sum_post   += (double)(wgt*psr[p].obsn[i].residual*psr[p].param[param_f].val[0]);
	  sumwt += wgt;
	  mean_post += (double)psr[p].obsn[i].residual;
	  count++;
	}
    }
    logdbg("textOutput %g %g %d",mean_pre,mean_post,count);
  mean_pre/=count;
  mean_post/=count;
  
  rms_pre = sqrt((sumsq_pre-sum_pre*sum_pre/sumwt)/sumwt)*1e3/psr[p].param[param_f].val[0]*1e3;
  rms_post = sqrt((sumsq_post-sum_post*sum_post/sumwt)/sumwt)*1e3/psr[p].param[param_f].val[0]*1e3;
    logdbg("textOutput %g %g %g %G %d",rms_pre,rms_post,sumsq_pre,sum_pre,count);


  psr[p].rmsPre  = rms_pre;
  psr[p].rmsPost = rms_post;
  psr[p].param[param_tres].val[0] = rms_post;
  psr[p].param[param_tres].paramSet[0] = 1;
  return sumwt;
}

/* ************************************************************************
   m2 - solves mass function for m2, using the Newton-Raphson method
   
   where:  mf = mass function
   m1 = primary mass
   si = sin(i) (i = inclination angle)
   
   solves: (m1+m2)^2 = (m2*si)^3 / mf
   
   returns -1 on error
   
   WVS Jan 2000
************************************************************************ */

double TNm2(longdouble mf, longdouble sini, longdouble m1)
{
  double guess = m1;
  double dx = 0.0;
  double eq = 0.0;
  double deq_dm2 = 0.0;
  int gi = 0;

  for (gi=0; gi<10000; gi++) {
    eq = pow(m1+guess,2) - pow(guess*sini,3) / mf;
    deq_dm2 = 2.0*(m1+guess) - 3.0 * pow(guess*sini,2) / mf;

    dx = eq / deq_dm2;
    guess -= dx;

    if (fabs (dx) <= fabs(guess)*1e-10)
      return guess;
  }
  fprintf (stdout,"m2: maximum iterations exceeded - %lf\n", fabs(dx/guess));
  return -1.0;
}

void TNprintGlitch(pulsar psr)
{
  double glep1z,glep2z,glepe;
  int iph;
  double fph;
  double dfof,edfof;
  double df1of1,edf1of1;

  iph = fortran_nint((double)psr.param[param_glph].val[0]);
  fph = (double)psr.param[param_glph].val[0]-iph;

  glep1z=(double)psr.param[param_glep].val[0]+TNdglep(psr,0,fph);
  if (fph >= 0)
    glep2z=(double)psr.param[param_glep].val[0]+TNdglep(psr,0,fph-1);
  else
    glep2z=(double)psr.param[param_glep].val[0]+TNdglep(psr,0,fph+1);

  //  glepe=ferr(NPAR1+(i-1)*NGLP+1)/
  glepe = (double)psr.param[param_glph].err[0]/(double)(fabs(psr.param[param_glf0].val[0]+psr.param[param_glf0d].val[0])*86400.0);

  printf("MJD for zero glitch phase = %.6f or %.6f, error = %g\n",glep1z,glep2z,glepe);

  dfof = (double)(psr.param[param_glf0].val[0]/psr.param[param_f].val[0]);
  edfof = (double)(psr.param[param_glf0].err[0]/psr.param[param_f].val[0]);
  printf("Delta f/f = %g +/- %g\n",dfof,edfof);

  df1of1 = (double)(psr.param[param_glf1].val[0]/psr.param[param_f].val[1]);
  edf1of1 = (double)(psr.param[param_glf1].err[0]/fabs(psr.param[param_f].val[1]));
  printf("Delta f1/f1 = %g +/- %g\n",df1of1,edf1of1);
}

double TNdglep(pulsar psr,int gn,double fph)
{
  double tds,plim,dph,t1;
  int niter;
  
  tds = psr.param[param_gltd].val[gn]*86400.0;
  niter=0;
  plim=1.0e-6;
  dph = 1000.0;
  t1 = -fph/(psr.param[param_glf0].val[gn]+psr.param[param_glf0d].val[gn]);
  do {
    dph = fph + psr.param[param_glf0].val[gn]*t1 + 0.5*psr.param[param_glf1].val[gn]*t1*t1;
    if (tds > 0.0) dph=dph+psr.param[param_glf0d].val[gn]*tds*(1.0-exp(-t1/tds));
    t1=t1-dph/(psr.param[param_glf0].val[gn]+psr.param[param_glf0d].val[gn]);
    niter++;
    if (niter>1000)
      {
	printf("WARNING: Glitch epoch convergence failed\n");
	return 0;
      }
  }while (fabs(dph) > plim);
  return t1/86400.0;
}
//
// Function to determine the error in the DM estimations and apply this to the 
// uncertainty on the TOAs
//
