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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TempoNest.h"

void getLinearPriors(pulsar *psr, double **Dmatrix, long double **LDpriors, double **Dpriors, int numtofit, int fitsig){

	double *linearParams = new double[numtofit];
	double *nonLinearParams = new double[numtofit];
	double *errorvec = new double[numtofit];
		
	nonLinearParams[0]=(double)psr[0].offset_e;
	linearParams[0]=0;
	for(int j=1;j<numtofit;j++){
		nonLinearParams[j]=(double)LDpriors[j][1];
		linearParams[j]=0;

	}
	TNIupdateParameters(psr,0,nonLinearParams,errorvec, linearParams);

	double *Drms=new double[numtofit];
	for(int i =0; i < numtofit; i++){
		Drms[i]=0;
		
		for(int j =0; j < psr->nobs; j++){
			Drms[i] += Dmatrix[j][i]*Dmatrix[j][i];
// 			printf("%i %i %g \n",i,j,Dmatrix[j][i]);
		}
		Drms[i]=sqrt(Drms[i]/psr->nobs);
	}

	for(int i =0; i < numtofit; i++){
		if(Dpriors[i][0] !=0 && Dpriors[i][1] != 0){
			if(Drms[i] != 0){
				Dpriors[i][0]=Dpriors[i][0]*linearParams[i];
				Dpriors[i][1]=Dpriors[i][1]*linearParams[i];
				//printf("DP: %i %g  %g\n",i,Dpriors[i][0], Dpriors[i][1]);
			}
			else if(Drms[i] == 0){
				Dpriors[i][0]=0;
				Dpriors[i][1]=0;
				//printf("DP: %i %g %g %g %g\n",i,fitsig*trms,Drms[i],Dpriors[i][0], Dpriors[i][1]);
			}

						
		}
		else {
			printf("DP: Prior %i is set to zero, marginalising?\n ", i);
		}

// 		printf("compare: %i %g %g \n",i,trms/Drms[i],linearParams[i]);
	}
}

void TNupdateDD(pulsar *psr,double val,double err,int pos, double &outval)
{
	
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
     outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}



void TNupdateBTJ(pulsar *psr,double val,double err,int pos,int arr, double &outval)
{
 if (pos==param_pb || pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot
     || pos==param_bpjph
     || pos==param_bpja1
     || pos==param_bpjec
     || pos==param_bpjom
     || pos==param_bpjpb
     )
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val*(SECDAY*365.25)*180.0/M_PI;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}



void TNupdateBT(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val*(SECDAY*365.25)*180.0/M_PI;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}



void TNupdateBTX(pulsar *psr,double val,double err,int pos,int k, double &outval)
{
  if (pos==param_fb)
    {
      outval = (val/powl(1.0e7,k+1));

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
	outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val*(SECDAY*365.25)*180.0/M_PI;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}



void TNupdateDDGR(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos==param_mtot)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;
 
    }
  else if (pos==param_pbdot || pos==param_xpbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val*(SECDAY*365.25)*180.0/M_PI;

    }
}



void TNupdateDDH(pulsar *psr,double val,double err,int pos, double &outval){
  if (pos==param_pb){
    outval = val/SECDAY;

  }else if( pos == param_a1 || pos == param_ecc || pos == param_t0 || 
            pos == param_sini || pos == param_m2 || pos == param_gamma || 
            pos == param_h3 || pos == param_stig ){
    outval = val;

  }else if (pos==param_om){
    outval = val*180.0/M_PI;

  }else if (pos==param_pbdot){
    outval = val;

  }else if (pos==param_a1dot){
    outval = val;

  }else if (pos==param_omdot){
    outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

  }
}



void TNupdateDDK(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}




void TNupdateDDS(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma || pos == param_shapmax)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}


void TNupdateELL1H(pulsar *psr,double val,double err,int pos, double &outval){
  if( pos == param_pb ){
    outval = val/SECDAY;

  }else if( pos == param_a1 || pos == param_eps1 || pos == param_eps2 ||
            pos == param_tasc || pos == param_sini || pos == param_m2 || 
            pos == param_eps1dot || pos == param_eps2dot || 
            pos == param_h3 || pos == param_h4 || pos == param_stig ){
    outval = val;

  }else if( pos == param_pbdot ){
    outval = val;

  }else if( pos == param_a1dot ){
    // JPWV 15 March 2010
    outval = val;

  }else if( pos == param_omdot ){
    // JPWV 15 March 2010
    outval = val;//*(SECDAY*365.25)*180.0/M_PI;

  }
}



void TNupdateELL1(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
  else if (pos==param_a1 || pos==param_eps1 || pos==param_eps2 || pos==param_tasc
	   || pos==param_sini || pos == param_m2 
       || pos==param_eps1dot || pos==param_eps2dot
       || pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_pbdot)
    {
     outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val*(SECDAY*365.25)*180.0/M_PI;

    }
}



void TNupdateMSS(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
    else if (pos==param_a1 || pos==param_ecc || pos==param_t0)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
     outval = val;

    }
  else if (pos==param_omdot)
    {
     outval = val*(SECDAY*365.25)*180.0/M_PI;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}


void TNupdateT2(pulsar *psr,double val,double err,int pos,int arr, double &outval){
  if (pos==param_pb || pos==param_a1 || pos==param_ecc || pos==param_t0 || 
      pos==param_sini || pos==param_m2 || pos == param_gamma || 
      pos==param_eps1 || pos==param_eps2 || pos==param_tasc ||
      pos == param_bpjph || pos==param_bpja1 || pos==param_bpjec || 
      pos==param_bpjom || pos == param_bpjpb || pos==param_shapmax){
    outval = val;

  }
  else if (pos==param_om || pos==param_kom || pos==param_kin)
    {
      outval = val*180.0/M_PI;

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot || pos == param_eps1dot || pos==param_eps2dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}


void TNupdateParameters(pulsar *psr,int p,double *val,double *error, double *outval)
{
  int i,j,k;
  logdbg("Updating parameters");
  outval[0] = val[0];

  j=1;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].paramSet[k]==1 && psr[p].param[i].fitFlag[k]==1 && (i!=param_start && i!=param_finish))
	    {
	      if (i==param_f) 
		{
		  if (k==0)
		    {

		        outval[j] = -val[j];	

		    }
		  else 
		    {
		      longdouble scale;
		      scale=1.0L;
		      if (k==2)      scale=1.0e9L;
		      else if (k>2 && k<10)  scale=1.0e18L;
		      else if (k>9) scale=1.0e23L;
		      
		      outval[j] =  - (psr[p].param[param_f].val[0]*(val[j]/pow(24.0*3600.0,k+1))/scale);

		    }
		}
	      else if (i==param_dm || i==param_px || i==param_fddc || i==param_fddi || i==param_dmassplanet || i==param_dmx)
		{
		  outval[j] = val[j];

		  // This is slow - should be a better approach
		  if (i==param_dm){
		    outval[j] = val[j];
		  }
		}
	      else if (i==param_dshk)
		{
		  outval[j]  = val[j];

		}
	      else if (i==param_pmrv)
		{
		  outval[j] = 10.0*val[j]*360.0*60.0*60.0/(2.0*M_PI);

		}
	      else if (i==param_glph) /* Glitch phase */
		{
		  outval[j] -= val[j];     

		}
	      else if (i==param_glf0d) /* Glitch */
		{
		  outval[j] -= val[j]; 

		}
	      else if (i==param_gltd) /* Glitch time delay */
		{
		  outval[j] -= val[j]; 

		}
	      else if (i==param_glf0) /* Glitch permanent pulse frequency increment */
		{
		  outval[j] -= val[j]; 

		}
	      else if (i==param_glf1) /* Glitch permanent pulse frequency deriv. increment */
		{
		  outval[j] -= val[j]; //*psr[p].param[param_f].val[0]; 

		}
	      else if (i==param_glf2) /* Glitch permanent pulse frequency second deriv. increment */
		{
		  outval[j] -= val[j]*1.0e-27; //*psr[p].param[param_f].val[0]; 

		}
 	      else if (i==param_telx || i==param_tely || i==param_telz)
		{
		  outval[j] -= val[j];

		}	      
	      else if (i==param_raj)
		{
		  char retstr[100];
		  outval[j] += val[j];

		  

		}
	      else if (i==param_decj)
		{
		  char retstr[100];
		 outval[j] += val[j];

		}
	      else if (i==param_pmra) /* Return in radian/sec */
		{
		  outval[j] += val[j]*180.0/M_PI*60.0*60.0*
		    1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val[0]);

		}
	      else if (i==param_pmdec) /* Return in radian/sec */
		{
		 outval[j] += val[j]*180.0/M_PI*60.0*60.0*1000.0*
		    SECDAY*365.25/24.0/3600.0;

		}	       
	      else if (i==param_wave_om) /* Whitening procedure using sinusoids */
		{
		  int k;
		  for (k=0;k<psr->nWhite;k++)
		    {
		      outval[j]  -= val[j]; 

		      outval[j] -= val[j]; 

		    }
		  if (psr->waveScale==2) // Ignore the non-frequency derivative terms
		    {
		      for (k=0;k<psr->nWhite;k++)
			{
			  printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
			  printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
			}
		      //		      j+=psr->nWhite*2; 
		    }
		  j--;
		}		  
	      else if (i==param_ifunc) 
		{
		  int k;
		  for (k=0;k<psr->ifuncN;k++)
		    {
		      outval[j] -= val[j];

		      j++;
		    }
		}		  
	      else if (i==param_clk_offs) 
		{
		  int k;
		  for (k=0;k<psr->clkOffsN;k++)
		    {
		      outval[j] -= val[j];

		      j++;
		    }
		}		  
	      else if (i==param_tel_dx) 
		{
		  int k;
		  printf("HERE --- hello\n");
		  if (psr[p].param[param_tel_dx].val[0] == -1)
		    {
		      printf("Updating with %g %g\n",val[j],error[j]);
		      outval[j] -= val[j]; //*SPEED_LIGHT/1000.0;

		    }
		  else if (psr[p].param[param_tel_dx].val[0] < 2)
		    {
		      for (k=0;k<psr->nTelDX;k++)
			{
			  outval[j] -= val[j];

			  j++;
			}
		    }
		  else
		    {
		      for (k=0;k<psr->nTelDX-1;k++)
			{
			 outval[j] -= val[j];

			  j++;
			}
		    }
		}		  
	      else if (i==param_tel_dy) 
		{
		  int k;
		  if (psr[p].param[param_tel_dy].val[0] == -1)
		    {
		     outval[j] -= val[j]; //*SPEED_LIGHT/1000.0;

		    }
		  else if (psr[p].param[param_tel_dy].val[0] < 2)
		    {
		      for (k=0;k<psr->nTelDY;k++)
			{
			  outval[j] -= val[j];

			  j++;
			}
		    }
		  else
		    {
		      for (k=0;k<psr->nTelDY-1;k++)
			{
			  outval[j] -= val[j];

			  j++;
			}
		    }

		      
		}		  
	      else if (i==param_tel_dz) 
		{
		  int k;
		  if (psr[p].param[param_tel_dz].val[0] == -1)
		    {
		    outval[j] -= val[j]; //*SPEED_LIGHT/1000.0;

		    }
		  else if (psr[p].param[param_tel_dz].val[0] < 2)
		    {
		      for (k=0;k<psr->nTelDZ;k++)
			{
			  outval[j] -= val[j];

			  j++;
			}
		    }
		  else
		    {
		      for (k=0;k<psr->nTelDZ-1;k++)
			{
			 outval[j] -= val[j];

			  j++;
			}
		    }
		}		  
	      else if (i==param_quad_ifunc_p) 
		{
		  int k;
		  for (k=0;k<psr->quad_ifuncN_p;k++)
		    {
		      outval[j] -= val[j];

		      j++;
		    }
		}		  
	      else if (i==param_quad_ifunc_c) 
		{
		  int k;
		  for (k=0;k<psr->quad_ifuncN_c;k++)
		    {
		      outval[j] -= val[j];

		      j++;
		    }
		}		  
	      else if (i==param_gwsingle)
		{
		  printf("%d GW SINGLE -- in here %g %g %g %g\n",j,val[j],error[j],val[j+1],error[j+1]);
		  outval[j] -= val[j];
		  j++;
		  outval[j] -= val[j];
		  j++;
		  outval[j] -= val[j];
		  j++;
		  outval[j] -= val[j];
		  printf("Now have: %g %g\n",psr[p].gwsrc_aplus_r,psr[p].gwsrc_across_r);
		}
	      else if (i==param_gwm_amp)
		{
		  outval[j] -= val[j];

		  j++;
		}
	      else if (i==param_dmmodel)
		{
		  for (k=0;k<psr[p].dmoffsDMnum;k++)
		    {
				  outval[j] += val[j];

				  j++;
			}
		 for (k=0;k<psr[p].dmoffsCMnum;k++){

				  outval[j] = val[j];

				  j++;
		    }
		  j--;
		}
	      else if (i==param_start)
		{
		}
	      else if (i==param_finish)
		{
		}
	      else if (strcmp(psr[p].binaryModel,"BT")==0)
		TNupdateBT(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"BTJ")==0)
		TNupdateBTJ(&psr[p],val[j],error[j],i,k,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"BTX")==0)
		TNupdateBTX(&psr[p],val[j],error[j],i,k,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"ELL1")==0)
		TNupdateELL1(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DD")==0)
		TNupdateDD(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDK")==0)
		TNupdateDDK(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDS")==0)
		TNupdateDDS(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDGR")==0)
		TNupdateDDGR(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"MSS")==0)
		TNupdateMSS(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"T2")==0)
		TNupdateT2(&psr[p],val[j],error[j],i,k,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDH")==0)
		TNupdateDDH(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"ELL1H")==0)
		TNupdateELL1H(&psr[p],val[j],error[j],i,outval[j]);
	      j++; /* Increment position in fit list */
	    }
	}
    }
  if (strcmp(psr[p].binaryModel,"DDGR")==0) 
    DDGRmodel(psr,0,0,-2);  /* Update GR parameters */	  
  
  logdbg("Updating jumps; nJumps = %d",psr[p].nJumps);
  /* Now check jumps */
  for (i=1;i<=psr[p].nJumps;i++)
    {
      logdbg("%d fitJump = %d",i,psr[p].fitJump[i]);
      if (psr[p].fitJump[i]==1)
	{
	  logdbg("%d Jump changes",i);
	  logdbg("value = %g",(double)val[j]);
	  logdbg("error = %g",(double)error[j]);
	  outval[j] += -val[j];

	  j++;
	}
      /*	      printf("Have jumps %g %g\n",(double)val[j],error[j][j]); */
    }
  logdbg("Complete updating parameters");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////AND NOW THE INVERSE OF ALL OF THESE///////////////////////////////////////////////////////////////



void TNIupdateDD(pulsar *psr,double val,double err,int pos, double &outval)
{
	
  if (pos==param_pb)
    {
      outval = val*SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
     outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}



void TNIupdateBTJ(pulsar *psr,double val,double err,int pos,int arr, double &outval)
{
 if (pos==param_pb || pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot
     || pos==param_bpjph
     || pos==param_bpja1
     || pos==param_bpjec
     || pos==param_bpjom
     || pos==param_bpjpb
     )
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val/((SECDAY*365.25)*180.0/M_PI);

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}



void TNIupdateBT(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val/((SECDAY*365.25)*180.0/M_PI);

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}



void TNIupdateBTX(pulsar *psr,double val,double err,int pos,int k, double &outval)
{
  if (pos==param_fb)
    {
      outval = (val*powl(1.0e7,k+1));

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
	outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val/((SECDAY*365.25)*180.0/M_PI);

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}



void TNIupdateDDGR(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val/SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos==param_mtot)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);
 
    }
  else if (pos==param_pbdot || pos==param_xpbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val/((SECDAY*365.25)*180.0/M_PI);

    }
}



void TNIupdateDDH(pulsar *psr,double val,double err,int pos, double &outval){
  if (pos==param_pb){
    outval = val/SECDAY;

  }else if( pos == param_a1 || pos == param_ecc || pos == param_t0 || 
            pos == param_sini || pos == param_m2 || pos == param_gamma || 
            pos == param_h3 || pos == param_stig ){
    outval = val;

  }else if (pos==param_om){
    outval = val/(180.0/M_PI);

  }else if (pos==param_pbdot){
    outval = val;

  }else if (pos==param_a1dot){
    outval = val;

  }else if (pos==param_omdot){
    outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

  }
}



void TNIupdateDDK(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val*SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}




void TNIupdateDDS(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val*SECDAY;

    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma || pos == param_shapmax)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}


void TNIupdateELL1H(pulsar *psr,double val,double err,int pos, double &outval){
  if( pos == param_pb ){
    outval = val*SECDAY;

  }else if( pos == param_a1 || pos == param_eps1 || pos == param_eps2 ||
            pos == param_tasc || pos == param_sini || pos == param_m2 || 
            pos == param_eps1dot || pos == param_eps2dot || 
            pos == param_h3 || pos == param_h4 || pos == param_stig ){
    outval = val;

  }else if( pos == param_pbdot ){
    outval = val;

  }else if( pos == param_a1dot ){
    // JPWV 15 March 2010
    outval = val;

  }else if( pos == param_omdot ){
    // JPWV 15 March 2010
    outval = val;//*(SECDAY*365.25)*180.0/M_PI;

  }
}



void TNIupdateELL1(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val*SECDAY;

    }
  else if (pos==param_a1 || pos==param_eps1 || pos==param_eps2 || pos==param_tasc
	   || pos==param_sini || pos == param_m2 
       || pos==param_eps1dot || pos==param_eps2dot
       || pos==param_a1dot)
    {
      outval = val;

    }
  else if (pos==param_pbdot)
    {
     outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val/((SECDAY*365.25)*180.0/M_PI);

    }
}



void TNIupdateMSS(pulsar *psr,double val,double err,int pos, double &outval)
{
  if (pos==param_pb)
    {
      outval = val*SECDAY;

    }
    else if (pos==param_a1 || pos==param_ecc || pos==param_t0)
    {
      outval = val;

    }
  else if (pos==param_om)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
     outval = val;

    }
  else if (pos==param_omdot)
    {
     outval = val/((SECDAY*365.25)*180.0/M_PI);

    }
  else if (pos==param_a1dot)
    {
      outval = val;

    }
}


void TNIupdateT2(pulsar *psr,double val,double err,int pos,int arr, double &outval){
  if (pos==param_pb || pos==param_a1 || pos==param_ecc || pos==param_t0 || 
      pos==param_sini || pos==param_m2 || pos == param_gamma || 
      pos==param_eps1 || pos==param_eps2 || pos==param_tasc ||
      pos == param_bpjph || pos==param_bpja1 || pos==param_bpjec || 
      pos==param_bpjom || pos == param_bpjpb || pos==param_shapmax){
    outval = val;

  }
  else if (pos==param_om || pos==param_kom || pos==param_kin)
    {
      outval = val/(180.0/M_PI);

    }
  else if (pos==param_pbdot)
    {
      outval = val;

    }
  else if (pos==param_a1dot || pos == param_eps1dot || pos==param_eps2dot)
    {
      outval = val;

    }
  else if (pos==param_omdot)
    {
      outval = val; /* *(SECDAY*365.25)*180.0/M_PI; */

    }
}


void TNIupdateParameters(pulsar *psr,int p,double *val,double *error, double *outval)
{
  int i,j,k;
  logdbg("Updating parameters");
  outval[0] = val[0];

  j=1;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].paramSet[k]==1 && psr[p].param[i].fitFlag[k]==1 && (i!=param_start && i!=param_finish))
	    {
	      if (i==param_f) 
		{
		  if (k==0)
		    {

		        outval[j] = val[j];	

		    }
		  else 
		    {
		      longdouble scale;
		      scale=1.0L;
		      if (k==2)      scale=1.0e9L;
		      else if (k>2 && k<10)  scale=1.0e18L;
		      else if (k>9) scale=1.0e23L;
		      
		      outval[j] =   val[j]/(psr[p].param[param_f].val[0]/pow(24.0*3600.0,k+1)/scale);

		    }
		}
	      else if (i==param_dm || i==param_px || i==param_fddc || i==param_fddi || i==param_dmassplanet || i==param_dmx)
		{
		  outval[j] = val[j];

		  // This is slow - should be a better approach
		  if (i==param_dm){
		    outval[j] = val[j];
		  }
		}
	      else if (i==param_dshk)
		{
		  outval[j]  = val[j];

		}
	      else if (i==param_pmrv)
		{
		  outval[j] = val[j]/(10.0*360.0*60.0*60.0/(2.0*M_PI));

		}
	      else if (i==param_glph) /* Glitch phase */
		{
		  outval[j] = val[j];     

		}
	      else if (i==param_glf0d) /* Glitch */
		{
		  outval[j] = val[j]; 

		}
	      else if (i==param_gltd) /* Glitch time delay */
		{
		  outval[j] = val[j]; 

		}
	      else if (i==param_glf0) /* Glitch permanent pulse frequency increment */
		{
		  outval[j] = val[j]; 

		}
	      else if (i==param_glf1) /* Glitch permanent pulse frequency deriv. increment */
		{
		  outval[j] = val[j]; //*psr[p].param[param_f].val[0]; 

		}
	      else if (i==param_glf2) /* Glitch permanent pulse frequency second deriv. increment */
		{
		  outval[j] = val[j]/(1.0e-27); //*psr[p].param[param_f].val[0]; 

		}
 	      else if (i==param_telx || i==param_tely || i==param_telz)
		{
		  outval[j] = val[j];

		}	      
	      else if (i==param_raj)
		{

		  outval[j] = val[j];

		  

		}
	      else if (i==param_decj)
		{

		 outval[j] = val[j];

		}
	      else if (i==param_pmra) /* Return in radian/sec */
		{
		  outval[j] = val[j]/(180.0/M_PI*60.0*60.0*1000.0*SECDAY*365.25/24.0/3600.0*cos(psr[p].param[param_decj].val[0]));

		}
	      else if (i==param_pmdec) /* Return in radian/sec */
		{
		 outval[j] = val[j]/(180.0/M_PI*60.0*60.0*1000.0*SECDAY*365.25/24.0/3600.0);

		}	       
	      else if (i==param_wave_om) /* Whitening procedure using sinusoids */
		{
		  int k;
		  for (k=0;k<psr->nWhite;k++)
		    {
		      outval[j]  = val[j]; 

		      outval[j] = val[j]; 

		    }
		  if (psr->waveScale==2) // Ignore the non-frequency derivative terms
		    {
		      for (k=0;k<psr->nWhite;k++)
			{
			  printf("Ignoring cos %g %g\n",val[j],error[j]); j++;
			  printf("Ignoring sin %g %g\n",val[j],error[j]); j++;
			}
		      //		      j+=psr->nWhite*2; 
		    }
		  j--;
		}		  
	      else if (i==param_ifunc) 
		{
		  int k;
		  for (k=0;k<psr->ifuncN;k++)
		    {
		      outval[j] = val[j];

		      j++;
		    }
		}		  
	      else if (i==param_clk_offs) 
		{
		  int k;
		  for (k=0;k<psr->clkOffsN;k++)
		    {
		      outval[j] = val[j];

		      j++;
		    }
		}		  
	      else if (i==param_tel_dx) 
		{
		  int k;
		  printf("HERE --- hello\n");
		  if (psr[p].param[param_tel_dx].val[0] == -1)
		    {
		      printf("Updating with %g %g\n",val[j],error[j]);
		      outval[j] = val[j]; //*SPEED_LIGHT/1000.0;

		    }
		  else if (psr[p].param[param_tel_dx].val[0] < 2)
		    {
		      for (k=0;k<psr->nTelDX;k++)
			{
			  outval[j] = val[j];

			  j++;
			}
		    }
		  else
		    {
		      for (k=0;k<psr->nTelDX-1;k++)
			{
			 outval[j] = val[j];

			  j++;
			}
		    }
		}		  
	      else if (i==param_tel_dy) 
		{
		  int k;
		  if (psr[p].param[param_tel_dy].val[0] == -1)
		    {
		     outval[j] = val[j]; //*SPEED_LIGHT/1000.0;

		    }
		  else if (psr[p].param[param_tel_dy].val[0] < 2)
		    {
		      for (k=0;k<psr->nTelDY;k++)
			{
			  outval[j] = val[j];

			  j++;
			}
		    }
		  else
		    {
		      for (k=0;k<psr->nTelDY-1;k++)
			{
			  outval[j] = val[j];

			  j++;
			}
		    }

		      
		}		  
	      else if (i==param_tel_dz) 
		{
		  int k;
		  if (psr[p].param[param_tel_dz].val[0] == -1)
		    {
		    outval[j] = val[j]; //*SPEED_LIGHT/1000.0;

		    }
		  else if (psr[p].param[param_tel_dz].val[0] < 2)
		    {
		      for (k=0;k<psr->nTelDZ;k++)
			{
			  outval[j] = val[j];

			  j++;
			}
		    }
		  else
		    {
		      for (k=0;k<psr->nTelDZ-1;k++)
			{
			 outval[j] = val[j];

			  j++;
			}
		    }
		}		  
	      else if (i==param_quad_ifunc_p) 
		{
		  int k;
		  for (k=0;k<psr->quad_ifuncN_p;k++)
		    {
		      outval[j] = val[j];

		      j++;
		    }
		}		  
	      else if (i==param_quad_ifunc_c) 
		{
		  int k;
		  for (k=0;k<psr->quad_ifuncN_c;k++)
		    {
		      outval[j] = val[j];

		      j++;
		    }
		}		  
	      else if (i==param_gwsingle)
		{
		  printf("%d GW SINGLE -- in here %g %g %g %g\n",j,val[j],error[j],val[j+1],error[j+1]);
		  outval[j] = val[j];
		  j++;
		  outval[j] = val[j];
		  j++;
		  outval[j] = val[j];
		  j++;
		  outval[j] = val[j];
		  printf("Now have: %g %g\n",psr[p].gwsrc_aplus_r,psr[p].gwsrc_across_r);
		}
	      else if (i==param_gwm_amp)
		{
		  outval[j] = val[j];

		  j++;
		}
	      else if (i==param_dmmodel)
		{
		  for (k=0;k<psr[p].dmoffsDMnum;k++)
		    {
				  outval[j] = val[j];

				  j++;
			}
		 for (k=0;k<psr[p].dmoffsCMnum;k++){

				  outval[j] = val[j];

				  j++;
		    }
		  j--;
		}
	      else if (i==param_start)
		{
		}
	      else if (i==param_finish)
		{
		}
	      else if (strcmp(psr[p].binaryModel,"BT")==0)
		TNIupdateBT(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"BTJ")==0)
		TNIupdateBTJ(&psr[p],val[j],error[j],i,k,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"BTX")==0)
		TNIupdateBTX(&psr[p],val[j],error[j],i,k,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"ELL1")==0)
		TNIupdateELL1(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DD")==0)
		TNIupdateDD(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDK")==0)
		TNIupdateDDK(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDS")==0)
		TNIupdateDDS(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDGR")==0)
		TNIupdateDDGR(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"MSS")==0)
		TNIupdateMSS(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"T2")==0)
		TNIupdateT2(&psr[p],val[j],error[j],i,k,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"DDH")==0)
		TNIupdateDDH(&psr[p],val[j],error[j],i,outval[j]);
	      else if (strcmp(psr[p].binaryModel,"ELL1H")==0)
		TNIupdateELL1H(&psr[p],val[j],error[j],i,outval[j]);
	      j++; /* Increment position in fit list */
	    }
	}
    }
  if (strcmp(psr[p].binaryModel,"DDGR")==0) 
    DDGRmodel(psr,0,0,-2);  /* Update GR parameters */	  
  
  logdbg("Updating jumps; nJumps = %d",psr[p].nJumps);
  /* Now check jumps */
  for (i=1;i<=psr[p].nJumps;i++)
    {
      logdbg("%d fitJump = %d",i,psr[p].fitJump[i]);
      if (psr[p].fitJump[i]==1)
	{
	  logdbg("%d Jump changes",i);
	  logdbg("value = %g",(double)val[j]);
	  logdbg("error = %g",(double)error[j]);
	  outval[j] = val[j];

	  j++;
	}
      /*	      printf("Have jumps %g %g\n",(double)val[j],error[j][j]); */
    }
  logdbg("Complete updating parameters");
}
