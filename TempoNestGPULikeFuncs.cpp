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

double *TNFactorialList = new double[21];
void *GPUglobalcontext;


extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise,  double **FNF, double *NFd, int N, int RF, int DMF, int F, int incRED, int incDM);
extern "C" void NewLRedMarginGPUWrapper_(void *context, double *TNDMVec, double *Freqs, double *ObsFreqs, double *powercoeff, double *resvec, double *BATvec, double *DMVec, double *Noise, int *SysFlags, int N, int RF,int DMF, int DMScatterCoeff, int GroupNoiseCoeff,  int D, int F, int T, int incRED, int incDM, int incBandNoise, int incGroupNoise, int numFitTiming, int numFitJumps, double *likevals, int incNGJitter, int numNGJitterEpochs);

void assignGPUcontext(void *context){
        GPUglobalcontext=context;
}



void store_factorial(){

	for(int k=0; k <=20; k++){
		TNFactorialList[k]=iter_factorial(k);
	}
}


void GPUothpl(int n,double x,double *pl){


        double a=2.0;
        double b=0.0;
        double c=1.0;
        double y0=1.0;
        double y1=2.0*x;
        pl[0]=1.0;
        pl[1]=2.0*x;



        for(int k=2;k<n;k++){

                double c=2.0*(k-1.0);
                double yn=(a*x+b)*y1-c*y0;
                pl[k]=yn;
                y0=y1;
                y1=yn;

        }



}


void LRedGPULogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	//printf("hereNM");
	clock_t startClock,endClock;

	double *EFAC;
	double *EQUAD;
	int pcount=0;

	int numfit=((MNStruct *)GPUglobalcontext)->numFitTiming + ((MNStruct *)GPUglobalcontext)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];

	for(int p=0;p<ndim;p++){

		Cube[p]=(((MNStruct *)GPUglobalcontext)->Dpriors[p][1]-((MNStruct *)GPUglobalcontext)->Dpriors[p][0])*Cube[p]+((MNStruct *)GPUglobalcontext)->Dpriors[p][0];
	}


	if(((MNStruct *)GPUglobalcontext)->doLinear==0){
	
		for(int p=0;p< ((MNStruct *)GPUglobalcontext)->numFitTiming + ((MNStruct *)GPUglobalcontext)->numFitJumps; p++){
			LDparams[p]=Cube[p]*(((MNStruct *)GPUglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GPUglobalcontext)->LDpriors[p][0]);
		}

		double phase=(double)LDparams[0];
		pcount++;
		for(int p=1;p<((MNStruct *)GPUglobalcontext)->numFitTiming;p++){
			((MNStruct *)GPUglobalcontext)->pulse->param[((MNStruct *)GPUglobalcontext)->TempoFitNums[p][0]].val[((MNStruct *)GPUglobalcontext)->TempoFitNums[p][1]] = LDparams[pcount];	
			pcount++;
		}
		for(int p=0;p<((MNStruct *)GPUglobalcontext)->numFitJumps;p++){
			((MNStruct *)GPUglobalcontext)->pulse->jumpVal[((MNStruct *)GPUglobalcontext)->TempoJumpNums[p]]= LDparams[pcount];
			pcount++;
		}
		
		if(((MNStruct *)GPUglobalcontext)->pulse->param[param_dmmodel].fitFlag[0] == 1){
			int DMnum=((MNStruct *)GPUglobalcontext)->pulse[0].dmoffsDMnum;
			for(int i =0; i < DMnum; i++){
				((MNStruct *)GPUglobalcontext)->pulse[0].dmoffsDM[i]=Cube[ndim-DMnum+i];
			}
		}
		
		
		fastformBatsAll(((MNStruct *)GPUglobalcontext)->pulse,((MNStruct *)GPUglobalcontext)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)GPUglobalcontext)->pulse,((MNStruct *)GPUglobalcontext)->numberpulsars,1);       /* Form residuals */
		
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			Resvec[o]=(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].residual+phase;
		}
	
	}
	else if(((MNStruct *)GPUglobalcontext)->doLinear==1){
	
		for(int p=0;p < numfit; p++){
			Fitparams[p]=Cube[p];
			pcount++;
		}
		
		double *Fitvec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];

		dgemv(((MNStruct *)GPUglobalcontext)->DMatrix,Fitparams,Fitvec,((MNStruct *)GPUglobalcontext)->pulse->nobs,numfit,'N');
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			Resvec[o]=((MNStruct *)GPUglobalcontext)->pulse->obsn[o].residual-Fitvec[o];
		}
		
		delete[] Fitvec;
	}

	if(((MNStruct *)GPUglobalcontext)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)GPUglobalcontext)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)GPUglobalcontext)->pulse->nobs; o1++){
				double time = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o1].bat;
				if(time > StepTime){
					Resvec[o1] += StepAmp;
				}
			}
		}
	}

	if(((MNStruct *)GPUglobalcontext)->numFitEFAC == 0){
		EFAC=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			EFAC[o]=1;
		}
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEFAC == 1){
		EFAC=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			EFAC[o]=Cube[pcount];
		}
		pcount++;
		
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEFAC > 1){
		EFAC=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int p=0;p< ((MNStruct *)GPUglobalcontext)->systemcount; p++){
			EFAC[p]=Cube[pcount];
			pcount++;
		}
	}				

	if(((MNStruct *)GPUglobalcontext)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
                for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
		}
		pcount++;
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEQUAD > 1){
                EQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
                for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
                        EQUAD[o]=pow(10.0,2*Cube[pcount]);
			pcount++;
                }
        }

  	
        int FitRedCoeff=2*(((MNStruct *)GPUglobalcontext)->numFitRedCoeff);
		int FitDMCoeff=2*(((MNStruct *)GPUglobalcontext)->numFitDMCoeff);

        int totCoeff=0;
        if(((MNStruct *)GPUglobalcontext)->incRED != 0)totCoeff+=FitRedCoeff;
        if(((MNStruct *)GPUglobalcontext)->incDM != 0)totCoeff+=FitDMCoeff;

        double *powercoeff=new double[totCoeff];
        for(int o=0;o<totCoeff; o++){
                powercoeff[o]=0;
        }

	double *WorkNoise=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];

	double tdet=0;
	double timelike=0;



	double *BATvec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];

	if(((MNStruct *)GPUglobalcontext)->whitemodel == 0){
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
				WorkNoise[o]=pow(((((MNStruct *)GPUglobalcontext)->pulse->obsn[o].toaErr)*pow(10.0,-6))*EFAC[((MNStruct *)GPUglobalcontext)->sysFlags[o]],2) + EQUAD[((MNStruct *)GPUglobalcontext)->sysFlags[o]];
				
				tdet=tdet+log(WorkNoise[o]);
				WorkNoise[o]=1.0/WorkNoise[o];
				timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	
				BATvec[o]=(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].bat;
	
	
		}
	}
	else if(((MNStruct *)GPUglobalcontext)->whitemodel == 1){
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			//	printf("Noise %i %g %g %g\n",m1,Noise[m1],EFAC[flagList[m1]],EQUAD);
				WorkNoise[o]=EFAC[((MNStruct *)GPUglobalcontext)->sysFlags[o]]*EFAC[((MNStruct *)GPUglobalcontext)->sysFlags[o]]*(pow(((((MNStruct *)GPUglobalcontext)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)GPUglobalcontext)->sysFlags[o]]);
				
				tdet=tdet+log(WorkNoise[o]);
				WorkNoise[o]=1.0/WorkNoise[o];
				timelike=timelike+pow(Resvec[o],2)*WorkNoise[o];
	
				BATvec[o]=(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].bat;
	
	
		}
	}

	double *NFd = new double[totCoeff];
        double **FNF=new double*[totCoeff];
        for(int i=0;i<totCoeff;i++){
                FNF[i]=new double[totCoeff];
        }


	double start,end;
	int go=0;
	for (int i=0;i<((MNStruct *)GPUglobalcontext)->pulse->nobs;i++)
	  {
	    if (((MNStruct *)GPUglobalcontext)->pulse->obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat)
		      start = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat;
		    if (end < (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat)
		      end = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat;
		  }
	      }
	  }

	double maxtspan=end-start;


	double *freqs = new double[totCoeff];

        double *DMVec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];
        double DMKappa = 2.410*pow(10.0,-16);
        int startpos=0;
        double freqdet=0;
	

        if(((MNStruct *)GPUglobalcontext)->incRED==2){
        
        	if(((MNStruct *)GPUglobalcontext)->incFloatRed == 0){
                for (int i=0; i<FitRedCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
            }
            else if(((MNStruct *)GPUglobalcontext)->incFloatRed >0){

                for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed ; i++){
                
                		int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[i]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

                        powercoeff[i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[i+FitRedCoeff/2]=powercoeff[i];
                        freqdet=freqdet+2*log(powercoeff[i]);
                        pcount++;
                }
                
                 for (int i=FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed; i<FitRedCoeff/2; i++){
                		//printf("Freq: %g \n", Cube[pcount]);
                		
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
                }
                
            }
                
                
                startpos=FitRedCoeff;

        }
        else if(((MNStruct *)GPUglobalcontext)->incRED==3){

			freqdet=0;
			
			for(int pl = 0; pl < ((MNStruct *)GPUglobalcontext)->numFitRedPL; pl ++){
			
                double redamp=Cube[pcount];
                pcount++;
                double redindex=Cube[pcount];
                pcount++;

                
   			   redamp=pow(10.0, redamp);

				for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed ; i++){

		                freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[i]/maxtspan;
		                freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
		                
		                double PLcomp=redamp*redamp*pow((freqs[i]*365.25),-1.0*redindex)/(maxtspan*24*60*60);

		                powercoeff[i]+= PLcomp;
		                powercoeff[i+FitRedCoeff/2]+= PLcomp;
				}
			}
				
				
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed ; i++){
				freqdet=freqdet+2*log(powercoeff[i]);
		//		printf("%i %g %g \n",i,powercoeff[i], freqdet);
			}

                 for (int i=FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed; i<FitRedCoeff/2; i++){
                                
                    //    Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitRedCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                startpos=FitRedCoeff;

        }
// 		printf("DM\n");
        double nlist[((MNStruct *)GPUglobalcontext)->incFloatDM][2];
	if(((MNStruct *)GPUglobalcontext)->incDM==2){
        
        	if(((MNStruct *)GPUglobalcontext)->incFloatDM == 0){

                for (int i=0; i<FitDMCoeff/2; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=((MNStruct *)GPUglobalcontext)->sampleFreq[startpos/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed+i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }
            }
            else if(((MNStruct *)GPUglobalcontext)->incFloatDM >0){

                for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatDM ; i++){
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        freqs[startpos+i]=((MNStruct *)GPUglobalcontext)->sampleFreq[startpos/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed +i]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                        pcount++;
                }
                
                for (int i=FitDMCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatDM; i<FitDMCoeff/2; i++){

                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];

                         pcount++;
                         
                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;
                        
						
                        
                        

                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);
                       
                }




                
            }            

                for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].freqSSB,2));
                }


        }
        else if(((MNStruct *)GPUglobalcontext)->incDM==3){
        
        
				for(int pl = 0; pl < ((MNStruct *)GPUglobalcontext)->numFitDMPL; pl ++){
			
		            double DMamp=Cube[pcount];
		            pcount++;
		            double DMindex=Cube[pcount];
		            pcount++;

		            DMamp=pow(10.0, DMamp);

					for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatDM ; i++){

				            freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[startpos/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed +i]/maxtspan;
		                    freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
		                    
		                     double PLcomp=DMamp*DMamp*pow((freqs[startpos+i]*365.25),-1.0*DMindex)/(maxtspan*24*60*60);
				            
				            powercoeff[startpos+i]+=PLcomp;
		                    powercoeff[startpos+i+FitDMCoeff/2]+=PLcomp;
					}
				}
			
				for (int i=0; i<FitDMCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatDM ; i++){
					freqdet=freqdet+2*log(powercoeff[startpos+i]);
//					printf("%i %g %g \n", i, powercoeff[startpos+i], freqdet);
				}
			



                 for (int i= FitDMCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatDM ; i<FitDMCoeff/2; i++){

						//Cube[pcount]=floor(Cube[pcount]);
                        freqs[startpos+i]=Cube[pcount]/maxtspan;
                        freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
                        pcount++;

                        int pnum=pcount;
                        double pc=Cube[pcount];
                        pcount++;


                        powercoeff[startpos+i]=pow(10.0,pc)/(maxtspan*24*60*60);
                        powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
                        freqdet=freqdet+2*log(powercoeff[startpos+i]);

                }

                for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
                        DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].freqSSB,2));
                }

        }
	LRedGPUWrapper_(freqs, Resvec, BATvec, DMVec, WorkNoise, FNF, NFd, ((MNStruct *)GPUglobalcontext)->pulse->nobs, FitRedCoeff, FitDMCoeff, totCoeff,((MNStruct *)GPUglobalcontext)->incRED, ((MNStruct *)GPUglobalcontext)->incDM);



        double **PPFM=new double*[totCoeff];
        for(int i=0;i<totCoeff;i++){
                PPFM[i]=new double[totCoeff];
                for(int j=0;j<totCoeff;j++){
                        PPFM[i][j]=0;
                }
        }

        for(int c1=0; c1<totCoeff; c1++){
                PPFM[c1][c1]=1.0/powercoeff[c1];
        }



        for(int j=0;j<totCoeff;j++){
                for(int k=0;k<totCoeff;k++){
                        PPFM[j][k]=PPFM[j][k]+FNF[j][k];
                }
        }
        
        double jointdet=0;
        double freqlike=0;
       double *WorkCoeff = new double[totCoeff];
       for(int o1=0;o1<totCoeff; o1++){
                WorkCoeff[o1]=NFd[o1];
        }



        int info=0;
        dpotrfInfo(PPFM, totCoeff, jointdet, info);
        dpotrs(PPFM, WorkCoeff, totCoeff);
        for(int j=0;j<totCoeff;j++){
                freqlike += NFd[j]*WorkCoeff[j];
        }
        lnew=-0.5*(((double)((MNStruct *)GPUglobalcontext)->pulse->nobs)*log(2.0*M_PI) + tdet+jointdet+freqdet+timelike-freqlike);


	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,200);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);

		
	}
 	//printf("Like: %g %g %g %g %g %g\n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
	 //printf("CPULIKE: %g %g %g %g %g %g \n", lnew, jointdet,tdet,freqdet,timelike,freqlike);

	delete[] EFAC;
	delete[] EQUAD;
	delete[] WorkNoise;
	delete[] powercoeff;
	delete[] Resvec;
	delete[] BATvec;
	delete[] NFd;
	delete[] freqs;
	delete[] DMVec;
	delete[] WorkCoeff;

	for (int j = 0; j < totCoeff; j++){
		delete[] PPFM[j];
	}
	delete[] PPFM;



	for (int j = 0; j < totCoeff; j++){
		delete[] FNF[j];
	}
	delete[] FNF;


}



double NewLRedMarginGPULogLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context){
	
	clock_t startClock,endClock;

	double **EFAC;
	double *EQUAD;
	int pcount=0;
	
	
   	
	int numfit=((MNStruct *)GPUglobalcontext)->numFitTiming + ((MNStruct *)GPUglobalcontext)->numFitJumps;
	long double LDparams[numfit];
	double Fitparams[numfit];
	double *Resvec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	int fitcount=0;
	
	pcount=0;


        int TimetoMargin=0;
        for(int i =0; i < ((MNStruct *)GPUglobalcontext)->numFitTiming+((MNStruct *)GPUglobalcontext)->numFitJumps; i++){
                if(((MNStruct *)GPUglobalcontext)->LDpriors[i][2]==1)TimetoMargin++;
        }

	for(int p=0;p< ((MNStruct *)GPUglobalcontext)->numFitTiming + ((MNStruct *)GPUglobalcontext)->numFitJumps; p++){
		if(((MNStruct *)GPUglobalcontext)->Dpriors[p][1] != ((MNStruct *)GPUglobalcontext)->Dpriors[p][0]){

			LDparams[p]=Cube[fitcount]*(((MNStruct *)GPUglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GPUglobalcontext)->LDpriors[p][0]);
			fitcount++;
		}
		else if(((MNStruct *)GPUglobalcontext)->Dpriors[p][1] == ((MNStruct *)GPUglobalcontext)->Dpriors[p][0]){
			LDparams[p]=((MNStruct *)GPUglobalcontext)->Dpriors[p][0]*(((MNStruct *)GPUglobalcontext)->LDpriors[p][1]) + (((MNStruct *)GPUglobalcontext)->LDpriors[p][0]);
		}


	}
	pcount=0;
	double phase=(double)LDparams[0];
	pcount++;
	for(int p=1;p<((MNStruct *)GPUglobalcontext)->numFitTiming;p++){
		((MNStruct *)GPUglobalcontext)->pulse->param[((MNStruct *)GPUglobalcontext)->TempoFitNums[p][0]].val[((MNStruct *)GPUglobalcontext)->TempoFitNums[p][1]] = LDparams[pcount];	
		pcount++;
	}
	for(int p=0;p<((MNStruct *)GPUglobalcontext)->numFitJumps;p++){
		((MNStruct *)GPUglobalcontext)->pulse->jumpVal[((MNStruct *)GPUglobalcontext)->TempoJumpNums[p]]= LDparams[pcount];
		pcount++;
	}
	if(TimetoMargin != ((MNStruct *)GPUglobalcontext)->numFitTiming+((MNStruct *)GPUglobalcontext)->numFitJumps){	
		fastformBatsAll(((MNStruct *)GPUglobalcontext)->pulse,((MNStruct *)GPUglobalcontext)->numberpulsars);       /* Form Barycentric arrival times */
		formResiduals(((MNStruct *)GPUglobalcontext)->pulse,((MNStruct *)GPUglobalcontext)->numberpulsars,1);       /* Form residuals */
	}

	
	for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
		Resvec[o]=(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].residual+phase;

	}
	
	
	pcount=fitcount;

	if(((MNStruct *)GPUglobalcontext)->incStep > 0){
		for(int i = 0; i < ((MNStruct *)GPUglobalcontext)->incStep; i++){
			double StepAmp = Cube[pcount];
			pcount++;
			double StepTime = Cube[pcount];
			pcount++;
			for(int o1=0;o1<((MNStruct *)GPUglobalcontext)->pulse->nobs; o1++){
				if(((MNStruct *)GPUglobalcontext)->sysFlags[o1] == 3){
					if(((MNStruct *)GPUglobalcontext)->pulse->obsn[o1].bat > StepTime){
						Resvec[o1] += StepAmp;
					}
				}
			}
		}
	}	
	
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get White Noise vector///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  

	double uniformpriorterm=0;

	if(((MNStruct *)GPUglobalcontext)->numFitEFAC == 0){
		EFAC=new double*[((MNStruct *)GPUglobalcontext)->EPolTerms];
		for(int n=1; n <=((MNStruct *)GPUglobalcontext)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)GPUglobalcontext)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
					EFAC[n-1][o]=1;
				}
			}
			else{
                for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
                   EFAC[n-1][o]=0;
                }
			}
		}
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEFAC == 1){
		EFAC=new double*[((MNStruct *)GPUglobalcontext)->EPolTerms];
		for(int n=1; n <=((MNStruct *)GPUglobalcontext)->EPolTerms; n++){
			
			EFAC[n-1]=new double[((MNStruct *)GPUglobalcontext)->systemcount];
			if(n==1){
				for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
					
					EFAC[n-1][o]=pow(10.0,Cube[pcount]);
					if(((MNStruct *)GPUglobalcontext)->EFACPriorType ==1) {uniformpriorterm += log(EFAC[n-1][o]);}
				}
				pcount++;
			}
			else{
                                for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){

                                        EFAC[n-1][o]=pow(10.0,Cube[pcount]);
                                }
                                pcount++;
                        }
		}
		
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEFAC > 1){
		EFAC=new double*[((MNStruct *)GPUglobalcontext)->EPolTerms];
		for(int n=1; n <=((MNStruct *)GPUglobalcontext)->EPolTerms; n++){
			EFAC[n-1]=new double[((MNStruct *)GPUglobalcontext)->systemcount];
			if(n==1){
				for(int p=0;p< ((MNStruct *)GPUglobalcontext)->systemcount; p++){
					EFAC[n-1][p]=pow(10.0,Cube[pcount]);
					if(((MNStruct *)GPUglobalcontext)->EFACPriorType ==1) {uniformpriorterm += log(EFAC[n-1][p]);}
					pcount++;
				}
			}
			else{
                                for(int p=0;p< ((MNStruct *)GPUglobalcontext)->systemcount; p++){
                                        EFAC[n-1][p]=pow(10.0,Cube[pcount]);
                                        pcount++;
                                }
                        }
		}
	}	

		

	if(((MNStruct *)GPUglobalcontext)->numFitEQUAD == 0){
		EQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			EQUAD[o]=0;
		}
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEQUAD == 1){
		EQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
                for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			EQUAD[o]=pow(10.0,2*Cube[pcount]);
			if(((MNStruct *)GPUglobalcontext)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount]));}
		}
		pcount++;
	}
	else if(((MNStruct *)GPUglobalcontext)->numFitEQUAD > 1){
		EQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){

			if(((MNStruct *)GPUglobalcontext)->includeEQsys[o] == 1){
				//printf("Cube: %i %i %g \n", o, pcount, Cube[pcount]);
				EQUAD[o]=pow(10.0,2*Cube[pcount]);
				if(((MNStruct *)GPUglobalcontext)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount])); }
				pcount++;
			}
			else{
				EQUAD[o]=0;
			}
			//printf("Equad? %i %g \n", o, EQUAD[o]);
		}
    	}
    

    double *SQUAD;
	if(((MNStruct *)GPUglobalcontext)->incShannonJitter == 0){
		SQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			SQUAD[o]=0;
		}
	}
	else if(((MNStruct *)GPUglobalcontext)->incShannonJitter == 1){
		SQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
                for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
			SQUAD[o]=pow(10.0,2*Cube[pcount]);
			if(((MNStruct *)GPUglobalcontext)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount]));}
		}
		pcount++;
	}
	else if(((MNStruct *)GPUglobalcontext)->incShannonJitter > 1){
        SQUAD=new double[((MNStruct *)GPUglobalcontext)->systemcount];
        for(int o=0;o<((MNStruct *)GPUglobalcontext)->systemcount; o++){
            SQUAD[o]=pow(10.0,2*Cube[pcount]);
	    if(((MNStruct *)GPUglobalcontext)->EQUADPriorType ==1) { uniformpriorterm +=log(pow(10.0,Cube[pcount])); }

			pcount++;
        }
    }
	double *ECORRPrior;
	if(((MNStruct *)GPUglobalcontext)->incNGJitter >0){
		double *ECorrCoeffs=new double[((MNStruct *)GPUglobalcontext)->incNGJitter];	
		for(int i =0; i < ((MNStruct *)GPUglobalcontext)->incNGJitter; i++){
			ECorrCoeffs[i] = pow(10.0, 2*Cube[pcount]);
			pcount++;
		}
    		ECORRPrior = new double[((MNStruct *)GPUglobalcontext)->numNGJitterEpochs];
		for(int i =0; i < ((MNStruct *)GPUglobalcontext)->numNGJitterEpochs; i++){
			ECORRPrior[i] = ECorrCoeffs[((MNStruct *)GPUglobalcontext)->NGJitterSysFlags[i]];
		}

		delete[] ECorrCoeffs;
	} 
	double *Noise;	
	double *BATvec;
	Noise=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	BATvec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	
	
	for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
		BATvec[o]=(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].bat;
	}
		
		
	if(((MNStruct *)GPUglobalcontext)->whitemodel == 0){
	
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			double EFACterm=0;
			double noiseval=0;
			double ShannonJitterTerm=0;
			
			
			if(((MNStruct *)GPUglobalcontext)->useOriginalErrors==0){
				noiseval=((MNStruct *)GPUglobalcontext)->pulse->obsn[o].toaErr;
			}
			else if(((MNStruct *)GPUglobalcontext)->useOriginalErrors==1){
				noiseval=((MNStruct *)GPUglobalcontext)->pulse->obsn[o].origErr;
			}


			for(int n=1; n <=((MNStruct *)GPUglobalcontext)->EPolTerms; n++){
				EFACterm=EFACterm + pow((noiseval*pow(10.0,-6))/pow(pow(10.0,-7),n-1),n)*EFAC[n-1][((MNStruct *)GPUglobalcontext)->sysFlags[o]];
			}	
			
			if(((MNStruct *)GPUglobalcontext)->incShannonJitter > 0){	
			 	ShannonJitterTerm=SQUAD[((MNStruct *)GPUglobalcontext)->sysFlags[o]]*((MNStruct *)GPUglobalcontext)->TobsInfo[o]/1000.0;
			}

			Noise[o]= 1.0/(pow(EFACterm,2) + EQUAD[((MNStruct *)GPUglobalcontext)->sysFlags[o]]+ShannonJitterTerm);

		}
		
	}
	else if(((MNStruct *)GPUglobalcontext)->whitemodel == 1){
	
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){

			Noise[o]=1.0/(EFAC[0][((MNStruct *)GPUglobalcontext)->sysFlags[o]]*EFAC[0][((MNStruct *)GPUglobalcontext)->sysFlags[o]]*(pow(((((MNStruct *)GPUglobalcontext)->pulse->obsn[o].toaErr)*pow(10.0,-6)),2) + EQUAD[((MNStruct *)GPUglobalcontext)->sysFlags[o]]));
		}
		
	}
	

//////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////Form the Power Spectrum//////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////  

	int FitRedCoeff=2*(((MNStruct *)GPUglobalcontext)->numFitRedCoeff);
	int FitDMCoeff=2*(((MNStruct *)GPUglobalcontext)->numFitDMCoeff);
	int FitBandNoiseCoeff=2*(((MNStruct *)GPUglobalcontext)->numFitBandNoiseCoeff);
	int FitGroupNoiseCoeff = 2*((MNStruct *)GPUglobalcontext)->numFitGroupNoiseCoeff;

	if(((MNStruct *)GPUglobalcontext)->incFloatDM != 0)FitDMCoeff+=2*((MNStruct *)GPUglobalcontext)->incFloatDM;
	if(((MNStruct *)GPUglobalcontext)->incFloatRed != 0)FitRedCoeff+=2*((MNStruct *)GPUglobalcontext)->incFloatRed;


        int totCoeff=0;
        if(((MNStruct *)GPUglobalcontext)->incRED != 0)totCoeff+=FitRedCoeff;
        if(((MNStruct *)GPUglobalcontext)->incDM != 0)totCoeff+=FitDMCoeff;

	if(((MNStruct *)GPUglobalcontext)->incBandNoise > 0)totCoeff += ((MNStruct *)GPUglobalcontext)->incBandNoise*FitBandNoiseCoeff;


	 if(((MNStruct *)GPUglobalcontext)->incGroupNoise > 0)totCoeff+= ((MNStruct *)GPUglobalcontext)->incGroupNoise*FitGroupNoiseCoeff;


	if(((MNStruct *)GPUglobalcontext)->incNGJitter >0)totCoeff+=((MNStruct *)GPUglobalcontext)->numNGJitterEpochs;


	double *powercoeff=new double[totCoeff];
	for(int o=0;o<totCoeff; o++){
		powercoeff[o]=0;
	}

	double priorterm=0;
	bool uniformprior=0;
	double start,end;
	int go=0;
	for (int i=0;i<((MNStruct *)GPUglobalcontext)->pulse->nobs;i++)
	  {
	    if (((MNStruct *)GPUglobalcontext)->pulse->obsn[i].deleted==0)
	      {
		if (go==0)
		  {
		    go = 1;
		    start = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat;
		    end  = start;
		  }
		else
		  {
		    if (start > (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat)
		      start = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat;
		    if (end < (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat)
		      end = (double)((MNStruct *)GPUglobalcontext)->pulse->obsn[i].bat;
		  }
	      }
	  }

	double maxtspan=1*(end-start);


	double *freqs = new double[totCoeff];

	double *DMVec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	double *ObsFreqs = new double[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	int *GroupNoiseGroups = new int[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	double DMKappa = 2.410*pow(10.0,-16);
	int startpos=0;
	double freqdet=0;
	double GWBAmpPrior=0;


	if(((MNStruct *)GPUglobalcontext)->incRED==2){


		for (int i=0; i<FitRedCoeff/2; i++){
			int pnum=pcount;
			double pc=Cube[pcount];
			freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[i]/maxtspan;
			freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];

			powercoeff[i]=pow(10.0,2*pc);
			powercoeff[i+FitRedCoeff/2]=powercoeff[i];
			freqdet=freqdet+2*log(powercoeff[i]);
			pcount++;
		}
		    
		startpos=FitRedCoeff;

	}
	else if(((MNStruct *)GPUglobalcontext)->incRED==3 || ((MNStruct *)GPUglobalcontext)->incRED==4){

		freqdet=0;
		
		if(((MNStruct *)GPUglobalcontext)->FitLowFreqCutoff == 1){
			double fLow = Cube[pcount];
			pcount++;

			double deltaLogF = 0.1;
			double RedMidFreq = 2.0;

			double RedLogDiff = log10(RedMidFreq) - log10(fLow);
			int LogLowFreqs = floor(RedLogDiff/deltaLogF);

			double RedLogSampledDiff = LogLowFreqs*deltaLogF;
			double sampledFLow = floor(log10(fLow)/deltaLogF)*deltaLogF;
			
			int freqStartpoint = 0;

			//printf("Freq Info: %g %g %g %i %g \n", fLow, log10(fLow), sampledFLow,  LogLowFreqs, RedLogSampledDiff);

			for(int i =0; i < LogLowFreqs; i++){
				((MNStruct *)GPUglobalcontext)->sampleFreq[freqStartpoint]=pow(10.0, sampledFLow + i*RedLogSampledDiff/LogLowFreqs);
				freqStartpoint++;
			//	printf("%i %g %g \n", freqStartpoint-1, log10(((MNStruct *)GPUglobalcontext)->sampleFreq[freqStartpoint-1]), ((MNStruct *)GPUglobalcontext)->sampleFreq[freqStartpoint-1]);

			}

			for(int i =0;i < FitRedCoeff/2-LogLowFreqs; i++){
				((MNStruct *)GPUglobalcontext)->sampleFreq[freqStartpoint]=i+RedMidFreq;
				freqStartpoint++;
			//	printf("%i %g %g \n", freqStartpoint-1, log10(((MNStruct *)GPUglobalcontext)->sampleFreq[freqStartpoint-1]), ((MNStruct *)GPUglobalcontext)->sampleFreq[freqStartpoint-1]);
			}

		}


		
		for(int pl = 0; pl < ((MNStruct *)GPUglobalcontext)->numFitRedPL; pl ++){
			double redamp=Cube[pcount];
			pcount++;
			double redindex=Cube[pcount];
			pcount++;
	
   			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
    			
                        double cornerfreq=0;
                        if(((MNStruct *)GPUglobalcontext)->incRED==4){
                                cornerfreq=pow(10.0, Cube[pcount])/Tspan;
                                pcount++;
                        }
	
			redamp=pow(10.0, redamp);
			if(((MNStruct *)GPUglobalcontext)->RedPriorType ==1) { uniformpriorterm +=log(redamp); }


			double Agw=redamp;
	
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed ; i++){
	
				freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[i]/maxtspan;
				freqs[startpos+i+FitRedCoeff/2]=freqs[startpos+i];
				
                                double rho=0;
                                if(((MNStruct *)GPUglobalcontext)->incRED==3){
                                        rho = (Agw*Agw/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-redindex))/(maxtspan*24*60*60);
                                }
                                if(((MNStruct *)GPUglobalcontext)->incRED==4){

                        rho = pow((1+(pow((1.0/365.25)/cornerfreq,redindex/2))),2)*(Agw*Agw/12.0/(M_PI*M_PI))/pow((1+(pow(freqs[i]/cornerfreq,redindex/2))),2)/(maxtspan*24*60*60)*pow(f1yr,-3.0);
                                }

 				powercoeff[i]+= rho;
 				powercoeff[i+FitRedCoeff/2]+= rho;

			}
		}
		
		int coefftovary=0;
		double amptovary=0.0;
		if(((MNStruct *)GPUglobalcontext)->varyRedCoeff==1){
			coefftovary=int(pow(10.0,Cube[pcount]))-1;
			pcount++;
			amptovary=pow(10.0,Cube[pcount])/(maxtspan*24*60*60);
			pcount++;

			powercoeff[coefftovary]=amptovary;
			powercoeff[coefftovary+FitRedCoeff/2]=amptovary;	
		}		
		
		double GWBAmp=0;
		if(((MNStruct *)GPUglobalcontext)->incGWB==1){
			GWBAmp=pow(10.0,Cube[pcount]);
			pcount++;
			//GWBAmpPrior=log(GWBAmp);
			uniformpriorterm += log(GWBAmp);
			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
			for (int i=0; i<FitRedCoeff/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed ; i++){
				double rho = (GWBAmp*GWBAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[i]*365.25,(-4.333))/(maxtspan*24*60*60);	
				powercoeff[i]+= rho;
				powercoeff[i+FitRedCoeff/2]+= rho;
			}
		}
		for (int i=0; i<FitRedCoeff/2; i++){
			freqdet=freqdet+2*log(powercoeff[i]);
		}

        startpos=FitRedCoeff;

    }


	if(((MNStruct *)GPUglobalcontext)->incsinusoid == 1){
		double sineamp=pow(10.0,Cube[pcount]);
		pcount++;
		double sinephase=Cube[pcount];
		pcount++;
		double sinefreq=pow(10.0,Cube[pcount])/maxtspan;
		pcount++;		
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			Resvec[o]-= sineamp*sin(2*M_PI*sinefreq*(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].bat + sinephase);
		}
	}





       if(((MNStruct *)GPUglobalcontext)->incDM==2){

			for (int i=0; i<FitDMCoeff/2; i++){
				int pnum=pcount;
				double pc=Cube[pcount];
				freqs[startpos+i]=((MNStruct *)GPUglobalcontext)->sampleFreq[startpos/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed+i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
	
				powercoeff[startpos+i]=pow(10.0,2*pc);
				powercoeff[startpos+i+FitDMCoeff/2]=powercoeff[startpos+i];
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
				
				if(((MNStruct *)GPUglobalcontext)->DMPriorType ==1) { uniformpriorterm += log(powercoeff[startpos+i])/2.0; }
				
				pcount++;
			}
           	 


			for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
				DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].freqSSB,2));
			}

		startpos=startpos+FitDMCoeff;

        }
        else if(((MNStruct *)GPUglobalcontext)->incDM==3){

		for(int pl = 0; pl < ((MNStruct *)GPUglobalcontext)->numFitDMPL; pl ++){
			double DMamp=Cube[pcount];
			pcount++;
			double DMindex=Cube[pcount];
			pcount++;
			
   			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
    			

			DMamp=pow(10.0, DMamp);
			if(((MNStruct *)GPUglobalcontext)->DMPriorType ==1) { uniformpriorterm += log(DMamp); }
			for (int i=0; i<FitDMCoeff/2; i++){
	
				freqs[startpos+i]=(double)((MNStruct *)GPUglobalcontext)->sampleFreq[startpos/2 - ((MNStruct *)GPUglobalcontext)->incFloatRed +i]/maxtspan;
				freqs[startpos+i+FitDMCoeff/2]=freqs[startpos+i];
				
 				double rho = (DMamp*DMamp)*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-DMindex))/(maxtspan*24*60*60);	
 				powercoeff[startpos+i]+=rho;
 				powercoeff[startpos+i+FitDMCoeff/2]+=rho;
			}
		}
		
		
		int coefftovary=0;
		double amptovary=0.0;
		if(((MNStruct *)GPUglobalcontext)->varyDMCoeff==1){
			coefftovary=int(pow(10.0,Cube[pcount]))-1;
			pcount++;
			amptovary=pow(10.0,Cube[pcount])/(maxtspan*24*60*60);
			pcount++;

			powercoeff[startpos+coefftovary]=amptovary;
			powercoeff[startpos+coefftovary+FitDMCoeff/2]=amptovary;	
		}	
			
		
		for (int i=0; i<FitDMCoeff/2; i++){
			freqdet=freqdet+2*log(powercoeff[startpos+i]);
		}


		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			DMVec[o]=1.0/(DMKappa*pow((double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].freqSSB,2));
		}

	startpos=startpos+FitDMCoeff;

    }




/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Subtract  DM Shape Events////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  

	if(((MNStruct *)GPUglobalcontext)->incDMShapeEvent != 0){
                for(int i =0; i < ((MNStruct *)GPUglobalcontext)->incDMShapeEvent; i++){

			int numDMShapeCoeff=((MNStruct *)GPUglobalcontext)->numDMShapeCoeff;

                        double EventPos=Cube[pcount];
			pcount++;
                        double EventWidth=Cube[pcount];
			pcount++;


			double *DMshapecoeff=new double[numDMShapeCoeff];
			double *DMshapeVec=new double[numDMShapeCoeff];
			for(int c=0; c < numDMShapeCoeff; c++){
				DMshapecoeff[c]=Cube[pcount];
				pcount++;
			}


			double *DMshapeNorm=new double[numDMShapeCoeff];
			for(int c=0; c < numDMShapeCoeff; c++){
				DMshapeNorm[c]=1.0/sqrt(sqrt(2.0*M_PI)*pow(2.0,c)*iter_factorial(c));
			}

			for(int k=0;k<((MNStruct *)GPUglobalcontext)->pulse->nobs;k++){	
				double time=(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[k].bat;

				double HVal=(time-EventPos)/(sqrt(2.0)*EventWidth);
				GPUothpl(numDMShapeCoeff,HVal,DMshapeVec);
				double DMsignal=0;
				for(int c=0; c < numDMShapeCoeff; c++){
					DMsignal += DMshapeNorm[c]*DMshapeVec[c]*DMshapecoeff[c]*DMVec[k];
				}

				  Resvec[k] -= DMsignal*exp(-0.5*pow((time-EventPos)/EventWidth, 2));
	
			}




		delete[] DMshapecoeff;
		delete[] DMshapeVec;
		delete[] DMshapeNorm;

		}
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add Yearly DM///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

    
	if(((MNStruct *)GPUglobalcontext)->yearlyDM == 1){
		double yearlyamp=pow(10.0,Cube[pcount]);
		pcount++;
		double yearlyphase=Cube[pcount];
		pcount++;
		for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
			Resvec[o]-= yearlyamp*sin((2*M_PI/365.25)*(double)((MNStruct *)GPUglobalcontext)->pulse->obsn[o].bat + yearlyphase)*DMVec[o];
		}
	}








/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add Band Noise//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



        if(((MNStruct *)GPUglobalcontext)->incBandNoise > 0){

		for(int b = 0; b < ((MNStruct *)GPUglobalcontext)->incBandNoise; b++){


 			double startfreq = ((MNStruct *)GPUglobalcontext)->FitForBand[b][0];
                        double stopfreq = ((MNStruct *)GPUglobalcontext)->FitForBand[b][1];
                        double BandScale = ((MNStruct *)GPUglobalcontext)->FitForBand[b][2];
                        int BandPriorType = ((MNStruct *)GPUglobalcontext)->FitForBand[b][3];


                        double Bandamp=Cube[pcount];
                        pcount++;
                        double Bandindex=Cube[pcount];
                        pcount++;

                        double Tspan = maxtspan;
                        double f1yr = 1.0/3.16e7;

                        Bandamp=pow(10.0, Bandamp);
                        if(BandPriorType == 1) { uniformpriorterm += log(Bandamp); }

			

			Bandamp=pow(10.0, Bandamp);
			for (int i=0; i<FitBandNoiseCoeff/2; i++){

				freqs[startpos+i]=((double)(i+1.0))/maxtspan;
				freqs[startpos+i+FitBandNoiseCoeff/2]=freqs[startpos+i];
				
				double rho = (Bandamp*Bandamp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-Bandindex))/(maxtspan*24*60*60);	
				powercoeff[startpos+i]+=rho;
				powercoeff[startpos+i+FitBandNoiseCoeff/2]+=rho;
			}
			
			
				
			for (int i=0; i<FitBandNoiseCoeff/2; i++){
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
			}

			 for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
				ObsFreqs[o]=((MNStruct *)GPUglobalcontext)->pulse->obsn[o].freq;
			}



			startpos=startpos+FitBandNoiseCoeff;
		}
    }



/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add Group Noise/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



        if(((MNStruct *)GPUglobalcontext)->incGroupNoise > 0){

		for(int i =0; i < ((MNStruct *)GPUglobalcontext)->pulse->nobs; i++){
			GroupNoiseGroups[i] = 0;
		}


		for(int g = 0; g < ((MNStruct *)GPUglobalcontext)->incGroupNoise; g++){

			int GrouptoFit=0;
			if(((MNStruct *)GPUglobalcontext)->FitForGroup[g][0] == -1){
				GrouptoFit = floor(Cube[pcount]);
				pcount++;
				//printf("Fitting nois to group %i \n", GrouptoFit);
			}
			else{
				GrouptoFit = ((MNStruct *)GPUglobalcontext)->FitForGroup[g][0];
				
			}		

			double GroupAmp=Cube[pcount];
			pcount++;
			double GroupIndex=Cube[pcount];
			pcount++;
		
			double Tspan = maxtspan;
			double f1yr = 1.0/3.16e7;
		

			GroupAmp=pow(10.0, GroupAmp);
			if(((MNStruct *)GPUglobalcontext)->DMPriorType ==1) { uniformpriorterm += log(GroupAmp); }

			for (int i=0; i<FitGroupNoiseCoeff/2; i++){

				freqs[startpos+i]=((double)(i+1.0))/maxtspan;
				freqs[startpos+i+FitGroupNoiseCoeff/2]=freqs[startpos+i];
			
				double rho = (GroupAmp*GroupAmp/12.0/(M_PI*M_PI))*pow(f1yr,(-3)) * pow(freqs[startpos+i]*365.25,(-GroupIndex))/(maxtspan*24*60*60);	
				powercoeff[startpos+i]+=rho;
				powercoeff[startpos+i+FitGroupNoiseCoeff/2]+=rho;
			}
		
		
			
			for (int i=0; i<FitGroupNoiseCoeff/2; i++){
				freqdet=freqdet+2*log(powercoeff[startpos+i]);
			}

			for(int i =0; i < ((MNStruct *)GPUglobalcontext)->pulse->nobs; i++){
				
				if(((MNStruct *)GPUglobalcontext)->sysFlags[i] == GrouptoFit){
					//printf("Fitting noise to group %i %i %i \n", g, i, FitGroupNoiseCoeff);
					GroupNoiseGroups[i] = g+1;
				}
			}


			startpos=startpos+FitGroupNoiseCoeff;
		}

    }




/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add ECORR Coeffs////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	if(((MNStruct *)GPUglobalcontext)->incNGJitter >0){
		for(int i =0; i < ((MNStruct *)GPUglobalcontext)->numNGJitterEpochs; i++){
			powercoeff[startpos+i] = ECORRPrior[i];
			freqdet = freqdet + log(ECORRPrior[i]);
		}
	}




	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Get Time domain likelihood//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////  


	double tdet=0;
	double timelike=0;

	for(int o=0; o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){
		timelike+=Resvec[o]*Resvec[o]*Noise[o];
		tdet -= log(Noise[o]);
	}
    
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////get TNDMVec////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////   

    
	
	double **TNDM=new double*[((MNStruct *)GPUglobalcontext)->pulse->nobs];
	for(int i=0;i<((MNStruct *)GPUglobalcontext)->pulse->nobs;i++){
		TNDM[i]=new double[TimetoMargin];
	}
	double *TNDMVec=new double[((MNStruct *)GPUglobalcontext)->pulse->nobs*TimetoMargin];
	
	if(TimetoMargin != ((MNStruct *)GPUglobalcontext)->numFitTiming+((MNStruct *)GPUglobalcontext)->numFitJumps){
	
		getCustomDMatrixLike(GPUglobalcontext, TNDM);
	
		for(int g=0;g<TimetoMargin; g++){
			for(int o=0;o<((MNStruct *)GPUglobalcontext)->pulse->nobs; o++){

				TNDMVec[g*((MNStruct *)GPUglobalcontext)->pulse->nobs + o]=TNDM[o][g];
			}
		}
	}
	

    
    

//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////Call GPU Code/////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    
    
    
	double *likevals=new double[2];
	int totalsize=TimetoMargin+totCoeff;

	NewLRedMarginGPUWrapper_(GPUglobalcontext, TNDMVec, freqs, ObsFreqs, powercoeff, Resvec, BATvec, DMVec, Noise, GroupNoiseGroups,  ((MNStruct *)GPUglobalcontext)->pulse->nobs, FitRedCoeff,FitDMCoeff, FitBandNoiseCoeff, FitGroupNoiseCoeff, TimetoMargin, totCoeff, totalsize, ((MNStruct *)GPUglobalcontext)->incRED,((MNStruct *)GPUglobalcontext)->incDM, ((MNStruct *)GPUglobalcontext)->incBandNoise, ((MNStruct *)GPUglobalcontext)->incGroupNoise, ((MNStruct *)GPUglobalcontext)->numFitTiming, ((MNStruct *)GPUglobalcontext)->numFitJumps, likevals, ((MNStruct *)GPUglobalcontext)->incNGJitter, ((MNStruct *)GPUglobalcontext)->numNGJitterEpochs);
    
    
//////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////calculate likelihood///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    
	double jointdet=likevals[0];
	double freqlike=likevals[1];
    
	
	double lnew=-0.5*(tdet+jointdet+freqdet+timelike-freqlike) + uniformpriorterm;

	if(isnan(lnew) || isinf(lnew)){

		lnew=-pow(10.0,20);
// 		printf("red amp and alpha %g %g\n",redamp,redalpha);
// 		printf("Like: %g %g %g \n",lnew,Chisq,covdet);
		
	}
	
	delete[] ObsFreqs;
	delete[] GroupNoiseGroups;
	delete[] DMVec;
	delete[] EFAC;
	delete[] EQUAD;
	delete[] SQUAD;
	delete[] powercoeff;
	delete[] freqs;
	delete[] Noise;
	delete[] Resvec;
	delete[] BATvec;
	delete[] likevals;
	for (int j = 0; j < ((MNStruct *)GPUglobalcontext)->pulse->nobs; j++){
		delete[]TNDM[j];
	}
	delete[]TNDM;
	delete[]TNDMVec;

	if(((MNStruct *)GPUglobalcontext)->incNGJitter >0){
		delete[] ECORRPrior;
	}

//	printf("GPUChisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);


	return lnew;

// 	if(isinf(lnew) || isinf(jointdet) || isinf(tdet) || isinf(freqdet) || isinf(timelike) || isinf(freqlike)){
 	//printf("Chisq: %g %g %g %g %g %g \n",lnew,jointdet,tdet,freqdet,timelike,freqlike);
// 	}

}




