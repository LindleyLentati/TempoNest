double  Template2DProfLike(int &ndim, double *Cube, int &npars, double *DerivedParams, void *context){


	int debug = 0;
        int pcount = 0;


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Timing Model////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

	double FoldingPeriod = ((MNStruct *)globalcontext)->ReferencePeriod;

	double *phases = new double[((MNStruct *)globalcontext)->numProfComponents]();
	//printf("Phase: %g \n",Cube[pcount] );
	phases[0] = Cube[pcount]*((MNStruct *)globalcontext)->ReferencePeriod;
	pcount++;

	for(int i = 1; i < ((MNStruct *)globalcontext)->numProfComponents; i++){

		double compphase = phases[0] + Cube[pcount]*((MNStruct *)globalcontext)->ReferencePeriod;
		double WrappedCompPhase  = -FoldingPeriod/2 + fmod(FoldingPeriod + fmod(compphase + FoldingPeriod/2, FoldingPeriod), FoldingPeriod);

		phases[i] = WrappedCompPhase;
		pcount++;
	}





	double *betas = new double[((MNStruct *)globalcontext)->numProfComponents]();
	for(int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++){
		betas[i] = (pow(10.0, Cube[pcount]))*((MNStruct *)globalcontext)->ReferencePeriod;
		pcount++;
	}

	int totalProfCoeff = 0;
	int *numcoeff= new int[((MNStruct *)globalcontext)->numProfComponents];
	for(int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++){
		numcoeff[i] =  floor(pow(10.0, Cube[pcount]));
		totalProfCoeff += numcoeff[i];
		pcount++;
	}
	if(debug == 1){printf("Total coeff %i \n", totalProfCoeff);}

	double STau = 0;
	if(((MNStruct *)globalcontext)->incProfileScatter > 0){
		STau = pow(10.0, Cube[pcount]);	
		pcount++;
	}

/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profile Params//////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

		   
	int Nbins  = (int)((MNStruct *)globalcontext)->ProfileInfo[0][1];

	double *ProfileBats = new double[Nbins];
	double pulsesamplerate = FoldingPeriod/Nbins;

	for(int b =0; b < Nbins; b++){
		 ProfileBats[b] = pulsesamplerate*(b-Nbins/2);
	}


/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Profiles////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////




	double *Betatimes = new double[Nbins];
	double **Hermitepoly =  new double*[Nbins];
	for(int i =0;i<Nbins;i++){Hermitepoly[i]=new double[totalProfCoeff]();}



    	int cpos = 0;
	for(int i = 0; i < ((MNStruct *)globalcontext)->numProfComponents; i++){

		for(int j =0; j < Nbins; j++){
			Betatimes[j] = ProfileBats[j] - phases[i];
			Betatimes[j] = (-FoldingPeriod/2 + fmod(FoldingPeriod + fmod(Betatimes[j]+FoldingPeriod/2, FoldingPeriod), FoldingPeriod))/betas[i];
			
		}
		//printf("BT0: %g %g \n", Cube[0], Betatimes[0]);
		for(int j =0; j < Nbins; j++){

			TNothplMC(numcoeff[i],Betatimes[j],Hermitepoly[j], cpos);
				
		}


		for(int k =0; k < numcoeff[i]; k++){

			double Bconst=(1.0/sqrt(betas[i]/((MNStruct *)globalcontext)->ReferencePeriod))/sqrt(pow(2.0,k)*sqrt(M_PI));
			for(int j =0; j < Nbins; j++){
				Hermitepoly[j][cpos+k]=Hermitepoly[j][cpos+k]*Bconst*exp(-0.5*Betatimes[j]*Betatimes[j]);
			}

		}


		cpos += numcoeff[i];
   	 }


	///////////////////////////////////////////Marginalise over arbitrary offset and absolute amplitude////////////////////////////////////////////////////////////


	double **M = new double*[Nbins];



	int Msize = totalProfCoeff+1;


	if(debug ==1){printf("Made Basis Vectors %i \n", Msize);}

	for(int i =0; i < Nbins; i++){
		M[i] = new double[Msize];

		M[i][0] = 1;


		for(int j = 0; j < totalProfCoeff; j++){
			M[i][j+1] = Hermitepoly[i][j];
		}
	  
	}

	int badprof = 0;
/*        for(int j = 0; j < totalProfCoeff; j++){
                double mean = 0;
                for(int b = 0; b < Nbins; b++){
                        mean += M[b][j+1];

                }
                mean = mean/Nbins;

                double std = 0;
                for(int b = 0; b < Nbins; b++){
                        std += (M[b][j+1]-mean)*(M[b][j+1]-mean);

                }
                std = sqrt(std/Nbins);
        	if(fabs(1.0-std) > 0.05){badprof = 1;}
	}*/
	double **MNM = new double*[Msize];
	for(int i =0; i < Msize; i++){
	    MNM[i] = new double[Msize];
	}

	double **TempMNM = new double*[Msize];
	for(int i =0; i < Msize; i++){
	    TempMNM[i] = new double[Msize];
	}

	double lnew = 0;
	double Chisq=0;
	double detN = 0;

	if(((MNStruct *)globalcontext)->incProfileScatter == 0){



		dgemm(M, M , MNM, Nbins, Msize,Nbins, Msize, 'T', 'N');


                int csum=1;
                for(int comp = 0; comp <  ((MNStruct *)globalcontext)->numProfComponents; comp++){
                        double largestoff = 0;
                        for(int compi = 0; compi < numcoeff[comp]; compi++){
                                for(int compj = compi+1; compj < numcoeff[comp]; compj++){
                                        if(fabs(MNM[csum+compi][csum+compj]) > largestoff){
                                                largestoff = fabs(MNM[csum+compi][csum+compj]);
//                                                printf("LO: %i %i %g \n", compi, compj, largestoff);
                                        }
                                }
                        }
                        csum+=numcoeff[comp];
                        if(largestoff > pow(10.0, -10)){badprof = 1;}
                }

		double priorval = 1000.0;

			MNM[0][0] += 1.0/(0.01*0.01);
			for(int j = 1; j < Msize; j++){
				//printf("Prior: %i %g %g \n", j, MNM[j][j], 1.0/priorval);
				MNM[j][j] += 1.0/(priorval*priorval);

			}



		if(debug ==1){printf("Calculating Like \n");}

		for(int fc=0; fc < ((MNStruct *)globalcontext)->numTempFreqs; fc++){
	

			double *dNM = new double[Msize]();
			double *TempdNM = new double[Msize]();
			double *NDiffVec = new double[Nbins]();
			double *shapevec = new double[Nbins]();


			for(int j = 0; j < Nbins; j++){
				NDiffVec[j] = (double)((MNStruct *)globalcontext)->TemplateChans[fc*Nbins + j];
			}
			dgemv(M,NDiffVec,dNM,Nbins,Msize,'T');

			for(int i =0; i < Msize; i++){
				TempdNM[i] = dNM[i];
				for(int j =0; j < Msize; j++){
					TempMNM[i][j] = MNM[i][j];
				}
			}


			int info=0;
			double Margindet = 0;
			dpotrfInfo(TempMNM, Msize, Margindet, info);
			dpotrs(TempMNM, TempdNM, Msize);


			if(debug ==1){printf("Size: %i %g \n", Msize, Margindet);}

		
			dgemv(M,TempdNM,shapevec,Nbins,Msize,'N');




			double TNoise = 0;
			double OneChisq = 0;
		

			for(int j =0; j < Nbins; j++){
				double datadiff =  (double)((MNStruct *)globalcontext)->TemplateChans[fc*Nbins + j] - shapevec[j];
				OneChisq += datadiff*datadiff;
			}

			TNoise = sqrt(OneChisq/((double)Nbins));
			detN += Nbins*log(TNoise*TNoise);
			Chisq += OneChisq/TNoise/TNoise;


			delete[] shapevec;
			delete[] NDiffVec;
			delete[] dNM;
			delete[] TempdNM;

			if(debug ==1 || debug ==2){printf("End Like: %.10g %g %g %g\n", lnew, detN, Chisq, OneChisq);}
		}

	}
	else{


		double **SM = new double*[Nbins];
		for(int i =0; i < Nbins; i++){
			SM[i] = new double[Msize];

			SM[i][0] = 1;
		}


		for(int fc=0; fc < ((MNStruct *)globalcontext)->numTempFreqs; fc++){



			double TFreq = ((MNStruct *)globalcontext)->TemplateFreqs[fc];
			double ScatterScale = pow(TFreq*pow(10.0,6),4)/pow(10.0,9.0*4.0);
			double STime = STau/ScatterScale;

			//printf("Scatter: %i %g %g %g \n", fc, TFreq, ScatterScale, STime);


			double *OneBasis = new double[Nbins]();

			fftw_plan plan;
			fftw_plan planback;
			fftw_complex *FFTBasis;



			FFTBasis = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nbins/2+1));

			plan = fftw_plan_dft_r2c_1d(Nbins, OneBasis, FFTBasis, FFTW_ESTIMATE);
			planback = fftw_plan_dft_c2r_1d(Nbins,FFTBasis,  OneBasis,  FFTW_ESTIMATE);


			for(int c = 0; c < totalProfCoeff; c++){


				for(int b = 0; b < Nbins; b++){
					OneBasis[b] = M[b][c+1]; 
				}

			
				fftw_execute(plan);

	
				for(int b = 0; b < Nbins/2 + 1; b++){

					double f = (b*1.0)/((MNStruct *)globalcontext)->ReferencePeriod;
					double w = 2.0*M_PI*f;
					double RScatter = 1.0/(w*w*STime*STime+1); 
					double IScatter = -w*STime/(w*w*STime*STime+1);

					//if(c==0){printf("Scattered Basis: %i %g %g %g %g \n", b,RScatter, IScatter, FFTBasis[b][0], FFTBasis[b][1]);}

					double RScattered = FFTBasis[b][0]*RScatter - FFTBasis[b][1]*IScatter;
					double IScattered = FFTBasis[b][1]*RScatter + FFTBasis[b][0]*IScatter;
					FFTBasis[b][0] = RScattered;
					FFTBasis[b][1] = IScattered;
				}

				fftw_execute(planback);

				for(int b = 0; b < Nbins; b++){
					//if(c==0){printf("Scattered Basis: %i %g %g \n", b, M[b][c+1], OneBasis[b]/1024);}
					 SM[b][c+1] = OneBasis[b]/1024; 
					
				}

				

			}

	  		fftw_destroy_plan ( plan );
			fftw_destroy_plan ( planback );

			fftw_free ( FFTBasis );

			dgemm(SM, SM , MNM, Nbins, Msize,Nbins, Msize, 'T', 'N');

			double priorval = 1000.0;

			MNM[0][0] += 1.0/(0.01*0.01);
			for(int j = 1; j < Msize; j++){
				//printf("Prior: %i %g %g \n", j, MNM[j][j], 1.0/priorval);
				MNM[j][j] += 1.0/(priorval*priorval);

			}


			double *dNM = new double[Msize]();
			double *TempdNM = new double[Msize]();
			double *NDiffVec = new double[Nbins]();
			double *shapevec = new double[Nbins]();


			for(int j = 0; j < Nbins; j++){
				NDiffVec[j] = (double)((MNStruct *)globalcontext)->TemplateChans[fc*Nbins + j];
			}
			dgemv(SM,NDiffVec,dNM,Nbins,Msize,'T');
			for(int i =0; i < Msize; i++){
				TempdNM[i] = dNM[i];
			}

			int info=0;
			double Margindet = 0;
			dpotrfInfo(MNM, Msize, Margindet, info);
			dpotrs(MNM, TempdNM, Msize);

	
			dgemv(SM,TempdNM,shapevec,Nbins,Msize,'N');
			//printf("A0: %g %g %g %g %g %g %g %g\n", TempdNM[0], TempdNM[1], dNM[0], dNM[1], MNM[0][0], MNM[0][1], MNM[1][0], MNM[1][1]);

			double TNoise = 0;
			double OneChisq = 0;
	

			for(int j =0; j < Nbins; j++){
				double datadiff =  (double)((MNStruct *)globalcontext)->TemplateChans[fc*Nbins + j] - shapevec[j];
				OneChisq += datadiff*datadiff;
			}

			TNoise = ((MNStruct *)globalcontext)->TemplateNoise[fc]; //sqrt(OneChisq/((double)Nbins));
			detN += Nbins*log(TNoise*TNoise);
			Chisq += OneChisq/TNoise/TNoise;

			
			

			delete[] shapevec;
			delete[] NDiffVec;
			delete[] dNM;
			delete[] TempdNM;

		}

		for (int j = 0; j < Nbins; j++){

		    delete[] SM[j];
		}

		delete[] SM;



	}


	lnew = -0.5*(Chisq+detN)-totalProfCoeff*((MNStruct *)globalcontext)->numTempFreqs;
	if(badprof == 1){lnew -= 10000.0;}


	delete[] Betatimes;

	for (int j = 0; j < Nbins; j++){
	    delete[] Hermitepoly[j];
	    delete[] M[j];
	}
	delete[] Hermitepoly;
	delete[] M;

	for (int j = 0; j < Msize; j++){
	    delete[] MNM[j];
		delete[] TempMNM[j];
	}
	delete[] MNM;
	delete[] TempMNM;
	
	
	delete[] ProfileBats;
	delete[] numcoeff;
	delete[] betas;
	delete[] phases;

	if(debug ==1 || debug ==2){printf("End Like: %.10g %g %g\n", lnew, detN, Chisq);}

	if(debug ==1){sleep(5);}

	return lnew;

}

