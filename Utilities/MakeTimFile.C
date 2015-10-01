#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <iterator>
#include <cstring>
#include <iostream>
#include <istream>
#include <fstream>


// psrchive stuff
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/FaradayRotation.h"
#include "Pulsar/PolnProfileStats.h"
#include "T2Observatory.h"
using namespace Pulsar;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////Program to Make Tim File from Archive/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc,char *argv[]){

	char *filename;
	char *TimName;
	char const *ObsCode;

	int providedTim=0;
	int providedObs=0;

  	// Check inputs
	for (int i=0;i<argc;i++){
		if (strcmp(argv[i],"-f")==0){
			filename=argv[i+1];
			printf("Opening File List: %s \n", filename);
		}
		else if (strcmp(argv[i],"-tim")==0){
			TimName=argv[i+1];
			printf("Writing to Tim File: %s \n", TimName);
			providedTim = 1;
		}
		else if (strcmp(argv[i],"-obscode")==0){
			ObsCode=argv[i+1];
			printf("Using Obs Code: %s \n", ObsCode);
			providedObs = 1;
		}
		else if (strcmp(argv[i],"-h")==0){
			printf("This program reads in a list of PSRChive files, and writes out a tim file that TempoNest can then use to do ToA estimation.\n");
			printf("It takes the following arguments:\n");
			printf("-h : This help info\n");
			printf("-f : A file list, with each archive to read in on a new line.\n");
			printf("-tim : The name of the output tim file (default new.tim)\n");
			printf("-obscode : The tempo2 observatory code for the observation (default obscode) \n");
			return 0;
			
		}
	}

	if(argc == 1){
		printf("No arguments provided: use -h option for help \n");
		return 0;
	}



	FILE *outputTim;
	if(providedTim == 0){
		outputTim = fopen("new.tim", "w");
		printf("Writing to Tim file: new.tim\n");
	}
	else{
		outputTim = fopen(TimName, "w");
	}

	if(providedObs == 0){
		ObsCode = "obscode";
	}

	fprintf(outputTim,"FORMAT 1\n");
	fprintf(outputTim,"MODE 1\n");

	int number_of_lines=0;
        std::ifstream checkfile;
        checkfile.open(filename);
        std::string line;
        while (getline(checkfile, line))
                ++number_of_lines;

	
        printf("number of Archives: %i \n",number_of_lines);


        checkfile.close();

	std::ifstream ProfileFiles;
	ProfileFiles.open(filename);

	for(int l = 0; l < number_of_lines; l++){

		std::string line;
		getline(ProfileFiles,line);
		printf("Opening File %i: %s \n", l, line.c_str());


		Reference::To< Archive > archive = Archive::load(line.c_str());
		if( !archive ){
			printf("Problem opening archive\n");			
			 return 0;	

		}

		int nsub = archive->get_nsubint();
		int nbins = 0;
		int nchans = 0;
		int npols = 0;
		double foldingperiod = 0;
		double inttime = 0;
		double centerfreq = 0;
		double chanfreq = 0;
		int intday = 0;
		double fracday = 0;
		int intsec = 0;
		double fracsecs = 0;
		bool isdedispersed = 0;

		if(nsub > 1){
			printf("Only one subint per archive supported currently: nsub %i \n", nsub);
			//return 0;
		}

		float *pvalues;
		for(int i=0; i < nsub; i++){
			Integration *subint = archive->get_Integration(i);

			nbins = subint->get_nbin();
			nchans = subint->get_nchan();
			npols = subint->get_npol();
			foldingperiod = subint->get_folding_period();
			inttime = subint->get_duration();
			centerfreq = subint->get_centre_frequency(); 
			MJD firstbin = subint->get_epoch();
			intday = firstbin.intday();
			fracday = firstbin.fracday();
			intsec = firstbin.get_secs();
			fracsecs = firstbin.get_fracsec();
			isdedispersed = subint->get_dedispersed();

			//if(nchans > 1){
			//	printf("Only one channel per archive supported currently: nchan %i\n", nchans);
				//return 0;
			//}


				
			//printf("archive info: %i %i %i %i %.16g %g %g %i %.16g %i %.16g\n", nsub, nbins, nchans, npols, foldingperiod, inttime, centerfreq, intday, fracday, intsec, fracsecs);
			//printf("dedispersed? %s \n", isdedispersed ? "true" : "false");
			double DMKappa = 2.410*pow(10.0,-16);
			double chanzerofreq = subint->get_centre_frequency(0);
			double DM = subint->get_dispersion_measure();
			double zerofreqdelay = DM/(DMKappa*pow((chanzerofreq)*pow(10.0,6), 2));
			double DMdelayDays=0;
			printf("DM? %g \n", DM);
			for(int j = 0; j < nchans; j++){
				int ipol = 0; //only want total Intensity
				Profile *prof = subint->get_Profile(ipol, j);
				pvalues = prof->get_amps();
				chanfreq = subint->get_centre_frequency(j);

				
				
				DMdelayDays = (DM/(DMKappa*pow((chanfreq)*pow(10.0,6), 2)) - zerofreqdelay)/24/60/60;
				double DMfracday = fracday + DMdelayDays;
				//printf("channel freq: %i %g %g\n", j, chanfreq, DMdelayDays*24*60*60);
				if(DMfracday < 0 || DMfracday >= 1){
					printf("DM delay has pushed day frac over a day.. fix this\n");
					return 0;
				}
			
		

				double Tobs = inttime;
				double ToAFreq = 0;
				if(isdedispersed){
					ToAFreq = centerfreq;
				}
				else{
					ToAFreq = chanfreq;
				}

				long double ProfileMJD = intday;   
				long double FirstBinSec = intsec+(long double)(fracsecs);  


				double oneflux = 0;	
				for(int b =0; b < nbins; b++){
					oneflux = oneflux + pvalues[b]; 
				}   

				char s[20];
				char s2[18];
				sprintf(s, "%.16g", DMfracday);

				int lastisaspace=0;
				for(int c =0; c < 18; c++){
					s2[c] = s[c+2];
					//printf("char: %i %c \n", i, s2[i]);
					if(s2[c]=='\0'){lastisaspace=1;}
				}
				//printf("%8s", (s[0] == '0' ? &s[1] : s));

				if(oneflux == 0){
					printf("Archive %i, Subint %i, Channel %i is all zeros, will not add to tim file \n", l, i, j);
				}
				else{

					printf("Archive details:\n");
					printf("Archive has %s been dedispersed, using %s frequency %g\n", isdedispersed ? "" : "not", isdedispersed ? "center" : "channel", ToAFreq);
					printf("SAT: %i.%s\n", intday, s2);
					printf("Tobs: %.5g\n", Tobs);
					fprintf(outputTim,"%s %.8f %i.%s 0.1 %s -tobs %.5g -nsub %i -nchan %i\n", line.c_str(), ToAFreq,  intday, s2, ObsCode, Tobs, nsub, nchans);


				}
			}

		}
	}

	fclose(outputTim);

	return 0;

}

