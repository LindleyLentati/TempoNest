#include <string.h>
#include <stdio.h>
#include "configfile.h"



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


void setupMNparams(int &IS, int &modal, int &ceff, int &nlive, double &efr){


	//IS: flag to use importance sampling, will be supported by upcoming release of multinest, at which point it should be set to 1
	IS=0;
	
	//modal: flag to allow multinest to search for multiple modes in the data. 1 = multimodal, 0 = single mode
	modal=1;

	//ceff: flag to set multinest to constant efficiency mode.  Adjusts sampling to maintain the efficiency set by the efr parameter. This is usefull for large dimensional problems (> 20dim), however the accuracy of the evidence suffers if importance sampling isn't used.
	ceff=0;

	//nlive: the number of live points used by multinest.  THe more you have the more it explores the parameter space, but the longer it takes to sample.
	//As a guide depending on the dimensionality of the problem (this is the numer of parameters sampled after discounting those to be marginalised over): 
	//Less than 5 dimensions: 100
	//Between 5 and 10 dimensions: 200
	//Between 10 and 50: 500
	//More than 50: 1000
	nlive=500;

	//efr: The sampling efficiency.  The Lower this number the more carefully multinest explores the parameter space, so the longer it takes.
	//The value depends on the dimensionality and the goal.
	//For parameter estimation it can be set higher than if an accurate Evidence value is required.
	//As a rough guide:
	//Less than 10 dimensions: 0.8 (parameter estimation), 0.3 (Evidence evaluation)
	//Between 10 and 20 dimensions: 0.3 (parameter estimation), 0.1 (Evidence evaluation)
	//Between 20 and 50: 0.1 (parameter estimation), 0.05 (Evidence evaluation)
	//More than 50: 0.05 (parameter estimation), 0.01 (Evidence evaluation)
	//These numbers may well be adjusted with experience, for parameter estimation of a single modal "blob", these may be much lower than required.
	//Also: In constant efficiency mode this must be set lower, to 0.05 for D<50 and 0.01 D>50.
	efr=0.1;



    // Use a configfile, if we can, to overwrite the defaults set in this file.
         try {
                 string strBuf;
                 strBuf = string("defaultparameters.conf");
                 ConfigFile parameters(strBuf);

		parameters.readInto(IS, "IS", IS);
		parameters.readInto(modal, "modal", modal);
		parameters.readInto(ceff, "ceff", ceff);
		parameters.readInto(nlive, "nlive", nlive);
		parameters.readInto(efr, "efr", efr);
 
	    } catch(ConfigFile::file_not_found oError) {
		printf("WARNING: parameters file '%s' not found. Using defaults.\n", oError.filename.c_str());
	    } // try

}
