#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "test_kinetic_energy.h"
#include "guided_hmc.h"

int test_kinetic_energy()
{
	const unsigned long ndim=100000;
	double * mtm;
	unsigned long i;
	double kin_eng,kin_eng_test;
	
	mtm=(double*)malloc(ndim*sizeof(double));
	
	kin_eng=0;
	for(i=0;i<ndim;++i)
	{
		mtm[i]=1.;
		kin_eng+=mtm[i]*mtm[i];
	}
	kin_eng*=0.5;
	
	kinetic_energy(ndim,mtm,&kin_eng_test);
	
	if(kin_eng != kin_eng_test)
	{
		free(mtm);
		return EXIT_FAILURE;
	}
	
	
	free(mtm);
	return EXIT_SUCCESS;
	
}


