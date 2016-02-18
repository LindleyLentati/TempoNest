#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

 
#include "mt19937.h"
#include "test_mt19937.h"


int test_random_uni()
{
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	double s=0;
	double s2=0;
	double u;
	double mu;
	double var;
	const unsigned long nrand=10000000;
	unsigned long i;
	ellipsis_mt19937_rng* rand;
	
	rand=(ellipsis_mt19937_rng*)malloc(sizeof(ellipsis_mt19937_rng));

	
	init_by_array(rand,init, length);
	
	/* find the mean and the variance */
	for(i=0;i<nrand;++i)
	{
		u=genrand_uniform(rand);
		s+=u;
		s2+=u*u;
	}
	mu=s/(double)nrand;
	var=(s2-s*s/(double)nrand)/(double)(nrand-1);
	
	/* mean should be 0.5 and the variance should be 1/12 */
	if(fabs(mu - 0.5)>1e-3 || fabs(var-1./12.)>1e-3)
	{
		free(rand);
		return EXIT_FAILURE;
	}

	free(rand);
	return EXIT_SUCCESS;
}

int test_random_norm()
{
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	double s=0;
	double s2=0;
	double u;
	double mu;
	double var;
	const unsigned long nrand=10000000;
	unsigned long i;
	ellipsis_mt19937_rng* rand;
	
	rand=(ellipsis_mt19937_rng*)malloc(sizeof(ellipsis_mt19937_rng));
	
	init_by_array(rand,init, length);
	
	/* find the mean and the variance */
	for(i=0;i<nrand;++i)
	{
		u=gerand_gauss(rand);
		s+=u;
		s2+=u*u;
	}
	mu=s/(double)nrand;
	var=(s2-s*s/(double)nrand)/(double)(nrand-1);
	
	/* mean should be 0 and the variance should be 1 */
	if(fabs(mu - 0)>1e-3 || fabs(var-1.)>1e-3)
	{
		free(rand);
		return EXIT_FAILURE;
	}

	free(rand);
	return EXIT_SUCCESS;
}

int test_rand_save_state()
{
	const unsigned nrand=10;
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	double* rands;
	unsigned i;
	ellipsis_mt19937_rng* rand;
	rands=(double*)malloc(nrand*sizeof(double));
	
	rand=(ellipsis_mt19937_rng*)malloc(sizeof(ellipsis_mt19937_rng));
	
	init_by_array(rand,init, length);

	/* store the first half */
	for(i=0;i<nrand/2;++i)
	{
		rands[i]=genrand_uniform(rand);
	}

	/* save the state */
	save_rand_state(rand,"myrand.dat");	
	
	
	/* store the 2nd half */
	for(i=0;i<nrand/2;++i)
	{
		rands[nrand/2+i]=genrand_uniform(rand);
	}
	
	/* read the state */
	read_rand_state(rand,"myrand.dat");

	
	/* check if the 2nd half is the same */
	for(i=0;i<nrand/2;++i)
	{
		if(rands[nrand/2+i]!=genrand_uniform(rand))
			return EXIT_FAILURE;
	}

	
	free(rand);
	return EXIT_SUCCESS;
	
}

