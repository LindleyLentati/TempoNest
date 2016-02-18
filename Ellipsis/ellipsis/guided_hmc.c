#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
 
#include "guided_hmc.h"
#include "mt19937.h"
#include "hanson.h"

/*Interface for calling from FORTRAN gnu, intel */
void run_guided_hmc_(int* num_dim,double* start_point,
	double* dim_scale_fact,int* max_steps,double* step_sizes,
	char* file_prefix_f,int* seed,int* resume,int* feedback_int,
	void (*neg_logpost)(int*,double*,double*,double*),
	void (*write_extract)(int*,double*,double*,double*) ,
	int* nburn,int* nsamp, int* doMaxLike)
{

	char* file_prefix_c;
	char * pch;
	int length;
	
	/* fortran string to c strig */
	pch=strchr(file_prefix_f,' ');
	length=(int) (pch-file_prefix_f);
	file_prefix_c=(char*)malloc((length+32)*sizeof(char));
	memcpy(file_prefix_c,file_prefix_f,length*sizeof(char));

	/* call the c interface */
	run_guided_hmc(*num_dim,start_point,*dim_scale_fact,*max_steps,step_sizes,
		file_prefix_c,*seed,*resume,*feedback_int,neg_logpost,write_extract,
		*nburn,*nsamp, *doMaxLike);

	
	free(file_prefix_c);	

	
}

/*Interface for calling from FORTRAN openf90 */
void run_guided_hmc__(int* num_dim,double* start_point,
	double* dim_scale_fact,int* max_steps,double* step_sizes,
	char* file_prefix_f,int* seed,int* resume,int* feedback_int,
	void (*neg_logpost)(int*,double*,double*,double*),
	void (*write_extract)(int*,double*,double*,double*) ,
	int* nburn,int* nsamp, int* doMaxLike)
{

	char* file_prefix_c;
	char * pch;
	int length;
	
	/* fortran string to c strig */
	pch=strchr(file_prefix_f,' ');
	length=(int) (pch-file_prefix_f);
	file_prefix_c=(char*)malloc((length+32)*sizeof(char));
	memcpy(file_prefix_c,file_prefix_f,length*sizeof(char));

	/* call the c interface */
	run_guided_hmc(*num_dim,start_point,*dim_scale_fact,*max_steps,step_sizes,
		file_prefix_c,*seed,*resume,*feedback_int,neg_logpost,write_extract,
		*nburn,*nsamp, *doMaxLike);

	
	free(file_prefix_c);	

	
}


/*Interface for C/C++ */
void run_guided_hmc(int num_dim,double* start_point,
	double dim_scale_fact,int max_steps,double* step_sizes,
	char* file_prefix,int seed,int resume,int feedback_int,
	void (*neg_logpost)(int*,double*,double*,double*),
	void (*write_extract)(int*,double*,double*,double*) ,
	int nburn,int nsamp, int doMaxLike)
{

	int resume_file_update_int;	
	char diag_file_name[128];
	char rand_file_name[128];
	FILE* diag_out_file;
	FILE* rand_out_file;
	int run_engine=1,run_ghs=1;
	int count;
	int iteration;
	int count_last_set;
	int iteration_last_set;
	double* momentum;
	double* proposal;
	double* grad;
	int i;
	ellipsis_mt19937_rng* rand;
	double uni_rand;
	int num_steps;
	double log_uni;
	double epsilon;
	double log_ratio;
	hanson_data* hdata;
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	
	printf("\n-------------------------------------------------\n");
	printf("          Guided Hamiltonian Sampler\n");
	printf("                Version %d.%d\n",GHS_VERSION_MAJOR,GHS_VERSION_MINOR);
	printf("  S. T. Balan, M. A. J. Ashdown & M. P. Hobson\n");
	printf("     Cavendish Laboratory, Cambridge, UK.\n");
	printf("-------------------------------------------------\n");


	if(doMaxLike == 1){
		printf("Running in Maximum Likelihood mode\n");
		max_steps = 1;
	}
	
	/* allocate memroy for random number generator and hanson's diagnostic*/
	rand=(ellipsis_mt19937_rng*)malloc(sizeof(ellipsis_mt19937_rng));
	hdata=(hanson_data*)malloc(sizeof(hanson_data));
	
	/* init Hanson data */
	init_hanson_diag_data(hdata,num_dim,0.8,1.2);
	init_by_array(rand,init, length);
	
	/* assign the file names */
	strcpy(diag_file_name,file_prefix);
	strcat(diag_file_name,".diag.dat");
	strcpy(rand_file_name,file_prefix);
	strcat(rand_file_name,".rand.dat");
	
	/* are we resuming from pervious state? */
	if(resume)
	{
		printf("\nResuming from previous state...\n");
		
		/* check if resume files exist files for output */
		diag_out_file=fopen(diag_file_name,"r");
		rand_out_file=fopen(rand_file_name,"r");
		
		if(diag_out_file == NULL || rand_out_file == NULL)
		{
			printf("No resume files foud...\n");
			printf("Starting sampling from scratch...\n");
			
			/* see if we can open files for writing */
			diag_out_file=fopen(diag_file_name,"w");
			rand_out_file=fopen(rand_file_name,"w");
			if(diag_out_file == NULL || rand_out_file == NULL)
			{
				printf("\nCannot open the diagnostic files for writing!");
				printf("\nPlease check the file_prefix.\n\n");
				run_ghs=0;
			}
			else
			{
				/* initialise the random number generator */
				init_genrand(rand,(unsigned long)seed);				
			}		
		}
		else
		{
			/* close the opened files*/
			fclose(diag_out_file);
			fclose(rand_out_file);

			/* read the diagnostic data */
			printf("Reading diagnostic data...\t");
			read_diag_data_from_file(hdata,diag_file_name);
			printf("Done!\n");

			/* has the algorithm converged already? */
			if(!hdata->keep_sampling)
			{
				printf("Algorithm has already converged!\n");
				run_ghs=0;
			}
			else
			{
		
				/* read the random number data */
				printf("Reading the random number state...\t");
				read_rand_state(rand,rand_file_name);
				printf("Done!\n");
		
				/* set the start point */
				for(i=0;i<num_dim;++i)
				{
					start_point[i]=hdata->start_point[i];
				}
			
			}		
		}		
	}
	else
	{
		printf("Starting sampling from scratch...\n");
		init_genrand(rand,seed);
		run_ghs=1;
	}	
	
		
	/* update interval for the diagnostic files */
	if(num_dim > 1000)
	{
		resume_file_update_int=100;
	}
	else
	{
		resume_file_update_int=500;
	}
	
	/* run guided hmc */
	if(run_ghs)
	{
		printf("\nNumber of dimensions in the posterior   = %u\n",num_dim);
		printf("Number of steps in the HMC              = %u\n",max_steps);
		printf("Dimensionality scale factor             = %f\n",dim_scale_fact);
		
		run_engine=1;
		count=hdata->num_ents;/*if resuming count should be equal to num_ents*/
		iteration=0;
		
		count_last_set=0;
		iteration_last_set=0;
		
		/* allocate memory*/
		momentum=(double*)malloc(num_dim*sizeof(double));
		proposal=(double*)malloc(num_dim*sizeof(double));
		grad=(double*)malloc(num_dim*sizeof(double));
		
		while(run_engine)
		{
			++iteration;
			++iteration_last_set;
			/* draw a momentum sample from N(0,1) */
			for(i=0;i<num_dim;++i)
			{
				momentum[i]=gerand_gauss(rand);
				proposal[i]=start_point[i];
			}
			
			/* randomise the trajectory */
			uni_rand=genrand_uniform(rand);
			num_steps=1+(int)( (1.-uni_rand)*((double)max_steps-1.) );
			epsilon=dim_scale_fact*uni_rand;
			
			/* Metropolis-Hastings */
			log_uni=log(genrand_uniform(rand));
			log_ratio=log_uni;
			
			/* evolve the sample in phase space */
			double log_lik=0;
			leapfrog(num_dim,proposal,grad,momentum,num_steps,
				epsilon,step_sizes,&log_ratio,neg_logpost,&log_lik);


			if(doMaxLike == 1){
				log_uni = 0;
			}
				
			/* accept / reject */
			if(-log_ratio > log_uni)
			{
				/*burn the first nburn samples*/
				if(count >= nburn)
				{
					/* calculate diagnostic */
					push_state(hdata,num_dim,proposal,grad);
				
					/* write extract */
					write_extract(&num_dim,proposal,&log_lik,grad);
				}
				
				for(i=0;i<num_dim;++i)
				{
					start_point[i]=proposal[i];
				}
				++count;
				++count_last_set;
				
				/* write feedback to the console */
				if(feedback_int >0 && count>0 
					&& count%feedback_int ==0)
				{
					/* if we have not passed the nburn, change the message*/
					if(count <= nburn)
					{
						printf("\nNumber of samples BURNED so far         = %d\n",
							count);
						printf("Acceptance rate (last %d samples)     = %f\n\n",
							feedback_int,
							(double)count_last_set/(double)iteration_last_set);
					}
					else
					{
						printf("\nNumber of samples DRAWN so far          = %d\n",
							hdata->num_ents);
						printf("Acceptance rate (last %d samples)     = %f\n\n",
							feedback_int,
							(double)count_last_set/(double)iteration_last_set);
					}
						
					count_last_set=0;
					iteration_last_set=0;
				}
				
				/* update resume files */
				if(count%resume_file_update_int==0)
				{
					/* write the RNG state */
					save_rand_state(rand,rand_file_name);
					
					/* write the dianostic file */
					write_diag_data_to_file(hdata,diag_file_name);
					
					/* keep sampling?  only if we need to run till convergence*/
					if(nsamp==0)
					{
						run_engine=hdata->keep_sampling;
					}					
				}
				
				/* force the sampler to stop if count > nburn+nsamp*/
				if(nsamp>0 && count >= nburn && count >= (nburn+nsamp))
				{
					/* write the RNG state */
					save_rand_state(rand,rand_file_name);
					
					/* write the dianostic file */
					write_diag_data_to_file(hdata,diag_file_name);
					
					/* keep sampling? */
					run_engine=hdata->keep_sampling;
					
					run_engine=0;	
				}
			}
		}
		
		/* check if the sampler stopped after hanson's check */
		if(hdata->keep_sampling)
		{
			printf("\nSampling stopped!\n");
			printf("Toal number of samples taken            = %d\n",
				hdata->num_ents);
			printf("Net acceptance rate                     = %f\n",
				(double)count/(double)iteration);
		}
		else
		{
			printf("\nSampling finished!\n");
			printf("Toal number of samples taken            = %d\n",
				hdata->num_ents);
			printf("Net acceptance rate                     = %f\n",
				(double)count/(double)iteration);
		}
		
		free(grad);
		free(proposal);
		free(momentum);
		/*free_hanson_diag_data(hdata);*/
	}
	else
	{
		printf("No sampling performed\n\n");
	}
	
	free_hanson_diag_data(hdata);
	free(rand);
	free(hdata);
}

void leapfrog(int num_dim,double* x,double *g,double *p,
	int num_steps,double epsilon,double* sigma,double *log_ratio,
	void (*neg_logpost)(int*,double*,double*,double*),double* llik)
{
	double ke,pe,h0,h1,dh;
	int i,j;
	
	/* Calculate the Hamiltonian at the initial position */
	pe=ke=dh=0;
	kinetic_energy(num_dim,p,&ke);
	neg_logpost(&num_dim,x,&pe,g);
	h0=ke+pe;
	
	for(i=0;i<num_steps;++i)
	{
		/* take a half leap in momnetum */
		/* and a full step in position space */
		for(j=0;j<num_dim;++j)
		{
			p[j]=p[j]-0.5*epsilon*sigma[j]*g[j];
			x[j]=x[j]+epsilon*sigma[j]*p[j];
		}
		
		neg_logpost(&num_dim,x,&pe,g);
		
		/* take a full leap in momentum */
		for(j=0;j<num_dim;++j)
		{
			p[j]=p[j]-0.5*epsilon*sigma[j]*g[j];
		}
		
		/* calculate the new Hamiltonian */
		kinetic_energy(num_dim,p,&ke);
		h1=ke+pe;
		dh=h1-h0;
		
		/* exit the loop if dh above rejection threshold */
		if(-dh <= *log_ratio)
		{
			break;
		}
		
	}
	/* if the break condition did not occur */
	/* log_ratio should be the dh at the end of the loop */
	*log_ratio=dh;
	*llik=-pe;
}

void kinetic_energy(int num_dim,double* p,double* ke)
{
	int i;
	*ke=0;
	for(i=0;i<num_dim;++i)
	{
		(*ke)+=p[i]*p[i];
	}
	(*ke)*=0.5;
}
 

