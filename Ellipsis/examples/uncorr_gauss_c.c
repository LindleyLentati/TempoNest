#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <guided_hmc.h>

FILE* test_gauss_outfile;

void write_gauss_ghs_extract(int* ndim,double* x,double* val,double* g);
void write_gauss_ghs_extract_with_logpostval(int* ndim,double* x,double* logpostval,double* g);
void nd_uncorr_gauss_neg_log_post(int* ndim,double* x,double* v,double* g);

int main()
{

	const int ndim=1000;
	double* st;
	double* stp_sz;
	void (*nlp)(int*,double*,double*,double*);
	void (*wrt_ext)(int*,double*,double*,double*);
	double scl_fct;
	char* fl_pfx="uncorr_gauss_c";
	int seed;
	int fb_int;
	int max_stp;
	int resume;
	char ext_file_name[128];
	int i;
	int nburn=5000;
	int nsamp=1000;
	int doMaxLike = 0;
	/*
	Note that if nsamp>0, the sampler will be stopped forcefully after nsamp samples
	have been taken. This does not mean that the algorithm will have converged.
	Check the Hanson values for convergence.
	*/
	
	printf("\n=================================================\n");
	printf("|    Example program in C: 1                    |\n");
	printf("|    Uncorrelated N-D Gaussian                  |"); 
	printf("\n=================================================\n");
	
	
	
	st=(double*)malloc(ndim*sizeof(double));
	stp_sz=(double*)malloc(ndim*sizeof(double));
	
	for(i=0;i<ndim;++i)
	{
		st[i]=0.;
		stp_sz[i]=1;
	}
	
	
	nlp=&nd_uncorr_gauss_neg_log_post;
	/*wrt_ext=&write_gauss_ghs_extract;*/
	wrt_ext=&write_gauss_ghs_extract_with_logpostval;
	
	scl_fct=1.;
	fb_int=1000;
	max_stp=10;
	resume=1;
	seed=1234;
	
	strcpy(ext_file_name,fl_pfx);
	strcat(ext_file_name,".extract.dat");
	
	if(resume==1)
	{
		test_gauss_outfile=fopen(ext_file_name,"a");
	}
	else
	{
		test_gauss_outfile=fopen(ext_file_name,"w");
	}

	run_guided_hmc(ndim,st,scl_fct,max_stp,stp_sz,fl_pfx,seed,resume,fb_int,nlp,wrt_ext,nburn,nsamp, doMaxLike);
		
	if(test_gauss_outfile) fclose(test_gauss_outfile);
	
	free(st);
	free(stp_sz);
	
	return 0;
}

void write_gauss_ghs_extract(int* ndim,double* x,double* val,double* g)
{
	int i;
	if(test_gauss_outfile != NULL)
	{
		for(i=0;i<*ndim-1;++i)
		{
			fprintf(test_gauss_outfile,"%e,",x[i]);
		}
		fprintf(test_gauss_outfile,"%e\n",x[*ndim-1]);
	}
	else
	{
		printf("\nERROR IN WRITING GHS EXTRACT!\n");
	}
}

void write_gauss_ghs_extract_with_logpostval(int* ndim,double* x,double* logpostval,double* g)
{
	int i;
	if(test_gauss_outfile != NULL)
	{
		for(i=0;i<*ndim;++i)
		{
			fprintf(test_gauss_outfile,"%e,",x[i]);
		}
		fprintf(test_gauss_outfile,"%e\n",*logpostval);
	}
	else
	{
		printf("\nERROR IN WRITING GHS EXTRACT!\n");
	}
}

void nd_uncorr_gauss_neg_log_post(int* ndim,double* x,double* v,double* g)
{
	int i;
	(*v)=0;
	for(i=0;i<*ndim;++i)
	{
		(*v)+=x[i]*x[i];
		g[i]=x[i];
	}
	(*v)*=0.5;
}


