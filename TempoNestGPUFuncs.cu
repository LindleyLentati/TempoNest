#include <stdio.h>
#include <stdlib.h>
#include "/usr/include/gsl/gsl_sf_gamma.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include <cula_lapack_device.h>
#include <cula_blas_device.h>


#define BLOCK_SIZE 16

double *GlobalGmat_d;
float *GlobalGmatFloat_d;

// Matrices are stored in row-major order:
// M(row, col) = *(M.elements + row * M.width + col)

double iter_factorialGPU(unsigned int n)
{
    double ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}


void checkStatus(culaStatus status)
{
    char buf[256];

    if(!status)
        return;

    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);

    culaShutdown();
    exit(EXIT_FAILURE);
}

void checkStatusCarryOn(culaStatus status, int &CarryOn)
{
    char buf[256];

    if(!status){
    	CarryOn=0;
        return;
        }

    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    //printf("%s\n", buf);
    CarryOn=1;
    
    culaShutdown();
    culaStatus status2;
	status2 = culaInitialize();
	
}


void checkCudaError(cudaError_t err)
{
    if(!err)
        return;

    printf("%s\n", cudaGetErrorString(err));

    culaShutdown();
    exit(EXIT_FAILURE);
}


// simple kernel function that adds two vectors
__global__ void vect_add(double *a, double *b, int N)
{
   int Bidx = blockIdx.x;
   //a[Bidx*N+Tidx] = a[Bidx*N+Tidx] + b[Bidx*N+Tidx]; 
   for(int i =0; i < N; i++){
   		a[Bidx*N+i] = a[Bidx*N+i] + b[Bidx*N+i]; 
   		}
}


// simple kernel function that calcs det of a matrix
__global__ void calc_det(double *a, double *det, int N)
{
	
	det[0]=0;
   for(int i =0; i < N; i++){
   		det[0]+=log(a[i*N+i]);
   		}
   		det[0]=det[0]*2;
   		
   		
}

// simple kernel function that calcs det of a matrix
__global__ void Floatcalc_det(float *a, double *det, int N)
{
	
	det[0]=0;
   for(int i =0; i < N; i++){
   		det[0]+=log(a[i*N+i]);
   		}
   		det[0]=det[0]*2;
   		
   		
}

// simple kernel function that calcs det of a matrix
__global__ void calc_detFloat(float *a, double *det, int N)
{
	
	det[0]=0;
   for(int i =0; i < N; i++){
   		det[0]+=log(a[i*N+i]);
   		}
   		det[0]=det[0]*2;
   		
   		
}

__global__ void Makecov(double *A_d, double *BATvec, double *NoiseVec, double *SpecParm, int Aheight, int Awidth) {

	// Each thread computes one element of C
	// by accumulating results into Cvalue
	

	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=SpecParm[1];
	double gwampsquared=SpecParm[0];
	double covconst=SpecParm[2];
	

	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	 __syncthreads();
	if(row >= Aheight || col >= Awidth) return;
	double timdiff= BATvec[row] - BATvec[col];	
	double tau=2.0*M_PI*fabs(timdiff);
	double covsum=0;

	for(int k=0; k <=5; k++){
	
		    double ret = 1;
			for(unsigned int i = 1; i <= 2*k; ++i){
				ret *= (double)i;
			}
    
			covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(ret*(2*k+1-modelalpha));

	}

	A_d[row * Awidth + col]=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);

	if(row==col){
		A_d[row * Awidth + col] += NoiseVec[row];
	}

}

__global__ void MatMulKernel(int Arow,int Acol,int Brow, int Bcol,double *A,double *B,double *C)
{

	int Crow=Arow;
	int Ccol=Bcol;
    double Ctemp = 0.0;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    
    __syncthreads();

	if(row < Arow && col < Bcol) {

   		Ctemp = A[row] * B[col * Brow + row];
						  //GGTest[col*N + row]


	   C[col*Crow+row] = Ctemp;
	}
	  // 
}



extern "C" void WhiteMarginGPUWrapper_(double *Noise, double *Res, double *likeInfo, int N, int G)
{

	double *dettemp;
	dettemp = (double*)malloc(sizeof(double));
	dettemp[0]=0;
	
	
	double *GRes;
	GRes = (double*)malloc(sizeof(double)*G);
	double *WorkingGRes;
	WorkingGRes = (double*)malloc(sizeof(double)*G);

	// declare GPU copies
	double *Res_d;
	double *Noise_d;
	double *dettemp_d;
	
	double *NG_d;
	double *GG_d;
	double *GRes_d;



    cudaError_t err;
    culaStatus status;

   // Allocate memory on GPU
  	 err = cudaMalloc( (void **)&dettemp_d, sizeof(double) );
	 checkCudaError(err);
 	 err = cudaMalloc( (void **)&Res_d, sizeof(double)*N );
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&Noise_d, sizeof(double)*N );
	 checkCudaError(err);
	 
 	 err = cudaMalloc( (void **)&NG_d, sizeof(double)*N*G);
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&GG_d, sizeof(double)*G*G);
	 checkCudaError(err);	 
   	 err = cudaMalloc( (void **)&GRes_d, sizeof(double)*G*G);
	 checkCudaError(err);


   // copy vectors from CPU to GPU
   	 err = cudaMemcpy( dettemp_d, dettemp, sizeof(double), cudaMemcpyHostToDevice );
 	 checkCudaError(err);
   	 err = cudaMemcpy(Noise_d, Noise, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
     err = cudaMemcpy( Res_d, Res, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

  	 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	 dim3 dimGrid;
	 dimGrid.x=(N + dimBlock.x - 1)/dimBlock.x;
	 dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;

	 MatMulKernel<<<dimGrid, dimBlock>>>(N,N,N, G,Noise_d,GlobalGmat_d,NG_d);
 	
 	double alpha=1.0;
 	double beta=0.0; 
 	status =  culaDeviceDgemm('T', 'N', G, G, N, alpha, GlobalGmat_d, N, NG_d, N, beta, GG_d, G);
	checkStatus(status);


 	 status = culaDeviceDpotrf('L', G, GG_d, G);
     checkStatus(status);

     calc_det<<< 1, 1 >>>( GG_d, dettemp_d, G);
     err = cudaMemcpy( dettemp, dettemp_d, sizeof(double), cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
     likeInfo[0]=dettemp[0];
	
	 status = culaDeviceDgemv('T', N, G, alpha, GlobalGmat_d, N, Res_d, 1, beta, GRes_d, 1);
     checkStatus(status);
 	 err = cudaMemcpy(GRes, GRes_d, sizeof(double)*G, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);

 	 status=culaDeviceDpotrs('L', G, 1, GG_d, G, GRes_d, G);
	 checkStatus(status);
	 err = cudaMemcpy(WorkingGRes, GRes_d, sizeof(double)*G, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
	 
	 double sum=0;
	 for(int i=0; i<G;i++){sum=sum+GRes[i]*WorkingGRes[i];}
	 likeInfo[1]=sum;
	 
	 cudaFree(dettemp_d);
 	 cudaFree(Res_d);
	 cudaFree(Noise_d);
	 cudaFree(NG_d);
	 cudaFree(GG_d);
	 cudaFree(GRes_d);

	 
  	 free(GRes); 
  	 free(dettemp);
  	 free(WorkingGRes);


   return;
}


extern "C" void vHRedMarginGPUWrapper_(double *Res, double *BatVec, double *NoiseVec, double *SpecInfo, double *likeInfo, double *FactorialList, int N, int G)
{


	cudaError_t err;
	culaStatus status;
	
	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=SpecInfo[1];
	double covconst=gsl_sf_gamma(1-modelalpha)*sin(0.5*M_PI*modelalpha);
	
	double gwamp=pow(10.0,SpecInfo[0]);
	double gwampsquared=gwamp*gwamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-modelalpha)))/(pow(flo,(modelalpha-1)));
	
	SpecInfo[0]=gwampsquared;
	SpecInfo[2]=covconst;
	
	double *Res_d;
	double *BatVec_d;
	double *NoiseVec_d;
	double *SpecInfo_d;
	double *CovMatrix_d;

  
  	 err = cudaMalloc( (void **)&Res_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&BatVec_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&NoiseVec_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&SpecInfo_d, sizeof(double)*3 );
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&CovMatrix_d, sizeof(double)*N*N );
	 checkCudaError(err);
	 

     err = cudaMemcpy( Res_d, Res, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy( BatVec_d, BatVec, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy( NoiseVec_d, NoiseVec, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy(SpecInfo_d, SpecInfo, sizeof(double)*3, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

 	 
  	 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	 dim3 dimGrid;//((G + dimBlock.x - 1) / dimBlock.x,(N + dimBlock.y - 1) / dimBlock.y);
	 dimGrid.x=(N + dimBlock.x - 1)/dimBlock.x;
	 dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;
	 
	 Makecov<<<dimGrid, dimBlock>>>(CovMatrix_d, BatVec_d, NoiseVec_d, SpecInfo_d, N,N);
	 

	double *dettemp;
	dettemp = (double*)malloc(sizeof(double));
	dettemp[0]=0;
	
	
	double *GRes;
	GRes =(double*)malloc(sizeof(double)*G);
	double *WorkingGRes;
	WorkingGRes =(double*)malloc(sizeof(double)*G);
	

	// declare GPU copies

	double *CG_d;
	double *GCG_d;
	double *GRes_d;
	double *dettemp_d;


   // Allocate memory on GPU
  	 err = cudaMalloc( (void **)&dettemp_d, sizeof(double) );
	 checkCudaError(err);

	 
	  	 err = cudaMalloc( (void **)&CG_d, sizeof(double)*N*G );
	 checkCudaError(err);
	  	 err = cudaMalloc( (void **)&GCG_d, sizeof(double)*G*G );
	 checkCudaError(err);
	  	 err = cudaMalloc( (void **)&GRes_d, sizeof(double)*G );
	 checkCudaError(err);

   // copy vectors from CPU to GPU
   	 err = cudaMemcpy( dettemp_d, dettemp, sizeof(double), cudaMemcpyHostToDevice );
 	 checkCudaError(err);
   	// err = cudaMemcpy(CovMatrix_d, CovMatrix, sizeof(double)*N*N, cudaMemcpyHostToDevice );
 	// checkCudaError(err);

 	 
	int carryOn=0;
  	double alpha=1.0;
 	double beta=0.0; 

 	status =  culaDeviceDsymm('L', 'U', N, G, alpha, CovMatrix_d, N, GlobalGmat_d, N, beta, CG_d, N);
	checkStatus(status);
	//printf("done first linalg\n");
  	status =  culaDeviceDgemm('T', 'N', G, G, N, alpha, GlobalGmat_d, N, CG_d, N, beta, GCG_d, G);
	checkStatus(status);

	 status = culaDeviceDgemv('T', N, G, alpha, GlobalGmat_d, N, Res_d, 1, beta, GRes_d, 1);
     	checkStatus(status);
 	 err = cudaMemcpy(GRes, GRes_d, sizeof(double)*G, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);	 
 	 //	printf("do chol\n");
		
 	 status = culaDeviceDpotrf('L', G, GCG_d, G);
     	checkStatusCarryOn(status,carryOn);

	calc_det<<< 1, 1 >>>( GCG_d, dettemp_d, G);
	err = cudaMemcpy( dettemp, dettemp_d, sizeof(double), cudaMemcpyDeviceToHost);
	checkCudaError(err);
     	likeInfo[0]=dettemp[0];
    // printf("det: %g \n",likeInfo[0]);
    
	if(carryOn == 1){
     	// printf("Bad chol\n");
     	 
		cudaFree(dettemp_d);
		cudaFree(Res_d);
		cudaFree(BatVec_d);
		cudaFree(NoiseVec_d);
		cudaFree(SpecInfo_d);
		cudaFree(CovMatrix_d);
		cudaFree(CG_d);
		cudaFree(GCG_d);
		cudaFree(GRes_d);
	
		free(dettemp);
		free(GRes);
		free(WorkingGRes);
		
		
		return;

  	 }

 	 status=culaDeviceDpotrs('L', G, 1, GCG_d, G, GRes_d, G);
	 checkStatus(status);
	 err = cudaMemcpy(WorkingGRes, GRes_d, sizeof(double)*G, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);

	 double sum=0;
	 for(int i=0; i<G;i++){sum=sum+GRes[i]*WorkingGRes[i];}
	 likeInfo[1]=sum;
	 //printf("like: %g \n",sum);
	 
	 cudaFree(dettemp_d);
 	 cudaFree(Res_d);
  	 cudaFree(BatVec_d);
 	 cudaFree(NoiseVec_d);
 	 cudaFree(SpecInfo_d);
	 cudaFree(CovMatrix_d);
 	 cudaFree(CG_d);
	 cudaFree(GCG_d);
	 cudaFree(GRes_d);

  	 free(dettemp);
   	 free(GRes);
  	 free(WorkingGRes);
  


   return;
}

extern "C" void vHRedGPUWrapper_(double *SpecInfo, double *BatVec, double *Res, double *NoiseVec, double *likeInfo, int N)
{

	cudaError_t err;
	culaStatus status;

	double secday=24*60*60;
	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double modelalpha=SpecInfo[1];
	double covconst=gsl_sf_gamma(1-modelalpha)*sin(0.5*M_PI*modelalpha);
	
	double gwamp=pow(10.0,SpecInfo[0]);
	double gwampsquared=gwamp*gwamp*(pow((365.25*secday),2)/(12*M_PI*M_PI))*(pow(365.25,(1-modelalpha)))/(pow(flo,(modelalpha-1)));
	
	SpecInfo[0]=gwampsquared;
	SpecInfo[2]=covconst;
	
	double *Res_d;
	double *BatVec_d;
	double *NoiseVec_d;
	double *SpecInfo_d;
	double *CovMatrix_d;

  
  	 err = cudaMalloc( (void **)&Res_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&BatVec_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&NoiseVec_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&SpecInfo_d, sizeof(double)*3 );
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&CovMatrix_d, sizeof(double)*N*N );
	 checkCudaError(err);
	 

     	 err = cudaMemcpy( Res_d, Res, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy( BatVec_d, BatVec, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy( NoiseVec_d, NoiseVec, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy(SpecInfo_d, SpecInfo, sizeof(double)*3, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

 	 
  	 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	 dim3 dimGrid;//((G + dimBlock.x - 1) / dimBlock.x,(N + dimBlock.y - 1) / dimBlock.y);
	 dimGrid.x=(N + dimBlock.x - 1)/dimBlock.x;
	 dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;
	 
	 Makecov<<<dimGrid, dimBlock>>>(CovMatrix_d, BatVec_d, NoiseVec_d, SpecInfo_d, N,N);
	 

	double *dettemp;
	dettemp = (double*)malloc(sizeof(double));
	dettemp[0]=0;
	
	

	double *WorkingRes;
	WorkingRes = (double*)malloc(sizeof(double)*N);

	// declare GPU copies
	double *dettemp_d;



   // Allocate memory on GPU
  	 err = cudaMalloc( (void **)&dettemp_d, sizeof(double) );
	 checkCudaError(err);

   // copy vectors from CPU to GPU
   	 err = cudaMemcpy( dettemp_d, dettemp, sizeof(double), cudaMemcpyHostToDevice );
 	 checkCudaError(err);

	 int carryOn=0;
 	 status = culaDeviceDpotrf('L', N, CovMatrix_d, N);
     	 checkStatusCarryOn(status,carryOn);
     


	 calc_det<<< 1, 1 >>>( CovMatrix_d, dettemp_d, N);
	 err = cudaMemcpy( dettemp, dettemp_d, sizeof(double), cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
     	 likeInfo[0]=dettemp[0];
     
          if(carryOn == 1){

     	 	 cudaFree(dettemp_d);
	 	 cudaFree(Res_d);
		 cudaFree(BatVec_d);
		 cudaFree(NoiseVec_d);
		 cudaFree(SpecInfo_d);
		 cudaFree(CovMatrix_d);
	  	 free(dettemp);
	  	 free(WorkingRes);
	  	 return;
  	 }
	

 	 status=culaDeviceDpotrs('L', N, 1, CovMatrix_d, N, Res_d, N);
	 checkStatus(status);
	 err = cudaMemcpy(WorkingRes, Res_d, sizeof(double)*N, cudaMemcpyDeviceToHost);

  	 checkCudaError(err);

	 double sum=0;
	 	
	for(int i=0; i<N;i++){
		sum=sum+Res[i]*WorkingRes[i];
	}

	 likeInfo[1]=sum;
	 
	 cudaFree(dettemp_d);
 	 cudaFree(Res_d);
  	 cudaFree(BatVec_d);
 	 cudaFree(NoiseVec_d);
 	 cudaFree(SpecInfo_d);
	 cudaFree(CovMatrix_d);

  	 free(dettemp);
  	 free(WorkingRes);
  	 
  


   return;
}


// simple kernel function that calculates the FMatrix
__global__ void make_fmatrix(double *FMatrix_d,double *Freqs_d, double *BATvec_d, int N,int F)
{

	int Bidx = blockIdx.x;
	
	for(int i=0;i<F/2;i++){
		
			FMatrix_d[i*N + Bidx]=cos(2*M_PI*Freqs_d[i]*BATvec_d[Bidx]);
			FMatrix_d[(i+F/2)*N + Bidx]=sin(2*M_PI*Freqs_d[i]*BATvec_d[Bidx]);
	}


}
__global__ void fastmake_fmatrix(double *FMatrix_d,double *Freqs_d, double *BATvec_d, int Aheight,int Awidth) {

	// Each thread computes one element of F
	// by accumulating results into Cvalue


	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	 __syncthreads();
	if(row >= Aheight || col >= Awidth) return;

	FMatrix_d[row * Awidth + col]=cos(2*M_PI*Freqs_d[col]*BATvec_d[row]);
	FMatrix_d[row*Awidth + col + Awidth/2]=sin(2*M_PI*Freqs_d[col]*BATvec_d[row]);

}



extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, int N, int F){

	double *Freqs_d;
	double *resvec_d;
	double *BATvec_d;
	double *Noise_d;
	
	double *FMatrix_d;
	double *NF_d;	
	double *FNF_d;
	double *NFd_d;
	
	double *FNFvec;
	FNFvec = (double*)malloc(sizeof(double)*F*F);

	cudaError_t err;
	culaStatus status;
      
  	 err = cudaMalloc( (void **)&Freqs_d, sizeof(double)*F );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&resvec_d, sizeof(double)*N );
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&BATvec_d, sizeof(double)*N );
	 checkCudaError(err);
   	 err = cudaMalloc( (void **)&Noise_d, sizeof(double)*N );
	 checkCudaError(err);
	 
   	 err = cudaMalloc( (void **)&FMatrix_d, sizeof(double)*N*F );
	 checkCudaError(err);
   	 err = cudaMalloc( (void **)&NF_d, sizeof(double)*N*F );
	 checkCudaError(err);
   	 err = cudaMalloc( (void **)&FNF_d, sizeof(double)*F*F );
	 checkCudaError(err);	 
   	 err = cudaMalloc( (void **)&NFd_d, sizeof(double)*F);
	 checkCudaError(err);	
	 
         // copy vectors from CPU to GPU
   	 err = cudaMemcpy( Freqs_d, Freqs, sizeof(double)*F, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
	 err = cudaMemcpy(resvec_d, resvec, sizeof(double)*N, cudaMemcpyHostToDevice );
	 checkCudaError(err);
   	 err = cudaMemcpy(BATvec_d, BATvec, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy( Noise_d, Noise, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

 	 
// 	 make_fmatrix<<< N, 1 >>>(FMatrix_d,Freqs_d,BATvec_d,N,F);
  	 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	 dim3 dimGrid;

	 dimGrid.x=(F + dimBlock.x - 1)/dimBlock.x;
	 dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;
	 fastmake_fmatrix<<<dimGrid, dimBlock>>>(FMatrix_d,Freqs_d,BATvec_d,N,F);

	 MatMulKernel<<<dimGrid, dimBlock>>>(N,N,N, F,Noise_d,FMatrix_d,NF_d);

 	 
	double alpha=1.0;
	double beta=0.0; 
	status =  culaDeviceDgemm('T', 'N', F, F, N, alpha, FMatrix_d, N, NF_d, N, beta, FNF_d, F);
	checkStatus(status);
	
 	 status = culaDeviceDgemv('T', N, F, alpha, NF_d, N, resvec_d, 1, beta, NFd_d, 1);
	 checkStatus(status);
	
	 err = cudaMemcpy(FNFvec, FNF_d, sizeof(double)*F*F, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
  	 	
	 err = cudaMemcpy(NFd, NFd_d, sizeof(double)*F, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);

	for(int f1=0;f1<F; f1++){
		for(int f2=0;f2<F; f2++){

			FNF[f2][f1]=FNFvec[f1*F + f2];
		}
	}

	cudaFree(Freqs_d);
	cudaFree(BATvec_d);
	cudaFree(Noise_d);
    	cudaFree(FMatrix_d);
    	cudaFree(NF_d);
	cudaFree(FNF_d);
	cudaFree(resvec_d);
	cudaFree(NFd_d);
	
	free(FNFvec);
		
}
 	 


extern "C" void LRedMarginGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *Noise, double **FNF, double *NFd, double *likeVals, int N, int F, int G){

	double *dettemp;
	dettemp = (double*)malloc(sizeof(double));
	dettemp[0]=0;
	
	double *dettemp_d;

	double *Freqs_d;
	double *resvec_d;
	double *BATvec_d;
	double *Noise_d;
	
		
	double *NG_d;
	double *GG_d;
	double *GNG_d;
	double *GNGd_d;
	
	double *GNGd;
	GNGd = (double*)malloc(sizeof(double)*N);

	
	double *FMatrix_d;
	double *NF_d;	
	double *FNF_d;
	double *NFd_d;
	
	double *FNFvec;
	FNFvec = (double*)malloc(sizeof(double)*F*F);

	cudaError_t err;
	culaStatus status;
    
  	 err = cudaMalloc( (void **)&dettemp_d, sizeof(double) );
	 checkCudaError(err);
      
  	 err = cudaMalloc( (void **)&Freqs_d, sizeof(double)*F );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&resvec_d, sizeof(double)*N );
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&BATvec_d, sizeof(double)*N );
	 checkCudaError(err);
   	 err = cudaMalloc( (void **)&Noise_d, sizeof(double)*N );
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&GNGd_d, sizeof(double)*N );
	 checkCudaError(err);

	 
   	 err = cudaMalloc( (void **)&NG_d, sizeof(double)*N*G);
	 checkCudaError(err);
	 err = cudaMalloc( (void **)&GG_d, sizeof(double)*G*G );
	 checkCudaError(err);
  	 err = cudaMalloc( (void **)&GNG_d, sizeof(double)*N*N );
	 checkCudaError(err);
	 
   	 err = cudaMalloc( (void **)&FMatrix_d, sizeof(double)*N*F );
	 checkCudaError(err);
   	 err = cudaMalloc( (void **)&NF_d, sizeof(double)*N*F );
	 checkCudaError(err);
   	 err = cudaMalloc( (void **)&FNF_d, sizeof(double)*F*F );
	 checkCudaError(err);	 
   	 err = cudaMalloc( (void **)&NFd_d, sizeof(double)*F);
	 checkCudaError(err);	
	 
   // copy vectors from CPU to GPU
   	 err = cudaMemcpy( Freqs_d, Freqs, sizeof(double)*F, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
	 err = cudaMemcpy(resvec_d, resvec, sizeof(double)*N, cudaMemcpyHostToDevice );
	 checkCudaError(err);
   	 err = cudaMemcpy(BATvec_d, BATvec, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 err = cudaMemcpy( Noise_d, Noise, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);
 	 
  	 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	 dim3 dimGrid;
	 dimGrid.x=(N + dimBlock.x - 1)/dimBlock.x;
	 dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;

	 MatMulKernel<<<dimGrid, dimBlock>>>(N,N,N, G,Noise_d,GlobalGmat_d,NG_d);
	 

 	 double alpha=1.0;
 	 double beta=0.0; 
 	 status =  culaDeviceDgemm('T', 'N', G, G, N, alpha, GlobalGmat_d, N, NG_d, N, beta, GG_d, G);
	 checkStatus(status);


 	 status = culaDeviceDpotrf('L', G, GG_d, G);
     	checkStatus(status);

	calc_det<<< 1, 1 >>>( GG_d, dettemp_d, G);
	err = cudaMemcpy( dettemp, dettemp_d, sizeof(double), cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
     	likeVals[0]=dettemp[0];
     
  	 status = culaDeviceDpotri('L', G, GG_d, G);
     	checkStatus(status);
     
  	 status =  culaDeviceDsymm('R', 'L', N, G, alpha, GG_d, G, GlobalGmat_d, N, beta, NG_d, N);
	 checkStatus(status);
	 
  	 status =  culaDeviceDgemm('N', 'T', N, N, G, alpha, NG_d, N, GlobalGmat_d, N, beta, GNG_d, N);
	 checkStatus(status);
	 
  	 status = culaDeviceDgemv('N', N, N, alpha, GNG_d, N, resvec_d, 1, beta, GNGd_d, 1);
     	checkStatus(status);
     

	 err = cudaMemcpy(GNGd, GNGd_d, sizeof(double)*N, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
  	 likeVals[1]=0;
  	 for(int i =0;i < N; i++){likeVals[1] += resvec[i]*GNGd[i];}

 	 make_fmatrix<<< N, 1 >>>(FMatrix_d,Freqs_d,BATvec_d,N,F);
 	 
   	 status =  culaDeviceDgemm('N', 'N', N, F, N, alpha, GNG_d, N, FMatrix_d, N, beta, NF_d, N);
	 checkStatus(status);
	 
  	 status =  culaDeviceDgemm('T', 'N', F, F, N, alpha, FMatrix_d, N, NF_d, N, beta, FNF_d, F);
	 checkStatus(status);
	
 	 status = culaDeviceDgemv('T', N, F, alpha, NF_d, N, resvec_d, 1, beta, NFd_d, 1);
     checkStatus(status);
	
	 err = cudaMemcpy(FNFvec, FNF_d, sizeof(double)*F*F, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);
  	 	
	 err = cudaMemcpy(NFd, NFd_d, sizeof(double)*F, cudaMemcpyDeviceToHost);
  	 checkCudaError(err);

		for(int f1=0;f1<F; f1++){
			for(int f2=0;f2<F; f2++){
 
				FNF[f2][f1]=FNFvec[f1*F + f2];
				//printf("GPUFNF: %i %i %g \n",f1,f2);
			}
		}

	cudaFree(dettemp_d);
	cudaFree(Freqs_d);
	cudaFree(BATvec_d);
	cudaFree(Noise_d);
    cudaFree(FMatrix_d);
    cudaFree(NF_d);
	cudaFree(FNF_d);
	cudaFree(resvec_d);
	cudaFree(NFd_d);
	
    cudaFree(NG_d);
    cudaFree(GG_d);
	cudaFree(GNG_d);
	cudaFree(GNGd_d);

	
	free(FNFvec);
	free(dettemp);
	free(GNGd);
	
}
 	 

 extern "C" void copy_floatgmat_(float *G, int N){

    cudaError_t err;

   // Allocate memory on GPU
	printf("copying G\n");


 	 err = cudaMalloc( (void **)&GlobalGmatFloat_d, sizeof(float)*N );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalGmatFloat_d, G, sizeof(float)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}


extern "C" void copy_gmat_(double *G, int N){

    cudaError_t err;

   // Allocate memory on GPU
	printf("copying G\n");


 	 err = cudaMalloc( (void **)&GlobalGmat_d, sizeof(double)*N );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalGmat_d, G, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}

