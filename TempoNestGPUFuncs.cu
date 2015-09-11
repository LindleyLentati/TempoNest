#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "dgesvd.h"
#include <cula_lapack_device.h>
#include <cula_blas_device.h>


#define BLOCK_SIZE 16

//double *GlobalGmat_d;
//double *GlobalStaticGmat_d;
//double *GlobalStaticUGmat_d;
//float *GlobalGmatFloat_d;
//double *GlobalStaticDmat_d;
double *GlobalTotalMatrix_d;

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

// simple kernel function that calcs det of a matrix
__global__ void calc_DiagLike(double *Vec, double *Noise, int N, double *val)
{
	
	val[0]=0;
	for(int i =0; i < N; i++){
   		val[0]+=Vec[i]*Vec[i]*Noise[i];
	}

   		
   		
}

// simple kernel function that calcs det of a matrix
__global__ void calc_DotLike(double *Vec1, double *Vec2, int N, double *val)
{
	
	val[0]=0;
	for(int i =0; i < N; i++){
   		val[0]+=Vec1[i]*Vec2[i];
   		//printf("GPU copy %i %g %g\n", i,Vec1[i],Vec2[i]);
	}

   		
   		
}

// simple kernel function that calcs det of a matrix
__global__ void copyvec(double *Vec1, double *Vec2, int N)
{
	
		int Bidx = blockIdx.x;
		 __syncthreads();
   		Vec1[Bidx]=Vec2[Bidx];
   		//printf("copy: %i %g\n",Bidx, Vec1[Bidx]);
   		
}

/*
__global__ void Makecov(double *A_d, double *BATvec, double *NoiseVec, double *SpecParm, int Aheight, int Awidth) {

	// Each thread computes one element of C
	// by accumulating results into Cvalue
	

	double LongestPeriod=1.0/pow(10.0,-5); //
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

	for(int k=0; k <=10; k++){
	
		    double ret = 1;
			for(unsigned int i = 1; i <= 2*k; ++i){
				ret *= (double)i;
			//	printf("Ret: %i %g \n",i,ret);
			}
    
			covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(ret*(2*k+1-modelalpha));
			//printf("covsum: %i %i %i %g \n",row,col,k,covsum);

	}

	A_d[row * Awidth + col]=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);

	if(row==col){
		A_d[row * Awidth + col] += NoiseVec[row];
	}

}


__global__ void MakeDMcov(double *A_d, double *BATvec, double *NoiseVec, double *DMVec, double *SpecParm, int Aheight, int Awidth) {

	// Each thread computes one element of C
	// by accumulating results into Cvalue
	

	double LongestPeriod=1.0/pow(10.0,-5);
	double flo=1.0/LongestPeriod;

	double gwampsquared=SpecParm[0];
	double modelalpha=SpecParm[1];
	double covconst=SpecParm[2];
	
	double dmampsquared=SpecParm[3];
	double dmmodelalpha=SpecParm[4];
	double dmcovconst=SpecParm[5];
	

	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	 __syncthreads();
	if(row >= Aheight || col >= Awidth) return;
	double timdiff= BATvec[row] - BATvec[col];	
	double tau=2.0*M_PI*fabs(timdiff);
	double covsum=0;
	double dmcovsum=0;

	for(int k=0; k <=5; k++){
	
		    double ret = 1;
			for(unsigned int i = 1; i <= 2*k; ++i){
				ret *= (double)i;
			}
    
			covsum=covsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(ret*(2*k+1-modelalpha));
			
			dmcovsum=dmcovsum+pow(-1.0,k)*(pow(flo*tau,2*k))/(ret*(2*k+1-dmmodelalpha));

	}
	
	double gwpart=0;
	if(SpecParm[0] !=0 )gwpart=gwampsquared*(covconst*pow((flo*tau),(modelalpha-1)) - covsum);
	
	double dmpart=0;
	if(SpecParm[3] !=0 )dmpart=dmampsquared*(dmcovconst*pow((flo*tau),(dmmodelalpha-1)) - dmcovsum)*DMVec[row]*DMVec[col];

	A_d[row * Awidth + col]= gwpart+dmpart;

	if(row==col){
		A_d[row * Awidth + col] += NoiseVec[row];
	}

	//printf("%i %i %g\n",row,col,A_d[row * Awidth + col]);

}

*/
__global__ void MatMulKernel(int Arow,int Acol,int Brow, int Bcol,double *A,double *B,double *C)
{

	int Crow=Arow;
	int Ccol=Bcol;
    double Ctemp = 0.0;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    
    __syncthreads();

	if(row < Arow && col < Bcol) {
		//if(row<32)printf("GPUNT: %i %i %g %g \n", row, col, B[col * Brow + row], A[row] );
   		Ctemp = A[row] * B[col * Brow + row];
						  //GGTest[col*N + row]


	   C[col*Crow+row] = Ctemp;
	}
	  // 
}


__global__ void SimpleDiagMatMulKernel(int N,int T,double *Noise_d,double *TMatrix_d,double *NTMatrix_d)
{

	
    
    __syncthreads();


		for(int i=0;i<T; i++){
			for(int j=0;j<N; j++){
				//if(i ==28)printf("GPU SDMMK %i %i %g %g \n",i,j,TMatrix_d[i*N + j],Noise_d[j]);
				NTMatrix_d[i*N + j]=TMatrix_d[i*N + j]*Noise_d[j];
			}
		}		
}


/*
extern "C" void WhiteMarginGPUWrapper_(void *context, double *TNDMVec, double *resvec, double *Noise, int N, int D, int NTime, int NJumps, double *likevals){



	double *resvec_d;
	double *Noise_d;


	double *DMatrix_d;
	double *NT_d;	
	double *TNT_d;
	double *NTd_d;


	cudaError_t err;
	culaStatus status;

	// allocate memory on GPU
	err = cudaMalloc( (void **)&resvec_d, sizeof(double)*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&Noise_d, sizeof(double)*N );
	checkCudaError(err);



	err = cudaMalloc( (void **)&DMatrix_d, sizeof(double)*N*D );
	checkCudaError(err);
	err = cudaMalloc( (void **)&NT_d, sizeof(double)*N*D );
	checkCudaError(err);
	err = cudaMalloc( (void **)&TNT_d, sizeof(double)*D*D );
	checkCudaError(err);	 
	err = cudaMalloc( (void **)&NTd_d, sizeof(double)*D);
	checkCudaError(err);	

	// copy vectors from CPU to GPU
	err = cudaMemcpy(resvec_d, resvec, sizeof(double)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
	err = cudaMemcpy( Noise_d, Noise, sizeof(double)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);

 	 

	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid;
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	


	if(D != NTime+NJumps){
	
		err = cudaMemcpy( DMatrix_d, TNDMVec, sizeof(double)*D*N, cudaMemcpyHostToDevice );
		checkCudaError(err);
	
		double *U_d;
		double *V_d;
		double *S_d;
	
		err = cudaMalloc( (void **)&U_d, sizeof(double)*N*N );
		checkCudaError(err);
		err = cudaMalloc( (void **)&V_d, sizeof(double)*D*D );
		checkCudaError(err);
		err = cudaMalloc( (void **)&S_d, sizeof(double)*D );
		checkCudaError(err);
	
	
		culaDeviceDgesvd('O','N', N, D, DMatrix_d, N, S_d, U_d, N, V_d, D);
	
		cudaFree(V_d);
		cudaFree(S_d);
		cudaFree(U_d);
		
		cudaDeviceSynchronize();	
	}
	else{
		DMatrix_d=GlobalStaticDmat_d;
		cudaDeviceSynchronize();	
	}
	


///////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Do the Algebra///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	
 	 
 	 
	//printf("entered 5 %i %i\n",T,N);
	
	
	dimGrid.x=(D + dimBlock.x - 1)/dimBlock.x;
	dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;	

	MatMulKernel<<<dimGrid, dimBlock>>>(N,N,N, D,Noise_d,DMatrix_d,NT_d);
	//SimpleDiagMatMulKernel<<<1,1>>>(N, T, Noise_d, TMatrix_d, NT_d);
	cudaDeviceSynchronize();

	double alpha=1.0;
	double beta=0.0; 
	

	status =  culaDeviceDgemm('T', 'N', D, D, N, alpha, DMatrix_d, N, NT_d, N, beta, TNT_d, D);
	checkStatus(status);

	cudaDeviceSynchronize();
	

	status = culaDeviceDgemv('T', N, D, alpha, NT_d, N, resvec_d, 1, beta, NTd_d, 1);
	checkStatus(status);
	
	cudaDeviceSynchronize();
	

	
	//printf("entered 6: %i \n", T);
	double *dettemp_d;
	double *tempval=new double[1];
	err = cudaMalloc( (void **)&dettemp_d, sizeof(double) );
	checkCudaError(err);

	int carryOn=0;
	status = culaDeviceDpotrf('L', D, TNT_d, D);
	cudaDeviceSynchronize();
	//printf("entered 6.25\n");
	
	checkStatusCarryOn(status,carryOn);
	
	cudaDeviceSynchronize();

	//printf("entered 6.5\n");
	calc_det<<< 1, 1 >>>( TNT_d, dettemp_d, D);
	err = cudaMemcpy( tempval, dettemp_d, sizeof(double), cudaMemcpyDeviceToHost);
	checkCudaError(err);
	likevals[0]=tempval[0];
	
	cudaDeviceSynchronize();
     
	if(carryOn == 1){


		cudaFree(resvec_d);
		cudaFree(Noise_d);



		cudaFree(DMatrix_d);
		cudaFree(NT_d);
		cudaFree(TNT_d);
		cudaFree(NTd_d);
		
		cudaFree(dettemp_d);
		delete(tempval);
		return;
	}

	//printf("entered 7\n");
	double *WorkVec_d;
	err = cudaMalloc( (void **)&WorkVec_d, sizeof(double)*D );
	checkCudaError(err);
	copyvec<<< D, 1 >>>(WorkVec_d, NTd_d, D);
	
	cudaDeviceSynchronize();
	

	status=culaDeviceDpotrs('L', D, 1, TNT_d, D, WorkVec_d, D);
	checkStatus(status);
	
	cudaDeviceSynchronize();
	
	double *freqlike_d;
	err = cudaMalloc( (void **)&freqlike_d, sizeof(double));
	checkCudaError(err);

	calc_DotLike<<< 1, 1 >>>(WorkVec_d, NTd_d, D, freqlike_d);
	
	cudaDeviceSynchronize();
	err = cudaMemcpy( tempval, freqlike_d, sizeof(double), cudaMemcpyDeviceToHost);
	checkCudaError(err);
	likevals[1]=tempval[0];

	


	cudaFree(resvec_d);
	cudaFree(Noise_d);

	cudaFree(DMatrix_d);
	cudaFree(NT_d);
	cudaFree(TNT_d);
	cudaFree(NTd_d);
	
	cudaFree(WorkVec_d);
	cudaFree(dettemp_d);
	delete(tempval);

}
*/


// simple kernel function that calculates the FMatrix
__global__ void make_fmatrix(double *TMatrix_d,double *Freqs_d, double *ObsFreqs_d, double *BATvec_d, double *DMVec_d, int *SysGroups_d, int *BandFreqs, int N,int RF, int DMF, int BandNoiseCoeff, int GroupNoiseCoeff,  int incRED, int incDM, int incBandNoise, int incGroupNoise, int ReplaceTMatrix, int TimetoMargin, int numTime, int numJumps, double *DMatrix_d)
{

	int Bidx = blockIdx.x;


	if(TimetoMargin != numJumps + numTime){
		for(int i=0;i<TimetoMargin;i++){
			TMatrix_d[i*N + Bidx]=DMatrix_d[i*N + Bidx];
		}
	}

	int startpos=0;
	if(incRED !=0){
		if(ReplaceTMatrix==0){
			for(int i=0;i<RF/2;i++){
				TMatrix_d[(TimetoMargin+i)*N + Bidx]=cos(2*M_PI*Freqs_d[i]*BATvec_d[Bidx]);
				TMatrix_d[(TimetoMargin+i+RF/2)*N + Bidx]=sin(2*M_PI*Freqs_d[i]*BATvec_d[Bidx]);
			}
		}
		startpos=RF;
	}

	
      if(incDM !=0){
		if(ReplaceTMatrix==0){
		        for(int i=0;i<DMF/2;i++){
				//if(Bidx==0)printf("D: %i %i %g %g \n", Bidx,i,1.0/Freqs_d[i], DMVec_d[Bidx]);
		                TMatrix_d[(TimetoMargin+startpos+i)*N + Bidx]=cos(2*M_PI*Freqs_d[startpos+i]*BATvec_d[Bidx])*DMVec_d[Bidx];
		                TMatrix_d[(TimetoMargin+startpos+i+DMF/2)*N + Bidx]=sin(2*M_PI*Freqs_d[startpos+i]*BATvec_d[Bidx])*DMVec_d[Bidx];
		        }
		}
		startpos=startpos+DMF;
       }

	if(incBandNoise > 0){

		
		for(int b = 0; b < incBandNoise; b++){	

			if(ReplaceTMatrix==0){
				int startfreq = BandFreqs[b*3+0];
				int stopfreq = BandFreqs[b*3+1];
				int BandScale = BandFreqs[b*3+2];


				for(int i=0;i<BandNoiseCoeff/2;i++){
					if(ObsFreqs_d[Bidx] > startfreq && ObsFreqs_d[Bidx] < stopfreq){
						TMatrix_d[(TimetoMargin+startpos+i)*N + Bidx]=cos(2*M_PI*Freqs_d[startpos+i]*BATvec_d[Bidx]);
						TMatrix_d[(TimetoMargin+startpos+i+BandNoiseCoeff/2)*N + Bidx]=sin(2*M_PI*Freqs_d[startpos+i]*BATvec_d[Bidx]);
					}
					else{
			
						TMatrix_d[(TimetoMargin+startpos+i)*N + Bidx]=0;
						TMatrix_d[(TimetoMargin+startpos+i+BandNoiseCoeff/2)*N + Bidx]=0;
		
					}
				}
			}
			startpos=startpos+BandNoiseCoeff;

		}


       }

	//printf("In GPU : %i \n", incGroupNoise); 
	if(incGroupNoise > 0){
		for(int g = 0; g < incGroupNoise; g++){

			for(int i=0;i<GroupNoiseCoeff/2;i++){
				//printf("GPU Groups %i %i \n", Bidx, SysGroups_d[Bidx]);
				if(SysGroups_d[Bidx] == g+1){
					TMatrix_d[(TimetoMargin+startpos+i)*N + Bidx]=cos(2*M_PI*Freqs_d[startpos+i]*BATvec_d[Bidx]);
					TMatrix_d[(TimetoMargin+startpos+i+GroupNoiseCoeff/2)*N + Bidx]=sin(2*M_PI*Freqs_d[startpos+i]*BATvec_d[Bidx]);
				}
				else{
		
					TMatrix_d[(TimetoMargin+startpos+i)*N + Bidx]=0;
					TMatrix_d[(TimetoMargin+startpos+i+GroupNoiseCoeff/2)*N + Bidx]=0;
	
				}
			}



			startpos=startpos+GroupNoiseCoeff;

		}

	}
}
/*


// simple kernel function that calculates the FMatrix
__global__ void add_EcorrToFMatrix(double *FMatrix_d, double *EMatrix_d, int Nobs, int FSize, int EpochSize)
{

	int Bidx = blockIdx.x;
        for(int i=0;i<EpochSize;i++){
                FMatrix_d[(FSize+i)*Nobs + Bidx]= EMatrix_d[i*Nobs + Bidx];
        }


}
*/
/*

// simple kernel function that calculates the FMatrix
__global__ void make_DMfmatrix(double *FMatrix_d,double *Freqs_d, double *BATvec_d, double *DMVec_d, int N,int F)
{

	int Bidx = blockIdx.x;
	
	for(int i=0;i<F/4;i++){
		
			FMatrix_d[i*N + Bidx]=cos(2*M_PI*Freqs_d[i]*BATvec_d[Bidx]);
			FMatrix_d[(i+F/4)*N  + Bidx]=cos(2*M_PI*Freqs_d[i]*BATvec_d[Bidx])*DMVec_d[Bidx];
			FMatrix_d[(i+F/2)*N + Bidx]=sin(2*M_PI*Freqs_d[i]*BATvec_d[Bidx]);
			FMatrix_d[(i+3*F/4)*N + Bidx]=sin(2*M_PI*Freqs_d[i]*BATvec_d[Bidx])*DMVec_d[Bidx];
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
*/

/*

extern "C" void LRedGPUWrapper_(double *Freqs, double *resvec, double *BATvec, double *DMVec, double *Noise, double **FNF, double *NFd, int N, int RF,int DMF, int F, int incRED, int incDM){

	double *Freqs_d;
	double *resvec_d;
	double *BATvec_d;
	double *Noise_d;
	double *DMVec_d;
	
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
     err = cudaMalloc( (void **)&DMVec_d, sizeof(double)*N );
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
     err = cudaMemcpy( DMVec_d, DMVec, sizeof(double)*N, cudaMemcpyHostToDevice );
     checkCudaError(err);

 	 
// 	 make_fmatrix<<< N, 1 >>>(FMatrix_d,Freqs_d,BATvec_d,N,F);
  	 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	 dim3 dimGrid;

	 dimGrid.x=(F + dimBlock.x - 1)/dimBlock.x;
	 dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;
	 //fastmake_fmatrix<<<dimGrid, dimBlock>>>(FMatrix_d,Freqs_d,BATvec_d,N,F);
 	 //make_fmatrix<<< N, 1 >>>(FMatrix_d,Freqs_d,BATvec_d,DMVec_d,N,RF,DMF, 0, incRED, incDM, 0);

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
	cudaFree(DMVec_d);
	free(FNFvec);
		
}
*/
/*
// simple kernel function that calculates the TMatrix
__global__ void make_Tmatrix(double *TMatrix_d, double *DMatrix_d, double *FMatrix_d, int N, int T, int D, int F)
{


 //   int row = blockIdx.y * blockDim.y + threadIdx.y;
  //  int col = blockIdx.x * blockDim.x + threadIdx.x;
     __syncthreads();
     
 //	int Bidx = blockIdx.x;
    int Bidx = blockIdx.y * blockDim.y + threadIdx.y;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

//	for(int i=0;i<T;i++){
		if(i<D){
			TMatrix_d[i*N + Bidx]=DMatrix_d[i*N + Bidx];
		}
		else{
			if(i==D)printf("i==D %i
			TMatrix_d[i*N + Bidx]=FMatrix_d[(i-D)*N + Bidx];
		}
//	}
	
	
	
	 	int Bidx = blockIdx.x;


	for(int i=0;i<T;i++){
		if(i<D){
			TMatrix_d[i*N + Bidx]=DMatrix_d[i*N + Bidx];
		}
		else{
			//if(i==D)printf("i==D %i %g \n", Bidx, FMatrix_d[(i-D)*N + Bidx]);
			TMatrix_d[i*N + Bidx]=FMatrix_d[(i-D)*N + Bidx];
		}
	}

}

*/
// simple kernel function that adds powercoeff to TNT
__global__ void addCoeffsKernel(int T, int D,int F,double *TNT_d, double *powercoeffs_d)
{


     __syncthreads();
    for(int i =0; i < F; i++){
    	
    	TNT_d[(D+i)*T+D+i]+=1.0/powercoeffs_d[i];
    }


}


extern "C" void NewLRedMarginGPUWrapper_(void *context, double *TNDMVec, double *Freqs, double *ObsFreqs, double *powercoeff, double *resvec, double *BATvec, double *DMVec, double *Noise, int *SysGroups, int N, int RF,int DMF, int BandNoiseCoeff, int GroupNoiseCoeff, int D, int F, int T, int incRED, int incDM, int incBandNoise, int incGroupNoise, int NTime, int NJumps, double *likevals, int incNGJitter, int numNGJitterEpochs, int *BandInfo, int ReplaceTMatrix){


	//printf("entered 1\n");
	double *Freqs_d;
	double *ObsFreqs_d;
	double *powercoeff_d;
	double *resvec_d;
	double *BATvec_d;
	double *Noise_d;
	double *DMVec_d;
	int *SysGroups_d;
	int *BandInfo_d;
	
//	double *FMatrix_d;
	double *DMatrix_d;
//	double *TMatrix_d;
	double *NT_d;	
	double *TNT_d;
	double *NTd_d;


	cudaError_t err;
	culaStatus status;

	// allocate memory on GPU
	err = cudaMalloc( (void **)&Freqs_d, sizeof(double)*F );
	checkCudaError(err);
        err = cudaMalloc( (void **)&ObsFreqs_d, sizeof(double)*N );
        checkCudaError(err);
	err = cudaMalloc( (void **)&powercoeff_d, sizeof(double)*F );
	checkCudaError(err);
	err = cudaMalloc( (void **)&resvec_d, sizeof(double)*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&BATvec_d, sizeof(double)*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&Noise_d, sizeof(double)*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&DMVec_d, sizeof(double)*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&SysGroups_d, sizeof(int)*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&BandInfo_d, sizeof(int)*3*incBandNoise);
	checkCudaError(err);


//	err = cudaMalloc( (void **)&FMatrix_d, sizeof(double)*N*F );
//	checkCudaError(err);

//	err = cudaMalloc( (void **)&TMatrix_d, sizeof(double)*N*T );
//	checkCudaError(err);


	err = cudaMalloc( (void **)&NT_d, sizeof(double)*N*T );
	checkCudaError(err);
	err = cudaMalloc( (void **)&TNT_d, sizeof(double)*T*T );
	checkCudaError(err);	 
	err = cudaMalloc( (void **)&NTd_d, sizeof(double)*T);
	checkCudaError(err);	

	// copy vectors from CPU to GPU
	err = cudaMemcpy( Freqs_d, Freqs, sizeof(double)*F, cudaMemcpyHostToDevice );
	checkCudaError(err);
        err = cudaMemcpy( ObsFreqs_d, ObsFreqs, sizeof(double)*N, cudaMemcpyHostToDevice );
        checkCudaError(err);
	err = cudaMemcpy( powercoeff_d, powercoeff, sizeof(double)*F, cudaMemcpyHostToDevice );
	checkCudaError(err);
	err = cudaMemcpy(resvec_d, resvec, sizeof(double)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
	err = cudaMemcpy(BATvec_d, BATvec, sizeof(double)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
	err = cudaMemcpy( Noise_d, Noise, sizeof(double)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
	err = cudaMemcpy( DMVec_d, DMVec, sizeof(double)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
	err = cudaMemcpy( SysGroups_d, SysGroups, sizeof(int)*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
 	err = cudaMemcpy( BandInfo_d, BandInfo, sizeof(int)*3*incBandNoise, cudaMemcpyHostToDevice );
	checkCudaError(err);	 

	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid;
	
/////////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the Design Matrix////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	

	//printf("entered 2\n");
	if(D != NTime+NJumps){

		err = cudaMalloc( (void **)&DMatrix_d, sizeof(double)*N*D );
		checkCudaError(err);
		err = cudaMemcpy( DMatrix_d, TNDMVec, sizeof(double)*D*N, cudaMemcpyHostToDevice );
		checkCudaError(err);
	
		double *U_d;
		double *V_d;
		double *S_d;
	
		err = cudaMalloc( (void **)&U_d, sizeof(double)*N*N );
		checkCudaError(err);
		err = cudaMalloc( (void **)&V_d, sizeof(double)*D*D );
		checkCudaError(err);
		err = cudaMalloc( (void **)&S_d, sizeof(double)*D );
		checkCudaError(err);
	
	
		culaDeviceDgesvd('O','N', N, D, DMatrix_d, N, S_d, U_d, N, V_d, D);
	
		cudaFree(V_d);
		cudaFree(S_d);
		cudaFree(U_d);
		
		cudaDeviceSynchronize();	
	}
	
   

///////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the F Matrix/////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	
	
	//printf("entered 3\n");
	make_fmatrix<<< N, 1 >>>(GlobalTotalMatrix_d,Freqs_d,ObsFreqs_d, BATvec_d,DMVec_d,SysGroups_d, BandInfo_d, N,RF,DMF, BandNoiseCoeff, GroupNoiseCoeff, incRED, incDM, incBandNoise, incGroupNoise, 				ReplaceTMatrix, D, NTime,NJumps,DMatrix_d);
 	 
	cudaDeviceSynchronize();



///////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Add ECORR Matrix/////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	
/*

	if(incNGJitter > 0){

		int NGCoeffStartPoint = RF+DMF+incGroupNoise*GroupNoiseCoeff+incBandNoise*BandNoiseCoeff;

		add_EcorrToFMatrix<<< N, 1 >>>(FMatrix_d, GlobalEMatrix_d, N, NGCoeffStartPoint, numNGJitterEpochs);
	 	 
		cudaDeviceSynchronize();

	}

*/
/*

///////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Form the T Matrix/////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	 


	//printf("entered 4\n");
	dimGrid.x=(T + dimBlock.x - 1)/dimBlock.x;
	dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y; 

	if(D != NTime+NJumps){
		make_Tmatrix<<<N,1>>>(TMatrix_d, DMatrix_d, FMatrix_d, N, T, D, F);
	}
	else{

		make_Tmatrix<<<N,1>>>(TMatrix_d, GlobalStaticDmat_d, FMatrix_d, N, T, D, F);
	}

	cudaDeviceSynchronize();
*/
///////////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////Do the Algebra///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// 	
 	 
 	 
	//printf("entered 5 %i %i\n",T,N);
	
	
	dimGrid.x=(T + dimBlock.x - 1)/dimBlock.x;
	dimGrid.y = (N + dimBlock.y - 1)/dimBlock.y;	

	MatMulKernel<<<dimGrid, dimBlock>>>(N,N,N, T,Noise_d,GlobalTotalMatrix_d,NT_d);
	//SimpleDiagMatMulKernel<<<1,1>>>(N, T, Noise_d, TMatrix_d, NT_d);
	cudaDeviceSynchronize();

	double alpha=1.0;
	double beta=0.0; 
	

	status =  culaDeviceDgemm('T', 'N', T, T, N, alpha, GlobalTotalMatrix_d, N, NT_d, N, beta, TNT_d, T);
	checkStatus(status);

	cudaDeviceSynchronize();
	

	status = culaDeviceDgemv('T', N, T, alpha, NT_d, N, resvec_d, 1, beta, NTd_d, 1);
	checkStatus(status);
	
	cudaDeviceSynchronize();
	
	
	dimGrid.x=(T + dimBlock.x - 1)/dimBlock.x;
	dimGrid.y = (T + dimBlock.y - 1)/dimBlock.y;	
	

	addCoeffsKernel<<<1,1>>>(T,D,F,TNT_d,powercoeff_d);

	cudaDeviceSynchronize();
	
	
	//printf("entered 6: %i \n", T);
	double *dettemp_d;
	double *tempval=new double[1];
	err = cudaMalloc( (void **)&dettemp_d, sizeof(double) );
	checkCudaError(err);

	int carryOn=0;
	status = culaDeviceDpotrf('L', T, TNT_d, T);
	cudaDeviceSynchronize();
	//printf("entered 6.25\n");
	
	checkStatusCarryOn(status,carryOn);
	
	cudaDeviceSynchronize();

	//printf("entered 6.5\n");
	calc_det<<< 1, 1 >>>( TNT_d, dettemp_d, T);
	err = cudaMemcpy( tempval, dettemp_d, sizeof(double), cudaMemcpyDeviceToHost);
	checkCudaError(err);
	likevals[0]=tempval[0];
	
	cudaDeviceSynchronize();
     
	if(carryOn == 1){

		cudaFree(Freqs_d);
		cudaFree(ObsFreqs_d);
		cudaFree(powercoeff_d);
		cudaFree(resvec_d);
		cudaFree(BATvec_d);
		cudaFree(Noise_d);
		cudaFree(DMVec_d);
		cudaFree(SysGroups_d);

//		cudaFree(FMatrix_d);
		if(D != NTime+NJumps){cudaFree(DMatrix_d);}
//		cudaFree(TMatrix_d);
		cudaFree(NT_d);
		cudaFree(TNT_d);
		cudaFree(NTd_d);
		cudaFree(BandInfo_d);
		
		cudaFree(dettemp_d);
		delete(tempval);
		return;
	}

	//printf("entered 7\n");
	double *WorkVec_d;
	err = cudaMalloc( (void **)&WorkVec_d, sizeof(double)*T );
	checkCudaError(err);
	copyvec<<< T, 1 >>>(WorkVec_d, NTd_d, T);
	
	cudaDeviceSynchronize();
	

	status=culaDeviceDpotrs('L', T, 1, TNT_d, T, WorkVec_d, T);
	checkStatus(status);
	
	cudaDeviceSynchronize();
	
	double *freqlike_d;
	err = cudaMalloc( (void **)&freqlike_d, sizeof(double));
	checkCudaError(err);

	calc_DotLike<<< 1, 1 >>>(WorkVec_d, NTd_d, T, freqlike_d);
	
	cudaDeviceSynchronize();
	err = cudaMemcpy( tempval, freqlike_d, sizeof(double), cudaMemcpyDeviceToHost);
	checkCudaError(err);
	likevals[1]=tempval[0];

	
	//printf("entered 8\n");
	cudaFree(Freqs_d);
	cudaFree(ObsFreqs_d);
	cudaFree(powercoeff_d);
	cudaFree(resvec_d);
	cudaFree(BATvec_d);
	cudaFree(Noise_d);
	cudaFree(DMVec_d);	
	cudaFree(SysGroups_d);

//	cudaFree(FMatrix_d);
	if(D != NTime+NJumps){cudaFree(DMatrix_d);}
//	cudaFree(TMatrix_d);
	cudaFree(NT_d);
	cudaFree(TNT_d);
	cudaFree(NTd_d);
	cudaFree(BandInfo_d);
	
	cudaFree(WorkVec_d);
	cudaFree(dettemp_d);
	cudaFree(freqlike_d);
	delete(tempval);
	
	//printf("entered 9\n");
		
}
 	 


extern "C" void copy_staticTmat_(double *T, int totalsize, int Nobs){

    cudaError_t err;

 	 err = cudaMalloc( (void **)&GlobalTotalMatrix_d, sizeof(double)*totalsize*Nobs );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalTotalMatrix_d, T, sizeof(double)*totalsize*Nobs, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}


/*

 	 

 extern "C" void copy_floatgmat_(float *G, int N){

    cudaError_t err;

   // Allocate memory on GPU
	//printf("copying G\n");


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
	//printf("copying G\n");


 	 err = cudaMalloc( (void **)&GlobalGmat_d, sizeof(double)*N );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalGmat_d, G, sizeof(double)*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}

extern "C" void copy_staticgmat_(double *G, int M, int N){

    cudaError_t err;

   // Allocate memory on GPU
	//printf("copying G\n");


 	 err = cudaMalloc( (void **)&GlobalStaticGmat_d, sizeof(double)*N*M );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalStaticGmat_d, G, sizeof(double)*N*M, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}

extern "C" void copy_staticumat_(double *G, int M, int N){

    cudaError_t err;

   // Allocate memory on GPU
	//printf("copying G\n");


 	 err = cudaMalloc( (void **)&GlobalStaticUGmat_d, sizeof(double)*M*N );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalStaticUGmat_d, G, sizeof(double)*M*N, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}

extern "C" void copy_staticdmat_(double **TNDM, double *TNDMVec, int N, int D){

    cudaError_t err;

   // Allocate memory on GPU
	//printf("copying G\n");
	
	err = cudaMalloc( (void **)&GlobalStaticDmat_d, sizeof(double)*N*D );
	checkCudaError(err);

	err = cudaMemcpy(GlobalStaticDmat_d, TNDMVec, sizeof(double)*D*N, cudaMemcpyHostToDevice );
	checkCudaError(err);
	
	double *U_d;
	double *V_d;
	double *S_d;
	
	err = cudaMalloc( (void **)&U_d, sizeof(double)*N*N );
	checkCudaError(err);
	err = cudaMalloc( (void **)&V_d, sizeof(double)*D*D );
	checkCudaError(err);
	err = cudaMalloc( (void **)&S_d, sizeof(double)*D );
	checkCudaError(err);
	
	
	culaDeviceDgesvd('O','N', N, D, GlobalStaticDmat_d, N, S_d, U_d, N, V_d, D);
	
	cudaFree(V_d);
	cudaFree(S_d);
	cudaFree(U_d);
	
    cudaDeviceSynchronize();
    

}

extern "C" void copy_staticECorrmat_(double *E, int EcorrSize, int Nobs){

    cudaError_t err;

 	 err = cudaMalloc( (void **)&GlobalEMatrix_d, sizeof(double)*EcorrSize*Nobs );
	 checkCudaError(err);

     // copy vectors from CPU to GPU
 	 err = cudaMemcpy( GlobalEMatrix_d, E, sizeof(double)*EcorrSize*Nobs, cudaMemcpyHostToDevice );
 	 checkCudaError(err);

   return;
}

*/

