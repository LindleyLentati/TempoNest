/* Routines based around the TEMPO2 use of the Cholesky algorithm    */
/* These have been put together by G. Hobbs, W. Coles and R. Shannon */

#include "cholesky.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
extern "C"{
#include "quadmath.h"
}


void TNcholDecompositionL(__float128 **a, int n, double &det, int &info)
{
   int i,j,k;
   __float128 sum;

   __float128 *p;
   p = (__float128 *) malloc(sizeof(__float128)*(n+1));
    
   info = 0;
   for (i=0;i<n;i++)
   {
	  for (j=i;j<n;j++)
	  {
	   
	    for (sum=a[i][j],k=i-1;k>=0;k--)
	      {
		sum-=a[i][k]*a[j][k]; 
		//	printf("%d %d %d %Lg %Lg %.3Le\n",i,j,k,a[i][k],a[j][k],sum);
	      }
		 if (i==j)
		 {
			//	      printf("Currently have %d %d %Lg\n",i,j,sum);
		   if (sum <= 0.Q)
			{
//			   printf("Here with %d %d %Qg %Qg\n",i,j,a[i][j],sum);
			   for (sum=a[i][j],k=i-1;k>=0;k--)
			   {
				  sum-=a[i][k]*a[j][k];
				  printf("Failed: %d %d %d %Qg %Qg %.3Qe\n",i,j,k,a[i][k],a[j][k],sum);
			   }
			   printf("Failed in long decomp- the matrix is not positive definite\n");
//			   exit(1);
			  	info = 1;
			}
			p[i] = sqrtq(sum);
			//	      printf("Currently have %d %d %Lg %g\n",i,j,sum,p[i]);
		 }
		 else
		 {
			a[j][i] = sum/p[i];
			/*	      if (j==120)
					  printf("j=120, setting %g %Lg %g %d\n",a[j][i],sum,p[i],i);
					  if (j==130)
					  printf("j=130, setting %g %Lg %g %d\n",a[j][i],sum,p[i],i);*/
		 }
	  }
   }

	det = 0;
	for(int i =0; i < n; i ++){
		a[i][i] = p[i];
		det += 2*log(p[i]);
	}
	free(p);
	
}

void T2cholDecompositionL(__float128 **a, int n, __float128 *p, int &info)
{
   int i,j,k;
   __float128 sum;


    info = 0;
   for (i=0;i<n;i++)
   {
	  for (j=i;j<n;j++)
	  {
	   
	    for (sum=a[i][j],k=i-1;k>=0;k--)
	      {
		sum-=a[i][k]*a[j][k]; 
		//	printf("%d %d %d %Lg %Lg %.3Le\n",i,j,k,a[i][k],a[j][k],sum);
	      }
		 if (i==j)
		 {
			//	      printf("Currently have %d %d %Lg\n",i,j,sum);
		   if (sum <= 0.Q)
			{
//			   printf("Here with %d %d %Qg %Qg\n",i,j,a[i][j],sum);
			   for (sum=a[i][j],k=i-1;k>=0;k--)
			   {
				  sum-=a[i][k]*a[j][k];
				  printf("Failed: %d %d %d %Qg %Qg %.3Qe\n",i,j,k,a[i][k],a[j][k],sum);
			   }
			   printf("Failed in long decomp- the matrix is not positive definite\n");
//			   exit(1);
			  	info = 1;
			}
			p[i] = sqrtq(sum);
			//	      printf("Currently have %d %d %Lg %g\n",i,j,sum,p[i]);
		 }
		 else
		 {
			a[j][i] = sum/p[i];
			/*	      if (j==120)
					  printf("j=120, setting %g %Lg %g %d\n",a[j][i],sum,p[i],i);
					  if (j==130)
					  printf("j=130, setting %g %Lg %g %d\n",a[j][i],sum,p[i],i);*/
		 }
	  }
   }
}

/**
 * UINV is a lower triangluar matrix.
 * Matricies are row-major order, i.e. uinv[r][c].
 *
 */
void cholesky_formUinvL(double **uinv,__float128** mm,int np, double &det, int &info){
   int i,j,k;
   //long double **mm;
   __float128 *cholpp;
   //double *cholp;
   
     

   cholpp = (__float128 *) malloc(sizeof(__float128)*(np+1));
   
  
   __float128 sum;
   
   //double *cholp  = (double *)malloc(sizeof(double)*(np+1));
   
   
   

     
   T2cholDecompositionL(mm,np,cholpp, info);

        det=0;
	for(int i=0;i<np;i++){
		det+=log(cholpp[i]);
	}

	det=det*2;

   for (i=0;i<np;i++)
     {
       mm[i][i] = 1.Q/cholpp[i];
       uinv[i][i] = (double) mm[i][i];
       for (j=0;j<i;j++)
	 uinv[j][i] = 0.0;
       for (j=i+1;j<np;j++)
	 {
	   sum=0.Q;
	   for (k=i;k<j;k++) sum-=mm[j][k]*mm[k][i];
	   mm[j][i]=sum/cholpp[j];
	   uinv[j][i] = (double) mm[j][i];
	   //uinv[i][j]=uinv[j][i];
	 }
     }




   
   free(cholpp);

   return;
}


