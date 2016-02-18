/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mt19937.h"

/* initializes mt[MT19937_N] with a seed */
void init_genrand(ellipsis_mt19937_rng* rng,unsigned long s)
{
    rng->mti=MT19937_N+1;
    rng->mt[0]= s & 0xffffffffUL;
    for (rng->mti=1; rng->mti<MT19937_N; rng->mti++) {
        rng->mt[rng->mti] = 
	    (1812433253UL * (rng->mt[rng->mti-1] ^ (rng->mt[rng->mti-1] >> 30)) + rng->mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        rng->mt[rng->mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(ellipsis_mt19937_rng* rng,unsigned long init_key[], int key_length)
{
    int i, j, k;
    rng->mti=MT19937_N+1;
    init_genrand(rng,19650218UL);
    i=1; j=0;
    k = (MT19937_N>key_length ? MT19937_N : key_length);
    for (; k; k--) {
        rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^ (rng->mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        rng->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=MT19937_N) { rng->mt[0] = rng->mt[MT19937_N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=MT19937_N-1; k; k--) {
        rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^ (rng->mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        rng->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=MT19937_N) { rng->mt[0] = rng->mt[MT19937_N-1]; i=1; }
    }

    rng->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(ellipsis_mt19937_rng* rng)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (rng->mti >= MT19937_N) { /* generate MT19937_N words at one time */
        int kk;

        if (rng->mti == MT19937_N+1)   /* if init_genrand() has not been called, */
            init_genrand(rng,5489UL); /* a default initial seed is used */

        for (kk=0;kk<MT19937_N-MT19937_M;kk++) {
            y = (rng->mt[kk]&UPPER_MASK)|(rng->mt[kk+1]&LOWER_MASK);
            rng->mt[kk] = rng->mt[kk+MT19937_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MT19937_N-1;kk++) {
            y = (rng->mt[kk]&UPPER_MASK)|(rng->mt[kk+1]&LOWER_MASK);
            rng->mt[kk] = rng->mt[kk+(MT19937_M-MT19937_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (rng->mt[MT19937_N-1]&UPPER_MASK)|(rng->mt[0]&LOWER_MASK);
        rng->mt[MT19937_N-1] = rng->mt[MT19937_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        rng->mti = 0;
    }
  
    y = rng->mt[rng->mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(ellipsis_mt19937_rng* rng)
{
    return (long)(genrand_int32(rng)>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(ellipsis_mt19937_rng* rng)
{
    return genrand_int32(rng)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(ellipsis_mt19937_rng* rng)
{
    return genrand_int32(rng)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(ellipsis_mt19937_rng* rng)
{
    return (((double)genrand_int32(rng)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(ellipsis_mt19937_rng* rng) 
{ 
    unsigned long a=genrand_int32(rng)>>5, b=genrand_int32(rng)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* added by Sree on 29-09-2012 */
/*-----------------------------*/

/* generate uniform [0,1] random numbers */
double genrand_uniform(ellipsis_mt19937_rng* rng)
{
	/* same sas genrand_real1(ellipsis_mt19937_rng* rng) */
	return genrand_int32(rng)*(1.0/4294967295.0); 
	/* divided by 2^32-1 */ 
}

/* generate normal [0,1] random numbers */
double gerand_gauss(ellipsis_mt19937_rng* rng)
{
	double u,v,r,c;
	u=2.*genrand_real1(rng)-1.;
	v=2.*genrand_real1(rng)-1.;
	r=u*u+v*v;
	if (r==0||r>1)
	{
		return gerand_gauss(rng);
	}
	c=sqrt(-2.*log(r)/r);
	return u*c;
}

/* save random number state */
void save_rand_state(ellipsis_mt19937_rng* rng,const char* filename)
{
	FILE *outfile;
	unsigned i;
	outfile=fopen(filename,"w");
	
	for(i=0;i<MT19937_N;++i)
	{
		fprintf(outfile,"%lu\n",rng->mt[i]);
	}
	fprintf(outfile,"%d\n",rng->mti);
	fclose(outfile);
}

/* read random number state */
void read_rand_state(ellipsis_mt19937_rng* rng,const char* filename)
{
	FILE *outfile;
	unsigned i,info=0;
	outfile=fopen(filename,"r");
	
	for(i=0;i<MT19937_N;++i)
	{
		info=fscanf(outfile,"%lu",&(rng->mt[i]));
	}
	info=fscanf(outfile,"%d",&(rng->mti));
	
	if(!info) printf("\nERROR in reading random number state\n");
	
	fclose(outfile);
}

