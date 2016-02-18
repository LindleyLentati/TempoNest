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

#ifndef ELLIPSIS_MT19937_H
#define ELLIPSIS_MT19937_H

/* Period parameters */  
#define MT19937_N 624
#define MT19937_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

typedef struct
{
	unsigned long mt[MT19937_N];
	int mti;
} ellipsis_mt19937_rng;

/* initializes mt[MT19937_N] with a seed */
void init_genrand(ellipsis_mt19937_rng* rng,unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(ellipsis_mt19937_rng* rng,unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(ellipsis_mt19937_rng* rng);

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(ellipsis_mt19937_rng* rng);

/* generates a random number on [0,1]-real-interval */
double genrand_real1(ellipsis_mt19937_rng* rng);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(ellipsis_mt19937_rng* rng);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(ellipsis_mt19937_rng* rng);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(ellipsis_mt19937_rng* rng);
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* added by Sree on 29-09-2012 */
/*-----------------------------*/

/* generate uniform [0,1] random numbers */
double genrand_uniform(ellipsis_mt19937_rng* rng);

/* generate normal [0,1] random numbers */
double gerand_gauss(ellipsis_mt19937_rng* rng);

/* save random number state */
void save_rand_state(ellipsis_mt19937_rng* rng,const char* filename);

/* read random number state */
void read_rand_state(ellipsis_mt19937_rng* rng,const char* filename);

#endif/*ELLIPSIS_MT19937_H*/
