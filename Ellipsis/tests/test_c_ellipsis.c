#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "test_mt19937.h"
#include "test_hanson.h"
#include "test_kinetic_energy.h"

int main(void)
{
	int ret=EXIT_SUCCESS;
	
	ret+=test_random_uni();
	ret+=test_random_norm();
	ret+=test_rand_save_state();	
	ret+=test_csv_IO();	
	ret+=test_hanson_diag();
	ret+=test_kinetic_energy();
	
	return ret;
}

