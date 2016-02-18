#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hanson.h"

int test_hanson_diag()
{
	const unsigned long ndim=10000;
	hanson_data* hdata;
	hdata=(hanson_data*)malloc(sizeof(hanson_data));

	init_hanson_diag_data(hdata,ndim,0.8,1.2);
	free_hanson_diag_data(hdata);
	free(hdata);
	return EXIT_SUCCESS;
}

int test_csv_IO()
{
	double x[10];
	double y[10];
	double z[10];	
	double mat[30];	
	int i,j;
	FILE* csvfile;
	
	/* write a csv file */
	csvfile=fopen("csv_file.dat","w");
	for(i=0;i<10;++i)
	{
		x[i]=(double)i;
		y[i]=(double)i*2;
		z[i]=(double)i*3;
	}
	write_csv_file(csvfile,10,x);
	write_csv_file(csvfile,10,y);
	write_csv_file(csvfile,10,z);
	
	fclose(csvfile);
	
	/* read the csv file */
	for(i=0;i<30;++i)
	{
		mat[i]=0;
	}
	csvfile=fopen("csv_file.dat","r");
	read_csv_file(csvfile,3,10,mat);
	
	/* test if they are equal */
	for(i=0;i<3;++i)
	{
		for(j=0;j<10;++j)
		{
			if(mat[i*(int)10+j] != (double)j*(i+1) )
				return EXIT_FAILURE;
		}
	}
	
	fclose(csvfile);
	
	
	return EXIT_SUCCESS;
}

