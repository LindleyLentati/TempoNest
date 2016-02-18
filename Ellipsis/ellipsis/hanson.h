#ifndef ELLIPSIS_HANSON_H
#define ELLIPSIS_HANSON_H

#include <stdio.h>

typedef struct
{
	int num_dims;
	int num_ents;
	int keep_sampling;
	double tol_min;
	double tol_max;
	long double* sum_x3dpsi;
	long double* sum_x2dpsi;
	long double* sum_xdpsi;
	long double* sum_dpsi;
	long double* sum_x;
	long double* sum_x2;
	double* hanson;
	double* start_point;
} hanson_data;

void write_csv_file(FILE* outfile,int num_dims,double* vect);
void write_LDcsv_file(FILE* outfile,int num_dims,long double* vect);
void read_csv_file(FILE* infile,int num_rows,
	int num_cols,double* matrix);

void init_hanson_diag_data(hanson_data* diag_data,int num_dims,
	double tol_min,double tol_max);
void free_hanson_diag_data(hanson_data* diag_data);
void push_state(hanson_data* diag_data,int num_dims,
	double* x, double* g);
void write_diag_data_to_file(hanson_data* diag_data,const char* file_name);
int read_diag_data_from_file(hanson_data* diag_data,const char* file_name);


#endif/*ELLIPSIS_HANSON_H*/
