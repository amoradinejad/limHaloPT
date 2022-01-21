#ifndef _SURVEY_SPECS_H_
#define _SURVEY_SPECS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int *make_1D_int_array(long size);
double 	*make_1Darray(long size);
double  **make_2Darray(long nrows, long ncolumns);
double 	*init_1Darray(long n,double xmin,double xmax);
double 	*loginit_1Darray(long n,double xmin,double xmax);
long 	count_lines_in_file(char *fname);
void 	return_arr(double a, double b, long n,double ** arr1);
void 	free_2Darray(double ** m);
	long count_cols_in_file(char *fname);

double *log10init_1Darray(long n, double inc, double xmin);


#endif