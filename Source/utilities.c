
/** @file utilities.c Documented basic utility functions used by other modules of the code. 
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 *  
 * In summary, the following functions can be called from other modules:
 * -# make_1Darray()  			dynamically allocates memory to a 1d array
 * -# make_2Darray()  			dynamically allocates memory to a 2d array
 * -# free_2Darray()  			free the memory allocated to a 2d array      
 * -# init_1Darray()  			initialize a 1d array with linear spacing
 * -# loginit_1Darray()         initialize a 1d array with natural-log spacing
 * -# log10init_1Darray()       initialize a 1d array with log10 spacing
 * -# count_lines_in_file()     count the number of lines of a file
 * -# count_cols_in_file()      count number of columns of a file
 * -# return_arr()              
 */


#include "header.h"



/**
 * Allocate memory to a 1d array of type double and length size
 * 
 * @param size      Input: length of the array 
 * @return a pointer to a 1d array
 */
double *make_1Darray(long size)
{
	double *array;
	long i;

	array = (double*) malloc(size*sizeof(double));
	if (array==NULL) exit(2);
	for(i=0L;i<size;i++) array[i]= 0.0;
	return array;	

}

/**
 * Allocate memory to a 1d array of type integer and length size
 * 
 * @param size      Input: length of the arrat 
 * @return a pointer to an integer type 1d array
 */
int *make_1D_int_array(long size)
{
	int *array;
	long i;

	array = (int*) malloc(size*sizeof(int));
	if (array==NULL) exit(2);
	for(i=0L;i<size;i++) array[i]= 0.0;
	
	return array;	

}

/**
 * Allocate memory to a 2d array of type double 
 * 
 * @param nrows      Input: number of rows of the output array
 * @param ncols      Input: number of columns of the output array
 * @return a double pointer to a double type 2d array
 */
 double **make_2Darray(long nrows, long ncols)
{   
	long size;
	long i,j;
	double **array;

	size = nrows*ncols;

	array = (double **) malloc(nrows*sizeof(double *));
	if(array==NULL) exit(2);

	array[0] = (double *) malloc(size*sizeof(double));
	if(array[0]==NULL) exit(2);

	for(i=1L;i<nrows;i++)
	{
		array[i] = array[0] + i * ncols;
	}


	for(i=0L; i<nrows; i++) 
	{
		for(j=0L;j<ncols;j++)
		{
			array[i][j] = 0.0;
		}
	}

	return array;
}


/**
 * Free the memory allocated to a 2d array
 * 
 * @param m      Input: double pointer to the elements of 2d array
 * 
 * @return void
 */
void free_2Darray(double ** m)
{

	free(m[0]);
	free(m);

	return;
}


/**
 * initialize a 1d array, with values in the range of [xmin,xmax] and evenely-space on linear scale
 * 
 * @param n     	Input: number of elements
 * @param xmin      Input: start point
 * @param xmax     Input: end point
 * 
 * @return a pointer to a double type 1d array, with values initialized
 */
double *init_1Darray(long n,double xmin,double xmax)
{
	double inc;
	double *A;
	long i;

	A = make_1Darray(n+1L);

	inc = (xmax-xmin)/n;

	for(i=0L;i<n+1L;i++)
	{
		A[i] = xmin + inc * i;
		//printf("%12.6e\n",A[i]);
	}

	return A;
	}


/**
 * initialize a 1d array, with values in the range of [xmin,xmax] and evenely-space on natural-log scale
 * 
 * @param n     	Input: number of elements
 * @param xmin      Input: start point
 * @param xmax     Input: end point
 * @return a pointer to a double type 1d array, with values initialized
 */
double *loginit_1Darray(long n,double xmin,double xmax)
{

	double inc;
	double *A;
	long i;

	A = make_1Darray(n);
	
	xmin = log(xmin);
	xmax = log(xmax);
	inc  = (xmax-xmin)/n;

	for(i=0L;i<n;i++)
	{
		A[i] = xmin + inc * i;
		A[i] = exp(A[i]); 
	}

	return A;
}


/**
 * initialize a 1d array, with values in the range of [xmin,xmax] and evenely-space on log10 scale
 * 
 * @param n     	  Input: number of elements
 * @param inc       Input: increment
 * @param xmin      Input: start point
 * @return a pointer to a double type 1d array, with values initialized
 */
double *log10init_1Darray(long n, double inc, double xmin)
{
	
	double *A, *logA;
	long i;

	double logxmin;
	logA = make_1Darray(n);
	A = make_1Darray(n);

	logxmin = log10(xmin);

	for(i=0L;i<n;i++)
	{
		logA[i] = logxmin + inc * i;
		A[i] = pow(10,logA[i]); 
	}

	free(logA);

	return A;
}



/**
 * Count the number of lines of a file
 * 
 * @param fname     	Input: filename
 * @return long integer value of nlines
 */
long count_lines_in_file(char *fname)
{
	FILE* ifp;
	long numlines = 0L;
	char c;
	int verbose = 1;
	char line[MAXL];
	
	ifp = fopen(fname, "r");
	if(ifp == NULL)
	{
    	fprintf( stderr, "Unable to open file \n" );
    	exit( 1 );
	}

	int numcoms = 0;
	while((c=fgetc(ifp)) != EOF)
	{ 
			if(c =='\n')
       			numlines++;
     
       		if(c =='#')
       			numcoms++;
       		    	
	}
	if (verbose > 3) printf("Number of lines: %ld Number of comment lines: %ld\n",numlines,numcoms);		
	
	int data_len = 	numlines-numcoms;

	fclose(ifp);

	return data_len;
}


/**
 * Count the number of columns of a file
 * 
 * @param fname     	Input: filename
 * @return long integer value of ncols
 */
long count_cols_in_file(char *fname)
{
	FILE* ifp;
	long nrows = 0L, ncols =1;
	char c;
	int verbose = 1;
	char line[MAXL];
	
	ifp = fopen(fname, "r");
	if(ifp == NULL)
	{
    	fprintf( stderr, "Unable to open file \n" );
    	exit( 1 );
	}

	int count = 0;
	while(fgets(line, sizeof line, ifp) )
	{ 
	 	if(*line == '#')  continue; 
	 	while((c=fgetc(ifp)) != EOF)
		{ 
			if(c =='\n')
	       		nrows++;
			if (nrows == 2 && c== '\t')
				ncols ++;
		}  
	} 	

	if (verbose > 3) printf("num_rows: %ld \t num_cols =%ld\n",nrows,ncols);
	fclose(ifp);
		
	return ncols;
}

