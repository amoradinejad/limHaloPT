
/** @file header.h
 * 
 */

#ifndef HEADER_H
#define HEADER_H

#define _GNU_SOURCE

#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_odeiv2.h>    
#include <gsl/gsl_roots.h>	 
#include <gsl/gsl_sf_expint.h>
#include <ctype.h>
#include "cuba.h"


#define PSC  		101L
#define ST 		    	102L
#define TK 		    	103L

#define GROWTH		104L
#define DERGROWTH		105L

#define NONLINEAR	   	106L
#define LINEAR		107L

#define GAUSSIAN		114L
#define NONGAUSSIAN	115L

#define INIT  		116L
#define LOCAL		117L
#define EQUILATERAL	118L
#define ORTHOGONAL	119L
#define QSF			120L
#define HS			121L
#define NGLOOP		122L
#define derNGLOOP		123L

#define QUADRATIC		124L
#define TIDE		125L
#define GAMMA		126L

#define LPOWER		127L
#define NLPOWER		128L
#define CDM 		90L
#define BA 			91L
#define TOT 		92L

#define CO10  		131L
#define CO21  		132L
#define CO32  		133L
#define CO43  		134L
#define CO54  		135L
#define CO65  		136L
#define CII 		137L

#define MATTER		138L
#define LINEMATTER	139L
#define LINE 		140L

#define DST			141L	
#define GFILTER		142L
#define BSPLINE		143L


#define TREE		144L	
#define LOOP		145L
#define WIR			146L
#define NOIR		147L

#define HALO 		148L


#define PS_KMIN   	1.e-5
#define PS_KMAX   	100.1 
#define PS_ZMAX   	14.

#define CLEANUP        1 
#define DO_NOT_EVALUATE -1.0
#define NPARS     6
#define MAXL 2000


/** 
 * List of limHaloPT header files
 */
/// \cond DO_NOT_DOCUMENT
#include "../Class/include/class.h"
/// \endcond
#include "cubature.h"
#include "Global_Structs.h"
#include "utilities.h"
#include "cosmology.h"
#include "survey_specs.h"
#include "line_ingredients.h"	
#include "wnw_split.h"
#include "IR_res.h"
#include "ps_halo_1loop.h"
#include "ps_line_pt.h"
#include "ps_line_hm.h"


/** 
 * Function declarations of main.c module
 */
void  initialize();
void 	cleanup();


/**
 * A structure passed to the integrators to hold the parameters fixed in the integration
 */
struct integrand_parameters
{
	double p1;
	double p2;
	double p3;
	double p4;
	double p5;
	double p6;
	double p7;
	double p8;
	double p9;
	double p10;
	double p11;
	long 	p12;
	long 	p13;
};

/**
 * Another structure passed to the integrators to hold the parameters fixed in the integration
 */
struct integrand_parameters2
{

	struct Cosmology *p1;
	struct Cosmology *p2;
	struct Cosmology *p3;

	double p4;
	double p5;
	double p6; 
	double p7;
	double p8;
	double p9;
	double p10;
	double p11;
	double p12;

	long p13;
	long p14;
	long p15;
	long p16;
	long p17;
	long p18;

	int p19;

	double *p20;
	size_t p22;
};
	

#endif







