
#ifndef _TEST_H_
#define _TEST_H_

#define _GNU_SOURCE

#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <omp.h>
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

#define PSC             101L
#define ST              102L
#define TK              103L

#define CO10            131L
#define CO21            132L
#define CO32            133L
#define CO43            134L
#define CO54            135L
#define CO65            136L
#define CII             137L

#define NPARS     6

extern struct globals gb;

#include "../Class/include/class.h"
#include "global_structs.h"
#include "cosmology.h"

typedef struct initialize_struct{
      double Mh_min;
      long   mode_mf;
      size_t ninterp;
      int    nlines;
      int    lines[7];
}initialize_struct;

void initialize(char *argv[], struct initialize_struct *init_struct);
void  cleanup();

int Cosmology_init(struct Cosmology *Cx, double pk_kmax, double pk_zmax, 
            int nlines, int * line_types, size_t npoints_interp, double M_min, long mode_mf);

double PS_line_HM(struct Cosmology *Cx, double k, double z, double M_min, long mode_mf,  long line_type, int line_id);
double PS_shot_HM(struct Cosmology *Cx, double k, double z, double M_min, double *input, long mode_mf, long line_type, int line_id);
double b22_ls(struct Cosmology *Cx, double z);

double PS_shot(struct Cosmology *Cx, double z, size_t line_id);
double Tbar_line(struct Cosmology *Cx, size_t line_id, double z);
void   line_bias(struct Line *Lx, double z, double *result);

double rhoc(struct Cosmology *Cx, double z);

#endif


