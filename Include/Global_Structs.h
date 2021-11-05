
/** @file Global_Structs.h
 * 
 */
#ifndef GLOBALSTRUCTS_H_
#define GLOBALSTRUCTS_H_


/**
 * Structure to store cosmology structure from CLASS code
*/
struct Class_Cosmology_Struct{

    struct precision                    pr;             /* for precision parameters */
    struct background                   ba;             /* for cosmological background */
    struct thermo                       th;             /* for thermodynamics */
    struct perturbs                     pt;             /* for source functions */
    struct transfers                    tr;             /* for transfer functions */
    struct primordial                   pm;             /* for primordial spectra */
    struct spectra                      sp;             /* for output spectra */
    struct nonlinear                    nl;             /* for non-linear spectra */
    struct lensing                      le;             /* for lensed spectra */
    struct output                       op;             /* for output files */
    ErrorMsg errmsg;                    /* for error messages */
};


/**
 * Structure that holds varioud quantities that need to be evaluated for a given choice of cosmological 
 * paramteres. This includes, the Class_Cosmology_Struct (initialized in cosmology.c), and 
 * Line Structure (initialized in line_ingredients.c). 
*/
struct Cosmology
{

     struct Class_Cosmology_Struct    ccs;
     struct Line                      **Lines;

     int                              NLines;
     long                             mode_nu;

     double cosmo_pars[6];
};



/**
 * Structure that holds the Line-related quantities, including the interpolators for first and second moments of the line luminosity 
 * and the linear and quadratic luminosity-weighted line biases.
*/
struct Line
{
      long                    LineType;
      int                     initialized;
      size_t                  npointsInterp; 

      double                  line_freq;                                         

      gsl_interp_accel        *mom1_accel_ptr;
      gsl_spline              *mom1_spline_ptr;
      gsl_interp_accel        *mom2_accel_ptr;
      gsl_spline              *mom2_spline_ptr;
      
      gsl_interp_accel        *b1_LW_accel_ptr;
      gsl_spline              *b1_LW_spline_ptr;
      gsl_interp_accel        *b2_LW_accel_ptr;
      gsl_spline              *b2_LW_spline_ptr;

};


/**
 * A global structure including the values of cosmological parmaeters, 2d interpolator of SFR, and names of various files.
*/
struct globals 
{
	double H0;
	double c;
	
	double As;
	double logAs;
	double ns;
	double h;
	double Omega_cdm;
	double Omega_b;
	double Omega_r;
	double Omega_lambda;
	double Omega_g;
	double Omega_nu;

	double b1;
	double sigFOG0;
	
	long Npars;
	double z_i;
	double rho;
	double mass;
	double kp;
	double ng;
	double volume;

	double kf;
	double h_m;

	double M_min;
	double M_max;
	double z_max;

	char project_home[FILENAME_MAX];
	char output_dir[FILENAME_MAX];
	char data_dir[FILENAME_MAX];
	char data_priors[FILENAME_MAX];


	// Min and max values
	double 				PS_kmin;
	double 				PS_kmax;

	// File names
	char			SFR_filename[FILENAME_MAX];
	char 			Planck_Fisher_filename[FILENAME_MAX];

	gsl_interp_accel    *logM_accel_ptr;
  	gsl_interp_accel    *z_accel_ptr;
  	gsl_spline2d        *logSFR_spline2d_ptr;	
};

#endif
 


















