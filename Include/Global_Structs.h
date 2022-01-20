
/** @file Global_Structs.h
 * 
 */
#ifndef GLOBALSTRUCTS_H_
#define GLOBALSTRUCTS_H_


/**
 * Structure to store cosmology structure from CLASS-v3.1 
 */
struct Class_Cosmology_Struct{ 

    struct precision                    pr;             /* for precision parameters */
    struct background                   ba;             /* for cosmological background */
    struct thermodynamics               th;             /* for thermodynamics */ 
    struct perturbations                pt;             /* for source functions */
    struct transfer                     tr;             /* for transfer functions */
    struct primordial                   pm;             /* for primordial spectra */
    struct harmonic                     hr;             /* for output spectra */
    struct fourier                      fo;             /* for non-linear spectra */
    struct lensing                      le;             /* for lensed spectra */
    struct distortions					sd;				/* for CMB distortions*/
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
    struct Line             **Lines;
    struct PS_xtr		    *PS_xtrapol;

    int                     NLines;
    long                    mode_model;

    double 			        cosmo_pars[NPARS];

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

struct PS_xtr
{
    int                     initialized;
	
    gsl_interp_accel        *z_accel_ptr;
    gsl_interp_accel        *logk_accel_ptr;
    gsl_spline2d            *logpkz_spline2d_ptr;
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

	long Npars; 
	
	char project_home[FILENAME_MAX];
	char output_dir[FILENAME_MAX];
	char data_dir[FILENAME_MAX];
	char data_priors[FILENAME_MAX];


	// Min and max values
	double  PS_kmin;
	double  PS_kmax;
	double  kmax_CLASS;

	// File names
	char			SFR_filename[FILENAME_MAX];
	char 			Planck_Fisher_filename[FILENAME_MAX];

	gsl_interp_accel    *logM_accel_ptr;
  	gsl_interp_accel    *z_accel_ptr;
  	gsl_spline2d        *logSFR_spline2d_ptr;	
};

#endif
 


















