
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


struct Cosmology
{

     struct Class_Cosmology_Struct    ccs;
     struct Line                      **Lines;

     int                              NLines;
     long                             mode_nu;

     double cosmo_pars[6];
};



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

	double fnl;
	double nu;

	double b1;
	double sigFOG0;

	double ns_min;
	double ns_max;

	double h_min ;
	double h_max;

	double Omega_cdm_min;
	double Omega_cdm_max;

	double Omega_b_min ;
	double Omega_b_max;
	
	double logAs_min;
	double logAs_max;

	double fnl_min;
	double fnl_max;
	
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


 


















