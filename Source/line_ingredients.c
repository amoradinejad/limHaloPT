
/** @file line_ingredients.c Documented line_ingredients module
 * 
 * This module includes functions that are needed for computing the line clustering and shot contributions.
 *  
 * Azadeh Moradinezhad Dizgah, September 6th 2021
 *
 * In summary, the following functions can be called from other modules:
 * -# Line_alloc_init()       allocate memory and initizlized the line structure which contains 4 interpolators for first and second mass moments and linear and quadratic line biases.
 * -# Line_free()             frees the memory allocated to line structure
 * -# Line_evaluate()         evaluates the interpolators initialized in Line_alloc_init()
 * -# mult_func()             computes the multiplicity function needed for computing the halo mass function
 * -# mass_func()             computes the halo mass finction. Three options are available, Press-Schecter, Sheth-Tormen, Tinker
 * -# mass_func_sims()        reads in the measured mass function on Hidden-Valley simulations by Farnik, and convert it to compare with the theoretical predictions 
 * -# halo_bias()             computes the halo biases assuming the above theoretical predictions of the halo mass function
 * -# logSFR_Behroozi_read()  reeds in the data file of Behroozi 2013 for SFR(M,z)
 * -# logSFR_alloc_init()     allocates memory for 2d interpolator of logSFR(M,z)
 * -# SFR_behroozi_free()     frees the memory allocated to logSFR interpolator
 * -# logSFR_Behroozi()       evaluates the logSFR_Behroozi interpolator
 * -# luminosity()            computes the line luminosity
 * -# mass_moment1()          computes the first mass moment
 * -# mass_moment2()          computes the first mass moment
 * -# bias_lum_weighted()     computes the luminosity-weighetd line bias
 * -# p_sig_shot()            computes the coefficient accounting for the scatter in L(M) in shot noise
 * -# p_sig_Tbar()            computes the coefficient accounting for the scatter in L(M) in mean brightness temprature
 * -# mean_intens()           compues the mean intensity of the line
 * -# Tbar_line()             compues the mean brightness temprature of the line
 * 
 */


#include "header.h"
struct globals gb;


/**
 * Allocate the memory and initialize the the line structure. This structure contains interpolators for computing the luminosity-weighted mass moments and line biases
 * For a given line defined with "line_type" variable, this function first computes the above four quantities for a wide range of redshifts. 
 * Next it iniialized 4 interpolators for these quantities, and store them in line structure.
 * 
 * @param Cx                Input: Cosmology structure
 * @param line_type         Inpute: name of the line to compute.
 *                              It can be set to CII, CO10, CO21, CO32, CO43, CO54, CO65 
 * @param npoints_interp    Input: number of interpolation points 
 * @param M_min             Input: minimum halo mass for mass integrals
 * @param mode_mf           Inpute: theoretical model of halo mass function to use. 
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * @return the total clustering line power spectrum, including the 1- and 2-halo term         
 */

struct Line * Line_alloc_init(struct Cosmology *Cx, long line_type, size_t npoints_interp, double M_min, long mode_mf) 
{

  struct Line *   ThisLine;

  // Allocation of the Line object
  ThisLine  = (struct Line *)malloc(sizeof(struct Line));

  // Attributes
  ThisLine->LineType  = line_type;

  int J=0;
  if(line_type == CII){
    ThisLine->line_freq = 1902. * pow(10.,9.); ///CII
    J = 0;
  }
  else
  { 
    if(line_type == CO10) 
      J = 1;
    else if(line_type == CO21) 
      J = 2;
    else if(line_type == CO32) 
      J = 3;
    else if(line_type == CO43)
      J = 4; 
    else if(line_type == CO54)
      J = 5; 
    else if(line_type == CO65)
      J = 6; 
    ThisLine->line_freq = J * 115.27 * pow(10.,9.);  ////in unit of Hz for CO
  }
  
  // Allocation of interpolation objects
  ThisLine -> mom1_accel_ptr      = gsl_interp_accel_alloc();
  ThisLine -> mom1_spline_ptr     = gsl_spline_alloc(gsl_interp_linear, npoints_interp);
  ThisLine -> mom2_accel_ptr      = gsl_interp_accel_alloc();
  ThisLine -> mom2_spline_ptr     = gsl_spline_alloc(gsl_interp_linear, npoints_interp);
  
  ThisLine -> b1_LW_accel_ptr     = gsl_interp_accel_alloc();
  ThisLine -> b1_LW_spline_ptr    = gsl_spline_alloc(gsl_interp_linear, npoints_interp);
  ThisLine -> b2_LW_accel_ptr     = gsl_interp_accel_alloc();
  ThisLine -> b2_LW_spline_ptr    = gsl_spline_alloc(gsl_interp_linear, npoints_interp); 

  // Initialization of the interpolating objects
  double *zz, *mom1, *mom2, *b1_LW, *b2_LW;
  double zmin =0.0, zmax = 13.;
  // Allocate temporary arrays 
  mom1    = make_1Darray(npoints_interp);
  mom2    = make_1Darray(npoints_interp);
  b1_LW   = make_1Darray(npoints_interp);
  b2_LW   = make_1Darray(npoints_interp);
  zz      = init_1Darray(npoints_interp,zmin,zmax);

  double bias_arr[npoints_interp][2];

  // Calculating values to be interpolated
  #pragma omp parallel
  #pragma omp for
  for(int i=0;i<npoints_interp;i++){
    mom1[i]  =  mass_moment1(Cx,zz[i],M_min,mode_mf,line_type);
    mom2[i]  =  mass_moment2(Cx,zz[i],M_min,mode_mf,line_type);
    bias_lum_weighted(Cx,zz[i],M_min,mode_mf,line_type,bias_arr[i]); 
    b1_LW[i] = bias_arr[i][0];
    b2_LW[i] = bias_arr[i][1];
    printf("hm %ld %d %12.6e %12.6e %12.6e %12.6e %12.6e \n", line_type, i, zz[i], mom1[i], mom2[i], b1_LW[i], b2_LW[i]);
    // printf("process %d done computing redshift %12.6e\n", omp_get_thread_num(),zz[i]);
  }              

  // Initialize interpolating object
  gsl_spline_init(ThisLine -> mom1_spline_ptr, zz, mom1, npoints_interp);
  gsl_spline_init(ThisLine -> mom2_spline_ptr, zz, mom2, npoints_interp);
  gsl_spline_init(ThisLine -> b1_LW_spline_ptr, zz, b1_LW,npoints_interp);
  gsl_spline_init(ThisLine -> b2_LW_spline_ptr, zz, b2_LW,npoints_interp);

  printf("for line %d mom1, mom2, b1 and b2 interpolators initialized\n", line_type); 
  // Freeing temporary arrays
  free(zz);
  free(mom1);
  free(mom2);
  free(b1_LW);
  free(b2_LW);

  ThisLine->initialized = _TRUE_;

  return ThisLine;

}


/**
 * Free the line structure 
 *  
 * @param Lx     Input: Pointer to line structure
 * @return the error status
 */

int Line_free(struct Line * Lx)
{
  gsl_interp_accel_free(Lx->mom1_accel_ptr);
  gsl_spline_free(Lx->mom1_spline_ptr);
  gsl_interp_accel_free(Lx->mom2_accel_ptr);
  gsl_spline_free(Lx->mom2_spline_ptr);
  gsl_interp_accel_free(Lx->b1_LW_accel_ptr);
  gsl_spline_free(Lx->b1_LW_spline_ptr);
  gsl_interp_accel_free(Lx->b2_LW_accel_ptr);
  gsl_spline_free(Lx->b2_LW_spline_ptr);

  free(Lx);

  return _SUCCESS_;
}


/**
 * Allocate the memory and initialize the the line structure. This structure contains interpolators for computing the luminosity-weighted mass moments and line biases
 * For a given line defined with "line_type" variable, this function first computes the above four quantities for a wide range of redshifts. 
 * Next it iniialized 4 interpolators for these quantities, and store them in line structure.
 * 
 * @param Lx                Input: Pointer to the line structure
 * @param zz                Input: this is an array with 4 elements to determine which of the 4 interpolators should be evaluated. 
 *                                 - If any of the elements are set to DO_NOT_EVALUATE, the quantitiy corresponding to that index is not computed. O
 *                                 - If any of the elements is set to z, the corresponding quantity would be evaluated at that redshift
 * @param res               Output: an array containing the results. The number of elements of this array depends on how the zz array is set.
 * @return the error status
 */
int Line_evaluate(struct Line * Lx, double *zz, double *res)
{

  if(zz[0] != DO_NOT_EVALUATE){
    res[0] = gsl_spline_eval(Lx->mom1_spline_ptr, zz[0], Lx->mom1_accel_ptr);
  }

  if(zz[1] != DO_NOT_EVALUATE){
    res[1] = gsl_spline_eval(Lx->mom2_spline_ptr, zz[1], Lx->mom2_accel_ptr);
  }

  if(zz[2] != DO_NOT_EVALUATE){
    res[2] = gsl_spline_eval(Lx->b1_LW_spline_ptr, zz[2], Lx->b1_LW_accel_ptr);
  }
  if(zz[3] != DO_NOT_EVALUATE){
    res[3] = gsl_spline_eval(Lx->b2_LW_spline_ptr, zz[3], Lx->b2_LW_accel_ptr);
  }
  return _SUCCESS_;
}



/**
 * Compute the multiplicity function needed to compute the halo mass function 
 * Three models are implemented: Press-Schechter, Sheth-Tormen and Tinker
 * see Pillepich et al arxiv: 0811.4176 for the expressions. 
 * 
 * @param sigma             Input: variance of matter fluctuations
 * @param mode_mf           Input: switch for setting the model of mass function, can be set to PSC, ST, TR
 * @return the multiplicity function  
 */

double mult_func(double sigma, long mode_mf)
{
	double f 	   = 0;
	double delta_c = 1.686;
	double nu      = delta_c/sigma;


	if(mode_mf == PSC)
		f = sqrt(2./M_PI)* nu * exp(-pow(nu,2.)/2.);
	else if(mode_mf == ST){
		double A = 0.3222;
		double a = 0.707;  ///  In Barkana & Loeb Rev a = 0.75
		double p = 0.3;
		
		f = A * sqrt(2.*a/M_PI) * nu* exp(-(a*pow(nu,2.)/2.)) * (1.+pow(pow(nu,-2.)/a,p));
	}
	else if(mode_mf == TR){
		//// Mass function of Tinker et al 2008 arxiv: 0803.2706 for Delta =200
		double A = 0.186;   
		double a = 1.47;
		double b = 2.57;
		double c = 1.19;	
		f        = A*(pow(sigma/b,-a) + 1.) *exp(-c/pow(sigma,2.));


		//// Mass function of Tinker et al 2010 arxiv:1001.3162 for Delta =200
		// double alpha = 0.368;
		// double beta  = 0.589;
		// double gamma = 0.864;
		// double phi   = -0.729;
		// double eta   = -0.243;
		// f 		       = alpha * (1.+pow(beta*nu,-2.* phi))*pow(nu,2.*eta)*exp(-gamma*pow(nu,2.)/2.);
	}	

	return f;
}	

/**
 * Compute the halo mass function for Press-Schechter, Sheth-Tormen and Tinker models
 * see Pillepich et al arxiv: 0811.4176 for the expressions. 
 * 
 * @param Cx            Input: Cosmology structure
 * @param M             Input: Halo mass function
 * @param z             Input: redshift
 * @param mode_mf       Input: switch for setting the model of mass function, can be set to PSC, ST, TR
 * @return the halo mass function 
 */
double mass_func(struct Cosmology *Cx, double M, double z, long mode_mf ) ///in unit of halos per Mpc^3 per solar mass, compared at z=0 with Murray etal https://arxiv.org/abs/1306.5140
{
	double f 	   = 0;
	double R       = R_scale(Cx,M);   
	double omega_m = Cx->cosmo_pars[3L] + Cx->cosmo_pars[4L];

	double rho_m   = omega_m* rhoc(Cx, 0.);
	double sigma   = sqrt(sig_sq(Cx, z, R));

	/////If calculating derivative of sigma numerically
	double Mmin 	 	    = (1-0.05)*M;
	double Mmax 	 	    = (1+0.05)*M;
	double Rmin 	 	    = R_scale(Cx, Mmin);
	double Rmax 	 	    = R_scale(Cx, Mmax);
	double sigma_min 	  = sqrtl(sig_sq(Cx, z, Rmin));
	double sigma_max 	  = sqrtl(sig_sq(Cx, z, Rmax));
 	double der_lnsig_dM = (log(sigma_max) - log(sigma_min))/(2.*0.05*M);

	f = - mult_func(sigma,mode_mf) * rho_m/M * der_lnsig_dM;  

	return f;

}

/**
 * Read in the measured mass function of Hidden-valey sims and build an interpolator for HMF(M) for a fixed redshift.
 * @param Cx            Input: Cosmology structure
 * @param M             Input: halo mass
 * @param z             Input: redshift
 * @param mode_mf       Input: switch for setting the model of mass function, can be set to PSC, ST, TR
 * @return the interpolated measured halo mass function       
 */

double mass_func_sims(struct Cosmology *Cx, double M, double z, long mode_mf) ///M in unit of M_sun and HMF in unit of #-of-halos/Mpc^3/M_sun
{
    ////NOTE: not to be used at z other than 2. This function is for testing purposes only. We test the theoretical predictions at z=2 if ST MF or measured MF are used. 

    char HMF_filename[FILENAME_MAX];
    sprintf(HMF_filename,"/Volumes/Data/Documents/Git/LIM_PS_HM/Output/HMF/hmf-z2.00.dat");
    int nlines = 0;
    nlines     = count_lines_in_file(HMF_filename);

    double *log10M_input = make_1Darray(nlines);
    double *hmf_input    = make_1Darray(nlines); 
    double *log10hmf_in  = make_1Darray(nlines);

    static gsl_interp_accel   *hmf_accel_ptr;
    static gsl_spline         *hmf_spline_ptr;

    static int first = 1;
    if(first ==1 )
    {
      FILE *HMF_file;
      HMF_file = fopen(HMF_filename,"r");

      char line[MAXL];
      int i =0;
      while(fgets(line, sizeof line, HMF_file) != NULL )
      { 
        if(*line == '#')  continue; 
        sscanf(line,"%lg %lg\n",&log10M_input[i],&hmf_input[i]); 
        log10hmf_in[i] = log10(1./log(10.)*hmf_input[i]);
        // printf("%12.6e %12.6e \n", log10M_input[i],log10hmf_in[i]);

        i += 1;   
      }
      fclose(HMF_file);

      hmf_accel_ptr  = gsl_interp_accel_alloc();
      hmf_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

      gsl_spline_init(hmf_spline_ptr,log10M_input,log10hmf_in,nlines);

      first = 0;
    }

    double log10M   = log10(M*gb.h);  //convert to Msun/h
    double log10HMF = gsl_spline_eval(hmf_spline_ptr,log10M,hmf_accel_ptr); //in unit of (h/Mpc)^3 
    double HMF      = pow(gb.h,3.)*pow(10.,log10HMF)/M; //convert to unit of (1/Mpc)^3

    free(log10M_input);
    free(hmf_input);
    free(log10hmf_in);

   return HMF;
}


/**
 * computes the halo biases for three mass functions, press-schecter, Sheth-Tormen, and Tinker mass functions
 * 
 * @param Cx            Input: Cosmology structure
 * @param M             Input: halo mass
 * @param z             Input: redshift
 * @param mode_mf       Input: switch for setting the model of mass function, can be set to PSC, ST, TR
 * @param bias_arr      Output: the output array containning linear and quadratic local-in-matter halo biases, and quadratic and cubic tidal biases
 * @return void    
 */

void halo_bias(struct Cosmology *Cx, double M, double z, long mode_mf, double *bias_arr)
{
	double delta_c = 1.686;
	double R       = R_scale(Cx, M);   
	double sigma   = sqrt(sig_sq(Cx, z, R));
	double nu      = delta_c/sigma ;

  double epsilon1 = 0., epsilon2 = 0., E1 = 0., E2 =0.;
  double b1 =0., b2 =0., bg2=0., btd = 0.;
  double alpha = 0., p = 0.;
  double a =0., b =0., c=0.;

  ///Note that for PSC and ST mass functions, same form of the biases can be assumed, with different coefficents.
  ///See astro-ph/0006319
	if(mode_mf == PSC){
      epsilon1  = (pow(nu,2.)-1.)/delta_c;
      epsilon2  = pow(nu/delta_c,2.)*(pow(nu,2.)-3.);
      b1        = 1. + epsilon1; 
      b2        = 2. * (1.-17./21.) * epsilon1 + epsilon2;
  }
	else if(mode_mf == ST){ 
      ///Assuming spherical collapse
      alpha     = 0.707;
      p         = 0.3 ;
      epsilon1  = (alpha*pow(nu,2.)-1.)/delta_c;
      epsilon2  = alpha*pow(nu/delta_c,2.)*(alpha*pow(nu,2.)-3.);
      E1        = (2.*p/delta_c)/(1.+pow(alpha*pow(nu,2.),p));
      E2        = E1 * ((1.+2.*p)/delta_c + 2.*epsilon1);  

      b1        = 1. + epsilon1 + E1;
      b2        = 2. * (1.-17./21.) * (epsilon1 + E1) + epsilon2 + E2;
    
    ////Assuming ellipsoidal collapse, astro-ph/9907024
  		// double a = 0.707;
  		// double b = 0.5;
  		// double c = 0.6;
  		// b1 = 1.+1./(sqrt(a) * delta_c) * (sqrt(a)*(a*pow(nu,2.)) + sqrt(a) *b * pow(a*pow(nu,2.),1.-c)\
  		//   - pow(a*pow(nu,2.),c)/(pow(a*pow(nu,2.),c)+b*(1.-c)*(1.-c/2.))); 
	}
	else if(mode_mf == TR){
		//// Bias of Tinker et al 2008 arxiv: 0803.2706 for Delta =200
      a = 0.707;   
		  b = 0.35;
		  c = 0.80;
		  b1 = 1.+1./(sqrt(a)*delta_c)*(sqrt(a)*(a*pow(nu,2.))+sqrt(a)*b*pow(a*pow(nu,2.),1.-c)\
		     - pow(a*pow(nu,2.),c)/(pow(a*pow(nu,2.),c)+b*(1.-c)*(1.-c/2.))); 
  	
		//// Bias of Tinker et al 2010 arxiv:1001.3162 for Delta =200
		// double Delta = 200.;  
		// double y = log10(Delta);
		// double A = 1. + 0.24 *y * exp(-pow(4./y,4.));
		// double B = 0.183 ;
		// double C = 0.019 + 0.107 *y + 0.19 *exp(-pow(4./y,4.));
		// double aa = 0.44 * y - 0.88;
		// double bb = 1.5;
		// double cc = 2.4;
		// b1 = 1. - A*pow(nu,aa)/(pow(nu,aa) + pow(delta_c,aa)) + B *pow(nu,bb) + C* pow(nu,cc);
	}

  bg2 = -2./7. * (b1-1.);
  btd =  23./42. * (b1-1.);   

  bias_arr[0] = b1;
  bias_arr[1] = b2;
  bias_arr[2] = bg2;
  bias_arr[3] = btd;

  // printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",nu,M*gb.h,b1,b2,bg2,btd);

	return;
}


/**
 * Read in the file for the star formation rate byy Behroozi et al 2013
 * 
 * @param z_arr         Output: pointer to an array of redshifts read from the file
 * @param logM_arr      Output: pointer to an array of halo masses read from the file
 * @param log10SFR      Output: pointer to an array of SFR read from the file
 * @return void    
 */

void logSFR_Behroozi_read(double *z_arr, double *logM_arr, double *log10SFR)
{
  extern struct globals gb;
  double  result = 0.0;
  long  i, j, l, m;
  FILE  *ifp;
  int   verbose =1;
  double M_min, M_max;

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t numlines = count_lines_in_file(gb.SFR_filename); 
  size_t num_z = 137 ;
  size_t num_M = numlines/num_z;

  double   *zp1, *log10Mh, *log10SM;

  zp1     = make_1Darray(numlines*sizeof(double));
  log10Mh = make_1Darray(numlines*sizeof(double));
  log10SM = make_1Darray(numlines*sizeof(double));

  // Open the file
  ifp = fopen(gb.SFR_filename,"r");
  
  if(ifp == NULL)
  {
    printf("Failed to open the file");
    exit(1);
  }
 
  //// Write a function that counts the columns of a file, read ead each line and parse the line by a given delimiter count_columns(filename, token) (return integer)
  for(i=0L;i<numlines;i++){
      fscanf(ifp,"%lg %lg %lg %lg \n", &zp1[i],&log10Mh[i],&log10SFR[i],&log10SM[i]);
      // printf("%12.6e %12.5e %12.6e \n",zp1[i],log10Mh[i],log10SFR[i]);
  }


  gb.M_min =  pow(10.,log10Mh[0]); 
  gb.M_max =  pow(10.,log10Mh[numlines-1]);  
  gb.z_max = zp1[num_z-1] - 1;

  //construct an array of unique values of log10Mh
  j =0;
  for(i=0L;i<numlines;i++){
    if(log10Mh[i+1] != log10Mh[i]){
        logM_arr[j] = log10Mh[i];

        j += 1;  
    }
  }

  for(l=0;l<num_z;l++){
    z_arr[l] = zp1[l] - 1. ;  
  }

  // for(i=0;i<numlines;i++){
  //     if(logSFR_arr[i] == -1000)
  //           logSFR_arr[i] = 0.;
  //     else 
  //           logSFR_arr[i] = log10SFR[i];
  // }



  free(zp1);
  free(log10Mh);
  free(log10SM);

  fclose(ifp);
    
  return;
}      



/**
 * Allocate memory and initialize the 2d interpolator for the star formation rate of Behroozi et al 2013 as a function of halo mass and redshift
 * 
 * @return the error status
 */

int logSFR_alloc_init()
{
  extern struct globals gb;
  
  size_t numlines = count_lines_in_file(gb.SFR_filename); 
  size_t num_z    = 137 ;
  size_t num_M    = numlines/num_z;

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gb.logM_accel_ptr       = gsl_interp_accel_alloc();
  gb.z_accel_ptr          = gsl_interp_accel_alloc();
  gb.logSFR_spline2d_ptr  = gsl_spline2d_alloc(T, num_M, num_z);

  double *z_arr      = make_1Darray(num_z*sizeof(double));
  double *logM_arr   = make_1Darray(num_M*sizeof(double));
  double *log10SFR   = make_1Darray(numlines*sizeof(double));
  double *logSFR_arr = make_1Darray(numlines*sizeof(double));

  logSFR_Behroozi_read(z_arr,logM_arr,log10SFR);

  int i = 0, l, m;
  for(l=0;l<num_M;l++)
     for(m=0;m<num_z;m++){
        gsl_spline2d_set(gb.logSFR_spline2d_ptr, logSFR_arr, l, m, log10SFR[i]);
        i += 1;
  }
  gsl_spline2d_init(gb.logSFR_spline2d_ptr, logM_arr, z_arr, logSFR_arr, num_M, num_z);
  free(logM_arr);
  free(z_arr);
  free(logSFR_arr);
  free(log10SFR);

  return _SUCCESS_;
}


/**
 * Free the memory allocated to the interpolators of star formation rate by Behroozi et al 2013
 * 
 * @return the error status
 */

int SFR_Behroozi_free()
{
  gsl_interp_accel_free(gb.logM_accel_ptr);
  gsl_interp_accel_free(gb.z_accel_ptr);
  gsl_spline2d_free(gb.logSFR_spline2d_ptr);

  return _SUCCESS_;
}



/**
 * Evaluate the SFR interpolator object for a given value of mass and redshift
 * 
 * @param logM          Input: log10 of halo mass
 * @param z             Input: redshift
 * @return log10SFR    
 */

double logSFR_Behroozi(double logM, double z)
{
  double m1, m2, result = 0, eval;
  if(pow(10.,logM)> gb.M_max || pow(10.,logM)<gb.M_min){
      result = -1000.;
      // printf("ERROR, choice of M is outside interpolation limits \n");
  }
  else if(pow(10.,logM)< gb.M_max && pow(10.,logM)>gb.M_min){
      if(z <= gb.z_max) 
         eval =  gsl_spline2d_eval(gb.logSFR_spline2d_ptr, logM, z, gb.logM_accel_ptr, gb.z_accel_ptr );
      else if(z>gb.z_max){
        m1 = gsl_spline2d_eval(gb.logSFR_spline2d_ptr, logM, gb.z_max, gb.logM_accel_ptr, gb.z_accel_ptr ) + 0.2943*(z-8.);
        m2 = 3.3847-0.2413*z;
        eval = GSL_MIN(m1,m2);
      }
      if(eval < 0. &&  fabs(eval) > 500.)
         result = -1000.;
      else 
      result = eval;  
  }

  return result;
}


/**
 * Compute the line specific luminosity in unit of solar luminosity 
 * For CO ladder, I am using the fits in Table 4 of ??? et al arXiv:1508.05102, while for CII we use Silva et al arXiv: 
 * 
 * @param M          Input: halo mass
 * @param z          Input: redshift
 * @param mode_lum   Inpute: which luminosity model, basically which line considered
 * @return line luminosity    
 */

double luminosity(double M, double z, long mode_lum)
{
  double f=0., L_CO = 0.;
  double delta_MF=0. , nu_CO10=0. , L_IR=0. ,Lprime_CO=0.;
  double a_CO  = 0., b_CO  = 0.;
  double a_CII = 0., b_CII = 0.;
  double logM  = log10(M);
  int J = 0;
  
  if(mode_lum == CII){

    a_CII = 0.8475;  ////M1
    b_CII = 7.2203;

    // a_LCII = 1.0000;   ////M2
    // b_LCII =  6.9647;

    // a_LCII = 0.8727;   ///M3
    // b_LCII = 6.7250;

    // a_LCII = 0.9231;   ///M4
    // b_LCII = 6.5234; 

    if(logSFR_Behroozi(logM, z) == -1000.)
      f = 0.;
    else 
      f = pow(10.,a_CII* logSFR_Behroozi(logM, z) + b_CII);

  }
  else{ 
            delta_MF = 1.;
            nu_CO10 =  115.27;  
            if(logSFR_Behroozi(logM, z) == -1000.)
                  L_IR = 0;
            else 
                  L_IR = pow(10.,logSFR_Behroozi(logM, z))/delta_MF * pow(10.,10.); //Kennicutt relation 1998
            
            if(mode_lum == CO10){
                  a_CO = 1.27;  /// a = 1.37 Charilli
                  b_CO =-1.0; /// b = -1.74
                  J = 1;
            }
            else if(mode_lum == CO21){
                  a_CO = 1.11;
                  b_CO = 0.6;
                  J = 2;
            }
            else if(mode_lum == CO32){
                  a_CO = 1.18;
                  b_CO = 0.1;
                  J = 3;
            }
            else if(mode_lum == CO43){
                  a_CO = 1.09;
                  b_CO = 1.2;
                  J =4;
            }
            else if(mode_lum == CO54){
                  a_CO = 1.05;
                  b_CO = 1.8;
                  J = 5;
            }
            else if(mode_lum == CO65){
                  a_CO = 1.04;
                  b_CO = 2.2;
                  J = 6;
            }
            Lprime_CO = pow(10.,(log10(L_IR)-b_CO)/a_CO); ///in unit of K km/s pc^2
            L_CO = 4.9 * pow(10.,-5.) * pow(J,3.)* Lprime_CO;   ///in unit of L_sun
            f = L_CO;
      }

  return f;
}


/**
 * Compute the first luminosityy-weighted mass moment. 
 * The function mass_moment1_integ() is the integrand and mass_moment1() compute the moment
 * 
 * @param Cx         Input: pointer to cosmology structure
 * @param z          Input: redshift
 * @param M_min      Input: minimum halo mass
 * @param mode_mf    Input: model of halo mass function to consider, PSC, ST, TR
 * @param mode_lum   Inpute: which luminosity model, basically which line considered
 * @return the first mass moment    
 */

int mass_moment1_integ(unsigned       nd,                   // Number of dimensions in the domain space -- number of dim we're integrating over
                      const double  *x,                               // The point at which the integrand is evaluated
                      void          *p,                               // Pointer to a structure that holds the parameters
                      unsigned      fdim,                             // Number of dimensions that the integrand return
                      double        *fvalue                           // Array of values of the integrand of dimension fdim
                      )
{
      double f      = 0;
      double result = 0;
      double M    = exp(x[0]);

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);   

      struct Cosmology *Cx  = pij.p1;
      double z              = pij.p4;
      long mode_mf          = pij.p13;
      long mode_lum         = pij.p14;

      result = M * mass_func(Cx, M, z, mode_mf) * luminosity(M, z, mode_lum);
      *fvalue = result;

      return 0;
}

double mass_moment1(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum)  /// in unit of M_sun/Mpc^3
{
      struct integrand_parameters2 par;

      unsigned fdim      = 1;                         // Dimensionality of the integrand function
      unsigned dim       = 1;                         // Dimensionality of the domain of integration
      double   xmin[1], xmax[1];                      // Integration limits: these are arrays of dimension dim
      unsigned maxEval = 0;                           // Maximum number of integrand evaluations (0 for none)
      double   AbsErr  = 0.0;                         // Required absolute error (0.0 for none)
      double   RelErr  = 1.0e-3;                      // Required relative error (1.0e-3)
      double   result  = 0.;                          // Final result
      double   error   = 0.;                          // Error estimate on the result
      error_norm norm = ERROR_INDIVIDUAL;

      xmin[0] = log(M_min);  ///In units of solar mass;
      xmax[0] = log(1.e16);   ///In units of solar mass

      // xmin[0] = log(pow(10.,9.280699860662846135)/gb.h);
      // xmax[0] = log(pow(10.,1.443569618680730571e+01)/gb.h);

      par.p1  = Cx;
      par.p4  = z;
      par.p13 = mode_mf;
      par.p14 = mode_lum;
      
      hcubature(fdim, mass_moment1_integ, &par, dim, xmin, xmax, maxEval, AbsErr, RelErr, norm, &result, &error);
      
      return result;    
}


/**
 * Compute the second luminosityy-weighted mass moment. 
 * The function mass_moment2_integ() is the integrand and mass_moment2() compute the moment
 * 
 * @param Cx         Input: pointer to cosmology structure
 * @param z          Input: redshift
 * @param M_min      Input: minimum halo mass
 * @param mode_mf    Input: model of halo mass function to consider, PSC, ST, TR
 * @param mode_lum   Inpute: which luminosity model, basically which line considered
 * @return the second mass moment    
 */
int mass_moment2_integ(unsigned           nd,               // Number of dimensions in the domain space -- number of dim we're integrating over
                      const double  *x,                               // The point at which the integrand is evaluated
                      void          *p,                               // Pointer to a structure that holds the parameters
                      unsigned      fdim,                             // Number of dimensions that the integrand return
                      double        *fvalue                           // Array of values of the integrand of dimension fdim
                      )
{
  double f      = 0;
  double result = 0;
  double M  = exp(x[0]);
 
  struct integrand_parameters2 pij;
  pij = *((struct integrand_parameters2 *)p);   

  struct Cosmology *Cx  = pij.p1;
  double z              = pij.p4;
  long mode_mf          = pij.p13;
  long mode_lum         = pij.p14;

  result = M * mass_func(Cx, M, z, mode_mf)* pow(luminosity(M, z, mode_lum),2.);

  *fvalue = result;

  return 0;
}


double mass_moment2(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum)  /// in unit of M_sun/Mpc^3
{
      struct integrand_parameters2 par;

      unsigned fdim    = 1;                           // Dimensionality of the integrand function
      unsigned dim     = 1;                           // Dimensionality of the domain of integration
      double  xmin[1], xmax[1];                       // Integration limits: these are arrays of dimension dim
      unsigned maxEval = 0;                           // Maximum number of integrand evaluations (0 for none)
      double   AbsErr  =0.0;                          // Required absolute error (0.0 for none)
      double   RelErr  =1.0e-3;                       // Required relative error (1.0e-3)
      double   result  =0.;                           // Final result
      double   error   =0.;                           // Error estimate on the result
      error_norm norm = ERROR_INDIVIDUAL;

      xmin[0] = log(M_min);  ///In units of solar mass;
      xmax[0] = log(1.e16);   ///In units of solar mass

      // xmin[0] = log(pow(10.,9.280699860662846135)/gb.h);
      // xmax[0] = log(pow(10.,1.443569618680730571e+01)/gb.h);


      par.p1  = Cx;
      par.p4  = z;
      par.p13 = mode_mf;
      par.p14 = mode_lum;
      
      hcubature(fdim, mass_moment2_integ, &par, dim, xmin, xmax, maxEval, AbsErr, RelErr, norm, &result, &error);
      
      return result;    
}

///// Moments of mass function needed to calculate  the mean intensity and shot noise. 
// static int mass_moment1_integ(const int *ndim,
//            const cubareal x[],
//            const int *ncomp,
//            cubareal ff[],
//            void *p)  
// {
//     double f      = 0;
//     double result = 0;

//     struct integrand_parameters2 pij;
//    	pij = *((struct integrand_parameters2 *)p);	
    
//     struct Cosmology *Cx = pij.p1;
//     double z             = pij.p4;
//     double logMmin       = pij.p5;
//     double logMmax       = pij.p6;
//     long   mode_mf       = pij.p13;
//     long   mode_lum      = pij.p14;

//     double logM          = (x[0] * (logMmax - logMmin) + logMmin); 
//     double cos           = (2. * x[1] - 1.); 
//     double M             = exp(logM);

//     double MF            = mass_func(Cx, M, z, mode_mf);
//     double lum           = luminosity(M, z, mode_lum);

//     ff[0] = 0.5 * 2. * (logMmax - logMmin) * M * MF * lum;;

// 	 return 0;

// }

// double mass_moment1(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum)  /// in unit of M_sun/Mpc^3
// {
//    struct integrand_parameters2 par;

//   double AbsErr = 0.0;        // Required absolute error (0.0 for none)
//   double RelErr = 1.e-2; 

//   par.p5  = log(M_min);  ///In units of solar mass;
//   par.p6  = log(1.e16);   ///In units of solar mass

//   par.p1  = Cx;
//   par.p4  = z;
//   par.p13 = mode_mf;
//   par.p14 = mode_lum;
  	
//   int ndim = 2,  ncomp = 1, nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 2.e8, key = 13;
//   int nregions, neval, fail;
//   double prob, result, error;
      
//   Cuhre(ndim,ncomp, mass_moment1_integ, &par, nvec,
//             RelErr, AbsErr, verbose | last, mineval, maxeval, key,
//             NULL, NULL, &nregions, &neval, &fail, &result, &error, &prob);
  

// 	return result;	
// }

// static int mass_moment2_integ(const int *ndim,
//            const cubareal x[],
//            const int *ncomp,
//            cubareal ff[],
//            void *p)  
// {

//   struct integrand_parameters2 pij;
//   pij = *((struct integrand_parameters2 *)p); 
  
//   struct Cosmology *Cx = pij.p1;
//   double z             = pij.p4;
//   double logMmin       = pij.p5;
//   double logMmax       = pij.p6;
//   long   mode_mf       = pij.p13;
//   long   mode_lum      = pij.p14;

//   double logM          = (x[0] * (logMmax - logMmin) + logMmin); 
//   double cos           = (2. * x[1] - 1.); 
//   double M             = exp(logM);

//   double MF            = mass_func(Cx, M, z, mode_mf);
//   double lum           = luminosity(M, z, mode_lum);

//   ff[0] = 0.5 * 2. * (logMmax - logMmin) * M * MF * pow(lum,2.);

//   return 0;

// }


// double mass_moment2(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum)  /// in unit of M_sun/Mpc^3
// {
// 	struct integrand_parameters2 par;

//   double AbsErr = 0.0;        // Required absolute error (0.0 for none)
//   double RelErr = 1.e-2; 

//   par.p5 = log(M_min);  ///In units of solar mass;
//   par.p6 = log(1.e16);   ///In units of solar mass

//   par.p1  = Cx;
//   par.p4  = z;
//   par.p13 = mode_mf;
//   par.p14 = mode_lum;
    
//   int    ndim = 2,  ncomp = 1, nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 2.e8, key = 13;
//   int    nregions, neval, fail;
//   double prob, result, error;
      
//   Cuhre(ndim,ncomp, mass_moment2_integ, &par, nvec,
//             RelErr, AbsErr, verbose | last, mineval, maxeval, key,
//             NULL, NULL, &nregions, &neval, &fail, &result, &error, &prob);
  

// 	return result;	
// }


// static int bias_lum_weighted_integ(const int *ndim,
//            const cubareal x[],
//            const int *ncomp,
//            cubareal ff[],
//            void *p)  
// {
   
//     struct integrand_parameters2 pij;
//     pij = *((struct integrand_parameters2 *)p);     

//     struct Cosmology *Cx = pij.p1;
//     double z             = pij.p4;
//     double logMmin       = pij.p6;
//     double logMmax       = pij.p7;
//     long   mode_mf       = pij.p13;
//     long   mode_lum      = pij.p14;

//     double logM          = (x[0] * (logMmax - logMmin) + logMmin); 
//     double cos           = (2. * x[1] - 1.); 
//     double M             = exp(logM);

//     double MF            = mass_func(Cx, M, z, mode_mf);
//     double lum           = luminosity(M, z, mode_lum);

//     double *bias_arr = make_1Darray(4);
//     halo_bias(Cx, M, z, mode_mf, bias_arr);
//     double b1 = bias_arr[0]; 
//     double b2 = bias_arr[1]; 
    
//     ff[0] = 0.5 * 2. * (logMmax - logMmin) * M * MF * b1 * lum;
//     ff[1] = 0.5 * 2. * (logMmax - logMmin) * M * MF * b2 * lum;
    
//     free(bias_arr);

//     // if(ff[0]  != 0. && ff[1] != 0.)
//     //   printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n", z, M,MF,lum,b1,b2,ff[0],ff[1]);
    
//     return 0;

// }

// void bias_lum_weighted(struct Cosmology *Cx, double z, double M_min, 
//                           long mode_mf, long mode_lum, double *result)  
// {
//   struct integrand_parameters2 par;

//   double AbsErr = 0.0;        // Required absolute error (0.0 for none)
//   double RelErr = 1.e-2; 

//   par.p6 = log(M_min);  ///In units of solar mass;
//   par.p7 = log(1.e16);   ///In units of solar mass

//   par.p1  = Cx;
//   par.p4  = z;
//   par.p13 = mode_mf;
//   par.p14 = mode_lum;
    
//   int    ndim = 2,  ncomp = 2, nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 2.e8, key = 13;
//   int    nregions, neval;
//   int    fail[ncomp];
//   double error[ncomp];
//   double prob[ncomp];

//   Cuhre(ndim,ncomp, bias_lum_weighted_integ, &par, nvec,
//             RelErr, AbsErr, verbose | last, mineval, maxeval, key,
//             NULL, NULL, &nregions, &neval, fail, result, error, prob); 

//   // for(int i=0;i<ncomp;i++)
//   //   printf("%d %12.6e %12.6e \n", i, result[i],error[i]);

//   return;    
// }


/**
 * Compute the luminosityy-weighted linear and quadratic line biases. The normalization of first mass moment is not included yet. 
 * The function bias_lum_weighted_integ() is the integrand and bias_lum_weighted() computes the bias
 * 
 * @param Cx         Input: pointer to cosmology structure
 * @param z          Input: redshift
 * @param M_min      Input: minimum halo mass
 * @param mode_mf    Input: model of halo mass function to consider, PSC, ST, TR
 * @param mode_lum   Inpute: which luminosity model, basically which line considered
 * @return un-normalized line bias   
 */

int bias_lum_weighted_integ(unsigned  nd,               // Number of dimensions in the domain space -- number of dim we're integrating over
                        const double  *x,                   // The point at which the integrand is evaluated
                        void          *p,               // Pointer to a structure that holds the parameters
                        unsigned      fdim,         // Number of dimensions that the integrand return
                        double        *fvalue           // Array of values of the integrand of dimension fdim
                        )
{
      double f      = 0;
      double result = 0;
      double M      = exp(x[0]);
       
      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);   

      struct Cosmology *Cx = pij.p1;
      double z             = pij.p4;
      long   mode_mf       = pij.p13;
      long   mode_lum      = pij.p14;

      double MF    = mass_func(Cx, M, z, mode_mf);
      double lum   = luminosity(M, z, mode_lum);

      double *bias_arr = make_1Darray(4);
      halo_bias(Cx, M, z, mode_mf, bias_arr);
      double b1 = bias_arr[0]; 
      double b2 = bias_arr[1]; 
    
      fvalue[0] = M * MF * b1 * lum;
      fvalue[1] = M * MF * b2 * lum;

      return 0;
}


void bias_lum_weighted(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum, double *result)  
{

      struct integrand_parameters2 par;

      unsigned fdim    = 2;                           // Dimensionality of the integrand function
      unsigned dim     = 1;                             // Dimensionality of the domain of integration
      double   xmin[1], xmax[1];                        // Integration limits: these are arrays of dimension dim
      unsigned maxEval = 0;                           // Maximum number of integrand evaluations (0 for none)
      double   AbsErr  = 0.0;                         // Required absolute error (0.0 for none)
      double   RelErr  = 1.0e-3;                      // Required relative error (1.0e-3)
      double   error[2];                          // Error estimate on the result
      error_norm norm = ERROR_INDIVIDUAL;

      xmin[0] = log(M_min);  ///In units of solar mass;
      xmax[0] = log(1.e16);   ///In units of solar mass

      // xmin[0] = log(pow(10.,9.280699860662846135)/gb.h);
      // xmax[0] = log(pow(10.,1.443569618680730571e+01)/gb.h);

      par.p1  = Cx;
      par.p4  = z;
      par.p13 = mode_mf;
      par.p14 = mode_lum;
      
      hcubature(fdim, bias_lum_weighted_integ, &par, dim, xmin, xmax, maxEval, AbsErr, RelErr, norm, result, error);
      
      return;    
}


/**
 * Model from Keating et al 2016 to account for the observed variation in halo activity, i.e. scatter in the L(M) relation
 * p_sig_shot replaces the f_duty in the shot-noise used in some LIM paper (ex. Lidz et al 2011).
 * p_sig_shot_integrand() is the integrand, and p_sig_shot() computes the scatter factor for the shot noise.
 * 
 * @param scatter    Input: variance of the log-scatter
 * @return the scatter coeff of the shot noise  
 */

double p_sig_shot_integrand(double x, void *par)
  {
    double f=0;
    double result = 0;
    double scatter;

    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);

    scatter = pij.p4 ;

    result = pow(10.,2.*x)/sqrtl(2.*M_PI*pow(scatter,2.))* exp(-pow(x/scatter,2.)/2.);

    return result;

  } 

double p_sig_shot(double scatter)
{
  extern struct globals gb;

  double result=0., error=0.;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

  struct integrand_parameters2 par; 

  double xmin = -10.;
  double xmax = 10.;

  gsl_function F;
  F.function = &p_sig_shot_integrand;
  F.params = &par;

  par.p4  = scatter;
 
  gsl_integration_qags(&F,xmin,xmax,0.0,1.0e-2,1000000,w,&result,&error);
  gsl_integration_workspace_free (w);
     
  return result;

}


/**
 * Model from Keating et al 2016 to account for the observed variation in halo activity, i.e. scatter in the L(M) relation
 * p_sig_Tbar replace the f_duty in the average brightness temprature used in some LIM paper (ex. Lidz et al 2011).
 * p_sig_Tbar_integrand() is the integrand, and p_sig_Tbar() computes the scatter factor for the mean brightness temprature.
 * 
 * @param scatter    Input: variance of the log-scatter
 * @return the scatter coeff of Tbar  
 */

double p_sig_Tbar_integrand(double x, void *par)
  {
    double f=0;
    double result = 0;
    double scatter;

    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);

    scatter = pij.p4 ;

    result = pow(10.,x)/sqrtl(2.*M_PI*pow(scatter,2.))* exp(-pow(x/scatter,2.)/2.);

    return result;

  } 

double p_sig_Tbar(double scatter)
{
  extern struct globals gb;

  double result=0., error=0.;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

  struct integrand_parameters2 par; 

  double xmin = -10.;
  double xmax = 10.;

  gsl_function F;
  F.function = &p_sig_Tbar_integrand;
  F.params = &par;

  par.p4  = scatter;
 
  gsl_integration_qags(&F,xmin,xmax,0.0,1.0e-2,1000000,w,&result,&error);
  gsl_integration_workspace_free (w);
     
  return result;

}

/**
 * Compute the linear and quadratic line biases, accounting ffor the normalization w.r.t. the first mass moment 
 *  
 * @param Lx        Input: Pointer to line structure
 * @param z         Input: Redshift
 * @param result    Input: a pointer to an array containing the results of b1_line and b2_line
 * @return void
 */

void line_bias(struct Line *Lx, double z, double *result)
{

  double zz[4] = {z,DO_NOT_EVALUATE,z,z};
  double res[4];

  Line_evaluate(Lx, zz, res);

  double mass_mom1 = res[0];
  double b1_lumW   = res[2]; 
  double b2_lumW   = res[3]; 

  result[0] = b1_lumW/mass_mom1;
  result[1] = b2_lumW/mass_mom1;
  

  return;
}


/**
 * Compute the line mean intensity in unit of erg Mpc^-2 Sr^-1
 * 
 * @param Cx        Input: Pointer to cosmology structure
 * @param line_id   Inpute: id of line of interest, an integer value
 * @param z         Input: Redshift
 * @return the line mean intensity
 */

double mean_intens(struct Cosmology *Cx, size_t line_id, double z)  
{
    ///Note: nu_J is the rest-frame emission frequency related to the observed frequency as nu_obs = nu_J/(1+z_J)
    /// For a CO transition from J-> J-1, the rest-frame frequency is nu_J = J nu_CO where nu_Co = 115 GHz. 

  double nu_line = Cx->Lines[line_id]->line_freq;
  double a = 1./(1.+z);
  double L_sun = 3.846 * 1.e33;  ///in unit of erg/s

  double zz[4] = {z,DO_NOT_EVALUATE,DO_NOT_EVALUATE,DO_NOT_EVALUATE};
  double res[4];
  Line_evaluate(Cx->Lines[line_id], zz, res);
  double mass_mom1 = res[0];

  double y_tilde = gb.c/nu_line * pow(1.+z,2.)/Hubble(Cx,z);  
  double fac = 1./(4.*M_PI) * y_tilde * pow(1.+z,-2.)*L_sun ;

  double f = fac * mass_mom1;

  return f;
}

  
/**
 * Compute the mean brightness temprature of CO in unit of microK, compared with Pullen et al and Lidz et al 2011
 *  
 * @param Cx        Input: Pointer to cosmology structure
 * @param line_id   Inpute: id of line of interest, an integer value
 * @param z         Input: Redshift
 * @return the line mean temprature assuming Rayleigh-Jeans limit
 */    

double Tbar_line(struct Cosmology *Cx, size_t line_id, double z) 
{   
  double f       = 0.;
  double k_B     = 1.38064852*1.e-16;  ///Boltzmann constant in unit of erg K^-1
  double fac     = 3.086*1.e19;   ////conversion factor from Mpc to km

  double nu_line = Cx->Lines[line_id]->line_freq;

  double sig_CO  = 0.37;
  double Tbar_scatter = p_sig_Tbar(sig_CO);

  f =  Tbar_scatter * pow(10.,6.)* pow(gb.c/nu_line,2.) *mean_intens(Cx,line_id,z) * pow(1.+z,2.)/(2.*k_B) * pow(fac,-2.);  ///factor of 10^6 is the conversion factor from K to microK
  
  return f;
}

