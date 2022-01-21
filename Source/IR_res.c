
/** @file IR_res.c Documented IR_res module 
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 * This module is computes the leading and next-to-leading IR-resummed matter power spectrum
 * The wiggle-nowiggle seperation is performed in wnw_split.c module. 
 *
 * In summary, the following functions can be called from other modules:
 * -# pm_IR_LO()
 * -# pm_IR_NLO()
 * -# IR_Sigma2()
 * -# pm_nowiggle()
 * -# pm_nowiggle_gfilter()
 * -# pm_nowiggle_bspline()
 * -# pm_nowiggle_dst()
 */

#include "header.h"
struct globals gb;


/**
 * Compute the leading-order IR-resummed matter power spectrum, ala Ivanovic et al
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber in unit of 1/Mpc. 
 * @param z            Input: redshift
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 * @return value of leading IR-ressumed power spectrum           
 */
double pm_IR_LO(struct Cosmology *Cx, double k, double z, long SPLIT)
{
    double kf0 = 1.e-4;
    static double sig2_LO = - 1.;
    if(sig2_LO == - 1.){
          sig2_LO = IR_Sigma2(Cx,z, kf0, SPLIT);
    }
    double p_nowiggle = pm_nowiggle(Cx, k, z, kf0, 0, SPLIT);
    double p_wiggle   = PS(Cx, k, z) - p_nowiggle;
    double sup        = exp(-k * k * sig2_LO);
    double f          = p_nowiggle + sup * p_wiggle;

    return f;

}

/**
 * Compute the next-to-leading-order IR-resummed matter power spectrum, ala Ivanovic et al
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber in unit of 1/Mpc. 
 * @param z            Input: redshift
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 * @return value of NL IR-ressumed power spectrum           
 */
double pm_IR_NLO(struct Cosmology *Cx, double k,  double z, long SPLIT)
{
    static double sig2_NLO = - 1.;
    double k0 = 1.e-4;

    if(sig2_NLO == - 1.){
      sig2_NLO = IR_Sigma2(Cx,z, k0, SPLIT);
    }

    double p_nowiggle = pm_nowiggle(Cx, k, z, k0, 0, SPLIT);
    double p_wiggle   = PS(Cx, k, z) - p_nowiggle;
    double sup        = exp(-k * k * sig2_NLO);

    double *pm_loops  =  make_1Darray(2);
    Compute_G_loops(Cx, k, z, WIR, MATTER, SPLIT,pm_loops);

    double p22_IR = pm_loops[0];
    double p13_IR = pm_loops[1];
    free(pm_loops);

    double pm_LO = p_nowiggle + sup * p_wiggle; 
    double f     = p_nowiggle + sup * p_wiggle * (1. + k * k * sig2_NLO) + p22_IR + p13_IR ;

    return f;
}



/**
 * Integrand to compute the suppression factor IR_sigma2 
 * 
 * @param x            Input: integration variable, k-values
 * @param par          Input: integration parameters
 * @return integrand to be used in IR_sigma2() function         
 */
double IR_Sigma2_integrand(double x, void *par)
{
    double result = 0;
    
    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);

    struct Cosmology *Cx = pij.p1;
    double bao_scale     = pij.p4;
    double z             = pij.p5;
    double k0           = pij.p6; 
    long   SPLIT         = pij.p13;

    double k_osc = 1./bao_scale;  /// BAO_scale = 110. Mpc/h.

    result = 1./(6.*M_PI*M_PI)*pm_nowiggle(Cx, x, z, k0, 0, SPLIT)* (1. - gsl_sf_bessel_j0(x/k_osc) + 2. * gsl_sf_bessel_j2(x/k_osc));;

    return result;

} 


/**
 * Compute the suppression factor IR_sigma2 
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param z            Input: redshift
 * @param k0          Input: first element of the k-array, used in normalization of EH no-wiggle spectrum
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 * @return value of IR resummation suppression factor          
 */
double IR_Sigma2(struct Cosmology *Cx, double z, double k0, long SPLIT)
{
    double result=0., error=0.;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    struct integrand_parameters2 par; 

    double kmin = 1.e-6;
    double kmax = 0.2;

    gsl_function F;
    F.function = &IR_Sigma2_integrand;
    F.params = &par;

    par.p1  = Cx;
    par.p4  = 110.;
    par.p5  = z;
    par.p6  = k0;
    par.p13 = SPLIT;

    gsl_integration_qags(&F,kmin,kmax,0.0,1.0e-3,1000000,w,&result,&error);
    gsl_integration_workspace_free(w);

    return result;
}



/**
 * Compute the no-wiggle componenet of the matter power spectrum
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber in unit of h/Mpc. 
 * @param z            Input: redshift
 * @param kf0          Input: first element of the k-array, used in normalization of EH no-wiggle spectrum
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 * @return value of no-wiggle power spectrum         
 */
double pm_nowiggle(struct Cosmology *Cx, double k, double z, double kf0, int cleanup, long SPLIT)
{
  double pm_nw = 0.;

  if(SPLIT == DST)
    pm_nw = pm_nowiggle_dst(Cx, k, z, cleanup);
  else if(SPLIT == GFILTER)
    pm_nw = pm_nowiggle_gfilter(Cx, k, z, cleanup);
  else if(SPLIT == BSPLINE){
    pm_nw = pm_nowiggle_bspline(Cx, k, z, cleanup);
  }

  return pm_nw;
}


/**
 * Compute the no-wiggle componenet of the matter power spectrum, reading in and interpolating the output of apython code which computed the broadband by fitting families of Bsplines (see Vlah et al 2015)
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber in unit of h/Mpc. 
 * @param z            Input: redshift
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 * @return value of no-wiggle power spectrum         
 */
double pm_nowiggle_bspline(struct Cosmology *Cx, double k, double z, int cleanup)
{

  static gsl_interp_accel   *pknw_accel_ptr;
  static gsl_spline         *pknw_spline_ptr;
  FILE *pksplit_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;
  if(first == 1){
    char pksplit_filename[FILENAME_MAX];
    sprintf(pksplit_filename,"Input/wnw_split/pk_Bspline.txt");

    int nlines = 0;
    nlines     = count_lines_in_file(pksplit_filename);

    double *k_in, *pk_nw, *log_k, *log_pknw;
    k_in     = make_1Darray(nlines);
    log_k    = make_1Darray(nlines);
    pk_nw    = make_1Darray(nlines);
    log_pknw = make_1Darray(nlines);
    
    pksplit_file = fopen(pksplit_filename,"r");

    if(pksplit_file==NULL){
      printf("Failed to open the file with k-values");
      exit(1);
    }

    char line[MAXL];
    int err;
    while(fgets(line, sizeof line, pksplit_file) != NULL )
    { 
        if(*line == '#')  continue; 
          for(i=0;i<nlines;i++){
            err         = fscanf(pksplit_file,"%lg %lg \n",&k_in[i],&pk_nw[i]); 
            log_k[i]    = log(k_in[i]);
            log_pknw[i] = log(pk_nw[i]/pow(2.*M_PI,3.)); 
          }
    }
    fclose(pksplit_file);

    pknw_accel_ptr  = gsl_interp_accel_alloc();
    pknw_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pknw_spline_ptr,log_k,log_pknw,nlines);

    kmin = k_in[0];
    kmax = k_in[nlines-1];

    free(k_in);
    free(log_k);
    free(pk_nw);
    free(log_pknw);
    
    first = 0;
  }   

 
  double pknw, logpknw, logk;

  if(k<kmin || k >kmax) {
      pknw = 0.;
  }
  else{
      logk    = log(k);
      logpknw = gsl_spline_eval(pknw_spline_ptr,logk,pknw_accel_ptr);
      pknw    = exp(logpknw);
  }


  
  double growth2 = pow(growth_D(Cx, k, z),2.);  
  double pknw_z = growth2 * pknw;


  if (cleanup == 1){
        gsl_interp_accel_free(pknw_accel_ptr);
        gsl_spline_free(pknw_spline_ptr);
  }

  
  return pknw_z;     
}


/**
 * Compute the no-wiggle componenet of the matter power spectrum, using Gaussian filter (see Vlah et al 2015)
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber in unit of h/Mpc. 
 * @param z            Input: redshift
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 * @return value of no-wiggle power spectrum         
 */
double pm_nowiggle_gfilter(struct Cosmology *Cx, double k, double z, int cleanup)
{
  double pm_nw = 0.;
  static gsl_interp_accel   *pknw_accel_ptr;
  static gsl_spline         *pknw_spline_ptr;
  FILE *pksplit_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;
  if(first == 1){

    int nlines       = 600;
    double *k_in     = loginit_1Darray(nlines, 1.e-4, 20.);
    double *log_k    = make_1Darray(nlines);
    double *pk_nw    = make_1Darray(nlines);
    double *log_pknw = make_1Darray(nlines);
    
    for(i=0;i<nlines;i++){
        pk_nw[i]    = pk_Gfilter_nw(Cx, k_in[i], k_in[0]);
        log_k[i]    = log(k_in[i]);
        log_pknw[i] = log(pk_nw[i]);
        // printf("nw %d %12.6e %12.6e \n",i, k_in[i],pk_nw[i]);
    }


    pknw_accel_ptr  = gsl_interp_accel_alloc();
    pknw_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pknw_spline_ptr,log_k,log_pknw,nlines);

    kmin = k_in[0];
    kmax = k_in[nlines-1];

    free(k_in);
    free(log_k);
    free(pk_nw);
    free(log_pknw);
    
    first = 0;
  }   

  double pknw, logpknw, logk;

  if(k<kmin || k >kmax) {
      pknw = 0.;
  }
  else{
      logk    = log(k);
      logpknw = gsl_spline_eval(pknw_spline_ptr,logk,pknw_accel_ptr);
      pknw    = exp(logpknw);
  }
  

  double growth2 = pow(growth_D(Cx, k, z),2.);
  double pknw_z = growth2 * pknw;


  if (cleanup == 1){
        gsl_interp_accel_free(pknw_accel_ptr);
        gsl_spline_free(pknw_spline_ptr);
  }


  return pknw_z;   
}


/**
 * Compute the no-wiggle componenet of the matter power spectrum, reading in and interpolating the output of apython code which computed the broadband by discrete sin-transform, See Hamann et al 2010.
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber in unit of h/Mpc. 
 * @param z            Input: redshift
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 * @return value of no-wiggle power spectrum         
 */
double pm_nowiggle_dst(struct Cosmology *Cx, double k, double z, int cleanup)
{
  double pm_nw = 0.;
    static gsl_interp_accel *pknw_accel_ptr;
  static gsl_spline         *pknw_spline_ptr;
  FILE *pksplit_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;
  if(first == 1){
    char pksplit_filename[FILENAME_MAX];
    sprintf(pksplit_filename,"Input/wnw_split/pk_dst.txt");

    int nlines = 0;
    nlines     = count_lines_in_file(pksplit_filename);

    double *k_in, *pk_nw, *log_k, *log_pknw;
    k_in     = make_1Darray(nlines);
    log_k    = make_1Darray(nlines);
    pk_nw    = make_1Darray(nlines);
    log_pknw = make_1Darray(nlines);
    
    pksplit_file = fopen(pksplit_filename,"r");

    if(pksplit_file==NULL){
      printf("Failed to open the file with k-values");
      exit(1);
    }

    char line[MAXL];
    int err;
    while(fgets(line, sizeof line, pksplit_file) != NULL )
    { 
        if(*line == '#')  continue; 
          for(i=0;i<nlines;i++){
            err         = fscanf(pksplit_file,"%lg %lg \n",&k_in[i],&pk_nw[i]); 
            log_k[i]    = log(k_in[i]);
            log_pknw[i] = log(pk_nw[i]/pow(2.*M_PI,3.)); 
            //printf("%12.6e %12.6e \n",log_k[i], log_pknw[i]);
          }
    }
    fclose(pksplit_file);

    pknw_accel_ptr  = gsl_interp_accel_alloc();
    pknw_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pknw_spline_ptr,log_k,log_pknw,nlines);

    kmin = k_in[0];
    kmax = k_in[nlines-1];

    free(k_in);
    free(log_k);
    free(pk_nw);
    free(log_pknw);
    
    first = 0;
  }   

 
  double pknw, logpknw, logk;

  if(k<kmin || k >kmax) {
      pknw = 0.;
  }
  else{
      logk    = log(k);
      logpknw = gsl_spline_eval(pknw_spline_ptr,logk,pknw_accel_ptr);
      pknw    = exp(logpknw);
  }
  


  double growth2 = pow(growth_D(Cx, k, z),2.);  
  double pknw_z = growth2 * pknw;

  if (cleanup == 1){
        gsl_interp_accel_free(pknw_accel_ptr);
        gsl_spline_free(pknw_spline_ptr);
  }
    
  return pknw_z;    
}


