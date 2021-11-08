/** @file cosmology.c  Documented cosmology module 
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *   
 * The first routine of this module initalizes the Cosmology structure, which is the main building block of this entire code. 
 * This structure includes two sub-structures: the CLASS cosmology structure and line structure. Once the CLASS cosmology is initialized,
 * various useful functions can be directly called from CLASS, example to compute matter power spectrum and transfer function, 
 * angular and comoving radii, growth factor and growth rate, variance of matter fluctuations and its derivative. 
 * Lastly, the module also includes various window functions and their derivatives. 
 * 
 * In summary, the following functions can be called from other modules:
 * 
 * -# Cosmology_init()              allocates memory to and initializes cosmology structure
 * -# Cosmology_free()              frees the memory allocated to cosmology structure
 * -# CL_Cosmology_initilize()      initializes the class cosmology structure
 * -# CL_Cosmology_free()           frees the class cosmology structure
 * -# PS()                          computes matter power spectrum calling class function
 * -# Transfer()                    computes matter transfer function calling class function
 * -# growth_D()                    computes the scale-dep growth factor
 * -# growth_f()                    computes the scale-dep growth rate dlnD(k,a)/dlna
 * -# scale_indep_growth_D()        computes the scale-indep growth factor using directly CLASS functions
 * -# scale_indep_growth_f()        computes the scale-indep growth rate dlnD(k,a)/dlna  using directly CLASS functions
 * -# Hubble()                      computes hubbble parameter using directly CLASS functions
 * -# angular_distance()            computes angular diamtere distance  using directly CLASS functions
 * -# comoving_radial_distance()    computes radial distance  using directly CLASS functions
 * -# sig_sq()                      computes variance of smoothed matter fluctuations
 * -# der_sig_sq()                  computes derivative of the variance of smoothed matter fluctuations w.r.t. smoothing scale
 * -# sigma0_sq()                   computes variance of unsmoothed matter fluctuations
 * -# rhoc()                        computes the critical density of the universe 
 * -# R_scale()                     computes the size of a spherical halo corresponding to a given mass at z=0
 * -# R_scale_wrong()               computes the size of a spherical halo corresponding to a given mass at a given redshift
 * -# window_rth()                  computes top-hat filter in real space
 * -# window_g()                    computes Gaussian window
 * -# window_kth()                  computes top-hat filter in Fourier space
 * -# derR_window_rth()             computes derivative of top-hat filter in real space w.r.t. smoothing scale 
 * -# derR_logwindow_g()            computes derivative of top-hat filter in Fourier space w.r.t. smoothing scale
 */

#include "header.h"
struct globals gb;


/**
 * Allocate memory and initialize the cosmology structure, which includes the CLASS cosmology structure and line strucrure 
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param pk_kmax           Input: kmax for computation of matter power spectrum by CLASS
 * @param pk_zmax           Input: zmax for computation of matter power spectrum by CLASS
 * @param nlines            Input: number of lines whose properties we want to compute 
 * @param line_type         Inpute: name of the line to compute.
 *                                  It can be set to CII, CO10, CO21, CO32, CO43, CO54, CO65 
 * @param npoints_interp    Input: number of points in redshift for interpolation of line properties
 * @param M_min             Input: minimum halo mass for mass integrals
 * @param mode_mf           Inpute: theoretical model of halo mass function to use. 
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * @return an integer if succeeded
 */

int Cosmology_init(struct Cosmology *Cx, double pk_kmax, double pk_zmax, 
                  int nlines, int * line_types, size_t npoints_interp, double M_min, long mode_mf)
{

  long mode_nu = Cx->mode_nu;

  CL_Cosmology_initilize(Cx,pk_kmax,pk_zmax);

  if (nlines > 0){
    Cx->NLines = nlines;

    // Allocating array of line objects
    Cx->Lines = (struct Line **)malloc(sizeof(struct Line *) * nlines);

    for(int i=0; i<nlines; i++){
      Cx->Lines[i] = Line_alloc_init(Cx, line_types[i], npoints_interp, M_min, mode_mf);
    }
  }

  return _SUCCESS_;
}


/**
 * Free the memory allocated to cosmology structure
 * 
 * @param Cx    Input: pointer to Cosmology structure
 * @return the error status
 */

int Cosmology_free(struct Cosmology *Cx)
{
  for(int i=0; i<(Cx->NLines); i++){
    Line_free(Cx->Lines[i]);
  }

  free(Cx->Lines);

  CL_Cosmology_free(Cx);
  Pk_dlnPk(Cx, 0.1, 0., CLEANUP);

  return _SUCCESS_;
}


/**
 * Allocate memory and initialize the CLASS cosmology structure
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param pk_kmax           Input: kmax for computation of matter power spectrum by CLASS
 * @param pk_zmax           Input: zmax for computation of matter power spectrum by CLASS
 * @return the error status
 */

int CL_Cosmology_initilize(struct Cosmology *Cx, double pk_kmax, double pk_zmax)
{

  ////To interpolate P(k,z) at  various values of (k,z), enter 'z_max_pk', the maximum value of z at which
 ////such interpolations are needed. 


    struct  file_content       pfc;
    int     counter = 0;
    int     size_max = 40;

    pfc.size        = size_max;
    pfc.name        = (FileArg*)malloc(sizeof(FileArg) * pfc.size);
    pfc.value       = (FileArg*)malloc(sizeof(FileArg) * pfc.size);
    pfc.filename    = (char*)malloc(sizeof(char) * FILENAME_MAX);
    pfc.read        = (short*)malloc(sizeof(short) * pfc.size);


    //printf("Initializing Class\n");

    sprintf(pfc.filename,"doesnt_matter");

    sprintf(pfc.name[counter],"accurate_lensing");
    sprintf(pfc.value[counter],"1");                   
    pfc.read[counter] = _TRUE_;
    counter++;


    sprintf(pfc.name[counter],"h");
    sprintf(pfc.value[counter],"%18.12e", Cx->cosmo_pars[2]);  ///h
    pfc.read[counter] = _TRUE_;
    counter++;

    //// Photon density: T_cmb
    sprintf(pfc.name[counter],"T_cmb");
    sprintf(pfc.value[counter],"2.7255");                   
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"m_ncdm");
    sprintf(pfc.value[counter],"0.1");                   
    pfc.read[counter] = _TRUE_;
    counter++;

    // sprintf(pfc.name[counter],"N_ncdm");
    // sprintf(pfc.value[counter],"%d", 3);               
    // pfc.read[counter] = _TRUE_;
    // counter++;

    // sprintf(pfc.name[counter],"m_ncdm");
    // sprintf(pfc.value[counter],"%18.12e, %18.12e, %18.12e", 0.1/3.,0.1/3.,0.1/3.);               
    // pfc.read[counter] = _TRUE_;
    // counter++;


    /////Density of photons equivalant to T_CMB = 2.7255
    // sprintf(pfc.name[counter],"Omega_g");
    // sprintf(pfc.value[counter],"%12.6e ", gb.Omega_g);                   
    // pfc.read[counter] = _TRUE_;
    // counter++;


    //DE equation of state
    sprintf(pfc.name[counter],"w"); 
    sprintf(pfc.value[counter],"-1");      ///Omega_b        
    pfc.read[counter] = _TRUE_;
    counter++;

    // Density of Baryons:  'Omega_b' 
    sprintf(pfc.name[counter],"Omega_b"); 
    sprintf(pfc.value[counter],"%18.12e", Cx->cosmo_pars[3]);      ///Omega_b        
    pfc.read[counter] = _TRUE_;
    counter++;


    // Number of ultra-relativistic species / massless neutrino density:  'N_eff' 
    sprintf(pfc.name[counter],"N_eff");
    sprintf(pfc.value[counter],"3.046");               
    pfc.read[counter] = _TRUE_;
    counter++;


   // Density of cdm (cold dark matter): 'Omega_cdm' 
    sprintf(pfc.name[counter],"Omega_cdm");     
    sprintf(pfc.value[counter],"%18.12e",Cx->cosmo_pars[4]);              
    pfc.read[counter] = _TRUE_;
    counter++;

 
    // Primordial Helium fraction 'YHe', e.g. 0.25; if set to 'BBN' or 'bbn', will be inferred from Big Bang Nucleosynthesis (default: set to 'BBN')
    sprintf(pfc.name[counter],"YHe");
    sprintf(pfc.value[counter],"0.24");                                                        
    pfc.read[counter] = _TRUE_;
    counter++;


    //reionization optical depth
    sprintf(pfc.name[counter],"tau_reio");
    sprintf(pfc.value[counter],"0.0568");                                                        
    pfc.read[counter] = _TRUE_;
    counter++;

    ///pivot scale in unit of 1/Mpc
    sprintf(pfc.name[counter],"k_pivot");        
    sprintf(pfc.value[counter],"%12.6e ",gb.kp);                                                
    pfc.read[counter] = _TRUE_;
    counter++;

    // scalar adiabatic perturbations: curvature power spectrum value at pivot scale, tilt at the same scale
    sprintf(pfc.name[counter],"ln10^{10}A_s");  
    sprintf(pfc.value[counter],"%18.12e", Cx->cosmo_pars[0]);  
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"n_s");
    sprintf(pfc.value[counter],"%18.12e",Cx->cosmo_pars[1]);           
    pfc.read[counter] = _TRUE_;
    counter++;


    // Output
    sprintf(pfc.name[counter],"output");
    sprintf(pfc.value[counter],"mPk, mTk");
    pfc.read[counter] = _TRUE_;
    counter++;

    // maximum k in P(k), 'P_k_max_h/Mpc' in units of h/Mpc or 'P_k_max_1/Mpc' in units of 1/Mpc. 
    // If scalar Cls are also requested, a minimum value is automatically imposed (the same as in scalar Cls computation) (default: set to 0.1h/Mpc)

    //sprintf(pfc.name[counter],"P_k_max_h/Mpc");
    sprintf(pfc.name[counter],"P_k_max_1/Mpc");
    sprintf(pfc.value[counter],"%12.6e", pk_kmax);
    pfc.read[counter] = _TRUE_;
    counter++;


    sprintf(pfc.name[counter],"z_max_pk");
    sprintf(pfc.value[counter],"%12.6e", pk_zmax);
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"background_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"thermodynamics_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"perturbations_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"bessels_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"transfer_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"primordial_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"spectra_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;

    sprintf(pfc.name[counter],"nonlinear_verbose");
    sprintf(pfc.value[counter],"0");
    pfc.read[counter] = _TRUE_;
    counter++;



    pfc.size = counter;

  ///////////////////////////////////
    // Calling CLASS 2.5.0
    ///////////////////////////////////

  if (input_init(&pfc ,&Cx->ccs.pr,&Cx->ccs.ba,&Cx->ccs.th,&Cx->ccs.pt,&Cx->ccs.tr,&Cx->ccs.pm,&Cx->ccs.sp,&Cx->ccs.nl,&Cx->ccs.le,&Cx->ccs.op, Cx->ccs.errmsg) == _FAILURE_) {
        printf("\n\nError running input_init\n=>%s\n",Cx->ccs.errmsg); 
        return _FAILURE_;
  }

  if (background_init(&Cx->ccs.pr,&Cx->ccs.ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",Cx->ccs.ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&Cx->ccs.pr,&Cx->ccs.ba,&Cx->ccs.th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",Cx->ccs.th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&Cx->ccs.pr,&Cx->ccs.ba,&Cx->ccs.th,&Cx->ccs.pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",Cx->ccs.pt.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&Cx->ccs.pr,&Cx->ccs.ba,&Cx->ccs.th,&Cx->ccs.pt,&Cx->ccs.nl,&Cx->ccs.tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",Cx->ccs.tr.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&Cx->ccs.pr,&Cx->ccs.pt,&Cx->ccs.pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",Cx->ccs.pm.error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(&Cx->ccs.pr,&Cx->ccs.ba,&Cx->ccs.th,&Cx->ccs.pt,&Cx->ccs.pm,&Cx->ccs.nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",Cx->ccs.nl.error_message);
    return _FAILURE_;
  }
 
  if (spectra_init(&Cx->ccs.pr,&Cx->ccs.ba,&Cx->ccs.pt,&Cx->ccs.pm,&Cx->ccs.nl,&Cx->ccs.tr,&Cx->ccs.sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",Cx->ccs.sp.error_message);
    return _FAILURE_;
  }

    //printf("Cosmology Successfully Initialized.\n");

    free(pfc.name);      
    free(pfc.value);     
    free(pfc.filename);  
    free(pfc.read);

    return _SUCCESS_;
}


/**
 *  Free the memory allocated to CLASS cosmology structure
 * 
 * @param Cx    Input: pointer to Cosmology structure
 * @return the error status
 */
int CL_Cosmology_free(struct Cosmology *Cx)
{
  
  if (spectra_free(&Cx->ccs.sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",Cx->ccs.sp.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&Cx->ccs.nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",Cx->ccs.nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&Cx->ccs.pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",Cx->ccs.pm.error_message);
    return _FAILURE_;
  }

 if (transfer_free(&Cx->ccs.tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",Cx->ccs.tr.error_message);
    return _FAILURE_;
  }

  //if (bessel_free(&Cx->ccs.bs) == _FAILURE_) {
  //  printf("\n\nError in primordial_free \n=>%s\n",Cx->ccs.bs.error_message);
   // return _FAILURE_;
 // }

  if (perturb_free(&Cx->ccs.pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",Cx->ccs.pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&Cx->ccs.th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",Cx->ccs.th.error_message);
    return _FAILURE_;
  }

  if (background_free(&Cx->ccs.ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",Cx->ccs.ba.error_message);
    return _FAILURE_;
  }
  
  return _SUCCESS_;

}


/**
 * Compute the matter power spectra (in unit of (Mpc)^3) as a function of k (in unit of 1/Mpc) and z, 
 * Setting the switch "mode", to LINEAR or NONLINEAR, we can compute the linear or nonlinear spectrum respectively.
 * 
 * The CLASS spectra_pk_at_k_and_z() and spectra_pk_nl_at_k_and_z, evaluate the matter power spectrum 
 * at a given value of k and z by interpolating in a table of all P(k)'s computed at this z 
 * by spectra_pk_at_z() (when kmin <= k <= kmax), or eventually by using directly the primordial 
 * spectrum (when 0 <= k < kmin): the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @param modes             Input: switch to decide whether to compute linear or nonlinear spectrum
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * 
 * @return the double value of matter power spectrum
 */
double Pk_dlnPk(struct Cosmology *Cx, double k, double z, int mode)
{

  double pk,pk_cb;
  double result;

  double *pk_ic    = (double*) calloc(Cx -> ccs.sp.ic_ic_size[Cx -> ccs.sp.index_md_scalars], sizeof(double));
  double *pk_cb_ic = (double*) calloc(Cx -> ccs.sp.ic_ic_size[Cx -> ccs.sp.index_md_scalars], sizeof(double)); 

  if(k< PS_KMIN || k>  PS_KMAX){
    printf("Error in PS: the requested value of k (%12.6e) exceeds the tabulation.\nReturning 0.0", k);
    return 0.0;
  } 
  else if (k>= PS_KMIN && k<= PS_KMAX){ /* make sure that the requested k-value is within the range set for CLASS*/
    if (mode == LPOWER){
      spectra_pk_at_k_and_z(&Cx -> ccs.ba, &Cx -> ccs.pm, &Cx -> ccs.sp, k, z, &pk, pk_ic, &pk_cb, pk_cb_ic); 
    }
    else if (mode == NLPOWER){
      spectra_pk_nl_at_k_and_z(&Cx->ccs.ba, &Cx->ccs.pm, &Cx->ccs.sp, k, z, &pk, &pk_cb); 
    }
  }
       
  result = pk;

  free(pk_ic);
  free(pk_cb_ic);

  return result;
} 


/**
 * Read in the linear power spectrum, used to set the initial conditions of Hidden-Valley sims
 *   
 * Input k is in unit of 1/Mpc. First convert it to h/Mpc, and also convert the final matter power spectrum in unit of (Mpc/h)^3
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @param mode              Input: switch to decide whether to evaluate the interpolator of the power spectrum or free the interpolator
 * @return the HV linear matter power spectrum
 */

double Pk_dlnPk_HV(struct Cosmology *Cx, double k, double z, int mode)
{

  static gsl_interp_accel   *pk_accel_ptr;
  static gsl_spline         *pk_spline_ptr;
  FILE *pk_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;
  if(first == 1){
    char pk_filename[FILENAME_MAX];
    sprintf(pk_filename,"/Volumes/Data/Documents/Git/LIM_PS_HM/Output/matter_L/pklin_1.0000_HVsims.txt");

    int nlines = 0;
    nlines     = count_lines_in_file(pk_filename);

    double *k_in, *pk_lin, *log_k, *log_pk;
    k_in    = make_1Darray(nlines);
    log_k   = make_1Darray(nlines);
    pk_lin  = make_1Darray(nlines);
    log_pk  = make_1Darray(nlines);
    pk_file = fopen(pk_filename,"r");

    if(pk_file==NULL){
      printf("Failed to open the file with k-values");
      exit(1);
    }


    char line[MAXL];
    int err;
    while(fgets(line, sizeof line, pk_file) != NULL )
    { 
      if(*line == '#')  continue; 
      for(i=0;i<nlines;i++){
        err = fscanf(pk_file,"%lg %lg\n",&k_in[i],&pk_lin[i]); 
        // printf("%12.6e %12.6e \n",k_in[i],pk_lin[i] );
        log_k[i] = log(k_in[i]);
        log_pk[i] = log(pk_lin[i]);
      }
    }
    fclose(pk_file);

    pk_accel_ptr  = gsl_interp_accel_alloc();
    pk_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pk_spline_ptr,log_k,log_pk,nlines);

    kmax = k_in[nlines-1];
    kmin = k_in[0];

    free(k_in);
    free(log_k);
    free(pk_lin);
    free(log_pk);
    
    first = 0;
  }   


  double kk = k/gb.h;
  double pk, log_pk, logk;
  if(kk<kmin || kk >kmax) {
    pk = 0.;
  }
  else{
    logk   = log(kk);
    log_pk = gsl_spline_eval(pk_spline_ptr,logk,pk_accel_ptr);
    pk     = exp(log_pk);
  }
    
  // double growth;  
  // if(z == 0.)
  //   growth = 1. ;
  // else if(z == 2.)
  //   growth = 0.41806802;
  // else if(z == 4.)
  //   growth = 0.2537648;

  double pkz, growth2;
  if(mode == LPOWER){
    growth2 = pow(growth_D(Cx, z),2.);
    pkz     = growth2 * pk/pow(gb.h,3.);  
  }  
  else if (mode == CLEANUP){
        gsl_interp_accel_free(pk_accel_ptr);
        gsl_spline_free(pk_spline_ptr);
  }


  return pkz;
    
}


/**
 * Compute the transfer function for different species depending on the switch "mode", which 
 * can be set to cdm, baryons or total matter transfer function.
 * 
 * CLASS function spectra_tk_at_k_and_z() routine evaluates the matter transfer functions at a given
 * value of k and z by interpolating in a table of all \f$ T_i(k,z)\f$'s
 * computed at this z by spectra_tk_at_z() (when kmin <= k <= kmax).
 * Returns an error when k<kmin or k > kmax. 
 *   
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @param mode              Input: switch to decide for which species we want to get the transfer function
 * @return the transfer function
 */

double Mk_dlnMk(struct Cosmology *Cx, double k, double z, int mode)
{
  double *tk;
  double t_cdm = 0.0;
  // double t_m = 0.0;

  class_alloc(tk,sizeof(double)* Cx -> ccs.sp.tr_size,Cx -> ccs.errmsg);


  if (k>= PS_KMIN && k<= PS_KMAX){
    if(mode == TRANS){
        spectra_tk_at_k_and_z(&Cx -> ccs.ba, &Cx -> ccs.sp, k , z, tk);  ////This class function calculates the linear transfer function for an arbitrary wavenumber and redshift
        t_cdm = -tk[Cx -> ccs.sp.index_tr_delta_cdm];  
        // t_m = -tk[Cx -> ccs.sp.index_tr_delta_tot];   
    }
    else if(mode == DER){
      t_cdm =0.;
    }         
  }

  free(tk);

  return t_cdm; 
  //return t_m; 


}


/** 
 * The integrand function passed to qags integrator to compute the variance of the matter density
 *   
 * @param x                 Input: integration variable
 * @param par               Input: integration parmaeters
 * @return value of the integrand 
 */
double sig_sq_integrand(double x, void *par)
  {
    double f=0;
    double result = 0;
    
    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);


    double k = exp(x);
    struct Cosmology *Cx = pij.p1;
    double z = pij.p4;
    double R = pij.p5;

    result = 1./(2.*pow(M_PI,2.))* pow(k,3) * pow(window_rth(k,R),2.)* Pk_dlnPk(Cx,k,z,LPOWER);


    return result;

  } 


/** 
 * Compute variance of smoothed matter density fluctuations. 
 * The function sig_sq_integrand() defines the integrand and sig_sq() computes the k-integral
 *   
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the variance
 */
double sig_sq(struct Cosmology *Cx, double z, double R)
{
  double result=0., error=0.;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000000);

  struct integrand_parameters2 par; 

  double kmin = log(gb.PS_kmin);
  double kmax = log(gb.PS_kmax);

  gsl_function F;
  F.function = &sig_sq_integrand;
  F.params = &par;

  par.p1  = Cx;
  par.p4  = z;
  par.p5  = R;

  gsl_integration_qags(&F,kmin,kmax,0.0,1.0e-3,10000000,w,&result,&error);
  gsl_integration_workspace_free (w);
    

  return result;

}


/** 
 * Compute the logarithmic derivative of the variance of smoothed matter density fluctuations w.r.t. smoothing scale
 *   
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the log-derivative of variance
 */

double der_lnsig_sq(struct Cosmology *Cx, double z, double R)
{

  double inc  = 0.01;
  double Rmin = (1.-inc)*R;
  double Rmax = (1.+inc)*R;

  double f =  (log(sig_sq(Cx,z,Rmax)) - log(sig_sq(Cx,z,Rmin)))/(2.*inc*R);

  return f; 

}


/** 
 * The integrand function passed to qags integrator to compute the variance of the unsmoothed matter density
 *   
 * @param x                 Input: integration variable
 * @param par               Input: integration parmaeters
 * @return value of the integrand 
 */
double sigma0_sq_integrand(double x, void *par)
{
  
  struct integrand_parameters2 pij;
  pij = *((struct integrand_parameters2 *)par);

  double k = exp(x);

  struct Cosmology *Cx = pij.p1;
  double z             = pij.p4;

  double result = 1./(2.*pow(M_PI,2.))* pow(k,3.) * Pk_dlnPk(Cx,k,z,LPOWER);

  return result;

} 


/** 
 * Compute variance of unsmoothed matter density fluctuations. 
 * The function sigma0_integrand() defines the integrand and sigma0_sq() computes the k-integral
 *   
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @return the unsmoothed variance
 */
double sigma0_sq(struct Cosmology *Cx, double z, double kmax)  ///kmax is in unit of 1/Mpc
{

	double result, error;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

	struct integrand_parameters2 par; 

	double kmin = log(0.0001*gb.h);

	gsl_function F;
	F.function = &sigma0_sq_integrand;
	F.params   = &par;

	par.p1 = Cx;
	par.p4 = z;

  gsl_integration_qags(&F,kmin,log(kmax),0.0,1.0e-4,1000000,w,&result,&error);
  gsl_integration_workspace_free (w);
     
  return result;

	
}

/**
 * Compute the growth factor D(k,z) which is scale-indep if mode_nu = NUM, and scale-dep if mode_nu = MASS
 * The scale-dep growth is calculated by taking the ratio of the transfer function at redshift z and zero.
 * The scale-indep growth is computed by CLASS directly
 * The switch "mode" can be set to CDM, BA, TOT to return the growth factor of 
 * cdm, baryon and total matter. 
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @return the growth factor, can be k-dep (ex. with nonzero neutrino mass)
 */

double growth_D(struct Cosmology *Cx, double z)
{ 

  double tau;
  int last_index; ///junk
  double * pvecback;

  pvecback = (double*) calloc(Cx -> ccs.ba.bg_size,sizeof(double));

  class_call(background_tau_of_z(&Cx -> ccs.ba,z,&tau),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);
  class_call(background_at_tau(&Cx -> ccs.ba,tau,Cx -> ccs.ba.long_info,Cx -> ccs.ba.inter_normal,&last_index,pvecback),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);

  double Dz = pvecback[Cx -> ccs.ba.index_bg_D];

  free(pvecback);

  return Dz;
}



/**
 * Compute the scale-dependant linear growth rate f(k,z) (i.e the velocity growth factor)
 * by taking numerical derivative of the scale_dep_growth_D() function f(k,a) = d ln D(k,a)/d ln a.
 * The switch "mode" can be set to CDM, BA, TOT to return the growth factor of the corresponding matter component.
 *
 * This is a useful function when constraining physics that induces scale-dependant growth 
 * such as massive neutrinos. 
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @return the growth rate, can be k-dep (ex. with nonzero neutrino mass)
 */

double growth_f(struct Cosmology *Cx, double z)
{ 

  double tau;
  int last_index; ///junk
  double * pvecback;

  pvecback = (double*) calloc(Cx -> ccs.ba.bg_size,sizeof(double));

  class_call(background_tau_of_z(&Cx -> ccs.ba,z,&tau),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);
  class_call(background_at_tau(&Cx -> ccs.ba,tau,Cx -> ccs.ba.long_info,Cx -> ccs.ba.inter_normal,&last_index,pvecback),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);

  double fz = pvecback[Cx -> ccs.ba.index_bg_f];


  free(pvecback);

  return fz;
}


/**
 * Compute the the hubble rate (exactly the quantity defined by CLASS as index_bg_H in the background module). 
 * This function is to a good approximation equal to Hubble(a,Cx) = gb.h*sqrt(Eofa(a,Cx)) 
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @return the hubble parameter
 */

double Hubble(struct Cosmology *Cx, double z)
{

  double tau;
  int last_index; ///junk
  double * pvecback;

  pvecback = (double*) calloc(Cx -> ccs.ba.bg_size,sizeof(double));
  
  class_call(background_tau_of_z(&Cx -> ccs.ba,z,&tau),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);
  class_call(background_at_tau(&Cx -> ccs.ba,tau,Cx -> ccs.ba.long_info,Cx -> ccs.ba.inter_normal,&last_index,pvecback),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);
 
  double H = gb.c*pvecback[Cx -> ccs.ba.index_bg_H];

  free(pvecback);

  return H;
}


/** 
 * Compute the angular diameter distance (exactly the quantity defined by 
 * CLASS as ba.index_bg_ang_distance in the background module). 
 *
 * luminosity distance d_L = (1+z) d_M
 * angular diameter distance d_A = d_M/(1+z)
 * where d_M is the transverse comoving distance, which is equal to comoving distance for flat cosmology 
 * and has a dependance on curvature for non-flat cosmologies, as described in lines 849 - 851
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @return  D_A
 */

double angular_distance(struct Cosmology *Cx, double z)
{
       
  double tau;
  int last_index; ///junk
  double * pvecback;

  pvecback = (double*) calloc(Cx -> ccs.ba.bg_size,sizeof(double));

  class_call(background_tau_of_z(&Cx -> ccs.ba,z,&tau),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);
  class_call(background_at_tau(&Cx -> ccs.ba,tau,Cx -> ccs.ba.long_info,Cx -> ccs.ba.inter_normal,&last_index,pvecback),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);

  double D_A = pvecback[Cx -> ccs.ba.index_bg_ang_distance];

  free(pvecback);

  return D_A;
}



/** 
 * Compute the comoving radial distance  
 *
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @return the double value D_c
 */

double comoving_radial_distance(struct Cosmology *Cx, double z)
{
  double tau;
  int last_index; ///junk
  double * pvecback;

  ///For a flat cosmology, comoving distance is equal to conformal distance. This pieace of code is how 
  ///the comving distance for flat and nonflat cases are computed. Chnage the expression of D_A below
  ///According to this if considering non-flat cosmology. 
  // if (pba->sgnK == 0) comoving_radius = pvecback[pba->index_bg_conf_distance];
  // else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(pba->K);
  // else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(-pba->K);

  pvecback = (double*) calloc(Cx -> ccs.ba.bg_size,sizeof(double));

  class_call(background_tau_of_z(&Cx -> ccs.ba,z,&tau),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);
  class_call(background_at_tau(&Cx -> ccs.ba,tau,Cx -> ccs.ba.long_info,Cx -> ccs.ba.inter_normal,&last_index,pvecback),
            Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);

  double D_c = pvecback[Cx -> ccs.ba.index_bg_conf_distance];

  free(pvecback);

  return D_c;

}

/** 
 * Compute the critical density in unit of M_sun/Mpc^3 
 *   
 * @param Cx                Input: pointer to Cosmology structure
 * @param z                 Input: redshift to compute the spectrum
 * @return the double value of rho_c
 */

double rhoc(struct Cosmology *Cx, double z)
{

  double a =1./(1.+z);
  double Hz2 = pow(Hubble(Cx, z),2.)*pow(10.,6.); ///E (a) = H(a)^2/H0^2
  double G = 6.67408*pow(10.,-11.)*3.24078*pow(10.,-23.);    ///G is in unit of m^3 kg^-1 s^-2, conversion factor from m to Mpc
  double f2 = 3.*Hz2/(8.*M_PI*G) *1./(1.99*pow(10.,30.));     ///To convert to solar mass

  return f2;
}


/** 
 * Compute the Lagrangian radius of halos in unit of  1/Mpc^3 , fixing z=0
 *   
 * @param Cx                Input: pointer to Cosmology structure
 * @param h_mass            Input: halo mass in unit of solar mass
 * @return R_s
 */

double R_scale(struct Cosmology *Cx, double M)  
{

  double f       = 0 ;
  double omega_m0 = (Cx->cosmo_pars[3L]+Cx->cosmo_pars[4L]); 
  double rho_m   = omega_m0 *rhoc(Cx, 0.);

  f = pow(3.*M/(4.*M_PI*rho_m),1./3.);

  return f;
}


/** 
 * Compute the comoving virial radius of halos in unit of  1/Mpc^3, which is defined as the radius at which 
 * the average density within this radius is Delta X rho_c 
 *   
 * @param Cx      Input: pointer to Cosmology structure
 * @param M       Input: halo mass in unit of solar mass
 * @return R_vir
 */

double R_vir(struct Cosmology *Cx, double M)  
{

  double rho_c = rhoc(Cx, 0);
  double Delta = 200.;
  
  double f = pow(3.*M/(4.*M_PI*rho_c*Delta),1./3.);

  return f;

}


/** 
 * Compute the cold dark matter concentration-mass relation
 *   
 * @param M      Input: halo mass in unit of solar mass
 * @param z      Input: redshift of interest
 * @return the cdm concentration
 */

double concentration_cdm(double M, double z)
{
  /*Farnik, where did you take this expression from?*/
  // double alpha = -0.061094;
  // double belta = 0.370894;
  // double B = 4.710019;

  // double result = B / pow(1.0 + z, belta) * pow(M/1.e14, alpha);

  //I will use the Duffy et al 2010 function
  double A      = 7.85;
  double alpha  = -0.081;
  double beta   = -0.71;
  double Mpivot =  2.e12/gb.h;   //in unit of M_sun
  double result = A*pow(M/Mpivot,alpha) * pow(1.+z, beta);

  // double mstar = 7.4e12;
  // double result = 9./(1.+z)*pow(M/mstar,-0.13);

  return result;
}


/** 
 * Compute the NFW halo profile in Fourier space, given by Eq. 3.7 of 2004.09515
 * The profile is normalized to unity at k->0, (see fig 3 of 1003.4740)
 *   
 * @param Cx     Input: pointer to Cosmology structure
 * @param k      Input: wavenumber in unit of 1/Mpc
 * @param M      Input: halo mass in unit of solar mass
 * @param z      Input: redshift of interest
 * @return the nfw profile
 */

double nfw_profile(struct Cosmology *Cx, double k, double M, double z)  
{
  
  double rvir  = R_vir(Cx, M); /* the co-moving virial radius_from_mass */
  double c     = concentration_cdm(M, z);
  double rs    = rvir / c;
  double f_fac = 1./(log(1.+c) - c/(1.+c));
  double rho_s = M/(4.*M_PI*pow(rs,3.)) * f_fac;  ///rho_s is computed by enforcing int dr r^2 u(r) = 1

  double part_1 = sin(k*rs) * (gsl_sf_Si((1.+c)*k*rs) - gsl_sf_Si(k*rs));
  double part_2 = cos(k*rs) * (gsl_sf_Ci((1.+c)*k*rs) - gsl_sf_Ci(k*rs));
  double part_3 = - sin(c*k*rs)/((1.+c)*k*rs);
  double trig   = part_1 + part_2 + part_3;
  
  double result = 4.*M_PI * pow(rs,3.) * rho_s/M * trig;  

  return result;

}

/** 
 * Fourier transform of top-hat window in real space
 * 
 * @param k                 Input: wavenumber in unit of 1/Mpc
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the window function
 */

double window_rth(double k, double R)
{
  double f = 0.0; 

  f = 3.0 * gsl_sf_bessel_j1(k * R)/(k * R);
  // f =  3.0 * (sin(k*R)/pow((k*R),3) - cos(k*R)/pow(k*R,2));

  return f;
}


/** 
 * Derivative w.r.t. smoothing scale of the Fourier transform of top-hat window in real space
 * 
 * @param k                 Input: wavenumber in unit of 1/Mpc
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the derivative of the window function
 */
double derR_window_rth(double k, double R)
{
  double f = 0.0; 

  f = 1./(k*pow(R, 2.)) *(1.5 * k * R * gsl_sf_bessel_j0(k * R) - 4.5 * gsl_sf_bessel_j1(k * R) - 1.5 * k * R *  gsl_sf_bessel_j2(k * R));
  //f = 3.*((3.*cos(k*R))/(pow(k,2.)*pow(R,3.)) - (3.*sin(K*R))/(pow(k,3.)*pow(R,4.)) + sin(k*R)/(k*pow(R,2.)))
  
  return f;
}

/** 
 * Top-hat window in Fourier space
 * 
 * @param k                 Input: wavenumber in unit of 1/Mpc
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the window function
 */
double window_kth(double k, double R)
{
  double f = 0.0;

  if (k<= 1./R)
    f = 1.0; 
  else 
    f =0.0;

  return f;
}


/** 
 * Gaussian window
 * 
 * @param k                 Input: wavenumber in unit of 1/Mpc
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the window function
 */
double window_g(double k, double R)
{
  double f =0.0;

  f = exp(-pow(k * R,2.0)/ 2.0); 

  return f; 
}


/** 
 * Derivative w.r.t smoothing scale of Gaussian window
 * 
 * @param k                 Input: wavenumber in unit of 1/Mpc
 * @param R                 Input: smoothing scale in unit of Mpc
 * @return the derivative of the window function
 */
double derR_logwindow_g(double k, double R)
{
  double f =0.0;

  f = - 0.2116 * pow(k,2.) * R ;

  return f; 
}







