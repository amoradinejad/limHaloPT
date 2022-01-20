
/** @file ps_line_pt.c Documented computation of Poisson shot noise and tree-level line power spectrum in real and redshift-space.
 *  
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 *  
 * NOTE TODO: Add the 1loop redshift-space power spectrum of the line. This requires implementing FFTLog, still in progress
 * For the moment we stick to the tree-level expression of line power spectrum in redshift-space. 
 *
 * In summary, the following functions can be called from other modules:
 * -# PS_tot_interloper()     computes power spectrum of interlopers (tree-level with RSD)
 * -# PS_line_real()          computes tree-level line power spectrum in real-space
 * -# PS_line_RSD()           computes the tree-level line power spectrum in redshift-space, as a function of wavenumber and angle w.r.t LOS
 * -# ps_line_multipoles()    computes the redshift-space multipoles of the line power spectrum
 * -# PS_shot()               computes the poisson shot noise
 */

#include "header.h"
struct globals gb;


/**
 * Compute the 3D power spectrum of interloper line power spectrum in unit of micro K^2 Mpc^3
 * 
 * @param Cx            Input: Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc. 
 * @param z_s           Input: redshift of the signal line
 * @param z_inter       Input: redshift of the interloper line
 * @param line_id       Inpute: id of the line to be considered. 
 * @return tree-level P_clust(k)       
 */
double PS_tot_interloper(struct Cosmology *Cx, double k, double mu, double z_s, double z_inter, size_t line_id)
{

  double Tave_line   = Tbar_line(Cx, line_id, z_inter);  ///to plot the power spectrum in units of micro K^2 Mpc^3
  double bias_arr[2];
  line_bias(Cx->Lines[line_id], z_inter, bias_arr);
  double b1z = bias_arr[0];

  ///// Alcock-Paczinski type effect due to interlopers////
  double sine        = sqrt((1.-pow(mu,2.)));
  double A_par       = Hubble(Cx,z_inter)/Hubble(Cx,z_s) * (1.+z_s)/(1.+z_inter);
  double A_perp      = comoving_radial_distance(Cx,z_s)/comoving_radial_distance(Cx,z_inter);
  double k_par       = k * mu *  A_par; 
  double k_perp      = k * sine * A_perp;
  double k_inter     = sqrt(pow(k_par,2.) + pow(k_perp,2.));
  double mu_inter    = k_par/k_inter;

  double sigz        = 0.001 * (1.+z_inter);   ////What is sigz for CO surveys? 
  double sig_FOG     = Cx->cosmo_pars[NPARS-1L] * sqrt(1.+z_inter);
  double sigv        = (1.+z_inter)*sqrt(pow(sig_FOG,2.)/2. + pow(gb.c *sigz,2.));
  double sup_fac     = exp(-pow(k_inter*mu_inter*sigv,2.)/pow(Hubble(Cx,z_inter),2.));
  double inter_vol   = pow(A_perp,2.) * A_par;
      
  double growth_rate = growth_f(Cx, k_inter, z_inter);
  double beta        = growth_rate/b1z;
  double a0_P        = pow(1.+beta*pow(mu_inter,2.),2.);

  double clust       = sup_fac * a0_P * pow(Tave_line,2.) * pow(b1z,2.) * PS(Cx, k_inter, z_inter);
  double shot        = PS_shot(Cx, line_id, z_inter);
  double tot         = inter_vol * (clust + shot);

  return tot; 
}


/**
 * Compute the real-space 3D power spectrum of emission lines in unit of micro K^2 Mpc^3
 * 
 * @param Cx            Input: Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc. 
 * @param z             Input: redshift
 * @param line_id       Inpute: id of the line to be considered. 
 * @return tree-level P_clust(k)       
 */
double PS_line_real(struct Cosmology *Cx, double k, double z, size_t line_id)
{

      double deltac = 1.686;
      double a      = 1./(1.+z); 

      double bias_arr[2];
      line_bias(Cx->Lines[line_id], z, bias_arr);
      double b1z = bias_arr[0];
      double Tave_line = Tbar_line(Cx, line_id, z);  
      double f = pow(Tave_line,2.)*pow(b1z,2.)*PS(Cx, k, z);

      return f; 
}

/**
 * Compute the redshift-space 3D power spectrum of emission lines in unit of micro K^2 Mpc^3 as a function of wavenumber and angle w.r.t. LOS
 * 
 * @param Cx            Input: Cosmology structure
 * @param Cx_ref        Input: Reference cosmology structure, needed for AP effect
 * @param k             Input: wavenumber in unit of 1/Mpc. 
 * @param mu            Inpute: angle w.r.t LOS
 * @param z             Input: redshift
 * @param line_id       Inpute: id of the line to be considered. 
 * @return tree-level P_clust(k,mu)    
 */
double PS_line_RSD(struct Cosmology *Cx, struct Cosmology *Cx_ref, double k, double mu, double z, size_t line_id)
{

      double f      = 0.;
      double deltac = 1.686;
      double a      = 1./(1.+z); 

      double bias_arr[2];
      line_bias(Cx->Lines[line_id], z, bias_arr);
      double b1z = bias_arr[0];
      double Tave_line = Tbar_line(Cx, line_id, z);  ///to plot the power spectrum in units of micro K^2 Mpc^3

      ///// Alcock-Paczinski effect: infering distances from redshifts and positions requires assuming a ref cosmology////
      double k_true  = k * pow((1.-pow(mu,2.))*pow(angular_distance(Cx_ref,z)/angular_distance(Cx,z),2.)\
                         + pow(mu,2.)* pow(Hubble(Cx,z)/Hubble(Cx_ref,z),2.),1./2.);    
      double mu_true = k/k_true * mu * Hubble(Cx,z)/Hubble(Cx_ref,z);

      double growth_rate = growth_f(Cx, k, z);
        
      double sigz    = 0.001 * (1.+z);   ////What is sigz for CO surveys? 
      double sig_FOG = Cx->cosmo_pars[NPARS-1L] * sqrt(1.+z);
      double sigv    = (1.+z)*sqrt(pow(sig_FOG,2.)/2. + pow(gb.c *sigz,2.));
      double sup_fac = exp(-pow(k_true*mu_true*sigv,2.)/pow(Hubble(Cx_ref,z),2.));
      double AP_fac  = pow(angular_distance(Cx_ref,z)/angular_distance(Cx,z),2.) * Hubble(Cx,z)/Hubble(Cx_ref,z);
            
      double beta    = growth_rate/b1z;
      double a0_P    = pow(1.+beta*pow(mu_true,2.),2.);

      f = AP_fac * sup_fac * a0_P * pow(Tave_line, 2.)*pow(b1z, 2.) * PS(Cx, k, z);
      
      return f; 
}


/**
 * The integrand function passed to hcubature integrator to compute the line power spectrum multipoles

 * @param nd         Input: Dimensionality of the domain of integration
 * @param x          Input: integration variable
 * @param p          Input: integration parmaeters
 * @param fdim       Input: Dimensionality of the integrand function
 * @param fvalue     Input: Array of values of the integrand of dimension fdim 
 * return the error status
 */
int ps_line_multipoles_integrand(unsigned    ndim,           // Number of dimensions in the domain space -- number of dim we're integrating over
                        const double         *x,             // The point at which the integrand is evaluated
                        void                 *p,             // Pointer to a structure that holds the parameters
                        unsigned             fdim,           // Number of dimensions that the integrand return
                        double               *fvalue         // Array of values of the integrand of dimension fdim
                        )
{
      struct integrand_parameters2 pij;
      unsigned i,j;

      pij = *((struct integrand_parameters2 *)p);

      struct Cosmology *Cx     = pij.p1;
      struct Cosmology *Cx_ref = pij.p2;
      double k                 = pij.p4;
      double z                 = pij.p5;
      size_t line_id           = pij.p22;
      int ell                  = pij.p19; 

      double mu     = x[0];
      double result =(2.* ell + 1.)/2.* gsl_sf_legendre_Pl(ell,mu) * PS_line_RSD(Cx, Cx_ref, k, mu, z, line_id);

      *fvalue =  result;
      
      return  0;
}


/**
 * Compute the multipole moments of redshift-space power spectrum of emission lines in unit of micro K^2 Mpc^3, integrating over the angle w.r.t LOS, weighted by 
 * 
 * @param Cx            Input:  Cosmology structure
 * @param Cx_ref        Input:  Reference cosmology structure, needed for AP effect
 * @param k             Input:  wavenumber in unit of 1/Mpc. 
 * @param z             Input:  redshift
 * @param line_id       Inpute: id of the line to be considered. 
 * @param ell           Inpute: the multipole
 * @return P_ell(k)       
 */
double ps_line_multipoles(struct Cosmology *Cx, struct Cosmology *Cx_ref, double k, double z, size_t line_id, int ell)
{
      struct integrand_parameters2 par;

      unsigned fdim = 1; 
      unsigned ndim = 1;                                    // Dimensionality of the domain of integration
      double xmin[1], xmax[1];                        // Integration limits: these are arrays of dimension dim
      size_t maxEval = 0;                                   // Maximum number of integrand evaluations (0 for none)
      double AbsErr=0.0;                                    // Required absolute error (0.0 for none)
      double RelErr=1.e-3;                            // Required relative error (1.0e-2)                               
      double result=0.;                         // Final result
      double error=0.;                          // Error estimate on the result
  
      error_norm norm = ERROR_INDIVIDUAL;

      xmin[0] = -1.;
      xmax[0] = 1.;
      
      par.p1  = Cx;
      par.p2  = Cx_ref;
      par.p4  = k;
      par.p5  = z;
      par.p22 = line_id;
      par.p19 = ell;

      hcubature(fdim, ps_line_multipoles_integrand, &par, ndim,xmin,xmax,maxEval,AbsErr,RelErr,norm,  &result, &error);   

      return result;
}

/**
 * Compute the Poisson shot noise in unit of micro K^2 Mpc^3
 * 
 * @param Cx            Input:  Cosmology structure
 * @param z             Input:  redshift
 * @param line_id       Inpute: id of the line to be considered. 
 * @return P_poisson    
 */
 double PS_shot(struct Cosmology *Cx, double z, size_t line_id)  ///in unit of micro K^2 Mpc^3
{
  double f =0.;
  double a = 1./(1.+z);

  double zz[4] = {DO_NOT_EVALUATE,z,DO_NOT_EVALUATE,DO_NOT_EVALUATE};
  double res[4];
  Line_evaluate(Cx->Lines[line_id], zz, res);
  double mass_mom2  = res[1];
  
  double k_B      = 1.38064852*1.e-16;  ///Boltzmann constant in unit of erg K^-1
  double len_conv = 3.086*1.e19;   ////conversion factor from Mpc to km
  double nu_line  = Cx->Lines[line_id]->line_freq;

  double L_sun  = 3.83 * 1.e33;  ///in unit of erg/s
  double fduty  = 0;
  double sig_CO = 0.37;
  double scatter_shot = p_sig_shot(sig_CO);
  
  double Hz  = Hubble(Cx,z);  //I get 2.036142e+02 , Farnik has 204.45512791569357
  // double Hz  = 204.45512791569357;  // 

  double fac = 1.e6 * 1./(8.*M_PI*k_B*Hz) * pow(gb.c/nu_line,3.) * pow(1.+z,2.)* pow(len_conv,-2.)*L_sun;
  
  f = scatter_shot * pow(fac,2.) * mass_mom2 ;

  return f;
}            


