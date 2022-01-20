
/** @file ps_line_hm.c Documented halo-model computation of line power spectrum, including clustering and stochastic contributions beyond Poisson limit 
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 * This module has two main functions: 
 *      - PS_line_HM() to compute clustering (1- and 2-halo terms). The 2-halo term, includes nonlinear corrections to halo power spectrum arising from nonlinearities of matter fluctuations and halo biases.
 *      - PS_shot_HM() to compute the stochastic contrubutions beyond Poisson shot noise (see arXiv:1706.08738)
 * 
 * The other functions in these modules are utilities for computing the above two main functions. 
 *  
 * In summary, the following functions can be called from other modules:
 * -# PS_line_HM()  computes the 1loop-HM line power spectrum 
 * -# PS_shot_HM()  computes the full line power spectrum, beyond poisson limit 
 * -# mhmc()        computes th corrections to mass integration of halo-model matter power spectrum
 * -# HM_1h2h()     performs the mass integraks for computing 1- and 2-halo terms of line-line, line-matter and matter-matter power spectra.
 * -# b22_ls()      computes the large-scale limit of P_b2b2 loop which behaves like a constant and so contributes to the shot noise. 
 */

#include "header.h"
struct globals gb;


/**
 * Compute the clustering contribution to the line power spectrum using halo-model. 
 * If nfw=1, the dependance of the power spectrum on the halo profile is neglected. Otherwise, NFW halo profile is assumed
 * 
 * @param Cx            Input: pointer to Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc. 
 * @param z             Input: redshift
 * @param M_min         Input: minimum halo mass for mass integrals
 * @param mode_mf       Inpute: theoretical model of halo mass function to use. 
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * @param line_type     Inpute: name of the line to compute.
 *                              It can be set to CII, CO10, CO21, CO32, CO43, CO54, CO65 
 * @param line_id       Inpute: id of the line to be considered. 
 * @return P_clust(k)        
 */
double PS_line_HM(struct Cosmology *Cx, double k, double z, double M_min, long mode_mf,  long line_type, int line_id)
{
      long cleanup = 0;

      double f = 0.;
      double a = 1./(1.+z);
      
      double k_B      = 1.38064852*1.e-16;  ///Boltzmann constant in unit of erg K^-1
      double len_conv = 3.086*1.e19;   ////conversion factor from Mpc to km
      double L_sun    = 3.83 * 1.e33;  ///in unit of erg/s
      
      int J = 0;
      double nu_line;
      if(line_type == CII){
        nu_line = 1902. * pow(10.,9.); ///CII
        J = 0;
      }
      else
      { 
        if(line_type == CO10){ 
          J = 1;
        }
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

        nu_line = J * 115.27 * pow(10.,9.);  ////in unit of Hz for CO
      }

      double scatter    = 0.37;
      double scatter_1h = p_sig_shot(scatter);
      
      double Hz  = Hubble(Cx, z);
      double fac = 1.e6 * 1./(8.*M_PI*k_B*Hz) * pow(gb.c/nu_line,3.) * pow(1.+z,2.)* pow(len_conv,-2.)*L_sun;

      double cs2      = 1.;
      double khat     = 1. * gb.h;

      double mom1 = mass_moment1(Cx,z,M_min,mode_mf,line_type);
      double Tbar = Tbar_line(Cx,line_id,z);  ///to plot the power spectrum in units of micro K^2 Mpc^3

      // Comoute constant linear and quadratic line biases. When nfw=1, these would be k-independant
      double *result = make_1Darray(3);
      HM_1h2h(Cx, k, z, M_min, mode_mf, line_type, LINE, result);  /// in unit of M_sun/Mpc^3
      double b1_line  = result[1]/mom1;  //the line line bias which is a mass integrated luminosity-HMF weighted halo bias
      double b2_line  = result[2]/mom1;  //the line quadratic bias which is a mass integrated luminosity-HMF weighted halo bias 
      double bg2_line =  -2./7.  * (b1_line-1.);
      double btd_line =  23./42. * (b1_line-1.);  
      // printf("biases %12.6e %12.6e %12.6e %12.6e %12.6e \n", z, b1_line, b2_line, bg2_line, btd_line);
      
      //Now compute the 1loop bias and matter corrections to the 2-halo term, using the biases computed above. 
      double * ps_hloops = make_1Darray(6);
      double * ps_mloops = make_1Darray(2);
      Compute_G_loops(Cx, k, z, WIR, MATTER, GFILTER, ps_mloops);
      Compute_G_loops(Cx, k, z, WIR, HALO, GFILTER, ps_hloops);

      double pb1b2       = pow(Tbar,2.)*b1_line * b2_line * ps_hloops[0];
      double pb1bg2      = pow(Tbar,2.)*2. * b1_line * bg2_line * ps_hloops[1];
      double pb22        = pow(Tbar,2.)* 0.25 * pow(b2_line, 2.) * ps_hloops[2];
      double pbg22       = pow(Tbar,2.)*pow(bg2_line, 2.)  * ps_hloops[3];
      double pb2bg2      = pow(Tbar,2.)*b2_line * bg2_line * ps_hloops[4];
      double pb1b3nl     = pow(Tbar,2.)*2. * b1_line * (bg2_line + 2./5. * btd_line) * ps_hloops[5];
      double pline_loops = pb1b2 + pb1bg2 + pb22 + pbg22 + pb2bg2 + pb1b3nl;

      double pm_lin_IR   = pm_IR_LO(Cx, k, z, GFILTER);
      double pm_1loop_IR = pm_IR_NLO(Cx, k, z, GFILTER);
      double pm_ct       = - 2. * cs2 * pow(k, 2.)/(1.+pow(k/khat,2.)) * pm_lin_IR;
      double pline_tot   = pow(Tbar,2.)* pow(b1_line, 2.) * (pm_1loop_IR + pm_ct) + pline_loops;

      //If you wanted to plot the linear and 22 and 13 matter loops seperately, uncomment the following 3 lines
      double p22   = pow(Tbar,2.)* pow(b1_line, 2.) * ps_mloops[0];
      double p13   = pow(Tbar,2.)* pow(b1_line, 2.) * ps_mloops[1];
      double plin  = pow(Tbar,2.)* pow(b1_line, 2.) * PS(Cx, k, z) ;
      
      //Now compute the total line power specturm, including the 1- and 2-halo terms.
      double pk_1h = scatter_1h * pow(fac,2.) * result[0];
      double pk_2h = pline_tot;
      double total = pk_1h + pk_2h;

      // printf("line loops %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",k, plin, p22, p13, pow(Tbar,2.)*pow(b1_line, 2.) * pm_1loop_IR, pow(Tbar,2.)*pow(b1_line, 2.)*pm_ct, pb1b2, pb1bg2, pb22, pbg22, pb2bg2, pb1b3nl) ;


      //Save the output to files. Note that I am multiplying the power spectrum componenet by h^3, since I want to plot the P(k) in unit of (Mpc/h)^3
      FILE *fp1, *fp2;
      char filename1[FILENAME_MAX], filename2[FILENAME_MAX];
      sprintf(filename1,"%s/line/ps/no_nfw/pk_hm_J%d_z%d.txt", gb.output_dir,J,(int)z);
      sprintf(filename2,"%s/line/ps/no_nfw/pk_2h_comps_J%d_z%d.txt", gb.output_dir,J,(int)z);
      fp1 = fopen(filename1, "ab");
      fp2 = fopen(filename2, "ab");
      fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e\n",k/gb.h, pk_1h*pow(gb.h,3.), pk_2h*pow(gb.h,3.), total*pow(gb.h,3.));
      fprintf(fp2,"%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",k/gb.h, plin*pow(gb.h,3.), p22*pow(gb.h,3.), p13*pow(gb.h,3.), pow(Tbar,2.)* pow(b1_line, 2.) * pm_1loop_IR*pow(gb.h,3.), pow(Tbar,2.)*pow(b1_line, 2.)*pm_ct*pow(gb.h,3.), pb1b2*pow(gb.h,3.), pb1bg2*pow(gb.h,3.), pb22*pow(gb.h,3.), pbg22*pow(gb.h,3.), pb2bg2*pow(gb.h,3.), pb1b3nl*pow(gb.h,3.)) ;

      //Close the output files
      fclose(fp1);
      fclose(fp2);

      //Free the memeory allocated to array.
      free(result);
      free(ps_hloops);
      free(ps_mloops);

      return total;
}


/**
 * Compute the shot noise contributions, including corrections beyond poisson limit (see 1706.08738 for more details)
 * If nfw=1, the dependance of the power spectrum on the halo profile is neglected. Otherwise, NFW halo profile is assumed
 * 
 * @param Cx            Input: Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc. 
 * @param z             Input: redshift
 * @param M_min         Input: minimum halo mass for mass integrals
 * @param input         inpute: an array of input values with 4 values, Tave_line, b1_line, pb22_ls, line_shot, rhom_bar
 * @param mode_mf       Inpute: theoretical model of halo mass function to use. 
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * @param line_type     Inpute: name of the line to compute.
 *                              It can be set to CII, CO10, CO21, CO32, CO43, CO54, CO65 
 * @return P_stoch(k)     
 */
double PS_shot_HM(struct Cosmology *Cx, double k, double z, double M_min, double *input, long mode_mf, long line_type)
{
      double linematter[1], matter[2], line[3];

      double a        = 1./(1.+z);  
      double k_B      = 1.38064852*1.e-16;  ///Boltzmann constant in unit of erg K^-1
      double len_conv = 3.086*1.e19;        ////conversion factor from Mpc to km
      double L_sun    = 3.83 * 1.e33;       ///in unit of erg/s
      
      int J = 0;
      double nu_line;
      if(line_type == CII){
        nu_line = 1902. * pow(10.,9.); ///CII
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
        nu_line = J * 115.27 * pow(10.,9.);  ////in unit of Hz for CO
      }

      double scatter    = 0.37;
      double scatter_1h = p_sig_shot(scatter);
      double scatter_2h = p_sig_Tbar(scatter);
      double Hz         = Hubble(Cx,z);
      double fac        = 1.e6 * 1./(8.*M_PI*k_B*Hz) * pow(gb.c/nu_line,3.) * pow(1.+z,2.)* pow(len_conv,-2.)*L_sun;

      ///Since the following quantities do not depend on k, I am  computing them once and pass them as input to this function
      double Tave_line = input[0];  ///to plot the power spectrum in units of micro K^2 Mpc^3
      double b1_line   = input[1];
      double pb22_ls   = input[2]; 
      double line_shot = input[3];
      double rhom_bar  = input[4];
   
      //If accounting for nfw in computing the stochastic terms, you should. replace the b1_line and pb22_ls above, with a k-dep one 
      //First, compute the line 1halo term. 
      HM_1h2h(Cx, k, z, M_min, mode_mf, line_type, LINE, line);  /// in unit of M_sun/Mpc^3
      double pk_1h    = scatter_1h * pow(fac,2.) * line[0];
     
      double mom1     = mass_moment1(Cx,z,M_min,mode_mf,line_type);
      // double b1_line  = line[1]/mom1;  //the line line bias which is a mass integrated luminosity-HMF weighted halo bias
      // double b2_line  = line[2]/mom1;  //the line line bias which is a mass integrated luminosity-HMF weighted halo bias
      // double pb22_ls  = 0.25 * pow(b2_line,2.) * pow(Tave_line,2.) * b22_ls(Cx,z);

      //Now compute the correction to the mass integrals of pkm_1h when setting lower mass limit
      double Ms   = 1.e6;  //Lower mass limit for the halo model integration
      double hm_corrs[2];
      mhmc(Cx,z,mode_mf,hm_corrs);
      // double nfw   = nfw_profile(Cx, k, Ms, z);
      double nfw = 1.;

      double p1h_cor = pow(nfw,2.) * Ms/rhom_bar * hm_corrs[0]; //correction to matter 1halo term, Eq. 2.27
      // double p2h_cor = nfw * hm_corrs[1]; //correction to matter 2halo term, we dont use it here

      //Next, compute the corrections to poisson noise from matter-matter and lin-matter 
      HM_1h2h(Cx, k, z, M_min, mode_mf, line_type, LINEMATTER, linematter);  /// in unit of M_sun/Mpc^3
      HM_1h2h(Cx, k, z, M_min, mode_mf, line_type, MATTER, matter);  /// in unit of M_sun/Mpc^3
   
      double pklm_1h  = -2. * scatter_2h * b1_line * Tave_line * fac/rhom_bar * linematter[0]; 
      double pkm_1h   = pow(b1_line * Tave_line, 2.) * (pow(rhom_bar,-2.) * matter[0] + p1h_cor);

      //Finally, put them all together to compute the stochastic power spectrum.
      double pk_stoch = line_shot + pklm_1h + pkm_1h + pb22_ls;

      printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",z, k/gb.h,  pk_1h * pow(gb.h,3.), line_shot * pow(gb.h,3.), pklm_1h * pow(gb.h,3.), pkm_1h* pow(gb.h,3.), pb22_ls* pow(gb.h,3.), pk_stoch* pow(gb.h,3.));

      // FILE *fp1;
      // char filename1[FILENAME_MAX];
      // sprintf(filename1,"%s/line/stoch/no_nfw/pk_stoch_hm_J%d_z%d.txt", gb.output_dir, J, (int)z);
      // fp1 = fopen(filename1, "ab");
      // fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",k/gb.h,  pk_1h * pow(gb.h,3.), line_shot * pow(gb.h,3.), pklm_1h * pow(gb.h,3.), pkm_1h* pow(gb.h,3.), pb22_ls* pow(gb.h,3.), pk_stoch* pow(gb.h,3.));
      // fclose(fp1);

      return pk_stoch;    
}


/**
 * The integrand function passed passed to Cuhre integration routine of CUBA library to compute the corrections to mass integration. 
 *   
 * When computing the matter power spectrum using halo-model, the mass integrations for 1- and 2-loop terms
 * get contributions from halos of all masses. For numerical computation, we need to impose a lower and upper integration 
 * limit. While the result of the integration are not sensitive to the upper bound (due to the fact that the mass function drops rapidly at high M_h)
 * the choice of the lower bound affects the results. We can compute the leading order corrections to the integral that are
 * accurate up to (k R_s)^2.  (see App. A of arXiv:1511.02231 for more details.)
 * 
 * @param ndim       Input: Dimensionality of the domain of integration
 * @param x          Input: An array of integration variables
 * @param ncomp      Input: Dimensionality of the integrand function
 * @param ff         Input: Array of values of the integrand of dimension fdim
 * @param p          Input: integration parmaeters
 * return the error status
 */
static int mhmc_integ(const int *ndim,
                       const cubareal x[],
                       const int *ncomp,
                       cubareal ff[],
                       void *p)         
{
    
      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);

      struct Cosmology *Cx = pij.p1;
      double z             = pij.p5;
      double logMmin       = pij.p6;
      double logMmax       = pij.p7;
      long mode_mf         = pij.p13;

      double logM    = (x[0] * (logMmax - logMmin) + logMmin); 
      double cos     = (2. * x[1] - 1.); 
      double M       = exp(logM);

      double *bias_arr = make_1Darray(4);
      halo_bias(Cx, M, z, mode_mf, bias_arr);
      double b1   = bias_arr[0];
      
      double MF  = mass_func(Cx, M, z, mode_mf);

      ff[0] = 0.5 * 2. * (logMmax - logMmin) * pow(M,2.) * MF ;
      ff[1] = 0.5 * 2. * (logMmax - logMmin) * pow(M,2.) * MF * b1;

      
      return 0;

} 



/**
 * Compute the corrections to mass integration of HM matter power spectrum. 
 *   
 * @param Cx            Input: Cosmology structure
 * @param z             Input: redshift
 * @param mode_mf       Inpute: theoretical model of halo mass function to use. 
 *                              It can  be set to Press-Schecter (PSC), sheth-Tormen (ST), Tinker (TR) 
 * @param results       Output: a 2d array of the integration results, 
 *                              - results[0]: correction to 1-halo term, 
 *                              - result[1]: correcrions to 2-halo term assuming linear halo bias
 * @return void
 */
void mhmc(struct Cosmology *Cx, double z, long mode_mf, double *result)
{
      extern struct globals gb;
      struct integrand_parameters2 par;

      double AbsErr = 0.0;        // Required absolute error (0.0 for none)
      double RelErr = 1.e-2;      // Required relative error, higher precision makes the code extremely slow

      par.p1  = Cx;
      par.p5  = z;
      par.p13 = mode_mf;
      par.p6 =  log(1.e6);  
      par.p7 =  log(1.e17);
      
      int ncomp = 2, ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 2.e8, key = 7;
      int nregions, neval;
      
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
      double integral[ncomp];

      Cuhre(ndim,ncomp, mhmc_integ, &par, nvec,
            RelErr, AbsErr, verbose | last, mineval, maxeval, key,
            NULL, NULL, &nregions, &neval, fail,integral, error, prob);
        
      double rhom_bar = (Cx->cosmo_pars[3] + Cx->cosmo_pars[4]) *rhoc(Cx, 0); 
      for(int i =0;i<2;i++){
          result[i] = 1. - 1./rhom_bar * integral[i];
      } 
      
      return;
}


/**
 * The integrand function passed passed to Cuhre integration routine to compute 1- and 2-halo integrals
 * 
 * @param ndim       Input: Dimensionality of the domain of integration
 * @param x          Input: An array of integration variables
 * @param ncomp      Input: Dimensionality of the integrand function
 * @param ff         Input: Array of values of the integrand of dimension fdim
 * @param p          Input: integration parmaeters
 * return the error status
 */
static int HM_1h2h_integ(const int *ndim,
                       const cubareal x[],
                       const int *ncomp,
                       cubareal ff[],
                       void *p)         
{
      double f      = 0;
      double result = 0;

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p); 

      struct Cosmology *Cx = pij.p1;
      double k             = pij.p4;
      double z             = pij.p5;
      double logMmin       = pij.p6;
      double logMmax       = pij.p7;
      long mode_mf         = pij.p13;
      long line_type       = pij.p14;
      long mode_hm         = pij.p15;

      double logM    = (x[0] * (logMmax - logMmin) + logMmin); 
      double cos     = (2. * x[1] - 1.); 
      double M       = exp(logM);

      double *bias_arr = make_1Darray(4);
      halo_bias(Cx, M, z, mode_mf, bias_arr);
      double b1   = bias_arr[0];
      double b2   = bias_arr[1];

      double MF  = mass_func(Cx, M, z, mode_mf);
      // double nfw = nfw_profile(Cx, k, M, z);
      double nfw = 1.;
      double lum = luminosity(M, z, line_type);
      
      ///we assume the profile of both matter and line are NFW
      if(mode_hm == LINE){
            ///integrand of line 1halo term 
            ff[0] = 0.5 * 2. * (logMmax - logMmin) * M * MF * pow(nfw,2.) * pow(lum,2.);
           ///integrand of 2halo term proportional to b1, the linear local-in-matter halo bias 
            ff[1] = 0.5 * 2. * (logMmax - logMmin) * M * MF * nfw * lum * b1;
            // // ///integrand of 2halo term proportional to b2, the quadratic local-in-matter bias 
            ff[2] = 0.5 * 2. * (logMmax - logMmin) * M * MF * nfw * lum * b2;           
      }      
      else if(mode_hm == LINEMATTER){
            /// integrand of 1halo term of line-matter cross-spectrum 
            ff[0] = 0.5 * 2. * (logMmax - logMmin) * pow(M,2.) * MF * pow(nfw,2.) * lum;   
      }
      else if(mode_hm == MATTER){
            ff[0] = 0.5 * 2. * (logMmax - logMmin) * pow(M,3.) * MF * pow(nfw,2.);
            ff[1] = 0.5 * 2. * (logMmax - logMmin) * pow(M,2.) * MF * nfw * b1;

      }      
      // printf("integrands %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",k,z,M,ff[0],ff[1],ff[2]);

      free(bias_arr);

      return 0;
}


/**
 * Compute the mass integrals needed for 1- and 1-halo line, line-matter and matter power spectrum
 * If nfw=1, the dependance of the power spectrum on the halo profile is neglected. Otherwise, NFW halo profile is assumed
 * 
 * @param Cx            Input: Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc. 
 * @param z             Input: redshift
 * @param M_min         Input: minimum halo mass for mass integrals
 * @param mode_mf       Input: theoretical model of halo mass function to use. 
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * @param line_type     Input: name of the line to compute.
 *                              It can be set to CII, CO10, CO21, CO32, CO43, CO54, CO65 
 * @param mode_hm       Input: a switch to decide whetehr to compute gthe mass integrations. It can be set to:
 *                              - LINE for line power spectrum, 
 *                              - LINEMATTER for line-matter cross spectrum 
 *                              - MATTER for matter power spectrum
 * @param results       Output: anarray of the integration results. Number of elements varies depending on mode_hm switch:
 *                              - 3 elements if mode_hm = LINE, 
 *                              - 1 element if mode_hm  = LINEMATTER
 *                              - 2 element if mode_hm  = MATTER
 *                              esults[0]: correction to 1-halo term, result[1]: correcrions to 2-halo term assuming linear halo bias
 * @return void
 */
void HM_1h2h(struct Cosmology *Cx, double k, double z, double M_min,
            long mode_mf, long line_type, long mode_hm, double *result)  /// in unit of M_sun/Mpc^3
{
      struct integrand_parameters2 par;

      int    ncomp =0;                               
      double AbsErr = 0.0;        // Required absolute error (0.0 for none)
      double RelErr = 1.e-2;      // Required relative error, higher precision makes the code extremely slow

      par.p1  = Cx;
      par.p4  = k;
      par.p5  = z;
      par.p13 = mode_mf;
      par.p14 = line_type;
      par.p15 = mode_hm;

      if(mode_hm == MATTER){ 
            par.p6 =  log(1.e6);  
            par.p7 =  log(1.e17);
      }
      else{
            par.p6  = log(M_min);
            par.p7  = log(1.e16);
      }

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 2.e8, key = 7;
      int nregions, neval;
      
      if(mode_hm == LINE)
            ncomp = 3; 
      else if(mode_hm == MATTER)
            ncomp = 2;
      else if(mode_hm == LINEMATTER)
            ncomp = 1;    

      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
    
      Cuhre(ndim,ncomp, HM_1h2h_integ, &par, nvec,
            RelErr, AbsErr, verbose | last, mineval, maxeval, key,
            NULL, NULL, &nregions, &neval, fail, result, error, prob);

      // for(int i=0;i<ncomp;i++)
      //     printf("res_HM: %d %12.6e %12.6e %12.6e \n", i, k, z, result[i]);


      return;  
}



/**
 * The integrand function passed passed to Cuhre integration routine to compute large-scale limit of b22 (shot-noise contribution)
 * 
 * @param ndim       Input: Dimensionality of the domain of integration
 * @param x          Input: An array of integration variables
 * @param ncomp      Input: Dimensionality of the integrand function
 * @param ff         Input: Array of values of the integrand of dimension fdim
 * @param p          Input: integration parmaeters
 * return the error status
 */
static int b22_ls_integrand(const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)  
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);   


      double logqmax = log(100.);
      double logqmin= log(1.e-4);

      double logq  = (x[0] * (logqmax - logqmin) + logqmin); 
      double cos   = (2. * x[1] - 1.); 
      double q     = exp(logq); 


      struct Cosmology *Cx = pij.p1;
      double z = pij.p4;

      double plin_IR_q = pm_IR_LO(Cx, q, z, GFILTER);

      ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_q;  

      return 0;
} 


/**
 * Compute the large-scale limit of P_b2b2 loop
 * 
 * @param Cx            Input: Cosmology structure
 * @param z             Input: redshift
 * @return b22_ls
 */
double b22_ls(struct Cosmology *Cx, double z)
{
      struct integrand_parameters2 par; 

      par.p1  = Cx;
      par.p4  = z;


      int ncomp =1;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail;
      double error;
      double prob;
      double result;
      int key = 13;
     

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e8;
      int nregions, neval;
      
      Cuhre(ndim,ncomp, b22_ls_integrand, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, &fail, &result, &error, &prob);
         
      // printf("b2_ls %12.6e %12.5e \n ",result,error)  ;
       
      return result;

}
