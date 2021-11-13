
/** @file wnw_split.c Documented wiggle-nowiggle split based on 3d Gaussian filter in linear k, 
 * and using the Eisentein-Hu wiggle-no wiggle template arXiv:astro-ph/9709112 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * The algorithm closely follows Ref. arXiv:1509.02120 by Vlah et al. (described in Appendix A)
 *
 * The following function will be called from other modules:
 * -# pk_Gfilter_nw()
 */



#include "header.h"

/**
 * Compute the nowiggle component of linear matter power spectrum using 3d Gaussian filter
 * Computing the nowiggle component involves calculating an integral (Eq. A3 of Vlah et al)
 * Below, pk_nw_integrand()is the corresponding integrand and pk_Gfilter_nw() is the integrator 
 * 
 * @param Cx            Input: pointer to Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc
 * @param k0            Input: smallest value of k, i.e. the largest scale
 * @return broadband component in unit of (Mpc)^3
 */
double pk_Gfilter_nw(struct Cosmology *Cx, double k, double k0)    
{ 
    struct integrand_parameters2 par; 
    par.p1 = Cx;
    par.p4 = k;
    par.p5 = k0;

    double a       = 0.25; 
    double logqmin =  log10(k) - 4.* a;
    double logqmax =  log10(k) + 4.* a;

    double result, error;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    gsl_function F;
    F.function = &pk_nw_integrand;
    F.params = &par;

    gsl_integration_qags(&F,logqmin, logqmax,0.0,1.0e-3,1000000,w,&result,&error);
    gsl_integration_workspace_free (w);
   
    double pk0     = Pk_dlnPk(Cx, k0, 0., LPOWER);
    double ratio   = 1./(sqrt(2.*M_PI)* a) * result; 
    double out     =  EH_PS_nw(Cx, k, k0, pk0) * ratio;

    // printf("%12.6e %12.6e %12.6e \n", k, out, Pk_dlnPk(Cx, k, z, LPOWER));

    return out;
}


/**
 * Integrand to compute the nowiggle matter power spectrum
 * 
 * @param x            Input: integration variable, k-values
 * @param par          Input: integration parameters
 * @return integrand to be used in pk_Gfilter_nw() function         
 */

double pk_nw_integrand(double x, void *par)  ///integration variable x = logq
{  
    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);

    struct Cosmology *Cx = pij.p1;
    double k             = pij.p4;
    double kf0           = pij.p5;

    double logq  = x;
    double q     = pow(10.,logq);
    double logk  = log10(k);
    double a     = 0.25; 

    double pkf0  = Pk_dlnPk(Cx,kf0, 0., LPOWER);
    double ratio = Pk_dlnPk(Cx,q, 0., LPOWER)/EH_PS_nw(Cx,q, kf0, pkf0);  
    double out   = ratio * exp(-1./(2.*pow(a,2.))*pow(logk-logq,2.));

    return out;
}        


/**
 * Compute the Eisentein-Hu approximate wiggle component of linear matter power spectrum
 *
 * @param Cx            Input: pointer to Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc
 * @param k0            Input: smallest value of k, i.e. the largest scale
 * @param pk0           Input: value of the power spectrum at the largest scale
 * @return P_w(k) in unit of (Mpc)^3
 */

double EH_PS_w( struct Cosmology *Cx, double k, double k0, double pk0)
{
    double h    = Cx->cosmo_pars[2];
    double ns   = Cx->cosmo_pars[1];
    double kh   = h*k;
    double kh0  = h*k0;
    // double norm = pk0/(pow(k0,ns)*pow(T(kh0,Cx),2.));
    // double out  = pow(k,ns) * pow(T(kh,Cx),2.) * norm;
    double norm = pk0/(pow(k0,ns)*pow(T(Cx,k0),2.));
    double out  = pow(k,ns) * pow(T(Cx,k),2.) * norm;

    return out;
}    


/**
 * Compute the Eisentein-Hu approximate nowiggle component of linear matter power spectrum
 *
 * @param Cx            Input: pointer to Cosmology structure
 * @param k           Input: wavenumber in unit of 1/Mpc
 * @param k0          Input: smallest value of k, i.e. the largest scale
 * @param pk0         Input: value of the power spectrum at the largest scale
 * @return P_nw(k) in unit of (Mpc)^3
 */

double EH_PS_nw( struct Cosmology *Cx, double k,double k0,double pk0) 
{
    double h    = Cx->cosmo_pars[2];
    double ns   = Cx->cosmo_pars[1];
    double kh   = h*k;
    double kh0  = h*k0;
    // double norm = p0/(pow(k0,ns)*pow(T(kh0,Cx),2.));
    // double out  = pow(k,ns) * pow(T0(kh,Cx),2.) * norm;
    double norm = pk0/(pow(k0,ns)*pow(T(Cx,k0),2.));
    double out  = pow(k,ns) * pow(T0(Cx,k),2.) * norm;

    return out;
}


/**
 * Compute the no-baryon transfer function given in Eq. 29 of EH ref
 *
 * @param Cx          Input: pointer to Cosmology structure
 * @param k           Input: wavenumber in unit of 1/Mpc
 * @return value of nor-baryon transfer fit
 */

double T0( struct Cosmology *Cx, double k)
{
    double h     = Cx->cosmo_pars[2];
    double ombh2 = pow(h,2.) * Cx->cosmo_pars[3];
    double omch2 = pow(h,2.) * Cx->cosmo_pars[4];
    double om0h2 = ombh2 + omch2;
    double om0   = om0h2/pow(h,2.);
    double theta = 2.728/2.7;    //OBBE-FIRAS value 

    double s     = 44.5 * log(9.83/om0h2)/sqrt(1.+10.*pow(ombh2,3./4.));    ///approximate sound speed given in Eq. (26) of EH
    double AG    = 1. - 0.328*log(431.*om0h2)*ombh2/om0h2 + 0.38*log(22.3*om0h2)*pow(ombh2/om0h2,2.);                                            
    double Gamma = om0 * h * (AG + (1. - AG)/(1.+pow(0.43*k*s,4.)));
    double q     = k/h *pow(theta,2.)/Gamma ;
    double L0    = log(2.*exp(1.) + 1.8 * q);
    double C0    = 14.2 + 731./(1.+62.5*q);
    double out   = L0/(L0+C0*pow(q,2.));

    return out;

}

/** 
 * Compute the total baryon+CDM transfer function
 *
 * @param Cx          Input: pointer to Cosmology structure
 * @param k           Input: wavenumber in unit of 1/Mpc
 * @return value of baryon+cdm transfer function
 */

double T(struct Cosmology *Cx, double k)
{
    double h     = Cx->cosmo_pars[2];
    double ombh2 = pow(h,2.) * Cx->cosmo_pars[3];
    double omch2 = pow(h,2.) * Cx->cosmo_pars[4];
    double om0h2 = ombh2 + omch2;
    double om0   = om0h2/pow(h,2.);
    double theta = 2.728/2.7;    //OBBE-FIRAS value 

    double HH0     = 1.e3*1.e2*h/299792458.;   ///H0 value devided by the speed of light
    double zeq     = 2.5e4 * om0h2 * pow(theta,-4.);
    double keq     = sqrt(2.*om0*pow(HH0,2.)*zeq);
    double k_silk  = 1.6*pow(ombh2,0.52)*pow(om0h2,0.73)*(1.+pow(10.4*om0h2,-0.95));  ////in 1/Mpc
    
    double beta_node = 8.41*pow(om0h2,0.435);
    double s         = 44.5 * log(9.83/om0h2)/sqrt(1.+10.*pow(ombh2,3./4.));    ///approximate sound speed given in Eq. (26) of EH
    double st        = s/pow(1.+ pow(beta_node/(k*s),3.),1./3.);

    double bb1   = 0.313 * pow(om0h2,-0.419) * (1.+0.607*pow(om0h2,0.674));
    double bb2   = 0.238 * pow(om0h2,0.223);
    double zd    = 1291. * pow(om0h2,0.251)/(1.+0.659*pow(om0h2,0.828))*(1.+bb1*pow(ombh2,bb2));
    double Rd    = 31.5 * ombh2*pow(theta,-4.)*pow(zd/1.e3,-1.);
    double y     = (1.+zeq)/(1.+zd);
    double G     = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)));

    double a1     = pow(46.9*om0h2,0.670)*(1.+pow(32.1*om0h2,-0.532));
    double a2     = pow(12.0*om0h2,0.424)*(1.+pow(45.*om0h2,-0.582));
    double alphac = pow(a1,-ombh2/om0h2) * pow(a2,-pow(ombh2/om0h2,3.));
    double b1     = 0.944*pow(1.+pow(458.*om0h2,-0.708),-1.);
    double b2     = pow(0.395 * om0h2,-0.0266);
    double betac  = pow(1. + b1*(pow(omch2/om0h2,b2)-1.),-1.);

    double alphab = 2.07*keq*s*pow(1.+Rd,-3./4.)*G;
    double betab  = 0.5 + ombh2/om0h2 + (3.-2.*ombh2/om0h2) * sqrt(pow(17.2*om0h2,2.)+1.); 
    double f      = 1./(1.+pow(k*s/5.4, 4.));

    double Tb = (Tt0(Cx,k,1.,1.)/(1.+pow(k*s/5.2,2.)) + alphab/(1.+pow(betab/(k*s),3.)) * exp(-pow(k/k_silk,1.4)))* gsl_sf_bessel_j0(k*st); ///Eq. 21 of EH ref.
    double Tc = f*Tt0(Cx,k,1.,betac) + (1.-f) *Tt0(Cx,k,alphac,betac);       ///Eq. 17 of EH ref                                    

    double out = ombh2/om0h2 * Tb + omch2/om0h2 * Tc;

    return out;
}    


/**
 * Compute the function defined in Eq. 19 of EH ref, which will be used to compute the fit for CDM transfer function in Eq. 17.
 *
 * @param Cx          Input: pointer to Cosmology structure
 * @param k           Input: wavenumber in unit of 1/Mpc. 
 * @param x1          Input: betac AM:WHAT WAS THIS VARIABLE???
 * @param x2          Input: betac AM:WHAT WAS THIS VARIABLE???
 * @return            value of the function
 */

double Tt0(struct Cosmology *Cx, double k, double x1, double x2) ///x1 = alphac, x2 = betac in Eq. 19.
{
    double h     = Cx->cosmo_pars[2];
    double ombh2 = pow(h,2.) * Cx->cosmo_pars[3];
    double omch2 = pow(h,2.) * Cx->cosmo_pars[4];
    double om0h2 = ombh2 + omch2;
    double theta = 2.728/2.7;    //OBBE-FIRAS value 

    double qq  = k*pow(theta,2.)*pow(om0h2,-1.);
    double C   = 14.2/x1+386./(1.+69.9*pow(qq,1.08));
    double out = log(exp(1.)+1.8*x2*qq)/(log(exp(1.)+1.8*x2*qq)+C*pow(qq,2.));
 
    return out;
}     

