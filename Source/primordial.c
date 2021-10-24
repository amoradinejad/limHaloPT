#include "header.h"
struct globals gb;
  

double f_s(int s, double mu_s) /////To be completed once Hayden added the results for higher spins
{
      double f =0;
      if(s==2)
            f = (-(985. - 664. * pow(mu_s,2.) + 16. * pow(mu_s,4.))/576.);
      else if(s==3)
            f = (4800519.0 - 2564556.0    * pow(mu_s,2.) + 96366.0* pow(mu_s,4.) - 64.0* pow(mu_s,6.))/460800.0 ;
      else if(s==4)     
            f = (-221842845.0 + 341268176.0* pow(mu_s,2.) - 39893856.0* pow(mu_s,4.) + 278784.0* pow(mu_s,6.) + 2816.0* pow(mu_s,8.))/12902400.0;
      
      return f;   
}


double Heaviside(double ks, double kl)
{
      double f = 0;
      double sq_fac = 0.1;
      if(kl/ks > sq_fac)
            f = 0;
      else if(kl/ks <= sq_fac )
            f = 1.;

      return f;
}     

double Bispec_zeta(double k1, double k2, double k3, int s, struct Cosmology *Cx, long mode_PNG )
{     
  
      double f = 0;
      double As   = exp(Cx->cosmo_pars[0L])*pow(10.,-10.);
      double ns   = Cx->cosmo_pars[1];
      double kt = k1 + k2 + k3;

      double mu_s, mu12, mu13, mu23;
      if(mode_PNG == HS || mode_PNG == QSF){
            mu_s = Cx->cosmo_pars[6];           
          mu12 = (pow(k3,2.) - pow(k1,2.) - pow(k2,2.))/(2.*k1*k2);
            mu13 = (pow(k2,2.) - pow(k1,2.) - pow(k3,2.))/(2.*k1*k3);
            mu23 = (pow(k1,2.) - pow(k2,2.) - pow(k3,2.))/(2.*k2*k3);
      
            if(mu12>1)
                  mu12 = 1;
            else if(mu12<-1)
                  mu12 = -1;

            if(mu13>1)
                  mu13 = 1;
            else if(mu13<-1)
                  mu13 = -1;

            if(mu23>1)
                  mu23 = 1;
            else if(mu23<-1)
                  mu23 = -1;
      }

      

      double b_angle = 0, b_osc =0, r_s=0, fs=0, phase=0, P_12=0, P_13=0, P_23 =0 ;
      double num = 0, denum =0, fac =0;

      if(mode_PNG == HS){
            fs        = f_s(s, mu_s);
            phase       = 0;
            P_12    = gsl_sf_legendre_Pl(s,mu12);
            P_13    = gsl_sf_legendre_Pl(s,mu13);
            P_23    = gsl_sf_legendre_Pl(s,mu23);
            b_angle = 1./pow(kt,2.*s+1.)\
                        *(P_12/(k3*pow(k1*k2,3.-s)*pow(k1*k2*k3,1.-ns)) * ((2.*s-1.)*((k1+k2)*kt + 2.*s*k1*k2) + pow(kt,2.))\
                      + P_13/(k2*pow(k1*k3,3.-s)*pow(k1*k2*k3,1.-ns)) * ((2.*s-1.)*((k1+k3)*kt + 2.*s*k1*k3) + pow(kt,2.))\
                      + P_23/(k1*pow(k2*k3,3.-s)*pow(k1*k2*k3,1.-ns)) * ((2.*s-1.)*((k2+k3)*kt + 2.*s*k2*k3) + pow(kt,2.)));  
            b_osc   = 1./pow(k1*k2,4.-ns) * pow(k1/k2,1.5) * cosl(mu_s* logl(k1/k2) + phase) * Heaviside(k2,k1)* P_12\
                      + 1./pow(k1*k2,4.-ns) * pow(k2/k1,1.5) * cosl(mu_s* logl(k2/k1) + phase) * Heaviside(k1,k2)* P_12\
                      + 1./pow(k1*k3,4.-ns) * pow(k1/k3,1.5) * cosl(mu_s* logl(k1/k3) + phase) * Heaviside(k3,k1)* P_13\
                      + 1./pow(k1*k3,4.-ns) * pow(k3/k1,1.5) * cosl(mu_s* logl(k3/k1) + phase) * Heaviside(k1,k3)* P_13\
                      + 1./pow(k2*k3,4.-ns) * pow(k2/k3,1.5) * cosl(mu_s* logl(k2/k3) + phase) * Heaviside(k3,k2)* P_23\
                      + 1./pow(k2*k3,4.-ns) * pow(k3/k2,1.5) * cosl(mu_s* logl(k3/k2) + phase) * Heaviside(k2,k3)* P_23;            
            if(s==2){
                  num         = pow(M_PI,4.5) * (1. + 4.*pow(mu_s,2.))*(9. + 4*pow(mu_s,2.));
                  denum       = 256.0*(1.+1./tanh(M_PI*mu_s));
                  fac   = sqrt(2.*(81.+4.*pow(mu_s,2.))/(mu_s*sinh(2.*M_PI*mu_s)));
            }
            else if(s==3)
            {
                  num         = - 5. *pow(M_PI,4.5)*pow((1. + 4.*pow(mu_s,2.))*(9. + 4*pow(mu_s,2.))*(25. + 4.*pow(mu_s,2.)),2.);
                  denum       = 688128.*(225. + 1036. * pow(mu_s,2.) + 2032. * pow(mu_s,4.) + 64.* pow(mu_s,6.)); 
                  fac   = sqrt(2.*(121.+4.*pow(mu_s,2.))/(mu_s*sinh(2.*M_PI*mu_s))) ;
            }
            else if(s==4){
                  num   = pow(M_PI,4.5)* pow((1. + 4.*pow(mu_s,2.))*(9. + 4*pow(mu_s,2.))*(25. + 4.*pow(mu_s,2.))*(49. + 4.*pow(mu_s,2.)),2.);
                  denum       = 31850496.*(11025. + 196144. *pow(mu_s,2.) + 24864. * pow(mu_s,4.) + 5376. * pow(mu_s,6.) + 256. * pow(mu_s,8.));
                  fac   = sqrt(2.*(169.+4.*pow(mu_s,2.))/(mu_s*sinh(2.*M_PI*mu_s))) ;
            }
            r_s =  num/denum * fs * fac ;       
            f   = b_angle + 0.5* r_s *b_osc ; 
      }
      else if(mode_PNG == QSF){
            double A = 8.*k1*k2*k3/pow(k1+k2+k3,3.);
            // f = 54./5.*sqrt(3)* pow(2.*M_PI*M_PI*As,2.)/gsl_sf_bessel_Ynu(mu_s,8./27.) * gsl_sf_bessel_Ynu(mu_s,A)/(pow(k1*k2*k3,3./2.)*pow(k1+k2+k3,3./2.));
            f = 2./5.*pow(gb.kp,2.*(1.-ns))* pow(2.*M_PI*M_PI*As,2.)*pow(3,1.5)/gsl_sf_bessel_Ynu(mu_s,8./27.) * gsl_sf_bessel_Ynu(mu_s,A)/(pow(k1*k2*k3,3./2.)*pow(k1+k2+k3,3./2.));

      }
      else if(mode_PNG == LOCAL)
            f     = 6./5.*pow(2.*M_PI*M_PI*As,2.)*(1./pow(k1*k2,4.-ns) +  1./pow(k1*k3,4.-ns) +  1./pow(k2*k3,4.-ns)); 
      else if(mode_PNG == EQUILATERAL)
            f     = 18./5.*pow(2.*M_PI*M_PI*As,2.)*(- 1./pow(k1*k2,4.-ns) - 1./pow(k1*k3,4.-ns) - 1./pow(k2*k3,4.-ns) - 2./pow(k1*k2*k3,2.*(4.-ns)/3.)\
                  + 1./(pow(k1,(4.-ns)/3.)*pow(k2,2.*(4.-ns)/3.)*pow(k3,4.-ns)) + 1./(pow(k1,(4.-ns)/3.)*pow(k3,2.*(4.-ns)/3.)*pow(k2,4.-ns))\
                  + 1./(pow(k2,(4.-ns)/3.)*pow(k1,2.*(4.-ns)/3.)*pow(k3,4.-ns)) + 1./(pow(k2,(4.-ns)/3.)*pow(k3,2.*(4.-ns)/3.)*pow(k1,4.-ns))\
                  + 1./(pow(k3,(4.-ns)/3.)*pow(k1,2.*(4.-ns)/3.)*pow(k2,4.-ns)) + 1./(pow(k3,(4.-ns)/3.)*pow(k2,2.*(4.-ns)/3.)*pow(k1,4.-ns))); 
      else if(mode_PNG == ORTHOGONAL)
            f     = 18./5.*pow(2.*M_PI*M_PI*As,2.)*(- 3./pow(k1*k2,4.-ns) - 3./pow(k1*k3,4.-ns) - 3./pow(k2*k3,4.-ns) - 8./pow(k1*k2*k3,2.*(4.-ns)/3.)\
                  + 3./(pow(k1,(4.-ns)/3.)*pow(k2,2.*(4.-ns)/3.)*pow(k3,4.-ns)) + 3./(pow(k1,(4.-ns)/3.)*pow(k3,2.*(4.-ns)/3.)*pow(k2,4.-ns))\
                  + 3./(pow(k2,(4.-ns)/3.)*pow(k1,2.*(4.-ns)/3.)*pow(k3,4.-ns)) + 3./(pow(k2,(4.-ns)/3.)*pow(k3,2.*(4.-ns)/3.)*pow(k1,4.-ns))\
                  + 3./(pow(k3,(4.-ns)/3.)*pow(k1,2.*(4.-ns)/3.)*pow(k2,4.-ns)) + 3./(pow(k3,(4.-ns)/3.)*pow(k2,2.*(4.-ns)/3.)*pow(k1,4.-ns))); 
      
      return f;
}



double Bispec_zeta_normalized(double k1, double k2, double k3, int s, struct Cosmology *Cx, long mode_PNG)
{
      double f    = 0;
      double As   = exp(Cx->cosmo_pars[0])*pow(10.,-10.);
      double ns   = Cx-> cosmo_pars[1];
      double fnl  = Cx-> cosmo_pars[5];

      double P_zeta     = 2.*M_PI*M_PI*As/pow(k1,4.-ns);
      double norm       = 5./18.*Bispec_zeta(k1,k1,k1,s,Cx,mode_PNG)/pow(P_zeta,2.);

      f     = fnl*Bispec_zeta(k1,k2,k3,s,Cx,mode_PNG)/norm;
      
      return f;

}

