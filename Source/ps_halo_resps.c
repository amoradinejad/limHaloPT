#include "header.h"
struct globals gb;

double Pk_rspace(double k, double z, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double alpha  = Cx->cosmo_pars[1]; 
	double fnl    = Cx->cosmo_pars[2];

	double deltac = 1.686;
	double bz  = 2. * deltac * (b1 - 1.);
	double result  = (pow(b1, 2.) + 2.*fnl*b1*bz/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	return result;		  
}

double Pk_zspace(double k, double mu, double z, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double alpha  = Cx->cosmo_pars[1];
	double fnl    = Cx->cosmo_pars[2];

	double deltac = 1.686;
	double a   = 1./(1.+z);
	double f   = growth(a,Cx, DERGROWTH); 
	double bz  = 2. * deltac * (b1 - 1.);
	double result  = (pow(b1 + f*pow(mu,2.), 2.) + 2. * fnl * bz * (b1 + f * pow(mu,2.))/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	return result;		  
}


/////Real- and redshit-sapce responses of the short-scale power spectrum to large-scale modes of density and gravitational potential
double resp_Pk_rspace(double k, double z, double R, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double alpha  = Cx->cosmo_pars[2]; 
	double Delta0 = Cx->cosmo_pars[3]; 
	double fnl    = Cx->cosmo_pars[4];

	double Phi0 =  pow(sig_sq(z, R, PHI, Cx),1./2.);

	double deltac = 1.686;
	double bz  = 2. * deltac * (b1 - 1.);
	double pg  = (pow(b1, 2.) + 2.*fnl*b1*bz/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	double result = pg + resp_Delta0_rspace(k,z,nbar,Cx) * Delta0\
			  + resp_phi_rspace(k,z,nbar,Cx) * Phi0;

	return result;		  
}

double resp_Pk_zspace(double k, double mu, double z, double R, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double alpha  = Cx->cosmo_pars[3]; 
	double Delta0 = Cx->cosmo_pars[4]; 
	double Delta2 = Cx->cosmo_pars[5];
	double fnl    = Cx->cosmo_pars[6];

	double Phi0 =  pow(sig_sq(z, R, PHI, Cx),1./2.);

	double deltac = 1.686;
	double a   = 1./(1.+z);
	double f   = growth(a,Cx, DERGROWTH); 
	double bz  = 2. * deltac * (b1 - 1.);
	double pg  = (pow(b1 + f*pow(mu,2.), 2.) + 2. * fnl * bz * (b1 + f * pow(mu,2.))/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	double result = pg + resp_Delta0_zspace(k,mu,z,nbar,Cx) * Delta0\
			  + resp_Delta2_zspace(k,mu,z,nbar,Cx) * Delta2\
			  + resp_phi_zspace(k,mu,z,nbar,Cx) * Phi0;

	return result;		  
}

/////Real-sapce response functions to large-scale density and gravitational potential modes
double resp_Delta0_rspace(double k, double z, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double b2     = Cx->cosmo_pars[1];
	double alpha  = Cx->cosmo_pars[2]; 
	double fnl    = Cx->cosmo_pars[4];

	double deltac = 1.686;
	double bz  = 2. * deltac * (b1 - 1.);
	double bzd = 2. * (deltac * (b2 + 13./21. * (b1 - 1.)) - b1 + 1.);
	double pg  = (pow(b1,2.) + 2. * fnl *  bz * b1/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	double result = (47./21. * pow(b1,2.) + 2. * b1 * b2 - pow(b1,2.)/3. * Pk_dlnPk(k,z,Cx,DER)) * Pk_dlnPk(k,z,Cx,POWER)\
		   + fnl * (26./21. * b1 * bz + 2. * b1 * bzd + 2. * b2 * bz - 2./3. * b1 * bz\
		   	* (Pk_dlnPk(k,z,Cx,DER) - Mk_dlnMk(k,z,Cx,DER))) * Pk_dlnPk(k,z,Cx,POWER)/Mk_dlnMk(k,z,Cx,TRANS)\
		   + b1*(1.+alpha)/nbar - 2. * b1 * pg;
	
	return result;	   
}

double resp_phi_rspace(double k, double z, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double b2     = Cx->cosmo_pars[1];
	double alpha  = Cx->cosmo_pars[2]; 
	double fnl    = Cx->cosmo_pars[4];

	double deltac = 1.686;
	double bz  = 2. * deltac * (b1 - 1.);
	double bzd = 2. * (deltac * (b2 + 13./21. * (b1 - 1.)) - b1 + 1.);
	double pg  = pow(b1,2.) * Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;
	
	double result   = fnl * (4. * pow(b1,2.) + 2. * b1 * bzd) * Pk_dlnPk(k,z,Cx,POWER) + fnl*bz*(1.+alpha)/nbar - 2.*fnl*bz*pg;
	
	return result;	   
}

/////Redshift-sapce response functions to large-scale density, tidal field and gravitational potential modes
double resp_Delta0_zspace(double k, double mu, double z, double nbar, struct Cosmology *Cx)
{	
	double b1     = Cx->cosmo_pars[0];
	double b2     = Cx->cosmo_pars[1];
	double alpha  = Cx->cosmo_pars[3]; 
	double fnl    = Cx->cosmo_pars[6];

	double deltac = 1.686;
	double a   = 1./(1.+z);
	double f   = growth(a,Cx, DERGROWTH); 
	double bz  = 2. * deltac * (b1 - 1.);
	double bzd = 2. * (deltac * (b2 + 13./21. * (b1 - 1.)) - b1 + 1.);
	double pg  = (pow(b1 + f*pow(mu,2.), 2.) + 2. * fnl * bz * (b1 + f * pow(mu,2.))/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	double result = (-1./3. * Pk_dlnPk(k,z,Cx,DER) * (1.+f*pow(mu,2.)) * pow(b1+f*pow(mu,2.),2.)\
		+ 1./21. * (b1+f*pow(mu,2.)) * (f*pow(mu,2.)*(42.*b1-7.*f+31.)+7.*b1*f+47.*b1+42.*b2+28.*pow(f,2.)*pow(mu,4.))) * Pk_dlnPk(k,z,Cx,POWER)\
		+ fnl * (2./3. * b1 * bz * (Mk_dlnMk(k,z,Cx,DER) - Pk_dlnPk(k,z,Cx,DER)) * (1.+f*pow(mu,2.)) * (b1+f*pow(mu,2.))\
		+ 2./21. *(f*pow(mu,2.) * (bz*(14.*f*pow(mu,2.)+ 5.) + 21.*bzd)\
		+ b1 * bz*(7.*f*(3.*pow(mu,2.)+1.) + 13.) + 21.*b2*bz + 21. * b1 * bzd))*Pk_dlnPk(k,z,Cx,POWER)/Mk_dlnMk(k,z,Cx,TRANS)\
		+ (b1+f/3.)*(1.+alpha)/nbar - 2.*(b1+f/3.)*pg;

	return result;		  		  
}

double resp_Delta2_zspace(double k, double mu, double z, double nbar, struct Cosmology *Cx)
{

	double b1     = Cx->cosmo_pars[0];
	double b2     = Cx->cosmo_pars[1];
	double bs2    = Cx->cosmo_pars[2];
	double alpha  = Cx->cosmo_pars[3]; 
	double fnl    = Cx->cosmo_pars[6];

	double deltac = 1.686;
	double a   = 1./(1.+z);
	double f   = growth(a,Cx, DERGROWTH); 
	double bz  = 2. * deltac * (b1 - 1.);
	double bzd = 2. * (deltac * (b2 + 13./21. * (b1 - 1.)) - b1 + 1.);
	double pg  = (pow(b1 + f*pow(mu,2.), 2.) + 2. * fnl * bz * (b1 + f * pow(mu,2.))/Mk_dlnMk(k,z,Cx,TRANS))* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	double result = (2./21.* (b1 + f *pow(mu,2.))*(pow(mu,2.)* (12.* b1 + 42.* bs2 - f *(7.* f + 8.))\
 	     		+ 7.* b1* f - 4.* b1 - 14.* bs2 + 4.* f * (7.* f + 6.)* pow(mu,4.))\
 	     		- 1./3. * Pk_dlnPk(k,z,Cx,DER)*((2.*f + 3.) * pow(mu,2.)-1.)*pow(b1+f*pow(mu,2.),2.)) * Pk_dlnPk(k,z,Cx,POWER)\
			+ fnl * (4./21. * bz *(pow(mu,2.)*(6.*b1+21.*bs2-4.*f) + 7.*b1*f -2.*b1 - 7.*bs2 + 2.*f*(7.*f + 6.)*pow(mu,4.))\
				+ 2./3. * (Mk_dlnMk(k,z,Cx,DER) - Pk_dlnPk(k,z,Cx,DER)) * bz*((2.*f + 3.) * pow(mu,2.) - 1.) * (b1 + f * pow(mu,2.)))*Pk_dlnPk(k,z,Cx,POWER)/Mk_dlnMk(k,z,Cx,TRANS)\
			+ 2./3. * f * (1.+alpha)/nbar - 4./3. * f * pg;
		
	return result;		  
}

double resp_phi_zspace(double k, double mu, double z, double nbar, struct Cosmology *Cx)
{
	double b1     = Cx->cosmo_pars[0];
	double b2     = Cx->cosmo_pars[1];
	double alpha  = Cx->cosmo_pars[3]; 
	double fnl    = Cx->cosmo_pars[6];

	double deltac = 1.686;
	double a   = 1./(1.+z);
	double f   = growth(a,Cx, DERGROWTH); 
	double bz  = 2. * deltac  * (b1 - 1.);
	double bzd = 2. * (deltac * (b2 + 13./21. * (b1 - 1.)) - b1 + 1.);
	double pg  = pow(b1 + f*pow(mu,2.), 2.)* Pk_dlnPk(k,z,Cx,POWER) + (1. + alpha)/nbar;

	double result = 2.*fnl*((b1+f*pow(mu,2.)) *(f*pow(mu,2.)*(bz+2.) + bzd + 2.*b1)) * Pk_dlnPk(k,z,Cx,POWER)\
			  + fnl*bz*(1.+alpha)/nbar - 2.*fnl*bz*pg;

	return result;
}




