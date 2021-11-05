
/** @file survey_specs.c Documented computation of some survey-related functions
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 * In summary, the following functions can be called from other modules:
 * -# shell_volume()            computes the comoving volume of a survey covering redshift up to z
 * -# kmin_val()                computes the fundumental k-mode of a given redshift shell
 * -# kmax_Brent_solver()       computes the kmax value such that kmax(z=0) = 0.15 h/Mpc
 *              
 */

#include "header.h"
struct globals gb;



/**
 * Compute the comoving volume of a survey covering redshift up to z
 * 
 * @param Cx            Input: Cosmology structure
 * @param z             Input: redshift
 * @param fsky          Input: sky-coverage of teh survey 
 * @return the comoving z-shell volume
 */

double shell_volume( struct Cosmology *Cx,double z, double fsky) // The volume of the shell in (Mpc)^3
{
    double f=0.; 

    f = 4. * fsky* M_PI/3. * pow(comoving_radial_distance(Cx,z),3.);

    return f;
}


/**
 * Compute the size of fundumental mode corresponding to the comoving volume enclosed in a redshift bin [zmin,zmax]
 * 
 * @param Cx            Input: Cosmology structure
 * @param zmin          Input: minimum redshift
 * @param zmin          Input: maximum redshift
 * @param fsky          Input: sky-coverage of teh survey 
 * @return kmin
 */

double kmin_val( struct Cosmology *Cx, double zmin, double zmax, double fsky)
{
    double f =0;

    f = 2. * M_PI *pow(3./(4. *M_PI) * (shell_volume(Cx,zmax,fsky)-shell_volume(Cx,zmin,fsky)),-1./3.); // kmin in units of 1/Mpc

    return f;

}


/**
 * Compute the maximum k-value used in Fisher forecast at each redshift bin. We follow Giannantonio et al. to for determining kmax, and use gsl Brent solver to solve for kmax in each redshift bin.
 * The goal is to compute the kmax such that at z=0, the variance of the matter fluctations has a fixed value, for instance 0,36. This can be achieved 
 * by solving Eq. 40 of Giannantonio: sigma^2(z) = int_kmin^kmax(z) dk k^2/(2pi^2) P_m(k,z) = 0.36. Instead of fixing sigma^2 to 0.36, I chose the variance 
 * such that kmax(z=0) = 0.15 h/Mpc. This corresponds to the variance of ~0.33 at z=0 . In the forecast, I additionally always impose kmax<0.3 h/Mpc cut. 
 * 
 * @param Cx            Input: Cosmology structure
 * @param z             Input: redshift of interest
 * @return kmax
 */

double kmax_Brent(double kmax, void *params)
{
    double f;
    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)params);

    f = sigma0_sq(pij.p1,pij.p4,kmax) - pij.p6;

    return f;   
}

double kmax_Brent_solver( struct Cosmology *Cx, double z)
{
    int status;

    int iter = 0, max_iter = 100000;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    double k_lo = 0.001*Cx->cosmo_pars[2], k_hi = 200.*Cx->cosmo_pars[2];

    gsl_function F;

    struct integrand_parameters2 par; 
    double R_s = R_scale(Cx,gb.h_m);


    par.p1 = Cx;
    par.p4  = z;
    par.p6 = 3.305751e-01;  //variance of unsmoothed density field such that kmax(z=0) = 0.15 h/Mpc
    ///in short paper we used, 3.631872e-01 ;  
        
    F.function = &kmax_Brent;
    F.params = &par;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, k_lo, k_hi);

    do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      k_lo = gsl_root_fsolver_x_lower (s);
      k_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (k_lo, k_hi,0., 1e-5);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return r;
}

