#ifndef COSMOLOGY_H
#define COSMOLOGY_H

int Cosmology_init(struct Cosmology *Cx, double pk_kmax, double pk_zmax, 
            int nlines, int * line_types, size_t npoints_interp, double M_min, long mode_mf);
int Cosmology_free(struct Cosmology *Cx);
int CL_Cosmology_initilize(struct Cosmology *Cx, double pk_kmax, double z_pk);
int CL_Cosmology_free(struct Cosmology *Cx);

double Pk_dlnPk(struct Cosmology *Cx, double k, double z , int mode);
double Mk_dlnMk(struct Cosmology *Cx, double k, double z, int mode);
double Pk_dlnPk_HV(struct Cosmology *Cx, double k, double z , int mode);

double sigma0_sq_integrand(double x, void *p);
double sigma0_sq(struct Cosmology *Cx, double z, double kmax);
double sig_sq_integrand(double x, void *par);
double sig_sq(struct Cosmology *Cx, double z, double R);
int der_sig_sq_integrand(unsigned  nd,       // Number of dimensions in the domain space -- number of dim we're integrating over
    const double  *x,           // The point at which the integrand is evaluated
    void          *p,           // Pointer to a structure that holds the parameters
    unsigned      fdim,           // Number of dimensions that the integrand return
    double        *fvalue         // Array of values of the integrand of dimension fdim
    );
double der_sig_sq(struct Cosmology *Cx, double z, double R);
double der_lnsig_sq(struct Cosmology *Cx, double z, double R);

double growth_D(struct Cosmology *Cx, double z);
double growth_f(struct Cosmology *Cx, double z);
double scale_indep_growth_D(struct Cosmology *Cx, double z);
double scale_indep_growth_f(struct Cosmology *Cx, double z);
double Hubble(struct Cosmology *Cx, double a);
double angular_distance(struct Cosmology *Cx, double z);
double comoving_radial_distance(struct Cosmology *Cx, double z);

double rhoc(struct Cosmology *Cx, double z);
double R_scale(struct Cosmology *Cx, double h_mass);

double R_vir(struct Cosmology *Cx, double h_mass);  
double c_cdm(double M, double z);
double nfw_profile(struct Cosmology *Cx, double k, double M, double z);

double window_kth(double k, double R);
double window_rth(double k, double R);
double derR_window_rth(double k, double R);
double window_g(double k, double R);
double derR_logwindow_g(double k, double R);


#endif