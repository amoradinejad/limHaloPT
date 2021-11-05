#ifndef LINE_INGREDIENT_H
#define LINE_INGREDIENT_H

double mass_func_sims(struct Cosmology *Cx, double M, double z, long mode_mf);
struct Line * Line_alloc_init(struct Cosmology *Cx, long line_type, size_t npoints_interp, double M_min, long mode_mf);
int Line_free(struct Line * Lx);
int Line_evaluate(struct Line * Lx, double *zz, double *res);

double mult_func(double sigma, long mode_mf);
double mass_func(struct Cosmology *Cx, double M, double z,long mode_mf);
void halo_bias(struct Cosmology *Cx, double M, double z, long mode_mf, double *bias_arr);

void   logSFR_Behroozi_read(double *z_arr, double *logM_arr, double *log10SFR);
int    logSFR_alloc_init();
double logSFR_Behroozi(double logM, double z);
int    SFR_Behroozi_free();
double luminosity(double M, double z, long mode_lum);

// static int mass_moment1_integ(const int *ndim,
//            const cubareal x[],
//            const int *ncomp,
//            cubareal ff[],
//            void *p) ;
// static int mass_moment2_integ(const int *ndim,
//            const cubareal x[],
//            const int *ncomp,
//            cubareal ff[],
//            void *p) ;
// static int bias_lum_weighted_integ(const int *ndim,
//            const cubareal x[],
//            const int *ncomp,
//            cubareal ff[],
//            void *p);
int mass_moment1_integ( unsigned     ndim,      // Number of dimensions in the domain space -- number of dim we're integrating over
        const double   *x,          // The point at which the integrand is evaluated
        void       *p,          // Pointer to a structure that holds the parameters
        unsigned     fdim,          // Number of dimensions that the integrand return
        double     *fvalue          // Array of values of the integrand of dimension fdim
        );
int mass_moment2_integ( unsigned     ndim,      // Number of dimensions in the domain space -- number of dim we're integrating over
        const double   *x,          // The point at which the integrand is evaluated
        void       *p,          // Pointer to a structure that holds the parameters
        unsigned     fdim,          // Number of dimensions that the integrand return
        double     *fvalue          // Array of values of the integrand of dimension fdim
        );

int bias_lum_weighted_integ(unsigned  nd,               // Number of dimensions in the domain space -- number of dim we're integrating over
                        const double  *x,                   // The point at which the integrand is evaluated
                        void          *p,               // Pointer to a structure that holds the parameters
                        unsigned      fdim,         // Number of dimensions that the integrand return
                        double        *fvalue           // Array of values of the integrand of dimension fdim
                        );
double mass_moment1(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum); 
double mass_moment2(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum) ;
void   bias_lum_weighted(struct Cosmology *Cx, double z, double M_min, long mode_mf, long mode_lum, double *result);  

void line_bias(struct Line *Lx, double z, double *result);
double mean_intens(struct Cosmology *Cx, size_t line_id, double z);
double Tbar_line(struct Cosmology *Cx, size_t line_id, double z);

double p_sig_shot_integrand(double x, void *par);
double p_sig_Tbar_integrand(double x, void *par);
double p_sig_shot(double sig_line);
double p_sig_Tbar(double sig_line);

#endif