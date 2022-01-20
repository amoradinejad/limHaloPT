#ifndef LINE_PS_H
#define LINE_PS_H

double PS_line( struct Cosmology *Cx, double k, double z, size_t line_id);
double PS_line_RSD(struct Cosmology *Cx, struct Cosmology *Cx_ref,double k, double mu, double z, size_t line_id);
double PS_shot(struct Cosmology *Cx, double z, size_t line_id);

double ps_line_multipoles(struct Cosmology *Cx, struct Cosmology *Cx_ref, double k, double z, size_t line_id, int ell);
int ps_line_multipoles_integrand(unsigned    ndim,           // Number of dimensions in the domain space -- number of dim we're integrating over
                        const double *x,              // The point at which the integrand is evaluated
                        void         *p,              // Pointer to a structure that holds the parameters
                        unsigned     fdim,            // Number of dimensions that the integrand return
                        double       *fvalue          // Array of values of the integrand of dimension fdim
                        );

#endif