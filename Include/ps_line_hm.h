#ifndef _LINE_PS_HM_H_
#define _LINE_PS_HM_H_

// double PS_line_HM_RSD(struct Cosmology *Cx, double k, double mu, double z, size_t line_id, struct Cosmology *Cx_ref);
double PS_line_HM(struct Cosmology *Cx, double k, double z, double M_min, long mode_mf,  long line_type, int line_id);
double PS_shot_HM(struct Cosmology *Cx, double k, double z, double M_min, double *input, long mode_mf, long line_type);


// void HM_1h2h_interp(double k, double z, size_t line_id, struct Cosmology *Cx, long mode_hm, double *result);

void HM_1h2h(struct Cosmology *Cx, double k, double z, double M_min,
            long mode_mf,long line_type, long mode_hm, double *result);
static int HM_1h2h_integ(const int *ndim,
                  const cubareal x[],
                  const int *ncomp,
                  cubareal ff[],
                  void *p);
void mhmc(struct Cosmology *Cx, double z, long mode_mf, double *result);
static int mhmc_integ(const int *ndim,
                  const cubareal x[],
                  const int *ncomp,
                  cubareal ff[],
                  void *p);


// double b22_ls_integrand(double x, void *par);
double b22_ls(struct Cosmology *Cx, double z);
static int b22_ls_integrand(const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p); 

#endif