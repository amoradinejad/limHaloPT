
#ifndef _LINE_PS_PT_H_
#define _LINE_PS_PT_H_

double PS_hh_G(struct Cosmology *Cx, double k, double z, double M, long mode, long IR_switch, long SPLIT, long mode_mf);
double PS_hh_PNG(struct Cosmology *Cx, double k, double z, double M, long mode, long IR_switch, long SPLIT, long mode_mf);

void Compute_G_loops(struct Cosmology *Cx, double k, double z, long IR_switch, long hm_switch, long SPLIT, double *result);
void Compute_PNG_loops(struct Cosmology *Cx, double k, double z, long IR_switch, long SPLIT, double *result);
static int G_loop_integrands(const int *ndim,
                              const cubareal x[],
                              const int *ncomp,
                              cubareal ff[],
                              void *p);
static int PNG_loop_integrands(const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p);

double F2_s(double k1,double k2,double mu);
double S2_s(double k1,double k2,double mu);
double F3_s(double k,double q, double mu);
double S2(double mu);
double F2(double k1,double k2,double mu);


#endif