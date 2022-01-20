
#ifndef IR_RES_H
#define IR_RES_H

double pm_nowiggle_dst(struct Cosmology *Cx, double k,double z, int mode);
double pm_nowiggle_gfilter(struct Cosmology *Cx, double k, double z, int mode);
double pm_nowiggle_bspline(struct Cosmology *Cx, double k, double z, int mode);

double pm_nowiggle(struct Cosmology *Cx, double k, double z, double kf0, int cleanup, long SPLIT);
double pm_IR_LO(struct Cosmology *Cx, double k, double z,  long SPLIT);
double pm_IR_NLO(struct Cosmology *Cx, double k, double z, long SPLIT);

double IR_Sigma2_integrand(double x, void *par);
double IR_Sigma2(struct Cosmology *Cx, double z,double kf0, long SPLIT);

#endif