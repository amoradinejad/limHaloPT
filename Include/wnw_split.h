
#ifndef _WNW_SPLIT_H_
#define _WNW_SPLIT_H_

double pk_nw_integrand(double x, void *par); 
double pk_Gfilter_nw( struct Cosmology *Cx, double k, double kf0); 
double EH_PS_w( struct Cosmology *Cx, double k, double k0, double p0);
double EH_PS_nw( struct Cosmology *Cx, double k,double k0,double p0);
double T0( struct Cosmology *Cx, double k);
double T( struct Cosmology *Cx, double k);
double Tt0( struct Cosmology *Cx, double k, double x1, double x2);

#endif