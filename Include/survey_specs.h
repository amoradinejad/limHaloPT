#ifndef SURVEY_SPECS_H
#define SURVEY_SPECS_H

double shell_volume( struct Cosmology *Cx,double z, double fsky);
double kmin_val(struct Cosmology *Cx, double zmin, double zmax, double fsky);
double kmin(struct Cosmology *Cx, double zmin, double zmax);

double kmax_Brent(double kmax, void *params);
double kmax_Brent_solver( struct Cosmology *Cx, double z);
      

#endif