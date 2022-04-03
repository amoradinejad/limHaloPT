
/** @file setup_teardown.c Documented setup_teardown module, including the initialize and clean functions that are called first and last thing in limHaloPT package
 *
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 */

#include "header.h"


/**
 * Initizlize the path to the required directories, set the values of cosmological parmaeters, and initialize the interpolator of the SFR(M,z) 
 * from tabulated data provided in gb.SFR_filename.  
 * 
 * The global structure "gb" have several elements to hold the paths to project source directory, input, and output folders, 
 * and values of cosmological parmaeters. 
 * 
 * @return void 
 */
void initialize()
{
      getcwd(gb.project_home , sizeof(gb.project_home));  //This gives the path to the source directory      
      sprintf(gb.output_dir,"%s/Output", gb.project_home);
      sprintf(gb.data_dir,"%s/Input", gb.project_home);
      sprintf(gb.SFR_filename,"%s/release-sfh_z0_z8_052913/sfr/sfr_release.dat",gb.data_dir);  

      logSFR_alloc_init();

      /**
       * Cosmological parameters corresponding to initial conditions of HiddenValley (HV) simulations
       */
      gb.c         = 2.99792458e5;  /// In units of km/s
      gb.h         = 0.677;
      gb.Omega_cdm = 0.11923/pow(gb.h,2.);  ////omega_cdm = Omega_cdm h^2 ;
      gb.Omega_b   = 0.02247/pow(gb.h,2.);   ///omega_b = Omega_b h^2;
      gb.Omega_r   = 0.0000910958;  ////radiation = photons + neutrinos
      gb.Omega_g   = 5.3956421715871286e-05;   ////photons, input for Class
      gb.Omega_nu  = 0.00;    ////neutrinos
      gb.ns        = 0.96824;
      gb.As        = 2.1085e-9;
      gb.logAs     = log(gb.As*pow(10.,10.));  ///3.0665
      gb.kp        = 0.05;  //in unit of Mpc^-1

      gb.sigFOG0 = 250.;      

      gb.PS_kmin = 1.e-4;
      gb.PS_kmax = 200.;
      // gb.kmax_CLASS = PS_KMAX * (1.-0.1);
      gb.kmax_CLASS = PS_KMAX;
      gb.Npars = NPARS;

      gb.PS_zmax   = 14.;
      gb.line_zmax = 14.;

      return;
}


/**
 * Free the memory allocated to cosmology structure and SFR interpolator
 * 
 * @return void 
 */
void cleanup(struct Cosmology *Cx)
{
      Cosmology_free(Cx);
      SFR_Behroozi_free();

      return;
}



