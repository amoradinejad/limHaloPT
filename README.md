# limHaloPT

Author: Azadeh Moradinezhad Dizgah
<br>

with contributions from Alberto Vallinotto and Farnik Nikakhtar 

Welcom to limHaloPT, a numerical package for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temprature fluctuations within halo-model framework. The current version of the code, is limited to real-space, and redshift-space distortions will be included in the next release. 

The extended halo model of line intensity power spectrum implemented in limHaloPT, combines the predictions of EFTofLSS for halo power spectrum with the standard halo model to account for the nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities (Infrared Resummation) in the 2-halo term. The deviations from Poisson shot noise on large scales are also computed within the halo model.

This package is released together with the following publication, [arXiv:2111.03717](https://arxiv.org/abs/2111.03717), where the prediction of the model are tested against new suite of simulated intensity (brightness temprature) maps of CO and [CII] lines. The mesheded fileds from MithraLIMSims will be publically avilable on [MithraLIMSims](http://cyril.astro.berkeley.edu/MithraLIMSims). The code to analyse the simulated maps is an extension of the toolkit used in analysing Hidden-Valley simulations, and is publically avilable on [LIM Analysis](https://github.com/farnikn/MithraLIMSims). As discussed in the manuscript above, the packages to compute the theory predictions and creating the simulated intensity maps can be straightforwardly extended to compute the power spectrum signal of other emission lines (emitted from star-froming galaxies), beside CO and [CII]. 

The full documentation of the code can be found on [Documentation](https://amoradinejad.github.io/limHaloPT/html/index.html). 
<br>


## Dependencies
The limHaloPT package calls various functions from [CLASS](https://github.com/lesgourg/class_public) Boltzman solver, including the matter power spectrum and transfer functions, growth factor etc. Therefore, you need to first download and compile CLASS code, create a " CLASS/lib/" folder and place the "libclass.a" in that folder. Furtehrmore, the loop calculations are performed with direct numerical integration, using routines of [CUBA](http://www.feynarts.de/cuba/) library. Lastly, the code heavily uses functions of [GSL](https://www.gnu.org/software/gsl/doc/html/) scientific library. Therfore, make sure that the two libraries are correctly linked to limHaloPT by making necassary modifcations to the makefile (placed in Source directory) of limHaloPT package. 

Note that in order to compute the luminosities of spectral lines, a model of star formation rate as a function of halo mass and redshift, SFR(M_h,z), should be assumed. Currently, the implemented model uses SFR(M_h,z) from [Behroozi et al. (2013)](https://arxiv.org/abs/1207.6105). The necassary input file, "sfr_release.dat", is included in "Input/release-sfh_z0_z8_052913/sfr/" subdirectory. <br>


## Structure of the package
limHaloPT consists of 10 main modules, which include the following categories of functions:
- main.c: This is the most external module, from which any function that you need is called. Depending on what quantities you want to calculate, you can modify the main() function in main.c module (as marked in the code). Two example calls to functions which compute the clustering and shot noise contributions is included in main.c module. After adding the function calls in this module, you need to re-compile the code and then run it. Further details of the modules and descriptions of the functions can be found in documentation of the code.
- ps_line_hm.c: This module includes various functions that are needed for computation of halo-model line power spectrum, including clustering and stochastic contributions beyond Poisson limit. At the moment, the redshift-space distrotions is not included in the halo model implementation. 
- ps_line_pt.c: This module includes various functions to compute the Poisson shot-noise and linear clustering components of line power spectrum, contributions of line interloper to the power spectrum signal. Both real and redshift-space linear power spectra are available.
- line_ingredients.c: This module includes all the necassary functions to compute the line power spectrum by the above two modules. Many of the functions in this module can be used to develope the extensions of limHaloPT to other line statistics beyond the power spectrum.
- ps_halo_1loop.c: This module includes the ingredients to compute the real-space 1loop contributions of the halo/galaxy power spectrum. Both loops due to gravitational evolution and primordial non-Gaussianity can be calculated. The computations of the loop integrals is performed with direct numerical integration routines from CUBA library.
- ir_res.c: This module includes functions to compute IR-resummed matter power spectrum at leading and next-to-leading order. 
- wnw_split.c: This module includes function to perform the splitting of matter power spectrum into its broadband and wiggle (BAO) contributions. 
- survey_specs.c: This module includes functions to calculate quantities related to the choices of maximum and minnimum scales that can be probed by a given survey.
- cosmology.c: This module includes function to compute the cosmological quantities at the level of background and perturbations. 
- read_input.c: This module contains two functions initialize() and clean() which are the first and last functions called by the main() function. Initialize() reads in the input values from an .ini files.
- setup_teardown.c: This module contains utilityy functions to read in an .ini file.
- utilities.c: This module contains some utility functions, for example to build a dynamically allocated 1-dimensional array. These utility functions are used by the rest of the modules. <br>


## Compilation 
- To compile, within the main directory of limHaloPT, type: "make" <br>
This would create an executable called "limHaloPT" in the same directory, which you will use to run the code. If you modified the code, you need to first do "make clean" before doing "make". 

The entire limHaloPT package was developed, compiled, and tested on Mac OS X, using gcc version 7.5.0 compiler. <br>


## Basic usage
- To run the code, in the main directory of limHaloPT, type  "./limHaloPT LCDM.ini"  

Currently, the main output of limHaloPT are the mean brightness temprature, linear and quadratic biases, clustering and shot noise contributions of 7 emission spectral lines, depending on the switch that you set in the ini file. The computation of these functions are performed within main.c module. An example of the ini file is provided (LCDM.ini). If you want to call any of the functions of limHaloPT, apart from those called by defacult, you neeed to add the function call to main.c, recompile the code by first cleaning the previous build using "make clean" and building the package again with "make". 

Example: Lets say you want to compute the mean brightness temprature, linear and quadratic biases for one or multiple emission lines, as a function of wavenumber at several redshifts. In the .ini file, you should set the switch "switch_Tbar" to be equal to "Tbar", and the switches of shot and power spectra blank, i.e. "switch_HMshot = " and "switch_HMshot = ". The relavent functions that are called within main.c are PS_line_HM() and line_bias(). For this function call, first two arrays for values of redshifts and wavenumbers for which the power spectrum is computed are created. the computed values for the three quantities is saved by default in an output file saved in "Output" directory. The number of emission lines to be included in the computation is set with "nlines" variable in the LCDM.ini file, while the name of the line is passed with line1, line2, ... variables. T The first argument of PS_line_HM, is the cosmology structure needed for performing any calculation in the code. See cosmology.c module for more details. The LCDM.ini file contains description of parmaeters that should be set in the .ini file.

So here is how PS_line_HM() function is called within main():
```
      int nz_mean     = gb.mean_nz;
      double *z_mean  = init_1Darray(nz_mean,0.0,gb.mean_zmax);

      double lbias_arr[2];
      double b1   = 0.;
      double b2   = 0.;
      double Tbar = 0.;

      FILE *fp1;
      char filename1[FILENAME_MAX];

      for(int i=0;i<nlines;i++){
            sprintf(filename1,"%s/Tbar_J%d.txt", gb.output_dir, JJ[i]);
            fp1    = fopen(filename1, "ab");
            fprintf(fp1,"%s %s %s %s %s \n","#", "z", "b1", "b2", "Tbar");
            for(int j=0;j<nz_mean;j++){
                  Tbar = Tbar_line(&Cx_ref, i, z_mean[j]);
                  line_bias(Cx_ref.Lines[i],z_mean[j],lbias_arr);
                  b1   = lbias_arr[0];
                  b2   = lbias_arr[1];
                  fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e \n",z_mean[j], b1, b2, Tbar); 
            }
            fclose(fp1);
      }
      free(z_mean);
```

By default, in additionn to having th eoutput in main() function for mean brightness temprature, biases, shot and clustering components of power spectrum, when computing the clustering and shot powers, individual loop contributions (in the former) and individual beyon-Possion contribution to the shot (in the latter), are also saved to output files which are stored in "Output" directory. <br>


<br>

## Attribution
You can use this package freely, provided that in your publication you cite the following paper
Moradinezhad, Nikakhtar, Keating, Castorina: [arXiv:2111.03717](https://arxiv.org/abs/2111.03717). Furthermore, since limHaloPT relies on CLASS Boltzman code, you should also cite at least this paper [arxiv:1104.2933](https://arxiv.org/abs/1104.2933) as required byy CLASS developers.<br>


## License
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Copyright 2021 Azadeh Moradinezhad Dizgah.<br>
limHaloPT is free software made available under the MIT License. For details see the [LICENSE](https://github.com/amoradinejad/limHaloPT/blob/d40a4a75188ae70f56ed76236d1fd9ee1aae312d/LICENSE) file.


## Contributing
To contribute, create a fork on github, make changes and commits, and submit a pull request on github. If you found any issues with the code, please get in touch by posting an issue on this github repository.