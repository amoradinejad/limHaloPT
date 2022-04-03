# limHaloPT

Author: Azadeh Moradinezhad Dizgah
<br>

with contributions from Farnik Nikakhtar 

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
- utilities.c: This module contains some utility functions, for example to build a dynamically allocated 1-dimensional array. These utility functions are used by the rest of the modules. <br>


## Compilation 
- To compile, type: make <br>

If you modified the code, you need to first do "make clean" before doing "make". The entire limHaloPT package was developed, compiled, and tested on Mac OS X, using gcc version 7.5.0 compiler. <br>


## Basic usage
In the current version of limHaloPT, to use the package, you need to modify the main.c module. After each modification, the package needs to be re-compiled before it can be run. In the future versions, you should be able to just set an .ini file to set which quantities to be computed and the values of the parmaeters to be used. 

Depending on what quantities you want to calculate, you should modify the main() function in main.c module (as marked in the code). Before calling any function within main.c, you may want to also change default values for some other initialization steps. This is also marked in the code. Two example calls to functions that compute the clustering and shot noise contributions aree included in main.c module. We descibe one of them here. 

Example: Lets say you want to compute the power spectrum of mean brightness temprature fluctuations for one or multiple emission lines, within halo-model, as a function of wavenumber, and at several redshifts. The relavent function that needs to be called is PS_line_HM(). First you should create two arrays for values of redshifts and wavenumbers for which you want to calculate the power spectrum. If you only want to calculate the power spectrum for a single wavenumber and redshift, you would pass just two numbers to PS_line_HM(). You should also define for which emission lines do you want to compute the power spectrum. The choice of the lines is done in the first part of the main.c module, where the number of the lines, their name and their corresponding index is set. For example to compute the power spectrum for CO(1-0) and [CII] you would set 
```
int nlines     = 2;
int lines[2]   = {CO10,CII};
int JJ[2]      = {1,0};
``` 

You would refer to each line in the example code below by its index inside the loop over i which extends to i=nlines. The name of the lines and their id is passed to PS_line_HM() function as the last two arguments. The first argument of PS_line_HM, is the cosmology structure needed for performing any calculation in the code. See cosmology.c module for more details. 

So here is how you would call the PS_line_HM() function:
```
      int nk     = 200;
      int nz     = 50;
      double *k  = loginit_1Darray(nk, 1.e-3, 2.);
      double *z  = init_1Darray(nz,0.0,11.);

      /** 
       * Compute the line clustering signal using halo model 
       */
      printf("Calculating halo-model line power spectrum\n");
      int line_id = 0;
      double ps_clust_hm = 0.;
      for(int i=0;i<nlines;i++){
            line_id = i;
            for(int j=0;j<nz;j++){
                  for(int l=0;l<nk;l++){
                        ps_clust_hm =PS_line_HM(&Cx_ref, k[l], z[j], M_min, mode_mf, lines[i], line_id);
                        printf("%d %12.6e %12.6e %12.6e \n", i, z[j], k[l], ps_clust_hm);

                  }
            }     
      }
```

Once you made your modifications to main.c module, you should do the following:
- Clean the previous make by "make clean"
- Re-compile the code with  "make" 
- Run the code by  "./limHaloPT" <br>


## Attribution
You can use this package freely, provided that in your publication you cite the following paper
Moradinezhad, Nikakhtar, Keating, Castorina: [arXiv:2111.03717](https://arxiv.org/abs/2111.03717). Furthermore, since limHaloPT relies on CLASS Boltzman code, you should also cite at least this paper [arxiv:1104.2933](https://arxiv.org/abs/1104.2933) as required byy CLASS developers.<br>


## License
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Copyright 2021 Azadeh Moradinezhad Dizgah.<br>
limHaloPT is free software made available under the MIT License. For details see the [LICENSE](https://github.com/amoradinejad/limHaloPT/blob/d40a4a75188ae70f56ed76236d1fd9ee1aae312d/LICENSE) file.


## Contributing
To contribute, create a fork on github, make changes and commits, and submit a pull request on github. If you found any issues with the code, please get in touch by posting an issue on this github repository.
