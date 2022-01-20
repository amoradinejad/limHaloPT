# limHaloPT

Author: Azadeh Moradinezhad Dizgah
<br>

with contributions from Farnik Nikakhtar 

Welcom to limHaloPT, a numerical package for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temprature fluctuations within halo-model framework. The current version of the code, is limited to real-space, and redshift-space distortions will be included in the next release. 

The extended halo model of line intensity power spectrum implemented in limHaloPT, combines the predictions of EFTofLSS for halo power spectrum with the standard halo model to account for the nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities (Infrared Resummation) in the 2-halo term. The deviations from Poisson shot noise on large scales are also computed within the halo model.

This package is released together with the following publication, [arXiv:2111.03717](https://arxiv.org/abs/2111.03717), where the prediction of the model are tested against new suite of simulated intensity (brightness temprature) maps of CO and [CII] lines. The mesheded fileds from MithraLIMSims will be publically avilable on [MithraLIMSims](http://cyril.astro.berkeley.edu/MithraLIMSims). The code to analyse the simulated maps is an extension of the toolskit used in analysin Hidden-Valley simulations, and is publically avialble on [LIM Analysis](https://github.com/farnikn/MithraLIMSims). As discussed in the manuscript above, the packages to compute the theory predictions and creating the simulated intensity maps can be straightforwardly extended to compute the power spectrum signal of other emission lines (emitted from star-froming galaxies), beside CO and [CII]. 
<br>

The full documentation of the code can be found on [Documentation](https://amoradinejad.github.io/limHaloPT/docs/html/index.html).

## Dependencies
The limHaloPT package calls various functions from [CLASS](https://github.com/lesgourg/class_public) Boltzman solver, including the matter power spectrum and transfer functions, growth factor etc. Therefore, you need to first download and compile CLASS code, and place the "libclass.a" file in the " CLASS/lib/" folder. Furtehrmore, the loop calculations are performed with direct numerical integration, using routines of [CUBA](http://www.feynarts.de/cuba/) library. Furthermore, the code heavily uses functions of [GSL](https://www.gnu.org/software/gsl/doc/html/) scientific library. Therfore, make sure that the two libraries are correctly linked to limHaloPT by making necassary modifcations to the makefile (placed in Source directory) of limHaloPT package. 
<br>

## Compilation and Usage
- To compile, type: make <br>
- To run, type: ./limHaloPT 

If you modified the code, you need to first do "make clean" before doing "make". Depending on what quantities you want to calculate, you can modify the main() function in main.c module (as marked in the code). As examples, I have included the calls to two functions to compute the clustering and shot noise contributions. 
<br>

## Attribution
You can use this package freely, provided that in your publication you cite the following paper
Moradinezhad, Nikakhtar, Keating, Castorina: [arXiv:2111.03717](https://arxiv.org/abs/2111.03717). Furthermore, since limHaloPT relies on CLASS Boltzman code, you should also cite at least this paper [arxiv:1104.2933](https://arxiv.org/abs/1104.2933) as required byy CLASS developers. 
<br>

## License
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Copyright 2021 Azadeh Moradinezhad Dizgah.<br>
limHaloPT is free software made available under the MIT License. For details see the [LICENSE](https://github.com/amoradinejad/limHaloPT/blob/d40a4a75188ae70f56ed76236d1fd9ee1aae312d/LICENSE) file.


