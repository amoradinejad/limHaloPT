# limHaloPT

Author: Azadeh Moradinezhad Dizgah


Welcom to limHaloPT, a numerical package for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temprature fluctuations within halo-model framework. The current version of the code, is limited to real-space. Redshift-space distortions will be included in the next release. 

The extended halo model of line intensity power spectrum implemented in limHaloPT, combines the predictions of EFTofLSS for halo power spectrum with the standard halo model to account for the nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities (Infrared Resummation) in the 2-halo term.
The deviations from Poisson shot noise on large scales is computed within halo model (see e.g. Guinzberg et al arXiv:1706.08738 for previous work in the context of halo/galaxy clustering).

The loop calculations, as well as integrations over halo masses are performed using integration routines of CUBA library (see https://arxiv.org/abs/hep-ph/0404043). Therfore, you need to compile this library and place the lib file in an Include directory that is linked to limHaloPT.

This package is released together with the following publication, arxiv:2111.XXXXX, where the prediction of the model are tested against new suite of simulated intensity (brightness temprature) maps of CO and [CII] lines. The mesheded fileds from MithraLIMSims are publically avilable on http://cyril.astro.berkeley.edu/MithraLIMSims.  As discussed in the paper, this code can be straightforwardly extended to compute the power spectrum signal of other emission lines (emitted from star-froming galaxies), beside CO and [CII].


## Compilation and Usage

To compile the package, type:
- make

To run the code, type:
- ./limHaloPT  

Depending on what quantities you want to calculate, you can modify the main() function in main.c module (as marked in the code). As examples, I have included the calls to two functions to compute the clustering and shot noise contributions. 

If you modified the code, you need to first do "make clean" before doing "make".


## Acknowledgment

You can use this package freely, provided that in your publication you cite the following paper: 
Moradinezhad & Nikakhtar & Keating & Castorina: arXiv:2111.XXX


## License

MIT Liscence

