# limHaloPT

Author: Azadeh Moradinezhad Dizgah


Welcom to limHaloPT, a numerical package for computing the power spectrum of line intensity fluctuations using a one-loop halo model, and the large-scale shot noise power spectrum byond Poisson approximations. 

This package is released together with the following publication, arxiv:2111.XXXXX, where the prediction of the model are tested against new suite of simulated intensity (brightness temprature) maps of CO and [CII] lines. The mesheded fileds from MithraLIMSims are publically avilable on http://cyril.astro.berkeley.edu/MithraLIMSims. 

The extended halo model of line intensity power spectrum implemented in MithraLIMSims, combines the predictions of EFTofLSS for halo power spectrum with the standard halo model to account for the nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter disctribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities (Infrared Resummation) in the 2-halo term.
The deviations from Poisson shot noise on large scales is computed within halo model (see e.g. Guinzberg et al arXiv:1706.08738 for previous work in the context of halo/galaxy clustering).

The loop calculations, as well as integrations over halo masses are performed using CUBA library. Therfore, you need to compile this library and place the lib file in an Include directory that is linked to limHaloPT.


## Compilation and Usage

To compile the package, type:
- make

To run the code, type:
- ./run.exe  



## Acknowledgment

If using this package in a publication, please cite the following paper: 
Moradinezhad & Nikakhtar & Keating & Castorina: arXiv:2111.XXX


## License

MIT Liscence

