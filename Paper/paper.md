---
title: 'limHaloPT: A Numerical Package for Accurate Modeling of Line Intensity Power Spectrum'
tags:
  - C language
  - cosmology
  - intensity mapping
authors:
  - name: Azadeh Moradinezhad Dizgah
    orcid: 0000-0001-8841-9989
    affiliation: 1
  - name: Alberto Vallinotto  
    affiliation: 2
  - name: Farnik Nikakhtar
    orcid: 0000-0002-3641-4366
    affiliation: 3
affiliations:
 - name: Departement of Theoretical Physics, University of Geneva, 24 quai Ernest Ansermet, 1211 Geneva 4, Switzerland
   index: 1
 - name: Independent Researcher, Italy
   index: 2  
 - name: Department of Physics and Astronomy, University of Pennsylvania, 209 S. 33rd St., Philadelphia, PA 19104, USA
   index: 3   
date: 10 November 2021
bibliography: paper.bib
---

&nbsp;
&nbsp;
&nbsp;



# Summary

`limHaloPT` is a modular numerical package, written in C, for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temperature fluctuations using the halo-model framework. This package is the first publicly available code, which combines the one-loop prediction of the halo power spectrum and the halo-model framework to model the power spectrum of emission lines originating from star-forming galaxies. Furthermore, the code includes routines to compute the stochastic contributions to the line power spectrum beyond the Poisson approximation. Several utility functions, e.g., for computing the theoretical halo mass functions, halo biases, one-loop halo power spectrum, are provided in the package, which can be used in contexts other than Line intensity mapping (LIM). This code is released together with a scientific publication [@MoradinezhadDizgah:2021dei], in which details of the implemented model and the comparison of model predictions against simulated intensity maps are presented. The current version of the code includes the first six rotational ladder of carbon monoxide, CO, and fine structure line of ionized carbon [CII], and the model of the power spectrum is limited to real space.


# Scientific Context and Statement of Need

Line intensity mapping is a novel technique to map the large-scale structure (LSS) of the Universe by measuring aggregate emission of the atomic and molecular emission lines from the unresolved sources [@Kovetz:2017agg]. Measurements of spatial fluctuations and frequency of the line provide a 3-dimensional map of the LSS, whose statistical properties capture a significant amount of information about astrophysics and cosmology. To fully exploit this rich data, accurate theoretical models of the signal and efficient numerical codes for evaluating the models are crucial. 

The modeling of the line signal is based on the halo-model framework and requires two main ingredients; modeling the relation between line luminosity and halo masses and modeling the relation between halo properties and the underlying dark matter distribution. Until now, the models used in the literature neglect the nonlinear effects in the latter relation and use the tree-level perturbation theory to relate the halo properties, and by extension line intensity properties, to the underlying dark matter distribution. Furthermore, the shot noise component of the line power spectrum is commonly addumed to be Poissonian. As of numerical implementation, the only publicly available code to compute the line power spectrum, [HaloGen](https://github.com/EmmanuelSchaan/HaloGen/tree/LIM) [@Schaan:2021hhy], is based on this simplified model.   

The extended halo model of line intensity power spectrum implemented in `limHaloPT` combines the predictions of EFTofLSS in Eulerian space for halo power spectrum [@Baumann:2010tm; @Carrasco:2012cv] with the standard halo model [@Seljak:2000gq; @Cooray:2002dia] to account for the nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities i.e., the Infrared Resummation [@Senatore:2014via; @Blas:2016sfa] in the 2-halo term. The deviations from Poisson shot noise on large scales are also computed within the halo model [@Ginzburg:2017mgf].

Recently, there has been a shift in the cosmology community in publicly releasing the packages developed by various groups to facilitate the follow-up research by the wider community, without the need of each research group to replicate the numerical tools previously developed by other groups. Great examples of this approach are the [CLASS](https://github.com/lesgourg/class_public) Boltzman code [@Blas2011], and the [nbodykit](https://nbodykit.readthedocs.io/en/latest/) toolkit for analysis of the LSS data from Nbody simulations and from galaxy surveys [@Hand:2017pqn]. In LIM, limHaloPT is the first package that includes detailed modeling of the line power spectrum. The modular structure of the package facilitates future extensions of the code to other LIM statistics, such as the line bispectrum, as well as embedding this code in a full likelihood analysis pipeline such as [CosmoSIS](https://bitbucket.org/joezuntz/cosmosis/wiki/Home) [@Zuntz:2014csq].   


# Dependencies

The `limHaloPT` package calls various functions from the [CLASS](https://github.com/lesgourg/class_public) Boltzmann solver, including the matter power spectrum, the transfer functions, the growth factor, etc. Therefore, prior to installation of `limHaloPT`, the CLASS code should be compiled and the static library of "libclass.a" should be placed in the "Class/lib" folder. Furthermore, the loop calculations are performed with direct numerical integration, using routines of [CUBA](http://www.feynarts.de/cuba/) library. Lastly, the code heavily uses functions of the [GSL](https://www.gnu.org/software/gsl/doc/html/) scientific library. Therefore, the two libraries should be correctly linked to `limHaloPT` by making necessary modifications to the makefile (placed in the main directory). 

In order to compute the luminosities of spectral lines, a model for star formation rate as a function of halo mass and redshift, SFR(M_h,z), should be assumed. Currently, the implemented model uses SFR(M_h,z) from [Behroozi et al. (2013)](https://arxiv.org/abs/1207.6105). The necessary input file, "sfr_release.dat", is included in "Input/release-sfh_z0_z8_052913/sfr/" subdirectory. 


# Usage

Currently, the main output of `limHaloPT` are the mean brightness temprature, the linear and the quadratic biases, the clustering and the shot noise contributions of seven emission spectral lines. The lines to compute are set using a switch in the ini file, an example of the ini file is provided in the package (see "LCDM.ini"). More detailed description of the main outputs is given in the code documentation on github repository. In addition to saving the output for mean brightness temprature, biases, shot and clustering components of power spectrum, when computing the clustering and shot powers, individual loop contributions (in the former) and individual beyond-Possion contribution to the shot (in the latter), are also saved to output files which are stored in "Output" directory, by default. 

In addition to using `limHaloPT` through "main.c" module which calls three specific functions, the package can also be exported as an external library to any other C code. Upon compilation of the package, a static library called "liblimHaloPT.a" is created and placed in the "lib" directory, which can then be linked to an external C code. Example of how to link an external code to "liblimHaloPT.a" static library can be found in the "Test" directory.



# Future Extensions 

Future releases will provide additional modules, for example, to include observational effects such as redshift-space distortions and the Alcock-Paczynski effect. Furthermore, we plan to extend this code to include modeling of other emission lines originating from star-forming galaxies, cross-correlations between different emission lines, and the bispectrum of line intensity fluctuations. 


# Acknowledgements

The research of A.M.D. is supported by the SNSF project "The  Non-Gaussian  Universe and  Cosmological Symmetries", project number:200020-178787. A.M.D., further acknowledges partial funding from Tomalla Foundation for Research in Gravity. 


# References

