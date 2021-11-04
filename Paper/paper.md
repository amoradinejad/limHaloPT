---
title: 'limHaloPT: A numerical Package for Accurate Modeling of Line Intensity Power spectrum'
tags:
  - C language
  - cosmology
  - intensity mapping
authors:
  - name: Azadeh Moradinezhad Dizgah
    orcid: 0000-0001-8841-9989
    affiliation: 1
affiliations:
 - name: Departement of Theoretical Physics, University of Geneva, 24 quai Ernest Ansermet, 1211 Geneva 4, Switzerland
   index: 1
date: November 2021
bibliography: paper.bib
---

&nbsp;
&nbsp;
&nbsp;


# Summary

`limHaloPT` is a modular numerical package, written in C, for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temperature fluctuations using the halo-model framework. This package is the first publically available code, which combines the one-loop prediction of the halo power spectrum and the halo-model framework to model the power spectrum of emission lines originating from star-forming galaxies. Furthermore, the code includes routines to compute the stochastic contributions to line power spectrum beyond the Poisson approximation. Several utility functions, e.g., for computing the theoretical halo mass functions, halo biases, one-loop halo power spectrum, are provided in the package, which can be used in contexts other than LIM. This code is released together with a scientific publication [@MoradinezhadDizgah:2021], in which details of the implemented model and the comparison of model predictions against simulated intensity maps are presented. The current version of the code includes the rotational ladder of carbon monoxide, CO, and fine structure line of ionized carbon [CII] and is limited to real space.


# Scientific Context and Statement of Need

Line intensity mapping (LIM) is a novel technique to map the large-scale structure (LSS) of the Universe by measuring aggregate emission of the atomic and molecular emission lines from the unresolved source [@Kovetz:2017agg]. Measurements of spatial fluctuations and frequency of the line provide a 3-dimensional map of the LSS, whose statistical properties capture a significant amount of information about astrophysics and cosmology. To fully exploit this rich data, accurate theoretical models of the signal and efficient numerical codes for evaluating the models are crucial. 

The modeling of the line signal is based on the halo-model framework and requires two main ingredients; modeling the relation between line luminosity and halo masses and modeling the relation between halo properties and the underlying dark matter distribution. Until now, the models used in the literature neglect the nonlinear effects in the latter relation and use the tree-level perturbation theory to relate the halo properties, and by extension line intensity properties, to the underlying dark matter distribution. As of numerical implementation, the only publicly available code to compute the line power spectrum, [HaloGen](https://github.com/EmmanuelSchaan/HaloGen/tree/LIM) [@Schaan:2021hhy], is based on this simplified model.   

The extended halo model of line intensity power spectrum implemented in `limHaloPT` combines the predictions of EFTofLSS in Eulerian space for halo power spectrum [@Baumann:2010tm; @Carrasco:2012cv] with the standard halo model [@Seljak:2000gq; @Cooray:2002dia] to account for nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities i.e., the Infrared Resummation [@Senatore:2014via; @Senatore:2017pbn; @Blas:2016sfa] in the 2-halo term. The deviations from Poisson shot noise on large scales are also computed within the halo model (see, e.g., [@Ginzburg:2017mgf]).

Recently, there has been a shift in the cosmology community in publicly releasing the packages developed by various groups to facilitate the follow-up research by the wider community, without the need of each research group to replicate the numerical tools previously developed by other groups. Great examples of this approach are [CLASS](https://github.com/lesgourg/class_public) Boltzman code [@Blas2011], and [nbodykit](https://nbodykit.readthedocs.io/en/latest/) [@Hand:2017pqn] toolkit for analysis LSS data from Nbody simulations and from galaxy surveys. In LIM, limHaloPT is the first package that includes detailed modeling of the line power spectrum. The modular structure of the package facilitates future extensions of the code to other LIM statistics, such as bispectrum, as well embedding this code in a full likelihood analysis pipeline such as [CosmoSIS](https://bitbucket.org/joezuntz/cosmosis/wiki/Home) [@Zuntz:2014csq].   


# Dependencies

The `limHaloPT' package calls various functions from [CLASS](https://github.com/lesgourg/class_public) Boltzman solver, including the matter power spectrum and transfer functions, growth factor, etc. Furthermore, the loop calculations are performed with direct numerical integration, using routines of [CUBA](http://www.feynarts.de/cuba/) and [Cubature](https://github.com/stevengj/cubature) libraries for C. The code also uses several functions of [GSL](https://www.gnu.org/software/gsl/doc/html/) scientific library. 


# Future Extensions 

Extension of the code to include other emission lines originating from star-forming galaxies is straightforward. Future releases will provide additional modules, for example, to include observational effects such as redshift-space distortions and the Alcock-Paczynski effect. Furthermore, we plan to extend this code to higher-order statistics of line intensity fluctuations. 


# Acknowledgements

The research of A.M.D. is supported by the SNSF project "The  Non-Gaussian  Universe and  Cosmological Symmetries", project number:200020-178787. A.M.D., further acknowledges partial funding from Tomalla Foundation for Research in Gravity. 


# References

