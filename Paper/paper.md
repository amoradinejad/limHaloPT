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
 - name: Departement of Theoretical Physics, Universityy of Geneva, 24 quai Ernest Ansermet, 1211 Geneva 4, Switzerland
   index: 1
date: November 2021
bibliography: paper.bib
---

# Summary

Line intensity mapping (LIM) is a novel technique to map the large-scale structure (LSS) of the Universe over a wide range of redshifts and scales, largely inaccessible to traditional galaxy surveys. Instead of resolving individual galaxies, LIM surveys measure aggregate emission of the atomic and molecular emission lines from the unresolved source. Measurements of spatial fluctuations and frequency of the line provide a 3-dimensional map of the LSS, which can be used to constrain astrophysics and cosmology. In particular, these measurements provide a detailed view of the underlying dark matter distribution, which carries a wealth of information about the origin, composition, and evolution of the Universe. As such, these observations offer a tremendous opportunity to test fundamental physics. To fully exploit this rich data, accurate theoretical models of the signal are crucial, which match the precision of the data. 


# Statement of need

limHaloPT, a numerical package for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temprature fluctuations within halo-model framework. The current version of the code, is limited to real-space, and redshift-space distortions will be included in the next release.

The extended halo model of line intensity power spectrum implemented in limHaloPT, combines the predictions of EFTofLSS for halo power spectrum with the standard halo model to account for the nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities (Infrared Resummation) in the 2-halo term. The deviations from Poisson shot noise on large scales is computed within the halo model.




# Acknowledgements

The author of this code is supported by the SNSF project "The  Non-Gaussian  Universe  and  Cosmological Symmetries", project number:200020-178787. 


# References
