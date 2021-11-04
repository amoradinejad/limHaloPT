---
title: 'limHaloPT: A numerical Package for Accurate Modeling of Line Intensity Power spectrum: one-loop halo model and non-Possionian shot'
tags:
  - C language
  - cosmology
  - intensity mapping
  - one-loop halo model
authors:
  - name: Azadeh Moradinezhad Dizgah ^[corresponding author]
    orcid: 0000-0001-8841-9989
    affiliation: 1
affiliations:
 - name: Departement of Theoretical Physics, University of Geneva, 24 quai Ernest Ansermet, 1211 Geneva 4, Switzerland
   index: 1
date: November 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---


# Summary

limHaloPT, is a modular numerical package, written in C language, for computing the clustering and shot-noise contributions to the power spectrum of line intensity/temprature fluctuations (emitted from within star-forming galaxies) using halo-model framework. This is the first publically available code, which combines the one-loop prediction of the halo power spectrum and the halo-model framework to model the line power spectrum. Furthermore, the code includes routines to compute the stochastic contributions to line power spectrum, beyond Poisson approximation. Several utility functions, e.g., to for computing the theoretical halo mass functions, halo biases, one-loop halo power spectrum, are provided in the package, which can be used in the contexts other than LIM. This code is realsed together with an accompanying scientific publications [@MoradinezhadDizgah:2021], in which detailes of the implemented model, as well as comparison of model predictions against simulated intensityy maps are presented. The current version of the code, is limited to real-space. Future releases will provide additional modules, for example to include observational effects such as redshift-space distortions and Alcock-Paczynski effect. Furthermore, we plan to extend this code to higher-order statistics of LIM. 


# Scientific context and Statement of need

Line intensity mapping (LIM) is a novel technique to map the large-scale structure (LSS) of the Universe by measuring aggregate emission of the atomic and molecular emission lines from the unresolved source [@Kovetz:2017agg]. Measurements of spatial fluctuations and frequency of the line provide a 3-dimensional map of the LSS, whose statistical properties capture significant amount of information about astrophysics and cosmology. To fully exploit this rich data, accurate theoretical models of the signal and efficient numerical codes for evaluating the models are crucial. The modeling of the line signal is based on halo-model framework, and require two assumptions; modeling of the relation between line luminosity and halo masses, and modeling the relation between halo properties and the underlying dark matter distribution. Until now, the models used in the literature neglect the nonlinear effects in the latter relation and are use the tree-level perturbbation theory to relate the halo properties, and by extension line intensity properties, to the underlying dark matter distribution. As of numerical implementation, the only publicaly avaialable code to compute the line power spectrum, HaloGen [@Schaan:2021hhy], is based on this simplified model.   

The extended halo model of line intensity power spectrum implemented in limHaloPT, combines the predictions of EFTofLSS (in Eulerian space, and assuming renormalized halo bias expansion) for halo power spectrum [@Baumann:2010tm; @Carrasco:2012cv] with the standard halo model [@Seljak:2000gq; @Cooray:2002dia] to account for nonlinear evolution of matter fluctuations and the nonlinear biasing relation between line intensity fluctuations and the underlying dark matter distribution in 2-halo term. Furthermore, the model includes the effect of large bulk velocities (Infrared Resummation) in the 2-halo term. The deviations from Poisson shot noise on large scales are also computed within the halo model (see e.g., [@Ginzburg:2017mgf]).

Recently, there has been a shift in the community in publicly releasing the packages developed by varous groups to facilitate the follow up research by the wider community, without the need of each research group to replicate the numerical tools previously developed by other groups. Great examples of this approch are CLASS Boltzman code and nbodykit toolkit for analysis LSS data from Nbody simulations and data from galaxy surveys. In the context of LIM, limHaloPT is the first package that includes detailed modeling of the line power spectrum.  


# Acknowledgements

The research of A.M.D. is supported by the SNSF project "The  Non-Gaussian  Universe  and  Cosmological Symmetries", project number:200020-178787. A.M.D., further acknowledges partial funding from Tomalla Foundation for research in Gravity. 


# References



