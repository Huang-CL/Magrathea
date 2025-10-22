---
title: 'Magrathea v2: planetary interior structure code in C++'
tags:
  - C++
  - astronomy
  - planets
  - exoplanets
  - interiors
authors:
  - name: David R. Rice
    orcid: 0000-0001-6009-8685
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Chenliang Huang
    orcid: 0000-0001-9446-6853
    affiliation: 2
  - name: Robert Royer
    orcid: 0000-0002-0900-0192
    affiliation: 3
  - name: Mangesh Daspute
    orcid: 
    affiliation: 4
  - name: Krishang Mittal
    orcid: 
    affiliation: 5
  - name: Baptiste Journaux
    orcid: 0000-0002-0957-3177
    affiliation: 6
  - name: Jason H. Steffen
    orcid: 0000-0003-2202-384
    affiliation: 3
  - name: Allona Vazan
    orcid: 0000-0001-9504-3174
    affiliation: 1
affiliations:
 - name: Astrophysics Research Center (ARCO), Department of Natural Sciences, The Open University of Israel, Raanana 4353701, Israel
   index: 1
 - name: Shanghai Astronomical Observatory, Chinese Academy of Sciences, Shanghai 200030, People’s Republic of China
   index: 2
 - name: Department of Physics and Astronomy, University of Nevada, Las Vegas, 4505 South Maryland Parkway, Las Vegas, NV 89154, US
   index: 3
 - name: Department of Physics, Ariel University, Ariel 40700, Israel
   index: 4
 - name: Department of Astronomy, University of Wisconsin–Madison, 475 N. Charter Street, Madison, WI, 53706, USA
   index: 5
 - name: Department of Earth and Space Science, University of Washington, Seattle, WA 98195, United States
   index: 6 

date: 23 October 2025
bibliography: paper.bib


# Summary

Magrathea is an open-source C++ code for modeling the internal structure of differentiated planets. The initial release, @Huang:2022, introduced the base solver, a modular framework for defining equations of state (EOSs) used within phase diagrams for each differentiated layer, and outlined a series of planned extensions. Many of those updates are now implemented here, alongside new physics, solvers, and usability improvements. The result is a more versatile platform that supports a winder range of planetary compositions, adds new tools for composition retrieval, and makes it easier for users to adapt the code to their own models.


# Statement of need

Constraining a planet’s composition is essential for understanding its formation and evolution. Observations of mass and radius alone are not sufficient, since many different interiors can yield the same bulk density. Interior structure solvers are therefore essential tools for exploring this degeneracy and constraining the possible distributions of iron, silicate, volatiles, and atmosphere. With observational programs routinely measuring the densities of small to large planets, researchers require codes with models that are transparent and flexible---able to adapt to our changing understanding of planet mineralogy. 

Magrathea is designed as such a platform. Rather than enforcing a fixed planet model, Magrathea provides a framework in which users can define their own phase diagrams, equations of state (EOSs), and thermal profiles. This adaptability has led to broad uptake: Magrathea has been used to generate mass–radius diagrams and infer interiors of observed planets [@Desai:2024; @Daspute:2025; @Rice:2025; @Taylor:2025], to connect theoretical composition models to observables [@Steffen:2025; @Dou:2024; @Childs:2023], and to incorporate new high-pressure equation of state (EOS) measurements into planetary modeling [@Huang:2021]. A list of other open source interior models can be found in @Acuna:2025. With the continued expansion of the physics and usability, Magrathea helps the community keep up with the growing precision of exoplanet observations and experimental constraints on planetary materials.

# Summary of the base code

The core solver of Magrathea is a one-dimensional, spherically symmetric integrator of the equations of hydrostatic equilibrium, mass continuity, and energy transport written in C++. For a user-defined planet consisting of up to four differentiated layers, the code integrates inward and outward solutions using a shooting-to-fitting-point method with adaptive Runge–Kutta–Fehlberg stepping. The solver returns the radius of the planet, the radii of each compositional boundary, and profiles of pressure, temperature, density, and phase as functions of enclosed mass. The hydrostatic integration routines are implemented in `hydro.cpp`. 

A key design choice is **modularity**:
- A large variety of EOS forms are supported in `EOS.cpp`, including Birch–Murnaghan, Vinet, Holzapfel, Keane, van der Waals gases, and tabulated.
- Parameters for EOSs are defined and stored in a large library (70+ EOSs) in `EOSlist.cpp`.  
- Phase diagrams for each layer define which material is used at a given P-T condition in `phase.cpp`.

MAGRATHEA offers **nine run modes** through human-readable `.cfg` files:  
1. **Full solver** takes masses for each layer and returns the planet’s radius and interior profiles.  
2. **Temperature-free solver** for isothermal interiors.  
3. **Two-layer mode** for rapid mass–radius curves.  
4. **Bulk mode** for ensembles of planets.  
5. **Composition finder**: determine an unknown layer mass to match observed M and R using a secant method.  
6. **On-the-fly EOS modification** for testing parameter uncertainties.  
7. **Iterated EOS modification** with two-layer solver.  
8. **Iterated EOS modification** with full solver.  
9. **MCMC composition retrieval** for probabilistic inference given mass, radius, and corresponding uncertainties.  

This modularity and range of modes makes Magrathea not just a solver but a platform for exploring interior models.


# Major updates in this version

Since the initial release [@Huang:2022], Magrathea has undergone expansions in physics, solvers, and usability. We summarize the key updates here.

**New physical models and materials**
- **Default Mantle:** Added upper-mantle polymorphs of Mg\(_2\)SiO\(_4\) (forsterite, wadsleyite, ringwoodite) [@Dorogokupets:2015], see \autoref{fig:phases}.
- **Default Hydrosphere:** Updated H\(_2\)O EOSs and phase boundaries with ices [@Journaux:2020], liquid, gas [@Wagner:2002], and supercritical [@Mazevet:2019, with 2021 entropy correction] largely inspired by the AQUA package [@Haldemann:2020], see \autoref{fig:phases}.
- **Additional Gas EOSs:** Including the solar-metalicity table for hydrogen/helium from [@Chabrier:2021] and van der Waals gases.
- **Carbon Mantles:** EOSs and phase diagrams for phases of carbon [@Benedict:2014] and silicon carbide [@Miozzi:2018], see \autoref{fig:phases}.
- **EOS library growth:** Dozens of additional EOSs: AQUA table [@Haldemann:2020], fcc- and bcc-iron [@Dorogokupets:2017], and the mantle materials from @Stixrude:2011.

**New functionality and solvers**
- **Composition finders:**  
  - A secant-method routine that determines the mass of a third unknown layer given a target mass, radius, and ratio between the other two layers looped over layer ratios and mass and radius posterior draws.
  - An Markov chain Monte Carlo based routine following @Rogers:2010 and @Dorn:2015 for probabilistic composition inference given mass, radius, and associated uncertainties with Metropolis–Hastings method.
- **Tabulated EOSs:** Support for tabulated $P$–$T$–$\rho$–$\nabla T_S$ EOS tables using bilinear interpolation.
- **Modular phase diagrams** Allow users to store multiple phase-diagram configurations and call them in the configuration file—for example, toggling between a silicate-based and carbon-based mantle phase diagram.

**Usability**
- **Input handling:** All input parameters migrated to `run/*.cfg` files with descriptive keys and documentation, improving reproducibility and scripting.
- **Parallelization:** Bulk runs and composition finder routines can exploit OpenMP in `compfind.cpp`, enabling execution with multiple threads across cores.
- **Diagnostics:** More informative error messages when solutions fail to converge.
- **Tutorial and documentation:** A guided set of examples and practice problems now resides in the `docs/` folder, and an online documentation site is available at [magrathea.readthedocs.io](https://magrathea.readthedocs.io).

Together, these updates make Magrathea a platform for both forward modeling and statistical inference of exoplanet interiors. By expanding the physics library, adding composition retrieval solvers, and improving usability, the code now enables a wider range of applications than in its initial release.  

Development of Magrathea is ongoing. Planned future expansions include building versatile methods for mixing materials, adding treatments of thermal evolution, and coupling to an atmosphere model. By continuing to integrate new experimental and theoretical results, Magrathea will remain a robust and adaptable tool for interpreting the increasing number and precision of exoplanet observations.

![New phase diagrams in the code._Left, default hydrosphere compiled from many sources---A: low pressure ice/liquid [@Journaux:2020], B: ice-VII [@Bezacier:2014,@Sotin:2007], C: ice-X [@Grande:2022], D: IAPWS-95 liquid/gas [@Wagner:2002], E: supercritical [@Brown:2018], F: van der Waals gas, G: supercritical [@Mazevet:2019]. Center, default mantle with lower pressure Mg\(_2\)SiO\(_4\) phases. Right, carbon phase diagram in dark and SiC phase diagram in light grey. On each plot is shown the P-T conditions inside a 100% composition planet of one Earth-mass with two or three different outer temperatures. The density inside of the planet is shown by each plot's colorbar and the radius of the planet is denoted in the legend. \label{fig:phases}](phase_panels.pdf)


# Acknowledgements

We acknowledge the many researchers using the Magrathea code base; we appreciate each contribution to science and the code. We thank Douglas Adams for the timeless stories compiled in the "Hitchikers Guide to the Galaxy" which inspired the name of the code. We acknowledge support from the College of Sciences and the Nevada Center for Astrophysics at the University of Nevada, Las Vegas. A.V. acknowledges support by ISF grants 770/21 and 773/21. C.H. is sponsored by Shanghai Pujiang Program (grant NO. 23PJ1414900).

# References
