---
title: 'Magrathea v2: A planetary interior modeling platform in C++'
tags:
  - C++
  - astronomy
  - planets
  - exoplanets
  - interiors
authors:
  - name: David R. Rice
    orcid: 0000-0001-6009-8685
    corresponding: true
    affiliation: "1, 5"
  - name: Chenliang Huang
    orcid: 0000-0001-9446-6853
    affiliation: 2
  - name: Robert Royer
    orcid: 0000-0002-0900-0192
    affiliation: 3
  - name: Mangesh Daspute
    affiliation: 4
  - name: Krishang Mittal
    affiliation: 5
  - name: Baptiste Journaux
    orcid: 0000-0002-0957-3177
    affiliation: 6
  - name: Jason H. Steffen
    orcid: 0000-0003-2202-3847
    affiliation: 3
  - name: Allona Vazan
    orcid: 0000-0001-9504-3174
    affiliation: 1
affiliations:
 - name: Astrophysics Research Center (ARCO), Department of Natural Sciences, The Open University of Israel, Raanana 4353701, Israel
   index: 1
 - name: Shanghai Astronomical Observatory, Chinese Academy of Sciences, Shanghai 200030, People’s Republic of China
   index: 2
 - name: Department of Physics and Astronomy, University of Nevada, Las Vegas, 4505 South Maryland Parkway, Las Vegas, NV 89154, USA
   index: 3
 - name: Department of Physics, Ariel University, Ariel 40700, Israel
   index: 4
 - name: Department of Astronomy, University of Wisconsin–Madison, 475 N. Charter Street, Madison, WI 53706, USA
   index: 5
 - name: Department of Earth and Space Science, University of Washington, Seattle, WA 98195, USA
   index: 6

date: 24 April 2026
bibliography: paper.bib
---

# Summary

Magrathea is an open-source C++ code for modeling the internal structure of differentiated planets. The initial release, @Huang:2022, introduced the base solver, a modular framework for defining equations of state (EOS) within phase diagrams for each layer, and a plan for expanding the code. Many of those updates are now implemented. Magrathea v2 is a more versatile platform that supports a wider range of compositions, adds new tools for composition retrieval, and makes it easier for users to adapt the code to their own models.

# Statement of need

Constraining a planet’s composition is essential for understanding its formation and evolution. Observations of mass and radius alone are not sufficient, since many different interiors can yield the same bulk density. Interior structure solvers are therefore essential tools for constraining possible compositions. With observational programs routinely measuring the densities of small to large planets, researchers require codes with models that are transparent and flexible which can adapt as our understanding of planet mineralogy changes.

Magrathea is designed as such a platform. Rather than enforcing a fixed planet model, Magrathea provides a framework in which users can define their own phase diagrams, equations of state (EOS), and thermal profiles. The goal is not one preferred interior model, but a fast and readable code base for building and testing many of them. With the continued expansion of the physics and usability in Version 2.0, Magrathea helps the community keep up with the growing precision of exoplanet observations and experimental constraints on planetary materials.

# State of the field

The field now includes several open tools for planet interior structure modeling, with broader summaries and comparisons given in @Acuna:2025 and @Baumeister:2025. Two open-source options are GASTLI [@Acuna:2025] which targets volatile-rich planets with coupled interior--atmosphere modeling and ExoPlex [@Unterborn:2023] which focuses on rocky-planet mantle minearology.

Magrathea gives users explicit control over phase diagrams and EOS choices in each differentiated layer allowing for diverse planet models from sub-Earth to Neptune-mass planets. Magrathea allows users to swap mineral physics assumptions quickly, run forward models fast enough for large sweeps, and test how interior assumptions propagate into inferred compositions. The combination of modularity and speed form the C++ implementation is the main reason we seek to build on Magrathea.

# Software design

The core solver of Magrathea is a one-dimensional, spherically symmetric integrator of the equations of hydrostatic equilibrium, mass continuity, temperature gradient, and equation of state. For a user-defined planet consisting of up to four differentiated layers, the code integrates inward and outward solutions using a shooting-to-fitting-point method with adaptive Runge--Kutta--Fehlberg stepping. The solver returns the radius of the planet, the radii of each compositional boundary, and profiles of pressure, temperature, density, and phase as functions of enclosed mass. Solving one planet takes about one second for most configurations.

A key design choice is modularity. A large variety of EOS forms are supported in `EOS.cpp`, including Birch--Murnaghan, Vinet, Holzapfel, Keane, ideal and van der Waals gases, Debye or Einstein thermal terms, and tabulated EOS. Parameters for each material are defined in a library of more than 70 EOS in `EOSlist.cpp`. Phase diagrams for each layer define which material is used at a given pressure--temperature condition in `phase.cpp`. Alternative phase diagrams are also stored in a library and can be selected at run time.

Magrathea separates the solver from the model. The code offers nine run modes through human-readable `.cfg` files, including a full four-layer solver, isothermal and two-layer modes, bulk runs, two/three-layer composition finder, on-the-fly EOS modification modes, and an MCMC retrieval mode. This range of models and use cases makes Magrathea not just a solver but a platform for exploring interior models.

# Major updates in this version

Since the initial release [@Huang:2022], Magrathea has undergone expansions in physics, solvers, and usability.

**New physical models and materials**

- **Default Mantle:** Added upper-mantle polymorphs of Mg$_2$SiO$_4$ (forsterite, wadsleyite, ringwoodite) [@Dorogokupets:2015], see \autoref{fig:phases}.
- **Default Hydrosphere:** Updated H$_2$O EOS and phase boundaries for ices [@Journaux:2020], liquid and gas [@Wagner:2002], and supercritical water [@Mazevet:2019], largely inspired by the AQUA package [@Haldemann:2020], see \autoref{fig:phases}.
- **Additional Gas EOS:** Including the solar-metallicity hydrogen/helium table from @Chabrier:2021 and van der Waals gases.
- **Carbon Mantles:** EOS and phase diagrams for phases of carbon [@Lowitzer:2006; @Benedict:2014] and silicon carbide [@Miozzi:2018], see \autoref{fig:phases}.
- **EOS library growth:** Including the AQUA table [@Haldemann:2020], fcc- and bcc-iron [@Dorogokupets:2017], and mantle materials from @Stixrude:2011.

**New functionality and solvers**

- **Composition finders:**  
  - A secant-method routine that determines the mass of a third unknown layer given a target mass, radius, and ratio between the other two layers, looped over layer ratios and mass--radius posterior draws.  
  - A Markov chain Monte Carlo routine following @Rogers:2010 and @Dorn:2015 for probabilistic composition inference given mass, radius, and associated uncertainties with the Metropolis--Hastings method.
- **Tabulated EOS:** Support for tabulated $P$--$T$--$\rho$--$\nabla T_S$ EOS tables using bilinear interpolation.
- **Modular phase diagrams:** Users can store multiple phase-diagram configurations and call them in the configuration file, for example switching between silicate and carbon mantle models without recompiling.

**Usability**

- **Input handling:** All input parameters were moved to `run/*.cfg` files with descriptive keys.
- **Parallelization:** Bulk runs and composition finder routines can exploit OpenMP in `compfind.cpp` enabling execution with multiple threads.
- **Diagnostics:** More informative error messages are returned when solutions fail to converge.
- **Tutorial and documentation:** A guided set of examples and practice problems resides in the `docs/` folder with online documentation at [magrathea.readthedocs.io](https://magrathea.readthedocs.io)

Together, these changes make Magrathea v2 a substantially expanded platform rather than a simple update to the default planet model. 

![New phase diagrams in the code. Left, updated default hydrosphere. Center, default mantle with lower-pressure Mg$_2$SiO$_4$ phases. Right, carbon and SiC mantle phase diagrams. On each plot are shown the pressure--temperature conditions inside a one Earth-mass planet with different outer temperatures. \label{fig:phases}](phase_panels.pdf)

# Research impact statement

Magrathea has been used in published work to generate mass--radius diagrams [@Taylor:2025], infer interiors of observed planets [@MacDonald:2022; @Desai:2024], connect theoretical composition from formation to observables [@Childs:2023; @Dou:2024; @Steffen:2025], and test how updated high-pressure EOS measurements and interior assumptions change inferred structures [@Huang:2021; @Lozovsky:2026]. The new composition finders are used and described in @Daspute:2025, @Rice:2025, and @Kroft:2026. 

Beyond its use in published research, Magrathea is already contributing to emerging intercomparison efforts in exoplanet interior modeling. @Schulze:2026 compares EOS and material choices across common rocky-planet models and shows that those choices can change inferred compositions at a level comparable to current observational uncertainties. In that context, open and modular tools like Magrathea help make differences in physical assumptions, EOS choices, and retrieval workflows easier to isolate and test.

Magrathea now covers a broader set of interior assumptions, supports faster composition inference workflows, and makes EOS-level sensitivity studies easier to perform. These are the main research impacts of this version.

# AI usage disclosure

Generative AI was not used to implement most of the updates described in this paper. Limited use of Generative AI was explored during development of the MCMC mode. All AI-assisted code changes were reviewed by the authors against the existing code base, and outputs were tested against other non-AI retrieval modes. Generative AI was also used in drafting and reorganizing this manuscript. All scientific claims and software descriptions were checked and revised by the authors against prior publications and the broader literature.

# Acknowledgements

We acknowledge the many researchers using the Magrathea code base; we appreciate each contribution to science and the code. We thank Douglas Adams for the timeless stories compiled in "The Hitchhiker's Guide to the Galaxy" which inspired the name of the code. We acknowledge support from the College of Sciences and the Nevada Center for Astrophysics at the University of Nevada, Las Vegas. A.V. acknowledges support by ISF grants 770/21 and 773/21. C.H. is sponsored by Shanghai Pujiang Program (grant No. 23PJ1414900).

# References
