---
title: 'Magrathea v3: planetary interior structure code in C++'
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
 - name: INSERT MANGESH
   index: 4
 - name: Department of Astronomy, University of Wisconsin–Madison, 475 N. Charter Street, Madison, WI, 53706, USA
   index: 5

date: 26 September 2025
bibliography: paper.bib


# Summary



# Statement of need

Understanding the internal structure of exoplanets is a key step in connecting mass–radius measurements to composition, formation, and evolution. Observational programs with *Kepler*, *TESS*, and JWST now routinely measure the densities of small planets, but interior degeneracies make inference strongly dependent on the assumptions of the chosen structural model. Researchers therefore require codes that are both transparent and adaptable to new physics, equations of state, and phase diagrams. 

Magrathea is designed as such a platform. Rather than providing only a single planet model, it allows users to construct and test their own interior prescriptions. Since its first publication [@Huang:2022], Magrathea has been adopted by more than a dozen groups for applications ranging from volatile-rich super-Earths to habitability analyses [@Rice:2025; @Childs:2023]. Recent work using Magrathea has quantified uncertainties in the interiors of the TRAPPIST-1 planets [@Rice:2025], demonstrating the importance of propagating both observational and model uncertainties. By releasing the code openly and continuing to expand its physical fidelity and usability, Magrathea provides a flexible, community-driven alternative to black-box or closed-source solvers.

# Summary of the base code

The base solver of Magrathea is a one-dimensional, spherically symmetric planet interior integrator written in C++. The code supports up to four fully differentiated layers—iron core, silicate mantle, water/ice hydrosphere, and atmosphere—with user-defined mass fractions. For a given configuration, Magrathea integrates the equations of hydrostatic equilibrium, mass conservation, temperature profile (isothermal or isentropic), and equation of state (EOS) using a shooting-to-fitting-point method with adaptive Runge–Kutta integration. The solver returns the radius of the planet, the radii of each compositional boundary, and profiles of pressure, temperature, density, and phase as functions of enclosed mass.

A distinctive feature of Magrathea is its modular EOS and phase diagram library. Each layer can host multiple phases, with transitions determined by pressure–temperature conditions. Users may select from a large library of tabulated or analytic EOSs, or add new ones via straightforward interfaces in the source. This flexibility makes Magrathea well suited both for forward modeling of specific exoplanets and for exploring how uncertainties in high-pressure physics affect planetary interiors.


# Major updates in this version

Since the initial release [@Huang:2022], Magrathea has undergone significant improvements in physics coverage, numerical methods, and usability. We summarize the key updates here.

**New physical models and materials**
- Inclusion of upper-mantle polymorphs of Mg\(_2\)SiO\(_4\) (forsterite, wadsleyite, ringwoodite) and phase transitions to bridgmanite and post-perovskite MgSiO\(_3\) [@Rice:2025].
- Expanded hydrosphere treatment with updated high-pressure ices and water phases, including supercritical water and recent experimental constraints [@Huang:2021].
- Support for advanced EOS formulations, including Chabrier EOS for hydrogen/helium, AQUA and SeaFreeze water models, and van der Waals gases.
- A carbon phase diagram module enabling models of carbon-rich interiors.

**New functionality and solvers**
- Modular phase diagram handling, allowing users to mix, replace, or extend layer definitions.
- Support for tabulated P–T–ρ–∇T EOS tables, in addition to analytic forms.
- New composition finder solvers:  
  - a secant-method routine that determines the mass of an unknown layer given a target mass and radius;  
  - an MCMC-based routine for probabilistic composition inference, usable with posterior samples from observations.
- Reordered and expanded run modes (now nine in total) to cover bulk planet input, two-layer approximations, EOS modification, and inverse retrieval problems.

**Usability and performance**
- Migration of all run parameters into `.cfg` input files, improving reproducibility and scripting.
- Improved error handling and convergence diagnostics at phase boundaries.
- Parallelization of bulk runs and composition finder routines with OpenMP, enabling efficient inference across thousands of mass–radius samples.
- More consistent default tolerances for the ODE integrators, reducing numerical artifacts while preserving speed.

Together, these updates transform Magrathea from a flexible forward solver into a broader platform for both forward modeling and statistical inference of exoplanet interiors. The new features have already been demonstrated in recent applications to the TRAPPIST-1 system [@Rice:2025], where the composition finder enabled quantification of model-dependent uncertainties in water mass fraction.


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge the many researchers using the Magrathea code base; we appreicate each contribution to science and the code. We thank Douglas Adams for the timeless stories compiled in the "Hitchikers Guide to the Galaxy" which inspired the name of the code. 

# References
