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

date: 1 October 2025
bibliography: paper.bib


# Summary

In our first publication, @Huang:2022, we listed a number of future updates which are now completed along with a host of other expansions to the code's available physics and features.

# Statement of need

Essential to understanding the formation and evolution of a planet, is constraining the composition of the planet. Observations of mass and radius alone are not sufficient, since many different interiors can yield the same density. Interior structure solvers are therefore essential tools for connecting observed properties to the underlying distribution of iron, silicate, volatiles, and atmosphere. With observational programs routinely measuring the densities of small to large planets, researchers require codes with models that are transparent and flexible---able to adapt to our changing understanding of planet minerology. 

Magrathea is designed as such a platform. Rather than enforcing a fixed planet model, Magrathea provides a framework in which users can define their own phase diagrams, equations of state (EOSs), and thermal structures. This adaptability has led to broad uptake: at least a dozen groups now use Magrathea. DESCRIBE A FEW WORKS. By continuing to expand its physical fidelity and usability, Magrathea enables the community to keep pace with rapidly improving exoplanet data and experimental constraints on planetary materials.

# Summary of the base code

The core solver of Magrathea is a one-dimensional, spherically symmetric integrator of the equations of hydrostatic equilibrium, mass continuity, and energy transport written in C++. For a user-defined planet consisting of up to four differentiated layers, the code integrates inward and outward solutions using a shooting-to-fitting-point method with adaptive Runge–Kutta–Fehlberg stepping. The solver returns the radius of the planet, the radii of each compositional boundary, and profiles of pressure, temperature, density, and phase as functions of enclosed mass.

A key design choice is **modularity**:
- Support for a large variety of EOS forms, including Birch–Murnaghan, Vinet, Holzapfel, Keane, van der Waals gases, and tabulated EOSs are implemented in `EOS.cpp`.
- EOSs parameters are defined and stored in a large library (70+ EOSs) in `EOSlist.cpp`.  
- Phase diagrams for each layer define which material is used at a given P-T condition in `phase.cpp`.
- The hydrostatic integration routines are implemented in `hydro.cpp`.  

MAGRATHEA offers **nine run modes** through human-readable `.cfg` files:  
1. **Full solver** define mass in each layer retrieve radius and interior conditions.  
2. **Temperature-free solver** for isothermal interiors.  
3. **Two-layer mode** for rapid mass–radius curves.  
4. **Bulk mode** for ensembles of planets.  
5. **Composition finder**: determine an unknown layer mass to match observed M and R using a secant method.  
6. **On-the-fly EOS modification** for testing parameter uncertainties.  
7. **Iterated EOS modification** with two-layer solver.  
8. **Iterated EOS modification** with full solver.  
9. **MCMC composition retrieval** for probabilistic inference given mass, radius, and corresponding uncertainties.  

This diversity makes Magrathea not just a solver but a *platform* for testing planetary interior hypotheses.


# Major updates in this version

Since the initial release [@Huang:2022], Magrathea has undergone major expansions in physics, solvers, and usability. We summarize the key updates here.

**New physical models and materials**
- Inclusion of upper-mantle polymorphs of Mg\(_2\)SiO\(_4\) (forsterite, wadsleyite, ringwoodite) [CITE] FIGURE.
- A new default hydrosphere treatment with updated H2O ices [@Journaux:2020], liquid, gas [IAPWS], and supercritical [Mazevet] largely inspired by the AQUA package [CITE] FIGURE.
- Additional gas EOSs CHAMBRIER solar-metalicity table for hydrogen/helium and van der Waals gases.
- EOS and phase diagram for phases of carbon and silicate carbide.
- A number of other materials including: AQUA table, fcc and bcc iron, the mantle materials from STIXRUDE

**New functionality and solvers**
- Modular phase diagram handling, allowing users to store and call in the run different compositions (e.g. swapping a magnesium silicate mantle for a carbon mantle).
- Support for tabulated P–T–ρ–∇T EOS tables using bilinear interpolation.
- New composition finder solvers:  
  - a secant-method routine that determines the mass of an unknown layer given a target mass, radius and ratio between the other layers;  
  - an MCMC-based routine for probabilistic composition inference, usable with posterior samples from observations.

**Usability and performance**
- Migration of all run parameters into `.cfg` input files, improving reproducibility and scripting.
- Improved error handling and convergence diagnostics at phase boundaries.
- Parallelization of bulk runs and composition finder routines with OpenMP, enabling efficient inference across thousands of mass–radius samples.
- More consistent default tolerances for the ODE integrators, reducing numerical artifacts while preserving speed.

Together, these updates transform Magrathea from a flexible forward solver into a broader platform for both forward modeling and statistical inference of exoplanet interiors. The new features have already been demonstrated in recent applications to the TRAPPIST-1 system [@Rice:2025], where the composition finder enabled quantification of model-dependent uncertainties in water mass fraction.

![A few new phase diagrams in the code._Left, default mantle with lower pressure phases from INSERT. Center, default hydrosphere compiled from many sources: A: INSERT, B: INSERT. Right, carbon phase diagram in dark and SiC phase diagram in light grey. On each the P-T conditions inside a 100% composition planet of X Earth-mass with three different tempaeratures. \label{fig:example}](figure.png)


We acknowledge the many researchers using the Magrathea code base; we appreicate each contribution to science and the code. We thank Douglas Adams for the timeless stories compiled in the "Hitchikers Guide to the Galaxy" which inspired the name of the code. 

# References
