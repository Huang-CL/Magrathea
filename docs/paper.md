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
date: 23 September 2025
bibliography: paper.bib


# Summary



# Statement of need

The code...

At least a dozen research groups are now or have used Magrathea, these include...

# Major updates in this version

Since publishing the first verison of the code along with CITE, the capabilities, useability, and performance of Magrathea has all been significantly improved and expanded. We separate here lists of updates to functionality and to the materials/physics available for the model. (MIGHT BE BETTER WORKDING THAN THIS) The updates to functionality:

- Input parameters for each mode are moved to input files with .cfg extensions
- Reordering of the now 9 modes see table
- Modular phase diagrams
- P-T-rho-dt/dp tables
- Secant Method composition solver
- MCMC composition solver

The updates to the model:
- Upper mantle phases
- Chambrier, AQUA, and SeaFreeze
- Supercritical Water
- Van der Waals
- Carbon Phase Diagram
- Updated Hydrosphere


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

We acknowledge contributions from 

# References
