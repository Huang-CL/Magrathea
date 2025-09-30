This directory contains example input files to run MAGRATHEA.

*inputcore.txt*, *inputmantle.txt* and *inputwater.txt* are input files for `mode3.cfg` single layer planets with masses from 0.1 to 4.0 Earth masses. Files generated with *massradiusinput.py*.

*inputternary.txt* is an input file for `mode3.cfg` containing 5151 planets of 1 Earth-mass with integer percentages of mass in each of three layers. 

*inputplanetMR.txt* is an input file for `mode4.cfg` containing mass and radii for a planet drawn from a posterior distribution from an observation.

*PosteriorsPPv.txt* input for mode 4 or 5 containing 1000 values of the reference volume, bulk modulus, and the derivative of the bulk modulus drawn from Gaussians using the mean and uncertainties from Sakai et al. 2016 for silicate post-perovskite.

*PosteriorsVinet.txt* input for mode 4 or 5 containing 1000 values of EOS parameters for Ice VII, Ice VII', and Ice X from the reported values and uncertainties in Grande et al. 2022.
