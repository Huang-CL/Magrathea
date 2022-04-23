This directory has plotting scripts written in python to help analyze the output of MAGRATHEA. Scripts are written for Python 3.6. Scripts may require matplotlib, numpy, astropy, python-ternary, regex, cmasher, and scipy.

quickdensityplot.py plots the results of result/StructureEarth.txt the output of input_mode=0. Outputs interior.pdf.

ternaryplot.py plots result/OneEarthTern.txt generated using run/inputternary.txt using input_mode=6 with 300 K surface tempearture and 0 K jumps between layers. Usings python-ternary package: https://github.com/marcharper/python-ternary. Outputs oneearthternary.pdf.
