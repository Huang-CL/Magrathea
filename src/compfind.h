#ifndef COMPFIND_H_
#define COMPFIND_H_

#include "hydro.h"

void multiplanet(vector<PhaseDgm> &Comp, vector<double> Tgap, int solver, vector<double> ave_rho, double P0, bool isothermal, string infile, string outfile);
// Find planet radii and radii of layers for an input file of planets
// Input file must specify planet mass and mass fractions of each layer


void compfinder(vector<PhaseDgm> &Comp, int findlayer, vector<int> layers, double minPMR, double maxPMR, float step, double rerr, int num_threads,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, string infile, string outfile);
// Find 3-layer solutions for samples of mass and radii
// Set which layer to find and keep other two in constant partial mass ratio
// Loops through multiple partial mass ratios





#endif
