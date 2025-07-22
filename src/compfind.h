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


//Structures for MCMC storage and iteration
struct MCMCRecord {
    double Mass;       
    double fCore;        
    double fMantle;      
    double fWater;      
    double fAtm;       
    double log_likelihood;  
    double RCore;
    double RMantle;
    double RWater;     
    double RPlanet;
};

struct Params {
    double Mass;
    double fCore;
    double fMantle;
    double fWater;
};

struct LikelihoodResult {
    double log_likelihood;
    double RCore;
    double RMantle;
    double RWater;
    double RPlanet; 
};


void mcmcsample(vector<PhaseDgm> &Comp, double MassPrior, double MUncPrior, double RadPrior, double RUncPrior,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, string outfile, int numchains, int steps, int numlayers);
// Find 2, 3 or 4-layer solutions for mass and radius with uncertainty

LikelihoodResult log_likelihood(vector<PhaseDgm> &Comp, double MassPrior, double MUncPrior, double RadPrior, double RUncPrior,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, double Mass, double fCore, double fMantle, double fWater);
//Find planet radius and likelihood calculation following Rogers & Seager 2010 ApJ

void metropolis_hastings(vector<MCMCRecord>& chain,vector<PhaseDgm> &Comp, double MassPrior, double MUncPrior, double RadPrior, double RUncPrior,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, int steps, int numlayers);
// MCMC method


#endif
