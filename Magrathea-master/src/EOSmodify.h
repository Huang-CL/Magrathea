#ifndef EOSMODIFY_H_
#define EOSMODIFY_H_

#include "hydro.h"

void twolayer(int index, double fraction, vector<double> &Mp, vector<double> &Rp, double P0, bool printmodel=false);
// calculate a mass-radius curve for two-layer planet with a constant mass fraction.
// index=0 stands for planet with no iron core.  Only have Si-mantle and ice crust.
// index=1 for planet with no Si-mantle.  2 for no water.
// fraction is the mass fraction of the inner layer.
// return a list of planet mass Mp, planet radius Rp and the length of array n

void twolayer(int absent_index, double fraction, double P0, int adjust_index, string eosfile, string outfile);
// calculate a mass-radius curve for two-layer planet with a constant mass fraction by giving the parameters of EOSs and phases boundaries of high pressure phases.
// Format of the input file
// 1. number of phases, a list of name of each phases. If the first phase in the list already exist in the database, the name should be the same as the phasetype string (not the name of EOS pointer) exclude parentheses part, e.g. Fe hcp, Si PPv, Water etc.
// 2. for each phases, how many parameter provided, the index of each parameter according to the index table.  The number of parameters has to be the same for all phases.
// Then a list of index with their default values.  The values has to be the same for all phases.
// 3. main body.  beginning pressure, parameters for phase 1, and so on.
// absent_index=0 stands for planet with no iron core.  Only have Si-mantle and ice crust.
// absent_index=1 for planet with no Si-mantle.  2 for no water.
// fraction is the mass fraction of the inner layer.
// adjust_index. The index of the phase whose EOS variation impact are going to be investigated. 
// using vector to input a list of V0, K0, and K0p of the ice VII EOS.
// calculating the planet model with those setting.
// Storage the result in a file, first row is the mass.  
// The following rows are the radii of planets of corresponding masses, with the given EOS.


void fullmodel(vector<PhaseDgm> &Comp, vector<double> M_Comp, vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, int adjust_index, string eosfile, string outfile);
// calculate a planet radius for complete planet model by giving the mass of each components, parameters of EOSs and phases boundaries of high pressure phases.
// Format of the input file
// 1. number of phases, a list of name of each phases. If the first phase in the list already exist in the database, the name should be the same as the phasetype string (not the name of EOS pointer) exclude parentheses part, e.g. Fe hcp, Si PPv, Water etc.
// 2. for each phases, how many parameter provided, the index of each parameter according to the index table.  The number of parameters has to be the same for all phases.
// Then a list of index with their default values.  The values has to be the same for all phases.
// 3. main body.  beginning pressure, parameters for phase 1, and so on.
// adjust_index. The index of the phase whose EOS variation impact are going to be investigated. 0:iron core, 1:mantle, 2, water layer. Do not support changing atmosphere.
// calculating the planet model with those setting.
// Storage the planet radius in a file



#endif
