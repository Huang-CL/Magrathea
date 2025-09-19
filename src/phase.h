#ifndef PHASE_H_
#define PHASE_H_

#include "EOS.h"
#include "EOSlist.h"

/*
  Available Phase Diagrams
  Core: "Fe_default" : hcp and Liquid iron
         "Fe_fccbcc" : Includes low pressure fcc and bcc iron
  Mantle: "Si_default" : Upper Mantle: Fo, Wds, Rwd, and liquid ; Lower Mantle: Brg, PPv
          "Si_simple" : Brg, PPv, and liquid 
          "PREM" : PREM tabulated mantle
          "C_simple": Graph, Diam
  Hydrosphere: "water_default" :  H2O Water/Ice boundaries primarily form Dunaeva et al. 2010
               "water_tabulated" : AQUA Haldemann et al. 2020 Tabulated Ice, Liquid, Vapor, Supercritical
  Atmosphere: "gas_default" : Ideal Gas: Isothermal for P<100 bar. Adiabatic ideal gas for P > 100 bar
              "HHe_tabulated" : H/He Gas: Isothermal ideal for P<100 bar, P>100 bar: tabulated real gas, Chabrier & Debras 2021 Y=0.275
*/

struct PhaseDgm
/*
  At lower pressure, using the input function to set the phase diagram and pick the correct phase.  This part of the phase diagram is fixed and can't be adjusted.
  At high pressure end, assuming the phase boundary is independent to the temperature.  The phase change pressures associated with phases at low ended are listed with the phase name.  The phase and the phase change pressure can be adjusted with a function.
*/
{
  PhaseDgm(string Comp_name, EOS* (*f)(double, double), int k=0, EOS** phase_name=NULL, double *start_pressure=NULL);
  PhaseDgm(const PhaseDgm &);

  ~PhaseDgm();

  string getname(){return Comp_type;}
  EOS* (*phase_lowP)(double P, double T); // P in cgs
  void set_phase_highP(int k, double *start_pressure, EOS** phase_name);
// start pressure is an array with dimension k-1.  The first phase will replace the one with the same name in the phase diagram.

  EOS* find_phase(double P, double T);
  // Pressure in microbar

  EOS* find_phase_boundary(double Pl, double Pu, double Tl, double Tu, bool inward, double &Po, double &To, double &rhoo, double &Pn, double &Tn, double &rhon); // Used when integrate adiabatic profile across the phase boundary.  Given the pressure in cgs, integration direction, the lower and upper pressure and temperature limit (at previous and next integral step depends on the integration direction).  Return the pressure, temperature, density at old (o) phase boundary and new (n) phase boundary.

private:
  string Comp_type;
  int n;			// number of high pressure phases.
  EOS **phase_list;		// an array of pointer to classes.
  double *start_P;	      // The phase change pressure at low end. In unit of GPa.  with the size of n-1
};
  
// modular conditonal phase diagrams for each layer
//Core
EOS* find_phase_Fe_default(double P, double T);
EOS* find_phase_Fe_fccbcc(double P, double T);

//Mantle
EOS* find_phase_Si_default(double P, double T);
EOS* find_phase_Si_simple(double P, double T);
EOS* find_phase_PREM(double P, double T);
EOS* find_phase_C_simple(double P, double T);

// Hydrosphere
EOS* find_phase_water_default(double P, double T);
EOS* find_phase_water_tabulated(double P, double T);
EOS* find_phase_water_legacy(double P, double T);

//Atmsophere
EOS* find_phase_gas_default(double P, double T);
EOS* find_phase_HHe_tabulated(double P, double T);

EOS* find_phase(double m, double MC, double MM, double MW, double MG, double P, double T, bool inward = false);
// given the accumulated mass (in g), P (in cgs) and T, return the corresponding phase.  If inward = true, return the outer component at the mass transition boundary
EOS* find_phase(double m, vector<PhaseDgm>& Comp, vector<double> M, double P, double T, bool inward = false);
// given the accumulated mass (in g), P (in cgs) and T, return the corresponding phase.  If inward = true, return the outer component at the mass transition boundary

struct phase_params
{
  double Pl, Pu, Tl, Tu;
  EOS* Phase;
  PhaseDgm* cmpn;
};

extern PhaseDgm core, core1, mant, mant1, mant2, mant3, water, water1, water2, atm, atm1;
// Atmosphere doesn't support self-consistent phase diagram.  Multiple components with each masses specified have to be constructed in order to use multi-phase atmosphere.

#endif	// PHASE_H_
