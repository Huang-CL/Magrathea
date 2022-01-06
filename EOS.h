#ifndef EOS_H_
#define EOS_H_

#include "define.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_deriv.h>

extern const double mu;
extern const double rho_eps_rel;
extern const double T_eps_rel;
extern const bool verbose;

struct EOS
// The pressure unit inside this package is GPa
{
  EOS();
  EOS(string phaseinput, double params[][2], int length);
  EOS(string phaseinput, string filename); // construction EOS from interpolate an input file
  EOS(string phaseinput, double (*f)(double P, double T), double (*g)(double rho, double T)=NULL); // construct EOS from external functions
  EOS(string phaseinput, double *Plist, double *rholist, int len_list); // construction EOS from interpolate an input pressure density list
  EOS(string phaseinput, double params[][2], double bparams[], int length, int blength); // construction EOS for RTpress
  ~EOS();
  
  void setphasename(string phaseinput){phasetype=phaseinput;}
  void modifyEOS(double params[][2], int length=3);  // modify the constructed EOS parameters
  void modifyEOS(int index, double value);	     // modify one value of the EOS
  void modify_extern_density(double (*f)(double P, double T)){density_extern = f;}
  void modify_extern_entropy(double (*g)(double rho, double T)){entropy_extern = g; thermal_type = 1;}
  void modify_dTdP(double (*h)(double P, double T)){dTdP = h; thermal_type = 2;} // If the dTdP is set, it will overwrite the entropy method. Pressure in cgs unit.
  
  double BM3(double rho);	// input rho in g/cm^3, return pressure in GPa, type 0
  double BM4(double rho);	// type 1
  double Vinet(double rho);	// type 2
  double Holzapfel(double rho);	// type 3
  double Keane(double rho);	// type 4
  // type 5, Water, low pressure ice EOS from Choukroun & Grasset
  // type 6, ideal gas law.
  // type 7, interpolate an input file
  // type 8-13 the same as 0-5 but for RTpress style. 8 BM3, 9 BM4, 10 Vinet, 11 Holzapfel, 12 Keane, 13 Choukroun
  
  double Pth(double rho, double T); // thermal pressure in GPa, and electron pressure or anharmonic if provided.
  double adiabatic_index();	    // get the adiabatic index for ideal gas.  Vibrational freedom is always ignored.
  double density(double P, double T, double rho_guess); // input P in cgs (microbar), at given temperature, return density in g/cm^3
  double (*density_extern)(double P, double T);		// using external function, input P, T in cgs, return density in g/cm^3
  double (*entropy_extern)(double rho, double T);	// using the external entropy function.
  void printEOS();					// print EOS table into a file named with phasename.txt
  string getEOS(){return phasetype;}		// get the type of EOS
  int geteqntype(){return eqntype;}		// get the type of EOS fitting eqn
  double getV0(){return V0;}
  double getmmol(){return mmol;}
  double getP0(){return P0;}
  double getT0(){return T0;}
  double (*dTdP)(double P, double T); // return the temperature gradient at given pressure and temperature point. Pressure in cgs unit.
  void DebyeT(double x, double &gamma, double &Theta);	   // return the Grueneisen parameter, Debye temperature or Einstein temperature according to Altshuler form.  If Theta0 is not available, a Debye temperature scaling factor is returned
  double entropy(double rho, double T); // Given the volume per mol and temperature, calculate the entropy over n*R, or P V^{7/5} / R for ideal gas.
  double pSpV_T(double V, double T);
  // partial S (entropy) partial V at constant T
  double pSpT_V(double V, double T);
  // partial S (entropy) partial T at constant V
  double pPprho_T(double rho, double T);
  // partial P partial rho at constant T in GPa / g/cm^3
  double pPpT_rho(double rho, double T);
  // partial P partial T at constant rho in GPa / K
  double dTdm(double m, double r, double rho, double T);
  // partial T partial enclosed mass
  double dTdP_S(double P, double T, double &rho_guess);
  // partial T partial P along isentrope in K / GPa, given pressure in GPa
  int getthermal(){return thermal_type;}	

  double fV (double V) {return 0.5*(pow(V0/V,2./3.)-1);}
// finite volumetric strain, take volume in cm^3 / mol
  double gamma0S (double V);
  // Grueneisen parameter along the reference adiabat (eq A.3), take volume in cm^3 / mol
  double bV(double V);
  // thermal coefficients b(V) (Eq.10) in eV/atom, take volume in cm^3 / mol
  double bVp(double V);
  // derivative of b(V) (Eq. B.2) eV mol/(cm^3 atom), take volume in cm^3 / mol
  double TOS(double V);
  // reference adiabat temperature profile (Eq. 7), take volume in cm^3 / mol
  double fT (double T) {return pow(T/T0, beta) -1;}
  // thermal deviation from the reference temperature (Eq. 9)
  double fTp (double T) {return beta/T0*pow(T/T0, beta-1);}
  // thermal deviation from the reference temperature (Eq. 9)
  double Cv (double V, double T);
  //total heat capacity in eV/atom/K (Eq. B.4), take volume in cm^3 / mol
  double Spot(double V, double T);
  // potential contribution of entropy in eV/atom/K, take volume in cm^3 / mol
  double gamma (double V, double T);
  // Grueneisen parameter (Eq. 17), take volume in cm^3 / mol
  double Press (double rho, double T);
  // pressure in GPa (Eq. 6, 13, 14), take density in g/cm^3
  double dTdV_S(double V, double T);
  // adiabatic temperature gradient in K mol/cm^3, take volume in cm^3 / mol
  double density(double V) {return mmol/V;}
  // return density in g/cm^3
  double volume(double rho){return mmol/rho;}
  // return volume in g/cm^3, take density in g/cm^3
  double density(double P1, double T1, double rho_guess, double P2, double &T2);
  // Given the pressure (cgs), temperature, density of the previous step, the pressure of the next step, return the temperature and density at the new pressure. This solver doesn't conserve the entropy well enough. Only used as an approximation in the first integration step from the core of the planet where dTdm has 0/0 limit.

private:
  string phasetype;
  int eqntype;
  double V0, K0, K0p, K0pp, mmol, P0, Cp, Theta0, gamma0, beta, gammainf, gamma0p, e0, g, T0;
  double at1, at2, at3, at4, ap1, ap2, ap3, ap4;
  int n, Z;
  bool Debye_approx;		       // Debye approximate or Einstein approximate.
  int thermal_type;		       // Indicates the thermal type of the phase.  0 indicates no temperature profile available, 1 indicates entropy method, 2 indicates the temperature gradient method, 3 indicates ideal gas, 4 indicates the EOS is fitted along the isentrope, 5 indicates no Theta0, 6 indicates has Theta 0 but no electron pressure, 7 indicates has electron pressure as well, type 8, RTpress style
  double *rhotable, *Ptable;	// density table in cgs, Ptable in GPa.
  int bn;				// number of indices of b
  double* b;				// fitted polynomial parameters of the thermal coefficients b(V) in eV/atom

  gsl_interp_accel *acc; // The gsl interpolation accelerator.
  gsl_spline *spline; // The gsl workspace objects
  int nline;
  
/*
  phasetype is the name of a phase. The comment about the EOS used for the phase should be in the parentheses separated by a space.
0.	eqntype. 8 for RTpress style.
1.	V0 in cm^3 / mol. For RTpress style, the input is in \AA^3 /atom and the number would be converted to cm^3/mol by the constructor
	For ice, 1 \AA^3 = N_A / (2*10^24) cm^3/mol = 0.3011 cm^3/mol
2.	K0 in GPa
3.	K0p
4.	K0pp in GPa ^-1
5.	mmol in g / mol, or mean molecular weight of gas, or in g / mol for RTpress style
6.	P0 (GPa) the minimum pressure.  The pressure correspond to V0
7.	Cp in 10^7 erg/g/K Heat capacity at constant pressure
8.	Theta0 (K), a fitting parameter to Einstein temperature or Debye temperature
9.	gamma0, a fitting parameter of Grueneisen parameter
10.	beta, a fitting parameter of Grueneisen parameter.  In RTpress style, it represents the "m" which stands for the power-law exponent in the thermal deviation term. theoretically expected value: 0.6
11.	gammainf, a fitting parameter of Grueneisen parameter
12.	gamma0p, volume derivative of the Grueneisen parameter
13.	e0 (10^-6 K^-1), electronic contribution to Helmholtz free energy
14.	g, is an electronic analogue of the Grueneisen parameter
15.	n is the number of atoms in the chemical formula of the compound.  Should have n*NA atoms within V.  The n of ideal gas is the number of atoms per molecule for the purpose of adiabatic index.  NOTE: n=2 for collinear molecules e.g. carbon dioxide!  Isothermal atmosphere can be achieved by setting n=0.
16.     Z is the atomic number (number of electron)
17.	T0, the reference temperature for the thermal pressure
18.	Debye_approx, whether use Debye approximation or Einstein approximation. Debye approximation is slower but more accurate at temperature lower than Debye/Einstein temperature.  Positive number for Debye, otherwise Einstein.
19.     thermal_type, indicates the thermal type of the phase.  0 indicates no temperature profile available, 1 indicates entropy method, 2 indicates the temperature gradient method.  The only method to set the gradient is using the modify_extern_dTdP function, 3 indicates ideal gas, 4 indicates the EOS is fitted along the isentrope, type 8 indicates RTpress style.
20-27.  at1-at4 & ap1 - ap4
 */
};

double P_EOS(double rho, void *params);
double dP_EOS(double rho, void *params);
void PdP_EOS(double rho, void *params, double *P, double *dP);

struct EOS_params
{
  vector<double> x;
  EOS* Phase;
};

#endif // EOS_H_