#ifndef HYDRO_H_
#define HYDRO_H_

#include "define.h"
#include "EOS.h"
#include "phase.h"

extern const double ode_eps_rel0;
extern const double ode_eps_rel1;
extern const double ode_eps_rel2;
extern const double R_eps_rel;
extern const double P_eps_rel;
extern const double fit_eps_rel;
extern int fit_iter;
extern int count_shoot;
extern int count_step;


struct hydro			// pressure, density of each species, and maybe temperature against radius
{
  hydro(double Rp, vector<PhaseDgm>& Comp, vector<double> Mass_Comp, vector<double> Tgap, double ode_tol, double P0, bool isothermal); // Given the test planet radius, a vector of components, their masses, and the temperature discontinuity between each gaps (the last temperature in the list is the planet equilibrium temperature). Bool isothermal determines whether use isothermal temperature profile or self-consistent calculation. The function integrates hydrostatic equation outside in.
  hydro(double Pc, double MCin, double MMin, double MWin, double P0, double Teq=300);
// Given the center pressure, Core (iron) mass, mantle (Si) mass, water mass, integrate hydrostatic equation inside out against mass.  The integration end when either the target mass reached or pressure smaller than P0.
  hydro(double Rp, double Pc, double Tc, vector<PhaseDgm>& Comp, vector<double> Mass_Comp, vector<double> Tgap, bool isothermal, double Mfit, double ode_tol, double P0, double &Ri, double &Pi, double &Ti, double &Ro, double &Po, double &To);
  // After finishing a round of iteration of outside-in integration and receive a reasonable value of Rp, Pc, and Tc, using these value as the initial guess to integrate equations from both surface and the center and meet at the fitting mass (in g).
  // The R, P, and T value at inner edge (from inside-out integration with subscript i) and outer edge (from outside-in integration with subscript o) at the connection point is passed out. 
  void print(string outfile, bool debug=false);
  double totalM(){return M.back();}		// return the total mass of a planet model
  double totalR(){return rb.back();}		// return the total radius of a planet model
  double getP(int layer){return P[layer];}		// return the pressure
  double getM(int layer){return M[layer];}		// return the enclosed mass, M[0] may not be 0 if the initial radius is too small.
  double getR(int layer){return rb[layer];}		// return the radius, rb[0] may not be 0 if the initial radius is too large.
  double getrho(int layer){return rho[layer];}		// return the density
  double getT(int layer){return T[layer];}			// return temperature
  int getsize(){return rb.size();}		// return the number of layers
  int getLayer_from_r(double r) const;	        // return the layer index (from 0, count from bottom) by given the radius in RE. rb(l)<=r*RE<rb(l+1)
  int getLayer_from_m(double m) const;	        // return the layer index (from 0, count from bottom) by given the mass in MEarth. M(l)<=m*MEarth<M(l+1)
  vector<double> getRs();
  vector<double> getTs();	// return the temperatures at the outer side of each component interfaces as well as planet surface 
  string checkdummy();		// Check if there is any Dummy EOS used in the profile
  void setstatus(int s);
  int getstatus(){return ode_status;}
  
private:
  vector<double> rb, P, M, T, rho; // all the value at grid points.
  vector<EOS*> Phaselist;
  vector<double> M_Comp, R_Comp;
  vector<PhaseDgm> Comp;
  double fitM;
  int ode_status;			// 0 legit result, 1 dummy EOS used, 2 not converged, 3 for two layer mode, when the code stops when the center pressure meets it accuracy target before the surface pressure reaches within 2%.
};

double R_hydro(double Rp, void *params); // given the planet test radius, output the radius where mass reaches 0.

hydro* Rloop(vector<PhaseDgm> &Comp, vector<double> Mass_Comp, vector<double> ave_rho, vector<double> Tgap, double P0, bool isothermal, double &Rp, double &Pc, double &Tc);			// First round of iteration.  Input: list of componsations, masses and typical densities of each layers, the temperature discontinuity between each gaps (the last temperature in the list is the planet equilibrium temperature), boolean isothermal determines whether use isothermal temperature profile or self-consistent calculation.  The function iterates planet radius in order to get zero mass at the center.  The suggested planet radius, central pressure and temperature are also returned.

int fitting_error(const gsl_vector *x, void *params, gsl_vector *f); // return the fitting error of M, P, and T given by the integrations from both end at half radius.

hydro* fitting_method(vector<PhaseDgm>& Comp, vector<double> M_Comp, vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal);	// Using the fitting method to integrate the model from both the center and the surface.  Meet at an intermediate mass.  The precedure requires an accurate initial condition, run a Rloop first to determine a good initial condition.

double P_hydro(double Pc, void *params);
// given the center pressure, output the difference between P0 and surface pressure.

hydro* getmass(double MC, double MM, double MW, double P0);	// Pick a default Pc value. 
// Return a pointer to a hydro object.  So the object lifetime beyond the scope of the function and the memory can be released elsewhere.

struct mass_params
{
  vector<double> x;
  vector<double> M;
  EOS* Phase;
};

struct loop_params
{
  vector<double> x;
  vector<PhaseDgm> Comp;
  vector<double> M;
  vector<double> Tgap;
  bool iso;
  hydro* model;
};

#endif // HYDRO_H_
