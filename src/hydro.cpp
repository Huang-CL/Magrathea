#include "hydro.h"
#include "EOSlist.h"
#include <gsl/gsl_integration.h>

int derivs_m(double x, const double y[], double dydx[], void * params)
// derivative equations
// two equations, x is mass, y0 is radius, y1 is P.
// Using mass instead of radius as variable because the integral upper limit is given with mass.
// Only works for constant room temperature and no gas.
{
  struct double_params *p=(struct double_params *) params;
  double rho= p->x[0];
  double MC = p->x[1];
  double MM = p->x[2];
  double MW = p->x[3];
  double Teq = p->x[4];
  
  if (y[1] < 0)
  {
    p->x[0] = -1;
    dydx[0] = -1/(4*pi*sq(y[0]));
    dydx[1] = -G*x/(4*pi*pow(y[0],4));

    return GSL_SUCCESS;
  }
  
  EOS* Phase=find_phase(x,MC,MM,MW,0,y[1],Teq);
  if (!Phase)
  {
    if (verbose)
      cout<<"Warning: Can't find the phase when conducting constant temperature integral at location: r="<<y[0]<<"cm, m="<<x/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, Component masses "<<MC<<' '<<MM<<' '<<MW<<endl;
    return GSL_EBADFUNC;
  }
  
  rho = Phase->density(y[1],Teq,rho);
  if (!gsl_finite(rho)) 
  {
    if (verbose)
      cout<<"Warning: Can't find the density for "<<Phase->getEOS()<<" when performing constant temperature integral at location: r="<<y[0]<<"cm, m="<<x/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, Component masses "<<MC<<' '<<MM<<' '<<MW<<endl;
    return GSL_EBADFUNC;
  }
    
  dydx[0] = 1/(4*pi*sq(y[0])*rho); // dr/dm = 1/(r pi r^2 rho)
  dydx[1] = -G*x/(4*pi*pow(y[0],4)); // dP/dm=-Gm/(4 pi r^4)
  p->x[0] = rho;

  if (!gsl_finite(dydx[0]) || !gsl_finite(dydx[1]))
  {
    if (verbose)
      cout<<"Warning: Derivative not finite. Can't find the density for "<<Phase->getEOS()<<" when performing constant temperature integral at location: r="<<y[0]<<"cm, m="<<x/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, Component masses "<<MC<<' '<<MM<<' '<<MW<<endl;
    return GSL_EBADFUNC;
  }

  return GSL_SUCCESS;
}

int derivs_m2(double x, const double y[], double dydx[], void * params)
// Three derivative equations used to integrate outside-in against mass.
// x is total mass minus the enclosed mass, y0 is radius, y1 is P, y2 is T.
// When integrate the function from the outer boundary, the enclosed mass is very large.  If the initial step size has to be very small due to the fluctuation of the density and the density of gas is very small, the mass change per step can be smaller than the precision of a double number.  It can crash the integration process.
// Using mass as variable because we need to stop by each density discontinuities, which is set by the mass.
{
  struct mass_params *p = (struct mass_params *) params;
  double rho = p -> x[0];
  int thermal = (int)p -> x[1];
  vector<double> M = p -> M;
  EOS *Phase = p -> Phase;
  double Mtot = accumulate(M.begin(), M.end(), 0.0) * ME;

  rho = Phase -> density(y[1], y[2], rho);
  if (!gsl_finite(rho)) 
  {
    if (verbose)
    {
      cout<<"Warning: Can't find density for "<<Phase->getEOS()<<" when performing outside-in integral against m at location: r="<<y[0]<<"cm, m="<<(Mtot-x)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<"K, Component masses ";
      for (int i=0; i < int(p -> M.size()); i++)
	cout<<p -> M[i]<<' ';
      cout<<endl;
    }
    return GSL_EBADFUNC;
  }

  dydx[0] = - 1 / (4 * pi * sq(y[0]) * rho); // dr/d(Mtot-m) = - 1 / (4 pi r^2 rho)
  dydx[1] = G * (Mtot - x) / (4 * pi * pow(y[0], 4)); // dP/d(Mtot-m)= G m / (4 pi r^4)
  if (thermal == 0 || thermal == 4)
    dydx[2] = 0;
  else if (thermal ==2)
    dydx[2] = Phase->dTdP(y[1], y[2])*dydx[1];
  else
    dydx[2] = -Phase->dTdm(Mtot-x, y[0], rho, y[1], y[2]);
  
  p -> x[0] = rho;
  
  if (!gsl_finite(dydx[0]) || !gsl_finite(dydx[1]) || !gsl_finite(dydx[2]))
  {
    if (verbose)
    {
      cout<<"Warning: Derivative not finite. Can't find density for "<<Phase->getEOS()<<" when performing outside-in integral against m at location: r="<<y[0]<<"cm, m="<<(Mtot-x)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<"K, Component masses ";
      for (int i=0; i < int(p -> M.size()); i++)
	cout<<p -> M[i]<<' ';
      cout<<endl;
    }
    return GSL_EBADFUNC;
  }

  return GSL_SUCCESS;
}

int derivs_m3(double x, const double y[], double dydx[], void * params)
// Three derivative equations used to integrate inside-out against mass.
// x is the enclosed mass, y0 is radius, y1 is P, y2 is T.
// When integrate the function from the inner boundary, the radius is very small.  If using the Mtot - m enclosed instead, the mass change per step can be smaller than the precision of a double number.  It can crash the integration process.
// Using mass as variable because we need to stop by each density discontinuities, which is set by the mass.
{
  struct mass_params *p=(struct mass_params *) params;
  double rho= p -> x[0];
  int thermal = (int)p -> x[1];
  vector<double> M = p -> M;
  EOS *Phase = p -> Phase;
  rho = Phase -> density(y[1], y[2], rho);

  if (!gsl_finite(rho)) 
  {
    if (verbose)
    {
      cout<<"Warning: Can't find density for "<<Phase->getEOS()<<" when performing inside-out integral against m at location: r="<<y[0]<<" m="<<x/ME<<" P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
      for (int i=0; i < int(p -> M.size()); i++)
	cout<<p -> M[i]<<' ';
      cout<<endl;
    }
    return GSL_EBADFUNC;
  }

  dydx[0] = 1 / (4 * pi * sq(y[0]) * rho); // dr/dm = 1 / (4 pi r^2 rho)
  dydx[1] = - G * x / (4 * pi * pow(y[0], 4)); // dP/dm= - G m / (4 pi r^4)
  if (thermal == 0 || thermal == 4)
    dydx[2] = 0;
  else if (thermal == 2)
    dydx[2] = Phase->dTdP(y[1], y[2])*dydx[1];
  else
    dydx[2] = Phase->dTdm(x, y[0], rho, y[1], y[2]);

  p -> x[0] = rho;
  
  if (!gsl_finite(dydx[0]) || !gsl_finite(dydx[1]))
  {
    if (verbose)
    {
      cout<<"Warning: Derivative not finite. Can't find density for "<<Phase->getEOS()<<" when performing inside-out integral against m at location: r="<<y[0]<<"cm, m="<<x/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<"K, Component masses ";
      for (int i=0; i < int(p -> M.size()); i++)
	cout<<p -> M[i]<<' ';
      cout<<endl;
    }
    return GSL_EBADFUNC;
  }

  return GSL_SUCCESS;
}

hydro::hydro(double Rp, vector<PhaseDgm> &Comp_in, vector<double> Mass_Comp, vector<double> Tgap, double ode_tol, double P0, bool isothermal):ode_status(2)
// Given the test planet radius, a vector of components, their masses, and the temperature discontinuity between each gaps (Temperature at the interior boundary minus exterior boundary.  The last temperature in the list is the planet equilibrium temperature). Bool isothermal determines whether use isothermal temperature profile or self-consistent calculation. The function integrates hydrostatic equation outside in.  The integration end when either the target mass reached or radius approaches 0.
// Integration inside out need to fine tune the center pressure.  The center pressure is way larger than the surface pressure.  The required accuracy for the center pressure and the integration process is higher than double.
// Integration outside-in will cause the same problem for the center mass.  It's hard to shoot the exact mass at the center.
// Integration outside in is more straight forward for temperature because only the temperature at outer boundary is set.
// Because the density has discontinuity at the shift between layers, the integration is operated within each layers.  And it's natural to choose mass as variable.
{
  if (Comp_in.size() != Mass_Comp.size() || Mass_Comp.size() != Tgap.size())
  {
    cout<<"Error. The length of input vectors doesn't match."<<endl;
    return;
  }
  double rhot, Tt, mmax;
  double hmax = ME*ode_tol;
  EOS *Phase, *new_Phase;
  int thermal;
  M_Comp = Mass_Comp;
  for (unsigned int i=0; i<Comp_in.size(); i++) 
    Comp.push_back(Comp_in[i]); 

  int n_Comp = M_Comp.size();
  double Mtot = accumulate(M_Comp.begin(), M_Comp.end(), 0.0) * ME;
  fitM=0;

  int status=0;
  R_Comp.assign(n_Comp-1, -1);			// reset the radius

  int i = n_Comp-1;
  Tt = Tgap.back();
  
  for ( ; M_Comp[i]*ME<1 ; i--) //Massless layer, skipped
  {
    R_Comp[i-1] = Rp;
    Tt += Tgap[i-1];
  }
  
  rb.resize(1, Rp);
  P.resize(1, P0);
  M.resize(1, Mtot);
  T.resize(1, Tt);
  Phase = Comp[i].find_phase(P[0], T[0]);

  if (!Phase)
  {
    if (verbose)
      cout<<"Warning: Can't find density for "<<Phase->getEOS()<<" at outer boundary in the first round of mode 0 integration with initial radius Rp="<<Rp<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    
    return;
  }
  
  Phaselist.resize(1, Phase);
  rhot = Phase -> density(P[0], T[0], -1);

  if (!gsl_finite(rhot)) //Failed to find the density
  {
    if (verbose)
      cout<<"Warning: Can't find density for "<<Phase->getEOS()<<" at outer boundary in the first round of mode 0 integration with initial radius Rp="<<Rp<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    
    return;
  }

  rho.resize(1, rhot);

  struct mass_params params = {{rho[0], 0.0}, M_Comp, Phase};
  gsl_odeiv2_system sys = {derivs_m2, NULL, 3, &params};
  const gsl_odeiv2_step_type *ODE = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(ODE, 3);
  
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1E-15, ode_tol); // eps_abs, eps_rel.  Control the error below eps_abs + eps_rel * y.  The eps_abs is the number precision (e.g. float).  Only important when y=0
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(3);

  double y[3] = {rb[0], P[0], T[0]};
  double m = Mtot - M[0];
  double h = pow(rb[0], 4) * P[0] / (10 * G * M[0]); // With the mass of step size, the pressure increases by about 1%.
  
  for (; i >= 0; i--)
  {    
    Phase = Comp[i].find_phase(P[0], T[0]);
    params.Phase = Phase;
    thermal = isothermal ? 0 : Phase->getthermal();
    thermal = Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.

    params.x[1] = (double)thermal;
    
    if (!Phase)
    {
      if (verbose)
      {
	cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" in mode 0 first round integration at location: r="<<y[0]<<"cm, m="<<(Mtot-m)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<"K, Component masses ";
	for (int j=0; j < int(M_Comp.size()); j++)
	  cout<<M_Comp[j]<<", ";
	cout<<endl;
      }
      rb.clear();
      P.clear();
      M.clear();
      T.clear();
      rho.clear();
      Phaselist.clear();
	
      gsl_odeiv2_evolve_free (e);
      gsl_odeiv2_control_free (c);
      gsl_odeiv2_step_free (s);

      return;
    }
    rhot = Phase -> density(y[1], y[2], -1);
    params.x[0] = rhot;

    rb.insert(rb.begin(), y[0]);
    P.insert(P.begin(), y[1]);
    M.insert(M.begin(), Mtot - m);
    T.insert(T.begin(), y[2]);
    rho.insert(rho.begin(), rhot);
    Phaselist.insert(Phaselist.begin(), Phase);
    
    mmax = accumulate(M_Comp.begin()+i, M_Comp.end(), 0.0) * ME;

    while (m < mmax && y[1] < 1E15)		// integrate to center.  If the mass is larger than 0 at the center, the gravity will diverge at the center.  Therefore cut the integration at 1E5 GPa.
    {
      status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &m, mmax, &h, y);
      if (status != GSL_SUCCESS)
      {
	if (status == GSL_FAILURE && y[0] > 1E-4 * Rp)
	{
	  if (verbose)
	    cout<<"Warning: In the first round integration of mode 0, at radius "<<y[0]<<"cm, mass "<<(Mtot - m)/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed because the required step size reaches machine precision.  This may caused by the EOS diverges at high pressure and high temperature, e.g. Bezacier's ice EOS. Solve this by changing the phase setting to guarantee Bezacier ice EOS is not ice EOS used for the highest ice pressure.  It may also caused by the discontinuity in density.  Try to resolve the error by reducing the ODE integration control accuracy requirement (larger ode_eps_rel0)."<<endl;
	}
	else if (status == GSL_FAILURE)
	  // If the radius is very small when the pressure is about to diverge, the mass to small that in each step may exceed the precision of float number.  In this case, the integration can be considered as finished.
	  break;
	else if (status == GSL_EBADFUNC)
	{
	  if (verbose)
	    cout<<"Warning: In the first round integration of mode 0, at radius "<<y[0]<<"cm, mass "<<(Mtot - m)/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed to calculate the density or derivative. "<<endl;
	  if (h/ME>1E-13)
	  {
	    h/=2;
	    gsl_odeiv2_step_reset(s);
	    gsl_odeiv2_evolve_reset(e);
	    continue;
	  }
	}
	else
	  if (verbose)
	    cout<<"Warning: In the first round integration of mode 0, at radius "<<y[0]<<"cm, mass "<<(Mtot - m)/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. Error code: "<<gsl_strerror (status)<<' '<<status<<endl;

	rb.clear();
	P.clear();
	M.clear();
	T.clear();
	rho.clear();
	Phaselist.clear();
	
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return;
      }

      new_Phase = Comp[i].find_phase(y[1], y[2]);
      rhot = new_Phase -> density(P[0], y[2], rhot);

      if (!new_Phase)
      {
	if (verbose)
	{
	  cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" in mode 0 first round integration at location: r="<<y[0]<<" m="<<(Mtot-m)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
	  for (int j=0; j < int(M_Comp.size()); j++)
	    cout<<M_Comp[j]<<", ";
	  cout<<endl;
	}
	rb.clear();
	P.clear();
	M.clear();
	T.clear();
	rho.clear();
	Phaselist.clear();

	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return;
      }

      if (new_Phase != Phase)
	// cross phase boundary but doesn't enter new composition
      {
	if ((M[0]-Mtot+m)>hmax)	// Last step size too large
	{
	  h = M[0]-Mtot+m;

	  do
	  {
	    if (h>hmax || status == GSL_FAILURE)
	      h/=2;
	  
	    if (new_Phase != Phase || status == GSL_FAILURE)
	    {
	      gsl_odeiv2_step_reset(s);
	      gsl_odeiv2_evolve_reset(e);
	
	      params.x[0] = rho[0];
	      y[0] = rb[0];
	      y[1] = P[0];
	      y[2] = T[0];
	      m = Mtot - M[0];
	    }
	    else
	    {
	      rb.insert(rb.begin(), y[0]);
	      P.insert(P.begin(), y[1]);
	      M.insert(M.begin(), Mtot - m);
	      Phaselist.insert(Phaselist.begin(), new_Phase);
	      T.insert(T.begin(), y[2]);
	      rho.insert(rho.begin(), rhot);
	    }

	    if (h<hmax/10)
	    {
	      if (verbose)
	      {
		cout<<"Warning: Failed in a phase transition for component: "<<Comp[i].getname()<<" in mode 0 first round integration at location: r="<<y[0]<<" m="<<(Mtot-m)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
		for (int j=0; j < int(M_Comp.size()); j++)
		  cout<<M_Comp[j]<<", ";
		cout<<endl;
	      }
	      rb.clear();
	      P.clear();
	      M.clear();
	      T.clear();
	      rho.clear();
	      Phaselist.clear();
	      gsl_odeiv2_evolve_free (e);
	      gsl_odeiv2_control_free (c);
	      gsl_odeiv2_step_free (s);

	      return;
	    }
	  
	    status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &m, h, y);

	    if (status == GSL_SUCCESS)
	    {
	      new_Phase = Comp[i].find_phase(y[1], y[2]);
	      rhot = new_Phase -> density(y[1], y[2], rhot);
	    }
	  } while (h>hmax);
	}
	thermal = isothermal ? 0 : new_Phase->getthermal();
	thermal = new_Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.
	params.Phase = new_Phase;
	Phase = new_Phase;
	params.x[1] = (double)thermal;
      }
        
      rb.insert(rb.begin(), y[0]);
      P.insert(P.begin(), y[1]);
      M.insert(M.begin(), Mtot - m);
      T.insert(T.begin(), y[2]);
      rho.insert(rho.begin(), rhot);
      Phaselist.insert(Phaselist.begin(), new_Phase);
      params.x[0] = rhot;
      h = min(h, 2*pow(y[0],3)*rhot);
      //cout<<y[0]/RE<<"RE "<<m/ME<<"ME "<<(mmax-m)/ME<<"ME "<<y[1]/1E10<<"GPa "<<params.x[0]<<"g/cm^3 "<<y[2]<<"K "<<h/ME<<"ME "<<Phaselist[0]->getEOS()<<" thermal_type "<<Phaselist[0]->getthermal()<<" params thermal "<<params.x[1]<<endl;
    }

    if (i)
    {
      R_Comp[i-1] = y[0];		// Get the size of each composition layer
      y[2] += Tgap[i-1];
    }

    while (i>1 && M_Comp[i-1]*ME<1)  //Massless layer, skipped
    {
      i--;
      R_Comp[i-1] = y[0]; 
      y[2] += Tgap[i-1];
    } 

    if ( y[1] >= 1E15 || status == GSL_FAILURE || (i==1 && M_Comp[0]*ME<1))		// Reaches the pressure limit, or code diverges almost at the center.  Exit iteration no matter which component current at.
    {
      for(int j = i-2; j >= 0; j-- )
	R_Comp[j] = R_Comp[j+1];
      break;
    }
    
    if (i)
    {
      gsl_odeiv2_step_reset(s);
      gsl_odeiv2_evolve_reset(e);
      h = pow(rb[0], 4) * P[0] / (10 * G * M[0]);
    }
  }
  count_shoot++;
  count_step+=rb.size();
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
}


hydro::hydro(double Pc, double MCin, double MMin, double MWin, double P0, double Teq):ode_status(2)
// Given the center pressure, Core (iron) mass, mantle (Si) mass, water mass, integrate hydrostatic equation inside out against mass.  The integration end when either the target mass reached or pressure smaller than P0.
{
  // set up the first grid
  double rhot;
  EOS *Phase;

  Comp.clear();
  Comp.push_back(Fe);
  Comp.push_back(Si);
  Comp.push_back(water);
  Comp.push_back(atm);
  
  M_Comp = {MCin, MMin, MWin, 0};
  double Mtot = accumulate(M_Comp.begin(), M_Comp.end(), 0.0) * ME;
    
  fitM = 0;

  rb.resize(1,0);
  P.resize(1,Pc);
  M.resize(1,0);
  T.resize(1,Teq);

  Phase = find_phase(10,Comp, M_Comp ,Pc,Teq);

  if (!Phase)
  {
    if (verbose)
      cout<<"Warning: Can't find the phase at center with mode 1 with Pc="<<Pc<<endl;

    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    return;
  }
  Phaselist.resize(1,Phase);
  rho.resize(1,Phase->density(Pc,Teq,-1));

  if (!gsl_finite(rho[0])) //Failed to find the density
  {
    if (verbose)
      cout<<"Warning: Can't find density for "<<Phase->getEOS()<<" at center with mode 1 with Pc="<<Pc<<endl;

    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    
    return;
  }
  
  // The ODE is singular at r=0, assume a small constant density ball. Using the Tylor expansion to get the second grid
  double R_T=1E4;		// The radius of center uniform density ball is 100 m.
  rb.push_back(R_T);
  P.push_back(Pc-2*pi*G*sq(rho[0]*R_T)/3.); // Delta_P = 2 pi G rho^2 Delta_r^2 / 3.
  M.push_back(4*pi*pow(R_T,3)*rho[0]/3.);
  T.push_back(Teq);
  Phase=find_phase(M[1], Comp, M_Comp, P[1], Teq);

  if (!Phase)
  {
    if (verbose)
      cout<<"Warning: Can't find the phase after first step with mode 1 with Pc="<<Pc<<endl;
    
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();

    return;
  }
  Phaselist.push_back(Phase);

  rho.push_back(Phase->density(P[1],Teq,rho[0]));

  if (!gsl_finite(rho.back())) //Failed to find the density
  {
    if (verbose)
      cout<<"Warning: Can't find density for "<<Phase->getEOS()<<" after first step with mode 1 with Pc="<<Pc<<endl;
    
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();

    return;
  }

  rhot=rho[1];
  struct double_params params = {{rhot, MCin, MMin, MWin, Teq}};

  gsl_odeiv2_system sys = {derivs_m, NULL, 2, &params};
  const gsl_odeiv2_step_type *ODE = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(ODE,2);
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-15, ode_eps_rel2); // eps_abs, eps_rel.  Control the error below eps_abs + eps_rel * y.  The eps_abs is the number precision (e.g. float).  Only important when y=0
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(2);

  double y[2] = {R_T, P[1]};
  double h=0.25*M[1], m=M[1];
  int status;
  R_Comp.assign(3, -1);			// reset the radius

  while (y[1] > P0 && m<Mtot)		// integrate to P0.
  {
    status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &m, Mtot, &h, y);
    if (status != GSL_SUCCESS)
    {
      if (status == GSL_FAILURE)
      {
	if (verbose)
	  cout<<"Warning: In the mode 1 integration, at radius "<<y[0]<<", mass "<<m/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed because the required step size reaches machine precision. This may caused by the EOS diverges at high pressure and high temperature, e.g. Bezacier's ice EOS. Solve this by changing the phase setting to guarantee Bezacier ice EOS is not ice EOS used for the highest ice pressure.  It may also caused by the discontinuity in density.  Try to resolve the error by reducing the ODE integration control accuracy requirement (larger ode_eps_rel2)."<<endl;
      }
      else if (status == GSL_EBADFUNC)
      {
	if (verbose)
	  cout<<"Warning: In the mode 1 integration, at radius "<<y[0]<<", mass "<<m/ME<<", pressure "<<y[1]/1E10<<"GPa, no solution can be found.  The final step size is "<<h/ME<< " MEarth. The solver failed to calculate the density or derivative."<<endl;
	if (h/ME>1E-13)
	{
	  h/=2;
	  gsl_odeiv2_step_reset(s);
	  gsl_odeiv2_evolve_reset(e);
	  continue;
	}
      }
      else
      {
	if (verbose)
	  cout<<"Warning: In the mode 1 integration, at radius "<<y[0]<<", mass "<<m/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, no solution can be found.  The final step size is "<<h/ME<<" MEarth. Error code: "<<gsl_strerror (status)<<' '<<status<<endl;
      }

      rb.clear();
      P.clear();
      M.clear();
      T.clear();
      rho.clear();
      Phaselist.clear();
      
      gsl_odeiv2_evolve_free (e);
      gsl_odeiv2_control_free (c);
      gsl_odeiv2_step_free (s);

      return;
    }
    
    rb.push_back(y[0]);
    P.push_back(y[1]);
    M.push_back(m);
    T.push_back(Teq);
    rho.push_back(params.x[0]);
    Phaselist.push_back(find_phase(m,Comp,M_Comp,y[1],Teq));
    
    if(R_Comp[0]<0 && m >= MCin*ME)	// Get the size of each composition layer
      R_Comp[0] = cbrt( pow(y[0],3.) - 3*(m-MCin*ME)/params.x[0]/4/pi );
    if(R_Comp[1]<0 && m >= (MCin+MMin)*ME)
      R_Comp[1] = cbrt( pow(y[0],3.) - 3*(m-(MCin+MMin)*ME)/params.x[0]/4/pi );
    if(R_Comp[2]<0 && m >= (MCin+MMin+MWin)*ME)
      R_Comp[2] = cbrt( pow(y[0],3.) - 3*(m-(MCin+MMin+MWin)*ME)/params.x[0]/4/pi );
  }

  if(R_Comp[0]<0)
    R_Comp[0] = rb.back();
  if(R_Comp[1]<0)
    R_Comp[1] = rb.back();
  if(R_Comp[2]<0)
    R_Comp[2] = rb.back();

  count_shoot++;
  count_step+=rb.size();

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
}


hydro::hydro(double Rp, double Pc, double Tc,  vector<PhaseDgm> &Comp_in, vector<double> Mass_Comp, vector<double> Tgap, bool isothermal, double Mfit, double ode_tol, double P0, double &Ri, double &Pi, double &Ti, double &Ro, double &Po, double &To):ode_status(2)
// After finishing a round of iteration of outside-in integration and receive a reasonable value of Rp, Pc, and Tc, using these values as the initial guess to integrate equations from both the surface, and the center, and meet two branches at the fitting mass Mfit (in g).  ode_tol is the ode relative tolerance 
// The R, P, and T value at inner edge (from inside-out integration with subscript i) and outer edge (from outside-in integration with subscript o) at the connection point is passed out. 
{
  fit_iter++;
  if (Comp_in.size() != Mass_Comp.size() || Mass_Comp.size() != Tgap.size())
  {
    if (verbose)
      cout<<"Warning. The length of input vectors doesn't match."<<endl;
    return;
  }

  double rhot, Tt, Tcor=0;
  //double Pt, mt;
  double hmax = ME*ode_tol;
  EOS *Phase, *new_Phase;
  int thermal;
  M_Comp = Mass_Comp;
  for (unsigned int i=0; i<Comp_in.size(); i++) 
    Comp.push_back(Comp_in[i]); 
  int n_Comp = M_Comp.size();
  double Mtot = accumulate(M_Comp.begin(), M_Comp.end(), 0.0) * ME;
  fitM = Mfit;
  
  int status=0;
  R_Comp.assign(n_Comp-1, -1);			// reset the radius

  // ---------------------------------
  // Integrate the inner part inside out

  int i = 0;
  Tt = Tc;

  for ( ; M_Comp[i]*ME<1 ; i++) //Massless layer at the center, skipped.  T should not be changed according to Tgap because Tc doesn't include consider Tgap of massless layer at the center.
  {
    R_Comp[i] = 0;
  }
  
  rb.resize(1, 0);
  P.resize(1, Pc);
  M.resize(1, 0);
  T.resize(1, Tt);
  Phase = Comp[i].find_phase(Pc, T[0]);
  if (!Phase)
  {
    if (verbose)
      cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" at the center in the second round of mode 0 integration with initial Pc="<<Pc/1E10<<"GPa, Tc="<<Tc<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    return;
  }
  Phaselist.resize(1, Phase);
  rhot = Phase -> density(P[0], T[0], -1);
  if (!gsl_finite(rhot)) //Failed to find the density
  {
    if (verbose)
      cout<<"Warning: Can't find the density for "<<Phase->getEOS()<<" at the center in the second round of mode 0 integration with initial Pc="<<Pc/1E10<<"GPa, Tc="<<Tc<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    return;
  }
  rho.resize(1, rhot);
  thermal = isothermal ? 0 : Phase->getthermal();
  thermal = Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.
  
// The ODE is singular at r=0, assume a small constant density ball. Using the Tylor expansion to get the second grid
  double M1 = min(1E14, 0.9*M_Comp[i]*ME);		// The mass of the center uniform density ball.  About a radius of 100m. 
  M.push_back(M1);
  rb.push_back(cbrt(3*M1/(4*pi*rho[0])));
  P.push_back(Pc - 2*pi*G*sq(rho[0]*rb[1])/3.); // Delta_P = 2 pi G rho^2 Delta_r^2 / 3.
  Phaselist.push_back(Phase);			// Assuming no phase changes in the center ~100 m size ball.

  if (thermal == 2)
  {
    Tt -= Phase->dTdP(P[0], Tc) * (P[0]-P[1]);
    rhot = Phase -> density(P[1], Tt, rhot);
  }
  else if (thermal == 0 || thermal == 4) // isothermal or fitted along isentrope
    rhot = Phase -> density(P[1], Tt, rhot);
  else
    rhot = Phase -> density(P[0], T[0], rho[0], P[1], Tt);

  if (!gsl_finite(rhot) || !gsl_finite(Tt)) //Failed to find the density or temperature
  {
    if (verbose)
      cout<<"Warning: Can't find the density or temperature after the first step in the second round of mode 0 integration with initial Pc="<<Pc/1E10<<"GPa, Tc="<<Tc<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    return;
  }

  rho.push_back(rhot);
  T.push_back(Tt);
  
  struct mass_params params = {{rhot, 0.0}, M_Comp, Phase};
  gsl_odeiv2_system sys = {derivs_m3, NULL, 3, &params};
  const gsl_odeiv2_step_type *ODE = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(ODE, 3);

  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-15, ode_tol); // eps_abs, eps_rel.  Control the error below eps_abs + eps_rel * y.  The eps_abs is the number precision (e.g. float).  Only important when y=0
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(3);

  double y[3] = {rb[1], P[1], T[1]};
  double m = M[1], mmax;
  double h = 0.25 * M[1];

  for ( ; i<n_Comp; i++)
  {
    Phase = Comp[i].find_phase(y[1], y[2]);
    params.Phase = Phase;
    thermal = isothermal ? 0 : Phase->getthermal();
    thermal = Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.

    params.x[1] = (double)thermal;
    
    if (!Phase)
    {
      if (verbose)
      {
	cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" in mode 0 second round inside-out integration at location: r="<<y[0]<<" m="<<m/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
	for (int j=0; j < int(M_Comp.size()); j++)
	  cout<<M_Comp[j]<<", ";
	cout<<endl;
      }
      rb.clear();
      P.clear();
      M.clear();
      T.clear();
      rho.clear();
      Phaselist.clear();
      gsl_odeiv2_evolve_free (e);
      gsl_odeiv2_control_free (c);
      gsl_odeiv2_step_free (s);

      return;
    }

    rhot = Phase -> density(y[1], Tt, -1);
    params.x[0] = rhot;
    
    rb.push_back(y[0]);
    P.push_back(y[1]);
    M.push_back(m);
    T.push_back(y[2]);
    rho.push_back(rhot);
    Phaselist.push_back(Phase);

    mmax = min(Mfit, accumulate(M_Comp.begin(), M_Comp.begin()+i+1, 0.0) * ME);


    while (m < mmax  && y[1]>P0)		// integrate from center to Mfit
    {
      status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &m, mmax, &h, y);

      if (status != GSL_SUCCESS)
      {
	if (status == GSL_FAILURE)
	{
	  if (verbose)
	    cout<<"Warning: In the second round of mode 0 inside-out integration, at radius "<<y[0]<<", mass "<<m/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed because the required step size reaches machine precision. This may caused by the EOS diverges at high pressure and high temperature, e.g. Bezacier's ice EOS. Solve this by changing the phase setting to guarantee Bezacier ice EOS is not ice EOS used for the highest ice pressure.  It may also caused by the discontinuity in density.  Try to resolve the error by reducing the ODE integration control accuracy requirement (larger ode_eps_rel1)."<<endl;
	}
	else if (status == GSL_EBADFUNC)
	{
	  if (verbose)
	    cout<<"Warning: In the second round of mode 0 inside-out integration, at radius "<<y[0]<<", mass "<<m/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed to calculate the density or derivative. "<<endl;
	  if (h/ME>1E-13)
	  {
	    h/=2;
	    gsl_odeiv2_step_reset(s);
	    gsl_odeiv2_evolve_reset(e);
	    continue;
	  }
	}
	else
	{
	  if (verbose)
	    cout<<"Warning: In the second round of mode 0 inside-out integration, at radius "<<y[0]<<", mass "<<m/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. Error code: "<<gsl_strerror (status)<<endl;
	}
	
	rb.clear();
	P.clear();
	M.clear();
	T.clear();
	rho.clear();
	Phaselist.clear();
	
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return;
      }

      new_Phase = Comp[i].find_phase(y[1], y[2]);
      rhot = new_Phase -> density(y[1], y[2], rhot);
      
      if (!new_Phase)
      {
	if (verbose)
	{
	  cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" in mode 0 second round inside-out integration at location: r="<<y[0]<<" m="<<m/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
	  for (int j=0; j < int(M_Comp.size()); j++)
	    cout<<M_Comp[j]<<", ";
	  cout<<endl;
	}
	rb.clear();
	P.clear();
	M.clear();
	T.clear();
	rho.clear();
	Phaselist.clear();
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return;
      }
      
      if (new_Phase != Phase) 	// cross phase boundary but doesn't enter new composition
      {
	if (m-M.back()>hmax) // the last step size is two large.
	{
	  h = m-M.back();

	  do
	  {
	    if (h>hmax || status == GSL_FAILURE)
	      h/=2;
	  
	    if (new_Phase != Phase || status == GSL_FAILURE)
	    {
	      gsl_odeiv2_step_reset(s);
	      gsl_odeiv2_evolve_reset(e);
	
	      params.x[0] = rho.back();
	      y[0] = rb.back();
	      y[1] = P.back();
	      y[2] = T.back();
	      m = M.back();
	    }
	    else
	    {
	      rb.push_back(y[0]);
	      P.push_back(y[1]);
	      M.push_back(m);
	      Phaselist.push_back(new_Phase);
	      T.push_back(y[2]);
	      rho.push_back(rhot);
	    }

	    if (h<hmax/10)
	    {
	      if (verbose)
	      {
		cout<<"Warning: Failed in a phase transition for component: "<<Comp[i].getname()<<" in mode 0 second round inside-out integration at location: r="<<y[0]<<" m="<<m/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
		for (int j=0; j < int(M_Comp.size()); j++)
		  cout<<M_Comp[j]<<", ";
		cout<<endl;
	      }
	      rb.clear();
	      P.clear();
	      M.clear();
	      T.clear();
	      rho.clear();
	      Phaselist.clear();
	      gsl_odeiv2_evolve_free (e);
	      gsl_odeiv2_control_free (c);
	      gsl_odeiv2_step_free (s);

	      return;
	    }
	  
	    status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &m, h, y);
	  
	    if (status == GSL_SUCCESS)
	    {
	      new_Phase = Comp[i].find_phase(y[1], y[2]);
	      rhot = new_Phase -> density(y[1], y[2], rhot);
	    }
	  } while (h>hmax);
	}
	thermal = isothermal ? 0 : new_Phase->getthermal();
	thermal = new_Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.
	params.Phase = new_Phase;
	Phase = new_Phase;
	params.x[1] = (double)thermal;
      }
      
      rb.push_back(y[0]);
      P.push_back(y[1]);
      M.push_back(m);
      Phaselist.push_back(new_Phase);
      T.push_back(y[2]);
      rho.push_back(rhot);
      params.x[0] = rhot;
    }

    if (m > (1-1E-15) * Mfit || y[1] <= P0)	// reaches the fitting point or the pressure approaches 0 before the fitting mass, need correction for the fitting values.
      break;
    
    if (i < n_Comp-1)
    {
      R_Comp[i] = y[0];
      y[2] -= Tgap[i];
      if (y[2] < 0)
      {
	Tcor = Tgap[n_Comp-1]-y[2];
	y[2] = Tgap[n_Comp-1];
      }
    } 

    while (M_Comp[i+1]*ME<1 && i<n_Comp-2)  //Massless layer, skipped
    {
      i++;
      R_Comp[i] = y[0];
      y[2] -= Tgap[i];
      if (y[2] < 0)
	y[2] = Tgap[n_Comp-1];
    }

    if (i < n_Comp-1)
    {
      gsl_odeiv2_step_reset(s);
      gsl_odeiv2_evolve_reset(e);
      h = pow(rb.back(), 4) * P.back() / (1000 * G * M.back());
    }    
  }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  int ninner = rb.size();	// layers of the inner integration
  
  Ri = rb.back();
  Pi = P.back();
  Ti = T.back() - Tcor;

  if (m < 0.999*Mfit && Pi <= P0) // The pressure approaches P0 before the fitting mass, need correction for the fitting values.
  {
    if (verbose)
    {
      cout<<"Warning: Starting pressure "<<Pc/1E10<<" GPa is too small for the inner part during fitting"<<endl;
      for (int j=0; j < int(M_Comp.size()); j++)
  	cout<<M_Comp[j]<<", ";
      cout<<Rp<<' '<<Pc<<' '<<Tc<<' '<<Mfit<<' '<<ode_tol;
      cout<<endl;
    }
    
    Ri *= cbrt(Mfit/M.back());
    Pi -= (Mfit-M.back()) * G * M.back() / (4 * pi * pow(rb.back(), 4));
    Ti = (T[ninner-2]*(M.back()-Mfit) + Ti*(Mfit-M[ninner-2])) / (M.back() - M[ninner-2]);
  }
// ---------------------------------
  // Integrate the outer part outside-in

  rb.insert(rb.begin()+ninner, Rp);
  P.insert(P.begin()+ninner, P0);
  M.insert(M.begin()+ninner, Mtot);

  i = n_Comp-1;
  Tt = Tgap.back();

  for ( ; M_Comp[i]*ME<1 ; i--) //Massless layer, skipped
  {
    R_Comp[i-1] = Rp;
    Tt += Tgap[i-1];
  }
    
  T.insert(T.begin()+ninner, Tt);
  Phase = Comp[i].find_phase(P.back(), Tt);
  if (!Phase)
  {
    if (verbose)
      cout<<"Warning: Can't find density at the surface in the second round of mode 0 integration with initial Rp="<<Rp<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();

    return;
  }
  Phaselist.insert(Phaselist.begin()+ninner, Phase);
  rhot = Phase -> density(P.back(), T.back(), -1);
  if (!gsl_finite(rhot)) //Failed to find the density
  {
    if (verbose)
      cout<<"Warning: Can't find density at the surface in the second round of mode 0 integration with initial Rp="<<Rp<<endl;
    rb.clear();
    P.clear();
    M.clear();
    T.clear();
    rho.clear();
    Phaselist.clear();
    
    return;
  }
  
  rho.insert(rho.begin()+ninner, rhot);

  params = {{rho[ninner], 0.0}, M_Comp, Phase};
  sys = {derivs_m2, NULL, 3, &params};
  s = gsl_odeiv2_step_alloc(ODE, 3);
  c = gsl_odeiv2_control_y_new(1e-15, ode_tol); // eps_abs, eps_rel.  Control the error below eps_abs + eps_rel * y.  The eps_abs is the number precision (e.g. float).  Only important when y=0
  e = gsl_odeiv2_evolve_alloc(3);

  y[0] = rb[ninner];
  y[1] = P[ninner];
  y[2] = T[ninner];
  m = Mtot - M[ninner];
  h = y[1] * pow(rb.back(), 4) / (1000 * G * M.back());	    // The step size is a 10^-4 of the pressure doubling radius.

  for ( ; i >= 0; i--)
  {
    Phase = Comp[i].find_phase(y[1], y[2]);
    params.Phase = Phase;
    thermal = isothermal ? 0 : Phase->getthermal();
    thermal = Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.

    params.x[1] = (double)thermal;
    
    if (!Phase)
    {
      if (verbose)
      {
	cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" in mode 0 second round outside-in integration at location: r="<<y[0]<<" m="<<(Mtot-m)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
	for (int j=0; j < int(M_Comp.size()); j++)
	  cout<<M_Comp[j]<<", ";
	cout<<endl;
      }
      rb.clear();
      P.clear();
      M.clear();
      T.clear();
      rho.clear();
      Phaselist.clear();
      gsl_odeiv2_evolve_free (e);
      gsl_odeiv2_control_free (c);
      gsl_odeiv2_step_free (s);

      return;
    }

    rhot = Phase -> density(y[1], y[2], -1);
    params.x[0] = rhot;
    
    rb.insert(rb.begin()+ninner, y[0]);
    P.insert(P.begin()+ninner, y[1]);
    M.insert(M.begin()+ninner, Mtot - m);
    T.insert(T.begin()+ninner, y[2]);
    rho.insert(rho.begin()+ninner, rhot);
    Phaselist.insert(Phaselist.begin()+ninner,Phase);
    
    mmax = min(Mtot - Mfit, accumulate(M_Comp.begin()+i, M_Comp.end(), 0.0) * ME);

    while (m < mmax && y[1] < 1E15)		// integrate from surface to Mfit
    {
      status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &m, mmax, &h, y);

      if (status != GSL_SUCCESS)
      {
	if (status == GSL_FAILURE)
	{
	  if (verbose)
	    cout<<"Warning: In the second round of mode 0 outside-in integration, at radius "<<y[0]<<", mass "<<(Mtot - m)/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed because the required step size reaches machine precision. This may caused by the EOS diverges at high pressure and high temperature, e.g. Bezacier's ice EOS. Solve this by changing the phase setting to guarantee Bezacier ice EOS is not ice EOS used for the highest ice pressure.  It may also caused by the discontinuity in density. Try to resolve the error by reducing the ODE integration control accuracy requirement (larger ode_eps_rel1)."<<endl;
	}
	else if (status == GSL_EBADFUNC)
	{
	  if (verbose)
	    cout<<"Warning: In the second round of mode 0 outside-in integration, at radius "<<y[0]<<", mass "<<(Mtot - m)/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. The solver failed to calculate the density or derivative. "<<endl;
	  if (h/ME>1E-13)
	  {
	    h/=2;
	    gsl_odeiv2_step_reset(s);
	    gsl_odeiv2_evolve_reset(e);
	    continue;
	  }
	}
	else
	{
	  if (verbose)
	    cout<<"Warning: In the second round of mode 0 outside-in integration, at radius "<<y[0]<<", mass "<<(Mtot - m)/ME<<"MEarth, pressure "<<y[1]/1E10<<"GPa, temperature "<<y[2]<<" no solution can be found.  The final step size is "<<h/ME<<" MEarth. Error code: "<<gsl_strerror (status)<<endl;
	}

	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	rb.clear();
	P.clear();
	M.clear();
	T.clear();
	rho.clear();
	Phaselist.clear();
	
	return;
      }

      new_Phase = Comp[i].find_phase(y[1], y[2]);
      rhot = new_Phase -> density(y[1], y[2], rhot);

      if (!new_Phase)
      {
	if (verbose)
	{
	  cout<<"Warning: Can't find the phase for component: "<<Comp[i].getname()<<" in mode 0 second round outside-in integration at location: r="<<y[0]<<" m="<<(Mtot-m)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
	  for (int j=0; j < int(M_Comp.size()); j++)
	    cout<<M_Comp[j]<<", ";
	  cout<<endl;
	}
	rb.clear();
	P.clear();
	M.clear();
	T.clear();
	rho.clear();
	Phaselist.clear();
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return;
      }
      
      if (new_Phase != Phase) // cross phase boundary but doesn't enter new composition
      {
	if ((M[ninner]-Mtot+m)>hmax)	//the last step size is two large.
	{
	  h = M[ninner]-Mtot+m;
	
	  do
	  {
	    if (h>hmax || status == GSL_FAILURE)
	      h/=2;
	  
	    if (new_Phase != Phase || status == GSL_FAILURE)
	    {
	      gsl_odeiv2_step_reset(s);
	      gsl_odeiv2_evolve_reset(e);
	
	      params.x[0] = rho[ninner];
	      y[0] = rb[ninner];
	      y[1] = P[ninner];
	      y[2] = T[ninner];
	      m = Mtot - M[ninner];
	    }
	    else
	    {
	      rb.insert(rb.begin()+ninner, y[0]);
	      P.insert(P.begin()+ninner, y[1]);
	      M.insert(M.begin()+ninner, Mtot - m);
	      Phaselist.insert(Phaselist.begin()+ninner, new_Phase);
	      T.insert(T.begin()+ninner, y[2]);
	      rho.insert(rho.begin()+ninner, rhot);
	    }

	    if (h<hmax/10)
	    {
	      if (verbose)
	      {
		cout<<"Warning: Failed in a phase transition for component: "<<Comp[i].getname()<<" in mode 0 second round outside-in integration at location: r="<<y[0]<<" m="<<(Mtot-m)/ME<<"MEarth, P="<<y[1]/1E10<<"GPa, T="<<y[2]<<" Component masses ";
		for (int j=0; j < int(M_Comp.size()); j++)
		  cout<<M_Comp[j]<<", ";
		cout<<endl;
	      }
	      rb.clear();
	      P.clear();
	      M.clear();
	      T.clear();
	      rho.clear();
	      Phaselist.clear();
	      gsl_odeiv2_evolve_free (e);
	      gsl_odeiv2_control_free (c);
	      gsl_odeiv2_step_free (s);

	      return;
	    }
	  
	    status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &m, h, y);

	    if (status == GSL_SUCCESS)
	    {
	      new_Phase = Comp[i].find_phase(y[1], y[2]);
	      rhot = new_Phase -> density(y[1], y[2], rhot);
	    }
	  } while (h>hmax);
	}
	thermal = isothermal ? 0 : new_Phase->getthermal();
	thermal = new_Phase->getthermal()==4 ? 4 : thermal; //for type 4, the isentrope temperature profile is enforced to the phase with this option even the isothermal option is chosen for the planet solver.
	params.Phase = new_Phase;
	Phase = new_Phase;
	params.x[1] = (double)thermal;
      }

      rb.insert(rb.begin()+ninner, y[0]);
      P.insert(P.begin()+ninner, y[1]);
      M.insert(M.begin()+ninner, Mtot - m);
      Phaselist.insert(Phaselist.begin()+ninner, new_Phase);      
      T.insert(T.begin()+ninner, y[2]);
      rho.insert(rho.begin()+ninner, rhot);
      params.x[0] = rhot;
    }

    if (m > (1-1E-15) * (Mtot - Mfit) || y[1] >= 1E15)	// reaches the fit point
      break;
    
    if (i)
    {
      R_Comp[i-1] = y[0];
      y[2] += Tgap[i-1];
    }

    while (M_Comp[i-1]*ME<1 && i>1)  //Massless layer, skipped
    {
      i--;
      R_Comp[i-1] = y[0]; 
      y[2] += Tgap[i-1];
    } 

    if (i)
    {
      gsl_odeiv2_step_reset(s);
      gsl_odeiv2_evolve_reset(e);
      h = pow(rb[ninner], 4) * P[ninner] / (1000 * G * M[ninner]);
    }
  }

  for (int j=0; j < n_Comp-1; j++) // Just in case some weird things occur near the fitting point
    if (R_Comp[j] < 0)
      R_Comp[j] = Ri;

  Ro = rb[ninner];
  Po = P[ninner];
  To = T[ninner];

  if ( m < 0.999*(Mtot-Mfit) && y[1] > 9E14) // The radius approachs 0 before the fitting mass, need correction for the fitting values.
  {
    if (verbose)
    {
      cout<<"Warning: Starting radius "<<Rp<<" cm is too small for the outer part during fitting"<<endl;
      for (int j=0; j < int(M_Comp.size()); j++)
	cout<<M_Comp[j]<<", ";
      cout<<Rp<<' '<<Pc<<' '<<Tc<<' '<<Mfit<<' '<<ode_tol<<endl;
    }

    Ro = cbrt(pow(Ro, 3) - 3 * (Mtot-m-Mfit) / (4 * pi));

    Po = (P[ninner+1]*(M[ninner]-Mfit) + Po*(Mfit-M[ninner+1])) / (M[ninner]-M[ninner+1]);
    To = (T[ninner+1]*(M[ninner]-Mfit) + To*(Mfit-M[ninner+1])) / (M[ninner]-M[ninner+1]);
  }

  // only keep the inside value at the fitting points to avoid non-monotonic order caused by the fitting error.
  rb.erase(rb.begin()+ninner);
  P.erase(P.begin()+ninner);
  M.erase(M.begin()+ninner);
  T.erase(T.begin()+ninner);
  rho.erase(rho.begin()+ninner);
  Phaselist.erase(Phaselist.begin()+ninner);

  count_shoot++;
  count_step+=rb.size();
  
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  if (rb.size() != P.size() || P.size() != M.size() || M.size() != T.size() || T.size() != rho.size() || rho.size() != Phaselist.size())
    cout<<"array sizes don't match "<<rb.size()<<' '<<P.size()<<' '<<M.size()<<' '<<T.size()<<' '<<rho.size()<<' '<<Phaselist.size()<<endl;
}


void hydro::print(string outfile, bool debug)
{
  ofstream fout(outfile.c_str(), ofstream::app);
  if(!fout)
  {
    cout<<"Error: Failed to open file "<<outfile<<endl;
    return;
  }
  double rtemp=0, Mtemp=0;
  int j=0;
  EOS *newPhase = Phaselist[0], *Phase=NULL;

  fout<<"Index\t Radius (earth)"<<"\t "<<"P (GPa)"<<"\t "<<"M (earth)"<<"\t "<<"Density (g cm^-3)"<<"\t "<<"T (K)"<<"\t "<<"Phase"<<endl;
  
  for(int i=0;i<int(rb.size());i++)
  {
    if (i==0 || rb[i]-rtemp > 1E-4*rb.back() || M[i] - Mtemp > 1E-4 *M.back() || i==int(rb.size())-1 || newPhase!=Phase)
    {
      if (newPhase!=Phase && j!=i-1 && Phase!=NULL && !debug) // make sure the previous step of the phase transition also get printed. Debug mode should have already output this
      {
	fout<<std::setprecision(8)<<i<<"\t "<<rb[i-1] / RE<<"\t "<<P[i-1] / 1E10<<"\t "<<M[i-1] / ME<<"\t "<<rho[i-1]<<"\t "<<T[i-1]<<"\t "<<Phase -> getEOS()<<endl;
      }
      Phase = newPhase;
      if (debug)
	fout<<i<<"\t "<<std::setprecision(16)<<rb[i] / RE<<"\t "<<P[i] / 1E10<<"\t "<<M[i] / ME<<"\t "<<rho[i]<<"\t "<<T[i]<<"\t "<<Phase -> getEOS()<<endl;
      else
	fout<<std::setprecision(8)<<i<<"\t "<<rb[i] / RE<<"\t "<<P[i] / 1E10<<"\t "<<M[i] / ME<<"\t "<<rho[i]<<"\t "<<T[i]<<"\t "<<Phase -> getEOS()<<endl;

      j = i;
      rtemp = rb[i];
      Mtemp = M[i];
    }
    else if (debug)
      fout<<i<<"\t "<<std::setprecision(16)<<rb[i] / RE<<"\t "<<P[i] / 1E10<<"\t "<<M[i] / ME<<"\t "<<rho[i]<<"\t "<<T[i]<<"\t "<<Phase -> getEOS()<<endl;
    if (i!=int(rb.size())-1)
      newPhase = Phaselist[i+1];
  }
  fout<<endl;
  fout.close();
}

int hydro::getLayer_from_r(double r) const
// return the layer index (from 0, count from bottom) by given the radius in RE. rb(l)<=r*RE<rb(l+1)
{
  r *= RE;
  int l=0,t=rb.size();
  int m;

  if(r>rb.back())
    return rb.size();
  else if(r<rb[0])
    return -1;

  while(t>l+1)
  {
    m=(l+t)>>1;
    if(rb[m]<=r)
      l=m;
    else
      t=m;
  }
  
  return l;
}

int hydro::getLayer_from_m(double m) const
// return the layer index (from 0, count from bottom) by given the mass in MEarth. M(l)<=m*MEarth<M(l+1)
{
  m *= ME;
  int l=0,t=M.size();
  int mid;

  if(m>M.back())
    return M.size();
  else if(m<M[0])
    return -1;

  while(t>l+1)
  {
    mid=(l+t)>>1;
    if(M[mid]<=m)
      l=mid;
    else
      t=mid;
  }
  
  return l;
}

vector<double> hydro::getRs()
// return the radii of core, mantle, and water layer, and the total radius in the unit of earth radii.
{
  if (!rb.empty())
  {
    vector<double> Rs(R_Comp);
    double scale = 1/RE;
    for(int i=0; i < int(Rs.size()); i++)
      Rs[i] *= scale;
    Rs.push_back(rb.back()*scale);
    return Rs;
  }
  else
  {
    cout<<"Error: Invalid planet model. Can't get boundary radius for ";
    for (int i=0; i < int(M_Comp.size()); i++)
      cout<<M_Comp[i]<<", ";
    cout<<endl;
    return {};			// return empty vector
  }
}
  
vector<double> hydro::getTs()
// return the temperatures at the outer side of each component interfaces as well as planet surface 
{
  if (!rb.empty())
  {
    vector<double> Rs = getRs();
    vector<double> Ts;
    int layer;
    for(int i=0; i < int(Rs.size()); i++)
    {
      layer = getLayer_from_r(Rs[i]);
      Ts.push_back(getT(layer));
    }
    return Ts;
  }
  else
  {
    cout<<"Error: Invalid planet model. Can't get temperature for ";
    for (int i=0; i < int(M_Comp.size()); i++)
      cout<<M_Comp[i]<<", ";
    cout<<endl;
    return {};			// return empty vector
  }
}


double R_hydro(double Rp, void *params)
// given the planet test radius, output the radius where mass reaches 0.
// p->x[2] > 0 means find a function value.  p->x[2] < 0 means planet structure integration was failed.
{
  struct loop_params *p = (struct loop_params *) params;
  double P0 = p -> x[2];
  vector<double> M = p -> M;
  vector<double> Tgap = p -> Tgap;
  bool isothermal = p -> iso;
  

  double Mtot = accumulate(M.begin(), M.end(), 0.0) * ME;

  if (Rp < 0)
    return Rp - Mtot;

  hydro temp(Rp, p -> Comp, M, Tgap, ode_eps_rel0, P0, isothermal);

  if (temp.getsize() == 0)	// can't find a solution
  {
    p -> x[1] = -1;
    p -> x[0] = numeric_limits<double>::quiet_NaN();
    return numeric_limits<double>::quiet_NaN();
  }
  if (temp.getM(0) < 1)	// The integration reaches the target mass before r reaches 0
  {
    p -> x[0] = temp.getR(0);
    p -> x[1] = temp.getP(0);
    return p -> x[0];
  }
  else  			// the integration reaches the center before reaches the target mass.  The integration is stopped at 1E5 GPa to prevent pressure diverge.  Assuming the remaining matter fixed at density 1 g cm^-3.  After the correction, the value should always negative because no material will have a density less than 1 g cm^3 at 1E5 GPa.  The smaller the remaining mass is, the closer to zero the value gets. 
  {
    p -> x[0] = -cbrt( 3 * temp.getM(0) / (4 * pi * 1) - pow(temp.getR(0), 3));
    p -> x[1] = temp.getP(0);

    if (p -> x[0] >= 0)
    {
      if (verbose)
	cout<<"Warning: The density of "<<find_phase(temp.getM(0), p -> Comp, M, p->x[1], temp.getT(0))->getEOS()<<" at pressure "<<p->x[1]/1E10<<"GPa and temperature "<<temp.getT(0)<<" K is "<<temp.getrho(0)<<" g cm^-3.  Such small density at high pressure may cause the program to be less efficient."<<endl;
      p -> x[0] = -cbrt( 3 * temp.getM(0) / (4 * pi * temp.getrho(0)) - pow(temp.getR(0), 3));
    }
    return p -> x[0];
  }
}

void hydro::setstatus(int s)
{
  if (s>3 || s<0)
  {
    cout<<"Error: Incorrect planet model status "<<s<<endl;
    exit(1);
  }
  else
    ode_status = s;
}


hydro* Rloop(vector<PhaseDgm> &Comp, vector<double> M_Comp, vector<double> ave_rho, vector<double> Tgap, double P0, bool isothermal, double &Rp, double &Pc, double &Tc)
// First round of iteration.  Input: list of componsations, masses and typical densities of each layers, the temperature discontinuity between each gaps (the last temperature in the list is the planet equilibrium temperature), boolean isothermal determines whether use isothermal temperature profile or self-consistent calculation.  The function iterates planet radius in order to get zero mass at the center.  The suggested planet radius, central pressure and temperature are also returned.
{
  int n_Comp = int(Comp.size());
  
  if (int(M_Comp.size()) != n_Comp || int(ave_rho.size()) != n_Comp || int(Tgap.size()) != n_Comp)
  {
    cout<<"Error: The array lengthes provided to Rloop (first round integration in mode 0) do not match."<<endl;
    return NULL;
  }

  Rp = Pc = Tc = 0;

  double Mtot = accumulate(M_Comp.begin(), M_Comp.end(), 0.0);
  for (int i = 0; i < n_Comp; i++)
  {
    if (i==3)
      Rp += M_Comp[i]*ME/(ave_rho[i]+sqrt(M_Comp[i]*Mtot)/20.);
    else
      Rp += M_Comp[i]*ME/ave_rho[i];
  }

  Rp = cbrt(3*Rp / (4*pi));			// initial guess of planet radius.

  int status;
  int iter = 0, max_iter = 100;

  double R_lo = 0.8 * Rp, R_hi = 1.2 * Rp, R_rsd = Rp;
  double success = -1;		// success > 0 means find a function value.  success < 0 means planet structure integration was failed.
  struct loop_params params = {{R_rsd, success, P0}, Comp, M_Comp, Tgap, isothermal, NULL};

  const gsl_root_fsolver_type *EQN = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (EQN);
  gsl_function F;

  F.function  = &R_hydro;
  F.params = &params;

  status = gsl_root_fsolver_set (s, &F, R_lo, R_hi);

  if (status == GSL_EINVAL && R_hi > R_lo)// Both shoots completed successfully and root is not straddled by the endpoints.
  {
    if (params.x[0] <= 0)	// R_hydro(R_hi, &params) called later. params can still reach its result.
      // Initial upper radius was too small
    {
      do
      {
	R_lo = R_hi;
	R_hi *= 1.4;
	iter++;

	if (params.x[1]<0)
	  // planet structure integration was failed.
	{
	  cout<<"Error: Failed to conduct the first round outside-in integration for ";
	  for (int i=0; i < n_Comp; i++)
	    cout<<M_Comp[i]<<", ";
	  cout<<endl<<"Can't find upper limit for the planet radius."<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
	if (iter >= max_iter || R_hi > 100*RE)
	{
	  cout<<"Error: Can't find planet radius upper limit for the first round outside-in iteration."<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
      }while (R_hydro(R_hi, &params) <= 0);
    }
      
    else
      // Initial lower radius was too large
    {
      do
      {
	R_hi = R_lo;
	R_lo *= 0.7;
	iter++;

	if (params.x[1]<0)
	  // planet structure integration was failed.
	{
	  cout<<"Error: Failed to conduct the first round outside-in integration for ";
	  for (int i=0; i < n_Comp; i++)
	    cout<<M_Comp[i]<<", ";
	  cout<<endl<<"Can't find lower limit for the planet radius."<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
	if (iter >= max_iter)
	{
	  cout<<"Error: Can't find planet radius lower limit within limited steps for the first round outside-in iteration."<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
      }while(R_hydro(R_lo, &params) >= 0);
    }
    gsl_root_fsolver_set (s, &F, R_lo, R_hi); // update correct endpoints that straddle target radius
  }
  else if (status == GSL_EINVAL)
  {
    cout<<"Error: The lower radius bound "<<R_lo<<" is larger than upper radius bound "<<R_hi<<" when conducting the first round outside-in integration for ";
    for (int i=0; i < n_Comp; i++)
      cout<<M_Comp[i]<<", ";
    cout<<endl;
    gsl_root_fsolver_free (s);
    return NULL;
  }
  else if (status != GSL_SUCCESS || params.x[1] < 0)
  {
    cout<<"Error: Failed to set up the radius solver for the first round outside-in integration for ";
    for (int i=0; i < n_Comp; i++)
      cout<<M_Comp[i]<<", ";
    cout<<endl;
    gsl_root_fsolver_free (s);
    return NULL;
  }

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);

    if (params.x[1]<0)
      // planet structure integration was failed.
    {
      cout<<"Error: Failed to conduct the first round outside-in integration for ";
      for (int i=0; i < n_Comp; i++)
	cout<<M_Comp[i]<<", ";
      cout<<endl;
      gsl_root_fsolver_free (s);
      return NULL;
    }

    Rp = gsl_root_fsolver_root (s);
    R_lo = gsl_root_fsolver_x_lower (s);
    R_hi = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_interval (R_lo, R_hi, 1E-10, R_eps_rel); // test radius precision
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status == GSL_CONTINUE)
  {
    cout<<"Error: Can't find outside-in planet model within maximum iteration "<<max_iter<<endl;
    gsl_root_fsolver_free (s);
    return NULL;
  }

  gsl_root_fsolver_free (s);

  Rp = R_hi;
  hydro *temp=new hydro(R_hi, Comp, M_Comp, Tgap, ode_eps_rel0, P0, isothermal); // Need the branch of solution that does not diverge at the center to get the central pressure and temperature

  Pc = temp -> getP(0) * sqrt(1 + temp -> getR(0) / Rp); // Make a correction to the pressure
  Tc = temp -> getT(0);

  return temp;
}

int fitting_error(const gsl_vector *x, void *params, gsl_vector *f)
// return the fitting error of R, P, and T given by the integrations from both end at the fitting point.
{
  struct loop_params *p = (struct loop_params *) params;
  double Mfit = p -> x[0];
  double ode_tol = p -> x[1];
  double P0 = p -> x[2];
  vector<double> M = p -> M;
  vector<double> Tgap = p -> Tgap;
  bool isothermal = p -> iso;
  if (p->model != NULL)
    delete p->model;
  
  const double Rp = gsl_vector_get(x, 0);
  const double Pc = gsl_vector_get(x, 1);
  const double Tc = gsl_vector_get(x, 2);

  double Ri, Pi, Ti, Ro, Po, To;

  p->model = new hydro(Rp, Pc, Tc, p->Comp, M, Tgap, isothermal, Mfit, ode_tol, P0, Ri, Pi, Ti, Ro, Po, To);

  if (p->model->getsize() == 0)	// can't find a solution
    return GSL_FAILURE;

  gsl_vector_set (f, 0, (Ri-Ro) / (Ri+Ro)); // set the mass different at the fitting point and try to normalize the difference
  gsl_vector_set (f, 1, (Pi-Po) / (Pi+Po));
  gsl_vector_set (f, 2, (Ti-To) / (Ti+To));

  return GSL_SUCCESS;
}


hydro* fitting_method(vector<PhaseDgm> &Comp, vector<double> M_Comp, vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal)	// Using the fitting method to integrate the model from both the center and the surface, and meet at an intermediate mass.  The precedure requires an accurate initial condition, run a Rloop first to determine a good initial condition.
{
  int n_Comp = int(M_Comp.size());
  if (int(M_Comp.size()) != n_Comp || int(ave_rho.size()) != n_Comp || int(Tgap.size()) != n_Comp)
  {
    cout<<"Error: The array lengthes provided to two branch fitting (second round integration in mode 0) do not match."<<endl;
    return NULL;
  }
  else if (n_Comp > 100 && verbose)
    cout<<"Warning: "<<n_Comp<<" layers of compositions are too many (over 100) and may lead to error."<<endl;

  double Mtot = accumulate(M_Comp.begin(), M_Comp.end(), 0.0) * ME;
  double Mfit = 0.2 * Mtot;
  double Mtemp;
  
  for (int i=0; i<n_Comp; i++)
    // avoid Mfit locates near the boundary of two components
  {
    Mtemp = accumulate(M_Comp.begin(), M_Comp.begin()+i+1, 0.0) * ME;
    if (fabs(Mtemp-Mfit)/Mtot < 1E-4)
    {
      for (int j=i+1; j<n_Comp; j++)
      {
	if (M_Comp[j]*ME > 5E-4*Mtot)
	{
	  Mfit += accumulate(M_Comp.begin()+i+1,M_Comp.begin()+j,0)*ME + 2E-4*Mtot;
	  break;
	}
	if (j == n_Comp-1 && verbose)
	  cout<<"Warning: "<<n_Comp<<" layers of compositions are too many (over 100) and may lead to error."<<endl;
      }
      break;
    }
    else if (Mtemp > Mfit)
      break;
  }
  
  if (Mtot < 1E-14 * ME && verbose)
    cout<<"Warning: The total planet mass "<<Mtot/ME<<" Earth mass is too small.  The code may not handle correctly."<<endl;

  double Rp, Pc, Tc;
  count_shoot=0;
  count_step=0;
  hydro *temp = Rloop(Comp, M_Comp, ave_rho, Tgap, P0, isothermal, Rp, Pc, Tc);
  
  if (!temp)
    return NULL;

  delete temp;

  const gsl_multiroot_fsolver_type *EQNS;
  gsl_multiroot_fsolver *s;

  int status;
  size_t iter, max_iter=15;
  double ode_tol = ode_eps_rel1;
  fit_iter = 0;
  const size_t n = 3;
  struct loop_params params = {{Mfit, ode_tol, P0}, Comp, M_Comp, Tgap, isothermal, NULL};
  gsl_multiroot_function f = {&fitting_error, n, &params};
  gsl_vector *x;
  EQNS = gsl_multiroot_fsolver_hybrids;

  for (int m=0; m<4; m++)
  {
    iter = 0;
    params.x[1] = ode_tol;

    x = gsl_vector_alloc (n);
  
    gsl_vector_set (x, 0, Rp);
    gsl_vector_set (x, 1, Pc);
    gsl_vector_set (x, 2, Tc);

    s = gsl_multiroot_fsolver_alloc (EQNS, 3);
    status = gsl_multiroot_fsolver_set (s, &f, x);

    if (status != GSL_SUCCESS)   
    {
      cout<<"Error: Multiroot solver failed to set up. The initial guesses of Rp="<<Rp<<" Pc="<<Pc<<" Tc="<<Tc<<". Component masses ";
      for (int i=0; i < n_Comp; i++)
	cout<<M_Comp[i]<<", ";
      cout<<endl;
      cout<<"Error code: "<<gsl_strerror (status)<<' '<<status<<endl;
      gsl_multiroot_fsolver_free (s);
      gsl_vector_free (x);
      return NULL;
    }
      
    do
    {
      iter++;

      status = gsl_multiroot_fsolver_iterate (s);

      if (status != GSL_SUCCESS) 
      {
	Rp = gsl_vector_get (s->x, 0);
	Pc = gsl_vector_get (s->x, 1);
	Tc = gsl_vector_get (s->x, 2);

	if (status == GSL_ENOPROG || status == GSL_ENOPROGJ || status == GSL_EBADFUNC) //Iteration is not making progress towards solution.  Exit the current iteration, decrease the integral tolerance and try again.
	{
	  gsl_multiroot_fsolver_free (s);
	  gsl_vector_free (x);
	  break;	
	}
	else
	{
	  cout<<"Error: Multiroot solver failed to iterate. The solution failed at Rp="<<Rp<<" Pc="<<Pc/1E10<<"GPa, Tc="<<Tc<<". Component masses ";
	  for (int i=0; i < n_Comp; i++)
	    cout<<M_Comp[i]<<", ";
	  cout<<endl;
	  cout<<"Error code: "<<gsl_strerror (status)<<' '<<status<<endl;
	  gsl_multiroot_fsolver_free (s);
	  gsl_vector_free (x);
	  return NULL;
	}
      }
    
      status = gsl_multiroot_test_residual (s->f, fit_eps_rel);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    if (status == GSL_SUCCESS)
    {
      Rp = gsl_vector_get (s->x, 0);
      Pc = gsl_vector_get (s->x, 1);
      Tc = gsl_vector_get (s->x, 2);
  
      gsl_multiroot_fsolver_free (s);
      gsl_vector_free (x);

      if (params.model->checkdummy() != "")
      {
	cout<<"Warning: P-T of the solution pass through a dummy EOS "<<params.model->checkdummy()<<". Result is not reliable. Component masses ";
	for (int i=0; i < n_Comp; i++)
	  cout<<M_Comp[i]<<", ";
	cout<<endl;
	params.model->setstatus(1);
      }
      else
	params.model->setstatus(0);
      return params.model;
    }
    else if (status == GSL_CONTINUE && m == 3)
    {   
      cout<<"Error: Can't find planet model using two side fitting method within maximum iterations "<<max_iter<<". Try stricter R_eps_rel.  Stricter ode tolerance (ode_eps_rel1) and less strict fitting tolerance (fit_eps_rel) may also help. The solution failed at Rp="<<Rp<<" Pc="<<Pc<<" Tc="<<Tc<<". Component masses ";
      for (int i=0; i < n_Comp; i++)
	cout<<M_Comp[i]<<", ";
      cout<<endl;
      cout<<"Error code: "<<gsl_strerror (status)<<' '<<status<<endl;
      gsl_multiroot_fsolver_free (s);
      gsl_vector_free (x);

      params.model->setstatus(2);
      return params.model;
    }
    else if (status == GSL_ENOPROG || status == GSL_ENOPROGJ || status == GSL_EBADFUNC || status == GSL_CONTINUE) //Iteration is not making progress towards solution
    {
      ode_tol /= 3.;
    }
    else
    {
      Rp = gsl_vector_get (s->x, 0);
      Pc = gsl_vector_get (s->x, 1);
      Tc = gsl_vector_get (s->x, 2);

      cout<<"Error: Multiroot solver failed to iterate. The solution failed at Rp="<<Rp<<" Pc="<<Pc/1E10<<"GPa Tc="<<Tc<<". Component masses ";
      for (int i=0; i < n_Comp; i++)
	cout<<M_Comp[i]<<", ";
      cout<<endl;
      cout<<"Error code: "<<gsl_strerror (status)<<' '<<status<<endl;
      gsl_multiroot_fsolver_free (s);
      gsl_vector_free (x);
      return NULL;
    }
  }

  cout<<"Error: Multiroot solver failed to making progress towards solution. The most plausible cause is that the fitting error tolerance is too strict compared to the ode error tolerance.  Try stricter ode tolerance (ode_eps_rel1) or less strict fitting tolerance (fit_eps_rel).  Try stricter R_eps_rel may also help.  \n The solution failed at Rp="<<Rp<<" Pc="<<Pc/1E10<<"GPa Tc="<<Tc<<". Component masses ";
  for (int i=0; i < n_Comp; i++)
    cout<<M_Comp[i]<<", ";
  cout<<endl;
  
  return NULL;
}

double P_hydro(double Pc, void *params)
// given the center pressure, output the difference between P0 and surface pressure.
{
  struct double_params *p = (struct double_params *) params;
  double MC = p->x[0];
  double MM = p->x[1];
  double MW = p->x[2];
  double Mtot = (MC+MM+MW)*ME;
  double P0 = p->x[4];
  int nlayer;			// number of layers
  
  if (Pc < 0)
    return Pc*1E10;
  
  hydro temp(Pc,MC,MM,MW,P0);
  nlayer = temp.getsize();
  
  if (nlayer == 0)	// can't find a solution
  {
    p -> x[4] = -1;
    return 0;
  }
  if (temp.totalM() >= Mtot)	// The integration reaches the target mass before P0.  The center pressure is too large.
    p -> x[3] = temp.getP(nlayer-1) - P0;
  else  			// the integration reaches P0 with less mass. The center pressure is too small.
    p -> x[3] = temp.getP(nlayer-1) - P0 - (Mtot-temp.totalM())*G*temp.totalM()/(4*pi*pow(temp.totalR(),4));

  p -> x[5] = 1;
  return p -> x[3];
}

hydro* getmass(double MC, double MM, double MW, double P0)
// Given the mass of core, mantle, and water, iterate center pressure to shot for a correct water layer mass.
{
  double Mtot = MC+MM+MW;
  double Pc = Mtot*1E12;
  if (Mtot < 1E-14 && verbose)
    cout<<"Warning: The total planet mass "<<Mtot<<" Earth mass is too small.  The code may not handle correctly."<<endl;

  int status, status_P;
  int iter = 0, max_iter = 100;

  double P_lo = 0.3*Pc, P_hi = 3*Pc, P_rsd = Pc;
  double success = -1;		// success > 0 means find a function value.  success < 0 means planet structure integration was failed.
  struct double_params params = {{MC, MM, MW, P_rsd, P0, success}};

  const gsl_root_fsolver_type *EQN = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (EQN);
  gsl_function F;

  F.function  = &P_hydro;
  F.params = &params;

  status = gsl_root_fsolver_set (s, &F, P_lo, P_hi);
  if (status == GSL_EINVAL && P_hi > P_lo)// root is not straddled by the endpoints.
  {
    if (params.x[3] <= 0)	// P_hydro(P_hi, &params) called later. params can still reach its result.
      // Initial upper central pressure was too small
    {
      do
      {
        P_lo = P_hi;
	P_hi *= 3;
	iter++;

	if (params.x[5]<0)
	  // planet structure integration was failed.
	{
	  cout<<"Error: Failed to conduct the mode 1 inside-out integration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}

	if (iter >= max_iter)
	{
	  cout<<"Error: Can't find central pressure upper limit within limited steps for mode 1 inside-out iteration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
      }while(P_hydro(P_hi,&params) <= 0);
    }
      
    else
      // Initial lower central pressure was too large
    {
      do
      {
	P_hi = P_lo;
	P_lo *= 0.3;
	iter++;

	if (params.x[5]<0)
	  // planet structure integration was failed.
	{
	  cout<<"Error: Failed to conduct the mode 1 inside-out integration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
	if (iter >= max_iter)
	{
	  cout<<"Error: Can't find central pressure lower limit within limited steps for mode 1 inside-out iteration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
	  gsl_root_fsolver_free (s);
	  return NULL;
	}
      }while(P_hydro(P_lo,&params) >= 0);
    }
    gsl_root_fsolver_set (s, &F, P_lo, P_hi); // update correct endpoints that straddle target central pressure.
  }
  else if (status == GSL_EINVAL)
  {
    cout<<"Error: The lower pressure bound "<<P_lo<<" is larger than upper pressure bound "<<P_hi<<" when conducting the mode 1 integration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
    gsl_root_fsolver_free (s);
    return NULL;
  }
  else if (status != GSL_SUCCESS)
  {
    cout<<"Error: Failed to set up the radius solver for the mode 1 integration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
    gsl_root_fsolver_free (s);
    return NULL;
  }
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);

    if (params.x[5]<0)
      // planet structure integration was failed.
    {
      cout<<"Error: Failed to conduct the mode 1 inside-out integration for (MC/MM/MW in ME) "<<MC<<' '<<MM<<' '<<MW<<'.'<<endl;
      gsl_root_fsolver_free (s);
      return NULL;
    }
    
    Pc = gsl_root_fsolver_root (s);
    P_lo = gsl_root_fsolver_x_lower (s);
    P_hi = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_interval (P_lo, P_hi, 1, P_eps_rel); // test pressure precision
    
    P_rsd = params.x[3];
    status_P = gsl_root_test_residual (P_rsd, P0*0.02); // test surface pressure precision
  }
  while (status == GSL_CONTINUE && status_P == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);
  
  if (status == GSL_CONTINUE && status_P == GSL_CONTINUE)
  {
    cout<<"Error: Can't find center pressure within maximum iterations "<<max_iter<<endl;
    return NULL;
  }

  hydro *temp=new hydro(Pc, MC, MM, MW, P0);
  
  if (abs(P_rsd) >= P0*0.02)
  {
    temp->setstatus(3);
    if(verbose)
      cout<<"Warning: The integral accuracy is not ideal for component masses "<<MC<<' '<<MM<<' '<<MW<<". A smaller ode_eps_rel2 and P_eps_rel are suggested."<<endl;
  }
  else
    temp->setstatus(0);
  return temp;
}

string hydro::checkdummy()
// Check if there is any Dummy EOS used in the profile
{
  for (int i=0; i<getsize(); i++)
    if (Phaselist[i]->getEOS().find("Dummy") != std::string::npos)
      return Phaselist[i]->getEOS();
  return "";
}  
