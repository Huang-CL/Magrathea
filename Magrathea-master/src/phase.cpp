#include "phase.h"
//Conditionals for Phase Diagrams in each layer begin LINE 249

PhaseDgm::PhaseDgm(string Comp_name, EOS* (*f)(double, double), int k, EOS** phase_name, double *start_pressure)
{
  Comp_type = Comp_name;
  phase_lowP = f;
  n = k;
  if(k > 0)
    phase_list = new EOS*[k];
  if(k > 1)
    start_P = new double[k-1];	// In GPa
  
  for(int i=0; i<n; i++)
  {
    phase_list[i] = phase_name[i];
    if(i != n-1)
      start_P[i] = start_pressure[i];
  }
}

PhaseDgm::PhaseDgm(const PhaseDgm &a)
{
  Comp_type = a.Comp_type;
  n = a.n;

  if(n > 0)
    phase_list = new EOS*[n];
  if(n > 1)
    start_P = new double[n-1];	// In GPa
  
  for(int i=0; i<n; i++)
  {
    phase_list[i] = a.phase_list[i];
    if(i != n-1)
      start_P[i] = a.start_P[i];
  }

  phase_lowP = a.phase_lowP;
}

PhaseDgm::~PhaseDgm()
{
  if(n>0)
    delete[] phase_list;
  if(n>1)
    delete[] start_P;
}

void PhaseDgm::set_phase_highP(int k, double *start_pressure, EOS** phase_name)
// start pressure is an array with dimension k-1.  The first phase will replace the one with the same name in the phase diagram.
{
  if(n > 1)
    delete[] start_P;
  if(n > 0)
    delete[] phase_list;
  n = k;
  phase_list = new EOS*[k];
  if(k > 1)
    start_P = new double[k-1];
  for(int i=0; i<n; i++)
  {
    phase_list[i] = phase_name[i];
    if(i != n-1)
      start_P[i] = start_pressure[i];
  }
}

EOS* PhaseDgm::find_phase(double P, double T)
// Pressure in microbar
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }
  EOS* phase_matched;
  string first_token, s;
  if(n == 0)
    return phase_lowP(P,T);
  else if(n == 1)
  {
    phase_matched = phase_lowP(P,T);
    s = phase_matched->getEOS();
    first_token = s.substr(0,s.find(" ("));
    if(first_token == phase_list[0]->getEOS())
      return phase_list[0];
    else
      return phase_matched;
  }
  else
  {
    int i = n-2;
    while(P < start_P[i]*1E10 && i > 0)
      i--;

    if(P >= start_P[i]*1E10)
      return phase_list[i+1];
    else
    {
      phase_matched = phase_lowP(P,T);
      s = phase_matched->getEOS();
      first_token = s.substr(0,s.find(" ("));
      if(first_token == phase_list[0]->getEOS())
        return phase_list[0];
      else
        return phase_matched;
    }
  }
}

double checkphase(double P, void *params)
{
  struct phase_params *p = (struct phase_params *) params;
  PhaseDgm* cmpn = p->cmpn;
  double Pu=p->Pu, Pl=p->Pl;
  
  if (P < Pl || P > Pu)
  {
    cout<<"Pressure "<<P/1E10<<" is outside the bound between "<<Pl<<" and "<<Pu<<" when solving the pressure of phase boundary for component "<<cmpn->getname()<<" from phase "<<p->Phase->getEOS()<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  
  double Tt;
  if (P == Pu)			// In case the output Tt is very close but not exact Tu or Tl.
    Tt = p->Tu;
  else if (P == Pl)
    Tt = p->Tl;
  else
    Tt = (p->Tl*(Pu-P)+p->Tu*(P-Pl))/(Pu-Pl);			// Linear interpolate temperature at P

  if (cmpn->find_phase(P, Tt) == p->Phase)
    return 1.;
  else
    return -1.;
}

EOS* PhaseDgm::find_phase_boundary(double Pl, double Pu, double Tl, double Tu, bool inward, double &Po, double &To, double &rhoo, double &Pn, double &Tn, double &rhon)
// Used when integrate adiabatic profile across the phase boundary.  Given the pressure in cgs, integration direction, the lower and upper pressure and temperature limit (at previous and next integral step depends on the integration direction).  Return the pressure, temperature, density at old (o) phase boundary and new (n) phase boundary.
{
  EOS* new_phase;
  double P1, P2;
  
  int iter = 0, max_iter = 100;
  int status;
  const gsl_root_fsolver_type *TPL = gsl_root_fsolver_bisection;

  gsl_root_fsolver *s = gsl_root_fsolver_alloc (TPL);
  gsl_function F;

  EOS* old_phase;
  if (inward)
    old_phase = find_phase (Pl, Tl);
  else
    old_phase = find_phase (Pu, Tu);
    
  struct phase_params params = {Pl, Pu, Tl, Tu, old_phase, this};
  F.function = &checkphase;
  F.params = &params;

  status = gsl_root_fsolver_set (s, &F, Pl, Pu);
  if (status == GSL_EINVAL)
  {
    cout<<"Error: The phase transition boundary is not between the lower pressure bound "<<Pl<<" and upper pressure bound "<<Pu<<", low temperature "<<Tl<<" and high temperature "<<Tu<<" when find the phase boundary of component "<<Comp_type<<" from phase "<<old_phase->getEOS()<<" to new phase ";

    if (inward)
      new_phase = find_phase (Pu, Tu);
    else
      new_phase = find_phase (Pl, Tl);
    
    gsl_root_fsolver_free (s);
    Po = numeric_limits<double>::quiet_NaN();
    To = numeric_limits<double>::quiet_NaN();
    rhoo = numeric_limits<double>::quiet_NaN();
    Pn = numeric_limits<double>::quiet_NaN();
    Tn = numeric_limits<double>::quiet_NaN();
    rhon = numeric_limits<double>::quiet_NaN();
    return NULL;
  }
  else if (status != GSL_SUCCESS)
  {
    cout<<"Error: Failed to set up the phase boundary of component "<<Comp_type<<" from phase "<<old_phase->getEOS()<<endl;
    gsl_root_fsolver_free (s);
    Po = numeric_limits<double>::quiet_NaN();
    To = numeric_limits<double>::quiet_NaN();
    rhoo = numeric_limits<double>::quiet_NaN();
    Pn = numeric_limits<double>::quiet_NaN();
    Tn = numeric_limits<double>::quiet_NaN();
    rhon = numeric_limits<double>::quiet_NaN();
    return NULL;
  }
  do
  {
    iter++;

    status = gsl_root_fsolver_iterate (s);

    P1 = gsl_root_fsolver_x_lower (s);
    P2 = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_interval (P1, P2, 1E-10, 1E-12);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status == GSL_CONTINUE)
  {
    cout<<"Error: Can't find the pressure of the phase boundary from "<<find_phase (Pl, Tl) -> getEOS()<<" and "<<find_phase (Pu, Tu) -> getEOS()<<" at pressure "<<Pl/1E10<<" GPa within maximum interation "<<max_iter<<endl;
    gsl_root_fsolver_free (s);
    Po = numeric_limits<double>::quiet_NaN();
    To = numeric_limits<double>::quiet_NaN();
    rhoo = numeric_limits<double>::quiet_NaN();
    Pn = numeric_limits<double>::quiet_NaN();
    Tn = numeric_limits<double>::quiet_NaN();
    rhon = numeric_limits<double>::quiet_NaN();
    return NULL;
  }

  if (inward)
  {
    Pn = P2;
    Po = P1;
  }
  else
  {
    Pn = P1;
    Po = P2;
  }

  Tn = (Tl*(Pu-Pn)+Tu*(Pn-Pl))/(Pu-Pl);			// Linear interpolate temperature at P
  To = (Tl*(Pu-Po)+Tu*(Po-Pl))/(Pu-Pl);			// Linear interpolate temperature at P
  new_phase = find_phase (Pn, Tn);
  rhon = new_phase->density (Pn, Tn, rhoo);
  rhoo = find_phase (Po, To) -> density(Po, To, rhoo);
  gsl_root_fsolver_free (s);
  return new_phase;
}

// Phase curve function in Dunaeva et al. 2010, input P in GPa, return T
double dunaeva_phase_curve(double P, double a, double b, double c, double d, double e)
{
  P *= 1E4;			// convert P from GPa to bar
  if (P <= 0)
  {
    cout<<"Error: P = "<<P<<" bar. Dunaeva phase curve can't handle negative pressure."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  return a + b*P + c*log(P) + d/P + e*sqrt(P);
}
 
// ========== Phase Diagram for Core  ================
// ---------------------------------
// Fe Default: hcp and Liquid iron
EOS* find_phase_Fe_default(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }

  P /= 1E10;			// convert microbar to GPa

  // Default Core
  if( T > 12.8*P + 2424 && T > 13.7*P + 2328)   // melting curve from Dorogokupets et al. 2017, Scientific Reports. fcc and hcp Fe melting curve.
    return Fe_liquid;
  else
    return Fe_hcp;             // use hcp Iron for all regions.
}

// ---------------------------------
// Fe fccbcc: Includes low pressure fcc and bcc iron
EOS* find_phase_Fe_fccbcc(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }

  P /= 1E10;			// convert microbar to GPa
	
  if(T>575+18.7*P+0.213*pow(P,2)-0.000817*pow(P,3) && P<98.5) // Anzellini et al. 2013
  {
    if(T<1991*pow((((P-5.2)/27.39)+1),1/2.38) && P<98.5)
      if(T<-41.1*P+1120)	
        return Fe_bcc;
      else
        return Fe_fcc;
    else
      return Fe_liquid;
  }
  else
  {
    if(T<3712*pow((((P-98.5)/161.2)+1),1/1.72))
      if(T<-61.2*P+1266)     // Dorogokupets et al. 2017
        return Fe_bcc;
      else	
        return Fe_hcp;
    else
      return Fe_liquid;
  }
}

// ========== Phase Diagram for Mantle  ================
// ---------------------------------
// Si Default: Upper Mantle: Fo, Wds, Rwd, and liquid ; Lower Mantle: Brg, PPv
EOS* find_phase_Si_default(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }

  P /= 1E10;			// convert microbar to GPa
  if(P > 112.5 + 7E-3*T)      // Phase transfer curve from Ono & Oganov 2005, Earth Planet. Sci. Lett. 236, 914
    return Si_PPv_Sakai;
  else if (T > 1830*pow(1+P/4.6, 0.33)) // Melting curve from Belonoshko et al. 2005 Eq. 2
    return Si_Liquid_Wolf;
  else if (P > 24.3+(-2.12E-4*T)+(-3.49E-7*pow(T, 2))) // Dorogokupets et al. 2015
    return Pv_Doro;
  else if (P > 8.69+6.05E-3*T)
    return Rwd;
  else if (P > 9.45+2.76E-3*T)
    return Wds;
  else
    return Fo;
}

// ---------------------------------
// Si Simplified: Brg, PPv, and liquid 
EOS* find_phase_Si_simple(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }

  P /= 1E10;			// convert microbar to GPa
  if(P > 112.5 + 7E-3*T)	// Phase transfer curve from Ono & Oganov 2005, Earth Planet. Sci. Lett. 236, 914
    return Si_PPv_Sakai;
  else if (T > 1830*pow(1+P/4.6, 0.33)) // Melting curve from Belonoshko et al. 2005 Eq. 2
    return Si_Liquid_Wolf;
  else
    return Si_Pv;
}

// ---------------------------------
// PREM tabulated mantle
EOS* find_phase_PREM(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }

  P /= 1E10;			// convert microbar to GPa	
  if(P > 3500)
    return Si_Seager;
  else if(P > 23.83)
    return Si_BM2fit;
  else
    return Si_PREM;
}
//-----------------------------------
// C Default
EOS* find_phase_C_simple(double P, double T)
{
  if (P <= 0 || T <= 0)
    return NULL;
  P /= 1E10;      // convert microbar to GPa
  
  if (P>=970.679+(-1.52854E-2*T)+(-5.72152E-7*pow(T,2))) // Transition from Benedict et al (2018)
    return BC8;
  else if (P>=1.949+(T+273)/400)  // Transition from Kennedy and Kennedy (1976)
    return Diam;
  else
    return Graph;
}
 
// ========== Phase Diagram for Hydrosphere  ================

// ---------------------------------
// H2O Water/Ice boundaries primarily form Dunaeva et al. 2010
EOS* find_phase_water_default(double P, double T)
// input P in cgs
{
  P /= 1E10;			// convert microbar to GPa
  double Tt1, Tt2;
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }
  if( P < 0.208566 )		// liquid water or Ice Ih (Dunaeva et al. 2010)
  {
    Tt1 = dunaeva_phase_curve(P, 273.0159, -0.0132, -0.1577, 0, 0.1516);
    if (!gsl_finite(Tt1))
    {
      cout<<"Error: Can't find the phase of H2O at P="<<P<<" GPa and T="<<T<<" K."<<endl;
      return NULL;
    }
    if(T > Tt1)
      return Water;
    else
      return IceIh;
  }

  else if( P < 0.6324 )		// liquid water or Ice II, III, V
  {
    Tt1 = dunaeva_phase_curve(P, 10.277, 0.0265, 50.1624, 0.5868, -4.3288);
    Tt2 = dunaeva_phase_curve(P, 5.0321, -0.0004, 30.9482, 1.0018, 0);
    if (!gsl_finite(Tt1) || !gsl_finite(Tt2))
    {
      cout<<"Error: Can't find the phase of H2O at P="<<P<<" GPa and T="<<T<<" K."<<endl;
      return NULL;
    }
    if( (P < 0.3501 && T > Tt1) || (P >= 0.3501 && T > Tt2))
      return Water;
    else
    {
      return Ice_Dummy;
    }
  }
  
  else if(P < 2.216)		// liquid water or Ice VI, VII
  {
    Tt1 = dunaeva_phase_curve(P, 4.2804, -0.0013, 21.8756, 1.0018, 1.0785);
    Tt2 = dunaeva_phase_curve(P, -47.8507, 0, -389.006, 0.9932, 28.8539);
    if (!gsl_finite(Tt1) || !gsl_finite(Tt2))
    {
      cout<<"Error: Can't find the phase of H2O at P="<<P<<" GPa and T="<<T<<" K."<<endl;
      return NULL;
    }
    if(T > Tt1)
      return Water;
    else if(T > Tt2)
      return IceVI_Bezacier;
      
    else
      return IceVII_Bezacier;
    //return IceVII_FFH2004;  // Examples of choosing other EOSs
    //return IceVII;
    //return IceZeng2013FFH;
  }

  else if(P < 5.10)		// liquid water or Ice VII.
  {
    if(T>1200)			// A dummy melting curve to avoid ice VII EOS extrapolated to temperature too high.
      return Water_sc_Mazevet;
    else
      return IceVII_Bezacier;
  }
  else if(P < 30.9)		// liquid water or Ice VII'.
  {
    if(T>1200)			// A dummy melting curve to avoid Ice VII EOS extrapolated to temperature too high.
      return Water_sc_Mazevet;
    else
      //return IceVIIp;         //Region of possible transitional Ice VII' reported in Grande not used in default
      return IceVII_Bezacier;
  }
  else if (P<700)				// Phase diagram becomes more uncertain over 30.9 GPa
  {						// use Ice X if T<2250 and P<700 otherwise supercritical similiar to AQUA 
    if(T<2250)
      return IceX;
    else
      return Water_sc_Mazevet;
  }
  else
    return Water_sc_Mazevet;
}


// ---------------------------------
// AQUA Haldemann et al. 2020 Tabulated Ice, Liquid, Vapor, Supercritical
EOS* find_phase_water_tabulated(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }

  P /= 1E10;			// convert microbar to GPa
  return H2O_AQUA;

  if (P < 62.5)
  {
    return H2O_AQUA;
  }
  else
    return Water_sc_Mazevet;

}

// ========== Phase Diagram for Atmosphere  ================
// ---------------------------------
// Ideal Gas: Isothermal for P<100 bar, the radiative/convective boundary (Nixon & Madhusudhan 2021). Adiabatic ideal gas for P > 100 bar
EOS* find_phase_gas_default(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }
  else if (P < 1E8)
    return Gas_iso;
  else
    //return vdW_H2;
    return Gas;
}

// ---------------------------------
// H/He Gas: Isothermal ideal for P<100 bar, P>100 bar: tabulated real gas, Chabrier & Debras 2021 Y=0.275
EOS* find_phase_HHe_tabulated(double P, double T)
{
  if (P <= 0 || T <= 0)
  {
    return NULL;
  }
  else if (P < 1E8)
    return Gas_iso;
  else
    return Gas_hhe;
}


PhaseDgm core("core", find_phase_Fe_default); //Phase Diagrams for Core
PhaseDgm core1("core1", find_phase_Fe_fccbcc); 
PhaseDgm mant("mantle", find_phase_Si_default); //Phase Diagram for Mantle
PhaseDgm mant1("mantle1", find_phase_Si_simple); 
PhaseDgm mant2("mantle2", find_phase_PREM); 
PhaseDgm mant3("mantle3", find_phase_C_simple); //Phase Diagram for Carbon Mantle
PhaseDgm water("water", find_phase_water_default); //Phase Diagram for Hydrosphere
PhaseDgm water1("water1", find_phase_water_tabulated);
PhaseDgm atm("atm", find_phase_gas_default); //Phase Diagram for Atmosphere
PhaseDgm atm1("atm1", find_phase_HHe_tabulated);

EOS* find_phase(double m, double MC, double MM, double MW, double MG, double P, double T, bool inward)
// given the accumulated mass (in g), P (in cgs) and T, return the corresponding phase
{
  if (inward)
  {
    if(m < MC*ME)			// iron core
      return core.find_phase(P,T);

    else if(m < (MM+MC)*ME)		// silicon mantle
      return mant.find_phase(P,T);

    else if(m < (MM+MC+MW)*ME)		// hydrosphere
      return water.find_phase(P,T);

    else				// return the outermost existing layer 
    {
      if(MG > 1E-18)
        return atm.find_phase(P,T);
      else if(MW > 1E-18)
        return water.find_phase(P,T);
      else if(MM > 1E-18)
        return mant.find_phase(P,T);
      else
        return core.find_phase(P,T);
    }
  }
  else
  {
    if(m <= MC*ME)			// iron core
      return core.find_phase(P,T);

    else if(m <= (MM+MC)*ME)		// silicon mantle
      return mant.find_phase(P,T);

    else if(m <= (MM+MC+MW)*ME)		// hydrosphere
      return water.find_phase(P,T);

    else				// return the outermost existing layer 
    {
      if(MG > 1E-18)
        return atm.find_phase(P,T);
      else if(MW > 1E-18)
        return water.find_phase(P,T);
      else if(MM > 1E-18)
        return mant.find_phase(P,T);
      else
        return core.find_phase(P,T);
    }
  }
}  


EOS* find_phase(double m, vector<PhaseDgm> &Comp, vector<double> M, double P, double T, bool inward)
{
  if (Comp.size() != M.size())
  {
    cout<<"Error. The length of input vectors doesn't match."<<endl;
    exit(1);
  }

  double mass_layers = 0;
  int n_comp = Comp.size();

  for (int i=0; i<n_comp; i++)
  {
    if (M[i]*ME<1)
      continue;
    mass_layers += M[i];
    if (inward)
    {
      if (m < mass_layers*ME)
        return Comp[i].find_phase(P, T);
    }
    else
      if (m <= mass_layers*ME)
        return Comp[i].find_phase(P, T);
  }
  // if nothing matched, return the outermost existing layer
  return Comp.back().find_phase(P, T);
}
