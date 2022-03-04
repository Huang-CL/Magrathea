#include "EOSlist.h"

/*
  phasetype is the name of a phase. The comment about the EOS used for the phase should be in the parentheses separated by a space.
0.	eqntype. 8-12 for RTpress style.
1.	V0 in cm^3 / mol.
	For ice, 1 \AA^3 = N_A / (2*10^24) cm^3/mol = 0.3011 cm^3/mol
2.	K0 in GPa
3.	K0p
4.	K0pp in GPa ^-1
5.	mmol in g / mol, or mean molecular weight of gas, or in g / mol for RTpress style
6.	P0 (GPa) the minimum pressure.  The pressure correspond to V0
7.	cp in 10^7 erg/g/K Heat capacity per mass at constant pressure
8.	Theta0 (K), a fitting parameter to Einstein temperature or Debye temperature
9.	gamma0, a fitting parameter of Grueneisen parameter
10.	beta, a fitting parameter of Grueneisen parameter.  In RTpress style, it represents the "m" which stands for the power-law exponent in the thermal deviation term.  Theoretically expected value: 0.6.
11.	gammainf, a fitting parameter of Grueneisen parameter
12.	gamma0p, volume derivative of the Grueneisen parameter
13.	e0 (10^-6 K^-1), electronic contribution to Helmholtz free energy
14.	g, is an electronic analogue of the Grueneisen parameter
15.	n is the number of atoms in the chemical formula of the compound.  Should have n*NA atoms within V.  The n of ideal gas is the number of atoms per molecule for the purpose of adiabatic index.  NOTE: n=2 for collinear molecules e.g. carbon dioxide!  Isothermal atmosphere can be achieved by setting n=0.
16.     Z is the atomic number (number of electron)
17.	T0, the reference temperature for the thermal pressure
18.	alpha, the coefficient of thermal expansion at a reference pressure P0 in 10^-6 K^-1
19.Debye_approx, whether use Debye approximation or Einstein approximation. Debye approximation is slower but more accurate at temperature lower than Debye/Einstein temperature.  Positive number for Debye, otherwise Einstein.
20.     thermal_type, indicates the thermal type of the phase.  0 indicates no temperature profile available, 1 indicates entropy method, 2 indicates the temperature gradient method.  The only method to set the gradient is using the modify_extern_dTdP function, 3 indicates ideal gas, 4 indicates the EOS is fitted along the isentrope, type 8 indicates RTpress style.
21-28.  at1-at4 & ap1 - ap4


For RTpress style of EOS, also need a _b array. They are fitted polynomial parameters of the thermal coefficients b(V) in erg/mol.  Convert eV/atom to erg/mol need to multiply eV_erg*n*NA. For example, for MgSiO3, 0.9821 eV/atom = 4.824E12 *0.9821 erg/mol = 4.738E12 erg/mol.
*/



// ==========  Iron ================

// ---------------------------------
// Liquid Iron,  Dorogokupets et al. 2017, Scientific Reports.

double Fe_liquid_array[][2] = {{0,2}, {1,7.957}, {2,83.7}, {3,5.97}, {5,mFe}, {6, 1E-4}, {8, 263}, {9, 2.033}, {10, 1.168}, {11,0}, {13,198}, {14, 0.884}, {15, 1}, {16,26}, {17, 1811}};

EOS *Fe_liquid = new EOS("Fe liquid (Dorogokupets)", Fe_liquid_array, sizeof(Fe_liquid_array)/2/sizeof(Fe_liquid_array[0][0]));

// -----------------------------------
// Liquid Iron 2, Anderson & Ahrens, 1994 JGR

double Fe_liquid2_array[][2] = {{0,1}, {1,7.95626}, {2,109.7}, {3,4.66}, {4,-0.043}, {5,mFe}, {15,1}, {16,26}};

EOS *Fe_liquid2 = new EOS("Fe liquid (Anderson)", Fe_liquid2_array, sizeof(Fe_liquid2_array)/2/sizeof(Fe_liquid2_array[0][0]));


// -----------------------------------
// Iron hcp, Bouchet et al. 2013, PRB 87, 094102

double Fe_hcp_array[][2] = {{0,3}, {1,6.29}, {2,253.844}, {3,4.719}, {5,mFe}, {8,44.574}, {9,1.408}, {10,0.826}, {11,0.827}, {13,212.1}, {14,1.891}, {15,1}, {16,26}};

EOS *Fe_hcp = new EOS("Fe hcp (Bouchet)", Fe_hcp_array, sizeof(Fe_hcp_array)/2/sizeof(Fe_hcp_array[0][0]));

// -----------------------------------
// Iron hcp, Dorogokupets et al. 2017, Scientific Reports.

double Fe_hcp2_array[][2] = {{0,2}, {1,6.8175}, {2,148.0}, {3,5.86}, {5,mFe}, {8, 227}, {9, 2.2}, {10, 0.01}, {11,0}, {13,126}, {14,-0.83}, {15,1}, {16,26}, {17, 298.15}};

EOS *Fe_hcp2 = new EOS("Fe hcp (Dorogokupets)", Fe_hcp2_array, sizeof(Fe_hcp2_array)/2/sizeof(Fe_hcp2_array[0][0]));

// -----------------------------------
// Iron hcp, Smith et al. 2018, Nature Astronomy. <4.0 Earth-mass suggested (Gruneisen determined from Fig. 3b)

double Fe_hcp3_array[][2] = {{0,2}, {1,mFe/8.43}, {2,177.7}, {3,5.64}, {5,mFe}, {8,322}, {9,2.09}, {10,1.01}, {11,0.0500}, {15,1}, {16,26}};

EOS *Fe_hcp3 = new EOS("Fe hcp (Smith)", Fe_hcp3_array, sizeof(Fe_hcp3_array)/2/sizeof(Fe_hcp3_array[0][0]));

// -----------------------------------
// Iron Dummy. Used to fill in phase space that no EOS provided.

EOS *Fe_Dummy = new EOS("Fe Dummy", Fe_hcp3_array, sizeof(Fe_hcp3_array)/2/sizeof(Fe_hcp3_array[0][0]));

// -----------------------------------
// Fe-7Si shock, Wicks et al. 2018, Science Advances.

double Fe_7Si_array[][2] = {{0,2}, {1,7.02}, {2,136.2}, {3,5.97}, {5,mFe*0.93+mSi*0.07}};

EOS *Fe_7Si = new EOS("Fe-7Si (Wicks)", Fe_7Si_array, sizeof(Fe_7Si_array)/2/sizeof(Fe_7Si_array[0][0]));

// -----------------------------------
// Fe-15Si shock, Wicks et al. 2018, Science Advances.

double Fe_15Si_array[][2] = {{0,2}, {1,6.784}, {2,227.9}, {3,4.74}, {5,mFe*0.85+mSi*0.15}};

EOS *Fe_15Si = new EOS("Fe-15Si (Wicks)", Fe_15Si_array, sizeof(Fe_15Si_array)/2/sizeof(Fe_15Si_array[0][0]));

// -----------------------------------
// Iron Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Fe_Seager = new EOS("Fe (Seager)", "./tabulated/iron.txt");

// ==========  Silicate ================

// ---------------------------------
// MgSiO3, perovskite, Shim & Duffy 2000, American Mineralogist

double Si_Pv_Shim_array[][2] = {{0,0}, {1,24.43}, {2,261}, {3,4}, {5,mMg+mSi+3*mO}, {8,1000}, {9,1.42}, {10,2}, {15,5}};

EOS *Si_Pv_Shim = new EOS("Si Pv (Shim)", Si_Pv_Shim_array, sizeof(Si_Pv_Shim_array)/2/sizeof(Si_Pv_Shim_array[0][0]));

// ---------------------------------
// MgSiO3, perovskite, Oganov & Ono 2004, Nature, GGA
// Ono & Oganov 2005, Earth Planet. Sci. Lett. 236, 914

double Si_Pv_array[][2] = {{0,2}, {1,25.206}, {2,230.05}, {3,4.142}, {5,mMg+mSi+3*mO}, {6, -11.2}, {8, 1000}, {9,1.506}, {10,7.02469}, {11,1.14821}, {15,5}};

EOS *Si_Pv = new EOS("Si Pv (Oganov)", Si_Pv_array, sizeof(Si_Pv_array)/2/sizeof(Si_Pv_array[0][0]));

// ---------------------------------
// MgSiO3, post-perovskite, Oganov & Ono 2004, Nature, GGA
// Ono & Oganov 2005, Earth Planet. Sci. Lett. 236, 914

double Si_PPv_array[][2] = {{0,2}, {1,25.239}, {2,199.96}, {3,4.541}, {5,mMg+mSi+3*mO}, {8, 1500}, {9,1.553}, {10,4.731}, {11,1.114}, {15,5}};

EOS *Si_PPv = new EOS("Si PPv (Oganov)", Si_PPv_array, sizeof(Si_PPv_array)/2/sizeof(Si_PPv_array[0][0]));

// ---------------------------------
// MgSiO3, post-perovskite, Sakai, Dekura, & Hirao, 2016, Scientific Reports

double Si_PPv_Sakai_array[][2] = {{0,4}, {1,24.73}, {2,203}, {3,5.35}, {5,mMg+mSi+3*mO}, {8,848}, {9,1.47}, {10,2.7}, {11,0.93}, {15,5}};

EOS *Si_PPv_Sakai = new EOS("Si PPv (Sakai)", Si_PPv_Sakai_array, sizeof(Si_PPv_Sakai_array)/2/sizeof(Si_PPv_Sakai_array[0][0]));

// ---------------------------------
// MgSiO3, liquid, Mosenfelder et al. (2009), Table 5, BM4LC, Journal of Geophysical Research: Solid Earth 38.202E-5 m^3/kg = 64.1 AA^3

double Si_liquid_array[][2] = {{0,1}, {1,64.1}, {2,24.7}, {3,9.2}, {4,-1.87}, {5,mMg+mSi+3*mO}, {9,0.37}, {10,-1.71}, {11, 0}, {15,5}, {17, 1673}, {20, 4}};

EOS *Si_liquid = new EOS("Si liquid (Mosenfelder)", Si_liquid_array, sizeof(Si_liquid_array)/2/sizeof(Si_liquid_array[0][0]));

// -----------------------------------
// MgSiO3, liquid, Wolf & Bower 2018, Table 1 S11, Using RTpress structure

double Si_Liquid_Wolf_array[][2] = {{0, 10}, {1, 38.99}, {2, 13.2}, {3, 8.238}, {5, mMg+mSi+3*mO}, {9, 0.1899}, {10, 0.6}, {12, -1.94}, {15, 5},  {17, 3000}};
double Si_Liquid_Wolf_b[] = {4.738E12, 2.97E12, 6.32E12, -1.4E13, -2.0E13};

EOS *Si_Liquid_Wolf = new EOS("Si liquid (Wolf)", Si_Liquid_Wolf_array, Si_Liquid_Wolf_b, sizeof(Si_Liquid_Wolf_array)/2/sizeof(Si_Liquid_Wolf_array[0][0]), sizeof(Si_Liquid_Wolf_b)/sizeof(Si_Liquid_Wolf_b[0]));

// -----------------------------------
// Silicate PREM mantle EOS in Appendix F.1 of Stacey & Davis 2008, used in Zeng 2016
EOS *Si_PREM = new EOS("Si (PREM)", "./tabulated/SiPREM.txt");


// -----------------------------------
// Silicate PREM BM2 extrapolation used in Zeng 2016
double Si_BM2fit_array[][2] = {{0,0}, {1,25.223}, {2,206},{3,4},{5,mMg+mSi+3*mO}};
EOS *Si_BM2fit = new EOS("Si (PREM, Zeng)", Si_BM2fit_array, sizeof(Si_BM2fit_array)/2/sizeof(Si_BM2fit_array[0][0]));

// -----------------------------------
// Silicate Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Si_Seager = new EOS("Si (Seager)", "./tabulated/silicate.txt");

// ---------------------------------
// Si Dummy. Used to fill in phase space that no EOS provided.

EOS *Si_Dummy = new EOS("Si Dummy", Si_BM2fit_array, sizeof(Si_BM2fit_array)/2/sizeof(Si_BM2fit_array[0][0]));

double dTdP_Si_Dummy (double P, double T)
// A temperature gradient that equals to the melting curve. Guarantee the temperature won't drop below the melting curve. 
{
  P /= 1E10;
  if (T > 1830*pow(1+P/4.6, 0.33))
    return 131.2826087*pow(1+P/4.6, -0.67)/1E10;
  else
  {
    cout<<"Error: The pressure "<<P<<" GPa and temperature "<<T<<" K are inconsistent with liquid silicate."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
}

// ==========  Water ================

// -----------------------------------
// Ice Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Ice_Seager = new EOS("Ice (Seager)", "./tabulated/water.txt");

// -----------------------------------
// Liquid water, ExoPlex, unkown source
double Water_ExoPlex_array[][2] = {{0,1}, {1,18.797}, {2,2.06}, {3,6.29}, {4,-0.9175}, {5,18.01528}, {7,4.184}};

EOS *Water_ExoPlex = new EOS("Water (ExoPlex)", Water_ExoPlex_array, sizeof(Water_ExoPlex_array)/2/sizeof(Water_ExoPlex_array[0][0]));

// -----------------------------------
// Liquid water, Valencia et al. 2007, ApJ, 656:545
double Water_array[][2] = {{0,0}, {1,18.047}, {2,2.18}, {5,18.01528}};

EOS *Water = new EOS("Water (Valencia)", Water_array, sizeof(Water_array)/2/sizeof(Water_array[0][0]));

// -----------------------------------
// Ice Ih, Feistel & Wagner 2006, Acuna et al. 2021

double IceIh_array[][2] = {{0,0}, {1,19.56}, {2,9.5}, {3,5.3}, {5,18.01528}, {7,1.913}};

EOS *IceIh = new EOS("Ice Ih", IceIh_array, sizeof(IceIh_array)/2/sizeof(IceIh_array[0][0]));
// -----------------------------------
// Ice Ih, ExoPlex, unkown source

double IceIh_ExoPlex_array[][2] = {{0,0}, {1,19.65}, {2,9.2}, {3,5.5}, {5,18.01528}, {7,4.184}};

EOS *IceIh_ExoPlex = new EOS("Ice Ih (ExoPlex)", IceIh_ExoPlex_array, sizeof(IceIh_ExoPlex_array)/2/sizeof(IceIh_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VI, ExoPlex, Bezacier et al. 2014

double IceVI_ExoPlex_array[][2] = {{0,0}, {1,14.17}, {2,14.01}, {3,4}, {5,18.01528}};

EOS *IceVI_ExoPlex = new EOS("Ice VI (ExoPlex)", IceVI_ExoPlex_array, sizeof(IceVI_ExoPlex_array)/2/sizeof(IceVI_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VI, Bezacier et al. 2014 & Tchijov et al. 2004

double IceVI_Bezacier_array[][2] = {{0,0}, {1,14.17}, {2,14.05}, {3,4}, {5,18.01528}, {7, 2.6}, {17, 300}, {18, 146.}};

EOS *IceVI_Bezacier = new EOS("Ice VI (Bezacier)", IceVI_Bezacier_array, sizeof(IceVI_Bezacier_array)/2/sizeof(IceVI_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII, ExoPlex, Bezacier et al. 2014

double IceVII_ExoPlex_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}};

EOS *IceVII_ExoPlex = new EOS("Ice VII (ExoPlex)", IceVII_ExoPlex_array, sizeof(IceVII_ExoPlex_array)/2/sizeof(IceVII_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VII, Bezacier et al. 2014 & Tchijov et al. 2004

double IceVII_Bezacier_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}, {7, 2.3}, {17, 300}, {18, 115.8}};

EOS *IceVII_Bezacier = new EOS("Ice VII (Bezacier)", IceVII_Bezacier_array, sizeof(IceVII_Bezacier_array)/2/sizeof(IceVII_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII, Zachary Grande

double IceVII_array[][2] = {{0,2}, {1,12.80}, {2,18.47}, {3,2.51}, {5,18.01528}};

EOS *IceVII = new EOS("Ice VII (Grande)", IceVII_array, sizeof(IceVII_array)/2/sizeof(IceVII_array[0][0]));

// -----------------------------------
// Ice VII', Zachary Grande

double IceVIIp_array[][2] = {{0,2}, {1,12.38}, {2,20.76}, {3,4.49}, {5,18.01528}};

EOS *IceVIIp = new EOS("Ice VII' (Grande)", IceVIIp_array, sizeof(IceVIIp_array)/2/sizeof(IceVIIp_array[0][0]));

// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 300K, apply 3rd Birch Murnaghan fit result in the paper to Vinet EOS.

double IceVII_FFH2004_array[][2] = {{0,2}, {1,12.4}, {2,21.1}, {3,4.4}, {5,18.01528}};

EOS *IceVII_FFH2004 = new EOS("Ice VII (FFH2004, Vinet)", IceVII_FFH2004_array, sizeof(IceVII_FFH2004_array)/2/sizeof(IceVII_FFH2004_array[0][0]));


// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 300K, fitting with 3rd order Vinet

double IceVII_FFH2004fit_array[][2] = {{0,2}, {1,12.42}, {2,19.84}, {3,4.99}, {5,18.01528}};

EOS *IceVII_FFH2004fit = new EOS("Ice VII (FFH2004fit, Vinet fit)", IceVII_FFH2004fit_array, sizeof(IceVII_FFH2004fit_array)/2/sizeof(IceVII_FFH2004fit_array[0][0]));


// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 300K, fitting with 3rd order Birch Murnaghan.

double IceVII_FFH2004BM_array[][2] = {{0,0}, {1,12.4}, {2,21.1}, {3,4.4}, {5,18.01528}};

EOS *IceVII_FFH2004BM = new EOS("Ice VII (FFH2004, BM)", IceVII_FFH2004BM_array, sizeof(IceVII_FFH2004BM_array)/2/sizeof(IceVII_FFH2004BM_array[0][0]));

// -----------------------------------
// Ice X, Hermann & Schwerdtfeger 2011, PRL 106, 187403

double IceX_HS_array[][2] = {{0,0}, {1,7.588}, {2,194.73}, {3,4}, {5,18.01528}};

EOS *IceX_HS = new EOS("Ice X (Hermann)", IceX_HS_array, sizeof(IceX_HS_array)/2/sizeof(IceX_HS_array[0][0]));

// -----------------------------------
// Ice X, Zachary Grande

double IceX_array[][2] = {{0,2}, {1,10.18}, {2,50.52}, {3,4.5}, {5,18.01528}};

EOS *IceX = new EOS("Ice X (Grande)", IceX_array, sizeof(IceX_array)/2/sizeof(IceX_array[0][0]));

// -----------------------------------
// Ice Dummy.  Used to fill in phase space that no EOS provided.

EOS *Ice_Dummy = new EOS("Ice Dummy", IceVI_ExoPlex_array, sizeof(IceVI_ExoPlex_array)/2/sizeof(IceVI_ExoPlex_array[0][0]));

double FMNR[2][13] = {{2.45, 2.5, 3.25, 3.5, 3.75, 4, 5, 6, 7, 9, 11, 13, 15}, {37.3, 43.6, 140, 183, 234, 290, 587, 966, 1440, 2703, 4405, 6416, 8893}};

double Zeng2013FFH(double P, double T)
{
  P/=1E10;			// convert to GPa

  return 18.01528 * (0.0805 + 0.0229*(1-exp(-0.0743*P)) + 0.1573*(1-exp(-0.0061*P)));

}

  
EOS *IceZeng2013FFH = new EOS("Ice (FFH 2004)", Zeng2013FFH);

EOS *IceZeng2013FMNR = new EOS("Ice (FMNR 2009)", FMNR[1], FMNR[0], 13);

double Gas_array[3][2] = {{0,6}, {5,3}, {15,2}};
EOS *Gas = new EOS("Ideal Gas", Gas_array, 3);

double Gas_iso_array[3][2] = {{0,6}, {5,3}, {15,0}};
EOS *Gas_iso = new EOS("Isothermal Ideal Gas", Gas_iso_array, 3);

double watervapor_array[3][2] = {{0,6}, {5,18}, {15,3}};
EOS *watervapor = new EOS("Water vapor", watervapor_array, 3);

// ==========  OTHER  ================
// -----------------------------------
// Gold, Heinz & Jeanloz 1983, J. Appl. Phys. (included for a Hitchiker's-related joke)
double Gold_array[][2] = {{0,0}, {1,10.215}, {2,166.65}, {3,5.4823}, {5,196.96657}, {8,170}, {9,2.95}, {12,1.7}, {15,1}};
EOS *Gold = new EOS("Gold", Gold_array, sizeof(Gold_array)/2/sizeof(Gold_array[0][0]));
// -----------------------------------
// Platinum, Matsui et al. 2009, J. Appl. Phys. (included for a Hitchiker's-related joke)
double Plat_array[][2] = {{0,2}, {1,10.03}, {2,273}, {3,5.20}, {5,195.084}, {8,230}, {9,2.7}, {10,1.1}, {15,1}};
EOS *Plat = new EOS("Plat", Plat_array, sizeof(Plat_array)/2/sizeof(Plat_array[0][0]));


// ============== An example on the format of dTdP function ==============
double dTdP_gas(double P, double T)
// return the adiabatic temperature gradient at given pressure and temperature point for ideal gas.
{
  if (P != 0)
    return 2.*T / (7.*P);
  else
  {
    cout<<"Error: Can't get adiabatic temperature gradient for diatomic gas at P=0."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
}
