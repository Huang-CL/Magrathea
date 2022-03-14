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
7.	Theta0 (K), a fitting parameter to Einstein temperature or Debye temperature
8.	gamma0, a fitting parameter of Grueneisen parameter
9.	beta, a fitting parameter of Grueneisen parameter.  In RTpress style, it represents the "m" which stands for the power-law exponent in the thermal deviation term.  Theoretically expected value: 0.6.
10.	gammainf, a fitting parameter of Grueneisen parameter
11.	gamma0p, volume derivative of the Grueneisen parameter
12.	e0 (10^-6 K^-1), electronic contribution to Helmholtz free energy
13.	g, is an electronic analogue of the Grueneisen parameter
14.	n is the number of atoms in the chemical formula of the compound.  Should have n*NA atoms within V.  The n of ideal gas is the number of atoms per molecule for the purpose of adiabatic index.  NOTE: n=2 for collinear molecules e.g. carbon dioxide!  Isothermal atmosphere can be achieved by setting n=0.
15.     Z is the atomic number (number of electron)
16.	T0, the reference temperature for the thermal pressure
17.     alpha0, the zeroth order coefficient of thermal expansion at a reference pressure P0 in 10^-6 K^-1
18.     alpha1, the first order coefficient of thermal expansion at a reference pressure P0 in 10^-6 K^-2
19.	xi, a power law index to describe the pressure effect of the coefficient of thermal expansion 
20.	cp_a in 10^7 erg/g/K Heat capacity per mass at constant pressure
21.	cp_b, fitting coefficient for specific heat capacity, in 10^7 erg/g/K^2
22.	cp_c, cp = cp_a + cp_b*T - cp_c/T^2. cp in 10^7 erg/g/K, cp_c in 10^7 erg*K/g
23.	Debye_approx, whether use Debye approximation or Einstein approximation. Debye approximation is slower but more accurate at temperature lower than Debye/Einstein temperature.  Positive number for Debye, otherwise Einstein.
24.     thermal_type, indicates the thermal type of the phase.  0 indicates no temperature profile available, 1 indicates entropy method, 2 indicates the temperature gradient method.  The only method to set the gradient is using the modify_extern_dTdP function, 3 indicates ideal gas, 4 indicates the EOS is fitted along the isentrope, type 8 indicates RTpress style.
25-32.  at1-at4 & ap1 - ap4


For RTpress style of EOS, also need a _b array. They are fitted polynomial parameters of the thermal coefficients b(V) in erg/mol.  Convert eV/atom to erg/mol need to multiply eV_erg*n*NA. For example, for MgSiO3, 0.9821 eV/atom = 4.824E12 *0.9821 erg/mol = 4.738E12 erg/mol.
*/



// ==========  Iron ================

// ---------------------------------
// Liquid Iron,  Dorogokupets et al. 2017, Scientific Reports.

double Fe_liquid_array[][2] = {{0,2}, {1,7.957}, {2,83.7}, {3,5.97}, {5,mFe}, {6, 1E-4}, {7, 263}, {8, 2.033}, {9, 1.168}, {10,0}, {12,198}, {13, 0.884}, {14, 1}, {15,26}, {16, 1811}};

EOS *Fe_liquid = new EOS("Fe liquid (Dorogokupets)", Fe_liquid_array, sizeof(Fe_liquid_array)/2/sizeof(Fe_liquid_array[0][0]));

// -----------------------------------
// Liquid Iron 2, Anderson & Ahrens, 1994 JGR

double Fe_liquid2_array[][2] = {{0,1}, {1,7.95626}, {2,109.7}, {3,4.66}, {4,-0.043}, {5,mFe}, {14,1}, {15,26}};

EOS *Fe_liquid2 = new EOS("Fe liquid (Anderson)", Fe_liquid2_array, sizeof(Fe_liquid2_array)/2/sizeof(Fe_liquid2_array[0][0]));


// -----------------------------------
// Iron hcp, Bouchet et al. 2013, PRB 87, 094102

double Fe_hcp_array[][2] = {{0,3}, {1,6.29}, {2,253.844}, {3,4.719}, {5,mFe}, {7,44.574}, {8,1.408}, {9,0.826}, {10,0.827}, {12,212.1}, {13,1.891}, {14,1}, {15,26}};

EOS *Fe_hcp = new EOS("Fe hcp (Bouchet)", Fe_hcp_array, sizeof(Fe_hcp_array)/2/sizeof(Fe_hcp_array[0][0]));

// -----------------------------------
// Iron hcp, Dorogokupets et al. 2017, Scientific Reports.

double Fe_hcp2_array[][2] = {{0,2}, {1,6.8175}, {2,148.0}, {3,5.86}, {5,mFe}, {7, 227}, {8, 2.2}, {9, 0.01}, {10,0}, {12,126}, {13,-0.83}, {14,1}, {15,26}, {16, 298.15}};

EOS *Fe_hcp2 = new EOS("Fe hcp (Dorogokupets)", Fe_hcp2_array, sizeof(Fe_hcp2_array)/2/sizeof(Fe_hcp2_array[0][0]));

// -----------------------------------
// Iron hcp, Smith et al. 2018, Nature Astronomy. <4.0 Earth-mass suggested (Gruneisen determined from Fig. 3b)

double Fe_hcp3_array[][2] = {{0,2}, {1,mFe/8.43}, {2,177.7}, {3,5.64}, {5,mFe}, {7,322}, {8,2.09}, {9,1.01}, {10,0.0500}, {14,1}, {15,26}};

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

double Si_Pv_Shim_array[][2] = {{0,0}, {1,24.43}, {2,261}, {3,4}, {5,mMg+mSi+3*mO}, {7,1000}, {8,1.42}, {9,2}, {14,5}};

EOS *Si_Pv_Shim = new EOS("Si Pv (Shim)", Si_Pv_Shim_array, sizeof(Si_Pv_Shim_array)/2/sizeof(Si_Pv_Shim_array[0][0]));

// ---------------------------------
// MgSiO3, perovskite, Oganov & Ono 2004, Nature, GGA
// Ono & Oganov 2005, Earth Planet. Sci. Lett. 236, 914

double Si_Pv_array[][2] = {{0,2}, {1,25.206}, {2,230.05}, {3,4.142}, {5,mMg+mSi+3*mO}, {6, -11.2}, {7, 1000}, {8,1.506}, {9,7.02469}, {10,1.14821}, {14,5}};

EOS *Si_Pv = new EOS("Si Pv (Oganov)", Si_Pv_array, sizeof(Si_Pv_array)/2/sizeof(Si_Pv_array[0][0]));

// ---------------------------------
// MgSiO3, post-perovskite, Oganov & Ono 2004, Nature, GGA
// Ono & Oganov 2005, Earth Planet. Sci. Lett. 236, 914

double Si_PPv_array[][2] = {{0,2}, {1,25.239}, {2,199.96}, {3,4.541}, {5,mMg+mSi+3*mO}, {7, 1500}, {8,1.553}, {9,4.731}, {10,1.114}, {14,5}};

EOS *Si_PPv = new EOS("Si PPv (Oganov)", Si_PPv_array, sizeof(Si_PPv_array)/2/sizeof(Si_PPv_array[0][0]));

// ---------------------------------
// MgSiO3, post-perovskite, Sakai, Dekura, & Hirao, 2016, Scientific Reports

double Si_PPv_Sakai_array[][2] = {{0,4}, {1,24.73}, {2,203}, {3,5.35}, {5,mMg+mSi+3*mO}, {7,848}, {8,1.47}, {9,2.7}, {10,0.93}, {14,5}};

EOS *Si_PPv_Sakai = new EOS("Si PPv (Sakai)", Si_PPv_Sakai_array, sizeof(Si_PPv_Sakai_array)/2/sizeof(Si_PPv_Sakai_array[0][0]));

// ---------------------------------
// MgSiO3, liquid, Mosenfelder et al. (2009), Table 5, BM4LC, Journal of Geophysical Research: Solid Earth 38.202E-5 m^3/kg = 64.1 AA^3

double Si_liquid_array[][2] = {{0,1}, {1,64.1}, {2,24.7}, {3,9.2}, {4,-1.87}, {5,mMg+mSi+3*mO}, {8,0.37}, {9,-1.71}, {10, 0}, {14,5}, {16, 1673}, {24, 4}};

EOS *Si_liquid = new EOS("Si liquid (Mosenfelder)", Si_liquid_array, sizeof(Si_liquid_array)/2/sizeof(Si_liquid_array[0][0]));

// -----------------------------------
// MgSiO3, liquid, Wolf & Bower 2018, Table 1 S11, Using RTpress structure

double Si_Liquid_Wolf_array[][2] = {{0, 10}, {1, 38.99}, {2, 13.2}, {3, 8.238}, {5, mMg+mSi+3*mO}, {8, 0.1899}, {9, 0.6}, {11, -1.94}, {14, 5},  {16, 3000}};
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
double Water_ExoPlex_array[][2] = {{0,1}, {1,18.797}, {2,2.06}, {3,6.29}, {4,-0.9175}, {5,18.01528}, {20,4.184}};

EOS *Water_ExoPlex = new EOS("Water (ExoPlex)", Water_ExoPlex_array, sizeof(Water_ExoPlex_array)/2/sizeof(Water_ExoPlex_array[0][0]));

// -----------------------------------
// Liquid water, Valencia et al. 2007, ApJ, 656:545
double Water_array[][2] = {{0,0}, {1,18.047}, {2,2.18}, {5,18.01528}};

EOS *Water = new EOS("Water (Valencia)", Water_array, sizeof(Water_array)/2/sizeof(Water_array[0][0]));

// -----------------------------------
// Ice Ih, Feistel & Wagner 2006, Acuna et al. 2021

double IceIh_array[][2] = {{0,0}, {1,19.56}, {2,9.5}, {3,5.3}, {5,18.01528}, {20,1.913}};

EOS *IceIh = new EOS("Ice Ih", IceIh_array, sizeof(IceIh_array)/2/sizeof(IceIh_array[0][0]));
// -----------------------------------
// Ice Ih, ExoPlex, unkown source

double IceIh_ExoPlex_array[][2] = {{0,0}, {1,19.65}, {2,9.2}, {3,5.5}, {5,18.01528}, {20,4.184}};

EOS *IceIh_ExoPlex = new EOS("Ice Ih (ExoPlex)", IceIh_ExoPlex_array, sizeof(IceIh_ExoPlex_array)/2/sizeof(IceIh_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VI, ExoPlex, Bezacier et al. 2014

double IceVI_ExoPlex_array[][2] = {{0,0}, {1,14.17}, {2,14.01}, {3,4}, {5,18.01528}};

EOS *IceVI_ExoPlex = new EOS("Ice VI (ExoPlex)", IceVI_ExoPlex_array, sizeof(IceVI_ExoPlex_array)/2/sizeof(IceVI_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VI, Bezacier et al. 2014 & Tchijov et al. 2004

double IceVI_Bezacier_array[][2] = {{0,0}, {1,14.17}, {2,14.05}, {3,4}, {5,18.01528}, {16, 300}, {17, 146.}, {20, 2.6}};

EOS *IceVI_Bezacier = new EOS("Ice VI (Bezacier)", IceVI_Bezacier_array, sizeof(IceVI_Bezacier_array)/2/sizeof(IceVI_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII, ExoPlex, Bezacier et al. 2014

double IceVII_ExoPlex_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}};

EOS *IceVII_ExoPlex = new EOS("Ice VII (ExoPlex)", IceVII_ExoPlex_array, sizeof(IceVII_ExoPlex_array)/2/sizeof(IceVII_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VII, Bezacier et al. 2014 & Tchijov et al. 2004

double IceVII_Bezacier_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}, {16, 300}, {17, 115.8}, {20, 2.3}};

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
// Ice VII, Frank, Fei & Hu 2004 original with temperature using thermal expansion representation and the heat capacity from Myint et al. 2017.

double IceVII_FFH2004T_array[][2] = {{0,0}, {1,12.4}, {2,21.1}, {3,4.4}, {5,18.01528}, {16,300}, {17, -420}, {18, 1.56}, {19, 1.1}, {20, 0}, {21, 4E-3}, {24,9}};

EOS *IceVII_FFH2004T = new EOS("Ice VII (FFH2004, thermal)", IceVII_FFH2004T_array, sizeof(IceVII_FFH2004T_array)/2/sizeof(IceVII_FFH2004T_array[0][0]));

// -----------------------------------
// Ice VII, Fei et al. 1993 modified by Sotin et al. 2007

double IceVII_Fei_array[][2] = {{0,0}, {1,12.3}, {2,23.9}, {3,4.2}, {5,18.01528}, {7, 1470}, {8, 1.2}, {9, 1}, {10, 0}, {14,3}};

EOS *IceVII_Fei = new EOS("Ice VII (Fei)", IceVII_Fei_array, sizeof(IceVII_Fei_array)/2/sizeof(IceVII_Fei_array[0][0]));

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

double Gas_array[3][2] = {{0,6}, {5,3}, {14,2}};
EOS *Gas = new EOS("Ideal Gas", Gas_array, 3);

double Gas_iso_array[3][2] = {{0,6}, {5,3}, {14,0}};
EOS *Gas_iso = new EOS("Isothermal Ideal Gas", Gas_iso_array, 3);

double watervapor_array[3][2] = {{0,6}, {5,18}, {14,3}};
EOS *watervapor = new EOS("Water vapor", watervapor_array, 3);

// ==========  OTHER  ================
// -----------------------------------
// Gold, Heinz & Jeanloz 1983, J. Appl. Phys. (included for a Hitchiker's-related joke)
double Gold_array[][2] = {{0,0}, {1,10.215}, {2,166.65}, {3,5.4823}, {5,196.96657}, {7,170}, {8,2.95}, {11,1.7}, {14,1}};
EOS *Gold = new EOS("Gold", Gold_array, sizeof(Gold_array)/2/sizeof(Gold_array[0][0]));
// -----------------------------------
// Platinum, Matsui et al. 2009, J. Appl. Phys. (included for a Hitchiker's-related joke)
double Plat_array[][2] = {{0,2}, {1,10.03}, {2,273}, {3,5.20}, {5,195.084}, {7,230}, {8,2.7}, {9,1.1}, {14,1}};
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
