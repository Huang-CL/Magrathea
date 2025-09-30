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
24.     thermal_type, indicates the thermal type of the phase.  0 indicates no temperature profile available, 1 indicates entropy method, 2 indicates the temperature gradient method.  The only method to set the gradient is using the modify_dTdP function, 3 indicates ideal gas, 4 indicates the EOS is fitted along the isentrope, type 8 indicates RTpress style.
25-32.  at1-at4 & ap1 - ap4


For RTpress style of EOS, also need a _b array. They are fitted polynomial parameters of the thermal coefficients b(V) in erg/mol.  Convert eV/atom to erg/mol need to multiply eV_erg*n*NA. For example, for MgSiO3, 0.9821 eV/atom = 4.824E12 *0.9821 erg/mol = 4.738E12 erg/mol.
*/


// ==========  Iron  ================

// ---------------------------------
// Liquid Iron, Dorogokupets et al. 2017, Scientific Reports.
// DEFAULT

double Fe_liquid_array[][2] = {{0,2}, {1,7.957}, {2,83.7}, {3,5.97}, {5,mFe}, {6, 1E-4}, {7, 263}, {8, 2.033}, {9, 1.168}, {10,0}, {12,198}, {13, 0.884}, {14, 1}, {15,26}, {16, 1811}};

EOS *Fe_liquid = new EOS("Fe liquid (Dorogokupets)", Fe_liquid_array, sizeof(Fe_liquid_array)/2/sizeof(Fe_liquid_array[0][0]));

// -----------------------------------
// Liquid Iron, Anderson & Ahrens, 1994 JGR

double Fe_liquid2_array[][2] = {{0,1}, {1,7.95626}, {2,109.7}, {3,4.66}, {4,-0.043}, {5,mFe}, {14,1}, {15,26}};

EOS *Fe_liquid2 = new EOS("Fe liquid (Anderson)", Fe_liquid2_array, sizeof(Fe_liquid2_array)/2/sizeof(Fe_liquid2_array[0][0]));

// -----------------------------------
// Alpha Iron (bcc), Dorogokupets et al. 2017, Scientific Reports 

double Fe_bcc_array[][2] = {{0,0}, {1,7.092}, {2,164.0}, {3,5.50}, {5,mFe}, {8,303}, {9,1.736}, {10,1.125}, {11,0}, {13,198}, {14,1.0}, {17,1043}};

EOS *Fe_bcc = new EOS("Fe bcc (Dorogokupets)", Fe_bcc_array, sizeof(Fe_bcc_array)/2/sizeof(Fe_bcc_array[0][0]));

// -----------------------------------
// Gamma Iron (fcc), Dorogokupets et al. 2017, Scientific Reports

double Fe_fcc_array[][2] = {{0,0}, {1,6.9285}, {2,146.2}, {3,4.67}, {5,mFe}, {8,222.5}, {9,2.203}, {10,0.01}, {11,0}, {13,198}, {14,0.5}};

EOS *Fe_fcc = new EOS("Fe fcc (Dorogokupets)", Fe_fcc_array, sizeof(Fe_fcc_array)/2/sizeof(Fe_fcc_array[0][0]));

// -----------------------------------
// Epsilon Iron (hcp), Smith et al. 2018, Nature Astronomy. (Gruneisen determined from fitting Fig. 3b)
// DEFAULT

double Fe_hcp_array[][2] = {{0,2}, {1,mFe/8.43}, {2,177.7}, {3,5.64}, {5,mFe}, {7,322}, {8,2.09}, {9,1.01}, {10,0.0500}, {14,1}, {15,26}};

EOS *Fe_hcp = new EOS("Fe hcp (Smith)", Fe_hcp_array, sizeof(Fe_hcp_array)/2/sizeof(Fe_hcp_array[0][0]));

// -----------------------------------
// Epsilon Iron (hcp), Bouchet et al. 2013, PRB 87, 094102

double Fe_hcp2_array[][2] = {{0,3}, {1,6.29}, {2,253.844}, {3,4.719}, {5,mFe}, {7,44.574}, {8,1.408}, {9,0.826}, {10,0.827}, {12,212.1}, {13,1.891}, {14,1}, {15,26}};

EOS *Fe_hcp2 = new EOS("Fe hcp (Bouchet)", Fe_hcp2_array, sizeof(Fe_hcp2_array)/2/sizeof(Fe_hcp2_array[0][0]));

// -----------------------------------
// Epsilon Iron (hcp), Dorogokupets et al. 2017, Scientific Reports.

double Fe_hcp3_array[][2] = {{0,2}, {1,6.8175}, {2,148.0}, {3,5.86}, {5,mFe}, {7, 227}, {8, 2.2}, {9, 0.01}, {10,0}, {12,126}, {13,-0.83}, {14,1}, {15,26}, {16, 298.15}};

EOS *Fe_hcp3 = new EOS("Fe hcp (Dorogokupets)", Fe_hcp3_array, sizeof(Fe_hcp3_array)/2/sizeof(Fe_hcp3_array[0][0]));

// -----------------------------------
// Iron-Silicate Alloy, 7 wt% silicate, Wicks et al. 2018, Science Advances.

double Fe_7Si_array[][2] = {{0,2}, {1,7.02}, {2,136.2}, {3,5.97}, {5,mFe*0.93+mSi*0.07}};

EOS *Fe_7Si = new EOS("Fe-7Si (Wicks)", Fe_7Si_array, sizeof(Fe_7Si_array)/2/sizeof(Fe_7Si_array[0][0]));

// -----------------------------------
// Iron-Silicate Alloy, 15 wt% silicate, Wicks et al. 2018, Science Advances.

double Fe_15Si_array[][2] = {{0,2}, {1,6.784}, {2,227.9}, {3,4.74}, {5,mFe*0.85+mSi*0.15}};

EOS *Fe_15Si = new EOS("Fe-15Si (Wicks)", Fe_15Si_array, sizeof(Fe_15Si_array)/2/sizeof(Fe_15Si_array[0][0]));

// -----------------------------------
// Iron, Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Fe_Seager = new EOS("Fe (Seager)", "./tabulated/iron.txt");

// -----------------------------------
// Iron Dummy, Used to fill in phase space that no EOS provided.

EOS *Fe_Dummy = new EOS("Fe Dummy", Fe_hcp3_array, sizeof(Fe_hcp3_array)/2/sizeof(Fe_hcp3_array[0][0]));

// ==========  Silicate  ================

// ---------------------------------
// Liquid Magnesium Silicate, MgSiO3, Mosenfelder et al. 2009, Journal of Geophysical Research: Solid Earth, Table 5, BM4LC38.202E-5 m^3/kg = 64.1 AA^3

double Si_liquid_array[][2] = {{0,1}, {1,64.1}, {2,24.7}, {3,9.2}, {4,-1.87}, {5,mMg+mSi+3*mO}, {8,0.37}, {9,-1.71}, {10, 0}, {14,5}, {16, 1673}, {24, 4}};

EOS *Si_liquid = new EOS("Si liquid (Mosenfelder)", Si_liquid_array, sizeof(Si_liquid_array)/2/sizeof(Si_liquid_array[0][0]));

// -----------------------------------
// Liquid Magnesium Silicate, MgSiO3, Wolf & Bower 2018, Table 1 S11, Using RTpress structure
// DEFAULT

double Si_Liquid_Wolf_array[][2] = {{0, 10}, {1, 38.99}, {2, 13.2}, {3, 8.238}, {5, mMg+mSi+3*mO}, {8, 0.1899}, {9, 0.6}, {11, -1.94}, {14, 5},  {16, 3000}};
double Si_Liquid_Wolf_b[] = {4.738E12, 2.97E12, 6.32E12, -1.4E13, -2.0E13};

EOS *Si_Liquid_Wolf = new EOS("Si liquid (Wolf)", Si_Liquid_Wolf_array, Si_Liquid_Wolf_b, sizeof(Si_Liquid_Wolf_array)/2/sizeof(Si_Liquid_Wolf_array[0][0]), sizeof(Si_Liquid_Wolf_b)/sizeof(Si_Liquid_Wolf_b[0]));

//----------------------------------------
// Forsterite/Olivine, Mg2Si04, Dorogokupets et al. 2015, Russ. Geol. Geophys.
double Fo_array[][2] = {{0,0}, {1,43.67}, {2,127.4}, {3,4.3}, {5,2*mMg+mSi+4*mO}, {7,949}, {8,1.066}, {9,2.225}, {10,0.0}, {14,7}};

EOS *Fo = new EOS("Fo/Ol (Dorogokupets)", Fo_array, sizeof(Fo_array)/2/sizeof(Fo_array[0][0]));

//----------------------------------------
// Wadsleyite, Mg2Si04, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Wds_array[][2] = {{0,0}, {1,40.54}, {2,169.0}, {3,4.14}, {5,2*mMg+mSi+4*mO}, {7,921}, {8,1.185}, {9,2.10}, {10,0.0}, {14,7}};

EOS *Wds = new EOS("Wds (Dorogokupets)", Wds_array, sizeof(Wds_array)/2/sizeof(Wds_array[0][0]));

//----------------------------------------
//Ringwoodite, Mg2Si04, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Rwd_array[][2] = {{0,0}, {1,39.5}, {2,187.4}, {3,3.98}, {5,2*mMg+mSi+4*mO}, {7,929}, {8,1.21}, {9,1.35}, {10,0.0}, {14,7}};

EOS *Rwd = new EOS("Rwd (Dorogokupets)", Rwd_array, sizeof(Rwd_array)/2/sizeof(Rwd_array[0][0]));

//----------------------------------------
// Akimotoite, MgSi03, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Akm_array[][2] = {{0,0}, {1,26.35}, {2,215.3}, {3,4.91}, {5,108.27}, {7,995}, {8,2.00}, {9,1.41}, {10,0.0}, {14,5}};

EOS *Akm = new EOS("Akm (Dorogokupets et al.)", Akm_array, sizeof(Akm_array)/2/sizeof(Akm_array[0][0]));

//--------------------------------------
// Forsterite/Olivine, Mg2SiO4, Sotin et al. 2007, Icarus, Using Duffy et al. 1995 & Bouhifd et al. 1996

double Fo_Sotin_array[][2] = {{0,0}, {1,(2*mMg+mSi+4*mO)/3.22}, {2,128}, {3,4.3}, {5,2*mMg+mSi+4*mO}, {14,7}, {16, 300}, {17,28.32}, {18, 0.00758}, {20,0.840}};

EOS *Fo_Sotin = new EOS("Fo/Ol (Sotin)", Fo_Sotin_array, sizeof(Fo_Sotin_array)/2/sizeof(Fo_Sotin_array[0][0]));

//---------------------------------------
// Enstatite/Orthopyroxene, MgSiO3, Sotin et al. 2007, Icarus, Using Vacher et al. 1998 & Anderson et al. 1991

double En_array[][2] = {{0,0}, {1,(mMg+mSi+3*mO)/3.215}, {2,111}, {3,7}, {5,mMg+mSi+3*mO}, {14,5}, {16, 300}, {17,28.6}, {18, 0.0072}, {20,0.840}};

EOS *En = new EOS("En/Opx (Sotin)", En_array, sizeof(En_array)/2/sizeof(En_array[0][0]));

//--------------------------------------
// Magnesiowustite, MgO, Sotin et al. 2007, Icarus,

double Mw_array[][2] = {{0,0}, {1,(mMg+mO)/3.584}, {2,157}, {3,4.4}, {5,mMg+mO}, {7,430}, {8, 1.45}, {9,3}, {14,2}};

EOS *Mw = new EOS("Magnesiowustite (Sotin)", Mw_array, sizeof(Mw_array)/2/sizeof(Mw_array[0][0]));

// ---------------------------------
// Bridgmanite/Perovskite, MgSiO3, Oganov & Ono 2004, Nature, GGA
// DEFAULT

double Si_Pv_array[][2] = {{0,2}, {1,25.206}, {2,230.05}, {3,4.142}, {5,mMg+mSi+3*mO}, {6, -11.2}, {7, 1000}, {8,1.506}, {9,7.02469}, {10,1.14821}, {14,5}};

EOS *Si_Pv = new EOS("Brg (Oganov)", Si_Pv_array, sizeof(Si_Pv_array)/2/sizeof(Si_Pv_array[0][0]));

// ---------------------------------
// Bridgmanite/Perovskite, MgSiO3, Shim & Duffy 2000, American Mineralogist

double Si_Pv_Shim_array[][2] = {{0,0}, {1,24.43}, {2,261}, {3,4}, {5,mMg+mSi+3*mO}, {7,1000}, {8,1.42}, {9,2}, {14,5}};

EOS *Si_Pv_Shim = new EOS("Brg (Shim)", Si_Pv_Shim_array, sizeof(Si_Pv_Shim_array)/2/sizeof(Si_Pv_Shim_array[0][0]));

//----------------------------------------
// Bridgmanite/Perovskite, MgSi03, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Pv_Doro_array[][2] = {{0,0}, {1,24.45}, {2,252.0}, {3,4.38}, {5, mMg+mSi+3*mO}, {7,943}, {8,1.70}, {9,3.00}, {10,0.0}, {14,5}};

EOS *Pv_Doro = new EOS("Pv (Dorogokupets)", Pv_Doro_array, sizeof(Pv_Doro_array)/2/sizeof(Pv_Doro_array[0][0]));

// ---------------------------------
// Post-Perovskite, MgSiO3, Sakai, Dekura, & Hirao, 2016, Scientific Reports
// DEFAULT

double Si_PPv_Sakai_array[][2] = {{0,4}, {1,24.73}, {2,203}, {3,5.35}, {5,mMg+mSi+3*mO}, {7,848}, {8,1.47}, {9,2.7}, {10,0.93}, {14,5}};

EOS *Si_PPv_Sakai = new EOS("Si PPv (Sakai)", Si_PPv_Sakai_array, sizeof(Si_PPv_Sakai_array)/2/sizeof(Si_PPv_Sakai_array[0][0]));

// ---------------------------------
// Post-Perovskite, MgSiO3, Oganov & Ono 2004, Nature, GGA

double Si_PPv_array[][2] = {{0,2}, {1,25.239}, {2,199.96}, {3,4.541}, {5,mMg+mSi+3*mO}, {7, 1500}, {8,1.553}, {9,4.731}, {10,1.114}, {14,5}};

EOS *Si_PPv = new EOS("Si PPv (Oganov)", Si_PPv_array, sizeof(Si_PPv_array)/2/sizeof(Si_PPv_array[0][0]));

//----------------------------------------
// Postperovskite, MgSi03, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double PPv_Doro_array[][2] = {{0,0}, {1,24.2}, {2,253.7}, {3,4.03}, {5,mMg+mSi+3*mO}, {7,943}, {8,1.67}, {9,2.22}, {10,0.0}, {14,5}};

EOS *PPv_Doro = new EOS("PPv (Dorogokupets)", PPv_Doro_array, sizeof(PPv_Doro_array)/2/sizeof(PPv_Doro_array[0][0]));

// -----------------------------------
// Silicate PREM mantle EOS in Appendix F.1 of Stacey & Davis 2008, used in Zeng 2016
EOS *Si_PREM = new EOS("Si (PREM)", "./tabulated/SiPREM.txt");

// -----------------------------------
// Silicate PREM BM2 extrapolation used in Zeng 2016
double Si_BM2fit_array[][2] = {{0,0}, {1,25.223}, {2,206},{3,4},{5,mMg+mSi+3*mO}};
EOS *Si_BM2fit = new EOS("Si (PREM, Zeng)", Si_BM2fit_array, sizeof(Si_BM2fit_array)/2/sizeof(Si_BM2fit_array[0][0]));

// -----------------------------------
// Silicate, Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Si_Seager = new EOS("Si (Seager)", "./tabulated/silicate.txt");

// ---------------------------------
// Si Dummy, Used to fill in phase space that no EOS provided.

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

// ==========  Water  ================

// -----------------------------------
// Water, Bollengier et al. (2019), as tabulated through SeaFreeze, Journaux et al. 2020, JGR Planets, 124, 1
//DEFAULT
EOS *Water_SF = new EOS("Water (Bollengier/SF)", "./tabulated/SFwater1.txt");

// -----------------------------------
// Water Liquid, IAPWS-R6-95, Wagner & PruB 2022, JPCRD, 31, 387
//DEFAULT
EOS *Water_IAPWS = new EOS("Water (IAPWS)", "./tabulated/water_IAPWS95.txt");

// -----------------------------------
// Liquid water, ExoPlex, unkown source
double Water_ExoPlex_array[][2] = {{0,1}, {1,18.797}, {2,2.06}, {3,6.29}, {4,-0.9175}, {5,18.01528}, {20,4.184}};

EOS *Water_ExoPlex = new EOS("Water (ExoPlex)", Water_ExoPlex_array, sizeof(Water_ExoPlex_array)/2/sizeof(Water_ExoPlex_array[0][0]));

// -----------------------------------
// Liquid water, Valencia et al. 2007, ApJ, 656:545
double Water_array[][2] = {{0,0}, {1,18.047}, {2,2.18}, {5,18.01528}};

EOS *Water = new EOS("Water (Valencia)", Water_array, sizeof(Water_array)/2/sizeof(Water_array[0][0]));

// -----------------------------------
// Water & Super Critical, Brown 2018, Fluid Phase Equilib., 463, 18, tabulated by SeaFreeze
// DEFAULT
EOS *Water_Brown = new EOS("Water SC (Brown)", "./tabulated/water_Brown.txt");

// -----------------------------------
// Supercritical water. Mazevet et al. 2019, A&A 621
// https://www.ioffe.ru/astro/H2O/index.html
// DEFAULT
EOS *Water_sc_Mazevet = new EOS("Water SC (Mazevet)", H2OSC);

// -----------------------------------
// Water/Ice, SeaFreeze tabulated, Journaux et al. 2020, JGR Planets, 124, 1
// DEFAULT
EOS *IceIh_SF = new EOS("Ice Ih (SeaFreeze)", "./tabulated/SFiceih.txt");

// -----------------------------------
// Water/Ice, SeaFreeze tabulated, Journaux et al. 2020, JGR Planets, 124, 1
// DEFAULT
EOS *IceII_SF = new EOS("Ice II (SeaFreeze)", "./tabulated/SFiceii.txt");

// -----------------------------------
// Water/Ice, SeaFreeze tabulated, Journaux et al. 2020, JGR Planets, 124, 1
// DEFAULT
EOS *IceIII_SF = new EOS("Ice III (SeaFreeze)", "./tabulated/SFiceiii.txt");

// -----------------------------------
// Water/Ice, SeaFreeze tabulated, Journaux et al. 2020, JGR Planets, 124, 1
// DEFAULT
EOS *IceV_SF = new EOS("Ice V (SeaFreeze)", "./tabulated/SFiceV.txt");

// -----------------------------------
// Water/Ice, SeaFreeze tabulated, Journaux et al. 2020, JGR Planets, 124, 1
// DEFAULT
EOS *IceVI_SF= new EOS("Ice VI (SeaFreeze)", "./tabulated/SFiceVI.txt");

// -----------------------------------
// Ice Ih, Feistel & Wagner 2006, Acuna et al. 2021
// DEFAULT

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
// DEFAULT

double IceVI_Bezacier_array[][2] = {{0,0}, {1,14.17}, {2,14.05}, {3,4}, {5,18.01528}, {16, 300}, {17, 146.}, {20, 2.6}};

EOS *IceVI_Bezacier = new EOS("Ice VI (Bezacier)", IceVI_Bezacier_array, sizeof(IceVI_Bezacier_array)/2/sizeof(IceVI_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII, Bezacier et al. 2014 & Tchijov et al. 2004
// DEFAULT

double IceVII_Bezacier_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}, {16, 300}, {17, 115.8}, {20, 2.3}};

EOS *IceVII_Bezacier = new EOS("Ice VII (Bezacier)", IceVII_Bezacier_array, sizeof(IceVII_Bezacier_array)/2/sizeof(IceVII_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII Isothermal, Bezacier et al. 2014

double IceVII_ExoPlex_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}};

EOS *IceVII_ExoPlex = new EOS("Ice VII (ExoPlex)", IceVII_ExoPlex_array, sizeof(IceVII_ExoPlex_array)/2/sizeof(IceVII_ExoPlex_array[0][0]));

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
// Ice X, Zachary Grande
// DEFAULT

double IceX_array[][2] = {{0,2}, {1,10.18}, {2,50.52}, {3,4.5}, {5,18.01528}};

EOS *IceX = new EOS("Ice X (Grande)", IceX_array, sizeof(IceX_array)/2/sizeof(IceX_array[0][0]));

// -----------------------------------
// Ice X, Hermann & Schwerdtfeger 2011, PRL 106, 187403

double IceX_HS_array[][2] = {{0,0}, {1,7.588}, {2,194.73}, {3,4}, {5,18.01528}};

EOS *IceX_HS = new EOS("Ice X (Hermann)", IceX_HS_array, sizeof(IceX_HS_array)/2/sizeof(IceX_HS_array[0][0]));

// -----------------------------------
// Ice, Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Ice_Seager = new EOS("Ice (Seager)", "./tabulated/water.txt");

// -----------------------------------
// Water/Ice/Supercritical/Vapor, AQUA tabulated, Haldermann et al. 2022, A&A 643, A105 
EOS *H2O_AQUA = new EOS("H2O (AQUA)", "./tabulated/AQUA.txt");

// -----------------------------------
// Water/Ice, SeaFreeze tabulated, Journaux et al. 2020, JGR Planets, 124, 1
EOS *H2O_SeaFreeze = new EOS("H2O (SeaFreeze)", "./tabulated/SeaFreeze.txt");

// -----------------------------------
// Ice Dummy,  Used to fill in phase space that no EOS provided.

EOS *Ice_Dummy = new EOS("Ice Dummy", IceVI_ExoPlex_array, sizeof(IceVI_ExoPlex_array)/2/sizeof(IceVI_ExoPlex_array[0][0]));

// -----------------------------------
// Modified EOSs to match Zeng 2013
double FMNR[2][13] = {{2.45, 2.5, 3.25, 3.5, 3.75, 4, 5, 6, 7, 9, 11, 13, 15}, {37.3, 43.6, 140, 183, 234, 290, 587, 966, 1440, 2703, 4405, 6416, 8893}};

double Zeng2013FFH(double P, double T, double rho_guess)
{
  P/=1E10;			// convert to GPa

  return 18.01528 * (0.0805 + 0.0229*(1-exp(-0.0743*P)) + 0.1573*(1-exp(-0.0061*P)));
}
  
EOS *IceZeng2013FFH = new EOS("Ice (FFH 2004)", Zeng2013FFH);

EOS *IceZeng2013FMNR = new EOS("Ice (FMNR 2009)", FMNR[1], FMNR[0], 13);

// =========  Atmosphere  ================

// -----------------------------------
// Adiabtic Ideal Gas, Parameter 5 can be changed for mean molecular weight of gas. 3 g/mol = mix of H/He
// DEFAULT
double Gas_array[3][2] = {{0,6}, {5,2}, {14,2}};

EOS *Gas = new EOS("Ideal Gas", Gas_array, 3);

// -----------------------------------
// van der Waals, H2, 
double vdW_H2_array[5][2] = {{0,6}, {5,2}, {14,2}, {33, 0.2452}, {34, 0.0265}};

EOS *vdW_H2 = new EOS("H2 vdW", vdW_H2_array, 5);

// -----------------------------------
// van der Waals, He, 
double vdW_He_array[5][2] = {{0,6}, {5,4}, {14,1}, {33, 0.0346}, {34, 0.0238}};

EOS *vdW_He = new EOS("He vdW", vdW_He_array, 5);

// -----------------------------------
// van der Waals, H2O, 
//Default Hydrosphere
double vdW_H2O_array[5][2] = {{0,6}, {5,18}, {14,3}, {33, 5.537}, {34, 0.0305}};

EOS *vdW_H2O = new EOS("H2O vdW", vdW_H2O_array, 5);

// -----------------------------------
// van der Waals, CH4, 
double vdW_CH4_array[5][2] = {{0,6}, {5,16}, {14,3}, {33, 2.303}, {34, 0.0431}};

EOS *vdW_CH4 = new EOS("CH4 vdW", vdW_CH4_array, 5);

// -----------------------------------
// van der Waals, NH3, 
double vdW_NH3_array[5][2] = {{0,6}, {5,17}, {14,3}, {33, 4.225}, {34, 0.0371}};

EOS *vdW_NH3 = new EOS("NH3 vdW", vdW_NH3_array, 5);

// -----------------------------------
// van der Waals, CO2, 
double vdW_CO2_array[5][2] = {{0,6}, {5,44}, {14,2}, {33, 3.658}, {34, 0.0429}};

EOS *vdW_CO2 = new EOS("CO2 vdW", vdW_CO2_array, 5);

// -----------------------------------
// Isothermal Ideal Gas
// DEFAULT
double Gas_iso_array[3][2] = {{0,6}, {5,3}, {14,0}};

EOS *Gas_iso = new EOS("Isothermal Ideal Gas", Gas_iso_array, 3);

// -----------------------------------
// Water Vapor Ideal Gas, same as adiabatic but 5 and 14 changed for water
double watervapor_array[3][2] = {{0,6}, {5,18}, {14,3}};

EOS *watervapor = new EOS("Water vapor", watervapor_array, 3);

// -----------------------------------
// Water Vapor, IAPWS-R6-95, Wagner & PruB 2022, JPCRD, 31, 387
// Default Hydrosphere
EOS *Water_Vap_IAPWS = new EOS("H20 Vapor (IAPWS)", "./tabulated/water_IAPWS95.txt");

// -----------------------------------
// H/He, Chabrier & Debras 2021 Apj, Y=0.275
EOS *Gas_hhe = new EOS("H/He (Chabrier)", "./tabulated/ChabrierHHe0275.txt");

// ==========  CARBON  ================

// Graphite (BME), Seager et al. 2007, ApJ
double Graph_array[][2] = {{0,0},{1,5.33822},{2,33.8},{3,8.9},{5,mC}};   

EOS *Graph = new EOS("Graph",Graph_array,sizeof(Graph_array)/2/sizeof(Graph_array[0][0]));

// -----------------------------------
// Diamond (Vinet), Benedict et al. 2018, Phys Rev B
// Uses the 1 terms from the double Debye model found in Benedict et al (2014)
double Diamond_array[][2] = {{0,2},{1,5.7304},{2,432.4},{3,3.793},{5,mC},{7,1887.8},{8,0.5836},{9,1},{10,0.499},{14,1},{15,6}};

EOS *Diam = new EOS("Diamond",Diamond_array,sizeof(Diamond_array)/2/sizeof(Diamond_array[0][0]));
// -----------------------------------
// BC8 (Vinet), Benedict et al. 2018, Phys  Rev B
// Uses the 1 terms from the double Debye model found in Benedict et al (2018)

double BC8_array[][2] = {{0,2},{1,3.75902},{2,221.2},{3,4.697},{5,mC},{7,2800.6},{8,0.561},{9,1},{10,0.449},{14,1},{15,6}};

EOS *BC8 = new EOS("BC8",BC8_array,sizeof(BC8_array)/2/sizeof(BC8_array[0][0]));


// Silicon Carbide B3 (Zinc Blende) - Vinet EOS for better high-pressure extrapolation
// Based on Miozzi et al. 2018 parameters but using Vinet equation
double SiC_B3_Vinet_array[][2] = {{0, 2}, {1, 12.47}, {2, 224}, {3, 4.08}, {5, 40.096}, {7, 1200}, {8, 1.06}, {10, 0}, {14, 2}, {17, 6.2}, {16, 300}};

EOS *SiC_B3_Vinet = new EOS("SiC B3 Vinet (Miozzi)", SiC_B3_Vinet_array, sizeof(SiC_B3_Vinet_array)/2/sizeof(SiC_B3_Vinet_array[0][0]));

// Silicon Carbide B1 (Rock Salt) - Vinet EOS for better high-pressure extrapolation  
// Based on Miozzi et al. 2018 parameters but using Vinet equation
// Added thermal parameters to match original B1 implementation
double SiC_B1_Vinet_array[][2] = {{0, 2}, {1, 9.90}, {2, 339}, {3, 3.06}, {5, 40.096}, {7, 1200}, {8, 0.50}, {9, 1.67}, {10, 0}, {14, 2}, {16, 300}};

EOS *SiC_B1_Vinet = new EOS("SiC B1 Vinet (Miozzi)", SiC_B1_Vinet_array, sizeof(SiC_B1_Vinet_array)/2/sizeof(SiC_B1_Vinet_array[0][0]));

// ==========  OTHER  ================

// -----------------------------------
// Gold, Heinz & Jeanloz 1983, J. Appl. Phys. (included for a Hitchiker's-related joke)
double Gold_BM3_array[][2] = {{0,0}, {1,10.215}, {2,166.65}, {3,5.4823}, {5,196.96657}, {7,170}, {8,2.95}, {11,1.7}, {14,1}};

EOS *Gold_BM3 = new EOS("Gold", Gold_BM3_array, sizeof(Gold_BM3_array)/2/sizeof(Gold_BM3_array[0][0]));

// -----------------------------------
// Gold, Matsui et al. 2010, J. Phys. (included for a Hitchiker's-related joke)
double Gold_array[][2] = {{0,2}, {1,10.215}, {2,167}, {3,6}, {5,196.96657}, {7,170}, {8,2.97}, {9,1.6}, {14,1}};

EOS *Gold = new EOS("Gold", Gold_array, sizeof(Gold_array)/2/sizeof(Gold_array[0][0]));

// -----------------------------------
// Platinum, Matsui et al. 2009, J. Appl. Phys. (included for a Hitchiker's-related joke)
double Plat_array[][2] = {{0,2}, {1,9.0902}, {2,273}, {3,5.20}, {5,195.084}, {7,230}, {8,2.7}, {9,1.1}, {14,1}};

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

double PH2OSC(double rho,double T)
// Mazevet et al. 2019, A&A, 621
// https://www.ioffe.ru/astro/H2O/index.html
{
  constexpr double UN_T6   = 0.3157746;             // ha/kB × 1e−6
  constexpr double C13     = 1.0 / 3.0;
  constexpr double Zmean   = 10.0 / 3.0;
  constexpr double CMImean = 18.0 / 3.0;
  constexpr double DENSCONV = 11.20587 * CMImean;   // g cm⁻³ → n_i  [au]
  constexpr double TnkCONV  = 8.31447e13 / CMImean; // n_i kT factor [1.0e-6 erg cm⁻³]
  constexpr double aW = 2.357;
  constexpr double bW = 340.8;
  // many extra params:
  constexpr double P1=2.35, P3=5.9,  P4=3.78, P5=17.0, P7=1.5, P8=0.09;
  constexpr double QW=0.00123797, PW=2.384, PQ=1.5;
  constexpr double Q4=4.0, Q1=0.4, Q2=90.0;

  const double T6   = T / 1.0e6;                   // K→MK
  const double TEMP = T6 / UN_T6;                  // to atomic units (Ha/kB)

  const double DENSI   = rho / DENSCONV;           // ion # density  [au]
  const double DENSMOL = DENSI / 3.0;              // pseudo‑molecule density

  // r_s (electron density parameter)
  const double RS  = pow(0.75 / pi / DENSI / Zmean, C13);
  const double GAME = 1.0 / (RS * TEMP);           // Γ_e
  const double GAMEsq = sqrt(GAME);
   
  //  1) super‑ionic / plasma contribution 
  const double ZNA    = 1.0 + P8 / RS / GAMEsq;
  const double ZNA1RS = -P8 / RS / GAMEsq;
  const double ZNA1G  =  0.5 * ZNA1RS;

  const double ZNB     = P1 * RS / ZNA;
  const double ZNB1RS  = ZNB * (1.0 - ZNA1RS / ZNA);
  const double ZNB1G   = -ZNB * ZNA1G / ZNA;

  const double ZNC   = 1.0 + P5 / GAME;
  const double RS4   = pow(RS, P4);
  const double ZNC4  = ZNC * sqrt(ZNC);       // ZNC^(1+0.5) = ZNC^1.5
  const double ZNE   = P3 * RS4 / ZNC4;

  const double ZNE1RS = P4 * ZNE;
  const double ZNE1G  = ZNE * P7 / ZNC * P5 / GAME;

  const double ZN     = 1.0 + ZNB + ZNE;
  const double ZN1RS  = ZNB1RS + ZNE1RS;
  const double ZN1G   = ZNB1G  + ZNE1G;

  const double ZN1R = (ZN1G - ZN1RS) / 3.0;

  const double ZEF  = Zmean / ZN;                 // effective ⟨Z⟩
  const double ZDR  = -ZN1R / ZN;

  const double DENSEF = DENSI * ZEF;

  //  electron‑gas EOS (ideal Fermi gas) 
  double CHI, FE, PE, UE, SE, CVE, CHITE, CHIRE;
  ELECNR(DENSEF, TEMP,
         CHI, FE, PE, UE, SE, CVE, CHITE, CHIRE);

  const double FNkTsi = FE * ZEF;

  const double FEDR  = PE * (1.0 + ZDR);

  const double FDR = FEDR * ZEF + FE * ZEF * ZDR;

  const double PnkTsi = FDR;

  //  2) non‑ideal molecular piece
  const double cW    = 1.0 + pow(QW / TEMP, PW);
  const double bWPQ  = bW * DENSMOL * sqrt(bW * DENSMOL); // (bW*ρ)^1.5

  const double FNkTmol = (-aW * DENSMOL / TEMP + bW * DENSMOL + bWPQ * cW / PQ) / 3.0;
  const double PnkTmol = (-aW * DENSMOL / TEMP + bW * DENSMOL + bWPQ * cW)      / 3.0;

  //  3) mix molecular + plasma via smooth switch YL/YH
  const double X   = Q4 * log(Q1 * rho + Q2 * T);
  const double X1R = Q4 * Q1 * rho   / (Q1 * rho + Q2 * T);

  const double YL  = FERMIF(X);
  const double YH  = 1.0 - YL;

  const double YH1X = YH * YL;

  const double YH1R = YH1X * X1R;

  const double PnkTni = PnkTmol * YL + PnkTsi * YH + (FNkTsi - FNkTmol) * YH1R;

  //  4) Ideal gas of pseudo‑molecules
  const double PnkTid = C13;

  //  5) total (non‑ideal + ideal)
  double PnkT,PGPa;
  PnkT = PnkTni + PnkTid;

  //  7) auxiliary conversions
  const double Tnk = TnkCONV * rho * T6;
  PGPa = PnkT * Tnk / 1.0e10;   // → GPa
  // PGPa = P/1E10 = PnkT * n k T / 1E10 = PnkT * rho * NA *k / MImean * T / 1E10
  return PGPa;
}


double dTdP_S_H2O_of_rho(double rho, double T)
// n_i number per volume
// N_i number per mass
// PnkT - pressure / n_i kT, where n_i=2*n_H+n_O=3*n_{H2O}
// FNkT - free energy / N_i kT (up to an additive constant)
// UNkT - internal energy / N_i kT
// CV - heat capacity per unit volume /kN_i
// CHIT - logarithmic derivative of pressure over temperature at constant density
// CHIR - logarithmic derivative of pressure over density at constant temperature
{
  constexpr double UN_T6   = 0.3157746;             // ha/kB × 1e−6
  constexpr double C13     = 1.0 / 3.0;
  constexpr double AUM     = 1822.88848;            // m_u / m_e
  constexpr double Zmean   = 10.0 / 3.0;
  constexpr double CMImean = 18.0 / 3.0;
  constexpr double DENSCONV = 11.20587 * CMImean;   // g cm⁻³ → n_i  [au]
  constexpr double TnkCONV  = 8.31447e13 / CMImean; // n_i k factor [1e−6 erg cm⁻³]
  constexpr double aW = 2.357;
  constexpr double bW = 340.8;
  constexpr double TCRIT = 0.00205;                 // 647.15 K in Hartree (k_B=1)
  // many extra params:
  constexpr double P1=2.35, P3=5.9,  P4=3.78, P5=17.0, P7=1.5, P8=0.09;
  constexpr double QW=0.00123797, PW=2.384, PQ=1.5;
  constexpr double Q4=4.0, Q1=0.4, Q2=90.0;
  constexpr double PC1=0.0069, PC2=0.0031, PC3=0.00558, PC4=0.019;
  constexpr double SBASE = 4.9;                     // entropy offset

  const double T6   = T / 1.0e6;                   // K→MK
  const double TEMP = T6 / UN_T6;                  // to atomic units (Ha/kB)

  const double DENSI   = rho / DENSCONV;           // ion # density  [au]
  const double DENSMOL = DENSI / 3.0;              // pseudo‑molecule density

  // r_s (electron density parameter)
  const double RS  = pow(0.75 / pi / DENSI / Zmean, C13);
  const double GAME = 1.0 / (RS * TEMP);           // Γ_e
  const double GAMEsq = sqrt(GAME);
   
  //  1) super‑ionic / plasma contribution 
  const double ZNA    = 1.0 + P8 / RS / GAMEsq;
  const double ZNA1RS = -P8 / RS / GAMEsq;
  const double ZNA1G  =  0.5 * ZNA1RS;
  const double ZNA2RS = -ZNA1RS;
  const double ZNA2RSG=  0.5 * ZNA2RS;
  const double ZNA2G  =  0.5 * ZNA2RSG;

  const double ZNB     = P1 * RS / ZNA;
  const double ZNB1RS  = ZNB * (1.0 - ZNA1RS / ZNA);
  const double ZNB1G   = -ZNB * ZNA1G / ZNA;

  const double ZNB2RS  = ZNB1RS * (1.0 - ZNA1RS / ZNA) - ZNB * ZNA2RS / ZNA
                          + ZNB * pow(ZNA1RS / ZNA, 2);
  const double ZNB2G   = -ZNB1G * ZNA1G / ZNA - ZNB * ZNA2G / ZNA
                          + ZNB * pow(ZNA1G / ZNA, 2);
  const double ZNB2RSG = -ZNB1RS * ZNA1G / ZNA - ZNB * ZNA2RSG / ZNA
                          + ZNB * ZNA1G * ZNA1RS / (ZNA*ZNA);

  const double ZNC   = 1.0 + P5 / GAME;
  const double RS4   = pow(RS, P4);
  const double ZNC4  = ZNC * sqrt(ZNC);       // ZNC^(1+0.5) = ZNC^1.5
  const double ZNE   = P3 * RS4 / ZNC4;

  const double ZNE1RS = P4 * ZNE;
  const double ZNE1G  = ZNE * P7 / ZNC * P5 / GAME;
  const double ZNE2RS = P4 * ZNE1RS;
  const double ZNE2G  = ZNE1G * (ZNE1G / ZNE + P5 / GAME / ZNC - 1.0);
  const double ZNE2RSG= P4 * ZNE1G;

  const double ZN     = 1.0 + ZNB + ZNE;
  const double ZN1RS  = ZNB1RS + ZNE1RS;
  const double ZN1G   = ZNB1G  + ZNE1G;
  const double ZN2RS  = ZNB2RS + ZNE2RS;
  const double ZN2G   = ZNB2G  + ZNE2G;
  const double ZN2RSG = ZNB2RSG+ ZNE2RSG;

  const double ZN1R = (ZN1G - ZN1RS) / 3.0;
  const double ZN1T = -ZN1G;
  const double ZN2R = (ZN2RS - 2.0 * ZN2RSG + ZN2G) / 9.0;
  const double ZN2T = ZN2G;
  const double ZN2RT= -(ZN2G - ZN2RSG) / 3.0;

  const double ZEF  = Zmean / ZN;                 // effective ⟨Z⟩
  const double ZDR  = -ZN1R / ZN;
  const double ZDT  = -ZN1T / ZN;
  const double ZDRR = -ZN2R / ZN + pow(ZN1R / ZN, 2);
  const double ZDTT = -ZN2T / ZN + pow(ZN1T / ZN, 2);
  const double ZDRT = ZN1R * ZN1T / (ZN*ZN) - ZN2RT / ZN;

  const double DENSEF = DENSI * ZEF;

  //  electron‑gas EOS (ideal Fermi gas) 
  double CHI, FE, PE, UE, SE, CVE, CHITE, CHIRE;
  ELECNR(DENSEF, TEMP,
         CHI, FE, PE, UE, SE, CVE, CHITE, CHIRE);

  const double FNkTsi = FE * ZEF;

  const double FEDR  = PE * (1.0 + ZDR);
  const double FEDT  = -UE + PE * ZDT;
  const double FEDRR = PE * (CHIRE - 1.0) * pow(1.0 + ZDR, 2) + PE * ZDRR;
  const double FEDRT = (PE * (CHITE - 1.0) + PE * (CHIRE - 1.0) * ZDT) *
                       (1.0 + ZDR) + PE * ZDRT;
  const double FEDTT = UE - CVE + PE * (CHITE - 1.0) * ZDT * 2.0
                       + PE * (CHIRE - 1.0) * pow(ZDT, 2)
                       + PE * ZDTT;

  const double FDR = FEDR * ZEF + FE * ZEF * ZDR;
  const double FDT = FEDT * ZEF + FE * ZEF * ZDT;

  const double FsiDRR = FEDRR * ZEF + 2.0 * FEDR * ZEF * ZDR
                         + FE * ZEF * (pow(ZDR, 2) + ZDRR);
  const double FsiDTT = FEDTT * ZEF + 2.0 * FEDT * ZEF * ZDT
                         + FE * ZEF * (pow(ZDT, 2) + ZDTT);
  const double FsiDRT = FEDRT * ZEF + FEDR * ZEF * ZDT + FEDT * ZEF * ZDR
                         + FE * ZEF * (ZDR * ZDT + ZDRT);

  const double PnkTsi = FDR;
  const double UNkTsi = -FDT;

  //  2) non‑ideal molecular piece
  const double cW    = 1.0 + pow(QW / TEMP, PW);
  const double cW1T  = -PW * pow(QW / TEMP, PW);
  const double cW2T  = -PW * cW1T;
  const double bWPQ  = bW * DENSMOL * sqrt(bW * DENSMOL); // (bW*ρ)^1.5

  const double FNkTmol = (-aW * DENSMOL / TEMP + bW * DENSMOL + bWPQ * cW / PQ) / 3.0;
  const double PnkTmol = (-aW * DENSMOL / TEMP + bW * DENSMOL + bWPQ * cW)      / 3.0;
  const double UNkTmol = -( aW * DENSMOL / TEMP + bWPQ * cW1T / PQ)              / 3.0;

  const double FmDRR = (-aW * DENSMOL / TEMP + bW * DENSMOL + bWPQ * cW * PQ) / 3.0;
  const double FmDTT = -( aW * DENSMOL / TEMP - bWPQ * cW2T / PQ) / 3.0;
  const double FmDRT = ( aW * DENSMOL / TEMP + bWPQ * cW1T) / 3.0;

  //  3) mix molecular + plasma via smooth switch YL/YH
  const double X   = Q4 * log(Q1 * rho + Q2 * T);
  const double X1R = Q4 * Q1 * rho   / (Q1 * rho + Q2 * T);
  const double X1T = Q4 * Q2 * T     / (Q1 * rho + Q2 * T);
  const double X2R = Q4 * Q1 * Q2 * rho * T / pow(Q1 * rho + Q2 * T, 2);
  const double X2T = X2R;
  const double X2RT= -X2R;

  const double YL  = FERMIF(X);
  const double YH  = 1.0 - YL;

  const double YH1X = YH * YL;
  const double YH2X = YH1X * (YL - YH);

  const double YH1R = YH1X * X1R;
  const double YH1T = YH1X * X1T;

  const double YH2R  = YH2X * pow(X1R, 2) + YH1X * X2R;
  const double YH2T  = YH2X * pow(X1T, 2) + YH1X * X2T;
  const double YH2RT = YH2X * X1R * X1T        + YH1X * X2RT;

  const double FNkTni = FNkTmol * YL + FNkTsi * YH;
  const double PnkTni = PnkTmol * YL + PnkTsi * YH + (FNkTsi - FNkTmol) * YH1R;
  const double UNkTni = UNkTmol * YL + UNkTsi * YH - (FNkTsi - FNkTmol) * YH1T;

  const double FDRR = FmDRR * YL + FsiDRR * YH + 2.0 * (PnkTsi - PnkTmol) * YH1R
    + (FNkTsi - FNkTmol) * YH2R;
  const double FDTT = FmDTT * YL + FsiDTT * YH - 2.0 * (UNkTsi - UNkTmol) * YH1T
    + (FNkTsi - FNkTmol) * YH2T;
  const double FDRT = FmDRT * YL + FsiDRT * YH + (PnkTsi - PnkTmol) * YH1T
    - (UNkTsi - UNkTmol) * YH1R + (FNkTsi - FNkTmol) * YH2RT;

  //  4) Ideal gas of pseudo‑molecules
  const double THLmol = sqrt(2 * pi / (18.0 * AUM * TEMP));
  const double FNkTid = (log(DENSMOL * pow(THLmol, 3)) - 1.0) / 3.0;
  const double PnkTid = C13;
  const double UNkTid = 0.5;

  //  5) total (non‑ideal + ideal)
  double PnkT,FNkT,UNkT,CV,CHIT,CHIR;
  FNkT = FNkTni + FNkTid;
  PnkT = PnkTni + PnkTid;		
  UNkT = UNkTni + UNkTid;

  CV   = UNkT - FDTT;				
  CHIR = FDRR / PnkT + 1.0;	
  CHIT = FDRT / PnkT + 1.0;	

  //  6) thermal corrections (high‑T & low‑T tweaks)
  const double TTC  = TEMP / TCRIT;
  const double TL2  = PC4 * TTC;
  const double ULB  = pow(TL2, 2) * sqrt(TL2);  // (TL2)^2.5
  const double ULB1 = 1.0 + ULB;
  const double FL   = log(ULB1 / ULB);
  const double UL   = 2.5 / ULB1;
  const double CVL  = UL * (1.0 - 1.5 * ULB) / ULB1;

  const double TTC2 = TTC * TTC;
  const double ULC1 = 1.0 + TTC2;
  const double FC   = (PC1 * log(ULC1 / TTC2) + PC2 * atan(TTC)) / TCRIT
    - PC3 / TEMP;
  const double UC   = ((2.0 * PC1 * TTC - PC2 * TTC2) / ULC1 - PC3) / TEMP;
  const double CVC  = 2.0 / TCRIT * (PC1 * (1.0 - TTC2) - PC2 * TTC) / (ULC1 * ULC1);

  FNkT += FL + FC - SBASE;
  UNkT += UL + UC;
  CV   += CVL + CVC;

  //  7) auxiliary conversions
  double nabla_ad = CHIT / (CV*CHIR/PnkT + sq(CHIT));				// d ln T / d ln P
  double dT_dP_S  = nabla_ad * 1E6 / (PnkT * TnkCONV * rho);
  // nabla_ad * T / P = nabla_ad * T / (PnkT * n k T) = nabla_ad * T / (PnkT * rho * NA *k / MImean * T)
  return dT_dP_S;
}

double H2OSC(double P, double T, double rho_guess)
// P in GPa
{
  if(!(rho_guess>0.5 && rho_guess<20))
    rho_guess = max(1 + 0.5*log(P/5E10),0.1);
  return density_solver(P,T,PH2OSC,rho_guess);
}

double dTdP_S_H2OSC(double P, double T, double &rho_guess)
{
  rho_guess = density_solver(P,T,PH2OSC,rho_guess);
  return dTdP_S_H2O_of_rho(rho_guess, T);
}
	

//  H2O SC Mazevet Helper Functions
//  ELECNR – ideal, non‑relativistic electron Fermi gas
void ELECNR(double DENSE, double TEMP,
            double &CHI, double &FEid, double &PEid, double &UEid,
            double &SEid, double &CVE, double &CHITE, double &CHIRE)
{
    static const double SQPI = sqrt(pi);
    const double CLE  = sqrt(2 * pi/ TEMP);
    const double FDENS= SQPI * std::pow(CLE, 3) * DENSE / 4.0;

    double X, XDF, XDFF;
    FINVER(FDENS, 1, X, XDF, XDFF);     // inverse Fermi integral
    CHI = X;

    double F12, F32;
    FERINT(X, 1, F12);
    FERINT(X, 2, F32);

    UEid = F32 / FDENS;
    PEid = UEid / 1.5;
    FEid = CHI - PEid;
    SEid = UEid - FEid;

    const double XDF32 = 1.5 * FDENS * XDF;
    CHITE = 2.5 - XDF32 / PEid;
    CHIRE = XDF32 * F12 / F32;
    CVE   = 1.5 * PEid * CHITE;
}

//  FERMIF: simple logistic Fermi function 1/(exp(x)+1)
double FERMIF(double X)
{
    if (X > 40.0)  return 0.0;
    if (X < -40.0) return 1.0;
    return 1.0 / (exp(X) + 1.0);
}


//  FERINT – Antia (1993) rational fits to F_q(x)
//            q = N-1/2 with N=0..3 (−½, ½, 3/2, 5/2)
void FERINT(double X, int N, double &F)
{
    if (N < 0 || N > 3) throw std::invalid_argument("FERINT: N out of range");

    // coefficient tables (flattened)
    static const double A[8][4] = {
      {1.71446374704454e7,  5.75834152995465e6,  4.32326386604283e4,  6.61606300631656e4},
      {3.88148302324068e7,  1.30964880355883e7,  8.55472308218786e4,  1.20132462801652e5},
      {3.16743385304962e7,  1.07608632249013e7,  5.95275291210962e4,  7.67259953168120e4},
      {1.14587609192151e7,  3.93536421893014e6,  1.77294861572005e4,  2.10427138842443e4},
      {1.83696370756153e6,  6.42493233715640e5,  2.21876607796460e3,  2.44325236813275e3},
      {1.14980998186874e5,  4.16031909245777e4,  9.90562948053293e1,  1.02589947781696e2},
      {1.98276889924768e3,  7.77238678539648e2,  1.0,                  1.0                },
      {1.0,                 1.0,                 0.0,                  0.0               }
    };

    static const double B[8][4] = {
      {9.67282587452899e6,  6.49759261942269e6,  3.25218725353467e4,  1.99078071053871e4},
      {2.87386436731785e7,  1.70750501625775e7,  7.01022511904373e4,  3.79076097261066e4},
      {3.26070130734158e7,  1.69288134856160e7,  5.50859144223638e4,  2.60117136841197e4},
      {1.77657027846367e7,  7.95192647756086e6,  1.95942074576400e4,  7.97584657659364e3},
      {4.81648022267831e6,  1.83167424554505e6,  3.20803912586318e3,  1.10886130159658e3},
      {6.13709569333207e5,  1.95155948326832e5,  2.20853967067789e2,  6.35483623268093e1},
      {3.13595854332114e4,  8.17922106644547e3,  5.05580641737527e0,  1.16951072617142e0},
      {4.35061725080755e2,  9.02129136642157e1,  1.99507945223266e-2, 3.31482978240026e-3}
    };

    // polynomial degrees per Antia table
    static const int LA[4] = {7,7,6,6};
    static const double C[12][4] = {
      {-4.46620341924942e-15,  4.85378381173415e-14,  2.80452693148553e-13,  8.42667076131315e-12},
      {-1.58654991146236e-12,  1.64429113030738e-11,  8.60096863656367e-11,  2.31618876821567e-09},
      {-4.44467627042232e-10,  3.76794942277806e-09,  1.62974620742993e-08,  3.54323824923987e-07},
      {-6.84738791621745e-08,  4.69233883900644e-07,  1.63598843752050e-06,  2.77981736000034e-05},
      {-6.64932238528105e-06,  3.40679845803144e-05,  9.12915407846722e-05,  1.14008027400645e-03},
      {-3.69976170193942e-04,  1.32212995937796e-03,  2.62988766922117e-03,  2.32779790773633e-02},
      {-1.12295393687006e-02,  2.60768398973913e-02,  3.85682997219346e-02,  2.39564845938301e-01},
      {-1.60926102124442e-01,  2.48653216266227e-01,  2.78383256609605e-01,  1.24415366126179e+00},
      {-8.52408612877447e-01,  1.08037861921488e+00,  9.02250179334496e-01,  3.18831203950106e+00},
      {-7.45519953763928e-01,  1.91247528779676e+00,  1.0,                 3.42040216997894e+00},
      { 2.98435207466372e+00,  1.0,                 0.0,                 1.0},
      { 1.0,                  0.0,                 0.0,                 0.0}
    };

    static const double D[12][4] = {
      {-2.23310170962369e-15,  7.28067571760518e-14,  7.01131732871184e-13,  2.94933476646033e-11},
      {-7.94193282071464e-13,  2.45745452167585e-11,  2.10699282897576e-10,  7.68215783076936e-09},
      {-2.22564376956228e-10,  5.62152894375277e-09,  3.94452010378723e-08,  1.12919616415947e-06},
      {-3.43299431079845e-08,  6.96888634549649e-07,  3.84703231868724e-06,  8.09451165406274e-05},
      {-3.33919612678907e-06,  5.02360015186394e-05,  2.04569943213216e-04,  2.81111224925648e-03},
      {-1.86432212187088e-04,  1.92040136756592e-03,  5.31999109566385e-03,  3.99937801931919e-02},
      {-5.69764436880529e-03,  3.66887808001874e-02,  6.39899717779153e-02,  2.27132567866839e-01},
      {-8.34904593067194e-02,  3.24095226486468e-01,  3.14236143831882e-01,  5.31886045222680e-01},
      {-4.78770844009440e-01,  1.16434871200131e+00,  4.70252591891375e-01,  3.70866321410385e-01},
      {-4.99759250374148e-01,  1.34981244060549e+00, -2.15540156936373e-02,  2.27326643192516e-02},
      { 1.86795964993052e+00,  2.01311836975930e-01,  2.34829436438087e-03,  0.0},
      { 4.16485970495288e-01, -2.14562434782759e-02,  0.0,                 0.0}
    };

    static const int LC[4] = {11,10,9,1};

    // choose branch (X<2 or >=2)
    if (X < 2.0) {
        double t = exp(X);
        double up=0.0, down=0.0;
        for (int i=LA[N]; i>=0; --i) {
            up   = up   * t + A[i][N];
            down = down * t + B[i][N];
        }
        F = t * up / down;
    } else {
        double t = 1.0 / (X*X);
        double up=0.0, down=0.0;
        for (int i=LC[N]; i>=0; --i) {
            up   = up   * t + C[i][N];
            down = down * t + D[i][N];
        }
        F = sqrt(X) * pow(X, N) * up / down;
    }
}

//  FINVER – Antia (1993) fits to inverse Fermi integrals X_q(f)
void FINVER(double F, int N, double &X, double &XDF, double &XDFF)
{
    if (N < 0 || N > 3) throw std::invalid_argument("FINVER: N out of range");
    if (F <= 0.0)       throw std::invalid_argument("FINVER: F must be >0");

    // shortened coefficient tables (same as Fortran)
    static const double A[6][4] = {
      {-1.570044577033e4, 1.999266880833e4, 1.715627994191e2, 2.138969250409e2},
      { 1.001958278442e4, 5.702479099336e3, 1.125926232897e2, 3.539903493971e1},
      {-2.805343454951e3, 6.610132843877e2, 2.056296753055e1, 1.0},
      { 4.121170498099e2, 3.818838129486e1, 1.0,                0.0},
      {-3.174780572961e1, 1.0,                0.0,                0.0},
      { 1.0,               0.0,                0.0,                0.0}
    };

    static const double B[7][4] = {
      {-2.782831558471e4,  1.771804140488e4,  2.280653583157e2,  7.10854551271e2},
      { 2.886114034012e4, -2.014785161019e3,  1.193456203021e2,  9.873746988121e1},
      {-1.274243093149e4,  9.130355392717e1,  1.167743113540e1,  1.067755522895e0},
      { 3.063252215963e3, -1.670718177489e0, -3.226808804038e-1,-1.182798726503e-2},
      {-4.225615045074e2,  0.0,              3.519268762788e-3,  0.0},
      { 3.168918168284e1,  0.0,              0.0,                0.0},
      {-1.008561571363e0,  0.0,              0.0,                0.0}
    };
    static const int LA[4] = {5,4,3,2};
    static const int LB[4] = {6,3,4,3};

    // C & D tables truncated (only first 7 rows used like in Antia fits)
    static const double C[7][4] = {
      { 2.206779160034e-8,  -1.277060388085e-2,-6.321828169799e-3, -3.312041011227e-2},
      {-1.437701234283e-6,   7.187946804945e-2,-2.183147266896e-2,  1.315763372315e-1},
      { 6.103116850636e-5,  -4.262314235106e-1,-1.057562799320e-1, -4.820942898296e-1},
      {-1.169411057416e-3,   4.997559426872e-1,-4.657944387545e-1,  5.099038074944e-1},
      { 1.814141021608e-2,  -1.285579118012e+0,-5.951932864088e-1,  5.495613498630e-1},
      {-9.588603457639e-2,  -3.930805454272e-1, 3.684471177100e-1,-1.498867562255e+0},
      { 1.0,                1.0,               1.0,                1.0}
    };

    static const double D[7][4] = {
      { 8.827116613576e-8,  -9.745794806288e-3,-4.381942605018e-3, -2.315515517515e-2},
      {-5.750804196059e-6,   5.485432756838e-2,-1.513236504100e-2,  9.198776585252e-2},
      { 2.429627688357e-4,  -3.299466243260e-1,-7.850001283886e-2, -3.835879295548e-1},
      {-4.601959491394e-3,   4.077841975923e-1,-3.407561772612e-1,  5.415026856351e-1},
      { 6.932122275919e-2,  -1.145531476975e+0,-5.074812565486e-1, -3.847241692193e-1},
      {-3.217372489776e-1,  -6.067091689181e-2,-1.387107009074e-1,  3.739781456585e-2},
      { 3.124344749296e+0,   0.0,               0.0,               -3.008504449098e-2}
    };

    static const int LD[4] = {6,5,5,6};

    // ------- branch F<4 ------------------------------------------------------
    if (F < 4.0) {
        double t = F;
        double up=0.0, up1=0.0, up2=0.0;
        double dn=0.0, dn1=0.0, dn2=0.0;
        for (int i=LA[N]; i>=0; --i) {
            up  = up  * t + A[i][N];
            if (i>=1) up1 = up1 * t + A[i][N] * i;
            if (i>=2) up2 = up2 * t + A[i][N] * i * (i-1);
        }
        for (int i=LB[N]; i>=0; --i) {
            dn  = dn  * t + B[i][N];
            if (i>=1) dn1 = dn1 * t + B[i][N] * i;
            if (i>=2) dn2 = dn2 * t + B[i][N] * i * (i-1);
        }
        X    = log(t * up / dn);
        XDF  = 1.0/t + up1/up - dn1/dn;
        XDFF = -1.0/(t*t) + up2/up - pow(up1/up,2) - dn2/dn + pow(dn1/dn,2);
        return;
    }

    // ------- branch F≥4 ------------------------------------------------------
    const double P = -1.0 / (0.5 + N);
    double t   = pow(F, P);
    double t1  = P * t / F;
    double t2  = P * (P - 1.0) * t / (F*F);
    double up=0.0, up1=0.0, up2=0.0;
    double dn=0.0, dn1=0.0, dn2=0.0;
    for (int i=6; i>=0; --i) {
        up  = up  * t + C[i][N];
        if (i>=1) up1 = up1 * t + C[i][N] * i;
        if (i>=2) up2 = up2 * t + C[i][N] * i * (i-1);
    }
    for (int i=LD[N]; i>=0; --i) {
        dn  = dn  * t + D[i][N];
        if (i>=1) dn1 = dn1 * t + D[i][N] * i;
        if (i>=2) dn2 = dn2 * t + D[i][N] * i * (i-1);
    }
    const double R  = up / dn;
    const double R1 = (up1 - up * dn1 / dn) / dn;
    const double R2 = (up2 - (2.0 * up1 * dn1 + up * dn2) / dn + 2.0 * up * pow(dn1/dn,2)) / dn;

    X    = R / t;
    double RT = (R1 - R / t) / t;
    XDF  = t1 * RT;
    XDFF = t2 * RT + pow(t1,2) * (R2 - 2.0*RT) / t;
}


