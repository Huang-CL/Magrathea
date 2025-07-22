#ifndef EOSLIST_H_
#define EOSLIST_H_

#include "EOS.h"

extern EOS *Fe_liquid, *Fe_liquid2, *Fe_fcc, *Fe_bcc, *Fe_hcp, *Fe_hcp2, *Fe_Seager, *Fe_hcp3, *Fe_7Si, *Fe_15Si, *Fe_Dummy, *Si_Pv_Shim, *Si_Pv, *Si_PPv, *Si_PPv_Sakai, *Si_PREM, *Si_BM2fit, *Si_Seager, *Si_Dummy, *Si_liquid, *Si_Liquid_Wolf, *Fo, *Wds, *Rwd, *Akm, *Pv_Doro, *PPv_Doro, *Fo_Sotin, *En, *Mw, *Ice_Seager, *H2O_AQUA, *H2O_SeaFreeze, *Water_ExoPlex, *Water, *Water_sc_Mazevet, *Water_sc_dummy, *IceIh, *IceIh_ExoPlex, *IceVI_ExoPlex, *IceVI_Bezacier, *IceVII_ExoPlex, *IceVII_Bezacier, *IceVII, *IceVIIp, *IceVII_FFH2004, *IceVII_FFH2004fit, *IceVII_FFH2004BM, *IceVII_Fei, *IceVII_FFH2004T, *IceX_HS, *IceX, *IceZeng2013FFH, *IceZeng2013FMNR, *Ice_Dummy, *Gas, *Gas_iso, *Gas_hhe, *watervapor, *vdW_H2, *vdW_He, *vdW_H2O, *vdW_CH4, *vdW_NH3, *vdW_CO2, *Gold, *Plat, *Graph, *Diam;

double dTdP_Si_Dummy (double P, double T);
// A temperature gradient that equals to the melting curve. Guarantee the temperature won't drop below the melting curve. 

double H2OSC(double P,double T, double rho_guess);

double dTdP_S_H2OSC(double P, double T, double &rho_guess);

void ELECNR(double DENSE, double TEMP,
            double &CHI, double &FEid, double &PEid, double &UEid,
            double &SEid, double &CVE, double &CHITE, double &CHIRE);

double FERMIF(double X);
void FERINT(double X, int N, double &F);
void FINVER(double F, int N, double &X, double &XDF, double &XDFF);


#endif
