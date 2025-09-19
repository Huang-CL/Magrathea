#include "define.h"
#include "EOS.h"
#include "EOSmodify.h"
#include "phase.h"
#include "hydro.h"
#include "compfind.h"
#include "parser.h"

struct timeval start_time, end_time;

//Global Parameter Defaults
double rho_eps_rel = 1E-11;	// relative error tolerance of the density solver
double T_eps_rel = 1E-11;	// relative error tolerance of the temperature solver for the adiabatic profile
double ode_eps_rel0 = 1E-7; // relative error tolerance for first round ode integrator (1E-7) used in mode 0
double ode_eps_rel1 = 1E-10; // relative error tolerance for second round ode integrator (1E-10) used in mode 0
int fit_iter;
double R_eps_rel = 2E-5; // relative error tolerance in mode 0 first round radius determination (5E-5).  Should be around sqrt(ode_eps_rel0).
double ode_eps_rel2 = 1E-10; // relative error tolerance for ode integrator (1E-10) used in mode 1
double P_eps_rel = 1E-10;	// relative error tolerance in mode 1 central pressure determination (1E-10).  Should not be more restrict than ode_eps_rel.
double fit_eps_rel = 1E-4; // the relative error tolerance at the fitting point in mode 0 round 2 (1E-4). Should be four orders of magnitudes larger than ode_eps_rel1.
vector<double> ave_rho = {15, 5, 2, 1E-3};// Assuming the density of the core is 15, mantle is 5, water is 2, and gas is 1E-3.
bool verbose = false;		  // Whether print warnings.
double P_surface = 1E5;		  // The pressure level that the broad band optical transit radius probes. (in microbar)
int count_shoot = 0;			  // used to count the total number of shootings per each solution
int count_step = 0;			  // used to count the sum of integral steps in all shooting results.


int main(int argc, char* argv[])
{
  gsl_set_error_handler_off();  //Dismiss the abortion from execution, which is designed for code testing.
  hydro* planet;
  ifstream fin;
  int n_settings;

  //Initialize all the parameters, will throw error later if missing
  int input_mode=0;
  string core_phasedgm="Fe_default";
  string mantle_phasedgm="Si_default";
  string hydro_phasedgm="water_default";
  string atm_phasedgm="gas_default";
  vector<double> Mcomp={0,0,0,0}; 
  vector<double> Tgap = {0,0,0,0};
  string outputfile=" "; 
  int layer_index=0;
  float mass_frac=0;
  float min_mass=0;
  float max_mass=0;
  float step_mass=0;
  string inputfile=" ";
  int solver=1;
  int findlayer=1;
  float minPMR=0;
  float maxPMR=0;
  float step=0;
  vector<int> layers={1,1,0,0};
  float rerr=0;
  float MassPrior=0;
  float MUncPrior=0;
  float RadPrior=0;
  float RUncPrior=0;
  int numlayers=2;
  int numchains=1;
  int chainsteps=10;

  Water_sc_Mazevet -> modify_dTdP(dTdP_S_H2OSC);

  //Parse input file
  try{
    Settings options(argv[1]);  // Will throw an exception if the file couldn't be opened.
    n_settings = options.LoadSettings(); // Load settings from the file and return the number of options found.
    std::cout << "Loaded " << n_settings << " settings.\n";
    //Error Tolerances
    rho_eps_rel = options.GetOptionDouble("rho_eps_rel");
    T_eps_rel = options.GetOptionDouble("T_eps_rel");	
    ode_eps_rel0 = options.GetOptionDouble("ode_eps_rel0");
    ode_eps_rel1 = options.GetOptionDouble("ode_eps_rel1");
    R_eps_rel = options.GetOptionDouble("R_eps_rel");
    ode_eps_rel2 = options.GetOptionDouble("ode_eps_rel2");
    P_eps_rel = options.GetOptionDouble("P_eps_rel");
    fit_eps_rel = options.GetOptionDouble("fit_eps_rel");
    //Global Run Options
    verbose = options.GetOptionBool("verbose");		  
    P_surface = options.GetOptionDouble("P_surface");
    vector<double> ave_rho = {options.GetOptionDouble("ave_rho_core"), options.GetOptionDouble("ave_rho_mantle"), options.GetOptionDouble("ave_rho_hydro"), options.GetOptionDouble("ave_rho_atm")};
    //Mode Options
    // Choose between the 8 available input_mode values:
    // 0: regular solver, 1: temperature-free solver, 2: two-layer solver, 
    // 3: bulk input mode with regular solver
    // 4: composition finder, finds third layer mass to match a mass and radius measurement
    // 5: modify a built-in EOS on they fly, 
    // 6: iterate over EOS modifications with two-layer solver, 5: iterate over EOS with regular solver
    input_mode = options.GetOptionDouble("input_mode");
    //Get Phase Diagrams for modes which require
    if (input_mode==0 or input_mode==3 or input_mode==4 or input_mode==7 or input_mode==8)
    {
      core_phasedgm=options.GetOptionString("core_phasedgm");
      mantle_phasedgm=options.GetOptionString("mantle_phasedgm");
      hydro_phasedgm=options.GetOptionString("hydro_phasedgm");
      atm_phasedgm=options.GetOptionString("atm_phasedgm");
    }
    //Get parameters specific to each mode
    switch (input_mode){
    case 0:
      Mcomp[0]=options.GetOptionDouble("mass_of_core");
      Mcomp[1]=options.GetOptionDouble("mass_of_mantle");
      Mcomp[2]=options.GetOptionDouble("mass_of_hydro");
      Mcomp[3]=options.GetOptionDouble("mass_of_atm");
      Tgap[0]=options.GetOptionDouble("temp_jump_3"); 
      Tgap[1]=options.GetOptionDouble("temp_jump_2");
      Tgap[2]=options.GetOptionDouble("temp_jump_1");
      Tgap[3]=options.GetOptionDouble("surface_temp");
      outputfile=options.GetOptionString("output_file");
      break;
    case 1:
      Mcomp[0]=options.GetOptionDouble("mass_of_core");
      Mcomp[1]=options.GetOptionDouble("mass_of_mantle");
      Mcomp[2]=options.GetOptionDouble("mass_of_hydro");
      outputfile=options.GetOptionString("output_file");
      break;
    case 2:
      layer_index=options.GetOptionDouble("layer_index");
      mass_frac=options.GetOptionDouble("mass_frac");
      min_mass=options.GetOptionDouble("min_mass");
      max_mass=options.GetOptionDouble("max_mass");
      step_mass=options.GetOptionDouble("step_mass");
      break;
    case 3:
      inputfile=options.GetOptionString("input_file");
      solver=options.GetOptionDouble("solver");
      Tgap[0]=options.GetOptionDouble("temp_jump_3"); 
      Tgap[1]=options.GetOptionDouble("temp_jump_2");
      Tgap[2]=options.GetOptionDouble("temp_jump_1");
      Tgap[3]=options.GetOptionDouble("surface_temp");
      outputfile=options.GetOptionString("output_file");      
      break;
    case 4:
      inputfile=options.GetOptionString("input_file");
      outputfile=options.GetOptionString("output_file");
      findlayer=options.GetOptionDouble("find_layer");
      layers[options.GetOptionDouble("layer_inner")-1]=1;
      layers[options.GetOptionDouble("layer_outer")-1]=1;
      minPMR=options.GetOptionDouble("PMR_min");
      maxPMR=options.GetOptionDouble("PMR_max");
      step=options.GetOptionDouble("PMR_step");
      rerr=options.GetOptionDouble("R_error");
      Tgap[0]=options.GetOptionDouble("temp_jump_3"); 
      Tgap[1]=options.GetOptionDouble("temp_jump_2");
      Tgap[2]=options.GetOptionDouble("temp_jump_1");
      Tgap[3]=options.GetOptionDouble("surface_temp");
      break;
    case 5:   
      break; //must be run from main.cpp
    case 6:
      inputfile=options.GetOptionString("input_file");
      outputfile=options.GetOptionString("output_file");
      layer_index=options.GetOptionDouble("layer_index");
      mass_frac=options.GetOptionDouble("mass_frac");   
      break;
    case 7:
      inputfile=options.GetOptionString("input_file");
      outputfile=options.GetOptionString("output_file");
      Mcomp[0]=options.GetOptionDouble("mass_of_core");
      Mcomp[1]=options.GetOptionDouble("mass_of_mantle");
      Mcomp[2]=options.GetOptionDouble("mass_of_hydro");
      Mcomp[3]=options.GetOptionDouble("mass_of_atm");
      Tgap[0]=options.GetOptionDouble("temp_jump_3"); 
      Tgap[1]=options.GetOptionDouble("temp_jump_2");
      Tgap[2]=options.GetOptionDouble("temp_jump_1");
      Tgap[3]=options.GetOptionDouble("surface_temp");
      break;
    case 8:
      MassPrior=options.GetOptionDouble("mass_prior");
      MUncPrior=options.GetOptionDouble("mass_unc");
      RadPrior=options.GetOptionDouble("radius_prior");
      RUncPrior=options.GetOptionDouble("radius_unc");
      numlayers=options.GetOptionDouble("num_layers");
      numchains=options.GetOptionDouble("num_chains");
      chainsteps=options.GetOptionDouble("chain_steps");
      outputfile=options.GetOptionString("output_file");
      Tgap[0]=options.GetOptionDouble("temp_jump_3"); 
      Tgap[1]=options.GetOptionDouble("temp_jump_2");
      Tgap[2]=options.GetOptionDouble("temp_jump_1");
      Tgap[3]=options.GetOptionDouble("surface_temp");
      break;
    }
  } catch (SettingsParserException& e) { // This will be triggered if an exception is thrown above.
    std::cout << "\nException: " << e.what() << "\n"; // e.what() prints the exception message.
    return EXIT_FAILURE;
  }

  //Set Phase Diagrams from string in config
  vector<PhaseDgm> Comp = {core, mant, water, atm};
  if (core_phasedgm=="Fe_default")
    Comp[0]=core;
  else if (core_phasedgm=="Fe_fccbcc")
    Comp[0]=core1;
  else
    cout<<"core_phasedgm does not exist, using default"<<endl;
  if (mantle_phasedgm=="Si_default")
    Comp[1]=mant;
  else if (mantle_phasedgm=="Si_simple")
    Comp[1]=mant1;
  else if (mantle_phasedgm=="PREM")
    Comp[1]=mant2;
  else if (mantle_phasedgm=="C_simple")
    Comp[1]=mant3;
  else
    cout<<"mant_phasedgm does not exist, using default"<<endl;
  if (hydro_phasedgm=="water_default")
    Comp[2]=water; 
  else if (hydro_phasedgm=="water_tabulated")
    Comp[2]=water1;
  else if (hydro_phasedgm=="water_legacy")
    Comp[2]=water2;
  else
    cout<<"hydro_phasedgm does not exist, using default"<<endl;
  if (atm_phasedgm=="gas_default")
    Comp[3]=atm;
  else if (atm_phasedgm=="HHe_tabulated")
    Comp[3]=atm1;
  else
    cout<<"atm_phasedgm does not exist, using default"<<endl;


  //START main code for the 9 input modes
	// 0: regular solver, 1: temperature-free solver, 2: two-layer solver, 
	// 3: bulk input mode with regular solver
	// 4: composition finder, finds third layer mass to match a mass and radius measurement
	// 5: modify a built-in EOS on they fly, 
	// 6: iterate over EOS modifications with two-layer solver, 7: iterate over EOS with regular solver
  // 8: MCMC mass fraction finder

  if (input_mode == 0)
  {
    double deltat;
    gettimeofday(&start_time,NULL);
    planet=fitting_method(Comp, Mcomp, Tgap, ave_rho, P_surface, false);
    cout<<"# of shots "<<count_shoot<<", # of total steps "<<count_step<<endl;
    gettimeofday(&end_time, NULL);
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;
    cout<<"run time "<<deltat<<'s'<<endl;
    if (!planet)
    {
      for (unsigned int i=0; i < Mcomp.size(); i++)
        cout<<Mcomp[i]<<", ";
      cout<<"\t No solution found."<<endl;
    }
    else
    {
      planet->print(outputfile, false); // Save the result in an asc file with this name.
      cout<<"Planet saved: "<<outputfile<<endl;
    }
    delete planet;
  }  

  else if(input_mode == 1)
  {
    double deltat;
    gettimeofday(&start_time,NULL);
    planet=getmass(Mcomp[0],Mcomp[1],Mcomp[2],P_surface);
    // Mass in Earth Masses of Core, Mantle, Hydrosphere
    gettimeofday(&end_time, NULL);
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;
    cout<<"run time "<<deltat<<'s'<<endl;
    if (planet)
    {
      planet->print(outputfile);
      cout<<"Planet saved: "<<outputfile<<endl;
    }
    // Save the result in an asc file with this name.
    delete planet;
  }
  
  else if(input_mode == 2)
  {
    vector<double> Mp,Rp;
    double deltat;
    for (float value = min_mass; value <= max_mass; value += step_mass) {
      Mp.push_back(value);
    }
    gettimeofday(&start_time,NULL);
    twolayer(layer_index,mass_frac,Mp,Rp,P_surface,true);
    for(int i=0; i < int(Mp.size()); i++)
      cout<<Mp[i]<<" \t"<<Rp[i]<<endl;
    gettimeofday(&end_time, NULL);
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;
    cout<<"run time "<<deltat<<'s'<<endl;
  }

  else if(input_mode == 3)
  {
    multiplanet(Comp, Tgap, solver, ave_rho, P_surface, false, inputfile, outputfile);
  }

  else if(input_mode == 4)
  {
    //Multi-threading is commented out by default for easy install
    //Must uncomment Line 2 in Makefile and all occurences of "pragma..." in comfind.cpp 
    int num_threads=0; 
    
    compfinder(Comp,findlayer,layers,minPMR,maxPMR,step,rerr,num_threads,Tgap,ave_rho,P_surface,false,inputfile,outputfile);
  }

  else if(input_mode == 5)
  {
    vector<double> Mp,Rp;
    double fraction = 0;
    
    int nPhases=3;
    EOS **highEOSs = new EOS*[nPhases];
    for(int i=0; i<nPhases; i++)
      highEOSs[i] = new EOS();
    highEOSs[0] -> setphasename("Ice VII");
    highEOSs[1] -> setphasename("Ice VII'");
    highEOSs[2] -> setphasename("Ice X");

    double IceVII_Grande[5][2]={{0.,0.},{1.,13.62},{2.,10.48},{3.,4.},{5.,18.01528}};
    double IceVIIp_Grande[5][2]={{0.,0.},{1.,12.32},{2.,22.36},{3.,4.},{5.,18.01528}};
    double IceX_Grande[5][2]={{0.,0.},{1.,11.0},{2.,38.6},{3.,4.},{5.,18.01528}};
    double start_P[2]={4.76,34.14};
    highEOSs[0] -> modifyEOS(IceVII_Grande,5);
    highEOSs[1] -> modifyEOS(IceVIIp_Grande,5);
    highEOSs[2] -> modifyEOS(IceX_Grande,5);

    water.set_phase_highP(nPhases, start_P, highEOSs);
    planet = getmass(0,0,0.0999968,P_surface);
  
    if (planet)
      planet->print("./result/Structure.txt");
  
    twolayer(0,fraction,Mp,Rp,P_surface);

    for(int i=0; i < int(Mp.size()); i++)
      cout<<Mp[i]<<" \t"<<Rp[i]<<endl;

    for(int i=0; i<nPhases; i++)
      delete highEOSs[i];
    delete[] highEOSs;
  }  

  else if(input_mode == 6)
  {
    double deltat;    
    gettimeofday(&start_time,NULL);
    twolayer(layer_index,mass_frac,P_surface,1,inputfile, outputfile);
    gettimeofday(&end_time, NULL);
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;
    cout<<"run time "<<deltat<<'s'<<endl;
  }
  
  else if(input_mode == 7)
  {
    double deltat;    
    gettimeofday(&start_time,NULL);
    fullmodel(Comp,Mcomp,Tgap,ave_rho,P_surface,false,2,inputfile,outputfile);
    gettimeofday(&end_time, NULL);
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;
    cout<<"run time "<<deltat<<'s'<<endl;
  }

  else if(input_mode == 8)
  {
    double deltat;    
    gettimeofday(&start_time,NULL);
    mcmcsample(Comp,MassPrior,MUncPrior,RadPrior,RUncPrior,Tgap,ave_rho,P_surface,false,outputfile,numchains,chainsteps,numlayers);
    gettimeofday(&end_time, NULL);
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;
    cout<<"run time "<<deltat<<'s'<<endl;
  }

  // ============================
  
  delete Fe_liquid;
  delete Fe_liquid2;
  delete Fe_fcc;
  delete Fe_bcc;
  delete Fe_hcp;
  delete Fe_hcp2;
  delete Fe_hcp3;
  delete Fe_Dummy;
  delete Fe_7Si;
  delete Fe_15Si;
  delete Fe_Seager;
  delete Si_Pv_Shim;
  delete Si_Pv;
  delete Si_PPv;
  delete Si_PPv_Sakai;
  delete Si_PREM;
  delete Si_BM2fit;
  delete Si_Seager;
  delete Si_liquid;
  delete Si_Liquid_Wolf;
  delete Si_Dummy;
  delete Fo;
  delete Wds;
  delete Rwd;
  delete Akm;
  delete Pv_Doro;
  delete PPv_Doro;
  delete Fo_Sotin;
  delete En;
  delete Mw;
  delete Ice_Seager;
  delete H2O_AQUA;
  delete H2O_SeaFreeze;
  delete Water_SF;
  delete Water_Brown;
  delete Water_IAPWS;
  delete Water_ExoPlex;
  delete Water;
  delete Water_sc_Mazevet;
  delete IceIh_SF;
  delete IceII_SF;
  delete IceIII_SF;
  delete IceV_SF;
  delete IceVI_SF;
  delete IceIh;
  delete IceIh_ExoPlex;
  delete IceVI_ExoPlex;
  delete IceVI_Bezacier;
  delete IceVII_ExoPlex;
  delete IceVII_Bezacier;
  delete IceVII;
  delete IceVIIp;
  delete IceVII_FFH2004;
  delete IceVII_FFH2004fit;
  delete IceVII_Fei;
  delete IceVII_FFH2004BM;
  delete IceX_HS;
  delete IceX;
  delete IceZeng2013FFH;
  delete IceZeng2013FMNR;
  delete Ice_Dummy;
  delete Gas;
  delete Gas_iso;
  delete Gas_hhe;
  delete watervapor;
  delete Water_Vap_IAPWS;
  delete Gold;
  delete Plat;
  delete Graph;
  delete Diam;
  delete vdW_H2;
  delete vdW_He;
  delete vdW_H2O;
  delete vdW_CH4;
  delete vdW_NH3;
  delete vdW_CO2; 

  return 0;
}
