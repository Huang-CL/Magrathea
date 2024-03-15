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

  string core_phasedgm="Fe_default";
  string mantle_phasedgm="Si_default";
  string hydro_phasedgm="water_default";
  string atm_phasedgm="gas_default";
  int input_mode=0;
  vector<double> Mcomp={0,0,0,0}; 
  vector<double> Tgap = {0,0,0,0};
  string outputfile=" ";  

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
    input_mode= options.GetOptionDouble("input_mode");
    //Get Phase Diagrams for modes which require
    if (input_mode==0)
    {
      core_phasedgm=options.GetOptionString("core_phasedgm");
      mantle_phasedgm=options.GetOptionString("mantle_phasedgm");
      hydro_phasedgm=options.GetOptionString("hydro_phasedgm");
      atm_phasedgm=options.GetOptionString("atm_phasedgm");
    }
    switch (input_mode){
      case 0:
        Mcomp[0]=options.GetOptionDouble("mass_of_core");
        Mcomp[1]=options.GetOptionDouble("mass_of_mantle");
        Mcomp[2]=options.GetOptionDouble("mass_of_hydro");
        Mcomp[3]=options.GetOptionDouble("mass_of_atm");
        Tgap[0]=options.GetOptionDouble("temp_jump_3"); 
        Tgap[1]=options.GetOptionDouble("temp_jump_2");
        Tgap[2]=options.GetOptionDouble("temp_jump_2");
        Tgap[3]=options.GetOptionDouble("surface_temp");
	outputfile=options.GetOptionString("output");
	break;
      case 1:
        Mcomp[0]=options.GetOptionDouble("mass_of_core");
        Mcomp[1]=options.GetOptionDouble("mass_of_mantle");
        Mcomp[2]=options.GetOptionDouble("mass_of_hydro");
        outputfile=options.GetOptionString("output");
        break;
    }
	
   } catch (SettingsParserException& e) { // This will be triggered if an exception is thrown above.
       std::cout << "\nException: " << e.what() << "\n"; // e.what() prints the exception message.
       return EXIT_FAILURE;
   }

  //Set Phase Diagrams
  vector<PhaseDgm> Comp = {core, mant, water, atm};
  if (input_mode==0){
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
    else
      cout<<"mant_phasedgm does not exist, using default"<<endl;
    if (hydro_phasedgm=="water_default")
      Comp[2]=water; 
    else if (hydro_phasedgm=="water_tabulated")
      Comp[2]=water1;
    else
      cout<<"hydro_phasedgm does not exist, using default"<<endl;
    if (atm_phasedgm=="gas_default")
      Comp[3]=atm;
    if (atm_phasedgm=="HHe_tabulated")
      Comp[3]=atm1;
    else
      cout<<"atm_phasedgm does not exist, using default"<<endl;
  }	

  if (input_mode == 0)
  {
    planet=fitting_method(Comp, Mcomp, Tgap, ave_rho, P_surface, false);
    cout<<"# of shots "<<count_shoot<<", # of total steps "<<count_step<<endl;
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
    planet=getmass(Mcomp[0],Mcomp[1],Mcomp[2],P_surface);
    // Mass in Earth Masses of Core, Mantle, Hydrosphere
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
    gettimeofday(&start_time,NULL);

    twolayer(0,0,Mp,Rp,P_surface,true);

    for(int i=0; i < int(Mp.size()); i++)
      cout<<Mp[i]<<" \t"<<Rp[i]<<endl;
    gettimeofday(&end_time, NULL);
  
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;

    cout<<"running time "<<deltat<<'s'<<endl;
  }

  else if(input_mode == 3)
  {
    vector<PhaseDgm> Comp = {core,mant,water,atm};
    vector<double> Tgap = {0,0,0,300}; // Using full temperature solver, Temperature gap between each layer and surface temperature.

    string infilename="./input/inputcore.txt";
    string outfilename="./result/coreplanets.txt";

    int solver=1; //Calculation mode. Mode 2 is twice faster, but only applicable for planet with no gas layer and whose temperature effect is not important

    
    multiplanet(Comp, Tgap, solver, ave_rho, P_surface, false, infilename, outfilename);
    

  }

  else if(input_mode == 4)
  {
    vector<PhaseDgm> Comp = {core, mant, water, atm};
    vector<double> Tgap = {0,0,0,300}; // Using full temperature solver, Temperature gap between each layer and surface temperature.

    string infilename="./input/inputplanetMR.txt";
    string outfilename="./result/outputcompfindplanet.txt";

    //Solver will hold 2 layers in constant Partial Mass Ratio (PMR) and find the mass of 3rd layer: 
    //PMR(%) is OMF/(IMF+OMF)*100, where OMG is outer-mass fraction, IMF is inner-mass fraction
    //Solver can be looped through steps in PMR
    int findlayer=3; //1 to find core, 2 for mantle, 3 for water, 4 for atmosphere mass fraction
    vector<int> layers={1,1,0,0}; //mark which layers {C,M,W,A} to hold in constant partial ratio
    double minPMR=67.0; //Minimum PMR(%) for iteration, must be multiple of 0.1
    double maxPMR=67.6; //Maximum PMR(%) for iteration, must be multiple of 0.1
    float step=0.2; // Step size for partial mass ratio, must be multiple of 0.1

    double rerr=0.001; // Error in the simulated radius to target radius

    //Multi-threading is commented out by default for easy install
    //Must uncomment Line 2 in Makefile and all occurences of "pragma..." in comfind.cpp 
    int num_threads=0; 
    
    compfinder(Comp,findlayer,layers,minPMR,maxPMR,step,rerr,num_threads,Tgap,ave_rho,P_surface,false,infilename, outfilename);

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
    double fraction = 1;
    int index = 0;
    
    gettimeofday(&start_time,NULL);

    twolayer(index,fraction,P_surface,1,"./run/PosteriorsPPv.txt", "./result/PPvMCMC.txt");
  
    gettimeofday(&end_time, NULL);
  
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;

    cout<<"running time "<<deltat<<'s'<<endl;
  }
  
  else if(input_mode == 7)
  {
    vector<PhaseDgm> Comp = {core, mant, water, atm};
    vector<double> Tgap = {0, 0, 0, 300};
    vector<double> Mcomp =  {0.29,0.69,1.02,0.001};

    fullmodel(Comp,Mcomp,Tgap,ave_rho,P_surface,false,2,"./run/PosteriorsVinet.txt", "./result/IceMCMCnew.txt");
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
  delete Water_ExoPlex;
  delete Water;
  delete Water_sc_dummy;
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
  delete Gold;
  delete Plat;

  return 0;
}
