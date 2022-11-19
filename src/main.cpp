#include "define.h"
#include "EOS.h"
#include "EOSmodify.h"
#include "phase.h"
#include "hydro.h"

struct timeval start_time, end_time;

const double rho_eps_rel = 1E-11;	// relative error tolerance of the density solver
const double T_eps_rel = 1E-11;	// relative error tolerance of the temperature solver for the adiabatic profile
const double ode_eps_rel0 = 1E-7; // relative error tolerance for first round ode integrator (1E-7) used in mode 0
const double ode_eps_rel1 = 1E-10; // relative error tolerance for second round ode integrator (1E-10) used in mode 0
int fit_iter;
const double R_eps_rel = 2E-5; // relative error tolerance in mode 0 first round radius determination (5E-5).  Should be around sqrt(ode_eps_rel0).
const double ode_eps_rel2 = 1E-10; // relative error tolerance for ode integrator (1E-10) used in mode 1
const double P_eps_rel = 1E-10;	// relative error tolerance in mode 1 central pressure determination (1E-10).  Should not be more restrict than ode_eps_rel.
const double fit_eps_rel = 1E-4; // the relative error tolerance at the fitting point in mode 0 round 2 (1E-4). Should be four orders of magnitudes larger than ode_eps_rel1.
vector<double> ave_rho = {15, 5, 2, 1E-3};// Assuming the density of the core is 15, mantle is 5, water is 2, and gas is 1E-3.
const bool verbose = false;		  // Whether print warnings.
const double P_surface = 1E5;		  // The pressure level that the broad band optical transit radius probes. (in microbar)
int count_shoot = 0;			  // used to count the total number of shootings per each solution
int count_step = 0;			  // used to count the sum of integral steps in all shooting results.

int main()
{
  gsl_set_error_handler_off();  //Dismiss the abortion from execution, which is designed for code testing.
  hydro* planet;
  ifstream fin;
  int input_mode=0;
  // Choose between the 8 available input_mode values:
  // 0: regular solver, 1: temperature-free solver, 2: two-layer solver, 
  // 3: modify a built-in EOS on they fly, 
  // 4: iterate over EOS modifications with two-layer solver, 5: iterate over EOS with regular solver
  // 6: bulk input mode with regular solver
  // 7: composition finder *in devopment*, secant method to find third layer mass to match a mass and radius measurement

  if (input_mode == 0)
  {
    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap = {0, 0, 0, 300};
    // The temperature of the outer boundary of the inner component minus the inner boundary of the outer component.  A positive number indicates temperature increases inward.  0 indicates the temperature is continuous at the boundary of components.  The last number is the planetary surface temperature.
    vector<double> Mcomp =  {1.0,0.5,0.1,0.00001}; // Mass in Earth Masses of Core, Mantle, Hydrosphere, Atmosphere
    planet=fitting_method(Comp, Mcomp, Tgap, ave_rho, P_surface, false);
    cout<<count_shoot<<' '<<count_step<<endl;
    if (!planet)
    {
      for (unsigned int i=0; i < Mcomp.size(); i++)
	cout<<Mcomp[i]<<", ";
      cout<<"\t No solution found."<<endl;
    }
    else
      planet->print("./result/Structure.txt", true); // Save the result in an asc file with this name.

    delete planet;
  }  

  else if(input_mode == 1)
  {
    planet=getmass(3,3,3,P_surface);
    // Mass in Earth Masses of Core, Mantle, Hydrosphere
    if (planet)
      planet->print("./result/Structure.txt");
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

  else if(input_mode == 4)
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
  
  else if(input_mode == 5)
  {
    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap = {0, 0, 0, 300};
    vector<double> Mcomp =  {0.29,0.69,1.02,0.001};

    fullmodel(Comp,Mcomp,Tgap,ave_rho,P_surface,false,2,"./run/PosteriorsVinet.txt", "./result/IceMCMCnew.txt");
  }

  else if(input_mode == 6)
  {
    double deltat;
    // read in a file with a table of mass fractions in the unit of earth's mass.  Calculate the radius of such planets and write in a file.

    /*
    Example input file
    Mass, fCore, fMantle, fWater
    2     0.2    0.4      0.4
    1.5   0.5    0.4      0.1
    */
    // Tgap and File names are prompted in command line for user entry

    string filename;
    cout<<"This function requires an input file that contains the planet mass in the unit of the Earth mass and the mass fraction of core, mantle and water.  The first column is the total planet mass in Earth unit.  The sum of the mass fraction should be less than or equal to 1.  The remaining mass would be put in the ideal gas layer.  The columns from the left to the right should be in planet mass, core, mantle, and water mass fraction.  The first row of the table would be considered as the header and WILL NOT be read in."<<endl<<"The please input the filename (with relative path from current directory):"<<endl;
    cin >> filename;
    fin.open(filename.c_str());
    if(!fin)
    {
      cout<<"ERROR: failed to open input planet list file, "<<filename<<" Exit."<<endl;
      return 1;
    }

    double *Mp, *fW, *fC, *fM, MC, MM, MW, MG;
    vector<double> Rs;
    string sline;
    int nline = 0;

    getline(fin,sline);
    streampos beginpos = fin.tellg();

    while(getline(fin,sline))
    {
      if(!sline.empty())
	nline++;
    }

    fin.clear();
    fin.seekg(beginpos);

    Mp = new double[nline];
    fC = new double[nline];
    fM = new double[nline];
    fW = new double[nline];

    for(int i=0; i<nline; i++)
      fin>>Mp[i]>>fC[i]>>fM[i]>>fW[i];

    fin.close();

    cout<<"Please input file name to store the output result.  Code will create a new file with the name provided.  If the file exists, any contents that existed will be discarded. \n file name:"<<endl;
    cin >> filename;
    ofstream fout(filename.c_str(), std::ofstream::trunc);
    if(!fout)
    {
      cout<<"ERROR: failed to open output file., "<<filename<<"  Exit"<<endl;
      return 1;
    }

    char mode='1';
    cout<<"Calculation mode. Mode 2 is twice faster, but only applicable for planet with no gas layer and whose temperature effect is not important."<<endl;
    cout<<"Choose the mode to use (1 or 2. default 1):";
    cin>>mode;

    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap;
    if (mode!='2')
    {
      double T1, T2, T3, T4;
      cout<<"Please input the temperature jump at the boundary of Fe core and Silicate mantle. Positive number if core is hotter than mantle, 0 if temperature is continuous at the boundary."<<endl;
      cin>>T1;
      Tgap.push_back(T1);
      cout<<"Please input the temperature jump at the boundary of Silicate mantle and ice shell."<<endl;
      cin>>T2;
      Tgap.push_back(T2);
      cout<<"Please input the temperature jump at the boundary of ice shell and gas layer."<<endl;
      cin>>T3;
      Tgap.push_back(T3);
      cout<<"Please input the planet surface temperature."<<endl;
      cin>>T4;
      Tgap.push_back(T4);
    }
    gettimeofday(&start_time,NULL);
  
    fout<<"MCore, MMantle, MWater, MGas, RCore, RMantle, RWater, RPlanet"<<endl;
    cout<<"Percentage completed:"<<endl;
    for(int i=0; i<nline; i++)
    {
      MC = fC[i]*Mp[i];
      MM = fM[i]*Mp[i];
      MW = fW[i]*Mp[i];
      MG = Mp[i] - (MW+MC+MM);

      if(MG < 0)
      {
	if(MG > -0.001*Mp[i])	// mass of gas smaller than 0 probably due to round error, fix the problem by reduce water mass and set gas mass to 0.
	{
	  MW += MG;
	  MG = 0;
	}
	else
	{
	  cout<<"Mass, fCore, fMantle, fWater "<<Mp[i]<<' '<<fC[i]<<' '<<fM[i]<<' '<<fW[i]<<" gives negative gas mass "<<MG<<".  nan will be in the output file.";
	  fout<<MC<<"\t "<<MM<<"\t "<<MW<<"\t "<<MG<<"\t nan"<<endl;
	  continue;
	}
      }
      if(mode == '2')
	planet = getmass(MC, MM, MW, P_surface);
      else
	planet = fitting_method(Comp, {MC, MM, MW, MG}, Tgap, ave_rho, P_surface, false);
      cout<<count_shoot<<' '<<count_step<<endl;
      if (!planet)
	fout<<MC<<"\t "<<MM<<"\t "<<MW<<"\t "<<MG<<"\t No solution found."<<endl;
      else
      {
	Rs = planet -> getRs();
	fout<<MC<<"\t "<<MM<<"\t "<<MW<<"\t "<<MG<<"\t ";
	for(int j=0; j < int(Rs.size()); j++)
	  fout<<Rs[j]<<"\t ";
	if (planet->getstatus() == 1)
	  fout<<" Dummy EOS used.";
	else if (planet->getstatus() == 2)
	  fout<<" Solution did not converge.";
	fout<<endl;
	delete planet;
      }

      if (100*(i+1) / nline > 100*i/nline) // progress bar
	cout<<100*(i+1)/nline<<'%'<<'\r'<<std::flush;
    }  
    cout<<endl;
    fout.close();

    gettimeofday(&end_time, NULL);
  
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;

    cout<<"running time "<<deltat<<'s'<<endl;

    delete[] Mp;
    delete[] fW;
    delete[] fC;
    delete[] fM;
  }

  else if(input_mode == 7)
  {
    // read in a file with a table of mass and radius posterior samples
    
    // Example input file
    // Mass Radius
    // 1.3617007672434112 1.1000571279817417
    // 1.3122419209981564 1.1101895856223976 
    
    string filename="./run/T1hpost.txt";
    fin.open(filename.c_str());
    if(!fin)
    {
      cout<<"ERROR: failed to open input planet list file, "<<filename<<" Exit."<<endl;
      return 1;
    }

    double *Mp, *Rtarg;
    vector<double> Rs;    
    string sline;
    int nline = 0;

    getline(fin,sline);
    streampos beginpos = fin.tellg();

    while(getline(fin,sline))
    {
      if(!sline.empty())
	    nline++;
    }

    fin.clear();
    fin.seekg(beginpos);

    Mp = new double[nline];
    Rtarg = new double[nline];

    for(int i=0; i<nline; i++)
      fin>>Mp[i]>>Rtarg[i];

    fin.close();
    filename="./result/T1houtput.txt";
    ofstream fout(filename.c_str(), std::ofstream::trunc);
    if(!fout)
    {
      cout<<"ERROR: failed to open output file., "<<filename<<"  Exit"<<endl;
      return 1;
    }

    fout<<"MPlanet, MCore, MMantle, MWater, RCore, RMantle, RWater, RPlanet, RPosterior"<<endl;
//    #pragma omp parallel for schedule(dynamic) num_threads(3) private(planet, Rs)   
    for(int i=0; i<nline; i++) // Loop for each posterior
    {
      int step=1; // Step size for core/manlte %
      double rerr=0.001; // Acceptable error in the simulated radius to target radius
      double MC, MM, MW, MG;
      vector<PhaseDgm> Comp = {Fe, Si, water, atm};
      vector<double> Tgap = {0,0,0,300}; // Using full temperature solver, Temperature gap between each layer and surface temperature.
      char final='n';     
      for(int j=0; j<101; j+=step){   // Loop for each core:mantle mass fraction
        double R1=0, R2=0;
        double fM0=j/100.0; 
        double fM=fM0;   // Constant Core:Mantle Ratio, start at 1, mantle % increases by j
        double fC=1-fM0; // Starts 100% Core
        double fW1=0.0, fW0=0.0, fW2=0.0;  // Starts 0% Water        
        do   // Finds water fraction to match posterior sampled radius
        {
          if (fW1==0) // First iteration
          {
            MC = fC*Mp[i];
            MM = fM*Mp[i];
            MW = fW0*Mp[i];
            MG = 0;

	    planet = fitting_method(Comp, {MC, MM, MW, MG}, Tgap, ave_rho, P_surface, false);
            if (!planet)
    	    {
              fW0=0.0001;  // No solution. Try increasing water fraction
              R1=0;
            }
            else
            {
	          Rs = planet -> getRs();
	          delete planet;
	          R1 = Rs[Rs.size()-1]; // Radius of planet
	          if (Rtarg[i]<R1)  
	          {
    	        if (j==0) // Posterior radius is less than 100% core planet
    	        {
    	          fout<<Mp[i]<<"\t "<<Rtarg[i]<<"\t Radius too small for mass."<<endl; 
    	          final='y';
	              j=101; // Move to next posterior
	              break;
	            }
	            else // Planet must have a larger core:mantle ratio
	            {
	              final='y';  // Don't keep this data point
	              j=101; // Move to the next posterior
	              break;  
	            }
	          }	                   
	        }	          
	        fW0=fW1;
	        fW1=fW1+0.001;  // Initial step size 0.001       
          }
          
          else  // Second iteration onward
          {
            
            fC=(1-fM0)-(1-fM0)*fW1; // Keep core:mantle ratio constant subtract water fraction proportionally from core and mantle
            fM=fM0-fM0*fW1;
            MC = fC*Mp[i];
            MM = fM*Mp[i];
            MW = fW1*Mp[i];
            MG = 0;

	    planet = fitting_method(Comp, {MC, MM, MW, MG}, Tgap, ave_rho, P_surface, false);

            if (!planet)
    	    {
    	      fW1=fW1+0.0001; // No solution. Try increasing water fraction
    	      R1=0;
    	    }  
            else
            {
	          Rs = planet -> getRs();
	          R2 = Rs[Rs.size()-1];
	          fW2=fW1-(fW1-fW0)*(R2-Rtarg[i])/(R2-R1);  // Secant Method
	          if (fW2<=0) // Secant method returned negative result
    	          {
    	          final='y';
	              j=101; // Move to next posterior
	              break;
	          }
	          fW0=fW1;
	          fW1=fW2;
	          R1=R2;
	          delete planet;
            }
          }       
        }
        while(abs(Rtarg[i]-R1)/Rtarg[i]>rerr);  // Find Radius with error
//        #pragma omp critical
        if (final!='y') // Don't record if the posterior radius is less than radius with 0% water
        {
          fout<<Mp[i]<<"\t "<<MC<<"\t "<<MM<<"\t "<<MW<<"\t ";
	      for(int k=0; k < int(Rs.size()); k++)
	        fout<<Rs[k]<<"\t ";
          fout<<Rtarg[i]<<endl;
        }  
      } // End of for loop for each core:mantle ratio
    }  // End of loop for each posterior sample

    fout.close();

    delete[] Mp;
    delete[] Rtarg;
  }

  // ============================
  
  delete Fe_liquid;
  delete Fe_liquid2;
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
  delete Ice_Seager;
  delete Water_ExoPlex;
  delete Water;
  delete IceIh_ExoPlex;
  delete IceVI_ExoPlex;
  delete IceVII_ExoPlex;
  delete IceVII;
  delete IceVIIp;
  delete IceVII_FFH2004;
  delete IceVII_FFH2004fit;
  delete IceVII_FFH2004BM;
  delete IceX_HS;
  delete IceX;
  delete IceZeng2013FFH;
  delete IceZeng2013FMNR;
  delete Ice_Dummy;
  delete Gas;

  return 0;
}
