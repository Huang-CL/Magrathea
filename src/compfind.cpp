//#include "EOSlist.h"
//#include "phase.h"
#include "compfind.h"


void multiplanet(vector<PhaseDgm> &Comp, vector<double> Tgap, int solver, vector<double> ave_rho, double P0, bool isothermal, string infile, string outfile)
// Find planet radii and radii of layers for an input file of planets
// Input file must specify planet mass and mass fractions of each layer
{
  struct timeval start_time, end_time;
  double deltat;
  hydro* planet;

  ifstream fin(infile.c_str());  //open input file
  if (!fin) {
    cout << "Error: cannot open " << infile << endl;
    exit(1);
  }

  vector<double> Mp, fC, fM, fW, Tsurf; //containers
  vector<double> Rs;
  string line;

  getline(fin, line); //skip header

  bool hasTsurf = false;
  size_t lineno  = 1;                // start after header
  while (getline(fin, line)) {
    ++lineno;
    if (line.find_first_not_of(" \t\r\n") == string::npos) continue; // blank

    stringstream ss(line);
    double mp, fc, fm, fw, ts;

    if (!(ss >> mp >> fc >> fm >> fw)) {           // need 4 numbers
      cout << "ERROR: malformed row " << lineno << endl;
      exit(1);
    }

    if (ss >> ts) {               // 5th number,  then we have Tsurf
      hasTsurf = true;
      Tsurf.push_back(ts);
    } else if (hasTsurf) {        // earlier rows had 5, this one only 4
      cout << "ERROR: inconsistent column count at row "
          << lineno << endl;
      exit(1);
    }

    Mp.push_back(mp);
    fC.push_back(fc);
    fM.push_back(fm);
    fW.push_back(fw);
  }

  fin.close();
  const int nline = static_cast<int>(Mp.size());

  cout << "Read " << nline << " planets (" << (hasTsurf?"with":"without")
     << " Tsurf column)\n";

  ofstream fout(outfile.c_str(), ofstream::trunc); //open output
  if (!fout) {
    cout << "ERROR: cannot open " << outfile << endl;
    exit(1);
  }

  gettimeofday(&start_time,NULL);
  fout<<"MCore\t MMantle\t MWater\t MGas\t RCore\t RMantle\t RWater\t RPlanet"<<endl;
  cout<<"Percentage completed:"<<endl;

  //#pragma omp parallel for schedule(dynamic) num_threads(3) private(planet, Rs)   
  for(int i=0; i<nline; i++)
  {
    double  MC, MM, MW, MG;
    MC = fC[i]*Mp[i];
    MM = fM[i]*Mp[i];
    MW = fW[i]*Mp[i];
    MG = Mp[i] - (MW+MC+MM);
    if(hasTsurf==true)
      Tgap[3]=Tsurf[i];

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

    if(solver == '2')
	    planet = getmass(MC, MM, MW, P0);
    else
	    planet = fitting_method(Comp, {MC, MM, MW, MG}, Tgap, ave_rho, P0, isothermal);
    
    //#pragma omp critical
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

}



void compfinder(vector<PhaseDgm> &Comp, int findlayer, vector<int> layers, double minPMR, double maxPMR, float step, double rerr, int num_threads, vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, string infile, string outfile)
// read in a file with a table of mass and radius posterior samples
// Example input file
// Mass Radius
// 1.3617007672434112 1.1000571279817417
// 1.3122419209981564 1.1101895856223976 
// calculate a mass-radius curve for two-layer planet with a constant mass fraction.
{

  struct timeval start_time, end_time;
  double deltat;
  ifstream fin(infile.c_str());
  if(!fin)
  {
    cout<<"Error: Failed to open the eos input file "<<infile<<endl;
    exit(1);
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

  ofstream fout(outfile.c_str(), std::ofstream::trunc);
  if(!fout)
  {
    cout<<"ERROR: failed to open output file., "<<outfile<<"  Exit"<<endl;
    exit(1);
  }
  fout<<"MPlanet\t MCore\t MMantle\t MWater\t MGas\t RCore\t RMantle\t RWater\t RPlanet\t RPosterior"<<endl;

  gettimeofday(&start_time,NULL);
  hydro *planet;
  //#pragma omp parallel for schedule(dynamic) num_threads(3) private(planet, Rs)   
  for(int i=0; i<nline; i++) // Loop for each posterior
  {
    char dontrecord='n';     
    for(int j=static_cast<int>(minPMR*10); j<static_cast<int>(maxPMR*10+1); j+=static_cast<int>(step*10)){   // Loop for each core:mantle mass fraction
      int iterations=0;
      double R1=0, R2=0;
      vector<double> Mcomp;
      double fouter0=j/10.0/100.0; 
      double fouter=fouter0;
      double finner=1-fouter0; // Starts 100% Core
      double funk1=0.0, funk0=0.0, funk2=0.0;  // Starts 0% Water        
      do   // Finds water fraction to match posterior sampled radius
      {
        iterations=iterations+1;
        if (funk1==0) // First iteration
        {
          int count_dumb=0;
          for(int k=0;k<int(layers.size());k++)
          {
            if(layers[k]==0)
              Mcomp.push_back(0.0);
            else
            {
              if(count_dumb==0)
              {
                Mcomp.push_back(finner*Mp[i]);
                count_dumb=count_dumb+1;
              }
              else
                Mcomp.push_back(fouter*Mp[i]);
            }
          }
          Mcomp[findlayer-1]=funk0*Mp[i];
          //cout<<Mcomp[0]<<' '<<Mcomp[1]<<' '<<Mcomp[2]<<' '<<Mcomp[3]<<endl; //Uncomment to watch solver find mass
	        planet = fitting_method(Comp, Mcomp, Tgap, ave_rho, P0, isothermal);
          if (!planet)
    	    {
            funk0=0.0001;  // No solution. Try increasing water fraction
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
    	          dontrecord ='y';
	              j=maxPMR*10+1; // Move to next posterior
	              break;
	            }
	            else // Planet must have a larger core:mantle ratio
	            {
	              dontrecord='y';  // Don't keep this data point
	              j=maxPMR*10+1; // Move to the next posterior
	              break;  
	            }
	          }	                   
	        }	          
	        funk0=funk1;
	        funk1=funk1+0.00001;  // Initial step size 0.00001       
        }
          
        else  // Second iteration onward
        {    
          finner=(1-fouter0)-(1-fouter0)*funk1; // Keep core:mantle ratio constant subtract water fraction proportionally from core and mantle
          fouter=fouter0-fouter0*funk1;
          
          int count_dumb=0;
          for(int k=0;k<int(layers.size());k++)
          {
            if(layers[k]==0)
              Mcomp[k]=0.0;
            else
            {
              if(count_dumb==0)
              {
                Mcomp[0]=finner*Mp[i];
                count_dumb=count_dumb+1;
              }
              else
                Mcomp[k]=fouter*Mp[i];
            }
          }
          Mcomp[findlayer-1]=funk1*Mp[i];
          //cout<<Mcomp[0]<<' '<<Mcomp[1]<<' '<<Mcomp[2]<<' '<<Mcomp[3]<<endl; //Uncomment to watch solver find mass
	        planet = fitting_method(Comp, Mcomp, Tgap, ave_rho, P0, isothermal);
          if (!planet)
    	    {
    	      funk1=funk1+0.0001; // No solution. Try increasing water fraction
    	      R1=0;
    	    }  
          else
          {
	          Rs = planet -> getRs();
	          R2 = Rs[Rs.size()-1];
	          funk2=funk1-(funk1-funk0)*(R2-Rtarg[i])/(R2-R1);  // Secant Method
	          if (funk2<=0) // Secant method returned negative result
    	      {
              dontrecord='y';
              j=maxPMR*10+1; // Move to next posterior
              break;
	          }
	          funk0=funk1;
	          funk1=funk2;
	          R1=R2;
	          delete planet;
          }
        }       
      }
      while(abs(Rtarg[i]-R1)/Rtarg[i]>rerr && iterations<30);  // Find Radius with error, break after 30 tries
      //#pragma omp critical
      if (dontrecord!='y') // Don't record if the posterior radius is less than radius with 0% water
      {
        fout<<Mp[i]<<"\t "<<Mcomp[0]<<"\t "<<Mcomp[1]<<"\t "<<Mcomp[2]<<"\t "<<Mcomp[3]<<"\t ";
	      for(int k=0; k < int(Rs.size()); k++)
	        fout<<Rs[k]<<"\t ";
        fout<<Rtarg[i]<<endl;
      }  
      else
        dontrecord='n';
    } // End of for loop for each core:mantle ratio
    if (100*(i+1) / nline > 100*i/nline) // progress bar
	    cout<<100*(i+1)/nline<<'%'<<'\r'<<std::flush;
  }  // End of loop for each posterior sample
  cout<<endl;
  fout.close();
  
  gettimeofday(&end_time, NULL);
  
  deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;

  cout<<"running time "<<deltat<<'s'<<endl;


  delete[] Mp;
  delete[] Rtarg;
}


void mcmcsample(vector<PhaseDgm> &Comp, double MassPrior, double MUncPrior, double RadPrior, double RUncPrior,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, string outfile, int numchains, int steps, int numlayers)
//Overall loop for each chain, may be parallized by uncommenting "#pragma omp" lines
{
  ofstream fout(outfile);
  fout << "Mass\t fCore\t fMantle\t fWater\t fAtm\t log_likelihood\t RCore\t RMantle\t RWater\t RPlanet"<<endl;
    //#pragma omp parallel for schedule(dynamic) num_threads(3)
    for (int i = 0; i < numchains; ++i)
    {
        vector<MCMCRecord> chain;
        metropolis_hastings(chain, Comp,
                            MassPrior, MUncPrior,
                            RadPrior,  RUncPrior,
                            Tgap, ave_rho,
                            P0, isothermal,
                            steps, numlayers);

        //-------------------------- summary print ---------------------------
        std::ostringstream summary;
        MCMCRecord last = chain.back();
        summary << "Chain " << i << " final sample after " << steps << " steps:\n"
                << "  Mass    = " << last.Mass    << "\n"
                << "  fCore   = " << last.fCore   << "\n"
                << "  fMantle = " << last.fMantle << "\n"
                << "  fWater  = " << last.fWater  << "\n"
                << "  fAtm    = " << last.fAtm    << "\n";

        //---------------------- build thread-local output -------------------
        std::ostringstream local;
        for (const auto& rec : chain)
        {
            local << rec.Mass        << '\t'
                  << rec.fCore       << '\t'
                  << rec.fMantle     << '\t'
                  << rec.fWater      << '\t'
                  << rec.fAtm        << '\t'
                  << rec.log_likelihood << '\t'
                  << rec.RCore       << '\t'
                  << rec.RMantle     << '\t'
                  << rec.RWater      << '\t'
                  << rec.RPlanet     << '\n';
        }
        //single critical section for ompenmp
        //#pragma omp critical
        {
            cout << summary.str();        // console output
            fout << local.str();               // serialised file write
        }
    } // end parallel for
    fout.close();
    cout << "MCMC chains saved to " << outfile << endl;
  } 

  void metropolis_hastings(vector<MCMCRecord>& chain,vector<PhaseDgm> &Comp, double MassPrior, double MUncPrior, double RadPrior, double RUncPrior,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, int steps, int numlayers)
  //Primary MCMC sampling and walker
  {
    /* ---------- common RNG ---------- */
    random_device rd;
    default_random_engine gen(rd());
    normal_distribution<double> prop_Mass(MassPrior, MUncPrior);
    normal_distribution<double> prop_dist(0.0, 0.1);
    uniform_real_distribution<double> uniform(0.0, 1.0);
    uniform_real_distribution<double> log_uniform(log(1e-5), log(0.98)); //for atmosphere choose log uniform
    normal_distribution<double> prop_logatm(0.0, 0.5); //Atmosphere step proposed

    /* ---------- helpers ---------- */
    // 1. Generate an initial Params object that obeys the layer‑count rules
    auto init_params = [&](int L) -> Params {
        if (L == 2) {
            double c = uniform(gen);
            return {MassPrior, c, 1.0 - c, 0.0};
        }
        if (L == 3) {
            double c = uniform(gen);
            double m = uniform(gen) * (1.0 - c);
            return {MassPrior, c, m, 1.0 - c - m};
        }
        if (L == 4) {
            double fAtm = std::exp(log_uniform(gen));
            double c    = uniform(gen) * (1.0 - fAtm);
            double m    = uniform(gen) * (1.0 - c - fAtm);
            return {MassPrior, c, m, 1.0 - c - m - fAtm};
        }
        throw std::invalid_argument("numlayers must be 2, 3 or 4");
    };

    // 2. Propose a new Params and return it together with the atmosphere fraction (needed only for L==4)
    auto propose = [&](const Params& cur, int L) -> std::pair<Params,double> {
        Params p = cur;
        p.Mass = std::max(0.01, prop_Mass(gen));   // keep mass > 0.01

        if (L == 2) {
            p.fCore = cur.fCore + prop_dist(gen);
            if (p.fCore < 0.0) p.fCore = 0.0;
            else if (p.fCore > 1.0) p.fCore = 1.0; //modify for prior constraints
            p.fMantle = 1.0 - p.fCore;
            p.fWater  = 0.0;
            return {p, 0.0};
        }

        if (L == 3) {
            do {
                p.fCore   = cur.fCore   + prop_dist(gen);
                p.fMantle = cur.fMantle + prop_dist(gen);
            } while (p.fCore < 0.0 || p.fMantle < 0.0 ||
                     p.fCore + p.fMantle > 1.0); //modify for prior constraints
            p.fWater = 1.0 - p.fCore - p.fMantle;
            return {p, 0.0};
        }

        /* L == 4 */
        double fAtm;
        do {
            p.fCore   = cur.fCore   + prop_dist(gen);
            p.fMantle = cur.fMantle + prop_dist(gen);
            fAtm      = std::exp(std::log(1.0 - cur.fCore - cur.fMantle - cur.fWater) +
                                 prop_logatm(gen));
            p.fWater  = 1.0 - p.fCore - p.fMantle - fAtm;
        } while (p.fCore < 0.0 || p.fMantle < 0.0 || p.fWater < 0.0); //modify for prior constraints
        return {p, fAtm};
    };

    /* ---------- initial state ---------- */
    Params current = init_params(numlayers);
    LikelihoodResult lr = log_likelihood(Comp, MassPrior, MUncPrior, RadPrior, RUncPrior,
                                         Tgap, ave_rho, P0, isothermal,
                                         current.Mass, current.fCore, current.fMantle, current.fWater);

    double logL_cur = lr.log_likelihood;
    double fAtm_cur = (numlayers == 4) ? 1.0 - current.fCore - current.fMantle - current.fWater : 0.0;

    chain.clear();
    chain.reserve(steps + 1);
    chain.push_back({current.Mass, current.fCore, current.fMantle, current.fWater,
                     fAtm_cur, logL_cur, lr.RCore, lr.RMantle, lr.RWater, lr.RPlanet});

    /* ---------- main MH loop ---------- */
    for (int s = 0; s < steps; ++s) {
        std::pair<Params,double> prop_pair = propose(current, numlayers);
        Params  prop      = prop_pair.first;
        //double  fAtm_prop = prop_pair.second; // (unused except for debugging)

        LikelihoodResult lr_prop = log_likelihood(Comp, MassPrior, MUncPrior, RadPrior, RUncPrior,
                                                  Tgap, ave_rho, P0, isothermal,
                                                  prop.Mass, prop.fCore, prop.fMantle, prop.fWater);

        double log_alpha = lr_prop.log_likelihood - logL_cur;

        if (std::log(uniform(gen)) < log_alpha) {   // accept
            current  = prop;
            lr       = lr_prop;
            logL_cur = lr_prop.log_likelihood;
            fAtm_cur = (numlayers == 4)
                         ? 1.0 - current.fCore - current.fMantle - current.fWater
                         : 0.0;
        }

        chain.push_back({current.Mass, current.fCore, current.fMantle, current.fWater,
                         fAtm_cur, logL_cur, lr.RCore, lr.RMantle, lr.RWater, lr.RPlanet});
    }
}

LikelihoodResult log_likelihood(vector<PhaseDgm> &Comp, double MassPrior, double MUncPrior, double RadPrior, double RUncPrior,  vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, double Mass, double fCore, double fMantle, double fWater) {
    // Find planet radius and calculate log likelihood
    hydro *planet;
    double  MC, MM, MW, MG;
    vector<double> Rs;
    MC = fCore*Mass;
    MM = fMantle*Mass;
    MW = fWater*Mass;
    MG = Mass - (MW+MC+MM);
    if(MG<1E-14)
      MG=0;
    planet = fitting_method(Comp, {MC,MM,MW,MG}, Tgap, ave_rho, P0, isothermal); //Full Solver
    if (!planet)
    {
      return {-1e9, 0, 0, 0, 0};
    }
    Rs = planet -> getRs();
    delete planet;

    // ln(L) = -0.5 [ ((M_model - M_obs)/sigM)^2 + ((R_model - R_obs)/sigR)^2 ]
    double chi2_mass = pow((Mass - MassPrior) / MUncPrior, 2);
    double chi2_rad  = pow((Rs.back() - RadPrior) / RUncPrior, 2);
    double lnL = -0.5 * (chi2_mass + chi2_rad);

    return {lnL, Rs[0], Rs[1], Rs[2], Rs.back()};
}
