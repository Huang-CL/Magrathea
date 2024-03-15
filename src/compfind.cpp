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

  ifstream fin(infile.c_str());
  if(!fin)
  {
    cout<<"Error: Failed to open the eos input file "<<infile<<endl;
    exit(1);
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

  ofstream fout(outfile.c_str(), std::ofstream::trunc);
  if(!fout)
  {
    cout<<"ERROR: failed to open output file., "<<outfile<<"  Exit"<<endl;
    exit(1);
  }
  gettimeofday(&start_time,NULL);
  fout<<"MCore\t MMantle\t MWater\t MGas\t RCore\t RMantle\t RWater\t RPlanet"<<endl;
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

    if(solver == '2')
	    planet = getmass(MC, MM, MW, P0);
    else
	    planet = fitting_method(Comp, {MC, MM, MW, MG}, Tgap, ave_rho, P0, isothermal);
    
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
    int interations=0;   
    for(int j=static_cast<int>(minPMR*10); j<static_cast<int>(maxPMR*10+1); j+=static_cast<int>(step*10)){   // Loop for each core:mantle mass fraction
      double R1=0, R2=0;
      vector<double> Mcomp;
      double fouter0=j/10.0/100.0; 
      double fouter=fouter0;
      double finner=1-fouter0; // Starts 100% Core
      double funk1=0.0, funk0=0.0, funk2=0.0;  // Starts 0% Water        
      do   // Finds water fraction to match posterior sampled radius
      {
        interations=interations+1;
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
          cout<<Mcomp[0]<<' '<<Mcomp[1]<<' '<<Mcomp[2]<<' '<<Mcomp[3]<<endl; //Uncomment to watch solver find mass
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
          cout<<Mcomp[0]<<' '<<Mcomp[1]<<' '<<Mcomp[2]<<' '<<Mcomp[3]<<endl; //Uncomment to watch solver find mass
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
      while(abs(Rtarg[i]-R1)/Rtarg[i]>rerr && interations<30);  // Find Radius with error, break after 30 tries
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

