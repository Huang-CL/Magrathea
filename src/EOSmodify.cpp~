#include "EOSlist.h"
#include "phase.h"
#include "EOSmodify.h"

void twolayer(int index, double fraction, vector<double> &Mp, vector<double> &Rp, double P0, bool printmodel)
// calculate a mass-radius curve for two-layer planet with a constant mass fraction.
// index=0 stands for planet with no iron core.  Only have Si-mantle and ice crust.
// index=1 for planet with no Si-mantle.  2 for no water.
// fraction is the mass fraction of the inner layer.
{
  double Mpt, MC, MM, MW;  

  hydro* planet;
  Mp.resize(0);
  Rp.resize(0);
  
  switch(index)
  {
  case 0:
    for(Mpt=0.1;Mpt<50;Mpt*=1.1)
    {
      MM=Mpt*fraction;
      MW=Mpt-MM;
      MC=0;
      planet=getmass(MC,MM,MW,P0);
      if(printmodel == true)
	planet->print("densitymap.txt");
      Mp.push_back(planet->totalM()/ME);
      Rp.push_back(planet->totalR()/RE);
      delete planet;
    }
    break;
  case 1:
    for(Mpt=0.1;Mpt<50;Mpt*=1.1)
    {
      MC=Mpt*fraction;
      MW=Mpt-MC;
      MM=0;
      planet=getmass(MC,MM,MW,P0);
      if(printmodel == true)
	planet->print("densitymap.txt");
      Mp.push_back(planet->totalM()/ME);
      Rp.push_back(planet->totalR()/RE);
      delete planet;
    }
    break;
  case 2:
    for(Mpt=0.1;Mpt<50;Mpt*=1.1)
    {
      MC=Mpt*fraction;
      MM=Mpt-MC;
      MW=0;
      planet=getmass(MC,MM,MW,P0);
      if(printmodel == true)
	planet->print("densitymap.txt");
      Mp.push_back(planet->totalM()/ME);
      Rp.push_back(planet->totalR()/RE);
      delete planet;
    }
    break;
  default:
    cout<<"Error: Incorrect skip layer index "<<index<<endl;
  };
}

void twolayer(int absent_index, double fraction, double P0, int adjust_index, string eosfile, string outfile)
/*
  Format of the input file

  line 1. number of phases, a list of name of each phases. If the first phase in the list already exist in the database, the name should be the same as the phasetype string (not the name of EOS pointer) exclude parentheses part, e.g. Fe hcp, Si PPv, Water etc.
  line 2. for each phases, how many parameter provided, the index of each parameter according to the index table.  The number of parameters has to be the same for all phases. Then a list of index with their default values.
  line 3 and follows. main body.  beginning pressure, parameters for phase 1, and so on.
 */
{
  PhaseDgm *adjust_cp;
  if (adjust_index == 0)
    adjust_cp = &Fe;
  else if (adjust_index == 1)
    adjust_cp = &Si;
  else if (adjust_index == 2)
    adjust_cp = &water;
  else
  {
    cout<<"Error: index of adjust component "<<adjust_index<<" is not between 0-2."<<endl;
    exit(1);
  }

  if (absent_index == adjust_index)
  {
    cout<<"Error: the phase, "<<adjust_cp->getname()<<", whose EOS impact are investigated is absent in the two layer model."<<endl;
    exit(1);
  }

  ifstream fin(eosfile.c_str());
  if(!fin)
  {
    cout<<"Error: Failed to open the eos input file "<<eosfile<<endl;
    exit(1);
  }

  ofstream fout(outfile.c_str());
  if(!fout)
  {
    cout<<"Error: Failed to open output file "<<outfile<<endl;
    return;
  }
  vector<double> Mp, Rp;
  int nPhases, length, param_index, i, j;
  bool first_iter = true;
  fin >> nPhases;
  EOS **highEOSs = new EOS*[nPhases];
  for(i=0; i<nPhases; i++)
    highEOSs[i] = new EOS();
  string phasename;
  getline(fin,phasename,'\t');
  for(i=0; i<nPhases; i++)
  {
    if(i==nPhases-1)
      getline(fin,phasename);
    else
      getline(fin,phasename,'\t');
    highEOSs[i] -> setphasename(phasename);
  }
  fin >> length;
  string indices;
  stringstream ss;
  double (*params)[2] = new double[length][2];
  getline(fin,indices);
  ss.str(indices);
  for(j=0; j<length; j++)
  {
    ss >> param_index;
    params[j][0] = (double) param_index;
  }
  int k;
  double default_k;
  while(ss>>k)
  {
    ss>>default_k;
    for(i=0; i<nPhases;i++)
      highEOSs[i] -> modifyEOS(k, default_k);
  }

  double *start_P=NULL;
  if(nPhases > 1)
    start_P = new double[nPhases-1];

  while(fin>>params[0][1])
  {
    for(j=1; j<length; j++)
      fin >> params[j][1];
    highEOSs[0] -> modifyEOS(params, length);

    if(nPhases > 1)
    {
      for(i=1; i<nPhases; i++)
      {
	fin >> start_P[i-1];
	for(j=0; j<length; j++)
	  fin >> params[j][1];
	highEOSs[i] -> modifyEOS(params, length);
      }
      adjust_cp->set_phase_highP(nPhases, start_P, highEOSs);
    }

    else
      adjust_cp->set_phase_highP(1, NULL, highEOSs);

    getline(fin,indices);
    twolayer(absent_index, fraction, Mp, Rp, P0);

    if(first_iter)
    {
      first_iter = false;
      for(int j=0; j < int(Mp.size()); j++)
	fout<<Mp[j]<<" \t";
      fout<<endl;
    }
    for(int j=0; j < int(Rp.size()); j++)
      fout<<Rp[j]<<" \t";
    fout<<endl;
  }
  fin.close();

  delete[] highEOSs;
  delete[] params;
  delete[] start_P;
}


void fullmodel(vector<PhaseDgm> &Comp, vector<double> M_Comp, vector<double> Tgap, vector<double> ave_rho, double P0, bool isothermal, int adjust_index, string eosfile, string outfile)
/*
  Format of the input file

  line 1. number of phases, a list of name of each phases. If the first phase in the list already exist in the database, the name should be the same as the phasetype string (not the name of EOS pointer) exclude parentheses part, e.g. Fe hcp, Si PPv, Water etc.
  line 2. for each phases, how many parameter provided, the index of each parameter according to the index table.  The number of parameters has to be the same for all phases. Then a list of index with their default values.
  line 3 and follows. main body.  beginning pressure, parameters for phase 1, and so on.
 */
{
  PhaseDgm *adjust_cp;
  if (adjust_index > 2)
  {
    cout<<"Error: index of adjust component "<<adjust_index<<" is not between 0-2."<<endl;
    adjust_cp=NULL;
    exit(1);
  }
  else
  {
    adjust_cp = &Comp[adjust_index];
  }
  
  ifstream fin(eosfile.c_str());
  if(!fin)
  {
    cout<<"Error: Failed to open the eos input file "<<eosfile<<endl;
    exit(1);
  }

  ofstream fout(outfile.c_str());
  if(!fout)
  {
    cout<<"Error: Failed to open output file "<<outfile<<endl;
    return;
  }

  int nPhases, length, param_index, i, j;
  hydro *planet;
  fin >> nPhases;
  EOS **highEOSs = new EOS*[nPhases];
  for(i=0; i<nPhases; i++)
    highEOSs[i] = new EOS();
  string phasename;
  getline(fin,phasename,'\t');
  for(i=0; i<nPhases; i++)
  {
    if(i==nPhases-1)
      getline(fin,phasename);
    else
      getline(fin,phasename,'\t');
    highEOSs[i] -> setphasename(phasename);
  }
  fin >> length;
  string indices;
  stringstream ss;
  double (*params)[2] = new double[length][2];
  getline(fin,indices);
  ss.str(indices);
  for(j=0; j<length; j++)
  {
    ss >> param_index;
    params[j][0] = (double) param_index;
  }
  int k;
  double default_k;
  while(ss>>k)
  {
    ss>>default_k;
    for(i=0; i<nPhases;i++)
      highEOSs[i] -> modifyEOS(k, default_k);
  }

  double *start_P=NULL;
  if(nPhases > 1)
    start_P = new double[nPhases-1];

  while(fin>>params[0][1])
  {
    for(j=1; j<length; j++)
      fin >> params[j][1];
    highEOSs[0] -> modifyEOS(params, length);

    if(nPhases > 1)
    {
      for(i=1; i<nPhases; i++)
      {
	fin >> start_P[i-1];
	for(j=0; j<length; j++)
	  fin >> params[j][1];
	highEOSs[i] -> modifyEOS(params, length);
      }
      adjust_cp->set_phase_highP(nPhases, start_P, highEOSs);
    }

    else
      adjust_cp->set_phase_highP(1, NULL, highEOSs);

    getline(fin,indices);
    planet = fitting_method(Comp, M_Comp, Tgap, ave_rho, P0, isothermal);

    fout<<planet->totalR()<<endl;
    delete planet;
  }
  fin.close();

  delete[] highEOSs;
  delete[] params;
  delete[] start_P;
}
