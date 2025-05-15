#include "EOS.h"

EOS::EOS():phasetype(""),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), temptable(NULL), adiabattable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
}

EOS::EOS(string phaseinput, double params[][2], int length):phasetype(phaseinput),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), temptable(NULL), adiabattable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  
  for(int i=0;i<length;i++)
  {
    switch(lround(params[i][0]))
    {
    case 0:
      eqntype=round(params[i][1]);
      break;
    case 1:
      V0=params[i][1];
      break;
    case 2:
      K0=params[i][1];
      break;
    case 3:
      K0p=params[i][1];
      break;
    case 4:
      K0pp=params[i][1];
      break;
    case 5:
      mmol=params[i][1];
      break;
    case 6:
      P0=params[i][1];
      break;
    case 7:
      Theta0=params[i][1];
      break;
    case 8:
      gamma0=params[i][1];
      break;
    case 9:
      beta=params[i][1];
      break;
    case 10:
      gammainf=params[i][1];
      break;
    case 11:
      gamma0p=params[i][1];
      break;
    case 12:
      e0=params[i][1];
      break;
    case 13:
      g=params[i][1];
      break;
    case 14:
      n=round(params[i][1]);
      break;
    case 15:
      Z=round(params[i][1]);
      break;
    case 16:
      T0=params[i][1];
      break;
    case 17:
      alpha0=params[i][1];
      break;
    case 18:
      alpha1=params[i][1];
      break;
    case 19:
      xi=params[i][1];
      break;
    case 20:
      cp_a=params[i][1];
      break;
    case 21:
      cp_b=params[i][1];
      break;
    case 22:
      cp_c=params[i][1];
      break;
    case 23:
      Debye_approx = params[i][1]>0 ? true : false;
      break;
    case 24:
      thermal_type = params[i][1];
      break;
    case 25:
      at1 = params[i][1];
      break;
    case 26:
      at2 = params[i][1];
      break;
    case 27:
      at3 = params[i][1];
      break;
    case 28:
      at4 = params[i][1];
      break;
    case 29:
      ap1 = params[i][1];
      break;
    case 30:
      ap2 = params[i][1];
      break;
    case 31:
      ap3 = params[i][1];
      break;
    case 32:
      ap4 = params[i][1];
      break;
    default:
      cout<<"Error: Incorrect index "<<round(params[i][0])<<" for EOS constructor "<<phasetype<<endl;
      exit(1);
    };
  }
  if (eqntype == 6)
  {
    if (!gsl_finite(mmol))	// default mean molecular weight of gas.  Mix of hydrogen and helium
      mmol = 2.3;
    if (n<0)		// number of atom
      n = 2;
  }

  if (eqntype == 6)		// ideal gas
    thermal_type = 3;
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(gamma0)) // thermal type not specified and have enough information to calculate using Grueneisen parameter
  {
    if (gsl_finite(Theta0))
    {
      if (gsl_finite(e0) && gsl_finite(g)) // Also have enough information to get Pe
	thermal_type = 7;
      else
	thermal_type = 6;
    }
    else
      thermal_type = 5;
  }
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(cp(300)) && gsl_finite(alpha(10,300)) && gsl_finite(T0)) // thermal type not specified and have enough information to calculate using thermal expansion
    thermal_type = 9;
  
  if (eqntype >= 8)		// RTpress EOS style
  {
    cout<<"Error: Incorrect equation type. Type "<<eqntype<<" is an index for RTpress style EOS. It cannot be initialized by this constructor."<<endl;
    exit(1);
  }
}

EOS::EOS(string phaseinput, string filename):phasetype(phaseinput),eqntype(7), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), temptable(NULL), adiabattable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  ifstream fin;
  string sline;
  fin.open(filename.c_str());
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  
  if(!fin)
  {
    if (verbose)
      cout<<"Warning: Failed to open input EOS file "<<filename<<" \n This MAY cause segmentation fault in the run time if this phase "<<phaseinput<<" is used in the runtime!!"<<endl;
    return;
  }

  nline = 0;
  getline(fin,sline);
  streampos beginpos=fin.tellg();
  stringstream stemp;
  bool previoustab=false;

  for(size_t i=0; i<sline.size()-1; i++){ //Checks columns for 2D or 3D table
    if(sline[i] == '\t')
      {
      if(!previoustab){		// two continuous tabs doesn't add new table column
        tabletype++;}    
        previoustab = true;
      } 
      else {
        previoustab = false;
      }   
    }  

  while(getline(fin,sline))
  {
   if(!sline.empty())
     nline++;
  }
  fin.clear();
  fin.seekg(beginpos);

  if(tabletype==4) //Pressure, Temperature, Density, Adiabatic Gradiant Table
  {
    tlen=0;
    float ptemp=0, ptemp1, ttemp, rhotemp,adiabattemp;
    for(int i=0; i<nline; i++){
      fin>>ptemp1>>ttemp>>rhotemp>>adiabattemp; 
      if(ptemp1!=ptemp && tlen>0)   //Table must be ordered by (i, j) for i in P for j in T
        break;
      else{
        tlen++;
        ptemp=ptemp1;
      }
    }
    fin.clear();
    fin.seekg(beginpos);
    if(nline % tlen != 0)
    {
      if (verbose)
        cout<<"Warning: The input EOS file "<<filename<<" is not rectangular"<<endl;
      return;
    }
      
    rhotable=new double[nline];
    Ptable=new double[nline/tlen];
    temptable=new double[tlen];
    adiabattable=new double[nline];
    thermal_type=2;
   
    fin>>Ptable[0]>>temptable[0]>>rhotable[0]>>adiabattable[0];
    fin.seekg(beginpos);
    for(int i=0; i<nline; i++){
      fin>>ptemp>>ttemp>>rhotable[i]>>adiabattable[i];
      if(i<tlen){
        temptable[i]=ttemp;     
      }
      else if(i%tlen==0){
        Ptable[i/tlen]=ptemp;
      }
    }
    fin.close();

    if(Ptable[1]<Ptable[0] || temptable[1]<temptable[0])
    {
      if (verbose)
        cout<<"Warning: Input EOS file "<<filename<<" is not ordered correctly, please refer to README"<<endl;
      return;
    }

    accP=gsl_interp_accel_alloc ();
    accT=gsl_interp_accel_alloc ();
    spline2drho = gsl_spline2d_alloc(gsl_interp2d_bilinear, tlen, nline/tlen);
    spline2dadi = gsl_spline2d_alloc(gsl_interp2d_bilinear, tlen, nline/tlen);
    gsl_spline2d_init(spline2drho,temptable,Ptable,rhotable,tlen,nline/tlen);
    gsl_spline2d_init(spline2dadi,temptable,Ptable,adiabattable,tlen,nline/tlen);

  }  
  
  else if(tabletype==2) //Pressure Density Table
  {
    rhotable=new double[nline];
    Ptable=new double[nline];

    // grl interpolation require strictly increasing x.  Reverse the array order if x is decreasing.
    fin>>rhotable[0]>>Ptable[0];
    fin>>rhotable[1]>>Ptable[1];
    fin.seekg(beginpos);

    if(Ptable[0]>Ptable[1])    
     for(int i=nline-1; i>=0; i--)
       fin>>rhotable[i]>>Ptable[i];
    else
      for(int i=0; i<nline; i++)
        fin>>rhotable[i]>>Ptable[i];
  
    fin.close();
  
    accP = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_steffen, nline);
    gsl_spline_init (spline, Ptable, rhotable, nline);
  }

  else
  {
    if (verbose)
      cout<<"Warning: Requires 2 column (P, rho) or 4 column (P,T,rho,dT/dP_S) table in input EOS file "<<filename<<" \n This MAY cause segmentation fault in the run time if this phase "<<phaseinput<<" is used in the runtime!!"<<endl;
    return;
  }

}


EOS::EOS(string phaseinput, double (*f)(double P, double T), double (*g)(double rho, double T)):phasetype(phaseinput),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), rhotable(NULL), Ptable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  density_extern=f;
  entropy_extern=g;
  dTdP = NULL;
  temptable = NULL;
  adiabattable = NULL;

  if (entropy_extern)
    thermal_type = 1;
  else
    thermal_type = 0;
}

EOS::EOS(string phaseinput, double *Plist, double *rholist, int len_list):phasetype(phaseinput),eqntype(7), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), bn(0), accT(NULL), spline2drho(NULL), spline2dadi(NULL), nline(len_list), tlen(0)
{
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  temptable = NULL;
  adiabattable = NULL;
  rhotable=new double[nline];
  Ptable=new double[nline];

  // grl interpolation require strictly increasing x.  Reverse the array order if x is decreasing.
  
  if(Plist[0]>Plist[1])    
    for(int i=0; i<nline; i++)
    {
      rhotable[i] = rholist[nline-1-i];
      Ptable[i]   = Plist[nline-1-i];
    }
  else
    for(int i=0; i<nline; i++)
    {
      rhotable[i] = rholist[i];
      Ptable[i]   = Plist[i];
    }
  
  accP = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_steffen, nline);
  gsl_spline_init (spline, Ptable, rhotable, nline);
}

EOS::EOS(string phaseinput, double params[][2], double bparams[], int length, int blength):phasetype(phaseinput),eqntype(8), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(8), rhotable(NULL), Ptable(NULL), bn(blength), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  // construction EOS for RTpress
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  temptable = NULL;
  adiabattable = NULL;

  for(int i=0;i<length;i++)
  {
    switch(lround(params[i][0]))
    {
    case 0:
      eqntype=round(params[i][1]);
      break;
// eqntype has to be equal or larger than 8
    case 1:
      V0=params[i][1];
      break;
    case 2:
      K0=params[i][1];
      break;
    case 3:
      K0p=params[i][1];
      break;
    case 4:
      K0pp=params[i][1];
      break;
    case 5:
      mmol=params[i][1];
      break;
    case 6:
      P0=params[i][1];
      break;
    case 7:
      Theta0=params[i][1];
      break;
    case 8:
      gamma0=params[i][1];
      break;
    case 9:
      beta=params[i][1];
      break;
    case 10:
      gammainf=params[i][1];
      break;
    case 11:
      gamma0p=params[i][1];
      break;
    case 12:
      e0=params[i][1];
      break;
    case 13:
      g=params[i][1];
      break;
    case 14:
      n=round(params[i][1]);
      break;
    case 15:
      Z=round(params[i][1]);
      break;
    case 16:
      T0=params[i][1];
      break;
    case 17:
      alpha0=params[i][1];
      break;
    case 18:
      alpha1=params[i][1];
      break;
    case 19:
      xi=params[i][1];
      break;
    case 20:
      cp_a=params[i][1];
      break;
    case 21:
      cp_b=params[i][1];
      break;
    case 22:
      cp_c=params[i][1];
      break;
    case 23:
      Debye_approx = params[i][1]>0 ? true : false;
      break;
    case 24:
      break;		// thermal_type has to be 8
    case 25:
      at1 = params[i][1];
      break;
    case 26:
      at2 = params[i][1];
      break;
    case 27:
      at3 = params[i][1];
      break;
    case 28:
      at4 = params[i][1];
      break;
    case 29:
      ap1 = params[i][1];
      break;
    case 30:
      ap2 = params[i][1];
      break;
    case 31:
      ap3 = params[i][1];
      break;
    case 32:
      ap4 = params[i][1];
      break;
    default:
      cout<<"Error: Incorrect index "<<round(params[i][0])<<" for EOS constructor "<<phasetype<<endl;
      exit(1);
    };
  }

  if (eqntype<8)
  {
    cout<<"Error: Incorrect equation type. Type "<<eqntype<<" cannot be initialized by the constructor for RTpress style EOS."<<endl;
    exit(1);
  }

  // fitted polynomial parameters of the thermal coefficients b(V) in erg/mol.  Convert eV/atom to erg/mol need to multiply eV_erg*n*NA. For example, for MgSiO3, 0.9821 eV/atom = 4.824E12 *0.9821 erg/mol = 4.738E12 erg/mol.
  b = new double[bn];
  for (int i=0; i<bn; i++)
    b[i] = bparams[i];
}


EOS::~EOS()
{
  if(rhotable)
    delete[] rhotable;
  if(Ptable)
    delete[] Ptable;
  if(temptable)
    delete[] temptable;
  if(adiabattable)
    delete[] adiabattable;
  if(spline)
    gsl_spline_free (spline);
  if(accP)
    gsl_interp_accel_free (accP);
  if(spline2drho)
    gsl_spline2d_free (spline2drho);
  if(spline2dadi)
    gsl_spline2d_free (spline2dadi);
  if(accT)
    gsl_interp_accel_free (accT);
  if(bn>0)
    delete[] b;
}

void EOS::modifyEOS(double params[][2], int length)  // modify the constructed EOS parameters
{
  if (density_extern)		// Warning for EOS modification can't overwrite external EOS functions.
  {
    if (verbose)
      cout<<"Warning: External EOS function is set for state "<<phasetype<<". The modification on the EOS fitting parameters won't change the EOS function.  Consider change the external function if really need to modify this EOS."<<endl;
    return;
  }

  int index;
  
  for(int i=0;i<length;i++)
  {
    index = lround(params[i][0]);
    
    switch(index)
    {
    case 0:
      eqntype=round(params[i][1]);
      break;
    case 1:
      V0=params[i][1];
      break;
    case 2:
      K0=params[i][1];
      break;
    case 3:
      K0p=params[i][1];
      break;
    case 4:
      K0pp=params[i][1];
      break;
    case 5:
      mmol=params[i][1];
      break;
    case 6:
      P0=params[i][1];
      break;
    case 7:
      Theta0=params[i][1];
      break;
    case 8:
      gamma0=params[i][1];
      break;
    case 9:
      beta=params[i][1];
      break;
    case 10:
      gammainf=params[i][1];
      break;
    case 11:
      gamma0p=params[i][1];
      break;
    case 12:
      e0=params[i][1];
      break;
    case 13:
      g=params[i][1];
      break;
    case 14:
      if(eqntype>=8)		// RTpress style
      {
	cout<<"Error: The number of atoms of a RTpress style EOS is not allowed to be modified."<<endl;
	exit(1);
      }
      n=round(params[i][1]);
      break;
    case 15:
      Z=round(params[i][1]);
      break;
    case 16:
      T0=params[i][1];
      break;
    case 17:
      alpha0=params[i][1];
      break;
    case 18:
      alpha1=params[i][1];
      break;
    case 19:
      xi=params[i][1];
      break;
    case 20:
      cp_a=params[i][1];
      break;
    case 21:
      cp_b=params[i][1];
      break;
    case 22:
      cp_c=params[i][1];
      break;
    case 23:
      Debye_approx = params[i][1]>0 ? true : false;
      break;
    case 24:
      thermal_type = params[i][1];
      break;
    case 25:
      at1 = params[i][1];
      break;
    case 26:
      at2 = params[i][1];
      break;
    case 27:
      at3 = params[i][1];
      break;
    case 28:
      at4 = params[i][1];
      break;
    case 29:
      ap1 = params[i][1];
      break;
    case 30:
      ap2 = params[i][1];
      break;
    case 31:
      ap3 = params[i][1];
      break;
    case 32:
      ap4 = params[i][1];
      break;

    default:
      cout<<"Error: Incorrect index "<<round(params[i][0])<<" for EOS constructor "<<phasetype<<endl;
      exit(1);
    };
  }

  if (eqntype == 6)		// ideal gas
    thermal_type = 3;
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(gamma0)) // thermal type not specified and have enough information to calculate using Grueneisen parameter
  {
    if (gsl_finite(Theta0))
    {
      if (gsl_finite(e0) && gsl_finite(g)) // Also have enough information to get Pe
	thermal_type = 7;
      else
	thermal_type = 6;
    }
    else
      thermal_type = 5;
  }
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(cp(300)) && gsl_finite(alpha(10,300)) && gsl_finite(T0)) // thermal type not specified and have enough information to calculate using thermal expansion
    thermal_type = 9;

  
  if (eqntype >= 8)		// RTpress EOS style
    thermal_type = 8;
}

void EOS::modifyEOS(int index, double value)	     // modify one value of the EOS
{
  if (density_extern)		// Warning for EOS modification can't overwrite external EOS functions.
  {
    if (verbose)
      cout<<"Warning: External EOS function is set for state "<<phasetype<<". The modification on the EOS fitting parameters won't change the EOS function.  Consider change the external function if really need to modify this EOS."<<endl;
    return;
  }

  switch(index)
  {
  case 0:
    eqntype=round(value);
    break;
  case 1:
    V0=value;
    break;
  case 2:
    K0=value;
    break;
  case 3:
    K0p=value;
    break;
  case 4:
    K0pp=value;
    break;
  case 5:
    mmol=value;
    break;
  case 6:
    P0=value;
    break;
  case 7:
    Theta0=value;
    break;
  case 8:
    gamma0=value;
    break;
  case 9:
    beta=value;
    break;
  case 10:
    gammainf=value;
    break;
  case 11:
    gamma0p=value;
    break;
  case 12:
    e0=value;
    break;
  case 13:
    g=value;
    break;
  case 14:
    if(eqntype>=8)		// RTpress style
    {
      cout<<"Error: The number of atoms of a RTpress style EOS is not allowed to be modified."<<endl;
      exit(1);
    }
    n=round(value);
    break;
  case 15:
    Z=round(value);
    break;
  case 16:
    T0=value;
    break;
  case 17:
    alpha0=value;
    break;
  case 18:
    alpha1=value;
    break;
  case 19:
    xi=value;
    break;
  case 20:
    cp_a=value;
    break;
  case 21:
    cp_b=value;
    break;
  case 22:
    cp_c=value;
    break;
  case 23:
    Debye_approx = value>0 ? true : false;
    break;
  case 24:
    thermal_type = value;
    break;
  case 25:
    at1 = value;
    break;
  case 26:
    at2 = value;
    break;
  case 27:
    at3 = value;
    break;
  case 28:
    at4 = value;
    break;
  case 29:
    ap1 = value;
    break;
  case 30:
    ap2 = value;
    break;
  case 31:
    ap3 = value;
    break;
  case 32:
    ap4 = value;
    break;
  
  default:
    cout<<"Error: Incorrect index "<<index<<" for EOS constructor "<<phasetype<<endl;
    exit(1);
  };

  if (eqntype == 6)		// ideal gas
    thermal_type = 3;
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(gamma0)) // thermal type not specified and have enough information to using Grueneisen parameter
  {
    if (gsl_finite(Theta0))
    {
      if (gsl_finite(e0) && gsl_finite(g)) // Also have enough information to get Pe
	thermal_type = 7;
      else
	thermal_type = 6;
    }
    else
      thermal_type = 5;
  }
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(cp(300)) && gsl_finite(alpha(10,300)) && gsl_finite(T0)) // thermal type not specified and have enough information to calculate using thermal expansion
    thermal_type = 9;
  
  if (eqntype >= 8)		// RTpress EOS style
  {
    if (!gsl_finite(beta))
      beta = 0.6;
    thermal_type = 8;
  }
}

double EOS::BM3(double rho)
// input rho in g/cm^3, return pressure in GPa
{
  if (!gsl_finite(V0) || !gsl_finite(K0))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for BM3 EOS."<<endl;
    exit(1);
  }
  if (!gsl_finite(K0p))
    K0p=4;			// reduced to BM2
  
  double V = mmol/rho;

  double eta = V0/V;
  
  return P0 + 1.5*K0 * (pow(eta,7./3.)-pow(eta,5./3.)) * (1+0.75*(K0p-4)*(pow(eta,2./3.)-1));
}


double EOS::BM4(double rho)
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p) || !gsl_finite(K0pp))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for BM4 EOS."<<endl;
    exit(1);
  }
  double V = mmol/rho;

  double eta = V0/V;
  return P0 + 1.5*K0 * (pow(eta,7./3.)-pow(eta,5./3.)) * (1 + 0.75*(K0p-4)*(pow(eta,2./3.)-1) + 0.375*sq(pow(eta,2./3.)-1)*(K0*K0pp+K0p*(K0p-7)+143./9.));
}


double EOS::Vinet(double rho)
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for Vinet EOS."<<endl;
    exit(1);
  }

  double V = mmol/rho;

  double eta = V0/V;
  return P0 + 3*K0 * pow(eta,2./3.) * (1-pow(eta,-1./3.)) * exp(1.5*(K0p-1)*(1-pow(eta,-1./3.)));
}

double EOS::Holzapfel(double rho)
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p) || Z<0)
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for Holzapfel EOS."<<endl;
    exit(1);
  }
  
  double V = mmol/rho;

  double x = pow(V/V0,1./3.);
  double c0= -log(3*K0/(1003.6*pow(Z/V0,5./3.)));
  double c2= 1.5*(K0p-3)-c0;

  return 3*K0 * pow(x,-5) * (1-x) * exp(c0*(1-x)) * (1+c2*x*(1-x));
}

double EOS::Keane(double rho)
// Keane 1954, Aust. J. Phys. 7, 322-333
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p) || !gsl_finite(gammainf))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for Keane EOS."<<endl;
    exit(1);
  }

  double V = mmol/rho;

  double y = V0/V;
  double Kinfp = 2*(gammainf+1./6.);

  return K0p*K0/sq(Kinfp)*(pow(y,Kinfp)-1) - (K0p-Kinfp)*K0/Kinfp*log(y);
}

double EOS::H2OSC(double rho, double T)
// Mazevet
{
  constexpr double PI      = 3.141592653;
  constexpr double TWO_PI  = 2.0 * PI;
  constexpr double UN_T6   = 0.3157746;             // ha/kB × 1e−6
  constexpr double C13     = 1.0 / 3.0;
  constexpr double AUM     = 1822.88848;            // m_u / m_e
  constexpr double Zmean   = 10.0 / 3.0;
  constexpr double CMImean = 18.0 / 3.0;
  constexpr double DENSCONV = 11.20587 * CMImean;   // g cm⁻³ → n_i  [au]
  constexpr double TnkCONV  = 8.31447e13 / CMImean; // n_i kT factor [erg cm⁻³]
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
  const double RS  = pow(0.75 / PI / DENSI / Zmean, C13);
  const double GAME = 1.0 / (RS * TEMP);           // Γ_e
  const double GAMEsq = sqrt(GAME);
   
  // -------------------------------------------------------------------------
  //  1) super‑ionic / plasma contribution (most of the algebra is straight
  //     from the Fortran; only operator precedence and pow→std::pow changed)
  // -------------------------------------------------------------------------
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

  //  electron‑gas EOS (ideal Fermi gas) ------------------------------------
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

    // ---------------------------------------------------------------------
    //  2) non‑ideal molecular piece
    // ---------------------------------------------------------------------
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

// ---------------------------------------------------------------------
    //  3) mix molecular + plasma via smooth switch YL/YH
    // ---------------------------------------------------------------------
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
     // ---------------------------------------------------------------------
    //  4) Ideal gas of pseudo‑molecules
    // ---------------------------------------------------------------------
    const double THLmol = sqrt(TWO_PI / (18.0 * AUM * TEMP));
    const double FNkTid = (log(DENSMOL * pow(THLmol, 3)) - 1.0) / 3.0;
    const double PnkTid = C13;
    const double UNkTid = 0.5;

    // ---------------------------------------------------------------------
    //  5) total (non‑ideal + ideal)
    // ---------------------------------------------------------------------
    double PnkT,FNkT,UNkT,CV,CHIT,CHIR,PGPa,USPEC;
    FNkT = FNkTni + FNkTid;
    PnkT = PnkTni + PnkTid;
    UNkT = UNkTni + UNkTid;

    CV   = UNkT - FDTT;
    CHIR = FDRR / PnkT + 1.0;
    CHIT = FDRT / PnkT + 1.0;

    // ---------------------------------------------------------------------
    //  6) thermal corrections (high‑T & low‑T tweaks)
    // ---------------------------------------------------------------------
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

        // ---------------------------------------------------------------------
    //  7) auxiliary conversions
    // ---------------------------------------------------------------------
    const double Tnk = TnkCONV * rho * T6;
    PGPa = PnkT * Tnk / 1.0e10;   // → GPa
    USPEC = UNkT * Tnk / rho;      // erg g⁻¹
  
    double CP = CV + (CHIT*CHIT)/CHIR;

    double nabla_ad = CHIT / (CHIR * CP);

    double dT_dP_S  = nabla_ad * T / PGPa;

  return PGPa;
}


void EOS::DebyeT(double x, double &gamma, double &Theta)  // return the Grueneisen parameter, Debye temperature or Einstein temperature according to Altshuler form.
// If Theta0 is not available, a Debye temperature scaling factor is returned
{
  if ((!gsl_finite(V0) || thermal_type < 4) && !(thermal_type==2 && gsl_finite(gamma0) && gsl_finite(Theta0))) // don't have thermal pressure data
  {
    cout<<"Error: Cannot calculate the Debye temperature for phase "<<phasetype<<".  Lack of necessary information."<<endl;
    gamma = numeric_limits<double>::quiet_NaN();
    Theta = numeric_limits<double>::quiet_NaN();
    return;
  }

  // set up default parameters for gammainf, beta, and n if not provided by the input data.
  
  if (!gsl_finite(gammainf))	// According to Al'tshuler et al. 1987, gammainf = 2/3 for all elements except alkali elments, for which gammainf = 0.5
    gammainf = 2./3.;
  if (!gsl_finite(beta))	// Altshuler form.
    beta = gamma0 / (gamma0-gammainf);

  gamma = gammainf + (gamma0-gammainf)*pow(x,beta);

  if (thermal_type == 5 || thermal_type == 4)
    Theta = pow(x,-gammainf) * exp((gamma0-gammainf)/beta*(1-pow(x,beta))); 
  else
    Theta = Theta0 * pow(x,-gammainf) * exp((gamma0-gammainf)/beta*(1-pow(x,beta))); // gammma = - d ln(Theta) / d ln(V)
}


double EOS::Pth(double V, double T)
// calculate the thermal pressure in GPa, ref. Bouchet et al. 2013 PRB 87, 094102, Shim & Duffy, 2000, American Mineralogist
// Or RTpress style Pth.
{
  if (!gsl_finite(V0) || thermal_type < 2 || thermal_type == 3)  // don't have thermal pressure data
    return 0;

  if ((thermal_type == 2 && !gsl_finite(gamma0)) || thermal_type == 9)
    return 0;

  if (gsl_finite(T) && getthermal() == 8) // RTpress style
  {
    double cf = 1E-10;
// conversion factor from erg/cm^3 to GPa
    double T_OS = TOS(V);
    double fTp_OS = fTp(T_OS);
    double bVpV = bVp(V);	// in  erg/cm^3 (microbar)

    return - cf*bVpV*fT(T) + cf*gamma0S(V)*(T-T0)/V*Cv(V,T_OS) + cf*bVpV/(beta-1)*(T*(fTp(T)-fTp_OS) - T0*(fTp(T0)-fTp_OS));
  }

  // set up default parameters for n if not provided by the input data.
  if (n<0)
    n = 1;

  double gamma, Theta, x = V/V0;
  DebyeT(x, gamma, Theta);

  if (thermal_type == 4 || thermal_type == 5 || (thermal_type == 2 && !gsl_finite(Theta0)))	// Only gamma available
    // ref. Ono & Oganov 2005, Eq. 3, Earth Planet. Sci. Lett. 236, 914
    return 3E-10*gamma*n*R*(T-T0)/V;
  
  double Eth = 3*n*R*T*gsl_sf_debye_3(Theta/T);
  double Eth0 = 3*n*R*T0*gsl_sf_debye_3(Theta/T0);

  if (thermal_type == 7  || (thermal_type==2 && gsl_finite(e0) && gsl_finite(g)))
    return 1E-10*gamma*(Eth-Eth0)/V + 1.5E-16*n*R*e0*pow(x,g)*g/V*(sq(T)-sq(T0));
  // 1E-16 = 1E-10 * 1E-6.  GPa -> microbar and e0 in 10^-6 K^-1
  else // don't have enough information to get Pe
    return 1E-10 * gamma*(Eth-Eth0)/V;	  // convert to GPa
}

double EOS::adiabatic_index()	    // get the adiabatic index for ideal gas.  Vibrational freedom is always ignored.
{
  if (eqntype != 6 || n<0)
  {
    cout<<"Error: "<<phasetype<<" is not ideal gas or number of atom per molecule unknown.  No adiabatic index applied."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  if (n == 1)			// monatomic ideal gas
    return 5./3.;
  else if (n == 2)		//  diatomic gas and collinear molecules e.g. carbon dioxide
    return 1.4;
  else if (n == 0)		// isothermal atmosphere
    return 1.;
  else			// polyatomic gas
    return 4./3.;
}
  
double EOS::density(double P, double T, double rho_guess)
// input P in cgs (microbar), return density in g/cm^3
{
  if(!gsl_finite(P) || !gsl_finite(T)) // Check if P, T, or rho_guess is infinite or nan due to some error.  Stop code to avoid further error.
  {
    if (verbose)
      cout<<"Warning: Request density for "<<phasetype<<" at infinite/nan value.  P="<<P/1E10<<" T="<<T<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  int status;
  
  if(P < 0 || P > 1E16)		// unrealistic pressure
    return numeric_limits<double>::quiet_NaN();

  else if(density_extern)
    return density_extern(P, T);
  
  else if(eqntype == 7)		// interpolate an input file
  {

    P /= 1E10;
    double rho;
    if(tabletype == 4) //Search for Rho in 3D table need Temp and Press
    {
      status = gsl_spline2d_eval_e(spline2drho, T, P, accT, accP, &rho);
      
      if(status == GSL_EDOM)
      {
        if (verbose)
	        cout<<"Warning: Pressure "<<P<<"GPa or Temperature "<<T<<"K is outside the tabulated range for "<<this->phasetype<<". The density at the end point is returned"<<endl;
        if(P < Ptable[0] && T < temptable[0])
	        return rhotable[0];
        else if(P>Ptable[nline/tlen-1] && T>temptable[tlen-1])
          return rhotable[nline-1];
        else if(P<Ptable[0])
        {
          gsl_spline2d_eval_e(spline2drho, T, Ptable[0], accT, accP, &rho);
          return rho;
        }
        else if(P>Ptable[nline/tlen-1])
        {
          gsl_spline2d_eval_e(spline2drho, T, Ptable[nline/tlen-1], accT, accP, &rho);
          return rho;
        }
        else if(T < temptable[0])
        {
          gsl_spline2d_eval_e(spline2drho, temptable[0], P, accT, accP, &rho);
          return rho;
        }
        else
        {
          gsl_spline2d_eval_e(spline2drho, temptable[tlen-1], P, accT, accP, &rho);
          return rho;
        }          
      }
      else	
        return rho;      
    }

    else
    {
      status = gsl_spline_eval_e(spline, P, accP, &rho);

      if(status == GSL_EDOM)
      {
        if (verbose)
	        cout<<"Warning: Pressure "<<P<<"GPa is outside the tabulated range for "<<this->phasetype<<". The density at the end point is returned"<<endl;
        if(P < Ptable[0])
	        return rhotable[0];
        else
	        return rhotable[nline-1];
      }
      else	
        return rho;
    }
  }

  else if(eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
    return P*mmol*mp/(kb*T);
  }

  else
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    if (eqntype >= 8 && (!gsl_finite(n)||!gsl_finite(gamma0)||!gsl_finite(gamma0p)||!gsl_finite(V0)||!gsl_finite(beta)||!gsl_finite(T0)||!(bn>0)))
      cout<<"Error: Don't have enough input parameters to calculate the density of "<<phasetype<<" using RTpress style EOS."<<endl;

    P /= 1E10;			// convert pressure from microbar to GPa
    // if no temperature information, T should be numeric_limits<double>::quiet_NaN()
    struct EOS_params params = {{P, T}, this};

    if(rho_guess < 0.5 || !gsl_finite(rho_guess) || dP_EOS(rho_guess, &params) < 0)	// rho_guess will be set to negative if it is unknown. Ideal gas doesn't need a rho_guess.
      // if rho_guess is too small, dP/drho can be negative, and the solver may be tricked to the unphysical branch of the solution.
    {
      rho_guess = density(V0) + P/1E3;
    }

    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *TPL = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (TPL);
    gsl_function_fdf FDF;

    double rho = rho_guess, rho0;

    FDF.f = &P_EOS;
    FDF.df = &dP_EOS;
    FDF.fdf = &PdP_EOS;
    FDF.params = &params;

    gsl_root_fdfsolver_set (s, &FDF, rho);

    do
    {
      iter++;

      status = gsl_root_fdfsolver_iterate (s);
      rho0 = rho;
      rho = gsl_root_fdfsolver_root (s);
      if (rho<0.95*rho0)// limit the step size of each iteration to increase stability.
      {
	rho = 0.95*rho0;
	gsl_root_fdfsolver_set (s, &FDF, rho);
      }
      else if (rho>1.05*rho0)
      {
	rho = 1.05*rho0;
	gsl_root_fdfsolver_set (s, &FDF, rho);
      }

      status = gsl_root_test_delta (rho, rho0, 1E-16, rho_eps_rel);
    }
    while (status == GSL_CONTINUE && gsl_finite(rho) && iter < max_iter);

    if (!gsl_finite(rho))
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P<<" GPa and temperature "<<T<<" K, initial guessed rho:"<<rho_guess<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<". Likely no solution exist for this physical condition under the EOS used."<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }
    else if (status == GSL_CONTINUE)
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P<<" GPa and temperature "<<T<<" K within maximum interation "<<max_iter<<", initial guessed rho:"<<rho_guess<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }

    gsl_root_fdfsolver_free (s);

    if(thermal_type == 9 && T>T0)	// thermal expansion
      // convert alpha to K^-1
      return rho*exp(-1E-6*pow(1+K0p*P/K0, -xi)*(alpha0*(T-T0)+0.5*alpha1*(sq(T)-sq(T0))));
    else
      return rho;
  }
}

void EOS::printEOS()
// Create a tabulated EOS at  ./tabulated/phasename.txt
// The table from a pressure of 0.1GPa (or P0 if it is larger) to 2000 GPa at the temperature T0 (default 300 K) of the EOS.
{
  string filename = "./tabulated/" + phasetype + ".txt";
  ofstream fout(filename.c_str(),ofstream::trunc);
  if(!fout)
  {
    cout<<"Error: Failed to open file "<<filename<<" to save the EOS data"<<endl;
    return;
  }
  fout<<"Pressure (GPa)\t Density (g/cm^3)"<<endl;
  double rho_guess = 1, P = max(0.1*1e10,P0*1e10);
  for(int i = 0; i <= 100; i++ )
  {
    if(P > 2000*1e10)
      break;
    rho_guess = density(P,T0,rho_guess);
    fout << P/1e10 << "\t " << rho_guess << endl;
    P *= 1.1;
  }
}

double EOS::entropy(double rho, double T)
// Given the density and temperature, calculate the entropy over n*R, or P V^{7/5} / R = T V^0.4 for ideal gas.
// At either constant P or constant V, entropy always increases with temperature (y decreases with temperature).
{
  double gamma;
  if (eqntype == 6)		// ideal gas,  S ~ nR log(T^(1/(gamma-1))*V) + const.  For better performance and more concise code, here returns T rho^{1-gamma}.
  {
    gamma = this->adiabatic_index();
    return T*pow(rho,1-gamma);
  }

  if (entropy_extern)
    return entropy_extern(rho, T);

  double Theta, V = mmol/rho, x = V/V0;
   
  if (!gsl_finite(V0) || thermal_type < 5 || thermal_type == 8 || thermal_type == 9) // don't have thermal pressure data, or outside the range of the EOS, return negative number
    return -1;
  
  DebyeT(x, gamma, Theta);
  double y = Theta/T;

  double Sth;
  if (Debye_approx == true)		// Debye model
  {
    Sth = 4*gsl_sf_debye_3(y) - 3*log(1-exp(-y));
  }
  else				// Einstein model
  {
    Sth = 3 * (y/(exp(y)-1) - log(1-exp(-y)));
  }

  if (!gsl_finite(e0) || !gsl_finite(g)) // don't have enough information to get Pe
    return Sth;
  else
  {
    return Sth + 3E-6*e0*pow(x,g)*T;
  }
}

double P_EOS(double rho, void *params)
// function of EOS, given a rho (in cgs), return the difference between pressure from EOS and target P (in GPa).  Let this function equals 0 to solve for the correct rho.
{
  struct EOS_params *p = (struct EOS_params *) params;

  double P = p->x[0];
  double T = p->x[1];
  EOS* Phase = p->Phase;
  if (P < Phase->getP0())
  {
    if (verbose)
      cout<<"Warning: Incorrect phase diagram. "<<Phase->getEOS()<<" doesn't exist at pressure "<<P<<" GPa, which is smaller than its transition pressure at "<<Phase->getP0()<<" GPa. The transition pressure applied."<<endl;
    return Phase->getP0() - P;
  }

  return Phase->Press(rho, T) - P;
}

double dP_EOS(double rho, void *params)
{
  gsl_function F;
  double result, abserr;

  F.function = &P_EOS;
  F.params = params;
  gsl_deriv_central(&F, rho, 1E-4, &result, &abserr);
  return result;
}

void PdP_EOS(double rho, void *params, double *P, double *dP)
{
  *P=P_EOS(rho,params);
  *dP=dP_EOS(rho,params);
}

double EOS::gamma0S(double V)
// Grueneisen parameter along the reference adiabat (eq A.3), take volume in cm^3 / mol
{
  double a1 = 6*gamma0;
  double a2 = -12*gamma0 + 36*sq(gamma0) -18*gamma0p;
  double f = fV(V);
  return (2*f + 1) * (a1 + a2*f) / 6 / (1 + a1*f + 0.5*a2*sq(f));
}

double EOS::bV(double V)
// thermal coefficients b(V) in erg/mol (Eq.10), take volume in cm^3 / mol
{
  double sum=0;
  double Vdev = (V/V0-1);
  for (int i=0; i<bn; i++)
    sum += b[i]*pow(Vdev,i);
  return sum;
}

double EOS::bVp(double V)
// derivative of b(V) in erg/cm^3 (microbar) (Eq. B.2), take volume in cm^3 / mol
{
  double sum=0;
  double Vdev = (V/V0-1);
  for (int i=1; i<bn; i++)
    sum += b[i]*i/V0*pow(Vdev,i-1);
  return sum;
}

double EOS::TOS(double V)
// reference adiabat temperature profile (Eq. 7), take volume in cm^3 / mol
{
  double a1 = 6*gamma0;
  double a2 = -12*gamma0 + 36*sq(gamma0) -18*gamma0p;
  double f = fV(V);
  return T0 * sqrt(1+a1*f+0.5*a2*sq(f));
}

double EOS::Cv (double V, double T)
//total heat capacity in erg/mol/K (Eq. B.4), take volume in cm^3 / mol
{
  return bV(V)*fTp(T) + 1.5*n*R;
}

double EOS::Spot(double V, double T)
// potential contribution of entropy in erg/mol/K, take volume in cm^3 / mol
{
  return bV(V)/(beta-1)*fTp(T);
}

double EOS::gamma (double V, double T)
// Grueneisen parameter, take volume in cm^3 / mol
{
  if (!gsl_finite(gamma0))
  {
    cout<<"Error: gamma0 of phase "<<phasetype<<" is not available.  Thermal type is "<<thermal_type<<". Can't calculate Grueneisen parameter at V="<<V<<" T="<<T<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  
  if (eqntype >= 8)
  {
    // RTpress Eq. 17
    double T_OS = TOS(V);
    return gamma0S(V) * Cv(V,T_OS) / Cv(V,T) + V * bVp(V) / bV(V) * (Spot(V,T)-Spot(V,T_OS)) / Cv(V,T);
  }
  
  if (!gsl_finite(gammainf))	// According to Al'tshuler et al. 1987, gammainf = 2/3 for all elements except alkali elments, for which gammainf = 0.5
    gammainf = 2./3.;
  if (!gsl_finite(beta))	// Altshuler form.
    beta = gamma0 / (gamma0-gammainf);

  return gammainf + (gamma0-gammainf)*pow(V/V0,beta);
}

double EOS::cp(double T)
// specific heat capacity in J/g/K at constant pressure
{
  if (!gsl_finite(cp_a) && cp_b==0 && cp_c==0)
    return numeric_limits<double>::quiet_NaN();

  if (!gsl_finite(cp_a))
    return cp_b*T - cp_c/sq(T);
  else
    return cp_a + cp_b*T - cp_c/sq(T);
}

double EOS::alpha (double P, double T)
// coefficient of thermal expansion in K^-1. Input P in GPa, T in K
{
  if (!gsl_finite(alpha0) && alpha1==0)
    return numeric_limits<double>::quiet_NaN();

  else if (K0p<0 || K0<0 || P<0)
  {
    if (verbose)
      cout<<"Warning: thermal expansion of phase "<<phasetype<<" is not available because some physical parameters are negative, which is nonphysical.  Thermal type is "<<thermal_type<<"."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  
  double alphaP0;		// alpha at P=0
  if (!gsl_finite(alpha0))
  {
    alphaP0 = alpha1*T;
  }
  else
  {
    if (T>T0)
      alphaP0 = alpha0 + alpha1*T;
    else				// avoid alpha becomes negative at low temperature
      alphaP0 = alpha0 + alpha1*T0;
  }
  return 1E-6 * alphaP0 * pow(1+K0p*P/K0, -xi);
}

double EOS::Press(double rho, double T)
// pressure in GPa (Eq. 6, 13, 14) in Wolf&Bower 2018, take density in g/cm^3. For thermal expansion representation, this return the pressure at T0.
{
  double P;
  double V = volume(rho);	// volume in cm^3/mol

  switch(geteqntype())
  {
  case 0:			// BM3
  case 8:
    P = BM3(rho);
    break;
  case 1:			// BM4
  case 9:
    P = BM4(rho);
    break;
  case 2:			// Vinet
  case 10:
    P = Vinet(rho);
    break;
  case 3:			// Holzapfel
  case 11:
    P = Holzapfel(rho);
    break;
  case 4:			// Keane
  case 12:
    P = Keane(rho);
    break;
  case 6:			// ideal gas
    P = rho*kb*T/(mmol*mp);
    return P;
  case 5:
    P=H2OSC(rho, T);
    break;
  default:
    cout<<"Error: No such EOS type "<<geteqntype()<<" used in "<<phasetype<<endl;
    P = -1;
    exit(1);
  };

  if (gsl_finite(T) && (getthermal() > 4 || getthermal() == 2)) // Pth = 0 if getthermal == 9
// thermal pressure
    P+= Pth(V,T);

  return P;
}

double S_T(double V, void *params)
// entropy at constant T, volume in cm^3/mol
{
  struct EOS_params *p = (struct EOS_params *) params;
  
  EOS* Phase = p->Phase;
  double rho = Phase->density(V);
  double T = p->x[0];

  return Phase->entropy(rho, T);
}

double S_V(double T, void *params)
// entropy at constant V, volume in cm^3/mol
{
  struct EOS_params *p = (struct EOS_params *) params;

  EOS* Phase = p->Phase;
  double rho = Phase->density(p->x[0]);
  
  return Phase->entropy(rho, T);
}

double EOS::pSpV_T(double V, double T)
// partial S (entropy) partial V at constant T
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{T}, this};
  
  F.function = &S_T;
  F.params = &params;
  gsl_deriv_central(&F, V, 1E-4, &result, &abserr);
  return result;
}

double EOS::pSpT_V(double V, double T)
// partial S (entropy) partial T at constant V
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{V}, this};
  
  F.function = &S_V;
  F.params = &params;
  gsl_deriv_central(&F, T, 1E-2, &result, &abserr);
  return result;
}


double EOS::dTdV_S(double V, double P, double T)
  // adiabatic temperature gradient in K mol/cm^3, take volume in cm^3 / mol, P in GPa
{
  if (thermal_type == 1)	// has external entropy
    // dT/dV_S = - (dS/dV_T) / (dS/dT_V)
  {
    return - pSpV_T(V,T) / pSpT_V(V,T);
  }
  
  if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    double gamma = adiabatic_index();
    return (1-gamma)*T/V;
  }

  if (thermal_type == 9)	// thermal expansion
  {
    if (!gsl_finite(cp(300)) || !gsl_finite(alpha(10,300)) || !gsl_finite(mmol))
    {
      cout<<"Error: Information of phase "<<phasetype<<" is not enough to calculate temperature gradient using the thermal expension method."<<endl;
      return 0;
    }
    double a=alpha(P,T);
    // cp to GPa cm^3 g^-1 K^-1
    return (a*T*K0)/(mmol*(sq(a)*T*K0*V/mmol-1E-3*cp(T)));
  }
  
  return -gamma(V,T)*T/V;
}

double P_T(double rho, void *params)
// pressure at constant T in GPa, density in g/cm^3
{
  struct EOS_params *p = (struct EOS_params *) params;
  
  EOS* Phase = p->Phase;
  double T = p->x[0];

  return Phase->Press(rho, T);
}

double P_rho(double T, void *params)
// pressure at constant rho in GPa, density in g/cm^3
{
  struct EOS_params *p = (struct EOS_params *) params;

  EOS* Phase = p->Phase;
  double rho = p->x[0];
  
  return Phase->Press(rho, T);
}

double EOS::pPprho_T(double rho, double T)
  // partial P partial rho at constant T in GPa / g/cm^3
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{T}, this};
  
  F.function = &P_T;
  F.params = &params;
  gsl_deriv_central(&F, rho, 1E-4, &result, &abserr);
  return result;
}

double EOS::pPpT_rho(double rho, double T)
// partial P partial T at constant rho in GPa / K
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{rho}, this};
  
  F.function = &P_rho;
  F.params = &params;
  gsl_deriv_central(&F, T, 1E-2, &result, &abserr);
  // 1E-2 is the step size in conducting derivative. A step size too small will increase the error due to machine error. A value too large may also increase the error and may even crash the code.
  return result;
}


double EOS::dTdm(double m, double r, double rho, double P, double T)
  // adiabatic temperature gradient in K/g, P in cgs
{
  if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      return 0;
    }

    double gamma = adiabatic_index();
    return - (gamma-1)*mp*mmol*G*m/(gamma*kb*rho*4*pi*pow(r,4));
  }

  else if (thermal_type == 9)	// thermal expension
  {
    if (!gsl_finite(cp(300)) || !gsl_finite(alpha(10,300)))
    {
      cout<<"Error: Information of phase "<<phasetype<<" is not enough to calculate temperature gradient using the thermal expension method."<<endl;
      return 0;
    }

    return -1E-7*(alpha(P/1E10,T)*T*G*m)/(4*pi*pow(r,4)*rho*cp(T));
  }
  
  double V = volume(rho);
  double dTdV = dTdV_S(V, P/1E10, T);
  if (r<1)		// At the center of the planet where dTdm has a 0/0 limit
  {
    if (m>400 && verbose)
      cout<<"Warning: At the center of of the planet when conducting the first step ODE integration, the material density is "<<m*3./4./pi<<"g/cm^3, which seems to be too high."<<endl;
    return 0;
  }
  return 1E-10*dTdV*G*m/(4*pi*pow(r,4)) / (rho/V*pPprho_T(rho,T) - dTdV*pPpT_rho(rho,T));
}

double EOS::dTdP_S(double P, double T, double &rho_guess)
// partial T partial P along isentrope in K / GPa, given pressure in GPa
{
  if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    double gamma = adiabatic_index();
    return - (gamma-1)*T/(gamma*P);
  }
  else if(eqntype==5)
  {
    cout<<"hi"<<endl;
    return 0;
  }
  else if(eqntype==7) //tabular P-T table
  {
    P /= 1E10;
    double adiabat;
    int status;
    status = gsl_spline2d_eval_e(spline2dadi, T, P, accT, accP, &adiabat);
    //return 0;
    if(status == GSL_EDOM)
    {
      if(P < Ptable[0] && T < temptable[0]) //return end point if P or T outside table bounds
	      return adiabattable[0];
      else if(P>Ptable[nline/tlen-1] && T>temptable[tlen-1])
        return adiabattable[nline-1];
      else if(P<Ptable[0])
      {
        gsl_spline2d_eval_e(spline2dadi, T, Ptable[0], accT, accP, &adiabat);
        return adiabat/1E10;
      }
      else if(P>Ptable[nline/tlen-1])
      {
        gsl_spline2d_eval_e(spline2dadi, T, Ptable[nline/tlen-1], accT, accP, &adiabat);
        return adiabat/1E10;
      }
      else if(T < temptable[0])
      {
        gsl_spline2d_eval_e(spline2dadi, temptable[0], P, accT, accP, &adiabat);
        return adiabat/1E10;
      }
      else
      {
        gsl_spline2d_eval_e(spline2dadi, temptable[tlen-1], P, accT, accP, &adiabat);
        return adiabat/1E10;
      }          
      }
    else	
      return adiabat/1E10;     
  }
  rho_guess = density(P*1E10, T, rho_guess);
  double V = volume(rho_guess);
  double dTdV = dTdV_S(V, P, T);

  return -V*dTdV/(rho_guess*pPprho_T(rho_guess,T));
}

double P_EOS_S(double rho, void *params)
// function used to solve volume and temperature along adiabatic temperature profile with known temperature gradient.  Read in rho, return the difference between pressure from EOS and target P (in GPa).  Let this function equals 0 to solve for the correct rho.
{
  struct EOS_params *p = (struct EOS_params *) params;

  double P2 = p->x[0];
  double T1 = p->x[1];
  double rho1 = p->x[2];
  double dTdV = p->x[3];
  EOS* Phase = p->Phase;
  double V = Phase->volume(rho);
  double V1 = Phase->volume(rho1);
  
  return Phase->Press(rho, T1+(V-V1)*dTdV) - P2;
}

double dP_EOS_S(double rho, void *params)
{
  gsl_function F;
  double result, abserr;

  F.function = &P_EOS_S;
  F.params = params;
  gsl_deriv_central(&F, rho, 1E-4, &result, &abserr);
  return result;
}

void PdP_EOS_S(double rho, void *params, double *P, double *dP)
{
  *P=P_EOS_S(rho,params);
  *dP=dP_EOS_S(rho,params);
}

double EOS::density(double P1, double T1, double rho, double P2, double &T2)
// Given the pressure (cgs), temperature, density of the previous step, the pressure of the next step, return the temperature and density at the new pressure.  This solver doesn't conserve the entropy well enough. Only used as an approximation in the first integration step from the core of the planet where dTdm has 0/0 limit.
{
  if( !gsl_finite(P1) || !gsl_finite(P2) || !gsl_finite(T1) || !gsl_finite(rho)) // Check if P, the guess of T and rho is infinite or nan due to some error.  Stop code to avoid further error.
  {
    if (verbose)
      cout<<"Warning: Request density for "<<phasetype<<" at infinite/nan value.  P="<<P2/1E10<<" T="<<T1<<" rho_guess="<<rho<<endl;
    T2 = numeric_limits<double>::quiet_NaN();
    return numeric_limits<double>::quiet_NaN();
  }

  if(P2 < 0)
  {
    T2 = P2;
    return P2;
  }

  if (!entropy_extern && thermal_type < 3)
    // Using external density function, but no external entropy function. assuming isothermal.
// don't have thermal pressure data, isothermal applied.
  {
    T2 = T1;
    return density(P2, T2, rho);
  }

  else if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    double gamma = adiabatic_index();
    T2 = T1 * pow(P1/P2, (1-gamma)/gamma);
    return P2*mmol*mp/(kb*T2);
  }

  else
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
	
    if(rho < 0.5 || !gsl_finite(rho))		// rho will be set to negative if it is unknown.
      rho = density(V0) + P2/1E13;
  
    P1 /= 1E10;			// convert pressure from microbar to GPa
    P2 /= 1E10;

    if (eqntype >= 8 && (!gsl_finite(n)||!gsl_finite(gamma0)||!gsl_finite(gamma0p)||!gsl_finite(V0)||!gsl_finite(beta)||!gsl_finite(T0)||!(bn>0)))
    {
      cout<<"Error: Don't have enough input parameters to calculate the density of "<<phasetype<<" using RTpress style EOS."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }
  
    else if (thermal_type >=4 && thermal_type <8 && (!gsl_finite(V0) || !gsl_finite(mmol)||!gsl_finite(gamma0)))
    {
      cout<<"Error: Don't have enough input parameters to calculate the density of "<<phasetype<<" using Debye temperature method."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }
  
    else if (thermal_type ==1 && !entropy_extern)
    {
      cout<<"Error: Don't have user defined external entropy function to calculate the density of "<<phasetype<<"."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }

    else if (thermal_type == 9 && (!gsl_finite(cp_a) || !gsl_finite(alpha0) || !gsl_finite(T0) || !gsl_finite(mmol)))
    {
      cout<<"Error: Information of phase "<<phasetype<<" is not enough to calculate temperature gradient using the thermal expension method."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }
    
    double dTdV = dTdV_S(volume(rho), P1, T1);

    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *TPL = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (TPL);
    gsl_function_fdf FDF;

    double rho2 = rho, rho1=rho;

    struct EOS_params params = {{P2, T1, rho, dTdV}, this};

    FDF.f = &P_EOS_S;
    FDF.df = &dP_EOS_S;
    FDF.fdf = &PdP_EOS_S;
    FDF.params = &params;
  
    int status;
    gsl_root_fdfsolver_set (s, &FDF, rho2);
    
    do
    {
      iter++;

      status = gsl_root_fdfsolver_iterate (s);
      rho1 = rho2;
      dTdV = dTdV_S(volume(rho1), P1, T1);
      params.x[3] = dTdV;
      rho2 = gsl_root_fdfsolver_root (s);
      if (rho2<0.95*rho1)// limit the step size of each iteration to increase stability.
      {
	rho2 = 0.95*rho1;
	gsl_root_fdfsolver_set (s, &FDF, rho2);
      }
      else if (rho2>1.05*rho1)
      {
	rho2 = 1.05*rho1;
	gsl_root_fdfsolver_set (s, &FDF, rho2);
      }

      T2 = T1 + mmol/sq(rho2)*(rho-rho2)*dTdV; // 
      status = gsl_root_test_delta (rho1, rho2, 1E-16, rho_eps_rel);
    }
    while (status == GSL_CONTINUE && gsl_finite(rho) && iter < max_iter);

    if (!gsl_finite(rho2))
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P2<<" GPa and temperature "<<T1<<" K, initial guessed density:"<<rho<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<". Likely no solution exist for this physical condition under the EOS used."<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }
    else if (status == GSL_CONTINUE)
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P2<<" GPa and temperature "<<T1<<" K within maximum interation "<<max_iter<<", initial guessed density:"<<rho<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }
    gsl_root_fdfsolver_free (s);
    return rho2;
  }
}

// =============================================================================
//  ELECNR – ideal, non‑relativistic electron Fermi gas
// =============================================================================


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

// =============================================================================
//  FERMIF: simple logistic Fermi function 1/(exp(x)+1)
// =============================================================================

double FERMIF(double X)
{
    if (X > 40.0)  return 0.0;
    if (X < -40.0) return 1.0;
    return 1.0 / (exp(X) + 1.0);
}

// =============================================================================
//  FERINT – Antia (1993) rational fits to F_q(x)
//            q = N-1/2 with N=0..3 (−½, ½, 3/2, 5/2)
// =============================================================================

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
    static const int LB[4] = {7,7,7,7};
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
    static const int LD[4] = {11,11,10,9};

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

// =============================================================================
//  FINVER – Antia (1993) fits to inverse Fermi integrals X_q(f)
// =============================================================================

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
