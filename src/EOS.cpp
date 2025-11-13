#include "EOS.h"

EOS::EOS():phasetype(""),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), a_vdW(numeric_limits<double>::quiet_NaN()),b_vdW(numeric_limits<double>::quiet_NaN()),Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), temptable(NULL), adiabattable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  shA = shB = shC = shD = shE = std::numeric_limits<double>::quiet_NaN();
  has_shomate = false;
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
}

EOS::EOS(string phaseinput, double params[][2], int length):phasetype(phaseinput),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), a_vdW(numeric_limits<double>::quiet_NaN()),b_vdW(numeric_limits<double>::quiet_NaN()),Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), temptable(NULL), adiabattable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  shA = shB = shC = shD = shE = std::numeric_limits<double>::quiet_NaN();
  has_shomate = false;
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
    case 33:
      a_vdW = params[i][1];
      break;
    case 34:
      b_vdW = params[i][1];
      break;
    case 35: shA = params[i][1]; has_shomate = true; break;
    case 36: shB = params[i][1]; has_shomate = true; break;
    case 37: shC = params[i][1]; has_shomate = true; break;
    case 38: shD = params[i][1]; has_shomate = true; break;
    case 39: shE = params[i][1]; has_shomate = true; break;
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

EOS::EOS(string phaseinput, string filename):phasetype(phaseinput),eqntype(7), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), a_vdW(numeric_limits<double>::quiet_NaN()),b_vdW(numeric_limits<double>::quiet_NaN()),Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), temptable(NULL), adiabattable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  shA = shB = shC = shD = shE = std::numeric_limits<double>::quiet_NaN();
  has_shomate = false;
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
    thermal_type=10;
   
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


EOS::EOS(string phaseinput, double (*f)(double P, double T, double rho_guess), double (*g)(double rho, double T)):phasetype(phaseinput),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), a_vdW(numeric_limits<double>::quiet_NaN()),b_vdW(numeric_limits<double>::quiet_NaN()),Debye_approx(false), rhotable(NULL), Ptable(NULL), bn(0), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  shA = shB = shC = shD = shE = std::numeric_limits<double>::quiet_NaN();
  has_shomate = false;
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

EOS::EOS(string phaseinput, double *Plist, double *rholist, int len_list):phasetype(phaseinput),eqntype(7), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), a_vdW(numeric_limits<double>::quiet_NaN()),b_vdW(numeric_limits<double>::quiet_NaN()),Debye_approx(false), thermal_type(0), bn(0), accT(NULL), spline2drho(NULL), spline2dadi(NULL), nline(len_list), tlen(0)
{
  shA = shB = shC = shD = shE = std::numeric_limits<double>::quiet_NaN();
  has_shomate = false;
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

EOS::EOS(string phaseinput, double params[][2], double bparams[], int length, int blength):phasetype(phaseinput),eqntype(8), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), a_vdW(numeric_limits<double>::quiet_NaN()),b_vdW(numeric_limits<double>::quiet_NaN()),Debye_approx(false), thermal_type(8), rhotable(NULL), Ptable(NULL), bn(blength), accP(NULL), accT(NULL), spline(NULL), spline2drho(NULL), spline2dadi(NULL), nline(0), tlen(0)
{
  shA = shB = shC = shD = shE = std::numeric_limits<double>::quiet_NaN();
  has_shomate = false;
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
    case 33:
      a_vdW = params[i][1];
      break;
    case 34:
      b_vdW = params[i][1];
      break;
    case 35: shA = params[i][1]; has_shomate = true; break;
    case 36: shB = params[i][1]; has_shomate = true; break;
    case 37: shC = params[i][1]; has_shomate = true; break;
    case 38: shD = params[i][1]; has_shomate = true; break;
    case 39: shE = params[i][1]; has_shomate = true; break;
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
    case 33:
      a_vdW = params[i][1];
      break;
    case 34:
      b_vdW = params[i][1];
      break;
    case 35: shA = params[i][1]; has_shomate = true; break;
    case 36: shB = params[i][1]; has_shomate = true; break;
    case 37: shC = params[i][1]; has_shomate = true; break;
    case 38: shD = params[i][1]; has_shomate = true; break;
    case 39: shE = params[i][1]; has_shomate = true; break;
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
  case 33:
    a_vdW = value;
    break;
  case 34:
    b_vdW = value;
    break;
  case 35: shA = value; has_shomate = true; break;
  case 36: shB = value; has_shomate = true; break;
  case 37: shC = value; has_shomate = true; break;
  case 38: shD = value; has_shomate = true; break;
  case 39: shE = value; has_shomate = true; break;
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


double EOS::vdW_gas(double rho, double T)
// Keane 1954, Aust. J. Phys. 7, 322-333
{
  double V = mmol/rho;
  
  if (!gsl_finite(a_vdW) || !gsl_finite(b_vdW))
    // ideal gas
  {
    return rho*kb*T/(mmol*mp);
  }

  return 1E-10*R*T/(V*(1-1E3*b_vdW/V))-1E2*a_vdW/sq(V);
}


void EOS::DebyeT(double x, double &gamma, double &Theta)  // return the Grueneisen parameter, Debye temperature or Einstein temperature according to Altshuler form.
// If Theta0 is not available, a Debye temperature scaling factor is returned
{
  if ((!gsl_finite(V0) || thermal_type < 4) && !(thermal_type==2) && !(gsl_finite(gamma0) && gsl_finite(Theta0))) // don't have thermal pressure data
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

double EOS::adiabatic_index() const
// get the adiabatic index for ideal gas.  Vibrational freedom is always ignored.
// Note: For van der Waals gas, the adiabatic index C_p/C_v is no longer a constant.
// But Cv is still relates to the internal degrees of freedom of the gas molecules (translational, rotational, vibrational).
// So we still have Cv = R / (gamma-1)
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

double EOS::gamma_shomate(double T) const
{
  // NIST/JANAF Shomate: Cp°(T) = A + B t + C t^2 + D t^3 + E / t^2,  t = T/1000
  // Returns gamma(T) = Cp / (Cp - R), with R molar gas constant.
  static const double Rm = 8.314462618; // J/mol/K
  if (!has_shomate || !gsl_finite(T) || T <= 0.0)
    return adiabatic_index(); // fallback to existing bucketed gamma

  const double t = T / 1000.0;
  const double Cp = shA + shB*t + shC*t*t + shD*t*t*t + shE/(t*t); // J/mol/K
  const double Cv = Cp - Rm;
  // Guard against pathological near-critical regions or bad coeffs:
  if (!gsl_finite(Cp) || !gsl_finite(Cv) || Cv < 1)
    return adiabatic_index();
  // fallback to existing bucketed gamma
  double g = Cp / Cv;
  if (!gsl_finite(g))
    return adiabatic_index();
  // fallback to existing bucketed gamma
  if (g < 1.01) g = 1.01;
  if (g > 1.67) g = 1.67;
  return g;
}
  
double EOS::density(double P, double T, double rho_guess)
// input P in cgs (microbar), return density in g/cm^3
{
  if(!gsl_finite(P) || !gsl_finite(T)) // Check if P or T is infinite or nan due to some error.  Stop code to avoid further error.
  {
    if (verbose)
      cout<<"Warning: Request density for "<<phasetype<<" at infinite/nan value.  P="<<P/1E10<<" T="<<T<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  int status;
  
  if(P < 0 || P > 1E16)		// unrealistic pressure
    return numeric_limits<double>::quiet_NaN();

  else if(density_extern)
    return density_extern(P, T, rho_guess);
  
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

  else if(eqntype == 6)
// Van der Waals EoS (P  + a(n/V)^2) (1 - b(n/V)) = (n/V) R T
    // a = 27 (R Tc)^2 / (64 Pc), b = R Tc / (8 Pc)
    // Tc and Pc available in "The properties of gases and liquids 5th edition" Poling, Prausnitz, O'Connel, Appendix A
    // a and b from CRC Handbook of Chemistry and Physics, section Fluid Properties, Lide & Haynes.
    // rho = (n/V)*mu
    // The equation have single solution in most scenario, but may have three solutions when P<Pc and T<Tc
    // It's too complicate to write the analytic solution for all situation, use root solver instead
  {
    if (!gsl_finite(mmol))
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      exit(1);
    }
    else if (!gsl_finite(a_vdW) || !gsl_finite(b_vdW))
      // ideal gas
      return P*mmol*mp/(kb*T);
    else			// Van der Waals EOS
    {
      if(P>R*T/b_vdW)
        // the volume is too close to b that the root solver would fail
        return 1E-3*mmol/b_vdW;

      P /= 1E10;			// convert pressure from microbar to GPa
      // if no temperature information, T should be numeric_limits<double>::quiet_NaN()

      struct EOS_params params = {{P, T}, this};

      if(rho_guess < 1E-5 || !gsl_finite(rho_guess) || dP_EOS(rho_guess, &params) < 0)	// rho_guess will be set to negative if it is unknown.
        rho_guess = 1E-3*mmol/b_vdW - 1E-4;
        // slightly less than the maximum density.
        // the 1E-4 guarantee the density would fall in the larger side, which would crash the code


      int iter = 0, max_iter = 200;
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
        if (rho<0.8*rho0)// limit the step size of each iteration to increase stability.
        {
          rho = 0.8*rho0;
          gsl_root_fdfsolver_set (s, &FDF, rho);
        }
        else if (rho>1.2*rho0)
        {
          rho = 1.2*rho0;
          gsl_root_fdfsolver_set (s, &FDF, rho);
        }

        status = gsl_root_test_delta (rho, rho0, 1E-16, rho_eps_rel);
      }
      while (status == GSL_CONTINUE && gsl_finite(rho) && iter < max_iter);

      if (!gsl_finite(rho))
      {
        if (verbose)
          cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P<<" GPa and temperature "<<T<<" K, initial guessed rho:"<<rho_guess<<endl;
      
        gsl_root_fdfsolver_free (s);
        return numeric_limits<double>::quiet_NaN();
      }
      else if (status == GSL_CONTINUE)
      {
        if (verbose)
          cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P<<" GPa and temperature "<<T<<" K within maximum interation "<<max_iter<<", initial guessed rho:"<<rho_guess<<endl;
      
        gsl_root_fdfsolver_free (s);
        return numeric_limits<double>::quiet_NaN();
      }

      gsl_root_fdfsolver_free (s);

      return rho;
    }
  }

  else
  {
    if (!gsl_finite(mmol))
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      exit(1);
    }

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
  if (eqntype == 6)
// ideal gas,  S ~ nR log(T^(1/(gamma-1))*V) + const.  For better performance and more concise code, here returns T rho^{1-gamma}.
    // van der Waals, returns T (V-b)^{gamma-1}.
  {
    double gamma = (has_shomate ? gamma_shomate(T) : this->adiabatic_index());
    if (!gsl_finite(a_vdW) || !gsl_finite(b_vdW))
      return T*pow(rho,1-gamma);
    else
    {
      double V = volume(rho);
      return T*pow(V-1E3*b_vdW, gamma-1);
    }
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
  case 6:			// van der Waals gas
    P = vdW_gas(rho, T);
    return P;
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
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      exit(1);
    }

    double gamma = (has_shomate ? gamma_shomate(T) : adiabatic_index());
    if (!gsl_finite(a_vdW) || !gsl_finite(b_vdW))
      return (1-gamma)*T/V;
    else
      return (1-gamma)*T/(V-1E3*b_vdW);
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
  if (eqntype == 6)		// gas
  {
    if (!gsl_finite(mmol))
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      exit(1);
    }

    double gamma = (has_shomate ? gamma_shomate(T) : adiabatic_index());
    if (!gsl_finite(a_vdW) || !gsl_finite(b_vdW))
      return -(gamma-1)*mp*mmol*G*m/(gamma*kb*rho*4*pi*pow(r,4));
    else
    {
      double V = volume(rho);
      return -(gamma-1)*T*G*m / (4*pi*pow(r,4)*(P*gamma +(gamma-2)*1E2*a_vdW/sq(V) + 2E5*a_vdW*b_vdW/pow(V,3)));
    }
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
  
  else if (thermal_type == 2) 	// External temperature gradient function
  {
    double dPdm =  -G*m/(4*pi*pow(r,4));
    return dPdm * dTdP(P, T, rho);
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
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      exit(1);
    }

    double gamma = (has_shomate ? gamma_shomate(T) : adiabatic_index());
    if (!gsl_finite(a_vdW) || !gsl_finite(b_vdW))
      return (gamma-1)*T/(gamma*P);
    else
    {
      double rho = density(P,T,rho_guess);
      double V = volume(rho);
      return (gamma-1)*T / (P*gamma +(gamma-2)*1E2*a_vdW/sq(V) + 2E5*a_vdW*b_vdW/pow(V,3));
    }
  }
  if(eqntype==7) //tabular P-T table
  {
    P /= 1E10;
    double adiabat;
    int status;
    status = gsl_spline2d_eval_e(spline2dadi, T, P, accT, accP, &adiabat);

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

  if (thermal_type == 2)				// External temperature gradient function
  {
    T2 = T1 + dTdP(P1, T1, rho)*(P2-P1);
    return density(P2, T2, rho);
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
    
    // Use dTdP_S for adiabatic temperature gradient calculation
    if(rho < 0.01 || !gsl_finite(rho))          // rho will be set to negative if it is unknown.
      rho = density(V0) + P2/1E13;
    
    T2 = T1 + dTdP_S(P1/1E10, T1, rho) * (P2/1E10 - P1/1E10);
    return density(P2, T2, rho);
  }
}

// Define a structure to hold the parameters
typedef struct
{
  double P;
  double T;
  double (*pressure_func)(double, double);
} density_params;

// Define the function to solve for the root
double fP_EOS(double rho, void *params)
{
  density_params *p = (density_params *)params;
  return p->pressure_func(rho, p->T) - p->P; // f(rho) = Pressure(rho, T) - P
}

// Define the derivative of the function using gsl_deriv_central
double fdP_EOS(double rho, void *params)
{
  // Use gsl_deriv_central to compute the derivative
  double result, abserr;
  gsl_function F;
  F.function = &fP_EOS; // Directly use f as the function
  F.params = params;
  gsl_deriv_central(&F, rho, 1e-4, &result, &abserr);

  return result;
}

// Define the function and its derivative together
void fdfP_EOS(double rho, void *params, double *f, double *df)
{
  density_params *p = (density_params *)params;
  *f = p->pressure_func(rho, p->T) - p->P; // f(rho) = Pressure(rho, T) - P
  *df = fdP_EOS(rho, params); // Use the df function defined above
}

double density_solver(double P, double T, double (*pressure_func)(double rho, double T), double rho_guess)
// Given the function pressure_func, which calculates pressure at (rho, T), use a root solver to find the density (rho) at the given pressure (P) in microbar and temperature (T).
{
  if(!gsl_finite(P) || !gsl_finite(T)) // Check if P or T is infinite or nan due to some error.  Stop code to avoid further error.
  {
    if (verbose)
      cout<<"Warning: Request density at infinite/nan value.  P="<<P/1E10<<" T="<<T<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  
  if(P < 0 || P > 1E16)		// unrealistic pressure
    return numeric_limits<double>::quiet_NaN();

  P /= 1E10;			// convert pressure from microbar to GPa

  density_params params[3] = {P, T, pressure_func};

  int status;
  int iter = 0, max_iter = 100;
  
  const gsl_root_fdfsolver_type *TPL = gsl_root_fdfsolver_newton;
  gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (TPL);
  gsl_function_fdf FDF;

  double rho = rho_guess, rho0;

  FDF.f = &fP_EOS;
  FDF.df = &fdP_EOS;
  FDF.fdf = &fdfP_EOS;
  FDF.params = &params;

  gsl_root_fdfsolver_set (s, &FDF, rho);

  do
  {
    iter++;

    status = gsl_root_fdfsolver_iterate (s);
    rho0 = rho;
    rho = gsl_root_fdfsolver_root (s);
    if (rho<0.8*rho0)// limit the step size of each iteration to increase stability.
    {
      rho = 0.8*rho0;
      gsl_root_fdfsolver_set (s, &FDF, rho);
    }
    else if (rho>1.2*rho0)
    {
      rho = 1.2*rho0;
      gsl_root_fdfsolver_set (s, &FDF, rho);
    }

    status = gsl_root_test_delta (rho, rho0, 1E-16, rho_eps_rel);
  }
  while (status == GSL_CONTINUE && gsl_finite(rho) && iter < max_iter);

  if (!gsl_finite(rho))
  {
    if (verbose)
      cout<<"Warning: Can't find the density at pressure "<<P<<" GPa and temperature "<<T<<" K, initial guessed rho:"<<rho_guess<<" for user set up function.  Likely no solution exist for this physical condition under the EOS used."<<endl;
      
    gsl_root_fdfsolver_free (s);
    return numeric_limits<double>::quiet_NaN();
  }
  else if (status == GSL_CONTINUE)
  {
    if (verbose)
      cout<<"Warning: Can't find the density at pressure "<<P<<" GPa and temperature "<<T<<" K within maximum interation "<<max_iter<<", initial guessed rho:"<<rho_guess<<endl;
      
    gsl_root_fdfsolver_free (s);
    return numeric_limits<double>::quiet_NaN();
  }

  gsl_root_fdfsolver_free (s);

  return rho;
}

