# MAGRATHEA #

![Header](plot/magratheahead.jpg)

**Excerpt from The Hitchhiker's Guide to the Galaxy, Page 634784, Section 5a, Entry: MAGRATHEA**

Planet interior structure code for astronomers, planetary scientists, mice, and more. 


## What is this repository for? ##

A 1D planet structure code written in C++ which considers the case of fully differentiated interiors. 
The code integrates the hydrostatic equation in order to shoot for the correct planet radius given the mass in each layer.
The code returns the pressure, temperature, density, phase, and radius at steps of enclosed mass.
The code supports 4 layers: core, mantle, hydrosphere, and atmosphere. Each layer has a phase diagram with equations of state (EoS) chosen for each phase.
The code was developed by [Chenliang Huang](https://huang-cl.github.io/), [David R. Rice](https://davidrrice.github.io/), and [Jason H. Steffen](https://www.jasonhsteffen.com/) at the Univerisity of Nevada, Las Vegas starting in 2017.
See a [list of works](citations.md) that use MAGRATHEA and [instructions on how to cite](CITATION.md). **If you don't see something on this ReadMe check our publication in this repository: [*MAGRATHEA.pdf*](MAGRATHEA_publication.pdf).** A tutorial and practice projects/problems are found in [Tutorial_Practice_Problems.pdf](Tutorial_Practice_Problems.pdf)

We encourage the community to contribute to and use MAGRATHEA for their interior modeling needs.

 <p align="center">
<img width = "600" src="plot/funplanets.png"/>
 </p>

## Prerequisite ##

[Install the GSL library](https://www.gnu.org/software/gsl/)(>= v2.0).  Download the compressed package from the [GNU ftp site](ftp://ftp.gnu.org/gnu/gsl/).  Extract the file and install the package following the instruction in the `INSTALL` file.  For the simplest case, 
 
    sudo ./configure
    sudo make
    sudo make install

should configure, build, and install the GSL package.  A few prerequisites, such as `g++`, `make`, `make-guile`, may need to be installed following the error messages throughout the process. 

On Ubuntu systems, the GSL package can also be installed from the Ubuntu repository using `sudo apt install libgsl27 libgsl-dev gsl-bin`. 

On Windows systems, we suggest using WSL and following the above instructions (some file paths within the code may break in Windows, please use WSL).

When running many simulations with bulk input or composition finder modes, OpenMP can be used to run individual planets in parallel, please contact us to get this set up in the code. 

## Quick Start ##

Installation:

1. Open a terminal. Navigate to where you wish to install the directory.
2. Clone the repository: `git clone https://github.com/Huang-CL/Magrathea.git`.
3. If the gsl library is not installed globally (under /usr/local/ or equivalent), edit `Makefile` of MAGRATHEA to include the actual path toward the gsl headers (e.g. `~/include` or `~/gsl/include`) and gsl library (e.g. `~/lib` or `~/gsl/lib`) following `-I` in `CFLAGS` and `-L` in `LDFLAGS`. The path of headers and libraries can be found using `gsl-config --cflags` and `gsl-config --libs`.  If the `gsl-config` command is not found, `gsl` may not have been installed properly.  _Note: The path after `-I` should end with `/include`, not `/include/gsl`._
4. Run `make -B` inside the code's directory. _After the first compilation use `make` to compile the code._

Running a planet:

5. Open `run/mode0.cfg` in a text editor.
6. Line 12-14: set mass of the core, mantle, hydrosphere, and atmosphere in Earth-masses.
7. Line 16-20: set surface temperature and any temperature discontinuities.
8. Line 21: set output file name and path.
9. In terminal, compile changed file with `make`.
10. Run MAGRATHEA with `./planet run/mode0.cfg`.
11. If an error message like "error while loading shared libraries: libgsl.so.23: cannot open shared object file: No such file or directory" is reported when running the code, add `export LD_LIBRARY_PATH=/usr/local/lib` (directory of GSL library files) to the `.bashrc` file, or add `setenv LD_LIBRARY_PATH /usr/local/lib` to the `.cshrc` file.  To apply the changes immediately, execute either `source ~/.bashrc` or `source ~/.cshrc`.

## Capability and Output ##

MAGRATHEA uses a shooting to a fitting-point method with a Runge-Kutta-Fehlberg method with an adaptive step-size. The user's supplied mass fractions determine the enclosed mass at each layer's boundary. Within a layer the enclosed mass at which the phase changes is determined by P-T conditions. When a phase changes, the solver backs up to find the exact location of the phase change. The fitting point is at 20% the mass of the planet. The solver integrates until the inner and outer branch of integration agree at the fitting point with a relative error less than 10<sup>-4</sup>.

### Input Modes ##

**There are 9 modes with different functionality described below. These modes can be used with the 9 config (.cfg) files in the `run` directory. The code is then ran with `./planet run/mode_.cfg`.**

### mode0.cfg ###
The basic capability: calculate the structure of a planet given the mass of each layer.

The config file requires defining the mass in the planet's core, mantle, hydrosphere, and atmosphere in the unit of Earth-mass.  Additionally the mode requires four temperatures, first the "surface temperature" in Kelvin which is the temperature at the top of the planet, where enclosed mass is equal to total mass. Then three temperatures which will be used for discontinuities in temperature at layer boundaries. The first being `temp_jump_1`, the size of temperature discontinuity between the atmosphere and hydrosphere, and the last being `temp_jump_3`, the size of temperature discontinuity between the core and the mantle. Setting these jumps to 0 will ensure layers in thermodynamic equilibrium. These parameters will be used in the `fitting_method()` function in `src/main.cpp`.

The function also requires a guess for the density in each layer from `ave_rho` and the surface pressure from `P_surface` located in the `#Global Run Options` section of the config file. Default surface pressure is 100 mbar. Lastly, the phase diagrams are chosen. This is detailed further in the _Phase Diagram_ section below.

The planet structure will be output as an ASCII file with pressure (GPa), interior mass (Earth mass), density (g cm<sup>-3</sup>), temperature (K), and phase of composition as a function of radius (Earth radius).

### mode1.cfg ###

Running the code with this input file will use the `getmass()` function which assumes an isothermal planet. This function runs much quicker than the full solver. User provides in the config file the mass of the core, mantle, and hydrospehre in Earth-masses. No atmosphere layer is available in this solver. Outputs interior conditions to a file as in mode 0.

### mode2.cfg ###

The fastest solver available using the `twolayer()` function. Using this input file will find radii for an array of planet masses each with only two layers using an isothermal inside-out shooting method.  Built to quickly make mass-radius curves with a constant mass ratio between the layers. 

The config file requires defining a `layer_index` where =0 is a planet with only mantle and hydrosphere (no core), =1 is only core and hydrosphere, and =2 is core and mantle. The file also requires defining `mass_fraction` which is the mass fraction of the inner layer. An array of masses in Earth-masses is defined by `min_mass`, `max_mass`, and `step_mass`. Planets will be solved from `min_mass` to `max_mass` in steps of `step_mass`. Returns radii for each mass in the terminal.


### mode3.cfg ###

This mode is for bulk inputs of planets using the solvers from mode 0 or mode 1. Requires a space separated file, where each row of the table lists the total mass in Earth-masses and fraction of mass in each layer for a planet.

Example input file:

    Mass  fCore  fMantle  fWater
    2     0.2    0.4      0.4
    1.5   0.5    0.39     0.1
	
Any remaining mass (1-fCore-fMantle-fWater) will be put into the atmosphere.  Example input files can also be found in the `input` directory.

In the config file, the user sets the path to the input file and the path to where the output file should be created. The full solver from mode 0 can be used by setting `solver=1` and the temperature-free solver can be used with `solver=2` (solver=1 is recommended). For `solver=1`, the surface temperature and temperature jumps are defined in the same way as mode 0. These parameters will be used in the `multiplanet()` function from `src/compfind.cpp`.

MAGRATHEA will generate an output file with mass of core, mantle, water, and atmosphere and the radius of the core, mantle, water, and planet for each line in the input file.  If the solution crosses the part of the phase diagram that the code has not been fully implemented, "Dummy EOS used" is added to the end of the row.  If the solver cannot find the solution that meets the required accuracy with the given number of iteration, "Solution did not converge" is added.  "No solution found" is also possible in the output if the solver failed to find a solution.

### mode4.cfg ###

This mode is our composition finder which finds the mass fraction of an unkown layer for a planet of given mass and radius. Requires a space separated file, where each row of the table lists a planet's total mass in Earth-masses and radius in Earth-radii.

Example input file:

    M (Earth-masses) 	R (Earth-radii)
    1.06658             1.09393
    1.02955             1.04153
	
Example input files can also be found in the `input` directory.

In the config file, the user sets the path to the input file and the path to where the output file should be created. Then they set which layer to find, `find_layer` and which two layers to set in constant **partial mass ratio**, `layer_outer` and `layer_inner`. The indexes are 1 for core, 2 for mantle, 3 for water/hydrosphere, and 4 for atmosphere. The index of `layer_inner` must be less than that of `layer_outer`. 

The **partial mass ratio (PMR%)** is given by the percentage ratio of the outer mass fraction to the total non-unkown mass fraction. In other words, OMF/(IMF+OMF)*100 where OMF is the outer mass fraction and IMF is the inner mass fraction. The user can then make an array of PMRs to loop over from `PMR_min` to `PMR_max` with a step of `PMR_step` these parameters must be between 0 and 100 and truncated to the tenths place. If `PMR_min` is equal to `PMR_max` the composition finder will find the unkown mass fraction for just one PMR for each mass and radius. Lastly, the user must define `R_error` which is error tolerance in simulated radius to target radius.

The surface temperature and temperature jumps are defined in the same way as mode 0. The above parameters will be used in the `compfinder()` function from `src/compfind.cpp`.

MAGRATHEA will use the full solver from mode 0 and for each mass and target radius in the input file and will first find the radius of a planet with no unkown mass and the rest of the mass distributed according to the PMR. A secant method is then used to find the unkown mass which is needed to match the planet radius. If a negative mass is needed an error is printed in the output file. An output file with mass of core, mantle, water, and atmosphere and the radius of the core, mantle, water, and planet and the target radius for each line in the input file and for each PMR. If the target radius is not matched within the error tolerance in 30 interations the finder moves to the next PMR or next line in the input file.

### mode5,6,7.cfg ###

These modes allow for changing an EoS in the model temporarily during a run (changing an EoS pernamently or adding an EoS discussed below). These modes can be used to measure how the uncertainty in a measurement affects a planet's radius. Mode 5 changes an EoS and uses the twolayer function. Mode 6 and Mode 7 iteratively change an EoS from an input file with the twolayer (Mode 2) and fullmodel (Mode 0) respectively.

For the moment, we save documenting how to create an input file for these modes for a later date. Examples of input files that change the EoSs are included in the input directory.

### mode8.cfg ###

This mode is our MCMC method which returns samples of interior solutions for a planet with a mass,mass uncertainty, radius, and radius uncertainty. These samples can be used to create trace and corner plots with `plot/mcmcsolutions.py`. This mode assumes gaussian priors and uniform priors on mass fractions from 0-1. This mode is under continued development; to add constraints please contact the Authors.

In the config file, the user sets the median and one sigma uncertainties for mass and radius. The user also designates the number of layers to solve---either 2 for core/mantle, 3 for core/mantle/hydrosphere, or 4 for core/mantle/hydrosphere/atmosphere. The user then defines the number of chains and the steps per chain to use for the planet. The surface temperature, temperature jumps, and phase diagrams are defined in the same way as mode 0. The above parameters will be used in the `mcmcsample()` function from `src/compfind.cpp`.

MAGRATHEA will use the full solver from mode 0 and step through the parameter space of mass and mass fractions for the given number of layers using a metropolis hastings algorithm. We suggest using at least 3 chains and 1000 steps per chain. The run time can be approximated as chain*steps seconds, but depends on the planet and the user's computer. The output is a tab separated file with the samples of mass and mass fraction of each layer with the calculated log-likelihood and radius of the planet.

## Build your own planet model ##

MAGRATHEA is built for model flexibility with transparent storage structures for equations of state (EoS). We make it easy for users to build a reproduceable model and cite the material measurements that have gone into the model.

The built-in EoSs for various planet building materials and phases are listed in the file `EOSlist.h`.  The detailed definition of each one can be found in `EOSlist.cpp`.


### Adding new equations of state ###

The code takes a new EoS in the format of (0) 3rd order Birch-Murnaghan, (1) 4th order Birch-Murnaghan, (2) Vinet, (3) Holzapfel, (4) Keane, (6) Ideal Gas, (7) Density-Pressure input table or Pressure-Temperature-Density-dT/dP_S input table from file for interpolation, (8-12) same as 0-4 in combination with RTPress. 

New phases should be listed into the file `EOSlist.h`, defined in `EOSlist.cpp`, and deleted in `main.cpp`.

Examples:

	// Epsilon Iron (hcp), Smith et al. 2018, Nature Astronomy. (Gruneisen determined from fitting Fig. 3b)
	// DEFAULT
	double Fe_hcp_array[][2] = {{0,2}, {1,mFe/8.43}, {2,177.7}, {3,5.64}, {5,mFe}, {7,322}, {8,2.09}, {9,1.01}, {10,0.0500}, {14,1}, {15,26}};
	EOS *Fe_hcp = new EOS("Fe hcp (Smith)", Fe_hcp_array, sizeof(Fe_hcp_array)/2/sizeof(Fe_hcp_array[0][0]));

	// Post-Perovskite, MgSiO3, Sakai, Dekura, & Hirao, 2016, Scientific Reports
	// DEFAULT
	double Si_PPv_Sakai_array[][2] = {{0,4}, {1,24.73}, {2,203}, {3,5.35}, {5,mMg+mSi+3*mO}, {7,848}, {8,1.47}, {9,2.7}, {10,0.93}, {14,5}};
	EOS *Si_PPv_Sakai = new EOS("Si PPv (Sakai)", Si_PPv_Sakai_array, sizeof(Si_PPv_Sakai_array)/2/sizeof(Si_PPv_Sakai_array[0][0]));

Analytical EoS are defined by an array with as many parameters as needed from the table below. For example, the Epsilon Iron EoS above has a bulk modulus of 177.7 GPa.

Index | Variable | Unit | Comment 
:---------: | :---------: | :----------: | -------------
0 | EOS formula type | | Index of the EoS formats indicated in parentheses [above](#adding-new-equations-of-state) 
1 |	V0 | cm<sup>3</sup> mol<sup>-1</sup> | Molar volume at reference point 
2 |	K0 | GPa | Bulk modulus 
3 |	K0' | | Pressure derivative of the bulk modulus. Default 4 
4 |	K0'' | GPa<sup>-1</sup> | Second pressure derivative 
5 |	m_mol | g mol<sup>-1</sup> | Molar mass 
6 |	P0 | GPa | The minimum pressure, corresponding to V0. Default 0
7 |	Theta0 | K | Fitting parameter of Einstein or Debye temperature. Default 1
8 |	gamma_0 | | Fitting parameter of Gruneisen parameter
9 | beta | | Fitting parameter of Gruneisen parameter
10 | gamma_inf | | Fitting parameter of Gruneisen parameter. Default 2/3
11 | gamma_0' | | Volume derivative of the Gruneisen parameter
12 | e0 | 10<sup>-6</sup> K<sup>-1</sup> | Electronic contribution to Helmholtz free energy. Default 0
13 | g | | Electronic analogue of the Gruneisen parameter
14 | n | | Number of atoms in the chemical formula. Default 1
15 |  Z | | Atomic number (number of electron) 
16 | T0 | K | Reference temperature for the thermal pressure. Default 300 
17 | alpha0 | 10<sup>-6</sup> K<sup>-1</sup> | The zeroth order coefficient of thermal expansion at a reference pressure P0 
18 | alpha1 | 10<sup>-6</sup> K<sup>-2</sup> | The first order coefficient of thermal expansion at a reference pressure P0 
19 | xi | | Power law index in the coefficient of thermal expansion. Default 0 
20 | c_p0 | 10<sup>7</sup> erg g<sup>-1</sup> K<sup>-1</sup> | Specific heat capacity at constant pressure
21 | c_p1 | 10<sup>7</sup> erg g<sup>-1</sup> K<sup>-2</sup> | Coefficient for specific heat capacity
22 | c_p2 | 10<sup>7</sup> erg g<sup>-1</sup> K  | Coefficient for specific heat capacity
23 | Debye approx |  | Positive number for Debye, otherwise Einstein 
24 |thermal type | | See Paper

Tabulated EoS can be read from files in the `tabulated` directory. A table can either have 2 or 4 columns. If 2, the user defines the density at each pressure and the code will assume the material is isothermal and interpolate the density from pressure alone. The table must monotonically increase in pressure. If 4 columns are provided, they must be pressure, temperature, density, and the adiabatic gradient (dT/dP_S). The code will linearly interpolate the density from both the temperature and pressure. The table must be rectangular in pressure and temperature with each temperature listed for a pressure before moving to the next pressure, and both must be increasing.

### Phase Diagrams ###

To pick desired materials/phases in each layer, change the corresponding return values of `find_phase_water_default`, `find_phase_Fe_default`,  `find_phase_Si_default`, or `find_phase_gas_default` in `phase.cpp` using conditionals to set the desired pressure and temperature where the material/phase will exist.

Example:

	// Fe Default: hcp and Liquid iron
	EOS* find_phase_Fe_default(double P, double T)
	{
 	 if (P <= 0 || T <= 0)
 	 {
  	  return NULL;
  	 }
  	 P /= 1E10;			// convert microbar to GPa
 	 // Default Core
 	 if( T > 12.8*P + 2424 && T > 13.7*P + 2328)   // melting curve from Dorogokupets et al. 2017, Scientific Reports. fcc and hcp Fe melting curve.
    	   return Fe_liquid;
 	 else
    	   return Fe_hcp;             // use hcp Iron for all regions.
	}

A number of predefined phase diagrams are available to test out which can be defined in the config files for the modes which will take changeable phase diagrams.

Layer | Name to call | Details
-------|:--------|-----------
Core | "Fe_default" | hcp and Liquid iron
&nbsp; | "Fe_fccbcc" | Includes low pressure fcc and bcc iron
Mantle |  "Si_default" | Upper Mantle: Fo, Wds, Rwd, and liquid ; Lower Mantle: Brg, PPv
&nbsp; | "Si_simple" | Brg, PPv, and liquid 
&nbsp; | "C_simple" | Graphite and Diamond, Carbon Mantle
&nbsp; | "SiC" | Silicon Carbide Mantle
&nbsp; | "PREM" | PREM tabulated mantle
Hydrosphere | "water_default" | H2O Water/Ice boundaries primarily form Dunaeva et al. 2010
&nbsp; | "water_tabulated" | AQUA Haldemann et al. 2020 Tabulated Ice, Liquid, Vapor, Supercritical
Atmosphere | "gas_default" | Ideal Gas: Isothermal for P<100 bar. Adiabatic ideal gas for P > 100 bar
&nbsp; | "HHe_tabulated" | H/He Gas: Isothermal ideal for P<100 bar, P>100 bar: tabulated real gas, Chabrier & Debras 2021 Y=0.275

### Adding new phase diagrams ###

Saving a phase diagram for repeated use and comparison with others is often helpful. There is a number of places where a new phase diagram needs to be added. First, define the new `find_phase_X_Y` in `phase.h` then create the function using conidtionals which return an EoS in `phase.cpp` folowing the example of other find_phase functions. Then later in `phase.cpp` pass the function into a phase diagram (`PhaseDgm`) and the layer to which it belongs. Give the `PhaseDgm` a new name (we suggest to continue to interate on our format) and add this new PhaseDgm near the end of `phase.h`. Lastly in `main.cpp` the parser must find the user defined string in the config file and pass it to the vector of phase diagrams in the `//Set Phase Diagrams` section. Write a new conditional for the correct layer and set the correct position in the vector to your new `PhaseDgm`.


### Useful unit conversions ###

The input parameters required to create an EOS in the program are typically in cgs units. However, the EoS measurement/calculation article may list their result in different units. Thus, the following unit conversions maybe helpful.

In the equations below, *m* is the number of formula units per unit cell. For example, cubic ice-VII (space group Pn3m) has two water molecules per unit cell, so *m=2*.
*n* is the number atoms per molecule formula. For example, MgSiO<sub>3</sub> has 5 atoms in a molecule, so *n=5*.

* 1 &#8491;<sup>3</sup>/cell = 10<sup>-24</sup>N<sub>A</sub>/m cm<sup>3</sup>/mol = 0.6022/m cm<sup>3</sup>/mol.
* 1 &#8491;<sup>3</sup>/atom = 10<sup>-24</sup>nN<sub>A</sub> cm<sup>3</sup>/mol = 0.6022n cm<sup>3</sup>/mol.
* 1 eV/atom = 1.602&times;10<sup>-12</sup>nN<sub>A</sub> erg/mol = 9.649&times;10<sup>11</sup>n erg/mol. 
* 1 GPa = 10<sup>10</sup>&micro;bar = 0.01 Mbar.

## Functions to obtain the calculated planetary parameters ##
The functions in the table below can be used to obtain the calculated planetary parameters after a successful solving of planetary structure.


Function | Output unit | Comment
-------|:--------:|-----------
`double totalM()` | g | return the total mass of a planet
`double totalR()` | cm | return the total radius of a planet
`int getLayer_from_r(double r)` |  | Input radius in the RE, return the layer index `l` (from 0, count from bottom), rb(l)<=r*R⊕<rb(l+1)
`int getLayer_from_m(double m)` |  | Input mass in the ME, return the layer index `l` (from 0, count from bottom), M(l)<=m*M⊕<M(l+1)
`double getP(int l)` | &micro;bar | return the pressure at layer `l`
`double getM(int l)` | g | return the pressure at layer `l`
`double getR(int l)` | cm | return the radius at layer `l`
`double getrho(int l)` | g/ cm<sup>3</sup> | return the ensity at layer `l`
`double getT(int l)` | K | return the temperature at layer `l`
`int getsize()` |  | return the total number of layers
`vector<double> getRs()` | R⊕ | Return the radii of core, mantle, and water layer, and the total radius in the unit of earth radii.
`vector<double> getTs()` | K | Return the temperatures at the outer side of each component interfaces as well as planet surface 

Example code:

    vector<double> Tgap = {0, 0, 0, 300};
    vector<double> Mcomp =  {1.0,0.5,0.1,0.00001}; 
    planet=fitting_method(Comp, Mcomp, Tgap, ave_rho, P_surface, false);
    if (planet)
    {
      int l = planet->getLayer_from_r(1);
      cout<<"layer:"<<l<<" P="<<planet->getP(l)<<"microbar R="<<planet->getR(l)/RE<<"REarth  rho="<<planet->getrho(l)<<"g/cm^3"<<endl;
      l = planet->getLayer_from_m(1);
      cout<<"layer:"<<l<<" P="<<planet->getP(l)<<"microbar M="<<planet->getM(l)/ME<<"MEarth  rho="<<planet->getrho(l)<<"g/cm^3"<<endl;
      vector<double> Rs = planet->getRs();
      vector<double> Ts = planet->getTs();
      for (int i=0; i<4; i++)
        cout<<planet->getLayer_from_r(Rs[i])<<' '<<Rs[i]<<"REarth "<<Ts[i]<<'K'<<endl;
    }

Example output:

	layer:605 P=1.21568e+11microbar R=0.998631REarth  rho=1.95656g/cm^3
	layer:547 P=1.5147e+12microbar M=0.999971MEarth  rho=11.8156g/cm^3
	547 0.726926REarth 630.97K
	602 0.978419REarth 460.504K
	682 1.07939REarth 300K
	752 1.12843REarth 300K

## Print EoS into a table ##

The code can print the pressure-density relation of a built-in EoS into an ASCII table by using the `printEOS` function of the EoS object. This can be a simple way to double check the EoS is set up correctly. The table covers the pressure range from 0.1 GPa (or P<sub>0</sub> if it is larger) to 2000 GPa at the temperature T<sub>0</sub> (default 300 K) of the EoS. The output table file is located at `./tabulated/phasename.txt`.


## Examples/Plots ##

Within the directories `input` and `result` are example input files and output files which are described in each directory's markdown (.md) file.

Python plotting scripts are included in the directory `plot`. Scripts are written for Python 3.6. Scripts may require matplotlib, numpy, astropy, python-ternary, regex, cmasher, and scipy.

### 3D Representations of Planets ###

Magrathea outputs can be imaged in the 3D open-source software Blender. Instructions and code are in https://github.com/DavidRRice/Blender-Magrathea.
 <p align="center">
<img width = "200" src="plot/planet1.png"/>
 </p>

## Don't Panic (FAQ) ##

_**Who do I talk to?**_

Find our emails on our websites:

Chenliang Huang, Shanghai Astronomical Observatory [website](https://huang-cl.github.io/)

David R. Rice, ARCO, Open University of Israel [website](https://davidrrice.github.io/)

_**Where is the EoS/functionality I want?**_

Open an issue with details of what you need. 

_**What if I need to model a planet, but I'm currently being tortured by Vogon poetry?**_

We are open to collaborations! Emails found on websites above.

_**Where do I find other codes for exoplanet modeling?**_

MAGRATHEA and many other useful exoplanet-related codes are archived on NASA's [Exoplanet Modeling and Analysis Center](https://emac.gsfc.nasa.gov/).

_**Where does the name MAGRATHEA come from?**_

Magrathea is a fictional planet in Douglas Adams's *The Hitchiker's Guide to the Galaxy*: 

"for all the richest and most successful merchants life inevitably became rather dull and niggly, and they began to imagine that this was therefore the fault of the worlds they'd settled on... And thus were created the conditions for a staggering new form of specialist industry: custom-made luxury planet building. The home of this industry was the planet Magrathea, where hyperspatial engineers sucked matter through white holes in space to form it into dream planets- gold planets, platinum planets, soft rubber planets with lots of earthquakes- all lovingly made to meet the exacting standards that the Galaxy's richest men naturally came to expect.

But so successful was this venture that Magrathea itself soon became the richest planet of all time and the rest of the Galaxy was reduced to abject poverty. And so the system broke down, the Empire collapsed, and a long sullen silence settled over a billion worlds...

Magrathea itself disappeared and its memory soon passed into the obscurity of legend. In these enlightened days of course, no one believes a word of it."

---

*MAGRATHEA is or has been supported by the [Nevada Center for Astrophysics](https://www.physics.unlv.edu/~bzhang/NCfA.html), the University of Nevada, Las Vegas's [Physics & Astronomy Department](https://www.physics.unlv.edu/) and [Star & Planet Formation Group](https://unlv-spfg.github.io/), the [Astrophysics Research Center](https://arco.org.il/) of the Open University of Israel, and the University of Arizona's [Lunar and Planetary Laboratory](https://www.lpl.arizona.edu/)*

*We want to thank the many users and supporters of this code. We greatly appreciate the contributions, discussions, and support of [Allona Vazan](https://www.openu.ac.il/en/personalsites/AllonaVazan.aspx), [Michael Lozovsky](https://michloz8.wixsite.com/michael-lozovsky), [Ashkan Salamat](https://nexcl.unlv.edu/the-team/ashkan-salamat), and [Piyush Puranik](https://github.com/preppie22)*
