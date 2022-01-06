# Magrathea #

(Excerpt from The Hitchhiker's Guide to the Galaxy, Page 634784, Section 5a, Entry: MAGRATHEA)

Planet interior structure code for astronomers, planetary scientists, mice, and more by Huang, Rice, and Steffen (2022).


## What is this repository for? ##

A planet structure code which considers the case of fully differentiated interiors.
The code integrates the hydrostatic equation in order to shoot for the correct planet radius givent the mass in each layer.
The code returns the pressure, temperature, density, phase, and radius at steps of enclosed mass.
The code support 4 layers: core, mantle, hydrosphere, and atmosphere. Each layer has a phase diagram with equations of state chosen for each phase.  
The code was developed by Chenliang Huang and David R. Rice at UNLV starting in 2017.


## How do I get set up? ##

To run the code, one needs to [install the gsl library](https://www.gnu.org/software/gsl/)(>= v2.0).  Download the compressed package from the [GNU ftp site](ftp://ftp.gnu.org/gnu/gsl/).  Extract the file and install the package following the instruction in the `INSTALL` file.  For the simplest case, 
 
    sudo ./configure
    sudo make
    sudo make install

should configure, build, and install the gsl package.  A few prerequisites, such as `g++`, `make`, `make-guile`, may need to be install following the error messages throughout the process. 

On Ubuntu system, the gsl package can also be installed from Ubuntu repository using 
`sudo apt-get install libgsl23 libgsl-dev gsl-bin`.

If the gsl library is not installed globally (under /usr/local/ or equivalent), the actual path toward the gsl headers (e.g. `~/include` or `~/gsl/include`) and gsl library (e.g. `~/lib` or `~/gsl/lib`) should be put in the `Makefile` following `-I` in `CFLAGS` and `-L` in `LDFLAGS`.  Note: The path after `-I` should end with `/include`, not `/include/gsl`.  The path of headers and libraries can be find with `gsl-config --cflags` and ` gsl-config --libs`.

On Windows system, download the [cyqwin terminal](https://www.cygwin.com/). Include package 'gsl' upon installation. We suggest including all packages in devel, science, math, and python along with a text editor. Directory will be found in: /cygdrive/c/Users/USERNAME/PATH_TO_MAGRATHEA.

For the first time, or working with a new computer, compile the code with `make -B`, or `make clean` then `make`.
Later on, compile the code with `make`.
The name of the compiled program is called planet by default.  To run the compiled program, run
`./planet`.

If an error message like "error while loading shared libraries: libgsl.so.23: cannot open shared object file: No such file or directory" is reported when running the code, add `export LD_LIBRARY_PATH=/usr/local/lib` (directory of gsl library files) to the `.bashrc` file, or add `setenv LD_LIBRARY_PATH /usr/local/lib` to the `.cshrc` file.

## Configuration ##

Magrathea.pdf (Huang, Rice, & Steffen 2022) contains the most compelete details for using the first version of the code.

### Pick phases ###

The build-in phases are listed in the a file `EOSlist.h`.  The detailed definition of each one can be found in `EOSlist.cpp`.

To pick desired phases, change the corresponding return values of `find_water_phase`, `find_Fe_phase`, or `find_Si_phase` in `phase.cpp` using conditionals to set the desired pressure and temperature according to the phase diagram.

### Adding new phases ###

The code takes new EOS in the format of (1) 3rd order Birch-Murnaghan, (2) 4th order Birch-Murnaghan, (3) Vinet, (4) Holzapfel, (5) Keane, (6) input table from file for interpolation, (7) user-defined function.
New phases should be listed into the file `EOSlist.h` and defined in `EOSlist.cpp`.


## Capability and Output ##

### Construct a planet ###
The basic capability. Calculate the structure of a planet given the mass of each layer.
Adjust the parameters of getmass function in `main.cpp`.  Four input parameters are the planet core mass, mantle mass, hydrosphere mass, and atmosphere mass in the unit of Earth mass. Temperature of the surface and discontinuities in temperature at layer boundaries can also be inputted.
The planet structure will be output as an ASCII file with pressure (GPa), interior mass (Earth mass), density (g cm$^{-3}$), temperature (K), and phase of composition as a function of radius (Earth radius).


## FAQ ##
Don't Panic.


### Who do I talk to? ###

Chenliang Huang, University of Arizona [website](https://www.lpl.arizona.edu/~huangcl/)

David R. Rice, University of Nevada, Las Vegas [website](https://www.physics.unlv.edu/~drice986/)
