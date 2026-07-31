# Installation and First Run

Magrathea is an open-source **1-D planetary interior structure solver** written in C++.  
It is designed for astronomers, planetary scientists, and students who want to model the internal structure of differentiated planets.

\---

## 🪐 Overview

Magrathea integrates the hydrostatic equilibrium equations to compute the radius, pressure, temperature, density, and phase structure of a planet consisting of up to **four layers**:

* **Core** – iron or alloy
* **Mantle** – silicate
* **Hydrosphere** – water/ice
* **Atmosphere** – ideal gas

Each layer uses a modular **phase diagram** and **equation of state (EOS)** that can be modified or replaced.  
The solver outputs the full radial structure and layer boundaries, allowing you to explore how composition or temperature affect planetary size and density.

\---

## ⚙️ 1. Prerequisites

Magrathea requires the **GNU Scientific Library (GSL ≥ 2.0)**.

### Linux / macOS

```bash
sudo apt install libgsl27 libgsl-dev gsl-bin     # Ubuntu / Debian
# or
brew install gsl                                 # macOS
```

To build from source:

```bash
tar -xzf gsl-X.Y.tar.gz
cd gsl-X.Y
./configure
make
sudo make install
```

### Windows

Use **Windows Subsystem for Linux (WSL)** and follow the same Linux instructions inside the terminal.  
(Direct compilation on native Windows is not supported because of path issues.)

\---

## 📦 2. Getting the Code

Clone from GitHub:

```bash
git clone https://github.com/Huang-CL/Magrathea.git
cd Magrathea
```

If GSL is not installed in a standard path (e.g. `/usr/local/`), edit the **Makefile**:

```makefile
CFLAGS += -I/path/to/gsl/include
LDFLAGS += -L/path/to/gsl/lib
```

To find the correct paths:

```bash
gsl-config --cflags
gsl-config --libs
```

\---

## 🧱 3. Building Magrathea

From the top-level Magrathea directory, run:

```bash
make -B
```

This creates an executable named **`planet`** in the project directory. After the first compilation, use `make` when C++ source files have changed.

Run the automated tests to check the installation and core functionality:

```bash
./planet --test
```

\---

## 🚀 4. Running Your First Planet

The default `run/mode0.cfg` provides an Earth-like example with a core and mantle and no hydrosphere or atmosphere.

1. From the top-level Magrathea directory, run:

```bash
./planet run/mode0.cfg
```

2. Inspect the end of the output file to find the calculated planet radius:

```bash
tail -2 result/StructureEarth.txt
```

To create a different planet, open `run/mode0.cfg` in a text editor and adjust these parameters:

|Parameter|Description|
|-|-|
|`mass\_of\_core`, `mass\_of\_mantle`, `mass\_of\_hydro`, `mass\_of\_atm`|Layer masses (Earth masses)|
|`surface\_temp`, `temp\_jump\_1–3`|Surface temperature and discontinuities between layers|
|`output\_file`|Output file name and path|

Save the file and rerun, from the top-level Magrathea directory:

```bash
./planet run/mode0.cfg
```

Changes to a `.cfg` file do not require recompilation. Run `make` only after changing C++ source files.

If you encounter an error like:

```
error while loading shared libraries: libgsl.so.23: cannot open shared object file
```

add the GSL library path:

```bash
export LD\_LIBRARY\_PATH=/usr/local/lib
source \~/.bashrc
```

\---

## 📊 5. Output

Magrathea writes an ASCII table (by default in `result/`) containing:

```
Pressure (GPa) | Enclosed Mass (M⊕) | Density (g cm⁻³) | Temperature (K) | Phase
```

and prints the total planet radius and radii of each compositional boundary.

\---

## 🧩 6. Next Steps: Other Run Modes

Magrathea includes **nine configurable modes** defined by the `.cfg` files in the `run/` directory:

|Mode|Purpose|
|-|-|
|**0**|Full 4-layer solver (default)|
|**1**|Isothermal, temperature-free solver|
|**2**|Two-layer mass–radius curves|
|**3**|Bulk solver for many planets|
|**4**|Composition finder (match observed M \& R)|
|**5–7**|EOS uncertainty \& iteration modes|
|**8**|MCMC composition inference (experimental)|

Try the [Practice Problems](user-guide/practice-problems.md) for guided examples, or download the [full tutorial handout](https://github.com/Huang-CL/Magrathea/blob/master/docs/Tutorial_Practice_Problems.pdf).

\---

## 📚 7. Learn More

* **Publication:** Huang et al. (2022), *MAGRATHEA: An open-source spherical symmetric planet interior structure code*, *MNRAS*
* **GitHub:** [Huang-CL/Magrathea](https://github.com/Huang-CL/Magrathea)
* **Citation info:** see [`CITATION.md`](https://github.com/Huang-CL/Magrathea/blob/paper/CITATION.md)

\---

*Developed by Chenliang Huang, David R. Rice, and Jason H. Steffen at the University of Nevada Las Vegas.*