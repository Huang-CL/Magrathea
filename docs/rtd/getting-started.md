# MAGRATHEA: Quick Installation and First Run

MAGRATHEA is an open-source **1-D planetary interior structure solver** written in C++.  
It is designed for astronomers, planetary scientists, and students who want to model the internal structure of differentiated planets.

---

## 🪐 Overview

MAGRATHEA integrates the hydrostatic equilibrium equations to compute the radius, pressure, temperature, density, and phase structure of a planet consisting of up to **four layers**:

- **Core** – iron or alloy  
- **Mantle** – silicate  
- **Hydrosphere** – water/ice  
- **Atmosphere** – ideal gas  

Each layer uses a modular **phase diagram** and **equation of state (EOS)** that can be modified or replaced.  
The solver outputs the full radial structure and layer boundaries, allowing you to explore how composition or temperature affect planetary size and density.

---

## ⚙️ 1. Prerequisites

MAGRATHEA requires the **GNU Scientific Library (GSL ≥ 2.0)**.

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

---

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

---

## 🧱 3. Building MAGRATHEA

Run:

```bash
make -B        # first build
make           # subsequent builds
```

This creates an executable named **`planet`** in the project directory.

---

## 🚀 4. Running a Test Planet

The simplest run mode is `mode0.cfg`, which solves for a planet’s radius and interior given its layer masses and temperatures.

1. Open `run/mode0.cfg` in a text editor.  
2. Adjust these key parameters:

| Lines | Parameter | Description |
|-------|------------|-------------|
| 12–14 | `mass_of_core`, `mass_of_mantle`, `mass_of_hydro`, `mass_of_atm` | Layer masses (Earth masses) |
| 16–20 | `surface_temp`, `temp_jump_1–3` | Surface temperature and discontinuities between layers |
| 21 | `output_file` | Output file name and path |

**Example**

```cfg
mass_of_core = 0.33
mass_of_mantle = 0.67
mass_of_hydro = 0.0
mass_of_atm = 0.0

surface_temp = 300
temp_jump_1 = 0
temp_jump_2 = 0
temp_jump_3 = 0

output_file = result/earthlike.txt
```

3. Compile and run:

```bash
make
./planet run/mode0.cfg
```

If you encounter an error like:

```
error while loading shared libraries: libgsl.so.23: cannot open shared object file
```

add the GSL library path:

```bash
export LD_LIBRARY_PATH=/usr/local/lib
source ~/.bashrc
```

---

## 📊 5. Output

MAGRATHEA writes an ASCII table (by default in `result/`) containing:

```
Pressure (GPa) | Enclosed Mass (M⊕) | Density (g cm⁻³) | Temperature (K) | Phase
```

and prints the total planet radius and radii of each compositional boundary.

---

## 🧩 6. Next Steps: Other Run Modes

MAGRATHEA includes **nine configurable modes** defined by the `.cfg` files in the `run/` directory:

| Mode | Purpose |
|------|----------|
| **0** | Full 4-layer solver (default) |
| **1** | Isothermal, temperature-free solver |
| **2** | Two-layer mass–radius curves |
| **3** | Bulk solver for many planets |
| **4** | Composition finder (match observed M & R) |
| **5–7** | EOS uncertainty & iteration modes |
| **8** | MCMC composition inference (experimental) |

For examples, open **`Tutorial_Practice_Problems.pdf`** in the repository.

---

## 📚 7. Learn More

- **Publication:** Huang et al. (2022), *MAGRATHEA: An open-source spherical symmetric planet interior structure code*, *MNRAS*  
- **GitHub:** [Huang-CL/Magrathea](https://github.com/Huang-CL/Magrathea)
- **Citation info:** see [`CITATION.md`](https://github.com/Huang-CL/Magrathea/blob/paper/CITATION.md)

---

_Developed by Chenliang Huang, David R. Rice, and Jason H. Steffen at the University of Nevada Las Vegas._
