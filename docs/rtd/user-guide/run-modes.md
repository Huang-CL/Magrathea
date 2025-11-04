# Run Modes and .cfg files

MAGRATHEA uses human-readable `.cfg` configuration files to control each planetary interior calculation found in `\run`.  
Each file specifies global solver parameters, numerical tolerances, and mode-specific inputs.

---

## Global Configuration Parameters

These options are shared across **all run modes** and should appear near the top of every `.cfg` file.

| Parameter | Type | Description |
|------------|------|--------------|
| `rho_eps_rel` | float | Relative error tolerance of the density solver (default `1e-11`) |
| `T_eps_rel` | float | Relative error tolerance of temperature solver for adiabatic profiles (`1e-11`) |
| `ode_eps_rel0` | float | ODE integrator tolerance (first pass) |
| `ode_eps_rel1` | float | ODE integrator tolerance (refinement) |
| `R_eps_rel` | float | Radius fitting tolerance (typ. `2e-5`) |
| `ode_eps_rel2` | float | ODE integrator tolerance for inside-out mode |
| `P_eps_rel` | float | Pressure fitting tolerance |
| `fit_eps_rel` | float | Fitting error tolerance for second-round iteration |
| `verbose` | bool | Print debug and warning messages |
| `P_surface` | float | Surface pressure (µbar) at which radius is defined |
| `ave_rho_core` | float | Typical core density (g/cm³, e.g. 15) |
| `ave_rho_mantle` | float | Typical mantle density (e.g. 5) |
| `ave_rho_hydro` | float | Typical water/ice density (e.g. 2) |
| `ave_rho_atm` | float | Typical atmospheric density (e.g. 1e-3) |

Example global block:

```cfg
rho_eps_rel = 1e-11
T_eps_rel = 1e-11
ode_eps_rel0 = 1e-7
ode_eps_rel1 = 1e-10
R_eps_rel = 2e-5
ode_eps_rel2 = 1e-10
P_eps_rel = 1e-10
fit_eps_rel = 1e-4
verbose = false
P_surface = 1e5
ave_rho_core = 15
ave_rho_mantle = 5
ave_rho_hydro = 2
ave_rho_atm = 0.001
```

---

## Overview

Each run mode is selected via:

```cfg
input_mode = [0–8]
```

Run modes differ in input requirements and outputs but share the same solver infrastructure from `hydro.cpp` and `phase.cpp`.

---

## Mode 0 — Full 4-Layer Solver (Default)

**Purpose:** Compute the full interior structure of a differentiated planet using self-consistent temperature and density profiles.

**Description:**  
Integrates the hydrostatic and thermal equations inward from the surface.  
Supports four compositional layers: core, mantle, hydrosphere, and atmosphere.

**Key configuration parameters:**

| Parameter | Description |
|------------|--------------|
| `mass_of_core` | Core mass in Earth masses |
| `mass_of_mantle` | Mantle mass |
| `mass_of_hydro` | Water/ice layer mass |
| `mass_of_atm` | Atmosphere mass |
| `surface_temp` | Surface temperature (K) |
| `temp_jump_1`–`temp_jump_3` | Temperature discontinuities between layers |
| `output_file` | Path for results file |
| `core_phasedgm`, `mantle_phasedgm`, `hydro_phasedgm`, `atm_phasedgm` | Phase diagram presets (e.g., `Fe_default`, `Si_default`, `water_default`, `gas_default`) |

**Example:**

```cfg
input_mode = 0
mass_of_core = 0.33
mass_of_mantle = 0.67
mass_of_hydro = 0.0
mass_of_atm = 0.0
surface_temp = 300
temp_jump_1 = 0
temp_jump_2 = 0
temp_jump_3 = 0
core_phasedgm = Fe_default
mantle_phasedgm = Si_default
hydro_phasedgm = water_default
atm_phasedgm = gas_default
output_file = result/earthlike.txt
```

**Output:** Pressure–density–temperature–phase profiles and total radius.

---

## Mode 1 — Temperature-Free Solver

**Purpose:** Compute structure assuming an **isothermal** profile (no thermal gradient).

**Parameters:** Same as Mode 0, except no temperature gradients are used.  
`isothermal = true` is implied.

**Output:** Planet radius and layer boundaries.

---

## Mode 2 — Two-Layer Mass–Radius Solver

**Purpose:** Generate a **mass–radius curve** for two-layer planets with fixed mass fractions.

| Parameter | Description |
|------------|--------------|
| `layer_index` | Which layer to exclude: 0=no core, 1=no mantle, 2=no water |
| `mass_frac` | Mass fraction of the inner layer |
| `min_mass`, `max_mass` | Range of planet masses (Earth masses) |
| `step_mass` | Step size (Earth masses) |
| `output_file` | Output file name |

**Example:**

```cfg
input_mode = 2
layer_index = 0
mass_frac = 0.33
min_mass = 0.1
max_mass = 10
step_mass = 0.1
output_file = result/twolayer_curve.txt
```

**Output:** Tabulated mass–radius pairs.

---

## Mode 3 — Bulk Planet List Solver

**Purpose:** Compute radii and layer boundaries for multiple planets from an input table.

| Parameter | Description |
|------------|--------------|
| `input_file` | Path to table with mass and layer fractions |
| `output_file` | Output file for all modeled planets |
| `solver` | Solver type (integer, typically `1` for default) |

**Example input table (`planets.txt`):**
```
Mass  f_core  f_mantle  f_water
1.0   0.33    0.67      0.0
5.0   0.25    0.70      0.05
```

**Configuration example:**
```cfg
input_mode = 3
input_file = planets.txt
solver = 1
output_file = result/bulk_output.txt
```

**Output:** Table of modeled layer radii for each planet.

---

## Mode 4 — Composition Finder

**Purpose:** Given observed **mass and radius**, infer the most likely internal composition.  
Searches layer fractions to match the measured radius within uncertainty.

| Parameter | Description |
|------------|--------------|
| `input_file` | Table of mass–radius samples |
| `find_layer` | Index of layer to solve for (1–3) |
| `layer_inner`, `layer_outer` | Which layers to keep constant |
| `PMR_min`, `PMR_max`, `PMR_step` | Range and step for partial mass ratios |
| `R_error` | Allowed fractional radius error |
| `output_file` | Results file |

**Example:**

```cfg
input_mode = 4
input_file = observed_samples.txt
find_layer = 3
layer_inner = 1
layer_outer = 2
PMR_min = 0.1
PMR_max = 0.9
PMR_step = 0.05
R_error = 0.02
output_file = result/composition_search.txt
```

**Output:** Mass fractions of each layer for observed M–R pairs.

---

## Mode 5 — On-the-Fly EOS Modification

**Purpose:** Change an EOS definition directly from an input file to test sensitivity.  
MAGRATHEA reads a custom EOS file and recomputes the planetary structure.

| Parameter | Description |
|------------|--------------|
| `eos_file` | Path to custom EOS definition |
| `output_file` | Output file name |

**Example:**

```cfg
input_mode = 5
eos_file = input/custom_water_eos.txt
output_file = result/customEOS_test.txt
```

---

## Mode 6 — Iterated EOS Variation (Two-Layer)

**Purpose:** Iterate EOS parameter modifications through the **two-layer solver (Mode 2)**.

| Parameter | Description |
|------------|--------------|
| `input_file` | File containing list of EOS variations |
| `output_file` | Output results file |
| `layer_index` | Which layer is excluded |
| `mass_frac` | Mass fraction of the inner layer |

**Example:**

```cfg
input_mode = 6
input_file = eos_variations.txt
output_file = result/mr_variations.txt
layer_index = 2
mass_frac = 0.5
```

**Output:** Radius curves for each EOS variant.

---

## Mode 7 — Iterated EOS Variation (Full Solver)

**Purpose:** Apply EOS modifications across **full multi-layer** planets (Mode 0 equivalent).

| Parameter | Description |
|------------|--------------|
| `input_file` | EOS variation table |
| `output_file` | Output results file |

**Example:**

```cfg
input_mode = 7
input_file = eos_table.txt
output_file = result/fullmodel_variations.txt
```

**Output:** Profiles or radii for each EOS configuration.

---

## Mode 8 — MCMC Composition Inference

**Purpose:** Retrieve posterior distributions of core, mantle, and water fractions consistent with observed **mass–radius uncertainties**.

| Parameter | Description |
|------------|--------------|
| `mass_prior`, `mass_unc` | Observed planet mass and uncertainty |
| `radius_prior`, `radius_unc` | Observed radius and uncertainty |
| `num_layers` | Number of layers in model (2–4) |
| `num_chains` | Number of parallel MCMC chains |
| `chain_steps` | Number of steps per chain |
| `output_file` | Output MCMC sample file |

**Example:**

```cfg
input_mode = 8
mass_prior = 1.00
mass_unc = 0.05
radius_prior = 1.10
radius_unc = 0.03
num_layers = 3
num_chains = 4
chain_steps = 5000
output_file = result/mcmc_samples.txt
```

**Output:** MCMC chain file with `log_likelihood`, `fCore`, `fMantle`, `fWater`, and `RPlanet` columns.

---

## Summary Table

| Mode | Description | Typical Output |
|------|--------------|----------------|
| 0 | Full 4-layer self-consistent solver | Profiles (`.txt`) |
| 1 | Isothermal solver | Radii only |
| 2 | Two-layer solver | Mass–radius curves |
| 3 | Bulk solver | Many planets |
| 4 | Composition finder | Layer fractions |
| 5 | EOS modification | Radius impact |
| 6 | Iterated EOS (2-layer) | Curves per EOS set |
| 7 | Iterated EOS (full) | Profiles per EOS set |
| 8 | MCMC inference | Posterior samples |

---

_See also: [Quick Installation and First Run](quickstart.md)_
