# Practice Problems

These short projects introduce MAGRATHEA's forward solver, bulk mode, equation-of-state controls, and composition finder. Complete the [Installation and First Run](../getting-started.md) guide before beginning.

The exercises below are shortened from the full [MAGRATHEA Tutorial and Practice Problems handout](https://github.com/Huang-CL/Magrathea/blob/master/docs/Tutorial_Practice_Problems.pdf), which includes additional background, teaching notes, and extensions.

:::{note}
Run all commands from the MAGRATHEA repository root. Use a different `output\_file` for each calculation because some modes append to an existing file.
:::

## Project 1: Match a Planet's Mass and Radius

**Goal:** Explore how the core-to-mantle ratio changes a rocky planet's radius.

Use `run/mode0.cfg` to model either:

|Planet|Mass|Target radius|Reference|
|-|-:|-:|-|
|LHS 1140 b|6.38 M⊕|1.635 R⊕|[Lillo-Box et al. (2020)](https://arxiv.org/abs/2010.06928)|
|HD 137496 b|4.04 M⊕|1.31 R⊕|[Azevedo Silva et al. (2022)](https://arxiv.org/abs/2111.08764)|

The LHS 1140 b values follow the original practice handout. Updated measurements may differ.

1. Set the total core and mantle mass equal to the planet mass. Set the hydrosphere and atmosphere masses to zero.
2. Begin with a core mass fraction of about 0.32. For LHS 1140 b:

```ini
   mass\_of\_core=2.0416
   mass\_of\_mantle=4.3384
   mass\_of\_hydro=0.0
   mass\_of\_atm=0.0
   ```

3. Set a unique output filename, then compile and run:

```bash
   make
   ./planet run/mode0.cfg
   ```

4. Read the planet radius from the end of the output file.
5. Change the core and mantle masses while keeping the total mass fixed. Repeat until the model radius approaches the target radius.
6. Report the core mass fraction:

   \[
\\mathrm{CMF} = \\frac{M\_{\\mathrm{core}}}{M\_{\\mathrm{planet}}}.
]

**Extension:** Add a hydrosphere or atmosphere and find additional structures that match the same mass and radius.

\---

## Project 2: Build a Mass--Radius Diagram

**Goal:** Calculate mass--radius curves for several fixed compositions.

Use `run/mode3.cfg` to model planets from 0.1 to 4 M⊕ with these compositions:

|Model|Core fraction|Mantle fraction|Water fraction|
|-|-:|-:|-:|
|Pure iron|1.000|0.000|0.000|
|Mercury-like|0.680|0.320|0.000|
|Equal core and mantle|0.500|0.500|0.000|
|Earth-like refractory ratio|0.325|0.675|0.000|
|Pure mantle|0.000|1.000|0.000|
|Pure water|0.000|0.000|1.000|

1. Create one bulk input file for each composition. Each row should contain:

```text
   Mass  fCore  fMantle  fWater
   ```

   Example files are available in `input/`; `input/massradiusinput.py` can also be adapted to generate a mass grid.

2. In `run/mode3.cfg`, set:

```ini
   solver=1
   surface\_temp=300
   temp\_jump\_1=0
   temp\_jump\_2=0
   temp\_jump\_3=0
   ```

3. Change `input\_file` and `output\_file` for each composition, then run:

```bash
   ./planet run/mode3.cfg
   ```

4. Plot planet radius against total mass for all six models.

**Extension:** Add measured planets or repeat the calculation with a thin atmosphere.

\---

## Project 3: Change a Planet Model

**Goal:** Test how material and EOS choices affect a water-rich planet.

Begin with a 1 M⊕ planet in `run/mode0.cfg`:

```ini
mass\_of\_core=0.165
mass\_of\_mantle=0.335
mass\_of\_hydro=0.500
mass\_of\_atm=0.0
surface\_temp=300
```

Run the default model and save its output. Then test the following variants separately.

### A. Tabulated Water EOS

Set:

```ini
surface\_temp=350
hydro\_phasedgm="water\_tabulated"
```

This selects the AQUA tabulated water EOS. Use a new output filename and rerun the model.

### B. Iron--Silicate Core

The EOS library includes `Fe\_7Si` and `Fe\_15Si` in `src/EOSlist.cpp`. In `src/phase.cpp`, temporarily replace the solid `Fe\_hcp` return value in `find\_phase\_Fe\_default()` with one of these alloys.

Recompile and rerun:

```bash
make
./planet run/mode0.cfg
```

### C. Post-Perovskite Bulk Modulus

`Si\_PPv\_Sakai` is defined in `src/EOSlist.cpp` with a bulk modulus of 203 GPa. Temporarily change its index-2 parameter to 250 GPa, recompile, and rerun.

:::{warning}
Make source-code experiments on a temporary branch or save the original files. Restore the default EOS values after completing the project.
:::

Compare the default, tabulated-water, iron--silicate-core, and modified post-perovskite models. Record:

* Planet radius
* Core--mantle boundary radius
* Pressure and temperature at the core--mantle boundary

**Extension:** Use `plot/quickdensityplot.py` to compare the interior profiles. Consider whether the large post-perovskite change produces physically unrealistic behavior.

See [Phase Diagrams](../models/phasediagrams.md) and the [EOS Library](../models/eos.md) for additional model options.

\---

## Project 4: Use the Composition Finder

**Goal:** Use `run/mode4.cfg` to find interior structures matching the LHS 1140 b values used in the original handout.

Create `input/LHS1140binput.txt`:

```text
M (Earth-masses)  R (Earth-radii)
6.38               1.635
```

Set the input and output paths in `run/mode4.cfg`.

### A. Dry Core--Mantle Solution

Find the mantle mass while allowing no water:

```ini
find\_layer=2
layer\_inner=1
layer\_outer=3
PMR\_min=0.0
PMR\_max=0.0
```

Run:

```bash
./planet run/mode4.cfg
```

Use the output to calculate the dry planet's core mass fraction. Compare it with Project 1.

### B. Maximum Water for an Iron--Water Planet

Find the hydrosphere mass while forcing the mantle mass to zero:

```ini
find\_layer=3
layer\_inner=1
layer\_outer=2
PMR\_min=0.0
PMR\_max=0.0
```

This gives an iron--water end-member solution.

### C. Water Fraction Across Core--Mantle Ratios

Keep the core and mantle in a range of partial mass ratios while solving for water:

```ini
find\_layer=3
layer\_inner=1
layer\_outer=2
PMR\_min=20.0
PMR\_max=<upper mantle PMR>
PMR\_step=1.0
```

Here,

\[
\\mathrm{PMR} = 100
\\frac{M\_{\\mathrm{outer}}}
{M\_{\\mathrm{inner}} + M\_{\\mathrm{outer}}}.
]

Run the model and plot water mass fraction against core mass fraction.

**Extension:** Add draws from the observational mass and radius uncertainties to the input file and repeat the calculation.

See [Run Modes and Configuration Files](run-modes.md) for a complete description of mode 4.

