# Phase Diagrams

This page explains **what phase diagrams do**, lists the **built‑in phase diagrams** (with the exact names you use in `.cfg` files), and shows **how to add new ones**. Phase diagrams are found in `phase.cpp` & `phase.h`.

---

## What a “phase diagram” is (in code)

In MAGRATHEA, each layer (core, mantle, hydrosphere, atmosphere) has a **phase diagram object** that chooses which **equation of state (EOS)** to use at a given pressure–temperature point.

- The class is `PhaseDgm`. It wraps:
  - a **low‑pressure selector** function like `find_phase_Si_default(P, T)` (pointer stored in `phase_lowP`), and
  - an optional **high‑pressure phase list** with fixed transition pressures `start_P` (in **GPa**).  
- The method `find_phase(P, T)` returns a pointer to the EOS to use. **Input pressure `P` is in cgs (microbar)** and is internally converted to **GPa**; `T` is in Kelvin.  
- When the structure integration crosses a phase boundary, the solver calls `find_phase_boundary(...)` to locate the precise transition.

> Practical takeaway: **you provide P (µbar) and T (K); the phase diagram returns the EOS** for that point. The hydro solver then integrates with that EOS.

---

## Built‑in phase diagrams and config names

Use these **strings** in your `.cfg` files (e.g., `core_phasedgm = Fe_default`) to select a phase diagram for that layer.

### Core

| Config name | Selector function | Summary |
|---|---|---|
| `Fe_default` | `find_phase_Fe_default(P,T)` | hcp Fe solid; liquid Fe above melting (Dorogokupets et al. 2017). |
| `Fe_fccbcc` | `find_phase_Fe_fccbcc(P,T)` | Includes **bcc**/**fcc** Fe fields at low P–T and liquid; switches to hcp at high P. |

### Mantle

| Config name | Selector function | Summary |
|---|---|---|
| `Si_default` | `find_phase_Si_default(P,T)` | Upper mantle: **Fo**, **Wds**, **Rwd**, liquid; Lower mantle: **Pv/Brg**, **PPv** with literature transition/melt curves. |
| `Si_simple` | `find_phase_Si_simple(P,T)` | Simplified: **Pv/Brg**, **PPv**, and liquid. |
| `PREM` | `find_phase_PREM(P,T)` | Tabulated PREM‑like mantle at low P; switches to high‑P fits above thresholds. |
| `C_simple` | `find_phase_C_simple(P,T)` | Carbon mantle: **Graphite → Diamond → BC8** transitions. |
| `SiC` | `find_phase_SiC(P,T)` | Silicon carbide mantle: **B3 (zinc blende)** → **B1 (rock salt)** transition with T‑dependence. |

### Hydrosphere (H₂O)

| Config name | Selector function | Summary |
|---|---|---|
| `water_default` | `find_phase_water_default(P,T)` | Comprehensive diagram using SeaFreeze fits, IAPWS gas/liquid, Brown (supercritical), Mazevet (supercritical), and high‑pressure ice (VI, VII, X) regions; includes a variety of triple lines. |
| `water_tabulated` | `find_phase_water_tabulated(P,T)` | Use **AQUA** tabulated EOS directly (Haldemann et al. 2020). |
| `water_legacy` | `find_phase_water_legacy(P,T)` | Simplified diagram following Dunaeva et al. (2010) with major phases (Ih, II, III, V, VI, VII, X, water). |

### Atmosphere

| Config name | Selector function | Summary |
|---|---|---|
| `gas_default` | `find_phase_gas_default(P,T)` | Ideal gas isothermal for **P < 100 bar**; ideal‑gas adiabat deeper. |
| `HHe_tabulated` | `find_phase_HHe_tabulated(P,T)` | Ideal isothermal <100 bar; **Chabrier & Debras (2021) H/He** real‑gas table deeper (Y≈0.275). |

**Example (mode 0):**
```cfg
core_phasedgm = Fe_default
mantle_phasedgm = Si_default
hydro_phasedgm = water_default
atm_phasedgm   = gas_default
```

---

## How selection works (mechanics)

1. **Low‑P region:** your selector function (e.g., `find_phase_Si_default`) returns an EOS based on P–T inequalities (melt curves, transitions, etc.).  
2. **Optional high‑P override:** you can attach a list of EOSs that take over above specified **start pressures** (in GPa). This is configured via `PhaseDgm::set_phase_highP(k, start_P, phase_list)` and is used by the EOS‑variation tools.  
3. `PhaseDgm::find_phase(P,T)` applies the high‑P overrides (if present) and otherwise calls the low‑P selector.  
4. During integration, when the profile steps across a boundary, `find_phase_boundary(...)` finds the bracketed transition pressure and interpolates the boundary temperature consistently.

Notes:
- The **pressure unit inside selectors is GPa** (they all do `P /= 1e10` to convert from µbar).  
- Always handle non‑physical inputs (`P <= 0` or `T <= 0`) by returning `NULL` (the built‑ins do this).

---

## Adding your own phase diagram

You can add a new diagram (say, a **metal‑alloy core** or a **volatile‑rich mantle**) in a few steps.

### 1) Declare a selector in `phase.h`
```cpp
// Header: pick a clear name & layer
EOS* find_phase_MyMantle(double P, double T);
```

### 2) Implement the selector in `phase.cpp`
Follow existing patterns: convert P→GPa, gate invalid inputs, then return EOS pointers based on P–T logic.

```cpp
EOS* find_phase_MyMantle(double P, double T)
{
  if (P <= 0 || T <= 0) return NULL;
  P /= 1E10; // µbar → GPa

  // Example logic: PPv above a T-dependent boundary; else Pv; liquid on melt curve
  if (T > 1830*pow(1+P/4.6, 0.33)) return Si_Liquid_Wolf; // placeholder melt
  if (P > 112.5 + 7E-3*T)         return Si_PPv_Sakai;    // PPv
  return Si_Pv;                                          // Pv/Brg
}
```

> Use EOS pointers already defined in `EOSlist.cpp/h`, or add new EOS first if needed.

### 3) Instantiate a `PhaseDgm`
At the bottom of `phase.cpp`, follow the existing pattern (one line per diagram).

```cpp
PhaseDgm mant5("mantle5", find_phase_MyMantle);
```

### 4) Make it selectable from `.cfg`
In the code that parses config (see the “Set Phase Diagrams” section in `main.cpp`), add a string option so users can write, for example:

```cfg
mantle_phasedgm = MyMantle   # maps to PhaseDgm mant5
```

> Tip: keep the **config string** short and human‑readable and the **PhaseDgm variable** name consistent (e.g., `mant5` above).

### 5) (Optional) Attach high‑P overrides
If you want your diagram to **switch EOS families** above fixed pressures, build an array of EOS pointers and a matching `start_P` list (GPa) and call:
```cpp
mant5.set_phase_highP(k, start_P, highEOSs);
```
This is how the EOS‑scanning utilities update phase diagrams on the fly.

---

## Gotchas and best practices

- **Units:** The solver passes **µbar** to the selector; convert to **GPa** at the top of your function (`P /= 1E10`). Temperatures are in **K**.  
- **Return existing EOS objects:** Use pointers defined in `EOSlist.cpp/h`. Avoid `new` inside selectors.  
- **Keep boundaries monotonic:** When you use polynomials or logs of P, make sure your functions behave well over the domain you expect.  
- **Safeguard extremes:** The built‑ins avoid extrapolating some EOS to very high T by switching to supercritical water EOSes or ideal gas. Mimic those safety checks for new materials.

---

## Quick reference: where things live

- **Phase selection logic:** `phase.cpp` (`find_phase_*` functions for each layer).  
- **Diagram objects:** `phase.cpp` bottom (`PhaseDgm core, mant, water, atm, ...`).  
- **Public API:** `phase.h` (declarations, `PhaseDgm` class, utility methods).  
- **EOS definitions:** `EOSlist.h/.cpp` (70+ built‑in EOS with literature references).  
- **Config mapping:** `main.cpp` (parses `*_phasedgm` strings and sets which `PhaseDgm` to use).

---

## Example: custom “dry super‑Earth” mantle

**Goal:** Pv→PPv transition with a higher‑than‑default boundary and a hotter melt curve.

```cpp
// phase.h
EOS* find_phase_Si_hotPPv(double P, double T);

// phase.cpp
EOS* find_phase_Si_hotPPv(double P, double T)
{
  if (P <= 0 || T <= 0) return NULL;
  P /= 1E10;

  if (T > 2000*pow(1+P/4.6, 0.33)) return Si_Liquid_Wolf;
  if (P > 130.0 + 8E-3*T)          return Si_PPv_Sakai;
  return Si_Pv;
}

// At bottom of phase.cpp:
PhaseDgm mant_hot("mantle_hotPPv", find_phase_Si_hotPPv);
```

**Use in config:**
```cfg
mantle_phasedgm = mantle_hotPPv
```

Replace the literal numbers with your preferred literature fits.

---

## Minimal `.cfg` snippet using built‑ins

```cfg
# Select phase diagrams per layer
core_phasedgm   = Fe_default
mantle_phasedgm = Si_default
hydro_phasedgm  = water_default
atm_phasedgm    = gas_default
```

That’s it — the solver will pick EOSs per layer and evolve P–T–ρ consistently across phase boundaries.
