# Run Modes

Document each solver mode here. For each mode include:

- **What it does** (one-liner).
- **Required config keys** (with units and defaults).
- **Inputs/Outputs** (filenames, columns).
- **Common pitfalls** (non-convergence, bad initial guesses, unit gotchas).

## Mode 0 — Hydrostatic Structure
- **Purpose:** Solve hydrostatic equilibrium for a specified composition/phase stack.
- **Key settings:** `mode = 0`, `layers`, `t_profile = isentropic|isothermal`, `phase_diagram` file, `eos_list` file.
- **Outputs:** `profile_*.txt` (r, P, T, ρ), total M and R.
- **Pitfalls:** EOS validity (P–T ranges), step size too large.

## Mode 1 — Composition Finder
- **Purpose:** Invert for layer compositions that hit target bulk (M, R, ρ̄, etc.).
- **Key settings:** `mode = 1`, `targets.mass`, `targets.radius`, `search.bounds`.
- **Outputs:** `solutions.csv` with best-fit mixtures.
- **Pitfalls:** Degeneracies; provide bounds/priors; check solution stability.

## Mode 3 — (placeholder)
Describe what Mode 3 does and its major switches.

## Mode 4 — (placeholder)
Describe what Mode 4 does and its major switches.

> Replace the placeholders with the exact names from your code/config parser.
