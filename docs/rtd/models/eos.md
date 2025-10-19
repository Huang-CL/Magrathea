# Equation of State (EOS) Library

Describe your materials, sources, and validity ranges. Include how to add a new EOS.

## EOS Catalog

| Material | Symbol | Phase(s) | Validity (P–T) | Source/Ref | Notes |
|---------:|:------:|---------:|----------------|-----------|-------|
| Forsterite | Fo | olivine | e.g., 1–150 GPa, 300–2500 K | Doe+2020 | --- |
| Enstatite | En | pyroxene | ... | ... | ... |
| Iron | Fe | bcc/hcp | ... | ... | ... |

> Fill with your actual materials and references.

## EOS Files and Formats
Explain where `EOSlist.*` lives, file format, units, and how the code looks up materials.

## Adding a New EOS
1. Add parameters to the EOS file or implement a new model class.
2. Document units and references.
3. Rebuild and verify with a small test (compare density at a few P–T points).

## Phase Diagrams
- How phase is selected vs P–T.
- How to override phase regions for experiments.
- Edge cases near boundaries.
