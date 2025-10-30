# Equation of State (EOS) Library

MAGRATHEA includes a comprehensive set of **Equations of State (EOS)** for solid, liquid, and gaseous planetary materials.  
This page describes all EOS entries defined in `EOSlist.cpp`, with their variable names, references, and notes.  

---

## Default EOS Set

These are the **default EOSs** used by the standard phase diagrams (`Fe_default`, `Si_default`, `water_default`, and `gas_default`).

- **Fe_hcp3** — *Epsilon (hcp) iron*, Dorogokupets et al. (2017, *Scientific Reports*). Default solid core EOS.  
- **Fe_liquid** — *Liquid iron*, Dorogokupets et al. (2017). Default molten core EOS.  
- **Si_Pv** — *Perovskite (Bridgmanite) MgSiO₃*, Oganov & Ono (2004, *Nature*). Default mantle EOS.  
- **Si_PPv_Sakai** — *Post‑Perovskite MgSiO₃*, Sakai et al. (2016, *Scientific Reports*). Default deep mantle EOS.  
- **Si_Liquid_Wolf** — *Liquid MgSiO₃*, Wolf & Bower (2018). Default silicate melt EOS.  
- **Water_IAPWS** — *Liquid water*, Wagner & Prub (2022, *JPCRD*). Default low‑pressure hydrosphere EOS.  
- **Water_Brown** — *Supercritical water*, Brown (2018, *Fluid Phase Equilibria*). Default high‑temperature fluid EOS.  
- **IceVI_Bezacier**, **IceVII_Bezacier**, **IceX** — *High‑pressure ice polymorphs*, Bezacier et al. (2014), Fei et al. (1993). Default ices.  
- **Gas** — *Adiabatic ideal gas*, default atmospheric EOS.  
- **Gas_iso** — *Isothermal ideal gas*, used for fixed‑temperature atmospheres.

---

## Iron EOS

- **Fe_liquid**, Dorogokupets et al. 2017 — Liquid Fe, primary melt EOS.  
- **Fe_liquid2**, Anderson & Ahrens 1994 — Alternate liquid Fe fit.  
- **Fe_bcc**, Dorogokupets et al. 2017 — α‑Fe (bcc).  
- **Fe_fcc**, Dorogokupets et al. 2017 — γ‑Fe (fcc).  
- **Fe_hcp**, Smith et al. 2018 — ε‑Fe (hcp), solid core EOS.  
- **Fe_hcp2**, Bouchet et al. 2013 — Alternate ε‑Fe model.  
- **Fe_hcp3**, Dorogokupets et al. 2017 — Default solid Fe EOS.  
- **Fe_7Si**, Wicks et al. 2018 — Fe‑7wt%Si alloy.  
- **Fe_15Si**, Wicks et al. 2018 — Fe‑15wt%Si alloy.  
- **Fe_Seager**, Seager et al. 2007 — Tabulated Fe EOS.  
- **Fe_Dummy** — Placeholder Fe EOS.

---

## Silicate and Mantle EOS

- **Si_liquid**, Mosenfelder et al. 2009 — Liquid MgSiO₃.  
- **Si_Liquid_Wolf**, Wolf & Bower 2018 — Default silicate melt (RTpress form).  
- **Fo**, **Wds**, **Rwd**, **Akm**, Dorogokupets et al. 2015 — Olivine, Wadsleyite, Ringwoodite, Akimotoite.  
- **Fo_Sotin**, **En**, **Mw**, Sotin et al. 2007 — Simplified mantle phases.  
- **Si_Pv**, Oganov & Ono 2004 — Perovskite (Bridgmanite).  
- **Si_Pv_Shim**, Shim & Duffy 2000 — Alternate perovskite.  
- **Pv_Doro**, Dorogokupets et al. 2015 — Bridgmanite variant.  
- **Si_PPv_Sakai**, Sakai et al. 2016 — Default post‑perovskite.  
- **Si_PPv**, Oganov & Ono 2004 — Alternate post‑perovskite.  
- **PPv_Doro**, Dorogokupets et al. 2015 — Post‑perovskite variant.  
- **Si_PREM**, Stacey & Davis 2008 — Tabulated PREM silicate.  
- **Si_BM2fit**, Zeng 2016 — Simplified PREM extrapolation.  
- **Si_Seager**, Seager et al. 2007 — Tabulated silicate.  
- **Si_Dummy** — Placeholder silicate EOS.

### Extended Silicate Library (Stixrude & Lithgow‑Bertelloni 2011)

Covers mineral families used for complex mantle models.

**Anorthite**, **Spinel**, **Fayalite**, **Fe_Wadsleyite**, **Fe_Ringwoodite**, **Enstatite**, **Ferrosilite**, **Diopside**, **Hedenbergite**, **HP_clinopyroxene**, **Mg_Akimotoite**, **Fe_Akimotoite**, **Pyrope**, **Mg_Majorite**, **Almandine**, **Grossular**, **Quartz**, **Coesite**, **Stishovite**, **Seifertite**, **Fe_Post_Perovskite**, **Fe_Perovskite**, **HP_Clinoferrosilite**, **Periclase**, **Wustite**, **Kyanite**, **Nepheline** — All use Vinet EOS fits from *Stixrude & Lithgow‑Bertelloni (2011)*.

---

## Water, Ice, and Hydrosphere EOS

- **Water_SF**, Bollengier et al. 2019 — SeaFreeze table.  
- **Water_IAPWS**, Wagner & Prub 2022 — IAPWS R6‑95 water.  
- **Water_ExoPlex**, unknown source — ExoPlex liquid water.  
- **Water**, Valencia et al. 2007 — Simplified water EOS.  
- **Water_Brown**, Brown 2018 — Supercritical water.  
- **Water_sc_Mazevet**, Mazevet et al. 2019 — A&A 621, high‑T plasma water.  
- **IceIh_SF**, **IceII_SF**, **IceIII_SF**, **IceV_SF**, **IceVI_SF**, Journaux et al. 2020 — Tabulated ice phases.  
- **IceIh**, Feistel & Wagner 2006 — Solid ice Ih.  
- **IceIh_ExoPlex**, ExoPlex variant.  
- **IceVI_ExoPlex**, Bezacier et al. 2014 — Ice VI.  
- **IceVI_Bezacier**, Bezacier et al. 2014 — Default ice VI EOS.  
- **IceVII_Bezacier**, Bezacier et al. 2014 — Default ice VII EOS.  
- **IceVII_FFH2004**, Frank, Fei & Hu 2004 — Ice VII, Vinet fit.  
- **IceVII_FFH2004fit/BM/T**, FFH2004 variations.  
- **IceVII_Fei**, Fei et al. 1993 — Ice VII model with melt curve.  
- **IceVIIp**, **IceVII**, **IceVII (Grande)** — Alternate high‑pressure ice fits.  
- **IceX**, Grande — Default high‑P ice.  
- **IceX_HS**, Hermann & Schwerdtfeger 2011 — PRL fit.  
- **Ice_Seager**, Seager et al. 2007 — Tabulated ice.  
- **H2O_AQUA**, Haldemann et al. 2022 — AQUA tabulated EOS.  
- **H2O_SeaFreeze**, Journaux et al. 2020 — Combined water/ice table.  
- **Ice_Dummy**, Placeholder EOS.  
- **IceZeng2013FFH/FMNR**, Zeng 2013 — Modified Zeng ice fits.

---

## Atmosphere and Gas EOS

- **Gas** — Ideal adiabatic gas (default).  
- **Gas_iso** — Isothermal ideal gas.  
- **vdW_H2**, **vdW_He**, **vdW_H2O**, **vdW_CH4**, **vdW_NH3**, **vdW_CO2** — Van der Waals EOS for major planetary gases.  
- **watervapor** — Ideal water vapor.  
- **Water_Vap_IAPWS** — IAPWS vapor EOS.  
- **Gas_hhe** — H/He mixture, Chabrier & Debras (2021).

---

## Carbon and Carbides

- **Graph**, Seager et al. 2007 — Graphite, Birch‑Murnaghan EOS.  
- **Diam**, Benedict et al. 2018 — Diamond, Vinet EOS.  
- **BC8**, Benedict et al. 2018 — BC8 carbon polymorph.  
- **SiC_B3_Vinet**, **SiC_B1_Vinet**, Miozzi et al. 2018 — Silicon carbide polymorphs.

---

## Calibration and Miscellaneous Metals

- **Gold_BM3**, Heinz & Jeanloz 1983 — Gold calibration EOS (BM3).  
- **Gold**, Matsui et al. 2010 — Alternate gold fit.  
- **Plat**, Matsui et al. 2009 — Platinum standard.  

---

## EOS Parameter Schema

Each EOS uses a numerical parameter list `{{index, value}}`, where key fields include:

| Index | Parameter | Meaning |
|-------|------------|----------|
| 1 | V₀ | Molar volume (cm³/mol) |
| 2 | K₀ | Bulk modulus (GPa) |
| 3 | K₀′ | Pressure derivative |
| 5 | M | Molar mass (g/mol) |
| 7 | Θ₀ | Debye/Einstein temperature (K) |
| 8–11 | γ parameters | Grüneisen model coefficients |
| 14 | n | Atoms per formula unit |
| 16–19 | α parameters | Thermal expansion coefficients |
| 20–22 | cₚ coefficients | Heat capacity model |
| 23 | Debye flag | +Debye, −Einstein |
| 33–34 | Van der Waals a,b | Gas constants |

Full definitions appear in the top comments of `EOSlist.cpp`.

---

## Implementation Notes

- All EOS return pressure–density–temperature relations compatible with the hydrostatic solver.  
- Most solids use **Birch–Murnaghan (BM3)** or **Vinet** forms with Debye/Einstein thermal corrections.  
- Liquids may use **RTpress** or **tabulated** equations.  
- Atmosphere and gas EOS use **ideal or van der Waals** forms.  
- Tabulated EOS load from the `./tabulated/` directory.  
- Units: pressure in **GPa**, density in **g cm⁻³**, temperature in **K**.

---

_This document summarizes all EOS definitions contained in `EOSlist.cpp`._
