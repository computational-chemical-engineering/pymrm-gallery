# Catalog B — Reactor models

Sections D–I of the taxonomy: the reactor archetypes. These pages import
closures from [catalog A](catalog-A-foundations.md).

Legend: **S** = structure code, **T** = tier, **P** = build priority
(see [taxonomy.md](taxonomy.md)).

---

## D — Fixed / packed bed reactors

### D1 — The model hierarchy (build these as one connected series)

The pedagogical and practical value here is the *ladder*: the same reaction on
the same bed, modelled at five increasing levels of detail, with the pymrm code
diffed between rungs. No single paper owns this hierarchy; Froment & Bischoff,
*Chemical Reactor Analysis and Design* is the canonical treatment.

| ID | Model | S | T | P |
|---|---|---|---|---|
| D1.1 | 1D pseudo-homogeneous plug flow | S2 | T0 | P1 |
| D1.2 | 1D pseudo-homogeneous + axial dispersion (Danckwerts BCs) | S4 | T0 | P1 |
| D1.3 | 1D heterogeneous (gas + solid phases, film resistance) | S7 | T0 | P1 |
| D1.4 | 1D heterogeneous + intraparticle profiles (nested particle model) | S8 | T0 | **P1** |
| D1.5 | 2D pseudo-homogeneous with radial dispersion + wall heat transfer | S6 | T0 | **P1** |
| D1.6 | 2D heterogeneous | S6+S8 | T1 | P2 |
| D1.7 | Particle-resolved CFD (reference only, not reproduced) | — | T1 | P3 |

D1.4 is the flagship: it is exactly what pymrm's Schur-complement coupling
(`S8`) is for, and it is the rung most textbooks describe but few codes make
easy. The existing teacher solution
`particle-model-coupled-to-column-model.ipynb` is the seed.

### D2 — Thermal behaviour, stability, runaway

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| D2.1 | Barkelew runaway diagram | *Chem. Eng. Prog. Symp. Ser.* 55 (1959) | Runaway boundary in N–S plane | S2 | T0 | P1 |
| D2.2 | Van Welsenaere–Froment criteria | *Chem. Eng. Sci.* 25 (1970) | Two explicit runaway criteria | S2 | T0 | **P1** |
| D2.3 | Morbidelli–Varma generalised sensitivity | *AIChE J* 28 (1982), *Chem. Eng. Sci.* (1988) | Parametric sensitivity criterion | S2 | T1 | P2 |
| D2.4 | Multiplicity & ignition/extinction | Hlavacek & Votruba; Uppal–Ray–Poore (CSTR) | Steady-state multiplicity, hysteresis | S1/S3 | T0 | **P1** |
| D2.5 | Vortmeyer–Schaefer equivalence | *Chem. Eng. Sci.* 29 (1974) | Axial dispersion ↔ two-phase heat model | S4/S7 | T1 | P2 |
| D2.6 | Wrong-way behaviour / creeping fronts | Mehta, Sams & Luss (1981) | Transient temperature excursion on cooling | S4 | T1 | P2 |
| D2.7 | Adiabatic reactor with recycle / autothermal | — | Ignition curves, multiplicity | S2+S10 | T1 | P2 |

**Gallery angle for D2.** D2.2 and D2.4 are outstanding pymrm demonstrations
because the *interesting* result (the runaway boundary, the hysteresis loop) is
a sweep over many solves — trivially parallel, visually striking, and something
the original papers could only sketch with a handful of points.

### D3 — Specific industrial fixed-bed cases (each = one page with plant/lab data)

| ID | Case | Reference anchor | T | P |
|---|---|---|---|---|
| D3.1 | Steam methane reformer tube | Xu & Froment (1989) + Froment reformer papers | T1 | **P1** |
| D3.2 | Ammonia synthesis converter (multi-bed, quench/intercooled) | Dyson & Simon (1968); Temkin | T0 | P1 |
| D3.3 | Methanol synthesis reactor (Lurgi/quasi-isothermal) | Graaf; Vanden Bussche & Froment | T1 | P1 |
| D3.4 | o-Xylene → phthalic anhydride multitubular | Froment (1967) | T0 | **P1** |
| D3.5 | Ethylene oxide multitubular | Borman & Westerterp (1995) | T1 | P1 |
| D3.6 | SO2 oxidation multi-bed | Classic; abundant data | T0 | P1 |
| D3.7 | Fixed-bed FT (ARGE-type) | Post et al.; Jess & Kern | T1 | P2 |
| D3.8 | Autothermal reformer / catalytic partial oxidation | De Groote & Froment (1996) | T1 | P2 |
| D3.9 | Hydrotreater trickle-bed (see also G) | Korsten & Hoffmann (1996) | T1 | P2 |
| D3.10 | Reverse-flow SO2/VOC oxidation reactor | Matros & Bunimovich, *Catal. Rev.* 38 (1996) | T1 | P2 |

---

## E — Fluidised beds

### E1 — Hydrodynamic foundations

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| E1.1 | Two-phase theory | Toomey & Johnstone (1952) | Excess gas goes to bubbles | S7 | T0 | P1 |
| E1.2 | Davidson–Harrison bubble | *Fluidised Particles* (1963) | Bubble rise, cloud, throughflow | S3 | T0 | P1 |
| E1.3 | Bubble growth / coalescence | Darton et al. (1977); Werther; Mori & Wen (1975) | d_b(h) | — | T1 | P1 |
| E1.4 | Sit–Grace interphase mass transfer | *Chem. Eng. Sci.* 36 (1981) | K_be | — | T1 | P1 |
| E1.5 | Solids mixing / dispersion | May (1959); van Deemter (1961) | Axial solids dispersion | S4 | T0 | P2 |
| E1.6 | Elutriation & TDH | Geldart; Wen & Chen | Entrainment above the bed | S1 | T1 | P2 |

### E2 — Fluidised bed reactor models

| ID | Model | Canonical reference | S | T | P |
|---|---|---|---|---|---|
| E2.1 | **Kunii–Levenspiel bubbling bed** | *IEC Proc. Des. Dev.* 7:481 (1968); *Fluidization Engineering* (1991) | S7 | T0 | **P1** |
| E2.2 | Kato–Wen bubble assemblage | *Chem. Eng. Sci.* 24 (1969) | S7 | T1 | P2 |
| E2.3 | Grace generalised three-phase | *Can. J. Chem. Eng.* 64 (1986) | S7 | T1 | P1 |
| E2.4 | Thompson–Bi–Grace generic model | *Chem. Eng. Sci.* 54 (1999) | S7 | T1 | P2 |
| E2.5 | Werther / Hilligardt–Werther | (1980s) | S7 | T1 | P2 |
| E2.6 | Turbulent-regime reactor model | Cui, Grace; Bi & Grace regime maps | S7 | T2 | P2 |
| E2.7 | CFB core–annulus | Berruti & Kalogerakis (1989); Pugsley & Berruti (1996); Senior & Brereton (1992) | S7 | T1 | P2 |
| E2.8 | CFB riser 1.5D with cluster | Various | S7 | T2 | P3 |
| E2.9 | Two-fluid / KTGF (reference only) | Lun et al. (1984); Gidaspow (1994) | — | T1 | P3 |
| E2.10 | Membrane-assisted fluidised bed | Deshmukh, van Sint Annaland & Kuipers (TU/e) | S7+H | T2 | **P1** |
| E2.11 | Spouted bed | Mathur & Epstein | S7 | T2 | P3 |
| E2.12 | Fluidised bed FT / MTO / polyolefin | Various | S7 | T2 | P2 |

**Data note (verified).** The NETL/PSRI *CFD Challenge Problems* release
experimental data for a 0.305 m × 16 m CFB riser and a 0.92 m × 9 m bubbling
bed, publicly at `mfix.netl.doe.gov`. This is the best open validation dataset
for section E and should anchor E2.7 and E2.9.

---

## F — Gas–liquid and slurry reactors

### F1 — Bubble column hydrodynamics and transfer

| ID | Model | Canonical reference | Predicts | T | P |
|---|---|---|---|---|---|
| F1.1 | Akita–Yoshida | *IEC Proc. Des. Dev.* 12 (1973), 13 (1974) | ε_g, k_L a, d_b | T0 | P1 |
| F1.2 | Hikita et al. | *Chem. Eng. J.* 20 (1980) | ε_g, k_L a incl. liquid properties | T1 | P1 |
| F1.3 | Wilkinson et al. | *AIChE J* 38 (1992) | ε_g with scale/pressure effects | T1 | P1 |
| F1.4 | **Krishna–Ellenberger two-bubble-class** | *AIChE J* 42:2627 (1996) | Small + large bubble holdup, churn-turbulent | T1 | **P1** |
| F1.5 | Vandu–Krishna | *Chem. Eng. Process.* (2004) | k_L a in churn-turbulent | T1 | P1 |
| F1.6 | Deckwer axial dispersion correlations | *Bubble Column Reactors* (1992) | D_ax,L and D_ax,G | T0 | P1 |
| F1.7 | Joshi–Sharma circulation cell | *Trans. IChemE* 57 (1979) | Liquid circulation pattern | T1 | P2 |
| F1.8 | Shah–Kelkar–Godbole–Deckwer review | *AIChE J* 28 (1982) | Consolidated design parameters | T0 | P1 |
| F1.9 | Regime transition (homogeneous→churn) | Krishna; Zahradník et al. | Transition ε_g, u_g | T1 | P1 |

### F2 — Gas–liquid reactor models

| ID | Model | Canonical reference | S | T | P |
|---|---|---|---|---|---|
| F2.1 | Axial dispersion two-phase bubble column | Deckwer | S7+S4 | T0 | **P1** |
| F2.2 | Two-bubble-class reactor model | Krishna, de Swart, Ellenberger | S7 | T1 | P1 |
| F2.3 | Slurry bubble column for FT | Maretto & Krishna, *Catal. Today* 52 (1999) | S7+S8 | T1 | **P1** |
| F2.4 | Krishna–Sie FT reactor selection | *Chem. Eng. Sci.* 55 (2000) | S7 | T1 | P2 |
| F2.5 | Airlift loop reactor | Chisti | S7 | T1 | P2 |
| F2.6 | Stirred tank G–L reactor | Van 't Riet; Bakker | S1+S7 | T1 | P1 |
| F2.7 | 2D circulation bubble column | (teacher exercise exists) | S6+S7 | T2 | P1 |
| F2.8 | Bubble column with PBE (bubble size distribution) | Luo & Svendsen; Prince & Blanch | S11 | T2 | P3 |

### F3 — Absorption with chemical reaction

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| F3.1 | Hatta number regimes | Hatta (1932) | Slow/fast/instantaneous regimes | S3 | T0 | **P1** |
| F3.2 | Van Krevelen–Hoftijzer enhancement | *Rec. Trav. Chim.* 67 (1948) | E for 2nd-order reaction | S3 | T0 | P1 |
| F3.3 | DeCoursey enhancement | *Chem. Eng. Sci.* 29 (1974) | Explicit E approximation | S3 | T1 | P1 |
| F3.4 | Danckwerts gas–liquid reaction theory | *Gas-Liquid Reactions* (1970) | Surface renewal + reaction | S4 | T0 | P1 |
| F3.5 | CO2–amine absorption | Versteeg & van Swaaij (1988); Danckwerts & Sharma | Rate + speciation | S3+S10 | T1 | **P1** |
| F3.6 | Rate-based (nonequilibrium stage) column | Krishnamurthy & Taylor, *AIChE J* 31 (1985) | Column profiles without stage efficiency | S7+S9 | T1 | P2 |
| F3.7 | Reactive distillation | Taylor & Krishna, *Chem. Eng. Sci.* 55 (2000) | Coupled reaction/separation | S7+S10 | T1 | P3 |

**Gallery angle for F3.** F3.1 is the ideal `S3` teaching page: the enhancement
factor asymptotes are analytical, the full film problem is a two-line pymrm
solve, and overlaying them shows exactly where each asymptote fails.

---

## G — Three-phase packed beds (trickle bed)

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| G1.1 | Larkins pressure drop | Larkins, White & Jeffrey, *AIChE J* 7 (1961) | Δp via Lockhart–Martinelli parameter | S3 | T0 | P1 |
| G1.2 | Lockhart–Martinelli | *Chem. Eng. Prog.* 45 (1949) | Two-phase Δp framework | S3 | T0 | P1 |
| G1.3 | Sáez–Carbonell relative permeability | *AIChE J* 31 (1985) | Δp + holdup, mechanistic | S3 | T1 | P1 |
| G1.4 | Attou–Boyer–Ferschneider | *Chem. Eng. Sci.* 54 (1999) | Momentum-based two-fluid Δp | S3 | T1 | P2 |
| G1.5 | Al-Dahhan–Duduković high-pressure | *Chem. Eng. Sci.* 50 (1995) | Δp, holdup, wetting at pressure | S3 | T1 | P1 |
| G1.6 | El-Hisnawi / Mills–Duduković wetting efficiency | (1982); *AIChE J* 27 (1981) | η_CE correlations | — | T1 | P1 |
| G1.7 | Wammes–Westerterp high-pressure hydrodynamics | *Chem. Eng. Sci.* (1991) | Regime + holdup at pressure (TU/e lineage) | S3 | T1 | P1 |
| G1.8 | Trickle-bed reactor with partial wetting | Ramachandran & Chaudhari (1983); Herskowitz & Smith, *AIChE J* 29 (1983) | Apparent rate with incomplete wetting | S8 | T1 | **P1** |
| G1.9 | Korsten–Hoffmann hydrotreater | *AIChE J* 42 (1996) | Industrial HDS trickle bed | S7 | T1 | P2 |
| G1.10 | Iliuta–Larachi multiphase Eulerian | (2000s) | Unified hydrodynamics | S7 | T2 | P3 |
| G1.11 | Periodic/liquid-induced operation | Silveston & Hanika | Performance under flow modulation | S4 | T2 | P3 |

---

## H — Membrane reactors and membrane separations

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| H1.1 | Sieverts law H2 permeation | Sieverts; exponent analysis by Caravella et al. (2010) | J_H2 ∝ (p^n − p^n), n≈0.5–1 | S7 | T0 | **P1** |
| H1.2 | Ward–Dao permeation mechanism | *J. Membr. Sci.* 153 (1999) | Which permeation step limits | S3 | T1 | P2 |
| H1.3 | Concentration polarisation in Pd membranes | Caravella, Barbieri & Drioli, *Sep. Purif. Technol.* (2009); *IJHE* (2015) | Loss of driving force along module | S6+S7 | T1 | **P1** |
| H1.4 | Itoh Pd membrane dehydrogenation reactor | *AIChE J* 33 (1987) | Equilibrium shift | S7 | T0 | P1 |
| H1.5 | Packed-bed membrane reactor | Tsotsis, Champagnie et al. (1992) | Conversion beyond equilibrium | S7 | T1 | P1 |
| H1.6 | Gallucci et al. membrane reactor review models | *Chem. Eng. Sci.* 92 (2013) (TU/e) | Consolidated MR modelling | S7 | T1 | **P1** |
| H1.7 | Wijmans–Baker solution–diffusion | *J. Membr. Sci.* 107:1 (1995) | Unified transport for dense membranes | S3 | T0 | **P1** |
| H1.8 | Robeson upper bound | *J. Membr. Sci.* 62 (1991), 320 (2008) | Permeability/selectivity trade-off | — | T0 | P1 |
| H1.9 | Maxwell–Stefan zeolite membrane | Krishna & van den Broeke (1995); Krishna & Baur (2003) | Loading-dependent multicomponent permeation | S9 | T1 | P2 |
| H1.10 | Perovskite oxygen transport (Wagner) | Bouwmeester & Burggraaf | Ambipolar O2 flux | S3 | T1 | P2 |
| H1.11 | Membrane-assisted autothermal reforming | van Sint Annaland group (TU/e) | Integrated H2 production | S7 | T2 | **P1** |
| H1.12 | Ammonia synthesis membrane reactor | This suite's `ammonia_synthesis_reactor` project | 2D axisymmetric, coupled p–T–c | S6+S7+S10 | T3 | *deferred* |
| H1.13 | RO/NF concentration polarisation + film theory | — | Flux decline, rejection | S3 | T1 | P2 |
| H1.14 | Electrodialysis / ion-exchange membrane | — | Nernst–Planck with electroneutrality | S10 | T2 | P3 |

**Gallery angle for H.** H1.12 is in-house work but **not yet published**, so it
is deferred under the published-work-only policy
([blueprint §9](blueprint.md#published-work-only-policy)) and returns when the
paper appears. The first wave for this section is therefore H1.1 (Sieverts
permeation — well-tabulated data, and the pressure exponent *n* is itself a
published controversy worth plotting against data) paired with H1.4 (Itoh's
membrane dehydrogenation reactor, AIChE J 1987 — conversion beyond the
equilibrium limit, which needs only H1.1's closure plus a 1D reactor model).

H1.3 remains the best "pymrm improves the original" candidate in this section:
concentration polarisation is usually estimated with a 1D film coefficient, and
pymrm can resolve the actual 2D boundary layer (`S6`) and show the error.

---

## I — Structured and intensified reactors

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| I1.1 | Monolith channel with washcoat | Hawthorn (1974); Groppi & Tronconi | Sh(z), washcoat diffusion, conversion | S6+S8 | T0 | **P1** |
| I1.2 | Catalytic converter transient | Oh & Cavendish, *IEC Prod. Res. Dev.* 21 (1982) | Light-off, cold start | S4+S7 | T0 | **P1** |
| I1.3 | Heck–Wei–Katzer monolith | *AIChE J* 22 (1976) | Monolith reactor performance | S6 | T0 | P1 |
| I1.4 | Groppi–Tronconi heat transfer in honeycombs | *AIChE J* / *Chem. Eng. Sci.* (1996–2000) | Conductive monolith thermal behaviour | S6 | T1 | P2 |
| I1.5 | Open-cell foam reactor | Boger & Heibel; Giani, Groppi & Tronconi (2005) | Δp, transfer in foams | S6 | T2 | P2 |
| I1.6 | Microreactor / Graetz with reaction | Hessel, Löwe | Isothermality, short contact time | S6 | T1 | P1 |
| I1.7 | Reverse-flow reactor | Matros & Bunimovich, *Catal. Rev.* 38 (1996) | Heat trapping, autothermal VOC oxidation | S4 | T1 | **P1** |
| I1.8 | Chemical looping combustion reactor | Adánez et al., *PECS* 38 (2012); Lyngfelt | Fuel/air reactor coupling | S7+S12 | T1 | P2 |
| I1.9 | Sorption-enhanced reforming | Hufton, Mayorga & Sircar, *AIChE J* 45 (1999) | Reaction + in-situ CO2 capture | S4+S12 | T1 | **P1** |
| I1.10 | Simulated moving bed (reactive) | Broughton & Gerhold (1961); triangle theory (Mazzotti, Storti, Morbidelli) | Countercurrent separation/reaction | S5+S7 | T1 | P2 |
| I1.11 | Electrified (Joule-heated) reforming | Wismann et al., *Science* 364 (2019) | Wall-heated compact reformer | S6 | T3 | P2 |
| I1.12 | Rotating packed bed / HiGee | Ramshaw | Intensified mass transfer | S3 | T2 | P3 |
| I1.13 | Catalytic plate/heat-exchanger reactor | — | Coupled endo/exothermic | S6+S7 | T2 | P2 |
