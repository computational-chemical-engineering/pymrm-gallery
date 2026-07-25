# Catalog C — Adjacent unit operations and neighbouring fields

Section J of the taxonomy. These are not reactors, but they are (a) modelled
with the *same* pymrm structures, (b) heavily cited, and (c) where a large
population of potential gallery users actually works. Including them roughly
doubles the gallery's audience for maybe 30% more effort, because most reuse
`S4`/`S5`/`S11` machinery already built for reactors.

Legend: **S** = structure code, **T** = tier, **P** = build priority.

---

## J1 — Adsorption, chromatography, PSA

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| J1.1 | Langmuir isotherm | Langmuir, *JACS* 40 (1918) | Monolayer equilibrium | S1 | T0 | P1 |
| J1.2 | Freundlich / Sips / Tóth / Dubinin | — | Heterogeneous-surface isotherms | S1 | T0 | P1 |
| J1.3 | BET | Brunauer, Emmett & Teller (1938) | Multilayer, surface area | S1 | T0 | P1 |
| J1.4 | IAST | Myers & Prausnitz, *AIChE J* 11 (1965) | Mixture equilibrium from pure isotherms | S10 | T0 | **P1** |
| J1.5 | Glueckauf LDF | *Trans. Faraday Soc.* 51 (1955) | Lumped intraparticle rate (15D/R²) | S4 | T0 | **P1** |
| J1.6 | Thomas / Bohart–Adams / Yoon–Nelson | (1944); (1920) | Analytical breakthrough curves | S4 | T0 | P1 |
| J1.7 | Rosen solution | *J. Chem. Phys.* 20 (1952) | Breakthrough with pore diffusion | S4 | T0 | P2 |
| J1.8 | General rate model | Guiochon; Gu | Full film + pore + surface resolution | S4+S8 | T1 | **P1** |
| J1.9 | Equilibrium theory / shock–wave | Rhee, Aris & Amundson | Front sharpening, wave interactions | S5 | T0 | P2 |
| J1.10 | van Deemter plate height | van Deemter, Zuiderweg & Klinkenberg, *Chem. Eng. Sci.* 5 (1956) | HETP vs velocity, A/B/C terms | S4 | T0 | **P1** |
| J1.11 | Skarstrom PSA cycle | (1960) | Cyclic steady state, purity/recovery | S4+S7 | T0 | P2 |
| J1.12 | TSA / VSA / VPSA cycles | Ruthven; Sircar | Cyclic separation performance | S4+S7 | T1 | P2 |
| J1.13 | SMB triangle theory | Storti, Mazzotti, Morbidelli (1993) | Operating region for complete separation | S5+S7 | T1 | P2 |
| J1.14 | Zeolite micropore diffusion | Kärger & Ruthven | Loading-dependent D, single-file | S9 | T1 | P2 |
| J1.15 | CO2 direct air capture sorbent bed | Modern (Sinha, Realff et al.) | Cyclic capture under dilute feed | S4 | T3 | P2 |

**Gallery angle for J1.** J1.5+J1.6+J1.8 form a natural triptych: analytical
breakthrough → LDF → full general rate model, all on the *same* breakthrough
dataset, showing exactly what each simplification costs. This is one of the
cleanest "ladder" stories in the catalog and reuses D1's structure.

## J2 — Crystallisation, particulate processes, population balances

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| J2.1 | Population balance equation | Hulburt & Katz, *Chem. Eng. Sci.* 19 (1964); Ramkrishna (2000) | Number density evolution | S11 | T0 | **P1** |
| J2.2 | MSMPR steady-state CSD | Randolph & Larson (1971, 1988) | Exponential CSD, growth/nucleation from slope | S11 | T0 | **P1** |
| J2.3 | Nucleation & growth kinetics | Mullin; Mersmann; Nývlt | B and G from supersaturation | S1 | T0 | P1 |
| J2.4 | Kumar–Ramkrishna fixed pivot / cell average | *Chem. Eng. Sci.* 51 (1996) | Conservative PBE discretisation | S11 | T1 | P2 |
| J2.5 | QMOM / DQMOM | McGraw (1997); Marchisio & Fox (2005) | Moment closure for PBE | S1 | T1 | P2 |
| J2.6 | Smoluchowski aggregation | (1917) | Coagulation kernel dynamics | S11 | T0 | P1 |
| J2.7 | Coulaloglou–Tavlarides | *Chem. Eng. Sci.* 32 (1977) | Drop breakage/coalescence in stirred tanks | S11 | T0 | P2 |
| J2.8 | Luo–Svendsen breakup | *AIChE J* 42 (1996) | Bubble/drop breakup kernel | S11 | T1 | P2 |
| J2.9 | Prince–Blanch coalescence | *AIChE J* 36 (1990) | Coalescence kernel | S11 | T1 | P2 |
| J2.10 | Granulation / Iveson–Litster regime map | *Powder Technol.* (1998) | Growth regimes | S11 | T2 | P3 |
| J2.11 | Attrition / breakage in beds | — | Fines generation | S11 | T2 | P3 |

**Note on S11.** A population balance with a growth term is mathematically a
convection equation along an internal coordinate. That means pymrm's TVD
machinery (`interp_cntr_to_stagg_tvd`) applies *directly* — size is just another
axis. This is an under-exploited fit and worth making explicit, because sharp
CSD fronts are exactly where naive PBE discretisations smear.

## J3 — Electrochemical reactors and devices

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| J3.1 | Butler–Volmer / Tafel | — | Charge-transfer kinetics | S1 | T0 | P1 |
| J3.2 | Nernst–Planck + electroneutrality | Newman & Thomas-Alyea, *Electrochemical Systems* | Ion transport with migration | S10 | T0 | **P1** |
| J3.3 | Porous electrode theory | Newman & Tobias, *JES* 109 (1962) | Distributed reaction in electrodes | S3+S10 | T0 | **P1** |
| J3.4 | **Doyle–Fuller–Newman (P2D)** | Doyle, Fuller & Newman, *JES* 140:1526 (1993) | Li-ion cell voltage/concentration | S8+S10 | T1 | **P1** |
| J3.5 | Single Particle Model / SPMe | Various reductions of J3.4 | Fast approximation, comparison page | S3 | T1 | P1 |
| J3.6 | Springer–Zawodzinski–Gottesfeld | *JES* 138 (1991) | PEM fuel cell membrane/water | S3+S10 | T1 | P1 |
| J3.7 | Bernardi–Verbrugge | *AIChE J* 37 (1991) / *JES* (1992) | PEFC transport | S3+S10 | T1 | P2 |
| J3.8 | SOFC electrochemical model | Achenbach (1994); Costamagna | Polarisation, thermal coupling | S6+S10 | T1 | P2 |
| J3.9 | Water electrolyser (PEM/alkaline) | Marangio et al.; Carmo review | Polarisation curve, efficiency | S3+S10 | T2 | P1 |
| J3.10 | Redox flow battery | Shah et al., *Electrochim. Acta* 53 (2008) | Cell + tank dynamics | S4+S10 | T2 | P2 |
| J3.11 | CO2 electroreduction cell | Weng, Bell & Weber, *PCCP* 20 (2018) | Selectivity, local pH, boundary layer | S3+S10 | T3 | **P1** |
| J3.12 | Electrodialysis / capacitive deionisation | — | Desalination performance | S10 | T3 | P3 |

**Gallery angle for J3.** J3.4 is one of the most-cited models in the whole
catalog and is *structurally* a reactor–particle coupling (`S8`) — the same Schur
machinery as D1.4. Demonstrating that pymrm solves DFN with the identical code
pattern as a packed-bed–pellet model is a strong statement about the library's
generality, and it reaches a large audience outside classical reaction engineering.

## J4 — Biochemical reactors

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| J4.1 | Monod | Monod (1949) | μ(S) | S1 | T0 | P1 |
| J4.2 | Haldane / Andrews substrate inhibition | Andrews (1968) | Inhibition at high S | S1 | T0 | P1 |
| J4.3 | Contois / Moser / Tessier | — | Alternative growth laws | S1 | T1 | P2 |
| J4.4 | Luedeking–Piret | *J. Biochem. Microbiol.* 1 (1959) | Growth/non-growth product formation | S1 | T0 | P1 |
| J4.5 | Herbert–Pirt maintenance | — | Substrate partition to maintenance | S1 | T0 | P1 |
| J4.6 | Michaelis–Menten / Briggs–Haldane | (1913); (1925) | Enzyme kinetics | S1 | T0 | P1 |
| J4.7 | Immobilised enzyme/biofilm pellet | Atkinson; Characklis | Diffusion-limited biocatalysis (= B1.1 with MM kinetics) | S3 | T0 | **P1** |
| J4.8 | **ASM1 activated sludge** | Henze et al., IWA (1987) | Wastewater plant dynamics | S1 | T1 | P1 |
| J4.9 | **ADM1 anaerobic digestion** | Batstone et al., IWA (2002) | Digester dynamics, biogas | S1+S10 | T1 | P1 |
| J4.10 | Biofilm reactor model | Wanner & Gujer (1986); Rittmann | Multispecies biofilm gradients | S4+S12 | T1 | P2 |
| J4.11 | Fed-batch fermentation optimisation | — | Feeding policy | S1 | T1 | P2 |
| J4.12 | Flux balance analysis | Palsson | Metabolic yields (constraint-based, not PDE) | — | T1 | P3 |

**Note.** J4.7 is literally the Thiele problem (B1.1) with Michaelis–Menten
kinetics. Cross-linking it makes the gallery's structure-code axis pay off
visibly: a bioengineer arriving at J4.7 discovers B1.1 and the whole `S3` family.

## J5 — Polymerisation reactors

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| J5.1 | Method of moments (free radical) | Ray, *Rev. Macromol. Chem.* 8 (1972) | Mn, Mw, PDI | S1 | T0 | **P1** |
| J5.2 | Flory–Schulz / most probable distribution | — | MWD baseline | S1 | T0 | P1 |
| J5.3 | Gel (Trommsdorff–Norrish) effect | Chiu, Carratt & Soong (1983) | Autoacceleration, runaway | S1 | T0 | P1 |
| J5.4 | Multigrain model | Floyd, Choi, Taylor & Ray (1986) | Ziegler–Natta particle growth | S8 | T1 | P2 |
| J5.5 | Polymer particle growth / broken-grain | — | Fragmentation and diffusion limitation | S8+S12 | T2 | P3 |
| J5.6 | Smith–Ewart emulsion polymerisation | (1948) | Particle number, rate | S1+S11 | T0 | P2 |
| J5.7 | Full MWD via PBE / Galerkin | Wulkow | Chain length distribution | S11 | T2 | P3 |
| J5.8 | Multi-site Ziegler–Natta deconvolution | Soares & Hamielec | MWD from site heterogeneity | S1 | T1 | P2 |

## J6 — Emerging / cross-cutting

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| J6.1 | Photoreactor LVRPA / six-flux | Cassano & Alfano, *Catal. Today* 58 (2000); Brucato & Rizzuti | Local photon absorption rate | S6 | T1 | P2 |
| J6.2 | Photocatalytic degradation kinetics | Langmuir–Hinshelwood + LVRPA | Rate vs irradiance | S3 | T2 | P2 |
| J6.3 | Non-thermal plasma reactor (DBD) | Bogaerts et al. | Conversion, energy efficiency | S1 | T3 | P3 |
| J6.4 | Microwave / induction-heated reactor | — | Non-uniform volumetric heating | S6 | T3 | P3 |
| J6.5 | Sonochemical reactor | — | Cavitation-driven rates | S1 | T3 | P3 |
| J6.6 | Digital-twin / data-driven hybrid | Physics-informed NN literature | Hybrid closure fitting | — | T3 | P3 |

---

## Coverage summary

Counts below are generated from the tables in these three files, not estimated.

| Section | Entries | P1 |
|---|---|---|
| A — transport closures | 43 | 30 |
| B — particle & non-catalytic | 29 | 17 |
| C — kinetics | 26 | 15 |
| D — fixed bed | 24 | 14 |
| E — fluidised bed | 18 | 7 |
| F — gas–liquid & slurry | 24 | 18 |
| G — trickle bed | 11 | 7 |
| H — membrane | 14 | 9 |
| I — structured/intensified | 13 | 6 |
| J — adjacent operations | 64 | 31 |
| **Total** | **266** | **154** |

266 catalogued models, of which 154 are marked P1 — but P1 means "good candidate
for an early page", not "must build 154 pages". A sustainable pace of one page
per week reaches a credible gallery (~50 pages spanning every section) in a year;
the full catalog is a multi-year community effort, which is the argument for
making outside contribution easy from day one.

Note that P1 density is highest in A (transport closures) because correlation
pages are cheap — load data, plot correlation, show deviation. Those pages are
individually less impressive than a coupled reactor model but they are what
makes the gallery *useful*, since a modeller picking closures needs exactly that
comparison and rarely finds it in one place.
