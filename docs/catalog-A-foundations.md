# Catalog A — Foundations: transport closures, particle models, kinetics

Sections A, B, C of the taxonomy. These are the building blocks that every
reactor page in [catalog B](catalog-B-reactors.md) imports. Many are one-line
correlations rather than differential models — they still deserve pages, because
a gallery page that *shows the correlation against the data it was fitted to* is
exactly what a modeller needs when deciding whether it applies to their case.

Legend: **S** = structure code, **T** = tier, **P** = build priority
(see [taxonomy.md](taxonomy.md)).

---

## A1 — Momentum / pressure drop in porous media

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| A1.1 | Ergun equation | Ergun, *Chem. Eng. Prog.* 48:89 (1952) | Δp in packed beds, viscous + inertial | S3 | T0 | P1 |
| A1.2 | Kozeny–Carman | Carman (1937) | Δp, creeping flow, basis of Ergun viscous term | S3 | T0 | P1 |
| A1.3 | Darcy–Forchheimer | — | Δp in foams, monoliths, general porous media | S3 | T0 | P1 |
| A1.4 | Eisfeld–Schnitzlein wall correction | *Chem. Eng. Sci.* 56 (2001) | Δp at low tube-to-particle ratio | S3 | T2 | P2 |
| A1.5 | Richardson–Zaki | *Trans. IChemE* 32 (1954) | Bed expansion / hindered settling | S1 | T0 | P1 |
| A1.6 | Wen–Yu minimum fluidisation | *Chem. Eng. Prog. Symp.* (1966) † | u_mf from particle properties | S3 ‡ | T0 | P1 |
| A1.7 | Geldart classification | *Powder Technol.* 7:285 (1973) | Fluidisation regime from ρ, d_p | — | T0 | P1 |
| A1.8 | Gidaspow / Syamlal–O'Brien / Wen–Yu drag | Gidaspow (1994); Syamlal & O'Brien (1989) § | Gas–solid interphase drag closures | S1 ‡ | T1 | P2 |
| A1.9 | EMMS drag | Li & Kwauk (1994) | Meso-scale-corrected drag, heterogeneous flow | S1 | T3 | P3 |

† `A1.6`'s page is built from Wen & Yu's *A.I.Ch.E. Journal* **12**(3) 610–612
communication, not from the Symposium Series paper listed here. The Symposium
paper is reference 23 of the communication and carries the derivation of eq. (1);
it is not on disk, and the page reconstructs the derivation instead, validating
the reconstruction by recovering both printed constants to under 0.1 %. It is
recorded as `origin_not_consulted` in `models.yaml`.

‡ `A1.6` and `A1.8` were catalogued `S1` and are built as `S3`, following `A1.1`
and `A1.7`: they are algebraic closures with no operator, grid or solve
(`pymrm_api: []`). Corrected 2026-08-02 (`A1.6`) and 2026-08-03 (`A1.8`). Note
that `check_metadata.py` compares `meta.yaml` against `models.yaml` only and
would not have caught the mismatch with this file.

§ `A1.8`'s page is built from the **MFIX Documentation Theory Guide**, Syamlal,
Rogers & O'Brien (1993), DOE/METC-94/1004,
[doi:10.2172/10145548](https://doi.org/10.2172/10145548), and the catalogue's two
citations are both wrong for what it delivers. The report's own reference list
gives the origin of its eq. (11) as Syamlal & O'Brien (1987), *"A Generalized
Drag Correlation for Multiparticle Systems," **Unpublished report***, so there is
**no origin paper to acquire** and the 1993 DOE report is the citable *published*
source. "Syamlal & O'Brien (1989)" is a different item in the same list — a
bubble-simulation paper, AIChE Symp. Ser. No. 270, **85**, 22–31. And the report
contains **neither** of the other two named drag closures: no Wen–Yu (the strings
do not occur) and no Gidaspow drag law or blend rule. The page therefore builds
Syamlal–O'Brien against Ergun (via `A1.1`) and Richardson & Zaki (via `A1.5`'s
source), and names the two documents that would complete the three-way
comparison — Gidaspow's 1994 monograph and the Wen & Yu Symposium Series paper.

**Gallery angle for A1.** One page overlaying Ergun, Kozeny–Carman,
Forchheimer, and the Eisfeld wall correction against a single Δp dataset, with
the wall-effect regime highlighted, is far more useful than four pages. Ergun's
original data are **not** tabulated: the 1952 paper contains no tables at all,
and no text layer either (`pdftotext` returns four bytes for four pages), so the
640 experiments survive only as scatter in its figures. The dataset on the built
page was digitised from Figures 5 and 7 off 600 dpi renders, as was every
constant quoted from the paper.

## A2 — Dispersion and residence time

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| A2.1 | Danckwerts boundary conditions | Danckwerts, *Chem. Eng. Sci.* 2:1 (1953) | Inlet/outlet closure for axial dispersion | S4 | T0 | P1 |
| A2.2 | Wehner–Wilhelm | *Chem. Eng. Sci.* 6 (1956) | Correct open/closed vessel BCs | S4 | T0 | P1 |
| A2.3 | Taylor–Aris dispersion | Taylor (1953); Aris (1956) | D_eff from laminar velocity profile | S6→S4 | T0 | P1 |
| A2.4 | Tanks-in-series | Levenspiel | RTD from N stirred tanks | S1 | T0 | P1 |
| A2.5 | Edwards–Richardson axial dispersion | *Chem. Eng. Sci.* 23 (1968) | Pe_ax in packed beds vs Re | S4 | T1 | P1 |
| A2.6 | Gunn dispersion correlations | *Chem. Eng. Sci.* 42 (1987) | Axial + radial Pe, wide Re range | S6 | T1 | P1 |
| A2.7 | Westerterp wave model | Westerterp, Dil'man & Kronberg, *AIChE J*/*Chem. Eng. Sci.* (1995) | Hyperbolic alternative to parabolic dispersion; no upstream signalling | S4 | T2 | P2 |
| A2.8 | Zwietering segregation | *Chem. Eng. Sci.* 11 (1959) | Micromixing bounds on conversion | S1 | T0 | P2 |
| A2.9 | Baldyga–Bourne engulfment | *Chem. Eng. Sci.* (1989) | Micromixing-limited selectivity | S1 | T1 | P2 |
| A2.10 | Compartment models from CFD | Bezzo & Macchietto (2004) | Reduced-order mixing from CFD fields | S7 | T2 | P3 |

**Gallery angle for A2.** A2.3 is the single best "pymrm improves on the
original" candidate in the whole catalog: solve the actual 2D Graetz problem
(`S6`) and show where the Taylor–Aris 1D effective coefficient breaks down at
short times / high Pe. The teacher solutions already contain a Taylor dispersion
exercise to start from.

## A3 — Interphase mass and heat transfer

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| A3.1 | Whitman two-film | *Chem. Metall. Eng.* 29 (1923) | Series film resistances | S3 | T0 | P1 |
| A3.2 | Higbie penetration | *Trans. AIChE* 31 (1935) | k_L ∝ √(D/t_c) | S4 | T0 | P1 |
| A3.3 | Danckwerts surface renewal | *Ind. Eng. Chem.* 43 (1951) | k_L with age distribution | S4 | T0 | P1 |
| A3.4 | Wakao–Funazkri | *Chem. Eng. Sci.* 33:1375 (1978) | Sh/Nu particle–fluid in packed beds | S3 | T1 | P1 |
| A3.5 | Ranz–Marshall / Frössling | (1952) / (1938) | Sh, Nu for single sphere | S3 | T0 | P1 |
| A3.6 | Calderbank–Moo-Young | *Chem. Eng. Sci.* 16 (1961) | k_L for bubbles/drops, small vs large | S3 | T0 | P1 |
| A3.7 | Van 't Riet | *IEC Proc. Des. Dev.* 18 (1979) | k_La in stirred tanks from P/V, u_g | — | T1 | P1 |
| A3.8 | Onda correlations | *J. Chem. Eng. Japan* 1 (1968) | k_L a, k_G a, wetted area in packings | — | T1 | P1 |
| A3.9 | Billet–Schultes | *Chem. Eng. Res. Des.* (1993, 1999) | Packing hydraulics + mass transfer | — | T1 | P2 |
| A3.10 | Rocha–Bravo–Fair | *IECR* 32 (1993), 35 (1996) | Structured packing hydraulics/MT | — | T1 | P2 |
| A3.11 | Dixon–Cresswell heat transfer params | *AIChE J* 25 (1979) | λ_er, h_w, and their interdependence | S6 | T1 | P2 |
| A3.12 | Yagi–Kunii effective conductivity | *AIChE J* 3 (1957) | λ_eff of packed beds, static + dynamic | S6 | T0 | P1 |
| A3.13 | Zehner–Bauer–Schlünder | (1970s) | Stagnant bed conductivity incl. radiation | S6 | T1 | P2 |
| A3.14 | Martin–Nilles wall heat transfer | *Chem. Ing. Tech.* 65 (1993) | h_w in tubular packed beds | S6 | T1 | P2 |
| A3.15 | Graetz–Nusselt problem | Graetz (1883) | Developing thermal/concentration BL in a tube | S6 | T0 | P1 |

## A4 — Multicomponent diffusion

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| A4.1 | Fick + Wilke mixture rule | Wilke, *Chem. Eng. Prog.* 46 (1950) | Effective D_i in a mixture | S3 | T0 | P1 |
| A4.2 | Maxwell–Stefan | Krishna & Wesselingh, *Chem. Eng. Sci.* 52:861 (1997) | Coupled fluxes, osmotic/reverse diffusion | S9 | T0 | P1 |
| A4.3 | Dusty gas model | Mason & Malinauskas (1983) | Combined bulk + Knudsen + viscous in pores | S9 | T1 | P1 |
| A4.4 | Knudsen / Bosanquet | — | Pore-size-limited diffusion | S3 | T0 | P1 |
| A4.5 | Fuller–Schettler–Giddings | *Ind. Eng. Chem.* 58 (1966) | Binary gas diffusivity estimation | — | T0 | P1 |
| A4.6 | Chapman–Enskog | — | Kinetic-theory binary diffusivity | — | T0 | P1 |
| A4.7 | Krishna zeolite/micropore MS | Krishna & van den Broeke (1995); Krishna & Baur (2003) | Loading-dependent MS diffusion in micropores | S9 | T1 | P2 |
| A4.8 | Stefan tube / Arnold diffusion | — | Classic MS validation experiment | S9 | T0 | P1 |
| A4.9 | Duncan–Toor three-bulb experiment | *AIChE J* 8 (1962) | Osmotic + reverse diffusion, MS proof case | S9 | T0 | **P1** |

**Gallery angle for A4.** A4.9 is the highest-value early page in the entire
gallery. The Duncan–Toor experiment is small, the data are published as a table,
Fick's law provably fails on it, and Maxwell–Stefan provably succeeds. pymrm
already has the machinery (`S9`) and the teacher book has a ternary-diffusion
exercise. It is the perfect "why this library exists" demonstration.

---

## B1 — Intraparticle diffusion–reaction

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| B1.1 | Thiele modulus / effectiveness factor | Thiele, *Ind. Eng. Chem.* 31:916 (1939); Damköhler; Zeldovich | η(φ), slab/cyl/sphere | S3 | T0 | **P1** |
| B1.2 | Aris generalised shape modulus | *Chem. Eng. Sci.* 6 (1957) | Shape-independent η via V/A | S3/S13 | T0 | P1 |
| B1.3 | Bischoff generalised modulus | *AIChE J* 11 (1965) | η for arbitrary kinetics | S3 | T0 | P1 |
| B1.4 | Weisz–Prater criterion | *Adv. Catal.* 6 (1954) | Observable test for diffusion limitation | S3 | T0 | P1 |
| B1.5 | Weisz–Hicks non-isothermal | *Chem. Eng. Sci.* 17:265 (1962) | η > 1, multiplicity, β and γ | S3 | T0 | **P1** |
| B1.6 | Prater relation | *Chem. Eng. Sci.* 8 (1958) | ΔT_max inside a pellet | S3 | T0 | P1 |
| B1.7 | Mears criteria | *IEC Proc. Des. Dev.* 10 (1971) | External/internal gradient screening | — | T0 | P1 |
| B1.8 | Wheeler pore model | *Adv. Catal.* 3 (1951) | Effective diffusivity from pore structure | S3 | T0 | P2 |
| B1.9 | Wakao–Smith random pore | *Chem. Eng. Sci.* 17 (1962) | D_eff for bidisperse pore systems | S3 | T1 | P2 |
| B1.10 | Feng–Stewart pore network | *IEC Fundam.* 12 (1973) | Multicomponent pore transport | S9 | T2 | P3 |
| B1.11 | Multicomponent pellet with MS/DGM | Ethanol dehydrogenation, SMR pellets | Coupled selectivity in pellets | S9+S3 | T1 | P2 |
| B1.12 | Contour-averaged non-spherical pellet | Peters (this suite's `particle_model` project) | η for arbitrary pellet shape via 1D reduction | S13 | T3 | *deferred* |
| B1.13 | Bimodal / macro–micro pellet | — | Two-scale intraparticle resistance | S8 | T2 | P2 |

**Gallery angle for B1.** B1.1 + B1.5 together are the flagship pages: the
analytical η(φ) is exact, the Weisz–Hicks multiplicity is a genuinely hard
nonlinear solve, and both are in the existing teacher solutions.

B1.12 is the group's own contribution but the manuscript is in second revision
(July 2026), so it is **deferred** under the published-work-only policy
([blueprint §9](blueprint.md#published-work-only-policy)) even though its data is
already deposited at 4TU with a DOI. Note the side effect: B1.12 was the only
planned `S13` page, so `construct_div(nu=callable)` currently has no
demonstration in the gallery. A published annular or variable-area pellet
problem would fill that gap in the meantime.

## B2 — Catalyst deactivation

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| B2.1 | Voorhies coking law | *Ind. Eng. Chem.* 37 (1945) | Coke ∝ t^n | S1 | T0 | P1 |
| B2.2 | Froment–Bischoff deactivation | *Chem. Eng. Sci.* 16 (1961), 17 (1962) | Coke-dependent activity, parallel/series | S4 | T0 | P1 |
| B2.3 | Levenspiel deactivation kinetics | (1972) | Order-based decay laws | S1 | T0 | P1 |
| B2.4 | Beeckman–Froment pore blockage | *IEC Fundam.* 18 (1979) | Percolation-type deactivation | S3 | T2 | P3 |
| B2.5 | Sintering (Ruckenstein–Pulvermacher) | *AIChE J* 19 (1973) | Metal area loss | S1 | T2 | P3 |
| B2.6 | Poisoning shell-progressive | Wheeler | Moving poison front in pellet | S12 | T1 | P2 |

## B3 — Non-catalytic gas–solid reaction

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| B3.1 | Shrinking core (Yagi–Kunii) | *5th Symp. Combust.* (1955) | Conversion with ash-layer/film/reaction control | S12 | T0 | **P1** |
| B3.2 | Grain model | Szekely & Evans, *Chem. Eng. Sci.* 25 (1970) | Structure-resolved conversion | S3+S12 | T0 | P1 |
| B3.3 | Random pore model | Bhatia & Perlmutter, *AIChE J* 26 (1980) | Rate maximum from pore overlap, ψ parameter | S3 | T0 | **P1** |
| B3.4 | Sohn law of additive reaction times | *Metall. Trans.* (1978) | Closed-form conversion across regimes | S1 | T1 | P1 |
| B3.5 | CaO carbonation deactivation | Grasa & Abanades, *IECR* 45 (2006) | Capacity decay over cycles | S1 | T1 | P1 |
| B3.6 | Char combustion (Field; Baum–Street) | (1967); (1971) | Burnout with diffusion/kinetic control | S12 | T1 | P2 |
| B3.7 | Kobayashi two-competing-rate devolatilisation | (1977) | Volatile yield vs heating rate | S1 | T1 | P2 |
| B3.8 | CPD / FG-DVC coal pyrolysis | Fletcher & Kerstein; Solomon | Network pyrolysis | S1 | T2 | P3 |
| B3.9 | Broido–Shafizadeh biomass pyrolysis | Di Blasi review, *PECS* 34 (2008) | Competing char/tar/gas paths | S1 | T1 | P1 |
| B3.10 | Oxygen-carrier reduction/oxidation | Abad & Adánez, *Chem. Eng. Sci.* (2007) | CLC particle conversion | S12 | T1 | P2 |

---

## C1 — Rate-law formalisms

| ID | Model | Canonical reference | Predicts | S | T | P |
|---|---|---|---|---|---|---|
| C1.1 | Langmuir–Hinshelwood–Hougen–Watson | Hougen & Watson (1943, 1947) | Surface-reaction-controlled rates | S1 | T0 | P1 |
| C1.2 | Eley–Rideal | — | Gas-phase attack on adsorbed species | S1 | T0 | P1 |
| C1.3 | Mars–van Krevelen | *Chem. Eng. Sci. Spec. Suppl.* 3 (1954) | Redox oxidation kinetics | S1 | T0 | P1 |
| C1.4 | Microkinetic / mean-field surface | Dumesic; Deutschmann (DETCHEM) | Elementary-step coverage-resolved rates | S1/S10 | T1 | P2 |
| C1.5 | Power-law + Arrhenius | — | Baseline empirical fit | S1 | T0 | P1 |
| C1.6 | Sabatier / volcano correlation | — | Activity vs binding energy | — | T1 | P3 |

## C2 — Named industrial kinetics (each is a page with its own data)

| ID | System | Canonical reference | Notes | T | P |
|---|---|---|---|---|---|
| C2.1 | Steam methane reforming | **Xu & Froment, *AIChE J* 35:88 (1989)** | The field standard; 3 reactions, LHHW; verified in searches as the most widely adopted SMR framework | T1 | **P1** |
| C2.2 | SMR alternative | Numaguchi & Kikuchi (1988) | The standard comparison to Xu–Froment — build as a companion page | T2 | P2 |
| C2.3 | Ammonia synthesis | Temkin & Pyzhev (1940); Dyson & Simon (1968) | Classic; Dyson–Simon is the usual engineering form | T0 | **P1** |
| C2.4 | Methanol synthesis | Graaf et al., *Chem. Eng. Sci.* 43 (1988) | Three-reaction LHHW | T1 | P1 |
| C2.5 | Methanol synthesis (alternative) | Vanden Bussche & Froment, *J. Catal.* 161:1 (1996) | CO2-based mechanism; direct comparison page with C2.4 | T1 | P1 |
| C2.6 | Water–gas shift | Moe (1962); Keiski et al. | HT and LT shift | T1 | P1 |
| C2.7 | Fischer–Tropsch | Yates & Satterfield, *Energy Fuels* 5 (1991); ASF distribution | Rate + chain-growth product distribution | T1 | P1 |
| C2.8 | Steam cracking of hydrocarbons | Sundaram & Froment, *Chem. Eng. Sci.* 32 (1977) | Radical network, molecular scheme | T1 | P2 |
| C2.9 | Hydrodesulfurisation | Broderick & Gates, *AIChE J* 27 (1981) | DBT hydrogenolysis/hydrogenation | T1 | P2 |
| C2.10 | o-Xylene → phthalic anhydride | Froment (1967); Calderbank | The classic hot-spot/runaway case | T0 | **P1** |
| C2.11 | Ethylene oxidation to EO | Borman & Westerterp, *IECR* 34 (1995) | Selectivity + hot spot (TU/e lineage) | T1 | P1 |
| C2.12 | SCR of NOx | Nova & Tronconi | Standard/fast SCR, NH3 inhibition | T1 | P2 |
| C2.13 | Three-way catalyst / CO oxidation | Voltz et al. (1973) | Self-inhibited CO oxidation; used in every converter model | T0 | P1 |
| C2.14 | Methanol-to-olefins / hydrocarbon pool | Gayubo, Bos et al. | Deactivating, autocatalytic | T2 | P2 |
| C2.15 | Oxidative coupling of methane | Stansch, Mleczko & Baerns (1997) | 10-step network | T2 | P2 |
| C2.16 | Propylene ammoxidation | — | Mars–van Krevelen industrial case | T2 | P3 |
| C2.17 | Methanation / Sabatier | — | Power-to-gas relevance | T2 | P1 |
| C2.18 | CO2 hydrogenation to methanol/DME | — | Modern, active field | T3 | P2 |
| C2.19 | Ethanol dehydrogenation | (in this suite's teacher exercises) | Already has an exercise to build from | T2 | P1 |
| C2.20 | Claus process | — | Industrial sulfur recovery | T2 | P3 |

**Gallery angle for C2.** Kinetics pages are the ones with the best data
situation: kinetic papers publish rate-vs-conversion tables and parity plots,
and a parity plot *is* the validation figure. C2.1, C2.3, C2.4/C2.5, C2.10 are
the strongest P1 candidates.
