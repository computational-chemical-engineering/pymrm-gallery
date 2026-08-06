# C1.3 — Mars–van Krevelen redox kinetics

The 1954 origin of redox (Mars–van Krevelen) kinetics, read from a real scan of
the paper itself — retiring the api-text route this case previously depended
on, which was the worst in the repository (four different glyphs for one
decimal point inside a single table row). Every numeric on this page comes off
cropped renders at the scan's **native 300 ppi**; the text layer was used for
nothing.

## What the page establishes

1. **What the paper prints, and what it does not.** Every number in the paper
   is in the SO₂ application. For the aromatic oxidations the mechanism is
   named for, the paper prints *figures only*: no rate table and no numeric
   k₁, k₂ or β for any aromatic at any temperature. The authors also disclaim
   mechanistic proof in the printed discussion ("agreement between a kinetic
   formula and experimental data is not a proof of the correctness of the
   assumed reaction mechanism" — their reply to J. M. Smith).
2. **The worked design reproduces.** The paper's two-stage SO₂ converter
   catalyst split along the optimum temperature curve, computed from printed
   inputs alone: **23.6 : 76.4** against the printed **23 : 77** (plant:
   27 : 73), by adaptive quadrature and a pymrm plug-flow solve agreeing to
   0.10 split-points. Nothing fitted in the chain; ≤ 0.4 points sensitivity
   across the whole printed feed-composition range.
3. **The printed intermediates hold.** E_app = 8.58 kcal/mol vs "about 9";
   van 't Hoff slope of eq. (8) = 22.6 vs the printed 23 kcal/mol; both k′
   tables reproduced (held-out 475–550 °C blocks: 3.80 % mean); the single
   bracketed k′ explained quantitatively (α′ = 0.968).
4. **The rate law earns its shape.** Equal-parameter discrimination on the
   reprinted Neumann data: MvK 1.0 conversion-points rms vs 2.9 (first order
   to equilibrium) and 3.8 (Boreskov 0.8-power).
5. **Four printed defects**, all proved from the paper's own neighbouring
   equations or columns: eq. (16)'s sign-inverted RHS; the α_eq-power
   mismatch between (13a) and (14); the 425 °C k′-average slip; Table 4's
   compositions identified as pyrite-roaster gas (ν = 1.386 fitted vs 11/8).

## Fit vs test

- **Fitted**: Neumann scale C and SO₂ fraction (on 325–465 °C blocks only);
  both Küster parameters (fit-only, conditions unprinted); one k per block in
  the discrimination; one stoichiometric ν on Table 4.
- **Test**: Neumann 475–550 °C blocks (nothing refitted); E_app vs "about 9";
  van 't Hoff vs 23; the catalyst split vs 23 : 77 (no fitted quantity
  anywhere in its chain).

## Data

Four CSVs with `columns:` sidecars, all transcribed at native resolution.
Tables 1–2 are experimental but **not the authors' own** — Küster (1904) and
Neumann (1928) reprinted with the authors' computed k′ columns;
`origin_not_consulted` in the sidecars. The paper states neither Neumann's gas
composition nor his velocity units; the sidecars and the page carry both
limitations.

- **Structures:** `S1` (1D steady plug flow)
- **Runtime:** ~7 s
- **Data tier:** 2 — printed tables (reprinted literature measurements)
- **Nothing digitised.** No figure point was read; Fig. 16's optimum curve is
  recomputed from the model, not traced.

## Source

Mars, P. and van Krevelen, D. W., *Oxidations carried out by means of vanadium
oxide catalysts*, Special Supplement to Chemical Engineering Science **3**
(1954) 41–59. doi:10.1016/S0009-2509(54)80005-4. Identity verified from the
document's own running head, title page and by-line (19 pp, complete, with the
printed conference discussion).
