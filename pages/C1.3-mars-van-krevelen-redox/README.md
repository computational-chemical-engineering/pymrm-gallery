# C1.3: Mars-van Krevelen redox kinetics

The 1954 origin of redox (Mars-van Krevelen) kinetics, read from a scan of the
paper itself. Every numeric on this page comes off cropped renders at the scan's
**native 300 ppi**; the text layer was used for nothing.

- **Structures:** `S1` (1D steady plug flow)
- **Reference:** Mars & van Krevelen (1954), Chem. Eng. Sci. Special Supplement 3, 41-59, doi:10.1016/S0009-2509(54)80005-4
- **Runtime:** about 45 s
- **Files:** `index.py` (the page, paired with `index.ipynb`), `model.py`
  (eq. (8) and `alpha_eq`, the eq. (14) fits, the three discrimination laws,
  `OptimumCurve`, `split_quad` and the pymrm `OptimumProfileConverter`),
  `test_model.py`

## What the page establishes

1. **What the paper prints, and what it does not.** Every number in the paper
   is in the SO₂ application. For the aromatic oxidations the mechanism is
   named for, the paper prints *figures only*: no rate table and no numeric
   k₁, k₂ or β for any aromatic at any temperature. The authors also disclaim
   mechanistic proof in the printed discussion ("agreement between a kinetic
   formula and experimental data is not a proof of the correctness of the
   assumed reaction mechanism", their reply to J. M. Smith).
2. **The worked design reproduces.** The paper's two-stage SO₂ converter
   catalyst split along the optimum temperature curve, computed from printed
   inputs alone: converter 1 takes **23.6 <!-- claim: split_converter1_pct_quad --> %** of
   the catalyst against the printed 23 : 77 (plant: 27 : 73), by adaptive quadrature and a pymrm
   plug-flow solve 0.08 <!-- claim: split_routes_gap_points --> split-points apart
   on 2000 cells. Nothing fitted in the chain. Over the six Table 4 feeds that can
   reach the converter-2 target at all the split runs from 22.47 <!-- claim: feed_split_min_pct -->
   to 23.73 <!-- claim: feed_split_max_pct -->; the seventh runs out of oxygen first.
3. **The printed intermediates hold.** E_app = 8.58 <!-- claim: Eapp_kcal_per_mol -->
   kcal/mol against "about 9"; van 't Hoff slope of eq. (8) =
   22.6 <!-- claim: vantHoff_dH_kcal_per_mol_SO3 --> against the printed 23 kcal/mol;
   both k′ tables reproduced (held-out 475-550 °C blocks:
   3.80 % <!-- claim: neumann_test_mean_dev --> mean).
4. **The rate law earns its shape.** Equal-parameter discrimination on the
   reprinted Neumann data, as the mean of the nine per-block rms: MvK 1.01 <!-- claim: discrimination_rms_mvk -->
   conversion points against 2.47 <!-- claim: discrimination_rms_boreskov -->
   (Boreskov's `(p_SO2/p_SO3)^0.8` law) and 2.85 <!-- claim: discrimination_rms_first_order -->
   (first order to equilibrium, MvK without its retardation term).
5. **Four printed defects**, all proved from the paper's own neighbouring
   equations or columns: eq. (16)'s sign-inverted RHS; the α_eq-power
   mismatch between (13a) and (14); the 425 °C k′-average slip; Table 4's
   compositions identified as pyrite-roaster gas (ν = 1.386 <!-- claim: table4_nu_bestfit -->
   fitted against 11/8).
6. **The one bracketed k′ is not explained by eq. (14).** Its row sits at
   α′ = 0.968 <!-- claim: bracketed_row_alpha_prime -->, but
   5 <!-- claim: table2_rows_more_sensitive_than_bracketed --> unbracketed rows sit
   closer to equilibrium, where eq. (14) amplifies a conversion error more.

## Fit and test

- **Fitted**: Neumann scale C and SO₂ fraction (on the 325-465 °C blocks only);
  both Küster parameters (fit only, conditions unprinted); one k per block in
  the discrimination; one stoichiometric ν on Table 4.
- **Test**: Neumann 475-550 °C blocks (nothing refitted); E_app against "about 9";
  van 't Hoff against 23; the catalyst split against 23 : 77 (no fitted quantity
  anywhere in its chain).

## Data

Four CSVs with `columns:` sidecars, all transcribed at native resolution.
Tables 1-2 are experimental but **not the authors' own**: Küster (1904) and
Neumann (1928) reprinted with the authors' computed k′ columns;
`origin_not_consulted` in the sidecars. The paper states neither Neumann's gas
composition nor his velocity units; the sidecars and the page carry both
limitations. **Nothing digitised**: no figure point was read; Fig. 16's optimum
curve is recomputed from the model, not traced.

## Corrections

- 2026-11-03: said "2.85 for first-order-to-equilibrium and 3.80 for Boreskov's
  0.8-power law" and "the latter failing worst in the 525-550 C blocks". The
  Boreskov fits of those two blocks ended on the upper end of a bounded search
  for k. At the minimum of the same sum of squares Boreskov's law gives
  2.47 <!-- claim: discrimination_rms_boreskov --> points rms (mean of the per-block
  rms) and fails worst in the 325 °C block.
- 2026-11-03: said the bracketed k′ was "explained quantitatively" and that the
  authors' decision to "bracket exactly this row and no other is what eq. (14)
  predicts". Eq. (14) moves the bracketed row's k′ by
  13.7 <!-- claim: bracketed_row_sensitivity_pct --> % for a half-point conversion
  error, while 5 <!-- claim: table2_rows_more_sensitive_than_bracketed --> unbracketed
  rows move more, up to 62.8 <!-- claim: table2_max_sensitivity_pct --> %. Its
  distance from its block does not single it out either:
  2 <!-- claim: table2_rows_further_from_block_than_bracketed --> unbracketed rows lie
  further from theirs. The page cannot say why it is bracketed.
- 2026-11-03: said "the result shown to move by at most 0.4 points across the whole
  printed feed-composition range". Only two of the seven feeds were computed. The
  feed richest in SO₂ cannot reach the converter-2 target; over the other six the split moves
  by 1.26 <!-- claim: feed_split_range_points --> split-points.
- 2026-11-03: said "the 81->97 % leg wants 445-490 C". The optimum temperature at
  the end of that leg is 401 <!-- claim: T_opt_097_C --> °C; at its start it is
  492 <!-- claim: T_opt_081_C --> °C.
- 2026-11-03: the k′-average checksum said "four more differ only by
  rounding-chain residues". The remaining blocks do: 500 °C is an exact rounding tie and
  2 <!-- claim: kavg_blocks_within_one_unit --> (465 and 550 °C) differ by at most
  one unit of the last digit.
- 2026-11-03: the discrimination said it fits "Neumann's 33 conversions"; it
  fits the 32 <!-- claim: discrimination_rows --> unbracketed ones.
