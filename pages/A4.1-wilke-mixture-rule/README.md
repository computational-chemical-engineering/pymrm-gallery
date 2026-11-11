# A4.1: Fick's law with Wilke's mixture rule for the effective diffusivity

The weighted-harmonic mixture rule

```
D'_A = (1 - y_A) / ( y_B/D_AB + y_C/D_AC + y_D/D_AD + ... )
```

printed in the abstract of Wilke (1950), reprinted forty-three years later as
eq. (6.1.14) of Taylor and Krishna with `(Wilke, 1950)` inline. Every reactor
model that carries one diffusion coefficient per species in a multicomponent gas
is using some version of it.

**The rule is exact when the other species are genuinely stagnant.** That is not
an approximation but an algebraic identity, so a page that reports agreement
there reports arithmetic. This page therefore goes where it is not exact, and
the paper itself supplies the place: eight ternary test examples, each solved
three ways, with the exact Stefan-Maxwell solution printed beside the rule as
Table 1.

## Headline

Against the exact solution of Wilke's own eight problems, computed here:

| | max error | median error |
|---|---|---|
| **Method 1**, the rule at the arithmetic-mean composition | 86.6 % <!-- claim: wilke_rule_max_err_vs_exact --> | 21.6 % <!-- claim: wilke_rule_median_err_vs_exact --> |
| **Method 2**, the *same* rule at a flux-weighted film composition | 12.1 % <!-- claim: method2_max_err_vs_exact --> | 0.18 % <!-- claim: method2_median_err_vs_exact --> |
| **null**, weighted *arithmetic* mean instead of harmonic | 212.7 % <!-- claim: arithmetic_null_max_err_vs_exact --> | 22.2 % <!-- claim: arithmetic_null_median_err_vs_exact --> |

Three readings of that table matter more than the numbers themselves.

- **The split is by flux share.** On the species whose flux is at least half the
  larger of the two the rule is 7.7 % <!-- claim: wilke_rule_median_err_majority_flux -->
  off at the median. In the examples whose dominant species has `|φ - 1|` below
  `0.05` (IV to VIII) it is never worse than 2.31 % <!-- claim: wilke_rule_max_err_dominant_phi_near_1 -->
  on that species. On the minority species it is 24.2 % <!-- claim: wilke_rule_min_err_minority_flux -->
  to 86.6 % <!-- claim: wilke_rule_max_err_minority_flux --> off, and the sign is
  not even consistent: +86.6 % <!-- claim: wilke_rule_err_IV_B --> on Example IV,
  -41.3 % <!-- claim: wilke_rule_err_VII_B --> on Example VII.
- **The harmonic mean is not what limits it.** Swapping it for the arithmetic
  mean that J. H. Arnold attributes to Hougen and Watson *in the paper's own
  printed Discussion* moves the median from 21.6 % <!-- claim: wilke_rule_median_err_vs_exact -->
  to 22.2 % <!-- claim: arithmetic_null_median_err_vs_exact --> while the worst case
  goes from 86.6 % <!-- claim: wilke_rule_max_err_vs_exact --> to 212.7 % <!-- claim: arithmetic_null_max_err_vs_exact -->.
  Keeping the harmonic mean and moving the **composition** takes the worst case to
  12.1 % <!-- claim: method2_max_err_vs_exact -->.
- **Wilke's own exact column is not uniformly good either.** Reproducing his
  Method-3 column against the page's exact solution, every cell but one is within
  2.31 % <!-- claim: table1_method3_repro_second_max -->, and Example VI's `N_A`
  is 10.42 % <!-- claim: table1_method3_repro_max_reldev --> out. That cell is
  the minority-flux species again, and it is his hand solution of the
  transcendental pair, not a reading error. It is reported at full weight; the
  page excludes nothing from the Method-2 or Method-3 reproduction.

## What pymrm adds

Not the answer: the exact problem is a 2 x 2 linear system with constant
coefficients and a closed-form solution, derived on the page. The finite-volume
film (`construct_grad`, `construct_div`, face compositions from
`compute_boundary_values`, `NumJac` and `pymrm.newton`) knows only
`N = -c [B]^-1 grad(y)` and `div N = 0`; on a nested ladder it converges onto the
closed form at observed orders of at least 1.99 <!-- claim: pymrm_grid_order_min -->
and at most 2.02 <!-- claim: pymrm_grid_order_max -->, with a worst flux error of
5.5e-6 <!-- claim: pymrm_vs_closed_form_reldiff_n320 --> at 320 cells. That is
the independent check on the closed-form derivation, and the formulation that
survives when the linearity does not (a reacting film, composition-dependent
diffusivities, a fourth species with its own flux). `pymrm.checks.find_roots`
finds every flux pair as a scanned root and reports how many there are: one per
problem, all 32 <!-- claim: examples_with_one_root_each_method --> of them.

## Source, and which Wilke it is

C. R. Wilke, *"Diffusional Properties of Multicomponent Gases"*, **Chemical
Engineering Progress 46**(2) 95–104, February 1950. No DOI; pre-DOI journal. The
file on disk is a microfilm scan whose text layer is among the two worst in this
library (182 characters extract from journal page 97 and the by-line comes back
as `C.R. WI E`), so **no character on this page comes from the text layer**.
`pdfimages -list` reports JPEG grayscale at 400 × 400 ppi native; everything was
rendered at that resolution, cropped and re-read at digit scale. Table 1 is
printed rotated 90° and was rotated losslessly before reading.

Identity confirmed from the document itself (title page, by-line, affiliation,
running feet `Vol. 46, No. 2` / `CHEMICAL ENGINEERING PROGRESS` /
`February, 1950`, pages 95–104, the 30-item Literature Cited and the signed
Discussion) and independently by Taylor and Krishna's reference list.

**This is not the other Wilkes.** It is not the Buddenberg–Wilke *viscosity*
mixture rule (*Ind. Eng. Chem.* **41**, 1345, 1949; *J. Chem. Phys.* **18**, 517,
1950) that Bird/Stewart/Lightfoot and Chapman & Cowling carry when they cite
"Wilke"; that rule appears in *this* paper too, as eqs. (39)–(43) and Table 2,
and is scoped out as a different catalogue case. It is not Wilke & Lee's
binary-diffusivity **estimation** correlation (*Ind. Eng. Chem.* **47**, 1253,
1955), a different paper, journal, year and equation. It is not Wilke–Chang.

**Second witness.** Taylor & Krishna (1993), *Multicomponent Mass Transfer*,
Wiley, read at its native 600 ppi: section 6.1.3 (printed p. 126) prints the
rule as eq. (6.1.14) with the attribution inline and states its condition;
section `8.6` (p. 204) records the pedigree; Example 8.5.1 (pp. 201–203),
"based on experiments conducted by Fairbanks and Wilke (1950) with a view to
assessing the validity of Wilke's effective diffusivity formula", is reproduced
in full. Note for `docs/papers-inventory.yaml`: the pedigree sentence is in
section `8.6`, not `8.7`, and the exact/linearised comparison is sub-section
`8.6.1` (p. 208) together with 8.5.5 (p. 203).

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext; edit `index.py`).
- `model.py`: the film equations, the three exact routes and Methods 1 and 2.
- `test_model.py`: fast checks of `model.py`, against the committed claims.
- `data/wilke-1950-table1.csv`: every cell of Table 1 (tier 6).
- `data/wilke-1950-printed-intermediates.csv`: the sample calculation, the
  Example-VI film-pressure-factor paragraph, the Example-IV diffusivity
  averages, and the paper's own accuracy claim (tier 6).
- `data/taylor-krishna-1993-example-8-5-1.csv`: the second witness's worked
  example (tier 6).

**No other page's dataset is loaded**, so the cross-page reconciliation rule does
not apply. No figure is digitised and no page image exists anywhere in this
directory.

## Three printed defects, proved from the paper's own numbers

1. **Two missing minus signs.** Table 1's Method-1 `N_B` for Examples I and II
   is printed without a minus, while the structurally identical Example III
   prints one. The magnitudes agree with Method 1 to within 0.06 % <!-- claim: method1_sign_defect_magnitude_dev -->
   , so the signs are lost, not
   the numbers.
2. **Example VI, Method 1, `N_A`.** The printed cell implies a film pressure
   factor of 379.9 <!-- claim: pf_A_from_printed_cell --> mm for A. The paper's
   own text on journal p. 99 states that quantity: *"In Example VI the film
   pressure factor for the diffusion of A as calculated under Method 1 is
   288 mm."* Solving Method 1 gives 287.4 <!-- claim: pf_A_exampleVI_method1 --> mm,
   1.2 <!-- claim: pf_A_exampleVI_half_units --> half-units of the printed last
   digit away, and the same paragraph's other three numbers round to their printed
   values. The paragraph is internally consistent; the table cell is the outlier.
3. **Eq. (22)'s first bracket.** Printed as `[1/D_AB - 1/D_AC]`; every other
   subscript in (22) is (21)'s under A and B exchanged with `D_AC -> D_BC`.
   Shown symbolically: (21) reproduces the Maxwell-Stefan balance identically,
   (22) does so only with `D_BC`, and (21)+(22) integrates to the paper's own
   eq. (23) only with `D_BC`. Eq. (23) is not in doubt: 8 <!-- claim: eq23_pairs_inside_1pct -->
   of the ten printed Method-3 pairs satisfy it inside `1 %`.

All three are **reported, never repaired**: the CSV stores every cell as
printed. So are two printed *typos*, quoted with `[sic]` rather than silently
corrected: journal p. 100 prints *"a maximum **derivation** of less than 10%"*
(crisp at native resolution, not a degraded "deviation"), and Taylor & Krishna's
printed p. 126 reads *"When species i diffuses through **a** n − 1 stagnant
gases"*.

One reference number in the paper is also wrong and is flagged rather than
silently mapped: Wilke cites "(25)" for both *"Sherwood (25) presents
differential equations"* and *"Gilliland (25) has integrated (21) and (22)"*,
but his ref. 25 is *Sherwood, T. K., Ibid., p. 11 (1937)*, the same book as
ref. 24, while Gilliland & Sherwood (1934) is his ref. **12**.

Two further cells are left unrepaired because the constraint does not decide
them. Example I's third Method-3 set, `(4.14, −3.93)`, misses eq. (23) by
-13.2 % <!-- claim: eq23_exI_set3_resid --> against a 0.73 % <!-- claim: eq23_exI_set3_round_bound -->
rounding bound; holding `N_A` the equation wants `N_B` = -3.723 <!-- claim: eq23_exI_set3_NB_required -->
and holding `N_B` it wants `N_A` = 4.302 <!-- claim: eq23_exI_set3_NA_required -->,
and nothing chooses between them. Example III's pair misses by +9.4 % <!-- claim: eq23_exIII_resid -->
against a 0.38 % <!-- claim: eq23_exIII_round_bound --> bound and no single-digit
repair reconciles it.

## The paper's own accuracy claim for Method 2

*"a maximum derivation [sic] of less than 10%"* holds on Wilke's own two printed
columns, 9.07 % <!-- claim: method2_max_dev_on_printed_columns -->, and still
holds when Example VII's alternative Method-3 reading is taken, 9.95 % <!-- claim: method2_max_dev_on_printed_alt_reading -->
. Against the exact
solution it fails: 12.1 % <!-- claim: method2_max_err_vs_exact -->.

## Two ambiguous glyphs, settled by arithmetic

Table 1's header multiplier and the exponent in the p. 101 line
`N_A = −0.105 × 10^?` are both **illegible at native resolution**, and the pixel
shape is not used. Both are fixed to `1e5` by arithmetic: Method 1 for Example I
recomputes onto the printed entry only with that multiplier, and the p. 101
sample calculation's own chain closes on the printed `−0.105` only with it.

## Reproduction of what the paper prints

The p. 101 sample calculation: all 9 <!-- claim: sample_calc_within_print --> printed
steps round to their printed values (largest relative deviation 0.29 % <!-- claim: sample_calc_max_reldev -->
). Table 1 against the printed cells:
Method 1 to 1.92 % <!-- claim: table1_method1_repro_max_reldev --> over the
cells that are not one of the three printed defects, Method 2 to 5.80 % <!-- claim: table1_method2_repro_max_reldev -->
and Method 3 to 10.42 % <!-- claim: table1_method3_repro_max_reldev -->
over all sixteen cells. Few cells round to the printed entry: 3 <!-- claim: table1_m1_cells_within_print --> of the thirteen Method-1 cells, 1 <!-- claim: table1_m2_cells_within_print --> of sixteen for Method 2 and 4 <!-- claim: table1_m3_cells_within_print --> of sixteen for Method 3, at the page's chosen R and T, which the paper does not print. The Example-IV diffusivity averages: eq. (25) gives 2.2222 <!-- claim: eq25_DppB_exampleIV -->
against a printed `2.22`, and the path integral
2.1571 <!-- claim: path_average_DppB_exampleIV --> against a printed `2.15` from a
graphical evaluation, 1.42 <!-- claim: path_average_DppB_exampleIV_half_units -->
half-units of the printed digit away. Taylor and Krishna's five intermediates
all round to their printed values, largest deviation 9.8e-5 <!-- claim: taylor_krishna_851_max_reldev -->
; their flux formula and Wilke's
eq. (14) at φ = 1 are the same expression, so that agreement is an identity.
The live number beside it: the rule frozen at one composition is off by -0.080 % <!-- claim: tk_frozen_composition_err -->
at the printed interface and by 31.0 % <!-- claim: stagnant_frozen_composition_worst_err -->
at the widest driving force,
with every other species stagnant.

## What the page cannot conclude

**Nothing measured.** Every dataset is tier 6: the authors' own computed
numbers. The paper's only experimental content is Figure 1 (from Fairbanks's
1948 M.S. thesis), which is a figure, is scoped out, and is not digitised. The
paper's own Comment section says validity of the Stefan-Maxwell equations *"has
not been established experimentally in multicomponent systems for other than the
stagnant gas case cited here"*, so the exact solution is used here as the
**reference**, and if those equations are wrong for these mixtures then so is
every error bar on the page. This is a model-versus-model reproduction, not a
validation.

**Nothing about these gases.** Table 1's own footnote: the diffusion
coefficients "are estimated and do not correspond to true values".

**Nothing beyond three components.** Every example is ternary with exactly one
stagnant species. Wilke's own extrapolation of Method 2 to more components is
quoted and is not tested.

## Related pages

`A4.2` reaches a related conclusion about scalar effective diffusivities from
Krishna & Wesselingh (1997) on the Duncan–Toor two-bulb measurements. It shares
no dataset, no solver and no number with this page, deliberately, so the two
are independent evidence rather than one argument told twice. `A4.3`, `A4.4`,
`A4.9` and `H1.9` are the rest of the multicomponent-diffusion cluster.

## Corrections

- 2026-11-10: said of the three routes to the exact answer that "the third
  shares nothing at all with either." The pymrm finite-volume film builds its
  face mobility from the same matrix `[B]` that defines the analytic routes'
  system; it shares that matrix and nothing downstream of it (integration
  against discretisation).
- 2026-11-10: printed "Every solve on this page is a bracketed" "root find on a fixed deterministic grid". The pymrm finite-volume film is a Newton solve from a linear start; the flux pairs of the exact problem, of Methods 1 and 2 and of Taylor and Krishna's example are scanned roots.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: the three datasets are committed CSVs in `data/`.
