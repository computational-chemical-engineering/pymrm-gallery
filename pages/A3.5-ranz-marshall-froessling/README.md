# A3.5: Ranz-Marshall / Frössling

`Nu = 2.0 + 0.60 Re^(1/2) Pr^(1/3)` and its mass-transfer twin, from Ranz, W. E.
& Marshall, W. R., Jr., **"Evaporation from Drops", Parts I and II**, *Chemical
Engineering Progress* **48**(3) 141–146 and **48**(4) 173–180 (1952). Pre-DOI, no
DOI exists; both parts read from 400 dpi renders of the scans on disk.

- **Structures:** `S3` (pointwise algebra, plus a 1-D spherical BVP)
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (the stagnant sphere, the Stefan problem and their closed forms),
  tested by `test_model.py`. The analysis of the paper's tables, the ladders and
  the break table stay in the notebook.
- **Runtime:** about 15 s.

## What the page shows

The correlation contains four numbers with four different statuses, and the page
never mixes them:

| constant | status |
|---|---|
| `2.0` | **theoretical**: Part I eq. (7), derived from the stagnant spherical field. Fitted to nothing. |
| `0.60` | **fitted**, to Figures 6, 7 and 9, which are Tables 1, 2, 3 and 4. |
| `1/2` on `Re` | **assumed**, from Frössling. |
| `1/3` on `Pr`/`Sc` | **assumed**, from Frössling. |

There is **no held-out set for the correlation anywhere in either part**, so every
rms quoted against those tables is a goodness of fit and is labelled `IN-SAMPLE`
where it is printed, with null baselines beside it. The intercept is the exception:
it was never fitted, so extrapolating the runs to `Re = 0` is a genuine test of it.

**Headline.** Free two-parameter extrapolations to `Re = 0`, each intercept with its
standard error from the fit's own residual scatter, assuming the straight line in `Re^(1/2)` holds down to `Re = 0`, and `t = (intercept - 2.0)/SE`:

| set | runs | intercept | vs 2.0 | t |
|---|---|---|---|---|
| heat, water (Table 1) | 18 | 2.0582 <!-- claim: t1_heat_intercept_free --> ± 0.0538 <!-- claim: t1_heat_intercept_se --> | +2.91 <!-- claim: heat_intercept_dev_pct --> % | +1.08 <!-- claim: heat_intercept_t --> |
| mass, benzene (Table 4) | 13 | 2.0523 <!-- claim: t4_benzene_intercept_free --> ± 0.1245 <!-- claim: t4_benzene_intercept_se --> | +2.61 <!-- claim: mass_intercept_dev_benzene_pct --> % | +0.42 <!-- claim: mass_intercept_t_benzene --> |
| mass, water (Table 1) | 15 | 1.7322 <!-- claim: t1_mass_intercept_free --> ± 0.0870 <!-- claim: t1_mass_intercept_se --> | -13.39 <!-- claim: mass_intercept_dev_water_pct --> % | -3.08 <!-- claim: mass_intercept_t_water --> |

Heat and benzene are inside the two-sided 5 % points of Student's t for their
degrees of freedom; water is beyond it. The still-air row says the same thing more
bluntly: it prints `N_Nu = 2.23` and `N_Nu' = 1.79` where eqs. (21) and (22) both
say 2.0, a ratio of 1.2458 <!-- claim: stillair_analogy_ratio --> where the
correlation demands 1.

**The authors saw this and named the cause.** Part II folio 173: "Data for mass
transfer show a steeper slope and a lower intercept, but the disagreement is always
less than 10 per cent at a given `N_Re`", followed by a list of contributing factors
that opens with "(1) inaccurate values of diffusivity". What the page adds is the
arithmetic: their bound is about the *curves*, where the fitted 0.60 absorbs most of
the offset, not about the intercept, where nothing was fitted; and the deficit is
species-specific: benzene, whose diffusivity came from Hirschfelder, Bird and Spotz
rather than from the authors' own Figure 5, lands where the *heat* data land.
`N_Nu' ∝ 1/D_v`, so the forced-convection intercept needs a diffusivity smaller by
1.155 <!-- claim: dv_factor_from_intercept --> and the still-air row one smaller by
1.117 <!-- claim: dv_used_over_dv_backed_out -->, the 0.204 sq.cm/sec back-out Part
II prints and calls "a low value compared with other methods of determination".

**Where the page does contradict the papers is on the conclusions.** Part I's
abstract claims the study "verified the simple expression for the Nusselt number at
zero Reynolds number"; folio 174 concludes that their calculated `D_v` "may be more
accurate than any reported in the literature". The tables decide against the
diffusivity claim: heat and benzene land on 2.0.

**A convention the papers never state for the tables.** The Reynolds column is
over-determined by the columns printed beside it. Leaving the reference temperature
`T_ref = T_a + w(T_d - T_a)` free, the residual is smallest at
`w` = 0.4996 <!-- claim: re_reference_temperature_weight --> (a root of the exact
derivative), the film temperature, with the fitted temperature exponent
-1.765 <!-- claim: re_law_T_exponent --> on kinetic theory's -1.7. It agrees with the
one reference temperature Part II prints, for the still-air run, to
0.15 <!-- claim: stillair_film_minus_stated_K --> K, and it predicts Table 3, which
prints no Reynolds column and so cannot enter the fit, to
0.82 <!-- claim: t3_heldout_rms_pct --> % rms against
1.77 <!-- claim: t3_heldout_rms_pct_free_stream --> % at the free stream.

## Agreement

| Check | Result |
|---|---|
| eq. (24)'s `0.54` and eqs. (21)'s `0.60` | imply `Pr` = 0.729 <!-- claim: pr_air_implied_by_eq24 -->: a definition, not a check |
| Table 1's coefficient through that `Pr` | `K` = 0.6020 <!-- claim: t1_heat_K_recovered -->, +0.34 <!-- claim: t1_heat_K_recovered_dev_pct --> % on the printed `0.60` |
| Reynolds exponent released | 0.4933 <!-- claim: t1_heat_Re_exponent_fitted --> ± 0.0079 <!-- claim: t1_heat_Re_exponent_se --> against the assumed `1/2` |
| eq. (24) against Table 1 (in-sample) | 1.72 <!-- claim: t1_heat_eq24_rms_pct --> % rms against a best null of 39.8 <!-- claim: t1_heat_null_best_rms_pct --> % |
| eq. (21) against Table 3 (in-sample) | 11.0 <!-- claim: t3_eq21_rms_pct --> % rms, 19.1 <!-- claim: t3_eq21_max_pct --> % worst, mean bias +6.4 <!-- claim: t3_eq21_bias_pct --> % ± 3.2 <!-- claim: t3_eq21_bias_se_pct --> % |
| Table 3 within-pair identities | inside each pair's printing allowance (worst 1.78 <!-- claim: t3_pair_abscissa_max_dev_pct --> % and 1.82 <!-- claim: t3_pair_nu_ratio_max_dev_pct --> %); eq. (21)'s within-run `N_Nu` ratio off by up to 21.9 <!-- claim: t3_pair_worst_nu_ratio_error_pct --> % |
| Table 4 run 11, settled by its ratio identity | `N_Nu'` in [9.73 <!-- claim: t4_run11_Nuprime_lower_bound -->, 10.06 <!-- claim: t4_run11_Nuprime_upper_bound -->], read 10.0; the ratios of runs 10-12 agree to 1.03 <!-- claim: t4_ratio_spread_runs_10_12_pct --> %, against 6.93 <!-- claim: t4_run11_alt_ratio_spread_pct --> % for the reading 10.6 |
| test No. 80 recomputed | 0.00042 <!-- claim: test80_printed_vs_recomputed_rel --> relative, inside half a unit of the printed last digit |
| pymrm sphere, the stagnant floor `2.0` | second order on nested geometric grids; Richardson limit below 1e-6 <!-- claim: pymrm_conduction_dev_from_2 --> from 2; a far field at ten radii shifts it by 11.1 <!-- claim: pymrm_conduction_b10_excess_pct --> % |
| pymrm Stefan flow against its closed form | 0.000962 <!-- claim: stefan_pymrm_vs_closed_rel --> at 600 cells, first order; with orders 1 and 2 removed, below 1e-8 <!-- claim: stefan_extrap_vs_closed_rel --> |
| `p_f` as the arithmetic mean | 0.41 <!-- claim: pf_arith_dev_at_y020_pct --> % out at `y = 0.20`; the log mean is exact |

## What pymrm adds

The operators of the two linear spherical problems and the flux reads. The `2.0` of
eq. (7) becomes a computation, with a far field of ten radii shifting it by
11.1 <!-- claim: pymrm_conduction_b10_excess_pct --> %, the correction a dense spray
needs and the correlation cannot express. The Stefan solve (diffusion plus the radial
convection it drives, with the total flow a root of the flux balance) settles which
average Part II's `p_f` is: the logarithmic mean, which makes the paper's
mass-transfer accounting exact for the isothermal stagnant film, and confirms the closed form
that sizes the effect Ranz said in 1993 the papers had underplayed: +1 % on `Sh` at a surface
mole fraction of 1.97 % <!-- claim: stefan_y_for_01pct --> and +10 % at
17.61 % <!-- claim: stefan_y_for_10pct -->. The fit/test split and the recovered
reference temperature, which carry the page, are NumPy and owe pymrm nothing.

## What this page does NOT do

- **It does not test the `0.60` or the two exponents.** The `0.60` was fitted to the
  very tables compared with, and every heat-transfer run is in air, so the `1/3`
  cannot be identified from these rows; Part II's own test of it is graphical, in
  Figure 12, which is not digitised.
- **It does not say which of the authors' two diffusivity statements is right.** It
  measures the size of the water deficit and the diffusivity factor that would close
  it, not its cause.
- **Nothing is digitised.** Figures 5, 12, 13 and 14 exist only as curves, so the
  free-convection exponent of eqs. (10) and (11) is never exercised and Table 6's
  calculated rates cannot be reproduced. Table 5 carries no transfer coefficient.
  Packed beds are `A3.4`.

## Corrections

- 2026-11-09: measured the water mass-transfer intercept against the replicate
  scatter, "so the last row is about 20 σ". The replicate pairs give the scatter of
  one run, not the uncertainty of an intercept extrapolated from a fit. The fit's
  own standard error is 0.0870 <!-- claim: t1_mass_intercept_se --> (15 runs, 13
  degrees of freedom), and the deficit is -3.08 <!-- claim: mass_intercept_t_water -->
  standard errors, beyond the two-sided 5 % point of Student's t. The water fit's
  residual scatter alone is 1.44 <!-- claim: t1_mass_fit_resid_rms_pct --> % rms.
- 2026-11-09: said of the authors' diffusivity that "the intercept says it is 1.117
  too large". The factor 1.117 <!-- claim: dv_used_over_dv_backed_out --> is the
  still-air row's `2.0/1.79`; the forced-convection intercept needs
  1.155 <!-- claim: dv_factor_from_intercept -->.

## Reading the sources

**No text layer at all** (`pdftotext` returns one byte per page on both parts), so
every digit came off a 400 dpi render (the files' native resolution; `pdfimages
-list` reports JPEG RGB tiles at 400 × 400 ppi, so rendering larger only
interpolates). Each numeric column was cropped on its own and magnified 2–8×, and
ambiguous glyphs were compared against unambiguous instances of both candidate
digits in the same column.

**Part I's page range, recorded both ways.** Part I's scan frames carry two folio
lines, one above the black edge of the leaf and one inside it at the foot, differing
by exactly two. Four of the six frames show both legibly (139/141, 140/142, 141/143,
144/146); on the fourth the upper line is cut off at the leaf edge and on the fifth
there is none, so those two pairs are inferred from the sequence rather than read.
Part II settles which belongs to the page: its first page has a blank top margin and
prints `Vol. 48, No. 4 | Chemical Engineering Progress | Page 173` at the **foot** of
the text block, so in this issue the folio is a footer and the upper line on a Part I
frame is the preceding leaf's footer caught in the same frame. Read that way Part I
is **141–146**, which is the universal citation and the one Ranz gives himself in
his 1993 ISI *Citation Classic*. All folios cited on the page are the ones printed
at the foot of the page the value appears on. Part II is 173–180 either way.

**Three unreadable cells and one digit settled by arithmetic**, both handled in the
open: an ink blot covers the `N_Nu'` entries of Table 1 runs 12–14, which are left
empty rather than guessed. Table 4 run 11's final glyph is degraded and the glyph's
counter position alone would suggest 10.6. It is settled instead by the table's own
arithmetic: both of Table 4's transfer numbers come from the same measured rate, so
their ratio depends only on the printed temperatures and pressure, and runs 10 and
12 (within 0.2 K, same pressure) bracket run 11's `N_Nu'`. With the first two glyphs
reading "10." and the last as wide as this column's round digits, the value is
**10.0**.

**Two files on disk are not this paper.**
`MISC-Ranz1993-citation-classic-commentary-CurrContents22.pdf` is a one-page 1993
reminiscence carrying none of the correlations, quoted once for Ranz's caveat on
radial convection, and `Charlesworth1960-evaporation-drops-dissolved-solids-AIChEJ6-9.pdf`
is a companion study on a different problem.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: all five datasets are committed CSVs in `data/`.
