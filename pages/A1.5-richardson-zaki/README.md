# A1.5, Richardson–Zaki

`u/u_t = ε^n` is one of the most-quoted expressions in particle technology, and it
is not what the paper writes. The paper writes `log V_c = n log ε + log V_i` and
then spends four pages establishing that **`V_i` is not `V_0`** in fluidisation.
The famous form is the sedimentation special case with the wall term dropped.

Richardson and Zaki (1954) is a **data paper**: eight printed tables, 66 <!-- claim: n_comparisons --> measured
exponents, two wall laws, and only then the correlation, as five expressions,
one for each of 5 <!-- claim: rz_windows --> Reynolds windows. This page transcribes all eight tables off
600 dpi renders, certifies the transcription against identities the tables carry,
costs the correlation against a null model window by window, decides between the
two wall laws on the authors' own numbers, and puts the flux function into a pymrm
conservation law.

- **Structures:** `S1` (the pymrm solve; the correlation itself has no PDE in it)
- **Data tier:** 2, printed tables of the authors' own laboratory measurements
- **Nothing here is digitised.** No figure point was read.

## Source, and which pagination is cited

**Richardson, J. F. and Zaki, W. N.**, *Sedimentation and fluidisation: Part I*,
Trans. Instn Chem. Engrs **32**, 35–53 (1954).

**The document on disk is the verbatim reprint** in *Trans IChemE* Vol. 75,
December 1997 Jubilee Supplement,
[doi:10.1016/S0263-8762(97)80006-8](https://doi.org/10.1016/S0263-8762(97)80006-8),
where the paper occupies **reprint pages S82–S100**. Every page carries the
original running footer `TRANS. INSTN CHEM. ENGRS, Vol. 32, 1954`, but **the
original page numbers 35–53 are nowhere printed on the reprint.** Nineteen
consecutive pages is consistent with 35–53; the page does not assert that range
as read. **Every page reference on the page and in the sidecars is a reprint
page, marked as such.**

**Origins cited but not consulted.** Figure 18, from which eqs. (35)–(39) are
fitted, carries Steinour's sedimentation slopes (ref. 4) and Lewis, Gilliland and
Bauer's fluidisation slopes (ref. 9) as well as the authors' own. Neither paper is
on disk. Nothing on this page separates their contribution from the authors' own,
because the paper does not.

## The correlation, as printed

```
log V_c = n log ε + log V_i                              (28)   the actual form

n = 4.65 + 19.5 d/D                        Re < 0.2      (33)
n = (4.35 + 17.5 d/D) Re^-0.03      0.2 < Re < 1         (37)
n = (4.45 + 18   d/D) Re^-0.1         1 < Re < 200       (38)
n = 4.45 Re^-0.1                    200 < Re < 500       (39)
n = 2.39                                  Re > 500       (34)

V_i = V_0                              sedimentation     (40)
log V_0 = log V_i + d/D                fluidisation      (41)
n = 2.7 K^0.16,  K = (π/6) d_s³/d_p³   non-spherical     (42)
```

Validity, stated by the authors: `Re` from 2e-4 to 7e3, `d/D` from 0 to 4e-2.
`Re` is built on the **terminal** velocity, so `n` does not move as the bed expands.

## What the page checks

| route | check | result |
|---|---|---|
| 1, measurement | eqs. (33)–(39) vs. the 58 <!-- claim: n_rows_nonoil --> non-oil measured slopes | mean **1.20** <!-- claim: corr_nonoil_mean_abs_pct --> %, worst 6.75 <!-- claim: corr_nonoil_max_abs_pct --> %, **but this is a FIT RESIDUAL**, see below |
| 1, measurement | the **wall coefficient α**, fitted separately on each half of the data | **sedimentation (Table I, the rows that can fail) -0.235 <!-- claim: wall_alpha_sed_tableI --> [-1.78 <!-- claim: wall_alpha_sed_tableI_ci_lo -->, +0.16 <!-- claim: wall_alpha_sed_tableI_ci_hi -->]; fluidisation +1.095 <!-- claim: wall_alpha_fluidisation --> [+0.83 <!-- claim: wall_alpha_flu_ci_lo -->, +1.50 <!-- claim: wall_alpha_flu_ci_hi -->]**, each interval excludes the other law: Table I's upper end is more than 0.3 <!-- claim: wall_tableI_gap_to_one --> below 1, fluidisation's lower end more than 0.3 <!-- claim: wall_flu_gap_to_zero --> above 0 |
| 1, measurement | the paired wall test, which removes `V_0` | slope **0.79** <!-- claim: wall_paired_slope --> [0.47 <!-- claim: wall_paired_ci_lo -->, 1.34 <!-- claim: wall_paired_ci_hi -->] on the **6** <!-- claim: wall_paired_resolvable --> of 10 <!-- claim: wall_paired_total --> pairs that can resolve it |
| 1, measurement | eq. (42) vs. Table VIII, with a null baseline | rms 1.59 <!-- claim: eq42_rms_pct --> % vs. a constant's 7.94 <!-- claim: eq42_null_rms_pct --> %, **but 1.60 <!-- claim: eq42_rms_pct_no_plates --> % vs. 2.66 <!-- claim: eq42_null_rms_pct_no_plates --> % once the plates are dropped** |
| 2, identity | eq. (39) at Re = 500 vs. the separately-stated eq. (34) | **-0.015** <!-- claim: crossover_Re500_pct --> %, but the two share Table VI's last row (`n0` = 2.39 <!-- claim: tableVI_top_n0 --> at `Re` = 489 <!-- claim: tableVI_top_Re -->), so this is a **transcription check**, not two independent fits meeting |
| 2, identity | Table VI × Table VII, `n0 = slope × \|(d/D)_{n→0}\|` | worst 0.635 <!-- claim: tableVI_VII_max_pct --> %, mean 0.196 <!-- claim: tableVI_VII_mean_pct --> % |
| 2, identity | branch continuity | **-13.9 <!-- claim: crossover_Re200_dD004_pct --> % jump at Re = 200 and d/D = 0.04**; -1.82 <!-- claim: crossover_Re02_pct --> % at Re = 0.2; +2.30 <!-- claim: crossover_Re1_pct --> % at Re = 1 |
| 2, identity | the 17.5 and 18 against the division displayed above them (4.35/0.25 and 4.45/0.25) | both **rounded up**, by +0.57 <!-- claim: coeff_175_vs_174_pct --> % and +1.12 <!-- claim: coeff_18_vs_178_pct --> % |
| 2, identity | `K` recomputed from the printed `d_s`, `d_p` and definition | **fails on 4** <!-- claim: K_identity_rows_bad --> rows, worst 45.6 <!-- claim: K_identity_max_pct --> % |
| 2, identity | transcription identities | locate the printed errors listed on the page; Stokes reproduces Table I to 0.67 <!-- claim: stokes_tableI_max_pct --> % and correctly **fails** on Table II, by 19 <!-- claim: stokes_tableII_settling_min_pct --> to 66 <!-- claim: stokes_tableII_settling_max_pct --> % |
| 4, digitised |, | **not used. Nothing on this page is digitised.** |

## The things a reader should take away

**The mean deviation of 1.20 <!-- claim: corr_nonoil_mean_abs_pct --> % is a fit residual and the page refuses to call it anything else.**
Eqs. (33)–(39) were fitted to these very slopes. The informative number is the
**null baseline** beside it, per Reynolds window:

| branch | correlation rms | constant-n rms |
|---|---|---|
| eq. (33) | 0.95 <!-- claim: eq33_rz_rms_pct --> % | 3.07 <!-- claim: eq33_null_rms_pct --> % |
| eq. (37), 5 <!-- claim: eq37_rows --> rows | 1.559 <!-- claim: eq37_rz_rms_pct --> % | 1.538 <!-- claim: eq37_null_rms_pct --> %, **decides nothing**; the gap reverses under two of the three sensitivities the page prints and splits under the third |
| eq. (38) | **2.11** <!-- claim: eq38_rz_rms_pct --> % | **15.80** <!-- claim: eq38_null_rms_pct --> %, this is where it earns its keep |
| eq. (34) |, | **it *is* the constant**; the mean of the measured slopes is 2.3900 <!-- claim: eq34_mean_slope --> |

Pooled over all non-oil rows, one constant n = 4.00 <!-- claim: null_pooled_constant_n --> gives 31.2 <!-- claim: null_pooled_rms_pct --> % rms against the
correlation's 1.78 <!-- claim: corr_nonoil_rms_pct --> %; quoting that alone would be misleading.

**The wall term is the one genuinely two-sided test on the page, once the rows
that cannot fail are taken out of it.** Eqs. (40) and (41) are the same statement
with α = 0 and α = 1, and both sides are printed columns. But **Table II's
intercept column *is* its `log10 V_0` column**: on its usable 14 <!-- claim: wall_alpha_sed_tableII_zero_rows --> rows, across two and
three tube diameters per group, the residual the fit sees is zero to the printed
digits, against Table I's much larger one. Those rows cannot show a wall
effect of any size while still narrowing the interval, fitted alone they give
**α = +0.009** <!-- claim: wall_alpha_sed_tableII --> [+0.000 <!-- claim: wall_alpha_sed_tableII_ci_lo -->, +0.012 <!-- claim: wall_alpha_sed_tableII_ci_hi -->]. So the headline sedimentation fit is **Table I
alone, α = -0.235** <!-- claim: wall_alpha_sed_tableI -->, which still excludes +1; the pooled fit
(**-0.064** <!-- claim: wall_alpha_sed_pooled --> [-0.41 <!-- claim: wall_alpha_sed_pooled_ci_lo -->, +0.08 <!-- claim: wall_alpha_sed_pooled_ci_hi -->], 3.9 <!-- claim: wall_pooled_ci_narrowing --> times narrower) is printed beside it and the exclusion
is explained. The verdict survives the correction; the precision does not.

The limits are stated: the fluidisation α is a refit of the data Fig. 21 was
drawn from, so recovering 1 is nearly guaranteed; the **sedimentation** α is the
informative half, because eq. (40) was justified on a different figure on which a
wall term would not have shown as a slope, and a fluidisation-sized wall effect
would have appeared here. It does not.

**And the coefficient 1 does not survive its own range.** Restricted to the `d/D`
range the two data sets share, the fluidisation α is **+2.267** <!-- claim: wall_alpha_flu_common_range --> [+1.20 <!-- claim: wall_alpha_flu_common_ci_lo -->, +3.4 <!-- claim: wall_alpha_flu_common_ci_hi -->] ,
an interval that *excludes* 1 (its lower end is more than 0.1 <!-- claim: wall_flu_common_gap_to_one --> above it), against sedimentation's **-0.101** <!-- claim: wall_alpha_sed_common_range_tableI --> [-1.20 <!-- claim: wall_alpha_sed_common_ci_lo -->, +0.21 <!-- claim: wall_alpha_sed_common_ci_hi -->]
on the Table I rows there; and α falls from **+2.162** <!-- claim: wall_alpha_flu_dD_below_005 --> below `d/D = 0.05` to
**+1.033** <!-- claim: wall_alpha_flu_dD_above_005 --> above it, which a single linear term in `d/D` cannot do. The two-law
verdict (a wall term in fluidisation, none in sedimentation) stands; the *linear
form* of eq. (41) does not hold across the range it was fitted on. Run 12, the
sign-corrupted cell, is a Table II row and so is not in the headline fit at all:
on the pooled fit it reads -0.436 <!-- claim: wall_alpha_sed_with_run12 --> as printed and **-0.067** <!-- claim: wall_alpha_sed_run12_sign_restored --> with the minus sign
restored, a printing error, not a fragility.

**The solids balance in the pymrm section cannot fail for almost any defect, and
the break table says so.** It stays at round-off, below 1e-12 <!-- claim: pymrm_solids_residual -->, for a wrong `n`, a wrong `V_0`,
an unstable CFL, a coarse grid, a flipped flux sign and a reversed upwind
direction; the only defect it catches is `nu = 1`. The check that *can* fail is
the L1 error against the exact shock, which converges with successive orders
within 0.8 <!-- claim: pymrm_order_min --> and 1.2 <!-- claim: pymrm_order_max --> and moves ×111 <!-- claim: brk_nu1_ratio --> for `nu = 1` and ×187 <!-- claim: brk_upwind_ratio --> for a reversed upwind
direction, but only ×1.5 <!-- claim: brk_n3_own_speed_ratio --> when `n` is changed and judged against its own shock speed.
**The solve tests the integration, not the physics.**

## What is not here

- **No experimental validation.** The data are the authors' measurements, but the
  correlation was fitted to them and there is no hold-out set. The only rows the
  authors exclude, the two oils at high `d/D`, are excluded for a reason that is
  **perfectly confounded** with the correlation's own stated `d/D` limit: inside
  the eq. (38) window every row above `d/D = 4e-2` is an oil row and there are 0 <!-- claim: nonoil_highdD_in_eq38_window -->
  non-oil rows there at all.
- **No claim about ε outside roughly 0.42 to 0.96.** The extreme printed abscissa
  ticks read on any figure are `log ε = -0.38` and `-0.02` (Fig. 12), printed
  axis labels, not digitised coordinates.
- **No result resting on Table V's `d_s` column**, which is the one column of the
  paper that is illegible at 600 dpi in this scan.
- **No mixed-size or mixed-shape result.** The paper's own is a single figure with
  no tabulated numbers.
- **No claim that pymrm improves the correlation.** It does not. The solve
  integrates the flux the correlation becomes downstream, and the break table
  measures exactly how little that shows.
- **No verdict on eq. (37).** Its rows cannot separate the correlation from a
  constant in either direction; the page prints the three sensitivities that show
  why, and does *not* say eq. (37) is worse than a constant.
- **No claim that `log V_0 − log V_i` equals `d/D`.** The page claims a wall term
  in fluidisation and none in sedimentation. On the shared `d/D` range the
  fluidisation interval excludes eq. (41)'s coefficient of 1.
- **No page range asserted from the document.** The reprint numbers the paper
  S82–S100 and `35–53` appears nowhere on it; the two `cite_data` lines in the
  notebook print the catalogue's range with a flag under each saying it was
  inherited, not read.

## Corrections

- 2026-11-02: said the transcription identities "locate five printed errors" and
  that the printed and recomputed `Re` "agree to better than 0.5 % everywhere
  except two rows named in section 1". The two differ by more than 0.5 % on 4 <!-- claim: Re_rows_above_half_pct --> fluidisation rows: the two Table III ballotini rows (water properties assumed), a Table IV
  ballotini row, and Table IV lead shot in oil, which is 3.80 <!-- claim: transcription_Re_max_pct --> % off and was
  not named. It is now listed with the other anomalies, which makes 6 <!-- claim: printed_errors_listed --> printed errors in the list: the `d/D` of Table I runs H, N and S, the lost
  minus sign of Table II run 12, the Table IV ballotini `V_0` cell, and that lead-shot row.
- 2026-11-02: said the five expressions cover "four Reynolds windows". There are
  5 <!-- claim: rz_windows --> windows: below 0.2, 0.2 to 1, 1 to 200, 200 to 500 and above 500.
- 2026-11-02: said the gap between eq. (37) and a constant on its five rows
  "reverses under all three sensitivities the page prints" ("under every
  sensitivity"). That overstated it: two of the three reverse it, and leave-one-out
  splits, the correlation being ahead in 2 <!-- claim: eq37_loo_rz_wins --> of the 5 <!-- claim: eq37_rows --> deletions and the constant in the rest.
- 2026-11-02: said the individual-row bootstrap "would have given intervals
  roughly half as wide" as the cluster bootstrap for the residual correlations.
  It gives 0.75 <!-- claim: resid_ci_width_ratio --> times the width on average.
- 2026-11-02: said Table II's usable rows carry "an identically zero residual" and
  that the intercept column agrees on "14 of them exactly". The 14 <!-- claim: wall_alpha_sed_tableII_zero_rows --> rows agree to the
  printed decimals of the intercept column, not exactly.
- 2026-11-02: said the solids balance is "exactly zero" for a wrong `n`, a wrong
  `V_0`, an unstable time step and a 25-cell grid. It is at round-off, not
  exactly zero for all of them: below 1e-12 <!-- claim: pymrm_solids_residual -->.
- 2026-11-02: gave the wall-coefficient bootstrap intervals from 2000
  resamples to two decimals, with upper ends +0.17 (Table I), +1.51
  (fluidisation) and +3.47 (common range). At 2000 resamples the second
  decimal moved with the random seed. From 200000 resamples the upper ends are
  +0.16 <!-- claim: wall_alpha_sed_tableI_ci_hi -->, +1.50
  <!-- claim: wall_alpha_flu_ci_hi --> and +3.4
  <!-- claim: wall_alpha_flu_common_ci_hi -->, and every exclusion stated with
  them still holds.

## Related

`A1.1` (Ergun, the fixed-bed end of the same ε(u) curve), `A1.6` (Wen and Yu ,
the minimum-fluidisation matching condition; **not built here**), `A1.7` (Geldart
, decides which fluidised-bed model applies at all), `A1.8` (fluidisation regime
map), `E1.1` (the two-phase theory, the gas-solid analogue; parked).
