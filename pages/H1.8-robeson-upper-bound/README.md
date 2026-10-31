# H1.8 — The Robeson upper bound, revisited

From Robeson, L. M., *The upper bound revisited*, *J. Membr. Sci.* **320**(1–2)
390–400 (2008),
[doi:10.1016/j.memsci.2008.04.030](https://doi.org/10.1016/j.memsci.2008.04.030).
The prior bounds are read from this paper's restatement of Robeson (1991) and
Robeson et al. (1994), which were not consulted.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the edited source, paired with `index.ipynb` (jupytext percent format) |
| `index.ipynb` | the page; runs in about 10 s |
| `model.py` | `P_bound`, `alpha_bound`, `excursion` and `envelope_fit`, the reusable functions |
| `test_model.py` | import, the round trip and sign convention, one stored claim via `load_claims` |
| `data/robeson-2008-table12-upper-bounds.csv` | Table 12, prior and present (k, n), `NA` structure kept |
| `data/robeson-2008-near-bound-points.csv` | Tables 1–11, the near-bound points |
| `data/robeson-2008-table13-knudsen-transition.csv` | Tables 13a/13b |
| `data/robeson-2008-printed-claims.csv` | the prose claims tested, verbatim |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Every numeric cell was read twice, from the born-digital text layer and from
digit-scale crops of 300 dpi renders: 382 <!-- claim: cells_total --> cells,
no disagreement. No figure was digitised.

## What the page is

A correlation audit of $P_i = k\,\alpha_{ij}^{\,n}$. There is no PDE, no grid
and no pymrm call in it; `structures` and `pymrm_api` are empty, deliberately.
Robeson's $(k,n)$ are not a fit: the line is drawn *"by eye"*, so nothing here
is a goodness-of-fit statement, and nothing is validation against measurement.

## The results

**Nine, not thirteen.** Only 9 <!-- claim: pairs_comparable --> of Table 12's
13 <!-- claim: pairs_in_table12 --> pairs have both a prior and a present
bound, and the paper's Conclusions classify exactly those nine.

**"Primarily the front factor" needs a reference.** $\Delta\log k$ is the
shift at $\alpha = 1$, where the membrane does not separate; that selectivity
lies inside the tabulated range of 5 <!-- claim: pairs_alpha1_inside_range -->
of the nine pairs and below every point of the others. At each pair's own
geometric-mean selectivity the front-factor share has median
0.98 <!-- claim: front_factor_share_median -->, range
0.59 <!-- claim: front_factor_share_min --> to
1.39 <!-- claim: front_factor_share_max -->; on the four references Robeson's
own Table 13 supplies it runs 0.5647 <!-- claim: v2_kshare_min --> to
1.0334 <!-- claim: v2_kshare_max -->. The two lines are near-parallel for
6 <!-- claim: pairs_near_parallel --> pairs and not for H2/N2
(0.34 <!-- claim: h2n2_spread_over_shift -->), H2/CH4
(0.93 <!-- claim: h2ch4_spread_over_shift -->) and H2/CO2
(2.05 <!-- claim: h2co2_spread_over_shift -->).

**Rank by selectivity gain, not by $k$.** The selectivity gain at fixed
permeability puts the paper's three "significant" pairs on top
(3 <!-- claim: dlogA_top3_hits --> of 3), with a third-to-fourth gap of
0.011 <!-- claim: dlogA_gap_min_over_eval_points --> to
0.111 <!-- claim: dlogA_gap_third_to_fourth --> decades depending on where it
is evaluated; $\Delta\log k$ gets 2 <!-- claim: dlogk_top3_hits --> of 3.

**The H2/CO2 bounds cross** at $\alpha^\star =$
37.47 <!-- claim: h2co2_crossing_alpha -->, inside Table 10's range, and at its
largest tabulated $\alpha$ (100.9 <!-- claim: h2co2_alpha_max -->) the 2008
bound sits 30.4 <!-- claim: h2co2_present_below_prior_at_alpha_max_pct --> %
below the 1994 one. "Inside the range" rests on
1 <!-- claim: h2co2_points_above_crossing --> tabulated point, which break row
(T1) removes.

**What "by eye" costs.** 19 <!-- claim: points_above_present_bound --> of the
117 <!-- claim: points_total --> points tabulated as close to the present
bound lie above it, by up to a factor 3.49 <!-- claim: excursion_max -->.

## Validation

- **V1**, Tables 13a/13b from Table 12: median deviation
  0.017 % <!-- claim: t13_identity_median_abs_rel -->,
  25 <!-- claim: t13_rows_inside_loose_band --> of 26 <!-- claim: t13_rows --> rows inside a generous
  printing band; 13a O2/N2 misses by -0.612 % <!-- claim: t13a_o2n2_rel -->.
- **V2** is a consistency identity, not a second route: its residual is the
  difference of two V1 residuals.
- **V3**, the second route: a covering line re-estimated from the 117 points.
  At each pair's geometric-mean selectivity it agrees on the size of the shift
  to 0.086 <!-- claim: envelope_offset_max_abs_wide_span --> decades on the
  well-levered pairs (up to 0.543
  <!-- claim: envelope_offset_max_abs_wide_span_max_over_alpha0 --> across
  references) and does not reproduce the front-factor share (median
  0.43 <!-- claim: envelope_implied_share_median -->).
- **V4**, the break table: single defects (sign and log conventions, glyph,
  digit and row errors in Tables 12, 1–11 and 13, two algebra slips), each
  claim recomputed under each; every claim fails at least one row.

## Corrections

- 2026-10-30: stated that $\alpha = 1$, where $\Delta\log k$ is the shift, lies
  "outside every dataset in the paper" and that "no pair's data reach" it, and
  that the envelope's $\Delta\log k$ is extrapolated "one to four decades
  outside the data". $\alpha = 1$ lies inside the tabulated range of
  5 <!-- claim: pairs_alpha1_inside_range --> of the nine pairs (He/CH4,
  He/CO2, H2/CO2, H2/CH4, He/N2) and below every point of the other four.
- 2026-10-30: stated that the envelope offset is "systematically positive".
  Its sign depends on the reference: at each pair's geometric-mean
  selectivity it is positive for 7 <!-- claim: envelope_offsets_positive -->
  pairs and negative for 4 <!-- claim: envelope_offsets_negative -->. What
  holds at every reference is the slope: the covering line is shallower than
  the printed one for 10 <!-- claim: envelope_dn_positive --> pairs.
- 2026-10-30: printed the V2 identity with the V1 residuals added; with
  $\rho = P_{k\alpha^n}/P_{\text{printed}} - 1$ they are subtracted.
- 2026-10-30: stated that the 117 points come from "roughly 65 primary
  papers"; they cite 44 <!-- claim: references_cited --> distinct references.
- 2026-10-30: stated that for He/H2 "the slope term contributes at most 3.0 %"
  of the shift; 3.0 % <!-- claim: heh2_spread_over_shift --> is how much the
  shift varies across the pair's range, and the slope term itself reaches
  3.1 % <!-- claim: heh2_slope_term_max_frac --> at the largest tabulated selectivity.
- 2026-10-30: the literal printing band of Table 13a's N2/CH4 row was computed
  from the product cell read as a number, which drops the printed trailing zero
  of `0.270`, giving that row's literal band as
  8.43 % <!-- claim: t13a_n2ch4_band_tight_float_read --> instead of
  0.92 % <!-- claim: t13a_n2ch4_band_tight -->; no count changes.
- 2026-10-30: stated that V1 would show a wrong digit in any transcribed
  (k, n) "at orders of magnitude above" its median deviation, and that V1 and
  V2 "genuinely pin the final printed digit of n(He/H2)". A last-digit error
  in n(He/H2) leaves all 4 <!-- claim: v1_heh2_rows_inside_tight_band -->
  He/H2 rows inside their literal printing bands; V1 detects misread glyphs,
  lost digits and dropped commas, and a last-digit error is caught only
  against the committed values.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/H1.8-robeson-upper-bound/index.py   # after editing index.py
python scripts/run_kernels.py pages/H1.8-robeson-upper-bound
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest -q pages/H1.8-robeson-upper-bound/test_model.py
```
