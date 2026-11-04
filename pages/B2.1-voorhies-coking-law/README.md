# B2.1 : Voorhies' coking law

**Catalog ID:** `B2.1` · **Section:** B · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1` · **Data tier:** 2 (tables printed in the paper, the
author's own laboratory measurements)

Reproduces Alexis Voorhies, Jr., "Carbon Formation in Catalytic Cracking",
*Industrial & Engineering Chemistry* **37**(4) 318–322 (1945),
doi:10.1021/ie50424a010, the origin of $C_c = A\,\theta^n$, the empirical law
every time-on-stream deactivation model descends from. The paper's sharper,
non-obvious claim is the one this page tests: carbon on catalyst at a given
residence time is, "within limits, independent of the hydrocarbon feed rate".

Built from the 300 ppi CCITT-G4 scan; every numeral read from cropped
native-resolution renders at digit scale (the exponents $n$ and coefficients
$A$ *are* the case). The by-line prints "ALEXIS VOORHIES, JR.", the Jr. the
catalogue omits is carried here. **No figure was digitised**: Figs. 1, 3, 4, 5
are figure-only and out of scope, so the time exponent $n$ enters only as the
author's printed fit and is *never tested against data on this page* (Tables I
and II sit at a single $\theta$).

## What the page shows

- **"Independent of feed rate" quantified, with standard errors.** Across ten
  groups (eight Table-I temperature blocks, a 2-fold feed-rate span; cetane and
  Decalin in Table II, a 4-fold span) the feed-rate exponent of carbon on catalyst
  is $m_{cat}$ = -0.116 <!-- claim: m_cat_tableI_worst_block --> to
  +0.034 <!-- claim: m_cat_group_max -->, closer to 0 (the paper's claim) than to +1 (the throughput hypothesis) by a factor of
  9.6 <!-- claim: m_cat_closer_to_zero_ratio_min --> in the least favourable group, while the same rows put carbon on feed at $m_{feed}$ ≈ −1 (-0.95 <!-- claim: m_feed_group_max -->
  to -1.18 <!-- claim: m_feed_group_min -->). Pooled per table, with a separate
  intercept for every group: Table II gives -0.042 <!-- claim: m_cat_tableII_pooled -->
  ± 0.029 <!-- claim: m_cat_tableII_pooled_se -->; Table I gives -0.050 <!-- claim: m_cat_tableI_pooled -->
  ± 0.018 <!-- claim: m_cat_tableI_pooled_se -->, which is not exactly zero, and
  +0.022 <!-- claim: m_cat_tableI_withT --> ± 0.022 <!-- claim: m_cat_tableI_withT_se -->
  once one temperature coefficient inside the blocks is allowed for (the
  temperature varies within a block and correlates with the feed rate). The paper's claim is at
  fixed temperature, so the second fit is its test. Every fit
  is tens of standard errors from +1 (the smallest distance, 35 <!-- claim: m_cat_pooled_from_one_se_min -->, is Table II's). The standard
  error exceeds the rounding floor of the table in 8 <!-- claim: groups_se_exceeds_floor -->
  of the ten groups, so the floor is the resolution of the print, not the error bar of
  a slope. Verdict: a feed-rate exponent of a few hundredths at most, incompatible with
  +1, indistinguishable from 0 once temperature is accounted for, and not "proved
  exactly zero".
- **The two carbon columns are one measurement, so each table carries one
  test, not two.** The implied catalyst-to-oil ratio
  $C_c/(C_f\,U\,\theta/60)$ is constant within every group to a worst-group CV of
  2.9 <!-- claim: WU_within_group_cv_max_pct --> %, at the 2-significant-figure rounding floor
  (a half-digit allowance of up to 5.9 <!-- claim: WU_half_digit_allowance_max_pct --> %
  in the worst row), consistent with both columns being computed from the single printed
  carbon determination (combustion of the discharged catalyst). The identity doubles as
  a per-row transcription guard, recovers a ratio the paper never prints for these units
  (per-system means 0.86 <!-- claim: WU_implied_cetane --> to 1.09 <!-- claim: WU_implied_synthetic -->,
  not the value printed for the fluid-derivation unit), and surfaces a drift of that
  ratio across the temperature blocks of 7.9 <!-- claim: WU_block_drift_natural_pct --> %
  (natural) and 8.7 <!-- claim: WU_block_drift_synthetic_pct --> % (synthetic), 4.6 <!-- claim: WU_block_drift_natural_over_se -->
  and 9.1 <!-- claim: WU_block_drift_synthetic_over_se --> standard errors of the within-block
  scatter (block means of four rows each, about six degrees of freedom), reported, not resolved.
- **The printed algebra chain (2)+(3)+(4)→(5), by two routes with independent
  solvers** (symbolic isolation in exact rationals; the same isolation root-found
  with `pymrm.checks.find_roots` on the physical range of a conversion and fitted by
  least squares): the conversion is 95.7145 <!-- claim: eq5_coeff_symbolic --> divided by U to the power 0.34130 <!-- claim: eq5_expU_symbolic --> and by θ to the power 0.19113 <!-- claim: eq5_expTheta_symbolic -->,
  printed as "96 / 0.34 / 0.19", every printed digit the correct rounding. The routes agree
  to below 1e-9 <!-- claim: eq5_route_gap_max -->. None of the six <!-- claim: eq5_misreadings_total -->
  leading-digit, transposed-glyph or exponent-scale misreadings tried rounds to the printed
  triple; the last digits are within the printed rounding (6 <!-- claim: eq5_last_digit_alternatives_matching -->
  of the 10 <!-- claim: eq5_last_digit_alternatives_total --> last-digit ±1 alternatives still round to it)
  and rest on the digit-scale crops, which were read independently twice, at transcription and
  at verification, in agreement. Bonus exact object: the chain implies instantaneous
  conversion = 0.809 <!-- claim: inst_over_avg_conversion --> × the period average, everywhere.
- **"Doubles per 190–200 °F" measured per catalyst, with standard errors** (goodness of fit,
  the claim was derived from Table I): natural 205.6 <!-- claim: doubling_dT_natural_F -->
  ± 7.1 <!-- claim: doubling_dT_natural_se_F --> °F, which is 0.8 <!-- claim: doubling_natural_vs_band_high_se -->
  standard errors from the band's upper end, so Table I does not contradict the band for the
  natural catalyst; synthetic 175.8 <!-- claim: doubling_dT_synthetic_F -->
  ± 4.8 <!-- claim: doubling_dT_synthetic_se_F --> °F, -3.0 <!-- claim: doubling_synthetic_vs_band_low_se -->
  standard errors from the band's lower end (negative: below it) on a straight line through its rows, which the rows reject (quadratic term -3.6 <!-- claim: doubling_quadratic_t_synthetic --> standard errors); on its four block means it is -1.5 <!-- claim: doubling_synthetic_vs_band_low_blockmean_se --> standard errors (2 degrees of freedom). The two
  catalysts differ by 29.8 <!-- claim: doubling_dT_split_F --> °F, 3.5 <!-- claim: doubling_dT_split_over_se -->
  standard errors of the difference on the 16 rows and 1.8 <!-- claim: doubling_dT_split_over_se_blockmean -->
  on the four block means each. A two-point route on the extreme blocks gives
  208.3 <!-- claim: doubling_dT_natural_twopoint_F --> and 174.0 <!-- claim: doubling_dT_synthetic_twopoint_F -->
  °F from the same rows.
- **A printed factor-2 defect settled from the paper's own two lines**: the
  Nomenclature's "$K$ = constant $(= A^2)$" is inconsistent with eq. (8),
  integration gives $K = A^2/2$ (symbolically, and by ODE integration, which
  overshoots `A theta^0.5` by exactly $\sqrt2$ = 1.4142 <!-- claim: K_printed_numeric_overshoot --> with the
  printed $K$). Consequence-free in the paper ($K$ appears nowhere else); reported, never repaired.
- **$A$ is local, $n$ is portable, measured from the paper's own two campaigns**: eq. (2)
  against Table I at the same nominal conditions differ by a factor 1.82 <!-- claim: A_campaign_gap_ratio -->
  in $A$ (different unit, fresh against regenerated batches), while the four printed exponents stay
  in 0.38 <!-- claim: n_printed_min -->–0.53 <!-- claim: n_printed_max -->, exactly the range Froment
  & Bischoff quote, reconciled against `B2.2`'s transcription in the notebook.
- **The "faulty point" remark made quantitative** (labelled inference): the printed
  "36 % to 41 %" pair implies via eq. (3) that the flagged run's carbon yield sat a factor
  1.46 <!-- claim: faulty_point_excess_ratio --> above the correlation line.

## Fit vs test

Nothing on this page fits a model: the least-squares fits (feed-rate exponents,
temperature coefficients, the power-law form of the numeric route) are the statistics of
the comparisons, each stated where it is made in the notebook. The slope and
doubling tests are consistency checks of the paper's printed claims against the paper's
own printed tables (the doubling claim is a goodness-of-fit check, Fig. 2 was drawn from
Table I); the eq.-(5) and eq.-(8) work is reproduction of printed algebra; eqs. (1), (2), (3),
(6), (7) are the author's fits to figure-only data, transcribed and never refit. The same
labels appear in the notebook, `meta.yaml`, and the data sidecars.

## Corrections

- 2026-11-04: "The printed band is a fair pooled summary and sits outside its own table's support for either catalyst individually" was wrong for the natural catalyst. Its interval, 205.6 <!-- claim: doubling_dT_natural_F --> ± 7.1 <!-- claim: doubling_dT_natural_se_F --> °F, is 0.8 <!-- claim: doubling_natural_vs_band_high_se --> standard errors from the band's upper end; the natural catalyst is also inside on its four block means, 0.4 <!-- claim: doubling_natural_vs_band_blockmean_se --> standard errors from the upper end. The synthetic catalyst is -3.0 <!-- claim: doubling_synthetic_vs_band_low_se --> standard errors from the band's lower end on a straight line through its rows, which the rows reject (quadratic term -3.6 <!-- claim: doubling_quadratic_t_synthetic --> standard errors), and -1.5 <!-- claim: doubling_synthetic_vs_band_low_blockmean_se --> on its four block means (2 degrees of freedom), so the tables do not show it outside the band either.
- 2026-11-04: "the floors ARE the verdict's error bar" was wrong: the rounding floor bounds the print only, and the standard error of the slope exceeds it in 8 <!-- claim: groups_se_exceeds_floor --> of the ten groups.
- 2026-11-04: "The implied-$W/U$ temperature drift (~5–8 % down across blocks) is unexplained." understated the drift: it is 7.9 <!-- claim: WU_block_drift_natural_pct --> % for the natural and 8.7 <!-- claim: WU_block_drift_synthetic_pct --> % for the synthetic catalyst. It is still unexplained.

## Files

- `index.py`, `index.ipynb`: the page, as a jupytext pair; edit `index.py`. Runs in
  about 12 seconds.
- `data/voorhies-1945-tableI.csv`, `-tableII.csv`: the two printed tables,
  with per-cell provenance and the glyph notes (the two `1.1` v/v/hr readings)
  in the sidecars.
- `data/voorhies-1945-printed-constants.csv`: every printed constant and
  quantitative prose claim, with locations.
- `agreement.json`: every number the page states, written by its Claims cell and
  checked by CI; every claim is failed by at least one row of the break table, or is
  named structural there.
- `meta.yaml`: page metadata, validated against `models.yaml`.

Cross-page: loads `froment-bischoff-1961-printed-claims.csv` from `B2.2` only
to reconcile the quoted exponent range (checked equal to this page's min and max;
no borrowed number reaches any reported result).
