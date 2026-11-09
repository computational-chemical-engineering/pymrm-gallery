# A1.7: Geldart's powder groups

The four letters everyone quotes come from two inequalities. This page recomputes
both from the expressions Geldart derives them from, tests the A/B one against
the size fractions he measured himself, measures what the minimum-fluidization
correlation underneath it is actually worth against those same measurements, and
shows that the B/D one is the Davidson cloud-existence condition in other symbols.

- **Structures:** `S3` (the single pymrm solve, E1.2's model copied unchanged; the
  classification itself has no PDE in it)
- **Runtime:** ~4 s
- **Data tier:** 2, a table printed in the paper, and the values in it are the
  author's own laboratory measurements.

## Source

**Geldart, D.**, *Types of gas fluidization*, Powder Technology **7**(5) 285–292
(1973), [doi:10.1016/0032-5910(73)80037-3](https://doi.org/10.1016/0032-5910(73)80037-3)
(the paper itself was read, not a reprint of it, on 600 dpi renders of journal
pages 285–292). The PDF text layer of this Elsevier scan renders the paper's own
title as "Types of Gas Fhidization" and was not used for anything.

**Origins cited but not consulted.** Equation (3), the minimum fluidization
velocity, carries Geldart's reference 11: *L. Davies and J. F. Richardson, Trans.
Inst. Chem. Engrs. 44 (1966) T293*. Not on disk; used exactly as Geldart prints
it. Same for equation (9) (Verloop and Heertjes, ref. 10) and equation (10)
(Oltrogge, ref. 36). The Davies–Taylor coefficient used in section 5 comes
from `E1.2`, not from this paper.

## The two inequalities

```
group A :  (rho_s - rho_f) d'   <= 225      eq. (6),  d' in um, rho in g/cm3
group D :  (rho_s - rho_f) d'^2 >= 1e6      eq. (8)
```

Recomputed from equations (5) and (7) with Geldart's own `g = 981`,
`mu = 1.8e-4`, `K_MB = 100`, `eps_0 = 0.4` and `d_B = 25 cm`, they are **229.36 <!-- claim: eq6_constant_recomputed -->**
and **1.016e6 <!-- claim: eq8_constant_recomputed -->**: his printed values are 1.94 <!-- claim: eq6_constant_dev_pct --> % and
1.59 <!-- claim: eq8_constant_dev_pct --> % low. Neither rounding changes any verdict in his own table: 0 <!-- claim: AB_rows_in_225_to_229_gap --> rows
fall between 225 and the recomputed constant, and 0 <!-- claim: BD_rows_in_1e6_gap --> between `1e6` and its recomputed
value, the largest `(rho_s-rho_f)d'^2` in Table 1 being 1.192e5 <!-- claim: BD_max_prod8_in_table1 -->.

## What the page checks, in the brief's order

| route | check | result |
|---|---|---|
| 3, measurement | A/B criterion vs. the measured `U_MB/U_0` of Table 1 | **2 <!-- claim: AB_discriminating_rows_agreeing --> of the 3 <!-- claim: AB_discriminating_rows --> rows that can discriminate** (whole table 20 <!-- claim: AB_rows_agreeing -->/21 <!-- claim: AB_rows_compared -->, against a null baseline of 19 <!-- claim: AB_rows_null_baseline -->) |
| 3, measurement | eq. (3) vs. the measured `U_0` column | implied constant spans 0.46 <!-- claim: eq3_K0eff_ratio_min -->–2.94 <!-- claim: eq3_K0eff_ratio_max --> × the printed `8e-4`, a factor 6.4 <!-- claim: eq3_K0eff_spread_factor -->, median 0.98 <!-- claim: eq3_K0eff_ratio_median -->; geometric mean 0.998 <!-- claim: eq3_K0eff_geomean --> with a standard error of 9.6 <!-- claim: eq3_K0eff_geomean_se_pct --> % |
| 3, text | the worked example on journal page 289 | worst deviation -1.94 <!-- claim: worked_sand_U0_dev_pct --> % (`U_0` of 100 µm sand), a coincidence of digits with the 1.94 <!-- claim: eq6_constant_dev_pct --> % above and an unrelated quantity |
| 2, identity | eqs. (6) and (8) recomputed from (5) and (7) | 1.94 <!-- claim: eq6_constant_dev_pct --> % and 1.59 <!-- claim: eq8_constant_dev_pct --> %, and no Table 1 row in either gap |
| 2, identity | solids balance linking `eps_MB` and `H_MB/H_0` | all rows give an admissible `eps_0`, 0.441 <!-- claim: eps0_inferred_min -->–0.570 <!-- claim: eps0_inferred_max --> |
| 2, identity | eqs. (9) and (10) placed relative to (6) | eq. (10) crosses XY inside the plot at 104 <!-- claim: eq10_crosses_XY_um --> µm; eq. (9) lies 2.4 <!-- claim: eq9_over_XY_size_ratio_bottom -->–8.0 <!-- claim: eq9_over_XY_size_ratio_top --> × to the right of XY across it |
| 2, cross-page | eq. (7) is E1.2's no-cloud condition | `0.711` vs. `1/sqrt(2)` moves the boundary by 0.55 <!-- claim: eq8_davies_taylor_vs_geldart_pct --> % |
| 2, cross-page | the copy of E1.2's model reproduces E1.2's stored cloud claims | relative gap at round-off, bounded by E1.2's own `rtol` of each claim |
| 4, digitised | none | **not used. Nothing on this page is digitised.** |

**The A/B score is reported the honest way round.** 19 <!-- claim: AB_rows_null_baseline --> of the 21 <!-- claim: AB_rows_compared --> comparable rows
are measured group A, so a null model that predicts group A for every powder
already scores 19 <!-- claim: AB_rows_null_baseline --> (Diakon 6 <!-- claim: AB_rows_null_baseline_diakon -->/8 <!-- claim: AB_rows_compared_diakon -->) against the criterion's 20 <!-- claim: AB_rows_agreeing --> (7 <!-- claim: AB_rows_agreeing_diakon -->/8 <!-- claim: AB_rows_compared_diakon -->): one row
of information. The only rows on which the two can differ are the three the
criterion calls group B, and there it is right 2 <!-- claim: AB_discriminating_rows_agreeing --> of 3 <!-- claim: AB_discriminating_rows -->. The margin also rests on a
**tie-break**: 2 <!-- claim: AB_tie_rows --> Diakon rows print `U_MB = U_0` exactly, eq. (4) is *displayed*
as `>= 1` but *applied in section 4.3* as `> 1`, and the strict form is used
because it is the only one that reproduces Geldart's own verdict on those rows.
Break table 2 costs the alternative at 18 <!-- claim: AB_rows_agreeing_nonstrict -->/21 <!-- claim: AB_rows_compared --> (Diakon 5 <!-- claim: AB_rows_agreeing_diakon_nonstrict -->/8 <!-- claim: AB_rows_compared_diakon -->, discriminating rows
0 <!-- claim: AB_discriminating_rows_agreeing_nonstrict -->/3 <!-- claim: AB_discriminating_rows -->, null baseline 21 <!-- claim: AB_rows_null_baseline_nonstrict -->/21 <!-- claim: AB_rows_compared -->). The last group A row sits 0.46 <!-- claim: AB_last_A_row_inside_pct --> % inside the
boundary and is measured group A (`U_MB/U_0` = 1.054 <!-- claim: AB_last_A_row_ratio -->). Half a unit in the last printed digit of its `rho_s` and `d'` moves its product by 0.69 <!-- claim: AB_last_A_row_rounding_pct --> %, more than that margin; at the upper rounding edge it is group B and the agreement is 19 <!-- claim: AB_rows_agreeing_upper_rounding --> of 21 <!-- claim: AB_rows_compared -->.

On the single row where eq. (6) and the strict ratio part company (the 220 µm
Diakon fraction, `(rho_s-rho_f)d'` = 259.3 <!-- claim: AB_disagree_row_prod6 -->, `U_MB/U_0` = 1.029 <!-- claim: AB_disagree_row_ratio -->), **Geldart's own prose sides with
eq. (6)**: section 4.3 puts that fraction in group B, as eq. (6) does, against the
mechanical strict ratio.

**Break tables.** Every claim has a break row built from a call that injects the
defect (a physical parameter, a unit, a misread printed constant, a wrong column
or a mis-transcribed cell), and each row records which claims fail. Four
blind spots are named rather than asserted away: `rho_f` leaves every count and
classification unchanged, the spread of the constant implied by eq. (3) catches a
wrong exponent of 3 but not of 1 (the median catches both), an error of `0.10` in one
`eps_MB` is invisible to the voidage check, and the whole-table A/B score is nearly
powerless because the classes are lopsided.

`eps_MB` and `H_MB/H_0` are described as **reported**, not measured: section 4.2
states the method for `U_0` and `U_MB` only and says nothing about those two.

## What is not here

- **No C/A boundary.** Geldart prints no expression for it; on his Figure 3 it is
  a hand-drawn shaded band.
- **No test of the B/D boundary against data.** Table 1 stops at its largest fraction and its
  densest powder; at that powder's density equation (8) does not begin until 817 <!-- claim: BD_eq8_onset_um_at_rho_s_1p5 --> µm.
- **No claim that pymrm improves the classification.** It does not. The one solve
  on the page answers a question the classification raises and does not answer
  (how big the cloud is on either side of the B/D line) and uses `E1.2`'s model,
  copied unchanged, to do it: R_c/R = 1.025 <!-- claim: cloud_RcR_sand_100um --> in a 100 µm sand and 3.23 <!-- claim: cloud_RcR_sand_600um --> in a 600 µm sand, and no dividing streamline in a 700 µm sand.

## Corrections

- 2026-11-09: break table 2's note said "The Diakon column is the one with power. 14 of the 22 fractions sit so" far inside group A that no plausible mis-reading of the boundary moves them. Of the rows with a measured `U_0`, 14 <!-- claim: AB_rows_at_most_150_comparable --> have a product of at most `150`, the lowest mis-read boundary in that table, and 15 <!-- claim: AB_rows_at_most_150 --> of all 22 <!-- claim: table1_rows --> rows do.
- 2026-11-09: said that "the eq. (3) SPREAD barely moves under a wrong exponent (6.35 -> 5.69) while its median moves by three orders of magnitude". The spread factor of the implied constant is 6.35 <!-- claim: eq3_K0eff_spread_factor --> as printed, 5.69 <!-- claim: eq3_exponent1_spread_factor --> for an exponent read as 1 and 47.3 <!-- claim: eq3_exponent3_spread_factor --> for 3, so a scatter catches the second misreading and not the first; the median moves by a factor 0.00998 <!-- claim: eq3_exponent1_median_factor --> and 81.8 <!-- claim: eq3_exponent3_median_factor -->, two orders of magnitude and not three.
- 2026-11-09: said that, with the median of the implied constants at 0.981 <!-- claim: eq3_K0eff_ratio_median --> times the printed `8e-4`, that value "is the right number for this collection of powders" to within the median's distance from 1. The geometric mean is 0.998 <!-- claim: eq3_K0eff_geomean --> times the printed value with a standard error of 9.6 <!-- claim: eq3_K0eff_geomean_se_pct --> %, so the data agree with `8e-4` but place the centre no better than that.
- 2026-11-09: said that "The eq. (6)/(8) constants move by 90 % to 10000 % under every plausible mis-reading of mu, g, K_MB, the 8e-4, eps_0 or d_B". A decade slip in `mu`, `K_MB`, the `8e-4` or `g` in SI moves at least one of the two constants by 90 <!-- claim: constants_slip_min_shift_pct --> % or more; `eps_0 = 0.5` and `d_B = 50 cm` move only the eq. (8) constant, by 25 <!-- claim: eq8_constant_shift_eps0_pct --> % and 41 <!-- claim: eq8_constant_shift_dB_pct --> %.
- 2026-11-09: described the range of the constant implied by eq. (3) for the rows "restricted to Diakon above 100 um, the only fractions with an exactly printed density". All 8 <!-- claim: rows_density_exact --> Diakon fractions have an exactly printed density; the 7 <!-- claim: eq3_diakon_rows_above_100um --> above 100 µm give 0.79 <!-- claim: eq3_K0eff_ratio_min_diakon_gt100 -->–1.16 <!-- claim: eq3_K0eff_ratio_max_diakon_gt100 -->, and all 8 give 0.79 <!-- claim: eq3_K0eff_ratio_min_diakon -->–1.66 <!-- claim: eq3_K0eff_ratio_max_diakon -->.
- 2026-11-09: said that the crossing of eq. (9) with line XY is "far BELOW the diagram rather than beyond its right edge". It is outside the labelled diagram in both coordinates: at 16310 <!-- claim: eq9_crosses_XY_um --> µm, 16.3 <!-- claim: eq9_cross_size_over_axis_max --> times the right end of the labelled abscissa, and at 0.0138 <!-- claim: eq9_crosses_XY_drho --> g/cm³, a factor 14.5 <!-- claim: eq9_cross_axis_min_over_drho --> below the lower end of the labelled ordinate.

## Related

`E1.2` (the Davidson bubble: owns the percolation solve and the cloud threshold,
and its `model.py` is copied here), `E2.1` (the bubbling-bed model, which presumes
a cloud), `A1.6` (Wen and Yu, the obvious replacement for equation (3); **not
built here**), `A1.1` (Ergun), `A1.8` (fluidisation regime map).
