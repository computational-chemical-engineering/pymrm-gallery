# A1.8: Two ways to close gas-solid drag

One scalar decides almost everything a two-fluid fluidised-bed simulation does,
and the MFIX Theory Guide says in one paragraph why there is no single answer to
it. This page builds the two families that paragraph names, divides them, and
puts each against a measurement it never saw.

- **Structures:** `S3`, algebraic closures. **No pymrm operator appears** except
  `pymrm.checks.find_roots`, which locates every crossing, threshold and extremum
  on its physical range, and the page says so where a reader would look for a
  solver, as `A1.1` and `A1.6` do.
- **Runtime:** about 50 s
- **Data tier:** 2, two tables printed in papers, neither closure fitted to
  either.

## Source

**Syamlal, M., Rogers, W. and O'Brien, T. J.**, *MFIX Documentation Theory
Guide*, Technical Note DOE/METC-94/1004 (DE94000087), Morgantown Energy
Technology Center, U.S. Department of Energy (December 1993),
[doi:10.2172/10145548](https://doi.org/10.2172/10145548), freely downloadable
from OSTI. Section 2.2.1, journal pages 10-11, read on 400 dpi renders, which is
the embedded bilevel page images' own native resolution.

**Do not use this file's text layer for any digit.** It renders `0.06 Re` as
`0.O6Re` with a capital letter O, ε_g^`2.65` as `_g-2'6s`, "Ergun (1952)" as
`Ergun (f952)`, and the equation numbers (12) and (16) as `(121` and `(1''`.

**The origin of eq. (11) is an unpublished report.** The report's own reference
list gives Syamlal and O'Brien (1987), *"A Generalized Drag Correlation for
Multiparticle Systems," Unpublished report*, so this DOE report is not a
convenient source for the terminal-velocity-to-drag conversion, it is the citable
**published** one. The catalogue's citation for the case, "Syamlal & O'Brien
(1989)", is a *different* item in the same list: "Computer Simulation of Bubbles
in a Fluidized Bed", AIChE Symposium Series No. 270, **85**, 22-31, a
bubble-simulation paper. Both were read off a 400 dpi render of report page 44.

**Origins cited but not consulted.** Garside & Al-Dibouni (1977) for eqs.
(12)-(14) and Dalla Valle (1948) for eq. (16); both are printed in full with
attribution in the 1993 report, which is the reprint route of `AGENTS.md`.

**Also read directly:** Richardson, J. F. and Zaki, W. N., *Sedimentation and
fluidisation: Part I*, Trans. Instn Chem. Engrs **32** (1954) 35-53, from the
verbatim Golden Jubilee reprint, Trans IChemE **75** (Dec 1997) S82-S100,
[doi:10.1016/S0263-8762(97)80006-8](https://doi.org/10.1016/S0263-8762(97)80006-8),
on 300 dpi renders. Its Table VI is shipped here. That paper is `A1.5`'s source;
`A1.5` is untouched and only one of its tables is used.

## Scope: read this before anything else

The catalogue asks for **Gidaspow / Syamlal-O'Brien / Wen-Yu**. The document on
disk carries **Syamlal-O'Brien complete and neither of the other two**:

- the strings "Wen and Yu" and "Wen & Yu" do **not occur anywhere** in the
  report;
- the report prints **no Gidaspow drag closure and no blend rule**. Gidaspow
  (1986), "Hydrodynamics of Fluidization and Heat Transfer: Supercomputer
  Modeling", Appl. Mech. Rev. **39**, 1-23, is cited exactly once in §2.2.1, for
  the observation that a packed-bed correlation must be supplemented at low
  solids fraction.

  The name Gidaspow runs through the introduction, the ill-posedness
  discussion, the solids-stress section and the heat-transfer section, and the
  report **adopts a Ding & Gidaspow (1990) expression as its own eq. (88)**
  (granular-energy transfer, §2.5.3) and a Syamlal & Gidaspow (1985)
  conductivity model in the heat-transfer section. What is load-bearing, and
  true, is only that no Gidaspow *drag* law and no blend rule appear anywhere.

Nothing about either is written here from memory. What replaces them is the pair
the report's own paragraph on journal page 10 points at, both on disk: **Ergun
(1952)**, via `A1.1`, and **Richardson & Zaki (1954)**. Two documents would
complete the case, and they are named on the page and in `models.yaml`:
Gidaspow's 1994 monograph, and the Wen & Yu Chem. Eng. Progr. Symposium Series
paper (**not** the A.I.Ch.E. Journal communication on disk, which is `A1.6`'s
source and contains no drag law).

## Contents

| file | what it is |
| --- | --- |
| `index.py`, `index.ipynb` | the page (edit `index.py`, then `jupytext --sync`); the paired notebook is executed and rendered |
| `data/mfix-1993-syamlal-obrien-constants.csv` | the **eleven** printed constants of eqs. (11)-(16), each tagged with its equation number |
| `data/richardson-zaki-1954-table6.csv` | Table VI: the tabulated velocity-voidage indices over the table's Reynolds range, plus the paper's own two log₁₀ columns. Each n₀ is a fitted slope, then extrapolated along a second fitted line, so "tabulated", not "measured" |
| `agreement.json` | the claims the notebook registers, written by its last cell and checked by CI |
| `meta.yaml` | page metadata |

Both CSVs have a `.meta.yaml` provenance sidecar. Three further datasets are
**borrowed**: `A1.1`'s Ergun constants and its marker file, and `A1.7`'s Geldart
table. The page's *The data* section lists what those pages have already
established about the rows it uses.

## What the page checks

| route | check | result |
|---|---|---|
| 2, identity | `V_rm(eps_g = 1) = 1` at every Re | residual below 1e-9 <!-- claim: check1_Vrm_unity_identity -->: the radicand of eq. (12) is a perfect square. **Blind to all three voidage exponents, to the `0.8` and to the switch voidage** (1^x = 1): 0 <!-- claim: check1_voidage_constants_seen --> of those five are seen, and 4 <!-- claim: check1_readings_seen --> of ten injected mis-readings |
| 2, identity | the force balance returns `v_slip = V_rm v_t` | residual below 1e-9 <!-- claim: check2_force_balance_identity --> over four decades of size. Catches 5 <!-- claim: check2_structural_readings_seen --> of five structural mis-readings of eq. (11), with residuals from 0.28 <!-- claim: check2_smallest_structural_residual --> to 23.1 <!-- claim: check2_largest_structural_residual -->; **blind to every constant of eqs. (12) and (16)**: 0 <!-- claim: check2_constant_readings_seen --> of four seen |
| 3, measurement | vs Richardson & Zaki's Table VI, no fitting | implied index biased **+3.2 <!-- claim: check3_index_bias_pct_eps045 --> %** (standard error 1.2 <!-- claim: check3_index_bias_se_pct_eps045 --> %) at eps_g = `0.45`, **+5.4 <!-- claim: check3_index_bias_pct_eps050 --> %** (1.1 <!-- claim: check3_index_bias_se_pct_eps050 --> %) at `0.50`, **+29.2 <!-- claim: check3_index_bias_pct_eps095 --> %** (2.9 <!-- claim: check3_index_bias_se_pct_eps095 --> %) at `0.95`, growing with voidage because the closed form is not a power law. In the viscous limit the index is the exponent of A plus one against their tabulated asymptote: **+10.5 <!-- claim: check3_viscous_index_vs_RZ_pct --> %**. **Conditional on a reading, see below** |
| 3, measurement | vs Geldart's measured minimum-fluidisation velocities of the eight Diakon cuts | **all three overpredict at the derived voidage**: Ergun **+42.1 <!-- claim: check4_ergun_bias_pct --> %**, Syamlal-O'Brien **+101.2 <!-- claim: check4_syamlal_obrien_bias_pct --> %**, Richardson-Zaki **+111.1 <!-- claim: check4_richardson_zaki_bias_pct --> %**, each on 8 <!-- claim: check4_rows_overpredicted_ergun --> of eight rows. No closure is validated; Ergun has the smallest absolute bias at 4 <!-- claim: check4_ergun_least_abs_bias_voidages --> of five voidages tried, each by at least 3.98 <!-- claim: check4_ergun_lead_min_over_se --> standard errors of the paired difference; Syamlal-O'Brien and Richardson-Zaki are separated by more than two standard errors at only 1 <!-- claim: check4_so_rz_order_resolved_voidages --> of the five, so no order of those two is claimed |
| 2, identity | Table VI against its own two log₁₀ columns | residuals below 2e-3 <!-- claim: check5_tableVI_logRe_residual --> and below 2e-3 <!-- claim: check5_tableVI_logn0_residual -->; a single-digit slip in either moves it to 1.00 <!-- claim: check5_slip_residual_n0 --> |
| 2, cross-page | `A1.6`'s and `A1.1`'s published numbers, recomputed by independent code | Ergun at eps_MB **+58.43 <!-- claim: check4_ergun_bias_pct_at_epsMB --> %**; at the settled eps₀ **+42.06 <!-- claim: check4_ergun_bias_pct --> %**; at a textbook voidage **+4.82 <!-- claim: check4_ergun_bias_pct_at_eps042 --> %**; Wen & Yu eq. (1) **-25.21 <!-- claim: check4_wen_yu_bias_pct --> %**, all to a relative gap below 1e-6 <!-- claim: reconcile_A16_wen_yu_gap -->; `A1.1`'s refit reproduced to a relative gap below 1e-9 <!-- claim: refit_vs_A11_k1_gap --> |
| 4, digitised | none | **nothing on this page is digitised.** The one digitised dataset it loads is `A1.1`'s, used only to bound Ergun's own constants |

**The Check-3 break table is swept over voidage.** Eq. (14) is branched, so the
evaluation voidage decides which of its constants are live: at `0.50` the dilute
branch is never reached and its exponent and the switch voidage move nothing in
any of the four checks, while the page's own headline, +29.2 <!-- claim: check3_index_bias_pct_eps095 --> % at `0.95`, is
computed **on** that branch. Swept over six voidages, 7 <!-- claim: check3_readings_invisible_at_eps050 --> of fifteen injected
mis-readings that the `0.50` column cannot see become visible: -15.3 <!-- claim: check3_break_dilute_exponent_2p65_to_1p65_pts_at_eps095 --> points at `0.95` for a
dilute exponent one lower, -41.1 <!-- claim: check3_break_dilute_branch_deleted_pts_at_eps095 --> for deleting the dilute branch, +17.5 <!-- claim: check3_break_switch_0p85_to_0p55_pts_at_eps060 --> at `0.60` and
+44.5 <!-- claim: check3_break_switch_0p85_to_0p95_pts_at_eps095 --> at `0.95` for moving the switch across them. The two constants of eq. (16) are
seen nowhere (2 <!-- claim: check3_readings_never_seen --> readings): the single-sphere drag curve cancels out of an index that
compares two hindered velocities to the same terminal velocity, which Check 4 does
exercise.

## Headline results

| | |
| --- | --- |
| the two closures at eps_g = `0.44`, viscous limit | **1.192 <!-- claim: ratio_viscous_limit_eps044 -->**: they agree to 19 <!-- claim: limit_gap_viscous_pct_eps044 --> % |
| the same, inertial limit | **0.957 <!-- claim: ratio_inertial_limit_eps044 -->**: 4 <!-- claim: limit_gap_inertial_pct_eps044 --> % |
| the same, worst point in between | **0.579 <!-- claim: worst_ratio_SO_over_Ergun_at_eps044 --> at Re_m = 4.11 <!-- claim: dip_Re_m_eps044 -->**, the root of the analytic derivative, i.e. a factor **1.73 <!-- claim: disagreement_factor_eps044 -->**: on Ergun's own abscissa that is Re/(1-ε) = **7.33 <!-- claim: dip_abscissa_ergun_eps044 -->**, and his 244 <!-- claim: ergun_markers_rows --> recovered markers start at **7.45 <!-- claim: ergun_markers_abscissa_min -->**, so the worst point sits 1.6 <!-- claim: dip_below_markers_pct_eps044 --> % *below* his data. At eps_g = `0.50` the dip is at 9.49 <!-- claim: dip_abscissa_ergun_eps050 -->, inside it |
| where the interior minimum stops being the lowest point of the ratio | above eps_g = 0.53 <!-- claim: dip_stops_being_lowest_eps -->, where the inertial limit becomes lower; the ratio has an interior minimum at 4 <!-- claim: dip_voidages_with_interior_minimum --> of the six voidages tabulated and it is the lowest point at 3 <!-- claim: dip_voidages_where_minimum_is_lowest --> of them |
| Ergun drag / isolated-sphere drag as eps_s → 0 | **7.9 <!-- claim: ergun_isolated_Re1_pct --> %** at Reynolds number one, **0.10 <!-- claim: ergun_isolated_Re0p01_pct --> %** at a hundredth of it, **5.79 <!-- claim: ergun_isolated_Re1e6_ratio -->×** at Re = 10⁶, overshooting above Re = 29 <!-- claim: ergun_isolated_crossover_Re --> |
| the eps_g = `0.85` branch switch in eq. (14) | **+0.049 <!-- claim: B_branch_jump_pct_at_switch --> %** in B; at most **0.049 <!-- claim: Vrm_branch_jump_max_pct --> %** in V_rm, 0.61 <!-- claim: B_jump_in_half_units_of_dilute_exponent --> half-units of the last printed digit of the dilute exponent. **It is not a discontinuity** |
| Dalla Valle's viscous limit | C_D·Re → `23.04` against Stokes' 24, **-4.00 <!-- claim: dalla_valle_vs_stokes_pct --> %** |
| the closures cross in the viscous limit at | eps_g = 0.507 <!-- claim: cross_viscous_lower_eps --> and 0.827 <!-- claim: cross_viscous_upper_eps -->, in the inertial limit at 0.428 <!-- claim: cross_inertial_eps -->; the viscous-limit ratio has its minimum 0.822 <!-- claim: viscous_ratio_minimum --> at 0.682 <!-- claim: viscous_ratio_minimum_eps --> |
| implied index vs Richardson & Zaki, eps_g = `0.45` → `0.95` | **+3.2 <!-- claim: check3_index_bias_pct_eps045 --> % → +29.2 <!-- claim: check3_index_bias_pct_eps095 --> %**, on a 15 <!-- claim: check3_row_spread_pts_eps045 -->-point row-to-row spread at `0.45`, **under the slip reading of V_rm**. Under the superficial reading it is **-26.3 <!-- claim: check3_alt_reading_bias_pct_eps045 --> % → -0.3 <!-- claim: check3_alt_reading_bias_pct_eps095 --> %**. Row by row the implied index exceeds every tabulated one from eps_g = 0.53 <!-- claim: check3_all_rows_above_from_eps --> |
| Geldart's eight cuts, at the derived eps₀ (0.441 <!-- claim: eps0_diakon_min --> to 0.498 <!-- claim: eps0_diakon_max -->) | Ergun **+42.1 <!-- claim: check4_ergun_bias_pct --> %** (standard error 5.3 <!-- claim: check4_ergun_bias_se_pct --> %), Syamlal-O'Brien **+101.2 <!-- claim: check4_syamlal_obrien_bias_pct --> %** (17.7 <!-- claim: check4_syamlal_obrien_bias_se_pct --> %), Richardson-Zaki **+111.1 <!-- claim: check4_richardson_zaki_bias_pct --> %** (10.1 <!-- claim: check4_richardson_zaki_bias_se_pct --> %) |
| the same, at a textbook voidage (`A1.6`'s choice) | Ergun **+4.8 <!-- claim: check4_ergun_bias_pct_at_eps042 --> %** |
| voidage each closure would need to be unbiased on those rows | Ergun **0.4147 <!-- claim: check4_eps_erasing_ergun_bias -->**, Richardson-Zaki 0.3708 <!-- claim: check4_eps_erasing_richardson_zaki_bias -->, Syamlal-O'Brien 0.3538 <!-- claim: check4_eps_erasing_syamlal_obrien_bias -->: all 3 <!-- claim: check4_closures_needing_voidage_below_table --> below the lowest settled voidage of the rows, Ergun by the least |
| Syamlal-O'Brien / Ergun on those rows | **0.97 <!-- claim: ratio_SO_over_Ergun_min_diakon --> to 1.68 <!-- claim: ratio_SO_over_Ergun_max_diakon -->**, with its largest value of 1.6812 <!-- claim: ratio_SO_over_Ergun_max_diakon --> at 263 <!-- claim: ratio_SO_over_Ergun_max_row_d_sv_um --> µm and 1.6717 <!-- claim: ratio_SO_over_Ergun_last_cut_diakon --> at 318 <!-- claim: ratio_SO_over_Ergun_last_row_d_sv_um --> µm: **not** monotone |
| Syamlal-O'Brien / Richardson-Zaki on those rows | **0.74 <!-- claim: ratio_SO_over_RZ_min_diakon --> to 1.12 <!-- claim: ratio_SO_over_RZ_max_diakon -->**: an *algebraic identity*, not a measured result, since v_t cancels out of it to a relative gap below 1e-12 <!-- claim: check4_vt_cancels_out_of_SO_over_RZ --> |
| Wen & Yu's empirical u_mf correlation, same rows | **-25.2 <!-- claim: check4_wen_yu_bias_pct --> %**, and it is the one comparison with no voidage in it at all |

## What is not here

- **Nothing about Gidaspow's closure or Wen-Yu drag.** The eps_g = `0.8`
  discontinuity a comparison of blended closures exists for is neither computed
  nor asserted. The only switch measured here is Garside & Al-Dibouni's at `0.85`,
  and it turns out to be continuous.
- **The inertial half of every closure is essentially untested against
  measurement.** Every Geldart row is a group A powder at a single-particle Re_t
  below 28 <!-- claim: check4_Re_t_max_diakon --> and a bed Reynolds number below 2 <!-- claim: check4_Re_mf_max_diakon -->. Table VI constrains eqs. (12)-(14)
  only, and Check 3 is exactly blind to the single-sphere drag curve.
- **No claim that Ergun is a reference standard, and no claim that any closure is
  validated.** `A1.6` shows Ergun is +58 <!-- claim: check4_ergun_bias_pct_at_epsMB --> % biased at Geldart's reported ε_MB and
  +4.8 <!-- claim: check4_ergun_bias_pct_at_eps042 --> % at a textbook voidage; both are recomputed here rather than quoted. What
  the page *does* support is the narrow statement that Ergun has the smallest
  absolute bias at 4 <!-- claim: check4_ergun_least_abs_bias_voidages --> of five voidages tried and needs the smallest change of voidage to
  be unbiased, with the absolute bias attributed to the voidage, which nobody
  measured, rather than to the closures.
- **Neither Check-4 ratio is measurement-facing.** v_t cancels identically out of
  Syamlal-O'Brien / Richardson-Zaki, so 0.74 <!-- claim: ratio_SO_over_RZ_min_diakon --> to 1.12 <!-- claim: ratio_SO_over_RZ_max_diakon --> is Check 3 restated at Geldart's eight
  (ε₀, Re_t) pairs and would read the same if the experiment had never been done.
  Geldart's U₀ enters only the per-cent columns.
- **Which of two readings of V_rm the report meant.** Eq. (11) forces the slip
  ratio; the report's own prose introducing V_rm implies the superficial one.
  They are one unit of index apart and the alternative **reverses the sign** of
  Check 3 (-26.3 <!-- claim: check3_alt_reading_bias_pct_eps045 --> %, -24.2 <!-- claim: check3_alt_reading_bias_pct_eps050 --> %, -0.3 <!-- claim: check3_alt_reading_bias_pct_eps095 --> %, viscous limit -11.0 <!-- claim: check3_alt_reading_viscous_index_vs_RZ_pct --> %). Both are printed
  and CI-tracked; the page adopts the one eq. (11) implements and says so. Garside &
  Al-Dibouni (1977), which would settle it, is not on disk.
- **The table lookup holds Table VI's index constant outside its Reynolds range.**
  4 <!-- claim: check4_context_rows_clipped --> of the context rows fall outside it; 0 <!-- claim: check4_diakon_rows_clipped --> of the eight Diakon rows the result is taken
  on do.
- **Table VI is not a raw measurement.** Each n₀ is a fitted slope, extrapolated
  along a second fitted line to d/D = 0, from liquid-solid runs. And Garside &
  Al-Dibouni fitted their relation to overlapping literature, so the comparison
  measures the cost of one closed form standing in for another, not two
  independent witnesses.
- **The switch voidage `0.85` is the one printed constant no data comparison
  settles** when it is displaced by a little: a small displacement is worth
  0.048 <!-- claim: check3_break_switch_small_displacement_max_pts --> points at most in the Check-3 table and is bounded by the +0.049 <!-- claim: B_branch_jump_pct_at_switch --> % continuity of
  eq. (14).
- **The non-Diakon Geldart rows** carry no printed sphericity, and
  Syamlal-O'Brien has no shape factor to give them. Whole-table numbers are
  printed as context, not as a result.

## Corrections

- 2026-11-09: placed the worst point of the closure ratio from a 400-point log grid
  of Re_m: "0.579 at Re_m = 4.09", on Ergun's abscissa "Re/(1-eps) = 7.30",
  "the worst point is 2.0 % below his data" and "at eps_g = 0.50 the dip is at 9.72".
  Those were the nearest grid node. The root of the analytic derivative of the
  ratio is at Re_m = 4.11 <!-- claim: dip_Re_m_eps044 -->, Re/(1-ε) = 7.33 <!-- claim: dip_abscissa_ergun_eps044 -->, so the worst point is
  1.6 <!-- claim: dip_below_markers_pct_eps044 --> % below the markers, and the dip at eps_g = `0.50` is at 9.49 <!-- claim: dip_abscissa_ergun_eps050 -->.
- 2026-11-09: said the dip stops being an interior minimum "above eps_g ≈ 0.50; above that they simply diverge". The ratio is not monotone above that: it has an interior minimum at 4 <!-- claim: dip_voidages_with_interior_minimum --> of the six voidages tabulated, and that minimum is the lowest point only up to eps_g = 0.53 <!-- claim: dip_stops_being_lowest_eps -->, where the inertial limit becomes lower.
- 2026-11-09: said "Ergun is least biased at every voidage tried". Of the five voidage readings it has the smallest absolute mean bias at 4 <!-- claim: check4_ergun_least_abs_bias_voidages -->, each by at least 3.98 <!-- claim: check4_ergun_lead_min_over_se --> standard errors of the paired difference; at the voidage Wen and Yu's eq. (2) forces on a sphere Richardson-Zaki's is smaller by 12.7 <!-- claim: check4_rz_lead_pts_eps0383 --> points (standard error 16.9 <!-- claim: check4_rz_lead_se_pts_eps0383 -->), so no closure is least biased there.
- 2026-11-09: said "The ORDERING is stable". Syamlal-O'Brien and Richardson-Zaki are separated by more than two standard errors of their paired difference at only 1 <!-- claim: check4_so_rz_order_resolved_voidages --> of the five voidage readings, so no order of those two is claimed.
- 2026-11-09: said "5 of 10 injected mis-readings leave it at the threshold". Reading `0.06` and `0.12` both a decimal place out leaves the radicand a perfect square, so that residual was round-off, not a catch; the identity sees 4 <!-- claim: check1_readings_seen --> of ten.
- 2026-11-09: said the closed forms of the two limits reproduce the full expression
  evaluated at "$Re = 10^{{-9}}$ and $10^{{7}}$ to". The full expression had been evaluated at those values of the slip
  velocity in m/s, not of Re_m. At the stated Re_m the relative gap is 0.0047 <!-- claim: limit_closed_form_gap --> at worst (the page now evaluates at the
  stated Re_m).

## Related

[`A1.1`](../A1.1-ergun-pressure-drop/) (the packed-bed branch and every Ergun
constant used here), [`A1.6`](../A1.6-wen-yu-minimum-fluidisation/) (the
empirical u_mf correlation, and the caveats that bound Check 4),
[`A1.7`](../A1.7-geldart-classification/) (the measured dataset and the derived
voidage), `A1.5` (Richardson-Zaki in its own right), `E1.2` and `E2.1`
(fluidised-bed models that consume u_mf rather than a drag law).
