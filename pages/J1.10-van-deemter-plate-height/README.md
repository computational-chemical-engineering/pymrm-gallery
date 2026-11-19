# J1.10: The van Deemter equation

**Catalog ID:** `J1.10` · **Section:** J · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S4` (1-D transient PDE) · **Data tier:** 6 (the source's own
four printed tables, transcribed; **no figure digitised**)

Reproduces J. J. van Deemter, F. J. Zuiderweg and A. Klinkenberg, "Longitudinal
diffusion and resistance to mass transfer as causes of nonideality in
chromatography", *Chemical Engineering Science* **5**(6) 271-289 (1956),
doi:10.1016/0009-2509(56)80003-1, the paper behind $H = A + B/u + Cu$.

Built from the native-resolution scan (confirmed with `pdfimages -list`). **The
text layer was used for nothing numeric**: it misreads the journal *volume* in the
document's very first line, extracting `Vol 6` where the page prints `Vol. 5`.
Every table cell, equation constant and quoted glyph was read on a digit-scale
crop of a render at the scan's resolution.

**The paper prints four tables, not two.** Table 1 (p. 277) is the feed-volume
experiment, Table 2 (p. 280) is Simpson & Wheaton's ion-exclusion data restated,
Table 3 (p. 281) is the plate heights derived from Table 2, and Table 4 (p. 286)
is the gas-liquid measurement that produces the three-term formula. All four are
transcribed; the two that carry the plate height are 3 and 4.

**No point was taken off any figure.** The ambiguous readings were settled by
arithmetic constraints the tables impose on themselves (the doubling series of
Table 2's velocities, the constant column length $nH$ and the particle diameter
$H/(H/d_p)$ of Table 3); the readings in Table 4 that Figs. 11 and 13 were
consulted for also follow from inside Table 4 alone, from its ascending velocity
column and from eq. (55), and the page shows it. Two printed *curve labels* on
Fig. 8 are read as text, and used only to identify which of Table 3's two printed
endpoints came from which route.

## Fit, not test, where it matters

**Eqs. (54) and (55) were fitted by the authors to the twelve Table 4 points
this page compares them against**, and the paper says so: "the values of the
constants in (54) are highly uncertain, owing to the relatively small number of
experimental data that were available". Every residual against Table 4 here is
a *goodness of fit*, is called one, and is reported beside a null baseline (the
best constant $H$, which every fit beats).

The same applies to $\lambda$ and $D_{II}$: the authors extracted them *from*
Table 3, so refitting Table 3 reproduces their extraction, not an independent
measurement of either constant.

## What the page shows

- **The 1956 algebra is sound, all four steps.** (34) to (38), (38) + (45, 46, 48)
  to (49) and its $C$, (38) + (41, 50, 51) to (52): each residual simplifies to
  exactly zero. The paper names the route for each but never prints the
  manipulation. A route the paper never takes, the Laplace transform of eqs.
  (29)/(30), whose residence-time cumulants give $H = \kappa_2/\kappa_1^2$,
  reproduces eq. (38) symbolically as well.

- **Eq. (38), computed without eq. (38).** Integrating eqs. (29)/(30) in time
  with pymrm and reading $H$ from the difference of the first two time moments
  at two interior stations agrees to 4.1e-5 <!-- claim: pymrm_H_vs_eq38_rel -->.
  The two-station difference cancels every boundary effect: swapping a Danckwerts
  inlet for a Dirichlet one moves $H$ by below 1e-5 <!-- claim: pymrm_bc_insensitivity -->
  relative (a bound, because it is a difference of two nearly equal plate
  heights), while the same swap read from one station is not cancelled.

- **A refinement study with a known answer.** First-order upwind adds $uh/2$ of
  numerical dispersion, so switching the limiter off must shift $H$ by *exactly
  $h$*. Measured coefficient: 0.99999 <!-- claim: pymrm_upwind_dispersion_coeff -->.
  With the van Leer correction the observed order of $H$, from successive
  differences on four nested grids, is at least 2.5 <!-- claim: pymrm_grid_order -->
  and at most 3.5 <!-- claim: pymrm_grid_order -->; bare upwind is first order.

- **The time march, measured.** The step is linear apart from the van Leer
  limiter, which the march takes as a correction deferred to the start of the
  step: one factorised matrix for the whole run. The same theta scheme with the
  limited flux implicit, solved by `newton` with the exact Jacobian of
  `ConvFluxTVD`, gives the same two-station $H$ to
  below 1e-5 <!-- claim: pymrm_deferred_vs_newton_H_rel -->. Over a 16-fold change of the
  step, $H$ moves below 1e-6 <!-- claim: pymrm_dt_swing_H -->; with the limiter
  off, Crank-Nicolson's order on the concentration at a fixed time departs from 2
  by below 0.1 <!-- claim: pymrm_dt_order_upwind_dev -->, and backward Euler is
  first order.

- **$D_{II}$ and $\lambda$ recovered from Table 3.** Refitting eq. (49) on the
  five varying-flow rows gives $\lambda$ = 2.99 <!-- claim: lambda_from_table3 -->
  against the printed 3 (0.24 % <!-- claim: lambda_from_table3_rel_dev -->) and
  $D_{II}$ = 1.34e-10 <!-- claim: DII_from_table3_m2s --> m² s⁻¹ against the
  printed "$1.3\times10^{-10}$" (3.7 % <!-- claim: DII_from_table3_rel_dev -->). The
  extraction is **almost independent of $D_I$**, which the paper never prints
  for that system: it moves only 1.8 % <!-- claim: DII_D_I_sensitivity --> as
  $D_I$ is swept twofold and then removed altogether.

- **The film thickness, three times.** 10.06 <!-- claim: df_celite_n_butane_um -->,
  10.57 <!-- claim: df_celite_iso_butane_um --> and
  8.66 <!-- claim: df_sterchamol_um --> µm from three different printed constant sets,
  against the paper's $\approx 10$, $\approx 10$ and 9 µm. Eq. (56)'s coefficient
  comes back at 3.333e-8 <!-- claim: eq56_coefficient --> against the printed
  "$3.3\times10^{-8}$".

- **Table 1's `calc` column, every value, and what it really pins.** The column
  depends on one number, $v\sqrt n = \Delta S_0/4$, and the printed values pin
  $\Delta S_0$ to the window from 166.0 <!-- claim: table1_dS0_window_lo_ml --> to
  167.4 <!-- claim: table1_dS0_window_hi_ml --> ml (edges in closed form from
  eq. (20)), above the paper's own "$\sim$ 165 ml", which is written with a tilde
  and is therefore not a defect, merely less precise than its own consequences.

- **Table 3 rebuilt from Table 2 by the paper's two routes.** Median deviation
  2.1 % <!-- claim: table3_n_median_rel_dev --> over the printed endpoints; split
  **by route**, the height route (median 1.6 % <!-- claim: table3_n_median_height_route -->)
  beats the width route (2.2 % <!-- claim: table3_n_median_width_route -->), as it
  must, because Table 2's widths are printed as integers and $n$ goes as $1/w^2$.
  The printed pair is in ascending order, so route and endpoint are not the same
  split: Fig. 8's own printed labels fix low = height for the 50-100 mesh row, and
  the page prints the rows where the two routes invert.

- **A printed equation that cannot be right.** Eq. (21) as printed carries
  $(a+\delta)/\delta$ where the derivation and the paper's *own* eqs. (24) and
  (25) require $(a+\delta)/a$: as printed it tends to 2
  <!-- claim: eq21_printed_a0_limit --> in the limit eq. (24) fixes at 4, and
  diverges where eq. (25) fixes it at $a+\sqrt{2\pi}$. Eq. (22) is not free, it is
  exactly the inflection condition of eq. (19), verified symbolically, so the
  denominator is the only thing left, and deriving $\Delta s$ from eq. (19)'s own
  inflection point returns the printed term with $a$ underneath at residual
  exactly zero. As printed, the width route to the plate number has a root for
  only 2 <!-- claim: eq21_printed_solvable_rows --> of Table 2's nine conditions.
  **Reported, not repaired**: quoted verbatim, shown to fail its own two limits,
  and the corrected form used with every dependent number labelled.

- **The velocity convention, stated rather than assumed.** The List of Symbols
  (p. 288) defines $u$ as *interstitial*; Tables 2 to 4 and Figs. 11 to 14 are
  headed *superficial*. Fig. 9 converts (abscissa $uF_Id_p$); the gas-liquid
  section does not appear to: the quoted $2\gamma D_I = 8\times10^{-6}$ m² s⁻¹ is
  numerically identical to eq. (54)'s iso-butane $B$ read in mm² s⁻¹, i.e. with
  $F_I = 1$, and implies $\gamma$ = 0.444 <!-- claim: gamma_iso_butane_superficial -->,
  below the paper's own stated range. $F_I$ is not printed for those columns, so
  **the ambiguity cannot be closed**; the page prints the band of $F_I$ each
  reading needs and notes that $2\lambda d_p$ is immune because the eddy term
  contains no $u$.

## What pymrm adds

Not the plate height: eq. (38) is exact for this model and the page proves it
twice. What it adds is a number for the paper's own *unquantified* condition.
Eqs. (33)/(34) hold "for locations $z$ much larger than both $2D/u$ and
$F_Iu/\alpha$", and how much larger is never said.

Because the cumulants of a linear column are exactly additive in $z$, $\kappa_3$
and $\kappa_4$ per unit length give the skewness and excess kurtosis of a column
of any length in closed form. **Both of the paper's groups fall as $1/z$; the
skewness falls only as $z^{-1/2}$.** At the ion-exclusion column that produced
Table 3 the groups are 2.4e-3 <!-- claim: group1_at_recovered_column --> and
3.4e-3 <!-- claim: group2_at_recovered_column -->, comfortably "$\ll 1$", and the
elution curve is still skewed by 0.14 <!-- claim: skewness_at_recovered_column -->.
Reaching a skewness of 0.05 needs a column 8.2 <!-- claim: z_for_skewness_over_column -->
times longer (in closed form, and checked by a root find).

All four cumulants come from the Laplace expansion, so every quoted shape number
is analytic; the pymrm run is the *check* on them, and it is reported with its
caveat. The PDE reproduces $\kappa_3$ to 0.21 % <!-- claim: pymrm_kappa3_vs_laplace_rel -->
and $\kappa_4$ to 1.5 % <!-- claim: pymrm_kappa4_vs_laplace_rel --> at the
reference settings, but neither is a converged discretisation error: the grid
error on the higher cumulants is positive and the time error negative, and the
page prints the $\Delta t$ table in which both deviations change sign. $H$ itself
has no such problem: it is grid-limited and barely moves with $\Delta t$.

## What the page cannot conclude

- **Nothing about whether eq. (53) is right.** Its constants were fitted to the
  points it is compared with. The nearest thing to a test the paper offers is the
  slope ratio, 1.448 <!-- claim: slope_ratio_theory --> predicted from the two
  printed $KF_I/F_{II}$ against the paper's "about 1.5" and eq. (54)'s "1.6",
  and that is one number.
- **Nothing about $\gamma$ or $d_f$ as physical quantities.** They are the
  paper's own inversions, redone, and the velocity convention moves them in
  opposite directions: $\gamma$ up by $1/F_I$, $d_f^2$ down by $F_I$, each by up
  to the factor the paper's own range of $\gamma$ allows (printed on the page).
  $\lambda$ is untouched by it: the eddy term carries no $u$.
- **Nothing about Simpson & Wheaton's measurements.** Their paper is not on disk
  and was not consulted; Table 2 is van Deemter et al.'s restatement.
- **Nothing about the ion-exclusion column's geometry.** Its length is printed
  nowhere; the 558 <!-- claim: table3_recovered_length_mm --> mm used here is
  recovered from Table 3 and is an inference.

## Files

```
index.py       the page source (jupytext percent format), paired with index.ipynb
index.ipynb    the executed page
model.py       Column (eqs. 29/30), the moment readings, eqs. (20) to (22) and the fits
test_model.py  fast checks of model.py
meta.yaml      metadata, printed defects, inferences, fit-vs-test labelling
data/          four tables + the printed scalars, each with a provenance sidecar
```

**Cite the source, not this page.**

## Corrections

- 2026-11-19: said that $H$ "moves 1.2 × 10⁻⁶ over 8× in $\Delta t$". The van Leer
  correction of that march read the inlet data of $c_{\text{in}} = 1$ at every
  step instead of the current inlet value (reproduced with the old code); with
  the current inlet value $H$ moves below 1e-6 <!-- claim: pymrm_dt_swing_H -->
  relative over a 16-fold range of the step.
- 2026-11-19: said that the time order of $\mu_1(z_2)$ was claimed "between 1.5 and 2.5" and that the orders "fall from about 2 at the coarse end to about 1". Both came from the same stale inlet data. With the current inlet value the successive differences of $\mu_1(z_2)$ change sign 1 <!-- claim: pymrm_dt_mu1_sign_changes --> time on the ladder, so the limited march has no regular order there, and the page claims a time order only with the limiter off.
