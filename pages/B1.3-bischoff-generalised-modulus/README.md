# B1.3 - Bischoff's generalised modulus: measuring how narrow the narrow band is

Bischoff (1965) builds the effectiveness-factor modulus from the integral of
the rate, which makes the large-modulus ends of the curves for *every*
reaction rate form coincide exactly, and claims the rest of the curves stay
within "about 15%" (orders one-half to three) or "about 30%" (zero order
included). This page measures those claims, computes the LHHW and
volume-change collapses the paper only redrew or tabulated, and audits the
construction out of sample on a finite cylinder against a 1998 competitor.

- **Structure:** `S3` (1D steady BVP; one 2-D audit)
- **Reference:** Bischoff (1965) AIChE J 11(2) 351-355,
  doi:10.1002/aic.690110229; secondary comparison Pan & Zhu (1998)
  Chem Eng Sci 53(5) 933-946, doi:10.1016/S0009-2509(97)00385-0
- **Source:** `index.py` (paired with `index.ipynb`); the rate forms, the three
  routes to eta and the cylinder in `model.py`, tested by `test_model.py`
- **Runtime:** about 70 s

## Results

Both printed claims hold, each under one reading, and the page pins which:
14.61 % <!-- claim: band_nhalf_to_3_max_relspread --> at
m = 1.270 <!-- claim: band_nhalf_to_3_argmax_m --> for orders one-half to three
(his "about 15%", as eta_max/eta_min - 1), and
40.55 % <!-- claim: band_incl_n0_at_kink_relspread --> max/min =
28.85 % <!-- claim: band_incl_n0_frac_of_upper --> of the upper curve at the
zero-order kink m = 1 exactly (his "about 30%" is the fraction-of-upper
reading). Extrema are root-found, and both headline bands are recomputed on
the closed forms, agreeing below 1e-9 <!-- claim: band_second_route_gap -->.
On the standard modulus the same curves peak at
120.6 % <!-- claim: std_modulus_band_incl_n0 --> spread and never collapse
(100 % <!-- claim: std_modulus_asympt_spread_incl_n0 --> at large modulus with
zero order, 63.3 % <!-- claim: std_modulus_asympt_spread_nhalf_to_3 --> without).
The LHHW family's band is 30.16 % <!-- claim: lhhw_band_max_relspread --> as
max/min and 23.17 % <!-- claim: lhhw_band_frac_of_upper --> as a fraction of the
upper curve, between the two simple-order bands in either reading; the
volume-change family stays within
3.18 % <!-- claim: volchange_dev_from_n1_curve --> of the first-order curve.
A dead-zone onset formula the paper does not print, m* = (n+1)/(1-n), is
derived from its eq. (14).

pymrm: the slab solve matches the quadrature of Bischoff's parametric solution
to 2.0e-6 <!-- claim: pymrm_slab_worst_rel --> on 400 cells over eight rate
forms, including a concentration-dependent diffusivity solved by one Newton
solve on the full residual. The dead zone is solved with a regularised
rate and its exact Jacobian and extrapolated in eps and in the grid to the
exact eta = 1/m, the gap below 1e-7 <!-- claim: deadzone_extrapolated_dev -->.
Without the regularisation, at tol 1e-12 with stagnation=True, newton does not
converge on any of the four grids tried; at newton's default tol the monotone
solve reports success on
4 <!-- claim: deadzone_plain_default_tol_converged --> of the four, on profiles
that go negative in the dead zone, so that success is not a converged solution.

Out of sample (finite cylinder Z = 1, second order, reference = 2-D pymrm):
the shape's own first-order curve read at Bischoff's generalised modulus
predicts eta to 0.030 <!-- claim: aris_bischoff_combined_cyl_n2_max_abs -->
absolute (4.4 % <!-- claim: aris_bischoff_combined_cyl_n2_max_rel --> relative,
inside the advertised band); reading Figure 1's slab curve costs
0.099 <!-- claim: slab_curve_estimate_cyl_n2_max_abs --> (the B1.2 shape band);
Pan & Zhu's 1998 "arbitrary kinetics, < 1.5%" polynomial misses by
0.090 <!-- claim: pan_zhu_n2_out_of_sample_max_abs -->,
6 <!-- claim: pan_zhu_n2_over_claimed_bound --> times its claimed bound.

## Printed defects found (reported, never repaired)

- **Bischoff Table 1, omega*Co = -0.25 / eq. (33):** prints 0.905; the
  closed-form eq. (33) gives 0.9107 <!-- claim: table1_eq33_at_wco_minus_0p25 -->,
  a gap of 0.0057 <!-- claim: table1_eq33_defect_cell_abs --> absolute,
  11 <!-- claim: table1_eq33_defect_halfunit_multiple --> times the half-unit
  rounding radius of the printed three-decimal cell. The other seven
  cells pin the transcription to within printed rounding.
- **Bischoff eqs. (27)-(28) (n = 1/2 elliptic reduction) as printed** are
  inconsistent with their own source eq. (20) by
  29 % <!-- claim: eq27_as_printed_min_rel_m --> to
  41 % <!-- claim: eq27_as_printed_worst_rel_m --> in m: the where-block is
  reused from the n = 2 case, but the y = sqrt(t) substitution changes the
  prefactor exponent and the phi argument. Amended in those two places
  (labelled inference) the formula agrees below
  1e-12 <!-- claim: eq27_amended_worst_rel_m -->, and the paper's own printed
  15% shows Figure 1 used the amended reduction, since the as-printed curve
  tops out at eta = 1/sqrt(2) and would give a band of about
  41 % <!-- claim: band_asprinted_nhalf_would_give --> at the kinetic end.
- **Bischoff eq. (14), denominator bracket exponent:** prints -1/2 where its
  own parents eqs. (8)+(13) force +1/2 (the numerator's inner bracket
  correctly prints -1/2, isolating the defect to one exponent). As printed
  the equation is dimensionally inconsistent and, in normalised variables,
  low by exactly the factor 1/(n+1), against the paper's own printed eq. (21).
  The page's Route B implements the forced +1/2 reading, labelled as such.
- **Pan & Zhu Table 5, mu2 at Z = 2.5:** prints 2.50; their eq. (30)
  gives 2.40 exactly, a gap of 0.10 <!-- claim: pan_mu2_z25_defect_abs -->
  (likely a copy of the Z = 2.0 cell).
- **Pan & Zhu eq. (36):** prints q3 = -54 - 405 mu1 Q1 + 9 mu2 Q2; the
  matching conditions they state force +405, and the printed sign makes
  their own Table 6 eta_a column negative.

## Data

**Provenance tier 6: nothing here is a measurement.** All three transcribed
tables are the source authors' own computed values, read from digit-scale
crops of native-resolution renders (300 ppi Bischoff, 600 ppi Pan & Zhu).
The exact eta(m) curves computed by the page are exported as
`eta-generalised-modulus.csv`.

## Corrections

- 2026-11-01: said that on the standard modulus the five curves with zero
  order keep a spread of sqrt(2) - 1 = 41.4 % at large modulus (here, in the
  page metadata and in the notebook's Reuse section). The ratio of the
  asymptotes is sqrt((3+1)/(0+1)) = 2, so the spread tends to
  100 % <!-- claim: std_modulus_asympt_spread_incl_n0 -->, as the curves
  computed at M = 1000 confirm.
- 2026-11-01: compared the LHHW band, 30.16 % as eta_max/eta_min - 1, with
  Bischoff's "about 30%", which this page reads as a fraction of the upper
  curve (28.85 % <!-- claim: band_incl_n0_frac_of_upper -->), and concluded that "again about the same" holds against
  that figure. Read the same way, the LHHW band is
  23.17 % <!-- claim: lhhw_band_frac_of_upper -->; in either reading it lies
  between the two simple-order bands, nearer the one that includes zero order.
- 2026-11-01: said the volume-change curves stay at least 10 % of the local
  band width from either edge of the Figure-1 envelope on m from 0.2 to 8.
  That is the local minimum at the kink m = 1; above m of about 2 the curves
  approach the upper edge, eta = 1/m, and at m = 8 they are within
  3.1e-6 <!-- claim: volchange_margin_at_m8_frac_of_band --> of the band width
  from it. They stay inside the envelope.
- 2026-11-01: said every extremum on the page was root-found. The first-order
  deviations on the cylinder were the largest of eight sampled moduli; root-found,
  the slab-curve reading errs by
  0.108 <!-- claim: slab_curve_estimate_cyl_1st_max_abs -->, not 0.107.
- 2026-11-01: said the as-printed n = 1/2 curve implies a band of 41.5 % at
  m = 0.05. That was 1/eta - 1, taking the upper curve as 1; with the other
  three orders the band there is
  41.3 % <!-- claim: band_asprinted_nhalf_would_give -->.
- 2026-11-01: the data sidecar of `eta-generalised-modulus.csv` gave its error
  as about 1e-10. The centre concentration was root-found to an absolute
  tolerance of 1e-12, which at large modulus (centre concentrations below
  1e-9) resolves nothing, and the LHHW and volume-change rate integrals lost
  digits to cancellation: rows were off by up to 2.3e-5 relative (LHHW, zeta* = 0.1,
  m = 5.8). The file is regenerated; its power-law rows differ from the closed
  forms by below 1e-9 <!-- claim: export_power_rows_vs_closed_forms -->.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        rate forms, the three routes to eta, the 2-D cylinder, Pan & Zhu
test_model.py   fast checks of model.py and one stored claim
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
data/           the printed tables and the exported curves, each with a sidecar
```
