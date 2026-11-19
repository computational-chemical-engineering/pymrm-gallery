# J4.6: Michaelis-Menten and Briggs-Haldane

**One rate law, two constants, and what data can decide between them.**

Section 1.5.1 of Froment, De Wilde & Bischoff derives the enzyme rate law twice.
Assume the complex decomposition is rate determining and the binding step
equilibrates, and you get eq. (1.5.1-8) with K_M = k_-1/k_1; the line after it
names this "the Michaelis-Menten equation". Assume only a pseudo steady state on
the complex ("the second approach, formulated by Briggs and Haldane") and you get
the identical expression with K_M = (k_-1 + k_2)/k_1, eq. (1.5.1-14). The two
constants differ by the factor 1 + k_2/k_-1, and nothing in the rate law bounds
it. This page asks what an experiment can do about that.

## What it finds

**No steady-state rate measurement can separate them, and this is exact.**
Substituting k_2 = r_m/C_E0, k_-1 = k_2/rho, k_1 = (k_-1 + k_2)/K_M into
eq. (1.5.1-13) returns eq. (1.5.1-15) identically and dr/drho = 0, in sympy. The
numerical sweep is a check on the code: with r_m and K_M fixed and rho = k_2/k_-1
swept over eight decades, the rate law from the elementary constants agrees with
eq. (1.5.1-15) to below 1e-14 <!-- claim: rate_law_rho_degeneracy_max_rel -->
relative on the nine swept rho, and to below 1e-14 <!-- claim: rate_law_rho_degeneracy_random_draws_max_rel -->
over 2000 random (rho, C_E0) draws. Raising k_2 by one percent with C_E0 left
alone moves it to 9.8e-3 <!-- claim: rate_law_break_k2_plus_1pct -->.

**A batch transient sees rho, but almost entirely where nobody samples.** With
the observables fixed, C_A(t) spreads by 0.9773 <!-- claim: rho_spread_induction_over_eps -->
eps across those eight decades, eps = C_E0/(c_A0 + K_M), and the maximum sits at
t = 2.011e-7 <!-- claim: rho_spread_induction_peak_min --> min, inside the
induction layer, whose time scale 1/[k_1(c_A0 + K_M)] runs from
2.09e-8 <!-- claim: induction_time_rho_lo_min --> min at rho = 1e-4 to
2.09e-4 <!-- claim: induction_time_rho_hi_min --> min at rho = 1e4. A uniform
181-point scan of the run reads that maximum 66.3 % <!-- claim: sampled_max_shortfall_uniform_scan -->
low; the page locates it as a root of the analytic time derivative after a
logarithmic and uniform scan. Froment's sentence that "only transient
experimentation (stopped flow, or relaxation ...) can help" is right, and this is
the number behind it.

**Over the window a slow experiment can see, the requirement is severe.** For
t >= 3 min (the table's first sample after zero) only 0.242 <!-- claim: rho_spread_over_eps_trace_limit -->
eps survives. Against the table's 5e-4 cmol/L print resolution that needs
eps = 2.069e-3 <!-- claim: eps_star_rho_visible_at_print_resolution --> with every
other constant fixed, and eps = 5.798e-2 <!-- claim: eps_star_discrimination_with_compensation -->
with r_m, K_M and C_E0 free to compensate. At that loading the quasi-steady-state
approximation is itself off by 88.9 <!-- claim: qssa_error_at_discrimination_eps_over_resolution -->
times the resolution against the free substrate C_A, and
38.8 <!-- claim: qssa_error_at_discrimination_eps_total_over_resolution --> times
against the total substrate C_A + C_A-E.

**How much of that is compensation.** The two epsilons are not the same
measurement: the first is a maximum over t in [3, 18] min of the rho-spread, the
second an RMS over a 25-point grid from t = 0 of a refitted misfit. The same RMS
on the same grid with nothing free gives eps = 3.185e-3 <!-- claim: eps_star_rms_25pt_grid_no_compensation -->.
Of the raw factor 28.0 <!-- claim: eps_star_raw_ratio_unlike_statistics -->,
1.539 <!-- claim: eps_star_statistic_and_window_factor --> is the change of
statistic and window and 18.20 <!-- claim: eps_star_compensation_penalty --> is
the price of compensation.

**On the actual data, rho is free.** Rawlings & Ekerdt's Exercise 9.15 gives
r_m = 0.162583 <!-- claim: mm_fit_rm_cmol_per_L_min --> cmol/(L min) and
K_M = 0.503055 <!-- claim: mm_fit_KM_cmol_per_L --> cmol/L at
4.02e-3 <!-- claim: mm_fit_rms_cmol_per_L --> cmol/L RMS, against
3.90e-2 <!-- claim: null_first_order_rms_cmol_per_L --> for the best
first-order null. Pinning rho anywhere in [1e-4, 1e4] and re-optimising the rest
moves the achievable sum of squares by 8.6 % <!-- claim: rho_profile_sse_span_fraction -->,
monotonically, and the two elementary-step parameters give
F = 1.90 <!-- claim: full_model_F_statistic --> against
F(2, 3) = 9.55 <!-- claim: full_model_F_crit_95 --> at 95 %.

## What it does not claim

It does not validate "the Michaelis-Menten mechanism". What the seven points
support is that a two-parameter saturating law beats the best one-parameter
simple order by 9.71 <!-- claim: null_first_order_rms_ratio --> times in RMS: a
claim about saturation, compatible with every mechanism in the two-parameter
family the page constructs, which has both of the book's derivations as end
members.

## Two printed defects, reported and not repaired

Both are settled from the book's own equations and were read on 300 ppi crops.

- **p. 24**, verbatim: *"the rate levels off and becomes zero order with respect
  to the reactant, rA = k1CE0."* [sic] The limit of the book's own eq. (1.5.1-8)
  is k_2 C_E0; its eq. (1.5.1-15) calls k_2 C_E0 "the maximum possible rate"; its
  Problem 2.4(b) asks the student to show that limit; and k_1 C_E0 has units of
  1/s, which a rate cannot have.
- **p. 25**, verbatim: *"the Lineweaver-Burke plot of 1/r versus CA."* [sic] 1/r
  is exactly affine in 1/C_A (below 1e-13 <!-- claim: lineweaver_correct_axes_nonlinearity -->
  of its range off a straight line) and strictly convex in C_A. The statistic is
  a maximum, so the line that minimises it is the Chebyshev line: no straight
  line comes closer than 31.7 % <!-- claim: lineweaver_printed_axes_nonlinearity_minimax -->
  of the curve's range over 0.05 to 1.0 cmol/L, in closed form with three-point
  equioscillation. The best least-squares line is further out and depends on a
  convention: 50.3 % <!-- claim: lineweaver_printed_axes_nonlinearity --> under
  the page's logarithmic weighting, 66.9 % <!-- claim: lineweaver_printed_axes_nonlinearity_uniform_measure -->
  under a uniform one. Levenspiel p. 615 prints the reciprocal axes and no
  surname at all.

The clean text layer drops the leading term of eq. (1.5.1-17) on p. 26 and makes
a correct equation look wrong; the render shows it is fine.

## Sources

**Read from:** Froment, G. F., De Wilde, J. & Bischoff, K. B., *Chemical Reactor
Analysis and Design*, 3rd edn, Wiley (2011), ISBN 978-0-470-56541-4, section
1.5.1, book pp. 23 to 26, eqs. (1.5.1-1) to (1.5.1-15).

**Origins, cited and not consulted:** Michaelis, L. & Menten, M. L., *Biochem.
Z.* 49, 333-369 (1913); Briggs, G. E. & Haldane, J. B. S., *Biochem. J.* 19,
338-339 (1925). Neither was read. The book prints no citation for either; the
searches are described in `meta.yaml` and the dataset sidecar.

**Data:** Rawlings, J. B. & Ekerdt, J. G., *Chemical Reactor Analysis and Design
Fundamentals*, 2nd edn (2025 printing), Exercise 9.15: seven c_S-versus-time
pairs, provenance unstated by the source. The table is not an exactly-rounded
model evaluation (maximum residual 13.87 <!-- claim: mm_fit_maxres_over_print_resolution -->
times the print resolution), and K_M sits within 0.61 % <!-- claim: mm_fit_KM_rel_dev_from_half -->
of the round value 0.5 cmol/L. Nothing the page concludes depends on it.

**Levenspiel**, *Chemical Reaction Engineering*, 3rd edn (1999), p. 615, is used
for the reciprocal plot's axes and for his claim that the direct
C_A-versus-tau fit "is direct, is less prone to fiddling, and is more reliable".
Measured here: the reciprocal construction run as the exercise instructs returns
K_M -46.7 % <!-- claim: lb_KM_bias_vs_nls --> and r_m -34.1 % <!-- claim: lb_rm_bias_vs_nls -->
off the direct fit, of which 14.0 % <!-- claim: lb_KM_bias_finite_difference_only -->
upward on K_M is the central-difference rate estimate alone.

## Files

| file | what it is |
|---|---|
| `index.py` / `index.ipynb` | the page, paired by jupytext |
| `model.py` | the closed form, the scipy reference, `max_gap`, and the pymrm marchers `BatchMassAction` and `BatchReduced` |
| `test_model.py` | fast checks of `model.py` and one stored claim |
| `data/rawlings-ekerdt-ex9.15-batch.csv` | the seven printed points, with a sidecar that leads on the unstated provenance |
| `data/froment-1.5.1-printed.csv` | the first fifteen of the section's eighteen numbered equations and nine prose claims, verbatim, defects flagged and not repaired |
| `agreement.json` | the claims in format 2, the prose literals and the break table |

## Validation summary

- ten symbolic identities close eqs. (1.5.1-2) to (1.5.1-15), both routes;
- two conservation balances in the pymrm mass-action marcher, both below
  1e-13 <!-- claim: pymrm_substrate_balance_max --> cmol/L, on a model that
  integrates all four fields independently;
- backward Euler first order, the observed order farthest from 1 over successive
  differences 0.996 <!-- claim: pymrm_be_time_order -->, Richardson against Radau
  4.19e-7 <!-- claim: pymrm_richardson_vs_radau --> cmol/L, LSODA against Radau
  below 1e-9 <!-- claim: lsoda_vs_radau_reference --> cmol/L;
- two second routes: the fit through the pymrm marcher instead of the Lambert W
  closed form (k and K agree to below 1e-5 <!-- claim: fit_second_route_k_rel_diff -->
  and below 1e-6 <!-- claim: fit_second_route_K_rel_diff --> relative), and rho*
  from an analytic sensitivity with its maximum root-found, 2.1934e-3 <!-- claim: rho_star_equilibrium_analytic -->
  against 2.1951e-3 <!-- claim: rho_star_equilibrium_route -->,
  0.080 % <!-- claim: rho_star_equilibrium_two_route_rel_diff --> apart;
- a defect-injection table: every non-structural claim fails at least one row,
  and each row injects its defect through an argument of the study that
  computes the claim.

## Reuse

Fit eq. (1.5.1-15), or its integral for a batch run, directly and report
(r_m, K_M). Do not report k_-1/k_1 or (k_-1 + k_2)/k_1 as though the data
distinguished them. If you linearise, the abscissa is 1/C_A; and consider not
linearising. If you check the quasi-steady state, the group is C_E0/(c_A0 + K_M):
comparing C_E0 with the substrate alone mispredicts the error by a factor
95.0 <!-- claim: qssa_group_spread_CE0_over_cA0 --> over the K_M swept here,
against 1.91 <!-- claim: qssa_group_spread_CE0_over_cA0_plus_KM --> for the
group. To put this rate law inside transport, take `BatchReduced.source` to
`J4.7` or `B1.1` unchanged.
