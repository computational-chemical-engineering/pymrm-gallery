# B2.3: Levenspiel's four deactivation classes

**Catalog ID:** `B2.3` · **Section:** B · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1`, `S5` · **Data tier:** 6 (model-defining equations
transcribed from the source; **the source publishes no data**)

Reproduces Octave Levenspiel, "Experimental Search for a Simple Rate Equation
to Describe Deactivating Porous Catalyst Particles", *Journal of Catalysis*
**25**(2) 265–272 (1972), doi:10.1016/0021-9517(72)90227-8, the paper that
sets out the four broad classes of catalyst deactivation (parallel, series,
side-by-side, independent) as $n$th-order kinetics with $n'$-order, $d$-order
decay, and asks *"whether experiments can be devised so as to give simply the
orders and rate constants of these equations"*.

**This is Levenspiel the author, not Levenspiel the book.** *Chemical Reaction
Engineering* is a different work (this paper cites its 2nd edition, Chap. 15,
as its own ref. 3) and remains outstanding as a textbook-canonical source; it
is what `A2.4` is catalogued to. Nothing here closes that.

Built from the 300 ppi CCITT-G4 scan. The OCR text layer is badly damaged
("detrrmination", "ordrrs", "particlrs", "parallrl", "rractor") and was used
for nothing: every symbol, subscript, prime and exponent was read from a
cropped native-resolution render at glyph scale, which matters because the
paper turns on telling $n$ from $n'$, $k$ from $k'$ from $k_d$ from $k'_d$, and
$C_{\rm A}$ from $C_{\rm R}$ from $C_{\rm P}$.

**There is no data in this source, and none is used.** All eight pages were
rendered at native resolution and read: there is no table anywhere, and
Figs. 2–6 are schematic test plots with unlabelled axes and no tick values
(Fig. 6's own caption: "A straight line on this plot indicates that the guessed
order of deactivation is correct"). No measured order is compared against,
because none is reported.

## What the page shows

- **All four batch-solids devices are exactly degenerate somewhere: five
  identity families, and each collapse is algebra, not a numerical near-miss.**
  The sharpest is the constant-flow mixed reactor: its performance relation
  makes $C_{\rm R} = \kappa\,C_{\rm A}^{\,n}\,a$ *exactly* at every instant, so
  **series deactivation of order $(n', d)$ is identical to parallel
  deactivation of order $(n n',\, d+n')$**, same observable, every time, any
  $n$. Over nine $(n, n', d)$ combinations the two differ by below
  1e-8 <!-- claim: degen_mixed_series_is_parallel --> at the production
  integrator tolerance (rtol 1e-9), and by below
  1e-11 <!-- claim: degen_mixed_series_is_parallel_rtol1e-13 --> at rtol 1e-13:
  the residual is the integrator's. The held-$C_{\rm A}$ device makes parallel,
  series and independent literally the same equation (below
  1e-12 <!-- claim: degen_constC_par_ser_ind -->, round-off), and at $d = 1$
  side-by-side joins them. Wherever the poison profile equals the reactant
  profile ($d = 1$, $n' = n$, $\sigma = \kappa/\lambda$), side-by-side **is**
  parallel with the same constants, in the batch, in the mixed reactor *and in
  the plug-flow bed* (below 1e-12 on the condition in each, 1.66e-2
  <!-- claim: degen_batch_sbs_off_condition --> off it). And with $\sigma$ left
  free, at $n = n' = d = 1$ the poison profile is a *power* of the reactant
  profile, so **side-by-side with any $\sigma$ is parallel with
  $n' \to \lambda\sigma/\kappa$**: below 1e-7
  <!-- claim: degen_sbs_is_parallel_famB_batch --> in the batch and below
  1e-8 <!-- claim: degen_sbs_is_parallel_famB_plug --> in the bed (integrator
  tolerances) over $\sigma$ = 0.2–3, against 1.2e-1
  <!-- claim: degen_sbs_is_parallel_famB_mixed --> and 3.5e-1
  <!-- claim: degen_sbs_is_parallel_famB_constC --> in the two well-mixed
  devices, where it correctly fails. $n' = 0$ collapses all four everywhere
  (0.0 exactly), which is the paper's own statement. **The bed is therefore
  *not* "the only device that is not exactly degenerate".** What the bed and
  the batch keep, and the two well-mixed devices do not, is a distinction
  between parallel, series and independent; what the bed alone has is a
  spatial coordinate.
- **Separability measured, not asserted: 64 fits, and the answer is
  negative.** Every class is fitted to every other class's response with all
  its parameters free ($\kappa$, $\lambda$, $d$, $n'$, and $\sigma$ for
  side-by-side), and the residual it cannot remove is the conversion resolution
  an experiment would need. At 0.2 conversion percentage points, **one of the
  twelve class pairs is separable in one of the four devices** (plug flow,
  parallel vs independent, 2.20e-3 <!-- claim: sep_plug_par_vs_ind -->).
  Parallel vs series in a bed is 8.07e-4 <!-- claim: sep_plug_par_vs_ser -->:
  real, but a tenth of a percentage point. **44 of the 64 entries have an exact
  identity behind them and the matrix reports min(fit, identity twin), not the
  fit**: on one of those (plug flow, independent truth, series candidate) the
  multistart alone falls into a wrong basin 1.31e-3
  <!-- claim: sep_optimiser_worst_overstatement --> short of a minimum that is
  exactly zero. The error is in the safe direction and moves no verdict, but
  printed as a fit it would have read as a measured separation.
- **The diagonal checks the fitter, and a profile measures what eq. (6) can
  return.** A class fitted to its own response has exact minimum zero. For
  parallel, series and independent the seed grid contains the truth, so those
  diagonal fits only evaluate it; the side-by-side fit starts off the truth and
  gets below 1e-9 <!-- claim: sep_stall_worst --> on every device, six orders
  below the instrument resolution. That the five constants of eq. (6) trade
  against each other is measured with a profile: holding $d$ anywhere in
  [0.5, 1.5] in either device, or $n'$ anywhere in [0.5, 1.5] (mixed flow), the
  other four refit to an RMS residual below 2e-4
  <!-- claim: sep_profile_d_mixed -->, a tenth of the conversion resolution,
  while +1 % on any single constant with the others held at the truth costs at
  least 3.9e-4 <!-- claim: sep_stall_1pct_cheapest_batch -->.
- **The discriminating measurement is the profile, not the effluent, and both
  inlet values are exact.** At the feed face $c = 1$, so $C_{\rm R} = 0$: the
  series inlet **never decays at all** ($a \equiv 1$), and the parallel inlet
  is $a = e^{-\lambda\theta}$ = 0.049787 <!-- claim: bed_a_inlet_end_par -->.
  The finite-volume first cell centre approaches each at first order, so the
  200-cell cell-centre reading is not reported as the inlet activity. The
  outlet activity is read on the outlet face (second order in space, observed
  orders 1.9985 <!-- claim: bed_a_outlet_space_order_par --> parallel and
  2.0017 <!-- claim: bed_a_outlet_space_order_ser --> series between the last
  grids) and extrapolated twice, in space and then in time,
  because the space extrapolation at a fixed step count is still time-limited
  by 9.8e-6 <!-- claim: bed_a_outlet_time_bias_par --> (parallel). The doubly
  extrapolated outlets, 0.19017 <!-- claim: bed_a_outlet_end_par --> and
  0.17711 <!-- claim: bed_a_outlet_end_ser -->, meet the quadrature route,
  which has a node at the outlet and shares no operator with the bed, to below
  1e-7 <!-- claim: bed_a_outlet_vs_quadrature_par -->, so the gradient
  reverses sign (+0.1404 <!-- claim: bed_profile_contrast_par --> vs
  -0.8229 <!-- claim: bed_profile_contrast_ser -->). One assay of the spent
  catalyst reads that off, where the outlet history is a whole campaign at a
  resolution nobody has. This is `B2.2`'s "the coke profile names the fouling
  mechanism", reached from the opposite direction, with `B2.2`'s own numbers
  carried (loaded from its claims): its *within-bed* parallel coke contrast at
  its operating point needs assays better than 10 %, a different quantity from
  the inlet activity used here, and the page says so.
- **The mixed-reactor degeneracy has an analytic boundary, and the measurement
  shows what it is worth.** The parallel-to-series relabelling needs
  $d_s = d - n'/n \ge 0$, and $d$ is an *order*, so below $d = n'/n$ no exact
  twin exists. Measured, with the twin evaluated directly where it exists: at
  or above the line the best impersonation is below 1e-8
  <!-- claim: mixed_boundary_worst_inside -->; below it the residual grows
  smoothly: below 1e-3 <!-- claim: mixed_boundary_np1_d0.75 --> at
  $d - n'/n = -0.25$, 4.0e-3 <!-- claim: mixed_boundary_np1_d0.25 --> at
  −0.75, and it crosses 2e-3 between −0.25 and −0.5. The line is a prediction
  about *exactness* and it holds; it is not a prediction about measurability,
  and the $(n' = 2, d = 1.5)$ row (half a unit below the line, still fitted to
  below 1e-4 <!-- claim: mixed_boundary_np2_d1.5 -->) shows why.
- **Levenspiel's recommended device does exactly what he says, and the page
  prices what that costs.** Holding $C_{\rm A}$ by lowering the flow rate gives
  $d$ cleanly from the shape of $\tau'(t)$ (eqs. (28)–(31) reproduce to
  round-off at $d$ = 0, 1, 2, 3) but *because* it decouples, a single run
  cannot name the class at all. Recovering the class needs a set of runs, and
  the span is root-found: with $k'_d$ known to 5 % the held levels must span
  0.204 <!-- claim: constC_span_needed_prec5pct --> in $C_{\rm A}/C_{\rm A0}$;
  at 2 %, 0.111 <!-- claim: constC_span_needed_prec2pct -->. **That 0.204 is an
  order, not a design constant**: it depends on two choices the paper never
  prints (where the level set starts and how many levels it has), and
  root-finding at 5 % over $C_{\rm LO}$ = 0.05–0.40 and 3–9 levels gives
  0.073 <!-- claim: constC_span_prec5pct_min --> to
  0.382 <!-- claim: constC_span_prec5pct_max -->, monotone in the starting
  level.
- **One knob, two constraints, with an exact blind spot at $d = 1$.** With a
  poison that is consumed, holding $C_{\rm A}$ does not hold $C_{\rm P}$. But
  the consumption rate goes as $a^{d-1}$, so at $d = 1$ it is *constant* and
  the Fig.-6 line stays straight however much poison is consumed (bias below
  1e-7 <!-- claim: constC_sbs_bias_at_d1_max --> for every $\sigma$ up to 20).
  At $d = 2$ the drift is real: the recovered order is
  1.383 <!-- claim: constC_sbs_recovered_d_base --> instead of 2 at the base
  $\sigma$ = 0.8, 5 % bias is reached at $\sigma$ =
  0.076 <!-- claim: constC_sbs_sigma_for_bias -->, and the worst over the
  scanned range is 49 % <!-- claim: constC_sbs_bias_max_scanned -->. **But the
  protocol priced here is not the one the paper prescribes.** Levenspiel's
  instruction for eq. (6) is "keep $C_{\rm P}$ constant", with $C_{\rm A}$ the
  conditional extra (p. 271; both sentences are rows of the claims CSV and both
  are quoted in full on the page). Under his own protocol (the $\sigma$ = 0
  row, poison in excess so $C_{\rm P}$ never moves) the order comes back
  unbiased to below 1e-7 <!-- claim: constC_sbs_bias_sig0 -->.
- **PRINTED DEFECT, reported not repaired:** eq. (24) prints
  $\ln\tau' = k_dt + \ln\!\big(\tfrac{1}{k_d}\ln\tfrac{C_{\rm A0}}{C_{\rm A}}\big)$
  where the paper's own eq. (21), one line above, requires $1/k$. The residual
  is **exactly $\ln(k_d/k)$**, and $1/k_d$ is dimensionally impossible inside
  $\ln\tau'$. It shifts the intercept of the Fig.-5 test plot, i.e. the
  recovered $k$, and leaves the slope and hence $k_d$ untouched.
- **The pymrm bed, and the paper's pseudo-steady assumption priced.** The
  finite-volume bed (`ConvFluxTVD` with the van Leer limiter, one `newton`
  solve per fluid profile, the outlet read on the extrapolated outlet face)
  meets the quadrature at second order in space on two frozen profiles and two
  reaction orders, and in the full march each axis is second order with the
  other axis's error cancelled (Heun orders from
  1.991 <!-- claim: bed_time_order_heun_min --> to
  1.998 <!-- claim: bed_time_order_heun_max -->; Euler from
  1.0007 <!-- claim: bed_time_order_euler_min --> to
  1.0035 <!-- claim: bed_time_order_euler_max -->). The
  same discrete equations solved by van Leer deferred correction agree to
  below 1e-11 <!-- claim: bed_route_gap_march -->. At production settings the
  bed lands 1.53e-6 <!-- claim: bed_vs_quadrature_production --> from the
  quadrature in conversion (largest over the observation times after the fresh
  bed), far below the resolution the verdicts are read at,
  which is what licenses the sweep to run on the quadrature. Retaining the fluid
  accumulation term the paper drops on p. 269, the cost is first order in
  $\varepsilon$ (slope 1.000 <!-- claim: pseudo_steady_cost_slope -->) and
  crosses a 0.2 % conversion resolution at $\varepsilon^*$ =
  4.30e-3 <!-- claim: pseudo_steady_eps_star -->, a run of about 4 minutes for
  a 1 s fluid residence time, which is where the paper draws the line in words
  ("Deactivation in the order of minutes or longer can use the fixed batch of
  solids", p. 267).
- **The printed structure proved.** 26 symbolic identities re-derive
  eqs. (12)–(32) and Fig. 6's slope and intercept from the predecessors the
  paper derives them from, with zero nonzero residuals; each device
  additionally meets the paper's own closed forms: eq. (15) to below 1e-7
  <!-- claim: eq15_maxdev --> (integrator tolerance), eq. (19) and
  eqs. (28)–(31) to round-off, and eq. (22) to 3.0e-6
  <!-- claim: eq22_maxdev_pymrm --> *on the pymrm bed* (its discretisation
  error).

## Fit vs test

Nothing on this page is fitted to any measurement and nothing is validated
against one: the source reports none. Every number is exact reproduction of
printed algebra, an exact algebraic identity between the printed rate forms, or
a computed property of those forms in a reactor. **Direction matters and is
stated on the page:** a separability number is the *minimum* misfit of a wrong
class, and an optimiser returns an upper bound on a minimum, so the
"not separable" results are proofs and the "separable" results mean "no better
fit was found". Where the two disagree the identity wins, in the matrix
itself: every entry is min(fit, identity twin), 44 of the 64 have a twin, and
the one entry where the fit alone stops short is printed side by side with what
the algebra gives.

## Files

- `index.py`: the page source (jupytext percent format), paired with
  `index.ipynb`.
- `model.py`: the four contactings (`run_batch`, `run_mixed`, `run_const_c`,
  the quadrature bed `PlugExact`) and the pymrm bed `Bed`; `test_model.py`
  checks them in seconds.
- `data/levenspiel-1972-printed-equations.csv`: the printed equations used,
  with the eq.-(24) defect carried as printed.
- `data/levenspiel-1972-printed-claims.csv`: the printed prose claims tested,
  quoted verbatim.
- `agreement.json`: the page's claims in format 2 (`docs/claims-design.md`),
  each registered once as a value, a bound or an exact result with how it is
  checked, and the break-table rows each claim must fail.

Runtime: about 11 minutes per run.

## Corrections

- 2026-08-09: stated that the bed is "the only device that is not exactly
  degenerate"; it is exactly degenerate too (side-by-side is parallel at
  $n = n' = d = 1$).
- 2026-10-29: stated the separability diagonal as an optimiser stall of 8.55e-6
  (batch) and 2.38e-5 (mixed flow), with nine entries where the multistart
  stopped more than 1e-6 short; those were artefacts of a one-sided
  finite-difference Jacobian, and the side-by-side diagonal is below 1e-9 on
  every device, with one such entry.
- 2026-10-29: stated three fit-only off-diagonal entries as 1.85e-5, 5.31e-5 and
  1.86e-5; they are 7.4e-6, 3.45e-5 and 1.39e-5.
- 2026-10-29: stated a $d = 1$ bias of 5.5e-7 in the Fig.-6 transform; it is
  below 1e-7.
- 2026-11-02: said the bed's Heun march reads an apparent time order of 3.16 at
  200 steps because its time error had fallen onto the 400-cell spatial floor.
  The cause was the read-off: at 25, 50 and 100 steps the 41 observation times
  are not all step times (at 100 steps 21 of the 41 are) and the conversion was interpolated between steps. With
  the space error cancelled the same ladder (400 cells, 25 to 200 steps) gives 3.28; on step counts
  whose grid holds every observation time the orders are 1.99 to 2.00.
- 2026-11-02: gave the slope of the pseudo-steady cost against $\varepsilon$ as
  1.009; the excess over 1 was the time error of the $\varepsilon = 0$
  reference (200 Heun steps). Against a reference of 1600 steps the slope is
  1.000; $\varepsilon^*$ is unchanged at 4.30e-3.
- 2026-11-02: said degeneracy 5 fails in the two well-mixed devices by 3e-2
  and 1e-1; it fails by 1.2e-1 (mixed flow) and 3.5e-1 (held $C_{\rm A}$).
