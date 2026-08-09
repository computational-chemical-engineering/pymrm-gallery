# B2.3 — Levenspiel's four deactivation classes

**Catalog ID:** `B2.3` · **Section:** B · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1`, `S5` · **Data tier:** 6 (model-defining equations
transcribed from the source; **the source publishes no data**)

Reproduces Octave Levenspiel, "Experimental Search for a Simple Rate Equation
to Describe Deactivating Porous Catalyst Particles", *Journal of Catalysis*
**25**(2) 265–272 (1972), doi:10.1016/0021-9517(72)90227-8 — the paper that
sets out the four broad classes of catalyst deactivation (parallel, series,
side-by-side, independent) as $n$th-order kinetics with $n'$-order, $d$-order
decay, and asks *"whether experiments can be devised so as to give simply the
orders and rate constants of these equations"*.

**This is Levenspiel the author, not Levenspiel the book.** *Chemical Reaction
Engineering* is a different work — this paper cites its 2nd edition, Chap. 15,
as its own ref. 3 — and remains outstanding as a textbook-canonical source; it
is what `A2.4` is catalogued to. Nothing here closes that.

Built from the 300 ppi CCITT-G4 scan. The OCR text layer is badly damaged
("detrrmination", "ordrrs", "particlrs", "parallrl", "rractor") and was used
for nothing: every symbol, subscript, prime and exponent was read from a
cropped native-resolution render at glyph scale — which matters because the
paper turns on telling $n$ from $n'$, $k$ from $k'$ from $k_d$ from $k'_d$, and
$C_{\rm A}$ from $C_{\rm R}$ from $C_{\rm P}$.

**There is no data in this source, and none is used.** All eight pages were
rendered at native resolution and read: there is no table anywhere, and
Figs. 2–6 are schematic test plots with unlabelled axes and no tick values
(Fig. 6's own caption: "A straight line on this plot indicates that the guessed
order of deactivation is correct"). No measured order is compared against,
because none is reported.

## What the page shows

- **All four batch-solids devices are exactly degenerate somewhere — five
  identity families, and each collapse is algebra, not a numerical near-miss.**
  The sharpest is the constant-flow mixed reactor: its performance
  relation makes $C_{\rm R} = \kappa\,C_{\rm A}^{\,n}\,a$ *exactly* at every
  instant, so **series deactivation of order $(n', d)$ is identical to parallel
  deactivation of order $(n n',\, d+n')$** — same observable, every time, any
  $n$; the printed 1.3e-9 across nine $(n, n', d)$ combinations is the
  *integrator tolerance*, not a separation (rtol 1e-13 gives 1.1e-13). The
  held-$C_{\rm A}$ device makes parallel, series and independent literally the
  same equation (8.9e-16), and at $d = 1$ side-by-side joins them. Wherever the
  poison profile equals the reactant profile — $d = 1$, $n' = n$,
  $\sigma = \kappa/\lambda$ — side-by-side **is** parallel with the same
  constants, in the batch, in the mixed reactor *and in the plug-flow bed*
  (0.0, 7.3e-15, 0.0 on the condition; 1.66e-2 off it). And with $\sigma$ left
  free, at $n = n' = d = 1$ the poison profile is a *power* of the reactant
  profile, so **side-by-side with any $\sigma$ is parallel with
  $n' \to \lambda\sigma/\kappa$** — 2.5e-9 in the batch and 1.9e-11 in the bed
  over $\sigma$ = 0.2–3, against 1.2e-1 and 3.5e-1 in the two well-mixed
  devices, where it correctly fails. $n' = 0$ collapses all four everywhere
  (0.0 exactly), which is the paper's own statement.
  **The bed is therefore *not* "the only device that is not exactly
  degenerate"** — an earlier version of this page said so on six surfaces and
  it was false. What the bed and the batch keep, and the two well-mixed devices
  do not, is a distinction between parallel, series and independent; what the
  bed alone has is a spatial coordinate.
- **Separability measured, not asserted — 64 fits, and the answer is
  negative.** Every class is fitted to every other class's response with all
  its parameters free ($\kappa$, $\lambda$, $d$, $n'$, and $\sigma$ for
  side-by-side), and the residual it cannot remove is the conversion resolution
  an experiment would need. At 0.2 conversion percentage points, **one of the
  twelve class pairs is separable in one of the four devices** (plug flow,
  parallel vs independent, 2.20e-3). Parallel vs series in a bed is 8.07e-4 —
  real, but a tenth of a percentage point. **44 of the 64 entries have an exact
  identity behind them and the matrix reports min(fit, identity twin), not the
  fit**: on nine of those the multistart alone stopped more than 1e-6 short of
  a minimum that is exactly zero, by up to 1.31e-3. All nine errors were in the
  safe direction and no verdict moved, but printed as fits they would have read
  as measured separations.
- **The diagonal is an optimiser stall, and the 20× floor built on it never
  bound.** A class fitted to its own response has exact minimum zero, so the
  8.55e-6 (batch) and 2.38e-5 (mixed) it stops at are solver stalls — at
  recovered parameters 3.6 % and 10.5 % from the truth, while a 1 % move along
  any single parameter costs 17–311× the stall. Twenty times the worst stall is
  below the instrument resolution on every device (worst case 0.24 of it), so
  the floor is gone and the thresholds are the instrument resolution and
  nothing else; the notebook prints the check.
- **The discriminating measurement is the profile, not the effluent — and both
  inlet values are exact.** At the feed face $c = 1$, so $C_{\rm R} = 0$: the
  series inlet **never decays at all** ($a \equiv 1$), and the parallel inlet
  is $a = e^{-\lambda\theta} = 0.049787$. The finite-volume first cell centre
  approaches each at first order (series errors 2.19e-2 / 1.11e-2 / 5.59e-3 /
  2.80e-3 at $m$ = 100/200/400/800), so the m = 200 cell-centre reading 0.9889
  is **not** reported as the inlet activity. Outlets, which have no closed
  form, are extrapolated **twice** — in space from the two finest grids and
  then in time, because the space extrapolation performed at a fixed step
  count is still time-limited by 9.8e-6 (parallel), an order more than the
  space extrapolation's own residual of 1.8e-7. The doubly-extrapolated
  outlets are 0.19017 and 0.17711, quoted to the five decimals the two
  refinements resolve, so the gradient reverses sign (+0.1404 vs −0.8229). The
  break table now carries both halves of that: halving the step count moves
  the space-only reading by 2.9e-5, and does **not** move the reported one.
  One assay of the spent catalyst reads
  that off, where the outlet history is a whole campaign at a resolution nobody
  has. This is `B2.2`'s "the coke profile names the fouling mechanism", reached
  from the opposite direction — with `B2.2`'s own caveat carried: it measured
  its *within-bed* contrast below a 10 %-assay floor at its operating point, a
  different quantity from the inlet activity used here, and the page says so.
- **The mixed-reactor degeneracy has an analytic boundary, and the measurement
  shows what it is worth.** The parallel-to-series relabelling needs
  $d_s = d - n'/n \ge 0$, and $d$ is an *order*, so below $d = n'/n$ no exact
  twin exists. Measured, with the twin evaluated directly where it exists: at
  or above the line the best impersonation is ≤ 1.5e-9 (nothing at all);
  below it the residual grows smoothly — 7.5e-4 at
  $d - n'/n = -0.25$, 4.0e-3 at −0.75 — and only crosses 2e-3 near −0.5. The
  line is a prediction about *exactness* and it holds; it is not a prediction
  about measurability, and the $(n' = 2, d = 1.5)$ row (half a unit below the
  line, still fitted to 4.16e-5) shows why.
- **Levenspiel's recommended device does exactly what he says, and the page
  prices what that costs.** Holding $C_{\rm A}$ by lowering the flow rate gives
  $d$ cleanly from the shape of $\tau'(t)$ — eqs. (28)–(32) reproduce to ~1e-16
  at $d$ = 0, 1, 2, 3 and for general $d$ — but *because* it decouples, a single
  run cannot name the class at all. Recovering the class needs a set of runs,
  and the span is root-found: with $k'_d$ known to 5 % the held levels must span
  0.204 in $C_{\rm A}/C_{\rm A0}$; at 2 %, 0.111. **That 0.204 is an order, not
  a design constant**: it depends on two choices the paper never prints — where
  the level set starts and how many levels it has — and root-finding at 5 %
  over $C_{\rm LO}$ = 0.05–0.40 and 3–9 levels gives 0.073 to 0.382, a factor
  of five, monotone in the starting level. The sweep is on the page. A narrow
  bracket of levels tells you nothing whichever way you choose it.
- **One knob, two constraints — with an exact blind spot at $d = 1$.** With a
  poison that is consumed, holding $C_{\rm A}$ does not hold $C_{\rm P}$. But
  the consumption rate goes as $a^{d-1}$, so at $d = 1$ it is *constant* and
  the Fig.-6 line stays perfectly straight however much poison is consumed
  (bias ≤ 5.5e-7 for every $\sigma$ up to 20). At $d = 2$ the drift is real:
  the recovered order is 1.383 instead of 2 at the base $\sigma$ = 0.8
  (−30.8 %), 5 % bias is reached at $\sigma$ = 0.076, and the worst over the
  scanned range is 49 %. **But the protocol priced here is not the one the
  paper prescribes.** Levenspiel's instruction for eq. (6) is "keep
  $C_{\rm P}$ constant", with $C_{\rm A}$ the conditional extra — "for
  side-by-side deactivation if $C_{\rm P}$ and, *if possible*, $C_{\rm A}$ are
  kept constant. **If $C_{\rm A}$ cannot be kept constant in side-by-side
  deactivation, analysis is still not particularly difficult**" (p. 271; both
  sentences are rows of the claims CSV and both are now quoted in full on the
  page). Under his own protocol — the $\sigma$ = 0 row, poison in excess so
  $C_{\rm P}$ never moves — the order comes back unbiased to 5.4e-9. What is
  sized here is the cost of holding $C_{\rm A}$ instead, which is what a single
  flow-rate knob forces; it is not a gap the author missed.
- **PRINTED DEFECT, reported not repaired:** eq. (24) prints
  $\ln\tau' = k_dt + \ln\!\big(\tfrac{1}{k_d}\ln\tfrac{C_{\rm A0}}{C_{\rm A}}\big)$
  where the paper's own eq. (21), one line above, requires $1/k$. The residual
  is **exactly $\ln(k_d/k)$**, and $1/k_d$ is dimensionally impossible inside
  $\ln\tau'$ (seconds, where $\tau'$ is g cat·sec/liter — both from the paper's
  Nomenclature). It shifts the intercept of the Fig.-5 test plot, i.e. the
  recovered $k$, and leaves the slope and hence $k_d$ untouched.
- **The pymrm bed, and the paper's pseudo-steady assumption priced.** The
  finite-volume bed meets an independent quadrature at observed order 2.01 in
  space (2.04 at $n$ = 2; 2.00 on a second frozen profile) and 1.97 in time —
  against 1.01 for Euler on the same refinement pair, which is what the
  time-order break row now compares (order against order; the earlier version
  set an order against an error and could not fail). The Heun row at nt = 200
  reads an apparent order 3.16 because its time error has fallen onto the
  m = 400 spatial floor, and the page says so.
  It lands 2.81e-6 from the quadrature in conversion at production settings — 711× below the
  resolution the verdicts are read at, which is what licenses the sweep to run
  on the quadrature. Retaining the fluid accumulation term the paper drops on
  p. 269, the cost is first order in $\varepsilon$ (slope 1.009) and crosses a
  0.2 % conversion resolution at $\varepsilon^*$ = 4.30e-3 — a run of about
  4 minutes for a 1 s fluid residence time, which is exactly where the paper
  draws the line in words ("Deactivation in the order of minutes or longer can
  use the fixed batch of solids", p. 267).
- **The printed structure proved.** 26 symbolic identities re-derive
  eqs. (12)–(32) and Fig. 6's slope and intercept from the predecessors the
  paper derives them from, with zero nonzero residuals; each device
  additionally meets the paper's own closed forms — eq. (15) to 8.6e-9,
  eq. (19) to 2.7e-15, eq. (22) to 5.46e-6 *on the pymrm bed*, and eqs. (28)–(31)
  to ~2e-16 on all four printed branches of $d$.

## Fit vs test

Nothing on this page is fitted to any measurement and nothing is validated
against one — the source reports none. Every number is exact reproduction of
printed algebra, an exact algebraic identity between the printed rate forms, or
a computed property of those forms in a reactor. **Direction matters and is
stated on the page:** a separability number is the *minimum* misfit of a wrong
class, and an optimiser returns an upper bound on a minimum — so the
"not separable" results are proofs and the "separable" results mean "no better
fit was found". Where the two disagree the identity wins, and it now wins *in
the matrix itself*: every entry is min(fit, identity twin), 44 of the 64 have a
twin, and the nine entries where the fit alone stopped short are printed side
by side with what the algebra gives.

## Files

- `build_page.py` — generates `index.ipynb` (run `python build_page.py`, then
  execute the notebook).
- `data/levenspiel-1972-printed-equations.csv` — the printed equations used,
  with the eq.-(24) defect carried as printed.
- `data/levenspiel-1972-printed-claims.csv` — the 17 printed prose claims
  tested, quoted verbatim.
- `agreement.json` — 217 metrics for CI regression checking; break-row coverage
  asserted key-for-key in the notebook (24 injected defects, every one of which
  moves its metric, plus three robustness rows that must *not* move and do
  not — the widened multistart, and each bed outlet re-extrapolated off the
  other step-count pair).
  63 metrics sit below `check_agreement.py`'s `ABS_FLOOR` = 1e-12 and are
  documented rather than protected. They are **not** all identities: the
  notebook classifies every one of them into 12 kinds — exact degeneracies,
  identity-backed matrix entries, three *counts* that carry the page's negative
  headline, a budget counter, a boolean, and five numerical ODE-vs-closed-form
  deviations — states for each what it cannot detect, and asserts that the
  above-floor companion it names is itself above the floor.

Runtime: about 5 minutes end to end (310 s measured; two consecutive executions gave 309.8 s and 307.2 s).
