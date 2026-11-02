# B2.2 Froment-Bischoff coke deactivation

**Catalog ID:** `B2.2` · **Section:** B · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S4`, `S5` · **Data tier:** 6 (the paper contains no
measurements; everything is reproduction, see below)

Reproduces G. F. Froment & K. B. Bischoff, "Non-steady state behaviour of
fixed bed catalytic reactors due to catalyst fouling", *Chem. Eng. Sci.*
**16**(3-4) 189-201 (1961), doi:10.1016/0009-2509(61)80030-4, the paper that
tied catalyst activity to the *deposited carbon* rather than to time on
stream, and showed that the coke lays down as a profile whose direction names
the fouling mechanism: descending for coke from the reactant (parallel),
ascending for coke from the product (consecutive).

Built from the 300 ppi CCITT-G4 scan of Part I; every numeral read from
cropped native-resolution renders (the text layer gives 800 °C where the page
prints 600 °C, and "0.88 to 058" where it prints 0.38 to 0.53). **Part II
(CES 17 (1962) 105) was not read**: no page image exists on disk, only
api-text of a class known to drop decimal points, and it contributes nothing.

## What the page shows

- **All three printed solutions re-derived and met by one pymrm bed.** The
  closed forms (32)-(33) and (36)-(37) are verified symbolically and against
  the paper's own eq.-(29) ODE route; the consecutive case gets an exact
  scalar-ODE reduction the paper does not contain. The pymrm finite-volume bed
  (van Leer TVD by `ConvFluxTVD`, each gas profile by `newton`, Heun in coking
  time) meets all three, worst error in y 1.7e-5 <!-- claim: ce_maxerr_y_n400 -->
  on 400 cells; PE and PH converge at second order, CE at orders rising, 1.73 to 1.93.
- **The paper's four approximations, each measured**: γ = 0 costs 0.29 % at
  the paper's own γ = 1e-3 and 5.9 % at its printed 0.02 bound; the four-term
  series is 22 % off at its claimed corner; the large-η patch *under-counts
  carbon by a growing amount* (1.16 in αc by αbη = 10), so the printed
  Fig. 3(b) is biased and the claim that the approximation "improves as η
  increases" is false in c; η ≈ t costs 2.2e-5 for the paper's own reactor.
- **A one-character typo in eq. (49) adjudicated** (αz printed where the
  integration of eq. (48) requires az; the readings differ by α/a, of order
  10^2).
- **The Fig.-4 slope claims re-examined** with the exponent computed exactly.
  The 1.0 limits reproduce; the two printed "slope reaches 0.5 at ~4 (CE) /
  ~10 (PH)" onsets do not (0.710 and 0.726 there; CE crosses 0.5 once, near
  13, with no plateau; 0.5 is PH's true asymptote but far beyond the window).
  The Voorhies comparison sharpens: every quoted experimental exponent is
  reachable by the consecutive mechanism, the upper three by
  parallel-hyperbolic, and none by parallel-exponential, whose exponent never
  drops below 1.
- **The profile diagnostic gets an amplitude.** At the paper's own operating
  point the ascending (consecutive) signature needs carbon assays of about
  21 % and fades only logarithmically (below a 20 % assay floor from
  αbη = 7.8 <!-- claim: diag_ce_floor20_crossing_T -->; the inlet anchor
  c(0) = 0 is exact at every time); the descending (parallel) signature is
  already below a 10 %-assay detection floor and fades as 1/t. Equal-weight
  mixed mechanisms cancel the profile *identically* while producing a smooth
  Voorhies curve, the paper's "balance each other" remark made exact.
- **Two exact objects the paper did not print**: the travelling hot-spot
  locus x* = ln(e^T - 1) (unit speed), and y = W/s for the hyperbolic case.

## Fit vs test

Nothing is fitted. The paper contains no measurements, so every agreement is
a reproduction of its mathematics or an internal-consistency check between
its printed statements, never a validation. The authors' own
order-of-magnitude constants (k₂°, α, K, their fit to their unpublished
data) fix only the ranges computed, exactly as in the paper. The five quoted
Voorhies exponents are other papers' measurements, transcribed from this
paper alone and used only in the comparison the paper itself makes.

## Corrections

- 2026-11-02: said the distance-averaged coke at αbη = 5 is 4.553 and 4.063
  (PE, aL = 1 and 2), 1.313 and 1.883 (CE). These were read at the point of a
  240-point logarithmic time grid nearest 5, which is 5.049. At 5 they are
  4.505 <!-- claim: spacetime_pe_cbar_aL1_T5 -->, 4.015, 1.304 and
  1.871 <!-- claim: spacetime_ce_cbar_aL2_T5 -->; the signs of the space-time
  test are unchanged.
- 2026-11-02: said the ascending (consecutive) signature "never fades". Its
  half-bed contrast falls from 0.606 at αbη = 5 to
  0.214 <!-- claim: diag_halfcontrast_ce_T1e3 --> at 10^3; only the inlet
  anchor c(0) = 0 is permanent. The same section said the signature survives assays of about
  22 % and drops below a 20 % floor "by t' ~ 20": the old text placed the crossing
  of the 20 % floor near t' = 20, but it is crossed at αbη = 7.8. The contrasts
  were averages over 400 nodes per half bed; that changes the assay needed at
  αbη = 5 from 21.5 % (shown as 22 %) to 21.4 %
  <!-- claim: diag_sigma_max_ce_base -->.
- 2026-11-02: said refining n alone at 720 coking steps gives order 1.93. That
  was an average over a ladder of 2.02, 2.08 and 1.69, whose last step met the
  time error of the 720 steps. With the time error cancelled (differences of outlet
  values between successive grids at 240 steps) the orders are 2.002 to 2.008.
- 2026-11-02: attributed the consecutive case's lower observed order (1.81) to
  the van Leer limiter. Unlimited second-order upwind gives lower orders
  still, and refining n alone gives the same errors; the largest error sits in
  the first cells, where the consecutive coke rises from zero, and the orders
  rise with refinement (1.93 from 400 to 800 cells).
- 2026-11-02: said the γ = 0 closure deviates by 5.87 % at the printed bound
  γ = 0.02, the largest over ten sampled times; the largest over the window,
  at αbη = 2.16, is 5.88 % <!-- claim: gamma_claim_dev_at_002 -->. And said
  the pymrm bed at γ = 0.02 agrees with the eq.-29 route to 1.1e-6, a
  comparison that interpolated the reference linearly onto the cells; at the
  cell centres the agreement is 1.3e-6 <!-- claim: pymrm_vs_eq29_gamma002 -->,
  the bed's discretisation error.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        the exact references, the paper's large-eta patch, the exact
                Voorhies exponent and profile contrast, and the pymrm bed
test_model.py   fast checks of model.py and two stored claims
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
data/           the printed Appendix-2 conditions and the paper's quantitative
                prose claims, each with a .meta.yaml provenance sidecar
```

Runtime: about 3 minutes end to end.
