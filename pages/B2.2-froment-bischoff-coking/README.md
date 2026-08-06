# B2.2 — Froment–Bischoff coke deactivation

**Catalog ID:** `B2.2` · **Section:** B · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S4`, `S5` · **Data tier:** 6 (the paper contains no
measurements; everything is reproduction — see below)

Reproduces G. F. Froment & K. B. Bischoff, "Non-steady state behaviour of
fixed bed catalytic reactors due to catalyst fouling", *Chem. Eng. Sci.*
**16**(3–4) 189–201 (1961), doi:10.1016/0009-2509(61)80030-4 — the paper that
tied catalyst activity to the *deposited carbon* rather than to time on
stream, and showed that the coke lays down as a profile whose direction names
the fouling mechanism: descending for coke from the reactant (parallel),
ascending for coke from the product (consecutive).

Built from the 300 ppi CCITT-G4 scan of Part I; every numeral read from
cropped native-resolution renders (the text layer gives 800 °C where the page
prints 600 °C, and "0.88 to 058" where it prints 0.38 to 0.53). **Part II
(CES 17 (1962) 105) was not read** — no page image exists on disk, only
api-text of a class known to drop decimal points — and contributes nothing.

## What the page shows

- **All three printed solutions re-derived and met by one pymrm bed.** The
  closed forms (32)–(33) and (36)–(37) are verified symbolically and against
  the paper's own eq.-(29) ODE route (3.0e-13); the consecutive case gets an
  exact scalar-ODE reduction the paper does not contain. The pymrm
  finite-volume bed (upwind + van Leer deferred correction, Heun in coking
  time) meets all three at observed order ~2 in both axes, max error 1.7e-5.
- **The paper's four approximations, each measured**: γ = 0 costs 0.29 % at
  the paper's own γ = 1e-3 and 5.9 % at its printed 0.02 bound; the four-term
  series is 22 % off at its claimed corner; the large-η patch *under-counts
  carbon by a growing amount* (1.16 in αc by αbη = 10) — so the printed
  Fig. 3(b) is biased and the claim that the approximation "improves as η
  increases" is false in c; η ≈ t costs 2.2e-5 for the paper's own reactor.
- **A one-character typo in eq. (49) adjudicated** (αz printed where the
  integration of eq. (48) requires az; the readings differ by α/a ≈ 200).
- **The Fig.-4 slope claims re-examined.** The 1.0 limits reproduce; the two
  printed "slope reaches 0.5 at ~4 (CE) / ~10 (PH)" onsets do not survive
  exact computation (0.710 and 0.726 there; CE first reaches 0.5 near 13 and
  has no plateau; 0.5 is PH's true asymptote but far beyond the window). The
  Voorhies comparison sharpens: every quoted experimental exponent is
  reachable by the consecutive mechanism, the upper three by
  parallel-hyperbolic, and none by parallel-exponential, whose exponent never
  drops below 1.
- **The profile diagnostic gets an amplitude.** At the paper's own operating
  point the ascending (consecutive) signature survives ~22 % carbon assays
  with only logarithmic decay of its half-bed contrast (below a 20 % assay
floor by t' ~ 20; the inlet anchor is exact at every t); the descending
(parallel) signature is already below a
  10 %-assay detection floor and fades as 1/t. Equal-weight mixed mechanisms
  cancel the profile *identically* while producing a perfectly smooth
  Voorhies curve — the paper's "balance each other" remark made exact.
- **Two exact objects the paper did not print**: the travelling hot-spot
  locus x* = ln(e^T − 1) (unit speed), and y = W/s for the hyperbolic case.

## Fit vs test

Nothing is fitted. The paper contains no measurements, so every agreement is
a reproduction of its mathematics or an internal-consistency check between
its printed statements — never a validation. The authors' own
order-of-magnitude constants (k₂°, α, K — their fit to their unpublished
data) fix only the ranges computed, exactly as in the paper. The five quoted
Voorhies exponents are other papers' measurements, transcribed from this
paper alone and used only in the comparison the paper itself makes.

## Files

- `build_page.py` — generates `index.ipynb` (run `python build_page.py`, then
  execute the notebook).
- `data/froment-bischoff-1961-*.csv` + `.meta.yaml` — the printed Appendix-2
  conditions and the paper's quantitative prose claims, with provenance.
- `agreement.json` — 55 metrics for CI regression checking.

Runtime: about 90 seconds end to end.
