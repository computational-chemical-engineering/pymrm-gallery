# A2.5 — Where the axial dispersion coefficient comes from

`A2.1` takes the axial-dispersion reactor as given. `A2.3` derives a dispersion
coefficient for a tube. This page is the experimental one: a pulse of argon in
air, two response curves, and the two-point moment method that turns them into
`D_L`.

- **Structures:** `S4` (1-D transient advection–dispersion)
- **Reference:** Edwards & Richardson (1968), Chem. Eng. Sci. 23(2) 109–123,
  doi:10.1016/0009-2509(68)87056-3
- **Runtime:** ~321 s

## Agreement

**Experimental**, against 55 marker positions digitised from Figure 9.

| Check | Result |
|---|---|
| Eq. (18) vs 55 markers | **11.75 %** mean absolute deviation, bias +0.003 %, over 3.7 decades of Re |
| null: Eq. (15), no radial-mixing term | 29.9 % with a −29.8 % bias; its ceiling of Pe = 2 sits below **30 of the 55** points |
| refit β to the digitised points (γ fixed) | **8.98** against the printed 9.7; block bootstrap 95 % interval [7.50, 9.94] |
| refit both parameters — a *different* estimator | γ = 0.702, β = 8.58; intervals [0.672, 0.745] and [6.90, 9.50] |
| Eq. (18)'s constants from Eq. (17)'s | 1.28 % and 0.31 %; one mis-read digit gives 14.8 % and 10.6 % |
| gas crossover Re, printed as 1.8 | **1.830** from printed values only — +1.7 %, not exact |
| Table 1 voidages vs 1 − bulk/material density | worst 0.0018 absolute over all eight rows |
| two-point moment inversion of a simulated experiment | **1.4e-05** relative at n = 3200, observed order 3.01 |
| numerical dispersion predicted before the runs | `u h/2 + u²Δt/2` accounts for the excess to **0.4 %** over 8 runs |
| two-zone model vs an independent closed form | 3.1e-05 at n = 2000 — and broken on purpose, ×173 and ×1058 |

Deviations are `(model − measured)/measured` everywhere on the page.

## The 11.75 % is a goodness of fit, not a prediction

γ and β were fitted by Edwards and Richardson to the very markers this page
digitises. What *is* independent:

- the **shape**. Eq. (15) has the right molecular branch and the right high-Re
  limit and still cannot exceed Pe = 2, where over half the data do. The
  maximum in the Péclet–Reynolds plot is the evidence for the third mechanism.
- the **refit**. Recovering β = 8.98 from a fresh digitisation tests the
  extraction and the page-image transcription of Eq. (17) against each other.
  No difference from the printed 9.7 is claimed: each estimator's bootstrap
  interval is quoted beside its own estimate, and whether 9.7 falls inside the
  γ-fixed one is decided in the second decimal by the block size.
- the **location of the maximum**: Pe = 3.186 at Re = 3.195, against the paper's
  own bound ("greater than 3"). Gunn (1993) puts a maximum at Re ≈ 4, but for
  his Figure 1's *combined* gas-phase set — Gunn & Pryce's data together with
  E&R's — not for E&R's alone.

## What the page adds

The measurement is an inverse problem, so the page simulates it. A pulse
propagated on a pymrm bed, sampled at two detectors 100 cm apart and inverted
with Eqs. (6) and (7), returns the input `D_L` to 1.4e-05 — and the numerical
dispersion of the cheaper schemes is **predicted before the runs**: `u h/2` for
upwind plus `u²Δt/2` for implicit Euler, good to 0.4 %. At n = 400 the bare
upwind scheme reports a dispersion coefficient 4.5 times the truth from response
curves that look perfectly smooth. In an inverse problem, truncation error does
not appear as a wiggle; it appears as a plausible parameter value.

The last section turns the authors' own shape diagnostic round. They report that
for fine particles the computed second response curve peaks *later* than the
measured one and attribute it to channels forming. A two-zone bed constrained to
their measured `D_L` = 0.503 cm²/s reproduces that sign when the zone holding
the smaller share of the void is the **slower** one — in 10 of 12 volume splits,
and in all 8 where the peak moves by more than 0.05 % of the transit time at
all. The two exceptions sit near an even split, where the displacement is 0.011 %
and effectively signless. The claim is therefore conditional: *when* such a bed
displaces the peak measurably, the direction E&R report is the one a slow
minority zone gives.

The *size* of that displacement bounds nothing. At the published θ₁ = 0.15
slice it is at most 0.36 % of the transit time, but 1.14 % with an 18.7 % shape
residual once the volume fraction is allowed to move at the same measured `D_L`.
What the page concludes is the either/or: a discrepancy large enough to see
implies either a slow region rather than a channel, or structure coarser than
any two-zone description tuned to 0.503 cm²/s can hold.

## Honest limits

- No particle-size test. Figure 9's five marker shapes were deliberately not
  extracted, because Eq. (18) contains no particle-size term, so the paper's "no
  consistent trend with particle size" is untested here.
- One gas pair only, so the ν in Re is untested. The authors say so themselves.
- The extraction was audited and corrected: four merged clusters re-fitted glyph
  by glyph, one row moved off a bare stretch of the dashed curve, one missed
  marker added, three glyphs still fused and listed in the sidecar with the
  measured effect of including the most nearly resolvable of them (11.75 % →
  11.62 %). The figure review is **pending and non-blocking**; the overlay is
  `queue_cases/A2.5/review/fig9-overlay.png` (git-ignored).
- The break table's two solver diagnostics are complementary but **not jointly
  sufficient**: reading the voidage off the adjacent row of Table 1 moves `u` by
  6.1 % and `D_L` by 6.9 % and moves no diagnostic on the page at all. They test
  the solver, not the parameter chain.
- One printed check does not come out. The liquid crossover Reynolds number
  3 × 10⁻⁴ needs ν ≈ 0.018 cm²/s where water at 20 °C is 0.0100. Stated, not
  repaired. The gas one comes out to 1.7 %, which is not "exactly".
- Figure 15's detector separation is **not printed**. The reconstruction assumes
  Figure 13's 100 cm; the page reports what the paper's other test section,
  21.3 cm, gives instead.
- The outlet boundary condition is a **measured** blind spot: with 80 cm of bed
  past the second detector, Dirichlet and zero-gradient outlets return identical
  numbers. That question belongs to `A2.1`.

## Rebuilding

```bash
python build_page.py                    # regenerate index.ipynb
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF needed: the three datasets are committed CSVs. `review/extract-fig9.py`
regenerates the Figure 9 dataset from a 600 dpi render of the page, which the
repository deliberately does not carry.
