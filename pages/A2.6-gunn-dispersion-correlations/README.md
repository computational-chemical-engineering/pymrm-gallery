# A2.6 — Gunn's dispersion correlations, and what the number 2 depends on

`A2.1` takes the axial-dispersion reactor as given. `A2.3` derives a dispersion
coefficient for an empty tube. `A2.5` measures one in a packed bed and fits two
constants to the measurement. This page is the fourth question: a stochastic model
with **one** empirical function in it, covering axial and radial dispersion,
spheres and cylinders, gases and liquids — and what its constants say about the
Péclet number 2 that everyone quotes.

- **Structures:** `S6` (2-D axisymmetric dispersion–reaction)
- **Reference:** Gunn, D. J. (1987), Chem. Eng. Sci. 42(2) 363–373,
  doi:10.1016/0009-2509(87)85066-2
- **Runtime:** ~9 s

## The headline

**Pe → 2 is a property of spheres, not of packed beds.** The high-Reynolds-number
limit of Gunn's eq. (1) is `Ud/D → 2p∞/(1−p∞)`, and eq. (3) gives
`p∞ = 0.17 + 0.33 = 0.50` exactly — so the limit is exactly 2. Eqs (4) and (5),
fitted to cylinder measurements that never see a sphere, give `p∞` = 0.46 and
0.37, hence **1.7037** and **1.1746**. The habit of writing `D = u dp / 2` at high
flow is a sphere habit; for rings it is 41 % wrong.

The page also states the width of the claim: two printed digits fix the limit only
to about ±4 %. What is not within rounding is the **1.70-fold spread between
shapes** — and Gunn's own Figure 2 supports it independently of the constants:
England and Gunn's measured cylinder Péclet groups sit between about 1.0 and 1.7
above Re = 100, visibly below the 2 the same paper quotes for spheres. Eqs (4) and
(5) were fitted over Re ≈ 1–400 and exp(−24/400) = 0.942, so `p∞` is interpolated,
not extrapolated.

## Agreement

| Check | Result |
|---|---|
| high-Re limit, by expansion vs direct evaluation of eq. (1) | **7.2e-08** on all three shapes; drops to a factor **1.4e7** if eq. (1)'s second term is dropped |
| eq. (8), a worked example with every input printed | 1.65 × 0.37 / 0.6 = **1.0175** vs the printed 1.02 (−0.25 %) |
| seven numbers the paper states in words | deviations **0.0 % to 11.4 %**, each inside the precision of the phrase |
| eq. (1) [1987] vs eq. (42) [1993] with Γ = 4(1−ε)α₁²/ε | **8.7e-11** — an identity, and labelled one |
| eq. (1)+(3) vs the 55 markers `A2.5` digitised | **20.3 %** mean absolute deviation, **+15.8 %** bias |
| E&R's own eq. (18) on the same rows, recomputed here | **11.75 %**, +0.003 % — reproduces `A2.5`'s published headline to 1e-05 |
| null: `p` a **free** constant (the strongest null) | 22.96 % at `p` = 0.0793 — only **2.7 points** worse |
| null: `p` frozen at 0.17 | 24.9 % |
| null: `p` frozen at 0.50 (the mixing-cell picture) | 67.6 % |
| the unprinted porosity, ε = 0.34 vs 0.40 | 26.6 % vs 17.5 % — **9.1 points**, 3.4× the cost of deleting p(Re) |
| the porosity these markers actually prefer | **16.4 % at ε = 0.443** — an interior minimum, not a monotone trend |
| 2-D bed vs Gunn's eq. (17), flat inlet | **7.5e-07** at n_z = 800, order 1.99 |
| 2-D bed vs a Bessel-mode closed form derived here | **2.9e-05** at (800, 80), order 2.00 — ×2646 on `nu=0`, ×252 on n_r = 5 |
| wall channel: 1.65 read as interstitial | **+5.46 %** unconverted with no radial mixing, **+1.16 %** with eq. (9)'s |
| zero-gradient outlet instead of eq. (18) | **+1.03 %**, and identical from L/d = 5 to 200 |
| **the same two, recomputed by `solve_bvp` collocation** | wall channel **1.4e-04** relative; outlet **1.0353 %** vs the bed's 1.0339 % |

Deviations are `(model − reference)/reference` everywhere on the page.

## Are these markers inside Gunn's fit set? Probably not

An earlier draft of this page asserted that Gunn fitted `p(Re)` to a body of data
including Edwards and Richardson's points, and called the 20.3 % "a goodness of fit
with partial overlap". The sources do not support that.

- Gunn 1987 p. 365: `p` was estimated from experiments "in which the effect of
  molecular diffusion was very small".
- Gunn 1993 p. 335 says what that meant — `ReSc ≫ 1`, i.e. **liquid-phase**
  measurements, "for beds of spheres, solid cylinders and hollow cylinders [2, 5]"
  (Gunn 1969 and Gunn 1987) — with gas-phase results used "**to compare**".
- These 55 markers are gas phase at Sc = 0.77: `ReSc` runs **0.0066 to 37.8**,
  median 1.40, with 26 below 1 and **47 of 55 below 10**. They fail that criterion.
- The 1987 paper attributes eqs (4) and (5) to England & Gunn (1970) and eq. (3) to
  nothing, and says of *these* points that they "are well supported by eq. (1)" —
  the language of comparison. Gunn 1993's Figure 1 lists them under comparison too.

So the comparison is **probably out of sample, though not provably held out**, and
the page says exactly that rather than either extreme. It is still not called a
validation.

Two limits on what those markers can resolve, both measured rather than assumed:

- **Deleting `p(Re)` costs 2.7 percentage points against the strongest null** — a
  *free* constant `p` = 0.0793 at 22.96 %, which beats freezing `p` at a printed
  plateau. The data span 61 % of `p`'s excursion, so they do exercise it, but not
  enough to identify *which* Reynolds-dependent function it is. At each model's own
  preferred porosity the margin is 4.6 points, so it is itself porosity-dependent.
- **The porosity, which the paper never prints, moves the headline by 9.1 points**
  over ε = 0.34 to 0.40 — 3.4× as much. It is **not** monotone: the residual has an
  interior minimum at ε = 0.443 (16.4 %) and rises again above it.

## What the page adds

**A number on Gunn's wall argument.** Eqs (7)–(8) conclude the interstitial
velocity in the wall region is within 2 % of the bulk, "suggesting that there may
be no significant differences". A 2-D bed says how much: reading Price's
superficial 1.65 as interstitial raises unconverted reactant by 5.46 % with no
radial mixing and 1.16 % with eq. (9)'s coefficient, so radial dispersion erodes
79 % of the channel; Gunn's corrected 1.0175 moves it by 0.001 %.

**A closed form for the 2-D problem.** Feeding the bed a zero-flux radial
eigenmode makes eq. (17) exact again with `k → k + D_R(β₁/R)²`. It is the only
check on the page that reads the radial operator — the flat-inlet check against
eq. (17) leaves `nu=0` at *exactly* its baseline value, so a page carrying only
that one could put a Cartesian divergence in the radial direction and report
machine precision.

**A correction to the obvious reading of "except where the reactor is long".** The
zero-gradient outlet's error at the *outlet* is a length-independent multiplicative
bias — 1.03 % at L/d = 5 and 1.03 % at L/d = 200. What length buys is the
*interior*: the same error 4 diameters in falls from 0.15 % to 9.4e-08 to 4.4e-16
at L/d = 20 — one unit in the last place — and to a true zero beyond. Over a
400-fold sweep in `kL/U` the outlet error is 0.911 to 0.998 times `λD/U`, so that
group is what to compute, not particle diameters.

**A second, independent computation of two of those numbers.** V8 recomputes the
outlet-condition error and the wall-channel cost with `scipy.integrate.solve_bvp` —
collocation on a different mesh, in a different solver, with the reference
asserting `status == 0` and its own rms residuals rather than being trusted. The
wall-channel cost agrees to **1.4e-04**, which is the only check on the page that
reads the velocity profile, the area weights and the mixing-cup weighting (V5a and
V5b both run a uniform profile). **The outlet error did not agree**, and the
disagreement was this page's error: it had been read at the last cell *centre*,
h/2 short of the outlet, and since eq. (18) and the zero-gradient outlet differ in
outlet gradient by construction, the offset did not cancel in their ratio. That
made the metric 11.4 % low, and no break row could have caught it — a break row
asks whether a number *moves*, not whether it is *right*. Read at the face it is
1.0339 % against the collocation 1.0353 %.

## Honest limits

- **The porosity is an assumption**, taken from a sentence on journal page 368 that
  is about tube diameter, not about any figure. It is the largest lever on every
  comparison here.
- **No liquid-phase data.** The Schmidt-group dependence over 1000-fold, which is
  what Gunn's Figure 1 exists to show and the single most valuable extension to
  this page, is untested. Jacques & Vermeulen (1957) or Miller & King (1966) would
  supply it.
- **The borrowed dataset carries `A2.5`'s open questions.** Three glyphs remain
  fused (including the most nearly resolvable, whose inclusion moves this page's
  headline 20.30 % → 20.65 %); a ±5 px ink-centroid systematic moves it to 19.27 %
  or 21.59 %; `A2.5`'s figure review is pending and non-blocking. **No new
  digitisation was done and none is proposed** — there is no maintainer available.
- **Two metrics have no break row and are labelled structural**: the cross-document
  identity against Gunn 1993's eq. (42), and eq. (17) substituted back into
  eq. (16). Several more are guarded by a family rather than an individual row, and
  four by a second independent computation.
- **The two closed-form solver checks cannot see a coefficient that is wrong
  everywhere**, measured: a 10 % error made *consistently* in `D` or `D_R` — in the
  operator, in `lam`, in `k_eff` and in the closed form together — moves them by
  0.0 % and 6.3 %. The same 10 % injected into the **operator alone** moves V5a
  ×1369 and V5b ×53.7 for `D`, and V5b ×445 for `D_R`. An operator/reference
  mismatch is caught loudly; only a shared error gets through.
- **One stated range does not come out**: the upper end of "0.2 < Re < 300" is not
  a crossing — the liquid/gas ratio decays asymptotically and is still 1.084 at
  Re = 300. Stated, not repaired.
- **Three printing inconsistencies in the source are recorded, not repaired**: the
  Figure 5 caption names "eq. (4)" where eqs (9) and (10) are required; the
  tortuosity for spheres is 1.4 in eq. (3) and 1.2 in eq. (10) — 17 %, described
  by the paper only as "small amounts"; and eq. (2) prints `16α₁⁴(1−ε)` where
  eq. (1), which it must reduce to, prints `16α₁⁴(1−ε)²`. Gunn 1993's eq. (42)
  confirms the square.
- **The reactor operating point is chosen, not printed.** The paper has no worked
  reactor example; every number of the operating point is on the page. Its
  Danckwerts inlet is written with the area-mean velocity rather than the local
  `u(r)`; the alternative is a break row worth 1.1 % of the wall-channel answer.

## Rebuilding

```bash
python build_page.py                    # regenerate index.ipynb
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF needed. The two Gunn datasets are committed CSVs; the 55 measured points
are loaded from `pages/A2.5-edwards-richardson-dispersion/data/` through
`load_data(..., page=...)`.
