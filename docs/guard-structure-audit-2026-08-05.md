# Guard-structure audit, 2026-08-05

Repository-wide sweep of all **41 built page directories** for the three structural
shapes that "the check that cannot fail" has been found to take. Each shape had
previously been found on exactly one page, by a verifier who happened to look;
none had been swept for.

> **Shape 1** — the diagnostic that differentiates the object it tests
> (`A4.7`'s `GAMMA_IDENTITY`). Two sides computed from one object, so agreement is
> algebraically guaranteed for any parameter values: it can detect a wrong
> *formula*, never a wrong *number*.
>
> **Shape 2** — every check evaluated at the same point, on a branched or
> piecewise model (`A1.8`'s four break tables, all on the dense branch). Printed
> constants on the untouched branch can be deleted outright.
>
> **Shape 3** — the refinement study that refines the wrong knob. (a) a knob that
> does not control resolution (`A4.7`'s geometric `dt`); (b) a page with both a
> steady and a transient result where only one is refined (`H1.9`); (c) a
> convergence claim with no observed order.

Method: read `build_page.py` (or `index.ipynb` where there is none — `A4.9`,
`B1.1`, `F3.1`) for every page, the code and not the prose; then inject the
defect and execute. `report_agreement` was stubbed in every run and
`git status --porcelain -- pages/` returned **0 entries** before, during and
after. Nothing in `pages/`, `models.yaml` or any other doc was modified.

## Result in one paragraph

**Fifteen pages are clean on all three shapes. Three more are the known
instances, and all three are confirmed already repaired.** Ten pages carry a
shape that reaches a published headline, a reported number or Reuse advice, and
thirteen carry one that is incidental or already disclosed in the page's own
prose. Two of the ten are the same defect propagated by copying, and two more are
the same defect propagated by a shared helper — **the dominant transmission
vector is page-to-page inheritance, not independent invention.** Shape 3(a), the
geometric-`dt` trap, produced **zero** new instances and is the one shape rare
enough not to need its own rule.

| outcome | count | pages |
|---|---|---|
| **clean on all three shapes** | **15** | `A1.5` `A1.6` `A2.1` `A2.5` `A3.4` `A4.3` `B1.4` `B3.1` `E1.2` `F1.4` `F3.1` `G1.8` `H1.1` `H1.4` `J4.8` |
| known instance, **confirmed repaired** | 3 | `A1.8` `A4.7` `H1.9` |
| carries a shape, **reaches a headline / reported number / Reuse** | 10 | `A2.3` `A4.4` `A4.9` `B1.1` `C2.1` `E2.1` `F2.3` `F3.5` `J1.5` `J3.4` |
| carries a shape, incidental or self-disclosed | 13 | `A1.1` `A1.7` `A4.2` `B1.2` `B1.6` `D2.1` `D2.2` `F1.3` `G1.7` `H1.7` `I1.2` `J3.1` `J3.5` |

### The correlate that explains almost all of it

**23 of the 41 pages carry no defect-injection table at all.** Nine of the ten
headline-reaching findings are on such a page; the tenth (`A4.4`) has a break
table that simply does not cover the constant in question. Ten of the fifteen
clean pages have one.

A second structural fact worth recording: `check_agreement.py` ignores any metric
below `ABS_FLOOR = 1e-12`. **74 of the 920 published metrics sit below that floor**
and are therefore invisible to CI entirely. Every Shape-1 identity in this report
is in that set. A page whose only guard is an identity has no CI coverage at all,
which is why these survive re-execution indefinitely.

---

## Findings, ranked by reach

### 1. `E2.1` — three of the page's checks are algebraic identities, and `meta.yaml` ranks them *above* the checks that work. CONFIRMED

`pages/E2.1-kunii-levenspiel-bubbling-bed/build_page.py:695-707`, `:684-693`,
`:499-506`. All three are in `agreement.json`; `meta.yaml:64-68` lists them under
**"The internal checks are stronger"** than the appendix comparisons.

**(i) `mass_balance_dev` = 8.106e-06, printed as `closure 0.001 %`.** The exchange
block `S` (`:382-387`) has column sums `−δγ_k K_r` exactly — the `K_bc` and `K_ce`
terms cancel pairwise between rows — so summing the discrete equations telescopes
the convective flux, and with `ve = umf − vc` (`:363`) the three velocities sum to
`u₀` identically. The equality is forced for **any** `γ_b, γ_c, γ_e, K_bc, K_ce,
K_r, δ, v_c`.

Measured by re-executing the page's own cells (baseline reproduces the published
`8.106057127982915e-06` and `6.3421228470936985e-09` to every digit):

| injected defect | `mass_balance_dev` | `eq57_limit_dev` | printed verdict |
|---|---|---|---|
| as published | 8.106e-06 | 6.342e-09 | `closure 0.001 %` |
| `K_bc` coefficient **4.5 → 3.5** (22 % error) | 9.307e-06 | 6.634e-09 | `closure 0.001 %` |
| `K_bc` coefficient 5.85 → 5.58 | 8.366e-06 | 6.405e-09 | `closure 0.001 %` |
| `K_ce` coefficient **6.78 → 8.76** (29 % error) | 1.461e-05 | 5.475e-09 | `closure 0.001 %` |
| `u_br` coefficient **0.711 → 0.171** | −2.035e-03 | 7.601e-09 | `closure 0.2 %` |

**(ii) `eq57_limit_dev` = 6.34e-09.** As `K_r → 0`, `Kf_closed` (`:333-339`)
collapses to `(1/(1−ε_mf))(u₀/u_br)K_m(γ_b+γ_c+γ_e)`; but `γ_e` is *defined*
(`:321`) as `(1−ε_mf)(1−δ)/δ − γ_c − γ_b`, so the bracket is `(1−ε_mf)(1−δ)/δ`
and, with `(1−δ)u_b = u_br`, everything cancels to `u₀/(u₀−u_mf)` identically.
`0.711`, `4.5`, `5.85`, `6.78`, `α`, `D`, `ε_mf`, `L_m`, `ε_m` all drop out. The
residual never leaves `[5.5e-9, 1.5e-8]` under any injection above — it is the
O(K_r) truncation from evaluating at `K_m = 1e-8`, not agreement.

**(iii) `ceiling_u0_spread` = 3.74e-16.** `K_bc L_f/u_b` contains no `u₀` at all,
so the spread is forced. The notebook prose derives this correctly (`:468-479`);
only `meta.yaml` mis-ranks it.

**What does have power, measured in the same runs:** the appendix comparisons.
`appendixC_max_dev` moves **4.61e-03 → 1.05e-01 (23×)** under the `K_bc 4.5 → 3.5`
injection. Those are the 0.46 %/0.57 %/1.8 % numbers `meta.yaml` ranks *below* the
identities.

**Also Shape 3(b).** The refinement study (`:658-682`, `upwind_order = 0.9993`)
runs `solve_distributed(convect=False)` exclusively — the degenerate branch where
two of three convective velocities are identically zero, the `ve < 0` recirculation
branch is never entered and the fixed-point loop never iterates. The page's *new*
result (`extension_max_change` = 0.1469, the "7 to 15 %" claim in *What pymrm
adds*) is the `convect=True` branch and has no refinement study — and
`build_page.py:784-786` explicitly borrows the frozen-branch convergence to
license it. Refined to n = 6400 the number is 0.14640 against the published
0.14686, so **the numbers are fine and the evidence for them is absent.**

**Failure scenario.** A reader takes "closure 0.001 %" as confirmation that
`K_bc = 4.5(…)` and `K_ce = 6.78(…)` were transcribed correctly from Kunii &
Levenspiel. A 22 % error in either leaves the verdict unchanged. The page's real
evidence is the appendix reproduction, which `meta.yaml` explicitly demotes.

---

### 2. `F3.5` — the penetration model is never time-refined, and its published headline is ~1 percentage point wrong. CONFIRMED

`pages/F3.5-co2-amine-absorption/build_page.py:1299-1360`. `pen_vs_film_max` is in
`agreement.json` (0.1048827394642553), is a **PRIMARY** validation bullet in
`meta.yaml:27`, and appears in `meta.yaml:47` under `adds` as *"the film-vs-
penetration sensitivity (up to ~10 %) the paper could not cheaply report"*.

The steady `Film` is grid-refined (`grid_delta`, n_x 400→800, 1.3e-5). The
`Penetration` subclass — backward Euler, `n_t = 160` on a geomspace grid — has no
refinement of `n_t`, `t0_frac`, `fac` or `n_x`, and it is the result carrying the
sharp front. Instead a `bias` measured on the *non-reacting* case is divided out
(`:1339`), on the assumption that the same relative time error occurs in the
reacting cases.

Measured by re-executing the page with only `n_t` changed (n_t = 160 reproduces
the published value to all 16 digits):

| `n_t` | quadrature bias | `pen_vs_film_max` |
|---|---|---|
| **160 (published)** | **+1.42 %** | **0.104883** |
| 320 | +0.72 % | 0.110322 |
| 640 | +0.37 % | 0.113011 |

The bias halves on every doubling — it is pure first-order backward-Euler time
error, not a quadrature property. Successive differences 5.44e-3 and 2.69e-3 give
an observed order of **1.01**; Richardson extrapolates to **≈ 0.116**. The
published 10.5 % is **10.6 % relative** from the converged 11.6 %, i.e. beyond
`check_agreement.py`'s 5 % tolerance had it ever been refined.

**Failure scenario.** A reader takes 10.5 % as the size of the film idealisation's
error. It is 11.6 %, and the individual unpromoted entry in the printed table
(0.61 %) is understated by roughly a factor of three. The *direction* is benign —
the difference is larger than claimed, so the conclusion that the reproduction is
not a film artefact survives — but the number itself carries ≈1 pp of error that
no check on the page can see. One `n_t = 320` run settles it.

---

### 3. `A2.3` — the page's only original claim rests on a single unrefined grid, and Reuse contradicts the result section. CONFIRMED

`pages/A2.3-taylor-aris-dispersion/build_page.py:453-504`. Every other result on
the page is refined (closure `n_r = 25…800` at `:378`, direct sim `n_z = 400…3200`
at `:397`). The "What pymrm adds" study that produces the τ threshold runs at
`s3 = Slug(PE_REF, n_z=1600, n_r=24)`, `dt = 0.0005`, once. **No metric from it
appears in `agreement.json`.**

| run | asymptotic `rate/target` | τ(±10 %) | τ(±5 %) | τ(±3 %) | τ(±2 %) |
|---|---|---|---|---|---|
| n_z=1600, dt=5e-4 (**the page's**) | 1.0223 | 0.140 | 0.175 | **0.195** | never |
| n_z=3200, dt=2.5e-4 | 1.0110 | 0.145 | 0.185 | **0.2125** | 0.2325 |

One refinement moves the ±3 % threshold by +9 % and **converts "never settles"
into 0.23 for the ±2 % band**. Which tolerances are attainable at all is set by
the mesh, not the physics: the +2.2 % numerical-dispersion floor makes
`|rate/target − 1| < band` satisfiable earlier, so the coarse grid biases the
settling time *low* and the reported τ values increase monotonically under
refinement.

Two secondary items on the same page:
- the caveat at `:503` quotes **+1.8 %**, which is the interval-averaged error of a
  *different* run (`dt=0.001`); the bias that actually contaminates the band test
  is the instantaneous asymptote, **+2.23 %**, leaving 0.8 % of margin on a ±3 %
  band;
- Reuse (`:532`) says "the τ ≳ 1 result" where the result section (`:506`) says
  τ ≳ 0.2. **Internal contradiction inside Reuse advice.**

**Failure scenario.** A reader takes τ ≳ 0.2 as the physical time after which the
lumped Taylor–Aris coefficient is defensible, and finds Reuse telling them τ ≳ 1.
The true converged threshold is higher than 0.2, and the ±2 % band the page
reports as unattainable is attainable.

---

### 4. `F2.3` — a printed exponent is provably inert on the reactor path, and Reuse advertises the function that contains it. CONFIRMED

`pages/F2.3-slurry-bubble-column-ft/build_page.py:213-218`:

```python
def hydrodynamics(U, eps_s, d_t=DT, rho_g=RHO_G):     # DT = 7.0 m  (line 159)
    ...
    eb = 0.3 * min(d_t, 1.0) ** -0.18 * max(U - udf, 1e-12) ** 0.58
```

`hydrodynamics` is called exactly once (`:257`) and with no `d_t` argument, so
`min(d_t, 1.0)` ≡ `min(7.0, 1.0)` ≡ `1.0`, and `1.0 ** -0.18` ≡ `1.0` **exactly,
for any exponent**. The `-0.18` can be deleted outright without moving
`conversion_U012`, `conversion_U040` or any printed number — bit-identical.

This is `A1.8`'s degenerate point recurring verbatim: an evaluation where
`1**x = 1` makes an exponent algebraically invisible. It is the second independent
instance in the repository.

The second copy of eq. 9 (`:303-306`, `eps_model`, hard-coded `0.10 ** -0.18`) is
the one compared against the 79 digitised markers — but at a single column
diameter, so what the data constrain is the product `C·0.10^-0.18 = 0.454`. The
`C = 0.300` vs `0.268` study tests that product, not the exponent. **So `-0.18` is
untested on one code path and inert on the other.**

**Reach: Reuse advice.** *"`hydrodynamics` is the closure other pages want. It
takes (U, ε_s) and returns the two-class split"*, and *"To adapt this to your
system … The reactor geometry is three constants at the top."* A reader who
changes the geometry to a 0.10 m column moves into the live branch, where the
factor becomes `0.10**-0.18 = 1.5136` — a **51 % change in ε_b**, governed by a
constant no printed number on this page can distinguish from any other value.

---

### 5. `A4.9` — two of four checks are one Newton-residual identity, blind to the diffusivity that carries the headline. CONFIRMED

`pages/A4.9-duncan-toor/index.ipynb`, cell 18. Prose (cell 17): *"Four checks. The
first three do not use the experimental data at all."*

Check 1 (`flux_nonuniformity`) recomputes `flux(x_ss, …)` — the same function whose
divergence `solve_capillary` drove to `tol=1e-12` with `construct_div(nu=0)`. On a
Cartesian grid that residual **is** `(N_{i+1} − N_i)/Δz = 0`, so check 1 is the
converged Newton residual re-divided by a scale. Check 2 (`mole_balance_error`) is
zero **iff** the flux is uniform, i.e. iff check 1 passes. Two of four checks are
one identity.

Measured by re-executing the page's own cells (baseline reproduces the published
metrics exactly):

| injection | `flux_nonuniformity` | `mole_balance_error` | mean dev vs data | N₂ crossing | final Δx(N₂) |
|---|---|---|---|---|---|
| **as published** | 1.577e-14 | 2.163e-16 | 0.59 % | 0.960 h | +0.0849 |
| `D(N₂–CO₂)` 16.8 → **61.8** e-6 (digit swap) | 1.927e-14 | 4.327e-16 | 2.28 % | **3.240 h** | **+0.0043** |
| `D(N₂–CO₂)` 16.8 → 1.68 e-6 (lost decade) | 2.155e-14 | 4.327e-16 | 1.30 % | 0.880 h | +0.2090 |
| `D(H₂–N₂)` 83.3 → 8.33 e-6 (lost decade) | 1.678e-14 | 1.731e-16 | 5.19 % | **never** | −0.1106 |
| `LA_EFF` 258.1 → 285.1 e2 | **1.577e-14** | **2.163e-16** | 0.80 % | 1.040 h | +0.0934 |
| `vol_2` 78.63 → 87.63 e-6 | **1.577e-14** | 3.273e-16 | 0.73 % | 1.000 h | +0.0902 |

A digit transposition in the one diffusivity that drives the whole phenomenon ends
the N₂ reversal **20× smaller and 2.3 h later**, and neither check moves off
machine precision. `LA_EFF` and `vol_2` slips leave check 1 *bit-identical*.

**What does have power:** check 4 alone, the comparison against the 28 digitised
Fig. 2 points — it moves 0.59 % → 2.28 %/1.30 %/5.19 % on the three diffusivity
slips. The page should say that check 4 is the only guard, rather than presenting
four.

**Same defect, second page.** `A4.2` carries the identical block —
`flux_nonuniformity` = 1.94e-14 and `mole_balance_error` = 2.16e-16 at
`build_page.py:1072-1080` — with the same algebra and the same blindness (a
10× error in `D(H₂–N₂)` flips `J2/J1` from −0.4675 to +0.3899 while flux
non-uniformity stays at machine zero). These are the *only* two pages in the
repository with those two metric names. `A4.2` borrows `A4.9`'s dataset; the
validation block travelled with it. **Impact on `A4.2` is incidental** — no
headline depends on it, and its `3b` identity *is* correctly labelled.

---

### 6. `J3.4` — the check ranked first of six detects none of the three error classes its own prose names. CONFIRMED

`pages/J3.4-doyle-fuller-newman/build_page.py:693-699` (prose) and `:708-712`
(code). `salt_conservation_drift` = 3.33e-16 is in `agreement.json` and is ranked
**first of six "in order of what they can catch"**:

> "1. **Salt conservation.** … Any sign error in the migration term, the
> transference number or the reaction coupling breaks it."

It is a telescoping identity. With `r_c = eps*(c−c_o)/dt + div @ n_salt` (`:493`)
and `n_salt[0] = n_salt[-1] = 0` (`:487`), `Σ_i dx_i (div n)_i = n[N] − n[0] = 0`
exactly, so `Σ eps·dx·c` is constant for **any** interior flux expression and any
parameters.

| injected defect | salt drift | V(u = 0.5) |
|---|---|---|
| as published | 3.33e-16 | 1.984998 V |
| migration sign flipped in `n_salt` | **4.44e-16** | 2.164254 V (+179 mV) |
| `(1−t⁰₊)` written as `t⁰₊` in `i2` and `n_salt` | **3.33e-16** | 2.091128 V (+106 mV) |
| salt diffusivity D 10× too small | **3.33e-16** | march dies at u = 0.19 |
| reaction coupling `ai` sign flipped | **3.33e-16** | does not discharge at all |

**All three error classes the prose names are invisible.** `J3.1` (`:1830-1832`)
and `J3.4`'s own *Reuse* section (`:957-960`) both describe this correctly as
structural; only the Validation ranking overclaims.

**Also Shape 3.** `Cell.march` (`:554`, `:573`) fixes `dt, dt_max =
2e-5·t_ref, 6e-3·t_ref` growing ×1.3 — a schedule that is a function of current
alone and **completely decoupled from `n_s`, `n_c`**. So the "numerical
independence" study (`:714-722`) holds the temporal error constant by
construction, and the page has no bound on it anywhere. Measured: the dt error is
~1.8× the grid error (0.256 mV vs 0.149 mV at I = 10) and still drifting
monotonically; at I = 20 it moves `u@1.7 V` by 9.5e-4 against 3.9e-4 for the grid.
Check 2 also reads V at u = 0.5, on the flat plateau, whereas every headline is
read at a voltage crossing on the near-vertical collapse.

**Same root cause, sibling page.** `J3.1` shares `Cell.march` (`:1373`) and the
same fixed schedule; its `grid_V_spread_mV`, `grid_mean_eta_spread_frac` and
`grid_peak_eta_ratio` are all in `agreement.json` and all hold dt constant. The
peak overpotential grows under **both** knobs (1.90× grid, 1.40× dt), and the
page's margin argument counts only the first. **`J3.5` fixed this locally** by
adding `dt_scale`, which scales `dt0` *and* `dt_max` — the correct fix — and the
two siblings did not inherit it.

**Impact:** moderate-low on conclusions (Δu ≈ 6e-4 against a 25.5 mV model–figure
gap), high on the check ranking.

---

### 7. `J1.5` — both refinement studies sit at the one τ the headline is not taken at. CONFIRMED

`pages/J1.5-ldf-breakthrough/build_page.py:270-284` (check 1, grid) and `:286-290`
(check 2, time step) both use `target = np.array([0.05])`. The reported
`solver_worst_rel_err` = 2.6465e-3 (`agreement.json`) is a max over `TAUS`
spanning 0.001…0.5 and is attained at **τ = 0.001** — the thin-shell regime the
page's whole argument is about (`ldf_max_rel_err_early` = 0.955) — which neither
study visits.

| τ = 0.001, dt = 2e-5 | n_r=25 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|
| rel err | 2.74e-2 | 4.03e-3 | 2.87e-3 | **2.646e-3** | 2.589e-3 | 2.574e-3 |

Grid refinement **saturates**: the error is essentially 100 % temporal. Halving dt
at n_r = 200, τ = 0.001 gives clean first order and jointly refining moves the
worst case 2.65e-3 → 6.5e-4.

Secondary (3c): `:287` prints "**time-step independence**" while the τ = 0.05
ladder falls 6.4e-4 → 3.3e-4 → 1.7e-4 → 7.6e-5, still ∝ dt and eightfold from
converged; and "the scheme should be second order in dr" (`:271`) is asserted with
printed ratios 3.21/2.38/1.55/1.16. Re-run at dt = 1e-6 the ratios are
**3.80/3.76/3.42** — the claim is true, just not established by the code as written.

**Impact: low.** Every *scientific* headline (`ldf_worst_abs_err` 0.170,
`ldf_max_rel_err_early` 0.955, `long_time_decay_constant` 9.8696) is analytic and
never touches the solver. Only the two solver metrics are affected.

---

### 8. `C2.1` — a fourth shape: the headline is a max over a set whose baseline scatter exceeds the defect. CONFIRMED

`pages/C2.1-xu-froment-smr/build_page.py:396-417`; `param_round_trip_worst_pct` =
1.4432 is in `agreement.json`, in `README.md:31-32` and in `meta.yaml:24`.

The Table 5 ⇄ Table 6 round trip is genuinely three-source and is **not** Shape 1.
But the published number is `rt.rel_pct.max()`, and the baseline per-row scatter
already runs to 1.4432 % (driven by `k3`), so any single-parameter slip landing
below that ceiling is invisible in the headline:

| perturbation | printed "worst deviation" |
|---|---|
| baseline | 1.4432 % |
| `K_CO` A `8.23e-5 → 8.32e-5` (digit transposition, +1.09 %) | **1.4432 %** (row goes 0.276 → 0.809 %) |
| `k2` A `1.955e6 → 1.965e6` (+0.51 %) | **1.4432 %** (row 0.290 → 0.798 %) |
| `k1` A `4.225e15 → 4.226e15` | **1.4432 %** |
| `k1` E `240.1 → 240.2` | **1.4432 %** |
| `k1` E `240.1 → 241.1` | 18.85 % |
| `k1` A lost decade | 90.13 % |
| `K_CH4` `T_ref 823 → 648` | 77.90 % |

**Can detect:** lost decimals, decade errors, the split reference temperature, E
slips ≳ 0.5 kJ/mol. **Cannot detect:** a digit transposition in *any* `A`, or a
last-digit slip in any `E`. The per-row table is printed in the notebook, so the
movement is visible to a reader — the defect is in the headline and the metric,
not in the page's own output.

This is a genuinely distinct shape from the three swept for, and it is worth
naming: **a max/worst-case summary is only as sharp as its own baseline scatter.**

Two smaller items on the same page. `K1 = exp(−26830/T + 30.114)` and
`K2 = exp(4400/T − 4.036)` are hard-coded (`:246-258`) and check 2 (`:419-434`)
compares only the **slopes** against Table 7 — a `+1` slip in the K1 intercept
multiplies K1 by e = 2.718 at every temperature and check 2 is unchanged. And the
Reuse advice (`:628-634`) tells the reader to scale `k1,k2,k3` by **1.225** and
**2.246** for the other two activity levels; neither number is used or checked
anywhere on the page.

---

### 9. `B1.1` — a check that restates a definition, and a second solver that only ever visits one β. CONFIRMED

`pages/B1.1-thiele-weisz-hicks/index.ipynb`, cell 22, billed as one of "four checks":

```python
# 4. Prater relation: the ignited branch should approach the bound beta.
dT_max = BETA_M * (1.0 - y_min_ignited)
```

There is **no temperature field anywhere on this page.** `wh_rate` (cell 15) has
the Prater relation already substituted, so ΔT/T_s ≡ β(1−y) *by construction*, and
`newton(..., callback=clip_approach(x, g, 0.0, None))` enforces y ≥ 0. The printed
number is β(1−y_min) ∈ [0, β] for **any** β, **any** γ, **any** grid and **any**
number of Newton iterations. It cannot fail. Not in `agreement.json` — prose reach
only.

**Shape 2.** Cell 16 draws the whole Weisz–Hicks diagram for
`BETAS = [0.0, 0.2, 0.4, 0.6, 0.8]` using the shooting reference alone, including
the printed per-β `eta_max` and `multivalued in phi` verdicts; pymrm runs at
`BETA_M = 0.6` only. Measured `eta_max` = 1.00 / 2.32 / 10.73 / 44.55 / 151.22.
Only the fourth is ever seen by a second method. *Note*: the `agreement.json`
headlines (`eta_ignited`, `phi_fold_lower/upper`) are taken at β = 0.6, which
check 3 *does* cross-validate, so the exposure is prose-level, not headline-level.

**Separate finding on the same page, and the more serious one.** The YAML front
matter (cell 0) describes the page as *"where η exceeds 1000"* and the intro says
*"a thousand times faster than its own surface conditions"*. **The largest η the
page computes anywhere is 151.22**, at β = 0.8, γ = 20 (measured above). No
printed number on the page supports the description that the gallery listing
carries. Weisz–Hicks reach η ~ 10³ at γ = 30, which this page never runs.

---

### 10. `A4.4` — a scaling check that is identically zero for any prefactor, and the prefactor reaches Reuse. PLAUSIBLE→CONFIRMED

`pages/A4.4-knudsen-bosanquet/build_page.py:404-408`, `:429-431`;
`eq85_scaling_dev` = **0.0** in `agreement.json`. The metric is built from three
*ratios* of `knudsen_D` to itself, so the prefactor `EPS_TAU/KN_DEN` cancels
identically:

| `KN_DEN` | `eq85_scaling_dev` | `D_K(100 nm, H₂)` |
|---|---|---|
| 3.0 (printed) | 0.000e+00 | 6.02e-05 |
| 2.0 (d₀ read as a radius) | **0.000e+00** | 9.03e-05 |
| 30.0 (factor 10) | **0.000e+00** | 6.02e-06 |

Nothing else on the page can see it either: every result is reported in the
dimensionless group `Kn = D_AB/D_AM` (`:679`). The unguarded prefactor does reach
print — the "Kn = 1 for H₂ in H₂–N₂: d₀ = … nm" line (`:402`), the `D_K(d₀ = 100
nm)` values (`:393`), the abscissa of the Bosanquet transition figure — and it is
load-bearing in the **first row of the Reuse table** (*"Your own pore size,
temperature and gas | `knudsen_D(d0, T, M)`; then `Kn = D_AB/D_K` decides whether
any of this matters"*).

The page's otherwise-exhaustive accounting of borrowed inputs (`:1155-1219`)
enumerates T, p and D_ij but **omits `KN_DEN` and `EPS_TAU`**. The check is not
vacuous — it does have power over `KN_EXP` and over the T- and d₀-linearity — only
prefactor-blind. `A4.3` measures and declares exactly this class of blind spot
(`:1508-1578`, *"NOTHING ON THIS PAGE PINS THE d_0/3"*); `A4.4` does not.

---

## Incidental and self-disclosed (13 pages)

| page | what it carries | why it is not ranked higher |
|---|---|---|
| `A1.1` | Shape 1 — check 6 prints `\|f_v(0) − k1\| = 1e-14` for `f_v(x) = k1 + k2·x` | line itself reads "exact by construction"; not in `agreement.json` |
| `A1.7` | Shape 1 — `eps0_identity_residual` = 5.6e-17 on rows selected *because* `H/H₀ == 1.000`, so the expression collapses to `eps_MB` identically. The break-table row that appears to give it power uses a selector frozen from the *unperturbed* table | page states it at `:846-850`: "structural … tests only that those two cells were transcribed as 1.000" |
| `A4.2` | Shape 1 — the copied `flux_nonuniformity` / `mole_balance_error` block (see finding 5) | no headline depends on them; the `3b` identity *is* correctly labelled. Separately, `n_z = 40` is never refined — I refined it: J2/J1 moves 3e-5 relative against a 6.7 % claim, benign |
| `B1.2` | Shape 3(c) — "with the expected **second-order** refinement"; the observed order is **2.81/2.85/2.88/2.89** | the claim is *conservative*, so no published number is wrong |
| `B1.6` | transient (`N_T = 100`) has no refinement study while the steady results are refined 200/400/800 | I ran it: N_T = 50→400 moves `transient_worst_violation` by 4e-4 relative; the published 0.2519 is right to four digits |
| `D2.1` | Shape 2 — `Tube.__init__`'s `simplified=True` branch (`:989-993`) is never instantiated and can be deleted | not in `agreement.json`; the simplified-law headline comes from `profile_ivp`/`blowup_z`, which *are* cross-checked LSODA vs Radau |
| `D2.2` | Shape 1 — `T_M² − a(T_M − T_w)` is identically zero for any `a`, `T_w`, since `T_M` is that root; plus 3(c) | annotated "(exact, from Eq. 8)"; not in `agreement.json`. The grid ladder runs at `p0 = 0.0164` rather than at `p_crit`, and `n_z = 1500` is not one of its levels — but check 4's grid-free phase-plane route bounds it at 0.182 % |
| `F1.3` | Shape 1 ×2 — `continuity_limit_departure` is `V_small/V_b → 1` from the same call chain (a 10× change in eq. 4's prefactor leaves it at 1.1e-5); plus a `brentq` round trip | the notebook says exactly this ("This is a structural identity"); the issue is that `agreement.json` ships the number bare beside seventeen that do have power |
| `G1.7` | Shape 1 ×2 — check 4 is literally `x == x` (`:525-529`, `rhs` is `lhs` with `Re80` inlined, printing 0.0e+00); plus an eq. 3 round trip | both self-labelled "(identity, as implemented)"; neither is published |
| `H1.7` | **not a shape — a wrong reported number.** Prose says a 50 µm film raises wall concentration 25 % and a 200 µm one 76 %; the page's own `film_wall_ratio` gives **+31 % and +194 %** at the plotted flux. `meta.yaml` says 23 %/70 %, so one of the two is also stale | reaches the `adds` block that feeds the gallery listing; the 200 µm figure is understated 2.8×, and the page's own left panel plots the correct 2.94 profile beneath the sentence |
| `I1.2` | **not a shape — a broken data path.** `PAGE = "I1.2-oh-cavendish-converter"` but the directory is `I1.2-catalytic-converter-lightoff`. The raw URL **404s**; the correct one returns 200. Every `load_data(..., page=PAGE)` fails on a fresh Colab VM, which is exactly what the `AGENTS.md` rule exists to prevent. The only such mismatch in the repo; `check_metadata.py` does not validate the constant | clean on all three shapes. Also: extending the case-2 grid ladder to n_x = 640 gives an observed order of **−0.74** on the finest pair (peak sampling snaps to cell positions), so "converged to about 1.5 K" is a defensible *spread bound* but not a demonstrated convergence — and the 1.5 K ≪ 5 K argument survives |
| `J3.1` | Shape 3 — shares `J3.4`'s dt-decoupled `Cell.march` (see finding 6). One mild redundancy: the arcsinh inversion (`:619-632`) is evaluated only at symmetric α, inside the region the page's own defect table shows is blind | Shapes 1 and 2 are already *measured and published* by the page — `defect_worst_at_symmetric_alpha_eq17_vs_eq30` = **0.0** names Shape 2's degenerate point outright. Headline survives (margin 15.7× → ~11×) |
| `J3.5` | Shape 1 — `coulomb_identity_dev` = 6.7e-16 presented as "Marquis et al.'s Eq. 32 … holding exactly"; it is a discrete conservation identity of the DFN residual (`du/dt ≡ I/Q` for any parameters, grid or dt) | filed under "Structural conservation"; nothing is concluded from it. Fix is one word of framing. The prior documented instance is fully disclosed, and `dt_scale` is the correct Shape-3 fix |

## Known instances — all three confirmed repaired

- **`A1.8`** (Shape 2). `EPS_PROBE = [0.45, 0.50, 0.60, 0.70, 0.85, 0.95]` and a
  swept break table now recover both deletable constants; `:1119-1145` states the
  defect openly rather than patching quietly. I reproduced the historical
  behaviour: at ε = 0.50 the dilute-branch perturbations move the bias by **exactly
  zero**, and my −15.338 / +17.487 match the published
  `check3_break_dilute_exponent_2p65_to_1p65_pts_at_eps095` and
  `check3_break_switch_0p85_to_0p55_pts_at_eps060` to all printed digits.
- **`A4.7`** (Shapes 1 and 3). The `dt` schedule is now
  `step = min(dt_max, max(dt0, dt_frac·fo), target−fo)` with `fine(n_r, s)` scaling
  all three of `dt0`, `dt_frac`, `dt_max`; `:420-424` names the old defect.
  `GAMMA_IDENTITY` is retained but explicitly labelled "**not** a guard here" and
  replaced by an injection table that moves the survey's conclusions
  (`extrema_survey_verdicts_moved` = 1.0, `worst_shift` = 0.265) while the identity
  stays at ~1e-7.
- **`H1.9`** (Shape 3(b)). Grid study is now under the transient, three-level dt
  study, observed orders computed, Richardson extrapolation, and the old coarse
  values retained as separate metrics beside the production ones (−0.81 % → −9.82 %,
  and the sign flip on the second metric). The withdrawn "backward Euler damps a
  maximum" explanation is deleted rather than reworded. A suspicion about an
  unreported occupancy-clamp branch was instrumented and found to fire **zero**
  times in the reported runs.

One residual on `H1.9`: the validation cell's heading reads *"[Gamma] and [B],
each by two routes that share no code"* and `gamma_numerical`'s docstring calls
itself *"The independent route"* — but `gamma_numerical` central-differences
`iast.inverse` on the same `mix` object that `inverse_and_gamma` belongs to, so for
Γ the two routes share everything. The cell body says so correctly two lines below
("both routes call the same IAST solve"). The heading contradicts its own body; the
`[B]` half of the claim is accurate (the GMRES route never forms `[B]`).

## Pages clean on all three shapes (15)

`A1.5` `A1.6` `A2.1` `A2.5` `A3.4` `A4.3` `B1.4` `B3.1` `E1.2` `F1.4` `F3.1`
`G1.8` `H1.1` `H1.4` `J4.8`

Four are worth naming as the standard the rest should be held to, because each
one *finds and publishes its own blind spot*:

- **`A4.3`** breaks the `d₀/3` prefactor in every route simultaneously and prints
  "NOTHING ON THIS PAGE PINS THE d_0/3" (`:1508-1578`). This is the exact defect
  `A4.4` carries unlabelled.
- **`G1.8`** titles a check "**The collapse — a check that cannot fail**",
  demonstrates it by corrupting the sphere term ×1.01/1.5/3.0 and showing the
  identity does not move, and reports that **4 of its injected defects would not be
  rejected** — one of them scoring *better* than the undamaged chain. It also
  evaluates at `1.0 ** 0.3` deliberately, to *measure* that at L_m ≈ 1 the second
  term carries only 13 % of the resistance and check 4 is "essentially blind" to
  it. That is Shape 2's degenerate point, found and published by the page itself.
- **`B1.4`** publishes its blind spots as named metrics —
  `Phi_identity_blind_wrong_nu` = 4.2e-12 while η is 26 % wrong,
  `Phi_identity_blind_three_cells` = 2.1e-15 — alongside what the identity *does*
  catch.
- **`A2.5`** publishes `parameter_chain_D_L_shift_unseen` = 0.069: all three
  diagnostics identical to five decimals under the wrong Table-1 row.

Two clean pages deserve a specific note. **`H1.1`**'s α_H check has a genuine
circularity — `r_i` was *selected* by requiring the identity to close — and the
code says so outright ("the check chose them, so it cannot catch them being
wrong"); because the selection is over a **discrete** pair of readings 2.4 % apart
against a 0.03 % residual, the check keeps its power over every other input.
**`H1.4`**'s headline sits on an algebraic ceiling and the page prints *"'Nothing
fitted' is true, but almost nothing could have made it fail"*, then measures it
(k·V_r over ×100 moves X by exactly 0) and switches to the quantity that varies.

## Which shapes turned out to be worth a rule

| shape | new instances | headline / Reuse reach | verdict |
|---|---|---|---|
| **1** — differentiates the object it tests | ~12 pages | **4** (`E2.1`, `A4.9`, `J3.4`, `A4.4`) | by far the most common; needs a rule |
| **2** — all checks at one point on a branched model | 4 (`F2.3`, `D2.1`, `B1.1`, + `A1.8` known) | **1** (`F2.3`) | rarer, but two of the four are the *identical* `x**exponent where x ≡ 1` degenerate point — worth one line, not its own rule |
| **3(a)** — geometric `dt`, refining `dt0` | **0** | 0 | **not worth a rule.** Already swept repo-wide on 2026-08-02 (commit `fb388f3`); this sweep confirms zero instances, and `J3.5` independently implemented the correct fix (`dt_scale` scaling `dt0` *and* `dt_max`) |
| **3(b)** — steady refined, transient not | 3 (`A2.3`, `F3.5`, `E2.1`, + `H1.9` known) | **3** | the highest-yield shape per instance; every instance reached a headline |
| **3(c)** — convergence claim with no observed order | 13 pages run a refinement study and print no order | 1 (`J1.5`) | mostly benign. Where I checked the assumed order it was **right**: `C2.1` assumed p = 1 and measures 1.001/1.002/0.998/1.000; `D2.2` assumed p = 1 and measures 0.983/0.991/0.996/0.997. Fold into the rule below rather than legislate separately |

**The shape not on the list.** `C2.1` (finding 8) is a fourth shape: a worst-case
headline whose own baseline scatter exceeds any plausible single-digit slip. It is
not Shape 1, 2 or 3, and nothing currently looks for it.

**And the transmission vector is inheritance, not invention.** Two of the ten
headline findings are one copied validation block (`A4.9` → `A4.2`) and two more
are one shared helper (`J3.4` / `J3.1`'s `Cell.march`, with `J3.5`'s fix not
inherited). `AGENTS.md` step 3 instructs builders to *"copy that page's directory
and substitute the physics"* — which is right, and which is exactly how a check
whose power depended on the *old* physics arrives as decoration in the new page.

## The one change

Add to `docs/agent-brief.md`, under **Make sure your validation can fail**,
replacing the three existing bullets with one requirement:

> **Every number you put in `agreement.json` needs a row in a break table.**
> For each published metric, inject one defect that metric is supposed to catch,
> re-run, and print the before/after beside it. A metric with no such row is not
> evidence — label it `structural` and say in the same sentence what it cannot
> detect, as `A4.3`, `B1.4` and `G1.8` do.
>
> Three cheap ways to find the row that matters:
> - **Perturb a printed constant by one digit**, not by a factor. A digit
>   transposition is the error that actually happens, and a worst-case summary is
>   only as sharp as its own baseline scatter — if your headline is a `max()`,
>   check that one perturbed row can exceed it.
> - **Evaluate at more than one point on a branched model**, and check no
>   exponent sits where its base is 1, 0 or a fixed constant.
> - **Refine every reported result at the condition it is reported at**, and print
>   the observed order from the sequence. If a page has both a steady and a
>   transient result, the transient needs its own study.
>
> **When you copy a page's directory, the break table does not travel — rebuild
> it.** A check's power depends on the physics, not on the code. `A4.9`'s flux
> uniformity is its Newton residual restated and it arrived on `A4.2` unchanged;
> `J3.4`'s `Cell.march` fixes `dt` as a function of current alone, so every grid
> study built on it holds the temporal error constant, and `J3.5`'s fix for that
> was never inherited by `J3.1`.

**Why this one.** Nine of the ten headline-reaching findings are on pages with no
defect-injection table; ten of the fifteen clean pages have one. The rule is the
single strongest predictor in the data, it subsumes 3(c) (an observed order is
just the refinement study's break row), and it converts the currently-invisible
class into something CI can see — today, 74 of 920 metrics sit below
`check_agreement.py`'s `ABS_FLOOR = 1e-12` and are not compared at all.

Two smaller changes worth making at the same time, neither a judgement call:

1. **`scripts/check_metadata.py` should assert `PAGE == <directory name>`.**
   `I1.2` has been shipping a `PAGE` constant that 404s on Colab, and no gate
   catches it.
2. **`meta.yaml`'s validation ranking is a claim and should be checked like one.**
   `E2.1` ranks three algebraic identities above the appendix comparisons that are
   its only real evidence — and that ranking, not the notebook prose, is what a
   reader filtering the gallery sees.

## Dispatch order

| priority | page | fix |
|---|---|---|
| 1 | `F3.5` | re-run at `n_t ≥ 320`; the published `pen_vs_film_max` is ~1 pp low and beyond CI tolerance |
| 1 | `E2.1` | re-rank `meta.yaml`; label the three identities; add a refinement study in the `convect=True` branch |
| 2 | `A2.3` | refine the τ study; fix the +1.8 % caveat; reconcile Reuse's τ ≳ 1 with the result's τ ≳ 0.2 |
| 2 | `F2.3` | the `-0.18` exponent is inert and untested — say so, or test it against `F1.4`'s multi-diameter data |
| 2 | `A4.9` | say that check 4 is the only guard; label checks 1 and 2 structural. Same label on `A4.2` |
| 2 | `J3.4` | demote `salt_conservation_drift` from first of six; add a dt study (port `J3.5`'s `dt_scale`), and to `J3.1` |
| 3 | `B1.1` | the "η exceeds 1000" description is unsupported — the page's max is 151 |
| 3 | `H1.7` | the 50/200 µm polarisation numbers are wrong in prose *and* differently wrong in `meta.yaml` |
| 3 | `C2.1` | print the per-row max alongside the headline, or report the median |
| 3 | `A4.4` | add `KN_DEN` / `EPS_TAU` to the borrowed-input accounting, `A4.3`-style |
| 3 | `J1.5` | refine at τ = 0.001; drop the unearned "time-step independence" wording |
| 4 | `I1.2` | fix the `PAGE` constant; add the assertion to `check_metadata.py` |
| 4 | — | one-line labels: `A1.1`, `A1.7`, `D2.2`, `F1.3`, `G1.7`, `J3.5`, `H1.9`'s heading |
