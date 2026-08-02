# A2.1 — The Danckwerts boundary conditions (covers A2.2)

A tubular reactor with axial dispersion needs two boundary conditions, and the
obvious pair is wrong. Danckwerts (1953) gave the right pair: a **flux** balance
at the inlet, which makes the concentration drop across the bed entrance —
discontinuously in the closed-vessel idealisation — and a **zero gradient** at
the outlet. Wehner & Wilhelm (1956) derived both by putting the bed between an
inlet pipe and an outlet pipe, and in their fuller model the inlet drop is
continuous: they say of Danckwerts' discontinuity that it is *"not necessarily
correct"*, and the page quotes and reconciles that.

- **Structure:** `S4` (1-D convection–dispersion–reaction)
- **References:** Danckwerts (1953), *Chem. Eng. Sci.* 2(1) 1–13,
  doi:10.1016/0009-2509(53)80001-1; Wehner & Wilhelm (1956), *Chem. Eng. Sci.*
  6(2) 89–93, doi:10.1016/0009-2509(56)80014-6
- **Runtime:** ~5 s
- **Data:** none. Tier 6, analytical; no figure on either paper is digitised.

## Why one page for two catalog entries

Wehner & Wilhelm's eq. (20) for the reaction section is algebraically identical
to Danckwerts' eq. (33), and they say so on their page 91. A separate `A2.2`
page would have to validate against the same closed form. What they add is the
*justification* of Danckwerts' conditions — Danckwerts himself wrote "Intuition
suggests that neither of these situations can arise" — and a boundary condition
and its justification are one result. The full argument is in
`queue_cases/A2.1.yaml` under `scope_decision`.

## The conditions, in pymrm's convention

`bc` is `a·∂c/∂n + b·c = d` with `n` the **outward** normal, so the same dict
means different physics at the two ends.

```python
D, u = 1.0 / Pe, 1.0
bc = (
    # inlet, Danckwerts eq. 31:  u c* = u c - D dc/dz  at z = 0.
    # n = -z, so dc/dn = -dc/dz  ->  (D/u) dc/dn + c = c*
    {"a": D / u, "b": 1.0, "d": c_star},
    # outlet, Danckwerts eq. 32:  dc/dz = 0 at z = L.
    # n = +z, so dc/dn = dc/dz   ->  a = 1, b = 0, d = 0
    {"a": 1.0, "b": 0.0, "d": 0.0},
)
```

## Agreement

| Check | Result |
|---|---|
| four transcriptions, two papers (D. eqs. 33/34, W&W eqs. 20/21/22) | agree to **1.1e-16**; one wrong digit lifts it to 2.9e-02 |
| pymrm vs Danckwerts eq. 34, van Leer deferred correction | **2.2e-07** at n = 1600, observed order **2.00** |
| the same with bare upwind | order **1.00**; the residual against `Pe_eff = 1/(1/Pe + 1/2n)` is order 2.00 |
| eight decades of Peclet number vs eq. 34 | worst **1.9e-04** (at Pe = 1e4, cell Peclet 6) |
| eq. 35 plug-flow limit / eq. 36 stirred-tank limit | 4.5e-08 / 2.2e-09 |
| eq. 37 asymptote, residual decay exponent | **2.00** over Pe ≥ 1e3, i.e. `R²/Pe` is the exact first correction (1.95 if the pre-asymptotic end is included) |
| W&W three sections, **no Danckwerts BC imposed anywhere** | reproduces the closed vessel to **8.0e-06**, invariant over Pe_a across five decades and Pe_c across four |
| harmonic vs arithmetic face average of `D` at the jump | an **order**, not a factor: 1.00 against 2.00, so the ratio is 84× / 172× / 347× / 696× at `n_b` = 100 / 200 / 400 / 800 |
| naive Dirichlet inlet vs its own closed form (Hulburt = eq.34 / eq.22) | 7.1e-05 over four decades |

## The defect table

| injected defect | rel. deviation from eq. 34 |
|---|---|
| none | 8.9e-07 |
| inlet BC: sign of `a` flipped | **1.97** |
| inlet BC: naive Dirichlet `c(0) = c*` | **0.497** |
| inlet BC: `a = D` with D taken as 1 | 0.356 |
| **outlet BC: sign of `a` flipped** | **8.9e-07 — bit-identical matrix and RHS** |
| outlet BC: Dirichlet `c(L) = 0` | 0.998 |
| `construct_div nu = 1` | 1.00 |
| reaction sign flipped | 26.1 |
| n = 10 cells | 0.006 (but 0.041 at Pe = 200) |

The outlet row is the point: because `b = 0` and `d = 0` there, the sign of `a`
is immaterial, so the outward-normal error is **undetectable at that end by any
test**. Not "the number does not move" — the page assembles the system for
`a = +1`, `-1` and `1e7` and finds max |ΔA| = max |Δb| = **0.0**, against 1.8e+03
for a flipped inlet sign. Set `d = 0.05` and the outlet sign matters again
(`f(exit)` 0.24692 vs 0.22193), so the because-clause is doing the work.
Verifying an outlet condition proves nothing about an inlet condition.

## What the page adds

The size of the mistake. `c(0) = c_in` over-predicts the outlet concentration by
179 % at Pe = 0.1, 50 % at the Peclet number Wehner and Wilhelm printed on their
own Figure 4, and 1 % by Pe = 200 — and as Pe → 0 it predicts *no conversion at
all* where the answer must be a stirred tank. Both papers had that algebra;
neither plotted it.

Also: bare upwind adds a numerical diffusivity `u·Δz/2` that is
indistinguishable from a lower Peclet number, so a 200-cell grid at Pe = 200 is
really running at Pe_eff = 133. A boundary-condition study on an upwind grid can
measure its own truncation error and call it physics.

## Honest limits

Tier 6, not experimental. Neither paper measures anything that tests these
boundary conditions; Danckwerts' one experimental figure belongs to his
*open*-vessel residence-time analysis, which he says himself cannot be carried
into the reactor problem. The demonstration that `df/dz = 0` emerges at the bed
exit reproduces Wehner & Wilhelm's chain of reasoning, not an independent proof
of its premise — imposing boundedness at the far downstream truncation is itself
a zero-gradient condition, applied one section away. Nor is the three-section
route independent *code*: it shares the whole operator stack and the solver with
the closed-vessel solve, and its condition at `z = 1` is discretely identical to
`bc_out = {a:1, b:0, d:0}`.

Two printed residuals are **structural** and labelled as such: the fore-section
flux being constant (only its *value* of 1.0 informs, as a domain-truncation
check) and the after section being flat, which the discrete equations force for
any `Pe_c` and any truncation length — cutting the after section from 20 decay
lengths to 0.2 leaves it just as flat and the bed unchanged to 1.8e-12. The
recovered `df/dz(1)` is 84 % the truncation error of the 3-point fit that
recovers it, measured by applying that fit to the exact eq. 33.

The residence-time half of Danckwerts 1953 — F- and C-diagrams, hold-back,
segregation, the laminar-pipe distribution — is not built here. `A2.4` is the
catalog entry for it.

## Rebuilding

```bash
python build_page.py                    # regenerate index.ipynb
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF and no dataset needed.
