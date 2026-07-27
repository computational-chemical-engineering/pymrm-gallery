# J1.5 — The linear driving force

Almost every adsorption, ion-exchange and chromatography model replaces
diffusion inside the particle with a single first-order relaxation at
`k = 15 D/r²`. This page solves the sphere that approximation stands in for.

- **Structures:** `S3` (1-D BVP), `S4` (1-D transient PDE)
- **Reference:** Glueckauf (1955), Trans. Faraday Soc. 51 1540–1551,
  doi:10.1039/TF9555101540
- **Runtime:** ~3 s

## Agreement

pymrm's sphere matches the exact series to **6.6e-5** relative at `n_r = 400`,
and to 2.6e-3 at worst across `τ = 0.001` to `0.5` — the worst case being the
earliest time, where the internal profile is steepest. Spatial refinement runs
into the backward-Euler time-error floor, which the time-step study then
measures directly.

## Where 15 comes from

Assume a parabolic profile `q = a + b(r'/r)²`. Then `q* − q̄ = (2/5)b` while the
surface flux gives `dq̄/dt = 6Db/r²`; eliminating `b` leaves exactly `15 D/r²`.

## What the page adds

The coefficient the LDF *would need* at each instant,

```
k_eff(τ) = (dq̄/dτ) / (q* − q̄)
```

computed in closed form from the series, so nothing is differentiated
numerically and nothing is formed as `1 − (1 − ε)`. It is not flat:

| regime | behaviour |
|---|---|
| `τ → 0` | diverges as `τ^(−1/2)` — uptake grows as a square root, and no exponential can follow that |
| `τ = 0.022` | passes through **15**, particle 44 % loaded |
| `τ → ∞` | settles on **π² = 9.8696**, recovered to five figures |

So 15 is exactly right at one instant. Glueckauf's stated domain — near
equilibrium, linear or mildly curved isotherms — is precisely the region where
the particle stays near there. Give it a steep front and the process happens
entirely where the LDF cannot follow, under-predicting uptake by tens of
percent, which is the `K_d > 3` failure he warns about.

## Honest limits

Provenance tier 6: the reference is the analytic series, not measurement.
Glueckauf's Tables 3–5 are not transcribed — the 1955 scan drops decimal points
wholesale, rendering 0.8647 as `8647` — and nothing is lost, since the exact
series is a stronger reference than his tabulated approximation to it.

The page treats a step in surface concentration with a linear isotherm. The
coupled case, where the surface value tracks a bulk the particle is itself
depleting, is the `S8` problem of `D1.4` and is not solved here.

## Rebuilding

```bash
python build_page.py && python ../../scripts/run_pages.py
```
