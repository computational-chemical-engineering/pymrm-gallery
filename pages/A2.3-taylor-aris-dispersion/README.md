# A2.3 — Taylor–Aris dispersion

Shear plus molecular diffusion produces spreading that looks exactly like
diffusion but is thousands of times faster — and, counter-intuitively, gets
*slower* as the solute diffuses faster.

- **Structures:** `S3` (1-D steady BVP), `S6` (2-D PDE)
- **Reference:** Taylor (1953), Proc. R. Soc. Lond. A 219(1137) 186–203,
  doi:10.1098/rspa.1953.0139
- **Runtime:** ~48 s

## Agreement

Two independent routes to Taylor's Eq. 25, `k = a²u₀²/(192 D)`:

| Route | Result |
|---|---|
| homogenisation closure (1-D, no convection) | 1.0e-4 relative at `n_r = 200`, converging **O(h²)** with ratio 4.00 |
| the same closure over Pe = 2 to 398 | worst deviation **0.016 %** |
| direct 2-D transient simulation | converges **O(Δz)** to 0.7 % at `n_z = 3200`; mass conserved to 1e-13 |

Against the paper: Taylor's worked capillary run recomputes to 0.04 % in `k`
and 0.38 % in the diffusivity it implies, and that value — 5.98e-6 cm²/s — lies
inside the 4.35e-6 to 1.50e-5 range Furth & Ullmann measured independently.

## Watch the velocity convention

`u₀` is the **centreline** velocity, as in Taylor's paper. The mean is `u₀/2`,
so the equivalent form is `k = a²ū²/(48 D)`. Using the mean velocity with 192
gives an answer four times too small — the commonest way to get this wrong.

## What the page adds

Taylor bounded when his result applies — radial equilibration takes about
`a²/(3.83²D)`, so wait "much longer than" that — but could not say how much
longer, because that needs the transient solution. Running it: the growth rate
settles permanently within 10 % of Eq. 25 by **2.1 equilibration times**, 5 % by
2.6 and 3 % by 2.9. So `τ ≳ 0.2`, considerably less demanding than the phrase
suggests.

## Honest limits

Data tier 2, and the comparison is indirect: Taylor's raw concentration profiles
are not transcribed, only his worked inference and the independently measured
range it is checked against. His second inferred value (Eq. 42) is left out
because the tube radius for those runs is not stated beside the calculation, and
back-solving it implies a different tube.

The direct simulation carries first-order upwind numerical dispersion — 1.8 % at
the resolution used for the timing sweep — so the tightest accuracy band sits at
the edge of what that run resolves. The closure, which has no convective term at
all, is the precise check.

## Rebuilding

```bash
python build_page.py                    # regenerate index.ipynb
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF needed: the dataset is a committed CSV.
