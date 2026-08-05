# A2.3 — Taylor–Aris dispersion

Shear plus molecular diffusion produces spreading that looks exactly like
diffusion but is thousands of times faster — and, counter-intuitively, gets
*slower* as the solute diffuses faster.

- **Structures:** `S3` (1-D steady BVP), `S6` (2-D PDE)
- **Reference:** Taylor (1953), Proc. R. Soc. Lond. A 219(1137) 186–203,
  doi:10.1098/rspa.1953.0139
- **Runtime:** ~110 s

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
longer, because that needs the transient solution. Running it, on the finest
mesh: the growth rate settles permanently within 10 % of Eq. 25 by **2.1
equilibration times**, 5 % by 2.7, 3 % by 3.1 and 2 % by 3.4. So `τ ≳ 0.2`,
considerably less demanding than the phrase suggests.

### That result was not refined until 2026-08-05, and it moved

Every other result on this page is refined — the closure over `n_r` = 25…800,
the direct simulation over `n_z` = 400…3200. The settling study, which is the
page's *only original claim*, ran once at `n_z` = 1600. It was not converged:

| `n_z` | dz | numerical-dispersion floor | τ(±10 %) | τ(±5 %) | τ(±3 %) | τ(±2 %) |
|---|---|---|---|---|---|---|
| 800 | 0.100 | +4.46 % | 0.1275 | 0.1550 | never | never |
| **1600** (published until 2026-08-05) | 0.050 | +2.23 % | 0.1375 | 0.1750 | **0.1950** | **never** |
| 3200 | 0.025 | +1.11 % | 0.1450 | 0.1850 | **0.2125** | **0.2325** |

Two things were wrong, and both ran the same way. **A coarse mesh flatters a
settling criterion**: first-order upwind lifts the whole growth-rate curve above
Taylor's value by a constant proportional to Δz, so a band is entered *earlier*
than it should be and every threshold is biased low — the ±3 % threshold by 9 %.
And **a band narrower than that bias can never be entered at all**, however long
the run, which is how "the ±2 % band never settles" came to be published as a
result about the physics when it was a statement about the mesh.

The floor halves exactly with dz (observed order 1.000, Richardson limit
−0.001 %), so it is entirely discretisation and the scheme converges onto
Taylor's constant. Refining `dt` (5e-4 → 2.5e-4) or `n_r` (24 → 48) instead
moves it by under 0.02 % — the study now refines the knob that actually controls
the answer. The thresholds were still rising on the last doubling, so read them
as **lower bounds**.

Also fixed: the Reuse section said `τ ≳ 1` where the result section said
`τ ≳ 0.2` — an internal contradiction, by a factor of five, in the advice a
reader acts on. The result section was right.

## Honest limits

Data tier 2, and the comparison is indirect: Taylor's raw concentration profiles
are not transcribed, only his worked inference and the independently measured
range it is checked against. His second inferred value (Eq. 42) is left out
because the tube radius for those runs is not stated beside the calculation, and
back-solving it implies a different tube.

The direct simulation carries first-order upwind numerical dispersion — 1.11 %
at the finest resolution used for the timing sweep, halving with every mesh
doubling — so the tightest accuracy band sits at the edge of what that run
resolves, and the ±1 % band cannot be reported at all. The closure, which has no
convective term, is the precise check.

Every published metric has a row in the notebook's defect-injection table
(section 5). `mass_conservation_err` is **structural** — every boundary is
no-flux and the scheme is conservative, so it stays at 1e-13 through a wrong
geometry and an unresolved radial profile alike. It is a regression guard on the
assembly, not evidence about the physics.

## Rebuilding

```bash
python build_page.py                    # regenerate index.ipynb
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF needed: the dataset is a committed CSV.
