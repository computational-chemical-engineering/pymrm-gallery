# A3.12 — Yagi–Kunii effective thermal conductivity

A packed bed conducts heat by several routes at once. Yagi & Kunii's Equation (15)
adds them — solid conduction, the thin gas film at the contact points, radiation
surface-to-surface, radiation void-to-void — and Equation (6) puts Ranz's lateral
mixing in parallel once gas flows. It is still the skeleton of every stagnant
conductivity correlation in use.

- **Structures:** `S6` (cylindrical `nu=1` radial conduction, for the apparatus)
- **Reference:** Yagi & Kunii (1957), *AIChE J* **3**(3) 373–381,
  doi:10.1002/aic.690030317
- **Runtime:** ~12 s

## Three kinds of result, kept apart

The paper contains theory, *other people's* data re-analysed, and the authors' own
measurements. The page labels each:

| | where | label |
|---|---|---|
| Equations (1)–(19) | pages 374–376 | reproduced, cross-checked against their own printed limits |
| Table 1 (twenty flow studies) and Table 3 (twenty stagnant rows) | pages 373, 377 | third-party measurement + the authors' reduction of it |
| iron spheres, porcelain, cement clinker, firebrick, Raschig rings, 0–1000 °C | **Figures 13–17, the plotted points** | **scoped out** — figures, and no maintainer is available for a digitisation review |
| the smooth **theoretical** lines in those same figures | Figures 13–17 | *read* — the authors' computed output, 17 curve-end values, used as a reproduction check |
| the endpoints of Figure 12's five temperature profiles | page 378 | *read* — 10 values off a labelled ordinate, to fix the Δt the apparatus ran |

## The finding

**One of Equation (15)'s five ingredients is untestable on the only table of numbers
in the paper — and it is the radiation term the paper's central claim rests on.**
Equation (15) needs `phi`, the effective contact-film thickness, which the paper
supplies only as Figures 9 and 10. With `phi` free per row, the solid-to-solid
coefficient `h_rs` enters only through the sum `1/phi + Dp·h_rs/kg` and is
algebraically degenerate with it:

| term deleted | rows `phi` can no longer match | worst shortfall | the `phi` it demands |
|---|---|---|---|
| solid–solid radiation, Eq. (7) | **0 of 14** | 0 % | 0.0075–0.054 — inside the baseline's own span |
| void–void radiation, Eq. (8) | 2 of 14 | 14.2 % | 0.0022–0.053 |
| solid conduction, `gamma(kg/ks)` | 0 of 14 | 0 % | **0.050–0.42** — 12 rows above every printed value, 2 off Figure 9's axis |
| *control:* numerator halved | 6 of 14 | 46.7 % | 0.0012–0.020 |

Deleting Eq. (7) entirely and re-solving for `phi` reproduces every printed
calculated value to **2.2e-16**, with a `phi` that never leaves the band the baseline
inversion already occupies. The control shows the test is not blind. **`gamma`, by
contrast, is identifiable** — it only looked degenerate because the search was allowed
to run to `phi = 1`, which is a "film" as thick as the particle. And the paper's
central conclusion — radiation dominant above 400 °C — concerns a regime **no row of
Table 3 reaches**; its hottest is 300 °C.

What breaks the degeneracy is temperature: `phi` is a geometric ratio, `h_rs` scales
as `T³`. That experiment exists — it is Figures 13–17 — and its measured points are
the part of the paper this page cannot read.

## Agreement

| check | result |
|---|---|
| Eq. (16) = Eq. (15) with `h_r = 0`, `gamma = 1` | 4.4e-16 over 20 000 draws |
| Eq. (14) = the `kg → 0` limit of Eq. (13) | residual vanishes as `O(kg)`, order **1.0000** |
| Eq. (13) rebuilt from the Eqs. (9)–(12) network | 8.9e-16 |
| Eq. (19)'s constant 0.743, from the closed form with `2*pi*l` | 0.743353 — the printed value to three digits |
| the same from a pymrm `nu=1` annulus solve carrying `ln` exactly | 0.744188, **+0.16 %** — the paper's own rounding of ln 10 to 2.3 |
| Table 4's `ks/kg = 4.48` from its own `ks` and the printed `kg(400 °C)` | −0.13 % |
| **Eq. (15) vs the ends of the paper's own printed theoretical curves, Figs. 13–17** | **mean 2.4 %, worst 4.2 % over 14 readings** |
| Table 3, 17 measured rows vs the paper's calculated column | mean \|dev\| **15.3 %**, bias −3.5 %, worst 52.0 % |
| Table 1's intercept vs Table 3's `Exp.`, five shared references | 8.8 % worst — but **exactly 0 % on three**, see below |

Three stated conclusions reproduced from printed parameters: radiation reaches half
of `ke0` at **461 °C** (iron spheres) and **265 °C** (insulating firebrick; 233 °C if
Table 4's printed `ks/kg` is held instead of its `ks`) against the paper's "higher
than 400 °C"; 0.18 mm cement clinker varies **9.5 %** over 200–1000 °C against
**154 %** for the 5 mm grade ("nearly constant for particles smaller than 0.35 mm");
randomly packed Raschig rings run **12.9–24.4 %** above regularly packed with no
fitted parameter ("about 20 % greater").

**The one comparison that is not what it looks like.** Table 3 daggers all five of the
shared references as "Extrapolated values, `N_ReM → 0`", so its `Exp.` column and
Table 1's intercept are *the same quantity* — the authors' own zero-flow extrapolation
of the same third-party flow experiments. Three of the five agree to 0.000000 %
because they are one number printed twice. The page keeps the comparison as a
transcription and internal-consistency check and labels it as one.

## Two printed defects, both proved from the paper's own results

- **Equation (17) is a factor of two low.** It sets `pi*l` where the steady
  cylindrical solution requires `2*pi*l`. The paper's own Eq. (19) constant, 0.743,
  proves it: the printed form gives **2.00095 ×** that. Anyone re-deriving the
  reduction from Eq. (17) would report every `ke0` twice too large. Eq. (19) itself
  is right, so no measurement in the paper is affected.
- **`gamma = 1.2 mm/9.0 mm` is printed as 0.0134**, where the division gives 0.1333.
  Table 4 carries the ratio, so the paper's own calculations are unaffected; using
  the decimal shifts `ke0/kg` by +9.2 % at 400 °C.

## What pymrm adds

The correlation is algebra and needs no solver. `nu=1` earns its place in the
*apparatus*:

- **The reduction formula is biased against its own model, by about +10 %.** For
  steady radial conduction the conductivity Eq. (19) returns is exactly the average of
  `ke0` uniform in **temperature**, while the abscissa is the **volume**-weighted mean
  temperature, which sits 0.703 of the way from the hot face to the cold one. At the
  five temperature drops Figure 12 actually shows — 195 to 260 °C — a bed obeying
  Eq. (15) exactly plots **+9.8 % to +11.5 %** above the Eq. (15) curve, and that
  offset is nearly constant from a mean bed temperature of 110 °C to one of 690 °C.
  Swept to Δt = 800 °C it reaches +21.7 %, but that is about three times the largest
  drop the paper ran and the page labels it an extrapolation.
- **"Nearly straight" is a measurement, and it cuts one way harder than the other.**
  A constant-conductivity annulus at this bore ratio departs from a straight line by
  13.7 % of the temperature drop; with `ke0(T)` from Eq. (15) it still departs by
  7.4–8.3 % at the drops Figure 12 shows — a factor of **1.7 to 1.9**, not the 4.5 a
  sweep to Δt = 600 °C would give. So the remark is strong evidence *against* a
  constant conductivity and only partial evidence *for* Eq. (15).

## Caveats

- **The measurement comparison is weak by construction** — see the finding above.
  Treat the radiation terms as *supported by evidence this page has not examined*,
  not as validated here.
- **The emissivity `p` is never printed in the paper.** `p = 1` is assumed. The
  curve-end comparison is the only quantitative evidence about it: consistent with
  `p = 1` (+2.4 % mean) and equally with `p = 0.9` (−1.3 %), inconsistent with
  `p ≲ 0.8` (−4.9 %). At `p = 1` both radiation coefficients collapse to the same
  `0.1952(T/100)³`, so every other metric on the page is blind to the emissivity
  structure of Eqs. (7) and (8); three break rows exist for exactly that.
- **Table 3's calculated column cannot be reproduced from printed material.** At
  `eps = 0.40` it implies `phi` spanning a factor of **7.6** — or 2.3 without the one
  pathological row — while Table 4 assigns that voidage the single value 0.034.
- **For insulating firebrick, Table 4's printed `ks/kg = 4.48` is the input the
  authors used, not its `ks = 0.20`.** Holding `ks` misses their own Figure 16 curve
  by 50–61 % at low temperature; holding the ratio matches to a few per cent.
- **The variable-`k` annulus solve is n-limited.** Converged and stable to four digits
  over `n = 100–800`; the fixed point does not reach tolerance at `n = 1200` without
  relaxation.
- **Ranz (1952) is cited, not consulted.** It is not on disk; nothing about it is
  written from memory.

## What the break table cannot detect

77 metrics, 25 injections, 77 of 77 moved, and two assertions tying every one of them
to the same `recompute()`. But every row perturbs an input and watches a number move,
so **none of them can catch a baseline that is wrong rather than insensitive** — which
is exactly how the three 0.000000 % agreements above survived an earlier draft. The
defence is the one number computed a second, independent way: Eq. (15) against the
ends of the paper's own printed theoretical curves.

## Reuse

`phi` is the parameter that will bite you. If you need a number, Table 4's printed
pairs are the defensible source: 0.034 at `eps = 0.40`, 0.040 at 0.43, 0.050 at 0.50,
0.060 at 0.54, which an exponential fit reproduces to under 3 %. If you re-use
Figures 13–17 as data, subtract about ten per cent for the reduction bias.
