# A3.15 — Graetz–Nusselt problem

The developing thermal boundary layer in a tube, solved as a 2-D PDE — and the
two amplitudes Graetz declared out of reach and calibrated out of his own
apparatus instead.

- **Structure:** `S6` (2-D PDE)
- **Reference:** Graetz (1882), *Ueber die Wärmeleitungsfähigkeit von
  Flüssigkeiten. Erste Abhandlung*, Annalen der Physik **254**(1) 79–94,
  doi:10.1002/andp.18822540106
- **Data:** tier 1 — 51 measured runs on six liquids, transcribed from the
  300 dpi page bitmaps
- **Runtime:** ~26 s

## What the paper does, and where it stops

Graetz's problem was not heat transfer; it was measuring the conductivity of
liquids with a flow calorimeter. To read his instrument he had to solve

$$2\alpha(1-r^2/R^2)\,\partial u/\partial z = a^2\,\big(u_{rr} + u_r/r\big),$$

which is the problem now named after him. He gets the eigenvalues, and then
stops:

> "Es wäre natürlich von wesentlichem Vortheil, wenn sich die $p_i$ bestimmt
> durch $R$ und $\mu_i$ ausdrücken liessen. Indess gehört dazu eine eingehende
> Untersuchung der Function $V(r,\beta,R)$, die wahrscheinlich nicht einfach
> ist … Für physikalische Zwecke lassen sich dagegen die $p_i$ sehr einfach als
> Constanten des Apparates experimentell bestimmen."

He fits $p_1$, $p_2$ **and** $k/c$ to three runs at three pressure heads.

## Agreement

Against what the paper prints:

| Check | Result |
|---|---|
| his power series for $V(\mu)=0$, vs exact rationals | 4 of 6 coefficients right to every printed digit; $\mu^8$ +0.58 %, $\mu^{10}$ ×2.65 |
| $\mu_1 = 2.7043$ | 2.4e-5 |
| the four $J_0$ zeros of his plug-flow case | worst 4.3e-5 |
| $\mu_2 = 6.50$ | **2.7 % low** — the 10-term series cannot place it |
| his calibration re-solved, $\log p_1$ | 2.3e-5 |
| his working formula vs his own $k/c$ column, 51 runs | 0.71 % mean |

Against limits he did not have, neither of which the series supplies:

| Check | Result |
|---|---|
| Lévêque constant $2/(\Gamma(4/3)\,9^{1/3}) = 1.076732$, derived | 2.9e-4, order 1.9 |
| $\mathrm{Nu}_\infty = \mu_1^2/2$ | 3.65679 |
| 2-D solve vs eigenfunction projection (no shared code) | $\mu_1$ 5.5e-6, $p_1$ 4.2e-4, $p_{1,\text{area}}$ 4.3e-4 |
| **observed order, axial / radial** | **1.94 / 2.12** |

The orders are the load-bearing pair. **The 5.5e-6 is not the accuracy of the
solve** — it is a partial cancellation of the axial and radial errors at
$n_z=600$, $n_r=128$: refining *either axis alone* makes it 7–8× worse (4.5e-5 at
600/256, 3.8e-5 at 1200/128), and the fit window and the domain length say the
same. §3b of the notebook measures all three knobs and states ~4.5e-5 as the
honest figure. No break-table row can catch a baseline that is right by accident,
which is why this is measured directly rather than left to the table.

## Two printed defects, both settled from the paper's own numbers

- **The exponent in the working formula on page 90 reads 2.4013**, where page 87
  derives 2.7043. Inverting the formula for all 51 runs decides it: 0.71 % mean
  deviation from his own $k/c$ column with 2.7043, 26.3 % with 2.4013. (2.401 is
  the brass-wall-corrected Bessel root four pages earlier — a plausible
  compositor's substitution, but the data settle it, not the guess.)
- **Two series coefficients are arithmetically wrong**, found against an
  independent exact-rational solution of his own recurrence, not against another
  printing. One caution, recorded on the page and in the data sidecar: the third
  digit-group of the $\mu^{10}$ coefficient is a **broken glyph** on the 300 ppi
  bitmap. At 16× it is a **3**, which gives the ×2.65; had it been a 9 the
  printed value would sit 2.4 % from exact and the ×2.65 claim would collapse to
  an ordinary rounding slip. The +0.58 % on $\mu^{8}$ does not depend on it.

## The result, and it is not the expected one

His calibration has three constants and three observations, so it fits exactly
and proves nothing. Substituting the computed amplitudes removes the fit — and
those three runs then **no longer agree on a single $k/c$: they spread by 15.6 %,
monotonically with flow rate**, against 0.77 % propagated from the half-division
reading error of the 1/10-degree thermometers he describes. Twenty sigma — on his
*stated* thermometry. That budget carries E₁ and E₂ only; adding the ¹/₅-degree
bath thermometer E₃ gives 17.9σ, adding a 1 % weighing/timing error 9.4σ, and
dropping the √5 averaging 8.0σ. **The honest headline is that range, 20 down to
8, and 20 is its optimistic end.** §8g prints it, and says plainly that a
break-table row which only *coarsens* an assumed source can never reveal an
*omitted* one.

**And it holds out of sample, which needs no error model at all.** Graetz sets
this test himself on page 89 — whether it matters which pressure head you use —
and asserts on page 90 that it passes. Run on every liquid he measured at more
than one head (48 of the 51 runs were never fitted to anything):

| liquid | runs | his fitted $p_i$ | computed $p_i$ |
|---|---|---|---|
| water | 13 | 3.6 % | **17.0 %** |
| copper sulfate | 11 | 4.8 % | **20.9 %** |
| alcohol | 7 | 7.2 % | **17.8 %** |

Same direction every time: the lowest head gives the highest $k/c$. Two further
checks confirm the 15.6 % is not an artefact of the method — 20 modes give the
identical three values (to 2.4e-15), and inverting through the 2-D solver's own
$\theta_{\rm cup}(x)$ with *no eigenfunction anywhere* gives 15.58 %.

Removing the fit moves every conductivity down by 16.6–24.5 % and destroys the
agreement he had with Lorberg's recomputation of Weber (his 0.0945 for water
against 0.0940). So the empirical calibration was absorbing a real systematic
effect in the apparatus that his own method could not see, and computing the
constants makes his numbers *worse*.

**The page does not diagnose the effect.** Three candidates are listed with what
each would need; two of the three move the answer the wrong way. See *What this
page cannot conclude*.

## Reuse

`solve_graetz` is a general axisymmetric convection–radial-diffusion solver:
`nu_r` for the geometry, `plug` for the velocity profile, `inv_pe2` for axial
dispersion. Copy the van Leer deferred correction **with its under-relaxation** —
without `omega < 1` the limiter switches and the iteration limit-cycles (measured
on the page: stalled at 2.0e-4 after 120 iterations, against 26 iterations to
9e-12 with omega = 0.6), silently making every number depend on the iteration
count.
The loop here raises rather than returning a non-converged field.

## Cite the source, not this page

Graetz, L., *Ueber die Wärmeleitungsfähigkeit von Flüssigkeiten. Erste
Abhandlung*, Annalen der Physik **254**(1) 79–94 (1882).
The catalogue records this as "Graetz (1883)"; the paper is signed *München,
15. Oct. 1882*. The *zweite Abhandlung*, Ann. Phys. **261** (1885), is a
different paper and was not read here.
