# B3.3 Bhatia & Perlmutter's random pore model, Part I (isothermal, kinetic control)

From Bhatia, S. K. & Perlmutter, D. D., *A random pore model for fluid–solid
reactions: I. Isothermal, kinetic control*, AIChE Journal **26**(3), 379–386
(May 1980), [doi:10.1002/aic.690260308](https://doi.org/10.1002/aic.690260308).
Read from the 8-page Wiley scan at its native 300 ppi (`pdfimages -list`:
CCITT-G4 bilevel, 300 × 300 ppi on every page); every equation and every
constant was transcribed from cropped native renders, and the text layer was
used only as a search index.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb` |
| `index.ipynb` | the executed notebook; runs in about 10 s |
| `model.py` | the random pore model, Petersen's model as reprinted, the Figure 7 comparison, the pymrm block inversions and march, and the paper's integral I |
| `test_model.py` | fast checks of `model.py` and one stored claim |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `agreement.json` | the claims, format 2, written by the last cell |
| `data/bhatia-1980-printed-scalars.csv` | the 47 numbers the paper prints outside its equations |
| `data/random-pore-reference.csv` | the page's own model evaluations, recomputed on every run |

## The scope decision, and why it is the first result

**The paper contains no table.** All eight pages were rendered at 300 ppi and
read column-block by column-block: eight figures, a NOTATION list, a
LITERATURE CITED list, nothing else tabular. The text layer grepped for
`table` returns two lines, both inside words (`suitable`, p. 385; a garbled
`surfaces`, p. 380).

The paper's only comparison with experiment is **Figure 8**, the char data of
Hashimoto et al. (1979), which is not on disk here. It is **not digitised**:
no point, curve or axis coordinate is taken from it or from any other figure.
**This page therefore does not establish that the random pore model describes
any real solid.** The paper prints its Figure 8 fit results in prose on p. 385
($S_o = 520$ m²/cm³ for both chars, $\psi = 6.9$ and $13.7$), and everything
derived from those three numbers is a prediction of the model, not a
validation of it.

## What the page establishes

- **The derivation chain, symbolically.** Fifteen `sympy` identities from
  eqs. (17)–(30) to (31)–(37), including that eq. (37)'s radical is the perfect
  square $(1+\psi\tau/2)^2$, so $S/S_o$ is $dX/d\tau$, and that Petersen's
  surface maximum sits at $\epsilon = 1/2$ for every structure, which is all
  eq. (34) says. Labelled structural.
- **The printed structural claims, as numbers.** The best-fit grain shape
  factor from the paper's own integral $I$ is 0.4905224
  <!-- claim: m_star_psi2 --> at $\psi = 2$ and 1 at $\psi = 0$ (printed
  $0.49 \le m \le 1$), and 0.6638850 <!-- claim: m_star_psi1 --> at $\psi = 1$
  (printed "very close to two-thirds"). $\tau_{1/2}$ at $\sigma = 100$ is within
  2.92 % <!-- claim: sigma100_max_rel --> of its plateau, and at
  $\sigma = 0.25$ spreads only 7.16 % <!-- claim: sigma025_psi_spread --> across
  $0 \le \psi \le 100$. The printed $X_M < 0.393$ is a truncation: $X_M$ tends to
  $1 - e^{-1/2}$ and passes 0.393 at $\psi$ = 1293
  <!-- claim: XM_exceeds_printed_psi -->.
- **Both special-case reductions the abstract claims.** Bhatia's eq. (40) is
  the published `B3.1`'s Yagi–Kunii eq. (6) at $\omega = 0$, to round-off, and
  `B3.1`'s committed reaction-control gap (loaded with `load_claims`) is
  reproduced on its own grid. The approach to Petersen (1957) is measured:
  under the convention p. 384 prints for curve B ($\epsilon_o$ and $L_o$ fixed),
  the conversion curves differ by at most 0.0218
  <!-- claim: fig7_max_dX_upto_75 --> up to $X = 0.75$ and first part by 0.02 at
  $X = 0.7423$ <!-- claim: fig7_X_at_two_points -->, the paper's "about 75 %".
  The other two pairings of the over-determined triple part at 0.7599
  <!-- claim: fig7_X_part_eps_S --> and 0.7687 <!-- claim: fig7_X_part_S_L -->.
  The two overlap laws, $1-e^{-V_E}$ and $3u^2-2u^3$, differ at order
  $V_E^{3/2}$, so no parameter choice makes the models identical.
- **One thing the paper's own numbers say that it does not.** Figure 7's curve C
  asks for a Petersen surface 4.596 <!-- claim: curveC_So_over_Smax --> times
  the largest any structure of its length can have: the triple is infeasible,
  which is why a parameter had to be freed. The curve-D length the paper never
  prints is 24.672 <!-- claim: curveD_L_over_LAB --> times curves A and B's.
- **The rate law marched.** Crank–Nicolson on eq. (37), a pymrm `newton` solve
  per step, reproduces eq. (32) to 9.8e-8 <!-- claim: march_cn_vs_eq32_max -->
  at 800 steps, observed order within 2 % of 2 <!-- claim: march_observed_order -->
  from successive differences.

## Three printed defects, reported and not repaired

1. **"This match at $m = 1$"** [*sic*, p. 383], where the preceding sentence,
   Figure 6 and the dashed line's label all say $m = 2/3$ at $\psi = 1$. At
   $\psi = 1$, $m = 1$ misfits 5.336 <!-- claim: m1_penalty_at_psi1 --> times
   worse than the optimum. Which repair was intended cannot be decided.
2. **The optimal-porosity worked example.** The printed $\epsilon_o = 0.1$ is not
   a root of the paper's own eq. (35); the root is 0.0852877
   <!-- claim: eq35_root_eps_o -->. The optimum is flat, so 0.1 costs only
   2.76e-5 <!-- claim: eq35_S_deficit_at_printed --> of the maximum surface.
3. **Park and Levenspiel (1976)** in the p. 380 text against **(1975)** in the
   LITERATURE CITED.

## The checks that can fail, and the one that cannot

Every claim has a break-table row that injects one defect through an argument
of `model.py` or a printed parameter, and fails the claim's own check; the
only exception is the count of symbolic identities, which is structural. Two
identities are narrower than they look, and the page says so: Petersen's two
routes for surface and conversion share algebra once $G = 3/(2u_o)$, so they
test the root of eq. (7) and the transcription of eq. (4), not eq. (6); and
$X(s = 2G/3 - 1) = 1$ holds for every root of eq. (7), so its break row
misreads eq. (7) itself. The `B3.1` identity at $\omega = 0$ is blind to all
three $\omega$ coefficients and to the value of $\gamma$; the two things it can
see (the exponent on $x$ and the ratio of the two $\gamma$ coefficients) each
have a row. Three candidate rows cannot exist, and the page prints why.

## Traps for anyone reusing this

- Eq. (7) has three real roots; only $G = 3/(2u_o)$ is consistent with eqs. (3),
  (4) and (6).
- $\psi$ carries $(1-\epsilon_o)$; eq. (32) carries $\psi\tau/4$, not $\psi\tau/2$.
- Eq. (3) is not monotone in $r_p$: two roots for any attainable $S$, none above
  $S_{\max}(L)$.
- Eq. (32) inverted from a start deep in the conversion loses its exponential
  against 1 and gives `newton` an exactly singular Jacobian; the page solves the
  log form.
- The closed form of $I$ loses digits to cancellation when the misfit is small.

## Scope, honestly

Kinetic control throughout, the paper's own title. Part II (non-isothermal
operation and diffusion control) is not on disk and nothing here depends on it.
Petersen (1957) and Avrami (1940) are read only as this paper reprints them.
`B3.2` is the neighbouring page where intergrain diffusion matters; `B3.1` is
the shrinking core this model collapses onto.

## Corrections

- 2026-11-02: said a Gauss-Legendre quadrature reproduces the closed form of
  $I$ "at observed second order in the node count", and called the 1e-12 left at
  400 nodes the quadrature error. Gauss-Legendre converges exponentially here and
  is below 1e-13 <!-- claim: I_gl400_vs_40digit_rel --> against a 40-digit
  quadrature at 400 nodes. The remaining gap to the closed form is the closed
  form's own cancellation error, below 1e-10
  <!-- claim: I_closed_vs_40digit_rel -->.
- 2026-11-02: said the best grain shape factor misfits by more than 8 % at
  $\psi = 2$. It misfits by 7.99 % <!-- claim: nrms_best_psi2 --> there and
  reaches 8 % at $\psi$ = 2.001 <!-- claim: nrms_8pct_psi -->.
- 2026-11-02: said the two routes to the best-fit $m$ agree to 1.4e-8. That gap
  was the resolution of a minimiser on the closed form of $I$, whose cancellation
  noise near the minimum limits it to about 1e-7. The root of $dI/dm$ agrees with
  a 40-digit root to below 1e-11 <!-- claim: m_star_two_routes -->.
- 2026-11-02: listed the printed $X_M < 0.393$ among the printed claims the
  page recovers. The supremum of eq. (33) is $1 - e^{-1/2}$, and $X_M$ passes
  0.393 at $\psi$ = 1293 <!-- claim: XM_exceeds_printed_psi -->; the printed
  bound is a truncation.
