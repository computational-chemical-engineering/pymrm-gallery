# I1.3 — Heck, Wei & Katzer 1976: what the 1-D monolith model is adequate for

The paper's two claims, tested rather than illustrated — Table 1's asymptotic
Nusselt numbers derived from scratch, the Nusselt singularity root-found and shown
to sit 8.7 % downstream of the cause the text gives it, and the 1-D/2-D adequacy
claim re-measured with the thresholds root-found instead of sampled.

- **Structure:** `S6` (2-D PDE)
- **Reference:** Heck, R. H., Wei, J. & Katzer, J. R., *Mathematical Modeling of
  Monolithic Catalysts*, AIChE Journal **22**(3) 477–484 (1976),
  doi:10.1002/aic.690220310
- **Data:** tier 6 — the paper's own printed values. One printed table (Table 1)
  and every stated number in running text, read off the 300 ppi page bitmaps
- **Runtime:** ~3.6 min

## The two claims

> "A two-dimensional model is shown to predict unusual behavior of the Nusselt
> number in the presence of rapid reaction. However, a simpler one-dimensional
> model is adequate for predicting monolith behavior."

Neither needs a figure. Table 1 anchors the first; the second is a comparison the
page can run itself.

## Scope: figures are out

Digitisation needs a maintainer review that is not available, so Figs. 1–13 are
used only for **printed characters inside their frames**. The notebook enumerates
every such row and **counts them from the CSV itself** rather than restating a
number here, because that list is what the scope promise is audited against: of the
31 transcribed rows, 15 name a figure, **11 were read inside a figure's frame**
(6 annotation boxes, 3 abscissa tick labels, 1 curve label, 1 case label) and 9 of
those 11 appear nowhere else in the paper. The three abscissa rows are Fig. 2
(0.20), Figs. 3–5 (0.25) and Fig. 11 (0.50) — the first and third because the
paper's two light-off windows are quoted for those two figures and they are not the
same length, the second because it is the other length the paper fixes for its 1-D
model; together they set the swept $L^{*}$ range. Each row is labelled with the
figure it sits in. No curve is traced. What that costs is
on the page: no profile shape, no light-off *position* against the paper's own,
and none of the conversion curves of Figs. 4–6 and 12 is checked. Every comparison
is against something the paper says in words, and there turned out to be a great
many of those.

## Table 1, derived rather than quoted

| cell | printed | this page | route |
|---|---|---|---|
| circle, $Nu_H$ | 4.364 | **48/11** = 4.363636 | closed form, derived |
| triangle, $Nu_H$ | 3.111 | **28/9** = 3.111111 | closed form, derived |
| square, $Nu_H$ | 3.608 | 3.607950 | pymrm FV **and** Fourier–Galerkin |
| circle, $Nu_T$ | 3.657 | 3.656793 | pymrm FV (= `A3.15`'s $\mu_1^2/2$ to 7 digits) |
| square, $Nu_T$ | 2.976 | **2.977523** | pymrm FV **and** Fourier–Galerkin |

Four agree to every printed digit. The fifth does not: `2.976` implies
[2.9755, 2.9765], which excludes 2.977523, so it is not rounding. **The page does
not attribute it** — Table 1 is credited to Shah and London (1971), which was not
consulted. (The paper's reference list prints that source as "**Shar**, R. K., and
A. L. London" where Table 1's own footnote says "Shah and London" — reported, not
repaired.)

Three cells are not reproduced, and one of them *cannot* be: the sinusoidal
geometry is printed as a glyph with no aspect ratio anywhere, so the paper does not
determine the problem.

## Claim 1: the singularity is real, the stated cause is not

Both wall-flux functions the paper prints give the same picture, and every landmark
is root-found:

| | Fig. 8's $-9.46+(46.0X)^8$ | Fig. 9's $-3+(40X)^8$ |
|---|---|---|
| flux zero, where $Nu$ passes through **zero** | 0.028789 | 0.028680 |
| $Nu$ pole, zero of $T_w-\bar T_G$ | 0.031300 | 0.031183 |
| pole downstream of flux zero | **+8.72 %** | **+8.73 %** |

The text says "as the wall flux passes through zero, infinite negative and positive
values for the Nusselt number occur". At the flux zero eq. (19)'s *numerator*
vanishes, so $Nu = 0$ there; the pole is the *denominator*'s zero and lies further
on, because the bulk temperature integrates the flux history and keeps falling
after the flux turns positive.

**This is a wrong summary, not a missed effect.** The paper's very next paragraph
gives the correct mechanism explicitly, in terms of the two signs. Its own
$X^{*} = 0.032$ and $0.033$ statements bracket the pole in (0.032, 0.033), and this
solve puts it 4.1 % below the centre of that bracket — and **not a grid effect on
any grid the paper could plausibly have run**: at its own ten radial points the pole
moves by 0.33 %, and on its *full* stated discretisation (ten radial points,
$\Delta X^{*} = 0.00025$, TVD off so the axial scheme is first order like a
1955-vintage marching code) it is 0.031241, 2.4 % short of 0.032. *How* short is
scheme-dependent — 0.031241 is this page's own radial operator on the paper's stated
counts, not the paper's 1955 code — so the page states the **margin** instead: the
two paper-count solves span 2.2 to 2.4 % short and sit within 0.3 % of converged,
while reaching 0.032 at all takes +2.6 %, about 8× the discretisation error the
paper's own stated grid produces here.

## Three meanings of one symbol, and a check value that measures the grid

$X^{*}$ is defined once, in the Notation, as $4x/(D\,Re\,Pr)$. It has to mean three
different things for the paper to be self-consistent, and each is settled here
against a computation:

| where | required scaling | how it was settled |
|---|---|---|
| Notation, eq. (9) | $4x/(D\,Re\,Pr)$ | the 1-D energy balance with $Nu = hD/k_G$ |
| eqs. (15)–(16), every figure abscissa | $x/(D\,Re\,Pr)$ | the paper's own printed $\bar T^{*}(0.25) = 0.0224$ |
| the two Grigull–Tratz correlations | $1000\,x/(D\,Re\,Pr)$ | reproduces this page's Graetz solve to 1.3 % |

That printed 0.0224 does more than settle the scaling. Converged, the answer is
**0.021147**; on the **ten radial points p. 480 says were used** it is 0.02158
(uniform) to 0.02245 (wall-clustered). So the number is not a check on the paper's
2-D model, it is a measurement of that model's radial discretisation error, and
**5.9 %** is its size — the same order as the 1-D/2-D difference the model is being
used to adjudicate. That is why every adequacy statement here is a *difference
between two models on one grid*, never an absolute.

## Claim 2: the paper's thresholds are sampled; these are root-found

The two quoted windows — 1-D "below 304 °C … between 304° and 343 °C" (p. 479,
about **Fig. 2**), 2-D "$T^o_G <$ 293 °C … 293 $<T^o_G<$ 343 °C" (p. 481, about
**Fig. 11**) — are read off the four inlet temperatures each figure happened to
run. They are also quoted for channels of *different length*: Fig. 2's abscissa is
labelled to 0.20 in $X/(D\,Re\,Pr)$ and Fig. 11's to 0.50. Run at each figure's own
length, this page's own models reproduce **both** statements and the apparent gap
between them:

| the paper says | at its figure's own length | this page |
|---|---|---|
| 304 °C does not light off (Fig. 2, ends at 0.20) | $G^{*}(304) = 0.2405$ | beyond the figure ✓ |
| 293 °C does not light off (Fig. 11, ends at 0.50) | $G^{*}(293) = 0.5767$ | beyond the figure ✓ |
| the two windows differ by −11 °C | 306.44 vs 295.61 °C | **−10.83 °C**, reproduced to 0.17 |

The gap is **decomposed rather than asserted**, because it is not purely length:
taking the *same* 1-D model out to Fig. 11's length gives **−13.42 °C** from the
length alone and **+2.59 °C** from the model difference at the common $L^{*}=0.50$.
Length dominates and is what fixes the sign; the model difference — +4.73 °C at
$L^{*}=0.20$ down to +2.59 at 0.50 — is what survives.

**A retraction.** An earlier draft of this page read the two windows as thresholds
on one channel and called the 11 °C a contradiction with the mechanism the authors
derive on p. 482 ("lower Nusselt number means earlier light-off"). **That was
wrong and it is withdrawn** — in the notebook, here, in `meta.yaml` and in the
`models.yaml` entry. It is a channel-length artefact: the paper's two numbers are
correct about their own figures, consistent with each other, and consistent with
its own mechanism. What survives is that the windows are *sampled*, and that on a
*common* length the 2-D sits above the 1-D.

Root-found on a common channel length, with identical parameters:

| $L^{*}$ | 1-D, $Nu=3.608$ | 1-D, Grigull–Tratz $Nu(X^{*})$ | 2-D | 2-D − 1-D |
|---|---|---|---|---|
| 0.20 | 306.44 | 308.89 | 311.17 | **+4.73** |
| 0.25 | 303.47 | 305.65 | 307.59 | **+4.12** |
| 0.30 | 300.89 | 302.87 | 304.56 | **+3.67** |

The 2-D model needs a *higher* inlet temperature at every length — the direction the
paper's own argument requires. The magnitude is length-dependent (about 2.8 °C per
0.05 of $L^{*}$), which is exactly why the paper's two windows cannot be compared
with each other.

**About half of that gap is channel shape, not dimensionality.** Split like-for-like
on the *same* threshold and the *same* length, because the paper's 1-D runs use the
*square*-channel Nusselt number of its own Table 1 against a *round* 2-D tube:

| at $L^{*} = 0.25$ | threshold | |
|---|---|---|
| 1-D, $Nu = 3.608$ (square — what the paper's 1-D uses) | 303.47 °C | |
| 1-D, $Nu = 4.364$ (round tube — shape matched to the 2-D) | 305.54 °C | shape **+2.07** |
| 2-D (round tube, entrance region resolved) | 307.59 °C | entrance + dimensionality **+2.05** |

An earlier draft measured the shape effect on the *upper* threshold (329.9 → 335.2,
+5.29 °C) and called it "as large as the whole gap" measured on the *lower* one.
That cross-threshold comparison is withdrawn; the like-for-like split above replaces
it, and the upper-threshold pair is still reported as a sensitivity of that
threshold.

**There is no 2-D inlet threshold to report.** With an infinite local Sherwood
number at $X^{*} = 0$, light-off always occurs at some $X^{*} > 0$, so the "strict"
inlet threshold a march returns is just the inlet temperature at which light-off
lands inside the first step: 378.9, 386.9, 395.0, 403.3 °C at $dG$ = 4e-4, 2e-4,
1e-4, 5e-5 — about +8.2 °C per halving, diverging. The metric carries its step size
in its name (`T_up_2d_first_step_dG1em4`) and has a $dG$ break row. An earlier draft
reported it as a converged `T_up_2d_strict`.

### "Less than a 2 °C difference in inlet gas temperature"

Measured as an equivalent inlet shift, by matching light-off positions:

| 2-D inlet | Grigull–Tratz 1-D (the Fig. 12 model, which the claim is about) | constant-$Nu$ 1-D (Figs. 2–6) |
|---|---|---|
| 310 °C | +2.16 | +4.53 |
| 316 °C | **+2.84** | **+5.71** |
| 330 °C | +5.13 | +10.11 |
| 340 °C | +6.78 | +15.46 |

The claim is the right size for the *procedure* it was made about — within about
1.5× at the inlet temperature of the figure that motivated it — and 2–5× too small
for the constant-$Nu$ model that produced the paper's parametric study. **Adequacy
belongs to the Fig. 12 procedure, not to one-dimensionality.**

## Two groups had to be reconstructed, and that is why the page states differences

The paper prints neither $M_G C_{PG}$ nor $\mathcal{D}$, so $\Delta T_{AD}$ = 425.2 K
and $\beta$ = 3.327e-3 come from ordinary gas properties. **Nothing is fitted to any
result of the paper** — and **neither reconstructed group is checked by anything the
paper prints**, which the page measures rather than glosses:

* The rate maximum root-finds at $C_w$ = 0.00261 against the stated 0.003 (one
  significant figure). But $C_w$ at the maximum barely depends on
  $\Delta T_{AD}$: root-found, it still rounds to 0.003 for **every**
  $\Delta T_{AD}$ from 0 up to 1474 K (minimum 0.002509, against a rounding edge of
  0.0025). What that check really verifies is the *transcription* of $k_a^{o}$ and
  $E_a/R$ — at $\Delta T_{AD} = 0$ the maximum is exactly $1/k_a$ = 0.00299 — and
  both were printed, not reconstructed. $\beta$ does not enter it at all.
* The implied fully lit wall at 741 °C has **no printed referent**: the paper
  states no adiabatic flame temperature in words anywhere in pp. 477–483 (checked
  page by page on the 300 ppi bitmaps; the ~750 °C plateau is a plotted curve, and
  figures are scoped out). It is a consequence of the reconstruction, not a test.

An earlier draft called these "two printed facts that check the reconstruction".
That is withdrawn.

What holds the reconstruction to account instead is §7f, which measures its cost:
±10 % in $\Delta T_{AD}$ and ×1.5 in $\beta$ move the
*absolute* upper threshold by **22.3 °C** and the 1-D/2-D *equivalent shift* by
**4.55 °C** — five times less sensitive in degrees, but not insensitive. What
survives every variation is the two statements this page makes: the constant-$Nu$
1-D model exceeds the 2 °C claim under all of them (smallest +3.96 °C) and the
Grigull–Tratz model is at or above it under all of them (smallest +1.84 °C).
**So no absolute threshold here should be read as the paper's number** — including
the 303.5 °C that lands within a degree of the printed 304. The notebook says so in
the cell that prints that table.

## A trap worth copying out

The 1-D light-off position is an integral along the cold branch. Parametrised by
$C_G$ it has a square-root endpoint, because $dC_G/dC_w$ vanishes at the fold, and
scipy's adaptive `quad` comes out up to **0.97 %** wrong on it. Parametrised by
$C_w$ the integrand *vanishes* at the fold and 100 Gauss points converge it. Had the
$C_G$ form been the reference, a 1 % error would have sat under every threshold on
the page and **no break row would have shown it** — which is why the two routes here
differ in parametrisation as well as in discretisation.

## Reuse

`march_2d` is the piece worth taking. The interior is *linear* in the wall flux, so
each implicit step needs exactly two linear solves (flux 0 and flux 1) and then a
**scalar** root-find for the wall value. A 2-D nonlinear reacting-wall problem
collapses to the same one-dimensional algebra as Fig. 1's mass-transfer line, the
fold structure becomes explicit, and light-off can be root-found rather than
stumbled over.

`Graetz` takes an **array-valued Neumann `d`**, so a wall flux that varies with the
axial coordinate costs nothing and the operator is still assembled once — that is
what makes the paper's Figs. 8–9 reproducible. Copy the van Leer deferred correction
*with* its `omega = 0.6`, as on `A3.15`.

## Cite the source, not this page

Heck, R. H., Wei, J. & Katzer, J. R., *Mathematical Modeling of Monolithic
Catalysts*, AIChE Journal **22**(3) 477–484 (1976), doi:10.1002/aic.690220310.

Shah and London (1971), Voltz et al. (1973) and Grigull and Tratz (1965) were
**not consulted**; all three are used exactly as Heck et al. reprint them.
