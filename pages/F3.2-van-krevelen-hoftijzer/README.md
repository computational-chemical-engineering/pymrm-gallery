# F3.2 — Van Krevelen–Hoftijzer: the chart replaced by a solve

The enhancement factor for a second-order gas–liquid reaction is an implicit
relation. Textbooks hand it over as a chart. Solving it instead prices the
approximation everywhere on that chart, and lets the assumption behind it be
measured rather than repeated.

- **Structure:** `S3` (1-D steady BVP)
- **Origin, cited but NOT consulted:** Van Krevelen, D. W. & Hoftijzer, P. J.
  (1948), *Kinetics of gas-liquid reactions. Part I. General theory*, Recueil
  des Travaux Chimiques des Pays-Bas **67**(7), 563–586,
  doi:10.1002/recl.19480670708. Not on disk. Both books on disk print only
  "67, 563 (1948)": the issue number, the end page and the DOI are a
  **Crossref record lookup** on that DOI, not a reading of the article.
- **Read, and the source of every equation:** Froment, De Wilde & Bischoff
  (2011), *Chemical Reactor Analysis and Design*, 3rd edn, Ch. 6 §§6.3.1–6.3.5
  (book pp. 326–341) — eqs. (6.3.2-11) and (6.3.5-1).
- **Read, as an independent check on the constants:** Levenspiel (1999),
  *Chemical Reaction Engineering*, 3rd edn, Ch. 23, Fig. 23.4 (book p. 530).
- **Runtime:** ~3.6 min — 218.6 s and 218.8 s on two consecutive executions on
  a 16-core box; `meta.yaml` declares 230 s

## The baseline, first, because everything here is an error against it

The reference film solution is **not a single solve**. It is a Richardson
extrapolation of the same graded film on n_x = 400 and 800, which the observed
order 2.0000 licenses. A single graded solve at n_x = 800 is grid-limited
exactly where this page's headline sits — its error, 1.9e−6 relative, is
*smaller* than the page's own pymrm-versus-collocation agreement, so a second
solver family cannot see it and only a convergence study can. It is the
"extrapolation switched off" row of the break table, and what that row moves is
**counted by the notebook** rather than asserted: **27 of the 42 reported
metrics** move by more than the coverage tolerance, **19 of them in their
fourth decimal**. An earlier version of this page said "nine reported numbers",
counted by hand; it was an undercount.

| the baseline, checked | |
|---|---|
| reference vs. scipy collocation, 20 chart points | 4.37e−8 |
| reference vs. a finer extrapolation (800, 1600) | 3.84e−8 |
| single n_x = 800 graded solve vs. collocation | 5.99e−6 — **137× worse** |
| first integral q·β(0) = 1 + q − F_A, never told to the discretisation | 6.4e−10 |

## Results

**The approximation, priced over the chart's own printed domain.** The error
surface has one interior extremum, root-found in both directions:

| | value |
|---|---|
| worst error | **−2.7955 %** at γ = 4.7279, q = 1.5217 (E_i = 2.5217) |
| worst over the ten *printed* curve labels only | −2.7294 % |
| the book's printed claim | "valid to within 10 percent" — margin **3.5772×** |
| sign | Van Krevelen–Hoftijzer **under**-estimates everywhere — and provably so |

The sign is a theorem, not a 612-point sample: β′ = (α′ + F_A)/q ≥ 0, so
β(ξ) ≥ β(0) everywhere; writing T(F) for the pseudo-first-order enhancement
with β frozen at β(0) = 1 − (F−1)/q, (6.3.5-1) says F_VKH = T(F_VKH) while the
true film sees more B than that, so F_film ≥ T(F_film); T decreases, so
F − T(F) increases and F_film ≥ F_VKH for every (γ, q). The map illustrates it.

**Neither of the book's two wordings of its assumption orders the error.**
Book p. 335, beside the chart, says the curves were computed "under the
assumption that B is only weakly depleted near the interface" — and adds "For
moderately fast reactions, this assumption was reasonably confirmed by more
rigorous computations." Book p. 340, five lines above (6.3.5-1) itself, says
they proceeded "by assuming that the concentration of B remains approximately
constant close to the interface". The first is about the *level* of depletion,
which β(0) measures; the second about the *constancy* of β where the reaction
happens, which β(0) does not measure at all — that is the rate-weighted mean
β̄ divided by β(0). Along the worst curve (q = 1.5217):

| | γ | β(0) | β̄/β(0) | VKH error |
|---|---|---|---|---|
| B half gone | 2.2036 | 0.5 | 1.208 | −1.8304 % |
| the worst point | 4.7279 | 0.1647 | **2.115** | **−2.7955 %** |
| B 99 % gone | 13.0378 | 0.01 | 19.20 | −1.5035 % |
| B essentially gone | 100 | 1.367e−11 | 3.684e+9 | −0.0373 % |
| B gone (γ = 1000) | 1000 | numerically zero | — | **−0.00037 %** |

β(0) falls monotonically, β̄/β(0) rises monotonically, and the error does
neither: it peaks in between and vanishes at both ends, because the
instantaneous ceiling F_A → 1 + q is a limit the approximation reproduces by
construction. **Froment's own sentence is qualified to *moderately fast*
reactions and, measured, it stands** — the worst point is inside that band and
errs by 2.8 %, well inside his printed 10 percent. What fails is the
unqualified reading, that weak depletion is *why* the approximation works.

At γ = 1000 β(0) is reported as numerically zero rather than as a bound with
digits: the solve returns a number of order 1e−58 whose sign is round-off.

**Levenspiel's three printed expansions, and what each expands.** The constants
check Froment's equation independently — except the last one:

| printed (Fig. 23.4, book p. 530) | measured |
|---|---|
| E ≅ 1 + M_H²/3 | within 1 % up to M_H = 0.7231 (root-found) |
| E = M_H(1 − (M_H−1)/2E_i), for E_i > 5M_H | −0.5366 % at the boundary, against +9.3004 % for E ≅ M_H |
| the same expansion of (6.3.5-1), with 2(E_i−1) | −0.7374 % — the two books' groups agree to 0.20 points there |
| E ≅ E_i − E_i²(E_i−1)/M_H², for E_i < M_H/5 | tracks **VKH** to 0.9953 of its deficit at M_H = 400; tracks the **film solve** to only 0.2286 |

**And the reason is exact, not empirical.** Putting F_A = 1 + q − d into
(6.3.5-1) gives d → q(1+q)²/γ² for large γ, which in Levenspiel's symbols *is*
E_i²(E_i−1)/M_H²: his correction is the algebraic asymptote of the
approximation. The film's own deficit is fixed by the first integral of the two
film equations, q·β(0) = 1 + q − F_A, so it *equals* q·β(0) — the interfacial
value of a reaction layer of thickness ~1/γ, exponentially small, and therefore
no power law at all. Its local exponent steepens octave by octave (−1.92,
−2.44, −3.34, **−5.34**) while the VKH deficit converges on the printed −2
(−1.9949 on the finest octave). At M_H = 800 the printed correction is
**44.24× too large** — and 44× of a 0.00134 % correction, so *leaving it out*
is there more accurate than applying it. **Two independent methods** — pymrm
finite volumes and scipy collocation, sharing no assembly — agree on that
deficit to 8.61e−5, the first of them also refined on a finer, differently
graded mesh. (That is two methods, one of them refined; it is not three
independent routes.)

**Reconciled with `F3.1`, which owns Hatta.** `F3.1` reports 2.076893 % as its
maximum VKH error. Its definition — a sampled max over E_i ∈ {5, 20, 100} on a
warm-started Ha sweep — recomputed here cold gives 2.076896 %, reproducing it
to 1.29e−6. That residual is **`F3.1`'s own grid error — not this page's
numerics, and not a warm-vs-cold difference either**, and both halves are
measured rather than asserted: run *cold* on `F3.1`'s own grid — uniform,
n_x = 400 — the same recomputation lands on 2.076893 %, reproducing `F3.1`'s
warm-started stored value to better than 1e−10 relative. The warm start
therefore contributes nothing to the residual, and what is left is the uniform
400-cell discretisation. This page's extrapolated reference is the accurate
side of it: it agrees with scipy's collocation to 4.37e−8. For scale, the
single graded n_x = 800 grid this page used before the baseline was
extrapolated lands 0.000214 points away, eighty times further. The stored
number is **loaded** from `F3.1`'s `agreement.json`, not retyped.

The gap to this page's headline is 0.7186 points, and it splits: root-finding γ
at the same three E_i buys **0.0016**, letting q move off {4, 19, 99} buys
**0.7170**. It is almost all q.

**What grading buys, and what it costs.** `F3.1`'s pseudo-first-order check
matches Ha/tanh(Ha) to 6.337e−3 on a uniform grid of 400 cells, "concentrated
at high Ha where the reaction layer is thinner than a cell". Run here it comes
out 7.088e−6 on the graded n_x = 800 grid — 893.9× tighter — but **that factor
is not what grading buys.** Split on this page's own three grids:

| | |
|---|---|
| uniform 400 → uniform 800 (refinement) | 7.163× |
| uniform 800 → graded 800 (grading, same n_x) | **124.8×** |
| product | 893.9× |

Crediting all 893.9× to grading, as an earlier version of this page did, is
wrong by a factor of seven. And at the moderate γ where the headline lives
grading is a **liability**: at the grid-study point (γ, q) = (10, 4) a graded
grid is **248.7× less accurate** than a uniform one of the same size, because
its outer cells are nearly eight times wider. Hence the extrapolated baseline,
which has the high-γ end and the fourth decimal at the same time.

## Data

**Provenance tier 6 — nothing here was measured.** The two CSVs are
transcriptions: the numeric constants, thresholds, axis ends and the two
printed wordings of the assumption, and the ten curve labels the charts share.
No curve is digitised and no point is read off either figure; every curve the
notebook plots is computed.

**The Froment file's equations do not survive text extraction**: every
Symbol-font operator becomes an unmappable Private-Use-Area glyph, so an
extracted equation looks complete and has no operators at all. Nothing here
came off that text layer. Its *equations and prose* are vector text at any
resolution — but **its figure is not**: `pdfimages -list` shows Fig. 6.3.2-1 on
PDF p. 373 as an embedded 536×771 greyscale JPEG at 148 ppi, so the 300 ppi
page render used for the curve labels was a 2× interpolation. The labels and
axis titles were re-read on the extracted raster at native resolution and are
unchanged. Levenspiel is a JBIG2 scan at a native 600 ppi and was rendered
there.

**One printed defect, reported and not repaired.** Levenspiel spells the first
author "van Krevelens" — with a terminal s — in both places he names him, the
caption of Fig. 23.4 and his chapter reference list. It is settled against
another document rather than against memory: Froment cites the same paper, same
volume 67, same page 563, and spells it *Van Krevelen*. (Crossref gives a third
capitalisation, "van Krevelen".) The two books' *citations* do not conflict —
Froment cites the 1948 paper, Levenspiel's caption credits a 1954 one, and his
reference list carries both.

## What this page cannot conclude

Nothing about absorbers. Both sides of every comparison are equations; the film
model itself is never tested against a measurement here. And the 1948 origin
was not read, so nothing here is evidence about what Van Krevelen and Hoftijzer
wrote — only about the equation two monographs attribute to them.
