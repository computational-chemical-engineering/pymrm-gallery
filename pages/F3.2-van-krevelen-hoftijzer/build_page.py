#!/usr/bin/env python3
"""Generate index.ipynb for page F3.2 (Van Krevelen-Hoftijzer enhancement).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Van Krevelen-Hoftijzer: the chart replaced by a solve, and the assumption behind it measured"
description: "The enhancement factor for a second-order gas-liquid reaction is an implicit relation that textbooks hand over as a chart. Solving the film equations instead prices the approximation everywhere on that chart - worst error -2.7955 % at gamma = 4.7279, q = 1.5217, comfortably inside the 10 percent Froment, De Wilde & Bischoff print for it - and tests the assumption they state for it, in both of the two wordings they print: that B is only weakly depleted near the interface (book p. 335) and that its concentration remains approximately constant close to the interface (book p. 340). Neither failure orders the error: at the worst point B at the interface is already down to 0.1647 of bulk and the rate-weighted mean of beta is 2.115 times its interfacial value, and where B is gone completely the approximation is exact to -0.00037 %. Levenspiel's three printed expansions of the same result check every constant independently, and one of them turns out to expand the approximation rather than the equations it approximates: his E = E_i - E_i^2(E_i-1)/M_H^2 is algebraically the asymptote of the Van Krevelen-Hoftijzer relation, which it tracks to 0.9953 of the deficit at M_H = 400, while the true deficit is q times beta(0) exactly and therefore exponentially small, so at M_H = 800 the printed correction is 44.24 times too big."
categories: [sec:F, struct:S3, tier:T0, data:tier6, phase:gas-liquid]
date: 2026-08-14
---

# Van Krevelen-Hoftijzer: the chart replaced by a solve, and the assumption behind it measured

**Catalog ID:** `F3.2` · **Structures:** `S3` (1-D steady BVP) · **Tier:** T0
"""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

A gas $A$ dissolving into a liquid that holds a reactant $B$ is absorbed faster
than it would be by dissolution alone. The ratio is the **enhancement factor**,
and for a genuinely second-order reaction it has no closed form: the two film
equations are coupled and nonlinear, and every design chart in the field
descends from one approximate solution of them.

That approximation is **Van Krevelen and Hoftijzer's (1948)**. This page does
two things with it that a chart cannot: it *solves* the implicit relation
instead of reading it off, and it *measures the assumption* the approximation
is built on.

### The source situation, stated plainly

**The 1948 origin was not consulted.** Van Krevelen, D. W. and Hoftijzer,
P. J., *Kinetics of gas-liquid reactions. Part I. General theory*, Recueil des
Travaux Chimiques des Pays-Bas **67**, 563 (1948), is not on disk and was not
read. Every equation on this page is transcribed from a monograph that
restates it with attribution, and the transcription is named equation by
equation below.

**The text read is Froment, De Wilde and Bischoff.** *Chemical Reactor
Analysis and Design*, 3rd edition, John Wiley & Sons (2011), ISBN
978-0-470-56541-4, Chapter 6 *Gas-Liquid Reactions*. It qualifies as a source
for a page *about this result* rather than merely a book that mentions it:
section 6.3.2 (book p. 331) says the diagram is "as first given by Van Krevelen
and Hoftijzer [1948]"; section 6.3.3 (book p. 335) says the intermediate curves
"were calculated by Van Krevelen and Hoftijzer [1948] under the assumption that
B is only weakly depleted near the interface"; and section 6.3.5 (book p. 340)
derives the approximation itself and prints its accuracy claim. Named,
attributed, derived, and its assumption stated.

**Levenspiel is a cross-check, not a second source.** *Chemical Reaction
Engineering*, 3rd edition, John Wiley & Sons (1999), Chapter 23, carries the
same result as a chart in a different notation with three asymptotic
expansions printed in its annotation boxes. Those expansions are used here to
check Froment's constants independently and to reconcile the two notations.
Nothing on this page adjudicates between the books.

### Every equation was read on a render, and here is why that matters

The Froment file is born-digital and its prose extracts cleanly, **but its
equations do not survive extraction and the failure is silent**. Every
Symbol-font operator becomes an unmappable Private-Use-Area glyph - `U+F02D`
for minus, `U+F03D` for equals, `U+F02B` for plus - which renders as nothing.
The extraction of the central equation of this page, (6.3.5-1), reads exactly

```
  '    1   FA   1

b D A C Ai
a DB C Bb
```

Seven characters have gone: `U+F067` (the leading $\gamma$), `U+F03D` (the
equals), `U+F067` again, `U+F02D` (a minus), `U+F028` (an opening bracket),
`U+F02D` (the second minus) and `U+F029`. The square root sign has no character
at all - it is drawn, not typeset - and the fraction bar is a line. What is left
looks like a complete equation and contains no operator whatsoever. Every
equation and every numeral used here was read on a 300 ppi render (Levenspiel:
600 ppi, his file's native resolution), cropped and enlarged.
"""))

# ---------------------------------------------------------- published model
cells.append(md(r"""## The published model

### The film equations (Froment, section 6.3.1)

A stagnant liquid film of thickness $y_L$ at the interface, everything beyond
it well mixed, and an irreversible reaction $aA + bB \rightarrow$ products.
With $\xi = y/y_L$, $\alpha = C_A/C_{Ai}$ and $\beta = C_B/C_{Bb}$,

$$
\frac{\mathrm{d}^2\alpha}{\mathrm{d}\xi^2} = \gamma^2\,\alpha\beta,
\qquad
\frac{\mathrm{d}^2\beta}{\mathrm{d}\xi^2} = \frac{\gamma^2}{q}\,\alpha\beta ,
$$

$$
\gamma = \frac{\sqrt{k D_A}}{k_L},
\qquad
q = \frac{a\,D_B\,C_{Bb}}{b\,D_A\,C_{Ai}} ,
$$

with $\alpha(0) = 1$ (saturated interface), $\alpha(1) = 0$ - Van Krevelen and
Hoftijzer "could set $C_{Ab} = 0$", book p. 340 - $\beta'(0) = 0$ ($B$ is
non-volatile) and $\beta(1) = 1$. The enhancement factor is the interfacial
flux divided by the purely diffusive one,

$$
F_A = \frac{N_A}{k_L (C_{Ai} - C_{Ab})}
= -\left.\frac{\mathrm{d}\alpha}{\mathrm{d}\xi}\right|_{\xi = 0} .
$$

$\gamma$ is Froment's symbol and his eq. (6.3.2-1) *et seq.*; he notes it "is
sometimes called the Hatta number" (book p. 328). $q$ is the group his
Fig. 6.3.2-1 carries as its right-hand axis title, and his eq. (6.3.3-10) makes
the instantaneous ceiling $F_A = 1 + q$.

### The approximation this page is about

**Pseudo-first-order**, exact when $\beta \equiv 1$ - Froment (6.3.2-11):

$$
F_A = \frac{\gamma}{\tanh\gamma}
$$

**Van Krevelen-Hoftijzer**, "entirely analogous with that obtained for a
pseudo-first-order reaction (6.3.2-11), but with $\gamma$ replaced by
$\gamma'$" - Froment (6.3.5-1), book p. 340, read on a 300 ppi crop enlarged
2x:

$$
\gamma' = \gamma\sqrt{1 - (F_A - 1)\,\frac{b\,D_A\,C_{Ai}}{a\,D_B\,C_{Bb}}}
\;=\; \gamma\sqrt{1 - \frac{F_A - 1}{q}} ,
\qquad
F_A = \frac{\gamma'}{\tanh\gamma'} .
$$

$F_A$ appears on both sides: this is the implicit relation, and solving it is
what the chart exists to avoid. The sentence that follows it in the book is the
claim this page measures: *"This approximate solution is valid to within 10
percent of the numerical solution."*

And the sentence three sections earlier, book p. 335, is the assumption this
page tests: the intermediate curves *"were calculated by Van Krevelen and
Hoftijzer [1948] under the assumption that B is only weakly depleted near the
interface. For moderately fast reactions, this assumption was reasonably
confirmed by more rigorous computations."*

### The same result in Levenspiel's notation

| Froment | Levenspiel | this page |
|---|---|---|
| $\gamma = \sqrt{kD_A}/k_L$ | $M_H = \sqrt{\mathcal{D}_A k C_B}/k_{Al}$ | `gamma` |
| $F_A$ | $E$ | `f_a` |
| $q = aD_BC_{Bb}/(bD_AC_{Ai})$ | $E_i - 1$, from $E_i = 1 + \mathcal{D}_BC_BH_A/(b\mathcal{D}_Ap_{Ai})$ | `q` |
| $F_A \to 1 + q$ | $E \to E_i$ | the plateau |

The two groups are the same objects written differently - $p_{Ai} = C_{Ai}/H_A$
turns Levenspiel's $E_i$ into Froment's $1 + q$ - **but the curve labels on the
two charts are not.** Froment's Fig. 6.3.2-1 labels its curves with $q$;
Levenspiel's Fig. 23.4 labels the same family with $E_i$. Both print the same
ten numerals. A curve labelled 1 plateaus at $F_A = 2$ in one book and at
$E = 1$ in the other, and the page computes what that costs.

Levenspiel's Fig. 23.4 also prints three asymptotic expansions in its
annotation boxes, and every constant in them is checked here against the
implicit relation above and against the film equations themselves:

$$
E \cong 1 + \frac{M_H^2}{3} + \dots
\qquad
E = M_H\left(1 - \frac{M_H - 1}{2E_i}\right) + \dots
\qquad
E \cong E_i - \frac{E_i^2(E_i-1)}{M_H^2} + \dots
$$
"""))

# ------------------------------------------------------ parameters and scope
cells.append(md(r"""## Parameters and assumptions

**Assumptions**, all of them the source's: steady state; a stagnant film with
no surface renewal; irreversible second-order kinetics, first order in each
reactant; constant diffusivities; $C_{Ab} = 0$, so no reaction in the bulk;
$B$ non-volatile. The problem is then dimensionless and has exactly two
parameters, $\gamma$ and $q$.

**Nothing here is fitted, and nothing here is measured.** Both sides of every
comparison on this page are equations: an approximation on one side and the
equations it approximates on the other. This is provenance tier 6. The page
therefore establishes nothing whatever about how well the film model describes
a real absorber - only how well one solution of it describes another.

**The search domain**, and each end taken from the row that actually carries
it. $\gamma \in [0.5, 1000]$ is Froment's own printed $\gamma$ axis,
transcribed. The $q$ range comes from the printed **curve labels** - $q$ is
what those labels are - so its upper end is the largest of them, and its lower
end goes one decade below the smallest, because $q$ is a physical ratio and not
a drawn curve. The two ends happen to coincide numerically with the $\gamma$
axis; they are looked up from different rows all the same. The worst point
turns out to lie in the interior of both, so widening the domain cannot move
it.

**Where this page stops, and what `F3.1` owns.** `F3.1` owns the Hatta number
and the regime map, and it already reports a maximum Van Krevelen-Hoftijzer
error of about 2.1 % over the three $E_i$ values it sweeps. That number is
loaded from `F3.1`'s own `agreement.json` in the Validation section, recomputed
under this page's definitions, and reconciled - not retyped. This page does not
restate the regime picture, and it does not touch DeCoursey, which is `F3.3`.
"""))

# ------------------------------------------------------------- colab + setup
cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import json
import sys
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

from pymrm import (NumJac, clip_approach, compute_boundary_values,
                   construct_div, construct_grad, newton)

RAW = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
       "pymrm-gallery/main")
if "google.colab" in sys.modules:
    urllib.request.urlretrieve(RAW + "/shared/gallery_utils.py", "gallery_utils.py")
else:
    for _p in (Path.cwd(), *Path.cwd().parents):
        if (_p / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(_p / "shared"))
            break
from gallery_utils import cite_data, load_data, load_meta, report_agreement

PAGE = "F3.2-van-krevelen-hoftijzer"
plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5, "axes.grid": True,
                     "grid.alpha": 0.25})
# Okabe-Ito, assigned in fixed order and never cycled
C_BLUE, C_ORANGE, C_GREEN = "#0072B2", "#D55E00", "#009E73"
C_PURPLE, C_YELLOW, C_GREY = "#CC79A7", "#E69F00", "0.45"

# DETERMINISM: nothing on this page is stochastic.  No sampling, no random
# initial guess, and no warm start - every film solve begins from a closed-form
# profile that depends on nothing but its own (gamma, q).  Two consecutive
# executions give an identical agreement.json.
np.set_printoptions(precision=8, suppress=False)'''))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

There is no experimental dataset here and there could not be: the case is an
approximation to a pair of differential equations. What the two CSVs carry is
**transcription** - the numbers and sentences the two books print around this
result - so that the notebook tests them instead of retyping them.

Both books were read on renders at their files' native resolutions, and the
page/PDF mapping is recorded in the sidecars so any reader can re-render
exactly what was read (Froment: book = PDF $-$ 42; Levenspiel: book = PDF
$-$ 16, both confirmed on printed running heads).

One correction to what the sidecars used to say, because it is the exact trap
the house rule exists for. The Froment file is born-digital and its equations
and prose are vector text at any resolution - **but its figure is not**.
`pdfimages -list` shows Fig. 6.3.2-1, on PDF p. 373, as an embedded
$536\times771$ greyscale JPEG at 148 ppi, so the 300 ppi page render the ten
curve labels were first read on was a 2x interpolation of a 148 ppi raster. The
raster was extracted and the labels, the four axis ends and both axis titles
re-read at native resolution enlarged 4x: they are unchanged. Only the claim
was wrong, not the data.
"""))

cells.append(code('''PRN = load_data("vkh-printed-claims.csv", page=PAGE).set_index("claim_id")
PMETA = load_meta("vkh-printed-claims.csv", page=PAGE)
LAB = load_data("enhancement-chart-labels.csv", page=PAGE)
LMETA = load_meta("enhancement-chart-labels.csv", page=PAGE)


def printed(key, field="value"):
    """A printed constant, looked up - never retyped in a cell."""
    return float(PRN.loc[key, field])


def quoted(key):
    return str(PRN.loc[key, "verbatim"])


print(cite_data(PMETA))
print(f"{len(PRN)} transcribed rows, {len(LAB)} chart labels\\n")
for k in ("froment_vkh_validity_pct", "froment_film_completion_gamma",
          "froment_bulk_reaction_gamma"):
    print(f'  {k:32s} {PRN.loc[k, "value"]:>8}   book p.{PRN.loc[k, "book_page"]}'
          f'  "{quoted(k)}"')
for k in ("lev_small_M_coeff", "lev_pfo_correction_denom", "lev_pfo_correction_offset",
          "lev_inst_correction_power", "lev_pfo_threshold", "lev_inst_threshold"):
    print(f'  {k:32s} {PRN.loc[k, "value"]:>8}   book p.{PRN.loc[k, "book_page"]}'
          f'  "{quoted(k)}"')

# ---- the domain of every search on this page, read off the printed axes -----
# gamma comes from the printed gamma AXIS; q comes from the printed curve
# LABELS, which is what q is drawn as - two different transcribed sources, and
# each is looked up from the row that actually carries it.
LABELS = [int(v) for v in LAB["label"]]
G_LO, G_HI = printed("froment_gamma_axis_lo"), printed("froment_gamma_axis_hi")
Q_HI = float(max(LABELS))
Q_LO = float(min(LABELS)) / 10.0   # one decade below the smallest printed label
print(f"\\nsearch domain: gamma in [{G_LO:g}, {G_HI:g}] (Froment's printed gamma axis),"
      f"  q in [{Q_LO:g}, {Q_HI:g}]")
print(f"printed curve labels, both books: {LABELS}")
print(f"  q's range comes from those labels, not from the gamma axis: the largest"
      f" printed label is {max(LABELS):g}"
      f"\\n  and the search goes one decade below the smallest, {min(LABELS):g}, because"
      f" q is a physical ratio")
print(f"  and not a drawn curve.  The worst point turns out to be interior in both"
      f" directions, so"
      f"\\n  widening the domain cannot move it.")'''))

cells.append(md(r"""### One printed defect, reported and not repaired

Levenspiel spells the first author **"van Krevelens"**, with a terminal *s*, in
both places he names him. It is settled against another document rather than
against memory: Froment cites the same paper - same journal, same volume 67,
same page 563 - and spells it *Van Krevelen*.

This is a defect in the cross-check source, not in the E1.1 target, and it
changes nothing computed here. The two books' *citations* do not conflict:
Froment cites the 1948 *Rec. Trav. Chim.* paper, Levenspiel's chart caption
credits a 1954 *Trans. I. Chem. E.* paper, and Levenspiel's own reference list
carries **both**. Those are two publications, not a discrepancy. Neither is on
disk.
"""))

cells.append(code('''for k in ("froment_surname_in_text", "froment_surname_in_reflist",
          "levenspiel_surname_in_caption", "levenspiel_surname_in_reflist"):
    print(f'{PRN.loc[k, "source"]:11s} book p.{PRN.loc[k, "book_page"]:>4}  "{quoted(k)}"')

FROMENT_SPELLING = "Van Krevelen"
LEV_SPELLING = "van Krevelens"
SPELLING_DIFFERS = (FROMENT_SPELLING.lower() not in
                    [quoted("levenspiel_surname_in_caption").split()[i] + " " +
                     quoted("levenspiel_surname_in_caption").split()[i + 1]
                     for i in range(len(quoted("levenspiel_surname_in_caption").split()) - 1)])
IN_FROMENT = sum(FROMENT_SPELLING in quoted(k)
                 for k in ("froment_surname_in_text", "froment_surname_in_reflist"))
IN_LEV = sum(LEV_SPELLING in quoted(k)
             for k in ("levenspiel_surname_in_caption", "levenspiel_surname_in_reflist"))
SAME_PAPER = all(s in quoted("levenspiel_surname_in_reflist") for s in ("67", "563", "1948"))
print(f"\\n'{FROMENT_SPELLING}' in {IN_FROMENT} of the 2 Froment rows;"
      f"  '{LEV_SPELLING}' in {IN_LEV} of the 2 Levenspiel rows.")
print(f"both reference lists carry volume 67, page 563, year 1948: {SAME_PAPER}")
assert IN_FROMENT == 2 and IN_LEV == 2 and SAME_PAPER
print("-> one paper, two spellings.  This page uses the spelling without the s,"
      "\\n   which is the one the E1.1 target and the catalogue use.")'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Three solvers, and they are kept apart on purpose.

1. **The film equations, in pymrm.** Finite volumes on a slab (`nu=0`) with a
   graded grid clustered at the interface, per-species boundary conditions in
   one pair of `{"a","b","d"}` dictionaries, and `newton` + `NumJac` on the
   pointwise reaction term. This is the "numerical solution" Froment's accuracy
   claim is made against.
2. **The same equations by collocation.** `scipy.integrate.solve_bvp`, adaptive
   mesh, fourth order, sharing no assembly, no Jacobian and no mesh with (1).
   Every headline that depends on the film solution is computed both ways.
3. **The implicit relation**, solved twice: for $F_A$ by Brent, and again for
   $\gamma'$ with $F_A$ eliminated - a different unknown and a different
   residual for the same pair of printed equations.

Two details that are not decoration. The grid is **graded**: with
$x_f = (\mathrm{e}^{s t} - 1)/(\mathrm{e}^{s} - 1)$ and $s = 8$ the reaction
layer near $\xi = 0$ is resolved at $\gamma = 1000$, where a uniform grid of
the same size is not - `construct_grad` and `construct_div` take arbitrary face
positions, so this costs one line. And `NumJac` is given `(n_x, 2)` and never a
bare `(n_x,)`: the reaction term is pointwise in the two fields, which is
exactly the default stencil, while a bare 1-D shape would declare every cell
coupled to every other and build a dense Jacobian.

**The reference is not a single solve.** Grading is not free: it buys the
high-$\gamma$ end and it *costs* accuracy at moderate $\gamma$, which is
exactly where this page's headline sits, because the outer cells of a graded
grid are several times wider than a uniform one's. Both effects are measured
below. So the baseline every error here is quoted against is a **Richardson
extrapolation of the same graded film on $n_x = 400$ and $n_x = 800$**,
$F(2n) + [F(2n) - F(n)]/3$, which the observed order 2.0000 licenses. A single
graded solve at $n_x = 800$ is *not* converged at the fourth decimal of the
headline; the break table's "extrapolation switched off" row is that defect,
and the extrapolated baseline is checked against collocation and against a
finer extrapolation in the Validation section.

Every solve is a **cold start** from a closed-form profile that depends only on
its own $(\gamma, q)$. Nothing on this page is continued from a neighbouring
solve, so no reported number can inherit a warm-start history.
"""))

cells.append(code('''XTOL = dict(xtol=1e-13, rtol=8.9e-16)


# ------------------------------------------- Froment (6.3.2-11) + (6.3.5-1)
def f_a_vkh(gamma, q, sign=-1.0, offset=1.0, denom="q", xtol=None):
    """F_A = gamma'/tanh(gamma') with gamma' from (6.3.5-1), solved for F_A.

    `sign`, `offset`, `denom` and `xtol` exist so that the break table can
    mistranscribe the printed equation or under-solve it; the page always calls
    this with the printed form and the tight tolerance.
    """
    d = q if denom == "q" else 1.0 + q
    tol = XTOL if xtol is None else dict(xtol=xtol, rtol=8.9e-16)

    def resid(F):
        s = np.sqrt(max(1.0 + sign * (F - offset) / d, 0.0))
        g = gamma * s
        return F - (g / np.tanh(g) if g > 1e-8 else 1.0 + g * g / 3.0)

    hi = 1.0 + q
    while resid(hi) < 0 and hi < 1e12:      # a broken form can exceed the ceiling
        hi *= 2.0
    return brentq(resid, 1.0, hi, **tol, maxiter=300)


def f_a_vkh_route_b(gamma, q):
    """The SAME two printed equations, solved for gamma' instead of for F_A.

    F = u/tanh u and F = 1 + q(1 - u^2/gamma^2) with u = gamma'; eliminating F
    leaves one scalar equation in u.  Different unknown, different residual.
    """
    def resid(u):
        return u / np.tanh(u) - (1.0 + q * (1.0 - (u / gamma) ** 2))

    u = brentq(resid, 1e-12, gamma, **XTOL, maxiter=300)
    return u / np.tanh(u)


# ------------------------------------------------------- the film, in pymrm
def faces(n, stretch):
    """Faces graded towards xi = 0, the interface, where the reaction layer is."""
    t = np.linspace(0.0, 1.0, n + 1)
    return t if stretch <= 1e-12 else np.expm1(stretch * t) / np.expm1(stretch)


class Film:
    """Stagnant film on xi in [0, 1], fields (alpha, beta) = (A, B)."""

    def __init__(self, n_x=800, stretch=8.0, bc_beta_bulk="dirichlet"):
        self.n_x, self.stretch = n_x, stretch
        self._memo = {}
        self.shape = (n_x, 2)
        self.x_f = faces(n_x, stretch)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        # a . dc/dn + b . c = d, with n the OUTWARD normal.
        right = ({"a": [[0.0, 0.0]], "b": [[1.0, 1.0]], "d": [[0.0, 1.0]]}
                 #    alpha = 0 (C_Ab = 0)     beta = 1 (bulk value)
                 if bc_beta_bulk == "dirichlet" else
                 {"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[0.0, 0.0]]})
        #        the break-table alternative: dbeta/dn = 0 at the bulk edge
        self.bc = ({"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[1.0, 0.0]]}, right)
        #            alpha = 1 (saturated)     dbeta/dn = 0 (B non-volatile)
        grad_mat, grad_bc = construct_grad(self.shape, self.x_f, self.x_c,
                                           self.bc, axis=0)
        div_mat = construct_div(self.shape, self.x_f, nu=0, axis=0)   # slab
        self.flux_mat, self.flux_bc = -grad_mat, -grad_bc
        self.jac_diff = div_mat @ self.flux_mat
        self.g_diff = (div_mat @ self.flux_bc).toarray().ravel()
        # F_A lives in the first face row of the flux operator
        self.iface = np.asarray(self.flux_mat[[0], :].todense()).ravel()
        self.iface_bc = float(self.flux_bc.toarray().ravel()[0])
        # (n_x, 2), never a bare (n_x,): the rate is pointwise in the fields
        self.numjac = NumJac((n_x, 2))

    def enhancement(self, u):
        return float(self.iface @ np.asarray(u).ravel() + self.iface_bc)

    def solve(self, gamma, q, frozen_b=False, tol=1e-10):
        def rate(u):
            # frozen_b is the break table's "no depletion of B at all": B is held
            # at its bulk value in the rate AND its own consumption is dropped,
            # so beta = 1 identically and the film IS pseudo-first-order.
            b = np.ones_like(u[:, 1]) if frozen_b else np.clip(u[:, 1], 0, None)
            r = np.clip(u[:, 0], 0, None) * b
            return np.stack([gamma**2 * r,
                             np.zeros_like(r) if frozen_b else (gamma**2 / q) * r],
                            axis=1)

        # COLD START from a closed form, not from a neighbouring solve: the
        # pseudo-first-order profile at this same gamma.
        g = min(gamma, 700.0)                      # sinh overflows past ~710
        u0 = np.stack([np.sinh(g * (1.0 - self.x_c)) / np.sinh(g),
                       np.ones_like(self.x_c)], axis=1)

        def residual(u):
            g_r, jac_r = self.numjac(rate, u)
            return (self.jac_diff @ u.reshape((-1, 1))
                    + self.g_diff.reshape((-1, 1))
                    + g_r.reshape((-1, 1))), self.jac_diff + jac_r

        res = newton(residual, u0, tol=tol, maxfev=200,
                     callback=lambda x, r: clip_approach(x, r, 0.0, None))
        assert res.success, f"film solve failed at gamma={gamma}, q={q}"
        return np.asarray(res.x).reshape(self.shape)

    def beta_interface(self, u):
        """beta at xi = 0 through compute_boundary_values, not off a cell centre."""
        v, _ = compute_boundary_values(u, self.x_f, self.x_c, self.bc[0], axis=0,
                                       bound_id=0)
        return float(np.asarray(v).ravel()[1])

    def beta_constancy(self, gamma, q, **kw):
        """Rate-weighted mean of beta over the reaction zone, over beta(0).

        The book states its assumption TWICE and the two statements are not the
        same test.  Book p. 335 is about the LEVEL of depletion at the
        interface, which beta(0) measures.  Book p. 340 - the page the equation
        itself is printed on - is about beta staying "approximately constant
        close to the interface", which beta(0) does not measure at all.  This
        is the second test: weight beta by the local rate alpha*beta, where the
        approximation is actually being made, and divide by beta(0).  It is 1
        exactly when beta does not vary over the reaction zone.
        """
        u = self.solve(gamma, q, **kw)
        a = np.clip(u[:, 0], 0.0, None)
        b = np.clip(u[:, 1], 0.0, None)
        w = np.diff(self.x_f) * a * b
        s, b0 = float(w.sum()), self.beta_interface(u)
        # A break row CAN make this undefined rather than large: the
        # zero-gradient-bulk row removes the supply of B entirely, so beta and
        # the rate go to zero everywhere and there is no reaction zone to
        # average over.  Return nan explicitly, so the row records no move on
        # this metric instead of raising a warning and dividing 0 by 0.
        if not (s > 0.0) or b0 == 0.0:
            return float("nan")
        return float((w * b).sum() / s) / b0

    def f_a(self, gamma, q, **kw):
        # A pure-function memo.  f_a depends on nothing but (gamma, q, kw) and
        # the operators built once in __init__, and the root-finders ask for the
        # same point repeatedly.  It changes no number - only the runtime.
        key = (gamma, q, tuple(sorted(kw.items())))
        if key not in self._memo:
            u = self.solve(gamma, q, **kw)
            self._memo[key] = (self.enhancement(u), self.beta_interface(u))
        return self._memo[key]


class Extrapolated:
    """Richardson extrapolation of the SAME film on two grids: the reference.

    The graded discretisation converges at order 2 - observed 2.0000 in the
    Validation section - so F(2n) + [F(2n) - F(n)]/3 removes its leading error
    term.  Everything this page reports about the approximation is measured
    against this, not against a single solve, because a single graded solve at
    n_x = 800 is wrong in the fourth decimal of the headline error.  The break
    table's "extrapolation switched off" row is exactly that defect.
    """

    def __init__(self, coarse, fine, order=2.0):
        self.coarse, self.fine = coarse, fine
        self.n_x, self.stretch = fine.n_x, fine.stretch
        self.w = 1.0 / (2.0**order - 1.0)

    def f_a(self, gamma, q, **kw):
        f_c, b_c = self.coarse.f_a(gamma, q, **kw)
        f_f, b_f = self.fine.f_a(gamma, q, **kw)
        return f_f + (f_f - f_c) * self.w, b_f + (b_f - b_c) * self.w

    def solve(self, gamma, q, **kw):
        """Profiles are drawn from the fine grid; only scalars are extrapolated."""
        return self.fine.solve(gamma, q, **kw)

    def beta_constancy(self, gamma, q, **kw):
        # A ratio of two integrals of the SAME profile, so the grid partly
        # cancels - but only partly, and least of all at the top of the gamma
        # range, where beta is a thin exponential and the rate weight sits in
        # the first few cells.  An earlier version of this page took this
        # diagnostic from the fine grid alone on the ground that the grid
        # "cancels to five digits"; that is true at the worst point and false
        # by orders of magnitude at gamma = 100.  The gap is PRINTED below
        # rather than asserted here, and the ratio is extrapolated exactly like
        # F_A: same two grids, same order.
        c = self.coarse.beta_constancy(gamma, q, **kw)
        f = self.fine.beta_constancy(gamma, q, **kw)
        return f + (f - c) * self.w


# ------------------------------------------------ the film, by collocation
def f_a_bvp(gamma, q, tol=1e-10, n0=200):
    """The same two-point problem, adaptive 4th-order collocation."""
    def fun(x, y):
        r = gamma**2 * y[0] * y[2]
        return np.vstack([y[1], r, y[3], r / q])

    def bc(ya, yb):
        return np.array([ya[0] - 1.0, yb[0], ya[3], yb[2] - 1.0])

    x = np.expm1(6.0 * np.linspace(0, 1, n0)) / np.expm1(6.0)
    y = np.zeros((4, x.size))
    y[0], y[1], y[2] = 1.0 - x, -1.0, 1.0
    s = solve_bvp(fun, bc, x, y, tol=tol, max_nodes=400000)
    assert s.status == 0, f"solve_bvp failed at gamma={gamma}, q={q}: {s.message}"
    return float(-s.y[1][0]), float(s.sol(0.0)[2])


# -------------------------------------------------------- root-find helpers
def argext_log(f, lo, hi, n_scan=13, h=1e-3, rootfind=True):
    """An interior extremum of f on [lo, hi], ROOT-FOUND from f' = 0.

    The scan only brackets; the value returned is f at a root of the central
    difference of f in ln x.  Returns (x, f(x), was_root_found).
    """
    xs = np.geomspace(lo, hi, n_scan)
    fs = np.array([f(x) for x in xs])
    i = int(np.argmax(np.abs(fs)))
    if not rootfind or i in (0, n_scan - 1):
        return float(xs[i]), float(fs[i]), False

    def dfd(lx):
        x = np.exp(lx)
        return f(x * np.exp(h)) - f(x * np.exp(-h))

    a, b = np.log(xs[i - 1]), np.log(xs[i + 1])
    if dfd(a) * dfd(b) >= 0:
        return float(xs[i]), float(fs[i]), False
    lx = brentq(dfd, a, b, **XTOL, maxiter=200)
    return float(np.exp(lx)), float(f(np.exp(lx))), True


def root_log(f, lo, hi):
    return float(np.exp(brentq(lambda lx: f(np.exp(lx)), np.log(lo), np.log(hi),
                               **XTOL, maxiter=200)))


def root_or_edge(f, lo, hi):
    """root_log where a root exists, and the nearer end where it does not.

    Only the break table ever reaches the fallback: the baseline's two
    depletion roots are genuine, and the notebook prints their residuals.
    """
    a, b = f(lo), f(hi)
    if a * b > 0:
        return float(hi if abs(b) < abs(a) else lo)
    return root_log(f, lo, hi)


_GRIDS = {}


def graded_film(n_x, stretch):
    """One Film object per (n_x, stretch), so its solve memo is shared."""
    return _GRIDS.setdefault((n_x, float(stretch)), Film(n_x=n_x, stretch=stretch))


FILM = graded_film(800, 8.0)                 # the fine grid of the pair
FILM_C = graded_film(400, 8.0)               # its coarse partner
REF = Extrapolated(FILM_C, FILM)             # THE BASELINE OF THIS PAGE
REF_FINE = Extrapolated(FILM, graded_film(1600, 8.0))   # only to check REF
print(f"film: n_x = {FILM.n_x}, grading s = {FILM.stretch:g},"
      f" first cell {FILM.x_f[1]:.3e} wide, last {1 - FILM.x_f[-2]:.3e}")
print(f"the reference is Richardson(n_x = {FILM_C.n_x}, {FILM.n_x}) on that grid,"
      f" order 2, NOT a single solve")
print(f"F_A at (gamma, q) = (10, 4):  pymrm n_x = 800 {FILM.f_a(10.0, 4.0)[0]:.10f}"
      f"   extrapolated {REF.f_a(10.0, 4.0)[0]:.10f}")
print(f"{'':29s}collocation {f_a_bvp(10.0, 4.0)[0]:.10f}"
      f"   VKH {f_a_vkh(10.0, 4.0):.10f}")'''))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 1. The chart, computed

Froment's Fig. 6.3.2-1 and Levenspiel's Fig. 23.4 are the same picture: the
implicit relation, solved once and drawn as a family of curves so that nobody
has to solve it again. Here it is solved - the curves below are computed from
the equations at the ten labels both books print, and nothing is traced.

The plateau check is the point of the notation table above: at the *same*
printed label, Froment's curve levels off at $1 + q$ and Levenspiel's at $E_i$.
"""))

cells.append(code('''GAMMA = np.geomspace(G_LO, G_HI, 60)
CURVES = {L: np.array([f_a_vkh(g, float(L)) for g in GAMMA]) for L in LABELS}
SOLVED = {L: np.array([REF.f_a(g, float(L))[0] for g in GAMMA[::6]])
          for L in (1, 10, 100, 1000)}

fig, ax = plt.subplots(figsize=(6.4, 5.0))
for L in LABELS:
    ax.loglog(GAMMA, CURVES[L], color=C_BLUE, lw=1.4)
    ax.annotate(f"{L}", (GAMMA[-1], CURVES[L][-1]), fontsize=7.5,
                xytext=(3, -2), textcoords="offset points", color=C_BLUE)
for L, mk in zip((1, 10, 100, 1000), ("o", "s", "^", "D")):
    ax.loglog(GAMMA[::6], SOLVED[L], mk, ms=4.5, mfc="none", mew=1.1,
              color=C_ORANGE, label=f"film solve, q = {L}" if L == 1 else None)
ax.loglog(GAMMA, GAMMA / np.tanh(GAMMA), color=C_GREY, lw=1.1, ls="--",
          label=r"$\\gamma/\\tanh\\gamma$  (6.3.2-11)")
ax.set_xlim(G_LO, G_HI * 1.5)
ax.set_ylim(printed("froment_FA_axis_lo"), printed("froment_FA_axis_hi"))
ax.set_xlabel(r"$\\gamma = \\sqrt{kD_A}/k_L$   ($M_H$ in Levenspiel's notation)")
ax.set_ylabel(r"$F_A$   ($E$ in Levenspiel's notation)")
ax.set_title("The enhancement chart, solved rather than read\\n"
             "curves: VKH (6.3.5-1) at the ten printed labels, as q",
             fontsize=9.5)
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout()
plt.show()

# ---- the plateau each book's own label implies -----------------------------
PLATEAU_FR = np.array([1.0 + L for L in LABELS])     # Froment: label = q
PLATEAU_LEV = np.array([float(L) for L in LABELS])   # Levenspiel: label = E_i
LABEL_RATIO = PLATEAU_FR / PLATEAU_LEV
CHART_LABEL_MAX_RATIO = float(LABEL_RATIO.max())
print("the same printed label, the two books' plateaux")
print(f'{"label":>7}{"Froment 1+q":>13}{"Levenspiel E_i":>16}{"ratio":>9}')
for L, a, b, r in zip(LABELS, PLATEAU_FR, PLATEAU_LEV, LABEL_RATIO):
    print(f"{L:7d}{a:13.0f}{b:16.0f}{r:9.4f}")
print(f"\\nworst at the smallest label: {CHART_LABEL_MAX_RATIO:.4g}x"
      f"   ({100 * (CHART_LABEL_MAX_RATIO - 1):.4g} %);"
      f" least at the largest: {LABEL_RATIO.min():.4f}")
print(f"the same abscissa is drawn over different ranges: Froment"
      f" {printed('froment_gamma_axis_lo'):g} to {printed('froment_gamma_axis_hi'):g},"
      f" Levenspiel {printed('lev_MH_axis_lo'):g} to {printed('lev_MH_axis_hi'):g}"
      f"\\n(this page searches Froment's, the narrower one)")
# the computed curve must actually reach the plateau its label implies
PLATEAU_ERR = max(abs(f_a_vkh(1e6, float(L)) / (1.0 + L) - 1) for L in LABELS)
print(f"computed VKH curves reach 1 + q at gamma = 1e6 to {PLATEAU_ERR:.2e} relative")'''))

cells.append(md(r"""### 2. Where the approximation is worst, root-found

The error $\varepsilon(\gamma, q) = F_A^{\mathrm{VKH}}/F_A^{\mathrm{film}} - 1$
is a smooth surface with one interior extremum. It is located by root-finding
$\partial\varepsilon/\partial\ln\gamma = 0$ at fixed $q$, and then
$\mathrm{d}\varepsilon^\*/\mathrm{d}\ln q = 0$ over $q$ - **not** by taking the
largest value on a sweep. The difference is not cosmetic: the same search with
the root-finds disabled is one of the rows in the break table below.
"""))

cells.append(code('''def err_pct(g, q, film=None, vkh=None, **kw):
    """100 (F_VKH/F_film - 1) at one point of the chart."""
    film = REF if film is None else film
    vkh = f_a_vkh if vkh is None else vkh
    F, _ = film.f_a(g, q, **kw)
    return 100.0 * (vkh(g, q) - F) / F


def block_worst(film=None, vkh=None, rootfind=True, n_scan=13, labels=None, **kw):
    """The worst point of the chart, and the depletion diagnostics there."""
    film = REF if film is None else film
    labels = LABELS if labels is None else labels

    def eps(g, q):
        return err_pct(g, q, film=film, vkh=vkh, **kw)

    def worst_at(q):
        return argext_log(lambda g: eps(g, q), G_LO, G_HI, n_scan,
                          rootfind=rootfind)[1]

    q_star, _, ok_q = argext_log(worst_at, Q_LO, Q_HI, 9, rootfind=rootfind)
    g_star, e_star, ok_g = argext_log(lambda g: eps(g, q_star), G_LO, G_HI,
                                      n_scan, rootfind=rootfind)
    out = {"vkh_worst_rel_error_pct": e_star,
           "vkh_worst_gamma": g_star,
           "vkh_worst_q": q_star,
           "froment_10pct_claim_margin":
               printed("froment_vkh_validity_pct") / abs(e_star),
           "beta_interface_at_worst": film.f_a(g_star, q_star, **kw)[1],
           # the p. 340 wording of the same assumption, measured where the
           # p. 335 wording is most badly violated
           "beta_constancy_at_worst": film.beta_constancy(g_star, q_star, **kw),
           "vkh_worst_on_printed_labels_pct":
               float(max((argext_log(lambda g: eps(g, float(L)), G_LO, G_HI,
                                     n_scan, rootfind=rootfind)[1]
                          for L in labels), key=abs)),
           "_root_found": bool(ok_q and ok_g)}

    def beta_at(g):
        return film.f_a(g, q_star, **kw)[1]

    g_half = root_or_edge(lambda g: beta_at(g) - 0.5, G_LO, G_HI)
    g_99 = root_or_edge(lambda g: beta_at(g) - 0.01, G_LO, G_HI)
    out.update({"gamma_at_half_depletion": g_half,
                "gamma_at_99pct_depletion": g_99,
                "vkh_error_at_half_depletion_pct": eps(g_half, q_star),
                "vkh_error_at_99pct_depletion_pct": eps(g_99, q_star),
                "vkh_error_at_full_depletion_pct": eps(G_HI, q_star),
                "beta_constancy_at_99pct_depletion":
                    film.beta_constancy(g_99, q_star, **kw),
                "_beta_constancy_at_half": film.beta_constancy(g_half, q_star, **kw),
                "_beta_at_full_depletion": beta_at(G_HI)})
    return out


W = block_worst()
G_STAR, Q_STAR = W["vkh_worst_gamma"], W["vkh_worst_q"]
print(f"worst error on the chart: {W['vkh_worst_rel_error_pct']:+.4f} %"
      f"  at gamma = {G_STAR:.4f}, q = {Q_STAR:.4f}"
      f"  (E_i = 1 + q = {1 + Q_STAR:.4f})")
print(f"  found by root-finding, not by sampling: {W['_root_found']}")
print(f"  worst over the ten PRINTED curve labels only:"
      f" {W['vkh_worst_on_printed_labels_pct']:+.4f} %")
print(f"  the book's printed claim is {printed('froment_vkh_validity_pct'):g} percent;"
      f" margin {W['froment_10pct_claim_margin']:.4f}x")
print(f"\\nVKH is an UNDER-estimate everywhere it was evaluated:"
      f" sign of the worst error {np.sign(W['vkh_worst_rel_error_pct']):+.0f}")
for q in (0.05, 0.1, 0.5, Q_STAR, 5.0, 50.0, 500.0):
    e = argext_log(lambda g: err_pct(g, q), G_LO, G_HI, 13)[1]
    print(f"   q = {q:8.4f}   worst over gamma {e:+8.4f} %")'''))

cells.append(md(r"""### 3. The assumption, measured - in both of the book's wordings

**The book states the assumption twice, and the two statements are not the same
test.** Book p. 335, beside the chart, says the intermediate curves *"were
calculated by Van Krevelen and Hoftijzer [1948] under the assumption that B is
only weakly depleted near the interface"* - and adds, in the same breath,
*"For moderately fast reactions, this assumption was reasonably confirmed by
more rigorous computations."* Book p. 340, five lines above (6.3.5-1) itself,
says they proceeded *"by assuming that the concentration of B remains
approximately constant close to the interface"*.

The first wording is about the **level** of depletion at the interface, which
$\beta(0)$ measures exactly. The second is about the **constancy** of $\beta$
over the zone where the reaction happens, which $\beta(0)$ does not measure at
all: it is measured here by the rate-weighted mean
$\bar\beta = \int \alpha\beta\cdot\beta\,\mathrm{d}\xi \big/
\int \alpha\beta\,\mathrm{d}\xi$ divided by $\beta(0)$, which is 1 when
$\beta$ is constant where it matters. Both are printed below, so the finding
does not depend on which wording is tested.

**Neither survives as an explanation of the error.** At the worst point $B$ at
the interface is already down to a sixth of its bulk value and $\bar\beta$ is
twice $\beta(0)$ - neither wording is nearly satisfied there - and yet as
$\gamma$ grows further, with $B$ consumed *completely* and $\bar\beta/\beta(0)$
running away without bound, the approximation becomes **exact**. The error
vanishes at both ends and peaks in between, and the two assumption measures are
monotone in $\gamma$, so no monotone function of either can order it.

The reason is visible in the profiles: at large $\gamma$ the reaction collapses
onto a plane and $F_A \to 1 + q$, a limit the approximation reproduces by
construction whatever it assumed on the way there. Note also what the page is
*not* saying: Froment's own claim on p. 335 is qualified to *moderately fast*
reactions, and the worst point found above sits inside that band - between his
printed $\gamma > 3$ line and the instantaneous ones - and still errs by less
than the 10 percent he prints. Measured, his qualified sentence stands. What
fails is the unqualified reading of it, that weak depletion is *why* the
approximation works.
"""))

cells.append(code('''# the two wordings, printed from the CSV rather than retyped here
for k in ("froment_assumption_p335", "froment_assumption_p340"):
    print(f'book p.{PRN.loc[k, "book_page"]}  "{quoted(k)}"\\n')
assert "weakly depleted" in quoted("froment_assumption_p335")
assert "reasonably confirmed" in quoted("froment_assumption_p335")
assert "remains approximately constant" in quoted("froment_assumption_p340")

print(f"at q = {Q_STAR:.4f}, the worst curve of the whole chart:")
print(f'{"":26s}{"gamma":>10}{"beta(0)":>12}{"beta_w/beta(0)":>16}'
      f'{"F_A film":>12}{"VKH error":>12}')
rows = [("B half gone", W["gamma_at_half_depletion"], 0.5,
         W["_beta_constancy_at_half"], W["vkh_error_at_half_depletion_pct"]),
        ("the worst point", G_STAR, W["beta_interface_at_worst"],
         W["beta_constancy_at_worst"], W["vkh_worst_rel_error_pct"]),
        ("B 99 % gone", W["gamma_at_99pct_depletion"], 0.01,
         W["beta_constancy_at_99pct_depletion"],
         W["vkh_error_at_99pct_depletion_pct"]),
        ("B gone (axis end)", G_HI, W["_beta_at_full_depletion"],
         None, W["vkh_error_at_full_depletion_pct"])]
for name, g, b, c, e in rows:
    print(f"  {name:24s}{g:10.4f}{b:12.4e}"
          f"{('-' if c is None else f'{c:.4g}'):>16}"
          f"{REF.f_a(g, Q_STAR)[0]:12.5f}{e:+11.5f} %")
print(f"\\n  p.335's test, beta(0), falls monotonically;  p.340's test,"
      f" beta_w/beta(0), rises monotonically;")
print(f"  the error does NEITHER - it peaks at gamma = {G_STAR:.4f} and vanishes"
      f" at both ends of the range.")
print(f"  beta(0) at the axis end is numerically zero: the solve returns"
      f" {W['_beta_at_full_depletion']:.1e}, whose")
print(f"  magnitude is far below anything this discretisation resolves and whose"
      f" sign is round-off, so the")
print(f"  page states the order and not a bound with digits - and the constancy"
      f" ratio, which divides by")
print(f"  it, is left out of that row rather than printed as a number.")
G_MID = 100.0
B_MID = REF.f_a(G_MID, Q_STAR)[1]
C_MID = REF.beta_constancy(G_MID, Q_STAR)
E_MID = err_pct(G_MID, Q_STAR)
print(f"\\n  One decade lower, where beta(0) is still meaningful:  at gamma ="
      f" {G_MID:g}, beta(0) = {B_MID:.3e}")
print(f"  and beta_w/beta(0) = {C_MID:.4g}, so p.340's constancy assumption is out by"
      f" {np.log10(C_MID):.0f} orders of")
print(f"  magnitude - and the error there is {E_MID:+.4f} %,"
      f" {abs(W['vkh_worst_rel_error_pct'] / E_MID):.0f} times SMALLER"
      f" than at the worst point.")
# The constancy ratio is EXTRAPOLATED like everything else, and it has to be:
# a ratio of two integrals of the same profile does NOT cancel the grid once
# beta is a thin exponential.  Both gaps are measured here, not asserted.
C_FINE_STAR = FILM.beta_constancy(G_STAR, Q_STAR)
C_FINE_MID = FILM.beta_constancy(G_MID, Q_STAR)
print(f"\\n  the constancy ratio is Richardson-extrapolated on the same two grids"
      f" as F_A, and needs to be:")
print(f"    at the worst point the fine n_x = {FILM.n_x} grid alone gives"
      f" {C_FINE_STAR:.6g} against {W['beta_constancy_at_worst']:.6g}"
      f"   ({abs(C_FINE_STAR / W['beta_constancy_at_worst'] - 1):.1e})")
print(f"    at gamma = {G_MID:g} it gives"
      f" {C_FINE_MID:.6g} against {C_MID:.6g}"
      f"   ({abs(C_FINE_MID / C_MID - 1):.1e})")
print(f"  so the grid all but cancels where the headline lives and does NOT at the"
      f" top of the range: the")
print(f"  one-grid shortcut this page used to take on this diagnostic costs"
      f" {abs(C_FINE_MID / C_MID - 1) / abs(C_FINE_STAR / W['beta_constancy_at_worst'] - 1):.0f}"
      f" times more at gamma = {G_MID:g}")
print(f"  than at the point it was justified on, which is why it is gone.")
print(f"\\nthe two depletion thresholds are ROOT-FOUND on beta(0), not sampled:"
      f"\\n  beta(0) = 0.5  at gamma = {W['gamma_at_half_depletion']:.6f}"
      f"   (residual {REF.f_a(W['gamma_at_half_depletion'], Q_STAR)[1] - 0.5:+.2e})"
      f"\\n  beta(0) = 0.01 at gamma = {W['gamma_at_99pct_depletion']:.6f}"
      f"   (residual {REF.f_a(W['gamma_at_99pct_depletion'], Q_STAR)[1] - 0.01:+.2e})")

# ---- the error map, and the depletion contours on top of it ----------------
GG = np.geomspace(G_LO, G_HI, 34)
QQ = np.geomspace(Q_LO, Q_HI, 18)
EPS = np.zeros((QQ.size, GG.size))
BET = np.zeros_like(EPS)
for i, q in enumerate(QQ):
    for j, g in enumerate(GG):
        F, b = FILM.f_a(g, q)
        EPS[i, j] = 100.0 * (f_a_vkh(g, q) - F) / F
        BET[i, j] = b

EPS_MAX = float(EPS.max())
VKH_ALWAYS_UNDER = bool(EPS_MAX < 0.0)
print(f"\\nover the whole {EPS.size}-point map the largest error is"
      f" {EPS_MAX:+.3e} %: VKH under-estimates everywhere it was evaluated"
      f" ({VKH_ALWAYS_UNDER})")
print("this one is a THEOREM and the map only illustrates it, so it is a check that"
      "\\ncould only fail by numerical accident:  alpha'' = gamma^2 alpha beta >= 0,"
      " so alpha' rises"
      "\\nfrom alpha'(0) = -F_A and beta' = (alpha' + F_A)/q >= 0, i.e. beta(xi) >="
      " beta(0) everywhere."
      "\\nWrite T(F) for the pseudo-first-order enhancement with beta frozen at"
      " 1 - (F-1)/q, which is"
      "\\nexactly beta(0); (6.3.5-1) says F_VKH = T(F_VKH).  The true film sees"
      " beta >= beta(0), hence"
      "\\nF_film >= T(F_film); T is strictly decreasing, so F - T(F) is strictly"
      " increasing and"
      "\\nF_film >= F_VKH for EVERY (gamma, q).  The map is drawn on the single"
      f" n_x = {FILM.n_x} graded"
      "\\ngrid, which is enough for a sign; every NUMBER reported comes from the"
      " extrapolated reference.")
print(f"the worst point sits at gamma = {G_STAR:.4f}, on the film side of the"
      f" book's own regime line - it prints"
      f' "{quoted("froment_film_completion_gamma")}" (book p.'
      f'{PRN.loc["froment_film_completion_gamma", "book_page"]}) and'
      f' "{quoted("froment_bulk_reaction_gamma")}"')

fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))
cf = axes[0].contourf(GG, QQ, EPS, levels=np.linspace(-3.0, 0.0, 13), cmap="viridis")
for _g, _ls in ((printed("froment_bulk_reaction_gamma"), ":"),
                (printed("froment_film_completion_gamma"), "-.")):
    if G_LO <= _g <= G_HI:
        axes[0].axvline(_g, color="w", lw=1.0, ls=_ls)
axes[0].contour(GG, QQ, BET, levels=[0.01, 0.5], colors=[C_ORANGE, C_YELLOW],
                linewidths=1.6)
axes[0].plot([G_STAR], [Q_STAR], "*", ms=13, color=C_ORANGE, mec="k", mew=0.6)
axes[0].set_xscale("log"); axes[0].set_yscale("log")
axes[0].set_xlabel(r"$\\gamma$"); axes[0].set_ylabel("$q = E_i - 1$")
axes[0].set_title("VKH error [%]; orange $\\\\beta(0)=0.01$, yellow $\\\\beta(0)=0.5$,\\n"
                  "star = the root-found worst point", fontsize=9)
fig.colorbar(cf, ax=axes[0])

for name, g, c in (("B half gone", W["gamma_at_half_depletion"], C_YELLOW),
                   ("worst point", G_STAR, C_ORANGE),
                   ("B 99 % gone", W["gamma_at_99pct_depletion"], C_BLUE)):
    u = FILM.solve(g, Q_STAR)
    axes[1].plot(FILM.x_c, u[:, 0], color=c, lw=2.0,
                 label=rf"$\\alpha$, {name} ($\\gamma$ = {g:.2f})")
    axes[1].plot(FILM.x_c, u[:, 1], color=c, lw=1.4, ls="--")
axes[1].set_xlim(0, 0.6)
axes[1].set_xlabel(r"$\\xi = y/y_L$")
axes[1].set_ylabel("dimensionless concentration")
axes[1].set_title(f"profiles at q = {Q_STAR:.4f}\\n"
                  r"solid $\\alpha$ (gas $A$), dashed $\\beta$ (reactant $B$)",
                  fontsize=9)
axes[1].legend(frameon=False, fontsize=7.5)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""### 4. Levenspiel's three expansions, and what each one expands

Levenspiel prints three "more precisely" corrections in the annotation boxes of
Fig. 23.4. They are an independent transcription of the same result by a
different author in a different notation, so every constant in them checks a
constant in Froment's equation - and the two of them that can be checked
against the *film equations* as well turn out to behave very differently.

- **Small $M_H$:** $E \cong 1 + M_H^2/3$. This is the Taylor series of
  $\gamma/\tanh\gamma$, so the constant 3 checks Froment's (6.3.2-11)
  directly. It is good to 1 % up to a root-found $M_H$ reported below.
- **Pseudo-first-order branch,** declared for $E_i > 5M_H$:
  $E = M_H\left(1 - (M_H-1)/2E_i\right)$. Expanding Froment's (6.3.5-1) in the
  same limit gives the same expression with $E_i - 1$ where Levenspiel prints
  $E_i$ - the difference between the two books' groups, second order in that
  branch. Both are evaluated against the film solve at the boundary of the
  declared domain.
- **Instantaneous branch,** declared for $E_i < M_H/5$:
  $E \cong E_i - E_i^2(E_i-1)/M_H^2$. **This one is the asymptote of the
  approximation, not of the equations,** and both halves of that statement are
  exact rather than empirical. Put $F_A = 1 + q - d$ into (6.3.5-1): for large
  $\gamma$ it gives $\gamma\sqrt{d/q} \to 1 + q$, i.e.
  $d \to q(1+q)^2/\gamma^2$, which in Levenspiel's symbols *is*
  $E_i^2(E_i-1)/M_H^2$. The film's own deficit, on the other hand, is fixed by
  the first integral of the two film equations, $q\,\beta(0) = 1 + q - F_A$: it
  equals $q\,\beta(0)$ identically, and $\beta(0)$ is the interfacial value of
  a reaction layer whose thickness goes as $1/\gamma$, so it is exponentially
  small and can follow no power law at all. The measured local exponent
  steepens over successive octaves, and the printed correction runs further and
  further above the truth.
"""))

cells.append(code('''DEFICIT_M = (50.0, 100.0, 200.0, 400.0, 800.0)


def block_lev(film=None, vkh=None, small_coeff=None, pfo_denom=None,
              pfo_offset=None, inst_power=None, deficit_film=None,
              bvp_tol=1e-10, **kw):
    """Levenspiel's three printed expansions, tested against both routes."""
    film = REF if film is None else film
    vkh = f_a_vkh if vkh is None else vkh
    c3 = printed("lev_small_M_coeff") if small_coeff is None else small_coeff
    c2 = printed("lev_pfo_correction_denom") if pfo_denom is None else pfo_denom
    c1 = printed("lev_pfo_correction_offset") if pfo_offset is None else pfo_offset
    cp = printed("lev_inst_correction_power") if inst_power is None else inst_power
    out = {}

    # (a) the small-M_H expansion, 1 % threshold ROOT-FOUND
    q_a = 4.0
    out["lev_small_M_1pct_threshold"] = root_log(
        lambda m: abs((1.0 + m**2 / c3) / film.f_a(m, q_a, **kw)[0] - 1) - 0.01,
        0.1, 3.0)

    # (b) the pseudo-first-order branch, AT the boundary E_i = 5 M_H it declares
    m_b = 10.0
    Ei_b = printed("lev_pfo_threshold") * m_b
    q_b = Ei_b - 1.0
    F_b = film.f_a(m_b, q_b, **kw)[0]
    out["lev_pfo_expansion_err_pct"] = 100.0 * (
        m_b * (1.0 - (m_b - c1) / (c2 * Ei_b)) / F_b - 1.0)
    out["vkh_form_pfo_expansion_err_pct"] = 100.0 * (
        m_b * (1.0 - (m_b - c1) / (c2 * q_b)) / F_b - 1.0)
    out["uncorrected_E_eq_M_err_pct"] = 100.0 * (m_b / F_b - 1.0)

    # (c) the instantaneous branch
    Ei_c = 20.0
    q_c = Ei_c - 1.0
    fd = film if deficit_film is None else deficit_film

    def pred(m):
        return Ei_c**cp * (Ei_c - 1.0) / m**2

    def dnum(m):
        return Ei_c - fd.f_a(m, q_c, **kw)[0]

    out["lev_inst_asymptote_ratio_vkh"] = (Ei_c - vkh(400.0, q_c)) / pred(400.0)
    out["lev_inst_asymptote_ratio_num"] = dnum(400.0) / pred(400.0)
    out["lev_inst_overstatement_factor"] = pred(800.0) / dnum(800.0)
    dn = [dnum(m) for m in DEFICIT_M]
    dv = [Ei_c - vkh(m, q_c) for m in DEFICIT_M]
    # abs() only so that a break row which drives a deficit negative still
    # returns a number instead of a nan; every baseline deficit is positive.
    ex = [float(np.log(abs(a[i + 1] / a[i]))
                / np.log(DEFICIT_M[i + 1] / DEFICIT_M[i]))
          for a in (dn, dv) for i in range(len(DEFICIT_M) - 1)]
    n_gap = len(DEFICIT_M) - 1
    out["inst_deficit_exponent_num_steepest"] = min(ex[:n_gap])
    out["inst_deficit_exponent_vkh_finest"] = ex[2 * n_gap - 1]
    out["_exponents_num"], out["_exponents_vkh"] = ex[:n_gap], ex[n_gap:]
    out["_deficit_num"], out["_deficit_vkh"] = dn, dv
    # TWO INDEPENDENT METHODS, one of them also refined - not three routes.
    # pymrm finite volumes and scipy's collocation share no assembly; the
    # middle entry is the same assembly on a finer, differently graded mesh, so
    # it prices the discretisation of the first and not its formulation.
    alt = graded_film(2 * fd.n_x, fd.stretch + 1.0)
    three = [dnum(800.0), Ei_c - alt.f_a(800.0, q_c, **kw)[0],
             Ei_c - f_a_bvp(800.0, q_c, tol=bvp_tol)[0]]
    out["inst_deficit_routes_max_rel_diff"] = max(three) / min(three) - 1.0
    out["_three"] = three
    # the EXACT reason the film deficit has no power law: the first integral of
    # the two film equations makes it q*beta(0) identically, and beta(0) at the
    # interface is exponentially small in gamma.
    out["_beta0_at_deficit"] = fd.f_a(800.0, q_c, **kw)[1]
    return out


L = block_lev()
print(f"(a) 1 + M_H^2/{printed('lev_small_M_coeff'):g} stays within 1 % of the film"
      f" solve up to M_H = {L['lev_small_M_1pct_threshold']:.4f}  (root-found)")
_m, _Ei = 10.0, printed("lev_pfo_threshold") * 10.0
print(f"\\n(b) at the declared boundary E_i = {printed('lev_pfo_threshold'):g} M_H,"
      f" with M_H = {_m:g} and E_i = {_Ei:g}:")
print(f"      E = M_H                              {L['uncorrected_E_eq_M_err_pct']:+8.4f} %")
print(f"      Levenspiel, printed 2 E_i            {L['lev_pfo_expansion_err_pct']:+8.4f} %")
print(f"      same expansion of (6.3.5-1), 2(E_i-1){L['vkh_form_pfo_expansion_err_pct']:+8.4f} %")
print(f"    -> the correction is worth"
      f" {abs(L['uncorrected_E_eq_M_err_pct'] / L['lev_pfo_expansion_err_pct']):.2f}x,"
      f" and the two books' groups differ by"
      f" {abs(L['lev_pfo_expansion_err_pct'] - L['vkh_form_pfo_expansion_err_pct']):.4f}"
      f" percentage points there")
print(f"\\n(c) instantaneous branch, E_i = 20:  deficit E_i - E against the printed"
      f" E_i^2(E_i-1)/M_H^2")
print(f'{"M_H":>8}{"film solve":>14}{"VKH":>14}{"printed":>14}{"printed/film":>14}')
for m, dn, dv in zip(DEFICIT_M, L["_deficit_num"], L["_deficit_vkh"]):
    p = 20.0**2 * 19.0 / m**2
    print(f"{m:8.0f}{dn:14.4e}{dv:14.4e}{p:14.4e}{p / dn:14.2f}")
print(f"  local exponent d ln(deficit)/d ln M_H, octave by octave")
print(f"    film solve : {[round(e, 4) for e in L['_exponents_num']]}")
print(f"    VKH        : {[round(e, 4) for e in L['_exponents_vkh']]}")
print(f"  the VKH deficit converges on the printed -2; the film solve's does not"
      f"\\n  the M_H = 800 deficit by TWO INDEPENDENT METHODS, the first also refined:"
      f"\\n    pymrm FV, extrapolated (n_x = {FILM_C.n_x}, {FILM.n_x}, s = 8)"
      f"  {L['_three'][0]:.6e}"
      f"\\n    pymrm FV, n_x = {2 * FILM.n_x}, s = 9 (the same assembly, refined)  "
      f"{L['_three'][1]:.6e}"
      f"\\n    scipy collocation (no assembly in common)               "
      f"{L['_three'][2]:.6e}"
      f"\\n  max relative spread {L['inst_deficit_routes_max_rel_diff']:.3e}"
      f" - far below the {L['lev_inst_overstatement_factor']:.1f}x it is being used to"
      f" establish")
print(f"\\n  WHY there can be no power law, exactly.  Subtracting the two film"
      f" equations gives"
      f"\\n  (beta - alpha/q)'' = 0, and the four boundary conditions close it:"
      f" q beta(0) = 1 + q - F_A"
      f"\\n  identically, so the film deficit E_i - E IS q beta(0) - and beta(0) is"
      f" the interfacial"
      f"\\n  value of a WKB-thin reaction layer, exponentially small in M_H."
      f"  At M_H = 800 the film"
      f"\\n  solve gives beta(0) = {L['_beta0_at_deficit']:.4e}, so"
      f" q beta(0) = {19.0 * L['_beta0_at_deficit']:.4e} against the deficit"
      f" {L['_three'][0]:.4e}.")
print(f"  The VKH deficit, by contrast, is a power law exactly: putting"
      f" F_A = 1 + q - d in (6.3.5-1)"
      f"\\n  gives gamma sqrt(d/q) -> 1 + q, i.e. d -> q(1+q)^2/gamma^2, which in"
      f" Levenspiel's symbols"
      f"\\n  IS E_i^2(E_i - 1)/M_H^2.  His 'more precisely' correction is the"
      f" ALGEBRAIC ASYMPTOTE OF THE"
      f"\\n  APPROXIMATION, not of the equations - which is why it tracks VKH to"
      f" {L['lev_inst_asymptote_ratio_vkh']:.4f} and the film"
      f"\\n  solve to only {L['lev_inst_asymptote_ratio_num']:.4f}.")
print(f"  And 44x too big is 44x of a small correction: at M_H = 800 the"
      f" uncorrected E = E_i is already"
      f"\\n  within {100 * L['_three'][0] / 20.0:.5f} % of the film solve, so"
      f" Levenspiel's correction is there WORSE than"
      f"\\n  leaving it out - it moves E away from the truth by"
      f" {100 * (20.0**2 * 19.0 / 800.0**2 - L['_three'][0]) / 20.0:.5f} percentage"
      f" points.")

fig, ax = plt.subplots(figsize=(6.0, 4.2))
MM = np.geomspace(30.0, 1000.0, 30)
ax.loglog(MM, [20.0 - f_a_vkh(m, 19.0) for m in MM], color=C_BLUE, lw=2.0,
          label="VKH (6.3.5-1)")
ax.loglog(MM, [20.0 - REF.f_a(m, 19.0)[0] for m in MM], color=C_ORANGE, lw=2.0,
          label="film solve (pymrm)")
ax.loglog(MM, 20.0**2 * 19.0 / MM**2, color=C_GREY, ls="--", lw=1.4,
          label=r"printed $E_i^2(E_i-1)/M_H^2$")
ax.set_xlabel("$M_H$"); ax.set_ylabel("$E_i - E$")
ax.set_title("What Levenspiel's instantaneous correction expands\\n"
             "($E_i$ = 20; the printed curve tracks the approximation)",
             fontsize=9.5)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
plt.show()'''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Six checks, none of which uses fitted data because there is none. Two of them
are about the **baseline itself**, because the whole page is a set of errors
quoted against it: the extrapolated reference is compared with an independent
solver family and with a finer extrapolation of itself, and the grid study is
run on the graded *and* the uniform grid so that what grading buys and what it
costs are both on the page.
"""))

cells.append(code('''def block_routes(film=None, ref=None, vkh=None, stretch=None, newton_tol=1e-10,
                 bvp_tol=1e-10, bc_beta_bulk="dirichlet", n_grid=(100, 200, 400, 800),
                 **kw):
    """Two solvers for the film, two for the implicit relation, and the grid."""
    film = FILM if film is None else film
    ref = REF if ref is None else ref
    vkh = f_a_vkh if vkh is None else vkh
    s = film.stretch if stretch is None else stretch
    pts = [(g, q) for q in (1.0, 4.0, 19.0, 99.0)
           for g in (1.0, 3.0, 10.0, 30.0, 100.0)]
    col = {p: f_a_bvp(*p, tol=bvp_tol) for p in pts}
    out = {
        "fv_vs_bvp_max_rel_diff": float(max(
            abs(film.f_a(g, q, tol=newton_tol, **kw)[0] / col[(g, q)][0] - 1.0)
            for g, q in pts)),
        "fv_vs_bvp_max_abs_diff_beta": float(max(
            abs(film.f_a(g, q, tol=newton_tol, **kw)[1] - col[(g, q)][1])
            for g, q in pts)),
        # the SAME comparison for the baseline this page actually reports
        "ref_vs_bvp_max_rel_diff": float(max(
            abs(ref.f_a(g, q, tol=newton_tol, **kw)[0] / col[(g, q)][0] - 1.0)
            for g, q in pts)),
        # and the baseline against a finer extrapolation of itself: the
        # refinement evidence for the reference, not for one grid of it
        "ref_self_refinement_max_rel_diff": float(max(
            abs(ref.f_a(g, q, tol=newton_tol, **kw)[0]
                / REF_FINE.f_a(g, q, tol=newton_tol, **kw)[0] - 1.0)
            for g, q in pts)),
        # the first integral of the two film equations, which no discretisation
        # is told about: q beta(0) = 1 + q - F_A.  Normalised by 1 + q because
        # at large gamma both sides are round-off.
        "first_integral_max_abs_resid": float(max(
            abs(q * ref.f_a(g, q, tol=newton_tol, **kw)[1]
                - (1.0 + q - ref.f_a(g, q, tol=newton_tol, **kw)[0])) / (1.0 + q)
            for g, q in pts)),
        "vkh_two_routes_max_rel_diff": float(max(
            abs(vkh(g, q) / f_a_vkh_route_b(g, q) - 1.0) for g, q in pts)),
        "vkh_pfo_limit_max_rel_error": float(max(
            abs(vkh(g, 1e10) / (g / np.tanh(g)) - 1.0)
            for g in np.geomspace(G_LO, 300.0, 40))),
        "vkh_instantaneous_limit_rel_error": float(abs(vkh(1e6, 19.0) / 20.0 - 1.0)),
    }
    g_ref, q_ref = 10.0, 4.0
    exact = f_a_bvp(g_ref, q_ref, tol=1e-11)[0]

    def errs_on(st):
        return [abs(Film(n_x=n, stretch=st, bc_beta_bulk=bc_beta_bulk)
                    .f_a(g_ref, q_ref, tol=newton_tol, **kw)[0] / exact - 1.0)
                for n in n_grid]

    errs, errs_u = errs_on(s), errs_on(0.0)
    out["fv_convergence_order"] = float(np.log2(errs[-2] / errs[-1]))
    out["fv_rel_error_n800"] = float(errs[-1])
    # WHAT GRADING COSTS, at the same n and against the same reference
    out["fv_rel_error_n800_uniform"] = float(errs_u[-1])
    out["_grid_errs"], out["_grid_errs_uniform"] = errs, errs_u
    return out


R = block_routes()
print("1. the film, two solver families, over 20 points of the chart")
print(f"   max |pymrm n_x = {FILM.n_x} / collocation - 1| on F_A :"
      f" {R['fv_vs_bvp_max_rel_diff']:.3e}")
print(f"   max |difference| on beta(0)                : "
      f"{R['fv_vs_bvp_max_abs_diff_beta']:.3e}")
print(f"\\n2. THE BASELINE ITSELF - Richardson(n_x = {FILM_C.n_x}, {FILM.n_x}),"
      f" which is what every")
print(f"   number on this page is quoted against")
print(f"   max |reference / collocation - 1|          : "
      f"{R['ref_vs_bvp_max_rel_diff']:.3e}")
print(f"   max |reference / a finer extrapolation - 1|: "
      f"{R['ref_self_refinement_max_rel_diff']:.3e}"
      f"   (Richardson {FILM.n_x}, {2 * FILM.n_x})")
print(f"   -> the reference is converged {R['fv_vs_bvp_max_rel_diff'] / R['ref_vs_bvp_max_rel_diff']:.0f}x"
      f" better than the single n_x = {FILM.n_x} solve it")
print(f"      is built from, which is the whole reason it exists")
print(f"\\n3. the first integral q beta(0) = 1 + q - F_A, which the discretisation"
      f" is never told")
print(f"   max |q beta(0) - (1 + q - F_A)| / (1 + q)  : "
      f"{R['first_integral_max_abs_resid']:.3e}")
print("\\n4. grid refinement at (gamma, q) = (10, 4), against the collocation solve")
print(f'{"":6s}{"n_x":>7}{"graded s = 8":>16}{"uniform":>14}{"graded/uniform":>17}')
for n, e, eu in zip((100, 200, 400, 800),
                    R["_grid_errs"], R["_grid_errs_uniform"]):
    print(f'{"":6s}{n:7d}{e:16.4e}{eu:14.4e}{e / eu:17.1f}')
print(f"   observed order {R['fv_convergence_order']:.4f}   (expected 2)")
print(f"   AT THIS POINT GRADING IS A LIABILITY: at the same n_x the graded grid is"
      f"\\n   {R['fv_rel_error_n800'] / R['fv_rel_error_n800_uniform']:.0f}x LESS"
      f" accurate than a uniform one, because its outer cells are"
      f"\\n   {(1 - FILM.x_f[-2]) / (1.0 / FILM.n_x):.1f}x wider."
      f"  Grading earns its place at the other end of the chart"
      f"\\n   (see the pseudo-first-order check below), and the extrapolated"
      f" reference is how"
      f"\\n   this page has both.")
print("\\n5. the implicit relation, two unknowns")
print(f"   max |solve-for-F_A / solve-for-gamma' - 1| : "
      f"{R['vkh_two_routes_max_rel_diff']:.3e}")
print("\\n6. the two limits (6.3.5-1) must reproduce")
print(f"   q -> inf gives gamma/tanh gamma, max rel : "
      f"{R['vkh_pfo_limit_max_rel_error']:.3e}")
print(f"   gamma -> inf gives 1 + q, rel            : "
      f"{R['vkh_instantaneous_limit_rel_error']:.3e}")'''))

cells.append(md(r"""### 5. Reconciling with `F3.1`, which owns Hatta

`F3.1` reports a maximum Van Krevelen-Hoftijzer error of about 2.1 %. That is a
different number from this page's headline and it should be: `F3.1` sweeps
three fixed $E_i$ values on a fixed $\mathrm{Ha}$ grid and takes the largest
sampled value, while this page lets $q$ vary continuously and root-finds the
extremum. How much of the gap is *root-finding* and how much is *letting $q$
move* is measured below rather than asserted.

The reconciliation is done by **loading `F3.1`'s own `agreement.json` and
recomputing its definition here**, not by retyping its number. The residual
that survives is then *attributed by measurement* rather than left to be
guessed at, and it belongs to `F3.1`, not to this page: running the same
recomputation **cold on `F3.1`'s own grid** - uniform, $n_x = 400$ - lands on
`F3.1`'s stored number, so the residual is that grid's discretisation error and
the warm start contributes nothing to it. That cold run is also what settles
the half of the sentence a reader would otherwise have to take on trust: this
is **not** a warm-start-versus-cold-start difference, which is the first thing
`F3.1`'s continued sweep suggests. This page's extrapolated reference is the
accurate side of the residual, and the single graded grid this page used before
the baseline was extrapolated is on the wrong side of it by two orders of
magnitude more.
"""))

cells.append(code('''def load_agreement(page):
    for _p in (Path.cwd(), *Path.cwd().parents):
        c = _p / "pages" / page / "agreement.json"
        if c.is_file():
            return json.loads(c.read_text())
    with urllib.request.urlopen(f"{RAW}/pages/{page}/agreement.json") as fh:
        return json.load(fh)


def block_f31(film=None, vkh=None, n_scan=13, **kw):
    """F3.1's OWN definitions - its three E_i, its two Ha grids - recomputed here."""
    film = REF if film is None else film
    ha = np.geomspace(0.05, 300.0, 55)
    v = float(max(abs(err_pct(h, e - 1.0, film=film, vkh=vkh, **kw))
                  for e in (5.0, 20.0, 100.0) for h in ha))
    out = {"f31_definition_vkh_max_error_pct": v}
    # the same three E_i, but with gamma ROOT-FOUND instead of sampled: this
    # splits the gap to this page's headline into "root-finding" and "letting q
    # move", which are not the same size at all.
    out["f31_definition_rootfound_pct"] = float(max(
        abs(argext_log(lambda g: err_pct(g, e - 1.0, film=film, vkh=vkh, **kw),
                       G_LO, G_HI, n_scan)[1])
        for e in (5.0, 20.0, 100.0)))
    # F3.1's other overlapping check: the pseudo-first-order solve against
    # Ha/tanh(Ha) over ITS Ha range.  Run on THREE single grids, because the
    # 893.9x this page used to credit to grading is really grading TIMES a
    # refinement: F3.1 uses n_x = 400 uniform and this page n_x = 800 graded.
    ha1 = np.geomspace(0.02, 200.0, 70)

    def pfo(f):
        return float(max(abs(f.f_a(h, 1.0, **{**kw, "frozen_b": True})[0]
                             / (h / np.tanh(h)) - 1.0) for h in ha1))

    n_half = max(film.n_x // 2, 8)
    out["pfo_max_rel_error_graded"] = pfo(graded_film(film.n_x, film.stretch))
    out["pfo_max_rel_error_uniform_same_n"] = pfo(graded_film(film.n_x, 0.0))
    out["pfo_max_rel_error_uniform_half_n"] = pfo(graded_film(n_half, 0.0))
    if F31 is not None:
        out["f31_stored_vs_recomputed_rel"] = abs(
            v / F31["metrics"]["vkh_max_error_pct"] - 1.0)
        out["f31_first_order_error_ratio"] = (
            F31["metrics"]["first_order_max_rel_error"]
            / out["pfo_max_rel_error_graded"])
    return out


try:
    F31 = load_agreement("F3.1-hatta-regimes")
except Exception as exc:                      # no checkout and no network
    F31 = None
    print(f"F3.1's agreement.json is not reachable ({type(exc).__name__});"
          " the reconciliation below is skipped rather than typed from memory.")

F = block_f31()
# the same block on the single graded grid: what this page reported before the
# baseline was extrapolated, kept so that shortcut's cost stays visible
F_PLAIN = block_f31(film=FILM)
# and the same block on F3.1's OWN grid - uniform, n_x = 400 - which is what
# ATTRIBUTES the residual instead of hand-waving it.  A COLD run on that grid
# reproduces F3.1's stored number, so the residual is neither the warm start
# nor this page's numerics: it is F3.1's discretisation.
FILM_F31 = graded_film(400, 0.0)
F_F31GRID = block_f31(film=FILM_F31)
print("what F3.1 states about the rows this page overlaps with, and whether it bites")
print("  * its VKH error is a max over E_i in {5, 20, 100} on 55 Ha points, warm")
print("    started along each E_i - so it is a SAMPLED max on a continued sweep.")
print("    Affects this page: the definitions differ, so the numbers must differ.")
print("  * its first-order solve matches Ha/tanh(Ha) over 0.02 <= Ha <= 200 at")
print("    n_x = 400 UNIFORM, with the error 'concentrated at high Ha where the")
print("    reaction layer is thinner than a cell'.  Affects this page directly:")
print("    it is the reason this page grades its grid AT ALL, and the same check")
print("    is run here on three grids so that grading and refinement are separated")
print("    instead of being credited to grading together.")
if F31 is not None:
    print(f"\\n  F3.1 agreement.json  first_order_max_rel_error ="
          f" {F31['metrics']['first_order_max_rel_error']:.4e}  (uniform, n_x = 400)")
    print(f"  the SAME check run here on a uniform n_x = {max(FILM.n_x // 2, 8)} grid ="
          f" {F['pfo_max_rel_error_uniform_half_n']:.4e}"
          f"  -> F3.1 reproduced to"
          f" {abs(F['pfo_max_rel_error_uniform_half_n'] / F31['metrics']['first_order_max_rel_error'] - 1):.1e}")
    print(f"  DECOMPOSING THE 893.9x THIS PAGE USED TO CREDIT TO GRADING ALONE:")
    print(f"    uniform n_x = {max(FILM.n_x // 2, 8)} -> uniform n_x = {FILM.n_x}"
          f"   (refinement, no grading) "
          f" {F['pfo_max_rel_error_uniform_half_n'] / F['pfo_max_rel_error_uniform_same_n']:8.2f}x")
    print(f"    uniform n_x = {FILM.n_x} -> graded n_x = {FILM.n_x}"
          f"     (grading, same n_x)   "
          f" {F['pfo_max_rel_error_uniform_same_n'] / F['pfo_max_rel_error_graded']:8.2f}x")
    print(f"    product                                              "
          f" {F['f31_first_order_error_ratio']:8.1f}x")
    print(f"  SO GRADING IS WORTH"
          f" {F['pfo_max_rel_error_uniform_same_n'] / F['pfo_max_rel_error_graded']:.0f}x"
          f" ON THIS CHECK, NOT {F['f31_first_order_error_ratio']:.0f}x - and it is"
          f" worth it here because")
    print(f"  the check runs to Ha = 200, where the reaction layer is thinner than a"
          f" uniform cell.")
    print(f"  At the moderate gamma this page's headline lives at, the same grading"
          f" COSTS a factor"
          f" {R['fv_rel_error_n800'] / R['fv_rel_error_n800_uniform']:.0f} (Validation,"
          f" check 4), which is why the")
    print(f"  baseline is an extrapolation of two graded grids rather than either grid"
          f" on its own.")
if F31 is not None:
    print(f"\\nF3.1 agreement.json  vkh_max_error_pct  = "
          f"{F31['metrics']['vkh_max_error_pct']:.6f}")
    print(f"recomputed here, same definition, cold  = "
          f"{F['f31_definition_vkh_max_error_pct']:.6f}"
          f"   (relative {F['f31_stored_vs_recomputed_rel']:.2e})")
    print(f"  the SAME recomputation, run COLD on F3.1's OWN grid (uniform,"
          f" n_x = {FILM_F31.n_x}) = "
          f"{F_F31GRID['f31_definition_vkh_max_error_pct']:.9f}")
    print(f"      -> reproduces F3.1's stored value to"
          f" {F_F31GRID['f31_stored_vs_recomputed_rel']:.1e} relative.")
    print(f"  -> SO THE RESIDUAL IS F3.1's OWN GRID ERROR, not this page's numerics"
          f" - and the cold run on")
    print(f"     F3.1's grid PROVES the other half too: the warm start is worth"
          f" nothing here, since cold on")
    print(f"     that grid lands on the warm-started stored number.  What is left,"
          f" {F['f31_stored_vs_recomputed_rel']:.2e}, is the")
    print(f"     uniform 400-cell discretisation; this page's extrapolated reference"
          f" is the accurate side of")
    print(f"     it, agreeing with scipy's collocation to"
          f" {R['ref_vs_bvp_max_rel_diff']:.2e} over 20 points of the chart.")
    print(f"     (For scale: this page's OWN pre-extrapolation shortcut, a single"
          f" n_x = {FILM.n_x} graded grid,")
    print(f"      lands {F_PLAIN['f31_definition_vkh_max_error_pct'] - F31['metrics']['vkh_max_error_pct']:+.6f}"
          f" away - {abs(F_PLAIN['f31_stored_vs_recomputed_rel'] / F['f31_stored_vs_recomputed_rel']):.0f}"
          f" times the residual left here.)")
    print(f"\\nthis page's root-found worst over the chart = "
          f"{abs(W['vkh_worst_rel_error_pct']):.6f}")
    print(f"  the same three E_i, gamma root-found instead of sampled ="
          f" {F['f31_definition_rootfound_pct']:.6f}")
    print(f"  -> of the"
          f" {abs(W['vkh_worst_rel_error_pct']) - F['f31_definition_vkh_max_error_pct']:.4f}"
          f" extra points in this page's headline,"
          f" {F['f31_definition_rootfound_pct'] - F['f31_definition_vkh_max_error_pct']:.4f}"
          f" come from root-finding gamma")
    print(f"     and"
          f" {abs(W['vkh_worst_rel_error_pct']) - F['f31_definition_rootfound_pct']:.4f}"
          f" from letting q move off {{4, 19, 99}}.  It is almost all q.")'''))

# --------------------------------------------------------------- break table
cells.append(md(r"""### The break table

Every metric this page reports needs something that moves it. Each row below
**recomputes** the metrics it claims to cover, under one injected defect, using
the same functions; the coverage map printed underneath is generated from those
measured moves and is never a hand-written claim. A static guard parses each
row's own source and rejects any metric returned as a typed constant - and the
cell below states what that guard cannot do, because a guard trusted past its
reach is worse than none.

Two rows are worth reading before the rest. **"The reference taken as one
graded solve"** is the defect this page actually shipped and a verification
caught; it is injected verbatim rather than quietly fixed. And **"uniform grid
instead of graded"** is *not* filed as a one-way defect any more: at high
$\gamma$ it is one, and on the metrics that carry the headline it agrees with
the extrapolated reference to within the move tolerance, which is why it stops
appearing as cover for *those*. It still buys sub-CI links on the other
moderate-$\gamma$ metrics and large ones at high $\gamma$; the split is printed
under the coverage map rather than described.
"""))

cells.append(code('''BREAKS, BREAK_FNS = [], []
MOVE_TOL = 1e-6


def brk(label, fn, note=""):
    try:
        got = fn()
    except Exception as exc:
        got = {"_raised": f"{type(exc).__name__}: {exc}"}
    BREAKS.append((label, got, note))
    BREAK_FNS.append((label, fn))


def _vkh_sign():
    v = lambda g, q: f_a_vkh(g, q, sign=+1.0)          # noqa: E731
    return {**block_worst(vkh=v), **block_lev(vkh=v), **block_f31(vkh=v),
            **{k: val for k, val in block_routes(vkh=v).items()
               if k.startswith("vkh_")}}


brk("(6.3.5-1) transcribed with + where it prints -", _vkh_sign,
    "the sign under the square root is the whole content of the correction:"
    " with it flipped gamma' exceeds gamma and the approximation over-shoots")


def _vkh_offset():
    v = lambda g, q: f_a_vkh(g, q, offset=0.0)         # noqa: E731
    return {**block_worst(vkh=v), **block_lev(vkh=v), **block_f31(vkh=v),
            **{k: val for k, val in block_routes(vkh=v).items()
               if k.startswith("vkh_")}}


brk("(6.3.5-1) transcribed as F_A where it prints (F_A - 1)", _vkh_offset,
    "the -1 is what makes gamma' -> gamma when there is no enhancement; without"
    " it the pseudo-first-order limit is wrong too, which is why the limit"
    " checks are metrics and not decoration")


def _vkh_denom():
    v = lambda g, q: f_a_vkh(g, q, denom="ei")        # noqa: E731
    return {**block_worst(vkh=v), **block_lev(vkh=v), **block_f31(vkh=v),
            **{k: val for k, val in block_routes(vkh=v).items()
               if k.startswith("vkh_")}}


brk("q read as E_i in (6.3.5-1) - the two books' groups swapped", _vkh_denom,
    "exactly the mistake the two charts invite, since they label the same"
    " curves with q and with E_i; it is worth more than the whole error the"
    " approximation makes")


def _frozen_b():
    return {**block_worst(frozen_b=True), **block_lev(frozen_b=True),
            **block_f31(frozen_b=True), **block_routes(frozen_b=True)}


brk("the depletion of B dropped from the FILM model (beta = 1 identically)",
    _frozen_b,
    "the row that shows every error metric on this page is a MEASUREMENT OF"
    " DEPLETION: take the depletion out of the film and the film becomes"
    " pseudo-first-order, while (6.3.5-1) goes on correcting for a depletion"
    " that no longer happens, so the reported error changes sign and size and"
    " beta(0) becomes 1")


def _coarse():
    f = Film(n_x=100, stretch=8.0)
    return {**block_worst(film=f), **block_lev(film=f, deficit_film=f),
            **block_f31(film=f),
            **block_routes(film=f, n_grid=(25, 50, 100, 200))}


brk("film solved on 100 cells instead of 800", _coarse,
    "the baseline against which a 2.8 % approximation error is quoted has to be"
    " better than 2.8 %; this prices how much better")


def _no_extrapolation():
    """The baseline as a SINGLE graded solve - what this page reported before."""
    return {**block_worst(film=FILM), **block_lev(film=FILM, deficit_film=FILM),
            **block_f31(film=FILM), **block_routes(ref=FILM)}


brk(f"the reference taken as one graded solve at n_x = {FILM.n_x}"
    " instead of extrapolated", _no_extrapolation,
    "THE ROW THIS PAGE EARNED THE HARD WAY.  A single graded solve at n_x = 800"
    " is grid-limited exactly where the headline lives: it moves the worst error"
    " into its fourth decimal, in the flattering direction, and it moves the"
    " margin against the printed 10 percent, the location of the extremum, the"
    " two depletion-curve errors, both Levenspiel expansion errors and the"
    " recomputation of F3.1's definition.  Its error in F_A is smaller than this"
    " page's own pymrm-versus-collocation agreement - the two are printed side"
    " by side in the Validation section - which is exactly why it stayed"
    " invisible until the baseline was extrapolated")


def _uniform():
    f = Film(n_x=800, stretch=0.0)
    return {**block_worst(film=f), **block_lev(film=f, deficit_film=f),
            **block_f31(film=f), **block_routes(film=f, ref=f, stretch=0.0)}


brk("uniform grid instead of graded, at the same n_x", _uniform,
    "NOT A ONE-WAY DEFECT, and the page says so rather than filing it as one."
    " At large gamma the reaction layer is thinner than a uniform cell - the"
    " defect F3.1 names in its own first-order check - and this row prices it:"
    " the pseudo-first-order error to Ha = 200 is two orders worse without"
    " grading.  But at the MODERATE gamma where this page's headline sits the"
    " uniform grid is the better one, and on the metrics that carry the"
    " headline - the worst error itself, the margin against the printed 10"
    " percent, the two depletion thresholds and the worst over the printed"
    " labels - it agrees with the extrapolated reference to within the move"
    " tolerance below, so it does not appear as cover for THOSE.  It is not"
    " cover-free, and an earlier wording of this note wrongly implied it was:"
    " on the other moderate-gamma metrics it sits just the other side of that"
    " tolerance and still buys links, and on the high-gamma metrics it buys"
    " large ones, which is the reason to keep it.  The split is printed under"
    " the coverage map.  The agreement is itself evidence:"
    " two different discretisations, one of them extrapolated, landing on the"
    " same fourth decimal")


def _bc():
    f = Film(n_x=800, stretch=8.0, bc_beta_bulk="neumann")
    return {**block_worst(film=f), **block_lev(film=f, deficit_film=f),
            **block_routes(film=f, bc_beta_bulk="neumann")}


brk("beta given a zero-gradient bulk edge instead of beta(1) = 1", _bc,
    "the bulk value of B is what q is defined against, so this is not a small"
    " perturbation: it removes the supply of B entirely")


def _sampled():
    return block_worst(rootfind=False)


brk("worst point SAMPLED on the scan grid instead of root-found", _sampled,
    "the reason every extremum on this page is a root: the sampled maximum is"
    " grid-limited, and both its location and its value move")


def _coarse_scan():
    return block_worst(n_scan=7)


brk("the bracketing scan halved to 7 points", _coarse_scan,
    "a root-find is only as good as its bracket; this shows the answer does not"
    " depend on the scan once the bracket holds")


def _lev_coeff():
    return block_lev(small_coeff=2.0)


brk("Levenspiel's 1 + M_H^2/3 transcribed with a 2", _lev_coeff,
    "the 3 is the Taylor coefficient of gamma/tanh gamma, so this row checks a"
    " constant of FROMENT's equation through LEVENSPIEL's transcription")


def _lev_denom():
    return block_lev(pfo_denom=1.0)


brk("Levenspiel's 2 E_i transcribed as E_i", _lev_denom, "")


def _lev_offset():
    return block_lev(pfo_offset=0.0)


brk("Levenspiel's (M_H - 1) transcribed as M_H", _lev_offset,
    "the -1 is worth little at M_H = 10 and everything at M_H = 1, which is why"
    " the expansion is evaluated at the boundary the book declares for it")


def _lev_power():
    return block_lev(inst_power=1.0)


brk("Levenspiel's E_i^2(E_i - 1) transcribed as E_i(E_i - 1)", _lev_power,
    "the exponent is what the instantaneous asymptote ratio tests; without this"
    " row that ratio would be an unchecked 0.995")


def _labels_as_ei():
    """Froment's labels read as E_i: the implied plateau becomes L, not 1 + L."""
    lab = [max(L - 1.0, 1.0) for L in LABELS]
    return {"chart_label_max_plateau_ratio": float(max(float(L) / L for L in LABELS)),
            **block_worst(labels=lab)}


brk("Froment's curve labels read as E_i, the way Levenspiel's are", _labels_as_ei,
    "the notation trap, priced: read both books' labels the same way and the"
    " plateau ratio collapses to 1, which is exactly the mistake the metric"
    " exists to record.  DECLARED LIMIT OF THIS ROW: max_L L/L is identically 1"
    " for ANY label set, so for that one metric this row is an identity and not"
    " a measurement; only the row below, which drops a transcribed label,"
    " actually tests the CSV.  chart_label_max_plateau_ratio is itself"
    " definitional - max_L (1+L)/L is 2 whenever the smallest printed label is 1"
    " - and is reported as a property of the two notations, not as a result")


def _labels_truncated():
    """The CSV's smallest curve label dropped."""
    lab = sorted(LABELS)[1:]
    return {"chart_label_max_plateau_ratio":
                float(max((1.0 + L) / L for L in lab)),
            **block_worst(labels=lab)}


brk("the smallest printed curve label dropped from the CSV", _labels_truncated,
    "the label-ratio metric is a property of the transcribed LABEL SET as much"
    " as of the two notations, and this is what a dropped row costs it")


def _vkh_loose():
    v = lambda g, q: f_a_vkh(g, q, xtol=1e-4)          # noqa: E731
    return {k: val for k, val in block_routes(vkh=v).items() if k.startswith("vkh_")}


brk("the implicit relation solved to xtol = 1e-4", _vkh_loose,
    "the limit checks on (6.3.5-1) sit at 1e-8 and 1e-10, which is Brent's"
    " tolerance and not physics; this is the row that moves them, and it is why"
    " they are reported as solver checks rather than as agreement")


def _newton_loose():
    return block_routes(newton_tol=1e-4)


brk("Newton stopped at an update of 1e-4", _newton_loose,
    "the second axis that carries error here is the solver tolerance, not a"
    " time step - this page has no time axis")


def _bvp_loose():
    return {k: v for k, v in block_routes(bvp_tol=1e-4).items()
            if k.startswith(("fv_vs", "ref_vs"))} | {
        "inst_deficit_routes_max_rel_diff":
            block_lev(bvp_tol=1e-4)["inst_deficit_routes_max_rel_diff"]}


brk("collocation tolerance loosened to 1e-4", _bvp_loose,
    "an independent route is only independent if it is also accurate; this"
    " prices the second solver's own error")

for label, got, note in BREAKS:
    print(f"* {label}")
    print(f"    {({k: v for k, v in got.items() if not k.startswith('_')})}")
    if note:
        print(f"    -> {note}")'''))

# ------------------------------------------------------------ agreement.json
cells.append(code('''METRICS = {
    # ---- where the approximation is worst, and the assumption behind it
    "vkh_worst_rel_error_pct": W["vkh_worst_rel_error_pct"],
    "vkh_worst_gamma": W["vkh_worst_gamma"],
    "vkh_worst_q": W["vkh_worst_q"],
    "vkh_worst_on_printed_labels_pct": W["vkh_worst_on_printed_labels_pct"],
    "froment_10pct_claim_margin": W["froment_10pct_claim_margin"],
    "beta_interface_at_worst": W["beta_interface_at_worst"],
    "gamma_at_half_depletion": W["gamma_at_half_depletion"],
    "gamma_at_99pct_depletion": W["gamma_at_99pct_depletion"],
    "vkh_error_at_half_depletion_pct": W["vkh_error_at_half_depletion_pct"],
    "vkh_error_at_99pct_depletion_pct": W["vkh_error_at_99pct_depletion_pct"],
    "vkh_error_at_full_depletion_pct": W["vkh_error_at_full_depletion_pct"],
    # ---- the OTHER wording of the same assumption, book p. 340
    "beta_constancy_at_worst": W["beta_constancy_at_worst"],
    "beta_constancy_at_99pct_depletion": W["beta_constancy_at_99pct_depletion"],
    # ---- Levenspiel's three printed expansions
    "lev_small_M_1pct_threshold": L["lev_small_M_1pct_threshold"],
    "lev_pfo_expansion_err_pct": L["lev_pfo_expansion_err_pct"],
    "vkh_form_pfo_expansion_err_pct": L["vkh_form_pfo_expansion_err_pct"],
    "uncorrected_E_eq_M_err_pct": L["uncorrected_E_eq_M_err_pct"],
    "lev_inst_asymptote_ratio_vkh": L["lev_inst_asymptote_ratio_vkh"],
    "lev_inst_asymptote_ratio_num": L["lev_inst_asymptote_ratio_num"],
    "lev_inst_overstatement_factor": L["lev_inst_overstatement_factor"],
    "inst_deficit_exponent_num_steepest": L["inst_deficit_exponent_num_steepest"],
    "inst_deficit_exponent_vkh_finest": L["inst_deficit_exponent_vkh_finest"],
    "inst_deficit_routes_max_rel_diff": L["inst_deficit_routes_max_rel_diff"],
    # ---- the numerics, including the baseline's own error
    "fv_vs_bvp_max_rel_diff": R["fv_vs_bvp_max_rel_diff"],
    "fv_vs_bvp_max_abs_diff_beta": R["fv_vs_bvp_max_abs_diff_beta"],
    "ref_vs_bvp_max_rel_diff": R["ref_vs_bvp_max_rel_diff"],
    "ref_self_refinement_max_rel_diff": R["ref_self_refinement_max_rel_diff"],
    "first_integral_max_abs_resid": R["first_integral_max_abs_resid"],
    "vkh_two_routes_max_rel_diff": R["vkh_two_routes_max_rel_diff"],
    "vkh_pfo_limit_max_rel_error": R["vkh_pfo_limit_max_rel_error"],
    "vkh_instantaneous_limit_rel_error": R["vkh_instantaneous_limit_rel_error"],
    "fv_convergence_order": R["fv_convergence_order"],
    "fv_rel_error_n800": R["fv_rel_error_n800"],
    "fv_rel_error_n800_uniform": R["fv_rel_error_n800_uniform"],
    # ---- the two charts, and the page next door
    "chart_label_max_plateau_ratio": CHART_LABEL_MAX_RATIO,
    "f31_definition_vkh_max_error_pct": F["f31_definition_vkh_max_error_pct"],
    "f31_definition_rootfound_pct": F["f31_definition_rootfound_pct"],
    "pfo_max_rel_error_graded": F["pfo_max_rel_error_graded"],
    "pfo_max_rel_error_uniform_same_n": F["pfo_max_rel_error_uniform_same_n"],
    "pfo_max_rel_error_uniform_half_n": F["pfo_max_rel_error_uniform_half_n"],
}
for _k in ("f31_stored_vs_recomputed_rel", "f31_first_order_error_ratio"):
    if _k in F:
        METRICS[_k] = F[_k]

# ---- coverage: GENERATED from the break rows' measured moves ---------------
COVERAGE, UNKNOWN = {}, []
for _label, _got, _note in BREAKS:
    if not isinstance(_got, dict):
        continue
    for _k, _v in _got.items():
        if _k.startswith("_"):
            continue
        if _k not in METRICS:
            UNKNOWN.append((_label, _k))
            continue
        _b = METRICS[_k]
        _rel = abs(float(_v) - _b) / max(abs(_b), abs(float(_v)), 1e-300)
        if _rel > MOVE_TOL:
            COVERAGE.setdefault(_k, []).append((_label, _rel))
assert not UNKNOWN, f"break rows recompute names that are not metrics: {UNKNOWN}"

# ---- how much the pre-extrapolation baseline actually moved, COUNTED --------
# The page used to say "nine reported numbers".  That was an undercount, and
# there is no reason to count them by hand: the "one graded solve" row IS the
# old baseline, so the movers are a set difference over METRICS.
_NOEXT = next(g for lbl, g, _ in BREAKS
              if lbl.startswith("the reference taken as one graded solve"))
_MOVERS = [(k, METRICS[k], float(v)) for k, v in sorted(_NOEXT.items())
           if k in METRICS
           and abs(float(v) - METRICS[k])
           / max(abs(METRICS[k]), abs(float(v)), 1e-300) > MOVE_TOL]
N_MOVED = len(_MOVERS)
N_MOVED_4DP = sum(1 for _k, _b, _v in _MOVERS if round(_b, 4) != round(_v, 4))
print(f"\\nwhat the pre-extrapolation baseline moved, counted rather than asserted:"
      f" {N_MOVED} of the {len(METRICS)} metrics")
print(f"move by more than the {MOVE_TOL:g} coverage tolerance, and {N_MOVED_4DP} of"
      f" them change in their fourth decimal:")
for _k, _b, _v in _MOVERS:
    print(f"   {_k:44s} {_v:>16.6g} -> {_b:<16.6g}"
          f"{'   4th decimal' if round(_b, 4) != round(_v, 4) else ''}")

# ---- and the rows must COMPUTE what they return ---------------------------
# A row that returns a typed constant records a move of exactly 1.0 whatever
# the reported value is, so its coverage links cannot fail: a hand-written
# coverage claim wearing the generator's clothes.  This guard parses each row's
# own source, plus one level of the notebook helpers it names, and rejects any
# agreement.json key bound to a numeric literal or to a local name that is only
# ever assigned one.  It cannot say whether the expression is the RIGHT one -
# that is what the measured moves are for.
import ast
import inspect
import textwrap


def _is_number(node):
    if isinstance(node, ast.Constant):
        return isinstance(node.value, (int, float)) and not isinstance(node.value, bool)
    return (isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub))
            and _is_number(node.operand))


def _key_bindings(src):
    tree, out = ast.parse(textwrap.dedent(src)), []
    for n in ast.walk(tree):
        if isinstance(n, ast.Dict):
            out += [(k.value, v) for k, v in zip(n.keys, n.values)
                    if isinstance(k, ast.Constant) and isinstance(k.value, str)]
        elif isinstance(n, ast.DictComp):
            # {k: -2.8 for k in ("metric_a", "metric_b")} - the keys are not
            # literals, so bind every constant string in the comprehension's
            # own iterables to its value expression.  A verification found this
            # form open when the rest were closed; no row used it.
            for gen in n.generators:
                out += [(c.value, n.value) for c in ast.walk(gen.iter)
                        if isinstance(c, ast.Constant) and isinstance(c.value, str)]
        elif isinstance(n, ast.Call) and getattr(n.func, "id", "") == "dict":
            out += [(kw.arg, kw.value) for kw in n.keywords if kw.arg]
        elif isinstance(n, ast.Assign):
            for t in n.targets:
                if (isinstance(t, ast.Subscript) and isinstance(t.slice, ast.Constant)
                        and isinstance(t.slice.value, str)):
                    out.append((t.slice.value, n.value))
    return out


def _sources(fn):
    src = textwrap.dedent(inspect.getsource(fn))
    seen = [src]
    for n in ast.walk(ast.parse(src)):
        if isinstance(n, ast.Name):
            obj = globals().get(n.id)
            if (inspect.isfunction(obj) and getattr(obj, "__module__", "") == "__main__"
                    and obj is not fn):
                try:
                    seen.append(textwrap.dedent(inspect.getsource(obj)))
                except (OSError, TypeError):
                    seen.append(f"# UNREADABLE SOURCE: {n.id}")
    return seen


def literal_metrics(fn):
    bad = []
    for src in _sources(fn):
        frozen = {t.id for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Assign)
                  for t in n.targets if isinstance(t, ast.Name) and _is_number(n.value)}
        bad += [k for k, v in _key_bindings(src) if k in METRICS
                and (_is_number(v) or (isinstance(v, ast.Name) and v.id in frozen))]
    return sorted(set(bad))


def _negative_control():
    """THE GUARD'S TEETH: a row that types its metrics instead of computing them."""
    return {"vkh_worst_rel_error_pct": -2.8, "lev_inst_overstatement_factor": 44.0,
            # the third form is a DICT COMPREHENSION, which walked past the
            # guard until a verification pointed at it
            **{k: 2.115 for k in ("beta_constancy_at_worst",)},
            "_diagnostic": len(LABELS)}


CAUGHT = literal_metrics(_negative_control)
LITERAL_ROWS = {lbl: bad for lbl, fn in BREAK_FNS for bad in [literal_metrics(fn)] if bad}
assert set(CAUGHT) == {"vkh_worst_rel_error_pct", "lev_inst_overstatement_factor",
                       "beta_constancy_at_worst"}, (
    "the literal guard no longer catches the row it was written for")
assert not LITERAL_ROWS, f"break rows returning typed constants: {LITERAL_ROWS}"
print(f"{len(BREAK_FNS)} break rows parsed; {len(LITERAL_ROWS)} return a metric as a"
      f" typed constant.  Negative control caught: {CAUGHT}")
print("WHAT THE GUARD CANNOT DO, stated so it is not over-trusted: it rejects a"
      " numeric LITERAL,"
      "\\nand a local name only ever bound to one, walking one level into the"
      " helpers a row names."
      "\\nIt sees literals in a dict display, in dict(...), in an out[k] = ..."
      " assignment, through a"
      "\\n{**d} merge, and - since a verification found this one open while the"
      " rest were closed - in a"
      "\\nDICT COMPREHENSION over literal key strings, which is the third form the"
      " negative control uses."
      "\\nWhat still passes: a comprehension over a NAME rather than a literal,"
      " float(-2.8), -2.8 * 1.0,"
      "\\nnp.float64(-2.8), and a literal returned two call levels down.  And it can"
      " never tell whether"
      "\\na computed expression is the RIGHT expression -"
      "\\nonly the measured moves below can, and only for rows that move something."
      "  All"
      f" {len(BREAK_FNS)} rows"
      "\\nhere perturb an input and re-enter block_worst / block_lev / block_f31 /"
      " block_routes.")

UNCOVERED = sorted(set(METRICS) - set(COVERAGE))
BEST = {k: max(r for _, r in v) for k, v in COVERAGE.items()}
print(f"\\n{len(METRICS)} metrics, {len(COVERAGE)} covered by at least one measured"
      f" move, {len(UNCOVERED)} uncovered: {UNCOVERED}")
print(f"{sum(len(v) for v in COVERAGE.values())} measured row-metric links from"
      f" {len(BREAKS)} rows")
CI_TOL = 0.05          # check_agreement.py's REL_TOL: a smaller move is invisible to CI
WEAK = sorted(k for k, r in BEST.items() if r < CI_TOL)
print(f"\\n{'metric':44s}{'rows':>6}{'best move':>12}")
for k in sorted(METRICS):
    n = len(COVERAGE.get(k, []))
    print(f"{k:44s}{n:6d}{BEST.get(k, 0.0):12.3e}"
          f"{'   < CI 5 %' if k in WEAK else ''}")
print(f"\\nweakest cover: {min(BEST.items(), key=lambda t: t[1])}")

# ---- what the relabelled uniform row still buys, measured not asserted -----
U_LABEL = "uniform grid instead of graded, at the same n_x"
U_LINKS = sorted(((r, k) for k, v in COVERAGE.items() for lbl, r in v
                  if lbl == U_LABEL), reverse=True)
U_REAL = [(r, k) for r, k in U_LINKS if r >= CI_TOL]
U_NOISE = [(r, k) for r, k in U_LINKS if r < CI_TOL]
U_BEST_FOR = [k for r, k in U_NOISE if r >= BEST[k]]
print(f"\\nthe relabelled uniform row, split by size: {len(U_REAL)} real links"
      f" (>= CI's {CI_TOL:.0%}) - the high-gamma")
print(f"metrics grading exists for, which is why the row is kept - and"
      f" {len(U_NOISE)} noise-scale links from")
print(f" {min(r for r, _ in U_NOISE):.1e} to {max(r for r, _ in U_NOISE):.1e},"
      f" invisible to CI and none of them any metric's best cover"
      f" ({not U_BEST_FOR}):")
print(f"   {', '.join(k for _, k in U_NOISE)}")
print("(that none of them is a best cover is not a check that can fail: the WEAK"
      " assertion above")
print(" already refuses any metric whose best cover is below CI's threshold.  It is"
      " printed because")
print(" the row's note used to say the relabelling removed these links, and it did"
      " not.)")
print("What the row buys NOTHING on is the metrics that carry the headline, where it"
      " agrees with the")
print("extrapolated reference to inside the move tolerance - which IS the claim its"
      " note makes:")
_U = {k: r for r, k in U_LINKS}
for _k in ("vkh_worst_rel_error_pct", "froment_10pct_claim_margin",
           "vkh_worst_on_printed_labels_pct", "gamma_at_half_depletion",
           "gamma_at_99pct_depletion"):
    print(f"   {_k:44s}"
          f"{'   still covered, ' + format(_U[_k], '.1e') if _k in _U else '   no link'}")
FLOOR = 1e-12
BELOW_FLOOR = sorted(k for k, v in METRICS.items() if abs(v) < FLOOR)
print(f"below check_agreement.py's ABS_FLOOR = {FLOOR:g}, so outside the regression"
      f" suite: {BELOW_FLOOR}")
print("  its above-floor companions, which come from the same two equations and"
      "\\n  ARE inside the suite: vkh_pfo_limit_max_rel_error,"
      " vkh_instantaneous_limit_rel_error")
assert not UNCOVERED, f"uncovered metrics: {UNCOVERED}"
assert not WEAK, f"metrics whose best cover is below CI's 5 %: {WEAK}"

report_agreement("F3.2", METRICS)'''))

# ------------------------------------------------------------- prose assertion
cells.append(code('''# Every number written in the markdown of this notebook, in meta.yaml, in
# README.md and in models_entry.yaml is checked here against the live
# computation.  THE NOTEBOOK FAILS TO EXECUTE if any of them drifts.
CLAIMS = [
    ("worst VKH error [%]", -2.7955, W["vkh_worst_rel_error_pct"], 5e-4),
    ("worst gamma", 4.7279, W["vkh_worst_gamma"], 5e-4),
    ("worst q", 1.5217, W["vkh_worst_q"], 5e-4),
    ("E_i at the worst point", 2.5217, 1 + W["vkh_worst_q"], 5e-4),
    ("beta(0) at the worst point", 0.1647, W["beta_interface_at_worst"], 5e-4),
    ("beta_w/beta(0) at half depletion", 1.208,
     W["_beta_constancy_at_half"], 5e-4),
    ("beta_w/beta(0) at the worst point", 2.115,
     W["beta_constancy_at_worst"], 5e-4),
    ("beta_w/beta(0) at 99 % depletion", 19.20,
     W["beta_constancy_at_99pct_depletion"], 5e-3),
    ("beta(0) at gamma = 100", 1.367e-11, B_MID, 5e-15),
    ("beta_w/beta(0) at gamma = 100", 3.684e9, C_MID, 5e6),
    ("the fine grid alone, same ratio, same point", 3.657e9, C_FINE_MID, 5e6),
    ("VKH error at gamma = 100 [%]", -0.0373, E_MID, 5e-5),
    ("how much smaller that error is", 75.0,
     abs(W["vkh_worst_rel_error_pct"] / E_MID), 0.5),
    ("margin on the printed 10 percent", 3.5772, W["froment_10pct_claim_margin"], 5e-4),
    ("worst on the printed labels [%]", -2.7294,
     W["vkh_worst_on_printed_labels_pct"], 5e-4),
    ("gamma at half depletion", 2.2036, W["gamma_at_half_depletion"], 5e-4),
    ("error at half depletion [%]", -1.8304,
     W["vkh_error_at_half_depletion_pct"], 5e-4),
    ("gamma at 99 % depletion", 13.0378, W["gamma_at_99pct_depletion"], 5e-4),
    ("error at 99 % depletion [%]", -1.5035,
     W["vkh_error_at_99pct_depletion_pct"], 5e-4),
    ("error at full depletion [%]", -0.00037,
     W["vkh_error_at_full_depletion_pct"], 5e-6),
    ("small-M 1 % threshold", 0.7231, L["lev_small_M_1pct_threshold"], 5e-4),
    ("Levenspiel pfo expansion [%]", -0.5366, L["lev_pfo_expansion_err_pct"], 5e-4),
    ("VKH-form pfo expansion [%]", -0.7374,
     L["vkh_form_pfo_expansion_err_pct"], 5e-4),
    ("E = M_H at the branch boundary [%]", 9.3004,
     L["uncorrected_E_eq_M_err_pct"], 5e-4),
    ("the two books' groups differ by [points]", 0.20,
     abs(L["lev_pfo_expansion_err_pct"] - L["vkh_form_pfo_expansion_err_pct"]), 5e-3),
    ("instantaneous ratio, VKH", 0.9953, L["lev_inst_asymptote_ratio_vkh"], 5e-4),
    ("instantaneous ratio, film", 0.2286, L["lev_inst_asymptote_ratio_num"], 5e-4),
    ("overstatement at M_H = 800", 44.24, L["lev_inst_overstatement_factor"], 5e-2),
    ("steepest measured exponent", -5.3385,
     L["inst_deficit_exponent_num_steepest"], 5e-4),
    # all four local exponents, not only the steepest: the other three are
    # printed in README.md and models_entry.yaml and were outside this guard.
    ("local exponent, octave 50 -> 100", -1.9162, L["_exponents_num"][0], 5e-4),
    ("local exponent, octave 100 -> 200", -2.4396, L["_exponents_num"][1], 5e-4),
    ("local exponent, octave 200 -> 400", -3.3385, L["_exponents_num"][2], 5e-4),
    ("local exponent, octave 400 -> 800", -5.3385, L["_exponents_num"][3], 5e-4),
    ("VKH exponent, finest octave", -1.9949,
     L["inst_deficit_exponent_vkh_finest"], 5e-4),
    ("two methods on the smallest deficit", 8.61e-5,
     L["inst_deficit_routes_max_rel_diff"], 5e-7),
    ("the deficit E = E_i is already within [%]", 0.00134,
     100.0 * L["_three"][0] / 20.0, 5e-6),
    ("pymrm n_x = 800 vs collocation", 5.99e-6, R["fv_vs_bvp_max_rel_diff"], 5e-8),
    ("the REFERENCE vs collocation", 4.37e-8, R["ref_vs_bvp_max_rel_diff"], 5e-10),
    ("the REFERENCE vs a finer extrapolation", 3.84e-8,
     R["ref_self_refinement_max_rel_diff"], 5e-10),
    ("how much better the reference is than one grid", 137.0,
     R["fv_vs_bvp_max_rel_diff"] / R["ref_vs_bvp_max_rel_diff"], 0.5),
    ("the single graded n_x = 800 solve's own error", 1.9e-6,
     R["fv_rel_error_n800"], 5e-8),
    ("the first integral, residual", 6.4e-10,
     R["first_integral_max_abs_resid"], 5e-11),
    ("grid order", 2.0000, R["fv_convergence_order"], 5e-4),
    ("chart label ratio", 2.0, CHART_LABEL_MAX_RATIO, 1e-12),
    ("F3.1 definition recomputed [%]", 2.0769,
     F["f31_definition_vkh_max_error_pct"], 5e-4),
    ("F3.1 definition recomputed, 6 dp [%]", 2.076896,
     F["f31_definition_vkh_max_error_pct"], 5e-6),
    ("beta(0) at gamma = 1000 is numerically zero", 1.0,
     float(abs(W["_beta_at_full_depletion"]) < 1e-40), 0.0),
    ("VKH under-estimates over the whole map", 1.0, float(VKH_ALWAYS_UNDER), 0.0),
    ("pseudo-first-order check, graded grid", 7.088e-6,
     F["pfo_max_rel_error_graded"], 5e-9),
    ("pseudo-first-order check, uniform, same n_x", 8.847e-4,
     F["pfo_max_rel_error_uniform_same_n"], 5e-7),
    ("what grading is worth on that check", 124.8,
     F["pfo_max_rel_error_uniform_same_n"] / F["pfo_max_rel_error_graded"], 5e-2),
    ("what the 400 -> 800 refinement is worth on it", 7.163,
     F["pfo_max_rel_error_uniform_half_n"]
     / F["pfo_max_rel_error_uniform_same_n"], 5e-3),
    ("what grading COSTS at the grid-study point", 248.7,
     R["fv_rel_error_n800"] / R["fv_rel_error_n800_uniform"], 5e-2),
    # the count README.md and models_entry.yaml quote for what the old
    # single-grid baseline moved - generated above, not counted by hand
    ("metrics the old baseline moves", 27.0, float(N_MOVED), 0.0),
    ("of those, in their fourth decimal", 19.0, float(N_MOVED_4DP), 0.0),
]
if F31 is not None:
    CLAIMS.append(("graded over F3.1's uniform grid", 893.9,
                   F["f31_first_order_error_ratio"], 5e-2))
    CLAIMS.append(("F3.1's stored first-order error", 6.337e-3,
                   F31["metrics"]["first_order_max_rel_error"], 5e-7))
    CLAIMS.append(("F3.1's stored VKH error [%]", 2.076893,
                   F31["metrics"]["vkh_max_error_pct"], 5e-6))
    CLAIMS.append(("F3.1's first-order check reproduced to", 6.3e-15,
                   abs(F["pfo_max_rel_error_uniform_half_n"]
                       / F31["metrics"]["first_order_max_rel_error"] - 1.0), 5e-16))
    CLAIMS.append(("the single-grid recomputation's offset from F3.1", 0.000214,
                   abs(F_PLAIN["f31_definition_vkh_max_error_pct"]
                       - F31["metrics"]["vkh_max_error_pct"]), 5e-7))
    CLAIMS.append(("the same recomputation on F3.1's OWN grid, cold [%]", 2.076893,
                   F_F31GRID["f31_definition_vkh_max_error_pct"], 5e-6))
    CLAIMS.append(("which reproduces F3.1's stored value to better than 1e-10", 1.0,
                   float(F_F31GRID["f31_stored_vs_recomputed_rel"] < 1e-10), 0.0))
if F31 is not None:
    CLAIMS.append(("F3.1 stored vs recomputed", 1.29e-6,
                   F["f31_stored_vs_recomputed_rel"], 5e-9))
    CLAIMS.append(("points bought by root-finding and by q", 0.7186,
                   abs(W["vkh_worst_rel_error_pct"])
                   - F["f31_definition_vkh_max_error_pct"], 5e-4))
    CLAIMS.append(("of which root-finding gamma", 0.0016,
                   F["f31_definition_rootfound_pct"]
                   - F["f31_definition_vkh_max_error_pct"], 5e-4))
    CLAIMS.append(("of which letting q move", 0.7170,
                   abs(W["vkh_worst_rel_error_pct"])
                   - F["f31_definition_rootfound_pct"], 5e-4))
BAD = [(n, w, g) for n, w, g, t in CLAIMS if not abs(w - g) <= t]
for n, w, g, t in CLAIMS:
    print(f"  {'ok ' if abs(w - g) <= t else 'DRIFT'} {n:42s} written {w:>12}"
          f"   computed {g:.10g}")
assert not BAD, f"prose has drifted from the computation: {BAD}"
print(f"\\nall {len(CLAIMS)} numbers written in prose match the live computation")'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**Not the equation.** Van Krevelen and Hoftijzer's relation is arithmetic in
one unknown and a bisection solves it; nothing in this repository is needed for
that, and the page says so.

**The thing the chart cannot give you: the reference solution.** The whole
point of an approximate enhancement factor is that the film equations are
coupled and nonlinear, so a designer has an approximation and no way to price
it. Those equations are 30 lines of pymrm - two fields, per-species boundary
conditions in one pair of dictionaries, a graded grid because the reaction
layer is thin, `newton` and `NumJac` for the pointwise product. Once they are
solved, three claims that have been printed and reprinted for decades become
measurements:

- the accuracy claim is checked rather than believed: $-2.7955$ % worst against
  the printed 10 percent, over the chart's whole printed domain, at a point
  found by root-finding and not by sampling;
- the stated assumption is checked rather than repeated - in **both** of the
  book's two wordings of it - and neither orders the error: at the worst point
  $B$ is at 0.1647 of bulk and $\bar\beta$ is 2.115 times $\beta(0)$, while
  where $B$ is completely gone, and the constancy assumption is violated by
  orders of magnitude, the approximation is exact to $-0.00037$ %;
- one of Levenspiel's asymptotic corrections is shown to expand the
  approximation rather than the equations, which is only visible if you have
  the equations' own solution to compare against - 44.24 times too large at
  $M_H = 800$, on a deficit two independent methods agree on. And the reason is
  exact, not empirical: the film's first integral makes the deficit
  $q\,\beta(0)$, which is exponentially small, while
  $E_i^2(E_i-1)/M_H^2$ is *algebraically* the large-$\gamma$ asymptote of
  (6.3.5-1) itself.

**And the grid is not a detail - in both directions.** `F3.1` notes that its
own first-order error is "concentrated at high Ha where the reaction layer is
thinner than a cell". Because `construct_grad` and `construct_div` take
arbitrary face positions, moving to a graded grid is one line. How much that is
worth is *measured* rather than asserted, and the measurement is split rather
than lumped: on `F3.1`'s own first-order check, over its own $\mathrm{Ha}$
range, grading at the same $n_x$ is worth **124.8x** and the refinement from
$n_x = 400$ to 800 a further 7.163x, which multiply to the 893.9x against the
6.337e-3 loaded from `F3.1`'s `agreement.json`. Crediting all 893.9x to grading
- as an earlier version of this page did - is wrong by a factor of seven.

And grading is a **liability** at the moderate $\gamma$ where this page's
headline sits: at the grid-study point a graded grid of the same size is 248.7
times *less* accurate than a uniform one, because its outer cells are nearly
eight times wider. That is why the baseline here is a Richardson extrapolation of two
graded grids rather than either grid alone. The single graded solve this page
first shipped was wrong in the fourth decimal of every moderate-$\gamma$
number - by less than its own pymrm-versus-collocation agreement, which is
precisely why a second solver family does not catch this class and a
convergence study does. It is the "extrapolation switched off" row of the break
table.

**What this page cannot conclude.** Nothing about absorbers. Both sides of
every comparison are equations, the film model itself is never tested against a
measurement here, and the 1948 origin was not read - so nothing here is
evidence about what Van Krevelen and Hoftijzer actually wrote, only about the
equation two monographs attribute to them.
"""))

# ------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

| Change | Where |
|---|---|
| A different reaction order | change `rate` in `Film.solve`; Froment book p. 341 prints Hikita and Asai's $m,n$ form |
| Reaction in the bulk, $C_{Ab} > 0$ | change `d` in the right-hand boundary dictionary; Froment's (6.3.2-10) is the pseudo-first-order version |
| Gas-side resistance | make the interface a Robin condition; Froment's (6.3.2-12) is the lumped result |
| Reversible reaction | add the product as a third field; the `NumJac` shape becomes `(n_x, 3)` |
| Surface renewal instead of a film | replace the steady BVP with a transient penetration solve (`S4`) - that is `A3.2`/`A3.3` |
| A column, not a film | this is the local closure inside `F2.1`/`F3.5` |

**Related pages:** `F3.1` Hatta number and the regimes (owns Hatta, and the
2.1 % reconciled above) · `F3.3` DeCoursey · `F3.5` CO2-amine absorption ·
`B1.1` Thiele modulus, the same `S3` equation with the boundary conditions
swapped · `A3.2`/`A3.3` the penetration and surface-renewal alternatives to the
film.

## References

**Origin, cited but not consulted.** Van Krevelen, D. W. and Hoftijzer, P. J.
(1948). *Kinetics of gas-liquid reactions. Part I. General theory.* Recueil des
Travaux Chimiques des Pays-Bas **67**(7), 563-586.
[doi:10.1002/recl.19480670708](https://doi.org/10.1002/recl.19480670708). Not
on disk; not read. **Where the extra bibliographic detail comes from:** both
books on disk print only "67, 563 (1948)". The issue number, the end page 586
and the DOI are a **Crossref record lookup** (the DOI resolved through
`api.crossref.org`, which returns volume 67, issue 7, pages 563-586, both
authors and the 1948 date), not a reading of the article - on a page whose
first section is about not writing anything it did not read, that distinction
is worth making explicit. Crossref also gives the author's name a third
capitalisation, "D. W. van Krevelen"; see the printed-defect section.

**Read, and the source of every equation here.** Froment, G. F., De Wilde, J.
and Bischoff, K. B. (2011). *Chemical Reactor Analysis and Design*, 3rd
edition. John Wiley & Sons. ISBN 978-0-470-56541-4. Chapter 6, sections 6.3.1,
6.3.2 (book pp. 328-332), 6.3.3 (book p. 335) and 6.3.5 (book pp. 340-341);
reference list book p. 365.

**Read, as an independent check on the constants.** Levenspiel, O. (1999).
*Chemical Reaction Engineering*, 3rd edition. John Wiley & Sons. ISBN
0-471-25424-X. Chapter 23, Fig. 23.4 and its annotation boxes (book p. 530);
chapter reference list (book p. 537).

**Named in the sources, not evaluated here.** Porter, K. E. (1966). *Trans.
Inst. Chem. Eng.* **44**, T25. Kishinevskii, M. Kh. *et al.* (1971). Alper, E.
(1973). Hikita, H. and Asai, S. (1964). All four are cited on Froment book
pp. 340-341 as explicit alternatives to the implicit relation; this page solves
the implicit one instead of replacing it, so none of them is computed.
"""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.kernelspec = {"display_name": "Python 3", "language": "python",
                          "name": "python3"}
nb.metadata.language_info = {"name": "python"}
nbf.write(nb, "index.ipynb")
print(f"wrote index.ipynb with {len(cells)} cells")
