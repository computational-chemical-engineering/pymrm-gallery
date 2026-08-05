#!/usr/bin/env python3
"""Generate index.ipynb for page A3.6 (Calderbank & Moo-Young 1961). Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "Calderbank & Moo-Young: one regime split, two powers of the diffusivity, and the arithmetic that does not close"
description: "The correlation that says k_L depends on nothing but the physical properties - written twice, once for rigid interfaces with Sc^-2/3 and once for mobile ones with Sc^-1/2. This page reproduces the authors' own two theoretical derivations of the small-bubble constant (0.28 and 0.34) from scratch, rebuilds the exponent split as a pymrm boundary layer and maps the transition between the two that the paper only sketches, and reports the two places where the paper's own printed arithmetic does not close: eq. (6) cannot produce eq. (7), and eq. (9) does not meet Chilton-Colburn where the paper says it does."
categories: [sec:A, struct:S3, struct:S4, tier:T0, data:tier6, phase:gas-liquid]
date: 2026-08-05
---

# Calderbank & Moo-Young: one regime split, two powers of the diffusivity

**Catalog ID:** `A3.6` · **Structures:** `S3` (1D steady BVP), `S4` (1D marched PDE) · **Tier:** T0

Calderbank & Moo-Young's result is one sentence long: in a dispersion whose
particles move freely under gravity, the continuous-phase mass-transfer
coefficient depends on **nothing but the physical properties** — not on bubble
size, not on slip velocity, not on agitator power. They write it twice, because
there are two regimes:

$$k_L\,N_{Sc}^{2/3} = 0.31\left(\frac{\Delta\rho\,\mu_c\,g}{\rho_c^{2}}\right)^{1/3}
\quad\text{(1)}\qquad\qquad
k_L\,N_{Sc}^{1/2} = 0.42\left(\frac{\Delta\rho\,\mu_c\,g}{\rho_c^{2}}\right)^{1/3}
\quad\text{(2)}$$

Eq. (1) is for bubbles below about 2½ mm, which "move in liquids as rigid
spheres under conditions of hindered flow in the boundary layer sense"; eq. (2)
is for bubbles above it, where "the conditions of unhindered flow envisaged by
HIGBIE obtained". **The entire content of the split is the exponent on the
Schmidt group** — $2/3$ against $1/2$, i.e. $k_L \propto D^{2/3}$ against
$k_L \propto D^{1/2}$.

That structure is what makes this case worth building, and it is also what makes
it dangerous to check: a single evaluation point sits on one branch and can say
nothing about the other. So this page does four separate things.

1. **Reproduces the two derivations the authors do themselves**, in the note
   added in proof, of the constant in eq. (1). Froessling combined with Allen's
   rise velocity gives them **0.28**; Friedlander combined with Stokes' law gives
   them **0.34**. Both are recomputed here from scratch and both come back.
   Neither is a fit.
2. **Rebuilds the exponent split in pymrm** as what it physically is — a
   concentration boundary layer over an interface that is either held still or
   allowed to move — and recovers $2/3$ and $1/2$ from the solved fields, with
   the magnitudes checked against two independent closed forms.
3. **Maps the transition**, which the paper describes qualitatively (its Fig. 3
   and its remarks on surface-active agents) and never models. One dimensionless
   group controls it, and the halfway point is computed.
4. **Reports two places where the paper's own arithmetic does not close.**
   Eq. (6) as printed cannot produce eq. (7); and eq. (9) does not meet
   Chilton–Colburn at $N_{Re} = 10{,}000$, nor within ±12 %. Both are shown from
   the paper's own printed constants, with every alternative reading printed
   beside the printed one, and neither is repaired.

**Nothing here tests the three fitted constants.** 0.31, 0.42 and 0.13 are
regression lines drawn through Figs 1, 2 and 6, and the only way to judge them is
against the points in those figures. That would be a goodness of fit, not a
validation, and it is scoped out — see *The data*."""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

### What the paper set out to do

By 1960 the Froessling and Higbie pictures of interphase mass transfer both
existed and neither had been tested on gas–liquid dispersions, for a simple
reason: nobody could measure interfacial area in a bubble swarm, so nobody could
separate $k_L$ from $k_L a$. Calderbank's earlier work supplied the missing
half — a light-transmission method for interfacial area, a light-reflection
method for optically dense dispersions, and a γ-ray transmission method for gas
hold-up, from which the Sauter mean bubble size follows. With area in hand, $k_L$
in a swarm becomes measurable.

The paper's finding is stated flatly on journal page 42:

> It will be seen that both correlations show that the liquid phase mass transfer
> coefficients are independent of bubble size and slip velocity and depend only
> on the physical properties of the system.

and the mechanism, on page 41:

> It was concluded that small "rigid-sphere" bubbles experienced friction drag
> causing hindered flow in the boundary layer sense, and that under these
> circumstances the mass transfer coefficient was proportional to the 2/3 power
> of the diffusion coefficient, as found by FROESSLING and others. For large
> bubbles (> 2½ mm diameter) form drag predominated and the conditions of
> unhindered flow envisaged by HIGBIE obtained, and as postulated by HIGBIE the
> mass transfer coefficient was proportional to the ½ power of the diffusion
> coefficient.

So the split is not a curve-fitting convenience. It is a claim about which of two
mechanisms is operating, and the discriminating evidence is a **power of $D$**.

### Where this page's numbers come from, and the citation trap on this file

Every constant, exponent and stated result used here was read off `pdftoppm`
renders of the Elsevier scan **at 300 ppi, which is that scan's native
resolution**: `pdfimages -list` reports every image on every page as CCITT-G4
bilevel at 300×300 ppi (each page is stored as a stack of 2250 × 115 px strips).
Rendering at 600 dpi would be 2× interpolation and would add nothing. Every
numeric was then **cropped and re-read at digit scale** before being written
down, because the 1961 Pergamon typesetting sets the decimal point as a *mid-dot*
and this scan drops most of them: `0·31` comes back from the text layer as
`0 31`, and the English abstract's copy of eq. (3) is OCR-ed as `0 18` where the
French and German abstracts of the same equation both print `0,13`.

**The text layer of this PDF states the wrong year and the wrong volume.** Its
first line comes back from `pdftotext` as

```
Chemcal   En@neenng Snence, 1981,   1'01        19, pp   39 to 54
```

— 1981, volume 19. Both are OCR damage on a pre-1980 scan. The header **on the
page image**, cropped at 300 ppi and enlarged with nearest-neighbour
interpolation, reads

> Chemical Engineering Science, 1961, Vol 16, pp 39 to 54  Pergamon Press Ltd,
> London  Printed in Great Britain

and the running footer on journal pages 44, 47, 49 and 51 independently reads
*"Chem Engng Sci Vol 16, Nos 1 and 2 December, 1961"* — of which the text layer
corrupts one to 1991, the same trap a third time. The publisher PII
`0009-2509(61)87005-X` encodes 61 and the page range 39–54 matches on both
readings. **The citation on this page was established from the article's own
printed header on an image**, not from the text layer, not from the filename and
not from the PII.

### What this paper prints as numbers, and what lives only in figures

| Printed as a number | Lives only in a figure |
|---|---|
| Eqs. (1), (2), (4): the three correlation constants | Figs 1, 2, 6: the several hundred points behind them |
| $N_{Sh} = N_{Nu} = 2.0 + 0.31\,N_{Ra}^{1/3}$ (page 44) | Fig. 3: the transition region, one system per curve |
| The complete packed-bed chain: eq. (5), Chu's $\Delta p$, eqs. (6), (7), (8) | Fig. 4: transition-region froth data |
| The complete pipe chain: Blasius, eqs. (9), Chilton–Colburn, and two claims about how they compare | Fig. 5: $k_La$ on sieve trays |
| Eq. (10) and the four competing suspension-power laws | — |
| The note added in proof: Allen's ¼, Friedlander's 0.89, and the resulting **0.28** and **0.34** | — |
| Baird's theoretical $K_L = 0.975\,D_L^{1/2}(g/d)^{1/4}$ | — |
| The Range of Variables table: $k_L$, $D_L$, $\Delta\rho$, $\rho_c$, $\mu_c$, $d$ per data set | — |

**No figure is digitised on this page.** The digitisation route needs a
maintainer review and none is available; more importantly it is not needed,
because everything checkable here is printed. The one thing the figures could
settle — whether 0.31, 0.42 and 0.13 are good fits — is therefore *not settled
here*, and the page says so rather than substituting something that looks like an
answer."""))

# ------------------------------------------------------- the published model
cells.append(md(r"""## The published model

### The two correlations, and the one identity joining them

With $N_{Sc} = \mu_c/(\rho_c D_L)$, $\Delta\rho$ the phase density difference and
$g$ gravity, the paper's eqs. (1) and (2) are

$$k_L\,N_{Sc}^{2/3} = \frac{h_c}{C_P\rho_c}N_{Pr}^{2/3}
  = 0.31\left(\frac{\Delta\rho\,\mu_c\,g}{\rho_c^{2}}\right)^{1/3},\qquad
k_L\,N_{Sc}^{1/2} = 0.42\left(\frac{\Delta\rho\,\mu_c\,g}{\rho_c^{2}}\right)^{1/3}.$$

The heat-transfer half of eq. (1) is the Chilton–Colburn analogue and is not
separately checked here; it carries the same constant.

On page 42 eq. (2) is also written in dimensionless form as
$N_{Sh} = 0.42\,N_{Sc}^{1/2}N_{Gr}^{1/3}$, and on page 44 eq. (1) appears in a
"more precise form" that adds the stagnant limit,

$$N_{Sh} = N_{Nu} = 2.0 + 0.31\,N_{Ra}^{1/3},$$

with $N_{Ra} = N_{Gr}N_{Sc}$. **These are the same statement written twice**, and
the notebook checks the equivalence symbolically — not because the algebra is in
doubt but because it is a *transcription* check that both readings survive
together. The `2.0` is the conduction limit of a sphere in a stagnant medium and
is reproduced with pymrm below.

### The turbulence branch, and the two chains built on it

Where the dispersed phase is *not* free to move under gravity — a fixed
submerged body, or particles held beyond the point of complete suspension — the
paper switches to Kolmogoroff local isotropy. Taking the turbulence Reynolds
number $N'_{Re} = \rho_c^{1/3}(P/v)^{1/3}d^{4/3}/\mu_c$ and requiring
$k_L \propto D_L^{2/3}$ with $k_L$ independent of $d$, it finds the exponents
$x = 1/3$, $y = 3/2$ and hence

$$k_L\,N_{Sc}^{2/3} = 0.13\left[\frac{(P/v)\,\mu_c}{\rho_c^{2}}\right]^{1/4}
\tag{4}$$

"with a standard deviation of 66 per cent" over a figure spanning $10^7$ in its
variables. Eq. (4) is then pushed through two closures the authors did not fit:

**Packed beds.** The flow work per unit volume of fluid is
$P/v = (\Delta p/L)\,G/(\rho_c\Sigma)$ (eq. 5), with $\Sigma$ the voidage; Chu's
correlation gives $\Delta p/L$; the two combine to a printed $P/v$; that goes into
eq. (4) to give eq. (6); and setting $\Sigma = 0.41$ gives eq. (7),
$j_d = 0.69\,N_{Re}^{-1/3}$, to be compared with Fallat's measured
$j_d = 0.626\,N_{Re}^{-0.322}$.

**Pipes.** The Blasius equation gives $\Delta p$; the same route gives eq. (9),
$j_{d,h} = 0.058\,N_{Re}^{-0.31}$, to be compared with Chilton–Colburn's
$j_{d,h} = 0.023\,N_{Re}^{-0.20}$.

Both chains are pure algebra on printed constants, so **both can be re-run
exactly**, and that is the highest-ranked validation available for this paper.

### What the authors prove about eq. (1) without fitting anything

The note added in proof (journal pages 52–53) contains the two results that make
this page possible. Correlation (1), they say, follows from Froessling's equation
combined with Allen's free-rise velocity

$$V = \tfrac14\left(\frac{\Delta\rho^{2}g^{2}}{\rho_c\mu_c}\right)^{1/3} d,$$

whereupon "a constant of 0·28 instead of 0·31, as proposed, is thereby obtained";
and it follows again from Stokes' law combined with Friedlander's Stokes-region
result $N_{Sh} = 0.89\,(dV/D_L)^{1/3}$, "when a constant of 0·34 is obtained".

**Those are predictions, not fits**, and they are the two numbers this page
reproduces from scratch."""))

# --------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

**Units.** The paper is CGS throughout: $k_L$ in cm/s, $D_L$ in cm²/s, densities
in g/cm³, viscosity in cP. Only one input to any calculation here is *not*
printed in the paper, and it is $g = 981$ cm/s²; it is supplied by the notebook
and flagged wherever it enters.

**Three constants come from outside this paper**, and each is named where it is
used rather than buried:

| Constant | Value | Where it comes from | What it is used for |
|---|---|---|---|
| Froessling's coefficient | 0.552 | Froessling (1938), the paper's ref. [4] — **not printed in this article** | reproducing the printed 0.28 |
| Stokes' law denominator | 18 | $V = g\Delta\rho d^2/(18\mu_c)$ | reproducing the printed 0.34 |
| Allen's drag law | $C_D = 10\,N_{Re}^{-1/2}$ | Allen (1900), the paper's ref. [49] | *checking* the paper's printed ¼, not replacing it |

The first is the only one that matters, and the page handles it in both
directions: it computes 0.28 forward from Froessling's published 0.552, and it
also **inverts** the paper's own two printed numbers (0.28 and ¼) to say what
coefficient they imply. The second statement uses nothing but this article.

**The boundary-layer model's assumptions**, stated before it is built:

- A flat interface of contact length $L$, with the solute saturated at the
  interface ($c = c^*$) and the bulk clean at entry. This is the *local* picture
  the paper's two mechanisms describe; it is not a bubble.
- A velocity field $u(y) = u_s + \dot\gamma\,y$ parallel to the interface, with
  $v = 0$. **This is divergence-free exactly**, so it is a genuine solution of
  continuity and not a boundary-layer approximation with a neglected term.
  $u_s = 0$ is the rigid, no-slip interface; $\dot\gamma = 0$ is the fully mobile
  one that translates with the surrounding fluid.
- Steady state, and streamwise diffusion neglected against streamwise
  convection — the standard high-Péclet boundary-layer ordering.

**What the model deliberately does not claim.** It is a *mechanism* model, not a
bubble model. It reproduces the two exponents and the interpolation between them;
it does not predict 0.31 or 0.42, because those constants absorb the bubble
geometry, the slip velocity and the swarm hydrodynamics that this flat-interface
picture throws away. Saying otherwise would be the overclaim this page exists to
avoid."""))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

**Tier 6.** This paper contains no individual measurement. Several hundred $k_L$
values, from the authors' own aerated mixing vessels and sieve and sintered-plate
columns and from twenty-two other published sources, appear **only as points on
Figs 1, 2 and 6**. What is printed is (a) the constants and exponents, (b) a
handful of statements the authors make about their own algebra, and (c) a table
of the *range* of each variable per data set.

Two CSVs are committed, and both are transcriptions of printed numbers.

`calderbank-mooyoung-1961-printed-numbers.csv` — every constant, exponent and
stated result used on this page, with the equation or sentence it was read from.
The `kind` column separates **fitted** constants (0.31, 0.42, 0.13) from
**predicted** ones (0.28, 0.34, 2.0) and from constants merely **quoted** from
elsewhere (Chu's 35.4, Blasius's 0.0396, Fallat's 0.626, Chilton–Colburn's 0.023,
Allen's ¼, Friedlander's 0.89, Baird's 0.975).

`calderbank-mooyoung-1961-range-of-variables.csv` — the two `This work` rows of
the Range of Variables table on journal pages 42–43. **These are the two sides of
the split**: the large-bubble row (eq. 2) and the small-bubble/rigid-sphere row
(eq. 1). Their diameter columns separate at the paper's own boundary, 0.20–0.80 cm
against 0.02–0.06 cm. They are *marginal* ranges — no cell is paired with any
other — so evaluating a correlation at the corners produces a **window** that
should bracket the printed $k_L$ window, and nothing sharper than that.

**One printed exponent in that table is wrong, and the page proves it from the
paper's own correlation rather than from plausibility.** The small-bubble row
prints $D_L = 3.4$–$17.8 \times 10^{-4}$ cm²/s. The demonstration and all three
candidate exponents are computed below.

### Cross-page rule

**No dataset from another gallery page is loaded here.** Both CSVs are new and
are transcribed from this paper's own pages, so the cross-page reconciliation
obligation in `AGENTS.md` does not arise. Three published pages are *related* and
are cited as context in *Reuse*, but none of their numbers is reused, reported or
depended on: `A3.1` (Whitman two-film), `A3.3` (Danckwerts surface renewal) and
`A3.4` (Wakao–Funazkri). Where this page's boundary-layer model overlaps `A3.3`'s
penetration element the overlap is named explicitly in *What pymrm adds*, so that
the two pages cannot be read as independent confirmations of each other."""))

cells.append(code(r"""# Colab: install pymrm if it is not already present.
try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm"""))

cells.append(code(r'''import sys, time, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(4):
        if (local / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(local / "shared")); break
        local = local.parent
    else:
        url = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
               "pymrm-gallery/main/shared/gallery_utils.py")
        urllib.request.urlretrieve(url, "gallery_utils.py")
        sys.path.insert(0, ".")

import numpy as np
import matplotlib.pyplot as plt
from gallery_utils import load_data, load_meta, report_agreement

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

T_START = time.time()
np.random.seed(20260805)          # nothing here is stochastic; seeded anyway
PAGE = "A3.6-calderbank-moo-young"

PN = load_data("calderbank-mooyoung-1961-printed-numbers.csv", page=PAGE).set_index("id")
RV = load_data("calderbank-mooyoung-1961-range-of-variables.csv", page=PAGE).set_index("row_id")

P_READS = set()

def p(key):
    """A printed value, by id. Nothing on this page types a constant twice."""
    P_READS.add(key)
    return float(PN.loc[key, "printed_value"])

print("Printed numbers transcribed from the 300 ppi renders:", len(PN), "rows")
print(PN[["journal_page", "equation", "quantity", "printed_value", "kind"]]
      .to_string(max_colwidth=52))
print()
print("Range of Variables, the two sides of the split:")
print(RV[["regime", "correlation", "kL_lo", "kL_hi", "D_L_lo_mant", "D_L_hi_mant",
          "D_L_exp_printed", "mu_c_lo", "mu_c_hi", "d_lo", "d_hi"]].to_string())
'''))

# -------------------------------------------------- the printed exponent bug
cells.append(md(r"""### The small-bubble row's diffusivity exponent

The test has a **control**: the same corner evaluation is applied to the
large-bubble row, whose printed exponent is $-6$, using eq. (2). If the method
works, the control must reproduce the large-bubble row's printed $k_L$ window
without any adjustment — and it must be able to fail, which is exactly what the
$-4$ reading demonstrates."""))

cells.append(code(r'''import itertools

G_CGS = 981.0     # cm/s^2 - the ONLY input here not printed in the paper

C1, C2, C4 = p("c_eq1"), p("c_eq2"), p("c_eq4")

def kL_eq1(drho, mu, rho, D):
    """Eq. (1): k_L Sc^(2/3) = 0.31 (drho mu g / rho^2)^(1/3), CGS."""
    return C1 * (drho * mu * G_CGS / rho**2)**(1/3) * (mu / (rho * D))**(-2/3)

def kL_eq2(drho, mu, rho, D):
    """Eq. (2): k_L Sc^(1/2) = 0.42 (drho mu g / rho^2)^(1/3), CGS."""
    return C2 * (drho * mu * G_CGS / rho**2)**(1/3) * (mu / (rho * D))**(-1/2)

def corner_window(row, corr, D_exp):
    r = RV.loc[row]
    Ds = (r.D_L_lo_mant * 10.0**D_exp, r.D_L_hi_mant * 10.0**D_exp)
    vals = [corr(a, b/100.0, c, d)      # cP -> poise
            for a, b, c, d in itertools.product((r.drho_lo, r.drho_hi),
                                                (r.mu_c_lo, r.mu_c_hi),
                                                (r.rho_c_lo, r.rho_c_hi), Ds)]
    return min(vals), max(vals)

# --- the control: large bubbles, printed exponent -6, eq. (2) -------------
lo2, hi2 = corner_window("this_work_large", kL_eq2, int(RV.loc["this_work_large", "D_L_exp_printed"]))
obs2 = (RV.loc["this_work_large", "kL_lo"], RV.loc["this_work_large", "kL_hi"])
print(f"CONTROL  large bubbles, eq. (2), D exponent as printed "
      f"({int(RV.loc['this_work_large','D_L_exp_printed'])}):")
print(f"   predicted window [{lo2:.5g}, {hi2:.5g}] cm/s   printed measured [{obs2[0]}, {obs2[1]}] cm/s")
print(f"   top   ratio {hi2/obs2[1]:.4f}   bottom ratio {lo2/obs2[0]:.4f}")

# --- the test: small bubbles, three candidate exponents, eq. (1) ----------
print("\nTEST  small bubbles and rigid spheres, eq. (1):")
small_windows = {}
for e in (-4, -5, -6):
    lo, hi = corner_window("this_work_small", kL_eq1, e)
    small_windows[e] = (lo, hi)
    obs = (RV.loc["this_work_small", "kL_lo"], RV.loc["this_work_small", "kL_hi"])
    tag = "  <- AS PRINTED" if e == int(RV.loc["this_work_small", "D_L_exp_printed"]) else ""
    print(f"   D x 10^{e:<3d} -> [{lo:.5g}, {hi:.5g}] cm/s   "
          f"top ratio {hi/obs[1]:8.4f}   bottom ratio {lo/obs[0]:8.4f}{tag}")

R_HI_PRINTED_EXP = small_windows[-4][1] / RV.loc["this_work_small", "kL_hi"]
R_HI_SIX         = small_windows[-6][1] / RV.loc["this_work_small", "kL_hi"]
print(f"\nThe printed -4 puts the correlation {R_HI_PRINTED_EXP:.1f}x above the row's own "
      f"printed k_L ceiling;\n-6 lands {abs(1-R_HI_SIX)*100:.1f} % below it. "
      f"-5 is ruled out too ({small_windows[-5][1]/RV.loc['this_work_small','kL_hi']:.2f}x), "
      "so the choice is unique.\nEVERYTHING BELOW USES -6, and the page says so.")
D_EXP_SMALL = -6
'''))

# --------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Two models, and they are deliberately small, because the paper's own algebra
carries most of the validation and the pymrm work exists to test the **one thing
the algebra cannot reach** — the mechanism behind the exponents.

### 1. The stagnant limit: a sphere in a quiescent medium (`S3`)

$$\frac{1}{r^{2}}\frac{\mathrm d}{\mathrm d r}\!\left(r^{2}\frac{\mathrm d c}{\mathrm d r}\right)=0,
\qquad c(R)=c^*,\quad c(R_\infty)=0,$$

solved with `construct_grad` + `construct_div(nu=2)`. $N_{Sh} = k_L d/D_L$ must
approach the **2.0** printed in the paper's more precise form of eq. (1) as
$R/R_\infty \to 0$. The geometry code `nu` is the thing under test here, and the
break table injects `nu=0` and `nu=1`.

### 2. The regime split: a concentration boundary layer (`S4`)

$$u(y)\,\frac{\partial c}{\partial x}=D_L\,\frac{\partial^{2}c}{\partial y^{2}},
\qquad u(y)=u_s+\dot\gamma\,y,$$

with $c(x,0)=c^*$, $c(0,y)=0$ and a far-field at $y_{\max}$. This is marched in
$x$ exactly as a transient problem is marched in time, with $u(y)$ playing the
role of the accumulation coefficient — built once with
`construct_coefficient_matrix`, never inside the loop.

- $u_s = 0$: rigid interface. The only length scale is
  $\delta_{\rm L}=(9D_LL/\dot\gamma)^{1/3}$, so $k_L \propto D_L^{2/3}$ — the
  **Lévêque/Froessling** branch, eq. (1).
- $\dot\gamma = 0$: mobile interface. The only length scale is
  $\delta_{\rm H}=\sqrt{D_LL/u_s}$, so $k_L\propto D_L^{1/2}$ — the **Higbie**
  branch, eq. (2).
- Both non-zero: the transition the paper's Fig. 3 shows and never models.

**The mean coefficient is taken from the discrete mass balance, not from the
boundary gradient.** Summing the discrete equation over all cells telescopes the
diffusive divergence to the two boundary faces, so the cumulative wall flux is
*identically* $\sum_i u_i V_i c_i$ at the exit plane. That is an exact statement
about the scheme, so it is a **structural** identity and is labelled one; it is
used because it is robust, not because it proves anything. The boundary-gradient
route is computed alongside it as a second, non-structural reading."""))

cells.append(code(r'''from scipy.sparse.linalg import splu, spsolve
from scipy.special import gamma as gammafn
from pymrm import construct_grad, construct_div, construct_coefficient_matrix


def sh_sphere(ratio=1e-4, n_r=800, nu=2, stretch=8.0, outer="dirichlet"):
    """Steady diffusion from a sphere of radius R out to R/ratio.  Sh = k d / D.

    Exact for this domain: Sh = 2 / (1 - ratio),  -> 2 as ratio -> 0.
    nu is GEOMETRY: 0 Cartesian, 1 cylindrical, 2 spherical.  2 is correct here.
    """
    R, Rout = 1.0, 1.0 / ratio
    w = np.linspace(0.0, 1.0, n_r + 1)
    r_f = R + (Rout - R) * np.expm1(stretch * w) / np.expm1(stretch)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n_r, 1)                      # (space, field) - never a bare (n,)
    # OUTWARD normal, a dc/dn + b c = d.
    #   r = R      : c = c* = 1                       -> a=0, b=1, d=1
    #   r = R_inf  : c = 0 (clean bulk)               -> a=0, b=1, d=0
    #                or no flux (closed shell)        -> a=1, b=0, d=0
    bc = ({"a": 0.0, "b": 1.0, "d": 1.0},
          {"a": 0.0, "b": 1.0, "d": 0.0} if outer == "dirichlet"
          else {"a": 1.0, "b": 0.0, "d": 0.0})
    div = construct_div(shape, r_f, nu=nu)
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc)
    c = spsolve((div @ grad).tocsc(), -np.asarray((div @ grad_bc).todense()).ravel())
    dcdr = float((grad @ c.reshape(-1, 1) + grad_bc)[0, 0])
    return -2.0 * R * dcdr                # d = 2R, D = 1, c* = 1


class BoundaryLayer:
    """u(y) dc/dx = D d2c/dy2 over an interface of contact length L.

        u(y) = u_s + gdot*y,  v = 0  (divergence-free EXACTLY, not an ordering)
        c(x, 0) = 1,  c(0, y) = 0,  no flux at y_max.

    u_s = 0        -> rigid  interface, Leveque,  k ~ D^(2/3)
    gdot = 0       -> mobile interface, Higbie,   k ~ D^(1/2)
    """

    def __init__(self, D=1e-9, u_s=0.0, gdot=1.0, L=1.0, n_y=400, n_x=400,
                 stretch=6.0, y_span=16.0, y_max=None):
        self.D, self.u_s, self.gdot, self.L = D, u_s, gdot, L
        d_lev = (9 * D * L / gdot)**(1/3) if gdot > 0 else np.inf
        d_hig = np.sqrt(D * L / u_s) if u_s > 0 else np.inf
        self.y_max = y_max if y_max is not None else y_span * min(d_lev, d_hig)
        self.n_y, self.n_x = n_y, n_x
        w = np.linspace(0.0, 1.0, n_y + 1)
        # exponential stretch: fine at the interface, coarse in the far field.
        # A fixed map of w, so refining n_y refines the WHOLE grid in proportion.
        self.y_f = self.y_max * np.expm1(stretch * w) / np.expm1(stretch)
        self.y_c = 0.5 * (self.y_f[:-1] + self.y_f[1:])
        self.shape = (n_y, 1)
        self.V = np.diff(self.y_f)
        self.div = construct_div(self.shape, self.y_f, nu=0)   # nu=0: Cartesian
        self.u = u_s + gdot * self.y_c

    def run(self, n_x=None, interface="dirichlet"):
        n_x = n_x or self.n_x
        # OUTWARD normal.  y = 0 : n points in -y.
        #   saturated interface  c = c*        -> a=0, b=1, d=1
        #   (break row) zero flux              -> a=1, b=0, d=0
        # y = y_max : n points in +y, zero flux -> a=1, b=0, d=0
        bc = ({"a": 0.0, "b": 1.0, "d": 1.0} if interface == "dirichlet"
              else {"a": 1.0, "b": 0.0, "d": 0.0},
              {"a": 1.0, "b": 0.0, "d": 0.0})
        grad, grad_bc = construct_grad(self.shape, self.y_f, self.y_c, bc)
        A = (self.div @ grad).tocsc() * self.D
        b = np.asarray((self.div @ grad_bc).todense()).ravel() * self.D
        U = construct_coefficient_matrix(self.u.reshape(-1, 1), self.shape).tocsc()
        dx = self.L / n_x
        lu = splu((U / dx - A).tocsc())        # constant operator: factorised ONCE
        c = np.zeros(self.n_y)
        kx = np.empty(n_x)
        for j in range(n_x):
            c = lu.solve(self.u * c / dx + b)
            kx[j] = -self.D * float((grad @ c.reshape(-1, 1) + grad_bc)[0, 0])
        # STRUCTURAL: backward Euler telescopes the diffusive divergence onto the
        # two faces, so the cumulative wall flux IS sum(u V c) at the exit plane.
        kbar_balance = float((self.u * self.V * c).sum()) / self.L
        kbar_gradient = float(kx.mean())       # independent of the balance route
        return kbar_balance, kbar_gradient, kx[-1], c


def leveque(D, gdot, L):
    """Mean k over 0..L for u = gdot*y.  Similarity: f'' + 3 eta^2 f' = 0,
    so f'(0) = -1/Gamma(4/3) and k(x) ~ x^(-1/3), whence kbar = (3/2) k(L)."""
    return 1.5 * D**(2/3) * (gdot / (9 * L))**(1/3) / gammafn(4/3)

def higbie(D, u_s, L):
    """Mean k over 0..L for u = u_s.  Penetration with contact time L/u_s."""
    return 2.0 * np.sqrt(D * u_s / (np.pi * L))

print("operators built;  sh_sphere(1e-4) =", f"{sh_sphere(1e-4):.8f}")
'''))

# ------------------------------------------------------------------ results
cells.append(md(r"""## Results

### 1. The two theoretical derivations the authors do themselves

Both are recomputed from scratch. The Froessling route needs one constant that
this article does not print, so it is done in **both directions**."""))

cells.append(code(r'''FROESSLING = 0.552     # Froessling (1938), ref [4] of this paper -- NOT printed here
STOKES_DEN = 18.0      # V = g drho d^2 / (18 mu)
ALLEN_CD   = 10.0      # Allen (1900), ref [49]: C_D = 10 Re^(-1/2)

# --- Froessling + Allen ---------------------------------------------------
# Sh = A Re^(1/2) Sc^(1/3)  =>  k Sc^(2/3) = A sqrt(nu V / d)
# Allen: V/d = C_A (g^2 drho^2 / (rho mu))^(1/3)   =>  nu V/d = C_A (g drho mu/rho^2)^(2/3)
C_ALLEN = p("c_allen")
c28 = FROESSLING * np.sqrt(C_ALLEN)
c28_printed = p("c_froessling_allen")
print("Froessling + Allen")
print(f"   forward : 0.552 * sqrt({C_ALLEN}) = {c28:.6f}   printed {c28_printed}   "
      f"rel {c28/c28_printed-1:+.4%}")
print(f"   inverse : the paper's OWN 0.28 and 1/4 imply a Froessling coefficient of "
      f"{c28_printed/np.sqrt(C_ALLEN):.4f}")
print(f"             against Froessling's published 0.552 -> "
      f"{(c28_printed/np.sqrt(C_ALLEN))/FROESSLING-1:+.2%}")
# and a check on Allen's own printed 1/4, from the drag law he proposed
C_ALLEN_FROM_CD = (2.0 / (3.0 * ALLEN_CD / 2.0))**(2/3)
print(f"   Allen's 1/4 recomputed from C_D = 10 Re^(-1/2): {C_ALLEN_FROM_CD:.4f} "
      f"({C_ALLEN_FROM_CD/C_ALLEN-1:+.2%} on the printed 1/4)")

# --- Friedlander + Stokes -------------------------------------------------
# Sh = 0.89 (dV/D)^(1/3);  Stokes V = g drho d^2/(18 mu)
#   => k Sc^(2/3) = 0.89 (1/18)^(1/3) (g drho mu / rho^2)^(1/3)
C_FRIED = p("c_friedlander")
c34 = C_FRIED * (1.0 / STOKES_DEN)**(1/3)
c34_printed = p("c_stokes_fried")
print("\nFriedlander + Stokes")
print(f"   {C_FRIED} * (1/18)^(1/3) = {c34:.6f}   printed {c34_printed}   "
      f"rel {c34/c34_printed-1:+.4%}")
print(f"   USES NOTHING FROM OUTSIDE THIS ARTICLE except the 18 of Stokes' law.")

print(f"\nThe two theoretical predictions bracket the fitted constant: "
      f"{c34:.3f} > {C1} > {c28:.3f}")
print(f"   fitted / Froessling+Allen = {C1/c28:.4f}, "
      f"fitted / Friedlander+Stokes = {C1/c34:.4f}")
'''))

cells.append(md(r"""### 2. The regime split as a number: where the two correlations cross

Eqs. (1) and (2) differ only in the Schmidt exponent, so they give the same
$k_L$ at exactly one Schmidt number,
$N_{Sc}^\star = (0.42/0.31)^{6}$. The paper notices the crossing in words —
"the mass transfer coefficients for large and small bubbles approach each other
at low values of the Schmidt group ... TOOR and MARCHELLO ... predict that the
boundary layer and penetration mechanisms of mass transfer lead to identical
results at low values of the Schmidt group" — but never prints the number."""))

cells.append(code(r'''SC_STAR = (C2 / C1)**6
print(f"N_Sc* = (0.42/0.31)^6 = {SC_STAR:.4f}")
def ratio_eq2_over_eq1(Sc):
    """k_L from eq. (2) divided by k_L from eq. (1), at the same properties."""
    return (C2 / C1) * Sc**(2/3 - 1/2)
for Sc in (1.0, SC_STAR, 100.0, 500.0, 3000.0, 15000.0):
    print(f"   Sc = {Sc:9.2f}   k_large / k_small = {ratio_eq2_over_eq1(Sc):6.3f}")

# the Schmidt window the paper's own two data sets actually occupy
for row in ("this_work_large", "this_work_small"):
    r = RV.loc[row]
    e = int(r.D_L_exp_printed) if row == "this_work_large" else D_EXP_SMALL
    Ds = (r.D_L_lo_mant * 10.0**e, r.D_L_hi_mant * 10.0**e)
    Sc = [ (b/100.0) / (c * d) for b, c, d in
           itertools.product((r.mu_c_lo, r.mu_c_hi), (r.rho_c_lo, r.rho_c_hi), Ds) ]
    print(f"   {r.regime:32s} Sc in [{min(Sc):.4g}, {max(Sc):.4g}]")
print(f"\nBoth data sets sit far above N_Sc* = {SC_STAR:.2f}, so the crossing is an "
      "extrapolation\nof the two fitted lines, not something either data set observes.")
'''))

cells.append(code(r'''fig, ax = plt.subplots(figsize=(6.6, 4.3))
Sc = np.logspace(0, 5, 400)
grp = 1.0     # k_L / (drho mu g / rho^2)^(1/3), i.e. the correlations reduced
ax.loglog(Sc, C1 * Sc**(-2/3), lw=2, label=r"eq. (1)  $0.31\,N_{Sc}^{-2/3}$  (rigid / hindered)")
ax.loglog(Sc, C2 * Sc**(-1/2), lw=2, label=r"eq. (2)  $0.42\,N_{Sc}^{-1/2}$  (mobile / unhindered)")
ax.axvline(SC_STAR, color="0.4", ls="--", lw=1)
ax.annotate(rf"$N_{{Sc}}^\star={SC_STAR:.2f}$", (SC_STAR, C1 * SC_STAR**(-2/3)),
            xytext=(12, 26), textcoords="offset points", fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8, color="0.4"))
for row, col in (("this_work_small", "tab:blue"), ("this_work_large", "tab:orange")):
    r = RV.loc[row]
    e = int(r.D_L_exp_printed) if row == "this_work_large" else D_EXP_SMALL
    Ds = (r.D_L_lo_mant * 10.0**e, r.D_L_hi_mant * 10.0**e)
    scs = [(b/100.0)/(c*d) for b, c, d in itertools.product(
        (r.mu_c_lo, r.mu_c_hi), (r.rho_c_lo, r.rho_c_hi), Ds)]
    ax.axvspan(min(scs), max(scs), color=col, alpha=0.10)
    ax.text(np.sqrt(min(scs)*max(scs)), 4e-4, r.regime.split(" and ")[0],
            ha="center", fontsize=8, color=col)
ax.set_xlabel(r"$N_{Sc}=\mu_c/\rho_c D_L$")
ax.set_ylabel(r"$k_L\,/\,(\Delta\rho\,\mu_c g/\rho_c^2)^{1/3}$")
ax.set_title("The regime split, and where the two branches would meet")
ax.legend(fontsize=8.5, loc="lower left"); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); plt.show()
print("Shaded bands are the Schmidt windows of the paper's own two data sets, from the\n"
      "Range of Variables table. The crossing lies outside both.")
'''))

cells.append(md(r"""### 3. The packed-bed chain, step by step

Every step is printed, so every step can be re-run. The first is a **symbolic**
identity: eq. (5) applied to Chu's correlation must give the $P/v$ the paper
prints, exponent for exponent."""))

cells.append(code(r'''import sympy as sp

mu, S, G, d, rho = sp.symbols("mu Sigma G d rho", positive=True)
CHU_C, CHU_A = p("c_chu"), p("a_chu")

# Chu, as printed on journal page 50
dpdL = sp.Rational(354, 10) * (mu*(1-S)/(d*G))**sp.Rational(44, 100) \
       * G**2 * (1-S) / (d*rho*S**3)
# eq. (5): P/v = (dp/L) G / (rho Sigma)     -- per unit volume of FLUID
Pv_derived = sp.simplify(dpdL * G / (rho*S))
# the P/v the paper prints immediately underneath
Pv_printed = sp.Rational(354, 10) * mu**sp.Rational(44, 100) * (1-S)**sp.Rational(144, 100) \
             * G**sp.Rational(256, 100) / (d**sp.Rational(144, 100) * rho**2 * S**4)
CHU_RESID = float(sp.simplify(Pv_derived/Pv_printed - 1))
print("eq. (5) + Chu  ->  the printed P/v :  symbolic residual =", CHU_RESID)
print("   derived:", sp.nsimplify(sp.powsimp(Pv_derived, force=True)))
print("   STRUCTURAL once both are transcribed correctly - it is a TRANSCRIPTION check,")
print("   and it is below check_agreement's ABS_FLOOR of 1e-12, so CI does not compare it.")
print("   The CI-protected version of the same statement is eq6_prefactor_rel_dev below,")
print("   which carries Chu's 35.4 and the 0.44 into a number of order one.\n")

# --- eq. (4) + that P/v -> eq. (6) ---------------------------------------
C6_printed = p("c_eq6")
C6_derived = C4 * CHU_C**0.25
print(f"eq. (6) prefactor:  0.13 * 35.4^(1/4) = {C6_derived:.6f}   printed {C6_printed}   "
      f"rel {C6_derived/C6_printed-1:+.4%}")

# the full j_d that the chain produces, symbolically
jd_chain = sp.powsimp(sp.simplify(
    sp.Rational(13, 100) * ((Pv_printed*mu)/rho**2)**sp.Rational(1, 4) * rho/G), force=True)
Re_sym = d*G/mu
jd_chain_coeff = sp.simplify(sp.powsimp(jd_chain * Re_sym**sp.Rational(36, 100), force=True))
print("\nj_d from the chain :", sp.nsimplify(jd_chain_coeff), "* Re^(-0.36)")
print("   -> the voidage enters as (1-Sigma)^0.36 * Sigma^(-1).")
print("   Eq. (6) AS PRINTED has ((1-Sigma)/Sigma)^(1/3): the SAME exponent on both.")
print("   Sigma^-1 is forced, because P/v carries Sigma^-4 and eq. (4) takes its 4th root.\n")

SIGMA = p("sigma_fallat")
jd_chain_num   = C6_derived * (1-SIGMA)**0.36 / SIGMA          # exponent -0.36
jd_printed_eq6 = C6_printed * ((1-SIGMA)/SIGMA)**(1/3)         # eq. (6) as printed
jd_sigma_inv   = C6_printed * (1-SIGMA)**(1/3) / SIGMA         # eq. (6) with Sigma^-1
C7 = p("c_eq7")
print(f"At Sigma = {SIGMA}:")
print(f"   eq. (6) exactly as printed      -> j_d = {jd_printed_eq6:.4f} Re^(-1/3)   "
      f"({jd_printed_eq6/C7-1:+.2%} on the printed eq. (7) constant {C7})")
print(f"   eq. (6) with Sigma^(-1)         -> j_d = {jd_sigma_inv:.4f} Re^(-1/3)   "
      f"({jd_sigma_inv/C7-1:+.2%})")
print(f"   the chain's own exponent -0.36  -> j_d = {jd_chain_num:.4f} Re^(-0.36)")
print(f"   PRINTED eq. (7)                 -> j_d = {C7:.4f} Re^(-1/3)")
print("\nNEITHER READING OF EQ. (6) REACHES EQ. (7). The Sigma^-1 form gets within "
      f"{abs(jd_sigma_inv/C7-1)*100:.1f} %;\nthe printed form is short by a factor "
      f"{C7/jd_printed_eq6:.2f}. Nothing is repaired here.")

EQ7_FROM_PRINTED = jd_printed_eq6/C7 - 1
EQ7_FROM_SIGINV  = jd_sigma_inv/C7 - 1
'''))

cells.append(code(r'''# What each candidate does against Fallat's OWN measured correlation, which is the
# thing eq. (7) exists to be compared with.
C_FAL, E_FAL, C8 = p("c_fallat"), p("e_fallat"), p("c_eq8")
print("j_d against Fallat's measured  j_d = 0.626 Re^(-0.322)")
print(f"{'Re':>8} {'Fallat':>10} {'eq.(7) 0.69':>13} {'eq.(8) 0.67':>13} "
      f"{'Sigma^-1':>13} {'chain -0.36':>13}")
rows = []
for Re in (100.0, 1000.0, 10000.0):
    fal = C_FAL * Re**E_FAL
    e7  = C7 * Re**(-1/3)
    e8  = C8 * Re**(-1/3)
    si  = jd_sigma_inv * Re**(-1/3)
    ch  = jd_chain_num * Re**(-0.36)
    rows.append((Re, fal, e7, e8, si, ch))
    print(f"{Re:8.0f} {fal:10.5f} {e7:8.5f}{e7/fal-1:+7.1%} {e8:8.5f}{e8/fal-1:+7.1%} "
          f"{si:8.5f}{si/fal-1:+7.1%} {ch:8.5f}{ch/fal-1:+7.1%}")
FAL_DEV_SIGINV = max(abs(r[4]/r[1]-1) for r in rows)
FAL_DEV_EQ7    = max(abs(r[2]/r[1]-1) for r in rows)
FAL_DEV_CHAIN  = max(abs(r[5]/r[1]-1) for r in rows)
print(f"\nWorst deviation from Fallat over the three decades:")
print(f"   printed eq. (7)        {FAL_DEV_EQ7:.1%}   <- the paper calls this agreement "
      f"'substantially complete'")
print(f"   eq. (6) with Sigma^-1  {FAL_DEV_SIGINV:.1%}")
print(f"   the chain, exp -0.36   {FAL_DEV_CHAIN:.1%}")
print("So eq. (7) IS close to Fallat -- but the printed eq. (6) does not produce eq. (7),")
print("and the reading that does is a little further from Fallat than eq. (7) is.")
'''))

cells.append(md(r"""### 4. The pipe chain, and the claim about it

Same machinery, different closure. **This is the control on the packed-bed
result**: if the pipe chain reproduces exactly, then eq. (4)'s constant, the
definition of $j_d$ and the flow-work argument are all sound, and whatever is
wrong upstream is localised to the voidage."""))

cells.append(code(r'''from scipy.optimize import brentq

C_BLA = p("c_blasius")
C9_printed, E9_printed = p("c_eq9"), p("e_eq9")
C9_derived = C4 * C_BLA**0.25
E9_exact = -(1/16 + 1/4)
print(f"eq. (9) prefactor : 0.13 * 0.0396^(1/4) = {C9_derived:.6f}   printed {C9_printed}"
      f"   rel {C9_derived/C9_printed-1:+.4%}")
print(f"eq. (9) exponent  : -(1/16 + 1/4) = {E9_exact:.4f}   printed {E9_printed}")
print("THE PIPE CHAIN REPRODUCES. That is what localises the packed-bed problem.\n")

C_CC, E_CC = p("c_cc"), p("e_cc")
RE_LO, RE_HI = p("re_lo_cc"), p("re_hi_cc")
RE_CLAIM, MAXDEV_CLAIM = p("re_agree_cc"), p("maxdev_cc")/100.0
print(f"The paper: eq. (9) 'agrees with' Chilton-Colburn at N_Re = {RE_CLAIM:.0f},")
print(f"           with a maximum deviation of +/- {MAXDEV_CLAIM:.0%} over "
      f"{RE_LO:.0f}-{RE_HI:.0f}.\n")
pipe = {}
for exp9, lab in ((E9_printed, "printed -0.31"), (E9_exact, "exact -0.3125")):
    dev = lambda Re: C9_derived*Re**exp9/(C_CC*Re**E_CC) - 1
    Rex = brentq(dev, 10.0, 1e9)
    d_lo, d_hi, d_claim = dev(RE_LO), dev(RE_HI), dev(RE_CLAIM)
    pipe[lab] = (Rex, max(abs(d_lo), abs(d_hi)), d_claim)
    print(f"   {lab}: they cross at N_Re = {Rex:.0f}")
    print(f"        deviation at {RE_LO:.0f}: {d_lo:+.2%}   at {RE_HI:.0f}: {d_hi:+.2%}"
          f"   at the claimed {RE_CLAIM:.0f}: {d_claim:+.2%}")
    print(f"        max |deviation| over the stated range: {max(abs(d_lo),abs(d_hi)):.2%}")
PIPE_CROSS_PRINTED = pipe["printed -0.31"][0]
PIPE_MAXDEV_PRINTED = pipe["printed -0.31"][1]
PIPE_DEV_AT_CLAIM = pipe["printed -0.31"][2]
dev9 = lambda Re: C9_derived*Re**E9_printed/(C_CC*Re**E_CC) - 1
inv9 = lambda Re: 1.0/(1.0 + dev9(Re)) - 1.0   # the paper's own direction, (CC - eq9)/eq9
PIPE_MAXDEV_PAPER_DIR = max(abs(inv9(RE_LO)), abs(inv9(RE_HI)))
print(f"\nThe crossing is at N_Re ~ {PIPE_CROSS_PRINTED:.0f}, not {RE_CLAIM:.0f} -- a factor "
      f"{RE_CLAIM/PIPE_CROSS_PRINTED:.2f}.\nAt {RE_CLAIM:.0f} the two differ by "
      f"{abs(PIPE_DEV_AT_CLAIM):.1%}, and the worst case over the stated range is "
      f"{PIPE_MAXDEV_PRINTED:.1%}\n(convention (eq9-CC)/CC; in the paper's own direction, "
      f"(CC-eq9)/eq9, it is {PIPE_MAXDEV_PAPER_DIR:.1%}),\nnot {MAXDEV_CLAIM:.0%}. The deviation "
      f"takes BOTH signs inside the stated range ({dev9(RE_LO):+.1%} at\n{RE_LO:.0f} to "
      f"{dev9(RE_HI):+.1%} at {RE_HI:.0f}, crossing zero at N_Re = {PIPE_CROSS_PRINTED:.0f}), so "
      "the paper's '+/-' is\nright about the shape; the size and the crossing are what its "
      "sentence gets wrong.\nReported, not repaired.")
'''))

cells.append(md(r"""### 5. Two more identities the paper's constants must satisfy

Eq. (10) — the power needed to just suspend solids — follows from equating
eqs. (1) and (4). The paper prints its exponent structure and leaves the constant
as `const.`; the constant is recovered here. And the "more precise form" of
eq. (1) on page 44 must be eq. (1) itself."""))

cells.append(code(r'''drho, g_, D_, Pv_ = sp.symbols("Delta_rho g D P_over_v", positive=True)
# eq. (1) and eq. (4), both as k_L Sc^(2/3):
lhs1 = sp.Rational(31, 100) * (drho*mu*g_/rho**2)**sp.Rational(1, 3)
lhs4 = sp.Rational(13, 100) * (Pv_*mu/rho**2)**sp.Rational(1, 4)
Pv_sol = sp.solve(sp.Eq(lhs1, lhs4), Pv_)[0]
Pv_sol = sp.powsimp(sp.simplify(Pv_sol), force=True)
print("Equating eq. (1) and eq. (4) for P/v gives")
print("   P/v =", Pv_sol)
# printed eq. (10):  P/v = const (g drho)^(4/3) mu^(1/3) / rho^(2/3)
Pv_10_shape = (g_*drho)**sp.Rational(4, 3) * mu**sp.Rational(1, 3) / rho**sp.Rational(2, 3)
_ratio = sp.simplify(sp.powsimp(Pv_sol / Pv_10_shape, force=True))
EQ10_RESID = float(max(abs(complex(_ratio.subs(sub)) / complex(_ratio.subs(
    {drho: 1, mu: 1, g_: 1, rho: 1})) - 1)
    for sub in ({drho: 2, mu: 3, g_: 5, rho: 7},
                {drho: 11, mu: 13, g_: 17, rho: 19},
                {drho: sp.Rational(1, 3), mu: 29, g_: sp.Rational(2, 7), rho: 31})))
EQ10_CONST = float(_ratio.subs({drho: 1, mu: 1, g_: 1, rho: 1}))
print(f"   printed eq. (10) form: P/v = const (g drho)^(4/3) mu^(1/3) rho^(-2/3)")
print(f"   the ratio to that shape is parameter-free to {EQ10_RESID:.3e} over three")
print(f"   decades of every symbol  ->  the exponent structure matches exactly.")
print(f"   the constant the paper leaves as 'const.' is (0.31/0.13)^4 = {EQ10_CONST:.4f}")
print(f"   cross-check, computed the other way: {(C1/C4)**4:.4f}")
print("   STRUCTURAL and below ABS_FLOOR; the CI-carried number is eq10_implied_const.\n")

# --- eq. (1) <-> Sh = 2.0 + 0.31 Ra^(1/3) --------------------------------
# Sh = k d/D ; Ra = Gr Sc = d^3 drho g /(mu D)
d_, Sc_ = sp.symbols("d N_Sc", positive=True)
kL_from_eq1 = sp.Rational(31, 100) * (drho*mu*g_/rho**2)**sp.Rational(1, 3) \
              * (mu/(rho*D_))**sp.Rational(-2, 3)
Sh_from_eq1 = sp.powsimp(sp.simplify(kL_from_eq1*d_/D_), force=True)
Sh_from_p44 = sp.Rational(31, 100) * (d_**3*drho*g_/(mu*D_))**sp.Rational(1, 3)
RA_RESID = float(sp.simplify(Sh_from_eq1/Sh_from_p44 - 1))
print("eq. (1) rewritten as a Sherwood number :", sp.nsimplify(Sh_from_eq1))
print("page 44's  0.31 N_Ra^(1/3)             :", sp.nsimplify(Sh_from_p44))
print("   residual =", RA_RESID, " -> the two printings agree, exponent for exponent.")
print("   STRUCTURAL once both are read correctly; it is a check on the TRANSCRIPTION")
print("   of two separately typeset equations, and the break table exercises it.")
'''))

cells.append(md(r"""### 6. pymrm: the stagnant floor, and the two exponents

The `2.0` of the page-44 form is the conduction limit of a sphere. It is
recovered from the spherical operator, and then the boundary layer is solved on
both sides of the split."""))

cells.append(code(r'''print("Sh for a sphere in a stagnant medium (nu = 2, spherical):")
ratios = np.array([1e-1, 1e-2, 1e-3, 1e-4])
sh = np.array([sh_sphere(r, n_r=800) for r in ratios])
for r, s in zip(ratios, sh):
    print(f"   R/R_inf = {r:7.0e}   Sh = {s:.8f}   exact 2/(1-R/R_inf) = {2/(1-r):.8f}"
          f"   rel {s/(2/(1-r))-1:+.2e}")
# extrapolate R/R_inf -> 0 : Sh = 2/(1-r) = 2(1 + r + ...), so Sh*(1-r) -> 2
SH_STAGNANT = float(sh[-1] * (1 - ratios[-1]))
SH_PRINTED = p("sh_floor")
print(f"\n   Sh(R/R_inf -> 0) = {SH_STAGNANT:.8f}   printed floor {SH_PRINTED}"
      f"   rel {SH_STAGNANT/SH_PRINTED-1:+.3e}")
print(f"   wrong geometry, nu = 0 (Cartesian)   -> Sh = {sh_sphere(1e-4, nu=0):.6f}")
print(f"   wrong geometry, nu = 1 (cylindrical) -> Sh = {sh_sphere(1e-4, nu=1):.6f}")
print("   so this number DOES see the geometry, which is the point of computing it.")
'''))

cells.append(code(r'''D_REF, L_REF, GDOT = 1e-9, 1.0, 1.0      # SI; only ratios matter

print("Boundary layer, both limits, against closed forms that share no code with pymrm:")
conv = {}
for lab, kw, exact in (("rigid  (u_s = 0, Leveque)",  dict(u_s=0.0, gdot=GDOT),
                        leveque(D_REF, GDOT, L_REF)),
                       ("mobile (gdot = 0, Higbie)",  dict(u_s=1e-3, gdot=0.0),
                        higbie(D_REF, 1e-3, L_REF))):
    ks = []
    for n in (200, 400, 800):
        kb, kg, _, _ = BoundaryLayer(D=D_REF, L=L_REF, n_y=n, n_x=n, **kw).run()
        ks.append(kb)
        print(f"   {lab:28s} n = {n:4d}  kbar = {kb:.8e}  rel {kb/exact-1:+.3e}")
    # first order in the x-march (backward Euler): Richardson on the last pair
    k_ext = 2*ks[-1] - ks[-2]
    order = np.log(abs(ks[0]-ks[1])/abs(ks[1]-ks[2]))/np.log(2.0)
    conv[lab] = (ks[-1]/exact - 1, k_ext/exact - 1, order)
    print(f"   {'':28s} observed order {order:.3f};  extrapolated rel "
          f"{k_ext/exact-1:+.3e}\n")
LEV_RAW, LEV_EXT, LEV_ORD = conv["rigid  (u_s = 0, Leveque)"]
HIG_RAW, HIG_EXT, HIG_ORD = conv["mobile (gdot = 0, Higbie)"]
'''))

cells.append(code(r'''def d_ln_k_d_ln_D(u_s, D0=D_REF, n=400, rel=0.15, gdot=GDOT, L=L_REF, **kw):
    """The exponent the paper reports: d ln k_L / d ln D_L, from solved fields."""
    Ds = D0*np.array([1/(1+rel), (1+rel)])
    ks = [BoundaryLayer(D=D, u_s=u_s, gdot=gdot, L=L, n_y=n, n_x=n, **kw).run()[0]
          for D in Ds]
    return float(np.log(ks[1]/ks[0]) / np.log(Ds[1]/Ds[0]))

d_lev = (9*D_REF*L_REF/GDOT)**(1/3)      # the rigid-branch length scale
def u_s_of(lam):   return lam * GDOT * d_lev

EXP_RIGID  = d_ln_k_d_ln_D(0.0)
EXP_MOBILE = d_ln_k_d_ln_D(u_s_of(1e4))
print(f"rigid  branch: d ln k / d ln D = {EXP_RIGID:.6f}   paper prints 2/3 = {2/3:.6f}"
      f"   rel {EXP_RIGID/(2/3)-1:+.2e}")
print(f"mobile branch: d ln k / d ln D = {EXP_MOBILE:.6f}   paper prints 1/2 = {0.5:.6f}"
      f"   rel {EXP_MOBILE/0.5-1:+.2e}")
print("\nBOTH LIMITS ARE FIXED BY DIMENSIONAL ANALYSIS once the velocity profile is set --")
print("there is only one length scale in each -- so these two numbers test the")
print("IMPLEMENTATION, not the physics. What they can catch is the wrong profile, and")
print("the break table shows they do: swapping the two profiles swaps the two exponents.")
'''))

cells.append(code(r'''LAMS = np.array([0.0, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 10.0, 30.0, 100.0])
EXPS = np.array([d_ln_k_d_ln_D(u_s_of(l)) for l in LAMS])
print("Lambda = u_s / (gdot * delta_Leveque)   ->   d ln k_L / d ln D_L")
for l, e in zip(LAMS, EXPS):
    print(f"   {l:8.3g}   {e:.6f}")

f_half = lambda l: d_ln_k_d_ln_D(u_s_of(l)) - 7/12      # halfway between 1/2 and 2/3
LAM_HALF = brentq(f_half, 0.05, 5.0, xtol=1e-4)
f_90 = lambda l: d_ln_k_d_ln_D(u_s_of(l)) - (0.5 + 0.1*(2/3 - 0.5))
LAM_90 = brentq(f_90, 0.5, 100.0, xtol=1e-3)
print(f"\nhalfway  (n = 7/12 = {7/12:.4f}) at Lambda = {LAM_HALF:.4f}")
print(f"within 10 % of the mobile value (n = {0.5+0.1*(2/3-0.5):.4f}) at Lambda = {LAM_90:.3f}")

fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.semilogx(np.where(LAMS > 0, LAMS, 1e-3), EXPS, "o-", lw=1.6, ms=4.5)
ax.axhline(2/3, color="tab:blue", ls="--", lw=1)
ax.axhline(0.5, color="tab:orange", ls="--", lw=1)
ax.axvline(LAM_HALF, color="0.4", ls=":", lw=1)
ax.text(1.2e-3, 2/3+0.004, r"eq. (1):  $k_L\propto D^{2/3}$  (rigid, Froessling)",
        color="tab:blue", fontsize=9)
ax.text(1.2e-3, 0.5+0.004, r"eq. (2):  $k_L\propto D^{1/2}$  (mobile, Higbie)",
        color="tab:orange", fontsize=9)
ax.annotate(rf"$\Lambda_{{1/2}}={LAM_HALF:.3f}$", (LAM_HALF, 7/12),
            xytext=(20, -22), textcoords="offset points", fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8, color="0.4"))
ax.set_xlabel(r"$\Lambda = u_s\,/\,(\dot\gamma\,\delta_{\rm Lev})$   "
              r"(interface mobility; 0 = rigid, $\infty$ = free)")
ax.set_ylabel(r"$\mathrm{d}\ln k_L/\mathrm{d}\ln D_L$")
ax.set_title("The transition the paper describes and never models")
ax.grid(alpha=0.3, which="both")
fig.tight_layout(); plt.show()
'''))

cells.append(md(r"""### 7. Eq. (2) against the theory it says it represents

The note added in proof says eq. (2) "represents the data of BAIRD obtained with
spherical cap bubbles of equivalent diameters up to 5 cm **as well as the
theoretical equation proposed by this worker**",
$K_L = 0.975\,D_L^{1/2}(g/d)^{1/4}$.

Those two cannot both be right everywhere, because **eq. (2) contains no bubble
diameter and Baird's equation does**. They agree on exactly one diameter for any
given liquid, and that diameter is computable from the two printed constants."""))

cells.append(code(r'''C_BAIRD = p("c_baird")
def d_agree(drho, mu, rho):
    """Diameter at which eq. (2) and Baird's theoretical equation coincide (cm)."""
    lhs = C2 * (drho*mu*G_CGS/rho**2)**(1/3) * (rho/mu)**0.5   # eq. (2) / D^(1/2)
    return G_CGS * (C_BAIRD/lhs)**4

r = RV.loc["this_work_large"]
print("Diameter at which eq. (2) and Baird's K_L = 0.975 D^(1/2)(g/d)^(1/4) coincide:")
cases = [("mu = 0.6 cP  (row minimum)", r.drho_lo, r.mu_c_lo/100.0, r.rho_c_hi),
         ("mu = 1 cP    (water)",       1.00,      0.010,           1.000),
         ("mu = 87 cP   (row maximum)", r.drho_hi, r.mu_c_hi/100.0, r.rho_c_hi)]
d_cross = {}
for lab, dr, m, rh in cases:
    dc = d_agree(dr, m, rh)
    d_cross[lab] = dc
    where = ("INSIDE" if r.d_lo <= dc <= r.d_hi else
             "below" if dc < r.d_lo else "above")
    print(f"   {lab:28s} d* = {dc:8.4f} cm = {dc*10:7.3f} mm   ({where} the row's "
          f"{r.d_lo}-{r.d_hi} cm)")
D_STAR_WATER = d_cross["mu = 1 cP    (water)"]
print(f"\nNot one of the three crossings falls inside the diameters the large-bubble runs")
print("actually covered: at low viscosity the crossing is below the smallest bubble and")
print("at high viscosity far above the largest. So eq. (2) and Baird's equation do not")
print("coincide anywhere in the tested window; they only bracket it.")

# how far apart they are at the extremes of the row's own diameter range
print("\nRatio  eq.(2) / Baird  across the row's own diameters (D_L = 1.9e-6 cm2/s):")
for lab, dr, m, rh in cases:
    for dd in (r.d_lo, r.d_hi, p("d_baird_max")):
        DL = r.D_L_lo_mant * 10.0**int(r.D_L_exp_printed)
        k2 = kL_eq2(dr, m, rh, DL)
        kb = C_BAIRD * DL**0.5 * (G_CGS/dd)**0.25
        print(f"   {lab:28s} d = {dd:5.2f} cm  ->  {k2/kb:6.3f}")
BAIRD_RATIO_5CM = (kL_eq2(1.00, 0.010, 1.000, 1.9e-6)
                   / (C_BAIRD * (1.9e-6)**0.5 * (G_CGS/p("d_baird_max"))**0.25))
print(f"\nAt Baird's largest bubble (5 cm) in water, eq. (2) stands "
      f"{BAIRD_RATIO_5CM:.2f}x above his theory.")
print("The claim that eq. (2) 'represents' it is therefore a statement about a window,")
print("not an identity, and the paper does not say which window.")
'''))

# --------------------------------------------------------------- validation
cells.append(md(r"""## Validation

### What is checked, and against what

| # | Check | Route A | Route B | Can it fail? |
|---|---|---|---|---|
| 1 | eq. (5) + Chu → the printed $P/v$ | sympy on the printed Chu form | the paper's own printed $P/v$ | yes — a mis-read exponent breaks it; **structural once both are right**, and below `ABS_FLOOR` |
| 2 | eq. (6) prefactor | $0.13\times35.4^{1/4}$ | printed 0.318 | yes |
| 3 | eq. (9) prefactor and exponent | $0.13\times0.0396^{1/4}$, $-1/16-1/4$ | printed 0.058, $-0.31$ | yes |
| 4 | eq. (7) from eq. (6) | two readings of the voidage group | printed 0.69 | **it does fail** |
| 5 | eq. (9) vs Chilton–Colburn | brentq on the two power laws | the printed claim | **it does fail** |
| 6 | 0.28 | Froessling × √Allen | printed 0.28 | yes |
| 7 | 0.34 | Friedlander × $18^{-1/3}$ | printed 0.34 | yes |
| 8 | $N_{Sh}\to 2$ | pymrm, `nu=2`, extrapolated | printed 2.0 | yes — `nu=0` gives 2.0e-4 |
| 9 | Lévêque magnitude | pymrm march | similarity closed form | yes |
| 10 | Higbie magnitude | pymrm march | penetration closed form | yes |
| 11 | exponents 2/3, 1/2 | pymrm, from solved fields | printed exponents | **only against a wrong profile** — fixed by dimensional analysis otherwise, and labelled |
| 12 | eq. (1) ⇄ page-44 Rayleigh form | sympy | sympy | yes on transcription; structural otherwise |
| 13 | eq. (10) exponents | sympy on (1) and (4) | printed eq. (10) | yes on transcription; structural otherwise |

### What perturbation testing cannot detect

**A break row establishes sensitivity; it never establishes correctness.** Every
row below perturbs an input and watches a number move, so no row can catch a
baseline that is *wrong* rather than *insensitive* — a discretisation error, a
metric read at the wrong place, or a coincidence. That is why checks 6, 7, 9 and
10 exist: each reaches a printed number by a route that shares no code with the
one it is compared against.

Specifically, **the break table cannot see**:

- a systematic error in the $x$-march, because it cancels out of a *ratio*
  between a perturbed and an unperturbed run. The convergence study and the
  Richardson extrapolation above are what bound it, and the closed-form
  comparisons are what confirm the extrapolated value;
- whether $\Lambda_{1/2}$ is being read at the right place. That is the `A2.6`
  failure mode, and it is covered by a dedicated row that computes the same
  transition from the **local** coefficient at $x=L$ instead of the mean;
- whether the fitted constants 0.31, 0.42 and 0.13 are any good. **Nothing on
  this page can see that**, because the evidence is in Figs 1, 2 and 6.

### Grid, march and domain independence"""))

cells.append(code(r'''print("Convergence of the boundary layer, all three knobs, on the rigid branch:")
base = dict(D=D_REF, u_s=0.0, gdot=GDOT, L=L_REF)
ex_lev = leveque(D_REF, GDOT, L_REF)
print(f"{'n_y':>6}{'n_x':>7}{'y_span':>9}{'stretch':>9}{'kbar':>16}{'rel err':>12}")
grid_rows = []
for n_y, n_x, span, st in ((400, 400, 16.0, 6.0), (800, 400, 16.0, 6.0),
                           (1600, 400, 16.0, 6.0), (400, 800, 16.0, 6.0),
                           (400, 1600, 16.0, 6.0), (400, 400, 24.0, 6.0),
                           (400, 400, 32.0, 6.0), (400, 400, 16.0, 8.0)):
    kb = BoundaryLayer(n_y=n_y, n_x=n_x, y_span=span, stretch=st, **base).run()[0]
    grid_rows.append((n_y, n_x, span, st, kb))
    print(f"{n_y:6d}{n_x:7d}{span:9.1f}{st:9.1f}{kb:16.8e}{kb/ex_lev-1:+12.3e}")
ny_spread = max(abs(r[4]/ex_lev-1) for r in grid_rows[:3]) - min(abs(r[4]/ex_lev-1) for r in grid_rows[:3])
print(f"\nRefining n_y alone moves the answer by {ny_spread:.2e} in relative error, while")
print("refining n_x moves it by an order of magnitude more: THE MARCH STEP IS THE ERROR,")
print("and it is what the Richardson extrapolation above removes. The y grid, the domain")
print("height and the stretching are all already converged at the production setting.")
DOMAIN_SENS = abs(grid_rows[6][4]/grid_rows[0][4] - 1)
print(f"Doubling the far-field height moves kbar by {DOMAIN_SENS:.2e} RELATIVE, i.e. it")
print("is already outside the layer at the production setting.")
'''))

cells.append(md(r"""### The break table

Every metric in `agreement.json` needs a row that moves it. Rows that move
nothing are kept and **labelled**, because an identity is worth having once it is
named as one."""))

cells.append(code(r"""lo1, hi1 = corner_window("this_work_small", kL_eq1, D_EXP_SMALL)
M = {
    # --- the paper's own algebra ---
    "chu_identity_resid":              abs(CHU_RESID),
    "eq6_prefactor_rel_dev":           C6_derived/C6_printed - 1,
    "eq9_prefactor_rel_dev":           C9_derived/C9_printed - 1,
    "eq7_from_printed_eq6_rel_dev":    EQ7_FROM_PRINTED,
    "eq7_from_sigma_inverse_rel_dev":  EQ7_FROM_SIGINV,
    "fallat_dev_sigma_inverse":        FAL_DEV_SIGINV,
    "fallat_dev_printed_eq7":          FAL_DEV_EQ7,
    "pipe_cross_Re":                   PIPE_CROSS_PRINTED,
    "pipe_dev_at_claimed_Re":          PIPE_DEV_AT_CLAIM,
    "pipe_maxdev_over_range":          PIPE_MAXDEV_PRINTED,
    "eq10_implied_const":              EQ10_CONST,
    "eq10_exponent_resid":             abs(EQ10_RESID),
    "ra_form_resid":                   abs(RA_RESID),
    # --- the two theoretical derivations ---
    "c28_rel_dev":                     c28/c28_printed - 1,
    "c34_rel_dev":                     c34/c34_printed - 1,
    "froessling_implied_by_paper":     c28_printed/np.sqrt(C_ALLEN),
    "allen_quarter_from_drag_law":     C_ALLEN_FROM_CD,
    # --- the split ---
    "crossover_Sc":                    SC_STAR,
    "eq2_window_hi_over_printed":      hi2/obs2[1],
    "eq2_window_lo_over_printed":      lo2/obs2[0],
    "eq1_window_hi_over_printed":      R_HI_SIX,
    "eq1_window_hi_printed_exponent":  R_HI_PRINTED_EXP,
    "baird_cross_d_water_cm":          D_STAR_WATER,
    "baird_ratio_at_5cm":              BAIRD_RATIO_5CM,
    # --- pymrm ---
    "sh_stagnant_extrap":              SH_STAGNANT,
    "sh_stagnant_rel_dev":             SH_STAGNANT/SH_PRINTED - 1,
    "leveque_rel_err_raw":             LEV_RAW,
    "leveque_rel_err_extrap":          LEV_EXT,
    "leveque_observed_order":          LEV_ORD,
    "higbie_rel_err_raw":              HIG_RAW,
    "higbie_rel_err_extrap":           HIG_EXT,
    "higbie_observed_order":           HIG_ORD,
    "exp_rigid":                       EXP_RIGID,
    "exp_mobile":                      EXP_MOBILE,
    "exp_split_gap":                   EXP_RIGID - EXP_MOBILE,
    "exp_mobile_rel_dev":              EXP_MOBILE/0.5 - 1,
    "lambda_half":                     LAM_HALF,
    "lambda_within_10pct_mobile":      LAM_90,
    "bl_domain_sensitivity":           DOMAIN_SENS,
}
print(len(M), "metrics defined. The break table below must move every one of them.")
"""))

cells.append(md(r"""### The break table

Every metric in `agreement.json` needs a row that moves it, so the table is
organised by **injected defect**, and each defect recomputes every metric it can
reach. A coverage check at the end names any metric no defect touches.

Rows that move nothing are kept and **labelled**, because an identity is worth
having once it is named as one."""))

cells.append(code(r"""# ---- helpers that recompute a metric under an injected defect --------------
def _bl_metrics(ns=(200, 400, 800), **kw):
    # re-run the boundary layer under an injected defect, at three march steps.
    # NOTE the pops happen ONCE, before the comprehension: doing them inside it
    # would apply the defect to the first grid only.
    u_s = kw.pop("u_s", 0.0)
    gdot = kw.pop("gdot", GDOT)
    face = kw.pop("interface", "dirichlet")
    return [BoundaryLayer(D=D_REF, u_s=u_s, gdot=gdot, L=L_REF,
                          n_y=n, n_x=n, **kw).run(interface=face)[0] for n in ns]

def _lev_set(ks, exact):
    ext = 2*ks[-1] - ks[-2]
    den = abs(ks[1]-ks[2])
    order = np.log(abs(ks[0]-ks[1])/den)/np.log(2.0) if den > 0 else np.nan
    return ks[-1]/exact - 1, ext/exact - 1, order

def _pipe(c9=None, e9=None, ccc=None, cce=None):
    c9, e9 = (C9_derived if c9 is None else c9), (E9_printed if e9 is None else e9)
    ccc, cce = (C_CC if ccc is None else ccc), (E_CC if cce is None else cce)
    dev = lambda Re: c9*Re**e9/(ccc*Re**cce) - 1
    try:
        Rex = brentq(dev, 1e-6, 1e12)
    except ValueError:
        Rex = np.nan
    return {"pipe_cross_Re": Rex,
            "pipe_dev_at_claimed_Re": dev(RE_CLAIM),
            "pipe_maxdev_over_range": max(abs(dev(RE_LO)), abs(dev(RE_HI)))}

def _fallat(c_eq6=None, sig=None, c_fal=None, e_fal=None, c7=None):
    c_eq6 = C6_printed if c_eq6 is None else c_eq6
    sig   = SIGMA      if sig   is None else sig
    c_fal = C_FAL      if c_fal is None else c_fal
    e_fal = E_FAL      if e_fal is None else e_fal
    c7    = C7         if c7    is None else c7
    printed = c_eq6*((1-sig)/sig)**(1/3)
    siginv  = c_eq6*(1-sig)**(1/3)/sig
    Res = (100.0, 1000.0, 10000.0)
    return {"eq7_from_printed_eq6_rel_dev": printed/c7 - 1,
            "eq7_from_sigma_inverse_rel_dev": siginv/c7 - 1,
            "fallat_dev_sigma_inverse":
                max(abs(siginv*R**(-1/3)/(c_fal*R**e_fal) - 1) for R in Res),
            "fallat_dev_printed_eq7":
                max(abs(c7*R**(-1/3)/(c_fal*R**e_fal) - 1) for R in Res)}

def _windows(c1=None, c2=None, d_exp=None):
    c1 = C1 if c1 is None else c1
    c2 = C2 if c2 is None else c2
    d_exp = D_EXP_SMALL if d_exp is None else d_exp
    f1 = lambda dr, mu, rho, D: c1*(dr*mu*G_CGS/rho**2)**(1/3)*(mu/(rho*D))**(-2/3)
    f2 = lambda dr, mu, rho, D: c2*(dr*mu*G_CGS/rho**2)**(1/3)*(mu/(rho*D))**(-1/2)
    a, b = corner_window("this_work_large", f2, -6)
    c, d = corner_window("this_work_small", f1, d_exp)
    e, g = corner_window("this_work_small", f1, -4)
    return {"eq2_window_hi_over_printed": b/obs2[1],
            "eq2_window_lo_over_printed": a/obs2[0],
            "eq1_window_hi_over_printed": d/RV.loc["this_work_small", "kL_hi"],
            "eq1_window_hi_printed_exponent": g/RV.loc["this_work_small", "kL_hi"]}

def _baird(cb=None, c2=None):
    cb = C_BAIRD if cb is None else cb
    c2 = C2 if c2 is None else c2
    lhs = c2*(1.00*0.010*G_CGS/1.000**2)**(1/3)*(1.000/0.010)**0.5
    DL = 1.9e-6
    k2 = c2*(1.00*0.010*G_CGS/1.000**2)**(1/3)*(0.010/(1.000*DL))**(-1/2)
    return {"baird_cross_d_water_cm": G_CGS*(cb/lhs)**4,
            "baird_ratio_at_5cm": k2/(cb*DL**0.5*(G_CGS/p("d_baird_max"))**0.25)}

def _chu_resid(a_chu=44):
    # symbolic residual of eq. (5)+Chu against the printed P/v, with Chu's
    # exponent perturbed: non-zero the moment the transcription is wrong
    dp = sp.Rational(354, 10)*(mu*(1-S)/(d*G))**sp.Rational(a_chu, 100) \
         * G**2*(1-S)/(d*rho*S**3)
    return abs(float(sp.simplify(sp.simplify(dp*G/(rho*S))/Pv_printed - 1)
                     .subs({mu: 2, S: sp.Rational(1, 3), G: 5, d: 7, rho: 11})))

def _eq10_resid(quarter=sp.Rational(1, 4)):
    l4 = sp.Rational(13, 100)*(Pv_*mu/rho**2)**quarter
    sol = sp.powsimp(sp.simplify(sp.solve(sp.Eq(lhs1, l4), Pv_)[0]), force=True)
    rat = sp.simplify(sp.powsimp(sol/Pv_10_shape, force=True))
    ref = complex(rat.subs({drho: 1, mu: 1, g_: 1, rho: 1}))
    return float(max(abs(complex(rat.subs(sub))/ref - 1)
                     for sub in ({drho: 2, mu: 3, g_: 5, rho: 7},
                                 {drho: 11, mu: 13, g_: 17, rho: 19})))

def _ra_resid(sc_exp=sp.Rational(2, 3)):
    k = sp.Rational(31, 100)*(drho*mu*g_/rho**2)**sp.Rational(1, 3)*(mu/(rho*D_))**(-sc_exp)
    rat = sp.simplify(sp.powsimp(k*d_/D_/Sh_from_p44, force=True))
    ref = complex(rat.subs({drho: 1, mu: 1, g_: 1, rho: 1, D_: 1, d_: 1}))
    return float(max(abs(complex(rat.subs(sub)) - ref)
                     for sub in ({drho: 2, mu: 3, g_: 5, rho: 7, D_: 2, d_: 3},)))

# ---- the defects -----------------------------------------------------------
DEFECTS = []
def defect(name, moves, note=""):
    DEFECTS.append((name, moves, note))

defect("Chu's constant 35.4 read as 354",
       {"eq6_prefactor_rel_dev": C4*354.0**0.25/C6_printed - 1})
defect("Chu's exponent 0.44 read as 0.50",
       {"chu_identity_resid": _chu_resid(50)})
defect("Blasius 0.0396 read as 0.396",
       {"eq9_prefactor_rel_dev": C4*0.396**0.25/C9_printed - 1,
        **_pipe(c9=C4*0.396**0.25)})
defect("eq. (4) constant 0.13 read as 0.18",
       {"eq6_prefactor_rel_dev": 0.18*CHU_C**0.25/C6_printed - 1,
        "eq9_prefactor_rel_dev": 0.18*C_BLA**0.25/C9_printed - 1,
        "eq10_implied_const": (C1/0.18)**4,
        **_pipe(c9=0.18*C_BLA**0.25)})
defect("eq. (4) exponent 1/4 read as 1/3",
       {"eq10_exponent_resid": _eq10_resid(sp.Rational(1, 3))})
defect("eq. (1) Schmidt exponent 2/3 read as 1/2",
       {"ra_form_resid": _ra_resid(sp.Rational(1, 2)),
        "crossover_Sc": np.inf})
defect("eqs (1) and (2) constants swapped",
       {"crossover_Sc": (C1/C2)**6, **_windows(c1=C2, c2=C1)})
defect("Allen's 1/4 read as 1/2",
       {"c28_rel_dev": FROESSLING*np.sqrt(0.5)/c28_printed - 1,
        "froessling_implied_by_paper": c28_printed/np.sqrt(0.5)})
defect("Allen's drag law C_D = 10 Re^-1/2 read as 24 Re^-1/2",
       {"allen_quarter_from_drag_law": (2.0/(3.0*24.0/2.0))**(2/3)})
defect("Stokes' 18 read as 24",
       {"c34_rel_dev": C_FRIED*(1/24.0)**(1/3)/c34_printed - 1})
defect("Friedlander's 0.89 read as 0.98",
       {"c34_rel_dev": 0.98*(1/STOKES_DEN)**(1/3)/c34_printed - 1})
defect("Fallat's voidage 0.41 read as 0.30", _fallat(sig=0.30),
       "the +0.0% on fallat_dev_printed_eq7 is OUT OF THIS DEFECT'S REACH - that metric reads only eq. (7) and Fallat's own constants, neither of which this row touches")
defect("eq. (6) constant 0.318 read as 0.618", _fallat(c_eq6=0.618),
       "the +0.0% on fallat_dev_printed_eq7 is OUT OF THIS DEFECT'S REACH - eq. (7) as printed does not read eq. (6)'s constant")
defect("Fallat's exponent -0.322 read as -0.222", _fallat(e_fal=-0.222),
       "the two +0.0% eq.(7)-derivation metrics are OUT OF THIS DEFECT'S REACH - neither reads Fallat's exponent")
defect("eq. (9) exponent -0.31 read as -0.25", _pipe(e9=-0.25))
defect("Chilton-Colburn 0.023 read as 0.032", _pipe(ccc=0.032))
defect("Baird's 0.975 read as 0.75", _baird(cb=0.75))
defect("small-bubble D exponent left at the printed -4",
       {"eq1_window_hi_over_printed": R_HI_PRINTED_EXP})
defect("sphere solved with nu = 0 (Cartesian)",
       {"sh_stagnant_extrap": sh_sphere(1e-4, nu=0)*(1-1e-4),
        "sh_stagnant_rel_dev": sh_sphere(1e-4, nu=0)*(1-1e-4)/SH_PRINTED - 1})
defect("sphere solved with nu = 1 (cylindrical)",
       {"sh_stagnant_extrap": sh_sphere(1e-4, nu=1)*(1-1e-4),
        "sh_stagnant_rel_dev": sh_sphere(1e-4, nu=1)*(1-1e-4)/SH_PRINTED - 1})
defect("sphere outer boundary closed (no flux)",
       {"sh_stagnant_extrap": sh_sphere(1e-4, outer="noflux")*(1-1e-4),
        "sh_stagnant_rel_dev": sh_sphere(1e-4, outer="noflux")*(1-1e-4)/SH_PRINTED - 1})

_raw, _ext, _ord = _lev_set(_bl_metrics(u_s=u_s_of(1e4), gdot=0.0), ex_lev)
defect("rigid branch run with the MOBILE profile",
       {"leveque_rel_err_raw": _raw, "leveque_rel_err_extrap": _ext,
        "leveque_observed_order": _ord, "exp_rigid": EXP_MOBILE,
        "exp_split_gap": EXP_MOBILE - EXP_MOBILE})
_exh = higbie(D_REF, 1e-3, L_REF)
_raw, _ext, _ord = _lev_set(_bl_metrics(u_s=0.0, gdot=GDOT), _exh)
defect("mobile branch run with the RIGID profile",
       {"higbie_rel_err_raw": _raw, "higbie_rel_err_extrap": _ext,
        "higbie_observed_order": _ord, "exp_mobile": EXP_RIGID,
        "exp_mobile_rel_dev": EXP_RIGID/0.5 - 1,
        "exp_split_gap": EXP_RIGID - EXP_RIGID})
_raw, _ext, _ord = _lev_set(_bl_metrics(u_s=0.0, gdot=GDOT, y_span=3.0), ex_lev)
defect("far field truncated at 3 delta instead of 16",
       {"leveque_rel_err_raw": _raw, "leveque_rel_err_extrap": _ext,
        "leveque_observed_order": _ord,
        "bl_domain_sensitivity": abs(
            BoundaryLayer(D=D_REF, u_s=0.0, gdot=GDOT, L=L_REF, n_y=400, n_x=400,
                          y_span=6.0).run()[0]
            / BoundaryLayer(D=D_REF, u_s=0.0, gdot=GDOT, L=L_REF, n_y=400, n_x=400,
                            y_span=3.0).run()[0] - 1)})
_raw, _ext, _ord = _lev_set(_bl_metrics(ns=(12, 24, 48)), ex_lev)
_rawh, _exth, _ordh = _lev_set(_bl_metrics(ns=(12, 24, 48), u_s=1e-3, gdot=0.0), _exh)
defect("order measured at n = 12/24/48, outside the asymptotic range",
       {"leveque_observed_order": _ord, "higbie_observed_order": _ordh,
        "leveque_rel_err_raw": _raw, "leveque_rel_err_extrap": _ext,
        "higbie_rel_err_raw": _rawh, "higbie_rel_err_extrap": _exth})

_kb_nf = BoundaryLayer(D=D_REF, u_s=0.0, gdot=GDOT, L=L_REF,
                       n_y=800, n_x=800).run(interface="noflux")[0]
defect("interface left flux-free instead of saturated",
       {"leveque_rel_err_raw": _kb_nf/ex_lev - 1, "leveque_rel_err_extrap": _kb_nf/ex_lev - 1},
       "a zero-flux interface transfers nothing: kbar collapses to 0")

def _lam(local=False, n=400):
    def expo(l):
        rel = 0.15
        out = []
        for D in (D_REF/(1+rel), D_REF*(1+rel)):
            b = BoundaryLayer(D=D, u_s=u_s_of(l), gdot=GDOT, L=L_REF, n_y=n, n_x=n).run()
            out.append(b[2] if local else b[0])
        return np.log(out[1]/out[0])/np.log((1+rel)**2)
    return (brentq(lambda l: expo(l) - 7/12, 0.05, 5.0, xtol=1e-4),
            brentq(lambda l: expo(l) - (0.5 + 0.1*(2/3-0.5)), 0.2, 100.0, xtol=1e-3))
_lh, _l90 = _lam(local=True)
defect("transition read from the LOCAL k at x = L, not the mean",
       {"lambda_half": _lh, "lambda_within_10pct_mobile": _l90})

# ---- rows kept and labelled because they do NOT move ----------------------
defect("sphere shell taken 10x deeper (1e-5 instead of 1e-4)",
       {"sh_stagnant_extrap": sh_sphere(1e-5, n_r=800)*(1-1e-5)},
       "STRUCTURAL by label, and it does drift -4.6e-4 relative: that residue is grid coarsening (n_r held while the domain deepens 10x), not the domain truncation the extrapolation removes")
defect("kbar from the boundary gradient instead of the mass balance",
       {"leveque_rel_err_raw":
        BoundaryLayer(D=D_REF, u_s=0.0, gdot=GDOT, L=L_REF, n_y=800, n_x=800).run()[1]
        / ex_lev - 1},
       "NEARLY STRUCTURAL: backward Euler conserves mass exactly, so the two routes "
       "agree by construction")

# ---- print ----------------------------------------------------------------
print(f"{'injected defect':<62}{'metric':<32}{'base':>12}{'broken':>12}{'move':>11}")
covered, n_rows, n_moved = set(), 0, 0
for nm, moves, note in DEFECTS:
    first = True
    for k, v in moves.items():
        base = M[k]
        rel = (v - base)/abs(base) if abs(base) > 1e-300 else np.nan
        moved = (not np.isfinite(v)) or abs(v - base) > 1e-9*max(1.0, abs(base))
        if moved and "STRUCTURAL" not in note.upper():
            covered.add(k); n_moved += 1
        n_rows += 1
        mv = ("  --" if not np.isfinite(rel)
              else (f"{rel:+.1%}" if abs(rel) < 20 else f"{rel:+.3g}x"))
        print(f"{(nm if first else ''):<62}{k:<32}{base:12.5g}"
              f"{(v if np.isfinite(v) else float('nan')):12.5g}{mv:>11}")
        first = False
    if note:
        print(f"{'':<62}   {note}")

print(f"\n{n_moved} of {n_rows} rows move their metric, over {len(DEFECTS)} injected defects.")
missing = [k for k in M if k not in covered]
assert not missing, f"metrics with no moving break row: {missing}"
print(f"COVERAGE: {len(covered)} of {len(M)} metrics have a defect row that moves them.")
if missing:
    print("NO ROW MOVES THESE, and each is labelled on the page:")
    for k in missing:
        print("   ", k)
else:
    print("Every metric in agreement.json is reachable by an injected defect.")
"""))

cells.append(md(r"""### Which printed constants the checks exercise

Computed from `p()`'s read log, not asserted — a claim of coverage that lives
only in prose is the "comment asserting a check that does not exist" defect
this repository documents."""))

cells.append(code(r"""unread = [k for k in PN.index if k not in P_READS]
print(f"{len(P_READS)} of {len(PN)} printed-number rows are read by p() somewhere on this page.")
print("NEVER READ BY ANY CHECK - transcribed and declared, exercised by nothing here:")
for k in unread:
    row = PN.loc[k]
    print(f"   {k:<18} {row['quantity']}  (printed {row['printed_value']})")
declared = {"sd_eq4", "fig6_span", "d_split"}
assert declared <= set(unread), (declared, unread)
print()
print("The three DECLARED unexercised constants are among them: the 66 % standard")
print("deviation and the 10^7 span are properties of Fig. 6, which is never read,")
print("and the 2.5 mm split enters no formula on this page. Note also that the")
print(f"large-bubble Range-of-Variables row starts at d = {RV.loc['this_work_large','d_lo']} cm "
      f"= {RV.loc['this_work_large','d_lo']*10:.1f} mm,")
print(f"BELOW the printed split of {PN.loc['d_split','printed_value']} cm - the boundary is "
      "corroborated by the diameter")
print("columns but tested by nothing here. Any remaining rows above are retyped as")
print("exact rationals in the sympy identity cells rather than read through p();")
print("a future CSV correction would not propagate to those cells.")
"""))

cells.append(code(r"""report_agreement("A3.6", M)

FLOOR = 1e-12
below = {k: v for k, v in M.items() if abs(v) < FLOOR}
print(f"\n{len(below)} of {len(M)} metrics are below check_agreement.py's ABS_FLOOR "
      f"of {FLOOR:g} and are\nTHEREFORE NOT COMPARED BY CI AT ALL:")
for k, v in below.items():
    print(f"   {k:<28} = {v:.3g}")
companions = {"chu_identity_resid": "eq6_prefactor_rel_dev, which carries Chu's 35.4",
              "eq10_exponent_resid": "eq10_implied_const",
              "ra_form_resid": "sh_stagnant_extrap and crossover_Sc"}
print("Each is an algebraic identity whose job is to certify a transcription, each HAS")
print("a break row above that makes it non-zero, and each has a companion metric that")
print("IS compared by CI:")
for k in below:
    print(f"   {k:<24} -> {companions[k]}")
print(f"\nRuntime so far: {time.time()-T_START:.1f} s")
"""))

cells.append(code(r"""from IPython.display import Markdown, display

# The headline sentences, with every number interpolated from what was computed
# above rather than typed.
display(Markdown(f'''
**Reproduced from the paper's own theory, nothing fitted.** Froessling + Allen
gives **{c28:.6f}** against a printed **{c28_printed}** ({c28/c28_printed-1:+.2%});
Friedlander + Stokes gives **{c34:.6f}** against a printed **{c34_printed}**
({c34/c34_printed-1:+.2%}). The two bracket the fitted **{C1}**.

**Reproduced from the paper's own algebra.** Eq. (5) + Chu &rarr; the printed
*P*/*v*: symbolic residual **{CHU_RESID:.1e}**. Eq. (6)'s prefactor
**{C6_derived:.6f}** against a printed **{C6_printed}**
({C6_derived/C6_printed-1:+.2%}). Eq. (9)'s prefactor **{C9_derived:.6f}**
against a printed **{C9_printed}** ({C9_derived/C9_printed-1:+.3%}), exponent
**{E9_exact:.4f}** against a printed **{E9_printed}**.

**Not reproduced, and reported rather than repaired.** At the paper's own
&Sigma; = {SIGMA}, eq. (6) as printed gives *j*<sub>d</sub> =
**{jd_printed_eq6:.4f}** Re<sup>-1/3</sup> ({EQ7_FROM_PRINTED:+.1%} on eq. (7)'s
printed **{C7}**) and the minimal &Sigma;<sup>-1</sup> repair of printed eq. (6)
gives **{jd_sigma_inv:.4f}** ({EQ7_FROM_SIGINV:+.1%}). Eq. (9) and
Chilton&ndash;Colburn cross at *N*<sub>Re</sub> = **{PIPE_CROSS_PRINTED:.0f}**,
not the printed **{RE_CLAIM:,.0f}**, and differ by up to
**{PIPE_MAXDEV_PRINTED:.1%}** over {RE_LO:,.0f}&ndash;{RE_HI:,.0f}, not
&plusmn;{MAXDEV_CLAIM:.0%} &mdash; though the deviation does take both signs inside the
range, so the paper's &plusmn; is right about the shape.

**Corrected from the paper's own correlation.** The small-bubble row's
*D*<sub>L</sub> exponent, printed as &minus;4, puts eq. (1)
**{R_HI_PRINTED_EXP:.1f}&times;** above that row's own printed *k*<sub>L</sub>
ceiling; &minus;6 lands **{abs(1-R_HI_SIX)*100:.1f}&nbsp;%** below it. The
control &mdash; the large-bubble row, printed &minus;6 &mdash; reproduces at
{hi2/obs2[1]:.3f} of its printed ceiling with no adjustment.

**pymrm.** Stagnant sphere **{SH_STAGNANT:.7f}** against the printed
**{SH_PRINTED}**. Boundary layer against the Leveque and Higbie closed forms at
**{LEV_EXT:+.1e}** and **{HIG_EXT:+.1e}** extrapolated ({LEV_RAW:+.1e} and
{HIG_RAW:+.1e} raw, observed orders {LEV_ORD:.2f} and {HIG_ORD:.2f}). Exponents
**{EXP_RIGID:.6f}** and **{EXP_MOBILE:.6f}** against the printed 2/3 and 1/2.

**New.** The two branches would cross at *N*<sub>Sc</sub> =
**{SC_STAR:.4f}**, far below either data set. The rigid-to-mobile transition is
halfway at &Lambda; = **{LAM_HALF:.4f}** and within 10&nbsp;% of the mobile
exponent by &Lambda; = **{LAM_90:.3f}**.
'''))
"""))

# --------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Three things, and two refusals.**

**(1) The transition is turned into a number.** The paper describes it three
times — Fig. 3's intermediate-size bubbles, the observation that deep pools drive
large-bubble coefficients down to rigid-sphere values, and the observation that
"surface active agents cause the high mass transfer coefficients for large
bubbles to be reduced to those obtained with rigid spheres" — and models it
nowhere. All three are the same statement: **the interface stops moving**. Once
that is written as $u(y) = u_s + \dot\gamma y$, one dimensionless group
$\Lambda = u_s/(\dot\gamma\,\delta_{\rm Lev})$ controls the whole thing, the
exponent $\mathrm d\ln k_L/\mathrm d\ln D_L$ interpolates monotonically from
$2/3$ to $1/2$, and the halfway point is computed above. The paper's two
correlations are the two ends of that one curve. This is the same instinct as
`A3.3`'s $\coth$ interpolation between film and renewal theory, applied to a
different pair of limits.

**(2) The paper's own arithmetic is re-run, and two steps do not close.** Eq. (6)
as printed cannot produce eq. (7) — the voidage enters as $\Sigma^{-1}$, not
$\Sigma^{-1/3}$, because $P/v$ carries $\Sigma^{-4}$ and eq. (4) takes its fourth
root — and even the corrected reading falls short of the printed 0.69. And
eq. (9) meets Chilton–Colburn near $N_{Re}\approx4500$, not 10 000, with a
worst-case deviation larger than the ±12 % claimed (the ± itself is right —
the deviation takes both signs inside the range). **What licenses
both findings is that the same machinery reproduces the pipe prefactor exactly**
($0.13\times0.0396^{1/4} = 0.05799$ against a printed 0.058) and the packed-bed
prefactor to 0.28 % — so eq. (4)'s constant, the flow-work argument and the
$j_d$ definition are all pinned down before anything is called wrong. That is the
`F2.3` order of operations: establish what is *not* free first.

**(3) A printed table exponent is corrected from the paper's own correlation.**
The small-bubble row's $D_L$ is printed as $\times10^{-4}$; eq. (1) evaluated on
that row's own property corners then stands 21× above the row's own printed
$k_L$ ceiling, while $\times10^{-6}$ lands within 2 % of it and $\times10^{-5}$
is ruled out at 4.6×. The test has a control — the large-bubble row, printed
$\times10^{-6}$, reproduces without adjustment — so the method can fail and
doesn't.

**The refusals.**

**This page cannot judge 0.31, 0.42 or 0.13.** They are regression lines through
Figs 1, 2 and 6; the only evidence about them is the scatter in those figures;
digitising them needs a maintainer review that is not available. Two things
*are* said about them without touching the figures: the paper's own two
theoretical routes bracket 0.31 (0.28 below, 0.34 above), and evaluating eqs. (1)
and (2) at the corners of the paper's own Range of Variables table produces $k_L$
windows that bracket the measured ones. **Neither is a validation.** The second
in particular is a *goodness of fit* — 0.42 was fitted to the very runs whose
$k_L$ range it is being compared against — and it is labelled as one everywhere
it appears. Its value is as a transcription check: a mis-read constant or
exponent would move it by orders, and it does not move.

**The pymrm boundary layer does not predict the constants.** It reproduces the
exponents and the interpolation between them, which is the mechanism claim. The
prefactors 0.31 and 0.42 absorb bubble shape, slip velocity, swarm hindrance and
the surface-area definition, none of which a flat interface of contact length $L$
contains. Claiming otherwise would be exactly the kind of quiet overclaim the
gallery's verifier brief exists to catch."""))

# --------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**Take the correlation.** Eqs. (1) and (2) are two lines of code and need no
solver. The rule for choosing between them, in the paper's own terms: bubbles
below ≈2½ mm, or any interface contaminated by surface-active material, or any
deep pool where a large bubble has had time to accumulate a skin → eq. (1);
clean, freshly formed bubbles above 2½ mm in a shallow pool → eq. (2). The
computed crossover $N_{Sc}^\star$ tells you when the choice stops mattering, and
the shaded bands on the first figure show that both of the paper's own data sets
sit far above it — so in a water-like liquid the two branches differ by nearly a
factor of four and you have to make the choice.

**Do not read eq. (2) as size-independent outside its window.** The three
crossings computed against Baird's theoretical equation all fall *outside* the
0.20–0.80 cm the large-bubble runs covered — below it in thin liquids, far above
it in thick ones — so the two expressions bracket the tested range rather than
coincide in it.

**Take the boundary-layer solver.** `BoundaryLayer` is the generic
"convection-diffusion in a thin layer over a moving or rigid wall" object. Change
`u_s`/`gdot` for a different interface mobility; give `u` any profile you like
(it is built once with `construct_coefficient_matrix` and never touched inside
the march); add a source term through `NumJac` if the solute reacts. It is `S4`
machinery marched in space rather than time, which is the same trick `S2` uses
for a plug-flow reactor.

**Take the exponent diagnostic.** `d_ln_k_d_ln_D` is the cheapest possible test
of which transport mechanism a numerical model is actually in, and it applies far
beyond bubbles: any time a model is claimed to be boundary-layer-controlled or
penetration-controlled, differentiate its output with respect to the diffusivity
and look at the number.

**Related pages, and what does *not* travel from them.** `A3.3`
(Danckwerts surface renewal) solves the same penetration limit this page's
mobile branch reduces to; the two are **not independent confirmations of each
other**, and this page's mobile-branch check is against the closed form
$2\sqrt{D u_s/\pi L}$, not against `A3.3`'s output. `A3.1` (Whitman two-film)
supplies the resistance-in-series bookkeeping that turns any $k_L$ from here into
an overall coefficient. `A3.4` (Wakao–Funazkri) is the packed-bed Sherwood
correlation that the eq. (6)→(7) branch of this paper is an early competitor to;
if you need packed-bed $k_L$, use that page, not eq. (7), for the reason
documented above. **No number, dataset or sentence from any of the three is
reused here.**

**What to be careful with.**

- Eq. (4)'s 66 % standard deviation is not a typo and is not a rounding of 6.6 %.
  A correlation with that scatter is an order-of-magnitude tool.
- Eq. (7) should not be used. The step that produces it does not reproduce, and
  Fallat's own measured correlation — which the paper prints — is available
  instead.
- The ±12 % claim about eq. (9) and Chilton–Colburn does not hold as stated; the
  recomputed crossing and worst case are printed above.
- The small-bubble row of the Range of Variables table has a wrong power of ten
  on $D_L$, and two further rows of the same column are unusable.
- Eq. (2) has no diameter in it while the theory it says it represents does, so
  "independent of bubble size" is a statement about a window whose edges are
  computed above and move by a factor of ~20 in diameter across the paper's own
  viscosity range."""))

cells.append(code(r'''print(f"Total runtime: {time.time()-T_START:.1f} s")
'''))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print("wrote index.ipynb with", len(cells), "cells")
