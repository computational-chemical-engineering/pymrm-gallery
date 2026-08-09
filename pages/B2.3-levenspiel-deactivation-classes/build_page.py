#!/usr/bin/env python3
"""Generate index.ipynb for page B2.3 (Levenspiel deactivation classes).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Levenspiel's four deactivation classes: which reactor can tell them apart"
description: "Levenspiel (1972) asks whether experiments can be devised that return the orders of an nth-order deactivating rate equation, for the four classes parallel, series, side-by-side and independent. The paper prints no data at all, so this page proves its printed structure (26 symbolic identities) and then measures the question it poses: all four classes solved in all four of its batch-solids contactings, and every class fitted to every other class's response - 64 fits. The answer is harder than the taxonomy suggests. FIVE exact degeneracy families reach ALL FOUR devices, the plug-flow bed included: in a well-mixed fluid series deactivation of order (n', d) IS parallel of order (n n', d+n') for any n; in the held-C_A device parallel, series and independent are literally the same equation; and at the base orders n = n' = d = 1, in a batch of fluid or a bed, side-by-side deactivation with poison group sigma IS parallel deactivation with n' -> lambda sigma / kappa, for EVERY sigma - the poison profile is a power of the reactant profile there, and it is NOT in the two well-mixed devices, where the same identity fails by 1e-1. What the bed and the batch keep - and the two well-mixed devices do not - is a distinction between parallel, series and independent; what the bed alone has is a spatial coordinate. At a 0.2-percentage-point conversion resolution NO device separates more than one of the twelve class pairs. What does discriminate is not the outlet history but the activity PROFILE: after the same run the bed's inlet sits at exactly a = 1 under series decay (no R has been made there) and at exactly a = e^(-lambda) = 0.0498 under parallel, a contrast of opposite sign that one spent-catalyst assay would see. Eq. (24) as printed is offset from the paper's own eq. (21) by exactly ln(k_d/k)."
categories: [sec:B, struct:S1, struct:S5, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-08
---

# Levenspiel's four deactivation classes: which reactor can tell them apart

**Catalog ID:** `B2.3` · **Structures:** `S1`, `S5` · **Tier:** T0

A deactivating catalyst has three things going on at once — temperature,
composition, and the activity of the pellet itself — and only two of them are
under the experimenter's control. Levenspiel's 1972 paper is about the third.
It writes the rate of a deactivating pellet as a product,

$$-r'_{\rm A} = f_1(\text{fluid conditions})\cdot f_2(\text{present activity}), \tag{1}$$

instantiates it as $n$th-order kinetics with $n'$-order, $d$-order decay in
**four broad classes** — parallel, series, side-by-side, independent (its
eqs. 4–7) — and then asks the question this page is about:

> *"We ask in particular whether experiments can be devised so as to give
> simply the orders and rate constants of these equations."* (p. 272)

That is a question about **experimental design**, and the paper answers it in
words: one particular contacting, it says, "alone allows decoupling of
concentration and activity effects". This page answers it with numbers.
It solves all four classes in all four of the paper's own batch-solids
contactings, and then does the only thing that can settle a claim of
discriminability: **fits every class to every other class's response** and
reports the conversion resolution an experiment would need to tell them apart.

The answer is harder than the taxonomy suggests, in two stages.

**First, the four devices are not merely poor at separating the classes — every
one of them is *exactly* degenerate somewhere**, by algebraic identities the
page proves symbolically and then confirms by subtraction. In any well-mixed
fluid at constant flow, series deactivation of order $(n', d)$ **is** parallel
deactivation of order $(n n',\, d+n')$, for any reaction order $n$: same
observable, every time. In the held-$C_{\rm A}$ device Levenspiel recommends,
parallel, series and independent are literally the same equation — decoupling
concentration from activity is exactly what destroys the class distinction
inside one run. And in a **batch of fluid or a plug-flow bed** — the two
devices that do keep parallel, series and independent apart — side-by-side
deactivation with poison group $\sigma$ **is** parallel deactivation with
$n' \to \lambda\sigma/\kappa$, for *every* $\sigma$, at the base orders
$n = n' = d = 1$. The bed is not the
exception to the degeneracies; it is one of the two devices whose degeneracy is
narrower.

**Second, the two devices that keep those three classes apart are still not
enough.** At a conversion resolution of 0.2 percentage points, the plug-flow
bed separates one of the twelve class pairs, and no other device separates any.
Parallel and series in a bed differ by 8.1e-4 in RMS conversion — real, but a
tenth of a percentage point.

What *does* discriminate is not the outlet history at all, and it is the one
thing only the bed has: an **activity profile**. After the same run, a bed
decaying in parallel ends with its inlet at $a = e^{-\lambda}$ and its outlet
near 0.19, while the same bed decaying in series ends with its inlet at
**exactly 1** — no R has been made at the feed face, so nothing has decayed
there — and its outlet near 0.18. Both inlet values are exact, not computed;
the contrast reverses sign, and one assay of the spent catalyst reads it off.
That is `B2.2`'s finding — the profile names the mechanism — arrived at from
the opposite direction."""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

**The source, precisely.** Octave Levenspiel, *"Experimental Search for a
Simple Rate Equation to Describe Deactivating Porous Catalyst Particles"*,
**Journal of Catalysis 25**(2) 265–272 (1972),
[doi:10.1016/0021-9517(72)90227-8](https://doi.org/10.1016/0021-9517(72)90227-8),
Department of Chemical Engineering, Oregon State University, Corvallis,
Oregon 97331; received September 28, 1971. Identified from the document's own
first page — running head "JOURNAL OF CATALYSIS 25, 265–272 (1972)", display
title, the single-author by-line "OCTAVE LEVENSPIEL", the affiliation block,
the received date and the abstract — on a native-resolution render. The scan is
CCITT-G4 bilevel at **300 ppi native**, 8 pages, complete.

**The text layer of this file is not usable and was not used.** Its OCR
mangles the abstract into "the experimental detrrmination of the ordrrs of
reaction and deactivation", "deactivation of catalyst particlrs", "parallrl",
"the various rractor types which may he usctf". Every symbol, subscript,
prime and exponent on this page was read from a cropped 300 ppi render at
glyph scale — which matters more here than usual, because the whole paper
turns on telling $n$ from $n'$, $k$ from $k'$ from $k_d$ from $k'_d$, and
$C_{\rm A}$ from $C_{\rm R}$ from $C_{\rm P}$.

**This is Levenspiel the author, not Levenspiel the book.** *Chemical
Reaction Engineering* is a different work (this paper cites its 2nd edition,
Chap. 15, as its own ref. 3); it remains outstanding as a textbook-canonical
source and is what `A2.4` is catalogued to. Nothing on this page closes that.

**Where it sits.** Third of a ladder. `B2.1` is
[Voorhies (1945)](../B2.1-voorhies-coking-law/) — coke as a power of
time-on-stream, $C_c = A\theta^n$, activity correlated with the clock.
`B2.2` is [Froment & Bischoff (1961)](../B2.2-froment-bischoff-coking/) —
activity correlated with the carbon actually deposited, which becomes a
*profile* along a bed whose direction names the fouling mechanism. This paper
takes the third step: activity as a state variable $a$ with its own $n'$-,
$d$-order rate law, and the question of what experiment could ever measure
those orders. The connection to `B2.2` turns out to be sharper than a shared
subject: the single device here that supports a spatial gradient is also the
one whose exact degeneracy is narrowest, and the measurement that actually
separates the classes in it is the shape of the deposit — which is `B2.2`'s
finding arrived at from the opposite direction.

**What the paper contains, and what it does not.** It contains 32 numbered
equations, six figures and no measurements. There is **no table anywhere in
the eight pages**; Figs. 2–6 are schematic test plots with unlabelled axes,
no tick values, and open circles drawn to illustrate a straight line
("A straight line on this plot indicates that the guessed order of
deactivation is correct", Fig. 6 caption). The paper's title word is
*Search*, and its Discussion states the subject plainly: whether experiments
*can be devised*. **No comparison against any measured order is attempted or
promised on this page, because the paper reports none.**"""))

# ----------------------------------------------------------- published model
cells.append(md(r"""## The published model

All equation numbers are the paper's. Every symbol below was read from a
cropped native-resolution render.

### The four classes (eqs. 4–7, p. 266)

Each class pairs the same $n$th-order reaction rate with a different driving
concentration for the decay. $a$ is the activity, defined in eq. (3) as the
ratio of the pellet's present rate to its rate when fresh, so $a(0) = 1$.

| eq. | class | scheme | rate of reaction | rate of decay |
|---|---|---|---|---|
| (4) | *parallel* | A→R, A→P↓ | $-r'_{\rm A} = k\,C_{\rm A}^{\,n} a$ | $-\dfrac{da}{dt} = k_d\,C_{\rm A}^{\,n'} a^{\,d}$ |
| (5) | *series* | A→R→P↓ | $-r'_{\rm A} = k\,C_{\rm A}^{\,n} a$ | $-\dfrac{da}{dt} = k_d\,C_{\rm R}^{\,n'} a^{\,d}$ |
| (6) | *side-by-side* | A→R, P→P↓ | $-r'_{\rm A} = k\,C_{\rm A}^{\,n} a$ | $-\dfrac{da}{dt} = k_d\,C_{\rm P}^{\,n'} a^{\,d}$ |
| (7) | *independent* | — | $-r'_{\rm A} = k\,C_{\rm A}^{\,n} a$ | $-\dfrac{da}{dt} = k_d\,a^{\,d}$ |

A fifth scheme (eq. 8, A→R with both A and R poisoning) is dismissed by the
paper in one line — $-da/dt = k_d(C_{\rm A}+C_{\rm R})^{n'}a^d$ and
"since $C_{\rm A} + C_{\rm R}$ remains constant for a specific feed, this type
of deactivation reduces to the simple-to-treat independent deactivation of
Eq. 7". Verified below.

The Nomenclature (p. 265) defines the shorthands that carry the rest of the
paper: $k' = k\,C_{\rm A}^{\,n}$, $k'_d = k_d\,C_{\rm A}^{\,n'}$,
$k'' = kW/V$, and the **weight-time**
$\tau' = W C_{\rm A0}/F_{\rm A0}$ (g cat·sec/liter of fluid), the capacity
factor that plays the role of space time for a fixed weight of catalyst.

### The four contactings (Fig. 1, p. 268)

Fig. 1's top row is the batch-solids family this paper is "primarily concerned
with": **(a)** batch of both solid and fluid, **(b)** plug flow of fluid,
**(c)** mixed flow of fluid, **(d)** recycle flow of fluid. (The bottom row —
raining solids and a fluidised bed — is for deactivation "in the order of
seconds or less" and is out of scope here, as is (d), which the paper itself
says "offers no particular advantage".) On top of the contacting pattern the
paper adds a second axis: the flow rate may be held **constant**, or lowered
continually so as to hold a concentration constant.

For the illustrative first-order, concentration-independent case (eq. 9,
$-r'_{\rm A} = kC_{\rm A}a$ with $-da/dt = k_d a$) the paper works each
contacting through to a straight-line test:

| contacting | performance relation | test plot |
|---|---|---|
| (a) batch fluid | $\ln\ln\dfrac{C_{\rm A}}{C_{\rm A\infty}} = \ln\ln\dfrac{C_{\rm A0}}{C_{\rm A\infty}} - k_d t$ (15) | Fig. 2 |
| (c) mixed, constant flow | $\ln\!\left(\dfrac{C_{\rm A0}}{C_{\rm A}}-1\right) = \ln(k\tau') - k_d t$ (19) | Fig. 3 |
| (b) plug, constant flow | $\ln\ln\dfrac{C_{\rm A0}}{C_{\rm A}} = \ln(k\tau') - k_d t$ (22) | Fig. 4 |
| (c) mixed, $C_{\rm A}$ held | $\ln\tau' = k_d t + \ln\dfrac{C_{\rm A0}-C_{\rm A}}{k\,C_{\rm A}}$ (23) | Fig. 5 |
| (b) plug, $C_{\rm A,out}$ held | $\ln\tau' = k_d t + \ln\!\left(\dfrac{1}{k_d}\ln\dfrac{C_{\rm A0}}{C_{\rm A}}\right)$ (24) | Fig. 5, shifted |

**Eq. (24) is transcribed exactly as printed, and it is wrong as printed.**
The paper derives it from its own eq. (21) one line above, and eq. (21) gives
$1/k$ where eq. (24) prints $1/k_d$. The discrepancy is settled symbolically
and dimensionally below; it is **reported, not repaired**, and the CSV carries
the printed form.

### The general-order analysis (eqs. 25–32, p. 271)

The paper's recommendation. Hold $C_{\rm A}$ constant by lowering the flow
rate; then $k' = kC_{\rm A}^{\,n}$ and $k'_d = k_dC_{\rm A}^{\,n'}$ are
constants, eqs. (4) collapse to

$$-r'_{\rm A} = k'a \tag{25}\qquad -\frac{da}{dt} = k'_d\,a^{\,d} \tag{26}$$

the mixed-flow performance relation becomes $\tau' = (C_{\rm A0}-C_{\rm A})/(k'a)$
(27), and integrating (26) for each $d$ gives the four printed straight lines
(28)–(31) and their general form

$$(\tau')^{\,d-1} = C_1 + C_2\,t, \qquad
C_1 = \left(\frac{C_{\rm A0}-C_{\rm A}}{k'}\right)^{d-1}, \quad
C_2 = C_1\,(d-1)\,k'_d \tag{32}$$

with $C_1$, $C_2$ printed on Fig. 6 as its intercept and slope. All of
(28)–(32) are re-derived symbolically below."""))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

**Everything is dimensionless, and every group is a ratio the paper itself
defines.** Writing $c = C_{\rm A}/C_{\rm A0}$, $r = 1-c = C_{\rm R}/C_{\rm A0}$
(constant density, A→R one-for-one), $p = C_{\rm P}/C_{\rm P0}$,
$s = w/W$ the fraction of the catalyst weight passed, and $\theta = t/t_{\rm obs}$
the time in units of the observation window:

| group | definition | meaning | base value |
|---|---|---|---|
| $\kappa$ | $k\,C_{\rm A0}^{\,n-1}\tau'$ | Damköhler number on fresh catalyst | 1.5 |
| $\lambda$ | $k_d\,C_{\rm ref}^{\,n'}\,t_{\rm obs}$ | decay rate × observation window | 3.0 |
| $n$ | — | reaction order (the paper's own worked case is $n=1$) | 1 |
| $n'$ | — | concentration order of the decay | 1 |
| $d$ | — | order of the decay in activity | 1 |
| $\sigma$ | see below | poison consumed per unit activity lost | 0.8 |
| $c_{\rm set}$ | — | the level $C_{\rm A}$ is held at in the constant-$C_{\rm A}$ device | 0.4 |

The base point puts fresh-catalyst conversion at 60 % (mixed) / 78 % (plug)
and runs until the activity has fallen by roughly an order of magnitude —
i.e. a run an experimenter would actually design. **Every conclusion below is
stated with the $(n, n', d)$ it was computed at**, and the ones that depend on
those orders are mapped, not asserted.

**$\sigma$ is not printed, and is declared here.** Eq. (6) reads
"P → P↓": the poison is *consumed* by depositing, so a poison mass balance is
unavoidable, but the paper never gives the stoichiometry. $\sigma$ is the
dimensionless group multiplying it — moles of P removed from the fluid per
unit of activity destroyed, scaled by $\tau'/C_{\rm P0}$. It is a
reconstruction, not a transcription; it is swept, and the $\sigma \to 0$ limit
(a poison in such excess that its concentration never moves) is reported
separately because it is exactly degenerate with independent deactivation.

**Assumptions carried from the source.**

- **Pseudo-steady fluid**, stated by the paper on p. 269: *"this and the
  following derivations for a batch of solids are based on the
  pseudo-steady-state assumption ... Since a batch of solids can only be used
  in experimentation if deactivation is not too rapid this assumption is
  reasonable."* Everything except the last section assumes it. The last
  section drops it and measures what it costs.
- **Isothermal, constant density, single reaction A→R** with no volume change
  ($\epsilon_{\rm A} = 0$), so $C_{\rm R} = C_{\rm A0} - C_{\rm A}$ everywhere.
  This is the paper's own setting: it prints $\epsilon_{\rm A}$ in the
  Nomenclature but never uses it.
- $a(0) = 1$ for the whole batch of solids ("The activity starts at unity",
  p. 266).

**Measurement resolution.** A discriminability claim is meaningless without
one, so every separability number below is reported as a **raw RMS deviation
in fractional conversion** — the resolution an experiment must have to refute
the wrong class. For labelling only, `RES_X` = 2e-3 (0.2 conversion
percentage points, a good laboratory chromatograph) and `RES_TAU` = 0.01
(1 % in flow rate) are used; the raw numbers are printed beside every verdict
so a reader with a different instrument can rescale.

**Direction of the separability claim, stated once.** $D(T\!\to\!C)$ is the
*minimum* misfit of class $C$ to a response generated by class $T$. An
optimiser returns an **upper bound** on that minimum. So a large $D$ means
"no fit this good was found" — separation could be worse than reported, never
better. The claims that run the other way — $D = 0$, *not* separable — are
therefore the safe ones, and every one of them is backed here by an exact
algebraic identity rather than by an optimiser."""))

# ------------------------------------------------------------ colab cell ----
cells.append(code(r'''try:
    import pymrm  # noqa: F401
except ImportError:
    %pip install -q pymrm pyyaml
'''))

# ------------------------------------------------------------- imports ------
cells.append(code(r'''"""Bootstrap, data, and the two ledgers (metrics M, break rows BREAKS)."""
import sys, pathlib, time
if "google.colab" in sys.modules:
    import urllib.request
    base = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
            "pymrm-gallery/main/shared/gallery_utils.py")
    urllib.request.urlretrieve(base, "gallery_utils.py")
else:
    for p in (pathlib.Path.cwd(), *pathlib.Path.cwd().parents):
        if (p / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(p / "shared")); break

import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "B2.3-levenspiel-deactivation-classes"
EQS    = load_data("levenspiel-1972-printed-equations.csv", page=PAGE)
EQS["eq"] = EQS["eq"].astype(str)
EQS = EQS.set_index("eq")
CLAIMS = load_data("levenspiel-1972-printed-claims.csv", page=PAGE).set_index("claim_id")
print(cite_data(load_meta("levenspiel-1972-printed-equations.csv", page=PAGE)))

def md_table(df, index=False):
    """DataFrame -> pipe table (no tabulate dependency; requirements.txt has none)."""
    d = df.reset_index() if index else df
    cols = [str(c) for c in d.columns]
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in d.iterrows():
        lines.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(lines)

M = {}          # agreement metrics, assembled across the page
BREAKS = []     # defect-injection rows: (metric, base, injection, result, note)
plt.rcParams.update({"figure.dpi": 110, "font.size": 9.5})

CLASSES = ("parallel", "series", "side-by-side", "independent")
SHORT = {"parallel": "par", "series": "ser", "side-by-side": "s-b-s",
         "independent": "ind"}
# Okabe-Ito, colour-vision-deficiency safe; fixed order used on every figure
COL = {"parallel": "#0072B2", "series": "#E69F00",
       "side-by-side": "#009E73", "independent": "#CC79A7"}
DEV = {"batch": "Fig. 1a  batch fluid",
       "plug": "Fig. 1b  plug flow, constant rate",
       "mixed": "Fig. 1c  mixed flow, constant rate",
       "constC": "Fig. 1c  mixed flow, $C_A$ held"}

RES_X   = 2.0e-3     # illustrative conversion resolution (0.2 percentage points)
RES_TAU = 1.0e-2     # illustrative resolution in ln(tau'), i.e. 1 % in flow rate
BASE = dict(kappa=1.5, lam=3.0, d=1.0, nprime=1.0, n=1.0, sig=0.8, c_set=0.4)
TH = np.linspace(0.0, 1.0, 41)         # observation times, theta = t/t_obs
print(f"{len(EQS)} printed equations and {len(CLAIMS)} printed claims loaded")
'''))

# -------------------------------------------------------------- the data ----
cells.append(md(r"""## The data

**There are none.** This paper publishes no measurement of any kind, and that
is a claim about the source, so here is the search that establishes it.

All eight pages (265–272) were rendered at the scan's native 300 ppi and read
as images. **There is no table in the document** — no `TABLE` head, no ruled
block, no column of numbers anywhere in the eight pages. (The two-column
Nomenclature on p. 265 is a list of symbol definitions, not data.) The numeral
density is 70–112 glyphs per page and is almost entirely equation numbers: the
only numerals in the body text that are not equation numbers, figure or
section cross-references, page numbers, dates or bibliographic detail are the
two timescales on p. 267 ("0.1–1 sec", "minutes") and the orders on p. 272,

> $n = 1;\quad n' = 0,1;\quad d = 1,3$

— and those are attributed to the paper's refs. 3 and 4 (Levenspiel's own
monograph and Kunii & Levenspiel's *Fluidization Engineering*) as the orders
*those* books show how to test, not as anything measured here. **Neither was
consulted**, and no number from either appears on this page.

Figs. 2–6 are **schematic**, not data. Their abscissa is an unlabelled $t$;
their ordinate is an algebraic expression, not a scale; there is not a single
tick value on any of the ten axes; and the open circles are drawn to
illustrate what a correct guess would look like. The Fig. 6 caption says so:
*"A straight line on this plot indicates that the guessed order of deactivation
is correct."* Nothing was digitised, and nothing could have been.

The Discussion states the paper's own subject in the same terms
(p. 272): *"We ask in particular whether experiments can be devised so as to
give simply the orders and rate constants of these equations."* The abstract's
damaged phrase "the experimental detrrmination of the ordrrs" means *how to
determine* them; the title word is **Search**.

**Consequence for this page.** Nothing here is fitted to measurement and
nothing here is validated against measurement. Every number is one of three
things, labelled wherever it appears: (i) exact reproduction of printed
algebra; (ii) an exact algebraic identity between the paper's own rate forms;
(iii) a computed property of those rate forms in a reactor. The two CSVs under
`data/` are transcriptions — the printed equations this page uses, and the
printed claims it tests, quoted verbatim; both counts are printed by the
bootstrap cell above rather than typed here.

**Data tier 6** (model-defining constants transcribed from the source, no
experimental dataset)."""))

cells.append(code(r'''"""The transcription, and the one thing it can be checked against: the paper's
own cross-references.  Each of eqs. (12)-(32) is derived by the paper from a
numbered predecessor, so the transcription is re-derived, line by line, in the
symbolic section below rather than merely retyped."""
display(Markdown("**The four classes as printed (eqs. 4-7):**\n\n" +
                 md_table(EQS.loc[["4", "5", "6", "7"],
                                  ["ascii_form", "describes"]], index=True)))
n_defect = int((EQS["status"] != "verified").sum())
print(f"\nrows flagged as printed defects in the transcription: {n_defect} "
      f"-> {list(EQS.index[EQS['status'] != 'verified'])}")
print(f"claims tested: {len(CLAIMS)};  of these, "
      f"{int(CLAIMS['tested_as'].str.startswith('not tested').sum())} explicitly out of scope")
'''))

# ------------------------------------------------- pymrm implementation -----
cells.append(md(r"""## PyMRM implementation

Three of the four contactings are lumped systems — a batch of fluid, or a
well-mixed fluid over a batch of solids — and are ODEs in $\theta$ (`S1`).
The fourth, plug flow of fluid over a batch of solids, is a one-dimensional
convection problem in the catalyst-weight coordinate $s = w/W$ with a
distributed activity $a(s,\theta)$ (`S5`), and that is the pymrm model:

$$\varepsilon\,\frac{\partial c}{\partial \theta} + \frac{\partial c}{\partial s}
   = -\kappa\,c^{\,n} a, \qquad
   \frac{\partial a}{\partial \theta} = -\lambda\,\phi^{\,n'} a^{\,d}$$

with $\phi = c$, $1-c$, $p$ or $1$ for the four classes, $c(0,\theta) = 1$ at
the feed, and

$$\varepsilon \;=\; \frac{\text{fluid residence time}}{\text{observation window}}$$

the small parameter the paper's pseudo-steady assumption sets to zero.

- **Transport assembled once** in `__init__`: `construct_convflux_upwind` with
  $v = 1$ on the weight coordinate, `construct_div` with `nu=0` (Cartesian —
  $s$ is a *weight* fraction, not a radius, so there is no area factor), and a
  van Leer deferred correction through `interp_cntr_to_stagg_tvd` for
  second-order accuracy on a convection-only equation.
- **Boundary conditions on the outward normal.** Inlet `{a:0, b:1, d:1}` is
  $c = 1$ (the feed); outlet `{a:1, b:0, d:0}` is $\partial c/\partial n = 0$,
  the standard outflow condition. The outlet reading is taken with
  `compute_boundary_values`, which — with the outflow bc left as pymrm builds
  it — is the **second-order** reading. (`B2.2` measured both configurations
  of this on the same operator stack and found the pairing matters: correcting
  the outflow face flux makes the *extrapolated* face second order and leaves
  `compute_boundary_values` first order, and vice versa. This page uses the
  uncorrected pairing, and the observed order is measured below, not assumed.)
- **`NumJac` is not used**: at frozen $a$ the fluid equation is linear for
  $n = 1$ and needs one sparse solve; for $n \ne 1$ the Newton Jacobian is the
  operator plus a diagonal, written directly. The activity is a *pointwise*
  ODE at each cell centre with no neighbour coupling, so it is marched with
  Heun (RK2) rather than through a Jacobian at all.
- **The pseudo-steady route ($\varepsilon = 0$)** and the **full transient
  ($\varepsilon > 0$)** are the same class, the same operators and the same
  boundary conditions — only the accumulation term differs. That is what makes
  the cost of the paper's assumption measurable rather than arguable.

**A second, independent route.** For the pseudo-steady bed at frozen $a(s)$
the fluid equation is a first-order cascade in $s$, so
$c(s) = [\,1-(1-n)\kappa\!\int_0^s a\,]^{1/(1-n)}$ exactly (and $e^{-\kappa\int a}$
for $n = 1$). That quadrature shares **no** assembly, no operator and no
solver with the finite-volume bed, and it is used as the reference the bed is
measured against — and, because it is 60× cheaper, as the engine of the
separability sweep once the bed has certified it."""))

cells.append(code(r'''"""The four reactor models.  All dimensionless (see Parameters)."""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_simpson
from scipy.optimize import brentq, least_squares, minimize_scalar

def _pow(x, e):
    """x**e, safe at x = 0 for fractional e."""
    return np.exp(e * np.log(np.clip(x, 1e-300, None)))

def _cascade(I, g, order):
    """Exact solution of dy/ds = -g y^order, y(0) = 1, given I(s) = int_0^s of
    the (known) coefficient profile."""
    if g == 0.0:
        return np.ones_like(I)
    if order == 1.0:
        return np.exp(-g * I)
    return _pow(np.clip(1.0 - (1.0 - order) * g * I, 0.0, None), 1.0 / (1.0 - order))

# --------------- the activity carried in the paper's own linearising variable
# Eq. (26) is -da/dt = k'_d a^d.  With  w = (a^(1-d) - 1)/(1-d)  (= ln a at
# d = 1) it becomes  dw/dt = -k'_d  for EVERY d - which is the substitution
# behind the paper's own eqs. (28)-(32).  Marching w instead of a removes a^d
# from every right-hand side, so d < 1 (finite-time extinction; the paper's
# printed d = 0 case) stops being a stiff kink: 6x fewer function evaluations
# at d = 0 and 25x at d = 0.5, with identical answers where both work.
TINY_D = 1e-8

def a_of_w(w, d):
    """Invert w -> a; a = 0 once the activity is used up (only possible d < 1)."""
    w = np.asarray(w, float)
    if abs(1.0 - d) < TINY_D:
        return np.exp(w)
    x = (1.0 - d) * w
    ok = x > -1.0
    return np.where(ok, np.exp(np.log1p(np.where(ok, x, 0.0)) / (1.0 - d)), 0.0)

def _extinction_event(d, idx=0):
    """Terminal event for the lumped reactors: with d < 1 the activity reaches
    ZERO at a finite time (k'_d t = 1/(1-d) in w), and integrating past it is
    both meaningless and where every stiff step of this page was being spent."""
    def ev(_, y):
        return 1.0 + (1.0 - d) * y[idx] if d < 1.0 - TINY_D else 1.0
    ev.terminal, ev.direction = True, -1
    return ev

def w_of_a(a, d):
    a = np.asarray(a, float)
    if abs(1.0 - d) < TINY_D:
        return np.log(a)
    return (_pow(a, 1.0 - d) - 1.0) / (1.0 - d)

# ---------------------------------------------------------- Fig. 1a ---------
def run_batch(cls, kappa, lam, d, nprime, n=1.0, sig=0.0, theta=TH, rtol=1e-9):
    """Batch of solids AND batch of fluid.  Observable: conversion X_A(theta)."""
    def rhs(_, y):
        c, a, p = max(y[0], 1e-300), float(a_of_w(y[1], d)), max(y[2], 1e-300)
        f = {"parallel": c, "series": 1.0 - c, "side-by-side": p}.get(cls, 1.0)
        dw = -lam * _pow(max(f, 0.0), nprime)
        return [-kappa * _pow(c, n) * a, dw, sig * dw * _pow(a, d)]
    s = solve_ivp(rhs, (0.0, theta[-1]), [1.0, 0.0, 1.0], t_eval=theta,
                  rtol=rtol, atol=1e-13, method="LSODA",
                  events=_extinction_event(d, idx=1))
    c, a = np.empty(theta.size), np.zeros(theta.size)   # dead catalyst: a = 0,
    k = s.y[0].size                                     # and c frozen
    c[:k], a[:k] = s.y[0], a_of_w(s.y[1], d)
    c[k:] = s.y[0][-1] if k else 1.0
    return 1.0 - np.clip(c, 0.0, 1.0), a

# ---------------------------------------------------------- Fig. 1c ---------
def _mixed_c(kappa, a, n):
    """Solve 1 - c = kappa c^n a (the mixed-flow performance relation, eq. 16)."""
    if a <= 0.0:
        return 1.0
    if n == 1.0:
        return 1.0 / (1.0 + kappa * a)
    return brentq(lambda c: 1.0 - c - kappa * _pow(c, n) * a, 1e-14, 1.0,
                  xtol=1e-15, rtol=8.9e-16)

def _mixed_p(g, nprime):
    """Solve 1 - p = g p^n' (the poison balance over a well-mixed reactor)."""
    if g <= 0.0:
        return 1.0
    if nprime == 1.0:
        return 1.0 / (1.0 + g)
    if nprime == 0.0:
        return max(1.0 - g, 0.0)
    return brentq(lambda p: 1.0 - p - g * _pow(p, nprime), 1e-14, 1.0,
                  xtol=1e-15, rtol=8.9e-16)

def run_mixed(cls, kappa, lam, d, nprime, n=1.0, sig=0.0, theta=TH, rtol=1e-9):
    """Batch of solids, mixed flow of fluid at CONSTANT flow rate."""
    def rhs(_, y):
        a = float(a_of_w(y[0], d))
        c = _mixed_c(kappa, a, n)
        f = (_mixed_p(sig * lam * _pow(a, d), nprime) if cls == "side-by-side"
             else {"parallel": c, "series": 1.0 - c}.get(cls, 1.0))
        return [-lam * _pow(max(f, 0.0), nprime)]
    s = solve_ivp(rhs, (0.0, theta[-1]), [0.0], t_eval=theta, rtol=rtol,
                  atol=1e-13, method="LSODA", events=_extinction_event(d))
    a = np.zeros(theta.size)                            # dead catalyst: a = 0
    a[:s.y[0].size] = a_of_w(s.y[0], d)
    return 1.0 - np.array([_mixed_c(kappa, ai, n) for ai in a]), a

# ------------------------------------------- Fig. 1c, flow rate changing ----
def run_const_c(cls, kappa, lam, d, nprime, c_set, n=1.0, sig=0.0, theta=TH,
                rtol=1e-11):
    """Mixed flow with the flow rate lowered continually so that C_A stays at
    c_set - the contacting the paper recommends.  Observable: tau'(theta),
    relative to its value on fresh catalyst."""
    r_set = 1.0 - c_set
    def tau_rel(a):
        return r_set / (kappa * _pow(c_set, n) * max(a, 1e-300))
    def rhs(_, y):
        a = float(a_of_w(y[0], d))
        f = (_mixed_p(sig * lam * _pow(a, d) * tau_rel(a), nprime)
             if cls == "side-by-side"
             else {"parallel": c_set, "series": r_set}.get(cls, 1.0))
        return [-lam * _pow(max(f, 0.0), nprime)]
    s = solve_ivp(rhs, (0.0, theta[-1]), [0.0], t_eval=theta, rtol=rtol,
                  atol=1e-14, method="LSODA", events=_extinction_event(d))
    a = np.zeros(theta.size)          # past extinction the flow rate needed to
    a[:s.y[0].size] = a_of_w(s.y[0], d)   # hold C_A is zero, i.e. tau' -> inf
    return np.array([tau_rel(ai) for ai in a]), a

# ---------------------------------------------------------- Fig. 1b ---------
class PlugExact:
    """Plug flow of fluid over a batch of solids, pseudo-steady, by quadrature.

    Shares nothing with the pymrm bed below: no operator, no grid stencil, no
    linear solve.  Cumulative Simpson is LINEAR in its integrand on a fixed
    grid, so its matrix is built once and every quadrature is one matvec."""
    def __init__(self, ns=21):
        assert ns % 2 == 1, "Simpson needs an odd node count"
        self.s = np.linspace(0.0, 1.0, ns); self.ns = ns
        self.Q = np.column_stack([cumulative_simpson(np.eye(ns)[i], x=self.s,
                                                     initial=0.0)
                                  for i in range(ns)])

    def profiles(self, a, kappa, lam, n, nprime, d, sig, cls):
        c = _cascade(self.Q @ a, kappa, n)
        p = (_cascade(self.Q @ _pow(a, d), sig * lam, nprime)
             if (cls == "side-by-side" and sig > 0.0) else np.ones_like(a))
        return c, p

    def run(self, cls, kappa, lam, d, nprime, n=1.0, sig=0.0, theta=TH, rtol=1e-8):
        def rhs(_, w):
            a = a_of_w(w, d)
            c, p = self.profiles(a, kappa, lam, n, nprime, d, sig, cls)
            f = {"parallel": c, "series": 1.0 - c,
                 "side-by-side": p}.get(cls, np.ones_like(a))
            return -lam * _pow(np.clip(f, 0.0, None), nprime)
        sol = solve_ivp(rhs, (0.0, theta[-1]), np.zeros(self.ns), t_eval=theta,
                        rtol=rtol, atol=1e-12, method="LSODA")
        A = a_of_w(sol.y, d)
        X = np.array([1.0 - self.profiles(A[:, j], kappa, lam, n, nprime, d,
                                          sig, cls)[0][-1]
                      for j in range(theta.size)])
        return X, A

PLUG = PlugExact(21)

def respond(config, cls, kappa, lam, d, nprime, sig, n=1.0, c_set=0.4, theta=TH):
    """The observable of one device: conversion, or ln tau' for the held-C_A one."""
    if config == "batch":
        return run_batch(cls, kappa, lam, d, nprime, n=n, sig=sig, theta=theta)[0]
    if config == "mixed":
        return run_mixed(cls, kappa, lam, d, nprime, n=n, sig=sig, theta=theta)[0]
    if config == "plug":
        return PLUG.run(cls, kappa, lam, d, nprime, n=n, sig=sig, theta=theta)[0]
    if config == "constC":
        return np.log(run_const_c(cls, kappa, lam, d, nprime, c_set, n=n,
                                  sig=sig, theta=theta)[0])
    raise ValueError(config)
'''))

cells.append(code(r'''"""The pymrm bed: one class, pseudo-steady (eps = 0) or full transient."""
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from pymrm import (construct_convflux_upwind, construct_div, compute_boundary_values,
                   interp_cntr_to_stagg_tvd, vanleer)

class Bed:
    """eps dc/dtheta + dc/ds = -kappa c^n a  on s in [0,1], s = w/W."""

    def __init__(self, m=200, tvd=True):
        self.m, self.tvd = m, tvd
        self.s_f = np.linspace(0.0, 1.0, m + 1)
        self.s_c = 0.5 * (self.s_f[:-1] + self.s_f[1:])
        # bc on the OUTWARD normal: inlet  {a:0,b:1,d:1} -> c = 1 (feed)
        #                           outlet {a:1,b:0,d:0} -> dc/dn = 0 (outflow)
        self.bc = ({"a": 0.0, "b": 1.0, "d": 1.0}, {"a": 1.0, "b": 0.0, "d": 0.0})
        conv, conv_bc = construct_convflux_upwind((m, 1), self.s_f, self.s_c,
                                                 self.bc, v=1.0)
        div = construct_div((m, 1), self.s_f, nu=0)   # nu=0: Cartesian weight coord
        self.div = div
        self.Lop = (div @ conv).tocsc()               # assembled ONCE
        self.bvec = np.asarray((div @ conv_bc).todense()).ravel()

    def _dcorr(self, c):
        """van Leer deferred correction of the first-order upwind flux."""
        if not self.tvd:
            return np.zeros(self.m)
        _, dg = interp_cntr_to_stagg_tvd(c.reshape(-1, 1), self.s_f, self.s_c,
                                         self.bc, 1.0, tvd_limiter=vanleer, axis=0)
        return np.asarray(self.div @ np.asarray(dg).reshape(-1, 1)).ravel()

    def _cascade_solve(self, coef, order, tol=1e-13, it=80):
        """Newton on  dy/ds + coef y^order = 0, y(0) = 1 (fluid or poison)."""
        y = np.ones(self.m)
        for _ in range(it):
            J = (self.Lop + diags(coef * order * _pow(y, order - 1.0))).tocsc()
            r = self.Lop @ y + self.bvec + coef * _pow(y, order) + self._dcorr(y)
            dy = splu(J).solve(-r)
            y = np.clip(y + dy, 1e-14, None)
            if np.max(np.abs(dy)) < tol:
                break
        return y

    def steady_c(self, kappa, a, n=1.0):
        return self._cascade_solve(kappa * a, n)

    def steady_p(self, g, a, d, nprime):
        return np.ones(self.m) if g == 0.0 else self._cascade_solve(g * _pow(a, d),
                                                                    nprime)

    def outlet(self, c):
        _, _, v_bc, _ = compute_boundary_values(c.reshape(-1, 1), self.s_f,
                                                self.s_c, self.bc, axis=0)
        return float(np.ravel(v_bc)[0])

    def phi(self, cls, c, p):
        return {"parallel": c, "series": 1.0 - c,
                "side-by-side": p}.get(cls, np.ones_like(c))

    def march(self, cls, kappa, lam, d, nprime, n=1.0, sig=0.0, nt=200, theta=TH,
              scheme="heun"):
        """Pseudo-steady (eps = 0): fluid solved at each instant, activity
        marched with Heun (RK2; 'euler' for the break table).  The marched
        variable is w = (a^(1-d)-1)/(1-d), the paper's own linearising
        substitution, so no a^d appears in the increment."""
        dt = float(theta[-1]) / nt
        def rate(w):
            a = a_of_w(w, d)
            c = self.steady_c(kappa, a, n)
            p = (self.steady_p(sig * lam, a, d, nprime) if cls == "side-by-side"
                 else np.ones(self.m))
            f = np.clip(self.phi(cls, c, p), 0.0, None)
            return -lam * _pow(f, nprime), c, a
        w = np.zeros(self.m)
        f1, c, a = rate(w)
        ts, Xs, Aall = [0.0], [1.0 - self.outlet(c)], [a.copy()]
        for j in range(nt):
            if scheme == "euler":
                w = w + dt * f1
            else:
                f2, _, _ = rate(w + dt * f1)
                w = w + 0.5 * dt * (f1 + f2)
            f1, c, a = rate(w)
            ts.append((j + 1) * dt); Xs.append(1.0 - self.outlet(c))
            Aall.append(a.copy())
        return np.interp(theta, ts, Xs), np.array(Aall).T, np.array(ts)

    def jac_sparsity(self):
        """Pattern of d(rhs)/dy for the transient state y = [c(0..m-1), a(0..m-1)].
        The c row reads c_{i-2..i+1} (upwind + van Leer) and a_i; the a row reads
        c_i and a_i.  Handing this to BDF turns a 2m x 2m dense numerical
        Jacobian into a handful of coloured differences."""
        from scipy.sparse import lil_matrix
        S = lil_matrix((2 * self.m, 2 * self.m), dtype=bool)
        for i in range(self.m):
            for j in range(max(0, i - 2), min(self.m, i + 2)):
                S[i, j] = True
            S[i, self.m + i] = True
            S[self.m + i, i] = True
            S[self.m + i, self.m + i] = True
        return S.tocsr()

    def transient(self, cls, kappa, lam, d, nprime, eps, n=1.0, sig=0.0,
                  theta=TH, rtol=1e-8):
        """The paper's pseudo-steady assumption DROPPED: eps dc/dtheta retained."""
        def rhs(_, y):
            c = np.clip(y[:self.m], 1e-14, None)
            a = a_of_w(y[self.m:], d)
            rc = -(self.Lop @ c + self.bvec + kappa * a * _pow(c, n)
                   + self._dcorr(c)) / eps
            p = (self.steady_p(sig * lam, a, d, nprime) if cls == "side-by-side"
                 else np.ones(self.m))
            f = np.clip(self.phi(cls, c, p), 0.0, None)
            return np.concatenate([rc, -lam * _pow(f, nprime)])
        c0 = self.steady_c(kappa, np.ones(self.m), n)
        sol = solve_ivp(rhs, (0.0, float(theta[-1])),
                        np.concatenate([c0, np.zeros(self.m)]), t_eval=theta,
                        rtol=rtol, atol=1e-11, method="BDF",
                        jac_sparsity=self._SP)
        return (np.array([1.0 - self.outlet(sol.y[:self.m, j])
                          for j in range(theta.size)]),
                a_of_w(sol.y[self.m:, :], d))

BED = Bed(m=200)
BED._SP = BED.jac_sparsity()
print(f"pymrm bed assembled: {BED.m} cells, "
      f"{BED.Lop.nnz} nonzeros in the convection-divergence operator")
'''))

# ------------------------------------------------------------------ results -
cells.append(md(r"""## Results

### 1. The printed structure, proved

A symbolic identity per printed step — the count is printed by the cell below,
not typed here. Each equation is re-derived from the numbered
predecessor the paper derives it from, so the transcription is checked line by
line rather than retyped — a mis-read subscript or prime surfaces as a nonzero
residual instead of as a plausible wrong equation. These are **identities**:
once the transcription is right they cannot fail, and they are labelled that
way in the coverage map. What they protect against is the transcription
itself, which on this scan is the real risk."""))

cells.append(code(r'''"""Symbolic verification of eqs. (8), (12)-(15), (17), (19), (21)-(23),
(25)-(32) and the slope/intercept printed on Fig. 6."""
t = sp.Symbol("t", nonnegative=True)
k, kd, kpp, kp, kdp = sp.symbols("k k_d k'' k' k_d'", positive=True)
CA0, CA, CAinf, CR, Q = sp.symbols("C_A0 C_A C_Ainf C_R Q", positive=True)
tau = sp.Symbol("tau'", positive=True)
npr = sp.Symbol("n'", nonnegative=True)
dd = sp.Symbol("d", positive=True)
aF = sp.Function("a")(t)
av = sp.Symbol("a", positive=True)
SYM = {}

# (8) -> (7): C_A + C_R is constant for a specific feed
SYM["eq8_reduces_to_eq7"] = sp.simplify(
    (kd * (CA + CR) ** npr).subs(CR, CA0 - CA) - kd * CA0 ** npr)

# (11) -> (12) -> (13) -> (14) -> (15): the batch-fluid chain
a12 = sp.dsolve(sp.Derivative(aF, t) + kd * aF, aF, ics={aF.subs(t, 0): 1}).rhs
SYM["eq12"] = sp.simplify(a12 - sp.exp(-kd * t))
CAt = sp.dsolve(sp.Derivative(sp.Function("C")(t), t)
                + kpp * sp.exp(-kd * t) * sp.Function("C")(t),
                sp.Function("C")(t),
                ics={sp.Function("C")(0): CA0}).rhs
SYM["eq13"] = sp.simplify(sp.log(CA0 / CAt) - kpp / kd * (1 - sp.exp(-kd * t)))
SYM["eq14"] = sp.simplify(sp.limit(sp.log(CA0 / CAt), t, sp.oo) - kpp / kd)
CAinf_e = CA0 * sp.exp(-kpp / kd)
SYM["eq15"] = sp.simplify(sp.expand_log(
    sp.log(sp.log(CAt / CAinf_e)) - sp.log(sp.log(CA0 / CAinf_e)) + kd * t,
    force=True))

# (16) -> (17) -> (18) -> (19): mixed flow at constant flow rate
XA = 1 - CA / CA0
SYM["eq17"] = sp.simplify(CA0 * (XA / (k * av * CA)) * k * av + 1 - CA0 / CA)
SYM["eq19"] = sp.simplify(sp.log((1 + k * sp.exp(-kd * t) * tau) - 1)
                          - (sp.log(k * tau) - kd * t))

# (20) -> (21) -> (22): plug flow at constant flow rate
X = sp.Symbol("X", nonnegative=True)
e20 = sp.integrate(1 / (k * av * CA0 * (1 - X)), (X, 0, XA))
SYM["eq21"] = sp.simplify(CA0 * e20 - sp.log(CA0 / CA) / (k * av))
tau21 = sp.log(CA0 / CA) / (k * sp.exp(-kd * t))
SYM["eq22"] = sp.simplify(sp.log(sp.log(CA0 / CA)) - (sp.log(k * tau21) - kd * t))

# (23): mixed flow with C_A held constant
tau18 = (CA0 / CA - 1) / (k * sp.exp(-kd * t))
SYM["eq23"] = sp.simplify(sp.log(tau18)
                          - (kd * t + sp.log((CA0 - CA) / (k * CA))))

# (25)-(27) and (28)-(31): the general-order chain at constant C_A
# (27): tau' = C_A0 (W/F_A0) with W/F_A0 = X_A/(-r'_A) and -r'_A = k' a  (eq. 25)
SYM["eq27"] = sp.simplify(CA0 * ((1 - CA / CA0) / (kp * av)) - (CA0 - CA) / (kp * av))
a_d = {0: 1 - kdp * t, 1: sp.exp(-kdp * t), 2: 1 / (1 + kdp * t),
       3: 1 / sp.sqrt(1 + 2 * kdp * t)}
for j, expr in a_d.items():
    SYM[f"eq26_solution_d{j}"] = sp.simplify(sp.diff(expr, t) + kdp * expr ** j)
    SYM[f"eq26_ic_d{j}"] = sp.simplify(expr.subs(t, 0) - 1)
taud = {j: Q / e for j, e in a_d.items()}          # Q = (C_A0 - C_A)/k' , eq. 27
SYM["eq28"] = sp.simplify(1 / taud[0] - (kp / (CA0 - CA)
                                         - kp * kdp / (CA0 - CA) * t)
                          ).subs(Q, (CA0 - CA) / kp)
SYM["eq28"] = sp.simplify(SYM["eq28"])
SYM["eq29"] = sp.simplify((sp.log(taud[1]) - (sp.log((CA0 - CA) / kp) + kdp * t)
                           ).subs(Q, (CA0 - CA) / kp))
SYM["eq30"] = sp.simplify((taud[2] - ((CA0 - CA) / kp + (CA0 - CA) * kdp / kp * t)
                           ).subs(Q, (CA0 - CA) / kp))
SYM["eq31"] = sp.simplify((taud[3] ** 2 - (((CA0 - CA) / kp) ** 2
                                           + ((CA0 - CA) / kp) ** 2 * 2 * kdp * t)
                           ).subs(Q, (CA0 - CA) / kp))

# (32) and Fig. 6, for GENERAL d - via the substitution the paper itself uses
SUB = {sp.Derivative(aF, t): -kdp * aF ** dd}
SYM["eq32_substitution"] = sp.simplify(sp.powsimp(
    sp.diff(aF ** (1 - dd), t).subs(SUB) - (dd - 1) * kdp, force=True))
SYM["fig6_slope_general_d"] = sp.simplify(sp.powsimp(
    sp.diff((Q / aF) ** (dd - 1), t).subs(SUB) - Q ** (dd - 1) * (dd - 1) * kdp,
    force=True))
SYM["fig6_intercept_general_d"] = sp.simplify(sp.powsimp(
    ((Q / aF) ** (dd - 1)).subs(aF, 1) - Q ** (dd - 1), force=True))

bad = [kk for kk, vv in SYM.items() if sp.simplify(vv) != 0]
print(f"symbolic identities checked: {len(SYM)};  nonzero residuals: {bad if bad else 'none'}")
M["printed_chain_identities"] = float(len(SYM))
M["printed_chain_nonzero"] = float(len(bad))
assert not bad, bad
'''))

cells.append(md(r"""### 2. Eq. (24) is wrong as printed, by exactly $\ln(k_d/k)$

The paper introduces eq. (24) with *"At any instant in the plug flow reactor
Eq. 21 applies. Also noting that $\tau'$ and $t$ are the two variables we
obtain on suitable rearrangement"*. Eq. (21) is
$\tau' = \dfrac{1}{k\,e^{-k_dt}}\ln\dfrac{C_{\rm A0}}{C_{\rm A}}$, so the
rearrangement is forced. Two independent arguments settle it, and the second
does not use the first:"""))

cells.append(code(r'''"""Eq. (24): printed vs required by the paper's own eq. (21) -- and the units."""
printed_24  = kd * t + sp.log(sp.log(CA0 / CA) / kd)     # AS PRINTED, p. 270
required_24 = kd * t + sp.log(sp.log(CA0 / CA) / k)      # from eq. (21)
gap = sp.simplify(sp.log(tau21) - printed_24)
M["eq24_residual_is_log_kd_over_k"] = float(sp.simplify(gap - sp.log(kd / k)) == 0)
M["eq24_repaired_residual"] = float(sp.simplify(sp.log(tau21) - required_24) == 0)
print(f"eq. (21) minus eq. (24) as printed  =  {gap}")
print(f"eq. (21) minus eq. (24) with 1/k    =  {sp.simplify(sp.log(tau21) - required_24)}")

# --- the dimensional argument, from the paper's own Nomenclature -------------
# k   : liter^n / (mol^(n-1) sec g cat)          -> [k]  = L^n mol^(1-n) s^-1 g^-1
# k_d : liter^n' / (sec mol^n')                  -> [k_d] = L^n' mol^-n' s^-1
# tau': g cat sec / liter of fluid
# At n = 1, n' = 0:  1/k = g s / L  = [tau']   while  1/k_d = s.
dims = pd.DataFrame(
    [["1/k   (eq. 21 requires)", "g cat sec / liter", "yes"],
     ["1/k_d (eq. 24 prints)",   "sec",               "no"],
     ["tau'  (the quantity)",    "g cat sec / liter", "-"]],
    columns=["quantity", "units at n = 1, n' = 0 (Nomenclature, p. 265)",
             "can sit inside ln tau'?"])
display(Markdown("**Eq. (24)'s printed argument is dimensionally impossible:**\n\n"
                 + md_table(dims)))

# --- and numerically, on the device eq. (24) is the test plot for ------------
_k, _kd, _CA0, _CA = 2.0, 0.7, 1.0, 0.35
off = np.log(np.log(_CA0 / _CA) / _kd) - np.log(np.log(_CA0 / _CA) / _k)
M["eq24_numeric_offset"] = float(off)
M["eq24_offset_expected"] = float(np.log(_k / _kd))
print(f"\nintercept as printed minus intercept required, at k = {_k}, k_d = {_kd}: "
      f"{off:.6f}  vs  ln(k/k_d) = {np.log(_k/_kd):.6f}")
print("The slope of Fig. 5 is unaffected; only the intercept - i.e. only the "
      "recovered value of k - is. REPORTED, NOT REPAIRED: the CSV and the table "
      "above carry the printed form.")
BREAKS.append(("eq24_numeric_offset", off, "read eq. (24) as if it printed 1/k",
               0.0, "the defect metric is exactly the printed-vs-required gap"))
'''))

cells.append(md(r"""### 3. The four classes in the four devices

Every class solved in every contacting, at the base orders
$(n, n', d) = (1, 1, 1)$. The curves differ — visibly. The rest of this page
is about how much of that difference survives a fit."""))

cells.append(code(r'''"""Figure 1: the observable of each device, for each class, at the base point."""
RESP = {}
for cfg in ("batch", "plug", "mixed", "constC"):
    for cls in CLASSES:
        RESP[(cfg, cls)] = respond(cfg, cls, BASE["kappa"], BASE["lam"], BASE["d"],
                                   BASE["nprime"], BASE["sig"], n=BASE["n"],
                                   c_set=BASE["c_set"])
fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.0), sharex=True)
for ax, cfg in zip(axes.ravel(), ("batch", "plug", "mixed", "constC")):
    for cls in CLASSES:
        ax.plot(TH, RESP[(cfg, cls)], color=COL[cls], lw=1.6, label=cls)
    ax.set_title(DEV[cfg], fontsize=9)
    ax.set_ylabel(r"$\ln\,\tau'/\tau'_0$" if cfg == "constC" else r"conversion $X_A$")
    ax.grid(alpha=0.25, lw=0.5)
axes[1, 0].set_xlabel(r"$\theta = t/t_{\rm obs}$")
axes[1, 1].set_xlabel(r"$\theta = t/t_{\rm obs}$")
axes[0, 0].legend(fontsize=8, frameon=False)
fig.suptitle("The four deactivation classes of eqs. (4)-(7), in the paper's own "
             "batch-solids contactings", fontsize=10)
fig.tight_layout(); plt.show()
for cfg in ("batch", "plug", "mixed", "constC"):
    for cls in CLASSES:
        M[f"span_{cfg}_{SHORT[cls]}"] = float(np.ptp(RESP[(cfg, cls)]))
print("observable spans (the signal a fit has to reproduce):")
print(pd.DataFrame({cfg: {cls: np.ptp(RESP[(cfg, cls)]) for cls in CLASSES}
                    for cfg in ("batch", "plug", "mixed", "constC")}).round(4)
      .to_string())
'''))

cells.append(md(r"""### 4. Five exact degeneracies — and they reach every device

Before any optimiser runs: **every one of the four devices collapses some pair
of classes onto another identically**, and each collapse is a one-line
algebraic fact about the device, not a numerical coincidence. Each is stated,
proved symbolically, and then confirmed by running the two "different" classes
and subtracting.

1. **$n' = 0$ collapses all four, everywhere.** With no concentration in the
   decay law, eqs. (4)–(6) *are* eq. (7). This is the paper's own statement
   (p. 267, *"When concentration independent (or $n' = 0$) then any type of
   batch-solids system may be used"*), and it is why the taxonomy has a
   fourth member at all.
2. **Mixed flow at constant rate: series *is* parallel, relabelled.** The
   performance relation (16) is algebraic, so
   $C_{\rm R}/C_{\rm A0} = \kappa\,(C_{\rm A}/C_{\rm A0})^{n} a$ *exactly* at
   every instant. Hence
   $C_{\rm R}^{\,n'} \propto C_{\rm A}^{\,n n'} a^{\,n'}$, and series with
   orders $(n', d)$ is **identical** to parallel with orders
   $(n\,n',\, d+n')$ and $k_d \to k_d \kappa^{n'}$. Every observable, at every
   time. A well-mixed fluid cannot tell the two apart — ever.
3. **Held-$C_{\rm A}$ mixed flow: everything is independent deactivation.**
   That is what "decoupling" means. With $C_{\rm A}$ (hence $C_{\rm R}$) held
   fixed, $k'_d = k_dC_{\rm A}^{n'}$ or $k_dC_{\rm R}^{n'}$ or $k_d$ are all
   just *constants*, and eq. (26) is the same equation for all three classes.
   At $d = 1$ side-by-side joins them, because the poison consumed per unit
   time then goes as $a^{d-1}$ = const and $C_{\rm P}$ settles at a fixed
   value: **all four classes are one equation there.**
   **The device Levenspiel recommends is the one in which, within a single
   run, the classes are least distinguishable.**
4. **Side-by-side *is* parallel wherever the poison profile equals the
   reactant profile.** The poison obeys the same first-order form as the
   reactant — $-dp/d(\cdot) = \sigma\lambda\,a^{\,d}p^{\,n'}$ against
   $-dc/d(\cdot) = \kappa\,a\,c^{\,n}$ — so at $d = 1$, $n' = n$ and
   $\sigma = \kappa/\lambda$ the two are the *same* differential equation with
   the same initial value: $p \equiv c$, identically, and eq. (6) becomes
   eq. (4) with the same constants. This holds in the batch, in the
   **plug-flow bed** and in the constant-flow mixed reactor alike — it is a
   property of the rate forms, not of the contacting.
5. **And with $\sigma$ free: side-by-side is parallel with a relabelled
   $n'$.** At $n = n' = d = 1$ the two profiles are powers of one another —
   $p = c^{\lambda\sigma/\kappa}$ in the bed (both are $\exp(-\text{const}
   \int a)$), and $c = p^{\kappa/(\lambda\sigma)}$ in a batch of fluid (where
   $p = 1-\sigma(1-a)$ is exactly affine in the activity). So the side-by-side
   driver is a *power of the parallel driver*, and side-by-side with **any**
   $\sigma$ is parallel with $n' \to \lambda\sigma/\kappa$, same
   $\kappa,\lambda,d$. Degeneracy 4 is the equal-parameter slice of this
   family ($\sigma = \kappa/\lambda$ makes the relabelled order 1 again).
   **This one holds in the batch and the bed and nowhere else** — in the two
   well-mixed devices the profiles are not powers of each other, and the
   measurement below shows the identity failing there by 3e-2 and 1e-1.

**So the plug-flow bed is not the one device without an exact degeneracy.** It
carries degeneracies 1, 4 and 5, exactly as the batch of fluid does. What
separates the bed and the batch from the two well-mixed devices is narrower and
still real: they are the two devices in which **parallel, series and
independent remain algebraically distinct**. What separates the bed from the
batch is narrower again, and is what §6 is about: it is the only contacting
with a spatial coordinate, so it is the only one whose activity carries a
*profile*.

Each identity is checked **on** its condition and **off** it."""))

cells.append(code(r'''"""The degeneracies: predicted from algebra, then measured by subtraction."""
rows = []

# (1) n' = 0 collapses the four classes in every device
for cfg in ("batch", "plug", "mixed", "constC"):
    ref = respond(cfg, "independent", 1.5, 3.0, 1.0, 0.0, 0.0, c_set=BASE["c_set"])
    worst = 0.0
    for cls in ("parallel", "series", "side-by-side"):
        y = respond(cfg, cls, 1.5, 3.0, 1.0, 0.0, 0.0, c_set=BASE["c_set"])
        worst = max(worst, float(np.max(np.abs(y - ref))))
    M[f"degen_nprime0_{cfg}"] = worst
    rows.append([f"n' = 0: all four classes identical", DEV[cfg], f"{worst:.2e}"])
# ... and the control: give the decay a concentration order and it stops holding
worst_np0 = 0.0
for cfg in ("batch", "plug", "mixed", "constC"):
    ref = respond(cfg, "independent", 1.5, 3.0, 1.0, 0.0, 0.0, c_set=BASE["c_set"])
    y = respond(cfg, "parallel", 1.5, 3.0, 1.0, 0.7, 0.0, c_set=BASE["c_set"])
    worst_np0 = max(worst_np0, float(np.max(np.abs(y - ref))))
M["degen_nprime0_off_condition"] = worst_np0
rows.append(["...the same comparison at n' = 0.7", "worst of the four",
             f"{worst_np0:.2e}"])

# (2) mixed flow: series(n', d) == parallel(n n', d + n'), lam -> lam kappa^n'
def _degen2(rtol):
    worst = 0.0
    for nn in (0.5, 1.0, 2.0):
        for npv, dv, lm, kp_ in ((1.0, 1.0, 3.0, 1.5), (0.5, 2.0, 2.0, 4.0),
                                 (2.0, 0.5, 1.0, 0.8)):
            ys = run_mixed("series", kp_, lm, dv, npv, n=nn, rtol=rtol)[0]
            yp = run_mixed("parallel", kp_, lm * kp_ ** npv, dv + npv, nn * npv,
                           n=nn, rtol=rtol)[0]
            worst = max(worst, float(np.max(np.abs(ys - yp))))
    return worst
worst2 = _degen2(1e-9)                 # the production tolerance
M["degen_mixed_series_is_parallel"] = worst2
M["degen_mixed_series_is_parallel_rtol1e-13"] = _degen2(1e-13)
rows.append(["series(n',d) == parallel(n n', d+n')", DEV["mixed"], f"{worst2:.2e}"])

# (3) held-C_A device: parallel == series == independent after relabelling lam,
#     and at d = 1 side-by-side joins them (consumption ~ a^(d-1) = const, so
#     C_P settles and its driver is a constant like the other three)
cs = BASE["c_set"]; worst3 = 0.0
for npv, dv in ((1.0, 1.0), (2.0, 3.0), (0.5, 0.0)):
    t1 = np.log(run_const_c("parallel", 1.5, 3.0, dv, npv, cs)[0])
    t2 = np.log(run_const_c("series", 1.5, 3.0 * (cs / (1 - cs)) ** npv, dv, npv, cs)[0])
    t3 = np.log(run_const_c("independent", 1.5, 3.0 * cs ** npv, dv, npv, cs)[0])
    worst3 = max(worst3, float(max(np.max(np.abs(t1 - t2)), np.max(np.abs(t1 - t3)))))
M["degen_constC_par_ser_ind"] = worst3
rows.append(["parallel == series == independent", DEV["constC"], f"{worst3:.2e}"])

# (4) side-by-side(sigma = kappa/lambda, n' = n, d = 1) == parallel, IDENTICAL
#     parameters, because the poison then solves the reactant's own equation.
#     Checked in every device that has one, not just the batch.
KA, LA = 1.5, 3.0
famA = {}
for cfg, runner in (("batch", run_batch), ("plug", PLUG.run),
                    ("mixed", run_mixed)):
    on = 0.0
    for nn, dv in ((1.0, 1.0), (2.0, 1.0), (1.5, 1.0)):     # n' = n on condition
        yp = runner("parallel", KA, LA, dv, nn, n=nn)[0]
        ys = runner("side-by-side", KA, LA, dv, nn, n=nn, sig=KA / LA)[0]
        on = max(on, float(np.max(np.abs(yp - ys))))
    famA[cfg] = on
    M[f"degen_sbs_is_parallel_famA_{cfg}"] = on
    rows.append(["side-by-side(sig=kappa/lam) == parallel, n'=n and d=1",
                 DEV[cfg], f"{on:.2e}"])
off_cond = 0.0
for nn, npv, dv in ((1.0, 1.0, 2.0), (1.0, 2.0, 1.0)):      # d != 1, or n' != n
    yp = run_batch("parallel", KA, LA, dv, npv, n=nn)[0]
    ys = run_batch("side-by-side", KA, LA, dv, npv, n=nn, sig=KA / LA)[0]
    off_cond = max(off_cond, float(np.max(np.abs(yp - ys))))
M["degen_batch_sbs_is_parallel_on_condition"] = famA["batch"]
M["degen_batch_sbs_off_condition"] = off_cond
rows.append(["...the same comparison OFF that condition (d = 2, or n' != n)",
             DEV["batch"], f"{off_cond:.2e}"])

# (5) sigma FREE: at n = n' = d = 1 the poison profile is a POWER of the
#     reactant profile, so side-by-side(sigma) == parallel(n' -> lam sig/kap).
#     Bed:   p = exp(-sig lam I), c = exp(-kap I)  ->  p = c^(sig lam/kap).
#     Batch: p = 1 - sig(1-a) and dln c/dln p = kap/(lam sig) -> c = p^(kap/lam sig).
#     True in those two devices ONLY; the control below is the other two.
SIG_FAM = (0.2, 0.5, 0.8, 1.5, 3.0)
famB, famB_rows = {}, []
for cfg in ("batch", "plug", "mixed", "constC"):
    worst = 0.0
    for sv in SIG_FAM:
        ys = respond(cfg, "side-by-side", KA, LA, 1.0, 1.0, sv, n=1.0,
                     c_set=BASE["c_set"])
        yp = respond(cfg, "parallel", KA, LA, 1.0, LA * sv / KA, 0.0, n=1.0,
                     c_set=BASE["c_set"])
        worst = max(worst, float(np.max(np.abs(ys - yp))))
    famB[cfg] = worst
    M[f"degen_sbs_is_parallel_famB_{cfg}"] = worst
famB_off = 0.0                       # off the n = n' = d = 1 condition, in the bed
for nn, npv, dv in ((1.0, 0.5, 1.0), (1.0, 2.0, 1.0), (1.0, 1.0, 2.0), (2.0, 1.0, 1.0)):
    ys = PLUG.run("side-by-side", KA, LA, dv, npv, n=nn, sig=BASE["sig"])[0]
    yp = PLUG.run("parallel", KA, LA, dv, npv * LA * BASE["sig"] / KA, n=nn)[0]
    famB_off = max(famB_off, float(np.max(np.abs(ys - yp))))
M["degen_sbs_famB_off_condition"] = famB_off
for cfg in ("batch", "plug"):
    rows.append([f"side-by-side(any sig) == parallel(n' = lam sig/kap), n=n'=d=1",
                 DEV[cfg], f"{famB[cfg]:.2e}"])
rows.append(["...the same identity in the two WELL-MIXED devices (it fails)",
             "mixed / constC",
             f"{famB['mixed']:.2e} / {famB['constC']:.2e}"])
rows.append(["...and in the bed OFF the n = n' = d = 1 condition", DEV["plug"],
             f"{famB_off:.2e}"])

display(Markdown("**Exact degeneracies, measured as max |difference| in the "
                 "device's observable:**\n\n"
                 + md_table(pd.DataFrame(rows, columns=["identity", "device",
                                                        "max |difference|"]))))
ratio = "infinitely" if famA["batch"] == 0.0 else f"{off_cond / famA['batch']:.1e}x"
print(f"\nControls. Degeneracy 4 off its condition (d = 2, or n' != n) separates "
      f"by {off_cond:.2e} - {ratio} larger than on it ({famA['batch']:.2e}). "
      f"Degeneracy 5 holds in the batch ({famB['batch']:.1e}) and the bed "
      f"({famB['plug']:.1e}) over sigma = {SIG_FAM[0]:g}-{SIG_FAM[-1]:g}, and "
      f"FAILS in the mixed reactor ({famB['mixed']:.1e}) and the held-C_A device "
      f"({famB['constC']:.1e}) - {famB['mixed']/max(famB['plug'],1e-300):.0e}x "
      f"larger - and fails in the bed itself off n = n' = d = 1 "
      f"({famB_off:.1e}). The identities are sharp, not approximate.")
print(f"\nDegeneracy 2 is reported at the production integrator tolerance "
      f"(rtol 1e-9): {worst2:.2e}. It is a TOLERANCE, not a separation - "
      f"re-integrating the same nine order combinations at rtol 1e-13 gives "
      f"{M['degen_mixed_series_is_parallel_rtol1e-13']:.2e}, i.e. the identity "
      f"is exact and the printed number measures the ODE solver.")
BREAKS.append(("degen_batch_sbs_is_parallel_on_condition", famA["batch"],
               "step off the n' = n, d = 1 condition (d = 2)", off_cond,
               "shows the identity is conditional, not a numerical artefact"))
BREAKS.append(("degen_sbs_is_parallel_famB_plug", famB["plug"],
               "evaluate the same sigma-family identity in the MIXED reactor",
               famB["mixed"],
               "degeneracy 5 is a property of the batch and the bed only"))
BREAKS.append(("degen_mixed_series_is_parallel", worst2,
               "tighten the integrator from rtol 1e-9 to 1e-13",
               M["degen_mixed_series_is_parallel_rtol1e-13"],
               "the metric measures the ODE tolerance; the identity itself is exact"))
'''))

cells.append(md(r"""### 5. Separability, measured

The only honest test of a discriminability claim is to **inject the wrong
class and see whether anything moves**. For every device and every pair
$(T, C)$ of classes, the response generated by $T$ is handed to a fitter that
is free to choose $C$'s parameters — $\kappa$, $\lambda$, $d$, $n'$ and, for
side-by-side, $\sigma$ — and the residual it cannot remove is reported.

$D(T\!\to\!C)$ is therefore **the RMS conversion resolution an experiment
needs before it can refute $C$**. Read it against your own instrument: 2e-3
is a good laboratory chromatograph, 1e-4 is not achievable.

**An optimiser bound is not a measurement, and where §4 supplies an exact twin
the twin is used instead.** Each of the five degeneracies names a *specific*
parameter vector at which some candidate reproduces some truth exactly — $n'=0$
for degeneracy 1, $(nn', d+n')$ with $\lambda\kappa^{n'}$ for degeneracy 2, and
so on. Those points are inside the fitter's own bounds and next to its own seed
grid, and yet the optimiser does not always find them: with the fit alone, a
number of entries whose exact minimum is **zero** came out well above it — the
count and the size are printed below, not typed here.
So every entry is evaluated **both** ways — the multistart fit, and the
identity twin where §4 predicts one — and the matrix reports the smaller,
flagged as an identity. An entry with no predicted twin is the fit alone and is
an upper bound; the table below says which is which.

**The diagonal's exact minimum is zero by construction**, since it is a class
fitted to its own response with every parameter free. Any non-zero diagonal is
therefore an **optimiser stall**, not a null and not a result, and the twin
(the truth itself) sets those entries to zero exactly. The size of the stall is
still worth reporting — it is what the fit alone would have claimed — and its
recovered parameters are printed below, because the interesting question is
*where* the optimiser stalled.

**A floor of twenty times that stall would be a natural extra safeguard, and it
is not applied here, because it never binds.** Twenty times the worst stall is
below the instrument resolution on all four devices, so such a threshold would
equal the instrument resolution exactly on every one of them. The check is
printed below rather than asserted in prose. Nothing that follows is claimed
inside the machinery's own noise — but the constraint that guarantees it is the
instrument, not the fitter, and the page says which."""))

cells.append(code(r'''"""The separability matrices.  Two-stage multistart: screen every seed cheaply,
then polish the best two.  Deterministic - a fixed seed list, no randomness."""
def _unpack(cls, x):
    kappa, lam, d = np.exp(x[0]), np.exp(x[1]), x[2]
    if cls == "independent":
        return kappa, lam, d, 0.0, 0.0
    if cls == "side-by-side":
        return kappa, lam, d, x[3], np.exp(x[4])
    return kappa, lam, d, x[3], 0.0

def _bounds(cls):
    # deliberately physical: d and n' are ORDERS in eqs. (4)-(7), kappa and lam
    # are a Damkoehler number and a decay-per-window.  Wider boxes only add
    # pathological corners where the ODE is expensive and no fit ever lands.
    lo, hi = [np.log(5e-2), np.log(5e-2), 0.0], [np.log(3e1), np.log(6e1), 4.0]
    if cls != "independent":
        lo.append(0.0); hi.append(3.0)
    if cls == "side-by-side":
        lo.append(np.log(1e-2)); hi.append(np.log(2e1))
    return np.array(lo), np.array(hi)

def _seeds(cls, kappa0, lam0):
    """A fixed, class-independent seed grid over the two ORDERS - the parameters
    the landscape is nonconvex in.  It contains the base point (d, n') = (1, 1)
    for every candidate alike, which is what makes the diagonal (the fitter
    recovering the truth) a meaningful null rather than a lucky start."""
    out = []
    for dv in (0.5, 1.0, 2.0, 3.0):
        base = [np.log(kappa0), np.log(lam0), dv]
        if cls == "independent":
            out.append(base)
        elif cls == "side-by-side":
            # sigma is a smooth scale, not an order: one seed, found by descent
            out += [base + [npv, np.log(1.0)] for npv in (0.5, 1.0, 2.0)]
        else:
            out += [base + [npv] for npv in (0.5, 1.0, 2.0)]
    return out

class _Budget(Exception):
    """Raised inside the residual to cap the cost of one fit."""

CAPPED = []          # entries whose evaluation budget bound

def fit_class(config, target, cls, kappa0=1.5, lam0=3.0, n=1.0, c_set=0.4,
              theta=TH, nfev=300, keep=2, rounds=5, seeds=None, budget=3500,
              tag=None):
    """Best fit of one class to one response.  Deterministic: a fixed seed grid
    over the two ORDERS (the parameters the landscape is nonconvex in), ranked
    by a direct residual evaluation, then the best `keep` polished with
    x_scale='jac' - which is what gets the diagonal of the matrix (the fitter
    recovering the truth) down to machine zero instead of 1e-6.

    `budget` caps the residual evaluations of ONE fit.  It binds only where the
    minimum is an exact degeneracy the optimiser is still chasing below 1e-9;
    the best value seen is kept, so a capped fit can only report a residual that
    is too LARGE.  Capped entries are listed, and the summary below refuses to
    call anything separable on one."""
    lo, hi = _bounds(cls)
    state = {"n": 0, "best": np.inf, "capped": False, "bx": None}
    def res(x):
        state["n"] += 1
        if state["n"] > budget:
            state["capped"] = True
            raise _Budget
        kk, ll, dv, npv, sg = _unpack(cls, x)
        try:
            y = respond(config, cls, kk, ll, dv, npv, sg, n=n, c_set=c_set,
                        theta=theta)
        except Exception:
            return np.full(target.size, 1e2)
        r = (y - target) if np.all(np.isfinite(y)) else np.full(target.size, 1e2)
        v_ = float(np.sqrt(np.mean(r ** 2)))
        if v_ < state["best"]:          # keep the ARGMIN, not just the minimum:
            state["best"] = v_          # the recovered parameters are a result
            state["bx"] = np.array(x, float)
        return r
    grid = [np.clip(np.asarray(x0, float), lo + 1e-9, hi - 1e-9)
            for x0 in (seeds or _seeds(cls, kappa0, lam0))]
    try:
        ranked = sorted(grid, key=lambda x0: float(np.sqrt(np.mean(res(x0) ** 2))))
    except _Budget:
        ranked = grid
    bx = ranked[0]
    for x0 in ranked[:keep]:
        x, prev = x0, np.inf
        for _ in range(rounds):
            # a trust-region solve can stall well short of the minimum; restart
            # it from its own answer while it keeps gaining.  This is what takes
            # the exact degeneracies from ~1e-6 to machine zero.
            try:
                r = least_squares(res, x, bounds=(lo, hi), xtol=1e-15, ftol=1e-15,
                                  gtol=1e-15, max_nfev=nfev, x_scale="jac")
            except _Budget:
                break
            except Exception:
                break
            v = float(np.sqrt(np.mean(r.fun ** 2)))
            x = r.x
            if v <= state["best"] + 1e-300:
                bx = r.x
            if v < 1e-13 or v > 0.9 * prev:
                break
            prev = v
        if state["capped"]:
            break
    if state["capped"] and tag is not None:
        CAPPED.append(tag)
    return state["best"], (bx if state["bx"] is None else state["bx"])

# ------------------------------------------------- the identity twins, sec. 4
def _in_bounds(cls, p):
    """Is this parameter vector inside the box the fitter is allowed to search?
    A twin the optimiser could not legally reach is NOT a fair replacement."""
    kk, ll, dv, npv, sg = p
    x = [np.log(kk), np.log(ll), dv]
    if cls != "independent":
        x.append(npv)
    if cls == "side-by-side":
        x.append(np.log(max(sg, 1e-300)))
    lo, hi = _bounds(cls)
    return bool(np.all(np.asarray(x) >= lo - 1e-12)
                and np.all(np.asarray(x) <= hi + 1e-12))

def twin_params(config, truth, cand, kap, lam, d, npv, sig, n, c_set):
    """The parameter vector at which `cand` reproduces `truth` EXACTLY, from the
    identities of section 4 - or None where no identity predicts one.  This is
    an algebraic prediction; the residual it produces is measured, not asserted."""
    if cand == truth:                                   # its own truth: exact
        return (kap, lam, d, npv, sig)
    if truth == "independent":                          # degeneracy 1: n' = 0
        return (kap, lam, d, 0.0, 1.0)
    if config == "constC" and (truth != "side-by-side" or d == 1.0):
        # degeneracy 3: every driver is a CONSTANT, so only lam' = lam f^n' shows
        f = {"parallel": c_set, "series": 1.0 - c_set,
             "independent": 1.0}.get(truth)
        if f is None:       # side-by-side at d = 1: C_P settles at a fixed value
            f = _mixed_p(sig * lam * (1.0 - c_set) / (kap * _pow(c_set, n)), npv)
        lam_eff = lam * _pow(f, npv)
        if cand == "parallel":
            return (kap, lam_eff / c_set, d, 1.0, 0.0)
        if cand == "series":
            return (kap, lam_eff / (1.0 - c_set), d, 1.0, 0.0)
        if cand == "independent":
            return (kap, lam_eff, d, 0.0, 0.0)
        return (kap, lam_eff, d, 0.0, 1.0)              # s-b-s with n' = 0
    if config == "mixed" and truth == "series" and cand == "parallel":
        return (kap, lam * _pow(kap, npv), d + npv, n * npv, 0.0)   # degeneracy 2
    if config == "mixed" and truth == "parallel" and cand == "series":
        ds = d - npv / n                                # the twin needs d_s >= 0
        return None if ds < 0.0 else (kap, lam * _pow(kap, -npv / n), ds,
                                      npv / n, 0.0)
    if truth == "parallel" and cand == "side-by-side":
        if d == 1.0 and npv == n:                       # degeneracy 4: p == c
            return (kap, lam, d, npv, kap / lam)
        if config in ("batch", "plug") and d == 1.0 and n == 1.0:
            return (kap, lam, 1.0, 1.0, kap * npv / lam)          # degeneracy 5
        return None
    if truth == "side-by-side" and cand == "parallel":
        if d == 1.0 and npv == n and sig == kap / lam:  # degeneracy 4
            return (kap, lam, d, npv, 0.0)
        if config in ("batch", "plug") and d == 1.0 and n == 1.0 and npv == 1.0:
            return (kap, lam, 1.0, lam * sig / kap, 0.0)          # degeneracy 5
        return None
    return None

def twin_residual(config, truth, cand, base, n, theta=TH):
    """RMS residual of the identity twin, or (None, None) if none is predicted."""
    p = twin_params(config, truth, cand, base["kappa"], base["lam"], base["d"],
                    base["nprime"], base["sig"], n, base["c_set"])
    if p is None or not _in_bounds(cand, p):
        return None, p
    y = respond(config, truth, base["kappa"], base["lam"], base["d"],
                base["nprime"], base["sig"], n=n, c_set=base["c_set"],
                theta=theta)
    yt = respond(config, cand, p[0], p[1], p[2], p[3], p[4], n=n,
                 c_set=base["c_set"], theta=theta)
    return float(np.sqrt(np.mean((yt - y) ** 2))), p

def sep_matrix(config, base=BASE, n=None, **kw):
    """D[i, j] = min(multistart fit, identity twin).  FIT records what the
    optimiser alone found; TWIN records the identity, where section 4 predicts
    one; BX records the parameters the optimiser recovered."""
    n = base["n"] if n is None else n
    D = np.zeros((4, 4)); FIT = np.zeros((4, 4)); TWIN = np.full((4, 4), np.nan)
    BX = {}
    for i, truth in enumerate(CLASSES):
        y = respond(config, truth, base["kappa"], base["lam"], base["d"],
                    base["nprime"], base["sig"], n=n, c_set=base["c_set"])
        for j, cand in enumerate(CLASSES):
            v, bx = fit_class(config, y, cand, base["kappa"], base["lam"], n=n,
                              c_set=base["c_set"], tag=(config, truth, cand),
                              **kw)
            FIT[i, j] = v; BX[(truth, cand)] = _unpack(cand, bx)
            tw, _ = twin_residual(config, truth, cand, base, n)
            if tw is not None:
                TWIN[i, j] = tw
            D[i, j] = v if tw is None else min(v, tw)
    return D, FIT, TWIN, BX

# 64 fits: 4 devices x 4 truths x 4 candidates.  No timing is printed - two
# executions of this notebook must produce byte-identical content.
SEP, SEP_FIT, SEP_TWIN, SEP_BX = {}, {}, {}, {}
for cfg in ("batch", "plug", "mixed", "constC"):
    SEP[cfg], SEP_FIT[cfg], SEP_TWIN[cfg], SEP_BX[cfg] = sep_matrix(cfg)
print(f"{4 * len(CLASSES) ** 2} fits done "
      f"({len(_seeds('side-by-side', 1.5, 3.0))} seeds for the 5-parameter "
      f"side-by-side candidate, {len(_seeds('parallel', 1.5, 3.0))} for the "
      f"others, no randomness anywhere)")
for cfg, D in SEP.items():
    lab = np.where(np.isfinite(SEP_TWIN[cfg]), "*", " ")
    df = pd.DataFrame([[f"{D[i, j]:.2e}{lab[i, j]}" for j in range(4)]
                       for i in range(4)],
                      index=[f"truth: {c}" for c in CLASSES],
                      columns=[f"fit {SHORT[c]}" for c in CLASSES])
    display(Markdown(f"**{DEV[cfg]}** - min(multistart fit, identity twin); "
                     f"`*` marks an entry whose minimum is an EXACT identity of "
                     f"section 4, not a measured separation (units: "
                     + ("ln tau'" if cfg == "constC" else "fractional conversion")
                     + ")\n\n" + md_table(df, index=True)))
    for i, tr in enumerate(CLASSES):
        for j, cd in enumerate(CLASSES):
            M[f"sep_{cfg}_{SHORT[tr]}_vs_{SHORT[cd]}"] = float(D[i, j])

# --- how much of the matrix is identity, and how far the optimiser alone was --
tw_all = np.concatenate([SEP_TWIN[c].ravel() for c in SEP])
ft_all = np.concatenate([SEP_FIT[c].ravel() for c in SEP])
has = np.isfinite(tw_all)
over = ft_all[has] - tw_all[has]
M["sep_identity_backed"] = float(has.sum())
M["sep_identity_twin_worst"] = float(np.max(tw_all[has]))
M["sep_optimiser_overstated"] = float(int((over > 1e-6).sum()))
M["sep_optimiser_worst_overstatement"] = float(np.max(over))
ov_rows = [[DEV[cfg], f"truth {tr}", f"fit {SHORT[cd]}",
            f"{SEP_FIT[cfg][i, j]:.2e}", f"{SEP_TWIN[cfg][i, j]:.2e}"]
           for cfg in SEP for i, tr in enumerate(CLASSES)
           for j, cd in enumerate(CLASSES)
           if np.isfinite(SEP_TWIN[cfg][i, j])
           and SEP_FIT[cfg][i, j] - SEP_TWIN[cfg][i, j] > 1e-6]
display(Markdown(
    f"**{int(M['sep_identity_backed'])} of 64 entries have an exact identity "
    f"behind them** (worst twin residual {M['sep_identity_twin_worst']:.1e}, "
    f"which is the ODE tolerance, not a separation). For "
    f"{int(M['sep_optimiser_overstated'])} of them the multistart alone stopped "
    f"more than 1e-6 short - these are the entries a fit-only matrix would have "
    f"printed as measured separations:\n\n"
    + md_table(pd.DataFrame(ov_rows, columns=["device", "truth", "candidate",
                                              "multistart fit alone",
                                              "identity twin"]))))
print(f"Every one of those errors is in the SAFE direction (an optimiser bound "
      f"is an upper bound on a minimum), and none of them crosses a threshold "
      f"below - the largest, {M['sep_optimiser_worst_overstatement']:.2e}, is "
      f"under the {RES_X:g} conversion resolution. But printed as fits they "
      f"would have read as measurements of a separation that is exactly zero.")
assert M["sep_optimiser_worst_overstatement"] < RES_X, \
    "an optimiser shortfall now crosses the instrument resolution - a verdict "\
    "would change, and the sentence above would be false"

# ------------------------------------------------------- the optimiser's stall
# The diagonal is a class fitted to its OWN response with every parameter free,
# so its exact minimum is identically zero: any non-zero diagonal is an
# optimiser stall, not a null and not a result.  The twin (the truth itself)
# sets the diagonal of the matrix above to zero; what is reported here is how
# far the multistart ALONE got, and where it stopped.
for cfg in SEP:
    M[f"sep_stall_{cfg}"] = float(np.max(np.diag(SEP_FIT[cfg])))
M["sep_stall_worst"] = float(max(M[f"sep_stall_{cfg}"] for cfg in SEP))
M["sep_stall_worst_not_sbs"] = float(max(
    SEP_FIT[cfg][i, i] for cfg in SEP for i, c in enumerate(CLASSES)
    if c != "side-by-side"))
stall_df = pd.DataFrame({DEV[cfg]: {c: f"{SEP_FIT[cfg][i, i]:.1e}"
                                    for i, c in enumerate(CLASSES)} for cfg in SEP})
display(Markdown("**The multistart's stall on the diagonal** - the same "
                 "machinery asked to recover the truth, whose exact minimum is "
                 "0 by construction:\n\n" + md_table(stall_df, index=True)))
print(f"worst stall over all devices and classes: {M['sep_stall_worst']:.2e}; "
      f"excluding the side-by-side candidate: {M['sep_stall_worst_not_sbs']:.2e} "
      f"- so the stall is a side-by-side phenomenon.")
assert M["sep_stall_worst"] < 1e-4, "the fitter cannot recover the truth at all"

# WHERE it stalled: the recovered parameters, which are the point of the claim.
TRU = np.array([BASE["kappa"], BASE["lam"], BASE["d"], BASE["nprime"], BASE["sig"]])
PNAMES = ("kappa", "lambda", "d", "n'", "sigma")
rec_rows, cost_rows = [], []
for cfg in ("batch", "mixed"):          # the two devices whose stall is > 1e-6
    rec = np.array(SEP_BX[cfg][("side-by-side", "side-by-side")])
    dev = np.abs(rec / TRU - 1.0)
    M[f"sep_stall_pardev_{cfg}"] = float(np.max(dev))
    rec_rows.append([DEV[cfg], f"{SEP_FIT[cfg][2, 2]:.2e}"]
                    + [f"{v:.4f}" for v in rec] + [f"{100*np.max(dev):.1f} %"])
    y0 = respond(cfg, "side-by-side", *TRU, n=BASE["n"], c_set=BASE["c_set"])
    costs = []
    for q in range(5):
        p = TRU.copy(); p[q] *= 1.01
        y1 = respond(cfg, "side-by-side", *p, n=BASE["n"], c_set=BASE["c_set"])
        costs.append(float(np.sqrt(np.mean((y1 - y0) ** 2))))
    M[f"sep_stall_1pct_cheapest_{cfg}"] = float(min(costs))
    M[f"sep_stall_1pct_cheapest_over_stall_{cfg}"] = float(
        min(costs) / SEP_FIT[cfg][2, 2])
    M[f"sep_stall_1pct_dearest_over_stall_{cfg}"] = float(
        max(costs) / SEP_FIT[cfg][2, 2])
    cost_rows.append([DEV[cfg]] + [f"{c:.2e}" for c in costs]
                     + [f"{min(costs)/SEP_FIT[cfg][2, 2]:.0f}x"])
display(Markdown("**Where the multistart stalled** (side-by-side fitted to its "
                 "own response; truth "
                 + ", ".join(f"{p} = {v:g}" for p, v in zip(PNAMES, TRU))
                 + "):\n\n"
                 + md_table(pd.DataFrame(rec_rows, columns=["device", "stall"]
                                         + list(PNAMES) + ["worst |rel. dev.|"]))
                 + "\n\n**and what a move along a single coordinate costs** "
                   "(+1 % in one parameter, RMS change in the response):\n\n"
                 + md_table(pd.DataFrame(cost_rows, columns=["device"]
                                         + [f"+1 % {p}" for p in PNAMES]
                                         + ["cheapest / stall"]))))
_rat = [M[f"sep_stall_1pct_{w}_over_stall_{c}"] for w in ("cheapest", "dearest")
        for c in ("batch", "mixed")]
print(f"Read the two tables together. A 1 % move along ANY single coordinate "
      f"costs {min(_rat):.0f}x to {max(_rat):.0f}x the stall, yet the stalled "
      f"point sits {100*M['sep_stall_pardev_batch']:.1f} % (batch) and "
      f"{100*M['sep_stall_pardev_mixed']:.1f} % (mixed flow) from the truth in "
      f"its worst coordinate. That is what 'the five constants of eq. (6) trade "
      f"against each other' means quantitatively: there is a valley along which "
      f"the response is far less sensitive than along any parameter axis, and "
      f"the trust-region solve stops in it. It is a stall - the exact minimum "
      f"is zero and the identity twin reaches it - but the reason it stalls is "
      f"a real non-identifiability of eq. (6)'s constants, and that is now a "
      f"printed result rather than an assertion.")
M["fits_capped"] = float(len(CAPPED))
print(f"\nfits whose evaluation budget bound: {len(CAPPED)} of 64"
      + (f" -> {CAPPED}" if CAPPED else ""))
'''))

cells.append(code(r'''"""Read the matrices: which pairs are separable, and at what resolution."""
res_of = {"constC": RES_TAU}
rows = []
for cfg, D in SEP.items():
    # The threshold is the INSTRUMENT RESOLUTION, and nothing else.  A floor of
    # 20 x the fitter's own stall would be the natural extra safeguard; the
    # check printed below shows it never binds on any device, so it is not
    # applied rather than carried as a decorative constraint.
    thr = res_of.get(cfg, RES_X)
    M[f"sep_threshold_{cfg}"] = float(thr)
    off = [(CLASSES[i], CLASSES[j], D[i, j]) for i in range(4) for j in range(4)
           if i != j]
    sep = [o for o in off if o[2] > thr]
    for o in sep:
        assert (cfg, o[0], o[1]) not in CAPPED, \
            f"separability claimed on a capped fit: {(cfg, o[0], o[1])}"
    best = max(off, key=lambda o: o[2])
    rows.append([DEV[cfg], f"{thr:.1e}", f"{len(sep)}/12", f"{best[2]:.2e}",
                 f"{best[0]} vs {best[1]}", f"{min(o[2] for o in off):.1e}"])
    M[f"sep_{cfg}_best_pair"] = float(best[2])
    M[f"sep_{cfg}_n_separable"] = float(len(sep))
display(Markdown(
    "**Which devices can separate anything at all?**\n\n"
    + md_table(pd.DataFrame(rows, columns=[
        "device", "threshold used", "pairs above it",
        "hardest-to-confuse pair (RMS)", "which pair", "easiest confusion"]))))
print(f"threshold = the instrument resolution: conversion {RES_X:g}, "
      f"ln tau' {RES_TAU:g}. Raw residuals are in the matrices above - rescale "
      f"at will.")

# The removed floor, priced: it never bound, and dropping it changes no verdict.
flr = pd.DataFrame([[DEV[cfg], f"{M[f'sep_stall_{cfg}']:.1e}",
                     f"{20 * M[f'sep_stall_{cfg}']:.1e}",
                     f"{res_of.get(cfg, RES_X):g}",
                     "no" if 20 * M[f"sep_stall_{cfg}"] < res_of.get(cfg, RES_X)
                     else "YES"] for cfg in SEP],
                   columns=["device", "multistart stall", "20 x stall",
                            "instrument resolution", "would the floor bind?"])
display(Markdown("**Does the fitter's own floor ever bind?**\n\n"
                 + md_table(flr)))
M["sep_floor_binds_anywhere"] = float(any(
    20 * M[f"sep_stall_{cfg}"] >= res_of.get(cfg, RES_X) for cfg in SEP))
M["sep_floor_worst_ratio"] = float(max(
    20 * M[f"sep_stall_{cfg}"] / res_of.get(cfg, RES_X) for cfg in SEP))
assert M["sep_floor_binds_anywhere"] == 0.0
print(f"Nowhere: the worst case is {M['sep_floor_worst_ratio']:.2f} of its "
      f"device's resolution, so every threshold above equals the instrument "
      f"resolution exactly and every count above is unchanged by dropping it. "
      f"A criterion that never binds is not a safeguard, and the page does not "
      f"present it as one.")
'''))

cells.append(md(r"""### 6. Why the bed is different: the profile

The plug-flow bed is the only one of the four contactings in which $C_{\rm A}$
and $C_{\rm R}$ vary *along* the catalyst, so it is the only one whose activity
is a **field** rather than a number. And the two classes lay down opposite
profiles: under **parallel** decay the inlet — richest in A — dies first; under
**series** decay the inlet stays fresh, because no R has been made yet, and
the damage starts at the outlet. That is `B2.2`'s "descending vs ascending"
signature reappearing one level of abstraction up, and it is what a well-mixed
fluid throws away.

**Both inlet values are exact, and neither needs the grid.** At $s = 0$ the
feed is imposed, so $c(0,\theta) = 1$ and $C_{\rm R}(0,\theta) = 0$ at every
time. Hence

- **parallel:** $\partial a/\partial\theta = -\lambda\,c^{n'}a^{d} = -\lambda a$
  at the feed face, so $a(0,\theta) = e^{-\lambda\theta}$ exactly;
- **series:** $\partial a/\partial\theta = -\lambda(1-c)^{n'}a^{d} = 0$ at the
  feed face for any $n' > 0$, so $a(0,\theta) \equiv 1$ exactly, at every time.

A finite-volume bed cannot print either number: its first unknown sits at the
cell *centre* $s = 1/2m$, and that value approaches the face value at **first
order** in $1/m$. The refinement is measured below, and the exact values are
what the page reports.

**The outlet has no closed form, so it is extrapolated — twice.** Once in
space, from the two finest grids (the cell-centre offset is first order, ratio
2, so the correction is the last change itself); and then once in **time**,
because the space extrapolation is performed at a fixed step count and Heun's
temporal error there turns out to be *larger* than what the space
extrapolation leaves behind. Both refinements are printed below, with their
observed orders, and the reported outlet is quoted only to the precision the
two of them jointly resolve. This is the third time in this one section that a
number has had to stop being a sample of whatever resolution the code was run
at; the closing note says what the general rule is.

**This, and not the conversion history, is where the bed's information
actually is.** The outlet histories of the two classes are 8.1e-4 apart in RMS
conversion after a full run — a tenth of a percentage point, and a hard
measurement. The inlet activities are more than an order of magnitude apart and
the *sign* of the gradient differs. Any assay of the spent catalyst — coke
burn-off, a crushed-and-retested activity, anything that resolves 10 % — reads
that inlet contrast straight off.

**One caveat carried from `B2.2`, which this page cites as agreement.** `B2.2`
measured its own *within-bed* coke contrast at its operating point and found
the descending (parallel) signature "already below a 10 %-assay detection floor
and fading as $1/t$". That is not a contradiction of what is claimed here,
because the two are different quantities: `B2.2`'s discriminator is the
*gradient* of the deposit along the bed, this page's is the **inlet activity
itself**, at one sampling point. But a
reader taking the "one assay is enough" advice to a real bed should size the
contrast they can actually reach, `B2.2`-style, before designing on it. The
recommendation this page ends on is therefore not about which reactor to run
but about **what to measure in it**."""))

cells.append(code(r'''"""Figure 2: activity profiles down the pymrm bed, parallel vs series."""
prof = {}
for cls in ("parallel", "series"):
    Xb, Ab, tsb = BED.march(cls, BASE["kappa"], BASE["lam"], BASE["d"],
                            BASE["nprime"], n=BASE["n"], nt=200, theta=TH)
    prof[cls] = (Xb, Ab)
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.2))
NT_PROF = 200                       # the march grid the profiles are recorded on
for ax, cls in zip(axes, ("parallel", "series")):
    for tt in (0.0, 0.25, 0.5, 1.0):
        j = int(round(tt * NT_PROF))
        ax.plot(BED.s_c, prof[cls][1][:, j], color=COL[cls],
                alpha=0.35 + 0.65 * tt, lw=1.6,
                label=rf"$\theta$ = {j / NT_PROF:.2f}")
    ax.set_title(f"{cls} deactivation (eq. {4 if cls=='parallel' else 5})",
                 fontsize=9)
    ax.set_xlabel("$s = w/W$  (fraction of the catalyst passed)")
    ax.set_ylabel("activity $a$"); ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.25, lw=0.5); ax.legend(fontsize=7.5, frameon=False)
fig.tight_layout(); plt.show()

# --- the inlet is EXACT, and is not read off the grid -----------------------
# At s = 0 the feed is imposed: c = 1, so C_R = 0.  Hence
#   parallel: da/dtheta = -lam c^n' a^d = -lam a  ->  a(0, theta) = exp(-lam t)
#   series:   da/dtheta = -lam (1-c)^n' a^d = 0   ->  a(0, theta) = 1, exactly
# The first finite-volume unknown sits at s = 1/2m, not at s = 0, and converges
# to these at FIRST order.  The exact values are what is reported; the grid
# study below is the evidence, and the outlet (no closed form) is extrapolated
# TWICE - once in space and once in time (see the second table).
M["bed_a_inlet_end_par"] = float(np.exp(-BASE["lam"] * TH[-1]))
M["bed_a_inlet_end_ser"] = 1.0
GRIDS = (100, 200, 400, 800)
NT_CONV = 100          # base step count of the refinement study.  NOTHING
                       # reported below is left a function of it: the outlet is
                       # extrapolated in nt as well, and the break table halves
                       # NT_CONV to show what that is worth.
INL, OUT = {}, {}      # (class, m, nt) -> end cell-centre activities, as FLOATS
def bed_ends(cls, m, nt):
    """a at the first and last cell centre at theta = 1.  Cached, and returned
    as floats - never as formatted strings, which quantise the very differences
    the refinement study exists to measure."""
    if (cls, m, nt) not in OUT:
        _, Ag, _ = Bed(m).march(cls, BASE["kappa"], BASE["lam"], BASE["d"],
                                BASE["nprime"], n=BASE["n"], nt=nt, theta=TH)
        INL[(cls, m, nt)] = float(Ag[0, -1])
        OUT[(cls, m, nt)] = float(Ag[-1, -1])
    return INL[(cls, m, nt)], OUT[(cls, m, nt)]

gr = []
for cls in ("parallel", "series"):
    ex = M[f"bed_a_inlet_end_{SHORT[cls]}"]
    errs, outs = [], []
    for m in GRIDS:
        ain, aout = bed_ends(cls, m, NT_CONV)
        errs.append(abs(ain - ex)); outs.append(aout)
        gr.append([cls, m, f"{ain:.6f}", f"{errs[-1]:.3e}",
                   "" if len(errs) == 1 else f"{np.log2(errs[-2]/errs[-1]):.2f}",
                   f"{outs[-1]:.6f}",
                   "" if len(outs) == 1 else f"{outs[-1]-outs[-2]:+.2e}"])
    M[f"bed_a_inlet_cellcentre_err_{SHORT[cls]}"] = float(errs[-1])
    M[f"bed_a_inlet_cellcentre_order_{SHORT[cls]}"] = float(
        np.log2(errs[-2] / errs[-1]))
display(Markdown(
    "**The bed's first and last cell centres, refined in SPACE** (the inlet "
    "column is measured against the EXACT face value; the outlet, which has no "
    "closed form, is Richardson-extrapolated from the last pair - and then "
    "extrapolated again, in time, in the table below):\n\n"
    + md_table(pd.DataFrame(gr, columns=[
        "class", "cells m", "a at s = 1/2m", "|error| vs the exact inlet",
        "observed order", "a at s = 1-1/2m", "change in it"]))))

# --- the INLET study is spatial; the OUTLET reading is not -------------------
# Both are checked on FLOATS, on BOTH classes, at the finest grid of the study.
MF = GRIDS[-1]
M["bed_inlet_time_sensitivity"] = float(max(
    abs(bed_ends(c, MF, 2 * NT_CONV)[0] - bed_ends(c, MF, NT_CONV)[0])
    for c in ("parallel", "series")))
M["bed_profile_time_sensitivity"] = float(max(
    abs(bed_ends(c, MF, 2 * NT_CONV)[1] - bed_ends(c, MF, NT_CONV)[1])
    for c in ("parallel", "series")))
print(f"observed order of the cell-centre error: "
      f"{M['bed_a_inlet_cellcentre_order_par']:.2f} (parallel), "
      f"{M['bed_a_inlet_cellcentre_order_ser']:.2f} (series) - first order, "
      f"which is what a 1/2m offset from the face demands.")
print(f"The INLET study is purely spatial: doubling nt at m = {MF} moves "
      f"either inlet reading by at most {M['bed_inlet_time_sensitivity']:.1e}, "
      f"{M['bed_a_inlet_cellcentre_err_par']/M['bed_inlet_time_sensitivity']:.0f}x "
      f"below the parallel cell-centre error it is set against and "
      f"{M['bed_a_inlet_cellcentre_err_ser']/M['bed_inlet_time_sensitivity']:.0f}x "
      f"below the series one.")
print(f"The OUTLET reading is NOT: the same doubling moves it by "
      f"{M['bed_profile_time_sensitivity']:.1e}, "
      f"{M['bed_profile_time_sensitivity']/M['bed_inlet_time_sensitivity']:.0f}x "
      f"more, and larger than what the space extrapolation leaves behind. So "
      f"the outlet is extrapolated a SECOND time, in nt.")

# --- second extrapolation: in TIME ------------------------------------------
# Richardson in space (first order, ratio 2) at each of three step counts, then
# Richardson in time (Heun is second order, ratio 2, so the correction is the
# last change over 3).  Reported to the precision the two extrapolations
# actually support, which is five decimals, not six.
def space_rich(cls, nt, pair=None):
    o = [bed_ends(cls, m, nt)[1] for m in (GRIDS[-2:] if pair is None else pair)]
    return o[-1] + (o[-1] - o[-2])

NT_SET = (NT_CONV // 2, NT_CONV, 2 * NT_CONV)
RR, tr = {}, []
for cls in ("parallel", "series"):
    R = [space_rich(cls, nt) for nt in NT_SET]
    RR[cls] = R
    dR = (R[1] - R[0], R[2] - R[1])
    M[f"bed_a_outlet_time_order_{SHORT[cls]}"] = float(np.log2(dR[0] / dR[1]))
    M[f"bed_a_outlet_spaceonly_{SHORT[cls]}"] = float(R[1])
    M[f"bed_a_outlet_end_{SHORT[cls]}"] = float(R[2] + dR[1] / 3.0)
    M[f"bed_a_outlet_altpair_{SHORT[cls]}"] = float(R[1] + dR[0] / 3.0)
    M[f"bed_a_outlet_time_bias_{SHORT[cls]}"] = float(
        M[f"bed_a_outlet_spaceonly_{SHORT[cls]}"]
        - M[f"bed_a_outlet_end_{SHORT[cls]}"])
    # how far the SPACE extrapolation itself can still be off: the two grid
    # pairs available disagree by this much, at the base step count
    M[f"bed_a_outlet_space_spread_{SHORT[cls]}"] = float(abs(
        space_rich(cls, NT_CONV, GRIDS[-3:-1]) - space_rich(cls, NT_CONV)))
    M[f"bed_a_outlet_residual_bound_{SHORT[cls]}"] = float(max(
        abs(dR[1] / 3.0), M[f"bed_a_outlet_space_spread_{SHORT[cls]}"]))
    M[f"bed_profile_contrast_{SHORT[cls]}"] = float(
        M[f"bed_a_outlet_end_{SHORT[cls]}"] - M[f"bed_a_inlet_end_{SHORT[cls]}"])
    for j, (nt, r) in enumerate(zip(NT_SET, R)):
        tr.append([cls, nt, f"{r:.7f}", "" if j == 0 else f"{r - R[j-1]:+.2e}",
                   "" if j < 2 else f"{M[f'bed_a_outlet_time_order_{SHORT[cls]}']:.2f}"])
    tr.append([cls, "extrapolated", f"{M[f'bed_a_outlet_end_{SHORT[cls]}']:.7f}",
               f"{dR[1]/3.0:+.2e}", ""])
display(Markdown(
    "**The space-extrapolated outlet, refined in TIME** (each row repeats the "
    "Richardson of the two finest grids above at that step count; Heun is "
    "second order, so the change falls by 4 per doubling and the remaining "
    "correction is that change over 3):\n\n"
    + md_table(pd.DataFrame(tr, columns=[
        "class", "steps nt", "space-Richardson outlet", "change in it",
        "observed order in nt"]))))
print(f"observed temporal order of the space-extrapolated outlet: "
      f"{M['bed_a_outlet_time_order_par']:.2f} (parallel), "
      f"{M['bed_a_outlet_time_order_ser']:.2f} (series) - second order, which "
      f"is what Heun demands and what the validation section measures "
      f"independently on the conversion history.")
print(f"WHAT THE SECOND EXTRAPOLATION IS WORTH. The space-Richardson ALONE, at "
      f"nt = {NT_CONV}, gives {M['bed_a_outlet_spaceonly_par']:.6f} for the "
      f"parallel outlet - {abs(M['bed_a_outlet_time_bias_par']):.1e} above the "
      f"limit, which is larger than what the space extrapolation itself leaves "
      f"({M['bed_a_outlet_space_spread_par']:.1e}, the disagreement between the "
      f"two grid pairs). A reading that is space-converged can still be "
      f"time-limited, and this one was.")
print(f"\nat theta = 1:  parallel  a(inlet) = {M['bed_a_inlet_end_par']:.6f} "
      f"(EXACT, = exp(-lambda)), a(outlet) = {M['bed_a_outlet_end_par']:.5f} "
      f"(space + time extrapolated, residual bound "
      f"{M['bed_a_outlet_residual_bound_par']:.1e})   -> descending")
print(f"               series    a(inlet) = {M['bed_a_inlet_end_ser']:.6f} "
      f"(EXACT, C_R = 0 at the feed face), a(outlet) = "
      f"{M['bed_a_outlet_end_ser']:.5f} (space + time extrapolated, residual "
      f"bound {M['bed_a_outlet_residual_bound_ser']:.1e})   -> ascending")
print(f"The outlets are quoted to five decimals and not six: the residual "
      f"bounds above are {max(M['bed_a_outlet_residual_bound_par'], M['bed_a_outlet_residual_bound_ser']):.1e} "
      f"at worst, which resolves the fifth decimal and not the sixth. Six was "
      f"printed here before, and the last two of them were wrong.")
# The defect this repairs becomes two rows per class - one that must move and
# one that must not.  Both classes, because both outlets are reported.
for cls in ("parallel", "series"):
    sh = SHORT[cls]
    BREAKS.append((f"bed_a_outlet_spaceonly_{sh}", M[f"bed_a_outlet_spaceonly_{sh}"],
                   f"halve NT_CONV to {NT_SET[0]} and report the "
                   f"space-Richardson alone, as the outlet was reported before",
                   float(RR[cls][0]),
                   "the space-only reading IS a function of the step count, "
                   "and nothing else on this page constrains that choice - "
                   "which is exactly why it is not what is reported"))
    BREAKS.append((f"bed_a_outlet_end_{sh}", M[f"bed_a_outlet_end_{sh}"],
                   f"read the time extrapolation off the nt = {NT_SET[0]}/"
                   f"{NT_SET[1]} pair instead of {NT_SET[1]}/{NT_SET[2]}",
                   M[f"bed_a_outlet_altpair_{sh}"],
                   "ROBUSTNESS (must NOT move): the reported outlet is the "
                   "limit of the march, so which pair the limit is read from "
                   "cannot matter - and after the second extrapolation it "
                   "does not"))
print(f"the two contrasts have opposite sign: "
      f"{M['bed_profile_contrast_par']:+.4f} vs {M['bed_profile_contrast_ser']:+.4f}")
print(f"For comparison, the production grid (m = {BED.m}) reads its first cell "
      f"centre at {prof['series'][1][0, -1]:.4f} (series) and "
      f"{prof['parallel'][1][0, -1]:.6f} (parallel). Neither is reported as the "
      f"inlet activity: a sampled cell centre standing in for an exactly known "
      f"limit is the wrong-baseline defect class this gallery keeps finding "
      f"(`B1.7`), and here the exact statement is the stronger one - the series "
      f"inlet does not merely stay high, it does not decay AT ALL.")
print(f"\nTHE SAME MISTAKE, THREE TIMES, ON THIS ONE SECTION. The inlet was "
      f"first a sampled cell centre where an exact limit exists; the outlet "
      f"was then a single-grid sample where a space extrapolation exists; and "
      f"the outlet was then a space extrapolation at one step count where a "
      f"TIME extrapolation exists. All three are one defect: a number read at "
      f"whatever resolution the code happened to be run at, and printed as if "
      f"it were the limit. The general rule this section now follows is that a "
      f"reported value must not be a function of any resolution the page "
      f"chose - every such knob is either extrapolated away or shown, by "
      f"halving it in the break table, not to matter - and that a value may "
      f"not be printed to more digits than its own refinement resolves.")
'''))

cells.append(md(r"""### 7. Where even the mixed reactor separates — and the boundary is analytic

Degeneracy 2 is **one-way**. Series with orders $(n', d)$ always has an exact
parallel twin at $(n n', d+n')$. Going the other way, a parallel truth
$(n'_p, d_p)$ has a series twin only if the required series order
$d_s = d_p - n'_p/n$ is **non-negative** — and $d$ is an *order* in the
paper's own rate forms, so a negative one is not a member of the family.

So there is a line, $d = n'/n$, below which the exact twin ceases to exist.
The prediction is sharp; what it buys is not. Below the line the best series
impersonation is no longer perfect, but it is not immediately *bad* either —
the residual grows smoothly with distance from the line, and only becomes
measurable at a realistic resolution well inside the region. Tested on both
sides:"""))

cells.append(code(r'''"""The predicted boundary d = n'/n in the constant-flow mixed reactor."""
rows = []
for npv, dv in ((1.0, 1.5), (1.0, 1.0), (1.0, 0.75), (1.0, 0.5), (1.0, 0.25),
                (2.0, 1.5)):
    y = respond("mixed", "parallel", BASE["kappa"], BASE["lam"], dv, npv, 0.0)
    D = fit_class("mixed", y, "series", BASE["kappa"], BASE["lam"],
                  tag=("mixed-boundary", npv, dv))[0]
    ds = dv - npv / BASE["n"]
    # where the twin exists it is evaluated directly, exactly as in section 5:
    # an optimiser bound must never stand in for an identity
    tw, _ = twin_residual("mixed", "parallel", "series",
                          dict(BASE, d=dv, nprime=npv, sig=0.0), BASE["n"])
    if tw is not None:
        D = min(D, tw)
    rows.append([npv, dv, f"{ds:+.2f}", "attainable" if ds >= 0 else "d < 0: none",
                 f"{D:.2e}{'*' if tw is not None else ''}",
                 "no" if D <= RES_X else "YES"])
    M[f"mixed_boundary_np{npv:g}_d{dv:g}"] = float(D)
display(Markdown(
    "**Mixed flow, constant rate: can a series model impersonate a parallel "
    "truth?**\n\n" + md_table(pd.DataFrame(rows, columns=[
        "n'", "d", "d - n'/n", "series twin", "D(parallel -> series)",
        f"separable at {RES_X:g}?"]))))
inside = [r for r in rows if float(r[2]) >= 0]
outside = [r for r in rows if float(r[2]) < 0]
M["mixed_boundary_worst_inside"] = float(max(M[f"mixed_boundary_np{r[0]:g}_d{r[1]:g}"]
                                             for r in inside))
M["mixed_boundary_best_outside"] = float(min(M[f"mixed_boundary_np{r[0]:g}_d{r[1]:g}"]
                                             for r in outside))
print(f"\nInside the region where the exact twin exists (d >= n'/n) the best "
      f"series impersonation is <= {M['mixed_boundary_worst_inside']:.1e} - i.e. "
      f"nothing. Outside it the residual grows smoothly with distance from the "
      f"line, from {M['mixed_boundary_np1_d0.75']:.1e} at d - n'/n = -0.25 to "
      f"{M['mixed_boundary_np1_d0.25']:.1e} at -0.75, and only crosses the "
      f"{RES_X:g} resolution around d - n'/n = -0.5.")
print("The line is a prediction about EXACTNESS, and it holds. It is not a "
      "prediction about measurability, and the (n' = 2, d = 1.5) row shows why: "
      f"it lies half a unit below the line yet is still fitted to "
      f"{M['mixed_boundary_np2_d1.5']:.1e}, because an inexact series twin is "
      "still an excellent one there. Crossing the boundary is necessary for the "
      "mixed reactor to separate the two classes, and nowhere near sufficient.")
BREAKS.append(("mixed_boundary_np1_d0.5", M["mixed_boundary_np1_d0.5"],
               "move d back across the predicted boundary to d = 1.5",
               M["mixed_boundary_np1_d1.5"],
               "the separability of the mixed reactor is switched by the "
               "analytic condition, in both directions"))
'''))

cells.append(md(r"""### 8. The recommended device, taken seriously

Levenspiel's recommendation (Discussion, point 2, p. 272) is unambiguous:

> *"the most useful reactor set up uses a batch of solids and mixed flow of
> fluid with changing flow rate so as to keep the concentration of the
> pertinent reaction component unchanged with time within the reactor. Thus by
> constantly lowering the flow rate in a basket or recycle reactor the activity
> and concentration dependencies can be decoupled and studied separately."*

Everything about that is right, and the page confirms both halves: the decay
order $d$ comes straight out of the shape of $\tau'(t)$ through eqs. (28)–(32)
**whatever $n'$ is** (that is the decoupling), and — from degeneracy 3 — the
class does not come out at all. The two facts are the same fact. So the
question the device leaves open is: *how many runs, over what span of held
concentration, does the class actually take?*

The recipe is the paper's own: $k'_d = k_dC_{\rm A}^{\,n'}$ under parallel and
$k_dC_{\rm R}^{\,n'}$ under series, so a set of runs at different held levels
gives $k'_d$ as a function of the level, and a power law in $C_{\rm A}$ is not
a power law in $C_{\rm A0}-C_{\rm A}$. Below, the span required is
**root-found** against a stated precision on $k'_d$, not read off a sweep — and
then swept over the two design choices it hides, because root-finding a
function of an arbitrary parameter does not make that parameter less arbitrary.

**The second half of this section prices a protocol the paper does not
prescribe, and says so.** For eq. (6) Levenspiel's own instruction is to hold
$C_{\rm P}$; holding $C_{\rm A}$ instead is the inverse of his priority. Both
sentences are in `data/levenspiel-1972-printed-claims.csv` and both are printed
below from it in full, hedge included. Under the author's own instruction — a
poison not consumed, so
$C_{\rm P}$ fixed — this page's own $\sigma = 0$ row shows **no bias at all**.
What is priced below is the cost of the *other* choice, and it is a real cost
only because one control knob cannot hold two concentrations."""))

cells.append(code(r'''"""How wide a set of held-C_A runs does it take to name the class?

At each held level c the device returns k'_d exactly (degeneracy 3 makes it
the ONLY thing it returns).  Under a parallel truth k'_d = lam c^n'; a series
reader fits ln k'_d = ln A + n'_s ln(1-c).  The relative misfit of that wrong
reading, over L levels spanning [c0, c0+dc], is the signal.  Root-found."""
L_LEVELS, C_LO = 5, 0.15
def wrong_class_misfit(dc, npv=BASE["nprime"], truth="parallel",
                       lo=None, levels=None):
    c = np.linspace(C_LO if lo is None else lo,
                    (C_LO if lo is None else lo) + dc,
                    L_LEVELS if levels is None else levels)
    kd_true = (c if truth == "parallel" else 1 - c) ** npv     # lam folded out
    x = np.log((1 - c) if truth == "parallel" else c)
    y = np.log(kd_true)
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.sqrt(np.mean((y - A @ coef) ** 2)))   # RMS in ln k'_d

for prec in (0.02, 0.05, 0.10):
    f = lambda dc: wrong_class_misfit(dc) - prec
    dc_star = brentq(f, 1e-4, 0.80, xtol=1e-10)
    M[f"constC_span_needed_prec{int(prec*100)}pct"] = float(dc_star)
    print(f"k'_d known to {prec*100:4.1f} %  ->  held-C_A levels must span "
          f"dc = {dc_star:.3f} in C_A/C_A0 "
          f"(i.e. {C_LO:.2f} to {C_LO+dc_star:.2f}), with {L_LEVELS} levels")
M["constC_misfit_at_span_0p5"] = wrong_class_misfit(0.5)
print(f"\nover a half-decade of held concentration (dc = 0.5) the wrong reading "
      f"misses by {100*M['constC_misfit_at_span_0p5']:.1f} % in k'_d - "
      f"comfortable; over dc = 0.05 it misses by "
      f"{100*wrong_class_misfit(0.05):.2f} %, which no rate constant is known to.")

# --- and the span is a property of TWO choices that are not the paper's -----
# C_LO (where the set of held levels starts) and L_LEVELS (how many there are)
# are design choices, not printed constants.  Root-finding a function of an
# arbitrary parameter does not make the parameter less arbitrary, so both are
# swept and the answer is reported as a RANGE, with the base cell marked.
LO_SCAN, LEV_SCAN = (0.05, 0.10, 0.15, 0.25, 0.40), (3, 5, 9)
span_grid, all_spans = [], []
for lo in LO_SCAN:
    row = [f"{lo:.2f}" + ("  <-- base" if lo == C_LO else "")]
    for lv in LEV_SCAN:
        v = brentq(lambda dc: wrong_class_misfit(dc, lo=lo, levels=lv) - 0.05,
                   1e-4, min(0.80, 0.98 - lo), xtol=1e-10)
        all_spans.append(v)
        row.append(f"{v:.3f}" + ("*" if (lo == C_LO and lv == L_LEVELS) else ""))
    span_grid.append(row)
M["constC_span_prec5pct_min"] = float(min(all_spans))
M["constC_span_prec5pct_max"] = float(max(all_spans))
display(Markdown(
    "**The required span is a property of two unprinted design choices** - "
    "where the set of held levels starts (C_LO) and how many there are "
    "(L_LEVELS). Span needed at 5 % precision on k'_d, root-found in every "
    "cell; `*` is the base case quoted above:\n\n"
    + md_table(pd.DataFrame(span_grid,
                            columns=["C_LO \\ levels"]
                            + [str(lv) for lv in LEV_SCAN]))))
print(f"Range over the sweep: {M['constC_span_prec5pct_min']:.3f} to "
      f"{M['constC_span_prec5pct_max']:.3f} - a factor of "
      f"{M['constC_span_prec5pct_max']/M['constC_span_prec5pct_min']:.1f}, and "
      f"monotone in the starting level. The base answer "
      f"{M['constC_span_needed_prec5pct']:.3f} is therefore an ORDER, not a "
      f"design constant: a set of held levels must span a few tenths of "
      f"C_A/C_A0, and where it starts matters as much as how wide it is.")
BREAKS.append(("constC_span_prec5pct_max", M["constC_span_prec5pct_max"],
               f"start the level set at C_A/C_A0 = {LO_SCAN[0]:g} instead of "
               f"{LO_SCAN[-1]:g}",
               M["constC_span_prec5pct_min"],
               "the required span is set by the unprinted C_LO as much as by "
               "the demanded precision"))
BREAKS.append(("constC_span_needed_prec5pct",
               M["constC_span_needed_prec5pct"],
               "demand 2 % precision on k'_d instead of 5 %",
               M["constC_span_needed_prec2pct"],
               "the required span is set by the precision, and moves with it"))
'''))

cells.append(code(r'''"""One knob, two constraints: side-by-side in the recommended device.

Holding C_A fixed sets the flow rate.  If the poison is CONSUMED (eq. 6 prints
"P -> P(down)"), C_P then drifts as the flow rate falls, so k'_d is not constant
and the eq.-(32) straight line is bent.  What the paper actually PRESCRIBES for
eq. (6) is printed below, from the claims CSV, in full."""
for cid in ("p271_keep_constant", "p271_side_by_side_hedge"):
    display(Markdown(f"> *\"{CLAIMS.loc[cid, 'quote']}\"* &nbsp;&nbsp;(p. "
                     f"{CLAIMS.loc[cid, 'page']}, `{cid}`)"))
print("Read together: for eq. (6) the author's instruction is to keep C_P "
      "constant, and holding C_A too is the conditional extra ('if possible'). "
      "He then hedges the failure of the SECOND condition, not the first: 'If "
      "C_A cannot be kept constant in side-by-side deactivation, analysis is "
      "still not particularly difficult.' What follows holds C_A and lets C_P "
      "drift - the inverse of that priority - because that is what a single "
      "flow-rate knob does when the poison is consumed. The author's own "
      "protocol is the sigma = 0 row below, and it is unbiased.")

def recovered_d(sig, d_true):
    """The order an experimenter would read off the paper's own Fig.-6 plot."""
    tau, _ = run_const_c("side-by-side", BASE["kappa"], BASE["lam"], d_true,
                         BASE["nprime"], BASE["c_set"], n=BASE["n"], sig=sig)
    def bend(dv):
        y = np.log(tau) if abs(dv - 1) < 1e-9 else (tau ** (dv - 1) - 1) / (dv - 1)
        A = np.vstack([np.ones_like(TH), TH]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        return float(np.sqrt(np.mean((y - A @ coef) ** 2)) / max(np.ptp(y), 1e-300))
    return float(minimize_scalar(bend, bounds=(0.05, 4.0), method="bounded",
                                 options={"xatol": 1e-10}).x)

# --- d = 1 is a blind spot of the device, and exactly so --------------------
# In a mixed reactor at held C_A the poison consumed per unit time goes as
# tau' x rate ~ (1/a) x a^d = a^(d-1).  At d = 1 that is CONSTANT: C_P settles
# at a fixed value and never drifts, so the eq.-(32) line stays perfectly
# straight no matter how much poison is consumed.  The over-determination has
# no cost at all at d = 1 - and only at d = 1.
SIG_SCAN = np.geomspace(1e-3, 2e1, 21)
bias1 = np.array([abs(recovered_d(sv, 1.0) - 1.0) for sv in SIG_SCAN])
M["constC_sbs_bias_at_d1_max"] = float(bias1.max())
print(f"d = 1: the recovered order is unbiased to {M['constC_sbs_bias_at_d1_max']:.1e} "
      f"for every sigma up to {SIG_SCAN[-1]:g} - an exact cancellation "
      f"(consumption rate ~ a^(d-1) = const at d = 1), not a small number.")

# --- d = 2, where the drift is real ---------------------------------------
D_TRUE = 2.0
M["constC_sbs_recovered_d_sig0"] = recovered_d(0.0, D_TRUE)
M["constC_sbs_recovered_d_base"] = recovered_d(BASE["sig"], D_TRUE)
bias2 = np.array([abs(recovered_d(sv, D_TRUE) / D_TRUE - 1.0) for sv in SIG_SCAN])
M["constC_sbs_bias_max_scanned"] = float(bias2.max())
target = next((t for t in (0.05, 0.02, 0.01) if bias2.max() > t), None)
assert target is not None, f"no bias target reached; max scanned {bias2.max():.3g}"
j = int(np.argmax(bias2 > target))
M["constC_sbs_bias_target"] = float(target)
M["constC_sbs_sigma_for_bias"] = float(brentq(
    lambda sv: abs(recovered_d(sv, D_TRUE) / D_TRUE - 1.0) - target,
    SIG_SCAN[j - 1], SIG_SCAN[j], xtol=1e-9))
M["constC_sbs_bias_sig0"] = float(abs(M["constC_sbs_recovered_d_sig0"] / D_TRUE
                                      - 1.0))
print(f"\nd = {D_TRUE:g}: recovered {M['constC_sbs_recovered_d_sig0']:.4f} at "
      f"sigma = 0 and {M['constC_sbs_recovered_d_base']:.4f} at "
      f"the base sigma = {BASE['sig']} - a bias of "
      f"{100*(M['constC_sbs_recovered_d_base']/D_TRUE - 1):+.2f} %.")
print(f"The sigma = 0 row IS the author's own protocol for eq. (6) - a poison "
      f"in such excess that C_P does not move is C_P kept constant - and it "
      f"recovers d to {M['constC_sbs_bias_sig0']:.1e} relative. The bias below "
      f"is the price of the OTHER protocol, the one a single flow-rate knob "
      f"forces on you if the poison is consumed and you choose to hold C_A.")
print(f"root-found (bracketed by the scan, not read off it): the bias reaches "
      f"{100*target:.0f} % once sigma exceeds {M['constC_sbs_sigma_for_bias']:.4f}; "
      f"the largest bias over sigma <= {SIG_SCAN[-1]:g} is "
      f"{100*M['constC_sbs_bias_max_scanned']:.1f} %.")
BREAKS.append(("constC_sbs_recovered_d_base", M["constC_sbs_recovered_d_base"],
               "set sigma = 0 (poison in such excess that C_P never moves)",
               M["constC_sbs_recovered_d_sig0"],
               "the bias IS the poison consumption; with none, d comes back exactly"))
BREAKS.append(("constC_sbs_bias_at_d1_max", M["constC_sbs_bias_at_d1_max"],
               f"move off d = 1, where the a^(d-1) cancellation holds, to d = {D_TRUE:g}",
               M["constC_sbs_bias_max_scanned"],
               "the blind spot is exactly at d = 1 and nowhere else"))
'''))

# --------------------------------------------------------------- validation -
cells.append(md(r"""## Validation

Four layers, in order of what they can catch.

1. **The pymrm bed against an independent exact route.** The pseudo-steady bed
   at frozen $a(s)$ is a first-order cascade, so the outlet is a quadrature —
   no grid, no operator, no linear solve in common with the finite-volume bed.
   Both refinement axes are swept separately.
2. **Every reactor against the paper's own closed forms.** Eqs. (15), (19),
   (22) and (28)–(31) are exact solutions of the very systems solved here, for
   the illustrative first-order case; they are the strongest available check
   and they are the paper's, not this page's.
3. **A defect-injection table**, with a row for every reported metric.
4. **The pseudo-steady assumption itself**, dropped and costed."""))

cells.append(code(r'''"""(1a) SPATIAL order, isolated: frozen a(s), bed outlet vs exact quadrature."""
from scipy.integrate import cumulative_simpson as _csimp
def exact_outlet(kappa, n, afun, ns=20001):
    s = np.linspace(0, 1, ns)
    return _cascade(np.concatenate(([0.0], _csimp(afun(s), x=s))), kappa, n)[-1]

rows = []
for label, afun in (("a = e^{-3s}", lambda s: np.exp(-3 * s)),
                    ("a = 0.2 + 0.8 s^2", lambda s: 0.2 + 0.8 * s ** 2)):
    for n_ in (1.0, 2.0):
        ref = exact_outlet(1.5, n_, afun); prev = None
        for m in (25, 50, 100, 200, 400):
            b = Bed(m)
            e = abs(b.outlet(b.steady_c(1.5, afun(b.s_c), n_)) - ref)
            rows.append([label, n_, m, f"{e:.3e}",
                         "" if prev is None else f"{np.log2(prev/e):.2f}"])
            prev = e
        M[f"bed_spatial_err_n{n_:g}_{'exp' if 'e^' in label else 'quad'}"] = float(e)
        M[f"bed_spatial_order_n{n_:g}_{'exp' if 'e^' in label else 'quad'}"] = \
            float(np.log2(float(rows[-2][3]) / e))
display(Markdown("**Spatial convergence of the pymrm bed against the quadrature "
                 "(activity frozen, so no time error can leak in):**\n\n"
                 + md_table(pd.DataFrame(rows,
                                          columns=["frozen a(s)", "n", "cells m",
                                                   "|outlet error|",
                                                   "observed order"]))))
print("compute_boundary_values with the outflow bc left as pymrm builds it: "
      "second order, at both reaction orders.")
'''))

cells.append(code(r'''"""(1b) TEMPORAL order, and the joint production setting."""
ref_par = PlugExact(161).run("parallel", BASE["kappa"], BASE["lam"], BASE["d"],
                             BASE["nprime"], n=BASE["n"], rtol=1e-11)[0]
NTS = (25, 50, 100, 200)
b400 = Bed(400)
rows, ERRT = [], {}
for scheme in ("heun", "euler"):
    prev = None; ERRT[scheme] = []
    for nt in NTS:
        X, _, _ = b400.march("parallel", BASE["kappa"], BASE["lam"], BASE["d"],
                             BASE["nprime"], n=BASE["n"], nt=nt, scheme=scheme)
        e = float(np.max(np.abs(X - ref_par))); ERRT[scheme].append(e)
        rows.append([scheme, nt, f"{e:.3e}",
                     "" if prev is None else f"{np.log2(prev/e):.2f}"])
        prev = e
M["bed_time_order"] = float(np.log2(ERRT["heun"][1] / ERRT["heun"][2]))
M["bed_time_order_euler"] = float(np.log2(ERRT["euler"][1] / ERRT["euler"][2]))
M["bed_time_order_euler_finest"] = float(np.log2(ERRT["euler"][2]
                                                 / ERRT["euler"][3]))
M["bed_time_order_heun_finest"] = float(np.log2(ERRT["heun"][2]
                                                / ERRT["heun"][3]))
display(Markdown("**Temporal convergence on the activity march (m = 400), "
                 "second-order Heun against first-order Euler:**\n\n"
                 + md_table(pd.DataFrame(rows, columns=["scheme", "steps nt",
                                                        "max |X error|",
                                                        "observed order"]))))
print(f"Heun is second order ({M['bed_time_order']:.2f} at the 50->100 pair) "
      f"and Euler is first ({M['bed_time_order_euler']:.2f} at the same pair, "
      f"and {M['bed_time_order_euler_finest']:.2f} at 100->200). "
      f"bed_time_order is read at 50->100 DELIBERATELY: at nt = 200 the Heun "
      f"time error ({ERRT['heun'][3]:.1e}) has fallen onto the m = 400 SPATIAL "
      f"floor (the frozen-a study above puts that at ~6e-7 to 2e-6), so the "
      f"apparent order there is {M['bed_time_order_heun_finest']:.2f} - an "
      f"artefact of the reference, not a super-convergence. Euler's error is "
      f"still {ERRT['euler'][3]/ERRT['heun'][3]:.0f}x above that floor at "
      f"nt = 200, which is why it keeps a clean order all the way down.")
Xprod, _, _ = BED.march("parallel", BASE["kappa"], BASE["lam"], BASE["d"],
                        BASE["nprime"], n=BASE["n"], nt=200)
M["bed_vs_quadrature_production"] = float(np.max(np.abs(Xprod - ref_par)))
M["bed_vs_quadrature_over_resolution"] = float(M["bed_vs_quadrature_production"]
                                               / RES_X)
print(f"production bed (m = {BED.m}, nt = 200) vs the quadrature: "
      f"{M['bed_vs_quadrature_production']:.2e} in conversion, "
      f"{1/M['bed_vs_quadrature_over_resolution']:.0f}x below the "
      f"{RES_X:g} resolution the separability verdicts are read at. "
      "That is what licenses the sweep to run on the quadrature.")
# and the same comparison for the class whose profile is the opposite way up
ref_ser = PlugExact(161).run("series", BASE["kappa"], BASE["lam"], BASE["d"],
                             BASE["nprime"], n=BASE["n"], rtol=1e-11)[0]
Xser, _, _ = BED.march("series", BASE["kappa"], BASE["lam"], BASE["d"],
                       BASE["nprime"], n=BASE["n"], nt=200)
M["bed_vs_quadrature_series"] = float(np.max(np.abs(Xser - ref_ser)))
print(f"same, series (ascending profile): {M['bed_vs_quadrature_series']:.2e}")
'''))

cells.append(code(r'''"""(2) Every device against the PAPER'S OWN closed forms, eqs. (15), (19), (22),
(28)-(31).  These are exact solutions of the systems being solved, printed in
1972; they share no code with anything above."""
kap, lam = 1.2, 2.5          # first-order reaction, first-order decay, n' = 0
rows = []

# eq. (15): batch fluid.  k'' t_obs = kap_b, k_d t_obs = lam
Xb, ab = run_batch("independent", kap, lam, 1.0, 0.0, n=1.0)
cb = 1.0 - Xb
c_inf = np.exp(-kap / lam)
lhs = np.log(np.log(cb / c_inf)); rhs = np.log(np.log(1.0 / c_inf)) - lam * TH
M["eq15_maxdev"] = float(np.max(np.abs(lhs - rhs)))
rows.append(["(15)", "batch fluid, Fig. 2", f"{M['eq15_maxdev']:.2e}"])

# eq. (19): mixed flow, constant rate
Xm, am = run_mixed("independent", kap, lam, 1.0, 0.0, n=1.0)
cm = 1.0 - Xm
M["eq19_maxdev"] = float(np.max(np.abs(np.log(1.0 / cm - 1.0)
                                       - (np.log(kap) - lam * TH))))
rows.append(["(19)", "mixed flow, Fig. 3", f"{M['eq19_maxdev']:.2e}"])

# eq. (22): plug flow, constant rate - the pymrm BED, not the quadrature
Xp, Ap, _ = BED.march("independent", kap, lam, 1.0, 0.0, n=1.0, nt=200)
cp = 1.0 - Xp
M["eq22_maxdev_pymrm"] = float(np.max(np.abs(np.log(np.log(1.0 / cp))
                                             - (np.log(kap) - lam * TH))))
rows.append(["(22)", "plug flow, Fig. 4 - solved on the pymrm bed",
             f"{M['eq22_maxdev_pymrm']:.2e}"])

# eqs. (28)-(31): the held-C_A device at d = 0, 1, 2, 3.  lam_c < 1 so that the
# d = 0 branch, where the activity reaches zero at k'_d t = 1, stays alive across
# the whole window - the branch is exercised, not stepped over.
cs, lam_c = BASE["c_set"], 0.8
Qd = (1.0 - cs) / (kap * cs)          # (C_A0 - C_A)/k' in tau'_ref units
for dv, form in ((0, "1/tau'"), (1, "ln tau'"), (2, "tau'"), (3, "(tau')^2")):
    tau, aa = run_const_c("independent", kap, lam_c, float(dv), 0.0, cs, n=1.0)
    if dv == 0:
        lhs, rhs = 1.0 / tau, (1.0 / Qd) * (1.0 - lam_c * TH)
    elif dv == 1:
        lhs, rhs = np.log(tau), np.log(Qd) + lam_c * TH
    elif dv == 2:
        lhs, rhs = tau, Qd * (1.0 + lam_c * TH)
    else:
        lhs, rhs = tau ** 2, Qd ** 2 * (1.0 + 2.0 * lam_c * TH)
    dev = float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs)))
    M[f"eq{28+dv}_reldev"] = dev
    M[f"a_end_d{dv}_branch"] = float(aa[-1])
    rows.append([f"({28+dv})", f"held C_A, d = {dv} ({form} linear in t)",
                 f"{dev:.2e}"])
display(Markdown("**The paper's own closed forms, met by the reactors solved "
                 "here** (max deviation; relative for eqs. 28-31)\n\n"
                 + md_table(pd.DataFrame(rows, columns=["printed eq.",
                                                        "device / branch",
                                                        "deviation"]))))
print("All four printed branches of the deactivation order (d = 0, 1, 2, 3) are "
      "exercised on a LIVE activity: a(theta = 1) = "
      + ", ".join(f"{M[f'a_end_d{j}_branch']:.3f}" for j in range(4))
      + " for d = 0, 1, 2, 3.  d = 0 is the branch on which the activity reaches "
        "zero in finite time (at k'_d t = 1); the window stops short of it.")
'''))

cells.append(code(r'''"""(4) The paper's pseudo-steady assumption (p. 269), dropped and costed.

eps = (fluid residence time)/(observation window).  The bed is re-solved with
eps dc/dtheta retained, on the same operators; the deviation from the
pseudo-steady answer is root-found against the conversion resolution."""
BEDT = Bed(m=80)                     # eps study: the O(eps) term, not the grid
BEDT._SP = BEDT.jac_sparsity()
X_PS, _, _ = BEDT.march("parallel", BASE["kappa"], BASE["lam"], BASE["d"],
                        BASE["nprime"], n=BASE["n"], nt=200)   # eps = 0 reference
def pseudo_steady_cost(eps, cls="parallel"):
    Xt, _ = BEDT.transient(cls, BASE["kappa"], BASE["lam"], BASE["d"],
                           BASE["nprime"], eps, n=BASE["n"])
    return float(np.max(np.abs(Xt - X_PS)))

eps_grid = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2])
cost = np.array([pseudo_steady_cost(e) for e in eps_grid])
M["pseudo_steady_cost_eps1e-3"] = float(cost[2])
M["pseudo_steady_cost_eps3e-2"] = float(cost[-1])
M["pseudo_steady_cost_slope"] = float(np.polyfit(np.log(eps_grid[:4]),
                                                 np.log(cost[:4]), 1)[0])
eps_star = brentq(lambda e: pseudo_steady_cost(np.exp(e)) - RES_X,
                  np.log(1e-4), np.log(3e-2), xtol=1e-3)
M["pseudo_steady_eps_star"] = float(np.exp(eps_star))
fig, ax = plt.subplots(figsize=(5.0, 3.1))
ax.loglog(eps_grid, cost, "o-", color="#0072B2", lw=1.5, ms=4)
ax.axhline(RES_X, color="0.4", ls="--", lw=1,
           label=f"conversion resolution {RES_X:g}")
ax.axvline(M["pseudo_steady_eps_star"], color="#D55E00", ls=":", lw=1.2,
           label=rf"$\varepsilon^*$ = {M['pseudo_steady_eps_star']:.2e}")
ax.set_xlabel(r"$\varepsilon$ = fluid residence time / observation window")
ax.set_ylabel("max error in $X_A$ from the\npseudo-steady assumption")
ax.grid(alpha=0.25, which="both", lw=0.5); ax.legend(fontsize=8, frameon=False)
fig.tight_layout(); plt.show()
print(f"cost is first order in eps (slope {M['pseudo_steady_cost_slope']:.2f}), "
      f"and crosses the {RES_X:g} conversion resolution at "
      f"eps* = {M['pseudo_steady_eps_star']:.2e}.")
print(f"Read as the paper reads it (p. 267): with a fluid residence time of 1 s, "
      f"the assumption is worth its own error bar once the run lasts longer than "
      f"{1.0/M['pseudo_steady_eps_star']:.0f} s = "
      f"{1.0/M['pseudo_steady_eps_star']/60:.0f} min - which is exactly the "
      f"boundary the paper draws in words: \"Deactivation in the order of "
      f"minutes or longer can use the fixed batch of solids.\"")
BREAKS.append(("pseudo_steady_cost_eps1e-3", M["pseudo_steady_cost_eps1e-3"],
               "raise eps 30x, to 3e-2", M["pseudo_steady_cost_eps3e-2"],
               "the assumption's cost is first order in eps and does move"))
'''))

cells.append(md(r"""### The break table

Every reported metric needs something that moves it. Where a metric is an
identity that cannot fail — and several here are, deliberately — it is named
as structural and paired with a number that can."""))

cells.append(code(r'''"""Defect injection: break something on purpose, watch the number move."""
kap, lam, dv, npv = BASE["kappa"], BASE["lam"], BASE["d"], BASE["nprime"]

# (a) wrong nu in construct_div: cylindrical instead of Cartesian weight coord
class BedNu(Bed):
    def __init__(self, m=200):
        super().__init__(m)
        conv, conv_bc = construct_convflux_upwind((m, 1), self.s_f, self.s_c,
                                                  self.bc, v=1.0)
        div = construct_div((m, 1), self.s_f, nu=1)      # WRONG: cylindrical
        self.div = div
        self.Lop = (div @ conv).tocsc()
        self.bvec = np.asarray((div @ conv_bc).todense()).ravel()
Xnu, _, _ = BedNu(200).march("parallel", kap, lam, dv, npv, n=BASE["n"], nt=200)
BREAKS.append(("bed_vs_quadrature_production", M["bed_vs_quadrature_production"],
               "construct_div(nu=1) - cylindrical, not the Cartesian weight coord",
               float(np.max(np.abs(Xnu - ref_par))), "wrong geometry flag"))

# (b) first-order upwind: drop the van Leer deferred correction
Xtvd, _, _ = Bed(200, tvd=False).march("parallel", kap, lam, dv, npv,
                                       n=BASE["n"], nt=200)
BREAKS.append(("bed_vs_quadrature_production", M["bed_vs_quadrature_production"],
               "tvd=False - plain first-order upwind",
               float(np.max(np.abs(Xtvd - ref_par))), "deferred correction removed"))

# (c) Euler instead of Heun on the activity - ORDER against ORDER, both
# observed on the same refinement pair.  Setting an order against an ERROR (the
# obvious way to write this row) cannot fail, because they are different
# quantities: "moves? yes" would be automatic and nothing would be tested.
BREAKS.append(("bed_time_order", M["bed_time_order"],
               "scheme='euler' on the activity march (same 50->100 pair, m=400)",
               M["bed_time_order_euler"],
               "observed order against observed order: 2 -> 1"))

# (d) the inlet boundary condition flipped to zero-gradient
class BedBC(Bed):
    def __init__(self, m=200):
        super().__init__(m)
        self.bc = ({"a": 1.0, "b": 0.0, "d": 0.0},   # WRONG: no feed imposed
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        conv, conv_bc = construct_convflux_upwind((m, 1), self.s_f, self.s_c,
                                                  self.bc, v=1.0)
        self.Lop = (self.div @ conv).tocsc()
        self.bvec = np.asarray((self.div @ conv_bc).todense()).ravel()
try:
    Xbc, _, _ = BedBC(200).march("parallel", kap, lam, dv, npv, n=BASE["n"], nt=40)
    dbc = float(np.max(np.abs(Xbc - ref_par)))
except Exception:
    dbc = float("inf")
BREAKS.append(("bed_vs_quadrature_production", M["bed_vs_quadrature_production"],
               "inlet bc {a:1,b:0,d:0} - feed concentration not imposed",
               dbc, "boundary condition on the outward normal. NOT independent "
                    "evidence from row (a): an un-fed bed and a bed with the "
                    "wrong nu both collapse to the same degenerate answer, so "
                    "the two rows return the same number"))

# (e) the four rate forms, swapped one for another (the wrong-class injection)
for cfg in ("plug", "mixed"):
    y_true = respond(cfg, "parallel", kap, lam, dv, npv, BASE["sig"],
                     c_set=BASE["c_set"])
    y_wrong = respond(cfg, "series", kap, lam, dv, npv, BASE["sig"],
                      c_set=BASE["c_set"])
    BREAKS.append((f"sep_{cfg}_par_vs_ser", M[f"sep_{cfg}_par_vs_ser"],
                   "use eq. (5) where eq. (4) is true, WITHOUT refitting",
                   float(np.sqrt(np.mean((y_wrong - y_true) ** 2))),
                   "the unfitted gap - the fitted one is the metric itself"))

# (f) transcription: a mis-read subscript in eq. (5), C_R read as C_A
y_true = respond("plug", "series", kap, lam, dv, npv, 0.0)
y_typo = respond("plug", "parallel", kap, lam, dv, npv, 0.0)
BREAKS.append(("printed_chain_nonzero", 0.0,
               "eq. (5) transcribed with C_A for C_R (one subscript)",
               float(np.sqrt(np.mean((y_typo - y_true) ** 2))),
               "reported as the response gap; the symbolic residual is the guard"))

# (g) the n' = 0 collapse, broken by giving the decay a concentration
BREAKS.append(("degen_nprime0_plug", M["degen_nprime0_plug"],
               "give the decay a concentration order n' = 0.7",
               M["degen_nprime0_off_condition"],
               "the collapse is exactly at n' = 0 and nowhere else"))

# (h) the fitter itself: double the multistart on the headline entry
y_true = respond("plug", "parallel", kap, lam, dv, npv, BASE["sig"],
                 c_set=BASE["c_set"])
big = _seeds("series", kap, lam) + [[np.log(kap * f), np.log(lam * g), dvv, npvv]
                                    for f in (0.5, 2.0) for g in (0.25, 4.0)
                                    for dvv in (0.5, 2.5) for npvv in (0.25, 2.5)]
D_big = fit_class("plug", y_true, "series", kap, lam, seeds=big, keep=6)[0]
BREAKS.append(("sep_plug_par_vs_ser", M["sep_plug_par_vs_ser"],
               f"multistart widened from {len(_seeds('series', kap, lam))} to "
               f"{len(big)} seeds", D_big,
               "ROBUSTNESS (must NOT move): a separability number is an UPPER "
               "bound, and this bounds the bound - a lower value here would mean "
               "the headline overstates the separation"))
M["sep_plug_par_vs_ser_wideseed"] = float(D_big)

# (i) the resolution the verdicts are read at - ONE ROW PER DEVICE, so that the
# three counts carrying the page's negative headline (batch, mixed, constC) each
# have a mover of their own.  The alternative resolution is chosen per device to
# be one that actually moves that device's count.
for cfg, alt in (("plug", 1e-4), ("batch", 1e-4), ("mixed", 1e-4)):
    BREAKS.append((f"sep_{cfg}_n_separable", M[f"sep_{cfg}_n_separable"],
                   f"read the {cfg} matrix at resolution {alt:g} instead of "
                   f"{M[f'sep_threshold_{cfg}']:g}",
                   float(sum(1 for i in range(4) for j in range(4)
                             if i != j and SEP[cfg][i, j] > alt)),
                   "the count is a function of the stated resolution, by design"))

# ... and the held-C_A count needs a PHYSICAL injection, because no resolution
# moves it: every pair in that device is exactly degenerate at the base point.
# Degeneracy 3's fourth member holds only at d = 1 (the poison consumed per unit
# time goes as a^(d-1)), so d = 2 is where side-by-side leaves the collapse.
# Where an identity of sec. 4 still applies the twin is used, so this costs
# three fits, not sixteen.
def constC_count_at(d_inj, res):
    b = dict(BASE, d=d_inj); cnt = 0
    for i, truth in enumerate(CLASSES):
        y = respond("constC", truth, b["kappa"], b["lam"], b["d"], b["nprime"],
                    b["sig"], n=BASE["n"], c_set=b["c_set"])
        for j, cand in enumerate(CLASSES):
            if i == j:
                continue
            v, _ = twin_residual("constC", truth, cand, b, BASE["n"])
            if v is None:
                v = fit_class("constC", y, cand, b["kappa"], b["lam"],
                              n=BASE["n"], c_set=b["c_set"],
                              tag=("constC-d2", truth, cand))[0]
            cnt += int(v > res)
    return float(cnt)
CONSTC_ALT_RES = 1e-4
BREAKS.append(("sep_constC_n_separable", M["sep_constC_n_separable"],
               f"move the held-C_A device to d = 2, where the a^(d-1) "
               f"cancellation that pins C_P fails, AND read at "
               f"{CONSTC_ALT_RES:g}",
               constC_count_at(2.0, CONSTC_ALT_RES),
               "COMPOUND, and deliberately so: at the page's own 1e-2 "
               "resolution this count stays 0 for every d tried, because the "
               "collapse of parallel/series/independent is exact at any d and "
               "only side-by-side ever leaves it"))

df = pd.DataFrame(BREAKS, columns=["metric", "as built", "injected defect",
                                   "with the defect", "note"])
df["kind"] = ["robustness" if str(nn).startswith("ROBUSTNESS") else "mover"
              for nn in df["note"]]
df["moves?"] = [("yes" if (not np.isfinite(b) or abs(b - a) >
                           max(1e-14, 1e-6 * max(abs(a), abs(b)))) else "no")
                for a, b in zip(df["as built"], df["with the defect"])]
display(Markdown("**Defect-injection table**\n\n" + md_table(df)))
movers = df[df["kind"] == "mover"]
assert (movers["moves?"] == "yes").all(), movers[movers["moves?"] != "yes"]
rob = df[df["kind"] == "robustness"]
for a, b in zip(rob["as built"], rob["with the defect"]):
    assert b >= 0.9 * a, f"robustness row improved: {a} -> {b}"
print(f"{len(movers)} injected defects, every one of which moves its metric; "
      f"{len(rob)} robustness row(s), which must NOT move and do not "
      f"(worst change {max((abs(b/a - 1) for a, b in zip(rob['as built'], rob['with the defect'])), default=0.0):.1e}).")
M["break_rows"] = float(len(movers))
M["break_robustness_rows"] = float(len(rob))
'''))

cells.append(code(r'''"""Coverage: which injected defect moves each reported metric, key for key."""
CHK_SYM = "sec. 1 symbolic identities"
COVERAGE = {}
for kk in M:
    if kk.startswith("span_"):
        COVERAGE[kk] = ("break (e)/(f): the response gap under a swapped rate "
                        "form is exactly a change of these spans")
    elif kk == "sep_plug_par_vs_ser_wideseed":
        COVERAGE[kk] = "break (h) is this metric"
    elif kk.startswith("sep_") and "_vs_" in kk:
        COVERAGE[kk] = ("break (e) directly (unfitted class swap) + (h) "
                        "(multistart width); entries marked * in the matrix are "
                        "identities of sec. 4, not measurements")
    elif kk.startswith("sep_") and kk.endswith("n_separable"):
        COVERAGE[kk] = (
            ("break (i-constC), which is COMPOUND and has to be: NO resolution "
             "moves this count, because after the twins every off-diagonal "
             "entry of this device is at double-precision zero. The row "
             "therefore moves the DEVICE (d = 2, where the a^(d-1) "
             "cancellation fails) as well as the resolution, and its own note "
             "says so")
            if kk.startswith("sep_constC") else
            ("break (i-" + kk.split("_")[1] + ") directly: this device's own "
             "matrix re-read at another resolution"))
    elif kk.startswith("sep_") and kk.endswith("best_pair"):
        COVERAGE[kk] = "break (e),(h) class (same fits)"
    elif kk.startswith("mixed_boundary"):
        COVERAGE[kk] = ("break: the boundary row pair itself - crossing "
                        "d = n'/n switches it by "
                        f"{M['mixed_boundary_np1_d0.5']/M['mixed_boundary_np1_d1.5']:.0e}"
                        " (that ratio is the executed row, not an estimate)")
    elif kk.startswith("bed_spatial"):
        COVERAGE[kk] = "break (a) nu=1, (b) tvd=False, (d) inlet bc"
    elif kk.startswith("bed_time"):
        COVERAGE[kk] = ("break (c) Euler for Heun, order against order on the "
                        "same refinement pair")
    elif kk.startswith("bed_vs_quadrature"):
        COVERAGE[kk] = "break (a),(b),(d) directly"
    elif kk.startswith("bed_a_inlet_cellcentre"):
        COVERAGE[kk] = ("the grid refinement in sec. 6 IS the injection: the "
                        "cell-centre error halves per refinement against the "
                        "exact face value, which is what makes the reported "
                        "inlet activity exact rather than sampled")
    elif kk == "bed_profile_time_sensitivity":
        COVERAGE[kk] = ("the nt doubling beside it, on floats and on BOTH "
                        "classes: it is the number that shows the outlet "
                        "reading is NOT purely spatial, which is why the "
                        "outlet carries a second extrapolation. Read on the "
                        "series class alone, and off a 6-decimal string, it "
                        "was 14x too small and hid a wrong printed value")
    elif kk == "bed_inlet_time_sensitivity":
        COVERAGE[kk] = ("the same nt doubling, on the inlet cell centres: the "
                        "control showing that the sec.-6 SPATIAL study and its "
                        "observed orders are uncontaminated by the step count")
    elif kk.startswith("bed_a_outlet_spaceonly"):
        COVERAGE[kk] = ("break: halve NT_CONV. This metric is the space "
                        "extrapolation ALONE, it does move with the step "
                        "count, and that is precisely why it is not the "
                        "reported outlet")
    elif kk.startswith("bed_a_outlet_end") or kk.startswith("bed_profile_contrast"):
        COVERAGE[kk] = ("break: this class's own ROBUSTNESS row - the time "
                        "extrapolation read off the other nt pair, which must "
                        "not move it; plus break (f), which inverts the "
                        "profile these numbers measure")
    elif kk.startswith("bed_a_outlet_altpair"):
        COVERAGE[kk] = ("it IS this class's robustness row value: the reported "
                        "outlet read off the other nt pair, and the row "
                        "asserts the two agree")
    elif kk.startswith("bed_a_outlet_time_order"):
        COVERAGE[kk] = ("the nt ladder beside it, which is what this order is "
                        "read from; break (c) prices the same knob on the "
                        "conversion history, where Heun's 1.97 becomes Euler's "
                        "1.01 on one refinement pair")
    elif (kk.startswith("bed_a_outlet_time_bias")
          or kk.startswith("bed_a_outlet_space_spread")
          or kk.startswith("bed_a_outlet_residual_bound")):
        COVERAGE[kk] = ("the two refinements beside it: this is the SIZE of "
                        "the error the second extrapolation removes, and the "
                        "bound that fixes how many decimals may be printed. "
                        "Halving NT_CONV moves the first of them")
    elif kk.startswith("bed_a_") or kk.startswith("bed_profile"):
        COVERAGE[kk] = ("break (f): swapping the subscript in eq. (5) inverts the "
                        "profile, which is what these numbers measure")
    elif kk.startswith("pseudo_steady"):
        COVERAGE[kk] = "break: the eps sweep is itself the injection (1e-3 -> 3e-2)"
    elif kk.startswith("constC_span"):
        COVERAGE[kk] = ("break (constC): the demanded precision on k'_d moves "
                        "it, and so does the unprinted starting level C_LO")
    elif kk.startswith("constC_sbs"):
        COVERAGE[kk] = "break (constC-sbs): sigma = 0 removes the bias entirely"
    elif kk == "constC_misfit_at_span_0p5":
        COVERAGE[kk] = "same family: the span is the knob and it is swept"
    elif kk.startswith("degen_nprime0"):
        COVERAGE[kk] = "break (g) directly (n' = 0.7)"
    elif kk == "degen_mixed_series_is_parallel_rtol1e-13":
        COVERAGE[kk] = ("the rtol refinement beside it: this metric and its "
                        "production-tolerance twin differ by the SOLVER "
                        "tolerance only, which is what shows the identity exact")
    elif kk == "degen_mixed_series_is_parallel":
        COVERAGE[kk] = ("STRUCT: an exact algebraic identity, and the printed "
                        "value is the ODE tolerance rather than a separation "
                        "(rtol 1e-9 -> 1e-13 moves it by four orders). Its "
                        "companion that CAN move is mixed_boundary_np1_d0.5, "
                        "the same identity evaluated where it does not hold")
    elif kk == "degen_constC_par_ser_ind":
        COVERAGE[kk] = ("STRUCT: exact by construction (k'_d is a constant). "
                        "Companion that runs through the device itself: "
                        "span_constC_par. constC_span_needed_prec5pct is the "
                        "PRICE of the degeneracy, but it is an analytic "
                        "least-squares argument and shares no code with "
                        "run_const_c, so it cannot stand as the detector")
    elif kk.startswith("degen_sbs_is_parallel_famA"):
        COVERAGE[kk] = ("STRUCT: degeneracy 4, exact wherever it holds. "
                        "Companion/mover: degen_batch_sbs_off_condition, the "
                        "same comparison off n' = n, d = 1")
    elif kk.startswith("degen_sbs_is_parallel_famB"):
        COVERAGE[kk] = ("degeneracy 5, exact in the batch and the bed only - "
                        "the mixed/constC entries of this same family ARE the "
                        "break row that moves it, and so is "
                        "degen_sbs_famB_off_condition")
    elif kk == "degen_sbs_famB_off_condition":
        COVERAGE[kk] = "the mover for degeneracy 5 inside the bed itself"
    elif kk == "degen_batch_sbs_is_parallel_on_condition":
        COVERAGE[kk] = "break: the off-condition row (d = 2) in the same table"
    elif kk == "degen_batch_sbs_off_condition":
        COVERAGE[kk] = "the mover for the row above"
    elif kk == "degen_nprime0_off_condition":
        COVERAGE[kk] = "break (g) is this metric"
    elif kk in ("printed_chain_identities", "printed_chain_nonzero"):
        COVERAGE[kk] = (f"STRUCT: {CHK_SYM} cannot fail once the transcription is "
                        "right - that is what they are for. Companion: break (f), "
                        "which prices one mis-read subscript in the response")
    elif kk.startswith("eq24"):
        COVERAGE[kk] = ("break: reading eq. (24) with 1/k sends the residual to 0 "
                        "exactly (first BREAKS row)")
    elif kk.startswith("eq15") or kk.startswith("eq19") or kk.startswith("eq22") \
            or kk.startswith("eq28") or kk.startswith("eq29") \
            or kk.startswith("eq30") or kk.startswith("eq31"):
        COVERAGE[kk] = ("break (a),(b),(c),(d) for eq. (22), which is solved on the "
                        "pymrm bed; the lumped ones are ODE-vs-closed-form and move "
                        "with any rate-form defect, e.g. break (f)")
    elif kk.startswith("a_end_d"):
        COVERAGE[kk] = ("branch-liveness guard for eqs. (28)-(31): moves with "
                        "lam_c and with d, and is what shows d = 0 is exercised "
                        "on a live activity rather than stepped over")
    elif kk == "fits_capped":
        COVERAGE[kk] = ("a COUNTER, not an identity: sep_matrix now passes a tag "
                        "to fit_class, so a capped sweep fit would be recorded "
                        "and the assertion beside it would fire. It is 0 because "
                        "no fit hits the budget, not because nothing can")
    elif kk.startswith("sep_threshold"):
        COVERAGE[kk] = (
            ("STRUCT for this device: a stated design constant, and re-reading "
             "the held-C_A matrix at another resolution moves NOTHING, because "
             "every off-diagonal entry of it is exactly 0. Break (i-constC) "
             "does re-read it, but only as half of a compound injection; what "
             "moves that count is the other half, d = 2")
            if kk.endswith("constC") else
            ("break (i-" + kk.split("_")[2] + "): the stated resolution moves "
             "the counts read at it; the threshold IS that resolution"))
    elif kk.startswith("sep_stall"):
        COVERAGE[kk] = ("the multistart's own shortfall on a problem whose "
                        "exact minimum is 0. Its knob is break (h)'s - the "
                        "seed count - but note that (h) is run as a ROBUSTNESS "
                        "row, on a separation that must NOT move; widening the "
                        "seeds can only lower a stall, never raise one, so the "
                        "stalls are upper bounds and the displacement table "
                        "beside them says where the fit stopped")
    elif kk.startswith("sep_identity") or kk.startswith("sep_optimiser"):
        COVERAGE[kk] = ("the fit-vs-identity comparison in sec. 5 is itself the "
                        "injection: it is the multistart alone, measured against "
                        "the algebra, entry by entry")
    elif kk.startswith("sep_floor"):
        COVERAGE[kk] = ("the printed floor check in sec. 5: 20 x the stall "
                        "against each device's resolution, on all four devices")
    elif kk in ("break_rows", "break_robustness_rows"):
        COVERAGE[kk] = "STRUCT: a count of the rows of the table above it"
    else:
        COVERAGE[kk] = "UNCOVERED"
missing = [kk for kk, v in COVERAGE.items() if v == "UNCOVERED"]
assert set(COVERAGE) == set(M), (set(M) - set(COVERAGE), set(COVERAGE) - set(M))
assert not missing, missing
print("break-row coverage (which injected defect moves each reported metric):")
for kk in M:
    print(f"  {kk:<44s} {COVERAGE[kk]}")
ABS_FLOOR = 1e-12
below = sorted(kk for kk, v in M.items() if abs(v) < ABS_FLOOR)

def floor_note(kk):
    """For a metric below ABS_FLOOR: what KIND of number it is, what it cannot
    detect, and an ABOVE-FLOOR companion metric that can.  Not every sub-floor
    number here is a degeneracy or an identity - some are counts and some are
    numerical deviations - and saying so is the point of this table."""
    if kk.startswith("degen_nprime0"):
        return ("exact degeneracy 1",
                "nothing: 0 by code construction, since the class enters only "
                "through f^n' and _pow(x, 0) returns 1.0",
                "degen_nprime0_off_condition")
    if (kk.startswith("degen_sbs_is_parallel_famA")
            or kk == "degen_batch_sbs_is_parallel_on_condition"):
        return ("exact degeneracy 4", "any change that keeps p == c",
                "degen_batch_sbs_off_condition")
    if kk.startswith("degen_sbs_is_parallel_famB"):
        return ("exact degeneracy 5", "any change that keeps p a power of c",
                "degen_sbs_famB_off_condition")
    if kk == "degen_constC_par_ser_ind":
        return ("exact degeneracy 3",
                "anything at all inside one held-C_A run: k'_d is a constant "
                "for three of the four classes by construction",
                "span_constC_par")
    if kk == "degen_mixed_series_is_parallel_rtol1e-13":
        return ("exact degeneracy 2, integrated tightly",
                "nothing: it measures the ODE tolerance",
                "mixed_boundary_np1_d0.5")
    if kk == "printed_chain_nonzero":
        return ("count of failed symbolic identities",
                "a transcription error that is self-consistent across the whole "
                "printed chain",
                "sep_plug_par_vs_ser")
    if kk == "fits_capped":
        return ("count of fits that hit the evaluation budget",
                "nothing about the physics; it is a cost counter, and it is now "
                "wired to fire (a tag is passed from sep_matrix)",
                "sep_optimiser_worst_overstatement")
    if kk == "sep_floor_binds_anywhere":
        return ("boolean: does 20 x the stall ever exceed a resolution",
                "nothing - it is 0 because the floor never binds, which is why "
                "the floor was removed",
                "sep_floor_worst_ratio")
    if kk.endswith("n_separable"):
        return ("count of separable class pairs - THE PAGE'S NEGATIVE HEADLINE, "
                "not an identity",
                "a separation smaller than the stated resolution",
                "sep_plug_best_pair" if kk.startswith("sep_plug")
                else ("span_constC_par" if kk.startswith("sep_constC")
                      else f"sep_{kk.split('_')[1]}_best_pair"))
    if kk.endswith("best_pair"):
        return ("largest off-diagonal entry of a matrix that is entirely "
                "identity-backed",
                "any confusion in a device where every pair is exactly "
                "degenerate", "span_constC_par")
    if kk.startswith("sep_stall"):
        return ("optimiser stall on a problem whose exact minimum is 0",
                "nothing on this device: the multistart reached the truth",
                "sep_stall_worst")
    if kk.startswith("sep_"):
        cfg = kk.split("_")[1]
        return ("min(fit, identity twin) where an identity of sec. 4 applies",
                "any separation, because the exact minimum of this entry is zero",
                "span_constC_par" if cfg == "constC"
                else f"sep_{cfg}_best_pair")
    if kk.startswith("eq") and (kk.endswith("reldev") or kk.endswith("maxdev")):
        return ("NUMERICAL deviation of an ODE solution from the paper's own "
                "closed form - not a symbolic identity",
                "an error that the closed form shares, i.e. a mis-transcription "
                "of the printed solution itself", "eq22_maxdev_pymrm")
    if kk.startswith("degen_"):
        return ("exact algebraic degeneracy between the printed rate forms",
                "any change that preserves the identity",
                "degen_nprime0_off_condition")
    if kk.startswith("mixed_boundary"):
        return ("identity twin inside d >= n'/n, optimiser bound outside it",
                "the loss of exactness at the boundary, which is what the row "
                "below the line measures", "mixed_boundary_np1_d0.5")
    if kk.startswith("constC"):
        return ("property of the held-C_A device, which is exactly degenerate",
                "any class distinction inside one run",
                "constC_sbs_recovered_d_base")
    if kk.startswith("bed_") or kk.startswith("span_") or kk.startswith("a_end"):
        return ("computed property of the bed or of a response span",
                "an error shared by the independent quadrature",
                "bed_vs_quadrature_production")
    return None

fl_rows, unnamed, bad_comp = [], [], []
for kk in below:
    note = floor_note(kk)
    if note is None:
        unnamed.append(kk); continue
    kind, blind, comp = note
    if comp not in M or abs(M[comp]) <= ABS_FLOOR:
        bad_comp.append((kk, comp, M.get(comp)))
        continue
    fl_rows.append([kk, kind, blind, f"{comp} = {M[comp]:.2e}"])
assert not unnamed, f"sub-floor metrics with no stated kind: {unnamed}"
assert not bad_comp, f"sub-floor metrics whose companion is not above floor: {bad_comp}"
kinds = sorted({r[1] for r in fl_rows})
M["subfloor_metrics"] = float(len(below))
M["subfloor_kinds"] = float(len(kinds))
print(f"\n{len(below)} of {len(M)} metrics sit below check_agreement's "
      f"ABS_FLOOR = {ABS_FLOOR:g} and are therefore documented rather than "
      f"protected by CI. They are NOT all identities: they fall into "
      f"{len(kinds)} kinds -")
for kd in kinds:
    print(f"  {sum(1 for r in fl_rows if r[1] == kd):3d}  {kd}")
display(Markdown("**Every sub-floor metric, what it cannot detect, and its "
                 "above-floor companion.** A companion is only useful if it "
                 "runs through the SAME code the sub-floor metric comes from, "
                 "so that a defect there moves it. The held-C_A family named "
                 "`constC_span_needed_prec5pct` until this was checked - a "
                 "numpy least-squares fit to analytic power laws that never "
                 "calls `run_const_c`, and so could not have detected a defect "
                 "in the device at all; those rows now name numbers that are "
                 "outputs of that device.\n\n"
                 + md_table(pd.DataFrame(fl_rows, columns=[
                     "metric", "kind", "what it cannot detect",
                     "above-floor companion (asserted > ABS_FLOOR)"]))))
COVERAGE["subfloor_metrics"] = COVERAGE["subfloor_kinds"] = (
    "bookkeeping of the table above; every entry in it is asserted to name a "
    "companion that is itself above the floor")
assert set(COVERAGE) == set(M), (set(M) - set(COVERAGE), set(COVERAGE) - set(M))
# wall clock is deliberately neither printed nor reported: it is machine
# dependent, and two executions of this page must give identical content.
report_agreement("B2.3", M)
'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**To the four rate equations, nothing** — eqs. (4)–(7) are four lines of
algebra, and three of the paper's four contactings are lumped systems that
need no discretisation at all. This page says so rather than dressing them up.
What the reimplementation adds is in three places, and only one of them is the
solver:

- **The question the paper asks becomes measurable.** "Which reactor can tell
  the classes apart" is answered by 64 fits, not by an argument, and it comes
  out the opposite way round from the taxonomy: **every one of the four devices
  is exactly degenerate somewhere**, and the reason in each case is a one-line
  identity the paper does not state. Two of the four collapse parallel onto
  series or onto independent, which destroys the taxonomy outright; the other
  two — the batch of fluid and the plug-flow bed — collapse only side-by-side
  onto parallel, and do so for *every* value of the poison group $\sigma$. The
  identity for the mixed reactor — series$(n',d) \equiv$ parallel$(nn',\,d+n')$
  — is the sharp form of the whole result, and it carries an analytic boundary,
  $d < n'/n$, below which the exact twin ceases to exist; the boundary is
  *predicted*, and the measurement then shows what it is and is not worth
  (necessary for the mixed reactor to separate the two classes, nowhere near
  sufficient).
- **And the identities are used, not just admired.** Each one names a parameter
  vector at which some candidate reproduces some truth exactly. Those vectors
  are evaluated directly and reported beside the fits — because on a ninth of
  the matrix the multistart, seeded next to the twin and free to reach it,
  stopped short and would have printed an exactly-zero separation as a measured
  one. That is the asymmetry this page is built on: a non-separability result
  is a proof, a separability result is only "no better fit was found".
- **And the honest negative.** At a 0.2-percentage-point conversion resolution,
  one of the twelve class pairs is separable in one of the four devices. A page
  that only reported where separation is easy would have reported nothing at
  all here — which is why the measurement was worth making.
- **The bed is where pymrm is load-bearing, and it changes the
  recommendation.** The plug-flow contacting is the only one with a spatial
  coordinate; it is *not* the only one free of exact degeneracy, as §4 shows —
  and its outlet history is still a hard measurement. Its
  *profile* is not: descending for parallel, ascending for series, contrasts of
  opposite sign, with both inlet values known in closed form
  ($a(0,\theta) = e^{-\lambda\theta}$ under parallel decay, and identically 1
  under series, because no R has been made at the feed face). The page's
  practical advice is therefore about what to
  measure, not which vessel to buy, and it is the same signature `B2.2` found
  in the coke. And because the finite-volume bed keeps the fluid accumulation
  term, the **pseudo-steady assumption the paper states on p. 269 can be
  dropped and priced**: first order in $\varepsilon$, crossing a 0.2 %
  conversion resolution at a run length that lands on the "minutes or longer"
  the paper writes in words on p. 267.
- **One printed defect, settled two ways.** Eq. (24)'s $1/k_d$ should be
  $1/k$: symbolically the gap is exactly $\ln(k_d/k)$, and dimensionally the
  printed argument is in seconds where $\tau'$ is in g cat·sec/liter. Reported,
  not repaired.

**What it cannot add.** There is no measurement in this paper, so nothing here
is validated in the sense of agreeing with an experiment — every number is
reproduction of printed algebra, an exact identity, or a computed property of
the printed rate forms. And a separability number is an upper bound produced
by an optimiser: the *non*-separability results are proofs, the separability
results are "no better fit was found"."""))

# ------------------------------------------------------------------- reuse --
cells.append(md(r"""## Reuse

**If you are designing an experiment to find a deactivation mechanism**, the
short version of this page is:

1. **Check $n'$ first.** If the decay has no concentration dependence, there
   is only one class and no experiment can find another. Everything below is
   about $n' > 0$.
2. **A well-mixed fluid throws the class away.** In any constant-flow mixed
   reactor, series deactivation of order $(n', d)$ is *exactly* parallel
   deactivation of order $(nn', d+n')$ — same conversion history, at every
   time, for any $n$. Only below $d = n'/n$ does the exact twin cease to exist,
   and even there the impersonation stays good until you are roughly half an
   order below the line.
3. **Side-by-side is not a separate class in a batch or in a bed.** The poison
   profile is a *power* of the reactant profile there, so eq. (6) with poison
   group $\sigma$ is eq. (4) with $n' \to \lambda\sigma/\kappa$ — exactly, for
   every $\sigma$, at $n = n' = d = 1$. No experiment in either device can
   distinguish them, and the plug-flow bed is no exception. (It fails in the
   two well-mixed devices, where the profiles are not powers of one another.)
4. **Levenspiel's recommended device does what he says and no more.** Holding
   $C_{\rm A}$ constant by lowering the flow rate gives you $d$ cleanly from
   the shape of $\tau'(t)$ — that is real, and eqs. (28)–(32) reproduce
   exactly. It does *not* give you the class: within one run parallel, series
   and independent are the same equation. Recovering the class needs a **set**
   of runs at different held levels, and the span matters: with $k'_d$ known
   to 5 % you need a span of order 0.1–0.4 in $C_{\rm A}/C_{\rm A0}$, and
   *where the set starts matters as much as how wide it is* — the swept range
   over plausible starting levels and level counts is a factor of five. A
   narrow bracket of levels tells you nothing.
5. **Run a bed, and assay the bed.** The gradient is the information — but it
   is in the *solid*, not in the effluent. Judged on conversion history alone
   the bed separates one class pair in twelve at a 0.2-percentage-point
   resolution; judged on the spent catalyst it separates parallel from series
   at a glance, because the profile's gradient reverses sign — and the inlet
   activity is exactly $e^{-\lambda\theta}$ under parallel decay against
   exactly 1 under series. One sample beats a whole campaign of outlet
   analyses. Before designing on it, size the contrast your assay can actually
   reach: `B2.2` measured its own *within-bed* signature below a 10 % detection
   floor at its operating point, and that is the same family of measurement.
6. **Hold $C_{\rm P}$, not $C_{\rm A}$, if the poison is consumed.** That is
   the paper's own instruction for eq. (6) and this page confirms it: with the
   poison in excess the deactivation order comes back exactly. Holding
   $C_{\rm A}$ instead — one knob, two constraints — bends the Fig.-6 line and
   biases $d$, except at $d = 1$ where the $a^{d-1}$ consumption rate makes the
   drift cancel exactly. The paper hedges the failure of the $C_{\rm A}$
   condition, not of the $C_{\rm P}$ one; this page sizes what that hedge
   costs.

**Structure reuse.** The `Bed` class is a one-dimensional convection problem
with a distributed, slowly-evolving source parameter: identical machinery to
[`B2.2`](../B2.2-froment-bischoff-coking/) (coke profile), and the same
`construct_convflux_upwind` + `construct_div` + `interp_cntr_to_stagg_tvd`
stack as `J1.5` (breakthrough) and `A2.1` (dispersion). The activity is
marched in $w = (a^{1-d}-1)/(1-d)$ — the paper's own linearising substitution,
eqs. (28)–(32) — which is worth copying: it removes $a^d$ from every
right-hand side and makes the finite-time extinction of $d < 1$ ordinary
rather than stiff.

**Where this sits in the ladder.** [`B2.1`](../B2.1-voorhies-coking-law/)
(coke keeps time), [`B2.2`](../B2.2-froment-bischoff-coking/) (coke keeps
place), this page (what an experiment can actually recover). `B2.4`
(Beeckman–Froment) is the next rung.

**Cite the source, not this page:** Levenspiel, O., *Experimental Search for a
Simple Rate Equation to Describe Deactivating Porous Catalyst Particles*,
Journal of Catalysis **25**(2) 265–272 (1972),
[doi:10.1016/0021-9517(72)90227-8](https://doi.org/10.1016/0021-9517(72)90227-8).
Note that this is Levenspiel the author; *Chemical Reaction Engineering* is a
different work and is not reproduced anywhere here.

**Fit/test, one line:** nothing on this page is fitted to any measurement and
nothing is validated against one — the paper reports no data. Every number is
exact reproduction of printed algebra, an exact algebraic identity between the
printed rate forms, or a computed property of those rate forms in a reactor."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python",
                   "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
