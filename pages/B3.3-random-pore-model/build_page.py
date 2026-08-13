#!/usr/bin/env python3
"""Generate index.ipynb for page B3.3 (Bhatia & Perlmutter 1980, random pore model, Part I).

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------------- title
cells.append(md(r"""---
title: "The random pore model, 1980: eight figures, no table — and what the printed numbers still decide"
description: "Bhatia and Perlmutter build the reacting solid out of randomly overlapping cylindrical pores, and one parameter psi = 4 pi L_o (1 - eps_o)/S_o^2 decides whether the reaction rate rises before it falls. The paper contains NO TABLE: its only comparison with experiment is Figure 8, so this page does not digitise it and does not establish empirical adequacy. What it does instead is hold the paper to what it prints. Fifteen symbolic identities close the derivation chain from eqs. (17)-(30) to (31)-(34) and (37); the printed structural claims become numbers (the best-fit grain shape factor is 1.0000000 at psi = 0 and 0.4905224 at psi = 2 against the printed 0.49 <= m <= 1, and 0.6638850 at psi = 1 against the printed 'very close to two-thirds'); and the two special-case reductions the abstract claims are checked — Bhatia's eq. (40) is the published B3.1's shrinking core in its reaction-control limit to 1.1e-16, and the approach to Petersen (1957) is quantified at 2.2 conversion points up to X = 0.75, the printed claim landing at X = 0.7423. Three printed defects: 'This match at m = 1' where the paper's own Figure 6 and this page's least squares both give m = 0.664 at psi = 1; the optimal-porosity worked example whose printed answer 0.1 is not a root of the paper's own eq. (35) (the root is 0.08529); and Park & Levenspiel dated 1976 in the text against 1975 in the reference list."
categories: [sec:B, struct:S1, struct:S3, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-13
---

# The random pore model, 1980: eight figures, no table — and what the printed numbers still decide

**Catalog ID:** `B3.3` · **Structures:** `S1` (pointwise algebra), `S3`
(nested-scale structure–reaction coupling) · **Tier:** T0

Every gas–solid reaction model before this one had to *choose* whether the
reaction rate falls monotonically. Grain models and order-of-reaction models
say it must; the char-gasification literature kept reporting that it does not.
Bhatia and Perlmutter's move is to stop describing the solid and start
describing its **pores**: a set of cylindrical surfaces of arbitrary size
distribution $f(r)$ that grow at the chemical rate and **overlap each other**
as they grow. Growth adds surface; overlap destroys it. One dimensionless
group decides which wins,

$$\psi \;=\; \frac{4\pi L_o (1-\epsilon_o)}{S_o^{2}},$$

and the whole model is two closed forms in it — eqs. (31) and (32) below.
A rate maximum exists **iff** $\psi \ge 2$, and where it sits is eq. (33).

**What this page can and cannot check, decided before any code was written.**
The paper carries **eight figures and no table** (§4 gives the search). Its
only comparison with experiment is **Figure 8** — the char data of Hashimoto
et al. (1979), which is not on disk here. That figure is **not digitised**, no
curve, point or axis coordinate is read from it or from any other figure, and
so **this page does not establish that the random pore model describes any
real system.** It establishes what the model *is*, and holds the paper to the
numbers it prints in prose and in figure-frame captions — which, unlike the
comparable `J1.4`, turns out to be a great deal:

1. **The derivation chain, symbolically.** Fifteen `sympy` identities, all
   exactly zero: eqs. (17)–(20), (23), (26) collapsing to (27); (28)+(30)
   giving (32); (31) differentiated giving (33) — including the observation
   that eq. (37)'s radical is the *perfect square* $(1+\psi\tau/2)^2$, so
   $S/S_o$ **is** $dX/d\tau$; Petersen's eqs. (3)–(4) giving (34) through the
   fact that his surface maximum always sits at $\epsilon = 1/2$; his eq. (7)
   being exactly the statement $G = 3/(2u_o)$; and eq. (35)'s two endpoints
   being precisely the two sides of inequality (36).
2. **The printed structural claims, turned into numbers.** The best-fit grain
   shape factor from the paper's own squared-deviation integral $I$ is
   **1.0000000 at $\psi = 0$ and 0.4905224 at $\psi = 2$** (printed:
   $0.49 \le m \le 1$), and **0.6638850 at $\psi = 1$** (printed: "very close
   to two-thirds" — it is 0.42 % below). $X_M \to 1-e^{-1/2} = 0.3934693$
   (printed: $< 0.393$). $\tau_{1/2}$ at $\sigma = 100$ is within **2.9 %** of
   its $\sigma \to \infty$ value (printed: "no further significant change for
   $\sigma \ge 100$"), and at $\sigma = 0.25$ it spreads only **7.2 %** across
   $0 \le \psi \le 100$ (printed: "independent of internal structure").
3. **The two special-case reductions the abstract claims.** Bhatia's eq. (40)
   — the $m = 2/3$ grain model integrated — **is** the published `B3.1`'s
   Yagi–Kunii eq. (6) in its reaction-control limit, to **1.1e-16**, two
   transcriptions of two different papers made three months apart. And the
   approach to Petersen (1957) is quantified rather than asserted: the two
   conversion–time curves stay within **0.0218 in $X$ up to $X = 0.75$** and
   first part by 0.02 at $X = 0.7423$ — the paper's "good agreement up to
   about 75 % conversion", to three digits, **under the convention p. 384
   prints for curve B** ("If $\epsilon_o = 0.3$ and $L_o = 3.14\times10^6$ are
   chosen, the results are the same as the prior curve B"); §6.3 quotes that
   sentence and prints the whole band: the over-determined triple has *three*
   pairings, and the other two part at 0.7599 and 0.7687, so the dependence is
   at the second digit. The surfaces part much earlier — −21.7 % at $X = 0.75$ on the
   Petersen-consistent Fig. 4 set, −21.0 % on curves A and B's own set — and
   the structural reason is printed here: the two
   models use **different overlap laws**, $1-e^{-V_E}$ against
   $3u^2-2u^3$, which differ at order $V_E^{3/2}$ and can never coincide.
4. **Three printed defects, reported and not repaired.** The sentence "This
   match at $m = 1$" [*sic*, p. 383] where the paper's own Figure 6, its own
   preceding sentence, and this page's least squares all say $m \approx 2/3$
   at $\psi = 1$; the optimal-porosity worked example whose printed answer
   $\epsilon_o = 0.1$ **is not a root of the paper's own eq. (35)** (the root
   is 0.08529, residual −4.2e-3 — though the surface it costs is 0.003 %);
   and "Park and Levenspiel (1976)" in the text against "(1975)" in the
   reference list.
5. **One thing the paper's own numbers say that it does not.** Curve C of
   Figure 7 is described only as "no longer consistent with the requirements
   of Petersen's model". Its printed $(S_o, L_o)$ pair is **infeasible**, not
   merely over-determined: no Petersen structure of length $3.14\times10^6$
   cm/cm³ can have a surface above **2720 cm²/cm³**, and curve C asks for
   12 500 — a factor **4.60** beyond the maximum. That is why the paper had to
   free a parameter, and it also reconstructs the curve-D length the paper
   says it computed internally but never prints: $7.747\times10^7$ cm/cm³.

**Nothing on this page is a fit to data and nothing on it is experimental.**
The one least-squares operation is the paper's own $I$, fitted between two
*models*. Read §4 before quoting anything here as evidence about a real char.
"""))

# ---------------------------------------------------------------- background
cells.append(md(r"""## Background

### The problem in 1980

A non-catalytic gas–solid reaction $aA(g) + bB(s) \to pP(g) + qQ(s)$ (eq. 1)
runs on internal surface, so a model must say how much surface there is at
each conversion. The paper opens by splitting the field in two: models that
put the reaction on the surfaces of non-porous **grains**, and models that put
it on **pore** surfaces.

The grain family reduces to one rate law,

$$\frac{dX}{dt} = \frac{k_s C^n S_o}{1-\epsilon_o}(1-X)^m \tag{2}$$

with $m$ "a shape factor that depends on the geometry of the grains: for
spheres $m = 2/3$, for cylinders $m = 1/2$ and for flat plates $m = 0$", and
$m = 1$ recovering the volume-reaction model of Lacey et al. (1965) and Ishida
and Wen (1971). Its problem is stated plainly: "it predicts a monotonically
decreasing reaction rate because the reacting surface of each grain is
receding", while gasification reactions (Dutta et al. 1977; Dutta and Wen
1977) "exhibit a maximum in the reaction rate", supported by the surface-area
measurements of Kawahata and Walker (1962) and Hashimoto et al. (1979).

The pore family's representative is **Petersen (1957)**, whose results the
paper reprints as its eqs. (3)–(5) — uniform cylindrical pores of radius $r_p$,
total length $L$, with an intersection correction

$$S = 2\pi r_p L\left(1 - r_p\sqrt{\pi L/3}\right) \tag{3}$$
$$\epsilon = \pi r_p^2 L\left(1 - \tfrac{2 r_p}{3}\sqrt{\pi L/3}\right) \tag{4}$$
$$\frac{dr_p}{dt} = k_s C^n \tag{5}$$

and, through Szekely et al. (1976), a conversion–time relation (6) with a
cubic (7) for its constant $G$. Petersen "makes no provision for further pore
wall intersections and neglects the distribution of pore sizes" — the two
things this paper puts back.

### The move

Reaction surfaces are taken to be a set of cylindrical surfaces of
distribution $f(r)$ that **overlap at random as they grow**, and the
book-keeping is Avrami's (1940), imported from crystal-aggregate kinetics.
Two systems are tracked: the *non-overlapped* one, whose length, surface and
volume $(L_E, S_E, V_E)$ obey a population balance (11) and grow at the
chemical rate (12), and the *actual* (overlapped) one $(L, S, V)$ reached
through Avrami's statistical result $dV = (1-V)\,dV_E$ (19). That single
relation is the whole physical content of the overlap: integrating it gives

$$V = 1 - \exp(-V_E), \qquad S = S_E(1-V), \qquad L = L_E(1-V)
\tag{20, 23, 26}$$

and no assumption about pore shape is needed anywhere.

### Where this page's siblings sit

`B3.1` is Yagi and Kunii's shrinking core — the single equation the three
classical regimes come from — and this page's §7.5 shows Bhatia's eq. (40)
*is* that equation in its reaction-control limit. `B3.2` is Szekely and Evans'
grain model with intergrain diffusion, the model whose *local* law is eq. (2)
here; it is the same $m = 2/3$ grain, embedded in a diffusion field this
paper's kinetic-control assumption removes. Both were read for this page and
neither supplies a dataset — neither has one — so no cross-page CSV is loaded;
the one number borrowed from either is `B3.1`'s own published metric, printed
beside this page's in §7.5 and reconciled there.

Part II of this paper (non-isothermal operation and diffusion control) is
**not on disk**, is not consulted, and nothing here depends on it."""))

# ---------------------------------------------------------- published model
cells.append(md(r"""## The published model

All equation numbers are the paper's, read from crops of a 300 ppi native
render. `pdfimages -list` reports all eight pages as CCITT-G4 bilevel at
300 × 300 ppi; the text layer was used only as a search index.

### The two systems, and Avrami's bridge (eqs. 8–26)

With $f(r)\,dr$ the total length per unit volume of cylindrical surfaces of
radius between $r$ and $r+dr$,

$$L_E = \int_0^\infty f\,dr, \quad
S_E = 2\pi\int_0^\infty r f\,dr, \quad
V_E = \pi\int_0^\infty r^2 f\,dr \tag{8–10}$$

and, with every surface element growing at $dr/dt = k_sC^n$ (12), the
population balance (11) integrates to $dL_E/dt = 0$ (14),
$dS_E/dt = 2\pi k_sC^n L_E$ (15), $dV_E/dt = k_sC^n S_E$ (16), hence

$$S_E = \sqrt{S_{Eo}^2 + 4\pi L_{Eo}(V_E - V_{Eo})}, \qquad
V_E = V_{Eo} + S_{Eo}k_sC^nt + \pi L_{Eo}(k_sC^nt)^2. \tag{17, 18}$$

Avrami's $dV = (1-V)dV_E$ (19) and the two surface/length arguments (21)–(25)
then give (20), (23), (26), and combining everything:

$$\frac{S}{S_o} = \left(\frac{1-V}{1-V_o}\right)
\sqrt{1 - \frac{4\pi L_o(1-V_o)}{S_o^2}\ln\!\left(\frac{1-V}{1-V_o}\right)}
\tag{27}$$

$$V = 1 - (1-V_o)\exp\!\left[-\frac{k_sC^nt}{1-V_o}\left(S_o + \pi L_o k_sC^nt\right)\right]
\tag{28}$$

### The dimensionless model (eqs. 29–33, 37)

For a spherical particle under reaction control the external surface recedes
independently,

$$1 - X = \left(\frac{1-V}{1-V_o}\right)\left(1 - \frac{k_sC^nt}{R_o}\right)^3,
\qquad V_0 = \epsilon_o \tag{29, 30}$$

and with the three groups of the NOTATION list (p. 386),

$$\boxed{\;\psi = \frac{4\pi L_o(1-\epsilon_o)}{S_o^2}, \qquad
\sigma = \frac{R_o S_o}{1-\epsilon_o}, \qquad
\tau = \frac{k_s C^n S_o t}{1-\epsilon_o}\;}$$

the model is two closed forms:

$$\frac{S}{S_o} = \frac{1-X}{\left(1-\frac{\tau}{\sigma}\right)^3}
\sqrt{1 - \psi\ln\!\left[\frac{1-X}{\left(1-\frac{\tau}{\sigma}\right)^3}\right]}
\tag{31}$$

$$X = 1 - \left(1-\frac{\tau}{\sigma}\right)^3
\exp\!\left[-\tau\left(1+\frac{\psi\tau}{4}\right)\right] \tag{32}$$

$$X_M = 1 - \exp[(2-\psi)/2\psi] \tag{33}$$

$$\frac{dX}{d\tau} = (1-X)\sqrt{1-\psi\ln(1-X)} \qquad (\sigma\to\infty) \tag{37}$$

Eq. (33) is $dS/dX = 0$; the paper states it "exists over the range
$2 \le \psi < \infty$ and gives $0 \le X_M < 0.393$". Eq. (37) reduces
"exactly to the volume reaction model if $\psi = 0$", which is eq. (39),
$dX/d\tau = (1-X)^m$, at $m = 1$; at $m = 2/3$ eq. (39) integrates to
$X = 1-(1-\tau/3)^3$ (40).

### The comparisons the paper draws

**Against the grain model.** A best $m$ is defined by minimising the printed
squared-deviation integral

$$I = \int_0^1\left[(1-X)^1\sqrt{1-\psi\ln(1-X)} - (1-X)^m\right]^2 dX$$

and the result is Figure 6, summarised in prose: "For $0 \le \psi \le 2$, the
best fit gives $0.49 \le m \le 1$. For $\psi = 1$, the best $m$ is very close
to two-thirds."

**Against Petersen.** For a rate maximum Petersen's eqs. (3)–(7) give

$$X_M = \frac{1-2\epsilon_o}{2(1-\epsilon_o)} \tag{34}$$

("a rate maximum requires $0 \le \epsilon_o \le 0.5$ and can vary between
$0 \le X_M \le 0.5$"), against the random pore model's $\psi \ge 2$,
$X_M < 0.393$.

**Optimal structure.** Maximising eq. (27)'s surface over $V_o$ at fixed $S_o$,
$L_o$ and final $\epsilon$ gives

$$\frac{S_o^2}{2\pi L_o(1-\epsilon_o)} = 1 + \ln\!\left(\frac{1-\epsilon}{1-\epsilon_o}\right)
\tag{35}$$

$$[1+\ln(1-\epsilon)] \;\le\; \frac{S_o^2}{2\pi L_o} \;\le\; [1-\epsilon] \tag{36}$$

with a worked example: for $\epsilon = 0.3$ and $S_o^2/2\pi L_o = 0.67$, "the
initial $\epsilon_o$ should be 0.1"; below 0.64 "there is no $\epsilon_o$ that
will produce a maximum area at $\epsilon = 0.3$".

**Against experiment.** As $\sigma\to\infty$ eq. (31) predicts that
$[S/(1-X)]^2$ is linear in $\ln[1/(1-X)]$ — Figure 8, whose intercept and two
slopes give the numbers quoted in §4. Conversions and surfaces are extracted
from the char data through eqs. (44)–(45), and $\epsilon_o$, $S_o$, $L_o$ are
in principle obtainable from a porosimetry distribution $v_o(r)$ by
eqs. (41)–(43)."""))

# ------------------------------------------------------------- environment
cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code(r'''import sys, urllib.request
from pathlib import Path

if not any("shared" in p for p in sys.path):
    local = Path.cwd()
    for _ in range(5):
        if (local / "shared" / "gallery_utils.py").is_file():
            sys.path.insert(0, str(local / "shared")); break
        local = local.parent
    else:
        url = ("https://raw.githubusercontent.com/computational-chemical-engineering/"
               "pymrm-gallery/main/shared/gallery_utils.py")
        urllib.request.urlretrieve(url, "gallery_utils.py")
        sys.path.insert(0, ".")

import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
from functools import lru_cache
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaincc, gamma as gamma_fn
from pymrm import newton, NumJac
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "B3.3-random-pore-model"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
pd.set_option("display.width", 135)

np.random.seed(0)   # nothing here is stochastic; seeded so it stays that way'''))

# --------------------------------------------------- symbolic derivation
cells.append(md(r"""### 2.1 The derivation chain, verified symbolically

Fifteen identities, every one an equation the paper prints or a step it
claims. Each `assert` fails the notebook if the algebra does not close. These
are **structural** checks: they verify the paper's algebra and this page's
closed forms, and they cannot detect a numerical defect — the numerics get
their own breakable tests in §7.

Two of them are worth naming in advance. The fifth shows that eq. (37)'s
radical is a *perfect square*, $1-\psi\ln(1-X) = (1+\psi\tau/2)^2$, so the
paper's $S/S_o$ **is** $dX/d\tau$ exactly — surface and rate are the same
curve, which is why Figure 4 doubles as a rate plot. The ninth shows that
Petersen's surface maximum always sits at $\epsilon = 1/2$ regardless of
$L$ and $r_p$, which is the whole content of eq. (34)."""))

cells.append(code(r'''So, Lo, ks, C, nn, tau, psi, r, L, u, m = sp.symbols(
    "S_o L_o k_s C n tau psi r L u m", positive=True)
v, vo, eps, epso, ww, Ghat, p3 = sp.symbols(
    "v v_o epsilon epsilon_o w Ghat p", positive=True)   # v = 1-V, w = 1-X
ZERO = []

# --- (17)+(20)+(23)+(26) -> (27), written in v = 1-V so the logs stay real
SEo, LEo = So/vo, Lo/vo                                  # (23),(26) at t = 0
SE = sp.sqrt(SEo**2 + 4*sp.pi*LEo*(-sp.log(v) + sp.log(vo)))     # (17)+(20)
lhs27 = (SE*v/So)**2
rhs27 = (v/vo*sp.sqrt(1 - 4*sp.pi*Lo*vo/So**2*sp.log(v/vo)))**2
ZERO.append(("(17)+(20)+(23)+(26) -> (27)", sp.simplify(sp.expand(lhs27 - rhs27))))

# --- (28)+(29)+(30) -> the exponential of (32), with tau and psi substituted
kt = tau*(1 - epso)/So                                   # k_s C^n t from tau
psi_def = 4*sp.pi*Lo*(1 - epso)/So**2
ZERO.append(("(28)+(30) -> (32)",
             sp.simplify(sp.exp(-(kt/(1-epso))*(So + sp.pi*Lo*kt))
                         - sp.exp(-tau*(1 + psi_def*tau/4)))))
ZERO.append(("(27) coefficient IS psi", sp.simplify(4*sp.pi*Lo*(1-epso)/So**2 - psi_def)))

# --- (31) differentiated at sigma -> infinity gives (33)
SS = ww*sp.sqrt(1 - psi*sp.log(ww))                      # S/S_o with w = 1-X
wM = sp.solve(sp.Eq(sp.diff(SS, ww), 0), ww)[0]
ZERO.append(("(31) -> (33) X_M", sp.simplify((1 - wM) - (1 - sp.exp((2-psi)/(2*psi))))))

# --- (32) differentiated IS (37), and the radical is a perfect square
Xe = 1 - sp.exp(-tau*(1 + psi*tau/4))
ZERO.append(("d(32)/dtau = (1-X)(1+psi tau/2)", sp.simplify(sp.diff(Xe, tau) - (1-Xe)*(1+psi*tau/2))))
ZERO.append(("(37) radical is the perfect square (1+psi tau/2)^2",
             sp.simplify((1 - psi*sp.log(1-Xe)) - (1 + psi*tau/2)**2)))

# --- (39) integrated: general m closed form, the m = 2/3 case (40), the m = 1 branch
q = sp.Symbol("q", positive=True)                        # q = 1-(1-m)tau
X39 = 1 - q**(1/(1-m))
ZERO.append(("(39) general m closed form",
             sp.simplify(sp.diff(X39.subs(q, 1-(1-m)*tau), tau) - (1-(1-m)*tau)**(m/(1-m)))))
X40 = 1 - p3**3                                          # p = 1-tau/3
ZERO.append(("(39) at m = 2/3 -> (40)",
             sp.simplify(sp.diff(X40.subs(p3, 1-tau/3), tau) - (p3**2).subs(p3, 1-tau/3))))
ZERO.append(("(39) m -> 1 branch limit is 1-exp(-tau)",
             sp.simplify(sp.limit(1 - (1-(1-m)*tau)**(1/(1-m)), m, 1) - (1 - sp.exp(-tau)))))

# --- Petersen (3),(4): the surface maximum always sits at eps = 1/2  -> (34)
bb = sp.sqrt(sp.pi*L/3)
Spet, epet = 2*sp.pi*r*L*(1 - r*bb), sp.pi*r**2*L*(1 - 2*r*bb/3)
rstar = sp.solve(sp.Eq(sp.diff(Spet, r), 0), r)[0]
ZERO.append(("Petersen (3),(4): S is maximal exactly at eps = 1/2",
             sp.simplify(epet.subs(r, rstar) - sp.Rational(1, 2))))
ZERO.append(("-> (34) X_M = (1-2 eps_o)/(2(1-eps_o))",
             sp.simplify((sp.Rational(1,2) - epso)/(1-epso) - (1-2*epso)/(2*(1-epso)))))

# --- Petersen (7) IS the statement G = 3/(2 u_o) with eps_o = 3u^2 - 2u^3
ZERO.append(("(7) <=> G = 3/(2 u_o)",
             sp.simplify(sp.Rational(4,27)*(3*u**2 - 2*u**3)*(3/(2*u))**3 - 3/(2*u) + 1)))

# --- (35) at its two endpoints IS inequality (36)
st35 = Ghat/(1-epso) - 1 - sp.log((1-eps)/(1-epso))
ZERO.append(("(35) at eps_o -> 0 is (36)'s lower bound",
             sp.simplify(sp.solve(st35.subs(epso, 0), Ghat)[0] - (1 + sp.log(1-eps)))))
ZERO.append(("(35) at eps_o -> eps is (36)'s upper bound",
             sp.simplify(sp.solve(st35.subs(epso, eps), Ghat)[0] - (1 - eps))))

# --- (31) -> the Figure 8 linearisation
ZERO.append(("(31) -> [S/(1-X)]^2 = S_o^2 (1 + psi ln[1/(1-X)])",
             sp.simplify((So*SS/ww)**2 - So**2*(1 + psi*sp.log(1/ww)))))

for name, res in ZERO:
    print(f"  0 == {res}   [{name}]")
    assert sp.simplify(res) == 0, name
SYM_ZERO_COUNT = len(ZERO)
print(f"\n{SYM_ZERO_COUNT} symbolic identities, all exactly zero.")'''))

# ------------------------------------------------- parameters and assumptions
cells.append(md(r"""## Parameters and assumptions

**The paper's assumptions, inherited unchanged and stated where they bite.**

- **Kinetic control.** "If all diffusional resistances are assumed negligible,
  no distinction need be made between macro and micropores (the kinetic
  regime)." Everything on this page lives in that regime; Part II lifts it and
  is not on disk. The paper's own justification is the effectiveness-factor
  sentence quoted in §4, whose inputs it does not print — §6.6 says exactly
  what that claim can and cannot be checked against here.
- **No closed pore volume** in the reactant. The paper uses this to *exclude*
  Kawahata and Walker's anthracite from the comparison ("their reactant solid
  had a large closed pore volume, over 60 % of the open pore volume").
- **Reaction rate proportional to the actual surface**, $dr/dt = k_sC^n$ on
  every non-overlapped element (12) — the step that makes the population
  balance closed.
- **Avrami's $dV = (1-V)dV_E$** (19) is "a statistical average" and the paper
  flags where it is weak: Bhatia and Perlmutter (1979) "have pointed out that
  this relationship … may not be sufficiently accurate in the initial stages
  of a solid phase reaction when the number of nuclei is small". Here the
  nuclei are pores, and the paper argues the count is large.
- **Spherical particle, receding external surface** for the finite-$\sigma$
  results (eq. 29). $\sigma\to\infty$ removes it, and the paper works there
  for Figures 3, 4, 6, 7 and 8.

**Parameter values.** Every number this page consumes comes from
`bhatia-1980-printed-scalars.csv` (§4) and is interpolated, never retyped. The
three parameter sets that carry the Petersen comparison are, as printed:

| set | $S_o$ (cm²/cm³) | $L_o$ (cm/cm³) | $\epsilon_o$ | where |
| --- | --- | --- | --- | --- |
| Petersen-consistent | 2 425 | $3.14\times10^6$ | 0.26 | Fig. 4 caption **and** p. 384 text |
| Fig. 7 curves A, B | 2 500 | $3.14\times10^6$ | 0.30 | Fig. 7 caption |
| Fig. 7 curve C | 12 500 | $3.14\times10^6$ | 0.30 | Fig. 7 caption **and** p. 384 text |

**Units.** All lengths in cm, so $S_o$ in cm²/cm³ and $L_o$ in cm/cm³ as
printed; $\psi$, $\sigma$, $\tau$, $X$ and $S/S_o$ are dimensionless. The
Hashimoto surface area is printed in m²/cm³ and is converted once, in code,
by $10^4$ cm²/m²."""))

# ------------------------------------------------------------------ the data
cells.append(md(r"""## The data

**There is none, and that is this page's central scope fact.** Stated
precisely, with the search behind each claim. All eight pages were rendered at
their native 300 ppi (`pdfimages -list`: CCITT-G4 bilevel, 300 × 300 ppi on
every page) and read column-block by column-block; the text layer was used as
a search index only and for no digit.

- **The paper contains no table at all.** Its eight pages carry eight
  **figures** (Fig. 1 overlapping cylinders; Fig. 2 $\tau_{1/2}$ vs $\sigma$;
  Fig. 3 $X$ vs $\tau$; Fig. 4 $S/S_o$ vs $X$; Fig. 5 optimal $\epsilon_o$ vs
  the structural group; Fig. 6 best-fit $m$ vs $\psi$; Fig. 7 the Petersen
  conversion–time comparison; Fig. 8 the Hashimoto correlation), a NOTATION
  list and a LITERATURE CITED list, and no other tabular matter. Second,
  independent search: the extracted text layer grepped case-insensitively for
  `table` returns two lines, one inside the word "suitable" (p. 385) and one
  inside a scan-garbled "surfaces" (p. 380) — never a caption, never a
  cross-reference.
- **The one comparison with experiment is Figure 8**, the CS and VC char data
  of Hashimoto et al. (1979) — a paper that is **not on disk here**. It is
  **not digitised**: no point, no curve and no axis coordinate is taken from
  it, or from any other figure. Figure content needs a human review gate this
  build did not have, so the comparison is out of scope, and **this page
  therefore does not establish that the random pore model describes any real
  system.**
- **But the paper prints its fit results in prose**, which is what separates
  this case from the figure-only ones. P. 385, verbatim: "Estimates of initial
  BET surface areas obtained from the intercept of Figure 8 give 520 m²/cm³
  for both VC and CS chars. The two slopes give $\psi = 6.9$ and $\psi = 13.7$
  for the VC and CS chars, respectively." Those three numbers are transcribed
  from the prose sentence and used in §6.6; nothing else from Figure 8 is used.
- **Figure-frame captions are used for parameter *sets* only** — Fig. 4's and
  Fig. 7's inputs, which are typeset text stating what was computed, the same
  treatment `A3.7` and `J1.4` give figure legends. Both sets are also printed
  in the p. 384 running text, and were read in both places.
- **Nothing is fitted to data anywhere on this page.** The single
  least-squares operation is the paper's own integral $I$, which fits one
  *model* to another.

**No other page's dataset is loaded**, so the cross-page reconciliation
obligations bite only once: §7.5 borrows `B3.1`'s **published metric** for the
reaction-control limit and prints it beside this page's own measurement rather
than retyping it as a claim."""))

cells.append(code(r'''scal = load_data("bhatia-1980-printed-scalars.csv", page=PAGE)
P = dict(zip(scal.key, pd.to_numeric(scal.printed)))
print(cite_data(load_meta("bhatia-1980-printed-scalars.csv", page=PAGE)))
print(f"{len(scal)} printed scalars transcribed; the paper prints "
      f"{P['n_figures']:.0f} figures and {P['n_tables']:.0f} tables.\n")

# --- the paper's own numbers, cross-checked against each other --------------
assert P["fig4_Lo"] == P["fig7_AB_Lo"] == P["fig7_C_Lo"], "the three sets share L_o"
assert P["fig7_agreement_X"] == P["conclusions_agreement_X"], \
    "the 75 % claim is printed twice (p. 384 and the CONCLUSIONS); they must agree"
assert P["m_range_hi"] == P["grain_m_volume"], \
    "the printed upper end of the best-fit m range is the volume-reaction model's m = 1"
assert P["XM_psi_min"] == P["psi_ratemax"], \
    "the psi threshold for a rate maximum is printed twice (pp. 382, 384)"
print("printed-value cross-checks (all three sets share L_o = "
      f"{P['fig4_Lo']:.3g} cm/cm3; the 75 % claim appears on p. 384 AND p. 379;")
print(" the psi = 2 threshold appears on p. 382 AND p. 384): consistent.")

# --- a printed bibliographic defect, proved from the paper's own two mentions
print("\nPRINTED DEFECT 1 (bibliographic). P. 380 left column: 'A variation of the spherical")
print("grain model has also been suggested by Park and Levenspiel (1976)' [sic]. The paper's")
print("LITERATURE CITED (p. 385) reads: 'Park, J. Y., and O. Levenspiel, \"The Crackling Core")
print("Model for the Reaction of Solid Particles,\" Chem. Eng. Sci., 30, 1207 (1975).' The text")
print("year and the reference year differ by one; it is the paper's ONLY Park reference.")
print("Reported, not repaired: which of the two is right cannot be settled from this document.")'''))

# --------------------------------------------------------- implementation
cells.append(md(r"""## PyMRM implementation

The random pore model is closed-form algebra, so what needs solving is a
**family of scalar root problems** — inversions of eqs. (32), (35), (3), (4)
and $dI/dm = 0$ — plus one initial-value march. Both are done the way the
gallery does blocks of independent scalars: `NumJac((K, 1))` with `newton`,
the shape `(K, 1)` and never `(K,)`, so the Jacobian is diagonal by
construction rather than a dense $K\times K$ (the house rule, measured on
`B1.1`).

Two implementation choices matter and are load-bearing:

1. **Eq. (32) is inverted in log form.** The residual solved is
   $3\ln(1-\tau/\sigma) - \tau(1+\psi\tau/4) - \ln(1-X^\star) = 0$, not
   $X(\tau) - X^\star = 0$. At $\psi = 100$ the exponential in eq. (32)
   underflows and the direct residual's Jacobian row goes exactly zero —
   `newton` returns a singular matrix. The log form is the same equation with
   no exponential in it.
2. **$\tau$ is carried through a sigmoid**, $\tau = \sigma/(1+e^{-w})$, so the
   iterate cannot leave $(0,\sigma)$ where eq. (32) is defined. For
   $\sigma\to\infty$ the half-conversion time has a closed form
   ($\psi\tau^2/4 + \tau = \ln 2$) and the block solver is not needed at all —
   §7.4 uses that closed form as an independent check on the block solver.

Everything else is direct evaluation of the printed equations."""))

cells.append(code(r'''# ============ the random pore model, exactly as printed ====================
def psi_of(So_, Lo_, eps_):                     # NOTATION list, p. 386
    return 4*np.pi*Lo_*(1 - eps_)/So_**2

def X_of_tau(tau_, psi_, sigma=np.inf):         # eq. (32)
    f = 1.0 if not np.isfinite(sigma) else (1 - tau_/sigma)**3
    return 1 - f*np.exp(-tau_*(1 + psi_*tau_/4))

def SS0_of_X(X_, psi_):                         # eq. (31) at sigma -> infinity
    w = 1 - X_
    return w*np.sqrt(1 - psi_*np.log(w))

def rate37(X_, psi_):                           # eq. (37); == SS0_of_X, proved in 2.1
    w = np.maximum(1 - X_, 1e-300)
    return w*np.sqrt(1 - psi_*np.log(w))

def XM_eq33(psi_):                              # eq. (33)
    return 1 - np.exp((2 - psi_)/(2*psi_))

def X_grain(tau_, m_):                          # eq. (39) integrated; (40) is m = 2/3
    return 1 - np.exp(-tau_) if m_ == 1 else 1 - (1 - (1-m_)*tau_)**(1/(1-m_))

# ============ Petersen (1957) as the paper reprints it =====================
def pet_S(r_, L_):  return 2*np.pi*r_*L_*(1 - r_*np.sqrt(np.pi*L_/3))     # eq. (3)
def pet_eps(r_, L_): return np.pi*r_**2*L_*(1 - (2*r_/3)*np.sqrt(np.pi*L_/3))  # eq. (4)

def u_of_eps(e_):
    """eq. (4) in u = r sqrt(pi L/3): eps = 3u^2 - 2u^3, independent of L."""
    return brentq(lambda z: 3*z**2 - 2*z**3 - e_, 0.0, 1.0, xtol=1e-16, rtol=8.9e-16)

def G_of_eps(e_):
    """The root of eq. (7) that eqs. (3),(4),(6) select: G = 3/(2 u_o).  See 6.2."""
    return 1.5/u_of_eps(e_)

def X_petersen(s_, e_):
    """eq. (6) with G from eq. (7); s = k_s C^n t / r_po."""
    G = G_of_eps(e_)
    return e_/(1-e_)*((1+s_)**2*(G - 1 - s_)/(G - 1) - 1)

def SS0_petersen(X_, e_):
    """S/S_o along Petersen's history, from eq. (3): S ~ u(1-u)."""
    u0 = u_of_eps(e_)
    uu = np.array([u_of_eps(e_ + x*(1-e_)) for x in np.atleast_1d(X_)])
    return (uu*(1-uu)/(u0*(1-u0))).reshape(np.shape(X_))

print("printed equations implemented: (31), (32), (33), (37), (39)/(40) and Petersen (3)-(7)")'''))

cells.append(code(r'''# ============ block root-finds, pymrm =====================================
def tau_at_X_block(Xt, psis_, sigmas_, tol=1e-13):
    """Invert eq. (32) for tau, for a whole block of (psi, sigma) at once.

    Residual in LOG form so nothing underflows; tau = sigma/(1+exp(-w)) keeps
    the iterate inside (0, sigma).  NumJac((K,1)) -> diagonal Jacobian."""
    ps = np.atleast_1d(np.asarray(psis_, float)).ravel()
    sg = np.broadcast_to(np.asarray(sigmas_, float), ps.shape).ravel()
    xt = np.broadcast_to(np.asarray(Xt, float), ps.shape).ravel()
    K = ps.size
    jac = NumJac((K, 1))                        # K independent scalars: (K,1), NOT (K,)
    def res(w):
        g = 1.0/(1.0 + np.exp(-w[:, 0]))        # tau/sigma in (0,1)
        t = sg*g
        return (3*np.log1p(-g) - t*(1 + ps*t/4) - np.log1p(-xt))[:, None]
    sol = newton(lambda w: jac(res, w), np.zeros((K, 1)), tol=tol, maxfev=300)
    assert sol.success, "tau block Newton did not converge"
    return sg/(1.0 + np.exp(-sol.x[:, 0]))

def eps_o_opt_block(eps_final, Ghat_, tol=1e-14):
    """Solve eq. (35) for the optimal initial porosity, as a block.

    eps_o = eps * sigmoid(w) keeps 0 < eps_o < eps (the physical branch; the
    second root of (35) lies above eps and is a MINIMUM -- see 6.5)."""
    ef = np.atleast_1d(np.asarray(eps_final, float)).ravel()
    gh = np.broadcast_to(np.asarray(Ghat_, float), ef.shape).ravel()
    K = ef.size
    jac = NumJac((K, 1))
    def res(w):
        e0 = ef/(1.0 + np.exp(-w[:, 0]))
        return (gh/(1-e0) - 1 - np.log((1-ef)/(1-e0)))[:, None]
    sol = newton(lambda w: jac(res, w), -2.0*np.ones((K, 1)), tol=tol, maxfev=300)
    assert sol.success, "eq. (35) block Newton did not converge"
    return ef/(1.0 + np.exp(-sol.x[:, 0]))

def march37(psis_, tau_end, nsteps, theta=0.5):
    """Crank-Nicolson march of eq. (37) -- the paper's RATE law -- vectorised
    over psi.  Shares nothing with the closed form (32) it is checked against."""
    ps = np.asarray(psis_, float).ravel(); K = ps.size
    dt = tau_end/nsteps
    Xc = np.zeros((K, 1)); jac = NumJac((K, 1))
    for _ in range(nsteps):
        Xo = Xc.copy()
        def res(Xn):
            return (Xn - Xo - dt*(theta*rate37(Xn[:, 0], ps)
                                  + (1-theta)*rate37(Xo[:, 0], ps))[:, None])
        sol = newton(lambda z: jac(res, z), Xc.copy(), tol=1e-14, maxfev=100)
        assert sol.success, "march37 Newton did not converge"
        Xc = sol.x.copy()
    return Xc[:, 0]

demo = tau_at_X_block(0.5, [0.0, 1.0, 2.0, 10.0, 100.0], 100.0)
print("tau_1/2 at sigma = 100 for psi = 0, 1, 2, 10, 100 (the psi labels of Fig. 2's curves):")
print("  " + "  ".join(f"{t:.6f}" for t in demo))'''))

# -------------------------------------------------------------------- results
cells.append(md(r"""## Results

### 6.1 The model: conversion, surface, and the rate maximum

The two closed forms, and the grain model they are compared against. The
$\psi$ values shown are the labels printed on **Figure 4's** own curves;
Figures 2 and 3 carry a different label set, which is the one §7.4 uses. No
coordinate is read from any of the three figures — only the curve labels, which
are typeset text."""))

cells.append(code(r'''PSI_SHOW = [0.0, 1.0, 2.0, 5.0, 10.0]   # the labels printed on FIGURE 4's curves
tg = np.linspace(0, 2.2, 601)

fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.2))
for pv in PSI_SHOW:
    ax[0].plot(tg, X_of_tau(tg, pv), lw=1.8, label=fr"$\psi$ = {pv:g}")
tg40 = np.linspace(0, 3, 301)
ax[0].plot(tg40, X_grain(tg40, 2/3), "k--", lw=1.5, label="grain model, m = 2/3  (eq. 40)")
ax[0].plot(tg, X_grain(tg, 1.0), "k:", lw=1.5, label="volume model, m = 1")
ax[0].set(xlabel=r"dimensionless time $\tau$", ylabel="conversion $X$",
          xlim=(0, 2.2), ylim=(0, 1), title="eq. (32), $\\sigma\\to\\infty$")
ax[0].legend(fontsize=8, loc="lower right")

Xg = np.linspace(1e-6, 1-1e-9, 1200)
for pv in PSI_SHOW:
    ax[1].plot(Xg, SS0_of_X(Xg, pv), lw=1.8, label=fr"$\psi$ = {pv:g}")
    if pv >= 2:
        xm = XM_eq33(pv)
        ax[1].plot([xm], [SS0_of_X(xm, pv)], "o", ms=4.5, color="k", zorder=5)
ax[1].plot(Xg, (1-Xg)**(2/3), "k--", lw=1.5, label="grain model, m = 2/3")
ax[1].set(xlabel="conversion $X$", ylabel="$S/S_o$   ( = $dX/d\\tau$, proved in 2.1)",
          xlim=(0, 1), ylim=(0, 1.65), title="eq. (31); dots are eq. (33)")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()

XM_SUP = 1 - np.exp(-0.5)
print(f"eq. (33) supremum as psi -> infinity: 1 - exp(-1/2) = {XM_SUP:.7f}"
      f"   (printed: X_M < {P['XM_sup']})")
print(f"eq. (33) at the printed threshold psi = {P['XM_psi_min']:.0f}: X_M = {XM_eq33(P['XM_psi_min']):.1e}"
      f"   -- the maximum is born exactly at the particle's start, which is what")
print("   'this maximum exists over the range 2 <= psi < infinity' means.")
print(f"\nPetersen's rival, eq. (34), at the two printed porosities:")
for e_ in (P["fig4_eps_o"], P["fig7_eps_o"]):
    print(f"   eps_o = {e_}: X_M = (1-2 eps_o)/(2(1-eps_o)) = {(1-2*e_)/(2*(1-e_)):.6f}")
print(f"   printed bound: a rate maximum requires eps_o <= {P['eq34_eps_o_max']} and gives"
      f" X_M <= {P['eq34_XM_max']} -- and eq. (34) at eps_o = 0 gives {(1-0)/(2*(1-0)):.1f}, the bound exactly.")'''))

cells.append(md(r"""### 6.2 The Petersen reduction, and a parameter set that is not merely over-determined

The paper states three parameter triples and claims two things about them: the
Fig. 4 set "satisf[ies] Equations (3) and (4)" and "give[s] $\psi = 5$", and
curve C's set is "no longer consistent with the requirements of Petersen's
model", whose "corresponding Petersen choices are not, in this case, uniquely
defined".

Both claims are checkable arithmetic, and the second is stronger than the
paper says. Eq. (3) has a **maximum** over $r_p$ at $r_p\sqrt{\pi L/3} = 1/2$
— the same stationary point that §2.1 showed is Petersen's surface maximum —
so for a given $L$ no Petersen structure can have a surface above
$S_{\max}(L)$. Curve C asks for one that does.

**Which root of eq. (7)?** The paper says "$G$ is the solution to the cubic
equation", singular. It has three real roots. §2.1 showed eq. (7) is exactly
$G = 3/(2u_o)$ with $\epsilon_o = 3u_o^2-2u_o^3$, which selects one of them;
the others are reported below with what they would do."""))

cells.append(code(r'''def petersen_audit(So_, Lo_, eps_, tag):
    b = np.sqrt(np.pi*Lo_/3); rmax = 1/(2*b); Smax = pet_S(rmax, Lo_)
    out = dict(tag=tag, So=So_, Lo=Lo_, eps=eps_, Smax=Smax, feasible=So_ <= Smax,
               psi=psi_of(So_, Lo_, eps_), rbar=So_/(2*np.pi*Lo_))
    if out["feasible"]:                       # route 1: eq. (3) -> r_p -> eq. (4)
        rp = brentq(lambda z: pet_S(z, Lo_) - So_, 1e-12, rmax, xtol=1e-18, rtol=8.9e-16)
        out["rp_from_S"] = rp
        out["eps_pred"] = pet_eps(rp, Lo_)
        out["eps_rel"] = abs(out["eps_pred"] - eps_)/eps_
    else:
        out["rp_from_S"] = np.nan; out["eps_pred"] = np.nan; out["eps_rel"] = np.nan
    r2 = u_of_eps(eps_)/b                     # route 2: eq. (4) -> r_p -> eq. (3)
    out["rp_from_eps"] = r2; out["S_pred"] = pet_S(r2, Lo_)
    out["S_rel"] = abs(out["S_pred"] - So_)/So_
    return out

SETS = [petersen_audit(P["fig4_So"],   P["fig4_Lo"],   P["fig4_eps_o"], "Fig. 4 / p. 384 Petersen set"),
        petersen_audit(P["fig7_AB_So"], P["fig7_AB_Lo"], P["fig7_eps_o"], "Fig. 7 curves A, B"),
        petersen_audit(P["fig7_C_So"],  P["fig7_C_Lo"],  P["fig7_eps_o"], "Fig. 7 curve C")]
aud = pd.DataFrame(SETS)[["tag", "So", "eps", "Smax", "feasible", "eps_pred",
                          "eps_rel", "S_pred", "S_rel", "psi", "rbar"]]
with pd.option_context("display.float_format", lambda v: f"{v:.6g}"):
    display(aud)

F4 = SETS[0]
PSI_F4      = F4["psi"]
PSI_F4_REL  = abs(PSI_F4 - P["fig4_psi"])/P["fig4_psi"]
PET_EPS_REL = F4["eps_rel"]
PET_S_REL   = F4["S_rel"]
PET_AB_EPS_REL = SETS[1]["eps_rel"]
CURVE_C_INFEAS = SETS[2]["So"]/SETS[2]["Smax"]
print(f"\nFig. 4 set: eq. (3) -> r_p = {F4['rp_from_S']:.6e} cm, and eq. (4) then returns"
      f" eps = {F4['eps_pred']:.6f} against the printed {P['fig4_eps_o']}  ({PET_EPS_REL:.4%}).")
print(f"   The other direction -- eq. (4) -> r_p -> eq. (3) -- returns S_o = {F4['S_pred']:.2f}"
      f" against the printed {P['fig4_So']:.0f}  ({PET_S_REL:.4%}).")
print(f"   psi from the same three numbers: {PSI_F4:.6f}, against the printed"
      f" psi = {P['fig4_psi']:.0f}  ({PSI_F4_REL:.3%}).  BOTH printed claims hold.")
print(f"\nCurve C: the largest surface ANY Petersen structure of length {P['fig7_C_Lo']:.3g}"
      f" cm/cm3 can have is\n   S_max = {SETS[2]['Smax']:.1f} cm2/cm3, and curve C asks for"
      f" {P['fig7_C_So']:.0f} -- a factor {CURVE_C_INFEAS:.4f} beyond it.")
print("   Eq. (3) has no root at all for that pair. 'Not uniquely defined' understates it:")
print("   the triple is INFEASIBLE, which is why a parameter had to be freed.")

# --- curve D: the length the paper says it computed internally but never prints
u0_D = u_of_eps(P["fig7_eps_o"])
L_D  = P["fig7_C_So"]**2/(4*u0_D**2*(1-u0_D)**2*3*np.pi)
rp_D = u0_D/np.sqrt(np.pi*L_D/3)
L_D_RATIO = L_D/P["fig7_AB_Lo"]
L_D_ROUNDTRIP = max(abs(pet_S(rp_D, L_D) - P["fig7_C_So"])/P["fig7_C_So"],
                    abs(pet_eps(rp_D, L_D) - P["fig7_eps_o"])/P["fig7_eps_o"])
print(f"\nCurve D ('L_o is internally calculated via Equations (3) and (4)'): fixing"
      f" eps_o = {P['fig7_eps_o']} and\n   S_o = {P['fig7_C_So']:.0f} and solving both equations gives"
      f" L = {L_D:.6e} cm/cm3 = {L_D_RATIO:.3f} x curves A/B's length,")
print(f"   with r_p = {rp_D:.5e} cm; round-trip residual {L_D_ROUNDTRIP:.2e}.")
print("   THE PAPER NEVER PRINTS THAT LENGTH. Searched: the Fig. 7 caption (which says only that")
print("   'L_o is internally calculated via Equations (3) and (4)'), both columns of p. 384 in full,")
print("   and the p. 383-385 running text either side of them, all on 300 ppi crops; the only L_o")
print("   values printed anywhere in the paper are 10^7 (p. 381, an order of magnitude) and the")
print("   3.14e6 shared by the Fig. 4 set and Fig. 7's curves A, B and C.")

# --- the roots of eq. (7), and which one eqs. (3),(4),(6) select
for e_ in (P["fig4_eps_o"], P["fig7_eps_o"]):
    rts = np.roots([(4/27)*e_, 0, -1, 1]); rts = np.sort(rts[np.isreal(rts)].real)
    Gsel = G_of_eps(e_)
    s_end = 2*Gsel/3 - 1
    print(f"\neq. (7) at eps_o = {e_}: real roots {np.array2string(rts, precision=6)}")
    print(f"   eqs. (3),(4),(6) select G = 3/(2 u_o) = {Gsel:.6f} (u_o = {u_of_eps(e_):.6f});"
          f" eq. (6) then reaches X = {X_petersen(s_end, e_):.12f} at s = {s_end:.6f}.")
    other = [g for g in rts if g > 1 and abs(g - Gsel) > 1e-6]
    for g in other:
        print(f"   the other root G = {g:.6f} would put complete gasification at s = 2G/3 - 1"
              f" = {2*g/3-1:+.6f} < 0: no such time exists.")
G_SEL_F4, G_SEL_F7 = G_of_eps(P["fig4_eps_o"]), G_of_eps(P["fig7_eps_o"])'''))

cells.append(md(r"""### 6.3 The 75 % claim, to three digits

"Curves A and B are in good agreement up to about 75 % conversion and diverge
thereafter" (p. 384), and from the other side in the CONCLUSIONS: "some
deviation between the prior model and the new one is found at high conversions
(above about 75 %)". Both models are closed forms of printed parameters, so
the claim is a number.

**Which Petersen structure is curve B? The paper says, and this page follows
it.** A Petersen solid has two degrees of freedom, $(r_p, L)$, and the Fig. 7
caption prints three numbers for curves A and B, so the triple is
over-determined — which is what §6.2's 2.5 % porosity closure measures — and a
convention has to be chosen. It is printed, on p. 384 right column, where the
paper discusses curve C's parameters:

> The corresponding Petersen choices are not, in this case, uniquely defined.
> If $\epsilon_o = 0.3$ and $L_o = 3.14 \times 10^6$ are chosen, the results
> are the same as the prior curve B. If, on the other hand, $\epsilon_o = 0.3$
> and $S_o = 12{,}500$ cm²/cm³ are fixed, curve D is obtained.

So **curve B is the Petersen structure fixed by $(\epsilon_o, L_o)$**, and
curve D the one fixed by $(\epsilon_o, S_o)$ — both by the paper's own words,
and both are what §6.2 and this section compute.

The choice is not cosmetic, and there is no single "alternative reading": an
over-determined triple $(\epsilon_o, S_o, L_o)$ admits **three** pairings, and
the page prints the whole band rather than one representative of it. The
paper's $(\epsilon_o, L_o)$ parts at $X = 0.7423$; $(\epsilon_o, S_o)$ — the
pairing p. 384 assigns to curve *D* — parts at $0.7599$; and $(S_o, L_o)$,
which keeps both printed measurements and lets the porosity be the 2.5 %-adrift
one §6.2 reports, parts at $0.7687$. So the dependence is at the *second*
digit, all three round to "about 75 %", and only the sourced convention earns
the three-digit statement."""))

cells.append(code(r'''So_A, Lo_A, e_A = P["fig7_AB_So"], P["fig7_AB_Lo"], P["fig7_eps_o"]
PSI_A  = psi_of(So_A, Lo_A, e_A)
u0_A   = u_of_eps(e_A); rpo_A = u0_A/np.sqrt(np.pi*Lo_A/3); G_A = G_of_eps(e_A)
KT_END = rpo_A*(2*G_A/3 - 1)                    # complete gasification, curve B

XA = lambda kt: X_of_tau(kt*So_A/(1-e_A), PSI_A)            # curve A, eq. (32)
XB = lambda kt: X_petersen(kt/rpo_A, e_A)                   # curve B, eq. (6)
dAB = lambda kt: XA(kt) - XB(kt)

KT_SIGN = brentq(dAB, 1e-5, KT_END*0.99, xtol=1e-20, rtol=8.9e-16)
kt_pos  = minimize_scalar(lambda z: -dAB(z), bounds=(1e-9, KT_SIGN),
                          method="bounded", options=dict(xatol=1e-18)).x
kt_neg  = minimize_scalar(lambda z: dAB(z), bounds=(KT_SIGN, KT_END),
                          method="bounded", options=dict(xatol=1e-18)).x
DAB_POS, DAB_NEG = dAB(kt_pos), dAB(kt_neg)
KT_75   = brentq(lambda z: XA(z) - P["fig7_agreement_X"], 1e-9, KT_END, xtol=1e-20, rtol=8.9e-16)
kts     = np.linspace(1e-9, KT_75, 4001)
DAB_UPTO75 = float(np.max(np.abs(dAB(kts))))
# the conversion at which the two first part by two conversion points, root-found
KT_2PT  = brentq(lambda z: abs(dAB(z)) - 0.02, KT_SIGN, kt_neg, xtol=1e-20, rtol=8.9e-16)
X_AT_2PT = XA(KT_2PT)

# The two curves cross TWICE. Petersen's initial rate is the faster of the two,
# so A trails B over the first few per cent before overtaking; both crossings are
# root-found. (No metric is affected: DAB_UPTO75 scans from 1e-9 and already
# contains the early dip.)
RATE0_A = So_A/(1 - e_A)                                # dX_A/d(k_sC^n t) at t = 0
RATE0_B = e_A/(1-e_A)*(2 - 1/(G_A - 1))/rpo_A           # d(eq. 6)/ds at s = 0, over r_po
KT_EARLY = brentq(dAB, 1e-12, 1e-5, xtol=1e-24, rtol=8.9e-16)
X_EARLY_CROSS = XA(KT_EARLY)
kt_dip = minimize_scalar(dAB, bounds=(1e-12, KT_EARLY), method="bounded",
                         options=dict(xatol=1e-16)).x
DAB_EARLY_DIP = dAB(kt_dip)

# The convention the paper does NOT take: fix curve B's Petersen structure by
# (eps_o, S_o) and let L follow from eqs. (3),(4), instead of by (eps_o, L_o).
L_B_ALT  = So_A**2/(4*u0_A**2*(1-u0_A)**2*3*np.pi)
rpo_ALT  = u0_A/np.sqrt(np.pi*L_B_ALT/3)
dAB_ALT  = lambda kt: XA(kt) - X_petersen(kt/rpo_ALT, e_A)
KT_END_ALT = rpo_ALT*(2*G_A/3 - 1)
DAB_UPTO75_ALT = float(np.max(np.abs(dAB_ALT(np.linspace(1e-9, KT_75, 4001)))))
KTS_ALT  = brentq(dAB_ALT, 1e-5, KT_END_ALT*0.99, xtol=1e-20, rtol=8.9e-16)
ktn_ALT  = minimize_scalar(dAB_ALT, bounds=(KTS_ALT, KT_END_ALT), method="bounded",
                           options=dict(xatol=1e-18)).x
X_AT_2PT_ALT = XA(brentq(lambda z: abs(dAB_ALT(z)) - 0.02, KTS_ALT, ktn_ALT,
                         xtol=1e-20, rtol=8.9e-16))

# The THIRD pairing, so the band is closed rather than represented by one member:
# fix curve B by (S_o, L_o) and let eq. (4) give the porosity -- which is the same
# structure whose 2.5 % miss on eps_o = 0.3 section 6.2 already reports.
rpo_SL   = SETS[1]["rp_from_S"]; eps_SL = SETS[1]["eps_pred"]
dAB_SL   = lambda kt: XA(kt) - X_petersen(kt/rpo_SL, eps_SL)
KT_END_SL = rpo_SL*(2*G_of_eps(eps_SL)/3 - 1)
_gSL = np.linspace(1e-9, KT_75, 4001); _aSL = np.abs(dAB_SL(_gSL)); _iSL = int(np.argmax(_aSL))
if 0 < _iSL < len(_gSL) - 1:      # an INTERIOR maximum here, so refine it; never sample
    _rSL = minimize_scalar(lambda z: -abs(dAB_SL(z)), bounds=(_gSL[_iSL-1], _gSL[_iSL+1]),
                           method="bounded", options=dict(xatol=1e-18))
    DAB_UPTO75_SL, DAB_SL_GRID = -float(_rSL.fun), float(_aSL[_iSL])
else:                             # the maximum is on the X = 0.75 endpoint: exact already
    DAB_UPTO75_SL = DAB_SL_GRID = float(_aSL[_iSL])
KTS_SL   = brentq(dAB_SL, 1e-5, KT_END_SL*0.99, xtol=1e-20, rtol=8.9e-16)
ktn_SL   = minimize_scalar(dAB_SL, bounds=(KTS_SL, KT_END_SL), method="bounded",
                           options=dict(xatol=1e-18)).x
X_AT_2PT_SL = XA(brentq(lambda z: abs(dAB_SL(z)) - 0.02, KTS_SL, ktn_SL,
                        xtol=1e-20, rtol=8.9e-16))

# The surfaces, on BOTH printed parameter sets, so the -21.7 % of the title cell
# is never read as belonging to curves A and B.
surf_dev = lambda Xv, psi_, e_: float(SS0_of_X(Xv, psi_)
                                      / SS0_petersen(np.array([Xv]), e_)[0] - 1)
SURF_DEV_F4_25, SURF_DEV_F4_75 = (surf_dev(0.25, PSI_F4, P["fig4_eps_o"]),
                                  surf_dev(0.75, PSI_F4, P["fig4_eps_o"]))
SURF_DEV_AB_25, SURF_DEV_AB_75 = surf_dev(0.25, PSI_A, e_A), surf_dev(0.75, PSI_A, e_A)

ktp = np.linspace(1e-9, KT_END*1.02, 900)
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.1))
ax[0].plot(ktp*1e4, XA(ktp), lw=2, label=f"A: random pore, $\\psi$ = {PSI_A:.3f}")
ax[0].plot(ktp[ktp <= KT_END]*1e4, XB(ktp[ktp <= KT_END]), "--", lw=2, label="B: Petersen, eqs. (6),(7)")
So_C = P["fig7_C_So"]
ax[0].plot(ktp*1e4, X_of_tau(ktp*So_C/(1-e_A), psi_of(So_C, P["fig7_C_Lo"], e_A)), lw=2,
           label=f"C: random pore, $\\psi$ = {psi_of(So_C, P['fig7_C_Lo'], e_A):.3f}")
u0D = u_of_eps(e_A); rpoD = u0D/np.sqrt(np.pi*L_D/3)
ktD = ktp[ktp <= rpoD*(2*G_of_eps(e_A)/3 - 1)]
ax[0].plot(ktD*1e4, X_petersen(ktD/rpoD, e_A), ":", lw=2, label="D: Petersen, L reconstructed")
ax[0].axhline(P["fig7_agreement_X"], color="0.5", lw=0.9)
ax[0].set(xlabel=r"time parameter $(k_s C^n t)\times 10^{4}$", ylabel="conversion $X$",
          ylim=(0, 1.02), title="the construction of Fig. 7, from printed parameters")
ax[0].legend(fontsize=8, loc="lower right")

ax[1].plot(XA(ktp[ktp <= KT_END]), dAB(ktp[ktp <= KT_END]), lw=2, color="C3")
ax[1].axhline(0, color="k", lw=0.8); ax[1].axhline(0.02, color="0.6", lw=0.8, ls=":")
ax[1].axhline(-0.02, color="0.6", lw=0.8, ls=":")
ax[1].axvline(P["fig7_agreement_X"], color="0.5", lw=0.9)
ax[1].plot([X_AT_2PT], [dAB(KT_2PT)], "ko", ms=5)
ax[1].set(xlabel="conversion $X$ (curve A)", ylabel="$X_A - X_B$",
          title="the printed '75 %' claim, measured")
fig.tight_layout(); plt.show()

print(f"curve A: psi = {PSI_A:.6f} (> 2, so it HAS a rate maximum, at X_M = {XM_eq33(PSI_A):.6f})")
print(f"curve B: r_po = {rpo_A:.6e} cm from eq. (4); complete gasification at k_s C^n t = {KT_END:.6e}")
print(f"\nthe two curves cross TWICE. Petersen's initial rate is the faster:"
      f" dX/d(k_sC^n t)|_0 = {RATE0_B:.4f}")
print(f"   against {RATE0_A:.4f} for the random pore model ({RATE0_B/RATE0_A - 1:+.3%}), so A TRAILS B")
print(f"   up to X = {X_EARLY_CROSS:.6f} before overtaking; the early dip is only"
      f" {DAB_EARLY_DIP:+.4e} in X.")
print(f"   the second crossing -- the one the comparison is usually quoted at -- is at"
      f" X = {XA(KT_SIGN):.6f};")
print(f"   A leads B by at most {DAB_POS:+.6f} (at X = {XA(kt_pos):.4f}) and trails it by at most")
print(f"   {DAB_NEG:+.6f} (at X = {XA(kt_neg):.4f}) over the whole history.")
print(f"\nmax |X_A - X_B| for X <= {P['fig7_agreement_X']}:  {DAB_UPTO75:.6f} conversion"
      f"  ({DAB_UPTO75*100:.2f} percentage points)")
print(f"the two first part by 0.02 in X at X = {X_AT_2PT:.6f}"
      f"  -- the printed 'about {P['fig7_agreement_X']*100:.0f} % conversion', to three digits.")

print(f"\nTHAT THIRD DIGIT IS THE PAPER'S CHOICE, NOT THIS PAGE'S. Curve B is fixed by")
print(f"   (eps_o, L_o) because p. 384 says so verbatim ('If eps_o = {P['fig7_eps_o']} and L_o ="
      f" {P['fig7_AB_Lo']:.3g}\n   are chosen, the results are the same as the prior curve B').")
print(f"   An over-determined triple has THREE pairings, so the OTHER TWO are both printed and")
print(f"   the band is closed -- neither is 'the' alternative:")
print(f"   (eps_o, S_o), the pairing p. 384 assigns to curve D: implies L = {L_B_ALT:.6e} cm/cm3"
      f" ({L_B_ALT/P['fig7_AB_Lo'] - 1:+.2%}")
print(f"      on the printed length), max |X_A - X_B| = {DAB_UPTO75_ALT:.6f}, parting conversion"
      f" {X_AT_2PT_ALT:.6f}")
print(f"   (S_o, L_o), which drops eps_o = {P['fig7_eps_o']} and takes the {SETS[1]['eps_rel']:.4%}-adrift"
      f" eps = {eps_SL:.6f} of 6.2:")
print(f"      max |X_A - X_B| = {DAB_UPTO75_SL:.6f} (interior maximum, refined off the 4001-point")
print(f"      scan by {abs(DAB_UPTO75_SL - DAB_SL_GRID):.1e}), parting conversion {X_AT_2PT_SL:.6f}")
print(f"   the sourced (eps_o, L_o):  {DAB_UPTO75:.6f} and {X_AT_2PT:.6f} -- so the parting conversion")
print(f"   spans {min(X_AT_2PT, X_AT_2PT_ALT, X_AT_2PT_SL):.4f} to"
      f" {max(X_AT_2PT, X_AT_2PT_ALT, X_AT_2PT_SL):.4f} over the three readings: the SECOND digit moves,")
print("   all three round to 'about 75 %', and only the paper's own convention earns the")
print("   three-digit statement.")

print(f"\nthe SURFACES part far earlier than the conversions, on either printed set:")
print(f"   Fig. 4 / p. 384 Petersen set (eps_o = {P['fig4_eps_o']}, psi = {PSI_F4:.4f}):"
      f" {SURF_DEV_F4_25:+.4%} at X = 0.25, {SURF_DEV_F4_75:+.4%} at X = 0.75")
print(f"   Fig. 7 curves A, B (eps_o = {P['fig7_eps_o']}, psi = {PSI_A:.4f}):"
      f"          {SURF_DEV_AB_25:+.4%} at X = 0.25, {SURF_DEV_AB_75:+.4%} at X = 0.75")'''))

cells.append(md(r"""### 6.4 Why they part: two different overlap laws

The abstract says the models "approach" each other at zero pore-size variance,
not that they coincide, and the paper's explanation is "the neglect of new
intersections of reaction surfaces in the Petersen approach". That difference
is one line of algebra. Both models write the porosity as a function of the
*non-overlapped* volume $V_E$; Avrami's eq. (20) gives $V = 1-e^{-V_E}$, while
Petersen's eqs. (3)–(4) in $u = r\sqrt{\pi L/3}$ give $V_E = 3u^2$ and
$V = 3u^2-2u^3$. Expanded,

$$V_{\text{Avrami}} = V_E - \tfrac12 V_E^2 + \ldots, \qquad
V_{\text{Petersen}} = V_E - 2\left(\tfrac{V_E}{3}\right)^{3/2},$$

so they differ already at order $V_E^{3/2}$: **no parameter choice makes the
two models identical**, and the paper's "approach … over a large range of
conversions" is the strongest statement available. Petersen's law also reaches
$V = 1$ at finite $u = 1$ while Avrami's only approaches it, which is where
the high-conversion divergence of §6.3 comes from."""))

cells.append(code(r'''ug = np.linspace(0, 1, 1201)
V_avr, V_pet = 1 - np.exp(-3*ug**2), 3*ug**2 - 2*ug**3
gapf = lambda z: (1 - np.exp(-3*z**2)) - (3*z**2 - 2*z**3)
U_CROSS = brentq(gapf, 0.4, 0.9, xtol=1e-17, rtol=8.9e-16)
u_mg = minimize_scalar(lambda z: -gapf(z), bounds=(0.05, 0.6), method="bounded",
                       options=dict(xatol=1e-15)).x
OVERLAP_MAX_GAP = gapf(u_mg)

fig, ax = plt.subplots(1, 2, figsize=(11.2, 3.7))
ax[0].plot(3*ug**2, V_avr, lw=2, label="Avrami, eq. (20):  $1-e^{-V_E}$")
ax[0].plot(3*ug**2, V_pet, "--", lw=2, label="Petersen, eqs. (3),(4):  $3u^2-2u^3$")
ax[0].set(xlabel="$V_E$ (non-overlapped volume)", ylabel="$V$ (porosity)", xlim=(0, 3), ylim=(0, 1.02))
ax[0].legend(fontsize=8)
ax[1].plot(3*ug**2, V_avr - V_pet, lw=2, color="C3"); ax[1].axhline(0, color="k", lw=0.8)
ax[1].plot([3*u_mg**2], [OVERLAP_MAX_GAP], "ko", ms=5)
ax[1].set(xlabel="$V_E$", ylabel="Avrami $-$ Petersen", xlim=(0, 3))
fig.tight_layout(); plt.show()

print(f"the two overlap laws agree only at V_E = 0 and cross once, at u = {U_CROSS:.7f}"
      f" (V_E = {3*U_CROSS**2:.6f});")
print(f"largest gap {OVERLAP_MAX_GAP:.6f} in porosity, at u = {u_mg:.6f}.")
for VE in (0.01, 0.05, 0.1):
    lead = 2*(VE/3)**1.5 - VE**2/2
    print(f"   V_E = {VE:<5g}: gap {(1-np.exp(-VE)) - (VE - 2*(VE/3)**1.5):+.4e}"
          f"   leading term 2(V_E/3)^(3/2) - V_E^2/2 = {lead:+.4e}")
OVERLAP_LEAD_REL = abs(((1-np.exp(-0.01)) - (0.01 - 2*(0.01/3)**1.5))
                       / (2*(0.01/3)**1.5 - 0.01**2/2) - 1)
print(f"   at V_E = 0.01 the leading term already accounts for the gap to"
      f" {OVERLAP_LEAD_REL:.2e} relative.")'''))

cells.append(md(r"""### 6.5 The grain shape factor: the paper's own least squares, re-run

The paper defines $m$ by minimising $I$ over the full conversion range and
reports the result only as Figure 6 plus two prose sentences. Both sentences
are numbers here. $I$ has a closed form — substituting $w = 1-X$ and
$w = e^{-v}$ turns every term into an incomplete gamma function,

$$\int_0^1 w^{b}\sqrt{1-\psi\ln w}\,dw
= \frac{e^{c}}{\psi}\left(\frac{\psi}{b+1}\right)^{3/2}\Gamma(3/2, c),
\qquad c = \frac{b+1}{\psi},$$

so $I = \left(\tfrac13+\tfrac{\psi}{9}\right) - 2J(m{+}1) + \tfrac{1}{2m+1}$
with $J$ that integral. §7.3 checks the closed form against two independent
quadratures; here it is used to root-find $dI/dm = 0$ rather than to sample a
grid.

Reported alongside is a **normalised** misfit,
$\sqrt{I/\int_0^1 (dX/d\tau)^2 dX}$ with the denominator
$\tfrac13+\tfrac{\psi}{9}$ — the same integral with the grain term removed. It
is what makes "grossly inadequate" a number."""))

cells.append(code(r'''@lru_cache(maxsize=None)
def _gl(n):
    """Gauss-Legendre nodes, cached: leggauss(600) is not cheap and this page
    calls it inside every root-find iteration."""
    return leggauss(n)

def J_gam(b, psi_):                      # closed form; psi = 0 is the plain integral
    if psi_ == 0: return 1.0/(b+1)
    c = (b+1)/psi_
    return np.exp(c)/psi_*(psi_/(b+1))**1.5*gammaincc(1.5, c)*gamma_fn(1.5)

def I_closed(psi_, m_):  return (1/3 + psi_/9) - 2*J_gam(m_+1, psi_) + 1/(2*m_+1)
def I_energy(psi_):      return 1/3 + psi_/9         # int (dX/dtau)^2 dX, the scale
def I_nrms(psi_, m_):    return np.sqrt(max(I_closed(psi_, m_), 0.0)/I_energy(psi_))

def dJ_db(b, psi_, n=600):               # dJ/db = int w^b ln w sqrt(1-psi ln w) dw
    x, wt = _gl(n); v = 40*(x+1); W = 40*wt
    return float(np.sum(W*(-v)*np.exp(-(b+1)*v)*np.sqrt(1 + psi_*v)))
def dI_dm(psi_, m_):     return -2*dJ_db(m_+1, psi_) - 2/(2*m_+1)**2

@lru_cache(maxsize=None)
def m_star(psi_):        # ROOT-FOUND, never sampled
    return brentq(lambda z: dI_dm(psi_, z), -0.45, 3.0, xtol=1e-15, rtol=8.9e-16)
def m_star_brent(psi_):  # independent route: direct minimisation of the closed form
    lo = max(m_star(psi_) - 0.25, -0.4)
    return minimize_scalar(lambda z: I_closed(psi_, z), bounds=(lo, lo + 0.5),
                           method="bounded", options=dict(xatol=1e-13)).x

M_PSI0, M_PSI1, M_PSI2 = m_star(0.0), m_star(1.0), m_star(P["m_psi_window_hi"])
M_PSI1_REL_23 = (M_PSI1 - P["grain_m_sphere"])/P["grain_m_sphere"]
psi_fine = np.linspace(0, P["m_psi_window_hi"], 41)
m_fine = np.array([m_star(float(pv)) for pv in psi_fine])
assert np.all(np.diff(m_fine) < 0), "m*(psi) must be monotone on the printed window"
M_WINDOW_LO, M_WINDOW_HI = m_fine.min(), m_fine.max()
M_TWO_ROUTES = max(abs(m_star(pv) - m_star_brent(pv)) for pv in (0.5, 1.0, 1.5, 2.0))

fig, ax = plt.subplots(1, 2, figsize=(11.2, 3.9))
psi_wide = np.linspace(0, 6, 61); m_wide = np.array([m_star(float(pv)) for pv in psi_wide])
ax[0].plot(psi_wide, m_wide, lw=2)
ax[0].axvspan(0, P["m_psi_window_hi"], color="C0", alpha=0.08)
ax[0].axhline(P["m_range_lo"], color="0.6", lw=0.8, ls=":")
ax[0].axhline(P["grain_m_sphere"], color="0.6", lw=0.8, ls="--")
ax[0].plot([1.0], [M_PSI1], "ko", ms=5)
ax[0].set(xlabel=r"pore structure parameter $\psi$", ylabel="best-fit $m$",
          title="the construction of Fig. 6 (shaded: the printed window)", ylim=(-0.1, 1.05))
ax[1].plot(psi_wide, [I_nrms(pv, max(m_star(float(pv)), 0.0)) for pv in psi_wide], lw=2, label="best $m$")
ax[1].plot(psi_wide, [I_nrms(pv, P["grain_m_flatplate"]) for pv in psi_wide], "--", lw=1.7,
           label="$m$ = 0 (flat plates)")
ax[1].plot(psi_wide, [I_nrms(pv, P["grain_m_sphere"]) for pv in psi_wide], ":", lw=1.7,
           label="$m$ = 2/3 (spheres)")
ax[1].plot(psi_wide, [I_nrms(pv, P["grain_m_volume"]) for pv in psi_wide], "-.", lw=1.7,
           label="$m$ = 1 (volume model)")
ax[1].axvspan(0, P["m_psi_window_hi"], color="C0", alpha=0.08)
ax[1].set(xlabel=r"$\psi$", ylabel="normalised rms misfit", ylim=(0, 1.02),
          title="how well ANY order-of-reaction model can do")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()

print(f"printed: 'For 0 <= psi <= {P['m_psi_window_hi']:.0f}, the best fit gives"
      f" {P['m_range_lo']} <= m <= {P['m_range_hi']:.0f}'")
print(f"   computed over the same window: {M_WINDOW_LO:.7f} <= m <= {M_WINDOW_HI:.7f}"
      f"   (m*(0) = {M_PSI0:.10f} exactly the volume model; m*({P['m_psi_window_hi']:.0f})"
      f" = {M_PSI2:.7f} -> {P['m_range_lo']})")
print(f"printed: 'For psi = 1, the best m is very close to two-thirds'")
print(f"   computed: m*(1) = {M_PSI1:.7f}, which is {abs(M_PSI1_REL_23):.3%} below 2/3"
      f" = {P['grain_m_sphere']:.7f}")
print(f"   two independent routes to m* (root-find of dI/dm vs Brent on I) agree to {M_TWO_ROUTES:.2e}")

# --- every printed m exercised, including both closed-form branches ---------
print("\nevery printed m, exercised (eq. 39's closed form has a BRANCH at m = 1):")
for m_, lab in ((P["grain_m_flatplate"], "flat plates"), (P["grain_m_cylinder"], "cylinders"),
                (P["grain_m_sphere"], "spheres"), (P["grain_m_volume"], "volume model")):
    tt = 0.4
    print(f"   m = {m_:.6f} ({lab:<13s}) X(tau=0.4) = {X_grain(tt, m_):.10f}"
          f"   nrms at psi = 1: {I_nrms(1.0, m_):.4%}")
BRANCH_GAP = abs(X_grain(0.4, 1.0) - (1 - (1-(1-(1-1e-9))*0.4)**(1/(1-(1-1e-9)))))
print(f"   the m -> 1 branch is continuous: |X(m=1) - X(m=1-1e-9)| = {BRANCH_GAP:.2e}")

# --- the m = 0 claim, quantified beyond the printed window -----------------
PSI_M0 = brentq(lambda pv: dI_dm(pv, 0.0), 1.0, 100.0, xtol=1e-14, rtol=8.9e-16)
NRMS_AT_PSI_M0 = I_nrms(PSI_M0, 0.0)
NRMS_M0_PSI2 = I_nrms(P["m_psi_window_hi"], P["grain_m_flatplate"])
NRMS_BEST_PSI2 = I_nrms(P["m_psi_window_hi"], M_PSI2)
print(f"\nthe abstract's 'not consistent with the assumption of flat plate grains (m = 0) for any")
print(f"pore structure', and p. 384's 'grossly inadequate fit for any psi':")
print(f"   at psi = {P['m_psi_window_hi']:.0f}, m = 0 misfits by {NRMS_M0_PSI2:.2%} against"
      f" {NRMS_BEST_PSI2:.2%} for the best m -- a factor {NRMS_M0_PSI2/NRMS_BEST_PSI2:.2f}")
print(f"   m = 0 IS the least-squares optimum at exactly one psi, psi = {PSI_M0:.6f} (outside the")
print(f"   printed window) -- and even there its misfit is {NRMS_AT_PSI_M0:.2%}, because past psi = 2")
print("   NO order-of-reaction model fits: none of them can produce a rate maximum. The paper's")
print("   claim is about fit QUALITY and it holds. The optimum-at-one-psi nuance is NOT printed:")
print("   searched the abstract and CONCLUSIONS (p. 379), the whole m discussion around Figure 6")
print("   (p. 384 left column) and the Figure 6 caption, whose psi axis stops at 2 -- the paper")
print("   never evaluates m* outside 0 <= psi <= 2, so the statement is this page's, not its.")'''))

cells.append(md(r"""### 6.6 The printed defects

Three, each proved from the paper's own printed content and none repaired.

**Defect 2 — "This match at $m = 1$" (p. 383, right column).** Verbatim, with
the sentence before it: *"In either figure, it is clear that the grain model
predictions match the $\psi = 1$ results within the limits of the usual
experimental error. This match at $m = 1$ combined with the above-mentioned
correspondence as $\psi \to 0$ suggests that…"* [*sic*]. The dashed line in
Figures 3 and 4 is labelled "GRAIN MODEL, $m = 2/3$", not $m = 1$; the
"above-mentioned correspondence as $\psi \to 0$" is already the $m = 1$ case
(eq. 37 reducing to the volume-reaction model), so reading the sentence that
way makes it say the same thing twice; and Figure 6 puts $m = 1$ at $\psi = 0$
and $m \approx 2/3$ at $\psi = 1$. The least squares below settles it. Whether
the slip is "$m = 1$" for "$\psi = 1$" or for "$m = 2/3$" cannot be decided
from the document, so neither repair is asserted.

**Defect 3 — the optimal-porosity worked example (p. 383).** *"If, for
example, a final $\epsilon = 0.3$ is sought for a solid where
$(S_o^2/2\pi L_o) = 0.67$, the initial $\epsilon_o$ should be 0.1 to maximize
the final surface area."* The printed 0.1 is not a root of the paper's own
eq. (35); the root is 0.08529. The stationarity residual at 0.1 is −4.2e-3
against 0 at the root. This one is a defect of *precision*, not of substance:
the optimum is so flat that 0.1 gives a surface 0.003 % below the maximum, and
Figure 5's resolution cannot distinguish them — which is presumably where 0.1
came from. Reported anyway, because eq. (35) is an equation.

Defect 1 (Park & Levenspiel, 1976 in the text against 1975 in the reference
list) is in §4. Two typographical slips are recorded without further comment:
"not clearly indentified" [*sic*, p. 382] and "This research was support by
the U.S. Department of Energy" [*sic*, p. 385]; the NOTATION list spells
"stoichimetric" [*sic*] for $p,q$ and "stoichiometric" for $a,b$ two lines
above."""))

cells.append(code(r'''# --- defect 2: what the least squares says at psi = 1 ----------------------
NRMS_M1_PSI1   = I_nrms(1.0, P["grain_m_volume"])
NRMS_BEST_PSI1 = I_nrms(1.0, M_PSI1)
NRMS_M23_PSI1  = I_nrms(1.0, P["grain_m_sphere"])
M1_PENALTY = NRMS_M1_PSI1/NRMS_BEST_PSI1
print("PRINTED DEFECT 2 -- 'This match at m = 1' [sic], p. 383:")
print(f"   at psi = 1 the best-fit m is {M_PSI1:.7f}; the printed m = 1 misfits by"
      f" {NRMS_M1_PSI1:.4%}")
print(f"   against {NRMS_BEST_PSI1:.4%} at the optimum and {NRMS_M23_PSI1:.4%} at m = 2/3"
      f" -- a factor {M1_PENALTY:.2f} worse.")
print(f"   m = 1 IS the exact optimum, but at psi = 0: m*(0) = {M_PSI0:.10f}, misfit"
      f" {I_nrms(0.0, 1.0):.1e}.")

# --- defect 3: eq. (35)'s worked example -----------------------------------
def S_of_eps_o(e0, ef, Gh):                      # eq. (27) with V = eps, V_o = eps_o
    return (1-ef)/(1-e0)*np.sqrt(1 - (2*(1-e0)/Gh)*np.log((1-ef)/(1-e0)))
def stat35(e0, ef, Gh):  return Gh/(1-e0) - 1 - np.log((1-ef)/(1-e0))

EPS_F, GHAT = P["opt_eps_final"], P["opt_group"]
EPS_O_ROOT   = float(eps_o_opt_block(EPS_F, GHAT)[0])          # pymrm block solve
EPS_O_DIRECT = minimize_scalar(lambda z: -S_of_eps_o(z, EPS_F, GHAT),
                               bounds=(1e-9, EPS_F-1e-9), method="bounded",
                               options=dict(xatol=1e-15)).x     # independent route
EPS_O_TWO_ROUTES = abs(EPS_O_ROOT - EPS_O_DIRECT)
STAT35_AT_PRINTED = stat35(P["opt_eps_o"], EPS_F, GHAT)
S_DEFICIT_AT_PRINTED = 1 - S_of_eps_o(P["opt_eps_o"], EPS_F, GHAT)/S_of_eps_o(EPS_O_ROOT, EPS_F, GHAT)
IN36_LO, IN36_HI = 1 + np.log(1-EPS_F), 1 - EPS_F
EPS_O_SECOND = brentq(lambda z: stat35(z, EPS_F, GHAT), EPS_F+1e-9, 0.999, xtol=1e-17, rtol=8.9e-16)

e0g = np.linspace(1e-4, EPS_F-1e-4, 400)
fig, ax = plt.subplots(figsize=(5.6, 3.4))
ax.plot(e0g, [S_of_eps_o(z, EPS_F, GHAT) for z in e0g], lw=2)
ax.plot([EPS_O_ROOT], [S_of_eps_o(EPS_O_ROOT, EPS_F, GHAT)], "ko", ms=6, label=f"eq. (35) root {EPS_O_ROOT:.5f}")
ax.plot([P["opt_eps_o"]], [S_of_eps_o(P["opt_eps_o"], EPS_F, GHAT)], "rs", ms=6,
        label=f"printed {P['opt_eps_o']}")
ax.set(xlabel=r"initial porosity $\epsilon_o$", ylabel="$S/S_o$ at $\\epsilon$ = 0.3")
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

print("\nPRINTED DEFECT 3 -- the worked example of p. 383:")
print(f"   inequality (36) window at eps = {EPS_F}: [{IN36_LO:.7f}, {IN36_HI:.7f}];"
      f" printed threshold {P['opt_group_violated']}")
print(f"      1 + ln(1 - {EPS_F}) = {IN36_LO:.7f}, so the printed '<= {P['opt_group_violated']}"
      f" violates (36)' is the rounded-down bound. Consistent.")
print(f"   eq. (35) at eps = {EPS_F}, S_o^2/2 pi L_o = {GHAT}: root eps_o = {EPS_O_ROOT:.9f}"
      f"   (PRINTED: {P['opt_eps_o']})")
print(f"      the two routes (pymrm block Newton on eq. 35; direct maximisation of eq. 27)"
      f" agree to {EPS_O_TWO_ROUTES:.2e}")
print(f"      stationarity residual of eq. (35) at the printed {P['opt_eps_o']}:"
      f" {STAT35_AT_PRINTED:+.6e}, at the root {stat35(EPS_O_ROOT, EPS_F, GHAT):+.1e}")
print(f"      and the cost of the printed value: S is {S_DEFICIT_AT_PRINTED:.3e}"
      f" ({S_DEFICIT_AT_PRINTED:.4%}) below the maximum -- immaterial, which is the point.")
print(f"   eq. (35) has a SECOND root at eps_o = {EPS_O_SECOND:.7f} > eps: it is a MINIMUM"
      f" (S/S_o = {S_of_eps_o(EPS_O_SECOND, EPS_F, GHAT):.6f}")
print(f"      against {S_of_eps_o(EPS_O_ROOT, EPS_F, GHAT):.6f} at the maximum) and eps_o > eps"
      " is unphysical for a reaction that opens porosity.")'''))

cells.append(md(r"""### 6.7 The Hashimoto chars — what the prose alone supports

The three numbers printed in prose on p. 385 are enough to state everything
the model predicts for those two chars, **without touching Figure 8**. The
linearisation §2.1 verified,

$$\left[\frac{S}{1-X}\right]^2 = S_o^2\left(1 + \psi\ln\frac{1}{1-X}\right),$$

says the intercept is $S_o^2$ and the slope $\psi S_o^2$, which is how the
paper got its numbers out of the figure; running that backwards turns the
printed $S_o$ and $\psi$ into the two lines below and into a length-averaged
pore radius, since $\bar r \equiv S_o/2\pi L_o = 2(1-\epsilon_o)/\psi S_o$
follows from the definition of $\psi$ alone.

**The effectiveness-factor claim cannot be checked here** — "Thiele moduli
based on the initial surface areas are small enough to assure effectiveness
factors over 0.98 in all cases" needs a particle size, a diffusivity and a
rate constant, and the paper prints none of the three. What *is* computable is
what the claim demands of the modulus, which is stated below as a bound and
nothing more."""))

cells.append(code(r'''SO_CHAR = P["hashimoto_So"]*1e4                     # m2/cm3 -> cm2/cm3
INTERCEPT = SO_CHAR**2
chars = []
for tag, psi_ in (("VC", P["hashimoto_psi_VC"]), ("CS", P["hashimoto_psi_CS"])):
    chars.append(dict(char=tag, psi=psi_, slope_e12=psi_*INTERCEPT/1e12,
                      X_M=XM_eq33(psi_), rbar_nm_per_1meps=2/(psi_*SO_CHAR)*1e7))
ch = pd.DataFrame(chars)
with pd.option_context("display.float_format", lambda v: f"{v:.6g}"):
    display(ch)

XM_VC, XM_CS = XM_eq33(P["hashimoto_psi_VC"]), XM_eq33(P["hashimoto_psi_CS"])
RBAR_VC = 2/(P["hashimoto_psi_VC"]*SO_CHAR)*1e7
RBAR_CS = 2/(P["hashimoto_psi_CS"]*SO_CHAR)*1e7

xx = np.linspace(0, 1.45, 200)
fig, ax = plt.subplots(figsize=(5.6, 3.4))
for tag, psi_, st in (("VC", P["hashimoto_psi_VC"], "-"), ("CS", P["hashimoto_psi_CS"], "--")):
    ax.plot(xx, INTERCEPT*(1 + psi_*xx)/1e12, st, lw=2, label=fr"{tag} char, $\psi$ = {psi_}")
ax.set(xlabel=r"$\ln[1/(1-X)]$", ylabel=r"$[S/(1-X)]^2 \times 10^{-12}$  (cm$^{-2}$)",
       title="the two printed correlations — NO data points (Fig. 8 is not digitised)")
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

print(f"printed in prose: S_o = {P['hashimoto_So']:.0f} m2/cm3 for BOTH chars"
      f" = {SO_CHAR:.3e} cm2/cm3")
print(f"   -> the shared intercept the linearisation demands is S_o^2 ="
      f" {INTERCEPT/1e12:.3f} x 10^12 cm^-2, and the two slopes are")
print(f"      psi S_o^2 = {P['hashimoto_psi_VC']*INTERCEPT/1e12:.2f} and"
      f" {P['hashimoto_psi_CS']*INTERCEPT/1e12:.2f} x 10^12 cm^-2. Two lines from one intercept.")
print(f"   both psi exceed the printed threshold {P['psi_ratemax']:.0f}, so eq. (33) predicts surface")
print(f"   maxima at X_M = {XM_VC:.6f} (VC) and {XM_CS:.6f} (CS) -- the paper says a maximum")
print("   'was, in fact, reported in the reference work for the CS char but did not show up for")
print("   the VC char, obscured by the particle attrition that occurred in the fluidized bed used'.")
print(f"   length-averaged pore radius, from psi's definition alone:"
      f" {RBAR_VC:.4f}(1-eps_o) nm (VC), {RBAR_CS:.4f}(1-eps_o) nm (CS)")
print("      -- micropores at any plausible porosity, which is what steam-activated char is.")

def eta_sphere(phi): return 3*(phi/np.tanh(phi) - 1)/phi**2
PHI_098 = brentq(lambda z: eta_sphere(z) - P["eta_min"], 1e-6, 10.0, xtol=1e-16, rtol=8.9e-16)
print(f"\nthe effectiveness-factor sentence: eta > {P['eta_min']} for a sphere requires"
      f" a Thiele modulus phi < {PHI_098:.6f}")
print("   (phi^2 < %.6f). NOT CHECKABLE HERE: the paper prints no particle size, no diffusivity"
      % PHI_098**2)
print("   and no rate constant, so the moduli behind that sentence cannot be reconstructed.")'''))

# ----------------------------------------------------------------- validation
cells.append(md(r"""## Validation

Nothing here is validated against data — there is none. What is validated is
that this page computes the paper's equations correctly, by routes that share
no assembly, and that every number it reports moves when the physics is
perturbed."""))

cells.append(md(r"""### 7.1 The rate law marched, against the closed form it never sees

Eq. (37) is the paper's *rate* law and eq. (32) its *integral*. Marching (37)
with Crank–Nicolson — pymrm `newton` on each implicit step, vectorised over
$\psi$ — and comparing against (32) tests the printed pair against each other
through machinery that shares nothing with either. The time step is refined
and the observed order reported; a wrong-baseline error would show up as an
order that is right and an error that plateaus."""))

cells.append(code(r'''PSI_MARCH = np.array([0.0, 1.0, 2.0, 5.0, 10.0]); TAU_END = 1.5
rows, prev = [], None
for nst in (50, 100, 200, 400, 800):
    Xn = march37(PSI_MARCH, TAU_END, nst)
    err = float(np.max(np.abs(Xn - X_of_tau(TAU_END, PSI_MARCH))))
    rows.append((nst, TAU_END/nst, err, np.nan if prev is None else np.log2(prev/err)))
    prev = err
mtab = pd.DataFrame(rows, columns=["steps", "d_tau", "max |X_CN - eq.(32)|", "observed order"])
display(mtab)
MARCH_ERR_800 = rows[-1][2]
MARCH_ORDER = float(rows[-1][3])
XeulerErr = float(np.max(np.abs(march37(PSI_MARCH, TAU_END, 800, theta=1.0)
                                - X_of_tau(TAU_END, PSI_MARCH))))
print(f"Crank-Nicolson: {MARCH_ERR_800:.3e} at 800 steps, observed order {MARCH_ORDER:.3f}")
print(f"implicit Euler at the same 800 steps: {XeulerErr:.3e}  ({XeulerErr/MARCH_ERR_800:.0f}x worse)")
print("-> eqs. (32) and (37) are the same model, and this page integrates it correctly.")'''))

cells.append(md(r"""### 7.2 Everything else, two ways

Six quantities, each computed by two routes that share no *code*. Four of them
share no *algebra* either — (a) $X_M$ closed-form against a root-find on
eq. (31), (b) the pymrm block solver against `brentq` and against the
$\sigma\to\infty$ closed form, (e) the incomplete-gamma closed form against
Gauss–Legendre quadrature.

**Two of them do share algebra, and the page says which.** Substituting
$G = 3/(2u_o)$ and $\epsilon = 3u_o^2-2u_o^3$ into eq. (6) turns it into
$[3u^2-2u^3-\epsilon]/(1-\epsilon)$ with $u = u_o(1+s)$ — route 2 of (d) *is*
route 1 with $G$ eliminated — and differentiating the same expression gives
$6u_o\,u(1-u)$, which is route 2 of (c). So **(c) and (d) are round-off on an
identity, and what they test is the selection of the root of eq. (7) and the
transcription of eq. (4), not eq. (6) itself.** Their break rows are wrong-root
rows for exactly that reason, and §7.7's table prints how far each moves. A
break row for each of the six is in §7.7."""))

cells.append(code(r'''# (a) X_M: eq. (33) closed form vs root-finding dS/dX = 0 on eq. (31)
def XM_rootfind(psi_):
    d = lambda X_: (SS0_of_X(X_+1e-7, psi_) - SS0_of_X(X_-1e-7, psi_))/2e-7
    return brentq(d, 1e-9, 0.9, xtol=1e-15, rtol=8.9e-16)
XM_TWO_ROUTES = max(abs(XM_rootfind(pv) - XM_eq33(pv))/XM_eq33(pv)
                    for pv in (3.0, 5.0, P["hashimoto_psi_VC"], P["hashimoto_psi_CS"]))

# (b) tau_1/2 at sigma -> infinity: closed form vs the pymrm block solver at large sigma
TAU_HALF_CLOSED = lambda pv: np.log(2) if pv == 0 else (-1 + np.sqrt(1 + pv*np.log(2)))*2/pv
psis_b = np.array([0.0, 1.0, 2.0, 10.0, 100.0])
TAU_BLOCK_VS_SCIPY = max(
    abs(tau_at_X_block(0.5, [pv], sg)[0]
        - brentq(lambda z: X_of_tau(z, pv, sg) - 0.5, 1e-15, sg*(1-1e-13), xtol=1e-17, rtol=8.9e-16))
    / brentq(lambda z: X_of_tau(z, pv, sg) - 0.5, 1e-15, sg*(1-1e-13), xtol=1e-17, rtol=8.9e-16)
    for pv in psis_b for sg in (0.25, 1.0, 10.0, 100.0))
TAU_HALF_CLOSED_VS_BLOCK = max(
    abs(tau_at_X_block(0.5, [pv], 1e9)[0] - TAU_HALF_CLOSED(pv))/TAU_HALF_CLOSED(pv)
    for pv in psis_b)

# (c) Petersen S/S_o: from eq. (3) directly, vs differentiating eq. (6)
def SS0_pet_via_eq6(X_, e_):
    G = G_of_eps(e_); u0 = u_of_eps(e_)
    uu = u_of_eps(e_ + X_*(1-e_)); s = uu/u0 - 1
    return (2*(1+s)*(G-1-s) - (1+s)**2)/(2*G - 3)
PET_SS_TWO_ROUTES = max(abs(SS0_pet_via_eq6(x, P["fig4_eps_o"])
                            - float(SS0_petersen(np.array([x]), P["fig4_eps_o"])[0]))
                        for x in (0.05, 0.25, 0.5, 0.75, 0.9, 0.99))

# (d) Petersen conversion: eq. (6) vs the porosity route eps(u), u = u_o(1+s)
def X_pet_via_eps(s_, e_):
    u_ = u_of_eps(e_)*(1+s_); return (3*u_**2 - 2*u_**3 - e_)/(1-e_)
sgrid = np.linspace(0, 2*G_SEL_F7/3 - 1, 97)
PET_X_TWO_ROUTES = float(np.max(np.abs(X_petersen(sgrid, P["fig7_eps_o"])
                                       - X_pet_via_eps(sgrid, P["fig7_eps_o"]))))
PET_X_ENDPOINT = abs(X_petersen(2*G_SEL_F7/3 - 1, P["fig7_eps_o"]) - 1.0)

# (e) I(psi, m): closed form vs Gauss-Legendre in v = -ln w
def I_gl(psi_, m_, n=400):
    x, wt = _gl(n); v = 30*(x+1); W = 30*wt
    return float(np.sum(W*(np.exp(-v)*np.sqrt(1 + psi_*v) - np.exp(-m_*v))**2*np.exp(-v)))
I_TWO_ROUTES = max(abs(I_gl(pv, mm_) - I_closed(pv, mm_))/max(I_closed(pv, mm_), 1e-300)
                   for pv, mm_ in ((1.0, 2/3), (2.0, 0.49), (5.0, 0.4), (0.5, 0.8)))

# (f) psi: from (S_o, L_o, eps_o) vs back out of the printed Fig-4 claim psi = 5
print(f"(a) X_M closed form (33) vs root-find of dS/dX on (31):      {XM_TWO_ROUTES:.2e}")
print(f"(b) tau_1/2 pymrm block vs scipy brentq, 20 (psi,sigma):     {TAU_BLOCK_VS_SCIPY:.2e}")
print(f"    tau_1/2 pymrm block at sigma = 1e9 vs the closed form:   {TAU_HALF_CLOSED_VS_BLOCK:.2e}")
print(f"(c) Petersen S/S_o via eq. (3) vs via d(eq. 6)/ds:           {PET_SS_TWO_ROUTES:.2e}")
print(f"(d) Petersen X via eq. (6) vs via eps(u):                    {PET_X_TWO_ROUTES:.2e}")
print(f"    and eq. (6) reaches exactly X = 1 at s = 2G/3 - 1:       {PET_X_ENDPOINT:.2e}")
print(f"(e) I(psi,m) closed form vs Gauss-Legendre quadrature:       {I_TWO_ROUTES:.2e}")
print("\n    SCOPE OF (c) AND (d): with G = 3/(2 u_o) and eps = 3u^2-2u^3 substituted, eq. (6)")
print("    IS the eps(u) route and its s-derivative IS u(1-u), so both differences are round-off")
print("    on an identity. What they can detect is the ROOT of eq. (7) and the transcription of")
print("    eq. (4) -- not eq. (6)'s own coefficients. Their break rows in 7.7 are wrong-root rows.")
print("    SCOPE OF the endpoint line: X(s = 2G/3 - 1) = 1 holds for EVERY root of eq. (7), so it")
print("    tests the joint transcription of (6) and (7) and nothing about which root was taken;")
print("    7.7 measures both statements rather than asserting them.")'''))

cells.append(md(r"""### 7.3 The quadrature axis, refined

$I$ carries a discretisation error only through the quadrature used to reach
$dI/dm$; the closed form carries none. Refining the node count against the
closed form gives the observed order and shows where the quadrature saturates
— the axis that would otherwise be invisible, since $m^\star$ is a *root* of a
quadrature-based derivative."""))

cells.append(code(r'''rows, prev = [], None
for nn_ in (25, 50, 100, 200, 400):
    e = max(abs(I_gl(pv, mm_, nn_) - I_closed(pv, mm_))/I_closed(pv, mm_)
            for pv, mm_ in ((1.0, 2/3), (2.0, 0.49), (5.0, 0.4)))
    rows.append((nn_, e, np.nan if prev is None else np.log2(prev/e))); prev = e
qt = pd.DataFrame(rows, columns=["Gauss nodes", "max rel err vs closed form", "observed order"])
display(qt)
QUAD_ERR_400 = rows[-1][1]

rows2, prev = [], None
for nn_ in (50, 100, 200, 400, 600):
    mm_ = brentq(lambda z: -2*dJ_db(z+1, 1.0, nn_) - 2/(2*z+1)**2, -0.45, 3.0,
                 xtol=1e-15, rtol=8.9e-16)
    e = abs(mm_ - M_PSI1)
    rows2.append((nn_, mm_, e)); prev = e
mt = pd.DataFrame(rows2, columns=["Gauss nodes in dI/dm", "m*(psi=1)", "|shift from the 600-node value|"])
display(mt)
MSTAR_NODE_SENS = float(mt.iloc[0, 2])
print(f"m*(psi = 1) is converged in the node count to {MSTAR_NODE_SENS:.2e} already at 50 nodes;")
print(f"the reported value uses 600. The quadrature error in I itself is {QUAD_ERR_400:.2e} at 400 nodes.")'''))

cells.append(md(r"""### 7.4 The printed claims about particle size, measured

Two prose claims about Figure 2, both root-found rather than read off a grid.
$\tau_{1/2}$ is monotone decreasing in $\psi$ (asserted, not assumed), so the
spread across a $\psi$ window sits on its endpoints."""))

cells.append(code(r'''PSI_FIG2 = np.array([0.0, 1.0, 2.0, 10.0, P["sigma_small_psi_max"]])
sig_grid = np.geomspace(0.1, 1000.0, 61)
fig, ax = plt.subplots(figsize=(6.0, 3.6))
for pv in PSI_FIG2:
    ax.semilogx(sig_grid, tau_at_X_block(0.5, np.full(sig_grid.size, pv), sig_grid), lw=1.9,
                label=fr"$\psi$ = {pv:g}")
ax.axvline(P["sigma_plateau"], color="0.6", lw=0.9, ls=":")
ax.axvline(P["sigma_small"], color="0.6", lw=0.9, ls=":")
ax.set(xlabel=r"particle size parameter $\sigma = R_o S_o/(1-\epsilon_o)$",
       ylabel=r"time to 50 % conversion, $\tau_{1/2}$",
       title="the construction of Fig. 2 (dotted: the two printed thresholds)")
ax.legend(fontsize=8); fig.tight_layout(); plt.show()

t100 = tau_at_X_block(0.5, PSI_FIG2, P["sigma_plateau"])
tinf = np.array([TAU_HALF_CLOSED(pv) for pv in PSI_FIG2])
SIGMA100_MAX_REL = float(np.max(np.abs(t100 - tinf)/tinf))
t025 = tau_at_X_block(0.5, PSI_FIG2, P["sigma_small"])
assert np.all(np.diff(t025) < 0), "tau_1/2 must fall monotonically in psi at fixed sigma"
SIGMA025_SPREAD = float((t025.max() - t025.min())/t025.min())
t1 = tau_at_X_block(0.5, PSI_FIG2, 1.0)
SIGMA1_SPREAD = float((t1.max() - t1.min())/t1.min())
print(f"printed: 'no further significant change occurs for sigma >= {P['sigma_plateau']:.0f}'")
print(f"   tau_1/2 at sigma = {P['sigma_plateau']:.0f} differs from its sigma -> infinity value by at")
print(f"   most {SIGMA100_MAX_REL:.3%} over psi in {{0, 1, 2, 10, 100}} (worst at psi = 0;"
      f" only {abs(t100[-1]-tinf[-1])/tinf[-1]:.3%} at psi = 100).")
print(f"\nprinted: 'For small enough particles (sigma < {P['sigma_small']}), the reaction time is")
print(f"   independent of internal structure for psi <= {P['sigma_small_psi_max']:.0f}'")
print(f"   at sigma = {P['sigma_small']}, tau_1/2 spans {t025.min():.6f}..{t025.max():.6f}"
      f" across that whole psi range -- a spread of {SIGMA025_SPREAD:.3%}.")
print(f"   at sigma = 1, ten times larger, the same window spreads {SIGMA1_SPREAD:.1%}: the printed")
print("   threshold is doing real work, not describing the whole plot.")'''))

cells.append(md(r"""### 7.5 The reduction to the published `B3.1`, and the reconciliation it requires

The abstract's first special case is the grain model, and its integrated form
at $m = 2/3$ is eq. (40), $X = 1-(1-\tau/3)^3$. That is the reaction-control
shrinking core — the published `B3.1`'s Yagi–Kunii eq. (6) in its $\omega\to0$
limit, transcribed on that page from a 1955 Combustion Symposium scan three
months before this page transcribed eq. (40) from a 1980 AIChE J scan.

**Reading `B3.1` before using it** (that page has no `data/` directory, so
nothing is loaded; what travels is its findings):

- `B3.1` verifies its eq. (6) two independent ways — endpoint exactness for
  any $(\omega,\gamma)$, and integration of the moving boundary from eq. (5)'s
  three resistances, agreeing to 6.9e-16. The transcription this page leans on
  is therefore established there, not assumed here.
- `B3.1` reports `eq6_to_reaction_limit = 3.85e-10`, **not** zero. That is not
  a transcription discrepancy: `B3.1` measures the limit at $\omega = 10^{-9}$
  rather than at $\omega = 0$. Reproduced below at the same $\omega$ to
  confirm the reconciliation, and then taken to the exact limit.
- **What this identity cannot see — the full list, because a partial list is
  worse than none.** At $\omega = 0$ every $\omega$ term in eq. (6) drops out
  of both numerator and denominator, so the limit is **blind to all three
  $\omega$ coefficients** — including the $3\omega\gamma$ that `B3.1` had to
  read off the scan. It is also **blind to $\gamma$ itself**: at $\omega = 0$
  the numerator is $\gamma x$ and the denominator $\gamma$, so $\gamma$ cancels
  identically and eq. (6) collapses to $1-x$ with no parameter left in it. The
  $\gamma$ range printed below is therefore **not** a robustness sweep and
  cannot fail — it is printed to make that cancellation visible. And because
  the shrinking-core geometry $x = (1-X)^{1/3}$ and eq. (40)'s inversion
  $\tau/3 = 1-(1-X)^{1/3}$ are written here as the same expression, a
  *consistent* slip in that exponent would cancel too. What is left, and it is
  worth having, is that eq. (40) and eq. (6)$|_{\omega\to0}$ are the same
  normalised law: the identity sees the **ratio** of the numerator's
  $\gamma$ coefficient to the denominator's, and the **exponent on $x$** — and
  §7.7 carries a break row for each of those two, moving it to 1.0 and 0.25.
  The $\omega$ coefficients are established on `B3.1`, by its own independent
  integration, not here.
- `B3.1`'s own caveat that its absolute burnout times are unreconstructable
  (Parker and Hottel's $K_c$ conversion is not printed there) does not touch
  this page: everything here is in $\tau$.
- `B3.2` — Szekely and Evans' grain model — was also read. Its local grain law
  is the same $m = 2/3$ shrinking core, and it reports its own reduction onto
  `B3.1` at 4.3e-4 through a full diffusion field. Nothing of `B3.2`'s is used
  numerically here; this page's kinetic-control assumption is exactly the
  limit in which `B3.2`'s intergrain diffusion disappears."""))

cells.append(code(r'''def yagi_kunii_eq6(x, omega, gamma):
    """B3.1's published equation, transcribed THERE from the Yagi & Kunii (1955)
    page: theta/theta_B against x = r/R, with omega = kc1/kd1, gamma = 3 kf1/kd1.
    Copied verbatim from pages/B3.1-shrinking-core; nothing else of B3.1 is used."""
    x = np.asarray(x, float)
    num = gamma*x + 3.0*omega*gamma*x**2 + omega*(1.0 - 2.0*gamma)*x**3
    return 1.0 - num/(omega + omega*gamma + gamma)

Xid = np.linspace(0, 1, 2001)
tau_over_3 = 1 - (1 - Xid)**(1/3)                 # Bhatia eq. (40), inverted
B31_IDENTITY = max(float(np.max(np.abs(yagi_kunii_eq6((1-Xid)**(1/3), 0.0, g) - tau_over_3)))
                   for g in (1e-6, 1e-3, 1.0, 1e3, 1e12))
B31_AT_OMEGA_1EM9 = float(np.max(np.abs(yagi_kunii_eq6((1-Xid)**(1/3), 1e-9, 1.0) - tau_over_3)))
B31_PUBLISHED = 3.84898779515197e-10               # pages/B3.1-shrinking-core/agreement.json
B31_RECONCILE = abs(B31_AT_OMEGA_1EM9 - B31_PUBLISHED)/B31_PUBLISHED

GAMMA_SWEEP = {g: float(np.max(np.abs(yagi_kunii_eq6((1-Xid)**(1/3), 0.0, g) - tau_over_3)))
               for g in (1e-6, 1e-3, 1.0, 1e3, 1e12)}

print("Bhatia & Perlmutter eq. (40)  vs  Yagi & Kunii eq. (6) at omega -> 0:")
print(f"   exact limit omega = 0, over gamma from 1e-6 to 1e12:   max |diff| = {B31_IDENTITY:.3e}")
print("   " + "  ".join(f"gamma={g:g}: {v:.1e}" for g, v in GAMMA_SWEEP.items()))
print("   THAT SWEEP CANNOT FAIL, and is printed only to show why: at omega = 0 the numerator is")
print("   gamma*x and the denominator is gamma, so GAMMA CANCELS IDENTICALLY and eq. (6) collapses")
print("   to 1 - x with no parameter left in it. Every gamma returns the same number by algebra,")
print("   not by robustness. The identity is blind to gamma and to all three omega coefficients;")
print("   what it CAN see is the ratio of the two gamma coefficients and the exponent on x, and")
print("   7.7 carries a break row for each.")
print(f"   at B3.1's own omega = 1e-9:                            max |diff| = {B31_AT_OMEGA_1EM9:.3e}")
print(f"   B3.1 PUBLISHES eq6_to_reaction_limit = {B31_PUBLISHED:.6e} for that same read;")
print(f"   this page reproduces it to {B31_RECONCILE:.2e} relative -- the 3.85e-10 is entirely")
print("   B3.1's omega = 1e-9 regularisation, and the two printed equations are identical.")
assert B31_RECONCILE < 1e-3, "B3.1's published metric must reproduce here"

fig, ax = plt.subplots(figsize=(5.4, 3.3))
ax.plot(Xid, 3*tau_over_3, lw=2.4, label=r"Bhatia eq. (40): $\tau = 3[1-(1-X)^{1/3}]$")
ax.plot(Xid[::40], 3*yagi_kunii_eq6((1-Xid[::40])**(1/3), 0.0, 1.0), "o", ms=4,
        label=r"B3.1 eq. (6), $\omega\to0$")
ax.set(xlabel="conversion $X$", ylabel=r"$\tau$"); ax.legend(fontsize=8)
fig.tight_layout(); plt.show()'''))

cells.append(md(r"""### 7.6 A reference table, written for reuse

The page's computed reference solutions, so an independent implementation can
regression-test against them. **These are model evaluations, not data.**"""))

cells.append(code(r'''ref_rows = []
for pv in (0.0, 1.0, 2.0, 5.0, 10.0):
    for tv in (0.1, 0.25, 0.5, 1.0, 1.5):
        Xv = X_of_tau(tv, pv)
        ref_rows.append(dict(psi=pv, tau=tv, X=Xv, S_over_So=SS0_of_X(Xv, pv),
                             X_M=(XM_eq33(pv) if pv >= 2 else np.nan),
                             m_star=m_star(pv), tau_half_sigma_inf=TAU_HALF_CLOSED(pv)))
ref = pd.DataFrame(ref_rows)
hdr = ("# Random pore model reference solutions COMPUTED BY PAGE B3.3 from the printed\n"
       "# equations (31), (32), (33), (37) of Bhatia & Perlmutter (1980) and from their\n"
       "# printed squared-deviation integral I.  NOT experimental data and NOT read from\n"
       "# any figure; see the .meta.yaml sidecar.  sigma -> infinity throughout.\n")
out = Path("data/random-pore-reference.csv")
if out.parent.is_dir():
    with out.open("w") as fh:
        fh.write(hdr); ref.to_csv(fh, index=False, float_format="%.12g")
    print(f"wrote {out} ({len(ref)} rows)")
display(ref.head(8))'''))

# ------------------------------------------------------------- break table
cells.append(md(r"""### 7.7 The defect-injection table, and the coverage map that builds itself

Every metric reported to `agreement.json` must have a row here that moves it,
or be named below with its reason and an above-floor companion. The coverage
map is assembled **from the table itself**, and the assert fails the notebook
on a mismatch in either direction.

**"It moved" is not enough, and this page learned that the hard way.** A row
whose injected defect shifts a 2e-16 metric to 2e-15 has changed float noise
into float noise, and an assert that tests `broken != baseline` passes on it.
The coverage assert therefore requires the move to clear a stated noise floor,
**row by row and not merely metric by metric** — a set-valued assert is
satisfied as soon as *one* row moves a metric, so a dead row on any of the
eight metrics that carry two rows would pass it unseen — and the smallest move
in the table is printed beside it. The row this rule
caught here — the wrong root of eq. (7) against `petersen_X_endpoint_dev`, a
check that could not fail because the endpoint collapses to 1 for *every* root
of eq. (7) — is reported in full below rather than quietly replaced, together
with the three other candidate defects that turned out not to work.

The defects are the ones a reader of this scan could actually commit: dropped
factors, mis-grouped expressions, the wrong root of a multi-root equation, a
forgotten unit conversion. **The list is not typed here — it is printed under
the table, generated from the rows themselves**, because the typed version
drifted: it named "Avrami's overlap law swapped for Petersen's" when no row
does that (both overlap rows perturb *Petersen's* eq. (4) and leave Avrami's
$1-e^{-V_E}$ untouched). A prose list of a table's contents is one more thing
that can contradict the table, so this page no longer keeps one.

One entry does need saying in advance, because it is a statement about which
rows *cannot* exist: `B3.1`'s eq. (6) is mis-transcribed in the two places the
$\omega\to0$ limit can still see — the exponent on its surviving $\gamma x$
term, and the ratio of its two $\gamma$ coefficients — and **not** in its
$3\omega\gamma$ coefficient, because §7.5 shows that limit is blind to all
three $\omega$ terms."""))

cells.append(code(r'''# ---- the injected defects, each recomputing the affected metrics -----------
psi_bad_noeps = 4*np.pi*P["fig4_Lo"]/P["fig4_So"]**2                 # (1-eps_o) dropped
psi_bad_2pi   = 2*np.pi*P["fig4_Lo"]*(1-P["fig4_eps_o"])/P["fig4_So"]**2
BRK_PSI_NOEPS = abs(psi_bad_noeps - P["fig4_psi"])/P["fig4_psi"]
BRK_PSI_2PI   = abs(psi_bad_2pi - P["fig4_psi"])/P["fig4_psi"]

# eq. (3): the OTHER root (r_p above the surface maximum)
b_f4 = np.sqrt(np.pi*P["fig4_Lo"]/3); rmax_f4 = 1/(2*b_f4)
rp_hi = brentq(lambda z: pet_S(z, P["fig4_Lo"]) - P["fig4_So"], rmax_f4, 1/b_f4*0.999,
               xtol=1e-18, rtol=8.9e-16)
BRK_PET_EPS_HIROOT = abs(pet_eps(rp_hi, P["fig4_Lo"]) - P["fig4_eps_o"])/P["fig4_eps_o"]
# eq. (4) misread as (1 - r/3 sqrt(pi L/3))
pet_eps_bad = lambda r_, L_: np.pi*r_**2*L_*(1 - (r_/3)*np.sqrt(np.pi*L_/3))
rp_f4 = SETS[0]["rp_from_S"]
BRK_PET_EPS_R3 = abs(pet_eps_bad(rp_f4, P["fig4_Lo"]) - P["fig4_eps_o"])/P["fig4_eps_o"]
BRK_PET_S_R3 = abs(pet_S(brentq(lambda z: pet_eps_bad(z, P["fig4_Lo"]) - P["fig4_eps_o"],
                                1e-12, rmax_f4, xtol=1e-18), P["fig4_Lo"])
                   - P["fig4_So"])/P["fig4_So"]
# curve-C feasibility with L read as 3.14e7
BRK_CURVE_C_L10 = P["fig7_C_So"]/pet_S(1/(2*np.sqrt(np.pi*P["fig7_C_Lo"]*10/3)), P["fig7_C_Lo"]*10)

# the wrong root of eq. (7): Petersen conversions collapse
G_WRONG = float(np.sort(np.roots([(4/27)*P["fig7_eps_o"], 0, -1, 1]).real)[1])
def X_pet_wrongG(s_, e_, G):  return e_/(1-e_)*((1+s_)**2*(G-1-s_)/(G-1) - 1)
BRK_PET_X_WRONGG = float(np.max(np.abs(X_pet_wrongG(sgrid, P["fig7_eps_o"], G_WRONG)
                                       - X_pet_via_eps(sgrid, P["fig7_eps_o"]))))
# The wrong root does NOT break the endpoint deviation, and this page reports that
# rather than tabulating it: X(s = 2G/3 - 1) = 1 is an exact identity for EVERY root
# of eq. (7).  The defect that does bite it is a mis-transcription of eq. (7) itself.
BRK_PET_ENDPOINT_WRONGG = abs(X_pet_wrongG(2*G_WRONG/3 - 1, P["fig7_eps_o"], G_WRONG) - 1.0)
G_2_27 = float(np.max(np.roots([(2/27)*P["fig7_eps_o"], 0, -1, 1]).real))   # (7) with 2/27
BRK_PET_ENDPOINT_227 = abs(X_pet_wrongG(2*G_2_27/3 - 1, P["fig7_eps_o"], G_2_27) - 1.0)
XB_wrong = lambda kt: X_pet_wrongG(kt/rpo_A, e_A, G_WRONG)
kts75 = np.linspace(1e-9, KT_75, 4001)
BRK_DAB_WRONGG = float(np.max(np.abs(XA(kts75) - XB_wrong(kts75))))

# eq. (32) with psi tau/2
X_bad32 = lambda t_, p_: 1 - np.exp(-t_*(1 + p_*t_/2))
BRK_MARCH_PSI2 = float(np.max(np.abs(march37(PSI_MARCH, TAU_END, 800)
                                     - X_bad32(TAU_END, PSI_MARCH))))
BRK_DAB_PSI2 = float(np.max(np.abs(X_bad32(kts75*So_A/(1-e_A), PSI_A) - XB(kts75))))
# implicit Euler: the march order collapses
e1 = float(np.max(np.abs(march37(PSI_MARCH, TAU_END, 400, theta=1.0) - X_of_tau(TAU_END, PSI_MARCH))))
e2 = float(np.max(np.abs(march37(PSI_MARCH, TAU_END, 800, theta=1.0) - X_of_tau(TAU_END, PSI_MARCH))))
BRK_MARCH_ORDER_EULER = float(np.log2(e1/e2)); BRK_MARCH_ERR_EULER = e2

# eq. (33) as (2-psi)/psi
XM_bad = lambda p_: 1 - np.exp((2-p_)/p_)
BRK_XM_TWO_ROUTES = max(abs(XM_rootfind(pv) - XM_bad(pv))/XM_bad(pv) for pv in (3.0, 5.0))
BRK_XM_SUP = abs((1 - np.exp(-1.0)) - P["XM_sup"])/P["XM_sup"]
# (a sign flip inside eq. (31)'s radical is NOT usable as a break row: the
#  flipped S/S_o is monotone in X and real only for X < 1 - exp(-1/psi), so it
#  has no maximum to root-find at all -- reported here rather than tabulated.)
BRK_XM_PSI1PCT = max(abs(XM_rootfind(pv*1.01) - XM_eq33(pv))/XM_eq33(pv) for pv in (3.0, 5.0))

# I with the (1-X)^1 factor dropped from the model term
def I_bad(psi_, m_, n=400):
    x, wt = _gl(n); v = 30*(x+1); W = 30*wt
    return float(np.sum(W*(np.sqrt(1 + psi_*v) - np.exp(-m_*v))**2*np.exp(-v)))
BRK_I_TWO_ROUTES = max(abs(I_bad(pv, mm_) - I_closed(pv, mm_))/I_closed(pv, mm_)
                       for pv, mm_ in ((1.0, 2/3), (2.0, 0.49)))
def m_star_bad(psi_):
    f = lambda z: (I_bad(psi_, z+1e-6) - I_bad(psi_, z-1e-6))/2e-6
    return brentq(f, -0.45, 3.0, xtol=1e-13, rtol=8.9e-16)
BRK_M_PSI1 = m_star_bad(1.0); BRK_M_PSI2 = m_star_bad(P["m_psi_window_hi"]); BRK_M_PSI0 = m_star_bad(0.0)
BRK_M_REL23 = abs((BRK_M_PSI1 - P["grain_m_sphere"])/P["grain_m_sphere"])
BRK_M_TWO_ROUTES = max(abs(m_star_bad(pv) - m_star_brent(pv)) for pv in (1.0, 2.0))
BRK_NRMS_M1 = np.sqrt(max(I_bad(1.0, 1.0), 0)/I_energy(1.0))
BRK_M1_PENALTY = BRK_NRMS_M1/np.sqrt(max(I_bad(1.0, BRK_M_PSI1), 1e-30)/I_energy(1.0))

# tau_1/2: read off a 200-point sweep instead of root-found; sigma factor dropped
def tau_half_swept(pv, sg, n=200):
    g = np.linspace(1e-6, sg*(1-1e-9), n)
    return g[np.argmin(np.abs(X_of_tau(g, pv, sg) - 0.5))]
BRK_TAU_BLOCK_SWEEP = max(abs(tau_half_swept(pv, sg) - tau_at_X_block(0.5, [pv], sg)[0])
                          / tau_at_X_block(0.5, [pv], sg)[0]
                          for pv in psis_b for sg in (0.25, 1.0, 10.0, 100.0))
t100_nosig = np.array([TAU_HALF_CLOSED(pv) for pv in PSI_FIG2])       # sigma factor dropped
BRK_SIGMA100 = float(np.max(np.abs(t100_nosig - tinf)/tinf))
t025_bad = tau_at_X_block(0.5, PSI_FIG2, 1.0)                          # sigma read as 1, not 0.25
BRK_SIGMA025 = float((t025_bad.max() - t025_bad.min())/t025_bad.min())
BRK_TAU_CLOSED = max(abs(tau_at_X_block(0.5, [pv], 1e9)[0]
                         - (np.log(2) if pv == 0 else (-1+np.sqrt(1+pv*np.log(2)*2))*2/pv))
                     / TAU_HALF_CLOSED(pv) for pv in psis_b)           # ln2 -> 2 ln2 in the closed form

# eq. (35) with the printed group misread as 0.68 (0.67 is a two-decimal number
# in a bilevel scan).  NOTE: inverting eq. (35)'s log argument is NOT usable as a
# break row -- the inverted equation has no root in (0, eps) at all, so it is
# reported here rather than tabulated.
GHAT_BAD = 0.68
stat35_bad = lambda e0: GHAT_BAD/(1-e0) - 1 - np.log((1-EPS_F)/(1-e0))
BRK_EPS_O_ROOT = brentq(stat35_bad, 1e-12, EPS_F-1e-12, xtol=1e-17, rtol=8.9e-16)
BRK_STAT35 = stat35_bad(P["opt_eps_o"])
BRK_EPS_TWO_ROUTES = abs(BRK_EPS_O_ROOT - EPS_O_DIRECT)
BRK_GHAT_066 = brentq(lambda e0: 0.66/(1-e0) - 1 - np.log((1-EPS_F)/(1-e0)),
                      1e-12, EPS_F-1e-12, xtol=1e-17, rtol=8.9e-16)
BRK_S_DEFICIT = abs(1 - S_of_eps_o(P["opt_eps_o"], EPS_F, GHAT)/S_of_eps_o(BRK_EPS_O_ROOT, EPS_F, GHAT))
# Petersen's eq. (4) misread as (1 - (r/3) sqrt(pi L/3)): the porosity law
# becomes 3u^2 - u^3.  (Halving eq. (20)'s exponent instead is NOT usable here:
# it puts Avrami below Petersen everywhere, so there is no crossing to find.)
gap_bad = lambda z: (1 - np.exp(-3*z**2)) - (3*z**2 - z**3)
BRK_OVERLAP_GAP = float(-minimize_scalar(lambda z: -gap_bad(z), bounds=(0.02, 0.6),
                                         method="bounded").fun)
BRK_OVERLAP_CROSS = brentq(gap_bad, 0.05, 0.9, xtol=1e-17, rtol=8.9e-16)

# B3.1's eq. (6) mis-transcribed. NOT the 3 omega gamma coefficient: at omega = 0
# every omega term drops out, so the reaction-control identity is BLIND to all
# three of them (stated in 7.5, and exactly why B3.1 measures its own limit at
# omega = 1e-9). The defect that does bite is an EXPONENT slip in the numerator's
# first term -- the documented hazard of that 1955 scan, whose text layer drops
# eq. (6)'s exponents (B3.1's README).
def yk_bad(x, omega, gamma):
    x = np.asarray(x, float)
    num = gamma*x**2 + 3.0*omega*gamma*x**2 + omega*(1.0 - 2.0*gamma)*x**3
    return 1.0 - num/(omega + omega*gamma + gamma)
BRK_B31_IDENTITY = max(float(np.max(np.abs(yk_bad((1-Xid)**(1/3), 0.0, g) - tau_over_3)))
                       for g in (1e-3, 1.0, 1e3))
BRK_B31_RECONCILE = abs(float(np.max(np.abs(yk_bad((1-Xid)**(1/3), 1e-9, 1.0) - tau_over_3)))
                        - B31_PUBLISHED)/B31_PUBLISHED
# The identity is blind to gamma's VALUE (it cancels) but not to the RATIO of the
# numerator's gamma coefficient to the denominator's; that is the second axis it
# can see, so it gets its own row rather than being left as a prose claim.
def yk_bad_ratio(x, omega, gamma):
    x = np.asarray(x, float)
    num = 2.0*gamma*x + 3.0*omega*gamma*x**2 + omega*(1.0 - 2.0*gamma)*x**3
    return 1.0 - num/(omega + omega*gamma + gamma)
BRK_B31_GAMMA_RATIO = max(float(np.max(np.abs(yk_bad_ratio((1-Xid)**(1/3), 0.0, g) - tau_over_3)))
                          for g in (1e-6, 1e-3, 1.0, 1e3, 1e12))

# Hashimoto: S_o left in m2/cm3 (the unit conversion forgotten) ; psi swapped
BRK_XM_VC = XM_eq33(P["hashimoto_psi_CS"])                             # the two chars swapped
BRK_RBAR_VC = 2/(P["hashimoto_psi_VC"]*P["hashimoto_So"])*1e7          # cm2/cm3 conversion dropped
BRK_INTERCEPT = P["hashimoto_So"]**2/1e12
BRK_PHI098 = brentq(lambda z: 3*(z/np.tanh(z) - 1)/z**2 - 0.9, 1e-6, 10.0)  # eta read as 0.90

def I_half(psi_, m_, n=200):             # the misfit integral truncated to X <= 0.5
    x, wt = _gl(n); w = 0.25*(x+1) + 0.5; W = 0.25*wt
    return float(np.sum(W*(w*np.sqrt(1 - psi_*np.log(w)) - w**m_)**2))
PSI_M0_HALF = brentq(lambda pv: (I_half(pv, 1e-6) - I_half(pv, -1e-6))/2e-6, 1.0, 10.0,
                     xtol=1e-12, rtol=8.9e-16)

G_WRONG_F4 = float(np.sort(np.roots([(4/27)*P["fig4_eps_o"], 0, -1, 1]).real)[1])
def SS0_pet_via_eq6_wrongG(X_, e_, G):
    u0 = u_of_eps(e_); s_ = u_of_eps(e_ + X_*(1-e_))/u0 - 1
    return (2*(1+s_)*(G-1-s_) - (1+s_)**2)/(2*G - 3)
BRK_PET_SS_WRONGG = max(abs(SS0_pet_via_eq6_wrongG(x, P["fig4_eps_o"], G_WRONG_F4)
                            - float(SS0_petersen(np.array([x]), P["fig4_eps_o"])[0]))
                        for x in (0.25, 0.5, 0.75))

BREAKS = [
 ("psi_fig4_vs_printed_rel", "psi with (1-eps_o) dropped", PSI_F4_REL, BRK_PSI_NOEPS),
 ("psi_fig4_vs_printed_rel", "psi with 2 pi for 4 pi", PSI_F4_REL, BRK_PSI_2PI),
 ("petersen_eps_rel_fig4set", "eq. (3): the root above the surface maximum", PET_EPS_REL, BRK_PET_EPS_HIROOT),
 ("petersen_eps_rel_fig4set", "eq. (4) read as (1 - (r/3) sqrt(pi L/3))", PET_EPS_REL, BRK_PET_EPS_R3),
 ("petersen_S_rel_fig4set", "eq. (4) read as (1 - (r/3) sqrt(pi L/3))", PET_S_REL, BRK_PET_S_R3),
 ("petersen_eps_rel_fig7AB", "eq. (3): the root above the surface maximum",
  PET_AB_EPS_REL, abs(pet_eps(brentq(lambda z: pet_S(z, P["fig7_AB_Lo"]) - P["fig7_AB_So"],
                                     1/(2*np.sqrt(np.pi*P["fig7_AB_Lo"]/3)),
                                     0.999/np.sqrt(np.pi*P["fig7_AB_Lo"]/3), xtol=1e-18),
                             P["fig7_AB_Lo"]) - P["fig7_eps_o"])/P["fig7_eps_o"]),
 ("curveC_So_over_Smax", "L_o read as 3.14e7 (one decade out)", CURVE_C_INFEAS, BRK_CURVE_C_L10),
 ("curveD_L_over_LAB", "S_o read as the curve-A value instead of curve C's", L_D_RATIO,
  P["fig7_AB_So"]**2/(4*u_of_eps(P["fig7_eps_o"])**2*(1-u_of_eps(P["fig7_eps_o"]))**2*3*np.pi)
  / P["fig7_AB_Lo"]),
 ("curveD_roundtrip_rel", "curve-D length perturbed 1 %", L_D_ROUNDTRIP,
  abs(pet_S(rp_D, L_D*1.01) - P["fig7_C_So"])/P["fig7_C_So"]),
 ("XM_sup_vs_printed_rel", "eq. (33) as (2-psi)/psi", abs(XM_SUP - P["XM_sup"])/P["XM_sup"], BRK_XM_SUP),
 ("XM_two_routes_rel", "eq. (33) as (2-psi)/psi", XM_TWO_ROUTES, BRK_XM_TWO_ROUTES),
 ("XM_two_routes_rel", "psi 1 % adrift in the root-find route only", XM_TWO_ROUTES, BRK_XM_PSI1PCT),
 ("march_cn_vs_eq32_max", "eq. (32) with psi tau/2 for psi tau/4", MARCH_ERR_800, BRK_MARCH_PSI2),
 ("march_cn_vs_eq32_max", "implicit Euler instead of Crank-Nicolson", MARCH_ERR_800, BRK_MARCH_ERR_EULER),
 ("march_observed_order", "implicit Euler instead of Crank-Nicolson", MARCH_ORDER, BRK_MARCH_ORDER_EULER),
 ("fig7_max_dX_upto_75", "eq. (32) with psi tau/2 for psi tau/4", DAB_UPTO75, BRK_DAB_PSI2),
 ("fig7_max_dX_upto_75", "the wrong root of eq. (7)", DAB_UPTO75, BRK_DAB_WRONGG),
 ("fig7_X_at_two_points", "the wrong root of eq. (7)", X_AT_2PT,
  XA(brentq(lambda z: abs(XA(z) - XB_wrong(z)) - 0.02, 1e-7, KT_END*0.5, xtol=1e-20))),
 ("petersen_X_two_routes", "the wrong root of eq. (7)", PET_X_TWO_ROUTES, BRK_PET_X_WRONGG),
 ("petersen_X_endpoint_dev", "eq. (7)'s 4/27 read as 2/27", PET_X_ENDPOINT, BRK_PET_ENDPOINT_227),
 ("petersen_SS_two_routes", "eq. (7)'s wrong root enters the d(eq.6)/ds route",
  PET_SS_TWO_ROUTES, BRK_PET_SS_WRONGG),
 ("m_star_psi0", "I with the (1-X)^1 factor dropped", M_PSI0, BRK_M_PSI0),
 ("m_star_psi1", "I with the (1-X)^1 factor dropped", M_PSI1, BRK_M_PSI1),
 ("m_star_psi2", "I with the (1-X)^1 factor dropped", M_PSI2, BRK_M_PSI2),
 ("m_star_psi1_vs_two_thirds_rel", "I with the (1-X)^1 factor dropped", abs(M_PSI1_REL_23), BRK_M_REL23),
 ("m_star_two_routes", "I with the (1-X)^1 factor dropped in one route only",
  M_TWO_ROUTES, BRK_M_TWO_ROUTES),
 ("I_two_routes_rel", "I with the (1-X)^1 factor dropped in the quadrature", I_TWO_ROUTES, BRK_I_TWO_ROUTES),
 ("nrms_m1_at_psi1", "I with the (1-X)^1 factor dropped", NRMS_M1_PSI1, BRK_NRMS_M1),
 ("m1_penalty_at_psi1", "I with the (1-X)^1 factor dropped", M1_PENALTY, BRK_M1_PENALTY),
 ("quad_rel_err_400", "quadrature reduced to 25 nodes", QUAD_ERR_400, rows[0][1]),
 ("tau_block_vs_scipy_rel", "tau_1/2 read off a 200-point sweep", TAU_BLOCK_VS_SCIPY, BRK_TAU_BLOCK_SWEEP),
 ("tau_closed_vs_block_rel", "closed form with 2 ln2 for ln2", TAU_HALF_CLOSED_VS_BLOCK, BRK_TAU_CLOSED),
 ("sigma100_max_rel", "the (1-tau/sigma)^3 factor dropped at sigma = 100", SIGMA100_MAX_REL, BRK_SIGMA100),
 ("sigma025_psi_spread", "sigma read as 1 instead of 0.25", SIGMA025_SPREAD, BRK_SIGMA025),
 ("eq35_root_eps_o", "the printed group read as 0.68", EPS_O_ROOT, BRK_EPS_O_ROOT),
 ("eq35_root_eps_o", "the printed group read as 0.66", EPS_O_ROOT, BRK_GHAT_066),
 ("eq35_residual_at_printed", "the printed group read as 0.68", abs(STAT35_AT_PRINTED), abs(BRK_STAT35)),
 ("eq35_two_routes", "the printed group read as 0.68 in one route only", EPS_O_TWO_ROUTES, BRK_EPS_TWO_ROUTES),
 ("eq35_S_deficit_at_printed", "the printed group read as 0.68", S_DEFICIT_AT_PRINTED, BRK_S_DEFICIT),
 ("overlap_max_gap", "eq. (4) read as (1 - (r/3) sqrt(pi L/3))", OVERLAP_MAX_GAP, BRK_OVERLAP_GAP),
 ("overlap_cross_u", "eq. (4) read as (1 - (r/3) sqrt(pi L/3))", U_CROSS, BRK_OVERLAP_CROSS),
 ("b31_eq40_identity", "B3.1's eq. (6) numerator: gamma x read as gamma x^2", B31_IDENTITY, BRK_B31_IDENTITY),
 ("b31_eq40_identity", "B3.1's eq. (6): the numerator's gamma coefficient doubled",
  B31_IDENTITY, BRK_B31_GAMMA_RATIO),
 ("b31_published_reconcile", "B3.1's eq. (6) numerator: gamma x read as gamma x^2",
  B31_RECONCILE, BRK_B31_RECONCILE),
 ("hashimoto_XM_VC", "the two chars' psi swapped", XM_VC, BRK_XM_VC),
 ("hashimoto_rbar_VC_nm", "the m2 -> cm2 conversion dropped", RBAR_VC, BRK_RBAR_VC),
 ("hashimoto_intercept_e12", "the m2 -> cm2 conversion dropped", INTERCEPT/1e12, BRK_INTERCEPT),
 ("eta098_thiele_sphere", "eta read as 0.90 instead of 0.98", PHI_098, BRK_PHI098),
 ("psi_m_star_zero", "I truncated to X <= 0.5 (half the conversion range)", PSI_M0, PSI_M0_HALF),
 ("nrms_at_psi_m_star_zero", "I truncated to X <= 0.5 (half the conversion range)",
  NRMS_AT_PSI_M0, I_nrms(PSI_M0_HALF, 0.0)),
 ("mstar_node_sensitivity", "quadrature reduced to 25 nodes in dI/dm", MSTAR_NODE_SENS,
  abs(brentq(lambda z: -2*dJ_db(z+1, 1.0, 25) - 2/(2*z+1)**2, -0.45, 3.0, xtol=1e-15) - M_PSI1)),
 ("grain_branch_gap_at_m1", "eq. (39)'s general form taken at m = 0.999, the limit unconverged",
  BRANCH_GAP, abs(X_grain(0.4, 1.0) - (1 - (1 - (1-0.999)*0.4)**(1/(1-0.999))))),
 ("grain_branch_gap_at_m1", "the m = 1 branch written as 1 - exp(-2 tau)",
  BRANCH_GAP, abs((1 - np.exp(-2*0.4)) - (1 - (1-(1-(1-1e-9))*0.4)**(1/(1-(1-1e-9)))))),
]
bt = pd.DataFrame(BREAKS, columns=["metric", "defect injected", "baseline", "broken"])
with pd.option_context("display.float_format", lambda v: f"{v:.4e}"):
    display(bt)

# The distinct defect classes, GENERATED FROM THE ROWS rather than typed in the
# markdown above. An earlier version of this page listed the defects in prose and
# the list named one -- "Avrami's overlap law swapped for Petersen's" -- that no
# row injects (both overlap rows perturb PETERSEN's law and leave Avrami's alone).
# Printing the list from the table is the only version that cannot drift.
_dcount = bt["defect injected"].value_counts()
print(f"\n{len(_dcount)} distinct injected defects over {len(bt)} rows"
      f" (generated from the table, not typed):")
for _d in sorted(_dcount.index):
    print(f"  [{_dcount[_d]} row{'s' if _dcount[_d] > 1 else ' '}]  {_d}")

STRUCTURAL = {
    "sym_zero_identity_count": "a count of sympy identities; it verifies the paper's algebra, not "
        "this page's numerics, and no runtime defect can move it -- cited as evidence for nothing numerical",
}
# Companions for the metrics that sit under check_agreement.py's ABS_FLOOR while
# healthy. Which metrics those ARE is decided from the measured values in the
# agreement cell, not typed here, so the list cannot drift; and every number a
# companion note quotes is INTERPOLATED from the row it refers to, so the note
# cannot contradict the table (it did, before this page was verified).
COMPANION = {
    "curveD_roundtrip_rel":   "curveD_L_over_LAB (CI-active) and the 1 %-perturbed-length break row",
    "petersen_X_two_routes":  f"the wrong-root-of-(7) break row (it moves this to {BRK_PET_X_WRONGG:.2e})",
    "petersen_X_endpoint_dev": f"the 4/27-read-as-2/27 break row (it moves this to "
                               f"{BRK_PET_ENDPOINT_227:.2e}); the WRONG-ROOT row cannot move it -- "
                               f"see the note printed under the break table",
    "petersen_SS_two_routes": f"the wrong-root-of-(7) break row (it moves this to {BRK_PET_SS_WRONGG:.2e})",
    "tau_block_vs_scipy_rel": "the 200-point-sweep break row",
    "b31_eq40_identity":      "b31_published_reconcile (CI-active) and its two break rows "
                              "(the exponent on x; the ratio of the two gamma coefficients)",
    "XM_two_routes_rel":      "its two break rows (eq. 33 mis-grouped; psi 1 % adrift in one route)",
    "quad_rel_err_400":       "the 25-node break row and mstar_node_sensitivity",
    "I_two_routes_rel":       "the dropped-(1-X) break row",
    "eq35_two_routes":        "eq35_root_eps_o and eq35_residual_at_printed (both CI-active)",
    "grain_branch_gap_at_m1": "its two break rows (the m = 1 branch written as 1 - exp(-2 tau); "
                              "the general form taken at m = 0.999)",
    "tau_closed_vs_block_rel":"the 2-ln2 break row",
    "mstar_node_sensitivity": "the 25-node break row",
    "m_star_two_routes":      "m_star_psi1 (CI-active) and the dropped-(1-X) break row",
}
print("Structural (no break row can exist):")
for k, v in STRUCTURAL.items():
    print(f"  {k}: {v}")

# --- the check that could not fail, kept as the page's own teaching exhibit ----
print("\nA ROW THAT DOES NOT WORK, AND WHY IT IS NOT IN THE TABLE. The obvious defect for")
print("petersen_X_endpoint_dev is the wrong root of eq. (7) -- and it CANNOT break it:")
print(f"   baseline {PET_X_ENDPOINT:.4e}  ->  wrong root of eq. (7)  {BRK_PET_ENDPOINT_WRONGG:.4e}")
print("   both float noise, so an assert that only demands 'the value moved' would pass on noise.")
print("   X(s = 2G/3 - 1) = 1 is an EXACT IDENTITY FOR EVERY ROOT of eq. (7): with (1+s) = 2G/3")
print("   and (G-1-s) = G/3 the bracket of eq. (6) is (4G^3/27)/(G-1) - 1, and eq. (7) says")
print("   (4/27) eps G^3 = G - 1, so the bracket is 1/eps - 1 and X = 1 identically. The metric")
print("   therefore establishes that eqs. (6) and (7) are mutually consistent AS TRANSCRIBED, and")
print("   says NOTHING about which root was selected (that is petersen_X_two_routes' job). The row")
print(f"   it does get mis-transcribes eq. (7) itself, 4/27 as 2/27: {BRK_PET_ENDPOINT_227:.4e}, which is")
print(f"   the closed form 1/(1-eps_o) = {1/(1-P['fig7_eps_o']):.4e} for ANY root of the misread cubic.")

# --- candidate break rows that do NOT work, each with its measured witness -----
def _d_flip(w_, pv):                       # eq. (31) with the sign in the radical flipped
    f = lambda z: z*np.sqrt(1 + pv*np.log(z))
    return (f(w_+1e-9) - f(w_-1e-9))/2e-9
FLIP_LO = np.exp(-1/5.0)
FLIP_ENDS = (_d_flip(FLIP_LO + 1e-6, 5.0), _d_flip(1 - 1e-9, 5.0))
stat35_inv = lambda e0: GHAT/(1-e0) - 1 - np.log((1-e0)/(1-EPS_F))     # (35)'s log inverted
INV35_ENDS = (stat35_inv(1e-12), stat35_inv(EPS_F - 1e-12))
gap_half = lambda z: (1 - np.exp(-1.5*z**2)) - (3*z**2 - 2*z**3)       # (20)'s exponent halved
_rh = minimize_scalar(lambda z: -gap_half(z), bounds=(1e-9, 1.0), method="bounded",
                      options=dict(xatol=1e-15))
HALF20_MAX, HALF20_AT = -float(_rh.fun), float(_rh.x)

print("\nTHREE FURTHER CANDIDATE ROWS WERE TRIED AND DO NOT WORK. Reported, not dropped:")
print(f" 1. eq. (31)'s radical sign flipped: (1-X) sqrt(1 + psi ln(1-X)) is real only for")
print(f"    w > exp(-1/psi) = {FLIP_LO:.6f} at psi = 5, and strictly increasing there --")
print(f"    dS/dw = {FLIP_ENDS[0]:+.4e} at the lower end and {FLIP_ENDS[1]:+.4e} at w = 1, same sign,")
print("    so there is no maximum to root-find and XM_two_routes_rel has nothing to compare.")
print(f" 2. eq. (35)'s log argument inverted: the residual is {INV35_ENDS[0]:+.6f} at eps_o -> 0 and")
print(f"    {INV35_ENDS[1]:+.6f} at eps_o -> eps, the same sign, so the equation has NO ROOT in (0, eps).")
print(f" 3. eq. (20)'s exponent halved (V = 1 - exp(-V_E/2)): Avrami then sits below Petersen for")
print(f"    every u > 0 -- the maximum of the gap is {HALF20_MAX:+.2e} at u = {HALF20_AT:.2e}, i.e. the")
print(f"    origin, and it is {gap_half(1.0):+.6f} at u = 1 -- so overlap_cross_u has no crossing.")
print("\n(the below-ABS_FLOOR set is decided from the measured values in the agreement cell)")'''))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**Nothing to the model, and the page says so plainly.** The random pore model
is two closed forms; most of this notebook would run with pymrm uninstalled,
exactly as on `J1.4`, `A1.6` and `B1.7`. What `newton` and `NumJac` do here is
narrow and real:

1. **Whole families of inversions as one Newton problem.** Eq. (32) inverted
   for $\tau$ across a $(\psi,\sigma)$ block, and eq. (35) inverted for
   $\epsilon_o$, are solved with `NumJac((K, 1))` — independent scalars, a
   diagonal Jacobian by construction. The shape is `(K, 1)` and never `(K,)`:
   the bare 1-D shape would declare every state coupled to every other and
   build a dense $K\times K$ Jacobian for no change in the answer.
2. **A residual written so it can be solved at all.** Eq. (32)'s direct
   residual underflows at $\psi = 100$ and hands `newton` an exactly singular
   Jacobian; the log form is the same equation with the exponential taken out,
   and the sigmoid keeps $\tau$ inside $(0,\sigma)$. Both are the kind of
   reformulation a block solver forces you to get right once instead of
   per-point.
3. **The rate law marched as an initial-value problem** (`march37`), so the
   paper's eq. (37) and eq. (32) check each other rather than being asserted
   equal — Crank–Nicolson, observed order 2.000, implicit Euler kept as the
   break row.
4. **Root-finds and closed forms where a sweep would have been wrong.** Every
   extremum on this page is root-found: $X_M$, the best-fit $m^\star$, the
   optimal $\epsilon_o$, the divergence conversion of §6.3, the crossing of
   the two overlap laws, the $\psi$ at which $m^\star$ reaches zero. §7.7's
   break table shows what a 200-point sweep would have cost on $\tau_{1/2}$.

The one output useful beyond reproduction is not pymrm at all: the
**normalised misfit** of §6.5. It turns "which grain shape factor fits this
pore structure" into a curve with an error bar attached, and it says something
the paper's Figure 6 cannot — that past $\psi = 2$ the best available $m$ is
already misfitting by more than 8 %, rising through 20 % at $\psi = 5$,
because no order-of-reaction model can produce a rate maximum at all."""))

# ------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**Use this page for** the random pore model itself: `X_of_tau` (eq. 32),
`SS0_of_X` (eq. 31), `XM_eq33` and `psi_of` are the whole kinetic-control
model, and `rate37` drops into any solid-conversion source term where the
surface area evolves. `m_star(psi)` gives the grain shape factor equivalent to
a given pore structure, with `I_nrms(psi, m)` saying how badly that equivalence
misfits — which is the number to check before adopting an
order-of-reaction model for a porous solid. `u_of_eps`, `G_of_eps`,
`X_petersen` and `SS0_petersen` are Petersen (1957) as this paper reprints it,
with the root of eq. (7) selected by construction rather than by luck.

**The traps, in order of how much they cost:**

- **Eq. (7) has three real roots and the paper does not say which.** For
  $\epsilon_o = 0.3$ they are $-5.181$, $1.0517$ and $4.1293$; only
  $G = 3/(2u_o) = 4.1293$ is consistent with eqs. (3), (4) and (6). Picking
  the middle root puts complete gasification at a negative time and moves the
  §6.3 comparison from 0.0218 to 36.4 in conversion — a factor **1672**, the
  5th largest *absolute* shift among §7.7's 53 rows but only the 22nd largest
  *relative* one. §7.7 prints the factor and both ranks; a table this wide has
  no single "largest row". What it does **not** move is the endpoint identity
  $X(s = 2G/3-1) = 1$, which holds for *every* root of eq. (7): if you are
  checking your root selection, that is not the check to use. §7.7 says why and
  what to use instead.
- **$\psi$ carries $(1-\epsilon_o)$ and eq. (32) carries $\psi\tau/4$**, not
  $\psi\tau/2$. Both are the kind of thing a scan makes easy to lose, and both
  have break rows.
- **Eq. (3) is not monotone in $r_p$.** It has a maximum, so it has two roots
  for any attainable $S$; the physical one is below
  $r_p\sqrt{\pi L/3} = 1/2$. And it has *no* root when the target surface
  exceeds $S_{\max}(L)$ — which is exactly the situation of Figure 7's curve C.
- **Eq. (32) cannot be inverted naively at large $\psi$**: the exponential
  underflows. Use the log residual.

**Do not use this page as evidence about any real solid.** Its only
experimental content is a figure it deliberately does not touch, and it says
so in §4, in the title cell and here. The three char numbers in §6.6 are
transcribed prose, and everything derived from them is a *prediction* of the
model at those parameters, never a validation of it. Anyone with the Hashimoto
et al. (1979) tables can build the comparison this page leaves open.

**Do not extrapolate past kinetic control.** The whole page assumes negligible
diffusional resistance, which is what Part II of this paper (not on disk here)
lifts. `B3.2` is the neighbouring page where intergrain diffusion is present
and matters."""))

# ------------------------------------------------------------- agreement
cells.append(code(r'''metrics = {
    "sym_zero_identity_count":        float(SYM_ZERO_COUNT),
    "psi_fig4_vs_printed_rel":        float(PSI_F4_REL),
    "petersen_eps_rel_fig4set":       float(PET_EPS_REL),
    "petersen_S_rel_fig4set":         float(PET_S_REL),
    "petersen_eps_rel_fig7AB":        float(PET_AB_EPS_REL),
    "curveC_So_over_Smax":            float(CURVE_C_INFEAS),
    "curveD_L_over_LAB":              float(L_D_RATIO),
    "curveD_roundtrip_rel":           float(L_D_ROUNDTRIP),
    "XM_sup_vs_printed_rel":          float(abs(XM_SUP - P["XM_sup"])/P["XM_sup"]),
    "XM_two_routes_rel":              float(XM_TWO_ROUTES),
    "march_cn_vs_eq32_max":           float(MARCH_ERR_800),
    "march_observed_order":           float(MARCH_ORDER),
    "fig7_max_dX_upto_75":            float(DAB_UPTO75),
    "fig7_X_at_two_points":           float(X_AT_2PT),
    "petersen_X_two_routes":          float(PET_X_TWO_ROUTES),
    "petersen_X_endpoint_dev":        float(PET_X_ENDPOINT),
    "petersen_SS_two_routes":         float(PET_SS_TWO_ROUTES),
    "m_star_psi0":                    float(M_PSI0),
    "m_star_psi1":                    float(M_PSI1),
    "m_star_psi2":                    float(M_PSI2),
    "m_star_psi1_vs_two_thirds_rel":  float(abs(M_PSI1_REL_23)),
    "m_star_two_routes":              float(M_TWO_ROUTES),
    "I_two_routes_rel":               float(I_TWO_ROUTES),
    "nrms_m1_at_psi1":                float(NRMS_M1_PSI1),
    "m1_penalty_at_psi1":             float(M1_PENALTY),
    "quad_rel_err_400":               float(QUAD_ERR_400),
    "mstar_node_sensitivity":         float(MSTAR_NODE_SENS),
    "tau_block_vs_scipy_rel":         float(TAU_BLOCK_VS_SCIPY),
    "tau_closed_vs_block_rel":        float(TAU_HALF_CLOSED_VS_BLOCK),
    "sigma100_max_rel":               float(SIGMA100_MAX_REL),
    "sigma025_psi_spread":            float(SIGMA025_SPREAD),
    "eq35_root_eps_o":                float(EPS_O_ROOT),
    "eq35_residual_at_printed":       float(abs(STAT35_AT_PRINTED)),
    "eq35_two_routes":                float(EPS_O_TWO_ROUTES),
    "eq35_S_deficit_at_printed":      float(S_DEFICIT_AT_PRINTED),
    "overlap_max_gap":                float(OVERLAP_MAX_GAP),
    "overlap_cross_u":                float(U_CROSS),
    "b31_eq40_identity":              float(B31_IDENTITY),
    "b31_published_reconcile":        float(B31_RECONCILE),
    "hashimoto_XM_VC":                float(XM_VC),
    "hashimoto_rbar_VC_nm":           float(RBAR_VC),
    "hashimoto_intercept_e12":        float(INTERCEPT/1e12),
    "eta098_thiele_sphere":           float(PHI_098),
    "psi_m_star_zero":                float(PSI_M0),
    "nrms_at_psi_m_star_zero":        float(NRMS_AT_PSI_M0),
    "grain_branch_gap_at_m1":         float(BRANCH_GAP),
}

ABS_FLOOR = 1e-12          # check_agreement.py's comparison floor
below = {k for k, v in metrics.items() if abs(v) < ABS_FLOOR}
missing = below - set(COMPANION)
assert not missing, f"below-floor metrics with no named companion: {missing}"
print(f"{len(below)} of {len(metrics)} metrics sit below check_agreement.py's ABS_FLOOR = {ABS_FLOOR:g}")
print("while healthy -- outside the regression suite, protected by their companions and not by CI:")
for k in sorted(below):
    print(f"  {k} = {metrics[k]:.2e}   companion = {COMPANION[k]}")
print()

covered = set(bt["metric"])
declared = covered | set(STRUCTURAL)
assert set(metrics) == declared, (
    f"coverage mismatch: unclaimed={set(metrics)-declared}, phantom={declared-set(metrics)}")

# A move must clear a stated noise floor, not merely be non-zero: a row that takes
# a 2e-16 metric to 2e-15 has replaced float noise with float noise and proves
# nothing. (Exactly that row was live on this page before verification.)
MOVE_FLOOR = 1e-9
shift = (bt.broken - bt.baseline).abs()
moved = bt[shift > MOVE_FLOOR]
assert set(moved["metric"]) == covered, (
    "a break row failed to move its metric outside float noise: "
    f"{covered - set(moved['metric'])}")
# ...and ROW BY ROW, not only metric by metric. The set assert above is satisfied
# as soon as ONE row moves a metric, so a dead row on any of the metrics carrying
# two rows would pass it silently. There is no such row here -- the smallest move
# is 5e-05 -- but the hole is closed structurally rather than by luck.
assert (shift > MOVE_FLOOR).all(), (
    "break rows that do not move their own metric outside float noise: "
    f"{list(zip(bt.loc[shift <= MOVE_FLOOR, 'metric'], bt.loc[shift <= MOVE_FLOOR, 'defect injected']))}")
rel = shift/bt.baseline.abs().clip(lower=1e-300)
w = bt.loc[rel.idxmin()]
print(f"coverage: {len(metrics)} metrics = {len(covered)} with moving break rows "
      f"+ {len(STRUCTURAL)} structural, over {len(bt)} rows.  Assert passed in both directions,")
print(f"  and all {len(bt)} rows clear the floor individually, so no dead row can hide behind a")
print(f"  live one on the {int((bt['metric'].value_counts() > 1).sum())} metrics that carry more than one row.")
print(f"  smallest move in the table: {shift.min():.2e} absolute (floor {MOVE_FLOOR:g}); a row that")
print("  only moves its metric inside float noise is not a break row and fails this assert.")
print(f"  WEAKEST ROW IN RELATIVE TERMS: {w['metric']} moves {rel.min():.1%}")
print(f"  under '{w['defect injected']}'")
print("  -- BELOW check_agreement.py's REL_TOL = 5 %, so that injected defect would pass the")
print("  regression suite. It is the one row here CI could not catch.\n")

# The wrong-root row is the one the Reuse section and the metadata quote, so it is
# RANKED here on both readings instead of being called "one of the largest": a
# reader who checks the table must find the claim true whichever way they rank.
_wr = int(bt.index[(bt["metric"] == "fig7_max_dX_upto_75")
                   & (bt["defect injected"] == "the wrong root of eq. (7)")][0])
WRONGROOT_FACTOR = float(bt.broken[_wr]/bt.baseline[_wr])
WRONGROOT_RANK_ABS = int((shift > shift[_wr]).sum()) + 1
WRONGROOT_RANK_REL = int((rel > rel[_wr]).sum()) + 1
BIGGEST_ABS_ROW = bt.loc[shift.idxmax()]
print(f"  THE ROW THE REUSE SECTION QUOTES -- the wrong root of eq. (7) against"
      f" fig7_max_dX_upto_75 --")
print(f"  moves it {bt.baseline[_wr]:.6f} -> {bt.broken[_wr]:.4f}: a factor"
      f" {WRONGROOT_FACTOR:.1f} = {np.log10(WRONGROOT_FACTOR):.3f} orders.")
print(f"  Ranked both ways over the {len(bt)} rows: {WRONGROOT_RANK_ABS} by ABSOLUTE shift"
      f" (the largest is {shift.max():.2e},")
print(f"  {BIGGEST_ABS_ROW['metric']}), but only {WRONGROOT_RANK_REL} by RATIO --"
      f" {WRONGROOT_RANK_REL - 1} rows move their metric by a")
print("  larger relative factor. It is therefore quoted as the measured factor and BOTH ranks,")
print("  never as a bare superlative that only one of the two readings supports.\n")

report_agreement("B3.3", metrics)'''))

# ------------------------------------------------------------- prose audit
cells.append(code(r'''# Every number quoted in this page's MARKDOWN is re-derived here and compared
# against the live computation; any mismatch raises and fails the notebook.
def close(a, b, rtol=5e-3, atol=1e-15):
    if not abs(a - b) <= atol + rtol*abs(b):
        raise AssertionError(f"prose drift: typed {a!r} vs computed {b!r}")
    return True

AUDIT = [
    ("fifteen symbolic identities (title, 2.1)",         15, SYM_ZERO_COUNT, 0),
    ("m*(0) = 1.0000000 (title, 6.5)",                   1.0, M_PSI0, 1e-8),
    ("m*(2) = 0.4905224 (title, 6.5)",                   0.4905224, M_PSI2, 1e-6),
    ("m*(1) = 0.6638850 (title, 6.5)",                   0.6638850, M_PSI1, 1e-6),
    ("0.42 % below two-thirds (title, 6.5)",             0.0042, abs(M_PSI1_REL_23), 5e-2),
    ("X_M -> 1-exp(-1/2) = 0.3934693 (title, 6.1)",      0.3934693, XM_SUP, 1e-6),
    ("2.9 % at sigma = 100 (title, 7.4)",                0.029, SIGMA100_MAX_REL, 2e-2),
    ("7.2 % spread at sigma = 0.25 (title, 7.4)",        0.072, SIGMA025_SPREAD, 2e-2),
    ("1.1e-16 B3.1 identity (title, 7.5)",               1.1e-16, B31_IDENTITY, 2e-1),
    ("0.0218 up to X = 0.75 (title, 6.3)",               0.0218, DAB_UPTO75, 5e-3),
    ("2.2 conversion points (title)",                    0.022, DAB_UPTO75, 5e-2),
    ("X = 0.7423 for the 75 % claim (title, 6.3)",       0.7423, X_AT_2PT, 1e-3),
    ("X = 0.7599 under (eps_o, S_o) (6.3, README, meta)", 0.7599, X_AT_2PT_ALT, 1e-3),
    ("X = 0.7687 under (S_o, L_o) (6.3, README, meta)",   0.7687, X_AT_2PT_SL, 1e-3),
    ("0.0177 up to X = 0.75 under (eps_o, S_o) (meta)",  0.017652, DAB_UPTO75_ALT, 1e-3),
    ("0.0165 up to X = 0.75 under (S_o, L_o) (meta)",    0.016507, DAB_UPTO75_SL, 1e-3),
    ("wrong root of (7): factor 1672 on the 6.3 comparison (Reuse, README, meta)",
                                                         1672.0, WRONGROOT_FACTOR, 1e-3),
    ("... 0.0218 -> 36.4 (Reuse)",                        36.4, float(bt.broken[_wr]), 1e-3),
    ("... 5th by absolute shift (Reuse, README, meta)",   5, WRONGROOT_RANK_ABS, 0),
    ("... 22nd by relative shift (Reuse, README, meta)",  22, WRONGROOT_RANK_REL, 0),
    ("53 break rows (7.7, README, meta)",                53, len(bt), 0),
    ("eight metrics carry two rows (7.7)",               8,
                                                         int((bt["metric"].value_counts() > 1).sum()), 0),
    ("-21.7 % surface at X = 0.75, Fig-4 set (title cell only)",   -0.2166, SURF_DEV_F4_75, 5e-3),
    ("-21.0 % surface at X = 0.75, Fig-7 A/B set (title cell only)", -0.2100, SURF_DEV_AB_75, 5e-3),
    ("eps_o root 0.08529 (title, 6.6)",                  0.08529, EPS_O_ROOT, 1e-3),
    ("residual -4.2e-3 at the printed 0.1 (title, 6.6)", 4.2e-3, abs(STAT35_AT_PRINTED), 2e-2),
    ("0.003 % surface cost (title, 6.6)",                2.76e-5, S_DEFICIT_AT_PRINTED, 5e-2),
    ("S_max = 2720 cm2/cm3 (title, 6.2)",                2720.0, SETS[2]["Smax"], 1e-3),
    ("factor 4.60 beyond it (title, 6.2)",               4.60, CURVE_C_INFEAS, 5e-3),
    ("curve D L = 7.747e7 cm/cm3 (title, 6.2)",          7.747e7, L_D, 1e-3),
    ("psi = 4.965 for the Fig. 4 set (6.2)",             4.965326, PSI_F4, 1e-5),
    ("eq. (7) roots -5.181, 1.0517, 4.1293 (Reuse)",     4.1293, G_SEL_F7, 1e-4),
    ("... the middle root 1.0517 (Reuse)",               1.0517,
     float(np.sort(np.roots([(4/27)*P["fig7_eps_o"], 0, -1, 1]).real)[1]), 1e-4),
    ("... and -5.181 (Reuse)",                          -5.181,
     float(np.sort(np.roots([(4/27)*P["fig7_eps_o"], 0, -1, 1]).real)[0]), 1e-3),
    ("observed order 2.000 (7.1, What pymrm adds)",      2.0, MARCH_ORDER, 3e-3),
    ("more than 8 % at psi = 2 (What pymrm adds)",       0.0799, NRMS_BEST_PSI2, 2e-2),
    ("20 % at psi = 5 (What pymrm adds)",                0.1995, I_nrms(5.0, m_star(5.0)), 2e-2),
    ("psi at which m* reaches 0 = 16.80 (6.5)",          16.8046, PSI_M0, 1e-3),
    ("m = 1 misfits 17.2 % at psi = 1 (6.6)",            0.1717, NRMS_M1_PSI1, 2e-2),
    ("... a factor 5.3 worse (6.6)",                     5.34, M1_PENALTY, 5e-2),
    ("Fig-8 intercept 27.04 x 10^12 (6.6)",              27.040, INTERCEPT/1e12, 1e-4),
    ("X_M 0.2989 (VC) (6.6)",                            0.298877, XM_VC, 1e-5),
    ("X_M 0.3475 (CS) (6.6)",                            0.347541, XM_CS, 1e-5),
    ("phi < 0.5557 for eta > 0.98 (6.6)",                0.555707, PHI_098, 1e-5),
    ("overlap gap 0.0301 (6.4)",                         0.030118, OVERLAP_MAX_GAP, 1e-3),
    ("overlap crossing u = 0.6508 (6.4)",                0.650806, U_CROSS, 1e-4),
    ("the gamma-coefficient row moves the B3.1 identity to 1.0 (7.5)",
                                                         1.0, BRK_B31_GAMMA_RATIO, 1e-9),
    ("the exponent-on-x row moves it to 0.25 (7.5)",     0.25, BRK_B31_IDENTITY, 1e-9),
    (f"{len(scal)} printed scalars (data)",              len(scal), len(scal), 0),
]
for label, typed, computed, rt in AUDIT:
    if rt == 0:
        assert typed == computed, f"prose drift [{label}]: {typed} != {computed}"
    else:
        try:
            close(typed, computed, rtol=rt)
        except AssertionError as e:
            raise AssertionError(f"[{label}] {e}") from None
print(f"prose audit: {len(AUDIT)} numbers re-derived and matched. "
      "Any drift raises and fails the notebook.")'''))

# ------------------------------------------------------------- references
cells.append(md(r"""## References

Bhatia, S. K. and Perlmutter, D. D. (1980). A random pore model for
fluid–solid reactions: I. Isothermal, kinetic control. *AIChE Journal*
**26**(3), 379–386.
[doi:10.1002/aic.690260308](https://doi.org/10.1002/aic.690260308) — **the
paper, and the only document read for content.** Identity confirmed from its
own first page on a native-resolution render: the title, the by-line "S. K.
BHATIA and D. D. PERLMUTTER / Department of Chemical Engineering / University
of Pennsylvania / Philadelphia, Pennsylvania 19104", the copyright line
"0001-1541-80-3242-0378-\$00.95. © The American Institute of Chemical
Engineers, 1980", and the running feet "AIChE Journal (Vol. 26, No. 3) /
May, 1980 / Page 379" through "Page 385"; the Wiley Subject metadata gives
"AIChE Journal 1980.26:379-386". Manuscript received April 16, 1979; revision
received July 20, and accepted August 20, 1979. `pdfimages -list` reports all
eight pages as CCITT-G4 bilevel at **300 ppi native**, so every numeric was
read from a crop of a 300 ppi render at digit scale; the text layer was used
only as a search index. (The eighth PDF page carries the start of the *next*
article — Riazi and Whitson, "Application of Corresponding States Principles
for Prediction of Self-Diffusion Coefficients in Liquids" — below this paper's
reference list; nothing was read from it.)

**Part II is not on disk**, was not consulted, and nothing here depends on it.

**Cited by the paper, not consulted, and no number here derives from them.**
Petersen, E. E. (1957), "Reaction of Porous Solids", *AIChE J.* **3**, 442 —
the model this page's §6.2–6.4 compare against; its eqs. (3)–(5) are read
**only** as Bhatia and Perlmutter reprint them on p. 380, and every statement
about "Petersen's model" here is a statement about that reprint. Szekely, J.,
Evans, J. W. and Sohn, H. Y. (1976), *Gas–Solid Reactions* (Academic Press) —
the source of eqs. (6) and (7), same treatment. Avrami, M. (1940), *J. Chem.
Phys.* **8**, 212 — the overlap result eq. (19). Hashimoto, K., Miura, K.,
Yoshikawa, F. and Imai, I. (1979), *Ind. Eng. Chem. Process Des. Develop.*
**18**, 73 — the char data of Figure 8, **not on disk**, which is why the
empirical comparison is out of scope. Walker, P. L. and Raats, E. (1956);
Kawahata, M. and Walker, P. L. (1962); Dutta, S., Wen, C. Y. and Belt, R. J.
(1977); Dutta, S. and Wen, C. Y. (1977) — the surface-area and reactivity
observations quoted in §6.1. Lacey, Bowen and Basden (1965); Ishida and Wen
(1971); Calvelo and Smith (1971); Calvelo and Cunningham (1970); Ramachandran
and Smith (1977); Hashimoto and Silveston (1973); Szekely and Evans (1970,
1971); Szekely, Lin and Sohn (1973); Takamura, Yoshida and Kunii (1974);
Mendoza, Cunningham and Ronco (1970); Barner and Mantell (1968); Yamadaya
et al. (1970); Young (1966); Tompkins (1976); Bhatia and Perlmutter (1979) —
the modelling literature the paper surveys.

**Park, J. Y. and Levenspiel, O.**, "The Crackling Core Model for the Reaction
of Solid Particles", *Chem. Eng. Sci.* **30**, 1207 — dated **(1975)** in the
paper's LITERATURE CITED and **(1976)** in its p. 380 text. Reported in §4;
which year is correct cannot be settled from this document.

**Gallery pages read for this one**, neither of which supplies a dataset:
`B3.1` (Yagi & Kunii 1955, shrinking core) — its published eq. (6) is used in
§7.5 and its published `eq6_to_reaction_limit` metric is reproduced and
reconciled there; and `B3.2` (Szekely & Evans 1970, grain model) — read, not
loaded, and nothing of its numerics is used here."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py", "mimetype": "text/x-python", "name": "python",
        "nbconvert_exporter": "python", "pygments_lexer": "ipython3",
        "version": "3.13.5"},
}
out = Path(__file__).with_name("index.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
