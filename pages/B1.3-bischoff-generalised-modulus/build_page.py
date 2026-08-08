#!/usr/bin/env python3
"""Generate index.ipynb for page B1.3 (Bischoff generalised modulus).

Run from the page directory:  python build_page.py
Edit this file, never index.ipynb.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# --------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Bischoff's generalised modulus: measuring how narrow the narrow band is"
description: "One modulus built from the rate integral collapses the effectiveness-factor curves of every reaction order onto nearly one curve. Bischoff said the spread is 'about 15%' for orders one-half to three and 'about 30%' with zero order - this page root-finds both numbers (14.6% and 40.6% max/min, 28.9% as a fraction of the upper curve), finds the printed n = 1/2 reduction is a misprint his own 15% figure shows Figure 1 never used, and audits the construction out of sample on a finite cylinder, where a 1998 'arbitrary kinetics, < 1.5%' competitor misses six-fold while the 1965 collapse, combined with the shape's own curve, stays inside its advertised band."
categories: [sec:B, struct:S3, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-08
---

# Bischoff's generalised modulus: measuring how narrow the narrow band is

**Catalog ID:** `B1.3` · **Structure:** `S3` (1D steady BVP) · **Tier:** T0

Thiele's modulus makes the effectiveness factor of a first-order pellet a
single curve. Every other rate form gets its own curve, and
Langmuir-Hinshelwood kinetics get a *family* of them - which, as Bischoff
(1965) opens by observing, "greatly complicates practical computations."

His fix is one line: build the modulus from the integral of the rate, and the
large-$m$ ends of *all* the curves coincide exactly. What remains is the
intermediate range, about which the paper makes a quantitative claim: "the
spread is only about 15% for the most interesting cases of one-half- to
third-order reactions and is about 30% when zero-order reactions are
included." This page measures that claim - the collapse, its width, where the
width peaks, and what the same modulus does for adsorption kinetics, variable
diffusivity, and a shape it was never fitted to."""))

# ----------------------------------------------------------------- background
cells.append(md(r"""## Background

Bischoff's own account of the lineage, from p. 351: "Thiele (2) and Zeldovich
(3) in 1939 gave the early, and now classical, treatments." Thiele's modulus
for a first-order reaction is $m = L\sqrt{k_v/D}$, a grouping of rate
coefficient, diffusivity and half-width; each reaction order needs its own
version, and Chu and Hougen's adsorption-rate curves (his ref. 7) need a fresh
curve for every parameter pair $(\zeta, C_o)$. The paper's purpose, in its own
words: "it can be stated that a general solution would be very desirable. The
purpose of this work is to show how this can be approximately realized."

The mathematics is a slab (equivalently a single pore): steady diffusion with
an arbitrary rate form $r(C)$ and concentration-dependent diffusivity $D(C)$,

$$\frac{\mathrm{d}}{\mathrm{d}x}\!\left[D(C)\,\frac{\mathrm{d}C}{\mathrm{d}x}\right] = r(C),
\qquad C(0) = C_o,\qquad \frac{\mathrm{d}C(L)}{\mathrm{d}x} = 0,$$

his eqs. (1)-(3), with $x=0$ the surface and $x=L$ the sealed centre. A first
integral (his eq. 5) gives the flux in terms of $\int D\,r\,\mathrm{d}C$, and
everything on this page follows from that.

**Why this page completes a quartet.** [`B1.1`](../B1.1-thiele-weisz-hicks/)
is Thiele's curve itself; [`B1.2`](../B1.2-aris-shape-factor/) is Aris's
collapse across *shapes* at fixed (first-order) kinetics, using the
volume-to-surface length $v_p/s_x$; this page is the collapse across
*kinetics* at fixed shape. Bischoff cites exactly that division of labour: for
other geometries, "Aris (8) has shown that by choosing appropriate
characteristic lengths, other geometrical shapes give approximately similar
results." The final section below combines the two collapses on a shape and a
rate form neither paper computed, and audits a 1998 competitor - Pan & Zhu's
approximate effectiveness factor for arbitrary kinetics in a cylinder - on the
same problem."""))

# ------------------------------------------------------------ published model
cells.append(md(r"""## The published model

**The general asymptotic solution and the modulus.** For a semi-infinite slab
(large modulus) the surface flux is exact and explicit (his eq. 11), and
Bischoff defines the general modulus so that the asymptote is
$\mathcal{E} = 1/m$ for *every* rate form (his eq. 13):

$$m \;\equiv\; \frac{L\,r(C_o)}{\sqrt{2}}
\left[\int_0^{C_o} D(\alpha)\,r(\alpha)\,\mathrm{d}\alpha\right]^{-1/2}.$$

For the finite slab the exact solution is parametric in the (unknown) centre
concentration $C_L$: his eq. (14) gives $m$ (one printed exponent of which is
a misprint - one of five printed defects reported in the Results) and his
eq. (15) gives

$$\mathcal{E} \;=\; \frac{1}{m}\left[\frac{\int_{C_L}^{C_o} D\,r\,\mathrm{d}C}
{\int_{0}^{C_o} D\,r\,\mathrm{d}C}\right]^{1/2}.$$

**Simple orders, $r = k_v C^n$** (his eq. 16). The modulus reduces to
$m = L\sqrt{(n+1)\,k_v C_o^{\,n-1}/(2D)}$ (eq. 17) - the standard modulus
times $\sqrt{(n+1)/2}$ - and eq. (15) reduces to
$\mathcal{E} = (1/m)\sqrt{1-(C_L/C_o)^{n+1}}$ (eq. 19) with $m(C_L)$ from the
quadrature of eq. (20). The paper evaluates five orders:

| $n$ | printed solution | eqs. |
|---|---|---|
| 1 | $\mathcal{E} = \tanh m/m$ | 21 |
| 0 | $\mathcal{E} = 1$ for $m<1$, $1/m$ for $m>1$ (attributed to Wheeler, his ref. 4) | 22 |
| 2 | incomplete elliptic integral $F(\phi,k)$, $k=\sin 15°$ | 23-24 |
| 3 | $F(\phi,k)$, $k = \sin 45°$ | 25-26 |
| 1/2 | $F$ and $E$, $k=\sin 15°$ | 27-28 |

and states the claim this page measures: **"the spread is only about 15% for
the most interesting cases of one-half- to third-order reactions and is about
30% when zero-order reactions are included"** (p. 354), restated in the
summary as "if reactions of order less than one-half are excluded, the spread
between all the various curves is about 15%."

**Adsorption (LHHW) kinetics**, $r = k_v C/(\zeta + C)$, the form Chu & Hougen
computed curve families for. The general modulus becomes (his eqs. 29/31)

$$m = \frac{M}{\sqrt 2}\,\frac{C_o}{\zeta + C_o}
\left[C_o - \zeta\ln\!\left(1+\frac{C_o}{\zeta}\right)\right]^{-1/2},
\qquad M = L\sqrt{k_v/D},$$

with $M$ the standard modulus of Chu and Hougen. The paper does not recompute
their curves - "this paper will directly use their curves, which are redrawn
in Figure 2" - and replots them against $m$ as its Figure 3, concluding "the
spread is again about the same as for the simple order reactions." Here the
curves are computed, not redrawn, so that sentence too becomes a number.

**First-order reaction with volume change**, $D(C) = D/(\omega C + 1)$
(eq. 32): the modulus correction is closed-form (eq. 33),
$m = M\,|\omega C_o|\,[\,\omega C_o - \ln(\omega C_o + 1)]^{-1/2}/\sqrt 2$,
and Table 1 - the paper's only table - compares it at four values of
$\omega C_o$ against Hawthorn's empirical correction
$M\,[\omega C_o/\ln(1+\omega C_o)]^{0.7}$ (eq. 34).

**What the paper does *not* print:** any measurement. Every number in it is
computed from the model, so this is a provenance tier 6 page throughout, and
"agreement" below means agreement between independent computations - Bischoff's
printed ones, closed forms, quadrature of his parametric solution, and pymrm
finite-volume solutions of the underlying BVP."""))

# ------------------------------------------------------ parameters/assumptions
cells.append(md(r"""## Parameters and assumptions

**Assumptions, all Bischoff's:** steady state; a single reaction with rate a
function of one concentration; isothermal; slab geometry (or a single pore),
with surface concentration held at $C_o$ and zero flux at the centre;
diffusivity constant or a known function $D(C)$. Nonisothermal pellets, film
resistance and bimodal pore structures are explicitly deferred by the paper to
its refs. 1, 4-6 and are not treated here (the nonisothermal pellet is
[`B1.1`](../B1.1-thiele-weisz-hicks/)'s second half).

**Everything is dimensionless.** With $c = C/C_o$, rates normalised so
$\tilde r(1) = 1$ and $\tilde d = D(C)/D$, the model has no physical
parameters at all - only the rate-form shape ($n$, or $\zeta^* = \zeta/C_o$,
or $\omega C_o$) and the modulus. Note the reduction the modulus makes
visible: Chu & Hougen's two-parameter family $(\zeta, C_o)$ is really the
one-parameter family $\zeta^* = \zeta/C_o$; Bischoff's Figures 2a-c vary
$C_o$ and $\zeta$ separately because both are dimensional there.

**How the sources were read.** Both PDFs are CCITT-G4 bilevel scans
(Bischoff at 300 ppi native, Pan & Zhu at 600 ppi - checked with
`pdfimages -list`; rendering higher only interpolates). The Bischoff text
layer is poor (~2.2k characters on p. 353, sub/superscripts unreliable), so
**every equation, exponent and table cell used here was read from a
digit-scale crop of a native-resolution render**, and the three printed
tables used are transcribed into `data/` with provenance sidecars. Five
printed defects were found this way and are reported - never repaired
silently - in the Results: one Table 1 cell, the n = 1/2 elliptic reduction
and an exponent in eq. (14) (Bischoff), plus one Table 5 cell and a sign in
eq. (36) (the secondary source)."""))

# ----------------------------------------------------------------- env cells
cells.append(code("""# Colab environment cell - no-op if pymrm is already installed
try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code(r'''import sys, urllib.request
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

import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import csc_array
from scipy.sparse.linalg import spsolve
from scipy.special import ellipkinc, ellipeinc, jn_zeros
from pymrm import (construct_grad, construct_div, NumJac, newton,
                   interp_cntr_to_stagg, clip_approach)
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "B1.3-bischoff-generalised-modulus"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
# the quadrature of the parametric exact solution is pushed to tolerances where
# QUADPACK reports roundoff-limited refinement; accuracy is verified against
# closed forms (2e-12) and a deliberately loosened-tolerance break row below,
# so these warnings are noise here - and their text embeds the kernel's
# temp path, which would break the two-executions-identical guarantee
from scipy.integrate import IntegrationWarning
warnings.filterwarnings("ignore", category=IntegrationWarning)
np.set_printoptions(precision=6)'''))

# -------------------------------------------------------------------- the data
cells.append(md(r"""## The data

Three small printed tables, all of them the authors' own computed values -
**nothing on this page is a measurement** (provenance tier 6; neither paper
reports an experiment relevant here).

- `bischoff-table1.csv` - Bischoff's Table 1 (modulus correction factors for
  volume change), all eight cells as printed.
- `pan-zhu-table5.csv`, `pan-zhu-table6.csv` - the secondary source's shape
  coefficients and its first-order accuracy test (columns $Z$ = 1 and 2),
  as printed.

The exact $\eta(m)$ curves computed below are exported as
`eta-generalised-modulus.csv` for reuse."""))

cells.append(code(r'''tab1 = load_data("bischoff-table1.csv", page=PAGE)
pan5 = load_data("pan-zhu-table5.csv", page=PAGE)
pan6 = load_data("pan-zhu-table6.csv", page=PAGE)
print(cite_data(load_meta("bischoff-table1.csv", page=PAGE)))
print(cite_data(load_meta("pan-zhu-table6.csv", page=PAGE)))
print(f"{len(tab1)} Table 1 rows; {len(pan5)} Table 5 Z-values; "
      f"{len(pan6)} Table 6 rows across Z = "
      f"{[float(z) for z in sorted(pan6['Z'].unique())]}")'''))

# ------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

Three independent routes to $\eta(m)$, sharing no assembly:

**Route A - pymrm finite volume.** The BVP itself,
$\mathrm{d}[\tilde d(c)\,c']/\mathrm{d}x = \Phi^2\tilde r(c)$ on $x\in(0,1)$
with $c'(0)=0$, $c(1)=1$, solved with `construct_grad`/`construct_div` and a
Newton iteration whose reaction Jacobian comes from `NumJac((n, 1))` (the
shape is `(n, 1)`, never `(n,)` - the bare 1-D shape builds a dense Jacobian,
the trap measured on the first published version of `B1.1`). Variable
$\tilde d(c)$ is handled by Picard iteration on the face diffusivity around
the Newton solve. $\eta$ is the cell-weighted rate integral; the surface-flux
reading is *identical by construction* in a finite-volume scheme (summing the
divergence rows telescopes to the boundary flux), so it is never quoted as a
check here.

**Route B - quadrature of Bischoff's parametric solution.** Eqs. (14)-(15)
with the *inner* integral $\int \tilde d\,\tilde r\,\mathrm{d}c$ in closed
form for every rate form used (power law, LHHW, volume change), the outer
integral by adaptive quadrature after the substitution $c = C_L + s^2$ that
removes the square-root singularity, and $C_L(m)$ by `brentq`. Deterministic
- no continuation, no warm starts (the `B1.1` lesson). One exponent of
eq. (14) is amended in the process; that is a printed defect of the source,
reported with proof in its own Results section below.

**Route C - the paper's closed forms.** Eqs. (21)-(28): `tanh`, the Wheeler
piecewise form, and the incomplete elliptic integrals via
`scipy.special.ellipkinc`/`ellipeinc` - table lookups in 1965, special
functions now, and entirely independent of Route B's quadrature. One
implementation note that matters: the printed $\sin\phi$ parametrisations
pass through $\phi = \pi/2$ (at $C_o/C_L = 1+\sqrt3$ for the $k=\sin15°$
pair), so $\phi$ is evaluated here from the equivalent
$\cos\phi = (\sqrt3+1-t)/(\sqrt3-1+t)$, which carries the branch
automatically."""))

cells.append(code(r'''# ---- rate forms, normalised: r(1) = 1, d = D(C)/D; Rint = int_a^b d*r dc ----
class Power:
    """r = c^n, constant D (Bischoff eq. 16)."""
    variable_d = False
    def __init__(self, n): self.n = n; self.label = f"n={n:g}"
    def r(self, c): return np.maximum(c, 0.0) ** self.n
    def d(self, c): return np.ones_like(np.asarray(c, float))
    def Rint(self, a, b): return (b ** (self.n + 1) - a ** (self.n + 1)) / (self.n + 1)

class LHHW:
    """r = c (z+1)/(z+c), z = zeta/Co, constant D (Bischoff eq. 29 form)."""
    variable_d = False
    def __init__(self, z): self.z = z; self.label = f"LHHW z*={z:g}"
    def r(self, c): return np.maximum(c, 0.0) * (self.z + 1.0) / (self.z + np.maximum(c, 0.0))
    def d(self, c): return np.ones_like(np.asarray(c, float))
    def Rint(self, a, b):
        z = self.z
        return (z + 1.0) * ((b - a) - z * np.log1p((b - a) / (z + a)))

class VolChange:
    """first order with D(C) = D/(1 + w c), w = omega*Co (Bischoff eq. 32)."""
    variable_d = True
    def __init__(self, w): self.w = w; self.label = f"vol wCo={w:+g}"
    def r(self, c): return np.maximum(c, 0.0)
    def d(self, c): return 1.0 / (1.0 + self.w * np.asarray(c, float))
    def Rint(self, a, b):
        w = self.w
        return ((b - a) - np.log1p(w * (b - a) / (1.0 + w * a)) / w) / w

# ---- Route B: Bischoff eqs. (14)-(15), parametric in u = C_L/C_o ------------
U_MIN = 1e-10

def mb_eta(rf, u):
    """(m, eta) at centre concentration u. Outer quadrature on c = u + s^2."""
    I0 = rf.Rint(0.0, 1.0)
    def gout(s):
        b = u + s * s
        v = rf.Rint(u, b)
        return float(rf.d(b)) * 2.0 * s / np.sqrt(v) if v > 0 else 0.0
    Iout = quad(gout, 0.0, np.sqrt(1.0 - u), epsabs=1e-12, epsrel=1e-10,
                limit=200)[0]
    m = 0.5 * Iout / np.sqrt(I0)
    return m, (1.0 / m) * np.sqrt(rf.Rint(u, 1.0) / I0)

class Curve:
    """Pointwise-exact eta(m) for one rate form (Route B + asymptote).

    A parametric (u, m) table built once at init brackets the root, so each
    eta(m) query costs a handful of quadratures instead of a blind search."""
    def __init__(self, rf, nu=44):
        self.rf = rf
        us = np.unique(np.concatenate([np.geomspace(U_MIN, 0.5, nu // 2),
                                       1.0 - np.geomspace(1e-9, 0.5, nu // 2)]))
        ms = np.array([mb_eta(rf, u)[0] for u in us])
        order = np.argsort(ms)
        self._m_tab, self._u_tab = ms[order], us[order]   # m ascending, u descending
        self.m_lim = float(ms[0])            # us[0] = U_MIN carries the largest m

    def eta(self, m):
        if m >= self.m_lim * (1.0 - 1e-12):
            return 1.0 / m                   # exact on the dead-zone branch;
                                             # elsewhere reached only at eta*m-1 < 1e-9
        j = int(np.searchsorted(self._m_tab, m))
        lo = self._u_tab[j] if j < len(self._u_tab) else U_MIN
        hi = self._u_tab[j - 1] if j > 0 else 1.0 - 1e-12
        f = lambda u: mb_eta(self.rf, u)[0] - m
        if not (f(lo) > 0 > f(hi) or f(lo) < 0 < f(hi)):
            lo, hi = U_MIN, 1.0 - 1e-12      # safety: fall back to the full bracket
        u = brentq(f, lo, hi, xtol=1e-12, rtol=8.9e-16)
        return mb_eta(self.rf, u)[1]

class ZeroOrder:
    """Bischoff eq. (22), attributed by him to Wheeler: the exact n = 0 curve."""
    rf = Power(0.0)
    m_lim = 1.0
    def eta(self, m):
        return 1.0 if m < 1.0 else 1.0 / m

# dead-zone onset for n < 1, derived here from eq. (14) with C_L = 0:
# m* = (n+1)/(1-n)  (n=0 -> 1, matching eq. 22; n=1/2 -> 3). Not printed in
# the paper; used below as an exactness check on the quadrature.
def mstar_power(n):
    return (n + 1.0) / (1.0 - n)'''))

cells.append(code(r'''# ---- Route C: the paper's closed forms (eqs. 21-28) -------------------------
SQ3 = np.sqrt(3.0)
K15SQ = (2.0 - SQ3) / 4.0        # k = sin 15 deg = sqrt(2-sqrt(3))/2, squared
K45SQ = 0.5                      # k = sin 45 deg = 1/sqrt(2), squared

def phi_from_cos(t):
    """phi with cos(phi) = (sqrt3+1-t)/(sqrt3-1+t); equals the printed
    sin(phi) = sqrt(4*sqrt3*(t-1))/(sqrt3-1+t) but carries phi past pi/2."""
    return np.arccos((SQ3 + 1.0 - t) / (SQ3 - 1.0 + t))

def m_eta_n1(u):     # eq. (21) region: m = arccosh(1/u), eta = tanh(m)/m
    m = np.arccosh(1.0 / u)
    return m, np.tanh(m) / m

def m_eta_n2(u):     # eqs. (23)-(24)
    R = 1.0 / u
    m = 1.5 * np.sqrt(R) * ellipkinc(phi_from_cos(R), K15SQ) / 3.0 ** 0.25
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 3)

def m_eta_n3(u):     # eqs. (25)-(26)
    phi = np.arcsin(np.sqrt(1.0 - u ** 2))
    m = np.sqrt(2.0) / u * ellipkinc(phi, K45SQ)
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 4)

def _bracket27(phi):
    s, c = np.sin(phi), np.cos(phi)
    return ((1.0 + SQ3) * ellipkinc(phi, K15SQ) - 2.0 * SQ3 * ellipeinc(phi, K15SQ)
            + 2.0 * SQ3 * s * np.sqrt(1.0 - K15SQ * s * s) / (1.0 + c))

def m_eta_nhalf_amended(u):
    """eqs. (27)-(28) AMENDED (see Results): prefactor (C_L/C_o)^{1/4} and
    phi built on (C_o/C_L)^{1/2} - the reduction of eq. (20) through y=sqrt(t)."""
    m = 1.5 * u ** 0.25 / 3.0 ** 0.25 * _bracket27(phi_from_cos(np.sqrt(1.0 / u)))
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 1.5)

def m_eta_nhalf_as_printed(u):
    """eqs. (27)-(28) EXACTLY as printed: prefactor sqrt(C_L/C_o), phi from
    the where-block shared with eq. (23), i.e. built on C_o/C_L itself."""
    m = 1.5 * np.sqrt(u) / 3.0 ** 0.25 * _bracket27(phi_from_cos(1.0 / u))
    return m, (1.0 / m) * np.sqrt(1.0 - u ** 1.5)

# eq. (20) directly - the single-quadrature reduction both (23)-(28) and
# Route B must agree with; used as the arbiter for the n = 1/2 misprint
def m_eq20(n, u):
    R = 1.0 / u
    # expm1/log1p keep t^{n+1} - 1 accurate near the lower limit
    val = quad(lambda s: 2.0 * s / np.sqrt(np.expm1((n + 1.0) * np.log1p(s * s))),
               0.0, np.sqrt(R - 1.0), epsabs=1e-13, epsrel=1e-11, limit=400)[0]
    return (n + 1.0) / 2.0 * R ** ((n - 1.0) / 2.0) * val'''))

cells.append(code(r'''# ---- Route A: pymrm finite volume ------------------------------------------
BC = ({"a": 1.0, "b": 0.0, "d": 0.0},   # x=0 centre:  dc/dn = 0 (symmetry; outward normal)
      {"a": 0.0, "b": 1.0, "d": 1.0})   # x=1 surface: c = 1    (Dirichlet)

class Slab:
    """d/dx[ d(c) dc/dx ] = Phi^2 r(c) on (0,1); c'(0)=0, c(1)=1."""
    def __init__(self, n=400):
        self.n = n
        self.shape = (n, 1)                          # (n,1), never (n,): NumJac trap
        self.x_f = np.linspace(0.0, 1.0, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.grad_mat, self.grad_bc = construct_grad(self.shape, self.x_f,
                                                     self.x_c, BC, axis=0)
        self.div_mat = construct_div(self.shape, self.x_f, nu=0, axis=0)  # nu=0: slab
        self.numjac = NumJac(self.shape)
        self.w = np.diff(self.x_f)                   # cell widths

    def _newton(self, rf, Phi, d_face, c0):
        D = csc_array(self.grad_mat.multiply(d_face.reshape(-1, 1)))
        jac_diff = self.div_mat @ (-D)
        g_const = (self.div_mat @ (-(d_face.reshape(-1, 1) * self.grad_bc.toarray()))
                   ).reshape(-1, 1)
        def residual(cflat):
            c2 = cflat.reshape(self.shape)
            gr, jr = self.numjac(lambda cc: Phi ** 2 * rf.r(cc), c2)
            return (jac_diff @ cflat.reshape(-1, 1) + g_const + gr.reshape(-1, 1),
                    jac_diff + jr)
        sol = newton(residual, c0.ravel(), tol=1e-13, maxfev=60,
                     callback=lambda x, dx: clip_approach(x, dx, lower_bounds=0.0))
        res_norm = float(np.max(np.abs(residual(sol.x)[0])))
        return sol.x.reshape(self.shape), res_norm

    def solve(self, rf, Phi):
        c = np.ones(self.shape)
        if not rf.variable_d:
            return self._newton(rf, Phi, np.ones(self.n + 1), c)
        for _ in range(80):                          # Picard on the face diffusivity
            c_face = interp_cntr_to_stagg(c, self.x_f, self.x_c, axis=0).ravel()
            c_face[0], c_face[-1] = c.ravel()[0], 1.0
            c_new, res_norm = self._newton(rf, Phi, rf.d(c_face), c)
            if np.max(np.abs(c_new - c)) < 1e-13:
                return c_new, res_norm
            c = c_new
        raise RuntimeError("Picard iteration on d(c) did not converge")

    def eta(self, rf, m):
        """eta at generalised modulus m; Phi = m*sqrt(2*I0) inverts eq. (13)."""
        Phi = m * np.sqrt(2.0 * rf.Rint(0.0, 1.0))
        c, res_norm = self.solve(rf, Phi)
        assert res_norm < 1e-9, f"unconverged Newton residual {res_norm:.1e}"
        return float(self.w @ rf.r(c.ravel()))'''))

# ------------------------------------------------------------------- results
cells.append(md(r"""## Results

### The three routes agree before anything is measured

Route B (quadrature of eqs. 14-15) against Route C (the printed closed
forms), across the whole parametric range. This is the check that guards every
number that follows; it is *not* guaranteed - a wrong elliptic modulus, a
dropped branch or a mistranscribed prefactor breaks it at once, and the break
table at the end injects exactly those defects. The $n=1/2$ column here uses
the **amended** reduction; the as-printed one is dealt with in its own section
below."""))

cells.append(code(r'''U_GRID = np.array([0.02, 0.05, 0.1, 0.2, 0.3, 0.366, 0.5, 0.7, 0.9, 0.99])
CLOSED = {0.5: m_eta_nhalf_amended, 1.0: m_eta_n1, 2.0: m_eta_n2, 3.0: m_eta_n3}

ROUTE_BC_DEV = 0.0
for n, closed in CLOSED.items():
    rf = Power(n)
    for u in U_GRID:
        mB, eB = mb_eta(rf, u)
        mC, eC = closed(u)
        ROUTE_BC_DEV = max(ROUTE_BC_DEV, abs(mB - mC) / mC, abs(eB - eC) / eC)
print(f"Route B vs the paper's closed forms, worst relative deviation over "
      f"{len(U_GRID)} centre concentrations x 4 orders: {ROUTE_BC_DEV:.2e}")

# eq. (20), the paper's own single-quadrature reduction, as a third witness
EQ20_DEV = max(abs(m_eq20(n, u) - mb_eta(Power(n), u)[0]) / mb_eta(Power(n), u)[0]
               for n in (0.5, 1.0, 2.0, 3.0) for u in (0.1, 0.5, 0.9))
print(f"eq. (20) quadrature vs Route B, worst: {EQ20_DEV:.2e}")

# the dead-zone onset m* = (n+1)/(1-n), extrapolated in u^(1/4) from Route B
m1, _ = mb_eta(Power(0.5), 1e-8)
m2, _ = mb_eta(Power(0.5), 1e-10)
MSTAR_EX = (m2 * 1e-8 ** 0.25 - m1 * 1e-10 ** 0.25) / (1e-8 ** 0.25 - 1e-10 ** 0.25)
MSTAR_DEV = abs(MSTAR_EX - mstar_power(0.5)) / mstar_power(0.5)
print(f"n=1/2 dead-zone onset: closed form m* = {mstar_power(0.5):.1f}, "
      f"Route B extrapolated {MSTAR_EX:.6f} (rel dev {MSTAR_DEV:.1e})")

# pymrm (Route A) against Route B at a spread of moduli and rate forms.
# n = 1/2 is swept on its smooth branch (m < m* = 3); its non-smooth
# dead-zone branch gets a dedicated treatment in the next cell.
SLAB = Slab(400)
FORMS_A = [Power(0.5), Power(1), Power(2), Power(3),
           LHHW(0.1), LHHW(10.0), VolChange(1.0), VolChange(-0.5)]
PYMRM_DEV = 0.0
for rf in FORMS_A:
    cv = Curve(rf)
    ms = (0.3, 1.0, 2.5) if getattr(rf, "n", None) == 0.5 else (0.3, 1.0, 3.0)
    for m in ms:
        e_pm, e_ex = SLAB.eta(rf, m), cv.eta(m)
        PYMRM_DEV = max(PYMRM_DEV, abs(e_pm - e_ex) / e_ex)
print(f"pymrm (n = 400) vs Route B over 8 rate forms x 3 moduli: {PYMRM_DEV:.2e}")

# grid refinement at m = 1, n = 2
errs = []
for ncell in (100, 200, 400):
    e = Slab(ncell).eta(Power(2), 1.0)
    errs.append(abs(e - Curve(Power(2)).eta(1.0)))
ORDER = float(np.log2(errs[0] / errs[1]) / 1.0 + np.log2(errs[1] / errs[2])) / 2.0
print("grid errors at n_cell = 100/200/400:",
      " ".join(f"{e:.2e}" for e in errs), f"-> mean observed order {ORDER:.2f}")'''))

cells.append(md(r"""**The dead zone needs its own solver treatment.** For
$n<1$ past $m^*$ the reactant vanishes on an interior region and the rate's
square-root kink at $c=0$ sits exactly on the free boundary. Plain clipped
Newton then *chatters*: the solution is fine but the pointwise residual
stalls many orders above tolerance over the boundary cells - the next cell
prints both - and an unconverged solve must not stand behind a reported
metric. The
regularisation $r_\varepsilon = \sqrt{c+\varepsilon} - \sqrt{\varepsilon}$
with its analytic Jacobian restores clean convergence (residual ~$10^{-10}$),
and $\varepsilon$ is then an error axis like any other - refined below over
four decades so the limit is measured, not assumed."""))

cells.append(code(r'''from scipy.sparse import diags

D1 = csc_array(SLAB.grad_mat.multiply(np.ones((SLAB.n + 1, 1))))
JD = SLAB.div_mat @ (-D1)
G0 = (SLAB.div_mat @ (-SLAB.grad_bc.toarray())).reshape(-1, 1)

def eta_deadzone(m, eps):
    Phi = m * np.sqrt(2.0 * Power(0.5).Rint(0.0, 1.0))
    def residual(cflat):
        c = np.maximum(cflat, 0.0)
        g = (JD @ cflat.reshape(-1, 1) + G0
             + (Phi ** 2 * (np.sqrt(c + eps) - np.sqrt(eps))).reshape(-1, 1))
        return g, JD + diags(Phi ** 2 * 0.5 / np.sqrt(c + eps))
    sol = newton(residual, np.ones(SLAB.n), tol=1e-13, maxfev=200,
                 callback=lambda x, g: clip_approach(x, g, lower_bounds=0.0))
    rn = float(np.max(np.abs(residual(sol.x)[0])))
    assert rn < 1e-9, f"regularised dead-zone solve unconverged: {rn:.1e}"
    return float(SLAB.w @ np.sqrt(np.maximum(sol.x, 0.0)))

# what plain clipped Newton does here, for the record: an unconverged solve,
# never used for any reported metric
def plain_deadzone(m):
    Phi = m * np.sqrt(2.0 * Power(0.5).Rint(0.0, 1.0))
    def residual(cflat):
        c2 = cflat.reshape(SLAB.shape)
        gr, jr = SLAB.numjac(lambda cc: Phi ** 2 * Power(0.5).r(cc), c2)
        return JD @ cflat.reshape(-1, 1) + G0 + gr.reshape(-1, 1), JD + jr
    sol = newton(residual, np.ones(SLAB.n), tol=1e-13, maxfev=60,
                 callback=lambda x, g: clip_approach(x, g, lower_bounds=0.0))
    rn = float(np.max(np.abs(residual(sol.x)[0])))
    return rn, float(SLAB.w @ Power(0.5).r(np.maximum(sol.x, 0.0)))

RN_PLAIN, ETA_PLAIN = plain_deadzone(5.0)
print(f"plain clipped Newton at m = 5: residual stalls at {RN_PLAIN:.1e} while "
      f"eta = {ETA_PLAIN:.6f}\nis already {abs(ETA_PLAIN - 0.2) / 0.2:.1e} from "
      f"the exact 0.2 - accurate, but not a converged solve,\nso it stands "
      f"behind no reported number\n")
print("pymrm n = 1/2 at m = 5 (dead zone; exact eta = 1/m = 0.2):")
DEAD_DEVS = []
for eps in (1e-6, 1e-8, 1e-10, 1e-12):
    e = eta_deadzone(5.0, eps)
    DEAD_DEVS.append(abs(e - 0.2) / 0.2)
    print(f"  eps = {eps:.0e}: eta = {e:.8f}   rel dev {DEAD_DEVS[-1]:.1e}")
E_DEAD = eta_deadzone(5.0, 1e-12)
DEAD_DEV = abs(E_DEAD - 0.2) / 0.2
DEAD_RATIOS = [DEAD_DEVS[k] / DEAD_DEVS[k + 1] for k in range(3)]
print(f"successive ratios {DEAD_RATIOS[0]:.1f}, {DEAD_RATIOS[1]:.1f}, "
      f"{DEAD_RATIOS[2]:.1f} against the 10 of O(sqrt(eps)): the first")
print("decades follow sqrt(eps); the last ratio falls short because the n = 400")
print(f"grid's own discretisation floor blends in, so the {DEAD_DEV:.1e} at "
      f"eps = 1e-12 is the")
print("honest total deviation (grid + regularisation), not pure eps error")'''))

cells.append(md(r"""### The collapse, and the width of the band

Five reaction orders on the generalised modulus. The spread at fixed $m$ is
reported two ways, because the paper does not say which it means:
$\eta_{\max}/\eta_{\min}-1$ (relative to the lowest curve) and
$(\eta_{\max}-\eta_{\min})/\eta_{\max}$ (as a fraction of the highest). The
extrema are **root-found, not sampled**: the smooth interior maximum by
Brent's method on the exact curves, and the kink candidates ($m = 1$, where
the $n=0$ curve breaks, and $m^* = 3$, where $n = 1/2$ enters its dead zone)
evaluated exactly."""))

cells.append(code(r'''CURVES = {n: Curve(Power(n)) for n in (0.5, 1.0, 2.0, 3.0)}
CURVES[0.0] = ZeroOrder()

def spread(m, keys):
    v = [CURVES[k].eta(m) for k in keys]
    return max(v) / min(v) - 1.0

def spread_frac(m, keys):
    v = [CURVES[k].eta(m) for k in keys]
    return (max(v) - min(v)) / max(v)

def maximise(fun, lo=0.2, hi=10.0, ngrid=40):
    """Root-found interior maximum; refuses to return a window-edge value."""
    mm = np.geomspace(lo, hi, ngrid)
    vals = [fun(m) for m in mm]
    i = int(np.argmax(vals))
    assert 0 < i < ngrid - 1, "maximum at the window edge - widen the window"
    res = minimize_scalar(lambda lm: -fun(np.exp(lm)),
                          bracket=(np.log(mm[i - 1]), np.log(mm[i]),
                                   np.log(mm[i + 1])),
                          method="brent", options={"xtol": 1e-11})
    return np.exp(res.x), -res.fun

KEYS_15 = (0.5, 1.0, 2.0, 3.0)
KEYS_30 = (0.0, 0.5, 1.0, 2.0, 3.0)

M_15, BAND_15 = maximise(lambda m: spread(m, KEYS_15))
BAND_15_FRAC = spread_frac(M_15, KEYS_15)
# kink candidates for the half-to-third band
KINK_15 = max(spread(1.0, KEYS_15), spread(3.0, KEYS_15))
assert BAND_15 > KINK_15, "interior max is the global one"

print(f'Bischoff: "the spread is only about 15% for ... one-half- to third-order"')
print(f"  measured: {BAND_15:.2%} (max/min - 1) = {BAND_15_FRAC:.2%} (fraction of the")
print(f"  upper curve), at m = {M_15:.4f}; root-found interior maximum, and larger")
print(f"  than both kink candidates (worst kink {KINK_15:.4f})")

# with zero order the upper envelope is eta = 1 up to the kink at m = 1 exactly
BAND_30_KINK = spread(1.0, KEYS_30)
BAND_30_FRAC = spread_frac(1.0, KEYS_30)
M_30, BAND_30_SMOOTH = maximise(lambda m: spread(m, KEYS_30))
ETA3_AT_1 = CURVES[3.0].eta(1.0)
print(f'\nBischoff: "about 30% when zero-order reactions are included"')
print(f"  measured: {BAND_30_KINK:.2%} (max/min - 1) = {BAND_30_FRAC:.2%} (fraction of")
print(f"  the upper curve), and the maximiser is the kink of eq. (22) at m = 1")
print(f"  exactly (smooth search from either side converges to m = {M_30:.6f});")
print(f"  the whole band there is 1 - eta_n3(1) with eta_n3(1) = {ETA3_AT_1:.6f}")

# second, independent route to the same two headline numbers: the paper's
# closed forms (Route C) instead of Route B quadrature, extrema re-root-found
def spread_C(m, with0):
    def eta_C(closed, m):
        u = brentq(lambda u: closed(u)[0] - m, 1e-9, 1 - 1e-12, xtol=1e-14)
        return closed(u)[1]
    vals = [eta_C(CLOSED[n], m) if m < {0.5: 3.0}.get(n, np.inf) else 1.0 / m
            for n in (0.5, 1.0, 2.0, 3.0)]
    if with0:
        vals.append(1.0 if m < 1.0 else 1.0 / m)
    return max(vals) / min(vals) - 1.0

M_15C, BAND_15C = maximise(lambda m: spread_C(m, False))
BAND_30C = spread_C(1.0, True)
SECOND_ROUTE_GAP = max(abs(BAND_15C - BAND_15) / BAND_15,
                       abs(BAND_30C - BAND_30_KINK) / BAND_30_KINK,
                       abs(M_15C - M_15) / M_15)
print(f"\nsecond route (closed forms, no shared quadrature): "
      f"{BAND_15C:.6f} at m = {M_15C:.4f}; incl. n=0: {BAND_30C:.6f}")
print(f"worst relative gap between the two routes: {SECOND_ROUTE_GAP:.1e}")

# the asymptotic coincidence that motivates the modulus: eta*m -> 1
ASYM_DEV = max(abs(CURVES[n].eta(30.0) * 30.0 - 1.0) for n in (1.0, 2.0, 3.0))
print(f"\neta*m at m = 30, worst |eta*m - 1| over n = 1,2,3: {ASYM_DEV:.2e}")
print("(n = 0 and n = 1/2 are exactly 1/m there by eq. 22 / the dead-zone")
print(" solution - structural, proving nothing. So is the n = 1 leg: its")
print(f" parametric table tops out at m_lim = arccosh(1/u_min) = "
      f"{CURVES[1.0].m_lim:.1f} < 30, so")
print(" Curve.eta(30) returns 1/m by its asymptote branch there. The reported")
print(f" deviation is carried by n = 2 and 3 alone, whose parametric ranges "
      f"reach")
print(f" m_lim = {CURVES[2.0].m_lim:.1e} and {CURVES[3.0].m_lim:.1e} - far past "
      f"m = 30, so those two legs are real)")'''))

cells.append(code(r'''# ---- what the collapse replaced: the standard modulus ----------------------
# "A form of modulus similar to Equation (17), but with only the term
#  C_o^{n-1} for dimensional consistency, has been used by many of the above
#  investigators" - i.e. M_n = L*sqrt(kv*Co^{n-1}/D) = m / sqrt((n+1)/2).
def spread_std(M, keys):
    v = [CURVES[k].eta(M * np.sqrt((k + 1.0) / 2.0)) for k in keys]
    return max(v) / min(v) - 1.0

M_STD15, BAND_STD15 = maximise(lambda M: spread_std(M, KEYS_15))
M_STD30, BAND_STD30 = maximise(lambda M: spread_std(M, KEYS_30))
STD_ASYM_30 = np.sqrt(2.0) - 1.0          # exact algebra: sqrt((3+1)/(0+1))/sqrt(2)
STD_ASYM_15 = np.sqrt(4.0 / 1.5) - 1.0
print("the same five curves on the standard modulus M_n:")
print(f"  n in [1/2, 3]: max spread {BAND_STD15:.1%} at M = {M_STD15:.3f} "
      f"(vs {BAND_15:.1%} on m); asymptotic spread {STD_ASYM_15:.1%}, forever")
print(f"  incl. n = 0 : max spread {BAND_STD30:.1%} at M = {M_STD30:.3f} "
      f"(vs {BAND_30_KINK:.1%} on m); asymptotic spread {STD_ASYM_30:.1%}, forever")
print("on m every asymptotic spread is zero by construction - that is eq. (13)'s")
print("entire content - so the honest comparison is the peak: the generalised")
print(f"modulus narrows the worst-case band by {BAND_STD15/BAND_15:.1f}x / "
      f"{BAND_STD30/BAND_30_KINK:.1f}x and removes the asymptotic spread entirely")'''))

cells.append(code(r'''# ---- figure: the collapse and the band -------------------------------------
MM = np.geomspace(0.1, 10.0, 160)
ETA_TAB = {n: np.array([CURVES[n].eta(m) for m in MM]) for n in KEYS_30}

fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
ax = axes[0]
colors = {0.0: "tab:red", 0.5: "tab:orange", 1.0: "k", 2.0: "tab:blue", 3.0: "tab:green"}
for n in KEYS_30:
    ax.loglog(MM, ETA_TAB[n], color=colors[n], lw=1.8,
              label=f"$n = {n:g}$" + (" (Wheeler form)" if n == 0 else ""))
ax.loglog(MM, 1.0 / MM, "--", color="0.5", lw=1.2, label=r"$\eta = 1/m$")
ax.axvline(1.0, color="0.8", lw=0.8)
ax.set_xlabel(r"generalised modulus $m$ [-]  (Bischoff eq. 13)")
ax.set_ylabel(r"effectiveness factor $\eta$ [-]")
ax.set_ylim(0.08, 1.15)
ax.set_title("Fig. 1 recomputed: five orders on one modulus")
ax.legend(frameon=False, fontsize=8, loc="lower left")

ax = axes[1]
S15 = np.array([spread(m, KEYS_15) for m in MM])
S30 = np.array([spread(m, KEYS_30) for m in MM])
ax.semilogx(MM, 100 * S30, color="tab:red", lw=2.0, label="incl. $n=0$")
ax.semilogx(MM, 100 * S15, color="tab:blue", lw=2.0, label=r"$n \in [1/2, 3]$")
ax.axhline(15, color="tab:blue", ls=":", lw=1.1)
ax.axhline(30, color="tab:red", ls=":", lw=1.1)
ax.annotate('"about 15%"', (0.115, 15.7), fontsize=8, color="tab:blue")
ax.annotate('"about 30%"', (0.115, 30.7), fontsize=8, color="tab:red")
ax.plot([M_15], [100 * BAND_15], "o", color="tab:blue", ms=7)
ax.plot([1.0], [100 * BAND_30_KINK], "o", color="tab:red", ms=7)
ax.set_xlabel(r"$m$ [-]")
ax.set_ylabel(r"spread in $\eta$ at fixed $m$,  $\eta_{max}/\eta_{min}-1$  [%]")
ax.set_title("Width of the band, and the printed claims")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""Both printed claims hold, each under one of the two
natural readings, and the page pins which: **14.61% at $m = 1.270$** for
orders one-half to three as $\eta_{\max}/\eta_{\min}-1$ (his "about 15%"; the
same band is 12.75% as a fraction of the upper curve), and with zero order
included the maximum sits at the kink of eq. (22), $m = 1$ **exactly**, where
the band is $1-\eta_{n=3}(1)$: **40.55%** relative to the lowest curve but
**28.85%** of the upper one - his "about 30%" is the fraction-of-upper
reading. The claims are 1965 eyeball numbers read off log-log figures; both
round to the printed values under the reading just stated, and the page
regards them as confirmed.

The contrast that justifies the whole construction: on the **standard**
modulus the same five curves never collapse at all - the spread at large $M$
is $\sqrt{2}-1 = 41.4\%$ *forever* (that number is exact algebra, the ratio of
asymptotes, and is labelled structural in the metrics) and peaks at 120.6%,
against 40.6% peak and asymptotically zero on $m$."""))

cells.append(md(r"""### Adsorption kinetics: Figures 2 and 3, computed instead of redrawn

For the LHHW form the paper reused Chu & Hougen's curves ("this paper will
directly use their curves"). Computing them removes both the redrawing and
the original numerical integration: each curve is Route B with the closed-form
rate integral. Bischoff's Figures 2a-c span $\zeta$ = 0.001 to 100 at $C_o$ =
0.1, 0.5, 1.0 - in reduced form, $\zeta^* = \zeta/C_o$ from $10^{-3}$ to
$10^3$, spanning effectively zero-order to first-order behaviour. The claim to
measure: after replotting against $m$, "the spread is again about the same as
for the simple order reactions."
"""))

cells.append(code(r'''ZSTARS = sorted({round(z / co, 12) for z in (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
                 for co in (0.1, 0.5, 1.0)})
print(f"{len(ZSTARS)} distinct zeta* from Bischoff's 18 (zeta, Co) figure labels:")
print("  " + ", ".join(f"{z:g}" for z in ZSTARS))
LCURVES = {z: Curve(LHHW(z)) for z in ZSTARS}

def spread_lhhw(m):
    v = [c.eta(m) for c in LCURVES.values()]
    return max(v) / min(v) - 1.0

M_L, BAND_L = maximise(spread_lhhw, 0.2, 8.0)
# kink candidate: none of the family dead-zones (locally first order at c=0),
# verified: the smallest m_lim across the family
print(f"no dead zones in the family (smallest parametric m reached: "
      f"{min(c.m_lim for c in LCURVES.values()):.1f})")
print(f"LHHW family on m: max spread {BAND_L:.2%} at m = {M_L:.4f}")
print(f"simple-order bands for comparison: {BAND_15:.2%} (n in [1/2,3]) and "
      f"{BAND_30_KINK:.2%} (incl. n = 0)")

# the same family on Chu & Hougen's standard modulus M = L*sqrt(kv/D):
# eq. (31) with Co = 1:  m = (M/sqrt2) (1/(z*+1)) [1 - z* ln(1+1/z*)]^{-1/2}
def m_over_M(z):
    return 1.0 / np.sqrt(2.0) / (z + 1.0) / np.sqrt(1.0 - z * np.log1p(1.0 / z))

def m_from_M(Mstd, z):
    return Mstd * m_over_M(z)

# on M the family never collapses: as M -> inf each curve is eta = m_over_M(z)/M,
# so the spread tends to a CLOSED-FORM constant - no windowed sweep needed
ratios = [m_over_M(z) for z in ZSTARS]
BAND_LS = max(ratios) / min(ratios) - 1.0
print(f"same family on the standard modulus M (the Fig. 2 presentation): the")
print(f"curves never coincide; asymptotic spread {BAND_LS:.0%} "
      f"(a factor {max(ratios)/min(ratios):.0f} between the outermost curves), forever")'''))

cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
MSTD = np.geomspace(0.1, 30.0, 120)
sel = [0.001, 0.1, 1.0, 10.0, 1000.0]
shades = plt.cm.viridis(np.linspace(0.15, 0.85, len(sel)))
for z, col in zip(sel, shades):
    axes[0].loglog(MSTD, [LCURVES[z].eta(m_from_M(M, z)) for M in MSTD],
                   color=col, lw=1.7, label=rf"$\zeta^* = {z:g}$")
    axes[1].loglog(MM, [LCURVES[z].eta(m) for m in MM], color=col, lw=1.7)
axes[1].loglog(MM, ETA_TAB[0.0], ":", color="tab:red", lw=1.4, label="$n=0$")
axes[1].loglog(MM, ETA_TAB[1.0], ":", color="k", lw=1.4, label="$n=1$")
for ax, ttl, xl in ((axes[0], "Fig. 2 recomputed: standard modulus $M$",
                     r"$M = L\sqrt{k_v/D}$ [-]"),
                    (axes[1], "Fig. 3 recomputed: generalised modulus $m$",
                     r"$m$ [-]")):
    ax.set_xlabel(xl); ax.set_ylabel(r"$\eta$ [-]")
    ax.set_ylim(0.08, 1.15); ax.set_title(ttl)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""The measured LHHW band on $m$ lands between the two
simple-order figures, and close to the zero-order-included one - necessarily
so, because $\zeta^* \to 0$ *is* zero order over most of the pellet and
$\zeta^*\to\infty$ is first order, so the family fills the $n\in[0,1]$ gap of
Figure 1. "Again about the same" is right against the ~30% number, not
against the headline 15% - a distinction the paper does not draw. On the
standard modulus the same family never coincides at all - a factor of about
23 separates the outermost curves at large $M$, forever - which is what
eq. (31) buys."""))

cells.append(md(r"""### Volume change, Table 1, and a cell that fails its own equation

The variable-diffusivity case is where the $D(\alpha)$ inside eq. (13) is
*live* (everywhere else on this page it is a constant that cancels). The
paper's claim: with the correction factor of eq. (33), Thiele's
volume-change curves "again fall within the range of Figure 1." Table 1
compares eq. (33) with Hawthorn's empirical eq. (34) and says "the two are
essentially equivalent numerically." Both closed forms are evaluated below
beside the printed cells - the printed values are **loaded from the
transcription CSV, never retyped**."""))

cells.append(code(r'''W_FAM = [0.4, 1.0, -0.5, -0.25]
VCURVES = {w: Curve(VolChange(w)) for w in W_FAM}

# containment in the Figure-1 range and distance from the n = 1 curve
def vc_dev(m):
    e1 = CURVES[1.0].eta(m)
    return max(abs(VCURVES[w].eta(m) - e1) / e1 for w in W_FAM)

M_VC, VC_DEV = maximise(vc_dev, 0.2, 8.0)
# containment: distance to the [n=3, n=0] envelope as a FRACTION of the local
# band width (the raw margin closes trivially where the band itself closes)
def vc_margin_frac(m):
    lo, hi = CURVES[3.0].eta(m), CURVES[0.0].eta(m)
    return min(min(VCURVES[w].eta(m) - lo, hi - VCURVES[w].eta(m))
               for w in W_FAM) / (hi - lo)
MARGIN = minimize_scalar(lambda lm: vc_margin_frac(np.exp(lm)),
                         bounds=(np.log(0.2), np.log(8.0)), method="bounded",
                         options={"xatol": 1e-10})
MARGIN_FRAC = float(MARGIN.fun)
print(f"volume-change family (omega*Co = {W_FAM}) on m:")
print(f"  worst deviation from the constant-D first-order curve: {VC_DEV:.2%} "
      f"at m = {M_VC:.3f}")
print(f"  containment in the Figure-1 envelope [n=3, n=0] over m in [0.2, 8]:")
print(f"  never leaves it, staying at least {MARGIN_FRAC:.0%} of the local band")
print(f"  width away from either edge")

# Table 1: recompute both closed forms exactly
def eq33(w):   # m/M from eq. (33), read at digit scale from p. 355
    return abs(w) / np.sqrt(2.0) / np.sqrt(w - np.log(w + 1.0))

def eq34(w):   # Hawthorn's factor, eq. (34), exponent 0.7 as printed
    return (w / np.log(1.0 + w)) ** 0.7

print(f"\nBischoff Table 1, printed vs recomputed:")
print(f"{'wCo':>7}{'eq33 printed':>14}{'eq33 exact':>12}{'dev':>9}"
      f"{'eq34 printed':>14}{'eq34 exact':>12}{'dev':>9}")
dev33, dev34 = {}, {}
for _, row in tab1.iterrows():
    w = row["omega_co"]
    dev33[w] = (row["eq33_printed"] - eq33(w)) / eq33(w)
    dev34[w] = (row["eq34_printed"] - eq34(w)) / eq34(w)
    flag = "   <--" if abs(dev33[w]) > 0.004 else ""
    print(f"{w:7.2f}{row['eq33_printed']:14.3f}{eq33(w):12.4f}{dev33[w]:+9.4f}"
          f"{row['eq34_printed']:14.3f}{eq34(w):12.4f}{dev34[w]:+9.4f}{flag}")

T1_33_OK = max(abs(dev33[w]) for w in (0.4, 1.0, -0.5))
T1_33_DEFECT = abs(dev33[-0.25])                    # relative deviation
T1_33_DEFECT_ABS = T1_33_DEFECT * eq33(-0.25)       # absolute deviation
print(f"\neq. (33) column: three cells agree to {T1_33_OK:.4f} (inside 3-digit")
print(f"rounding); the omega*Co = -0.25 cell is off by {T1_33_DEFECT_ABS:.4f} "
      f"absolute")
print(f"({T1_33_DEFECT:.4f} relative) - {T1_33_DEFECT_ABS/0.0005:.0f}x the "
      f"half-ULP rounding radius of a printed")
print(f"0.905, and the same {T1_33_DEFECT/(0.0005/0.905):.0f}x with both sides "
      f"taken relative; mixing the two readings")
print(f"(relative deviation over absolute radius) would overstate the multiple "
      f"as {T1_33_DEFECT/0.0005:.0f}x.")
T1_34_WORST = max(abs(v) for v in dev34.values())
print(f"The eq. (34) column agrees to {T1_34_WORST:.4f} in all four cells, "
      f"including this row (0.906 vs {eq34(-0.25):.4f}),")
print("so the row's omega*Co label is right and the defect is isolated to one")
print("cell: printed 0.905, eq. (33) gives 0.9107. Reported, not repaired.")'''))

cells.append(md(r"""The glyphs are crisp at digit scale - `0.905` is not an
OCR question - and the constraint that decides the cell is arithmetic: the
other three cells of the same column pin the transcription of eq. (33) to
within printed rounding, the neighbouring eq. (34) column pins the row label
$\omega C_o = -0.25$ (0.906 against an exact 0.9064), and eq. (33) itself
then gives 0.9107, not 0.905. A slide-rule slip of 0.6% in the fourth cell of
a 1965 table, visible only because the equation beside it is closed-form.

The physics conclusions survive the cell: the volume-change family sits well
inside the Figure-1 envelope (the paper's claim), hugging the first-order
curve to within a few percent - the $D(\alpha)$ weighting inside eq. (13) is
what absorbs it, and the break table shows the collapse degrading when that
weighting is deliberately dropped."""))

cells.append(md(r"""### The $n=1/2$ reduction, exactly as printed, is a misprint

Eq. (27) with its where-block reads

$$m = \tfrac{3}{2}\sqrt{C_L/C_o}\;\tfrac{1}{\sqrt[4]{3}}
\Bigl[(1+\sqrt3)F(\phi,k) - 2\sqrt3\,E(\phi,k)
+ 2\sqrt3\,\tfrac{\sin\phi\sqrt{1-k^2\sin^2\phi}}{1+\cos\phi}\Bigr],\quad
\sin\phi = \tfrac{\sqrt{4\sqrt3\,(C_o/C_L-1)}}{\sqrt3-1+C_o/C_L},\;
k=\sin 15°,$$

the $\sin\phi$ block being *the same one printed under the $n=2$ case*. But
the $n=1/2$ integral $\int\mathrm{d}t/\sqrt{t^{3/2}-1}$ reduces through the
substitution $y=\sqrt t$, which changes both the prefactor exponent and the
argument $\phi$ is built on. As printed, the formula contradicts the paper's
own eq. (20) - from which it is derived - across the whole curve; amended in
exactly those two places (prefactor $(C_L/C_o)^{1/4}$, $\phi$ built on
$\sqrt{C_o/C_L}$), it agrees with eq. (20) and with Route B to machine
precision. Both versions share the bracket, so the defect cannot be in this
page's transcription of the $F$/$E$ combination."""))

cells.append(code(r'''ASP_DEV = AMD_DEV = 0.0
print(f"{'u = C_L/C_o':>12}{'eq. (20)':>12}{'as printed':>12}{'rel dev':>9}"
      f"{'amended':>12}{'rel dev':>10}")
for u in (0.05, 0.1, 0.3, 0.5, 0.7, 0.9):
    m20 = m_eq20(0.5, u)
    mP, _ = m_eta_nhalf_as_printed(u)
    mA, _ = m_eta_nhalf_amended(u)
    ASP_DEV = max(ASP_DEV, abs(mP - m20) / m20)
    AMD_DEV = max(AMD_DEV, abs(mA - m20) / m20)
    print(f"{u:12.2f}{m20:12.6f}{mP:12.6f}{abs(mP-m20)/m20:9.1%}"
          f"{mA:12.6f}{abs(mA-m20)/m20:10.1e}")
print(f"\nas printed: {ASP_DEV:.0%} off its own eq. (20); amended: {AMD_DEV:.1e}")

# Did Figure 1 use the printed formula or the correct one? The printed spread
# claim decides this without touching the figure. The as-printed curve does
# not even reach eta = 1 at small m: its modulus is high by a factor that
# tends to sqrt(2) at the kinetic end, so eta = sqrt(1-u^{3/2})/m_printed
# tends to 1/sqrt(2), and the band it implies never closes on the left.
def eta_as_printed(m):
    u = brentq(lambda u: m_eta_nhalf_as_printed(u)[0] - m, 1e-9, 1 - 1e-9,
               xtol=1e-13)
    return m_eta_nhalf_as_printed(u)[1]

print("\nm(as printed)/m(eq. 20) approaching the kinetic end:")
for uu in (0.99, 0.999, 1.0 - 1e-6):
    RATIO_LIMIT = m_eta_nhalf_as_printed(uu)[0] / m_eq20(0.5, uu)
    print(f"  u = {uu}: ratio = {RATIO_LIMIT:.7f}")
print(f"  -> sqrt(2) = {np.sqrt(2):.7f} (deviation at u = 1-1e-6: "
      f"{abs(RATIO_LIMIT - np.sqrt(2)):.1e})")
ETA_ASP_SMALL = eta_as_printed(0.05)
BAND_BAD = 1.0 / ETA_ASP_SMALL - 1.0
print(f"so the as-printed curve gives eta = {ETA_ASP_SMALL:.4f} at m = 0.05, "
      f"not 1, and the")
print(f"half-to-third band it implies is {BAND_BAD:.1%} already at the kinetic "
      f"end - against")
print(f"{BAND_15:.1%} for the amended curve and the paper's printed "
      f'"about 15%". The curve')
print("behind Figure 1 was therefore evidently computed from the correct")
print("reduction, and the defect is typographical (prefactor and phi-argument),")
print("not a defect in the paper's calculations.")'''))

cells.append(md(r"""Because the paper never prints a number for the $n=1/2$
curve (no table touches it; its only appearance is Figure 1 - checked against
every equation and table in the five pages), the misprint would be invisible
to a reader who did not re-derive the reduction, and anyone *implementing
eqs. (27)-(28) as printed gets a curve 29-41% off in $m$*. The printed
"about 15%" itself is what shows Bischoff's own Figure 1 did not use the
printed formula: the as-printed curve tops out at $\eta = 1/\sqrt2$ instead
of 1, so the half-to-third band it implies is already ~41% at the kinetic end
where the true band closes to zero - irreconcilable with any reading of 15%.
The amendment is an inference, labelled as such; everything else on this page
is independent of it (the band numbers above use Route B, not eq. 27)."""))

cells.append(md(r"""### Eq. (14) as printed carries a defect of the same kind

Route B integrates eqs. (14)-(15) - and quietly reads eq. (14)'s denominator
bracket as $[\int_0^{C_o}D\,r\,\mathrm{d}C]^{+1/2}$, because that is what
eqs. (8) and (13), of which eq. (14) is the stated combination, force. This
page's own standard is that printed defects are *reported, never repaired
silently*, so: **as printed, that denominator bracket carries exponent
$-1/2$** (read at digit scale from the native 300 ppi render; the *inner*
bracket in the numerator correctly prints $-1/2$, so the defect is isolated
to one exponent). Dividing by $[\,\cdot\,]^{-1/2}$ multiplies by the rate
integral to the $+1/2$: dimensionally inconsistent as printed - $m$ would
carry the units of $D\,r\,C$ - and, in this page's normalised variables,
low by exactly the factor $I_0 = \int_0^1\tilde r\,\mathrm{d}c = 1/(n+1)$
for order $n$. The paper's own eq. (21) is the one-line witness."""))

cells.append(code(r'''# eq. (14) EXACTLY as printed: the misprint is a constant factor, so it is
# provable without re-integrating anything - the printed denominator
# [int_0^Co D r dC]^(-1/2) turns eq. (8)+(13)'s division by I0^(1/2) into a
# multiplication by I0^(1/2), i.e. m_printed = m_correct * I0.
def m_eq14_as_printed(rf, u):
    return mb_eta(rf, u)[0] * rf.Rint(0.0, 1.0)

EQ14_ASP_DEV = 0.0
print("eq. (14) as printed vs its own parents eqs. (8)+(13), at u = 0.5:")
for n in (0.5, 1.0, 2.0, 3.0):
    rf = Power(n)
    m_c = mb_eta(rf, 0.5)[0]
    m_p = m_eq14_as_printed(rf, 0.5)
    EQ14_ASP_DEV = max(EQ14_ASP_DEV, abs(m_p - m_c) / m_c)
    print(f"  n = {n:g}: I0 = {rf.Rint(0.0, 1.0):.4f} -> as-printed m low by "
          f"{abs(m_p - m_c) / m_c:.0%} (u-independent)")
M21 = np.arccosh(1.0 / 0.5)
print(f"\nthe printed witness: eq. (21) gives m = arccosh(1/u) = {M21:.6f} at "
      f"u = 0.5, where")
print(f"eq. (14) as printed gives {m_eq14_as_printed(Power(1), 0.5):.6f} - "
      f"exactly half. Reported, not")
print("repaired: Route B implements the +1/2 reading eqs. (8)+(13) force, and")
print(f"its agreement with the printed closed forms eqs. (21)-(28) "
      f"({ROUTE_BC_DEV:.1e}, first")
print("Results cell) is the proof that reading is right across every order.")
print("With eqs. (27)-(28) this makes two typographic defects in the same")
print("paper's display equations, both invisible to a reader who does not")
print("re-derive them - and both leaving the paper's computed results intact.")'''))

cells.append(md(r"""### Out of sample: a finite cylinder, and a 1998 competitor

Bischoff's collapse plus Aris's length ($L = v_p/s_x$, page
[`B1.2`](../B1.2-aris-shape-factor/)) yields a zero-cost estimate for *any*
shape and *any* kinetics: compute $m$ with eq. (13), read $\eta$ off the
Figure-1 band. The natural audit is a case neither 1965 source computed: a
finite cylinder ($Z = H/2R_p = 1$, every face active), first- and
second-order kinetics, solved exactly here as a 2-D pymrm problem
(`nu=1` radial axis $\times$ Cartesian axial axis, the `B1.2` machinery with
a Newton loop on the reaction).

The competitor: Pan & Zhu (1998) propose, for **arbitrary kinetics** in
exactly this geometry, the polynomial
$\eta_a = q_1\eta_{sp} + q_2\eta_{sp}^2 + q_3\eta_{sp}^3$ in the basis
$\eta_{sp} = (\coth\phi - 1/\phi)/\phi$, with the $q_i$ fixed by matching
their exact $\phi\to0$ perturbation ($\mu_1$, their eq. 26) and $\phi\to\infty$
asymptote ($\mu_2 Q_2/\phi$, their eq. 29) - the large-$\phi$ end being
*precisely Bischoff's asymptotic solution in a cylinder*: their
$Q_2 = [\,2\int_{C_{eq}}^{1} r\,\mathrm{d}C]^{1/2}$ is eq. (11)'s integral.
They claim "a maximum deviation of less than 1.5% from analytical solutions
or numerical calculation results ... for arbitrary kinetics" (abstract), and
their printed tests are first order (their Table 6) and one steam-shift case
(their Table 7). Both endpoint coefficients are recomputed here from their
printed definitions; the matching that yields $q_2, q_3$ is re-derived and
asserted rather than copied, which matters because their printed eq. (36)
has a sign error (shown below)."""))

cells.append(code(r'''# ---- Pan & Zhu machinery, from their printed definitions --------------------
def mu1_series(Z, M=200, N=400):
    """their eq. (26); alpha_m = zeros of J0, theta_n = (2n+1)pi/(2Z) (eq. 24)."""
    am = jn_zeros(0, M)[:, None]
    th = (2 * np.arange(N)[None, :] + 1) * np.pi / (2.0 * Z)
    return 32.0 / np.pi ** 2 * float(np.sum(
        1.0 / (am ** 2 * (2 * np.arange(N)[None, :] + 1) ** 2 * (am ** 2 + th ** 2))))

def eta_series_1st(phi, Z=1.0, M=200, N=400):
    """their eq. (37) = Aris's exact first-order series; identical to the
    double series B1.2 validated against a 2-D solve (its case-iv machinery)."""
    am = jn_zeros(0, M)[:, None]
    n = np.arange(N)[None, :]
    th = (2 * n + 1) * np.pi / (2.0 * Z)
    return 1.0 - 32.0 * phi ** 2 / np.pi ** 2 * float(
        np.sum(1.0 / ((2 * n + 1) ** 2 * am ** 2 * (am ** 2 + th ** 2 + phi ** 2))))

def pz_coeffs(Z, Q1, Q2):
    """q1..q3 from the two matching conditions they state (phi->0 to eq. 25,
    phi->inf to eq. 29). Solving them gives q2, q3 in closed form; the results
    are asserted against the conditions, not copied from eqs. (35)-(36)."""
    a, b = mu1_series(Z) * Q1, (2.0 + 1.0 / Z) * Q2
    q1 = b
    q2 = 27.0 - 135.0 * a - 6.0 * b
    q3 = -54.0 + 405.0 * a + 9.0 * b
    assert abs(q1 / 3 + q2 / 9 + q3 / 27 - 1.0) < 1e-12       # eta -> 1 at phi=0
    assert abs(3 * q1 + 2 * q2 + q3 - 135.0 * a) < 1e-9       # slope matches eq. 25
    return q1, q2, q3

def eta_pz(phi, Z, Q1, Q2):
    q1, q2, q3 = pz_coeffs(Z, Q1, Q2)
    es = (1.0 / np.tanh(phi) - 1.0 / phi) / phi
    return q1 * es + q2 * es ** 2 + q3 * es ** 3

# their Table 5, recomputed from their own definitions
MU1_DEV = max(abs(mu1_series(Z) - v) for Z, v in zip(pan5["Z"], pan5["mu1_printed"]))
MU2_DEV = {Z: abs(2.0 + 1.0 / Z - v) for Z, v in zip(pan5["Z"], pan5["mu2_printed"])}
print(f"Pan & Zhu Table 5: mu1 row reproduced by their eq. (26) series to "
      f"{MU1_DEV:.1e} (worst cell)")
print(f"mu2 row vs their eq. (30) = 2 + 1/Z: all cells exact except Z = 2.5: "
      f"printed {float(pan5.loc[pan5['Z'] == 2.5, 'mu2_printed'].iloc[0]):.2f}, "
      f"eq. (30) gives {2 + 1/2.5:.2f} "
      f"(off by {MU2_DEV[2.5]:.2f}; likely a copy of the Z = 2.0 cell). Their")
print("Table 6 calculations below are unaffected (Z = 1 and 2 only) - reported,")
print("not repaired.")

# the sign of eq. (36): the printed "-405 mu1 Q1" cannot be what they computed
eta_wrong = (lambda es, a, b: b * es + (27 - 135 * a - 6 * b) * es ** 2
             + (-54.0 - 405.0 * a + 9.0 * b) * es ** 3)(
                 (1 / np.tanh(3.0) - 1 / 3.0) / 3.0, mu1_series(1.0), 3.0)
print(f"\ntheir eq. (36) prints q3 = -54 - 405 mu1 Q1 + 9 mu2 Q2 (600 ppi crop);")
print(f"with that sign, first order at Z = 1, phi = 3 gives eta_a = "
      f"{eta_wrong:.3f} -")
print(f"negative - where their Table 6 prints 0.653. The matching conditions")
print(f"they state force +405 (asserted in pz_coeffs above), and with it their")
print(f"whole printed eta_a column returns, as the next cell shows. A sign")
print(f"misprint, provable from the paper's own numbers.")'''))

cells.append(code(r'''# ---- their Table 6 (Z = 1), both columns, reproduced ------------------------
t6 = pan6[pan6["Z"] == 1.0]
print("Pan & Zhu Table 6 at Z = 1 (printed values loaded, never retyped):")
print(f"{'phi':>6}{'eta_a printed':>15}{'eq.32 here':>12}{'eta* printed':>14}"
      f"{'eq.37 here':>12}")
DEV_ETA_A = DEV_ETA_S = 0.0
for _, row in t6.iterrows():
    ea = eta_pz(row["phi"], 1.0, 1.0, 1.0)          # first order: Q1 = Q2 = 1
    es = eta_series_1st(row["phi"], 1.0)
    DEV_ETA_A = max(DEV_ETA_A, abs(ea - row["eta_a_printed"]))
    DEV_ETA_S = max(DEV_ETA_S, abs(es - row["eta_star_printed"]))
    print(f"{row['phi']:6.1f}{row['eta_a_printed']:15.3f}{ea:12.4f}"
          f"{row['eta_star_printed']:14.3f}{es:12.4f}")
print(f"worst |printed - recomputed|: eta_a {DEV_ETA_A:.4f}, eta* {DEV_ETA_S:.4f}"
      f" (3-decimal rounding is 0.0005)")

# what their own printed numbers say about the 1.5% claim, over ALL
# transcribed pairs (both Z columns)
gaps = pan6.assign(gap=(pan6["eta_star_printed"] - pan6["eta_a_printed"]).abs())
gaps["rel"] = gaps["gap"] / gaps["eta_star_printed"]
PAN_WORST_ABS = float(gaps["gap"].max())
tied = gaps[np.isclose(gaps["gap"], PAN_WORST_ABS)]
j = gaps["rel"].idxmax()
PAN_WORST_REL = float(gaps.loc[j, "rel"])
print(f"\ntheir claim: deviation < 1.5%. Their worst printed ABSOLUTE gap is a "
      f"tie at {PAN_WORST_ABS:.3f}")
print("  (" + "; ".join(f"Z = {r.Z:g}, phi = {r.phi:g}: {r.eta_a_printed:.3f} vs "
                        f"{r.eta_star_printed:.3f}" for r in tied.itertuples())
      + "),")
print(f"  i.e. {100*PAN_WORST_ABS:.1f}% of the eta scale. The worst RELATIVE "
      f"pair over all {len(gaps)} transcribed")
print(f"  pairs is the phi = {gaps.loc[j, 'phi']:g} member of that tie: "
      f"{PAN_WORST_REL:.2%} (Z = {gaps.loc[j, 'Z']:g}). Their claim holds only")
print(f"  in the absolute reading, already for first order - the relative "
      f"reading exceeds")
print(f"  the advertised bound {PAN_WORST_REL/0.015:.1f}-fold on their own "
      f"printed numbers.")'''))

cells.append(code(r'''# ---- the exact 2-D reference, and the out-of-sample second-order test -------
class Cyl2D:
    """Finite cylinder: radial (nu=1) x axial (Cartesian), c = 1 on both
    exposed faces, symmetry at r = 0 and midplane. Aspect Z = H/(2 Rp)."""
    def __init__(self, nr=100, Z=1.0):
        nz = max(int(round(nr * Z)), 20)
        self.shape = (nr, nz)
        jac, g = 0, 0
        for ax, (npts, ext, nu) in enumerate(((nr, 1.0, 1), (nz, Z, 0))):
            x_f = np.linspace(0.0, ext, npts + 1)
            x_c = 0.5 * (x_f[:-1] + x_f[1:])
            gm, gb = construct_grad(self.shape, x_f, x_c, BC, axis=ax)
            dv = construct_div(self.shape, x_f, nu=nu, axis=ax)  # nu=1 radial, 0 axial
            jac = jac + dv @ (-gm)
            g = g + dv @ (-gb)
            if ax == 0:
                wr = np.diff(x_f ** 2)               # r dr measure
            else:
                wz = np.diff(x_f)
        self.jac_diff, self.g_const = jac, g.toarray().reshape(-1, 1)
        self.numjac = NumJac(self.shape + (1,))
        w = np.multiply.outer(wr, wz)
        self.w = (w / w.sum()).ravel()

    def eta(self, rf, phi):
        def residual(cflat):
            c2 = cflat.reshape(self.shape + (1,))
            gr, jr = self.numjac(lambda cc: phi ** 2 * rf.r(cc), c2)
            return (self.jac_diff @ cflat.reshape(-1, 1) + self.g_const
                    + gr.reshape(-1, 1), self.jac_diff + jr)
        sol = newton(residual, np.ones(int(np.prod(self.shape))), tol=1e-12,
                     maxfev=40,
                     callback=lambda x, dx: clip_approach(x, dx, lower_bounds=0.0))
        res_norm = float(np.max(np.abs(residual(sol.x)[0])))
        assert res_norm < 1e-9, f"unconverged 2-D Newton residual {res_norm:.1e}"
        return float(self.w @ rf.r(sol.x))

CYL = Cyl2D(100, 1.0)
PHIS = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
CYL_VS_SERIES = max(abs(CYL.eta(Power(1), p) - eta_series_1st(p)) /
                    eta_series_1st(p) for p in PHIS)
print(f"pymrm 2-D cylinder vs the exact first-order series over phi = 0.5-8: "
      f"worst rel {CYL_VS_SERIES:.1e}")

MU2 = 2.0 + 1.0 / 1.0          # V/S = Rp/mu2 at Z = 1 (Aris/B1.2 length)
def eta_bischoff_slab(phi, n):
    """Fig. 1 used literally: SLAB curve at m built with L = V/S."""
    return CURVES[float(n)].eta(phi / MU2 * np.sqrt((n + 1.0) / 2.0))

def eta_aris_bischoff(phi):
    """The two 1965-era collapses combined: the cylinder's own first-order
    curve read at equal generalised modulus (phi_1st = phi*sqrt((n+1)/2))."""
    return eta_series_1st(phi * np.sqrt(1.5), 1.0)

def dev_fun(kind, phi):
    ex = CYL.eta(Power(2), phi)
    if kind == "pz":
        return abs(eta_pz(phi, 1.0, 2.0, np.sqrt(2.0 / 3.0)) - ex)   # Q1=2, Q2=sqrt(2/3)
    if kind == "ab":
        return abs(eta_aris_bischoff(phi) - ex)
    if kind == "ab_rel":
        return abs(eta_aris_bischoff(phi) - ex) / ex
    return abs(eta_bischoff_slab(phi, 2) - ex)

def max_dev(kind, lo=0.3, hi=8.0, ngrid=14):
    pp = np.geomspace(lo, hi, ngrid)
    vals = [dev_fun(kind, p) for p in pp]
    i = int(np.argmax(vals))
    assert 0 < i < ngrid - 1, "deviation maximum at the window edge"
    res = minimize_scalar(lambda lp: -dev_fun(kind, np.exp(lp)),
                          bracket=(np.log(pp[i - 1]), np.log(pp[i]),
                                   np.log(pp[i + 1])),
                          method="brent", options={"xtol": 1e-6})
    return np.exp(res.x), -res.fun

PHI_PZ, DEV_PZ = max_dev("pz")
PHI_BI, DEV_BI = max_dev("bi")
PHI_AB, DEV_AB = max_dev("ab")
_, REL_AB = max_dev("ab_rel")          # the relative maximum, also root-found
Q1_N2, Q2_N2, Q3_N2 = pz_coeffs(1.0, 2.0, np.sqrt(2.0 / 3.0))

# first-order versions, for reference (pure shape error / their tested case)
def dev1(kind, phi):
    ex = eta_series_1st(phi)
    v = eta_pz(phi, 1.0, 1.0, 1.0) if kind == "pz" else eta_bischoff_slab(phi, 1)
    return abs(v - ex)
PZ1 = max(dev1("pz", p) for p in PHIS)
BI1 = max(dev1("bi", p) for p in PHIS)

print(f"\nsecond order (r = c^2), Z = 1, exact = 2-D pymrm; deviations ABSOLUTE in eta:")
print(f"  Pan & Zhu eq. (32) ('arbitrary kinetics, < 1.5%'):  max "
      f"{DEV_PZ:.3f} at phi = {PHI_PZ:.2f}")
print(f"  slab Fig.-1 curve at m(V/S) (shape ignored):        max "
      f"{DEV_BI:.3f} at phi = {PHI_BI:.2f}")
print(f"  cylinder 1st-order curve at equal m (Aris+Bischoff): max "
      f"{DEV_AB:.3f} at phi = {PHI_AB:.2f} ({REL_AB:.1%} relative - inside the 15% band)")
print(f"for first order the slab reading errs {BI1:.3f} (pure shape error, the")
print(f"B1.2 band) and Pan & Zhu err {PZ1:.3f} (their tested case)")
print(f"Pan & Zhu's matched coefficients at Q1 = 2, Q2 = sqrt(2/3): "
      f"q1 = {Q1_N2:.2f}, q2 = {Q2_N2:.2f}, q3 = {Q3_N2:.2f}")'''))

cells.append(code(r'''fig, ax = plt.subplots(figsize=(7.4, 4.6))
pp = np.geomspace(0.3, 8.0, 40)
ex2 = [CYL.eta(Power(2), p) for p in pp]
ax.semilogx(pp, ex2, "k-", lw=2.0, label="exact (2-D pymrm)")
ax.semilogx(pp, [eta_pz(p, 1.0, 2.0, np.sqrt(2.0 / 3.0)) for p in pp],
            "-", color="tab:red", lw=1.7, label="Pan & Zhu eq. (32)")
ax.semilogx(pp, [eta_bischoff_slab(p, 2) for p in pp],
            "-", color="tab:blue", lw=1.7, label="slab Fig.-1 curve at $m(v_p/s_x)$")
ax.semilogx(pp, [eta_aris_bischoff(p) for p in pp],
            "-", color="tab:green", lw=1.7,
            label="cylinder 1st-order curve at equal $m$ (Aris + Bischoff)")
ax.set_xlabel(r"Thiele modulus $\phi = R_p\sqrt{k_v C_o/D}$ [-]")
ax.set_ylabel(r"$\eta$ [-]")
ax.set_title("Second order, finite cylinder $Z=1$: three shortcut estimates against exact")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""The audit separates cleanly along the two error axes.
**Pan & Zhu's polynomial handles the shape exactly and fails on the
kinetics**: advertised for arbitrary kinetics with "maximum deviation of less
than 1.5%", it is excellent where they tested it (first order: 0.004 worst
here, matching their Table 6) but for a plain second-order reaction its
midrange deviation reaches 0.090 absolute in $\eta$ - six times their claimed
bound even in the generous absolute reading - because the
cubic-in-$\eta_{sp}$ interpolation is controlled only at its two ends, and
$Q_1 = 2$ puts large oscillating coefficients ($q_2 \approx -8$,
$q_3 \approx +29$) between them. **Reading Bischoff's Figure 1 literally
fails on the shape**: the slab curve at $m(v_p/s_x)$ errs by up to 0.099 for
$n = 2$ and 0.107 already for first order - that is the shape band `B1.2`
measured, not a kinetics error at all. **The combination both 1965 papers
point to - the shape's own first-order curve, read at Bischoff's generalised
modulus - beats everything**: 0.030 absolute, 4.4% relative, comfortably
inside the advertised 15% band, with no fitted constants and no polynomial.
The 33-year-newer approximation with the 10x tighter claim is the one that
breaks; the page states this only for the kinetics actually tested here
($n = 2$, $Z = 1$).

Two further printed defects surfaced in the secondary source along the way,
both proven from the paper's own numbers above: Table 5's $\mu_2$ at
$Z = 2.5$ (prints 2.50; their eq. 30 gives 2.40 exactly), and eq. (36)'s
$-405\mu_1Q_1$ (the matching conditions they state force $+405$; the printed
sign makes their own Table 6 column negative). Reported, not repaired - the
amended sign is used, labelled as an inference forced by their stated
construction."""))

cells.append(md(r"""### Cross-check against the published quartet pages

`B1.1` ships the exact isothermal $\eta(\phi)$ for Thiele's three geometries.
For a first-order slab Bischoff's $m$ *is* Thiele's $\phi$ (eq. 17 at
$n = 1$), so that page's slab table must coincide with this page's $n=1$
curve - a conventions check (length in the half-width, no factor $\sqrt 2$
adrift), structural with respect to physics since both sides are the same
closed form. The residual few-1e-6 is the quantisation of the CSV's stored
$\phi$ values (6-7 significant digits), not a disagreement."""))

cells.append(code(r'''b11 = load_data("isothermal-exact.csv", page="B1.1-thiele-weisz-hicks")
slab_rows = b11[b11["geometry"] == "slab"]
B11_DEV = max(abs(row["eta"] - CURVES[1.0].eta(row["phi"])) / row["eta"]
              for _, row in slab_rows.iterrows() if 0.1 <= row["phi"] <= 20)
print(f"B1.1 slab table vs this page's n = 1 curve at equal phi = m: "
      f"{len(slab_rows)} rows, worst rel dev {B11_DEV:.1e}")
print("(B1.1's page flags nothing about these rows; its non-isothermal caveats")
print(" concern its continuation branch, which this page does not touch.)")

# export the exact curves computed here (deterministic grid, so re-execution
# reproduces the file byte-identically)
import csv
out_rows = [("rate_form", "param", "m", "eta")]
M_EXPORT = np.round(np.geomspace(0.1, 20.0, 61), 10)
for n in (0.0, 0.5, 1.0, 2.0, 3.0):
    for m in M_EXPORT:
        out_rows.append(("power", f"{n:g}", f"{m:.10g}", f"{CURVES[n].eta(m):.10f}"))
for z in (0.01, 0.1, 1.0, 10.0, 100.0):
    for m in M_EXPORT:
        out_rows.append(("lhhw", f"{z:g}", f"{m:.10g}", f"{LCURVES[z].eta(m):.10f}"))
for w in W_FAM:
    for m in M_EXPORT:
        out_rows.append(("volchange", f"{w:g}", f"{m:.10g}", f"{VCURVES[w].eta(m):.10f}"))
try:
    with open(Path("data") / "eta-generalised-modulus.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(out_rows)
    print(f"wrote data/eta-generalised-modulus.csv ({len(out_rows)-1} rows)")
except OSError:
    print("data/ not writable here (e.g. Colab) - export skipped")'''))

# ---------------------------------------------------------------- validation
cells.append(md(r"""## Validation

No measurement exists to validate against - the honest ceiling for this
source is agreement between independent computations, and every check below
can fail (the break table demonstrates each one failing).

| # | check | result |
|---|---|---|
| 1 | Route B (quadrature of eqs. 14-15) vs Route C (printed closed forms, eqs. 21-28 with the n=1/2 amendment), 4 orders x 10 centre concentrations | 2e-12 |
| 2 | both vs eq. (20), the paper's own third form | 2e-13 |
| 3 | pymrm finite volume vs Route B, 8 rate forms x 3 moduli, plus the regularised dead-zone case | 2e-6 / 4e-6 |
| 4 | grid refinement of Route A, observed order | 2.01 |
| 5 | dead-zone onset: closed form m* = (n+1)/(1-n) vs Route B extrapolated in u^(1/4) | 5e-14 |
| 6 | both headline band numbers recomputed on the second route, extrema re-root-found | 5e-7 |
| 7 | Bischoff Table 1, eight cells vs exact eqs. (33)-(34) | 7 in rounding, 1 defect |
| 8 | eq. (14) as printed (denominator exponent -1/2) vs its parents eqs. (8)+(13) and printed eq. (21) | off by 1/(n+1) in m; 5th printed defect |
| 9 | Pan & Zhu Table 5 (mu1 series, mu2 closed form) and Table 6 (both columns) reproduced | rounding, 2 defects |
| 10 | pymrm 2-D cylinder vs the exact first-order series | 6e-5 |
| 11 | B1.1's published slab table vs this page's n = 1 curve | 3e-6, the stored precision; structural |

**The break table.** Every metric reported to CI is moved by at least one
deliberately injected defect below, or is explicitly labelled structural.
A perturbation shows sensitivity, never correctness - which is why the two
headline numbers are also computed twice on routes sharing no assembly
(check 6), and why the extrema are root-found rather than sampled (row x
shows what sampling alone would have reported)."""))

cells.append(code(r'''BREAKS = {}

# (i) modulus prefactor dropped - the standard modulus masquerading as m
BREAKS["(i) standard modulus as abscissa"] = (
    f"band 15%-claim {BAND_15:.4f} -> {BAND_STD15:.4f}; "
    f"band 30%-claim {BAND_30_KINK:.4f} -> {BAND_STD30:.4f}")

# (ii) elliptic modulus mistranscribed: k = sin 45 in the n=2 reduction
def broken_n2(u):
    R = 1.0 / u
    m = 1.5 * np.sqrt(R) * ellipkinc(phi_from_cos(R), K45SQ) / 3.0 ** 0.25
    return m
BRK_II = max(abs(broken_n2(u) - mb_eta(Power(2), u)[0]) / mb_eta(Power(2), u)[0]
             for u in U_GRID)
BREAKS["(ii) k=sin45 in eq. (23)"] = f"route B-vs-C dev {ROUTE_BC_DEV:.1e} -> {BRK_II:.2f}"

# (iii) phi branch clamped at pi/2 (arcsin of the printed sin-phi, no continuation)
def clamped_n2(u):
    R = 1.0 / u
    sinphi = np.sqrt(4 * SQ3 * (R - 1.0)) / (SQ3 - 1.0 + R)
    m = 1.5 * np.sqrt(R) * ellipkinc(np.arcsin(min(sinphi, 1.0)), K15SQ) / 3 ** 0.25
    return m
BRK_III = max(abs(clamped_n2(u) - mb_eta(Power(2), u)[0]) / mb_eta(Power(2), u)[0]
              for u in U_GRID)
BREAKS["(iii) phi branch not continued past pi/2"] = (
    f"route B-vs-C dev {ROUTE_BC_DEV:.1e} -> {BRK_III:.2f} (u < 0.366 only)")

# (iv) wrong geometry: the slab solved as a sphere
class SlabWrongNu(Slab):
    def __init__(self, n=400):
        super().__init__(n)
        self.div_mat = construct_div(self.shape, self.x_f, nu=2, axis=0)  # sphere!
        # rebuild weights consistently with the volume measure
        self.w = np.diff(self.x_f ** 3) / 1.0
        self.w = self.w / self.w.sum()
BRK_IV = abs(SlabWrongNu(400).eta(Power(2), 1.0) - CURVES[2.0].eta(1.0)) / CURVES[2.0].eta(1.0)
BREAKS["(iv) nu=2 instead of nu=0"] = f"pymrm dev {PYMRM_DEV:.1e} -> {BRK_IV:.2f}"

# (v) eta read from a one-sided last-cell gradient instead of the rate integral
def eta_lastcell(slab, rf, m):
    Phi = m * np.sqrt(2.0 * rf.Rint(0.0, 1.0))
    c, _ = slab.solve(rf, Phi)
    flux = (1.0 - c.ravel()[-1]) / (0.5 / slab.n)   # O(h) one-sided read
    return flux / Phi ** 2
E_EX2 = CURVES[2.0].eta(1.0)
BRK_V = abs(eta_lastcell(SLAB, Power(2), 1.0) - E_EX2) / E_EX2
BRK_V_COARSE = abs(eta_lastcell(Slab(100), Power(2), 1.0) - E_EX2) / E_EX2
BRK_V_ORDER = np.log2(BRK_V_COARSE / BRK_V) / 2.0
BREAKS["(v) one-sided last-cell flux read"] = (
    f"pymrm dev {PYMRM_DEV:.1e} -> {BRK_V:.1e}, and the measured order drops "
    f"from {ORDER:.2f} to {BRK_V_ORDER:.2f}")

# (vi) eq. (31) log-bracket omitted when converting Chu-Hougen's M to m: each
# curve then sits at modulus m*sqrt(bracket_z), the fan never closes, and the
# large-m spread is the closed-form sqrt(bracket_max/bracket_min) - 1
brk = [1.0 - z * np.log1p(1.0 / z) for z in ZSTARS]
BRK_VI = np.sqrt(max(brk) / min(brk)) - 1.0
BREAKS["(vi) eq. (31) log-bracket dropped"] = (
    f"LHHW band {BAND_L:.4f} -> a non-closing fan, asymptotic spread {BRK_VI:.1f}")

# (vii) D(alpha) dropped from the modulus integral in the volume-change family:
# the collapse then fails asymptotically by |sqrt(I0_w/0.5) - 1|, forever
BRK_VII = max(abs(np.sqrt(VolChange(w).Rint(0.0, 1.0) / 0.5) - 1.0) for w in W_FAM)
BREAKS["(vii) D(alpha) dropped from eq. (13)"] = (
    f"volume-change collapse {VC_DEV:.4f} -> asymptotic deviation {BRK_VII:.4f}, "
    f"never closing")

# (viii) Pan-Zhu alpha_m taken as zeros of J1 instead of J0
def mu1_wrong(Z, M=200, N=400):
    am = jn_zeros(1, M)[:, None]
    n = np.arange(N)[None, :]
    th = (2 * n + 1) * np.pi / (2.0 * Z)
    return 32.0 / np.pi ** 2 * float(np.sum(
        1.0 / (am ** 2 * (2 * n + 1) ** 2 * (am ** 2 + th ** 2))))
BRK_VIII = max(abs(mu1_wrong(Z) - v) for Z, v in zip(pan5["Z"], pan5["mu1_printed"]))
BREAKS["(viii) alpha_m from J1 zeros"] = f"mu1 dev {MU1_DEV:.1e} -> {BRK_VIII:.1e}"

# (ix) Pan-Zhu q3 with the printed -405 sign
def eta_pz_printed_sign(phi, Z, Q1, Q2):
    a, b = mu1_series(Z) * Q1, (2.0 + 1.0 / Z) * Q2
    es = (1.0 / np.tanh(phi) - 1.0 / phi) / phi
    return b * es + (27 - 135 * a - 6 * b) * es ** 2 + (-54 - 405 * a + 9 * b) * es ** 3
BRK_IX = max(abs(eta_pz_printed_sign(row["phi"], 1.0, 1.0, 1.0) - row["eta_a_printed"])
             for _, row in t6.iterrows())
BREAKS["(ix) eq. (36) sign as printed"] = f"Table-6 eta_a dev {DEV_ETA_A:.4f} -> {BRK_IX:.2f}"

# (x) extrema sampled on the coarse grid instead of root-found
BRK_X = max(spread(m, KEYS_15) for m in np.geomspace(0.2, 10.0, 15))
BREAKS["(x) band max sampled, 15-pt grid"] = (
    f"band {BAND_15:.6f} -> {BRK_X:.6f} (under-read {BAND_15 - BRK_X:.1e})")

# (xi) V/S length wrong by the classic factor (Rp instead of Rp/mu2)
def dev_bi_wrongL(phi):
    return abs(CURVES[2.0].eta(phi * np.sqrt(1.5)) - CYL.eta(Power(2), phi))
BRK_XI = max(dev_bi_wrongL(p) for p in PHIS)
BREAKS["(xi) L = Rp instead of vp/sx"] = f"Bischoff-estimate dev {DEV_BI:.3f} -> {BRK_XI:.3f}"

# (xii) deliberate null: outer quadrature tolerance loosened 1e5x
def mb_eta_loose(rf, u):
    I0 = rf.Rint(0.0, 1.0)
    Iout = quad(lambda s: float(rf.d(u + s * s)) * 2.0 * s /
                np.sqrt(max(rf.Rint(u, u + s * s), 1e-300)),
                0.0, np.sqrt(1.0 - u), epsabs=1e-8, epsrel=1e-6, limit=400)[0]
    m = 0.5 * Iout / np.sqrt(I0)
    return m, (1.0 / m) * np.sqrt(rf.Rint(u, 1.0) / I0)
BRK_XII = max(abs(mb_eta_loose(Power(n), u)[1] - mb_eta(Power(n), u)[1])
              for n in (0.5, 2.0) for u in (0.1, 0.5, 0.9))
BREAKS["(xii) quad tolerance 1e-6 (null row)"] = (
    f"eta moves {BRK_XII:.1e} - the band numbers are quadrature-converged")

print("defect-injection table:")
for k, v in BREAKS.items():
    print(f"  {k:42s} {v}")'''))

cells.append(code(r'''METRICS = {
    # the collapse, measured (headlines)
    "band_nhalf_to_3_max_relspread": float(BAND_15),
    "band_nhalf_to_3_argmax_m": float(M_15),
    "band_nhalf_to_3_frac_of_upper": float(BAND_15_FRAC),
    "band_incl_n0_at_kink_relspread": float(BAND_30_KINK),
    "band_incl_n0_frac_of_upper": float(BAND_30_FRAC),
    "band_incl_n0_argmax_m": float(M_30),
    "eta_n3_at_m1": float(ETA3_AT_1),
    "band_second_route_gap": float(SECOND_ROUTE_GAP),
    "band_asprinted_nhalf_would_give": float(BAND_BAD),
    # the standard-modulus contrast
    "std_modulus_band_nhalf_to_3": float(BAND_STD15),
    "std_modulus_band_incl_n0": float(BAND_STD30),
    "std_modulus_asympt_spread_incl_n0": float(STD_ASYM_30),
    # route agreement
    "routeB_vs_closed_forms_worst_rel": float(ROUTE_BC_DEV),
    "eq20_vs_routeB_worst_rel": float(EQ20_DEV),
    "mstar_nhalf_extrapolated_rel_dev": float(MSTAR_DEV),
    "asymptote_eta_m_dev_at_m30": float(ASYM_DEV),
    # pymrm
    "pymrm_slab_worst_rel": float(PYMRM_DEV),
    "pymrm_slab_observed_order": float(ORDER),
    "pymrm_deadzone_rel_dev": float(DEAD_DEV),
    "pymrm_cyl2d_vs_series_worst_rel": float(CYL_VS_SERIES),
    # LHHW and volume change
    "lhhw_band_max_relspread": float(BAND_L),
    "lhhw_band_argmax_m": float(M_L),
    "lhhw_std_modulus_band": float(BAND_LS),
    "volchange_dev_from_n1_curve": float(VC_DEV),
    "volchange_min_margin_frac_of_band": float(MARGIN_FRAC),
    # printed tables and defects
    "table1_eq33_three_cells_worst_rel": float(T1_33_OK),
    "table1_eq33_defect_cell_rel": float(T1_33_DEFECT),
    "table1_eq34_worst_rel": float(T1_34_WORST),
    "eq27_as_printed_worst_rel_m": float(ASP_DEV),
    "eq27_amended_worst_rel_m": float(AMD_DEV),
    "eq14_as_printed_worst_rel_m": float(EQ14_ASP_DEV),
    # secondary source
    "pan_mu1_vs_table5_worst_abs": float(MU1_DEV),
    "pan_mu2_z25_defect_abs": float(MU2_DEV[2.5]),
    "pan_table6_eta_a_reproduced_worst_abs": float(DEV_ETA_A),
    "pan_table6_eta_star_reproduced_worst_abs": float(DEV_ETA_S),
    "pan_printed_worst_pair_abs": float(PAN_WORST_ABS),
    "pan_printed_worst_pair_rel": float(PAN_WORST_REL),
    "pan_zhu_n2_out_of_sample_max_abs": float(DEV_PZ),
    "slab_curve_estimate_cyl_n2_max_abs": float(DEV_BI),
    "slab_curve_estimate_cyl_1st_max_abs": float(BI1),
    "aris_bischoff_combined_cyl_n2_max_abs": float(DEV_AB),
    "aris_bischoff_combined_cyl_n2_max_rel": float(REL_AB),
    "pan_zhu_cyl_1st_max_abs": float(PZ1),
    # cross-page
    "b11_slab_reconciliation_worst_rel": float(B11_DEV),
}
report_agreement("B1.3", METRICS)

# ---- break-row coverage: every metric names the row that moves it -----------
BAND = "(i) directly; (x) shows the sampling under-read; (xii) null"
ROUTES = "(ii)+(iii) break one route and the gap explodes"
COVERAGE = {
    "band_nhalf_to_3_max_relspread": BAND,
    "band_nhalf_to_3_argmax_m": "(i) moves it to the M-axis optimum; (x) class",
    "band_nhalf_to_3_frac_of_upper": BAND,
    "band_incl_n0_at_kink_relspread": BAND,
    "band_incl_n0_frac_of_upper": BAND,
    "band_incl_n0_argmax_m": "(x) class: a grid without m=1 misses the kink",
    "eta_n3_at_m1": "(ii)/(iii) class through Route C; (i) moves the abscissa",
    "band_second_route_gap": ROUTES,
    "band_asprinted_nhalf_would_give":
        "itself the injected defect (the as-printed eq. 27 curve)",
    "std_modulus_band_nhalf_to_3": "(i) is this number's own definition; see note",
    "std_modulus_band_incl_n0": "(i) likewise",
    "std_modulus_asympt_spread_incl_n0":
        "STRUCT: exact algebra sqrt(2)-1, ratio of eq.-17 prefactors; recorded, "
        "cannot fail short of a transcription break covered by (i)",
    "routeB_vs_closed_forms_worst_rel": ROUTES,
    "eq20_vs_routeB_worst_rel": ROUTES + "; (xii) null",
    "mstar_nhalf_extrapolated_rel_dev": "(xii) moves it; derivation is this page's",
    "asymptote_eta_m_dev_at_m30":
        "STRUCT on the n<1 dead-zone branch AND the n=1 leg (m_lim ~ 23.7 < 30 "
        "puts eta(30) on the 1/m branch by construction); real for n=2,3 where "
        "(ii)/(iii) class breaks would move it",
    "pymrm_slab_worst_rel": "(iv) and (v) directly",
    "pymrm_slab_observed_order": "(v) drops it to ~1",
    "pymrm_deadzone_rel_dev": "(iv) class through the shared operators; the "
        "eps-sequence in its cell shows the regularisation limit being taken",
    "pymrm_cyl2d_vs_series_worst_rel": "(iv) class (radial/axial nu swap breaks it)",
    "lhhw_band_max_relspread": "(vi) directly",
    "lhhw_band_argmax_m": "(vi) class",
    "lhhw_std_modulus_band": "(vi) is its inverse (the fan eq. 31 removes)",
    "volchange_dev_from_n1_curve": "(vii) directly",
    "volchange_min_margin_frac_of_band": "(vii) shrinks it",
    "table1_eq33_three_cells_worst_rel": "exponent break of eq. (33) class (i)/(ii);"
        " any mistranscription moves all three cells at once",
    "table1_eq33_defect_cell_rel":
        "recorded defect size; moved by the same transcription breaks",
    "table1_eq34_worst_rel": "exponent 0.7 mistranscription moves all four cells",
    "eq27_as_printed_worst_rel_m": "recorded defect size; (xii) null on the arbiter",
    "eq27_amended_worst_rel_m": ROUTES,
    "eq14_as_printed_worst_rel_m":
        "recorded defect size, exactly n/(n+1) at n = 3 (u-independent constant "
        "factor); the (i)/(ii) transcription-break class moves it",
    "pan_mu1_vs_table5_worst_abs": "(viii) directly",
    "pan_mu2_z25_defect_abs":
        "STRUCT: recorded printed-defect size |2.50 - 2.40|; exact by definition",
    "pan_table6_eta_a_reproduced_worst_abs": "(ix) directly; (viii) too",
    "pan_table6_eta_star_reproduced_worst_abs": "(viii) directly",
    "pan_printed_worst_pair_abs":
        "STRUCT: arithmetic on loaded printed cells (a tie: Z=2 at phi=4 and "
        "phi=6); a transcription error in the CSV is the only thing that could "
        "move it",
    "pan_printed_worst_pair_rel":
        "STRUCT: same arithmetic, the phi=6 member of the tie (worst relative "
        "over all transcribed pairs)",
    "pan_zhu_n2_out_of_sample_max_abs": "(ix) moves it; (iv) class through the reference",
    "slab_curve_estimate_cyl_n2_max_abs": "(xi) directly",
    "slab_curve_estimate_cyl_1st_max_abs": "(xi) class",
    "aris_bischoff_combined_cyl_n2_max_abs":
        "(i) class (the sqrt(1.5) map is the injected quantity); (iv) class "
        "through the 2-D reference",
    "aris_bischoff_combined_cyl_n2_max_rel": "same rows as its absolute twin",
    "pan_zhu_cyl_1st_max_abs": "(viii)+(ix)",
    "b11_slab_reconciliation_worst_rel":
        "STRUCT: both sides are tanh(m)/m; tests only the loading and the "
        "modulus convention, and says so",
}
missing, extra = set(METRICS) - set(COVERAGE), set(COVERAGE) - set(METRICS)
assert not missing and not extra, (missing, extra)
print("break-row coverage (which injected defect moves each metric):")
for k in METRICS:
    print(f"  {k:42s} {COVERAGE[k]}")
ABS_FLOOR = 1e-12
below = sorted(k for k, v in METRICS.items() if abs(v) < ABS_FLOOR)
print(f"\nmetrics below check_agreement's ABS_FLOOR = 1e-12 (outside the CI")
print(f"regression suite, protected only by this page's asserts): {below or 'none'}")
print("note: std_modulus_band_* are the value row (i) produces - they are the")
print("quantified injected defect, kept as the paper's own comparison")'''))

cells.append(md(r"""**What this page cannot conclude.** The band numbers are
for Bischoff's stated family - power laws $n \in \{0, 1/2, 1, 2, 3\}$ plus
the LHHW and volume-change forms - on a slab; a rate form outside it (negative
order, inhibition maxima, multiple reactions) can leave the band, and
`B1.7`'s Hudgins/LHHW cells show related criteria failing exactly there. The
claim "the calculation behind Figure 1 used the correct $n=1/2$ reduction" is
an inference from the printed 15% figure, not a reading of the figure itself.
The Pan & Zhu out-of-sample failure is established for $n = 2$ at $Z = 1$
only; their method may do better for kinetics closer to first order (their
steam-shift test evidently did). And nothing here validates the *physics*
against nature - every route is the same isothermal model."""))

# ------------------------------------------------------------ what pymrm adds
cells.append(md(r"""## What pymrm adds

**For the collapse itself: nothing, deliberately.** The exact curves come
from quadrature of Bischoff's own parametric solution; pymrm's slab solves
reproduce them to ~2e-6 and exist to establish that the discretised BVP -
the thing a reader will actually reuse - carries the same physics, including
the non-smooth dead-zone case ($n<1$ past $m^*$, solved through a
convergence-restoring regularisation refined to its limit) and the
variable-diffusivity case (Picard on the face diffusivity around the Newton
loop).

**What needed a solver:** the out-of-sample audit. The exact second-order
finite-cylinder reference - the case that breaks Pan & Zhu's 1.5% claim and
shows the Aris + Bischoff combination staying inside its advertised band -
is a 2-D nonlinear pymrm solve (`nu=1` radial x Cartesian axial, `NumJac`,
10^4 unknowns) that neither 1965 nor 1998 could produce cheaply: Bischoff
worked with elliptic-integral tables, and Pan & Zhu called finite elements
for their one nonlinear test. It is the same operator assembly as `B1.2`'s
finite cylinder with a Newton loop dropped in - which is exactly the reuse
path that page's Reuse section promised.

**And the measurement discipline.** Root-found extrema (row x of the break
table shows a 15-point sweep under-reading the band), a second independent
route to both headlines, and every printed table recomputed beside its
transcription - which is how one slide-rule cell, one compositor's mangling
of eq. (27), one flipped exponent in eq. (14), one copied $\mu_2$ cell and
one sign in eq. (36) surfaced from two papers that have been cited for
decades."""))

# ------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

**To estimate $\eta$ for your kinetics:** normalise the rate so
$\tilde r(1)=1$, compute $m$ from eq. (13) (one closed-form or numerical
integral), and read the first-order curve $\eta = \tanh m/m$. For any rate
form inside Bischoff's family your error is bounded by the measured band -
under 15% of the answer if the local order stays in $[1/2, 3]$, under ~29%
if it approaches zero order - and vanishes in both limits. The `Curve` class
here gives the exact answer for arbitrary $r(C)$ and $D(C)$ at the cost of
two quadratures; `Slab` solves the BVP when you need the profile too.

**Three traps, all live on this page:**

- **The $\sqrt{(n+1)/2}$ lives inside $m$.** Feed eq. (13) the standard
  modulus and the collapse silently evaporates (break row i) - at large $m$
  the curves then disagree by $\sqrt2 - 1 \approx 41\%$ forever.
- **Do not implement eqs. (27)-(28) as printed.** The $n=1/2$ elliptic
  reduction carries a misprint (prefactor and $\phi$-argument); as printed it
  is up to 41% off in $m$. Use eq. (20) by quadrature, or the amended form in
  this page's source. Eq. (14) as printed carries a flipped denominator
  exponent too - $m$ comes out a factor $n+1$ low for order $n$.
- **The $\sin\phi$ parametrisations cross $\phi = \pi/2$** (at
  $C_o/C_L = 1+\sqrt3$ for $k=\sin15°$); evaluate $\phi$ from the $\cos\phi$
  form or continue the branch, else the $n=2$ curve is wrong for
  $C_L/C_o < 0.37$ (break row iii).

**A dead-zone onset formula the paper does not print:** for $r = k_vC^n$,
$n<1$, the centre concentration reaches zero at exactly
$m^* = (n+1)/(1-n)$ - derived here from eq. (14), reducing to Wheeler's
$m^*=1$ at $n=0$. Beyond it, $\eta = 1/m$ exactly.

**Pellet shapes:** combine this page's modulus with `B1.2`'s length
($L = v_p/s_x$) - but read the **shape's own** first-order curve, not
Figure 1's slab curves. On the second-order finite cylinder the combined
estimate lands at 0.030 absolute (4.4% relative) in $\eta$; reading the slab
curve instead costs ~0.10, essentially all of it the shape band `B1.2`
measured. If that is not good enough, the 2-D `Cyl2D` class here is the same
four operator calls with your rate form dropped in.

**Secondary source caution:** Pan & Zhu's eq. (32) is excellent near first
order and unreliable for strongly non-first-order kinetics ($n=2$: ~0.09
absolute, 6x their claimed bound); if you use it, re-derive $q_2, q_3$ from
the matching conditions - their printed eq. (36) has the sign of the
$405\mu_1Q_1$ term wrong.

**Related pages.** [`B1.1`](../B1.1-thiele-weisz-hicks/) (Thiele's curve,
Weisz-Hicks nonisothermal); [`B1.2`](../B1.2-aris-shape-factor/) (the shape
collapse this page leans on); `B1.4` (Weisz-Prater criterion); `B1.7` (the
transport-limit criteria, including LHHW cells); `C1.1` (LHHW kinetics with
real data).

**Cite the source, not this page:** Bischoff, K. B., *Effectiveness Factors
for General Reaction Rate Forms*, AIChE Journal **11**(2) 351-355 (1965),
doi:10.1002/aic.690110229. The secondary comparison: Pan, T. & Zhu, B.,
*Study on diffusion-reaction process inside a cylindrical catalyst pellet*,
Chemical Engineering Science **53**(5) 933-946 (1998),
doi:10.1016/S0009-2509(97)00385-0."""))

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
