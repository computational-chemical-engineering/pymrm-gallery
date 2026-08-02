#!/usr/bin/env python3
"""Generate index.ipynb for page G1.8. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

# ---------------------------------------------------------------- front matter
cells.append(md(r"""---
title: "Trickle-bed reactor with partial catalyst wetting"
description: "A partially wetted pellet fed from two sides at once, and a printed figure whose legend is off by one row — reconstructed, quantified, and shown both ways."
categories: [sec:G, struct:S8, tier:T1, data:tier6, phase:gas-liquid-solid]
date: 2026-08-02
---

# Trickle-bed reactor with partial catalyst wetting

**Catalog ID:** `G1.8` · **Structures:** `S8` (diffusion–reaction in a particle
with mixed external resistances) · **Tier:** T1

In a trickle bed the liquid does not cover the packing. Part of every pellet is
wetted and part is dry, and a **volatile** reactant can reach the dry part
straight from the gas — bypassing the gas–liquid film and the liquid–solid film
that stand between it and the wetted part. Partial wetting is then not a
penalty at all: it is a *shortcut*, and Herskowitz & Smith's Figure 6 measures
how large a one.

This page reproduces that figure from the paper's own equations, and in doing so
reports **two printed defects** in the source: a cross-reference to the wrong
table, and a legend whose rows are misaligned by one against the curves they
label. The second one is why the figure looks unreproducible until it is fixed.

**This is not experimental validation.** Figure 6 contains no data — the four
lines on it are the authors' own computed model output. Reproducing them tests a
transcription and a chain of algebra. The page is **tier 6** and says so
everywhere."""))

# ------------------------------------------------------------------ background
cells.append(md(r"""## Background

Herskowitz & Smith (1983) is a review of trickle-bed reactors, and its
partial-wetting section is the piece the gallery wants. The physical picture:

- The external surface of a catalyst pellet is split into a **wetted** fraction
  $f_e$ (covered by flowing liquid) and a **dry** fraction $1-f_e$ (covered at
  most by a stagnant film, and in contact with the gas).
- A limiting reactant that is **non-volatile** can only arrive through the
  liquid. Partial wetting then costs external area, and the effect on the
  overall effectiveness factor is modest — the paper says so explicitly for
  $0.6 < f_e < 1$.
- A limiting reactant that **is volatile** can arrive at the dry surface
  directly from the gas. The gas-side external resistance is normally far
  smaller than the gas–liquid *plus* liquid–solid resistance in series, so the
  dry patch is the *better* supplied one. Partial wetting then **raises** the
  overall effectiveness factor, and the more strongly the more diffusion-limited
  the pellet is.

Figure 6 quantifies the second case in the limit where the gas-side resistance
is negligible altogether. Its ordinate is

$$
\chi = \frac{\eta_o - \eta_o(f_e=1)}{\eta_o(f_e=1)} \tag{19}
$$

the *fractional gain* in overall effectiveness factor relative to a fully wetted
pellet, plotted against the Thiele modulus $\phi$ for four liquid mass
velocities $L_m$. Everything on that figure is positive and rises with $\phi$:
at $\phi \approx 45$ the gain reaches an order of magnitude.

That is the practical warning the section is built around. If you measure a rate
in a trickle bed and vary the liquid flow to test for mass-transfer limitation,
the two mechanisms move in **opposite** directions — lowering $L_m$ dries more
surface (helping a volatile reactant) while starving the wetted surface
(hurting it). The paper's own words: doing so "may lead to pitfalls."

### The gallery's structural claim

`S8` here is the classical pellet problem of `B1.1` with a *mixed* external
boundary condition: the same interior equation, two different Robin conditions
on two parts of one surface. The paper handles that by superposition of two
whole-surface problems (its Eq. 6, from Ramachandran & Smith 1979), which is an
**approximation** — and the page is careful to say that nothing here tests it."""))

# ------------------------------------------------------------- published model
cells.append(md(r"""## The published model

Read from a 600 dpi render of the paper (the PDF text layer mangles every
equation on these pages). Page numbers are journal pages.

**The particle problem** (page 4). For an isothermal, irreversible, $m$-th order
reaction in a porous particle,

$$
\nabla^2 C - \phi^2 C^m = 0, \qquad
\phi = \frac{V_p}{S_x}\sqrt{\frac{\rho_p k_v C_L^{m-1}}{D_e}}, \qquad
C = C_i/C_L \tag{3}
$$

with $\nabla$ non-dimensionalised on the **characteristic length $V_p/S_x$** —
for a sphere, $R/3$. The external surface is divided into a wetted part $S_w$
and a "dry" part $S_D$, carrying different Robin conditions:

$$
\frac{\partial C}{\partial n} = \alpha_{gLs}\,(C^*_L - C) \ \ \text{on } S_w,
\qquad
\frac{\partial C}{\partial n} = \alpha_{gs}\,(1 - C) \ \ \text{on } S_D
\tag{4a, 4b}
$$

with $C^*_L = C_{L,b}/C_L$, and $\alpha_{gLs}$ the gas–liquid and liquid–solid
resistances in series (Eq. 4d). The overall effectiveness factor is the volume
average of the local rate, $\eta_o = V_p^{-1}\int_{V_p} C^m\,dV$ (Eq. 2a).

**Table 1, sphere row** (page 4). For $m = 1$, superposing a fully wetted and a
fully dry particle in the proportion $f_e : 1-f_e$ (Eq. 6) gives a closed form:

$$
\eta_o = \frac{f_e C^*_L}
              {\dfrac{\phi^2}{\alpha_{gLs}} + \dfrac{3\phi^2}{3\phi\coth 3\phi - 1}}
       + \frac{1 - f_e}
              {\dfrac{\phi^2}{\alpha_{gs}} + \dfrac{3\phi^2}{3\phi\coth 3\phi - 1}}
\tag{Table 1}
$$

**The two closures Figure 6 uses** (page 8):

$$
f_e = 0.77\,L_m^{0.1}
\tag{20}
$$

$$
\frac{1}{\alpha_{gLs}} = \frac{1}{1.05\,L_m^{0.3}} + \frac{1}{6.91\,L_m^{0.6}}
\tag{21}
$$

Eq. 21 is stated to be Eq. 4d with the Goto & Smith (1975) and Specchia et al.
(1978) correlations substituted at the operating conditions printed in the
figure caption. The figure's own two stated closing assumptions are
$\alpha_{gs} \to \infty$ ("very large") and $C^*_L = 1.0$.

### Printed defect 1: the wrong table is cited

Page 8 reads: *"The overall effectiveness factor is calculated from the
approximate solution for a spherical particle given in **Table 2**."* Table 2
holds the pressure-drop constants $\beta$ and $\gamma$. The spherical solution is
the **third row of Table 1**, page 4. Checked on 600 dpi renders of both pages
and against the paper's notation list. This is a typographic slip with no effect
on any number; it is recorded because a reader following the citation finds
nothing usable.

### The chain collapses

With $\alpha_{gs}\to\infty$ the dry term loses its external resistance and
becomes the plain sphere effectiveness factor
$\eta_s(\phi) = (3\phi\coth 3\phi - 1)/(3\phi^2)$; with $C^*_L = 1$ the wetted
term is exactly $\eta_o(f_e=1)$. So Table 1 reads
$\eta_o = f_e\,\eta_o(f_e{=}1) + (1-f_e)\,\eta_s$, and Eq. 19 collapses to

$$
\chi = (1 - f_e)\,\eta_s(\phi)\,\frac{\phi^2}{\alpha_{gLs}}
\tag{collapse}
$$

a product of one factor per assumption: how much surface is dry, how deep the
pellet is, and how large the resistance the dry patch bypasses. This is used
below as a *reading aid and a root-finding target*; it is an algebraic identity,
not evidence, and the page labels it as such."""))

# ---------------------------------------------------- parameters & assumptions
cells.append(md(r"""## Parameters and assumptions

Everything numeric on this page traces to one of three places, and the page
never mixes them:

| | source | what it is |
|---|---|---|
| $f_e$, $\alpha_{gLs}$ | Eqs. 20, 21, page 8 | printed closures, transcribed from a 600 dpi render |
| $\eta_o$ | Table 1 sphere row, page 4 | printed closed form, transcribed the same way |
| $\alpha_{gs}\to\infty$, $C^*_L = 1$ | figure caption + page 8 text | stated assumptions |
| the four curve positions | digitised, `data/` | the **only** extracted quantity |
| $L_m \approx 10$ for the fourth curve | reconstructed on this page | **labelled reconstruction**, see below |

The caption's conditions ($Sc_L = 100$, $d_p = 1.5\times10^{-3}$ m,
$\varepsilon_p = 0.5$, $\varepsilon_B = 0.4$, $\tau = 2$,
$\mu_L = 8\times10^{-4}$ kg/(m·s)) are what produced the constants in Eq. 21.
They do not enter this page again: Eq. 21 already contains them, and
re-deriving it would mean transcribing two further correlations from papers not
on disk. That is a deliberate limit and it is listed as a blind spot at the end.

Assumptions inherited, all of them the paper's:

1. **First order**, $m = 1$ — Table 1 exists only for $m=1$.
2. **Isothermal.** The paper separately argues $T_s - T_f < 0.5$ K under the
   relevant conditions.
3. **Eq. 6 superposition**: the partially wetted particle is a weighted average
   of two whole-surface particles. This is the paper's approximation to a
   genuinely multi-dimensional boundary-value problem, and *nothing on this
   page tests it*.
4. **$\alpha_{gs}$ literally infinite.** Section "Validation" measures how much
   that matters — quite a lot."""))

# ------------------------------------------------------------ environment cell
cells.append(md(r"""### Environment"""))

cells.append(code("""try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml"""))

cells.append(code('''import sys, urllib.request
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
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from IPython.display import Markdown, display

from pymrm import construct_grad, construct_div, NumJac, newton, clip_approach
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "G1.8-trickle-bed-partial-wetting"
np.seterr(all="ignore")
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# Deviation convention, used EVERYWHERE on this page and nowhere reversed:
#     dev = (model - figure) / figure
def dev(model, figure):
    return (np.asarray(model, float) - np.asarray(figure, float)) / np.asarray(figure, float)

COL = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]   # one colour per drawn curve'''))

# --------------------------------------------------------------- paper algebra
cells.append(md(r"""### The paper's algebra, transcribed

Four functions, one per printed equation. They contain no fitted quantity."""))

cells.append(code('''def f_e(L_m):
    """Eq. 20, page 8:  f_e = 0.77 L_m^0.1  (wetting efficiency)."""
    return 0.77 * np.asarray(L_m, float) ** 0.1


def alpha_gLs(L_m):
    """Eq. 21, page 8:  1/a_gLs = 1/(1.05 L_m^0.3) + 1/(6.91 L_m^0.6).

    Two external resistances in series (gas-liquid, then liquid-solid), so the
    RECIPROCALS add - which is why the equation is printed for 1/a_gLs.
    """
    L = np.asarray(L_m, float)
    return 1.0 / (1.0 / (1.05 * L ** 0.3) + 1.0 / (6.91 * L ** 0.6))


def eta_sphere(phi):
    """Classical sphere effectiveness factor on the V/S length scale.

    Table 1's sphere row contains 3 phi^2 / (3 phi coth 3 phi - 1), which is
    1/eta_sphere: with L = V/S = R/3 the Thiele modulus on the RADIUS is 3 phi.
    """
    lam = 3.0 * np.asarray(phi, float)
    return (lam / np.tanh(lam) - 1.0) / (3.0 * np.asarray(phi, float) ** 2)


def eta_o_table1(phi, L_m, alpha_gs=np.inf, C_star_L=1.0):
    """Table 1, sphere row (page 4) - NOT Table 2, which the paper cites."""
    inv_eta = 1.0 / eta_sphere(phi)
    a = alpha_gLs(L_m)
    fe = f_e(L_m)
    wet = fe * C_star_L / (phi ** 2 / a + inv_eta)
    dry = (1.0 - fe) / (phi ** 2 / alpha_gs + inv_eta)
    return wet + dry


def chi_table1(phi, L_m, alpha_gs=np.inf, C_star_L=1.0):
    """Eq. 19 evaluated on Table 1's sphere row."""
    # denominator: the same row at f_e = 1 and C*_L = 1, i.e. the wetted term alone
    eta_full = 1.0 / (phi ** 2 / alpha_gLs(L_m) + 1.0 / eta_sphere(phi))
    return eta_o_table1(phi, L_m, alpha_gs, C_star_L) / eta_full - 1.0


def chi_collapsed(phi, L_m):
    """The closed form for a_gs -> inf and C*_L = 1:  (1-f_e) eta_s phi^2 / a_gLs.

    ALGEBRAIC IDENTITY with chi_table1(..., alpha_gs=inf). Agreement between the
    two is guaranteed and is used below only as a transcription check on the
    collapse, never as evidence about the physics.
    """
    return (1.0 - f_e(L_m)) * eta_sphere(phi) * phi ** 2 / alpha_gLs(L_m)'''))

cells.append(code('''# --- check the transcription of Eq. 20 against the numbers printed IN the figure
L_legend = np.array([0.50, 1.0, 2.0, 7.0])
fe_printed = np.array([0.72, 0.77, 0.83, 0.94])
fe_eq20 = f_e(L_legend)

print("Eq. 20 against the f_e column printed inside Figure 6")
print("  L_m    printed   Eq. 20    rounded   match")
for L, p, c in zip(L_legend, fe_printed, fe_eq20):
    print(f"  {L:<6.2f} {p:<9.2f} {c:<9.4f} {round(c, 2):<9.2f} {'OK' if abs(round(c,2)-p) < 5e-3 else 'MISMATCH'}")

# and the collapse, as the identity it is
_phi = np.logspace(0, 1.7, 25)
identity = float(np.max(np.abs(chi_collapsed(_phi, 2.0) / chi_table1(_phi, 2.0) - 1.0)))
print(f"\\ncollapse vs full Table 1 chain (algebraic identity): max rel diff {identity:.2e}")'''))

# -------------------------------------------------------------------- the data
cells.append(md(r"""## The data

**One dataset, and it is not measurements.** Figure 6 carries four drawn lines
and no markers whatsoever; the lines are the authors' own model output. The CSV
holds, per drawn curve, a straight line fitted to it in $(\log\phi, \log\chi)$:
its slope, its position expressed as $\chi$ at $\phi = 10$, the traced $\phi$
range, and the rms residual of the fit.

Three consequences the page holds to throughout:

1. **Tier 6, never "validated".** Reproducing a computed curve tests
   transcription and algebra. It cannot test the model against nature.
2. **Away from the fit's own parametrisation there are no raw points.** Any
   deviation quoted at a $\phi$ other than through these four fitted lines is
   *model against a fitted line*, and is labelled that way.
3. **The reduction to a straight line is itself lossy** — quantified in
   Validation, where the model's own best straight-line fit turns out to have an
   rms of the same order as the digitisation's.

The maintainer reviewed Figure 6 in a private decision artifact on 2026-08-02
and endorsed the conclusion drawn from these four fits. What was reviewed was
**the reading** — the figure alongside the reconstruction below — not a
marker-by-marker audit, which Figure 6 does not admit. That review is what
licenses the use of these rows; the numbers in the CSV are unchanged from the
2026-07-29 extraction."""))

cells.append(code('''fig6 = load_data("herskowitz-smith-1983-fig6.csv", page=PAGE)
meta6 = load_meta("herskowitz-smith-1983-fig6.csv", page=PAGE)
print(cite_data(meta6))
print()
print(fig6.to_string(index=False))
print()
print("estimated error:", meta6["acquisition"]["estimated_error"])

CHI10 = fig6["chi_at_phi_10"].to_numpy()          # position of each drawn curve
SLOPE = fig6["slope_dlogchi_dlogphi"].to_numpy()  # its fitted log-log slope
PHI_LO = float(fig6["phi_min"].iloc[0])
PHI_HI = float(fig6["phi_max"].iloc[0])
RMS_FIT = fig6["fit_rms_decades"].to_numpy()


def line_chi(phi, i):
    """The i-th DIGITISED FITTED LINE evaluated at phi. Not raw pixels."""
    return CHI10[i] * (np.asarray(phi, float) / 10.0) ** SLOPE[i]'''))

# -------------------------------------------------------- pymrm implementation
cells.append(md(r"""## PyMRM implementation

The closed form above is an answer, not a model. The model is Eq. 3 with the
Robin conditions 4a/4b, and that is what pymrm solves — one spherical
finite-volume BVP per external condition, combined by the paper's Eq. 6.

Two conventions carry all the risk here, and both are about the **length scale**:

- The paper's $\phi$ is built on $V_p/S_x = R/3$, so on a grid $u = r/R$ the
  interior equation is $\nabla_u^2 C = (3\phi)^2 C^m$.
- The paper's $\partial C/\partial n$ in Eqs. 4a/4b is *also* per unit $V_p/S_x$,
  so $\partial C/\partial n = \tfrac13\,\partial C/\partial u$ and the Robin
  condition in grid units is $\partial C/\partial u + 3\alpha C = 3\alpha C^*$.

Both factors of 3 are asserted, not printed, so both are broken on purpose in
Validation to show the check notices.

House conventions followed:

- **Outward normal.** At $u=0$ the outward normal points into the pellet, so
  symmetry is `{a:1, b:0, d:0}`. At $u=1$ the Robin condition
  $\partial C/\partial n = \alpha(C^*-C)$ becomes
  $a\,\partial C/\partial u + bC = d$ with `a=1, b=3α, d=3αC*`; the
  $\alpha\to\infty$ limit is the Dirichlet `{a:0, b:1, d:C*}`.
- **Constant operators assembled once**, cached per boundary condition, never
  inside Newton.
- **Layout `(n_u, 1)`** — spatial axis first, single field last. `NumJac` couples
  the last axis in full, so a bare `(n_u,)` would declare every cell coupled to
  every other and build a dense Jacobian. `axes_diagonals` is *not* passed: at
  `ndims = 1` it would be reinterpreted as absolute indices and leave the
  Jacobian with no diagonal, and the source term here is pointwise anyway.
- `nu=2` in `construct_div` is the spherical geometry factor."""))

cells.append(code('''class WettedSphere:
    """Isothermal spherical pellet, Eq. 3 with a single Robin surface condition.

    Grid u = r/R on [0, 1]. The paper's length scale is V_p/S_x = R/3, so
        d/dn|paper = (1/3) d/du       and       lap|paper = (1/9) lap_u,
    which is where every factor of 3 below comes from. `n_scale` exists only so
    that the Validation section can set it wrong on purpose.
    """

    def __init__(self, n_u=800, nu=2, n_scale=3.0):
        self.n_u, self.shape, self.nu, self.n_scale = n_u, (n_u, 1), nu, n_scale
        self.u_f = np.linspace(0.0, 1.0, n_u + 1)              # faces
        self.u_c = 0.5 * (self.u_f[:-1] + self.u_f[1:])        # centres
        self.div = construct_div(self.shape, self.u_f, nu=nu)  # nu: 0 slab, 1 cyl, 2 sphere
        vf = self.u_f ** (nu + 1) / (nu + 1.0)                 # cell volumes for this nu
        self.dv, self.v_tot = np.diff(vf), vf[-1]
        self.numjac = NumJac(self.shape)   # (n_u, 1): last axis is the FIELD -> diagonal block
        self._ops = {}

    def operators(self, alpha, C_star):
        key = (float(alpha), float(C_star))
        if key not in self._ops:
            s = self.n_scale
            # --- boundary conditions, OUTWARD normal --------------------------
            #   u = 0 (centre) : symmetry,                      dC/dn = 0
            #   u = 1 (surface): dC/dn = alpha (C* - C) with n in units of V/S,
            #                    i.e. (1/s) dC/du + alpha C = alpha C*
            left = {"a": 1.0, "b": 0.0, "d": 0.0}
            if np.isinf(alpha):
                right = {"a": 0.0, "b": 1.0, "d": C_star}          # alpha -> inf
            else:
                right = {"a": 1.0, "b": s * alpha, "d": s * alpha * C_star}
            g, gb = construct_grad(self.shape, self.u_f, self.u_c, (left, right))
            lap = self.div @ g
            lap_bc = (self.div @ gb).toarray().reshape(-1, 1)
            self._ops[key] = (lap, lap_bc)
        return self._ops[key]

    def solve(self, phi, alpha, C_star=1.0, m=1.0, maxfev=60):
        lap, lap_bc = self.operators(alpha, C_star)
        lam2 = (self.n_scale * phi) ** 2                        # lap_u C = (3 phi)^2 C^m
        # clip at 0: the maximum principle puts C in [0, max(1, C*)], and a Newton
        # iterate dipping negative would make C**m complex for non-integer m
        src = lambda y: -lam2 * np.clip(y, 0.0, None) ** m

        def residual(y):
            y = y.reshape(self.shape)
            g_s, j_s = self.numjac(src, y)
            return (lap @ y.reshape((-1, 1)) + lap_bc
                    + np.asarray(g_s).reshape((-1, 1))), lap + j_s

        hi = max(1.0, C_star)
        r = newton(residual, np.full((self.n_u, 1), min(1.0, C_star)),
                   maxfev=maxfev, tol=1e-13,
                   callback=lambda x, g: clip_approach(x, g, 0.0, hi))
        y = r.x.ravel()
        res, _ = residual(y)
        return y, float(np.max(np.abs(res)))

    def eta(self, y, m=1.0):
        """eta = V^-1 int C^m dV  (Eq. 2a)."""
        return float(np.sum(np.clip(y, 0.0, None) ** m * self.dv) / self.v_tot)

    # ---------------------------------------------------------------- the model
    def chi(self, phi, L_m, m=1.0, alpha_gs=np.inf, C_star_L=1.0):
        """Eq. 19 from three finite-volume solves and the Eq. 6 superposition."""
        a, fe = float(alpha_gLs(L_m)), float(f_e(L_m))
        eta_wet = self.eta(self.solve(phi, a, C_star_L, m)[0], m)   # wetted path
        eta_dry = self.eta(self.solve(phi, alpha_gs, 1.0, m)[0], m)  # dry path
        eta_ref = self.eta(self.solve(phi, a, 1.0, m)[0], m)         # f_e = 1, C*_L = 1
        eta_o = fe * eta_wet + (1.0 - fe) * eta_dry                  # Eq. 6
        return eta_o / eta_ref - 1.0


sph = WettedSphere(n_u=800)
_y, _r = sph.solve(10.0, float(alpha_gLs(2.0)))
print(f"pymrm sphere, n_u = {sph.n_u}: eta = {sph.eta(_y):.8f}, "
      f"max |Newton residual| = {_r:.2e}")'''))

cells.append(md(r"""### The pymrm solve against Table 1

The finite-volume route and Table 1's closed form share **nothing** — one
discretises the BVP, the other is a formula — so this comparison can fail, and
Validation shows exactly what makes it fail. It is the check that the two length
scales and the sphere row are transcribed right."""))

cells.append(code('''PHI = np.logspace(np.log10(PHI_LO), np.log10(PHI_HI), 60)

devs = {}
for L in (0.5, 2.0, 7.0, 10.0):
    a = float(alpha_gLs(L))
    num = np.array([sph.eta(sph.solve(p, a)[0]) for p in PHI])
    ana = 1.0 / (PHI ** 2 / a + 1.0 / eta_sphere(PHI))       # Table 1 wetted path
    devs[L] = np.abs(num / ana - 1.0)
num_d = np.array([sph.eta(sph.solve(p, np.inf)[0]) for p in PHI])
dev_dirichlet = np.abs(num_d / eta_sphere(PHI) - 1.0)

ETA_MAXDEV = float(max(max(v.max() for v in devs.values()), dev_dirichlet.max()))
print("pymrm finite volume vs Table 1 sphere row, over the digitised phi window")
for L, v in devs.items():
    print(f"  Robin, a_gLs(L_m={L:<5}) = {float(alpha_gLs(L)):5.3f}:  max |dev| = {v.max():.2e}")
print(f"  Dirichlet (a_gs -> inf)          :  max |dev| = {dev_dirichlet.max():.2e}")
print(f"  worst anywhere                   :  {ETA_MAXDEV:.2e}   (at phi = {PHI[int(np.argmax(dev_dirichlet))]:.1f})")'''))

# ----------------------------------------------------------------- the finding
cells.append(md(r"""## Results

### The figure does not reproduce as printed

Evaluate the chain at the legend's four $L_m$ values and compare with where the
four drawn curves actually sit."""))

cells.append(code('''chi_at10_printed = chi_collapsed(10.0, L_legend)
d_printed = dev(chi_at10_printed, CHI10)

print("Reading the legend as printed - chi at phi = 10, dev = (model - figure)/figure")
print("  L_m     f_e      a_gLs    model     figure    dev")
for L, cm, cf, d in zip(L_legend, chi_at10_printed, CHI10, d_printed):
    print(f"  {L:<7.2f} {float(f_e(L)):<8.3f} {float(alpha_gLs(L)):<8.3f} "
          f"{cm:<9.4f} {cf:<9.4f} {d*100:+8.1f} %")'''))

cells.append(code('''display(Markdown(rf"""
The top curve lands at **{d_printed[0]*100:+.1f} %** with nothing fitted — the
transcription, the axis calibration and the whole chain are right. The other
three are out by **{d_printed[1]*100:+.0f} %**, **{d_printed[2]*100:+.0f} %** and
**{d_printed[3]*100:+.0f} %**, and not even monotonically. That is not a model
that is slightly wrong; it is a model that is right for one curve out of four.
"""))'''))

cells.append(md(r"""### What the curve *spacing* says

The **gaps between adjacent curves**, in decades, are the same measurement as
the four positions above with one degree of freedom removed — a common
multiplicative offset in $\chi$. They are therefore **not a second, independent
witness**, and this page does not offer them as one. What they do rule out is a
specific rival: that the printed reading is right and the figure's $\chi$ axis is
mis-calibrated by a constant factor. That is worth ruling out, because check 4 in
Validation turns out to be weak against exactly that error class.

What the gaps need, and what they do not. A gap in decades is invariant to where
the $\chi$ axis *starts* and to the abscissa entirely, so it is visible without
knowing where either axis begins. It is **not** calibration-free: it scales
linearly with the ordinate's decades-per-pixel, which the sidecar reports as
415.5 px/decade with residuals under 0.006 decades over $\chi$ = 0.01–20. That
one number is what the gaps rest on.

The figure's gaps are also **$\phi$-dependent**, because the four fitted lines
are not exactly parallel — their slopes differ by 3 %. The model's gaps are not
$\phi$-dependent at all, since the chain factorises. Every gap quoted below is
evaluated at **$\phi = 10$**, and the drift across the digitised window is
printed alongside so a reader can see how much that choice matters."""))

cells.append(code('''gaps_fig = np.log10(CHI10[:-1] / CHI10[1:])     # the drawn lines, at phi = 10

def gaps_line(phi):
    """The DRAWN lines' gaps at an arbitrary phi. They drift with phi because the
    four fitted slopes differ by 3 %; the model's gaps do not drift at all."""
    c = np.array([line_chi(phi, i) for i in range(4)])
    return np.log10(c[:-1] / c[1:])

def gaps_for(Ls):
    c = chi_collapsed(10.0, np.asarray(Ls, float))
    return np.log10(c[:-1] / c[1:])

g_printed = gaps_for([0.5, 1.0, 2.0, 7.0])
g_shifted = gaps_for([0.5, 2.0, 7.0, 10.0])

print("Gaps between adjacent curves, in decades, evaluated at phi = 10")
print("                        1->2    2->3    3->4")
print(f"  figure                {gaps_fig[0]:6.3f}  {gaps_fig[1]:6.3f}  {gaps_fig[2]:6.3f}")
print(f"  L_m = 0.5,1,2,7       {g_printed[0]:6.3f}  {g_printed[1]:6.3f}  {g_printed[2]:6.3f}   <- the legend")
print(f"  L_m = 0.5,2,7,10      {g_shifted[0]:6.3f}  {g_shifted[1]:6.3f}  {g_shifted[2]:6.3f}   <- shifted by one row")
print()
print(f"  the figure's 2->3 gap, {gaps_fig[1]:.3f}, is the legend model's 3->4 gap, {g_printed[2]:.3f}")
print()
print("The figure's gaps are phi-dependent (the fitted lines are not exactly parallel):")
for p in (PHI_LO, 2.0, 10.0, PHI_HI):
    gl = gaps_line(p)
    print(f"  phi = {p:<5.4g}            {gl[0]:6.3f}  {gl[1]:6.3f}  {gl[2]:6.3f}")
GAP_DRIFT = float(np.max(np.abs(gaps_line(PHI_HI) - gaps_line(PHI_LO))))
GAP_MISS_SHIFT = float(np.max(np.abs(g_shifted - gaps_fig)))
GAP_MISS_PRINT = float(np.max(np.abs(g_printed - gaps_fig)))
print(f"\\n  drift over the whole window: at most {GAP_DRIFT:.3f} decades")'''))

cells.append(code('''display(Markdown(rf"""
The legend's third gap, **{g_printed[2]:.3f}** decades, is the figure's
**second** gap, **{gaps_fig[1]:.3f}** — and the figure's other two,
{gaps_fig[0]:.3f} and {gaps_fig[2]:.3f}, appear nowhere in the printed model's
gap sequence. Only one of the three gaps is shared between the two readings, so
this is not a tidy displacement of a sequence; it is the pattern a legend column
offset by one row produces, seen from the outside, and it is visible without
knowing where either axis *starts*.

How much does it matter that the gaps are read at $\\phi = 10$? Across the whole
digitised window the figure's gaps drift by up to **{GAP_DRIFT:.3f} decades** —
*larger* than the {GAP_MISS_SHIFT:.3f} decades by which the shifted reading
misses, so it is the drift, not the shifted reading's residual, that sets how
precisely these gaps can be quoted at all. But it is
{GAP_MISS_PRINT/GAP_DRIFT:.0f} times smaller than the {GAP_MISS_PRINT:.3f}
decades by which the printed reading misses. Where the gaps are evaluated
therefore cannot change which reading wins.
"""))'''))

cells.append(md(r"""### Printed defect 2: the legend is misaligned by one row

Solve, per drawn curve, for the $L_m$ that would put it exactly where it sits.
This is a one-dimensional root find on the collapsed chain at $\phi = 10$; the
map $L_m \mapsto \chi$ is monotone decreasing, so the root is unique."""))

cells.append(code('''L_rec = np.array([brentq(lambda x, c=c: float(chi_collapsed(10.0, x)) - c,
                         1e-3, 1e4, xtol=1e-12, rtol=1e-14) for c in CHI10])

# how tightly does chi pin L_m?  d log chi / d log L_m
h = 1e-6
sens = np.array([(np.log(chi_collapsed(10.0, L * (1 + h))) - np.log(chi_collapsed(10.0, L * (1 - h))))
                 / (np.log(L * (1 + h)) - np.log(L * (1 - h))) for L in L_rec])
chi_err_pct = 100.0 * (10.0 ** RMS_FIT - 1.0)          # the fit rms, expressed in chi
L_tol = chi_err_pct / np.abs(sens)

L_shift = np.array([0.50, 2.0, 7.0, 10.0])             # 3 printed values reused + 1 reconstruction

print("Reconstructing L_m from each drawn curve's position at phi = 10")
print("  curve  legend   reconstructed   shifted reading   diff     digitisation tol.  is it a test?")
for i in range(4):
    is_test = "yes - vs a PRINTED value" if i < 3 else "no  - vs a rounded reconstruction"
    print(f"  {i+1:<6} {L_legend[i]:<8.2f} {L_rec[i]:<15.4f} {L_shift[i]:<17.2f} "
          f"{(L_rec[i]/L_shift[i]-1)*100:+6.1f} %  +/- {L_tol[i]:4.1f} %          {is_test}")
print()
print("  (tolerance = the fit's own rms residual, "
      f"{RMS_FIT.min():.3f}-{RMS_FIT.max():.3f} decades = {chi_err_pct.min():.1f}-{chi_err_pct.max():.1f} % in chi,")
print("   propagated through d log chi / d log L_m = "
      f"{sens[0]:.2f}, {sens[1]:.2f}, {sens[2]:.2f}, {sens[3]:.2f})")
print()
print("  Curve 4 carries no printed L_m, so its row compares the reconstruction with the")
print("  round number this page adopts for it. That is a rounding residual, not evidence.")'''))

cells.append(code('''within = np.abs(L_rec[:3] / L_shift[:3] - 1) * 100 < L_tol[:3]
display(Markdown(rf"""
**The reconstruction returns the paper's own printed numbers, in a different
order.** Curves 1–3 come back as
{L_rec[0]:.3f}, {L_rec[1]:.3f} and {L_rec[2]:.3f} — that is the legend's
0.50, 2.0 and 7.0 to within {abs(L_rec[0]/L_shift[0]-1)*100:.1f} %,
{abs(L_rec[1]/L_shift[1]-1)*100:.1f} % and
{abs(L_rec[2]/L_shift[2]-1)*100:.1f} %, all inside what the digitisation itself
can resolve ({L_tol[0]:.0f} %, {L_tol[1]:.0f} %, {L_tol[2]:.0f} %).
{'All three' if within.all() else 'Only some'} of the testable curves land inside
their tolerance. Curve 4 has no printed value to be tested against, so its
{abs(L_rec[3]/L_shift[3]-1)*100:.1f} % is the residual of rounding
{L_rec[3]:.2f} to 10 and carries no evidential weight either way.

So the legend's **second row, $L_m = 1.0$, belongs to no drawn curve**, and the
fourth drawn curve — at $L_m = {L_rec[3]:.2f}$, nearest round value **10** — has
no legend row. From its second row down, the legend block is offset by one
position against the curves it labels.

Two things this does **not** establish. It cannot say whether the slip is in the
legend (right curves, wrong labels) or in the plotting (right labels, curves
computed for the wrong $L_m$): both produce this figure exactly. And the value
**10** for the fourth curve is a *reconstruction* — the reconstruction says
{L_rec[3]:.2f}, and 10 is the nearest round number, which the paper's own text
uses as the practical bound ("operating at $L_m < 10$ kg/(m²)(s)"). It is
labelled as reconstructed wherever it appears, and it is the only number on this
page that is not printed somewhere in the paper.
"""))'''))

cells.append(md(r"""### Both readings, side by side

In both panels the four **thick pale lines are the drawn curves** — identically
placed left and right, because they are the figure. The **thin lines are the
pymrm solve**, coloured to match the curve the reading assigns them to.

The thing to look at in the left panel is not that the thin lines are wrong, but
*where they land*: the model at $L_m = 2$ (green) sits almost exactly on the pale
curve the legend calls $L_m = 1$, and the model at $L_m = 7$ (purple) on the one
the legend calls $L_m = 2$. Each model curve has found a drawn curve — just not
its own. That displacement by one position is the whole finding, and the right
panel is the same picture with it undone."""))

cells.append(code('''chi_num_shift = {L: np.array([sph.chi(p, L) for p in PHI]) for L in L_shift}
chi_num_print = {L: np.array([sph.chi(p, L) for p in PHI]) for L in L_legend}

fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)
for ax, (title, curves, labels) in zip(
        axes,
        [("As printed: $L_m$ = 0.50, 1.0, 2.0, 7.0", chi_num_print, L_legend),
         ("Shifted by one row: $L_m$ = 0.50, 2.0, 7.0, 10*", chi_num_shift, L_shift)]):
    for i in range(4):
        ax.loglog(PHI, line_chi(PHI, i), color=COL[i], lw=3.2, alpha=0.30, solid_capstyle="round")
    for i, L in enumerate(labels):
        ax.loglog(PHI, curves[L], color=COL[i], lw=1.5,
                  label=rf"$L_m$ = {L:g}" + ("*" if L == 10.0 else ""))
    ax.set_xlabel(r"Thiele modulus  $\\phi$")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
axes[0].set_ylabel(r"$\\chi = [\\eta_o - \\eta_o(f_e{=}1)]\\,/\\,\\eta_o(f_e{=}1)$")
fig.suptitle("Herskowitz & Smith (1983) Figure 6 - thick pale lines are the four drawn curves,\\n"
             "thin lines are the pymrm solve of Table 1's model.  * = reconstructed, not printed",
             fontsize=9)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""### Deviations, both readings

The comparison away from $\phi = 10$ is model against the **fitted straight
line**, since the CSV holds no raw points. Convention throughout:
$(\text{model} - \text{figure})/\text{figure}$."""))

cells.append(code('''PHI_TAB = np.array([2, 3, 5, 7, 10, 15, 20, 30], float)
dev_shift = np.array([dev([sph.chi(p, L) for p in PHI_TAB], line_chi(PHI_TAB, i))
                      for i, L in enumerate(L_shift)])
dev_print = np.array([dev([sph.chi(p, L) for p in PHI_TAB], line_chi(PHI_TAB, i))
                      for i, L in enumerate(L_legend)])

print("Deviation (%) against the fitted drawn lines - SHIFTED reading")
print("  phi   " + "".join(f"  L_m={L:<7g}" for L in L_shift))
for j, p in enumerate(PHI_TAB):
    print(f"  {p:<5.0f} " + "".join(f"{dev_shift[i, j]*100:+10.1f}" for i in range(4)))
print("\\nDeviation (%) against the fitted drawn lines - AS PRINTED")
print("  phi   " + "".join(f"  L_m={L:<7g}" for L in L_legend))
for j, p in enumerate(PHI_TAB):
    print(f"  {p:<5.0f} " + "".join(f"{dev_print[i, j]*100:+10.1f}" for i in range(4)))

SHIFT_WORST = float(np.max(np.abs(dev_shift)))
SHIFT_WORST_123 = float(np.max(np.abs(dev_shift[:3])))
PRINT_WORST = float(np.max(np.abs(dev_print)))
d10_shift = dev([sph.chi(10.0, L) for L in L_shift], CHI10)

# how much of curve 4's excess is the rounding of 10.18 to 10?
dev_c4_exact = dev([sph.chi(p, L_rec[3]) for p in PHI_TAB], line_chi(PHI_TAB, 3))
print("\\ncurve 4 with the reconstructed L_m = "
      f"{L_rec[3]:.2f} instead of the rounded 10:")
print("  phi   " + "".join(f"{p:8.0f}" for p in PHI_TAB))
print("  dev   " + "".join(f"{v*100:+8.1f}" for v in dev_c4_exact))
print(f"  worst {np.max(np.abs(dev_c4_exact))*100:.1f} % against "
      f"{np.max(np.abs(dev_shift[3]))*100:.1f} % at L_m = 10 exactly, i.e. most of curve 4's")
print("  excess is the rounding, not a failure of the shifted reading.")
C4_EXACT_WORST = float(np.max(np.abs(dev_c4_exact)))'''))

cells.append(code('''display(Markdown(rf"""
At $\\phi = 10$ the shifted reading gives
**{d10_shift[0]*100:+.1f} %, {d10_shift[1]*100:+.1f} %, {d10_shift[2]*100:+.1f} %**
and **{d10_shift[3]*100:+.1f} %**, against
{d_printed[0]*100:+.0f} %, {d_printed[1]*100:+.0f} %, {d_printed[2]*100:+.0f} %,
{d_printed[3]*100:+.0f} % as printed. Over $\\phi$ = 2 to 30 the shifted reading
stays inside **{SHIFT_WORST_123*100:.0f} %** for the three curves that use only
printed $L_m$ values, and inside **{SHIFT_WORST*100:.0f} %** including the
reconstructed fourth at the round $L_m = 10$; as printed the worst is
**{PRINT_WORST*100:.0f} %**.

Curve 4 is the loose one for two reasons, and neither is the shifted reading.
About **{(1 - C4_EXACT_WORST/np.max(np.abs(dev_shift[3])))*100:.0f} %** of its
excess is the rounding of {L_rec[3]:.2f} to 10: at the
reconstructed value its worst deviation falls from
{np.max(np.abs(dev_shift[3]))*100:.1f} % to **{C4_EXACT_WORST*100:.1f} %**. The
rest is the slope residual of the next section — curve 4 is the steepest drawn
line ({SLOPE[3]:.3f}), so it is where the model's slope deficit shows most. The
page keeps the round 10 in the figure, because a fitted $L_m$ would make the
comparison circular, and reports both numbers.
"""))'''))

cells.append(code('''DIG = float(chi_err_pct.mean())          # the fits' own rms, expressed as % in chi

fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.4), sharey=True)
for ax, (d, labels, title) in zip(axes, [
        (dev_print, L_legend, "As printed"),
        (dev_shift, L_shift, "Shifted by one row")]):
    for i, L in enumerate(labels):
        ax.semilogx(PHI_TAB, d[i] * 100, "o-", color=COL[i], ms=4, lw=1.3,
                    label=rf"$L_m$ = {L:g}" + ("*" if L == 10.0 else ""))
    ax.axhline(0, color="k", lw=1.0)
    ax.axhspan(-DIG, DIG, color="tab:green", alpha=0.10)
    ax.set_xlabel(r"$\\phi$")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
axes[0].set_ylabel("deviation  (model - figure)/figure  [%]")
axes[0].set_yscale("symlog", linthresh=20)
axes[1].text(2.2, DIG + 2.0, f"digitisation rms, ~{DIG:.1f} % in chi",
             fontsize=8, color="tab:green")
fig.tight_layout()
plt.show()'''))

# --------------------------------------------------------- residual slope tilt
cells.append(md(r"""### What the shift does *not* fix

The deviations under the shifted reading are not flat. They run **positive at
low $\phi$ and negative at high $\phi$**, with the same shape on all four
curves. That is a slope difference, and it survives the reassignment."""))

cells.append(code('''lp = np.log10(PHI)
A = np.vstack([lp, np.ones_like(lp)]).T
lc = np.log10(chi_collapsed(PHI, 2.0))
(model_slope, model_int), *_ = np.linalg.lstsq(A, lc, rcond=None)
model_fit_rms = float(np.sqrt(np.mean((A @ [model_slope, model_int] - lc) ** 2)))

# chi_collapsed(phi, L) = (1-f_e(L)) * [eta_s(phi) phi^2] / a_gLs(L): the phi-dependence
# is a single factor, so the log-log SLOPE is identical for every L_m.
slopes_all_L = [float(np.linalg.lstsq(A, np.log10(chi_collapsed(PHI, L)), rcond=None)[0][0])
                for L in L_shift]

def local_slope(phi):
    d = 1e-6
    return float((np.log(chi_collapsed(phi * (1 + d), 2.0)) - np.log(chi_collapsed(phi * (1 - d), 2.0)))
                 / (np.log(phi * (1 + d)) - np.log(phi * (1 - d))))

print("The model is not a straight line on these axes.")
print("  phi        1.07     2       5      10      20     45.5")
print("  local slope " + "".join(f"{local_slope(p):7.3f}" for p in (1.07, 2, 5, 10, 20, 45.5)))
print()
print(f"  best straight-line fit of the model over phi = {PHI_LO}-{PHI_HI}:")
print(f"    slope {model_slope:.4f}   rms {model_fit_rms:.4f} decades")
print(f"    (identical for all four L_m: {['%.4f' % s for s in slopes_all_L]})")
print(f"  the four digitised lines:")
print(f"    slopes {['%.4f' % s for s in SLOPE]}   mean {SLOPE.mean():.4f}")
print(f"    fit rms {['%.4f' % r for r in RMS_FIT]} decades")
SLOPE_GAP = float(SLOPE.mean() - model_slope)'''))

cells.append(code('''# The dominant uncertainty on this gap is NOT the fit rms - it is the ABSCISSA
# calibration, which this page's own digitisation never varied.  A log-log slope is
# (decades in chi)/(decades in phi), so it scales linearly with the assumed px/decade
# on the phi axis.  The CSV used 421.9 px/decade, read off the four labelled verticals
# (phi = 1, 5, 10, 50) - and those four do NOT sit on a common logarithmic ruler:
# a least-squares fit through them leaves +0.013 decades at phi = 10, i.e. that
# gridline is misplaced in the artwork by about 4 px.  So 421.9 is not the only
# defensible value, and the other three verticals support a slightly smaller one.
CAL_CSV = 421.9                       # px/decade in phi, the value behind this CSV
CAL_ALT = {                           # other readings of the same four verticals
    "4-point least squares":       418.6,
    "3-point, phi = 10 dropped":   417.8,
    "phi = 1 -> 10 baseline only": 426.0,
}
DEC = np.log10(PHI_HI / PHI_LO)

def e2e(gap):
    """A slope gap, accumulated over the digitised window, as a percentage in chi."""
    return (10.0 ** (gap * DEC) - 1.0) * 100.0

print("Sensitivity of the slope gap to the abscissa calibration")
print("  px/decade  reading of the phi verticals    mean slope   gap      end-to-end")
gaps_by_cal = {}
for lab, cal in [("as digitised (this CSV)", CAL_CSV)] + list(CAL_ALT.items()):
    s = SLOPE.mean() * cal / CAL_CSV          # slope scales linearly with px/decade
    gaps_by_cal[lab] = s - model_slope
    print(f"  {cal:8.1f}  {lab:<30} {s:8.4f}  {s-model_slope:+7.3f}  {e2e(s-model_slope):+8.0f} %")

# READ, not computed here: an independent 600 dpi re-digitisation of the same figure
# (adversarial verification pass, 2026-08-02) with its own gridline detection, its own
# calibration and 620-640 traced columns per curve.
INDEP_SLOPE_MEAN = 1.1028     # its directly measured mean fitted slope, at 418.6 px/dec
INDEP_GAP_MIN = 0.007         # smallest gap it obtained by RE-FITTING under a narrower
                              # baseline (re-fitting also moves the fit window, which
                              # rescaling this CSV's slopes cannot capture)
INDEP_GAP = INDEP_SLOPE_MEAN - model_slope
print(f"\\n  independent re-trace, measured directly at 418.6 px/decade: mean slope "
      f"{INDEP_SLOPE_MEAN:.4f}")
print(f"  the same calibration applied to THIS CSV by rescaling:      "
      f"{SLOPE.mean()*418.6/CAL_CSV:.4f}")
print("  -> the two digitisations differ by calibration alone, not by tracing.")

SLOPE_GAP_LO = float(min(min(gaps_by_cal.values()), INDEP_GAP_MIN))
SLOPE_GAP_HI = float(max(gaps_by_cal.values()))
print(f"\\n  defensible range for the gap: {SLOPE_GAP_LO:+.3f} to {SLOPE_GAP_HI:+.3f}"
      f"  ({e2e(SLOPE_GAP_LO):+.0f} % to {e2e(SLOPE_GAP_HI):+.0f} % end to end)")
print("  The SIGN survives every one of them. The SIZE does not.")'''))

cells.append(code('''end_to_end = e2e(SLOPE_GAP)
display(Markdown(rf"""
**The residual tilt, stated honestly.** The model's log–log slope is not a
single number: it falls from {local_slope(1.07):.2f} at
$\\phi = {PHI_LO}$ to {local_slope(45.5):.2f} at $\\phi = {PHI_HI}$, because
$\\eta_s\\phi^2 \\to \\phi$ only asymptotically. Its best straight-line fit over
the digitised window has slope **{model_slope:.3f}**. The four drawn curves,
reduced to straight lines by the same kind of fit, measure
**{SLOPE[0]:.3f}, {SLOPE[1]:.3f}, {SLOPE[2]:.3f}, {SLOPE[3]:.3f}** — mean
{SLOPE.mean():.3f}. The drawn lines are the steeper ones.

**The sign of the gap is solid; its size is not.** As digitised it is
{SLOPE_GAP:+.3f} in slope, {end_to_end:+.0f} % end to end over the
{DEC:.2f} decades of the window. But the dominant uncertainty is not the fit rms
— it is the **abscissa calibration**, and the abscissa is the weak axis of this
figure. Its four labelled verticals do not sit on a common logarithmic ruler:
the $\\phi = 10$ gridline is about 4 px off the line the other three define, so
the CSV's 421.9 px/decade sits roughly 1 % above the best-supported value. An
independent 600 dpi re-digitisation of the same figure, calibrating the abscissa
at 418.6 px/decade, measures a mean slope of **{INDEP_SLOPE_MEAN:.4f}** — a gap
of **{INDEP_GAP:+.3f}** ({e2e(INDEP_GAP):+.0f} % end to end) — and that is this
CSV's {SLOPE.mean():.4f} rescaled by exactly the ratio of the two calibrations.
The two traces agree; only their rulers differ. Across the defensible
calibrations the gap runs from **{SLOPE_GAP_LO:+.3f}** to **{SLOPE_GAP_HI:+.3f}**
({e2e(SLOPE_GAP_LO):+.0f} % to {e2e(SLOPE_GAP_HI):+.0f} % end to end). **The page
quotes the range, not a single number.**

**What does not defend it.** All four measured slopes lie above
{model_slope:.3f}, but that is *not* a scatter argument and this page does not
use it as one: all four were traced against **one** abscissa calibration, so a
common-mode error in that calibration moves all four together. Four numbers that
look like four samples are one measurement.

**It is unexplained, and it is not the legend misalignment.** The collapsed chain
factorises as $(1-f_e)\\,\\alpha_{{gLs}}^{{-1}}$ times $\\eta_s(\\phi)\\phi^2$, so
$L_m$ moves a curve up and down and **cannot change its slope at all** — every
$L_m$ gives the identical {model_slope:.3f}. Reassigning curves therefore has
exactly zero effect on this residual, and it would be there under any reading.

A second, smaller caveat, quantified rather than assumed away: the model's own
best straight line leaves an rms of **{model_fit_rms:.3f} decades**, comparable
to the digitisation's {RMS_FIT.min():.3f}–{RMS_FIT.max():.3f}. So the reduction
of each drawn curve to a straight line discards curvature of the same size as the
model's, and the CSV cannot distinguish a straight line from this model's curve.
The discrepancy is reported as open, with its size as a range.
"""))'''))

# ------------------------------------------------------------------ validation
cells.append(md(r"""## Validation

Ranked before any code was written, highest available first. Every check below
is followed by a **deliberate break**: a defect it is supposed to catch is
injected and the number is measured again. A check whose number does not move is
decoration, and this page names the ones that do not.

| # | check | what it is | can it fail? |
|---|---|---|---|
| 1 | Eq. 20 vs the $f_e$ printed inside the figure | printed intermediate | yes |
| 2 | pymrm finite volume vs Table 1's sphere row | independent routes | yes |
| 3 | grid convergence of the pymrm solve | discretisation | yes |
| 4 | the shift reconstruction | reconstruction vs printed values | yes |
| 5 | the collapse vs the full Table 1 chain | **algebraic identity** | **no** |
| 6 | $\alpha_{gs}\to\infty$ sensitivity | assumption audit | n/a — a sweep |

### 1–3, and what breaks them"""))

cells.append(code('''def check_eq20(coef=0.77, expo=0.1):
    """max |Eq.20 - printed f_e| over the four legend rows."""
    return float(np.max(np.abs(coef * L_legend ** expo - fe_printed)))


def check_pymrm_vs_table1(n_u=800, nu=2, n_scale=3.0, L=2.0, phis=None):
    """max relative deviation of the finite-volume eta from Table 1's sphere row."""
    s = WettedSphere(n_u=n_u, nu=nu, n_scale=n_scale)
    a = float(alpha_gLs(L))
    ps = PHI if phis is None else np.asarray(phis, float)
    num = np.array([s.eta(s.solve(p, a)[0]) for p in ps])
    ana = 1.0 / (ps ** 2 / a + 1.0 / eta_sphere(ps))
    return float(np.max(np.abs(num / ana - 1.0)))


def check_reconstruction(coef=0.77, expo=0.1, c21=(1.05, 0.3, 6.91, 0.6), chi_scale=1.0):
    """The three reconstructed L_m that the shift hypothesis predicts to be 0.5/2/7.

    Returns the max relative distance from those printed values. `chi_scale`
    corrupts the ORDINATE calibration instead of the algebra.
    """
    def a_g(L):
        return 1.0 / (1.0 / (c21[0] * L ** c21[1]) + 1.0 / (c21[2] * L ** c21[3]))

    def chi(phi, L):
        return (1.0 - coef * L ** expo) * eta_sphere(phi) * phi ** 2 / a_g(L)

    out = []
    for c, target in zip(CHI10[:3] * chi_scale, [0.5, 2.0, 7.0]):
        try:
            r = brentq(lambda x: chi(10.0, x) - c, 1e-4, 1e5, xtol=1e-12, rtol=1e-14)
        except ValueError:
            return np.inf
        out.append(abs(r / target - 1.0))
    return float(max(out))


BASE = {"eq20": check_eq20(), "pymrm": check_pymrm_vs_table1(), "recon": check_reconstruction()}
print("baselines:", {k: f"{v:.3e}" for k, v in BASE.items()})'''))

cells.append(code('''rows = []


def add(label, eq20=None, pymrm=None, recon=None, note=""):
    rows.append((label, eq20, pymrm, recon, note))


add("(baseline)", BASE["eq20"], BASE["pymrm"], BASE["recon"], "")

# --- defects that should move check 1 (Eq. 20 transcription) -----------------
add("Eq. 20 exponent 0.1 -> 0.01", eq20=check_eq20(expo=0.01),
    recon=check_reconstruction(expo=0.01), note="mis-read superscript")
add("Eq. 20 exponent 0.1 -> 0.7", eq20=check_eq20(expo=0.7),
    recon=check_reconstruction(expo=0.7), note="mis-read superscript")
add("Eq. 20 coefficient 0.77 -> 0.87", eq20=check_eq20(coef=0.87),
    recon=check_reconstruction(coef=0.87), note="mis-read digit")

# --- defects that should move check 2 (the pymrm/Table 1 pair) ---------------
add("nu = 0 (slab) not sphere", pymrm=check_pymrm_vs_table1(nu=0), note="wrong geometry")
add("nu = 1 (cylinder) not sphere", pymrm=check_pymrm_vs_table1(nu=1), note="wrong geometry")
add("Robin scale 3 -> 1 (n per R)", pymrm=check_pymrm_vs_table1(n_scale=1.0),
    note="length scale of dC/dn")
add("Robin scale 3 -> 9", pymrm=check_pymrm_vs_table1(n_scale=9.0), note="length scale of dC/dn")
add("n_u = 10", pymrm=check_pymrm_vs_table1(n_u=10), note="under-resolved grid")
add("n_u = 50", pymrm=check_pymrm_vs_table1(n_u=50), note="under-resolved grid")

# --- defects that should move check 4 (the reconstruction) -------------------
add("Eq. 21: 6.91 -> 6.31", recon=check_reconstruction(c21=(1.05, 0.3, 6.31, 0.6)),
    note="mis-read digit")
add("Eq. 21: 6.91 -> 6.71", recon=check_reconstruction(c21=(1.05, 0.3, 6.71, 0.6)),
    note="mis-read digit")
add("Eq. 21: 6.91 -> 7.91", recon=check_reconstruction(c21=(1.05, 0.3, 7.91, 0.6)),
    note="mis-read digit - IMPROVES the metric")
add("Eq. 21: 1.05 -> 1.95", recon=check_reconstruction(c21=(1.95, 0.3, 6.91, 0.6)),
    note="mis-read digit")
add("Eq. 21 exponent 0.6 -> 0.5", recon=check_reconstruction(c21=(1.05, 0.3, 6.91, 0.5)),
    note="mis-read superscript, 2nd term")
add("Eq. 21 exponent 0.6 -> 0.7", recon=check_reconstruction(c21=(1.05, 0.3, 6.91, 0.7)),
    note="mis-read superscript, 2nd term")
add("Eq. 21 exponents 0.3/0.6 -> 0.5/0.5", recon=check_reconstruction(c21=(1.05, 0.5, 6.91, 0.5)),
    note="mis-read superscripts")
add("chi axis calibration x1.05", recon=check_reconstruction(chi_scale=1.05),
    note="ordinate scale error, not a transcription defect")
add("chi axis calibration x0.95", recon=check_reconstruction(chi_scale=0.95),
    note="ordinate scale error, not a transcription defect")

w = max(len(r[0]) for r in rows)
print(f"{'injected defect':<{w}}  {'check 1':>10}  {'check 2':>10}  {'check 4':>10}   note")
print("-" * (w + 50))
for lab, a, b, c, n in rows:
    fa = f"{a:10.3e}" if a is not None else " " * 10
    fb = f"{b:10.3e}" if b is not None else " " * 10
    fc = f"{c:10.3e}" if c is not None else " " * 10
    print(f"{lab:<{w}}  {fa}  {fb}  {fc}   {n}")'''))

cells.append(code('''moved = {
    "check 1 (Eq. 20 vs printed $f_e$)": ([r[1] for r in rows[1:] if r[1] is not None], BASE["eq20"]),
    "check 2 (pymrm vs Table 1)": ([r[2] for r in rows[1:] if r[2] is not None], BASE["pymrm"]),
    "check 4 (the reconstruction)": ([r[3] for r in rows[1:] if r[3] is not None], BASE["recon"]),
}
lines = []
for k, (vals, b) in moved.items():
    ratios = sorted(v / b for v in vals if np.isfinite(v))
    if ratios[0] < 1.0:
        lines.append(f"- **{k}**: baseline {b:.2e}; the injected defects move it by factors "
                     f"{ratios[0]:.2f}x to {ratios[-1]:.0f}x — and the smallest is **below 1**, "
                     f"i.e. one defect makes the metric look *better* than the undamaged chain.")
    else:
        lines.append(f"- **{k}**: baseline {b:.2e}; every injected defect moved it, "
                     f"by factors {ratios[0]:.0f}x to {ratios[-1]:.0f}x.")

# --- check 4's real limits, measured rather than asserted --------------------
L_TOL_MIN = float(L_tol[:3].min())
recon_rows = [r for r in rows[1:] if r[3] is not None and np.isfinite(r[3])
              and not r[0].startswith("chi axis")]
unrejected = sorted((r for r in recon_rows if r[3] * 100 < L_TOL_MIN), key=lambda r: r[3])
best = unrejected[0]
# WHY: at L_m ~ 1 the first term of Eq. 21 carries almost all of 1/alpha_gLs
inv1, inv2 = 1.0 / (1.05 * 1.0 ** 0.3), 1.0 / (6.91 * 1.0 ** 0.6)
SHARE_FIRST = inv1 / (inv1 + inv2)
s_up = next(r[3] for r in rows if r[0] == "chi axis calibration x1.05")
s_dn = next(r[3] for r in rows if r[0] == "chi axis calibration x0.95")
CHECK4_BEST_DEFECT = float(best[3])
lines.append(
    f"- **but check 4 is only as sharp as the digitisation.** Its tolerance is "
    f"{L_TOL_MIN:.1f}–{L_tol[:3].max():.1f} % in $L_m$, set by the fits' own rms, against a "
    f"baseline of {BASE['recon']*100:.2f} %. **{len(unrejected)}** of the injected defects stay "
    f"under even the tightest of those bars and so would **not** be rejected: "
    + ", ".join(f"`{r[0]}` ({r[3]*100:.2f} %)" for r in unrejected) + ". The first of them, "
    f"`{best[0]}`, comes in at **{best[3]*100:.2f} %** — *better* than the undamaged "
    f"{BASE['recon']*100:.2f} %, so a larger number here is not even monotone in the size of "
    f"the defect.")
lines.append(
    f"- **The reason is structural, not statistical.** At $L_m \\\\approx 1$ the "
    f"$1.05\\\\,L_m^{{0.3}}$ term carries **{SHARE_FIRST*100:.0f} %** of "
    f"$1/\\\\alpha_{{gLs}}$, leaving the $6.91\\\\,L_m^{{0.6}}$ term only "
    f"{(1-SHARE_FIRST)*100:.0f} % — so corrupting the second term barely moves $\\\\chi$ at all. "
    f"Check 4 catches errors in Eq. 20 and in Eq. 21's **first** term; against Eq. 21's "
    f"**second** term it is essentially blind, and the page claims nothing more for it.")
lines.append(
    f"- **It is also weak against an ordinate-scale error.** Rescaling the whole $\\\\chi$ "
    f"calibration by ±5 % leaves it at {s_up*100:.1f} % / {s_dn*100:.1f} % — inside the band, "
    f"so unrejected. That is precisely the error class the curve-spacing test above *does* "
    f"cover, and it is the honest reason to keep that test even though it is not an "
    f"independent witness.")
display(Markdown("Every injected defect moves the check it is aimed at — but moving is not "
                 "rejecting:\\n\\n" + "\\n".join(lines)))'''))

cells.append(md(r"""### 3. Grid convergence

The pymrm/Table 1 comparison is only a transcription check if the
discretisation error is smaller than what it is meant to detect. Refine and
measure the order."""))

cells.append(code('''ns = [50, 100, 200, 400, 800, 1600]
errs = [check_pymrm_vs_table1(n_u=n, phis=[PHI_HI]) for n in ns]
print("worst-case grid convergence, at the largest phi in the window "
      f"(phi = {PHI_HI}, the thinnest boundary layer)")
print("  n_u     max |dev|    order")
prev = None
for n, e in zip(ns, errs):
    o = "" if prev is None else f"{np.log2(prev / e):5.2f}"
    print(f"  {n:<7} {e:.3e}   {o}")
    prev = e
GRID_ORDER = float(np.log2(errs[-2] / errs[-1]))
GRID_ORDER_MIN = float(np.min([np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]))
print(f"\\n  observed order rises from {GRID_ORDER_MIN:.2f} to {GRID_ORDER:.2f} as the surface")
print("  layer becomes resolved - consistent with a second-order scheme whose error at the")
print("  coarse end is dominated by the boundary layer, not by the interior stencil.")'''))

cells.append(md(r"""### 5. The collapse — a check that cannot fail

`chi_collapsed` and `chi_table1(..., alpha_gs=inf)` agree to machine precision.
They are **the same algebra written twice**: substituting $\alpha_{gs}=\infty$
and $C^*_L=1$ into Table 1's row and cancelling gives the collapsed form
line for line. It therefore cannot detect a mis-read coefficient in Table 1, a
wrong length scale, or anything physical — a defect in the sphere row moves both
sides identically. It is kept because it confirms the *collapse step*, which the
page's root finding relies on, and for no other reason.

The measurement below shows exactly that: a deliberately corrupted sphere row
leaves the identity untouched while the pymrm comparison (check 2) blows up."""))

cells.append(code('''def _corrupt(fac):
    """Multiply the sphere term 3 phi^2/(3 phi coth 3phi - 1) by `fac`."""
    def eta_bad(phi):
        return eta_sphere(phi) / fac
    return eta_bad


def identity_with_corrupt_sphere(fac):
    """check 5 (the collapse identity) evaluated with the sphere term corrupted."""
    eta_bad = _corrupt(fac)

    def full(phi, L):
        ie = 1.0 / eta_bad(phi)
        a, fe = float(alpha_gLs(L)), float(f_e(L))
        return (fe / (phi ** 2 / a + ie) + (1 - fe) * eta_bad(phi)) / (1.0 / (phi ** 2 / a + ie)) - 1

    def coll(phi, L):
        return (1 - float(f_e(L))) * eta_bad(phi) * phi ** 2 / float(alpha_gLs(L))

    p = np.logspace(0, 1.6, 20)
    return float(np.max(np.abs(coll(p, 2.0) / full(p, 2.0) - 1.0)))


def pymrm_vs_corrupt_sphere(fac, L=2.0):
    """check 2 evaluated against the SAME corrupted sphere term. Measured, not asserted."""
    eta_bad = _corrupt(fac)
    a = float(alpha_gLs(L))
    num = np.array([sph.eta(sph.solve(p, a)[0]) for p in PHI])
    ana = 1.0 / (PHI ** 2 / a + 1.0 / eta_bad(PHI))
    return float(np.max(np.abs(num / ana - 1.0)))


print("A corrupted sphere term, put through both checks")
print("  sphere term x   check 5 (identity)   check 2 (pymrm vs Table 1)")
for fac in (1.0, 1.01, 1.5, 3.0):
    print(f"  {fac:<15.2f} {identity_with_corrupt_sphere(fac):<20.2e} {pymrm_vs_corrupt_sphere(fac):.2e}")
print("\\n  -> check 5 does not move at all; check 2 moves by the size of the defect.")
print("     The identity is blind to the very error it looks like it is testing.")
IDENTITY = identity_with_corrupt_sphere(1.0)'''))

cells.append(md(r"""### 6. Is $\alpha_{gs}\to\infty$ safe?

The paper says only "very large". Every number on this page assumes literally
infinite. That assumption is doing real work, and the sweep says how much."""))

cells.append(code('''print("chi under a FINITE gas-solid Biot number, relative to the infinite limit, L_m = 2")
print("  a_gs      phi=2      phi=10     phi=30")
for ags in (10, 30, 100, 300, 1000, 3000):
    r = [chi_table1(p, 2.0, alpha_gs=ags) / chi_collapsed(p, 2.0) - 1 for p in (2.0, 10.0, 30.0)]
    print(f"  {ags:<9} " + "".join(f"{v*100:+9.2f}%" for v in r))
AGS_1000_AT30 = float(chi_table1(30.0, 2.0, alpha_gs=1000.0) / chi_collapsed(30.0, 2.0) - 1)
AGS_100_AT30 = float(chi_table1(30.0, 2.0, alpha_gs=100.0) / chi_collapsed(30.0, 2.0) - 1)
print(f"\\nfor scale, alpha_gLs itself is only {float(alpha_gLs(0.5)):.2f} to "
      f"{float(alpha_gLs(10.0)):.2f} over the four curves")'''))

cells.append(code('''display(Markdown(rf"""
"Very large" has to mean *very*: at $\\phi = 30$ even
$\\alpha_{{gs}} = 1000$ still costs **{abs(AGS_1000_AT30)*100:.1f} %** in
$\\chi$, and $\\alpha_{{gs}} = 100$ — already about
{100/float(alpha_gLs(2.0)):.0f} times $\\alpha_{{gLs}}$ — costs
**{abs(AGS_100_AT30)*100:.0f} %**. The reconstruction above is
therefore **conditional on $\\alpha_{{gs}} = \\infty$ exactly**, and would return
different $L_m$ values under any finite value. The page states this as a
condition of its finding, not as a caveat buried in the code.
"""))'''))

cells.append(md(r"""### The alternative explanations, and why they lose

Before concluding the legend is misaligned, the two hypotheses that would rescue
the printed reading are given their best chance. Each is fitted *per curve*, so
each has four free parameters against four curves and must fit perfectly — the
question is whether what it needs is physically admissible."""))

cells.append(code('''print("ALT A - a finite alpha_gs rescuing the printed reading (fitted per curve at phi = 10)")
alt_a = []
for L, c in zip(L_legend, CHI10):
    try:
        r = brentq(lambda a: float(chi_table1(10.0, L, alpha_gs=a)) - c, 1e-3, 1e9,
                   xtol=1e-10, rtol=1e-14)
        alt_a.append(r)
        print(f"  L_m = {L:<5} needs alpha_gs = {r:10.1f}")
    except ValueError:
        alt_a.append(np.nan)
        print(f"  L_m = {L:<5} no admissible root")

print("\\nALT B - a corrected f_e per curve (printed reading kept)")
alt_b = []
for L, c in zip(L_legend, CHI10):
    r = brentq(lambda f: (1 - f) * float(eta_sphere(10.0)) * 100.0 / float(alpha_gLs(L)) - c,
               -5.0, 1 - 1e-12)
    alt_b.append(r)
    print(f"  L_m = {L:<5} needs f_e = {r:.4f}   (Eq. 20 gives {float(f_e(L)):.4f}, "
          f"figure prints {fe_printed[list(L_legend).index(L)]:.2f})")
alt_a, alt_b = np.array(alt_a), np.array(alt_b)

# ALT B pressed harder: could ANY single power law c L^n supply those f_e values?
# A power law has a CONSTANT log-log slope, so measure the required one per interval.
lb, lf = np.log(L_legend), np.log(alt_b)
ALTB_SLOPES = np.diff(lf) / np.diff(lb)
N_PL = float((lf[-1] - lf[0]) / (lb[-1] - lb[0]))    # pinned to the two extreme values,
C_PL = float(alt_b[0] / L_legend[0] ** N_PL)         # which is the kindest single power law
pl = C_PL * L_legend ** N_PL
ALTB_PL_MISS = (pl / alt_b - 1.0)
print("\\nIs the required f_e(L_m) a power law at all?")
print("  L_m                     " + "".join(f"{L:9.2f}" for L in L_legend))
print("  required f_e            " + "".join(f"{v:9.4f}" for v in alt_b))
print("  d log f_e / d log L_m       " + "".join(f"{v:9.3f}" for v in ALTB_SLOPES)
      + "   <- must be CONSTANT for a power law")
print(f"  best power law {C_PL:.3f} L^{N_PL:.3f}" + "".join(f"{v:9.4f}" for v in pl))
print("  it misses by            " + "".join(f"{v*100:+8.1f}%" for v in ALTB_PL_MISS))
ALTB_PL_WORST = float(np.max(np.abs(ALTB_PL_MISS)))'''))

cells.append(code('''display(Markdown(rf"""
**ALT A fails on consistency.** One $\\alpha_{{gs}}$ is stated for the whole
figure, but the four curves need
{alt_a[0]:.0f}, {alt_a[1]:.1f}, {alt_a[2]:.1f} and {alt_a[3]:.1f} — spread over
a factor of {np.nanmax(alt_a)/np.nanmin(alt_a):.0f}, not monotone in $L_m$, and
three of the four far too small to be called "very large" when
$\\alpha_{{gLs}}$ itself is only {float(alpha_gLs(0.5)):.2f}–{float(alpha_gLs(7.0)):.2f}.

**ALT B fails on two counts, and the first one is decisive on its own.** The
$f_e$ values it needs — {alt_b[0]:.4f}, {alt_b[1]:.4f}, {alt_b[2]:.4f},
{alt_b[3]:.4f} — **contradict the $f_e$ column printed inside the figure**
({fe_printed[0]:.2f}, {fe_printed[1]:.2f}, {fe_printed[2]:.2f},
{fe_printed[3]:.2f}), which Eq. 20 already reproduces to the printed precision.
The paper would have to disagree with its own legend.

It also fails on **curvature**. It is tempting to argue from the *correction* to
$f_e$ ({(alt_b[0]-float(f_e(0.5)))*100:+.1f},
{(alt_b[1]-float(f_e(1.0)))*100:+.1f}, {(alt_b[2]-float(f_e(2.0)))*100:+.1f},
{(alt_b[3]-float(f_e(7.0)))*100:+.1f} percentage points) being non-monotone, but
that proves nothing: a non-monotone *difference* between two monotone functions
is perfectly ordinary, and the required $f_e$ values themselves *are* monotone.
What actually excludes a power law is that a power law has a **constant** log–log
slope, and the required values do not: theirs run
**{ALTB_SLOPES[0]:.3f} → {ALTB_SLOPES[1]:.3f} → {ALTB_SLOPES[2]:.3f}**, strongly
concave. Pin the kindest single power law to the two extreme required values and
it is ${C_PL:.3f}\\,L_m^{{{N_PL:.3f}}}$, which misses the middle two by
**{ALTB_PL_MISS[1]*100:+.1f} %** and **{ALTB_PL_MISS[2]*100:+.1f} %**. So no
rewriting of Eq. 20's coefficient or exponent reaches these values; it would take
a different correction per curve.

**The shift wins because it needs nothing.** It introduces no free parameter for
three of the four curves: it reuses 0.50, 2.0 and 7.0 exactly as printed, and
only the fourth curve needs a value the paper does not print. The curve positions
return those three to within {max(abs(L_rec[:3]/L_shift[:3]-1))*100:.1f} % in
$L_m$, and the curve spacings — *the same measurement with a common
multiplicative offset in $\\chi$ removed, not a second witness* — come to
{gaps_fig[0]:.3f}/{gaps_fig[1]:.3f}/{gaps_fig[2]:.3f} against the shifted
reading's {g_shifted[0]:.3f}/{g_shifted[1]:.3f}/{g_shifted[2]:.3f}.
"""))'''))

cells.append(md(r"""#### The general argument: $\chi$ separates

Both rivals above are fitted, and a fitted rival can always be answered with
"then fit something else". There is a stronger statement available, and it needs
no fitting at all.

Every variant of Eqs. 20–21 — and every unstated-parameter hypothesis that keeps
the paper's structure — makes $\chi$ **separable**:

$$
\chi(\phi, L_m) = g(L_m)\,h(\phi), \qquad
g = \frac{1 - f_e(L_m)}{\alpha_{gLs}(L_m)}, \qquad h = \eta_s(\phi)\,\phi^2
$$

All of the $L_m$ dependence sits in $g$. A gap between two curves, in decades, is
then $\log_{10} g(L_i)/g(L_j)$ **whatever $h$ is** — so the figure's three
measured gaps, with nothing fitted and no model evaluated, *force* the value of
$d\log g/d\log L_m$ over whichever $L_m$ intervals a reading assigns to them."""))

cells.append(code('''def g_of_L(L):
    """The whole L_m dependence of chi, under the paper's own closures."""
    return (1.0 - f_e(L)) / alpha_gLs(L)

def forced_dlogg(Ls):
    """What the FIGURE's gaps require of d log g / d log L_m, under a reading."""
    Ls = np.asarray(Ls, float)
    return -gaps_fig / np.log10(Ls[1:] / Ls[:-1])

def model_dlogg(Ls):
    Ls = np.asarray(Ls, float)
    return np.diff(np.log10(g_of_L(Ls))) / np.diff(np.log10(Ls))

print("chi = g(L_m) h(phi), so the figure's gaps ALONE fix d log g / d log L_m")
for name, Ls in (("AS PRINTED", [0.5, 1.0, 2.0, 7.0]), ("SHIFTED   ", [0.5, 2.0, 7.0, 10.0])):
    f_, m_ = forced_dlogg(Ls), model_dlogg(Ls)
    ivs = "  ".join(f"[{Ls[i]:g},{Ls[i+1]:g}]".rjust(9) for i in range(3))
    print(f"\\n  {name}  L_m = {', '.join('%g' % x for x in Ls)}")
    print(f"    interval                 {ivs}")
    print("    forced by the figure  " + "".join(f"{v:11.3f}" for v in f_))
    print("    the paper's own g     " + "".join(f"{v:11.3f}" for v in m_))
    print("    relative difference   " + "".join(f"{(a/b-1)*100:+10.1f}%" for a, b in zip(f_, m_)))

FORCED_PRINT, FORCED_SHIFT = forced_dlogg([0.5, 1.0, 2.0, 7.0]), forced_dlogg([0.5, 2.0, 7.0, 10.0])
MODEL_SHIFT = model_dlogg([0.5, 2.0, 7.0, 10.0])
DLOGG_SHIFT_WORST = float(np.max(np.abs(FORCED_SHIFT / MODEL_SHIFT - 1.0)))
print(f"\\n  A power law in L_m would need this slope CONSTANT. Read as printed the figure")
print(f"  demands {FORCED_PRINT[0]:.3f} -> {FORCED_PRINT[1]:.3f} -> {FORCED_PRINT[2]:.3f}:")
print("  not constant, and not even monotone - steeper in the middle than at either end.")'''))

cells.append(code('''display(Markdown(rf"""
Read **as printed**, the figure's own gaps demand
$d\\log g/d\\log L_m$ = **{FORCED_PRINT[0]:.3f}, {FORCED_PRINT[1]:.3f},
{FORCED_PRINT[2]:.3f}** over $[0.5, 1]$, $[1, 2]$ and $[2, 7]$. That is not
constant — so no product of powers of $L_m$ can supply it, which disposes of
every rewriting of Eq. 20 or Eq. 21 in one line. Worse, it is not even
**monotone**: it steepens to {FORCED_PRINT[1]:.3f} in the middle interval and
then flattens to {FORCED_PRINT[2]:.3f}. A power law would give a **constant**;
the paper's own closures give a **monotone steepening**
({model_dlogg([0.5,1.0,2.0,7.0])[0]:.3f} → {model_dlogg([0.5,1.0,2.0,7.0])[1]:.3f}
→ {model_dlogg([0.5,1.0,2.0,7.0])[2]:.3f} over these same three intervals); the
printed reading needs neither.

Read **shifted**, the *same three gaps* demand
**{FORCED_SHIFT[0]:.3f}, {FORCED_SHIFT[1]:.3f}, {FORCED_SHIFT[2]:.3f}** over
$[0.5, 2]$, $[2, 7]$ and $[7, 10]$, against the model's
**{MODEL_SHIFT[0]:.3f}, {MODEL_SHIFT[1]:.3f}, {MODEL_SHIFT[2]:.3f}** — agreement
to within **{DLOGG_SHIFT_WORST*100:.0f} %** on the worst interval, with nothing
fitted.

This is the strongest statement the page can make, because it does not depend on
which rival is proposed. Any hypothesis that leaves $\\chi$ separable — a finite
$\\alpha_{{gs}}$, a corrected $f_e$, different constants in Eq. 21, a different
correlation altogether — is fixed by these three numbers, and the printed reading
asks them for a shape no such $g$ has. **The printed reading is not recoverable
by any reparameterisation; only by a different assignment of $L_m$ to curves.**
(ALT A is the one escape from separability, since a finite $\\alpha_{{gs}}$
couples $\\phi$ and $L_m$ — and it is refuted above on its own terms.)
"""))'''))

cells.append(md(r"""### What none of this can detect

Named as claims the page does **not** make:

- **Eq. 6's superposition is untested.** Replacing a mixed-boundary sphere by a
  weighted average of two whole-surface spheres is the paper's approximation
  (from Ramachandran & Smith 1979), and every route on this page — closed form,
  pymrm solve, reconstruction — is built on top of it. If it is wrong, all of
  them are wrong together. Testing it needs a 2-D axisymmetric solve of one
  sphere with two boundary patches, which is not attempted here.
- **The constants inside Eq. 21 are untested.** They came from the Goto & Smith
  and Specchia correlations evaluated at the caption's conditions. Neither
  correlation is on disk; the page uses Eq. 21 as printed and cannot check the
  substitution.
- **The mechanism of the misalignment is untested.** See above: legend-vs-plot
  cannot be distinguished from the figure.
- **A non-converged solve would be invisible to check 2 at $m = 1$**, because
  the $m=1$ problem is linear and `newton` converges in a single step —
  `maxfev = 1` is not a meaningful defect here. The check's power against
  non-convergence is therefore untested at $m=1$; the printed Newton residual is
  what guards that, and it is reported next to every solve.
- **No claim is made that the model is correct.** It reproduces curves the
  authors computed from it. Figure 6 contains no measurements, and neither does
  this page."""))

cells.append(code('''# the one defect check 2 cannot see at m = 1, measured rather than asserted
s_lin = WettedSphere(n_u=200)
_, r_full = s_lin.solve(10.0, float(alpha_gLs(2.0)), maxfev=60)
_, r_one = s_lin.solve(10.0, float(alpha_gLs(2.0)), maxfev=1)
e_full = s_lin.eta(s_lin.solve(10.0, float(alpha_gLs(2.0)), maxfev=60)[0])
e_one = s_lin.eta(s_lin.solve(10.0, float(alpha_gLs(2.0)), maxfev=1)[0])
print(f"m = 1 is linear: maxfev=60 -> eta {e_full:.10f} (residual {r_full:.1e})")
print(f"                 maxfev=1  -> eta {e_one:.10f} (residual {r_one:.1e})")
print(f"                 difference {abs(e_one/e_full-1):.2e}  <- so 'maxfev=1' is not a defect at m=1")'''))

# ------------------------------------------------------------- what pymrm adds
cells.append(md(r"""## What pymrm adds

**Honestly: for Figure 6 itself, nothing.** The figure is reproduced by
evaluating a closed form the paper prints. The pymrm solve returns the same
numbers to the discretisation error measured above, and it exists on this page
as an *independent route* — the thing that makes "we transcribed Table 1
correctly" a claim that can fail rather than an assertion.

Three things it does add, all of them beyond where Table 1 stops.

**1. The closed form is first-order only; the pellet is not.** Table 1 exists
because $\nabla^2 C = \phi^2 C$ is linear. Real hydrogenations are not first
order in the volatile reactant, and there the partial-wetting gain has no closed
form at all. The finite-volume route does not care."""))

cells.append(code('''print("chi at L_m = 2 for reaction orders the closed form cannot reach")
print("  phi      m=1        m=1.5      m=2       (m=2)/(m=1)")
rows_m = []
for p in (2.0, 5.0, 10.0, 20.0, PHI_HI):
    r = [sph.chi(p, 2.0, m=mm) for mm in (1.0, 1.5, 2.0)]
    rows_m.append(r)
    print(f"  {p:<8.4g} {r[0]:<10.4f} {r[1]:<10.4f} {r[2]:<10.4f} {r[2]/r[0]:.3f}")
M_RATIO_HI = rows_m[-1][2] / rows_m[-1][0]
M_RATIO_LO = rows_m[0][2] / rows_m[0][0]'''))

cells.append(code('''display(Markdown(rf"""
The correction **changes sign across the range**: at $\\phi = 2$ a second-order
reaction gains {(M_RATIO_LO-1)*100:.0f} % *more* from partial wetting than a
first-order one, while at $\\phi = {PHI_HI}$ it gains
{abs(M_RATIO_HI-1)*100:.0f} % *less*. So reading Figure 6 for a second-order
hydrogenation at high $\\phi$ **over**states the partial-wetting gain by about
{abs(M_RATIO_HI-1)*100:.0f} %, and at low $\\phi$ understates it by about
{(M_RATIO_LO-1)*100:.0f} %. Table 1 cannot say this — the closed form does not
exist for $m \\neq 1$, and the pymrm route needs no change beyond an argument.
"""))'''))

cells.append(md(r"""**2. The assumptions become sweeps.** $\alpha_{gs}$, $C^*_L$ and the reaction
order are arguments, not hard-wired limits, so the audit in Validation is three
lines rather than a re-derivation.

**3. The pellet is the same object as `B1.1`.** Same interior equation, same
operators, different surface condition — which is the gallery's structural
claim for `S8` made concrete: `B1.1`'s Dirichlet pellet is the
$\alpha \to \infty$ member of this family, and the page checks that limit
explicitly above.

What pymrm does **not** add here: it does not lift Eq. 6's superposition, which
would need a 2-D axisymmetric sphere with two boundary patches. That is a real
extension and it is not on this page."""))

# ------------------------------------------------------------------ agreement
cells.append(md(r"""## Agreement"""))

cells.append(code('''metrics = {
    # the finding
    "L_rec_curve1": L_rec[0], "L_rec_curve2": L_rec[1],
    "L_rec_curve3": L_rec[2], "L_rec_curve4": L_rec[3],
    # chi at phi = 10, both readings, dev = (model - figure)/figure
    "chi10_dev_shifted_worst": float(np.max(np.abs(d10_shift))),
    "chi10_dev_printed_worst": float(np.max(np.abs(d_printed))),
    "chi_dev_shifted_worst_2to30": SHIFT_WORST,
    "chi_dev_shifted_worst_2to30_printed_Lm_only": SHIFT_WORST_123,
    "chi_dev_printed_worst_2to30": PRINT_WORST,
    # the unexplained residual - size quoted as a range, dominated by the abscissa
    "slope_model_fitted": model_slope,
    "slope_figure_mean": float(SLOPE.mean()),
    "slope_gap": SLOPE_GAP,
    "slope_gap_lo_over_abscissa_calibrations": SLOPE_GAP_LO,
    "slope_gap_hi_over_abscissa_calibrations": SLOPE_GAP_HI,
    "model_straightline_fit_rms_decades": model_fit_rms,
    # the separability argument
    "dlogg_forced_shifted_worst_rel_err": DLOGG_SHIFT_WORST,
    # declared limits of check 4
    "recon_check_best_unrejected_defect": CHECK4_BEST_DEFECT,
    "eq21_first_term_resistance_share_at_Lm1": SHARE_FIRST,
    # ALT B, pressed on curvature rather than monotonicity
    "altb_best_powerlaw_worst_miss": ALTB_PL_WORST,
    # the gap statistic's own phi-drift
    "figure_gap_drift_over_window_decades": GAP_DRIFT,
    # numerics
    "pymrm_vs_table1_maxdev": ETA_MAXDEV,
    "grid_order": GRID_ORDER,
    "curve4_dev_worst_at_reconstructed_Lm": C4_EXACT_WORST,
    "eq20_max_abs_err_vs_printed": BASE["eq20"],
    "collapse_identity": IDENTITY,
}
report_agreement("G1.8", metrics)'''))

cells.append(code('''display(Markdown(rf"""
**Summary.** Reading Figure 6's legend as printed, the paper's own chain misses
three of its four curves by up to
{np.max(np.abs(d_printed))*100:.0f} %. Reading the legend as offset by one row
from its second entry down — so that the drawn curves are
$L_m$ = 0.50, 2.0, 7.0 and ≈10 — brings all four to within
{np.max(np.abs(d10_shift))*100:.1f} % at $\\phi = 10$ and
{SHIFT_WORST*100:.0f} % over $\\phi$ = 2–30, with **nothing fitted** for the
three curves that reuse printed $L_m$ values.

Two further readings of the same figure support it, and neither is an
independent witness. The curve **spacings** — the same four positions with a
common multiplicative offset in $\\chi$ removed — are matched by the shifted
reading to {GAP_MISS_SHIFT:.3f} decades, while the printed reading misses them by
up to {GAP_MISS_PRINT:.3f}. And because $\\chi$ **separates** into $g(L_m)h(\\phi)$,
those same three gaps force $d\\log g/d\\log L_m$ without any model at all: read
as printed they demand {FORCED_PRINT[0]:.3f}/{FORCED_PRINT[1]:.3f}/{FORCED_PRINT[2]:.3f},
which is non-monotone and so unreachable by any product of powers of $L_m$; read
shifted they demand {FORCED_SHIFT[0]:.3f}/{FORCED_SHIFT[1]:.3f}/{FORCED_SHIFT[2]:.3f}
against the model's {MODEL_SHIFT[0]:.3f}/{MODEL_SHIFT[1]:.3f}/{MODEL_SHIFT[2]:.3f}.
No reparameterisation rescues the printed reading; only a reassignment does.

One residual survives and is reported as unexplained: the drawn lines are
steeper than the model in log–log slope ({SLOPE.mean():.3f} against
{model_slope:.3f} as digitised). Its **sign** is robust, its **size** is not —
the abscissa calibration, not the fit rms, dominates, and across defensible
calibrations the gap spans {SLOPE_GAP_LO:+.3f} to {SLOPE_GAP_HI:+.3f}. $L_m$
cannot produce it either way, because the chain factorises.

The pymrm finite-volume solve reproduces Table 1's sphere row to
{ETA_MAXDEV:.1e} at $n_u$ = {sph.n_u}, at an observed order between
{GRID_ORDER_MIN:.1f} and {GRID_ORDER:.1f}. Deviation is
(model − figure)/figure throughout. Tier 6: the "figure" is the authors' own
computed output, not a measurement.
"""))'''))

# ----------------------------------------------------------------------- reuse
cells.append(md(r"""## Reuse

| you want | change |
|---|---|
| a different pellet geometry | `WettedSphere(nu=0)` slab, `nu=1` cylinder — Table 1 prints those rows too, and the $\eta$ term becomes $\phi/\tanh\phi$ or $\phi I_0(2\phi)/I_1(2\phi)$ |
| a different reaction order | `sph.chi(phi, L_m, m=...)` — no closed form needed |
| a finite gas-side resistance | `sph.chi(..., alpha_gs=...)`; see the sweep in Validation |
| a volatile reactant not at saturation | `C_star_L` on `solve`/`chi` |
| a different wetting correlation | replace `f_e`; Eq. 20 is a one-liner |
| different mass-transfer correlations | replace `alpha_gLs`; Eq. 21 already bakes in the caption's conditions, so a different system needs Eq. 4d re-evaluated |
| the fully wetted pellet | `alpha_gs = alpha_gLs` and `f_e = 1` — that is `B1.1` |

**Related pages.** `G1.7` (the hydrodynamics of the same reactor, including the
hold-up that sets $f_e$), `B1.1` (the Dirichlet pellet this generalises),
`B1.4` (the observable form of the same Thiele modulus), `F3.1` (the
gas–liquid film that $\alpha_{gLs}$ lumps).

## References

Herskowitz, M. and Smith, J. M. (1983). Trickle-bed reactors: a review.
*AIChE Journal* **29**(1) 1–18. doi:10.1002/aic.690290102 — the source, read
from a 600 dpi render of a lawfully accessed copy. Table 1 (page 4),
Eqs. 19–21 and Figure 6 (page 8).

Cited by that paper and used only through it, not consulted here:
Ramachandran, P. A. and Smith, J. M. (1979) for the Eq. 6 superposition;
Goto, S. and Smith, J. M. (1975) and Specchia, V. et al. (1978) for the
mass-transfer correlations condensed into Eq. 21; Mills, P. L. and Dudukovic,
M. P. (1979) and Herskowitz, M. (1981) for the wetting-efficiency work the
section reviews.

**Two printed defects reported by this page**, both checked against 600 dpi
renders: the page-8 reference to "Table 2" for the spherical solution, which is
in Table 1; and the misalignment of Figure 6's legend block against its curves.
Neither affects the model — only a reader trying to use the figure."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print(f"wrote index.ipynb with {len(cells)} cells")
