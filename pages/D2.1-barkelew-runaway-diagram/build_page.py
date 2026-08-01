#!/usr/bin/env python3
"""Generate index.ipynb for page D2.1. Run from the page directory."""
from __future__ import annotations

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "Barkelew's runaway diagram"
description: "The 1959 chart that collapses the runaway boundary of a cooled tubular reactor onto two dimensionless groups - rebuilt from the paper that reprints it, and tested where a similarity reduction is supposed to fail."
categories: [sec:D, struct:S2, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-01
---

# Barkelew's runaway diagram

**Catalog ID:** `D2.1` · **Structures:** `S2` (plug flow with reaction) · **Tier:** T0

A cooled tubular reactor runs away when the heat it makes outruns the heat the
wall can take. Barkelew integrated that problem a great many times and found
that the boundary between safe and unsafe operation collapses onto **two
dimensionless groups**, so that one chart covers every first-order exothermic
reaction. This page rebuilds the chart, and then asks the question a collapse
always invites: how well does it actually collapse?"""))

cells.append(md(r"""## Background

[`D2.2`](../D2.2-van-welsenaere-froment-runaway/) is the other half of this
story. Van Welsenaere and Froment (1970) located the runaway boundary from
**intrinsic geometric features** of the trajectory in the $p$–$T$ plane, and
turned them into closed-form inequalities on the inlet conditions — explicit
formulae in dimensional variables, one reactor at a time.

Barkelew (1959) did the opposite. He integrated the reactor equations over a
large parameter range and *correlated* the results, choosing groups in which
the whole family of solutions falls onto one picture. There is no formula at
the end of it; there is a chart, and using it takes trial and error. What it
buys is generality: the same four curves serve any tube, any coolant
temperature, any first-order exothermic chemistry.

That difference — an explicit criterion against a similarity reduction — is why
this is a separate page and not a section of `D2.2`. A collapse is a claim that
a two-parameter family really is two-parameter, and claims like that are
testable, including where they fail. `D2.2` has no such claim to test.

**Terminology, as both papers are careful about it.** With gas and solid at the
same temperature an ideal tubular reactor is always stable in the strict sense;
there is no bifurcation and no multiplicity here, unlike
[`B1.1`](../B1.1-thiele-weisz-hicks/). "Runaway" means *parametric sensitivity*:
a region where a one percent change in an inlet variable moves the hot spot by
a hundred degrees."""))

cells.append(md(r"""## The published model

### Where this comes from, and where it does not

**Barkelew's paper was not consulted.** *Chem. Engng Prog. Symp. Ser.* 55 (1959)
37 is pre-DOI, is not on disk, and has no open-access route. Everything on this
page is read from **Van Welsenaere & Froment (1970)**, whose Section 5 is
explicitly a comparison with Barkelew: it reprints his diagram as their Fig. 9,
says what his axes mean, states the one numerical property of his critical
locus, and applies their own criterion in his notation. That is the reprint
route this repository uses for `F1.3` and `B1.6`, and the distinction is kept in
the metadata: `reference` is Barkelew, `reference_read_from` is Van Welsenaere &
Froment.

The consequence matters for reading this page. **Everything below is Barkelew's
method as Van Welsenaere and Froment describe and use it.** Where his own paper
would have given definitions, theirs gives sentences, and the reconstruction of
his groups from those sentences is set out explicitly so that it can be checked
rather than believed.

### The reactor

The starting point is the same model as `D2.2` — their Eqs. 3–4. One dimension,
pseudo-homogeneous, constant wall temperature, one irreversible pseudo-first-order
reaction:

$$
\frac{\mathrm{d}p}{\mathrm{d}z} = -A\,p\,\mathrm{e}^{-a/T+b},
\qquad
\frac{\mathrm{d}T}{\mathrm{d}z} = B\,p\,\mathrm{e}^{-a/T+b} - C\,(T-T_w),
$$

$$
A = \frac{M P \rho_b}{\rho_g} p_B^0, \qquad
B = \frac{(-\Delta H)\rho_b}{c_p} p_B^0, \qquad
C = \frac{2U}{c_p R}, \qquad z = \frac{z'}{u},
$$

with $p = p^0$, $T = T_0 = T_w$ at $z = 0$. The axial coordinate is a contact
time, not a length.

### Barkelew's simplification

Van Welsenaere and Froment write that Barkelew *"modified the temperature
dependence of the rate coefficient in such way that it may be characterized by
only one constant instead of two"*. The two constants are $a$ and $b$; the
modification that removes one of them is the exponential approximation to the
Arrhenius factor about the wall temperature,

$$
\mathrm{e}^{-a/T+b} \;\longrightarrow\; k_w\,\mathrm{e}^{\tau},
\qquad
k_w = \mathrm{e}^{-a/T_w+b},
\qquad
\boxed{\;\tau = \frac{a\,(T-T_w)}{T_w^{2}}\;}
$$

which is what makes $\tau$ the dimensionless temperature on the axes of his
diagram."""))

cells.append(md(r"""### What is forced by the printed text, and what is not

Section 5 gives $\tau$ a name and uses it, but never writes it down, so the
boxed definition is a **reconstruction**. Most of it is forced by what is
printed. One step is not, and the page says which is which rather than calling
the whole thing pinned.

**Forced — the shape of the reduced system, and the reference temperature.**
The two sentences Section 5 gives in words are

> *"$N/S$ is the ratio of the rate of heat transfer per unit volume at $\tau = 1$
> to the rate of heat generation per unit volume at $\tau = 0$. The parameter $S$
> corresponds to the value $\tau$ would reach if the reactants were fed at the
> wall temperature and if the reactor were operated adiabatically."*

The first fixes the cooling term as $N\tau$ and the generation as
$S\,y\,\mathrm{e}^{\tau}$ — at $\tau = 1$ the transfer rate is $N$, at $\tau=0$
with $y=1$ the generation rate is $S$ — and a cooling term proportional to
$\tau$ forces $\tau \propto (T - T_w)$ with $\tau(T_w) = 0$. The second, and
only the second, names the wall temperature explicitly: it feeds the reactants
at it, so with $N = 0$ and $\mathrm{d}\tau/\mathrm{d}y = -S$ the adiabatic
$\tau$ is $S$ itself, and their own Eq. 32 turns that into
$S = a(\Delta T)_{ad}/T_w^2$.

**Forced too, by a third printed statement, is that the exponent is $\tau$ and
not a multiple of it.** Section 5 states that their own first criterion becomes
$\tau_m = 1$ in these variables. With a rate $\propto \mathrm{e}^{\lambda\tau}$
the maxima curve is $y \propto \tau\,\mathrm{e}^{-\lambda\tau}$ and its peak
sits at $\tau = 1/\lambda$, so their printed $\tau_m = 1$ forces $\lambda = 1$.

**Not forced — the temperature at which the exponent is linearised.** Nothing in
the paper says it. $a(T-T_w)/T_w^2$ is the tangent to $\ln k$ at the wall, but
the tangent at their own critical hot spot, $a(T-T_w)/T_M^2$, and the chord
between the two, $a(T-T_w)/(T_w T_M)$, satisfy every sentence quoted above
equally well: all three vanish at $T_w$, all three are linear in $T - T_w$, and
each comes with its own consistent $k_w$. That is a Frank–Kamenetskii
convention, not a constraint from the paper.

It is, however, **testable**, and this page tests it twice — once as a
deliberate defect in the break table, and once against the blow-up positions of
their Fig. 10. The wall tangent wins both times.

Substituting and rescaling the contact time by $\zeta = A k_w z$ leaves **two
parameters and nothing else**:

$$
\frac{\mathrm{d}y}{\mathrm{d}\zeta} = -y\,\mathrm{e}^{\tau},
\qquad
\frac{\mathrm{d}\tau}{\mathrm{d}\zeta} = S\,y\,\mathrm{e}^{\tau} - N\,\tau,
\qquad
y(0)=1,\; \tau(0)=0,
$$

with $y = p/p^0$ and

$$
S = \frac{a\,(\Delta T)_{ad}}{T_w^{2}} = \frac{a\,B\,p^0}{A\,T_w^{2}},
\qquad
N = \frac{C}{A\,k_w}.
$$

$N/S = C T_w^2 / (a B p^0 k_w)$ is exactly the ratio of the two rates their first
sentence describes. Note what has happened to the inlet pressure: $S \propto p^0$
and $N$ does not contain $p^0$ at all.

### The criterion

Barkelew plots $\tau_m/S$ against $N/S$ at constant $S$, where $\tau_m$ is the
hot spot. Section 5:

> *"The curves of constant $S$ have an envelope. Above the tangent to the
> envelope $\tau_m$ changes rapidly with $N/S$. Below that point it does not.
> The points of tangency therefore correspond to critical conditions."*

So the criterion is geometric: on each $S$-curve, the point where it touches the
envelope of the whole family. And then the one number Section 5 prints:

> *"It was observed that for all $S$-values $\tau_m$ at the point of tangency was
> very close to $1{\cdot}275$."*

That is the claim this page is built to test, and it is a **similarity claim** —
a statement that one number characterises the critical condition for every $S$.

### Where the two papers meet

Section 5 also carries Van Welsenaere and Froment's own criterion into
Barkelew's variables. Their first criterion puts the critical trajectory through
the maximum of the maxima curve; with the exponential rate law the maxima curve
is $y_m \propto \tau\,\mathrm{e}^{-\tau}$, whose peak is at $\tau = 1$ for any
choice of reference temperature. Hence:

> *"With Barkelew's dimensionless notation the first criterion leads to
> $\tau_m = 1$, so that $\tau_m/S = 1/S$. As indicated in Fig. 9 the corresponding
> point on each $S$-curve lies a bit to the right of the points of tangency
> between the $S$-curve and the envelope. Consequently, the first criterion
> proposed in this paper is slightly more conservative than Barkelew's."*

Two testable statements: the ordering (right of, for every $S$), and — through
Fig. 8, where they mark Barkelew's predictions against their own — that
*"Barkelew's predictions agree extremely well with these based upon the present
methods"*. Their Fig. 8 shows that comparison as marker positions only, so this
page does not read it; it recomputes **both** sides."""))

cells.append(md(r"""## Parameters and assumptions

The dimensionless system needs no parameters at all — that is the point of it.
Parameters enter only where the page leaves the diagram and asks what
Barkelew's criterion predicts for a specific tube, which is where it can be
compared with `D2.2`. For that the base case of their Section 1 is used
unchanged.

**Assumptions carried, all of them the papers'.** Plug flow; no axial
dispersion of heat or mass; constant wall temperature; constant physical
properties; one irreversible reaction, first order in the limiting reactant with
the second reactant in large excess; gas and solid at the same temperature. On
top of these, Barkelew's own: the Arrhenius factor replaced by an exponential in
$\tau$, which is accurate near $T_w$ and increasingly wrong above it — and
wrong in the unsafe direction, as the results below show.

**One trap in getting the numbers out**, inherited from `D2.2`. This is a 1970
scan and the Elsevier full-text API returns the publisher's OCR of it, which
discards the mid-dot decimal separator: `1275` for 1.275, `001651` for 0.01651,
`19837` for 19.837. The API is excellent for prose — it is how Section 5 was
located and read — and useless for numbers. Every number here was read from a
600 dpi render of the printed page."""))

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

import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.sparse import SparseEfficiencyWarning
from scipy.sparse.linalg import MatrixRankWarning
from pymrm import construct_convflux_upwind, construct_div, NumJac, newton
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "D2.1-barkelew-runaway-diagram"
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

cells.append(md(r"""## The data

Two files. The **parameters** are the base case of Section 1, needed only for
the dimensional comparison. The **Barkelew file** is everything Section 5 prints
about the diagram, plus the Section 6 values used to prove the parameter reading
is right.

**These are not measurements.** Barkelew's chart correlates his own numerical
integrations; Van Welsenaere and Froment's numbers are their own closed-form
extrapolations and Runge–Kutta runs. Neither paper contains an experiment, so
this page is provenance tier 6 — validated against a published reference
solution, in the same sense as [`B1.1`](../B1.1-thiele-weisz-hicks/) and
[`D2.2`](../D2.2-van-welsenaere-froment-runaway/), and not against experiment.

**No curve was digitised into these files.** The `label` and `axis` rows below
are text printed inside the figure frames. The markers on Fig. 9 — the crosses
at the tangency points and the circles at $\tau_m/S = 1/S$ — were deliberately
not measured: that would be a figure digitisation, and the printed 1.275 makes
it unnecessary. Four positions *were* measured off Fig. 10, for the one question
nothing printed can settle; they are quarantined in their own cell near the end
of the page, they are not in either CSV, and they are awaiting a maintainer's
visual review."""))

cells.append(code('''par = load_data("van-welsenaere-froment-1970-parameters.csv", page=PAGE)
bk = load_data("van-welsenaere-froment-1970-barkelew.csv", page=PAGE)
bk_meta = load_meta("van-welsenaere-froment-1970-barkelew.csv", page=PAGE)

P = dict(zip(par.symbol, par.value))
V = dict(zip(bk.quantity, bk.value))

print(par.to_string(index=False))
print()
print(bk[["quantity", "value", "unit", "kind", "source"]].to_string(index=False))
print(f"\\n{cite_data(bk_meta)}")
print("origin of the result (not consulted): "
      + " ".join(f"{bk_meta['reprint_of']['authors'][0]} "
                 f"({bk_meta['reprint_of']['year']}), "
                 f"{bk_meta['reprint_of']['container']}".split()))
print(f"provenance tier {bk_meta['provenance_tier']['tier']}")'''))

cells.append(md(r"""## PyMRM implementation

Three pieces, and they are deliberately kept separate because they fail
independently.

1. **The dimensionless Barkelew tube**, solved with pymrm. Convection in $\zeta$
   with a pointwise source — `construct_convflux_upwind` then `construct_div`,
   the same two calls as `D2.2` and [`C2.1`](../C2.1-xu-froment-smr/), with the
   state $(y, \tau)$. This produces the diagram.
2. **The same system in the phase plane**, $\mathrm{d}\tau/\mathrm{d}y$, solved
   by adaptive quadrature. It never forms a grid, and its right-hand side is a
   **separate transcription** of the reduced model — not a call into the
   finite-volume one — so it is the independent route against which the pymrm
   discretisation is measured. This is the pattern `D2.2` used for its critical
   inlet pressure.
3. **The dimensional tube**, with the true Arrhenius factor or Barkelew's
   exponential, for the comparison of their Fig. 3 with their Fig. 10.

Four implementation choices are worth stating.

*The domain is chosen, not fixed.* The hot spot sits at $\zeta \approx 0.1$ for
$S = 32$ and $\zeta \approx 0.5$ for $S = 4$, so one fixed domain either misses
the peak or wastes its cells on the tail. A cheap 200-cell pass locates the
peak, the fine solve then covers three times that. This is arithmetic on the
grid, not a continuation: nothing is warm-started and the answer does not depend
on the order the cases are solved in.

*The integration is capped at $\tau = 20$.* Deep inside the runaway region
$\tau_m \to S$ and $\mathrm{e}^{\tau}$ overflows for the larger $S$. Every
$\tau$ this page reports is below 13, so the cap only makes the useless region
cheap; where it binds, the value returned is the cap and is never used as a
result.

*Root-finds are bracketed and the bracket is asserted.* The tangency search is
bracketed by the two hot spots $\tau_m = 1.7$ and $\tau_m = 0.7$, chosen to
straddle it without reference to the paper's 1.275, and `brentq` is given a sign
change it has verified. No continuation chain enters any reported number.

*The tangency itself has a step size, and it is converged rather than chosen.*
The envelope is located from a difference of $\tau_m/S$ in $\ln S$ over a step
$h$; that difference has an $O(h^2)$ truncation error, so the raw value at any
fixed $h$ is biased. The pair $(h, h/2)$ is Richardson-extrapolated, and the
next section measures what that is worth. The extrapolated value is a pure
function of $S$, so it is memoised — arithmetic reuse, not a warm start."""))

cells.append(code('''A = P["M"] * P["P"] * P["rho_b"] / P["rho_g"] * P["p_B0"]
B = P["minus_dH"] * P["rho_b"] / P["c_p"] * P["p_B0"]
a_exp, b_exp = P["a"], P["b"]


def C_of_R(R):
    # 2U/(c_p R): c_p is printed VOLUMETRIC, so no rho_g here.
    return 2.0 * P["U"] / (P["c_p"] * R)


C_BASE = C_of_R(P["R"])


def kexp(T):
    """The true Arrhenius factor of their Eqs. 3-4."""
    return np.exp(b_exp - a_exp / np.asarray(T, float))


def T_M(Tw):
    """Their Eq. 8: the critical hot spot of the first criterion."""
    return 0.5 * (a_exp - np.sqrt(a_exp * (a_exp - 4.0 * Tw)))


# --- Barkelew's two groups, from the sentences quoted above ------------------
def N_group(Tw, R=None):
    """C/(A k_w): cooling capacity over reaction rate. Independent of p0."""
    return (C_BASE if R is None else C_of_R(R)) / (A * kexp(Tw))


def S_group(p0, Tw):
    """a (dT)_ad / Tw^2, with (dT)_ad = (B/A) p0 from their Eq. 32."""
    return a_exp * B * p0 / (A * Tw ** 2)


def p0_of_S(S, Tw):
    return S * A * Tw ** 2 / (a_exp * B)


def check5(cp_fac=1.0):
    """Check 5, as a number: worst |dev| over four printed derived quantities.

    `cp_fac` exists for the break test - multiplying the printed volumetric c_p
    by rho_g is the trap this page inherits from D2.2, and c_p enters only
    through C, hence only through ln K.
    """
    Cc = 2.0 * P["U"] / (P["c_p"] * cp_fac * P["R"])
    rows = [
        ("t_w  (Eq. 21b)", a_exp / 625.0, V["t_w_at_Tw625"]),
        ("ln K (Eq. 22)", float(np.log(A / Cc * np.exp(b_exp - 20.0))), V["lnK_at_Tw625"]),
        ("T_M  (Eq. 8), K", float(T_M(625.0)), V["T_M_at_Tw625"]),
        ("(dT)_ad (Eq. 32) at p0=0.0125 atm, K", B / A * 0.0125, V["dT_ad_at_p0_0125"]),
    ]
    df = pd.DataFrame(rows, columns=["quantity", "recomputed", "printed"])
    df["dev %"] = (df.recomputed - df.printed).abs() / df.printed.abs() * 100
    return df


print("The parameter reading, against four numbers the paper prints and that")
print("were not used to obtain it (a single lost decimal breaks all four):")
chk = check5()
print(chk.to_string(index=False, float_format=lambda v: f"{v:.4f}"))'''))

cells.append(code('''TAU_CAP = 20.0


def tau_max_phase(N, S, cap=TAU_CAP, rtol=1e-9, expo=np.exp, cool=1.0, order=1.0):
    """Hot spot from the phase-plane form dtau/dy = -S + N tau /(y^n e^tau).

    Independent variable is the reactant fraction y, so this route never forms
    a zeta grid at all. This `rhs` is a SEPARATE transcription of the reduced
    model from `BarkelewTube.reaction` below; the two share no code, which is
    what gives the cross-check of them power. `expo`, `cool` and `order` exist
    only for the break tests; the model is order = 1, cool = +1, expo = exp.
    """
    def rhs(y, t):
        return [-S + cool * N * t[0] / (y ** order * expo(t[0]))]

    def flat(y, t):                      # dtau/dy = 0: the hot spot
        return S * y ** order * expo(t[0]) - cool * N * t[0]

    def hot(y, t):                       # useless-region cap
        return t[0] - cap

    flat.direction, flat.terminal = -1.0, True
    hot.direction, hot.terminal = 1.0, True
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        sol = solve_ivp(rhs, [1.0, 1e-12], [0.0], rtol=rtol, atol=1e-12,
                        events=[flat, hot], method="LSODA")
    if sol.t_events[0].size:
        return float(sol.y_events[0][0][0])
    if sol.t_events[1].size:
        return cap
    return float(sol.y[0].max())


class BarkelewTube:
    """The dimensionless tube, state (n_zeta, 2) = [y, tau].

    `cool` and `s_fac` exist only for the break tests: they inject a defect
    into THIS route's reaction term and nowhere else, which is how the page
    measures whether the cross-check against the quadrature can see the model.
    """

    def __init__(self, N, S, zeta_end, n_z, nu=0, cool=1.0, s_fac=1.0):
        self.N, self.S = float(N), float(S)
        self.cool, self.s_fac = float(cool), float(s_fac)
        self.z_f = np.linspace(0.0, zeta_end, n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.shape = (n_z, 2)
        # outward normal, so both dicts read a.dx/dn + b.x = d
        # inlet : y = 1, tau = 0        -> a=0, b=1, d=(1, 0)
        # outlet: dx/dn = 0, outflow    -> a=1, b=0, d=0
        self.bc = ({"a": 0.0, "b": 1.0, "d": np.array([1.0, 0.0])},
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        self.u = np.tile(np.array([1.0, 0.0]), (n_z, 1))
        self.diverged = False
        # v = 1 by construction: zeta is a scaled contact time, so
        # d(v x)/dzeta = source.  nu = 0: contact time is Cartesian.
        cm, cbc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                            self.bc, v=1.0, axis=0)
        dm = construct_div(self.shape, self.z_f, nu=nu, axis=0)
        self.jac_const = dm @ cm
        self.g_const = dm @ cbc
        self.numjac = NumJac(self.shape)   # pointwise source: last axis only

    def reaction(self, u):
        y, t = u[..., 0], u[..., 1]
        r = y * np.exp(np.clip(t, -50.0, 50.0))
        return np.stack([-r, self.s_fac * self.S * r - self.cool * self.N * t],
                        axis=-1)

    def residual(self, u):
        g_r, j_r = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_r.reshape((-1, 1))
        return g, self.jac_const - j_r

    def solve(self, maxfev=200):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", MatrixRankWarning)
            warnings.simplefilter("ignore", SparseEfficiencyWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            try:
                res = newton(self.residual, self.u, maxfev=maxfev)
            except (ValueError, RuntimeError):
                self.diverged = True
                return None
        self.u = res.x.reshape(self.shape)
        self.diverged = not np.all(np.isfinite(self.u))
        return res


def tau_max_pymrm(N, S, n_z=1600, nu=0, cool=1.0, s_fac=1.0, return_profile=False):
    """Two passes: 200 cells to find the peak, then n_z cells around it."""
    ze = 1.0
    for _ in range(8):
        coarse = BarkelewTube(N, S, ze, 200, nu=nu, cool=cool, s_fac=s_fac)
        coarse.solve()
        if coarse.diverged:
            return (np.inf, None) if return_profile else np.inf
        i = int(np.argmax(coarse.u[:, 1]))
        if i < 180:
            break
        ze *= 2.0
    fine = BarkelewTube(N, S, max(3.0 * coarse.z_c[i], 1e-3), n_z, nu=nu,
                        cool=cool, s_fac=s_fac)
    fine.solve()
    if fine.diverged:
        return (np.inf, None) if return_profile else np.inf
    tm = float(fine.u[:, 1].max())
    return (tm, fine) if return_profile else tm'''))

cells.append(code('''def u_at_tau(S, target, tm=tau_max_phase, lo=1e-4, hi=40.0):
    """N/S at which the hot spot equals `target`; tau_m falls monotonically."""
    return brentq(lambda u: tm(u * S, S) - target, lo, hi, xtol=1e-10)


def tangency_at_h(S, tm=tau_max_phase, h=0.04):
    """Barkelew's critical point at ONE step size h, unextrapolated.

    The envelope of a one-parameter family is the locus where d/dS vanishes at
    fixed abscissa, so the tangency on curve S is the N/S at which the curves
    for S e^{+h} and S e^{-h} cross. That two-point difference is a centred
    difference in ln S, so it carries an O(h^2) truncation error. Bracketed by
    the hot spots tau_m = 1.7 and 0.7, which are neutral round numbers
    straddling the answer; the sign change is asserted, not assumed.
    """
    Sp, Sm = S * np.exp(h), S * np.exp(-h)
    f = lambda u: tm(u * Sp, Sp) / Sp - tm(u * Sm, Sm) / Sm
    lo, hi = u_at_tau(S, 1.7, tm), u_at_tau(S, 0.7, tm)
    if f(lo) * f(hi) >= 0:
        return None, None                     # no envelope: the check collapses
    ur = brentq(f, lo, hi, xtol=1e-8)
    return ur, tm(ur * S, S)


H_TAN = 0.01          # coarse member of the Richardson pair (h, h/2)


def _tangency_rich(S, tm, h):
    """Richardson extrapolation of the O(h^2) two-point tangency to h -> 0."""
    u1, t1 = tangency_at_h(S, tm, h)
    if u1 is None:
        return None, None
    u2, t2 = tangency_at_h(S, tm, 0.5 * h)
    if u2 is None:
        return None, None
    return (4.0 * u2 - u1) / 3.0, (4.0 * t2 - t1) / 3.0


@lru_cache(maxsize=None)
def _tangency_cached(S, h):
    return _tangency_rich(S, tau_max_phase, h)


def tangency(S, tm=None, h=H_TAN):
    """The converged tangency: abscissa N/S and ordinate tau_m.

    `tm=None` means the page's own model, whose value is a pure function of
    (S, h) and is therefore memoised. Passing a `tm` is how the break tests
    substitute a defective right-hand side; those are not cached.
    """
    if tm is None:
        return _tangency_cached(float(S), float(h))
    return _tangency_rich(S, tm, h)


S_PRINTED = [V["S_curve_1"], V["S_curve_2"], V["S_curve_3"], V["S_curve_4"]]
S_WIDE = [3., 4., 6., 8., 12., 16., 24., 32., 48., 64., 100., 150., 200.]
print("the four curves Barkelew's diagram is drawn with:", S_PRINTED)'''))

cells.append(md(r"""## Results

### First: the step size the whole page rests on

Every number below is a tangency, and a tangency is located from a finite
difference. Before quoting any of them to four figures, here is what the step
size does. The table gives $\tau_m$ at the tangency at five step sizes and the
two Richardson extrapolations that can be formed from the last three."""))

cells.append(code('''import time
t_start = time.time()

H_STUDY = [0.04, 0.02, 0.01, 0.005, 0.0025]
S_STUDY = S_PRINTED + [200.0]
raw = {S: {h: tangency_at_h(S, h=h)[1] for h in H_STUDY} for S in S_STUDY}

step = pd.DataFrame({"S": S_STUDY})
for h in H_STUDY:
    step[f"h={h:g}"] = [raw[S][h] for S in S_STUDY]
step["Richardson(.01,.005)"] = [(4 * raw[S][0.005] - raw[S][0.01]) / 3 for S in S_STUDY]
step["Richardson(.005,.0025)"] = [(4 * raw[S][0.0025] - raw[S][0.005]) / 3 for S in S_STUDY]
conv = {S: step["Richardson(.005,.0025)"][i] for i, S in enumerate(S_STUDY)}
step["error at h=0.04, %"] = [(raw[S][0.04] - conv[S]) / conv[S] * 100
                              for S in S_STUDY]
print("tau_m at the envelope tangency, against the step size h in ln S")
print(step.to_string(index=False, float_format=lambda v: f"{v:9.6f}"))
resid = max(abs(step["Richardson(.005,.0025)"][i] - step["Richardson(.01,.005)"][i])
            / conv[S] * 100 for i, S in enumerate(S_STUDY))
err = np.array([[abs(raw[S][h] - conv[S]) / conv[S] * 100 for h in H_STUDY]
                for S in S_STUDY])
order = np.polyfit(np.log(H_STUDY), np.log(err.mean(axis=0)), 1)[0]
print(f"\\nobserved order of the truncation error: {order:.2f}  (O(h^2) expected)")
print(f"the two Richardson estimates differ by at most {resid:.2e} % of the value,")
print(f"so the extrapolated tangency is converged well past the four figures quoted.")

mad_h = {h: np.mean([abs(raw[S][h] - V["tau_m_at_tangency"]) / V["tau_m_at_tangency"] * 100
                     for S in S_PRINTED]) for h in H_STUDY}
mad_conv = np.mean([abs(conv[S] - V["tau_m_at_tangency"]) / V["tau_m_at_tangency"] * 100
                    for S in S_PRINTED])
drift_h = {h: (raw[200.0][h] - raw[S_PRINTED[0]][h]) / raw[S_PRINTED[0]][h] * 100
           for h in H_STUDY}
drift_conv = (conv[200.0] - conv[S_PRINTED[0]]) / conv[S_PRINTED[0]] * 100
print(f"\\nmean |dev| from {V['tau_m_at_tangency']} over the four S-curves Barkelew drew:")
print(f"  at h = 0.04 (a plausible fixed choice) : {mad_h[0.04]:.3f} %")
print(f"  extrapolated to h -> 0                 : {mad_conv:.3f} %")
print(f"the fall in tau_m from S = {S_PRINTED[0]:.0f} to S = 200, this page's main result:")
print(f"  at h = 0.04                            : {drift_h[0.04]:+.3f} %")
print(f"  extrapolated to h -> 0                 : {drift_conv:+.3f} %")
print(f"  so a fixed h = 0.04 INFLATES the drift by "
      f"{abs(drift_h[0.04] - drift_conv):.3f} points, "
      f"{abs((drift_h[0.04] - drift_conv) / drift_conv) * 100:.0f} % of it.")
print(f"[{time.time() - t_start:.1f} s]")'''))

cells.append(md(r"""Two things to take from that table, and the second is the
reason it is printed at all.

The truncation error is $O(h^2)$ as advertised, and the pair
$(h, h/2) = (0.01, 0.005)$ extrapolates to a value that a second pair
$(0.005, 0.0025)$ confirms to far better than the digits quoted. **Every
$\tau_m$ from here on is that extrapolated value**, so the tangency is reported
as the exactly-defined mathematical object it is, not as an artefact of a step
size nobody varied.

And the bias is **not** neutral with respect to what this page concludes. It
grows monotonically with $S$ — from under a tenth of a percent at $S = 4$ to
more than half a percent at $S = 200$ — which is the *same direction* as the
drift in $\tau_m$ that is this page's main new result. An unconverged grid would
therefore have inflated the conclusion: at $h = 0.04$ the fall from $S = 4$ to
$S = 200$ reads about half a percentage point steeper than it is. The drift
survives convergence with about 95 % of its size intact, but a page that had not
checked would not have known that.

### The diagram

Their Fig. 9, recomputed: $\tau_m/S$ against $N/S$ for the four printed
$S$-values, over the printed axis ranges. The **envelope** and the tangency
points are computed, not read; the open circles are Van Welsenaere and
Froment's own first criterion in these coordinates, $\tau_m/S = 1/S$."""))

cells.append(code('''# tangency and the tau_m = 1 point, on every S
tan_u, tan_t, one_u = {}, {}, {}
for S in S_WIDE:
    u, t = tangency(S)
    tan_u[S], tan_t[S] = u, t
    one_u[S] = u_at_tau(S, 1.0)

# the S-curves themselves
NS_MIN, NS_MAX = V["fig9_NS_min"], V["fig9_NS_max"]
TS_MAX = V["fig9_tauS_max"]

curves = {}
for S in S_PRINTED:
    u_hi = NS_MAX
    u_lo = max(NS_MIN, u_at_tau(S, TS_MAX * S) if TS_MAX * S < TAU_CAP else NS_MIN)
    uu = np.linspace(u_lo, u_hi, 90)
    curves[S] = (uu, np.array([tau_max_phase(u * S, S) for u in uu]) / S)

# The envelope IS the locus of tangency points: each S-curve touches it exactly
# where d/dS vanishes, so no separate construction is needed - and a pointwise
# minimum over a truncated S-range would be an artefact of where it was
# truncated, not the envelope.
env_u = np.array([tan_u[S] for S in S_WIDE])
env_v = np.array([tan_t[S] / S for S in S_WIDE])
env_spline = PchipInterpolator(env_u, env_v)
env_uu = np.linspace(env_u.min(), min(env_u.max(), NS_MAX), 300)

fig, ax = plt.subplots(figsize=(7.2, 5.0))
for S, (uu, vv) in curves.items():
    ax.plot(uu, vv, "k-", lw=1.6)
    j = int(np.argmin(np.abs(vv - 0.385)))
    ax.text(uu[j] + 0.035, 0.378, f"S = {S:.0f}", fontsize=9, va="top", ha="left")
ax.plot(env_uu, env_spline(env_uu), "--", color="tab:red", lw=1.4,
        label="envelope (the locus of tangency points)")
ax.plot([tan_u[S] for S in S_PRINTED], [tan_t[S] / S for S in S_PRINTED],
        "x", color="tab:red", ms=9, mew=2,
        label="tangency = Barkelew's critical point")
ax.plot([one_u[S] for S in S_PRINTED], [1.0 / S for S in S_PRINTED],
        "o", mfc="none", mec="tab:blue", ms=8, mew=1.6,
        label=r"Van Welsenaere-Froment criterion 1, $\\tau_m/S = 1/S$")
ax.set(xlabel="$N/S$", ylabel=r"$\\tau_m/S$", xlim=(NS_MIN, NS_MAX),
       ylim=(0.0, TS_MAX), title="Barkelew's diagram, recomputed (their Fig. 9)")
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
plt.show()
print(f"[{time.time() - t_start:.1f} s]")'''))

cells.append(md(r"""### Does it collapse?

The whole value of the chart is that the tangency points share one $\tau_m$. The
paper says 1.275, for all $S$. Here is $\tau_m$ at the tangency over a range of
$S$ six times wider than the one Barkelew drew, at the converged step size."""))

cells.append(code('''tab = pd.DataFrame({
    "S": S_WIDE,
    "(N/S) tangency": [tan_u[S] for S in S_WIDE],
    "tau_m tangency": [tan_t[S] for S in S_WIDE],
    "dev from 1.275 %": [(tan_t[S] - V["tau_m_at_tangency"]) / V["tau_m_at_tangency"] * 100
                         for S in S_WIDE],
    "(N/S) at tau_m=1": [one_u[S] for S in S_WIDE],
    "gap to the right %": [(one_u[S] - tan_u[S]) / tan_u[S] * 100 for S in S_WIDE],
})
tab["drawn by Barkelew"] = ["yes" if S in S_PRINTED else "" for S in S_WIDE]
print(tab.to_string(index=False, float_format=lambda v: f"{v:8.4f}"))

drawn = tab[tab["drawn by Barkelew"] == "yes"]
mad_drawn = drawn["dev from 1.275 %"].abs().mean()
mad_wide = tab["dev from 1.275 %"].abs().mean()
print(f"\\nover the four S-curves Barkelew drew : mean |dev| {mad_drawn:.3f} %, "
      f"tau_m from {drawn['tau_m tangency'].min():.4f} to {drawn['tau_m tangency'].max():.4f}")
print(f"over S = {min(S_WIDE):.0f} to {max(S_WIDE):.0f}          : mean |dev| {mad_wide:.3f} %, "
      f"tau_m from {tab['tau_m tangency'].min():.4f} to {tab['tau_m tangency'].max():.4f}")
print(f"the gap is positive for all {len(tab)} values of S : "
      f"{bool((tab['gap to the right %'] > 0).all())}"
      f"   (range {tab['gap to the right %'].min():.1f} to "
      f"{tab['gap to the right %'].max():.1f} %)")'''))

cells.append(md(r"""$\tau_m$ at the tangency is not a constant. It rises to a
maximum near $S \approx 4$, at the very left edge of the range Barkelew drew,
and falls monotonically thereafter — so 1.275 is a good summary of his four
curves and a steadily worse one outside them. The panel below shows the drift,
and shows how far outside his range it has to be pushed before the collapse
stops being usable.

This is not a contradiction of anything the paper asserts. Section 5 says *"very
close to"*, which already declines to claim a constant, and the crosses printed
on their Fig. 9 could not have resolved a change of 0.06 in $\tau_m$ across
$S = 4$–32 anyway. What the computation supplies is the size of "very close to",
and the fact that it has a sign and a direction.

The second statement in Section 5 survives intact: the point $\tau_m/S = 1/S$
lies to the right of the tangency for **every** $S$ tested, so Van Welsenaere
and Froment's first criterion is more conservative than Barkelew's everywhere,
not merely on the four curves they could see. The size of that margin is not in
either paper — it collapses from a quarter of $N/S$ at $S = 4$ to a percent at
$S = 200$."""))

cells.append(code('''fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
ax = axes[0]
ax.semilogx(tab.S, tab["tau_m tangency"], "o-", color="k", lw=1.6, ms=5)
ax.axhline(V["tau_m_at_tangency"], color="tab:red", ls="--", lw=1.4,
           label=f"their stated {V['tau_m_at_tangency']}")
ax.axvspan(min(S_PRINTED), max(S_PRINTED), color="tab:blue", alpha=0.10,
           label="the range Barkelew drew")
ax.set(xlabel="$S$", ylabel=r"$\\tau_m$ at the tangency",
       title="the collapse, and where it drifts")
ax.legend(fontsize=8)

ax = axes[1]
ax.semilogx(tab.S, tab["gap to the right %"], "o-", color="tab:blue", lw=1.6, ms=5)
ax.axhline(0.0, color="k", lw=1.0)
ax.axvspan(min(S_PRINTED), max(S_PRINTED), color="tab:blue", alpha=0.10)
ax.set(xlabel="$S$", ylabel=r"$(N/S)_{\\tau_m=1}$ over tangency, %",
       title="how much more conservative criterion 1 is")
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""### What the two criteria predict for a real tube

$N$ depends only on the wall temperature and the tube radius; $S$ carries the
inlet pressure. So Barkelew's criterion, read backwards, gives a critical inlet
pressure for any $T_w$ — which is what the circled crosses on their Fig. 8 are.
The comparison below recomputes both sides: Barkelew's tangency locus with his
exponential rate law, against Van Welsenaere and Froment's first criterion
back-integrated with the **true** Arrhenius law, the quantity their Example 1
prints as 0.01651 atm."""))

cells.append(code('''# Critical S for a given N, from the computed tangency locus.
Ns_tan = np.array([tan_u[S] * S for S in S_WIDE])
S_of_N = PchipInterpolator(np.log(Ns_tan), np.log(S_WIDE))
Ns_one = np.array([one_u[S] * S for S in S_WIDE])
S_of_N_one = PchipInterpolator(np.log(Ns_one), np.log(S_WIDE))


def p0_crit_vwf(Tw, R=None, Bs=1.0, Cs=1.0):
    """Their criterion 1 with the TRUE rate law: back-integrate their Eq. 5
    from the top of the maxima curve until T falls to Tw, i.e. to the inlet.

    `Bs` and `Cs` scale B and C for the blindness tests further down; the model
    is Bs = Cs = 1.
    """
    Cc = (C_BASE if R is None else C_of_R(R)) * Cs
    Bb = B * Bs
    TM = T_M(Tw)
    pM = (TM - Tw) / ((Bb / Cc) * kexp(TM))          # their Eq. 7 at T = T_M

    def dTdp(p, T):                                  # their Eq. 5
        return [-Bb / A + (Cc / A) * (T[0] - Tw) / (p * kexp(T[0]))]

    def hit_Tw(p, T):
        return T[0] - Tw
    hit_Tw.terminal, hit_Tw.direction = True, -1.0
    sol = solve_ivp(dTdp, [pM, pM + 1.0], [TM], rtol=1e-10, atol=1e-12,
                    events=hit_Tw, method="LSODA")
    return float(sol.t_events[0][0])


Tw_grid = np.arange(600.0, 701.0, 5.0)
rows = []
for Tw in Tw_grid:
    N = float(N_group(Tw))
    Sb = float(np.exp(S_of_N(np.log(N))))
    S1 = float(np.exp(S_of_N_one(np.log(N))))
    rows.append((Tw, N, Sb, p0_of_S(Sb, Tw), p0_of_S(S1, Tw), p0_crit_vwf(Tw)))
cmp = pd.DataFrame(rows, columns=["Tw", "N", "S_crit",
                                  "p0 Barkelew", "p0 crit-1 (exp law)",
                                  "p0 VWF exact (true law)"])
cmp["dev Barkelew vs VWF %"] = ((cmp["p0 Barkelew"] - cmp["p0 VWF exact (true law)"])
                                / cmp["p0 VWF exact (true law)"] * 100)
show = cmp[cmp.Tw.isin([600., 625., 650., 675., 700.])]
print(show.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))

i625 = int(np.argmin(np.abs(cmp.Tw.values - 625.0)))
print(f"\\nat Tw = 625 K, their Example 1 prints (p0)_cr,1 = {V['p0_crit_1_at_Tw625']} atm")
print(f"  this page's back-integration of the same criterion : "
      f"{cmp['p0 VWF exact (true law)'][i625]:.5f} atm  "
      f"({abs(cmp['p0 VWF exact (true law)'][i625] - V['p0_crit_1_at_Tw625']) / V['p0_crit_1_at_Tw625'] * 100:.2f} %)")
print(f"  Barkelew's criterion, recomputed                   : "
      f"{cmp['p0 Barkelew'][i625]:.5f} atm  "
      f"({cmp['dev Barkelew vs VWF %'][i625]:+.2f} % against it)")
print(f"\\nover Tw = 600-700 K: mean |dev| "
      f"{cmp['dev Barkelew vs VWF %'].abs().mean():.3f} %, "
      f"worst {cmp['dev Barkelew vs VWF %'].abs().max():.3f} %, "
      f"signed range {cmp['dev Barkelew vs VWF %'].min():+.2f} to "
      f"{cmp['dev Barkelew vs VWF %'].max():+.2f} %")

dev_of_Tw = lambda Tw: (p0_of_S(float(np.exp(S_of_N(np.log(float(N_group(Tw)))))), Tw)
                        - p0_crit_vwf(Tw))
Tw_cross = brentq(dev_of_Tw, 600.0, 700.0, xtol=1e-6)
print(f"the two methods cross at Tw = {Tw_cross:.1f} K "
      f"(bracketed: {dev_of_Tw(600.0):+.2e} at 600 K, {dev_of_Tw(700.0):+.2e} at 700 K)")'''))

cells.append(md(r"""#### The two approximations, separated

Section 5's comparison mixes two different things: Barkelew's *criterion* (the
envelope tangency instead of $\tau_m = 1$) and Barkelew's *rate law* (the
exponential instead of the true Arrhenius factor). Their Fig. 8 puts his values
next to theirs and the text says the two *"agree extremely well"*, but a reader
cannot tell how much of whatever difference remains is which — the paper does
say his values were obtained *"using the modified rate equation"*, so it is not
hiding the mixture, it simply does not decompose it.

Computing the intermediate case — their criterion, his rate law — separates
them exactly, because the critical $p^0$ factorises into the two steps:

$$
\frac{p^0_{\text{Barkelew}}}{p^0_{\text{VWF, true}}}
= \underbrace{\frac{p^0_{\tau_m=1,\ \exp}}{p^0_{\text{VWF, true}}}}_{\text{rate law alone}}
\times
\underbrace{\frac{p^0_{\text{Barkelew}}}{p^0_{\tau_m=1,\ \exp}}}_{\text{criterion alone}}
$$

The second factor is the ratio of the critical $S$ at the tangency to the
critical $S$ at $\tau_m = 1$ at the same $N$: a **purely dimensionless**
quantity, so it does not inherit the reconstruction at all. Only the first
factor does."""))

cells.append(code('''dec = pd.DataFrame({"Tw": cmp.Tw})
dec["rate law alone %"] = (cmp["p0 crit-1 (exp law)"]
                           / cmp["p0 VWF exact (true law)"] - 1.0) * 100
dec["criterion alone %"] = (cmp["p0 Barkelew"]
                            / cmp["p0 crit-1 (exp law)"] - 1.0) * 100
dec["total %"] = cmp["dev Barkelew vs VWF %"]
dec["product of the two %"] = ((1 + dec["rate law alone %"] / 100)
                               * (1 + dec["criterion alone %"] / 100) - 1) * 100
print("what the Barkelew-vs-VWF disagreement is made of")
print(dec[dec.Tw.isin([600., 625., 650., 675., 700.])].to_string(
    index=False, float_format=lambda v: f"{v:9.4f}"))

opp = bool((dec["rate law alone %"] * dec["criterion alone %"] < 0).all())
smaller = bool((dec["total %"].abs() < dec["rate law alone %"].abs()).all()
               & (dec["total %"].abs() < dec["criterion alone %"].abs()).all())
print(f"\\nrate law alone : {dec['rate law alone %'].min():+.2f} to "
      f"{dec['rate law alone %'].max():+.2f} % (conservative everywhere)")
print(f"criterion alone: {dec['criterion alone %'].min():+.2f} to "
      f"{dec['criterion alone %'].max():+.2f} % (the other way everywhere)")
print(f"total          : {dec['total %'].min():+.2f} to {dec['total %'].max():+.2f} %")
print(f"opposite sign at every one of the {len(dec)} wall temperatures: {opp}")
print(f"total smaller in magnitude than BOTH components everywhere : {smaller}")
print(f"the product of the two reproduces the total to "
      f"{(dec['product of the two %'] - dec['total %']).abs().max():.2e} points, "
      f"as it must - the decomposition is exact, not a fit.")'''))

cells.append(code('''fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
ax = axes[0]
ax.plot(cmp.Tw, cmp["p0 VWF exact (true law)"], "k-", lw=2.0,
        label="VWF criterion 1, true rate law (exact)")
ax.plot(cmp.Tw, cmp["p0 Barkelew"], "--", color="tab:red", lw=1.8,
        label="Barkelew's tangency, exponential rate law")
ax.plot(cmp.Tw, cmp["p0 crit-1 (exp law)"], ":", color="tab:blue", lw=1.8,
        label=r"VWF criterion 1 in Barkelew's variables ($\\tau_m=1$)")
ax.plot([625.0], [V["p0_crit_1_at_Tw625"]], "*", color="k", ms=13,
        label="their Example 1, printed")
ax.set(xlabel="$T_w$ (K)", ylabel="critical $p^0$ (atm)", xlim=(600, 700),
       yscale="log", title="the runaway boundary, three ways")
ax.legend(fontsize=8)

ax = axes[1]
ax.plot(dec.Tw, dec["rate law alone %"], ":", color="tab:blue", lw=1.8,
        label="rate law alone (exponential vs true)")
ax.plot(dec.Tw, dec["criterion alone %"], "-.", color="tab:green", lw=1.8,
        label=r"criterion alone (tangency vs $\\tau_m = 1$)")
ax.plot(dec.Tw, dec["total %"], "-", color="tab:red", lw=2.0,
        label="total, Barkelew vs VWF exact")
ax.axhline(0.0, color="k", lw=1.0)
ax.axvline(Tw_cross, color="k", ls=":", lw=1.0)
ax.annotate(f"cross at {Tw_cross:.1f} K", (Tw_cross, 8.0), fontsize=8,
            ha="center")
ax.set(xlabel="$T_w$ (K)", ylabel="deviation, %", xlim=(600, 700),
       title="what the disagreement is made of")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()'''))

cells.append(md(r"""The two methods cross. Barkelew's chart is *more*
conservative than the true-rate criterion at low wall temperature — it allows a
lower inlet pressure — and *less* conservative at high, and the decomposition
says why. The exponential rate law on its own is conservative everywhere and by
more than the total discrepancy; the tangency criterion sits above $\tau_m = 1$
and pushes the other way by a growing amount. The two have opposite sign at
every wall temperature and the total is smaller in magnitude than either
component at all of them. That is the honest reading of *"agree extremely
well"*: the agreement is better than either approximation taken alone, and it is
better for a reason neither paper states.

### How much earlier the simplified rate law runs away

Their Fig. 3 and their Fig. 10 are the same four cases — $p^0 = 0.017$ atm,
$T_w = 625, 626, 627, 628$ K — with the true rate law and with Barkelew's. The
text says the modified equation *"leads to runaway well before the true rate
equation"*. This is what that costs."""))

cells.append(code('''class Tube:
    """The dimensional tube of their Eqs. 3-4, on the contact time z = z'/u.

    `simplified=True` swaps the Arrhenius factor for Barkelew's exponential.
    """

    def __init__(self, p0, Tw, length=1.0, n_z=2000, simplified=False, R=None):
        self.p0, self.Tw, self.simplified = p0, Tw, simplified
        self.C = C_BASE if R is None else C_of_R(R)
        self.kw = float(kexp(Tw))
        self.z_f = np.linspace(0.0, length / P["u"], n_z + 1)
        self.z_c = 0.5 * (self.z_f[:-1] + self.z_f[1:])
        self.shape = (n_z, 2)
        # inlet: p = p0, T = Tw   |   outlet: pure outflow (outward normal)
        self.bc = ({"a": 0.0, "b": 1.0, "d": np.array([p0, Tw])},
                   {"a": 1.0, "b": 0.0, "d": 0.0})
        self.u = np.tile(np.array([p0, Tw]), (n_z, 1))
        self.diverged = False
        cm, cbc = construct_convflux_upwind(self.shape, self.z_f, self.z_c,
                                            self.bc, v=1.0, axis=0)
        dm = construct_div(self.shape, self.z_f, nu=0, axis=0)
        self.jac_const = dm @ cm
        self.g_const = dm @ cbc
        self.numjac = NumJac(self.shape)

    def k(self, T):
        if self.simplified:
            return self.kw * np.exp(np.clip(
                a_exp * (T - self.Tw) / self.Tw ** 2, -50.0, 50.0))
        return np.exp(np.clip(b_exp - a_exp / T, -50.0, 50.0))

    def reaction(self, u):
        p, T = u[..., 0], u[..., 1]
        r = p * self.k(T)
        return np.stack([-A * r, B * r - self.C * (T - self.Tw)], axis=-1)

    def residual(self, u):
        g_r, j_r = self.numjac(self.reaction, u)
        g = self.g_const + self.jac_const @ u.reshape((-1, 1)) - g_r.reshape((-1, 1))
        return g, self.jac_const - j_r

    def solve(self, maxfev=200):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                res = newton(self.residual, self.u, maxfev=maxfev)
            except (ValueError, RuntimeError):
                self.diverged = True
                return None
        self.u = res.x.reshape(self.shape)
        self.diverged = not np.all(np.isfinite(self.u))
        return res


def profile_ivp(p0, Tw, simplified, length=1.0, T_stop=None):
    """The same reactor marched as an initial-value problem, stopping at the
    top of the paper's own temperature axis."""
    T_stop = V["fig10_T_max"] if T_stop is None else T_stop
    kw = float(kexp(Tw))
    kf = ((lambda T: kw * np.exp(a_exp * (T - Tw) / Tw ** 2)) if simplified
          else (lambda T: np.exp(b_exp - a_exp / T)))

    def rhs(z, w):
        r = w[0] * kf(w[1])
        return [-A * r, B * r - C_BASE * (w[1] - Tw)]

    def hot(z, w):
        return w[1] - T_stop
    hot.terminal, hot.direction = True, 1.0
    sol = solve_ivp(rhs, [0.0, length / P["u"]], [p0, Tw], rtol=1e-10,
                    atol=1e-12, events=hot, method="LSODA", dense_output=True)
    zz = np.linspace(0.0, sol.t[-1], 800)
    return zz * P["u"], sol.sol(zz)[1], (float(sol.t_events[0][0] * P["u"])
                                         if sol.t_events[0].size else None)


p0_fig = V["fig10_p0"]
Tw_fig = [V[f"fig10_Tw_{i}"] for i in (1, 2, 3, 4)]

fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4), sharey=True)
rows = []
for Tw, col in zip(Tw_fig, ["tab:blue", "tab:green", "tab:orange", "tab:red"]):
    zt, Tt, _ = profile_ivp(p0_fig, Tw, False)
    zs, Ts, zblow = profile_ivp(p0_fig, Tw, True)
    axes[0].plot(zt, Tt, color=col, lw=1.7, label=f"$T_w$ = {Tw:.0f} K")
    axes[1].plot(zs, Ts, color=col, lw=1.7, label=f"$T_w$ = {Tw:.0f} K")
    m = Tube(p0_fig, Tw, n_z=2000, simplified=False)
    m.solve()
    rows.append((Tw, float(m.u[:, 1].max()) if not m.diverged else np.inf,
                 float(Tt.max()), zblow))
for ax, ttl in zip(axes, ["true rate law (their Fig. 3)",
                          "Barkelew's exponential (their Fig. 10)"]):
    ax.set(xlabel="$z'$ (m)", xlim=(0.0, V["fig10_z_max"]),
           ylim=(600.0, V["fig10_T_max"]), title=ttl)
    ax.legend(fontsize=8, loc="upper left")
axes[0].set_ylabel("$T$ (K)")
axes[0].text(0.55, 620, f"$p^0$ = {p0_fig} atm", fontsize=9)
fig.tight_layout()
plt.show()

f3 = pd.DataFrame(rows, columns=["Tw", "hot spot, pymrm", "hot spot, marched",
                                 "z' at 800 K, exponential law (m)"])
f3["dev %"] = ((f3["hot spot, pymrm"] - f3["hot spot, marched"]).abs()
               / f3["hot spot, marched"] * 100)
print(f3.to_string(index=False, float_format=lambda v: f"{v:9.4f}"))
print("\\nWith the true rate law every case has a bounded hot spot inside the tube.")
print("With Barkelew's exponential every one of the same four reaches the top of")
print("the paper's temperature axis before the end of the tube - which is what")
print("their Fig. 10 shows, and why his chart is conservative.")'''))

cells.append(md(r"""### The one printed object that can discriminate the reconstruction

Everything checked so far lives either inside the reduced system — which
contains no reference temperature at all — or compares two criteria against each
other. None of it can settle the one step of the reconstruction that the paper
does not print: the temperature at which the exponent is linearised.

Their **Fig. 10** can, because it is drawn with the simplified rate law, and the
position at which each curve leaves the top of the printed 800 K axis depends on
that choice directly. The previous cell already computed those four positions.
Comparing them is the decisive test, and it is the only one available.

> **Review status: not yet page evidence.** The four positions on the printed
> figure are *curve positions read off a scan*, not values the paper prints, so
> under this repository's rules they are a figure digitisation and need the
> maintainer's visual review before they can count. They were measured by
> `queue_cases/D2.1/review/extract_figure10.py`, which writes a numbered overlay
> for that review; they are **not** in either committed CSV and nothing else on
> this page depends on them. The case file carries a non-blocking follow-up
> asking for the review. Until it happens, read the next table as an indication.
> The reconstruction stands on the printed sentences either way — this only
> makes it decisive rather than merely consistent."""))

cells.append(code('''# ---------------------------------------------------------------------------
# QUARANTINE: the only numbers on this page taken off a curve rather than out
# of the text. Measured on a 600 dpi render of journal page 1512 by
# queue_cases/D2.1/review/extract_figure10.py, which fits each of the six marks
# (two frame edges, four curves) linearly over 72 clean rows and extrapolates to
# the top frame row where T = 800 K actually is; worst fit residual 1.5 px. The
# frame there is x = 575.5 px (z' = 0) and x = 2154.3 px (z' = 1.0 m, the
# printed axis limits), the four curves cross it at x = 1220.0, 1321.5, 1435.6
# and 1768.2 px, and the printed 0.5 tick lands 0.03 % of full scale from the
# frame midpoint. AWAITING MAINTAINER VISUAL REVIEW of the overlay that script
# writes - not in any CSV, no sidecar, and no other cell reads this dict.
FIG10_MEASURED_PENDING_REVIEW = {628.0: 0.4082, 627.0: 0.4725,
                                 626.0: 0.5448, 625.0: 0.7555}


def blowup_z(Tw, k_wall, slope, p0=None, length=None, T_stop=None):
    """z' at which T reaches the top of the printed axis, for a rate law
    k(T) = k_wall * exp(slope * (T - Tw)). Radau, so that the candidates that
    make the problem violently stiff are integrated on the same footing."""
    p0 = V["fig10_p0"] if p0 is None else p0
    length = V["fig10_z_max"] if length is None else length
    T_stop = V["fig10_T_max"] if T_stop is None else T_stop

    def rhs(z, w):
        r = w[0] * k_wall * np.exp(np.clip(slope * (w[1] - Tw), -50.0, 50.0))
        return [-A * r, B * r - C_BASE * (w[1] - Tw)]

    def hot(z, w):
        return w[1] - T_stop
    hot.terminal, hot.direction = True, 1.0
    sol = solve_ivp(rhs, [0.0, length / P["u"]], [p0, Tw], rtol=1e-10,
                    atol=1e-12, events=hot, method="Radau")
    z = float(sol.t_events[0][0] * P["u"]) if sol.t_events[0].size else None
    return z, float(sol.y[1].max())


def candidates(Tw):
    """The linearisations of ln k that every printed sentence permits.

    Each is ln k(T) ~ ln k(T_r) + (T - T_r)*a/T_r^2 rewritten as k_w,eff e^tau
    with tau = a (T - T_w)/D, so all of them vanish at the wall and are linear
    in T - T_w, as Section 5 requires. Only D differs.
    """
    TM = float(T_M(Tw))
    kw = float(kexp(Tw))
    return {
        "tangent at T_w  (this page)": (kw, a_exp / Tw ** 2),
        "tangent at T_M": (float(kexp(TM)) * np.exp(-a_exp * (TM - Tw) / TM ** 2),
                           a_exp / TM ** 2),
        "chord T_w to T_M": (kw, a_exp / (Tw * TM)),
        "square dropped, a(T-T_w)/T_w": (kw, a_exp / Tw),
    }


names = list(candidates(625.0))
rows, hot_rows = [], []
for Tw in Tw_fig[::-1]:                     # 628, 627, 626, 625
    row = {"Tw": Tw, "figure (pending review)": FIG10_MEASURED_PENDING_REVIEW[Tw]}
    hot = {"Tw": Tw}
    for nm, (kw, sl) in candidates(Tw).items():
        z, Tmax = blowup_z(Tw, kw, sl)
        row[nm] = z if z is not None else np.nan
        hot[nm] = Tmax
    rows.append(row)
    hot_rows.append(hot)
fig10 = pd.DataFrame(rows)
hottest = pd.DataFrame(hot_rows)
print("z' (m) at which T reaches the printed 800 K axis, p0 = "
      f"{V['fig10_p0']} atm; NaN = no runaway within the printed "
      f"{V['fig10_z_max']:.0f} m")
print(fig10.to_string(index=False, float_format=lambda v: f"{v:8.4f}"))

devs = pd.DataFrame({"Tw": fig10.Tw})
for nm in names:
    devs[nm] = (fig10[nm] - fig10["figure (pending review)"]) \\
        / fig10["figure (pending review)"] * 100
print("\\ndeviation from the measured positions, %")
print(devs.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))
print()
for nm in names:
    n_ok = int(devs[nm].notna().sum())
    if n_ok == 0:
        print(f"  {nm:30s} no runaway in any of the four cases")
    else:
        print(f"  {nm:30s} {n_ok}/4 run away, |dev| {devs[nm].abs().min():.1f} "
              f"to {devs[nm].abs().max():.1f} %")
print("\\nhighest T reached inside the printed tube, K (800 = ran away)")
print(hottest.to_string(index=False, float_format=lambda v: f"{v:8.2f}"))
lsoda_z = list(f3["z' at 800 K, exponential law (m)"])[::-1]   # 628 -> 625
radau_z = list(fig10["tangent at T_w  (this page)"])
print("\\nNot an integrator artefact: the LSODA march that drew the figure above")
print("  gives " + ", ".join(f"{z:.4f}" for z in lsoda_z) + " m")
print("  Radau here gives " + ", ".join(f"{z:.4f}" for z in radau_z) + " m")
print(f"  worst difference {max(abs(a - b) for a, b in zip(lsoda_z, radau_z)):.1e} m")'''))

cells.append(md(r"""The wall tangent reproduces all four printed positions,
three of them to under one percent and the fourth to 2.4 %. Every alternative
that the printed sentences permit fails *qualitatively*, not marginally.
Linearising at $T_M$, or on the chord from $T_w$ to $T_M$, weakens the
temperature dependence enough that most of the four cases never reach 800 K at
all inside the printed 1 m of tube — three of four for the $T_M$ tangent, two of
four for the chord — and the ones that do survive are off by 29 to 59 %.
Dropping the square goes the other way and makes the exponent so steep that the
reaction cannot outrun the wall cooling at all: the hottest of the four gets
about 1.5 K above the wall and nothing runs away.

So the reconstruction is not merely consistent with the printed text — subject
to the maintainer's review of those four positions, it is the only member of its
family that reproduces the paper's own figure."""))

cells.append(md(r"""## Validation

**Ranked before any code was written.** The paper offers no experiment, so the
highest available route is *a stated numerical result in the text* — and
Section 5 has exactly one, $\tau_m = 1.275$ at the tangency, plus one stated
identity, $\tau_m = 1$ for the other criterion. **No figure was digitised into
the page's data**; the only things taken from Figs. 9 and 10 into the CSVs are
the labels printed inside them, and the four Fig. 10 positions above are
quarantined and awaiting review.

| # | check | what it tests | what it cannot see |
|---|---|---|---|
| 1 | $\tau_m$ at the envelope tangency against their 1.275 | the exponential temperature dependence, and the choice of $N/S$ and $\tau_m/S$ as the collapsing groups | the dimensional definitions of $N$ and $S$ — the reduced system does not contain them — and, as measured below, all but a little of the reaction order |
| 2 | the $\tau_m = 1$ point lies right of the tangency, every $S$ | the sign and ordering claim of Section 5 | magnitude; it is a sign test |
| 3 | Barkelew's critical $p^0$ against the criterion 1 back-integration | the dimensional group definitions and the reference temperature | **nothing about $B$**: it is invariant under any rescaling of $B$ down to floating point, structurally and not by luck, hence blind to $(-\Delta H)$, $\rho_b$ and $p_B^0$ — measured below. It is also nearly blind to $C$, including the $c_p\rho_g$ trap, which check 5 catches instead |
| 4 | pymrm against the phase-plane quadrature | the operator assembly, the discretisation — **and the reduced model as the finite-volume route transcribes it**: the two routes do *not* share a `reaction` function, and a flipped cooling sign or a 10 % error in `BarkelewTube.reaction` alone sends this check to `inf`, measured below | a defect written identically into *both* transcriptions, and everything dimensional — it lives entirely inside the reduced system |
| 5 | the parameter reading against four printed derived quantities | a lost decimal point anywhere in the table, and the $c_p\rho_g$ trap | — |

Checks 1 and 3 are the ones that carry weight, and they are **blind to different
things**: 1 lives entirely inside the dimensionless system and cannot see a
wrong $N$ or $S$; 3 is built on those definitions but cannot see anything that
cancels between its two routes. The table below measures both claims rather than
asserting them."""))

cells.append(code('''# Check 4: the pymrm reactor against the phase-plane quadrature, at the four
# tangency points, plus the observed order of convergence.
print("4. pymrm (finite volume in zeta) vs phase-plane quadrature (adaptive, no grid)")
conv_rows = []
for S in S_PRINTED:
    N = tan_u[S] * S
    ref = tau_max_phase(N, S)
    errs = []
    for n in (200, 400, 800, 1600):
        errs.append(abs(tau_max_pymrm(N, S, n_z=n) - ref) / ref * 100)
    order4 = np.polyfit(np.log([200, 400, 800, 1600]), np.log(errs), 1)[0]
    conv_rows.append((S, ref, *errs, -order4))
conv4 = pd.DataFrame(conv_rows, columns=["S", "tau_m reference", "err 200 %",
                                         "err 400 %", "err 800 %", "err 1600 %",
                                         "observed order"])
print(conv4.to_string(index=False, float_format=lambda v: f"{v:9.4f}"))
cross_mean = conv4["err 1600 %"].mean()
print(f"\\nmean |dev| at n_zeta = 1600: {cross_mean:.3f} %, "
      f"worst {conv4['err 1600 %'].max():.3f} %; first-order upwind, so O(h) is "
      f"the expected order and {conv4['observed order'].mean():.2f} is what comes out.")'''))

cells.append(md(r"""### Make sure the checks can fail

Each check is now given a defect it is *presented* as guarding against, and one
it is not. A number that does not move is decoration; the point of printing this
table is to say which cells are which — and one row below is printed precisely
*because* it cannot fail."""))

cells.append(code('''def check1(tm=None, Ss=S_PRINTED):
    """mean |dev| of tau_m at the tangency from their 1.275, over the four
    S-curves Barkelew drew. Returns None when the envelope ceases to exist."""
    d = []
    for S in Ss:
        _, t = tangency(S, tm)
        if t is None:
            return None
        d.append(abs(t - V["tau_m_at_tangency"]) / V["tau_m_at_tangency"] * 100)
    return float(np.mean(d))


def check3(T2=None, kk=None, tm=None, Bs=1.0, Cs=1.0):
    """mean |dev| of Barkelew's critical p0 from the criterion 1 exact value.

    T2(Tw) is the squared reference temperature in tau = a (T - Tw)/T2, and
    kk(Tw) the rate coefficient at the wall under that same linearisation.
    Bs, Cs scale B and C in BOTH routes, which is how the blindness rows are
    measured.
    """
    T2 = (lambda Tw: Tw ** 2) if T2 is None else T2
    kk = (lambda Tw: float(kexp(Tw))) if kk is None else kk
    us = []
    for S in S_WIDE:
        u, _ = tangency(S, tm)
        if u is None:
            return None
        us.append(u)
    S_of = PchipInterpolator(np.log(np.array(us) * np.array(S_WIDE)), np.log(S_WIDE))
    d = []
    for Tw in Tw_grid:                    # the same 600-700 K grid as above
        N = C_BASE * Cs / (A * kk(Tw))
        Sc = float(np.exp(S_of(np.log(N))))
        p0 = Sc * A * T2(Tw) / (a_exp * B * Bs)
        pv = p0_crit_vwf(Tw, Bs=Bs, Cs=Cs)
        d.append(abs(p0 - pv) / pv * 100)
    return float(np.mean(d))


def check4(n_z=1600, nu=0, cool=1.0, s_fac=1.0):
    d = []
    for S in S_PRINTED:
        N = tan_u[S] * S
        ref = tau_max_phase(N, S)
        d.append(abs(tau_max_pymrm(N, S, n_z=n_z, nu=nu, cool=cool, s_fac=s_fac)
                     - ref) / ref * 100)
    return float(np.mean(d))


def attempt(fn):
    """Run a defect and report HOW it fails, not merely that it did."""
    try:
        v = fn()
    except ValueError:
        return "brentq bracketing fails (ValueError)"
    except RuntimeError:
        return "solver raises (RuntimeError)"
    if v is None:
        return "no sign change: the envelope is gone"
    if not np.isfinite(v):
        return "the pymrm solve diverges (inf)"
    return f"{v:.3f} %"


lin = lambda t: 1.0 + t          # exp(tau) -> 1 + tau
T_M2 = lambda Tw: T_M(Tw) ** 2
k_at_TM = lambda Tw: float(np.exp(b_exp - a_exp / T_M(Tw))
                           * np.exp(-a_exp * (T_M(Tw) - Tw) / T_M(Tw) ** 2))

base1, base3, base4 = check1(), check3(), check4()
base5 = float(check5()["dev %"].max())
base = {1: base1, 3: base3, 4: base4, 5: base5}
print(f"baseline   check 1 {base1:.3f} %   check 3 {base3:.3f} %   "
      f"check 4 {base4:.3f} %   check 5 {base5:.3f} %\\n")

c3_TM = check3(T2=T_M2, kk=k_at_TM)

defects = [
    ("exp(tau) -> 1 + tau in the reduced ODE",
     lambda: check1(lambda N, S, **kw: tau_max_phase(N, S, expo=lin)), 1),
    ("cooling term sign flipped, -N tau -> +N tau",
     lambda: check1(lambda N, S, **kw: tau_max_phase(N, S, cool=-1.0)), 1),
    ("reaction second order in y instead of first",
     lambda: check1(lambda N, S, **kw: tau_max_phase(N, S, order=2.0)), 1),
    ("linearise the rate about T_M instead of T_w",
     lambda: c3_TM, 3),
    ("cooling sign flipped in BarkelewTube.reaction ONLY",
     lambda: check4(cool=-1.0), 4),
    ("S 10 % wrong in BarkelewTube.reaction ONLY",
     lambda: check4(s_fac=1.1), 4),
    ("nu = 1 (cylindrical) in construct_div",
     lambda: check4(nu=1), 4),
    ("n_zeta = 25 cells",
     lambda: check4(n_z=25), 4),
    ("c_p multiplied by rho_g (the page's own headline trap)",
     lambda: float(check5(cp_fac=P["rho_g"])["dev %"].max()), 5),
]
out = [(name, which, f"{base[which]:.3f} %", attempt(fn))
       for name, fn, which in defects]
brk = pd.DataFrame(out, columns=["injected defect", "check", "before", "after"])
print("defects each check is presented as guarding against")
print(brk.to_string(index=False))'''))

cells.append(code('''# ---- the row that cannot fail, printed as such ------------------------------
# tau = a(T - Tw)/Tw with the square dropped enters check 3 ONLY through the
# factor T2(Tw) in p0 = Sc*A*T2(Tw)/(a*B). kk, N and the interpolated Sc are
# untouched, so p0 is multiplied by exactly 1/Tw and the deviation is
# 1 - p0/(Tw * p_ref) BY CONSTRUCTION. It is arithmetic, not evidence.
sq_measured = check3(T2=lambda Tw: Tw)
sq_predicted = float(np.mean(np.abs(cmp["p0 Barkelew"] / cmp.Tw
                                    - cmp["p0 VWF exact (true law)"])
                             / cmp["p0 VWF exact (true law)"] * 100))
print("a row that CANNOT fail, printed so it is not mistaken for evidence")
print(pd.DataFrame([{
    "injected defect": "tau = a(T-T_w)/T_w, the square dropped",
    "check": 3,
    "before": f"{base3:.3f} %",
    "after": f"{sq_measured:.3f} %",
    "predicted from the baseline alone, no ODE re-run": f"{sq_predicted:.3f} %",
}]).to_string(index=False))
print(f"the two agree to {abs(sq_measured - sq_predicted):.2e} points, because "
      f"they are the same arithmetic:")
print(f"  p0 is multiplied by exactly 1/T_w ~ 1/{cmp.Tw.mean():.0f}, so the "
      f"deviation is forced to be about {100 * (1 - 1 / cmp.Tw.mean()):.2f} %")
print("  whatever the model, the groups or the criterion happen to be. This row")
print("  is a units check on one factor; it tests nothing about the reconstruction.")
print("  The row that DOES test the reconstruction is the T_M linearisation above,")
print(f"  {base3:.3f} % -> {c3_TM:.3f} %, and the Fig. 10 comparison earlier on")
print("  the page.")'''))

cells.append(code('''# ---- what each check is blind to, measured -----------------------------------
blind = []
for name, fn, which, claim in [
    ("B x 1.1, i.e. (-dH), rho_b or p_B0 10 % wrong",
     lambda: check3(Bs=1.1), 3, "check 3 sees the heat of reaction"),
    ("C x 1.2, i.e. U 20 % wrong",
     lambda: check3(Cs=1.2), 3, "check 3 sees the heat transfer coefficient"),
    ("c_p multiplied by rho_g (C /= 1.293)",
     lambda: check3(Cs=1.0 / P["rho_g"]), 3, "check 3 catches the c_p trap"),
    ("reaction second order in y instead of first",
     lambda: check1(lambda N, S, **kw: tau_max_phase(N, S, order=2.0)), 1,
     "the 1.275 pins the reaction order"),
]:
    v = fn()
    blind.append((name, which, f"{base[which]:.4f} %", f"{v:.4f} %",
                  f"{abs(v - base[which]) / base[which] * 100:.3g} %", claim))
bl = pd.DataFrame(blind, columns=["injected defect", "check", "before", "after",
                                  "relative move", "claim this page does NOT make"])
print("what the checks are blind to - measured, and named as claims not made")
print(bl.to_string(index=False))
print()
print("Check 3 does not move at all under a rescaling of B, and that is")
print("structural, not luck: p0(Barkelew) is proportional to 1/B, and")
print("substituting p = q/B removes B from their Eq. 5 and from p_M as well, so")
print("p0(VWF) is proportional to 1/B too. The ratio cannot move; what is left")
print("in the table is floating point. Check 3 therefore says nothing whatever")
print("about (-dH), rho_b or p_B0.")
print()
print("It is also nearly blind to C, including the c_p x rho_g trap this page")
print("inherits from D2.2 - which lands inside anything a reader would call")
print("agreement. That trap is caught by CHECK 5, not check 3: the recomputed")
print(f"ln K goes from {float(check5()['recomputed'][1]):.4f} to "
      f"{float(check5(cp_fac=P['rho_g'])['recomputed'][1]):.4f} against the "
      f"printed {V['lnK_at_Tw625']},")
print("which is why the page says which check is not silent about it.")
print()
print("And what check 1 cannot see: both group defects change only the map from")
print("the dimensional parameters to (N, S). check 1 never reads that map - its")
print("only inputs are S and the reduced right-hand side - so those two defects")
print(f"leave it the IDENTICAL computation, {base1:.4f} %. A page quoting the")
print("1.275 agreement as evidence that the groups are right would be quoting a")
print("number that cannot see a wrong group.")
print()
print("and what the nu = 1 defect actually does to the pymrm solve:")
tm_nu1, prof_nu1 = tau_max_pymrm(tan_u[S_PRINTED[0]] * S_PRINTED[0], S_PRINTED[0],
                                 n_z=400, nu=1, return_profile=True)
print(f"  tau_max = {tm_nu1:.3e}: the profile collapses to tau = 0 everywhere,")
print("  because a cylindrical divergence has a 1/r factor and the first face")
print("  of this grid sits at zeta = 0. The 100 % is a degenerate solve, not a")
print("  saturated metric.")'''))

cells.append(md(r"""Read the tables as follows.

**Check 1 has power over the form of the reduced system, and less of it than it
looks.** Replacing $\mathrm{e}^{\tau}$ by $1+\tau$ destroys the envelope
altogether — the $S$-curves stop having one, so there is no sign change to
bracket and no criterion — and flipping the sign of the cooling term breaks the
bracketing search one step earlier still, which is why the two failures are
reported differently: one is "no envelope", the other is a `ValueError` out of
`brentq`. But making the reaction second order in $y$ leaves check 1 *better*
than the baseline. That is not noise; it is the table saying that
$\tau_m \approx 1.275$ at the tangency is a property of the **exponential
temperature dependence**, not of the reaction order. A later cell measures how
little the order matters.

And moving the linearisation from $T_w$ to $T_M$ leaves check 1 not merely
unchanged but *the identical computation* — the reduced system contains no
reference temperature and no dimensional parameter at all.

**Check 3 sees the reference temperature**, and that is its job: the $T_M$
linearisation moves it about fourfold. But it does not move at all under a
rescaling of $B$ — $B$ cancels out of both of its routes analytically — and it
barely moves under $C$. Those are not weaknesses to be glossed; they are the
reason check 5 is on the page at all, and check 5 is what catches the
$c_p\rho_g$ trap, which check 3 leaves inside anything a reader would call
agreement.

**Check 4 sees more than the discretisation.** The two routes transcribe the
reduced model independently — `tau_max_phase` writes $\mathrm{d}\tau/\mathrm{d}y$
in its own local `rhs`, `BarkelewTube.reaction` writes the $\zeta$-form
separately — so injecting a defect into the finite-volume side alone is a real
test, and it fails loudly: both the flipped cooling sign and a 10 % error in $S$
send the pymrm solve to `inf`. What it cannot see is a defect written the same
way into both, and anything dimensional. The `nu = 1` row reads exactly 100 %
because the solve becomes degenerate rather than merely inaccurate: a
cylindrical divergence carries a $1/r$ factor and the first face of this grid is
at $\zeta = 0$, so the whole profile collapses to zero. The last cell prints
that value so the 100 % is not mistaken for a saturated metric.

**One row is printed because it cannot fail.** Dropping the square multiplies
$p^0$ by exactly $1/T_w$ and nothing else, so its 99.8 % is forced by
arithmetic — as the cell demonstrates by reproducing it from the baseline
numbers without re-running a single ODE. It is a scale check on one factor and
it is evidence for nothing.

None of the checks is redundant, and none of them is sufficient. That is why all
of them are printed."""))

cells.append(code('''# How much of the 1.275 belongs to the exponential and how much to the
# first-order kinetics? Re-derive the tangency for three reaction orders.
ord_t, ord_u = [], []
for n in (0.5, 1.0, 2.0):
    tm_n = lambda N, S, n=n, **kw: tau_max_phase(N, S, order=n)
    tg = [tangency(S, tm_n) for S in S_PRINTED]
    ord_t.append([n] + [t for _, t in tg])
    ord_u.append([n] + [u for u, _ in tg])
cols = ["reaction order"] + [f"S={S:.0f}" for S in S_PRINTED]
print("tau_m at the tangency (the ordinate of Barkelew's critical point)")
print(pd.DataFrame(ord_t, columns=cols).to_string(
    index=False, float_format=lambda v: f"{v:7.4f}"))
print("\\nN/S at the tangency (its abscissa)")
print(pd.DataFrame(ord_u, columns=cols).to_string(
    index=False, float_format=lambda v: f"{v:7.4f}"))
tv = np.array(ord_t)[:, 1:]
uv = np.array(ord_u)[:, 1:]
print(f"\\nOver reaction orders 0.5 to 2 and S = 4 to 32 the tangency hot spot stays")
print(f"inside {tv.min():.3f} to {tv.max():.3f}, a spread of {100*(tv.max()-tv.min())/tv.min():.1f} %"
      f" - no wider than the spread with")
print(f"S alone. The abscissa is a different matter: at S = 4 it moves from "
      f"{uv[2,0]:.3f} to {uv[0,0]:.3f},")
print(f"a factor of {uv[0,0]/uv[2,0]:.2f}. So their 1.275 characterises the exponential")
print("temperature dependence and says almost nothing about the reaction order,")
print("while the chart itself is order-specific. Check 1 should be read that way")
print("and no further.")'''))

cells.append(code('''report_agreement(PAGE, {
    "tau_m_tangency_mad_pct_drawn": float(mad_drawn),
    "tau_m_tangency_mad_pct_wide": float(mad_wide),
    "tau_m_tangency_S4": float(tan_t[S_PRINTED[0]]),
    "tau_m_tangency_S8": float(tan_t[S_PRINTED[1]]),
    "tau_m_tangency_S16": float(tan_t[S_PRINTED[2]]),
    "tau_m_tangency_S32": float(tan_t[S_PRINTED[-1]]),
    "tau_m_tangency_S200": float(tan_t[200.0]),
    "barkelew_vs_vwf_mean_dev_pct": float(cmp["dev Barkelew vs VWF %"].abs().mean()),
    "barkelew_vs_vwf_worst_dev_pct": float(cmp["dev Barkelew vs VWF %"].abs().max()),
    "rate_law_alone_dev_pct_600K": float(dec["rate law alone %"].iloc[0]),
    "rate_law_alone_dev_pct_700K": float(dec["rate law alone %"].iloc[-1]),
    "criterion_alone_dev_pct_600K": float(dec["criterion alone %"].iloc[0]),
    "criterion_alone_dev_pct_700K": float(dec["criterion alone %"].iloc[-1]),
    "p0_crit_vwf_at_625K": float(cmp["p0 VWF exact (true law)"][i625]),
    "p0_crit_barkelew_at_625K": float(cmp["p0 Barkelew"][i625]),
    "cross_method_mean_pct": float(cross_mean),
})
print(f"\\ntotal notebook compute: {time.time() - t_start:.1f} s")'''))

cells.append(md(r"""## What pymrm adds

**The envelope is computed rather than drawn.** Barkelew's construction is a
tangency between a curve and the envelope of a family, and in 1959 both had to
be drawn and the tangency judged by eye — which is why Section 5 can only report
that $\tau_m$ there was *"very close to"* one number, and why the crosses on
their Fig. 9 sit where a draughtsman put them. Solving the reduced system on
demand makes the envelope the *locus of tangency points* — one bracketed root
per $S$, at a step size extrapolated until the answer no longer depends on it —
so the tangency $\tau_m$ becomes a quantity with a value at every $S$ rather
than a constant asserted for four of them. That is what makes the drift visible
at all, and the step-size study is what makes the drift trustworthy: the naive
fixed step biases it in exactly the direction that would have flattered the
conclusion.

**The two criteria can be put in one picture in dimensional variables, and then
taken apart.** Section 5's comparison is qualitative in the text and graphical
in Fig. 8, and mixes rate laws — Barkelew's values use his exponential, theirs
use the true Arrhenius factor. Computing the intermediate case (their criterion,
his rate law) factorises the difference exactly into a rate-law part and a
criterion part, and the answer is that they have opposite signs at every wall
temperature and partly cancel. Neither paper states that, and it is why "agree
extremely well" is true.

**And the range is no longer set by what fits on a sheet of paper.** Barkelew
drew four curves; the collapse can be tested over a range six times wider at no
cost, which is where a similarity reduction is supposed to be interrogated.

**What this page does not add.** It does not improve on either method. Both are
cheap; the reactor here is a two-equation ODE and even the full sweep is
seconds. Nothing here is validated against a measurement, and neither source
paper contains one."""))

cells.append(md(r"""## Reuse

The dimensionless system is the *generic* first-order exothermic tube, so the
diagram is reusable as-is: compute $N = C/(A k_w)$ and $S = a(\Delta T)_{ad}/T_w^2$
for the tube in hand and read the critical $S$ off the tangency locus. Nothing
in `BarkelewTube` is specific to this chemistry.

Three ways to extend it:

- **A different reaction order.** Change `reaction` to $y^n \mathrm{e}^{\tau}$.
  The validation section already does this for $n = \tfrac12, 1, 2$, and the
  answer is that the tangency $\tau_m$ barely moves — the chart is far more
  general than "first order" suggests, though the *position* of the tangency in
  $N/S$ does shift. The rest of the machinery is untouched.
- **Keep the true Arrhenius factor.** Replace $\mathrm{e}^{\tau}$ by
  $\exp[\tau/(1+\tau/\gamma)]$ with $\gamma = a/T_w$ and the family gains a
  third parameter. The collapse then holds only asymptotically as
  $\gamma \to \infty$, and the machinery here measures how fast.
- **Add axial dispersion.** Build a `construct_grad`/`construct_div` diffusion
  operator, add it to `jac_const`, and change the inlet to a Danckwerts
  condition. The boundary-condition dicts already use the outward normal, so
  only `a` and `b` change.

For the explicit closed-form criteria in dimensional variables, and the full
comparison of their two criteria over the operating plane, see
[`D2.2`](../D2.2-van-welsenaere-froment-runaway/). This page and that one are
two answers to the same question, and they disagree by a few percent in a
direction that depends on the wall temperature."""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata.kernelspec = {"display_name": "Python 3", "language": "python",
                          "name": "python3"}
nb.metadata.language_info = {"name": "python", "version": "3.11"}
nbf.write(nb, "index.ipynb")
print("wrote index.ipynb with", len(cells), "cells")
