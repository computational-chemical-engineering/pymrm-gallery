#!/usr/bin/env python3
"""Generate index.ipynb for page B3.2. Run from the page directory."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = []

cells.append(md(r"""---
title: "The grain model: a solid made of shrinking cores, and the exact penetration law inside it"
description: "Szekely and Evans replaced the shrinking-core model's empirical rate constant with structure: porosity, grain size, two diffusivities and a rate constant that no longer depends on the solid. The page rebuilds their grain model with pymrm, derives the exact closed-form penetration law their quasi-steady equations imply, and shows the model collapsing onto the published B3.1 shrinking core at the exposed face."
categories: [sec:B, struct:S3, struct:S12, tier:T0, data:tier6, phase:gas-solid]
date: 2026-08-05
---

# The grain model: a solid made of shrinking cores, and the exact penetration law inside it

**Catalog ID:** `B3.2` · **Structures:** `S3` (1-D reaction–diffusion BVP) + `S12`
(moving boundary) · **Tier:** T0

A porous ore pellet is not a solid lump with a sharp reaction front. Szekely and
Evans' 1970 paper models it as what it is under the microscope: a packing of
small dense **grains**, each one a shrinking core, all of them fed by gas
diffusing through the interstices. Reaction then happens in a **diffuse zone**
whose thickness the shrinking-core picture cannot even express — and the rate
"constant" measured by fitting a shrinking-core model stops being a property of
the chemistry and becomes a property of the microstructure.

Their motivation is stated in the abstract: the shrinking core model buries
structural parameters "into an empirical reaction rate constant", whereas here
the reaction is described "in terms of the porosity, grain size, gas phase and
solid state diffusivities and a heterogeneous reaction rate constant, **which is
now independent of structure**"."""))

cells.append(md(r"""## Background

The classical description of a gas–solid reaction with a moving boundary is the
shrinking core model — the gallery's published
[`B3.1`](../B3.1-shrinking-core/): three resistances in series and a sharp
interface. Szekely and Evans open by naming its two shortcomings (their §1):
the reacted and unreacted zones "are not necessarily separated by a sharp
boundary, in fact under certain conditions the reaction is found to occur in a
diffuse zone, the thickness of which may be of the same order as the size of
the reacting specimen"; and the fitted rate constants "must be dependent on the
solid structure … no unique reaction rate constant can be assigned to a given
ore or sinter".

The paper proposes two structural models — a **pore model** (parallel
cylindrical pores) and a **grain model** (a packing of spherical grains) — and
shows they behave almost identically. This page builds the **grain model**,
which became the canonical one: the solid is spherical grains of radius $R_g$;
gas diffuses through the interstices with an effective diffusivity $D'$; each
grain is consumed as a little shrinking core with product-layer diffusivity $D$
and first-order surface rate constant $k$.

Two honesty notes the page carries from the source. First, the paper's own
validation claim is deliberately modest: "By taking 'reasonable values' for
these parameters the model was found to reproduce the **general trends**
exhibited by the experimental data of other investigators." That is a trend
claim, not a fit, and nothing on this page presents it as a validated
prediction. Second, this Part I treats a **semi-infinite** solid with a planar
exposed face — the finite spherical pellet is Part II (*Chem. Eng. Sci.* **26**
(1971) 1901), which is not available to this gallery as a readable scan and is
not used here."""))

cells.append(md(r"""## The published model

All equations below were transcribed from renders of the scan at its native
300 ppi (CCITT-G4 bilevel; rendering higher is interpolation). The text layer of
this 1970 Pergamon scan is unusable for mathematics — it returns Eq. 21's
denominator as `R R ~ _1 P [ _ D DR$+k` and drops the mid-dot decimal points
from the figure captions — so nothing here comes from OCR. The same discipline,
and the same paper era, as `B3.1`.

**Porosity** (Eq. 20): each volume element $L'^3$ holds $N$ grains of radius
$R_g$,

$$S_g = 1 - \tfrac{4}{3}\pi \frac{R_g^3}{L'^3}N .$$

**Grain consumption** (Eq. 21, with its printed constraint): for a grain
surrounded by gas at concentration $C'_{p}$, quasi-steady diffusion through its
product shell in series with the surface reaction gives

$$-\frac{\mathrm{d}R'}{\mathrm{d}t}
 = \frac{C'_{p}}
        {\rho\left[\dfrac{R'}{D}-\dfrac{R'^2}{D R_g}+\dfrac{1}{k}\right]};
 \qquad 0 \le R' \le R_g ,$$

where $R'$ is the unreacted-core radius and $\rho$ the molar density of the
solid reactant. The initial condition is $R' = R_g$ (Eq. 22).

**Gas phase** (Eq. 27, the continuum form of the row-by-row balance Eq. 24):

$$D'\,\frac{\partial^2 C'_p}{\partial y^2}
 = \frac{3R'(y)\,(1-S_g)\,D\;C'_p(y)}
        {R_g^{\,3}\left[1-\dfrac{R'(y)}{R_g}+\dfrac{D}{R'(y)\,k}\right]} ,$$

with $C'_p = C_0$ at $y=0$ (Eq. 28) and $C'_p \to 0$ as $y\to\infty$ (Eq. 29).
Note the equation is **quasi-steady**: gas holdup in the pores is neglected
against the solid inventory. Note also that Eqs. 21 and 27 print the *same*
series resistance in two different algebraic arrangements — a transcription slip
in either one would make gas consumption and solid consumption inconsistent,
which is exactly what the mass-balance check below is able to catch.

**Extent of reaction** (Eq. 34): the *equivalent penetration*, the depth of a
fully-reacted layer containing the same amount of reacted solid,

$$\mathrm{E.P.} = \int_{y=0}^{\infty} \frac{R_g^3 - R'^3(y)}{R_g^3}\,\mathrm{d}y .$$

**Effective diffusivity** (Appendix 2): $D' = P D_p / F_T$ with the porosity
substituted for $P$ and a tortuosity factor $F_T = 2.75$.

### The dimensionless form

With
$$\xi = \frac{R'}{R_g},\qquad
\psi = \frac{C'_p}{C_0},\qquad
\eta = \frac{y}{\ell},\qquad
\theta = \frac{t}{\tau},$$

$$\ell = \sqrt{\frac{D' R_g}{3(1-S_g)k}},\qquad
\tau = \frac{\rho R_g}{k C_0},\qquad
\boxed{\;g = \frac{k R_g}{D}\;}$$

the whole model collapses to **one parameter** $g$, the ratio of product-shell
to surface-reaction resistance at the grain scale:

$$\frac{\partial^2\psi}{\partial\eta^2}
   = \frac{\xi^2\,\psi}{1+g\,\xi(1-\xi)},\qquad
  \frac{\partial \xi}{\partial\theta} = -\frac{\psi}{1+g\,\xi(1-\xi)}
  \;\;(\xi \ge 0),\qquad
  E \equiv \frac{\mathrm{E.P.}}{\ell} = \int_0^\infty (1-\xi^3)\,\mathrm{d}\eta .$$

### The exact solution the quasi-steady equations imply

Because the gas equation is quasi-steady and the kinetics first order, the model
is integrable — a consequence the source pays for but does not use. Define the
local **exposure** $\zeta(\eta,\theta) = \int_0^\theta \psi\,\mathrm{d}\theta'$.
The grain ODE integrates pointwise to

$$G(\xi) \equiv (1-\xi) + \frac{g}{6}\left(1-3\xi^2+2\xi^3\right) = \zeta ,$$

so the local core radius is a function of local exposure alone, and every grain
is fully consumed when its exposure reaches $\theta^\* = G(0) = 1+g/6$.
Integrating the gas equation in time gives a *single autonomous BVP for the
exposure field*,

$$\frac{\partial^2\zeta}{\partial\eta^2} = \frac{x(\zeta)}{3},
\qquad x(\zeta) = 1-\xi^3(\zeta),\qquad \zeta(0,\theta)=\theta,\ \ \zeta(\infty)=0,$$

whose first integral yields a **closed-form penetration law**:

$$\boxed{\;E(\theta) = \sqrt{6\,X(\theta)}\;},\qquad
X(\theta) = \int_0^\theta x(s)\,\mathrm{d}s ,$$

with $X$ itself a closed-form polynomial in the surface core radius
$\xi_\theta = G^{-1}(\theta)$ (substituting $s = G(\xi)$ turns the integral
into $\int (1-\xi^3)(1+g\xi-g\xi^2)\,\mathrm{d}\xi$). The gas profile is exact
too: $\psi = \sqrt{X(\zeta)/X(\theta)}$. Both limits fall out:

* **early / kinetic**: $E \to 3\theta$ — the exponential initial profile
  (the grain-model analogue of the paper's Eq. 18) eating every grain in place;
* **late / diffusion through the reacted layer**:
  $E \to \sqrt{6(\theta-\Delta)}$ with a lag
  $\Delta = \tfrac14 + \tfrac{g}{30}$ — the sharp-interface shrinking core of
  the *whole specimen*, which is the classical picture the paper set out to
  generalise.

This closed form is derived here from Eqs. 21 + 27; it is the second,
independent route against which the pymrm solve is judged, and it is a close
relative of the additive-reaction-time law Sohn later built on this model
(planned page `B3.4`, which cites this paper).

### The reduction to the published `B3.1`

At the exposed face $\psi \equiv 1$, so the surface grains see constant gas —
exactly `B3.1`'s single shrinking core with no film resistance. And indeed
$G(\xi)/\theta^\*$ is **term for term** Yagi & Kunii's Eq. 6 in its no-film
limit, with `B3.1`'s groups $\omega = k_{c1}/k_{d1}$, $k_{d1}=12\mathbb{D}/D_p$
evaluating to $\omega = g/6$:

$$\frac{\theta}{\theta^\*} =
\frac{(1-\xi) + \omega\,(1-3\xi^2+2\xi^3)}{1+\omega},\qquad \omega=\frac{g}{6}.$$

Two different papers, read four decades apart into two different gallery pages,
print the same function — that identity, and the simulated surface-grain history
landing on it, are checks that cannot be satisfied by accident."""))

cells.append(md(r"""## Parameters and assumptions

**The paper's assumptions** (its §2, numbered as printed): (1) no gas-film
resistance at the exposed face; (2) semi-infinite solid, macroscopically
one-dimensional; (3) first-order irreversible reaction, equimolar or dilute so
bulk flow is negligible; (4) isothermal; (5) the initial structure persists —
no sintering, swelling or pore closure. Assumption (5) is the one the authors
themselves flag as the most serious limitation. Additional to these:
quasi-steady gas phase, and the cavities as "perfect mixers" (their Eq. 23–24
picture, taken to the continuum in Eq. 27).

**What the paper prints as numbers** — all read off native-resolution crops:

* Figure 7 (grain-model base case): $D_p = 8.0$ cm²/s, $R_g = 5.0\times10^{-3}$
  cm, $D = 2.0\times10^{-4}$ cm²/s, $k = 20.0$ cm/s, porosity 25 %.
* Table 1 (grain-model ranges): $R_g$ = 1–500 μ; $S_g$ = 12.5, 25.0, 50.0 %;
  $D$ = 2×10⁻⁷ … 2×10⁻¹ cm²/s; $D_p$ = 4.0–10.0 cm²/s;
  $k$ = 0.005, 0.02, 0.2, 20.0 cm/s (the Fig. 19 sweep).
* Appendix 2: $F_T = 2.75$; Fig. 5's inset prints $\rho = 0.07$ g mol/cm³ for
  the pore-model contour computation; the text brackets practical solids as
  $0.01 < \rho < 0.1$ g mol/cm³.

**What the paper never prints: $C_0$.** Every E.P.–time result is plotted
against $t/\rho$, and the bulk gas concentration used in the computations is
stated nowhere. Since time enters the model only as $kC_0t/(\rho R_g)$, the
absolute abscissa of Figs. 6–19 cannot be reconstructed from the printed text,
and this page does not attempt it — the same refusal as `B3.1`'s missing unit
conversion. Everything runs in the model's own dimensionless groups, where the
paper's printed parameter sets translate cleanly: the four printed rate
constants of the Fig. 19 sweep are, at the base-case $R_g$ and $D$, exactly
$g = 0.125,\ 0.5,\ 5,\ 500$, and the Fig. 7 base case is $g=500$. Those four
$g$ values are used throughout this page, so every curve below is anchored to a
printed parameter set. (One deliberately *chosen* value appears in the
dimensional closure check below, and is labelled there.)

**Fitted: nothing.** No constant on this page is adjusted to make any
comparison work; the only least-squares operation is a diagnostic straight-line
fit to the page's *own computed* late-time asymptote, with a stated window."""))

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
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.integrate import quad
from scipy.optimize import brentq
from pymrm import construct_grad, construct_div, NumJac, newton
from gallery_utils import report_agreement

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})'''))

cells.append(md(r"""## The data

**There is none, and the split is precise.** What the paper *prints as numbers*
is listed above: figure-caption parameter sets, Table 1's ranges, $F_T$, and
the closed forms (Eqs. 19–34). Every computed *result* — all the E.P.–time
curves of Figs. 5–19 — exists **only as figures**, and the paper contains no
experimental measurements at all (its comparison with "the experimental data of
other investigators" is the trend statement quoted above, made against
literature it cites, not against tables it prints).

No figure is digitised for this page. Two reasons, in order: the validation
available *without* any figure — an exact closed form plus reduction to the
independently published `B3.1` — ranks higher than a digitised curve in this
gallery's hierarchy; and a figure route would require a maintainer review gate
for which no maintainer is currently available. Since $C_0$ is unprinted, a
digitised ordinate could in any case only be compared up to an unknown scale
factor on the time axis.

This page is therefore provenance **tier 6** by necessity, like `B3.1`. Nothing
below is fitted to anything; *fit vs test* reduces to "everything is test"."""))

cells.append(md(r"""## PyMRM implementation

Two independent routes to the same model.

**Route 1 — the exact closed form** derived above: a polynomial, one scalar
cubic inversion $G^{-1}$, and a square root. No grid, no time march, no linear
algebra.

**Route 2 — the paper's own formulation marched with pymrm**, as Szekely and
Evans did numerically in 1969 (their Appendix 1: a two-point BVP coupled to the
interface ODEs). Here:

* the gas BVP is discretised with `construct_grad` / `construct_div` on a
  finite-volume grid, with **`nu=0`** — the sample scale is a Cartesian
  semi-infinite slab. (`nu=2` would be a *spherical pellet*, which is Part II's
  finite-pellet geometry, not this paper's. The grains are spherical, but their
  sphericity enters through the shrinking-core closure, not through a pymrm
  operator.) The semi-infinite domain is truncated at $\eta = L$ with $L$ far
  beyond the reaction zone; the far tail decays like $e^{-\eta}$, and the
  truncation insensitivity is demonstrated in the break table.
* boundary conditions use the **outward normal** ($a\,\partial\psi/\partial n +
  b\,\psi = d$): at $\eta=0$, `{"a": 0, "b": 1, "d": 1}` imposes $\psi = 1$; at
  $\eta=L$, `{"a": 0, "b": 1, "d": 0}` imposes $\psi = 0$. Both are Dirichlet,
  so the sign of $a$ is moot here — the comment states the physical equation
  anyway, per house rule.
* the grain field $\xi$ and the gas field $\psi$ are solved **coupled** at each
  time level with `newton`, the Jacobian assembled as (analytic sparse Laplacian)
  + (`NumJac((n, 2))` for the pointwise source terms). The state layout is
  spatial axis first, fields last, `(n, 2)` — and note the house trap: a bare
  1-D shape `(n,)` would declare the *spatial* axis fully coupled and build a
  dense Jacobian.
* Eq. 21's printed constraint $0 \le R' \le R_g$ — grains stop existing when
  consumed — is imposed as a complementarity condition,
  $\min(\text{trapezoidal residual},\ \xi) = 0$: either the ODE holds, or the
  grain is exhausted and $\xi = 0$. The time scheme is the trapezoidal rule;
  its observed order is *measured* below rather than asserted, because the
  exhaustion clamp is a kink the classical order theory does not cover."""))

cells.append(code(r'''# ---------------------------------------------------------------- exact route
def G_of_xi(xi, g):
    """Grain exposure-conversion law: G(xi) = zeta. Identical, term for term,
    to Yagi & Kunii eq. (6) (page B3.1) with no film and omega = g/6."""
    return (1.0 - xi) + (g / 6.0) * (1.0 - 3.0 * xi**2 + 2.0 * xi**3)


def _P(xi, g):
    """Antiderivative of (1 - xi^3)(1 + g xi - g xi^2)."""
    return (xi + g * xi**2 / 2 - g * xi**3 / 3 - xi**4 / 4 - g * xi**5 / 5
            + g * xi**6 / 6)


def xi_of_zeta(zeta, g):
    """Invert G on [0, 1] (monotone). Scalar."""
    ts = 1.0 + g / 6.0
    if zeta <= 0.0:
        return 1.0
    if zeta >= ts:
        return 0.0
    return brentq(lambda x: G_of_xi(x, g) - zeta, 0.0, 1.0,
                  xtol=1e-15, rtol=8.9e-16)


def X_of_theta(theta, g):
    """X(theta) = integral_0^theta x(s) ds, closed form."""
    ts = 1.0 + g / 6.0
    if theta >= ts:
        return _P(1.0, g) + (theta - ts)          # P(0) = 0
    # max() guards the P(1) - P(xi) cancellation at theta -> 0 (roundoff can
    # return a tiny negative); callers stay above theta ~ 1e-5, where the
    # relative cancellation error is ~1e-6 of X and negligible for every use.
    return max(_P(1.0, g) - _P(xi_of_zeta(theta, g), g), 0.0)


def E_exact(theta, g):
    """Equivalent penetration E.P./ell, exact: E = sqrt(6 X(theta))."""
    th = np.atleast_1d(np.asarray(theta, float))
    return np.sqrt(6.0 * np.array([X_of_theta(t, g) for t in th]))


def eta_of_zeta(zeta, theta, g):
    """Exact spatial coordinate of the exposure level zeta at time theta."""
    v, _ = quad(lambda s: 1.0 / np.sqrt((2.0 / 3.0) * X_of_theta(s, g)),
                zeta, theta, limit=200)
    return v


def yagi_kunii_eq6(x, omega, gamma):
    """B3.1's published equation, transcribed from the Yagi & Kunii page:
    theta/theta_B against x = r/R, groups omega = kc1/kd1, gamma = 3 kf1/kd1."""
    x = np.asarray(x, float)
    num = gamma * x + 3.0 * omega * gamma * x**2 + omega * (1.0 - 2.0 * gamma) * x**3
    return 1.0 - num / (omega + omega * gamma + gamma)


# --- the printed anchors (all read off native 300 ppi crops) -----------------
FIG7 = dict(Dp=8.0, Rg=5.0e-3, D=2.0e-4, k=20.0, Sg=0.25, FT=2.75)  # + Table 1
K_SWEEP = [0.005, 0.02, 0.2, 20.0]           # Fig. 19 / Table 1, cm/s
G_SWEEP = [kk * FIG7["Rg"] / FIG7["D"] for kk in K_SWEEP]

Deff = FIG7["Sg"] * FIG7["Dp"] / FIG7["FT"]                    # D' = Sg Dp / FT
ell = np.sqrt(Deff * FIG7["Rg"] / (3 * (1 - FIG7["Sg"]) * FIG7["k"]))
print("printed Fig. 7 base case ->  D' = Sg*Dp/FT = %.4f cm^2/s" % Deff)
print("length scale ell = sqrt(D' Rg / (3(1-Sg)k)) = %.4e cm" % ell)
print("g = k Rg / D for the printed k-sweep {0.005, 0.02, 0.2, 20} cm/s:",
      [f"{gg:g}" for gg in G_SWEEP])
print("theta* = 1 + g/6 (surface burnout):",
      [f"{1 + gg / 6:g}" for gg in G_SWEEP])

# --- identity: G(xi)/theta* IS B3.1's eq. 6 with no film, omega = g/6 --------
xi_t = np.linspace(0.0, 1.0, 1001)
worst_id = 0.0
for gg in G_SWEEP:
    om = gg / 6.0
    lhs = G_of_xi(xi_t, gg) / (1.0 + om)
    rhs_analytic = ((1 - xi_t) + om * (1 - 3 * xi_t**2 + 2 * xi_t**3)) / (1 + om)
    rhs_b31 = yagi_kunii_eq6(xi_t, om, 1e12)      # published form, film -> 0
    worst_id = max(worst_id, np.abs(lhs - rhs_analytic).max(),
                   np.abs(lhs - rhs_b31).max())
print(f"\nG(xi)/theta* vs B3.1 eq. 6 (no film, omega = g/6): "
      f"max |diff| = {worst_id:.2e} over the four printed g")
print("   -> STRUCTURAL: an algebraic identity between two page transcriptions.")
print("      It proves the two papers print the same function (and guards the")
print("      transcription), but it cannot fail once both are typed correctly;")
print("      it also sits below the CI comparison floor (ABS_FLOOR = 1e-12).")
print("      CI skips only while BOTH baseline and current sit under the floor,")
print("      so a regression LIFTING it above 1e-12 would fail CI - what is")
print("      unprotected is drift beneath the floor. (The 4.4e-13 itself is")
print("      finite-gamma truncation at gamma = 1e12, not roundoff.) The check with power")
print("      is the SIMULATED surface-grain history against it, further down.")'''))

cells.append(code(r'''# ------------------------------------------------- route 2: the pymrm march
class GrainModel:
    """Szekely-Evans grain model, semi-infinite 1-D, quasi-steady gas phase.

        diff * psi'' = sink_pref * xi^2 psi / (1 + g xi(1-xi))     (Eq. 27)
        dxi/dt      = -rate_pref * psi / (1 + g xi(1-xi)), xi >= 0 (Eq. 21)

    Dimensionless: all prefactors 1. Dimensional (cgs): diff = D',
    sink_pref = 3(1-Sg)k/Rg, rate_pref = k/(rho Rg), c_surf = C0.
    `defect` deliberately injects known errors for the break table.
    """

    def __init__(self, n=600, L=35.0, g=0.5, diff=1.0, sink_pref=1.0,
                 rate_pref=1.0, c_surf=1.0, defect=None):
        self.n, self.g, self.defect = n, g, defect
        self.diff, self.sink_pref, self.rate_pref = diff, sink_pref, rate_pref
        self.x_f = np.linspace(0.0, L, n + 1)
        self.x_c = 0.5 * (self.x_f[:-1] + self.x_f[1:])
        self.dx = np.diff(self.x_f)
        shape = (n, 1)
        # Outward-normal BCs, a*dpsi/dn + b*psi = d:
        #   eta = 0:  0*dpsi/dn + 1*psi = c_surf   -> psi(0) = c_surf  (Eq. 28)
        #   eta = L:  0*dpsi/dn + 1*psi = 0        -> psi(L) = 0       (Eq. 29,
        #             the truncated form of psi -> 0 as y -> infinity)
        d0 = 0.9 * c_surf if defect == "surface_conc_0p9" else c_surf
        bc_R = ({"a": 1.0, "b": 0.0, "d": 0.0} if defect == "bc_far_neumann"
                else {"a": 0.0, "b": 1.0, "d": 0.0})
        bc = ({"a": 0.0, "b": 1.0, "d": d0}, bc_R)
        grad, grad_bc = construct_grad(shape, self.x_f, self.x_c, bc)
        # nu = 0: Cartesian slab (semi-infinite SAMPLE; the grains' sphericity
        # lives in the closure, not in the divergence operator). nu = 2 would
        # be Part II's spherical pellet and is one of the injected defects.
        nu = 2 if defect == "nu2" else 0
        div = construct_div(shape, self.x_f, nu=nu)
        self.grad, self.grad_bc = grad, grad_bc
        self.Lop = (div @ (-diff * grad)).tocsr()      # -(diff * psi')' term
        self.bvec = np.asarray((div @ (-diff * grad_bc)).todense()).ravel()
        # Embed the psi-Laplacian in the interleaved (n,2) state [psi, xi]:
        N2 = 2 * n
        Lc = self.Lop.tocoo()
        self.A = sp.coo_matrix((Lc.data, (2 * Lc.row, 2 * Lc.col)),
                               shape=(N2, N2)).tocsc()
        self.b = np.zeros(N2)
        self.b[0::2] = self.bvec
        # Pointwise source: NumJac((n, 2)) couples the LAST axis (the two
        # fields) per cell. A bare (n,) here would couple the spatial axis in
        # full and build a dense n x n Jacobian - the trap measured on B1.1.
        self.numjac = NumJac((n, 2))
        self.t = 0.0
        self.xi = np.ones(n)
        self.psi = self.solve_psi_linear(self.xi)
        self.hist = []
        self.snapshots = {}
        self.I_gas = 0.0
        self.q_prev = self.influx()
        self.newton_iters = []

    # --- closure terms ------------------------------------------------------
    def _denom(self, xi):
        xp = np.clip(xi, 0.0, 1.0)
        if self.defect == "denom_no_shield":
            return 1.0 + self.g * xp                 # WRONG: drops (1 - xi)
        return 1.0 + self.g * xp * (1.0 - xp)

    def sink_coeff(self, xi):
        xp = np.clip(xi, 0.0, 1.0)
        area = xp**3 if self.defect == "sink_area_cubed" else xp**2
        return self.sink_pref * area / self._denom(xi)

    def rate(self, psi, xi):
        """-dxi/dt, unclamped (the NCP supplies the clamp)."""
        if self.defect == "rate_denom_mismatch":
            xp = np.clip(xi, 0.0, 1.0)
            return self.rate_pref * psi / (1.0 + self.g * xp)  # WRONG vs sink
        return self.rate_pref * psi / self._denom(xi)

    # --- residual for newton ------------------------------------------------
    def _source(self, u):
        psi, xi = u[:, 0], u[:, 1]
        r = np.empty_like(u)
        r[:, 0] = self.sink_coeff(xi) * psi
        be = (xi - self.xi_old
              + 0.5 * self.dt * (self.rate(psi, xi) + self.r_old))
        r[:, 1] = np.minimum(be, xi)     # complementarity: Eq. 21's 0 <= xi
        return r

    def _residual(self, u):
        s, js = self.numjac(self._source, u)
        F = self.A @ u.reshape(-1, 1) + self.b.reshape(-1, 1) + s.reshape(-1, 1)
        return F, self.A + js

    # --- solves -------------------------------------------------------------
    def solve_psi_linear(self, xi):
        """Quasi-steady gas profile for a frozen grain field (linear)."""
        Amat = (self.Lop + sp.diags(self.sink_coeff(xi))).tocsc()
        return sp.linalg.spsolve(Amat, -self.bvec)

    def influx(self):
        """Boundary flux from the scheme's own gradient operator."""
        f = -self.diff * (self.grad @ self.psi.reshape(-1, 1)
                          + self.grad_bc.todense())
        return float(np.asarray(f).ravel()[0])

    def march(self, t_end, dt, store_every=1, snap_at=()):
        self.dt = dt
        snap_at = sorted(snap_at)
        u = np.column_stack([self.psi, self.xi])
        nst = int(round(t_end / dt))
        for i in range(nst):
            self.xi_old = self.xi.copy()
            self.r_old = np.where(self.xi_old > 0.0,
                                  self.rate(self.psi, self.xi_old), 0.0)
            res = newton(self._residual, u, tol=1e-12, maxfev=30)
            assert res.success, f"newton failed at step {i} (t={self.t:.3f})"
            self.newton_iters.append(res.nit)
            u = res.x.reshape(self.n, 2)
            self.psi = u[:, 0].copy()
            self.xi = np.maximum(u[:, 1], 0.0)
            u[:, 1] = self.xi
            self.t += dt
            q = self.influx()
            self.I_gas += 0.5 * dt * (self.q_prev + q)   # trapezoid in time
            self.q_prev = q
            if (i + 1) % store_every == 0:
                self.hist.append((self.t, self.EP(), self.xi_surface(),
                                  self.front()))
            while snap_at and self.t >= snap_at[0] - 0.5 * dt:
                self.snapshots[snap_at.pop(0)] = (self.psi.copy(),
                                                  self.xi.copy())
        return self

    # --- observables --------------------------------------------------------
    def EP(self):
        """Eq. 34, the equivalent penetration."""
        return float(np.sum((1.0 - np.clip(self.xi, 0, 1)**3) * self.dx))

    def xi_surface(self):
        """xi extrapolated to eta = 0 (linear, from the first two centres)."""
        x0, x1 = self.x_c[0], self.x_c[1]
        return float(self.xi[0] + (self.xi[0] - self.xi[1]) * x0 / (x1 - x0))

    def front(self):
        """Depth of the completely-reacted zone (paper: 'zone of complete
        reduction'): the outermost face below which every grain is consumed."""
        alive = np.nonzero(self.xi > 0.0)[0]
        return float(self.x_f[alive[0]] if alive.size else self.x_f[-1])

    def solid_consumed(self):
        return self.EP() * self.sink_pref / (3.0 * self.rate_pref)


print("GrainModel ready (route 2)")'''))

cells.append(md(r"""## Results

The four printed $g$ values of the Fig. 19 sweep, run both ways: exact closed
form (lines) and the pymrm march (markers, two of the four). Dashed envelopes
are the two classical limits the grain model interpolates — the kinetic law
$E = 3\theta$ and the specimen-scale shrinking core
$E = \sqrt{6(\theta-\Delta)}$."""))

cells.append(code(r'''runs = {}
for gg in (0.5, 5.0):
    runs[gg] = GrainModel(n=600, L=35.0, g=gg).march(
        20.0, 0.01, store_every=10,
        snap_at=(0.25, 1.0, 2.0) if gg == 0.5 else ())

theta_plot = np.logspace(-2.3, np.log10(20.0), 220)
COL = {0.125: "tab:blue", 0.5: "tab:orange", 5.0: "tab:green", 500.0: "tab:red"}

fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.5))
for gg in G_SWEEP:
    ax[0].loglog(theta_plot, E_exact(theta_plot, gg), lw=1.9, color=COL[gg],
                 label=f"exact, $g$ = {gg:g}")
for gg, m in runs.items():
    th = np.array([h[0] for h in m.hist])
    Ep = np.array([h[1] for h in m.hist])
    ax[0].loglog(th[::10], Ep[::10], "o", ms=4, mfc="none", color=COL[gg],
                 label=f"pymrm, $g$ = {gg:g}")
ax[0].loglog(theta_plot, 3 * theta_plot, "k--", lw=1.1, alpha=0.6,
             label=r"kinetic limit $E=3\theta$ (slope 1)")
th_l = np.linspace(0.35, 20, 100)
ax[0].loglog(th_l, np.sqrt(np.maximum(6 * (th_l - 0.25 - 0.5 / 30), 0)), "k:",
             lw=1.6, alpha=0.9,
             label=r"diffusion limit, $g$ = 0.5 (slope $\frac{1}{2}$)")
ax[0].set(xlabel=r"$\theta = kC_0t/(\rho R_g)$",
          ylabel=r"$E = \mathrm{E.P.}/\ell$", ylim=(1e-2, 30),
          title="the printed Fig. 19 sweep, dimensionless:\n"
                r"$g = kR_g/D$ = 0.125, 0.5, 5, 500")
ax[0].legend(fontsize=8, loc="upper left")
ax[0].text(1.5, 0.09, "the $g \\leq 5$ curves nearly collapse:\n"
           "in these units the structure lives\n"
           "almost entirely in $\\ell$ and $\\tau$;\n"
           "$g = 500$ alone lags visibly",
           fontsize=8, ha="left")

for gg, m in runs.items():
    th = np.array([h[0] for h in m.hist])
    Ep = np.array([h[1] for h in m.hist])
    rel = np.abs(Ep - E_exact(th, gg)) / E_exact(th, gg)
    ax[1].semilogy(th, rel, lw=1.6, color=COL[gg], label=f"$g$ = {gg:g}")
ax[1].set(xlabel=r"$\theta$", ylabel="|pymrm $-$ exact| / exact",
          title="route 2 against route 1\n(n = 600, $\\Delta\\theta$ = 0.01)")
ax[1].legend(fontsize=9)
fig.tight_layout(); plt.show()

hl = {}
for gg, m in runs.items():
    th = np.array([h[0] for h in m.hist])
    Ep = np.array([h[1] for h in m.hist])
    sel = th >= 0.1
    hl[gg] = float(np.max(np.abs(Ep - E_exact(th, gg))[sel]
                          / E_exact(th, gg)[sel]))
    it = np.array(runs[gg].newton_iters)
    print(f"g = {gg:g}: max rel |E_pymrm - E_exact| over 0.1 <= theta <= 20 "
          f"= {hl[gg]:.2e}   (newton iterations/step {it.min()}-{it.max()})")
print("\nThe two routes share no code: one is a finite-volume march through")
print("(psi, xi) with a Newton solve per step, the other a polynomial and one")
print("scalar root find. They DO share the model statement itself - what that")
print("agreement cannot detect is discussed under Validation.")'''))

cells.append(code(r'''# ------- profiles: pymrm markers against the exact parametric solution ------
g05 = 0.5
m = runs[g05]
fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.3))
prof_err_psi, prof_err_xi = 0.0, 0.0
for i, (th, colr) in enumerate(((0.25, "tab:blue"), (1.0, "tab:orange"),
                                (2.0, "tab:red"))):
    psi_n, xi_n = m.snapshots[th]
    zs = np.linspace(min(th, 1e-2) * 1e-3, th, 300)[1:]
    etas = np.array([eta_of_zeta(z, th, g05) for z in zs])
    psis = np.sqrt([X_of_theta(z, g05) / X_of_theta(th, g05) for z in zs])
    xis = np.array([xi_of_zeta(z, g05) for z in zs])
    show = m.x_c <= 6.6
    ax[0].plot(etas, psis, lw=1.8, color=colr, label=rf"$\theta$ = {th:g}")
    ax[0].plot(m.x_c[show][::12], psi_n[show][::12], "o", ms=3.5, mfc="none",
               color=colr)
    ax[1].plot(etas, 1 - xis**3, lw=1.8, color=colr,
               label=rf"$\theta$ = {th:g}")
    ax[1].plot(m.x_c[show][::12], 1 - np.clip(xi_n[show][::12], 0, 1)**3, "o",
               ms=3.5, mfc="none", color=colr)
    sel = etas <= 12
    prof_err_psi = max(prof_err_psi, np.abs(
        np.interp(etas[sel], m.x_c, psi_n) - psis[sel]).max())
    prof_err_xi = max(prof_err_xi, np.abs(
        np.interp(etas[sel], m.x_c, xi_n) - xis[sel]).max())
ax[0].set(xlabel=r"$\eta$", ylabel=r"$\psi = C_p'/C_0$", xlim=(0, 8),
          title="gas profile (lines exact, circles pymrm), $g$ = 0.5")
ax[1].set(xlabel=r"$\eta$", ylabel=r"local conversion $x = 1-\xi^3$",
          xlim=(0, 8),
          title="solid conversion: the DIFFUSE zone\nthe shrinking-core "
                "picture replaces with a front")
for a in ax:
    a.legend(fontsize=9)
fig.tight_layout(); plt.show()
print(f"profile agreement (0 <= eta <= 12): max|psi| dev {prof_err_psi:.1e}, "
      f"max|xi| dev {prof_err_xi:.1e}")
print("(limited by linear interpolation across the corner at the")
print(" complete-conversion front, not by the integral quantities; E itself")
print(" converges at second order, measured under Validation)")'''))

cells.append(code(r'''# ------- the diffuse zone, quantified ---------------------------------------
# The paper's central qualitative claim: reaction occurs in a zone whose
# thickness "may be of the same order as the size of the reacting specimen".
# The exact solution makes that a number: the zone spans exposures
# 0 < zeta < theta*, and its width (between local conversion 5% and 95%)
# tends to a CONSTANT w(g) while the penetration keeps growing like sqrt(6 theta).
def zone_width(gg, theta):
    z_lo = G_of_xi((1 - 0.05)**(1 / 3), gg)     # x = 5 %
    z_hi = G_of_xi((1 - 0.95)**(1 / 3), gg)     # x = 95 %
    return eta_of_zeta(z_lo, theta, gg) - eta_of_zeta(z_hi, theta, gg)

print("g        theta*      w(4 theta*)   w(8 theta*)   [zone width, units of ell]")
wg = {}
for gg in G_SWEEP:
    ts = 1 + gg / 6
    w4, w8 = zone_width(gg, 4 * ts), zone_width(gg, 8 * ts)
    wg[gg] = w8
    print(f"{gg:7.3f}  {ts:9.3f}   {w4:10.4f}   {w8:10.4f}")
print("\nThe width is already time-independent at 4 theta* - a travelling")
print("reaction zone of fixed thickness behind a sqrt-in-time front.")

w500_cm = wg[500.0] * ell
print(f"\nFor the printed Fig. 7 base case (g = 500):")
print(f"  zone width = {wg[500.0]:.1f} ell = {w500_cm:.3f} cm,")
print(f"  using only printed constants (Dp, Rg, D, k, Sg, FT) - no C0 enters,")
print(f"  because ell and w are both C0-free. Szekely & Evans' diffuse-zone")
print(f"  claim, made quantitative: a ~{w500_cm*10:.0f} mm reaction zone in a")
print(f"  specimen whose penetration Fig. 7 plots on a half-centimetre axis.")

# ------- E.P. and the zone of complete reaction (Fig. 7's two curves) -------
fig, ax = plt.subplots(figsize=(7.4, 4.4))
for gg, m in runs.items():
    th = np.array([h[0] for h in m.hist])
    Ep = np.array([h[1] for h in m.hist])
    Yf = np.array([h[3] for h in m.hist])
    ax.plot(th, Ep, lw=1.9, color=COL[gg], label=f"E.P., $g$ = {gg:g}")
    ax.plot(th, Yf, lw=1.4, ls="--", color=COL[gg],
            label=f"complete-reaction zone, $g$ = {gg:g}")
ax.set(xlabel=r"$\theta$", ylabel=r"depth / $\ell$",
       title="the two curves Fig. 7 draws, dimensionless:\nequivalent "
             "penetration vs zone of complete reaction (pymrm)")
ax.legend(fontsize=9)
fig.tight_layout(); plt.show()
print("The gap between each pair of curves is the diffuse zone's inventory;")
print("it tends to the constant lag Delta = 1/4 + g/30 (checked below), which")
print("is why a shrinking-core fit to grain-model data returns a rate constant")
print("contaminated by structure - the paper's opening complaint.")

# ------- where the model departs from each classical limit ------------------
print("\nWhere the grain model departs from each limit (exact route, 5% bands):")
print("g        kinetic E=3theta good to   sqrt-law within 5% from   window ratio")
dep = {}
for gg in G_SWEEP:
    Delta = 0.25 + gg / 30
    th_e = brentq(lambda t: E_exact(t, gg)[0] / (3 * t) - 0.95, 1e-5, 1e5,
                  xtol=1e-12)
    th_l = brentq(lambda t: np.sqrt(6 * (t - Delta)) / E_exact(t, gg)[0] - 0.95,
                  Delta * (1 + 1e-12), 1e7, xtol=1e-10)
    dep[gg] = (th_e, th_l)
    print(f"{gg:7.3f}   theta < {th_e:8.4f}          theta > {th_l:9.3f}"
          f"          {th_l / th_e:9.0f}x")
print("\nSmall g: the two classical laws nearly meet - the grain model is only")
print("needed in a factor-~3 window of theta. Large g (the printed Fig. 7 base")
print("case): the window spans more than four orders of magnitude in time, and")
print("neither textbook law describes most of the conversion history.")'''))

cells.append(md(r"""## Validation

Six checks. Each one names what it can and cannot catch, every reported metric
appears in the break table below, and everything is deterministic — fixed
grids, direct solves, bracketed root finds, no continuation.

1. **Route 2 against route 1** (above): the finite-volume march against the
   closed form, over the whole conversion history, two values of $g$. This is
   the headline. What it *cannot* catch: a mis-transcription of the model
   statement itself, which both routes inherit — that risk is covered by check
   2 and by the native-render transcription protocol.
2. **Reduction to the published `B3.1`**: the *simulated* surface-grain
   history $\xi(0,\theta)$ against Yagi & Kunii's Eq. 6 as published on the
   `B3.1` page (no film, $\omega = g/6$) — a different paper, transcribed
   independently, validated there to 6.9e-16 against its own integration. The
   algebraic identity between the two printed functions is separately recorded
   (and labelled structural); *this* check exercises the simulation.
3. **The initial gas profile** against the exact $\sinh$ solution on the
   truncated domain — the grain-model analogue of the paper's printed Eq. 18.
   Blind, by construction, to everything about the march (it is $\theta = 0$)
   and to the area exponent (at $\xi = 1$, $\xi^2 = \xi^3$): the break table
   shows exactly that blindness.
4. **Both axes refined, orders observed**: grid at fixed $\Delta\theta$
   against the exact solution; $\Delta\theta$ at fixed grid against a
   same-grid small-step reference (so the temporal order is not hidden under
   the spatial error floor). The step is uniform — no geometric ramp — so the
   knob refined is the knob measured.
5. **Mass balance**, split honestly into its two parts: the *spatial* identity
   (boundary influx = cell-sink sum at fixed time) telescopes by construction
   of `construct_div` and is reported as structural; the *temporal* closure
   (time-integrated influx vs solid consumed) is a real check of the
   march + clamp consistency, converges at the scheme's order, and is the one
   check that catches the paper-specific hazard of Eqs. 21/27 printing one
   resistance in two algebraic forms.
6. **Dimensional closure**: the same class run in cgs units on the printed
   Fig. 7 base case, against $\ell\,E_\mathrm{exact}(t/\tau)$ — testing the
   $\ell,\tau,g$ reduction algebra end to end. $C_0$ is not printed in the
   paper, so a **chosen** illustrative $C_0 = 1.0\times10^{-5}$ g mol/cm³
   (ideal gas at roughly ambient pressure and furnace temperature) is used and
   labelled; the check's verdict is independent of that choice, because both
   routes receive it."""))

cells.append(code(r'''# ---- check 2: simulated surface grain vs B3.1's published closed form ------
print("2. surface-grain history vs B3.1 eq. 6 (no film, omega = g/6)")
surf_err = {}
for n in (150, 300, 600):
    mm = GrainModel(n=n, L=30.0, g=0.5).march(1.0, 0.0025, store_every=40)
    w = 0.0
    for th, _, xs, _ in mm.hist:
        if xs > 1e-3:
            om = 0.5 / 6
            th_pred = yagi_kunii_eq6(np.clip(xs, 0, 1), om, 1e12) * (1 + om)
            w = max(w, abs(th_pred - th) / (1 + om))
    surf_err[n] = w
    o = (f"   order = {np.log2(surf_err[n // 2] / w):.2f}" if n > 150 else "")
    print(f"   n = {n:4d}   max |theta(eq.6) - theta| / theta* = {w:.3e}{o}")
print("   second-order convergence onto ANOTHER PAPER's closed form - the")
print("   grain model reduces to the published shrinking core at the exposed")
print("   face, as it must. (xi is extrapolated to eta = 0; the residual error")
print("   is the extrapolation + march error, which is what refines away.)\n")

# ---- check 3: initial profile vs exact sinh --------------------------------
print("3. theta = 0 gas profile vs exact sinh((L-eta))/sinh(L)")
init_err = {}
for n in (150, 300, 600):
    mm = GrainModel(n=n, L=35.0, g=0.5)
    exact0 = np.sinh(35.0 - mm.x_c) / np.sinh(35.0)
    init_err[n] = float(np.abs(mm.psi - exact0).max())
    o = (f"   order = {np.log2(init_err[n // 2] / init_err[n]):.2f}"
         if n > 150 else "")
    print(f"   n = {n:4d}   max |psi - exact| = {init_err[n]:.3e}{o}")
print("   (the paper prints the pore-model twin of this profile as its Eq. 18)")'''))

cells.append(code(r'''# ---- check 4: both axes, observed orders -----------------------------------
print("4a. grid refinement, Delta_theta = 2e-3 fixed, theta_end = 2, vs exact")
Eex2 = E_exact(2.0, 0.5)[0]
grid_err = []
for n in (75, 150, 300, 600):
    mm = GrainModel(n=n, L=30.0, g=0.5).march(2.0, 2e-3)
    e = abs(mm.EP() - Eex2) / Eex2
    grid_err.append((n, e))
    o = (f"   order = {np.log2(grid_err[-2][1] / e):.2f}"
         if len(grid_err) > 1 else "")
    print(f"   n = {n:4d}   rel err = {e:.3e}{o}")
grid_order = float(np.log2(grid_err[-2][1] / grid_err[-1][1]))

print("\n4b. time-step refinement, n = 300 fixed, vs same-grid reference at")
print("    Delta_theta = 6.25e-4 (isolates the temporal error; against the")
print("    exact solution the spatial floor ~4e-5 would mask it)")
ref = GrainModel(n=300, L=30.0, g=0.5).march(2.0, 6.25e-4).EP()
dt_err = []
for dt in (0.08, 0.04, 0.02, 0.01, 0.005):
    mm = GrainModel(n=300, L=30.0, g=0.5).march(2.0, dt)
    e = abs(mm.EP() - ref) / ref
    dt_err.append((dt, e))
    o = (f"   order = {np.log2(dt_err[-2][1] / e):.2f}"
         if len(dt_err) > 1 else "")
    print(f"   dt = {dt:.3f}   rel err = {e:.3e}{o}")
dt_order = float(np.log2(dt_err[-2][1] / dt_err[-1][1]))
print("   theta_end = 2 crosses the surface burnout theta* = 1.083, so these")
print("   orders INCLUDE the complementarity clamp firing - the trapezoidal")
print("   rule's second order survives it, measured, not assumed.")'''))

cells.append(code(r'''# ---- check 5: mass balance, split into its two parts -----------------------
print("5a. spatial identity at fixed time (STRUCTURAL - telescoping sum):")
mm = GrainModel(n=300, L=30.0, g=0.5)
f = -(mm.grad @ mm.psi.reshape(-1, 1) + mm.grad_bc.todense())
f = np.asarray(f).ravel()
sink_sum = float(np.sum(mm.sink_coeff(mm.xi) * mm.psi * mm.dx))
spatial_gap = abs((f[0] - f[-1]) - sink_sum)
print(f"   (influx - outflux) - sum(sink * dx) = {spatial_gap:.2e}")
print("   This CANNOT fail while the operators come from construct_grad /")
print("   construct_div - the divergence telescopes. It is bookkeeping, kept")
print("   because it proves the flux post-processing reads the same operator")
print("   the solve uses; it detects nothing about the physics, and it is")
print("   below the CI floor. The check with power is 5b.\n")

print("5b. temporal closure: integral of boundary influx vs solid consumed")
mb = []
for dt in (0.04, 0.02, 0.01, 0.005):
    mm = GrainModel(n=300, L=30.0, g=0.5).march(2.0, dt)
    e = abs(mm.I_gas - mm.solid_consumed()) / mm.solid_consumed()
    mb.append((dt, e))
    o = f"   order = {np.log2(mb[-2][1] / e):.2f}" if len(mb) > 1 else ""
    print(f"   dt = {dt:.3f}   rel gap = {e:.3e}{o}")
mass_bal = mb[-1][1]
print("   Converges at the scheme's order -> the march consumes gas and solid")
print("   consistently, INCLUDING through grain-exhaustion events. This is the")
print("   check that would catch transcribing Eq. 21 and Eq. 27 with different")
print("   resistance denominators (the two equations print the same physics in")
print("   two algebraic forms) - demonstrated in the break table.\n")

# ---- check 6: dimensional closure on the printed Fig. 7 base case ----------
print("6. dimensional run, printed Fig. 7 constants, CHOSEN C0 = 1.0e-5 mol/cm3")
C0 = 1.0e-5              # chosen, NOT printed in the paper; labelled above
rho = 0.07               # g mol/cm3, the value printed in Fig. 5's inset
g500 = FIG7["k"] * FIG7["Rg"] / FIG7["D"]
tau = rho * FIG7["Rg"] / (FIG7["k"] * C0)
# Domain: at g = 500 the exposure profile at t = 5 tau reaches eta ~ 20
# (exact-route quadrature), so L = 30 ell. A first draft of this cell used
# L = 6 ell and THIS CHECK RETURNED 48% - domain truncation, caught here
# precisely because the far-field hazard is invisible to interior residuals
# (see the break table's far-boundary row). The check can fail, and did.
# Step: at g = 500 the start-up has a fast internal timescale ~ tau/g (the
# shell resistance 1 + g*xi*(1-xi) grows ~g-fold over the first instants), so
# the step must be finer than for the g = O(1) runs. Swept below so the
# claim is measured, not asserted.
dim_sweep = []
for dtf in (0.008, 0.004, 0.002):
    mdim = GrainModel(n=600, L=30.0 * ell, g=g500, diff=Deff,
                      sink_pref=3 * (1 - FIG7["Sg"]) * FIG7["k"] / FIG7["Rg"],
                      rate_pref=FIG7["k"] / (rho * FIG7["Rg"]), c_surf=C0)
    mdim.march(5.0 * tau, dtf * tau, store_every=max(1, int(0.1 / dtf)))
    t_d = np.array([h[0] for h in mdim.hist])
    EP_d = np.array([h[1] for h in mdim.hist])
    EP_pred = ell * E_exact(t_d / tau, g500)
    sel = t_d / tau >= 0.1
    dim_sweep.append((dtf, float(np.max(np.abs(EP_d - EP_pred)[sel]
                                        / EP_pred[sel]))))
dim_err = dim_sweep[-1][1]
print(f"   ell = {ell:.4e} cm, tau = {tau:.3f} s, g = {g500:g}")
for i, (dtf, e) in enumerate(dim_sweep):
    o = (f"   order = {np.log2(dim_sweep[i-1][1] / e):.2f}" if i else "")
    print(f"   dt = {dtf:.3f} tau   max rel |E.P. - ell*E_exact(t/tau)| "
          f"over 0.1 <= t/tau <= 5 = {e:.2e}{o}")
print("   (the maximum sits in the early fast transient and refines away -")
print("    order ~2 on the finest pair; the first pair is pre-asymptotic and")
print("    also compares maxima on different stored time grids. The error at")
print("    t = 5 tau itself is an order smaller than the reported maximum)")
print("   The march never sees ell, tau or theta - it runs in cm and seconds -")
print("   so this closes the reduction algebra end to end. The verdict is")
print("   independent of the chosen C0 (both routes receive it). It is also a")
print("   check that has already fired once: a first draft truncated the")
print("   domain at 6 ell and this comparison returned 48%, the far-field")
print("   hazard no interior residual can see (see the code comment).\n")

# ---- late-time lag: the fitted asymptote vs Delta = 1/4 + g/30 -------------
print("late-time asymptote: fit of (theta - E^2/6) on 12 <= theta <= 20")
lag_err = {}
for gg, m in runs.items():
    th = np.array([h[0] for h in m.hist])
    Ep = np.array([h[1] for h in m.hist])
    sel = th >= 12.0
    lag_fit = float(np.mean(th[sel] - Ep[sel]**2 / 6.0))
    Delta = 0.25 + gg / 30
    lag_err[gg] = abs(lag_fit - Delta) / Delta
    print(f"   g = {gg:g}: fitted lag {lag_fit:.4f} vs exact 1/4 + g/30 = "
          f"{Delta:.4f}   rel {lag_err[gg]:.1e}")
print("   (a diagnostic fit to this page's own computed curve - stated window,")
print("    nothing from the paper is fitted)")'''))

cells.append(code(r'''# ---- the break table -------------------------------------------------------
# Every reported metric, recomputed under every injected defect (short config:
# n = 300, L = 30, theta_end = 2, dt = 0.01, g = 0.5). "no move" rows are the
# finding, not a failure - each is explained below the table.
def metrics_under(defect):
    out = {}
    mm = GrainModel(n=300, L=30.0, g=0.5, defect=defect)
    exact0 = np.sinh(30.0 - mm.x_c) / np.sinh(30.0)
    out["init_profile"] = float(np.abs(mm.psi - exact0).max())
    mm.march(2.0, 0.01, store_every=5)
    th = np.array([h[0] for h in mm.hist])
    Ep = np.array([h[1] for h in mm.hist])
    sel = th >= 0.1
    Ee = E_exact(th, 0.5)
    out["ep_vs_exact"] = float(np.max(np.abs(Ep - Ee)[sel] / Ee[sel]))
    w = 0.0
    for t, _, xs, _ in mm.hist:
        if xs > 1e-3:
            om = 0.5 / 6
            w = max(w, abs(yagi_kunii_eq6(np.clip(xs, 0, 1), om, 1e12)
                           * (1 + om) - t) / (1 + om))
    out["surface_b31"] = w
    # denominator guarded: under the nu=2 defect the r=0 face has zero area,
    # the domain seals, and both gas intake and solid consumption are ~0
    sc = mm.solid_consumed()
    out["mass_balance"] = abs(mm.I_gas - sc) / max(sc, abs(mm.I_gas), 1e-12)
    return out

DEFECTS = [
    (None,                  "baseline (short config)"),
    ("nu2",                 "spherical divergence (nu=2) for the slab"),
    ("sink_area_cubed",     "sink area ~ xi^3 instead of xi^2"),
    ("denom_no_shield",     "shell resistance g*xi for g*xi*(1-xi)"),
    ("rate_denom_mismatch", "Eq.21 and Eq.27 denominators transcribed differently"),
    ("surface_conc_0p9",    "surface concentration 10% low"),
    ("bc_far_neumann",      "far boundary zero-flux instead of psi=0"),
]
rows = {name: metrics_under(d) for d, name in DEFECTS}
keys = ["ep_vs_exact", "surface_b31", "mass_balance", "init_profile"]
base = rows["baseline (short config)"]
print(f"{'defect':<50}" + "".join(f"{k:>15}" for k in keys))
for name, r in rows.items():
    line = f"{name:<50}"
    for k in keys:
        flag = " " if name.startswith("baseline") or r[k] > 3 * base[k] else "*"
        line += f"{r[k]:>14.2e}{flag}"
    print(line)
print("(* = metric did NOT move by >3x: that metric is blind to that defect)\n")

print("What the no-move cells establish:")
print(" - far-boundary defect moves nothing: at L = 30 the profile has decayed")
print("   to ~e^-20, so the truncated Dirichlet and Neumann conditions are")
print("   equivalent. That is the domain-truncation statement, demonstrated")
print("   rather than asserted - and a warning that NO check on this page")
print("   could catch a wrong far-field BC on an adequately long domain.")
print(" - init_profile is blind to the area exponent (xi^2 = xi^3 at xi = 1)")
print("   and to every rate/march defect (it is a theta = 0 check), which is")
print("   why it is never cited as evidence for the transient physics.")
print(" - surface_b31 CATCHES the surface-concentration defect (1.1e-01 vs")
print("   1.8e-03 baseline, a 63x move in the table above): eq. 6's predicted")
print("   theta is compared against the ABSOLUTE march time, so a rescaled")
print("   clock is exactly what it sees. An earlier draft claimed the check")
print("   compared shape only - contradicted by its own table, and corrected.")
print(" - mass_balance separates the defects by CONSISTENCY, not correctness:")
print("   it is blind to the shared-denominator defect (sink and rate change")
print("   together, so the wrong model conserves its own mass - note it even")
print("   IMPROVES there, the shared denominator being smoother), and blind")
print("   to the surface-concentration row for the same reason; it fires on")
print("   the two defects that make gas and solid inconsistent with each")
print("   other - the sink-area exponent and the Eq.21-vs-Eq.27 denominator")
print("   mismatch. (Under nu2 the r=0 face has zero area, the domain seals,")
print("   nothing is consumed, and the metric saturates at 1 by its guarded")
print("   definition - a symptom of the broken geometry, not a balance test.)")
print("\nWhat NO check on this page can catch: a mis-statement of the model")
print("shared by both routes AND by B3.1's closed form - i.e. an error in the")
print("common shrinking-core closure itself. That is bounded by B3.1's own")
print("independent validation (6.9e-16 against a resistance integration) and")
print("by the native-render transcription shown in 'The published model'.")'''))

cells.append(code(r'''report_agreement("B3.2", {
    # headline: the pymrm march against the exact closed form
    "ep_vs_exact_g05_relmax": hl[0.5],
    "ep_vs_exact_g5_relmax": hl[5.0],
    # reduction to the published B3.1 page (simulated surface grain, n = 600)
    "surface_b31_relmax_n600": surf_err[600],
    # exact sinh profile at theta = 0, n = 600
    "init_profile_relmax_n600": init_err[600],
    # observed orders, both axes (last refinement pair)
    "grid_order_observed": grid_order,
    "dt_order_observed": dt_order,
    # temporal mass-balance closure at dt = 0.005
    "mass_balance_rel_dt5e3": mass_bal,
    # dimensional closure on the printed Fig. 7 base case
    "dimensional_closure_relmax": dim_err,
    # late-time lag vs Delta = 1/4 + g/30
    "lag_fit_rel_g05": lag_err[0.5],
    "lag_fit_rel_g5": lag_err[5.0],
    # quantitative results this page adds (regression-protected numbers)
    "zone_width_ell_g500": wg[500.0],
    "zone_width_cm_fig7base": w500_cm,
    "theta_early5_g05": dep[0.5][0],
    "theta_late5_g05": dep[0.5][1],
    # structural identity, recorded but BELOW the CI floor (ABS_FLOOR=1e-12):
    # kept for the record, not regression protection
    "b31_identity_structural": worst_id,
})
import time as _time
print("\ntotal wall-clock so far is printed by the executed notebook's cells;")
print("all solves above are deterministic (fixed grids, direct sparse solves,")
print("bracketed root finds, no warm-start continuation anywhere).")'''))

cells.append(md(r"""## What pymrm adds

**To the semi-infinite grain model's penetration curve: nothing, and the page
proves it.** The exact law $E = \sqrt{6X(\theta)}$ — derived here from the
paper's own Eqs. 21 and 27 via the exposure transformation — replaces the
numerical integration Szekely and Evans ran in 1969, for their exact
formulation (quasi-steady gas, first-order kinetics). Their Appendix 1
computer solution was, for this quantity, avoidable. The closed form is this
page's contribution; it is a near relative of results in the later
Sohn–Szekely lineage (the additive-reaction-time law, planned `B3.4`), and no
novelty is claimed for it.

**What the pymrm march adds is everything the transformation cannot survive.**
The exposure trick dies the moment the gas equation gains a holdup term, the
kinetics stop being first order, the structure evolves (the authors' own
assumption 5), or the geometry becomes a finite pellet (Part II) — while the
finite-volume march above extends to each of those by changing a coefficient
or `nu`. The page's value is having the two side by side: an exact law to
stand on, and the general machine validated against it at second order in both
axes before any generalisation is attempted.

**What the analysis adds beyond the paper**: the diffuse zone the paper could
only plot is quantified — a travelling reaction zone of constant width
$w(g)\,\ell$ (0.295 cm for the printed base case, from printed constants
alone, no $C_0$ required), a closed-form lag $\Delta = 1/4 + g/30$ separating
the grain model from the specimen-scale shrinking core forever, and the 5 %
departure map showing the window where neither textbook limit works growing
from a factor of ~3 in time at $g = 0.125$ to more than four orders of
magnitude at the paper's own base case $g = 500$.

**Honesty**: nothing here is validated against measurement — the source
contains none. The paper's own claim for the model is reproducing "general
trends" with "reasonable values", and this page inherits exactly that
epistemic status, plus internal mathematical consistency."""))

cells.append(md(r"""## Reuse

**`E_exact(theta, g)` is the whole semi-infinite grain model** — one parameter
in, penetration history out — with `G_of_xi`, `X_of_theta`, `xi_of_zeta` and
the exact profiles as by-products. `B3.4` (Sohn's law of additive reaction
times, whose paper cites this one) should import these as its exact reference.
Anyone fitting shrinking-core rate constants to porous-solid data can use
$\Delta = 1/4 + g/30$ to see how much structure their "rate constant" has
absorbed.

**`GrainModel` is the reusable pymrm piece**: a 1-D quasi-steady field coupled
to a pointwise internal-state ODE with a complementarity constraint — the
shape of every deactivation, sorption-front and product-layer problem in
section B. To make it a finite spherical pellet (Part II, `B3.4`'s geometry):
`nu=2`, a symmetry condition `{"a": 1, "b": 0, "d": 0}` at the centre, and the
film condition of `B3.1` at the surface.

**Three traps, all demonstrated above rather than asserted.** (i) `nu` is the
*sample* geometry, not the grain geometry — the grains' sphericity is in the
closure, and `nu=2` here moves the headline metric by four orders of
magnitude. (ii) `NumJac((n, 2))`, never a bare `(n,)` — the 1-D-shape dense
Jacobian trap measured on `B1.1`. (iii) The far boundary of a truncated
semi-infinite domain is *undetectable* by any interior check once $L$ is
adequate — choose $L$ from the physics ($E(\theta_\mathrm{end})$ plus a tail
allowance), because no residual will warn you.

**Reading this paper.** Native resolution is 300 ppi CCITT-G4; render at 300
and **crop before reading numerics** — the mid-dot decimals ($2\!\cdot\!0
\times 10^{-4}$) vanish at page scale and the text layer is unusable for
mathematics. Part II (*CES* **26** (1971) 1901) has no readable scan in this
collection; its Elsevier api-text drops decimal points and must not be used
for constants.

**Related pages.** [`B3.1`](../B3.1-shrinking-core/) (the closure, and the
model this one reduces to at the exposed face), `B3.4` (planned: the
closed-form bridge between the two pictures), `B1.1` (the porous-catalyst
counterpart where the solid is not consumed), `B2.2` (deactivation fronts,
same mathematical shape).

**Cite the source, not this page:** Szekely, J. and Evans, J. W., *A
structural model for gas–solid reactions with a moving boundary*, Chemical
Engineering Science **25**(6), 1091–1107 (1970),
doi:10.1016/0009-2509(70)85053-9."""))

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
