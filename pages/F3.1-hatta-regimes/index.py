# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ---
# title: "Hatta number and the gas–liquid reaction regimes"
# description: "Where reaction and diffusion compete in a liquid film: the exact first-order solution, and how far the Van Krevelen–Hoftijzer and DeCoursey enhancement approximations can be trusted."
# categories: [sec:F, struct:S3, tier:T0, data:tier6, phase:gas-liquid]
# date: 2026-07-26
# ---
#
# # Hatta number and the gas–liquid reaction regimes
#
# **Catalog ID:** `F3.1` · **Structures:** `S3` (1D steady BVP) · **Tier:** T0
#
# A gas dissolving into a reacting liquid can behave in three completely different
# ways depending on one number. Get the regime wrong and an absorber is sized
# against the wrong resistance entirely. This page solves the stagnant-film
# equations with pymrm, checks the solve against the exact pseudo-first-order
# result, and measures how far the Van Krevelen-Hoftijzer and DeCoursey
# approximations depart from the film solution, and where.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.27.10"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "F3.1-hatta-regimes"  # equals the page directory

import os, sys, urllib.request
from pathlib import Path
_shared = next((p / "shared" for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "shared" / "gallery_bootstrap.py").is_file()), None)
if _shared is None:
    _shared = (Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
               / "pymrm-gallery" / GALLERY_REF)
    _shared.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        f"pymrm-gallery/{GALLERY_REF}/shared/gallery_bootstrap.py",
        _shared / "gallery_bootstrap.py")
sys.path.insert(0, str(_shared))
sys.modules.pop("gallery_bootstrap", None)
from gallery_bootstrap import setup
setup(PAGE, PYMRM, GALLERY_REF);  # the semicolon stops the cell echoing the path

# %%
# Imports. Nothing beyond pymrm and the Colab defaults.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import eye_array
from scipy.sparse.linalg import splu
from IPython.display import Markdown, display
from pymrm.checks import observed_orders
from gallery_utils import Claims, load_claims
from model import Film, argext_log, decoursey, err_pct, vkh

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# When a gas dissolves into a liquid that consumes it, reaction and diffusion
# compete inside a thin film at the interface. The **Hatta number** compares the
# two:
#
# $$
# \mathrm{Ha} = \frac{\text{reaction rate in the film}}{\text{diffusion rate through the film}} .
# $$
#
# Three regimes follow, and they call for different equipment:
#
# - $\mathrm{Ha} \ll 1$: **slow reaction.** Nothing reacts in the film; the
#   dissolved gas reaches the bulk and reacts there. Performance scales with
#   *liquid holdup*, so a bubble column beats a packed tower.
# - $\mathrm{Ha} \gg 1$: **fast reaction.** Everything is consumed within the
#   film and the bulk never sees it. Performance scales with *interfacial area*,
#   so the ranking reverses.
# - $\mathrm{Ha}$ large **and** the co-reactant depleted: **instantaneous**
#   regime, where the reaction collapses to a plane inside the film and the rate
#   is set by how fast the co-reactant can be supplied.
#
# The **enhancement factor** $E$ measures how much faster absorption is than
# physical dissolution alone. For a pseudo-first-order reaction it has an exact
# closed form. For the second-order case that actually occurs in practice, CO₂
# into an amine, for instance, it does not, and the field runs on two
# approximations: Van Krevelen and Hoftijzer's implicit expression (1948) and
# DeCoursey's explicit one (1974).
#
# This page reproduces the exact case, then asks how much the approximations cost.

# %% [markdown]
# ## The published model
#
# **Film model.** A stagnant liquid film of thickness $\delta$ at the interface;
# everything beyond it is well mixed. With $\xi = x/\delta$, dissolved gas $a =
# A/A_i$ and co-reactant $b = B/B_b$, an irreversible reaction $A + \nu B \to$
# products gives
#
# $$
# \frac{\mathrm{d}^2 a}{\mathrm{d}\xi^2} = \mathrm{Ha}^2\, a\, b , \qquad
# \frac{\mathrm{d}^2 b}{\mathrm{d}\xi^2} = \frac{\mathrm{Ha}^2}{q}\, a\, b ,
# $$
#
# $$
# \mathrm{Ha} = \frac{\delta\sqrt{k_2 B_b/D_A}}{1} = \frac{\sqrt{k_2 B_b D_A}}{k_L},
# \qquad
# q = \frac{D_B B_b}{\nu D_A A_i}, \qquad E_i = 1 + q ,
# $$
#
# with $a(0)=1$ (saturated interface), $a(1)=0$ (depleted bulk), $b(1)=1$, and
# $\mathrm{d}b/\mathrm{d}\xi|_0 = 0$ because the co-reactant is non-volatile.
# $E_i$ is the **instantaneous** enhancement, the ceiling set by co-reactant
# supply.
#
# The enhancement factor is the interfacial flux normalised by the purely
# diffusive one,
# $E = -\left.\mathrm{d}a/\mathrm{d}\xi\right|_{\xi=0}$.
#
# **Exact, pseudo-first-order** ($b \equiv 1$, co-reactant in large excess):
#
# $$
# E = \frac{\mathrm{Ha}}{\tanh \mathrm{Ha}} \;\longrightarrow\;
# \begin{cases} 1 & \mathrm{Ha}\ll 1 \\ \mathrm{Ha} & \mathrm{Ha}\gg 1 \end{cases}
# $$
#
# **Van Krevelen–Hoftijzer (1948)**, implicit in $E$:
#
# $$
# E = \frac{\mathrm{Ha}\sqrt{(E_i-E)/(E_i-1)}}{\tanh\!\left(\mathrm{Ha}\sqrt{(E_i-E)/(E_i-1)}\right)}
# $$
#
# **DeCoursey (1974)**, explicit:
#
# $$
# E = -\frac{\mathrm{Ha}^2}{2(E_i-1)} +
# \sqrt{\frac{\mathrm{Ha}^4}{4(E_i-1)^2} + \frac{E_i\,\mathrm{Ha}^2}{E_i-1} + 1}
# $$
#
# DeCoursey derived his relation for Danckwerts' surface-renewal model, not for
# the film (the title of the 1974 paper says so). With the co-reactant in excess,
# $E_i \to \infty$, it becomes $\sqrt{1+\mathrm{Ha}^2}$, Danckwerts' own
# pseudo-first-order result, while the film gives $\mathrm{Ha}/\tanh\mathrm{Ha}$.
# Van Krevelen-Hoftijzer becomes $\mathrm{Ha}/\tanh\mathrm{Ha}$ exactly.
#
# | Symbol | Code | Meaning |
# |---|---|---|
# | $a, b$ | `u[:, 0]`, `u[:, 1]` | dimensionless concentrations in the film |
# | $\mathrm{Ha}$ | `ha` | Hatta number |
# | $E_i$ | `e_i` | instantaneous enhancement, $1+q$ |
# | $E$ | `Film.enhancement(u)` | enhancement factor from the interfacial flux |
# | VKH, DeCoursey | `vkh(ha, e_i)`, `decoursey(ha, e_i)` | the two approximations, in `model.py` |

# %% [markdown]
# ## Parameters and assumptions
#
# **Assumptions:** steady state; stagnant film (no surface renewal); irreversible
# second-order kinetics; constant diffusivities; the bulk fully depleted in
# dissolved gas, $a(1)=0$; co-reactant non-volatile.
#
# The problem is dimensionless: only $\mathrm{Ha}$ and $E_i$ matter.
#
# **What serves as "data" here.** As with `B1.1`, the sources are analytical:
# Hatta's regimes and the first-order enhancement are exact, and the Van
# Krevelen–Hoftijzer and DeCoursey expressions are published *approximations*, not
# measurements. This page is provenance tier 6: validated against exact solutions
# and limiting cases, not against experiment. The interesting comparison is
# therefore not "does pymrm match the data" but "how far do the standard
# approximations depart from the equation they approximate", which is the question
# a practitioner actually faces.

# %%
# The grids and parameter sets of the page (assumed, dimensionless).
N_X = 400                                   # film cells, uniform
HA1 = np.geomspace(0.02, 200.0, 70)         # pseudo-first-order sweep
HA2 = np.geomspace(0.05, 300.0, 55)         # second-order sweep, and the extremum scan
E_I = (5.0, 20.0, 100.0)                    # instantaneous enhancement factors swept
HA_ORDER = 5.0                              # Hatta number of the grid study
LADDER = (100, 200, 400, 800)               # nested grids of the grid study
HA_ORDER2 = 2.0                             # Hatta number of the second grid study
LADDER2 = (400, 800, 1600, 3200)            # nested grids of the second grid study
N_FINE = 1600                               # film that checks the worst points
E_I_EXCESS = 1.0e6                          # co-reactant in large excess
HA_EXCESS = (0.5, 2.0, 10.0)                # Hatta numbers of the excess check

# %% [markdown]
# ## The data
#
# None. The comparison is between the film equations, their exact
# pseudo-first-order solution and two published approximations, all computed
# here; no dataset is loaded. `F3.2` prices Van Krevelen-Hoftijzer on its own
# reference and reconciles with this page; its numbers are loaded with
# `load_claims` in the Validation section, never retyped.

# %% [markdown]
# ## PyMRM implementation
#
# The film is a slab (`nu=0`) on $\xi \in [0,1]$ with a uniform grid, fields on
# the last axis, shape `(n_x, n_c)`. The two-component case carries per-species
# boundary conditions in the same `{"a","b","d"}` dictionaries; the coefficients
# are arrays over the component axis, which is how pymrm expresses "Dirichlet
# for $a$, Neumann for $b$" at the same face without any special casing. The
# diffusion operator is assembled once per film. With $b \equiv 1$ the problem
# is linear and takes one sparse solve; the second-order problem is solved by
# `newton` with a `NumJac((n_x, 2))` Jacobian of the pointwise rate, raising
# unless it converges. $E$ is read on the interface face from the face-flux
# operator, the flux the solve transports. The `Film` class is in `model.py`:
#
# ```python
# # a dc/dn + b c = d, n the OUTWARD normal (n = -xi at xi = 0, +xi at 1)
# # xi = 0, interface: a = 1 (saturated) and -db/dxi = 0 (B non-volatile)
# left = {"a": [[0.0, 1.0]], "b": [[1.0, 0.0]], "d": [[1.0, 0.0]]}
# # xi = 1, bulk edge: a = 0 (depleted bulk) and b = 1 (bulk value)
# right = {"a": [[0.0, 0.0]], "b": [[1.0, 1.0]], "d": [[0.0, 1.0]]}
# ```
#
# The extremum helper `argext_log` scans the whole $\mathrm{Ha}$ window and
# refines the largest $|f|$ as the single root of its central difference
# (`pymrm.checks.find_roots`); a maximum on a window end is reported as "at or
# beyond the end". Van Krevelen-Hoftijzer is solved for $E$ as the single root
# on $[1, E_i]$ the same way.

# %%
# The studies. Every number the page states comes out of one of these
# functions, and every break-table row calls them with one argument changed.
EPS = np.finfo(float).eps


def first_order_study(n_x=N_X, flux_read="face"):
    """Pseudo-first-order film over HA1 against the exact Ha/tanh(Ha)."""
    film = Film(n_x, 1, flux_read=flux_read)
    E = np.array([film.solve_first_order(h)[1] for h in HA1])
    exact = HA1 / np.tanh(HA1)
    rel = np.abs(E - exact) / exact
    return {"first_order_max_rel_error": float(rel.max()),
            # where the error sits: the largest error below Ha = 20 over the largest
            "first_order_worst_below_ha20_share": float(rel[HA1 < 20].max() / rel.max()),
            # reaction-layer thickness 1/Ha at the largest Ha, in cells of 1/n_x
            "reaction_layer_cells_at_ha_max": float(n_x / HA1[-1]),
            "_E": E, "_rel": rel, "_film": film}


def cond_inf(film, ha):
    """cond_inf of the linear film matrix, exact: the matrix is an M-matrix, so
    its inverse is non-negative and ||A^-1||_inf is the largest entry of A^-1 1."""
    A = film.jac_diff + ha**2 * eye_array(film.n_x, format="csc")
    inv_one = splu(A).solve(np.ones(film.n_x))
    if inv_one.min() <= 0.0:
        raise RuntimeError("film matrix is not an M-matrix; the cond_inf identity fails")
    return float(abs(A).sum(axis=1).max() * inv_one.max())


def order_study(ladder=LADDER, flux_read="face", ha=HA_ORDER, expected=3.0, suffix=""):
    """Observed order of E at Hatta number `ha` from successive differences.

    The order farthest from `expected` is returned as "convergence_order" + suffix,
    the round-off companion as "order_roundoff_margin" + suffix.
    """
    films = {n: Film(n, 1, flux_read=flux_read) for n in ladder}
    st = observed_orders(lambda n: films[n].solve_first_order(ha)[1], ladder)
    o = np.asarray(st["orders"], float)
    v = np.asarray(st["values"], float)
    # round-off companion: the finest successive difference over the round-off
    # floor of the direct solve, eps cond_inf(A) |E|, on the finest grid
    floor = EPS * cond_inf(films[ladder[-1]], ha) * abs(v[-1])
    return {"convergence_order" + suffix: float(o[np.argmax(np.abs(o - expected))]),
            "order_roundoff_margin" + suffix: float(abs(v[-1] - v[-2]) / floor),
            "_orders": o, "_values": v}


def order_study_ha2(ladder=LADDER2, flux_read="face"):
    """The grid study at Ha = HA_ORDER2, held against second order."""
    return order_study(ladder=ladder, flux_read=flux_read, ha=HA_ORDER2, expected=2.0,
                       suffix="_ha2")


def sweep_study(e_is=E_I, n_x=N_X, denom="q", q_power=2):
    """Second-order film over HA2 at each E_i: the largest sampled errors (%)."""
    film = Film(n_x, 2)
    res = {}
    for e_i in e_is:
        res[e_i] = np.array([(film.solve_second_order(h, e_i)[1], vkh(h, e_i, denom),
                              decoursey(h, e_i, q_power)) for h in HA2])
    return {"vkh_max_error_pct": float(max(
                np.max(100 * np.abs(r[:, 1] - r[:, 0]) / r[:, 0]) for r in res.values())),
            "decoursey_max_error_pct": float(max(
                np.max(100 * np.abs(r[:, 2] - r[:, 0]) / r[:, 0]) for r in res.values())),
            "_res": res}


def worst_study(name, e_is=E_I, n_x=N_X, refine=True, denom="q", q_power=2):
    """Largest |error| of one approximation over HA2 and E_I, refined in Ha.

    Raises when the overall maximum sits on a window end (at or beyond it)."""
    film = Film(n_x, 2)
    approx = ((lambda h, e: vkh(h, e, denom)) if name == "vkh"
              else (lambda h, e: decoursey(h, e, q_power)))
    table = {e: argext_log(lambda h, e=e: err_pct(h, e, film, approx),
                           HA2[0], HA2[-1], n_scan=len(HA2), refine=refine)
             for e in e_is}
    e_w = max(table, key=lambda e: abs(table[e][1]))
    x, f, interior = table[e_w]
    if refine and not interior:
        raise RuntimeError(f"{name}: largest error at the window end Ha = {x:g}")
    return {f"{name}_worst_error_pct": f, f"{name}_worst_ha": x,
            f"{name}_worst_ha_over_ei": x / e_w, "_e_i": e_w, "_table": table}


def worst_grid_study(name, n_x=N_X):
    """Relative change of the worst-point error when the film is refined to N_FINE."""
    w = worst_study(name, n_x=n_x)
    approx = vkh if name == "vkh" else decoursey
    fine = err_pct(w[f"{name}_worst_ha"], w["_e_i"], Film(N_FINE, 2), approx)
    return {f"{name}_worst_grid_rel_change": abs(fine / w[f"{name}_worst_error_pct"] - 1.0)}


def pfo_limit_study(reference="film"):
    """DeCoursey with the co-reactant in excess (E_i = 1e12) against the
    pseudo-first-order film result (or, as a row, Danckwerts' sqrt(1 + Ha^2))."""
    ref = ((lambda h: h / np.tanh(h)) if reference == "film"
           else (lambda h: np.sqrt(1.0 + h * h)))
    x, f, interior = argext_log(lambda h: 100 * (decoursey(h, 1e12) / ref(h) - 1.0),
                                HA2[0], HA2[-1], n_scan=len(HA2))
    return {"decoursey_pfo_limit_max_error_pct": abs(f), "_ha": x, "_interior": interior}


def limits_study(a_interface=1.0, e_i_large=20.0):
    """E at the two ends of HA2: E_i = 20 at Ha = 0.05, E_i = e_i_large at Ha = 300."""
    film = Film(N_X, 2, a_interface=a_interface)
    h0, h1 = HA2[0], HA2[-1]
    lo = film.solve_second_order(h0, 20.0)[1]
    hi = film.solve_second_order(h1, e_i_large)[1]
    return {"limit_E_small_Ha": lo,
            # the small-Ha limit as a rate: (E - 1)/(Ha^2/3), 1 for b = 1
            "small_ha_rate_over_closed_form": (lo - 1.0) / (h0**2 / 3.0),
            "limit_E_large_Ha": hi}


def excess_study(e_i=E_I_EXCESS):
    """Second order with B in excess against first order: exact, and same grid."""
    f2, f1 = Film(N_X, 2), Film(N_X, 1)
    dev, same = [], []
    for h in HA_EXCESS:
        E2 = f2.solve_second_order(h, e_i)[1]
        dev.append(abs(E2 - h / np.tanh(h)) / (h / np.tanh(h)))
        same.append(abs(E2 / f1.solve_first_order(h)[1] - 1.0))
    return {"excess_reduces_to_first_order": float(max(dev)),
            "excess_vs_first_order_same_grid": float(max(same))}

# %% [markdown]
# ## Results
#
# ### The exact case
#
# Pseudo-first-order, where $E = \mathrm{Ha}/\tanh\mathrm{Ha}$ exactly. Both
# asymptotes are visible: $E \to 1$ (no enhancement) and $E \to \mathrm{Ha}$.

# %%
FO = first_order_study()
film1 = FO["_film"]
E_exact = HA1 / np.tanh(HA1)

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
axes[0].loglog(HA1, E_exact, "k-", lw=2, label=r"exact $\mathrm{Ha}/\tanh\mathrm{Ha}$")
axes[0].loglog(HA1[::3], FO["_E"][::3], "o", color="tab:blue", ms=6,
               mfc="white", mew=1.6, label="pymrm")
axes[0].loglog(HA1, np.ones_like(HA1), ":", color="0.5", lw=1.2)
axes[0].loglog(HA1, HA1, ":", color="0.5", lw=1.2)
axes[0].text(0.05, 1.25, "slow: $E\\to1$", fontsize=9)
axes[0].text(30, 8, "fast: $E\\to\\mathrm{Ha}$", fontsize=9)
axes[0].set_xlabel(r"Hatta number, $\mathrm{Ha}$ [-]")
axes[0].set_ylabel("enhancement factor, $E$ [-]")
axes[0].set_title("Pseudo-first-order")
axes[0].legend(frameon=False)

for h, c in zip([0.3, 1.0, 3.0, 10.0], plt.cm.viridis(np.linspace(0.1, 0.85, 4))):
    axes[1].plot(film1.x_c, film1.solve_first_order(h)[0], color=c, lw=2,
                 label=rf"$\mathrm{{Ha}}$ = {h}")
axes[1].set_xlabel(r"position in the film, $\xi$ [-]")
axes[1].set_ylabel("$a$ [-]")
axes[1].set_title("Where the gas is consumed")
axes[1].legend(frameon=False)
fig.tight_layout(); plt.show()

display(Markdown(
    f"Largest relative error against the exact solution over the sweep: "
    f"{FO['first_order_max_rel_error']:.3e}, concentrated at high Ha: the "
    f"largest error below Ha = 20 is {FO['first_order_worst_below_ha20_share']:.1e} "
    f"of it. At Ha = {HA1[-1]:g} the reaction layer, of thickness 1/Ha, spans "
    f"{FO['reaction_layer_cells_at_ha_max']:.0f} cells of the {N_X}-cell film."))

# %% [markdown]
# ### Second order, and what the approximations cost
#
# Now the co-reactant is depleted inside the film. $E$ can no longer exceed the
# instantaneous ceiling $E_i$, and neither published formula is exact. Each film
# solve starts cold from the pseudo-first-order profile at its own Ha.

# %%
SW = sweep_study()
res = SW["_res"]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
cols = plt.cm.plasma(np.linspace(0.15, 0.75, len(E_I)))
for e_i, c in zip(E_I, cols):
    r = res[e_i]
    axes[0].loglog(HA2, r[:, 0], color=c, lw=2.2, label=rf"$E_i$ = {e_i:.0f}")
    axes[0].loglog(HA2, r[:, 1], "--", color=c, lw=1.4)
    axes[0].loglog(HA2, r[:, 2], ":", color=c, lw=1.6)
    axes[0].axhline(e_i, color=c, lw=0.7, ls="-.", alpha=0.5)
    axes[1].semilogx(HA2, 100 * (r[:, 1] - r[:, 0]) / r[:, 0], "--", color=c, lw=1.6,
                     label=rf"VKH, $E_i$ = {e_i:.0f}")
    axes[1].semilogx(HA2, 100 * (r[:, 2] - r[:, 0]) / r[:, 0], ":", color=c, lw=2,
                     label=rf"DeCoursey, $E_i$ = {e_i:.0f}")
axes[1].semilogx(HA2, 100 * (np.sqrt(1 + HA2**2) * np.tanh(HA2) / HA2 - 1), "-",
                 color="0.6", lw=1.2, label=r"$\sqrt{1+\mathrm{Ha}^2}$ vs film, $b=1$")
axes[0].loglog(HA2, HA2 / np.tanh(HA2), color="0.6", lw=1.2,
               label="pseudo-first-order")
axes[0].set_xlabel(r"$\mathrm{Ha}$ [-]"); axes[0].set_ylabel("$E$ [-]")
axes[0].set_title("solid pymrm · dashed VKH · dotted DeCoursey")
axes[0].legend(frameon=False, fontsize=8)
axes[1].axhline(0, color="k", lw=0.8)
axes[1].set_xlabel(r"$\mathrm{Ha}$ [-]")
axes[1].set_ylabel("approximation error [%]")
axes[1].set_title("Departure from the full solution")
axes[1].legend(frameon=False, fontsize=7, ncol=2)
fig.tight_layout(); plt.show()

# the largest |error| per E_i, refined in Ha from a scan of the whole window
WV, WD = worst_study("vkh"), worst_study("decoursey")
rows = []
for e_i in E_I:
    for lab, w in (("VKH", WV), ("DeCoursey", WD)):
        x, f, interior = w["_table"][e_i]
        fine = err_pct(x, e_i, Film(N_FINE, 2), vkh if lab == "VKH" else decoursey)
        rows.append((lab, f"{e_i:g}", f"{f:+.3f}",
                     f"{x:.4g}" if interior else f"at or beyond {x:g}",
                     f"{abs(fine / f - 1):.0e}"))
print(pd.DataFrame(rows, columns=["approximation", "E_i", "worst error [%]", "at Ha",
                                  f"rel. change on {N_FINE} cells"]).to_string(index=False))

# %%
# The concentration profiles behind the three regimes.
film2 = Film(N_X, 2)
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.7), sharey=True)
for ax, h, tag in zip(axes, [0.3, 3.0, 50.0],
                      ["slow: reaction happens in the bulk",
                       "fast: consumed inside the film",
                       "instantaneous: reaction plane"]):
    u, E = film2.solve_second_order(h, 20.0)
    ax.plot(film2.x_c, u[:, 0], color="tab:blue", lw=2.2, label="$a$ (gas)")
    ax.plot(film2.x_c, u[:, 1], color="tab:red", lw=2.2, label="$b$ (co-reactant)")
    ax.set_title(f"Ha = {h:g},  E = {E:.2f}\n{tag}", fontsize=9)
    ax.set_xlabel(r"$\xi$ [-]")
axes[0].set_ylabel("dimensionless concentration [-]")
axes[0].legend(frameon=False)
fig.tight_layout(); plt.show()

# %% [markdown]
# ### Where each approximation is worst, and why

# %%
PL = pfo_limit_study()
share = abs(WD["decoursey_worst_error_pct"]) / PL["decoursey_pfo_limit_max_error_pct"]
display(Markdown(
    f"Van Krevelen-Hoftijzer is worst at E_i = {WV['_e_i']:g}: "
    f"{WV['vkh_worst_error_pct']:+.4f} % at Ha = {WV['vkh_worst_ha']:.3f}, that is "
    f"Ha = {WV['vkh_worst_ha_over_ei']:.2f} E_i, in the transition from the fast to "
    f"the instantaneous regime. DeCoursey is worst at E_i = {WD['_e_i']:g}: "
    f"{WD['decoursey_worst_error_pct']:+.4f} % at Ha = {WD['decoursey_worst_ha']:.3f}, "
    f"Ha = {WD['decoursey_worst_ha_over_ei']:.4f} E_i, which is the transition from "
    f"the slow to the fast regime and far from the instantaneous one. With the "
    f"co-reactant in excess DeCoursey becomes Danckwerts' sqrt(1 + Ha^2), and that "
    f"differs from the film's Ha/tanh(Ha) by up to "
    f"{PL['decoursey_pfo_limit_max_error_pct']:.3f} % (at Ha = {PL['_ha']:.3f}); "
    f"DeCoursey's worst second-order error is {share:.4f} of that. Its largest "
    f"error on this page is the difference between the film model and the "
    f"surface-renewal model it was derived for, not a second-order effect."))

# %% [markdown]
# ## Validation
#
# Checks against exact results and limits, grid studies at two Hatta numbers, a split of each
# approximation's worst error from the film's discretisation, and the
# reconciliation with `F3.2`. None uses fitted data.

# %%
# 1. First order against the exact closed form (above).
print(f"1. max rel. error vs Ha/tanh(Ha), {N_X} cells : {FO['first_order_max_rel_error']:.3e}")

# 2. Grid study at Ha = HA_ORDER: orders from successive differences.
OS = order_study()
print(f"2. E at Ha = {HA_ORDER:g} on {LADDER} cells: orders "
      + ", ".join(f"{o:.4f}" for o in OS["_orders"])
      + f"; round-off margin of the finest difference {OS['order_roundoff_margin']:.1e}")
OS2 = order_study_ha2()
print(f"   E at Ha = {HA_ORDER2:g} on {LADDER2} cells: orders "
      + ", ".join(f"{o:.4f}" for o in OS2["_orders"])
      + f"; round-off margin of the finest difference {OS2['order_roundoff_margin_ha2']:.1e}")

# 3. Limits of the second-order solution at the ends of the sweep.
LI = limits_study()
print(f"3. E at Ha = {HA2[0]:g}, E_i = 20      : {LI['limit_E_small_Ha']:.6f}"
      f"   (E - 1)/(Ha^2/3) = {LI['small_ha_rate_over_closed_form']:.5f}")
print(f"   E at Ha = {HA2[-1]:g}, E_i = 20      : {LI['limit_E_large_Ha']:.4f}"
      f"   ({100 * (1 - LI['limit_E_large_Ha'] / 20.0):.2f} % under E_i)")

# 4. Second order reduces to pseudo-first order when B is in large excess.
EX = excess_study()
print(f"4. E_i = {E_I_EXCESS:.0e}: max rel. dev. from Ha/tanh(Ha) "
      f"{EX['excess_reduces_to_first_order']:.3e}; from the first-order solve on "
      f"the same grid {EX['excess_vs_first_order_same_grid']:.3e}")

# 5. The worst points refined against a finer film: the approximation errors
#    are not the 400-cell film's discretisation.
GV, GD = worst_grid_study("vkh"), worst_grid_study("decoursey")
print(f"5. worst-point error, {N_X} against {N_FINE} cells: VKH "
      f"{GV['vkh_worst_grid_rel_change']:.1e}, DeCoursey "
      f"{GD['decoursey_worst_grid_rel_change']:.1e} relative")


def share_study(q_power=2, reference="film"):
    """DeCoursey's worst error over its pseudo-first-order (film vs renewal) limit."""
    w = WD if q_power == 2 else worst_study("decoursey", q_power=q_power)
    return {"decoursey_worst_share_of_pfo_limit": abs(w["decoursey_worst_error_pct"])
            / pfo_limit_study(reference)["decoursey_pfo_limit_max_error_pct"]}


# 6. Reconciliation with F3.2, which root-finds the same three E_i on its own
#    graded, extrapolated reference. Loaded, never retyped.
F32 = load_claims("F3.2-van-krevelen-hoftijzer")


def f32_study(w=None):
    w = WV if w is None else w
    return {"vkh_worst_vs_f32_rootfound_rel":
            abs(abs(w["vkh_worst_error_pct"]) / F32["f31_definition_rootfound_pct"] - 1.0)}


F32S = f32_study()
print(f"6. F3.2, root-found on its reference, same three E_i: "
      f"{F32['f31_definition_rootfound_pct']:.4f} %; here {abs(WV['vkh_worst_error_pct']):.4f} %;"
      f" relative gap {F32S['vkh_worst_vs_f32_rootfound_rel']:.0e}")
print(f"   F3.2 also finds this page's sampled VKH maximum within "
      f"{F32['f31_stored_vs_recomputed_rel']:.0e} relative on its reference, and states it")
print("   is a sampled maximum on 55 Ha points at three E_i: the refined values above")
print("   are the extrema; the sampled ones are kept as the sweep's own numbers.")

# %%
# Register every number the page states.
claims.value("first_order_max_rel_error", FO["first_order_max_rel_error"], rtol=1e-6)
claims.bound("first_order_worst_below_ha20_share", FO["first_order_worst_below_ha20_share"],
             max=0.01, basis="chosen", note="the error grows with Ha: below Ha = 20 "
             "it stays under 1 % of the largest one")
claims.value("reaction_layer_cells_at_ha_max", FO["reaction_layer_cells_at_ha_max"], rtol=1e-12)
claims.bound("convergence_order", OS["convergence_order"], min=2.9, max=3.1, basis="chosen",
             note="observed order of the interface flux at Ha = 5 on the nested ladder "
             "100 to 800 cells, a pre-asymptotic range; the order farthest from 3 is held")
claims.bound("order_roundoff_margin", OS["order_roundoff_margin"], min=100.0, basis="chosen",
             note="finest successive difference over eps cond_inf(A) |E| of the direct solve")
claims.bound("convergence_order_ha2", OS2["convergence_order_ha2"], min=1.9, max=2.1,
             basis="chosen", note="observed order of the interface flux at Ha = 2 on the "
             "nested ladder 400 to 3200 cells, the second order of the scheme; the order "
             "farthest from 2 is held")
claims.bound("order_roundoff_margin_ha2", OS2["order_roundoff_margin_ha2"], min=3.0,
             basis="chosen", note="finest successive difference over eps cond_inf(A) |E| "
             "of the direct solve, at Ha = 2 on 3200 cells; cond_inf is the worst-case "
             "amplification, so this floor is an upper bound on the round-off")
claims.value("limit_E_small_Ha", LI["limit_E_small_Ha"], rtol=1e-6)
claims.value("small_ha_rate_over_closed_form", LI["small_ha_rate_over_closed_form"], rtol=1e-6)
claims.value("limit_E_large_Ha", LI["limit_E_large_Ha"], rtol=1e-6)
claims.value("excess_reduces_to_first_order", EX["excess_reduces_to_first_order"], rtol=1e-3)
claims.value("excess_vs_first_order_same_grid", EX["excess_vs_first_order_same_grid"], rtol=1e-3)
claims.value("vkh_max_error_pct", SW["vkh_max_error_pct"], rtol=1e-6)
claims.value("decoursey_max_error_pct", SW["decoursey_max_error_pct"], rtol=1e-6)
for _k in ("vkh_worst_error_pct", "decoursey_worst_error_pct"):
    claims.value(_k, (WV if _k.startswith("vkh") else WD)[_k], rtol=1e-6)
for _k in ("vkh_worst_ha", "vkh_worst_ha_over_ei", "decoursey_worst_ha",
           "decoursey_worst_ha_over_ei"):
    claims.value(_k, (WV if _k.startswith("vkh") else WD)[_k], rtol=1e-4)
claims.value("decoursey_pfo_limit_max_error_pct", PL["decoursey_pfo_limit_max_error_pct"],
             rtol=1e-6)
claims.value("decoursey_worst_share_of_pfo_limit",
             share_study()["decoursey_worst_share_of_pfo_limit"], rtol=1e-6)
claims.bound("vkh_worst_grid_rel_change", GV["vkh_worst_grid_rel_change"], max=1e-5, basis="chosen",
             note=f"discretisation of the {N_X}-cell film at the worst point, "
             f"against {N_FINE} cells")
claims.bound("decoursey_worst_grid_rel_change", GD["decoursey_worst_grid_rel_change"],
             max=1e-5, basis="chosen", note=f"discretisation of the {N_X}-cell film at the worst point, "
             f"against {N_FINE} cells")
claims.bound("vkh_worst_vs_f32_rootfound_rel", F32S["vkh_worst_vs_f32_rootfound_rel"],
             max=1e-4, basis="chosen", note="F3.2's graded, extrapolated reference against this page's "
             "uniform 400-cell film, at the same three E_i")
print(f"{len(claims.metrics)} claims registered")

# %% [markdown]
# ### Break table
#
# Each row injects one defect through an argument of the study that computes the
# claims it targets, and re-runs that study. Two rows make an approximation exact
# by construction: Van Krevelen-Hoftijzer at $E_i = 10^6$, and DeCoursey's
# excess limit compared with Danckwerts' $\sqrt{1+\mathrm{Ha}^2}$ instead of the
# film. `reaction_layer_cells_at_ha_max` is arithmetic on the grid and the sweep,
# moved only by the grid.

# %%
def take(study, *keys):
    """A row's claims from one study call; a study that raises gives nan."""
    try:
        out = study()
        return {k: out[k] for k in keys}
    except Exception as exc:                 # noqa: BLE001 - a failed study is a result
        print(f"  study raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
        return {k: np.nan for k in keys}


W_KEYS = ("_worst_error_pct", "_worst_ha", "_worst_ha_over_ei")
ROWS = [
    ("pseudo-first-order sweep on n_x = 200", lambda: take(
        lambda: first_order_study(n_x=200), "first_order_max_rel_error",
        "reaction_layer_cells_at_ha_max")),
    ("interface flux read from the first cell centre", lambda: {
        **take(lambda: first_order_study(flux_read="cell"),
               "first_order_worst_below_ha20_share"),
        **take(lambda: order_study(flux_read="cell"), "convergence_order"),
        **take(lambda: order_study_ha2(flux_read="cell"), "convergence_order_ha2")}),
    ("order ladder moved to 800 to 6400 cells (past the third-order range)", lambda: take(
        lambda: order_study(ladder=(800, 1600, 3200, 6400)), "convergence_order",
        "order_roundoff_margin")),
    ("Ha = 2 order ladder moved to 800 to 6400 cells (into round-off)", lambda: take(
        lambda: order_study_ha2(ladder=(800, 1600, 3200, 6400)),
        "order_roundoff_margin_ha2")),
    ("interface held at a(0) = 1.01", lambda: take(
        lambda: limits_study(a_interface=1.01), "limit_E_small_Ha",
        "small_ha_rate_over_closed_form")),
    ("E_i = 19 at the large-Ha end (q 5 % low)", lambda: take(
        lambda: limits_study(e_i_large=19.0), "limit_E_large_Ha")),
    ("E_i = 1e3 in the excess check", lambda: take(
        lambda: excess_study(e_i=1e3), "excess_reduces_to_first_order",
        "excess_vs_first_order_same_grid")),
    ("VKH divided by E_i instead of E_i - 1", lambda: {
        **take(lambda: sweep_study(denom="ei"), "vkh_max_error_pct"),
        **take(lambda: (lambda w: {**w, **f32_study(w)})(worst_study("vkh", denom="ei")),
               *("vkh" + k for k in W_KEYS), "vkh_worst_vs_f32_rootfound_rel")}),
    ("VKH at E_i = 1e6 only (exact by construction)", lambda: take(
        lambda: sweep_study(e_is=(1e6,)), "vkh_max_error_pct")),
    ("DeCoursey with q under Ha^4 instead of q^2", lambda: {
        **take(lambda: sweep_study(q_power=1), "decoursey_max_error_pct"),
        **take(lambda: share_study(q_power=1), "decoursey_worst_share_of_pfo_limit"),
        **take(lambda: worst_study("decoursey", q_power=1),
               "decoursey_worst_error_pct")}),
    ("worst points read off the 55-point scan, not refined", lambda: {
        **take(lambda: worst_study("vkh", refine=False), *("vkh" + k for k in W_KEYS)),
        **take(lambda: worst_study("decoursey", refine=False),
               *("decoursey" + k for k in W_KEYS))}),
    ("DeCoursey limit against sqrt(1 + Ha^2) (exact by construction)", lambda: {
        **take(lambda: pfo_limit_study(reference="renewal"),
               "decoursey_pfo_limit_max_error_pct"),
        **take(lambda: share_study(reference="renewal"),
               "decoursey_worst_share_of_pfo_limit")}),
    ("worst-point film on n_x = 50", lambda: {
        **take(lambda: worst_grid_study("vkh", n_x=50), "vkh_worst_grid_rel_change"),
        **take(lambda: worst_grid_study("decoursey", n_x=50),
               "decoursey_worst_grid_rel_change")}),
]

claims.define_rows([r[0] for r in ROWS])
BREAK = []
for row, run in ROWS:
    for name, out in run().items():
        caught_ = claims.breaks(row, name, out)
        base_ = claims.metrics[name]["value"]
        if claims.metrics[name]["kind"] == "bound":   # held to its bound, not its digits
            base_, out = f"(here {base_:.1e})", f"(here {float(out):.1e})"
        BREAK.append((row, name, base_, out, caught_))

bt = pd.DataFrame(BREAK, columns=["injected defect", "metric", "baseline", "with defect",
                                  "caught"])
pd.set_option("display.max_colwidth", 60)
print(bt.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
missing = [n for n in claims.metrics if not claims.metrics[n].get("structural")
           and n not in set(bt.metric)]
print(f"rows that fail no claim: {idle or 'none'}; claims with no row: {missing or 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# **For the pseudo-first-order case, nothing**, and as on `B1.1`, that is the
# point. The closed form is exact; the film solve reproduces it to the error
# printed under the first figure, with the observed orders of the grid studies, and
# that is what licenses the second-order result, where no closed form exists.
#
# **For the second-order case, an error bar on the two standard formulas, at
# three values of $E_i$.** The film equations are solved with pymrm operators
# and `newton` on a `NumJac((n_x, 2))` Jacobian, every solve checked; the worst
# point of each approximation is refined from a scan of the whole $\mathrm{Ha}$
# window by `pymrm.checks.find_roots`, and the film's own discretisation is
# shown to be a small part of it. The measurement says *why* each formula is
# worst where it is: Van Krevelen-Hoftijzer in the fast-to-instantaneous
# transition, DeCoursey at the slow-to-fast transition, where nearly all of its
# error is the difference between the film model and the surface-renewal model
# it was derived for. `F3.2` maps Van Krevelen-Hoftijzer over the whole chart
# with a continuously varying $q$.
#
# The full solve is not expensive: a few hundred cells and a Newton iteration
# per point.

# %% [markdown]
# ## Reuse
#
# | Change | Where |
# |---|---|
# | Reversible reaction / equilibrium | extend `rate` and add the product as a third component |
# | Surface renewal instead of a film | replace the steady BVP with a transient penetration solve (`S4`) |
# | Non-depleted bulk, $a(1) = a_b > 0$ | change `d` in the right-hand boundary dictionary |
# | Gas-side resistance | make the interface a Robin condition instead of Dirichlet |
# | Couple to a column | this is the local closure inside `F2.1`/`F3.5` |
#
# DeCoursey's relation belongs with a surface-renewal model; compared with a
# film solution, most of its error at the slow-to-fast transition is that model
# difference, so do not read it as a second-order correction to the film.
#
# **Related pages:** `F3.2` Van Krevelen–Hoftijzer · `F3.3` DeCoursey ·
# `F3.5` CO₂–amine absorption · `B1.1` Thiele (the same `S3` equation with the
# boundary conditions swapped) · `A4.9` Duncan–Toor.
#
# ## References
#
# Hatta, S. (1932). *On the absorption velocity of gases by liquids.* Technological
# Reports of Tôhoku University **10**, 119–128.
#
# Van Krevelen, D. W. and Hoftijzer, P. J. (1948). *Kinetics of gas-liquid
# reactions. Part I. General theory.* Recueil des Travaux Chimiques des Pays-Bas
# **67**(7), 563–586. [doi:10.1002/recl.19480670708](https://doi.org/10.1002/recl.19480670708)
#
# DeCoursey, W. J. (1974). *Absorption with chemical reaction: development of a new
# relation for the Danckwerts model.* Chemical Engineering Science **29**(9),
# 1867–1872. [doi:10.1016/0009-2509(74)85003-7](https://doi.org/10.1016/0009-2509(74)85003-7)
#
# Danckwerts, P. V. (1970). *Gas-Liquid Reactions.* McGraw-Hill.

# %% [markdown]
# ### Claims and the prose
#
# Every number the README states is checked below against the claim it quotes,
# so a number in the text cannot drift from the computation.

# %%
SAYS = {
    "first_order_max_rel_error": ["6.3e-3"],
    "reaction_layer_cells_at_ha_max": ["2"],
    "excess_reduces_to_first_order": ["5.4e-6"],
    "vkh_max_error_pct": ["2.1"],
    "decoursey_max_error_pct": ["8.7"],
    "vkh_worst_error_pct": ["-2.08"],
    "vkh_worst_ha": ["8.23"],
    "vkh_worst_ha_over_ei": ["1.65"],
    "decoursey_worst_error_pct": ["8.73"],
    "decoursey_worst_ha": ["1.43"],
    "decoursey_worst_ha_over_ei": ["0.0143"],
    "decoursey_pfo_limit_max_error_pct": ["8.81"],
    "decoursey_worst_share_of_pfo_limit": ["0.991"],
    "convergence_order": ["between 2.9 and 3.1"],
    "convergence_order_ha2": ["between 1.9 and 2.1"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
print(f"prose literals: {sum(len(v) for v in SAYS.values())} checked against their claims")

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
