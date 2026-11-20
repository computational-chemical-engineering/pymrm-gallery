# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
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
# title: "Duncan–Toor ternary gas diffusion"
# description: "The two-bulb experiment where Fick's law provably fails: nitrogen diffuses uphill, stops against its own gradient, then reverses."
# categories: [sec:A, struct:S3, struct:S9, tier:T0, data:tier4, phase:gas]
# date: 2026-07-25
# ---
#
# # Duncan–Toor ternary gas diffusion
#
# **Catalog ID:** `A4.9` · **Structures:** `S3` (1D steady BVP), `S9` (implicit
# multicomponent flux) · **Tier:** T0
#
# In a ternary mixture a species can reverse direction while its own gradient
# still points the same way, and then stop dead while that gradient is at its
# largest. Fick's law cannot produce either. This page reproduces the experiment
# that made it concrete, and the Maxwell–Stefan solve that explains it.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.20.3"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A4.9-duncan-toor"  # equals the page directory

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
# Imports. Everything beyond numpy, scipy, pandas and matplotlib is pymrm.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from pymrm.checks import find_roots, observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import (TwoBulbCell, NEWTON_TOL, fick_difference, toor_difference,
                   wilke_deff)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# Fick's law says a species moves down its own concentration gradient. For binary
# mixtures that is essentially always true. For three or more species it is not,
# because each species is dragged by every other one.
#
# Duncan and Toor (1962) built the cleanest possible demonstration: two bulbs of
# gas joined by a thin capillary, held at 35.2 °C and 1 atm. Bulb 1 was charged
# with hydrogen and nitrogen, bulb 2 with nitrogen and carbon dioxide. Hydrogen
# and carbon dioxide then exchange through the capillary, and nitrogen, which
# starts only slightly apart between the two bulbs (the difference is printed in
# "The data"), is dragged along by them.
#
# Four runs were made from the same initial charge, each stopped at a different
# time (1, 2, 4 and 12 h) and analysed by gas chromatography. Read in sequence,
# the nitrogen trajectory shows:
#
# 1. **ordinary diffusion**: nitrogen initially moves down its own small
#    gradient, as Fick would have it;
# 2. **osmotic diffusion**: the gradient reaches zero (marked *point O* in the
#    original figure) but the flux does not stop;
# 3. **uphill (reverse) diffusion**: nitrogen keeps going, now *against* its own
#    gradient, building a difference larger than it started with (the ratio
#    measured is printed in "The data");
# 4. **a diffusion barrier**: the flux finally reaches zero at the moment the
#    gradient is at its largest;
# 5. and only then ordinary downhill diffusion again.
#
# Steps 2 to 4 are impossible under Fick's law with any positive diffusivity.
# The experiment is small, the geometry is simple, and the failure is
# unambiguous, which makes it a clean demonstration of why a multicomponent
# transport model is needed.

# %% [markdown]
# ## The published model
#
# **Maxwell–Stefan.** For an ideal gas mixture of $n$ species at uniform
# temperature and pressure, the driving force on species $i$ is balanced by
# friction with every other species:
#
# $$
# \nabla x_i \;=\; \sum_{j \neq i} \frac{x_i N_j - x_j N_i}{c_t\,\mathrm{D}_{ij}}
# $$
#
# with $x_i$ the mole fraction, $N_i$ the molar flux, $c_t$ the total molar
# concentration, and $\mathrm{D}_{ij}$ the Maxwell–Stefan diffusivity of the
# $i\!-\!j$ pair. The capillary is closed and isobaric, so
#
# $$
# \sum_{i=1}^{n} N_i = 0 .
# $$
#
# **Matrix form.** Eliminating species $n$ gives $n-1$ independent equations that
# can be inverted for the fluxes,
#
# $$
# (N) = -\,c_t\,[B]^{-1}\,(\nabla x), \qquad
# B_{ii} = \frac{x_i}{\mathrm{D}_{in}} + \sum_{k \neq i}\frac{x_k}{\mathrm{D}_{ik}},
# \qquad
# B_{ij} = -x_i\!\left(\frac{1}{\mathrm{D}_{ij}} - \frac{1}{\mathrm{D}_{in}}\right).
# $$
#
# $[B]$ depends on composition, so the flux of each species depends on the
# gradients of *all* of them. That coupling is the entire phenomenon: the
# off-diagonal terms let $N_i$ be nonzero when $\nabla x_i = 0$.
#
# **Capillary and bulbs.** The capillary volume is about 0.4 % of a bulb volume
# (computed in "PyMRM implementation"), so it is quasi-steady: at each instant
# $\nabla\!\cdot\!N_i = 0$, and $N_i$ is uniform along the capillary. The bulbs
# are well mixed:
#
# $$
# V_1 c_t \frac{\mathrm{d}x_{i,1}}{\mathrm{d}t} = -A\,N_i , \qquad
# V_2 c_t \frac{\mathrm{d}x_{i,2}}{\mathrm{d}t} = +A\,N_i .
# $$
#
# | Symbol | Code | Meaning |
# |---|---|---|
# | $x_i$ | `x[:, i]` | mole fraction in the capillary, cell-centred |
# | $x_{i,1}, x_{i,2}$ | `xb1`, `xb2` | bulb compositions (the boundary values) |
# | $N_i$ | `TwoBulbCell.flux` | molar flux at capillary faces |
# | $[B]$ | `ms_b_matrix` | Maxwell–Stefan resistance matrix |
# | $\mathrm{D}_{ij}$ | `D_MS` | MS pair diffusivities |
# | $c_t$ | `c_t` | total molar concentration |
# | $V_1, V_2$, $L/A$ | `VOL_1`, `VOL_2`, `LA_EFF` | bulb volumes, effective length over area |

# %% [markdown]
# ## Parameters and assumptions
#
# **Assumptions:** ideal gas; isothermal and isobaric; capillary quasi-steady
# (justified below); bulbs perfectly mixed; no viscous or thermal diffusion.
#
# All parameters below were checked against the paper. Two are worth singling
# out:
#
# - The paper does **not** use the geometric length/area ratio of the capillary.
#   It measured an *effective* $(L/A)$ from seven binary diffusion runs (its
#   Table 1; both values are in the next cell) because a residual concentration
#   gradient exists inside the bulbs, so the effective diffusion path is longer
#   than the capillary. It states that Figure 2 is computed with the effective
#   value, and this page follows that.
# - The initial compositions are the $t=0$ points of Figure 2 itself, read from
#   the dataset. Nitrogen does **not** start equal in the two bulbs: it starts
#   apart and then crosses over. That crossing is the point of the experiment.

# %%
# Physical constants and cell geometry, Duncan and Toor (1962).
T = 308.35                 # K, water bath, 35.2 +/- 0.01 C
p = 101.325e3              # Pa, 1 atm
R_GAS = 8.314462           # J/(mol K)
c_t = p / (R_GAS * T)      # mol/m3, total molar concentration

LENGTH = 85.9e-3           # m, capillary, 8.59 cm
VOL_1 = 77.99e-6           # m3, bulb 1, charged with H2 + N2
VOL_2 = 78.63e-6           # m3, bulb 2, charged with N2 + CO2
# Effective length/area ratio measured from seven binary runs (Table 1), used for
# Fig. 2 in place of the geometric value. Only the ratio enters the model.
LA_EFF = 258.1e2           # 1/m, effective
LA_GEOM = 251.1e2          # 1/m, geometric (quoted for comparison only)

# Species: 0 = H2, 1 = N2, 2 = CO2 (the closing species of [B]).
SPECIES = ["H$_2$", "N$_2$", "CO$_2$"]
D_MS = np.zeros((3, 3))                    # Maxwell-Stefan pair diffusivities, m2/s
D_MS[0, 1] = D_MS[1, 0] = 83.3e-6          # H2 - N2  (0.833 cm2/s)
D_MS[0, 2] = D_MS[2, 0] = 68.0e-6          # H2 - CO2 (0.680 cm2/s)
D_MS[1, 2] = D_MS[2, 1] = 16.8e-6          # N2 - CO2 (0.168 cm2/s)

T_END = 16 * 3600.0        # s, the span of Fig. 2
N_Z = 40                   # capillary cells
N_STEPS = 400              # backward-Euler steps over T_END
CELL = dict(d_ms=D_MS, length=LENGTH, la_eff=LA_EFF, vol_1=VOL_1, vol_2=VOL_2, c_t=c_t)

print(f"c_t = {c_t:.2f} mol/m3; time step {T_END / N_STEPS / 3600:.2f} h")
print(f"effective L/A {LA_EFF / 100:.1f} 1/cm against geometric {LA_GEOM / 100:.1f} 1/cm: "
      f"{LA_EFF / LA_GEOM - 1:.2%} longer")

# %% [markdown]
# ## The data
#
# Bulb compositions versus time, digitised from Figure 2 of the original article.
# The extraction was programmatic (600 dpi render, morphological opening to
# isolate the markers, axis calibration from the gridlines) and is validated by
# three independent checks recorded in the provenance sidecar: composition
# closure, monotonicity of all six series, and the closed-cell exchange balance.
#
# Two markers of the figure overlap another series and could not be separated;
# they are absent rather than guessed (the count of points is printed below).
# The initial bulb compositions are the $t = 0$ rows, normalised to sum to one.

# %%
data = load_data("duncan-toor-1962-run1.csv", page=PAGE)
meta = load_meta("duncan-toor-1962-run1.csv", page=PAGE)
IDX = {"H2": 0, "N2": 1, "CO2": 2}


def bulb_composition(bulb, time_h):
    """Composition of one bulb at one measured time, from the CSV rows, normalised."""
    rows = data[(data["bulb"] == bulb) & (data["time_h"] == time_h)]
    if len(rows) != 3:
        raise ValueError(f"bulb {bulb} at {time_h} h has {len(rows)} species rows, not 3")
    x = np.zeros(3)
    for _, r in rows.iterrows():
        x[IDX[r["species"]]] = r["x"]
    return x / x.sum()


X1_0, X2_0 = bulb_composition(1, 0), bulb_composition(2, 0)
N_POINTS = len(data)
exp_n2 = data[data["species"] == "N2"].pivot(index="time_h", columns="bulb", values="x").dropna()
exp_dx = (exp_n2[1] - exp_n2[2])                       # measured N2 difference, bulb 1 - bulb 2
EXP_DX_RATIO = float(np.max(np.abs(exp_dx)) / abs(exp_dx.loc[0]))
claims.exact("n_points", N_POINTS, structural=True, note="rows of the CSV")
claims.exact("n2_diff_initial", float(X1_0[1] - X2_0[1]), structural=True,
             note="normalised t = 0 rows of the CSV")
claims.exact("exp_n2_diff_ratio", EXP_DX_RATIO, structural=True,
             note="largest measured |dx_N2| over the t = 0 value, CSV rows")

print(cite_data(meta))
print(f"\n{N_POINTS} points, {data['time_h'].nunique()} times, "
      f"estimated error {meta['acquisition']['estimated_error'].strip()}")
print(f"initial N2 difference between the bulbs: {X1_0[1] - X2_0[1]:+.4f}")
print(f"measured: largest |N2 difference| / initial = {EXP_DX_RATIO:.2f}")
data.head()


# %% [markdown]
# ## PyMRM implementation
#
# The model lives in `model.py` (`TwoBulbCell`). The capillary is discretised
# with `construct_grad` and `construct_div`; the reusable points:
#
# - **The bulb compositions are boundary values that change every step.** The
#   `shapes_d` argument of `construct_grad` returns separate constant matrices
#   that multiply the boundary values, so all operators are assembled once.
# - **The flux is implicit** (structure `S9`): at every face $[B]$ is built from
#   the face composition and a small linear system is solved; `np.linalg.solve`
#   broadcasts over all faces at once.
# - **`NumJac` needs a neighbour stencil.** The function handed to `NumJac` is the
#   whole residual $\nabla\cdot N(x)$, which reads both neighbours of a cell
#   through the face fluxes, so `NumJac((n_z, 2), axes_diagonals=[0],
#   axes_blocks=[-1])`. Every `newton` call checks `success` and raises.
# - **Time.** The bulbs are marched by backward Euler. Each step is a fixed-point
#   iteration on the two bulb compositions, one steady capillary solve per pass,
#   until the bulb fractions change by less than `1e-14`; the march raises if
#   that is not met. The same bulb equations are also integrated by `solve_ivp`
#   (DOP853, `rtol=1e-11`) with the capillary solved at every evaluation: this
#   route is exact in time to its tolerance, gives a continuous trajectory for
#   the root finds of point O and of the barrier, and measures the
#   backward-Euler time error on the same grid.

# %%
ODE_RTOL, ODE_ATOL = 1e-11, 1e-14     # solve_ivp tolerances of the time-exact route


def full(y):
    """Bulb compositions (n_c,) from the solve_ivp state [x1_H2, x1_N2, x2_H2, x2_N2]."""
    return np.append(y[:2], 1.0 - y[:2].sum()), np.append(y[2:], 1.0 - y[2:].sum())


def deviations(ts, h1, h2):
    """|model - measured| at every CSV row; the measured times are whole steps of ts."""
    dt, out = ts[1] - ts[0], []
    for _, r in data.iterrows():
        k = int(round(r["time_h"] * 3600.0 / dt))
        if abs(k * dt - r["time_h"] * 3600.0) > 1e-9 * dt:
            raise ValueError(f"t = {r['time_h']} h is not a step time")
        out.append((h1 if r["bulb"] == 1 else h2)[k, IDX[r["species"]]] - r["x"])
    return np.abs(np.asarray(out))


def regimes(c, sol, t_end_h=T_END / 3600.0):
    """Point O and the diffusion barrier on the time-exact route.

    Point O is the root of dx_N2(t) = x_N2,1 - x_N2,2. The barrier is the root of
    d(dx_N2)/dt, which is the N2 flux term of the bulb equations, evaluated with a
    capillary solve (the analytic derivative, not a difference of the trajectory).
    Both scans cover the whole span; a root count other than one gives nan.
    """
    dx = lambda th: float(sol.sol(th * 3600.0)[1] - sol.sol(th * 3600.0)[3])

    def n2_rate(th):
        r1, r2 = c.bulb_rates(*full(sol.sol(th * 3600.0)))
        return float(r1[1] - r2[1])

    ro = find_roots(dx, 0.0, t_end_h, n_scan=64)
    rb = find_roots(n2_rate, 0.0, t_end_h, n_scan=64)
    t_o = ro["roots"][0] if len(ro["roots"]) == 1 else np.nan
    t_b = rb["roots"][0] if len(rb["roots"]) == 1 else np.nan
    extreme = dx(t_b) if np.isfinite(t_b) else np.nan
    return dict(t_o=t_o, t_b=t_b, extreme=extreme, msg_o=ro["message"], msg_b=rb["message"])


def refined_max(f, lo, hi, n_scan=401):
    """Largest value of f on [lo, hi]: a scan, then a bounded maximisation around the
    best scan point. Raises when the best point is a domain end (at or beyond it)."""
    grid = np.linspace(lo, hi, n_scan)
    vals = np.array([f(g) for g in grid])
    k = int(np.argmax(vals))
    if k in (0, n_scan - 1):
        raise ValueError("maximum at the end of the window")
    res = minimize_scalar(lambda g: -f(g), bounds=(grid[k - 1], grid[k + 1]),
                          method="bounded", options={"xatol": 1e-10})
    return float(max(-res.fun, vals[k]))


def mole_balance(c, h1, h2):
    """Largest change of a species inventory V1 x1 + V2 x2 over the run, over the total."""
    tot0 = c.vol_1 * X1_0 + c.vol_2 * X2_0
    return float(np.max(np.abs(c.vol_1 * h1[-1] + c.vol_2 * h2[-1] - tot0)) / tot0.sum())


def study(n_z=N_Z, nu=0, **over):
    """Every data-facing number of the page for one set of constants.

    The baseline and every constant-perturbing break row run this same function;
    `over` replaces entries of CELL (d_ms, la_eff, vol_2, ...), so a perturbed
    constant enters the solve and all the post-processing.
    """
    c = TwoBulbCell(**{**CELL, **over}, n_z=n_z, nu=nu)
    ts, h1, h2, _ = c.march(X1_0, X2_0, T_END, N_STEPS)
    sol = c.integrate(X1_0, X2_0, T_END, rtol=ODE_RTOL, atol=ODE_ATOL)
    dev = deviations(ts, h1, h2)
    reg = regimes(c, sol)
    geo = (c.d_ms, c.area, c.length, c.vol_1, c.vol_2)
    t_exp = exp_dx.index.values * 3600.0
    k_exp = np.rint(t_exp / (ts[1] - ts[0])).astype(int)
    fick = fick_difference(t_exp, X1_0, X2_0, *geo)[:, 1]
    toor = toor_difference(t_exp, X1_0, X2_0, *geo)[:, 1]

    def toor_gap(th):
        y = sol.sol(th * 3600.0)
        return abs(toor_difference(th * 3600.0, X1_0, X2_0, *geo)[0, 1] - (y[1] - y[3]))

    try:
        toor_max = refined_max(toor_gap, 0.0, T_END / 3600.0)
    except ValueError:                       # largest gap at the end of the span: no interior maximum
        toor_max = np.nan
    return dict(
        mean_abs_deviation=float(dev.mean()), max_abs_deviation=float(dev.max()),
        rmse=float(np.sqrt(np.mean(dev**2))),
        point_o_time_h=reg["t_o"], barrier_time_h=reg["t_b"],
        n2_extreme_ratio=abs(reg["extreme"]) / abs(X1_0[1] - X2_0[1]),
        ms_n2_dev=float(np.mean(np.abs((h1[k_exp, 1] - h2[k_exp, 1]) - exp_dx.values))),
        fick_n2_dev=float(np.mean(np.abs(fick - exp_dx.values))),
        toor_n2_dev=float(np.mean(np.abs(toor - exp_dx.values))),
        toor_max_dev=toor_max, toor_dev_frac=toor_max / abs(reg["extreme"]),
        flux_nonuniformity=c.flux_nonuniformity(X1_0, X2_0),
        mole_balance_error=mole_balance(c, h1, h2),
        capillary_volume_fraction=c.capillary_volume_fraction,
        _reg=reg, _run=(ts, h1, h2), _sol=sol, _cell=c)


BASE = study()
REG, ode, cell = BASE["_reg"], BASE["_sol"], BASE["_cell"]
t_s, x_bulb1, x_bulb2 = BASE["_run"]
t_h = t_s / 3600.0
print(f"backward Euler: {N_STEPS} steps; solve_ivp: {ode.nfev} evaluations, "
      "each a steady capillary solve")
if not (np.isfinite(REG["t_o"]) and np.isfinite(REG["t_b"])):
    raise RuntimeError(f"expected one point O and one barrier: {REG['msg_o']}; {REG['msg_b']}")
print(f"capillary volume / bulb 1 volume = {BASE['capillary_volume_fraction']:.2%}: quasi-steady")

# %% [markdown]
# ## Results
#
# Lines are pymrm, markers are the digitised measurements. Nothing was fitted: the
# diffusivities, geometry and effective $(L/A)$ all come from the paper, and the
# initial compositions are the $t=0$ measurements.

# %%
COL = {"H2": "tab:blue", "N2": "tab:red", "CO2": "tab:green"}
KEY = ["H2", "N2", "CO2"]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
for ax, xb, b in zip(axes, (x_bulb1, x_bulb2), (1, 2)):
    for i, sp in enumerate(KEY):
        ax.plot(t_h, xb[:, i], color=COL[sp], lw=2, label=SPECIES[i])
        sub = data[(data["bulb"] == b) & (data["species"] == sp)]
        ax.plot(sub["time_h"], sub["x"], "o", color=COL[sp], ms=7,
                mfc="white", mew=1.8, zorder=5, label="_nolegend_")
    ax.set_xlabel("time, $t$ [h]")
    ax.set_title(f"Bulb {b}" + (": charged H$_2$ + N$_2$" if b == 1
                                else ": charged N$_2$ + CO$_2$"))
axes[0].set_ylabel("mole fraction, $x$ [-]")
axes[0].legend(frameon=False, loc="center right")
fig.suptitle("Duncan–Toor two-bulb cell: pymrm Maxwell–Stefan (lines) "
             "vs experiment (markers)")
fig.tight_layout()
plt.show()

# %%
# The nitrogen signal is the whole point: its difference between the bulbs.
# Point O and the barrier are roots on the time-exact route (regimes() above).
dx_n2 = x_bulb1[:, 1] - x_bulb2[:, 1]
T_O, T_B, DX_B = REG["t_o"], REG["t_b"], REG["extreme"]
claims.value("point_o_time_h", T_O, rtol=1e-6, unit="h", note="root of dx_N2(t), solve_ivp route")
claims.value("barrier_time_h", T_B, rtol=1e-6, unit="h",
             note="root of the N2 flux term d(dx_N2)/dt, solve_ivp route")
claims.value("n2_extreme_ratio", BASE["n2_extreme_ratio"], rtol=1e-6,
             note="|dx_N2| at the barrier over |dx_N2| at t = 0")
claims.value("capillary_volume_fraction", BASE["capillary_volume_fraction"], rtol=1e-9)

fig, ax = plt.subplots(figsize=(8, 4.6))
ax.plot(t_h, dx_n2, color="tab:red", lw=2.2, label="pymrm")
ax.plot(exp_dx.index, exp_dx.values, "o", color="tab:red", ms=8,
        mfc="white", mew=1.8, zorder=5, label="experiment")
ax.axhline(0.0, color="k", lw=0.8)
ax.plot(T_O, 0.0, "s", color="k", ms=7, zorder=6)
ax.plot(T_B, DX_B, "D", color="k", ms=7, zorder=6)
ax.annotate(f"point O: osmotic diffusion\nzero gradient, non-zero flux\n"
            f"t = {T_O:.2f} h",
            xy=(T_O, 0.0), xytext=(1.6, -0.052),
            arrowprops=dict(arrowstyle="->", lw=1.1), fontsize=9)
ax.annotate(f"diffusion barrier\nlargest gradient, zero flux\n"
            f"t = {T_B:.2f} h",
            xy=(T_B, DX_B), xytext=(8.6, 0.045),
            arrowprops=dict(arrowstyle="->", lw=1.1), fontsize=9)
ax.annotate("uphill (reverse) diffusion", xy=(4.0, 0.005), fontsize=9,
            ha="center", style="italic")
ax.set_xlabel("time, $t$ [h]")
ax.set_ylabel(r"$x_{\mathrm{N_2},1} - x_{\mathrm{N_2},2}$ [-]")
ax.set_title("Nitrogen: four regimes in one experiment")
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
plt.show()

print(f"N2 difference at t = 0      : {X1_0[1] - X2_0[1]:+.4f}")
print(f"crosses zero (point O) at   : {T_O:.2f} h   ({REG['msg_o']} on 0 to 16 h)")
print(f"extreme (diffusion barrier) : {DX_B:+.4f} at t = {T_B:.2f} h   ({REG['msg_b']})")
print(f"ratio |extreme| / |initial| : {BASE['n2_extreme_ratio']:.2f} "
      f"(measured, largest of the sampled times: {EXP_DX_RATIO:.2f})")

# %% [markdown]
# Nitrogen starts apart, closes that gap, and then keeps going, ending up about
# twice as far apart in the *opposite* sense (the printed ratio) before slowly
# relaxing. At the marked barrier its flux is zero while its gradient is at its
# largest.

# %% [markdown]
# ## Validation
#
# The checks are **not** of equal weight.
#
# **1. The comparison with the digitised points of Figure 2 is the only guard on
# the physics.** It is the one number here that moves when a transcribed
# diffusivity, the effective $(L/A)$ or a bulb volume is wrong; the break table
# below shows how far.
#
# **2. Grid refinement** on the time-exact `solve_ivp` route, so the time error is
# removed rather than mixed in. **3. Time-step refinement** of backward Euler on
# the production grid, with its error measured against the time-exact route on
# the same grid, so the shared spatial error cancels. Both bound the numerical
# error and say nothing about the parameters.
#
# **4. Flux uniformity and the closed-cell mole balance are one structural
# identity.** `solve_capillary` drives $\nabla\cdot N$ to the Newton tolerance
# with `construct_div(nu=0)`; on a Cartesian grid that residual *is*
# $(N_{i+1}-N_i)/\Delta z$, so flux uniformity is the converged Newton residual
# divided by a scale. With the flux uniform, $N_0 = N_L$ and the two bulb updates
# cancel, so the mole balance holds whenever uniformity does. Neither depends on
# the diffusivities, $(L/A)$ or the volumes; they catch an unconverged solve and
# a wrong geometry operator, and the break table shows both.

# %%
# 1. Agreement with the measurements (no fitted parameter).
for _k, _what in (("mean_abs_deviation", "mean"), ("max_abs_deviation", "largest"),
                  ("rmse", "root-mean-square")):
    claims.value(_k, BASE[_k], rtol=1e-6, note=f"{_what} |model - measured| over the CSV rows")
print(f"1. |deviation| vs experiment over {N_POINTS} points: mean "
      f"{100 * BASE['mean_abs_deviation']:.2f} mole %, largest "
      f"{100 * BASE['max_abs_deviation']:.2f} mole %, rms {100 * BASE['rmse']:.2f} mole %")
print("   The paper reports 0.45 mole % average deviation of its own Maxwell-Stefan")
print("   predictions and quotes an experimental error of 2.6 mole %; the digitisation")
print("   error stated in the sidecar is part of the figure above.")

# 2. Grid refinement, time removed: the same bulb equations by solve_ivp on
#    nested grids, N2 in bulb 1 at 4, 8, 12 and 16 h.
NZS = (10, 20, 40, 80, 160)
T_OBS_H = np.array([4.0, 8.0, 12.0, 16.0])   # whole steps of every ladder below


def grid_study(face="mean", rtol=ODE_RTOL):
    """Observed spatial orders, the spatial error at N_Z and the round-off margin."""
    sols = {nz: TwoBulbCell(**CELL, n_z=nz, face=face).integrate(
        X1_0, X2_0, T_END, rtol=rtol, atol=1e-3 * rtol) for nz in NZS}
    st = observed_orders(lambda nz: sols[nz].sol(T_OBS_H * 3600.0)[1], NZS)
    orders = np.array(st["orders"])
    q = {nz: sols[nz].sol(t_s)[1] for nz in (N_Z, NZS[-2], NZS[-1])}
    extrap = q[NZS[-1]] + (q[NZS[-1]] - q[NZS[-2]]) / 3.0    # Richardson, theory order 2
    fine = np.abs(sols[NZS[-1]].sol(T_OBS_H * 3600.0)[1] - sols[NZS[-2]].sol(T_OBS_H * 3600.0)[1])
    return dict(orders=orders, worst=float(orders.flat[np.argmax(np.abs(orders - 2.0))]),
                err=float(np.max(np.abs(q[N_Z] - extrap))),
                margin=float(np.min(fine) / (rtol * np.max(np.abs(q[NZS[-1]])))))


GS = grid_study()
claims.bound("grid_observed_order", GS["worst"], min=1.8, max=2.2, basis="theory",
             note="order farthest from 2 over n_z 10 to 160 and N2 in bulb 1 at 4, 8, 12 and 16 h")
claims.bound("grid_order_roundoff_margin", GS["margin"], min=1e3, basis="tolerance",
             note="smallest finest successive difference over rtol * |x|")
claims.value("grid_error_nz40", GS["err"], rtol=1e-3,
             note="max over the step times of |x_N2,1(n_z=40) - Richardson value|")
print(f"\n2. grid ladder n_z = {NZS}, solve_ivp rtol {ODE_RTOL:g}; observed orders "
      "(one per triple, max norm over 4, 8, 12, 16 h):")
print(np.array2string(GS["orders"], precision=3))
print(f"   min {GS['orders'].min():.3f}, max {GS['orders'].max():.3f}; round-off margin "
      f"{GS['margin']:.1e}; spatial error at n_z = {N_Z}: {GS['err']:.2e} in mole fraction")

# 3. Time-step refinement of backward Euler on the production grid.
NSS = (100, 200, 400, 800)


def time_study(theta=1.0, coupling_tol=1e-14):
    """Observed time orders, the time error at N_STEPS and the round-off margin."""
    runs = {}
    for ns in NSS:
        if ns == N_STEPS and theta == 1.0 and coupling_tol == 1e-14:
            runs[ns] = (t_s, x_bulb1)                     # the production march
        else:
            ts_, h1_, _, _ = cell.march(X1_0, X2_0, T_END, ns, theta=theta,
                                        coupling_tol=coupling_tol)
            runs[ns] = (ts_, h1_)

    def q(ns):
        ts_, h1_ = runs[ns]
        k = T_OBS_H * 3600.0 / (ts_[1] - ts_[0])
        if np.max(np.abs(k - np.rint(k))) > 1e-9:
            raise ValueError(f"the ladder times are not whole steps of n_steps = {ns}")
        return h1_[np.rint(k).astype(int), 1]

    st = observed_orders(q, NSS)
    orders = np.array(st["orders"])
    (t4, h4), (_, h8) = runs[N_STEPS], runs[2 * N_STEPS]
    return dict(orders=orders, worst=float(orders.flat[np.argmax(np.abs(orders - 1.0))]),
                err=float(np.max(np.abs(h4[:, 1] - ode.sol(t4)[1]))),
                sens=float(np.max(np.abs(h4[:, 1] - h8[::2, 1]))),
                margin=float(np.min(np.abs(q(N_STEPS) - q(2 * N_STEPS)))
                             / max(coupling_tol, NEWTON_TOL)))


TS = time_study()
claims.bound("timestep_observed_order", TS["worst"], min=0.9, max=1.1, basis="theory",
             note="order farthest from 1 over n_steps 100 to 800 and N2 in bulb 1 at 4, 8, 12 and 16 h")
claims.bound("time_order_roundoff_margin", TS["margin"], min=1e3, basis="tolerance",
             note="smallest finest successive difference over the solver tolerance")
claims.value("timestep_sensitivity", TS["sens"], rtol=1e-4,
             note="max N2 change in bulb 1 on halving dt from the production step")
claims.value("timestep_error_ns400", TS["err"], rtol=1e-4,
             note="max |backward Euler - solve_ivp| for N2 in bulb 1, same grid")
claims.value("time_over_grid_error_ratio", TS["err"] / GS["err"], rtol=1e-3)
claims.value("dev_over_grid_error", BASE["mean_abs_deviation"] / GS["err"], rtol=1e-3)
print(f"\n3. time ladder n_steps = {NSS}, observed orders:")
print(np.array2string(TS["orders"], precision=3))
print(f"   min {TS['orders'].min():.3f}, max {TS['orders'].max():.3f}; round-off margin "
      f"{TS['margin']:.1e}")
print(f"   time error at {N_STEPS} steps, against solve_ivp on the same grid: {TS['err']:.2e}; "
      f"change on halving dt: {TS['sens']:.2e}")
print(f"   time error / spatial error = {TS['err'] / GS['err']:.0f}; mean model-vs-data "
      f"deviation / spatial error = {BASE['mean_abs_deviation'] / GS['err']:.0f}")


def ode_gap(rtol_a=ODE_RTOL, rtol_b=ODE_RTOL / 10):
    """Largest difference of the solve_ivp route between two tolerances, at the step times."""
    ya = ode if rtol_a == ODE_RTOL else cell.integrate(X1_0, X2_0, T_END, rtol=rtol_a, atol=1e-3 * rtol_a)
    yb = cell.integrate(X1_0, X2_0, T_END, rtol=rtol_b, atol=1e-3 * rtol_b)
    return float(np.max(np.abs(ya.sol(t_s) - yb.sol(t_s))))


ODE_GAP = ode_gap()
claims.bound("ode_route_tol_gap", ODE_GAP, max=1e-9, basis="tolerance",
             note="solve_ivp at rtol 1e-11 against 1e-12, all four bulb fractions")
print(f"   solve_ivp route, rtol {ODE_RTOL:g} against {ODE_RTOL / 10:g}: below 1e-9 "
      f"(here {ODE_GAP:.1e})")

# 4. STRUCTURAL: flux uniformity at t = 0 and the closed-cell mole balance.
claims.bound("flux_nonuniformity", BASE["flux_nonuniformity"], max=1e-12, basis="tolerance",
             note="the converged Newton residual over a scale, newton tol 1e-12")
claims.bound("mole_balance_error", BASE["mole_balance_error"], max=1e-13, basis="identity",
             note="N_0 = N_L makes the bulb updates cancel")
print(f"\n4. flux non-uniformity along the capillary: below 1e-12 "
      f"(here {BASE['flux_nonuniformity']:.1e})   [STRUCTURAL]")
print(f"   species mole-balance error over 16 h  : below 1e-13 "
      f"(here {BASE['mole_balance_error']:.1e})   [STRUCTURAL]")

# %% [markdown]
# ### Break table: what each check sees
#
# Five constants were read off the paper and drive everything: three
# Maxwell–Stefan pair diffusivities, the effective $(L/A)$ and the bulb volumes.
# Each is perturbed by one digit, the error that happens when reading a scan, and
# the page's own `study()` runs with nothing reimplemented. Further rows inject
# the error classes the numerical and structural checks are for: an unconverged
# capillary solve, the wrong geometry in `construct_div`, a first-order face
# composition for $[B]$, a second-order time step in place of backward Euler,
# and loose solver tolerances under the refinement ladders.

# %%
def with_d(i, j, value):
    """D_MS with one pair diffusivity replaced."""
    d = D_MS.copy()
    d[i, j] = d[j, i] = value
    return d


CONST_ROWS = [
    ("D(N2-CO2) 16.8 -> 61.8 (digit swap)", dict(d_ms=with_d(1, 2, 61.8e-6))),
    ("D(N2-CO2) 16.8 -> 1.68 (lost decade)", dict(d_ms=with_d(1, 2, 1.68e-6))),
    ("D(H2-N2) 83.3 -> 8.33 (lost decade)", dict(d_ms=with_d(0, 1, 8.33e-6))),
    ("LA_EFF 258.1 -> 285.1 (digit swap)", dict(la_eff=285.1e2)),
    ("vol_2 78.63 -> 87.63 (digit swap)", dict(vol_2=87.63e-6)),
]
ROW_UNCONVERGED = "capillary newton stopped after 1 iteration, unchecked"
ROW_NU = "construct_div(nu=1), cylindrical geometry"
ROW_FACE = "[B] at the left-cell composition on interior faces"
ROW_THETA = "trapezoidal step (theta = 0.5) in the time study"
ROW_COUPLING = "bulb coupling tolerance 1e-6 in the time study"
ROW_ODE = "solve_ivp rtol 1e-6 in the grid study and the second route"
ROW_TOOR_V2 = (f"vol_2 {VOL_2 * 1e6:.2f} -> {CONST_ROWS[4][1]['vol_2'] * 1e6:.2f} "
               "in Toor's closed form only (all D_ij equal)")
claims.define_rows([r[0] for r in CONST_ROWS]
                   + [ROW_UNCONVERGED, ROW_NU, ROW_FACE, ROW_THETA, ROW_COUPLING, ROW_ODE,
                      ROW_TOOR_V2])

STUDY_KEYS = [k for k in BASE if not k.startswith("_")]
INJ = {}
for label, over in CONST_ROWS:
    try:
        INJ[label] = study(**over)
    except Exception as exc:                 # noqa: BLE001 - a failed study is a result
        print(f"  {label}: study raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
        INJ[label] = {k: np.nan for k in STUDY_KEYS}


def flux_count(maxfev=30, check=True):
    """Flux non-uniformity of the five constant rows (counted below 1e-10, far from both round-off and an unconverged solve)."""
    vals = [TwoBulbCell(**{**CELL, **over}, n_z=N_Z).flux_nonuniformity(
        X1_0, X2_0, maxfev=maxfev, check=check) for _, over in CONST_ROWS]
    return vals


FLUX_VALS = flux_count()
claims.count("n_injections_flux_within_bound", FLUX_VALS, 1e-10, above=False, min_margin=0.5,
             note="constant rows whose flux non-uniformity stays below 1e-10")

# What pymrm adds: Fick and Toor's linearisation against the measurements.
for _k in ("ms_n2_dev", "fick_n2_dev", "toor_n2_dev"):
    claims.value(_k, BASE[_k], rtol=1e-6, note="mean |dx_N2 - measured| over the measured times")
claims.value("toor_max_dev", BASE["toor_max_dev"], rtol=1e-5,
             note="max over t of |Toor - solve_ivp| in dx_N2, refined maximisation")
claims.value("toor_dev_frac", BASE["toor_dev_frac"], rtol=1e-5, note="toor_max_dev over |dx_N2| at the barrier")


def toor_gap_equal_d(vol_2_toor=None):
    """Largest |Toor - solve_ivp| in dx_N2 over the run with every D_ij equal to D(N2-CO2).

    Then [B] = I / D whatever the composition, so Toor's frozen [B] is exact and
    the gap is the integrator tolerance. `vol_2_toor` replaces vol_2 in Toor's
    closed form only, a defect the solve_ivp route does not share.
    """
    d = np.where(np.eye(3, dtype=bool), 0.0, D_MS[1, 2])
    c = TwoBulbCell(**{**CELL, "d_ms": d}, n_z=N_Z)
    sol = c.integrate(X1_0, X2_0, T_END, rtol=ODE_RTOL, atol=ODE_ATOL)
    ts = np.linspace(0.0, T_END, 401)
    y = sol.sol(ts)
    v2 = c.vol_2 if vol_2_toor is None else vol_2_toor
    toor = toor_difference(ts, X1_0, X2_0, c.d_ms, c.area, c.length, c.vol_1, v2)[:, 1]
    return float(np.max(np.abs(toor - (y[1] - y[3]))))


TOOR_EQUAL_D = toor_gap_equal_d()
claims.bound("toor_gap_equal_d", TOOR_EQUAL_D, max=1e-10, basis="tolerance",
             note="max over 401 times of |Toor - solve_ivp| in dx_N2, all D_ij equal: [B] constant")


def nu1_mole_balance():
    """Mole balance of a march with nu = 1; its bulb coupling converges slowly, so more passes."""
    c = TwoBulbCell(**CELL, n_z=N_Z, nu=1)
    _, h1, h2, _ = c.march(X1_0, X2_0, T_END, N_STEPS, max_coupling=300)
    return mole_balance(c, h1, h2)


_face = {}


def face_study():
    if not _face:
        _face.update(grid_study(face="upwind"))
    return _face


_theta = {}


def theta_study():
    if not _theta:
        _theta.update(time_study(theta=0.5))
    return _theta


ROWS = [(label, {k: (lambda k=k, label=label: INJ[label][k]) for k in STUDY_KEYS})
        for label, _ in CONST_ROWS]
ROWS += [
    (ROW_UNCONVERGED, {
        "flux_nonuniformity": lambda: cell.flux_nonuniformity(X1_0, X2_0, maxfev=1, check=False),
        "n_injections_flux_within_bound": lambda: int(np.sum(
            np.asarray(flux_count(maxfev=1, check=False)) < 1e-10))}),
    (ROW_NU, {
        "flux_nonuniformity": lambda: TwoBulbCell(**CELL, n_z=N_Z, nu=1).flux_nonuniformity(X1_0, X2_0),
        "mole_balance_error": nu1_mole_balance}),
    (ROW_FACE, {
        "grid_observed_order": lambda: face_study()["worst"],
        "grid_error_nz40": lambda: face_study()["err"],
        "dev_over_grid_error": lambda: BASE["mean_abs_deviation"] / face_study()["err"],
        "time_over_grid_error_ratio": lambda: TS["err"] / face_study()["err"]}),
    (ROW_THETA, {
        "timestep_observed_order": lambda: theta_study()["worst"],
        "timestep_sensitivity": lambda: theta_study()["sens"],
        "timestep_error_ns400": lambda: theta_study()["err"],
        "time_over_grid_error_ratio": lambda: theta_study()["err"] / GS["err"]}),
    (ROW_COUPLING, {
        "time_order_roundoff_margin": lambda: time_study(coupling_tol=1e-6)["margin"]}),
    (ROW_ODE, {
        "grid_order_roundoff_margin": lambda: grid_study(rtol=1e-6)["margin"],
        "ode_route_tol_gap": lambda: ode_gap(1e-6, 1e-7)}),
    (ROW_TOOR_V2, {
        "toor_gap_equal_d": lambda: toor_gap_equal_d(vol_2_toor=CONST_ROWS[4][1]["vol_2"])}),
]

BREAK = []
for row, studies in ROWS:
    for name, fn in studies.items():
        try:
            out = fn()
        except Exception as exc:                 # noqa: BLE001 - a failed study is a result
            print(f"  {row}: {name} raised {type(exc).__name__}: {str(exc)[:70]}; records nan")
            out = np.nan
        caught_ = claims.breaks(row, name, out)
        base_ = claims.metrics[name]["value"]
        if claims.metrics[name]["kind"] == "bound":   # round-off digits differ between machines
            base_, out = f"(here {base_:.1e})", f"(here {float(out):.1e})"
        BREAK.append((row, name, base_, out, caught_))

bt = pd.DataFrame(BREAK, columns=["injected defect", "metric", "baseline", "with defect", "caught"])
pd.set_option("display.max_colwidth", 60)
pd.set_option("display.max_rows", 300)
print(bt.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
idle = [r for r in bt["injected defect"].unique()
        if not bt[(bt["injected defect"] == r) & bt.caught].shape[0]]
print(f"\n{len(ROWS)} rows defined, {bt['injected defect'].nunique()} in the table; "
      f"rows that fail no claim: {idle or 'none'}")
missing = [k for k in claims.metrics if not claims.metrics[k].get("structural")
           and not bt[(bt.metric == k) & bt.caught].shape[0]]
print(f"claims no row fails: {missing or 'none'}")
_swap = INJ[CONST_ROWS[0][0]]
print(f"\nThe D(N2-CO2) digit swap moves the mean deviation from "
      f"{100 * BASE['mean_abs_deviation']:.2f} to {100 * _swap['mean_abs_deviation']:.2f} mole % and "
      f"point O from {T_O:.2f} to {_swap['point_o_time_h']:.2f} h, while flux uniformity stays "
      f"below 1e-12 (here {_swap['flux_nonuniformity']:.1e}).")

# %% [markdown]
# ## What pymrm adds
#
# The Maxwell–Stefan equations and the two-bulb balance are classical, and Duncan
# and Toor compared Maxwell–Stefan predictions with these measurements
# themselves. pymrm supplies the discretisation: constant `construct_grad` and
# `construct_div` operators whose boundary values (the bulb compositions) enter
# through `shapes_d` matrices, a `NumJac` Jacobian for the implicit flux, and a
# checked `newton` solve per capillary state. Against the two simpler models the
# page compares, it adds the following.
#
# **Against Fick's law, a qualitative correction.** With effective (Wilke)
# diffusivities each species relaxes on its own and monotonically towards
# equilibrium, so the nitrogen difference can only shrink towards zero. It can
# never cross zero, let alone end up larger than it started with the opposite
# sign. The experiment does both.

# %%
t_fine = np.linspace(0.0, T_END, 1601)
geo = (D_MS, cell.area, LENGTH, VOL_1, VOL_2)
dx_fick = fick_difference(t_fine, X1_0, X2_0, *geo)
d_eff = wilke_deff(0.5 * (X1_0 + X2_0), D_MS)

fig, ax = plt.subplots(figsize=(8, 4.4))
ax.plot(t_h, dx_n2, color="tab:red", lw=2.2, label="pymrm (Maxwell–Stefan)")
ax.plot(t_fine / 3600.0, dx_fick[:, 1], "--", color="tab:gray", lw=2,
        label="Fick with Wilke effective $D$")
ax.plot(exp_dx.index, exp_dx.values, "o", color="tab:red", ms=8,
        mfc="white", mew=1.8, zorder=5, label="experiment")
ax.axhline(0.0, color="k", lw=0.8)
ax.set_xlabel("time, $t$ [h]")
ax.set_ylabel(r"$x_{\mathrm{N_2},1} - x_{\mathrm{N_2},2}$ [-]")
ax.set_title("Fick's law cannot cross zero; the experiment does")
ax.legend(frameon=False)
fig.tight_layout()
plt.show()

print("Wilke effective diffusivities at the mean initial composition [1e-6 m2/s]: "
      + ", ".join(f"{sp} {v * 1e6:.2f}" for sp, v in zip(IDX, d_eff)))
print(f"mean |error| in the N2 difference over the {len(exp_dx)} measured times: "
      f"Fick {100 * BASE['fick_n2_dev']:.2f} mole %, Maxwell-Stefan "
      f"{100 * BASE['ms_n2_dev']:.2f} mole %")

# %% [markdown]
# **Against Toor's linearisation, a quantitative correction.** The original
# analysis froze $[B]$ at the average composition, which turns the problem into a
# matrix exponential. That made the phenomenon explicable and it captures every
# regime. What it cannot do is let $[B]$ vary along the capillary, where the
# composition swings across the ternary diagram while the pair diffusivities
# differ by a factor of about five. The gap is measured against the time-exact
# route, as the largest difference over the run (a refined maximisation).

# %%
dx_lin = toor_difference(t_fine, X1_0, X2_0, *geo)
fig, ax = plt.subplots(figsize=(8, 4.4))
ax.plot(t_h, dx_n2, color="tab:red", lw=2.2, label="pymrm: composition-dependent $[B]$")
ax.plot(t_fine / 3600.0, dx_lin[:, 1], "--", color="tab:purple", lw=2,
        label="Toor linearisation: $[B]$ frozen at the mean")
ax.plot(exp_dx.index, exp_dx.values, "o", color="tab:red", ms=8,
        mfc="white", mew=1.8, zorder=5, label="experiment")
ax.axhline(0.0, color="k", lw=0.8)
ax.set_xlabel("time, $t$ [h]")
ax.set_ylabel(r"$x_{\mathrm{N_2},1} - x_{\mathrm{N_2},2}$ [-]")
ax.set_title("Effect of resolving the composition dependence of $[B]$")
ax.legend(frameon=False)
fig.tight_layout()
plt.show()

print(f"largest |Toor - full| in the N2 difference: {BASE['toor_max_dev']:.4f}, "
      f"{BASE['toor_dev_frac']:.1%} of |dx_N2| at the barrier")
print(f"with every D_ij equal ([B] constant): largest |Toor - full| below 1e-10 "
      f"(here {TOOR_EQUAL_D:.1e}), so the gap is the composition dependence of [B]")
print(f"mean |error| against the measured N2 difference: linearised "
      f"{100 * BASE['toor_n2_dev']:.2f} mole %, full {100 * BASE['ms_n2_dev']:.2f} mole %")

# %% [markdown]
# ## Reuse
#
# **To adapt this page to your system**, change in order:
#
# | Change | Where |
# |---|---|
# | Different species or more of them | `SPECIES`, `D_MS` and the dataset; `TwoBulbCell` takes any `n_c` |
# | Different cell geometry | `LENGTH`, `VOL_1`, `VOL_2`, `LA_EFF` |
# | Non-ideal mixture | multiply $[B]^{-1}$ by the thermodynamic factor $[\Gamma]$ in `TwoBulbCell.flux` |
# | Porous medium instead of an open capillary | add Knudsen terms to $[B]$: this becomes the dusty gas model (`A4.3`) |
# | Add reaction | add a source term to the residual of `TwoBulbCell.solve_capillary` |
#
# `ms_b_matrix` and `TwoBulbCell.flux` are the reusable core and are independent
# of the geometry: any page tagged `struct:S9` can use the same pattern. Keep the
# flux-uniformity and mole-balance checks only as what they are: a convergence
# and geometry check, not a check on the parameters.
#
# **Related pages:** `A4.2` Maxwell–Stefan vs Fick · `A4.3` dusty gas model ·
# `A4.8` Stefan tube · `B1.11` multicomponent pellet · `H1.9` zeolite membrane
# permeation.
#
# ## References
#
# Duncan, J. B. and Toor, H. L. (1962). *An experimental study of three component
# gas diffusion.* AIChE Journal **8**(1), 38–41.
# [doi:10.1002/aic.690080112](https://doi.org/10.1002/aic.690080112)
#
# Krishna, R. and Wesselingh, J. A. (1997). *The Maxwell–Stefan approach to mass
# transfer.* Chemical Engineering Science **52**(6), 861–911.
#
# Taylor, R. and Krishna, R. (1993). *Multicomponent Mass Transfer.* Wiley.

# %%
# Prose literals of this notebook, README.md, meta.yaml and models.yaml, each
# checked against its claim; then the claims report.
SAYS = {
    "capillary_volume_fraction": ["0.4 %"],
    "mean_abs_deviation": ["0.59 %"], "max_abs_deviation": ["1.79 %"],
    "point_o_time_h": ["0.96"], "barrier_time_h": ["7.50"], "n2_extreme_ratio": ["2.0"],
    "n_points": ["28"],
    "fick_n2_dev": ["9.1 %"], "ms_n2_dev": ["1.4 %"], "toor_dev_frac": ["5.0 %"],
    "toor_gap_equal_d": ["below 1e-10"], "exp_n2_diff_ratio": ["between 2.2 and 2.3"],
    "n_injections_flux_within_bound": ["5"],
    "grid_observed_order": ["between 1.8 and 2.2"],
    "timestep_observed_order": ["between 0.9 and 1.1"],
    "grid_error_nz40": ["1.6e-6"], "timestep_error_ns400": ["3.0e-4"],
    "time_over_grid_error_ratio": ["183"], "dev_over_grid_error": ["3.6e3"],
    "flux_nonuniformity": ["below 1e-12"], "mole_balance_error": ["below 1e-13"],
}
for _k, _lits in SAYS.items():
    for _lit in _lits:
        claims.says(_lit, _k)
print(f"prose literals: {sum(len(v) for v in SAYS.values())} checked against their claims with claims.says")
claims.report();   # the semicolon keeps the notebook from echoing the payload
