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
# ---

# %% [markdown]
# ---
# title: "Gallucci et al. (2013): membrane or bed-to-wall film, which limits the hydrogen flux"
# description: "The review's statements on who limits the H2 flux through a Pd membrane (reactor), the membrane or the bed-to-wall film, put on one scale: the flux elasticity of a Sieverts membrane behind a film, decided from the review's own figures and printed numbers, with a pymrm membrane module for the area statement."
# categories: [sec:H, struct:S7, tier:T1, data:tier4, phase:gas]
# date: 2026-11-28
# ---
#
# # Gallucci et al. (2013): membrane or bed-to-wall film, which limits the hydrogen flux
#
# **Catalog ID:** `H1.6` · **Structures:** `S7` (retentate and permeate coupled by the membrane flux) · **Tier:** T1
#
# A Pd membrane in a packed or fluidized bed sees hydrogen through two resistances in series:
# the bed-to-wall film that brings H2 from the catalyst bed to the wall, and the membrane
# itself (Sieverts' law). Gallucci, Fernandez, Corengia and van Sint Annaland (2013) review
# membranes and membrane reactors for hydrogen production. On pp. 54 to 57 they state, with
# two reproduced figures and a few printed numbers, that thin membranes made the film the
# limiting step. The review prints no model. This page puts those statements on one scale,
# the **flux elasticity**: the fraction of a small relative change in the flux that a
# relative change in the membrane permeance (or the film conductance) produces. It decides
# each statement from the review's own numbers where they suffice, gives the deciding
# number where they do not, and solves a membrane module with pymrm for the statement on
# membrane area.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "main"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "H1.6-gallucci-membrane-film-limitation"  # equals the page directory

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
# Imports. pymrm contributes the module solve (operators, BlockLayout, newton) and the
# checks (find_roots, observed_orders); the figure conversions are closed forms.
import contextlib
import io
import math
import re
from dataclasses import dataclass, field, replace

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from pymrm.checks import find_roots, observed_orders
from scipy.optimize import minimize_scalar
from gallery_utils import Claims, assert_check_values, cite_data, load_data, load_meta
from model import (MembraneModule, a_film, area_elasticity, centre_depletion, cupmix_depletion,
                   depletion_vacuum, elasticities, eps_k_from_eta, eps_q_bound, eps_q_vacuum,
                   film_share, fit_profile, lambda_star, local_flux, module_number,
                   n_membrane_vacuum, pinch, recovery_at, share_elasticities, theta_star,
                   wall_factor, wall_root, depletion_with_permeate, theta_star_linear)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])


def fill(text, **v):
    """Fill ``<<name:fmt>>`` in ``text`` from ``v``."""
    def sub(m):
        name, _, fmt = m.group(1).partition(":")
        return format(v[name], fmt)
    return re.sub(r"<<([^<>]+)>>", sub, text)


def say(text, **v):
    """Prose whose numbers are written by code, so they cannot drift."""
    display(Markdown(fill(text, **v)))

# %% [markdown]
# ## Background
#
# A thin Pd film permeates hydrogen fast enough that the gas next to it runs short of H2: the
# membrane then sees a lower H2 partial pressure than the bed. The review states this for
# fixed-bed membrane reactors (p. 54, quoted from the page image):
#
# > "as long as the hydrogen flux through the membrane is a limiting step, the effect of
# > external mass transfer limitations such as the limitations to hydrogen transport between
# > the bulk of the catalytic bed (where hydrogen is produced) and the membrane wall (where
# > hydrogen is recovered) can be neglected. However, by increasing the membrane flux, the
# > external mass transfer limitations became limiting and determine the extent of membrane
# > area." (S3)
#
# It backs this with two figures. Fig. 15 (p. 55, reproduced from Tiemersma et al., 2006)
# shows radial H2 profiles at three membrane permeabilities:
#
# > "In their paper the authors claimed that with the membrane available in 2006 the effect of
# > mass transfer limitations was quite limited. However, by increasing their experimental
# > membrane flux (increase foreseen in a couple of years of membrane developments) the
# > external mass transfer limitation could become a limiting step (as indicated in the
# > figure)." (S1, p. 55)
#
# Fig. 16 (p. 56, from Peters et al., 2008) splits the H2 pressure drop into depletion,
# concentration polarisation, Pd membrane and support:
#
# > "As indicated in Fig. 16, it is evident that the extent of bed-to-wall mass transfer
# > limitations (also known as concentration polarization) is the limiting factor for the
# > permeation of hydrogen through the membrane and thus it determines the amount of membrane
# > area required for a given hydrogen recovery." (S2, p. 55)
#
# For fluidized beds the review reports an effectiveness factor on the Sieverts flux (p. 57):
# "The authors suggested that this effectiveness factor could possibly be due fine dust
# generated in the bed or to other mass transfer resistances. [...] could in principle be
# higher than 0.9 (as high as 0.98 according with Chen et al., 2007)." (S5)
#
# The statements are qualitative. "Limiting" needs a scale; the share of the pressure drop is
# the one Fig. 16 draws, but for a Sieverts membrane it is not the share of control over the
# flux. This page uses the flux elasticity, says what number makes each statement true or
# false, and computes that number where the review prints enough to do so.
#
# Not covered here: reaction and kinetics, the pressure drop and particle size (p. 55),
# temperature control, micro-reactors, the other tables (counts of papers), sweep gas. The
# exponent of Sieverts' law is H1.1's subject; concentration polarisation along a laminar
# tube is H1.3's; the equilibrium shift is H1.4's.

# %% [markdown]
# ## The published model
#
# The review prints no equation; it names "Sieverts' equation" once (p. 57). The page's model
# is the standard resistances-in-series pair, scaled by the feed H2 partial pressure
# $p_{in}$:
#
# | id | relation | symbol in code |
# |---|---|---|
# | E1 | membrane $J = Q_s(\sqrt{p_w} - \sqrt{p_p})$ (Sieverts, $n = 1/2$) | `local_flux` |
# | E2 | film $J = k_p(p_b - p_w)$ (linear, isothermal, no high-flux factor) | `local_flux` |
# | E4 | $u^2 + \beta u = \pi + \beta\sqrt{\theta}$, $u = \sqrt{p_w/p_{in}}$, $\pi = p_b/p_{in}$, $\theta = p_p/p_{in}$, $\beta = Q_s/(k_p\sqrt{p_{in}})$ | `wall_root` |
# | E5 | $j = J/(Q_s\sqrt{p_{in}}) = (\pi - \theta)/(\beta + u + \sqrt\theta) = u - \sqrt\theta = (\pi - u^2)/\beta$ | `local_flux(form=...)` |
# | E6 | $\varepsilon_Q = \partial\ln J/\partial\ln Q_s = 2u/(2u + \beta)$, $\varepsilon_k = \beta/(2u + \beta)$, $\varepsilon_Q + \varepsilon_k = 1$ | `elasticities` |
# | E7 | vacuum permeate, $x^2 = p_w/p_b$: $\varepsilon_Q = 2x^2/(1 + x^2)$; $\varepsilon_Q > 1/2$ iff $p_p/p_b < \theta^*_b = (3x^2 - 1)^2/(4x^2)$ | `eps_q_vacuum`, `theta_star` |
# | E8 | vacuum, fixed bulk, permeability times $\lambda$: $\beta_b = \lambda\beta_{b1}$, $\beta_{b1} = d_1/\sqrt{1 - d_1}$; film limiting at $\lambda^* = (2/\sqrt3)/\beta_{b1}$ | `lambda_star`, `depletion_vacuum` |
# | E9 | shares $c$ (film), $m$ (membrane), $\sigma$ (support) of the total drop: $\varepsilon_Q \le 1/(1 + n c/m) \le 1/(1 + c/(2m))$ for any interface pressure and support share | `share_elasticities`, `eps_q_bound` |
# | E10 | effectiveness factor $\eta$ as a series film: $x = y + \eta(1 - y)$, $\beta_b = (1 - x^2)/(x - y)$, $\varepsilon_k = \beta_b/(2x + \beta_b)$ | `eps_k_from_eta` |
# | E11 | module: $df/d\xi = -N j(\pi(f))$, $dg/d\xi = +N j$, $\pi(f) = f(1 + \phi)/(f + \phi)$, $N = Q_s\sqrt{p_{in}}A/F_{H_2,in}$, film units $a = N/\beta$ | `MembraneModule` |
# | E12 | $N(R^*) = \int_{1 - R^*}^1 df/j$ | `module_number` |
# | E13 | limits: $\beta \to 0$ (membrane), $\beta \to \infty$ ($a_{film}$, closed form), pinch $R_{max} = 1 - \theta\phi/(1 + \phi - \theta)$ | `n_membrane_vacuum`, `a_film`, `pinch` |
# | E14 | area elasticity $e_A = -d\ln a/d\ln\beta = \int(\varepsilon_Q/j)df / \int(1/j)df$ | `area_elasticity` |
# | E15 | Fig. 15 profile $1 - D r^2$ normalised at the first marker $r_0 = 1/24$; wall over last-marker depletion $(1 - r_0^2)/(r_l^2 - r_0^2)$, $r_l = 23/24$; cup-mix depletion $d_{1D} = (D/2)/(1 - D/2)$ | `wall_factor`, `centre_depletion`, `cupmix_depletion` |
#
# "Limiting" on this page: the film limits when $\varepsilon_k > 1/2$, i.e. when a small relative
# increase of the film conductance raises the flux more than the same increase of the permeance. Under vacuum permeate
# that happens only when the film takes two thirds of the drop (E7), not half.

# %% [markdown]
# ## Parameters and assumptions
#
# - **Test condition TC** (module): feed H2 fraction and pressure of the review's WGS test
#   (p. 46), permeate pure H2 at 1 bar (**assumption**: the review prints no permeate
#   pressure), $\beta = 1$, $N = 2$, target recovery $R^* = 0.8$. Isothermal, isobaric, plug
#   flow, no reaction.
# - **Fig. 15**: the plotted ratio is a weight fraction; the page uses it as the
#   partial-pressure ratio. With H2 depleted at the wall the mean molar mass rises there, so
#   the weight-fraction depletion overstates the partial-pressure depletion; every
#   bound and threshold of R2 is worded so that a smaller depletion only strengthens it. The
#   permeate pressure, bulk pressure and permeability are not printed.
# - **Fig. 16**: only the ratio of the film and membrane shares enters (E9). The categories
#   are not defined in the review; the source is a permeation cell with sweep gas
#   (Caravella et al., 2009, p. 623), so "bed-to-wall" is the review's label.
# - **Sieverts exponent** $n = 1/2$ (H1.1 owns $n$); $n = 1$ is shown as the sensitivity.

# %%
PN = load_data("gallucci-2013-printed-numbers.csv", page=PAGE, dtype=str, keep_default_na=False)
T2 = load_data("gallucci-2013-table2.csv", page=PAGE, dtype=str, keep_default_na=False)
RD15 = load_data("gallucci-2013-fig15-reads.csv", page=PAGE)
RD16 = load_data("gallucci-2013-fig16-reads.csv", page=PAGE)
for f_ in ("gallucci-2013-printed-numbers.csv", "gallucci-2013-table2.csv",
           "gallucci-2013-fig15-reads.csv", "gallucci-2013-fig16-reads.csv"):
    print(cite_data(load_meta(f_, page=PAGE)))


def read_u(name):
    """The read uncertainty ``u = ...`` stated in a figure sidecar."""
    return float(re.search(r"u = ([0-9.]+)", load_meta(name, page=PAGE)["acquisition"]["estimated_error"]).group(1))


def printed(pn, item):
    """A printed entry, as the string printed."""
    return pn.loc[pn["item"] == item, "printed"].iloc[0]


def num(pn, item):
    """A printed number as a float."""
    return float(printed(pn, item))


U15 = read_u("gallucci-2013-fig15-reads.csv")     # per point, Fig. 15 sidecar
U16 = read_u("gallucci-2013-fig16-reads.csv")     # per boundary, Fig. 16 sidecar
R0, RL = 1.0 / 24.0, 23.0 / 24.0                  # first and last marker, r/R = (i - 1/2)/12
THETA_EX = 0.8                                    # permeate-to-bed H2 ratio of the worked example (R2)
P_PERM_BAR = 1.0                                  # bar, permeate pure H2 (page assumption)
BETA_TC, N_TC, R_STAR = 1.0, 2.0, 0.8             # TC: conductance ratio, membrane number, recovery
BETAS = np.logspace(-2, 2, 41)                    # R4 sweep
NZ0 = 100                                         # coarsest of the three module levels (x 2, x 4)
say("TC: H2 fraction <<y>> % at <<p>> bar (p. 46), permeate <<pp:g>> bar (assumed), so "
    "$\\theta$ = <<th:.5f>> and $\\phi$ = <<ph:.5f>>. Read uncertainties: Fig. 15 <<u15>> per point, "
    "Fig. 16 <<u16>> per boundary (sidecars).",
    y=printed(PN, "p46_wgs_H2_pct"), p=printed(PN, "p46_wgs_P_bar"), pp=P_PERM_BAR,
    th=P_PERM_BAR / (num(PN, "p46_wgs_H2_pct") / 100 * num(PN, "p46_wgs_P_bar")),
    ph=(100 - num(PN, "p46_wgs_H2_pct")) / num(PN, "p46_wgs_H2_pct"), u15=U15, u16=U16)

# %% [markdown]
# ## The data
#
# Four files, each with a provenance sidecar. **Fig. 15 reads** (digitised from the 200 ppi
# embedded image, $u$ per point): the three curves at $z/L = 0.10$; markers 1 to 4 of the
# two triangle curves are hidden behind or merged with the squares and are not read.
# **Fig. 16 reads**: the segment boundaries of the four stacked bars ($u$ per boundary).
# **Table 2** (p. 53): every entry as the printed string. **Printed numbers**: p. 46 (SINTEF
# permeances, WGS feed), the Table 2 footnote, p. 57 ($\eta$). The `check_values.csv` rows
# are computed by a second route from the same equations, not read from the review.

# %% [markdown]
# ## PyMRM implementation
#
# `model.MembraneModule` solves E11 on $n_z$ cells in $\xi = A/A_{tot}$.
#
# - **Unknowns**: retentate $f$ and permeate $g$, both $(n_z, 1)$, side by side on the field
#   axis with `pymrm.BlockLayout`; all order one.
# - **Operators**, assembled once: `construct_convflux_upwind` (unit velocity, flow-rate form)
#   and `construct_div(nu=0)` per stream; first-order upwind, so each cell is backward Euler
#   in $\xi$.
# - **Boundary conditions**, outward normal: inlet $f = 1$, $g = 0$ (`{a: 0, b: 1, d: 1}` and
#   `d: 0`); outlet `{"outflow": True}`; outlet values read with `compute_boundary_values`.
# - **Source**: $N j(f)$ leaves the retentate and enters the permeate in the same cell; its
#   exact derivative $N\,dj/df$ enters through `construct_coefficient_matrix`, placed in the
#   two blocks; the Jacobian is sparse and lower block-bidiagonal.
# - **Area for a recovery**: $N$ is one bordering unknown with the row $g(\xi = 1) = R^*$,
#   solved together with the cells by `pymrm.newton`; the same Jacobian then gives the exact
#   discrete $dN/d\beta$, hence $e_A$ on the grid.
# - **Every module number is extrapolated** over $n_z$ = 100, 200, 400 with
#   `pymrm.checks.observed_orders`; thresholds are roots from `pymrm.checks.find_roots`.
# - **Second route** (`module_number`, `area_elasticity`, `recovery_at`): scipy quadrature of
#   E12 and E14, the closed forms of E13.

# %%
@dataclass(frozen=True)
class Cfg:
    """Everything the study reads; every break row changes one field."""
    pn: object = field(default_factory=lambda: PN)
    t2: object = field(default_factory=lambda: T2)
    rd15: object = field(default_factory=lambda: RD15)
    rd16: object = field(default_factory=lambda: RD16)
    k15: float = 2.0             # exponent of the Fig. 15 radial profile (E15)
    n_mem: float = 0.5           # Sieverts exponent in the Fig. 16 bound (E9)
    u15: float = U15             # Fig. 15 read uncertainty per point
    u16: float = U16             # Fig. 16 read uncertainty per boundary
    p_perm: float = P_PERM_BAR   # bar, permeate (assumption)
    beta_tc: float = BETA_TC
    n_tc: float = N_TC
    r_star: float = R_STAR
    nz0: int = NZ0
    source: str = "cell"         # module source term ("trapezoid" only as a defect)
    beta_small: float = 1e-6     # V2: the general route near beta -> 0
    beta_large: float = 1e6      # V2: the general route near beta -> infinity
    n_large: float = 50.0        # V2: the module near the pinch
    v4_read: str = "outlet"      # V4: the read of the ladder ("inlet" only as a defect)
    ea_weight: str = "elasticity"  # E14 weight of the second route ("share" only as a defect)
    n_route2: float = 0.5        # membrane exponent of the second route in V5 (a defect when not 1/2)


C0 = Cfg()
CACHE = {}


def tc(c):
    """Test condition from the printed WGS feed and the assumed permeate."""
    y = num(c.pn, "p46_wgs_H2_pct") / 100.0
    p_in = y * num(c.pn, "p46_wgs_P_bar")
    return c.p_perm / p_in, (1.0 - y) / y


def extrap(fun, nz0):
    """Richardson over three levels ``nz0 x (1, 2, 4)``: (extrapolated, error estimate, order)."""
    o = observed_orders(lambda n: np.asarray(fun(n), dtype=float), (nz0, 2 * nz0, 4 * nz0))
    return o["extrapolated"], o["error_estimate"], o["orders"][-1]


def module(c, beta, nz):
    theta, phi = tc(c)
    return MembraneModule(theta, phi, beta, nz, source=c.source)


def area_solve(c, beta, nz):
    """(N, e_A) on one grid at recovery ``c.r_star``, cached."""
    theta, phi = tc(c)
    key = ("area", theta, phi, beta, nz, c.r_star, c.source)
    if key not in CACHE:
        m = module(c, beta, nz).solve_area(c.r_star)
        CACHE[key] = (m.N, m.area_elasticity())
    return CACHE[key]


def recovery_solve(c, N, nz, beta=None):
    theta, phi = tc(c)
    beta = c.beta_tc if beta is None else beta
    key = ("rec", theta, phi, beta, nz, N, c.source)
    if key not in CACHE:
        m = module(c, beta, nz).solve_recovery(N)
        CACHE[key] = (m.recovery(), m.balance_error(), m)
    return CACHE[key]


def one_root(f, lo, hi, log=True, n_scan=24):
    hit = find_roots(f, lo, hi, n_scan=n_scan, log=log)
    if len(hit["roots"]) != 1:
        raise ValueError(f"expected one root: {hit['message']}")
    return hit["roots"][0]


# ---------------------------------------------------------------- R1, R2: Fig. 15
CURVES = ("experimental", "doubled", "quadrupled")
TAG = {"experimental": "1x", "doubled": "2x", "quadrupled": "4x"}


def f15_chain(v_last, k):
    """From one last-marker read to (D, d_1D) of the profile (E15)."""
    D = centre_depletion(1.0 - v_last, R0, RL, k)
    return D, cupmix_depletion(D, k)


def part_f15(c, M):
    rd = c.rd15
    r = (rd["cell"].to_numpy() - 0.5) / 12.0
    resid = []
    for cv in CURVES:
        sel = (rd["curve"] == cv).to_numpy() & (rd["cell"].to_numpy() > 1)
        D, res = fit_profile(r[sel], rd["rel_H2_wt_fraction"].to_numpy()[sel], R0, c.k15)
        M[f"r1_D_fit_{TAG[cv]}"] = D
        resid.append(res / c.u15)
    resid = np.concatenate(resid)
    M["_r1_resid_u"] = resid
    M["r1_rms_u"] = float(np.sqrt(np.mean(resid ** 2)))
    M["r1_max_u"] = float(np.max(np.abs(resid)))
    M["_r1_n_read"] = len(resid)
    M["r1_wall_factor"] = wall_factor(R0, RL, c.k15)
    last = {cv: float(rd.loc[(rd["curve"] == cv) & (rd["cell"] == 12), "rel_H2_wt_fraction"].iloc[0])
            for cv in CURVES}
    j_ref = {"2D": 0, "1D": 1}
    for cv in CURVES:
        t = TAG[cv]
        D, d1 = f15_chain(last[cv], c.k15)
        M[f"r1_D_last_{t}"], M[f"r1_d1D_{t}"] = D, d1
        M[f"r1_D_gap_u_{t}"] = (D - M[f"r1_D_fit_{t}"]) / (c.u15 * wall_factor(R0, RL, c.k15))
        for ref, d in (("2D", D), ("1D", d1)):
            x2 = 1.0 - d
            M[f"r2_epsQ_vac_{ref}_{t}"] = eps_q_vacuum(x2)
            M[f"r2_theta_star_{ref}_{t}"] = theta_star(x2)
            if t != "2x":
                M[f"r2_theta_star_n1_{ref}_{t}"] = theta_star_linear(x2)
                hh = 1e-6
                M[f"r2_theta_star_u_{ref}_{t}"] = abs(
                    theta_star(1.0 - f15_chain(last[cv] + hh, c.k15)[j_ref[ref]])
                    - theta_star(1.0 - f15_chain(last[cv] - hh, c.k15)[j_ref[ref]])) / (2 * hh) * c.u15
    # proportionality from the 1x curve (E8), each reference
    h = 1e-6
    for ref, j in (("2D", 0), ("1D", 1)):
        d1 = f15_chain(last["experimental"], c.k15)[j]
        M[f"r2_lambda_star_{ref}"] = lambda_star(d1)
        for cv, lam in (("doubled", 2.0), ("quadrupled", 4.0)):
            def pred(v1):
                dd = f15_chain(v1, c.k15)[j]
                return depletion_vacuum(lam * dd / math.sqrt(1.0 - dd))
            d_obs = f15_chain(last[cv], c.k15)[j]
            dp = pred(last["experimental"])
            s_pred = abs(pred(last["experimental"] + h) - pred(last["experimental"] - h)) / (2 * h) * c.u15
            s_obs = abs(f15_chain(last[cv] + h, c.k15)[j] - f15_chain(last[cv] - h, c.k15)[j]) / (2 * h) * c.u15
            M[f"r2_pred_{ref}_{TAG[cv]}"] = dp
            M[f"r2_gap_units_{ref}_{TAG[cv]}"] = (d_obs - dp) / math.hypot(s_pred, s_obs)
            if ref == "2D":           # the same prediction with the permeate of the 1x deciding ratio
                th_s = M["r2_theta_star_2D_1x"]

                def pred_t(v1, lam=lam):
                    return depletion_with_permeate(f15_chain(v1, c.k15)[0], lam, th_s)
                dpt = pred_t(last["experimental"])
                s_t = abs(pred_t(last["experimental"] + h) - pred_t(last["experimental"] - h)) / (2 * h) * c.u15
                gap_t = (d_obs - dpt) / math.hypot(s_t, s_obs)
                if cv == "quadrupled":
                    M["r2_pred_2D_4x_thstar"], M["r2_gap_units_2D_4x_thstar"] = dpt, gap_t
                else:
                    M["r2_gap_units_2D_2x_thstar"] = gap_t
    # the film of the 1x curve at a permeate ratio below the deciding one (centre reference)
    u1 = math.sqrt(1.0 - M["r1_D_last_1x"])
    beta1 = (1.0 - u1 * u1) / (u1 - math.sqrt(THETA_EX))
    M["r2_ex_epsk_2D_1x"] = beta1 / (2.0 * u1 + beta1)
    M["r2_ex_flux_2D_1x"] = (u1 - math.sqrt(THETA_EX)) / (1.0 - math.sqrt(THETA_EX))
    # direct control from the curves: at a fixed centre value the wall flux scales as D
    def sD(cv):
        return abs(f15_chain(last[cv] + h, c.k15)[0] - f15_chain(last[cv] - h, c.k15)[0]) / (2 * h) * c.u15
    D1, D4 = M["r1_D_last_1x"], M["r1_D_last_4x"]
    M["r2_fd_control_1x_4x"] = math.log(D4 / D1) / math.log(4.0)
    M["r2_fd_control_u"] = math.hypot(sD("quadrupled") / D4, sD("experimental") / D1) / math.log(4.0)
    M["r2_centre_ratio_half"] = 4.0 ** 0.5 / (D4 / D1)      # centre(4x)/centre(1x) for a control of 1/2


# ---------------------------------------------------------------- R3: Fig. 16
PRESSURES = (10, 15, 20, 26)


def f16_bound(dep, cp, pd, n):
    return eps_q_bound(cp - dep, pd - cp, n)


def epsk_min(cc, mm, sigma, n):
    """Minimum film elasticity over the interface ratio ``s`` in ``[sigma/(sigma + m), 1)``
    (permeate H2 pressure not negative): a scan refined by a bounded minimisation."""
    def f(s):
        return share_elasticities(cc, mm, sigma, s, n)[1]
    ss = np.linspace(sigma / (sigma + mm), 1.0 - 1e-6, 2001)
    v = np.array([f(s) for s in ss])
    i = int(np.argmin(v))
    r = minimize_scalar(f, bounds=(ss[max(i - 1, 0)], ss[min(i + 1, len(ss) - 1)]), method="bounded",
                        options={"xatol": 1e-12})
    return float(min(r.fun, v[i]))


def part_f16(c, M):
    margins = []
    for P in PRESSURES:
        row = c.rd16.loc[c.rd16["P_feed_bar"] == P].iloc[0]
        dep, cp, pd = float(row["depletion_top"]), float(row["cp_top"]), float(row["pd_top"])
        cc, mm = cp - dep, pd - cp
        b = f16_bound(dep, cp, pd, c.n_mem)
        h = 1e-6
        grad = [(f16_bound(*(np.array([dep, cp, pd]) + h * e), c.n_mem)
                 - f16_bound(*(np.array([dep, cp, pd]) - h * e), c.n_mem)) / (2 * h) for e in np.eye(3)]
        ub = c.u16 * math.sqrt(sum(g * g for g in grad))
        M[f"r3_c_{P}"], M[f"r3_m_{P}"], M[f"r3_cm_{P}"] = cc, mm, cc / mm
        M[f"r3_bound_{P}"] = b
        M[f"r3_bound_n1_{P}"] = eps_q_bound(cc, mm, 1.0)
        M[f"r3_ubound_{P}"] = ub
        M[f"r3_epsk_min_{P}"] = epsk_min(cc, mm, float(row["support_top"]) - pd, c.n_mem)
        M[f"r3_margin_u_{P}"] = (0.5 - b) / ub
        margins.append((0.5 - b) / ub)
    M["_r3_margins"] = np.array(margins)
    M["r3_n_below_half"] = int(np.sum(np.array(margins) > 0))


# ---------------------------------------------------------------- R6 input: p. 57
def part_p57(c, M):
    for item, t in (("p57_eta_low", "lo"), ("p57_eta_high", "hi")):
        eta = num(c.pn, item)
        M[f"r6_epsk_vac_{t}"] = eps_k_from_eta(eta, 0.0)
        M[f"r6_epsk_y1_{t}"] = 1.0 - eta          # the y -> 1 limit of E10, closed form
        M[f"_r6_epsk_y05_{t}"] = eps_k_from_eta(eta, 0.5)


# ---------------------------------------------------------------- R5: Table 2 and p. 46
READINGS = ("strict", "two", "one")


def half_unit(s, reading):
    """Half a unit of the last significant digit of a printed entry.

    ``"strict"``: every printed digit is significant (130 means 129.5 to 130.5).
    ``"two"``: trailing zeros of an integer are placeholders, but every entry carries at
    least two significant figures (130 and 1,500 have half-units 5 and 50, 20 has 0.5).
    ``"one"``: every trailing zero is a placeholder (20 means 15 to 25).
    """
    s = s.replace(",", "").replace("~", "")
    if "." in s or reading == "strict":
        return 0.5 * 10.0 ** -(len(s.split(".")[1]) if "." in s else 0)
    zeros = len(s) - len(s.rstrip("0"))
    if reading == "two":
        zeros = min(zeros, max(len(s) - 2, 0))
    return 0.5 * 10.0 ** zeros


def interval(s, reading):
    """The interval a printed entry stands for (a range ``a–b`` widened by half units)."""
    s = s.replace(",", "").replace("~", "")
    if "–" in s:
        a, b = s.split("–")
        return float(a) - half_unit(a, reading), float(b) + half_unit(b, reading)
    h = half_unit(s, reading)
    return float(s) - h, float(s) + h


def t2_rows(t2):
    """Rows printing thickness, permeance and permeability."""
    ok = [(r["thickness_um"] != "–") and (r["permeance_1e-8"] != "–") and (r["permeability_1e-13"] != "–")
          for _, r in t2.iterrows()]
    return t2.loc[ok]


def t2_consistent(row, reading):
    """Permeance (1e-8) x thickness (1e-6 m) / 10 against permeability (1e-13), with half units."""
    t0, t1 = interval(row["thickness_um"], reading)
    q0, q1 = interval(row["permeance_1e-8"], reading)
    p0, p1 = interval(row["permeability_1e-13"], reading)
    return (t0 * q0 / 10.0 <= p1) and (t1 * q1 / 10.0 >= p0)


def part_t2(c, M):
    rows = t2_rows(c.t2)
    M["_r5_n_rows"] = len(rows)
    for rd_ in READINGS:
        M[f"r5_n_consistent_{rd_}"] = int(sum(t2_consistent(r, rd_) for _, r in rows.iterrows()))
    M["_r5_out_two"] = [r["material"] for _, r in rows.iterrows() if not t2_consistent(r, "two")]
    M["_r5_out_strict"] = [r["material"] for _, r in rows.iterrows() if not t2_consistent(r, "strict")]
    zif = c.t2.loc[c.t2["material"] == "ZIF-8/Al2O3"].iloc[0]
    M["r5_zif8_ratio"] = (float(zif["thickness_um"]) * float(zif["permeance_1e-8"]) / 10.0
                          / float(zif["permeability_1e-13"]))
    pd_rows = c.t2.iloc[:7]                       # the Pd-based films at the top of the table
    q = [float(s.replace(",", "")) for s in pd_rows["permeance_1e-8"]]
    M["_r5_n_pd"] = len(q)
    M["r5_pd_span"] = max(q) / min(q)
    # SINTEF triple: text Q_s (p. 46) against Table 2's permeance at 1 bar, two readings
    sint = c.t2.loc[c.t2["institute"].str.startswith("SINTEF")].iloc[0]
    m_s, e_s = printed(c.pn, "p46_Qs_pure_mantissa"), num(c.pn, "p46_Qs_pure_exponent")
    qs = float(m_s) * 10.0 ** e_s
    qs_lo = (float(m_s) - half_unit(m_s, "two")) * 10.0 ** e_s
    p1 = num(c.pn, "t2_footnote_pH2_bar") * 1e5
    pe_hi = interval(sint["permeance_1e-8"], "two")[1] * 1e-8
    pe = float(sint["permeance_1e-8"].replace(",", "")) * 1e-8
    los = []
    for tag, fac in (("secant", 1.0), ("tangent", 2.0)):
        M[f"r5_sintef_ratio_{tag}"] = qs / (fac * pe * math.sqrt(p1))
        los.append(qs_lo / (fac * pe_hi * math.sqrt(p1)))
        M[f"r5_sintef_ratio_lo_{tag}"] = los[-1]
    M["_r5_sintef_los"] = np.array(los)
    M["r5_n_readings_inconsistent"] = int(np.sum(np.array(los) > 1.0))


# ---------------------------------------------------------------- R4: module
def part_local(c, M):
    theta, phi = tc(c)
    b = c.beta_tc
    M["r4_u_in"] = wall_root(1.0, theta, b)
    M["r4_j_in"] = local_flux(1.0, theta, b)
    M["r4_share_in"] = film_share(1.0, theta, b)
    M["r4_epsQ_in"], M["_r4_epsk_in"] = elasticities(1.0, theta, b)
    M["r4_R_max"] = pinch(theta, phi)[1]
    M["r4_N80_quad"] = module_number(c.r_star, theta, phi, b)
    M["r4_rec_N2_quad"] = recovery_at(c.n_tc, theta, phi, b)
    M["r4_N80_mem"] = module_number(c.r_star, theta, phi, 0.0)
    M["r4_N80_over_mem"] = M["r4_N80_quad"] / M["r4_N80_mem"]
    M["tc_r_star"] = c.r_star
    M["r4_N80_mem_vac"] = n_membrane_vacuum(c.r_star, phi)
    M["r4_a80_film"] = a_film(c.r_star, theta, phi)
    M["r4_beta_cross"] = M["r4_N80_mem"] / M["r4_a80_film"]
    M["r4_additivity"] = M["r4_N80_quad"] / (M["r4_N80_mem"] + b * M["r4_a80_film"]) - 1.0
    w = c.ea_weight
    M["r4_eA_quad"] = area_elasticity(c.r_star, theta, phi, b, weight=w)
    M["r4_beta_half_quad"] = one_root(lambda x: area_elasticity(c.r_star, theta, phi, x, weight=w) - 0.5,
                                      1e-2, 1e2)
    for rs in (0.5, 0.8, 0.9):
        for tname, th in (("vac", 0.0), ("tc", theta)):
            M[f"r4_beta_half_R{int(100 * rs)}_{tname}"] = one_root(
                lambda x: area_elasticity(rs, th, phi, x, weight=w) - 0.5, 1e-2, 1e2)


def part_sweep(c, M):
    """R4 by pymrm: N(beta) for R* on three levels, extrapolated, against the quadrature."""
    theta, phi = tc(c)
    Ns, ests, gaps = [], [], []
    for b in BETAS:
        x, e, _ = extrap(lambda n: area_solve(c, b, n)[0], c.nz0)
        q = module_number(c.r_star, theta, phi, b, n=c.n_route2)
        Ns.append(x)
        ests.append(abs(e) / x)
        gaps.append(abs(x - q) / q)
    M["_r4_sweep_N"] = np.array(Ns)
    M["r4_sweep_max_est"] = max(ests)
    M["r4_sweep_max_gap"] = max(gaps)
    M["r4_sweep_gap_over_est"] = max(g / e for g, e in zip(gaps, ests))


def part_tcmod(c, M):
    """pymrm values at TC, each with its grid error (rule 36)."""
    b = c.beta_tc
    x, e, o = extrap(lambda n: area_solve(c, b, n)[0], c.nz0)
    M["r4_N80_pymrm"], M["r4_N80_pymrm_est"], M["_r4_N80_order"] = x, abs(e), o
    x, e, _ = extrap(lambda n: area_solve(c, b, n)[1], c.nz0)
    M["r4_eA_pymrm"], M["r4_eA_pymrm_est"] = x, abs(e)
    M["v5_eA_gap"] = abs(M["r4_eA_pymrm"] - area_elasticity(c.r_star, *tc(c), b, weight=c.ea_weight))

    def bh(n):
        return one_root(lambda z: area_solve(c, z, n)[1] - 0.5, 1e-2, 1e2, n_scan=12)
    x, e, _ = extrap(bh, c.nz0)
    M["r4_beta_half_pymrm"], M["r4_beta_half_pymrm_est"] = x, abs(e)


def part_v(c, M):
    theta, phi = tc(c)
    # V1, V3 on the finest module level at TC
    rec, bal, m = recovery_solve(c, c.n_tc, 4 * c.nz0)
    forms = m.cell_flux_forms()
    M["v1_flux_forms_gap"] = float(np.max(np.abs(forms - forms[0])) / np.max(np.abs(forms[0])))
    M["v3_balance"] = bal
    # V2 limits by the general route
    M["v2_mem_vac_gap"] = abs(module_number(c.r_star, 0.0, phi, c.beta_small)
                              / n_membrane_vacuum(c.r_star, phi) - 1.0)
    M["v2_film_gap"] = abs(module_number(c.r_star, theta, phi, c.beta_large) / c.beta_large
                           / a_film(c.r_star, theta, phi) - 1.0)
    r50 = recovery_solve(c, c.n_large, 4 * c.nz0)[0]
    M["v2_pinch_gap"] = abs(r50 - pinch(theta, phi)[1])
    # V4 ladder at N = 2, beta = 1
    ns = [100 * 2 ** k for k in range(6)]
    vals = [recovery_solve(c, c.n_tc, n)[2].inlet_permeate() if c.v4_read == "inlet"
            else recovery_solve(c, c.n_tc, n)[0] for n in ns]
    o = observed_orders(lambda n: vals[ns.index(n)], ns)
    M["_v4_orders"] = np.array(o["orders"])
    M["v4_order_min"], M["v4_order_max"] = float(np.min(o["orders"])), float(np.max(o["orders"]))
    M["r4_rec_N2_pymrm"], M["r4_rec_N2_pymrm_est"] = o["extrapolated"], abs(o["error_estimate"])
    M["v4_rich_gap_over_est"] = abs(o["extrapolated"] - CV_TABLE["TC_recovery_N2_beta1"]) / abs(o["error_estimate"])
    M["v4_finest_diff"] = abs(vals[-1] - vals[-2])


CV_TABLE = dict(zip(load_data("check_values.csv", page=PAGE)["name"],
                    load_data("check_values.csv", page=PAGE)["value"]))
PARTS = {"f15": part_f15, "f16": part_f16, "p57": part_p57, "t2": part_t2, "local": part_local,
         "sweep": part_sweep, "tcmod": part_tcmod, "v": part_v}


def study(c, parts=tuple(PARTS)):
    """Every number the page reports, for one configuration."""
    M = {}
    for p in parts:
        PARTS[p](c, M)
    return M


M = study(C0)
print(f"study: {len(CACHE)} module solves cached")

# %% [markdown]
# Before any claim, the page's functions are asserted against check values computed by a
# second route (closed forms, quadrature, root finding; no pymrm and no page code;
# `data/check_values.meta.yaml`), at TC and at the test inputs the sidecar lists.

# %%
TH, PHI = tc(C0)
CV = {"TC_u_wall_inlet": M["r4_u_in"], "TC_j_inlet": M["r4_j_in"], "TC_film_share_inlet": M["r4_share_in"],
      "TC_epsQ_inlet": M["r4_epsQ_in"], "TC_epsk_inlet": M["_r4_epsk_in"], "TC_R_max_pinch": M["r4_R_max"],
      "TC_N_R80_beta1": M["r4_N80_quad"], "TC_recovery_N2_beta1": M["r4_rec_N2_quad"],
      "TC_f_out_N2_beta1": 1.0 - M["r4_rec_N2_quad"], "TC_mean_j_N2_beta1": M["r4_rec_N2_quad"] / N_TC,
      "TC_N_R80_membrane_limit": M["r4_N80_mem"], "TC_a_R80_film_limit": M["r4_a80_film"],
      "TC_N_R80_membrane_limit_vacuum": M["r4_N80_mem_vac"], "TC_beta_cross_R80": M["r4_beta_cross"],
      "TC_area_additivity_defect_R80": M["r4_additivity"], "TC_area_elasticity_R80": M["r4_eA_quad"],
      "TC_beta_half_area_R80": M["r4_beta_half_quad"],
      "P57_epsk_eta0p90_vacuum": eps_k_from_eta(0.90, 0.0), "P57_epsk_eta0p90_y0p5": eps_k_from_eta(0.90, 0.5),
      "P57_epsk_eta0p98_vacuum": eps_k_from_eta(0.98, 0.0),
      "F16_epsQ_from_shares_test": share_elasticities(0.80, 0.12, 0.01, 0.5)[0],
      "F16_epsk_from_shares_test": share_elasticities(0.80, 0.12, 0.01, 0.5)[1],
      "F16_epsQ_upper_bound_test": eps_q_bound(0.80, 0.12),
      "F15_wall_factor": wall_factor(R0, RL), "F15_D_centre_test": centre_depletion(0.30, R0, RL),
      "F15_d_1D_test": cupmix_depletion(centre_depletion(0.30, R0, RL)),
      "F15_epsQ_vacuum_test": eps_q_vacuum(1.0 - cupmix_depletion(centre_depletion(0.30, R0, RL))),
      "F15_theta_star_test": theta_star(1.0 - cupmix_depletion(centre_depletion(0.30, R0, RL))),
      "F15_lambda_star_test": lambda_star(cupmix_depletion(centre_depletion(0.30, R0, RL)))}
with contextlib.redirect_stdout(io.StringIO()):
    DEV = assert_check_values(CV, page=PAGE)
say("Check values: <<n:d>> of <<t:d>> rows agree, each inside its own tolerance; the largest relative "
    "deviation is <<w:.1e>>.", n=len(DEV), t=len(CV), w=max(abs(v) for v in DEV.values()))

# %% [markdown]
# ## Results
#
# ### R1. Fig. 15 is a parabola
#
# Each curve is fitted with the one-parameter profile of E15 (normalised at the first
# marker, which every curve passes at 1). The residuals are in units of the read uncertainty.

# %%
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
rr = np.linspace(R0, 1.0, 200)
for cv, mk in zip(CURVES, ("<", "^", "s")):
    d_ = RD15[RD15["curve"] == cv]
    D_ = M[f"r1_D_fit_{TAG[cv]}"]
    ax[0].plot(d_["r_over_R"], d_["rel_H2_wt_fraction"], mk, mfc="none", label=f"{cv} (read)")
    ax[0].plot(rr, (1 - D_ * rr ** 2) / (1 - D_ * R0 ** 2), "-", lw=1)
    sel = d_["cell"] > 1
    ax[1].plot(d_["r_over_R"][sel], ((1 - D_ * d_["r_over_R"][sel] ** 2) / (1 - D_ * R0 ** 2)
                                    - d_["rel_H2_wt_fraction"][sel]) / U15, mk, mfc="none")
ax[0].set(xlabel="r/R", ylabel="relative H$_2$ wt fraction", title="Fig. 15 reads and E15 fits")
ax[0].legend(fontsize=8)
ax[1].axhspan(-1, 1, color="0.9")
ax[1].set(xlabel="r/R", ylabel="fit minus read, units of u", title="residuals")
plt.tight_layout()
plt.show()
n_in = int(np.sum(np.abs(M["_r1_resid_u"]) < 1.0))
say("Of <<n>> points read (squares 2 to 12, triangles 5 to 12; the triangles' markers 1 to 4 are hidden "
    "behind the squares), "
    "<<k>> lie within one read uncertainty of their parabola; RMS residual <<rms:.2f>> u, largest "
    "<<mx:.2f>> u. Fitted centre-to-wall depletions $D$: <<a:.4f>> (1x), <<b:.4f>> (2x), <<c:.4f>> (4x); "
    "from the last marker alone (E15, factor <<wf:.6f>>): <<a2:.4f>>, <<b2:.4f>>, <<c2:.4f>>, within "
    "<<g:.2f>> read units of the fit. The parabola licenses the lumped film of E3: the cup-mix (1-D) "
    "wall depletions are <<e1:.4f>>, <<e2:.4f>>, <<e3:.4f>>.",
    n=M["_r1_n_read"], k=n_in, rms=M["r1_rms_u"], mx=M["r1_max_u"], a=M["r1_D_fit_1x"],
    b=M["r1_D_fit_2x"], c=M["r1_D_fit_4x"], wf=M["r1_wall_factor"], a2=M["r1_D_last_1x"],
    b2=M["r1_D_last_2x"], c2=M["r1_D_last_4x"],
    g=max(abs(M[f"r1_D_gap_u_{t}"]) for t in ("1x", "2x", "4x")),
    e1=M["r1_d1D_1x"], e2=M["r1_d1D_2x"], e3=M["r1_d1D_4x"])

# %% [markdown]
# ### R2. S1 on Fig. 15: what decides "quite limited"
#
# Under vacuum permeate the membrane's share of control $\varepsilon_Q$ follows from the wall
# depletion alone (E7); a permeate pressure only lowers it. So the vacuum value is an upper
# bound of membrane control, and $\theta^*_b$ is the permeate-to-bed H2 ratio above which the
# film limits ($\varepsilon_Q < 1/2$). Two references for the bulk: the 2-D centre
# ($x^2 = 1 - D$) and the 1-D cup-mix mean ($x^2 = 1 - d_{1D}$).

# %%
rows_ = ["| curve | ref. | vacuum $\\varepsilon_Q$ | $\\theta^*_b$ | $\\pm$ one read uncertainty | $\\theta^*_b$ at $n$ = 1 |",
         "|---|---|---|---|---|---|"]
for cv in CURVES:
    for ref in ("2D", "1D"):
        rows_.append(f"| {cv} | {ref} | {M[f'r2_epsQ_vac_{ref}_{TAG[cv]}']:.4f} | "
                     f"{M[f'r2_theta_star_{ref}_{TAG[cv]}']:.4f} | "
                     f"{M.get(f'r2_theta_star_u_{ref}_{TAG[cv]}', float('nan')):.4f} | "
                     f"{M.get(f'r2_theta_star_n1_{ref}_{TAG[cv]}', float('nan')):.4f} |")
display(Markdown("\n".join(rows_)))
say("Membrane control under vacuum is at least <<lo:.3f>> on every curve and reference, so the film "
    "of Fig. 15 limits only through the permeate. The film is not the limiting step at 1x unless the "
    "permeate H2 pressure exceeds <<t1a:.3f>> $\\pm$ <<u1a:.3f>> (centre; <<n1a:.3f>> at $n$ = 1) or "
    "<<t1b:.3f>> $\\pm$ <<u1b:.3f>> (cup-mix; <<n1b:.3f>> at $n$ = 1) of the bed's; \"could become a "
    "limiting step\" at 4x needs a ratio above <<t4a:.3f>> $\\pm$ <<u4a:.3f>> (centre; <<n4a:.3f>> at "
    "$n$ = 1) or <<t4b:.3f>> $\\pm$ <<u4b:.3f>> (cup-mix; <<n4b:.3f>> at $n$ = 1); each $\\pm$ is "
    "one read uncertainty on the last marker, propagated. The review prints neither ratio, so Fig. 15 "
    "alone cannot decide S1; these ratios are the thresholds for \"a limiting step\". Under vacuum "
    "the film holds <<ek1:.3f>> of the control at 1x (centre): \"quite limited\" is a statement of "
    "that size, which the 1/2 threshold does not test. A film that is not limiting is not necessarily "
    "a small one: at a permeate ratio of <<tx:.1f>> the film holds <<exk:.3f>> of the control at 1x "
    "and the flux is <<exf:.3f>> of the flux without the film. A smaller partial-pressure depletion "
    "than the weight-fraction read (A1) raises the thresholds.",
    lo=min(M[f"r2_epsQ_vac_{r}_{t}"] for r in ("2D", "1D") for t in ("1x", "2x", "4x")),
    t1a=M["r2_theta_star_2D_1x"], t1b=M["r2_theta_star_1D_1x"], t4a=M["r2_theta_star_2D_4x"],
    t4b=M["r2_theta_star_1D_4x"], u1a=M["r2_theta_star_u_2D_1x"], u1b=M["r2_theta_star_u_1D_1x"],
    u4a=M["r2_theta_star_u_2D_4x"], u4b=M["r2_theta_star_u_1D_4x"], n1a=M["r2_theta_star_n1_2D_1x"],
    n1b=M["r2_theta_star_n1_1D_1x"], n4a=M["r2_theta_star_n1_2D_4x"], n4b=M["r2_theta_star_n1_1D_4x"],
    ek1=1.0 - M["r2_epsQ_vac_2D_1x"], tx=THETA_EX, exk=M["r2_ex_epsk_2D_1x"], exf=M["r2_ex_flux_2D_1x"])
say("Proportionality (E8, vacuum, fixed bulk): from the 1x curve the film would take two thirds of the "
    "drop at a permeability <<la:.2f>> (centre) or <<lb:.2f>> (cup-mix) times Tiemersma's. The predicted "
    "4x depletion is <<p4:.4f>> (centre) against <<o4:.4f>> read, read minus prediction <<g4:+.1f>> "
    "combined read units; 2x: <<g2:+.1f>> units (cup-mix: <<h2:+.1f>> at 2x, <<h4:+.1f>> at 4x). The 2x "
    "curve follows the scaled film; at 4x the read is deeper by the gap stated above. A "
    "permeate pressure does not cause that: at fixed bulk it makes the predicted depletion shallower, "
    "<<p4t:.4f>> at the 1x deciding ratio (<<g4t:+.1f>> units from the read; 2x: <<g2t:+.1f>>). At a "
    "fixed centre value the wall flux scales as $D$, so the three curves give the control directly: "
    "$\\ln(D_4/D_1)/\\ln 4$ = <<fd:.3f>> $\\pm$ <<sfd:.3f>>. A control below 1/2 would need the "
    "centre H2 fraction at 4x to be <<cr:.3f>> of its 1x value or less; the normalised figure does not "
    "print the centre values. For orientation, the seven Pd-based films at the top of Table 2 span a "
    "factor <<sp:.1f>> in permeance.",
    la=M["r2_lambda_star_2D"], lb=M["r2_lambda_star_1D"], p4=M["r2_pred_2D_4x"], o4=M["r1_D_last_4x"],
    g4=M["r2_gap_units_2D_4x"], g2=M["r2_gap_units_2D_2x"], h2=M["r2_gap_units_1D_2x"],
    p4t=M["r2_pred_2D_4x_thstar"], g4t=M["r2_gap_units_2D_4x_thstar"], g2t=M["r2_gap_units_2D_2x_thstar"],
    fd=M["r2_fd_control_1x_4x"], sfd=M["r2_fd_control_u"], cr=M["r2_centre_ratio_half"],
    h4=M["r2_gap_units_1D_4x"], sp=M["r5_pd_span"])

# %% [markdown]
# ### R3. S2 on Fig. 16: "the limiting factor"
#
# From the shares of the drop alone, E9 bounds membrane control for any Pd/support interface
# pressure and any support share; only the ratio $c/m$ of the film (concentration
# polarisation) and membrane shares enters. The depletion share drops out.

# %%
fig, ax = plt.subplots(figsize=(7, 3.8))
xs = np.arange(len(PRESSURES))
segs = [("depletion", RD16["depletion_top"], np.zeros(4), "#9fa8da"),
        ("concentration polarisation", RD16["cp_top"], RD16["depletion_top"], "#8e3a5f"),
        ("Pd membrane", RD16["pd_top"], RD16["cp_top"], "#fff3b0"),
        ("support", RD16["support_top"], RD16["pd_top"], "#b2ebf2")]
for lab, top, bot, col in segs:
    ax.bar(xs, top - bot, bottom=bot, color=col, edgecolor="k", width=0.4, label=lab)
ax.set(xticks=xs, xticklabels=[str(p) for p in PRESSURES], xlabel="feed pressure (bar)",
       ylabel="share of the H$_2$ pressure drop (read)")
ax2 = ax.twinx()
ax2.plot(xs, [M[f"r3_bound_{P}"] for P in PRESSURES], "ko-", label="bound on $\\varepsilon_Q$ ($n$ = 1/2)")
ax2.plot(xs, [M[f"r3_bound_n1_{P}"] for P in PRESSURES], "k^--", mfc="none", label="bound, $n$ = 1")
ax2.axhline(0.5, color="r", lw=1)
ax2.set(ylim=(0, 1), ylabel="membrane control $\\varepsilon_Q$")
ax.legend(fontsize=7, loc="lower left")
ax2.legend(fontsize=7, loc="lower right")
plt.title("Fig. 16 redrawn from the reads, with the bound of E9")
plt.tight_layout()
plt.show()
rows_ = ["| feed (bar) | $c$ | $m$ | $c/m$ | bound $n$ = 1/2 | bound $n$ = 1 | margin from 1/2 (read units) |",
         "|---|---|---|---|---|---|---|"]
for P in PRESSURES:
    rows_.append(f"| {P} | {M[f'r3_c_{P}']:.4f} | {M[f'r3_m_{P}']:.4f} | {M[f'r3_cm_{P}']:.3f} | "
                 f"{M[f'r3_bound_{P}']:.4f} | {M[f'r3_bound_n1_{P}']:.4f} | {M[f'r3_margin_u_{P}']:.1f} |")
display(Markdown("\n".join(rows_)))
say("Membrane control is at most <<bmax:.3f>> at every feed pressure. Pressures with the bound below 1/2: <<n:d>> of <<t:d>>, the closest "
    "<<mg:.1f>> propagated read uncertainties away. Membrane control below 1/2 does not by itself make "
    "the film limiting. The film itself keeps at least <<ekm:.3f>> of the control when the permeate H2 "
    "pressure is not negative ($s \\ge \\sigma/(\\sigma + m)$); without that limit the support could "
    "take the control. On the elasticity scale S2's \"the limiting factor\" holds for the device of "
    "Fig. 16 with no unprinted number. It is a permeation cell with "
    "sweep gas (Caravella et al., 2009, p. 623), so the statement is about that cell, not a packed bed.",
    bmax=max(M[f"r3_bound_{P}"] for P in PRESSURES),
    n=M["r3_n_below_half"], t=len(PRESSURES), mg=float(np.min(M["_r3_margins"])),
    ekm=min(M[f"r3_epsk_min_{P}"] for P in PRESSURES))

# %% [markdown]
# ### R4. S3 at module scale (pymrm)
#
# At TC, the membrane area for 80 % recovery in film units $a = N/\beta$, from the membrane
# limit ($a \to N_{mem}/\beta$) to the film limit ($a \to a_{film}$). The area elasticity
# $e_A$ is the area-weighted mean of the local $\varepsilon_Q$ (E14): above
# $\beta_{1/2}$ ($e_A = 1/2$) a better membrane buys less area than a better film.

# %%
th_, ph_ = tc(C0)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
aq = np.array([module_number(R_STAR, th_, ph_, b) / b for b in BETAS])
ax[0].loglog(BETAS, aq, "-", color="0.4", label="quadrature (E12)")
ax[0].loglog(BETAS, M["_r4_sweep_N"] / BETAS, "o", ms=3, mfc="none", label="pymrm, extrapolated")
ax[0].loglog(BETAS, M["r4_N80_mem"] / BETAS, "--", lw=1, label="membrane limit $N_{mem}/\\beta$")
ax[0].axhline(M["r4_a80_film"], ls=":", color="k", lw=1, label="film limit $a_{film}$")
ax[0].set(xlabel="$\\beta$", ylabel="area in film units $a$", title="area for 80 % recovery at TC")
ax[0].legend(fontsize=7)
eA = [area_elasticity(R_STAR, th_, ph_, b) for b in BETAS]
ax[1].semilogx(BETAS, eA, "-")
ax[1].axhline(0.5, color="r", lw=1)
ax[1].axvline(M["r4_beta_half_quad"], color="0.5", ls="--", lw=1)
ax[1].set(xlabel="$\\beta$", ylabel="$e_A$", title="area elasticity (E14)")
plt.tight_layout()
plt.show()
say("At TC and $\\beta$ = <<b:g>>: $N$ = <<n:.6f>> (pymrm, extrapolated; grid error of the "
    "<<nz:d>>-cell value <<ne:.1e>>; quadrature <<nq:.6f>>), $e_A$ = <<e:.5f>> (pymrm, grid error <<ee:.1e>>). "
    "Limits: $N_{mem}$ = <<nm:.5f>> (vacuum <<nv:.5f>>), $a_{film}$ = <<af:.5f>>; the asymptotes cross at "
    "$\\beta$ = <<bc:.4f>>, and at $\\beta$ = 1 the area differs from their sum by <<ad:+.4f>> (relative). $e_A$ = 1/2 at "
    "$\\beta_{1/2}$ = <<bh:.5f>> (pymrm, grid error <<bhe:.1e>>; quadrature <<bhq:.6f>>). The inlet "
    "already shows the point of the elasticity scale: the film takes <<sh:.3f>> of the drop there "
    "but the membrane keeps <<eq:.3f>> of the control.",
    b=BETA_TC, n=M["r4_N80_pymrm"], nz=4 * NZ0, ne=M["r4_N80_pymrm_est"], nq=M["r4_N80_quad"],
    e=M["r4_eA_pymrm"], ee=M["r4_eA_pymrm_est"], nm=M["r4_N80_mem"], nv=M["r4_N80_mem_vac"],
    af=M["r4_a80_film"], bc=M["r4_beta_cross"], ad=M["r4_additivity"], bh=M["r4_beta_half_pymrm"],
    bhe=M["r4_beta_half_pymrm_est"], bhq=M["r4_beta_half_quad"], sh=M["r4_share_in"], eq=M["r4_epsQ_in"])
rows_ = ["| $R^*$ | $\\beta_{1/2}$, vacuum permeate | $\\beta_{1/2}$, TC permeate |", "|---|---|---|"]
for rs in (50, 80, 90):
    rows_.append(f"| {rs / 100:.1f} | {M[f'r4_beta_half_R{rs}_vac']:.4f} | {M[f'r4_beta_half_R{rs}_tc']:.4f} |")
display(Markdown("\n".join(rows_)))
say("The second half of S3 holds as a statement about the area: above $\\beta_{1/2}$ the film sets most "
    "of it. The first half (\"can be neglected\") does not follow from membrane control: at $\\beta$ = "
    "<<b:g>> the membrane holds $e_A$ = <<e:.3f>> while the area is <<rn:.3f>> times the "
    "membrane-limit area. The threshold sits near $\\beta$ = 1 for every recovery and permeate in the table "
    "(<<lo:.3f>> to <<hi:.3f>>): the film determines the area once its conductance is below the "
    "membrane's at the feed, $k_p p_{in} < Q_s\\sqrt{p_{in}}$, roughly. Pinch: recovery cannot exceed "
    "<<rm:.5f>> with a 1 bar permeate.",
    b=BETA_TC, e=M["r4_eA_pymrm"], rn=M["r4_N80_over_mem"],
    lo=min(M[f"r4_beta_half_R{r}_{t}"] for r in (50, 80, 90) for t in ("vac", "tc")),
    hi=max(M[f"r4_beta_half_R{r}_{t}"] for r in (50, 80, 90) for t in ("vac", "tc")), rm=M["r4_R_max"])

# %% [markdown]
# ### R5. Printed numbers: Table 2 and the SINTEF triple
#
# Table 2 prints thickness, permeance and permeability; permeance times thickness is the
# permeability (units $10^{-8}$ times $\mu$m equals $10^{-13}$ times 10). Each row is
# compared with the half units of every entry it combines, from the strings.

# %%
say("Rows consistent within their printed precision: <<n:d>> of <<t:d>> when trailing zeros of an integer "
    "are placeholders but every entry carries at least two significant figures (1,500 stands for 1,450 "
    "to 1,550, 20 for 19.5 to 20.5); the row outside is <<o2>>. With every digit significant: <<s:d>> of "
    "<<t:d>> (outside: <<os>>); with every trailing zero a placeholder (20 standing for 15 to 25): <<a:d>> "
    "of <<t:d>>. For ZIF-8/Al2O3 (McCarthy et al., 2010) thickness times permeance is <<z:.3f>> of the "
    "printed permeability; only the one-figure reading of its thickness admits it. Which entry is off is "
    "not decidable from the table.",
    n=M["r5_n_consistent_two"], s=M["r5_n_consistent_strict"], a=M["r5_n_consistent_one"],
    t=M["_r5_n_rows"], z=M["r5_zif8_ratio"], o2=", ".join(M["_r5_out_two"]), os=", ".join(M["_r5_out_strict"]))
say("SINTEF (Peters et al., 2011a): the text gives $Q_s$ = <<m>> $\\times 10^{<<e:.0f>>}$ mol m$^{-2}$ "
    "s$^{-1}$ Pa$^{-0.5}$ (p. 46); Table 2 gives a permeance of <<pt>> $\\times 10^{-8}$ mol m$^{-2}$ "
    "s$^{-1}$ Pa$^{-1}$ at 1 bar. At $n$ = 1/2 the text value is <<rs:.3f>> times the table's (secant "
    "reading, $Q_s$ = permeance $\\times\\sqrt{1\\,\\mathrm{bar}}$) or <<rt:.3f>> times (tangent reading); "
    "the lowest ratios the printed digits allow are <<ls:.3f>> and <<lt:.3f>>, so <<k:d>> of 2 readings "
    "are inconsistent. The page does not say which number is wrong; a flux exponent other than 1/2 "
    "between 1 and 25 bar would also do it (H1.1).",
    m=printed(PN, "p46_Qs_pure_mantissa"), e=num(PN, "p46_Qs_pure_exponent"),
    pt=T2.loc[T2["institute"].str.startswith("SINTEF"), "permeance_1e-8"].iloc[0],
    rs=M["r5_sintef_ratio_secant"], rt=M["r5_sintef_ratio_tangent"], ls=M["r5_sintef_ratio_lo_secant"],
    lt=M["r5_sintef_ratio_lo_tangent"], k=M["r5_n_readings_inconsistent"])

# %% [markdown]
# ### R6. One scale
#
# The membrane's share of control from each piece of evidence. These are different membranes
# and devices: the table places them on one scale, it does not compare packed and fluidized
# beds. The p. 57 effectiveness factor is read "as if the whole deficit were a series film"
# (the review names fine dust or other mass transfer resistances).

# %%
rows_ = ["| evidence | membrane control $\\varepsilon_Q$ | note |", "|---|---|---|"]
for cv in CURVES:
    rows_.append(f"| Fig. 15, {cv} | at most {M[f'r2_epsQ_vac_2D_{TAG[cv]}']:.3f} (centre), "
                 f"{M[f'r2_epsQ_vac_1D_{TAG[cv]}']:.3f} (cup-mix) | film limits if $p_p/p_b$ > "
                 f"{M[f'r2_theta_star_2D_{TAG[cv]}']:.3f} (centre) |")
for P in PRESSURES:
    rows_.append(f"| Fig. 16, {P} bar | at most {M[f'r3_bound_{P}']:.3f} | any interface, any support |")
for t, item in (("lo", "p57_eta_low"), ("hi", "p57_eta_high")):
    rows_.append(f"| p. 57, $\\eta$ = {printed(PN, item)} | {1 - M[f'r6_epsk_vac_{t}']:.3f} (vacuum) to "
                 f"{1 - M[f'r6_epsk_y1_{t}']:.3f} (permeate near bed) | as a series film |")
display(Markdown("\n".join(rows_)))

# %% [markdown]
# ### The finding

# %%
say("**On the elasticity scale the review's own numbers support S2 and leave S1 open.** Fig. 16: membrane "
    "control is bounded by <<b10:.3f>>, <<b15:.3f>>, <<b20:.3f>> and <<b26:.3f>> at 10, 15, 20 and 26 bar, below "
    "1/2 at <<n:d>> of <<t:d>> pressures with at least <<mg:.1f>> read uncertainties to spare, and the film "
    "itself keeps at least <<ekm:.3f>> of the control: the film is \"the limiting factor\" there. Fig. 15: "
    "the film is the limiting step only if the permeate-to-bed H2 ratio exceeds <<t1:.3f>> $\\pm$ <<u1:.3f>> "
    "(1x) or <<t4:.3f>> $\\pm$ <<u4:.3f>> (4x) (centre reference, one read uncertainty; <<t1b:.3f>> and "
    "<<t4b:.3f>> cup-mix); the review prints no permeate pressure. At a fixed centre value the three curves "
    "give a control of <<fd:.3f>> between 1x and 4x, so S1 stays open only through the centre values the "
    "normalised figure hides. Module (S3): the area stops responding mainly to the permeance once "
    "$\\beta$ exceeds <<bh:.4f>> (grid error <<bhe:.1e>>) at TC. Table 2: <<c2:d>> of <<nr:d>> rows "
    "consistent at two significant figures (<<cs:d>> with every digit significant, <<c1:d>> at one figure).",
    b10=M["r3_bound_10"], b15=M["r3_bound_15"], b20=M["r3_bound_20"], b26=M["r3_bound_26"],
    n=M["r3_n_below_half"], t=len(PRESSURES), mg=float(np.min(M["_r3_margins"])),
    t1=M["r2_theta_star_2D_1x"], t4=M["r2_theta_star_2D_4x"], t1b=M["r2_theta_star_1D_1x"],
    t4b=M["r2_theta_star_1D_4x"], bh=M["r4_beta_half_pymrm"], bhe=M["r4_beta_half_pymrm_est"],
    u1=M["r2_theta_star_u_2D_1x"], u4=M["r2_theta_star_u_2D_4x"], fd=M["r2_fd_control_1x_4x"],
    ekm=min(M[f"r3_epsk_min_{P}"] for P in PRESSURES),
    c2=M["r5_n_consistent_two"], nr=M["_r5_n_rows"], cs=M["r5_n_consistent_strict"], c1=M["r5_n_consistent_one"])

# %% [markdown]
# ## Validation
#
# - **V1, local flux**: the three forms of E5 at every cell of the finest TC solve.
# - **V2, limits**: the quadrature route at $\beta$ = 1e-6 with vacuum permeate against the
#   closed-form membrane limit; at $\beta$ = 1e6 against the closed-form film limit; the module
#   at $N$ = 50 against the pinch recovery.
# - **V3, conservation**: $f + g = 1$ in every cell (an identity of the scheme).
# - **V4, order**: recovery at $N$ = 2, $\beta$ = 1 on 100 to 3200 cells; observed orders of
#   successive triples, Richardson value against the check value within its estimate.
# - **V5, second route**: at every $\beta$ of R4 the extrapolated pymrm area against the
#   quadrature; $e_A$ by the discrete sensitivity against E14.
# - **V6, check values**: above, before the first result.

# %%
say("V1: largest relative gap between the flux forms <<v1:.1e>> (round-off). V2: <<a:.1e>> (membrane "
    "limit, $O(\\beta)$), <<b:.1e>> (film limit, $O(1/\\beta)$), <<c:.1e>> (pinch). V3: largest "
    "$|f + g - 1|$ <<v3:.1e>>. V4: orders <<o1:.4f>> to <<o2:.4f>> (first-order upwind), extrapolated "
    "recovery <<r:.9f>> with estimate <<re:.1e>>, <<g:.1e>> estimates from the check value; finest "
    "difference <<fd:.1e>>. V5: extrapolated areas within <<sg:.1e>> of the quadrature at every $\\beta$ "
    "(<<so:.3f>> of the grid estimate at most; largest estimate <<se:.1e>>), $e_A$ within <<eg:.1e>>.",
    v1=M["v1_flux_forms_gap"], a=M["v2_mem_vac_gap"], b=M["v2_film_gap"], c=M["v2_pinch_gap"],
    v3=M["v3_balance"], o1=M["v4_order_min"], o2=M["v4_order_max"], r=M["r4_rec_N2_pymrm"],
    re=M["r4_rec_N2_pymrm_est"], g=M["v4_rich_gap_over_est"], fd=M["v4_finest_diff"],
    sg=M["r4_sweep_max_gap"], so=M["r4_sweep_gap_over_est"], se=M["r4_sweep_max_est"], eg=M["v5_eA_gap"])

# %% [markdown]
# ### The claims

# %%
RTOL = 1e-6
STRUCTURAL = {"v1_flux_forms_gap": "algebraic identity of the three forms of E5 (round-off)",
              "v3_balance": "H2 identity of the conservative scheme (round-off)"}
BOUNDS = {"v1_flux_forms_gap": (None, 1e-13, "identity", STRUCTURAL["v1_flux_forms_gap"]),
          "v3_balance": (None, 1e-12, "identity", STRUCTURAL["v3_balance"]),
          "v2_mem_vac_gap": (None, 1e-5, "theory", "O(beta) at beta = 1e-6"),
          "v2_film_gap": (None, 1e-5, "theory", "O(1/beta) at beta = 1e6"),
          "v2_pinch_gap": (None, 1e-6, "theory", "approach to the pinch is geometric in N"),
          "v4_order_min": (0.9, 1.1, "theory", "first-order upwind (backward Euler in xi)"),
          "v4_order_max": (0.9, 1.1, "theory", "first-order upwind (backward Euler in xi)"),
          "v4_rich_gap_over_est": (None, 1.0, "chosen", "Richardson value within its own estimate"),
          "v4_finest_diff": (1e3 * 3200 * 2.2e-16, None, "theory", "finest difference at least 1e3 x round-off (n_z eps)"),
          "r4_sweep_gap_over_est": (None, 0.1, "chosen", "extrapolated gap at most a tenth of the grid estimate"),
          "r4_sweep_max_gap": (None, 1e-5, "chosen", "second-order remainder after Richardson"),
          "v5_eA_gap": (None, 1e-5, "chosen", "second-order remainder after Richardson")}
TIGHT = {k: 1e-8 for k in ("r4_N80_pymrm", "r4_eA_pymrm", "r4_beta_half_pymrm", "r4_rec_N2_pymrm")}
COUNTS = {"r3_n_below_half": ("_r3_margins", 0.0, 3.0, "margins in read units, at least 3 (chosen)"),
          "r5_n_readings_inconsistent": ("_r5_sintef_los", 1.0, 0.05, "lowest ratio at least 5 % above 1 (chosen)")}
EXACT = tuple(f"r5_n_consistent_{r_}" for r_ in READINGS) + ("tc_r_star",)


def register(M):
    """Register every reported number with its kind."""
    for name, v in M.items():
        if name.startswith("_"):
            continue
        if name in BOUNDS:
            lo, hi, basis, note = BOUNDS[name]
            claims.bound(name, v, min=lo, max=hi, basis=basis, structural=name in STRUCTURAL, note=note)
        elif name in COUNTS:
            src, thr, mm, note = COUNTS[name]
            claims.count(name, M[src], thr, min_margin=mm, relative=False, note=note)
        elif name in EXACT:
            claims.exact(name, v)
        else:
            claims.value(name, v, rtol=TIGHT.get(name, RTOL))
    n_in_ = int(np.sum(np.abs(M["_r1_resid_u"]) < 1.0))
    claims.count("r1_n_within_u", -np.abs(M["_r1_resid_u"]), -1.0, min_margin=0.2, relative=False,
                 note="points within one read uncertainty of their fit; margin in read units")
    assert n_in_ == claims.metrics["r1_n_within_u"]["value"]
    claims.exact("r1_n_read", M["_r1_n_read"], structural=True, note="number of points read (data size)")
    claims.exact("r5_n_rows", M["_r5_n_rows"], structural=True, note="rows printing all three entries (data)")
    claims.exact("r5_n_pd_rows", M["_r5_n_pd"], structural=True, note="the Pd-based films at the top of Table 2 (a selection)")


register(M)

# %% [markdown]
# ### The break table
#
# Every non-structural number registered above has a row that moves it. Each row injects one
# defect through `Cfg` (a figure read, a printed entry, the profile, the exponent, the
# assumed permeate, the module parameters or the discretisation) and reruns the parts of the
# study it reaches. Structural, with no row: the flux-form identity and the H2 balance
# (round-off identities), and the data sizes.

# %%
def printed_row(item, new, parts):
    """A break row that mistypes one printed entry."""
    pn = PN.copy()
    i = pn.index[pn["item"] == item][0]
    old = pn.at[i, "printed"]
    pn.at[i, "printed"] = new
    return f"printed {item} {old} -> {new}", replace(C0, pn=pn), parts


def t2_row(material, col, new, parts):
    """A break row that mistypes one Table 2 entry."""
    t2 = T2.copy()
    i = t2.index[t2["material"] == material][0]
    old = t2.at[i, col]
    t2.at[i, col] = new
    return f"Table 2 {material} {col} {old} -> {new}", replace(C0, t2=t2), parts


def f15_row(curve, cell, new, parts):
    """A break row that moves one Fig. 15 read."""
    rd = RD15.copy()
    i = rd.index[(rd["curve"] == curve) & (rd["cell"] == cell)][0]
    old = rd.at[i, "rel_H2_wt_fraction"]
    rd.at[i, "rel_H2_wt_fraction"] = new
    return f"Fig. 15 {curve} marker {cell} read {old:.4f} -> {new}", replace(C0, rd15=rd), parts


def f16_row(label, col, P, delta, parts):
    """A break row that shifts Fig. 16 boundary reads (``P`` None: every pressure)."""
    rd = RD16.copy()
    sel = rd.index if P is None else rd.index[rd["P_feed_bar"] == P]
    rd.loc[sel, col] = rd.loc[sel, col] + delta
    return label, replace(C0, rd16=rd), parts


DEFECTS = [
    ("Fig. 15 profile exponent 2 -> 1 (E15)", replace(C0, k15=1.0), ("f15",)),
    f15_row("quadrupled", 12, 0.60, ("f15",)),
    f15_row("experimental", 12, 0.90, ("f15",)),
    f15_row("doubled", 12, 0.80, ("f15",)),
    ("Fig. 15 read uncertainty doubled", replace(C0, u15=2 * U15), ("f15",)),
    f16_row("Fig. 16 Pd share at 15 bar 0.121 -> 0.45 (Pd top moved)", "pd_top", 15,
            0.45 - float(RD16.loc[RD16["P_feed_bar"] == 15, "pd_top"].iloc[0]
                         - RD16.loc[RD16["P_feed_bar"] == 15, "cp_top"].iloc[0]), ("f16",)),
    f16_row("Fig. 16 CP tops all read 0.02 low", "cp_top", None, -0.02, ("f16",)),
    ("Sieverts exponent n = 1 in the Fig. 16 bound", replace(C0, n_mem=1.0), ("f16",)),
    ("Fig. 16 read uncertainty doubled", replace(C0, u16=2 * U16), ("f16",)),
    printed_row("p57_eta_low", "0.8", ("p57",)),
    printed_row("p57_eta_high", "0.95", ("p57",)),
    t2_row("ZIF-8/Al2O3", "thickness_um", "25", ("t2",)),
    t2_row("Pd–Au–Pt", "permeability_1e-13", "135", ("t2",)),
    t2_row("Pd–Ag/PSS", "permeance_1e-8", "3,000", ("t2",)),
    t2_row("ZIF-8/Al2O3", "permeability_1e-13", "60", ("t2",)),
    printed_row("p46_Qs_pure_mantissa", "4.74", ("t2",)),
    ("permeate 1 -> 2 bar (assumption)", replace(C0, p_perm=2.0), ("local", "sweep", "tcmod", "v")),
    ("beta at TC doubled", replace(C0, beta_tc=2.0), ("local", "tcmod", "v")),
    ("target recovery 0.8 -> 0.85", replace(C0, r_star=0.85), ("local", "sweep", "tcmod", "v")),
    ("membrane number at TC 2 -> 2.2", replace(C0, n_tc=2.2), ("local", "v")),
    printed_row("p46_wgs_H2_pct", "50.0", ("local", "sweep", "tcmod", "v")),
    ("module levels halved (n_z 50, 100, 200)", replace(C0, nz0=50), ("sweep", "tcmod", "v")),
    ("trapezoidal source in the module (second order)", replace(C0, source="trapezoid"), ("sweep", "tcmod", "v")),
    ("V2 general route at beta 1e-3 and 1e3", replace(C0, beta_small=1e-3, beta_large=1e3), ("v",)),
    ("V2 module at N = 10 instead of 50", replace(C0, n_large=10.0), ("v",)),
    ("V4 ladder read at the inlet face (fixed by the bc)", replace(C0, v4_read="inlet"), ("v",)),
    ("E14 weighted by the share of the drop, not eps_Q", replace(C0, ea_weight="share"), ("local", "tcmod")),
    ("V5 quadrature route with n = 1 (linear membrane)", replace(C0, n_route2=1.0), ("sweep",)),
]
CONTROLS = [("control: clean rerun of the figure, table and local parts", C0, ("f15", "f16", "p57", "t2", "local"))]
claims.define_rows([r_ for r_, _, _ in CONTROLS + DEFECTS])
BREAK, CRASH = [], {}
for row, c_, parts in CONTROLS + DEFECTS:
    Mb = {}
    for part in parts:                # one part at a time, so a crash costs only its own claims
        try:
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                PARTS[part](c_, Mb)
        except (ValueError, RuntimeError, ZeroDivisionError, FloatingPointError, KeyError) as e:
            CRASH[f"{row} [{part}]"] = f"{type(e).__name__}: {e}"
    if "_r1_resid_u" in Mb:
        Mb["r1_n_within_u"] = int(np.sum(np.abs(Mb["_r1_resid_u"]) < 1.0))
    for name in claims.metrics:
        if claims.metrics[name].get("structural"):
            continue
        x_ = Mb.get(name, None)
        if x_ is None or (isinstance(x_, (float, np.floating)) and not np.isfinite(x_)):
            continue                  # only finite values are recorded
        BREAK.append((row, name, claims.breaks(row, name, x_)))

# %%
# The break table: one line per row, with the claims it fails.
by_row = {}
for r_, m_, c_ in BREAK:
    by_row.setdefault(r_, []).append(c_)
_bad = [r_ for r_ in by_row if r_.startswith("control") and any(by_row[r_])]
assert not _bad, f"a control row failed a claim: {_bad}"
print(f"{'injected defect':<62} {'recorded':>8} {'failed':>7}")
for r_, cs in by_row.items():
    print(f"{r_[:62]:<62} {len(cs):8d} {sum(cs):7d}")
for r_, e_ in CRASH.items():
    print(f"crashed part: {r_}: {e_}")
print(f"{len(CRASH)} of {sum(len(p_) for _, _, p_ in CONTROLS + DEFECTS)} row parts raised "
      f"(each row still records the claims of its other parts)")
_caught = {m_ for r_, m_, c_ in BREAK if c_}
_open = sorted(k_ for k_ in claims.metrics if k_ not in _caught
               and not claims.metrics[k_].get("structural"))
_n_struct = sum(1 for k_ in claims.metrics if claims.metrics[k_].get("structural"))
print(f"\n{len(by_row)} rows; {len(_caught)} of {len(claims.metrics)} claims are failed by at "
      f"least one row; structural: {_n_struct}; open: {', '.join(_open) or 'none'}")

# %% [markdown]
# ## What pymrm adds
#
# The figure verdicts (R1 to R3, R6) and the Table 2 checks (R5) are closed forms and reads:
# pymrm contributes nothing to them beyond `pymrm.checks.find_roots` where a threshold is a
# root. For the area statement (R4) pymrm solves the membrane module: retentate and permeate
# as one finite-volume system (`BlockLayout`, `construct_convflux_upwind`, `construct_div`,
# the exact source derivative through `construct_coefficient_matrix`), the area as a bordering
# unknown in one `pymrm.newton` solve, the outlet read by `compute_boundary_values`, and the
# exact discrete $dN/d\beta$ from the converged Jacobian. `pymrm.checks.observed_orders`
# supplies every grid error and `find_roots` the threshold $\beta_{1/2}$. The module itself
# has a quadrature solution (E12), the second route of R4. On this page the module gives no
# number the quadrature does not.

# %%
say("Stated by the review and not tested here: that the SINTEF work confirms Tiemersma (p. 55, "
    "\"completely confirmed by the work carried out at SINTEF\"; one is a reactor model, the other a "
    "permeation cell), and that fluidized beds reduce the bed-to-wall limitation (p. 55). Not claimed "
    "here: whether Tiemersma's film was limiting.")

# %% [markdown]
# ## Reuse
#
# - To decide whether a film limits a Sieverts membrane, compute $\varepsilon_Q$ of E6 (or its
#   bound E9 from a split of the pressure drop), not the share of the drop.
# - `MembraneModule` takes $\theta$, $\phi$, $\beta$; add a reaction source or axial
#   dispersion through the same operators. For a sweep or counter-current permeate see H1.4.
# - Related pages: H1.1 (Sieverts exponent), H1.3 (concentration polarisation along a tube),
#   H1.4 (membrane reactor, equilibrium shift), H1.7 (solution-diffusion), D1.1 (fixed-bed
#   models).

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
