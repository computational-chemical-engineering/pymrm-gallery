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
# title: "Heck, Wei & Katzer 1976: what the one-dimensional monolith model is adequate for"
# description: "The paper's two claims, tested: Table 1's asymptotic Nusselt numbers derived from scratch, the Nusselt singularity root-found downstream of the cause the text gives it, and the 1-D/2-D adequacy claim re-measured with root-found thresholds."
# categories: [sec:I, struct:S6, tier:T0, data:tier6, phase:gas-solid]
# date: 2026-08-15
# ---
#
# # Heck, Wei & Katzer 1976: what the one-dimensional monolith model is adequate for
#
# **Catalog ID:** `I1.3` · **Structure:** `S6` (2-D PDE) · **Tier:** T0
#
# The abstract makes two claims and the page tests both:
#
# > "A two-dimensional model is shown to predict unusual behavior of the Nusselt
# > number in the presence of rapid reaction. However, a simpler one-dimensional
# > model is adequate for predicting monolith behavior."
#
# *Unusual* compared with what: that is Table 1, the paper's only printed table,
# and this page derives five of its eight cells (two of them in closed form)
# instead of quoting them. *Adequate* for what, and where it stops: that is a
# 1-D versus 2-D comparison the page runs itself, and the paper's own numbers for
# it are **sampled inlet temperatures, not root-found thresholds**.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.01"  # replaced by the release tag at publication
PYMRM = ">=3.3.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "I1.3-heck-wei-katzer-monolith"  # equals the page directory

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
import re
import warnings
from dataclasses import replace

import numpy as np
import pandas as pd
from IPython.display import display_html
from scipy.optimize import brentq, minimize_scalar
from gallery_utils import Claims, load_claims, load_data
from model import (Channel, Graetz, Radial, Voltz, closed_form_nusselt, correlation,
                   cross, duct, fitted_order, local_nusselt, observed_order, parse_flux,
                   reconstruct_groups, richardson, square_spectral, wall_flux)

claims = Claims(PAGE.split("-")[0])


def show(styler):
    """Display a styled table as HTML only, with a pinned uuid, so reruns are identical."""
    display_html(styler.to_html(), raw=True)


# %% [markdown]
# ## Background
#
# A monolith converter is a bundle of parallel ceramic channels with platinum on
# the walls. CO oxidation over Pt is strongly exothermic *and* inhibited by CO, so
# the rate is **negative order** in CO at high CO: the rate curve rises, peaks and
# falls as the wall concentration increases. Put that curve against a straight
# mass-transfer line and you get one or three intersections, and the jump from the
# outer to the inner one is **light-off**.
#
# The question the paper asks is a modelling question, not a chemistry one. A 1-D
# model lumps the cross-section and closes the wall fluxes with a Nusselt and a
# Sherwood number. A 2-D model resolves the radial profiles and *computes* the
# Nusselt number from them. Young and Finlayson (1974) had reported that the
# Nusselt number behaves strangely when the reaction is fast, and concluded the
# expensive 2-D model was necessary. Heck, Wei and Katzer confirm the strange
# behaviour and then argue the 1-D model is good enough anyway.
#
# Both halves are checkable **without touching a figure**, because the paper
# states the numbers that matter in running text.
#
# `I1.2` (Oh & Cavendish 1982) is the neighbouring page: the same reactor, the
# same Voltz kinetics family, transient rather than steady, and reproduced against
# a printed table. This page shares no data with it and loads none of its files.

# %% [markdown]
# ## The published model
#
# Read off the 300 ppi page bitmaps (`pdfimages -list` reports 300 ppi CCITT-G4
# bilevel on all eight pages, so 300 is native). Nothing below came from the PDF
# text layer. The printed constants are listed from the data file in *The data*.
#
# ### The kinetics (p. 478), from Voltz et al. (1973)
#
# $$r=\frac{k_r\,C\,C_{\mathrm{O_2}}}{[1+k_a C]^{2}},\qquad
# k_r=k_r^{o}\exp\!\left[-\frac{E_r}{R\,(T_w+273)}\right],\qquad
# k_a=k_a^{o}\exp\!\left[-\frac{E_a}{R\,(T_w+273)}\right]$$
#
# with $C$ a **mole fraction** (Notation, p. 483) and $E_a < 0$ (adsorption). The
# paper writes no oxygen balance, so $C_{\mathrm{O_2}}$ is held at its inlet
# value here and stoichiometric depletion is a break-table row.
#
# ### The one-dimensional model (eqs. 8 to 14)
#
# $$\frac{dC_G}{dX^{\ast}}=\frac{Sh}{Le}(C_w-C_G),\qquad
# \frac{dT_G}{dX^{\ast}}=Nu\,(T_w-T_G),$$
#
# closed by $-r\Delta H = h(T_w-T_G)$ and $r=\frac{\rho_G k_m}{M_G}(C_G-C_w)$, and
# reduced with the adiabatic relation to a single algebraic wall condition. With
# $Le\cong1$ and $Sh\cong Nu$ the paper's eq. (14) is
#
# $$T_w \cong T_G^{o} + \Delta T_{\mathrm{AD}}\left(1-\frac{C_w}{C_G^{o}}\right),
# \qquad \Delta T_{\mathrm{AD}}=\frac{\Delta H\,C_G^{o}}{M_G C_{PG}} .$$
#
# **Equation (14) is not printed that way.** Verbatim, it reads
#
# > $T_w \cong T_G{}^{o} + \Delta T_{AD} = \left(1 - \dfrac{C_w}{C_G{}^{o}}\right)$
# > (14)
#
# with an "=" where a multiplication is required. As printed it sets a temperature
# equal to a dimensionless number, so it cannot be right; the form above is what
# eq. (13) gives when $Le\cong1$ and $Sh\cong Nu$ are imposed as the text
# instructs. That is a derivation, labelled an inference. Reported, not repaired.
#
# ### The two-dimensional model (eqs. 15 to 19)
#
# $$\frac{(1-R^{2})}{2}\frac{\partial T_G}{\partial X^{\ast}}
# =\frac{1}{R}\frac{\partial}{\partial R}\left(R\frac{\partial T_G}{\partial R}\right),
# \qquad
# \frac{(1-R^{2})}{2}\frac{\partial C_G}{\partial X^{\ast}}
# =\frac{1}{Le\,R}\frac{\partial}{\partial R}\left(R\frac{\partial C_G}{\partial R}\right),$$
#
# $$r\Delta H=\frac{2k_G}{D}\left.\frac{\partial T_G}{\partial R}\right|_w,
# \qquad r=\frac{2\rho_G\mathcal{D}}{M_G D}\left.\frac{\partial C_G}{\partial R}\right|_w,
# \qquad
# Nu=\frac{2\,\partial T/\partial R|_w}{T_w-\bar T_G}. $$
#
# Two things about this block carry everything below.
#
# **First, $X^{\ast}$ does not mean the same thing in eq. (9) and eq. (15).** The
# Notation defines $X^{\ast}=4x/(D\,Re\,Pr)$. Equation (9) requires exactly that (it is
# the 1-D energy balance with $Nu=hD/k_G$). Equation (15) requires
# $X^{\ast}=x/(D\,Re\,Pr)$, a factor of four smaller, and the figure abscissae are
# labelled `X/(D RE PR)` without the four. The paper's own printed check value
# settles which the 2-D results use (Validation, first section). A **third**
# scaling, $1000\,x/(D\,Re\,Pr)$, is required by the Grigull-Tratz correlations
# the paper reprints on p. 482, shown against a solve further down.
#
# **Second, eqs. (17) and (18) cannot both have the sign they are printed with.**
# The wall is a *sink* for CO and a *source* of heat, so $\partial C_G/\partial R|_w$
# and $\partial T_G/\partial R|_w$ must have opposite signs; as printed both
# equations set a positive $r$ proportional to a positive gradient. The same applies
# to eq. (10) against eq. (12): $\Delta H$ cannot be positive in one and negative in
# the other. These are sign-convention slips with an obvious intent, reported, not
# repaired: the physically consistent signs are used.

# %% [markdown]
# ## Parameters and assumptions
#
# Everything the paper prints is in `data/heck-1976-printed-values.csv`, with the
# journal page and the location for each row: running text, a displayed equation, a
# figure caption, or a character printed **inside** a figure's frame (an annotation
# box, an abscissa tick label, a curve label or a case label). **No curve was
# digitised anywhere on this page.**
#
# ### What the paper does not print
#
# The reacting model needs two dimensional groups the paper never states:
#
# * $\Delta T_{\mathrm{AD}} = \Delta H\,C_G^{o}/(M_G C_{PG})$: $\Delta H$ and
#   $C_G^{o}$ are printed, $M_G C_{PG}$ is not;
# * $\beta \equiv \rho_G k_m/M_G = (P/R_g T)\,Sh\,\mathcal{D}/D$: $Sh$ and $D$
#   are printed, the molar density and the diffusivity are not.
#
# Neither can be recovered from the paper without circularity, so both are
# **reconstructed from ordinary gas properties** (an exhaust gas taken as air, the
# CO-air diffusivity scaled with temperature, atmospheric pressure, properties at
# the inlet temperature of Fig. 1; the values are printed below). The 1-D and the 2-D model use *the same two numbers*, so the comparison
# between them, which is what the page is about, is far less sensitive to the
# reconstruction than either model alone; the Validation section measures how much
# less.
#
# ### Scope: figures are out
#
# Figures 1 to 13 carry the paper's profiles. Comparing against them would need
# digitisation, which needs a maintainer review that is not available, so **it is
# scoped out**: nothing is taken from a figure except characters *printed* in it,
# and every such row is listed and counted below. The page cannot check any profile
# shape, any light-off *position* against the paper's own, or the conversion curves
# of Figs. 4 to 6 and 12. What it can do instead is check every number the paper
# states in words.

# %%
T1 = load_data("heck-1976-table1-nusselt.csv", page=PAGE)
PV = load_data("heck-1976-printed-values.csv", page=PAGE)
P = {r.quantity: r.value for r in PV.itertuples()}         # printed values, by name


def pf(name):
    """A printed value as a float."""
    return float(P[name])


def t1(geom, col):
    """A Table 1 cell as printed."""
    return float(T1.loc[T1.geometry == geom, col].iloc[0])


# every numeric printed value is registered as read: the page quotes them
for r in PV.itertuples():
    try:
        claims.exact(f"paper_{r.quantity}", float(r.value), structural=True,
                     note="read from heck-1976-printed-values.csv")
    except ValueError:
        pass                                  # the four printed expressions
for r in T1.itertuples():
    for col in ("Nu_H", "Nu_T"):
        claims.exact(f"paper_table1_{r.geometry.replace(' ', '_')}_{col}",
                     float(getattr(r, col)), structural=True,
                     note="read from heck-1976-table1-nusselt.csv")

# the reconstruction (model.reconstruct_groups holds the gas properties)
CG0, NU_SQ, T_IN = pf("CG_CO_inlet"), pf("Nu_used_1D"), pf("TG_inlet_fig1")
NU_CI = t1("circular", "Nu_H")                # round-tube value for the shape split
DTAD, BETA, B0 = reconstruct_groups(pf("delta_H") * 1e6, CG0, pf("channel_diameter") * 1e-3,
                                    NU_SQ, T_IN + 273.0)
KIN = Voltz(pf("k_r0"), pf("k_a0"), pf("Er_over_R"), pf("Ea_over_R"), CG0,
            pf("CG_O2_inlet"))
CH = Channel(KIN, DTAD, B0)
print(f"reconstructed: dT_AD = {DTAD:.2f} K from the printed dH = {P['delta_H']} MJ/kmol "
      f"and C_G0 = {CG0}")
print(f"               beta  = {BETA:.4e} kmol/m2 s per unit mole fraction, "
      f"(P/RgT) Sh D_AB / D with Sh = {NU_SQ}")
import inspect
GAS = {k: v.default for k, v in inspect.signature(reconstruct_groups).parameters.items()
       if k in ("M_G", "C_PG", "P", "R_g", "D_AB300", "n_T")}
print("gas properties: " + ", ".join(f"{k} = {v:g}" for k, v in GAS.items())
      + "  (kg/kmol, J/kg K, Pa, J/kmol K, m2/s at 300 K, exponent of T)")
for k, v in GAS.items():
    claims.exact(f"input_{k}", float(v), structural=True, note="model.reconstruct_groups default")
claims.exact("deferred_omega", 0.6, structural=True, note="Graetz.solve default relaxation")
claims.exact("sweep_L_margin", 0.30, structural=True, note="chosen sweep end past the paper's lengths")
claims.value("dTad_reconstructed", DTAD, rtol=1e-9)
claims.value("beta_reconstructed", BETA, rtol=1e-9)

# %% [markdown]
# ## The data
#
# Two files, both tier 6 (the paper's own printed values, no measurement).
#
# * `heck-1976-table1-nusselt.csv`: Table 1, all eight cells, the paper's **only**
#   printed table.
# * `heck-1976-printed-values.csv`: every constant and every stated numerical
#   result, with its journal page.
#
# Table 1 is not the authors' computation: its footnote attributes it to Shah and
# London (1971), which was **not consulted**. Two of its cells have exact closed
# forms, derived here. The sinusoidal row cannot be checked from this paper at
# all: the geometry column prints only a glyph, and the paper states no aspect
# ratio or wall profile.
#
# The cell below lists every row that names a figure and **counts them from the
# CSV itself**, so the scope statement above can be audited without diffing the
# file.

# %%
show(T1.style.set_uuid("i13_t1").hide(axis="index"))
IN_FRAME = (("annotation", "annotation box"), ("tick label", "abscissa tick labels"),
            ("curve label", "curve label"), (r"case \d+ label", "case label"))
OUTSIDE = (("caption", "figure caption"), ("running text", "running text"),
           ("displayed", "displayed equation"))


def _kinds(w, pats):
    return ", ".join(lab for pat, lab in pats if re.search(pat, w)) or "-"


prov = pd.DataFrame([(r.quantity, r.where, _kinds(r.where, IN_FRAME), _kinds(r.where, OUTSIDE))
                     for r in PV.itertuples() if re.search(r"Fig", r.where)],
                    columns=["row", "where (as recorded in the CSV)",
                             "read INSIDE the figure frame", "also printed outside it"])
show(prov.style.set_uuid("i13_prov").hide(axis="index"))
inside = prov["read INSIDE the figure frame"] != "-"
n_in = int(inside.sum())
n_only = int((inside & (prov["also printed outside it"] == "-")).sum())
by_kind = {lab: int(prov["read INSIDE the figure frame"].str.contains(lab, regex=False).sum())
           for _, lab in IN_FRAME}
print(f"{len(PV)} printed constants and stated results; {len(T1)} rows in Table 1.")
print(f"{len(prov)} rows name a figure; {n_in} were read inside a figure's frame, {n_only} of "
      f"them\nprinted nowhere else in the paper; the other {len(prov) - n_in} are in running "
      "text or a caption.")
print("By kind: " + ", ".join(f"{lab} {n}" for lab, n in by_kind.items()) + ".")
show(PV[["quantity", "value", "unit"]].style.set_uuid("i13_pv").hide(axis="index"))
for name, n in (("data_rows", len(PV)), ("table1_rows", len(T1)), ("rows_naming_a_figure", len(prov)),
                ("rows_inside_a_frame", n_in), ("rows_only_in_a_frame", n_only)):
    claims.exact(name, n, structural=True, note="counted from the CSV")
for lab, n in by_kind.items():
    claims.exact("rows_" + lab.replace(" ", "_"), n, structural=True, note="counted from the CSV")

# %% [markdown]
# ## PyMRM implementation
#
# Three solvers, in `model.py`, all built from `pymrm` operators that are
# assembled once per solve.
#
# 1. **`duct(...)`**: fully developed flow and heat transfer in a duct
#    cross-section, for Table 1. Cartesian 2-D for the square (quarter domain,
#    symmetry planes as zero-gradient boundaries), cylindrical 1-D for the circle
#    (`nu=1` in `construct_div`). The constant-wall-temperature case is the
#    $u$-weighted eigenproblem, solved by shift-invert `eigsh`.
# 2. **`Graetz`**: the developing problem in a round tube on the whole $(X, R)$
#    domain, shape `(n_x, n_r)`, with van Leer face values in $X$. The axial
#    direction carries convection only, so the outlet takes no condition:
#    `{"outflow": True}`. The radial wall takes an imposed flux that varies along
#    $X$ (an array-valued `d`) or a fixed temperature, and wall values come from
#    `compute_boundary_values`. The production route is deferred correction on one
#    LU of the upwind operator, under-relaxed with $\omega = 0.6$. An independent
#    route, `solve_newton`, solves the same discrete equations with the exact van
#    Leer Jacobian of `ConvFluxTVD` and `newton(damping="monotone")`; the
#    Validation section compares the two.
# 3. **`Channel`**: the reacting monolith channel, 1-D and 2-D on a common axial
#    coordinate $G = x/(D\,Re\,Pr)$ and on the same reconstructed groups. With
#    $Le = 1$ the 2-D field is linear in the wall flux $q$, so the response to a
#    unit flux is computed once per step length and each backward-Euler step is one
#    linear solve (on a matrix from `construct_coefficient_matrix` and the radial
#    operators) and then a *scalar* root for the wall value: the same algebra as the
#    1-D mass-transfer line, with an effective $\beta$. The cold branch exists
#    while the bulk value lies above the fold, $r'(C_w) = -\beta$; light-off is
#    root-found **inside the crossing step**, on the length of that step, not
#    reported as the step itself.
#
# Every switch a break-table row needs is a **function argument**, so each row
# re-runs the computation with the defect in place.

# %% [markdown]
# ## Results
#
# ### 1. Table 1, derived rather than quoted
#
# The paper's only table is the "usual" behaviour against which the Nusselt number
# is called unusual. Two of its cells have exact closed forms; the square pair
# does not, and is computed twice on discretisations that share nothing.

# %%
nu_c_exact, nu_t_exact = closed_form_nusselt()
print(f"circular duct,   constant wall flux: Nu_H1 = {nu_c_exact} = {float(nu_c_exact):.6f}"
      f"   (Table 1 prints {t1('circular', 'Nu_H')})")
print(f"equilateral tri, constant wall flux: Nu_H1 = {nu_t_exact} = {float(nu_t_exact):.6f}"
      f"   (Table 1 prints {t1('equilateral triangle', 'Nu_H')})")
print("Both from exact polynomial solutions of the cross-section problem; nothing is quoted.")
claims.exact("nu_h1_circle_exact", float(nu_c_exact), structural=True,
             note="48/11, a symbolic result; no solver argument produces it")
claims.exact("nu_h1_triangle_exact", float(nu_t_exact), structural=True,
             note="28/9, a symbolic result")


def printed_relerr(exact, printed):
    return abs(exact - printed) / exact


claims.value("nu_h1_circle_printed_relerr",
             printed_relerr(float(nu_c_exact), t1("circular", "Nu_H")), rtol=1e-6)
claims.value("nu_h1_triangle_printed_relerr",
             printed_relerr(float(nu_t_exact), t1("equilateral triangle", "Nu_H")), rtol=1e-6)

# %%


def duct_ladders(**kw):
    """Circle and square Nusselt numbers and fRe, Richardson-extrapolated over four rungs."""
    lad_c = [duct(n, geom="circle", **kw) for n in (64, 128, 256, 512)]
    lad_s = [duct(n, geom="square", **{k: v for k, v in kw.items() if k != "nu_r"})
             for n in (32, 64, 128, 256)]
    out = {}
    for tag, lad in (("circle", lad_c), ("square", lad_s)):
        for q in ("Nu_H1", "Nu_T", "fRe"):
            out[f"{tag}_{q}"] = richardson([d[q] for d in lad])
    return out


DU = duct_ladders()
print(f"circle  fRe -> {DU['circle_fRe']:.6f}   (Poiseuille: 16)")
print(f"square  fRe -> {DU['square_fRe']:.6f}")
A315 = load_claims("A3.15-graetz-nusselt")["nu_infinity"]
tab1 = pd.DataFrame([
    ("circular", "Nu_H", t1("circular", "Nu_H"), DU["circle_Nu_H1"], float(nu_c_exact)),
    ("circular", "Nu_T", t1("circular", "Nu_T"), DU["circle_Nu_T"], np.nan),
    ("square", "Nu_H", t1("square", "Nu_H"), DU["square_Nu_H1"], np.nan),
    ("square", "Nu_T", t1("square", "Nu_T"), DU["square_Nu_T"], np.nan),
    ("equilateral triangle", "Nu_H", t1("equilateral triangle", "Nu_H"), np.nan,
     float(nu_t_exact))],
    columns=["geometry", "column", "printed", "pymrm (extrapolated)", "closed form"])
ref_col = tab1["pymrm (extrapolated)"].fillna(tab1["closed form"])
tab1["rel. dev."] = abs(tab1["printed"] - ref_col) / ref_col
tab1["half a printed unit, rel."] = 5e-4 / ref_col
show(tab1.style.set_uuid("i13_table1").hide(axis="index").format(
    {"printed": "{:.3f}", "pymrm (extrapolated)": "{:.6f}", "closed form": "{:.6f}",
     "rel. dev.": "{:.2e}", "half a printed unit, rel.": "{:.2e}"}))
n_agree = int((tab1["rel. dev."] <= tab1["half a printed unit, rel."]).sum())
print(f"{n_agree} of {len(tab1)} reproduced cells agree with every printed digit.")
print(f"circular Nu_T against A3.15's committed Graetz limit {A315:.8f}: relative gap "
      f"(here {abs(DU['circle_Nu_T'] / A315 - 1):.1e}), held below 1e-6.")


def table1_claims(du, nu_h_circle_exact=float(nu_c_exact), col_h="Nu_H", col_t="Nu_T"):
    """The Table 1 claims from a set of duct ladders (columns as read from the CSV)."""
    return {
        "nu_h1_square_pymrm": du["square_Nu_H1"], "nu_t_square_pymrm": du["square_Nu_T"],
        "nu_t_circle_pymrm": du["circle_Nu_T"], "nu_h1_circle_pymrm": du["circle_Nu_H1"],
        "fre_square": du["square_fRe"],
        "fre_circle_vs_16": abs(du["circle_fRe"] / 16.0 - 1.0),
        "nu_h1_circle_solve_vs_exact": printed_relerr(nu_h_circle_exact, du["circle_Nu_H1"]),
        "nu_h1_square_printed_relerr": printed_relerr(du["square_Nu_H1"], t1("square", col_h)),
        "nu_t_square_printed_relerr": printed_relerr(du["square_Nu_T"], t1("square", col_t)),
        "nu_t_circle_printed_relerr": printed_relerr(du["circle_Nu_T"], t1("circular", col_t)),
        "nu_t_circle_vs_A315": abs(du["circle_Nu_T"] / A315 - 1.0),
        "table1_cells_within_printed_digits": int(sum(
            abs(t1(g, col) - v) <= 5e-4 for g, col, v in (
                ("circular", col_h, nu_h_circle_exact), ("circular", col_t, du["circle_Nu_T"]),
                ("square", col_h, du["square_Nu_H1"]), ("square", col_t, du["square_Nu_T"]),
                ("equilateral triangle", col_h, float(nu_t_exact))))),
    }


TC = table1_claims(DU)
for k, v in TC.items():
    if k == "fre_circle_vs_16":
        claims.bound(k, v, max=1e-6, basis="theory", note="Poiseuille fRe = 16, Richardson")
    elif k == "nu_h1_circle_solve_vs_exact":
        claims.bound(k, v, max=1e-7, basis="theory", note="Richardson-extrapolated against 48/11")
    elif k == "nu_t_circle_vs_A315":
        claims.bound(k, v, max=1e-6, basis="chosen",
                     note="another page's committed value, independent route")
    elif k == "table1_cells_within_printed_digits":
        assert v == n_agree
        claims.exact(k, v, margin=float(np.min(abs(tab1["rel. dev."] - tab1["half a printed unit, rel."])
                                               / tab1["half a printed unit, rel."])),
                     min_margin=0.05, note="margin: smallest distance to the half-unit edge, relative")
    elif k.endswith("printed_relerr"):
        claims.value(k, v, rtol=1e-5)
    else:
        claims.value(k, v, rtol=1e-6)

# %%
SPEC_N = (12, 16, 20, 24)


def spectral_ladder(**kw):
    return [square_spectral(N, **kw) for N in SPEC_N]


SPEC = spectral_ladder()
print("Fourier-sine Galerkin, square duct, Nu_T (converges from below as observed here): "
      + ", ".join(f"N={N} {b:.6f}" for N, (_, b) in zip(SPEC_N, SPEC)))
d_spec = np.diff([b for _, b in SPEC])
print(f"pymrm finite volume, extrapolated: Nu_H1 = {DU['square_Nu_H1']:.6f}, "
      f"Nu_T = {DU['square_Nu_T']:.6f}")
two_routes = abs(SPEC[-1][0] - DU["square_Nu_H1"]) / DU["square_Nu_H1"]
print(f"""
The spectral sequence is not a bound: the sine-Galerkin eigenproblem is a
Rayleigh-Ritz restriction, so it bounds the eigenvalue from above, while the
truncated integral of u converges from below, and their product has no one-sided
guarantee. What the sequence shows is that it is still increasing at N = {SPEC_N[-1]} (last step
{d_spec[-1]:+.1e}),
where it has already reached {SPEC[-1][1]:.6f}, above the interval
[{t1('square', 'Nu_T') - 5e-4:.4f}, {t1('square', 'Nu_T') + 5e-4:.4f}] the printed {t1('square', 'Nu_T'):.3f} allows. The two routes
agree to {two_routes:.1e} relative on Nu_H1 at N = {SPEC_N[-1]}.""")
for N, (_, b) in zip(SPEC_N, SPEC):
    claims.value(f"nu_t_square_spectral_N{N}", b, rtol=1e-6)
assert np.all(d_spec > 0)
claims.value("nu_t_square_spectral_last_increment", d_spec[-1], rtol=1e-4,
             note="N = 20 to 24; every increment of the ladder is positive")
claims.value("nu_h1_square_two_routes_reldiff", two_routes, rtol=1e-4)
claims.value("nu_t_square_spectral_above_interval", SPEC[-1][1] - (t1("square", "Nu_T") + 5e-4),
             rtol=1e-4, note="N = 24 value minus the upper edge of the printed interval")

# %% [markdown]
# **Five of the eight cells are reproduced**, and the count of those that agree
# with every printed digit is printed above. The square duct at constant wall
# temperature does not: it is printed **2.976** where both routes give **2.9775**.
# `2.976` implies the interval [2.9755, 2.9765], which excludes 2.9775, so this is
# not a rounding difference. It is also not a defect this page can pin on Heck et
# al.: Table 1 is attributed to Shah and London (1971), which was **not
# consulted**, so the discrepancy could be in the 1971 tabulation or in the
# transcription of it, and the page does not choose. (The paper's own reference
# list prints that source as "**Shar**, R. K., and A. L. London" where Table 1's
# footnote says "Shah and London": reported, not repaired.)
#
# **Three cells are not reproduced.** The equilateral triangle at constant wall
# temperature needs a $u$-weighted eigenproblem on a non-tensor domain, which this
# page does not build. Both sinusoidal cells are **under-determined by the
# paper**: the geometry column prints a glyph and no aspect ratio or wall profile
# anywhere, so there is no problem to solve.
#
# ### 2. The unusual Nusselt number, and what actually causes it
#
# The paper demonstrates the anomaly on the *non-reacting* Graetz problem with an
# imposed wall flux that rises through zero, and it prints both the flux functions
# and the landmark positions, in running text, in Fig. 8's curve labels and in
# Fig. 9's caption (printed characters, no curve read; the two functions are read
# from the CSV below):
#
# * "the wall flux rises rapidly from a negative value and passes through zero near
#   $X^{\ast}=0.029$. As the wall flux passes through zero, infinite negative and
#   positive values for the Nusselt number occur."
# * "When $X^{\ast}=0.032$ ... the average gas temperature $T_G^{\ast}$ is higher than
#   $T_w^{\ast}$; thus a large negative value results ... At $X^{\ast}=0.033$,
#   $(\partial T_G^{\ast}/\partial R)_w$ and $(T_w^{\ast}-\bar T_G^{\ast})$ are both
#   positive."
#
# So the paper states one location for the singularity in the first quotation and
# brackets a different one in the second. That is decidable by computation.

# %%
QW = {"Fig. 8": parse_flux(P["Qw_fig8"]), "Fig. 9": parse_flux(P["Qw_fig9"])}
X_SIGN = (pf("X_sign_negative_Nu"), pf("X_sign_positive_Nu"))   # the p. 481 sign statements
GZ_PROD = dict(X=0.05, n_x=6000, n_r=200, r_grade=1.0)


def anomaly(coef, grid=GZ_PROD, **solve_kw):
    """Flux zero, Nusselt zero and pole of one printed wall flux, all root-found."""
    Q = wall_flux(*coef)
    s = Graetz(**grid).solve(("flux", Q), **solve_kw)
    zero_q = brentq(Q, 1e-4, 0.05)
    pole = cross(s["x"], s["Tw"] - s["cup"])
    nu = local_nusselt(s)
    pre = s["x"] < pole                    # before the pole Nu is continuous
    return dict(s=s, zero_q=zero_q, zero_nu=cross(s["x"][pre], nu[pre]), pole=pole,
                pole_area=cross(s["x"], s["Tw"] - s["area"]),
                nu_at_zero_q=float(np.interp(zero_q, s["x"], nu)))


AN = {nm: anomaly(coef) for nm, coef in QW.items()}
for nm, a in AN.items():
    print(f"{nm}  Qw = 0 at X* = {a['zero_q']:.6f}   (text: 'near {P['Qw_zero_stated']}')")
    print(f"        Nu there (here {a['nu_at_zero_q']:+.1e}), so it passes through zero, not "
          "infinity")
    print(f"        Nu = 0 at X* = {a['zero_nu']:.6f}, the same point, root-found on Nu")
    print(f"        Nu pole at X* = {a['pole']:.6f} (mixing-cup mean), "
          f"{a['pole_area']:.6f} (area mean)")
    print(f"        the pole sits {100 * (a['pole'] / a['zero_q'] - 1):.2f} % downstream of the "
          "flux zero\n")
a9 = AN["Fig. 9"]

# the paper's own ten radial points, reused below
G10 = dict(X=0.05, n_x=3200, n_r=int(pf("radial_points_2D")), r_grade=1.0)
pole_10r = anomaly(QW["Fig. 9"], G10)["pole"]

# normalisation: the problem is linear in T, so halving Q_w leaves the pole
gh = Graetz(X=0.05, n_x=3200, n_r=40, r_grade=1.0)
a, b, p = QW["Fig. 9"]
sh1, sh2 = gh.solve(("flux", wall_flux(a, b, p))), gh.solve(
    ("flux", lambda X: 0.5 * wall_flux(a, b, p)(X)))
pole_gap = abs(cross(sh1["x"], sh1["Tw"] - sh1["cup"]) - cross(sh2["x"], sh2["Tw"] - sh2["cup"]))
print(f"Q_w and Q_w/2 on one grid: poles differ by (here {pole_gap:.1e}), below 1e-12.")


def anomaly_claims(an, p10):
    a8, a9 = an["Fig. 8"], an["Fig. 9"]
    return {"fluxzero_X_fig8": a8["zero_q"], "fluxzero_X_fig9": a9["zero_q"],
            "nu_zero_X_fig8": a8["zero_nu"], "nu_zero_X_fig9": a9["zero_nu"],
            "nu_pole_X_fig8": a8["pole"], "nu_pole_X_fig9": a9["pole"],
            "nu_pole_area_X_fig8": a8["pole_area"], "nu_pole_area_X_fig9": a9["pole_area"],
            "nu_pole_over_fluxzero_fig8": a8["pole"] / a8["zero_q"],
            "nu_pole_over_fluxzero_fig9": a9["pole"] / a9["zero_q"],
            "nu_pole_downstream_fig8": a8["pole"] / a8["zero_q"] - 1.0,
            "nu_pole_downstream_fig9": a9["pole"] / a9["zero_q"] - 1.0,
            "nu_pole_below_bracket_centre_pct": 100 * (1 - a9["pole"] / np.mean(X_SIGN)),
            "nu_pole_X_fig9_n_r10": p10,
            "pole_n_r10_dev_pct": 100 * abs(p10 / a9["pole"] - 1),
            "pole_n_r10_short_of_0p032_pct": 100 * (1 - p10 / X_SIGN[0])}


ANC = anomaly_claims(AN, pole_10r)
for k, v in ANC.items():
    claims.value(k, v, rtol=1e-6 if "pct" not in k else 1e-5)
for nm in AN:
    tag = nm.replace("Fig. ", "fig")
    claims.bound(f"nu_zero_minus_fluxzero_{tag}", abs(AN[nm]["zero_nu"] - AN[nm]["zero_q"]),
                 max=1e-7, basis="identity", structural=True,
                 note="Nu's numerator is the imposed flux; the gap is the interpolation of Nu")
    claims.bound(f"nu_at_fluxzero_{tag}", abs(AN[nm]["nu_at_zero_q"]), max=1e-3,
                 basis="identity", structural=True, note="interpolated Nu at the flux zero")
claims.bound("pole_normalisation_gap", pole_gap, max=1e-12, basis="identity", structural=True,
             note="the problem is linear in T")

# %%
s = a9["s"]


def sign_rows(s):
    rows = []
    for xq in (0.025,) + X_SIGN:
        i = int(np.argmin(abs(s["x"] - xq)))
        d, q = s["Tw"][i] - s["cup"][i], s["q"][i]
        rows.append((s["x"][i], q, d, 2 * q / d))
    return pd.DataFrame(rows, columns=["X*", "Qw = dT/dR|w", "Tw - Tbar", "Nu"])


SG = sign_rows(s)
SG["paper says"] = ["cooling, Nu positive and ordinary",
                    "dT/dR|w > 0 but Tbar > Tw, so LARGE NEGATIVE Nu",
                    "dT/dR|w > 0 and Tw > Tbar, so large positive Nu"]
show(SG.style.set_uuid("i13_signs").hide(axis="index").format(
    {"X*": "{:.4f}", "Qw = dT/dR|w": "{:+.4f}", "Tw - Tbar": "{:+.5f}", "Nu": "{:+.3f}"}))


def sign_claims(sg):
    return {"fig9_sign_0p032_matches_paper": bool(sg["Tw - Tbar"].iloc[1] < 0),
            "fig9_sign_0p033_matches_paper": bool(sg["Tw - Tbar"].iloc[2] > 0)}


SC = sign_claims(SG)
for k, v in SC.items():
    claims.exact(k, v, margin=float(abs(SG["Tw - Tbar"].iloc[1 if "032" in k else 2])),
                 min_margin=0.01, note="margin: |Tw - Tbar| at that cell")
print(f"""The X* = {X_SIGN[1]} statement is reproduced; the X* = {X_SIGN[0]} statement is not. The pole is
root-found at X* = {a9['pole']:.6f}, so by X* = {X_SIGN[0]} the solve has already passed it and
(Tw - Tbar) is positive. The paper's two sentences bracket the pole in ({X_SIGN[0]}, {X_SIGN[1]});
this solve puts it {ANC['nu_pole_below_bracket_centre_pct']:.1f} % below the centre of that bracket. On the paper's own
{G10['n_r']} radial points the same solver gives X* = {pole_10r:.6f}, {ANC['pole_n_r10_dev_pct']:.2f} % from the converged
value and still {ANC['pole_n_r10_short_of_0p032_pct']:.1f} % short of {X_SIGN[0]}.""")

# %% [markdown]
# **The cause the text gives is not the cause.** At the zero of the wall flux the
# Nusselt number of eq. (19) has a vanishing *numerator*, so it passes through
# **zero**, not infinity, which is what the computation shows. The pole is the zero
# of the *denominator* $T_w-\bar T_G$, and it lies **8.7 %** further downstream in
# both of the paper's flux functions. The two are close but distinct, and they must
# be: the bulk temperature integrates the flux history, so it keeps falling for a
# while after the flux turns positive.
#
# This is not an effect the authors missed. The very next paragraph, the Fig. 9
# discussion quoted above, gives the correct mechanism in terms of the signs of
# the gradient and of $(T_w^{\ast}-\bar T_G^{\ast})$. What is wrong is the one-sentence
# summary, and the bracket that paragraph implies for the pole.
#
# **The anomaly itself is real and reproduces.** It is also
# *normalisation-invariant*: the problem is linear in $T$, so scaling $Q_w$ scales
# the field and leaves the pole where it is. Whatever convention the paper used for
# the magnitude of $Q_w$, which it never states, the landmark positions are
# unchanged.
#
# ### 3. "A simpler one-dimensional model is adequate": for what, and until where
#
# The reacting comparison needs the two reconstructed groups. **Nothing below is
# fitted to any result of the paper.** The paper prints one fact that touches the
# reconstruction, the position of the rate maximum on the adiabatic locus, and
# the next cell measures how little it constrains.

# %%
def rate_claims(ch, scan=True):
    """The rate maximum on the adiabatic locus at Fig. 1's inlet, and what it constrains.

    ``scan=False`` skips the range of dT_AD that still rounds to the printed value.
    """
    def cw_at_max(d):
        c = replace(ch, dTad=d)
        return float(minimize_scalar(lambda x: -c.r_of(x, T_IN), bounds=(1e-6, c.cg0),
                                     method="bounded", options=dict(xatol=1e-12)).x)
    tgt = pf("Cw_at_rate_maximum")
    half_unit = 0.5 * 10.0 ** np.floor(np.log10(tgt))     # the paper prints one figure
    cmax = cw_at_max(ch.dTad)
    ka = ch.kin.ka0 * np.exp(-ch.kin.ear / (T_IN + 273.0))
    out = {"rate_max_Cw": cmax, "rate_max_Cw_printed_absdiff": abs(cmax - tgt),
           "rate_max_at_dTad0": cw_at_max(0.0), "inverse_ka_at_inlet": 1.0 / ka,
           "lit_wall_T_C": T_IN + ch.dTad, "rate_max_half_unit": half_unit}
    if scan:
        d_hi = brentq(lambda d: cw_at_max(d) - (tgt + half_unit), 0.5 * ch.dTad,
                      8.0 * ch.dTad, xtol=1e-6)
        lo = minimize_scalar(cw_at_max, bracket=(0.1 * ch.dTad, 0.5 * ch.dTad, ch.dTad),
                             method="brent", options=dict(xtol=1e-10))
        out.update(rate_max_dTad_upper_K=d_hi, rate_max_min_Cw=float(lo.fun),
                   rate_max_min_at_dTad_K=float(lo.x),
                   rate_max_min_clearance=float(lo.fun) - (tgt - half_unit))
    return out


RC = rate_claims(CH)
tgt = pf("Cw_at_rate_maximum")
print(f"""The paper states 'the rate reaches a maximum at C_w = {tgt}' for T_G0 = {T_IN:.0f} C. The
maximum of eq. (5) on the eq. (14) locus is at C_w = {RC['rate_max_Cw']:.5f}, which rounds to the
one significant figure printed.

How much does that check? C_w at the maximum barely depends on dT_AD: it stays inside
[{tgt - RC['rate_max_half_unit']:g}, {tgt + RC['rate_max_half_unit']:g}], so still rounds to {tgt}, for every dT_AD from 0 up to {RC['rate_max_dTad_upper_K']:.0f} K
(root-found); its minimum over that range is {RC['rate_max_min_Cw']:.6f} at dT_AD = {RC['rate_max_min_at_dTad_K']:.0f} K,
{RC['rate_max_min_clearance']:.1e} inside the lower edge. So it does not check dT_AD. It checks the transcription
of the kinetics: at dT_AD = 0 the maximum is at {RC['rate_max_at_dTad0']:.5f} = 1/k_a(T_G0) = {RC['inverse_ka_at_inlet']:.5f},
and k_a0 and E_a/R are printed. beta does not enter it at all.

The same dT_AD puts a fully lit wall at T_G0 + dT_AD = {RC['lit_wall_T_C']:.0f} C. The paper prints no
adiabatic flame temperature in words anywhere in pp. 477-483 (the plateau exists
only as a plotted curve), so this is a consequence of the
reconstruction, not a test of it.""")
for k, v in RC.items():
    if k == "rate_max_half_unit":
        claims.exact(k, v, structural=True, note="half the last printed place of 0.003")
    else:
        claims.value(k, v, rtol=1e-5 if k != "rate_max_min_clearance" else 1e-3)

# %% [markdown]
# #### What "light-off threshold" has to mean
#
# The paper quotes two windows, and each describes a different figure:
#
# * 1-D, p. 479, about **Fig. 2**: "For inlet gas temperatures below **304 °C**, the
#   monolith temperature does not rise significantly ... between 304° and
#   **343 °C**, the reaction lights-off at some point in the channel."
# * 2-D, p. 481, about **Fig. 11**: "no reaction light-off for low inlet gas
#   temperature (case 1), $T^{o}_G <$ **293 °C**; light-off in the tube for
#   intermediate inlet gas temperatures (case 2 and 3),
#   293 $< T^{o}_G <$ **343 °C**."
#
# **Neither is a root-found threshold.** Each window is quoted from whichever of
# the four sampled cases of its figure happened to fall on the right side.
#
# **And the two are quoted for channels of different length.** "Lights off at some
# point in the channel" is a statement about a channel of some length, and the two
# figures do not share one: Fig. 2's abscissa is labelled to **0.20** in
# `X/(D RE PR)`, Fig. 11's to **0.50** (tick labels, transcribed as printed
# characters). A longer channel lights off at a lower inlet temperature, so a
# window read off the longer channel sits lower for that reason alone. The
# computation below runs each statement at the length of the figure it is made
# about, then decomposes the gap between the two windows into a length term and a
# model term.
#
# Both thresholds are root-found, and both are defined without reference to any
# figure:
#
# * **upper**, light-off *at the inlet*: the inlet temperature at which the cold
#   branch has already folded at $C_G = C_G^{o}$. Length-independent.
# * **lower**, light-off leaves the channel: the inlet temperature at which the
#   light-off position $G^{\ast}$ equals the channel length $L^{\ast}$. This one needs a
#   length, which the paper never states, so $L^{\ast}$ is swept over the lengths its
#   own 1-D figures fix (0.20 for Fig. 2, 0.25 for Figs. 3 to 5) with 0.30 as a
#   margin past both.
#
# The function below computes every reacting result of this section and of the
# Validation section from one `Channel`; a break row re-runs it with a defect.

# %%
DG, GMAX, GMAX_LONG = 1e-4, 0.40, 0.90
claims.exact("axial_step_dG", DG, structural=True, note="the axial step of every 2-D result")
L_F2, L_F3, L_F11 = pf("fig2_abscissa_end"), pf("fig3_abscissa_end"), pf("fig11_abscissa_end")
L_SWEEP = (L_F2, L_F3, 0.30)
T_F2, T_F11 = pf("lightoff_window_1D_low"), pf("lightoff_window_2D_low")
T_CASE4 = pf("lightoff_window_2D_inlet_case")
T_EQ = (310.0, T_IN, 320.0, 330.0, 340.0)
T_ROUTE = (305.0, 310.0, T_IN, 322.0, 328.0)
DG_LADDER = (4e-4, 2e-4, 1e-4, 5e-5)
CLAIM_C = pf("adequacy_claim_degC")
NU_GT = correlation(P["Nu_GrigullTratz_Qw"])        # constant wall flux, X* = 1000 G
NUF_GT = lambda G: NU_GT(max(1000.0 * G, 1e-9))     # noqa: E731
NUF_CONST = lambda G: NU_SQ                          # noqa: E731
RECON = (("baseline", 1.00, 1.00), ("dT_AD +10 %", 1.10, 1.00), ("dT_AD -10 %", 0.90, 1.00),
         ("beta x 1.5", 1.00, 1.50), ("beta / 1.5", 1.00, 1 / 1.5))
_MEMO = {}


def _fin(g):
    return g if np.isfinite(g) else 10.0


def rtag(lab):
    lab = lab.replace("+", "plus").replace("-", "minus").replace("/", "over")
    return re.sub(r"[^0-9A-Za-z]+", "_", lab).strip("_")


def reacting(ch, dG1=DG, dG2=DG, nu_r=1, nu_scale=1.0):
    """Every reacting result of the page for one channel, as ``{claim: value}``.

    ``nu_scale`` multiplies every 1-D Nusselt number (a break-table defect).
    """
    nu_sq, nu_ci = nu_scale * NU_SQ, nu_scale * NU_CI
    gt = (lambda G: nu_scale * NUF_GT(G)) if nu_scale != 1.0 else NUF_GT

    def g1q(T, Sh=nu_sq, c=ch):
        return _fin(c.lightoff_quad(T, Sh))

    def g1m(T, nu=gt, c=ch, dG=dG1):
        key = ("1d", c.kin, c.dTad, c.b0, T, nu, dG)
        if key not in _MEMO:
            _MEMO[key] = _fin(c.march_1d(T, nu, dG, GMAX))
        return _MEMO[key]

    def g2(T, Gmax=GMAX, c=ch, dG=dG2):
        key = ("2d", c.kin, c.dTad, c.b0, T, Gmax, dG, nu_r)
        if key not in _MEMO:
            _MEMO[key] = _fin(c.march_2d(T, dG, Gmax, nu_r=nu_r))
        return _MEMO[key]

    def root(f, lo=280.0, hi=380.0, xtol=1e-6):
        if not f(lo) * f(hi) < 0:          # a defect can move a threshold out of range
            lo, hi = 150.0, 700.0
        return brentq(f, lo, hi, xtol=xtol)

    m = {}
    # the two windows at their figures' own lengths, and the decomposition
    m["G_lightoff_1d_at_window_low"] = g1q(T_F2)
    m["G_lightoff_2d_at_window_low"] = g2(T_F11, GMAX_LONG)
    t1f2 = root(lambda T: g1q(T) - L_F2)
    t2f11 = root(lambda T: g2(T, GMAX_LONG) - L_F11)
    t1f11 = root(lambda T: g1q(T) - L_F11)
    m.update(T_low_1d_at_fig2_length=t1f2, T_low_2d_at_fig11_length=t2f11,
             T_low_1d_at_fig11_length=t1f11, T_windows_length_artefact_C=t2f11 - t1f2,
             T_windows_vs_paper_gap_C=abs((t2f11 - t1f2) - (T_F11 - T_F2)),
             T_windows_length_term_C=t1f11 - t1f2, T_windows_model_term_C=t2f11 - t1f11)
    # upper thresholds
    for tag, nu in (("square", nu_sq), ("circle", nu_ci)):
        m[f"T_up_1d_{tag}"] = root(lambda T: ch.fold_CG(T, ch.b0 * nu) - ch.cg0, 280.0, 420.0,
                                   xtol=1e-8)
    m["T_up_1d_shape_effect"] = m["T_up_1d_circle"] - m["T_up_1d_square"]
    first = [root(lambda T: ch.first_step_margin(T, d, nu_r=nu_r), 280.0, 520.0)
             for d in DG_LADDER]
    for d, v in zip(DG_LADDER, first):
        m[f"T_up_2d_first_step_dG{d:.0e}".replace("-0", "m")] = v
    m["T_up_2d_first_step_divergence_C"] = first[-1] - first[-2]
    m["T_up_2d_first_step_mean_per_halving_C"] = float(np.mean(np.diff(first)))
    m["G_lightoff_2d_at_case4"] = g2(T_CASE4)
    m["G_lightoff_2d_at_case4_over_fig11"] = m["G_lightoff_2d_at_case4"] / L_F11
    # lower thresholds on common lengths
    for L in L_SWEEP:
        tag = f"L{round(100 * L):03d}"
        a = root(lambda T: g1q(T) - L)
        b = root(lambda T: g1m(T) - L)
        c = root(lambda T: g2(T) - L)
        m.update({f"T_low_1d_const_{tag}": a, f"T_low_1d_GT_{tag}": b, f"T_low_2d_{tag}": c,
                  f"T_low_2d_minus_1d_const_{tag}": c - a, f"T_low_2d_minus_1d_GT_{tag}": c - b})
    m["T_low_1d_const_slope_per_005"] = ((m["T_low_1d_const_L020"] - m["T_low_1d_const_L030"])
                                         / ((L_SWEEP[-1] - L_SWEEP[0]) / 0.05))
    m["T_low_2d_L025_dG5em5"] = root(lambda T: g2(T, dG=DG / 2) - L_F3)
    # shape against entrance region and dimensionality, at L* = 0.25
    t_ci = root(lambda T: g1q(T, nu_ci) - L_F3)
    m.update(T_low_1d_circle_L025=t_ci, T_low_shape_effect_L025=t_ci - m["T_low_1d_const_L025"],
             T_low_entrance_dim_L025=m["T_low_2d_L025"] - t_ci)
    m["T_low_shape_share_pct"] = 100 * m["T_low_shape_effect_L025"] / m[
        "T_low_2d_minus_1d_const_L025"]
    # the paper's measure of adequacy: an equivalent inlet shift
    for T2 in T_EQ:
        g = g2(T2)
        tag = "at316" if T2 == T_IN else f"at{T2:.0f}"
        m[f"equiv_G2d_{tag}"] = g
        m[f"equiv_shift_const_{tag}"] = T2 - root(lambda T: g1q(T) - g, 280.0, 420.0)
        m[f"equiv_shift_GT_{tag}"] = T2 - root(lambda T: g1m(T) - g, 280.0, 420.0)
    m["equiv_shift_const_Nu_at316"] = m.pop("equiv_shift_const_at316")
    sc = [m["equiv_shift_const_Nu_at316" if T2 == T_IN else f"equiv_shift_const_at{T2:.0f}"]
          for T2 in T_EQ]
    sg = [m[f"equiv_shift_GT_at{'316' if T2 == T_IN else f'{T2:.0f}'}"] for T2 in T_EQ]
    m.update(equiv_shift_const_max=max(sc), equiv_shift_GT_max=max(sg),
             equiv_shift_const_over_claim_at316=m["equiv_shift_const_Nu_at316"] / CLAIM_C,
             equiv_shift_const_over_claim_min=min(sc) / CLAIM_C,
             equiv_shift_const_over_claim_max=max(sc) / CLAIM_C,
             equiv_shift_GT_over_claim_at316=m["equiv_shift_GT_at316"] / CLAIM_C,
             equiv_shift_GT_over_claim_min=min(sg) / CLAIM_C,
             equiv_shift_GT_over_claim_max=max(sg) / CLAIM_C,
             equiv_shift_GT_minus_claim_at316=m["equiv_shift_GT_at316"] - CLAIM_C)
    # the light-off position of the 2-D model against its axial step
    for d in DG_LADDER:
        m[f"G_lightoff_2d_316_dG{d:.0e}".replace("-0", "m")] = g2(T_IN, dG=d)
    gl = [m[f"G_lightoff_2d_316_dG{d:.0e}".replace("-0", "m")] for d in DG_LADDER]
    m["G_lightoff_2d_316_dG_order"] = observed_order(np.abs(np.diff(gl)))[-1]
    m["G_lightoff_2d_316_richardson"] = richardson(gl, order=1.0)
    m["G_lightoff_2d_316_dG1em4_vs_richardson_pct"] = 100 * (gl[2] / m[
        "G_lightoff_2d_316_richardson"] - 1)
    # two routes to the 1-D light-off position, and the C_G parametrisation
    n_giveup = 0
    for T in T_ROUTE:
        a = g1m(T, nu=NUF_CONST if nu_scale == 1.0 else (lambda G: nu_sq))
        b = g1q(T)
        with warnings.catch_warnings(record=True) as wlog:
            warnings.simplefilter("always")       # scipy's subdivision warning, counted
            c = ch.lightoff_quad_CG(T, nu_sq)
        n_giveup += any("subdivisions" in str(w.message) for w in wlog)
        m[f"route_march_minus_quad_steps_{T:.0f}"] = (a - b) / dG1
        m[f"route_quad_CG_relerr_{T:.0f}"] = abs(c - b) / b
        m[f"route_quad_{T:.0f}"] = b
    m["lightoff_march_vs_quad_max_steps"] = max(abs(m[f"route_march_minus_quad_steps_{T:.0f}"])
                                                for T in T_ROUTE)
    m["quad_CG_endpoint_error"] = max(m[f"route_quad_CG_relerr_{T:.0f}"] for T in T_ROUTE)
    m["quad_CG_rows_out_of_subdivisions"] = n_giveup
    # what the reconstruction moves
    rows = []
    for lab, f_dT, f_b in RECON:
        cv = replace(ch, dTad=ch.dTad * f_dT, b0=ch.b0 * f_b)
        tu = root(lambda T: cv.fold_CG(T, cv.b0 * nu_sq) - cv.cg0, 250.0, 500.0)
        g = g2(T_IN, c=cv)
        t1c = root(lambda T: g1q(T, c=cv) - g, 250.0, 500.0)
        t1v = root(lambda T: g1m(T, c=cv) - g, 250.0, 500.0)
        rows.append((lab, tu, g, T_IN - t1c, T_IN - t1v))
    for lab, tu, g, s1, s2 in rows:
        tag = rtag(lab)
        m.update({f"recon_{tag}_T_up": tu, f"recon_{tag}_G2d": g,
                  f"recon_{tag}_shift_const": s1, f"recon_{tag}_shift_GT": s2})
    tu_all = [r[1] for r in rows]
    s1_all, s2_all = [r[3] for r in rows], [r[4] for r in rows]
    m.update(recon_spread_upper_threshold_C=max(tu_all) - min(tu_all),
             recon_spread_equiv_shift_C=max(s1_all) - min(s1_all),
             recon_min_equiv_shift_const=min(s1_all), recon_min_equiv_shift_GT=min(s2_all),
             recon_rows_GT_shift_above_claim=int(sum(v > CLAIM_C for v in s2_all)),
             recon_rows_const_shift_above_claim=int(sum(v > CLAIM_C for v in s1_all)))
    return m


R = reacting(CH)
print(f"reacting model: {len(R)} results")
EXACT_R = ("quad_CG_rows_out_of_subdivisions", "recon_rows_GT_shift_above_claim",
           "recon_rows_const_shift_above_claim")


def register_reacting(m):
    """Register the reacting results: temperatures to 1e-5, the rest to 1e-4 (1e-3 when small)."""
    for k, v in m.items():
        if k in EXACT_R:
            continue
        if k == "G_lightoff_2d_316_dG_order":
            claims.bound(k, v, min=0.8, max=1.2, basis="theory",
                         note="backward Euler in G is first order")
        elif abs(v) > 100.0:
            claims.value(k, v, rtol=1e-5)
        elif abs(v) < 1.0 and k.endswith("_C"):
            claims.value(k, v, rtol=1e-3)
        else:
            claims.value(k, v, rtol=1e-4)
    for kind in ("GT", "const"):
        gaps = [m[f"recon_{rtag(lab)}_shift_{kind}"] - CLAIM_C for lab, _, _ in RECON]
        claims.exact(f"recon_rows_{kind}_shift_above_claim", m[f"recon_rows_{kind}_shift_above_claim"],
                     margin=float(np.min(np.abs(gaps))), min_margin=0.05,
                     note="margin in degC to the 2 degC claim")
    claims.exact("quad_CG_rows_out_of_subdivisions", m["quad_CG_rows_out_of_subdivisions"],
                 structural=True, note="scipy quad warnings, counted")


register_reacting(R)

# %%
print(f"""Fig. 2  (1-D) runs to L* = {L_F2:.2f}. The paper: below {T_F2:.0f} C 'the monolith temperature does
        not rise significantly'. This model puts that case's light-off at G* = {R['G_lightoff_1d_at_window_low']:.4f},
        {'beyond' if R['G_lightoff_1d_at_window_low'] > L_F2 else 'inside'} the figure, so the statement holds.
Fig. 11 (2-D) runs to L* = {L_F11:.2f}. The paper: 'no reaction light-off' below {T_F11:.0f} C. This model
        puts that case's light-off at G* = {R['G_lightoff_2d_at_window_low']:.4f}, {'beyond' if R['G_lightoff_2d_at_window_low'] > L_F11 else 'inside'} the figure, so that holds too.

Root-found at each figure's own length:
  1-D, Nu = {NU_SQ}, L* = {L_F2:.2f} (Fig. 2)    {R['T_low_1d_at_fig2_length']:7.2f} C
  2-D, round tube, L* = {L_F11:.2f} (Fig. 11)   {R['T_low_2d_at_fig11_length']:7.2f} C

The page's two models, on identical parameters, put the 2-D threshold {R['T_windows_length_artefact_C']:+.2f} C from
the 1-D one when each is read at its own figure's length; the paper's two windows
differ by {T_F11 - T_F2:+.0f} C, reproduced to {R['T_windows_vs_paper_gap_C']:.2f} C. Taking the same 1-D model out to Fig. 11's
length separates the two effects:

  length alone, 1-D from L* = {L_F2:.2f} to {L_F11:.2f}  {R['T_windows_length_term_C']:+8.2f} C
  model difference at the common L* = {L_F11:.2f}   {R['T_windows_model_term_C']:+8.2f} C
  ------------------------------------------------------
  the gap between the two windows           {R['T_windows_length_artefact_C']:+8.2f} C

The length term dominates and fixes the sign. Both of the paper's statements are
correct about their own figures and consistent with each other once the lengths are
accounted for; the model term has the sign of the mechanism the paper derives on
p. 482 ('lower Nusselt number means earlier light-off').""")
claims.exact("window_1d_case_beyond_fig2", bool(R["G_lightoff_1d_at_window_low"] > L_F2),
             margin=float(R["G_lightoff_1d_at_window_low"] - L_F2), min_margin=1e-3)
claims.exact("window_2d_case_beyond_fig11", bool(R["G_lightoff_2d_at_window_low"] > L_F11),
             margin=float(R["G_lightoff_2d_at_window_low"] - L_F11), min_margin=1e-3)

# %%
print(f"""Upper threshold, light-off at the inlet:
  1-D, Nu = {NU_SQ} (square)   {R['T_up_1d_square']:7.2f} C
  1-D, Nu = {NU_CI} (circle)   {R['T_up_1d_circle']:7.2f} C     shape effect {R['T_up_1d_shape_effect']:+.2f} C

The paper's 1-D bracket for this, from its running text alone, is ({T_F2:.0f}, {pf('lightoff_window_1D_high'):.0f}] C.

The 2-D model has no such threshold. A round tube has an infinite local Sherwood
number at X* = 0, so light-off always occurs at some X* > 0; what a march can report
is the inlet temperature at which light-off enters the first axial step, a number
set by the step:""")
for d in DG_LADDER:
    k = f"T_up_2d_first_step_dG{d:.0e}".replace("-0", "m")
    print(f"  dG = {d:.1e}   light-off inside the first step from {R[k]:7.2f} C"
          + ("   (the dG used everywhere else)" if d == DG else ""))
print(f"""
  {R['T_up_2d_first_step_mean_per_halving_C']:+.1f} C per halving on average, the last halving {R['T_up_2d_first_step_divergence_C']:+.2f} C: it diverges as
  dG -> 0, as the explanation requires, and is reported as the value at its step.

Fig. 11's case 4 ({T_CASE4:.0f} C, 'light-off at the tube inlet') lights off at G* = {R['G_lightoff_2d_at_case4']:.5f} here:
nonzero, {R['G_lightoff_2d_at_case4'] / L_F11:.3f} of the abscissa of its figure, so the paper's 'at the inlet' is
visual, not a threshold.""")

# %%
low = pd.DataFrame([(L, R[f"T_low_1d_const_L{round(100 * L):03d}"],
                     R[f"T_low_1d_GT_L{round(100 * L):03d}"], R[f"T_low_2d_L{round(100 * L):03d}"],
                     R[f"T_low_2d_minus_1d_const_L{round(100 * L):03d}"],
                     R[f"T_low_2d_minus_1d_GT_L{round(100 * L):03d}"]) for L in L_SWEEP],
                   columns=["L*", "1-D const Nu", "1-D Grigull-Tratz Nu(X*)", "2-D",
                            "2-D minus 1-D const", "2-D minus 1-D G-T"])
show(low.style.set_uuid("i13_lowthr").hide(axis="index").format(
    {"L*": "{:.2f}", "1-D const Nu": "{:.2f}", "1-D Grigull-Tratz Nu(X*)": "{:.2f}",
     "2-D": "{:.2f}", "2-D minus 1-D const": "{:+.2f}", "2-D minus 1-D G-T": "{:+.2f}"}))
pos = low[["2-D minus 1-D const", "2-D minus 1-D G-T"]].to_numpy()
claims.count("lower_threshold_2d_above_both_1d", pos.ravel(), 0.0, min_margin=0.5,
             relative=False, note="margin in degC")
print(f"""The 2-D model needs a higher inlet temperature than either 1-D model in
{int((pos > 0).sum())} of {pos.size} comparisons, at every channel length (the direction p. 482 requires).
The magnitude depends on L*: the 1-D threshold falls {R['T_low_1d_const_slope_per_005']:.1f} C per 0.05 of L*,
which is why the paper's two windows, read off figures whose abscissae differ by a
factor of {L_F11 / L_F2:.1f}, cannot be compared with each other as thresholds. The 2-D minus 1-D
column shrinks with length, +{R['T_low_2d_minus_1d_const_L020']:.2f} C at L* = {L_F2:.2f} and +{R['T_windows_model_term_C']:.2f} C at L* = {L_F11:.2f}.

None of these absolute temperatures is the paper's number: the {R['T_low_1d_const_L025']:.1f} C of the first
column is not evidence for the {T_F2:.0f} C it lands near. Two of the model's groups are
reconstructed (the Validation section moves this column by tens of degrees), and the
paper's {T_F2:.0f} is a sampled case. Only differences between columns, which share those
groups, are claimed.""")

# %% [markdown]
# #### Shape or dimensionality, split like-for-like on one threshold
#
# The paper's 1-D runs use the *square*-channel Nusselt number of its own Table 1,
# while its 2-D model is a *round* tube. Two different things are therefore bundled
# into "the 1-D/2-D difference", and they separate by running the 1-D model at the
# round-tube Nusselt number, on the same threshold and the same length.

# %%
print(f"""at L* = {L_F3:.2f}, all three on identical parameters:
  1-D, Nu = {NU_SQ} (square, what the paper's 1-D uses)  {R['T_low_1d_const_L025']:7.2f} C
  1-D, Nu = {NU_CI} (round tube, shape matched to 2-D)   {R['T_low_1d_circle_L025']:7.2f} C
  2-D  (round tube, entrance region resolved)        {R['T_low_2d_L025']:7.2f} C

  channel shape alone               {R['T_low_shape_effect_L025']:+.2f} C
  entrance region + dimensionality  {R['T_low_entrance_dim_L025']:+.2f} C
  -------------------------------------------
  the whole gap                     {R['T_low_2d_minus_1d_const_L025']:+.2f} C

{R['T_low_shape_share_pct']:.0f} % of the "1-D/2-D difference" here is the Nusselt number of a square channel
compared with a round one. The paper's own runs carry the same mismatch: its Figs. 2
to 6 use {NU_SQ} and its Fig. 11 is a cylinder.""")

# %%
eq = pd.DataFrame([(T2, R["equiv_G2d_at316" if T2 == T_IN else f"equiv_G2d_at{T2:.0f}"],
                    R["equiv_shift_const_Nu_at316" if T2 == T_IN else f"equiv_shift_const_at{T2:.0f}"],
                    R["equiv_shift_GT_at316" if T2 == T_IN else f"equiv_shift_GT_at{T2:.0f}"])
                   for T2 in T_EQ],
                  columns=["2-D inlet T (C)", "2-D light-off G*", "shift, const-Nu 1-D (C)",
                           "shift, Grigull-Tratz 1-D (C)"])
show(eq.style.set_uuid("i13_equiv").hide(axis="index").format(
    {"2-D inlet T (C)": "{:.0f}", "2-D light-off G*": "{:.5f}",
     "shift, const-Nu 1-D (C)": "{:+.2f}", "shift, Grigull-Tratz 1-D (C)": "{:+.2f}"}))
print(f"""The paper: the 1-D/2-D difference is 'equivalent to less than a {CLAIM_C:.0f} C difference in inlet
gas temperature', said of the Grigull-Tratz 1-D model of its Fig. 12. Measured:

                     shift at {T_IN:.0f} C   / claim at {T_IN:.0f} C   / claim over {T_EQ[0]:.0f} to {T_EQ[-1]:.0f} C
  Grigull-Tratz 1-D: {R['equiv_shift_GT_at316']:+6.2f} C       {R['equiv_shift_GT_over_claim_at316']:5.2f}             {R['equiv_shift_GT_over_claim_min']:.2f} to {R['equiv_shift_GT_over_claim_max']:.2f}
  constant-Nu 1-D:   {R['equiv_shift_const_Nu_at316']:+6.2f} C       {R['equiv_shift_const_over_claim_at316']:5.2f}             {R['equiv_shift_const_over_claim_min']:.2f} to {R['equiv_shift_const_over_claim_max']:.2f}

The Grigull-Tratz model is {R['equiv_shift_GT_minus_claim_at316']:+.2f} C over the claim at {T_IN:.0f} C; the constant-Nu model,
which produced most of the paper's figures, is {R['equiv_shift_const_over_claim_at316']:.1f} times the claim at {T_IN:.0f} C and
{R['equiv_shift_const_over_claim_min']:.1f} to {R['equiv_shift_const_over_claim_max']:.1f} times over the range. 'Adequate' is a property of the Fig. 12 procedure,
not of one-dimensionality.""")

# %% [markdown]
# #### Where the difference comes from
#
# The three sources are separable, and two of them are not about dimensionality:
#
# 1. **Channel shape.** The paper's 1-D runs use the square-channel Nusselt
#    number; the 2-D model is a round tube, whose asymptote is higher. The split
#    above measures it on the same threshold and length as the gap itself.
# 2. **The thermal entrance region.** The 2-D model starts with an infinite local
#    Nusselt number and relaxes towards the asymptote. A constant-$Nu$ 1-D model
#    never has that, which is why the paper's Fig. 12 procedure (the
#    Grigull-Tratz correlation before light-off) closes most of the gap.
# 3. **Dimensionality proper**, what is left once 1 and 2 are matched.
#
# With the entrance region put back in, the two models differ by a few degrees of
# inlet temperature, close to the paper's conclusion. What does not survive is the
# way the conclusion is stated, as a property of *one-dimensional models*, when the
# model that carries the paper's own parametric study (Figs. 2 to 6, constant
# $Nu$) is several times further off than the claim allows.

# %% [markdown]
# ## Validation
#
# ### Which X* the paper's 2-D results use, decided by its own printed number
#
# The Notation defines $X^{\ast} = 4x/(D\,Re\,Pr)$; eq. (15) requires
# $X^{\ast} = x/(D\,Re\,Pr)$. The paper prints one number that separates them: for the
# constant-wall-temperature Graetz problem with $\bar T_G^{\ast} = 1$ at the inlet and
# $T_w^{\ast} = 0$, "at the point $X^{\ast} = 0.25$, $\bar T_G^{\ast}$ is 0.0224". The same
# number decides whether the overbar is the mixing-cup or the area mean.

# %%
TBAR = pf("Tbar_at_X_0p25")


def tbar(half=0.5, X=0.30, n_r=200, r_grade=2.0, mean="cup"):
    s = Graetz(X=X, n_x=1600, n_r=n_r, half=half, x_grade=3.0, r_grade=r_grade).solve(
        ("temp", 0.0), inlet=1.0)
    return float(np.interp(0.25, s["x"], s[mean]))


TB = {"Tbar_X025_eq15": tbar(), "Tbar_X025_notation": tbar(half=2.0, X=0.40),
      "Tbar_X025_area": tbar(mean="area")}
TB["Tbar_X025_eq15_over_printed"] = TB["Tbar_X025_eq15"] / TBAR
TB["Tbar_X025_notation_over_printed"] = TB["Tbar_X025_notation"] / TBAR
TB["Tbar_X025_area_over_printed"] = TB["Tbar_X025_area"] / TBAR
print(f"""eq. (15) as printed,     mixing cup: {TB['Tbar_X025_eq15']:.6f}  ({TB['Tbar_X025_eq15_over_printed']:.2f} x the printed {TBAR})
Notation, 4x/(D Re Pr),  mixing cup: {TB['Tbar_X025_notation']:.6f}  ({TB['Tbar_X025_notation_over_printed']:.2f} x)
eq. (15), area mean:                 {TB['Tbar_X025_area']:.6f}  ({TB['Tbar_X025_area_over_printed']:.2f} x)
So the 2-D results use eq. (15)'s scaling and the overbar is the mixing-cup mean.""")
LAD_R = ((10, 1.0), (10, 2.0), (20, 1.0), (40, 1.0), (80, 2.0), (200, 2.0))
lad = [tbar(n_r=n, r_grade=g) for n, g in LAD_R]
print("radial grid, Tbar_cup(X* = 0.25)  (the paper states ten radial points):")
print("\n".join(f"  n_r = {n:3d}, grading {g:.1f}   {v:.6f}" for (n, g), v in zip(LAD_R, lad)))
TB.update({f"Tbar_X025_nr{n}_g{g:.0f}": v for (n, g), v in zip(LAD_R, lad)})
TB["Tbar_X025_converged"] = lad[-1]
TB["Tbar_X025_paper_grid_error"] = TBAR / lad[-1] - 1.0
TB["Tbar_X025_nr10_clustered_error"] = lad[1] / lad[-1] - 1.0
inside10 = min(lad[0], lad[1]) < TBAR < max(lad[0], lad[1])
print(f"""
Converged: {lad[-1]:.6f}. On ten radial points the answer lies between {min(lad[0], lad[1]):.6f} (uniform)
and {max(lad[0], lad[1]):.6f} (wall-clustered), and the printed {TBAR} {'lies' if inside10 else 'does not lie'} inside that interval. The
wall-clustered ten-point solve is {100 * TB['Tbar_X025_nr10_clustered_error']:.1f} % above converged and the printed value
{100 * TB['Tbar_X025_paper_grid_error']:.1f} % above it. So the printed number measures the radial discretisation error of
the paper's 2-D model, an error of the same order as the 1-D/2-D difference the model
is used to adjudicate. That is why every adequacy statement here is a difference
between two models on one grid.""")
for k, v in TB.items():
    claims.value(k, v, rtol=1e-6 if "error" not in k else 1e-5)
claims.exact("Tbar_printed_inside_ten_point_interval", inside10,
             margin=float(min(abs(TBAR - lad[0]), abs(TBAR - lad[1])) / TBAR), min_margin=1e-3)

# %% [markdown]
# ### Refinement: every axis that carries error

# %%
def square_orders(n_ref=384):
    """Observed order at 48 -> 96 of each axis, the other axis at the reference count."""
    ref = duct(n_ref, n_ref)
    out = {}
    for axis in (0, 1):
        eH, eT = [], []
        for n in (24, 48, 96):
            d = duct(n, n_ref) if axis == 0 else duct(n_ref, n)
            eH.append(abs(d["Nu_H1"] - ref["Nu_H1"]))
            eT.append(abs(d["Nu_T"] - ref["Nu_T"]))
        out[f"order_square_NuH1_axis{axis}"] = observed_order(eH)[1]
        out[f"order_square_NuT_axis{axis}"] = observed_order(eT)[1]
    return out


SQ_ORD = square_orders()
print("square duct, each axis refined against a (384, 384) reference, order at 48 -> 96:")
print("  " + ", ".join(f"{k} {v:.2f}" for k, v in SQ_ORD.items()))
print("The square is symmetric under swapping the axes, so the two axes must agree; the "
      "independent\npair of axes on this page is the Graetz solve below.")
for k, v in SQ_ORD.items():
    claims.bound(k, v, min=1.9, max=2.3, basis="theory",
                 note="second-order finite volumes; the reference error alone lifts 2 to about 2.07")

# %%
QW9 = wall_flux(*QW["Fig. 9"])


def pole_of(n_x, n_r, tvd=True, half=0.5, X=0.05, **kw):
    s = Graetz(X=X, n_x=n_x, n_r=n_r, half=half, r_grade=1.0).solve(("flux", QW9), tvd=tvd, **kw)
    return cross(s["x"], s["Tw"] - s["cup"])


REF_POLE = pole_of(6400, 320)
AX_NS, RAD_NS = (50, 100, 200, 400), (10, 20, 40, 80)
AX_BASE = pole_of(6400, 40)


def radial_errors(ref, wall_read="face"):
    return [abs(pole_of(3200, n, wall_read=wall_read) - ref) for n in RAD_NS]


rad_err = radial_errors(REF_POLE)


def axial_order(tvd=True):
    errs = [abs(pole_of(n, 40, tvd=tvd) - AX_BASE) for n in AX_NS]
    return fitted_order(AX_NS, errs), errs


ax_ord, ax_err = axial_order()
GR = {"nu_pole_X_fig9_refined": REF_POLE, "order_pole_radial": fitted_order(RAD_NS, rad_err),
      "order_pole_axial": ax_ord, "pole_axial_over_radial_err_coarsest": ax_err[0] / rad_err[0]}
print(f"reference pole (n_x = 6400, n_r = 320): X* = {REF_POLE:.7f}")
print(f"radial ladder n_r = {', '.join(map(str, RAD_NS))}: least-squares order "
      f"{GR['order_pole_radial']:.2f}")
print(f"axial ladder n_x = {', '.join(map(str, AX_NS))}: least-squares order "
      f"{GR['order_pole_axial']:.2f}")
print(f"""The axial ladder is noisy and its slope comes out above 2. That is not third-order
accuracy: the pole is located on this axis, so the locator's own interpolation error
is inside the ladder, and the axial errors are {GR['pole_axial_over_radial_err_coarsest']:.2f} times the radial ones at the
coarsest rung, so it is not the limiting axis. With the TVD correction off the same
fit gives first order (break table).""")


def paper_scheme(tvd=False):
    """The paper's stated counts: ten radial points and its axial step."""
    return pole_of(int(round(0.05 / pf("axial_step_2D"))), int(pf("radial_points_2D")), tvd=tvd)


p_pg = paper_scheme()
lo_pg, hi_pg = min(p_pg, pole_10r), max(p_pg, pole_10r)
GR.update(nu_pole_X_fig9_paper_scheme=p_pg,
          pole_paper_scheme_dev_pct=100 * (p_pg / REF_POLE - 1),
          pole_paper_counts_short_min_pct=100 * (1 - hi_pg / X_SIGN[0]),
          pole_paper_counts_short_max_pct=100 * (1 - lo_pg / X_SIGN[0]),
          pole_paper_counts_max_dev_pct=100 * max(abs(p_pg / REF_POLE - 1),
                                                  abs(pole_10r / REF_POLE - 1)),
          pole_to_0p032_pct=100 * (X_SIGN[0] / REF_POLE - 1))
GR["pole_reach_over_grid_error"] = GR["pole_to_0p032_pct"] / GR["pole_paper_counts_max_dev_pct"]
print(f"""
The paper's stated discretisation is ten radial points and dX* = {pf('axial_step_2D')}; with the TVD
correction off as well (a first-order marching scheme) the pole is at X* = {p_pg:.6f},
{GR['pole_paper_scheme_dev_pct']:+.1f} % from converged and {100 * (1 - p_pg / X_SIGN[0]):.1f} % short of {X_SIGN[0]}. That is this page's radial operator on
the paper's counts, not the paper's code, so the conclusion is a margin: the two
paper-count solves lie {GR['pole_paper_counts_short_min_pct']:.1f} to {GR['pole_paper_counts_short_max_pct']:.1f} % short of {X_SIGN[0]} and within {GR['pole_paper_counts_max_dev_pct']:.1f} % of converged,
while reaching {X_SIGN[0]} takes {GR['pole_to_0p032_pct']:+.1f} %, {GR['pole_reach_over_grid_error']:.0f} times that discretisation error.""")
for k, v in GR.items():
    if k == "order_pole_radial":
        claims.bound(k, v, min=1.9, max=2.3, basis="theory",
                     note="second order radially")
    elif k == "order_pole_axial":
        claims.bound(k, v, min=1.6, max=3.5, basis="theory",
                     note="van Leer, located on this axis, so noisy above 2")
    else:
        claims.value(k, v, rtol=1e-6 if "pct" not in k and "over" not in k else 1e-4)

# %% [markdown]
# ### Two solvers for the same discrete Graetz equations
#
# The production route is deferred correction (the van Leer correction in the
# residual, one LU of the upwind operator, under-relaxed). An independent route
# solves the same equations by Newton with the exact van Leer Jacobian from
# `ConvFluxTVD`. Each Newton iteration factorises the whole Jacobian, so on the
# production grid it costs several times the deferred
# correction's single factorisation, which is why the comparison runs on a
# coarser grid.

# %%
GV = Graetz(X=0.05, n_x=1000, n_r=40, r_grade=1.0)
s_dc, s_nw = GV.solve(("flux", QW9)), GV.solve_newton(("flux", QW9))
s_dc1 = GV.solve(("flux", QW9), omega=1.0, accept=np.inf)
NW = {"graetz_newton_vs_deferred_field": float(np.max(np.abs(s_dc["T"] - s_nw["T"]))),
      "graetz_newton_vs_deferred_pole": abs(cross(s_dc["x"], s_dc["Tw"] - s_dc["cup"])
                                            - cross(s_nw["x"], s_nw["Tw"] - s_nw["cup"]))}
print("n_x = 1000, n_r = 40, Fig. 9's flux:")
print(f"  deferred correction, omega = 0.6: {s_dc['nit']} sweeps")
print(f"  Newton with ConvFluxTVD, damping='monotone': {s_nw['nit']} iterations")
print(f"  largest field difference (here {NW['graetz_newton_vs_deferred_field']:.1e}), pole difference "
      f"(here {NW['graetz_newton_vs_deferred_pole']:.1e}), both below 1e-9")
print(f"  without under-relaxation (omega = 1) the correction still changes by "
      f"{s_dc1['incr']:.1e} after {s_dc1['nit']} sweeps")
claims.bound("graetz_newton_vs_deferred_field", NW["graetz_newton_vs_deferred_field"],
             max=1e-9, basis="tolerance", note="both solvers stop at 1e-10")
claims.bound("graetz_newton_vs_deferred_pole", NW["graetz_newton_vs_deferred_pole"],
             max=1e-9, basis="tolerance", note="both solvers stop at 1e-10")
claims.exact("graetz_deferred_sweeps", int(s_dc["nit"]))
claims.exact("graetz_newton_iterations", int(s_nw["nit"]))
claims.value("graetz_omega1_last_change", s_dc1["incr"], rtol=1e-3)
claims.exact("graetz_omega1_sweeps", int(s_dc1["nit"]))

# %% [markdown]
# ### The Grigull-Tratz correlation against the solve, and a third X*
#
# The paper reprints two correlations from Grigull and Tratz (1965) and uses the
# constant-flux one as the Nusselt number of its Fig. 12 one-dimensional model.
# They are printed without a definition of the $X^{\ast}$ they use, and it is neither
# of the two already in play. The correlation is read from the CSV as printed.

# %%
s_gt = Graetz(X=0.06, n_x=4000, n_r=200, r_grade=2.0, x_grade=3.0).solve(
    ("flux", lambda X: -np.ones_like(X)))                     # constant wall flux
NU_NUM = local_nusselt(s_gt)       # q < 0 and (Tw - Tbar) < 0, so positive
G_GT = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05)


def gt_ratio(scale):
    return [NU_GT(scale * g) / float(np.interp(g, s_gt["x"], NU_NUM)) for g in G_GT]


def gt_dev(scale=1000.0):
    return float(max(abs(r - 1.0) for r in gt_ratio(scale)))


def gt_factor(scale):
    return float(min(gt_ratio(scale)))


GTC = {"grigull_tratz_max_reldev": gt_dev(), "grigull_tratz_X_eq_G_min_factor": gt_factor(1.0),
       "grigull_tratz_X_eq_4G_min_factor": gt_factor(4.0)}
print(f"""Over G = {G_GT[0]} to {G_GT[-1]} the correlation reproduces this page's solve to {100 * GTC['grigull_tratz_max_reldev']:.2f} % with
X* = 1000 x/(D Re Pr), and overpredicts by at least a factor {GTC['grigull_tratz_X_eq_4G_min_factor']:.1f} (X* = 4x/(D Re Pr)) or
{GTC['grigull_tratz_X_eq_G_min_factor']:.1f} (X* = x/(D Re Pr)) under the paper's own two meanings. So the paper carries three
definitions of one symbol: 4x/(D Re Pr) in the Notation and eq. (9), x/(D Re Pr) in
eqs. (15)-(16) and on every figure abscissa, and 1000 x/(D Re Pr) in the reprinted
correlations. The page uses each equation with the scaling it requires.""")
for k, v in GTC.items():
    claims.value(k, v, rtol=1e-6)

# %% [markdown]
# ### Light-off positions: two routes in 1-D, and the axial step in 2-D
#
# The constant-$Nu$ 1-D light-off position is computed two ways that share no
# discretisation and no parametrisation: an explicit march of $dC_G/dG$, with the
# position root-found inside the step that crosses the fold, and a Gauss-Legendre
# quadrature of $dG/dC_w$ along the cold branch. A third route, the same integral
# parametrised by $C_G$ with scipy's adaptive `quad`, has a square-root endpoint,
# because $dC_G/dC_w$ vanishes at the fold. For the 2-D model the light-off
# position is followed over four axial steps.

# %%
rt = pd.DataFrame([(T, R[f"route_quad_{T:.0f}"], R[f"route_march_minus_quad_steps_{T:.0f}"],
                    R[f"route_quad_CG_relerr_{T:.0f}"]) for T in T_ROUTE],
                  columns=["T_G0 (C)", "G* quadrature in C_w", "march minus quadrature (steps)",
                           "quad in C_G, rel. error"])
show(rt.style.set_uuid("i13_routes").hide(axis="index").format(
    {"T_G0 (C)": "{:.0f}", "G* quadrature in C_w": "{:.6f}",
     "march minus quadrature (steps)": "{:+.3f}", "quad in C_G, rel. error": "{:.2e}"}))
print(f"""The march (dG = {DG:.0e}) and the quadrature agree to within {R['lightoff_march_vs_quad_max_steps']:.2f} axial steps: the
explicit march's own first-order error, which the dG row of the break table scales.
The C_G parametrisation is up to {100 * R['quad_CG_endpoint_error']:.2f} % out; scipy's quad ran out of its 600
subdivisions on {R['quad_CG_rows_out_of_subdivisions']} of the {len(T_ROUTE)} rows. Had that form been the reference, an error of that
size would have sat under every threshold on this page.""")
gl = [R[f"G_lightoff_2d_316_dG{d:.0e}".replace("-0", "m")] for d in DG_LADDER]
print(f"""
2-D light-off position at {T_IN:.0f} C, root-found in the crossing step:
""" + "\n".join(f"  dG = {d:.1e}   G* = {g:.6f}" for d, g in zip(DG_LADDER, gl)) + f"""
observed order {R['G_lightoff_2d_316_dG_order']:.2f} (backward Euler: 1); first-order Richardson on the last two
steps {R['G_lightoff_2d_316_richardson']:.6f}, so the dG = {DG:.0e} value is {R['G_lightoff_2d_316_dG1em4_vs_richardson_pct']:+.3f} % from it. The lower 2-D
threshold at L* = {L_F3:.2f} moves from {R['T_low_2d_L025']:.3f} to {R['T_low_2d_L025_dG5em5']:.3f} C when the step is halved.
Every 2-D result on this page is the value at dG = {DG:.0e}.""")

# %% [markdown]
# ### What the reconstruction moves
#
# Two groups were reconstructed. The thresholds depend on them strongly; the
# difference between the two models, which is what the page concludes about, much
# less. The table re-runs the comparison at the paper's Fig. 1 inlet temperature
# with each group varied.

# %%
sens = pd.DataFrame([(lab,) + tuple(R[f"recon_{rtag(lab)}_{c}"]
                                    for c in ("T_up", "G2d", "shift_const", "shift_GT"))
                     for lab, _, _ in RECON],
                    columns=["reconstruction", "1-D upper threshold (C)",
                             f"2-D light-off G* at {T_IN:.0f} C", "shift, const-Nu 1-D (C)",
                             "shift, G-T 1-D (C)"])
show(sens.style.set_uuid("i13_sens").hide(axis="index").format(
    {"1-D upper threshold (C)": "{:.2f}", f"2-D light-off G* at {T_IN:.0f} C": "{:.5f}",
     "shift, const-Nu 1-D (C)": "{:+.2f}", "shift, G-T 1-D (C)": "{:+.2f}"}))
print(f"""A +/-10 % change in dT_AD and a factor 1.5 either way in beta move the absolute upper
threshold over {R['recon_spread_upper_threshold_C']:.1f} C and the constant-Nu equivalent shift over {R['recon_spread_equiv_shift_C']:.2f} C.

Against the paper's {CLAIM_C:.0f} C: the constant-Nu 1-D shift exceeds it in {R['recon_rows_const_shift_above_claim']} of {len(RECON)} variations
(smallest {R['recon_min_equiv_shift_const']:+.2f} C); the Grigull-Tratz 1-D shift exceeds it in {R['recon_rows_GT_shift_above_claim']} of {len(RECON)}
(smallest {R['recon_min_equiv_shift_GT']:+.2f} C). So the statement that the constant-Nu model misses the claim does
not depend on the reconstruction, while whether the Grigull-Tratz model meets it does.

What no perturbation of this kind can test is whether the model structure is the
paper's: constant oxygen, Le = 1 and Sh = Nu are the paper's assumptions, and the
first is a break-table row.""")

# %% [markdown]
# ### The boundary read
#
# Both tube models take their wall value from `compute_boundary_values`, not
# from the last cell centre. For a Neumann wall the two differ at first order in the
# near-wall spacing. On the Graetz problem the wall value sits in the denominator of
# the Nusselt number, so the pole moves when it is read from the last cell; the
# cell below re-runs the pole that way on the production grid and on the paper's
# ten radial points, and compares the two reads on a test profile of the reacting
# channel's radial grid.

# %%
def cell_read_shifts(r_grade=1.0):
    out = {}
    for tag, nx, nr in (("prod", GZ_PROD["n_x"], GZ_PROD["n_r"]), ("nr10", G10["n_x"], N_R10)):
        g = Graetz(X=0.05, n_x=nx, n_r=nr, r_grade=r_grade)
        sf, sc = g.solve(("flux", QW9)), g.solve(("flux", QW9), wall_read="cell")
        out[f"pole_cell_read_shift_{tag}"] = (cross(sc["x"], sc["Tw"] - sc["cup"])
                                             - cross(sf["x"], sf["Tw"] - sf["cup"]))
    return out


N_R10 = int(pf("radial_points_2D"))
CS = cell_read_shifts()
print(f"pole read from the last cell instead of the face: {CS['pole_cell_read_shift_prod']:+.2e} on "
      f"{GZ_PROD['n_x']} x {GZ_PROD['n_r']} cells, {CS['pole_cell_read_shift_nr10']:+.2e} on the paper's {N_R10} radial points")
for k, v in CS.items():
    claims.value(k, v, rtol=1e-4)


def boundary_reads(grade=2.0):
    rad = Radial(48, grade=grade)
    c_test = np.exp(-3.0 * (1.0 - rad.r_c)).reshape(-1, 1)      # a smooth test profile
    return {f"boundary_read_vs_lastcell{'' if q else '_q0'}":
            abs(rad.wall_value(c_test, q) - c_test[-1, 0]) for q in (-0.7, 0.0)}


BR = boundary_reads()
sd = Graetz(X=0.05, n_x=200, n_r=40).solve(("temp", 0.37), inlet=1.0)
dirich = float(np.max(abs(sd["Tw"] - 0.37)))
print(f"test profile exp(-3(1-R)) on the reacting channel's 48 graded cells, face against last cell: "
      f"{BR['boundary_read_vs_lastcell_q0']:.2e} at q = 0, "
      f"{BR['boundary_read_vs_lastcell']:.2e} at q = -0.7")
print(f"Dirichlet wall at 0.37: largest deviation of the read value over the tube "
      f"(here {dirich:.1e}), below 1e-14")
for k, v in BR.items():
    claims.value(k, v, rtol=1e-6)
claims.bound("boundary_read_dirichlet_exact", dirich, max=1e-14, basis="identity",
             structural=True, note="the boundary reconstruction returns a Dirichlet value")

# %% [markdown]
# ### The defect-injection table
#
# Each row breaks one thing on purpose, **re-runs** the computation behind its
# claims with the defect in place, and hands the result to `claims.breaks`, which
# records whether the claim's own check catches it. The cell after the table prints
# the uncaught pairs.

# %%
BREAK = []


def brk(row, values):
    for name, v in values.items():
        BREAK.append((row, name, claims.breaks(row, name, v)))


# Table 1
PRINTED = ("nu_h1_square_printed_relerr", "nu_t_square_printed_relerr", "nu_t_circle_printed_relerr",
           "table1_cells_within_printed_digits")
du_p = duct_ladders(plug=True)
brk("plug velocity instead of Poiseuille", table1_claims(du_p))
brk("plug velocity instead of Poiseuille",
    {"nu_h1_square_two_routes_reldiff": abs(SPEC[-1][0] - du_p["square_Nu_H1"]) / du_p["square_Nu_H1"]})
brk("unweighted eigenproblem (drop u)", table1_claims(duct_ladders(weight_eigen=False)))
brk("nu=0 (Cartesian) in the radial div", table1_claims(duct_ladders(nu_r=0)))
row = "Table 1 columns swapped (Nu_T read as Nu_H)"
brk(row, {k: v for k, v in table1_claims(DU, col_h="Nu_T", col_t="Nu_H").items() if k in PRINTED})
brk(row, {"nu_h1_circle_printed_relerr": printed_relerr(float(nu_c_exact), t1("circular", "Nu_T")),
          "nu_h1_triangle_printed_relerr": printed_relerr(float(nu_t_exact),
                                                          t1("equilateral triangle", "Nu_T"))})
sp_w = spectral_ladder(weight_eigen=False)
brk("unweighted eigenproblem (drop u)",
    {f"nu_t_square_spectral_N{N}": b for N, (_, b) in zip(SPEC_N, sp_w)}
    | {"nu_t_square_spectral_last_increment": sp_w[-1][1] - sp_w[-2][1],
       "nu_t_square_spectral_above_interval": sp_w[-1][1] - (t1("square", "Nu_T") + 5e-4)})

# %%
N_R10 = int(pf("radial_points_2D"))
row = "the paper's own 10 radial points"
AN10 = {nm: anomaly(coef, dict(GZ_PROD, n_r=N_R10)) for nm, coef in QW.items()}
brk(row, {k: v for k, v in anomaly_claims(AN10, pole_10r).items() if "n_r10" not in k})
ref10 = pole_of(6400, N_R10)
brk(row, {"nu_pole_X_fig9_refined": ref10, "pole_to_0p032_pct": 100 * (X_SIGN[0] / ref10 - 1)})
for nm in QW:
    a, b, p = QW[nm]
    tag = nm.replace("Fig. ", "fig")
    row = f"{nm} flux read with exponent 6 for 8"
    an6 = dict(AN, **{nm: anomaly((a, b, 6.0))})
    p10_6 = anomaly((a, b, 6.0), G10)["pole"] if nm == "Fig. 9" else pole_10r
    other = "fig8" if tag == "fig9" else "fig9"
    brk(row, {k: v for k, v in anomaly_claims(an6, p10_6).items()
              if other not in k and (tag == "fig9" or "fig8" in k)})
row = "half = 2 (the Notation's X*)"
AN2 = {nm: anomaly(coef, dict(X=0.20, n_x=2500, n_r=60, r_grade=1.0, half=2.0))
       for nm, coef in QW.items()}
brk(row, {k: v for k, v in anomaly_claims(AN2, pole_10r).items() if "n_r10" not in k})
s_neg = Graetz(**GZ_PROD).solve(("flux", lambda X: -QW9(X)))
brk("wall flux with the sign of eq. (17) as printed", sign_claims(sign_rows(s_neg)))
tb2 = {"cup": tbar(half=2.0), "area": tbar(half=2.0, mean="area")}
brk(row, {"Tbar_X025_eq15": tb2["cup"], "Tbar_X025_eq15_over_printed": tb2["cup"] / TBAR,
          "Tbar_X025_area": tb2["area"], "Tbar_X025_area_over_printed": tb2["area"] / TBAR})
tbn = tbar(X=0.40)
brk("eq. (15)'s coefficient for the Notation's run",
    {"Tbar_X025_notation": tbn, "Tbar_X025_notation_over_printed": tbn / TBAR})
row = "the paper's own 10 radial points"
lad10 = [tbar(n_r=N_R10, r_grade=g) for _, g in LAD_R]
brk(row, {f"Tbar_X025_nr{n}_g{g:.0f}": v for (n, g), v in zip(LAD_R, lad10) if n != N_R10}
    | {"Tbar_X025_converged": lad10[-1], "Tbar_X025_paper_grid_error": TBAR / lad10[-1] - 1.0,
       "Tbar_X025_nr10_clustered_error": lad10[1] / lad10[-1] - 1.0})
v_g = tbar(n_r=N_R10, r_grade=2.0)
brk("ten-point grid graded for the uniform run",
    {"Tbar_X025_nr10_g1": v_g,
     "Tbar_printed_inside_ten_point_interval": bool(min(v_g, lad[1]) < TBAR < max(v_g, lad[1]))})
p10g = anomaly(QW["Fig. 9"], dict(G10, r_grade=2.0))["pole"]
lo_g, hi_g = min(p_pg, p10g), max(p_pg, p10g)
dev_g = 100 * max(abs(p_pg / REF_POLE - 1), abs(p10g / REF_POLE - 1))
brk("ten-point grid graded for the uniform run",
    {k: v for k, v in anomaly_claims(AN, p10g).items() if "n_r10" in k}
    | {"pole_paper_counts_short_min_pct": 100 * (1 - hi_g / X_SIGN[0]),
       "pole_paper_counts_short_max_pct": 100 * (1 - lo_g / X_SIGN[0]),
       "pole_paper_counts_max_dev_pct": dev_g, "pole_reach_over_grid_error": GR["pole_to_0p032_pct"] / dev_g})
v_u = tbar(n_r=N_R10, r_grade=1.0)
brk("ten-point grid not graded",
    {"Tbar_X025_nr10_g2": v_u, "Tbar_X025_nr10_clustered_error": v_u / lad[-1] - 1.0,
     "Tbar_printed_inside_ten_point_interval": bool(min(lad[0], v_u) < TBAR < max(lad[0], v_u))})
ax_off, ax_err_off = axial_order(tvd=False)
brk("TVD correction off", {"order_pole_axial": ax_off,
                           "pole_axial_over_radial_err_coarsest": ax_err_off[0] / rad_err[0]})
p_pg_t = paper_scheme(tvd=True)
lo_t, hi_t = min(p_pg_t, pole_10r), max(p_pg_t, pole_10r)
dev_t = 100 * max(abs(p_pg_t / REF_POLE - 1), abs(pole_10r / REF_POLE - 1))
brk("paper's counts with the TVD correction on",
    {"nu_pole_X_fig9_paper_scheme": p_pg_t,
     "pole_paper_scheme_dev_pct": 100 * (p_pg_t / REF_POLE - 1),
     "pole_paper_counts_short_min_pct": 100 * (1 - hi_t / X_SIGN[0]),
     "pole_paper_counts_short_max_pct": 100 * (1 - lo_t / X_SIGN[0]),
     "pole_paper_counts_max_dev_pct": dev_t,
     "pole_reach_over_grid_error": GR["pole_to_0p032_pct"] / dev_t})
s_dc3 = GV.solve(("flux", QW9), maxit=3, accept=np.inf)
brk("deferred correction stopped after 3 sweeps",
    {"graetz_newton_vs_deferred_field": float(np.max(np.abs(s_dc3["T"] - s_nw["T"]))),
     "graetz_newton_vs_deferred_pole": abs(cross(s_dc3["x"], s_dc3["Tw"] - s_dc3["cup"])
                                           - cross(s_nw["x"], s_nw["Tw"] - s_nw["cup"])),
     "graetz_deferred_sweeps": int(s_dc3["nit"])})
brk("Newton started from the inlet value instead of the upwind solution",
    {"graetz_newton_iterations": int(GV.solve_newton(("flux", QW9), guess="inlet")["nit"])})
s_dc1b = GV.solve(("flux", QW9), omega=0.8, accept=np.inf)
brk("omega = 0.8 for the unrelaxed run", {"graetz_omega1_last_change": s_dc1b["incr"],
                                         "graetz_omega1_sweeps": int(s_dc1b["nit"])})
brk("X* = G in Grigull-Tratz", {"grigull_tratz_max_reldev": gt_dev(1.0)})
brk("X* = 1000 G for the paper's own two meanings",
    {"grigull_tratz_X_eq_G_min_factor": gt_factor(1000.0),
     "grigull_tratz_X_eq_4G_min_factor": gt_factor(1000.0)})
brk("reacting channel's radial grid not graded", boundary_reads(grade=1.0))
brk("wall-graded Graetz radial grids", cell_read_shifts(r_grade=2.0))
brk("square reference grid one rung too coarse (192)", square_orders(192))
brk("radial reference grid one rung too coarse (160)",
    {"order_pole_radial": fitted_order(RAD_NS, radial_errors(pole_of(6400, 160)))})
brk("wall value read from the last cell", {"order_pole_radial": fitted_order(
    RAD_NS, radial_errors(pole_of(6400, 320, wall_read="cell"), wall_read="cell"))})

# %%
# the reconstruction and the reacting model
groups = (pf("delta_H") * 1e6, CG0, pf("channel_diameter") * 1e-3, NU_SQ, T_IN + 273.0)
dT_m, _, b0_m = reconstruct_groups(*groups, molar=False)
_, beta_n, _ = reconstruct_groups(*groups, n_T=0.0)
row = "C_PG read as molar (M_G dropped)"
ch_m = replace(CH, dTad=dT_m)
g_m1, g_m2 = ch_m.lightoff_quad(T_F2, NU_SQ), ch_m.march_2d(T_F11, DG, GMAX_LONG)
brk(row, {"dTad_reconstructed": dT_m, "G_lightoff_1d_at_window_low": g_m1,
          "G_lightoff_2d_at_window_low": g_m2, "window_1d_case_beyond_fig2": bool(g_m1 > L_F2),
          "window_2d_case_beyond_fig11": bool(g_m2 > L_F11), "lit_wall_T_C": T_IN + dT_m})
brk("diffusivity left at 300 K", {"beta_reconstructed": beta_n})
brk("E_a sign flipped (+961)", rate_claims(replace(CH, kin=replace(KIN, ear=-KIN.ear)), scan=False))
CH_O2 = replace(CH, kin=replace(KIN, o2="depleted"))
brk("oxygen depleted stoichiometrically", rate_claims(CH_O2))


def low_count(m):
    return int(sum(m[f"T_low_2d_minus_1d_{k}_L{round(100 * L):03d}"] > 0
                   for L in L_SWEEP for k in ("const", "GT")))


for row, kw in (("oxygen depleted stoichiometrically", dict(ch=CH_O2)),
                ("nu=0 (Cartesian) in the 2-D radial div", dict(ch=CH, nu_r=0)),
                ("dG x 10 (1e-3) in the 1-D march", dict(ch=CH, dG1=1e-3)),
                ("1-D Nusselt number on the radius (x 2)", dict(ch=CH, nu_scale=2.0))):
    Rb = reacting(**kw)
    brk(row, Rb | {"lower_threshold_2d_above_both_1d": low_count(Rb)})
brk("quad subdivision limit 600 -> 30",
    {"quad_CG_endpoint_error": float(max(abs(CH.lightoff_quad_CG(T, NU_SQ, limit=30)
                                             / CH.lightoff_quad(T, NU_SQ) - 1.0)
                                         for T in T_ROUTE))})

tbl = pd.DataFrame(BREAK, columns=["row", "claim", "caught"])
summary = tbl.groupby("row", sort=False).agg(pairs=("claim", "size"), caught=("caught", "sum"))
print(summary.to_string())
print(f"{len(summary)} rows, {len(tbl)} (row, claim) pairs, {int(tbl.caught.sum())} caught")

# %% [markdown]
# A row leaves a pair inside the claim's check wherever its defect does not enter
# the quantity, and the rows that re-run a whole computation report every claim of
# it. The four rows that re-run the reacting model: the 2-D slab leaves every 1-D
# and quadrature result where it was; the 1-D step leaves every 2-D and quadrature
# result; the doubled 1-D Nusselt number leaves every 2-D result; oxygen depletion
# leaves $1/k_a$, the lit-wall temperature (which contains only $\Delta T_{AD}$)
# and the counts, whose margins it does not cross. On Table 1 the unweighted
# eigenproblem leaves the constant-flux numbers and fRe, which contain no
# eigenvalue, and the circle's geometry factor does not reach the square. The rows
# that change the Graetz grid or its axial coefficient leave the flux zeros, which
# are roots of the printed polynomials, and the Nusselt zeros, which equal them by
# construction up to interpolation. The paper's counts
# with the TVD correction on leave the claims set by the uniform ten-point solve,
# which is the larger of the two paper-count poles. The flipped activation energy
# leaves the lit-wall temperature, and the structural
# half unit of the printed rate maximum and the scipy subdivision count move under
# no row. The cell below lists every uncaught pair and the claims, if any, that no
# row catches.

# %%
unc = tbl[~tbl.caught]
print(f"{int(tbl.caught.sum())} of {len(tbl)} pairs caught; {len(unc)} uncaught, by row:")
print("\n".join(f"  {r} ({int((unc.row == r).sum())}): " + ", ".join(unc.claim[unc.row == r])
                for r in dict.fromkeys(unc.row)))
unbroken = [n for n, c in claims.metrics.items() if not c.get("structural")
            and not any(caught for r, n2, caught in BREAK if n2 == n)]
print(f"claims without a catching row: {len(unbroken)}" + (": " + ", ".join(unbroken) if unbroken else ""))
n_struct = sum(1 for c in claims.metrics.values() if c.get("structural"))
print(f"{len(claims.metrics)} claims, {n_struct} structural")

# %% [markdown]
# ## What pymrm adds
#
# **What is pymrm's here.** Every field on this page is discretised with pymrm
# operators assembled once: `construct_grad` and `construct_div` for the duct
# cross-sections (Cartesian and cylindrical through `nu`) and for the radial
# diffusion of both tube models, `construct_convflux_upwind` and
# `interp_cntr_to_stagg_tvd` for the van Leer convection along the tube, with
# `{"outflow": True}` as the outlet a convection-only direction needs. The wall
# flux that varies along the tube is an array-valued boundary `d`. In the reacting
# channel the wall flux enters through `shapes_d`, whose boundary vector per unit
# flux equals the plain boundary vector of the dictionary with `d = 1`; what makes
# each light-off step one linear solve is that the interior is linear in that flux,
# not the API. `compute_boundary_values` supplies every wall value, which is what the kinetics
# are evaluated at. `ConvFluxTVD` with `newton(damping="monotone")` provides an
# independent solver for the same discrete Graetz equations, which converges in a
# handful of iterations and agrees with the deferred correction to the solvers'
# tolerance.
#
# **What is not pymrm's.** The closed forms (sympy), the spectral square-duct route
# (numpy and scipy), every root, the quadrature of the light-off integral and the
# analysis of the paper are this page's own work on top of scipy; they would read
# the same with any discretisation library. The physics is entirely the paper's.
#
# **What the page establishes about the paper:**
#
# * **Table 1 stops being a citation.** Five of its eight cells are derived, two
#   in closed form, and one printed cell (the square duct at constant wall
#   temperature) disagrees with two converged routes by more than its printed
#   digits allow; the page does not attribute it, because the cited 1971 source was
#   not consulted.
# * **The Nusselt singularity is root-found and its cause corrected.** At the flux
#   zero the Nusselt number passes through zero; the pole is the zero of
#   $T_w - \bar T_G$, downstream of it in both printed flux functions. The paper's
#   next paragraph gives the correct mechanism, so this is a wrong summary rather
#   than a missed effect.
# * **The 1-D/2-D thresholds are root-found instead of sampled, and on a common
#   length.** The paper's two windows are quoted at different channel lengths; run
#   at each figure's own length, the models reproduce both, and the gap between
#   them decomposes into a dominant length term and a model term of the sign the
#   paper's mechanism requires. About half of the remaining 1-D/2-D difference on a
#   common length is channel shape, not dimensionality.
# * **"Adequate" is given a number and a scope.** The "less than 2 °C" claim is
#   measured as an equivalent inlet shift for both 1-D models, at the inlet
#   temperature of Fig. 1 and over the range of the paper's cases, and printed as
#   a ratio to the claim in the Results.
# * **Three meanings of $X^{\ast}$ are separated**, each against a computation, and
#   the printed check value that settles one of them measures the paper's own
#   radial grid.
#
# ### What this page cannot conclude
#
# * **Nothing about any figure.** Figures are used only for characters printed
#   inside their frames, each listed and counted in *The data*. No profile shape,
#   no light-off position against the paper's own, and none of the conversion
#   curves of Figs. 4 to 6 and 12 is checked.
# * **The channel length is not printed.** The lower threshold needs one, so it is
#   swept and every conclusion is checked over the sweep; no single "the" lower
#   threshold is claimed.
# * **There is no 2-D inlet light-off threshold**; the number reported for it is
#   the value at its axial step and diverges as the step shrinks.
# * **Two dimensional groups are reconstructions**, so no absolute threshold here
#   should be read as the paper's number; only differences between models on the
#   same parameters are claimed.
# * **The square-duct $Nu_T$ discrepancy is not attributed**, and three of Table
#   1's cells are not reproduced (one cannot be).
# * **The transient model (eqs. 1 to 4, Fig. 7) is not built.** Its only
#   quantitative statement, the time by which the transient is almost complete,
#   needs the solid heat capacity and the monolith dimensions, none of which is
#   printed. `I1.2` does transient monolith light-off against a printed
#   table.

# %% [markdown]
# ## Reuse
#
# `Graetz` (in `model.py`) is a general developing-flow solver for a round tube:
# pick the wall condition with `("flux", f)` or `("temp", value)` and the axial
# coefficient with `half`. The wall flux may be any function of the axial
# coordinate, because pymrm takes an array-valued `d` in the boundary dictionary.
# On a large grid keep the deferred correction **with its under-relaxation**
# ($\omega = 0.6$; without it the correction does not settle, as printed above),
# and use `solve_newton` (`ConvFluxTVD`, `newton(damping="monotone")`) where the
# grid is small enough to refactorise every iteration or the profile is not as
# smooth as this one (pitfall P10).
#
# `Channel.march_2d` is the piece worth copying for a reacting wall. The interior
# is linear in the wall flux, so the response to a unit flux is computed once per
# step length, and each implicit step needs one linear solve and a **scalar** root
# for the wall value. A 2-D nonlinear boundary problem becomes the
# same algebra as Fig. 1's mass-transfer line with an effective $\beta$, the fold
# is explicit, and light-off is root-found inside its step. It needs $Le = 1$ (one
# field); with $Le \ne 1$ temperature and concentration are two fields with the
# field axis last, and the step becomes a two-unknown wall problem.
#
# Read wall values with `compute_boundary_values`, never from the last cell centre:
# read from the last cell, the Graetz pole moves as printed in the boundary-read
# check, by far more on a ten-point radial grid than on the production grid, and
# the radial convergence order drops to about one (break table).
#
# Related pages: [`A3.15`](../A3.15-graetz-nusselt/) for the Graetz eigenvalue
# route to the same circular-duct limit, [`I1.2`](../I1.2-catalytic-converter-lightoff/)
# for transient monolith light-off, [`D2.2`](../D2.2-van-welsenaere-froment-runaway/)
# for a threshold located as a root on a single-valued branch.
#
# ## Cite the source, not this page
#
# Heck, R. H., Wei, J. & Katzer, J. R., *Mathematical Modeling of Monolithic
# Catalysts*, AIChE Journal **22**(3) 477-484 (1976), doi:10.1002/aic.690220310.
#
# Table 1 is attributed by the paper to Shah and London (1971), Stanford University
# Dept. of Mechanical Engineering Tech. Rept. No. 75, **not consulted here**. The
# kinetics are Voltz, Morgan, Liederman & Jacob, *Ind. Eng. Chem. Prod. Res.
# Develop.* **12**, 295 (1973), **not consulted** either; the four constants are
# used exactly as Heck et al. print them, the reprint route of `AGENTS.md`. The
# two Nusselt correlations are Grigull & Tratz, *Intern. J. Heat Mass Transfer*
# **8**, 669 (1965), not consulted; they are used as the paper reprints them.

# %%
# prose literals of the notebook text
for lit, name in [
    ("304", "paper_lightoff_window_1D_low"), ("343", "paper_lightoff_window_1D_high"),
    ("293", "paper_lightoff_window_2D_low"), ("0.029", "paper_Qw_zero_stated"),
    ("0.20", "paper_fig2_abscissa_end"), ("0.50", "paper_fig11_abscissa_end"),
    ("0.25", "paper_fig3_abscissa_end"), ("0.0224", "paper_Tbar_at_X_0p25"),
    ("2.976", "paper_table1_square_Nu_T"), ("2.9775", "nu_t_square_pymrm"),
    ("2", "paper_adequacy_claim_degC"),
]:
    claims.says(lit, name)
claims.says("8.7 %", "nu_pole_downstream_fig9")
claims.says("8.7 %", "nu_pole_downstream_fig8")
claims.says("0.6", "deferred_omega")
claims.says("0.30", "sweep_L_margin")

# %%
# meta.yaml (validation, adds, agreement, caveats) and the models.yaml entry
for lit, name in [
    ("3.607950", "nu_h1_square_pymrm"), ("3.608", "paper_table1_square_Nu_H"),
    ("2.977523", "nu_t_square_pymrm"), ("2.976", "paper_table1_square_Nu_T"),
    ("below 1e-6", "nu_t_circle_vs_A315"), ("8.7 %", "nu_pole_downstream_fig9"),
    ("8.7 %", "nu_pole_downstream_fig8"), ("2.2", "pole_paper_counts_short_min_pct"),
    ("2.4", "pole_paper_counts_short_max_pct"), ("0.3", "pole_paper_counts_max_dev_pct"),
    ("0.0224", "paper_Tbar_at_X_0p25"), ("5.9 %", "Tbar_X025_paper_grid_error"),
    ("0.021147", "Tbar_X025_converged"), ("1.26 %", "grigull_tratz_max_reldev"),
    ("below 1e-9", "graetz_newton_vs_deferred_field"),
    ("below 1e-9", "graetz_newton_vs_deferred_pole"),
    ("0.96", "G_lightoff_2d_316_dG_order"), ("0.93", "lightoff_march_vs_quad_max_steps"),
    ("0.97 %", "quad_CG_endpoint_error"), ("-10.83", "T_windows_length_artefact_C"),
    ("-13.42", "T_windows_length_term_C"), ("+2.59", "T_windows_model_term_C"),
    ("+4.13", "T_low_2d_minus_1d_const_L025"), ("0.25", "paper_fig3_abscissa_end"),
    ("+2.85", "equiv_shift_GT_at316"), ("316", "paper_TG_inlet_fig1"),
    ("2.86", "equiv_shift_const_over_claim_at316"), ("2.26", "equiv_shift_const_over_claim_min"),
    ("7.74", "equiv_shift_const_over_claim_max"), ("1.43", "equiv_shift_GT_over_claim_at316"),
    ("1.09", "equiv_shift_GT_over_claim_min"), ("3.41", "equiv_shift_GT_over_claim_max"),
    ("0.028680", "fluxzero_X_fig9"), ("0.028789", "fluxzero_X_fig8"),
    ("0.029", "paper_Qw_zero_stated"), ("0.003", "paper_Cw_at_rate_maximum"),
    ("0.00261", "rate_max_Cw"), ("+5.72", "equiv_shift_const_Nu_at316"),
    ("2", "paper_adequacy_claim_degC"), ("22.3", "recon_spread_upper_threshold_C"),
    ("4.55", "recon_spread_equiv_shift_C"), ("303.5", "T_low_1d_const_L025"),
    ("304", "paper_lightoff_window_1D_low"), ("0.20", "paper_fig2_abscissa_end"),
    ("0.30", "sweep_L_margin"), ("0.50", "paper_fig11_abscissa_end"),
    ("1e-4", "axial_step_dG"), ("0.032", "paper_X_sign_negative_Nu"), ("0.033", "paper_X_sign_positive_Nu"),
    ("33", "data_rows"), ("15", "rows_naming_a_figure"), ("11", "rows_inside_a_frame"),
]:
    claims.says(lit, name)

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
