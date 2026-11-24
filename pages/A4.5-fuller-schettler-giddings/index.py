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
# title: "Fuller-Schettler-Giddings: a correlation and the difference between fitting and predicting"
# description: "The most-used estimate of binary gas diffusivity, rebuilt from its own printed tables. Its famous 4.3 % is a goodness of fit to the measurements the diffusion volumes were fitted to; held out one organic compound at a time, the same scheme does worse."
# categories: [sec:A, struct:S3, tier:T0, data:tier5, phase:gas]
# date: 2026-08-05
# ---
#
# # Fuller-Schettler-Giddings binary gas diffusivities
#
# **Catalog ID:** `A4.5` · **Structures:** `S3` (scalar diffusion) · **Tier:** T0
#
# $$
# D_{AB} \;=\; \frac{1.00\times10^{-3}\,T^{1.75}\,(1/M_A + 1/M_B)^{1/2}}
#                   {p\left[(\Sigma_A v_i)^{1/3} + (\Sigma_B v_i)^{1/3}\right]^{2}}
# $$
#
# A temperature, two molar masses and a table of "diffusion volumes" give a
# binary gas diffusivity. The estimate is in most process simulators and
# transport textbooks, and the accuracy quoted for it is the average error the
# paper prints, `4.32` per cent.
#
# This page establishes what that number is. Fuller, Schettler and Giddings
# obtained the diffusion volumes by nonlinear least squares **on the same
# measurements they then report the average error against**, and they say so:
# *"Despite the large sample, this is not a totally objective test since the
# parameters were obtained from these particular data."* So it is a **goodness
# of fit, not a validation**.
#
# The paper contains no held-out data, but it prints enough to construct a
# held-out test for one half of the method: the atomic-increment scheme that
# gives $\Sigma v$ for organic vapours is additive, so a compound can be removed
# from the fitting set and still be predicted from the increments the other
# compounds fix. Doing that one compound at a time is the page's main result.
# The other half of the method, the simple molecules with one fitted $\Sigma v$
# each, admits no held-out test, and the page says so.
#
# Everything comes from the printed page images of the paper. Nothing is taken
# from the later versions of the diffusion-volume table.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.24"  # replaced by the release tag at publication
PYMRM = ">=3.4.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "A4.5-fuller-schettler-giddings"  # equals the page directory

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
# Imports. Packages beyond pymrm and the Colab defaults (numpy, scipy, pandas,
# matplotlib, pyyaml) are installed here with subprocess, not %pip.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pymrm.checks import observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import ORG, fsg_params, sigma_v, D_fsg, VolumeFit, Tube, flux_quadrature

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
pd.set_option("display.width", 120)
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# **The paper.** Fuller, E. N., Schettler, P. D. and Giddings, J. C.,
# *"A New Method for Prediction of Binary Gas-Phase Diffusion Coefficients"*,
# **Industrial and Engineering Chemistry 58(5), 18-27 (1966)**,
# [doi:10.1021/ie50677a007](https://doi.org/10.1021/ie50677a007), read in full
# from the PDF; its title block, the ACS download stamp and the PDF's DOI
# metadata agree with the catalogue record.
#
# **Where it starts.** From the hard-sphere kinetic theory result the paper
# prints as its eq. (1) in "the Chapman-Enskog form", which carries a $T^{3/2}$
# temperature dependence and needs a collision diameter $\sigma_{AB}$. The paper
# states two objections to using it directly:
#
# 1. *"The first limitation is the 3/2 power temperature dependence; observed
#    values usually lie in the range 1.6 to 1.8."* Section 5 tests that sentence
#    against the paper's own data.
# 2. Few $\sigma$ are tabulated, and those apply over narrow temperature ranges.
#
# Arnold (1930) estimated $\sigma$ from additive Le Bas atomic volumes. Fuller
# *et al.* keep the additivity and drop the reading of the parameters as
# volumes: they fit a general form (their eq. 2) with a free temperature
# exponent $b$, free exponents $\alpha_1,\alpha_2,\alpha_3$ on the two sums and
# one free additive parameter $v_i$ per atom and per simple molecule, to the
# measured diffusivities of their Table III. The exponents came out close to the
# hard-sphere values and were fixed at $1/3, 1/3, 2$; the paper reports $b$ as
# "1.749 ± 0.013" and rounds it to 1.75. What survives is eq. (4) above and the
# table of $v_i$.
#
# **Why that matters for what the page may claim.** The $v_i$ are fitted
# parameters, not measured volumes, and the paper says what follows:
#
# > *"The present method can be no more reliable than the diffusion data used in
# > its development. Any systematic errors in the data would be reflected in the
# > atomic diffusion volumes."*
#
# **The erratum that is not on disk.** A published erratum exists,
# [doi:10.1021/ie50680a601](https://doi.org/10.1021/ie50680a601), *"Errata: A
# New Method for Predicting Gas Diffusivities"*, *Ind. Eng. Chem.* **58**
# (1966), and it has **not been consulted**. The article on disk contains no
# correction notice. An erratum to this paper would most plausibly correct the
# atomic diffusion volumes, the content of the correlation. This page uses the
# 1966 table exactly as printed. A different set of increments (C 15.9, H 2.31,
# O 6.11, N 4.54) is widely reproduced for this correlation; it is in no
# document on disk, is not asserted here to be the erratum's content, and
# appears only as an injected defect in the break table.
#
# **What this page is not.** `A4.6` (Chapman-Enskog, from Chapman & Cowling)
# prints measured $D_{12}$ pairs; setting the kinetic-theory result against this
# empirical one on one measured axis is the natural comparison and `A4.6`
# should own it. No number on this page comes from Chapman & Cowling.

# %% [markdown]
# ## The published model
#
# **Eq. (2), the general form fitted** (journal p. 20):
#
# $$
# D_{AB} = \frac{C\,T^{b}\,(1/M_A + 1/M_B)^{1/2}}
#               {p\left[(\Sigma_A v_i)^{\alpha_1} + (\Sigma_B v_i)^{\alpha_2}\right]^{\alpha_3}}
# $$
#
# **Eq. (3), the objective** (journal p. 20). The paper prints
# $\phi' = \sum_{i=1}^{N}[Y_i - E_i/Y_i]^2$, a typographic slip: the definition
# list below it reads *"sum of the squares of the relative differences between
# observed and estimated values"*, so
#
# $$
# \phi' = \sum_{i=1}^{N}\left(\frac{Y_i - E_i}{Y_i}\right)^{2},
# $$
#
# with $Y_i$ the observed $D_{AB}$ and $E_i$ the estimate. The paper's reason for
# the relative form is that the data span three orders of magnitude while the
# fractional experimental error is roughly constant.
#
# **Eq. (4), the published correlation**, obtained by fixing
# $\alpha_1=\alpha_2=1/3$, $\alpha_3=2$, $C=1.00\times10^{-3}$, $b=1.75$ and
# re-optimising the $v_i$: the equation at the top of the page, with $D_{AB}$ in
# cm²/s, $T$ in K, $p$ in atm and $M$ in g/mol. $C$ was chosen *"so that the
# atomic diffusion volumes would roughly correspond to Le Bas volumes"*.
#
# **Eq. (5), the error definition** of Table III, used everywhere on this page:
#
# $$
# \%\ \text{error} = \frac{D_{\text{calcd.}} - D_{\text{obsd.}}}{D_{\text{obsd.}}}\times 100
# $$
#
# | symbol | meaning | code |
# |---|---|---|
# | $C$, $b$ | prefactor, temperature exponent | `PAR["C"]`, `PAR["b"]` |
# | $v_i$ | Table I volumes and increments | `PAR["v_atom"]`, `PAR["v_ring"]`, `PAR["v_simple"]` |
# | $M_A$, $M_B$ | molar masses | `model.mw`, `PAR["mw_simple"]` |
# | $\phi'$ | eq. (3) | `VolumeFit.resid` |
#
# Equations 1 to 5 were read from 300 ppi renders of the page images, not from
# the text layer, which mangles the equations and, as the data section shows,
# also alters digits inside the data table.

# %% [markdown]
# ## Parameters and assumptions
#
# **Table I, the diffusion volumes**, transcribed from a 300 ppi crop of the page
# image (journal page 21). Parentheses in the printed table mark a value *"based
# on only a few data points"* and are carried into the CSV as `few_points_flag`.
# The atomic and structural increments are summed for organic vapours; the
# simple molecules each carry one fitted $\Sigma v$, the paper's reason being
# that *"typical organic substances tend to obey additivity relations much
# better than substances such as those listed in the lower part."*
#
# Two things the paper does not state were settled from its own printed error
# column. Both are inferences, and the rows they rest on are counted out of the
# independent tests in Section 1.
#
# - **Which molecules get the ring term.** Table I calls it "aromatic or
#   heterocyclic rings". Applied to the non-aromatic heterocycles as well
#   (piperidine, tetrahydrothiophene, ethylene oxide), it leaves those rows far
#   from their printed errors; withheld, they reproduce. The term as the authors
#   applied it covers **aromatic rings only**. `PAR["ring_extra"]` lists species
#   given a ring term beyond their aromatic rings; it is empty as published, and
#   the break table fills it.
# - **The molar mass of D₂.** With the isotopic mass, 4.028, the D₂ rows move
#   away from their printed errors; with $M = 4.00$ they reproduce. Section 1
#   prints both alternatives.
#
# Molar masses are modern IUPAC values; the paper prints none. Every published
# constant sits in the one dictionary `PAR` that every function takes.

# %%
VOL = load_data("fuller-1966-table1-diffusion-volumes.csv", page=PAGE)
PAR = fsg_params(VOL)          # C = 1.00e-3, b = 1.75, Table I volumes, molar masses
print(VOL.to_string(index=False))
n_vol = len(PAR["v_atom"]) + 1 + len(PAR["v_simple"])
print(f"\n{len(PAR['v_atom'])} atomic increments + 1 ring increment + "
      f"{len(PAR['v_simple'])} simple molecules = {n_vol} fitted volumes")
print(f"spot check: H2-N2 at 273.0 K and 1 atm, D = {D_fsg('H2', 'N2', 273.0, PAR):.4f} cm2/s; "
      f"benzene Sigma v = {sigma_v('BENZENE', PAR):.2f}")

# %% [markdown]
# ## The data
#
# Two datasets, both transcribed from the paper's printed tables and owned by
# this page; nothing is digitised from a figure.
#
# - `fuller-1966-table1-diffusion-volumes.csv`: Table I, the correlation's
#   content (printed above).
# - `fuller-1966-table3-diffusion-data.csv`: **Table III in full**, each
#   measured binary diffusivity with the literature reference number the paper
#   cites, the temperature, the observed $D_{AB}$ *"reduced to 1 atm. pressure
#   (cm²/sec.)"* and the FSG percentage error the authors printed beside it.
#   These are measurements (tier 5) collected by Fuller *et al.* from the
#   literature. The flags `tracer` (rows marked `(TR)`) and `starred` (force
#   constants estimated from critical properties, relevant only to other
#   methods' columns) are kept so the table is complete.
#
# **Only the FSG error column was transcribed.** Table III also carries eight
# other methods' errors; they are not in the CSV.
#
# **How the table was extracted.** The scan is 300 ppi native. Rows were located
# by word bounding boxes of the text layer and column ranges, then every row was
# checked against eq. (4): recomputed from $T$, the species and Table I, it must
# reproduce the printed FSG error, so a wrong digit larger than the clean-row
# scatter cannot hide. Rows that failed were re-read on a crop of the page image
# and corrected; the column `image_corrected` names the field, and the cell
# below counts them. A slip below that scatter can hide: the sidecar
# `data/fuller-1966-table3-diffusion-data.meta.yaml` records one such slip in
# the error column's last digit, found by re-reading the crop, and the audit
# sample of further rows read on crops.
#
# **What the paper says about its own data:**
#
# - *"Accurate experimental diffusion measurements are difficult to make... The
#   typical amount of error is probably of the order of 5 %, as can be seen from
#   the fact that observations reported in the literature for identical systems
#   often differ by this amount or more."*
# - *"Even considering all available data, no completely adequate test of
#   correlation methods appears possible."*
# - The eight methods compared were evaluated *"for all the experimental points
#   used in fixing parameters for the present correlation"*. Every comparison in
#   Tables II and III is therefore in sample for FSG and out of sample for the
#   other eight.
#
# **A printed inconsistency in Table III.** The two summary rows at its foot,
# `AVE. ERROR (PER CENT)` and `STANDARD DEVIATION (PER CENT)`, read
# `4.32 21.11 16.57 10.67 6.38 6.76 8.86 14.02 13.89` and
# `6.71 24.34 19.24 15.41 8.72 8.91 10.95 23.73 17.67`. Those nine averages
# follow **Table II's** method order (FSG, Gilliland, Andrussow, Arnold, HBS,
# Mod-HBS, Chen-Othmer, Slattery-Bird, Othmer-Chen), not Table III's own column
# headings (FSG, ANDR, HBS, CH-OT, OT-CH, GILL, ARN, MOD-HBS, SL-BD). One of the
# two orderings is wrong. The page uses only the FSG column, whose identity
# Section 1 settles, so it does not resolve which; anyone reading the other
# eight columns should settle the ordering first.

# %%
DATA = load_data("fuller-1966-table3-diffusion-data.csv", page=PAGE)
META3 = load_meta("fuller-1966-table3-diffusion-data.csv", page=PAGE)
print(cite_data(META3))
SPECIES = set(DATA.component_a) | set(DATA.component_b)
n_sys = DATA.groupby(["component_a", "component_b"]).ngroups
n_fixed = int((DATA.image_corrected.fillna("") != "").sum())
print(f"\n{len(DATA)} rows, {n_sys} distinct binary systems, {len(SPECIES)} species")
print(f"temperature {DATA.temperature_K.min():.1f} to {DATA.temperature_K.max():.1f} K, "
      f"D_obs {DATA.D_obs_cm2_s.min():.4f} to {DATA.D_obs_cm2_s.max():.4f} cm2/s "
      f"({DATA.D_obs_cm2_s.max() / DATA.D_obs_cm2_s.min():.0f}-fold)")
print(f"rows corrected against the page image: {n_fixed}")
print("the paper states '153 different binary systems (340 actual data points)'")
claims.exact("n_rows", len(DATA))
claims.exact("n_systems", n_sys, note="the paper states 153")
claims.exact("n_rows_image_corrected", n_fixed, note="count of the image_corrected column")
display(DATA.head(8))

# %% [markdown]
# ## PyMRM implementation
#
# A correlation is not a transport model, so the pymrm content here is what the
# correlation is for: putting $D_{AB}(T)$ into a transport calculation and
# measuring what its exponent costs.
#
# The problem is steady one-dimensional diffusion of a dilute species through a
# stagnant gas in a tube whose wall temperature varies linearly from $T_0$ to
# $T_1$, the situation in which the exponent $b$ matters rather than cancelling.
# In mole fractions,
#
# $$
# \frac{\mathrm{d}}{\mathrm{d}z}\!\left(-\Gamma(z)\,\frac{\mathrm{d}x}{\mathrm{d}z}\right) = 0,
# \qquad
# \Gamma(z) = c_t(z)\,D_{AB}(T(z)) = \frac{p}{RT}\,D_{\mathrm{ref}}\!\left(\frac{T}{T_\mathrm{ref}}\right)^{b},
# $$
#
# so $\Gamma \propto T^{\,b-1}$: the diffusivity rises with temperature, the
# total concentration falls, and the two do not cancel.
#
# `model.Tube` assembles `construct_grad`, `construct_div` and
# `construct_coefficient_matrix` once per grid. **The problem is linear in $x$**:
# one direct sparse solve (`spsolve`), no Newton iteration. $\Gamma$ at every
# face, the two end faces included, is $\Gamma$ evaluated at the face
# temperature, one rule for all faces. The boundary conditions are Dirichlet at
# both ends, written on the outward normal in the code ($n=-z$ at $z=0$, $n=+z$
# at $z=L$; with $a=0$ the sign of the normal does not enter). The flux is read
# at $z=0$ from the boundary gradient `compute_boundary_values` returns.
#
# **The validation route shares no code with the solver.** The flux is constant
# along the tube and
#
# $$
# N = \frac{x_0 - x_L}{\int_0^L \mathrm{d}z / \Gamma(z)},
# $$
#
# evaluated by adaptive quadrature (`scipy.integrate.quad`), which forms no grid
# and calls no pymrm. With $\Gamma$ at the face temperature on every face, end
# faces included, the scheme is second order; it is not the trapezoidal rule of
# this integral, and the quadrature shares only $\Gamma(T)$ with it.
#
# The case: benzene vapour through air in a 10 cm tube from 300 K to 600 K at
# 1 atm, $x_0 = 0.10$ and $x_L = 0$, with $D$ from eq. (4). Two ways of pinning
# $D$: at the cold end ($T_\mathrm{ref} = T_0$) and at the geometric mean
# temperature $\sqrt{T_0 T_1}$, called mid-range below.

# %%
L_TUBE, T_COLD, T_HOT = 10.0, 300.0, 600.0   # cm, K, K: page choices
T_MID = float(np.sqrt(T_COLD * T_HOT))       # K, geometric mean, the mid-range anchor
N_GRID = 400                                 # cells of the reported solve
LADDER = (25, 50, 100, 200, 400, 800)        # nested grids for the order study


def tube_case(T_ref, b=1.75, T_hot=T_HOT, par=PAR):
    """Keyword arguments of Tube and flux_quadrature with D from eq. (4) at T_ref."""
    return dict(L=L_TUBE, T0=T_COLD, T1=T_hot, D_ref=D_fsg("AIR", "BENZENE", T_ref, par),
                T_ref=T_ref, b=b)


case = tube_case(T_MID)
sol = Tube(N_GRID, **case).solve()
N_quad, quad_err = flux_quadrature(**case)
print(f"D(benzene-air) = {D_fsg('AIR', 'BENZENE', T_COLD, PAR):.4f} cm2/s at {T_COLD:.0f} K, "
      f"{case['D_ref']:.4f} at {T_MID:.1f} K")
print(f"flux, pymrm on {N_GRID} cells: {sol['flux']:.6e} mol/(cm2 s); quadrature "
      f"{N_quad:.6e}, its own relative error estimate {quad_err:.1e}")

# %% [markdown]
# ## Results
#
# ### 1. Every printed FSG error recomputed, and the rows that do not reproduce
#
# Eq. (4) with the Table I volumes, on every row, against the error the authors
# printed beside it. A row is called irreproducible when the two differ by more
# than half a percentage point; the cell prints the margin of that threshold.
# The two rows it finds were re-read at high magnification on the page image and
# are as printed. They are named below and excluded by identity from the
# statistics of the reproducible rows.

# %%
THRESH_PP = 0.5          # pp: |recomputed - printed| above this is irreproducible
NAMED = [("HE", "WATER", 352.5), ("CO2", "ETHYLENEOXIDE", 298.0)]   # (A, B, T in K)
NAMED_MASK = np.logical_or.reduce(
    [((DATA.component_a == a) & (DATA.component_b == b) & (DATA.temperature_K == T)).values
     for a, b, T in NAMED])
INFERRED = ["PIPERIDINE", "TETRAHYDROTHIOPHENE", "ETHYLENEOXIDE", "D2"]  # rows of the two inferences


def data_study(drop_page=None):
    """Counts of the dataset; drop_page removes one printed page of Table III."""
    d = DATA if drop_page is None else DATA[DATA.pdf_page != drop_page]
    inferred = d.component_a.isin(INFERRED) | d.component_b.isin(INFERRED)
    return dict(n_rows=len(d), n_systems=d.groupby(["component_a", "component_b"]).ngroups,
                n_rows_image_corrected=int((d.image_corrected.fillna("") != "").sum()),
                n_rows_inferred=int(inferred.sum()),
                n_rows_d2=int(((d.component_a == "D2") | (d.component_b == "D2")).sum()))


def recompute(par=PAR):
    """Eq. (5) error of eq. (4) on every row of Table III, per cent."""
    D = np.array([D_fsg(a, b, T, par) for a, b, T in
                  zip(DATA.component_a, DATA.component_b, DATA.temperature_K)])
    return 100 * (D - DATA.D_obs_cm2_s.values) / DATA.D_obs_cm2_s.values


HETERO = ("PIPERIDINE", "TETRAHYDROTHIOPHENE", "ETHYLENEOXIDE")   # non-aromatic heterocycles


def repro_study(par=PAR):
    delta = recompute(par) - DATA.err_fsg_pct.values
    ad = np.abs(delta)
    bad = ad > THRESH_PP
    inferred = (DATA.component_a.isin(INFERRED) | DATA.component_b.isin(INFERRED)).values
    hetero = (DATA.component_a.isin(HETERO) | DATA.component_b.isin(HETERO)).values
    return dict(n_hetero_reproduced=int((hetero & ~bad).sum()),
                n_rows_irreproducible=int(bad.sum()), n_rows_reproduced=int((~bad).sum()),
                irreproducible_are_named=bool(np.array_equal(bad, NAMED_MASK)),
                max_abs_delta_pp=float(ad[~NAMED_MASK].max()),
                median_abs_delta_pp=float(np.median(ad[~NAMED_MASK])),
                n_free_rows_failing=int((bad & ~inferred).sum()),
                _margin=float(np.abs(ad - THRESH_PP).min() / THRESH_PP), _delta=delta)


REP = repro_study()
DATA["err_calc"] = DATA.err_fsg_pct + REP["_delta"]
DATA["delta"] = REP["_delta"]
print(f"{REP['n_rows_reproduced']} of {len(DATA)} rows reproduce the printed error to "
      f"|delta| <= {REP['max_abs_delta_pp']:.3f} pp (median {REP['median_abs_delta_pp']:.4f} pp); "
      f"threshold {THRESH_PP} pp, nearest row to it at {REP['_margin']:.2f} of the threshold")
print(f"{REP['n_rows_irreproducible']} rows do not; they are the two named rows: "
      f"{REP['irreproducible_are_named']}")
display(DATA.loc[NAMED_MASK, ["component_a", "component_b", "ref_number", "temperature_K",
                              "D_obs_cm2_s", "err_fsg_pct", "err_calc", "delta"]]
        .rename(columns={"err_fsg_pct": "printed", "err_calc": "recomputed"}))

def inference_study(par=PAR):
    """Each inference's alternative applied to par, on the rows it concerns (named rows out)."""
    out = {}
    for tag, alt, sp in (("ring", dict(par, ring_extra=HETERO), HETERO),
                         ("d2_isotopic", dict(par, mw_simple={**par["mw_simple"], "D2": 4.028}),
                          ("D2",))):
        m = (DATA.component_a.isin(sp) | DATA.component_b.isin(sp)).values & ~NAMED_MASK
        d = recompute(alt)[m] - DATA.err_fsg_pct.values[m]
        out.update({f"{tag}_delta_min_pp": float(d.min()), f"{tag}_delta_max_pp": float(d.max()),
                    f"_{tag}_n": int(m.sum())})
    return out


INF = inference_study()
for tag, label in (("ring", "ring term on the non-aromatic heterocycles"),
                   ("d2_isotopic", "M(D2) = 4.028")):
    print(f"{label}: {INF[f'_{tag}_n']} rows, delta {INF[f'{tag}_delta_min_pp']:+.2f} to "
          f"{INF[f'{tag}_delta_max_pp']:+.2f} pp")
    for k in ("min", "max"):
        claims.value(f"{tag}_delta_{k}_pp", INF[f"{tag}_delta_{k}_pp"], rtol=1e-6, unit="pp")
d2 = ((DATA.component_a == "D2") | (DATA.component_b == "D2")).values
print(f"M(D2) = 4.00, as used: {d2.sum()} rows, delta {DATA.delta[d2].min():+.3f} to "
      f"{DATA.delta[d2].max():+.3f} pp")

DS = data_study()
n_inf = DS["n_rows_inferred"]
print(f"\nrows the two inferences rest on: {n_inf}; free tests: {len(DATA) - n_inf}, of which "
      f"{REP['n_free_rows_failing']} irreproducible")
claims.exact("n_rows_d2", DS["n_rows_d2"])
print(f"rows with a non-aromatic heterocycle: {REP['n_hetero_reproduced']} reproduce without "
      f"the ring term")
claims.exact("n_hetero_reproduced", REP["n_hetero_reproduced"])
for k in ("n_rows_irreproducible", "n_rows_reproduced", "n_free_rows_failing"):
    claims.exact(k, REP[k], margin=REP["_margin"], min_margin=0.1,
                 note=f"|delta| against {THRESH_PP} pp; margin relative to it")
claims.exact("irreproducible_are_named", REP["irreproducible_are_named"])
claims.exact("n_rows_inferred", n_inf)
claims.value("max_abs_delta_pp", REP["max_abs_delta_pp"], rtol=1e-6, unit="pp",
             note="over the rows other than the two named ones")
claims.value("median_abs_delta_pp", REP["median_abs_delta_pp"], rtol=1e-6, unit="pp")

# %% [markdown]
# **Route B, a check of the same rows that uses no volumes, no molar masses and
# no evaluation of eq. (4).** Within one binary system from one literature
# reference, eq. (4) makes $D_\text{calc}$ scale exactly as $T^{1.75}$, so the
# printed temperature, the printed observed $D$ and the printed error alone must
# satisfy that identity between consecutive temperatures. Half a unit in the last
# printed digit of $D_\text{obs}$ (four significant figures) and of the error, on
# both rows of a ratio, moves it by about 1e-3 at most; that is the bound used.
# The exponent recovered from the printed calculated values is printed with and
# without the ratio that touches a named row.

# %%
ROUND_ID = 1e-3          # rounding bound of one ratio of printed values


def routeB_study(T_offset=0.0, b_id=1.75, by_ref=True):
    d = DATA.assign(T=DATA.temperature_K + T_offset,
                    Dc=DATA.D_obs_cm2_s * (1 + DATA.err_fsg_pct / 100))
    keys = ["component_a", "component_b"] + (["ref_number"] if by_ref else [])
    r, t, touch, n_series = [], [], [], 0
    for _, s in d.groupby(keys):
        s = s.sort_values("T")
        if s["T"].nunique() < 2:
            continue
        n_series += 1
        for i, j in zip(s.index[:-1], s.index[1:]):
            if s.at[i, "T"] == s.at[j, "T"]:
                continue
            r.append(s.at[j, "Dc"] / s.at[i, "Dc"])
            t.append(s.at[j, "T"] / s.at[i, "T"])
            touch.append(NAMED_MASK[i] or NAMED_MASK[j])
    r, t, touch = np.array(r), np.array(t), np.array(touch)
    dev = np.abs(r / t**b_id - 1)
    slope = lambda m: float(np.sum(np.log(r[m]) * np.log(t[m])) / np.sum(np.log(t[m])**2))
    beyond = dev > ROUND_ID
    return dict(routeB_n_ratios=len(r), routeB_n_series=n_series,
                routeB_recovered_exponent=slope(np.ones_like(touch)),
                routeB_recovered_exponent_excl=slope(~touch),
                routeB_median_dev=float(np.median(dev)),
                routeB_max_dev_excl_anomaly=float(dev[~touch].max()),
                routeB_n_beyond_rounding=int(beyond.sum()),
                routeB_beyond_are_named=bool(np.array_equal(beyond, touch)),
                _margin=float(np.abs(dev - ROUND_ID).min() / ROUND_ID), _dev=dev, _touch=touch)


RB = routeB_study()
print(f"{RB['routeB_n_ratios']} consecutive-temperature ratios in {RB['routeB_n_series']} "
      f"(system, reference) series; {int(RB['_touch'].sum())} touches a named row")
print(f"deviation from the T^1.75 identity: median {RB['routeB_median_dev']:.1e}, largest "
      f"away from the named rows {RB['routeB_max_dev_excl_anomaly']:.1e}, on the named row "
      f"{RB['_dev'][RB['_touch']].max():.3f}")
print(f"ratios beyond the rounding bound {ROUND_ID:g}: {RB['routeB_n_beyond_rounding']}, "
      f"exactly those touching a named row: {RB['routeB_beyond_are_named']}")
print(f"b recovered from the printed calculated values: {RB['routeB_recovered_exponent']:.5f} "
      f"with every ratio, {RB['routeB_recovered_exponent_excl']:.5f} without the ratios "
      f"touching a named row")
claims.exact("routeB_n_ratios", RB["routeB_n_ratios"])
claims.exact("routeB_n_series", RB["routeB_n_series"])
claims.value("routeB_recovered_exponent", RB["routeB_recovered_exponent"], rtol=1e-6,
             note="all ratios")
claims.value("routeB_recovered_exponent_excl", RB["routeB_recovered_exponent_excl"], rtol=1e-6,
             note="ratios touching a named row excluded")
claims.bound("routeB_median_dev", RB["routeB_median_dev"], max=ROUND_ID, basis="theory",
             note="rounding of the printed D_obs and error")
claims.bound("routeB_max_dev_excl_anomaly", RB["routeB_max_dev_excl_anomaly"], max=ROUND_ID,
             basis="theory", note="rounding of the printed D_obs and error")
claims.exact("routeB_n_beyond_rounding", RB["routeB_n_beyond_rounding"],
             margin=RB["_margin"], min_margin=0.1)
claims.exact("routeB_beyond_are_named", RB["routeB_beyond_are_named"])

# %% [markdown]
# ### 2. Table II and the printed aggregates at the foot of Table III
#
# Table II's category memberships are not defined in the paper; a row is
# assigned to a category when either component is that species. Each printed
# aggregate is compared in half-units of its last printed digit (one half-unit
# is the most rounding can explain). The point counts and means of Table II are
# this page's reading of the printed table.

# %%
# Table II as printed: category, species, points, mean |error| % (FSG column)
CATS = [("H2", ("H2",), "59", "5.8"), ("He", ("HE",), "108", "3.9"),
        ("N2", ("N2",), "71", "3.8"), ("O2", ("O2",), "23", "2.3"),
        ("Ar", ("AR",), "58", "3.3"), ("Air", ("AIR",), "53", "5.7"),
        ("NH3, H2O", ("NH3", "WATER"), "41", "4.4"), ("Organic", tuple(ORG), "154", "5.0")]
FOOT = {"ave": "4.32", "tail": "7.4", "sd": "6.71"}   # Table III foot; Table II tail
PHI_PRINTED_MIN = 1.4762      # the paper's constrained minimum of eq. (3), b = 1.75
N_PAR = n_vol + 1             # Table I rows plus b


def half_units(x, printed):
    dec = len(printed.split(".")[1]) if "." in printed else 0
    return abs(x - float(printed)) / (0.5 * 10.0**-dec)


ERR_DIGITS = int(max(len(x.split(".")[1]) if "." in x else 0 for x in
                     load_data("fuller-1966-table3-diffusion-data.csv", page=PAGE,
                               dtype=str).err_fsg_pct))      # decimals of the printed errors


def aggregates_study(err_scale=1.0, err_basis="obs", n_par=N_PAR, tail=10.0, member="either",
                     err_digits=ERR_DIGITS):
    e = DATA.err_fsg_pct.values * err_scale
    if err_basis == "calc":                      # (calc - obs)/calc in place of eq. (5)
        e = e / (1 + e / 100)
    rows = []
    for name, sp, n_p, m_p in CATS:
        a, b = DATA.component_a.isin(sp).values, DATA.component_b.isin(sp).values
        mask = (a | b) if member == "either" else a
        rows.append((name, int(mask.sum()), int(n_p), np.abs(e[mask]).mean(), m_p))
    T2 = pd.DataFrame(rows, columns=["category", "n", "n printed", "mean |error| %", "printed"])
    T2["half-units"] = [half_units(x, p) for x, p in zip(T2["mean |error| %"], T2.printed)]
    ave, pct = np.abs(e).mean(), 100 * np.mean(np.abs(e) > tail)
    sd = np.sqrt(np.sum(e**2) / (len(e) - n_par))
    phi = np.sum((e / 100)**2)
    hu = list(T2["half-units"]) + [half_units(ave, FOOT["ave"]), half_units(pct, FOOT["tail"]),
                                   half_units(sd, FOOT["sd"])]
    return dict(ave_error_pct=ave, pct_over_10=pct, residual_sd_pct=sd,
                phi_prime_printed_column=phi, phi_prime_gap_printed_min=PHI_PRINTED_MIN - phi,
                aggregates_within_print=int(np.sum(np.array(hu) <= 1.0)),
                table2_counts_exact=int((T2.n == T2["n printed"]).sum()),
                aggregate_worst_halfunits=float(max(hu)),
                aggregate_input_rounding_halfunits=max(
                    0.5 * 10.0**-err_digits / (0.5 * 10.0**-len(p.split(".")[1]))
                    for _, _, _, _, p in rows),
                _margin=float(np.abs(np.array(hu) - 1.0).min()), _T2=T2, _hu=hu, _e=e)


AG = aggregates_study()
print(AG["_T2"].round(3).to_string(index=False))
e = AG["_e"]
print(f"\nAVE. ERROR {AG['ave_error_pct']:.4f} (printed {FOOT['ave']}), more than 10 %: "
      f"{AG['pct_over_10']:.2f} (printed {FOOT['tail']}), n = {len(e)}")
print(f"STANDARD DEVIATION (printed {FOOT['sd']}): sample sd of the signed errors "
      f"{np.std(e, ddof=1):.3f}, root mean square {np.sqrt(np.mean(e**2)):.3f}, residual sd "
      f"sqrt(sum e^2/(N - p)) with N = {len(e)}, p = {N_PAR}: {AG['residual_sd_pct']:.4f} "
      f"({len(e) - N_PAR} degrees of freedom)")
worst = AG["_T2"].category[np.argmax(AG["_T2"]["half-units"])]
print(f"aggregates within half a unit of the printed digit: {AG['aggregates_within_print']} of "
      f"{len(AG['_hu'])}; the other, {worst}, at {AG['aggregate_worst_halfunits']:.3f} half-units. "
      f"Half a unit of the last printed error digit moves a mean by at most "
      f"{0.5 * 10.0**-ERR_DIGITS:g} pp, {AG['aggregate_input_rounding_halfunits']:.1f} "
      f"half-units of a one-decimal mean")
print(f"Table II point counts equal: {AG['table2_counts_exact']} of {len(CATS)}")
claims.value("aggregate_worst_halfunits", AG["aggregate_worst_halfunits"], rtol=1e-6,
             note="half-units of the printed digit, worst of the 11 aggregates")
claims.value("aggregate_input_rounding_halfunits", AG["aggregate_input_rounding_halfunits"],
             rtol=1e-6, note="half a unit of the last printed error digit, over half a unit of "
             "the coarsest printed mean digit; every error moves a mean by at most this")
claims.value("ave_error_pct", AG["ave_error_pct"], rtol=1e-6, unit="%",
             note="mean |error| of the printed column, n = 340")
claims.value("pct_over_10", AG["pct_over_10"], rtol=1e-6, unit="%", note="n = 340")
claims.value("residual_sd_pct", AG["residual_sd_pct"], rtol=1e-6, unit="%",
             note="sqrt(sum e^2/(N - p)), N = 340, p = 28, 312 degrees of freedom")
claims.exact("aggregates_within_print", AG["aggregates_within_print"], margin=AG["_margin"],
             min_margin=0.01, note="8 Table II means, AVE. ERROR, >10 % tail, STANDARD DEVIATION")
claims.exact("table2_counts_exact", AG["table2_counts_exact"])

# %% [markdown]
# The residual standard deviation needs a comment, because the paper does not
# print its definition. Neither the sample standard deviation of the signed
# errors nor their root mean square gives the printed value; the residual
# standard deviation with the printed parameter count does: every row of Table I
# plus the exponent $b$. The degrees-of-freedom form is an inference; the
# parameter count it needs is printed.
#
# ### 3. The objective function does not match the table
#
# $\phi'$ can be computed three ways, and none lands where the paper says its
# minimum is: the paper prints "1.4707" as the least-squares minimum with $b$
# free and "1.4762" as the constrained minimum at $b = 1.75$. Eq. (4) is a
# special case of eq. (2), so a minimum of eq. (2) cannot lie above a value
# reachable inside eq. (4). The cell measures the three values, the margin of
# each verdict, and whether the free refit is the minimum (several starts, the
# gradient at the solution, the bounds).

# %%
FIT = VolumeFit(DATA.component_a, DATA.component_b, DATA.temperature_K, DATA.D_obs_cm2_s, PAR)
N_STARTS = 8             # random starts, log-normal spread 0.3 about the published volumes
D_ROUND = 5e-4           # half a unit in the 4th significant figure of D_obs, relative


def fit_study(par=PAR, b=None, objective="relative", tol=1e-14, lower=None, seed=1966):
    F = FIT if par is PAR else VolumeFit(DATA.component_a, DATA.component_b,
                                         DATA.temperature_K, DATA.D_obs_cm2_s, par)
    lo = F.lower if lower is None else lower
    phi = lambda x: float(np.sum(F.resid(x, b=b)**2))
    start = np.maximum(F.published, lo + 1e-9)
    f = F.fit(start=start, b=b, objective=objective, tol=tol, lower=lo)
    rng = np.random.default_rng(seed)
    phis = []
    for _ in range(N_STARTS):
        s = np.maximum(F.published * np.exp(rng.normal(0.0, 0.3, F.published.size)), lo + 1e-9)
        phis.append(phi(F.fit(start=s, b=b, objective=objective, tol=tol, lower=lo).x))
    r_pub = F.resid(F.published, b=b)
    return dict(phi_prime_published_volumes=phi(F.published), phi_prime_free_refit=phi(f.x),
                phi_prime_gap_published=PHI_PRINTED_MIN - phi(F.published),
                refit_simple_volume_max_rel_dev=float(np.max(np.abs(f.x[7:] / F.published[7:] - 1))),
                refit_n_at_bound=int(np.sum(f.active_mask != 0)),
                refit_start_spread=float(np.max(np.abs(np.array(phis) - phi(f.x)))),
                refit_optimality=float(f.optimality),
                _round=float(np.sum(2 * np.abs(r_pub) * D_ROUND)), _fit=f, _F=F,
                _bound_gap=float(np.min(f.x - lo)))


FS = fit_study()
fx = FS["_fit"].x
e_round = float(np.sum(2 * np.abs(e / 100) * 5e-5))   # half a unit (0.005 pp) in every printed error
print(pd.DataFrame({"phi'": ["paper, least squares (b free)", "paper, constrained (b = 1.75)",
                             "Table III's printed error column", "eq. (4), printed Table I volumes",
                             f"free refit of all {n_vol} volumes, b = 1.75"],
                    "value": ["1.4707", "1.4762", f"{AG['phi_prime_printed_column']:.4f}",
                              f"{FS['phi_prime_published_volumes']:.4f}",
                              f"{FS['phi_prime_free_refit']:.4f}"]}).to_string(index=False))
print(f"\ngap to 1.4762: printed column {AG['phi_prime_gap_printed_min']:.4f} against at most "
      f"{e_round:.1e} from rounding of the printed errors; published volumes "
      f"{FS['phi_prime_gap_published']:.4f} against at most {FS['_round']:.1e} from rounding of D_obs")
print(f"free refit: {N_STARTS} random starts reach the same phi' to below 1e-10; "
      f"first-order optimality below 1e-06; volumes on a bound: "
      f"{FS['refit_n_at_bound']} (nearest {FS['_bound_gap']:.2f} above it)")
FIX = dict(zip(FS["_F"].names, fx))
print(f"refit against printed: simple molecules within {100 * FS['refit_simple_volume_max_rel_dev']:.1f} % "
      f"(median {100 * np.median(np.abs(fx[7:] / FIT.published[7:] - 1)):.1f} %), e.g. He "
      f"{FIX['HE']:.3f} vs {PAR['v_simple']['HE']}, SF6 {FIX['SF6']:.2f} vs {PAR['v_simple']['SF6']}; "
      f"atomic C {FIX['C']:.2f} vs {PAR['v_atom']['C']}, H {FIX['H']:.2f} vs {PAR['v_atom']['H']}")
claims.value("phi_prime_printed_column", AG["phi_prime_printed_column"], rtol=1e-6)
claims.value("phi_prime_published_volumes", FS["phi_prime_published_volumes"], rtol=1e-6)
claims.value("phi_prime_free_refit", FS["phi_prime_free_refit"], rtol=1e-6)
assert e_round < 2e-3 and FS["_round"] < 2e-2
claims.bound("phi_prime_gap_printed_min", AG["phi_prime_gap_printed_min"], min=2e-3,
             basis="theory", note="rounding of the printed errors moves phi' by less")
claims.bound("phi_prime_gap_published", FS["phi_prime_gap_published"], min=2e-2,
             basis="theory", note="rounding of the printed D_obs moves phi' by less")
claims.value("refit_simple_volume_max_rel_dev", FS["refit_simple_volume_max_rel_dev"], rtol=1e-4)
claims.exact("refit_n_at_bound", FS["refit_n_at_bound"], margin=FS["_bound_gap"],
             min_margin=1e-3, note="margin: distance of the nearest volume to its bound")
claims.bound("refit_start_spread", FS["refit_start_spread"], max=1e-10, basis="tolerance")
claims.bound("refit_optimality", FS["refit_optimality"], max=1e-6, basis="tolerance",
             note="least_squares first-order optimality")

# %% [markdown]
# All three computed values lie below the paper's constrained minimum, by far
# more than rounding of the printed digits can move them, so **the printed
# $\phi'$ figures do not belong to the published table.** Why is not resolvable
# from what is printed, and the page does not guess.
#
# The refit recovers the printed simple-molecule volumes closely and moves the
# atomic increments more, which is what strongly correlated parameters do; it is
# the strongest single sign that the model, the data and the objective on this
# page are the paper's.
#
# ### 4. The held-out test: what the additive scheme predicts for a compound it never saw
#
# This is the number the page exists for. The simple molecules each carry their
# own fitted $\Sigma v$, so removing one from the fitting set removes the
# parameter that predicts it: **those points admit no held-out test**. The
# organic-vapour points are different: $\Sigma v$ for an organic is built from
# the atomic increments, which the *other* organics fix. Leave-one-compound-out
# over every organic vapour, refitting all volumes each time on the remaining
# data, tests the additivity principle out of sample. The excess of held-out over
# in-sample error is a mean over compounds, so its standard error is taken over
# compounds (points of one compound are not independent).

# %%
def loco_study(par=PAR, leak=False, min_points=1, b=None):
    F = FIT if par is PAR else VolumeFit(DATA.component_a, DATA.component_b,
                                         DATA.temperature_K, DATA.D_obs_cm2_s, par)
    full = F.fit(b=b)
    rows, at_bound = [], []
    for c in sorted(s for s in SPECIES if s in ORG):
        held = ((DATA.component_a == c) | (DATA.component_b == c)).values
        if held.sum() < min_points:
            continue
        f = F.fit(mask=None if leak else ~held, b=b, tol=1e-12)
        if np.any(f.active_mask != 0):
            at_bound.append((c, [F.names[i] for i in np.flatnonzero(f.active_mask)]))
        err = lambda x: 100 * float(np.mean(np.abs(F.resid(x, b=b)[held])))
        rows.append((c, int(held.sum()), err(f.x), err(full.x), err(F.published)))
    T = pd.DataFrame(rows, columns=["compound", "n", "held_out_%", "in_sample_refit_%",
                                    "in_sample_published_%"])
    w = T.n.values
    d = (T["held_out_%"] - T["in_sample_refit_%"]).values
    gap, K = float(np.average(d, weights=w)), len(T)
    se = float(np.sqrt(np.sum(w**2 * (d - gap)**2)) / w.sum() * np.sqrt(K / (K - 1)))
    return dict(loco_held_out_pct=float(np.average(T["held_out_%"], weights=w)),
                loco_in_sample_refit_pct=float(np.average(T["in_sample_refit_%"], weights=w)),
                loco_in_sample_published_pct=float(np.average(T["in_sample_published_%"], weights=w)),
                loco_gap_pp=gap, loco_gap_se_pp=se, loco_gap_t=gap / se,
                loco_n_compounds=K, loco_n_points=int(w.sum()), loco_fits_at_bound=len(at_bound),
                _T=T, _at_bound=at_bound)


LO = loco_study()
LT = LO["_T"]
print(f"{LO['loco_n_compounds']} organic vapours, {LO['loco_n_points']} points, one refit each")
print(f"mean |error| on those points: published volumes in sample {LO['loco_in_sample_published_pct']:.2f} %, "
      f"full refit in sample {LO['loco_in_sample_refit_pct']:.2f} %, held out "
      f"{LO['loco_held_out_pct']:.2f} %")
print(f"held out minus in sample: {LO['loco_gap_pp']:.2f} pp, standard error over the "
      f"{LO['loco_n_compounds']} compounds {LO['loco_gap_se_pp']:.2f} pp "
      f"(t = {LO['loco_gap_t']:.1f}, {LO['loco_n_compounds'] - 1} degrees of freedom); relative "
      f"{100 * (LO['loco_held_out_pct'] / LO['loco_in_sample_refit_pct'] - 1):.0f} %")
print(f"refits with a volume on its bound: {LO['loco_fits_at_bound']}: {LO['_at_bound']}")
gain = LT.assign(gain=LT["held_out_%"] - LT["in_sample_refit_%"]).nlargest(4, "gain")
print("largest excess: " + ", ".join(f"{r.compound} {r[3]:.1f} -> {r[2]:.1f} %"
                                    for r in gain.itertuples(index=False)))
display(LT.sort_values("held_out_%", ascending=False).head(12).round(2))
claims.value("loco_held_out_pct", LO["loco_held_out_pct"], rtol=1e-4, unit="%",
             note="point-weighted over 153 points, 42 compounds")
claims.value("loco_in_sample_refit_pct", LO["loco_in_sample_refit_pct"], rtol=1e-4, unit="%")
claims.value("loco_in_sample_published_pct", LO["loco_in_sample_published_pct"], rtol=1e-6, unit="%")
claims.value("loco_gap_pp", LO["loco_gap_pp"], rtol=1e-3, unit="pp")
claims.value("loco_gap_se_pp", LO["loco_gap_se_pp"], rtol=1e-3, unit="pp",
             note="clustered by compound, 41 degrees of freedom")
claims.bound("loco_gap_t", LO["loco_gap_t"], min=2.0, basis="chosen",
             note="the excess is called a degradation only beyond two standard errors")
claims.exact("loco_n_compounds", LO["loco_n_compounds"])
claims.exact("loco_n_points", LO["loco_n_points"])
claims.exact("loco_fits_at_bound", LO["loco_fits_at_bound"])

# %% [markdown]
# Predicting a compound the fit never saw costs more than fitting it, by more
# than two standard errors. Set against the paper's headline average, which is
# over all points and entirely in sample, the honest statement for a *new*
# organic vapour built from these increments is the held-out figure. The cost is
# largest where a compound is the only carrier of a structural feature. In some
# refits a volume goes to its lower bound; the printed list names them.
#
# ### 5. The temperature exponent, measured from the paper's own data
#
# The paper says observed exponents *"usually lie in the range 1.6 to 1.8"* and
# fits "1.749 ± 0.013". Series in Table III with four or more temperatures give
# a log-log slope of the **observed** $D$ against $T$, with its standard error,
# and no correlation involved. A series is placed inside or outside the band
# only when its slope is more than two standard errors from the edge.

# %%
BAND = (1.6, 1.8)            # the paper's "1.6 to 1.8"
CI_WIDTH = 2 * 0.013         # full width of the paper's "1.749 ± 0.013"


def exponent_study(min_temps=4, band=BAND, T_offset=0.0, source="obs", se_scale=1.0):
    rows = []
    D = DATA.D_obs_cm2_s * ((1 + DATA.err_fsg_pct / 100) if source == "calc" else 1.0)
    for k, s in DATA.groupby(["component_a", "component_b", "ref_number"]):
        if s.temperature_K.nunique() < min_temps:
            continue
        f = stats.linregress(np.log(s.temperature_K + T_offset), np.log(D[s.index]))
        rows.append((f"{k[0]}-{k[1]}", int(k[2]), len(s), s.temperature_K.min(),
                     s.temperature_K.max(), f.slope, se_scale * f.stderr))
    E = pd.DataFrame(rows, columns=["pair", "ref", "n", "T_min", "T_max", "b_obs", "se"])
    lo, hi = band
    dist = np.minimum(np.abs(E.b_obs - lo), np.abs(E.b_obs - hi))
    sure_in = (E.b_obs - 2 * E.se >= lo) & (E.b_obs + 2 * E.se <= hi)
    sure_out = (E.b_obs + 2 * E.se < lo) | (E.b_obs - 2 * E.se > hi)
    i, j = E.b_obs.idxmin(), E.b_obs.idxmax()
    return dict(exponent_n_series=len(E), exponent_median=float(E.b_obs.median()),
                exponent_min=float(E.b_obs[i]), exponent_max=float(E.b_obs[j]),
                exponent_in_band=int(((E.b_obs >= lo) & (E.b_obs <= hi)).sum()),
                exponent_certain_in=int(sure_in.sum()), exponent_certain_out=int(sure_out.sum()),
                exponent_spread_over_ci=float((E.b_obs[j] - E.b_obs[i]) / CI_WIDTH),
                exponent_extremes_z=float((E.b_obs[j] - E.b_obs[i]) / np.hypot(E.se[i], E.se[j])),
                _margin=float(dist.min()), _margin2=float(np.abs(dist - 2 * E.se).min()), _E=E)


EX = exponent_study()
EXP = EX["_E"]
print(EXP.round(4).to_string(index=False))
print(f"\nband {BAND[0]} to {BAND[1]}: {EX['exponent_in_band']} of {EX['exponent_n_series']} slopes inside; "
      f"by more than two standard errors {EX['exponent_certain_in']} inside, "
      f"{EX['exponent_certain_out']} outside, "
      f"{EX['exponent_n_series'] - EX['exponent_certain_in'] - EX['exponent_certain_out']} undecided")
print(f"median {EX['exponent_median']:.3f}; range {EX['exponent_min']:.3f} "
      f"({EXP.pair[EXP.b_obs.idxmin()]}) to {EX['exponent_max']:.3f} ({EXP.pair[EXP.b_obs.idxmax()]}), "
      f"{EX['exponent_extremes_z']:.1f} combined standard errors apart; range over the paper's "
      f"interval width {EX['exponent_spread_over_ci']:.1f}")
claims.exact("exponent_n_series", EX["exponent_n_series"])
for k in ("exponent_median", "exponent_min", "exponent_max"):
    claims.value(k, EX[k], rtol=1e-6)
claims.exact("exponent_in_band", EX["exponent_in_band"], margin=EX["_margin"], min_margin=1e-6,
             note="slopes in [1.6, 1.8]; margin against round-off, the statistics are in "
                  "exponent_certain_in and exponent_certain_out")
for k in ("exponent_certain_in", "exponent_certain_out"):
    claims.exact(k, EX[k], margin=EX["_margin2"], min_margin=1e-6, note="two standard errors")
claims.value("exponent_spread_over_ci", EX["exponent_spread_over_ci"], rtol=1e-6)
claims.value("exponent_extremes_z", EX["exponent_extremes_z"], rtol=1e-6)

fig, ax = plt.subplots(figsize=(7, 3.6))
ax.axvspan(*BAND, color="0.88", label="paper: 1.6 to 1.8")
ax.axvline(1.749, color="C3", lw=1.6, label="fitted b = 1.749")
ax.axvline(1.5, color="C0", ls="--", lw=1.4, label="hard sphere, 3/2")
ax.errorbar(EXP.b_obs, np.arange(len(EXP)), xerr=2 * EXP.se, fmt="o", ms=4, color="C1",
            capsize=2, label="slope, two standard errors")
ax.set_xlabel("measured d ln D / d ln T"); ax.set_yticks([])
ax.set_title("Temperature exponent of the multi-temperature series of Table III")
ax.legend(fontsize=8, loc="lower right"); fig.tight_layout(); plt.show()

# %% [markdown]
# Most slopes fall in the paper's band, but most series are too short or too
# narrow in temperature to place a slope inside or outside it. The per-system
# exponents do differ, by many standard errors at the extremes, so the paper's
# "± 0.013" on $b$ is the precision of a single global exponent fitted to all
# points, not the uncertainty in the temperature dependence of any one gas pair.
#
# ### 6. What the exponent costs in a transport calculation
#
# The tube of the implementation section. The exponent is moved from 1.75 to the
# smallest and largest slopes of Section 5 and to the hard-sphere 3/2, with $D$
# pinned at the cold end (300 K) or at the mid-range temperature. Each change is
# the difference of two pymrm solves on the same grid; the quadrature route
# computes the same changes independently.

# %%
def tube_study(n=N_GRID, ladder=LADDER, nu=0, face_rule="face", bc_right=None, b=1.75,
               b_pymrm=None, T_hot=T_HOT, T_anchor=T_COLD, b_lo=None, b_hi=None, b_hs=1.5):
    b_lo = EX["exponent_min"] if b_lo is None else b_lo
    b_hi = EX["exponent_max"] if b_hi is None else b_hi
    opts = dict(nu=nu, face_rule=face_rule, bc_right=bc_right)

    def pymrm_flux(nn, T_ref, bb):
        c = tube_case(T_ref, bb, T_hot)
        if b_pymrm is not None:                  # a defect in the pymrm route only
            c["b"] = b_pymrm
        return Tube(nn, **c, **opts).solve()
    quad_flux = lambda T_ref, bb: flux_quadrature(**tube_case(T_ref, bb, T_hot))[0]
    N_q, s = quad_flux(T_MID, b), pymrm_flux(n, T_MID, b)
    st = observed_orders(lambda nn: pymrm_flux(nn, T_MID, b)["flux"], ladder)
    v = st["values"]
    out = dict(flux_pymrm_vs_quadrature=abs(s["flux"] / N_q - 1),
               flux_grid_order=float(st["orders"][-1]),
               flux_grid_order_min=float(np.min(st["orders"])),
               flux_grid_order_max=float(np.max(st["orders"])),
               flux_roundoff_margin=float(abs(v[-1] - v[-2]) / abs(v[-1])
                                          / (ladder[-1]**2 * np.finfo(float).eps)),
               flux_richardson_vs_quadrature=float(abs(st["extrapolated"] / N_q - 1)),
               flux_spread_structural=s["spread"] / abs(s["flux"]),
               flux_scale_identity=abs(Tube(n, **{**tube_case(T_MID, b, T_hot),
                                                  "D_ref": 1.043 * tube_case(T_MID, b, T_hot)["D_ref"]},
                                            **opts).solve()["flux"] / s["flux"] - 1.043),
               _orders=st["orders"])
    gaps = []
    for tag, T_a in (("cold", T_anchor), ("mid", T_MID)):
        base, q_base = pymrm_flux(n, T_a, b)["flux"], quad_flux(T_a, b)
        for key, bb in (("min", b_lo), ("max", b_hi), ("15", b_hs)):
            ch = 100 * (pymrm_flux(n, T_a, bb)["flux"] / base - 1)
            out[f"flux_change_b_{key}_anchor_{tag}"] = ch
            gaps.append(abs(ch - 100 * (quad_flux(T_a, bb) / q_base - 1)))
    out["flux_change_routes_gap_pp"] = max(gaps)
    return out


TS = tube_study()
print(f"pymrm against quadrature at {N_GRID} cells: {TS['flux_pymrm_vs_quadrature']:.2e} relative")
print(f"observed orders on {LADDER[0]} to {LADDER[-1]} cells: "
      + ", ".join(f"{p:.3f}" for p in TS["_orders"])
      + f"; finest difference over the round-off floor n^2 eps: {TS['flux_roundoff_margin']:.0f}")
print(f"Richardson limit against quadrature: {TS['flux_richardson_vs_quadrature']:.1e}; "
      f"flux spread over the faces below 1e-10 of the flux; "
      f"D scaled by 1.043 scales the flux by 1.043 to below 1e-10")
print(f"\nflux change, % (b from 1.75 to {EX['exponent_min']:.3f} / {EX['exponent_max']:.3f} / 1.5):")
for tag, T_a in (("cold", T_COLD), ("mid", T_MID)):
    print(f"  D pinned at {T_a:.1f} K: " + " / ".join(
        f"{TS[f'flux_change_b_{k}_anchor_{tag}']:+.2f}" for k in ("min", "max", "15")))
print(f"pymrm and quadrature changes agree to {TS['flux_change_routes_gap_pp']:.0e} pp")
claims.value("flux_pymrm_vs_quadrature", TS["flux_pymrm_vs_quadrature"], rtol=1e-3,
             note="relative, 400 cells, D pinned mid-range")
claims.bound("flux_grid_order", TS["flux_grid_order"], min=1.9, max=2.1, basis="theory",
             note="finest triple 200-400-800")
claims.bound("flux_grid_order_min", TS["flux_grid_order_min"], min=1.8, max=2.1, basis="theory",
             note="second order approached from below; the coarsest triple is lowest")
claims.bound("flux_grid_order_max", TS["flux_grid_order_max"], min=1.9, max=2.1, basis="theory")
claims.bound("flux_roundoff_margin", TS["flux_roundoff_margin"], min=100, basis="theory",
             note="finest successive difference over n^2 eps")
claims.bound("flux_richardson_vs_quadrature", TS["flux_richardson_vs_quadrature"], max=1e-8,
             basis="tolerance")
claims.bound("flux_spread_structural", TS["flux_spread_structural"], max=1e-10,
             basis="identity", structural=True,
             note="a conservative finite-volume divergence: one flux through every face")
claims.bound("flux_scale_identity", TS["flux_scale_identity"], max=1e-10, basis="identity",
             structural=True, note="the problem is linear in Gamma")
for tag in ("cold", "mid"):
    for k in ("min", "max", "15"):
        claims.value(f"flux_change_b_{k}_anchor_{tag}", TS[f"flux_change_b_{k}_anchor_{tag}"],
                     rtol=1e-4, unit="%",
                     note=f"D pinned at {'300 K' if tag == 'cold' else 'sqrt(300 x 600) K'}")
claims.bound("flux_change_routes_gap_pp", TS["flux_change_routes_gap_pp"], max=1e-4,
             basis="tolerance", note="discretisation at 400 cells")

fig, ax = plt.subplots(figsize=(6.6, 3.2))
for b_, ls in ((EX["exponent_min"], "--"), (1.75, "-"), (EX["exponent_max"], ":")):
    t = Tube(200, **tube_case(T_COLD, b_))
    ax.plot(t.z_c, t.solve()["x"], ls, lw=1.8, label=f"b = {b_:.2f}")
ax.set_xlabel("z (cm)"); ax.set_ylabel("benzene mole fraction")
ax.set_title("Benzene through air, 300 K to 600 K, D pinned at 300 K")
ax.legend(fontsize=9); fig.tight_layout(); plt.show()

# %% [markdown]
# **The exponent is nearly free if $D$ is pinned in the middle of the
# temperature range, and costs several per cent if it is pinned at one end.**
# The paper reports the exponent's confidence interval, not what it does to an
# answer. The asymmetry is a property of the problem: the integrand $1/\Gamma$ is
# pinned at whatever temperature $D$ is pinned at.
#
# A scale error in $D$, by contrast, transfers one for one: the problem is linear
# in $\Gamma$, so a relative error in $D$ is the same relative error in the flux
# (the cell above checks it for a factor 1.043). The fit error and the held-out
# error of Sections 2 and 4 carry over to the flux unchanged.

# %% [markdown]
# ## Validation
#
# Four checks, each followed by what it cannot see.
#
# **(a) The correlation against the printed error column** (Section 1): one
# exact equation constrains three transcribed numbers per row. *Cannot detect:*
# an error in the reference number, a molar mass wrong in the same way in the
# paper's arithmetic and here, or the identity of Table III's other columns.
#
# **(b) Route B, the $T^{1.75}$ identity within a series**: only the printed
# $T$, $D_\text{obs}$ and error, sharing no code with (a). *Cannot detect:* an
# error common to every temperature of one series, or any error in a system
# measured at a single temperature.
#
# **(c) The printed aggregates** (Section 2). *Cannot detect:* compensating
# errors in single rows; an average cannot see the two irreproducible rows,
# which is why (a) and (b) exist.
#
# **(d) pymrm against adaptive quadrature** (Section 6), second order with a
# Richardson limit on the quadrature value. *Cannot detect:* a wrong exponent, a
# wrong diffusivity or a wrong molar mass, because both routes read the same
# $\Gamma(T)$; the break table shows it with a row that changes $b$ in both. The
# cell-to-cell spread of the flux is an identity of the conservative scheme and
# is registered as structural.
#
# A break row measures sensitivity, never correctness: every row below shows a
# claim failing when something is broken, none shows that the unbroken value is
# right. A bias shared by the paper's arithmetic and this reimplementation would
# pass every row; that is why the reproduction is called a check on the
# transcription, not on the physics.
#
# ### The break table
#
# Each row reruns one or more studies with one defect, injected through a call
# argument or a published constant. A study that cannot be evaluated under a
# defect records `nan` for its claims, which fails them. Each study run with no
# defect must return the registered values; the cell checks that first.

# %%
STUDIES = dict(data=data_study, repro=repro_study, inference=inference_study, routeB=routeB_study,
               agg=aggregates_study, fit=fit_study, loco=loco_study,
               exp=exponent_study, tube=tube_study)
BASES = {key: f() for key, f in STUDIES.items()}
for key, base in BASES.items():
    for name, x in base.items():
        if name not in claims.metrics:
            continue
        x0 = claims.metrics[name]["value"]
        if claims.metrics[name]["kind"] == "exact":
            assert x == x0, (key, name, x, x0)
        else:
            assert abs(x - x0) <= 1e-9 * abs(x0) + 1e-15, (key, name, x, x0)
STUDY_CLAIMS = {key: [k for k in base if k in claims.metrics and
                      not claims.metrics[k].get("structural")] for key, base in BASES.items()}

with_par = lambda **kw: dict(PAR, **kw)
ALT_SET = with_par(v_atom={**PAR["v_atom"], "C": 15.9, "H": 2.31, "O": 6.11, "N": 4.54})
SWAP = with_par(v_simple={**PAR["v_simple"], "N2": PAR["v_simple"]["O2"], "O2": PAR["v_simple"]["N2"]})
C16 = with_par(v_atom={**PAR["v_atom"], "C": 16.0})
B170 = with_par(b=1.70)
D2V = with_par(v_simple={**PAR["v_simple"], "D2": 7.6})
RING0 = np.r_[FIT.lower[:6], 0.0, FIT.lower[7:]]
NEUMANN = {"a": 1.0, "b": 0.0, "d": 0.0}         # -Gamma dx/dn = 0 at z = L: no flux out
ROWS = [
    ("Table III page 6 not read", dict(data=dict(drop_page=int(DATA.pdf_page.min())))),
    ("Table I: C 16.5 read as 16.0", dict(repro=dict(par=C16), inference=dict(par=C16),
                                          fit=dict(par=C16), loco=dict(par=C16))),
    ("Table I: Sigma v(D2) 6.7 read as 7.6", dict(repro=dict(par=D2V), inference=dict(par=D2V))),
    ("the alternative increments C 15.9, H 2.31, O 6.11, N 4.54 (in no document on disk)",
     dict(repro=dict(par=ALT_SET))),
    ("Sigma v(N2) and Sigma v(O2) swapped", dict(repro=dict(par=SWAP))),
    ("b 1.75 read as 1.749, the paper's unrounded fit", dict(repro=dict(par=with_par(b=1.749)))),
    ("eq. (4) prefactor 1.00e-3 read as 1.01e-3", dict(repro=dict(par=with_par(C=1.01e-3)))),
    ("ring term on piperidine, tetrahydrothiophene and ethylene oxide",
     dict(repro=dict(par=with_par(ring_extra=("PIPERIDINE", "TETRAHYDROTHIOPHENE", "ETHYLENEOXIDE"))))),
    ("M(D2) 4.00 read as 4.028",
     dict(repro=dict(par=with_par(mw_simple={**PAR["mw_simple"], "D2": 4.028})))),
    ("temperature column read 273.15 K high", dict(routeB=dict(T_offset=273.15),
                                                   exp=dict(T_offset=273.15))),
    ("identity exponent 1.75 read as 1.70 (route B only)", dict(routeB=dict(b_id=1.70))),
    ("series formed without the reference number", dict(routeB=dict(by_ref=False))),
    ("printed error column read 3 % high", dict(agg=dict(err_scale=1.03))),
    ("errors relative to D_calc in place of D_obs", dict(agg=dict(err_basis="calc"))),
    ("parameter count 28 read as 27 (b not counted)", dict(agg=dict(n_par=N_PAR - 1))),
    ("tail threshold 10 % read as 10.5 %", dict(agg=dict(tail=10.5))),
    ("Table II category by the first component only", dict(agg=dict(member="first"))),
    ("printed errors read as rounded to one decimal", dict(agg=dict(err_digits=1))),
    ("b 1.75 read as 1.70 in eq. (4) and every refit", dict(fit=dict(par=B170), loco=dict(par=B170))),
    ("absolute residuals in place of eq. (3)'s relative ones", dict(fit=dict(objective="absolute"))),
    ("refit tolerance 1e-14 raised to 1e-3", dict(fit=dict(tol=1e-3))),
    ("ring increment bounded below by 0", dict(fit=dict(lower=RING0))),
    ("held-out compound left in its own fit (in sample by construction)", dict(loco=dict(leak=True))),
    ("compounds with a single point skipped", dict(loco=dict(min_points=2))),
    ("series with three temperatures admitted", dict(exp=dict(min_temps=3))),
    ("band 1.6 to 1.8 read as 1.65 to 1.85", dict(exp=dict(band=(1.65, 1.85)))),
    ("slopes of the printed calculated D in place of the observed (1.75 by construction)",
     dict(exp=dict(source="calc"))),
    ("slope standard errors halved", dict(exp=dict(se_scale=0.5))),
    ("tube grid 400 cells read as 50", dict(tube=dict(n=50))),
    ("construct_div nu 0 read as 1 (cylindrical)", dict(tube=dict(nu=1))),
    ("outlet Dirichlet replaced by zero flux", dict(tube=dict(bc_right=NEUMANN))),
    ("face Gamma from the low-z cell (first order)", dict(tube=dict(face_rule="left_cell"))),
    ("b 1.75 read as 1.50 in the pymrm route only (one route only)", dict(tube=dict(b_pymrm=1.5))),
    ("b 1.75 read as 1.50 in both tube routes", dict(tube=dict(b=1.5))),
    ("hot end 600 K read as 450 K", dict(tube=dict(T_hot=450.0))),
    ("cold anchor 300 K read as 350 K", dict(tube=dict(T_anchor=350.0))),
    ("exponent extremes replaced by the paper's band 1.6 and 1.8",
     dict(tube=dict(b_lo=1.6, b_hi=1.8))),
    ("order ladder 6400 to 204800 cells", dict(tube=dict(ladder=(6400, 12800, 25600, 51200,
                                                                 102400, 204800)))),
]
claims.define_rows([row for row, _ in ROWS])

BREAK, REFUSED, OUTS = [], [], {}
for row, studies in ROWS:
    for key, kw in studies.items():
        try:
            out = OUTS[row, key] = STUDIES[key](**kw)
        except Exception as exc:                      # the study has no result under this defect
            REFUSED.append((row, key, f"{type(exc).__name__}: {exc}"[:120]))
            out = {}
        for name in STUDY_CLAIMS[key]:
            x = out.get(name, np.nan)
            if isinstance(x, (float, np.floating)) and not np.isfinite(x):
                x = np.nan
            BREAK.append((row, name, claims.breaks(row, name, x)))

bt = pd.DataFrame(BREAK, columns=["row", "claim", "caught"])
caught = bt[bt.caught].groupby("row", sort=False).claim.apply(list)
printed = list(bt.row.unique())
print(f"{len(ROWS)} rows defined, {len(printed)} printed; claims each row fails:")
for row in printed:
    names = caught.get(row, [])
    print(f"  {row}:")
    print("      " + (", ".join(names) if names else "none"))
assert len(printed) == len(ROWS)
print("\nrows whose study has no result (nan for its claims):")
for row, key, why in REFUSED:
    print(f"  {row} [{key}]: {why}")
missing = [n for n in claims.metrics
           if not claims.metrics[n].get("structural") and not bt[(bt.claim == n) & bt.caught].shape[0]]
print(f"\nclaims no row fails: {missing if missing else 'none'}")
both = OUTS["b 1.75 read as 1.50 in both tube routes", "tube"]["flux_pymrm_vs_quadrature"]
print(f"pymrm against quadrature with b = 1.50 in both routes: {both:.2e} "
      f"(at b = 1.75: {claims['flux_pymrm_vs_quadrature']:.2e})")

# %% [markdown]
# The row that changes $b$ in both tube routes leaves the pymrm and quadrature
# routes agreeing at the discretisation level (the last line above): that check
# reads the same $\Gamma$ twice and cannot see a wrong exponent. Its value moves
# only because the discretisation error depends on $b$; the flux changes, which
# depend on the base exponent, are what fail.

# %% [markdown]
# ## What pymrm adds
#
# Not much to the correlation, and something specific to using it.
#
# **What pymrm does not add.** Eq. (4) is four arithmetic operations. The
# reproduction of the printed tables, the leave-one-compound-out
# cross-validation and the measured temperature exponents are pandas, numpy and
# scipy; pymrm is not involved.
#
# **What it does add.** The question a correlation cannot answer is what its
# error does to an answer, and that needs a transport solve: one linear sparse
# solve of a variable-coefficient diffusion operator assembled from
# `construct_grad`, `construct_div` and `construct_coefficient_matrix`, with the
# flux read at the boundary by `compute_boundary_values`. On this page it gives
# one result the paper does not state: the cost of the exponent depends on where
# $D$ is pinned, a few tenths of a per cent mid-range against several per cent
# at the cold end (Section 6). A scale error in $D$ needs no solve at all: it
# carries over to the flux one for one. The same three operator calls serve
# whether $\Gamma$ comes from FSG, Chapman-Enskog or a measurement.

# %% [markdown]
# ## Reuse
#
# **The correlation.** `D_fsg(A, B, T, PAR, p)` and `sigma_v(name, PAR)` in
# `model.py`, with `PAR = fsg_params(volumes)` from
# `data/fuller-1966-table1-diffusion-volumes.csv`. Use the 1966 volumes as
# printed, with four cautions:
#
# 1. **Quote the paper's average error as a fit, not an accuracy.** For an
#    organic vapour outside the original set, the leave-one-out estimate of
#    Section 4 is the better guide, and a compound carrying a structural feature
#    no other compound carries does worse (the largest excesses are printed
#    there).
# 2. **For a simple molecule not in Table I there is no method**, and no
#    held-out estimate exists for those in it, because each has its own fitted
#    parameter. Do not transfer the organics' held-out figure to them.
# 3. **The 1966 increments are not the ones in most handbooks.** The break table
#    shows that the widely reproduced set destroys this page's reproduction.
#    Whichever set you use, say which, and cite the document you read it from.
# 4. **The 1966 erratum, `10.1021/ie50680a601`, has not been consulted.** Anyone
#    with access should check the volumes against it.
#
# **The dataset.** `data/fuller-1966-table3-diffusion-data.csv` holds the
# measured binary diffusivities with the correlation's own error beside each.
# For any estimate other than FSG it is an out-of-sample benchmark, because only
# FSG was fitted to it. The two named rows of Section 1 carry a printed FSG error
# eq. (4) cannot reproduce; their measurements are unaffected and both are kept.
#
# **The transport model.** `model.Tube` is a variable-coefficient 1-D diffusion
# solve with Dirichlet ends; replace `gamma` for any $\Gamma(z)$.
# `flux_quadrature` is its independent reference for a smooth coefficient.
#
# **Related pages.** `A4.6` (Chapman-Enskog) is the first-principles route to
# the same quantity and should carry the comparison of the two on one measured
# axis. `A4.2` (Maxwell-Stefan against Fick), `A4.3` (dusty gas), `A4.4`
# (Knudsen and Bosanquet) and `A4.9` (Duncan-Toor) all use a binary diffusivity
# of the kind this page estimates.

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
