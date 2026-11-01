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
# title: "Aris's shape factor, and the shapes he could not test"
# description: "Dividing volume by external surface area makes the effectiveness factor almost shape-independent. Aris measured 'almost' on the three shapes that have closed-form answers. Eleven more shapes widen the band by about a quarter, and the sphere, which his table makes look like the extreme case, is not the floor."
# categories: [sec:B, struct:S3, struct:S13, tier:T0, data:tier6, phase:gas-solid]
# date: 2026-07-31
# ---
#
# # Aris's shape factor, and the shapes he could not test
#
# **Catalog ID:** `B1.2` · **Structures:** `S3` (1D steady BVP), `S13`
# (multi-dimensional / non-standard geometry) · **Tier:** T0
#
# Thiele solved the catalyst pellet for a slab and for a sphere and got two
# different curves. Aris (1957) showed that they are nearly the same curve drawn
# against the wrong variable: use $v_p/s_x$, the ratio of particle volume to
# external surface area, as the characteristic length and they very nearly
# collapse.
#
# "Very nearly" is the entire content of the paper, and Aris quantified it on the
# only three shapes whose answer he could write down: flat plate, long cylinder,
# sphere. Of the general case he wrote that estimating the divergence "would be an
# excessively difficult task, equivalent to solving the problem completely."
#
# This page evaluates the shapes his table leaves out: the spherical shell he gave
# in closed form, the finite cylinder he left as an unevaluated series (from a long
# rod to a flat disc), a square rod and a cube. One pymrm class solves all of
# them, and each solve is checked against an independent series or closed form.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.11.01.4"  # replaced by the release tag at publication
PYMRM = ">=3.3.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "B1.2-aris-shape-factor"  # equals the page directory

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
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from scipy.optimize import minimize_scalar
from scipy.special import iv, jn_zeros
from pymrm.checks import find_roots, observed_orders
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import (EXPOSED, Pellet, slab, long_cylinder, sphere,  # this page's model.py
                   spherical_shell, finite_cylinder, square_rod, cube)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])
BREAK = []
FLUX_GAPS = {}      # |eta from the surface flux - eta as volume average|, per solve


def brk(row, metric, broken):
    """One injected defect: ``broken`` is what claim ``metric`` becomes with it."""
    caught = claims.breaks(row, metric, broken)
    BREAK.append((row, metric, claims[metric], broken, caught))
    return caught


def rel(a, b):
    return a / b - 1


def eta_both(pel, lam, key=None):
    """Solve once; eta as the volume average, and the surface-flux gap recorded."""
    chi = pel.solve(lam)
    v = pel.eta(lam, chi)
    if key is not None:
        FLUX_GAPS[key] = abs(pel.eta_flux(lam, chi) - v)
    return v, chi


# %% [markdown]
# ## Background
#
# A porous pellet with a first-order reaction in it obeys
#
# $$\nabla^2 c = \lambda^2 c \quad\text{in } R, \qquad c = c_s \ \text{ on } S,
# \qquad \lambda^2 = \frac{k\sigma}{D},$$
#
# and the effectiveness factor is the volume average of the normalised
# concentration, $\eta = \langle c/c_s\rangle$. Everything about the pellet enters
# through the shape of $R$ and the single group $\lambda$, which has units of
# 1/length. To build a dimensionless modulus you must choose a length, and the
# choice is arbitrary until you demand that different shapes give the same
# answer.
#
# Aris's argument is short. For all three classical shapes the strong-diffusion
# asymptote is $\eta \to m/(\lambda a)$ with $m = 1, 2, 3$ for plate, cylinder and
# sphere, and in each case $a/m$ is exactly $v_p/s_x$. So the choice
# $a = v_p/s_x$ makes the three curves coincide *exactly* at both ends. What is
# left to check is the middle.
#
# **Why this page is not [`B1.1`](../B1.1-thiele-weisz-hicks/).** That page plots
# the collapse and reports the residual spread between Thiele's three geometries,
# three curves that all have closed forms. The question here is the one Aris
# declined: what happens for shapes that have no closed form. Corners, edges,
# hollow interiors, finite aspect ratios.

# %% [markdown]
# ## The published model
#
# Aris's generalised modulus, and the three closed forms it produces (his
# eqs. 15-17):
#
# $$\Lambda \;=\; \lambda\,\frac{v_p}{s_x} \;=\; \left(\frac{v_p^2\,k\sigma}{s_x^2\,D}\right)^{1/2}$$
#
# | shape | characteristic length $a$ | $v_p/s_x$ | $\eta(\Lambda)$ | eq. |
# |---|---|---|---|---|
# | flat plate | half-thickness | $a$ | $\tanh\Lambda/\Lambda$ | 15 |
# | cylinder, ends sealed | radius | $a/2$ | $I_1(2\Lambda)/\bigl(\Lambda\,I_0(2\Lambda)\bigr)$ | 16 |
# | sphere | radius | $a/3$ | $(3\Lambda\coth 3\Lambda - 1)/(3\Lambda^2)$ | 17 |
#
# Two further results in the paper carry more weight than the table.
#
# **The general asymptote, his eq. 18.** For *any* shape, as $\lambda\to\infty$
# only a thin surface layer reacts, $\chi \simeq e^{-\lambda\xi}$ normal to the
# surface, and
#
# $$\eta \;=\; \frac{1}{\lambda^2 v_p}\oint_S \frac{\partial\chi}{\partial n}\,\mathrm{d}s
# \;\longrightarrow\; \frac{s_x}{\lambda\,v_p} \;=\; \frac{1}{\Lambda}.$$
#
# That holds for every shape with a piecewise-smooth boundary. It is the reason
# the collapse works at all, and it is testable on shapes Aris never wrote down.
#
# **The spherical shell.** Active catalyst between radii $pa$ and $a$, with only
# the outer surface exposed, so $s_x = 4\pi a^2$ and
# $v_p = 4\pi a^3(1-p^3)/3$. Setting $(1-p^3)\Lambda' = 3(1-p)\Lambda$, so that
# $\Lambda' = \lambda a(1-p)$ is built on the shell *thickness*,
#
# $$\eta = \frac{3(1-p)}{(1-p^3)\Lambda'^2}\cdot
# \frac{(1-p)^2\bigl(\Lambda' - \tanh\Lambda'\bigr) + p\,\Lambda'^2\tanh\Lambda'}
#      {(1-p)\tanh\Lambda' + p\,\Lambda'} .$$
#
# A one-parameter family joining the sphere ($p=0$) to the flat plate ($p=1$), and
# the paper says so explicitly, which makes it a free check on the transcription.
#
# **The finite cylinder, his case (iv)**: radius $a$, flat ends *porous* rather
# than sealed, given only as a double Fourier-Bessel series,
#
# $$\eta = 1 - \frac{32}{\pi^2}\sum_{m=1}^{\infty}\sum_{n=0}^{\infty}
# \frac{1}{(2n+1)^2 j_m^2}\cdot
# \frac{\lambda^2a^2}{\lambda^2a^2 + j_m^2 + (2n+1)^2p^2\pi^2},$$
#
# with $j_m$ the zeros of $J_0$. He never evaluates it. **And $p$ is ambiguous**:
# the nomenclature defines it as radius:length, the running text says "length
# $= 2a/p$", and the two differ by the factor of 2 to which the axial eigenvalue
# $(2n+1)\pi/\mathrm{length}$ is most sensitive. The results section settles it
# with a two-dimensional solve instead of guessing.
#
# **Mixtures, his eqs. 19-20.** For a bed of mixed sizes,
# $\eta = \sum_i f_i\,\eta(\Lambda_i)$, and only in the strong-diffusion limit may
# a mean size be taken first, $\bar\Lambda = (\sum_i f_i/\Lambda_i)^{-1}$.

# %% [markdown]
# ## Parameters and assumptions
#
# **Assumptions**, all Aris's: steady state; first-order irreversible kinetics;
# constant effective diffusivity; isothermal; no external film resistance, so
# $c = c_s$ everywhere on the exposed surface; and the pellet small enough that
# $c_s$ is uniform around it. His section 5 adds a film resistance in series,
# $E^{-1} = \eta^{-1} + \eta_f^{-1}$; that is not needed here and is not
# reproduced.
#
# **There are no physical parameters.** The problem is $\nabla^2\chi =
# \lambda^2\chi$ with $\chi = 1$ on the exposed boundary, and the only inputs are
# the shape and $\Lambda$. Every number on this page is dimensionless.
#
# **How the paper was obtained, and what could not be.** Aris (1957) is not open
# access. It was retrieved through the Elsevier article retrieval API under an
# institutional subscription (PII `0009250957850283`). The 1957 scan's text layer
# is the publisher's OCR of a mid-dot-decimal typesetting and is unusable for
# numbers: it renders Table 1's 0.698 as `898` and 0.432 as `482`. Chemical
# Engineering Science **reprinted the paper verbatim in 1995** for its Golden
# Jubilee (50(24) 3899-3903, PII `0009250996818197`), and that re-typesetting OCRs
# cleanly. Both were retrieved, and every transcribed number is one on which the
# two agree, or where the disagreement is stated.
#
# The entitlement returns a **one-page preview PDF only** (page 262 of the 1957
# printing, page 3899 of the 1995 reprint), and Table 1 is on page 265 / 3901.
# The article-object endpoint holds no page images for either PII. So the
# gallery's usual discipline of reading numbers off a 600 dpi render is not
# available for Table 1. What stands in for it is the paper's own algebra, below.

# %% [markdown]
# ## The data
#
# Three small tables of the author's own computed values. **They are not
# measurements**: Aris (1957) reports no experiments, and this is a provenance
# tier 6 page. Reproducing them tests that his equations have been transcribed
# correctly and that the numerics here agree with his; it is not a test against
# nature.

# %%
tab1 = load_data("aris1957-table1.csv", page=PAGE)
shell_tab = load_data("aris1957-spherical-shell.csv", page=PAGE)
mix = load_data("aris1957-mixture-example.csv", page=PAGE)
for _f in ("aris1957-table1.csv", "aris1957-spherical-shell.csv",
           "aris1957-mixture-example.csv"):
    print(cite_data(load_meta(_f, page=PAGE)))
print(f"{len(tab1)} Table 1 entries, {tab1['shape'].nunique()} shapes, "
      f"Lambda {tab1['Lambda'].min()} to {tab1['Lambda'].max()}")
print(f"{len(shell_tab)} spherical-shell entries at Lambda = {shell_tab['Lambda'].unique()}")
print(f"mixture example: {int((mix['quantity'] == 'Lambda').sum())} sizes, "
      f"{int(mix['quantity'].str.startswith(('mean', 'eta')).sum())} printed results")

# %% [markdown]
# ## PyMRM implementation
#
# **The references first.** Eqs. 15-17, the shell expression and the two
# eigenfunction series are written out here, independently of the solver. The
# finite-cylinder series takes the aspect ratio $\tau$ (half-length / radius) as a
# separate argument from $p$, so that the two readings of $p$ can be tested
# against each other instead of one being assumed. Both series are truncated;
# their truncation error is measured in the validation section by doubling the
# number of terms.

# %%
def eta_plate(Lam):
    """Aris eq. 15, flat plate."""
    return np.tanh(Lam) / Lam


def eta_cylinder(Lam):
    """Aris eq. 16, long cylinder with sealed ends."""
    return iv(1, 2 * Lam) / (Lam * iv(0, 2 * Lam))


def eta_sphere(Lam):
    """Aris eq. 17, sphere."""
    return (3 * Lam / np.tanh(3 * Lam) - 1) / (3 * Lam**2)


def eta_shell(Lam, p, on_thickness=True):
    """Aris's spherical-shell expression; p = inner:outer radius ratio.

    on_thickness=False builds Lambda' on the outer radius instead of the shell
    thickness; it exists only for the break table.
    """
    if p <= 0.0:
        return eta_sphere(Lam)
    if p >= 1.0:
        return eta_plate(Lam)
    Lp = 3.0 * (1.0 - p if on_thickness else 1.0) * Lam / (1.0 - p**3)
    pre = 3.0 * (1.0 - p) / ((1.0 - p**3) * Lp**2)
    num = (1 - p) ** 2 * (Lp - np.tanh(Lp)) + p * Lp**2 * np.tanh(Lp)
    return pre * num / ((1 - p) * np.tanh(Lp) + p * Lp)


def eta_fincyl_series(Lam, p, tau, M=400, N=3000):
    """Aris case (iv), exactly as printed, truncated at M radial and N axial terms.

    p enters only through the axial term, as printed. tau = half-length/radius
    fixes the geometry, and with it v_p/s_x = a tau/(2 tau + 1). The paper's two
    definitions of p correspond to tau = 1/(2p) (nomenclature) and tau = 1/p
    (running text).
    """
    lam_a = Lam * (2.0 * tau + 1.0) / tau
    jm = jn_zeros(0, M)[:, None]
    n = np.arange(N)[None, :]
    weight = 1.0 / (jm**2 * (2 * n + 1) ** 2)          # sums to pi^2/32 over all m, n
    axial = (2 * n + 1) ** 2 * p**2 * np.pi**2
    return 1.0 - 32.0 / np.pi**2 * float(np.sum(
        weight * lam_a**2 / (lam_a**2 + jm**2 + axial)))


J0_ZEROS = jn_zeros(0, 12800)


def eta_fincyl_axial(Lam, p, tau, M=6400):
    """Case (iv) with the axial sum done in closed form, M radial terms.

    Uses sum_n 1/((2n+1)^2 (A + (2n+1)^2 c^2)) = (pi^2/8 - pi tanh(pi a/2)/(4a))/A
    with a = sqrt(A)/c, c = p pi. Needed where lambda a is large: the double
    series then needs thousands of radial terms.
    """
    lam_a = Lam * (2.0 * tau + 1.0) / tau
    j = J0_ZEROS[:M]
    A = lam_a**2 + j**2
    a = np.sqrt(A) / (p * np.pi)
    S = (np.pi**2 / 8 - np.pi * np.tanh(np.pi * a / 2) / (4 * a)) / A
    return 1.0 - 32.0 / np.pi**2 * float(np.sum(lam_a**2 / j**2 * S))


def eta_box_series(Lam, half, va, K=200):
    """Rectangular rod (2 half-widths) or box (3), by eigenfunction expansion.

    half : half-widths in units of a;  va : (v_p/s_x)/a for that box;
    K terms per direction. The 3-D case is summed one index at a time.
    """
    lam = Lam / va
    m = 2 * np.arange(K) + 1
    w = 8.0 / (m**2 * np.pi**2)                       # sums to 1 over all m
    b = [(m * np.pi / (2 * h)) ** 2 for h in half]
    if len(half) == 2:
        mu = b[0][:, None] + b[1][None, :]
        return 1.0 - float(np.sum(np.outer(w, w) * lam**2 / (lam**2 + mu)))
    ww = np.outer(w, w)
    bb = b[1][:, None] + b[2][None, :]
    total = sum(w[i] * float(np.sum(ww * lam**2 / (lam**2 + b[0][i] + bb)))
                for i in range(K))
    return 1.0 - total


# %% [markdown]
# **The solver.** One class, `Pellet` in this page's `model.py`, covers every
# shape. The equation is linear, so there is no Newton loop: one sparse direct
# solve of $(A_\mathrm{diff} + \lambda^2 A_\mathrm{react})\,\chi = -g$ per
# $\lambda$, with the diffusion operator $A_\mathrm{diff}$ (from `construct_grad`
# and `construct_div`), its boundary vector $g$, and the unit reaction operator
# $A_\mathrm{react}$ (from `construct_coefficient_matrix`) all assembled once in
# `__init__`.
#
# Axis 0 may carry a geometry factor `nu` (`construct_div(nu=1)` cylindrical,
# `nu=2` spherical); every other axis is Cartesian. That single switch produces
# the slab, the long cylinder, the sphere, the spherical shell (`nu=2` on a domain
# starting at $r = pa$ instead of 0), the axisymmetric finite cylinder (`nu=1`
# radial plus Cartesian axial), the square rod and the cube. Each axis runs from a
# symmetry or sealed face to an exposed face:
#
# ```python
# # a dchi/dn + b chi = d with the OUTWARD normal n
# SEALED  = {"a": 1.0, "b": 0.0, "d": 0.0}   # dchi/dn = 0: symmetry plane, or a sealed face
# EXPOSED = {"a": 0.0, "b": 1.0, "d": 1.0}   # chi = 1: face exposed to the bulk
# ```
#
# $\eta$ has two routes in the class. `eta` is the volume average of $\chi$ with
# the cell volumes of the same $x^{\nu}$ measure `construct_div` uses (Aris's
# eq. 13). `eta_flux` integrates the face gradient over the exposed faces, read
# with `compute_boundary_values`, and divides by $\lambda^2 v_p$ (the first form of
# eq. 18). For a conservative finite-volume scheme these are **the same sum**
# rearranged: the divergence telescopes to the boundary flux, so the two differ
# only by the residual of the linear solve. They are not independent, and the
# page uses the volume average throughout; the gap between them is registered
# as a check on the solve and on how the boundary is read. Every comparison that
# tests the discretisation is against the references above.
#
# The grid may be graded toward the exposed faces with `q > 1` (faces at
# $1-(1-u)^q$ on uniform $u$, a nested map). It is used where the reacting layer,
# $O(1/\lambda)$ thick, is thin: the sweep over $\Lambda$ and the strong-diffusion
# asymptote.

# %% [markdown]
# ## Results
#
# ### Aris's three shapes, and his Table 1

# %%
CLOSED = {"plate": eta_plate, "cylinder": eta_cylinder, "sphere": eta_sphere}


def table1_diffs(closed=CLOSED):
    return np.array([r.eta_printed - closed[r.shape](r.Lambda) for r in tab1.itertuples()])


def lambda_needed(printed, fn=eta_cylinder):
    """Every Lambda at which fn reaches the printed value (scan and refine)."""
    return find_roots(lambda x: fn(x) - printed, 0.5, 40.0, n_scan=200, log=True)


d1 = table1_diffs()
HALF_UNIT = 5e-4    # half a unit of the printed third decimal
FLAG = 2e-3         # four half-units: what an arithmetic slip, not a rounding, exceeds
print("Aris Table 1 against his own eqs. 15-17")
print(f"{'shape':10s}{'Lambda':>8}{'printed':>10}{'eq.15-17':>11}{'diff':>10}")
for r, d in zip(tab1.itertuples(), d1):
    mark = "   <-- slip" if abs(d) > FLAG else ("   (last digit)" if abs(d) > HALF_UNIT else "")
    print(f"{r.shape:10s}{r.Lambda:8.1f}{r.eta_printed:10.3f}{r.eta_printed - d:11.5f}{d:+10.5f}{mark}")

N_SLIP = claims.count("table1_entries_off_more_than_0p002", np.abs(d1), FLAG, min_margin=0.1)
N_ROUND = claims.count("table1_entries_not_rounding_to_printed", np.abs(d1), HALF_UNIT,
                       min_margin=0.01)
ok = np.abs(d1) <= FLAG
claims.value("table1_worst_19_entries", np.abs(d1[ok]).max())
_ps5 = (tab1["shape"] != "cylinder") & (tab1["Lambda"] >= 5)
claims.value("table1_plate_sphere_worst_Lambda_5_10", np.abs(d1[_ps5.to_numpy()]).max())

CYL_SLIP = {}
for L in (5.0, 10.0):
    pr = float(tab1[(tab1["shape"] == "cylinder") & (tab1["Lambda"] == L)]["eta_printed"].iloc[0])
    roots = lambda_needed(pr)
    CYL_SLIP[L] = (pr, rel(pr, eta_cylinder(L)), roots["roots"][0] / L, roots["message"])
    claims.value(f"cyl_printed_high_pct_Lambda{L:g}", 100 * CYL_SLIP[L][1])
    claims.value(f"cyl_printed_needs_Lambda_factor_{L:g}", CYL_SLIP[L][2])
print(f"\n{len(d1) - N_SLIP} of {len(d1)} entries agree to {claims['table1_worst_19_entries']:.5f}; "
      f"{len(d1) - N_ROUND} of them round to the printed digits.")
print(f"plate and sphere at Lambda = 5 and 10: worst error "
      f"{claims['table1_plate_sphere_worst_Lambda_5_10']:.5f}")
print("the Lambda eq. 16 would need to give the printed cylinder values (find_roots):")
for L, (pr, hi, fac, msg) in CYL_SLIP.items():
    print(f"  printed {pr:.3f} at Lambda = {L:4.1f}: {hi:+.1%} above eq. 16; "
          f"needs {fac:.3f} x Lambda ({msg})")

# %% [markdown]
# Nineteen of the twenty-one entries come back to within one unit of the third
# decimal, and seventeen of them round to the printed digits exactly. The two that
# do not come back at all are both in the **cylinder** column, at the two largest
# moduli, and both are printed *high*.
#
# Three things narrow what that can be:
#
# - the **plate and sphere columns are right at those same two moduli**, so
#   nothing global has gone wrong with the last two rows;
# - the **first five rows of the cylinder column are right**, so the column really
#   is eq. 16 evaluated at the stated $\Lambda$, not some other quantity;
# - **no single rescaling of $\Lambda$ rescues them.** If the column had used a
#   differently defined modulus (a diameter for a radius, say) one factor would
#   fit both rows. The printed values need two different ones, each the only root
#   in the scanned range.
#
# The same digits appear in the 1957 printing and in the 1995 reprint, so this is
# not an OCR artefact either. What is left is two arithmetic slips at the end of
# one column.
#
# **But that is an inference, not a reading.** The two cells have never been seen
# on a page image, only through two text layers that agree; the entitlement stops
# at page 262 and no page-image route exists. A maintainer holding the full PDF
# can close this, and until then it stays open. Reported, not repaired: nothing
# on this page depends on the two cells.

# %% [markdown]
# ### The spherical shell: sphere to plate in one parameter

# %%
p_print = shell_tab["p"].to_numpy()
eta_print = shell_tab["eta_printed"].to_numpy()


def shell_formula_vs_printed(printed=eta_print, **kw):
    return float(np.abs(np.array([eta_shell(1.0, p, **kw) for p in p_print]) - printed).max())


def shell_pymrm(p, bc_low=None, n=400):
    if p == 0.0:
        pel, va = sphere(n)
    elif p == 1.0:
        pel, va = slab(n)
    else:
        pel, va = spherical_shell(p, n) if bc_low is None else spherical_shell(p, n, bc_low=bc_low)
    return eta_both(pel, 1.0 / va, key=None if bc_low is not None else f"shell p={p}")[0]


print("Aris's spherical-shell table at Lambda = 1")
print(f"{'p':>6}{'printed':>10}{'formula':>10}{'pymrm':>10}")
SHELL_PYMRM = {}
for p, pr in zip(p_print, eta_print):
    SHELL_PYMRM[p] = shell_pymrm(p)
    print(f"{p:6.2f}{pr:10.3f}{eta_shell(1.0, p):10.5f}{SHELL_PYMRM[p]:10.5f}")

claims.value("shell_formula_vs_printed", shell_formula_vs_printed())
claims.value("shell_formula_p075", eta_shell(1.0, 0.75), rtol=1e-4)
claims.value("shell_pymrm_vs_formula",
             max(abs(v - eta_shell(1.0, p)) for p, v in SHELL_PYMRM.items()))
SHELL_INNER = [p for p in p_print if 0 < p < 1]
SHELL_OPEN = {p: shell_pymrm(p, bc_low=EXPOSED) for p in SHELL_INNER}   # inner face exposed
claims.value("shell_exposed_inner_face_dev",
             max(abs(v - eta_shell(1.0, p)) for p, v in SHELL_OPEN.items()))

# the limits: how the expression approaches eqs. 17 and 15
EPS = 10.0 ** -np.arange(2, 6)


def limit_orders(**kw):
    """Orders of approach in p (to the sphere) and in 1 - p (to the plate)."""
    d_s = np.array([abs(eta_shell(1.0, e, **kw) - eta_sphere(1.0)) for e in EPS[:2]])
    d_p = np.array([abs(eta_shell(1.0, 1 - e, **kw) - eta_plate(1.0)) for e in EPS])
    return float(np.log10(d_s[0] / d_s[1])), float(np.log10(d_p[-2] / d_p[-1])), d_s, d_p


ord_s, ord_p, d_s, d_p = limit_orders()
claims.value("shell_sphere_limit_order", ord_s)
claims.value("shell_plate_limit_order", ord_p)
claims.value("shell_plate_gap_1mp_1e-5", d_p[-1])
print(f"\nformula vs printed : {claims['shell_formula_vs_printed']:.5f}")
print(f"pymrm vs formula   : {claims['shell_pymrm_vs_formula']:.1e}  (n = 400)")
print(f"the same with the inner face exposed by mistake: {claims['shell_exposed_inner_face_dev']:.2f}")
print("p -> 0, |shell - eq. 17| at p = " + ", ".join(f"{e:.0e}: {d:.2e}" for e, d in zip(EPS, d_s))
      + f"   order {ord_s:.3f} in p")
print("p -> 1, |shell - eq. 15| at 1-p = " + ", ".join(f"{e:.0e}: {d:.2e}" for e, d in zip(EPS, d_p))
      + f"   order {ord_p:.4f} in 1-p")

# %% [markdown]
# Six distinct terms, transcribed from a mangled scan, and they reproduce the five
# printed values *and* both stated limits: the expression approaches eq. 17 as
# $p^3$ and eq. 15 linearly in $1-p$, smoothly, not by a coincidence at one point.
#
# It also settles a disagreement between the two printings: at $p = 0.75$ the 1957
# OCR gives `728` and the 1995 reprint 0.723. The expression evaluates to 0.7227,
# so the reprint is right and the 1957 reading is corrupt.

# %%
claims.says("0.7227", "shell_formula_p075")

# %% [markdown]
# ### The finite cylinder, and which reading of *p* is right
#
# Case (iv) is the one shape in the paper with no evaluated answer, only a double
# series, and its series and its running text disagree about the aspect ratio by
# a factor of two. A two-dimensional solve decides it: build the cylinder from its
# geometry, compute $\eta$, and see which reading of the series matches.

# %%
FC_P, FC_L, FC_NR = (0.25, 0.5, 1.0, 2.0), (0.5, 1.0, 2.0, 5.0), 160
READINGS = {"nomenclature": lambda p: 1.0 / (2.0 * p), "running text": lambda p: 1.0 / p}


def fincyl_gaps(nu_r=1, n_r=FC_NR, record=False):
    """Worst |pymrm / series - 1| per reading of p, over FC_P x FC_L."""
    worst, table = {}, []
    for label, tau_of in READINGS.items():
        gaps = []
        for p in FC_P:
            tau = tau_of(p)
            pel, va = finite_cylinder(tau, n_r=n_r)
            if nu_r != 1:                                 # break table: wrong geometry factor
                pel = Pellet(pel.shape, (1.0, tau), (nu_r, 0))
            for L in FC_L:
                v, _ = eta_both(pel, L / va, key=f"fincyl {label} p={p} L={L}" if record else None)
                s = eta_fincyl_series(L, p, tau)
                gaps.append(abs(rel(s, v)))         # the series' error, relative to the solve
                if L == 1.0:
                    table.append((p, label, 2 * tau, pel.shape, v, s, gaps[-1]))
        worst[label] = max(gaps)
    return worst, table


FC_WORST, _tab = fincyl_gaps(record=True)
print("finite cylinder: 2-D pymrm against the series exactly as printed (Lambda = 1 shown)")
print(f"{'p':>6}{'reading':>14}{'length/radius':>15}{'grid':>12}{'pymrm':>10}{'series':>10}{'rel':>9}")
for p, label, l_a, shp, v, s, g in _tab:
    print(f"{p:6.2f}{label:>14}{l_a:15.2f}{str(shp):>12}{v:10.5f}{s:10.5f}{g:9.1e}")
claims.value("fincyl_series_p_as_nomenclature", FC_WORST["nomenclature"])
claims.value("fincyl_series_p_as_running_text", FC_WORST["running text"])
print(f"\nworst over p in {FC_P} and Lambda in {FC_L}:")
print(f"  p as radius:length (the nomenclature, length = a/p): {FC_WORST['nomenclature']:.1e}")
print(f"  p as in the running text (length = 2a/p)            : {FC_WORST['running text']:.1%}")

# %% [markdown]
# The series is right as printed, provided $p$ is the **radius:length** ratio its
# own nomenclature defines. Taking the running text's "length $= 2a/p$" at face
# value puts it out by tens of percent. The 2-D solve settles this without having
# to guess which the author meant; an independent derivation of the series would
# settle it too, but the paper never evaluates its series, so it cannot settle
# it itself. The remaining gap under the nomenclature
# reading is the grid error of the 2-D solve: the validation section refines the
# worst case on a ladder and its extrapolated limit lands on the series.

# %%
display(Markdown(
    f"Under the running-text reading the worst gap is "
    f"{100 * claims['fincyl_series_p_as_running_text']:.0f} %; under the nomenclature "
    f"reading {claims['fincyl_series_p_as_nomenclature']:.1e}, on {FC_NR} radial cells."))

# %% [markdown]
# ### Shapes with corners
#
# The square rod has the same $v_p/s_x = a/2$ as the long circular cylinder, and
# the cube has the same $v_p/s_x = a/3$ as the sphere. At equal $\Lambda$ they are
# therefore a direct test of Aris's proposition with everything except *shape*
# held fixed. Both are checked against an independent eigenfunction expansion
# before being used.

# %%
ROD_L = (0.2, 1.0, 5.0)
rod, rod_va = square_rod(n=200)
ROD = {L: eta_both(rod, L / rod_va, key=f"rod L={L}")[0] for L in ROD_L}
claims.value("square_rod_vs_series",
             max(abs(rel(ROD[L], eta_box_series(L, (1, 1), rod_va, 800))) for L in ROD_L))

CUBE_NS = (8, 16, 32)
CUBE_VA = 1.0 / 3.0
CUBE_SOL = {}


def cube_eta(n, lam=1.0 / CUBE_VA):
    pel, _ = cube(n)
    v, chi = eta_both(pel, lam, key=f"cube n={n}")
    CUBE_SOL[n] = (pel, chi)
    return v


cube_lad = observed_orders(cube_eta, CUBE_NS)
cube_ref = eta_box_series(1.0, (1, 1, 1), CUBE_VA, 200)
eta_cube32 = cube_lad["values"][-1]
claims.value("cube_vs_series_n32", abs(rel(eta_cube32, cube_ref)))
claims.bound("cube_order_8_16_32", cube_lad["orders"][-1], min=2.0, max=3.0, basis="chosen",
             note="second-order scheme; 8/16/32 is pre-asymptotic (1-D ladders: 2.0 at Lambda <= 1, "
                  "up to 3.0 where the layer is thin)")
claims.bound("cube_extrapolated_vs_series", abs(rel(cube_lad["extrapolated"], cube_ref)),
             max=2e-5, basis="chosen", note="Richardson limit of 8/16/32 against the series, "
             "below the n = 32 error (7.5e-5)")
claims.bound("cube_ladder_finest_diff", abs(rel(cube_lad["values"][-2], cube_lad["values"][-1])),
             min=1e-7, basis="chosen", note="ten times the flux-identity bound (1e-8)")

lam20 = 20.0 / CUBE_VA
cube24, _ = cube(24)
eta20 = {24: cube24.eta(lam20), 32: CUBE_SOL[32][0].eta(lam20)}
claims.value("cube_n24_n32_gap_Lambda20", abs(rel(eta20[24], eta20[32])))

print(f"square rod, 2-D pymrm (n = 200) vs its series : {claims['square_rod_vs_series']:.1e}")
print("cube, 3-D pymrm at Lambda = 1 against its series:")
for n, v in zip(CUBE_NS, cube_lad["values"]):
    print(f"   n = {n:2d}^3  eta {v:.6f}  error {rel(v, cube_ref):+.2e}")
print(f"   order from successive differences {cube_lad['orders'][-1]:.2f}; "
      f"Richardson limit off the series by {rel(cube_lad['extrapolated'], cube_ref):+.1e}")
print(f"cube at Lambda = 20: n = 24 and n = 32 differ by "
      f"{claims['cube_n24_n32_gap_Lambda20']:.0%}: the O(1/lambda) reacting layer "
      "is not resolved on these uniform 3-D grids")

eq_cyl, eq_va = finite_cylinder(1.0, n_r=200)     # length 2a, the 'equilateral' cylinder
eta_eq_cyl = eta_both(eq_cyl, 1.0 / eq_va, key="cylinder l=2a L=1")[0]


def order_label(values):
    """Shapes at v_p/s_x = a/3 and Lambda = 1, highest eta first, and the smallest gap."""
    s = sorted(values.items(), key=lambda kv: -kv[1])
    return " > ".join(k for k, _ in s), min(a[1] - b[1] for a, b in zip(s, s[1:]))


EQ_A3 = {"sphere": eta_sphere(1.0), "cylinder l=2a": eta_eq_cyl, "cube": eta_cube32}
_lab, _gap = order_label(EQ_A3)
claims.exact("order_at_va_a3_Lambda1", _lab, margin=_gap, min_margin=1e-3)
print("\nEqual v_p/s_x, different shape, at Lambda = 1:")
print(f"  v_p/s_x = a/2   long circular cylinder {eta_cylinder(1.0):.4f}"
      f"   long square rod {ROD[1.0]:.4f}")
print("  v_p/s_x = a/3   " + "   ".join(f"{k} {v:.4f}" for k, v in EQ_A3.items()))
print(f"  ordering: {_lab}  (smallest gap {_gap:.4f})")

# %% [markdown]
# **Where the cube's curve comes from.** The reacting layer thins as
# $1/\lambda$, and the $24^3$ and $32^3$ uniform grids disagree badly at
# $\Lambda = 20$, as printed above. So the sweep below takes the cube from its
# eigenfunction series at every $\Lambda$; the 3-D pymrm ladder checks that series
# at $\Lambda = 1$ only. The flat plate, long cylinder, sphere and the four shells come
# from their closed forms; the square rod and the five finite cylinders are pymrm
# solves on grids graded toward the exposed faces, each compared with its series
# at every swept $\Lambda$.

# %% [markdown]
# ### How wide is the band, really?
#
# Fourteen shapes: Aris's three, four spherical shells, five finite cylinders from
# a long rod to a flat disc, the square rod and the cube.

# %%
Lams = np.geomspace(0.05, 20.0, 28)
SWEEP_Q = 2.0
TAUS = (2.0, 1.0, 0.5, 0.25, 0.125)


def cyl_name(tau):
    return f"cylinder l/a={2 * tau:g}"


curves = {"flat plate": eta_plate(Lams), "long cylinder": eta_cylinder(Lams),
          "sphere": eta_sphere(Lams)}
for p in (0.25, 0.5, 0.75, 0.9):
    curves[f"shell p={p}"] = np.array([eta_shell(L, p) for L in Lams])
refs = dict(curves)
rod_q, _ = square_rod(n=200, q=SWEEP_Q)
curves["square rod"] = np.array([rod_q.eta(L / rod_va) for L in Lams])
refs["square rod"] = np.array([eta_box_series(L, (1, 1), rod_va, 800) for L in Lams])
curves["cube"] = refs["cube"] = np.array([eta_box_series(L, (1, 1, 1), CUBE_VA, 200) for L in Lams])
for tau in TAUS:
    pel, va = finite_cylinder(tau, q=SWEEP_Q)
    curves[cyl_name(tau)] = np.array([pel.eta(L / va) for L in Lams])
    refs[cyl_name(tau)] = np.array([eta_fincyl_series(L, 1 / (2 * tau), tau) for L in Lams])

SOLVED = ["square rod"] + [cyl_name(t) for t in TAUS]
sweep_err = {k: np.abs(curves[k] / refs[k] - 1) for k in SOLVED}
claims.value("sweep_pymrm_vs_series_worst", max(e.max() for e in sweep_err.values()))
print("pymrm solves in the sweep against their series, worst over all 28 Lambda:")
for k, e in sweep_err.items():
    print(f"   {k:22s} {e.max():.1e} at Lambda = {Lams[e.argmax()]:.3g}")


def disc_uniform_worst(n_r=120):
    """The flattest disc of the sweep on a uniform grid: worst error over the sweep."""
    pel, va = finite_cylinder(TAUS[-1], n_r=n_r, q=1.0)
    e = np.abs(np.array([pel.eta(L / va) for L in Lams]) / refs[cyl_name(TAUS[-1])] - 1)
    return float(e.max()), float(Lams[e.argmax()])


_du = disc_uniform_worst()
claims.value("sweep_disc_uniform_grid_worst_pct", 100 * _du[0])
print(f"   the same {cyl_name(TAUS[-1])} on a uniform grid: {_du[0]:.1e} at Lambda = {_du[1]:.3g}")

names = list(curves)


def cube_lowest_margins(R=refs):
    """cube - lowest other reference, per Lambda (negative: the cube is lowest)."""
    others = np.array([R[k] for k in R if k != "cube"])
    return R["cube"] - others.min(axis=0)


def others_minus_plate(R=refs):
    return float((np.array([R[k] for k in R if k != "flat plate"]).max(axis=0)
                  - R["flat plate"]).max())


N_CUBE_LOW = claims.count("sweep_points_cube_lowest", cube_lowest_margins(), 0.0,
                          above=False, min_margin=1e-6)
def cyl2a_between_margins(R=refs):
    """min(sphere - cylinder l=2a, cylinder l=2a - cube), per Lambda (positive: between)."""
    c = R[cyl_name(1.0)]
    return np.minimum(R["sphere"] - c, c - R["cube"])


N_CYL_BETWEEN = claims.count("sweep_points_cyl2a_between_sphere_cube", cyl2a_between_margins(),
                             0.0, min_margin=1e-5)
claims.bound("sweep_others_minus_plate_max", others_minus_plate(), max=0.0, basis="chosen",
             note="the plate is the upper edge of the band at every swept Lambda (references)")
print(f"\non the references: the cube is the lowest curve at {N_CUBE_LOW} of {len(Lams)} Lambda, "
      f"the cylinder of length 2a between sphere and cube at {N_CYL_BETWEEN}, "
      f"the plate the highest at all (largest excess of another shape over it "
      f"{claims['sweep_others_minus_plate_max']:.1e})")


# The band's edges are the plate (eq. 15) and the cube (series), so its widest
# point is found on those two curves, not read off the 28 sweep points.
def widest(upper, lower, lo=0.3, hi=3.0):
    r = minimize_scalar(lambda L: -(upper(L) - lower(L)), bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-6})
    return -r.fun, r.x


def eta_cube_series(L, va=CUBE_VA):
    return eta_box_series(L, (1, 1, 1), va, 200)


band_all, L_band = widest(eta_plate, eta_cube_series)
band_aris, L_aris = widest(eta_plate, eta_sphere)
below_abs, L_below = widest(eta_sphere, eta_cube_series)
r_rel = minimize_scalar(lambda L: eta_cube_series(L) / eta_sphere(L), bounds=(0.3, 20.0),
                        method="bounded", options={"xatol": 1e-6})
below_rel, L_below_rel = 1 - r_rel.fun, r_rel.x
claims.value("band_all_shapes", band_all)
claims.value("band_all_Lambda", L_band)
claims.value("band_aris_plate_sphere_max", band_aris)
claims.value("band_aris_Lambda", L_aris)
claims.value("band_ratio_all_over_aris", band_all / band_aris)
claims.value("cube_below_sphere", below_abs)
claims.value("cube_below_sphere_pct_there", 100 * below_abs / eta_sphere(L_below))
claims.value("cube_below_sphere_max_pct", 100 * below_rel)
claims.value("cube_below_sphere_max_pct_Lambda", L_below_rel)
print(f"widest band over all 14 shapes : {band_all:.4f} at Lambda = {L_band:.3f} (plate to cube)")
print(f"Aris's plate-sphere band       : {band_aris:.4f} at Lambda = {L_aris:.3f}"
      f"   -> the full band is {band_all / band_aris:.3f} x wider")
print(f"deepest excursion below sphere : {below_abs:.4f} at Lambda = {L_below:.3f} "
      f"({claims['cube_below_sphere_pct_there']:.2f} % of the sphere's eta there); "
      f"relative, {below_rel:.2%} at Lambda = {L_below_rel:.2f}")

# %%
lo, hi = np.array([curves[k] for k in names]).min(0), np.array([curves[k] for k in names]).max(0)
fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
ax = axes[0]
ax.fill_between(Lams, lo, hi, color="tab:orange", alpha=0.30, lw=0, label="all 14 shapes")
ax.fill_between(Lams, eta_sphere(Lams), eta_plate(Lams), color="tab:blue",
                alpha=0.30, lw=0, label="Aris's plate-sphere band")
for nm, col, ls in (("flat plate", "k", "-"), ("sphere", "tab:blue", "-"),
                    ("cylinder l/a=2", "tab:green", "--"), ("cube", "tab:red", "-")):
    ax.semilogx(Lams, curves[nm], color=col, ls=ls, lw=1.8, label=nm)
for shp, mk in (("plate", "s"), ("cylinder", "^"), ("sphere", "o")):
    sub = tab1[tab1["shape"] == shp]
    ax.semilogx(sub["Lambda"], sub["eta_printed"], mk, ms=6, mfc="white", mec="0.15",
                mew=1.3, ls="none", label="Aris Table 1" if shp == "plate" else "_nolegend_")
ax.set_xlabel(r"generalised modulus $\Lambda = \lambda\,v_p/s_x$ [-]")
ax.set_ylabel(r"effectiveness factor $\eta$ [-]")
ax.set_ylim(0.0, 1.06)
ax.set_title("The collapse, and what it hides")
ax.legend(frameon=False, fontsize=8, loc="lower left")
ax = axes[1]
ax.semilogx(Lams, hi - lo, "-", color="tab:orange", lw=2.2, label="all 14 shapes")
ax.semilogx(Lams, eta_plate(Lams) - eta_sphere(Lams), "-", color="tab:blue", lw=2.2,
            label="plate minus sphere (Aris)")
ax.plot([L_band], [band_all], "o", color="tab:orange", ms=7)
ax.plot([L_aris], [band_aris], "o", color="tab:blue", ms=7)
ax.set_xlabel(r"$\Lambda$ [-]")
ax.set_ylabel(r"spread in $\eta$ at fixed $\Lambda$ [-]")
ax.set_title("Width of the band")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
plt.show()

# %%
display(Markdown(
    f"The collapse is real: over nearly three decades of $\\Lambda$ every shape sits "
    f"inside a band that never exceeds {band_all:.3f} in $\\eta$ and closes to nothing "
    f"at both ends. Its upper edge is the flat plate at every swept modulus and its lower "
    f"edge the cube at {N_CUBE_LOW} of the {len(Lams)}. But two things about it are not what the paper's three "
    f"shapes suggest.\n\n"
    f"**The band is about a quarter wider than plate-to-sphere.** Aris's plate-sphere "
    f"spread peaks at {band_aris:.4f} ($\\Lambda$ = {L_aris:.2f}). Adding corners takes the "
    f"peak to {band_all:.4f} ($\\Lambda$ = {L_band:.2f}), {band_all / band_aris - 1:.0%} wider.\n\n"
    f"**The sphere is not the floor.** Table 1 makes it look like the extreme case, "
    f"lowest $\\eta$ at every $\\Lambda$, and the paper's Leva argument reinforces the "
    f"impression. But at fixed $v_p/s_x$ the cube sits up to {below_abs:.4f} below the "
    f"sphere (at $\\Lambda$ = {L_below:.2f}, {claims['cube_below_sphere_pct_there']:.1f} % of the "
    f"sphere's $\\eta$ there; the largest relative gap is {below_rel:.1%}, at "
    f"$\\Lambda$ = {L_below_rel:.1f}), and the finite cylinder of length $2a$ sits between "
    f"them at {N_CYL_BETWEEN} of the {len(Lams)} swept moduli. All three have exactly the same $v_p/s_x = a/3$."))

# %% [markdown]
# ### Two different comparisons, two different answers
#
# Aris ends section 3 with: "of all catalyst particles of given volume, diffusion
# limitation is most severe in a spherical one." That is a statement about the
# Leva shape factor $\phi = (36\pi v_p^2)^{1/3}/s_x \le 1$: at **given volume**
# the sphere has the largest $v_p/s_x$, hence the largest $\Lambda$, hence the
# lowest $\eta$. It is correct, and it is a different comparison from the one
# Table 1 invites.

# %%
SHAPES = {                                    # (v_p, s_x) in units of a
    "sphere":              (4 * np.pi / 3, 4 * np.pi),
    "cube":                (8.0, 24.0),
    "cylinder, length 2a": (2 * np.pi, 6 * np.pi),
    "cylinder, length 4a": (4 * np.pi, 10 * np.pi),
    "cylinder, length a":  (np.pi, 4 * np.pi),
}
LAM_STAR = 3.0                                # modulus of the equal-volume sphere, illustrative

print(f"{'shape':22s}{'v_p/s_x':>10}{'Leva phi':>10}{'Lambda at equal volume':>25}")
for nm, (v, s) in SHAPES.items():
    phi = (36 * np.pi * v**2) ** (1 / 3) / s
    print(f"{nm:22s}{v / s:10.4f}{phi:10.4f}{phi * LAM_STAR:25.3f}")
print("\nAt equal PARTICLE VOLUME every non-spherical particle sits at a smaller")
print("Lambda than the equal-volume sphere, hence a HIGHER eta: Aris's claim.")
print(f"\nAt equal v_p/s_x = a/3 and Lambda = 1: {claims['order_at_va_a3_Lambda1']}")

# %% [markdown]
# Both statements hold and they are not in conflict:
#
# - **At equal volume** the sphere is the worst case, because it has the least
#   surface per unit volume and therefore the largest $\Lambda$. This is Aris's
#   claim, and the Leva shape factor proves it.
# - **At equal $\Lambda$**, that is after the collapse has already been applied,
#   the sphere is the *best* of the three shapes that share $v_p/s_x = a/3$. The
#   cube has corners and edges, the finite cylinder edges only, the sphere neither,
#   and $\eta$ falls in that order: near an edge or a corner reactant arrives
#   from two or three faces at once.
#
# The second is the one that matters when you use the shortcut, and it is the one
# Table 1 cannot show, because all three of its shapes are smooth.

# %% [markdown]
# ### Eq. 18: the asymptote that makes the collapse work
#
# $\eta\Lambda \to 1$ as $\Lambda\to\infty$, for *every* shape. This is the
# paper's one rigorous general result, and unlike Table 1 it can be tested on
# shapes with edges. $\eta\Lambda$ is extrapolated to $1/\Lambda = 0$ by a straight
# line through $\Lambda$ = 10, 20, 50, for the exact references and for pymrm
# solves on graded grids, so that the error of the fit and the error of the grid
# are seen separately.

# %%
BIG = np.array([5.0, 10.0, 20.0, 50.0, 100.0])
FIT = slice(1, 4)
ASYM_SHAPES = {   # name: (pymrm pellet maker, its reference)
    "square rod": (lambda q: square_rod(n=200, q=q),
                   lambda L: eta_box_series(L, (1, 1), 0.5, 1600)),
    "cylinder l=2a": (lambda q: finite_cylinder(1.0, n_r=200, q=q),
                      lambda L: eta_fincyl_axial(L, 0.5, 1.0)),
    "shell p=0.5": (lambda q: spherical_shell(0.5, n=400, q=q),
                    lambda L: eta_shell(L, 0.5)),
}


def extrapolated(etaL, fit=FIT):
    return float(np.polyfit(1.0 / BIG[fit], etaL[fit], 1)[1])


def asym_pymrm(q=2.0):
    out = {}
    for nm, (make, _) in ASYM_SHAPES.items():
        pel, va = make(q)
        out[nm] = np.array([pel.eta(L / va) for L in BIG]) * BIG
    return out


asym_ref = {"sphere (eq. 17)": eta_sphere(BIG) * BIG, "flat plate (eq. 15)": eta_plate(BIG) * BIG}
for nm, (_, ref) in ASYM_SHAPES.items():
    asym_ref[nm] = np.array([ref(L) for L in BIG]) * BIG
asym = asym_pymrm()


def asym_metrics(A, R=asym_ref, fit=FIT):
    ext_dev = max(abs(extrapolated(v, fit) - 1) for v in A.values())
    grid50 = max(abs(rel(A[k][3], R[k][3])) for k in A)
    return ext_dev, grid50


claims.value("eq18_cyl_reference_truncation_Lambda50",
             abs(rel(eta_fincyl_axial(50.0, 0.5, 1.0), eta_fincyl_axial(50.0, 0.5, 1.0, 12800))))
claims.bound("fincyl_axial_vs_double_series", abs(rel(eta_fincyl_axial(1.0, 0.5, 1.0),
             eta_fincyl_series(1.0, 0.5, 1.0))), max=1e-6, basis="chosen",
             note="the two forms of case (iv) at Lambda = 1, l = 2a; both truncated")
print(f"case (iv), closed-form axial sum against the double series at Lambda = 1: "
      f"{claims['fincyl_axial_vs_double_series']:.0e}; its truncation at Lambda = 50 "
      f"(6400 against 12800 radial terms): {claims['eq18_cyl_reference_truncation_Lambda50']:.0e}")
claims.value("eq18_reference_extrapolated_worst_dev",
             max(abs(extrapolated(v) - 1) for v in asym_ref.values()))
_ext, _g50 = asym_metrics(asym)
claims.value("eq18_extrapolated_worst_dev", _ext)
claims.value("eq18_pymrm_worst_rel_err_Lambda50", _g50)

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for nm, v in asym_ref.items():
    ax.plot(1.0 / BIG, v, "-", lw=1.2, alpha=0.6)
for nm, v in asym.items():
    ax.plot(1.0 / BIG, v, "o", ms=5, label=f"{nm} (pymrm)")
ax.axhline(1.0, color="k", lw=1.0, ls="--")
ax.set_xlabel(r"$1/\Lambda$ [-]")
ax.set_ylabel(r"$\eta\,\Lambda$ [-]")
ax.set_xlim(-0.008, 0.21)
ax.set_title(r"Aris eq. 18: $\eta\Lambda \to 1$ for every shape (lines: references)")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
plt.show()

print(f"{'eta*Lambda':30s}" + "".join(f"{L:>9.0f}" for L in BIG) + "   extrapolated")
for nm, v in asym_ref.items():
    print(f"{nm + ' (reference)':30s}" + "".join(f"{a:9.5f}" for a in v) + f"   {extrapolated(v):.5f}")
for nm, v in asym.items():
    print(f"{nm + ' (pymrm)':30s}" + "".join(f"{a:9.5f}" for a in v) + f"   {extrapolated(v):.5f}")
print(f"\nworst deviation of the extrapolated value from 1: references "
      f"{claims['eq18_reference_extrapolated_worst_dev']:.1e} (the straight-line fit), pymrm "
      f"{claims['eq18_extrapolated_worst_dev']:.1e}; pymrm grid error at Lambda = 50 up to "
      f"{claims['eq18_pymrm_worst_rel_err_Lambda50']:.1e}")

# %% [markdown]
# Every shape heads to $\eta\Lambda = 1$, the corner shapes included. On the
# references the straight-line extrapolation lands on 1 to within its own fit
# error (the curvature in $1/\Lambda$; the cylinder reference is the closed-form
# axial sum, truncated far below that). The pymrm values extrapolate a little
# further from 1 because they carry a grid error at $\Lambda = 50$, where
# $\lambda h$ near the exposed faces is no longer small even on graded grids; that
# grid error, not the asymptote, is most of their deviation.

# %% [markdown]
# ### The mixture rule
#
# Aris's section 4: for a bed of mixed sizes the effectiveness must be averaged,
# not the size. He gives one worked example to show how badly the shortcut
# fails.

# %%
f_i = mix.loc[mix["quantity"] == "fraction", "value"].to_numpy()
Lam_i = mix.loc[mix["quantity"] == "Lambda", "value"].to_numpy()
PRINTED = {q: mix.loc[mix["quantity"] == q, ["value", "printing"]].to_numpy()
           for q in ("mean_modulus", "eta_at_mean_modulus", "eta_mixture")}


def mixture(mean=lambda f, L: 1.0 / np.sum(f / L), curve=eta_plate):
    Lam_bar = mean(f_i, Lam_i)                             # his eq. 20
    return Lam_bar, curve(Lam_bar), float(np.sum(f_i * curve(Lam_i)))   # and eq. 19


Lam_bar, eta_from_mean, eta_correct = mixture()
claims.value("mixture_mean_modulus_dev", abs(Lam_bar - PRINTED["mean_modulus"][0, 0]))
claims.value("mixture_eta_at_mean_dev", abs(eta_from_mean - PRINTED["eta_at_mean_modulus"][0, 0]))
claims.value("mixture_eta_eq19", eta_correct)
print(f"eq. 20, mean modulus  : {Lam_bar:.4f}   (printed {PRINTED['mean_modulus'][0, 0]})")
print(f"eta at that mean      : {eta_from_mean:.4f}   (printed {PRINTED['eta_at_mean_modulus'][0, 0]})")
print(f"eq. 19, averaging eta : {eta_correct:.4f}   (printed "
      + ", ".join(f"{v} in the {pr} printing" for v, pr in PRINTED["eta_mixture"]) + ")")

# %% [markdown]
# The first two reproduce the printed values to their last digit, which pins down
# both the formula and the fact that the example uses the flat-plate curve. The
# third does not, and the two printings disagree with each other as well as with
# eq. 19, so no reading of it can be established; it is left as printed. The point
# of the example, that taking the mean size first overstates $\eta$, holds against
# eq. 19 whatever the third printed number is.

# %% [markdown]
# ## Validation
#
# None of these checks compares the solver with itself.
#
# | # | check |
# |---|---|
# | 1 | eqs. 15-17 against Aris's Table 1 (above) |
# | 2 | the shell expression against his five printed values and its two stated limits (above) |
# | 3 | pymrm 1-D against eqs. 15-17 and the shell expression, with a grid ladder |
# | 4 | pymrm 2-D finite cylinder against the printed double series, both readings of $p$, with a grid ladder on the worst case |
# | 5 | pymrm 2-D rod and 3-D cube against eigenfunction expansions, with a 3-D grid ladder (above) |
# | 6 | the truncation of every series, by doubling its terms |
# | 7 | eq. 18, $\eta\Lambda \to 1$, on five shapes including two with edges (above) |
# | 8 | eqs. 19-20 against his worked mixture example (above) |
#
# **Check 2 is the one the paper pays for.** The shell expression has six distinct
# terms and had to be transcribed from a mangled scan. It must reduce to eq. 17 at
# $p = 0$ and to eq. 15 at $p = 1$, and it must reproduce five printed values in
# between. A mis-read coefficient would show in the limits or in the values.
#
# **Check 4 is the one that could not be settled from the page at all.** The
# printed series and the printed geometry disagree by a factor of two in the aspect
# ratio, and since the paper never evaluates the series there is nothing on the
# page to break the tie. An independent solution of the actual shape breaks it.

# %%
PYMRM_1D_L = (0.1, 1.0, 5.0, 20.0)


def pymrm_1d(nu_sphere=2):
    out = {}
    for nm, (pel, va), fn, Ls in (
            ("slab", slab(), eta_plate, PYMRM_1D_L),
            ("long cylinder", long_cylinder(), eta_cylinder, PYMRM_1D_L),
            ("sphere", (Pellet((400,), (1.0,), (nu_sphere,)), 1.0 / 3.0), eta_sphere, PYMRM_1D_L),
            ("shell p=0.5", spherical_shell(0.5), lambda L: eta_shell(L, 0.5), PYMRM_1D_L[:3])):
        out[nm] = max(abs(rel(eta_both(pel, L / va, key=f"1-D {nm} L={L}" if nu_sphere == 2 else None)[0],
                              fn(L))) for L in Ls)
    return out


P1D = pymrm_1d()
claims.value("pymrm_1d_worst_rel", max(P1D.values()))
for nm, v in P1D.items():
    print(f"  pymrm 1-D {nm:14s} worst relative error {v:.1e}  (n = 400)")

# a grid ladder on the sphere, both eta routes from the same solves
SPH_NS, SPH_L = (100, 200, 400, 800, 1600), 1.0
_sph = {}


def sphere_eta(n, read=None):
    if n not in _sph:
        pel, va = sphere(n)
        _sph[n] = (pel, pel.solve(SPH_L / va), va)
    pel, chi, va = _sph[n]
    FLUX_GAPS[f"sphere ladder n={n}"] = abs(pel.eta_flux(SPH_L / va, chi) - pel.eta(SPH_L / va, chi))
    return pel.eta(SPH_L / va, chi) if read is None else pel.eta_flux(SPH_L / va, chi, read=read)


sph_lad = observed_orders(sphere_eta, SPH_NS)
claims.bound("sphere_order_Lambda1", sph_lad["orders"][-1], min=1.9, max=2.1, basis="theory",
             note="second-order finite volumes, layer resolved at Lambda = 1")
claims.bound("sphere_ladder_finest_diff", abs(rel(sph_lad["values"][-2], sph_lad["values"][-1])),
             min=1e-7, basis="chosen", note="ten times the flux-identity bound (1e-8), "
             "which bounds what the solve residual does to eta")
claims.bound("sphere_extrapolated_vs_eq17", abs(rel(sph_lad["extrapolated"], eta_sphere(SPH_L))),
             max=1e-8, basis="chosen", note="Richardson limit of the ladder against eq. 17")
print(f"\nsphere ladder at Lambda = {SPH_L:g}, n = {SPH_NS}: orders "
      + ", ".join(f"{o:.3f}" for o in sph_lad["orders"])
      + f"; Richardson limit against eq. 17 {rel(sph_lad['extrapolated'], eta_sphere(SPH_L)):+.1e}")

# %%
# The worst case of check 4 (p = 2 under the nomenclature, tau = 1/4) on a ladder
# that refines both axes together, ending one level above the production grid.
FC_LAD_P, FC_LAD_L, FC_LAD_K = 2.0, 5.0, (1, 2, 4, 8)
FC_LAD_TAU = READINGS["nomenclature"](FC_LAD_P)
_fc = {}


def fincyl_eta(k, read=None):
    if k not in _fc:
        pel, va = finite_cylinder(FC_LAD_TAU, n_r=40 * k, n_z=10 * k)
        _fc[k] = (pel, pel.solve(FC_LAD_L / va), va)
    pel, chi, va = _fc[k]
    if read is None:
        FLUX_GAPS[f"fincyl ladder k={k}"] = abs(pel.eta_flux(FC_LAD_L / va, chi) - pel.eta(FC_LAD_L / va, chi))
        return pel.eta(FC_LAD_L / va, chi)
    return pel.eta_flux(FC_LAD_L / va, chi, read=read)


fc_lad = observed_orders(fincyl_eta, FC_LAD_K, ratio=2.0)
fc_ref = eta_fincyl_series(FC_LAD_L, FC_LAD_P, FC_LAD_TAU)
claims.bound("fincyl_order_finest", fc_lad["orders"][-1], min=2.0, max=3.0, basis="chosen",
             note="second-order scheme, pre-asymptotic where the layer is thin (Lambda = 5)")
claims.bound("fincyl_extrapolated_vs_series", abs(rel(fc_lad["extrapolated"], fc_ref)),
             max=1e-5, basis="chosen", note="Richardson limit against the printed series; "
             "the series truncation is measured separately")
claims.bound("fincyl_ladder_finest_diff", abs(rel(fc_lad["values"][-2], fc_lad["values"][-1])),
             min=1e-7, basis="chosen", note="ten times the flux-identity bound (1e-8)")


def fincyl_all_extrapolated():
    """Richardson limit on n_r = 80, 160, 320 (n_z in proportion) for all 16 cases of
    check 4 under the nomenclature reading: worst gap to the series, with the volume
    average and with the two-point cell read (break table)."""
    worst = {"face": 0.0, "cell": 0.0}
    for p in FC_P:
        tau = READINGS["nomenclature"](p)
        n_z = finite_cylinder(tau, n_r=FC_NR)[0].shape[1]
        va = tau / (2.0 * tau + 1.0)
        pels = {k: finite_cylinder(tau, n_r=FC_NR * k // 4, n_z=n_z * k // 4)[0] for k in (2, 4, 8)}
        for L in FC_L:
            chis = {k: pel.solve(L / va) for k, pel in pels.items()}
            s_ref = eta_fincyl_series(L, p, tau)
            for read in worst:
                f = (lambda k: pels[k].eta(L / va, chis[k])) if read == "face" else \
                    (lambda k: pels[k].eta_flux(L / va, chis[k], read="cell"))
                lad = observed_orders(f, (2, 4, 8), ratio=2.0)
                worst[read] = max(worst[read], abs(rel(lad["extrapolated"], s_ref)))
    return worst


FC_ALL = fincyl_all_extrapolated()
claims.bound("fincyl_extrapolated_vs_series_all16", FC_ALL["face"], max=1e-5, basis="chosen",
             note="Richardson limits on 80/160/320 radial cells, all 16 cases of check 4")
print(f"all 16 cases, Richardson limit on n_r = 80/{FC_NR}/320 against the series: worst "
      f"{FC_ALL['face']:.1e}")
print(f"finite cylinder p = {FC_LAD_P:g}, Lambda = {FC_LAD_L:g}, grids "
      + ", ".join(str(_fc[k][0].shape) for k in FC_LAD_K))
print("   error against the series: "
      + ", ".join(f"{rel(v, fc_ref):+.2e}" for v in fc_lad["values"]))
print("   orders " + ", ".join(f"{o:.2f}" for o in fc_lad["orders"])
      + f"; Richardson limit against the series {rel(fc_lad['extrapolated'], fc_ref):+.1e}")

# %%
# Series truncation: each series against itself with twice the terms, at the
# settings the page uses and over every case each comparison covers.
def trunc_fincyl(scale=1):
    M, N = 400 // scale, 3000 // scale
    return max(abs(rel(eta_fincyl_series(L, p, READINGS["nomenclature"](p), M, N),
                       eta_fincyl_series(L, p, READINGS["nomenclature"](p), 2 * M, 2 * N)))
               for p in FC_P for L in FC_L)


def trunc_sweep(scale=1):
    """Worst relative truncation of the sweep's series curves, at the largest Lambda."""
    L = Lams[-1]
    out = [abs(rel(eta_box_series(L, (1, 1, 1), CUBE_VA, 200 // scale),
                   eta_box_series(L, (1, 1, 1), CUBE_VA, 400 // scale))),
           abs(rel(eta_box_series(L, (1, 1), 0.5, 800 // scale),
                   eta_box_series(L, (1, 1), 0.5, 1600 // scale)))]
    out += [abs(rel(eta_fincyl_series(L, 1 / (2 * t), t, 400 // scale, 3000 // scale),
                    eta_fincyl_series(L, 1 / (2 * t), t, 800 // scale, 6000 // scale))) for t in TAUS]
    return max(out)


claims.value("fincyl_series_truncation_worst", trunc_fincyl())
claims.value("cube_series_truncation_Lambda1",
             abs(rel(cube_ref, eta_box_series(1.0, (1, 1, 1), CUBE_VA, 400))))
claims.value("rod_series_truncation_worst",
             max(abs(rel(eta_box_series(L, (1, 1), 0.5, 800), eta_box_series(L, (1, 1), 0.5, 1600)))
                 for L in ROD_L))
claims.value("sweep_series_truncation_Lambda20", trunc_sweep())
for k in ("fincyl_series_truncation_worst", "cube_series_truncation_Lambda1",
          "rod_series_truncation_worst", "sweep_series_truncation_Lambda20"):
    print(f"   {k:36s} {claims[k]:.1e}")

# %%
# The two eta routes: the same sum rearranged, so their gap is the solve residual.
claims.bound("eta_flux_vs_volume_worst", max(FLUX_GAPS.values()), max=1e-8, basis="identity",
             note="finite-volume conservation; the gap is the residual of the direct solve")
_kmax = max(FLUX_GAPS, key=FLUX_GAPS.get)
print(f"surface-flux route against the volume average, worst over {len(FLUX_GAPS)} solves: "
      f"below 1e-8 (here {claims['eta_flux_vs_volume_worst']:.0e}, {_kmax})")

# %%
display(Markdown(
    f"Checks 3 to 6 in numbers. In 1-D, pymrm reproduces eqs. 15-17 and the shell to "
    f"{claims['pymrm_1d_worst_rel']:.1e} at $n = 400$ (worst over $\\Lambda$ up to 20; "
    f"the error is largest where the layer is thinnest), and the sphere ladder at "
    f"$\\Lambda = 1$ is second order ({claims['sphere_order_Lambda1']:.3f}). The finite-"
    f"cylinder gap under the nomenclature reading is grid error: on the worst case the "
    f"ladder converges at order {claims['fincyl_order_finest']:.2f} (pre-asymptotic, the "
    f"layer is thin at $\\Lambda = 5$) and its Richardson limit lies "
    f"{claims['fincyl_extrapolated_vs_series']:.0e} from the series (over all 16 cases, "
    f"Richardson limits on 80/160/320 radial cells lie within "
    f"{claims['fincyl_extrapolated_vs_series_all16']:.0e}). The cube's 8/16/32 "
    f"ladder has order {claims['cube_order_8_16_32']:.2f} and a limit "
    f"{claims['cube_extrapolated_vs_series']:.0e} from its series. Every series is "
    f"truncated far below the gaps it is used to judge: at most "
    f"{max(claims[k] for k in ('fincyl_series_truncation_worst', 'cube_series_truncation_Lambda1', 'rod_series_truncation_worst')):.0e} "
    f"in the validation comparisons, and {claims['sweep_series_truncation_Lambda20']:.0e} "
    f"at $\\Lambda = 20$ in the sweep, against a pymrm-to-series gap there of up to "
    f"{claims['sweep_pymrm_vs_series_worst']:.1e}. The eq. 18 cylinder reference is "
    f"truncated by {claims['eq18_cyl_reference_truncation_Lambda50']:.0e} at "
    f"$\\Lambda = 50$, against a line-fit deviation of "
    f"{claims['eq18_reference_extrapolated_worst_dev']:.1e}."))

# %% [markdown]
# **What is *not* validated here.** Nothing on this page is compared with a
# measurement, and no version of it could be: Aris (1957) reports no experiments.
# Every number above is either the author's own arithmetic or a mathematical
# identity. That is provenance tier 6, and it is the honest ceiling for this
# source.
#
# The widened band is a statement about the *idealised* first-order isothermal
# problem, with a uniform surface concentration and no external film. A real
# pellet with a film resistance, a distribution of pore sizes or non-first-order
# kinetics will depart from the collapse for reasons that have nothing to do with
# shape: `B1.3` (Bischoff) is the kinetics half of that story and Aris's own
# section 5 is the film half. Neither is reproduced here.
#
# And the cube is the numerically weakest shape on the page: its curve in the
# sweep is the eigenfunction series at every $\Lambda$, not a pymrm solve. The
# series is checked against the 3-D solve at $\Lambda = 1$, one modulus, not
# across the sweep.

# %% [markdown]
# ### The break table
#
# Each row re-runs the affected computation with one thing deliberately wrong,
# injected through an argument, an input or a formula, and hands the result to
# `claims.breaks`, which records whether that claim's own check fails.

# %%
# ---- Table 1
_lab = "cylinder column evaluated with the modulus on the radius (eq. 16 at 2 Lambda)"
_c2 = {**CLOSED, "cylinder": lambda L: eta_cylinder(2 * L)}
_d2 = table1_diffs(_c2)
brk(_lab, "table1_worst_19_entries", np.abs(_d2[np.abs(_d2) <= FLAG]).max())
brk(_lab, "table1_entries_off_more_than_0p002", int(np.sum(np.abs(_d2) > FLAG)))
brk(_lab, "table1_entries_not_rounding_to_printed", int(np.sum(np.abs(_d2) > HALF_UNIT)))
for L in (5.0, 10.0):
    brk(_lab, f"cyl_printed_high_pct_Lambda{L:g}", 100 * rel(CYL_SLIP[L][0], eta_cylinder(2 * L)))
    brk(_lab, f"cyl_printed_needs_Lambda_factor_{L:g}",
        lambda_needed(CYL_SLIP[L][0], lambda x: eta_cylinder(2 * x))["roots"][0] / L)
_lab = "sphere column evaluated with the modulus on the radius (eq. 17 at Lambda / 3)"
_d3 = table1_diffs({**CLOSED, "sphere": lambda L: eta_sphere(L / 3)})
brk(_lab, "table1_plate_sphere_worst_Lambda_5_10", np.abs(_d3[_ps5.to_numpy()]).max())

# ---- the shell
_lab = _lab_shell = "shell Lambda' built on the outer radius instead of the thickness"
brk(_lab, "shell_formula_vs_printed", shell_formula_vs_printed(on_thickness=False))
brk(_lab, "shell_formula_p075", eta_shell(1.0, 0.75, on_thickness=False))
_o = limit_orders(on_thickness=False)
brk(_lab, "shell_sphere_limit_order", _o[0])
brk(_lab, "shell_plate_limit_order", _o[1])
brk(_lab, "shell_plate_gap_1mp_1e-5", _o[3][-1])
_ocr = shell_tab["eta_1957_ocr"].to_numpy()
brk("shell table with the 1957 OCR reading where it differs (p = 0.75)", "shell_formula_vs_printed",
    shell_formula_vs_printed(np.where(np.isnan(_ocr), eta_print, _ocr)))
brk("shell's inner face exposed by mistake (EXPOSED at r = pa)", "shell_pymrm_vs_formula",
    claims["shell_exposed_inner_face_dev"])
brk(_lab_shell, "shell_exposed_inner_face_dev",
    max(abs(v - eta_shell(1.0, p, on_thickness=False)) for p, v in SHELL_OPEN.items()))

# ---- 1-D and the ladders
brk("sphere assembled with nu = 1 in construct_div", "pymrm_1d_worst_rel",
    max(pymrm_1d(nu_sphere=1).values()))
_cell = observed_orders(lambda n: sphere_eta(n, read="cell"), SPH_NS)
_lab = "eta read from a two-point gradient at the last cell centre"
brk(_lab, "sphere_order_Lambda1", _cell["orders"][-1])
brk(_lab, "sphere_extrapolated_vs_eq17", abs(rel(_cell["extrapolated"], eta_sphere(SPH_L))))
_fcc = observed_orders(lambda k: fincyl_eta(k, read="cell"), FC_LAD_K, ratio=2.0)
brk(_lab, "fincyl_order_finest", _fcc["orders"][-1])
brk(_lab, "fincyl_extrapolated_vs_series", abs(rel(_fcc["extrapolated"], fc_ref)))
brk(_lab, "fincyl_extrapolated_vs_series_all16", FC_ALL["cell"])
_cc = observed_orders(lambda n: CUBE_SOL[n][0].eta_flux(3.0, CUBE_SOL[n][1], read="cell"), CUBE_NS)
brk(_lab, "cube_order_8_16_32", _cc["orders"][-1])
brk(_lab, "cube_extrapolated_vs_series", abs(rel(_cc["extrapolated"], cube_ref)))
_p32, _x32 = CUBE_SOL[32]
brk(_lab, "eta_flux_vs_volume_worst", abs(_p32.eta_flux(3.0, _x32, read="cell") - _p32.eta(3.0, _x32)))
_lab = "finest ladder level built with the previous level's cells (refinement skipped)"
brk(_lab, "sphere_ladder_finest_diff", abs(rel(sphere_eta(SPH_NS[-2]), sphere_eta(SPH_NS[-2]))))
brk(_lab, "fincyl_ladder_finest_diff", abs(rel(fincyl_eta(FC_LAD_K[-2]), fincyl_eta(FC_LAD_K[-2]))))
brk(_lab, "cube_ladder_finest_diff", abs(rel(cube_lad["values"][-2], CUBE_SOL[16][0].eta(3.0))))

# ---- the finite cylinder, rod and cube
_lab = "finite cylinder assembled with nu = 0 on the radial axis"
_w0, _ = fincyl_gaps(nu_r=0)
brk(_lab, "fincyl_series_p_as_nomenclature", _w0["nomenclature"])
brk(_lab, "fincyl_series_p_as_running_text", _w0["running text"])
_rod1 = Pellet((200, 200), (1.0, 1.0), (1, 0))
brk("square rod assembled with nu = 1 on axis 0", "square_rod_vs_series",
    max(abs(rel(_rod1.eta(L / rod_va), eta_box_series(L, (1, 1), rod_va, 800))) for L in ROD_L))
brk("cube grid coarsened to n = 16", "cube_vs_series_n32", abs(rel(cube_lad["values"][1], cube_ref)))
_lam_bad = 20.0 / 1.0        # modulus taken on the half-width: v_p/s_x = a instead of a/3
brk("cube at Lambda = 20 with the modulus on the half-width (lambda a = 20, not 60)",
    "cube_n24_n32_gap_Lambda20", abs(rel(cube24.eta(_lam_bad), CUBE_SOL[32][0].eta(_lam_bad))))
_cube_a2 = CUBE_SOL[32][0].eta(1.0 / 0.5)
brk("cube taken at v_p/s_x = a/2 instead of a/3", "order_at_va_a3_Lambda1",
    order_label({**EQ_A3, "cube": _cube_a2})[0])

# ---- the sweep and the band
brk("sweep grids uniform (q = 1) for the flat disc, l/a = 0.25", "sweep_pymrm_vs_series_worst",
    max(_du[0], claims["sweep_pymrm_vs_series_worst"]))
brk("uniform disc grid with n_r halved to 60", "sweep_disc_uniform_grid_worst_pct",
    100 * disc_uniform_worst(n_r=60)[0])
_lab = "cube series built on v_p/s_x = a/2 instead of a/3"
_refs_b = {**refs, "cube": np.array([eta_cube_series(L, 0.5) for L in Lams])}
brk(_lab, "sweep_points_cube_lowest", int(np.sum(cube_lowest_margins(_refs_b) < 0)))
brk(_lab, "sweep_points_cyl2a_between_sphere_cube", int(np.sum(cyl2a_between_margins(_refs_b) > 0)))
_bb = widest(eta_plate, lambda L: eta_cube_series(L, 0.5))
brk(_lab, "band_all_shapes", _bb[0])
brk(_lab, "band_all_Lambda", _bb[1])
brk(_lab, "band_ratio_all_over_aris", _bb[0] / band_aris)
_sb = widest(eta_sphere, lambda L: eta_cube_series(L, 0.5))
brk(_lab, "cube_below_sphere", _sb[0])
brk(_lab, "cube_below_sphere_pct_there", 100 * _sb[0] / eta_sphere(_sb[1]))
_rb = minimize_scalar(lambda L: eta_cube_series(L, 0.5) / eta_sphere(L), bounds=(0.3, 20.0),
                      method="bounded", options={"xatol": 1e-6})
brk(_lab, "cube_below_sphere_max_pct", 100 * (1 - _rb.fun))
brk(_lab, "cube_below_sphere_max_pct_Lambda", _rb.x)
_lab = "sphere with the modulus on the radius (eq. 17 at Lambda / 3)"
_ab = widest(eta_plate, lambda L: eta_sphere(L / 3))
brk(_lab, "band_aris_plate_sphere_max", _ab[0])
brk(_lab, "band_aris_Lambda", _ab[1])
brk("plate with the modulus on the full thickness (eq. 15 at 2 Lambda)",
    "sweep_others_minus_plate_max", others_minus_plate({**refs, "flat plate": eta_plate(2 * Lams)}))
_lab = "series with half the terms"
brk(_lab, "fincyl_series_truncation_worst", trunc_fincyl(scale=2))
brk(_lab, "cube_series_truncation_Lambda1",
    abs(rel(eta_box_series(1.0, (1, 1, 1), CUBE_VA, 100), cube_ref)))
brk(_lab, "rod_series_truncation_worst",
    max(abs(rel(eta_box_series(L, (1, 1), 0.5, 400), eta_box_series(L, (1, 1), 0.5, 800)))
        for L in ROD_L))
brk(_lab, "sweep_series_truncation_Lambda20", trunc_sweep(scale=2))

# ---- eq. 18 and the mixture
brk("series with half the terms", "eq18_cyl_reference_truncation_Lambda50",
    abs(rel(eta_fincyl_axial(50.0, 0.5, 1.0, 3200), eta_fincyl_axial(50.0, 0.5, 1.0))))
brk("closed-form axial sum with p doubled", "fincyl_axial_vs_double_series",
    abs(rel(eta_fincyl_axial(1.0, 1.0, 1.0), eta_fincyl_series(1.0, 0.5, 1.0))))
_ext1, _g501 = asym_metrics(asym_pymrm(q=1.0))
brk("eq. 18 grids uniform (q = 1)", "eq18_extrapolated_worst_dev", _ext1)
brk("eq. 18 grids uniform (q = 1)", "eq18_pymrm_worst_rel_err_Lambda50", _g501)
brk("eq. 18 straight line fitted through Lambda = 5, 10, 20", "eq18_reference_extrapolated_worst_dev",
    max(abs(extrapolated(v, slice(0, 3)) - 1) for v in asym_ref.values()))
_m = mixture(mean=lambda f, L: float(np.sum(f * L)))
brk("mean modulus as the arithmetic mean of the sizes", "mixture_mean_modulus_dev",
    abs(_m[0] - PRINTED["mean_modulus"][0, 0]))
brk("mean modulus as the arithmetic mean of the sizes", "mixture_eta_at_mean_dev",
    abs(_m[1] - PRINTED["eta_at_mean_modulus"][0, 0]))
brk("eq. 19 averaged on the sphere curve", "mixture_eta_eq19", mixture(curve=eta_sphere)[2])

# %%
print(f"{'injected defect':<74} {'claim':<40} {'clean':>11} {'broken':>11}  caught")
for name, metric, clean, broken, caught in BREAK:
    is_bound = claims.metrics[metric]["kind"] == "bound"
    shown = f"(here {clean:.1e})" if is_bound else (f"{clean:11.4g}" if np.isscalar(clean) and not isinstance(clean, str) else str(clean)[:11])
    bshow = str(broken)[:11] if isinstance(broken, str) else f"{broken:11.4g}"
    print(f"{name[:74]:<74} {metric[:40]:<40} {shown:>11} {bshow:>11}  "
          f"{'yes' if caught else 'NO'}")
_caught = {m for _, m, _, _, c in BREAK if c}
_structural = sorted(k for k, c in claims.metrics.items() if c.get("structural"))
_open = sorted(k for k in claims.metrics if k not in _caught and k not in _structural)
print(f"\n{len(BREAK)} rows; {len(_caught)} of {len(claims.metrics)} claims are failed by at "
      f"least one row; structural: {', '.join(_structural) or 'none'}; "
      f"open: {', '.join(_open) or 'none'}")
_uncaught = [(n, m) for n, m, _, _, c in BREAK if not c]
print("rows that do not fail their claim: "
      + ("; ".join(f"{n} -> {m}" for n, m in _uncaught) or "none"))

# %% [markdown]
# **A break table measures sensitivity, never correctness.** Nothing in it can
# catch a wrong baseline; that is what the references are for. Which quantities
# have a second route, and which do not:
#
# - the 1-D pymrm solves: eqs. 15-17 and the shell expression (closed forms);
# - the 2-D finite cylinder: Aris's double series; the 2-D square rod and the 3-D
#   cube: their eigenfunction series (written on this page, sharing no code with
#   the solver);
# - the series themselves: their own truncation, by doubling the terms;
# - $\eta$ from the surface flux: **not** a second route. It is the volume average
#   rearranged by the conservation of the scheme, and their gap is the solve
#   residual;
# - the band's edges: the plate's closed form and the cube's series; the ordering
#   of the shapes is read from the references, not from the grids.

# %% [markdown]
# ## What pymrm adds
#
# **For Aris's three shapes: nothing.** Eqs. 15-17 are exact and take one line
# each. The 1-D pymrm solves here exist only to establish that the code producing
# the new results reproduces the known ones. [`B1.1`](../B1.1-thiele-weisz-hicks/)
# already plots the collapse for those three shapes.
#
# **What it adds is a solver for every other shape with one class.** Aris's
# method, solve $\nabla^2\chi = \lambda^2\chi$ with $\chi = 1$ on the boundary and
# average $\chi$, is completely general. Every shape on this page happens to be
# separable, so each also has a series (case (iv) as Aris printed it, unevaluated;
# the rod and the cube written here), and those series are the references. What
# the solver adds is that the same few lines handle every shape without deriving a
# series for each, and that its answer is independent of those derivations, which
# is what settles the reading of case (iv). The same pymrm operators (`construct_grad`,
# `construct_div` with `nu` on axis 0, `construct_coefficient_matrix` for the
# reaction) handle:
#
# - a **spherical shell**, `nu=2` on a radial domain starting at $r = pa$ rather
#   than 0, a variable-area 1-D problem with a sealed inner face;
# - an **axisymmetric finite cylinder**, `nu=1` radial plus Cartesian axial, which
#   is where the printed series' factor-of-two ambiguity gets resolved;
# - a **square rod** and a **cube**, which have corners and no closed form, only
#   eigenfunction series.
#
# `compute_boundary_values` reads the surface flux for the second form of $\eta$,
# and `pymrm.checks.observed_orders` and `find_roots` run the grid ladders and the
# root scans. The linear solves are scipy's `spsolve`; the band's maximum is found
# with scipy's `minimize_scalar`.
#
# **What it does not do here.** The band's edges are the plate's closed form and
# the cube's eigenfunction series, at every swept $\Lambda$; the 3-D solve's part is
# the check of that series at $\Lambda = 1$. The rod and the five finite cylinders
# in the sweep are pymrm solves, and they lie inside the band; for each of them a
# series also exists, which is what their accuracy is measured against.
#
# **The practical consequence.** Use $v_p/s_x$ and read $\eta$ off the plate
# curve, and your error is bounded by the band width, worst near $\Lambda = 1$,
# and one-sided for every shape here: the true $\eta$ is lower. Read it off the
# *sphere* curve believing that to be the conservative choice, and you can still
# be several percent optimistic for a pellet with corners, which is most
# extrudates.

# %% [markdown]
# ## Reuse
#
# **`Pellet` (in `model.py`) is the reusable piece.** Per-axis `n`, `extent` and
# `nu`, plus an `origin` for hollow shapes. Every geometry on this page is one
# call:
#
# ```python
# Pellet((400,), (1.0,), (2,))                  # sphere
# Pellet((400,), (1.0,), (2,), origin=0.5)      # spherical shell, p = 0.5
# Pellet((160, 160), (1.0, 1.0), (1, 0))        # finite cylinder, half-length = radius
# Pellet((32, 32, 32), (1.0,) * 3, (0, 0, 0))   # cube
# ```
#
# To put kinetics other than first order on it, replace the
# `lambda**2 * jac_react` term with a `NumJac((*shape, 1))` Jacobian of the rate
# and a Newton solve; the diffusion operator is unchanged. That is the step `B1.3`
# takes, and Bischoff's generalised modulus is built on this same $v_p/s_x$
# length.
#
# **The trap, and it is silent.** The low end of every axis carries
# `{"a": 1, "b": 0, "d": 0}`, which with the **outward** normal means
# $\partial\chi/\partial n = 0$: symmetry for a solid shape, a sealed face for
# the shell. Mark the shell's inner face `EXPOSED` by mistake and you get a
# perfectly convergent solution to a different problem, a shell open on both
# sides. Nothing warns you, and the $p \to 0$ limit does not catch it (the open
# shell also tends to the sphere there); the solve against the shell expression
# at the printed $p$ does, by the deviation printed in the shell section.
#
# **A note on cost and resolution.** The reacting layer is $O(1/\lambda)$ thick.
# On the uniform $24^3$ and $32^3$ grids the cube is not resolved at
# $\Lambda = 20$; for the 2-D shapes, grading toward the exposed faces
# (`q = 2`) keeps the sweep within the printed error of its series up to
# $\Lambda = 20$ on the same number of cells. Measure on a grid ladder before
# trusting a value where $\lambda h$ is not small.
#
# **Getting this paper.** Aris (1957) is closed access and the 1957 scan OCRs
# badly. Chemical Engineering Science reprinted it verbatim in 1995 as
# **50(24) 3899-3903**, doi:`10.1016/0009-2509(96)81819-7`, PII
# `0009250996818197`, and that re-typesetting has a clean text layer. Worth trying
# for any pre-1970 CES classic: the Golden Jubilee issue reprinted a dozen of
# them.
#
# **Related pages.** [`B1.1`](../B1.1-thiele-weisz-hicks/) (Thiele's exact
# solutions and the non-isothermal pellet); `B1.3` (Bischoff: the same length
# scale, arbitrary kinetics); `B1.4` (Weisz-Prater); `B1.12` (contour-averaged
# non-spherical pellets, deferred pending publication); `J4.7` (the same equation
# in an enzyme pellet).
#
# **Cite the source, not this page:** Aris, R., *On shape factors for irregular
# particles. I. The steady state problem. Diffusion and reaction*, Chemical
# Engineering Science **6**(6) 262-268 (1957),
# doi:10.1016/0009-2509(57)85028-3.

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
