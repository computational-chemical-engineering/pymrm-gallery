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
# title: "Sieverts-law hydrogen permeation through Pd membranes"
# description: "Why the exponent is one half, what it buys you, and when it stops mattering. The permeance Itoh quotes is derived rather than quoted, and the same numbers show his membrane was 18x larger than his duty, which is why his measurement cannot see it."
# categories: [sec:H, struct:S3, struct:S5, tier:T0, data:tier2, phase:gas]
# date: 2026-08-01
# ---
#
# # Sieverts-law hydrogen permeation through Pd membranes
#
# **Catalog ID:** `H1.1` · **Structures:** `S3` (1-D steady BVP: the wall), `S5`
# (1-D convection with a nonlinear sink: the channel) · **Tier:** T0
#
# Hydrogen crosses a dense palladium wall in a way nothing else does: it
# dissociates on the surface, dissolves as *atoms*, diffuses, and recombines on
# the far side. Because two atoms come from one molecule, the dissolved
# concentration follows the **square root** of the gas pressure, and the flux
# follows the difference of two square roots rather than the difference of two
# pressures. That single exponent is the whole content of Sieverts' law, and it is
# what every palladium membrane model, including the reactor on page
# [`H1.4`](../H1.4-itoh-membrane-dehydrogenation/), is built on.
#
# This page is about the closure, not about a reactor. It does four things:
#
# 1. **derives** the permeance that Itoh (1987) prints as a formula, by solving
#    the wall as a radial diffusion problem in pymrm, and lands on the value he
#    prints, the one number in the source anything here can be checked against;
# 2. shows what the exponent ½ *does*: it makes a permeator strip its feed
#    completely in a **finite** length against a clean permeate, which an exponent
#    of 1 never does, and it makes the required driving mole fraction fall as the
#    *square* of the membrane's spare capacity;
# 3. turns that into a criterion for **when permeation is the limiting
#    resistance**, and finds Itoh's membrane 18× from being it, which is the
#    condition that puts `H1.4`'s reactor on its fast-permeation ceiling and
#    leaves the permeance nothing to move;
# 4. shows what a departure from ½ means, by putting a finite surface step in
#    series with the wall and measuring the apparent exponent it produces.
#
# **Nothing on this page is compared with a permeation measurement, because the
# source does not print one.** The validation section says exactly what each check
# can and cannot fail on, and breaks every one of them on purpose.

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "gallery-2026.10.31.5"  # replaced by the release tag at publication
PYMRM = ">=3.2.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "H1.1-sieverts-permeation"  # equals the page directory

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
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from pymrm import NewtonBreakdown, continuation
from gallery_utils import load_data, load_meta, cite_data, Claims
from model import (sieverts_c, wall_solve, alpha_from_wall, Permeator, G_half, exact_half,
                   exact_quadrature, exact_lin, q_series, n_apparent, root_log)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# Palladium's selectivity for hydrogen is old news (Itoh calls it "a phenomenon
# that has been known for more than 100 years"), and the mechanism is a chain of
# five steps in series:
#
# 1. transport of H₂ through the gas to the surface;
# 2. dissociative adsorption, H₂ + 2\* ⇌ 2H\*;
# 3. dissolution of adsorbed H into the metal lattice;
# 4. Fickian diffusion of dissolved atomic H across the wall;
# 5. the reverse of 3, 2 and 1 on the permeate side.
#
# When steps 2 and 3 are fast on both faces, the dissolved concentration at each
# surface is set by *equilibrium* with the local gas. Write that equilibrium for
# the dissociative reaction ½H₂(g) ⇌ H(dissolved): its equilibrium constant is
#
# $$
# K = \frac{a_\mathrm{H}}{\sqrt{p_\mathrm{H_2}}}
# \qquad\Longrightarrow\qquad
# c = C_0\sqrt{p/P_0},
# $$
#
# which is **Sieverts' law**: the concentration of dissolved hydrogen goes as the
# square root of the gas pressure, with $C_0$ the concentration in equilibrium
# with a reference pressure $P_0$. Feed that as a boundary condition to step 4 and
# the flux comes out as a difference of square roots. The exponent ½ is therefore
# not empirical: it is a direct consequence of the molecule splitting in two.
#
# The exponent is also the diagnostic. If the *surface* steps are not fast (thin
# membranes, poisoned or contaminated surfaces, low temperature), then the
# dissolved concentration is no longer in equilibrium with the gas, and the flux
# stops following a square root. In the extreme where the surface controls, the
# flux becomes proportional to the pressure *difference* itself, exponent 1.
# Anything between the two says the resistances are comparable. Section 5 builds
# exactly that series model and measures the apparent exponent it produces.

# %% [markdown]
# ## The published model
#
# ### The source, and what it actually carries
#
# The permeation closure is read from **Itoh, N. (1987), "A membrane reactor using
# palladium", *AIChE Journal* 33(9) 1576-1578**
# ([doi:10.1002/aic.690330921](https://doi.org/10.1002/aic.690330921)), the same
# 3-page communication that page `H1.4` reproduces. That paper is on disk and was
# read; the classics it cites were not, so this page states carefully what it
# takes from where.
#
# **Itoh names the law and uses it.** His sentence, from the Analysis section:
#
# > The permeation rate of hydrogen gas through the palladium membrane, $Q_H$, was
# > assumed to obey the **half-power pressure law** (Bohmholdt and Wicke, 1967),
# > i.e., it is proportional to the difference between **the roots** of the
# > hydrogen partial pressures in the reaction and separation sides, $p_H$ and
# > $p'_H$:
#
# and then Eq. 1, printed exactly as:
#
# $$
# Q_H \;=\; \alpha_H\left(\sqrt{\frac{p_H}{P_o}} - \sqrt{\frac{p'_H}{P_o}}\right),
# \qquad
# \alpha_H \;=\; \frac{2\pi l_o}{\ln(r_o/r_i)}\,D\,C_o
# \tag{1}
# $$
#
# followed by: *"$P_o$ is 1.013 × 10⁵ Pa and $D$ is Fick's diffusion coefficient
# of hydrogen dissolved in the palladium. $C_o$ is defined as a standard
# concentration of hydrogen in the palladium, which is at equilibrium with the
# gas-phase hydrogen at pressure $P_o$. The $\alpha_H$ is defined as a permeation
# rate constant of hydrogen gas. When literature values, i.e., 9.23 × 10⁻¹⁰ m²
# s⁻¹ for $D$ (Nagamoto and Inoue, 1985) and 1,280 mol m⁻³ for $C_o$ (Sieverts and
# Danz, 1936) at 473 K, were used, $\alpha_H$ calculated was 4.47 × 10⁻⁵
# mol s⁻¹."*
#
# **Three things follow, and each constrains what this page may claim.**
#
# - **The exponent is attributed to Bohmholdt and Wicke (1967), not to Sieverts.**
#   Sieverts and Danz (1936) is cited by Itoh for exactly one number: the
#   solubility constant $C_o = 1{,}280$ mol/m³ at 473 K, the amplitude of
#   $c = C_0\sqrt{p/P_0}$. The catalog entry for this case reads "Sieverts;
#   exponent analysis by Caravella et al. (2010)". Caravella is not on disk and is
#   **not** used here; the half of that citation this page can honour is Itoh's
#   own, and the page is written to it.
# - **Itoh prints no permeation measurement.** No flux, no permeance measured on
#   this membrane, no pressure sweep. $\alpha_H$ is *computed* from two literature
#   constants, and that computed value is the one number in the source against
#   which anything here can be checked.
# - **Itoh prints no temperature dependence.** $D$ and $C_o$ are quoted at 473 K
#   only, and $\alpha_H$ is a single number. Nothing about the activation energy
#   of the permeance can be said from this source, and this page says nothing
#   about it.
#
# **How the paper was read.** The PDF's page images are 300 dpi bilevel CCITT
# scans, so `pdftoppm -r 600` is a 2× upsample carrying no extra information;
# every character quoted above was read on the *native* bitmap. All three marginal
# superscripts in the constants sentence (−10, −3, −5) read cleanly, and the −10
# and −3 on the same line serve as unambiguous comparison glyphs for the −5.
# **No misprint was found anywhere in the permeation law.** The two printing
# problems in this paper (Eq. 5's sign and the $K_P$ expression) are elsewhere,
# and are resolved on page `H1.4`.
#
# ### Where Eq. 1 comes from
#
# Itoh states Eq. 1; he does not derive it. It is the exact steady solution of
# radial diffusion through a cylindrical shell with Sieverts equilibrium at each
# face. Inside the wall, dissolved atomic hydrogen obeys
#
# $$
# \frac{1}{r}\frac{\mathrm{d}}{\mathrm{d}r}\!\left(r\,D\,\frac{\mathrm{d}c}{\mathrm{d}r}\right) = 0,
# \qquad
# c(r_i) = C_0\sqrt{p_H/P_0},\quad c(r_o) = C_0\sqrt{p'_H/P_0},
# $$
#
# whose solution is **logarithmic in $r$** (not linear), and whose total molar
# throughput over a tube of length $l_0$ is
#
# $$
# Q_H = \frac{2\pi l_0 D}{\ln(r_o/r_i)}\bigl(c(r_i) - c(r_o)\bigr)
#     = \underbrace{\frac{2\pi l_0 D C_0}{\ln(r_o/r_i)}}_{\alpha_H}
#       \left(\sqrt{p_H/P_0} - \sqrt{p'_H/P_0}\right).
# $$
#
# Section 1 solves that BVP in pymrm rather than quoting the closed form, and
# compares the resulting $\alpha_H$ with Itoh's printed 4.47 × 10⁻⁵ mol/s. Two
# things are worth noticing before the numbers arrive: the geometry enters only
# through $\ln(r_o/r_i)$, and **the square roots enter only through the boundary
# conditions**. The exponent lives entirely in steps 2 and 3 of the mechanism; the
# diffusion step is linear.

# %% [markdown]
# ## Parameters and assumptions
#
# | Quantity | Value | Source |
# |---|---|---|
# | $T$ | 473 K | stated |
# | $P_0$ | 1.013 × 10⁵ Pa | stated, below Eq. 1 |
# | $D$ (H in Pd, 473 K) | 9.23 × 10⁻¹⁰ m²/s | stated; Itoh cites Nagamoto & Inoue (1985) |
# | $C_0$ (473 K) | 1280 mol/m³ | stated; Itoh cites **Sieverts & Danz (1936)** |
# | $l_0$ | 0.14 m | stated ("140 mm long") |
# | wall thickness | 200 µm | stated |
# | $r_i,\;r_o$ | 8.5 mm, 8.7 mm | **reconstructed**, see below |
# | $\alpha_H$ | 4.47 × 10⁻⁵ mol/s | stated; *computed by Itoh from Eq. 1* |
# | exponent $n$ | ½ | stated ("half-power pressure law") |
#
# **The radii are a reconstruction, and this page depends on it, so it rests on
# two independent readings, not one.**
#
# The running text says the palladium tube was "200 µm thick, **17.0 mm OD**, and
# 140 mm long". Read that way, $\alpha_H$ deviates by -2.4 % from the printed
# value; reading 17.0 mm as the **ID** gives -0.03 %. That comparison
# **discriminates between two candidate readings of one printed sentence**, but
# it is *not* an independent test of the geometry, because the radii it selects
# are its own output. Section 1a and validation check 1 say so plainly.
#
# The independent witness is **Figure 1**, which dimensions the same tube. It
# prints four dimension *values* (`0.2`, `17ϕ`, `28ϕ` and `140`) with the unit
# `[mm]`; those are printed numbers, read on the native bitmap like any table
# entry. *Which* diameter `17ϕ` labels depends on where its arrowheads terminate,
# and that is an **interpretation of the drawing**, so it is set out here rather
# than asserted. Read on the 300 dpi bitmap of journal page 1577 (page-image rows
# quoted; upper membrane line = rows 262-266, lower = 344-348, upper shell band =
# 231-238, lower = 369-377):
#
# | dimension | arrow column | upper tip | lower tip | what it terminates on |
# |---|---|---|---|---|
# | `0.2` | 705 | row 259, pointing **down** | row 267, pointing **up** | the two tips face *each other* across the membrane band from outside: the **wall thickness** |
# | `17ϕ` | 754 | row 268, pointing **up** | row 342, pointing **down** | the tips face *apart* and land on the **inner** faces of the two membrane lines: the **bore** |
# | `28ϕ` | 809 | line starts at row 239, under the shell band | row 368, pointing **down** | the **inner** face of the lower shell wall; the lower tip is crisp, the upper arrowhead's point is the least clean of the six (narrowest pixel row 255, not 239), so this row is context, not evidence |
#
# **Read that way, Figure 1 dimensions the palladium tube's *inner* diameter as
# 17 mm** ($r_i = 8.5$ mm, $r_o = 8.7$ mm) independently of $\alpha_H$, and in
# contradiction to the text's "17.0 mm OD". Only the arrowhead *placement* is
# used: the drawing is schematic, with the 0.2 mm wall drawn about 5 px where
# 17 mm is 74 px, so no length is measured off it. The wall is dimensioned
# *separately* as `0.2`, which is the internal reason `17` cannot include it.
#
# **And $\alpha_H$ fixes only the ratio $r_o/r_i$.** At the printed value that
# ratio is 1.023522, which $(r_i,\text{wall}) = (4.25\ \text{mm}, 100\ \mu m)$,
# $(8.5\ \text{mm}, 200\ \mu m)$ and $(17\ \text{mm}, 400\ \mu m)$ all satisfy
# identically. **The printed 200 µm wall is load-bearing**: it is what converts a
# ratio into a radius. Section 1a puts numbers on all of this.
#
# **Assumptions.** Isothermal; the wall is uniform, dense and defect-free; $D$ is
# independent of hydrogen content (a real approximation, since $D$ in Pd does
# depend on loading, but it is Itoh's, and no loading dependence is printed anywhere
# in the source); surface steps are fast enough that Sieverts equilibrium holds at
# both faces (section 5 is precisely the study of relaxing this); and the channel
# model of section 2 assumes isothermal, isobaric plug flow with the permeate-side
# hydrogen mole fraction prescribed.

# %% [markdown]
# ## The data
#
# Everything Itoh prints about permeation, and nothing else. There is no table in
# the paper; the values live in Eq. 1 and in the sentence beneath it, and all of
# them were transcribed from the native page bitmaps.
#
# **Data tier.** The values themselves are tier 2: numbers printed in a paper.
# The *validation* situation is tier 6: the source contains no permeation
# measurement, so the only external target on this page is $\alpha_H$, a number
# Itoh **computed** rather than measured. Every other check below is an internal
# identity or a closed-form limit. Nothing here is validated against experiment,
# and the page never says it is.

# %%
vals = load_data("itoh-1987-permeation.csv", page=PAGE)
meta = load_meta("itoh-1987-permeation.csv", page=PAGE)
print(vals.to_string(index=False, columns=["quantity", "value", "unit", "kind"]))
print(f"\n{cite_data(meta)}")


def stated(offset=0):
    """The data file's values by name. ``offset`` reads each value that many rows
    further down the file (a break-row handle: a misaligned read)."""
    v = vals.value.to_numpy(dtype=float)
    return {q: v[(i + offset) % len(v)] for i, q in enumerate(vals.quantity)}


S = stated()
T        = S["temperature"]              # K
P0       = S["P_0"]                      # Pa
D_PD     = S["D_H_in_Pd"]                # m2/s, H in Pd at 473 K
C_0      = S["C_0"]                      # mol/m3, Sieverts reference concentration
ALPHA_H  = S["alpha_H"]                  # mol/s, printed permeation rate constant
L0       = S["l_0"]                      # m, tube length
DELTA    = S["wall_thickness"]           # m
RI       = S["r_i"]                      # m, reconstructed
RO       = S["r_o"]                      # m, reconstructed
N_EXP    = S["exponent_n"]               # the printed exponent, 1/2
UC0      = S["u_C0"]                     # mol/s, cyclohexane feed of the reported run
YC0      = S["feed_cyclohexane_fraction"]
VA0      = S["v_A0"]                     # mol/s, argon purge of the reported run
X_MEAS   = S["measured_conversion"]      # the conversion Itoh reports for that run
RI_OD, RO_OD = RI - DELTA, RO - DELTA    # m, the radii if "17.0 mm" is the outer diameter
assert N_EXP == 0.5, "this page implements the half-power law that the source states"
print(f"\nT = {T:.0f} K, P_0 = {P0:.3e} Pa, printed exponent n = {N_EXP}")

# %% [markdown]
# ## PyMRM implementation
#
# Two models, deliberately kept separate because they answer different questions.
# Both live in this page's `model.py`.
#
# ### 1. The wall: a 1-D steady BVP (`S3`)
#
# Radial diffusion of dissolved atomic hydrogen across the cylindrical shell,
# $r \in [r_i, r_o]$, with Dirichlet values fixed by Sieverts equilibrium at each
# face (`wall_solve`). `construct_grad` builds the gradient (with the boundary
# contribution in `grad_bc`), `construct_coefficient_matrix` puts $D$ on the
# faces, and `construct_div` builds the divergence with **`nu=1`**: $r$ is a
# radial coordinate, so the cell volumes carry the $r$ weight. The system is
# linear, so one `spsolve` finishes it; no Newton is needed. A wall whose $D$
# varies with radius takes $D$ evaluated at the face positions; a jump in $D$
# between cells would need the harmonic face mean instead.
#
# The permeance is then *read off the solution*, not evaluated from a formula: the
# molar rate through the cylindrical face at radius $r$ is
# $Q(r) = 2\pi r\, l_0 \times (-D\,\partial c/\partial r)$, which the steady
# solution makes independent of $r$ (`alpha_from_wall`).
#
# Boundary conditions use the **outward normal** ($a\,\partial c/\partial n + b\,c
# = d$). Both ends are Dirichlet here ($a{=}0$, $b{=}1$, $d = c_{i}$ or $c_{o}$),
# so the sign convention does not bite, but the flux sign does: the face flux is
# $-D\,\partial c/\partial r$ in the $+r$ direction, so an outward (positive)
# hydrogen flow needs $c_i > c_o$.

# %%
def alpha_wall(r_i=None, r_o=None, D=None, l_0=None, n_r=200, nu=1):
    """alpha_H DERIVED from the pymrm wall solve, with this page's constants."""
    return alpha_from_wall(RI if r_i is None else r_i, RO if r_o is None else r_o,
                           D_PD if D is None else D, C_0, L0 if l_0 is None else l_0,
                           n_r=n_r, nu=nu)


def alpha_closed(r_i=None, r_o=None, D=None):
    """Eq. 1's closed form, for comparison only."""
    r_i, r_o = (RI if r_i is None else r_i), (RO if r_o is None else r_o)
    return 2 * np.pi * L0 * (D_PD if D is None else D) * C_0 / np.log(r_o / r_i)

# %% [markdown]
# ### 1a. The wall solution, and Itoh's permeance
#
# The profile is logarithmic, the flux is constant across the wall, and the
# permeance that comes out is Eq. 1's $\alpha_H$: derived, not quoted.

# %%
def wall_numbers(D=None, delta=None, nu=1):
    """Everything section 1a states about the wall, for given constants."""
    D = D_PD if D is None else D
    delta = DELTA if delta is None else delta
    ro = RI + delta
    r_c, c_w, r_f, Q_f = wall_solve(sieverts_c(1.0, C_0), 0.0, RI, ro, D, L0, n_r=200, nu=nu)
    a_fv, a_cf = Q_f[0], alpha_closed(RI, ro, D)
    c_log = C_0 * (1.0 - np.log(r_c / RI) / np.log(ro / RI))
    ratio = np.exp(2 * np.pi * L0 * D * C_0 / ALPHA_H)

    def dev_of_ri(r):
        return (alpha_from_wall(r, r + delta, D, C_0, L0, nu=nu) - ALPHA_H) / ALPHA_H
    return dict(
        r_c=r_c, c_w=c_w, c_log=c_log, Q_f=Q_f, alpha_fv=a_fv,
        fv_vs_closed=(a_fv - a_cf) / a_cf,
        dev_printed=(a_fv - ALPHA_H) / ALPHA_H,
        profile_vs_log=np.max(np.abs(c_w - c_log)) / C_0,
        flux_spread=np.ptp(Q_f) / np.mean(Q_f),
        dev_od=(alpha_from_wall(RI - delta, RI, D, C_0, L0, nu=nu) - ALPHA_H) / ALPHA_H,
        ratio_locus=ratio,
        # alpha_H rises strictly with r_i at a fixed wall, so each root is unique
        ri_lo=root_log(lambda r: dev_of_ri(r) + 0.01, 5e-3, RI),
        ri_hi=root_log(lambda r: dev_of_ri(r) - 0.01, RI, 20e-3),
        thin_10um=(alpha_from_wall(RI, ro - 10e-6, D, C_0, L0, nu=nu) - ALPHA_H) / ALPHA_H,
        thick_10um=(alpha_from_wall(RI, ro + 10e-6, D, C_0, L0, nu=nu) - ALPHA_H) / ALPHA_H)


W = wall_numbers()
ALPHA_FV, ALPHA_CF = W["alpha_fv"], alpha_closed()
print("wall solve at unit driving force (p_H = P_0, p'_H = 0), n_r = 200")
print(f"  alpha_H from the pymrm solve   : {ALPHA_FV:.6e} mol/s")
print(f"  alpha_H closed form, Eq. 1     : {ALPHA_CF:.6e} mol/s "
      f"({W['fv_vs_closed']:+.1e} relative)")
print(f"  alpha_H PRINTED by Itoh        : {ALPHA_H:.3e} mol/s "
      f"({W['dev_printed'] * 100:+.3f} % from the solve)")
print(f"  profile vs the analytic log    : max |dc|/C_0 below 1e-10 (here {W['profile_vs_log']:.0e})")
print(f"  flux across the {len(W['Q_f'])} faces      : relative spread "
      f"below 1e-10 (here {W['flux_spread']:.0e}); div(J) = 0 is what was solved")

print("\nwhich of the two readings of '17.0 mm' did Itoh substitute into Eq. 1?")
print(f"  17.0 mm read as ID (r_i = {RI * 1e3:.1f}, r_o = {RO * 1e3:.1f} mm): "
      f"{W['dev_printed'] * 100:+.2f} % vs printed")
print(f"  17.0 mm read as OD (r_i = {RI_OD * 1e3:.1f}, r_o = {RO_OD * 1e3:.1f} mm): "
      f"{W['dev_od'] * 100:+.2f} % vs printed")
print("  -> the ID reading, decisively. But read that for what it is: it PICKS")
print("     BETWEEN TWO READINGS of one printed sentence. It is not a test of the")
print("     geometry, because the geometry is its own output. The independent")
print("     witness is Figure 1's dimension arrows (see Parameters above), which")
print("     land on the INNER faces of the membrane lines and say the same thing.")

print("\nhow much geometry does the alpha_H identity actually pin down?")
print(f"  Eq. 1 constrains only the RATIO: alpha_H = {ALPHA_H:.3e} <=> "
      f"r_o/r_i = {W['ratio_locus']:.6f} exactly.")
print("\n".join(f"    (r_i, wall) = ({ri * 1e3:>5.2f} mm, {wall * 1e6:>3.0f} um) -> "
                f"alpha_H = {alpha_wall(ri, ri + wall):.4e} mol/s"
                for ri, wall in ((RI / 2, DELTA / 2), (RI, DELTA), (2 * RI, 2 * DELTA))))
print("  All three sit on the same locus. Only the PRINTED 200 um wall turns the")
print("  ratio into a radius, so that number is load-bearing.")
print(f"  Given that wall, agreement degrades to 1 % at r_i = {W['ri_lo'] * 1e3:.3f} mm "
      f"and {W['ri_hi'] * 1e3:.3f} mm,")
print(f"  a window of {(W['ri_lo'] - RI) * 1e6:+.0f}/{(W['ri_hi'] - RI) * 1e6:+.0f} um, i.e. "
      f"+-{(W['ri_hi'] - W['ri_lo']) / (2 * RI) * 100:.1f} % in r_i. Thickness bites harder:")
print("\n".join(f"    wall {(DELTA + dw) * 1e6:.0f} um at r_i = {RI * 1e3:.1f} mm -> "
                f"{W[key] * 100:+.2f} % vs printed"
                for key, dw in (("thin_10um", -10e-6), ("thick_10um", +10e-6))))

fig, ax = plt.subplots(1, 2, figsize=(11.2, 3.9))
xd = (W["r_c"] - RI) / (RO - RI)
ax[0].plot(xd, W["c_w"] / C_0, lw=2.5, label="pymrm solve")
ax[0].plot(xd, W["c_log"] / C_0, "k--", lw=1.2, label="analytic, log in $r$")
ax[0].plot(xd, 1 - xd, ":", color="tab:red", lw=1.6, label="a linear (slab) profile")
ax[0].set(xlabel="fractional depth into the wall, $(r - r_i)/(r_o - r_i)$", ylabel="$c/C_0$",
          title="dissolved H across the wall: log and linear are\n"
                "indistinguishable at this thickness")
ax[0].legend(fontsize=8)
pr = np.linspace(0, 1, 200)
ax[1].plot(pr, np.sqrt(pr), lw=2.5, label=r"Sieverts, $n=\frac{1}{2}$")
ax[1].plot(pr, pr, ":", color="tab:red", lw=1.8, label="$n = 1$")
ax[1].set(xlabel="$p_H/P_0$", ylabel=r"driving term $(p_H/P_0)^n$",
          title="the exponent, drawn")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()
print(f"The wall is {DELTA / RI * 100:.1f} % of the tube radius, so the log profile is barely")
print("distinguishable from a straight line, which is exactly why the alpha_H")
print("identity has almost no power against a wrong geometry. Validation check 2")
print("measures that rather than assuming it.")

# %% [markdown]
# ### 2. The channel: 1-D convection with a nonlinear sink (`S5`)
#
# A tube of the same membrane carrying hydrogen diluted in an inert, with the
# permeate side held at a fixed hydrogen mole fraction $y'$ (the large-purge or
# vacuum limit; the *coupled* two-channel problem, where the permeate composition
# is itself an unknown, is page `H1.4`). Writing $w = u_H/u_H^0$ for the hydrogen
# flow scaled on its inlet value, $\beta = u_{\text{inert}}/u_H^0$, and $L$ for
# dimensionless axial position,
#
# $$
# \frac{\mathrm{d}w}{\mathrm{d}L}
# = -\,\Pi\left[\left(\frac{w}{w+\beta}\right)^{n} - (y')^{n}\right],
# \qquad w(0) = 1,
# \qquad
# \boxed{\;\Pi \equiv \frac{\alpha_H}{u_H^0}\sqrt{\frac{P_T}{P_0}}\;}
# $$
#
# with $n = \tfrac12$ the printed exponent. **$\Pi$ is the whole problem**: the
# membrane's permeation capacity at full driving force, divided by the hydrogen
# throughput it has to handle. Everything in sections 3 and 4 is a statement about
# $\Pi$.
#
# In pymrm this is `Permeator`: `construct_convflux_upwind` with $v = 1$ (the
# equation is already in flow-rate form) and `construct_div` with **`nu=0`** ($L$
# is axial, not radial), assembled once, plus the pointwise sink linearised with
# `NumJac` on the state shape `(n_z, 1)` (**not** `(n_z,)`: a bare 1-D shape makes
# `NumJac`'s default stencil declare every cell coupled to every other and build a
# dense Jacobian).
#
# Newton does not iterate on $w$ but on $s = \sqrt{w}$. Substituting $w = s^2$
# into the same cell equations changes nothing in the discrete solution, and turns
# the half-power sink into
#
# $$
# \Pi\left(\frac{w}{w+\beta}\right)^{1/2} = \Pi\,\frac{s}{\sqrt{s^2+\beta}},
# $$
#
# which is smooth and linear through $s = 0$, where the sink in $w$ has an
# infinite slope. The Jacobian is the constant convection operator times
# $\mathrm{diag}(2s)$, minus `NumJac` of the composed sink (`Permeator.residual_s`).
# `pymrm.newton` solves it from the flat state $s = 1$ with `damping="monotone"`,
# the lower bound $s \ge 0$ and a step tolerance of 1e-10 on the order-one $s$;
# $w = s^2$ is then non-negative by construction. Section 4a shows why the
# change of variable is needed.
#
# Boundary conditions, outward normal: inlet $w = 1$ (Dirichlet, $a{=}0$,
# $b{=}1$, $d{=}1$); the outlet is convection only, so nothing is imposed there and
# it takes the library's pure-outflow marker `{"outflow": True}`. First-order
# upwind with that outlet is exactly backward Euler in $L$, and the exit face
# carries the last cell's value; outlet values are read with
# `compute_boundary_values`.

# %%
NEWTON_ITERS = []   # Newton iterations of every channel solve on the page that converged


def channel(Pi, beta=None, **kw):
    """A Permeator solved in s from the flat state; raises when Newton does not converge."""
    p = Permeator(Pi, BETA if beta is None else beta, **kw)
    p.solve()
    NEWTON_ITERS.append(p.result.nit)
    return p

# %% [markdown]
# ## Results
#
# ### 3. What the exponent ½ actually does
#
# Two consequences, both exact, both easy to check, and both invisible if the law
# is only ever written down.
#
# **(a) A ½-law strips its feed completely in a finite length, *if the permeate
# side is clean*.** Everything in this subsection is derived at $y' = 0$, and that
# condition is not cosmetic: with any permeate loading the tube strips only to the
# pressure-matched floor $\beta y'/(1-y')$, the approach to it is asymptotic, and
# the finite extinction length disappears entirely. Validation check 5b measures
# exactly that. With $y' = 0$, separating the equation gives $G(w_0) - G(w) = \Pi L$ with
# $G(w) = \sqrt{w(w+\beta)} + \beta\,\operatorname{arcsinh}\sqrt{w/\beta}$, and
# $G(0) = 0$. So $w$ reaches **exactly zero** at the finite length
#
# $$
# L_\mathrm{ext} = \frac{G(w_0)}{\Pi},
# $$
#
# because $\int_0 \mathrm{d}w/\sqrt{w}$ converges. With $n = 1$ the same integral
# is $\int_0 \mathrm{d}w/w$, which diverges: an exponent of 1 gives exponential
# decay that never terminates. The capacity at which a unit-length tube *just*
# achieves complete stripping is therefore $\Pi_\mathrm{crit} = G(1,\beta)$, a
# number, not a fit.
#
# **(b) Spare capacity is worth its square.** In steady operation the membrane
# must export a duty $q$, so the driving mole fraction settles at
# $y_H = (q/\alpha_H)^{1/n}$. With $n = \tfrac12$ that is $\Pi^{-2}$; with $n = 1$
# it is $\Pi^{-1}$. A membrane with 50× spare capacity needs a driving mole
# fraction 2500× smaller under a ½-law and only 50× smaller under a linear law.
# That is the whole reason a palladium membrane can hold a reactor's hydrogen
# partial pressure down at the $10^{-4}$ level, **against a clean permeate**. If
# the permeate is loaded, as in Itoh's own argon-purged reactor, what the reaction
# side reaches is set by the permeate composition instead, and the driving
# difference goes to zero rather than to $\Pi^{-2}$. The cell below computes both
# numbers for his tube, because they differ by a factor of 19 and are easy to
# confuse.

# %%
def itoh_numbers(nu_H2=3.0, y_C0=None, n=None):
    """Section 3's numbers for Itoh's tube from printed values only.

    ``nu_H2`` (hydrogen per cyclohexane), ``y_C0`` (feed fraction) and ``n``
    (exponent) are arguments so that the break table can change them.
    """
    y_C0 = YC0 if y_C0 is None else y_C0
    n = N_EXP if n is None else n
    n_ar = (1.0 - y_C0) / y_C0
    u_h0 = nu_H2 * UC0                        # mol/s of H2 to export
    beta = (UC0 * n_ar + UC0) / u_h0          # Ar + benzene over H2
    pi = ALPHA_H / u_h0                       # P_T = P_0, so the sqrt factor is 1
    pi_crit = G_half(1.0, beta)
    y_req = (1.0 / pi) ** (1.0 / n)
    q_both = nu_H2 * X_MEAS * UC0 / (UC0 * (1.0 + n_ar) + VA0)
    y_both = q_both / (1.0 + q_both)
    return dict(u_h0=u_h0, beta=beta, Pi=pi, Pi_crit=pi_crit, L_ext=pi_crit / pi,
                spare=pi / pi_crit, y_req=y_req, y_lin=1.0 / pi, y_both=y_both,
                both_over_req=y_both / y_req, q_both=q_both,
                square_example=((1.0 / pi) / (1.0 / (50 * pi))) ** (1.0 / n))


IT = itoh_numbers()
BETA, PI_ITOH, PI_CRIT = IT["beta"], IT["Pi"], IT["Pi_crit"]
print("Itoh's membrane, USED AS A SEPARATOR ON A CLEAN PERMEATE (y' = 0), for the")
print("duty his reaction sets. This is a hypothetical duty placed on his tube, NOT")
print("a description of his reactor: see the loaded-permeate block below.")
print(f"  hydrogen to export  u_H0 = 3 u_C0 = {IT['u_h0']:.3e} mol/s")
print(f"  inert carried along        beta  = {BETA:.3f}")
print(f"  permeation capacity ratio  Pi    = alpha_H / u_H0 = {PI_ITOH:.1f}")
print(f"  critical capacity for L_ext = 1  : Pi_crit = G(1, beta) = {PI_CRIT:.3f}")
print(f"  -> the hydrogen would be exhausted at L_ext = {IT['L_ext']:.4f} of the tube,")
print(f"     i.e. the membrane is {IT['spare']:.1f}x larger than the duty needs.")
print(f"  driving mole fraction that duty needs, n = 1/2 : {IT['y_req']:.3e}")
print(f"                                         n = 1   : {IT['y_lin']:.3e}"
      f"   ({IT['y_lin'] / IT['y_req']:.0f}x larger)")
print(f"  halving alpha_H moves the n = 1/2 requirement only to "
      f"{(1 / (0.5 * PI_ITOH)) ** 2:.2e}.")
print(f"  50x spare capacity lowers the n = 1/2 requirement {IT['square_example']:.0f}x.")

print("\nItoh's REACTOR is not that separator: his permeate side is LOADED.")
print(f"  argon purge v_A0 = {VA0:.2e} mol/s, measured conversion X = {X_MEAS}")
print("  at the co-current fast-permeation limit both sides carry the SAME H2")
print(f"  mole fraction  y = q/(1+q),  q = 3 X u_C0/(F_0 + v_A0) = {IT['q_both']:.4e}")
print(f"  -> y = {IT['y_both']:.3e}, which is {IT['both_over_req']:.0f}x the {IT['y_req']:.2e}")
print("     the clean-permeate separator would need, and the DRIVING difference")
print("     tends to zero, because that is what the ceiling means. Nor is the")
print(f"     hydrogen stripped in {IT['L_ext'] * 100:.1f} % of the length: with y' > 0 there is no")
print("     finite extinction length at all (check 5b).")
print(f"  What survives, and it is the useful part: Pi = {PI_ITOH:.1f} vs Pi_crit = {PI_CRIT:.2f}")
print(f"  means permeation is nowhere near limiting ({IT['spare']:.1f}x spare),")
print("  which is the condition that puts H1.4's model on its co-current")
print("  fast-permeation ceiling, and once there the conversion is a function of")
print("  K_p and the purge split alone, so alpha_H has no observable left to move.")

# %%
# --- profiles: the finite extinction of a 1/2-law vs the tail of a 1-law ----
# Compared at equal flux AT THE INLET, so the exponent is the only difference:
#   alpha_n y0^n = alpha_H y0^(1/2)   =>   Pi_n = Pi y0^(1/2 - n)
Y0 = 1.0 / (1.0 + BETA)
PI_DEMO = 0.90 * PI_CRIT                     # L_ext = 1.111: just past the outlet


def demo_profiles(pi_lin=None):
    """Outlet hydrogen of the 1/2-law and the equal-inlet-flux 1-law channels."""
    pi_lin = PI_DEMO * Y0 ** (0.5 - 1.0) if pi_lin is None else pi_lin
    p_half = channel(PI_DEMO, n=N_EXP, n_z=1600)
    p_lin = channel(pi_lin, n=1.0, n_z=1600)
    return p_half, p_lin, pi_lin


p_half, p_lin, PI_LIN = demo_profiles()
Lg = np.linspace(0.0, 1.6, 200)
ex_half, ex_lin = exact_half(Lg, BETA, PI_DEMO), exact_lin(Lg, BETA, PI_LIN)
L_EXT_DEMO = G_half(1.0, BETA) / PI_DEMO

fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.0))
for a, logy in zip(ax, (False, True)):
    a.plot(Lg, ex_half, lw=2.4, color="tab:blue", label=r"exact, $n=\frac{1}{2}$")
    a.plot(Lg, ex_lin, lw=2.0, color="tab:red", label="exact, $n=1$ (same inlet flux)")
    a.plot(p_half.z_c[::40], p_half.w.ravel()[::40], "o", ms=4, color="tab:blue",
           mfc="w", label="pymrm (cell values)")
    a.plot(p_lin.z_c[::40], p_lin.w.ravel()[::40], "s", ms=4, color="tab:red", mfc="w")
    a.axvline(L_EXT_DEMO, color="0.5", ls=":", lw=1.4)
    a.axvline(1.0, color="0.75", lw=1.0)
    a.set(xlabel="dimensionless length $L$", ylabel="$w = u_H/u_H^0$")
    if logy:
        a.set_yscale("log"); a.set_ylim(1e-12, 3)
        a.set_title(r"same, logarithmic: only the $\frac{1}{2}$-law lands")
    else:
        a.set_title("hydrogen left in the tube")
    a.legend(fontsize=8)
ax[0].text(L_EXT_DEMO + 0.02, 0.55, "$L_{ext}$", fontsize=9, color="0.4")
ax[0].text(1.015, 0.5, "tube outlet", fontsize=8, color="0.55", rotation=90)
fig.tight_layout(); plt.show()
W_HALF_OUT, W_LIN_OUT = p_half.outlet_value(), p_lin.outlet_value()
print(f"Pi = {PI_DEMO:.3f} (= 0.90 Pi_crit), so the exact extinction length is")
print(f"L_ext = {L_EXT_DEMO:.4f}, just past the outlet. Hydrogen left at L = 1 (exit face):")
print(f"  n = 1/2 : pymrm {W_HALF_OUT:.6e}, exact {exact_half(1.0, BETA, PI_DEMO)[0]:.6e}")
print(f"  n = 1   : pymrm {W_LIN_OUT:.6e}, exact {exact_lin(1.0, BETA, PI_LIN)[0]:.6e}"
      f"   ({W_LIN_OUT / W_HALF_OUT:.1f}x more)")
print("Both laws are normalised to the SAME flux at the inlet, so the exponent is")
print("the only difference between them. As the hydrogen thins out, sqrt(y) > y,")
print("so the 1/2-law keeps pulling where the linear law fades, and it reaches")
print(f"zero at the finite L_ext = {L_EXT_DEMO:.4f}, while the n = 1 solution obeys")
print("w + beta ln w = 1 - Pi L, whose logarithm sends w -> 0 only as L -> inf.")

# %% [markdown]
# ### 4. When is permeation the limiting resistance?
#
# $\Pi$ answers this on its own, and the answer is a comparison of two lengths.
# **Again on a clean permeate side**: $\Pi_\mathrm{crit}$ exists because
# $\int_0\mathrm{d}w/\sqrt{w}$ converges, which needs $y' = 0$; with $y' > 0$ the
# extinction length is infinite and this map's right-hand axis has no critical
# point on it. Under that condition, a tube of dimensionless length 1 exports all
# of its hydrogen if and only if $\Pi \ge \Pi_\mathrm{crit} = G(1,\beta)$; the
# fraction of the tube the membrane actually needs is
# $L_\mathrm{ext} = \Pi_\mathrm{crit}/\Pi$. Below $\Pi_\mathrm{crit}$ the membrane
# is the bottleneck and the outlet still carries hydrogen; far above it, the
# membrane finishes its job in a short entrance region and **nothing downstream is
# sensitive to $\alpha_H$ at all**.
#
# Itoh's membrane sits at $\Pi \approx$ 51 against $\Pi_\mathrm{crit} \approx$ 2.8.
# That ratio (not the 5.5 % stripping length, which belongs to the clean-permeate
# separator and not to his loaded reactor) is why page `H1.4` found its measured
# conversion completely insensitive to the permeance: permeation is 18× from being
# the limiting resistance, which is the condition that puts that model on its
# co-current fast-permeation ceiling.
#
# **A numerical note: why the solve is in $s$.** The discrete cell equation
#
# $$
# \frac{w_i - w_{i-1}}{h} + \Pi\left(\frac{w_i}{w_i+\beta}\right)^{1/2} = 0
# $$
#
# is continuous and strictly increasing in $w_i$, negative at $w_i = 0$ and
# positive at $w_i = w_{i-1}$, so a **unique non-negative root exists at every cell
# and every $\Pi$**. Past the front it decays *quadratically per cell*, underflows
# to exactly 0 within a handful of cells, and $w = 0$ satisfies the cell equation
# thereafter. Iterating on $w$ itself, `pymrm.newton` from the flat state stops
# without converging (`success` is False, "stopped on a bound") as soon as the
# extinction point enters the domain. The cause is the bound, not the problem:
# near $w = 0$ the sink's slope grows like $w^{-1/2}$, the linearised step of every
# tail cell overshoots past $w = 0$, and `newton`'s bound rule then holds the step
# length at its floor of 0.5, so the remaining error only halves per iteration.
# In $s = \sqrt{w}$ the sink is linear near zero, no tail step overshoots, and the
# same `newton` call converges from the flat state at every $\Pi$ (section 4a
# measures both).
#
# The sweep below solves every point in $s$ from the flat state, and compares each
# with the **march** (`Permeator.march_guess`), which solves the same bidiagonal
# cell equations one cell at a time with a scalar root finder: an independent route
# to the same discrete solution. It also records where the comparison solve in
# $w$ stalls.
#
# The map is a statement about the *law*: the reaction, the purge and the
# equilibrium of `H1.4` are all absent.

# %%
def explicit_march(p):
    """A DEFECT for the break table: the channel marched by explicit (forward) Euler.

    ``w_i = w_(i-1) - h Pi driving(w_(i-1))`` overshoots below zero at the front;
    once a value is negative the march stops and holds it.
    """
    h, w = np.diff(p.z_f), np.empty(p.shape[0])
    prev = 1.0
    for i in range(len(h)):
        if prev >= 0.0:
            prev = prev - h[i] * p.Pi * p.driving(prev)
        w[i] = prev
    return w.reshape(p.shape)


def sweep(n=None, variable="s", march=None, pis=None, with_w=True, cmp="w"):
    """Outlet hydrogen over a range of Pi, by the solve in s and by the march.

    ``variable`` is the unknown of the library solve (from the flat state; the
    march is used as its initial guess only where that fails), ``march`` the
    second route (default ``Permeator.march_guess``). ``with_w`` also records
    whether the comparison solve in w converges from the flat state.
    """
    n = N_EXP if n is None else n
    pis = np.logspace(np.log10(0.2), np.log10(60.0), 30) if pis is None else pis
    out, needed, resid, minw, dev, w_ok = [], [], [], [], [], []
    for pi in pis:
        p = Permeator(pi, BETA, n=n, n_z=800)
        wm = p.march_guess() if march is None else march(p)
        res = p.solve(variable=variable, require=False)      # flat start, inspected
        used_march = not res.success
        if used_march:
            res = p.solve(w0=wm, variable=variable)          # raises if this fails too
        NEWTON_ITERS.append(res.nit)
        r, _ = p.residual(p.w)
        out.append(p.outlet_value()); needed.append(used_march)
        resid.append(np.abs(r).max()); minw.append(min(p.w.min(), wm.min()))
        dev.append(np.max(np.abs(p.w - wm)))
        if with_w:
            w_ok.append(Permeator(pi, BETA, n=n, n_z=800).solve(variable=cmp,
                                                                require=False).success)
    out = np.array(out)
    ex = np.array([exact_half(1.0, BETA, pi)[0] for pi in pis])
    return dict(pis=pis, out=out, n_marched=int(np.sum(needed)),
                min_w=float(min(minw)), max_resid=float(max(resid)),
                maxdev=float(np.max(np.abs(out - ex))), march_dev=float(max(dev)),
                w_ok=np.array(w_ok, dtype=bool))


SW = sweep()
PI_EX = np.logspace(np.log10(0.2), np.log10(60.0), 300)
out_ex = np.array([exact_half(1.0, BETA, pi)[0] for pi in PI_EX])

W_FAIL = ~SW["w_ok"]
assert np.array_equal(W_FAIL, SW["pis"] > PI_CRIT), "the solve in w fails exactly above Pi_crit"
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.0))
ax[0].semilogx(PI_EX, 1 - out_ex, "k--", lw=1.2, label="exact")
ax[0].semilogx(SW["pis"], 1 - SW["out"], "o", ms=4.5, label="pymrm, solved in $s$ from the flat state")
ax[0].semilogx(SW["pis"][W_FAIL], 1 - SW["out"][W_FAIL], "s", ms=8, mfc="none",
               color="tab:red", label="where the solve in $w$ stalls")
ax[0].axvline(PI_CRIT, color="tab:red", ls=":", lw=1.6)
ax[0].text(PI_CRIT * 1.1, 0.22, r"$\Pi_{crit}$", color="tab:red", fontsize=10)
ax[0].set(xlabel=r"capacity ratio $\Pi = \alpha_H/u_H^0$",
          ylabel="fraction of H$_2$ recovered (exit face)", ylim=(0, 1.05),
          title="membrane-limited on the left, flow-limited on the right")
ax[0].legend(fontsize=8, loc="lower right")
ax[0].plot([PI_ITOH], [1 - exact_half(1.0, BETA, PI_ITOH)[0]], "k*", ms=13, zorder=5)
pig = np.logspace(np.log10(PI_CRIT), np.log10(200), 200)
ax[1].loglog(pig, PI_CRIT / pig, lw=2.2, label=r"$L_{ext}/L$ needed")
ax[1].loglog(pig, (1 / pig) ** 2, lw=2.0, color="tab:green",
             label=r"driving $y_H$, $n=\frac{1}{2}$")
ax[1].loglog(pig, 1 / pig, ":", color="tab:red", lw=1.8, label=r"driving $y_H$, $n=1$")
ax[1].plot([PI_ITOH], [PI_CRIT / PI_ITOH], "k*", ms=15, zorder=5)
ax[1].annotate("Itoh's tube", (PI_ITOH, PI_CRIT / PI_ITOH),
               textcoords="offset points", xytext=(-16, 13), fontsize=9)
ax[1].set(xlabel=r"capacity ratio $\Pi$", ylabel="fraction of length  /  mole fraction",
          title="spare capacity is worth its square")
ax[1].legend(fontsize=8)
fig.tight_layout(); plt.show()
print(f"finite-volume sweep over the WHOLE range, Pi = {SW['pis'][0]:.2f} to {SW['pis'][-1]:.1f}")
print(f"  solved in s from the flat state: {len(SW['pis']) - SW['n_marched']} of {len(SW['pis'])} points "
      "converged with no other help;")
print(f"  the march agrees with every one to below 1e-9 (here {SW['march_dev']:.0e}), and min w over")
print(f"  both routes is {SW['min_w']:.1e}; the s solutions satisfy the cell equations in w")
print(f"  to below 1e-10 (here {SW['max_resid']:.0e}).")
print(f"  the comparison solve in w stalls at the {int(W_FAIL.sum())} points above Pi_crit, and only there.")
print(f"  worst |pymrm - exact| outlet fraction over the swept range: {SW['maxdev']:.2e}")
print("  So the fixed grid does place the front; the solve in w could not find it.")
print()
print(f"Pi_crit = {PI_CRIT:.3f} at beta = {BETA:.3f}, for a clean permeate side.")
print(f"At Itoh's Pi = {PI_ITOH:.1f} the same tube run as a clean-permeate separator")
print(f"would need {IT['L_ext'] * 100:.1f} % of its length and a driving mole fraction of")
print(f"{IT['y_req']:.2e}. His actual reactor loads the permeate and sits on the")
print(f"fast-permeation ceiling instead, at y = {IT['y_both']:.2e} on BOTH sides. Either")
print(f"way permeation is {IT['spare']:.0f}x from limiting, so nothing measured at the")
print("outlet of such a tube constrains alpha_H. That is a property of the")
print("operating point, not of the measurement technique.")

# %% [markdown]
# #### 4a. The failure mode, measured
#
# Worth spelling out, because "the solver stopped converging" and "the problem has
# no solution" are very different statements and only one of them is true here.
# The table runs the *same* `Permeator` and the *same* `pymrm.newton` call (monotone
# damping, lower bound 0, tolerance 1e-10) from the flat state at five capacity
# ratios from just below $\Pi_\mathrm{crit}$ up to Itoh's own point, once in $w$
# and once in $s = \sqrt{w}$, and compares the solution in $s$ with the march. For
# the solve in $w$ it also reports the fraction of iterations whose step length
# was held at 0.5, the floor of `newton`'s bound rule. Two more rows: both forms at
# 12800 cells, and `pymrm.continuation` in $\Pi$ from $0.5\,\Pi_\mathrm{crit}$ to
# 60 in each form.

# %%
def failure_mode(main="s", cmp="w", n=None):
    """Section 4a's runs. ``main`` and ``cmp`` are the unknowns of the two solves."""
    n = N_EXP if n is None else n
    rows = {}
    for ratio in (0.95, 1.05, 3.0, 10.0, PI_ITOH / PI_CRIT):
        pi = ratio * PI_CRIT
        a = Permeator(pi, BETA, n=n, n_z=800)
        ra = a.solve(variable=cmp, require=False)
        lam = np.array([h["damping"] for h in ra.history])
        b = Permeator(pi, BETA, n=n, n_z=800)
        rb = b.solve(variable=main, require=False)
        wm = b.march_guess()
        res_m, _ = b.residual(wm)
        rows[ratio] = dict(cmp_ok=ra.success, cmp_nit=ra.nit,
                           lam_half=float(np.mean(np.isclose(lam, 0.5))),
                           ok=rb.success, nit=rb.nit,
                           dev=float(np.max(np.abs(b.w - wm))) if rb.success else np.nan,
                           march_resid=float(np.abs(res_m).max()), march_min=float(wm.min()))
        if rb.success:
            NEWTON_ITERS.append(rb.nit)
    fine = {}
    for var in (main, cmp):
        q = Permeator(1.05 * PI_CRIT, BETA, n=n, n_z=12800)
        rq = q.solve(variable=var, require=False)
        fine[var] = (rq.success, rq.nit)
    routes = {}
    for var in (main, cmp):
        p = Permeator(0.5 * PI_CRIT, BETA, n=n, n_z=800)
        res_fun = p.residual_s if var == "s" else p.residual

        def fun(x, pi, p=p, res_fun=res_fun):
            p.Pi = pi
            return res_fun(x)
        # the secant predictions overshoot below 0 in the extinction tail;
        # newton projects each onto its bound and warns every time
        with warnings.catch_warnings(), np.errstate(invalid="ignore"):
            warnings.filterwarnings("ignore", message="newton: .* outside the bounds",
                                    category=RuntimeWarning)
            cont = continuation(fun, p.flat_guess(),
                                np.logspace(np.log10(0.5 * PI_CRIT), np.log10(60.0), 40),
                                newton_kwargs=dict(tol=1e-10, damping="monotone",
                                                   bounds=Permeator.BOUNDS))
        routes[var] = dict(ok=cont.success, last=cont.last_param / PI_CRIT, msg=cont.message)
    return rows, fine, routes


FM, FINE, ROUTES = failure_mode()
print("the same pymrm.newton call from the flat state (n_z = 800), in w and in s = sqrt(w)")
out = [f"   {'Pi/Pi_crit':>10}   {'in w':<42}{'in s':<22}{'s vs march'}"]
for ratio, r in FM.items():
    sw = (f"converged, {r['cmp_nit']} its" if r["cmp_ok"] else
          f"not converged, lam = 0.5 in {r['lam_half'] * 100:.0f} % of its")
    ss = f"converged, {r['nit']} its" if r["ok"] else "not converged"
    out.append(f"   {ratio:>10.2f}   {sw:<42}{ss:<22}max |w_s - w_march| < 1e-9 "
               f"(here {r['dev']:.0e})")
print("\n".join(out))
print(f"   12800 cells, 1.05 Pi_crit: in s {'converged' if FINE['s'][0] else 'not converged'}"
      f" ({FINE['s'][1]} its), in w {'converged' if FINE['w'][0] else 'not converged'}.")
MARCH_WORST_RESID = max(r["march_resid"] for r in FM.values())
MARCH_MIN_W = min(r["march_min"] for r in FM.values())
LAM_HALF_MIN = min(r["lam_half"] for k, r in FM.items() if k > 1.0)
print(f"   the march itself satisfies the assembled cell equations to below 1e-11 "
      f"(here {MARCH_WORST_RESID:.0e}),")
print(f"   with min w = {MARCH_MIN_W:.1e}.")
print()
for var in ("w", "s"):
    r = ROUTES[var]
    print(f"   continuation in Pi, solved in {var}: "
          + ("reaches every value up to 60" if r["ok"] else
             f"stops at {r['last']:.4f} Pi_crit, with the message\n"
             f"      '{r['msg']}'"))
print("   In w the continuation message names a fold; the branch has none (w is single")
print("   valued in Pi, and the same call in s follows it through).")

# %%
def tail_numbers(n=None, l_ext=None):
    """The extinction tail from the march at Itoh's Pi, and the front on three grids."""
    n = N_EXP if n is None else n
    l_ext = PI_CRIT / PI_ITOH if l_ext is None else l_ext
    wm = Permeator(PI_ITOH, BETA, n=n, n_z=800).march_guess().ravel()
    i0, iz = int(np.argmax(wm < 1e-3)), int(np.argmax(wm == 0.0))
    front = {}
    for nz in (800, 3200, 12800):
        q = Permeator(PI_ITOH, BETA, n=n, n_z=nz).march_guess().ravel()
        zf = np.linspace(0.0, 1.0, nz + 1)[int(np.argmax(q == 0.0))]
        front[nz] = zf - l_ext
    return dict(wm=wm, i0=i0, iz=iz, tail=iz - i0, n_zero=int((wm == 0.0).sum()),
                front=front,
                front_order=np.log(abs(front[800] / front[3200])) / np.log(4.0))


TL = tail_numbers()
print("the extinction tail, from the march at Itoh's Pi (n_z = 800)")
print("   " + "  ".join(f"{x:.2e}" for x in TL["wm"][TL["i0"] - 1:TL["iz"]]))
print(f"   Each cell squares the one before it, so {TL['tail']} cells after w drops")
print(f"   below 1e-3 it has underflowed to exactly 0, and stays 0 for the")
print(f"   remaining {TL['n_zero']} of {TL['wm'].size} cells. A Newton step in w linearises")
print("   the sink at each of these cells, where its slope grows like w^(-1/2); the")
print("   step overshoots past 0, which is where the bound rule takes over.")
print()
print("   the front itself is placed to O(h), and refining shows it:")
print(f"   {'n_z':>7}{'first face with w = 0':>24}{'L - L_ext':>13}")
print("\n".join(f"   {nz:>7}{off + PI_CRIT / PI_ITOH:>24.5f}{off:>13.5f}"
                for nz, off in TL["front"].items()))
print(f"   against the exact L_ext = {PI_CRIT / PI_ITOH:.5f}. The offset is the width of")
print("   the underflow tail, a fixed number of CELLS, so it shrinks like h")
print(f"   (observed order {TL['front_order']:.2f} from 800 to 3200 cells).")

# %% [markdown]
# ### 5. What a departure from ½ means
#
# Sieverts' law assumes the surface steps are fast enough to keep the dissolved
# hydrogen in equilibrium with the gas at each face. Put a **finite** surface step
# in series with the wall and that assumption breaks. The model needs no physics
# beyond the mechanism already written down: with $\theta = c/C_0$ the scaled
# dissolved concentration and $\pi = p/P_0$ the scaled pressure, and writing $a$
# for the surface exchange conductance (mol/s at unit driving) and $\alpha_H$ for
# the wall conductance,
#
# $$
# q = a\left(\pi_1 - \theta_1^2\right)
#   = \alpha_H\left(\theta_1 - \theta_2\right)
#   = a\left(\theta_2^2 - \pi_2\right),
# \qquad
# \lambda \equiv \frac{\alpha_H}{a}.
# $$
#
# The three expressions are the near surface, the wall and the far surface
# carrying the same flux; the **squares** are the dissociative stoichiometry,
# which is where the ½ came from in the first place. Both limits are analytic:
# $\lambda \to 0$ (fast surface) returns $q = \alpha_H(\sqrt{\pi_1} -
# \sqrt{\pi_2})$, Eq. 1 exactly; $\lambda \to \infty$ gives $q =
# \tfrac{a}{2}(\pi_1 - \pi_2)$, linear in pressure. In between there is no power
# law at all, only a **local** apparent exponent $n_\mathrm{app} = \partial \ln
# q/\partial \ln \pi_1$ that depends on the pressure at which it is measured.
#
# Two things this construction is, and is not. It **is** a derivation from the
# mechanism, with one dimensionless group and no fitted constants. It is **not**
# in Itoh, who assumes the exponent and never tests it, and it is compared with no
# measurement here. Read the numbers below as "what a measured exponent would
# imply", never as a result about real palladium.

# %%
def series_numbers(pi_read=0.5, surface_power=2.0, r_area=None):
    """Section 5's numbers: lambda for given apparent exponents, and its scale."""
    lam_at = {tgt: root_log(lambda lam: n_apparent(pi_read, lam, surface_power=surface_power)
                            - tgt, 1e-6, 1e8)
              for tgt in (0.55, 0.60, 0.75, 0.90)}
    area = 2 * np.pi * (RI if r_area is None else r_area) * L0
    return dict(lam_at=lam_at, area=area, flux_per_area=ALPHA_H / area,
                lim_sieverts=abs(q_series(0.5, 1e-12, surface_power=surface_power)
                                 - np.sqrt(0.5)) / np.sqrt(0.5),
                lim_linear=abs(q_series(0.5, 1e12, surface_power=surface_power)
                               - 0.5 / (2 * 1e12)) / (0.5 / (2 * 1e12)))


LAMS = np.logspace(-4, 4, 160)
n_at = {pi: np.array([n_apparent(pi, lam) for lam in LAMS]) for pi in (0.5, 0.05)}
SE = series_numbers()
LAM_AT = SE["lam_at"]

print("apparent exponent produced by a surface step in series with the wall")
out = [f"  {'lambda = alpha_H/a':<24}{'n_app at p/P0=0.5':>19}{'at p/P0=0.05':>15}"]
for lam in (1e-4, 1e-2, 0.1, 1.0, 10.0, 1e3):
    out.append(f"  {lam:<24g}{n_apparent(0.5, lam):>19.4f}{n_apparent(0.05, lam):>15.4f}")
print("\n".join(out))
print("  limits (analytic): lambda -> 0 gives exactly 1/2, lambda -> inf gives 1")
print("\n".join(f"  n_app = {tgt:.2f} at p/P0 = 0.5 needs lambda = {lam:.3g}"
                for tgt, lam in LAM_AT.items()))

print("\nputting a scale on lambda for Itoh's tube, from printed numbers only:")
print(f"  membrane area 2 pi r_i l_0 = {SE['area']:.3e} m^2")
print(f"  wall conductance alpha_H   = {ALPHA_H:.2e} mol/s"
      f"  = {SE['flux_per_area']:.3e} mol/m^2/s at unit driving")
print(f"  so lambda = 1 means a surface exchange rate of {SE['flux_per_area']:.2e} mol/m^2/s,")
print("  and a 10 % shift of the exponent (n_app = 0.55) appears already at")
print(f"  lambda = {LAM_AT[0.55]:.3g}, i.e. a surface only {1 / LAM_AT[0.55]:.0f}x faster than the wall.")
print("  Itoh prints no surface kinetics, so his membrane CANNOT be placed on")
print("  this axis from this source. The chart says what a measured exponent")
print("  would imply; it says nothing about what palladium does.")

fig, ax = plt.subplots(figsize=(7.0, 4.2))
for pi, style in zip((0.5, 0.05), ("-", "--")):
    ax.semilogx(LAMS, n_at[pi], style, lw=2.2, label=f"measured at $p/P_0$ = {pi}")
ax.axhline(0.5, color="0.5", ls=":", lw=1.2)
ax.axhline(1.0, color="0.5", ls=":", lw=1.2)
ax.text(25, 0.525, "Sieverts: bulk-diffusion controlled", fontsize=8, color="0.35")
ax.text(1.5e-4, 0.955, "surface controlled", fontsize=8, color="0.35")
ax.set(xlabel=r"$\lambda = \alpha_H / a$   (wall conductance / surface conductance)",
       ylabel=r"apparent exponent $n_{app}$", ylim=(0.45, 1.05),
       title="the exponent is a resistance ratio, and it is not a constant")
ax.legend(fontsize=8, loc="center right")
fig.tight_layout(); plt.show()
print("The two curves differ, and that is the useful part: with a surface step in")
print("series the flux is NOT a power law, so a fitted exponent depends on the")
print("pressure range it was fitted over. Note also which way thickness pushes:")
print("alpha_H = 2 pi l D C_0 / ln(r_o/r_i) RISES as the wall thins, while the")
print("surface conductance a does not, so lambda = alpha_H/a rises and THINNER")
print("membranes drift toward n = 1. That is the standard reason very thin Pd")
print("films stop obeying a square root, and it drops out of this construction")
print("without having been put in.")

# %% [markdown]
# ## Validation
#
# Six checks. Exactly **one** of them compares with a number from the source, and
# that number is itself a calculation rather than a measurement: Itoh reports no
# permeation experiment. The rest are internal identities and closed-form limits.
# Every agreement number is registered as a claim below, and the defect-injection
# table at the end re-runs each computation with one thing broken on purpose,
# because an agreement is worth exactly what it would cost to break.
#
# **Deviation convention throughout: (model − reference)/reference.**

# %%
print("check 1 - the one number the source prints: alpha_H = 4.47e-5 mol/s")
print("   Derived from Eq. 1's own inputs (D, C_0, l_0, r_i, r_o) by a pymrm")
print("   radial solve, NOT by evaluating Eq. 1's formula.")
print(f"     pymrm solve {ALPHA_FV:.5e}  vs printed {ALPHA_H:.3e}"
      f"   -> {W['dev_printed'] * 100:+.3f} %")
print("   The printed value carries three significant figures, so it pins the")
print(f"   derivation only to about +-0.1 %; the deviation is {abs(W['dev_printed']) * 100:.3f} %.")
print()
print("   READ THIS CHECK EXACTLY. It tests the transcription of D, C_0, l_0 and")
print("   the 200 um wall, the structure of Eq. 1, and the radial assembly,")
print("   GIVEN the geometry. It does NOT test the geometry, because the radii")
print("   were SELECTED by requiring this same identity to close: the check chose")
print("   them, so it cannot catch them being wrong. What it can do, and does, is")
print("   discriminate the two candidate readings of the printed '17.0 mm' by")
print(f"   {abs(W['dev_od']) * 100:.1f} % against a {abs(W['dev_printed']) * 100:.2f} % residual. "
      "The independent evidence about the")
print("   geometry is Figure 1's dimension arrows (Parameters, above), which are")
print("   read from the drawing and never from alpha_H.")
print("   It also cannot test the physics of the law: nothing here was measured.")

# %%
def wall_sensitivity(r_o=None, thick_ratio=3.0):
    """Check 2: how much a wrong nu or a coarse grid moves alpha_H on this wall."""
    r_o = RO if r_o is None else r_o
    base = alpha_wall(r_o=r_o)
    a1, a0 = alpha_wall(RI, thick_ratio * RI, nu=1), alpha_wall(RI, thick_ratio * RI, nu=0)
    return dict(thin=abs(alpha_wall(r_o=r_o, nu=0) - base) / base * 100,
                thick=abs(a0 - a1) / a1 * 100,
                grid3=abs(alpha_wall(r_o=r_o, n_r=3) - base) / base)


print("check 2 - breaking the wall solve on purpose")
print("   Each row injects one defect and reports the alpha_H that comes out.")
print()
rows = [("none (baseline: nu=1, n_r=200)", ALPHA_FV),
        ("nu=0: Cartesian cell volumes", alpha_wall(nu=0)),
        ("nu=2: spherical cell volumes", alpha_wall(nu=2)),
        ("n_r = 3 (a deliberately coarse grid)", alpha_wall(n_r=3)),
        ("radii from the '17.0 mm OD' reading", alpha_wall(RI_OD, RO_OD)),
        ("wall 10x thicker (r_o = r_i + 2 mm)", alpha_wall(RI, RI + 10 * DELTA)),
        ("D read as 9.23e-9 (one digit slipped)", alpha_wall(D=10 * D_PD))]
out = [f"   {'injected defect':<44}{'alpha_H [mol/s]':>17}{'vs printed':>13}"]
out += [f"   {lab:<44}{a:>17.4e}{(a - ALPHA_H) / ALPHA_H * 100:>+12.2f} %" for lab, a in rows]
print("\n".join(out))
BREAK1 = dict(rows)
WS = wall_sensitivity()
_pc = lambda k: (BREAK1[k] - ALPHA_H) / ALPHA_H * 100
print()
print("   Read this honestly. Deviations in this table are quoted against the")
print("   PRINTED alpha_H; the percentages quoted below are against the BASELINE")
print(f"   SOLVE, which itself sits {W['dev_printed'] * 100:+.3f} % from printed. The nu=0 row is the")
print(f"   same defect either way: {(BREAK1['nu=0: Cartesian cell volumes'] - ALPHA_FV) / ALPHA_FV * 100:+.2f} % "
      f"vs the baseline, {_pc('nu=0: Cartesian cell volumes'):+.2f} % vs printed.")
print(f"   * A wrong `nu` moves alpha_H by only {WS['thin']:.2f} % on this wall,")
print(f"     because the wall is {DELTA / RI * 100:.1f} % of the tube radius and a log profile that")
print("     shallow is nearly a straight line. Check 1 is therefore WEAK against")
print("     geometry: it would catch a factor-10 slip in D instantly and a wrong")
print("     nu only marginally. On a wall with r_o/r_i = 3 the same defect moves")
print(f"     the answer by {WS['thick']:.0f} %, so the assembly is right: it is the")
print("     thin wall, not the code, that makes the identity insensitive here.")
_e5 = abs(alpha_wall(n_r=5) - ALPHA_CF) / ALPHA_CF
print(f"   * The n_r = 3 row looks invisible ({WS['grid3']:.0e} relative), but that is a")
print(f"     coincidence: the error changes sign between 3 and 5 cells ({_e5:.1e} at")
print("     n_r = 5). The trend with n_r is what check 3 measures.")
print("   * What check 1 DOES catch decisively is a wrong constant: a slipped")
print(f"     digit in D moves it by {_pc('D read as 9.23e-9 (one digit slipped)'):+.0f} %, "
      f"a 10x wall by {_pc('wall 10x thicker (r_o = r_i + 2 mm)'):+.0f} %. It does NOT catch")
print("     a wrong pair of radii, because it is what chose them; the radii row")
print("     in this table is the discrimination of two readings, not a break.")
print("     And it constrains only r_o/r_i, so the printed 200 um wall carries")
print(f"     the geometry (section 1a puts the window at "
      f"+-{(W['ri_hi'] - W['ri_lo']) / (2 * RI) * 100:.1f} % in r_i).")

# %%
def wall_ladder(nu=1, ns=(5, 10, 20, 40, 80, 160, 320)):
    """Check 3: alpha_H against the closed form on a grid ladder.

    The noise of a solve is the relative spread of its face fluxes, which the
    exact solution makes equal: round-off of the solve, not discretisation.
    """
    errs, noise = {}, {}
    for n in ns:
        _, _, _, Q_f = wall_solve(sieverts_c(1.0, C_0), 0.0, RI, RO, D_PD, L0, n_r=n, nu=nu)
        errs[n] = abs(Q_f[0] - ALPHA_CF) / ALPHA_CF
        noise[n] = np.ptp(Q_f) / abs(np.mean(Q_f))
    return errs, noise


print("check 3 - grid convergence of the wall solve against the closed form")
WL, WN = wall_ladder()
ns = list(WL)
out = [f"     n_r = {ns[0]:>4}    rel. error {WL[ns[0]]:.3e}"]
out += [f"     n_r = {n:>4}    rel. error {WL[n]:.3e}    order {np.log2(WL[m] / WL[n]):.2f}"
        for m, n in zip(ns[:-1], ns[1:])]
print("\n".join(out))
WALL_ORDER = np.log2(WL[80] / WL[160])
WALL_ERR_OVER_NOISE = WL[160] / WN[160]
print(f"   Second order ({WALL_ORDER:.2f} from 80 to 160 cells, where the error is at least")
print(f"   1e3 times the solve's own round-off (here {WALL_ERR_OVER_NOISE:.0e})), as a centred flux on a")
print("   smooth solution should be. Beyond 320 cells the error reaches that round-off")
print("   and stops falling. The two routes share no code (one is an spsolve on")
print("   assembled operators, the other is a logarithm), but the check still cannot")
print("   fail for the right reason: it measures the discretisation, not the physics.")


def graded_wall(nu=1, gain=1.0):
    """Check 4: a wall whose D rises linearly by the factor 1 + gain across it."""
    def d_prof(r):
        return D_PD * (1.0 + gain * (r - RI) / (RO - RI))
    _, _, _, Q_f = wall_solve(sieverts_c(1.0, C_0), 0.0, RI, RO, D_PD, L0, n_r=400, nu=nu,
                              D_profile=d_prof)
    q_ref = 2 * np.pi * L0 * C_0 / quad(lambda r: 1.0 / (r * d_prof(r)), RI, RO, limit=200)[0]
    return dict(Q=Q_f[0], Q_ref=q_ref, dev=(Q_f[0] - q_ref) / q_ref,
                uniform_gap=(ALPHA_CF - q_ref) / q_ref)


print()
print("check 4 - a wall the closed form cannot handle")
print("   The point of solving rather than quoting. Let D vary across the wall")
print("   (hypothetical: no such profile is printed anywhere). The reference is")
print("   then the resistance integral Q = 2 pi l0 dc / int dr/(r D(r)), by")
print("   quadrature, a genuinely independent route.")
GR = graded_wall()
print(f"     pymrm                {GR['Q']:.6e} mol/s")
print(f"     quadrature reference {GR['Q_ref']:.6e} mol/s   -> {GR['dev']:+.2e} relative")
print(f"     (a uniform-D closed form would give {ALPHA_CF:.4e}, "
      f"{GR['uniform_gap'] * 100:+.1f} % off)")
GR0 = graded_wall(nu=0)
print(f"     break test, same graded D with nu=0: {GR0['Q']:.6e} mol/s  ({GR0['dev'] * 100:+.2f} %)")
print("     - the quadrature reference does move when the assembly is broken.")

# %%
def channel_ladder(pi, ns, n=None, ref_pi=None, noise_n=None):
    """Max |w_fv - w_exact| on a grid ladder, the exact solution at ``ref_pi``.

    Each grid is solved in s from the flat state. Also returns the noise of the
    finest solve, its distance to the march of the same discrete system (two
    routes to one discrete root; ``noise_n`` changes the march's exponent, a
    break-table defect), and the march's count of cells at exactly w = 0.
    """
    n = N_EXP if n is None else n
    ref_pi = pi if ref_pi is None else ref_pi
    errs, iters, zeros = {}, {}, {}
    for nz in ns:
        p = Permeator(pi, BETA, n=n, n_z=nz)
        p.solve()
        NEWTON_ITERS.append(p.result.nit)
        errs[nz] = float(np.max(np.abs(p.w.ravel() - exact_half(p.z_c, BETA, ref_pi))))
        iters[nz] = p.result.nit
        zeros[nz] = int(np.sum(p.march_guess() == 0.0))
    q = Permeator(pi, BETA, n=n if noise_n is None else noise_n, n_z=ns[-1])
    noise = float(np.max(np.abs(q.march_guess() - p.w)))
    order = np.log2(errs[ns[-2]] / errs[ns[-1]])
    return dict(errs=errs, iters=iters, zeros=zeros, order=order, last=p,
                over_noise=errs[ns[-1]] / max(noise, 1e-16), noise=noise)


print("check 5 - the channel against its exact solution, and its grid ladder")
print("   Pi = 0.90 Pi_crit, so the hydrogen is not exhausted inside the domain.")
print("   The reference is the inverted antiderivative G, which shares no code")
print("   with the finite-volume assembly.")
print()
CL = channel_ladder(PI_DEMO, (100, 200, 400, 800, 1600, 3200))
out = [f"   {'n_z':>6}{'max |w_fv - w_exact|':>24}{'order':>9}{'Newton its':>12}"]
prev = None
for nz, e in CL["errs"].items():
    o = "" if prev is None else f"{np.log2(prev / e):>9.2f}"
    out.append(f"   {nz:>6}{e:>24.3e}{o:>9}{CL['iters'][nz]:>12}")
    prev = e
print("\n".join(out))
print(f"   measured order over the last refinement: {CL['order']:.2f}  (first-order upwind);")
print(f"   the finest error is at least 1e3 (here {CL['over_noise']:.0e}) times the gap to the march")
print("   of the same system (at least 1e3 is required for the order to mean anything).")


def extinction_invariant(n=None):
    """max |L + G(w(L))/Pi - L_ext| / L_ext over the cells, n_z = 3200."""
    p = channel(PI_DEMO, n=N_EXP if n is None else n, n_z=3200)
    L_front = p.z_c + G_half(p.w.ravel(), BETA) / PI_DEMO
    return float(np.max(np.abs(L_front - L_EXT_DEMO)) / L_EXT_DEMO)


DEV_FRONT = extinction_invariant()
print()
print("   the extinction length, extracted from the pymrm solution:")
print("     L + G(w(L))/Pi must equal L_ext at every cell (outflow outlet, cell values).")
print(f"     exact L_ext = {L_EXT_DEMO:.6f}; the worst cell deviates by {DEV_FRONT:.2e} relative.")

print()
print("   check 5a - the same ladder ABOVE Pi_crit, where the extinction point is")
print(f"   inside the domain. Pi = {PI_ITOH:.1f} (Itoh's own), solved in s from the flat state.")
CH = channel_ladder(PI_ITOH, (200, 400, 800, 1600, 3200))
assert CH["last"].w.min() >= 0.0, "the solution must be non-negative"
out = [f"   {'n_z':>6}{'max |w_fv - w_exact|':>24}{'order':>9}{'march cells at 0':>18}"]
prev = None
for nz, e in CH["errs"].items():
    o = "" if prev is None else f"{np.log2(prev / e):>9.2f}"
    out.append(f"   {nz:>6}{e:>24.3e}{o:>9}{CH['zeros'][nz]:>18}")
    prev = e
print("\n".join(out))
MARCH_ZEROS_3200 = CH["zeros"][3200]
print(f"   first order ({CH['order']:.2f}) and every value non-negative, so the fixed grid")
print("   resolves the free boundary to O(h). The march count of cells at exactly 0 encodes")
print("   the front position plus the width of the float64 underflow tail; it is")
print("   deterministic. (Counting exact zeros of a Newton solution would count round-off.)")

print()
print("   check 5b - the same channel with a loaded permeate (y' > 0), against")
print("   an independent quadrature reference:")


def loaded(yp, ref_yp=None, length=1.0):
    """Check 5b at one permeate loading; the reference at ``ref_yp``."""
    p = channel(PI_DEMO, y_perm=yp, n_z=3200, length=length)
    idx = np.arange(199, 3200, 400)
    ref = exact_quadrature(p.z_c[idx], BETA, PI_DEMO, y_perm=yp if ref_yp is None else ref_yp)
    w_stop = BETA * yp / (1 - yp)
    return dict(err=float(np.max(np.abs(p.w.ravel()[idx] - ref))), out=p.outlet_value(),
                floor=w_stop, above=(p.outlet_value() - w_stop) / w_stop)


LD = {yp: loaded(yp) for yp in (1e-4, 1e-2)}
print("\n".join(f"     y' = {yp:.0e}:  max |fv - quadrature| = {d['err']:.2e};  outlet "
                f"{d['out']:.5e}  vs the flux-free floor {d['floor']:.5e}"
                for yp, d in LD.items()))
print("     With a loaded permeate the 1/2-law no longer strips to zero: it")
print("     strips to the pressure-matched floor beta y'/(1 - y'), and reaches")
print("     it only asymptotically. The finite extinction length of section 3 is")
print("     specific to a clean permeate side.")

# %%
def channel_breaks(n_z=800):
    """Check 6: the channel at Pi_demo with one defect each, against the exact G solution."""
    defects = [("none (baseline)", {}),
               ("sink sign flipped (+ instead of -)", dict(sink_sign=+1.0)),
               ("mole fraction replaced by the flow w itself", dict(use_mole_fraction=False)),
               ("exponent n = 1 instead of 1/2", dict(n=1.0)),
               ("nu = 1 in construct_div (radial, not axial)", dict(nu=1)),
               ("outlet forced to w = 0 (Dirichlet, not outflow)", dict(outlet="dirichlet"))]
    res = {}
    for label, kw in defects:
        p = Permeator(PI_DEMO, BETA, n_z=n_z, **kw)
        try:
            r = p.solve(require=False)
        except NewtonBreakdown as exc:                # no Newton step can be formed
            res[label] = dict(state="breakdown", note=str(exc).split(":")[0])
            continue
        if r.success:
            res[label] = dict(state="converged", err=float(
                np.max(np.abs(p.w.ravel() - exact_half(p.z_c, BETA, PI_DEMO)))))
        else:
            res[label] = dict(state="not converged", note=r.message)
    base = res["none (baseline)"]["err"]
    conv = {k: v["err"] / base for k, v in res.items()
            if v["state"] == "converged" and k != "none (baseline)"}
    return res, conv, base


print("check 6 - breaking the channel on purpose")
print("   Reference: the exact G solution. The baseline is the O(h) error of the")
print("   n_z = 800 row above.")
print()
B2, B2_RATIO, E_BASE = channel_breaks()
out = [f"   {'injected defect':<50}{'max |w_fv - exact|':>20}"]
for label, v in B2.items():
    shown = f"{v['err']:>20.3e}" if v["state"] == "converged" else f"   ({v['state']}: {v['note']})"
    out.append(f"   {label:<50}{shown}")
print("\n".join(out))
N_BREAK_CONV = len(B2_RATIO)
N_BREAK_BREAKDOWN = sum(v["state"] == "breakdown" for v in B2.values())
RATIO_MIN = min(B2_RATIO.values())
print()
print(f"   {N_BREAK_CONV} of the five defects converged to a wrong ANSWER, and each moved the")
print(f"   error by at least {RATIO_MIN:.0f}x the baseline (w for the mole fraction "
      f"{B2_RATIO['mole fraction replaced by the flow w itself']:.0f}x, n = 1 "
      f"{B2_RATIO['exponent n = 1 instead of 1/2']:.0f}x,")
print(f"   nu = 1 {B2_RATIO['nu = 1 in construct_div (radial, not axial)']:.0f}x, sink sign "
      f"{B2_RATIO['sink sign flipped (+ instead of -)']:.0f}x); {N_BREAK_BREAKDOWN} stops Newton with a singular")
print("   Jacobian (NewtonBreakdown; the outlet-Dirichlet variant).")
print("   Check 5 therefore has real power against the sink, the driving")
print("   force, the exponent, the geometry and the outlet condition. It has NO")
print("   power against the value of alpha_H itself: alpha_H enters both routes")
print("   only through Pi, so any error in it cancels exactly. Only check 1")
print("   touches alpha_H, and only as a transcription.")

print()
print("   the series-resistance model of section 5 at its two analytic limits:")
print("\n".join(f"     {name:<28} exact {lim:.6e}   model {q_series(0.5, lam):.6e}   "
                f"({(q_series(0.5, lam) - lim) / lim:+.1e})"
                for lam, lim, name in ((1e-12, np.sqrt(0.5), "lambda -> 0  : sqrt(pi)"),
                                       (1e12, 0.5 / (2 * 1e12), "lambda -> inf: pi/(2 lambda)"))))
print("     The gaps are the next-order terms of the series model at lambda = 1e-12")
print("     and 1e12, not round-off. Both limits come out of the same solver, so this")
print("     is an algebraic consistency check on the series algebra, reported as the")
print("     identity it is: it cannot fail on physics.")

# %% [markdown]
# ### Claims
#
# Every number the page states is registered below with how it is checked
# (`docs/claims-design.md`): a **value** against the committed baseline, a
# **bound** at run time and in CI, an **exact** result by equality. The page is
# deterministic (no number moves between BLAS kernels), so values carry tight
# tolerances; residuals and round-off effects are bounds, and so are convergence
# orders, each with a companion bound that its finest error is at least 1e3 times
# the noise of the solve.

# %%
MIN_ITERS, MAX_ITERS = min(NEWTON_ITERS), max(NEWTON_ITERS)
print(f"{len(NEWTON_ITERS)} channel solves in s on the page, every one from the flat state: "
      f"{MIN_ITERS} to {MAX_ITERS} Newton iterations")
C = claims
# the data file
C.value("alpha_H_printed", ALPHA_H, rtol=1e-9, unit="mol/s")
C.value("D_H_in_Pd", D_PD, rtol=1e-9, unit="m2/s")
C.value("C_0", C_0, rtol=1e-9, unit="mol/m3")
C.value("wall_thickness", DELTA, rtol=1e-9, unit="m")
# check 1 and section 1a
C.value("alpha_H_dev_vs_printed_pct", W["dev_printed"] * 100, rtol=1e-3)
C.value("alpha_H_od_reading_dev_pct", W["dev_od"] * 100, rtol=1e-3)
C.value("alpha_fv_vs_closed_form_rel", abs(W["fv_vs_closed"]), rtol=1e-2,
        note="discretisation error at n_r = 200, second order")
C.bound("wall_profile_vs_log_maxrel", W["profile_vs_log"], max=1e-10, basis="chosen",
        note="cell values against the log profile at n_r = 200, near round-off")
C.bound("wall_flux_spread_rel", W["flux_spread"], max=1e-10, basis="tolerance",
        note="div(J) = 0 solved by a direct sparse LU: round-off of the solve")
C.value("ratio_locus_ro_over_ri", W["ratio_locus"], rtol=1e-8)
C.value("wall_over_radius_pct", DELTA / RI * 100, rtol=1e-6)
C.value("ri_window_lo_mm", W["ri_lo"] * 1e3, rtol=1e-5)
C.value("ri_window_hi_mm", W["ri_hi"] * 1e3, rtol=1e-5)
C.value("ri_window_halfwidth_pct", (W["ri_hi"] - W["ri_lo"]) / (2 * RI) * 100, rtol=1e-3)
C.value("wall_190um_dev_pct", W["thin_10um"] * 100, rtol=1e-3)
C.value("wall_210um_dev_pct", W["thick_10um"] * 100, rtol=1e-3)
# check 2
C.value("wall_nu0_break_thin_pct", WS["thin"], rtol=1e-3)
C.value("wall_nu0_break_thick_pct", WS["thick"], rtol=1e-3)
C.value("wall_grid_n3_rel", WS["grid3"], rtol=1e-2)
# checks 3 and 4
C.value("wall_grid_rel_err_n320", WL[320], rtol=1e-2)
C.value("wall_grid_rel_err_n5", WL[5], rtol=1e-3)
C.bound("wall_grid_order", WALL_ORDER, min=1.95, max=2.05, basis="theory",
        note="centred flux on a smooth solution, 80 to 160 cells")
C.bound("wall_grid_error_over_noise_n160", WALL_ERR_OVER_NOISE, min=1e3, basis="chosen",
        note="the order is fitted to error, not to round-off of the solve")
C.value("graded_D_fv_vs_quadrature_rel", abs(GR["dev"]), rtol=1e-2)
C.value("graded_D_uniform_gap_pct", GR["uniform_gap"] * 100, rtol=1e-3)
# section 3
C.value("Pi_itoh", PI_ITOH, rtol=1e-6)
C.value("Pi_crit", PI_CRIT, rtol=1e-6)
C.value("Pi_over_Pi_crit_itoh", IT["spare"], rtol=1e-6)
C.value("L_ext_over_tube_itoh_clean_permeate", IT["L_ext"], rtol=1e-6)
C.value("y_H_required_half_power_clean_permeate", IT["y_req"], rtol=1e-6)
C.value("y_H_both_sides_itoh_loaded_permeate", IT["y_both"], rtol=1e-6)
C.value("y_both_over_y_required", IT["both_over_req"], rtol=1e-6)
C.value("spare_capacity_50x_square", IT["square_example"], rtol=1e-6)
C.value("demo_outlet_half_n1600", W_HALF_OUT, rtol=1e-4)
C.value("demo_outlet_lin_over_half", W_LIN_OUT / W_HALF_OUT, rtol=1e-4)
# section 4 and 4a
C.exact("sweep_points_needing_march", SW["n_marched"])
C.count("sweep_points_above_pi_crit", SW["pis"], PI_CRIT, min_margin=0.05,
        note="the solve in w fails at exactly these points (asserted)")
C.exact("sweep_points_w_form_stalls", int(W_FAIL.sum()))
C.value("channel_sweep_maxdev", SW["maxdev"], rtol=1e-3)
C.bound("sweep_min_w", SW["min_w"], min=0.0, basis="chosen",
        note="the lower bound of the solve; w < 0 is outside the power law's domain")
C.bound("sweep_max_residual", SW["max_resid"], max=1e-10, basis="tolerance",
        note="newton tol 1e-10 on s; the residual is in w per unit L, h = 1/800")
C.bound("marched_worst_residual", MARCH_WORST_RESID, max=1e-11, basis="tolerance",
        note="the march's scalar roots in the assembled cell equations: round-off")
C.bound("marched_min_w", MARCH_MIN_W, min=0.0, basis="chosen",
        note="the implicit cell equation has its root in [0, w_(i-1)]")
C.exact("w_form_flat_converges_at_itoh_pi", FM[PI_ITOH / PI_CRIT]["cmp_ok"])
C.bound("w_form_step_at_half_fraction_min", LAM_HALF_MIN, min=0.4, basis="chosen",
        note="share of iterations at the bound rule's floor 0.5, above Pi_crit")
C.exact("s_form_flat_converges_everywhere",
        all(r["ok"] for r in FM.values()) and FINE["s"][0] and SW["n_marched"] == 0)
C.exact("w_form_flat_converges_n12800", FINE["w"][0])
C.bound("s_form_vs_march_maxabs", max(max(r["dev"] for r in FM.values()), SW["march_dev"]),
        max=1e-9, basis="tolerance", note="two routes to one discrete root, newton tol 1e-10")
C.exact("continuation_reaches_60", ROUTES["w"]["ok"])
C.bound("continuation_last_pi_over_crit", ROUTES["w"]["last"], min=0.9, max=1.1,
        basis="chosen", note="where continuation in w stops")
C.exact("continuation_s_reaches_60", ROUTES["s"]["ok"])
C.exact("tail_cells_to_zero", TL["tail"])
C.exact("march_cells_at_zero_n800", TL["n_zero"])
C.bound("front_offset_order", TL["front_order"], min=0.9, max=1.1, basis="theory",
        note="the underflow tail is a fixed number of cells")
# section 5
for tgt in (0.55, 0.60, 0.75, 0.90):
    C.value(f"lambda_at_napp_{tgt:.2f}", LAM_AT[tgt], rtol=1e-6)
C.value("surface_rate_lambda1_mol_m2_s", SE["flux_per_area"], rtol=1e-6, unit="mol/m2/s")
C.value("series_limit_sieverts_rel", SE["lim_sieverts"], rtol=1e-3)
C.value("series_limit_linear_rel", SE["lim_linear"], rtol=1e-3)
# check 5
C.value("channel_fv_vs_exact_maxabs_n3200", CL["errs"][3200], rtol=1e-3)
C.bound("channel_convergence_order", CL["order"], min=0.95, max=1.05, basis="theory",
        note="first-order upwind (backward Euler in L)")
C.bound("channel_finest_error_over_noise", CL["over_noise"], min=1e3, basis="chosen",
        note="the order is fitted to error, not to the solver tolerance")
C.value("extinction_length_rel_err", DEV_FRONT, rtol=1e-3)
C.value("channel_highPi_fv_vs_exact_maxabs_n3200", CH["errs"][3200], rtol=1e-3)
C.bound("channel_highPi_convergence_order", CH["order"], min=0.95, max=1.05, basis="theory",
        note="first-order upwind, approached from below")
C.bound("channel_highPi_finest_error_over_noise", CH["over_noise"], min=1e3, basis="chosen",
        note="the order is fitted to error, not to the solver tolerance")
C.exact("march_cells_at_zero_n3200", MARCH_ZEROS_3200)
C.value("channel_yperm_1e-04_fv_vs_quad", LD[1e-4]["err"], rtol=1e-3)
C.value("channel_yperm_1e-02_fv_vs_quad", LD[1e-2]["err"], rtol=1e-3)
C.bound("loaded_outlet_above_floor_min_rel", min(d["above"] for d in LD.values()), min=0.1,
        basis="chosen", note="the outlet stays above the floor beta y'/(1 - y')")
# check 6
C.value("channel_break_min_ratio", RATIO_MIN, rtol=1e-3)
C.value("channel_break_ratio_exponent", B2_RATIO["exponent n = 1 instead of 1/2"], rtol=1e-3)
C.value("channel_break_ratio_sink_sign", B2_RATIO["sink sign flipped (+ instead of -)"], rtol=1e-3)
C.value("channel_break_ratio_nu1", B2_RATIO["nu = 1 in construct_div (radial, not axial)"],
        rtol=1e-3)
C.value("channel_break_ratio_mole_fraction",
        B2_RATIO["mole fraction replaced by the flow w itself"], rtol=1e-3)
C.exact("channel_breaks_converged", N_BREAK_CONV)
C.exact("channel_breaks_breakdown", N_BREAK_BREAKDOWN)
# the solver
C.bound("newton_iters_min", MIN_ITERS, min=3, basis="chosen",
        note="every solve on the page starts from the flat state")
C.bound("newton_iters_max", MAX_ITERS, max=30, basis="chosen",
        note="the solve in s from the flat state takes about 20 iterations at most here")
print(f"{len(C.metrics)} claims registered")

# %% [markdown]
# ### The defect-injection table
#
# Each row breaks one thing on purpose, **re-runs** the computation behind its
# claims with the defect in place, and hands the result to `claims.breaks`,
# which records whether the claim's own check catches it. Every claim must be
# failed by at least one row, and the last cell of the page raises otherwise.

# %%
BREAK = []


def brk(row, metric, value):
    """One injected defect: ``value`` is what claim ``metric`` becomes with it."""
    BREAK.append((row, metric, value, claims.breaks(row, metric, value)))


# the data file read one row off
row = "data file read one row off"
S1 = stated(offset=1)
for metric, key in [("alpha_H_printed", "alpha_H"), ("D_H_in_Pd", "D_H_in_Pd"), ("C_0", "C_0"),
                    ("wall_thickness", "wall_thickness")]:
    brk(row, metric, S1[key])

# the wall
row = "D read as 9.23e-9 (one digit slipped)"
W10 = wall_numbers(D=10 * D_PD)
for metric, key, f in [("alpha_H_dev_vs_printed_pct", "dev_printed", 100),
                       ("alpha_H_od_reading_dev_pct", "dev_od", 100),
                       ("ratio_locus_ro_over_ri", "ratio_locus", 1),
                       ("ri_window_lo_mm", "ri_lo", 1e3), ("ri_window_hi_mm", "ri_hi", 1e3),
                       ("wall_190um_dev_pct", "thin_10um", 100),
                       ("wall_210um_dev_pct", "thick_10um", 100)]:
    brk(row, metric, W10[key] * f)
brk(row, "ri_window_halfwidth_pct", (W10["ri_hi"] - W10["ri_lo"]) / (2 * RI) * 100)
row = "nu = 0 in the wall: Cartesian cell volumes"
W0 = wall_numbers(nu=0)
brk(row, "alpha_fv_vs_closed_form_rel", abs(W0["fv_vs_closed"]))
brk(row, "wall_profile_vs_log_maxrel", W0["profile_vs_log"])
brk(row, "wall_flux_spread_rel", W0["flux_spread"])
WL0, WN0 = wall_ladder(nu=0, ns=(80, 160, 320))
brk(row, "wall_grid_rel_err_n320", WL0[320])
brk(row, "wall_grid_order", np.log2(WL0[80] / WL0[160]))
brk(row, "wall_grid_rel_err_n5", wall_ladder(nu=0, ns=(5,))[0][5])
brk(row, "graded_D_fv_vs_quadrature_rel", abs(GR0["dev"]))
row = "wall 10x thicker (r_o = r_i + 2 mm)"
WS10 = wall_sensitivity(r_o=RI + 10 * DELTA)
brk(row, "wall_nu0_break_thin_pct", WS10["thin"])
brk(row, "wall_over_radius_pct", 10 * DELTA / RI * 100)
brk(row, "wall_grid_n3_rel", WS10["grid3"])
brk("thick comparison wall at r_o/r_i = 2 instead of 3", "wall_nu0_break_thick_pct",
    wall_sensitivity(thick_ratio=2.0)["thick"])
WLf, WNf = wall_ladder(ns=(640, 1280))
brk("order fitted on 640 to 1280 cells, where the error meets round-off",
    "wall_grid_error_over_noise_n160", WLf[1280] / WNf[1280])
brk("graded D read as a tripling instead of a doubling", "graded_D_uniform_gap_pct",
    graded_wall(gain=2.0)["uniform_gap"] * 100)

# section 3
row = "hydrogen per cyclohexane taken as 1 instead of 3"
IT1 = itoh_numbers(nu_H2=1.0)
for metric, key in [("Pi_itoh", "Pi"), ("L_ext_over_tube_itoh_clean_permeate", "L_ext"),
                    ("y_H_required_half_power_clean_permeate", "y_req"),
                    ("y_H_both_sides_itoh_loaded_permeate", "y_both"),
                    ("y_both_over_y_required", "both_over_req"),
                    ("Pi_over_Pi_crit_itoh", "spare"), ("Pi_crit", "Pi_crit")]:
    brk(row, metric, IT1[key])
brk("linear law (n = 1) in the spare-capacity argument", "spare_capacity_50x_square",
    itoh_numbers(n=1.0)["square_example"])
row = "1-law compared at the same Pi, not at the same inlet flux"
_, _pl, _ = demo_profiles(pi_lin=PI_DEMO)
brk(row, "demo_outlet_lin_over_half", _pl.outlet_value() / W_HALF_OUT)
row = "exponent 0.55 instead of 1/2 in the solved channel"
brk(row, "demo_outlet_half_n1600",
    channel(PI_DEMO, n=0.55, n_z=1600).outlet_value())

# section 4
row = "the library solve in w and the comparison in s (forms swapped)"
SWx = sweep(variable="w", cmp="s")
brk(row, "sweep_points_needing_march", SWx["n_marched"])
brk(row, "sweep_points_w_form_stalls", int((~SWx["w_ok"]).sum()))
FMx, FINEx, ROUTESx = failure_mode(main="w", cmp="s")
brk(row, "w_form_flat_converges_at_itoh_pi", FMx[PI_ITOH / PI_CRIT]["cmp_ok"])
brk(row, "w_form_step_at_half_fraction_min",
    min(r["lam_half"] for k, r in FMx.items() if k > 1.0))
brk(row, "w_form_flat_converges_n12800", FINEx["s"][0])
brk(row, "s_form_flat_converges_everywhere",
    all(r["ok"] for r in FMx.values()) and FINEx["w"][0] and SWx["n_marched"] == 0)
brk(row, "continuation_reaches_60", ROUTESx["s"]["ok"])
brk(row, "continuation_last_pi_over_crit", ROUTESx["s"]["last"])
brk(row, "continuation_s_reaches_60", ROUTESx["w"]["ok"])
row = "Pi_crit taken from the undiluted feed (y_C0 = 1)"
brk(row, "sweep_points_above_pi_crit",
    int(np.sum(SW["pis"] > itoh_numbers(y_C0=1.0)["Pi_crit"])))
row = "exponent 0.55 instead of 1/2 in the solved channel"
SW55 = sweep(n=0.55)
brk(row, "channel_sweep_maxdev", SW55["maxdev"])
TL55 = tail_numbers(n=0.55)
brk(row, "tail_cells_to_zero", TL55["tail"])
brk(row, "march_cells_at_zero_n800", TL55["n_zero"])
brk(row, "march_cells_at_zero_n3200",
    int(np.sum(Permeator(PI_ITOH, BETA, n=0.55, n_z=3200).march_guess() == 0.0)))
CL55 = channel_ladder(PI_DEMO, (1600, 3200), n=0.55)
brk(row, "channel_fv_vs_exact_maxabs_n3200", CL55["errs"][3200])
brk(row, "extinction_length_rel_err", extinction_invariant(n=0.55))
CH55 = channel_ladder(PI_ITOH, (1600, 3200), n=0.55)
brk(row, "channel_highPi_fv_vs_exact_maxabs_n3200", CH55["errs"][3200])
row = "front compared with L_ext of the undiluted feed (y_C0 = 1)"
IT_u = itoh_numbers(y_C0=1.0)
brk(row, "front_offset_order", tail_numbers(l_ext=IT_u["Pi_crit"] / PI_ITOH)["front_order"])
row = "the march by explicit (forward) Euler instead of the implicit cell equation"
brk(row, "sweep_min_w", sweep(march=explicit_march, with_w=False)["min_w"])
brk(row, "marched_min_w", min(float(explicit_march(Permeator(r * PI_CRIT, BETA, n_z=800)).min())
                              for r in FM))
row = "the n = 0.55 march read on the n = 1/2 system"
_pm = Permeator(PI_ITOH, BETA, n_z=800)
_r55, _ = _pm.residual(Permeator(PI_ITOH, BETA, n=0.55, n_z=800).march_guess())
brk(row, "marched_worst_residual", float(np.abs(_r55).max()))
brk(row, "sweep_max_residual", float(np.abs(_r55).max()))
_pl = Permeator(0.95 * PI_CRIT, BETA, n_z=800, outlet="zero_gradient")
_pl.solve()
brk("solve in s with a zero-gradient outlet against the march", "s_form_vs_march_maxabs",
    float(np.max(np.abs(_pl.w - Permeator(0.95 * PI_CRIT, BETA, n_z=800).march_guess()))))

# section 5
row = "exponent read at p/P0 = 0.05 instead of 0.5"
SE05 = series_numbers(pi_read=0.05)
for tgt in (0.55, 0.60, 0.75, 0.90):
    brk(row, f"lambda_at_napp_{tgt:.2f}", SE05["lam_at"][tgt])
brk("membrane area from the OD-reading radius", "surface_rate_lambda1_mol_m2_s",
    series_numbers(r_area=RI_OD)["flux_per_area"])
row = "surface step without the dissociative square (power 1)"
brk(row, "series_limit_sieverts_rel",
    abs(q_series(0.5, 1e-12, surface_power=1.0) - np.sqrt(0.5)) / np.sqrt(0.5))
brk(row, "series_limit_linear_rel",
    abs(q_series(0.5, 1e12, surface_power=1.0) - 0.5 / (2 * 1e12)) / (0.5 / (2 * 1e12)))

# checks 5 and 6
row = "exact reference evaluated at a Pi 1 % high"
CLr = channel_ladder(PI_DEMO, (1600, 3200), ref_pi=1.01 * PI_DEMO)
brk(row, "channel_convergence_order", CLr["order"])
row = "ladder 1600, 2400 read as a factor 2"
brk(row, "channel_highPi_convergence_order",
    channel_ladder(PI_ITOH, (1600, 2400))["order"])
row = "second route of the noise estimate solved with exponent 0.55"
brk(row, "channel_finest_error_over_noise",
    channel_ladder(PI_DEMO, (1600, 3200), noise_n=0.55)["over_noise"])
brk(row, "channel_highPi_finest_error_over_noise",
    channel_ladder(PI_ITOH, (1600, 3200), noise_n=0.55)["over_noise"])
row = "permeate loading of the reference ten times the solved one"
brk(row, "channel_yperm_1e-04_fv_vs_quad", loaded(1e-4, ref_yp=1e-3)["err"])
brk(row, "channel_yperm_1e-02_fv_vs_quad", loaded(1e-2, ref_yp=1e-1)["err"])
brk("tube 20 times longer in the floor test", "loaded_outlet_above_floor_min_rel",
    min(loaded(yp, length=20.0)["above"] for yp in (1e-4, 1e-2)))
row = "check 6 run on 3200 cells instead of 800"
B2f, B2f_ratio, _ = channel_breaks(n_z=3200)
brk(row, "channel_break_min_ratio", min(B2f_ratio.values()))
brk(row, "channel_break_ratio_sink_sign", B2f_ratio["sink sign flipped (+ instead of -)"])
brk(row, "channel_break_ratio_nu1", B2f_ratio["nu = 1 in construct_div (radial, not axial)"])
brk(row, "channel_break_ratio_mole_fraction",
    B2f_ratio["mole fraction replaced by the flow w itself"])
brk(row, "channel_break_ratio_exponent", B2f_ratio["exponent n = 1 instead of 1/2"])
row = "zero-gradient outlet in place of w = 0 at the exit"
_rz = Permeator(PI_DEMO, BETA, n_z=800, outlet="zero_gradient").solve(require=False)
brk(row, "channel_breaks_breakdown", N_BREAK_BREAKDOWN - int(_rz.success))
brk(row, "channel_breaks_converged", N_BREAK_CONV + int(_rz.success))

# the solver
_e = Permeator(PI_ITOH, BETA, n_z=800)
brk("a start from an almost empty tube (w = 1e-6) instead of the flat state", "newton_iters_max",
    _e.solve(w0=np.full(_e.shape, 1e-6)).nit)
_m = Permeator(PI_ITOH, BETA, n_z=800)
brk("a start from the march instead of the flat state", "newton_iters_min",
    _m.solve(w0=_m.march_guess()).nit)

width = max(len(r) for r, *_ in BREAK)
out = [f"{'defect':<{width}}  {'claim':<46}{'broken':>13}  caught"]
for r, metric, value, caught in BREAK:
    v = repr(value)
    if isinstance(value, (float, np.floating)):
        v = f"(here {value:.1e})" if 0 < abs(value) < 1e-9 else f"{value:.4g}"
    out.append(f"{r:<{width}}  {metric:<46}{v:>13}  {'yes' if caught else 'no'}")
print("\n".join(out))

# %% [markdown]
# ## What pymrm adds
#
# **For the wall Itoh actually had, nothing.** A 200 µm uniform shell has a
# closed-form permeance, and the pymrm solve reproduces it to 1e-9 relative. That
# is the honest headline: section 1 is a *derivation and a transcription check*,
# not an improvement.
#
# What the solve buys is everything the closed form cannot express, and the page
# demonstrates one case: **a wall whose diffusivity varies with position**
# (check 4). Change one argument and the same assembly returns a permeance no
# $\ln(r_o/r_i)$ can give: the route to composite membranes, Pd-Ag gradients,
# loading-dependent $D$, or a supported film whose support is a second resistance
# in series. The reference there is a quadrature of the resistance integral, not
# the code being tested.
#
# For the channel, what the solve buys is a **lesson in choosing the unknown,
# with a number on it**. The finite-volume/Newton form reproduces the exact
# solution at first order and breaks under every defect the break table injects.
# Iterated in $w$, `pymrm.newton` stops converging exactly where the physics gets
# interesting, at $\Pi > \Pi_\mathrm{crit}$, where the hydrogen runs out inside the
# domain, and the tempting conclusion, that a free boundary needs a
# front-tracking method (`S12`), is **wrong here**. The discrete cell equation has
# a unique non-negative root at every cell and every $\Pi$; what fails is the
# iteration in $w$, whose tail steps all cross the bound $w = 0$. Iterated in
# $s = \sqrt{w}$, the same `newton` call converges from the flat state over the
# whole range, including Itoh's own operating point, and `continuation` in $\Pi$
# follows the branch through $\Pi_\mathrm{crit}$. A cell-by-cell march of the same
# equations agrees to below 1e-9, and the front is resolved at first order
# (section 4a, check 5a). The page keeps the stalled form visible because the
# difference between "no solution" and "an unknown in which Newton cannot reach
# the solution" is the whole reuse value.
#
# The genuinely new quantitative content is the pair of numbers this page puts on
# the law: $\Pi_\mathrm{crit} = G(1,\beta)$, the capacity ratio below which the
# membrane is the bottleneck *on a clean permeate side*, and the
# $y_H \sim \Pi^{-1/n}$ scaling that turns spare capacity into driving force.
# Neither is in the source. Together they *explain*, rather than restate, `H1.4`'s
# finding that its measured conversion cannot see $\alpha_H$: at $\Pi \approx$ 51
# against $\Pi_\mathrm{crit} \approx$ 2.8, permeation is 18× from being the
# limiting resistance, which is precisely the condition (with equilibrium
# kinetics) that puts `H1.4`'s model on its co-current fast-permeation ceiling,
# and on that ceiling the conversion depends on $K_p$ and the purge split alone,
# so the permeance has no observable left to move.

# %% [markdown]
# ## Reuse
#
# **Copy this page when your problem is a flux law that lives in a boundary
# condition.** The structure (a linear transport problem across a layer, with a
# *nonlinear equilibrium* setting the values at its faces) is not specific to
# hydrogen in palladium. The same two models cover:
#
# - oxygen through a mixed-conducting perovskite, where the surface equilibrium is
#   a different power of $p_{O_2}$;
# - solution-diffusion in a polymer film (`H1.7`), where the face concentrations
#   come from a sorption isotherm instead of a square root;
# - any supported catalyst layer where an adsorption isotherm sets the
#   concentration at the pellet surface (`B1.1`).
#
# The two pieces to change are (i) `sieverts_c`, the map from gas pressure to the
# concentration at the face, which `wall_solve` takes as its boundary values, and
# (ii) the driving force in `Permeator.driving` (the exponent `n`, the permeate
# loading `y_perm`). The operators do not change at all.
#
# **Three traps this page pays for.**
#
# - `NumJac` on a single-field 1-D problem needs shape `(n, 1)`. With `(n,)` the
#   last axis is *space*, the default stencil declares every cell coupled to every
#   other, and the Jacobian comes out dense: same answer, orders of magnitude
#   more time.
# - The permeance identity of check 1 looks like a strong check and is not, for
#   two separate reasons. On a wall this thin a wrong `nu` costs about a percent
#   (check 2 measures that instead of asserting it, and shows the identity has real
#   teeth once the wall is thick), so **measure your check's sensitivity on the
#   geometry you actually have**. And an input that was *fitted* to a check cannot
#   then be validated by it: the radii here were chosen by closing this identity,
#   so its residual says nothing about them. **When a reconstruction is
#   load-bearing, go and find a second, independent reading of the source**: here
#   Figure 1's dimension arrows, which turn a confession into a corroboration.
# - A sink with an exponent below 1 drives its variable to zero in *finite* space.
#   Iterated in that variable, `newton` from a flat start (and `continuation`)
#   stops once the zero is inside your domain, and `continuation`'s message then
#   names a fold that does not exist. **Do not read that as the discrete problem
#   having no solution**: it has one, uniquely, at every cell. The sink's slope
#   grows without bound toward zero, every tail cell's Newton step overshoots past
#   the bound, and the bound rule holds the step length at 0.5. Solve for
#   $s = \sqrt{w}$ instead: the sink becomes linear near zero, the bound
#   $s \ge 0$ keeps $w = s^2$ non-negative by construction, and the operators stay
#   the same, multiplied by $\mathrm{diag}(2s)$ in the Jacobian
#   (`Permeator.residual_s`). Reach for front tracking (`S12`) only after that.
#
# **Related pages.** [`H1.4`](../H1.4-itoh-membrane-dehydrogenation/) is this law
# inside a reactor, from the same paper; read its validation section for what a
# reactor measurement can and cannot say about a permeance.
# [`H1.7`](../H1.7-solution-diffusion/) is the polymeric counterpart, where the
# transport step is the same and the sorption step is not.
#
# ### Sources
#
# - **Read, and the source of every equation and number here:** Itoh, N. (1987),
#   "A membrane reactor using palladium", *AIChE Journal* **33**(9) 1576-1578,
#   [doi:10.1002/aic.690330921](https://doi.org/10.1002/aic.690330921): Eq. 1 and
#   the constants beneath it.
# - **Cited by Itoh for the half-power law; not consulted here:** Bohmholdt, G. &
#   Wicke, E. (1967), *Z. physik. Chem., Neue Folge* **56**, 133.
# - **Cited by Itoh for $C_0 = 1280$ mol/m³; not consulted here:** Sieverts, A. &
#   Danz, W. (1936), "Solubilities of D₂ and H₂ in Palladium", *Z. physik. Chem.
#   (B)* **34**, 158.
# - **Cited by Itoh for $D = 9.23\times10^{-10}$ m²/s; not consulted here:**
#   Nagamoto, H. & Inoue, H. (1985), *Chem. Eng. Commun.* **34**, 315.
#
# The catalog entry for this case also names Caravella et al. (2010) for the
# empirical pressure exponent. That paper is not on disk and is **not** used
# anywhere on this page; nothing here is a claim about measured exponents in real
# palladium membranes. A page that tests the exponent against permeation data
# would need such a source, and that is recorded as a follow-up on the case
# file.

# %%
# prose numbers: the Markdown cells of this notebook (README.md tags its own)
for lit, name in [
    ("18", "Pi_over_Pi_crit_itoh"), ("-2.4", "alpha_H_od_reading_dev_pct"),
    ("-0.03", "alpha_H_dev_vs_printed_pct"), ("1.023522", "ratio_locus_ro_over_ri"),
    ("2500", "spare_capacity_50x_square"), ("19", "y_both_over_y_required"),
    ("51", "Pi_itoh"), ("2.8", "Pi_crit"), ("5.5 %", "L_ext_over_tube_itoh_clean_permeate"),
    ("1e-9", "alpha_fv_vs_closed_form_rel"), ("below 1e-9", "s_form_vs_march_maxabs"),
]:
    claims.says(lit, name)

# %%
# meta.yaml (validation, adds, agreement, caveats) and the models.yaml entry
for lit, name in [
    ("4.47e-5", "alpha_H_printed"), ("-0.032", "alpha_H_dev_vs_printed_pct"),
    ("-2.4", "alpha_H_od_reading_dev_pct"), ("-0.03", "alpha_H_dev_vs_printed_pct"),
    ("1.16", "wall_nu0_break_thin_pct"), ("45", "wall_nu0_break_thick_pct"),
    ("1.023522", "ratio_locus_ro_over_ri"),
    ("8.417", "ri_window_lo_mm"), ("8.589", "ri_window_hi_mm"),
    ("2.9e-7", "graded_D_fv_vs_quadrature_rel"), ("2.4e-4", "channel_fv_vs_exact_maxabs_n3200"),
    ("1.4e-4", "extinction_length_rel_err"), ("51.4", "Pi_itoh"),
    ("4.9e-3", "channel_highPi_fv_vs_exact_maxabs_n3200"), ("3014", "march_cells_at_zero_n3200"),
    ("2e-4", "channel_yperm_1e-02_fv_vs_quad"), ("217", "channel_break_min_ratio"), ("340", "channel_break_ratio_exponent"),
    ("2898", "channel_break_ratio_sink_sign"), ("1028", "channel_break_ratio_nu1"),
    ("1.2e-6", "series_limit_sieverts_rel"), ("8.4e-7", "series_limit_linear_rel"),
    ("1e-9", "alpha_fv_vs_closed_form_rel"), ("-30.6", "graded_D_uniform_gap_pct"),
    ("2.84", "Pi_crit"), ("18", "Pi_over_Pi_crit_itoh"),
    ("5.5 %", "L_ext_over_tube_itoh_clean_permeate"),
    ("3.8e-4", "y_H_required_half_power_clean_permeate"),
    ("7.2e-3", "y_H_both_sides_itoh_loaded_permeate"), ("19", "y_both_over_y_required"),
    ("9", "tail_cells_to_zero"), ("16", "sweep_points_w_form_stalls"),
    ("below 1e-10", "sweep_max_residual"), ("7.2e-7", "wall_grid_rel_err_n5"),
    ("217", "channel_break_ratio_mole_fraction"), ("below 1e-9", "s_form_vs_march_maxabs"), ("1280", "C_0"), ("9.23e-10", "D_H_in_Pd"),
    ("2.4", "wall_over_radius_pct"), ("at least 1e3", "wall_grid_error_over_noise_n160"),
]:
    claims.says(lit, name)

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
