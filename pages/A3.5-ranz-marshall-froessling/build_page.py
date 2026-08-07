#!/usr/bin/env python3
"""Generate index.ipynb for page A3.5 (Ranz-Marshall / Froessling).

Design note: every number that appears in a markdown cell is computed here, at
build time, by executing the notebook's own code cells in order and then
substituting the results into the prose. Placeholders are written «name» so that
LaTeX braces are never touched. Nothing in the prose is typed by hand.

Run from the page directory:  python build_page.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell

# =====================================================================
# The notebook's code cells, in order. Each is a string; build_page.py
# execs them sequentially into one namespace so the prose can quote the
# results.
# =====================================================================

CODE: list[str] = []

CODE.append('''try:
    import pymrm
except ImportError:
    %pip install -q pymrm pyyaml''')

CODE.append('''import sys, urllib.request
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
from scipy.integrate import quad
from scipy.optimize import brentq, curve_fit, minimize_scalar
from scipy.sparse.linalg import spsolve
from pymrm import construct_grad, construct_div, construct_convflux_upwind
from gallery_utils import load_data, load_meta, cite_data, report_agreement

PAGE = "A3.5-ranz-marshall-froessling"
np.random.seed(0)                      # nothing here is stochastic; pinned anyway
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
M = {}                                 # every reported metric lands here
print("environment ready")''')

# ------------------------------------------------------------------ data
CODE.append('''t1 = load_data("ranz-marshall-1952-table1.csv", page=PAGE)
t2 = load_data("ranz-marshall-1952-table2.csv", page=PAGE)
t3 = load_data("ranz-marshall-1952-table3.csv", page=PAGE)
t4 = load_data("ranz-marshall-1952-table4.csv", page=PAGE)
K  = load_data("ranz-marshall-1952-printed-constants.csv", page=PAGE).set_index("quantity")

print(cite_data(load_meta("ranz-marshall-1952-table1.csv", page=PAGE)), "\\n")

def kp(name):
    """A printed constant, by name, as a float. Never retype one of these."""
    return float(K.loc[name, "value"])

# The still-air row is not a forced-convection run and never joins a regression.
still   = t1[t1.run == 19].iloc[0]
t1f     = t1[t1.run != 19].copy()
t1f["sqRe"] = np.sqrt(t1f.Re.to_numpy())
t4       = t4.copy(); t4["sqRe"] = np.sqrt(t4.Re.to_numpy())
mass_ok  = t1f.Nu_prime.notna().to_numpy()

print(f"Table 1: {len(t1f)} forced-convection runs, Re {t1f.Re.min():g} to {t1f.Re.max():g}, "
      f"plus one still-air row")
print(f"         mass-transfer column usable on {mass_ok.sum()} of {len(t1f)} "
      f"(ink damage on runs {list(t1f.run[~mass_ok])})")
print(f"Table 2: {len(t2)} runs, air {t2.air_temp_C.min():g}-{t2.air_temp_C.max():g} C")
print(f"Table 3: {len(t3)} rows, air {t3.air_temp_C.min():g}-{t3.air_temp_C.max():g} C, "
      f"abscissa printed")
print(f"Table 4: {len(t4)} benzene runs, Re {t4.Re.min():g} to {t4.Re.max():g}")
print(f"still air: N_Nu = {still.Nu_apparent:g}, N_Nu' = {still.Nu_prime:g}, "
      f"P = {still.press_mmHg:g} mm Hg")''')

# ------------------------------------------------------- printed constants
CODE.append('''A_TH  = kp("correlation_intercept")            # 2.0   theoretical
K_FIT = kp("correlation_coefficient_K")         # 0.60  fitted
N_RE  = kp("exponent_Re")                       # 1/2   assumed
M_PR  = kp("exponent_Pr_Sc")                    # 1/3   assumed
K_AQ  = kp("aqueous_collapsed_coefficient")     # 0.54  derived, eq. (24)

print(K[["value", "role"]].to_string(), "\\n")
print("roles:", dict(K.role.value_counts()))

# The three printed constants are not independent: eq. (24) collapses eqs. (21)
# for air, so it fixes the Prandtl number the authors used.
#
# THIS IS A DEFINITION, NOT A CHECK. PR_IMPLIED is (0.54/0.60)^3 and nothing
# independent is compared with it -- the page's own rule forbids importing a
# property value, so nothing can be. It is reported because every Prandtl number
# used downstream is this one, and a reader is entitled to see where it came from.
PR_IMPLIED = (K_AQ / K_FIT) ** (1.0 / M_PR)
M["pr_air_implied_by_eq24"] = PR_IMPLIED
print(f"\\neq. (24)'s {K_AQ:g} against eqs. (21)'s {K_FIT:g} with the printed 1/3 exponent")
print(f"  => the Prandtl number of air the authors used: {PR_IMPLIED:.5f}")
print(f"  (a DEFINITION, not an agreement: it is (0.54/0.60)^3 and this page never")
print(f"   compares it with a tabulated Prandtl number, because importing one is")
print(f"   exactly what the page refuses to do.)")

def nu_corr(Re, Pr_or_Sc=None, A=A_TH, Kc=K_FIT, n=N_RE, m=M_PR):
    """Eqs. (21)/(22).  Pass Pr_or_Sc=None to use the aqueous collapsed form."""
    Re = np.asarray(Re, float)
    if Pr_or_Sc is None:
        return A + K_AQ * Re ** n
    return A + Kc * np.asarray(Pr_or_Sc, float) ** m * Re ** n

def nu_from_abscissa(X, A=A_TH, Kc=K_FIT):
    """Eq. (21) evaluated on the paper's own printed abscissa N_Re^(1/2) N_Pr^(1/3)."""
    return A + Kc * np.asarray(X, float)

def pct(model, data):
    model, data = np.asarray(model, float), np.asarray(data, float)
    return 100.0 * (model - data) / data

def rms_pct(model, data):
    return float(np.sqrt(np.mean(pct(model, data) ** 2)))''')

# ----------------------------------------------- transcription / Re check
CODE.append("""# The Reynolds column is over-determined by the columns beside it: the tables
# define N_Re = D_p v rho / mu and print D_p, v, air temperature and pressure, so
# Re/(D_p v P) is a function of temperature alone. Fit
#     ln( Re / (D_p v P) ) = c0 + c1 ln T_ref
# and look at the residuals. A misread digit cannot hide there. But WHICH T_ref?
# Neither part states the convention for the TABLES. (Part II states one for the
# single still-air run -- "an average temperature of 290 deg K." -- and that
# statement is reconciled with the recovered value at the end of this cell.) So the
# reference temperature is left as an unknown weight w between the air and the drop:
#     T_ref = T_air + w (T_drop - T_air),   w = 0 free stream, w = 1/2 film, w = 1 drop.
def re_frame(df, dp_col="D_p_cm"):
    d = df[(df.Re > 0) & df.drop_temp_C.notna()]
    return pd.DataFrame({"D": d[dp_col], "v": d.air_vel_cm_s, "Ta": d.air_temp_C,
                         "Td": d.drop_temp_C, "P": d.press_mmHg, "Re": d.Re})

RE_ALL = pd.concat([re_frame(t1).assign(tab=1), re_frame(t2).assign(tab=2),
                    re_frame(t4).assign(tab=4)], ignore_index=True)
yv = np.log(RE_ALL.Re / (RE_ALL.D * RE_ALL.v * RE_ALL.P)).to_numpy()

def re_fit(w, frame=None, y=None):
    F = RE_ALL if frame is None else frame
    Y = yv if y is None else y
    T = (F.Ta + w * (F.Td - F.Ta) + 273.15).to_numpy()
    X = np.column_stack([np.ones(len(Y)), np.log(T)])
    c, *_ = np.linalg.lstsq(X, Y, rcond=None)
    r = Y - X @ c
    se = float(np.sqrt(np.sum(r ** 2) / (len(Y) - 2) * np.linalg.inv(X.T @ X)[1, 1]))
    return float(c[1]), se, r

print(f"{len(RE_ALL)} printed Reynolds numbers (every row with a measured drop")
print(f"temperature) against one rho/mu = C P T_ref^n law.\\n")
print(f"  {'w':>6}{'reference temperature':>24}{'exponent n':>13}{'+-':>8}"
      f"{'rms':>9}{'worst':>9}")
for w, lab in ((0.0, "air (free stream)"), (1/3, "one-third film"),
               (0.5, "film (mean)"), (1.0, "drop surface")):
    n_, se_, r_ = re_fit(w)
    print(f"  {w:6.3f}{lab:>24}{n_:13.3f}{se_:8.3f}"
          f"{100*np.sqrt(np.mean(r_**2)):8.3f} %{100*np.max(np.abs(r_)):8.3f} %")

# Minimised, not sampled.
opt = minimize_scalar(lambda w: float(np.sqrt(np.mean(re_fit(w)[2] ** 2))),
                      bounds=(0.0, 1.0), method="bounded",
                      options={"xatol": 1e-10})
W = float(opt.x)
n_best, se_best, r_best = re_fit(W)
M["re_reference_temperature_weight"] = W
M["re_law_T_exponent"]    = n_best
M["re_law_resid_rms_pct"] = float(100 * np.sqrt(np.mean(r_best ** 2)))
M["re_law_resid_max_pct"] = float(100 * np.max(np.abs(r_best)))
n_air, se_air, r_air = re_fit(0.0)
M["re_law_resid_rms_pct_free_stream"] = float(100 * np.sqrt(np.mean(r_air ** 2)))
M["re_law_T_exponent_free_stream"]    = n_air

lo = brentq(lambda w: np.sqrt(np.mean(re_fit(w)[2]**2))
            - 1.05*np.sqrt(np.mean(r_best**2)), 0.0, W, xtol=1e-10)
hi = brentq(lambda w: np.sqrt(np.mean(re_fit(w)[2]**2))
            - 1.05*np.sqrt(np.mean(r_best**2)), W, 1.0, xtol=1e-10)
M["re_reference_weight_band"] = float(hi - lo)

print(f"\\n  MINIMISED (bounded Brent, not swept): w = {W:.4f}, i.e. the arithmetic")
print(f"  mean of the air and drop temperatures to {abs(W-0.5)/0.5*100:.2f} %.")
print(f"  There the exponent is {n_best:+.3f} +- {se_best:.3f}, which is what kinetic")
print(f"  theory gives for air (mu ~ T^0.7 and rho ~ P/T give n ~ -1.7); at the free")
print(f"  stream it is {n_air:+.3f} +- {se_air:.3f}, {abs((n_air+1.7)/se_air):.0f} standard")
print(f"  errors away from it. The residual falls from "
      f"{M['re_law_resid_rms_pct_free_stream']:.3f} % to {M['re_law_resid_rms_pct']:.3f} %.")
print(f"  The band of w within 5 % of the best residual is [{lo:.3f}, {hi:.3f}],")
print(f"  width {M['re_reference_weight_band']:.3f} -- so the film rule is identified,")
print(f"  but a 0.45 or a 0.55 rule would not be distinguishable from it.")
print(f"\\n  THE PAPERS DO NOT STATE THIS CONVENTION FOR THE TABLES; it is recovered")
print(f"  here from the printed columns alone, and it matters downstream.")
print(f"  Re-evaluating every row's Reynolds number at the free stream instead of")
print(f"  the film shifts it by")
_ratio = ((RE_ALL.Ta + 273.15) / (RE_ALL.Ta + W*(RE_ALL.Td - RE_ALL.Ta) + 273.15)
          ).to_numpy() ** n_best
_nu_now = nu_corr(RE_ALL.Re.to_numpy())
_nu_alt = nu_corr(RE_ALL.Re.to_numpy() / _ratio)
M["re_freestream_worst_Re_pct"] = float(100 * np.max(np.abs(1/_ratio - 1)))
M["re_freestream_worst_Nu_pct"] = float(100 * np.max(np.abs(_nu_alt/_nu_now - 1)))
print(f"  up to {M['re_freestream_worst_Re_pct']:.1f} %, and the N_Nu that eq. (24)")
print(f"  then predicts by up to {M['re_freestream_worst_Nu_pct']:.1f} %. See Reuse.")
for tab in (1, 2, 4):
    sel = (RE_ALL.tab == tab).to_numpy()
    print(f"    Table {tab}: {sel.sum():2d} rows, rms "
          f"{100*np.sqrt(np.mean(r_best[sel]**2)):.3f} %, worst "
          f"{100*np.max(np.abs(r_best[sel])):.3f} %")
worst_i = int(np.argmax(np.abs(r_best)))
print(f"  worst single row: Table {RE_ALL.tab[worst_i]}, v = {RE_ALL.v[worst_i]:g} cm/s, "
      f"Re = {RE_ALL.Re[worst_i]:g}, off by {100*r_best[worst_i]:+.2f} %")
print(f"  (Table 2's runs 1-4 are excluded throughout: they print no drop temperature,")
print(f"   the table's own note saying adiabatic saturation was assumed instead.)")

# A NEARLY fit-free version of the same check. Table 1's runs 6-14 share a pressure
# and sit inside 0.7 C of each other in air temperature, so their nine Re/v values
# are very nearly one number. NOT exactly one number, and this page had that wrong
# in an earlier draft: their DROP temperatures spread the film temperature over
# about 1 K, which at the fitted exponent is a real effect of its own, computed
# below. And the quantity being bounded is a nine-row peak-to-peak, so the bound it
# has to be compared against is a peak-to-peak bound, not one row's half-width.
iso  = t1f[(t1f.press_mmHg == 738)]
r_v  = (iso.Re / iso.air_vel_cm_s).to_numpy()
M["t1_iso_Re_over_v_spread_pct"] = float(100.0 * np.ptp(r_v) / r_v.mean())

T_iso = (iso.air_temp_C + W * (iso.drop_temp_C - iso.air_temp_C) + 273.15).to_numpy()
rm_iso = T_iso ** n_best
M["t1_iso_expected_spread_from_T_pct"] = float(100.0 * np.ptp(rm_iso) / rm_iso.mean())

def _half(x):                      # rounding half-width of a printed decimal
    s = f"{x:g}"
    return 0.5 * 10.0 ** (-(len(s.split(".")[1]) if "." in s else 0))
_hw = np.array([100 * (_half(r.Re) / r.Re + _half(r.air_vel_cm_s) / r.air_vel_cm_s)
                for _, r in iso.iterrows()])
M["t1_iso_rounding_ptp_bound_pct"] = float(np.sort(_hw)[-1] + np.sort(_hw)[-2])

print(f"\\nTable 1 runs {iso.run.min()}-{iso.run.max()} (all 738 mm Hg, "
      f"{iso.air_temp_C.min()}-{iso.air_temp_C.max()} C air): Re/v = {r_v.mean():.5f}")
print(f"  observed peak-to-peak spread            {M['t1_iso_Re_over_v_spread_pct']:6.3f} %")
print(f"  film temperature moves over these rows  "
      f"{T_iso.min():.2f}-{T_iso.max():.2f} K, which at n = {n_best:.3f} is")
print(f"  a real spread in rho/mu of              "
      f"{M['t1_iso_expected_spread_from_T_pct']:6.3f} %")
print(f"  rounding alone allows, peak-to-peak,    "
      f"{M['t1_iso_rounding_ptp_bound_pct']:6.3f} %  "
      f"(worst single row +-{_hw.max():.3f} %)")
print(f"  So this is NOT a reference-temperature-free check: the temperature effect")
print(f"  alone is LARGER than the whole observed spread, and the spread sits well")
print(f"  inside temperature plus rounding. What it is, is a gross-transcription")
print(f"  detector -- a misread digit in Re or v shows up as a single outlier, and")
print(f"  the break table exercises exactly that.")

# WHAT THE PAPER ITSELF SAYS ABOUT THE REFERENCE TEMPERATURE. Part II states the
# convention once, for the still-air run: D_v = 0.204 sq.cm./sec. "at an average
# temperature of 290 deg K." Table 1's still-air row prints the air and the drop
# temperature that go with it, so the two can be reconciled -- and must be, because
# both numbers are already loaded on this page.
T_STATED = kp("Dv_backout_temperature_K")
P_STATED = kp("Dv_backout_pressure_mmHg")
T_still_film = float(still.air_temp_C + 0.5 * (still.drop_temp_C - still.air_temp_C)
                     + 273.15)
M["stillair_film_minus_stated_K"]  = float(T_still_film - T_STATED)
print(f"\\n  AND THE PAPER STATES IT ONCE ITSELF. Part II gives D_v = "
      f"{kp('Dv_backed_out_cm2_s'):g} sq.cm/sec")
print(f"  'at an average temperature of {T_STATED:g} deg K. and a pressure of "
      f"{P_STATED:g} mm. Hg.'")
print(f"  Table 1's still-air row prints air {still.air_temp_C:g} C and drop "
      f"{still.drop_temp_C:g} C at {still.press_mmHg:g} mm Hg.")
print(f"    free stream        {still.air_temp_C + 273.15:7.2f} K")
print(f"    film (w = 1/2)     {T_still_film:7.2f} K   <- the authors' 'average'")
print(f"    drop surface       {still.drop_temp_C + 273.15:7.2f} K")
print(f"    stated             {T_STATED:7.2f} K,  difference from the film "
      f"{M['stillair_film_minus_stated_K']:+.2f} K")
print(f"  and the stated pressure equals the row's own printed pressure exactly.")
print(f"  So w = {W:.4f} is not only a minimisation: for the one run where the")
print(f"  authors say which average they took, they took the film temperature.")""")

# --------------------------------------------- Table 3 as a held-out test of w
CODE.append('''# AND THE FIT HAS AN OUT-OF-SAMPLE TEST IT DID NOT ASK FOR. Table 3 prints NO
# Reynolds column at all -- only the abscissa N_Re^(1/2) N_Pr^(1/3) -- so not one of
# its rows can enter the fit above, which used the 36 rows of Tables 1, 2 and 4.
# But it prints D_p, v, the air and drop temperatures, and its footnote prints the
# pressure, so the recovered law can PREDICT its abscissa and be scored against the
# printed one. And it does so 80 K beyond where the law was fitted.
P3 = kp("table3_pressure_mmHg")

def t3_abscissa(w, frame=None):
    """Predict Table 3's printed abscissa from the rho/mu law refitted at this w."""
    d3 = t3 if frame is None else frame
    c = np.polyfit(np.log((RE_ALL.Ta + w*(RE_ALL.Td - RE_ALL.Ta) + 273.15).to_numpy()),
                   yv, 1)                                   # [slope, intercept]
    T3 = (d3.air_temp_C + w*(d3.drop_temp_C - d3.air_temp_C) + 273.15).to_numpy()
    Re3 = (np.exp(c[1]) * T3 ** c[0]
           * d3.D_p_corresponding_cm.to_numpy() * d3.air_vel_cm_s.to_numpy() * P3)
    return np.sqrt(Re3) * PR_IMPLIED ** M_PR

X3p = t3.Re_half_Pr_third.to_numpy()

def t3_rms(w, frame=None):
    return float(np.sqrt(np.mean((100*(t3_abscissa(w, frame)/X3p - 1)) ** 2)))

T_fit = (RE_ALL.Ta + W*(RE_ALL.Td - RE_ALL.Ta) + 273.15).to_numpy()
T_t3  = (t3.air_temp_C + W*(t3.drop_temp_C - t3.air_temp_C) + 273.15).to_numpy()
print(f"reference temperatures the law was FITTED on : "
      f"{T_fit.min():.1f} - {T_fit.max():.1f} K  ({len(RE_ALL)} rows)")
print(f"reference temperatures Table 3 asks it about : "
      f"{T_t3.min():.1f} - {T_t3.max():.1f} K  ({len(t3)} rows, NONE in the fit)")
print(f"  the two ranges do not overlap; the extrapolation is "
      f"{T_t3.max() - T_fit.max():.0f} K beyond the fit.\\n")

M["t3_heldout_rms_pct"]              = t3_rms(W)
M["t3_heldout_rms_pct_free_stream"]  = t3_rms(0.0)
print(f"  {'w':>6}{'reference temperature':>24}{'Table 3 abscissa rms':>22}{'worst':>10}")
for w, lab in ((0.0, "air (free stream)"), (W, "the RECOVERED value"),
               (1.0, "drop surface")):
    r_ = 100*(t3_abscissa(w)/X3p - 1)
    print(f"  {w:6.3f}{lab:>24}{np.sqrt(np.mean(r_**2)):19.3f} %"
          f"{np.max(np.abs(r_)):9.3f} %")

# Minimised on Table 3 ALONE -- rows the fit never saw.
M["t3_heldout_w"] = float(minimize_scalar(t3_rms, bounds=(0.0, 1.0),
                                          method="bounded",
                                          options={"xatol": 1e-10}).x)
print(f"\\n  and the w that best reproduces Table 3 ON ITS OWN is "
      f"{M['t3_heldout_w']:.4f},")
print(f"  against the {W:.4f} recovered from the other three tables -- rows the fit")
print(f"  never saw, at temperatures beyond every one it did see. This is the")
print(f"  held-out test the page otherwise does not have, and it is a test of the")
print(f"  reference-temperature convention only: the 0.60 was fitted to Table 3")
print(f"  (Figure 9), so nothing here is out of sample for the CORRELATION.")''')

# ------------------------------------------------------------- heat fits
CODE.append('''# HEAT TRANSFER, Table 1. The correlation's 0.60 was fitted to these rows
# (Figs. 6 and 7), so what follows is a goodness of fit -- EXCEPT for the
# intercept, which is theoretical and was never fitted to anything.
sq, Nu = t1f.sqRe.to_numpy(), t1f.Nu_corrected.to_numpy()

a_free, b_free = np.linalg.lstsq(np.column_stack([np.ones_like(sq), sq]), Nu, rcond=None)[0]
b_fixed = float(np.sum(sq * (Nu - A_TH)) / np.sum(sq ** 2))     # intercept held at 2.0

M["t1_heat_intercept_free"] = float(a_free)
M["t1_heat_slope_free"]     = float(b_free)
M["t1_heat_slope_at_A2"]    = b_fixed
M["t1_heat_K_recovered"]    = b_fixed / PR_IMPLIED ** M_PR      # back to eqs. (21)'s K

print("Table 1, 18 forced-convection runs, N_Nu (corrected) against N_Re^(1/2)")
print(f"  free two-parameter fit : N_Nu = {a_free:.4f} + {b_free:.4f} Re^(1/2)")
print(f"  intercept held at {A_TH:g}   : N_Nu = {A_TH:g} + {b_fixed:.4f} Re^(1/2)")
print(f"  eq. (24) as printed    : N_Nu = {A_TH:g} + {K_AQ:g} Re^(1/2)")
print(f"  so the coefficient these rows give, taken back through Pr^(1/3), is "
      f"K = {M['t1_heat_K_recovered']:.4f} against the printed {K_FIT:g} "
      f"({100*(M['t1_heat_K_recovered']/K_FIT - 1):+.2f} %)")

# Is the assumed 1/2 exponent supported? Fit it.
f_free = lambda R, Kk, q: A_TH + Kk * R ** q
(Kq, q_heat), _ = curve_fit(f_free, t1f.Re.to_numpy(), Nu, p0=[K_AQ, N_RE])
M["t1_heat_Re_exponent_fitted"] = float(q_heat)
print(f"\\n  Reynolds exponent released: q = {q_heat:.4f} against the assumed "
      f"{N_RE:g} ({100*(q_heat/N_RE - 1):+.2f} %)")

# Goodness of fit, and the null baselines it has to beat.
M["t1_heat_eq24_rms_pct"] = rms_pct(nu_corr(t1f.Re), Nu)
M["t1_heat_eq24_max_pct"] = float(np.max(np.abs(pct(nu_corr(t1f.Re), Nu))))
M["t1_heat_freefit_rms_pct"] = rms_pct(a_free + b_free * sq, Nu)
nulls = {"constant N_Nu (the mean)": np.full_like(Nu, Nu.mean()),
         "the intercept alone, N_Nu = 2.0": np.full_like(Nu, A_TH),
         "no intercept, N_Nu = 0.54 Re^(1/2)": K_AQ * sq}
M["t1_heat_null_best_rms_pct"] = min(rms_pct(v, Nu) for v in nulls.values())

print(f"\\n  eq. (24) against the 18 rows : rms {M['t1_heat_eq24_rms_pct']:.3f} %, "
      f"worst {M['t1_heat_eq24_max_pct']:.3f} %   [IN-SAMPLE]")
print(f"  free two-parameter fit       : rms {M['t1_heat_freefit_rms_pct']:.3f} %")
for name, v in nulls.items():
    print(f"  null: {name:36s} rms {rms_pct(v, Nu):6.2f} %")''')

# ------------------------------------------------------------- mass fits
CODE.append('''# MASS TRANSFER. Two species, and they do not agree about the intercept.
sqm, Nup = sq[mass_ok], t1f.Nu_prime.to_numpy()[mass_ok]
am_free, bm_free = np.linalg.lstsq(np.column_stack([np.ones_like(sqm), sqm]), Nup,
                                   rcond=None)[0]
bm_fixed = float(np.sum(sqm * (Nup - A_TH)) / np.sum(sqm ** 2))

sqb, Nub = t4.sqRe.to_numpy(), t4.Nu_prime.to_numpy()
ab_free, bb_free = np.linalg.lstsq(np.column_stack([np.ones_like(sqb), sqb]), Nub,
                                   rcond=None)[0]
bb_fixed = float(np.sum(sqb * (Nub - A_TH)) / np.sum(sqb ** 2))

M["t1_mass_intercept_free"]      = float(am_free)
M["t1_mass_slope_free"]          = float(bm_free)
M["t4_benzene_intercept_free"]   = float(ab_free)
M["t4_benzene_slope_free"]       = float(bb_free)
M["mass_intercept_dev_water_pct"]   = float(100 * (am_free / A_TH - 1))
M["mass_intercept_dev_benzene_pct"] = float(100 * (ab_free / A_TH - 1))
M["heat_intercept_dev_pct"]         = float(100 * (a_free / A_TH - 1))

print("N_Nu' extrapolated to N_Re = 0 by the SAME free two-parameter fit,")
print(f"against the theoretical intercept {A_TH:g}:\\n")
print(f"  heat,    water,  Table 1 ({len(sq):2d} runs) : {a_free:.4f}  "
      f"({M['heat_intercept_dev_pct']:+.2f} %)")
print(f"  mass,    water,  Table 1 ({len(sqm):2d} runs) : {am_free:.4f}  "
      f"({M['mass_intercept_dev_water_pct']:+.2f} %)")
print(f"  mass,  benzene,  Table 4 ({len(sqb):2d} runs) : {ab_free:.4f}  "
      f"({M['mass_intercept_dev_benzene_pct']:+.2f} %)")

# The Schmidt exponent is the only exponent a second species can touch.
M["sc_ratio_benzene_over_water"] = float((bb_fixed / bm_fixed) ** (1.0 / M_PR))
print(f"\\n  slopes at fixed intercept: water {bm_fixed:.4f}, benzene {bb_fixed:.4f}")
print(f"  through the assumed 1/3 that is a Schmidt-number ratio of "
      f"{M['sc_ratio_benzene_over_water']:.3f} (benzene/water in air)")

# THE ONE DEGRADED GLYPH IN TABLE 4, SETTLED BY THE TABLE'S OWN ARITHMETIC.
# Both printed transfer numbers of Table 4 come from the SAME measured evaporation
# rate -- N_Nu from the heat balance, N_Nu' from eq. (16) -- so their ratio is
# rate-independent and depends only on the printed air temperature, drop temperature
# and pressure. Runs 10, 11 and 12 are printed within 0.2 K and at one pressure, so
# their three ratios have to be the same number. That is a transcription check no
# property value can touch, and it is the check that decides run 11's last digit.
trio  = t4[t4.run.isin([10, 11, 12])]
ratio = (trio.Nu_prime / trio.Nu_apparent).to_numpy()
M["t4_ratio_spread_runs_10_12_pct"] = float(100 * np.ptp(ratio) / ratio.mean())
print(f"\\nTable 4 runs 10-12, air {trio.air_temp_C.min():g}-{trio.air_temp_C.max():g} C,"
      f" drop {trio.drop_temp_C.min():g}-{trio.drop_temp_C.max():g} C, "
      f"{trio.press_mmHg.iloc[0]:g} mm Hg:")
for _, r in trio.iterrows():
    print(f"  run {int(r.run):2d}: N_Nu = {r.Nu_apparent:4.1f}, N_Nu' = "
          f"{r.Nu_prime:5.1f}, ratio {r.Nu_prime/r.Nu_apparent:.4f}")
print(f"  spread {M['t4_ratio_spread_runs_10_12_pct']:.2f} % across the three")

# The bracket runs 10 and 12 put on run 11, from the printed rounding alone.
def _rbracket(row):
    return ((row.Nu_prime - 0.05) / (row.Nu_apparent + 0.05),
            (row.Nu_prime + 0.05) / (row.Nu_apparent - 0.05))
b10, b12 = _rbracket(t4[t4.run == 10].iloc[0]), _rbracket(t4[t4.run == 12].iloc[0])
R_LO, R_HI = max(b10[0], b12[0]), min(b10[1], b12[1])
nu11 = float(t4.loc[t4.run == 11, "Nu_apparent"].iloc[0])
M["t4_run11_Nuprime_upper_bound"] = float((nu11 + 0.05) * R_HI)
print(f"  their rounding brackets intersect in [{R_LO:.4f}, {R_HI:.4f}], which puts")
print(f"  run 11's N_Nu' in [{(nu11-0.05)*R_LO:.2f}, "
      f"{M['t4_run11_Nuprime_upper_bound']:.2f}] against a printed N_Nu of {nu11:g}.")
print(f"  The first two glyphs read '10.' unambiguously and the last is 17 px wide,")
print(f"  the width of this column's round digits (16-18) and not of its 1s (10-12),")
print(f"  so the value is 10.0 -- and NOT the 10.6 an earlier draft of this page read")
print(f"  off the glyph's counter position, which would need a ratio of "
      f"{10.6/nu11:.4f},")
print(f"  {100*((10.6/nu11)/R_HI - 1):.1f} % outside a bracket built from a 1 % "
      f"rounding budget.")

# What the discarded reading would have done, printed because the drift is SMALL.
Nub_alt = Nub.copy(); Nub_alt[t4.run.to_numpy() == 11] = 10.6
ab_alt, bb_alt = np.linalg.lstsq(np.column_stack([np.ones_like(sqb), sqb]), Nub_alt,
                                 rcond=None)[0]
print(f"\\n  Had 10.6 been kept: intercept {ab_free:.4f} -> {ab_alt:.4f}, slope "
      f"{bb_free:.4f} -> {bb_alt:.4f},")
print(f"  a drift of {100*abs(ab_alt/ab_free - 1):.1f} % -- UNDER "
      f"check_agreement.py's 5 % tolerance, so no regression test would ever have")
print(f"  caught it. The ratio spread above would have gone "
      f"{M['t4_ratio_spread_runs_10_12_pct']:.2f} % -> "
      f"{100*np.ptp([ratio[0], 10.6/nu11, ratio[2]])/np.mean([ratio[0], 10.6/nu11, ratio[2]]):.2f} %.")
print(f"  That is why the spread is reported as a metric of its own and carries its")
print(f"  own break row: it is the only number on this page that moves enough to")
print(f"  fail CI when that digit is read wrong.")''')

# ---------------------------------------------- still air, analogy, noise
CODE.append('''# The paper's OWN measurement at N_Re = 0, from Fig. 8's D_p^2-vs-time slope
# through eqs. (19) and (20). These are the numbers the correlation's 2.0 is
# supposed to be.
nu0, nup0 = float(still.Nu_apparent), float(still.Nu_prime)
M["stillair_Nu_over_theory"]      = nu0 / A_TH
M["stillair_Nuprime_over_theory"] = nup0 / A_TH
M["stillair_analogy_ratio"]       = nu0 / nup0

print(f"still air, measured : N_Nu = {nu0:g}  ({100*(nu0/A_TH-1):+.1f} % on {A_TH:g})")
print(f"                      N_Nu' = {nup0:g}  ({100*(nup0/A_TH-1):+.1f} % on {A_TH:g})")
print(f"\\neqs. (21) and (22) assert N_Nu = N_Nu' = {A_TH:g} at N_Re = 0, i.e. a ratio of"
      f" exactly 1.")
print(f"The paper's own still-air pair gives {M['stillair_analogy_ratio']:.4f}.")

# What the diffusivity would have to be. N_Nu' is inversely proportional to D_v,
# so this is a statement about D_v and the paper makes it itself.
DV_BACK = kp("Dv_backed_out_cm2_s")
M["dv_used_over_dv_backed_out"] = A_TH / nup0
print(f"\\nN_Nu' ~ 1/D_v, so forcing N_Nu' = {A_TH:g} on the still-air run needs a D_v "
      f"smaller by {A_TH/nup0:.4f}.")
print(f"That is the {DV_BACK:g} sq.cm/sec the paper prints -- and calls "
      f"'a low value compared with other methods of determination'.")
print(f"Implied D_v behind Table 1's mass column: {DV_BACK*A_TH/nup0:.5f} sq.cm/sec "
      f"(an inference from two printed numbers, not a printed value).")

# How big is 10 %? The only repeatability either part offers.
reps = []
for df, pairs, cols in ((t1, [(8, 9)], ["evap_rate_ml_s_x1e5", "Nu_apparent",
                                        "Nu_corrected", "Nu_prime"]),
                        (t2, [(3, 4), (6, 7)], ["evap_rate_ml_s_x1e5", "Nu_apparent",
                                                "Nu_corrected", "Nu_prime"])):
    for i, j in pairs:
        ri, rj = df[df.run == i].iloc[0], df[df.run == j].iloc[0]
        for c in cols:
            if pd.notna(ri[c]) and pd.notna(rj[c]) and ri[c] != 0:
                reps.append((c, 100 * abs(ri[c] - rj[c]) / (0.5 * (ri[c] + rj[c]))))
M["replicate_worst_pct"] = float(max(v for _, v in reps))
M["replicate_mean_pct"]  = float(np.mean([v for _, v in reps]))
print(f"\\nReplicate pairs (Table 1 runs 8/9; Table 2 runs 3/4 and 6/7), "
      f"{len(reps)} paired quantities:")
print(f"  mean disagreement {M['replicate_mean_pct']:.2f} %, worst "
      f"{M['replicate_worst_pct']:.2f} %")
print(f"  so the mass-side intercept deficit is about "
      f"{abs(M['mass_intercept_dev_water_pct'])/M['replicate_mean_pct']:.0f} "
      f"times the run-to-run scatter, and the heat-side excess about "
      f"{abs(M['heat_intercept_dev_pct'])/M['replicate_mean_pct']:.0f} times.")''')

# ------------------------------------------------------------- Table 3
CODE.append('''# TABLE 3: the only high-temperature transfer data in either part, and the only
# table that prints the correlation's abscissa, so it can be evaluated with no
# property data whatsoever.
X3, Nu3 = t3.Re_half_Pr_third.to_numpy(), t3.Nu.to_numpy()

M["t3_eq21_rms_pct"]  = rms_pct(nu_from_abscissa(X3), Nu3)
M["t3_eq21_max_pct"]  = float(np.max(np.abs(pct(nu_from_abscissa(X3), Nu3))))
M["t3_eq21_bias_pct"] = float(np.mean(pct(nu_from_abscissa(X3), Nu3)))

# and with the coefficient refitted on Table 1's room-temperature runs alone
M["t3_from_t1_rms_pct"] = rms_pct(nu_from_abscissa(X3, Kc=M["t1_heat_K_recovered"]), Nu3)

print(f"eq. (21) on Table 3's own printed abscissa ({len(X3)} rows, air "
      f"{t3.air_temp_C.min():g}-{t3.air_temp_C.max():g} C):")
print(f"  rms {M['t3_eq21_rms_pct']:.2f} %, worst {M['t3_eq21_max_pct']:.2f} %, "
      f"mean bias {M['t3_eq21_bias_pct']:+.2f} %")
print(f"  same rows with K refitted on Table 1 alone (K = "
      f"{M['t1_heat_K_recovered']:.4f}): rms {M['t3_from_t1_rms_pct']:.2f} %")
print(f"  against Table 1's own in-sample rms of {M['t1_heat_eq24_rms_pct']:.3f} % "
      f"-- a factor {M['t3_eq21_rms_pct']/M['t1_heat_eq24_rms_pct']:.1f}")

# Two identities inside Table 3 that need no properties at all. Its nine rows are
# FOUR motion-picture runs read at two drop diameters each (pairs A-D) plus one
# single row (E), which the pairing loop below skips; within a pair everything but
# D_p is held fixed.
print("\\nWithin-pair identities (same run, two diameters):")
print(f"  {'pair':>5}  {'X ratio':>8} {'sqrt(Dp) ratio':>15} {'dev':>7}   "
      f"{'Nu ratio':>8} {'-dDp2 ratio':>12} {'dev':>7}   {'Nu ratio predicted':>19}")
dev_x, dev_n, pair_rows = [], [], []
for g, d in t3.groupby("motion_picture_pair"):
    if len(d) != 2:
        continue
    hi, lo = d.iloc[0], d.iloc[1]
    rx  = hi.Re_half_Pr_third / lo.Re_half_Pr_third
    rdp = np.sqrt(hi.D_p_corresponding_cm / lo.D_p_corresponding_cm)
    rn  = hi.Nu / lo.Nu
    rd  = hi.neg_dDp2_dtau_cm2_s_x1e4 / lo.neg_dDp2_dtau_cm2_s_x1e4
    rpred = nu_from_abscissa(hi.Re_half_Pr_third) / nu_from_abscissa(lo.Re_half_Pr_third)
    dev_x.append(100*(rx/rdp - 1)); dev_n.append(100*(rn/rd - 1))
    pair_rows.append((g, rn, float(rpred), 100.0*(rn/float(rpred) - 1)))
    print(f"  {g:>5}  {rx:8.4f} {rdp:15.4f} {100*(rx/rdp-1):+6.2f} %   "
          f"{rn:8.4f} {rd:12.4f} {100*(rn/rd-1):+6.2f} %   {float(rpred):19.4f}")

M["t3_pair_abscissa_max_dev_pct"] = float(np.max(np.abs(dev_x)))
M["t3_pair_nu_ratio_max_dev_pct"] = float(np.max(np.abs(dev_n)))
# The N_Nu ratio the correlation predicts against the one the table prints, stated
# as a per-cent difference and NOT as a ratio of logarithms: pair B's abscissa moves
# only from 2.0 to 2.1, so a log ratio would divide by a two-figure rounding.
M["t3_pair_worst_nu_ratio_error_pct"] = float(max(abs(d) for *_, d in pair_rows))
M["t3_pairs_measurement_steeper"] = float(sum(1 for *_, d in pair_rows if d > 0))
print(f"\\n  the two identities hold to {M['t3_pair_abscissa_max_dev_pct']:.2f} % and "
      f"{M['t3_pair_nu_ratio_max_dev_pct']:.2f} % (2-figure printing allows about 2.5 %),")
print(f"  so the table is self-consistent; what disagrees is the correlation. The")
print(f"  N_Nu ratio inside a run, measured against what eq. (21) predicts:")
for g, rn, rp, d in pair_rows:
    print(f"    pair {g}: measured {rn:.4f}, eq. (21) {rp:.4f}, {d:+.1f} %")
print(f"  worst {M['t3_pair_worst_nu_ratio_error_pct']:.1f} %, and the measurement is the")
print(f"  steeper of the two in {int(M['t3_pairs_measurement_steeper'])} of "
      f"{len(pair_rows)} pairs. Two-figure printing allows about 1.3 % on a")
print(f"  predicted ratio and about 3.4 % on a measured one.")''')

# ------------------------------------------------------------------ figure
CODE.append('''# Each species goes on the abscissa its OWN free fit implies, which is the only way
# to put two Schmidt numbers on one axis without importing a diffusivity that the
# papers print only as a figure. The intercepts are unaffected by that choice.
Xh  = PR_IMPLIED ** M_PR * sq                 # heat, water: the paper's own Pr
Xw  = (bm_free / K_FIT) * sqm                 # mass, water
Xb  = (bb_free / K_FIT) * sqb                 # mass, benzene

fig, ax = plt.subplots(1, 3, figsize=(14.6, 4.5))
xs = np.linspace(0, 18, 200)

ax[0].plot(xs, nu_from_abscissa(xs), "k-", lw=2, label=r"eq. (21): $2.0+0.60\,X$")
ax[0].axhline(A_TH, color="0.6", ls=":", lw=1.5, label=r"theoretical limit $2.0$")
ax[0].plot(Xh, Nu, "o", ms=5, label="Table 1, water, 20-25 C air")
ax[0].plot(X3, Nu3, "s", ms=6, mfc="none", color="tab:red",
           label="Table 3, water, 85-221 C air")
ax[0].plot(0, nu0, "*", ms=15, color="tab:orange", label="still air, measured")
ax[0].set(xlabel=r"$N_{Re}^{1/2} N_{Pr}^{1/3}$", ylabel=r"$N_{Nu}$",
          title="heat transfer", xlim=(0, 16), ylim=(0, 12))
ax[0].legend(fontsize=8, loc="upper left")

ax[1].plot(xs, nu_from_abscissa(xs), "k-", lw=2, label=r"eq. (22): $2.0+0.60\,X$")
ax[1].axhline(A_TH, color="0.6", ls=":", lw=1.5)
ax[1].plot(Xw, Nup, "o", ms=5, color="tab:blue", label="Table 1, water")
ax[1].plot(Xb, Nub, "^", ms=6, mfc="none", color="tab:green", label="Table 4, benzene")
ax[1].plot(xs, am_free + K_FIT * xs, "--", lw=1.4, color="tab:blue",
           label="free fit, water")
ax[1].plot(xs, ab_free + K_FIT * xs, "--", lw=1.4, color="tab:green",
           label="free fit, benzene")
ax[1].plot(0, nup0, "*", ms=15, color="tab:orange", label="still air, measured")
ax[1].set(xlabel=r"$N_{Sc}^{1/3} N_{Re}^{1/2}$ (each species on its own fitted scale)",
          ylabel=r"$N_{Nu}'$", title="mass transfer", xlim=(0, 18), ylim=(0, 14))
ax[1].legend(fontsize=8, loc="upper left")

ax[2].plot(xs, nu_from_abscissa(xs), "k-", lw=2)
ax[2].axhline(A_TH, color="0.6", ls=":", lw=1.5)
ax[2].plot(Xh, Nu, "o", ms=5, color="tab:purple", alpha=.6)
ax[2].plot(Xw, Nup, "o", ms=5, color="tab:blue")
ax[2].plot(Xb, Nub, "^", ms=6, mfc="none", color="tab:green")
ax[2].plot(xs, a_free  + K_FIT * xs, "--", lw=1.4, color="tab:purple")
ax[2].plot(xs, am_free + K_FIT * xs, "--", lw=1.4, color="tab:blue")
ax[2].plot(xs, ab_free + K_FIT * xs, "--", lw=1.4, color="tab:green")
ax[2].plot(0, nu0,  "*", ms=15, color="tab:orange")
ax[2].plot(0, nup0, "*", ms=15, color="tab:orange")
for v, c, lab in ((a_free, "tab:purple", "heat, water"),
                  (am_free, "tab:blue", "mass, water"),
                  (ab_free, "tab:green", "mass, benzene")):
    ax[2].plot(0, v, "x", ms=10, mew=2, color=c, label=f"{lab}: {v:.3f}")
ax[2].set(xlabel=r"$X$", ylabel=r"$N_{Nu}$ or $N_{Nu}'$",
          title="the intercept, magnified", xlim=(-0.4, 4), ylim=(1.4, 5))
ax[2].legend(fontsize=8, loc="lower right", title="extrapolated to $X=0$",
             title_fontsize=8)
fig.tight_layout(); plt.show()

print("Right-hand panel: the three free two-parameter fits extrapolated to X = 0.")
print("The stars are the paper's own still-air measurements, which are separate")
print("experiments and are NOT in any of the three fits.")''')

# ------------------------------------------------- pymrm: the theoretical 2
CODE.append('''# WHERE THE 2.0 COMES FROM. Part I obtains it from eqs. (1), (2), (5) and (6)
# with the fluid at rest: steady conduction (or diffusion) out of a sphere into an
# unbounded stagnant medium. That is a one-dimensional spherical BVP, so pymrm can
# be asked the same question directly.
#
# construct_div(..., nu=2) supplies the spherical geometry; the grid is geometric
# because the answer is set by the shell nearest the drop and the far field has to
# reach many drop radii.
def nu_shell(n, ratio, D=1.0, Ca=1.0, Cb=0.0, nu=2, geometric=True):
    a, b = 1.0, float(ratio)
    r_f = (a * (b / a) ** np.linspace(0.0, 1.0, n + 1) if geometric
           else np.linspace(a, b, n + 1))
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n, 1)
    #  a dC/dn + b C = d, outward normal.  Both ends Dirichlet, so a = 0.
    bc = ({"a": 0.0, "b": 1.0, "d": Ca},        # r = a : drop surface, saturated
          {"a": 0.0, "b": 1.0, "d": Cb})        # r = b : far field
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc)
    div = construct_div(shape, r_f, nu=nu)
    lhs = div @ (-D * grad)
    rhs = -(div @ (-D * grad_bc)).toarray().ravel()
    C = spsolve(lhs.tocsc(), rhs)
    g = (grad @ C.reshape(-1, 1) + grad_bc).ravel()      # 2nd-order face gradient
    return float(2.0 * a * (-D * g[0]) / (D * (Ca - Cb)))

RATIO = 1.0e4
print("grid convergence at b/a =", f"{RATIO:.0e}")
ns, vals = [50, 100, 200, 400, 800], []
for n in ns:
    vals.append(nu_shell(n, RATIO)); print(f"  n = {n:4d}   N_Nu = {vals[-1]:.10f}")
errs = np.abs(np.array(vals) - 2.0 / (1.0 - 1.0 / RATIO))
orders = np.log(errs[:-1] / errs[1:]) / np.log(2.0)
M["pymrm_conduction_grid_order"] = float(orders[-1])
print(f"  observed orders {np.array2string(orders, precision=2)} -> "
      f"{M['pymrm_conduction_grid_order']:.2f}")

# Richardson in h at each domain size, then extrapolate the domain to infinity.
def nu_h0(ratio, n=400):
    c, f = nu_shell(n, ratio), nu_shell(2 * n, ratio)
    return f + (f - c) / 3.0                                # 2nd order

ratios = np.array([1e3, 3e3, 1e4, 3e4, 1e5])
nuv    = np.array([nu_h0(r) for r in ratios])
p      = np.polyfit(1.0 / ratios, nuv, 1)                   # N_Nu = Nu_inf + c a/b
M["pymrm_conduction_Nu_infinity"] = float(p[1])
M["pymrm_conduction_dev_from_2"]  = float(abs(p[1] - A_TH))
print(f"\\n  b/a -> infinity (linear in a/b, Richardson in h first): "
      f"N_Nu = {p[1]:.9f}")
print(f"  the printed theoretical value is {A_TH:g}; difference "
      f"{M['pymrm_conduction_dev_from_2']:.2e}")

# SECOND, INDEPENDENT ROUTE. The same limit is a series of shell resistances,
# integral(dr / 4 pi r^2 k) from a to infinity, which adaptive quadrature can do
# on the semi-infinite interval directly. It shares no assembly, no grid and no
# linear solve with the finite-volume calculation above.
Rtot, qerr = quad(lambda r: 1.0 / (4.0 * np.pi * r ** 2), 1.0, np.inf,
                  epsabs=1e-14, epsrel=1e-14)
nu_quad = 2.0 / (4.0 * np.pi * 1.0 * Rtot)
M["nu_quadrature_route"]        = float(nu_quad)
M["nu_two_routes_rel_diff"]     = float(abs(nu_quad - p[1]) / A_TH)
print(f"\\n  quadrature of the resistance integral : N_Nu = {nu_quad:.12f} "
      f"(quad error estimate {qerr:.1e})")
print(f"  pymrm finite volume, extrapolated     : N_Nu = {p[1]:.12f}")
print(f"  the two routes differ by {M['nu_two_routes_rel_diff']:.2e} relative")''')

# ------------------------------------------------ pymrm: radial convection
CODE.append('''# RANZ'S OWN CAVEAT. In 1993 he wrote that a serious limitation of the papers,
# "unfortunately not emphasized in the original papers", was the crude accounting
# for "the radial convection at a spherical boundary caused by diffusion at the
# same boundary". The object that does that accounting is p_f, defined in Part
# II's notation only as the "average value of (pi - p_A) across transfer path"
# -- and which average is not said.
#
# Solve the actual problem: species A diffusing out of the sphere through stagnant
# B, so the mixture itself moves radially. Diffusion by construct_grad, the radial
# convection by construct_convflux_upwind, geometry by construct_div(nu=2), and
# the total molar flow W closed by a scalar root-find.
def sh_stefan(n, ratio, ya, yb=0.0, D=1.0, c=1.0):
    a, b = 1.0, float(ratio)
    r_f = a * (b / a) ** np.linspace(0.0, 1.0, n + 1)
    r_c = 0.5 * (r_f[:-1] + r_f[1:])
    shape = (n, 1)
    bc = ({"a": 0.0, "b": 1.0, "d": ya}, {"a": 0.0, "b": 1.0, "d": yb})
    grad, grad_bc = construct_grad(shape, r_f, r_c, bc)
    div = construct_div(shape, r_f, nu=2)

    def flux_at_surface(W):
        v_f = (W / (4.0 * np.pi * c)) / r_f ** 2            # molar-average velocity
        conv, conv_bc = construct_convflux_upwind(shape, r_f, r_c, bc,
                                                  v=v_f.reshape(-1, 1))
        lhs = div @ (conv - c * D * grad)
        rhs = -(div @ (conv_bc - c * D * grad_bc)).toarray().ravel()
        yv  = spsolve(lhs.tocsc(), rhs)
        g   = (grad @ yv.reshape(-1, 1) + grad_bc).ravel()
        return 4.0 * np.pi * a ** 2 * (v_f[0] * c * ya - c * D * g[0])

    W = brentq(lambda W: flux_at_surface(W) - W, 1e-12, 100.0, xtol=1e-14, rtol=8.9e-16)
    return float(W / (2.0 * np.pi * a * c * D * (ya - yb))), float(W)

def sh_closed(ratio, ya, yb=0.0):
    """Exact for the spherical stagnant film."""
    return 2.0 * np.log((1 - yb) / (1 - ya)) / ((1 - 1 / ratio) * (ya - yb))

RS = 1.0e5
print(f"{'y_surface':>10}{'pymrm Sh':>14}{'closed form':>14}{'rel':>10}"
      f"{'Sh/2 (the correction)':>24}")
ys, sh_num, sh_ex = [0.01, 0.05, 0.10, 0.20], [], []
for ya in ys:
    s, _ = sh_stefan(600, RS, ya); e = sh_closed(RS, ya)
    sh_num.append(s); sh_ex.append(e)
    print(f"{ya:10.2f}{s:14.6f}{e:14.6f}{abs(s/e-1):10.1e}{s/2.0:24.4f}")
M["stefan_pymrm_vs_closed_rel"] = float(np.max(np.abs(np.array(sh_num) /
                                                     np.array(sh_ex) - 1)))
print(f"\\n  worst pymrm-vs-closed-form disagreement {M['stefan_pymrm_vs_closed_rel']:.1e} "
      f"(the grid study below gives the observed order)")

gr = [sh_stefan(n, RS, 0.20)[0] for n in (300, 600, 1200)]
ge = sh_closed(RS, 0.20)
M["stefan_grid_order"] = float(np.log(abs(gr[1]-ge) / abs(gr[2]-ge)) / np.log(2.0))
print(f"  n = 300/600/1200 at y = 0.20: relative errors "
      f"{[f'{abs(g/ge-1):.2e}' for g in gr]}, observed order over the finest pair "
      f"{M['stefan_grid_order']:.2f}")
print(f"  (upwind convection is first order and this approaches it from below; the")
print(f"   closed form is the authority and the solve is the check on it.)")

# Now the question the notation leaves open. Ranz's N_Nu' carries a factor p_f/pi.
# If p_f is the LOGARITHMIC mean of (pi - p_A) it removes the radial-convection
# enhancement exactly and returns the theoretical 2.0; if it is the arithmetic
# mean it does not.
def pf_log(y):   return y / np.log(1.0 / (1.0 - y))        # /pi
def pf_arith(y): return 1.0 - 0.5 * y                       # /pi

yy = np.array(ys)
back_log   = np.array(sh_num) * pf_log(yy) * (1 - 1/RS)   # pymrm Sh, not the closed form
back_arith = np.array(sh_ex) * pf_arith(yy) * (1 - 1/RS)
pf_log_dev = float(np.max(np.abs(back_log / A_TH - 1)))
M["pf_arith_dev_at_y020_pct"] = float(100 * abs(back_arith[-1] / A_TH - 1))
print(f"\\n  pymrm Sh x p_f/pi, p_f the LOG mean -> {np.array2string(back_log, precision=8)}")
print(f"  (put the CLOSED FORM in place of the pymrm Sh and that product is exactly")
print(f"   {A_TH:g} for every y -- an algebraic identity, printed here and deliberately")
print(f"   NOT reported as a metric because it cannot fail. The number above is the")
print(f"   pymrm solve being asked the same question, and it can.)")
print(f"  Sh x p_f/pi with p_f the ARITHMETIC mean -> "
      f"{np.array2string(back_arith, precision=6)}")
print(f"  so the log mean is the one that makes Ranz's N_Nu' equal {A_TH:g} at rest "
      f"(to {pf_log_dev:.1e}),")
print(f"  and that {pf_log_dev:.1e} is NOT a second, independent number: it is")
print(f"  identically the pymrm-vs-closed-form disagreement above, because the log")
print(f"  mean cancels the closed form exactly. It is reported once, under that name.")
print(f"  and the arithmetic mean is already {M['pf_arith_dev_at_y020_pct']:.2f} % out "
      f"at y = {ys[-1]:g}.")

# How large is the effect at all? Root-found, not swept.
def corr_excess(y): return -np.log(1.0 - y) / y - 1.0       # Sh/2 - 1
for tgt in (0.01, 0.05, 0.10):
    yroot = brentq(lambda y: corr_excess(y) - tgt, 1e-12, 0.95, xtol=1e-14)
    M[f"stefan_y_for_{int(tgt*100):02d}pct"] = float(yroot)
    print(f"  radial convection raises Sh by {100*tgt:4.0f} % once the surface mole "
          f"fraction reaches y = {yroot:.6f}")
print("\\n  Placing this paper's runs on that axis needs vapour pressures, which")
print("  Part II prints only as Figure 14 -- out of scope. The correction is stated")
print("  as a function of y and not evaluated at their conditions.")''')

# ------------------------------------------------------ printed arithmetic
CODE.append('''# PRINTED ARITHMETIC, checked and not repaired.
c_sens = kp("test80_correction_sensible_heat")
c_rad  = kp("test80_correction_radiation")
c_cap  = kp("test80_correction_capillary")
n_app  = kp("test80_Nu_apparent")
n_cor  = kp("test80_Nu_corrected")
recomputed = n_app * (1 - c_sens - c_rad - c_cap)
M["test80_printed_vs_recomputed_rel"] = float(abs(recomputed - n_cor) / n_cor)
print(f"Part I, test No. 80:  {n_app:g} (1 - {c_sens:g} - {c_rad:g} - {c_cap:g}) = "
      f"{recomputed:.6f}")
print(f"  printed: {n_cor:g}.  Relative difference "
      f"{M['test80_printed_vs_recomputed_rel']:.2e} -- rounding only, the printed "
      f"value is right to its last digit.")
print(f"  (the three corrections together remove "
      f"{100*(c_sens+c_rad+c_cap):.1f} % of the apparent N_Nu, and the largest of "
      f"them, sensible heat, {100*c_sens:.1f} %.)")

# Eq. (24) against eqs. (21): consistent only for one Prandtl number.
print(f"\\nPart II eq. (24) vs eqs. (21): {K_AQ:g} = {K_FIT:g} x Pr^(1/3) requires "
      f"Pr = {PR_IMPLIED:.5f} for air.")
print(f"  Table 1's own 18 rows, intercept held at {A_TH:g}, give a coefficient of "
      f"{b_fixed:.4f} against that {K_AQ:g}: "
      f"{100*(b_fixed/K_AQ-1):+.2f} %.")
print(f"  THIS IS NOT A SECOND RESULT. K_recovered = b_fixed x K/K_AQ identically,")
print(f"  so b_fixed/K_AQ - 1 IS K_recovered/K - 1, which the heat section already")
print(f"  reported as {100*(M['t1_heat_K_recovered']/K_FIT - 1):+.2f} %. It is printed"
      f" again here because a")
print(f"  reader arriving at this section should not have to reconstruct it, and it")
print(f"  is deliberately NOT reported as a second metric -- the same rule that kept")
print(f"  pf_logmean_returns_2_rel out of the agreement file.")''')

# --------------------------------------------------------------- break table
CODE.append('''# BREAK TABLE. Every metric reported below needs at least one row here that
# moves it. Rebuilt for this page's physics -- none of it travels from another.
def heat_fit_intercept(nu_col, sq_col=None, include_still=False):
    s = t1f.sqRe.to_numpy() if sq_col is None else sq_col
    v = np.asarray(nu_col, float)
    if include_still:
        s = np.r_[s, 0.0]; v = np.r_[v, float(still.Nu_apparent)]
    return np.linalg.lstsq(np.column_stack([np.ones_like(s), s]), v, rcond=None)[0]

def re_broken(df1, w=None):
    """Re-run the whole reference-temperature analysis on a corrupted Table 1."""
    R = pd.concat([re_frame(df1).assign(tab=1), re_frame(t2).assign(tab=2),
                   re_frame(t4).assign(tab=4)], ignore_index=True)
    yy = np.log(R.Re / (R.D * R.v * R.P)).to_numpy()
    if w is None:
        w = float(minimize_scalar(
            lambda ww: float(np.sqrt(np.mean(re_fit(ww, R, yy)[2] ** 2))),
            bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-10}).x)
    n_, _, r_ = re_fit(w, R, yy)
    return {"w": w, "n": n_, "rms": float(100*np.sqrt(np.mean(r_**2))),
            "max": float(100*np.max(np.abs(r_)))}

rows = []
def brk(what, metric, base, broken):
    f = abs(broken / base) if base not in (0.0,) else np.inf
    rows.append((what, metric, base, broken,
                 f"x{f:.3g}" if np.isfinite(f) else "-"))

# --- transcription / Re law
t1_bad = t1.copy(); t1_bad.loc[t1_bad.run == 1, "Re"] = 195.1     # 159.1 misread
B = re_broken(t1_bad)
brk("Table 1 run 1 Re read 195.1 instead of 159.1", "re_law_resid_max_pct",
    M["re_law_resid_max_pct"], B["max"])
brk("Table 1 run 1 Re read 195.1 instead of 159.1", "re_law_resid_rms_pct",
    M["re_law_resid_rms_pct"], B["rms"])
brk("Table 1 run 1 Re read 195.1 instead of 159.1", "re_reference_temperature_weight",
    M["re_reference_temperature_weight"], B["w"])
brk("Table 1 run 1 Re read 195.1 instead of 159.1", "re_law_T_exponent",
    M["re_law_T_exponent"], B["n"])
t1_dp = t1.copy(); t1_dp["D_p_cm"] = 0.954       # decimal point lost on every row
Bd = re_broken(t1_dp)
brk("Table 1 drop diameter read 0.954 instead of 0.0954", "re_law_resid_max_pct",
    M["re_law_resid_max_pct"], Bd["max"])
t1_dt = t1.copy(); t1_dt["drop_temp_C"] = t1_dt["air_temp_C"]   # drop temps discarded
Bt = re_broken(t1_dt)
brk("Table 1 drop temperatures replaced by the air temperatures",
    "re_reference_temperature_weight", M["re_reference_temperature_weight"], Bt["w"])
brk("Table 1 drop temperatures replaced by the air temperatures",
    "re_law_resid_rms_pct", M["re_law_resid_rms_pct"], Bt["rms"])
brk("reference temperature forced to the free stream (w = 0)",
    "re_law_resid_rms_pct", M["re_law_resid_rms_pct"],
    M["re_law_resid_rms_pct_free_stream"])
brk("reference temperature forced to the free stream (w = 0)", "re_law_T_exponent",
    M["re_law_T_exponent"], M["re_law_T_exponent_free_stream"])
brk("free-stream residual taken over Table 1 only", "re_law_resid_rms_pct_free_stream",
    M["re_law_resid_rms_pct_free_stream"],
    float(100*np.sqrt(np.mean(r_air[(RE_ALL.tab == 1).to_numpy()]**2))))
brk("free-stream exponent refitted with the drop temperatures discarded",
    "re_law_T_exponent_free_stream", M["re_law_T_exponent_free_stream"], Bt["n"])
brk("worst-case film-vs-free-stream shift taken on Table 1 alone",
    "re_freestream_worst_Re_pct", M["re_freestream_worst_Re_pct"],
    float(100*np.max(np.abs(1/_ratio[(RE_ALL.tab == 1).to_numpy()] - 1))))
brk("worst-case film-vs-free-stream N_Nu shift taken on Table 1 alone",
    "re_freestream_worst_Nu_pct", M["re_freestream_worst_Nu_pct"],
    float(100*np.max(np.abs((_nu_alt/_nu_now)[(RE_ALL.tab == 1).to_numpy()] - 1))))
brk("the 5 % band on w taken at 50 % instead", "re_reference_weight_band",
    M["re_reference_weight_band"],
    float(brentq(lambda w: np.sqrt(np.mean(re_fit(w)[2]**2))
                 - 1.5*np.sqrt(np.mean(r_best**2)), W, 1.0)
          - brentq(lambda w: np.sqrt(np.mean(re_fit(w)[2]**2))
                   - 1.5*np.sqrt(np.mean(r_best**2)), 0.0, W)))
iso_bad = iso.copy(); iso_bad.iloc[0, iso_bad.columns.get_loc("air_vel_cm_s")] = 92.5
rv2 = (iso_bad.Re / iso_bad.air_vel_cm_s).to_numpy()
brk("Table 1 run 6 velocity read 92.5 instead of 95.2", "t1_iso_Re_over_v_spread_pct",
    M["t1_iso_Re_over_v_spread_pct"], float(100 * np.ptp(rv2) / rv2.mean()))
_rm_fs = T_iso ** n_air
brk("the nine rows' temperature effect evaluated at the free-stream exponent",
    "t1_iso_expected_spread_from_T_pct", M["t1_iso_expected_spread_from_T_pct"],
    float(100 * np.ptp(_rm_fs) / _rm_fs.mean()))
brk("rounding bound taken as the worst single row's half-width, not a peak-to-peak",
    "t1_iso_rounding_ptp_bound_pct", M["t1_iso_rounding_ptp_bound_pct"],
    float(_hw.max()))

# --- the reference temperature the paper states, and Table 3 held out
brk("the authors' 'average temperature' read as the free stream, not the film",
    "stillair_film_minus_stated_K", M["stillair_film_minus_stated_K"],
    float(still.air_temp_C + 273.15 - T_STATED))
brk("reference temperature forced to the free stream (w = 0)", "t3_heldout_rms_pct",
    M["t3_heldout_rms_pct"], M["t3_heldout_rms_pct_free_stream"])
brk("reference temperature forced to the drop surface (w = 1)", "t3_heldout_rms_pct",
    M["t3_heldout_rms_pct"], t3_rms(1.0))
brk("free-stream Table 3 residual evaluated at the recovered w instead",
    "t3_heldout_rms_pct_free_stream", M["t3_heldout_rms_pct_free_stream"],
    M["t3_heldout_rms_pct"])
_t3_flat = t3.copy(); _t3_flat["drop_temp_C"] = _t3_flat["air_temp_C"]
brk("Table 3's drop temperatures discarded and replaced by its air temperatures",
    "t3_heldout_w", M["t3_heldout_w"],
    float(minimize_scalar(lambda w: t3_rms(w, _t3_flat), bounds=(0.0, 1.0),
                          method="bounded", options={"xatol": 1e-10}).x))

# --- heat side
brk("apparent N_Nu used instead of corrected", "t1_heat_intercept_free",
    M["t1_heat_intercept_free"], float(heat_fit_intercept(t1f.Nu_apparent)[0]))
brk("apparent N_Nu used instead of corrected", "t1_heat_eq24_rms_pct",
    M["t1_heat_eq24_rms_pct"], rms_pct(nu_corr(t1f.Re), t1f.Nu_apparent))
brk("still-air row wrongly swept into the regression", "t1_heat_intercept_free",
    M["t1_heat_intercept_free"],
    float(heat_fit_intercept(t1f.Nu_corrected, include_still=True)[0]))
brk("Reynolds exponent 1/2 -> 1/3", "t1_heat_eq24_rms_pct", M["t1_heat_eq24_rms_pct"],
    rms_pct(A_TH + K_AQ * t1f.Re.to_numpy() ** (1/3), Nu))
brk("eq. (24) used with K = 0.60 (Pr^(1/3) forgotten)", "t1_heat_eq24_max_pct",
    M["t1_heat_eq24_max_pct"],
    float(np.max(np.abs(pct(A_TH + K_FIT * sq, Nu)))))
brk("null baselines evaluated on Table 3 instead of Table 1",
    "t1_heat_null_best_rms_pct", M["t1_heat_null_best_rms_pct"],
    float(min(rms_pct(np.full_like(Nu3, Nu3.mean()), Nu3),
              rms_pct(np.full_like(Nu3, A_TH), Nu3),
              rms_pct(K_AQ * np.sqrt(X3 / PR_IMPLIED ** M_PR), Nu3))))
brk("Reynolds exponent released with the intercept held at the still-air 2.23",
    "t1_heat_Re_exponent_fitted", M["t1_heat_Re_exponent_fitted"],
    float(curve_fit(lambda R, Kk, q: nu0 + Kk * R ** q, t1f.Re.to_numpy(), Nu,
                    p0=[K_AQ, N_RE])[0][1]))
brk("Prandtl exponent 1/3 -> 1/2 when inverting eq. (24)", "pr_air_implied_by_eq24",
    M["pr_air_implied_by_eq24"], float((K_AQ / K_FIT) ** 2))

# --- mass side
Nup_guess = t1f.Nu_prime.to_numpy().copy()
Nup_guess[~mass_ok] = [3.32, 3.05, 2.96]     # the blotted cells, filled by guess
sqg = t1f.sqRe.to_numpy()
brk("the three ink-damaged N_Nu' cells filled in by guesswork",
    "t1_mass_intercept_free", M["t1_mass_intercept_free"],
    float(np.linalg.lstsq(np.column_stack([np.ones_like(sqg), sqg]), Nup_guess,
                          rcond=None)[0][0]))
brk("benzene: apparent N_Nu used where N_Nu' is meant", "t4_benzene_intercept_free",
    M["t4_benzene_intercept_free"],
    float(np.linalg.lstsq(np.column_stack([np.ones_like(sqb), sqb]),
                          t4.Nu_apparent.to_numpy(), rcond=None)[0][0]))
brk("Table 4 run 11 read 10.6 -- BARELY MOVES; the ratio identity is what excludes it",
    "t4_benzene_intercept_free", M["t4_benzene_intercept_free"], float(ab_alt))
# THE ROW THE OLD BREAK TABLE DID NOT HAVE. Reading run 11 as 10.6 moves the benzene
# intercept by only 1.5 %, i.e. under check_agreement.py's tolerance, so no metric on
# the old page would have failed. The ratio spread does fail, by a factor of seven.
_r_bad = np.array([ratio[0], 10.6 / nu11, ratio[2]])
brk("Table 4 run 11's N_Nu' read 10.6 instead of 10.0",
    "t4_ratio_spread_runs_10_12_pct", M["t4_ratio_spread_runs_10_12_pct"],
    float(100 * np.ptp(_r_bad) / _r_bad.mean()))
brk("Table 4 run 12's N_Nu' read 10.0 instead of 10.9",
    "t4_ratio_spread_runs_10_12_pct", M["t4_ratio_spread_runs_10_12_pct"],
    float(100 * np.ptp([ratio[0], ratio[1], 10.0 / 8.7])
          / np.mean([ratio[0], ratio[1], 10.0 / 8.7])))
brk("run 11's bracket taken from run 1 (ratio 1.00) instead of runs 10 and 12",
    "t4_run11_Nuprime_upper_bound", M["t4_run11_Nuprime_upper_bound"],
    float((nu11 + 0.05) * _rbracket(t4[t4.run == 1].iloc[0])[1]))
brk("water and benzene slopes exchanged", "sc_ratio_benzene_over_water",
    M["sc_ratio_benzene_over_water"], float((bm_fixed / bb_fixed) ** (1 / M_PR)))
brk("still-air N_Nu' read 7.9 instead of 1.79", "stillair_analogy_ratio",
    M["stillair_analogy_ratio"], float(nu0 / 7.9))
brk("replicate pairs taken from runs 8 and 10 instead of 8 and 9",
    "replicate_worst_pct", M["replicate_worst_pct"],
    float(100 * abs(t1.loc[t1.run == 8, "Nu_corrected"].iloc[0]
                    - t1.loc[t1.run == 10, "Nu_corrected"].iloc[0])
          / (0.5 * (t1.loc[t1.run == 8, "Nu_corrected"].iloc[0]
                    + t1.loc[t1.run == 10, "Nu_corrected"].iloc[0]))))

# --- Table 3
brk("Table 3 abscissa raised to 2/3 (the two assumed exponents exchanged)",
    "t3_eq21_rms_pct", M["t3_eq21_rms_pct"],
    rms_pct(nu_from_abscissa(X3 ** (2/3)), Nu3))
brk("Table 3 predicted with the mass coefficient instead of the heat one",
    "t3_from_t1_rms_pct", M["t3_from_t1_rms_pct"],
    rms_pct(nu_from_abscissa(X3, Kc=bm_fixed / PR_IMPLIED ** M_PR), Nu3))
brk("Table 3 row 7 given row 9's diameter (0.101 for 0.088)",
    "t3_pair_abscissa_max_dev_pct", M["t3_pair_abscissa_max_dev_pct"],
    float(max(M["t3_pair_abscissa_max_dev_pct"],
              abs(100 * ((t3.Re_half_Pr_third[6] / t3.Re_half_Pr_third[7])
                         / np.sqrt(t3.D_p_corresponding_cm[8]
                                   / t3.D_p_corresponding_cm[7]) - 1)))))
brk("Table 3 N_Nu of row 1 read 4.0 instead of 7.0", "t3_pair_nu_ratio_max_dev_pct",
    M["t3_pair_nu_ratio_max_dev_pct"],
    float(max(abs(100 * ((4.0 / t3.Nu[1])
                         / (t3.neg_dDp2_dtau_cm2_s_x1e4[0]
                            / t3.neg_dDp2_dtau_cm2_s_x1e4[1]) - 1)),
              M["t3_pair_nu_ratio_max_dev_pct"])))
brk("Table 3 worst pair dropped", "t3_pair_worst_nu_ratio_error_pct",
    M["t3_pair_worst_nu_ratio_error_pct"],
    float(np.sort([abs(d) for *_, d in pair_rows])[-2]))
brk("sign convention flipped when counting which side is steeper",
    "t3_pairs_measurement_steeper", M["t3_pairs_measurement_steeper"],
    float(sum(1 for *_, d in pair_rows if d < 0)))

# --- pymrm
brk("construct_div nu = 1 (cylindrical) instead of 2", "pymrm_conduction_Nu_infinity",
    M["pymrm_conduction_Nu_infinity"], nu_shell(400, RATIO, nu=1))
brk("domain stopped at b/a = 10 and reported as the limit",
    "pymrm_conduction_dev_from_2", M["pymrm_conduction_dev_from_2"],
    float(abs(nu_shell(400, 10.0) - A_TH)))
brk("grid coarsened to n = 25 at b/a = 1e4", "pymrm_conduction_grid_order",
    M["pymrm_conduction_grid_order"],
    float(np.log(abs(nu_shell(12, RATIO) - 2/(1-1/RATIO))
                 / abs(nu_shell(25, RATIO) - 2/(1-1/RATIO))) / np.log(25/12)))
brk("resistance integral truncated at 10 radii instead of infinity",
    "nu_quadrature_route", M["nu_quadrature_route"],
    float(2.0 / (4*np.pi*quad(lambda r: 1/(4*np.pi*r**2), 1.0, 10.0)[0])))
brk("radial convection switched off (v = 0) at y = 0.20",
    "stefan_pymrm_vs_closed_rel", M["stefan_pymrm_vs_closed_rel"],
    float(abs(nu_shell(600, RS) / sh_closed(RS, 0.20) - 1)))
brk("p_f taken as the arithmetic mean", "stefan_pymrm_vs_closed_rel",
    M["stefan_pymrm_vs_closed_rel"], float(np.max(np.abs(back_arith / A_TH - 1))))
brk("Stefan grid order taken over the coarse pair n = 300/600 instead",
    "stefan_grid_order", M["stefan_grid_order"],
    float(np.log(abs(gr[0]-ge) / abs(gr[1]-ge)) / np.log(2.0)))
brk("p_f arithmetic-mean deviation evaluated at y = 0.05 instead of 0.20",
    "pf_arith_dev_at_y020_pct", M["pf_arith_dev_at_y020_pct"],
    float(100 * abs(sh_closed(RS, 0.05) * pf_arith(0.05) * (1 - 1/RS) / A_TH - 1)))
brk("Stefan threshold solved for a 2 % excess instead of 1 %",
    "stefan_y_for_01pct", M["stefan_y_for_01pct"],
    float(brentq(lambda y: corr_excess(y) - 0.02, 1e-12, 0.95)))
brk("Stefan threshold solved for a 7 % excess instead of 5 %",
    "stefan_y_for_05pct", M["stefan_y_for_05pct"],
    float(brentq(lambda y: corr_excess(y) - 0.07, 1e-12, 0.95)))
brk("Stefan threshold solved for a 15 % excess instead of 10 %",
    "stefan_y_for_10pct", M["stefan_y_for_10pct"],
    float(brentq(lambda y: corr_excess(y) - 0.15, 1e-12, 0.95)))
brk("the two routes to 2.0 compared at b/a = 10 rather than in the limit",
    "nu_two_routes_rel_diff", M["nu_two_routes_rel_diff"],
    float(abs(nu_shell(400, 10.0) - M["nu_quadrature_route"]) / A_TH))

# --- printed arithmetic and the derived quantities
brk("capillary correction dropped from test No. 80",
    "test80_printed_vs_recomputed_rel", M["test80_printed_vs_recomputed_rel"],
    float(abs(n_app * (1 - c_sens - c_rad) - n_cor) / n_cor))
brk("still-air N_Nu read 3.23 instead of 2.23 (leading digit)",
    "stillair_Nu_over_theory", M["stillair_Nu_over_theory"], 3.23 / A_TH)
brk("still-air N_Nu' read 1.19 instead of 1.79", "stillair_Nuprime_over_theory",
    M["stillair_Nuprime_over_theory"], 1.19 / A_TH)
brk("D_v back-out done against the heat number 2.23 instead of N_Nu'",
    "dv_used_over_dv_backed_out", M["dv_used_over_dv_backed_out"], A_TH / nu0)
brk("water mass slope taken from the free fit instead of the fixed-intercept fit",
    "t1_mass_slope_free", M["t1_mass_slope_free"], bm_fixed)
brk("heat slope taken with the intercept forced to zero", "t1_heat_slope_free",
    M["t1_heat_slope_free"],
    float(np.sum(sq * Nu) / np.sum(sq ** 2)))
brk("heat coefficient recovered with Pr = 1 instead of the paper's",
    "t1_heat_K_recovered", M["t1_heat_K_recovered"], b_fixed)
brk("benzene slope taken from the apparent N_Nu column", "t4_benzene_slope_free",
    M["t4_benzene_slope_free"],
    float(np.linalg.lstsq(np.column_stack([np.ones_like(sqb), sqb]),
                          t4.Nu_apparent.to_numpy(), rcond=None)[0][1]))
brk("water intercept deviation taken against 2.23 instead of the theoretical 2.0",
    "mass_intercept_dev_water_pct", M["mass_intercept_dev_water_pct"],
    float(100 * (am_free / nu0 - 1)))
brk("benzene intercept deviation taken against 1.79", "mass_intercept_dev_benzene_pct",
    M["mass_intercept_dev_benzene_pct"], float(100 * (ab_free / nup0 - 1)))
brk("heat intercept deviation taken against the still-air 2.23",
    "heat_intercept_dev_pct", M["heat_intercept_dev_pct"],
    float(100 * (a_free / nu0 - 1)))
brk("replicate pairs taken from Table 1 runs 8 and 11 instead of 8 and 9",
    "replicate_mean_pct", M["replicate_mean_pct"],
    float(np.mean([100 * abs(t1.loc[t1.run == 8, c].iloc[0]
                             - t1.loc[t1.run == 11, c].iloc[0])
                   / (0.5 * (t1.loc[t1.run == 8, c].iloc[0]
                             + t1.loc[t1.run == 11, c].iloc[0]))
                   for c in ("evap_rate_ml_s_x1e5", "Nu_apparent", "Nu_corrected",
                             "Nu_prime")])))
brk("Table 3 bias computed on the log of N_Nu", "t3_eq21_bias_pct",
    M["t3_eq21_bias_pct"],
    float(100 * np.mean(np.log(nu_from_abscissa(X3) / Nu3))))
brk("Table 3 worst row dropped", "t3_eq21_max_pct", M["t3_eq21_max_pct"],
    float(np.sort(np.abs(pct(nu_from_abscissa(X3), Nu3)))[-2]))
brk("Re residual rms taken over Table 1 only", "re_law_resid_rms_pct",
    M["re_law_resid_rms_pct"],
    float(100*np.sqrt(np.mean(r_best[(RE_ALL.tab == 1).to_numpy()] ** 2))))
Nu_bad = Nu.copy(); Nu_bad[t1f.run.to_numpy() == 18] = 6.56   # 9.56 misread
brk("Table 1 run 18's N_Nu read 6.56 instead of 9.56", "t1_heat_freefit_rms_pct",
    M["t1_heat_freefit_rms_pct"],
    rms_pct(np.polyval(heat_fit_intercept(Nu_bad)[::-1], sq), Nu_bad))
brk("Table 1 run 18's N_Nu read 6.56 instead of 9.56", "t1_heat_intercept_free",
    M["t1_heat_intercept_free"], float(heat_fit_intercept(Nu_bad)[0]))
brk("Table 1 run 18's N_Nu read 6.56 instead of 9.56", "t1_heat_slope_free",
    M["t1_heat_slope_free"], float(heat_fit_intercept(Nu_bad)[1]))
brk("Table 1 run 18's N_Nu read 6.56 instead of 9.56", "t1_heat_slope_at_A2",
    M["t1_heat_slope_at_A2"],
    float(np.sum(sq * (Nu_bad - A_TH)) / np.sum(sq ** 2)))
brk("Table 1 run 18's N_Nu read 6.56 instead of 9.56", "t1_heat_K_recovered",
    M["t1_heat_K_recovered"],
    float(np.sum(sq * (Nu_bad - A_TH)) / np.sum(sq ** 2)) / PR_IMPLIED ** M_PR)
brk("Table 1 run 18's N_Nu read 6.56 instead of 9.56", "heat_intercept_dev_pct",
    M["heat_intercept_dev_pct"], float(100 * (heat_fit_intercept(Nu_bad)[0] / A_TH - 1)))
Nup_bad = Nup.copy(); Nup_bad[-1] = 4.48                     # 9.48 misread
_ab, _bb = np.linalg.lstsq(np.column_stack([np.ones_like(sqm), sqm]), Nup_bad,
                           rcond=None)[0]
brk("Table 1 run 18's N_Nu' read 4.48 instead of 9.48", "t1_mass_slope_free",
    M["t1_mass_slope_free"], float(_bb))
brk("Table 1 run 18's N_Nu' read 4.48 instead of 9.48", "t1_mass_intercept_free",
    M["t1_mass_intercept_free"], float(_ab))
brk("Table 1 run 18's N_Nu' read 4.48 instead of 9.48", "mass_intercept_dev_water_pct",
    M["mass_intercept_dev_water_pct"], float(100 * (_ab / A_TH - 1)))
brk("uniform radial grid instead of geometric at b/a = 1e4",
    "pymrm_conduction_grid_order", M["pymrm_conduction_grid_order"],
    float(np.log(abs(nu_shell(200, RATIO, geometric=False) - 2/(1-1/RATIO))
                 / abs(nu_shell(400, RATIO, geometric=False) - 2/(1-1/RATIO)))
          / np.log(2.0)))
brk("Table 3 predicted with half the Table 1 coefficient", "t3_from_t1_rms_pct",
    M["t3_from_t1_rms_pct"],
    rms_pct(nu_from_abscissa(X3, Kc=0.5 * M["t1_heat_K_recovered"]), Nu3))
brk("intercept held at the still-air 2.23 instead of the theoretical 2.0",
    "t1_heat_slope_at_A2", M["t1_heat_slope_at_A2"],
    float(np.sum(sq * (Nu - nu0)) / np.sum(sq ** 2)))

BT = pd.DataFrame(rows, columns=["what is broken", "metric", "baseline",
                                 "broken value", "factor"])
# A row that moves a metric by less than check_agreement.py's 5 % tolerance would
# not be caught by CI either, so it does not count as coverage.
BT["moves_past_CI"] = [abs(bv - ba) > 0.05 * max(abs(ba), 1e-12)
                       for ba, bv in zip(BT.baseline, BT["broken value"])]
pd.set_option("display.width", 200, "display.max_colwidth", 62)
print(BT.to_string(index=False, float_format=lambda v: f"{v:.6g}"))''')

# --------------------------------------------------------- coverage + report
CODE.append('''# COVERAGE MAP, asserted key for key against what is reported.
covered = set(BT.metric)
strong  = set(BT.metric[BT.moves_past_CI])
missing = sorted(set(M) - covered)
weak    = sorted(set(M) - strong)
extra   = sorted(covered - set(M))
print(f"{len(M)} metrics reported; {len(covered)} have a break row; {len(strong)} have")
print(f"one that moves them past check_agreement.py's 5 % tolerance.")
print(f"{int(BT.moves_past_CI.sum())} of {len(BT)} rows clear that tolerance. The other "
      f"{len(BT) - int(BT.moves_past_CI.sum())} are")
print(f"kept but are NOT all deliberate 'barely moves' exhibits: exactly one says so in")
print(f"its own label (Table 4 run 11), and the rest are ordinary rows whose metric")
print(f"happens to be insensitive to that particular perturbation and is covered by")
print(f"another row. The label in the 'what is broken' column is what marks intent.")
assert not extra, f"break table names a metric that is not reported: {extra}"
assert not missing, f"metric without a break row: {missing}"
assert not weak, f"metric with no break row that moves it past CI tolerance: {weak}"

FLOOR = 1e-12
below = {k: v for k, v in M.items() if abs(v) < FLOOR}
print(f"\\nmetrics below check_agreement.py's ABS_FLOOR = {FLOOR:g} "
      f"(outside the regression suite): {below if below else 'none'}")

report_agreement(PAGE, M)''')


# =====================================================================
# Execute the code cells to obtain every quoted number.
# =====================================================================
def _run() -> dict:
    import matplotlib
    matplotlib.use("Agg")
    ns: dict = {}
    for src in CODE:
        if src.lstrip().startswith("try:\n    import pymrm"):
            continue                                   # magic-only Colab cell
        exec(compile(src, "<cell>", "exec"), ns)       # noqa: S102
    return ns


NS = _run()
M = NS["M"]


def f(x, n=2):
    return f"{x:.{n}f}"


V = {
    "n_t1": str(len(NS["t1f"])),
    "n_mass": str(int(NS["mass_ok"].sum())),
    "n_t4": str(len(NS["t4"])),
    "n_t3": str(len(NS["t3"])),
    "re_min": f"{NS['t1f'].Re.min():g}",
    "re_max": f"{NS['t1f'].Re.max():g}",
    "pr_implied": f(M["pr_air_implied_by_eq24"], 3),
    "heat_a": f(M["t1_heat_intercept_free"], 4),
    "heat_b": f(M["t1_heat_slope_free"], 4),
    "heat_b_fixed": f(M["t1_heat_slope_at_A2"], 4),
    "heat_K": f(M["t1_heat_K_recovered"], 4),
    "heat_K_pct": f"{100 * (M['t1_heat_K_recovered'] / 0.60 - 1):+.2f}",
    "re_sigma_air": f"{abs((M['re_law_T_exponent_free_stream'] + 1.7) / NS['se_air']):.0f}",
    "heat_q": f(M["t1_heat_Re_exponent_fitted"], 4),
    "heat_rms": f(M["t1_heat_eq24_rms_pct"], 2),
    "heat_max": f(M["t1_heat_eq24_max_pct"], 2),
    "heat_null": f(M["t1_heat_null_best_rms_pct"], 1),
    "mass_a": f(M["t1_mass_intercept_free"], 4),
    "mass_dev": f"{M['mass_intercept_dev_water_pct']:+.1f}",
    "benz_a": f(M["t4_benzene_intercept_free"], 4),
    "benz_dev": f"{M['mass_intercept_dev_benzene_pct']:+.1f}",
    "heat_dev": f"{M['heat_intercept_dev_pct']:+.1f}",
    "sc_ratio": f(M["sc_ratio_benzene_over_water"], 2),
    "still_ratio": f(M["stillair_analogy_ratio"], 4),
    "still_ratio_pct": f(100 * (M["stillair_analogy_ratio"] - 1), 1),
    "dv_factor": f(M["dv_used_over_dv_backed_out"], 3),
    "rep_mean": f(M["replicate_mean_pct"], 2),
    "rep_worst": f(M["replicate_worst_pct"], 2),
    "sigma_mass": f"{abs(M['mass_intercept_dev_water_pct']) / M['replicate_mean_pct']:.0f}",
    "t3_rms": f(M["t3_eq21_rms_pct"], 1),
    "t3_max": f(M["t3_eq21_max_pct"], 1),
    "t3_bias": f"{M['t3_eq21_bias_pct']:+.1f}",
    "t3_from_t1": f(M["t3_from_t1_rms_pct"], 1),
    "t3_factor": f(M["t3_eq21_rms_pct"] / M["t1_heat_eq24_rms_pct"], 1),
    "t3_pair_err": f(M["t3_pair_worst_nu_ratio_error_pct"], 1),
    "t3_steeper": str(int(M["t3_pairs_measurement_steeper"])),
    "n_pairs": str(int(len(NS["pair_rows"]))),
    "t3_pair_x": f(M["t3_pair_abscissa_max_dev_pct"], 2),
    "t3_pair_n": f(M["t3_pair_nu_ratio_max_dev_pct"], 2),
    "re_n": str(len(NS["RE_ALL"])),
    "re_exp": f(M["re_law_T_exponent"], 2),
    "re_w": f"{M['re_reference_temperature_weight']:.4f}",
    "re_exp_air": f(M["re_law_T_exponent_free_stream"], 2),
    "re_rms_air": f(M["re_law_resid_rms_pct_free_stream"], 2),
    "re_band": f"{M['re_reference_weight_band']:.3f}",
    "re_t1_rms": f"{100*np.sqrt(np.mean(NS['r_best'][(NS['RE_ALL'].tab == 1).to_numpy()]**2)):.3f}",
    "re_t1_max": f"{100*np.max(np.abs(NS['r_best'][(NS['RE_ALL'].tab == 1).to_numpy()])):.3f}",
    "re_fs_re": f(M["re_freestream_worst_Re_pct"], 1),
    "re_fs_nu": f(M["re_freestream_worst_Nu_pct"], 1),
    "re_rms": f(M["re_law_resid_rms_pct"], 2),
    "re_max": f(M["re_law_resid_max_pct"], 2),
    "iso_spread": f(M["t1_iso_Re_over_v_spread_pct"], 3),
    "nu_inf": f"{M['pymrm_conduction_Nu_infinity']:.9f}",
    "nu_dev": f"{M['pymrm_conduction_dev_from_2']:.1e}",
    "nu_order": f(M["pymrm_conduction_grid_order"], 2),
    "nu_routes": f"{M['nu_two_routes_rel_diff']:.1e}",
    "stefan_rel": f"{M['stefan_pymrm_vs_closed_rel']:.1e}",
    "stefan_order": f(M["stefan_grid_order"], 2),
    "pf_log": f"{NS['pf_log_dev']:.1e}",
    "b10_pct": f(100 * abs(NS["nu_shell"](400, 10.0) / 2.0 - 1), 1),
    "pf_arith": f(M["pf_arith_dev_at_y020_pct"], 2),
    "y1": f(100 * M["stefan_y_for_01pct"], 2),
    "y5": f(100 * M["stefan_y_for_05pct"], 2),
    "y10": f(100 * M["stefan_y_for_10pct"], 2),
    "test80": f"{M['test80_printed_vs_recomputed_rel']:.1e}",
    "eq24_pct": f"{100 * (M['t1_heat_K_recovered'] / 0.60 - 1):.2f}",
    "iso_T_spread": f(M["t1_iso_expected_spread_from_T_pct"], 3),
    "iso_round": f(M["t1_iso_rounding_ptp_bound_pct"], 3),
    "still_film_K": f"{NS['T_still_film']:.2f}",
    "still_dK": f"{M['stillair_film_minus_stated_K']:+.2f}",
    "T_stated": f"{NS['T_STATED']:g}",
    "t3ho_rms": f(M["t3_heldout_rms_pct"], 2),
    "t3ho_fs": f(M["t3_heldout_rms_pct_free_stream"], 2),
    "t3ho_drop": f(NS["t3_rms"](1.0), 2),
    "t3ho_w": f"{M['t3_heldout_w']:.4f}",
    "t3ho_Tlo": f"{NS['T_t3'].min():.0f}",
    "t3ho_Thi": f"{NS['T_t3'].max():.0f}",
    "fitTlo": f"{NS['T_fit'].min():.0f}",
    "fitThi": f"{NS['T_fit'].max():.0f}",
    "t4_spread": f(M["t4_ratio_spread_runs_10_12_pct"], 2),
    "t4_spread_bad": f"{100*np.ptp(NS['_r_bad'])/np.mean(NS['_r_bad']):.2f}",
    "t4_bracket_hi": f(M["t4_run11_Nuprime_upper_bound"], 2),
    "t4_bracket_lo": f"{(NS['nu11']-0.05)*NS['R_LO']:.2f}",
    "t4_run11": f"{float(NS['t4'].loc[NS['t4'].run == 11, 'Nu_prime'].iloc[0]):.1f}",
    "t3ho_gap": f"{NS['T_t3'].max() - NS['T_fit'].max():.0f}",
    "benz_alt": f(NS["ab_alt"], 4),
    "benz_drift": f"{100*abs(NS['ab_alt']/M['t4_benzene_intercept_free'] - 1):.1f}",
    "n_metrics": str(len(M)),
    "n_break": str(len(NS["BT"])),
    "n_strong": str(int(NS["BT"].moves_past_CI.sum())),
    "n_weak": str(len(NS["BT"]) - int(NS["BT"].moves_past_CI.sum())),
    "n_weak_other": str(len(NS["BT"]) - int(NS["BT"].moves_past_CI.sum()) - 1),
}


def sub(text: str) -> str:
    return re.sub(r"«(\w+)»", lambda m: V[m.group(1)], text)


# =====================================================================
# The notebook: markdown interleaved with the code cells above.
# =====================================================================
cells: list = []
ci = iter(CODE)


def C():
    cells.append(code(next(ci)))


cells.append(md(sub(r"""---
title: "Ranz-Marshall: the only constant this paper tests is the 2.0, and the authors' own water data miss it"
description: "Nu = 2 + 0.6 Re^(1/2) Pr^(1/3) as Ranz and Marshall published it in 1952, taken apart constant by constant. The 0.60 is fitted to the very tables any comparison uses; the 1/2 and 1/3 are Froessling's, assumed and not testable from the transcribed tables; the 2.0 is the one piece of theory, and the paper's own data extrapolate to «heat_a» for heat, «benz_a» for benzene mass transfer and «mass_a» for water mass transfer. The authors saw the gap and named its likeliest cause; this page puts a number on it, against their own replicate scatter and against the claim in their abstract to have verified the intercept. pymrm recomputes the 2.0 from the spherical BVP it comes from, and puts a number on the radial-convection factor Ranz said in 1993 the papers had underplayed."
categories: [sec:A, struct:S3, tier:T0, data:tier3, phase:gas-liquid]
date: 2026-08-07
---

# Ranz-Marshall: the only constant this paper tests is the 2.0

**Catalog ID:** `A3.5` · **Structures:** `S3` (pointwise algebra, plus a 1-D spherical BVP) · **Tier:** T0

$$N_{Nu} = 2.0 + 0.60\,N_{Pr}^{1/3}N_{Re}^{1/2},
\qquad N_{Nu}' = 2.0 + 0.60\,N_{Sc}^{1/3}N_{Re}^{1/2}$$

Two of the most-used equations in chemical engineering, and they contain four
numbers with four different epistemic statuses:

| constant | status in the papers |
|---|---|
| `2.0` | **theoretical.** Part I derives it (eq. 7) from the conduction/diffusion field around a sphere at rest. Nothing was fitted. |
| `0.60` | **fitted**, to Figures 6, 7 and 9 — that is, to Tables 1, 2, 3 and 4, the same tables anyone would compare against. |
| `1/2` | **assumed.** Froessling's boundary-layer form, carried over. |
| `1/3` | **assumed.** Same — and not identifiable from Tables 1–4, every run of which is in air. Part II tests it graphically in Figure 12, which this page does not digitise. |

So a page that reports "the correlation reproduces the paper's data to a few per
cent" has reported a goodness of fit and nothing else. This page keeps the four
apart and asks what the data can actually decide.

**What comes out.** The intercept is the only testable constant, and the paper's
data test it three separate ways — an outcome the paper reports in words rather than
in numbers. Extrapolated to
$N_{Re}=0$ by a free two-parameter fit, Table 1's «n_t1» heat-transfer runs give
**«heat_a»** («heat_dev» %), Table 4's «n_t4» benzene mass-transfer runs give
**«benz_a»** («benz_dev» %) — the same answer to within the noise — and Table 1's
water mass-transfer runs give **«mass_a»**, «mass_dev» %, which is «sigma_mass»
times the run-to-run scatter the paper's own replicate pairs show. The still-air
row says the same thing more bluntly: it prints $N_{Nu}=2.23$ and $N_{Nu}'=1.79$
for a situation in which eqs. (21) and (22) both say 2.0, a ratio of
**«still_ratio»** where the correlation demands exactly 1.

**The authors saw this, and said so — qualitatively.** Part II, folio 173, first
column:

> "The points for heat transfer extrapolate close to the theoretical minimum value.
> Experimental values obtained for $N_{Nu}$ at $N_{Re} = 0$ were slightly greater
> than 2.0, and the small difference is attributed to free convection. … **Data for
> mass transfer show a steeper slope and a lower intercept**, but the disagreement
> is always less than 10 per cent at a given $N_{Re}$ and is usually much less than
> 5 per cent. Factors which may have contributed to errors in the data and anomalies
> in the heat- and mass-transfer analogy were: **(1) inaccurate values of
> diffusivity**; (2) $p_{Ai}$ may not have been the saturation vapor pressure; (3)
> the partial pressure of water vapor in the air may not have been zero at low air
> rates …"

So the direction, and the cause this page finds most consistent with the numbers,
are the authors' own. **What is added here is the arithmetic they did not do**, and
it changes what the sentence means:

- the deficit is **«mass_dev» %** on the intercept, not "less than 10 per cent" —
  their bound is about the *curves* at a given $N_{Re}$, where the fitted 0.60
  absorbs most of it, not about the extrapolated intercept;
- «sigma_mass» times their own replicate scatter, so it is not "small";
- it is species-specific: benzene, whose diffusivity came from Hirschfelder, Bird
  and Spotz rather than from the authors' own Figure 5, lands where the *heat* data
  land; only water misses;
- $N_{Nu}'\propto 1/D_v$, so closing the water gap needs a diffusivity smaller by
  «dv_factor» — which is exactly the 0.204 sq.cm/sec Part II prints, and calls "a
  low value compared with other methods of determination".

And that puts the numbers against two claims the papers do make. Part I's abstract
says the study "confirmed the analogy between heat and mass transfer at low Reynolds
numbers, and **verified the simple expression for the Nusselt number at zero Reynolds
number**" — the paper's own still-air pair has an analogy ratio of «still_ratio».
Folio 174 says the correlation's success carries "an implication that the calculated
value of the diffusivity of water vapor in air may be **more accurate than any
reported in the literature**" — while the same data, read at the intercept, say that
diffusivity is «dv_factor» too large. The gap is not an oversight of the authors'.
The published *conclusions* drawn from it are what these numbers sit against.

**And pymrm.** The 2.0 is a boundary-value problem, so it is computed here rather
than quoted: a spherical shell with `construct_div(nu=2)`, refined in grid and in
domain size, giving «nu_inf» — «nu_routes» relative from the exact answer, which
adaptive quadrature of the resistance integral supplies. The same solver,
with `construct_convflux_upwind` added, answers the question Ranz raised in 1993
about the radial-convection factor — and shows which of the two readings of the
paper's own definition of $p_f$ is the right one.""")))

cells.append(md(r"""## Background

Ranz and Marshall published "Evaporation from Drops" in two parts in *Chemical
Engineering Progress* 48 in 1952, out of W. R. Marshall's spray-drying programme
at Wisconsin. One drop, about a millimetre across, hangs from a microburet in an
upward air stream; liquid is fed in fast enough to hold the diameter constant,
and the feed rate is the evaporation rate. A 0.5-mil thermocouple inside the drop
gives its temperature. That is the whole experiment, and it is why the paper is
still read: it is a *single-particle* measurement in a field that otherwise had
to infer transfer coefficients from beds and sprays.

**The page-range discrepancy, recorded both ways.** Part I's scan frames carry two
folio lines — one above the black edge of the leaf, one inside it at the foot —
differing by exactly two. Four of the six frames show both legibly (139/141,
140/142, 141/143 and 144/146); on the fourth frame the upper line is cut off at the
leaf edge and on the fifth there is no upper line at all, so those two pairs are
inferred from the sequence rather than read. Part II settles which belongs to the
page: its first page has a *blank*
top margin and prints `Vol. 48, No. 4 | Chemical Engineering Progress | Page 173`
at the **foot** of the text block. In this issue the folio is a footer. The line
at the top of a Part I frame is therefore the preceding leaf's footer, caught in
the same frame, and **Part I's own folios are 141–146** — which is the universal
citation and the one Ranz himself gives in his 1993 ISI *Citation Classic*
("Chem. Eng. Progr. 48:141-6; 173-80, 1952"). The 139–144 reading on file comes
from reading the upper line as a running head. **Every page reference on this
page is the folio printed at the foot of the page the value appears on.** Part II
is 173–180 either way.

**Neither file has a text layer at all** — `pdftotext` returns one byte per page —
so every digit here was read off 400 dpi renders (the files' native resolution:
`pdfimages -list` reports JPEG RGB tiles at 400 × 400 ppi, so rendering larger
only interpolates), with each numeric column cropped and magnified on its own.

**Two files on disk are not this paper** and are named here so nobody re-opens
them hoping: `MISC-Ranz1993-citation-classic-commentary-CurrContents22.pdf` is a
one-page 1993 reminiscence carrying none of the correlations, and
`Charlesworth1960-evaporation-drops-dissolved-solids-AIChEJ6-9.pdf` is a
companion study about crust formation in drops containing dissolved solids. The
1993 file is quoted once on this page, for a caveat only its author could give."""))

cells.append(md(sub(r"""## The published model

**Part I, eq. (7) — the theory.** With the fluid at rest ($v_\theta = v_r = 0$),
eqs. (1), (2), (5) and (6) collapse to

$$N_{Nu} = N_{Nu}' = 2.0. \tag{7}$$

That is steady conduction (or diffusion) from a sphere into an unbounded stagnant
medium, and it is the only piece of the correlation that is derived.

**Part I, eqs. (8) and (9) — the form, from Froessling.**

$$N_{Nu}' = 2.0 + K_1 (N_{Sc})^m (N_{Re})^n, \qquad
  N_{Nu}  = 2.0 + K_2 (N_{Pr})^p (N_{Re})^q$$

with, in the paper's own words, "$K_1 = K_2$, $p = m = 1/3$, and $q = n = 1/2$."
Ranz went to Pittsburgh to read the only copy of Froessling's 1938 paper in the
United States and translated it himself (he tells that story in the 1993
commentary, not in the papers); the exponents are Froessling's, adopted here, not
measured here.

**Part II, eqs. (21) and (22) — the constant.**

$$N_{Nu} = 2.0 + 0.60\,N_{Pr}^{1/3}N_{Re}^{1/2}, \qquad
  N_{Nu}' = 2.0 + 0.60\,N_{Sc}^{1/3}N_{Re}^{1/2}$$

**Part II, eq. (24) — and a constant the papers do not flag as such.**

$$h = h_c = \frac{k}{D_p}\left(2.0 + 0.54\,N_{Re}^{1/2}\right),
\qquad\text{"restricted to aqueous drops".}$$

`0.54` is `0.60` times $Pr^{1/3}$ for air, so the three printed numbers together
pin the Prandtl number the authors used: $(0.54/0.60)^3 = $ **«pr_implied»**. That is
a **definition, not a check** — nothing independent is compared with it, and nothing
can be, because this page's own rule forbids importing a property value. It is worth
having anyway: it means the whole heat-transfer analysis can be done without
importing one.

**How the Nusselt numbers were measured** (Part I, eqs. 16, 19, 20). At finite
velocity from the feed rate; in still air from the slope of $D_p^2$ against time,

$$N_{Nu} = -\tfrac14\frac{\lambda_v\rho_l}{\Delta t\,k}\frac{\mathrm{d}D_p^2}{\mathrm{d}\tau},
\qquad
N_{Nu}' = -\tfrac14\frac{\rho_l M_m p_f}{\Delta p_A D_v \rho_{gm}}\frac{\mathrm{d}D_p^2}{\mathrm{d}\tau}.$$

The $p_f$ in that second expression is the object Ranz singled out forty-one
years later, and Part II's notation defines it only as the "average value of
$(\pi - p_A)$ across transfer path" — without saying *which* average. The
validation section settles it.""")))

cells.append(md(r"""## Parameters and assumptions

**Theirs.** Quasi-steady transfer around an isolated sphere; the drop surface at
a single temperature $t_i$ and saturated at $p_{Ai} = p_{Ai}(t_i)$; the transfer
path of uniform temperature and velocity; the drop spherical (Part II notes the
diameter along the capillary axis was about 8 % greater than across it, and $D_p$
is the arithmetic mean of the two); radiation, sensible heat of the feed and
conduction along the capillary removed as corrections (Part I eq. 18); no
thermal diffusion or Dufour effect ("negligible in this study").

**Ours.** Nothing added. Every transfer number used here is the paper's own
$N_{Nu}$, $N_{Nu}'$ and $N_{Re}$ column; no property of air, water or benzene is
imported from outside the two papers, and the one Prandtl number that appears is
back-calculated from the paper's own eq. (24). The pymrm sections use a unit
diffusivity and a unit sphere, because both answers are dimensionless.

**What is out of scope, and why.** Figures 5 (transport properties), 12 (the
five-decade master correlation), 13 (free convection) and 14 (vapour pressure over
NH₄NO₃ solutions) carry content that exists only as curves. No maintainer figure
review is available for this case, so **no point is read off any of them** and no
number on this page comes from one. Figure 12 is *named* twice below, once for what
its printed legend lists and once for the sentence on folio 174 that comments on it;
that is reading type, not reading a curve, and nothing numerical follows from it.
Three consequences are stated rather than worked around: the
free-convection exponent $1/4$ of eqs. (10) and (11) is **never exercised on this
page**, because the only zero-Reynolds datum in any table is the single still-air
row; Table 6's calculated evaporation rates cannot be reproduced, because
they need the vapour-pressure curve of Figure 14; and the paper's own test of the
Prandtl exponent, which lives entirely in Figure 12, is neither reproduced nor
contradicted here. Table 5 (drop temperature
against Reynolds number) is not used either — it carries no transfer coefficient."""))

cells.append(md(sub(r"""## The data

Five files, all transcribed from the papers' own tables — **tier 3, table
transcription, no digitisation anywhere.**

| file | what it is | rows |
|---|---|---|
| `…-table1.csv` | Part II Table 1, water drops in dry air, folio 174 | 18 forced + 1 still-air |
| `…-table2.csv` | Part II Table 2, water drops in 66–90 °C air, folio 175 | 9 |
| `…-table3.csv` | Part II Table 3, water drops in 85–221 °C air, from the motion-picture record, folio 175 | 9 |
| `…-table4.csv` | Part II Table 4, **benzene** drops, folio 176 | 13 |
| `…-printed-constants.csv` | every constant either part prints, with its role | 24 |

Four things about these rows that the analysis has to respect.

1. **The 0.60 was fitted to all of them.** Figures 6, 7 and 9 are exactly Tables
   1–4. **There is no held-out set for the correlation anywhere in either part**,
   and this page never claims one: where it needs an out-of-range test it refits the
   coefficient on Table 1's room-temperature runs alone and then applies it to
   Table 3. (There *is* a held-out set for the reference-temperature convention —
   Table 3 again, which prints no Reynolds column and so cannot enter that fit. It
   is used below, and its scope is stated where it is used.)
2. **Three $N_{Nu}'$ cells are unreadable.** An ink blot on the microfilm covers
   the mass-transfer entries of Table 1 runs 12, 13 and 14. They are left empty,
   not guessed; the mass-transfer fits use the «n_mass» rows that remain, and the
   break table shows what filling them in by eye would do.
3. **One digit in Table 4 is settled by arithmetic, not by the glyph.** Run 11's
   $N_{Nu}'$ is degraded on the film. It is recorded as **«t4_run11»**, and the reason is
   in the mass-transfer cell below: both of Table 4's transfer numbers come from the
   same measured rate, so their ratio depends only on the printed temperatures and
   pressure, and the neighbouring rows bracket run 11 tightly enough to exclude
   every other reading. An earlier draft read it 10.6 off the glyph's counter
   position and was wrong.
4. **Table 4 has no corrected $N_{Nu}$ column at all**, because for benzene the
   sensible-heat term is about 20 % of the heat flow and the authors did not
   attempt one. Apparent and corrected Nusselt numbers are never mixed on this
   page.

The Reynolds column is the transcription check, and it turns out to be more than
that. The tables define $N_{Re} = D_p v_o\rho/\mu$ and print $D_p$, $v_o$, the air
and drop temperatures and the pressure beside it, so it is over-determined:
«re_n» printed Reynolds numbers must lie on one $\rho/\mu(T,P)$ law. **Neither part
says at what temperature the air properties in those columns were evaluated** — the
one place either part states a reference temperature is for the single still-air run,
and that statement is reconciled below with what the columns give — so the next cell
leaves $T_{\rm ref}$ as an unknown and lets the 36 rows answer it.""")))

C()   # colab
C()   # imports
C()   # data
cells.append(md("""### The constants, with their roles"""))
C()   # constants

cells.append(md(sub(r"""### The transcription check

A single two-parameter law $\rho/\mu = C\,P\,T_{\rm ref}^{\,n}$ against every printed
Reynolds number in Tables 1, 2 and 4 that has a measured drop temperature, with the
reference temperature written $T_{\rm ref} = T_a + w\,(T_d - T_a)$ and $w$ left free.
Then a nearly fit-free companion: nine rows of Table 1 sit at one pressure inside
0.7 °C of each other in air temperature, so their nine $Re/v$ values are very nearly
one number — and the cell says how nearly, because *nearly* is the whole content of
the check.""")))
C()   # Re check

cells.append(md(sub(r"""**Minimising the residual over $w$** — by bounded Brent, not by sweeping — puts the
reference temperature at $w = $ **«re_w»**, the arithmetic mean of the air and the
drop temperatures. That is not a small preference. At the film temperature the fitted
exponent is **«re_exp»** — what kinetic theory gives for air, since
$\mu\sim T^{0.7}$ and $\rho\sim P/T$ make $\rho/\mu\sim P\,T^{-1.7}$ — and the
residual over all «re_n» rows is «re_rms» %; at the free stream the exponent is
«re_exp_air», «re_sigma_air» standard errors away from $-1.7$, and the residual is
«re_rms_air» %. **Table 1's own 18 rows come back to «re_t1_rms» % rms and
«re_t1_max» % worst** once the film temperature is used: eighteen four-figure
Reynolds numbers reproduced from four other printed columns by a two-parameter law.

**And the papers state it once themselves, for the one run where they had to.**
Part II's still-air paragraph gives $D_v = 0.204$ sq.cm/sec "at an average
temperature of «T_stated»° K. and a pressure of 741 mm. Hg." Table 1's still-air row
prints the air and drop temperatures that go with it, and their mean is
**«still_film_K» K** — «still_dK» K from the stated value, with the stated pressure
equal to the row's own printed pressure. So for the single run whose reference
temperature the authors write down, the average they took *is* the film temperature.
The recovered $w$ is not only a minimisation; it agrees with the one printed
statement of the convention in either part.

That leaves the *tables'* convention genuinely unstated — no part says at what
temperature the Reynolds columns of Tables 1, 2 and 4 were evaluated — and it
carries into anyone's reuse of the correlation: evaluating $N_{Re}$ at the free
stream instead of the film would move it by up to «re_fs_re» % on these very rows,
and the $N_{Nu}$ that eq. (24) then predicts by up to «re_fs_nu» %.

The identification has a width and the page states it: the band of $w$ within 5 % of
the best residual is «re_band» wide, so a $0.45$ or a $0.55$ rule would not be
distinguishable from $1/2$. What *is* distinguishable is the free stream, and the
drop surface.

**A correction this page made to itself.** An earlier draft called the $Re/v$ check
"fit-free and needing no reference temperature", on the grounds that nine rows at one
pressure and one air temperature must give one number. They do not. Their *drop*
temperatures spread the film temperature over about a kelvin, which at the fitted
exponent is a **«iso_T_spread» %** effect on $\rho/\mu$ — on its own larger than the
whole «iso_spread» % spread observed. And the observed quantity is a nine-row
peak-to-peak, so the bound it must be judged against is a peak-to-peak bound,
«iso_round» % from rounding alone, not one row's half-width. Read correctly the
check passes comfortably and is what it always was: a **gross-transcription
detector**, not a reference-temperature-free identity.

### The convention has an out-of-sample test, on rows the fit could not see

Table 3 prints no Reynolds column at all, so none of its rows can enter the fit
above. It does print $D_p$, $v$, both temperatures and (in its footnote) the
pressure, so the recovered law can be asked to *predict* the abscissa Table 3 does
print — at reference temperatures «t3ho_Tlo»–«t3ho_Thi» K against the
«fitTlo»–«fitThi» K the law was fitted on. The two ranges do not overlap.""")))
C()   # table 3 held out

cells.append(md(sub(r"""**«t3ho_rms» % rms** at the recovered $w$, against **«t3ho_fs» %** at the free
stream and **«t3ho_drop» %** at the drop surface; and the $w$ that best reproduces
Table 3 *on its own* is **«t3ho_w»**, against the «re_w» the other three tables give.
That is out-of-sample confirmation of the reference-temperature convention on nine
rows the fit never saw, «t3ho_gap» K beyond every temperature it did see — and it
simultaneously confirms Table 3's transcription.

Two limits on what that buys, stated because the page is careful about fit-versus-test
everywhere else. It is a test of the **reference-temperature convention only**: the
0.60 was fitted to Table 3 through Figure 9, so nothing here is out of sample for the
*correlation*. And it inherits the assumption that $N_{Pr}$ for air is the constant
«pr_implied» that eq. (24) implies, which is the same assumption the abscissa itself
carries.""")))

cells.append(md(r"""## PyMRM implementation

The algebra above needs no solver. What does need one is eq. (7) — the `2.0` —
which is a one-dimensional spherical boundary-value problem, and the radial
convection Ranz flagged in 1993, which is the same problem with a convective flux
added. Both are built with `construct_grad`, `construct_div(nu=2)` and
`construct_convflux_upwind`, on a geometric radial grid because the answer is set
by the shell nearest the drop while the far field has to reach many drop radii.

Those two cells are further down, after the results, because the results are what
the paper is about. The implementation notes worth having in advance:

- `construct_div(..., nu=2)` is the spherical geometry. Getting it wrong is the
  first row of the break table.
- Both boundaries are Dirichlet, so `{"a": 0, "b": 1, "d": …}` at each end and the
  outward-normal sign convention never bites.
- The surface flux is read from the pymrm face gradient (`grad @ C + grad_bc`),
  which is second order at a Dirichlet boundary — a hand-written one-sided
  difference would be first order and would masquerade as a physics error.
- The Stefan problem's total molar flow is closed by a **scalar root-find**
  (`brentq` on "the flux the solve produces equals the flux it was given"), not
  by a sweep."""))

cells.append(md(sub(r"""## Results

### Heat transfer — the fit, and the part of it that was not free

The 0.60 was fitted to these rows, so the rms below is a goodness of fit and is
labelled one. The intercept is the exception, and it is reported separately.""")))
C()   # heat

cells.append(md(sub(r"""Three things to take from that.

The coefficient these 18 rows give, taken back through the paper's own Prandtl
number, is **«heat_K»** against the printed 0.60, «heat_K_pct» %. Table 1 alone
reproduces the constant the authors fitted to four tables.

Released, the Reynolds exponent comes out **«heat_q»** against Froessling's
assumed 1/2. So of the two assumed exponents, the one on $N_{Re}$ *is* supported
by these data. The one on $N_{Pr}$ **cannot be identified from the transcribed
tables**: every heat-transfer run in Tables 1, 2 and 3 is in air, so
$N_{Pr}^{1/3}$ is a constant absorbed into the 0.60. That is stated here, not
buried — it is the exponent a reader is most likely to assume was measured.

Part II does test it, and this page does not. Figure 12's legend, under **HEAT
TRANSFER DATA**, lists `AIR TO WATER (PRESENT WORK)`, `STEEL TO OIL (17)`,
`STEEL TO WATER (17)` and `STEEL TO AIR (17)` — Kramers's spheres in oil and in
water, i.e. three decades of non-air Prandtl numbers, plotted on the
$N_{Pr}^{1/3}N_{Re}^{1/2}$ abscissa; and folio 174 concludes that "the over-all
success of this correlation gives credence to the values of the exponents on
$N_{Pr}$, $N_{Sc}$, and $N_{Re}$". That test is **graphical**, it lives in a figure
this page deliberately does not digitise, and so it is neither reproduced nor
contradicted here. The honest statement is: *untested by anything transcribed on
this page; tested graphically in the paper's own Figure 12.*

And the fit is worth something: «heat_rms» % rms against a best null of
«heat_null» %.

### Mass transfer — where the intercept stops working""")))
C()   # mass

cells.append(md(sub(r"""**This is the page's headline.** The same free two-parameter fit, applied to
three sets of runs from the same apparatus:

| set | runs | intercept at $N_{Re}=0$ | deviation from the theoretical 2.0 |
|---|---|---|---|
| heat, water (Table 1) | «n_t1» | **«heat_a»** | «heat_dev» % |
| mass, benzene (Table 4) | «n_t4» | **«benz_a»** | «benz_dev» % |
| mass, water (Table 1) | «n_mass» | **«mass_a»** | **«mass_dev» %** |

Two of the three land on the same answer, about «heat_dev» % above a constant the
paper derived rather than fitted — and they are the two that share nothing except
the apparatus: a heat measurement and a mass measurement, water and benzene, the
authors' own conductivity and someone else's diffusivity. The third misses by
«mass_dev» %, and the next cell shows the paper's own still-air measurement missing
it by very nearly the same amount in the same direction.

**The one degraded digit in Table 4, and why it is not decided by looking at it.**
The cell above settles run 11's $N_{Nu}'$ from Table 4's own arithmetic rather than
from the film. Both printed transfer numbers of that table are reduced from the same
measured evaporation rate, so their ratio is rate-independent and depends only on the
printed air temperature, drop temperature and pressure — and runs 10, 11 and 12 are
printed within 0.2 K of one another at one pressure. Runs 10 and 12 therefore bracket
run 11's $N_{Nu}'$ in **[«t4_bracket_lo», «t4_bracket_hi»]**; the first two glyphs
read "10." and the last is as wide as this column's round digits and not as narrow as
its 1s; so the value is **«t4_run11»**. The three ratios then agree to «t4_spread» %.

This matters as a method point, not only as a digit. Reading it 10.6 — which an
earlier draft did, from the counter position of a three-pixel white hole — moves the
benzene intercept to «benz_alt», a drift of «benz_drift» %. That is *under*
`check_agreement.py`'s 5 % tolerance, so no regression test on this page would have
caught it. The ratio spread would have gone «t4_spread» % → **«t4_spread_bad» %**,
which is why it is now reported as a metric in its own right and carries its own
break row.

The Schmidt exponent is worth a sentence of its own. Benzene is the only second
species in either part, so it is the only thing that touches the $1/3$ at all:
the two fitted slopes map onto a Schmidt-number ratio of **«sc_ratio»**
(benzene/water in air). That is a *consequence* of assuming the exponent, not a
test of it — closing it would need the two diffusivities, which the papers print
only in Figure 5. One species pair cannot separate a coefficient from an
exponent, and this page does not pretend otherwise.

### The still-air row, and the analogy the correlation asserts""")))
C()   # still air

cells.append(md(sub(r"""At $N_{Re}=0$ eqs. (21) and (22) are the same equation, so the paper's own
still-air pair should have a ratio of exactly 1. It is **«still_ratio»** —
«still_ratio_pct» % apart — against a run-to-run repeatability of «rep_mean» %
mean and «rep_worst» % worst, measured on the paper's three replicate pairs.

The paper explains the heat side itself: at $N_{Re}=0$ free convection is present
and lifts $N_{Nu}$ above 2.0. That explanation survives here, and gains support it
did not have — the *forced*-convection runs extrapolate to «heat_a», so the excess
really is a still-air effect and not a bias in the heat measurement.

The mass side is not a still-air effect either: the forced-convection runs carry it
too («mass_a»). Here the authors give the cause but not the size. Their list of
"factors which may have contributed to errors in the data and anomalies in the heat-
and mass-transfer analogy" opens with "inaccurate values of diffusivity", and every
ingredient needed to price that is printed. $N_{Nu}'\propto 1/D_v$, so the
diffusivity would have to be smaller by «dv_factor»; Part II prints exactly that
number (0.204 sq.cm/sec at «T_stated» K and 741 mm Hg) as what the still-air data
give *if* $N_{Nu}'=2.0$ is imposed, and calls it low against other determinations;
and Part I says its calculated water-vapour diffusivities are "approximately 10 per
cent lower than that indicated by the S.T.P. value given by the International
Critical Tables".

What is not in either part is the arithmetic — and with it, the conflict. The
authors' own bound on the anomaly is "always less than 10 per cent at a given
$N_{Re}$", which is a statement about the two *curves*, where the fitted 0.60 absorbs
most of the offset; at the intercept, where nothing was fitted, the gap is
«mass_dev» %, or «sigma_mass» times their replicate scatter. And where they conclude
on folio 174 that the correlation's success implies their calculated $D_v$ "may be
more accurate than any reported in the literature", the intercept says that same
$D_v$ is «dv_factor» too large. Both readings cannot stand; this page does not
claim to say which is right, only that the paper's own tables decide against the
stronger of the two published claims.

### Table 3: the high-temperature rows, and what they cost the correlation""")))
C()   # table 3

cells.append(md(sub(r"""Table 3 is inside Figure 9 and so inside the fit, and it still comes out at
**«t3_rms» % rms and «t3_max» % worst**, with a «t3_bias» % bias — a factor
«t3_factor» on Table 1's «heat_rms» %. Refitting the coefficient on Table 1's
room-temperature runs alone and predicting Table 3 changes almost nothing
(«t3_from_t1» %), which says the disagreement is not about the constant.

The within-pair identities say where it *is*. Table 3's nine rows are **four**
motion-picture runs read at two drop diameters each (pairs A–D) **plus one single
row** (E, which the pairing loop skips), so within a pair the air velocity, air
temperature and drop temperature are held fixed and only $D_p$ moves. Both identities that follow — the abscissa scaling as $D_p^{1/2}$, and
$N_{Nu}$ scaling as the printed $-\mathrm{d}D_p^2/\mathrm{d}\tau$ — hold to
«t3_pair_x» % and «t3_pair_n» %, well inside what two-figure printing allows. The
table is self-consistent. What fails is the correlation: inside a single run, at
fixed everything-else, the measured $N_{Nu}$ ratio
differs from what $2.0+0.60X$ predicts by up to **«t3_pair_err» %**, with the
measurement the steeper of the two in «t3_steeper» of «n_pairs» pairs — against
roughly 1.3 % on a predicted ratio and 3.4 % on a measured one from the two-figure
printing alone.""")))
C()   # figure

cells.append(md(sub(r"""## Validation

### 1. The theoretical 2.0, recomputed — and checked against the exact answer

Part I's eq. (7) is asserted, not solved, and it is the only constant in the
correlation worth solving for. pymrm solves the spherical BVP. The
series-resistance integral $\int_a^\infty \mathrm{d}r/4\pi r^2 k$ is then evaluated
by adaptive quadrature, sharing no grid, no assembly and no linear solve — but it is
worth being exact about what that second number is, because "computed two ways" is
easy to oversell: $2/(4\pi\int_1^\infty \mathrm{d}r/4\pi r^2)$ is **identically 2**
for any working quadrature. It is the *exact* answer, obtained numerically. So what
follows is a finite-volume solve measured against the analytic limit, not two
independent estimates that happen to agree — which is the stronger of the two things
to have, and the one worth reporting.""")))
C()   # pymrm conduction

cells.append(md(sub(r"""Second order on the grid («nu_order» observed), linear in $a/b$ on the domain, and
the refined-and-extrapolated finite-volume result sits «nu_routes» relative from the
exact 2. **That is the one number on this page checked against something other than
itself, and it is what a break table cannot give**: every row of the table below
perturbs an input and watches a number move, which shows sensitivity and never
correctness.

Note what the domain refinement is for. Stopping at $b/a = 10$ and calling the
answer 2.0 would be wrong by «b10_pct» %, and no grid refinement would reveal it — the
break table carries that row precisely because it is the failure mode a
grid-convergence study is blind to.

### 2. Radial convection, and which average $p_f$ is

Ranz, 1993:

> "A serious limitation on extended applications, unfortunately not emphasized in
> the original papers, has always been crude accounting (i.e., the use of the
> $D_f$ factor) for the radial convection at a spherical boundary caused by
> diffusion at the same boundary."

There is no $D_f$ in the notation of either part. The object that does that
accounting is $p_f$, and Part II's notation (folio 180) defines it only as the
"average value of $(\pi - p_A)$ across transfer path" — the identification is an
inference and is labelled one here. Which average it is has a right answer, and the
solver gives it.""")))
C()   # pymrm stefan

cells.append(md(sub(r"""The logarithmic mean returns Ranz's $N_{Nu}'$ to 2.0 at rest (to «pf_log», which is
the pymrm-vs-closed-form disagreement over again and not a second number — the log
mean cancels the closed form identically, so it is reported once under its own
name); the arithmetic mean does not, and is «pf_arith» % out by $y = 0.20$.
So $p_f$ **must** be the log mean, and read that way the paper's mass-transfer
accounting for radial convection is *exact*, not crude — for the isothermal
spherical stagnant film. What it does not do is appear anywhere in the heat
balance, which is where a Stefan flow also carries enthalpy; that asymmetry is the
most likely reading of the 1993 remark, and it is offered as a reading, not as a
result.

The size of the effect: the radial convection lifts $Sh$ by 1 % at a surface mole
fraction of «y1» %, 5 % at «y5» % and 10 % at «y10» %, each root-found rather
than read off a sweep. Placing this paper's runs on that axis needs vapour
pressures that exist only in Figure 14, so it is not done.

### 3. The paper's own arithmetic""")))
C()   # printed arithmetic

cells.append(md(sub(r"""No printed defect found. Test No. 80 recomputes to «test80» relative, i.e. the
printed 5.57 is correct to its last digit; and eq. (24)'s 0.54 is within
«eq24_pct» % of what Table 1's own 18 rows give with the intercept held at 2.0.
Both were checked because both *could* have failed — this section reports that
they did not, rather than omitting a check that came out clean.

That second number is **not a second result**. $K_{\rm recovered} \equiv b_{\rm fixed}
\times 0.60/0.54$ identically, so "eq. (24)'s 0.54 against Table 1" and the heat
section's "$K = $ «heat_K» against the printed 0.60" are one statement written two
ways. It is printed in both places because a reader arriving at either section should
not have to reconstruct it, and it is reported **once** in `agreement.json` — the same
rule that kept `pf_logmean_returns_2_rel` out of the metrics below.

### 4. Defect injection

«n_break» rows against «n_metrics» metrics. Every one of the «n_metrics» has at
least one row that moves it **past `check_agreement.py`'s own 5 % tolerance** — a
row that moves a number by less than CI would notice is not coverage, and the
assertion in the next cell enforces that distinction («n_strong» of «n_break» rows
clear it). No metric on this page falls below `ABS_FLOOR = 1e-12`, so all «n_metrics»
are inside the regression suite. The table is rebuilt for this page's physics;
nothing in it travelled from another page.""")))
C()   # break table

cells.append(md(sub(r"""### 5. Coverage, asserted key for key

Every metric in `agreement.json` must appear in the break table above, and the
assertion fails the build if one does not.""")))
C()   # coverage

cells.append(md(sub(r"""**What the break table cannot do**, stated because it is the defect this
repository keeps finding. Every row perturbs an input and watches a number move.
That establishes sensitivity. It cannot catch a baseline that is wrong by
accident — a value quantised by round-off, limited by a grid, read at the wrong
place, or right for the wrong reason. Four defences are used here instead:
the finite-volume 2.0 is measured **against the exact answer**, not against a
second estimate («nu_routes»); every threshold and every limit is
**root-found or extrapolated**, never sampled (the Stefan thresholds by `brentq`,
the conduction limit by Richardson in $h$ and then in $a/b$); the datasets
themselves are checked against **over-determinations the papers cannot avoid** (the
Reynolds column against the columns beside it, «re_rms» % rms over «re_n» rows; and
Table 4's $N_{Nu}'/N_{Nu}$ ratio against its own temperature columns); and the
recovered reference temperature is checked **out of sample**, on Table 3's nine rows
at «t3ho_Tlo»–«t3ho_Thi» K, and against the one place the papers state the
convention themselves.

**A worked example of the failure mode, from this page's own history.** Table 4 run
11's degraded digit was first read 10.6. Every metric the old page reported moved by
less than «benz_drift» % — under CI's tolerance — so the break table, the coverage
assertion and `check_agreement.py` all passed on a wrong number. What catches it is
not a perturbation but a **constraint**: the row's $N_{Nu}'/N_{Nu}$ ratio must match
its neighbours', and it did not, by «t4_spread_bad» % against «t4_spread» %. That
ratio is now a reported metric with its own break row. A break table tells you which
of your numbers are load-bearing; only a constraint tells you whether one is wrong.

«n_weak» of the «n_break» rows do not clear the 5 % tolerance. **One** of them is
kept precisely because it barely moves and says so in its own label — reading Table 4
run 11's digit the discarded way. The other «n_weak_other» are ordinary rows whose
metric happens to be insensitive to that particular perturbation and which are
covered elsewhere; they are not, and are not claimed to be, deliberate
"barely moves" exhibits.

## What pymrm adds

**Not the correlation.** Eqs. (21) and (22) are algebra and need no solver; this
page says so rather than dressing an evaluation up as a simulation. Four things
are added, and two of them need no solver at all.

1. **The `2.0` becomes a computation instead of a citation.** Part I asserts eq.
   (7); pymrm solves the spherical BVP it comes from, on a refined grid and to an
   extrapolated domain, and gets «nu_inf» — «nu_routes» relative from the exact
   answer, which the same cell obtains by quadrature. That matters practically, not
   decoratively: the same solve says a
   finite far field of ten radii shifts the "2" by «b10_pct» %, which is the correction
   anyone applying Ranz-Marshall inside a dense spray or a packed bed needs and
   which the correlation itself cannot express. (`A3.4`, Wakao-Funazkri, is the
   page for packed beds; this one does not go there.)

2. **A number on the radial-convection factor Ranz said the papers underplayed.**
   Adding `construct_convflux_upwind` to the same spherical solve turns "which
   average is $p_f$?" — a genuine ambiguity in Part II's own notation — into a
   settled question: the logarithmic mean, to «pf_log»; the arithmetic mean is
   «pf_arith» % out at $y = 0.20$. And the correction is quantified where the
   papers only gesture at it: 1 % at $y = $ «y1» %, 10 % at «y10» %.

3. **The fit/test split, which is the actual result and needed no solver.**
   Separating the theoretical 2.0 from the fitted 0.60 and the assumed exponents,
   and then extrapolating three independent sets of runs to $N_{Re}=0$, is what
   turns "the correlation fits the data to a few per cent" into "two of the three
   land together at about «heat_dev» % and the third misses by «mass_dev» %, in the
   direction the authors' own first-listed error term predicts". The authors state
   the direction — "Data for mass transfer show a steeper slope and a lower
   intercept" — and name the cause — "(1) inaccurate values of diffusivity". What is
   new here is the **size**, per species, against their own replicate scatter, and
   the resulting tension with the two conclusions they *did* draw from it: the
   abstract's "verified the simple expression for the Nusselt number at zero Reynolds
   number", and folio 174's implication that their $D_v$ "may be more accurate than
   any reported in the literature".

4. **A convention recovered from the tables, and then tested out of sample.** The
   Reynolds column is over-determined by the columns printed beside it; asking it at
   what temperature the air properties were evaluated returns $w = $ «re_w», the film
   temperature. Two things then corroborate it independently: Part II's own statement
   of a reference temperature for the still-air run («T_stated» K against the row's
   film temperature of «still_film_K» K), and Table 3 — nine rows with no Reynolds
   column, which therefore cannot have entered the fit — predicted to «t3ho_rms» %
   against «t3ho_fs» % at the free stream, «t3ho_Tlo»–«t3ho_Thi» K.

## Reuse

**Use `Nu = 2 + 0.6 Re^(1/2) Pr^(1/3)` for an isolated sphere or drop in a moving
fluid at 0 < Re < 200.** That is what the paper's own abstract claims ("a Reynolds
number range of 0 to 200"), and it covers Tables 1–3; note that Table 4's last
benzene run is printed at $N_{Re} = 220$, so the tables reach slightly past the
stated range. **Part I's abstract** — not Part II — says the results "could be
extrapolated with remarkable accuracy five times beyond the experimental range of
Reynolds numbers"; Part II puts it as "remarkably accurate even when extrapolated
five times beyond the experimental range in which they were determined" (folio 174).
Either way that claim rests on Figure 12, which is not transcribed here, so this
page neither supports nor contradicts it.

**Five warnings, all of them measured above.**

- **Evaluate the air properties at the film temperature, the mean of the air and
  drop temperatures.** Neither part states this for the tables; it is recovered on
  this page from the printed Reynolds column, at $w = $ «re_w» with the fitted
  temperature exponent landing on kinetic theory's $-1.7$, and it is corroborated
  both by the one reference temperature Part II does print («T_stated» K against a
  film temperature of «still_film_K» K) and out of sample on Table 3. Using the
  free-stream temperature instead moves $N_{Re}$ by up to «re_fs_re» % on the paper's
  own rows and the predicted $N_{Nu}$ by up to «re_fs_nu» %, and the error grows with
  the air-to-drop temperature difference, so it is worst exactly where drying
  calculations live.
- **The Prandtl exponent is not testable from anything transcribed here — but the
  paper does test it, graphically.** Every run in Tables 1–4 is in air, so the $1/3$
  cannot be identified from these rows. Figure 12 carries Kramers's steel-to-oil and
  steel-to-water points on the $N_{Pr}^{1/3}N_{Re}^{1/2}$ abscissa and folio 174
  credits the exponents on that basis. That figure is not digitised here, so if your
  fluid is not air the $1/3$ reaches you supported by a plot this page did not check,
  not by nothing — and the coefficient and the exponent are still not separable from
  these tables.
- **Expect roughly 10 % at high air temperature.** Table 3 (85–221 °C) sits
  «t3_rms» % rms from the correlation *while inside its own fit set*, against
  «heat_rms» % for the room-temperature runs.
- **The heat and mass forms are not interchangeable at low Reynolds number.** The
  paper's own still-air pair differs by «still_ratio_pct» % where the two
  equations demand identity.
- **If you use the mass-transfer form, watch which diffusivity you pair it with.**
  $N_{Nu}'$ was reduced from rate data using the authors' own $D_v$; adopting a
  different one rescales the whole correlation, and the «mass_dev» % intercept
  deficit found here is exactly that sensitivity showing up in their water data.
  The authors list "inaccurate values of diffusivity" first among the causes of the
  anomaly themselves, and separately conclude that their $D_v$ may be the most
  accurate available; those two statements pull in opposite directions and this page
  cannot settle which is right.

**Do not use this page for**: free convection (the $N_{Gr}^{1/4}$ form of eqs. 10
and 11 is never exercised here — say so if you cite it); drops containing
dissolved or suspended solids (Part II's second half, whose calculated rates need
Figure 14); or packed beds, which is `A3.4`.

**Datasets are reusable on their own.** All five carry `columns:` blocks and
provenance sidecars, and all five are table transcriptions, so nothing in them
depends on a figure reading or on this page's argument.""")))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                             "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
Path("index.ipynb").write_text(nbf.writes(nb), encoding="utf-8")
print(f"wrote index.ipynb: {len(cells)} cells ({sum(c.cell_type == 'code' for c in cells)} code)")
print(f"quoted {len(V)} interpolated values; {len(M)} metrics; {len(NS['BT'])} break rows")
