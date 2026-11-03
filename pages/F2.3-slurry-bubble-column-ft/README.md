# F2.3: Slurry bubble column for Fischer–Tropsch synthesis

From Maretto, C. & Krishna, R., *Modelling of a bubble column slurry reactor for
Fischer–Tropsch synthesis*, *Catalysis Today* **52**(2–3) 279–289 (1999),
[doi:10.1016/S0920-5861(99)00082-6](https://doi.org/10.1016/S0920-5861(99)00082-6).

- **Structures:** `S7` (multiphase coupling), `S8` (reactor plus transport)
- **Runtime:** about 20 s.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source (jupytext percent format), paired with `index.ipynb`. Edit this, not the notebook |
| `index.ipynb` | the executed page |
| `model.py` | `hydrodynamics` (eqs. 6 to 11), `kinetics` and `rate` (eqs. 1 to 3), `Operating` (one operating point), `SlurryColumn` (the reactor on pymrm operators, one `newton` solve with an exact Jacobian) and `reactor_ode` (the independent ODE route) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | every number the page states, with how CI checks it |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/maretto-krishna-1999-fig2.csv` | Figure 2: 115 <!-- claim: data_rows --> rows, 95 <!-- claim: data_resolved --> resolved markers and 20 <!-- claim: data_cluster --> detections in the merged cluster at the origin, reviewed, plus its provenance sidecar |

## Two printed constants that had to be corrected

Both are resolved against the paper's own reported results, and the page shows
the alternatives.

1. **Eq. 2's prefactor** is printed as 8.8533 × 10³ mol/(s·kg_cat·bar²). That
   makes the rate at the feed 1.6e4 <!-- claim: rate_feed_printed --> mol
   CO/(s·kg_cat) and the conversion 99.6 % <!-- claim: conversion_a+3_CO_U012 -->
   at U = 0.12 and 99.9 % <!-- claim: conversion_a+3_CO_U040 --> at 0.40.
   Of the four decades tried, only 10⁻³ comes near both reported conversions.
   **Used as 8.8533 × 10⁻³.**
2. **Eq. 1's rate is the CO rate**, despite the subscript `CO+H₂`: eqs. 19 and 20
   use it with coefficients 1 and 2. As a syngas rate the model gives
   72.8 % <!-- claim: conversion_a-3_syngas_U012 --> and
   24.3 % <!-- claim: conversion_a-3_syngas_U040 --> against the reported 96 % and
   63 %.

## The contraction, and which conversion eq. 21 means

Eq. 21 shrinks the large-bubble velocity with the syngas conversion χ. Read as the
molar conversion of the large-bubble gas, `(U − U_df)/(Uⁱ − U_df) = (1 + α)/(1 + α r)`
with `r` the syngas concentration relative to the feed; read as `χ = 1 − r`, it is
`1 + α(1 − r)`. The page uses the molar reading: 93.3 % <!-- claim: conversion_U012 -->
at U = 0.12 and 65.7 % <!-- claim: conversion_U040 --> at 0.40. The concentration
reading gives 93.1 % <!-- claim: conversion_concentration_reading_U012 --> at 0.12,
farther from the 96.3 % <!-- claim: paper_implied_conversion_U012 --> the paper's
productivity implies than the molar 93.3 %, and 63.6 % <!-- claim: conversion_concentration_reading_U040 -->
at 0.40, close to the 63.5 % <!-- claim: paper_implied_conversion_U040 --> the paper's
productivity implies there. But it evaluates eq. 21 at a conversion that is 16.2 % <!-- claim: concentration_reading_chi_gap_U040 -->
(points of conversion) below the one that leaves the column. Both are on the page.

## What is validated, and what is only reproduced

- **Hydrodynamics: validated against measurement.** 5.9 % <!-- claim: holdup_mad_eps_s_0 -->,
  6.3 % <!-- claim: holdup_mad_eps_s_016 --> and 4.9 % <!-- claim: holdup_mad_eps_s_035 -->
  mean deviation at ε_s = 0, 0.16, 0.35 from 79 <!-- claim: n_used_total --> markers above 2 U_df, one set of
  constants, nothing fitted.
- **Reactor conversions: reproduced, not validated.** The comparison is against
  the authors' simulation output. No commercial FT slurry column data are public.
  At U = 0.12 the model is below the conversion the paper's productivity implies by
  3.0 % <!-- claim: gap_implied_U012 --> (points of conversion) and the page does not
  identify why. The dropped liquid outflow term is not it: 0.1 m/s of liquid adds
  0.36 % <!-- claim: liquid_outflow_gain_ul01 --> (points), and that is an upper
  bound, because the conversion counts syngas that leaves dissolved in the liquid as
  converted.
- **Numerics.** The pymrm conversion converges at first order to an independent
  ODE route and is 3.9e-6 <!-- claim: pymrm_vs_ode_n1600 --> from it at 1600 cells.
  The syngas balance closes to round-off on every solve; it is an identity of the
  conservative scheme and cannot tell the two readings of eq. 21 apart.

## Three traps in the pymrm assembly

- Write the large-bubble balance for the molar flux, `q = s(state) y`, and convect
  `q` with a constant unit-velocity operator. The advective form `(U − U_df) dc/dz`
  drops `c d(U − U_df)/dz` and breaks the balance.
- The outlet needs `{"outflow": True}`, and the exit is read with
  `compute_boundary_values`. With no condition the outlet is closed to convection:
  at U = 0.12 `newton` still converges, to a conversion that looks plausible.
- Eq. 21 depends on the state, so the face velocity must be the upstream cell's;
  a cell-centred velocity array interpolated to the faces admits less than the
  feed at the inlet face.

## The constant nothing on this page pins

Eq. 9 is `ε_b = C · min(D_T, 1)^-0.18 · (U − U_df)^0.58`. At the reactor's 7 m the
cap makes the base 1, so the exponent is inert: deleting it leaves the conversion
unchanged. At the 0.10 m column of Figure 2 it is live, but every marker is at
that one diameter, so the data constrain only `C · D_T^n` =
0.454 <!-- claim: c_dt_product_0p10m -->: rescale `C` and any exponent reproduces the figure. In
between it is live and untested: eq. 9 at 0.10 m gives ε_b
1.514 <!-- claim: eps_b_ratio_0p10m_to_7m --> times the 7 m value and moves the U = 0.40
conversion to 60.6 % <!-- claim: conversion_U040_dt010 -->.

## Relationship to F1.4

Eq. 9 here has coefficient **0.3**; `F1.4`'s eq. 19 has **0.268**. The earlier
paper fitted gas–liquid systems, this one slurries above ε_s = 0.16. The page
tests both on the same markers: they are not interchangeable.

## Corrections

- 2026-11-03: stated the conversions as 93.1 % and 63.8 % (63.7 % in the text). The
  cell-centred velocity array was extrapolated to the inlet face, which admitted
  less large-bubble gas than the conversion counted as feed. With the face fixed and
  the same reading (χ = 1 − c/c_in) the model gives 93.1 % <!-- claim: conversion_concentration_reading_U012 -->
  and 63.6 % <!-- claim: conversion_concentration_reading_U040 -->. Separately, the
  page now reads χ as the molar conversion, giving 93.3 % <!-- claim: conversion_U012 -->
  and 65.7 % <!-- claim: conversion_U040 -->: a change of reading, not a fix.
- 2026-11-03: said "The likely candidate is the liquid outflow term" for the gap at
  low velocity. A liquid outflow of 0.1 m/s, as fast as the gas, adds
  0.36 % <!-- claim: liquid_outflow_gain_ul01 --> (points of conversion); 1 m/s adds
  1.7 % <!-- claim: liquid_outflow_gain_ul1 -->.
- 2026-11-03: said "The paper's own 96 % conversion at U = 0.12 m/s requires
  1.8 × 10⁻²" mol/(s·kg_cat). That matches no rate the page computes: the mean CO
  rate the 96 % needs is 5.4e-3 <!-- claim: rate_needed_U012 -->, the rate at the
  feed with the corrected prefactor is 1.6e-2 <!-- claim: rate_feed_corrected -->,
  and the model reaches 96 % only with the prefactor scaled by
  1.83 <!-- claim: rate_factor_for_96pct_U012 --> (the one root on the bracket 1 to 4).
  Also dropped: the unsourced statement that a commercial cobalt catalyst runs near
  1e-2.
- 2026-11-03: said "reproduced: the reactor is chemically controlled, as the paper
  states" for the kLa sensitivity test. A 3-fold reduction changes the productivity
  by -6.8 % <!-- claim: kla_third_change --> (doubling a: +56.9 % <!-- claim: a_double_change -->,
  paper +60 %): the ranking holds, and the reactor is mostly, not fully, kinetically
  controlled.
- 2026-11-03: the README and `meta.yaml` counted 99 digitised markers; the CSV has
  115 <!-- claim: data_rows --> rows, 95 <!-- claim: data_resolved --> resolved markers
  and 20 <!-- claim: data_cluster --> cluster detections.
- 2026-11-03: said "Leaving it `None` makes the pure-convection matrix singular, and
  the failure is silent" about a missing outlet condition. The pure-convection part is
  true; only the reactor case was wrong. In the reactor, with the exchange sink, the
  matrix is not singular; the outlet is closed to convection and the solve converges
  at U = 0.12 to a wrong conversion.
