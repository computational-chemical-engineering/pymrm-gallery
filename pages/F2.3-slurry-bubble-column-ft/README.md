# F2.3 — Slurry bubble column for Fischer–Tropsch synthesis

From Maretto, C. & Krishna, R., *Modelling of a bubble column slurry reactor for
Fischer–Tropsch synthesis*, *Catalysis Today* **52**(2–3) 279–289 (1999),
[doi:10.1016/S0920-5861(99)00082-6](https://doi.org/10.1016/S0920-5861(99)00082-6).

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 80 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `data/maretto-krishna-1999-fig2.csv` | 99 markers digitised from Figure 2, reviewed |
| `meta.yaml` | page metadata, validated against `models.yaml` |

## Two printed constants that had to be corrected

Both are stated on the page with the alternatives shown, and both were resolved
against **the paper's own reported results**, not by guesswork.

1. **Eq. 2's prefactor** is printed as 8.8533 × 10³ mol/(s·kg_cat·bar²). That
   gives an intrinsic rate of 1.6 × 10⁴ mol/(s·kg_cat); a commercial cobalt
   catalyst runs near 10⁻². The paper's own 96 % conversion at *U* = 0.12 m/s
   requires 1.8 × 10⁻², and 8.8533 × 10⁻³ delivers 1.6 × 10⁻². **Used as 10⁻³.**
2. **Eq. 1's rate is the CO rate**, despite the subscript `CO+H₂`. As a syngas
   rate the model gives 71.5 %/24.1 % against the reported 96 %/63 %; as a CO
   rate it gives 93.1 %/63.7 %.

What makes the diagnosis safe is that the reactor is **chemically controlled** —
the paper's own kLa sensitivity result is reproduced first — so nothing else is
free to absorb the error.

## What is validated, and what is only reproduced

Keeping these apart is the point of the page's structure.

- **Hydrodynamics: validated against measurement.** 5.9 %, 6.3 % and 4.9 % mean
  deviation at ε_s = 0, 0.16, 0.35 over 79 digitised markers, one set of
  constants, nothing fitted.
- **Reactor conversions: reproduced, not validated.** The comparison is against
  the authors' simulation output. No commercial FT slurry column data is public.

## Two traps in the pymrm assembly

- The outlet **needs** a real outflow condition. Leaving it `None` makes the
  pure-convection matrix singular, and the failure is silent — a rank-deficient
  solve returns a plausible-looking profile.
- With a varying velocity, discretise d(*Uc*)/d*z* as the **divergence of the
  flux**, not *U*·d*c*/d*z*. The latter loses the gas contraction, worth 65 % of
  the volumetric flow at full conversion.

## Relationship to F1.4

Eq. 9 here has coefficient **0.3**; `F1.4`'s Eq. 19 has **0.268**. The earlier
paper fitted gas–liquid systems, this one slurries above ε_s = 0.16. The page
tests both on the same points: they are not interchangeable.

## The constant nothing on this page pins

Eq. 9 is `ε_b = C · min(D_T, 1)^-0.18 · (U − U_df)^0.58`, and the `-0.18` is
**untestable from this page in both directions.** Section 5 of the notebook
measures it rather than asserting it.

- **On the reactor path the exponent is provably inert.** The diameter
  dependence is capped at D_T = 1 m and this reactor is 7 m across, so the base
  of the power is exactly 1 and `1**x = 1` for every x. Deleting `-0.18`, or
  replacing it with `-5`, leaves ε_b, ε_total and both published conversions
  **bit-identical**. That also means the metric recording it is exactly 0.0,
  which puts it under `check_agreement.py`'s `ABS_FLOOR = 1e-12` — CI never
  compares it, and it is published as a record, not as a guard.
- **On the figure path the exponent is live but inseparable.** At the 0.10 m
  column of Figure 2 the factor is 1.5136, and perturbing the exponent alone
  moves the mean holdup deviation from 4.9–6.3 % to as much as 38 %. But all 79
  markers are at that one diameter, so what they constrain is the **product**
  `C · D_T^n = 0.454`. Rescale `C` to hold the product fixed and any exponent
  reproduces the figure to the last bit. The `C = 0.300` vs `0.268` study tests
  the product, not the exponent.
- **In between it is both live and untested.** A 0.10 m column gives ε_b **51 %**
  larger than the 7 m reactor's at the same (U, ε_s), and moves the U = 0.40 m/s
  conversion by 4.8 points. Anything at or above 1 m is indistinguishable from
  7 m. If your column is under a metre, check the exponent against a source that
  measured more than one diameter.

Two other things section 5 measures: a transposed digit in the rate prefactor
(`8.8533e-3` read as `8.5833e-3`) moves the U = 0.40 conversion by 2.2 % and the
U = 0.12 conversion by 0.2 % — **both below the 5 % relative tolerance CI
applies**, so neither conversion metric would flag it; and forcing the outlet to
a Dirichlet zero raises the operator's error against `exp(-kz/U)` from 1.4e-2 to
9.0e-2, which is the measurement behind the outflow-condition warning above.
