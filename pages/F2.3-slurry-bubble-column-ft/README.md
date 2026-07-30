# F2.3 — Slurry bubble column for Fischer–Tropsch synthesis

From Maretto, C. & Krishna, R., *Modelling of a bubble column slurry reactor for
Fischer–Tropsch synthesis*, *Catalysis Today* **52**(2–3) 279–289 (1999),
[doi:10.1016/S0920-5861(99)00082-6](https://doi.org/10.1016/S0920-5861(99)00082-6).

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 2.5 min |
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
