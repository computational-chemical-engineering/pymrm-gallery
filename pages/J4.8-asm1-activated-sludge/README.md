# J4.8 — ASM1 activated sludge model

From Henze, M., Grady, C. P. L. Jr, Gujer, W., Marais, G. v. R. & Matsuo, T.,
*A general model for single-sludge wastewater treatment systems*, Water Research
**21**(5), 505–515 (1987), doi:10.1016/0043-1354(87)90058-3. This is the
abbreviated report of IAWPRC Scientific and Technical Report No. 1.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 58 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/henze-1987-table1-components.csv` | the 13 components, names and units |
| `data/henze-1987-table1-stoichiometry.csv` | the 8 × 13 Petersen matrix, as expressions |
| `data/henze-1987-table1-rates.csv` | the 8 process rate expressions |
| `data/henze-1987-table4-parameters.csv` | typical values at 20 °C |

Each CSV has a `.meta.yaml` provenance sidecar.

## Why this model

The model **is** a table. Reproducing it is the page, and the table is
**self-checking**: the paper declares its components in COD, in nitrogen and in
molar charge, and those declarations force three conservation relations on every
row. Either they close or they do not, and there is no tolerance to negotiate.

## What the checks establish

- **Charge closes on all 8 processes, exactly**, as a symbolic identity in
  Y_H, Y_A, f_P, i_XB, i_XP — not at a sample parameter set.
- **Nitrogen closes on all 8**, but only after adding a fourteenth component
  S_N2. Without it, process 2 leaves `(1-Y_H)/(2.86*Y_H)` g N per g of biomass
  grown anoxically — 0.16 to 0.41 over Table 4's yield range — with nowhere to
  go. ASM2 and ASM3 add S_N2 for exactly this reason.
- **COD closes on 6 of 8.** The two residuals are `(Y_H-1)/(1001*Y_H)` and
  `-1/(700*Y_A)`, and those are *exactly* the rounding of 20/7 to 2.86 and of
  32/7 to 4.57. Put the exact fractions in and all 24 balances are identically
  zero.
- **Ten printed coefficients are recovered by solving continuity for them.**
  Eight come back symbolically identical; the other two differ by the same
  rounding. A single mis-read digit anywhere in the oxygen, nitrate, ammonium,
  alkalinity or X_ND columns would break this.
- **The paper's equation (3) is a free cross-check.** It writes out
  r₂ = ν₂₁ρ₁ + ν₂₂ρ₂ + ν₂₇ρ₇ in full, typeset separately from the matrix.
  Rebuilding it from the transcribed matrix column and the transcribed rate rows
  reproduces it exactly.

## The key insight, in one line

The paper's two conversion factors are two differences on one linear scale.
Reference nitrogen at −III (ammonium and organic amine, which is where dichromate
COD leaves it), then per gram of N the theoretical oxygen demand is 0 for NH₄⁺,
−12/7 for N₂ and −32/7 for NO₃⁻. So

- 32/7 − 12/7 = 20/7 = 2.857… — the printed **2.86**, nitrate to N₂
- 32/7 − 0 = 32/7 = 4.571… — the printed **4.57**, ammonium to nitrate

That is why a *single* COD continuity can cover both denitrification and
nitrification, and it is the point the whole page turns on.

## The reactor, and what it is for

A plug-flow lane with an unaerated head end and internal mixed-liquor recycle
(the Modified Ludzack–Ettinger arrangement), because the paper says its
hydrolysis model exists to produce "realistic space-time … electron acceptor
profiles". `construct_convflux_upwind` + `construct_div` for transport,
`NumJac(shape)` + `newton` for the steady state.

Its job is a second, independent form of the same check: at steady state the
flux of a conserved quantity through the lane can only change by what crosses
the boundary. Nitrogen and charge close to 1e-12 and 3e-14 g m⁻³ d⁻¹ against an
oxygen demand of 652; COD leaks 0.2072, which the two stoichiometric residuals
predict to 6e-11 relative, and which drops to 7e-12 with the exact fractions.
The same defect shows up symbolically and numerically, and both vanish together.

## Not decided by inference

**Table 4 prints twenty parameter rows** against the nineteen the paper says the
model has. The extra is a second ammonification coefficient: `k_a = 0.016`
m³(g cell COD·day)⁻¹ in the right place in the kinetic block, and `k_A = 0.08`
m³ COD (g·day)⁻¹ appended after `b_A` under a symbol defined nowhere. Both
subscripts were re-read at 8× magnification. Both rows are in the CSV; the
reactor is run at both and the difference reported. No kinetic parameter enters
any continuity check, so nothing in the validation depends on the choice.

## Traps met while building this

- **The recycle fixed point contracts at R/(1+R)** — about 57 iterations at
  R = 2. Over-relaxing by (1+R), the obvious Newton acceleration, diverges: the
  fast components (S_S, S_O, S_NH) have essentially no memory of the inlet and
  their contraction factor is nothing like R/(1+R).
- **Take the boundary fluxes from the convection operator**, not from the end
  cell values. With a zero-gradient outlet the face value is extrapolated, and
  using `c[-1]` puts a spurious 1e-3 into an otherwise machine-precision balance.
- **Put a cell face on the aeration boundary.** Letting the anoxic/aerobic step
  fall inside a cell makes the discretisation error jump around with n and the
  convergence study unreadable.
- **The PDF text layer is unusable** for this paper's tables: 1/Y_H comes out as
  `I / r~` and (4.57−Y_A)/Y_A as `4 . 5 7 - Y4 / Y,~`. Journal pages 506, 507 and
  513 were rendered at 600 dpi and read as images.

## No experimental data, deliberately

The abbreviated report contains no worked example, no calibrated plant and no
dataset. Its one data figure is reproduced *from* Ekama et al. (1986) to
illustrate a measurement procedure. Re-simulating the model and calling the
result data would be circular, so the page says there is none and lets the
continuity checks carry the weight — which is the stronger test anyway, because
they admit no free parameter.
