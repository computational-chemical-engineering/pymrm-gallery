# J4.8: ASM1 activated sludge model

From Henze, M., Grady, C. P. L. Jr, Gujer, W., Marais, G. v. R. & Matsuo, T.,
*A general model for single-sludge wastewater treatment systems*, Water Research
**21**(5), 505-515 (1987), doi:10.1016/0043-1354(87)90058-3. This is the
abbreviated report of IAWPRC Scientific and Technical Report No. 1.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page source, paired with `index.ipynb` (`jupytext --sync`); edit this |
| `index.ipynb` | the executed page |
| `model.py` | the symbolic matrix, the composition vectors, the rates and the lane (`Lane`) |
| `test_model.py` | fast offline checks of `model.py` |
| `agreement.json` | the page's claims, written by its last cell and checked by CI |
| `meta.yaml` | page metadata, validated against `models.yaml` |
| `data/henze-1987-table1-components.csv` | the components, names and units |
| `data/henze-1987-table1-stoichiometry.csv` | the Petersen matrix, as expressions |
| `data/henze-1987-table1-rates.csv` | the process rate expressions |
| `data/henze-1987-table4-parameters.csv` | typical values at 20 °C |

Each CSV has a `.meta.yaml` provenance sidecar.

## Why this model

The model **is** a table. Reproducing it is the page, and the table is
**self-checking**: the paper declares its components in COD, in nitrogen and in
molar charge, and those declarations force three conservation relations on every
row. Either they close or they do not, and there is no tolerance to negotiate.

## What the checks establish

- **Charge closes on every process, exactly**, as a symbolic identity in
  Y_H, Y_A, f_P, i_XB, i_XP, not at a sample parameter set.
- **Nitrogen closes on every process**, but only after adding a fourteenth
  component S_N2. Without it, process 2 leaves `(1-Y_H)/(2.86*Y_H)` g N per g of
  biomass grown anoxically, 0.16 <!-- claim: n_loss_no_n2_min --> to
  0.41 <!-- claim: n_loss_no_n2_max --> over Table 4's yield range, with nowhere
  to go. ASM2 and ASM3 add S_N2 for exactly this reason.
- **COD closes on 6 <!-- claim: cod_closed_count --> processes, all but two.**
  The two residuals are `(Y_H - 1)/(1001*Y_H)` <!-- claim: cod_residual_p2_symbolic -->
  and `-1/(700*Y_A)` <!-- claim: cod_residual_p3_symbolic -->, and those are
  *exactly* the rounding of 20/7 to 2.86 and of 32/7 to 4.57. Put the exact
  fractions in and every balance is identically zero.
- **Printed coefficients are recovered by solving continuity for them.**
  Of the ten tried, 8 <!-- claim: coefficients_recovered_identically --> come back
  symbolically identical; the other 2 <!-- claim: coefficients_recovered_up_to_rounding -->
  differ by the same rounding. A single mis-read digit anywhere in the oxygen,
  nitrate, ammonium, alkalinity or X_ND columns would break this.
- **The paper's equation (3) is a free cross-check.** It writes out
  r₂ = ν₂₁ρ₁ + ν₂₂ρ₂ + ν₂₇ρ₇ in full, typeset separately from the matrix.
  Rebuilding it from the transcribed matrix column and the transcribed rate rows
  reproduces it exactly.

## The key insight, in one line

The paper's two conversion factors are two differences on one linear scale.
Reference nitrogen at -III (ammonium and organic amine, which is where dichromate
COD leaves it), then per gram of N the theoretical oxygen demand is 0 for NH₄⁺,
-12/7 for N₂ and -32/7 for NO₃⁻. So

- 32/7 - 12/7 = 20/7 = 2.857..., the printed **2.86**, nitrate to N₂
- 32/7 - 0 = 32/7 = 4.571..., the printed **4.57**, ammonium to nitrate

That is why a *single* COD continuity can cover both denitrification and
nitrification, and it is the point the whole page turns on.

## The reactor, and what it is for

A plug-flow lane with an unaerated head end and internal mixed-liquor recycle
(the Modified Ludzack-Ettinger arrangement), because the paper says its
hydrolysis model exists to produce "realistic space-time ... electron acceptor
profiles". `construct_convflux_upwind` + `construct_div` for transport with an
`{"outflow": True}` exit, `NumJac((n, 14))` + `newton` (bounded below by zero)
for the steady state, and the recycle inside the same Newton system: the inlet
value enters through `shapes_d` and the outlet value is the exit-face row of the
convection operator, so one solve replaces the fixed-point iteration on the
inlet. That iteration is kept as a second route: it needs
57 <!-- claim: picard_iterations --> lane solves, its updates shrinking by
0.65 <!-- claim: picard_contraction --> per iterate, close to R/(1+R) for a
passive component, and agrees with the one-solve route to below
1e-8 <!-- claim: recycle_routes_gap --> g m⁻³.

Its job is to carry the same check through the
discretisation: at steady state the
flux of a conserved quantity through the lane can only change by what crosses
the boundary. Nitrogen and charge close to below
1e-8 <!-- claim: reactor_nitrogen_balance --> and below
1e-8 <!-- claim: reactor_charge_balance --> g m⁻³ d⁻¹ at every grid tried,
against an oxygen demand of 652 <!-- claim: lane_otr --> g m⁻³ d⁻¹; COD leaks
0.2071 <!-- claim: reactor_cod_leak --> g m⁻³ d⁻¹, which the two stoichiometric
residuals predict to below 1e-8 <!-- claim: reactor_cod_leak_vs_prediction -->
relative, and which falls to below 1e-8 <!-- claim: reactor_cod_leak_exact_factors -->
g m⁻³ d⁻¹ with the exact fractions. For a conservative scheme at a converged
solve this is an identity: the balances test the operators, the solve and the
numeric matrix, not the transcription a second time.

Grid refinement on n = 100 to 1600: effluent ammonia and soluble organic
nitrogen converge at first order, observed orders from successive differences
at least 0.9 <!-- claim: lane_order_min --> and at most
1.1 <!-- claim: lane_order_max -->; the other outlet values are not yet in the
asymptotic range on that ladder (orders below 0.8 <!-- claim: lane_order_preasymptotic_max -->),
and total inorganic nitrogen on the sweeps' 100-cell grid differs from the
1600-cell value by 0.6 % <!-- claim: tn_grid_change -->.

## What the alkalinity column is for

Charge continuity alone says denitrification returns exactly half the alkalinity
nitrification destroys, a ratio of 2 <!-- claim: alk_nit_den_ratio --> with no
parameter in it. On the lane (an operating point of ours: the paper gives none,
100 cells per lane), total inorganic nitrogen is lowest at an unaerated fraction
of 0.50 <!-- claim: tn_optimum_fraction -->, 8.73 <!-- claim: tn_at_optimum -->
against 26.85 <!-- claim: tn_no_anoxic --> g N m⁻³ with no anoxic zone, a
67 % <!-- claim: tn_reduction_at_optimum --> reduction. The raw-water alkalinity
needed to stay above the paper's 50 g CaCO₃ m⁻³ limit is
3.37 <!-- claim: alk_need_f0 --> mol m⁻³ with no anoxic zone,
1.48 <!-- claim: alk_need_f30 --> at 30 % unaerated and
0.26 <!-- claim: alk_need_f50 --> at 50 %, reaching zero at an unaerated
fraction of 0.566 <!-- claim: alk_need_zero_fraction --> only because
nitrification is by then incomplete.

## Not decided by inference

**Table 4 prints one parameter row more** than the nineteen the paper says the
model has. The extra is a second ammonification coefficient: `k_a = 0.016`
m³(g cell COD·day)⁻¹ in the right place in the kinetic block, and `k_A = 0.08`
m³ COD (g·day)⁻¹ appended after `b_A` under a symbol defined nowhere. Both
subscripts were re-read at 8× magnification. Both rows are in the CSV; the
reactor is run at both. The larger value moves soluble organic nitrogen by
-79 % <!-- claim: ka_snd_rel_change --> and raises effluent inorganic nitrogen by
21 % <!-- claim: ka_tn_rel_change -->, while soluble nitrogen in total moves by
-0.5 % <!-- claim: ka_soluble_n_rel_change -->: the nitrogen moves from the
organic pool into nitrate. No kinetic parameter enters any continuity check, so
nothing in the validation depends on the choice.

## Traps met while building this

- **A recycle as a fixed point on the inlet value contracts at about R/(1+R)**
  per iterate; placing the recycle inside the Newton system removes the
  iteration.
- **The rate expressions clip negative concentrations**, so Newton is
  bounded below by zero, `bounds=(0, None)`.
- **The PDF text layer is unusable** for this paper's tables: 1/Y_H comes out as
  `I / r~` and (4.57-Y_A)/Y_A as `4 . 5 7 - Y4 / Y,~`. Journal pages 506, 507 and
  513 were rendered at 600 dpi and read as images.

## No experimental data, deliberately

The abbreviated report contains no worked example, no calibrated plant and no
dataset. Its one data figure is reproduced *from* Ekama et al. (1986) to
illustrate a measurement procedure. Re-simulating the model and calling the
result data would be circular, so the page says there is none and lets the
continuity checks carry the weight, which is the stronger test anyway, because
they admit no free parameter.

## Corrections

- 2026-11-19. The page said of the two ammonification values "The choice moves
  soluble organic nitrogen appreciably and effluent" ... "inorganic nitrogen very
  little". Effluent inorganic nitrogen rises by
  21 % <!-- claim: ka_tn_rel_change --> at the larger value; what changes
  little is soluble nitrogen in total
  (-0.5 % <!-- claim: ka_soluble_n_rel_change -->).
- 2026-11-19. The page said "using `c[-1]` puts a spurious 1e-3 into an
  otherwise machine-precision balance." Under its zero-gradient exit, nitrogen
  and charge read from the last cell still close to below
  1e-8 <!-- claim: zg_last_cell_balance --> g m⁻³ d⁻¹, since every face carries
  the same flux of a conserved quantity; only COD, which has a source along the
  lane, moves, by 0.51 <!-- claim: zg_last_cell_cod_shift --> g m⁻³ d⁻¹ with the
  exact fractions.
- 2026-11-19. The page said the COD leak "is grid independent, which"
  ... "numerical error would not be". It moves by
  0.19 % <!-- claim: reactor_cod_leak_grid_change --> between 75 and 300 cells:
  a discretisation-dependent value, small against the leak itself.
