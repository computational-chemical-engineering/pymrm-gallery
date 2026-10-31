# A2.1 — The Danckwerts boundary conditions (covers A2.2)

A tubular reactor with axial dispersion needs two boundary conditions, and the
obvious pair is wrong. Danckwerts (1953) gave the right pair: a **flux** balance
at the inlet, which makes the concentration drop across the bed entrance
(discontinuously in the closed-vessel idealisation), and a **zero gradient** at
the outlet. Wehner & Wilhelm (1956) derived both by putting the bed between an
inlet pipe and an outlet pipe, and in their fuller model the inlet drop is
continuous: they say of Danckwerts' discontinuity that it is *"not necessarily
correct"*, and the page quotes and reconciles that.

- **Structure:** `S4` (1-D convection, dispersion and reaction)
- **References:** Danckwerts (1953), *Chem. Eng. Sci.* 2(1) 1-13,
  doi:10.1016/0009-2509(53)80001-1; Wehner & Wilhelm (1956), *Chem. Eng. Sci.*
  6(2) 89-93, doi:10.1016/0009-2509(56)80014-6
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`ClosedVessel`, the closed vessel with the two dicts, and
  `ThreeSection`, Wehner & Wilhelm's fore, reaction and after sections as one
  domain), tested by `test_model.py`. The closed forms stay in the notebook.
- **Runtime:** about 7 s.
- **Data:** none. Tier 6, analytical; no figure on either paper is digitised.

## Why one page for two catalog entries

Wehner & Wilhelm's eq. (20) for the reaction section is algebraically identical
to Danckwerts' eq. (33), and they say so on their page 91. A separate `A2.2`
page would have to validate against the same closed form. What they add is the
*justification* of Danckwerts' conditions (Danckwerts himself wrote "Intuition
suggests that neither of these situations can arise"), and a boundary condition
and its justification are one result. The full argument is in
`queue_cases/A2.1.yaml` under `scope_decision`.

## The conditions, in pymrm's convention

`bc` is `a·∂c/∂n + b·c = d` with `n` the **outward** normal, so the same dict
means different physics at the two ends. `pymrm.describe_bc(bc, z_f)` prints the
equation each dict imposes.

```python
D, u = 1.0 / Pe, 1.0
bc = (
    # inlet, Danckwerts eq. 31:  u c* = u c - D dc/dz  at z = 0.
    # n = -z, so dc/dn = -dc/dz  ->  (D/u) dc/dn + c = c*
    {"a": D / u, "b": 1.0, "d": c_star},
    # outlet, Danckwerts eq. 32:  dc/dz = 0 at z = L.
    # n = +z, so dc/dn = dc/dz   ->  a = 1, b = 0, d = 0
    {"a": 1.0, "b": 0.0, "d": 0.0},
)
```

## The numerics

Van Leer TVD convection with `ConvFluxTVD` and `newton(..., damping="monotone")`
(pitfall P10). The problem is linear apart from the limiter, so Newton starts
from the first-order upwind solution with the residual written as the correction
to it; written in full, the residual applies the diffusion entries to order-one
fractions, and Newton stops at a round-off floor that moves with the BLAS kernel. The same equations by deferred correction are the second
solver route. The exit is read on the outlet face with `compute_boundary_values`
(pitfall P8).

## Agreement

| Check | Result |
|---|---|
| the transcriptions of two papers (D. eqs. 33/34, W&W eqs. 20/21/22) | agree to round-off; one wrong digit lifts the gap to 2.9e-2 <!-- claim: transcription_one_digit_misread --> |
| pymrm vs Danckwerts eq. 34, van Leer, outlet face | **3.9e-8** <!-- claim: closed_form_rel_err_n1600 --> at n = 1600, observed order 1.97, approaching 2 from below |
| the same, read at the last cell centre | 2.2e-7 <!-- claim: closed_form_rel_err_n1600_last_cell --> at n = 1600, larger on every rung |
| bare upwind | order 1.00; the residual against `Pe_eff = 1/(1/Pe + dz/2)` is order 2.00 |
| Peclet sweep, 1e-4 to 1e4, vs eq. 34 | worst **3.6e-7** <!-- claim: pe_sweep_worst_rel_err --> on the face; 1.9e-4 <!-- claim: pe_sweep_worst_last_cell --> at the last cell centre |
| eq. 35 plug-flow limit / eq. 36 stirred-tank limit (closed form) | 4.5e-8 <!-- claim: plug_flow_limit_eq35 --> / 2.2e-9 <!-- claim: stirred_tank_limit_eq36 --> |
| eq. 37 asymptote, residual decay exponent | **2.00** over Pe ≥ 1e3, so `R²/Pe` is the exact first correction (1.95 <!-- claim: eq37_residual_order_full_range --> if the pre-asymptotic end is included) |
| W&W three sections, **no Danckwerts BC imposed anywhere** | reproduces the closed vessel to **7.8e-6** <!-- claim: three_section_vs_eq33 -->, with Pe_a and Pe_c swept |
| harmonic vs arithmetic face average of `D` at the jump | an **order**, not a factor: 1.00 against 2.00, so the ratio is 84 <!-- claim: face_average_ratio_n100 --> / 172 <!-- claim: face_average_ratio_n200 --> / 347 <!-- claim: face_average_ratio_n400 --> / about 700 at `n_b` = 100 / 200 / 400 / 800 |
| naive Dirichlet inlet vs its own closed form (Hulburt = eq.34 / eq.22) | 9.4e-7 <!-- claim: dirichlet_vs_hulburt_closed_form --> over four decades |

The `n_b = 800` ratio sits on a 32800-cell uniform grid, where round-off alone
(two algebraically identical assemblies of the same operator, or the solver
route) moves the fore-section error by up to about 0.2 %; the last digit of that
one entry is not resolved.

## The defect table

| injected defect | rel. deviation from eq. 34 |
|---|---|
| none | 1.5e-7 <!-- claim: closed_form_rel_err_n800 --> |
| inlet BC: sign of `a` flipped | **1.97** <!-- claim: defect_inlet_sign_flip --> |
| inlet BC: naive Dirichlet `c(0) = c*` | **0.497** <!-- claim: defect_naive_dirichlet_inlet --> |
| inlet BC: `a = D` with D taken as 1 | 0.356 <!-- claim: defect_inlet_a_one --> |
| **outlet BC: sign of `a` flipped** | **1.5e-7** <!-- claim: defect_outlet_sign_flip -->, bit-identical matrix and RHS |
| outlet BC: Dirichlet `c(L) = 0` | 1.00 |
| `construct_div nu = 1` | 1.00 |
| reaction sign flipped | 26.1 <!-- claim: defect_reaction_sign --> |
| n = 10 cells | 6.3e-5 <!-- claim: defect_coarse_n10 --> (but 3.2e-3 <!-- claim: undergrid_pe200_vanleer --> at Pe = 200) |

The outlet row is the point: because `b = 0` and `d = 0` there, the sign of `a`
is immaterial, so the outward-normal error is **undetectable at that end by any
test**. Not "the number does not move": the page assembles the system for
`a = +1`, `-1` and `1e7` and finds max |ΔA| = max |Δb| = **0.0**, against 1.8e3
<!-- claim: inlet_sign_system_change --> for a flipped inlet sign. Set `d = 0.05`
and the outlet sign matters again (`f(exit)` 0.2470 <!-- claim: outlet_flux_exit_a_plus -->
against 0.2219 <!-- claim: outlet_flux_exit_a_minus -->), so the because-clause is
doing the work. Verifying an outlet condition proves nothing about an inlet
condition.

## What the page adds

The size of the mistake. `c(0) = c_in` over-predicts the outlet concentration by
179 % <!-- claim: naive_f_error_pe0p1 --> at Pe = 0.1, 50 % <!-- claim: dirichlet_inlet_error_at_ww_fig4 -->
at the Peclet number Wehner and Wilhelm printed on their own Figure 4, and 1 %
<!-- claim: naive_f_error_pe200 --> at Pe = 200, and as Pe → 0 it predicts *no
conversion at all* where the answer must be a stirred tank. Both papers had that
algebra; neither plotted it.

Also: bare upwind adds a numerical diffusivity `u·Δz/2` that is
indistinguishable from a lower Peclet number, so a 200-cell grid at Pe = 200 is
really running at Pe_eff = 133 <!-- claim: upwind_pe_eff_at_pe200_n200 -->. A
boundary-condition study on an upwind grid can measure its own truncation error
and call it physics. And at high Peclet number the outlet boundary layer is
thinner than a cell: at Pe = 1e4 the solve's exit, on the outlet face, agrees with
eq. 34 to 6.2e-8 <!-- claim: pe_sweep_rel_err_at_pe1e4 -->, while the last cell
value lies 1.9e-4 <!-- claim: pe_sweep_last_cell_at_pe1e4 --> from it.

## Honest limits

Tier 6, not experimental. Neither paper measures anything that tests these
boundary conditions; Danckwerts' one experimental figure belongs to his
*open*-vessel residence-time analysis, which he says himself cannot be carried
into the reactor problem. The demonstration that `df/dz = 0` emerges at the bed
exit reproduces Wehner & Wilhelm's chain of reasoning, not an independent proof
of its premise: imposing boundedness at the far downstream truncation is itself
a zero-gradient condition, applied one section away. Nor is the three-section
route independent *code*: it shares the whole operator stack and the solver with
the closed-vessel solve. Its condition at `z = 1` shares the zero diffusive flux
of `bc_out = {a:1, b:0, d:0}` but not the convective flux: the closed vessel's
outlet face carries `(9 c_N - c_(N-1))/8`, the three-section bed exit carries
`c_N`, and the two differ by 1.0e-6 <!-- claim: closed_outlet_face_vs_last_cell -->
relative at n = 800.

Four identities of the discrete equations are claimed as bounds: the after
section being flat, at its full length and cut short (cutting it from 20 decay
lengths to 0.2 leaves it just as flat); the bed being independent of that
length; the last bed cell equalling the first after-section cell; and the outlet
sign being absent from the assembled system. None is structural: letting the
reaction run on past `z = 1` breaks the first three, and an outlet `d = 0.05`
with `a` flipped breaks the last. The fore-section flux is
constant by construction; only its value informs, as a domain-truncation check.
The recovered `df/dz(1)` is 84 % <!-- claim: grad_at_exit_estimator_share --> the
truncation error of the 3-point fit that recovers it, measured by applying that
fit to the exact eq. 33.

The residence-time half of Danckwerts 1953 (F- and C-diagrams, hold-back,
segregation, the laminar-pipe distribution) is not built here. `A2.4` is the
catalog entry for it.

## Reproduce

```bash
cd pages/A2.1-danckwerts-boundary-conditions
jupytext --sync index.py      # index.py is the source; index.ipynb is paired with it
jupyter nbconvert --to notebook --execute --inplace index.ipynb
GALLERY_OFFLINE=1 python -m pytest -q test_model.py
```

No PDF and no dataset needed.

## Corrections

- 2026-10-31: stated that the pymrm solve agrees with eq. 34 to 2.2e-7 at
  n = 1600, to 1.9e-4 at worst over the Peclet sweep "at Pe = 1e4 where the cell
  Peclet number reaches 6", that a 10-cell grid costs 0.6 % (4.1 % at Pe = 200),
  that the Dirichlet-inlet solve matches Hulburt's closed form to 7.1e-5, and that
  the unbroken model deviates by 8.9e-7 at n = 800. Each read the exit at the last
  cell centre, half a cell inside the bed (pitfall P8). The solve's exit, on the
  outlet face, agrees to 3.9e-8 <!-- claim: closed_form_rel_err_n1600 -->, worst
  3.6e-7 <!-- claim: pe_sweep_worst_rel_err --> over the sweep and
  6.2e-8 <!-- claim: pe_sweep_rel_err_at_pe1e4 --> at Pe = 1e4; the 1.9e-4
  <!-- claim: pe_sweep_last_cell_at_pe1e4 --> there is the last cell's distance
  from the exit, which is large because the outlet layer is thinner than a cell.
  The 10-cell grid costs 0.0063 % <!-- claim: defect_coarse_n10 --> (0.32 %
  <!-- claim: undergrid_pe200_vanleer --> at Pe = 200); Hulburt is matched to
  9.4e-7 <!-- claim: dirichlet_vs_hulburt_closed_form -->, and the unbroken model
  deviates by 1.5e-7 <!-- claim: closed_form_rel_err_n800 -->.
- 2026-10-31: stated, as a trap found while building the page, that a convection
  outlet left as `None` makes the matrix singular while the solve still returns a
  plausible profile. The assembled matrix is regular, and the profile rises
  towards the outlet instead of flattening.
- 2026-10-31: stated that the three-section problem's effective condition at
  `z = 1` is discretely identical to the closed vessel's `bc_out = {a:1, b:0, d:0}`.
  Only the zero diffusive flux is identical. The closed vessel's outlet face
  carries `(9 c_N - c_(N-1))/8` and the three-section bed exit carries `c_N`; at
  n = 800 they differ by 1.0e-6 <!-- claim: closed_outlet_face_vs_last_cell -->
  relative.
