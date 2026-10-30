# B1.1 + B1.5 — Thiele modulus and the Weisz-Hicks non-isothermal pellet

How much of a catalyst pellet actually works. Thiele (1939) answered it exactly
for the isothermal case; Weisz & Hicks (1962) showed that once the pellet can
heat up, the question has up to three answers at once and the pellet can run
more than a hundred times faster than its own surface conditions suggest.

- **Structures:** `S3` (1D steady BVP), `S10` (constrained continuation)
- **References:** Thiele (1939) I&EC 31(7) 916-920, doi:10.1021/ie50355a027;
  Weisz & Hicks (1962) Chem Eng Sci 17(4) 265-275, doi:10.1016/0009-2509(62)85005-2
- **Runtime:** about 1 min
- **Source:** `index.py` (paired with `index.ipynb`); model code in `model.py`,
  tested by `test_model.py`

## Results

Isothermal eta matches the exact closed form to 2.2e-4 <!-- claim: isothermal_max_rel_error_phi_le_30 -->
for phi <= 30 across slab, cylinder and sphere, with second-order convergence
in all three (the slab's observed 2.96 at phi = 5 is pre-asymptotic: its h^2
coefficient is small there; at phi = 1 it shows 2.0). At phi = 100 the error grows to 6.5e-3 <!-- claim: isothermal_max_rel_error_all_phi -->
because the surface layer is thinner than a cell; refining to n_r = 1600 brings
it to 1.9e-5 <!-- claim: isothermal_err_phi100_n1600 -->.

Non-isothermal: at beta = 0.6, gamma = 20 the fold spans
0.2973 <!-- claim: phi_fold_lower --> < phi < 0.5719 <!-- claim: phi_fold_upper -->,
and phi = 0.4123 <!-- claim: phi_probe --> has three steady states, a
38 <!-- claim: ignited_over_extinguished -->x spread at a single Thiele modulus.
Reference and pymrm agree on all three branches, including the unstable middle
one:

| branch       | shooting reference (root-found) | pymrm (Newton at fixed phi, n_r = 400) |
|--------------|--------------------|--------------------|
| extinguished | 1.1731 <!-- claim: eta_extinguished --> | 1.1731 <!-- claim: eta_extinguished_pymrm --> |
| unstable     | 6.2055 <!-- claim: eta_unstable --> | 6.2054 <!-- claim: eta_unstable_pymrm --> |
| ignited      | 44.450 <!-- claim: eta_ignited --> | 44.450 <!-- claim: eta_ignited_pymrm --> |

The branch is traced by `pymrm.continuation` in the centre concentration, with
phi^2 as an unknown, so it is single-valued and has no folds in the swept
parameter. The folds are located on the shooting reference, refined from a scan.

The page's own warm-start loop traces the same branch to below 1e-9 in phi^2. The
rate is not clipped at zero and the solves pass no bounds (pymrm pitfall P9).
Reading the folds off 900 samples and interpolating the states instead is
computed beside the refined values: it shifts the states by up to 3e-3, which
would look like a discretisation error.

## Data

**Provenance tier 6 - not experimental.** Thiele (1939) is analytical and
Weisz-Hicks (1962) is computational; neither reports measurements. Validation is
against exact solutions and an independent shooting reference. Experimental
effectiveness-factor data exists elsewhere and would make a good companion page.

## Corrections

- 2026-10-30: stated that eta exceeds 1000 and the pellet runs a thousand times
  faster than its surface conditions; the largest eta on the diagram is 151, at
  beta = 0.8.
- 2026-10-30: stated that the middle branch is unreachable by Newton however the
  guess is chosen; Newton at fixed phi converges there from a guess close to it.
- 2026-10-30: stated that the discretisation is pinned to 1e-4; it is 2.2e-4.
- 2026-10-30: quoted a 0.3 % pymrm-vs-reference gap at the probe modulus (6.225
  against 6.207); that was linear interpolation, and the difference is about 1e-5.
- 2026-10-30: printed one mean convergence order, 2.35; each geometry is second
  order (the slab's 2.96 at phi = 5 is pre-asymptotic).
- 2026-10-30: stated the upper fold as 0.57179, read off samples; refined, it is
  0.57186.
