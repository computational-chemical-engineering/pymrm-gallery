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

## Refresh of 2026-10

- New layout: `index.py`, `model.py`, `test_model.py`, bootstrap cell, pymrm
  pin 3.2.0, every number registered with `Claims` (49 claims, 55 break rows,
  each re-running its computation).
- The page's own warm-start loop is replaced by `pymrm.continuation`; the old
  loop stays on the page as a comparison (same branch to 4e-12 in phi^2). The
  library routine in pymrm 3.1.0 stopped on the old sweep, whose two segments
  both contained 0.1 (as 1 - 0.9 and as 0.1, 3e-17 apart): its secant predictor
  divided by that step. Reported from this page and fixed in pymrm 3.2.0, which
  traces the old sweep through; the sweep still drops the repeated value.
- The rate is no longer clipped at zero and the solves no longer pass
  `bounds=(0, inf)` (pymrm pitfall P9); no number changed.
- The folds are refined instead of read off 900 samples (upper fold
  0.57179 before, 0.57186 now), and the steady states at the probe are
  root-found and solved by Newton instead of interpolated. The old table
  (reference 1.177 / 6.207 / 44.450, pymrm 1.173 / 6.225 / 44.450) suggested a
  discretisation error of 0.3 %; it was interpolation, and the true difference
  is about 1e-5. The page computes both methods side by side.
- Correction (2026-10-30), also on the page: the earlier version said eta
  exceeds 1000 and the pellet runs a thousand times faster (the page's largest
  eta is 151, at beta = 0.8); that the middle branch is unreachable by Newton
  however the guess is chosen (from a guess close to it, taken off the traced
  branch, Newton at fixed phi converges there in a few iterations); that the
  discretisation is pinned to 1e-4 (it is 2.2e-4); it quoted the 0.3 %
  interpolation gap above as a pymrm-vs-reference difference; and it printed
  one averaged order, 2.35, which hid the per-geometry orders (all second
  order; the slab's 2.96 at phi = 5 is pre-asymptotic).
- The secant ratio the old sweep produces, 2.4e14 (exact from the sweep's
  values), is shown with the result that pymrm 3.2.0 traces through it.

## Data

**Provenance tier 6 - not experimental.** Thiele (1939) is analytical and
Weisz-Hicks (1962) is computational; neither reports measurements. Validation is
against exact solutions and an independent shooting reference. Experimental
effectiveness-factor data exists elsewhere and would make a good companion page.
