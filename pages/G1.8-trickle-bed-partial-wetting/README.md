# G1.8: Trickle-bed reactor with partial catalyst wetting

Herskowitz, M. and Smith, J. M. (1983). *Trickle-bed reactors: a review.*
AIChE Journal **29**(1) 1–18. doi:10.1002/aic.690290102

Reproduces **Figure 6**, the gain in overall effectiveness factor when part of
a catalyst pellet is dry and the limiting reactant is volatile, from the
paper's own Table 1 (sphere row, page 4) and Eqs. 19–21 (page 8), with
$\alpha_{gs}\to\infty$ and $C^*_L = 1$ as stated.

- **Source:** `index.py` (paired with `index.ipynb`); model code in `model.py`
  (`Closures` for Eqs. 20 and 21, the Table 1 closed form, the root finders
  `reconstruct_L` and `alpha_gs_for`, and the pymrm pellet `WettedSphere`),
  tested by `test_model.py`.
- **Runtime:** about ten seconds.
- `data/`: four fitted lines, not markers, with their provenance sidecar (which
  also records the abscissa calibrations and the independent re-digitisation
  used for the slope range).

## What this page is

**Tier 6.** Figure 6 contains no measurements; its four lines are the authors'
own computed model output. Reproducing them tests a transcription and a chain of
algebra. Nothing on this page is experimental validation.

## The headline

The figure **does not reproduce as printed**. At $\phi = 10$ the paper's own
chain lands at +1.0 <!-- claim: chi10_dev_printed_curve1_pct --> %,
+69 <!-- claim: chi10_dev_printed_curve2_pct --> %,
+314 <!-- claim: chi10_dev_printed_curve3_pct --> % and
+152 <!-- claim: chi10_dev_printed_curve4_pct --> % of where the four drawn
curves sit: right for one curve out of four.

Solving for the $L_m$ that puts each drawn curve where it is gives
**0.508 <!-- claim: L_rec_curve1 -->, 2.046 <!-- claim: L_rec_curve2 -->,
7.061 <!-- claim: L_rec_curve3 --> and 10.18 <!-- claim: L_rec_curve4 -->**.
Three of those are the legend's own printed 0.50, 2.0 and 7.0, reused in a
different order. Reading the four curves as $L_m$ = 0.50, 2.0, 7.0 and ≈10
brings all four to +1.0 <!-- claim: chi10_dev_shifted_curve1_pct --> %,
+1.9 <!-- claim: chi10_dev_shifted_curve2_pct --> %,
+1.6 <!-- claim: chi10_dev_shifted_curve3_pct --> % and
+6.4 <!-- claim: chi10_dev_shifted_curve4_pct --> % at $\phi = 10$, with
**nothing fitted** for the first three.

**The legend block is misaligned by one row against its curves**, from the
second row down: its $L_m = 1.0$ entry belongs to no drawn curve, and the fourth
drawn curve has no entry. Endorsed by the gallery maintainer on 2026-08-02 from
a private decision artifact. The page shows both readings side by side.

Two further readings of the same figure support it. **Neither is an independent
witness**: both are the same four measured curve positions, re-expressed.

1. The three inter-curve **gaps**, in decades, at $\phi = 10$: the shifted
   reading matches them to 0.020 <!-- claim: gap_miss_shifted_decades -->
   decades, the printed reading misses by up to
   0.389 <!-- claim: gap_miss_printed_decades -->, and the gaps drift by
   0.031 <!-- claim: figure_gap_drift_over_window_decades --> decades across the
   window. This rules out "printed reading right, $\chi$ calibration wrong", and
   it needs the ordinate's decades-per-pixel.
2. Because $\chi = g(L_m)\,h(\phi)$ **separates**, those same gaps force
   $d\log g/d\log L_m$ with no model at all. As printed the forced values are
   non-monotone, so no product of powers of $L_m$ supplies them; shifted they
   agree with the paper's own closures within
   5.4 % <!-- claim: dlogg_forced_shifted_worst_rel_err -->. No
   reparameterisation rescues the printed reading; only a reassignment does.

## The residual the shift does not fix

The four drawn lines measure log–log slopes
**1.098 <!-- claim: slope_figure_min --> to 1.131 <!-- claim: slope_figure_max -->**
(mean 1.111 <!-- claim: slope_figure_mean -->); the model's best straight-line
fit over the same window is **1.077 <!-- claim: slope_model_fitted -->**. $L_m$
cannot produce that: the chain factorises into an $L_m$-dependent prefactor
times a $\phi$-dependent shape, so every $L_m$ gives the identical slope.

The **sign** is robust; the **size** is not. The dominant uncertainty is the
abscissa calibration, not the fit rms: the four labelled $\phi$ verticals do not
sit on a common log ruler, and an independent 600 dpi re-digitisation differs
from this CSV by exactly the ratio of the two calibrations. Across defensible
calibrations the gap spans **+0.007 <!-- claim: slope_gap_lo_over_abscissa_calibrations -->
to +0.046 <!-- claim: slope_gap_hi_over_abscissa_calibrations -->**
(+3 <!-- claim: slope_gap_lo_e2e_pct --> % to
+19 <!-- claim: slope_gap_hi_e2e_pct --> % end to end). **Reported as
unexplained, with its size as a range.**

## Two printed defects reported

1. **Page 8 cites "Table 2"** for the approximate spherical solution. Table 2
   holds the pressure-drop constants $\beta$ and $\gamma$; the sphere row is in
   **Table 1**, page 4. Checked on 600 dpi renders of both pages and against the
   notation list.
2. **Figure 6's legend is misaligned by one row** against the curves, as above.

## Validation

Every check is followed by a deliberate break that re-runs it with the defect
injected; the notebook's defect-injection table lists every row.

- Eq. 20 against the $f_e$ printed inside the figure: within
  0.005 <!-- claim: eq20_max_abs_err_vs_printed -->.
- The pymrm finite-volume pellet against Table 1's sphere row: within
  3.1e-4 <!-- claim: pymrm_vs_table1_maxdev --> over the digitised window,
  plus a surface-flux identity read on the boundary face. The scheme is second
  order (2.00 <!-- claim: grid_order_finest_phi2 --> at $\phi = 2$); at
  $\phi = 45.5$ the observed rate reaches 2.84 <!-- claim: grid_order --> and
  falls on finer grids, a pre-asymptotic excess.
- The reconstruction check has a real blind spot and the page states it. Its
  tolerance is 3.8 <!-- claim: L_rec_tol_curve3_pct --> to
  11 <!-- claim: L_rec_tol_curve1_pct --> % in $L_m$, and
  4 <!-- claim: recon_check_n_unrejected --> injected defects of Eq. 21's second
  term stay under even the tightest bar; one, `6.91 → 7.91`, makes the metric
  look *better* (1.72 % <!-- claim: recon_check_best_unrejected_defect -->
  against the undamaged 2.30 % <!-- claim: recon_check_baseline -->). At
  $L_m \approx 1$ the $1.05\,L_m^{0.3}$ term carries
  87 % <!-- claim: eq21_first_term_resistance_share_at_Lm1 --> of
  $1/\alpha_{gLs}$. The check is also weak against a ±5 % ordinate-scale error,
  which is why the gap test is kept.
- The collapse identity does not move under any defect: kept and labelled
  powerless.

## Rebuilding the page

```bash
# from the repository root
jupytext --sync pages/G1.8-trickle-bed-partial-wetting/index.py   # after editing index.py
python scripts/run_kernels.py pages/G1.8-trickle-bed-partial-wetting
GALLERY_OFFLINE=1 PYTHONPATH=shared python -m pytest -q pages/G1.8-trickle-bed-partial-wetting/test_model.py
```

## Corrections

- 2026-10-30: stated that the observed grid order rises from 2.1 to 2.8 as the
  surface layer resolves, consistent with a second-order scheme whose coarse-grid
  error is dominated by the boundary layer. The observed rate at $\phi = 45.5$
  rises from 2.13 <!-- claim: grid_order_coarsest_phimax --> to a plateau of
  2.84 <!-- claim: grid_order --> ($n_u$ = 400 to 1600) and then falls on finer
  grids, to 2.66 <!-- claim: grid_order_finest_phimax --> between $n_u$ = 3200
  and 6400; at $\phi = 2$ the order is 2.00 <!-- claim: grid_order_finest_phi2 -->.
  The scheme is second order; the 2.8 is a pre-asymptotic rate.
- 2026-10-30: stated that the Newton residual is reported next to every solve
  and guards against a non-converged solve. It was printed for one solve, and no
  solve checked convergence; every solve now raises when `newton` does not
  converge (all converged).
- 2026-10-30: the Reuse table said `WettedSphere(nu=0)` gives the slab and
  `nu=1` the cylinder. With the sphere's length ratio 3 they miss Table 1's slab
  and cylinder rows by tens of per cent; with `n_scale` 1 and 2 they agree to
  1.1e-5 <!-- claim: reuse_slab_vs_table1 --> and
  9.6e-5 <!-- claim: reuse_cylinder_vs_table1 -->. It also called the fully
  wetted pellet `B1.1`; `B1.1` is its $\alpha \to \infty$ (Dirichlet) member.
- 2026-10-30: the data sidecar gave Eq. 20 at $L_m = 7$ as 0.936; it is 0.935.
  It also listed the re-digitisation's smallest refit gap as a measurement; it
  is a derived model result and is now labelled so.
