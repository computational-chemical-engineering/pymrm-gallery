# A4.9 — Duncan-Toor ternary gas diffusion

Two bulbs of gas joined by a capillary. Nitrogen starts 0.06 mole fraction apart,
closes that gap, then keeps going *against* its own gradient until it is twice as
far apart the other way, and finally stops dead at the moment its gradient is
largest. Fick's law can produce none of that.

- **Structures:** `S3` (1D steady BVP), `S9` (implicit multicomponent flux)
- **Reference:** Duncan & Toor (1962), AIChE J 8(1) 38-41, doi:10.1002/aic.690080112
- **Runtime:** ~48 s

## Agreement

0.59 mole % mean absolute deviation over 28 digitised points, nothing fitted.
The paper reports 0.45 mole % for its own Maxwell-Stefan predictions and quotes
2.6 mole % experimental error; about 0.5 mole % of our figure is digitisation
error rather than model error.

Point O (osmotic diffusion) at 0.96 h, diffusion barrier at 7.52 h.

## The comparison with the data is the *only* guard

The page used to open its Validation section with *"Four checks. The first three
do not use the experimental data at all"*, which reads as a strength and is the
opposite of one. Two of those three were the same algebraic identity:
`solve_capillary` drives `div(N)` to `tol=1e-12` with `construct_div(nu=0)`, so on
a Cartesian grid the flux-uniformity check **is** the converged Newton residual;
and with the flux uniform, `N_0 = N_L`, so the mole balance is zero iff the
uniformity check passes.

Neither can see a wrong number. Measured on the page's own `simulate()`:

| injected defect | flux non-uniformity | mole balance | vs data | N₂ crossing | final Δx(N₂) |
| --- | --- | --- | --- | --- | --- |
| as published | 1.58e-14 | 2.16e-16 | 0.59 % | 1.00 h | +0.0849 |
| `D(N₂–CO₂)` 16.8 → 61.8 (digit swap) | 1.93e-14 | 4.33e-16 | **2.28 %** | **3.28 h** | **+0.0043** |
| `D(N₂–CO₂)` 16.8 → 1.68 | 2.16e-14 | 4.33e-16 | **1.30 %** | 0.92 h | +0.2090 |
| `D(H₂–N₂)` 83.3 → 8.33 | 1.68e-14 | 1.73e-16 | **5.19 %** | **never** | −0.1106 |
| `LA_EFF` 258.1 → 285.1 | *bit-identical* | *bit-identical* | 0.80 % | 1.08 h | +0.0934 |
| `vol_2` 78.63 → 87.63 | *bit-identical* | 3.27e-16 | 0.73 % | 1.04 h | +0.0902 |

A digit transposition in the one diffusivity that drives the whole phenomenon ends
the reversal 20× smaller and 2.3 h later, and both residuals stay at machine
precision. They are kept and labelled `STRUCTURAL`, because they do catch what
nothing else here catches — an unconverged Newton solve (1.7e-3) and a wrong `nu`
(0.98).

Two further facts worth recording. Both metrics sit **below
`check_agreement.py`'s `ABS_FLOOR = 1e-12`**, so CI does not compare them at all.
And the flux-uniformity metric *improves* on a coarser grid: from `n_z` 10 to 80
the actual error falls 88-fold while the residual rises from 2.2e-15 to 4.2e-14.
It is round-off in a linear solve, not accuracy.

**This block is shared with `A4.2`**, which borrowed this page's dataset and
copied the validation with it — the only two pages in the repository carrying
`flux_nonuniformity` and `mole_balance_error`. Both are now labelled. A check's
power depends on the physics, not on the code, so it does not survive a copy.

## Both numerical knobs are now refined

`n_z = 40` had never been refined. Grid: observed order 2.06, spatial error
1.5e-6 in mole fraction — about 3900× below the model-vs-data deviation. Time:
observed order 1.15 (backward Euler), 2.2e-4 at the production `n_steps = 400`.
Time is the dominant numerical error, 147× the grid's, and both are far below the
0.59 mole % the comparison reports.

## Data

Digitised from Figure 2 by programmatic marker extraction, validated by closure,
monotonicity and the closed-cell exchange balance. See
`data/duncan-toor-1962-run1.meta.yaml` for the method and its limits.

All model parameters were verified against the paper text. Note the paper uses an
*effective* (L/A) of 258.1 1/cm measured from binary runs, not the geometric
251.1 1/cm.
