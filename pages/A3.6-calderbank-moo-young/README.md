# A3.6: Calderbank & Moo-Young, continuous-phase heat and mass transfer in dispersions

Calderbank, P. H. & Moo-Young, M. B., *The continuous phase heat and mass-transfer
properties of dispersions*, **Chemical Engineering Science 16** (1961) 39-54,
doi:[10.1016/0009-2509(61)87005-X](https://doi.org/10.1016/0009-2509(61)87005-X).

> **Citation trap on this file.** The PDF's text layer prints the header as
> *"1981 … Vol 19"* and corrupts one running footer to 1991. The citation above
> was established from the article's **own printed header on a 300 ppi page
> image**, cropped and enlarged with nearest-neighbour interpolation:
> *"Chemical Engineering Science, 1961, Vol 16, pp 39 to 54 Pergamon Press Ltd,
> London"*. Three running footers independently read *"Vol 16, Nos 1 and 2
> December, 1961"*, and the PII `0009-2509(61)87005-X` encodes 61.

- **Structures:** `S3` (1-D steady BVP), `S4` (1-D marched PDE)
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (the stagnant sphere, the boundary layer and its two closed forms),
  tested by `test_model.py`. The paper's algebra, the ladders and the break table
  stay in the notebook.
- **Runtime:** about 35 s.

## What the page shows

Two correlations for the continuous-phase mass-transfer coefficient in a
dispersion whose particles move freely under gravity, differing only in the
Schmidt exponent: `2/3` for rigid interfaces (hindered flow, Froessling) and `1/2`
for mobile ones (unhindered flow, Higbie). The page

1. **reproduces the paper's own two theoretical derivations** of eq. (1)'s
   constant: 0.276000 <!-- claim: c28_value --> against a printed `0.28`
   (Froessling and Allen) and 0.339599 <!-- claim: c34_value --> against a
   printed `0.34` (Friedlander and Stokes), neither fitted;
2. **rebuilds the exponent split in pymrm** as a concentration boundary layer
   over an interface that is held still or allowed to move, a linear march
   checked against the Lévêque and Higbie closed forms, each grid axis on its own
   ladder;
3. **maps the transition** the paper describes and never models: one group
   `Λ = u_s/(γ̇ δ_Lev)` controls it, halfway at
   `Λ` = 0.2721 <!-- claim: lambda_half -->;
4. **reports two places where the paper's own arithmetic does not close** (eq. (6)
   cannot produce eq. (7), and eq. (9) does not meet Chilton-Colburn at
   `N_Re = 10,000` nor within `±12 %`) and **corrects a printed table exponent** from
   the paper's own correlation, with a control that passes.

## Agreement

| Check | Result |
|---|---|
| Froessling + Allen against the printed `0.28` | -1.43 % <!-- claim: c28_rel_dev -->; the paper's own `0.28` and ¼ imply a Froessling coefficient of 0.5600 <!-- claim: froessling_implied_by_paper --> |
| Friedlander + Stokes against the printed `0.34` | -0.12 % <!-- claim: c34_rel_dev --> |
| Allen's printed ¼ against his drag law `C_D = 10 Re^(-1/2)` | 0.2610 <!-- claim: allen_quarter_from_drag_law -->, +4.4 % <!-- claim: allen_quarter_rel_dev --> |
| eq. (5) + Chu → the printed `P/v`, eq. (10) exponents, page 44's Rayleigh form | three exact identities, each below 1e-12 <!-- claim: chu_identity_resid --> |
| eq. (6) prefactor `0.13 × 35.4^(1/4)` | 0.317098 <!-- claim: eq6_prefactor_derived --> against `0.318`, -0.28 % <!-- claim: eq6_prefactor_rel_dev --> |
| eq. (7) from eq. (6) at `Σ = 0.41` | **does not close**: as printed -48.0 % <!-- claim: eq7_from_printed_eq6_rel_dev -->, with `Σ^-1` -5.7 % <!-- claim: eq7_from_sigma_inverse_rel_dev --> |
| eq. (9) prefactor `0.13 × 0.0396^(1/4)`, the control | 0.0579919 <!-- claim: eq9_prefactor_derived --> against `0.058`, -0.014 % <!-- claim: eq9_prefactor_rel_dev --> |
| eq. (9) against Chilton-Colburn | **does not hold as stated**: crossing at `N_Re` = 4480 <!-- claim: pipe_cross_Re -->, not 10,000; largest deviation over 2000-20,000 is 15.2 % <!-- claim: pipe_maxdev_over_range --> |
| small-bubble `D_L` printed `× 10^-4` | eq. (1) stands 21.2 <!-- claim: eq1_window_hi_printed_exponent --> times above the row's own `k_L` ceiling; `× 10^-6` gives 0.983 <!-- claim: eq1_window_hi_over_printed --> of it |
| the control: large-bubble row, eq. (2), as printed | 0.818 <!-- claim: eq2_window_lo_over_printed --> and 0.896 <!-- claim: eq2_window_hi_over_printed --> of the printed window ends, no adjustment |
| pymrm sphere, the stagnant floor `2.0` | -3.34e-5 <!-- claim: sh_stagnant_rel_dev --> at 800 cells, second order; the Richardson limit below 1e-6 <!-- claim: sh_richardson_rel_dev --> from 2 |
| pymrm boundary layer against Lévêque and Higbie | -1.43e-4 <!-- claim: leveque_rel_err_raw --> and -1.60e-4 <!-- claim: higbie_rel_err_raw --> at 800 cells and 800 steps; first order in the step, second in the cells; the two-axis limit below 1e-7 <!-- claim: leveque_rel_err_xy_extrap --> |
| exponents `d ln k / d ln D` in the two limits | 2/3 and 1/2 to round-off: fixed by dimensional analysis, an implementation check |

Deviations are `(computed − printed)/printed`.

## What pymrm adds

The operators of the two linear problems, the spherical shell and the marched
layer, and the flux reads; with them the transition becomes a number. The paper
describes the transition three times (Fig. 3's intermediate bubbles, the decline
of large-bubble coefficients in deep pools, and the effect of surface-active
agents) and models it nowhere. Written as `u(y) = u_s + γ̇ y`, one group `Λ`
controls it, the exponent falls at every step of the scan from 2/3 to 1/2, it is
halfway at `Λ` = 0.2721 <!-- claim: lambda_half --> and within `10 %` of the mobile end
by `Λ` = 1.507 <!-- claim: lambda_within_10pct_mobile -->. The paper's two
correlations are the two ends of one curve. The paper's algebra, which carries the
rest of the page, is NumPy and sympy and owes pymrm nothing.

## What this page does NOT do

- **It does not judge `0.31`, `0.42` or `0.13`.** Those are regression lines through
  Figs 1, 2 and 6. The evidence about them is the scatter in those figures, and no
  figure is digitised here.
- **The corner evaluations of eqs. (1) and (2) against the Range of Variables
  table are a goodness of fit, not a validation**: those constants were fitted to
  the very runs whose `k_L` ranges they are compared with. What they can see is a
  wrong power of ten.
- **The boundary layer does not predict the constants.** It reproduces the
  exponents and the interpolation; `0.31` and `0.42` absorb bubble shape, slip
  velocity and swarm hindrance, none of which a flat interface contains. Its two
  pure limits are fixed by dimensional analysis, so the exponents test the
  implementation rather than the physics.
- **It shares a limit with `A3.3`.** The mobile branch reduces to the penetration
  problem Danckwerts' page solves; this page checks against the closed form, not
  against `A3.3`'s output. No number, dataset or sentence from `A3.1`, `A3.3` or
  `A3.4` is reused here.

## Corrections

- 2026-11-09: said Fallat's measured j-factor was compared "over three decades".
  The comparison runs from `N_Re = 100` to `10,000`, two decades; over that range
  eq. (7) is within 4.6 % <!-- claim: fallat_dev_printed_eq7 --> of Fallat, the
  `Σ^-1` reading within 6.4 % <!-- claim: fallat_dev_sigma_inverse --> and the
  chain's own exponents within 28.0 % <!-- claim: fallat_dev_chain -->.
- 2026-11-09: listed "doubling the far field by 8.4e-7" among the boundary
  layer's error sources. That study held the cell count fixed, so the doubled
  height doubled every cell. With the cell error removed by extrapolation on a cell
  ladder at each height, doubling the far field moves the mean coefficient of
  either branch by less than 1e-8 <!-- claim: bl_domain_effect -->.
- 2026-11-09: said the boundary-layer march gave "Observed orders 1.11 and 1.09
  in the march step". Those orders were read on a ladder that refined the cells
  together with the step. Refined alone, the step gives orders 0.9994
  <!-- claim: leveque_x_order_min --> to 0.9997 <!-- claim: leveque_x_order_max -->
  on the Lévêque branch and 0.9997 <!-- claim: higbie_x_order_min --> to 0.9998
  <!-- claim: higbie_x_order_max --> on the Higbie branch, the first order of
  backward Euler.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```

No PDF needed: both datasets are committed CSVs in `data/`.
