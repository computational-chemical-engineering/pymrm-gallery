# I1.3 — Heck, Wei & Katzer 1976: what the 1-D monolith model is adequate for

The paper's two claims, tested: Table 1's asymptotic Nusselt numbers derived from
scratch, the Nusselt singularity root-found downstream of the cause the text gives
it, and the 1-D/2-D adequacy claim re-measured with root-found thresholds.

- **Structure:** `S6` (2-D PDE)
- **Reference:** Heck, R. H., Wei, J. & Katzer, J. R., *Mathematical Modeling of
  Monolithic Catalysts*, AIChE Journal **22**(3) 477-484 (1976),
  doi:10.1002/aic.690220310
- **Source of the page:** `index.py` (paired with `index.ipynb`); model code in
  `model.py` (`duct` and `square_spectral` for Table 1, `Graetz` for the developing
  tube problem with a deferred-correction and a `ConvFluxTVD`/`newton` solver,
  `Voltz` and `Channel` for the reacting 1-D and 2-D channel with light-off
  root-found inside its axial step), tested by `test_model.py`.
- **Data:** tier 6, the paper's own printed values: Table 1 and every number stated
  in running text, read off the 300 ppi page bitmaps, each file with its
  `.meta.yaml` sidecar.
- **Runtime:** about 10 minutes, half of it the defect-injection table, which
  re-runs the reacting model four times with one defect each.

## Scope: figures are out

Figures are used only for characters printed inside their frames (annotation
boxes, abscissa tick labels, one curve label, one case label). The notebook lists
and counts every such row from the CSV. No curve is traced, no profile shape or
light-off position is compared with the paper's own, and none of the conversion
curves of Figs. 4 to 6 and 12 is checked.

## Table 1, derived rather than quoted

| cell | printed | this page | route |
|---|---|---|---|
| circle, $Nu_H$ | 4.364 | 48/11 | closed form |
| triangle, $Nu_H$ | 3.111 | 28/9 | closed form |
| square, $Nu_H$ | 3.608 | 3.607950 <!-- claim: nu_h1_square_pymrm --> | pymrm FV and Fourier-Galerkin |
| circle, $Nu_T$ | 3.657 | 3.656793 <!-- claim: nu_t_circle_pymrm --> | pymrm FV, equal to `A3.15`'s Graetz limit |
| square, $Nu_T$ | 2.976 | 2.977523 <!-- claim: nu_t_square_pymrm --> | pymrm FV and Fourier-Galerkin |

Four agree with every printed digit. The fifth does not: `2.976` implies
[2.9755, 2.9765], which excludes the computed value, and the page does not
attribute it, because Table 1's source (Shah and London 1971) was not consulted.

## Claim 1: the singularity is real, the stated cause is not

At the zero of the imposed wall flux the Nusselt number of eq. (19) has a
vanishing numerator, so it passes through zero. The pole is the zero of
$T_w - \bar T_G$, 8.7 % <!-- claim: nu_pole_downstream_fig9 --> further
downstream in both printed flux functions. The paper's next paragraph gives the
correct mechanism, so this is a wrong one-sentence summary, not a missed effect.
Its sign statements bracket the pole in (0.032, 0.033); the solve puts it
4.1 <!-- claim: nu_pole_below_bracket_centre_pct --> % below the centre of that
bracket, and the paper's own stated discretisation does not reach the bracket.

## Three meanings of $X^{*}$

The Notation's $4x/(D\,Re\,Pr)$ (eq. 9), eq. (15)'s $x/(D\,Re\,Pr)$ (every figure
abscissa) and the Grigull-Tratz correlations' $1000\,x/(D\,Re\,Pr)$, each settled
against a computation. The printed $\bar T^{*}(0.25) = 0.0224$ settles the second
and is itself 5.9 % <!-- claim: Tbar_X025_paper_grid_error --> above the converged
value: it measures the paper's ten-point radial grid.

## Claim 2: the paper's thresholds are sampled; these are root-found

The paper's two light-off windows are quoted from sampled cases at two different
channel lengths (Fig. 2 to 0.20, Fig. 11 to 0.50). Root-found at each figure's
length the models reproduce both statements, and the gap between them,
-10.83 <!-- claim: T_windows_length_artefact_C --> degC, is
-13.42 <!-- claim: T_windows_length_term_C --> degC of length and
+2.59 <!-- claim: T_windows_model_term_C --> degC of model. On a common length the
2-D model needs a higher inlet temperature than either 1-D model,
+4.13 <!-- claim: T_low_2d_minus_1d_const_L025 --> degC over the constant-Nu 1-D
at $L^{*} = 0.25$, of which +2.07 <!-- claim: T_low_shape_effect_L025 --> degC is
channel shape (square-channel Nu against a round tube).

The "less than a 2 degC difference in inlet gas temperature", as the ratio of the
equivalent inlet shift to the claim: for the Grigull-Tratz 1-D model of Fig. 12,
1.43 <!-- claim: equiv_shift_GT_over_claim_at316 --> at a 316 degC inlet and
1.09 <!-- claim: equiv_shift_GT_over_claim_min --> to
3.41 <!-- claim: equiv_shift_GT_over_claim_max --> over 310 to 340 degC; for the
constant-Nu model of the paper's parametric figures,
2.86 <!-- claim: equiv_shift_const_over_claim_at316 --> at 316 degC and
2.26 <!-- claim: equiv_shift_const_over_claim_min --> to
7.74 <!-- claim: equiv_shift_const_over_claim_max --> over 310 to 340 degC.

Two dimensional groups ($\Delta T_{AD}$ and $\beta$) are reconstructed from
ordinary gas properties, so no absolute threshold on the page is the paper's
number; only differences between models on the same parameters are claimed.

## Corrections

- 2026-10-31: stated that the "less than 2 degC" claim is 2 to 5 times too small for
  the constant-Nu 1-D model; the upper end is wrong. Over 310 to 340 degC the shift
  is 2.3 <!-- claim: equiv_shift_const_over_claim_min --> to
  7.7 <!-- claim: equiv_shift_const_over_claim_max --> times the claim. The notebook's
  "roughly 2.0x" was the ratio to the Grigull-Tratz shift, not to the claim; at
  316 degC the ratio to the claim is 2.9 <!-- claim: equiv_shift_const_over_claim_at316 -->.
- 2026-10-31: stated that the Grigull-Tratz 1-D shift is at or above 2 degC under
  every reconstruction variation, with the smallest +1.84 degC; that value is below
  the claim, and the shift exceeds it in 4 of 5 variations (smallest
  +1.85 <!-- claim: recon_min_equiv_shift_GT --> degC).
- 2026-10-31: `meta.yaml` stated that plausible variations of the reconstructed
  groups move the 1-D/2-D difference by under a degree; the constant-Nu equivalent
  shift spans 4.55 <!-- claim: recon_spread_equiv_shift_C --> degC
  (3.96 <!-- claim: recon_min_equiv_shift_const --> to
  8.52 <!-- claim: recon_beta_over_1_5_shift_const -->).

## Cite the source, not this page

Heck, R. H., Wei, J. & Katzer, J. R., *Mathematical Modeling of Monolithic
Catalysts*, AIChE Journal **22**(3) 477-484 (1976), doi:10.1002/aic.690220310.

Shah and London (1971), Voltz et al. (1973) and Grigull and Tratz (1965) were
not consulted; all three are used exactly as Heck et al. reprint them.
