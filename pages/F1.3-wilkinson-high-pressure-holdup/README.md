# F1.3 Wilkinson high-pressure bubble-column gas holdup

The correlation is Wilkinson, P. M., Spek, A. P. and van Dierendonck, L. L.,
*Design parameters estimation for scale-up of high-pressure bubble columns*,
*AIChE Journal* **38**(4) 544–554 (1992),
[doi:10.1002/aic.690380408](https://doi.org/10.1002/aic.690380408).

**That paper was not available.** The equations used here are its reprint as
Eqs. 1–4 of Krishna, R. and Ellenberger, J., *AIChE Journal* **42**(9) 2627–2634
(1996), [doi:10.1002/aic.690420923](https://doi.org/10.1002/aic.690420923), each
verified on a 600 dpi render of that paper. The page says so, and so does the
reference block.

## Contents

| file | what it is |
| --- | --- |
| `index.py`, `index.ipynb` | the page, as a jupytext pair; edit `index.py`. Runs in about 4 s |
| `data/krishna-ellenberger-1996-table2-liquids.csv` | Table 2's liquids, transcribed from a 600 dpi render |
| `agreement.json` | every number the page states, written by its Claims cell and checked by CI |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Two further datasets are loaded **cross-page** from
`pages/F1.4-krishna-ellenberger-holdup/data/` and are not copied here: the
digitised Figure 11 markers and the parameters file holding the gas densities,
Reilly's *B* and the Table 3 deviations. Read
`krishna-ellenberger-1996-fig11.meta.yaml` on that page before using the `gas`
column: a maintainer review found the marker positions trustworthy and the
shape-derived gas labels not, and nothing on this page depends on series
identity.

## What this page is, next to F1.4

`F1.4` tests Krishna and Ellenberger's replacement correlation, Eq. 19, against
their Figure 11, and uses Wilkinson only as the foil, evaluating its
rise-velocity equations at the plotted abscissa and never touching **Eq. 2**, its
transition equation, at all. This page runs the correlation as a correlation,
from operating conditions, and Eq. 2 turns out to be where everything happens.
The two pages share one deviation number, Wilkinson's +64 % <!-- claim: wilkinson_eps_b_bias -->,
which the page checks against `F1.4`'s committed value; this page retypes Eqs. 3
and 4 from the same source and applies them to the same rows, so that agreement is
a **transcription cross-check between the pages, not an independent computation**.

## Headline result

**This is the authors' own finding, made specific.** Krishna and Ellenberger
write that "the Wilkinson correlation severely underpredicts the values of the
voidage and gas velocity through the dense phase", and compare Eq. 2 against
their measured `ε_df` themselves. At the gas densities they tested, Wilkinson's
Eq. 2 gives a transition holdup of 0.0028 <!-- claim: eps_trans_wilkinson_air_tetradecane -->
for air–tetradecane and 1.5e-8 <!-- claim: eps_trans_wilkinson_helium --> for helium,
which is 1.7 <!-- claim: decades_low_air --> decades below the measurement at air,
1.4 <!-- claim: decades_low_argon --> at argon and 6.8 <!-- claim: decades_low_helium --> at helium.
The homogeneous regime effectively does not exist in the correlation at atmospheric
density, so the small-bubble population is empty.

**The gap does not close at higher density.** Reilly's Eq. 8 exceeds Eq. 2 at
every density: the ratio has one stationary point, a minimum of
2.85 <!-- claim: eq8_over_eq2_min_ratio --> near 29 <!-- claim: rho_G_at_eq8_over_eq2_min --> kg/m³,
so the curves never cross. And the target is not fixed: Eq. 2 crosses a *fixed* 0.15
at 14 <!-- claim: rho_G_where_eq2_reaches_0p15 --> kg/m³, but the paper states that increasing
gas density significantly increases the dense-phase voidage, and the three points
read off its Figure 6(a) scale as ρ_G^0.30 <!-- claim: fig6_eps_df_exponent -->, putting the
measured level at that density near 0.32 <!-- claim: measured_level_at_rho_0p15 --> where
Eq. 2 would still be a factor 2.1 <!-- claim: factor_low_at_rho_0p15 --> low. **That density is a
lower bound on where Eq. 2 could become adequate, not a crossover.**

That is one defect. The other is separate and the page separates them: at the
measured excess velocity Eq. 2 does not enter at all, and `ε_b` still comes out
64 % <!-- claim: wilkinson_eps_b_bias --> high on average, and high at every marker, which is Eqs. 3 and 4
jointly. `V_b` is 39 % <!-- claim: V_b_shortfall_all_air --> low on average (the uniform shortfall that gives the mean bias; the mean of each marker's own shortfall is a different statistic, 37 % <!-- claim: V_b_shortfall_marker_mean_all_air -->), and `V_b`'s first term
*is* Eq. 3's `V_small`, 46 % <!-- claim: V_small_fraction_of_V_b_mean --> of `V_b` on average
(23 % <!-- claim: V_small_fraction_of_V_b_min -->–77 % <!-- claim: V_small_fraction_of_V_b_max -->,
the high end at the lowest excess velocity), while Eq. 4's second term is the larger one at
42 <!-- claim: eq4_second_term_larger_points --> of the 63 <!-- claim: n_points --> markers, above an excess
of 0.072 <!-- claim: eq4_terms_equal_excess --> m/s. Run end to end from a physical velocity,
which adds Eq. 2's transition error, the bias at air is +106 % <!-- claim: end_to_end_bias_air -->, and on
individual markers Wilkinson's `ε_b` moves by up to 138 % <!-- claim: abscissa_effect_max_all_markers -->
across the four gas densities (128 % <!-- claim: abscissa_effect_max_air_markers --> at air)
when the transition velocity is included, against at most 15 % <!-- claim: abscissa_convention_max_effect -->
at the fixed velocities 0.15 and 0.30 m/s. The air-density figure assumes every marker is air, which the
unreliable gas labels force; across the four gases it is +47 % <!-- claim: wilkinson_eps_b_bias_all_helium -->
to +79 % <!-- claim: wilkinson_eps_b_bias_all_SF6 -->. The two errors have opposite signs in the
*total* holdup and 68 % <!-- claim: total_holdup_error_cancellation --> cancels, which is why Table 3 makes
Wilkinson 7.8 <!-- claim: table3_eps_b_wilkinson_over_replacement --> times worse than its replacement on
`ε_b` and level with it on `ε`.

Meanwhile Reilly's Eq. 8, which Krishna and Ellenberger recommend instead, passes
its own authors' stated data ceiling just past SF₆, turns non-monotone at
15 <!-- claim: rho_G_at_reilly_U_trans_peak --> kg/m³ and exceeds 1 at
70 <!-- claim: rho_G_where_reilly_reaches_1 --> kg/m³. So the two closures fail in opposite
directions rather than being complementary, and above the SF₆ density **neither** has
support in this source.

What carries the transcription: every group in Eqs. 3 and 4 comes out exactly
dimensionless by exponent arithmetic, while both transition equations do not
(their constants carry units: feed them mPa·s and mN/m and Eq. 2 returns
2.95e-153 <!-- claim: eq2_value_with_mPa_s_inputs -->), and the source's own Notation list on p. 2634
declares σ in N/m, μ_L in Pa·s and ρ_G in kg/m³, so SI is the paper's stated
convention and not an inference. The correlation is also C¹ across its own regime
boundary, which is a structural identity, **not** a transcription check.

## Corrections

- 2026-11-03: said Eq. 4's own second term "dominates only at the high-velocity end".
  The second term is larger than the first at 42 <!-- claim: eq4_second_term_larger_points --> of the
  63 <!-- claim: n_points --> markers, above an excess of 0.072 <!-- claim: eq4_terms_equal_excess --> m/s.
- 2026-11-03: said the abscissa convention "is worth at most 15 % in `ε_b`", that this
  "closes `F1.4`'s open caveat: its +64 % is not an artefact" and that running end to end
  "multiplies" the two defects "giving +67 to +88 %". The 15 % holds at the six tabulated
  conditions (0.15 and 0.30 m/s, three gases), where the excess velocity is large. On the
  markers, which sit at small excess velocity, the transition velocity matters: `ε_b` is up to
  138 % <!-- claim: abscissa_effect_max_all_markers --> higher run end to end across the four gas densities, and the mean bias at air is
  +106 % <!-- claim: end_to_end_bias_air --> (+99 % <!-- claim: end_to_end_bias_air_fig6_U_df --> with the
  `U_df` read off Figure 6) against +64 % <!-- claim: wilkinson_eps_b_bias --> at the plotted abscissa. The caveat is
  not closed; the old +67 to +88 % mixed the air-density bias with ratios from other gases and velocities.
- 2026-11-03: said Table 3 makes Wilkinson "ten times worse than its replacement on `ε_b`".
  The ratio of Table 3's two deviations is 7.8 <!-- claim: table3_eps_b_wilkinson_over_replacement --> times.
- 2026-11-03: said the three points that give `ε_df ~ ρ_G^0.30` were read "off its
  Figure 6(b)". `ε_df` is panel (a) of Figure 6 and `U_df` is panel (b), as the caption
  states.
- 2026-11-03: the column-diameter table listed implied biases at three more diameters
  (0.174, 0.19 and 0.38 m, +81 %, +84 % and +108 %) that are in no dataset of the
  repository and had no recorded source; the text said "six" diameters while listing
  five. They are removed. The 0.63 m case, quoted from `F1.4`'s text as the largest column
  of Table 3's runs, stays: the implied bias there is +128 % <!-- claim: implied_bias_at_0p63m -->.
