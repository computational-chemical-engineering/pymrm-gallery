# A4.6: Chapman–Enskog binary diffusion

A binary diffusivity that is *derived* from the intermolecular potential rather
than fitted to diffusion data, and the question that makes it worth building:
how much of its accuracy is the kinetic theory, and how much is the
two-parameter fit underneath it.

**Source.** Chapman, S. and Cowling, T. G. (1970), *The Mathematical Theory of
Non-Uniform Gases: An Account of the Kinetic Theory of Viscosity, Thermal
Conduction and Diffusion in Gases*, 3rd edn, Cambridge University Press
(Cambridge Mathematical Library reissue 1990, reprinted 1993). This is the
**origin**, not a reprint: the primary account of the result the case is named
for, by one of the two people it is named after, so there is no
`reference_read_from`.

The PDF is a complete CCITT-G4 bilevel scan with an ABBYY FineReader text layer.
The prose layer is good; the numbers are not: the typesetting uses a raised
mid-dot decimal separator that the OCR drops, and leading zeros come back as the
letter "o". Every numeral on the page was read off a **cropped render at the
scan's native resolution**, twice.

## Files

- `index.py` / `index.ipynb`: the page (paired by jupytext).
- `model.py`: `ChapmanEnskog` (eqs. 14.2, 1 and 4; 14.21, 1-6), `CollisionTable`
  (Table 6 as splines), the 12,6 collision integrals by quadrature
  (`omega_11_22`), and `DiffusionCell`, the closed Loschmidt tube solved with
  pymrm.
- `test_model.py`: the Lorentz limit against the stored value, a turning point,
  a constant-D cell against its series and one composition-dependent cell.
- `data/chapman-cowling-1970-table22.csv`: the **measured D₁₂ at S.T.P.** with
  the book's own literature references (**tier 2**).
- `data/chapman-cowling-1970-table6.csv`: the collision integrals.
- `data/chapman-cowling-1970-table11.csv`: molecular weights, viscosities and
  viscosity diameters.
- `data/chapman-cowling-1970-table17.csv`: the Lennard-Jones force constants,
  in both of the book's two published versions.

**No figure is digitised anywhere on this page, and no dataset is borrowed from
another page.**

## The result

Over the 39 <!-- claim: d12_n_pairs --> pairs of Table 22 for which Table 17
supplies 12,6 constants for both gases, the first approximation (eq. 14.2, 4)
predicts the **measured** D₁₂ with a **mean absolute deviation of 3.65 <!-- claim: d12_mean_abs_dev_pct --> %**, a bias of **-2.74 <!-- claim: d12_bias_pct --> %** (standard error 0.59 <!-- claim: d12_bias_se_pct --> %), an r.m.s. of 4.57 <!-- claim: d12_rms_pct --> %
and a worst case of **-10.06 <!-- claim: d12_worst_dev_pct --> %** (N₂–CO₂);
29 <!-- claim: d12_within_5pct --> pairs deviate by less than 5 % and
38 <!-- claim: d12_within_10pct --> by less than 10 %. That is what §14.4's
"in general agree fairly well" is worth, and the book prints no such number
anywhere.

The comparison is **out of sample**: the force constants were fitted to pure-gas
*viscosity* (p. 235, explicitly), and no diffusion measurement enters the
prediction at any point. Table 23's constants, which *are* derived from D₁₂ data,
are deliberately not used.

**But most of that accuracy is the parameter fit, not the theory.** The book
prints a second set of force constants for the same molecules, obtained mainly
from virial coefficients, and warns that other authors' values "often differ
appreciably". Running the *identical* calculation on them gives **10.23 <!-- claim: d12_other_params_subset_mean_abs_pct --> %** where the viscosity
constants give **4.13 <!-- claim: d12_visc_params_subset_mean_abs_pct --> %** on
the same 27 <!-- claim: d12_subset_n_pairs --> pairs, 2.48 <!-- claim: d12_other_over_visc_ratio --> times worse; the paired difference of
the absolute deviations is 6.11 <!-- claim: d12_subset_diff_pp --> percentage
points with a standard error of 1.23 <!-- claim: d12_subset_diff_se_pp -->.
Nothing about the kinetic theory changed between those two numbers. So the page
reports the headline as what it is: *the 12,6 model transfers a two-parameter
fit from pure-gas viscosity to binary diffusion with a few per cent penalty, and
from equilibrium virial data with about ten per cent.*

The **null model** shows the collision integrals are nonetheless doing real work:
the same expression with 𝒲⁽¹⁾(1) deleted and Table 11's viscosity diameters gives
**17.24 <!-- claim: d12_rigid_sphere_mean_abs_pct --> %**, biased -17.24 <!-- claim: d12_rigid_sphere_bias_pct --> %, every pair under-predicted
(39 <!-- claim: d12_rigid_sphere_n_under --> of 39), 4.7 <!-- claim: d12_rigid_over_126_ratio --> times worse. That bias is predicted
independently by the book's own statement that the viscosity half-sum exceeds
the diffusion σ₁₂ "by about 10 per cent", reproduced here at 10.46 <!-- claim: sigma_ratio_mean_pct --> % over all 50 rows.

## Validation

**The headline is computed a second, independent way.** Ω⁽¹'¹⁾\* and Ω⁽²'²⁾\* are
recomputed from the 12,6 potential by quadrature of the classical deflection
angle, with no table consulted. It reproduces Table 6 to **0.18 <!-- claim: omega11_quad_vs_t6_max_pct --> %** (Ω⁽¹'¹⁾\*) and 0.30 <!-- claim: omega22_quad_vs_t6_half_max_pct --> % (Ω⁽²'²⁾\*) at worst over all
sixteen rows, both worst cases on the kT/ε = 100 row, where refining any of the
quadrature's four knobs moves it by at most 1e-3 <!-- claim: quad_refinement_at_worst_row_pct --> %: that residual is a difference
between two numerical integrations rather than an error in this one. At the
lowest rows the potential orbits and the production quadrature is not converged: at
kT/ε = 0.5 its gap to the table, 0.174 <!-- claim: omega_gap_at_T05_pct --> %, is its own
error (doubling n_E and n_b together moves it by 0.211 <!-- claim: quad_joint_change_at_T05_pct --> %);
refined, it agrees to 0.037 <!-- claim: omega_gap_at_T05_refined_pct --> %. At kT/ε = 0.3
the gap is 0.037 <!-- claim: omega_gap_at_T03_pct --> %, the same doubling moves it by
0.091 <!-- claim: quad_joint_change_at_T03_pct --> % and the refined gap is
0.097 <!-- claim: omega_gap_at_T03_refined_pct --> %. Substituting the quadrature for Table 6 entirely changes every
predicted D₁₂ by at most 0.1 <!-- claim: d12_two_routes_max_pct --> % and moves
the headline to 3.655 <!-- claim: d12_quadrature_route_mean_abs_pct --> %; on the
Table 6 rows that bracket the pairs the quadrature's own error is up to 0.048 <!-- claim: quad_refinement_pair_range_pct --> %, which limits how finely this
route tests the table. A break table could not have established this: it shows
the answer is *sensitive* to Table 6, never that Table 6 was read *correctly*.

**One reading had to be reconstructed.** Table 6's third column header carries
a fraction whose glyph is destroyed in the scan. It is **½, not ⅓**, established
twice from printed material: eq. (9.8, 7) with the rigid-sphere normalisation of
eq. (10.2, 1) (largest relative residual 2.9e-4 <!-- claim: t6_A_identity_max_rel -->, 1.00 <!-- claim: t6_A_identity_max_over_rounding --> of its rounding allowance), and
the independent quadrature above. Read as ⅓, both checks fail (rows of the
break table). It changes no diffusivity.

**Four printed-table round trips**, each against its rounding allowance (half a
unit of the last printed digit of every entry involved, carried through the
relation). Table 22's σ₁₂ from its D₁₂ through eq. (14.2, 1): worst 0.00509 <!-- claim: t22_sigma12_from_D_max_abs -->, inside its allowance
(at most 1 <!-- claim: t22_sigma12_max_over_rounding -->), *and this pins the
state*: at 1 atm the mean signed residual lies at most 2 <!-- claim: t22_sigma12_bias_z --> standard errors from zero in magnitude, while at 1 bar it is
34.8 <!-- claim: t22_sigma12_bias_z_1bar --> standard errors from zero with
every one of the 50 <!-- claim: t22_sigma12_1bar_positive --> residuals
positive, so "S.T.P." is established from the book's own arithmetic rather than
assumed. Table 22's half-sum from **Table 11's** σ: worst 0.0075 <!-- claim: t22_halfsum_from_t11_max_abs -->, which is 1.36 <!-- claim: t22_halfsum_max_over_rounding --> of its allowance, on the O₂–Xe row;
from Table 17's σ instead it is wrong by 0.740 <!-- claim: t22_halfsum_from_t17_max_abs -->. Table 11's σ from its μ and M
through eq. (12.1, 6): largest relative residual 8.7e-4 <!-- claim: t11_sigma_from_mu_max_rel -->, rising to 8.8e-3 <!-- claim: t11_sigma_from_mu_no_corr_max_rel --> without the 1.016 correction;
the largest, D₂'s, is 2.07 <!-- claim: t11_sigma_from_mu_max_over_rounding --> of
its allowance. Table 17's rₘ against 2^(1/6)σ: worst
0.0052 <!-- claim: t17_rm_max_abs -->, inside its allowance (at most 1 <!-- claim: t17_rm_max_over_rounding -->). Plus Table 6's C column against
eq. (14.4, 1): r.m.s. at most 0.01 <!-- claim: t6_C_vs_dlnW_rms --> on the
log-log slope.

**Four printed numerical results** of the second approximation, with the limits
taken in closed form: eq. (14.3, 2) at five mass ratios to at most 1e-12 <!-- claim: eq1432_max_rel -->; [D₁₂]₂/[D₁₂]₁ = 1.083333 <!-- claim: lorentz_D2_over_D1 --> against the printed "1.083"; 8.3333 <!-- claim: lorentz_variation_pct --> % against "8⅓ per cent"; and Kihara's B = ¾
giving 11.1111 <!-- claim: kihara_variation_pct --> % against "1/9, or 11·1 per
cent". Section 14.3's bound "unlikely to exceed about 6 per cent" holds: the
largest end-to-end swing over the pairs is 3.28 <!-- claim: composition_swing_max_pct --> %, and the largest variation over the
whole composition range 3.39 <!-- claim: composition_variation_max_pct --> %
(the same He–Kr ends, measured on the light end's value); for all 37 <!-- claim: composition_heavier_rises --> pairs of unequal mass D₁₂ rises with the
concentration of the heavier molecule.

**Two "wrong" choices score better than the published ones**, and the page says
so: a geometric-mean σ₁₂ gives 3.26 <!-- claim: d12_mean_abs_geometric_sigma_pct --> %
and a 1 bar reference 3.10 <!-- claim: d12_mean_abs_1bar_pct --> %. Neither is
adopted: eq. (14.2, 1) prints the arithmetic mean, and the σ₁₂ round trip pins
1 atm. The headline therefore sits on a shallow optimum and is quoted to one
decimal at most.

## What pymrm adds

The closed-cell solver: `construct_grad` and `construct_div` for the zero-flux
tube, `construct_boundary_value_matrices` for the boundary face compositions,
`construct_coefficient_matrix` for the face diffusivity and its derivative in
the exact Jacobian, `newton` per backward-Euler step with every step checked,
`pymrm.checks.find_roots` for the inversion of the constant-D series and
`pymrm.checks.observed_orders` for the refinement ladders. The Chapman–Enskog
algebra and the quadrature are the page's own.

§14.32 states that "the variation of D₁₂ with the proportions of a mixture is
ignored" when the data are discussed, but the apparatus behind many Table 22
entries is a closed tube whose composition runs from 0 to 1 along its own length.
Solving that cell with D₁₂ = [D₁₂]₂(x₁) and inverting the constant-D closed form
for the apparent D, exactly as an experimenter would, gives **+2.53 <!-- claim: cell_apparent_excess_pct_HeXe --> %** above [D₁₂]₁ for He–Xe,
+1.80 <!-- claim: cell_apparent_excess_pct_H2N2 --> % for H₂–N₂ and **+0.31 <!-- claim: cell_apparent_excess_pct_N2CO --> % for the equal-mass control N₂–CO**,
each extrapolated to zero cell size and zero step (the production level is
within 0.01 <!-- claim: cell_excess_discretisation_pp --> percentage points; the
observed orders lie of at least 1.9 <!-- claim: cell_grid_order --> and at most 2.1 <!-- claim: cell_grid_order --> in cell
size and of at least 0.9 <!-- claim: cell_time_order --> and at most 1.1 <!-- claim: cell_time_order --> in step). It is
*not* the equimolar second approximation (+2.68 <!-- claim: cell_equimolar_excess_pct_HeXe --> %, +1.89 <!-- claim: cell_equimolar_excess_pct_H2N2 --> %, +0.31 <!-- claim: cell_equimolar_excess_pct_N2CO --> %): the cell falls short by 0.15 <!-- claim: cell_shortfall_pp_HeXe --> and 0.09 <!-- claim: cell_shortfall_pp_H2N2 --> percentage points for the first two, the
shortfall growing with the mass ratio at both steps
(2 <!-- claim: cell_shortfall_falls_with_mass_ratio --> of 2) and vanishing for the
control. A transient cell reports a weighted average over a composition
history, not a value at a stated composition. The page states plainly that this
is a model of the apparatus, not the apparatus, and is not a correction to be
applied to Table 22.

## Scope

Binary, dilute gases at S.T.P. The Fuller–Schettler–Giddings **empirical**
correlation for the same quantity is `A4.5`; the page argues that the joint
comparison (one measured-D₁₂ axis with a first-principles prediction and an
empirical correlation overlaid) should be built **there**, since `A4.5` must
transcribe its own source anyway and can load this page's Table 22 cross-page.
Multicomponent diffusion is `A4.2`/`A4.9`, pore transport `A4.3`/`A4.4`.

**This monograph does not source `A4.1`.** Its only Wilke is Buddenberg & Wilke
(1949) on the *viscosity* of a mixture, cited once in §12.43; the 1950 diffusion
mixture rule is not stated, named or attributed anywhere.

## Corrections

- 2026-11-16: "worst 0.0075 x 1e-8 cm, consistent with rounding" said of Table
  22's half-sum rebuilt from Table 11's σ. On the O₂–Xe row the residual is
  1.36 <!-- claim: t22_halfsum_max_over_rounding --> of its rounding allowance
  (half a unit of the printed half-sum plus half the half units of the two
  printed σ), so that row is not explained by rounding; the exact half-sum, 4.2875 <!-- claim: t22_halfsum_worst_exact -->, rounds to 4.29, not the printed 4.28 <!-- claim: t22_halfsum_worst_printed -->.

## Rebuilding

```bash
jupytext --sync index.py                 # keep index.py and index.ipynb paired
python ../../scripts/run_kernels.py .     # execute under the three BLAS kernels
```
