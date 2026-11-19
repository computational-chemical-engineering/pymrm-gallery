# J5.3: The gel effect without a switch

**A model paper with printed constants, no printed data, and two claims about
itself that can be checked without any data at all.**

Chiu, Carratt and Soong (1983) describe the Trommsdorff-Norrish gel effect with
one unbranched constitutive equation,

```
1/k_t = 1/k_t^0 + theta_t * P / exp[2.3 phi_m/(A + B phi_m)]        (eq. 31)
```

a reaction resistance and a mass-transfer resistance in series, and they
criticise the alternatives on book p. 351 for switching diffusion control on at
a prescribed critical conversion *"in a somewhat ad hoc fashion"*.

**The experimental comparison lives only in Figures 3 to 9**, whose markers are
Marten and Hamielec's data, so that half of the case is out of scope and
**nothing here is digitised**. What the paper prints in numerals is two tables
and its stated results about its own output, and that is what this page tests.
**The page therefore establishes nothing about whether the model describes
PMMA**, and says so wherever it matters.

## What it finds

**The paper's own onset construction, computed instead of read.** Book p. 354
defines the onset as the intersection of a line extrapolated from *"the sloped
region"* of the `log k_t` curve with *"the initial value"*, and states one onset
conversion at each of 50, 70 and 90 °C (`data/chiu1983-stated-results.csv`).
Reconstructed as a least-squares line over the decades `L_0-2` to `L_0-4`, the
model gives 0.25516 <!-- claim: onset_x_50C -->, 0.34838 <!-- claim: onset_x_70C --> and 0.44034 <!-- claim: onset_x_90C -->, off
the printed values by -0.00484 <!-- claim: onset_dev_50C -->, -0.00162 <!-- claim: onset_dev_70C --> and -0.00966 <!-- claim: onset_dev_90C -->,
with the increase with temperature reproduced.

**Those are at the higher initiator loading of Table II and the page says so.**
The paper quotes one onset per temperature while Figures 10 to 12 plot both
loadings and the sentence names neither. **The same construction on the other
loading of the same figures misses by up to** 0.03126 <!-- claim: onset_dev_other_loading_max -->.

**The limit on that number is the sentence, not the arithmetic**: the
defensible readings of *"the sloped region"* and *"the initial value"* that the
notebook lists spread the answer by 0.09009 <!-- claim: onset_definition_band_max -->, while the resampling grid moves it by
1.311e-05 <!-- claim: onset_resample_floor_max_rel --> relative (so five
decimals are quoted, not six) and a completely independent integration
reproduces it to below 1e-10 <!-- claim: onset_two_route_max_rel -->
relative, a number that prices the attracting manifold `lambda_0(x)`, not the
accuracy of either integration. The notebook prints the one reading the band
does not cover, the steepest-slope tangent, and the measured reason it is kept
out: its tangent point sits at the right edge of this page's conversion window
at every condition.

**The paper's own correlations, recomputed from its own Table II.** Figure 13's
Arrhenius plots give 34.908161 <!-- claim: dE_theta_t_0258_kcal --> and
34.516600 <!-- claim: dE_theta_t_01548_kcal --> kcal/mol for `theta_t`, the first
2.67 % <!-- claim: dE_theta_t_rel_vs_printed --> from the printed value, and
27.771564 <!-- claim: dE_theta_p_kcal --> for `theta_p`, 0.82 % <!-- claim: dE_theta_p_rel_vs_printed --> from the printed value. Figure 15's
*"All three points fall on the same straight line"* holds to 8.819e-04 <!-- claim: A_linear_in_Tgp2_max_rel --> relative, against 1.535e-02 <!-- claim: A_null_linear_in_T_max_rel --> for the same three points fitted
linearly in `T`, 17.4057 <!-- claim: A_null_gain --> times worse, with
neither form nested inside the other. The `T_gp` that would make the three `A`
values exactly collinear is 111.42857 <!-- claim: Tgp_implied_C --> °C,
inside the interval the printed rounding of `A` admits, so the three cells
cannot tell it from Table I's value.

**What that ratio is evidence *for* is concavity in `T`, not the quadratic
form**: the gain falls to 3.5347 <!-- claim: A_gain_min_at_Tgp_band_ends --> at
an end of that admissible `T_gp` interval, other non-nested two-parameter
regressors beat linear-in-`T` as well (the notebook prints them), and with three
points and one residual degree of freedom what all of them read is the second
difference of `A`, -0.007000 <!-- claim: A_second_difference -->.

**"Only a small number of adjustable parameters", counted and then exercised.**
The notebook counts the named parameters and the distinct fitted numbers of
Table II. The paper's own Figures 13 and 15 compress them, reproducing `A` to
0.0882 % <!-- claim: compress_A_max_rel --> and `theta_p` to 8.3695 % <!-- claim: compress_theta_p_max_rel --> but `theta_t` only to 18.8281 % <!-- claim: compress_theta_t_max_rel -->, which moves the onset conversion by
0.01633 <!-- claim: onset_compressed_shift_max -->. The compression is true
as a procedure and lossy as an identity.

**"In a somewhat ad hoc fashion", tested rather than repeated.** Eq. 31 and 32
carry no branch anywhere. The nearest thing to an onset knob is `theta_t`, and
it moves the onset by -0.21735 <!-- claim: onset_per_decade_theta_t --> per
decade, where a prescribed critical conversion moves it one for one by
construction. The counterfactual settles the paper's own argument: the same
model switched on at a critical conversion is concave **down** over its
pre-switch portion (largest second derivative -1.391e-06 <!-- claim: switch_d2x_max_pre -->) and jumps 0.760194 <!-- claim: switch_logkt_jump --> of a decade in `log k_t` at the switch, while
the published model is concave **up** from `x = 0.30643` <!-- claim: concave_up_lower_x_max --> to the sharp rise, every onset inside
that window, the smallest margin 0.129 <!-- claim: concave_onset_margin_min -->. **One qualification the paper does
not make:** below that lower bound the published model's curve is concave *down*
too (initiator depletion), so the sentence is true of the region it means and
not of the whole pre-gel curve.

**The quasi-steady-state assumption, at every condition instead of one.** The
paper's *"more than a factor of 2"* is 3.5176 <!-- claim: qssa_lam0_ratio_50C -->
at the condition it plots, with the conversion history still agreeing to
0.038611 <!-- claim: qssa_conv_max_dev_50C --> in `x`; the notebook prints
the factor at the other conditions, where it does not always exceed 2.

**Five identities proved rather than asserted.** Eq. 24 and eq. 26 re-derived
symbolically from eq. 21 (exactly zero); eq. 22's *"Since r_D >> r_m"* shown to
cost **exactly r_m/r_D**, symbolically and then a second time through the flux a
pymrm spherical-shell BVP actually computes, within 1.858e-05 <!-- claim: eq22_rm_over_rD_max_dev --> over five ratios on the finest shell;
the moment equations 12 to 17 reproduced from eq. 4a, 4b and 5 summed over a
seeded population to a relative residual below 1e-14 <!-- claim: moment_identity_ktc_max_rel -->, including the `k_tc` terms, which
Table I sets to zero and which therefore **no condition the paper reports ever
exercises**.

**And a structural point the paper's prose does not make.** Its *"limiting
conversion at long times"* is reproduced as an **ordering** (0.95535 <!-- claim: xlim_50C --> at 50 °C, rising with temperature, the smallest step 0.011165 <!-- claim: xlim_monotone_margin --> at a stated rate threshold), but this model
has **no limiting conversion in the mathematical sense**: as `x -> 1` the
Fujita-Doolittle factor returns to 1, `1/k_p -> 1/k_p^0 + theta_p lambda_0`
stays finite, and `dx/dt` vanishes only through `(1-x)`. Integrated far enough
every condition creeps on towards full conversion and is still rising.

## Printed features, reported and not repaired

- **Figure 3's legend prints `0.01584`** where Table II and the legends of
  Figures 4, 5 and 8 and the captions of Figures 6 and 8 print `0.01548`, the
  last two digits transposed. Five printings against one; the CSV carries
  `0.01548` and the odd one out is recorded.
- **Table I gives the first-order `k_d` the units `L/min`**, where eq. 1 makes
  them `1/min`. The CSV carries the printed string.
- **Book p. 353 reads *"The presence model"*** [sic].
- Eq. 31 and 32 print the base-ten conversion as **`2.3`**, not ln(10). The model
  is run with the printed `2.3` and the alternative is priced: 0.00116 <!-- claim: ln10_onset_shift --> in onset conversion.

## The file straddles three articles

PDF page 1 opens with the References and Notes of the **preceding** article (a
Zambelli stereochemistry paper); PDF page 10 opens the **next** one
(Macromolecules 1983, 16, 357-359, Leonard O. Moore), which has a **Table I of
its own** about the chain length of normal alkanes. Every numeral used here was
located on a page whose running head names *this* article before it was read, and
then read on a 300 ppi crop at digit scale.

## Files

- `index.py`, paired with `index.ipynb`: the page
- `model.py`: the model (eq. 10-17 with 29-32, the backward-Euler marcher, the
  shell BVP, the moment identity and the onset construction); `test_model.py`
  checks it
- `data/chiu1983-table1-rate-constants.csv`: Table I, book p. 352
- `data/chiu1983-table2-model-parameters.csv`: Table II, book p. 353
- `data/chiu1983-stated-results.csv`: the stated numerical results, each with
  the sentence it was read from

## Corrections

- 2026-11-19: meta.yaml called the seed of the moment-identity population "the volume and first page of the article". It is the year and the first page; the volume is 16.
- 2026-11-19: the Reuse section said "Three things on this page" generalise and then listed four. It now says four.
- 2026-11-19: the break table labelled the coarse shell row "BVP refined only to n = 80 instead of n = 640". The baseline ladder ends at n = 1280.
