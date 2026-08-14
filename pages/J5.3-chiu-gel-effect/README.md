# J5.3 — The gel effect without a switch

**A model paper with printed constants, no printed data, and two claims about
itself that can be checked without any data at all.**

Chiu, Carratt and Soong (1983) describe the Trommsdorff–Norrish gel effect with
one unbranched constitutive equation,

```
1/k_t = 1/k_t^0 + theta_t * P / exp[2.3 phi_m/(A + B phi_m)]        (eq. 31)
```

— a reaction resistance and a mass-transfer resistance in series — and they
criticise the alternatives on book p. 351 for switching diffusion control on at
a prescribed critical conversion *"in a somewhat ad hoc fashion"*.

**The experimental comparison lives only in Figures 3–9**, whose markers are
Marten & Hamielec's data, so that half of the case is out of scope and **nothing
here is digitised**. What the paper prints in numerals is two tables and six
stated results about its own output, and that is what this page tests. **The
page therefore establishes nothing about whether the model describes PMMA**, and
says so wherever it matters.

## What it finds

**The paper's own onset construction, computed instead of read.** Book p. 354
defines the onset as the intersection of a line extrapolated from *"the sloped
region"* of the `log k_t` curve with *"the initial value"*, and states 0.26,
0.35 and 0.45 at 50, 70 and 90 °C. Reconstructed as a least-squares line over
the decades `L_0-2` to `L_0-4`, the model gives **0.25516, 0.34838 and
0.44034** — low by 0.00484, 0.00162 and **0.00966** — with the increase with
temperature reproduced.

**Those are at `I_0 = 0.0258` mol/L and the page says so.** The paper quotes one
onset per temperature while Figures 10–12 plot both loadings and the sentence
names neither; 0.0258 is the loading listed first in all three legends and the
one Figures 16–18 plot, and that is the whole basis for the choice. **The same
construction on the 0.01548 curve of the same figures misses by up to
0.03126**, 3.24× the headline miss.

**The limit on that number is the sentence, not the arithmetic**: seven defensible
readings — six of *"the sloped region"*, one of *"the initial value"* — spread
the answer by **0.09009**, **9.33×** the largest miss, while the resampling grid
moves it by **1.311e-05** relative (so five decimals are quoted, not six) and a
completely independent integration reproduces it to **9.08e-12** — a number that
prices the attracting manifold `lambda_0(x)`, not the accuracy of either
integration. The band is reported beside the value, and the page says what the
band does **not** cover: the steepest-slope tangent, the most literal reading of
all, gives 0.44783, 0.47680 and 0.51866, and is excluded because its tangent
point sits at the right edge of *this page's* conversion window at all six
conditions. That is the only ground on which anything is kept out, and it is
applied to what is kept in: the widest-swinging reading in the band, the
4-6-decade window (0.29957, 0.41999, 0.49516), passes the same test — its fit
mask ends inside `x = 0.86` at four of the six conditions, including the one at
which the band is widest.

**The paper's own correlations, recomputed from its own Table II.** Figure 13's
Arrhenius plots give **34.908161** and **34.516600** kcal/mol against the printed
**34** (+2.67 % and +1.52 %) and **27.771564** against the printed **28**
(−0.82 %). Figure 15's *"All three points fall on the same straight line"* holds
to **8.819e-04** relative — against **1.535e-02** for the same three points fitted
linearly in `T`, **17.4057×** worse, with neither form nested inside the other,
so the direction of that ratio is evidence. The `T_gp` that would make the three
`A` values exactly collinear is **111.42857 °C** against Table I's 114, but the
printed rounding of `A` alone admits **[102.2222, 128.0000] °C** — the three cells
cannot tell them apart, and the page reports the interval.

**What that ratio is evidence *for* is concavity in `T`, not the quadratic
form**, and the page says so with its own numbers: the gain is 3.6124× and
3.5347× at the two ends of that admissible `T_gp` interval and diverges inside
it, and two other non-nested two-parameter regressors also beat linear-in-`T`
(`A` against `1/T`, 5.5926×; `A` against `(T − 130)^2`, 3.2607×). With three
points and one residual degree of freedom, what all of them read is the second
difference of `A`, **−0.007000**.

**"Only a small number of adjustable parameters", counted and then exercised.**
Four named parameters; **13 distinct fitted numbers** over **6 conditions**. The
paper's own Figures 13 and 15 compress those 13 to **8**, reproducing `A` to
0.0882 % and `theta_p` to 8.3695 % but `theta_t` only to **18.8281 %** — which
moves the onset conversion by **0.01633**, **1.691×** the distance by which this
page's reconstruction misses the printed onsets. The compression is true as a
procedure and lossy as an identity.

**"In a somewhat ad hoc fashion", tested rather than repeated.** Eq. 31 and 32
carry no branch anywhere. The nearest thing to an onset knob is `theta_t`, and it
moves the onset by **−0.21735 per decade** — it takes **4.601 decades** to buy
one unit of onset, where a prescribed critical conversion buys it 1:1 by
construction. The counterfactual settles the paper's own argument: the same model
switched on at `x_c = 0.3` is concave **down** over 100 % of its pre-switch
portion and jumps **0.760194** of a decade in `log k_t` at the switch, while the
published model is concave **up** from `x = 0.30643` to the sharp rise, every
onset inside that window by at least **0.12933**. **One qualification the paper
does not make:** below that lower bound the published model's curve is concave
*down* too — initiator depletion — so the sentence is true of the region it means
and not of the whole pre-gel curve.

**The quasi-steady-state assumption, at six conditions instead of one.** The
paper's *"more than a factor of 2"* is **3.5176×** at the condition it plots, with
the conversion history still agreeing to **0.038611** in `x`. Across all six the
factor runs 1.4889 to 4.8466 and **exceeds 2 at three and fails at three**.

**Five identities proved rather than asserted.** Eq. 24 and eq. 26 re-derived
symbolically from eq. 21 (exactly zero); eq. 22's *"Since r_D ≫ r_m"* shown to
cost **exactly r_m/r_D**, symbolically and then a second time through the flux a
pymrm spherical-shell BVP actually computes — **1.858e-05** over five ratios at
`n = 6400`, converging onto `r_m/r_D` at observed order **1.8953**; the moment
equations 12–17 reproduced from eq. 4a, 4b and 5 summed over a
seeded population to **2.10e-16** — including the `k_tc` terms, which Table I sets
to zero and which therefore **no condition the paper reports ever exercises**.

**One retraction, about what a number *was*.** Until 2026-08-14 this page priced
eq. 22 as `1 − K_22/K_exact` and reported the result, 1.53e-16, as a pymrm BVP
result and as one of its independent routes. It was neither: both factors are
closed-form scalars computed inside the shell function, so the expression is the
identity `1 − r_m(1/r_m − 1/r_D) = r_m/r_D` in floating point, and it returns the
same bits with the geometry set to Cartesian, with two cells, and with no solve
at all. The notebook now demonstrates that beside the corrected number.

**And a structural point the paper's prose does not make.** Its *"limiting
conversion at long times"* is reproduced as an **ordering** (0.95535, 0.98173,
0.99289 at a stated rate threshold), but this model has **no limiting conversion
in the mathematical sense**: as `x → 1` the Fujita–Doolittle factor returns to 1,
`1/k_p → 1/k_p^0 + theta_p lambda_0` stays finite, and `dx/dt` vanishes only
through `(1−x)`. Integrated far enough every condition creeps to 0.97897, 0.99850
and 0.99998 and is still rising.

## Printed features, reported and not repaired

- **Figure 3's legend prints `0.01584`** where Table II and the legends of
  Figures 4, 5 and 8 and the captions of Figures 6 and 8 print `0.01548` — the
  last two digits transposed. Five printings against one; the CSV carries
  0.01548 and the odd one out is recorded.
- **Table I gives the first-order `k_d` the units `L/min`**, where eq. 1 makes
  them `1/min`. The CSV carries the printed string.
- **Book p. 353 reads *"The presence model"*** [sic].
- Eq. 31 and 32 print the base-ten conversion as **2.3**, not ln(10). The model
  is run with the printed 2.3 and the alternative is priced: 0.00116 in onset
  conversion.

## The file straddles three articles

PDF page 1 opens with the References and Notes of the **preceding** article (a
Zambelli stereochemistry paper); PDF page 10 opens the **next** one
(Macromolecules 1983, 16, 357–359, Leonard O. Moore), which has a **Table I of
its own** about the chain length of normal alkanes. Every numeral used here was
located on a page whose running head names *this* article before it was read, and
then read on a 300 ppi crop at digit scale.

## Files

- `index.ipynb` — the page (runtime ≈ 190 s)
- `build_page.py` — regenerates it
- `data/chiu1983-table1-rate-constants.csv` — Table I, book p. 352
- `data/chiu1983-table2-model-parameters.csv` — Table II, book p. 353
- `data/chiu1983-stated-results.csv` — the six stated numerical results, each
  with the sentence it was read from
