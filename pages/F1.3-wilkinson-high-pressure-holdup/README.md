# F1.3 — Wilkinson high-pressure bubble-column gas holdup

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
| `index.ipynb` | the page; runs in about 4 s |
| `build_page.py` | regenerates `index.ipynb` from source. Edit this, not the notebook |
| `data/krishna-ellenberger-1996-table2-liquids.csv` | Table 2's five liquids, transcribed from a 600 dpi render |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Two further datasets are loaded **cross-page** from
`pages/F1.4-krishna-ellenberger-holdup/data/` and are not copied here: the 63
digitised Figure 11 markers and the parameters file holding the gas densities,
Reilly's *B* and the Table 3 deviations. Read
`krishna-ellenberger-1996-fig11.meta.yaml` on that page before using the `gas`
column — a maintainer review found the marker positions trustworthy and the
shape-derived gas labels not, and nothing on this page depends on series
identity.

## Regenerating

```bash
python build_page.py                    # rewrite index.ipynb
python ../../scripts/check_metadata.py  # then execute it via scripts/run_pages.py
```

## What this page is, next to F1.4

`F1.4` tests Krishna and Ellenberger's replacement correlation, Eq. 19, against
their Figure 11, and uses Wilkinson only as the foil — evaluating its
rise-velocity equations at the plotted abscissa and never touching **Eq. 2**, its
transition equation, at all. This page runs the correlation as a correlation,
from operating conditions, and Eq. 2 turns out to be where everything happens.
There is no overlap in the validation. The two pages share one deviation number,
Wilkinson's +63.9 %; this page retypes Eqs. 3 and 4 from the same source and
applies them to the same 63 rows, so that agreement is a **transcription
cross-check between the pages, not an independent computation**. This page also
closes the abscissa-convention caveat `F1.4` had to leave open.

## Headline result

**This is the authors' own finding, made specific.** Krishna and Ellenberger
write that "the Wilkinson correlation severely underpredicts the values of the
voidage and gas velocity through the dense phase", and compare Eq. 2 against
their measured `ε_df` themselves. At the gas densities they tested, Wilkinson's
Eq. 2 gives a transition holdup of 0.0028 for air–tetradecane against a measured
0.139, and 1.5 × 10⁻⁸ against 0.089 for helium — 1.7 decades low at air, 1.4 at
argon and 6.8 at helium. The homogeneous regime effectively does not exist in the
correlation at atmospheric density, so the small-bubble population is empty.

**The gap does not close at higher density.** Reilly's Eq. 8 exceeds Eq. 2 at
every density from 0.1 to 1000 kg/m³ — the ratio bottoms out at 2.85 near
29 kg/m³ and rises again, so the curves never cross. And the target is not fixed:
Eq. 2 crosses a *fixed* 0.15 at 14 kg/m³, but the paper states that increasing
gas density significantly increases the dense-phase voidage, and the three points
read off its Figure 6(b) scale as ρ_G^0.30, putting the measured level near
14 kg/m³ at about 0.32 — where Eq. 2 would still be a factor 2.1 low. **14 kg/m³
is a lower bound on where Eq. 2 could become adequate, not a crossover.**

That is one defect. The other is separate and the page separates them: at the
measured excess velocity Eq. 2 does not enter at all, and `ε_b` still comes out
64 % high on all 63 points, which is Eqs. 3 and 4 jointly — `V_b` is 39 % low,
and `V_b`'s first term *is* Eq. 3's `V_small`, 46 % of `V_b` on average (23–77 %,
and 77 % at the lowest excess velocity). Running the correlation end to end
multiplies the two, giving +67 to +88 %. The +64 % assumes all 63 points are air,
which the unreliable gas labels force; across the four gases it is +47 to +79 %,
a ±16 pp span. The two errors have opposite signs in the *total* holdup and about
two-thirds cancel, which is why Table 3 makes Wilkinson ten times worse than its
replacement on `ε_b` and level with it on `ε`.

Meanwhile Reilly's Eq. 8, which Krishna and Ellenberger recommend instead, passes
its own authors' stated data ceiling of 0.32 just past SF₆ at 6.7 kg/m³, turns
non-monotone at 15 kg/m³ and exceeds 1 at 70. So the two closures fail in
opposite directions rather than being complementary, and above 6.7 kg/m³
**neither** has support in this source.

What carries the transcription: every group in Eqs. 3 and 4 comes out exactly
dimensionless by exponent arithmetic, while both transition equations do not
(their constants carry units — feed them mPa·s and the answer is 10⁻¹⁵³), and the
source's own Notation list on p. 2634 declares σ in N/m, μ_L in Pa·s and ρ_G in
kg/m³, so SI is the paper's stated convention and not an inference. The
correlation is also C¹ across its own regime boundary, which is a structural
identity — **not** a transcription check. An earlier draft fitted the rate of
approach to that limit, "recovered" Eq. 4's printed 0.757 as 0.7566 and offered
it as one; the fit only reads back a constant typed into the function being
fitted, and the claim has been withdrawn.
