# E1.4: Sit and Grace (1981), bubble interaction and interphase mass transfer

Sit and Grace, *Chem. Eng. Sci.* **36**, 327-335 (1981),
[doi:10.1016/0009-2509(81)85012-9](https://doi.org/10.1016/0009-2509(81)85012-9), split the
single-bubble interphase coefficient of a fluidized bed into a throughflow term and a
penetration term (their eqns 1, 2), postulate that bubble interaction raises only the
throughflow term, and from that derive eqns (3), (4) for freely bubbling beds. The paper calls
its explanation "semi-quantitative". The page tests whether the paper's own printed numbers
(Table 2, the coefficients of interacting and isolated bubbles, Table 1) decide the postulate and whether eqns (3), (4)
follow from them. pymrm contributes two small solves (the penetration problem behind the
diffusive term, and the tracer balance of a growing bubble) and the checks; every quantity
also has a closed form.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `model.py` | eqns (1), (2) as two terms, potential-flow surface renewal, the pymrm penetration solve and the tracer march |
| `test_model.py` | the solves against their closed forms, two check values and one stored claim |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/sit-grace-1981-printed-numbers.csv` | every number of the text and equations the page uses, as printed |
| `data/sit-grace-1981-table1.csv` | Table 1, the tracer accounting of a vertical pair, as printed |
| `data/sit-grace-1981-table2.csv` | Table 2, the five studies, as printed |
| `data/sit-grace-1981-fig4.csv` | Fig. 4, the front bubble's nose height in the approach period, digitised |
| `data/sit-grace-1981-fig9.csv` | Fig. 9, tracer concentration of the front, back and single bubble, digitised |
| `data/check_values.csv` | the page's functions at test conditions, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. The page loads one claim of
[E1.2](../E1.2-davidson-bubble/), the through-flow of a two-dimensional Davidson bubble.

## Things to know before reusing this

- The throughflow term of eqns (1), (2) is Murray's, 0.5 <!-- claim: r1_md_2d --> of
  Davidson's in two dimensions; a model that uses Davidson's throughflow and adds this
  paper's enhancement applies it to a different base.
- Neither the voidage at minimum fluidization nor the bubble velocity is printed for any
  row; every statement about Table 2 is made over a bracket of both. Over every corner,
  throughflow-only enhancement orders 6 <!-- claim: r2_conc_T_min --> of
  6 <!-- claim: r2_n_pairs --> decidable pairs as printed. A free two-part fit leaves the
  diffusive enhancement undecided: with Pereira's row its excess is between
  -0.033 <!-- claim: r2_fit4_b_min --> and 0.040 <!-- claim: r2_fit4_b_max -->.
- At the paper's own conditions, eqn (1)'s throughflow term alone is
  1.0486 <!-- claim: r3_kT_over_iso --> times the isolated coefficient the paper prints,
  so eqn (1) cannot give it.
- Table 1 implies K = 11.34 <!-- claim: r4_K_front --> cm/s for the front bubble in the
  approach period under the page's reading of the transfer surface; two of its cells
  disagree with their products by 6.5 <!-- claim: r4_hu_0.6_rear --> and
  23.8 <!-- claim: r4_hu_1.1_resultant --> rounding allowances, and the page does not decide
  which entry is right.
- Eqn (3)'s constant follows from the paper's enhancement and coalescing fraction within
  0.40 <!-- claim: r5_hu_eqn3 --> half units of its last digit; those two inputs are not
  tested by any printed number, and nothing here tests eqns (3), (4) against a freely
  bubbling bed.
