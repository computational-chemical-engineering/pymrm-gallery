# A2.8 — Zwietering's micromixing bounds

**Catalog ID:** `A2.8` · **Section:** A · **Tier:** T0 · **Priority:** P2 ·
**Structures:** `S1`, `S3`, `S5` · **Data tier:** 6 (no measurement exists)

Reproduces Th. N. Zwietering, *The degree of mixing in continuous flow systems*,
**Chemical Engineering Science 11(1) 1–15 (1959)**,
[doi:10.1016/0009-2509(59)80068-3](https://doi.org/10.1016/0009-2509(59)80068-3).

> **The page range is not 1–11.** The article is cited that way essentially
> everywhere. Its last page carries the printed folio **15**, holds eqs. (II, 11)
> and (II, 12) and the three-item reference list, and the printed folios run 1 to
> 15 without a gap. Read off a 300 dpi render, not from metadata.

## What the page shows

For a reaction of any order other than one, knowing the residence-time
distribution is **not enough** to fix the conversion. Zwietering's two extremes —
Danckwerts' complete segregation (mixing as late as possible) and his own state of
maximum mixedness (as early as possible) — bracket every reactor with that RTD.

The page reproduces all 44 conversions and all 8 degrees of segregation he
printed, with **zero fitted parameters**, and then measures the thing the bounds
exist for: how wide the bracket is. It is wide, and it widens with reaction rate —
a factor **2.12** between the two bounds at $kc_0\tau = 50$ for a two-tank RTD.
Both of Zwietering's *intermediate* reactors sit well inside the band, at 0.23 and
0.70–0.85 of its width, so the bracketing check is a check that can fail.

## Two findings about the paper

1. **The conversion bodies of Tables 1 and 2 are interchanged with respect to
   their captions.** Proved from three closed forms Zwietering printed himself —
   eq. (42), derived from eq. (41)'s $2/\tau$ and therefore unambiguously the
   two-vessel chain, reproduces Table 2's row (captioned *three* vessels) to
   4.1e-04 and Table 1's only to 3.3e-02; eqs. (38) and (39)+(40) do the same for
   the other two algebraic rows. Scored over all 44 values the interchanged
   reading gives a worst deviation of 4.5e-03 against 3.5e-02. The $J$ column does
   **not** move with the body, which is what makes it an interchange rather than a
   relabelling.
2. **One printed value is not reproduced, and is not repaired.** The three-tank
   degree of segregation under maximum mixedness computes to 0.080296 against the
   printed 0.0831. The two-tank companion reproduces to four decimals and
   Appendix I's variance identity closes to 6.1e-16. Reported as unresolved.

Also measured: Zwietering's own hand-integration error. Every *algebraic* row in
both tables is inside his three-decimal rounding; the only row outside it is the
maximum-mixedness row he had to integrate by hand for the two-tank RTD, which
carries a **+1.10 %** bias with all four printed deviations of one sign.

## What pymrm contributes

Eq. (31) is normally solved as an initial-value problem integrated backwards from
infinite life expectancy. Written conservatively,

```
div(v gamma) = E(x) - K w(x) gamma^p ,    v = -w(x) ,    x = lambda/tau
```

it becomes an ordinary convection–reaction boundary-value problem with a decaying
velocity field and a distributed side feed — `construct_convflux_upwind` plus
`construct_div` plus a pointwise `NumJac` source — so the life-expectation
coordinate behaves like any other spatial axis and the van Leer deferred
correction applies unchanged at second order.

The four em dashes in Zwietering's maximum-mixedness rows are filled in. He
integrated eq. (31) only for $K = 5, 10, 20, 30$.

## Validation summary

| Check | Result |
|---|---|
| 44 printed conversions, nothing fitted | worst 4.5e-03, rms 8.3e-04, 40/44 inside 0.001 |
| the same, read as captioned | worst 3.5e-02, 0/44 inside 0.001 |
| pymrm finite volume vs gridless DOP853 integration | 8.6e-07 over 12 cases |
| grid refinement | order 2.00 (van Leer), 1.00 (bare upwind) |
| domain truncation | decay 8.44 per unit X; 4.1e-15 at production X = 6 |
| exponential RTD → stirred-tank root (his eq. II, 12) | 3.5e-14 |
| first order → both bounds collapse onto (1+K/n)^-n | 3.5e-14, two independent code paths |
| plug-flow limit (300 tanks) | bracket 4.0 % wide, down from 83 % at 2 tanks |
| bracket reverses through order one | +8.3e-03 at p = 1.1, −7.4e-03 at p = 0.9 |
| 8 printed degrees of segregation | 7 to 4.3e-05 (four exact rationals); 1 unresolved |
| defect injection, 15 rows | all physics defects ≥ 2.5e-02, all discretisation defects ≥ 3.1e-05 |

**Labelled structural:** Appendix I's variance identity (6.1e-16) — his appendix
proves it must hold, and the two sides share `alpha_P`, so it tests the mean
consistency of `alpha_P` and not its shape.

**Labelled untested:** the far-field boundary condition. Replacing the root by
`gamma = 1` moves the answer by 5.9e-07 — the grid error, i.e. nothing. Shorten
the domain to X = 0.5 and the same defect is worth 2.1e-02, which is what shows
the insensitivity belongs to the domain length rather than to the check.

**Outside CI:** six metrics fall below `check_agreement.py`'s `ABS_FLOOR` of
1e-12 and are named on the page as unprotected by the regression suite.

## Files

```
index.ipynb        generated by build_page.py - edit the builder, not the notebook
build_page.py      prose and code in one reviewable file
meta.yaml          page metadata
agreement.json     39 metrics, written by report_agreement
data/zwietering-1959-conversions.csv          44 printed values of ce/c0
data/zwietering-1959-segregation-degree.csv    8 printed values of J
```

Runtime about 37 s.

## Related pages

[`A2.1`](../A2.1-danckwerts-boundary-conditions/) (Danckwerts boundary
conditions), [`A2.3`](../A2.3-taylor-aris-dispersion/) (Taylor–Aris),
[`A2.5`](../A2.5-edwards-richardson-dispersion/) (measured dispersion),
[`A2.6`](../A2.6-gunn-dispersion-correlations/) (dispersion correlations). All
four end with an RTD; this page is about what an RTD still leaves undetermined.
`A2.4` (tanks in series) is **not** built, so the Erlang distributions used here
have no page of their own — they are analytic and verified in place. No sibling
page's dataset is loaded.
