# D2.2 — Parametric sensitivity and runaway in a fixed bed

Two 1970 criteria for when a cooled tubular reactor stops being controllable,
rebuilt and then swept over the whole operating plane the original could only
sample.

- **Structures:** `S2` (plug flow with reaction)
- **Reference:** Van Welsenaere & Froment (1970), Chem Eng Sci 25(10) 1503–1516,
  doi:10.1016/0009-2509(70)85073-4
- **Runtime:** ~10 s

## Agreement

**0.054 % mean absolute deviation over all 30 numbers printed in their
Section 6**, worst 0.40 %, nothing fitted. Split by what they test: 0.051 % over
the 25 closed-form criterion values, 0.069 % over the 5 numerical-integration
values.

Separately, the pymrm reactor solve and an independent phase-plane
back-integration agree on the critical inlet pressure to **0.18 %** over
600–700 K. Those two numbers are the only pair on the page that does not involve
the paper.

## Provenance: tier 6, not experimental

The paper contains **no experimental data of any kind.** Every value here is the
authors' own computed result — partly closed-form extrapolation, partly their
fourth-order Runge–Kutta. The page is validated against a published reference
solution in the same sense as `B1.1` and `F3.1`. Do not describe it as validated
against experiment.

## Two traps

1. **`c_p` is printed as kcal/m³·°C — it is already volumetric.** Multiplying by
   `rho_g` "to fix the units" changes `C = 2U/(c_p R)` by a factor of 1.293 and
   moves the runaway boundary with no other symptom. `rho_g` appears only in `A`.
2. **The Elsevier full-text API is useless for the numbers in this paper.** It
   returns the publisher's OCR of the 1970 scan, which discards the mid-dot
   decimal separator: `R = 00125 m` for 0.0125 m, `b = 19837` for 19.837,
   `001353` for 0.01353 atm. Excellent for prose, unusable for parameters.
   Everything numeric here was read from a 600 dpi page render.

The check that the reading is right costs nothing and is in the notebook: four
quantities the paper prints (ln K, t_w, T_M and the Example 1 bracket) each
depend on the whole parameter table, and none was used to obtain it.

## What the page adds

Their Figs. 6 and 7 are graphical solutions of one implicit equation (Eq. 20),
included explicitly so no computer was needed. Solving it directly turns both
into checks. Their Fig. 8, drawn continuously rather than at a handful of
back-integrated points, shows that the second criterion has the narrower bracket
and the more conservative critical value — as claimed — but that the **first**
criterion's midpoint is the better estimate. The paper compares them only at
625 K, where the two are indistinguishable.

## Rebuilding the page

```bash
python build_page.py                    # regenerate index.ipynb from the builder
python ../../scripts/run_pages.py       # execute it
python ../../scripts/check_agreement.py # metric regression check
```

No PDF is needed to rebuild: both datasets are committed CSVs.
