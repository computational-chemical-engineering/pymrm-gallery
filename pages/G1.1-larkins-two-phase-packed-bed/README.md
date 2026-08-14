# G1.1 — Larkins two-phase concurrent flow in packed beds

**Source.** Larkins, R. P., White, R. R. and Jeffrey, D. W. (1961),
*Two-Phase Concurrent Flow in Packed Beds*, AIChE Journal **7**(2) 231–239,
[doi:10.1002/aic.690070213](https://doi.org/10.1002/aic.690070213). Read on
300 ppi renders — the scan's native resolution. Book page = PDF page + 230.
**PDF page 1 opens with the previous article's nomenclature**; Larkins et al.
start below the rule on the same page.

## Why this paper got a page

It prints a **worked example**. Table 3 (book p. 234) is eight runs of measured
two-phase data on ⅜-in. Raschig rings; Table 4 (book p. 235) is the authors'
own arithmetic on the identically numbered rows. Every constant that turns one
into the other is printed somewhere in the nine pages, so the chain can be
redone cell by cell — which is the highest validation rank the builder brief
offers short of held-out data.

## What the page establishes

| | |
|---|---|
| liquid Reynolds column reproduced | 0.174 % worst, 0.147 % r.m.s. |
| δ_l (Ergun, α = 266, β = 2.33) | 0.218 % worst, 0.176 % r.m.s. |
| δ_g | 0.416 % worst |
| δ_lg (a measurement, corrected by two densities) | 0.639 % worst, 0.260 % r.m.s. |
| **run 5 does not close** | **−5.99 %**, and it is an inconsistency between the paper's two tables |
| air viscosity recovered, linear in T | 0.0259 % worst residual over seven rows |
| that line extrapolated to Table 2's own 80 °F | 0.819 % from the printed 0.0192 cP (0.672 % once the offset measured on the liquid rows is divided out) |
| average-pressure claim, one phase, horizontal | exact — proved, measured at 1.16e-13 |
| the same with the geometric mean of the same terminals | 1.12e-05 |
| one-section recipe over a 4 ft bed, two phases | 0.332 % on the pressure drop |
| bed length at which it reaches 1 %, **root-found**, on Table 3's largest-friction-loss row | 5.6012 ft — inside the paper's own 7 ft column |
| the same root-find over the other two-phase rows | up to 466.9537 ft; 1 of the 6 inside the 7 ft column, 4 inside the 40 ft reactor |
| sectional march, observed order | 1.9999 (the midpoint rule) |
| extrapolated pymrm vs adaptive quadrature (the reference) | 7.5e-10 relative — 6.6 orders below the recipe error it measures |
| the same pair on the *other* outlet read | 1.1e-09 — so the choice of read cannot reach a reported number |
| what the outlet face transports, off the operator | (9 P_N − P_{N−1})/8, i.e. `compute_boundary_values` — **not** v·P_N |
| what the outlet reconstruction does to the last step | shortens it from h to (8/9)h, displacing P_N by h·S/9 — asserted |

## The run-5 discrepancy, stated and not repaired

Run 5 has no gas, so R_l = 100 %, so Eq. (21)'s mixture density equals the
liquid density and the manometer correction of book p. 236 cancels exactly:
Eq. (20) reduces to δ_lg = ΔP/ΔL. **Table 3 prints `2.321`; Table 4 prints
`2.469`.** The run's liquid rate is corroborated twice inside Table 4 (it
returns the printed N_Re = 2,739 and, through Eq. (17), the printed
δ_l = 2.467), and the same procedure closes the other seven rows — run 160
included, at +0.134 %.

**Run 160 is the mirror single-phase row, not the same reduction with the
phases swapped.** With R_l = 0 nothing cancels: the manometer correction still
applies in full and Eq. (21) still contributes ρ_g, and dropping the correction
alone would move that row by +22.29 %. It corroborates the *ingredients* of
run 5's reduction, not the reduction itself.

**Both readings are weighed on the paper's own scale.** Taking Table 3's cell
makes the two-phase parameter 0.9408, −5.92 % off the unity the paper says the
parameter approaches as χ → ∞ (book p. 238) — a sentence that is asymptotic and
about the Figure 7 cloud, whose scatter the same paper puts at a 13 % standard
deviation with 87 % of points inside ±20 %. On that band −5.92 % is comfortably
inside, and run 160 at χ → 0 is itself +0.90 % off unity. Against it, run 5's
2.469 sits 0.081 % from Ergun's δ_l — 11.5 times closer than run 160's 0.93 %.

**That second arm is the same statistic as the first, and it rests on a named
hypothesis.** Run 5 has no gas, so |δ_lg/δ_l − 1| *is* Table 4's own two-phase
parameter minus one — the quantity the χ → ∞ check reads, and the notebook
asserts the identity — so this is one datum read in two directions. And against
the paper's 13 % scatter *alone* a 0.081 % residual is the **more** probable of
the two numbers, not the less: it is evidence only under the hypothesis that
Table 4's cell was **back-computed** from Eq. (17) rather than measured, with
run 160 the only calibration for "too close". **The page does not choose**, and
the CSV sidecars forbid reusing that row as a measurement.

## A second printed discrepancy: Eq. (23) is printed twice

Book p. 237 prints the liquid-saturation fit as
`log_10 R_l = - 0.774 + 0.525 (log_10 chi) - 0.109 (log_10 chi)^2` and, inside
Figure 7 on the facing column, as
`log_10 R_l = -0.744 + 0.525 log_10 X - 0.109 (log_10 X)^2`. Every other
coefficient agrees, and so does the whole of Eq. (22), which Figure 7 also
prints. The difference is worth 7.15 % in R_l. Against the six Table 4 rows
inside the declared range the Figure 7 reading is **better in r.m.s.**
(18.33 % against 19.55 %) and **worse on the paper's own "within 20 %" count**
(3 rows outside against 2) — so the printed data do not settle it and the page
adopts neither. Both readings put the ceiling below a full column (0.7730
against 0.7214) and both give a negative discriminant at R_l = 1.

Figures 3, 5 and 6 carry the same curves but print only the labels
`EQUATION (22)` and `EQUATION (23)`, not the coefficients; each was cropped and
read.

## Fit versus test

Eqs. (22) and (23) were fitted to the body of data these eight rows are a
sample of, so comparing them here is a **goodness of fit, not a test**, and
that is said in the notebook, in `meta.yaml`, in `../models_entry.yaml`, in
the CSV sidecars and in the case yaml. Every agreement number is printed beside
a null that uses no correlation: Eq. (22) beats the best single constant by
4.05× in log units; Eq. (23)'s null is 117.9 % r.m.s. against its own 19.55 %.

**And the 4.05× is printed with the number that does not flatter it.** Two of
those eight rows are the χ → 0 and χ → ∞ limits, where Eq. (22) returns exactly
1 for *any* constants, so the fit is not free there while the null is charged
in full. On the six rows that can discriminate the gain is **2.48×** and the
r.m.s. **14.79 %** — the eight-row headline is 63 % larger, and 14.79 % is the
figure to set against the paper's own 13 %. Eq. (23)'s null needs no such
split.

## Scope

**Nothing is digitised.** The paper's nine figures carry the 600+ points behind
the correlations, none of them is tabulated, and the brief ranks a worked
example above a digitised figure. So the page says nothing about the scatter in
Figure 7, the foaming deviations of Figure 8, or the Lockhart–Martinelli
comparison of Figure 9 beyond the words the paper prints.

**Reference 4 — Larkins's 1959 thesis — is not on disk.** It is the support the
paper gives for the average-pressure claim the page tests, so what is shown
there, and under what assumptions, is unknown here.

**The column is a design calculation, not a reproduction.** No run in Table 3 is
a seven-foot bed; the demonstration uses run 63's own rates and temperature with
the inlet at the printed 90 lb/sq.in. gauge equipment limit, and is validated
against three independent routes rather than against a measurement.

**Run 63 is the extreme row and the page says so.** It carries the largest
friction loss of the eight, which is the condition the paper's own step 1 names
("if the bed is of great length *or the friction loss is large*"), so the
5.6012 ft is a worst case. The same root-find over every two-phase row spans
5.6012 ft to 466.9537 ft; `Cfg.demo_run` selects the row and a break row moves
it, so the evaluation point is not hard-wired.

## Files

```
index.ipynb   nine sections, executed clean (about 27 s)
build_page.py regenerates index.ipynb
meta.yaml     page metadata
data/larkins1961-table3-measured.csv       Table 3, measured
data/larkins1961-table4-calculated.csv     Table 4, the authors' arithmetic
data/larkins1961-printed-constants.csv     every constant, with its book page
```

Each CSV carries a `.meta.yaml` sidecar with a `columns:` block and a
`READ_THIS_FIRST:` note.

## Regenerate

```bash
python build_page.py
python -c "import nbformat; from nbclient import NotebookClient; \
nb=nbformat.read('index.ipynb',as_version=4); \
NotebookClient(nb,timeout=1800,kernel_name='python3', \
resources={'metadata':{'path':'.'}}).execute(); \
nbformat.write(nb,'index.ipynb')"
```
