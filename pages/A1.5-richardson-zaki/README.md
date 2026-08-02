# A1.5 — Richardson–Zaki

`u/u_t = ε^n` is one of the most-quoted expressions in particle technology, and it
is not what the paper writes. The paper writes `log V_c = n log ε + log V_i` and
then spends four pages establishing that **`V_i` is not `V_0`** in fluidisation.
The famous form is the sedimentation special case with the wall term dropped.

Richardson and Zaki (1954) is a **data paper**: eight printed tables, 66 measured
exponents, two wall laws, and only then the correlation, as *five* expressions
covering four Reynolds windows. This page transcribes all eight tables off 600 dpi
renders, certifies the transcription against four identities the tables carry,
costs the correlation against a null model window by window, decides between the
two wall laws on the authors' own numbers, and puts the flux function into a pymrm
conservation law.

- **Structures:** `S1` (the pymrm solve; the correlation itself has no PDE in it)
- **Runtime:** ~15 s
- **Data tier:** 2 — printed tables of the authors' own laboratory measurements
- **Nothing here is digitised.** No figure point was read.

## Source, and which pagination is cited

**Richardson, J. F. and Zaki, W. N.**, *Sedimentation and fluidisation: Part I*,
Trans. Instn Chem. Engrs **32**, 35–53 (1954).

**The document on disk is the verbatim reprint** in *Trans IChemE* Vol. 75,
December 1997 Jubilee Supplement,
[doi:10.1016/S0263-8762(97)80006-8](https://doi.org/10.1016/S0263-8762(97)80006-8),
where the paper occupies **reprint pages S82–S100**. Every page carries the
original running footer `TRANS. INSTN CHEM. ENGRS, Vol. 32, 1954`, but **the
original page numbers 35–53 are nowhere printed on the reprint.** Nineteen
consecutive pages is consistent with 35–53; the page does not assert that range
as read. **Every page reference on the page and in the sidecars is a reprint
page, marked as such.**

**Origins cited but not consulted.** Figure 18, from which eqs. (35)–(39) are
fitted, carries Steinour's sedimentation slopes (ref. 4) and Lewis, Gilliland and
Bauer's fluidisation slopes (ref. 9) as well as the authors' own. Neither paper is
on disk. Nothing on this page separates their contribution from the authors' own,
because the paper does not.

## The correlation, as printed

```
log V_c = n log ε + log V_i                              (28)   the actual form

n = 4.65 + 19.5 d/D                        Re < 0.2      (33)
n = (4.35 + 17.5 d/D) Re^-0.03      0.2 < Re < 1         (37)
n = (4.45 + 18   d/D) Re^-0.1         1 < Re < 200       (38)
n = 4.45 Re^-0.1                    200 < Re < 500       (39)
n = 2.39                                  Re > 500       (34)

V_i = V_0                              sedimentation     (40)
log V_0 = log V_i + d/D                fluidisation      (41)
n = 2.7 K^0.16,  K = (π/6) d_s³/d_p³   non-spherical     (42)
```

Validity, stated by the authors: `Re` from 2e-4 to 7e3, `d/D` from 0 to 4e-2.
`Re` is built on the **terminal** velocity, so `n` does not move as the bed expands.

## What the page checks, in the brief's order

| route | check | result |
|---|---|---|
| 1, measurement | eqs. (33)–(39) vs. the 58 non-oil measured slopes | mean **1.20 %**, worst 6.75 % — **but this is a FIT RESIDUAL**, see below |
| 1, measurement | the **wall coefficient α**, fitted separately on each half of the data | **sedimentation (Table I, the rows that can fail) −0.235 [−1.78, +0.17]; fluidisation +1.095 [+0.83, +1.51]** — each interval excludes the other law |
| 1, measurement | the paired wall test, which removes `V_0` | slope **0.79 [0.47, 1.34]** on the **6 of 10** pairs that can resolve it |
| 1, measurement | eq. (42) vs. Table VIII, with a null baseline | rms 1.59 % vs. a constant's 7.94 % — **but 1.60 % vs. 2.66 % once the plates are dropped** |
| 2, identity | eq. (39) at Re = 500 vs. the separately-stated eq. (34) | **−0.015 %** — but the two share Table VI's last row (`n0 = 2.39` at `Re = 489`), so this is a **transcription check**, not two independent fits meeting |
| 2, identity | Table VI × Table VII, `n0 = slope × \|(d/D)_{n→0}\|` | worst 0.635 %, mean 0.196 %, 11 rows |
| 2, identity | branch continuity | **−13.9 % jump at Re = 200 and d/D = 0.04**; −1.82 % at Re = 0.2; +2.30 % at Re = 1 |
| 2, identity | the 17.5 and 18 against the division displayed above them | 17.400 and 17.800, **both rounded up**, +0.57 % and +1.12 % |
| 2, identity | `K` recomputed from the printed `d_s`, `d_p` and definition | **fails on 4 of 7 rows**, worst −45.6 % |
| 2, identity | four transcription identities | locate **five printed errors**; Stokes reproduces Table I to 0.67 % and correctly **fails** on Table II |
| 4, digitised | — | **not used. Nothing on this page is digitised.** |

## The three things a reader should take away

**The 1.20 % is a fit residual and the page refuses to call it anything else.**
Eqs. (33)–(39) were fitted to these very slopes. The informative number is the
**null baseline** beside it, per Reynolds window:

| branch | rows | correlation rms | constant-n rms |
|---|---|---|---|
| eq. (33) | 20 | 0.95 % | 3.07 % |
| eq. (37) | 5 | 1.559 % | 1.538 % — **decides nothing**; the gap reverses under all three sensitivities the page prints |
| eq. (38) | 23 | **2.11 %** | **15.80 %** — this is where it earns its keep |
| eq. (39) | 4 | 1.87 % | 3.42 % |
| eq. (34) | 6 | — | **it *is* the constant**; the mean of the six slopes is 2.3900 |

**The wall term is the one genuinely two-sided test on the page — once the rows
that cannot fail are taken out of it.** Eqs. (40) and (41) are the same statement
with α = 0 and α = 1, and both sides are printed columns. But **Table II's
intercept column *is* its `log10 V_0` column**: on all 15 rows, across two and
three tube diameters per group, the residual the fit sees is 0.0000 or 0.0004
decades, against Table I's −0.030 to +0.010. Those 14 rows cannot show a wall
effect of any size while still narrowing the interval — fitted alone they give
**α = +0.009 [+0.000, +0.012]**. So the headline sedimentation fit is **Table I
alone, α = −0.235 [−1.78, +0.17]**, which still excludes +1; the pooled fit
(**−0.064 [−0.41, +0.08]**, 4.0× narrower) is printed beside it and the exclusion
is explained. The verdict survives the correction; the precision does not.

The limits are stated: the fluidisation α is a refit of the data Fig. 21 was
drawn from, so recovering 1 is nearly guaranteed; the **sedimentation** α is the
informative half, because eq. (40) was justified on a different figure on which a
wall term would not have shown as a slope, and a fluidisation-sized wall effect
would have appeared here. It does not.

**And the coefficient 1 does not survive its own range.** Restricted to the `d/D`
range the two data sets share, the fluidisation α is **+2.267 [+1.20, +3.47]** —
an interval that *excludes* 1 — against sedimentation's **−0.101 [−1.22, +0.21]**
on the Table I rows there; and α falls from **+2.162** below `d/D = 0.05` to
**+1.033** above it, which a single linear term in `d/D` cannot do. The two-law
verdict (a wall term in fluidisation, none in sedimentation) stands; the *linear
form* of eq. (41) does not hold across the range it was fitted on. Run 12, the
sign-corrupted cell, is a Table II row and so is not in the headline fit at all:
on the pooled fit it reads −0.436 as printed and **−0.067 with the minus sign
restored**, against −0.064 excluded — a printing error, not a fragility.

**The solids balance in the pymrm section cannot fail, and the break table says
so.** It returns exactly zero for a wrong `n`, a wrong `V_0`, an unstable CFL and
a 16-fold-too-coarse grid; the only defect it catches is `nu = 1`. The check that
*can* fail is the L1 error against the exact shock, which converges at order
1.000 and moves ×111 for `nu = 1` and ×187 for a reversed upwind direction — but
only ×1.5 when `n` is changed and judged against its own shock speed. **The solve
tests the integration, not the physics.**

## What is not here

- **No experimental validation.** The data are the authors' measurements, but the
  correlation was fitted to them and there is no hold-out set. The only rows the
  authors exclude — the two oils at high `d/D` — are excluded for a reason that is
  **perfectly confounded** with the correlation's own stated `d/D` limit: inside
  the eq. (38) window every row above `d/D = 4e-2` is an oil row and there are no
  non-oil rows there at all.
- **No claim about ε outside roughly 0.42 to 0.96.** The extreme printed abscissa
  ticks read on any figure are `log ε = −0.38` and `−0.02` (Fig. 12) — printed
  axis labels, not digitised coordinates.
- **No result resting on Table V's `d_s` column**, which is the one column of the
  paper that is illegible at 600 dpi in this scan.
- **No mixed-size or mixed-shape result.** The paper's own is a single figure with
  no tabulated numbers.
- **No claim that pymrm improves the correlation.** It does not. The solve
  integrates the flux the correlation becomes downstream, and the break table
  measures exactly how little that shows.
- **No verdict on eq. (37).** Its five rows cannot separate the correlation from
  a constant in either direction; the page prints the three sensitivities that
  show why, and does *not* say eq. (37) is worse than a constant.
- **No claim that `log V_0 − log V_i` equals `d/D`.** The page claims a wall term
  in fluidisation and none in sedimentation. On the shared `d/D` range the
  fluidisation interval excludes eq. (41)'s coefficient of 1.
- **No page range asserted from the document.** The reprint numbers the paper
  S82–S100 and `35–53` appears nowhere on it; the two `cite_data` lines in the
  notebook print the catalogue's range with a flag under each saying it was
  inherited, not read.

## Related

`A1.1` (Ergun — the fixed-bed end of the same ε(u) curve), `A1.6` (Wen and Yu —
the minimum-fluidisation matching condition; **not built here**), `A1.7` (Geldart
— decides which fluidised-bed model applies at all), `A1.8` (fluidisation regime
map), `E1.1` (the two-phase theory, the gas-solid analogue; parked).
