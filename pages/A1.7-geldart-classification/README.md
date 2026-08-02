# A1.7 — Geldart's powder groups

The four letters everyone quotes come from two inequalities. This page recomputes
both from the expressions Geldart derives them from, tests the A/B one against
the 22 size fractions he measured himself, measures what the minimum-fluidization
correlation underneath it is actually worth against those same measurements, and
shows that the B/D one is the Davidson cloud-existence condition in other symbols.

- **Structures:** `S3` (the single pymrm solve; the classification itself has no PDE in it)
- **Runtime:** ~4 s
- **Data tier:** 2 — a table printed in the paper, and the values in it are the
  author's own laboratory measurements.

## Source

**Geldart, D.**, *Types of gas fluidization*, Powder Technology **7**(5) 285–292
(1973), [doi:10.1016/0032-5910(73)80037-3](https://doi.org/10.1016/0032-5910(73)80037-3)
— the paper itself was read, not a reprint of it, on 600 dpi renders of journal
pages 285–292. The PDF text layer of this Elsevier scan renders the paper's own
title as "Types of Gas Fhidization" and was not used for anything.

**Origins cited but not consulted.** Equation (3), the minimum fluidization
velocity, carries Geldart's reference 11: *L. Davies and J. F. Richardson, Trans.
Inst. Chem. Engrs. 44 (1966) T293*. Not on disk; used exactly as Geldart prints
it. Same for equation (9) (Verloop and Heertjes, ref. 10) and equation (10)
(Oltrogge, ref. 36). The Davies–Taylor coefficient 0.711 used in section 5 comes
from `E1.2`, not from this paper.

## The two inequalities

```
group A :  (rho_s - rho_f) d'   <= 225      eq. (6),  d' in um, rho in g/cm3
group D :  (rho_s - rho_f) d'^2 >= 1e6      eq. (8)
```

Recomputed from equations (5) and (7) with Geldart's own `g = 981`,
`mu = 1.8e-4`, `K_MB = 100`, `eps_0 = 0.4` and `d_B = 25 cm`, they are **229.36**
and **1.0159e6** — his printed values are 1.94 % and 1.59 % low. Neither rounding
changes any verdict in his own table: no row falls in the 225–229.36 gap, and
none in the 1e6–1.0159e6 gap either, the largest `(rho_s-rho_f)d'^2` in Table 1
being 1.192e5.

## What the page checks, in the brief's order

| route | check | result |
|---|---|---|
| 3, measurement | A/B criterion vs. the measured `U_MB/U_0` of Table 1 | **2 of the 3 rows that can discriminate** (whole table 20/21, against a null baseline of 19/21) |
| 3, measurement | eq. (3) vs. the measured `U_0` column | implied constant spans 0.46–2.94 × the printed 8e-4, median 0.981 |
| 3, text | the worked example on journal page 289 | worst 1.94 % over four velocities — a coincidence of digits with the row below, and an unrelated quantity |
| 2, identity | eqs. (6) and (8) recomputed from (5) and (7) | +1.94 %, +1.59 % — and no Table 1 row in either gap |
| 2, identity | solids balance linking `eps_MB` and `H_MB/H_0` | all 22 rows give an admissible `eps_0`, 0.441–0.570 |
| 2, identity | eqs. (9) and (10) placed relative to (6) | eq. (10) crosses XY inside the plot at 104 µm; eq. (9) lies 2.4–8.0× to the right of XY across it |
| 2, cross-page | eq. (7) is E1.2's no-cloud condition | 0.711 vs. 1/√2 moves the boundary by 0.55 % |
| 4, digitised | — | **not used. Nothing on this page is digitised.** |

**The A/B score is reported the honest way round.** 19 of the 21 comparable rows
are measured group A, so a null model that predicts group A for every powder
already scores 19/21 (Diakon 6/8) against the criterion's 20/21 (7/8) — one row
of information. The only rows on which the two can differ are the three the
criterion calls group B, and there it is right 2 of 3. The margin also rests on a
**tie-break**: two Diakon rows print `U_MB = U_0` exactly, eq. (4) is *displayed*
as `>= 1` but *applied in section 4.3* as `> 1`, and the strict form is used
because it is the only one that reproduces Geldart's own verdict on those rows.
Break table 2 costs the alternative at 18/21 (Diakon 5/8, discriminating rows
0/3, null baseline 21/21).

On the single row where eq. (6) and the strict ratio part company — the 220 µm
Diakon fraction — **Geldart's own prose sides with eq. (6)**: section 4.3 puts
that fraction in group B, as eq. (6) does, against the mechanical
`U_MB/U_0 = 1.029`.

Three deliberate-break tables measure what each check can and cannot catch, and
four blind spots are named rather than asserted away: `rho_f` moves nothing
anywhere, the eq. (3) *spread* is nearly blind to a wrong exponent while its
median is not, a 0.10 error in one `eps_MB` is invisible to the voidage check,
and the whole-table A/B score is nearly powerless because the classes are
lopsided.

`eps_MB` and `H_MB/H_0` are described as **reported**, not measured: section 4.2
states the method for `U_0` and `U_MB` only and says nothing about those two.

## What is not here

- **No C/A boundary.** Geldart prints no expression for it — on his Figure 3 it is
  a hand-drawn shaded band.
- **No test of the B/D boundary against data.** Table 1 stops at 318 µm and
  1.5 g/cm³; at that density equation (8) does not begin until 817 µm.
- **No claim that pymrm improves the classification.** It does not. The one solve
  on the page answers a question the classification raises and does not answer —
  how big the cloud is on either side of the B/D line — and reuses `E1.2`'s
  operator to do it.

## Related

`E1.2` (the Davidson bubble — owns the percolation solve and the cloud threshold),
`E2.1` (the bubbling-bed model, which presumes a cloud), `A1.6` (Wen and Yu —
the obvious replacement for equation (3); **not built here**), `A1.1` (Ergun),
`A1.8` (fluidisation regime map).
