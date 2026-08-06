# A4.5 — Fuller–Schettler–Giddings binary gas diffusivities

The most-used empirical estimate of binary gas-phase diffusivity, rebuilt from
its own printed tables — and the difference between the number it is quoted for
and the number it earns.

**Source.** Fuller, E. N., Schettler, P. D. & Giddings, J. C., *"A New Method
for Prediction of Binary Gas-Phase Diffusion Coefficients"*, **Ind. Eng. Chem.
58(5) 18–27 (1966)**,
[doi:10.1021/ie50677a007](https://doi.org/10.1021/ie50677a007). On disk as
`Fuller1966-diffusion-volumes-IEC58-18.pdf`, read in full. `pdfimages -list` gives **300 × 300 ppi
native**, so every render is `pdftoppm -r 300`; 600 dpi would interpolate.
Equations 1–5 and all 27 diffusion volumes were read off cropped page images,
and every one of the 340 data rows was verified against the physics and, where
it failed, re-read on a crop.

**The 1966 erratum, `10.1021/ie50680a601`, is NOT on disk and was NOT
consulted.** It is one page at the same publisher and it corrects exactly the
quantity this page is about. The 1966 table is used as printed; no "corrected"
table from memory has been applied.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page.
- `data/fuller-1966-table1-diffusion-volumes.csv` — Table I, the 27 fitted
  diffusion volumes. Tier 6 (fitted constants, not measurements).
- `data/fuller-1966-table3-diffusion-data.csv` — **Table III in full: 340
  measured binary diffusivities** over 90–1200 K, from 153 systems, each with
  its literature reference number and the FSG error the authors printed beside
  it. Tier 5 — measurements, but collected from the literature and every one of
  them inside the correlation's own fitting set.

## The result

**The famous 4.3 % is a goodness of fit, not a validation.** The 27 diffusion
volumes were fitted by nonlinear least squares to the same 340 measurements the
4.3 % is reported against, and the authors say so: *"Despite the large sample,
this is not a totally objective test since the parameters were obtained from
these particular data."* The page labels it that way in the notebook and in all
four metadata files.

**The paper contains no held-out data, but it prints enough to construct a
held-out test for half the method.** The atomic-increment scheme is additive, so
an organic vapour can be deleted from the fitting set entirely and still
predicted from the increments the other compounds fix. Leave-one-compound-out
over all 42 organic vapours, refitting all 27 volumes each time:

| on the 153 organic-vapour points | mean \|error\| |
|---|---|
| published volumes, in sample (Table II prints 5.0) | 5.03 % |
| own full refit, in sample | 5.04 % |
| **leave-one-compound-out, genuinely held out** | **6.61 %** |

So for an organic vapour that was *not* in the 1966 set, the honest expectation
is nearer **6.6 %** than 4.3 % — and much worse where a compound is the only
carrier of a structural feature (thiophene 3.4 → 13.7 %, methane 5.0 → 11.8 %,
benzene 5.8 → 11.4 %).

**The other half of the method admits no held-out test at all.** Each of the
twenty simple molecules carries its own fitted Σv, so deleting one deletes the
parameter that predicts it. Those 187 points cannot be cross-validated by any
means, and the 6.6 % must not be transferred to them. The page says this rather
than quietly averaging over it.

## Validation (tier 5 data; the correlation comparison is a fit, not a test)

**Reproduction of the printed tables** — this checks the transcription and the
implementation of eq. (4), not the physics, and is described that way. 338 of
340 printed FSG errors to within **0.043 percentage points** (median 0.0038),
inside the table's own rounding. Twelve printed aggregates fall out of the same
transcription without being fitted to: all eight of Table II's FSG category
means, the `AVE. ERROR` 4.3211 against a printed 4.32, the >10 % tail 7.35
against 7.4, the system count 153 against 153, and the `STANDARD DEVIATION`
6.7101 against 6.71 — the last reproduced exactly as √(Σe²/(N−p)) with p = 28,
a definition the paper never states but whose parameter count is printed
(27 rows of Table I plus b).

**A second route sharing no code.** Within one system from one reference, eq. (4)
forces D to scale exactly as Tᵇ — an identity involving only the printed
temperature, observed D and printed error, with no volumes, no molar masses and
no evaluation of eq. (4). Median deviation **3.1e-5** over 142 ratios in 51
series, worst 1.0e-4. It also recovers **b = 1.7492** from the printed table
alone.

**Two printed rows cannot be reproduced** and are reported rather than dropped:
He–H₂O at 352.5 K (printed −5.33, recomputed +0.21) and CO₂–ethylene oxide at
298.0 K (printed 15.09, recomputed 14.27). Both re-read at 3× magnification and
confirmed as printed; the first is flagged independently by the Tᵇ identity.

**The objective function does not match the published table.** φ′ from Table
III's own error column is 1.4048, from eq. (4) with the printed volumes 1.3996,
and a free refit reaches 1.3701 — all below the paper's stated constrained
minimum of **1.4762**. Since eq. (4) is a special case of eq. (2), a minimum of
eq. (2) cannot lie above a value reachable inside eq. (4). Reported, not
resolved.

**The paper's own claim about temperature exponents tested against its own
data**: 20 of the 26 multi-temperature series fall in the stated 1.6–1.8 band,
median 1.704 against the fitted 1.749. But the per-system exponents span
1.58–1.95 — **14× wider than the ±0.013 confidence interval on b**, which is a
statement about a global fit and not about any one gas pair.

**The transport solve** is validated against an adaptive quadrature of 1/Γ that
never forms a grid and never calls pymrm: 3.6e-7 relative at n = 400, grid order
**2.002**, with the quadrature's own error 1.1e-14.

**Break table: 14 injected defects.** C 16.5 → 16.0 moves the reproduction from
0.043 to 1.65 pp; a widely-circulated alternative increment set (C 15.9, H 2.31, O 6.11, N 4.54 — not in any document on disk) to 1.99; swapping Σv(N₂)/Σv(O₂)
to 3.48; b → 1.749 to 0.82; M(D₂) → 4.028 moves the fourteen D₂ rows 0.005 → 0.335;
and the raw text-layer parse before the image corrections gives **1520 pp**. On
the transport side n = 5 moves the flux deviation by 6.7e3, `nu = 1` by 2.8e6, a
wrong outlet BC by 2.8e6. **Two rows are labelled structural and do not move**:
a wrong exponent applied to *both* routes (the check reads the same Γ twice, so
it is blind to a wrong b, a wrong D and a wrong molar mass), and the
cell-to-cell flux spread, which is conservative by construction at ~1e-20 and
also sits below `check_agreement.py`'s `ABS_FLOOR = 1e-12`.

## What pymrm adds

Not much to the correlation — eq. (4) is four arithmetic operations and the
reproduction, the cross-validation and the exponent measurement are all
numpy/scipy. What it adds is the step a correlation cannot take: turning an
error in D into an error in an answer. For benzene through air across 300–600 K,
moving b across the range measured from the paper's own data costs **−5.9 % to
+7.3 %** in the flux if D is anchored at 300 K and **−0.2 % to +0.1 %** if it is
anchored mid-range — the exponent is nearly free or expensive depending entirely
on where you pin D, which the paper does not say. A *scale* error in D is
structural: the problem is linear in Γ, so 4.3 % in D is exactly 4.3 % in the
flux and nothing is learned by computing it.

## Scope

Two printed defects are recorded: the two irreproducible FSG entries, and the
fact that the `AVE. ERROR`/`STANDARD DEVIATION` rows at the foot of Table III
are in **Table II's** method order rather than Table III's own column order.
Only the FSG column is transcribed; anyone using the other eight should settle
that ordering first.

`A4.6` (Chapman–Enskog, from Chapman & Cowling) is built separately and prints
~45 measured D₁₂ pairs in its Table 22. The first-principles result set against
this empirical one on a single measured-D₁₂ axis is the natural comparison and
**`A4.6` should own it**. It is not attempted here and no number on this page
comes from Chapman & Cowling.

Runtime ~4 s.
