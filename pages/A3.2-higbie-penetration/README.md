# A3.2 - Higbie penetration theory

A liquid surface is exposed to a gas for a definite length of time, absorbs by
unsteady diffusion the whole while, and is then swept away. Nothing reaches a
steady state and no film thickness appears anywhere, so the coefficient is
`k_L = 2 sqrt(D / pi t_exp)` and `k_L ~ D^(1/2)` against film theory's `k_L ~ D`.

## Where the model comes from

R. Higbie, *Trans. AIChE* **31**, 365-389 (1935) is the origin. It is pre-DOI,
it is not on disk, and **it was not consulted**. Everything here is read from
**Bird, Stewart & Lightfoot, *Transport Phenomena*, 2nd edn (Wiley, 2002),
section 18.5** (book pp. 558-561), which derives the result, names the
"penetration model", attributes it to Higbie in the same paragraph and prints the
full citation in a footnote that also credits him with having "provided the basis
for the 'penetration model' of mass transfer". Every equation number on the
page is BSL's, and every digit was read on 300 dpi renders at the scan's native
resolution - its ABBYY text layer destroys equations and was used only to find
pages.

That one sentence and that footnote are **all** BSL prints about the 1935 paper,
so nothing on the page is attributed to it beyond what BSL prints and - where the
Froment cross-check is being described - what Froment prints about it (his
"Higbie's uniform age"). The two geometries the contact time is read off - the
falling film (`t_exp = L/v_max`) and the rising bubble (`t_exp = D/v_t`) - are
BSL's section 18.5 and BSL's Example 18.5-1; the page does not claim they are
Higbie's.

**Froment, De Wilde & Bischoff, 3rd edn (Wiley, 2011), section 6.4** is a
cross-check only: an independent statement of the same constant, the closed-form
enhancement factor for Higbie's uniform surface age, and Table 6.4.2.1.
Levenspiel names and cites Higbie but neither derives nor tests him, and is not
used.

## What the page establishes

- **Problem 18A.4 reproduced two ways, each scoped to what it tests.**
  Eq. 18.5-18 gives 0.27351 g-mol/hr against the printed 0.273 (+0.188 %); a
  pymrm march of Eq. 18.5-11, read as the outlet mixing cup and touching no error
  function, gives the same to 4.6e-4 %. The two routes are *not* independent end
  to end: the wetted width, the 3/2, the solubility conversion and the units
  enter both and cancel out of the ratio, so the 4.6e-4 % tests the dimensionless
  `Phi(Lambda)` against `sqrt(4 Lambda / pi)` and nothing else. The **+0.188 %
  against the printed 0.273 is what tests the prefactor.** Break rows show the
  invariance directly. The +0.188 % is a hair *outside* three-figure rounding -
  0.273513 rounds to 0.274 - and the page says so: two inputs 18A.4 does not
  print are each worth about that much (rho = 0.999 gives 0.273239, M(Cl2) = 71.0
  gives 0.273147, both rounding to 0.273). The gap is not theirs alone, though:
  only
  1.28e-5 g-mol/hr (0.0047 %) lies outside the window, the printed inputs are
  three-figure too - half a unit in the last place of D is +/-0.199 % in W_A by
  itself, the whole +/-0.183 % window - and truncation rather than rounding would
  put 0.273513 inside. Neither unprinted value is adopted.
- **A printed defect in the pinned source.** Problem 18A.7(a) says "Use
  Eq. 18.5-20" and prints Eq. 18.5-19's answer. The two differ by exactly
  sqrt(3). Reported, not repaired.
- **The one empirical test on disk.** For the same bubble the penetration
  prediction carries no adjustable constant, `t_exp = D/v_t`, and comes out
  **12.0 % below** Hammerton & Garner's measured `k_c`. Creeping flow would be
  49.2 % below.
- **The three-way comparison, taken.** Froment's Table 6.4.2.1 is the only place
  in either book with film, surface renewal and penetration side by side. Its
  surface-renewal column reproduces everywhere; its film column is high by one
  unit in the last printed place at two of six; **its penetration column
  disagrees with the formula printed at its own head at all six**, and two cells
  are impossible on the printed digits alone (0.94 where an enhancement factor
  cannot be below 1, and 10.39 above the same row's surface-renewal 10.05).
- **The penetration approximation, decomposed.** BSL justifies Eq. 18.5-11 with
  two arguments as though they were one. Solving Eq. 18.5-7 - which BSL sets up
  and then declines to solve - shows they are worth a factor 4810 apart: at
  Lambda = 0.1 the velocity profile costs 1.92 % and the finite film 4.0e-6. The
  1 % threshold is `Lambda* = 0.0564`, i.e. penetration depth below 0.475 delta.

## What it does not do

It does **not** test the sqrt(D) exponent. One measured coefficient at one
diffusivity cannot resolve an exponent, and the page carries the search behind
that statement rather than asserting it. The exponent question stands exactly
where [`A3.3`](../A3.3-danckwerts-surface-renewal/) left it: film theory and
surface renewal each fit a single datum exactly, one free constant each, zero
residual each. What is new is that on this datum penetration theory has no free
constant at all - so it is the only one of the three that can be wrong, and it
is, by 12 %.

The measurement is quoted from BSL, not consulted; BSL prints no error bar for
it and none is invented. The film thickness of Problem 18A.4's column is not
printed and is not inferred.

## Files

```
index.ipynb    the page
build_page.py  regenerates index.ipynb
meta.yaml      page metadata
data/          three CSVs of printed numbers, each with a provenance sidecar
agreement.json 43 metrics, written by gallery_utils.report_agreement
```

Runtime about 40 s. Nothing is stochastic; two executions give identical
content and an identical `agreement.json`.
