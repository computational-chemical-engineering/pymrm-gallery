# A4.1 — Fick's law with Wilke's mixture rule for the effective diffusivity

The weighted-harmonic mixture rule

```
D'_A = (1 - y_A) / ( y_B/D_AB + y_C/D_AC + y_D/D_AD + ... )
```

printed in the abstract of Wilke (1950), reprinted forty-three years later as
eq. (6.1.14) of Taylor and Krishna with `(Wilke, 1950)` inline. Every reactor
model that carries one diffusion coefficient per species in a multicomponent gas
is using some version of it.

**The rule is exact when the other species are genuinely stagnant.** That is not
an approximation but an algebraic identity, so a page that reports agreement
there reports arithmetic. This page therefore goes where it is not exact — and
the paper itself supplies the place: eight ternary test examples, each solved
three ways, with the exact Stefan–Maxwell solution printed beside the rule as
Table 1.

## Headline

Against the exact solution of Wilke's own eight problems, computed here:

| | max error | median error |
|---|---|---|
| **Method 1** — the rule at the arithmetic-mean composition | **86.6 %** | **21.6 %** |
| **Method 2** — the *same* rule at a flux-weighted film composition | 12.1 % | 0.2 % |
| **null** — weighted *arithmetic* mean instead of harmonic | 212.7 % | 22.2 % |

Two readings of that table matter more than the numbers themselves.

- **The split is by flux share.** On the species carrying the majority of the
  net flux the rule is 7.7 % off at the median and never worse than 2.4 % in
  Examples IV, V, VI, VII and VIII — every example whose dominant species has
  φ ≈ 1. On the minority species it is 24–87 % off, and
  the sign is not even consistent (+87 % on Example IV, −41 % on Example VII).
- **The harmonic mean is not what limits it.** Swapping it for the arithmetic
  mean that J. H. Arnold attributes to Hougen and Watson *in the paper's own
  printed Discussion* moves the median by 0.6 points while making the worst case
  two and a half times worse (86.6 % → 212.7 %). Keeping the harmonic mean and
  moving the **composition** takes the worst case from 86.6 % to 12.1 %.
- **Wilke's own exact column is not uniformly good either.** Reproducing his
  Method-3 column against the page's exact solution, fifteen of the sixteen
  cells are within 2.31 % — and Example VI's `N_A` is 10.42 % out. That cell is
  the minority-flux species again, and it is his hand solution of the
  transcendental pair, not a reading error. It is reported at full weight; the
  page excludes nothing from the Method-2 or Method-3 reproduction metrics.

## Source, and which Wilke it is

C. R. Wilke, *"Diffusional Properties of Multicomponent Gases"*, **Chemical
Engineering Progress 46**(2) 95–104, February 1950. No DOI; pre-DOI journal. The
file on disk is a microfilm scan whose text layer is among the two worst in this
library — 182 characters extract from journal page 97 and the by-line comes back
as `C.R. WI E` — so **no character on this page comes from the text layer**.
`pdfimages -list` reports JPEG grayscale at 400 × 400 ppi native; everything was
rendered at that resolution, cropped and re-read at digit scale. Table 1 is
printed rotated 90° and was rotated losslessly before reading.

Identity confirmed from the document itself (title page, by-line, affiliation,
running feet `Vol. 46, No. 2` / `CHEMICAL ENGINEERING PROGRESS` /
`February, 1950`, pages 95–104, the 30-item Literature Cited and the signed
Discussion) and independently by Taylor and Krishna's reference list.

**This is not the other Wilkes.** It is not the Buddenberg–Wilke *viscosity*
mixture rule (*Ind. Eng. Chem.* **41**, 1345, 1949; *J. Chem. Phys.* **18**, 517,
1950) that Bird/Stewart/Lightfoot and Chapman & Cowling carry when they cite
"Wilke" — that rule appears in *this* paper too, as eqs. (39)–(43) and Table 2,
and is scoped out as a different catalogue case. It is not Wilke & Lee's
binary-diffusivity **estimation** correlation (*Ind. Eng. Chem.* **47**, 1253,
1955), a different paper, journal, year and equation. It is not Wilke–Chang.

**Second witness.** Taylor & Krishna (1993), *Multicomponent Mass Transfer*,
Wiley, read at its native 600 ppi: section 6.1.3 (printed p. 126) prints the
rule as eq. (6.1.14) with the attribution inline and states its condition;
section 8.6 (p. 204) records the pedigree; Example 8.5.1 (pp. 201–203),
"based on experiments conducted by Fairbanks and Wilke (1950) with a view to
assessing the validity of Wilke's effective diffusivity formula", is reproduced
in full. Note for `docs/papers-inventory.yaml`: the pedigree sentence is in
section **8.6**, not 8.7, and the exact/linearised comparison is sub-section
8.6.1 (p. 208) together with 8.5.5 (p. 203).

## Files

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page.
- `data/wilke-1950-table1.csv` — every cell of Table 1 (tier 6).
- `data/wilke-1950-printed-intermediates.csv` — the sample calculation, the
  Example-VI film-pressure-factor paragraph, the Example-IV diffusivity
  averages, and the paper's own accuracy claim (tier 6).
- `data/taylor-krishna-1993-example-8-5-1.csv` — the second witness's worked
  example (tier 6).

**No other page's dataset is loaded**, so the cross-page reconciliation rule does
not apply. No figure is digitised and no page image exists anywhere in this
directory.

## Three printed defects, proved from the paper's own numbers

1. **Two missing minus signs.** Table 1's Method-1 `N_B` for Examples I and II
   is printed without a minus, while the structurally identical Example III
   prints one. Recomputing Method 1 gives −5.80 and −0.141: the magnitudes agree
   to 0.06 %, so the signs are lost, not the numbers.
2. **Example VI, Method 1, `N_A`.** The printed −0.104 implies a film pressure
   factor of 380 mm for A. The paper's own text on journal p. 99 states that
   quantity: *"In Example VI the film pressure factor for the diffusion of A as
   calculated under Method 1 is 288 mm."* Solving Method 1 gives −0.193 and
   **287.4 mm** — 0.21 % from the printed 288 — and the same paragraph's other
   three numbers (75, 72.5, 147.5 mm) reproduce exactly. The paragraph is
   internally consistent; the table cell is the outlier, and it happens to sit
   within 1 % of the Method-2 entry in the next column.
3. **Eq. (22)'s first bracket.** Printed as `[1/D_AB − 1/D_AC]`; every other
   subscript in (22) is (21)'s under A↔B with `D_AC → D_BC`. Shown symbolically:
   (21) reproduces the Maxwell–Stefan balance identically, (22) does so only
   with `D_BC`, and — decisively — (21)+(22) integrates to the paper's own
   eq. (23) only with `D_BC`. Eq. (23) is not in doubt: eight of the ten printed
   Method-3 pairs satisfy it inside 1 %.

All three are **reported, never repaired**: the CSV stores every cell as
printed. So are two printed *typos*, quoted with `[sic]` rather than silently
corrected: journal p. 100 prints *"a maximum **derivation** of less than 10%"*
(crisp at 400 ppi, not a degraded "deviation"), and Taylor & Krishna's printed
p. 126 reads *"When species i diffuses through **a** n − 1 stagnant gases"*.

One reference number in the paper is also wrong and is flagged rather than
silently mapped: Wilke cites "(25)" for both *"Sherwood (25) presents
differential equations"* and *"Gilliland (25) has integrated (21) and (22)"*,
but his ref. 25 is *Sherwood, T. K., Ibid., p. 11 (1937)* — the same book as
ref. 24 — while Gilliland & Sherwood (1934) is his ref. **12**.

Two further cells are left unrepaired because the constraint does not decide
them. Example I's third Method-3 set, `(4.14, −3.93)`, misses eq. (23) by
−13.2 % against a 0.73 % rounding bound; holding `N_A` the equation wants
`N_B = −3.723` and holding `N_B` it wants `N_A = 4.302`, and nothing chooses
between them. Example III's pair misses by +9.4 % against a 0.38 % bound and no
single-digit repair reconciles it at all.

## Two ambiguous glyphs, settled by arithmetic

Table 1's header multiplier and the exponent in the p. 101 line
`N_A = −0.105 × 10^?` are both **illegible at 400 ppi** — no amount of zoom
resolves either, and the pixel shape is not used. Both are fixed to `1e5` by
arithmetic: Method 1 for Example I recomputes to 5.799e-5 against a printed
5.80, and the p. 101 sample calculation's own nine-step chain closes on
−1.049e-6 = −0.105 × 1e-5.

## What the page cannot conclude

**Nothing measured.** Every dataset is tier 6 — the authors' own computed
numbers. The paper's only experimental content is Figure 1 (from Fairbanks's
1948 M.S. thesis), which is a figure, is scoped out, and is not digitised. The
paper's own Comment section says validity of the Stefan–Maxwell equations *"has
not been established experimentally in multicomponent systems for other than the
stagnant gas case cited here"* — so the exact solution is used here as the
**reference**, and if those equations are wrong for these mixtures then so is
every error bar on the page. This is a model-versus-model reproduction, not a
validation.

**Nothing about these gases.** Table 1's own footnote: the diffusion
coefficients "are estimated and do not correspond to true values".

**Nothing beyond three components.** Every example is ternary with exactly one
stagnant species. Wilke's own extrapolation of Method 2 to more components is
quoted and is not tested.

## Related pages

`A4.2` reaches a related conclusion about scalar effective diffusivities from
Krishna & Wesselingh (1997) on the Duncan–Toor two-bulb measurements. It shares
no dataset, no solver and no number with this page — deliberately, so the two
are independent evidence rather than one argument told twice. `A4.3`, `A4.4`,
`A4.9` and `H1.9` are the rest of the multicomponent-diffusion cluster.
