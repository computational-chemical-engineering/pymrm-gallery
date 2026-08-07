# J1.3 — BET 1938: the theory that yields the surface, against the theory that has to be given it

Almost everybody meets the BET equation as a line drawn through five points on a
$p/v(p_0-p)$ plot. That line is a **two-parameter fit over a fifth of an
isotherm's range**, and it is supposed to be good. The theory Brunauer, Emmett
and Teller were arguing against could fit the same isotherms, and had been doing
so for nine years.

So this page is organised the way the 1938 paper is. Section I of that paper is a
**quantitative refutation of a rival** — the polarization theory of de Boer,
Zwicker and Bradley — carried out before a single line of BET algebra appears.
That refutation is the page's spine, because it is the only part of the argument
where two theories are made to disagree about something measurable.

**Source.** Brunauer, S., Emmett, P. H. and Teller, E. (1938). *Adsorption of
gases in multimolecular layers*. J. Am. Chem. Soc. **60**(2), 309–319,
[doi:10.1021/ja01269a023](https://doi.org/10.1021/ja01269a023). Received
November 19, 1937. **It is the only document read.** Identity confirmed from its
own title page on a native-resolution render: the running head "Feb., 1938 …
309", the contribution line "[Contribution from the Bureau of Chemistry and Soils
and George Washington University]", the title and the by-line "By Stephen
Brunauer, P. H. Emmett and Edward Teller".

`pdfimages -list` reports every page as **CCITT-G4 bilevel at 300 ppi native**,
so pages were rendered at 300 ppi — a larger render would only interpolate — and
**every numeric was cropped out and re-read at digit scale**. The PDF's text
layer was not used for any digit. That discipline changed one reading: at page
scale the Fig. 1 legend appears to name catalyst "358" where Table I has 958; at
crop scale it reads 958 and the two agree.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page. Runtime ~6 s.
- `data/bet-1938-table1-nitrogen-90K.csv` — Table I, twelve nitrogen isotherms.
- `data/bet-1938-table2-e1el-gases.csv` — Table II, four gases on two adsorbents.
- `data/bet-1938-table3-silica-gel.csv` — Table III, seven gases on one silica gel.
- `data/bet-1938-table4-so2-silica-gel.csv` — Table IV, SO₂ at six temperatures.
- `data/bet-1938-table5-charcoal.csv` — Table V, eight isotherms on charcoal.
- `data/bet-1938-printed-claims.csv` — every scalar in prose, in a footnote or in
  a figure's typeset constant block that the page checks (91 entries).

**No other page's dataset is loaded**, so none of the cross-page reconciliation
obligations apply. **No curve is digitised, no figure is reproduced, and no page
image is committed anywhere.** Where a figure is used it is a block of
hand-lettered *constants* sitting inside the plot frame — Fig. 3 because it is
the only place in the paper where $v_\mathrm{m}$ and $c$ are printed for the
same isotherm, Fig. 4 because it is the only place a $v_\mathrm{m}$ is printed as
*calculated* rather than fitted.

## What the page establishes

**1. The rival theory fails on its own constants — factor 96.8.** The paper's
chain is reproduced input by input: $8.357-8.823 = -0.466$ exactly,
$d = -0.466/(1+11.1\times0.029) = -0.352523$ against the printed $-0.35$,
$C = d\,(\alpha/r^3) = -0.010223$ against the printed $-0.01$, and
$K_1 = C^2 = 1.045\times10^{-4}$ against the printed $10^{-4}$. Bradley's fitted
$k = 0.615$ needs $C = 0.989$ instead — **96.8 times larger in the dipole ratio
and 9361 times in the energy it transmits**. His $k$ is also outside the range
$k \le 1/2$ in which the polarization recursion has a decaying solution at all.

**2. …and it fails without using the rival's closed form.** Eq. (1b) is quoted by
BET from a paper nobody here has read, so nothing load-bearing rests on it: the
recursion $\mu_i = k(\mu_{i-1}+\mu_{i+1})$ is assembled and the decay ratio
*measured* off the solution. At the computed $k$ it returns eq. (1b) to
$8.5\times10^{-14}$ and is independent of where the ladder is truncated (spread
$0$ over $N = 100$–$600$). At Bradley's $k$ the spread is **4.48** and the
moments change sign along the ladder — there is no decaying film to find.

**3. The BET evidence, split into what was fitted and what was not.** The
isotherm agreement is a fit and is labelled one everywhere; the page never
reports it as a metric. What is *not* a fit is that seven gases, each with its
own fitted $v_\mathrm{m}$, must give the same surface once converted through
cross-sections computed from bulk densities. **A 46.9 % spread in
$v_\mathrm{m}$ collapses to 10.6 % in area** on silica gel (factor 4.42) and
18.6 % → 7.9 % on charcoal. The null baseline is printed beside the result every
time. The polarization theory cannot be asked this question: it has no
$v_\mathrm{m}$.

**4. And the agreement that was mostly guaranteed — from above.** The paper's
best-looking corroboration is $v_\mathrm{m}$ against point B, claimed on journal
page 315 for the **twelve** isotherms of Table I ("the two seldom differing by as
much as 10 %") and extended here to the eighteen non-butane pairs of Tables I
and III. True — and much weaker than it looks. This page derives in closed form,
and confirms by an independent root-find, that the BET inflection point can never
exceed

$$\frac{v(x_\mathrm{infl})}{v_\mathrm{m}} \;\le\; \frac{2}{\sqrt3} = 1.1547\ldots
\qquad\text{at } c = 27+15\sqrt3,\;\; x = 3\sqrt3-5,$$

**for every $c$ whatever**. So an eye-read landmark at or below the inflection
cannot be *high* by more than 15.5 %. **The bound is one-sided** — nothing stops
such a landmark sitting arbitrarily *low*, and on this paper's own Table III
butane isotherm point B is 51.7 % below $v_\mathrm{m}$ (point B / $v_\mathrm{m}$
= 0.4828). Over the eighteen non-butane pairs mean(inflection / point B) =
1.148812, i.e. the inflection sits 14.881 % **above** point B, equivalently point
B sits 12.954 % below the inflection. That ceiling is not in the 1938 paper and
is the page's own result.

## One printed statement does not survive

Proved in the F2.3 order — **pin what is not free first**. The paper never
defines its "±". Table II's four gases pin it exactly, and the two halves are
*not* rounded alike: the printed 840 ± 50, 650 ± 55, 1460 ± 120 and 1900 ± 30 are
the **midrange rounded to the nearest 10 cal** (4 of 4 exact, worst residual
4 cal) and the **half-range rounded to the nearest 5** (4 of 4; argon's 55 is not
a multiple of 10, and CO₂'s 122.5 is an exact tie printed as 120). Applied to
Table I:

| basis | midrange | half-range | rounds to | printed |
|---|---|---|---|---|
| all twelve rows, as the sentence says | 824.5 | 86.5 | 820 ± 85 (± 90 on a nearest-10 half-width) | 840 ± 70 |
| eleven rows, `Cr2O3 gel` (738) omitted | 841.5 | 69.5 | 840 ± 70 (both rules) | 840 ± 70 |

Verbatim on the crop, and the ellipsis spans a **sentence boundary**: "For
nitrogen $E_1-E_\mathrm{L}$ is uniformly 840 ± 70 cal. Since $E_\mathrm{L}$ is
about 1330 cal., $E_1$ is therefore 2170 ± 70 cal. for nitrogen on all twelve
adsorbents." The first sentence is about "the last column of Table I", which has
twelve rows, so the substance is unchanged — but the page is reporting an
authors' error and shows the quotation at full length for that reason.

840 ± 70 is the *eleven*-row number. Two rows fall outside the band as printed,
and of the twelve possible single-row omissions that one is the **only** one that
reproduces it — under **both** half-width rounding rules, which the page checks
rather than assumes. $E_1 = 2170 \pm 70$ inherits the same basis. **Reported, not
repaired**; the page does not claim to know which the authors intended.

Four smaller ones, all reported and none repaired: footnote 10's
`(8.357 - 8.823/r^3)` is mis-grouped (it must be $(8.357-8.823)/r^3$, which gives
the $-0.466$ used two lines later exactly); Table IV's $E_1$ column is 5 cal out
on one row, at a rounding half-way point; catalyst 954's nitrogen
$v_\mathrm{m}$ disagrees by 15.6 % between Table I and Fig. 4, with **no sample
mass printed for Table I's row**, so the page prints both readings and adopts
neither; and Table V's two argon rows repeat one $v_\mathrm{m}$ (215.5) and one
solid-packing surface (746) across a 13 K interval, where the only other gas
measured at both temperatures does not.

A fifth is quantitative rather than typographic. Recovering the molecular
cross-sections independently from Tables III and V agrees to 0.6631 % at worst —
but that maximum is **not** the printed rounding, and the page no longer says it
is. $V_\mathrm{STP}$ and $N_\mathrm{A}$ cancel in the ratio, so the rounding
envelope follows from the printed quantisation alone: for the pair that sets the
maximum (CO at −183 °C, liquid packing) it is 0.221 %, a factor 3.0 smaller than
what is observed. Thirteen of the fourteen pairs sit inside their own envelopes;
that one does not, all four of its inputs were re-read at digit scale, and the
cause is **unknown**.

## What pymrm is doing here, honestly

Nothing to the BET equation. Eqs. (26), (28), (A), (B) and (E) are closed forms;
there is no grid, no time step and no transport anywhere, and most of the page
would run with pymrm uninstalled — exactly as on `A1.6` and `A1.1`. `newton` and
`NumJac` earn their place in three narrow ways: a route to the isotherm that
shares no algebra with it (the layer equilibria solved for $\{s_i\}$ and eq. (15)
summed term by term, exercised on **every branch the paper uses** — $n = 1, 2, 5,
6, 7$ and $\infty$ — agreeing to $6.7\times10^{-16}$ and *catching* a mistyped
$(n+1)$); the same treatment applied to the rival's recursion; and root-finds
where a sweep would have been wrong. Sampling the ceiling on a 40-point log grid
instead of root-finding it moves the ceiling only in the **sixth** significant
figure (1.1546960 against 1.1547005, $3.9\times10^{-6}$ relative) but moves
$c^\ast$ by **1.43 %** (52.2335 against 52.9808 — the grid has a node at
52.2335). The maximum is quadratic, so its *location* degrades a thousandfold
faster than its *value*; it is the location a reader would quote, and that is the
break row.

## Caveats worth reading before reusing anything

- **Tier 6 throughout.** Not one number in any table is a measurement; every
  entry is a constant the authors derived from an isotherm. Reproducing one is
  reproduction, never validation.
- **The rival is judged in its opponent's words.** De Boer & Zwicker (1929),
  de Boer (1931, 1932) and Bradley (1936) are read *through* this paper's
  restatement. None is on disk; none was consulted. If you need to cite the
  polarization theory, go to the originals.
- **$c$ is a reconstruction wherever it appears**, inverted from
  $E_1-E_\mathrm{L}$ through footnote 16 — but the inversion is the authors' own
  operation run backwards. Journal page 313 says the energies were obtained
  *from* $c$ ("From $c$ one can obtain an approximate value for
  $E_1-E_\mathrm{L}$"), with footnote 16's prefactor already set to 1, so the
  prefactor is **not measurable anywhere in this paper** and the inversion
  returns the authors' fitted $c$ to the 0.28 % the integer $E$ column costs.
  What the paper does exhibit is one isotherm — nitrogen on catalyst 954 at
  77.3 K — carrying two $c$ a **factor 2.24** apart: 156.7 fitted directly in
  Fig. 3, against 350.36 implied by the 900 cal Fig. 4 carries over from the
  −183 °C fit for its *calculated* curve. §8 applies that factor as a
  sensitivity: the inflection ceiling is untouched (it is a maximum over all $c$,
  so no rescaling could touch it — a structural insensitivity, not evidence), and
  **the single-point rule of §7 does not survive it** — seven of Table I's twelve
  rows then breach the paper's 5 % and the smallest $c$ falls to 27.58, below the
  threshold 38. The page claims §7 on the paper's own numbers and says plainly
  that the margin is not robust.
- **The molecular cross-sections behind Tables III and V are recovered**, not
  read: Emmett & Brunauer (1937) is not on disk. The recovery is checked across
  the two tables independently and agrees to 0.66 % (see the fifth defect above
  for what that number is and is not). The circularity risk is bounded with a
  measurement rather than an argument: a σ set tuned to flatten the silica-gel
  areas exactly leaves the charcoal areas spread by 14.38 % against 7.37 % for
  the printed set on the same six rows — a factor 1.95, real but modest — and on
  the butane-inclusive seven-row basis the test has no power at all (31.28 %
  against 31.04 %). Both are printed.
- **Scoped out, not closed with outside data:** Fig. 6 (needs vapour pressures
  the paper does not print) and the Darco statement (refers to an isotherm that
  appears nowhere in the paper). Footnote 18's eqs. (C) and (D) are transcribed
  and deliberately not exercised, because the paper never evaluates them either.
- **This page is not a rate.** Adsorption breakthrough and the linear-driving-force
  law are `J1.5`; a BET equilibrium drops into that page's isotherm slot.
