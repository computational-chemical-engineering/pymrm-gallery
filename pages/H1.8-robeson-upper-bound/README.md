# H1.8 — The Robeson upper bound, revisited

**Source.** Lloyd M. Robeson, *"The upper bound revisited"*, Journal of Membrane
Science **320**(1–2) 390–400 (2008), doi:`10.1016/j.memsci.2008.04.030`, Lehigh
University. Born-digital PDF; every numeric cell used on the page was read a
second time from cropped 300 dpi renders at digit scale and compared against the
text layer. The two readings agree on all 382 **numeric cells** (44 + 234 + 104,
counted in the notebook; cells, not numerals — the digit count is several times
larger). The hazard specific to this paper is comma thousands separators inside
table cells (`1,073,700`, `30,967,000`, `5,369,140`).

**What the page is.** A correlation audit of $P_i = k\,\alpha_{ij}^{\,n}$. There
is no PDE, no grid, and no pymrm operator in it — `structures` and `pymrm_api`
are both empty, deliberately, and the page says so in *PyMRM implementation*.

## Fit, test, or neither

**Neither.** Robeson's $(k,n)$ are not a fit: p. 391 states the line is
*"determined empirically ('by eye')"*. So nothing here is a goodness-of-fit
statement, and the 19 tabulated points above the line are not fit residuals —
they are what "by eye" amounts to, measured. Three labels are used throughout:

- **Transcription** — Tables 1–13 as printed.
- **Internal consistency** — V1 and V2 below. These test *this page's*
  arithmetic and transcription. They cannot test Robeson's, because Table 12 is
  an input on both sides. **V2 is a consistency identity, not a second route**:
  its residual is the difference of the two corresponding V1 residuals, and the
  notebook measures that rather than asserting the opposite.
- **Independent estimate** — V3's covering-line envelope, a different objective
  on a subset of the data, offered as a diagnostic and **not** as a correction.
  It is nevertheless the page's **only genuine second route** to Result 1's
  shift, because it is the only computation here that does not read the present
  bound out of Table 12. **Its scope is the PRESENT line only**: both routes
  read the prior block from the same 2008 restatement, and both are evaluated
  at the same reference selectivity $\alpha_0$ drawn from Tables 1–11, so
  nothing on this page gives the prior block or $\alpha_0$ a second route.

**Nothing on this page is validation against measurement.** Every number is
Robeson's: his compiled points, his hand-drawn lines, his derived tables.

## Headline results

| | |
|---|---|
| pairs in Table 12 / pairs supporting a prior-vs-present comparison | **13 / 9** |
| front-factor share of the shift at each pair's own geometric-mean $\alpha$ | median **0.98**, range 0.59–1.39 |
| … the same share on the four references *Robeson's own Table 13* supplies | **0.56–1.03** |
| pairs whose two bounds are near-parallel (spread/shift < 0.20) | **6 of 9** |
| paper's three *"significant"* pairs recovered by selectivity gain at fixed $P$ | **3 of 3**, gap 0.111 decades at $P_{gm}$ |
| … the same, at every evaluation point inside the range all nine pairs share | **3 of 3**; gap 0.011–0.111 decades |
| second route (envelope present line): shift residual, 7 well-levered pairs, **at each pair's geometric-mean $\alpha$** | **≤ 0.086 decades**, sign **9 of 9** |
| … the same residual and sign count at the *worst* of four in-data references | **0.543 decades** (at each pair's $\alpha_{\max}$); sign **8 of 9** (at each pair's $\alpha_{\min}$) |
| … the same route on Result 3's ranking, at each pair's $P_{gm}$ | **3 of 3** |
| … swept over the same seven references the *primary* route is held to | **2 of 3** at two of them, against **3 of 3** at all seven for the primary route |
| … the same route's front-factor share | median **0.43**, range **−0.03 to 1.55**, *not* reproduced |
| … by $\Delta\log k$ / by permeability gain at fixed $\alpha$ | **2 of 3** / **2 of 3** |
| pairs whose bounds cross inside their own tabulated $\alpha$ range | **1** — H2/CO2 at $\alpha^\star = 37.47$ |
| 2008 H2/CO2 bound relative to the 1994 one at $\alpha = 100.9$ | **−30.4 %** (below it) |
| tabulated "close to the bound" points lying **above** the bound | **19 of 117** (16.2 %), $r$ up to 3.49 |
| Tables 13a/13b reproduced from Table 12 | median 0.017 %, **25 of 26** inside a generous printing band |
| V2 shift cross-check, four printed comparisons | agree to **4.9e-4 decades** |

## The four claims, in one paragraph each

**Nine, not thirteen.** H2/O2 and He/O2 have a prior bound and `NA` present;
CO2/N2 and N2/CH4 have `NA` prior and a present bound. The notebook asserts
that the set of pairs Robeson's own Conclusions classify equals those nine
exactly, and that the eleven data tables match the eleven present bounds — so
the structure is checked, not assumed.

**"Primarily the front factor" needs a reference.** $\Delta\log P(\alpha) =
\Delta\log k + \Delta n\log\alpha$, so $\Delta\log k$ *is* the shift at
$\alpha=1$ — outside every dataset in the paper. The reference-free version of
the claim is that the two lines are near-parallel, which holds for six of nine
pairs and fails for **three** — H2/N2 (0.34), H2/CH4 (0.93) and H2/CO2 (2.05).
Only one of the three is mitigated by the paper: for H2/CO2 it says the shift is
*"primarily a slight slope change"*. For H2/CH4 and H2/N2 it says only that the
shift is modest, and attributes no slope change — searched, every occurrence of
`slope` in the full text. For He/H2, the paper's headline pair, it is airtight:
over the $\alpha$ range its own twelve points span, the slope term contributes
at most 3.0 % of a 1.83-decade shift.

**$\Delta\log k$ is the wrong number to rank pairs by.** $|n|$ runs from 0.79
to 5.8 in Table 12, so equal front-factor moves buy wildly unequal separation.
Ranked by selectivity gain at fixed permeability, the paper's three
*"significant"* pairs are the top three with a 0.111-decade gap; ranked by
$\Delta\log k$ or by permeability gain, CO2/CH4 — which the paper calls
*modest* — interleaves with them. Both of Robeson's statements survive: the
abstract's is about the **shape** of the move, the Conclusions' ranking is
about its **size**.

**The H2/CO2 bounds cross at $\alpha = 37.47$**, inside the range Table 10
spans, and Table 10 has **exactly one** point out there — at $\alpha = 100.9$,
where the next tabulated value down is 23.1. The crossing itself is a property
of four printed numbers; the *"inside its own range"* half of the claim rests
on that single point, and a break row that mis-reads its decimal point takes
the result to zero. Both halves are stated wherever the result appears. The paper prints the cause (a slope
change, and a warning that the low-permeability end may be skewed) but not the
consequence. Searched for any statement of a crossing: the words
cross/crossing/intersect occur only inside the polymer name *crosslinked* and
in two reference titles.

## Validation

- **V1** — Tables 13a/13b reproduce from Table 12. Three nested checks: the
  $(M_j/M_i)^{1/2}$ column to 1.5e-4 from **integer** molecular weights (and
  only to 3.9e-3 from IUPAC weights, which *identifies* the convention); the
  product column to 1.6e-3; the permeability column to a median 0.017 % over
  26 rows, with propagated printing bands computed two ways.
- **V2 — a consistency identity, NOT a second route.** Tables 13a and 13b each
  evaluate He/H2 and H2/CO2 at one selectivity on **both** bounds, and the
  ratio of the printed permeabilities agrees with the decomposition to
  4.9e-4 decades. But Robeson's transition permeability *is* $k\alpha^n$ on the
  same two lines V1 tests, so that residual **is the difference of the two
  corresponding V1 residuals** — the notebook prints both and they agree to
  4e-16 decades. V2 cannot fail unless V1 fails on those four rows. Its real
  power is over *this page's* decomposition code: the natural-log break row
  moves it and moves nothing else. The same cell shows that the obvious second
  candidate — the shift read off the data as a ratio of geometric-mean
  excursions — is also an identity, equal to the decomposition to 5e-16
  decades on all nine pairs. **Given Table 12, Result 1's decomposition is a
  rearrangement and no arithmetic on Table 12 can corroborate it.**
- **V3 — the second route, and the only one.** A covering-line envelope
  re-estimated from Tables 1–11 by exact active-set enumeration (no optimiser,
  no seed). Substituting it for the printed present line gives an
  envelope-implied shift whose residual against the decomposition is exactly
  V3's offset at $\alpha_0$. **Scope: the present line only** — both routes
  read the prior block from the same restatement and both use the same
  $\alpha_0$.

  *At each pair's geometric-mean selectivity*, the two routes agree on the
  **sign** for 9 of 9 pairs and on the **magnitude** to 0.086 decades on the
  seven pairs whose points span more than 1.5 decades of selectivity (0.47 on
  O2/N2, whose points span 0.82).

  **Those two numbers are themselves reference-dependent, and the page sweeps
  them rather than quoting one.** The residual carries $\log\alpha_0$
  explicitly, so it moves for exactly the reason the front-factor share does.
  Across each pair's own minimum, geometric-mean, median and maximum tabulated
  selectivity — all four inside the pair's own data — the wide-span worst case
  runs **0.086 → 0.543 decades** (a factor 6.3) and the median over all nine
  **0.038 → 0.315**. The geometric mean, the reference used, is the kindest of
  the four. The sign count is **8 of 9 at each pair's $\alpha_{\min}$**, where
  He/N2's envelope-implied shift is −0.012 decades against a decomposed +0.202.
  Quoting "≤ 0.086 decades" without its reference is the same error this page
  exists to warn about, one section later.

  The two routes do **not** agree on the front-factor share — median 0.43
  against 0.98, range −0.03 to 1.55 — which is this page's central caveat
  arriving from the other side. The −0.03 (He/N2) is explained in words on the
  page: the envelope's front factor moved *down* (12,239 barrers against the
  prior 12,500) while the bound moved *up* by 0.299 decades, so the whole shift
  is carried by the slope term.

  Result 3's ranking also gets a second route out of the same substitution, and
  **it corroborates less strongly than the primary route, which the page states
  rather than glossing**: 3 of 3 at the same reference and at five of the seven
  references the primary route is swept over, but **2 of 3 at each pair's own
  $P_{\text{median}}$ and $P_{\max}$**, where He/H2 falls out of the top three
  — against 3 of 3 at all seven for the primary route. Weak where the lever arm
  is short, and the measurement says which pairs those are.

  What the selection of the points does *not* limit: the route's detecting
  power. $k_{\text{env}}$ never touches Table 12, so the residual responds to a
  vertical error in the printed present line exactly 1:1 — break row (DD)
  divides $k$(He/CH4) by ten and moves that pair's residual 0.07341 → 1.07341,
  measured.
- **V4** — fourteen break rows, three of them perturbing a single cell of
  Tables 1–11 (an earlier version perturbed only Table 12, leaving the 117
  points with no break row at all). The coverage map is **built from the
  measured mover list of each of the 58 reported metrics**, not written by
  hand, and the clean row recomputes every metric by a separate code path and
  is asserted equal to the reported values (bit-identical on every key).
  That is a second implementation of the page's **assembly**: `rerun_all`
  reads no narrative result frame, but it calls the same shared primitives
  — `P_bound`, `alpha_bound`, `excursion`, `ulp`, `sqrtM` and `envelope_fit`,
  the most intricate code here — which are written **once**, not twice. The
  page says so rather than claiming a second implementation of everything.

## Printed defects, reported and never repaired

1. Conclusions, p. 399: *"Significant shifts … were observed for He/CH4,
   He/CO2 and He/H2. The shifts in the upper bound for He/CO2, He/CO2 and
   He/H2 were due to aliphatic fluorocarbon polymers"* [sic] — He/CO2 twice,
   He/CH4 dropped. Table 12 confirms He/CH4 belongs (5,002 → 19,800).
2. p. 392: *"the one data point above the present upper bound"* for O2/N2;
   Table 1 has two above the printed line.
3. Table 13a's O2/N2 row misses $k\alpha^n$ by −0.612 %, against a worst-case
   propagated band of ±0.279 %. Backing the printed permeability out needs
   $\alpha = 1.4714$ where the table prints 1.473.
4. Table 2 prints `Polypyrrole` where Tables 1/5/6 print `Polypyrrolone` for
   the same material from the same reference; Table 12's prior block prints
   `He/CO2` without the final subscript.
5. The author repairs a compiled CH4 permeability in prose (80.2 → 180.2) but
   not in Table 9, which still carries $\alpha = 1.9$, consistent with 80.2 and
   not with 180.2. That point is the largest excursion in its table.

## Out of scope, with the search that establishes it

Fig. 12's $-1/n$ versus gas-diameter correlation and Freeman's eq. (2) need
inputs the 2008 paper does not print. Searched the full text layer for the
string `diameter`: **17 occurrences** (8 singular, 9 plural, counted as
occurrences in the `pdftotext -layout` extraction of all eleven pages, not as
matching lines), every one of them prose, a caption or eq. (1); the only numeric
diameters anywhere in the paper are PTMSP's pore size, 0.9–1.2 nm. No verdict is
offered on either. No figure was digitised anywhere on this page and none needed
to be.

## Reuse

`excursion(P, alpha, k, n)` returns $r$; $r>1$ means above the bound. Take
$(k,n)$ from the `present` block. Four cautions, in the order they bite: nine
pairs not thirteen; never quote a front-factor share without its reference
selectivity (**0.56–1.03 across Robeson's own Table 13 references, 0.59–1.39
across the nine pairs at this page's reference — the two ranges are not
interchangeable**); rank by selectivity gain, not by $k$; and the bound applies to
**homogeneous polymer films only** — the paper explicitly excludes
heterogeneous, surface-modified, mixed-matrix, carbon molecular-sieve and
thermally-rearranged membranes, so a Pd (`H1.1`) or zeolite (`H1.9`) membrane
beating it is not news.

The transport model that produces a $(P,\alpha)$ pair in the first place is
`H1.7`; that page *is* a pymrm page. This one supplies the empirical ceiling
it lives under, and does not derive it.

## Files

```
build_page.py   -> index.ipynb        (executed clean; two runs give identical content
                                       and an identical agreement.json)
meta.yaml
agreement.json
data/robeson-2008-table12-upper-bounds.csv          + .meta.yaml
data/robeson-2008-near-bound-points.csv             + .meta.yaml
data/robeson-2008-table13-knudsen-transition.csv    + .meta.yaml
data/robeson-2008-printed-claims.csv                + .meta.yaml
```

Runtime ≈ 8.7 s (measured: two consecutive `NotebookClient` executions in
place, 8.66 s and 8.66 s wall including kernel startup; four further runs in
an isolated copy of the directory fall between 8.92 s and 9.26 s). The break
table recomputes all 58 metrics under 14 rows, which is where the runtime
goes.
