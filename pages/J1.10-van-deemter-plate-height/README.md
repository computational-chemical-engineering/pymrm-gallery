# J1.10 — The van Deemter equation

**Catalog ID:** `J1.10` · **Section:** J · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S4` (1-D transient PDE) · **Data tier:** 6 (the source's own
four printed tables, transcribed; **no figure digitised**)

Reproduces J. J. van Deemter, F. J. Zuiderweg and A. Klinkenberg, "Longitudinal
diffusion and resistance to mass transfer as causes of nonideality in
chromatography", *Chemical Engineering Science* **5**(6) 271–289 (1956),
doi:10.1016/0009-2509(56)80003-1 — the paper behind $H = A + B/u + Cu$.

Built from the 300 ppi CCITT-G4 scan (native resolution, confirmed with
`pdfimages -list`). **The text layer was used for nothing numeric**: it
misreads the journal *volume* in the document's very first line, extracting
`Vol 6` where the page prints `Vol. 5`. Every table cell, equation constant and
quoted glyph was read on a 2× digit-scale crop of a 300 ppi render.

**The paper prints four tables, not two** — a correction to this repository's
own papers inventory. Table 1 (p. 277) is the feed-volume experiment, Table 2
(p. 280) is Simpson & Wheaton's ion-exclusion data restated, Table 3 (p. 281)
is the plate heights derived from Table 2, and Table 4 (p. 286) is the
gas–liquid measurement that produces the three-term formula. All four are
transcribed; the two that carry the plate height are 3 and 4.

**No point was taken off any figure**, and the figure consultations that did
happen are counted exactly: five ambiguous readings were settled by arithmetic
constraints the tables impose on themselves (two in Table 2, three in Table 3);
three more, spanning five cells of Table 4, were settled with Figs. 11 and 13 —
and the page shows that all three follow from inside Table 4 alone, from its
ascending velocity column and from eq. (55). Two printed *curve labels* on
Fig. 8 ($n = 240$, $n = 277$) are read as text, and used only to identify which
of Table 3's two printed endpoints came from which route.

## Fit, not test — where it matters

**Eqs. (54) and (55) were fitted by the authors to the twelve Table 4 points
this page compares them against**, and the paper says so: "the values of the
constants in (54) are highly uncertain, owing to the relatively small number of
experimental data that were available". Every residual against Table 4 here is
a *goodness of fit*, is called one, and is reported beside a null baseline (the
best constant $H$, which the fits beat by a factor 3.3 to 6.2 in rmse).

The same applies to $\lambda$ and $D_{II}$: the authors extracted them *from*
Table 3, so refitting Table 3 reproduces their extraction, not an independent
measurement of either constant.

## What the page shows

- **The 1956 algebra is sound, all four steps.** (34) → (38), (38) + (45,46,48)
  → (49) and its $C$, (38) + (41,50,51) → (52): each residual is exactly zero.
  The paper names the route for each — which equations go into which — but never
  prints the manipulation, and none of them is checked in the source. A fourth route the
  paper never takes — the Laplace transform of eqs. (29)/(30), whose
  residence-time cumulants give $H = \kappa_2/\kappa_1^2$ — reproduces eq. (38)
  symbolically as well.

- **Eq. (38), computed without eq. (38).** Integrating eqs. (29)/(30) in time
  with pymrm and reading $H$ from the difference of the first two time moments
  at two interior stations agrees to **4.1 × 10⁻⁵**. The two-station difference
  cancels every boundary effect — measured: swapping a Danckwerts inlet for a
  Dirichlet one moves $H$ by 1.5 × 10⁻⁶.

- **A refinement study with a known answer.** First-order upwind adds $uh/2$ of
  numerical dispersion, so switching the limiter off must shift $H$ by *exactly
  $h$*. Measured coefficient: **0.99999**. With the van Leer deferred correction
  the observed order is 3.00; bare, 1.00.

- **$D_{II}$ and $\lambda$ recovered from Table 3.** Refitting eq. (49) on the
  five varying-flow rows gives $\lambda = 2.99$ against the printed 3 (0.24 %)
  and $D_{II} = 1.34\times10^{-10}$ m² s⁻¹ against the printed
  $1.3\times10^{-10}$ (3.7 %). The extraction turns out to be **almost
  independent of $D_I$** — which the paper never prints for that system — moving
  only 1.8 % as $D_I$ is swept 2× and then removed altogether.

- **The film thickness, three times.** 10.06, 10.57 and 8.66 µm from three
  different printed constant sets, against the paper's $\approx 10$, $\approx 10$
  and 9 µm. Eq. (56)'s coefficient comes back at 3.333 × 10⁻⁸ against the printed
  3.3 × 10⁻⁸, and the Sterchamol ratio is *sharper* with this page's $d_f$
  (2.56 × 10⁻⁵) than with the paper's rounded 9 µm (2.37 × 10⁻⁵) against a
  printed 2.5 × 10⁻⁵.

- **Table 1's `calc` column, all eleven values, and what it really pins.** The
  column depends on one number, $v\sqrt n = \Delta S_0/4$, and the eleven printed
  values pin $\Delta S_0$ to **166.0–167.4 ml** — 0.6 to 1.5 % above the paper's
  own "$\sim$ 165 ml", which is written with a tilde and is therefore not a
  defect, merely less precise than its own consequences.

- **Table 3 rebuilt from Table 2 by the paper's two routes.** Median deviation
  2.1 % over eighteen printed endpoints; split **by route**, the height route
  (median 1.6 %) beats the width route (2.2 %), as it must, because Table 2's
  widths are printed as integers and $n$ goes as $1/w^2$. The printed pair is in
  ascending order, so route and endpoint are not the same split: Fig. 8's own
  printed labels ($n = 240$ height, $n = 277$ width) fix low = height for the
  50–100 mesh row, and in one row of nine the two routes invert.

- **A printed equation that cannot be right.** Eq. (21) as printed carries
  $(a+\delta)/\delta$ where the derivation and the paper's *own* eqs. (24) and
  (25) require $(a+\delta)/a$: as printed it gives 2 in the limit eq. (24) fixes
  at 4, and diverges where eq. (25) fixes it at $a+\sqrt{2\pi}$. Eq. (22) is not
  free — it is exactly the inflection condition of eq. (19), verified
  symbolically — so the denominator is the only thing left, and deriving
  $\Delta s$ from eq. (19)'s own inflection point returns the printed term with
  $a$ underneath at residual exactly zero. As printed, the width route to the
  plate number admits no root at all for seven of Table 2's nine conditions. **Reported, not repaired**: quoted verbatim, shown to fail its own
  two limits, and the corrected form used with every dependent number labelled.

- **The velocity convention, stated rather than assumed.** The List of Symbols
  (p. 288) defines $u$ as *interstitial*; Tables 2–4 and Figs. 11–14 are headed
  *superficial*. Fig. 9 converts (abscissa $uF_Id_p$); the gas–liquid section
  does not appear to — the quoted $2\gamma D_I = 8\times10^{-6}$ m² s⁻¹ is
  numerically identical to eq. (54)'s iso-butane $B$ read in mm² s⁻¹, i.e. with
  $F_I = 1$, and implies $\gamma = 0.444$, just below the paper's own stated
  0.5–1.0. $F_I$ is not printed for those columns, so **the ambiguity cannot be
  closed**; the page reports the band of $F_I$ each reading needs (0.50–1.00) and
  notes that $2\lambda d_p$ is immune because the eddy term contains no $u$.

## What pymrm adds

Not the plate height — eq. (38) is exact for this model and the page proves it
twice. What it adds is a number for the paper's own *unquantified* condition.
Eqs. (33)/(34) hold "for locations $z$ much larger than both $2D/u$ and
$F_Iu/\alpha$", and how much larger is never said.

Because the cumulants of a linear column are exactly additive in $z$, $\kappa_3$
and $\kappa_4$ per unit length give the skewness and excess kurtosis of a column
of any length in closed form. **Both of the paper's groups fall as $1/z$; the
skewness falls only as $z^{-1/2}$.** At the ion-exclusion column that produced
Table 3 the groups are 2.4 × 10⁻³ and 3.4 × 10⁻³ — comfortably "$\ll 1$" — and
the elution curve is still skewed by **0.14**. Reaching a skewness of 0.05 needs
a column 8.2× longer (root-found, not swept).

All four cumulants come from the Laplace expansion, so every quoted shape number
is analytic; the pymrm run is the *check* on them, and it is reported with its
caveat. The PDE reproduces $\kappa_3$ to 0.21 % and $\kappa_4$ to 1.5 % at the
reference settings, but neither is a converged discretisation error: the grid
error on the higher cumulants is positive and the time error negative, and the
page prints the $\Delta t$ table in which both deviations pass through zero.
$H$ itself has no such problem — it is grid-limited at observed order 3 and
moves 1.2 × 10⁻⁶ over 8× in $\Delta t$.

## What the page cannot conclude

- **Nothing about whether eq. (53) is right.** Its constants were fitted to the
  points it is compared with. The nearest thing to a test the paper offers is the
  slope ratio — 1.448 predicted from the two printed $KF_I/F_{II}$ against the
  paper's "about 1.5" and eq. (54)'s 1.6 — and that is one number.
- **Nothing about $\gamma$ or $d_f$ as physical quantities.** They are the
  paper's own inversions, redone, and the velocity convention moves them in
  opposite directions: $\gamma$ up by $1/F_I$, $d_f^2$ down by $F_I$, each by up
  to a factor 2. $\lambda$ is untouched by it — the eddy term carries no $u$.
- **Nothing about Simpson & Wheaton's measurements.** Their paper is not on disk
  and was not consulted; Table 2 is van Deemter et al.'s restatement.
- **Nothing about the ion-exclusion column's geometry.** Its length is printed
  nowhere; the 558 mm used here is recovered from Table 3 and is an inference.

## Files

```
index.ipynb   the page (executes clean in ~50 s)
build_page.py regenerates index.ipynb
meta.yaml     metadata, printed defects, inferences, fit-vs-test labelling
data/         four tables + 56 printed scalars, each with a provenance sidecar
```

**Cite the source, not this page.**
