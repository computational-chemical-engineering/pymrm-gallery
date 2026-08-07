# A2.4 — Tanks in series

**Catalog ID:** `A2.4` · **Section:** A · **Tier:** T0 · **Priority:** P1 ·
**Structures:** `S1`, `S3`, `S5` · **Data tier:** 6 (no measurement exists)

Built from Octave Levenspiel, *Chemical Reaction Engineering*, **3rd edition,
John Wiley & Sons (1999), ISBN 0-471-25424-X** — **Chapter 14 "The
Tanks-In-Series Model", sections 14.1 and 14.2**, book pp. 321–338, plus **book
p. 303 in Chapter 13**, which Fig. 14.7's own annotation names.

> **Origin not consulted.** The RTD of $N$ equal tanks, its means and variances,
> and the $F$ curve of Fig. 14.7 are attributed *by this book* to **MacMullin,
> R. B. and Weber, M., Jr., Trans. AIChE 31, 409 (1935)**. That paper is not on
> disk and was not read. Every equation on the page was read from the 1999
> monograph, on 600 dpi renders — the scan's native resolution. Nothing came from
> the text layer, which is excellent for prose and destroys every equation and
> subscript in the chapter (`pdftotext` returns Eq. 3 as `7. = (-&) NN e-tNli`).

## What the page shows

The tanks-in-series model has one adjustable number, and Chapter 14 fits it three
times. Reproducing an RTD from an $N$ fitted to that same RTD would be a goodness
of fit, so the page does not do it. It fits $N$ where Levenspiel fits it — from
**two moments**, in worked example E14.4 — and then asks what the fit is worth on
quantities held out of it.

**The headline is a limit, not an agreement.** E14.4 reads $\bar t = 60$ s and
$\sigma^2 = 900$ s² off a pair of tracer curves and concludes "4 tanks". A chain
of five *unequal* tanks of 9, 9, 9, 9 and 24 s has

$$\textstyle\sum t_i = 60\ \text{s},\qquad \sum t_i^2 = 4(81)+576 = 900\ \text{s}^2$$

— exactly the two numbers the method uses. Three things have to be said carefully
about what follows, and the page says all three:

- the difference is **23.2 % in outlet concentration**, first order, and only
  **0.16 %** in *conversion* (99.334 % against 99.488 %). Over the whole sweep
  the relative gap in conversion never exceeds 0.38 %;
- 23.2 % is the value at $k\bar t = 10$, and it is **monotone** — 1.53 % at
  $k\bar t = 2$, 69.7 % at $k\bar t = 50$ — so it is where the sweep stops, not a
  maximum over it;
- it is **one chain**, hence a lower bound. Over *every* chain matching the same
  two moments the supremum is closed form (the extremal chain has at most two
  distinct tank sizes, by a two-line Lagrange argument, and the bound is
  $\ln(1+k\sqrt S) + k(T-\sqrt S)$): **83.1 %** at $k\bar t = 10$, approached by
  one 30 s tank followed by 30 s of plug flow. E14.4's two numbers pin the outlet
  concentration only to within a factor of 5.9 there.

The exit curves differ by 6.9 % of the peak, and the second-order gap is 6.1 %.
Two moments do not determine a conversion, and the page measures by how much.

## Five printed defects, each proved from the book's own arithmetic

0. **An inverted modelling rule — the consequential one.** Fig. 14.7 (p. 327)
   carries, inside its own axes box, *"When $N > 50$ the curve becomes
   symmetrical in which case use fig. 13-11 with $N = \tfrac12(\mathbf{D}/uL)$"*.
   The ratio is upside down. The book proves it against itself twice: its own
   Fig. 13.12 (p. 303) prints $\sigma^2_\theta = 2(\mathbf{D}/uL)$ and Ch. 14
   prints $\sigma^2_\theta = 1/N$, which give $N = \tfrac12(uL/\mathbf{D})$; and
   under that corrected form the annotation's threshold $N>50$ is exactly
   $\mathbf{D}/uL < 0.01$, the printed validity limit of the very figure it
   points at. Read literally the rule sends a reader at $N = 64$ to
   $\mathrm{Pe} = 0.0078$ instead of 128 — 138 % wrong in conversion at
   $k\bar t = 2$, against 0.053 % de-inverted. Reported, not repaired.
1. **A factor of two.** Section 14.2 (p. 328) prints, for equal final conversion,
   $V_{N\,\text{tanks}}/V_p = 1 + k\bar t_i = 1 + k\bar t/(2N)$. With the book's
   own $\bar t = N\bar t_i$ those differ by two. Expanding the book's own Eq. 9
   against plug flow gives $x/\ln(1+x) = 1 + x/2 + \cdots$ with $x = k\bar t_i$,
   so the **middle** expression is the wrong one; the right-hand form and the
   companion equal-volume formula $1+(k\bar t)^2/(2N)$ are correct as printed.
2. **A transposed digit.** Fig. E14.4b (p. 334) labels its top axis tick
   0.105 s⁻¹, on an axis ticked 0, 0.005, 0.010 at equal spacing. The book's own
   printed $E(t) = 3.2922\times10^{-6}t^3e^{-0.0667t}$ peaks at **0.0149138** s⁻¹
   (at $t = 44.9775$ s), so the tick is 7.04× the whole curve. It should read
   0.015. Note which coefficients: the *unrounded* 4⁴/(3!·60⁴) and 4/60 give
   0.0149361 and a ratio of 7.03, and the tick belongs to the printed curve.
3. **"$L = 272$ moles."** E14.2 (p. 331). The quantity is a distance along the
   Ohio River, and the Comment two lines below reads "$L \le 272$ **miles**". The
   value 272 is exactly right; only the unit is a slip.
4. **"$\Delta(\sigma^2) = 1000 - 100 = 900$ s."** E14.4 (p. 333). A variance
   carrying the unit of a time — the same class as (3), and the value is right.

A sixth item is a typesetting ambiguity rather than a defect: Eq. 11's prefactor
denominator is set as $(N-1)!\,\bar t_N$ with $N$ **lowered**, where $\bar t^{\,N}$
is required for the kernel to normalise. The book prints the raised form itself
two pages later in E14.4, which settles the reading.

None of the six is repaired in the transcription; the CSV carries the glyphs as
printed.

## What pymrm adds

- **The identity, and the correction of the chapter's own conversion rule.** $N$
  equal tanks in series *is* `construct_convflux_upwind` + `construct_div` on
  $N$ cells, exactly — so the model's numerical diffusion is $u L/(2N)$ and the
  equivalent Péclet number is $\mathrm{Pe} = 2N$, which is Fig. 14.7's printed
  annotation the right way up (defect 0 above). With the rule de-inverted, the
  two models' first-order conversions agree to 0.71 % at worst for $N\ge16$,
  $k\bar t\le2$, and disagree by 44 % at $N=2$, $k\bar t=5$. Matching on the
  closed-vessel variance or on the conversion instead of on the small-deviation
  variance shifts Pe by an $O(1)$ offset (1.01 and 1.34 at $N=64$), so the three
  rules agree to $O(1/N)$ relatively and the page says so.
- **A caution about pymrm itself.** There is no pure-outflow boundary condition
  in the library. The nearest, zero-gradient, reconstructs the exit face to second
  order — right for a discretised PDE, wrong when the cell count is a physical
  parameter. It costs **5.9 % at $N=2$ and 3.8 % at $N=4$** and decays as
  $N^{-0.87}$, i.e. it is indistinguishable from ordinary discretisation error and
  is not it. The page suppresses the boundary flux and adds the outflow as a
  `construct_coefficient_matrix` sink on the last tank.

## Reproduce

```bash
cd queue_cases/A2.4/page
python build_page.py          # regenerate index.ipynb from source
jupyter nbconvert --to notebook --execute --inplace index.ipynb
```

Runtime ≈ 8 s, all of it in the dispersion root-finds. Two consecutive executions
give identical content and an identical `agreement.json`; nothing here is
stochastic and the seed is set anyway.

## What this page cannot conclude

- **Nothing is validated against experiment.** Chapter 14 contains no measurement
  of any kind. Every comparison is against the book's arithmetic or an analytical
  limit.
- **Independence is assumed, not tested.** Levenspiel's own footnote flags laminar
  flow as a case where it fails; the additivity check would pass anyway.
- **The equal-tank assumption cuts both ways.** The page cannot say the vessel *is*
  four equal tanks either — only that the two chains are indistinguishable by the
  two-moment method.
- **Fig. 14.3's $E_{\theta,\text{inf}}$ annotation is not settled.** The two
  inflection points have different heights (0.4267 and 0.6943 of the maximum at
  $N=4$) and the figure draws one level. Both are reported and neither is scored.
- **E14.3 is scoped out.** Its data live only in the curve of Fig. E14.3a, and no
  curve is traced anywhere on this page.
- **26 of the 81 metrics are below `check_agreement.py`'s `ABS_FLOOR`** and are
  outside CI entirely — several of them because they reproduce a printed integer
  exactly. Each is named on the page with an above-floor companion that the
  companion string must actually *name*, asserted in code; an earlier draft had
  four whose companion strings named no metric, so their assert ran over an empty
  set and proved nothing.
- **The break table's associations are declared, not measured.** Each row
  recomputes the one metric it is pinned to, and that row's baseline is asserted
  to equal the reported value; the wider mover list is a judgement with every
  name checked against `agreement.json`. The coverage map is built from those
  lists rather than written by hand.
