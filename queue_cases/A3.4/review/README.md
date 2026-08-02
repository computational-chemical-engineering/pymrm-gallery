# Review overlays — deliberately not committed

The PNGs in this directory are extraction overlays drawn on the source figure,
so they **are** the copyrighted page image. They are git-ignored.

They exist to be inspected once, by the maintainer, in a private artifact — the
same basis on which every dataset sidecar in this repository states that the
figure image is *not* reproduced and only the extracted numbers are published as
facts. Delete them once the review is answered.

| file | what to look at |
| --- | --- |
| `A3.4-fig2-recentred-overlay.png` | Figure 2, all 79 glyphs. **Green** = shape-fitted centre with the fitted outline; **orange x** joined by a line = where the 2026-07-30 crosshair was; **red** = found on 2026-08-02 |
| `A3.4-fig2-recentred-audit.png` | the same 79 at 2×, one panel each, numbered as `marker_id` in the CSV |
| `A3.4-fig2-axis-calibration.png` | unchanged from 2026-07-30; the three printed curves recomputed onto the drawn ones |
| `A3.4-fig3-extraction-overlay.png` | Figure 3, the 182 glyphs of the shipped erasure-band setting, outlines only |
| `A3.4-fig3-candidate-audit.png` | the same 182 at 2×, **worst fit first** — the failures are at the top left |

---

# Technique: locating a marker by fitting its shape

*Written for promotion into `docs/handoff.md`. It grew out of a maintainer note
on this case — "mostly symbols are identified well, but often the marker is off
centre; it might be helpful to try to fit the shape and then take the center of
the shape" — and it replaced the ink-density method `C2.1`, `F1.4` and the first
`A3.4` pass all used. The code is in the session scratchpad as `glyphfit.py`,
`setup2.py`, `setup3.py`, `overlay.py`; worth promoting into `scripts/` now that
a fifth figure has needed it.*

## The problem with a density method

Every detector this repository has used so far scores a *position*: box-filtered
ink density, minus whatever a locally straight structure could explain, with
non-maximum suppression. It answers "is there compact ink here?" and returns the
argmax. Three things follow, and all three were measured on `A3.4` Figure 2:

1. **The argmax of ink density is not the centre of the glyph.** For an open
   circle the density peak sits wherever the stroke happens to be thickest; for
   a triangle it drifts towards the base; for a glyph clipped by a curve it
   drifts away from the clip. Median displacement here was 3.3 px on a ~22 px
   glyph, with a 90th percentile of 5.7 px — measured on the 65 markers that
   were actually re-centred. Averaging in the 16 that fell back to their
   crosshair, and so have zero displacement by construction, understates it as
   2.3 px; do not do that.
2. **A density peak has no size and no shape,** so nothing stops the detector
   firing twice inside one glyph: 11 of the 81 crosshairs on this figure were a
   second or third mark on a glyph another crosshair already had. Do not repeat
   the explanation the first draft of this note gave — that "non-maximum
   suppression at 8 px does not prevent it". Eight of the twelve merged pairs
   were 1.4 to 5.0 px apart, and nothing 1.4 px apart survives an 8 px
   suppression radius, so either that pass did not apply NMS as documented or
   those crosshairs came from the hand audit. The **count** is a measurement;
   the **cause** is not known, and a shape fit exposes the double count either
   way.
3. **It cannot separate merged glyphs,** because two overlapping markers make
   one density plateau with one argmax.

## The method

Score a *template*, not a position. For a glyph of shape `s` and bounding-box
size `d`, placed with its symmetry centre at sub-pixel `(y, x)`:

```
score = <ink over the glyph body>  -  <ink over a thin outer margin annulus>
```

both averages weighted by a `care` mask that is zero wherever the ink is known
not to be a marker. Report the fitted `(y, x)`.

Four details, each of which cost a debugging cycle:

- **The reference point is the symmetry centre of the bounding box, not the area
  centroid.** A drafter centres a triangle on its data point by eye, and so does
  every plotting library; the area centroid of a triangle is a third of the way
  down and would be wrong by ~3 px.
- **The margin annulus has to be WIDE.** It ran from 1.16 to 1.62 times the
  glyph size here. The first version derived it by dilating the body by a couple
  of pixels, giving a 1-px ring — and with a 1-px ring the size is unconstrained,
  so every one of the 130 Figure 3 fits collapsed onto the smallest template
  that would sit inside a blob and scored 1.00. All-perfect scores are the
  symptom; a diagnosis of "every marker is a small triangle" is the tell.
- **Optimise properly.** Discrete search over (shape, size, integer offset) using
  `sliding_window_view` and `tensordot`, then Nelder–Mead on `(dy, dx, log d)`
  with an **explicit `initial_simplex`**. SciPy's default simplex perturbs each
  coordinate by 5 % *of its value*, which is 0.00025 when the start is 0, so a
  polish started at the origin converges instantly and reports zero displacement.
- **Cache the rasterised templates** on a quantised `(shape, size, stroke, dy,
  dx)` key. Rasterisation is the only expensive step; 260 markers fit in ~30 s.

Detection is then **matching pursuit in rounds**: take every peak of the coverage
map above threshold, fit it, keep it if it passes the gates, subtract the ink it
explains, recompute. Merged glyphs are explained one at a time. Round-based
batching is ~50× cheaper than one-at-a-time and gave the same answer here.

## Two gates, and why the obvious one is not enough

`score` normalises by the *visible* support, so it is unbounded above when most
of the template is masked away: a sliver of surviving ink beside an erased curve
scores 1.0. That single defect produced 75 false positives in one pass, all of
them triangles sitting on the Ranz–Marshall dashes. The fix is a second gate:

```
evidence = <cared-for ink under the body> / max(<cared-for body area>, f · <full body area>)
```

with `f ≈ 0.45`. The floor means a template can lose up to 55 % of its support to
masking and still be judged fairly, and no more than that. Set `f = 1` and every
marker lying on a curve is rejected; drop the floor entirely and the sliver
problem returns.

## Erasing what is not a marker: three mechanisms, not one

`docs/handoff.md` already records that curves whose equations are printed can be
erased at their computed position. That is necessary but not sufficient, and the
three mechanisms are worth keeping separate because they fail differently.

1. **Straight ink, by morphology.** Binary opening with a long thin line element,
   unioned over 12–18 orientations. Removes long curves, annotation leader lines
   and the axis frame. A glyph 18–26 px across cannot contain a 37–61 px straight
   run, so glyph ink survives — but a *chain* of merged markers along a drawn
   line can, which is why the length has to be raised when the data crowd the
   curve (37 px on Figure 2, 61 px on Figure 3).
2. **Frame and ticks, by connectivity.** Ink connected to the frame, *clipped to
   a band* of ~30 px from it. Connectivity alone is a trap: the eq. 12 curve runs
   into the right-hand frame, so a plain connected-component erasure deletes the
   curve, everything merged with it and about 40 % of the markers. Clipping to a
   band fixes it. Ticks are ~26 px long and morphology (1) never sees them.
3. **Dashed curves, by computed position, and only by that.** The
   Ranz–Marshall dashes on this figure are 12–20 px of compact ink: they are
   *indistinguishable from a marker by any shape or straightness test*, and the
   only thing that separates them is that you know where the curve is. This is
   the one mechanism that can destroy a real marker, so use the narrowest band
   that works and say so in the sidecar.

**Erase as "don't care", never as "no ink".** Set the mask to zero and let the
score renormalise. Subtracting the ink instead makes every marker that touches a
curve score as a partial glyph, and here that was 25 % of them.

## Measure the improvement, do not assert it

Keep both coordinates in the shipped CSV — old crosshair and new fitted centre,
one row per position ever recorded — so the comparison is a paired sample the
page can compute rather than a claim in a sidecar. On `A3.4` Figure 2 that gave:

| | along Re (columns) | along Sh (rows) |
| --- | --- | --- |
| mean displacement | −0.01 ± 2.90 px | +0.78 ± 2.43 px |
| standard errors from zero | 0.0 | 2.6 |
| in data units | +0.01 ± 1.65 % | −0.51 ± 1.64 % |

on the **65 markers that actually moved** (median distance 3.3 px, 90th
percentile 5.7 px, worst 9.9 px). **Verdict: the column component is consistent
with zero; the row component is a small systematic, not noise.** A component 2.6
standard errors from zero is a bias, and the first draft of this note called it
"random" — the right statement is its size, about half a per cent in Sh, an order
of magnitude below the scatter of the data. Refitting the correlation on exactly
the same 81 markers moved α from 1.168 to 1.152 and β from 0.5550 to 0.5572,
under a tenth of one standard error either way. The off-centre crosshairs were
adding scatter and a bias too small to matter, and nothing on the page turned on
them.

That is a *negative* result for the headline motivation, and it is the reason to
measure rather than assert. The method still earned its place, for the reason
nobody predicted: **two crosshairs on one glyph produce two fits at the same
point**, which is how a shape fit detects a double count. Density thresholding
cannot, at any threshold.

## Two lessons that are not about centring at all

Both came out of adversarial review of this case on 2026-08-02, both are general,
and both are the kind of thing that is invisible until someone attacks the page.

**A control that fails in the direction of your result is a symptom, not a
confirmation.** Figure 3 prints its own correlation inside the plot, so its ink
is an object whose fitted parameters are known in advance — 1.1 and 0.6. Under
the flat calibration it comes back α high and β low (2.4 % and 0.9 % on the
tracer now shipped; the first pass's own tracer said 3.9 % and 1.1 %), and the
first pass wrote in the sidecar that "the axes are therefore calibrated by the
paper's own printed function". But α high and β low were *exactly* the signs of
the disagreement the page was claiming about the markers. Tracing the offset
column by column, instead of averaging it, showed a tilt of about 6 px across
the plot — and that tilt was the 5.6 px left/right decade-tick skew the same
sidecar already recorded as "~0.5 % of systematic axis error" and had never
propagated into β. Rules that follow:

- if a figure prints a curve whose equation is given, **trace it and fit it**;
  the residual is your calibration error, free of charge;
- report the residual **against position**, never as a single mean. A tilt has
  mean zero and destroys the exponent;
- calibrate a log axis from **both** sets of decade ticks separately. If the two
  sides disagree, the render is rotated and the origin is a function of the
  other coordinate. Averaging them is right in the middle and wrong at the ends;
- ship the traced curve as a dataset, so the control runs on the page instead of
  being quoted from a sidecar.

**Sweep the erasure width and ship every marker set.** The band erased around a
computed curve is the only mechanism that can destroy a real marker, and its
half-width is a free parameter. On this figure, sweeping it 0 → 14 px moved the
marker count 210 → 126 and the free α 1.296 → 1.440 — a spread larger than the
regression's own standard error. One setting plus a sentence of prose is not a
disclosure; it is how a sensitivity claim gets made about the wrong knob (the
case file here
said "the fit is insensitive, α moved only 1.371 to 1.376", which was true of the
detection *threshold* and false of the *band* at the same marker count). Put a
`band_hw` column in the CSV, ship the whole sweep, and let the page compute the
systematic. It costs six pursuit runs and about two minutes.

And when you do choose a setting, justify it on something measurable. Here,
turning the band off gains 35 centres that sit a median 3.6 px from the drawn
line and score a median 0.684 against 0.793 for the rest — evidence that they are
surviving line fragments rather than rescued markers. That is a judgement, it is
printed as one, and the sweep is what it costs if it is wrong.

## When to reach for it

Worth it when markers are ~15 px or larger, filled or open, and merge — which is
most dense scatter plots in pre-1980 journals. Not worth it when the markers are
tiny relative to the line width, when a curve can carry the series identity
instead (`F2.3`), or when the model under test reads no column the shape would
supply (`F1.4`). And ask the `F1.4` question first: *which columns does the model
you are testing actually read?* Here it read positions only, so the eleven fitted
shapes are recorded as advisory and no series label is carried.
