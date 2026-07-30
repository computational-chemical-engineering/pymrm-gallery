# E2.1 — Figure 9 extraction, staged for review

**The page does not use any of this.** `queue_cases/E2.1/page/` is built and runs
clean as a tier-6 page validated against the paper's three worked appendices.
What follows is the *optional* upgrade, and it needs one visual judgement that a
maintainer has to make.

No page image is committed. Run `extract_figure9.py` against your own copy of the
paper to regenerate `fig9_overlay.png` and `fig9_contact.png` locally:

```bash
python extract_figure9.py ~/papers/pymrm-gallery/i260028a001.pdf /tmp/e21
```

## What Figure 9 is

Conversion `1 − X_A` (log, 1.0 to 0.04) against `𝒦_m = (1−ε_m)L_m K_r/u_0`
(linear, 0 to 8), for Kobayashi & Arai's ozone-decomposition bed —
`d_t` = 20 cm, `L_m` = 34 cm, `u_mf` = 2.1 cm/s. Three series, distinguished only
by glyph, each fitted by the authors with its own bubble size:

| glyph | `u_0` [cm/s] | `u_0/u_mf` | fitted `d_b` |
|---|---|---|---|
| ⊙ dotted circle | 6.6 | 3.14 | 3.7 cm |
| ● filled circle | 9.9 | 4.71 | 4.2 cm |
| ○ open circle | 13.2, 16.5, 20 | 6.28, 9.55 | 5.0 cm |

This is the paper's sharpest test: each series is at fixed `u_0` and `L_m` with
`𝒦_m` varied by changing `K_r` alone.

## Axis calibration — settled, and not the question

On the 600 dpi render of PDF page 9, sub-image `rows 200:1700, cols 2550:4900`:

```
row = 236.660 − 750.532 · log10(1 − X_A)     residuals ≤ 5 px over 15 printed ticks
col = 764.544 + 121.475 · 𝒦_m               residuals ≤ 2 px over  9 printed ticks
```

The frame corners recover 1.007 and 0.0394 against the nominal 1.0 and 0.04, so
the axis limits are exactly as labelled and the frame is the axis. This part is
not in doubt.

## What the detector finds

37 candidates from a ring template (outer radius 11–15.5 px, hole 4–9.5 px,
centre dot ≤ 2.5 px — measured on an isolated glyph). They split three ways.

**15 clean markers at `𝒦_m` > 1.8.** Isolated, and the three glyphs separate
unambiguously on interior ink: filled has hole ≈ 1.0, dotted has hole 0.08–0.20
with a solid centre, open has hole < 0.14 and no centre.

**6 false positives**, all glyphs from the annotation block ("Calcd. by Eqs. 54
and 56", "d_b = 5.0 cm", "cm"). They are removed by a principled filter rather
than by hand: no marker can sit above the topmost printed model curve, so
anything with `1 − X_A > 1.05 ·` (the `d_b` = 5.0 curve) is flagged
`suspect_text`. The filter catches all six and nothing else.

**16 unresolvable markers at `𝒦_m` < 1.8.** They overlap into a connected black
mass; see patches 0–15 of the contact sheet. This is not a classifier problem —
the ink is merged, and individual centres cannot be located. They would have to
be carried as `resolved = no`, the way `F2.3` carries its origin cluster.

## The result, if it were used

Assigning `d_b` by glyph exactly as the paper does, and computing `1 − X_A` from
equations 49, 50 and 54:

```
  glyph   d_b    𝒦_m    data    model    (model−data)/data
   open   5.0   2.21   0.3913   0.4307    +10.1 %
 filled   4.2   2.36   0.3576   0.3484     −2.6 %
   open   5.0   2.42   0.3856   0.4172     +8.2 %
   open   5.0   2.70   0.3903   0.4012     +2.8 %
 dotted   3.7   2.93   0.2796   0.2621     −6.3 %
   open   5.0   2.94   0.3752   0.3889     +3.6 %
 filled   4.2   3.25   0.3160   0.2975     −5.9 %
 dotted   3.7   3.33   0.2626   0.2425     −7.7 %
   open   5.0   3.76   0.3118   0.3543    +13.7 %
   open   5.0   3.84   0.3542   0.3513     −0.8 %
 dotted   3.7   4.20   0.2122   0.2087     −1.7 %
   open   5.0   4.67   0.3113   0.3244     +4.2 %
 filled   4.2   4.89   0.2127   0.2400    +12.9 %
 dotted   3.7   7.00   0.1350   0.1466     +8.6 %
 dotted   3.7   7.50   0.1451   0.1392     −4.0 %

n = 15, mean absolute deviation 6.2 %, bias +2.3 %
```

**The classification is internally corroborated.** The three glyph groups occupy
three non-overlapping bands — dotted 0.135–0.280, filled 0.213–0.358, open
0.311–0.391 — and never cross. If the glyph classification were wrong the groups
would interleave. Shape and position agree here, which is the opposite of the
`F1.4` and `F2.3` situation.

## The judgement being asked for

1. **Is the right-hand subset publishable on its own?** 15 of ~31 real markers,
   with 16 reported as unresolvable. Precedent exists (`F2.3` carries 20 of 99 as
   `resolved = no`), but that page kept the majority.

2. **Does it establish anything the appendices do not?** `d_b` was fitted per
   series *to these very points*, so 6.2 % measures the authors' own fit residual
   plus digitisation error, not a prediction. The one genuinely predictive claim
   the subset does test is the paper's statement that "the progression of bubble
   sizes is in the right direction" — and the three bands do come out ordered
   3.7 < 4.2 < 5.0 with velocity, with no fitting on our side.

3. **If yes**, `E2.1` gains an experimental tier and the page needs a `The data`
   section rewritten, a second CSV with a sidecar, and one paragraph in
   `What pymrm adds` distinguishing *validated against measurement* from
   *reproduced against the authors' fit*. If no, the page stays exactly as built.

Recommendation: **weak yes on (1) and (2)**, on the strength of the bubble-size
ordering rather than the 6.2 %. But this is a visual call and the page was not
built to depend on it.
