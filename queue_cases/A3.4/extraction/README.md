# The A3.4 figure extraction, as it was run

These five modules regenerate both datasets in `../page/data/` and both review
overlays in `../review/` from the source PDF. They are kept here rather than in
the session scratchpad because the technique is meant for promotion into
`scripts/` — see `../review/README.md` for the write-up, and
`docs/handoff.md` for where it belongs.

| module | what it is |
| --- | --- |
| `glyphfit.py` | the method: parametric glyph templates, the body-minus-margin score, sub-pixel fitting, matching pursuit, the straight-ink opening. **Figure-independent** — this is the part worth promoting |
| `setup2.py` | Figure 2: crop, calibration, what to erase and how |
| `setup3.py` | Figure 3: the same, plus a **column-dependent** least-squares calibration from the printed decade ticks, and `trace_line`, which traces the drawn correlation's own ink as a calibration control |
| `overlay.py` | numbered overlays and 2× contact sheets for review |
| `emit.py` | the driver: writes all three CSVs and all four overlays, sweeps the Figure 3 erasure band, and applies the two hand-audit decisions (`FIG2_AUDIT_REJECT`, `FIG2_FIT_OK`) |

## Running it

Render the two figure crops first — the PDF is not in the repository and never
will be:

```bash
pdftoppm -r 600 -f 8 -l 9 -png ~/papers/pymrm-gallery/Wakao1978-particle-to-fluid-transfer-CES33-1375.pdf p
python - <<'EOF'
from PIL import Image
Image.open("p-08.png").convert("L").crop((2296, 4200, 4592, 5800)).save("renders/fig2-full.png")
Image.open("p-09.png").convert("L").crop(( 400, 2450, 2250, 3760)).save("renders/fig3-full.png")
EOF

A34_RENDERS=renders python emit.py
```

`emit.py` reads the existing `wakao-funazkri-1978-fig2.csv` for the 2026-07-30
crosshairs and re-derives everything else, so it is idempotent: run it twice and
you get the same file.

## The two hand decisions, and where they live

Nothing else in this pipeline is hand-tuned per marker, but two things are, and
both are constants at the top of `emit.py`:

- `FIG2_AUDIT_REJECT = [7, 9]` — two of the eleven markers the pursuit proposed
  on Figure 2 were rejected on sight: one sits on the "Single spheres" leader
  line, one where that leader crosses the eq. (11) line. They are indices into
  the pursuit's output order, so **they are only valid for these settings**. If
  any detection parameter changes, re-audit before trusting them.
- `FIG2_FIT_OK = 0.50` — below this fit score a marker keeps its 2026-07-30
  crosshair instead of the fitted centre. **Sixteen** of eighty-one did, and
  they are flagged `recentred = False` in the CSV. Their displacement is zero by
  construction, so any statistic about how far the re-centring moved a marker
  must exclude them: on the 65 that moved the median displacement is 3.29 px,
  against 2.28 px if all 81 rows are averaged together.

## The two things this pipeline gets wrong if you run it naively

Both were found by adversarial review on 2026-08-02 and both are fixed here.

**Figure 3's ordinate is not a single scale.** The left- and right-hand decade
ticks disagree by 5.6 px because the render is rotated by 0.203°, and averaging
them — which the first pass did — is right in the middle of the plot and wrong
at both ends. Since the error is a *tilt* it feeds straight into β. `calibrate()`
now fits each side separately and `y_of_row(row, col)` takes the column;
`y_of_row_flat` is kept only so the page can show what changed. The symptom was
already in the sidecar and had been read as confirmation: fitting the drawn
correlation's own ink returned α 2.4 % high and β 0.9 % low, in the same
direction as the disagreement the page was claiming about the markers. A control
that fails in the direction of your result is a symptom, not a confirmation.

**The erasure band is not a detail; it is the dominant systematic.** `emit.py`
sweeps `BAND_HW = (0, 3, 5, 7, 10, 14)` and ships every marker set, because
across that sweep the free α runs 1.296 to 1.440 — a half-width of 0.072 against
an ordinary standard error of 0.052 and a cluster-bootstrap one of 0.132, i.e.
larger than the error bar the page was quoting. Shipping one
setting and describing the sensitivity in prose is how a "the fit is insensitive"
claim gets made about the wrong knob.
