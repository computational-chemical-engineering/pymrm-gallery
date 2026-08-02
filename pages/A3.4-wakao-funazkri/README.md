# A3.4 — Wakao–Funazkri particle-to-fluid mass and heat transfer

`Sh = 2 + 1.1 Sc^(1/3) Re^0.6`, from Wakao, N. & Funazkri, T., *Chemical
Engineering Science* **33**(10) 1375–1384 (1978),
[doi:10.1016/0009-2509(78)85120-3](https://doi.org/10.1016/0009-2509(78)85120-3).

The paper measures nothing. It takes thirty-five years of published packed-bed
mass-transfer data, observes that the coefficients were extracted with the wrong
bed model, and extracts them again with the right one. The page rebuilds that
re-analysis.

> **Status: staged, ready.** The page ships with its figure results disclosed as
> pending re-confirmation. Figure 2's marker identifications were confirmed on
> 2026-08-02 — *"mostly symbols are identified well, but often the marker is off
> centre"* — and the re-centring that followed, the markers it added, and the
> whole of Figure 3 are batched for the next overlay review. Nothing waits on it.
> **Figure 3's axis calibration changed on 2026-08-02** — a column-dependent row
> origin, correcting a 5.6 px decade-tick skew — so the overlays are not the ones
> any earlier review saw. See `../../A3.4.yaml` and `../review/`.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 5 s |
| `build_page.py` | regenerates `index.ipynb` from source. Edit this, never the notebook |
| `data/wakao-funazkri-1978-fig2.csv` | Figure 2, gas phase at Sc = 0.6: 79 distinct glyphs in 90 rows. One row per marker position ever recorded, so the effect of the re-centring is a paired sample; group by `marker_id` |
| `data/wakao-funazkri-1978-fig3.csv` | Figure 3, liquid phase: (Sh−2)/Sc^(1/3) vs Re — **the data α and β were fitted on**. Six marker sets, one per erasure-band half-width; select `band_hw == 5` for the shipped 182, and use the rest as the systematic uncertainty |
| `data/wakao-funazkri-1978-fig3-line.csv` | the correlation Figure 3 draws inside itself, traced from its own ink. Not data — a calibration control whose answer is known to be 1.1 and 0.6 |
| `data/wakao-funazkri-1978-parameters.csv` | every constant read off 600 dpi page renders, plus both figures' axis calibrations, so the page can convert pixels itself |
| `meta.yaml` | page metadata |

Each CSV has a `.meta.yaml` provenance sidecar. **Read both figure sidecars
before using those files** — they record what each extraction is known to miss,
what the review has been asked, and (for Figure 3) which direction the known
recall loss biases a refit.

Both figures were extracted by fitting the marker **shape** and taking the centre
of the fitted shape. The method, its four traps, and the measured effect are
written up in `../review/README.md`, for promotion into `docs/handoff.md`.

## Regenerating

```bash
python build_page.py                    # rewrite index.ipynb
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb = nbformat.read("index.ipynb", as_version=4)
NotebookClient(nb, timeout=1800, kernel_name="python3",
               resources={"metadata": {"path": "."}}).execute()
nbformat.write(nb, "index.ipynb")
EOF
```

## Three things to know before reusing this

**The reusable object is the bed model, not the correlation.** `theta_pymrm`
solves eq. (4) — steady 1-D convection–diffusion–reaction with Danckwerts
boundary conditions, the `S3` structure — written purely in `(N, Pe_L)`, so it
applies unchanged to a first-order catalytic bed, an adsorber, or any closed
vessel with a linear sink. `theta_eq7` is its analytical solution and makes a
free regression test.

**`Sh_of_theta` is the piece worth lifting.** Every packed-bed transfer
coefficient in the literature carries an implicit axial dispersion assumption,
and that function converts between assumptions. Pre-1978 correlations and modern
ones are usually not the same quantity.

**Use eq. (12) with eq. (2), never on its own.** Conclusion 3 of the paper is
explicit: the correlation was extracted under a particular dispersion
coefficient and is only self-consistent when paired with it. Using it with
`D_ax = 0`, or with the inert-tracer coefficient of eq. (3), puts back the error
it was built to remove.

## Headline result

The whole paper reduces to one ratio. Eq. (2) and the `Pe = 2` assumption it
replaces share the same convective term, so

```
D_ax(eq. 2) / D_ax(Pe = 2) = 1 + 40 / (Sc·Re)
```

which doubles the dispersion coefficient below `Sc·Re = 40`. For a gas that is
`Re < 67` — most of the range the old correlations were fitted over. For a
liquid it is `Re < 0.04`, below anything the paper accepts. That single line
predicts both halves of the paper's first conclusion, and the numbers confirm
it: the liquid-phase correction never exceeds 1.33 % over `3 ≤ Re ≤ 10⁴`.

Pushing Petrovic and Thodos' own correlation through the bed model twice — once
at their `Pe = 2`, once with eq. (2) — returns eq. (12) to **1.8 % on average
and 4.5 % at worst** over `100 ≤ Re ≤ 900`, with nothing fitted and almost no
dependence on the assumed bed height. The paper never prints that comparison.

## And the constants, against the data they were fitted on

Figure 3 plots `(Sh−2)/Sc^(1/3)` directly, so the correlation reduces there to a
plain power law with no Schmidt number in it. **Hold β at the printed 0.6 and
the 182 digitised markers return `α = 1.1000`,** against the printed 1.1 — and
that number stays inside 1.098–1.118 across every erasure-band setting and both
axis calibrations, and inside 1.007–1.118 once the estimator is varied too — 
within 9 % of the printed value and never more than 2 % above it. The free α over
the same choices runs 1.05 to 1.47, up to a third high. **That contrast is the
result.** The β-fixed α is not a fixed point either: least squares in linear `y`
gives 1.047 where the log metric gives 1.100, and restricting to `Re ≥ 50` gives
1.02. The page prints both of those columns rather than claiming a stability it
does not have.

**The page does not claim the printed constants are wrong.** An earlier version
did, on a free log-space fit and its ordinary standard errors. That claim is
withdrawn and the page now shows why: fitting the same points unweighted in
linear `y` — an equally defensible loss function, and the paper gives no reason
to prefer one — moves α from 1.34 to 1.16; correcting a 5.6 px skew in the
ordinate calibration moved it further; the
erasure band alone spans 1.296 to 1.440; and a cluster bootstrap over the
Reynolds bins gives a 95 % interval of 1.04 to 1.53, which contains 1.1. What
survives is the β-fixed α above, and a qualitative excess in the markers below
`Re ≈ 10` — visible on the printed page, and half anticipated by the paper's own
rejection of liquid data below `Re = 3` for natural convection.

## The calibration control, and why it is on the page

Figure 3 prints its own correlation across the middle of the plot. That curve is
known in advance to be exactly `1.1 Re^0.6`, so tracing its ink and fitting a
power law to it measures the axis calibration against an object with no
scientific content. It came back `α = 1.126, β = 0.5949` — high in α and low in
β, the same signs as the disagreement being claimed about the markers. That is a
symptom, not a confirmation, and the cause turned out to be a 5.6 px skew
between the left- and right-hand decade ticks: the page render is rotated by
0.2°. With a column-dependent row origin the control returns `1.107 / 0.598`.
The trace is shipped as a dataset so that the control runs on the page.
