# A3.7 — van 't Riet's k_L a correlation for stirred gas-liquid vessels

`k_L A = 2.6e-2 (P/V)^0.4 v_s^0.5` for water and
`k_L A = 2.0e-3 (P/V)^0.7 v_s^0.2` for strong ionic solutions, from
van 't Riet, K., *Ind. Eng. Chem. Process Des. Dev.* **18**(3) 357–364 (1979),
[doi:10.1021/i260071a001](https://doi.org/10.1021/i260071a001).

**It is a review, and the page is built accordingly.** Roughly ninety papers are
cited and five figures compile other people's measurements. Two things in it are
van 't Riet's own — eqs. 8 and 9 — and they are what the page is about. Every
number the page reads carries an `origin` field saying whose result it is.

> **Status: staged, ready.** Runs in about 10 s. Nothing waits on a maintainer:
> the figures are deliberately **not** digitised, and the page claims nothing
> that would need them.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page |
| `build_page.py` | regenerates `index.ipynb`. **Edit this, never the notebook.** The physics lives in two module-level strings that are emitted into the notebook *and* `exec`'d by the builder, so every number in the prose is computed from the same source the notebook runs |
| `data/vant-riet-1979-correlations.csv` | eqs. 8 and 9: prefactors, exponents, validity windows, accuracy claim |
| `data/vant-riet-1979-printed-numbers.csv` | every other number the review prints as text, each tagged with **whose result it is** |
| `data/vant-riet-1979-figure-legends.csv` | the legend tables printed inside Figures 1–5 — **printed text, not extracted curve data** |
| `meta.yaml` | page metadata |

Each CSV has a provenance sidecar with a `columns:` block. Read
`vant-riet-1979-figure-legends.meta.yaml` before using that file: it says
exactly what the legends are and are not.

## Four things to know before reusing this

**Use the branch that matches the liquid.** The two correlations differ by up to
6.0× over the printed window and are not interchangeable at any accuracy the
review claims. The review's own boundary is "ca. 10 g of NaCl/L".

**Do not trust the ionic `v_s` exponent to two figures.** Setting it to zero and
re-fitting the prefactor stays inside the review's own 20–40 % accuracy over the
span of Figure 4's two drawn correlation lines (26.4 %, 0.0045–0.047 m/s) and
over the wider span its data occupy (35.2 %, 0.0023–0.047 m/s) alike — and the
review says in the sentence introducing eq. 9 that the 0.2 was imported from
Zlokarnik rather than fitted.

**The validity windows are not the same window.** `500 < P/V < 10 000` W/m³ for
both; up to 2600 L for eq. 8 and 2–4400 L for eq. 9. Eq. 8's own validity
statement prints no lower volume bound, so the CSV cell is empty rather than
guessed — the running text one sentence above eq. 8 does say the compiled
volumes "may range between 2 and 2600 L", and the sidecar records both facts.

**If you are going to *measure* `k_L A` rather than predict it, size the gas
depletion first — and use the right group.** `N_G = (1-ε) k H_t/(v_s H)`, not
the review's `τ_G · k_L A`. They differ by `(1-ε)/(ε H)`, which is 0.300 for
oxygen in water at ε = 0.1.

## Headline results

**Figure 5's topology, from the printed constants.** Pure-water curves cannot
cross pure-water curves and ionic cannot cross ionic, so all nine candidate
crossings are pure-vs-ionic. Exactly **one** falls inside the plotted window —
the 4 cm/s water curve meeting the 0.5 cm/s ionic curve at **P/V = 826.5 W/m³**,
which is what the printed figure shows. The check reads all six transcribed
constants at once, and no figure coordinate is measured: only the crossing count
and the identity of the pair, both legible without extracting a position. The
window's left edge, 400 W/m³, is Figure 5's printed axis label `4.10²` —
transcribed text, declared as such — and the page prints its sensitivity: the
count is 1 for any left edge in (328.0, 826.5] W/m³, and also 1 at the
correlations' printed validity bound of 500 W/m³ (at 300 it would be 2).

**Which constants the data exercise.** With the prefactor re-fitted at the window
centre in every case:

| term | cost of deleting it | inside the review's own 20–40 %? |
| --- | --- | --- |
| eq. 8 `(P/V)^0.4` | 82.1 % | no — load-bearing |
| eq. 8 `v_s^0.5` | 73.6 % | no — load-bearing |
| eq. 9 `(P/V)^0.7` | 185.3 % | no — load-bearing |
| **eq. 9 `v_s^0.2`** | **26.4 % over Figure 4's drawn lines; 35.2 % over its data span** | **yes — not resolved** |
| the branching itself | 83.3 % under-prediction | no — the single most load-bearing feature |

**An internal inconsistency in the review.** It expects ions to raise `k_L A` by
2–10×, argued from bubble mechanism and Zlokarnik and **not** from eqs. 8 and 9.
The ratio of the two fitted correlations spans **1.30× to 5.98×** and lands
inside 2–10× on 83.6 % of the printed box — falling below the 2× floor in the
high-`v_s`, low-`P/V` corner. Reported unresolved; what would settle it is in
Figures 3 and 4, which are not digitised.

**What the well-mixed-gas assumption costs.** Eq. 2 — the method most of the data
come from — assumes a uniform gas composition. The review warns about this four
times and never computes anything. At the centre of the printed `P/V` window, at
the lowest velocity Figure 5 draws, in a 2600 L vessel, the page's own transient
gassing-out experiment shows eq. 2 under-reporting `k_L A` by **19.1 %** on the
pure-water branch (the quasi-steady closed form gives 16.0 % — optimistic,
because it assumes `τ_G ≪ 1/k_L A`, exactly the condition being violated) and by
**45.4 %** on the ionic one at the reported fit window. On the ionic branch the
response is genuinely non-exponential (`τ_G` = 29.8 s against `1/(kφ)` = 12.0 s),
so what eq. 2 returns is window-dependent — 34.1–56.5 % across the fit windows
the page sweeps, with the quasi-steady 45.6 % inside that range — and the
near-coincidence of the (0.2, 0.8)-window value with the quasi-steady φ is a
property of the window, not corroboration; the page says so where the two
numbers meet. Either way the bias is the same size as the entire accuracy the
review claims. Over the vessel volumes its own legends print, the quasi-steady
ionic figure runs 6.5 % at 2.5 L to 52.2 % at 5100 L, which is its sentence
about larger vessel diameters with numbers on it. And across the review's own
printed aspect-ratio window the tall vessel is the worst case — 59.3 % ionic
bias at T/H = 0.5 (H_t/T = 2) against 38.0 % at T/H = 1.5 — which is the
review's `H_t/T ≫ 1` warning with numbers on it.

**This is a self-consistency statement, not a correction.** If the underlying
measurements were biased low by φ, eq. 9 as fitted already contains that bias.
What is established is that eq. 9 taken at face value predicts conditions under
which the method that produced most of its data does not measure what it is
assumed to measure.

## What is deliberately not here

**No figure is digitised.** Figures 1–5 hold the only measurements in the review;
extracting them needs a maintainer overlay review and none is available. So the
page is a **reproduction, not a validation**, the data tier is 6, and no fit
quality for eqs. 8 and 9 is computed — because the data they were fitted to are
in exactly those figures.

**Three inputs are not printed in the review** and all three are marked: the
Henry constant `H` (the load-bearing one, swept 20–40), the gas holdup `ε` (which
cancels except through `1-ε`, swept 0.05–0.20) and the aspect ratio, taken at the
centre of the review's own printed `0.5 < T/H < 1.5` with both extremes swept at
constant vessel volume — the sweep prints the implied volume alongside, because
an earlier version computed the tank diameter where the height belongs, swept
0.325–8.775 m³ vessels instead of 2.6, and inverted the trend. The corrected
sweep supports the review's warning instead of contradicting it, and the buggy
formula is kept as a break-table row.

## Two printed defects, recorded and not repaired

1. The sentence under eq. 7 says *"where N = stirrer diameter"*; the
   Nomenclature gives `N` = stirrer speed and `D` = stirrer diameter.
2. **Four of the five figure legends cite the drawn correlation by an equation
   number one lower than the text does — and the fifth cites none at all.**
   Figures 1–3 print `CORRELATION EQUATION (7)` for eq. 8's line, and Figure 4
   prints `CORRELATION EQUATION (8)` twice for eq. 9's, on the very figure
   eq. 9 is derived from; Figure 5's legend names the liquids and cites no
   equation number. Eq. 7 carries no numerical constants and cannot be drawn at
   all. The offset is a column of the legend dataset — empty for Figure 5 — so
   it can be checked from the data.

## A pymrm note

**pymrm's Neumann outflow boundary extrapolates to the face rather than taking
the upwind cell.** A hand-written outlet flux `v_s·C_N` misses the operator's own
by 8.9e-3 at n = 8 and 7.8e-5 at n = 800, so a mass balance written that way
fails to close by ~1e-4 and looks like a physics error. The page's balance uses
the operator's face fluxes and closes at 1.4e-10. Worth promoting into
`docs/handoff.md`.

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
