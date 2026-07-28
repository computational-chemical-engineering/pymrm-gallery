# F1.4 — Krishna–Ellenberger two-bubble-class gas holdup

Large-bubble gas holdup in a churn-turbulent bubble column, from Krishna &
Ellenberger, *AIChE Journal* **42**(9) 2627–2634 (1996),
[doi:10.1002/aic.690420923](https://doi.org/10.1002/aic.690420923).

The page tests Eq. 19 — a correlation containing no fluid property at all —
against the 63 markers of the paper's Figure 11, and against the Wilkinson
correlation it replaced.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 3 s |
| `build_page.py` | regenerates `index.ipynb` from source. Edit this, not the notebook |
| `data/krishna-ellenberger-1996-fig11.csv` | 63 marker positions digitised from Figure 11 |
| `data/krishna-ellenberger-1996-parameters.csv` | Table 2 properties, gas densities, Reilly's *B*, the Table 3 deviations |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it. Read
`krishna-ellenberger-1996-fig11.meta.yaml` before using the `gas` column — it
records a maintainer review of the extraction and the reason most rows are
deliberately left `unassigned`.

## Regenerating

```bash
python build_page.py                    # rewrite index.ipynb
python ../../scripts/check_metadata.py  # then execute it via scripts/run_pages.py
```

## Two things to know before reusing this

**This is a correlation, not a discretised model.** Eqs. 8, 19 and 20 are
algebra; no pymrm operator appears on the page, and the page says so. The three
functions `eq19`, `wilkinson` and `reilly_transition` take SI arguments and are
standalone — lift them into a reactor model. `F2.3` is the page that does.

**The gas labels are incomplete on purpose.** The four series on Figure 11
differ only by marker shape, drawn as overlapping open outlines. Automatic shape
classification failed in the dense band and a review confirmed it, so only the
SF₆ group is labelled. This costs nothing for the main test — Eq. 19 has no
gas-density term, so all 63 points can be used without labels — and the SF₆
group alone carries the independence test, being the density extreme — see the
headline result below for how that test has to be framed.

## Headline result

Eq. 19 sits through the points with a mean absolute deviation of 13.8 % and a
bias of +2.8 %, comparable to the δ = 0.16 the authors report over their full
1,735-run set. Wilkinson deviates 4.6× as far and is biased 64 % high.

For gas independence, note that comparing group biases does not work on this
figure — the SF₆ points and the rest occupy disjoint velocity windows, so a
group difference is confounded with velocity. The page tests it by extrapolation
instead: a free power law fitted on the 51 helium/air/argon points predicts the
12 SF₆ points, 37× denser and unseen by the fit, with a bias of −7.7 % and the
same scatter as the correlation itself. Wilkinson requires 20–27 %.
