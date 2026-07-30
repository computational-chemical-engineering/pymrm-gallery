# A3.4 — Wakao–Funazkri particle-to-fluid mass and heat transfer

`Sh = 2 + 1.1 Sc^(1/3) Re^0.6`, from Wakao, Kaguei & Funazkri, *Chemical
Engineering Science* **33**(10) 1375–1384 (1978),
[doi:10.1016/0009-2509(78)85120-3](https://doi.org/10.1016/0009-2509(78)85120-3).

The paper measures nothing. It takes thirty-five years of published packed-bed
mass-transfer data, observes that the coefficients were extracted with the wrong
bed model, and extracts them again with the right one. The page rebuilds that
re-analysis.

> **Status: staged, not published.** The Figure 2 digitisation is waiting on a
> maintainer review. See `../../A3.4.yaml` and `../review/`.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 25 s |
| `build_page.py` | regenerates `index.ipynb` from source. Edit this, never the notebook |
| `data/wakao-funazkri-1978-fig2.csv` | 81 marker positions digitised from Figure 2 (gas phase, Sc = 0.6) |
| `data/wakao-funazkri-1978-parameters.csv` | every constant, read off 600 dpi page renders |
| `meta.yaml` | page metadata |

Each CSV has a `.meta.yaml` provenance sidecar. **Read
`wakao-funazkri-1978-fig2.meta.yaml` before using that file** — it records what
the extraction is known to miss and what the review has been asked.

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
