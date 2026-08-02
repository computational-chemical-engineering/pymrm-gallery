# G1.8 — Trickle-bed reactor with partial catalyst wetting

Herskowitz, M. and Smith, J. M. (1983). *Trickle-bed reactors: a review.*
AIChE Journal **29**(1) 1–18. doi:10.1002/aic.690290102

Reproduces **Figure 6** — the gain in overall effectiveness factor when part of
a catalyst pellet is dry and the limiting reactant is volatile — from the
paper's own Table 1 (sphere row, page 4) and Eqs. 19–21 (page 8), with
$\alpha_{gs}\to\infty$ and $C^*_L = 1$ as stated.

## What this page is

**Tier 6.** Figure 6 contains no measurements; its four lines are the authors'
own computed model output. Reproducing them tests a transcription and a chain of
algebra. Nothing on this page is experimental validation.

## The headline

The figure **does not reproduce as printed**. At $\phi = 10$ the paper's own
chain lands at +1.0 %, +69.2 %, +314.3 % and +151.8 % of where the four drawn
curves sit — right for one curve out of four.

Solving for the $L_m$ that puts each drawn curve where it is gives **0.508,
2.046, 7.061 and 10.18**. Three of those are the legend's own printed 0.50, 2.0
and 7.0, reused in a different order. Reading the four curves as
$L_m$ = 0.50, 2.0, 7.0 and ≈10 brings all four to +1.0 %, +1.9 %, +1.6 % and
+6.4 % at $\phi = 10$, with **nothing fitted** for the first three.

**The legend block is misaligned by one row against its curves**, from the
second row down: its $L_m = 1.0$ entry belongs to no drawn curve, and the fourth
drawn curve has no entry. Endorsed by the gallery maintainer on 2026-08-02 from
a private decision artifact. The page shows both readings side by side.

Two further readings of the same figure support it. **Neither is an independent
witness** — both are the same four measured curve positions, re-expressed:

1. The three inter-curve **gaps**, in decades (at $\phi = 10$; they drift by up
   to 0.031 decades across the window): 0.415/0.609/0.394 on the figure, against
   0.412/0.611/0.374 for the shifted reading and 0.191/0.221/0.611 as printed.
   This is the position measurement with one degree of freedom — a common
   multiplicative offset in $\chi$ — removed, so what it rules out is
   "printed reading right, $\chi$ calibration wrong". It needs the ordinate's
   decades-per-pixel and is not calibration-free.
2. Because $\chi = g(L_m)\,h(\phi)$ **separates**, those same gaps force
   $d\log g/d\log L_m$ with no model at all. As printed they demand
   −1.380/−2.024/−0.725 — non-monotone, so unreachable by any product of powers
   of $L_m$. Shifted they demand −0.690/−1.120/−2.546 against the model's
   −0.684/−1.122/−2.414, within 5 %. No reparameterisation rescues the printed
   reading; only a reassignment does.

## The residual the shift does not fix

The four drawn lines measure log–log slopes **1.098–1.131** (mean 1.111); the
model's best straight-line fit over the same window is **1.077**. $L_m$ cannot
produce that — the chain factorises into an $L_m$-dependent prefactor times a
$\phi$-dependent shape, so every $L_m$ gives the identical slope.

The **sign** is robust; the **size** is not. The dominant uncertainty is the
abscissa calibration, not the fit rms: the four labelled $\phi$ verticals do not
sit on a common log ruler ($\phi = 10$ is ~4 px off), and an independent 600 dpi
re-digitisation measures a mean slope of 1.1028 against this CSV's 1.1115 —
exactly the ratio of the two calibrations. Across defensible calibrations the gap
spans **+0.007 to +0.046** (+3 % to +19 % end to end). **Reported as unexplained,
with its size as a range.**

## Two printed defects reported

1. **Page 8 cites "Table 2"** for the approximate spherical solution. Table 2
   holds the pressure-drop constants $\beta$ and $\gamma$; the sphere row is in
   **Table 1**, page 4. Checked on 600 dpi renders of both pages and against the
   notation list.
2. **Figure 6's legend is misaligned by one row** against the curves, as above.

## Validation

Ranked before any code was written, and every check followed by a deliberate
break with the number re-measured:

| check | baseline | breaks by |
|---|---|---|
| Eq. 20 vs the $f_e$ printed inside the figure | 4.7e-03 | 25×–436× |
| pymrm finite volume vs Table 1's sphere row | 7.7e-06 | 1.2e3×–2.6e5× |
| the $L_m$ reconstruction | 2.3 % | 0.75×–1302× (see caveat) |
| the collapse identity | 1.8e-15 | **does not move at all** — kept and labelled powerless |

The reconstruction check has a real blind spot and the page states it. Its
tolerance is 3.8–11.4 % in $L_m$, and **four** injected defects stay under even
the tightest of those bars: Eq. 21's `6.91 → 7.91` (1.72 %), `6.91 → 6.71`
(2.71 %) and its exponent `0.6 → 0.5` / `0.6 → 0.7` (3.33 % / 3.43 %). The first
of these makes the metric look *better* than the undamaged 2.30 %. The reason is
structural: at $L_m \approx 1$ the $1.05\,L_m^{0.3}$ term carries **87 %** of
$1/\alpha_{gLs}$, so corrupting the $6.91\,L_m^{0.6}$ term barely moves $\chi$.
The check catches Eq. 20 and Eq. 21's first term; against Eq. 21's second term it
is essentially blind. It is also weak against an ordinate-scale error (±5 %
leaves it at 6.3 % / 10.6 %), which is why the gap test is kept.

## Files

```
index.ipynb        the page (generated — edit build_page.py, not this)
build_page.py      the generator
meta.yaml          page metadata
data/herskowitz-smith-1983-fig6.csv        four fitted lines, not markers
data/herskowitz-smith-1983-fig6.meta.yaml  provenance + both review rounds
```

Runtime ≈ 9 s.

## Regenerate

```bash
source ~/Code/pymrm_suite/.venv/bin/activate
cd <this directory>
python build_page.py
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb = nbformat.read("index.ipynb", as_version=4)
NotebookClient(nb, timeout=1800, kernel_name="python3",
               resources={"metadata": {"path": "."}}).execute()
nbformat.write(nb, "index.ipynb")
EOF
```

## Note for the integrator

`G1.8` **already has a record in `models.yaml`** with `status: planned`
(slug `trickle-bed-partial-wetting`). Upgrade it **in place** using
`../models_entry.yaml`; appending a second block makes `check_metadata.py`
report a duplicate id.
