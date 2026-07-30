# Handoff — state as of 2026-07-31

Start here if you are picking this up fresh. Read this, then
[`AGENTS.md`](../AGENTS.md), then [`pdf-findings.md`](pdf-findings.md).

## Built and live

**https://computational-chemical-engineering.github.io/pymrm-gallery/**

13 pages, 14 published catalog entries, both CI workflows green.

| Page | What it shows | Validation |
|---|---|---|
| `A4.9` Duncan–Toor ternary diffusion | Osmotic → uphill → diffusion barrier | **experimental** — 0.59 mole % vs the paper's own 0.45%, 28 digitised points |
| `C2.1` Xu–Froment steam reforming | The most-used SMR kinetics, against the runs they were fitted to | **experimental** — 0.0017 in conversion (2.7%) over 61 digitised points, nothing fitted |
| `B1.1`+`B1.5` Thiele + Weisz–Hicks | η(φ), and 3 steady states at one φ | 2.2e-4 vs exact; both methods agree on all three branches |
| `D2.2` Van Welsenaere–Froment runaway | Two criteria for the runaway boundary, swept over the operating plane | 0.054% over all 30 numbers in their Section 6; two independent methods agree to 0.18% |
| `F3.1` Hatta regimes | Enhancement factor, 3 regimes | 6.3e-3 vs exact; VKH good to 2.1%, DeCoursey to 8.7% |
| `A2.3` Taylor–Aris dispersion | Homogenisation closure, and when the lumped coefficient becomes defensible | 1.0e-4 vs Eq. 25 at n_r=200, O(h²); Taylor's own capillary run to 0.04% |
| `J1.5` LDF breakthrough | What the linear-driving-force constant actually is | 6.6e-5 vs the exact series at n_r=400 |
| `F1.4` Krishna–Ellenberger holdup | A correlation with no fluid property in it, against the figure it was fitted to | **experimental** — 13.8% mean deviation, +2.8% bias over 63 digitised points |
| `H1.7` Wijmans–Baker solution–diffusion | Two constants predict the third — and the figure cannot test the prediction | **experimental** — A and B fitted to 8 points, rejection predicted for the other 4 |
| `B3.1` Yagi–Kunii shrinking core | The one equation all three textbook regimes come from, plus a map of where they are safe | 6.9e-16 against an independent derivation; all three limits to 2.4e-8 |
| `F2.3` Maretto–Krishna FT slurry column | Plug-flow large bubbles over a well-mixed slurry, and two printed constants that stop it working | **experimental** holdup — 5–6% over 79 points; conversions 93.1/63.8% vs the paper's 96/63% |

Status counts live in `models.yaml`: 14 published, 15 planned, 2 deferred
(`H1.12`, `B1.12` — unpublished manuscripts, see the published-work-only policy
in [`blueprint.md §9`](blueprint.md#published-work-only-policy)).

**Three of eight pages are validated against experiment** (`A4.9`, `C2.1`,
`F1.4`), and `A2.3` against Taylor's own worked capillary run. `D2.2` is tier 6
by necessity — its source paper contains no measurements at all. `A3.4` is the
next chance to move the ratio, and it needs figure digitisation.

## Papers available

Nine PDFs at **`~/papers/pymrm-gallery/`** (priorities 1 and 2 complete).
Inventory with per-file text-layer quality is in
[`pdf-findings.md`](pdf-findings.md).

## Elsevier full text — use this instead of OCR

An institutional API key is stored at **`~/.config/elsevier/apikey`** (mode 600,
outside the repo; `.gitignore` blocks credential filenames — never commit it).

The **DOI endpoint 404s** on older Chemical Engineering Science articles because
their DOIs contain parentheses. **Use the PII endpoint**, and note the PII is
recoverable straight from the ScienceDirect filename:
`1-s2.0-`**`0009250970850734`**`-main.pdf`.

```bash
K=$(cat ~/.config/elsevier/apikey)
curl -sS -H "X-ELS-APIKey: $K" -H "Accept: text/plain" \
  "https://api.elsevier.com/content/article/pii/0009250970850734"
```

Verified working on Van Welsenaere & Froment (47 kB) and Wakao & Funazkri
(45 kB) — both of which had the *worst* local OCR of the set.

**Correction, 2026-07-26 — do not trust the API text for numbers.** An earlier
version of this file claimed the API "returns properly encoded text, so
exponents and subscripts survive". That is wrong for pre-1980 scans. What comes
back is the *publisher's* OCR of the same scan, and on Van Welsenaere & Froment
it drops decimal points wholesale:

| API text | Truth (600 dpi page image) |
|---|---|
| `R = 00125 m` | *R* = 0.0125 m |
| `M = 2948 kg/kmole` | *M* = 29.48 kg/kmole |
| `(p°),,, = 001353 atm` | (*p*⁰)ₗ = 0.01353 atm |
| `b = 19837` | *b* = 19.837 |
| `t,„ = 21 -818` | *t*ₘ = 21.818 |

The 1970 typesetting uses a mid-dot decimal separator, which the OCR discards.
So the API is excellent for **prose** — section structure, what each figure
shows, which assumptions are stated — and useless for **parameters**. Read
those from a 600 dpi page render, the same as for a Wiley scan. Treat any
integer with an implausible magnitude as a lost decimal point and go to the
image; never "fix" it by inference.

**Other limits:** Elsevier only. Xu & Froment and Krishna & Ellenberger are
Wiley and have no API route at all. And full text never yields figure data —
plots still need digitising. Be polite with request volume; this is not a
bulk-download tool.

## Extraction cost, revised

| Page | Paper | Extraction | Value |
|---|---|---|---|
| `D2.2` | Van Welsenaere & Froment | API for prose, **page image for every number** | Runaway criteria; sweep-based figures are striking |
| `A3.4` | Wakao & Funazkri | API for prose, page image for numbers | Sh–Re dataset, but likely a scatter *figure* → digitise |
| `E2.1` | Kunii & Levenspiel | good OCR (12.4k chars/page) | The canonical fluidised-bed model |
| `I1.2` | Oh & Cavendish | good OCR (9.0k) | Converter light-off, `S4`+`S7` |
| ~~`F1.4`~~ | ~~Krishna & Ellenberger~~ | ~~tables clean, holdup in figures~~ | **done 2026-07-28** |
| ~~`C2.1`~~ | ~~Xu & Froment~~ | ~~hard~~ | **done 2026-07-26** |
| ~~`D2.2`~~ | ~~Van Welsenaere & Froment~~ | ~~API + page image~~ | **done 2026-07-27** |

## Recommended next moves

1. **`J3.4`** — the largest build left, and now **specified rather than just
   staged**. The paper is transcribed into the sidecar's `model_inventory:`
   block, including one finding that changes the build: the paper's own Eq. 26
   gives S_c = 1.0e-4, so **solid-phase diffusion is negligible** and Appendix
   B's superposition machinery comes off the critical path. Three items remain to
   extract (Appendix A's conductivity polynomial and salt diffusivity, Eq. 16's
   open-circuit potential, Eq. 17's Butler–Volmer). It can only ever be a
   reference-solution page — Figure 2 contains no measurements. **`G1.8` is
   BLOCKED** on a question for the maintainer; see below. (`F2.3` — **done
   2026-07-30**.)
   `F2.3` can lift Eq. 19 and Reilly's Eq. 8 from the `F1.4` page as validated SI
   functions. `J3.4` is a full P2D battery model, the largest build in the queue,
   and can only ever be a reference-solution page. (`B3.1` — **done 2026-07-29**.)
2. **`E2.1` Kunii & Levenspiel** — best text layer of the set, and the canonical
   fluidised-bed model.
3. **`A3.4` Wakao & Funazkri** — the Sh–Re dataset. Elsevier, so the API gives
   the prose, but read every number off the page (see the correction above).
   Note `F2.3` consumes `F1.4`: Eq. 19 and Reilly's Eq. 8 are already
   implemented and validated on the `F1.4` page as standalone SI functions, so
   lift them rather than rewriting.

### When a printed constant is wrong, prove it from the paper's own results

`F2.3` needed two corrections before it would run, and the method for
establishing them generalises.

Eq. 2's rate prefactor is printed as 8.8533e3 mol/(s kg_cat bar²), which gives an
intrinsic rate 10⁶ larger than any cobalt catalyst. Eq. 1's rate is labelled
`R_CO+H2` but behaves as a CO rate. Neither could be fixed by fitting — that would
have made the whole comparison circular.

**What made the diagnosis safe was establishing chemical control first.** The
paper reports that a 10-fold rise or 3-fold fall in kLa is negligible; reproducing
that *before* touching the kinetics proves the mass-transfer correlations are not
free to absorb a rate error. Only then does the conversion comparison isolate the
kinetics, so each correction becomes a discrete choice between stated
alternatives — 10³ vs 10⁻³, syngas vs CO — with the paper's own reported
conversions selecting between them. The page prints the alternatives and what each
gives.

*Order matters: pin down what is NOT free before claiming a constant is wrong.*

Two pymrm traps recorded on that page, both silent failures:

- A convection outlet left as `None` makes the matrix **singular**, and a
  rank-deficient solve still returns a plausible-looking profile.
- With varying velocity, discretise `d(Uc)/dz` as the divergence of the flux.
  `U dc/dz` loses the gas contraction — 65 % of the volumetric flow here.

### `G1.8` is blocked — a model that closes for one curve out of four

Everything needed to reproduce Herskowitz & Smith's Figure 6 was found and each
piece verified on a 600 dpi render: **Table 1**'s sphere row (page 8 cites
"Table 2", which actually holds the pressure-drop constants — an error in the
paper), Eq. 19 for χ, Eq. 20 for *f*_e, Eq. 21 for α_gLS, with α_gs → ∞ and
*C*\*_L = 1 as stated.

Three independent checks pass. Eq. 20 reproduces the four *f*_e values printed
inside the figure exactly. Table 1's sphere expression reduces to the classical
sphere effectiveness factor when *f*_e = 1, which fixes φ on the *V*/*S* length.
The chain collapses to χ = (1 − *f*_e)·η_s(φ)·φ²/α_gLS, whose log-log slope is
1.03 against 1.10–1.13 measured.

**And then it reproduces only the lowest curve.** L_m = 0.50 comes out at 3.767
against 3.752 digitised at φ = 10 — 0.4 %, nothing fitted. L_m = 1.0, 2.0 and 7.0
are out by 1.69×, 4.03× and 2.34×. The line positions are not in doubt: an
overlay of the fits sits on the printed curves and the legend rows align with the
curve endpoints.

Nothing reconciles it cleanly. Matching each curve needs α_gs = 2462, 16.2, 4.7,
10.3 — no pattern, all contradicting "very large". Fitting *f*_e instead needs
0.721, 0.864, 0.958, 0.974 against the printed 0.72, 0.77, 0.83, 0.94, a
correction that grows with L_m.

So either the upper three curves used something unstated, or the figure is wrong.
**Do not publish the page claiming agreement.** Asserting a 1983 AIChE review is
in error is not a call to make from arithmetic alone; it is a question for the
maintainer. Full detail is in the staged sidecar's `model_closure_attempt:` block.

### Working the catalogue with parallel agents

`queue_cases/` holds one YAML per catalogued case; `scripts/case_queue.py` and
`scripts/dashboards.py` drive it, and `scripts/find_papers.py` does acquisition.
Agents build into `queue_cases/<ID>/page/` and touch nothing shared, so several
run at once; an integrator merges finished pages into `pages/`, adds the
`models.yaml` entry, and runs the gates.

**Three things that went wrong in the first batch, all worth avoiding again.**

**Never `git add -A` after an agent runs.** Review overlays are drawn on the
source page image, so they *are* the copyrighted figure — three were committed to
this public repo before being caught, contradicting the redistribution basis every
sidecar states. `queue_cases/*/review/*.png` is now git-ignored, but the gate is
human: look at what an agent produced before staging it.

**A DOI resolved from a terse citation is usually wrong.** CrossRef returns
something confident for any query: "Carman (1937)" matched a 2025 paper citing
Kozeny–Carman, and "Danckwerts" matched a 1968 re-derivation. Title-word overlap
is not enough — require the publication year to agree, and mark anything
auto-resolved as unverified so a wrong DOI never sends the maintainer after the
wrong paper.

**Filenames carry no metadata.** Half the PDFs on disk are named by publisher PII
(`i260028a001.pdf`) or an export id, so automatic matching missed them and one
case was reported as needing a paper the maintainer had already supplied.
`docs/papers-on-disk.yaml` maps catalogue ID to filename by hand and is consulted
first; add a line whenever a PDF arrives.

**Open access will not carry section A.** Of the first cases checked, most 1937–75
classics have neither DOI nor open copy — they predate DOIs. Expect the papers
list to be dominated by old classics and to work well only for recent sections.

### Verify an equation read off a page image before building on it

`B3.1` is the template. Its two governing equations were read from 600 dpi
renders because the scan's text layer mangles them (`theta_B` as `0B`, Eq. 6's
exponents dropped). A page-image reading is a transcription and needs checking
like any other. Two checks settled it:

- **Endpoint identities that only hold for the right coefficients.** Eq. 6 must
  give exactly 0 at *r*/*R* = 1 and exactly 1 at *r*/*R* = 0 for *any* parameter
  values, which requires its numerator to collapse to its denominator term by
  term. A mis-read coefficient breaks this immediately.
- **An independent derivation.** Integrating the moving boundary from Eq. 5's
  three resistances, without looking at Eq. 6, reproduced it to 6.9e-16 over six
  decades — and recovered the factor 3 on the film term and the 12 in
  *k*_d1 = 12D/D_p on the way.

Look for both before writing the notebook. This is the same instinct as the
`C2.1` Table 5/6 round trip, generalised: *a transcription you can only read once
should be checked against something you can derive.*

**And refuse to invent what the paper omits.** `B3.1` cannot produce absolute
burnout times, because Parker & Hottel's correlation is printed for the specific
combustion rate `K_c` and the unit conversion to `k_c1` is not given. The page
works entirely in the dimensionless groups instead, which costs nothing, and says
so. Do not reconstruct a missing unit conversion by inference.

### The batched figure review, 2026-07-29 — read this before digitising anything

Five figures were digitised, put to the maintainer as **one** review artifact
(source figure ⇄ overlay toggle, closed questions plus a free-text box per
figure), and all five came back answered in a single pass. That is the workflow
to repeat: batching costs nothing extra and turns five review round-trips into
one. The artifact is private and no page image enters the repo.

What came back, and what it changed:

- **`H1.7` Wijmans & Baker Fig. 5** — all correct. Page now built and published.
- **`F2.3` Maretto & Krishna Fig. 2** — *"On the eps_s=0.35 line you detect a few
  circles and diamonds. These are just squares almost on top of other squares."*
  Fixed, see the lesson below. Also: report the unresolved cluster flagged rather
  than drop it, and the two circles floating above the top curve are real data.
- **`J3.4` Doyle–Fuller–Newman Fig. 2** — confirmed **no experimental data at
  all**; reference-solution page only. One print speck removed.
- **`G1.8`** — **switch to Figure 6, do not digitise Figure 2.** Fig. 2 needs the
  identity of ~45 overlapping markers and the correlation under test depends on
  particle diameter, so the F1.4 shortcut does not apply. Fig. 6 is the model.
- **`B3.1`** — agreed it is analytic; **no figure needed, no review gate.**

Staged extractions with their review verdicts are in `docs/staged-data/`.

**Series identity can come from position instead of shape.** This is the F1.4
lesson in a second form and it is the reusable one. On `F2.3` the template
matcher picked the wrong *shape* wherever markers overlapped — two overlapping
squares read as a circle. But the three series each follow their own curve and
never cross, so identity was reassigned by fitting `eps = a + b·log U + c·log²U`
per series and moving every marker to the nearest curve, iterated to a fixed
point. That moved exactly the 10 markers the reviewer had described, without
touching the detector. **Ask what carries the series identity in the figure —
shape, position, or a curve — before trying to improve shape recognition.**

**Check whether the figure can resolve what you are testing.** On `H1.7` the
rejection panel appears to validate the model; closing the model shows it
predicts a 0.31 percentage-point rise, which is 4 px on a figure whose curve is
6 px thick. Agreement there is not evidence. Compute the predicted effect in
*pixels* before quoting an agreement.

**Reusable extraction code** lives in the session scratchpad, not the repo:
`markers.py` (shape-template matching pursuit, open and filled), `curves.py`
(column tracing with an order-preserving DP so adjacent curves cannot swap), and
`overlay.py`. Worth promoting into `scripts/` when a fourth figure needs them.

### What `F1.4` settled, and what it did not

Figure 11 is digitised (63 markers) and the page is live. Three things are worth
carrying forward.

**The figure was salvaged by dropping labels, not by fixing the classifier.**
The four series differ only by marker shape, and shape recognition failed in the
dense band — the maintainer review confirmed it. Rather than curate a subset,
every row except the SF₆ group is `unassigned`. This cost nothing, because
**Eq. 19 contains no gas-density term**: the correlation can be tested against
positions alone, using all 63 points. Labels are needed only for the
gas-independence claim, and the SF₆ group alone carries that, being the density
extreme. *Generalise this before spending hours on a classifier: ask which
columns the model you are testing actually reads.*

**A group comparison can be confounded even when the groups look clean.** The
first version of the gas-independence test compared the SF₆ points' bias against
the rest — until a check showed the two groups occupy *disjoint* velocity
windows on this figure, every SF₆ point below 0.044 m/s and every other point
above 0.051. The difference measured velocity, not density. Replaced by an
extrapolation test: fit a free power law on the lighter gases only, predict the
SF₆ points it never saw. Before reporting a between-group difference, check the
groups overlap in every other variable.

**Deviation direction is not cosmetic.** The first draft computed the two
correlations' deviations in opposite senses — measured/model for one,
model/measured for the other. At 14 % scatter the reciprocal differs materially:
the mean moved 13.3 → 13.8 %, the bias −0.2 → +2.8 %, and the headline
gas-density number 3.5 → 5.0 % (before that test was replaced outright). Fix a
single convention,
(model − measured)/measured, state it on the page, and use it everywhere.

Figures 7 (column diameter) and 9 (liquid properties) are still undigitised, so
the diameter and liquid-property independence claims remain untested.

## Hard-won lessons — read before building

The first two are also stored as user memories and will load automatically.

- **Never fabricate data.** If a dataset cannot be obtained, mark it
  `status: placeholder`, keep the page `status: planned`, and say so on the page.
  Re-simulating the model and presenting it as data is circular and would
  discredit the gallery.
- **Reported numbers must be deterministic.** On `B1.1` the fold was first
  located on a warm-start continuation curve; CI re-executed on another machine
  and got η_ignited 44.45 against 36.15 locally. Locate features on a smooth
  deterministic reference (analytical, or a shooting solve) and rank turning
  points by prominence.
- **Validation catches what inspection does not.** Three defects in `B1.1` were
  found by the monotonicity and closure checks, not by reading the code: a sign
  error, a sweep that missed a solution branch entirely, and a comparison that
  reported 89% deviation where the truth was 1e-5.
- **`check_agreement.py` failing is a signal, not a nuisance.** Update the
  baseline only when the model genuinely changed, and say why in the commit.
- **Look for a check the *paper* pays for.** `C2.1` gained three that cost
  nothing: Tables 5 and 6 are related by the Arrhenius form, so recomputing one
  from the other tests the page-image reading *and* the split reference
  temperature at once; the van 't Hoff slope of an equilibrium correlation is a
  reaction enthalpy, so it can be checked against the paper's own table; and two
  figures plotting the same runs must pair up point-for-point across
  independently fitted axes. Look for these before writing the notebook — they
  are worth more than any amount of code review.
- **Two independent methods beat one method plus a tolerance.** On `D2.2` the
  critical inlet pressure is computed twice: by bisection on the pymrm reactor,
  and by an adaptive quadrature of the trajectory in the *p*–*T* phase plane
  that never forms the reactor grid. They agree to 0.18 % across the whole
  operating range, and that number is the only one on the page that does not
  involve the paper. Look for a second route to the same quantity — it is often
  cheap, and it catches discretisation error that grid refinement alone will
  flatter.
- **Marker extraction is per-figure, not a recipe.** Morphological opening
  worked on Duncan & Toor because its markers are solid and much larger than the
  curves. Xu & Froment's are ~20 px glyphs on ~10 px curves and needed a
  different method: local ink density minus the largest value explainable by a
  locally straight structure (grey-scale opening with a long *line* element,
  maximised over six orientations). Take the maximum over orientations, or every
  steep near-origin curve section is flagged as a marker. Then audit every
  candidate visually at 600 dpi — on `C2.1` roughly a quarter of the automatic
  candidates were curvature artefacts with no glyph at the crosshair.

## Working commands

```bash
cd ~/Code/pymrm_suite/pymrm-gallery
source ~/Code/pymrm_suite/.venv/bin/activate

python scripts/check_metadata.py            # metadata + provenance validation
python scripts/check_metadata.py --report-missing   # catalog IDs not yet in models.yaml
python scripts/run_pages.py                 # execute every page notebook (slow: ~minutes)
python scripts/check_agreement.py           # metric regression check

export QUARTO_PYTHON=$(which python)
quarto render                               # build the site into _site/
```

Pages are generated from builder scripts rather than hand-edited JSON.
`pages/C2.1-xu-froment-smr/build_page.py` is the reference example: a list of
`nbformat` markdown/code cells written to `index.ipynb`, so prose and code stay
in one reviewable file. Nine required sections, listed in
[`AGENTS.md`](../AGENTS.md).
