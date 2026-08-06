# Builder playbook

You build one gallery page for one catalogue case. Read this file completely,
then: `AGENTS.md`, `docs/agent-brief.md` (Builder brief), the case's
`queue_cases/<ID>.yaml`, and the source's entry in `docs/papers-inventory.yaml`
(identity, native ppi, text-layer warnings, known traps — do not re-derive
these). Skim `docs/handoff.md`'s header for anything newer than this file.

## Non-negotiables (each one is a defect class found on a published page)

1. **Read the source on native-resolution crops.** Check `pdfimages -list`
   yourself (disk holds 150–600 ppi). Render at native; crop every numeric and
   read it at digit scale. Never trust the text layer for a digit: it silently
   alters digits into plausible wrong numbers, prints European thousands
   separators, drops mid-dot decimals, and sometimes opens with the *previous*
   article. Verify the file's identity from its own title page even though the
   inventory names it.
2. **Fit vs test, labelled everywhere.** If the correlation was fitted to the
   data you compare against, the agreement is a goodness of fit — say so in the
   notebook AND meta.yaml AND models_entry.yaml AND README.md AND the case
   yaml. Compute a null baseline beside any headline agreement.
3. **Every `agreement.json` metric needs a break row that moves it.** Where
   impossible, label the metric structural and state what it cannot detect.
   Assert the coverage map against `agreement.json` key-for-key. Metrics below
   `ABS_FLOOR = 1e-12` are outside CI while both sides stay under it — name
   them and give each an above-floor companion.
4. **Compute at least one headline a second, independent way** (closed form,
   different solver family, quadrature — sharing no assembly). A break row
   shows sensitivity, never correctness; wrong-baseline defects (float
   round-off, discretisation floor, cancellation, wrong read location,
   grid-limited extremum) are invisible to every perturbation. **Root-find
   extrema and thresholds; never report a sampled max or a swept crossing** —
   two verified pages had "corrected" numbers that were still grid maxima.
5. **Every number in prose is printed by the code**, interpolated, never typed.
   Same for metadata where feasible; anything typed must match printed output.
6. **Evaluate branched/piecewise models on every branch and at the switch.**
   A check living on one branch lets constants on the other be deleted
   unnoticed. Exercise every printed exponent somewhere it is live; if one is
   inert or unexercised by the data, say so explicitly.
7. **Refine every axis that carries error** (grid AND time step, through any
   front/event), report observed orders. The unmeasured axis has repeatedly
   been the larger error.
8. **Boundary reads**: outlet values via `compute_boundary_values` — but know
   the scope: with a zero-gradient outflow bc left as-is that read is 2nd
   order; hand-written `v*C_N` last-cell reads are O(h) and fail mass balances
   while looking like physics errors (handoff has both configurations).
9. **Cross-page CSVs**: loading another page's dataset means reading that
   page's findings; never retype a number that is a row in a CSV you loaded —
   print it beside yours and reconcile.
10. **Figures**: content living only in figures is scoped out (no maintainer
    available) — or, only if load-bearing, extracted by the compute-and-erase
    method documented on `pages/A2.5-edwards-richardson-dispersion/`. Never
    trace a curve you can compute. No page images committed anywhere.
11. **Printed defects in the source are reported, never repaired** — prove
    them from the paper's own numbers (the F2.3 order: pin what is NOT free
    first). Repairs-by-inference must be labelled as inferences.
12. **Nothing fabricated, nothing from memory.** Origins you did not read go
    under `origin_not_consulted`/`reference_read_from` per AGENTS.md. For
    monograph-sourced (canonical-class) cases, the E1.1 test: the book must
    state, attribute and carry the result, not merely mention the topic.
13. Seed anything stochastic. Two consecutive executions must give identical
    content, figures and `agreement.json` (raw .ipynb bytes never match —
    nbformat mints cell IDs; compare content).

## Mechanics

- Copy the closest published page's directory and substitute the physics —
  current best exemplars: pointwise algebra `A1.6`/`A3.1`; 1-D transient
  `A2.8`/`B2.2`; 2-D `A3.15`; nested-scale `B3.2`; kinetics-fit `C1.1`;
  criterion audit `B1.7`; review-sourced `A4.7`/`A3.7`. **The break table does
  not travel when you copy — rebuild it for your physics.**
- Deliverables in `queue_cases/<ID>/page/`: `build_page.py` → `index.ipynb`
  (executed clean), `meta.yaml`, `README.md`, `data/*.csv` each with a
  `.meta.yaml` sidecar carrying a `columns:` block. Nine sections in order
  (AGENTS.md). Plus `queue_cases/<ID>/models_entry.yaml` — check whether
  models.yaml already has the case (planned entry → say "upgrade in place";
  else "append"; state which, correctly).
- Leave the case `ready`. Nothing into `pages/`. No git writes.
- Keep runtime modest and declare it truthfully in meta.yaml.
- Blocked? Park with a complete `resume:` block; never wait for the maintainer.

## Report format (keep it dense)

Status (`ready`/parked) · what the source prints (numbers vs figures) · the
scope decision if any · headline numbers · break-table summary + declared
blind spots · the second independent computation and what it can catch ·
printed defects found · what the page cannot conclude · runtime · anything
you disagreed with in these instructions, argued with numbers.
