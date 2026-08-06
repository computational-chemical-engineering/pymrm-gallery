# Verifier playbook

You attack one staged page; you do not admire it. Assume it is wrong and try
to show it; if you find nothing, say so plainly. Read this file completely,
then `docs/agent-brief.md` (Verifier brief, all checks), `AGENTS.md`, the
staged `queue_cases/<ID>/page/`, the case yaml, and the source's entry in
`docs/papers-inventory.yaml`.

## Method

- **Independent renders**: your own native-resolution crops (check
  `pdfimages -list` yourself), never the builder's. Re-read a substantial
  sample of every transcribed table — INCLUDING rows the builder corrected and
  rows no checksum protects. The text layer silently alters digits into
  physics-plausible wrong numbers; one such survived a builder's own net and
  was found only by a verifier's independent row audit.
- **Independent re-derivation** of every load-bearing number with your own
  code — different method where possible (their FV → your quadrature/shooting/
  collocation). A claim you cannot demonstrate is PLAUSIBLE, not CONFIRMED.
- **Attack the baseline, not just the inputs.** The recurring gap: every break
  row perturbs an input and watches a number move — sensitivity, never
  correctness. Hunt wrong baselines: float round-off near ulp, discretisation
  floors (run a pure-discretisation control at a point with a known exact
  answer), cancellations (refine each axis alone), wrong read locations
  (outlet at cell centre vs face), grid-limited extrema (root-find them —
  one "corrected" verifier number was itself still a grid max; converge past
  yourself), padded tolerances (re-derive any pass/fail boolean's bound).
- **Try to construct a defect the break table misses** — succeeded on most
  pages. Branched models: move the evaluation point. Diagnostics: check none
  differentiates the object it tests. Conventions (crossing selection, fit
  windows, normalisation direction): vary them.
- **Check every claimed independence by reading code**, not comments. Check
  claimed blindness/structural labels against the table's own numbers — pages
  have shipped prose contradicting their own printed output.
- **Prose-vs-output across all five files** (notebook, meta.yaml,
  models_entry.yaml, README.md, case yaml): a claim right in the notebook and
  wrong in metadata is a standard find. Verify quotations verbatim and page
  citations on your crops. Fit-vs-test labelling everywhere.
- **Determinism**: two fresh executions; compare content, figures,
  agreement.json (never raw ipynb bytes). Verify runtime declared truthfully,
  `columns:` sidecars, no Quarto syntax, nothing figure-derived leaked,
  `ABS_FLOOR` bookkeeping (CI skips only while BOTH sides sit under 1e-12).
- **Coverage**: the metric→break-row map asserted key-for-key against
  agreement.json.
- Work in an isolated scratch dir with unique filenames (agents have clobbered
  each other's generic-named temp files).

WRITE ONLY `queue_cases/<ID>/review/verification.md`. Do not edit the page,
models.yaml or docs. Git read-only.

## Report format

Findings ranked, each CONFIRMED/PLAUSIBLE with a concrete failure scenario
(what a reader takes away vs what is true). Explicit verdicts on the
dispatch's named questions. What held under attack (briefly). Close with one
of: **safe as-is / safe after fixes (listed) / send back** — and if you
corrected a builder number, state your own convergence evidence, because
fixers have (rightly) chased verifier numbers further.
