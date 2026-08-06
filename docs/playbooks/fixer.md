# Fixer playbook

You apply a verification's fixes to one staged page. Read this file, then the
FULL `queue_cases/<ID>/review/verification.md`, then `AGENTS.md` and the
builder playbook's non-negotiables (they bind you too).

## Rules of engagement

- **Verify each finding yourself before changing anything.** The verifier is
  not automatically right; fixers have correctly overruled verifiers with
  measurements several times (a bound the verifier "corrected" that was still
  a grid max; a consequence refuted by the page's own cancellation; a proposed
  fix that would have left 40 % of the contamination). If you disagree, argue
  with numbers in your report and do what the numbers say.
- **Do not weaken what survived.** The verification's "held under attack" list
  is load-bearing; keep those results at full strength. If a finding
  *strengthens* the page (an offered upgrade), take it only if your own
  measurement confirms it — and decline offered upgrades your measurement
  refutes.
- Numbers changed by a fix propagate to every file that quotes them — the
  standard set is notebook prose, meta.yaml, models_entry.yaml, README.md,
  case yaml, sidecars. Grep for the old value.
- A defect the break table missed becomes a new break row (inject it
  verbatim) or a named blind-spot entry — the bug you just fixed is the
  page's best teaching exhibit.
- Extrema/thresholds: root-find, show refinement evidence, never sample.
- Every number in prose printed by code; wording fixes must not introduce a
  typed number.
- Rebuild via `build_page.py`, re-execute, confirm two consecutive runs give
  identical content/figures/agreement.json. Run `scripts/check_metadata.py`.
  Leave the case `ready`; nothing into `pages/`; no git writes; touch nothing
  outside `queue_cases/<ID>*`.

## Report format

Per finding: what changed and the new numbers. Where you overruled the
verifier, the measurement. New/changed agreement metrics. Determinism
confirmation. Runtime delta if significant.
