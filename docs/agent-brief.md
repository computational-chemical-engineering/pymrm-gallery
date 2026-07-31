# Agent briefs

Two roles: a **builder** takes one case from `queue_cases/` to a staged page, and
a **verifier** attacks that page before it is published. Dispatch prompts should
point here rather than restate it, so the rules stay in one place.

The gallery optimises for **being right**, not for throughput. A missing page
costs nothing; a page with a fabricated number costs the credibility of every
other page.

---

## Builder brief

### Read first, in this order

1. `AGENTS.md` — the page contract and house style. Follow it exactly.
2. `docs/handoff.md` — accumulated lessons. Several will apply to any case.
3. `docs/pdf-findings.md` — per-paper extraction traps and acquisition routes.
4. The closest published page under `pages/`, as a model.

### Choose the validation before choosing the model

**This is the single most important decision, and it is made before any code.**

Rank the ways this paper can be checked, and take the highest available:

1. **A worked example with printed intermediates.** Best of all — `E2.1`
   reproduced three appendices to 0.46 % and never touched a figure.
2. **An internal identity the paper must satisfy.** Continuity, a limit the model
   collapses to, two tables related by a transformation. `J4.8` proved ASM1's
   continuity as symbolic identities; `B3.1` derived eq. 6 independently and
   matched to 6.9e-16.
3. **A stated numerical result in the text.** `I1.2` hit all six numbers of
   Table III; `F2.3` used the two conversions the authors quote.
4. **A digitised figure.** *Last resort.* It is the only route that needs a human
   review, so it is the only route that stalls a case.

If (1)–(3) exist, **do not digitise a figure at all**. Say in the page that the
data are the paper's own tabulated values and the tier is 6.

### Digitisation, when unavoidable

- Ask **what carries series identity** — shape, position, or a curve — before
  building a classifier. On `F1.4` the correlation had no term distinguishing the
  series, so labels were decoration; on `F2.3` identity came from which curve a
  marker sat on, not its matched shape.
- Curves whose equations are printed **calibrate the axes exactly**, and can be
  erased at their computed position before marker detection.
- Compute the effect you are testing **in pixels**. If it is smaller than the
  line width, the figure cannot test it and agreement there is not evidence.
- Always produce a numbered overlay for review, in `queue_cases/<ID>/review/`.
  These are **git-ignored** — they are the copyrighted page image.

### Make sure your validation can fail

The most common defect found in review is **an agreement that is algebraically
guaranteed**, presented as evidence. It always looks like the page's strongest
result. Four of six pages carried one in the 2026-07-31 batch — see
[`handoff.md`](handoff.md#the-check-that-cannot-fail).

Before putting an agreement number on a page:

- **Do the two routes share code?** Same assembly or same operator means you are
  testing arithmetic, not physics.
- **Break something on purpose** — a flipped sign, the wrong `nu`, a mismatched
  boundary condition — and confirm the number *moves*. If it does not, the check
  is decoration. Do this; it is cheap and nothing else substitutes for it.
- **Is the residual structural?** Conservation checks are often exact by
  construction.

Keeping such a check is fine — label it as the identity it is, and say what it
cannot detect.

**Interpolate computed values into your prose; never retype them.** Four pages in
that same batch stated a number the notebook contradicted two cells above.

### Rules that do not bend

- **Never fabricate data.** No synthetic points, no textbook restatement passed
  off as the source's. If a dataset cannot be obtained, halt.
- **Never repair a mangled number by inference.** Read constants off 600 dpi
  renders (`pdftoppm -r 600 -f N -l N -png f.pdf out`) and *look at the image*.
- **Reconstruction is allowed; fabrication is not.** The line: every input traces
  to something printed. `J3.4`'s conductivity was not printed, so it was obtained
  by inverting two other published equations — and the page says so.
- **Distinguish validated from reproduced.** Measurement is validation. The
  authors' own simulation output is reproduction. Never blur them.
- **When a printed constant looks wrong, prove it from the paper's own results**,
  and print the alternatives. `F2.3` established chemical control *first*, so
  nothing else could absorb a rate error, before concluding two constants were
  mis-set.

### Halt rather than guess

Halt on: a constant you cannot read; a figure you are not confident in; a stated
result you cannot reproduce and cannot explain; a scope call.

Halting is a **good outcome** — but only if the case is resumable. Write into
`queue_cases/<ID>.yaml`:

```yaml
status: needs-input          # or needs-paper
blocker: {kind, question, detail, artifacts}
resume:
  staged_at: queue_cases/<ID>/page/
  established: what is already done, so nobody re-derives it
  answer_changes: {each plausible answer -> exactly what changes}
  files_to_touch: [...]
  do_not_redo: ...
```

The question must be **specific and answerable without reading the paper**.

### Where to work

Only `queue_cases/<ID>/page/` and `queue_cases/<ID>.yaml` (+ `review/`). Never
`models.yaml`, `docs/`, `scripts/`, or any existing `pages/*`. **Never run git.**

Also emit `queue_cases/<ID>/models_entry.yaml` — the block to splice into
`models.yaml`, with `id, slug, title, description, section, tier, priority,
structures, phases, predicts, pymrm_api, reference{authors,year,container,doi},
data{tier,method,status}, related`. Keep `slug` and `title` identical to your
`meta.yaml`; mismatches are the most common integration failure.

### Before reporting ready

Regenerate and execute:

```
python build_page.py
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb=nbformat.read("index.ipynb",as_version=4)
NotebookClient(nb,timeout=1800,kernel_name="python3",resources={"metadata":{"path":"."}}).execute()
nbformat.write(nb,"index.ipynb"); print("OK")
EOF
```

Keep runtime under ~5 minutes. Report the status you set, the key numbers, which
checks passed, and anything the integrator must know.

---

## Verifier brief

You are given a **staged page and its source paper**. Your job is to attack the
page, not to admire it. Assume it is wrong and try to show it.

The failure this exists to catch is not sloppiness — it is **plausible,
confident agreement that is actually circular or coincidental**. That has
happened in this repository: a gas-density independence test compared two groups
that had no overlap in the abscissa, so it measured velocity, not density, and it
looked clean.

Check, in order:

1. **Traceability.** Every number on the page must come from a printed source, a
   computation shown, or an explicitly labelled reconstruction. Any number you
   cannot trace is a finding. This catches real defects: `J3.5` asserted a
   utilisation of 0.31 that appears in no source — not in the cited page, not in
   the original, and not in its own output, which gives 0.4006.
2. **Powerless checks.** For every agreement number, ask whether it *can* fail:
   do the two routes share code, is the residual structural, does the number move
   when you deliberately break something it should catch? Break it and see. This
   is now the most common finding in the repository — details and four worked
   examples in [`handoff.md`](handoff.md#the-check-that-cannot-fail).
3. **Prose against output.** Re-read every number in the markdown against what
   the cells print. Drift between the two has been found on four pages.
4. **Borrowed claims.** If the page reuses a sentence, number or dataset from
   another page, verify it against the *source*, not the sibling page — a wrong
   band on the published `F1.4` was found exactly this way.
5. **Circularity.** Is any "agreement" using a quantity that was fitted to the
   thing it is being compared against? Was a constant tuned until the comparison
   worked? Watch for a knob fitted to one identity that also moves the quantity
   used to validate — `F3.5`'s scalar did both.
6. **Confounding.** For any group comparison, do the groups overlap in every
   other variable? If not, the difference may be the confound.
7. **Resolving power.** Is the claimed effect larger than the figure's line width
   or the digitisation error? If not, agreement is not evidence. Ask the same of
   a headline agreement: on `H1.4` the +0.13 % did not move at all when the rate
   constant was varied 100-fold, because the answer sat on an algebraic ceiling.
8. **Direction and convention.** Deviations defined consistently? Boundary
   conditions on the outward normal? Reciprocals at >5 % scatter are not
   interchangeable.
9. **The honesty of the caveats.** Does the page overstate? Does it call
   reproduction "validation"? Does it report evidence that points *against* its
   own conclusion — `F3.5` had two such numbers and presented them as support.

Report findings ranked by severity, each with a concrete failure scenario. If you
find nothing, say so plainly — do not invent findings to look useful.
