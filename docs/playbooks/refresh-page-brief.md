# Briefs for the refresh sweep

The sweep (`docs/refresh-plan.md`) refreshes one page per cycle:

1. an `implement` agent on Opus refreshes the page in both lanes from brief A;
2. if a published statement or a metric beyond its tolerance changed, a `think`
   agent (Opus, effort xhigh) reviews it from brief B, without editing;
3. the review's findings go back to the same implementer (SendMessage to its id),
   quoted with the reviewer's numbers, followed by the same acceptance checks;
4. the coordinator spot-checks and commits one page per commit (section "Commit").

Fill in `<PAGE>` (directory name), `<ID>`, and the page-specific lane B notes. Keep
the lesson list whole: every item on it was a defect a review found on an earlier page.

## A. Refresh brief (implement, model opus, background)

```
Refresh the gallery page pages/<PAGE> in /home/eajfpeters/Code/pymrm_suite/pymrm-gallery,
branch <BRANCH>, onto the new page layout and pymrm 3.3.0, in both lanes. Do not
commit, do not `git add`, do not `git rm` (delete files with plain `rm`): other agents
may work in this checkout and share the git index. Touch only pages/<PAGE>/, the <ID>
entry in models.yaml, and new rows appended at the end of docs/refresh-learnings.md
(re-read each shared file right before editing it).

Read first: AGENTS.md (whole file, in particular "Adding a page", "A page shows current
pymrm, not history", "Honesty requirements", the data rules), docs/page-infrastructure.md,
docs/claims-design.md, templates/index.py, templates/bootstrap_cell.py,
docs/refresh-learnings.md, and a refreshed page as worked example: pages/A2.3-*/ and
pages/A2.5-*/ (index.py, model.py, test_model.py, README.md Corrections). The pymrm
rules: ~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/references/pitfalls.md
(complete) and style-guide.md.

Lane A (infrastructure): index.py from templates/index.py, paired index.ipynb via
jupytext; bootstrap cell verbatim with GALLERY_REF="main", PAGE the directory name,
PYMRM from meta.yaml; meta.yaml pymrm_min_version "3.3.0" and realistic
runtime_seconds; model.py with NumPy-style docstrings (classes and functions only, no
plotting, no data loading); test_model.py (import, one fast solve, one stored value via
load_claims, runs offline in seconds); every number the page states is a Claims entry,
last code cell `claims.report();`; agreement.json format 2; build_page.py removed;
README `## Corrections` list for any published statement that was wrong (dated today,
keep existing dates).

Lane B (pymrm 3.3.0): <what the page does now, and the candidate change: ConvFluxTVD
with pymrm.newton (damping="monotone", stagnation=True, raise when not sol.success)
instead of deferred correction; construct_coefficient_matrix or BlockLayout instead of
hand indices; {"outflow": True} and compute_boundary_values instead of reading the last
cell>. Whichever route you choose, another may stay only as an independent validation
route described by what it computes, with its measured agreement. Use the operators,
assemble constant operators once, spatial axes first, NumJac((n,1)) never (n,). Comment
every bc with the physical equation (outward normal). Published numbers may change when
you show and verify why.

Lessons from earlier pages (apply all):
- Refining a read-off can expose a same-sign discretisation error; a spread is not an
  error bound unless the ladder is monotone; measure every "converged" claim (print the
  ladder), over EVERY case a statement covers, never a hand-picked subset.
- A ratio of two computed values uses matched, refined grids for both.
- Close residuals exactly where possible; NumJac is not exact.
- Do not call two routes "identical" or independent when they share a factorisation or
  iteration; measure a second permc_spec if you claim solver independence.
- When you say what breaks an identity, test the alternatives one at a time (scheme vs
  time integrator vs limiter) before naming the cause.
- Print at a precision round-off supports; round-off-sized values are bounds with
  `(here ...)`.
- Check each identity is real (linearity can make a sweep metric structural;
  zero-gradient upwind ends leak mass). Label structural claims as such.
- Counts on the page are printed by code. Break rows must inject the defect itself
  (through an argument or physical parameter), and their labels are built from the
  calls. Every metric needs a break row that moves it.
- Never retype CSV values or another page's claims: compute or load them, and
  reconcile printed values against the CSV.
- Each correction must be verified (on pymrm 2.3.1 if needed to reproduce the old
  statement); no pymrm version names, no history narration on the page, in comments,
  meta.yaml or models.yaml.
- Compare on the same basis. "What pymrm adds" must be truthful. No em dashes anywhere.

Python: a venv from requirements-lock.txt (pymrm 3.3.0, jupytext included). Run long
commands in the background.

Acceptance (report each):
1. `python scripts/run_kernels.py pages/<PAGE>`: every claim passes, printed text
   identical under the three kernels.
2. `python scripts/check_agreement.py --base <BASE>`: no regressions; list every metric
   that moved, old and new value, and why.
3. scripts/check_claims.py, `pytest tests/`, and test_model.py offline (GALLERY_OFFLINE=1).
4. ~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/scripts/check_notebook_math.py
   on index.ipynb and check_model.py on model.py.
5. `python scripts/check_metadata.py` OK.
6. `jupytext --sync` leaves index.py unchanged.

Report at most 400 words: files changed, lane B choice and its measured justification,
moved metrics, changed published statements (each with evidence), claim and break-row
counts as printed, acceptance results, runtime.
```

`<BASE>` is the branch the baselines are compared against; it was
`refresh/pymrm-3.2.0` through the transport sweep, `origin/main` works once the page's
old agreement.json is on main.

## B. Review brief (think, background, no edits)

Give it: the page directory and the git commands to see the old version
(`git show HEAD:pages/<PAGE>/build_page.py`, `git diff HEAD -- ...`), the rules to read
(AGENTS.md sections, docs/claims-design.md), the venv, the warning not to re-execute
index.ipynb in place, the implementer's report restated as facts, and this checklist
in priority order:

1. Each correction: is the new statement true and the old one really wrong? Is the
   evidence a refined ladder over every case it covers?
2. Do claims and break rows test what their labels say? Any "converged", "identical",
   "agrees" without a measurement, or digits beyond round-off?
3. Is the second route independent, and is "What pymrm adds" truthful?
4. Retyped numbers, history narration, em dashes.
5. pymrm conventions in model.py.

Ask for a numbered list of defects (file:line, what, evidence, fix, must-fix or
optional), one line "commit-ready yes or no", at most 600 words.

## Commit

While any agent works in the checkout:

- `git diff --cached --stat` must list nothing you did not stage, or commit with
  `git commit -- <paths>`.
- Stage the page directory, `git rm` the old build_page.py, and this page's part of
  models.yaml and docs/refresh-learnings.md. When another page's agent has edited the
  same files, use `python scripts/stage_page_part.py <ID> <row,row,...>` (an agent
  has no interactive `git add -p`).
- Spot-check before committing: `git status`, the README Corrections, `grep -n '—'`,
  a grep for history words (earlier, previously, used to, old, version names),
  `scripts/check_metadata.py`.
- Message: `<ID> on the new page layout and pymrm 3.3.0; published corrections`, with
  the moved metrics and the corrections listed in the body.
