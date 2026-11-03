# Brief templates

Fill in the placeholders and write the brief to the scratchpad; the dispatch prompt is
then one line: "Read and follow the brief at <path>. It is self-contained. Return the
report it asks for." Placeholders:

| placeholder | meaning |
| --- | --- |
| `<PAGE>` | page directory name, e.g. `F2.3-slurry-bubble-column-ft` |
| `<ID>` | catalog ID, e.g. `F2.3` |
| `<WT>` | the checkout the agent works in (a worktree of the sweep branch) |
| `<BRANCH>` | the branch checked out there |
| `<VENV>` | venv built from requirements-lock.txt (pymrm 3.4.0, jupytext) |
| `<SCRATCH>` | the coordinator's scratchpad directory |
| `<id>` | lower-case short ID for file names, e.g. `f23` |
| `<ROWS>` | the learnings row numbers the page added |

## A. Refresh, model-solving lane (`implement`, `model: opus`)

```
Refresh the gallery page pages/<PAGE> in the git checkout <WT>, branch <BRANCH>, onto
the new page layout and pymrm 3.4.0, in both lanes. Do not commit, do not `git add`, do
not `git rm` (delete files with plain `rm`): other agents may share the git index. Touch
only pages/<PAGE>/, the <ID> entry in models.yaml, and new rows appended at the end of
docs/refresh-learnings.md (re-read each shared file right before editing it).

Read first: .claude/skills/gallery-refresh/references/refresh-checklist.md (apply every
rule; the review will check against it), AGENTS.md (whole file), docs/page-infrastructure.md,
docs/claims-design.md, templates/index.py, templates/bootstrap_cell.py, and a refreshed
page of the same lane as worked example (<two recent pages, e.g. F2.3, J1.5>). The pymrm
rules: ~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/references/pitfalls.md
(complete), style-guide.md and api-map.md (the current API: which call does what), and the
plugin exemplar closest to this page's structure (conventions/exemplars/: dispersion_reactor,
steady_pellet, steady_tvd_channel for ConvFluxTVD and {"outflow": True}, assembly_styles;
model-patterns/exemplars/: membrane_reactor, darcy_gas_reactor, nested_particle_bed for
BlockLayout, construct_coefficient_matrix and coupled domains). Search docs/refresh-learnings.md for rows on this page's
section or structure; do not read all of it.

Lane A (infrastructure): index.py from templates/index.py, paired index.ipynb via
jupytext; bootstrap cell verbatim with GALLERY_REF="main", PAGE the directory name,
PYMRM from meta.yaml; meta.yaml pymrm_min_version "3.4.0" and realistic runtime_seconds;
model.py with NumPy-style docstrings (classes and functions only, no plotting, no data
loading); test_model.py (import, one fast solve, one stored value via load_claims, runs
offline in seconds); every number the page states is a Claims entry, last code cell
`claims.report();`; agreement.json format 2; build_page.py removed; README
`## Corrections` list for any published statement that was wrong (dated today, keep
existing dates).

Lane B (pymrm 3.4.0): <what the old page computes, with build_page.py line numbers of
the solver, root finds and outlet reads; the candidate changes>. Whichever route you
choose, another may stay only as an independent validation route described by what it
computes, with its measured agreement. Use the operators, assemble constant operators
once, spatial axes first, NumJac((n,1)) never (n,). Comment every bc with the physical
equation (outward normal). Published numbers may change when you show and verify why.
If the problem is linear, say so and do not add newton for show.

Python: <VENV>. Run long commands in the background and wait with `kill -0 PID` or the
harness notification, never a `pgrep -f` loop. Old outputs: `git show
origin/main:pages/<PAGE>/index.ipynb`; reproduce every old published value you change on
the current venv first (old code from build_page.py).

Progress file (required): <SCRATCH>/<id>_progress.md after every checkpoint
(measurements, decisions, corrections found, what is left). If you receive "pause",
finish the command in hand, update the progress file and stop.

check_model.py rewrites agreement.json: copy it aside before running it and restore it.
A format-1 baseline on main shows drift FAILs for moved metrics before the commit; list
each with its cause. Keep metric names from the old agreement.json where the quantity
is the same. If the baseline on main is already format 2, removed break labels need
docs/loosened.yaml entries.

Acceptance (report each):
1. `python scripts/run_kernels.py pages/<PAGE>`: every claim passes, printed text
   identical under the three kernels.
2. `python scripts/check_agreement.py --base origin/main`: list every metric that moved,
   old and new value, and why.
3. scripts/check_claims.py, `pytest tests/`, and test_model.py offline (GALLERY_OFFLINE=1).
4. ~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/scripts/check_notebook_math.py
   on index.ipynb and check_model.py on model.py.
5. `python scripts/check_metadata.py` OK.
6. `jupytext --sync` leaves index.py unchanged.

Report at most 400 words: files changed, lane B choice and its measured justification,
moved metrics, changed published statements (each with evidence), claim and break-row
counts as printed, acceptance results, runtime.
```

## A-light. Refresh, light lane (`implement`, Sonnet)

Same as A with these differences:

- Replace Lane A's model.py and test_model.py sentence with: "No model.py and no
  test_model.py unless a check script requires them; if one does, stop and report which
  check and why instead of inventing a model."
- Lane B: "No pymrm solve (<what the page computes>). Every root find with a hand
  bracket (brentq, fsolve; build_page.py lines <...>) becomes pymrm.checks.find_roots on
  the physical range, raising unless exactly one root is found; every extremum is a root
  of an analytic derivative; an argmax over data rows stays, naming the row. Do not add
  other pymrm calls. Keep 'What pymrm adds' honest: say plainly that pymrm contributes
  nothing to the correlation, if so." Name the sibling pages whose claims it may load
  (never retype another page's claims).
- Acceptance drops test_model.py and check_model.py.

## B. Review (`gallery-reviewer`)

```
Review the refreshed gallery page <ID> (<title>) before it is committed. Advice only.
Scratch: <SCRATCH>/review_<id>/ (create it; notes in NOTES.md).
Checkout: <WT> (branch <BRANCH>). Venv: <VENV>. Page: pages/<PAGE>/ (<lane>).
Old version: `git show origin/main:pages/<PAGE>/build_page.py`, old README, meta.yaml,
index.ipynb on origin/main. Implementer notes: <SCRATCH>/<id>_progress.md. Other
uncommitted changes: models.yaml <ID> entry, learnings rows <ROWS>.

Implementer's findings to verify: <the report restated as facts, with numbers>.

What to check, in this order:
1. <the riskiest judgement call: a reversed conclusion, a changed reading, a cited page>
2. Each Correction: real error or change of definition / fair rounding; old wording
   quoted exactly; what is true stated with claim markers.
3. <lane-specific numerics: routes, orders, roots, break rows>
4. The checklist, as a final pass.

Cheap checks you may run: <list, with runtime>.
Output (at most 400 words): numbered items for a Sonnet fix round, each MUST-FIX or
optional, with file:line, what is wrong, and the exact fix (replacement wording where
wording is the issue). Then one verdict line per question asked above.
```

## C. Fix round (`gallery-fixer`)

```
<ID> fix round (from the review, <date>).
Checkout: <WT> (branch <BRANCH>). Venv: <VENV>. Page: pages/<PAGE>/.
Review notes: <SCRATCH>/review_<id>/NOTES.md. Implementer notes: <SCRATCH>/<id>_progress.md.
Progress file: <SCRATCH>/<id>_fix_progress.md. Touch only the <ID> page, the models.yaml
<ID> entry and learnings rows <ROWS>. Correction dates stay <date>.
<anything that must not change, e.g. claims another page loads>

Items (line numbers approximate):
1. MUST-FIX. <file:line>: <what is wrong>. <exact fix or replacement wording>.
...

Acceptance: the standard list in your agent definition, for pages/<PAGE>
(<lane>; expected non-OK check_agreement lines: <...>).
Report at most 250 words: one line per item, acceptance results, new claim values.
```
