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
| `<FACTS>` | `<SCRATCH>/<id>_facts.md`, written by `scripts/page_facts.py <ID>` |
| `<REF>` | `<SCRATCH>/<id>_old`, the reference run of the old page (`accept_page.py --make-ref`) |
| `<WORKED>` | ONE refreshed page of the same lane, from the facts file |
| `<EXEMPLAR>` | the plugin exemplar closest to the page, and the part of it that matters |
| `<ROWS>` | the learnings row numbers the page added |

The plugin exemplars: `~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/exemplars/`
(dispersion_reactor, steady_pellet, steady_tvd_channel for ConvFluxTVD and
`{"outflow": True}`, assembly_styles, refinement_orders, correlation_roots,
continuation_s_curve, cstr_steady_states) and `.../model-patterns/exemplars/`
(membrane_reactor, darcy_gas_reactor, nested_particle_bed for BlockLayout,
construct_coefficient_matrix and coupled domains, slurry_bubble_column).

## The working-style block (every brief carries it verbatim)

```
Working style. Every call pays again for everything read so far, so:
- Read by line range (Read with offset and limit, or `sed -n 'a,bp'`), after a grep
  that finds the range. Never cat a whole file of more than about 200 lines; never
  read a file a second time; never print a notebook.
- Look rules up when you need them: `grep -n '^#' <file>` for the headings, then read
  the one section. The files: AGENTS.md, docs/claims-design.md,
  docs/page-infrastructure.md, and in ~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/
  conventions/references/: api-map.md (which call does what), pitfalls.md,
  style-guide.md, numerics.md. The learnings rows: `grep -n '^| <row> |'
  docs/refresh-learnings.md`.
- Change files with Edit. Write only creates a new file; never rewrite a whole file
  to change part of it.
- Run anything longer than a minute in the background with output to a log file, and
  wait with `kill -0 PID` or the harness notification, never a `pgrep -f` loop. Read
  only `tail` or a grep of the log.
- Check with `python scripts/accept_page.py pages/<PAGE> --quick --only <steps>`
  while iterating (steps: sync, kernels, agreement, claims, plugin, metadata, stable,
  dashes, lint, outputs, port). Run the full acceptance once, at the end, and open its
  log only at a FAIL, by line range.
- Keep the progress file short: measurements, decisions, what is left.
```

## A1. Port, model-solving lane (`implement`, Sonnet)

```
Port the gallery page pages/<PAGE> in the git checkout <WT>, branch <BRANCH>, onto the
new page layout WITHOUT changing what it computes. A second agent does the numerics,
the claims and the wording afterwards, starting from your progress file. Do not
commit, do not `git add`, do not `git rm` (delete files with plain `rm`): other agents
may share the git index. Touch only pages/<PAGE>/.

Read: <FACTS> (the old page's definitions, solver lines and metrics);
templates/index.py and templates/bootstrap_cell.py; AGENTS.md section "Adding a page"
only; pages/<WORKED>/index.py and model.py for the layout only, by line range. The
old page is pages/<PAGE>/build_page.py: its notebook cells are the strings passed to
`cells.append(...)`; read it in ranges.

<the working-style block>

Do:
1. index.py from templates/index.py, carrying the old cells in order: prose as
   Markdown cells, code as code cells. Bootstrap cell verbatim from
   templates/bootstrap_cell.py with GALLERY_REF = "main", PAGE = "<PAGE>",
   PYMRM = ">=3.4.0,<4". Pair it: `jupytext --set-formats ipynb,py:percent index.py`.
2. model.py: the classes and functions the old cells define that solve or evaluate
   the model, moved unchanged apart from a NumPy-style docstring each; no plotting,
   no data loading. index.py imports them.
3. meta.yaml: pymrm_min_version "3.4.0"; runtime_seconds measured.
4. test_model.py: import, one fast solve, compare one value with the page's
   agreement.json metric of the same name at rtol 1e-9; offline, seconds. (A2 moves
   it to `load_claims`.)
5. Keep the old last cell's `report_agreement("<ID>", {...})` call unchanged, so every
   old metric can be compared exactly with the reference run of the old page.
   Remove build_page.py with plain `rm`.

Do not change any solver call, grid, tolerance, bc, parameter, number or sentence; do
not add Claims, break rows, Corrections or README text. If a cell cannot move without
a change (a path, an import), make the smallest change and list it.

Acceptance (report each):
1. `python scripts/accept_page.py pages/<PAGE> --quick --only sync,kernels,port,stable --ref <REF>`:
   every line PASS. A port FAIL means a metric differs from the old page's run on the
   same kernel: find the cell that changed. Each LOOK line under port (a printed line
   of the old page missing or new) is a moved cell or the bootstrap, never a changed
   result; say which.
2. `python scripts/accept_page.py pages/<PAGE> --only claims`: PASS (test_model.py).
3. `python scripts/accept_page.py pages/<PAGE> --only plugin,dashes,metadata`: do NOT
   fix what these find (old em dashes, history words, pitfalls in the old code, the
   missing Claims): copy every FAIL line into the progress file for A2.

Progress file (required): <SCRATCH>/<id>_port_progress.md, after every step: what
moved where (old build_page.py lines to index.py cells and model.py functions), any
smallest change made, the acceptance lines, runtime. If you receive "pause", finish
the command in hand, update the progress file and stop.

Report at most 200 words: files, the mapping in one line per model.py function, the
acceptance lines as printed, runtime.
```

## A2. Numerics and claims, model-solving lane (`implement`, `model: opus`)

```
Refresh the gallery page pages/<PAGE> in the git checkout <WT>, branch <BRANCH>, in
lane B (pymrm 3.4.0) and the Claims. A port agent has already moved the page onto the
new layout without changing what it computes; its notes are
<SCRATCH>/<id>_port_progress.md and the page reproduces the old metrics through
`report_agreement`. Do not commit, do not `git add`, do not `git rm` (delete files
with plain `rm`): other agents may share the git index. Touch only pages/<PAGE>/, the
<ID> entry in models.yaml, and new rows appended at the end of
docs/refresh-learnings.md (re-read each shared file's end right before editing it).

Read first, and nothing else up front:
1. .claude/skills/gallery-refresh/references/refresh-checklist.md, all of it. Apply
   every rule; the review checks against it.
2. <FACTS>, and the port notes.
3. pages/<PAGE>/index.py and model.py, by range as you work on them.
4. pages/<WORKED>/: its Claims cell and break table in index.py, and its README
   Corrections, as the worked example.
5. <EXEMPLAR>.
6. ~/Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/references/pitfalls.md,
   complete (you rewrite solver code).
7. When the facts file lists borrowed data (`load_data(..., page=<another page>)`):
   AGENTS.md section "Data rules", and the source page's README findings on the rows
   you use.
Everything else (AGENTS.md, claims-design, page-infrastructure, api-map,
style-guide, numerics, learnings rows) by grep when a question needs it.

<the working-style block>

Lane B: <what the page computes, with old build_page.py line numbers from the facts
file for the solver, root finds and outlet reads; the candidate changes:
ConvFluxTVD with pymrm.newton (damping="monotone", stagnation=True, raise when not
sol.success) instead of deferred correction; construct_coefficient_matrix or
BlockLayout instead of hand indices; {"outflow": True} and compute_boundary_values
instead of reading the last cell; find_roots instead of a bracket>. Whichever route
you choose, another may stay only as an independent validation route described by
what it computes, with its measured agreement. Use the operators, assemble constant
operators once, spatial axes first, NumJac((n,1)) never (n,). Comment every bc with
the physical equation (outward normal). Published numbers may change when you show
and verify why; reproduce every old published value you change first (the old code
is in `git show origin/main:pages/<PAGE>/build_page.py`). If the problem is linear,
say so and do not add newton for show.

Claims: replace `report_agreement` by Claims format 2. Every number the page states
is a Claims entry, metric names kept from the old agreement.json where the quantity
is the same, the last code cell `claims.report();`, every non-structural claim with a
break row (checklist rules 1 to 11); test_model.py's stored value through
`load_claims`. README `## Corrections` for any published statement that was wrong
(dated today, keep existing dates); "What pymrm adds"; the models.yaml entry.

Python: <VENV>. check_model.py rewrites agreement.json; accept_page.py restores it.
A format-1 baseline on main shows drift FAILs for moved metrics before the commit;
list each with its cause. If the baseline on main is already format 2, removed break
labels need docs/loosened.yaml entries.

Progress file (required): <SCRATCH>/<id>_progress.md after every checkpoint
(measurements, decisions, corrections found, what is left). If you receive "pause",
finish the command in hand, update the progress file and stop.

Acceptance: `python scripts/accept_page.py pages/<PAGE>` (full, in the background),
once after the last edit. Every step PASS, except agreement, which is LOOK before the
commit when the baseline is format 1 (each of this page's lines explained: moved or
replaced metrics); every LOOK line under lint checked (tagged, quoted or an
identifier). Then the checklist's
"Before you report" items 2 to 4.

Report at most 400 words: files changed, lane B choice and its measured justification,
moved metrics, changed published statements (each with evidence), the accept_page
summary as printed, runtime.
```

When the port was not run (a page whose cells cannot be ported without numerics),
give A2 with this paragraph in place of the port sentence: "No port was run: first
move the page onto the new layout as brief A1 describes (<the A1 Do list>), check
that the old metrics are reproduced, then continue." This is the old single-agent
refresh; use it only when the split does not fit, and say so in the cost log.

## A-light. Refresh, light lane (`implement`, Sonnet)

```
Refresh the gallery page pages/<PAGE> in the git checkout <WT>, branch <BRANCH>, onto
the new page layout and pymrm 3.4.0. Do not commit, do not `git add`, do not `git rm`
(delete files with plain `rm`). Touch only pages/<PAGE>/, the <ID> entry in
models.yaml, and new rows appended at the end of docs/refresh-learnings.md (re-read
each shared file's end right before editing it).

Read first, and nothing else up front: the refresh checklist
(.claude/skills/gallery-refresh/references/refresh-checklist.md, all of it); <FACTS>;
templates/index.py and templates/bootstrap_cell.py; pages/<WORKED>/index.py (layout,
Claims cell, break table) and its README Corrections; the old build_page.py by range;
when the facts file lists borrowed data, AGENTS.md section "Data rules" and the
source page's README findings on the rows you use. Everything else by grep when a
question needs it.

<the working-style block>

Layout: index.py from templates/index.py, paired index.ipynb via jupytext; bootstrap
cell verbatim with GALLERY_REF="main", PAGE the directory name, PYMRM ">=3.4.0,<4";
meta.yaml pymrm_min_version "3.4.0" and measured runtime_seconds; every number the
page states is a Claims entry, last code cell `claims.report();`; agreement.json
format 2; build_page.py removed; README `## Corrections` for any published statement
that was wrong (dated today, keep existing dates). No model.py and no test_model.py
unless a check script requires them; if one does, stop and report which check and
why instead of inventing a model.

Lane B: no pymrm solve (<what the page computes>). Every root find with a hand
bracket (brentq, fsolve; build_page.py lines <...>) becomes pymrm.checks.find_roots on
the physical range, raising unless exactly one root is found; every extremum is a
root of an analytic derivative; an argmax over data rows stays, naming the row. Do
not add other pymrm calls. Keep "What pymrm adds" honest: say plainly that pymrm
contributes nothing to the correlation, if so. Sibling pages whose claims it may
load (never retype another page's claims): <...>.

Progress file, pause rule and acceptance as in A2 (accept_page.py full, once at the
end; it skips test_model.py and check_model.py when the page has none).

Report at most 400 words, as in A2.
```

## B. Review (`gallery-reviewer`)

```
Review the refreshed gallery page <ID> (<title>) before it is committed. Advice only.
Scratch: <SCRATCH>/review_<id>/ (create it; notes in NOTES.md).
Checkout: <WT> (branch <BRANCH>). Venv: <VENV>. Page: pages/<PAGE>/ (<lane>).
Old version: `git show origin/main:pages/<PAGE>/build_page.py`, old README, meta.yaml,
index.ipynb on origin/main. Facts: <FACTS>. Implementer notes:
<SCRATCH>/<id>_progress.md. Other uncommitted changes: models.yaml <ID> entry,
learnings rows <ROWS>.

Implementer's findings to verify: <the report restated as facts, with numbers>.

What to check, in this order:
1. <the riskiest judgement call: a reversed conclusion, a changed reading, a cited page>
2. Each Correction: real error or change of definition / fair rounding; old wording
   quoted exactly; what is true stated with claim markers.
3. <lane-specific numerics: routes, orders, roots, break rows>
4. The checklist, as a final pass.

Cheap checks you may run: `python scripts/accept_page.py pages/<PAGE> --only
dashes,lint,outputs` (seconds, changes no file); <others, with runtime>.
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

Acceptance: the standard acceptance in your agent definition, for pages/<PAGE>
(<lane>; expected non-OK agreement lines: <...>).
Report at most 250 words: one line per item, the accept_page summary, new claim values.
```
