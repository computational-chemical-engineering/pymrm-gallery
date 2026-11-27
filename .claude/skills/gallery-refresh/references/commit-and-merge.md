# Commit, PR and merge

## Staging when two pages share files

models.yaml and docs/refresh-learnings.md hold every page's edits. With one page
uncommitted, `git add -A pages/<PAGE> models.yaml docs/refresh-learnings.md`. With two:

```
git add -A pages/<PAGE>
python scripts/stage_page_part.py <ID> <row,row,...>   # this page's models.yaml entry and learnings rows
git diff --cached --stat                                # nothing you did not mean to stage
```

## Spot-checks before the commit (coordinator, not delegated)

- Read the README `## Corrections` yourself. Each entry: old wording in quotes and
  exact, what is true stated, every number with a `<!-- claim: name -->` marker. A
  phrase like "under those assumptions" must have its assumptions inside the same entry.
- `python scripts/accept_page.py pages/<PAGE> --only sync,agreement,claims,stable,dashes,contract,lint`
  (seconds): every step PASS, agreement PASS or LOOK with only the known pre-commit
  cases below, every LOOK line accounted for. Run it again after any edit of yours to
  index.py.
- After a fix round, read every prose sentence the fix diff added or changed against
  the page output: `git diff -U0 <before-fix> -- pages/<PAGE>/README.md
  pages/<PAGE>/meta.yaml pages/<PAGE>/index.py`, markdown lines only (`# ` in
  index.py). Count the new unsupported sentences in the cost-log row; if five fix
  rounds in a row add none, this step may be dropped.
- `--only contract`: `claims.report()` executed, check values asserted first, and
  `status: published` set by the branch (LOOK: the coordinator sets it).
- The claim count, break-table line and break-row count: `--only outputs` prints
  them from the executed notebook and agreement.json without running anything. The
  rows line must read PASS (as many recorded as defined); LOOK means the page does
  not declare its rows with `claims.define_rows`, which a refreshed page must, that
  no more than one row catches a claim, or that a row only records a crash (say so
  on the page, or make the row move a number).

## Learnings rows (coordinator, before the commit)

Every page closes with its learnings decided, also a new page and a second review.
Take the "Learnings" lines of the build or refresh report and of each review, keep
those the evidence supports, reserve row numbers with `gq.py rows <ID> <n>`, and append
them to `docs/refresh-learnings.md` with a target: a pymrm file (`numerics`,
`style-guide`, `verifier`, `pitfall`, `library`, ...) or `gallery: <where>`. Apply a
gallery row in the same commit and mark it `done`. If nothing qualifies, write
"learnings: none" in the page's cost-log notes. Then run
`python scripts/learnings_status.py --pymrm ~/Code/pymrm_suite/pymrm`: exit status 1
means ten rows wait for pymrm (flush them, refresh-plan section 5) or a row pymrm
already cites is still marked open (mark it).

## Known checker behaviour

- **Format-1 baseline on main**: before the commit, check_agreement compares against
  HEAD's format-1 file and shows drift FAILs or "metric disappeared" for moved or
  replaced metrics. After the commit the page shows `MIGRATED` and 0 regressions. No
  docs/loosened.yaml entry is needed or possible for these.
- **Format-2 baseline on main**: a removed break label or a loosened tolerance needs a
  docs/loosened.yaml entry with a reason (see the B2.3 follow-up commit).
- **check_model.py** rewrites agreement.json: back it up and restore it.

## Commit message

```
<ID> on the new page layout and pymrm 3.4.1; published corrections

<lane>: index.py/ipynb pair, Claims format 2 (<n> claims, <m> break rows, none open),
build_page.py removed, <model.py and test_model.py | no model.py>. <solver in one
sentence, runtime old -> new>.

<the main finding, if the page's conclusion changed>

Corrections (README, <date>):
1. <old statement in short>: <what is true>.
...

Moved metric(s): <name old -> new (why)>.

Reviewed by an independent think pass (<k> must-fix, <j> optional), <applied / which
dropped and why>.

<attribution lines>
```

## PR and merge

- Push batches: one PR for several pages, so one CI run (about 2.5 h) covers them. Do
  not push to the branch of a PR whose CI is running unless the commit belongs there;
  wait for the merge and push the next batch as a new PR.
- PR body: one bullet per page (lane, claims, the main finding, corrections count),
  then the shared checks (check_agreement 0 regressions, check_claims 0 problems,
  tests pass).
- Merge with `gh pr merge <n> --merge` when all four checks (metadata, notebooks
  locked, latest, portability) are green; the sweep owner has standing authorisation.
