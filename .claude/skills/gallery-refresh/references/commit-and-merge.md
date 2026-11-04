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
- `grep -n '—'` over the page files, the models.yaml entry and the new learnings rows.
- After any edit to index.py: `jupytext --sync`, then `scripts/check_claims.py`.
- `scripts/check_agreement.py --base origin/main` after staging: 0 regressions, or the
  known pre-commit cases below.
- Printed break-table counts, from the executed notebook outputs (the fixer may not
  report them).

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
<ID> on the new page layout and pymrm 3.4.0; published corrections

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
