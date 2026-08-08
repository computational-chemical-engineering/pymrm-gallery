# Coordinator playbook — the one-case loop, and how to halt/resume cheaply

The session-level protocol. Everything an agent needs is in the repo; the
coordinator's job is dispatch, integration, and honest relay — not re-teaching.

## Token discipline (why this file exists)

Measured on the 2026-08-02/06 run: each page cost ~0.7–1.2M subagent tokens,
dominated by (a) 800-word bespoke dispatch prompts restating standing rules —
now replaced by the playbooks + the 10-line template in `README.md` — and
(b) coordinator context growing monotonically across multi-day sessions. The
fixes:

- **ONE CASE AT A TIME, START TO FINISH.** Maintainer's instruction, 2026-08-07,
  and it overrides the older "~5 concurrent agents" rule everywhere it appears.
  Take one case all the way — builder → verifier → fix → integrate → commit —
  and only then start the next. Do not dispatch a second builder while the
  first case is unfinished.
  **Why, precisely:** the lean playbooks cut cost *per dispatch*, but the API
  limit is on usage *rate*, so N parallel builders burn N× as fast and a limit
  hit destroys N cases instead of one. On 2026-08-07 five concurrent builders
  were killed mid-task simultaneously; the wave published nothing and left five
  cases holding unverified half-built scripts, two with notebooks staler than
  their own build scripts. Serial work does not make a page cheaper — it caps
  the blast radius at one case, and a case finished and committed can never be
  lost to a limit.
- **Integrate and commit each case before starting the next.** Committed work
  survives a limit; staged work on disk survives only if its `resume:` block
  says honestly how far to trust it.
- **One session = as many complete cases as the budget allows.** End the
  session when a case is finished and the next would start; starting fresh
  costs nothing because all state is on disk. Never leave a case half-built to
  begin another.
- **Dispatch via the template** in `playbooks/README.md` + per-case notes from
  `papers-inventory.yaml` and the case yaml. Do not restate the playbooks.
- **Set `model:` per the tier table** in `playbooks/README.md`.
- Park blocked cases with a `resume:` block; never wait on the maintainer.
- **Relay reports compactly**; do not re-quote whole reports back to agents.

## Re-verify a fix that adds load-bearing material (2026-08-08)

The loop is builder → verifier → fixer → integrate, and a fixer that merely
applies a list needs no second pass. But when the fix **replaces or adds
something the page's conclusions rest on**, dispatch a SECOND verifier scoped
to the new material only — tell it the first pass is settled and name what
changed. Found on `H1.8`: the verifier killed the page's second route as a
tautology (Robeson's Table 13 permeability column *is* `k·α^n`, so its ratio
at fixed α is identically `10^(Δlog k + Δn·log α)` — the "independent" residual
was the difference of two identity residuals). The fixer built a genuine
replacement route, and that new route arrived carrying **the same defect class
the first pass had just caught elsewhere on the page**: its headline residuals
were themselves reference-dependent and quoted bare on three of five surfaces,
and its ranking claim failed at two of seven evaluation points. A defect class
does not stay fixed where it was found — a fix written against one instance
reproduces it in the new material. Scope the second pass tightly (it cost
~40 % of a full pass) and have it re-derive only what moved.

Two corollaries from the same case. **The fixer is allowed to overrule the
verifier and did, upward** — independently re-measuring 26 metrics with false
coverage attributions where the verifier said 22, and the second pass settled
it in the fixer's favour to the digit. And **the page's own printed output can
contradict its coverage map**: `H1.8`'s map claimed a break row that the page's
own break table two cells above showed to be a non-mover. Generate the map from
measured movers (the `A4.1` mechanism) — an assert over key *sets* cannot see it.

## Integration (per ready+verified+fixed case)

`python scripts/integrate_case.py <ID> <BEFORE_ID>` does: existing-entry
check, splice (via splice_entry), status flips (page meta, case yaml with
page: path = slug), copy to `pages/<ID>-<slug>/`, then check_metadata,
run_pages --changed, check_agreement. Review its output; on the (expected)
"NEW ... no committed baseline" line, proceed. Then:
`git pull --rebase` (stash in-flight case yamls around it — the maintainer
edits layout on GitHub) → commit named paths only (NEVER `git add -A`; check
the staged diff for deletions — a pure splice deletes nothing) → push.

## Halt (any time, cheaply)

1. Let running agents finish (builders write to disk; killing verifiers loses
   their reading). If an agent dies on a connection error, RESUME it via
   SendMessage to its id — context survives. Session-limit deaths need fresh
   dispatches.
2. Integrate what is ready; leave `ready`-but-unverified cases staged — the
   next session dispatches their verifiers from the queue state alone.
3. Update `docs/handoff.md` header (counts, gates) and next-moves; regenerate
   `scripts/dashboards.py` (0 base64 payloads) and republish BOTH dashboards
   to the maintainer's existing artifact URLs (pinned in handoff /
   papers-requested; never mint new URLs).
4. Final gates, pull, push, tree clean. Done — nothing lives only in context.

## Resume (fresh session, minimal context)

1. Read `docs/handoff.md` (header + next moves) and `git log --oneline -10`.
2. `python -c "..."` over `queue_cases/*.yaml` for status counts; anything
   `ready` gets a verifier first (quality over throughput — verifier on every
   ready page, no exceptions).
3. Pick the NEXT SINGLE CASE from `unclaimed` by tier/priority; consult
   `papers-inventory.yaml` for per-source traps; dispatch via the template.
4. Auto-memory carries the cross-session lessons; the playbooks carry the
   agent-facing ones. If a lesson is new, add it to the playbook (one place),
   not to future prompts.

## Standing cautions

- models.yaml is contested state: check for an existing entry before splicing
  (duplicate-id check is the backstop); never commit files a live agent owns.
- `.ipynb` raw bytes never reproduce (cell IDs, stdout chunking) — compare
  content; revert pure-metadata churn on already-committed pages.
- The maintainer's two dashboard URLs and the decisions/papers pages are the
  interface to the human: keep them current at every halt, and put decision
  requests THERE with the artifact/PDF named (never make them search).
