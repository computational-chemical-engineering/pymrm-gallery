# Coordinator playbook — the wave loop, and how to halt/resume cheaply

The session-level protocol. Everything an agent needs is in the repo; the
coordinator's job is dispatch, integration, and honest relay — not re-teaching.

## Token discipline (why this file exists)

Measured on the 2026-08-02/06 run: each page cost ~0.7–1.2M subagent tokens,
dominated by (a) 800-word bespoke dispatch prompts restating standing rules —
now replaced by the playbooks + the 10-line template in `README.md` — and
(b) coordinator context growing monotonically across multi-day sessions. The
fixes:

- **One session = one wave.** A wave: pick 4–6 unclaimed cases → dispatch
  builders (playbook prompts) → verifier per ready page → fix (sonnet for
  wording-only; strong model when numbers change; NEVER downgrade verifiers)
  → integrate → wind down → END THE SESSION. Do not carry a second wave in the
  same context; starting fresh costs nothing because all state is on disk.
- **Dispatch via the template** in `playbooks/README.md` + per-case notes from
  `papers-inventory.yaml` and the case yaml. Do not restate the playbooks.
- **Set `model:` per the tier table** in `playbooks/README.md`.
- **Cap concurrent agents at ~5**, one case each. Park blocked cases with a
  `resume:` block; never wait on the maintainer.
- **Relay reports compactly**; do not re-quote whole reports back to agents.

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
3. Pick the next wave from `unclaimed` by tier/priority; consult
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
