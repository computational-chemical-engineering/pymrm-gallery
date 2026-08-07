# Prompt for the next session

Copy everything in the block below as the first message of a fresh session.

---

```
Continue the pymrm-gallery. Repo ~/Code/pymrm_suite/pymrm-gallery, venv ~/Code/pymrm_suite/.venv.

Operate per docs/playbooks/coordinator.md: read docs/handoff.md header + git log -5,
then work ONE CASE AT A TIME — builder → verifier → fix → integrate → commit — and
only start the next once that one is pushed. Do NOT run builders in parallel; that
is what lost a whole wave to a usage limit on 2026-08-07. Dispatch agents with the
10-line template in docs/playbooks/README.md (playbooks carry the rules — do not
restate them) and set model: per its tier table (verifiers never downgraded).
Per-source traps are in docs/papers-inventory.yaml; integrate with
scripts/integrate_case.py; pull --rebase before every push.

Next cases, in order: A3.5, E1.1, A4.1, A2.4, J1.3 — each has an UNVERIFIED PARTIAL
BUILD on disk from the 2026-08-07 wave, whose builders were all killed by an API
session limit. Read each case's resume: block first. Do NOT re-execute the staged
build_page.py scripts and publish the result: A2.4's and J1.3's notebooks are older
than their own scripts, and A2.4's builder had found an eq. (4) scaling bug it never
finished fixing. Treat every staged digit as unread. Verifier on every ready page.
Park blocked cases with a resume: block — never wait for me.

Republish dashboards at wind-down to the SAME URLs (never mint new ones):
  papers    https://claude.ai/code/artifact/99b52225-ce54-4487-99ec-3a420f2ac4ad  (favicon 📄)
  decisions https://claude.ai/code/artifact/237909a0-3432-4b82-8cfe-15b05019d036  (favicon ❓)

Halt = coordinator.md's checklist: integrate what's ready, update handoff,
republish dashboards, gates green, push, end.
```

---

## Why each part is there

**One wave per session** is the token-discipline rule from
[`playbooks/coordinator.md`](playbooks/coordinator.md). All state lives in the repo,
so ending a session costs nothing and a second wave in the same context costs a lot.

**"Never block on me"** is the maintainer's standing correction. A case needing input
is parked — with a `resume:` block recording what is established, what each plausible
answer changes, and which files to touch — and the next case starts at once. Assume no
answers are coming; the maintainer reads the dashboards when they choose.

**Dispatch policy** is the main throughput lever, and it is about avoiding the figure
review gate rather than servicing it faster. Ranked: a worked example with printed
intermediates beats an internal identity beats a stated numerical result beats a
digitised figure. Only the last needs a human.

**The verifier** stays on the maintainer's explicit instruction. It roughly doubles the
cost per page and is aimed at the one failure this repository has actually suffered:
plausible, confident agreement that turned out to be circular.

**The favicons are recorded** because a changed favicon reads as a different page — the
maintainer finds these tabs by their icon. They were unrecorded until 2026-08-07 and had
to be chosen; keep them stable from here.

## State at the time of writing (2026-08-07, after the failed wave)

| | |
|---|---|
| Published pages | 58 directories, 59 catalogue entries, 68 models.yaml entries |
| Gates | `check_agreement` 58 pages / 0 regressions; `check_metadata` OK, 0 warnings |
| Unclaimed (dispatchable) | 58 |
| Needing a paper | 140 |
| Covered / deferred | 6 / 3 |
| Carrying unverified partial builds | `A3.5`, `E1.1`, `A4.1`, `A2.4`, `J1.3` |

Counts differ between page directories and catalogue entries because `B1.1` covers
`B1.5`, `A2.1` covers `A2.2`, and `C2.10` covers `D3.4`.

**The 2026-08-07 wave published nothing.** All five builders were killed mid-task by the
same account-level API session limit, not by anything about the cases. The full
post-mortem — what each left on disk, why neither staged notebook can be trusted, and
why the five stay `unclaimed` rather than `in-progress` — is in
[`handoff.md`](handoff.md) under Recommended next moves. The lesson worth carrying:
session-limit deaths cannot be resumed via `SendMessage` and need fresh dispatches,
whereas connection-error deaths can be resumed with context intact.

## Two decisions still carried

Neither blocks work — ask once, record the answer when it comes, continue.

1. Purging three copyrighted overlays from git history (needs a force-push to a public
   repo).
2. The Elsevier API key at `~/.config/elsevier/apikey` is blocked by the Bash permission
   classifier — the highest-leverage unblock for section A.
