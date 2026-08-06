# Playbooks — the standing instructions for gallery agents

Why this directory exists: through 2026-08-06 every agent dispatch carried an
800-word bespoke prompt restating the same rules, and every builder re-derived
the same scaffolding. That is the single largest avoidable token cost in the
workflow. The rules now live HERE, once; a dispatch prompt is ~10 lines that
points at the right playbook and adds only what is case-specific.

## Dispatch template (coordinator: copy, fill, send)

```
You are a <ROLE> agent for the pymrm-gallery.
Repo /home/eajfpeters/Code/pymrm_suite/pymrm-gallery, venv ../.venv.
READ FIRST: docs/playbooks/<role>.md — it is your full standing instruction set.
CASE: <ID> — <one-line title>
SOURCE: <file in ~/papers/pymrm-gallery/, systematic name> (native <N> ppi)
CASE-SPECIFIC NOTES:
- <trap or scope note from the inventory / queue file, 1-3 bullets>
Report back per the playbook's report format.
```

## Model tiering (coordinator: set `model:` on the Agent call)

| Role | Model | Why |
|---|---|---|
| Builder | session default (strong) | creative modelling + transcription judgement; errors here cost a verify+fix cycle |
| Verifier | session default (strong) | the quality guarantee; independent re-derivation is the point |
| Fixer, numbers change | session default (strong) | touching computed results |
| Fixer, wording/metadata only | `sonnet` | mechanical edits against an explicit list |
| Mapper / library curation | `sonnet` | read-title-page + bookkeeping; the checklist carries the judgement |
| Dashboard/doc chores | `sonnet` | mechanical |

Never downgrade the verifier. If in doubt about a fixer, look at whether any
`agreement.json` value will change: yes → strong, no → sonnet.

## The other half of the savings

- `scripts/integrate_case.py` — one command replaces the coordinator's manual
  splice/flip/copy/gates sequence per case.
- `scripts/splice_entry.py` — shape-aware models.yaml splice (read its docstring).
- Builders copy the closest published page (AGENTS.md rule) — that IS the code
  template; the playbook names current best exemplars per structure.
- Session hygiene: one session = one wave (see coordinator.md). All state lives
  in the repo + auto-memory, so ending a session costs nothing.
