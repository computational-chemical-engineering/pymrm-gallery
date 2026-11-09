#!/usr/bin/env python3
"""Token use of the agents a coordinating session started, per agent and per page.

    python scripts/agent_tokens.py ~/.claude-tue/projects/<project>/<session-id>
    python scripts/agent_tokens.py <session dir> --page A2.6
    python scripts/agent_tokens.py <session dir> --since 2026-11-09

Reads ``<session dir>/subagents/agent-*.jsonl`` and their ``.meta.json`` (agent type and
the description given at launch) and the session's own ``<session-id>.jsonl`` beside
the directory. Per agent: calls, prompt tokens (uncached + cache read + cache write),
cache writes, output, and the context of its last call. ``--page`` keeps the agents
whose description names the catalog ID and adds a per-page total, which is the number
the refresh cost log (docs/refresh-cost-log.md) records.

Cost is roughly calls times mean context: a file read early is paid again on every
later call, which is why the per-agent context is printed.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def usage(path: Path, since: str | None = None) -> dict:
    """Sum the usage of one transcript, counting each API response once."""
    seen, calls = set(), 0
    prompt = write = out = last = 0
    first = None
    for line in path.open(errors="ignore"):
        try:
            d = json.loads(line)
        except ValueError:
            continue
        ts = d.get("timestamp", "")
        if since and ts and ts < since:
            continue
        m = d.get("message")
        if not isinstance(m, dict) or not m.get("usage"):
            continue
        key = (m.get("id"), d.get("requestId"))
        if key in seen:
            continue
        seen.add(key)
        u = m["usage"]
        ctx = (u.get("input_tokens", 0) + u.get("cache_read_input_tokens", 0)
               + u.get("cache_creation_input_tokens", 0))
        calls += 1
        prompt += ctx
        write += u.get("cache_creation_input_tokens", 0)
        out += u.get("output_tokens", 0)
        last = ctx
        first = first or ts
    return {"calls": calls, "prompt": prompt, "write": write, "out": out,
            "last": last, "start": (first or "")[:16]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("session", type=Path, help="the session directory (holds subagents/)")
    ap.add_argument("--page", help="catalog ID, e.g. A2.6: only agents naming it")
    ap.add_argument("--since", help="ISO date or time; ignore calls before it")
    args = ap.parse_args()

    rows = []
    for f in sorted((args.session / "subagents").glob("agent-*.jsonl")):
        meta = f.with_suffix(".meta.json")
        info = json.loads(meta.read_text()) if meta.is_file() else {}
        desc = info.get("description", "")
        if args.page and not re.search(rf"(?<![\w.]){re.escape(args.page)}(?![\w.]*\d)",
                                       desc):
            continue
        u = usage(f, args.since)
        if u["calls"]:
            rows.append((u["start"], info.get("agentType", "?"), desc, u))

    print(f"{'start':16} {'agent':17} {'calls':>5} {'prompt M':>8} {'write M':>7} "
          f"{'out k':>6} {'last k':>6}  description")
    tot = {"calls": 0, "prompt": 0, "write": 0, "out": 0}
    for start, kind, desc, u in sorted(rows):
        print(f"{start:16} {kind[:17]:17} {u['calls']:5d} {u['prompt'] / 1e6:8.1f} "
              f"{u['write'] / 1e6:7.2f} {u['out'] / 1e3:6.0f} {u['last'] / 1e3:6.0f}  "
              f"{desc[:50]}")
        for k in tot:
            tot[k] += u[k]
    print(f"{'agents total':34} {tot['calls']:5d} {tot['prompt'] / 1e6:8.1f} "
          f"{tot['write'] / 1e6:7.2f} {tot['out'] / 1e3:6.0f}")
    main_log = args.session.with_suffix(".jsonl")
    if main_log.is_file() and not args.page:
        u = usage(main_log, args.since)
        print(f"{'coordinator':34} {u['calls']:5d} {u['prompt'] / 1e6:8.1f} "
              f"{u['write'] / 1e6:7.2f} {u['out'] / 1e3:6.0f} {u['last'] / 1e3:6.0f}")


if __name__ == "__main__":
    main()
