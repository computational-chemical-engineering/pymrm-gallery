#!/usr/bin/env python3
"""Collect the facts a refresh brief needs about one page, from the old version.

    python scripts/page_facts.py A3.4 > <scratch>/a34_facts.md
    python scripts/page_facts.py A3.4 --base origin/main

Run it with the gallery venv (it needs PyYAML).

The coordinator writes this file once and the refresh brief points to it, so neither
the coordinator nor the refresh agent reads whole files to find the solver calls, the
published metrics or the relevant learnings rows. Printed, as Markdown:

- the models.yaml entry's title, structure codes and lane hint (model-solving when
  build_page.py calls a pymrm operator or solver, light otherwise);
- the old build_page.py: length, every def and class, and the lines with solves, root
  finds, extrema, interpolation, outlet reads and data loads, with line numbers;
- data the old page borrows from another page (``load_data(..., page=...)``), which
  brings the AGENTS.md data rules into the brief;
- the published metrics of the old agreement.json;
- refreshed pages that share a structure code (candidates for the one worked
  example), with their lane;
- the learnings rows that name this page or a page sharing a structure code, as row
  number, pages and the first words of the finding.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SOLVER = re.compile(r"newton\(|solve_ivp|spsolve|splu|brentq|fsolve|root_scalar|"
                    r"minimize|find_roots|argmax|argmin|np\.interp|\[-1\]|"
                    r"construct_|NumJac|ConvFlux|BlockLayout|load_data|read_csv|"
                    r"report_agreement|curve_fit|least_squares")
MODEL_SOLVING = re.compile(r"construct_(grad|div|convflux)|NumJac|newton\(|ConvFlux")


def show(path: str, base: str) -> str:
    proc = subprocess.run(["git", "show", f"{base}:{path}"], cwd=ROOT,
                          capture_output=True, text=True)
    return proc.stdout if proc.returncode == 0 else ""


def entries(text: str) -> dict[str, dict]:
    return {str(e["id"]): e for e in (yaml.safe_load(text) or {}).get("models", [])}


def field(entry: dict, name: str) -> str:
    v = entry.get(name, "")
    return ", ".join(map(str, v)) if isinstance(v, list) else str(v or "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("id", help="catalog ID, e.g. A3.4")
    ap.add_argument("--base", default="origin/main")
    args = ap.parse_args()

    models = entries(show("models.yaml", args.base))
    entry = models.get(args.id)
    if not entry:
        raise SystemExit(f"{args.id}: no models.yaml entry on {args.base}")
    page = field(entry, "page").rstrip("/")
    structs = set(re.findall(r"S\d+", field(entry, "structures")))
    build = show(f"{page}/build_page.py", args.base)
    lane = "model-solving" if MODEL_SOLVING.search(build) else "light"

    print(f"# Page facts: {args.id} ({page}, on {args.base})\n")
    print(f"- title: {field(entry, 'title')}")
    print(f"- structures: {', '.join(sorted(structs)) or 'none'}")
    print(f"- lane hint: {lane} (from the calls in build_page.py; the coordinator decides)")

    print(f"\n## Old build_page.py ({len(build.splitlines())} lines)\n")
    if not build:
        print("No build_page.py on the base: the page is already refreshed.")
    defs, calls = [], []
    for i, ln in enumerate(build.splitlines(), 1):
        if re.match(r"\s*(def|class) ", ln):
            defs.append(f"{i}: {ln.strip()[:100]}")
        elif SOLVER.search(ln) and not ln.lstrip().startswith("#"):
            calls.append(f"{i}: {ln.strip()[:110]}")
    if defs:
        print("Definitions:\n```\n" + "\n".join(defs) + "\n```")
    if calls:
        print("Solves, root finds, extrema, outlet reads, data loads:\n```\n"
              + "\n".join(calls) + "\n```")

    borrowed = sorted(set(re.findall(r"load_data\([^)]*page\s*=\s*[\"']([^\"']+)", build))
                      - {page.rsplit("/", 1)[-1], args.id})
    print("\n## Borrowed data\n")
    print("From " + ", ".join(borrowed) + ": read AGENTS.md \"Data rules\" and the source "
          "page's findings on the rows used." if borrowed else "None.")

    agr = show(f"{page}/agreement.json", args.base)
    print("\n## Published metrics (old agreement.json)\n")
    try:
        data = json.loads(agr)
        if data.get("format") == 2 or "claims" in data:
            print("Already format 2.")
        for k, v in (data.get("metrics") or {}).items():
            print(f"- `{k}` = {v!r}")
    except ValueError:
        print("No agreement.json.")

    print("\n## Refreshed pages sharing a structure code\n")
    shared = []
    for cid, e in models.items():
        p = field(e, "page").rstrip("/")
        their = set(re.findall(r"S\d+", field(e, "structures")))
        if cid == args.id or not p or not structs & their:
            continue
        if show(f"{p}/index.py", args.base):
            has_model = bool(show(f"{p}/model.py", args.base))
            shared.append((cid, p, has_model))
            print(f"- {cid} `{p}` ({'model.py' if has_model else 'no model.py'}), "
                  f"structures {field(e, 'structures')}")
    if not shared:
        print("None.")

    print("\n## Learnings rows\n")
    ids = {cid for cid, _, _ in shared}
    learn = show("docs/refresh-learnings.md", args.base)
    own, near = [], []
    for ln in learn.splitlines():
        cells = [c.strip() for c in ln.split("|")]
        if len(cells) < 6 or not cells[1].isdigit():
            continue
        named = set(re.findall(r"`([A-J]\d+\.\d+)`", cells[3]))
        if args.id in named:
            own.append(f"- row {cells[1]} ({cells[3]}): {cells[4]}")
        elif named & ids:
            near.append(cells[1])
    print("Naming this page:" if own else "None name this page.")
    print("\n".join(own))
    if near:
        print(f"\nNaming a page that shares a structure code ({len(near)}; read one with "
              f"`grep -n '^| <row> |' docs/refresh-learnings.md` when its rule applies): "
              + ", ".join(near))


if __name__ == "__main__":
    main()
