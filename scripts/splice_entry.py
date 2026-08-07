"""Splice a queue_cases/<ID>/models_entry.yaml block into models.yaml.

Usage: python scripts/splice_entry.py <CASE_ID> <BEFORE_ID>   # new entry
       python scripts/splice_entry.py --replace <CASE_ID>     # upgrade in place

Inserts the entry immediately before `  - id: <BEFORE_ID>`. Handles the three
shapes entries arrive in (list item indented 2, list item at column 0, bare
mapping with no dash), flips a ready/in-progress status to published, and
validates the result parses with the case present exactly once.

`--replace` is for a case that ALREADY has a `planned` block in models.yaml:
it swaps that block for the built one instead of appending beside it. Use it
when integrate_case.py stops with "upgrade-in-place needed", then run the rest
of that script's steps (status flips, copy to pages/, gates) yourself.

Hard-won caveats (2026-08-05/06):
- An existing entry is now detected rather than trusted to you: splice() refuses
  to append over one and points at --replace. Appending anyway creates a
  duplicate id, which check_metadata.py catches (C2.10 did this); E1.1 was the
  second case to need the in-place path, which is why it is a mode now.
- Some entries lack a top-level `status:` altogether (B3.2 had only
  data.status) - add `status: published` and `page:` to the entry file first.
- After staging, diff the staged models.yaml: a pure splice DELETES NOTHING;
  unexpected deletions mean a concurrent agent's work was swept in.
"""
import pathlib, re, sys, yaml

def entry_block(case_id):
    """The case's models_entry.yaml, normalised to a `  - id:` list item."""
    ep = pathlib.Path(f"queue_cases/{case_id}/models_entry.yaml")
    raw = ep.read_text().splitlines(True)
    start = next(k for k, l in enumerate(raw)
                 if re.match(rf"^\s*-?\s*id:\s*{re.escape(case_id)}\s*$", l))
    body = raw[start:]
    first = body[0]
    ind = len(first) - len(first.lstrip())
    has_dash = first.lstrip().startswith("- ")
    if has_dash and ind == 2:
        block = "".join(body)
    elif has_dash and ind == 0:
        block = "".join(("  " + l if l.strip() else l) for l in body)
    else:  # bare mapping, no dash
        out = ["  - " + first.lstrip()]
        for l in body[1:]:
            out.append(("    " + l[ind:]) if l.strip() else l)
        block = "".join(out)
    block = re.sub(r"^(\s*)status:\s*(ready|in-progress).*$",
                   r"\1status: published", block, count=1, flags=re.M)
    if not block.endswith("\n"):
        block += "\n"
    block += "\n"
    return block

def _check(case_id, slug_expected, what):
    """models.yaml must parse with the case present exactly once and published."""
    d = yaml.safe_load(pathlib.Path("models.yaml").read_text())
    e = d["models"]
    hits = [m for m in e if m["id"] == case_id]
    if len(hits) != 1:
        raise SystemExit(f"{case_id}: {len(hits)} entries after {what} — expected 1")
    a = hits[0]
    assert a["status"] == "published", a["status"]
    if slug_expected:
        assert a["slug"] == slug_expected, a["slug"]
    print(f"{case_id}: {what}, {len(e)} entries, "
          f"status={a['status']}, slug={a['slug']}, page={a.get('page')}")

def splice(case_id, before_id, slug_expected=None):
    """Insert the entry immediately before `  - id: <before_id>`."""
    block = entry_block(case_id)
    mp = pathlib.Path("models.yaml")
    mt = mp.read_text().splitlines(True)
    if any(l.rstrip("\n") == f"  - id: {case_id}" for l in mt):
        raise SystemExit(f"{case_id} already in models.yaml — use --replace "
                         "(upgrade in place), not a splice")
    idx = next(k for k, l in enumerate(mt) if l.rstrip("\n") == f"  - id: {before_id}")
    mp.write_text("".join(mt[:idx]) + block + "".join(mt[idx:]))
    _check(case_id, slug_expected, f"spliced before {before_id}")

def replace(case_id, slug_expected=None):
    """Upgrade a `planned` entry IN PLACE: swap the existing block for the built one.

    Kept separate from splice() on purpose — appending over an existing id is
    the C2.10 duplicate-id incident, and integrate_case.py refuses that case and
    sends you here. The block runs from its `  - id:` line to the next one at
    the same indent (or the end of the list).
    """
    block = entry_block(case_id)
    mp = pathlib.Path("models.yaml")
    mt = mp.read_text().splitlines(True)
    starts = [k for k, l in enumerate(mt) if l.rstrip("\n") == f"  - id: {case_id}"]
    if len(starts) != 1:
        raise SystemExit(f"{case_id}: {len(starts)} existing blocks — expected exactly 1")
    start = starts[0]
    end = next((k for k in range(start + 1, len(mt))
                if re.match(r"^  - id: ", mt[k])), None)
    if end is None:  # last entry in the list: stop at the first non-indented line
        end = next((k for k in range(start + 1, len(mt))
                    if mt[k].strip() and not mt[k].startswith("    ")), len(mt))
    old_status = next((m.group(1) for m in
                       (re.match(r"^    status:\s*(\S+)", l) for l in mt[start:end]) if m),
                      "?")
    print(f"{case_id}: replacing lines {start+1}-{end} "
          f"(old status={old_status}, {end-start} lines) "
          f"with {len(block.splitlines())} lines")
    mp.write_text("".join(mt[:start]) + block + "".join(mt[end:]))
    _check(case_id, slug_expected, "upgraded in place")

if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 2 and a[0] == "--replace":
        replace(a[1])
    elif len(a) == 2:
        splice(a[0], a[1])
    else:
        sys.exit(__doc__)
