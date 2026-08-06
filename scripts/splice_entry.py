"""Splice a queue_cases/<ID>/models_entry.yaml block into models.yaml.

Usage: python scripts/splice_entry.py <CASE_ID> <BEFORE_ID>

Inserts the entry immediately before `  - id: <BEFORE_ID>`. Handles the three
shapes entries arrive in (list item indented 2, list item at column 0, bare
mapping with no dash), flips a ready/in-progress status to published, and
validates the result parses with the case present exactly once.

Hard-won caveats (2026-08-05/06):
- CHECK FOR AN EXISTING ENTRY FIRST: `grep -c "id: <CASE_ID>" models.yaml`.
  Some models_entry files say "upgrade the planned entry in place" - appending
  then creates a duplicate id, which check_metadata.py catches (C2.10 did this).
- Some entries lack a top-level `status:` altogether (B3.2 had only
  data.status) - add `status: published` and `page:` to the entry file first.
- After staging, diff the staged models.yaml: a pure splice DELETES NOTHING;
  unexpected deletions mean a concurrent agent's work was swept in.
"""
import pathlib, re, sys, yaml

def splice(case_id, before_id, slug_expected=None):
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
    mp = pathlib.Path("models.yaml")
    mt = mp.read_text().splitlines(True)
    idx = next(k for k, l in enumerate(mt) if l.rstrip("\n") == f"  - id: {before_id}")
    mp.write_text("".join(mt[:idx]) + block + "".join(mt[idx:]))
    d = yaml.safe_load(mp.read_text())
    e = d["models"]
    a = [m for m in e if m["id"] == case_id][0]
    assert a["status"] == "published", a["status"]
    if slug_expected:
        assert a["slug"] == slug_expected, a["slug"]
    print(f"{case_id}: spliced before {before_id}, {len(e)} entries, "
          f"status={a['status']}, slug={a['slug']}, page={a.get('page')}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    splice(sys.argv[1], sys.argv[2])
