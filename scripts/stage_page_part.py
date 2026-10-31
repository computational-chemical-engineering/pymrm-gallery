#!/usr/bin/env python3
"""Stage one page's part of the two files that every page refresh edits.

    python scripts/stage_page_part.py A2.8 70,71,72,73

When two pages are refreshed at once, models.yaml and docs/refresh-learnings.md
hold changes of both, and an agent has no interactive ``git add -p``. This stages
HEAD's text with only the named page's models.yaml entry swapped in and only the
named learnings rows appended, and leaves the working tree as it is.
"""

import re
import subprocess
import sys


def head(path):
    return subprocess.run(["git", "show", f"HEAD:{path}"], capture_output=True,
                          text=True, check=True).stdout


def stage(path, text):
    sha = subprocess.run(["git", "hash-object", "-w", "--stdin"], input=text,
                         capture_output=True, text=True, check=True).stdout.strip()
    subprocess.run(["git", "update-index", "--cacheinfo", f"100644,{sha},{path}"],
                   check=True)


def entry(text, page_id):
    return re.search(rf"^  - id: {re.escape(page_id)}\n.*?(?=^  - id: |\Z)", text,
                     re.S | re.M)


def main():
    page_id, rows = sys.argv[1], sys.argv[2].split(",")
    old, new = head("models.yaml"), open("models.yaml").read()
    mo, mn = entry(old, page_id), entry(new, page_id)
    stage("models.yaml", old[:mo.start()] + mn.group(0) + old[mo.end():])

    row = re.compile(r"\| (\d+) \|")
    lines = [ln for ln in open("docs/refresh-learnings.md").read().splitlines()
             if row.match(ln) and row.match(ln).group(1) in rows]
    if len(lines) != len(rows):
        sys.exit(f"found {len(lines)} of {len(rows)} learnings rows")
    stage("docs/refresh-learnings.md",
          head("docs/refresh-learnings.md").rstrip("\n") + "\n" + "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
