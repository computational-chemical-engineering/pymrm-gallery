"""Tests for scripts/learnings_status.py."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("learnings_status", ROOT / "scripts" / "learnings_status.py")
ls = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ls)

TABLE = """| # | Date | Page | Finding | Evidence | Target | Status |
|---|---|---|---|---|---|---|
| 1 | d | `A1.1` | a `\\|x\\|` term | e | `numerics` | open |
| 2 | d | `A1.2` | f | e | `gallery` | open |
| 3 | d | `A1.3` | f | e | checklist 8 | open |
| 4 | d | `A1.4` | f | e | `style-guide` | merged (pymrm #51) |
"""


def test_rows_keep_escaped_pipes_and_split_targets():
    t = ls.rows(TABLE)
    assert t[1] == ("`numerics`", "open") and t[4][1] == "merged (pymrm #51)"
    assert [n for n in t if ls.is_open(t[n][1]) and ls.targets_pymrm(t[n][0])] == [1]


def test_whole_file_parses():
    t = ls.rows((ROOT / "docs" / "refresh-learnings.md").read_text(encoding="utf-8"))
    assert len(t) > 300 and all(s for _, s in t.values())
