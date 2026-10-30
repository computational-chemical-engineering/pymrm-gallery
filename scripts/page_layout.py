"""The layout of a migrated page, shared by check_metadata.py and release_ref.py.

A page is *migrated* when it has an ``index.py`` or when a code cell of its
notebook imports ``gallery_bootstrap`` or assigns ``GALLERY_REF = ...`` at the
start of a line. Detection is deliberately broad: a migrated page then has to
pass every check (the exact ``GALLERY_REF = `` line, the template cell, the
pairing), so a reformatted bootstrap cannot escape them or the release rewrite.
The bootstrap is ``templates/bootstrap_cell.py`` (docs/page-infrastructure.md,
section 2) and must be the first code cell. Unmigrated pages are left alone.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_CELL = ROOT / "templates" / "bootstrap_cell.py"
NAMES = ("GALLERY_REF", "PYMRM", "PAGE")
MARKER = re.compile(r"^GALLERY_REF = ", re.MULTILINE)
BOOTSTRAP_IMPORT = re.compile(
    r"^[ \t]*(?:from[ \t]+gallery_bootstrap[ \t]+import|import[ \t]+gallery_bootstrap)\b",
    re.MULTILINE)
LOCK = ROOT / "requirements-lock.txt"
LITERAL = re.compile(
    r'^(?P<name>GALLERY_REF|PYMRM|PAGE) = "(?P<value>[^"\n]*)"[ \t]*(?P<comment>#[^\n]*)?[ \t]*$',
    re.MULTILINE)


def source(cell: dict) -> str:
    """Source of a notebook cell as one string (nbformat JSON keeps a list)."""
    src = cell.get("source", "")
    return "".join(src) if isinstance(src, list) else src


def code_sources(nb: dict) -> list[str]:
    """Sources of the code cells of a notebook given as a JSON dict."""
    return [source(c) for c in nb.get("cells", []) if c.get("cell_type") == "code"]


def bootstrap_position(nb: dict) -> int | None:
    """Index among the code cells of the first one that looks like the bootstrap:
    it imports ``gallery_bootstrap`` or assigns ``GALLERY_REF``."""
    for i, src in enumerate(code_sources(nb)):
        if MARKER.search(src) or BOOTSTRAP_IMPORT.search(src):
            return i
    return None


def read_notebook(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def is_migrated(page_dir: Path) -> bool:
    page_dir = Path(page_dir)
    if (page_dir / "index.py").is_file():
        return True
    nb_path = page_dir / "index.ipynb"
    return nb_path.is_file() and bootstrap_position(read_notebook(nb_path)) is not None


def migrated_pages(root: Path = ROOT) -> list[Path]:
    """Page directories under ``root/pages`` whose notebook has the bootstrap cell."""
    pages = Path(root) / "pages"
    return [p for p in sorted(pages.iterdir()) if p.is_dir() and is_migrated(p)] \
        if pages.is_dir() else []


def literals(src: str) -> dict[str, list[str]]:
    """Every value assigned to each of the three bootstrap literals in ``src``."""
    found: dict[str, list[str]] = {n: [] for n in NAMES}
    for m in LITERAL.finditer(src):
        found[m["name"]].append(m["value"])
    return found


def normalise(src: str) -> str:
    """Cell text with the three literal values replaced by placeholders.

    Spacing before a trailing comment is normalised too, since aligning the
    comments depends on the length of the values. Trailing whitespace is dropped.
    """
    def repl(m: re.Match) -> str:
        line = f'{m["name"]} = "<{m["name"]}>"'
        return f"{line}  {m['comment']}" if m["comment"] else line

    text = LITERAL.sub(repl, src)
    return "\n".join(line.rstrip() for line in text.strip("\n").splitlines())


def template_cell() -> str:
    """The canonical bootstrap cell, with placeholders, as stored in templates/."""
    return TEMPLATE_CELL.read_text(encoding="utf-8")


def fill(template: str, gallery_ref: str, pymrm: str, page: str) -> str:
    """The template cell with its three placeholders filled."""
    return (template.replace('"<GALLERY_REF>"', f'"{gallery_ref}"')
            .replace('"<PYMRM>"', f'"{pymrm}"')
            .replace('"<PAGE>"', f'"{page}"'))


def locked_pymrm(lock: Path = LOCK) -> str | None:
    """The ``pymrm==`` version pinned in requirements-lock.txt."""
    m = re.search(r"^pymrm==(\S+)", Path(lock).read_text(encoding="utf-8"), re.MULTILINE)
    return m[1] if m else None


def pin_for(min_version: str) -> str:
    """The bootstrap's pymrm specifier for a meta.yaml ``pymrm_min_version``."""
    floor = str(min_version).strip()
    major = int(floor.split(".")[0])
    return f">={floor},<{major + 1}"
