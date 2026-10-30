"""Tests for the migrated-page checks of scripts/check_metadata.py
(docs/page-infrastructure.md, sections 3 and 4)."""

import shutil
import sys
from pathlib import Path

import jupytext
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import check_metadata as cm  # noqa: E402
import page_layout  # noqa: E402

PAGE = "X0.0-demo"
META = {"id": "X0.0", "pymrm_min_version": "3.1.0"}
TEMPLATE_PY = ROOT / "templates" / "index.py"
TEMPLATE_PAGE = 'PAGE = "X0.0-replace-me"'


def make_page(tmp_path, edit=lambda text: text, sync=True, formats=True):
    """A migrated fixture page: index.py from the template, index.ipynb paired."""
    page = tmp_path / PAGE
    page.mkdir()
    text = TEMPLATE_PY.read_text(encoding="utf-8").replace(TEMPLATE_PAGE, f'PAGE = "{PAGE}"')
    (page / "index.py").write_text(text, encoding="utf-8")
    nb = jupytext.read(page / "index.py")
    if not formats:
        nb.metadata.pop("jupytext", None)
    jupytext.write(nb, page / "index.ipynb")
    edited = edit(text)
    (page / "index.py").write_text(edited, encoding="utf-8")
    if sync and edited != text:
        nb = jupytext.read(page / "index.py")
        if not formats:
            nb.metadata.pop("jupytext", None)
        jupytext.write(nb, page / "index.ipynb")
    return page


def errors(page, meta=META, ref="main"):
    return cm.check_migrated(page, meta, ref)


def test_fixture_page_passes(tmp_path):
    assert errors(make_page(tmp_path)) == []


def test_released_tree_passes_with_its_tag(tmp_path):
    tag = "gallery-2026.10.30"
    page = make_page(tmp_path, lambda t: t.replace('GALLERY_REF = "main"', f'GALLERY_REF = "{tag}"'))
    assert errors(page, ref=tag) == []
    assert any("GALLERY_REF" in e for e in errors(page))


def test_aligned_comments_still_match_the_template(tmp_path):
    page = make_page(tmp_path, lambda t: t.replace(
        'PYMRM = ">=3.1.0,<4"  #', 'PYMRM = ">=3.1.0,<4"           #'))
    assert errors(page) == []


@pytest.mark.parametrize("edit, meta, needle", [
    (lambda t: t.replace(">=3.1.0,<4", ">=3.0.0,<4"), META, "PYMRM"),
    (lambda t: t, {**META, "pymrm_min_version": "3.2"}, "requires '>=3.2,<4'"),
    (lambda t: t, {"id": "X0.0"}, "no pymrm_min_version"),
    (lambda t: t.replace(f'PAGE = "{PAGE}"', 'PAGE = "X0.0-other"'), META, "directory name"),
    (lambda t: t.replace('GALLERY_REF = "main"', 'GALLERY_REF = "dev"'), META, "expected 'main'"),
    (lambda t: t.replace("setup(PAGE, PYMRM, GALLERY_REF)",
                         "setup(PAGE, PYMRM, GALLERY_REF)\nimport numpy"), META, "templates/bootstrap_cell.py"),
    (lambda t: t.replace("from gallery_bootstrap import setup",
                         'from gallery_bootstrap import setup\nPAGE = "x"'), META, "PAGE 2 times"),
    (lambda t: t.replace("# Bootstrap: the pymrm", "# %%\nx = 1\n\n# %%\n# Bootstrap: the pymrm"),
     META, "must be the first"),
])
def test_each_failure_is_reported(tmp_path, edit, meta, needle):
    errs = errors(make_page(tmp_path, edit), meta)
    assert any(needle in e for e in errs), errs


def test_index_py_out_of_step_is_reported(tmp_path):
    page = make_page(tmp_path, lambda t: t + "\n# %%\nprint('only in index.py')\n", sync=False)
    errs = errors(page)
    assert any("differ" in e for e in errs), errs


def test_missing_pairing_metadata_is_reported(tmp_path):
    errs = errors(make_page(tmp_path, formats=False))
    assert any("jupytext.formats" in e for e in errs), errs


def test_migrated_page_without_index_py_is_reported(tmp_path):
    page = make_page(tmp_path)
    (page / "index.py").unlink()
    errs = errors(page)
    assert any("without index.py" in e for e in errs), errs


def test_migrated_page_with_build_page_is_reported(tmp_path):
    page = make_page(tmp_path)
    (page / "build_page.py").write_text("# old generator\n")
    errs = errors(page)
    assert any("build_page.py" in e for e in errs), errs


def test_reformatted_bootstrap_is_still_detected_and_reported(tmp_path):
    """No `GALLERY_REF = ` line: detection by the import, then the marker is required."""
    page = make_page(tmp_path, lambda t: t.replace('GALLERY_REF = "main"', 'GALLERY_REF="main"'))
    nb = page_layout.read_notebook(page / "index.ipynb")
    assert page_layout.bootstrap_position(nb) == 0
    (page / "index.py").unlink()
    assert page_layout.is_migrated(page)                     # by the import alone
    errs = errors(page)
    assert any("GALLERY_REF = " in e for e in errs), errs


def test_index_py_alone_marks_a_page_migrated(tmp_path):
    page = tmp_path / "X0.2-py-only"
    page.mkdir()
    shutil.copy(ROOT / "templates" / "page-skeleton" / "index.ipynb", page / "index.ipynb")
    assert not page_layout.is_migrated(page)
    (page / "index.py").write_text("# %%\nx = 1\n")
    assert page_layout.is_migrated(page)
    errs = errors(page)
    assert any("no bootstrap cell" in e for e in errs), errs


def test_pin_the_lock_does_not_satisfy_is_reported(tmp_path):
    """meta.yaml says 2.3 and the cell matches it, but the lock has pymrm 3.1.0."""
    page = make_page(tmp_path, lambda t: t.replace(">=3.1.0,<4", ">=2.3,<3"))
    meta = {**META, "pymrm_min_version": "2.3"}
    errs = errors(page, meta)
    assert len(errs) == 1 and "does not satisfy" in errs[0], errs
    assert cm.check_migrated(page, meta, "main", lock_version="2.5.0") == []


def test_lock_version_is_read():
    assert page_layout.locked_pymrm() is not None


def test_migrated_pages_detects_only_the_bootstrap(tmp_path):
    (tmp_path / "pages").mkdir()
    shutil.move(make_page(tmp_path), tmp_path / "pages" / PAGE)
    old = tmp_path / "pages" / "X0.1-old"
    old.mkdir()
    shutil.copy(ROOT / "templates" / "page-skeleton" / "index.ipynb", old / "index.ipynb")
    assert page_layout.migrated_pages(tmp_path) == [tmp_path / "pages" / PAGE]


def test_templates_agree_with_each_other():
    """templates/index.py carries the canonical cell and is paired with index.ipynb."""
    nb = page_layout.read_notebook(ROOT / "templates" / "index.ipynb")
    assert page_layout.bootstrap_position(nb) == 0
    cell = page_layout.code_sources(nb)[0]
    assert page_layout.normalise(cell) == page_layout.normalise(page_layout.template_cell())
    assert cm.check_paired(ROOT / "templates") == []


def test_pin_for():
    assert page_layout.pin_for("3.1.0") == ">=3.1.0,<4"
    assert page_layout.pin_for("2.3") == ">=2.3,<3"
