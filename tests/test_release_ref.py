"""Tests for scripts/release_ref.py (docs/page-infrastructure.md, section 1)."""

import sys
from pathlib import Path

import nbformat
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import page_layout  # noqa: E402
import release_ref  # noqa: E402

TAG = "gallery-2026.10.30"


def cell_text(gallery_ref="main"):
    return page_layout.fill(page_layout.template_cell(), gallery_ref, ">=3.1.0,<4", "X0.0-demo")


def write_nb(path, *code):
    nb = nbformat.v4.new_notebook()
    nb.cells = [nbformat.v4.new_markdown_cell('GALLERY_REF = "main" in prose is not code')]
    nb.cells += [nbformat.v4.new_code_cell(c) for c in code]
    nbformat.write(nb, path)


def test_notebook_replaced_once(tmp_path):
    path = tmp_path / "index.ipynb"
    write_nb(path, cell_text(), "x = 1")
    release_ref.rewrite_notebook(path, TAG)
    nb = nbformat.read(path, as_version=4)
    assert page_layout.literals(nb.cells[1].source)["GALLERY_REF"] == [TAG]
    assert "# replaced by the release tag" in nb.cells[1].source
    assert nb.cells[0].source.startswith('GALLERY_REF = "main"')     # markdown untouched


@pytest.mark.parametrize("cells", [
    ("x = 1",),                                   # zero
    (cell_text(), 'GALLERY_REF = "main"'),         # two, in two cells
    (cell_text() + '\nGALLERY_REF = "main"\n',),   # two, in one cell
])
def test_notebook_fails_unless_exactly_one(tmp_path, cells):
    path = tmp_path / "index.ipynb"
    write_nb(path, *cells)
    before = path.read_text()
    with pytest.raises(ValueError, match="exactly one"):
        release_ref.rewrite_notebook(path, TAG)
    assert path.read_text() == before


def test_py_replaced_once(tmp_path):
    path = tmp_path / "index.py"
    path.write_text("# %%\n" + cell_text() + "\n# %%\nx = 1\n")
    release_ref.rewrite_py(path, TAG)
    assert page_layout.literals(path.read_text())["GALLERY_REF"] == [TAG]


@pytest.mark.parametrize("text", [
    "# %%\nx = 1\n",
    "# %%\n" + cell_text() + '\n# %%\nGALLERY_REF = "main"\n',
    "# %%\n" + cell_text(TAG) + "\n",            # already released: zero "main"
])
def test_py_fails_unless_exactly_one(tmp_path, text):
    path = tmp_path / "index.py"
    path.write_text(text)
    with pytest.raises(ValueError, match="exactly one"):
        release_ref.rewrite_py(path, TAG)
    assert path.read_text() == text
