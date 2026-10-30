"""Tests for shared/gallery_bootstrap.py and offline mode in gallery_utils
(docs/page-infrastructure.md, section 2)."""

import ast
import subprocess
import sys
import types
import urllib.error
import venv
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "shared"))
import gallery_bootstrap as gb  # noqa: E402
import gallery_utils as gu  # noqa: E402

PAGE = next(p.name for p in sorted((ROOT / "pages").iterdir()) if p.is_dir())


# -- standard library only ----------------------------------------------------

def test_bootstrap_imports_only_the_standard_library():
    tree = ast.parse((ROOT / "shared" / "gallery_bootstrap.py").read_text(encoding="utf-8"))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            assert node.level == 0, "no relative imports"
            names.add(node.module.split(".")[0])
    foreign = names - set(sys.stdlib_module_names) - {"__future__"}
    assert not foreign, f"gallery_bootstrap imports non-stdlib modules: {foreign}"


def test_setup_runs_in_a_fresh_venv_without_pymrm(tmp_path):
    """A bare interpreter (no pip, no pymrm, no numpy) imports the bootstrap, and
    offline mode refuses to install pymrm instead of fetching it."""
    env_dir = tmp_path / "venv"
    venv.create(env_dir, with_pip=False)
    py = env_dir / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    code = (f"import sys; sys.path.insert(0, {str(ROOT / 'shared')!r}); "
            "from gallery_bootstrap import setup; "
            f"setup({PAGE!r}, '>=3.1.0,<4', 'main')")
    run = subprocess.run([str(py), "-c", code], capture_output=True, text=True,
                         env={"GALLERY_OFFLINE": "1", "PATH": ""}, cwd=tmp_path)
    assert run.returncode != 0
    assert "RuntimeError" in run.stderr and "GALLERY_OFFLINE" in run.stderr, run.stderr


# -- versions -------------------------------------------------------------------

@pytest.mark.parametrize("a, b", [
    ("3.1.0", "3.2.0rc1"), ("3.1.0rc1", "3.1.0"), ("3.1.0a1", "3.1.0b1"),
    ("3.1.0b2", "3.1.0rc1"), ("3.1.0.dev1", "3.1.0a1"), ("3.1.0", "3.1.0.post1"),
    ("2.9.9", "3.0"), ("3.10", "3.11.0"), ("3.9", "3.10"),
])
def test_version_order(a, b):
    assert gb.parse_version(a) < gb.parse_version(b)


def test_trailing_zeros_are_equal():
    assert gb.parse_version("3.1") == gb.parse_version("3.1.0")


@pytest.mark.parametrize("version, ok", [
    ("3.1.0", True), ("3.1", True), ("3.2.0rc1", True), ("3.9.9", True),
    ("3.1.0.post1", True), ("3.1.0rc1", False), ("3.1.0.dev1", False),
    ("3.0.2", False), ("4.0.0", False), ("4.0.0rc1", False),
])
def test_satisfies_admits_prereleases(version, ok):
    assert gb.satisfies(version, ">=3.1.0,<4") is ok


def test_unsupported_clause_raises():
    with pytest.raises(ValueError):
        gb.satisfies("3.1.0", "~=3.1")


# -- pymrm --------------------------------------------------------------------

def _no_network(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("network access attempted")
    monkeypatch.setattr(gb.urllib.request, "urlretrieve", refuse)
    monkeypatch.setattr(gb.subprocess, "check_call", refuse)


def test_newer_pymrm_warns_and_is_kept(monkeypatch):
    _no_network(monkeypatch)
    monkeypatch.setattr(gb, "_installed", lambda dist="pymrm": "4.0.1")
    with pytest.warns(UserWarning, match="newer"):
        gb._ensure_pymrm(">=3.1.0,<4")


def test_older_loaded_pymrm_asks_for_restart(monkeypatch):
    _no_network(monkeypatch)
    monkeypatch.setattr(gb, "_installed", lambda dist="pymrm": "3.0.2")
    monkeypatch.setitem(sys.modules, "pymrm", types.ModuleType("pymrm"))
    with pytest.raises(RuntimeError, match="Restart"):
        gb._ensure_pymrm(">=3.1.0,<4")


def test_older_unloaded_pymrm_is_reinstalled(monkeypatch):
    calls = []
    monkeypatch.delenv("GALLERY_OFFLINE", raising=False)
    monkeypatch.setattr(gb, "_installed", lambda dist="pymrm": "3.0.2")
    monkeypatch.delitem(sys.modules, "pymrm", raising=False)
    monkeypatch.setattr(gb.subprocess, "check_call", lambda cmd: calls.append(cmd))
    gb._ensure_pymrm(">=3.1.0,<4")
    assert calls and calls[0][-1] == "pymrm>=3.1.0,<4"


def test_older_pymrm_in_a_checkout_is_never_reinstalled(monkeypatch):
    _no_network(monkeypatch)
    monkeypatch.setattr(gb, "_installed", lambda dist="pymrm": "3.0.2")
    monkeypatch.delitem(sys.modules, "pymrm", raising=False)
    with pytest.raises(RuntimeError, match="checkout"):
        gb._ensure_pymrm(">=3.1.0,<4", checkout=True)
    with pytest.raises(RuntimeError, match="checkout"):
        gb.setup(PAGE, ">=3.1.0,<4", "main")


def test_cache_dir_falls_back_to_home(monkeypatch, tmp_path):
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    assert gb.cache_dir("main") == tmp_path / ".cache" / "pymrm-gallery" / "main"


def test_missing_pymrm_offline_raises(monkeypatch):
    _no_network(monkeypatch)
    monkeypatch.setenv("GALLERY_OFFLINE", "1")
    monkeypatch.setattr(gb, "_installed", lambda dist="pymrm": None)
    with pytest.raises(RuntimeError, match="GALLERY_OFFLINE"):
        gb._ensure_pymrm(">=3.1.0,<4")


# -- files ----------------------------------------------------------------------

@pytest.fixture
def clean_state(monkeypatch):
    monkeypatch.setattr(sys, "path", list(sys.path))
    monkeypatch.setattr(gb, "_installed", lambda dist="pymrm": "3.1.0")
    monkeypatch.setenv("GALLERY_REF", "unset")     # restored (removed) at teardown
    monkeypatch.delenv("GALLERY_REF")
    for name in ("model", "gallery_utils"):
        if name in sys.modules:
            monkeypatch.setitem(sys.modules, name, sys.modules[name])


def test_setup_in_checkout_puts_paths_and_fetches_nothing(monkeypatch, clean_state):
    _no_network(monkeypatch)
    stale = types.ModuleType("model")
    stale.__file__ = str(ROOT / "pages" / "some-other-page" / "model.py")
    monkeypatch.setitem(sys.modules, "model", stale)
    folder = gb.setup(PAGE, ">=3.1.0,<4", "gallery-2026.10.30")
    assert folder == ROOT / "pages" / PAGE
    assert sys.path[0] == str(ROOT / "pages" / PAGE)
    assert sys.path[1] == str(ROOT / "shared")
    assert "model" not in sys.modules
    assert gb.os.environ["GALLERY_REF"] == "gallery-2026.10.30"


def test_setup_without_checkout_offline_raises(monkeypatch, clean_state):
    _no_network(monkeypatch)
    monkeypatch.setenv("GALLERY_OFFLINE", "1")
    monkeypatch.setattr(gb, "_checkout_root", lambda page: None)
    with pytest.raises(FileNotFoundError, match="GALLERY_OFFLINE"):
        gb.setup(PAGE, ">=3.1.0,<4", "main")


def test_setup_without_checkout_downloads_to_temp(monkeypatch, clean_state, tmp_path):
    monkeypatch.delenv("GALLERY_OFFLINE", raising=False)
    monkeypatch.setattr(gb, "_checkout_root", lambda page: None)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.chdir(tmp_path)
    urls = []

    def fake(url, dest):
        urls.append(url)
        if url.endswith("/model.py"):
            raise urllib.error.HTTPError(url, 404, "Not Found", None, None)
        Path(dest).write_text("# fetched\n")

    monkeypatch.setattr(gb.urllib.request, "urlretrieve", fake)
    folder = gb.setup("X0.0-demo", ">=3.1.0,<4", "gallery-2026.10.30")
    assert folder == tmp_path / "cache" / "pymrm-gallery" / "gallery-2026.10.30" / "X0.0-demo"
    assert (folder / "gallery_utils.py").is_file() and not (folder / "model.py").exists()
    assert sys.path[0] == str(folder)
    assert all("/gallery-2026.10.30/" in u for u in urls)
    assert not (tmp_path / "gallery_utils.py").exists()   # never the working directory


# -- gallery_utils ----------------------------------------------------------------

def test_raw_url_follows_gallery_ref(monkeypatch):
    monkeypatch.delenv("GALLERY_OFFLINE", raising=False)
    monkeypatch.setenv("GALLERY_REF", "gallery-2026.10.30")
    url = gu._resolve("missing.csv", "X0.0-demo")
    assert url.endswith("/pymrm-gallery/gallery-2026.10.30/pages/X0.0-demo/data/missing.csv")
    monkeypatch.delenv("GALLERY_REF")
    assert "/pymrm-gallery/main/pages/" in gu._resolve("missing.csv", "X0.0-demo")


def test_offline_resolve_raises_for_missing_file(monkeypatch):
    monkeypatch.setenv("GALLERY_OFFLINE", "1")
    with pytest.raises(FileNotFoundError, match="GALLERY_OFFLINE"):
        gu._resolve("missing.csv", "X0.0-demo")
    with pytest.raises(FileNotFoundError, match="GALLERY_OFFLINE"):
        gu.load_meta("missing.csv", page="X0.0-demo")
    with pytest.raises(FileNotFoundError, match="GALLERY_OFFLINE"):
        gu.load_claims("X0.0-no-such-page")


def test_offline_resolve_still_finds_local_files(monkeypatch):
    monkeypatch.setenv("GALLERY_OFFLINE", "1")
    page = next(p for p in sorted((ROOT / "pages").iterdir())
                if (p / "data").is_dir() and any((p / "data").glob("*.csv")))
    csv = next((page / "data").glob("*.csv"))
    assert gu._resolve(csv.name, page.name) == csv


def test_numpy_simd_pin_matches_run_pages():
    src = (ROOT / "scripts" / "run_pages.py").read_text()
    assert f'NPY_PIN = "{gb.NPY_PIN}"' in src
    assert f'NPY_PIN = "{gb.NPY_PIN}"' in (ROOT / "scripts" / "run_kernels.py").read_text()


@pytest.mark.parametrize("machine, preset, loaded, expect", [
    ("x86_64", None, False, True),
    ("aarch64", None, False, False),     # the names exist on x86 only
    ("x86_64", "", False, False),        # a value already set is kept
    ("x86_64", None, True, False),       # too late: numpy reads it at import
])
def test_numpy_simd_pin_only_when_it_can_act(monkeypatch, machine, preset, loaded, expect):
    import platform
    monkeypatch.setattr(platform, "machine", lambda: machine)
    monkeypatch.setenv("NPY_DISABLE_CPU_FEATURES", "x")
    if preset is None:
        monkeypatch.delenv("NPY_DISABLE_CPU_FEATURES")
    else:
        monkeypatch.setenv("NPY_DISABLE_CPU_FEATURES", preset)
    modules = dict(sys.modules)
    modules.pop("numpy", None)
    if loaded:
        modules["numpy"] = types.ModuleType("numpy")
    monkeypatch.setattr(sys, "modules", modules)
    assert gb._pin_numpy_simd() is expect
    if expect:
        assert gb.os.environ["NPY_DISABLE_CPU_FEATURES"] == gb.NPY_PIN
