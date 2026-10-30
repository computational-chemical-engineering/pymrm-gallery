"""Environment setup for a migrated gallery page (docs/page-infrastructure.md, section 2).

The first code cell of a migrated page imports this module and calls
:func:`setup`. It runs before numpy or pymrm may exist, so it uses the Python
standard library only (``tests/test_gallery_bootstrap.py`` enforces that).

``setup(page, pymrm_spec, ref)``

1. makes sure the installed pymrm satisfies ``pymrm_spec``;
2. sets ``os.environ["GALLERY_REF"]`` so ``gallery_utils`` builds its URLs at ``ref``;
3. puts ``shared/`` and the page folder on ``sys.path`` in a checkout, or downloads
   ``gallery_utils.py`` and the page's ``model.py`` at ``ref`` into the user's
   cache folder otherwise (never the working directory, never shared /tmp);
4. with ``GALLERY_OFFLINE=1`` fetches nothing and raises instead.
"""

from __future__ import annotations

import importlib.metadata
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request
import warnings
from pathlib import Path

REPO_RAW = "https://raw.githubusercontent.com/computational-chemical-engineering/pymrm-gallery"

__all__ = ["setup", "parse_version", "satisfies", "offline", "cache_dir"]


# ---------------------------------------------------------------------------
# Versions (``packaging`` may be absent on a bare VM)
# ---------------------------------------------------------------------------

_VERSION = re.compile(
    r"^\s*v?(?P<release>\d+(?:\.\d+)*)"
    r"(?:[-_.]?(?P<pre_l>a|b|c|rc|alpha|beta|pre|preview)[-_.]?(?P<pre_n>\d*))?"
    r"(?:[-_.]?(?:post|rev|r)[-_.]?(?P<post>\d*))?"
    r"(?:[-_.]?dev[-_.]?(?P<dev>\d*))?"
    r"(?:\+[a-z0-9.]+)?\s*$",
    re.IGNORECASE,
)
_PRE_RANK = {"a": 0, "alpha": 0, "b": 1, "beta": 1, "c": 2, "rc": 2, "pre": 2, "preview": 2}
_INF = float("inf")


def parse_version(text: str) -> tuple:
    """Sortable key of a PEP 440 version string.

    Pre-releases sort below their final release (``3.1.0rc1 < 3.1.0``) and above
    the previous one (``3.2.0rc1 > 3.1.0``); ``.devN`` sorts below, ``.postN``
    above. Trailing zeros of the release do not matter (``3.1 == 3.1.0``).
    """
    m = _VERSION.match(text)
    if m is None:
        raise ValueError(f"not a version: {text!r}")
    release = [int(x) for x in m["release"].split(".")]
    while len(release) > 1 and release[-1] == 0:
        release.pop()
    if m["pre_l"]:
        pre = (_PRE_RANK[m["pre_l"].lower()], int(m["pre_n"] or 0))
    elif m["dev"] is not None and m["post"] is None:
        pre = (-1, 0)                   # 3.1.0.dev1 sorts before 3.1.0a1
    else:
        pre = (_INF, 0)                 # a final release
    post = int(m["post"] or 0) if m["post"] is not None else -1
    dev = int(m["dev"] or 0) if m["dev"] is not None else _INF
    return (tuple(release), pre, post, dev)


def _cmp(a: str, b: str) -> int:
    ka, kb = parse_version(a), parse_version(b)
    return (ka > kb) - (ka < kb)


def satisfies(version: str, spec: str) -> bool:
    """True when ``version`` meets every clause of ``spec`` (e.g. ``">=3.1.0,<4"``).

    Pre-releases are admitted, as ``packaging`` does with ``prereleases=True``:
    ``3.2.0rc1`` satisfies ``>=3.1.0,<4``. Supported operators: ``>= > <= < == !=``.
    """
    for clause in filter(None, (c.strip() for c in spec.split(","))):
        m = re.match(r"^(>=|<=|==|!=|>|<)\s*(\S+)$", clause)
        if m is None:
            raise ValueError(f"unsupported version clause {clause!r} in {spec!r}")
        op, ref = m.groups()
        c = _cmp(version, ref)
        ok = {">=": c >= 0, ">": c > 0, "<=": c <= 0, "<": c < 0,
              "==": c == 0, "!=": c != 0}[op]
        kv, kr = parse_version(version), parse_version(ref)
        if op == "<" and ok and kr[1] == (_INF, 0) and kv[0] == kr[0]:
            ok = False      # PEP 440: <4 excludes 4.0.0rc1 unless 4 itself is a pre-release
        if op == ">" and ok and kr[2] == -1 and kv[0] == kr[0] and kv[1] == kr[1]:
            ok = False      # PEP 440: >3.1 excludes 3.1.post1
        if not ok:
            return False
    return True


def _above_spec(version: str, spec: str) -> bool:
    """True when ``version`` fails ``spec`` only through an upper bound."""
    lowers = [c for c in spec.split(",") if c.strip().startswith(">")]
    return satisfies(version, ",".join(lowers)) and not satisfies(version, spec)


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------

def offline() -> bool:
    """True when ``GALLERY_OFFLINE`` forbids fetching (run_pages.py, CI model tests)."""
    return os.environ.get("GALLERY_OFFLINE", "").strip().lower() in ("1", "true", "yes")


def _installed(dist: str = "pymrm") -> str | None:
    try:
        return importlib.metadata.version(dist)
    except importlib.metadata.PackageNotFoundError:
        return None


def _pip_install(spec: str) -> None:
    if offline():
        raise RuntimeError(
            f"pymrm{spec} is needed and GALLERY_OFFLINE is set, so it is not installed.")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", f"pymrm{spec}"])
    importlib.invalidate_caches()


def _ensure_pymrm(spec: str, checkout: bool = False) -> None:
    have = _installed()
    if have is None:
        _pip_install(spec)
        return
    if satisfies(have, spec):
        return
    if _above_spec(have, spec):
        warnings.warn(
            f"pymrm {have} is newer than this page was checked with ({spec}); "
            f"keeping it.", stacklevel=3)
        return
    if checkout:
        raise RuntimeError(
            f"pymrm {have} is installed but this page needs pymrm{spec}. This is a "
            f"gallery checkout, so the bootstrap does not reinstall pymrm into your "
            f"development environment: install a matching pymrm (see "
            f"requirements-lock.txt) and restart the kernel.")
    if "pymrm" in sys.modules:
        raise RuntimeError(
            f"pymrm {have} is loaded but this page needs pymrm{spec}. Restart the "
            f"runtime (kernel) and run this cell first.")
    _pip_install(spec)


def _checkout_root(page: str) -> Path | None:
    """The gallery checkout this module lives in, when it holds ``pages/<page>``."""
    here = Path(__file__).resolve().parent
    if here.name == "shared" and (here.parent / "pages" / page).is_dir():
        return here.parent
    return None


def cache_dir(ref: str) -> Path:
    """``$XDG_CACHE_HOME/pymrm-gallery/<ref>``, or ``~/.cache/pymrm-gallery/<ref>``.

    A per-user folder: a predictable path in a shared ``/tmp`` placed at the
    front of ``sys.path`` would let another user plant modules.
    """
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "pymrm-gallery" / ref


def _download_dir(page: str, ref: str) -> Path:
    return cache_dir(ref) / page


def _fetch(url: str, dest: Path, required: bool) -> bool:
    try:
        urllib.request.urlretrieve(url, dest)  # noqa: S310 - fixed, trusted host
    except urllib.error.HTTPError as exc:
        if exc.code == 404 and not required:
            return False
        raise
    return True


def _drop_foreign(name: str, folder: Path) -> None:
    """Forget module ``name`` when it was loaded from somewhere other than ``folder``."""
    mod = sys.modules.get(name)
    if mod is None:
        return
    origin = getattr(mod, "__file__", None)
    if origin is None or Path(origin).resolve().parent != folder.resolve():
        del sys.modules[name]


def _front_of_path(*folders: Path) -> None:
    for folder in reversed(folders):
        s = str(folder)
        while s in sys.path:
            sys.path.remove(s)
        sys.path.insert(0, s)


def setup(page: str, pymrm_spec: str, ref: str) -> Path:
    """Prepare the interpreter for gallery page ``page`` at gallery ref ``ref``.

    Parameters
    ----------
    page : str
        Page directory name, e.g. ``"A2.4-tanks-in-series"``.
    pymrm_spec : str
        Version specifier of pymrm, e.g. ``">=3.1.0,<4"``.
    ref : str
        Gallery git ref the page was published with (``"main"`` before release).

    Returns
    -------
    pathlib.Path
        The folder ``model`` is imported from: the page folder in a checkout,
        the download folder ``cache_dir(ref) / page`` otherwise.

    Raises
    ------
    RuntimeError
        When pymrm is below the specifier in a checkout (never reinstalled
        there) or already loaded (restart needed), or must be installed while
        ``GALLERY_OFFLINE`` is set.
    FileNotFoundError
        When there is no checkout and ``GALLERY_OFFLINE`` is set.
    """
    root = _checkout_root(page)
    _ensure_pymrm(pymrm_spec, checkout=root is not None)
    os.environ["GALLERY_REF"] = ref

    if root is not None:
        page_dir = root / "pages" / page
        _front_of_path(page_dir, root / "shared")
        _drop_foreign("model", page_dir)
        _drop_foreign("gallery_utils", root / "shared")
        return page_dir

    if offline():
        raise FileNotFoundError(
            f"no gallery checkout holds pages/{page}, and GALLERY_OFFLINE is set, "
            f"so gallery_utils.py and model.py are not fetched.")
    dest = _download_dir(page, ref)
    dest.mkdir(parents=True, exist_ok=True)
    base = f"{REPO_RAW}/{ref}"
    _fetch(f"{base}/shared/gallery_utils.py", dest / "gallery_utils.py", required=True)
    if not _fetch(f"{base}/pages/{page}/model.py", dest / "model.py", required=False):
        (dest / "model.py").unlink(missing_ok=True)
    _front_of_path(dest)
    _drop_foreign("model", dest)
    _drop_foreign("gallery_utils", dest)
    return dest
