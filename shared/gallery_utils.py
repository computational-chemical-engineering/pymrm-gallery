"""Shared helpers for pymrm-gallery pages.

Kept deliberately small. Every function must work identically when the notebook
runs locally, in CI, and on a fresh Google Colab VM with no repository checkout.
"""

from __future__ import annotations

import io
import json
import os
from pathlib import Path

import numpy as np

REPO_RAW = "https://raw.githubusercontent.com/computational-chemical-engineering/pymrm-gallery"


def _raw_base() -> str:
    """Raw URL of the gallery at ``GALLERY_REF`` (default ``main``), read per call.

    The bootstrap of a migrated page sets ``GALLERY_REF`` to the release tag the
    page was published with, so its data come from the same tag.
    """
    return f"{REPO_RAW}/{os.environ.get('GALLERY_REF', 'main')}"


def _offline() -> bool:
    """True when ``GALLERY_OFFLINE`` forbids fetching (set by run_pages.py and CI)."""
    return os.environ.get("GALLERY_OFFLINE", "").strip().lower() in ("1", "true", "yes")

__all__ = [
    "in_colab",
    "load_data",
    "load_meta",
    "cite_data",
    "rmse",
    "nrmse",
    "report_agreement",
    "Claims",
    "ClaimsFailed",
    "load_claims",
    "parse_literal",
]


def in_colab() -> bool:
    """True when running inside Google Colab."""
    return "google.colab" in str(os.environ.get("PYTHONPATH", "")) or (
        "COLAB_RELEASE_TAG" in os.environ
    )


def _page_dir() -> Path:
    """Directory of the page being executed.

    Quarto sets ``execute-dir: project``, so resolve relative to this file's
    grandparent when a page-local path cannot be determined.
    """
    return Path.cwd()


def _resolve(name: str, page: str | None):
    """Return a local path if it exists, otherwise a raw GitHub URL.

    The URL points at ``GALLERY_REF`` (default ``main``). With ``GALLERY_OFFLINE``
    set, a file that is not local raises instead.
    """
    candidates = []
    if page:
        candidates.append(Path(__file__).resolve().parents[1] / "pages" / page / "data" / name)
    candidates.append(_page_dir() / "data" / name)
    candidates.append(_page_dir() / name)
    for c in candidates:
        if c.is_file():
            return c
    if page is None:
        raise FileNotFoundError(
            f"{name!r} not found locally and no `page=` given, so the raw URL "
            f"cannot be constructed. Pass page='<catalog-id>-<slug>'."
        )
    if _offline():
        raise FileNotFoundError(
            f"{name!r} not found locally for page {page!r}, and GALLERY_OFFLINE is "
            f"set, so it is not fetched from GitHub.")
    return f"{_raw_base()}/pages/{page}/data/{name}"


def load_data(name: str, page: str | None = None, **kwargs):
    """Load a page dataset as a pandas DataFrame.

    Resolves to the local file when the repository is checked out, and to the
    raw GitHub URL otherwise (Colab). ``kwargs`` go to :func:`pandas.read_csv`.
    """
    import pandas as pd

    src = _resolve(name, page)
    return pd.read_csv(src, comment="#", **kwargs)


def load_meta(name: str, page: str | None = None) -> dict:
    """Load the provenance sidecar for a dataset (``<name>.meta.yaml``)."""
    import yaml

    src = _resolve(f"{Path(name).stem}.meta.yaml", page)
    if isinstance(src, Path):
        text = src.read_text(encoding="utf-8")
    else:
        from urllib.request import urlopen

        with urlopen(src) as fh:  # noqa: S310 - fixed, trusted host
            text = fh.read().decode("utf-8")
    return yaml.safe_load(text)


def cite_data(meta: dict) -> str:
    """One-line human-readable provenance string for printing under a figure."""
    s = meta.get("source", {})
    a = meta.get("acquisition", {})
    authors = "; ".join(s.get("authors", [])) if isinstance(s.get("authors"), list) else s.get("authors", "")
    bits = [f"{authors} ({s.get('year', 'n.d.')})", s.get("container", "")]
    if s.get("doi"):
        bits.append(f"doi:{s['doi']}")
    if a.get("method") == "digitised":
        bits.append(f"digitised from {a.get('figure', 'figure')}")
    elif a.get("method"):
        bits.append(str(a["method"]))
    return " — ".join(b for b in bits if b)


def rmse(model, data) -> float:
    """Root-mean-square error between model and data arrays."""
    model, data = np.asarray(model, float), np.asarray(data, float)
    return float(np.sqrt(np.mean((model - data) ** 2)))


def nrmse(model, data) -> float:
    """RMSE normalised by the range of the data (dimensionless)."""
    data = np.asarray(data, float)
    span = np.ptp(data)
    return rmse(model, data) / span if span > 0 else np.nan


def report_agreement(page_id: str, metrics: dict, write: bool = True) -> dict:
    """Print agreement metrics and write ``agreement.json`` for CI regression checks.

    CI compares this file against a stored baseline so that a change in pymrm
    cannot silently degrade a page while the site still renders cached output.
    """
    payload = {"page": page_id, "metrics": {k: float(v) for k, v in metrics.items()}}
    width = max(len(k) for k in metrics) if metrics else 0
    print(f"Agreement metrics for {page_id}:")
    for k, v in payload["metrics"].items():
        print(f"  {k:<{width}} = {v:.5g}")
    if write:
        try:
            Path("agreement.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        except OSError:
            pass  # read-only environment (e.g. some Colab paths); metrics still printed
    return payload


# ---------------------------------------------------------------------------
# Claims: every published number declares what it is (docs/claims-design.md)
# ---------------------------------------------------------------------------

KINDS = ("value", "bound", "exact")
BASES = ("identity", "tolerance", "theory", "chosen")
_UPPER = ("below", "under", "less than", "at most")
_LOWER = ("above", "over", "more than", "at least")
_NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


class ClaimsFailed(AssertionError):
    """Raised by :meth:`Claims.report` listing every failed check at once."""


def parse_literal(text: str) -> dict:
    """Parse a number as written in prose, for :meth:`Claims.says`.

    Grammar: an optional bound phrase (``below``, ``under``, ``less than``,
    ``at most``, ``above``, ``over``, ``more than``, ``at least``, or
    ``between L1 and L2``), then a number with an optional exponent and an
    optional trailing ``%``. Returns ``op`` (``"eq"``, ``"max"``, ``"min"`` or
    ``"between"``), ``nums`` (floats, already multiplied by the ``%`` scale),
    ``unit`` (one unit in the last written place of the last number, scaled the
    same way) and ``integer`` (a bare integer with no decimal point or exponent).
    """
    import re

    s = " ".join(text.strip().split())
    low = s.lower()
    m = re.fullmatch(rf"between\s+({_NUM})\s*(%?)\s+and\s+({_NUM})\s*(%?)", low)
    if m:
        pct = 0.01 if (m.group(2) or m.group(4)) else 1.0
        a, b = float(m.group(1)) * pct, float(m.group(3)) * pct
        return {"op": "between", "nums": [min(a, b), max(a, b)],
                "unit": _unit(m.group(3)) * pct, "integer": False}
    op = "eq"
    for words, kind in ((_UPPER, "max"), (_LOWER, "min")):
        for w in words:
            if low.startswith(w + " "):
                op, low = kind, low[len(w) + 1:].strip()
                break
        if op != "eq":
            break
    m = re.fullmatch(rf"({_NUM})\s*(%?)", low)
    if not m:
        raise ValueError(f"cannot parse the prose literal {text!r}")
    pct = 0.01 if m.group(2) else 1.0
    integer = op == "eq" and pct == 1.0 and re.fullmatch(r"[+-]?\d+", m.group(1)) is not None
    return {"op": op, "nums": [float(m.group(1)) * pct],
            "unit": _unit(m.group(1)) * pct, "integer": integer}


def _unit(num: str) -> float:
    """One unit in the last written place: ``0.7186`` -> 1e-4, ``3.9e-5`` -> 1e-6."""
    mant, _, exp = num.lower().partition("e")
    decimals = len(mant.partition(".")[2])
    return 10.0 ** (int(exp or 0) - decimals)


def _env() -> dict:
    """The arithmetic a run was produced with, recorded beside its metrics."""
    import importlib.metadata as md

    env = {}
    for pkg in ("pymrm", "numpy", "scipy"):
        try:
            env[pkg] = md.version(pkg)
        except md.PackageNotFoundError:
            env[pkg] = None
    env["openblas_coretype"] = os.environ.get("OPENBLAS_CORETYPE", "default")
    env["threads"] = os.environ.get("OPENBLAS_NUM_THREADS", "default")
    return env


def _native(x):
    """JSON-safe scalar: numpy numbers become Python ones; labels stay strings."""
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, str):
        return x
    return float(x)


class Claims:
    """Register every number a page publishes, with how it is to be checked.

    Kinds (``docs/claims-design.md``):

    ``value``
        A physical, data or discretisation result above round-off. CI compares
        it with the committed baseline: ``|x - x0| <= rtol*|x0| + atol``.
    ``bound``
        A residual, closure, identity, two-route gap, solver-tolerance effect or
        a convergence order. Checked at run time and in CI as
        ``min <= x <= max``; the value itself is information only. Every bound
        states its ``basis``: ``identity``, ``tolerance``, ``theory`` or
        ``chosen``.
    ``exact``
        An integer, boolean or label. CI requires equality. When it was decided
        by comparing values with a threshold, record the margin
        (:meth:`count`, or ``margin=``/``min_margin=``).

    Every claim that is not ``structural`` needs at least one break-table row
    (:meth:`breaks`) that fails the claim's own check, and :meth:`report`
    raises otherwise.

    Examples
    --------
    >>> claims = Claims("B1.6")
    >>> claims.value("vs_shooting_reference_n800", 3.944e-5)
    >>> claims.breaks("reference off by 1 %", "vs_shooting_reference_n800", 1.2e-2)
    >>> claims.says("3.94e-5", "vs_shooting_reference_n800")
    >>> payload = claims.report(write=False)  # doctest: +SKIP
    """

    def __init__(self, page: str):
        self.page = page
        self.metrics: dict[str, dict] = {}
        self.quoted: dict[str, list[str]] = {}
        self._rows: dict[str, list[tuple[str, bool]]] = {}
        self._problems: list[str] = []

    # -- registration ---------------------------------------------------------

    def _add(self, name: str, entry: dict) -> None:
        if name in self.metrics:
            raise ValueError(f"claim {name!r} registered twice")
        v = entry["value"]
        if isinstance(v, float) and not np.isfinite(v):
            self._problems.append(f"{name}: value is {v}")
        self.metrics[name] = {k: w for k, w in entry.items() if w is not None}

    def value(self, name, x, *, rtol=0.05, atol=0.0, unit=None, note=None):
        """Register a value; CI holds it to ``rtol*|x0| + atol`` of the baseline."""
        x = float(x)
        if x == 0.0 and atol <= 0.0:
            raise ValueError(f"{name}: a value of exactly 0 needs atol > 0; "
                             "a quantity that is zero by construction is a bound")
        self._add(name, {"kind": "value", "value": x, "rtol": float(rtol),
                         "atol": float(atol) or None, "unit": unit, "note": note})

    def bound(self, name, x, *, min=None, max=None, basis, structural=False,
              unit=None, note=None):
        """Register a bound ``min <= x <= max``, checked at run time and in CI."""
        if min is None and max is None:
            raise ValueError(f"{name}: a bound needs min, max or both")
        if basis not in BASES:
            raise ValueError(f"{name}: basis must be one of {BASES}, not {basis!r}")
        if basis == "chosen" and not note:
            raise ValueError(f"{name}: a chosen bound needs a note saying why")
        x = float(x)
        self._add(name, {"kind": "bound", "value": x,
                         "min": None if min is None else float(min),
                         "max": None if max is None else float(max),
                         "basis": basis, "structural": True if structural else None,
                         "unit": unit, "note": note})
        if not _in_bound(x, self.metrics[name]):
            self._problems.append(f"{name}: {x:.6g} violates {_bound_text(self.metrics[name])}")

    def exact(self, name, x, *, margin=None, min_margin=None, structural=False, note=None):
        """Register an exact result. Give ``margin`` and ``min_margin`` when ``x``
        was decided against a threshold, so a result inside the noise raises."""
        if (margin is None) != (min_margin is None):
            raise ValueError(f"{name}: give margin and min_margin together")
        self._add(name, {"kind": "exact", "value": _native(x),
                         "margin": None if margin is None else float(margin),
                         "min_margin": None if min_margin is None else float(min_margin),
                         "structural": True if structural else None, "note": note})
        if margin is not None and not margin >= min_margin:
            self._problems.append(f"{name}: margin {margin:.3g} to its threshold is below "
                                  f"min_margin {min_margin:.3g}, so the result is decided "
                                  "by noise")

    def count(self, name, values, threshold, *, min_margin, above=True,
              relative=None, note=None):
        """Register how many ``values`` lie above (or below) ``threshold``.

        The margin is the smallest distance of any value to the threshold,
        relative to ``|threshold|`` unless the threshold is 0 or ``relative`` is
        False. Returns the count.
        """
        v = np.asarray(values, dtype=float).ravel()
        n = int(np.sum(v > threshold) if above else np.sum(v < threshold))
        if relative is None:
            relative = threshold != 0
        gap = np.abs(v - threshold)
        margin = float(gap.min() / (abs(threshold) if relative else 1.0)) if v.size else np.inf
        self.exact(name, n, margin=margin, min_margin=min_margin, note=note)
        return n

    # -- break table and prose ------------------------------------------------

    def breaks(self, row: str, name: str, broken) -> bool:
        """Record a break-table row: the value claim ``name`` takes with the
        defect injected. Returns whether the claim's own check catches it."""
        if name not in self.metrics:
            raise KeyError(f"break row {row!r} names unknown claim {name!r}")
        caught = not _passes(self.metrics[name], self.metrics[name]["value"], broken)
        self._rows.setdefault(name, []).append((row, caught))
        return caught

    def says(self, literal: str, name: str) -> bool:
        """Check a number as written in the prose against claim ``name``.

        A number passes when within one unit of its last written place of the
        computed value. A bound phrase (``below 1e-8``) must hold for every run
        that passes the claim's check: for a bound, at least as loose as the
        declared bound; for a value, its whole tolerance band inside the phrase.
        """
        if name not in self.metrics:
            raise KeyError(f"prose literal {literal!r} names unknown claim {name!r}")
        c = self.metrics[name]
        lit = parse_literal(literal)
        ok = _says_ok(lit, c)
        self.quoted.setdefault(name, []).append(literal)
        if not ok:
            self._problems.append(f"{name}: prose says {literal!r}, computed {c['value']!r}"
                                  + (f" with {_bound_text(c)}" if c["kind"] == "bound" else ""))
        return ok

    def __getitem__(self, name):
        return self.metrics[name]["value"]

    # -- output ---------------------------------------------------------------

    def report(self, write: bool = True, require_breaks: bool = True) -> dict:
        """Print the claims, write ``agreement.json`` (format 2), raise on failure.

        ``require_breaks=False`` is for a page migrated before it has a break
        table; the file records it and ``scripts/check_claims.py`` lists the page.
        """
        problems = list(self._problems)
        unbroken = [name for name, c in self.metrics.items() if not c.get("structural")
                    and not any(caught for _, caught in self._rows.get(name, []))]
        if require_breaks:
            problems += [f"{name}: no break row fails its check (add one, or mark the "
                         "claim structural and say why)" for name in unbroken]
        payload = {
            "page": self.page,
            "format": 2,
            "env": _env(),
            "metrics": self.metrics,
            "quoted": self.quoted,
            "breaks": {n: [r for r, caught in rows if caught] for n, rows in self._rows.items()},
        }
        if not require_breaks:
            # A page migrated before its break table: visible in the file and
            # listed by scripts/check_claims.py until the table exists.
            payload["breaks_required"] = False
            payload["unbroken"] = unbroken
        width = max((len(k) for k in self.metrics), default=0)
        print(f"Claims for {self.page}:")
        for k, c in self.metrics.items():
            v = c["value"]
            shown = f"{v:.5g}" if isinstance(v, float) else repr(v)
            if c["kind"] == "bound":
                # a bound publishes its limits; the value is this run's, and its
                # digits differ between machines, so it is shown as such
                shown = f"(here {v:.1e})"
            rule = {"value": lambda: f"rtol {c['rtol']:g}" + (f", atol {c['atol']:g}" if "atol" in c else ""),
                    "bound": lambda: _bound_text(c) + f" [{c['basis']}]",
                    "exact": lambda: "exact" + (f", margin {c['margin']:.3g}" if "margin" in c else "")}[c["kind"]]()
            print(f"  {k:<{width}} = {shown:<16} {c['kind']:<5}  {rule}")
        if write:
            try:
                Path("agreement.json").write_text(json.dumps(payload, indent=2) + "\n",
                                                  encoding="utf-8")
            except OSError:
                pass  # read-only environment (e.g. some Colab paths); claims still printed
        if problems:
            raise ClaimsFailed(f"{len(problems)} claim check(s) failed for {self.page}:\n  "
                               + "\n  ".join(problems))
        return payload


def _in_bound(x, c) -> bool:
    lo, hi = c.get("min"), c.get("max")
    return np.isfinite(x) and (lo is None or x >= lo) and (hi is None or x <= hi)


def _bound_text(c) -> str:
    lo, hi = c.get("min"), c.get("max")
    if lo is not None and hi is not None:
        return f"bound [{lo:g}, {hi:g}]"
    return f"bound <= {hi:g}" if hi is not None else f"bound >= {lo:g}"


def _passes(c: dict, x0, x) -> bool:
    """Whether ``x`` passes claim ``c`` against the reference value ``x0``."""
    if c["kind"] == "value":
        return bool(np.isfinite(x) and abs(float(x) - float(x0))
                    <= c["rtol"] * abs(float(x0)) + c.get("atol", 0.0))
    if c["kind"] == "bound":
        return bool(_in_bound(float(x), c))
    return _native(x) == _native(x0)


def _says_ok(lit: dict, c: dict) -> bool:
    x = c["value"]
    op, nums, unit = lit["op"], lit["nums"], lit["unit"]
    if c["kind"] == "exact":
        if isinstance(x, str):
            return False
        x = float(x)
        if op == "between":
            return nums[0] <= x <= nums[1]
        return {"eq": x == nums[0], "max": x < nums[0], "min": x > nums[0]}[op]
    if op == "eq":
        return abs(nums[0] - float(x)) <= unit
    if c["kind"] == "bound":
        lo, hi = c.get("min"), c.get("max")
        if op == "max":
            return hi is not None and nums[0] >= hi
        if op == "min":
            return lo is not None and nums[0] <= lo
        return lo is not None and hi is not None and nums[0] <= lo and nums[1] >= hi
    # a value quoted with a phrase: the whole tolerance band must sit inside it
    band = c["rtol"] * abs(x) + c.get("atol", 0.0)
    lo_band, hi_band = x - band, x + band
    if op == "max":
        return hi_band < nums[0]
    if op == "min":
        return lo_band > nums[0]
    return nums[0] <= lo_band and hi_band <= nums[1]


def load_claims(page: str) -> dict:
    """Another page's committed metrics as ``{name: value}``, either format.

    ``page`` is the page directory name (``"F3.1-..."``). Resolves the local
    checkout first and the raw GitHub URL otherwise. A comparison of a live
    number with one committed on another machine is a ``bound`` with a measured
    basis, never a tight equality.
    """
    for base in (Path.cwd(), *Path.cwd().parents, Path(__file__).resolve().parents[1]):
        p = base / "pages" / page / "agreement.json"
        if p.is_file():
            data = json.loads(p.read_text(encoding="utf-8"))
            break
    else:
        if _offline():
            raise FileNotFoundError(
                f"pages/{page}/agreement.json not found locally, and GALLERY_OFFLINE "
                f"is set, so it is not fetched from GitHub.")
        from urllib.request import urlopen

        with urlopen(f"{_raw_base()}/pages/{page}/agreement.json") as fh:  # noqa: S310
            data = json.load(fh)
    return {k: (v["value"] if isinstance(v, dict) else v) for k, v in data["metrics"].items()}
