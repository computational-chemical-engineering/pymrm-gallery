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

RAW_BASE = "https://raw.githubusercontent.com/computational-chemical-engineering/pymrm-gallery/main"

__all__ = [
    "in_colab",
    "load_data",
    "load_meta",
    "cite_data",
    "rmse",
    "nrmse",
    "report_agreement",
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
    """Return a local path if it exists, otherwise a raw GitHub URL."""
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
    return f"{RAW_BASE}/pages/{page}/data/{name}"


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
