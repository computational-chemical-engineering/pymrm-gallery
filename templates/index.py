# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ---
# title: "Model name (Author, year)"
# description: "One sentence: what the page shows and why it matters."
# categories: [sec:X, struct:S0, tier:T0, data:tier0]
# date: 2026-01-01
# ---
#
# # Model name
#
# **Catalog ID:** `X0.0` · **Structures:** `S0` · **Tier:** T0
#
# One paragraph: the phenomenon, the published model, and what this page does
# with it.
#
# <!-- Page skeleton for the refreshed layout (docs/page-infrastructure.md).
# Copy this file to pages/<catalog-id>-<slug>/index.py, set the three literals of
# the bootstrap cell below (meta.yaml pymrm_min_version is the pymrm== version of
# requirements-lock.txt), then create the paired notebook with
#     jupytext --set-formats ipynb,py:percent pages/<page>/index.py
# and keep the two in step with `jupytext --sync`. Edit index.py, not the
# notebook. The model classes go in model.py next to it (NumPy-style
# docstrings, no plotting, no data loading), tested by test_model.py. -->

# %%
# Bootstrap: the pymrm this page was checked with, and the gallery files it was
# published with. Works in a checkout, in CI and on a fresh Colab VM.
GALLERY_REF = "main"  # replaced by the release tag at publication
PYMRM = ">=3.1.0,<4"  # equals meta.yaml pymrm_min_version, up to the next major
PAGE = "X0.0-replace-me"  # equals the page directory

import os, sys, urllib.request
from pathlib import Path
_shared = next((p / "shared" for p in (Path.cwd(), *Path.cwd().parents)
                if (p / "shared" / "gallery_bootstrap.py").is_file()), None)
if _shared is None:
    _shared = (Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
               / "pymrm-gallery" / GALLERY_REF)
    _shared.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        "https://raw.githubusercontent.com/computational-chemical-engineering/"
        f"pymrm-gallery/{GALLERY_REF}/shared/gallery_bootstrap.py",
        _shared / "gallery_bootstrap.py")
sys.path.insert(0, str(_shared))
sys.modules.pop("gallery_bootstrap", None)
from gallery_bootstrap import setup
setup(PAGE, PYMRM, GALLERY_REF)

# %%
# Imports. Packages beyond pymrm and the Colab defaults (numpy, scipy, pandas,
# matplotlib, pyyaml) are installed here with subprocess, not %pip.
import numpy as np
import matplotlib.pyplot as plt
from gallery_utils import load_data, load_meta, cite_data, Claims
# from model import MyModel     # the page's model.py, put on sys.path by setup()

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
claims = Claims(PAGE.split("-")[0])

# %% [markdown]
# ## Background
#
# The physical problem, the question the original authors asked, and why the
# model matters.

# %% [markdown]
# ## The published model
#
# The equations as published, with a table mapping each symbol to its name in
# the code.

# %% [markdown]
# ## Parameters and assumptions
#
# Every parameter with its unit and source; every assumption stated.

# %%
# Parameters, one per line, with unit and source in a comment.

# %% [markdown]
# ## The data
#
# Where the data come from (each dataset has a `data/*.meta.yaml` sidecar). A
# dataset borrowed from another page: list that page's findings about the rows
# used (AGENTS.md, "Data rules").

# %%
# data = load_data("author-year-figN.csv", page=PAGE)
# meta = load_meta("author-year-figN.csv", page=PAGE)
# print(cite_data(meta))

# %% [markdown]
# ## PyMRM implementation
#
# The discretisation and the reusable choices: boundary conditions (outward
# normal, physical equation in a comment), operators assembled once, NumJac
# stencil.

# %%
# model = MyModel(...)
# result = model.solve()

# %% [markdown]
# ## Results

# %%
# Figures against the data.

# %% [markdown]
# ## Validation
#
# At least one check that does not use the data: an analytical limit, a
# conservation balance, grid or time-step independence, or a physical bound.

# %%
# Register every number the page states, with how it is checked
# (docs/claims-design.md), and a break-table row that moves it:
# claims.value("nrmse_fig3", err, rtol=0.05)
# claims.breaks("rate constant 10 % high", "nrmse_fig3", err_broken)
# claims.says("0.031", "nrmse_fig3")

# %% [markdown]
# ## What pymrm adds
#
# Truthful: if the page only reproduces the original, say so.

# %% [markdown]
# ## Reuse
#
# What to change to adapt the model, and the related pages by catalog ID.

# %%
# Claims: prints every registered number, writes agreement.json for CI, and
# raises when a claim fails its own check or has no break row.
claims.report();   # the semicolon keeps the notebook from echoing the payload
