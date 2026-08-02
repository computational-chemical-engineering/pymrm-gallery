"""Numbered overlays and contact sheets for review."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, Rectangle

import glyphfit as G


def outline(ax, g, color, lw=1.0):
    sh, s, y, x = g["shape"], g["size"], g["y"], g["x"]
    if sh == "plus":
        ax.add_patch(Rectangle((x - s / 2, y - s / 2), s, s, fill=False,
                               ec=color, lw=lw))
    elif sh == "circle":
        ax.add_patch(Circle((x, y), s / 2, fill=False, ec=color, lw=lw))
    else:
        p = G.POLY[sh] * s
        ax.add_patch(Polygon(np.stack([p[:, 0] + x, -p[:, 1] + y], 1),
                             closed=True, fill=False, ec=color, lw=lw))


def full_overlay(img, glyphs, path, title, dpi=170, show_old=True, numbers=True):
    h, w = img.shape
    fig, ax = plt.subplots(figsize=(w / dpi, h / dpi + 0.5), dpi=dpi)
    ax.imshow(img, cmap="gray", interpolation="nearest")
    for i, g in enumerate(glyphs):
        col = "tab:red" if g.get("origin", "").startswith("new") else "tab:green"
        outline(ax, g, col, 1.1)
        ax.plot([g["x"]], [g["y"]], "+", color=col, ms=4, mew=0.7)
        if show_old and g.get("old_x") is not None:
            ax.plot([g["old_x"]], [g["old_y"]], "x", color="tab:orange", ms=3.5,
                    mew=0.7)
            ax.plot([g["old_x"], g["x"]], [g["old_y"], g["y"]], "-",
                    color="tab:orange", lw=0.5)
        if numbers:
            ax.text(g["x"] + g["size"] / 2 + 2, g["y"] - g["size"] / 2 - 2, str(i),
                    color=col, fontsize=3.4, ha="left", va="bottom")
    ax.set_title(title, fontsize=6)
    ax.set_axis_off()
    fig.tight_layout(pad=0.1)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def contact_sheet(img, glyphs, path, ncol=10, half=26, zoom=3, dpi=150,
                  title="", show_old=True):
    n = len(glyphs)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 1.05, nrow * 1.12), dpi=dpi)
    axes = np.atleast_2d(axes)
    for k in range(nrow * ncol):
        ax = axes[k // ncol, k % ncol]
        ax.set_axis_off()
        if k >= n:
            continue
        g = glyphs[k]
        y, x = int(round(g["y"])), int(round(g["x"]))
        sub = G.window(img, y, x, half, fill=1.0)
        ax.imshow(sub, cmap="gray", interpolation="nearest",
                  extent=(x - half - .5, x + half + .5, y + half + .5, y - half - .5))
        col = "tab:red" if g.get("origin", "").startswith("new") else "tab:green"
        outline(ax, g, col, 0.9)
        ax.plot([g["x"]], [g["y"]], "+", color=col, ms=5, mew=0.8)
        if show_old and g.get("old_x") is not None:
            ax.plot([g["old_x"]], [g["old_y"]], "x", color="tab:orange", ms=5, mew=0.8)
        ax.set_xlim(x - half, x + half)
        ax.set_ylim(y + half, y - half)
        ax.set_title(f"{k} {g['shape'][:5]} {g['score']:.2f}", fontsize=3.6, pad=1)
    fig.suptitle(title, fontsize=6)
    fig.tight_layout(pad=0.15, rect=(0, 0, 1, 0.985))
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
