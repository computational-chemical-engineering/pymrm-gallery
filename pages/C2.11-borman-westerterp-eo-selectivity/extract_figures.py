#!/usr/bin/env python3
"""Regenerate data/borman-1995-fig6.csv from the source PDF.

The PDF is not in the repository and must not be (AGENTS.md). Point this at your
own lawfully obtained copy:

    python extract_figures.py <your copy>.pdf

Requires poppler's ``pdftoppm`` on PATH, plus numpy and pillow.

Method
------
Figure 6 (journal page 56, PDF page 8) shows three smooth curves the authors
calculated with model III: selectivity to ethene oxide against the partial
pressure of CO2 (solid), ethene oxide (dashed) and water (dash-dot). There are
no markers. The page is rendered at 400 dpi in grey and the figure cropped.

Axes: the frame lines are the longest ink rows and columns of the crop. The
calibration is a least-squares line through the frame edges and every tick mark
found on the inner side of the left axis (one per selectivity point, 65 to 77)
and of the bottom axis (one per 0.01 bar, 0 to 0.16); a tick is a short ink run
3 px inside the frame, no wider than 3 px, within a tenth of a unit of where
the frame alone puts it. The largest tick residual is printed.

Curves: every 6th pixel column is scanned between the frame lines for ink runs
no thicker than 10 px (a line, not a letter); the centre row of each run is one
point. Columns are kept only where the number of runs identifies the curves by
their order on the two flanks of the crossing region: three runs left of
0.045 bar (EO on top, then W, then C), three runs between 0.085 and 0.119 bar (C,
W, EO), two runs between 0.122 and 0.132 bar (C, W; EO ends at 0.12) and one run
right of 0.152 bar (C). Excluded: 18 px next to each frame line (the ticks on all
four sides reach about 12 px into the plot), the crossing
region 0.045 to 0.085 bar, and the columns of the labels "EO" (0.018 to 0.038
bar) and "C", "W" (0.132 to 0.152 bar). A column where a dash gap or the label
changes the count is skipped, not guessed.

Reading band: half the median line thickness plus the largest y tick residual
plus the largest x tick residual times the steepest slope of a quadratic
through each curve's points, all in
selectivity points, rounded up to 0.01; printed for the sidecar.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

DPI = 400
PAGE = 8
CROP = (440, 1760, 1560, 2440)       # (col_lo, row_lo, col_hi, row_hi) in the render
Y_LO, Y_HI = 65.0, 77.0               # selectivity at the bottom and top frame, %
X_LO, X_HI = 0.0, 0.16                # partial pressure at the left and right frame, bar
STEP = 6
MAX_RUN = 10
AXIS_MARGIN = 18                      # px kept clear of every frame line (ticks reach 12 px)
EXCLUDE = [(0.018, 0.038), (0.045, 0.085), (0.132, 0.152)]
FLANKS = [  # (x_lo, x_hi, runs, curve names top to bottom)
    (0.0, 0.045, 3, ("EO", "W", "C")),
    (0.085, 0.119, 3, ("C", "W", "EO")),
    (0.122, 0.132, 2, ("C", "W")),
    (0.152, 0.161, 1, ("C",)),
]


def render(pdf: Path, tmp: Path) -> np.ndarray:
    stem = tmp / "p"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(PAGE), "-l", str(PAGE),
                    "-gray", "-png", str(pdf), str(stem)], check=True)
    png = next(tmp.glob("p*.png"))
    return np.asarray(Image.open(png).convert("L"))


def runs(mask: np.ndarray) -> list[tuple[float, int]]:
    """Centre index and length of every run of True."""
    e = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    return [(0.5 * (a + b - 1), int(b - a)) for a, b in zip(e[::2], e[1::2])]


def frame(ink: np.ndarray) -> tuple[float, float, float, float]:
    """Left, right, top and bottom frame lines (centres, px)."""
    h, w = ink.shape
    rows = [c for c, _ in runs(ink.sum(1) > 0.6 * w)]
    top, bottom = rows[0], rows[-1]
    span = ink[int(top) + 3:int(bottom) - 2].sum(0)
    cols = [c for c, _ in runs(span > 0.6 * (bottom - top - 5))]
    cols = [c for c in cols if c > 100]          # the rotated axis title lies left of 100 px
    return cols[0], cols[-1], top, bottom


def calibrate(ink, left, right, top, bottom):
    """Least-squares pixel-to-value lines through the frame and the tick marks."""
    # y: ticks on the inner side of the left axis, one per selectivity point
    px_y, val_y = [top, bottom], [Y_HI, Y_LO]
    per = (bottom - top) / (Y_HI - Y_LO)
    for c, n in runs(ink[int(top) + 3:int(bottom) - 2, int(left) + 3]):
        r = c + int(top) + 3
        v = Y_HI - (r - top) / per
        if n <= 3 and abs(v - round(v)) < 0.1:
            px_y.append(r)
            val_y.append(float(round(v)))
    # x: ticks on the inner side of the bottom axis, one per 0.01 bar
    px_x, val_x = [left, right], [X_LO, X_HI]
    per_x = (right - left) / ((X_HI - X_LO) / 0.01)
    for c, n in runs(ink[int(bottom) - 3, int(left) + 5:int(right) - 4]):
        col = c + int(left) + 5
        v = (col - left) / per_x
        if n <= 3 and abs(v - round(v)) < 0.1:
            px_x.append(col)
            val_x.append(0.01 * round(v))
    by = np.polyfit(px_y, val_y, 1)
    bx = np.polyfit(px_x, val_x, 1)
    res_y = np.abs(np.polyval(by, px_y) - val_y)
    res_x = np.abs(np.polyval(bx, px_x) - val_x)
    return bx, by, len(px_x) - 2, len(px_y) - 2, float(res_x.max()), float(res_y.max())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", type=Path)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent / "data" / "borman-1995-fig6.csv")
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as t:
        img = render(args.pdf, Path(t))
    c0, r0, c1, r1 = CROP
    ink = img[r0:r1, c0:c1] < 128
    left, right, top, bottom = frame(ink)
    bx, by, nx, ny, rx, ry = calibrate(ink, left, right, top, bottom)
    pts, widths = [], []
    for col in range(int(left) + AXIS_MARGIN, int(right) - AXIS_MARGIN + 1, STEP):
        x = float(np.polyval(bx, col))
        if any(lo <= x <= hi for lo, hi in EXCLUDE):
            continue
        found = [(c + int(top) + AXIS_MARGIN, n)
                 for c, n in runs(ink[int(top) + AXIS_MARGIN:int(bottom) - AXIS_MARGIN, col])]
        if any(n > MAX_RUN for _, n in found):
            continue
        for lo, hi, k, names in FLANKS:
            if lo <= x < hi and len(found) == k:
                for name, (row, n) in zip(names, found):
                    pts.append((name, x, float(np.polyval(by, row))))
                    widths.append(n)
    per_point = 1.0 / abs(by[0])                 # px per selectivity point
    half_width = 0.5 * float(np.median(widths)) / per_point
    slope = 0.0
    for name in ("C", "EO", "W"):
        xy = np.array([(x, y) for n, x, y in pts if n == name])
        d = np.polyder(np.polyfit(xy[:, 0], xy[:, 1], 2))   # smooth slope, not pixel noise
        slope = max(slope, float(np.max(np.abs(np.polyval(d, xy[:, 0])))))
    band = np.ceil((half_width + ry + rx * slope) * 100.0) / 100.0
    with open(args.out, "w") as f:
        f.write("curve,P_bar,selectivity_pct\n")
        for name, x, y in pts:
            f.write(f"{name},{x:.5f},{y:.3f}\n")
    count = {n: sum(p[0] == n for p in pts) for n in ("C", "EO", "W")}
    print(f"frame px: left {left:.1f} right {right:.1f} top {top:.1f} bottom {bottom:.1f}")
    print(f"ticks used: x {nx}, y {ny}; largest residual x {rx:.2e} bar, y {ry:.3f} points")
    print(f"{len(pts)} points {count}; median line {np.median(widths):.0f} px "
          f"({per_point:.1f} px per point); steepest slope {slope:.0f} points/bar; "
          f"reading band {band:.2f} points")


if __name__ == "__main__":
    main()
