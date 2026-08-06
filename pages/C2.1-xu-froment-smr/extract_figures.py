#!/usr/bin/env python3
"""Regenerate data/xu-froment-1989-conversion.csv from the source PDF.

The PDF is not in the repository and must not be (AGENTS.md). Point this at your
own lawfully obtained copy:

    python extract_figures.py ~/papers/pymrm-gallery/"Xu1989-methane-steam-reforming-kinetics-AIChEJ35-88.pdf"

Requires poppler's ``pdftoppm`` on PATH, plus numpy, scipy and pillow.

Method, and why it is not the Duncan-Toor method
------------------------------------------------
Only the experimental markers are wanted. The smooth curves on Figures 2 and 3
are the authors' own model predictions, and extracting them would turn the
page's validation into model-against-model.

At 600 dpi the curves are ~10 px thick and the markers are multi-stroke glyphs
(x, +, open square, open circle, open and filled triangle) only ~20 px across,
drawn on top of the curves. That is too similar in scale for the morphological
opening that isolates Duncan & Toor's solid markers. Instead a marker is treated
as a *local excess* of ink over what the curve alone puts into the same
neighbourhood:

    D = ink area inside a WIN x WIN box (uniform filter)
    B = max over six orientations of the grey-scale opening of D with a long
        LINE element -- the largest value at that point explainable by an
        extended, locally straight structure
    E = D - B -> bumps no straight line can explain == markers

Taking the maximum over orientations is the essential step. Near the origin all
four curves are steep; a vertical line element reproduces a vertical curve, so B
tracks D there and no false peaks appear. A single horizontal element would flag
the entire steep section as marker.

The automatic pass is not trusted on its own. Curvature of the predicted curves
still produces excess peaks with no glyph under them, and a strong glyph
produces a second peak on the curve beside it. ACCEPT_F2 / ACCEPT_F3 below are
the candidate indices that survived a visual audit of every candidate at 600
dpi, and EXTRA_F2 / EXTRA_F3 are two glyph pairs the split still merged. Those
lists are the reason this script is deterministic rather than a heuristic: re-run
it and you get the audited set back. Change WIN, LINE, THRESH or SEP and the
candidate numbering changes, so the audit would have to be redone.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import (uniform_filter, grey_opening, label, maximum_filter)

DPI = 600
PAGE_FIGURES = 2           # journal page 89 is PDF page 2

# Crop boxes in the 600 dpi page render: (col_lo, row_lo, col_hi, row_hi)
CROP = {"fig2": (2380, 260, 4560, 1800), "fig3": (2380, 4020, 4560, 5450)}
# Region searched for markers, and text blocks to keep out of it, in crop pixels
SEARCH = {"fig2": (200, 1870, 150, 1418), "fig3": (200, 1870, 150, 1382)}
EXCLUDE = {"fig2": [(340, 990, 180, 300)], "fig3": []}      # fig2 legend block

WIN, LINE, THRESH, SEP, RAD = 25, 121, 120.0, 21, 15
ANGLES = (0, 30, 60, 90, 120, 150)

# Axis calibration from the printed tick marks: (value, pixel) pairs.
TICKS = {
    "fig2": {"x": [(0.10, 571.5), (0.20, 975.5), (0.30, 1376.5), (0.40, 1775.5)],
             "y": [(0.150, 208), (0.125, 408), (0.100, 613),
                   (0.075, 813), (0.050, 1021), (0.025, 1216)]},
    "fig3": {"x": [(0.10, 567.5), (0.20, 972.0), (0.30, 1375.0), (0.40, 1772.0)],
             "y": [(0.125, 123.5), (0.100, 374.5), (0.075, 627.5),
                   (0.050, 885.5), (0.025, 1139.0)]},
}

# Audited candidate indices (see the module docstring): first pass over the
# raw candidates, then a second pass over the split result, which renumbers.
ACCEPT = {
    "fig2": [3, 4, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 24, 25,
             26, 29, 30, 31, 32, 34, 35, 37, 38, 39, 40],
    "fig3": [3, 4, 5, 6, 7, 8, 10, 11, 13, 15, 17, 18, 19, 20, 22, 24, 28, 29,
             30, 35, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48],
}
REFINED_KEEP = {
    "fig2": [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 21, 22, 23,
             27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
    "fig3": [0, 1, 2, 3, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 18, 19, 22, 23, 24,
             29, 31, 34, 35, 36, 37, 38, 39, 40, 41, 42],
}
# Glyph pairs the split still merged, added by hand from the 600 dpi image.
EXTRA = {"fig2": [(1576, 700)],   # 'x' just above the open square at (1567, 723)
         "fig3": [(771, 846)]}    # filled triangle below the square at (772, 822)
SEP_FIG = {"fig2": 21, "fig3": 17}
CENTROID_RAD = {"fig2": RAD, "fig3": 11}

# Curve labels, assigned by conversion band. Gaps between bands are 0.0062 or
# wider, against a digitisation error of 0.0006.
BANDS = {"fig2": [(0.000, 0.022, 773), (0.022, 0.052, 798),
                  (0.052, 0.100, 823), (0.100, 0.200, 848)],
         "fig3": [(0.000, 0.021, 773), (0.021, 0.046, 798),
                  (0.046, 0.085, 823), (0.085, 0.200, 848)]}
QUANTITY = {"fig2": "x_CH4", "fig3": "x_CO2"}


# --------------------------------------------------------------------- imaging
def render(pdf: Path, out_dir: Path) -> dict[str, np.ndarray]:
    stem = out_dir / "page"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(PAGE_FIGURES),
                    "-l", str(PAGE_FIGURES), "-png", str(pdf), str(stem)],
                   check=True)
    hits = sorted(out_dir.glob("page-*.png"))
    if not hits:
        sys.exit("pdftoppm produced no output")
    page = Image.open(hits[0]).convert("L")
    if page.size != (4800, 6600):
        sys.exit(f"unexpected page size {page.size}; the crop boxes assume "
                 f"(4800, 6600) at {DPI} dpi")
    return {name: np.array(page.crop(box)) for name, box in CROP.items()}


def line_footprint(length: int, angle_deg: float) -> np.ndarray:
    n = length if length % 2 else length + 1
    fp = np.zeros((n, n), bool)
    c = n // 2
    t = np.linspace(-c, c, 4 * n)
    a = np.deg2rad(angle_deg)
    fp[np.clip(np.rint(c - t * np.sin(a)).astype(int), 0, n - 1),
       np.clip(np.rint(c + t * np.cos(a)).astype(int), 0, n - 1)] = True
    return fp


def excess_map(grey: np.ndarray) -> np.ndarray:
    ink = (grey < 128).astype(np.float32)
    dens = uniform_filter(ink, WIN, mode="constant") * WIN * WIN
    bg = None
    for a in ANGLES:
        op = grey_opening(dens, footprint=line_footprint(LINE, a))
        bg = op if bg is None else np.maximum(bg, op)
    return dens - bg


# ------------------------------------------------------------------ detection
def candidates(excess: np.ndarray, name: str) -> list[tuple[float, float]]:
    c0, c1, r0, r1 = SEARCH[name]
    mask = np.zeros(excess.shape, bool)
    mask[r0:r1, c0:c1] = excess[r0:r1, c0:c1] >= THRESH
    lab, n = label(mask)
    out = []
    for i in range(1, n + 1):
        sel = lab == i
        if sel.sum() < 8:
            continue
        rr, cc = np.nonzero(sel)
        w = excess[rr, cc]
        x, y = float((cc * w).sum() / w.sum()), float((rr * w).sum() / w.sum())
        if any(a <= x <= b and c <= y <= d for a, b, c, d in EXCLUDE[name]):
            continue
        out.append((x, y))
    return sorted(out)


def centroid(excess: np.ndarray, x: float, y: float, rad: int = RAD):
    sub = excess[int(y) - rad:int(y) + rad + 1, int(x) - rad:int(x) + rad + 1].copy()
    sub[sub < 0] = 0
    if sub.sum() == 0:
        return float(x), float(y)
    rr, cc = np.mgrid[int(y) - rad:int(y) + rad + 1, int(x) - rad:int(x) + rad + 1]
    return float((cc * sub).sum() / sub.sum()), float((rr * sub).sum() / sub.sum())


def split(excess: np.ndarray, x: float, y: float, sep: int, half=40, thresh=110):
    """Local maxima at least `sep` apart, to separate stacked glyph pairs."""
    r0, c0 = int(y) - half, int(x) - half
    sub = excess[r0:r0 + 2 * half + 1, c0:c0 + 2 * half + 1]
    peaks = np.argwhere((sub == maximum_filter(sub, size=sep)) & (sub >= thresh))
    kept = []
    for r, c in peaks:
        if all((r - a) ** 2 + (c - b) ** 2 > sep ** 2 for a, b in kept):
            kept.append((r, c))
    return [(c0 + c, r0 + r) for r, c in kept]


def markers(excess: np.ndarray, name: str) -> list[tuple[float, float]]:
    cand = candidates(excess, name)
    sep = SEP_FIG[name]
    pts: list[tuple[float, float]] = []
    for i in ACCEPT[name]:
        # The visual audit was carried out on candidate centres quantised to
        # 0.1 px, and `split` centres its search window on int(centre), so the
        # quantisation has to be reproduced or a window can shift by one pixel
        # and admit an extra curvature peak that was never audited.
        cx, cy = round(cand[i][0], 1), round(cand[i][1], 1)
        for mx, my in split(excess, cx, cy, sep=sep):
            x, y = centroid(excess, mx, my)
            if all((x - a) ** 2 + (y - b) ** 2 > (sep * 0.8) ** 2 for a, b in pts):
                pts.append((x, y))
    pts.sort()
    if len(pts) <= max(REFINED_KEEP[name]):
        sys.exit(f"{name}: split produced {len(pts)} points, fewer than the "
                 f"audited indices expect; the detector settings have changed "
                 f"and the visual audit must be redone")
    pts = [pts[i] for i in REFINED_KEEP[name]]
    for x, y in EXTRA[name]:
        pts.append(centroid(excess, x, y, rad=CENTROID_RAD[name]))
    return sorted(pts)


# ---------------------------------------------------------------- calibration
def fit(pairs):
    v = np.array([p[0] for p in pairs], float)
    px = np.array([p[1] for p in pairs], float)
    (m, c), *_ = np.linalg.lstsq(np.vstack([v, np.ones_like(v)]).T, px, rcond=None)
    return m, c, float(np.abs(np.vstack([v, np.ones_like(v)]).T @ [m, c] - px).max())


def to_data(pts, name):
    mx, cx, ex = fit(TICKS[name]["x"])
    my, cy, ey = fit(TICKS[name]["y"])
    print(f"  {name}: x tick residual {ex:.1f} px = {ex / mx:.5f}; "
          f"y tick residual {ey:.1f} px = {abs(ey / my):.5f}")
    rows = []
    for x, y in pts:
        wf, val = (x - cx) / mx, (y - cy) / my
        T = next(t for lo, hi, t in BANDS[name] if lo <= val < hi)
        rows.append((T, wf, QUANTITY[name], val))
    return sorted(rows)


# ------------------------------------------------------------------ self-check
def self_check(rows):
    ok = True
    f2 = [r for r in rows if r[2] == "x_CH4"]
    f3 = [r for r in rows if r[2] == "x_CO2"]
    for tag, sub in (("x_CH4", f2), ("x_CO2", f3)):
        for T in sorted({r[0] for r in sub}):
            v = [r[3] for r in sub if r[0] == T]
            drops = sum(1 for i in range(1, len(v)) if v[i] < v[i - 1] - 0.004)
            if drops:
                ok = False
            print(f"  {tag} {T} K: {len(v)} points, {drops} drop(s) beyond scatter")
    used, dwf, viol = set(), [], 0
    for T, wf, _, xc in f2:
        cand = [(abs(wf - w), i) for i, (t, w, _, _) in enumerate(f3)
                if t == T and i not in used]
        if not cand:
            continue
        d, i = min(cand)
        if d > 0.004:
            continue
        used.add(i)
        dwf.append(d)
        if f3[i][3] > xc + 0.0012:
            viol += 1
    print(f"  cross-figure pairing: {len(dwf)}/{len(f2)} paired, "
          f"max |dW/F| {max(dwf):.4f}, mean {np.mean(dwf):.4f}")
    print(f"  carbon closure violations (x_CO2 > x_CH4 + 0.0012): {viol}")
    if len(dwf) != len(f2) or viol:
        ok = False
    return ok


HEADER = """# Xu & Froment (1989), AIChE J 35(1) 88-96, Figures 2 and 3.
# Experimental markers only; the plotted curves are the authors'
# own model and were deliberately not extracted.
# p_t = 10 bar, H2O/CH4 = 3.0, H2/CH4 = 1.25 (Table 1).
# W_F is W/F0_CH4 in g_cat h / mol_CH4.
# Regenerate with extract_figures.py; provenance in the .meta.yaml sidecar.
temperature_K,W_F,quantity,value,figure"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pdf", type=Path)
    ap.add_argument("-o", "--out", type=Path,
                    default=Path(__file__).with_name("data")
                    / "xu-froment-1989-conversion.csv")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        crops = render(args.pdf, Path(tmp))

    rows = []
    print("axis calibration")
    for name in ("fig2", "fig3"):
        e = excess_map(crops[name])
        pts = markers(e, name)
        print(f"  {name}: {len(pts)} markers")
        rows += to_data(pts, name)

    print("self-checks")
    ok = self_check(rows)

    lines = [HEADER] + [f"{T},{wf:.4f},{q},{v:.5f},{2 if q == 'x_CH4' else 3}"
                        for T, wf, q, v in rows]
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {args.out} ({len(rows)} rows)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
