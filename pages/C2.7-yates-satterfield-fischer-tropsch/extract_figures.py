#!/usr/bin/env python3
"""Regenerate the digitised CSVs of page C2.7 from the source PDF.

The PDF is not in the repository and must not be (AGENTS.md). Point this at your
own lawfully obtained copy:

    python extract_figures.py <your copy>.pdf

Requires poppler's ``pdftoppm`` on PATH, plus numpy, scipy and pillow.

Method
------
Pages 170 and 172 (PDF pages 3 and 5) are rendered at 300 dpi, the native
resolution of the scan. Ink is every pixel darker than 128 of 255.

Calibration. Each figure frame is its axis limits (0 to 0.6 by 0 to 1.2 MPa in
Figures 1 and 2, 0 to 1.5 by 0 to 1.5 in Figure 8). The scan is turned by
about 0.3 degree, so each axis is an affine function of both pixel
coordinates, fitted by least squares to the centroids of the printed tick
marks (each tick sits at a known value of its own axis, every 0.1 on every
axis). The rms residual of the fit is the calibration
error, written to the sidecar.

Markers. Inside the frame (4 px clear of it) the ink is split into connected
components. A component of one marker's size is a marker, located at the
centroid of its ink. A larger component is a cluster: it is matched by a union
of copies of an isolated marker's ink (the template), placed at integer pixel
offsets, choosing the copies that minimise the number of unmatched pixels, and
each copy's centre is the template centroid at that offset. Figure 1 has one
cluster of three circles, Figure 2 three pairs of triangles; the 23rd point of
Figure 2 (Table V footnote: 23 at 240 C) is not visible, so it lies under one
of the 22 visible triangles (recorded as missing in the sidecar).

Figure 8 carries the 45 degree line through the markers: a binary opening
with a disk of radius 4 px removes the line (3 px wide) and the open 230 C
circles, and keeps the filled triangles and squares. Only markers that stand
alone after the opening are recorded with a symbol; the clump below 0.3
mmol/(min g) is recorded as one unassigned row. The two squares that stand alone
touch the line, so their centres come from a second opening with a 10 px square
element, which keeps a square whole and removes the line.
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

DPI = 300
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

# frame lines (pixel columns left, right; rows top, bottom) in the 300 dpi render
FRAME = {1: (1615.5, 2179.0, 173.5, 836.5), 2: (1614.0, 2177.5, 1058.5, 1724.0),
         8: (1600.5, 2171.5, 1061.0, 1726.5)}
LIMITS = {1: (0.6, 1.2), 2: (0.6, 1.2), 8: (1.5, 1.5)}
PDF_PAGE = {1: 3, 2: 3, 8: 5}
SINGLE = {1: (160, 220), 2: (110, 145)}     # ink area of one marker (px)


def render(pdf, page, tmp):
    """Render one PDF page at DPI and return the ink mask."""
    stem = Path(tmp) / f"p{page}"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(page), "-l", str(page),
                    "-png", str(pdf), str(stem)], check=True)
    png = sorted(Path(tmp).glob(f"p{page}*.png"))[0]
    return np.asarray(Image.open(png).convert("L")) < 128


def _runs(profile, thresh):
    """(start, stop) of the runs where ``profile >= thresh``."""
    on = np.concatenate([[False], profile >= thresh, [False]])
    edges = np.flatnonzero(np.diff(on.astype(int)))
    return list(zip(edges[::2], edges[1::2]))


def ticks(ink, fig):
    """Tick centroids (u, v, axis, value) of one figure (every 0.1 on every axis).

    The ticks touch the frame, so they are found in 11 px bands just inside it,
    as runs of at most 4 px of ink in the band's column (or row) sums, 20 px
    clear of the corners (where the diagonal of Figure 8 crosses the bands).
    """
    x0, x1, y0, y1 = FRAME[fig]
    xmax, ymax = LIMITS[fig]
    step = 0.1
    out = []
    X0, X1, Y0, Y1 = (int(round(v)) for v in FRAME[fig])
    for r0 in (Y1 - 14, Y0 + 3):                          # bottom and top bands
        band = ink[r0:r0 + 11, X0 + 20:X1 - 20]
        for a, b in _runs(band.sum(0), 8):
            if b - a <= 4:
                cols = np.arange(a, b)
                u = X0 + 20 + float((band[:, a:b].sum(0) * cols).sum() / band[:, a:b].sum())
                val = round((u - x0) / (x1 - x0) * xmax / step) * step
                out.append((u, r0 + 5.0, "x", val))
    for c0 in (X0 + 3, X1 - 14):                          # left and right bands
        band = ink[Y0 + 20:Y1 - 20, c0:c0 + 11]
        for a, b in _runs(band.sum(1), 8):
            if b - a <= 4:
                rows = np.arange(a, b)
                v = Y0 + 20 + float((band[a:b].sum(1) * rows).sum() / band[a:b].sum())
                val = round((y1 - v) / (y1 - y0) * ymax / step) * step
                out.append((c0 + 5.0, v, "y", val))
    return out


def calibrate(tk):
    """Affine maps pixel (u, v) -> axis values, fitted to the ticks; rms residuals."""
    maps, rms = {}, {}
    for ax in ("x", "y"):
        rows = [(u, v, val) for u, v, a, val in tk if a == ax]
        A = np.array([[1.0, u, v] for u, v, _ in rows])
        b = np.array([val for _, _, val in rows])
        c, *_ = np.linalg.lstsq(A, b, rcond=None)
        maps[ax] = c
        rms[ax] = float(np.sqrt(np.mean((A @ c - b) ** 2)))
    return maps, rms, len(tk)


def to_axes(maps, u, v):
    return (maps["x"] @ [1.0, u, v], maps["y"] @ [1.0, u, v])


def components(ink, fig, pad=4):
    x0, x1, y0, y1 = (int(round(v)) for v in FRAME[fig])
    sub = ink[y0 + pad:y1 - pad, x0 + pad:x1 - pad]
    lab, n = ndimage.label(sub)
    comps = []
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        m = lab[sl] == k
        comps.append((m, sl[1].start + x0 + pad, sl[0].start + y0 + pad))
    return comps


def fit_union(mask, tmpl, k):
    """Offsets of k template copies whose union best matches ``mask``."""
    th, tw = tmpl.shape
    H, W = mask.shape[0] + 2 * th, mask.shape[1] + 2 * tw
    big = np.zeros((H, W), bool)
    big[th:th + mask.shape[0], tw:tw + mask.shape[1]] = mask
    cands = []
    for i in range(H - th + 1):
        for j in range(W - tw + 1):
            if (big[i:i + th, j:j + tw] & tmpl).sum() >= 0.8 * tmpl.sum():
                c = np.zeros((H, W), bool)
                c[i:i + th, j:j + tw] = tmpl
                cands.append(((i - th, j - tw), c.ravel()))
    flat = big.ravel()
    best = None
    for combo in itertools.combinations(range(len(cands)), k):
        u = np.logical_or.reduce([cands[c][1] for c in combo])
        cost = int((u ^ flat).sum())
        if best is None or cost < best[0]:
            best = (cost, [cands[c][0] for c in combo])
    return best


def design(ink, fig):
    """Markers of Figure 1 or 2: rows (u, v, cluster, n_in_cluster)."""
    lo, hi = SINGLE[fig]
    comps = [c for c in components(ink, fig) if c[0].sum() >= 60]
    singles = [c for c in comps if lo <= c[0].sum() <= hi and max(c[0].shape) <= 17]
    tmpl = singles[0][0]
    ty, tx = np.where(tmpl)
    tcy, tcx = ty.mean(), tx.mean()
    rows, cid = [], 0
    for m, cx, cy in comps:
        if any(m is s[0] for s in singles):
            ys, xs = np.where(m)
            rows.append((cx + xs.mean(), cy + ys.mean(), "", 1))
        elif m.sum() > hi and m.shape[0] > 4 and m.shape[1] > 4:
            k = int(round(m.sum() / (0.85 * tmpl.sum())))
            k = min(max(k, 2), 3)
            cost, offs = fit_union(m, tmpl, k)
            cid += 1
            for (i, j) in offs:
                rows.append((cx + j + tcx, cy + i + tcy, f"c{fig}.{cid}", k))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    pdf = Path(ap.parse_args().pdf)
    with tempfile.TemporaryDirectory() as tmp:
        pages = {p: render(pdf, p, tmp) for p in (3, 5)}
    lines = ["T_C,P_CO_MPa,P_H2_MPa,u_px,v_px,cluster,n_in_cluster"]
    for fig, T in ((1, 220), (2, 240)):
        ink = pages[PDF_PAGE[fig]]
        maps, rms, nt = calibrate(ticks(ink, fig))
        rows = design(ink, fig)
        print(f"Figure {fig}: {len(rows)} markers; {nt} ticks, rms x {rms['x']:.4f} y {rms['y']:.4f} MPa")
        for u, v, cl, k in sorted(rows, key=lambda r: (r[0], r[1])):
            pco, ph2 = to_axes(maps, u, v)
            lines.append(f"{T},{pco:.4f},{ph2:.4f},{u:.1f},{v:.1f},{cl},{k}")
    (DATA / "yates-satterfield-1991-fig1-fig2-design.csv").write_text("\n".join(lines) + "\n")
    # Figure 8: filled markers that stand alone after removing the line
    ink = pages[PDF_PAGE[8]]
    maps, rms, nt = calibrate(ticks(ink, 8))
    print(f"Figure 8: {nt} ticks, rms x {rms['x']:.4f} y {rms['y']:.4f}")
    x0, x1, y0, y1 = (int(round(v)) for v in FRAME[8])
    sub = ink[y0 + 4:y1 - 4, x0 + 4:x1 - 4].copy()
    sub[1080 - y0 - 4:1230 - y0 - 4, 1640 - x0 - 4:1900 - x0 - 4] = False   # legend
    yy, xx = np.mgrid[-4:5, -4:5]
    op = ndimage.binary_opening(sub, structure=(xx ** 2 + yy ** 2) <= 16)
    lab, n = ndimage.label(op)
    # a 10 px square element keeps the squares whole and removes what the disk left of
    # the line where it touches them; triangles (14 px) cannot hold it
    sq = ndimage.binary_opening(sub, structure=np.ones((10, 10), bool))
    out = ["symbol,R_exp,R_pred,u_px,v_px,area_px,status"]
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        ys, xs = np.where(lab[sl] == k)
        if SYMBOL.get(k) == "square":
            ys, xs = np.where(sq[sl] & (lab[sl] == k))
        u, v = xs.mean() + sl[1].start + x0 + 4, ys.mean() + sl[0].start + y0 + 4
        area, h, w = len(xs), np.ptp(ys) + 1, np.ptp(xs) + 1
        rexp, rpred = to_axes(maps, u, v)
        out.append(f"{SYMBOL.get(k, '?')},{rexp:.4f},{rpred:.4f},{u:.1f},{v:.1f},{area},{STATUS.get(k, '?')}")
        print(k, f"{rexp:.4f} {rpred:.4f} area {area} {w}x{h}")
    (DATA / "yates-satterfield-1991-fig8-parity.csv").write_text("\n".join(out) + "\n")


# symbol and status of each Figure 8 component (labels in scan order, top to bottom),
# assigned by eye on an 8x zoom of the render
# (inverted triangle 240 C, square 220 C). "pair": two triangles merged, one row at the
# pair's centroid; "clump": squares and triangles below 0.3 merged, unassigned.
SYMBOL = {**{k: "triangle" for k in range(1, 16)}, 16: "square", 17: "square", 18: "mixed"}
STATUS = {**{k: "isolated" for k in (1, 2, 3, 4, 5, 8, 9, 11, 13, 14, 15, 16, 17)},
          **{k: "pair" for k in (6, 7, 10, 12)}, 18: "clump"}

if __name__ == "__main__":
    main()
