"""Regenerate the digitised figure data of page C2.5 from the source PDF.

The PDF is not in the repository and must not be (AGENTS.md). Point this at
your copy:

    python extract_figures.py ~/papers/pymrm-gallery/VandenBussche1996-methanol-synthesis-WGS-kinetics-JCatal161-1.pdf

Figures 5 to 8 of Vanden Bussche and Froment (1996) are embedded as 1-bit
scans (307 to 369 ppi). They are extracted at native resolution with
poppler's ``pdfimages`` (no resampling), calibrated on their own tick marks,
and read by column scans of the black pixels. Requires ``pdfimages`` on PATH,
numpy and pillow.

Outputs (in ``data/``):

- ``vandenbussche-1996-fig5b.csv``: the temperature trace of Fig. 5b.
- ``vandenbussche-1996-fig8.csv``: the simulated curve of Fig. 8.
- ``vandenbussche-1996-figure-readings.csv``: plateaus of Figs. 5, 6 and 7,
  inlet readings of Fig. 7a, the 15-bar CO2 minimum of Fig. 6c and the lowest
  of Klier's markers at 2 % CO2 in Fig. 8.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

# Embedded image index (pdfimages order) and its pixel size, checked on load.
IMAGES = {5: (1534, 515), 6: (1408, 972), 7: (1525, 526), 8: (1734, 1044)}

# Axis calibration: (pixel, value) pairs of tick marks found on each panel
# (outer or inner ticks, centre of the tick stroke) and the frame columns that
# carry the ends of the z/L axis (0 and 1 are the only labelled x ticks of
# Figs. 5 to 7; Fig. 5's unlabelled middle ticks are not at 0.5 and are not used).
CAL = {
    "5a": {"x": [(97.5, 0.0), (704.5, 1.0)], "y": [(443.5, 0.0), (271.0, 2.0), (98.0, 4.0)]},
    "5b": {"x": [(866.5, 0.0), (1525.5, 1.0)], "y": [(396.0, 500.0), (257.0, 525.0), (114.0, 550.0)]},
    "6a": {"x": [(86.5, 0.0), (654.5, 1.0)], "y": [(409.5, 0.0), (238.0, 1.2), (67.0, 2.4)]},
    "6b": {"x": [(817.5, 0.0), (1398.5, 1.0)], "y": [(410.5, 490.0), (310.5, 510.0)]},
    "6c": {"x": [(460.5, 0.0), (1008.5, 1.0)], "y": [(905.5, 2.0), (743.5, 2.5), (580.0, 3.0)]},
    "7a": {"x": [(104.5, 0.0), (727.0, 1.0)], "y": [(442.5, 180.0), (292.5, 230.0), (142.5, 280.0)]},
    "7b": {"x": [(867.5, 0.0), (1517.0, 1.0)], "y": [(443.0, 0.0), (304.5, 1.0), (166.5, 2.0)]},
    "8": {"x": [(150.5, 0.0), (411.5, 5.0), (669.5, 10.0), (929.5, 15.0), (1190.5, 20.0),
                (1450.5, 25.0), (1710.5, 30.0)],
          "y": [(914.5, 0.0), (801.0, 10.0), (689.0, 20.0), (576.0, 30.0), (464.0, 40.0),
                (352.0, 50.0), (240.0, 60.0), (128.0, 70.0), (14.5, 80.0)]},
}
PANEL_IMAGE = {"5a": 5, "5b": 5, "6a": 6, "6b": 6, "6c": 6, "7a": 7, "7b": 7, "8": 8}
# Row window inside each panel's frame (frame lines and outside excluded).
ROWS = {"5a": (62, 440), "5b": (62, 451), "6a": (60, 405), "6b": (66, 407), "6c": (585, 902),
        "7a": (57, 439), "7b": (58, 439), "8": (18, 909)}


def fit(pairs):
    """Piecewise-linear value(pixel) through the tick pairs (the scans are not exactly
    linear: Fig. 5b's two tick intervals differ by 4 px), extrapolated with the end
    segments, and the largest residual of a single straight line through the ticks."""
    p, v = np.array(sorted(pairs), dtype=float).T
    slope, icpt = np.polyfit(p, v, 1)

    def value(q):
        q = np.asarray(q, dtype=float)
        lo = v[0] + (q - p[0]) * (v[1] - v[0]) / (p[1] - p[0])
        hi = v[-1] + (q - p[-1]) * (v[-1] - v[-2]) / (p[-1] - p[-2])
        return np.where(q < p[0], lo, np.where(q > p[-1], hi, np.interp(q, p, v)))
    return value, float(np.max(np.abs(slope * p + icpt - v)))


def runs(col):
    """Runs of black pixels in a column: list of (start, end) indices, inclusive."""
    idx = np.flatnonzero(col)
    if idx.size == 0:
        return []
    cut = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate([[idx[0]], idx[cut + 1]])
    ends = np.concatenate([idx[cut], [idx[-1]]])
    return list(zip(starts, ends))


class Panel:
    def __init__(self, images, key):
        self.key = key
        self.a = images[PANEL_IMAGE[key]]
        self.xv, self.xres = fit(CAL[key]["x"])
        self.yv, self.yres = fit(CAL[key]["y"])
        xp = np.array(CAL[key]["x"])[:, 0]
        self.x_of = lambda v: np.interp(v, np.array(CAL[key]["x"])[:, 1], xp)  # value -> pixel
        self.r0, self.r1 = ROWS[key]
        yp, yv = np.array(CAL[key]["y"], dtype=float).T
        self.unit = abs((yv[-1] - yv[0]) / (yp[-1] - yp[0]))   # value per pixel

    def column_runs(self, c):
        return [(s + self.r0, e + self.r0) for s, e in runs(self.a[self.r0:self.r1, c])]

    def band(self, lo, hi, window, max_run=8):
        """Run centres inside the value window [window] over x in [lo, hi]."""
        c0, c1 = int(np.ceil(self.x_of(lo))), int(np.floor(self.x_of(hi)))
        vals = []
        for c in range(c0, c1 + 1):
            for s, e in self.column_runs(c):
                v = float(self.yv((s + e) / 2.0))
                if e - s + 1 <= max_run and min(window) <= v <= max(window):
                    vals.append(v)
        return np.array(vals)


def load_images(pdf, tmp):
    subprocess.run(["pdfimages", "-png", str(pdf), str(Path(tmp) / "img")], check=True)
    out = {}
    for i, size in IMAGES.items():
        im = Image.open(Path(tmp) / f"img-{i:03d}.png")
        if im.size != size:
            sys.exit(f"image {i}: size {im.size}, expected {size}: not this PDF")
        out[i] = np.array(im.convert("L")) < 128
    return out


def plateau(rows, panel, fig, curve, quantity, lo, hi, window, unit):
    v = panel.band(lo, hi, window)
    if v.size < 10:
        sys.exit(f"{fig} {curve}: only {v.size} pixel runs in the plateau window")
    rows.append((fig, curve, quantity, lo, hi, float(np.mean(v)), float(np.std(v)),
                 int(v.size), panel.unit, unit))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pdf", type=Path)
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        images = load_images(args.pdf, tmp)
    P = {k: Panel(images, k) for k in CAL}
    for k, p in P.items():
        print(f"panel {k}: {p.unit:.5g} per pixel, tick residuals x {p.xres:.3g}, y {p.yres:.3g}")

    # Fig. 5b: one solid curve; the topmost run in each column inside the frame
    p = P["5b"]
    trace = []
    for c in range(int(p.x_of(0.0)) + 4, int(p.x_of(1.0)) - 3):
        rr = p.column_runs(c)
        if rr:
            s, e = rr[0]
            trace.append((float(p.xv(c)), float(p.yv((s + e) / 2.0)), int(e - s + 1)))
    with open(DATA / "vandenbussche-1996-fig5b.csv", "w") as f:
        f.write("zL,T_K,run_px\n")
        for z, T, w in trace:
            f.write(f"{z:.5f},{T:.3f},{w}\n")

    # Fig. 8: the simulated line, thin runs only (the markers are wider than 8 px)
    p = P["8"]
    curve = []
    marker_cols = [(p.x_of(2.0), 16), (p.x_of(20.0), 16)]
    for c in range(int(p.x_of(0.6)), int(p.x_of(20.9)) + 1):
        if any(abs(c - mc) <= w for mc, w in marker_cols):
            continue
        rr = [(s, e) for s, e in p.column_runs(c) if e - s + 1 <= 8 and s > p.r0 + 2]
        rr = [(s, e) for s, e in rr if (e + s) / 2 > 120]     # legend text above 71 %
        if len(rr) == 1:
            s, e = rr[0]
            curve.append((float(p.xv(c)), float(p.yv((s + e) / 2.0)), int(e - s + 1)))
    with open(DATA / "vandenbussche-1996-fig8.csv", "w") as f:
        f.write("x_CO2_pct,X_MeOH_pct,run_px\n")
        for x, X, w in curve:
            f.write(f"{x:.4f},{X:.3f},{w}\n")

    rows = []
    # Fig. 5: plateaus (z/L 0.45 to 0.75, clear of the labels)
    for curve_, win in (("CO", (2.6, 3.1)), ("MeOH", (2.25, 2.5)), ("CO2", (1.95, 2.2)),
                        ("H2O", (0.9, 1.2))):
        plateau(rows, P["5a"], "5a", curve_, "y_molpct", 0.45, 0.75, win, "mol%")
    plateau(rows, P["5b"], "5b", "base", "T_K", 0.45, 0.95, (540, 560), "K")
    # Fig. 6: plateaus at 15, 30, 50 bar (z/L 0.5 to 0.75)
    for curve_, win in (("50 bar", (2.2, 2.5)), ("30 bar", (1.6, 1.9)), ("15 bar", (0.95, 1.2))):
        plateau(rows, P["6a"], "6a", curve_, "y_MeOH_molpct", 0.5, 0.75, win, "mol%")
    for curve_, win in (("50 bar", (544, 556)), ("30 bar", (530, 541)), ("15 bar", (512, 523))):
        plateau(rows, P["6b"], "6b", curve_, "T_K", 0.5, 0.75, win, "K")
    # Fig. 7: plateaus (z/L 0.8 to 0.97 in a, 0.5 to 0.8 in b). The inlet ends of Fig. 7a are
    # not read: near z/L = 0 the five curves rise steeply and their runs merge.
    for curve_, win in (("280 C", (300, 312)), ("240 C", (285, 295)), ("220 C", (277, 284)),
                        ("200 C", (268, 276))):
        plateau(rows, P["7a"], "7a", curve_, "T_C", 0.8, 0.97, win, "C")
    for curve_, lo, win in (("200 C", 0.74, (2.65, 2.9)), ("220 C", 0.5, (2.25, 2.45)),
                            ("240 C", 0.5, (1.85, 2.05)), ("280 C", 0.5, (1.1, 1.3))):
        plateau(rows, P["7b"], "7b", curve_, "y_MeOH_molpct", lo, 0.80, win, "mol%")
    # Fig. 6c: the minimum of the solid 15-bar CO2 curve (topmost run for z/L 0.06 to 0.35)
    p = P["6c"]
    best = None
    for c in range(int(p.x_of(0.06)), int(p.x_of(0.35)) + 1):
        rr = p.column_runs(c)
        if rr:
            s, e = rr[0]
            v = float(p.yv((s + e) / 2.0))
            if best is None or v < best[1]:
                best = (float(p.xv(c)), v)
    rows.append(("6c", "15 bar", "y_CO2_min_molpct", 0.06, 0.35, best[1], p.unit, 1, p.unit,
                 "mol%"))
    rows.append(("6c", "15 bar", "zL_at_y_CO2_min", 0.06, 0.35, best[0], 0.0, 1,
                 1.0 / (CAL["6c"]["x"][1][0] - CAL["6c"]["x"][0][0]), "-"))
    # Fig. 8: Klier's markers at 2 % CO2 form one blob; its lowest marker centre is the
    # blob's bottom less half the height of an isolated marker (the one at 10 % CO2)
    p = P["8"]
    a = p.a

    def blob_rows(xc, r_lo, r_hi):
        c = int(round(p.x_of(xc)))
        sub = a[r_lo:r_hi, c - 3:c + 4].all(axis=1)          # 7 px wide and black: marker body
        idx = np.flatnonzero(sub) + r_lo
        return idx.min(), idx.max()

    t10, b10 = blob_rows(10.0, 400, 520)
    t2, b2 = blob_rows(2.0, 100, 200)
    half = (b10 - t10) / 2.0
    rows.append(("8", "Klier 2 %", "X_lowest_marker_pct", 2.0, 2.0,
                 float(p.yv(b2 - half)), p.unit, 1, p.unit, "%"))
    rows.append(("8", "Klier 2 %", "X_highest_marker_pct", 2.0, 2.0,
                 float(p.yv(t2 + half)), p.unit, 1, p.unit, "%"))
    with open(DATA / "vandenbussche-1996-figure-readings.csv", "w") as f:
        f.write("figure,curve,quantity,zL_lo,zL_hi,value,sd,n_px,value_per_px,unit\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]:.4f},{r[6]:.4f},{r[7]},"
                    f"{r[8]:.5f},{r[9]}\n")
            print(r[:3], round(r[5], 4), round(r[6], 4), r[7])
    print(f"fig5b: {len(trace)} points; fig8: {len(curve)} points")


if __name__ == "__main__":
    main()
