"""Raster part of data/digitise.py: Figures 4c, 8c and 10f of Vandu and Krishna (2004).

Called by digitise.py with the PDF path. `pdfimages -png` extracts the stencil bitmaps at
their native 600 ppi (ink is True). For each panel:

- the axes are the longest horizontal and vertical ink runs of the panel; each tick is the
  run of ink columns (rows) just inside the axis, refined within +-9 px of its nominal place
  between two labelled ticks; the calibration is a least-squares line through every tick
  found, and its largest residual (px) is printed;
- every marker of a series is a peak of the fraction of the series' legend marker (the
  template, cut from the same panel's legend) that is covered by ink; a peak counts when
  that fraction is at least COVER and no other template scores higher within half a marker;
- centres (D1): circle and square the box centre of the template; triangle the vertex
  centroid, h/6 from the box centre toward the base;
- each thick horizontal line (the paper's drawn plateau) is the set of rows with a large
  ink count, at least LINE_ROWS rows thick; its value is the mean row.
"""
import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.signal import fftconvolve

HERE = Path(__file__).parent
COVER = 0.86          # fraction of a template's ink that must be present at a marker (default)
TICK_W, TICK_OFF = 5, 4.0   # px: widest tick, largest offset of a tick from its nominal place
LINE_ROWS, LINE_GAP = 6, 60  # px: thinnest drawn plateau line, widest gap (a marker) inside one

# page, image index on that page, panel crop (x0, x1, y0, y1) native px, ticks
# (first nominal px, value, last nominal px, value, steps), series: (shape, legend centre px);
# umin: markers below this U are not read (overlapping clusters); cover: per-series COVER where
# the plotted marker differs from the legend's; bands: the value band a series is searched in, set
# by eye, so that a half-hidden marker is not given to a series plotted far away
PANELS = {
    "fig4c": dict(page=4, index=1, crop=(2280, 3430, 0, 1170),
                  x=(2497, 0.0, 3376.5, 0.4, 20), y=(1006, 0.0, 82, 0.8, 20),
                  series={"water": ("tri_down", (2807, 140)), "tetradecane": ("circle", (2807, 202)),
                          "paraffin": ("square", (2807, 263)), "tellus": ("tri_up", (2807, 325))},
                  legend=(2712, 3205, 70, 392), umin=0.03,
                  cover={"paraffin": 0.75},       # the plotted squares are a little smaller than the legend's
                  bands={"water": (0.35, 0.8), "tetradecane": (0.35, 0.8), "paraffin": (0.1, 0.35),
                         "tellus": (0.0, 0.1)}),
    "fig8c": dict(page=6, index=0, crop=(2800, 3921, 0, 1235),
                  x=(2915, 0.0, 3871, 0.4, 20), y=(1073, 0.0, 69, 0.30, 30),
                  series={"tellus": ("circle", (3115, 527)), "slurry25": ("tri_down", (3115, 617))},
                  legend=(3013, 3752, 450, 815), cover={}, bands={"tellus": (0.0, 0.1), "slurry25": (0.15, 0.3)}),
    "fig10f": dict(page=7, index=0, crop=(2800, 3941, 1300, 2681),
                   x=(2875.5, 0.0, 3854, 0.4, 20), y=(2509, 0.0, 1481.5, 0.10, 25),
                   series={"tellus_D0.10": ("circle", (3249, 2229)), "tellus_D0.38": ("tri_up", (3249, 2304))},
                   legend=None, cover={"tellus_D0.10": 0.80}, bands={}),
}


def image(pdf, page, index, tmp):
    subprocess.run(["pdfimages", "-png", "-f", str(page), "-l", str(page), str(pdf), f"{tmp}/p{page}"],
                   check=True)
    return np.array(Image.open(f"{tmp}/p{page}-{index:03d}.png")).astype(bool)


def longest_run(v):
    best = cur = 0
    for b in v:
        cur = cur + 1 if b else 0
        best = max(best, cur)
    return best


def axes(a, x0, y0):
    """Row of the x axis (lowest of the longest runs) and column of the y axis (leftmost)."""
    rr = np.array([longest_run(r) for r in a])
    cc = np.array([longest_run(c) for c in a.T])
    return int(np.where(rr >= 0.95 * rr.max())[0].max()) + y0, int(np.where(cc >= 0.95 * cc.max())[0].min()) + x0


def ticks(a, spec, along, axis_pos, x0, y0):
    """Refined tick centres (px, value) of one axis; nominal places that fail are skipped."""
    p0, v0, p1, v1, n = spec
    out, miss = [], []
    for i in range(n + 1):
        p = p0 + (p1 - p0) * i / n
        lo = int(round(p)) - 9
        if along == "x":
            band = a[axis_pos - y0 - 12:axis_pos - y0 - 3, lo - x0:lo - x0 + 19].all(axis=0)
        else:
            band = a[lo - y0:lo - y0 + 19, axis_pos - x0 + 4:axis_pos - x0 + 12].all(axis=1)
        lab, k = ndi.label(band)
        c = lo + np.where(lab == 1)[0].mean() if k == 1 else np.nan
        if k != 1 or (lab == 1).sum() > TICK_W or abs(c - p) > TICK_OFF:     # a marker on the tick
            miss.append(round(v0 + (v1 - v0) * i / n, 4))
            continue
        out.append((c, v0 + (v1 - v0) * i / n))
    return out, miss


def fit(t):
    p, v = np.array(t).T
    b, c = np.polyfit(p, v, 1)
    return (lambda q: b * q + c), float(np.max(np.abs((v - (b * p + c)) / b)))


def template(a, centre, shape):
    """The legend marker containing ``centre``: its ink patch and its D1 centre offset."""
    lab, _ = ndi.label(a)
    cx, cy = centre
    win = lab[cy - 30:cy + 31, cx - 30:cx + 31]
    ids, cnt = np.unique(win[win > 0], return_counts=True)
    k = ids[np.argmax(cnt)]
    sl = ndi.find_objects((lab == k).astype(int))[0]
    patch = (lab[sl] == k)
    h, w = patch.shape
    oy = (h - 1) / 2 + {"tri_up": h / 6, "tri_down": -h / 6}.get(shape, 0.0)   # image y is down
    return patch.astype(float), ((w - 1) / 2, oy)


def markers(a, temps, region, limits):
    """Peaks of the template cover, one series per peak (the highest-scoring template)."""
    best = None
    for name, (patch, off) in temps.items():
        cover = fftconvolve(a.astype(float), patch[::-1, ::-1], mode="valid") / patch.sum()
        best = (cover, name, patch.shape, off) if best is None else best
        temps[name] = (patch, off, cover)
    found = []
    for name, (patch, off, cover) in temps.items():
        h, w = patch.shape
        cmin, (ylo, yhi) = limits[name]
        peak = (cover == ndi.maximum_filter(cover, size=(h // 2, w // 2))) & (cover >= cmin)
        rows_ok = (np.arange(cover.shape[0]) + off[1] >= ylo) & (np.arange(cover.shape[0]) + off[1] <= yhi)
        peak &= rows_ok[:, None]
        for yy, xx in zip(*np.nonzero(peak)):
            cy, cx = yy + off[1], xx + off[0]
            if region is not None and region[0] <= cx <= region[1] and region[2] <= cy <= region[3]:
                continue
            rival = max((c2[min(max(int(cy - o2[1]), 0), c2.shape[0] - 1),
                            min(max(int(cx - o2[0]), 0), c2.shape[1] - 1)]
                         for n2, (p2, o2, c2) in temps.items() if n2 != name), default=0.0)
            found.append((name, cx, cy, float(cover[yy, xx]), float(rival)))
    keep = []
    for f in sorted(found, key=lambda f: -f[3]):          # one marker per place
        if all(np.hypot(f[1] - k[1], f[2] - k[2]) > 18 for k in keep):
            keep.append(f)
    return keep


def lines(a, x_axis_row, y_axis_col, x0, y0, min_frac=0.25):
    """Thick horizontal lines (>= 3 rows with a long run), excluding the x axis."""
    w = a.shape[1]
    rows = [r for r in range(a.shape[0]) if r + y0 < x_axis_row - 6 and a[r].sum() >= min_frac * w]
    groups, out = [], []
    for r in rows:
        if groups and r == groups[-1][-1] + 1:
            groups[-1].append(r)
        else:
            groups.append([r])
    for g in groups:
        if len(g) < LINE_ROWS:
            continue                        # a legend frame or an axis, not a drawn plateau
        cols = np.nonzero(a[g].all(axis=0))[0]
        segs = np.split(cols, np.nonzero(np.diff(cols) > LINE_GAP)[0] + 1)
        seg = max(segs, key=len)
        out.append((float(np.mean(g)) + y0, int(seg[0]) + x0, int(seg[-1]) + x0, len(g)))
    return out


def panel(pdf, name, spec, tmp):
    a_full = image(pdf, spec["page"], spec["index"], tmp)
    x0, x1, y0, y1 = spec["crop"]
    a = a_full[y0:y1, x0:x1]
    xrow, ycol = axes(a, x0, y0)
    tx, mx = ticks(a, spec["x"], "x", xrow, x0, y0)
    ty, my = ticks(a, spec["y"], "y", ycol, x0, y0)
    fx, rx = fit(tx)
    fy, ry = fit(ty)
    temps = {}
    for s, (shape, c) in spec["series"].items():
        temps[s] = template(a_full, c, shape)
    legend = spec["legend"]
    if legend is None:                      # the legend frame: the box around the legend markers
        cs = [c for _, c in spec["series"].values()]
        legend = (min(c[0] for c in cs) - 60, x1, min(c[1] for c in cs) - 60, max(c[1] for c in cs) + 60)
    sub = a_full[y0:y1, x0:x1]
    inv = lambda v: (v - fy(0.0)) / (fy(1.0) - fy(0.0)) - y0          # value -> row in the crop
    limits = {}
    for s in spec["series"]:
        lo, hi = spec["bands"].get(s, (-np.inf, np.inf))
        r = sorted((inv(lo), inv(hi))) if np.isfinite(lo) else (-np.inf, np.inf)
        limits[s] = (spec["cover"].get(s, COVER), r)
    found = markers(sub, temps, (legend[0] - x0, legend[1] - x0, legend[2] - y0, legend[3] - y0), limits)
    rows = []
    for s, cx, cy, cov, rival in sorted(found, key=lambda f: (f[0], f[1])):
        X, Y = cx + x0, cy + y0
        if X <= ycol + 8 or Y >= xrow - 8 or fx(X) < spec.get("umin", 0.0):
            continue
        rows.append(("marker", s, spec["series"][s][0], fx(X), fy(Y), X, Y, cov, rival))
    for yy, xa, xb, th in lines(sub, xrow, ycol, x0, y0):
        # the series whose markers above 0.08 m/s have the nearest median value draws this line
        med = {s: np.median([r[4] for r in rows if r[0] == "marker" and r[1] == s and r[3] > 0.08])
               for s in spec["series"]}
        s = min(med, key=lambda k: abs(med[k] - fy(yy)))
        rows.append(("line", s, f"line {th} px thick", fx(xa), fy(yy), xa, yy, 0.0, 0.0))
        rows.append(("line", s, f"line {th} px thick", fx(xb), fy(yy), xb, yy, 0.0, 0.0))
    print(f"{name}: axis row {xrow}, column {ycol}; ticks x {len(tx)} (missing {mx}), y {len(ty)} "
          f"(missing {my}); tick residual x {rx:.2f} px, y {ry:.2f} px; px per unit x {1 / (fx(1) - fx(0)):.1f}, "
          f"y {1 / (fy(1) - fy(0)):.1f}; markers " + ", ".join(
              f"{s} {sum(r[1] == s for r in rows)}" for s in spec["series"]) +
          f"; lines {sum(r[0] == 'line' for r in rows) // 2}")
    return rows


def main(pdf, only=None):
    with tempfile.TemporaryDirectory() as tmp:
        for name, spec in PANELS.items():
            if only and name not in only:
                continue
            rows = panel(pdf, name, spec, tmp)
            with open(HERE / f"vandu-krishna-2004-{name}.csv", "w", newline="") as fh:
                fh.write(f"# Vandu and Krishna (2004) Fig. {name[3:]}: kLa/eps_G against U, read from the native "
                         "600 ppi bitmap by data/digitise_raster.py; see the sidecar\n")
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(["kind", "series", "shape", "U_m_s", "P_1_s", "px_x", "px_y", "cover", "rival"])
                for r in rows:
                    w.writerow([r[0], r[1], r[2], f"{r[3]:.5f}", f"{r[4]:.5f}", f"{r[5]:.1f}", f"{r[6]:.1f}",
                                f"{r[7]:.3f}", f"{r[8]:.3f}"])


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2:] or None)
