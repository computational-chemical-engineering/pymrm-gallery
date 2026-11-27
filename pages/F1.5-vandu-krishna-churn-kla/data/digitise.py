"""Extraction of the plotted values of Vandu and Krishna (2004), Figures 2, 5, 7, 9 and 11
(vector paths) and Figures 4, 8 and 10 (600 ppi bilevel bitmaps).

The PDF is not in the repository. Pass a local copy of the publisher PDF:

    python data/digitise.py /path/to/Vandu2004.pdf

Vector figures. Each page content stream is decompressed with `mutool show -b` and its
path operators (m, l, c, re and the paint operators) are parsed. Marker centres
(D1): a circle is the centre of the box of its Bezier control points, a square the
centre of its `re` rectangle, a triangle the mean of its three vertices, a diamond the
mean of its four. Figures 7 and 11 join their markers with polylines: a polyline
vertex is the data point. The axes are calibrated piecewise linearly between adjacent
labelled ticks (in log10 on a log axis), with the tick positions read from the axis
paths of the same stream; the script asserts that every tick it uses is a vertex of
its axis path. Coordinates are PDF points (y up), exact to the 0.001 pt the stream
prints.

Raster figures. `pdfimages -png` extracts the stencil bitmaps at their native 600 ppi
(no resampling). Marker centres are found by matching each marker shape against a
template cut from the same panel's legend (circle and square: the box centre;
triangle: the vertex centroid, h/6 from the box centre toward the base). Each axis is
calibrated by a least-squares line through its labelled tick marks, whose pixel
positions are listed below; the residuals are printed and quoted in the sidecars.

Outputs (next to this script): vandu-krishna-2004-vector.csv,
vandu-krishna-2004-fig2a-curves.csv and the raster CSVs. No page image is stored.
"""
import csv
import math
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent

# ---------------------------------------------------------------------------------------
# vector figures: page content objects, ticks (position, value), plot boxes, legends
# ---------------------------------------------------------------------------------------
VECTOR = {
    "fig5": dict(obj=33, log=(True, True),
                 xt=[(142.252, 1e2), (190.396, 1e3), (238.54, 1e4), (286.786, 1e5)],
                 yt=[(620.683, 1e-2), (692.899, 1e-1), (765.217, 1.0)],
                 legend=(149.0, 246.5, 629.0, 693.5),
                 shapes={"tri_down": "water", "circle": "tetradecane", "square": "paraffin",
                         "tri_up": "tellus"}),
    "fig9": dict(obj=66, log=(True, True),
                 xt=[(144.088, 1e2), (192.232, 1e3), (240.376, 1e4), (288.622, 1e5)],
                 yt=[(614.815, 1e-2), (687.134, 1e-1), (759.451, 1.0)],
                 legend=(151.5, 268.0, 620.0, 663.5),
                 shapes={"tri_down": "slurry25", "circle": "tellus"}),
    "fig7": dict(obj=52, log=(False, False),
                 xt=[(150.540, 0.0), (179.406, 0.05), (208.272, 0.10), (237.240, 0.15), (266.107, 0.20),
                     (295.074, 0.25)],
                 yt=[(620.382, 0.0), (644.454, 0.1), (668.627, 0.2), (692.599, 0.3), (716.771, 0.4),
                     (740.844, 0.5), (764.915, 0.6)],
                 legend=(154.0, 274.5, 624.0, 657.5),
                 polylines={3: "water", 4: "tetradecane"}),      # vertex counts: Table 4's eps_s levels
    "fig11": dict(obj=83, log=(False, True),
                  xt=[(153.785, 0.0), (198.972, 0.2), (244.056, 0.4), (289.242, 0.6)],
                  yt=[(619.786, 1e-2), (692.104, 1e-1), (764.320, 1.0)],
                  legend=(196.0, 261.5, 628.0, 661.5),
                  polylines={4: "water", 2: "tellus"}),          # Table 4's column diameters per liquid
    "fig2b": dict(obj=16, log=(False, False),
                  xt=[(372.469, 0.0), (406.799, 1.0), (441.131, 2.0), (475.562, 3.0), (509.895, 4.0)],
                  yt=[(167.061, 0.0), (191.204, 0.2), (215.247, 0.4), (239.289, 0.6), (263.330, 0.8),
                      (287.372, 1.0), (311.414, 1.2)],
                  box=(372.0, 530.0, 160.0, 320.0), legend=(380.0, 457.0, 246.5, 306.5),
                  shapes={"circle": "tetradecane", "square": "paraffin", "tri_down": "water",
                          "tri_up": "tellus", "diamond": "ethanol"}),
    "fig2a": dict(obj=16, log=(False, False),
                  xt=[(180.697, 0.0), (216.76, 2.0), (252.822, 4.0), (288.883, 6.0), (325.048, 8.0)],
                  yt=[(167.061, 0.0), (195.484, 0.2), (223.703, 0.4), (252.023, 0.6), (280.343, 0.8),
                      (308.664, 1.0)],
                  box=(180.0, 326.0, 160.0, 320.0), legend=(243.0, 320.5, 180.5, 240.5)),
}
LINE_MIN_PT = 5.0          # a stroked subpath shorter than this is a tick or a marker detail


def stream(pdf, obj):
    """The decompressed content stream of one PDF object."""
    return subprocess.run(["mutool", "show", "-b", str(pdf), str(obj)], capture_output=True,
                          check=True).stdout.decode("latin-1")


def paths(text):
    """(paint operator, [subpath, ...]) with each subpath a list of (op, x, y, extra)."""
    text = re.sub(r"\bBT\b.*?\bET\b", " ", text, flags=re.S)       # text objects carry no paths
    out, sub, subs, num = [], [], [], []
    for tok in text.split():
        if re.fullmatch(r"-?\d*\.?\d+", tok):
            num.append(float(tok))
            continue
        if tok == "m":
            if sub:
                subs.append(sub)
            sub = [("m", num[-2], num[-1], None)]
        elif tok == "l":
            sub.append(("l", num[-2], num[-1], None))
        elif tok == "c":
            sub.append(("c", num[-2], num[-1], tuple(num[-6:])))
        elif tok == "re":
            if sub:
                subs.append(sub)
            subs.append([("re", num[-4], num[-3], tuple(num[-4:]))])
            sub = []
        elif tok in ("S", "s", "f", "f*", "F", "B", "B*", "b", "b*", "n"):
            if sub:
                subs.append(sub)
            if subs:
                out.append((tok, subs))
            sub, subs = [], []
        num = []
    return out


def box_of(sub):
    xs, ys = [], []
    for op, x, y, e in sub:
        if op == "re":
            x0, y0, w, h = e
            xs += [x0, x0 + w]
            ys += [y0, y0 + h]
        else:
            xs.append(x)
            ys.append(y)
            if op == "c":
                xs += [e[0], e[2]]
                ys += [e[1], e[3]]
    return min(xs), max(xs), min(ys), max(ys)


def inside(b, region):
    return region[0] <= b[0] and b[1] <= region[1] and region[2] <= b[2] and b[3] <= region[3]


def marker(sub):
    """(shape, centre x, centre y) of a filled marker subpath, D1 rules."""
    ops = "".join(op[0] for op, *_ in sub)
    b = box_of(sub)
    if ops == "r":
        return "square", (b[0] + b[1]) / 2, (b[2] + b[3]) / 2
    if "c" in ops:
        return "circle", (b[0] + b[1]) / 2, (b[2] + b[3]) / 2
    v = [(x, y) for op, x, y, _ in sub]
    if len(v) > 1 and v[-1] == v[0]:
        v = v[:-1]
    if len(v) == 3:
        ys = sorted(y for _, y in v)
        apex_low = abs(ys[0] - ys[1]) > abs(ys[1] - ys[2])        # PDF y is up
        cx, cy = sum(x for x, _ in v) / 3, sum(y for _, y in v) / 3
        return ("tri_down" if apex_low else "tri_up"), cx, cy
    if len(v) == 4:
        return "diamond", sum(x for x, _ in v) / 4, sum(y for _, y in v) / 4
    return "other", (b[0] + b[1]) / 2, (b[2] + b[3]) / 2


def calibrate(p, ticks, log):
    """Piecewise linear between adjacent ticks (log10 on a log axis); the end intervals extend."""
    for i in range(len(ticks) - 1):
        (a, va), (b, vb) = ticks[i], ticks[i + 1]
        if p <= b or i == len(ticks) - 2:
            if log:
                la, lb = math.log10(va), math.log10(vb)
                return 10.0 ** (la + (p - a) * (lb - la) / (b - a))
            return va + (p - a) * (vb - va) / (b - a)


def check_ticks(subpaths, spec):
    """Every tick position used must be a vertex of a stroked axis path of the stream."""
    verts = [(x, y) for sub in subpaths for _, x, y, _ in sub]
    for axis, ticks in (("x", spec["xt"]), ("y", spec["yt"])):
        for pos, _ in ticks:
            k = 0 if axis == "x" else 1
            assert any(abs(v[k] - pos) < 5e-4 for v in verts), (axis, pos)


def vector(pdf):
    rows, curves, notes = [], [], []
    for fig, spec in VECTOR.items():
        allp = paths(stream(pdf, spec["obj"]))
        lx, ly = spec["log"]
        x0, x1 = spec["xt"][0][0], spec["xt"][-1][0]
        y0, y1 = spec["yt"][0][0], spec["yt"][-1][0]
        region = spec.get("box", (x0 - 6, x1 + 12, y0 - 6, y1 + 12))
        stroked = [s for op, subs in allp if op in ("S", "s") for s in subs if s[0][0] in ("m",)]
        check_ticks([s for s in stroked if inside(box_of(s), (x0 - 1, x1 + 1, y0 - 4, y1 + 4))], spec)
        cal = lambda x, y: (calibrate(x, spec["xt"], lx), calibrate(y, spec["yt"], ly))
        if "shapes" in spec:
            for op, subs in allp:
                if op != "f*":
                    continue
                for s in subs:
                    b = box_of(s)
                    if not inside(b, region) or inside(b, spec["legend"]) or b[1] - b[0] < 3.0:
                        continue
                    shape, cx, cy = marker(s)
                    if shape not in spec["shapes"]:
                        continue
                    vx, vy = cal(cx, cy)
                    rows.append((fig, spec["shapes"][shape], shape, cx, cy, vx, vy))
        lines = []
        for s in stroked:
            v = [(x, y) for op, x, y, _ in s if op in ("m", "l")]
            if len(v) != len(s) or len(v) < 2 or v[0] == v[-1]:
                continue
            length = sum(math.dist(a, b) for a, b in zip(v, v[1:]))
            b = box_of(s)
            if length < LINE_MIN_PT or not inside(b, region) or inside(b, spec["legend"]):
                continue
            if abs(b[2] - b[3]) < 1e-9 and abs(b[2] - y0) < 1e-3 or abs(b[0] - b[1]) < 1e-9 and abs(b[0] - x0) < 1e-3:
                continue                                           # an axis line itself
            lines.append(v)
        if fig in ("fig5", "fig9"):
            assert len(lines) == 1, (fig, len(lines))
            for end, (px, py) in (("left", lines[0][0]), ("right", lines[0][-1])):
                vx, vy = cal(px, py)
                rows.append((fig, "line_" + end, "line end", px, py, vx, vy))
            notes.append(f"{fig}: drawn line with {len(lines[0])} vertices")
        elif fig in ("fig7", "fig11"):
            assert sorted(len(v) for v in lines) == sorted(spec["polylines"]), (fig, [len(v) for v in lines])
            for v in lines:
                name = spec["polylines"][len(v)]
                for px, py in v:
                    vx, vy = cal(px, py)
                    rows.append((fig, name, "polyline vertex", px, py, vx, vy))
            # D1 evidence: distance of each vertex to the nearest triangle centroid and box centre
            cents = []
            for op, subs in allp:
                for s in subs:
                    if op == "f*" and not inside(box_of(s), spec["legend"]) and marker(s)[0].startswith("tri"):
                        b = box_of(s)
                        cents.append((marker(s)[1:], ((b[0] + b[1]) / 2, (b[2] + b[3]) / 2)))
            near = [(p, min(cents, key=lambda c: math.dist(p, c[0]))) for v in lines for p in v if cents]
            near = [(p, c) for p, c in near if math.dist(p, c[0]) < 2.0]     # vertices drawn with a triangle
            dc = [math.dist(p, c[0]) for p, c in near]
            db = [math.dist(p, c[1]) for p, c in near]
            if dc:
                notes.append(f"{fig}: vertex to triangle centroid {min(dc):.3f} to {max(dc):.3f} pt, "
                             f"to box centre {min(db):.3f} to {max(db):.3f} pt")
        elif fig == "fig2a":
            joined = []
            for v in lines:                     # a curve split over two paths continues at its last vertex
                for c in joined:
                    if c[-1] == v[0]:
                        c.extend(v[1:])
                        break
                else:
                    joined.append(list(v))
            assert len(joined) == 5, len(joined)
            for i, v in enumerate(sorted(joined, key=lambda c: c[len(c) // 3][1])):
                for px, py in v:
                    vx, vy = cal(px, py)
                    curves.append((f"curve{i + 1}", px, py, vx, vy))
    return rows, curves, notes


def write_csv(name, header, comment, rows, fmt):
    with open(HERE / name, "w", newline="") as fh:
        fh.write(f"# {comment}\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([f(x) for f, x in zip(fmt, r)])


def main(pdf):
    rows, curves, notes = vector(pdf)
    g = lambda x: f"{x:.15g}"
    s = str
    write_csv("vandu-krishna-2004-vector.csv", ["figure", "series", "shape", "x_pt", "y_pt", "value_x", "value_y"],
              "Vandu and Krishna (2004), plotted values of Figs. 2b, 5, 7, 9, 11 from the PDF vector paths "
              "(data/digitise.py; see the sidecar)", rows, (s, s, s, g, g, g, g))
    write_csv("vandu-krishna-2004-fig2a-curves.csv", ["curve", "x_pt", "y_pt", "t_s", "ratio"],
              "Vandu and Krishna (2004), Fig. 2a drawn response curves (polyline vertices; data/digitise.py; "
              "curves numbered from the slowest)", curves, (s, g, g, g, g))
    print(f"vector: {len(rows)} rows, Fig. 2a: {len(curves)} vertices")
    for n in notes:
        print(n)
    if "--vector-only" not in sys.argv:
        import digitise_raster
        digitise_raster.main(pdf)


if __name__ == "__main__":
    main(sys.argv[1])
