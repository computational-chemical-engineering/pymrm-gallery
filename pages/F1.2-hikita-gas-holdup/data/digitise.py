"""Digitised readings of Hikita et al. (1980), Figures 2 and 3, and their calibration.

The readings below are marker and line centres in page pixels of the PDF pages
rendered at 600 dpi (`pdftoppm -r 600 -gray -png`), twice the native 300 ppi of the
scan; pixel errors in the sidecars are quoted in native pixels (half of these). Each
axis is calibrated by a least-squares affine map from the printed ticks on both
frames (log10 of the tick value against the pixel column and row), which absorbs
the small skew of the scan. Run from the page directory:

    python data/digitise.py

It writes hikita-1980-fig2.csv and hikita-1980-fig3.csv and prints the tick-fit
residuals that the sidecars quote. No page image is stored in the repository.
"""
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent

# ---- Figure 2 (journal p. 63): eps against uG (cm/s), log-log, four gases in water.
# Ticks: (value, px_x, px_y). Abscissa ticks on the bottom frame (row 4583) and the top
# frame (row 3986); ordinate ticks on the left frame (column 513) and the right (1780).
F2_X = [(v, x, 4583) for v, x in ((4, 583.5), (6, 691.5), (8, 765.5), (10, 824.5), (20, 1002.5),
                                  (40, 1183.5), (60, 1288.5), (80, 1362.5), (100, 1421.5), (200, 1600.5))] + \
       [(v, x, 3986) for v, x in ((4, 582.5), (6, 689.5), (8, 759.5), (10, 821.5), (20, 999.5),
                                  (40, 1181.5), (60, 1285.5), (80, 1359.5), (100, 1417.5), (200, 1595.5))]
F2_Y = [(v, 513, y) for v, y in ((0.6, 3989), (0.4, 4094.5), (0.2, 4274.5), (0.1, 4451.5),
                                 (0.08, 4511.5), (0.06, 4583))] + \
       [(v, 1780, y) for v, y in ((0.4, 4085.5), (0.2, 4264.5), (0.1, 4446.5), (0.08, 4502.5))]
F2_MARKERS = {
    "air_water": [(598, 4475), (757, 4385), (847, 4322), (957, 4260), (1040, 4212), (1162, 4150)],
    "H2N2_11_water": [(783, 4394), (847, 4340), (945, 4285), (1020, 4235)],
    "H2_water": [(805, 4389), (875, 4350), (990, 4288), (1057, 4250)],
    "H2N2_51_water": [(791, 4414), (900, 4357), (1007, 4277), (1072, 4240), (1132, 4206)],
}
F2_SHAPE = {"air_water": "open circle", "H2N2_11_water": "open diamond", "H2_water": "open square",
            "H2N2_51_water": "filled diamond"}
# The two drawn lines: centre row of the line at a pixel column.
F2_LINES = {"upper": [(575, 4492), (640, 4452.5), (700, 4416.5), (1200, 4123.5), (1220, 4111.5)],
            "lower": [(720, 4448.5), (740, 4436.5), (1100, 4225.5), (1180, 4178.5)]}

# ---- Figure 3 (journal p. 64): (eps)obs against (eps)cal of eq. 8, log-log, fitted set.
F3_X = [(v, x, 1861) for v, x in ((0.08, 516.5), (0.09, 578.5), (0.1, 634.5), (0.2, 992.5), (0.3, 1199.5),
                                  (0.4, 1351.5), (0.5, 1465.5), (0.6, 1559.5), (0.7, 1639.0), (0.8, 1712.5),
                                  (0.9, 1776.5))] + \
       [(v, x, 660.5) for v, x in ((0.08, 512.5), (0.09, 574.5), (0.1, 629.5), (0.2, 986.5), (0.3, 1196.5),
                                   (0.4, 1347.5), (0.5, 1461.5), (0.6, 1557.5), (0.7, 1637.5), (0.8, 1709.5),
                                   (0.9, 1773.5))]
F3_Y = [(v, 441, y) for v, y in ((0.5, 761.5), (0.4, 876.5), (0.3, 1028.5), (0.2, 1239.5), (0.1, 1598.5),
                                 (0.09, 1652.5), (0.08, 1714.0), (0.07, 1784.5))] + \
       [(v, 1830, y) for v, y in ((0.5, 757.5), (0.4, 872.5), (0.3, 1024.5), (0.2, 1234.0), (0.1, 1595.5),
                                  (0.09, 1650.0), (0.08, 1713.5), (0.07, 1780.5))]
# Marker centres by shape; the system follows from the Fig. 3 legend (Table 2's
# fig3_marker column). "unassigned": the shape could not be decided (open up against
# open down triangles in the dense band; filled markers half hidden).
F3_MARKERS = {
    "filled square": [(581, 1633), (1065, 1130)],
    "cross": [(717, 1565), (808, 1441), (938, 1307), (1182, 1070)],
    "open square": [(714, 1481), (825, 1378), (1013, 1188), (1157, 1075), (1227, 992)],
    "filled diamond": [(748, 1520), (890, 1392), (995, 1247), (1083, 1183), (1160, 1115)],
    "plus": [(763, 1498), (860, 1385), (982, 1273), (1045, 1183)],
    "open circle": [(790, 1483), (867, 1327), (943, 1268), (1025, 1217), (1120, 1120), (1257, 992)],
    "filled triangle": [(760, 1443)],
    "filled circle": [(789, 1415), (877, 1298), (1258, 957)],
    "open diamond": [(878, 1363), (968, 1250), (1075, 1163)],
}
F3_UNASSIGNED = [(735, 1473, "filled; half hidden"), (767, 1475, "open triangle"), (822, 1410, "open triangle"),
                 (843, 1357, "filled"), (903, 1340, "open triangle"), (922, 1338, "open triangle"),
                 (925, 1265, "filled"), (945, 1247, "filled"), (980, 1222, "open triangle"),
                 (985, 1195, "filled"), (1098, 1115, "open triangle"), (1122, 1088, "filled"),
                 (1218, 1012, "filled"), (1312, 908, "filled; probably the n-butanol down triangle"),
                 (1368, 871, "open down triangle; probably methanol")]
# The drawn line of Fig. 3: centre row of the line (a run of about 10 px at 600 dpi) at
# columns where no marker touches it.
F3_LINE = [(620, 1608.5), (660, 1566.5), (700, 1528.5), (1300, 922.5), (1340, 880), (1380, 839.5), (1420, 799.5)]


def fit_axis(ticks):
    """Least-squares log10(value) = a + b px_x + c px_y; returns the coefficients and the
    residuals in native pixels along the axis."""
    v, x, y = (np.array(c, float) for c in zip(*ticks))
    A = np.column_stack([np.ones_like(x), x, y])
    coef, *_ = np.linalg.lstsq(A, np.log10(v), rcond=None)
    along = np.hypot(coef[1], coef[2])            # decades per 600-dpi pixel
    res = (A @ coef - np.log10(v)) / along / 2    # native pixels
    return coef, res, 1 / along


def to_value(coef, x, y):
    return 10 ** (coef[0] + coef[1] * x + coef[2] * y)


def system_of(shape):
    import csv
    with open(HERE / "hikita-1980-table2-properties.csv", newline="") as fh:
        rows = list(csv.DictReader(r for r in fh if not r.startswith("#")))
    return {r["fig3_marker"]: r["key"] for r in rows if r["fig3_marker"]}[shape]


def main():
    cx2, rx2, dx2 = fit_axis(F2_X)
    cy2, ry2, dy2 = fit_axis(F2_Y)
    cx3, rx3, dx3 = fit_axis(F3_X)
    cy3, ry3, dy3 = fit_axis(F3_Y)
    for tag, r, d in (("fig2 x", rx2, dx2), ("fig2 y", ry2, dy2), ("fig3 x", rx3, dx3), ("fig3 y", ry3, dy3)):
        print(f"{tag}: {len(r)} ticks, {d:.1f} px/decade at 600 dpi, max |residual| {np.abs(r).max():.2f} native px")
    head2 = ("# Hikita et al. (1980) Figure 2 (p. 63): gas holdup against superficial gas velocity, four gases in "
             "water at 20 C, D 10 cm. Digitised by data/digitise.py; see the sidecar.\n")
    with open(HERE / "hikita-1980-fig2.csv", "w") as fh:
        fh.write(head2 + "kind,series,shape,uG_cm_s,eps,px_x,px_y\n")
        for key, pts in F2_MARKERS.items():
            for x, y in pts:
                fh.write(f"marker,{key},{F2_SHAPE[key]},{to_value(cx2, x, y):.4f},{to_value(cy2, x, y):.5f},{x},{y}\n")
        for key, pts in F2_LINES.items():
            for x, y in pts:
                fh.write(f"line,{key},,{to_value(cx2, x, y):.4f},{to_value(cy2, x, y):.5f},{x},{y}\n")
    head3 = ("# Hikita et al. (1980) Figure 3 (p. 64): observed against calculated (eq. 8) gas holdup, the fitted "
             "systems. Digitised by data/digitise.py; see the sidecar.\n")
    with open(HERE / "hikita-1980-fig3.csv", "w") as fh:
        fh.write(head3 + "kind,system,shape,eps_cal,eps_obs,px_x,px_y,note\n")
        for shape, pts in F3_MARKERS.items():
            for x, y in pts:
                fh.write(f"marker,{system_of(shape)},{shape},{to_value(cx3, x, y):.5f},{to_value(cy3, x, y):.5f},"
                         f"{x},{y},\n")
        for x, y, note in F3_UNASSIGNED:
            fh.write(f"marker,unassigned,,{to_value(cx3, x, y):.5f},{to_value(cy3, x, y):.5f},{x},{y},{note}\n")
        for x, y in F3_LINE:
            fh.write(f"line,y=x,,{to_value(cx3, x, y):.5f},{to_value(cy3, x, y):.5f},{x},{y},\n")


if __name__ == "__main__":
    main()
