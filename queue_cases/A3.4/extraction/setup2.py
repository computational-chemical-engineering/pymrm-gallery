"""Figure 2 (page 1382): image, calibration, care/allow masks."""
import os

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

S = os.environ.get("A34_RENDERS", "renders")
PATH = f"{S}/fig2-full.png"

# calibration recorded with the 2026-07-30 extraction and confirmed by the three
# printed curves; NOT re-derived here (the maintainer took it as established).
COL0, PXDEC_X = 442.0, 405.2          # Re = 1 at COL0
ROW0, PXDEC_Y = 1225.25, 337.67       # Sh = 1 at ROW0
SC = 0.6
FRAME = dict(left=440.5, right=1657.5, top=211.5, bottom=1223.5)

# text blocks to exclude, (row0, row1, col0, col1), read off a 100-px grid overlay
TEXT_BOXES = [
    (262, 645, 505, 1215),     # the 8-line legend, glyph column included
    (695, 765, 665, 895),      # "Sc = 0 6"
    (795, 870, 510, 960),      # "Single spheres"
    (1090, 1170, 780, 1620),   # "Petrovic and Thodos  Eqn(11)"
]


def load():
    return 1.0 - np.asarray(Image.open(PATH).convert("L"), float) / 255.0


def re_of_col(c):
    return 10.0 ** ((np.asarray(c, float) - COL0) / PXDEC_X)


def col_of_re(re):
    return COL0 + PXDEC_X * np.log10(np.asarray(re, float))


def sh_of_row(r):
    return 10.0 ** ((ROW0 - np.asarray(r, float)) / PXDEC_Y)


def row_of_sh(sh):
    return ROW0 - PXDEC_Y * np.log10(np.asarray(sh, float))


def eq12(re):
    return 2.0 + 1.1 * SC ** (1 / 3) * np.asarray(re, float) ** 0.6


def eq9(re):
    return 2.0 + 0.6 * SC ** (1 / 3) * np.asarray(re, float) ** 0.5


def eq11(re, eps_b=0.40):
    return (0.357 / eps_b) * np.asarray(re, float) ** (-0.359) * \
        np.asarray(re, float) * SC ** (1 / 3)


def curve_mask(shape, halfwidth, which=(eq12, eq9, eq11)):
    """True on pixels within `halfwidth` px of the named printed curves."""
    m = np.zeros(shape, bool)
    cols = np.arange(int(FRAME["left"]), int(FRAME["right"]) + 1)
    re = re_of_col(cols)
    for f in which:
        rows = row_of_sh(f(re))
        good = np.isfinite(rows) & (rows > FRAME["top"]) & (rows < FRAME["bottom"])
        m[np.round(rows[good]).astype(int), cols[good]] = True
    return ndi.binary_dilation(m, ndi.generate_binary_structure(2, 2),
                               iterations=halfwidth)


def tick_mask(ink, reach=30):
    """Frame lines and ticks: ink connected to the frame, within `reach` px of
    it.  Clipping to `reach` keeps the erasure from following a chain of merged
    markers away from the frame, which a plain connected-component erasure does
    (the eq. 12 curve runs into the right-hand frame and drags the whole dense
    band with it)."""
    h, w = ink.shape
    b = ink > 0.5
    lab, _ = ndi.label(b, ndi.generate_binary_structure(2, 2))
    seed = np.zeros((h, w), bool)
    for r in (int(FRAME["top"]), int(FRAME["bottom"])):
        seed[r - 2:r + 3, int(FRAME["left"]):int(FRAME["right"])] = True
    for c in (int(FRAME["left"]), int(FRAME["right"])):
        seed[int(FRAME["top"]):int(FRAME["bottom"]), c - 2:c + 3] = True
    ids = np.unique(lab[seed & b])
    touch = np.isin(lab, ids[ids > 0])
    zone = np.zeros((h, w), bool)
    zone[:int(FRAME["top"]) + reach, :] = True
    zone[int(FRAME["bottom"]) - reach:, :] = True
    zone[:, :int(FRAME["left"]) + reach] = True
    zone[:, int(FRAME["right"]) - reach:] = True
    return touch & zone


def clean2(ink, long_len=37, short_len=23, curve_hw=12, dash_hw=9, reach=30,
           n_ang=12):
    """Everything on this figure that is not a marker, in three pieces.

    1. Locally straight ink over `long_len` px at any orientation: the eq. (12)
       and eq. (11) curves and the two annotation leaders.  A glyph is 18-26 px
       across and cannot contain a 37-px straight run, so glyphs survive.
    2. Ink connected to the frame and lying within `reach` px of it: the frame
       and every tick.  Ticks are only ~26 px long, so (1) does not see them.
    3. A band of half-width `dash_hw` about the COMPUTED position of eq. (9).
       The Ranz-Marshall curve is dashed, and its dashes are 12-20 px of compact
       ink - indistinguishable from a marker by any shape test, which is why
       they must be removed by position.  Plus straight ink over `short_len`
       within `curve_hw` of any computed curve, for the rest.

    Only (3) can destroy a real marker, and only one lying on the dashed curve.
    """
    import glyphfit as G
    b = ink > 0.5
    L = G.line_ink(ink, long_len, n_ang).astype(bool)
    S = G.line_ink(ink, short_len, n_ang).astype(bool)
    near = curve_mask(ink.shape, curve_hw)
    dash = curve_mask(ink.shape, dash_hw, which=(eq9,))
    remove = (L | tick_mask(ink, reach) | dash | (S & near)) & b
    return np.clip(b.astype(float) - remove, 0.0, 1.0), remove


def clean(ink, line_len=31, n_ang=12):
    """Remove every locally straight structure: the three printed curves, the
    two annotation leaders, the frame and its ticks.  A glyph 18-26 px across
    cannot contain a `line_len`-px straight run, so glyph ink survives.  The
    curves are removed as the ink they are, at the position the ink is in, but
    only after their computed positions have confirmed that this is what the
    straight ink IS - see calibration_check in the sidecar."""
    import glyphfit as G
    L = G.line_ink(ink, line_len, n_ang)
    return np.clip((ink > 0.5).astype(float) - L, 0.0, 1.0), L


def masks(ink, curve_hw=9):
    """care: 1 where ink is trustworthy evidence; allow: where a centre may sit."""
    h, w = ink.shape
    care = np.ones((h, w), float)
    care[:int(FRAME["top"]) + 4, :] = 0.0
    care[int(FRAME["bottom"]) - 4:, :] = 0.0
    care[:, :int(FRAME["left"]) + 4] = 0.0
    care[:, int(FRAME["right"]) - 4:] = 0.0
    # The frame, the ticks and the three curves are straight ink and are removed
    # by clean(); nothing is masked for them here, because a band mask wide
    # enough to cover a curve also swallows the markers sitting on it.
    if curve_hw:
        care[curve_mask((h, w), curve_hw)] = 0.0

    allow = np.zeros((h, w), bool)
    allow[int(FRAME["top"]) + 12:int(FRAME["bottom"]) - 12,
          int(FRAME["left"]) + 12:int(FRAME["right"]) - 12] = True
    for r0, r1, c0, c1 in TEXT_BOXES:
        allow[r0:r1, c0:c1] = False
    return care, allow


def trace_curves(ink, long_len=37, n_ang=18, window=20, max_width=8):
    """Offset of each printed curve's INK from its computed position, in px.

    Same one-thin-run-per-column rule as `setup3.trace_line`: a column counts
    only where the straight ink near the computed curve is a single run no wider
    than `max_width`, which is what keeps markers sitting on the curve out of the
    measurement.  Negative means the drawn ink is ABOVE the computed position,
    i.e. the calibration puts the curve too low.  Returns {name: (n, mean, sd)}.
    """
    import glyphfit as G
    straight = G.line_ink(ink, long_len, n_ang).astype(bool) & (ink > 0.5) \
        & ~tick_mask(ink)
    cols = np.arange(int(FRAME["left"]) + 25, int(FRAME["right"]) - 25)
    out = {}
    for name, f in (("eq12", eq12), ("eq11", eq11), ("eq9", eq9)):
        comp = row_of_sh(f(re_of_col(cols)))
        d = []
        for c, rc in zip(cols, comp):
            if not np.isfinite(rc) or rc < FRAME["top"] + 22 or rc > FRAME["bottom"] - 22:
                continue
            lo = int(rc) - window
            idx = np.where(straight[lo:lo + 2 * window + 1, c])[0]
            if len(idx) == 0:
                continue
            groups, cur = [], [idx[0]]
            for i in idx[1:]:
                if i - cur[-1] <= 1:
                    cur.append(i)
                else:
                    groups.append(cur); cur = [i]
            groups.append(cur)
            if len(groups) != 1 or len(groups[0]) > max_width:
                continue
            d.append(lo + float(np.mean(groups[0])) - rc)
        d = np.asarray(d)
        out[name] = (len(d), float(d.mean()) if len(d) else np.nan,
                     float(d.std()) if len(d) else np.nan)
    return out


SHAPES = ["circle", "square", "tri_up", "tri_down", "tri_right", "plus"]
SIZES = [18, 20, 22, 24, 26]
STROKE, R, SEARCH, GUARD = 0.26, 16, 5, 0.25


def prepare():
    """Binary ink, `care` (0 = no usable evidence), `allow` (centre may sit)."""
    ink = load()
    care0, allow = masks(ink, curve_hw=0)
    _, removed = clean2(ink, long_len=37, short_len=23, curve_hw=12, dash_hw=9, n_ang=18)
    return (ink > 0.5).astype(float), care0 * (1.0 - removed), allow, removed
