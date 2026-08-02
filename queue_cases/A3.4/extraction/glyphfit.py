"""Shape-fitted marker extraction.

The centring method, in one sentence: every marker is located by fitting a
parametric glyph template (the shapes the figure's series actually use) to the
ink at sub-pixel resolution, and the reported coordinate is the CENTRE OF THE
FITTED SHAPE, not a local maximum of ink density.

Score of a template of size s placed with its symmetry centre at (y, x):

    score = <ink over the glyph body> - <ink over a thin outer margin>

both averages weighted by a `care` mask that is zero on pixels belonging to the
frame, the ticks and the curves whose equations the paper prints (erased at
their COMPUTED position, never traced).  The margin term is what makes the score
peak at the centre: sliding the template off-centre both loses body ink and
gains margin ink, so the optimum is a genuine registration of the shape rather
than the centroid of a blob.

Fitting is a discrete search over (shape, size, integer offset) followed by a
Nelder-Mead polish of (dy, dx, log size) with an explicit initial simplex.

Detection is matching pursuit: take the best-scoring (shape, position, size),
accept it, remove the ink it explains, repeat.  Merged glyphs are therefore
explained one at a time instead of collapsing into a single centroid.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from matplotlib.path import Path
from numpy.lib.stride_tricks import sliding_window_view
from scipy import ndimage as ndi
from scipy.optimize import minimize
from scipy.signal import fftconvolve

SS = 3          # supersampling factor for rasterisation


# --------------------------------------------------------------------------
# glyph geometry.  Vertices live in a unit box: the SYMMETRY CENTRE of the
# glyph is the origin and `size` is the bounding-box width in pixels.  A hand
# drafter centres a symbol on its data point by eye, so the symmetry centre --
# not the area centroid -- is the right model of where the datum is.
# --------------------------------------------------------------------------
def _poly(n, start=90.0):
    a = np.deg2rad(start + np.arange(n) * 360.0 / n)
    return np.stack([np.cos(a), np.sin(a)], 1)


def _star(n=5, inner=0.45, rot=0.0):
    a = np.deg2rad(90 + rot + np.arange(2 * n) * 180.0 / n)
    r = np.where(np.arange(2 * n) % 2 == 0, 1.0, inner)
    return np.stack([r * np.cos(a), r * np.sin(a)], 1)


def _norm(p):
    p = np.asarray(p, float)
    return p / (2.0 * np.abs(p).max())


POLY = {
    "circle":    _norm(_poly(48)),
    "square":    _norm(_poly(4, start=45.0)),
    "diamond":   _norm(_poly(4)),
    "tri_up":    _norm(_poly(3, start=90.0)),
    "tri_down":  _norm(_poly(3, start=-90.0)),
    "tri_left":  _norm(_poly(3, start=180.0)),
    "tri_right": _norm(_poly(3, start=0.0)),
    "star":      _norm(_star()),
    "hexag":     _norm(_poly(6)),
    "wide_dia":  _norm(_poly(4) * np.array([1.5, 1.0])),
}
ALL_SHAPES = tuple(POLY) + ("plus",)


def _inside(poly, yy, xx, scale=1.0):
    return Path(poly * scale).contains_points(
        np.stack([xx.ravel(), yy.ravel()], 1)).reshape(xx.shape)


@lru_cache(maxsize=200000)
def render(shape, size, filled=True, stroke=0.26, R=16, dy=0.0, dx=0.0):
    """Rasterise one glyph.  Returns (body, margin) on a (2R+1)^2 grid."""
    g = (np.arange((2 * R + 1) * SS) + 0.5) / SS - 0.5 - R
    yy, xx = np.meshgrid(g - dy, g - dx, indexing="ij")
    # GAP: dead zone just outside the glyph, absorbing ink spread and the
    # drafter's line weight.  FAR: outer edge of the margin annulus.  The
    # annulus must be several pixels WIDE or it cannot constrain the fitted
    # size, and every fit collapses onto the smallest template that will sit
    # inside the blob (measured, 2026-08-02: with a 1-px annulus all 130
    # Figure 3 fits came back as under-sized triangles scoring 1.00).
    GAP, FAR = 1.16, 1.62
    if shape == "plus":
        arm, hw = 0.5 * size, 0.5 * stroke * size
        body = ((np.abs(xx) <= arm) & (np.abs(yy) <= hw)) | \
               ((np.abs(yy) <= arm) & (np.abs(xx) <= hw))
        inner = ((np.abs(xx) <= arm * GAP) & (np.abs(yy) <= hw * GAP + 0.10 * size)) | \
                ((np.abs(yy) <= arm * GAP) & (np.abs(xx) <= hw * GAP + 0.10 * size))
        margin = ((np.abs(xx) <= arm * FAR) & (np.abs(yy) <= arm * FAR)) & ~inner
    else:
        poly = POLY[shape]
        out = _inside(poly, yy, xx, size)
        body = out if filled else (out & ~_inside(poly, yy, xx, size * (1 - 2 * stroke)))
        margin = _inside(poly, yy, xx, size * FAR) & ~_inside(poly, yy, xx, size * GAP)
    n = 2 * R + 1
    body = body.reshape(n, SS, n, SS).mean((1, 3))
    margin = margin.reshape(n, SS, n, SS).mean((1, 3))
    body.flags.writeable = False
    margin.flags.writeable = False
    return body, margin


def _key(v, q=0.05):
    return round(round(v / q) * q, 4)


def tmpl(shape, size, filled, stroke, R, dy=0.0, dx=0.0):
    return render(shape, _key(size, 0.25), filled, _key(stroke, 0.01), R,
                  _key(dy, 0.05), _key(dx, 0.05))


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def window(img, y, x, R, fill=0.0):
    h, w = img.shape
    out = np.full((2 * R + 1, 2 * R + 1), fill, float)
    y0, x0 = int(round(y)) - R, int(round(x)) - R
    sy0, sx0 = max(0, y0), max(0, x0)
    sy1, sx1 = min(h, y0 + 2 * R + 1), min(w, x0 + 2 * R + 1)
    if sy1 > sy0 and sx1 > sx0:
        out[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = img[sy0:sy1, sx0:sx1]
    return out


def score_at(ink, care, y, x, shape, size, filled, stroke, R, margin_w=1.0):
    body, marg = tmpl(shape, size, filled, stroke, R, y - round(y), x - round(x))
    I = window(ink, y, x, R)
    W = window(care, y, x, R)
    wb, wm = (W * body).sum(), (W * marg).sum()
    if wb < 1e-6 or wm < 1e-6:
        return -1.0
    return float((W * body * I).sum() / wb - margin_w * (W * marg * I).sum() / wm)


def body_cover(img, care, y, x, shape, size, filled, stroke, R):
    """Plain weighted coverage of the glyph body by `img` (no margin term)."""
    body, _ = tmpl(shape, size, filled, stroke, R, y - round(y), x - round(x))
    I = window(img, y, x, R)
    W = window(care, y, x, R)
    wb = (W * body).sum()
    return float((W * body * I).sum() / wb) if wb > 1e-6 else 0.0


def evidence(ink, care, y, x, shape, size, filled, stroke, R, floor=0.45):
    """Cared-for ink under the glyph body, normalised by the FULL body area.

    Unlike `score_at` this does not shrink its denominator when part of the
    glyph is masked away, so it cannot be inflated by a sliver of surviving ink
    next to an erased curve -- the failure mode that filled the first detection
    pass with triangles sitting on the Ranz-Marshall dashes.
    """
    body, _ = tmpl(shape, size, filled, stroke, R, y - round(y), x - round(x))
    I = window(ink, y, x, R)
    W = window(care, y, x, R)
    den = max((W * body).sum(), floor * body.sum())
    return float((W * body * I).sum() / den) if den > 1e-6 else 0.0


def fit_glyph(ink, care, y0, x0, shapes, sizes, filled=False, stroke=0.26, R=16,
              margin_w=1.0, search=8, polish=True, guard=0.35):
    """Discrete (shape, size, integer offset) search, then sub-pixel polish."""
    S = int(search)
    n = 2 * R + 1
    pi = window(ink, y0, x0, R + S)
    pc = window(care, y0, x0, R + S)
    Wi = sliding_window_view(pi, (n, n))
    Wc = sliding_window_view(pc, (n, n))
    A = Wi * Wc
    best = (-9.0, None, None, 0, 0)
    for sh in shapes:
        for size in sizes:
            body, marg = tmpl(sh, size, filled, stroke, R)
            nb = np.tensordot(A, body, axes=([2, 3], [0, 1]))
            db = np.tensordot(Wc, body, axes=([2, 3], [0, 1]))
            nm = np.tensordot(A, marg, axes=([2, 3], [0, 1]))
            dm = np.tensordot(Wc, marg, axes=([2, 3], [0, 1]))
            with np.errstate(invalid="ignore", divide="ignore"):
                sc = nb / np.maximum(db, 1e-6) - margin_w * nm / np.maximum(dm, 1e-6)
            sc[db < guard * body.sum()] = -9.0
            sc[dm < guard * marg.sum()] = -9.0
            i = int(np.argmax(sc))
            iy, ix = np.unravel_index(i, sc.shape)
            if sc[iy, ix] > best[0]:
                best = (float(sc[iy, ix]), sh, size, int(iy) - S, int(ix) - S)
    s0, sh, size, dy, dx = best
    if sh is None:          # nothing measurable here (all support masked away)
        return dict(y=float(y0), x=float(x0), shape=shapes[0],
                    size=float(sizes[len(sizes) // 2]), score=-9.0)
    y, x = float(y0 + dy), float(x0 + dx)
    if not polish:
        return dict(y=y, x=x, shape=sh, size=size, score=s0)

    def neg(p):
        ddy, ddx, ls = p
        if abs(ddy) > 2.5 or abs(ddx) > 2.5 or abs(ls) > 0.22:
            return 9.0
        return -score_at(ink, care, y + ddy, x + ddx, sh, size * np.exp(ls),
                         filled, stroke, R, margin_w)

    simplex = np.array([[0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [0.0, 0.8, 0.0],
                        [0.0, 0.0, 0.10]])
    r = minimize(neg, simplex[0], method="Nelder-Mead",
                 options=dict(initial_simplex=simplex, xatol=0.03, fatol=2e-4,
                              maxiter=200))
    ddy, ddx, ls = r.x
    return dict(y=y + float(ddy), x=x + float(ddx), shape=sh,
                size=float(size * np.exp(ls)), score=float(-r.fun))


# --------------------------------------------------------------------------
# matching pursuit
# --------------------------------------------------------------------------
def line_ink(ink, length=31, n_ang=12, thr=0.5):
    """Ink explainable by a locally STRAIGHT structure: the union over
    orientations of a binary opening with a long thin line element.  Curves,
    leader lines, frame and ticks survive it; a glyph 20-26 px across cannot
    contain a 31-px straight run, so it does not."""
    b = ink > thr
    out = np.zeros_like(b)
    h = length // 2
    for a in np.linspace(0, np.pi, n_ang, endpoint=False):
        se = np.zeros((length, length), bool)
        t = np.linspace(-h, h, 4 * length)
        rr = np.clip(np.round(h + t * np.sin(a)).astype(int), 0, length - 1)
        cc = np.clip(np.round(h + t * np.cos(a)).astype(int), 0, length - 1)
        se[rr, cc] = True
        out |= ndi.binary_opening(b, se)
    return out.astype(float)


def coverage_map(res, care, shape, size, filled, stroke, R, margin_w):
    body, marg = tmpl(shape, size, filled, stroke, R)
    rc = res * care
    nb = fftconvolve(rc, body[::-1, ::-1], mode="same")
    db = fftconvolve(care, body[::-1, ::-1], mode="same")
    nm = fftconvolve(rc, marg[::-1, ::-1], mode="same")
    dm = fftconvolve(care, marg[::-1, ::-1], mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        s = nb / np.maximum(db, 1e-6) - margin_w * nm / np.maximum(dm, 1e-6)
    s[db < 0.35 * body.sum()] = -9.0
    s[dm < 0.35 * marg.sum()] = -9.0
    return s


def subtract(res, y, x, shape, size, filled, stroke, R):
    body, _ = tmpl(shape, size, filled, stroke, R, y - round(y), x - round(x))
    y0, x0 = int(round(y)) - R, int(round(x)) - R
    ys = slice(max(0, y0), min(res.shape[0], y0 + 2 * R + 1))
    xs = slice(max(0, x0), min(res.shape[1], x0 + 2 * R + 1))
    bb = body[ys.start - y0:ys.stop - y0, xs.start - x0:xs.stop - x0]
    res[ys, xs] *= (1.0 - bb)


def pursue_batch(ink, care, allow, shapes, sizes, filled=False, stroke=0.26,
                 R=16, margin_w=1.0, thresh=0.45, keep_thresh=0.60, min_sep=11.0,
                 rounds=5, seeds=(), search=5, verbose=True, nl_ink=None,
                 nl_min=0.30, guard=0.35, ev_min=0.0, ev_floor=0.45):
    """Matching pursuit in rounds: all peaks of one round are fitted, kept and
    subtracted together, then the maps are recomputed.  Equivalent in spirit to
    one-at-a-time pursuit but ~50x cheaper."""
    res = ink.copy()
    found = [dict(g) for g in seeds]
    for g in found:
        subtract(res, g["y"], g["x"], g["shape"], g["size"], filled, stroke, R)
    size0 = sizes[len(sizes) // 2]
    for rnd in range(rounds):
        smap = None
        for sh in shapes:
            m = coverage_map(res, care, sh, size0, filled, stroke, R, margin_w)
            smap = m if smap is None else np.maximum(smap, m)
        smap[~allow] = -9.0
        loc = ndi.maximum_filter(smap, size=int(min_sep))
        peaks = np.argwhere((smap >= loc) & (smap > thresh))
        order = np.argsort(-smap[peaks[:, 0], peaks[:, 1]])
        added = 0
        for iy, ix in peaks[order]:
            if any((iy - f["y"]) ** 2 + (ix - f["x"]) ** 2 < (0.8 * min_sep) ** 2
                   for f in found):
                continue
            g = fit_glyph(ink, care, float(iy), float(ix), shapes, sizes, filled,
                          stroke, R, margin_w, search=search, guard=guard)
            if any((g["y"] - f["y"]) ** 2 + (g["x"] - f["x"]) ** 2 < min_sep ** 2
                   for f in found):
                continue
            g["res_score"] = score_at(res, care, g["y"], g["x"], g["shape"],
                                      g["size"], filled, stroke, R, margin_w)
            g["evidence"] = evidence(ink, care, g["y"], g["x"], g["shape"],
                                     g["size"], filled, stroke, R, ev_floor)
            if g["evidence"] < ev_min:
                continue
            if nl_ink is not None:
                g["nl_cover"] = body_cover(nl_ink, care, g["y"], g["x"],
                                           g["shape"], g["size"], filled, stroke, R)
            if g["score"] < keep_thresh or g.get("nl_cover", 1.0) < nl_min:
                continue
            g["round"] = rnd
            found.append(g)
            subtract(res, g["y"], g["x"], g["shape"], g["size"], filled, stroke, R)
            added += 1
        if verbose:
            print(f"    round {rnd}: +{added} -> {len(found)} markers")
        if added == 0:
            break
    return found, res


def pursue(ink, care, allow, shapes, sizes, filled=False, stroke=0.26, R=16,
           margin_w=1.0, thresh=0.40, keep_thresh=0.42, min_sep=10.0,
           max_n=400, verbose=False, seeds=()):
    """Matching pursuit.  `seeds` are centres already accepted (removed first)."""
    res = ink.copy()
    found = []
    for g in seeds:
        found.append(dict(g))
        subtract(res, g["y"], g["x"], g["shape"], g["size"], filled, stroke, R)
    blocked = np.zeros(ink.shape, bool)
    size0 = sizes[len(sizes) // 2]
    while len(found) < max_n:
        best_s, best_sh, best_yx = -9.0, None, None
        for sh in shapes:
            m = coverage_map(res, care, sh, size0, filled, stroke, R, margin_w)
            m[~allow] = -9.0
            m[blocked] = -9.0
            i = int(np.argmax(m))
            if float(m.flat[i]) > best_s:
                best_s = float(m.flat[i]); best_sh = sh
                best_yx = np.unravel_index(i, m.shape)
        if best_s < thresh:
            break
        y0, x0 = float(best_yx[0]), float(best_yx[1])
        g = fit_glyph(ink, care, y0, x0, shapes, sizes, filled, stroke, R,
                      margin_w, search=4)
        far = all((g["y"] - f["y"]) ** 2 + (g["x"] - f["x"]) ** 2 > min_sep ** 2
                  for f in found)
        rs = score_at(res, care, g["y"], g["x"], g["shape"], g["size"],
                      filled, stroke, R, margin_w)
        g["res_score"] = rs
        if far and g["score"] >= keep_thresh:
            found.append(g)
            if verbose and len(found) % 20 == 0:
                print(f"    {len(found)} markers, last {g['shape']} "
                      f"score {g['score']:.3f}")
        else:
            blocked[max(0, int(y0) - 3):int(y0) + 4,
                    max(0, int(x0) - 3):int(x0) + 4] = True
        subtract(res, g["y"], g["x"], g["shape"], g["size"], filled, stroke, R)
    return found, res
