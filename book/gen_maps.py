#!/usr/bin/env python3
"""Generate the star-atlas plates that follow the catalogue.

Three main plates cover Dec +52..-52 in 10-hour slices of RA overlapping by
two hours, plus two half-disc plates for the south polar cap (Dec -45..-90).

Everything is exported as PDF so the plates stay vector art inside the book.
The stars are the very same BSC5 selection the catalogue table lists, so a dot
on a chart always has a row in the table.

starplot draws RA increasing to the left (the sky-chart convention) with Dec
increasing upward.  The book wants RA bottom-to-top and Dec left-to-right, so
the main plates are included in the LaTeX rotated 90 degrees *clockwise*
(\\includegraphics[angle=-90]), which maps left->up and up->right.

Usage:
    python gen_maps.py                 # all plates into book/maps/
    python gen_maps.py --only map1     # just one, while iterating
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_catalog import BSC5, parse_bsc5  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BOOK = ROOT / "book"
MAPS = BOOK / "maps"

CONSTELLATIONS = ROOT / "constellations.0.3.3.parquet"
HIPSTARS = ROOT / "stars.bigksy.0.1.3.mag11.parquet"

VMAX = 5.5          # same cut as the catalogue
DEC_BAND = 52.0     # main plates: +/- this
DEC_POLAR = -45.0   # polar cap starts here (7 deg of overlap with the band)

# marker size in points: linear in magnitude, so the brightest stars read as
# obviously larger without swamping the plate
MAG_FAINT = 6.2     # notional zero-size magnitude
SIZE_SCALE = 2.2


def plates() -> list[dict]:
    """The five plates, RA in degrees (may exceed 360 to express a wrap)."""
    # 11h per plate on an 8h step.  Miller stretches declination ~10%, so 11h
    # (not the 10h a linear projection would give) is what actually fills A4
    # once the plate is rotated; the extra hour is free overlap.
    out = []
    for i, ra0 in enumerate((0, 120, 240), start=1):
        out.append({
            "name": f"map{i}",
            "kind": "band",
            "ra_min": ra0, "ra_max": ra0 + 165,
            "dec_min": -DEC_BAND, "dec_max": DEC_BAND,
        })
    # south cap, split along the 0h-12h diameter
    out.append({"name": "polar_a", "kind": "polar",
                "ra_min": 180, "ra_max": 360,
                "dec_min": -90.0, "dec_max": DEC_POLAR})
    out.append({"name": "polar_b", "kind": "polar",
                "ra_min": 0, "ra_max": 180,
                "dec_min": -90.0, "dec_max": DEC_POLAR})
    return out


def ra_label(hours: float) -> str:
    """Gridline label for RA, given in hours: 8h, or 8h30 off the hour."""
    total = round(hours * 60)
    h, m = (total // 60) % 24, total % 60
    return f"{h}h" if m == 0 else f"{h}h{m:02d}"


def dec_label(deg: float) -> str:
    return f"{deg:+.0f}°"


# Klein-style frame metrics, in points measured outward from the map edge.
BAND_PT = 9.0        # gap between the two lines of the double map border
LABELPAD_PT = 7.0    # clearance from the outer border line to the labels

# Line weights, chosen for how they land once the plate is scaled onto the
# page (a plate is reduced roughly 3.1x, so these divide by about that).
LW_BORDER = 1.6      # ~0.50 pt printed
LW_TICK = 1.3        # ~0.40 pt printed
LW_GRID = 1.2        # ~0.38 pt printed
LW_BORDERS = 1.2     # ~0.38 pt printed
LW_FIGURES = 1.0     # ~0.32 pt printed

# Polar plates: gap from the rim to an RA label, and from the straight edge to
# a Dec label.  Both in points, so they stay isotropic on a non-square plate.
POLAR_LABEL_PT = 62.0
POLAR_DEC_PT = 26.0

# Colours for the on-screen PDF.  --mono gives the black-and-white set a
# printed edition wants.  The frame, ticks, labels and stars stay black in
# both: they are furniture, not data.
COLOURS = {
    "grid":     "#C9564B",   # light red
    "figures":  "#000000",   # constellation lines stay black
    "borders":  "#8A76C4",   # light violet
    "names":    "#2E8B57",   # constellation names
}
MONO = {"grid": "#666666", "figures": "#000000", "borders": "#888888",
        "names": "#000000"}

NAME_FS_FACTOR = 0.75   # constellation names, relative to the gridline labels


def draw_constellation_names(chart, p: dict, colour: str) -> None:
    """Name each constellation whose label anchor falls on this plate.

    A name is kept only if its whole box lies inside the map, so nothing runs
    over the border.  The test uses the real rendered text extent rather than
    an estimate from the character count, which would be wrong for the wide
    names ("Camelopardalis") that are exactly the ones at risk.
    """
    import cartopy.crs as ccrs
    ax = chart.ax
    fig = ax.figure

    df = pd.read_parquet(CONSTELLATIONS)[["name", "ra", "dec"]]
    fs = (chart.style.gridlines.label.font_size
          * getattr(chart, "scale", 1.0) * NAME_FS_FACTOR)

    placed = []
    for name, ra, dec in zip(df.name, df.ra, df.dec):
        if not in_fov(float(ra), float(dec), p):
            continue
        lon = -(float(ra) % 360.0)
        if lon < -180.0:
            lon += 360.0
        x, y = ax.projection.transform_point(lon, float(dec),
                                             ccrs.PlateCarree())
        placed.append(ax.text(x, y, str(name), color=colour, fontsize=fs,
                              ha="center", va="center", clip_on=False,
                              zorder=1004))

    if not placed:
        return 0

    fig.canvas.draw()
    try:
        renderer = fig.canvas.get_renderer()
    except AttributeError:
        renderer = fig._get_renderer()
    inv = ax.transAxes.inverted()

    if p["kind"] == "polar":
        px = 0.0 if p["opens_right"] else 1.0

        def inside(ax_x, ax_y):
            # the axes is exactly the half-disc's box, so the rim is the
            # ellipse with semi-axes 1.0 in x and 0.5 in y about the pole
            side = (ax_x >= px) if p["opens_right"] else (ax_x <= px)
            r = np.hypot((ax_x - px) / 1.0, (ax_y - 0.5) / 0.5)
            return side and r <= 1.0
    else:
        def inside(ax_x, ax_y):
            return 0.0 <= ax_x <= 1.0 and 0.0 <= ax_y <= 1.0

    kept = 0
    for t in placed:
        bb = t.get_window_extent(renderer=renderer)
        corners = [inv.transform((bx, by))
                   for bx, by in ((bb.x0, bb.y0), (bb.x1, bb.y0),
                                  (bb.x0, bb.y1), (bb.x1, bb.y1))]
        if all(inside(cx, cy) for cx, cy in corners):
            kept += 1
        else:
            t.remove()
    return kept


def draw_klein_frame(chart, p: dict, dec_locs) -> None:
    """Frame a band plate the way Klein's atlas does.

    From the map outwards: the map border as a double line, a tick band between
    those two lines carrying a mark every degree (the ticks stay inside the
    band and never cross the map), then the RA/Dec labels.  The plate's own
    bounding box then hugs the labels, and the page margins supply the rest.
    """
    import cartopy.crs as ccrs
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D

    ax = chart.ax
    fw, fh = ax.figure.get_size_inches()
    pos = ax.get_position()
    w_pt, h_pt = pos.width * fw * 72.0, pos.height * fh * 72.0
    fx, fy = (lambda v: v / w_pt), (lambda v: v / h_pt)

    ra_min, ra_max = float(p["ra_min"]), float(p["ra_max"])

    def xf(ra: float) -> float:
        """Axes x-fraction of an RA.  Miller is linear in longitude, and
        starplot draws RA increasing to the left, hence the 1 - ..."""
        return 1.0 - (ra - ra_min) / (ra_max - ra_min)

    y0, y1 = ax.get_ylim()

    def yf(dec: float) -> float:
        _, y = ax.projection.transform_point(0.0, dec, ccrs.PlateCarree())
        return (y - y0) / (y1 - y0)

    def add(artist):
        # cartopy clips anything added to a GeoAxes against the map boundary,
        # so frame parts go on the figure and borrow the axes transform
        ax.figure.add_artist(artist)
        artist.set_clip_on(False)
        artist.set_clip_path(None)

    # --- double map border -------------------------------------------------
    for lw, pad in ((LW_BORDER, 0.0), (LW_BORDER, BAND_PT)):
        add(Rectangle((-fx(pad), -fy(pad)),
                      1 + 2 * fx(pad), 1 + 2 * fy(pad),
                      transform=ax.transAxes, fill=False,
                      linewidth=lw, edgecolor="black", zorder=1002))

    # --- one tick per degree, living only inside the band ------------------
    for ra in range(int(np.ceil(ra_min)), int(np.floor(ra_max)) + 1):
        x = xf(float(ra))
        if not 0.0 <= x <= 1.0:
            continue
        for base, sgn in ((0.0, -1.0), (1.0, 1.0)):
            add(Line2D([x, x],
                       [base, base + sgn * fy(BAND_PT)],
                       transform=ax.transAxes, linewidth=LW_TICK,
                       color="black", zorder=1002))
    for dec in range(int(np.ceil(p["dec_min"])), int(np.floor(p["dec_max"])) + 1):
        y = yf(float(dec))
        if not 0.0 <= y <= 1.0:
            continue
        for base, sgn in ((0.0, -1.0), (1.0, 1.0)):
            add(Line2D([base, base + sgn * fx(BAND_PT)], [y, y],
                       transform=ax.transAxes, linewidth=LW_TICK,
                       color="black", zorder=1002))

    # --- labels, outside the band -----------------------------------------
    fs = chart.style.gridlines.label.font_size * getattr(chart, "scale", 1.0)
    kw = dict(transform=ax.transAxes, fontsize=fs, clip_on=False, zorder=1003)
    off = BAND_PT + LABELPAD_PT
    for ra in range(int(ra_min), int(ra_max) + 1, 15):
        x = xf(float(ra))
        if not 0.0 <= x <= 1.0:
            continue
        lab = ra_label(ra / 15.0)
        ax.text(x, 1 + fy(off), lab, ha="center", va="bottom", **kw)
        ax.text(x, -fy(off), lab, ha="center", va="top", **kw)
    for dec in dec_locs:
        y = yf(float(dec))
        if not 0.0 <= y <= 1.0:
            continue
        lab = dec_label(dec)
        ax.text(-fx(off), y, lab, ha="right", va="center", **kw)
        ax.text(1 + fx(off), y, lab, ha="left", va="center", **kw)



def draw_polar_labels(chart, p: dict, opens_right: bool, dec_locs) -> None:
    """Label a half-disc plate by hand.

    starplot's own gridline labels only come out on one of the two halves, and
    on that one it also stacks several RA labels on the pole where every
    meridian converges.  Placing them ourselves gives the same scheme as the
    band plates -- RA round the rim, Dec down the straight edge -- plus a
    single -90 at the pole.  starplot maps RA to longitude as lon = -RA
    (starplot.utils.ra_to_lon), which is what makes this projectable.
    """
    import cartopy.crs as ccrs
    ax = chart.ax
    inv = ax.transAxes.inverted()

    def axes_xy(ra_deg: float, dec: float):
        lon = -(ra_deg % 360.0)
        if lon < -180.0:
            lon += 360.0
        x, y = ax.projection.transform_point(lon, dec, ccrs.PlateCarree())
        return np.asarray(inv.transform(ax.transData.transform((x, y))))

    # A half-disc plate is about 1:2, so one unit of axes-x is a much shorter
    # physical distance than one unit of axes-y.  Offsetting labels by a fixed
    # axes fraction therefore pushed the sideways-pointing ones (the middle of
    # the page) barely half as far from the rim as the ones near the ends.
    # Work in points and convert per-axis so every label clears the arc equally.
    fw, fh = ax.figure.get_size_inches()
    pos = ax.get_position()
    w_pt, h_pt = pos.width * fw * 72.0, pos.height * fh * 72.0

    def offset_axes(direction, pts):
        """`pts` points along `direction`, returned in axes fractions."""
        dx, dy = direction[0] * w_pt, direction[1] * h_pt
        norm = np.hypot(dx, dy)
        if norm == 0:
            return np.zeros(2)
        return np.array([dx / norm * pts / w_pt, dy / norm * pts / h_pt])

    # match the size starplot uses for the band plates' gridline labels
    fs = chart.style.gridlines.label.font_size * getattr(chart, "scale", 1.0)
    kw = dict(transform=ax.transAxes, fontsize=fs, clip_on=False, zorder=1001)
    pole = axes_xy(0.0, -90.0)
    out = -1.0 if opens_right else 1.0     # away from the disc, along x

    # RA round the rim, pushed radially outward by a constant physical gap
    ra0, ra1 = int(p["ra_min"]), int(p["ra_max"])
    for ra in range(ra0, ra1 + 1, 15):
        q = axes_xy(float(ra), p["dec_max"])
        if ra in (ra0, ra1):
            # the two ends sit on the straight edge, not the arc
            up = 1.0 if q[1] > 0.5 else -1.0
            chart.ax.text(*(q + offset_axes((0.0, up), POLAR_LABEL_PT)),
                          ra_label(ra / 15.0), ha="center",
                          va="bottom" if up > 0 else "top", **kw)
            continue
        chart.ax.text(*(q + offset_axes(q - pole, POLAR_LABEL_PT)),
                      ra_label(ra / 15.0), ha="center", va="center", **kw)

    # Dec down both straight edges, and -90 at the pole itself
    for ra_edge in (ra0, ra1):
        for dec in dec_locs:
            q = axes_xy(float(ra_edge), dec)
            chart.ax.text(q[0] + out * POLAR_DEC_PT / w_pt, q[1], dec_label(dec),
                          ha="right" if opens_right else "left",
                          va="center", **kw)
    chart.ax.text(pole[0] + out * POLAR_DEC_PT / w_pt, pole[1], dec_label(-90),
                  ha="right" if opens_right else "left", va="center", **kw)


def ra_in(ra: float, lo: float, hi: float) -> bool:
    """RA membership that tolerates a window running past 360 degrees."""
    ra = ra % 360.0
    if hi <= 360.0:
        return lo <= ra <= hi
    return ra >= lo or ra <= hi - 360.0


def in_fov(ra: float, dec: float, p: dict) -> bool:
    return (p["dec_min"] <= dec <= p["dec_max"]
            and ra_in(ra, p["ra_min"], p["ra_max"]))


def load_hip_positions() -> dict[int, tuple[float, float]]:
    # column projection trips this file's dataset schema, so read then subset
    df = pd.read_parquet(HIPSTARS)[["hip", "ra", "dec"]].dropna(subset=["hip"])
    return {int(h): (float(r), float(d))
            for h, r, d in zip(df.hip, df.ra, df.dec)}


def load_constellation_lines() -> list[tuple[int, int]]:
    df = pd.read_parquet(CONSTELLATIONS, columns=["star_hip_lines"])
    pairs = []
    for lines in df.star_hip_lines:
        if lines is None:
            continue
        for pair in lines:
            pairs.append((int(pair[0]), int(pair[1])))
    return pairs


def marker_size(mag: float) -> float:
    return max(0.6, (MAG_FAINT - mag) * SIZE_SCALE)


def build_plate(p: dict, stars, hip_pos, lines, out_dir: Path,
                colours: dict = COLOURS) -> Path:
    from starplot import (LineStyle, MapPlot, MarkerStyle, ObjectStyle,
                          PathStyle, Stereographic, Miller)

    clip = None
    if p["kind"] == "band":
        centre = ((p["ra_min"] + p["ra_max"]) / 2.0) % 360.0
        projection = Miller(center_ra=centre)
    else:
        projection = Stereographic(center_ra=0.0, center_dec=-90.0)
        # Half-disc: out along one radius, round the Dec limit, back along the
        # other radius to the pole.  Without this the plate is drawn as a bare
        # rectangle and the grid spills past the rim into the corners.
        from shapely.geometry import Polygon
        arc = np.linspace(p["ra_min"], p["ra_max"], 361)
        pts = ([(p["ra_min"], -90.0)]
               + [(float(a), p["dec_max"]) for a in arc]
               + [(p["ra_max"], -90.0)])
        clip = Polygon(pts)

    chart = MapPlot(
        projection=projection,
        ra_min=p["ra_min"], ra_max=p["ra_max"],
        dec_min=p["dec_min"], dec_max=p["dec_max"],
        clip_path=clip,
        resolution=3000,
    )

    if p["kind"] == "polar":
        # clip_path filters objects but leaves the frame rectangular; cartopy's
        # set_boundary is what actually cuts the plate to a half-disc.  The
        # extent is exactly the half-disc bounding box (1:2), but which edge
        # the pole sits on depends on the RA half -- 180-360 puts it left,
        # 0-180 puts it right -- so locate it rather than assume.
        import cartopy.crs as ccrs
        import matplotlib.path as mpath
        px, _ = chart.ax.projection.transform_point(0.0, -90.0,
                                                    ccrs.PlateCarree())
        x0, x1 = chart.ax.get_xlim()
        opens_right = ((px - x0) / (x1 - x0)) < 0.5
        p["opens_right"] = opens_right

        # keep the traversal direction the same for both halves, otherwise the
        # reversed winding leaves the rim spine unstroked on the mirrored plate
        t = np.linspace(-np.pi / 2, np.pi / 2, 361)
        if opens_right:
            verts = np.column_stack([np.cos(t), 0.5 + 0.5 * np.sin(t)])
        else:
            t = t[::-1]
            verts = np.column_stack([1.0 - np.cos(t), 0.5 + 0.5 * np.sin(t)])
        chart.ax.set_boundary(mpath.Path(verts, closed=True),
                              transform=chart.ax.transAxes)
        # cartopy strokes the boundary spine only for one winding direction, so
        # draw the rim ourselves and both halves come out looking the same
        chart.ax.plot(np.append(verts[:, 0], verts[0, 0]),
                      np.append(verts[:, 1], verts[0, 1]),
                      transform=chart.ax.transAxes,
                      color="black", linewidth=0.8, clip_on=False, zorder=1000)

    # Labels in the manner of Klein's atlas: RA every hour along the RA axis,
    # Dec every 10 degrees along the Dec axis, on both edges.  They rotate with
    # the plate, so on the printed page they read along their own axis.
    grid_kw = dict(
        style=PathStyle(line=LineStyle(style="dotted", width=LW_GRID, alpha=0.9,
                                       color=colours["grid"])),
        labels=True,
        ra_formatter_fn=ra_label,
        dec_formatter_fn=dec_label,
    )
    polar_dec_locs = list(range(-80, -44, 10))
    band_dec_locs = list(range(-50, 51, 10))
    if p["kind"] == "band":
        grid_kw["ra_locations"] = list(range(0, 361, 15))          # 1h steps
        grid_kw["dec_locations"] = band_dec_locs
        grid_kw["labels"] = False       # placed by draw_klein_frame instead
    else:
        grid_kw["ra_locations"] = list(range(0, 360, 15))
        grid_kw["dec_locations"] = polar_dec_locs
        grid_kw["labels"] = False        # placed by draw_polar_labels instead
    chart.gridlines(**grid_kw)
    if p["kind"] == "polar":
        draw_polar_labels(chart, p, p["opens_right"], polar_dec_locs)
    else:
        draw_klein_frame(chart, p, band_dec_locs)
    chart.constellation_borders(
        style=LineStyle(style="dashed", width=LW_BORDERS, alpha=0.9,
                        color=colours["borders"]))

    # Constellation figures: a segment is drawn only when *both* endpoints fall
    # inside this plate, so no line is left dangling off the edge.
    seg_style = PathStyle(line=LineStyle(style="solid", width=LW_FIGURES,
                                         alpha=1.0, color=colours["figures"]))
    named = draw_constellation_names(chart, p, colours["names"])

    drawn = 0
    for h1, h2 in lines:
        a, b = hip_pos.get(h1), hip_pos.get(h2)
        if not a or not b:
            continue
        if in_fov(*a, p) and in_fov(*b, p):
            chart.line(coordinates=[a, b], style=seg_style)
            drawn += 1

    plotted = 0
    for s in stars:
        if not in_fov(s.ra_deg, s.de_deg, p):
            continue
        chart.marker(
            ra=s.ra_deg, dec=s.de_deg,
            style=ObjectStyle(marker=MarkerStyle(
                size=marker_size(s.vmag), fill="full", color="black")),
        )
        plotted += 1

    out = out_dir / f"{p['name']}.pdf"
    chart.export(str(out))
    print(f"  {p['name']}: {plotted} stars, {drawn} constellation segments, "
          f"{named} names -> {out.name}", file=sys.stderr)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="build a single plate by name")
    ap.add_argument("--out", type=Path, default=MAPS)
    ap.add_argument("--mono", action="store_true",
                    help="black-and-white plates for a printed edition")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    stars = parse_bsc5(BSC5, VMAX, 90.0)
    print(f"stars available: {len(stars)}", file=sys.stderr)
    hip_pos = load_hip_positions()
    lines = load_constellation_lines()
    print(f"HIP positions: {len(hip_pos)}, constellation segments: {len(lines)}",
          file=sys.stderr)

    for p in plates():
        if args.only and p["name"] != args.only:
            continue
        build_plate(p, stars, hip_pos, lines, args.out,
                    MONO if args.mono else COLOURS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
