#!/usr/bin/env python3
"""
Genera una carta celeste de una constelación a partir de un CSV en B1875.0.

Uso:
    python3 gen_atlas.py <input.csv> <CONST_3LETRAS> <output.png> [modo]

Donde:
    <input.csv>  CSV con columnas (name, x, y, z, mag) en frame B1875.0 (FK4),
                 al estilo de results/cat1875/*.csv.  No usar las variantes
                 .const.csv (esquema distinto).
    <CONST>      Código IAU de 3 letras, sin distinguir mayúsculas
                 (Oct, Ori, CMa, ...).
    <output.png> Ruta de salida (PNG).
    [modo]       Modo de etiqueta para estrellas con mag <= 6.0 (por defecto: full).
                   full   – nombre completo ("GC 1234", "bet CMa")
                   number – solo el número final ("1234"); si no hay número, nombre completo
                   mag    – magnitud × 10 como entero ("54" para mag 5.4)

El script convierte internamente las coordenadas B1875.0 -> ICRS y dibuja con
starplot:
  - los bordes IAU de la constelación con líneas a trazos cortos,
  - una grilla RA/Dec cada 10 grados con líneas punteadas,
  - todas las estrellas del CSV (ninguna se descarta), con radio escalado por
    magnitud y rotuladas según el modo elegido para EPSILON < mag <= 6.0.

Primera ejecución: starplot descarga su catálogo + efemérides (~50 MB) en su
directorio de datos.

Ejemplos:
    python3 gen_atlas.py results/cat1875/gc.csv CMa cma_gc.png number
    python3 gen_atlas.py results/cat1875/ppm.csv CMa cma_ppm.png mag
"""

import csv
import os
import re
import sys

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

from starplot import (
    CollisionHandler,
    Constellation,
    LabelStyle,
    LineStyle,
    MapPlot,
    MarkerStyle,
    ObjectStyle,
    PathStyle,
    Star,
    Stereographic,
    _,
)


LABEL_MODE = "full"   # "full" | "number" | "mag"

EPSILON = 1e-6        # mag <= EPSILON treated as the "missing" placeholder
MAG_LIMIT = 10.0      # magnitudes >= this share the smallest visible radius
LABEL_MAG_MAX = 6.0   # only stars brighter than this get a label
PLACEHOLDER_MAG = 9.0 # value substituted for the mag=0 placeholder


def ra_label(r):
    """Format RA gridline label (r in hours).  Show minutes only when non-zero."""
    total_minutes = round(r * 60)
    h = total_minutes // 60
    m = total_minutes % 60
    return f"{h}h" if m == 0 else f"{h}h{m:02d}"


def rec_to_radec_deg(x, y, z):
    """(x, y, z) on the unit sphere -> (RA, Dec) in degrees.  Mirrors trig.cpp."""
    ra_deg = np.degrees(np.arctan2(y, x)) % 360.0
    dec_deg = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
    return ra_deg, dec_deg


def label_for(name, mag, mode):
    if mode == "mag":
        effective_mag = PLACEHOLDER_MAG if mag <= EPSILON else mag
        if effective_mag > LABEL_MAG_MAX:
            return None
        return str(int(round(effective_mag * 10)))
    if mag <= EPSILON or mag > LABEL_MAG_MAX:
        return None
    if mode == "full":
        return name
    if mode == "number":
        m = re.search(r"(\d+)$", name)
        return m.group(1) if m else name
    raise ValueError(f"LABEL_MODE desconocido: {mode}")


def load_csv(path):
    """Return list of dicts {name, x, y, z, mag} from the input CSV."""
    with open(path, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)
    if not rows:
        sys.exit(f"Error: CSV vacío: {path}")
    header = rows[0]
    if len(header) < 5:
        sys.exit(f"Error: el CSV debe tener al menos 5 columnas (name, x, y, z, mag): {path}")
    out = []
    for r in rows[1:]:
        out.append({
            "name": r[0],
            "x": float(r[-4]),
            "y": float(r[-3]),
            "z": float(r[-2]),
            "mag": float(r[-1]),
        })
    return out


def b1875_xyz_to_icrs_radec(stars):
    """Batch-convert B1875 (x,y,z) records to ICRS (ra_deg, dec_deg) arrays."""
    xs = np.array([s["x"] for s in stars])
    ys = np.array([s["y"] for s in stars])
    zs = np.array([s["z"] for s in stars])
    ra1875, dec1875 = rec_to_radec_deg(xs, ys, zs)
    coords_b1875 = SkyCoord(
        ra=ra1875 * u.deg,
        dec=dec1875 * u.deg,
        frame="fk4",
        equinox="B1875.0",
    )
    icrs = coords_b1875.icrs
    return icrs.ra.deg, icrs.dec.deg


def compute_fov(boundary):
    """From a constellation boundary polygon, decide the chart FoV + projection.

    Returns (ra_min, ra_max, dec_min, dec_max, projection).  The boundary may
    be a normal polygon or one that wraps RA/encloses a pole.
    """
    ra_lo, dec_lo, ra_hi, dec_hi = boundary.bounds
    pad = 5.0

    south_polar = dec_lo <= -89.5
    north_polar = dec_hi >= 89.5
    ra_full = (ra_hi - ra_lo) > 350.0  # boundary spans most of the RA circle

    if south_polar or (ra_full and dec_hi < 0):
        # Polar south: Stereographic centered on the south pole with a
        # rectangular frame (shows RA/Dec labels on all four edges).
        return (
            0.0, 360.0,
            -90.0, min(90.0, dec_hi + pad),
            Stereographic(center_ra=0.0, center_dec=-90.0),
        )
    if north_polar or (ra_full and dec_lo > 0):
        return (
            0.0, 360.0,
            max(-90.0, dec_lo - pad), 90.0,
            Stereographic(center_ra=0.0, center_dec=90.0),
        )

    # Normal case: rectangular bbox + 5 deg pad on each side.
    ra_min = ra_lo - pad
    ra_max = ra_hi + pad
    dec_min = max(-90.0, dec_lo - pad)
    dec_max = min(90.0, dec_hi + pad)
    center_ra = (ra_lo + ra_hi) / 2.0
    center_dec = (dec_lo + dec_hi) / 2.0
    return ra_min, ra_max, dec_min, dec_max, Stereographic(
        center_ra=center_ra, center_dec=center_dec,
    )


def draw_constellation_lines(chart, iau_id, line_style):
    """Draw constellation stick figure by fetching HIP pairs and using chart.line().

    Using chart.line() instead of chart.constellations() bypasses starplot's
    60°-RA wrap-detection heuristic, which incorrectly routes lines the long
    way around for polar constellations (e.g. Oct) whose stars are physically
    adjacent but far apart in RA.
    """
    constellation = Constellation.get(iau_id=iau_id)
    hip_pos = {}
    for hip in constellation.star_hip_ids:
        try:
            s = Star.get(hip=hip)
            hip_pos[hip] = (s.ra, s.dec)
        except Exception:
            pass

    style = PathStyle(line=line_style)
    for hip1, hip2 in constellation.star_hip_lines:
        p1 = hip_pos.get(hip1)
        p2 = hip_pos.get(hip2)
        if p1 and p2:
            chart.line(coordinates=[p1, p2], style=style)


def usage_and_exit():
    print(__doc__, file=sys.stderr)
    sys.exit(1)


def main():
    args = sys.argv[1:]
    if len(args) not in (3, 4):
        usage_and_exit()
    input_csv, const_arg, output_png = args[:3]
    label_mode = args[3] if len(args) == 4 else LABEL_MODE

    if label_mode not in ("full", "number", "mag"):
        sys.exit(f"Error: modo desconocido '{label_mode}'. Usar full, number o mag.")
    if not os.path.isfile(input_csv):
        sys.exit(f"Error: archivo no encontrado: {input_csv}")
    if not output_png.lower().endswith(".png"):
        sys.exit("Error: la salida debe ser .png")

    iau_id = const_arg.strip().lower()
    constellation = Constellation.get(iau_id=iau_id)
    if constellation is None:
        sys.exit(f"Error: constelación IAU desconocida: {const_arg}")

    stars = load_csv(input_csv)
    ra_icrs, dec_icrs = b1875_xyz_to_icrs_radec(stars)

    ra_min, ra_max, dec_min, dec_max, projection = compute_fov(constellation.boundary)

    chart = MapPlot(
        projection=projection,
        ra_min=ra_min, ra_max=ra_max,
        dec_min=dec_min, dec_max=dec_max,
        # Never suppress labels: allow overlaps with markers, other labels,
        # and clipped positions; plot at the fallback position on failure.
        point_label_handler=CollisionHandler(
            allow_marker_collisions=True,
            allow_label_collisions=True,
            allow_clipped=True,
            plot_on_fail=True,
        ),
    )

    chart.gridlines(
        style=PathStyle(line=LineStyle(style="dotted", width=1)),
        ra_locations=list(range(0, 360, 10)),
        dec_locations=list(range(-80, 90, 10)),
        ra_formatter_fn=ra_label,
    )

    # Force all four frame borders visible — cartopy hides left/right for polar
    # stereographic projections by default.
    for side in ("left", "right", "top", "bottom"):
        chart.ax.spines[side].set_visible(True)
        chart.ax.spines[side].set_linewidth(0.5)

    chart.constellation_borders(style=LineStyle(style="dashed", width=1))

    # Draw constellation stick figure manually (one chart.line() per HIP pair)
    # to avoid starplot's 60°-RA wrap heuristic, which breaks polar constellations.
    draw_constellation_lines(
        chart, iau_id, LineStyle(style="solid", width=1, alpha=0.4),
    )

    for s, ra, dec in zip(stars, ra_icrs, dec_icrs):
        mag = s["mag"]
        effective_mag = PLACEHOLDER_MAG if mag <= EPSILON else mag
        radius_pt = max(2.0, 2.0 * (MAG_LIMIT - effective_mag))
        label = label_for(s["name"], mag, label_mode)
        chart.marker(
            ra=ra, dec=dec,
            style=ObjectStyle(
                marker=MarkerStyle(size=radius_pt, fill="full", color="black"),
                label=LabelStyle(font_size=14),
            ),
            label=label,
        )

    chart.export(output_png)
    print(f"Escrito: {output_png}")


if __name__ == "__main__":
    main()
