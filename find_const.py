#!/usr/bin/env python3
"""
Script para encontrar la constelación dadas coordenadas ecuatoriales celestes.

Modos:
  --2000 RA_h RA_m RA_s Dec_d Dec_m Dec_s
      Coordenadas J2000.0 (ICRS).
  --1875 RA_h RA_m RA_s Dec_d Dec_m Dec_s
      Coordenadas B1875.0 (FK4).
  --csv archivo.csv
      Procesa un CSV cuyas últimas 4 columnas son x, y, z y un valor
      numérico (p. ej. "mag" o "dist"). Las columnas x, y, z se asumen
      en B1875.0 (FK4) y se reemplazan por la constelación; el resto se
      conserva. La salida se escribe junto al original con sufijo
      ".const" (p. ej. gilliss.const.csv).

Ejemplos:
  python3 find_const.py --2000 11 22 05.29 -24 46 39.8
  python3 find_const.py --1875 11 18 30.00 -24 30 00.0
  python3 find_const.py --csv results/doubles/gilliss.csv
  python3 find_const.py --csv results/cat1875/gc.csv
"""

import csv
import os
import sys

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


def make_skycoord(ra, dec, frame_kind):
    """Build a SkyCoord in the requested equinox/frame.

    `ra` and `dec` may be strings ("0h47m...", "-70d10m...") or astropy
    Quantities — SkyCoord handles both.
    """
    if frame_kind == 'j2000':
        return SkyCoord(ra, dec, frame='icrs')
    if frame_kind == 'b1875':
        return SkyCoord(ra, dec, frame='fk4', equinox='B1875.0')
    raise ValueError(f"frame desconocido: {frame_kind}")


def parse_six_coord_args(args):
    if len(args) != 6:
        print("Error: se requieren 6 argumentos (RA: h m s, Dec: d m s)", file=sys.stderr)
        usage_and_exit()
    ra_string = f"{args[0]}h{args[1]}m{args[2]}s"
    dec_string = f"{args[3]}d{args[4]}m{args[5]}s"
    return ra_string, dec_string


def find_constellation_hmsdms(ra_string, dec_string, frame_kind):
    coord = make_skycoord(ra_string, dec_string, frame_kind)
    return coord.get_constellation(short_name=True)


def rec_to_radec_deg(x, y, z):
    """Convert rectangular unit-sphere coords to RA/Dec in degrees.

    Mirrors trig.cpp:rec2sph (ra = atan2(y, x); dec = atan2(z, sqrt(x^2+y^2))).
    """
    ra_deg = np.degrees(np.arctan2(y, x)) % 360.0
    dec_deg = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
    return ra_deg, dec_deg


def process_csv(input_path):
    if not os.path.isfile(input_path):
        print(f"Error: archivo no encontrado: {input_path}", file=sys.stderr)
        sys.exit(1)

    base, ext = os.path.splitext(input_path)
    output_path = f"{base}.const{ext}"

    with open(input_path, newline='') as fin:
        rows = list(csv.reader(fin))

    if not rows:
        print("Error: CSV vacío", file=sys.stderr)
        sys.exit(1)

    header = rows[0]
    if len(header) < 4:
        print("Error: el CSV debe tener al menos 4 columnas (..., x, y, z, valor)", file=sys.stderr)
        sys.exit(1)

    new_header = header[:-4] + ['const'] + header[-1:]
    data_rows = rows[1:]

    with open(output_path, 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(new_header)

        if not data_rows:
            print(f"Escrito (sin filas): {output_path}")
            return

        xs = np.array([float(r[-4]) for r in data_rows])
        ys = np.array([float(r[-3]) for r in data_rows])
        zs = np.array([float(r[-2]) for r in data_rows])

        ra_deg, dec_deg = rec_to_radec_deg(xs, ys, zs)
        coords = make_skycoord(ra_deg * u.deg, dec_deg * u.deg, 'b1875')
        consts = coords.get_constellation(short_name=True)

        for row, const in zip(data_rows, consts):
            writer.writerow(row[:-4] + [const] + row[-1:])

    print(f"Escrito: {output_path}")


def usage_and_exit():
    print(__doc__, file=sys.stderr)
    sys.exit(1)


def main():
    args = sys.argv[1:]
    if not args:
        usage_and_exit()

    flags = {'--2000', '--1875', '--csv'}
    present = {a for a in args if a in flags}
    rest = [a for a in args if a not in flags]

    if present == {'--csv'}:
        if len(rest) != 1:
            usage_and_exit()
        process_csv(rest[0])
        return

    if present == {'--2000'}:
        ra_string, dec_string = parse_six_coord_args(rest)
        print(find_constellation_hmsdms(ra_string, dec_string, 'j2000'))
        return

    if present == {'--1875'}:
        ra_string, dec_string = parse_six_coord_args(rest)
        print(find_constellation_hmsdms(ra_string, dec_string, 'b1875'))
        return

    usage_and_exit()


if __name__ == "__main__":
    main()
