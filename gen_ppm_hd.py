#!/usr/bin/env python3
"""Build cat/ppm_hd.csv: the PPM -> HD identifications that PPM itself omits.

The VizieR merge in cat/ppm.txt leaves the HD field (bytes 109-114) empty for
248279 of its 468861 records, and empty for *all* 275 records of the bright
stars supplement (PPM 400001-400321), which also carry no SAO and no DM.  Since
read_ppm.cpp projects every old-catalogue match onto HD through that field, a
blank there silently drops the star from every results/cross/cross_*_hd.csv --
zeta Sculptoris (PPM 400001 = HD 224990) is matched by eight old catalogues and
appears in none of them.

This script recovers the missing identifications from the Bright Star Catalogue
(book/ybsc5.txt), which carries HD number, J2000 position and V magnitude for
its 9110 stars.  That bounds what can be recovered: the bright end, which is
also the only end the catalogue-and-atlas book uses.  Fainter PPM stars without
HD stay without HD -- no offline source here can give them one.

Matching rules, chosen to mimic what PPM does in the records it *does* fill in:

  * an edge joins a BSC5 star to a blank-HD PPM record within MAX_SEP arcsec,
    and is refused when both magnitudes are visual and differ by more than
    MAX_DMAG -- the same cutoffs, and for the same reason, as
    likelihood/cross_likelihood.py.  PPM magnitudes only count as visual where
    Flag5 is 'V' or the record is in the supplement; elsewhere the field is
    photographic and says nothing about V;
  * each edge is weighted (theta/sigma_pos)^2 + (dm/sigma_m)^2, the log-
    likelihood of cross_likelihood.py, and edges are settled cheapest first
    with each side used once.  Position alone is not enough: the supplement
    puts a close pair at the mean of the two components, so the BSC5 star at
    0.00" from PPM 400225 is alpha-2 Her (V 5.39) while the magnitude on the
    record, 3.6, plainly means alpha-1 (V 3.48).  Weighing both together picks
    the component PPM meant, and stops two PPM records over one double from
    ending up crossed;
  * when a pair collapses onto *one* PPM record, the loser keeps getting no
    note -- exactly what happens today for the pairs PPM itself fills in, where
    PPM 44721 is gamma-1 And, HD 12533, and gamma-2 = 12534 goes unnamed;
  * an HD already spoken for by another PPM record is refused, so the mapping
    stays injective and no cross_*_hd.csv can gain a duplicate;
  * a PPM record that already has an HD is left alone -- this file only fills
    blanks, it never overrides the catalogue.

Usage:  python gen_ppm_hd.py [--max-sep 20] [--check]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
PPM = ROOT / "cat" / "ppm.txt"
BSC5 = ROOT / "book" / "ybsc5.txt"
OUT = ROOT / "cat" / "ppm_hd.csv"

# The supplement sits at PPM 400001-400321; read_ppm.cpp already singles out
# that range when reading magnitudes (they are always visual there).
SUPPL = (400001, 400321)

# Byte offsets, 0-based half-open.  cat/ppm.txt is the VizieR layout shifted one
# column right in the identifier field, exactly as read_ppm.cpp reads it.
PPM_ID = (1, 7)
PPM_MAG = (19, 23)
PPM_RA = ((27, 29), (30, 32), (33, 39))
PPM_DE = (41, (42, 44), (45, 47), (48, 53))
PPM_HD = (108, 114)
PPM_FLAG5 = 130

# Edge weight, after likelihood/cross_likelihood.py.  SIGMA_POS is tighter than
# the 30" used there because both sides are modern catalogues: what it has to
# absorb is not observing error but the few arcsec by which the supplement's
# mean place for a double sits off each component.  DMAG_MISS is the flat
# penalty charged when the PPM magnitude is photographic or absent, so a
# position-only edge still loses to one that agrees in both.
SIGMA_POS = 5.0     # arcsec
SIGMA_MAG = 0.5     # mag
MAX_DMAG = 3.0      # hard cutoff, both magnitudes visual
DMAG_MISS = 2.0


def read_ppm() -> list[tuple[int, float, float, str, float]]:
    """(number, RA deg, Dec deg, HD or "", visual mag or nan) per PPM record.

    The magnitude is returned only where read_ppm.cpp would treat it as visual:
    Flag5 == 'V', or the record belongs to the bright stars supplement.
    """
    out = []
    for line in PPM.open():
        if len(line) < 95:
            continue
        num = line[slice(*PPM_ID)].strip()
        if not num.isdigit():
            continue
        try:
            h, m, s = (line[slice(*f)] for f in PPM_RA)
            ra = (float(h) + float(m) / 60 + float(s) / 3600) * 15
            sign, d, am, asec = PPM_DE
            de = float(line[slice(*d)]) + float(line[slice(*am)]) / 60 \
                + float(line[slice(*asec)]) / 3600
        except ValueError:
            continue
        if line[sign] == "-":
            de = -de
        n = int(num)
        mag = np.nan
        if line[PPM_FLAG5:PPM_FLAG5 + 1] == "V" or SUPPL[0] <= n <= SUPPL[1]:
            try:
                mag = float(line[slice(*PPM_MAG)])
            except ValueError:
                pass
        out.append((n, ra, de, line[slice(*PPM_HD)].strip(), mag))
    return out


def read_bsc5() -> list[tuple[str, str, float, float, float]]:
    """(HD, name, V, RA deg, Dec deg) for every BSC5 star with a position."""
    out = []
    for line in BSC5.open():
        hd = line[25:31].strip()
        if not hd or not line[75:77].strip():
            continue
        try:
            v = float(line[102:107])
        except ValueError:
            v = 99.0
        ra = (float(line[75:77]) + float(line[77:79]) / 60
              + float(line[79:83]) / 3600) * 15
        de = float(line[84:86]) + float(line[86:88]) / 60 + float(line[88:90]) / 3600
        if line[83] == "-":
            de = -de
        out.append((hd, line[4:14].strip(), v, ra, de))
    return out


def unit(ra: np.ndarray, de: np.ndarray) -> np.ndarray:
    r, d = np.radians(ra), np.radians(de)
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-sep", type=float, default=20.0,
                    help="match radius in arcsec (default 20)")
    ap.add_argument("--check", action="store_true",
                    help="report and exit without writing")
    args = ap.parse_args()

    ppm = read_ppm()
    bsc = read_bsc5()
    print(f"PPM records: {len(ppm)};  BSC5 stars with HD and position: {len(bsc)}",
          file=sys.stderr)

    # HD numbers PPM already uses: an HD there is not ours to hand out again.
    taken = {p[3].lstrip("0") for p in ppm if p[3]}
    blank = [p for p in ppm if not p[3]]
    print(f"PPM records with an empty HD field: {len(blank)}", file=sys.stderr)

    # candidate BSC5 stars: those PPM does not already place somewhere
    free = [s for s in bsc if s[0].lstrip("0") not in taken]
    print(f"BSC5 stars PPM does not already carry: {len(free)}", file=sys.stderr)

    PV = unit(np.array([p[1] for p in blank]), np.array([p[2] for p in blank]))
    BV = unit(np.array([s[3] for s in free]), np.array([s[4] for s in free]))
    chord = 2 * np.sin(np.radians(args.max_sep / 3600) / 2)
    pairs = cKDTree(PV).query_ball_tree(cKDTree(BV), chord)

    # every surviving edge, then greedy on the weight: the most convincing
    # pairing is settled first, so two PPM records over a double cannot swap.
    edges = []
    dropped_dmag = 0
    for j, hits in enumerate(pairs):
        for i in hits:
            sep = np.degrees(2 * np.arcsin(
                np.linalg.norm(PV[j] - BV[i]) / 2)) * 3600
            dm = abs(blank[j][4] - free[i][2])
            if dm != dm:                        # nan: PPM magnitude not visual
                dm = DMAG_MISS
            elif dm > MAX_DMAG:
                dropped_dmag += 1
                continue
            w = (sep / SIGMA_POS) ** 2 + (dm / SIGMA_MAG) ** 2
            edges.append((w, sep, j, i))
    edges.sort()

    cand: dict[int, tuple[str, str, float, float]] = {}
    used: set[str] = set()
    rejected_taken = collisions = 0
    for _w, sep, j, i in edges:
        num = blank[j][0]
        hd, name, v = free[i][0], free[i][1], free[i][2]
        if num in cand:                 # a nearer BSC5 star already took it
            collisions += 1
            continue
        if hd in used:                  # this HD already placed on another PPM
            rejected_taken += 1
            continue
        cand[num] = (hd, name, v, sep)
        used.add(hd)
    print(f"edges refused on magnitude (>{MAX_DMAG:g} mag): {dropped_dmag}",
          file=sys.stderr)

    supp = {n for n in cand if SUPPL[0] <= n <= SUPPL[1]}
    n_supp_total = sum(1 for p in ppm if SUPPL[0] <= p[0] <= SUPPL[1])
    print(f"identifications recovered: {len(cand)}"
          f"  ({len(supp)} of the {n_supp_total} supplement records,"
          f" {len(cand) - len(supp)} elsewhere)", file=sys.stderr)
    print(f"fainter components of a pair left over: {collisions};"
          f"  HD already used: {rejected_taken}", file=sys.stderr)

    if args.check:
        return 0

    with OUT.open("w") as f:
        f.write("# PPM -> HD identifications missing from cat/ppm.txt.\n")
        f.write("# Recovered from book/ybsc5.txt (BSC5) by position, "
                f"match radius {args.max_sep:g}\"; see gen_ppm_hd.py.\n")
        f.write("# Read by read_ppm.cpp, which uses it only where the PPM "
                "record's own HD field is blank.\n")
        f.write("ppm,hd,name,vmag,dist\n")
        for num in sorted(cand):
            hd, name, v, sep = cand[num]
            f.write(f"{num},{hd},{name or '-'},{v:.2f},{sep:.2f}\n")
    print(f"wrote {OUT} with {len(cand)} identifications", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
