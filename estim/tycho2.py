#!/usr/bin/env python3
"""Tycho-2 photometric catalogue: build a compact cache and query it.

Build once from the raw Tycho-2 distribution:

    python3 tycho2.py --build

then in other scripts:

    from tycho2 import Tycho2
    cat = Tycho2()
    idx, sep_arcsec, iso_arcsec = cat.match(ra_deg, dec_deg)

Field layout (pipe-delimited, from the CDS ReadMe for catalog I/259):
main catalog (cat/tyc2.txt), 0-indexed after ``line.split('|')``:
  [0]  "TYC1 TYC2 TYC3" (whitespace-separated header block)
  [2]  RAmdeg   mean R.A., ICRS, epoch J2000
  [3]  DEmdeg   mean Dec., ICRS, epoch J2000
  [17] BTmag
  [19] VTmag
  [20] e_VTmag
supplement (cat/tyc2_suppl.txt), no proper motion / mean-epoch fields:
  [0]  "TYC1 TYC2 TYC3"
  [2]  RAdeg
  [3]  DEdeg
  [11] BTmag
  [13] VTmag
  [14] e_VTmag

Johnson conversion (gen_tycho2.cpp:517-527), applied only when BT is present:
  V   = VT - 0.090 * (BT - VT)
  B-V = 0.850 * (BT - VT)
Otherwise V = VT and B-V is unknown (NaN). Valid only for BT-VT < 1.8 (B-V < 1.53).
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_TXT = HERE.parent / "cat" / "tyc2.txt"
DEFAULT_SUPPL = HERE.parent / "cat" / "tyc2_suppl.txt"
DEFAULT_NPY = HERE.parent / "cat" / "tyc2_phot.npy"

DTYPE = np.dtype([
    ("ra", "f8"),
    ("dec", "f8"),
    ("V", "f4"),
    ("BV", "f4"),
    ("e_VT", "f4"),
    ("tyc1", "i4"),
    ("tyc2", "i4"),
    ("tyc3", "i2"),
])


def _f(s, default=0.0):
    s = s.strip()
    return float(s) if s else default


def _parse_main_line(line):
    f = line.split("|")
    if len(f) < 21:
        return None
    vt = _f(f[19])
    if vt <= 0.0:
        return None
    bt = _f(f[17])
    if bt > 0.0:
        v = vt - 0.090 * (bt - vt)
        bv = 0.850 * (bt - vt)
    else:
        v = vt
        bv = np.nan
    e_vt = _f(f[20], np.nan)
    ra = _f(f[2], np.nan)
    dec = _f(f[3], np.nan)
    if not (0.0 <= ra <= 360.0) or not (-90.0 <= dec <= 90.0):
        return None
    head = f[0].split()
    tyc1, tyc2, tyc3 = int(head[0]), int(head[1]), int(head[2])
    return (ra, dec, v, bv, e_vt, tyc1, tyc2, tyc3)


def _parse_suppl_line(line):
    f = line.split("|")
    if len(f) < 15:
        return None
    vt = _f(f[13])
    if vt <= 0.0:
        return None
    bt = _f(f[11])
    if bt > 0.0:
        v = vt - 0.090 * (bt - vt)
        bv = 0.850 * (bt - vt)
    else:
        v = vt
        bv = np.nan
    e_vt = _f(f[14], np.nan)
    ra = _f(f[2], np.nan)
    dec = _f(f[3], np.nan)
    if not (0.0 <= ra <= 360.0) or not (-90.0 <= dec <= 90.0):
        return None
    head = f[0].split()
    tyc1, tyc2, tyc3 = int(head[0]), int(head[1]), int(head[2])
    return (ra, dec, v, bv, e_vt, tyc1, tyc2, tyc3)


def build(txt_path=DEFAULT_TXT, suppl_path=DEFAULT_SUPPL, out_path=DEFAULT_NPY):
    t0 = time.time()
    rows = []
    n_lines = 0
    n_kept = 0
    with open(txt_path, "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            n_lines += 1
            rec = _parse_main_line(line)
            if rec is not None:
                rows.append(rec)
                n_kept += 1
            if n_lines % 500000 == 0:
                print(f"  main: {n_lines} lines, {n_kept} kept ({time.time()-t0:.0f}s)",
                      file=sys.stderr)
    print(f"main catalog: {n_lines} lines, {n_kept} with usable VT", file=sys.stderr)

    n_suppl = 0
    n_suppl_kept = 0
    if Path(suppl_path).exists():
        with open(suppl_path, "r", encoding="ascii", errors="replace") as fh:
            for line in fh:
                n_suppl += 1
                rec = _parse_suppl_line(line)
                if rec is not None:
                    rows.append(rec)
                    n_suppl_kept += 1
        print(f"supplement: {n_suppl} lines, {n_suppl_kept} with usable VT", file=sys.stderr)
    else:
        print(f"supplement not found at {suppl_path}, skipping", file=sys.stderr)

    arr = np.array(rows, dtype=DTYPE)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(out_path, arr)
    print(f"wrote {len(arr)} stars to {out_path} ({time.time()-t0:.0f}s total)", file=sys.stderr)
    return arr


def _radec_to_xyz(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    cd = np.cos(dec)
    return np.stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)], axis=-1)


def _chord_to_arcsec(chord):
    # chord length between unit vectors -> angular separation
    return np.degrees(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))) * 3600.0


class Tycho2:
    def __init__(self, npy_path=DEFAULT_NPY):
        self.data = np.load(npy_path)

    def __len__(self):
        return len(self.data)

    def designation(self, i):
        """Standard TYC designation, third component 1-3. Part of the common
        ``starcat`` interface -- the Gaia backends carry no identifiers and synthesise a
        coordinate-based one instead."""
        rec = self.data[i]
        return f"TYC {rec['tyc1']}-{rec['tyc2']}-{rec['tyc3']}"

    def box(self, ra_min, ra_max, dec_min, dec_max):
        """Return the sub-array of stars within a RA/Dec box (deg). Handles RA wraparound."""
        d = self.data
        dec_mask = (d["dec"] >= dec_min) & (d["dec"] <= dec_max)
        if ra_min <= ra_max:
            ra_mask = (d["ra"] >= ra_min) & (d["ra"] <= ra_max)
        else:
            ra_mask = (d["ra"] >= ra_min) | (d["ra"] <= ra_max)
        return d[dec_mask & ra_mask]

    def cone(self, ra_deg, dec_deg, radius_deg):
        """Sub-array of stars within radius_deg of a single point, via a padded box + exact cut."""
        pad = radius_deg / max(np.cos(np.radians(dec_deg)), 1e-6)
        sub = self.box(ra_deg - pad, ra_deg + pad, dec_deg - radius_deg, dec_deg + radius_deg)
        if len(sub) == 0:
            return sub
        xyz = _radec_to_xyz(sub["ra"], sub["dec"])
        center = _radec_to_xyz(np.array([ra_deg]), np.array([dec_deg]))[0]
        chord = np.linalg.norm(xyz - center, axis=1)
        sep_arcsec = _chord_to_arcsec(chord)
        return sub[sep_arcsec <= radius_deg * 3600.0]

    def match(self, ra_deg, dec_deg, radius_arcsec=8.0, isolation_arcsec=20.0):
        """Match arrays of field-star (ra, dec) [deg] against the catalogue.

        Returns three arrays, one per input star:
          idx           index into self.data of the nearest catalogue star, or -1 if
                        none within radius_arcsec
          sep_arcsec    separation to that star (NaN if idx == -1)
          nn2_arcsec    separation from the *matched* catalogue star to its own
                        nearest neighbour in the catalogue (NaN if idx == -1); compare
                        against isolation_arcsec to flag blends
        """
        from scipy.spatial import cKDTree

        ra_deg = np.atleast_1d(np.asarray(ra_deg, dtype="f8"))
        dec_deg = np.atleast_1d(np.asarray(dec_deg, dtype="f8"))
        n = len(ra_deg)
        idx_out = np.full(n, -1, dtype="i8")
        sep_out = np.full(n, np.nan, dtype="f8")
        nn2_out = np.full(n, np.nan, dtype="f8")
        if n == 0:
            return idx_out, sep_out, nn2_out

        ra_min, ra_max = ra_deg.min(), ra_deg.max()
        dec_min, dec_max = dec_deg.min(), dec_deg.max()
        pad = max(radius_arcsec, isolation_arcsec) / 3600.0
        pad_ra = pad / max(np.cos(np.radians(max(abs(dec_min), abs(dec_max), 1e-6))), 1e-6)
        sub = self.box(ra_min - pad_ra, ra_max + pad_ra, dec_min - pad, dec_max + pad)
        if len(sub) == 0:
            return idx_out, sep_out, nn2_out

        # map sub-array positions back to indices in self.data
        d = self.data
        dec_mask = (d["dec"] >= dec_min - pad) & (d["dec"] <= dec_max + pad)
        if ra_min - pad_ra <= ra_max + pad_ra:
            ra_mask = (d["ra"] >= ra_min - pad_ra) & (d["ra"] <= ra_max + pad_ra)
        else:
            ra_mask = (d["ra"] >= ra_min - pad_ra) | (d["ra"] <= ra_max + pad_ra)
        sub_global_idx = np.nonzero(dec_mask & ra_mask)[0]
        sub = d[sub_global_idx]

        cat_xyz = _radec_to_xyz(sub["ra"], sub["dec"])
        tree = cKDTree(cat_xyz)

        # nearest neighbour of each catalogue star to itself (k=2, skip self)
        if len(sub) > 1:
            self_chord, self_ii = tree.query(cat_xyz, k=2)
            cat_nn_arcsec = _chord_to_arcsec(self_chord[:, 1])
        else:
            cat_nn_arcsec = np.full(len(sub), np.inf)

        field_xyz = _radec_to_xyz(ra_deg, dec_deg)
        chord, ii = tree.query(field_xyz, k=1)
        sep_arcsec = _chord_to_arcsec(chord)

        matched = sep_arcsec <= radius_arcsec
        idx_out[matched] = sub_global_idx[ii[matched]]
        sep_out[matched] = sep_arcsec[matched]
        nn2_out[matched] = cat_nn_arcsec[ii[matched]]
        return idx_out, sep_out, nn2_out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="store_true", help="build the .npy cache from the raw txt")
    ap.add_argument("--txt", default=str(DEFAULT_TXT))
    ap.add_argument("--suppl", default=str(DEFAULT_SUPPL))
    ap.add_argument("-o", "--out", default=str(DEFAULT_NPY))
    ap.add_argument("--stats", action="store_true", help="print summary stats of the cache")
    args = ap.parse_args()

    if args.build:
        build(args.txt, args.suppl, args.out)
    if args.stats or not args.build:
        cat = Tycho2(args.out)
        d = cat.data
        n_bv = np.sum(~np.isnan(d["BV"]))
        print(f"{len(d)} stars, {n_bv} with B-V, V range [{d['V'].min():.2f}, {d['V'].max():.2f}]")


if __name__ == "__main__":
    main()
