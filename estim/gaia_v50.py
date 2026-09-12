#!/usr/bin/env python3
"""Reader for ASTAP's local star databases (V50, D80, ... in .1476 / .290 form).

V50 is a Gaia DR3 extract carrying **Johnson V plus B-V**, which makes it a drop-in
alternative to Tycho-2 as the photometric reference. This module exposes the same
surface as ``tycho2.Tycho2`` -- ``box()``, ``cone()``, ``match()`` and a ``.data``
structured array -- so the rest of the pipeline can switch backends without caring.

    from gaia_v50 import GaiaV50
    cat = GaiaV50()
    idx, sep_arcsec, nn2_arcsec = cat.match(ra_deg, dec_deg)

Unlike tycho2.py there is **no --build step**. A full conversion would be ~190 million
stars (several GB); instead each 5.14-degree area file is decoded on demand in a few
milliseconds and cached. ``match()`` loads whatever areas the query needs and appends
them to ``self.data``, so indices handed out earlier stay valid.

Format decoded from ``unit_star_database.pas`` in the ASTAP source (github.com/
CanardConfit/ASTAP) -- layout comments at lines 39-181, record types at 140-154, decode
at 2954-2996, declination boundaries at 2080-2119, area/cell arithmetic at 2323-2688:

* **Header, 110 bytes, no magic.** Bytes 0-107 ASCII description; byte 108 = version;
  byte 109 = record size (5 or 6). ``n_records = (filesize - 110) / record_size``.
* **Sky partition: 36 declination rings** (not HEALPix), cell counts per ring in
  ``CELLS_PER_RING`` summing to 1476. Ring step 90/17.5 = 5.142857 deg with half-height
  polar caps. Cell index ``= int(ra_deg * n_cells / 360)``. File ``<db>_<RR><CC>.1476``.
* **Differential compression.** A marker record (``ra7=ra8=ra9=0xFF``) sets the running
  group state: ``dec9 = dec7 - 128`` (signed high byte of Dec) and
  ``mag = (dec8 - 16)/10``. Ordinary records carry only the low bytes and inherit both.
* **Decode.** ``ra_deg = (b0 | b1<<8 | b2<<16) * 360/16777215``;
  ``dec_deg = ((dec9<<16) | (b4<<8) | b3) * 90/8388607``; magnitude from the current
  marker; if record size is 6, ``colour = int8(b5)``.
* **Colour byte semantics depend on header byte 108**: version 2 means (B-V)*50,
  anything else means (BP-RP)*10 (converted here via ASTAP's own polynomials).
  ``-128`` is the only "unknown" sentinel -- byte 0 is a genuine, merely rare colour bin.

Two properties that matter for photometry and are easy to miss:

* **V is quantised to 0.1 mag** (uniform quantisation error 0.029 mag), while **B-V is
  quantised to 0.02 mag**. Colour is stored five times finer than the magnitude.
* **B-V is capped at +-2.54** by the signed byte. Stars redder than that come back as
  unknown rather than clipped, which is the safe failure but does mean carbon stars get
  no reference colour.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
DEFAULT_DB_PATH = Path(os.environ.get("ASTAP_DB_PATH", "/usr/local/opt/astap"))
DEFAULT_DB_NAME = "v50"

# Cell counts per declination ring, south to north (sums to 1476).
CELLS_PER_RING = [1, 3, 9, 15, 21, 27, 33, 38, 43, 48, 52, 56, 60, 63, 65, 67, 68, 69,
                  69, 68, 67, 65, 63, 60, 56, 52, 48, 43, 38, 33, 27, 21, 15, 9, 3, 1]
N_RINGS = len(CELLS_PER_RING)

# Ring edges in degrees. Copied verbatim from ASTAP's dec_boundaries1476 table
# (unit_star_database.pas:2080-2119) rather than computed from the 5.142857 deg step:
# the table's 8-decimal literals do not round-trip exactly, and a computed equator edge
# lands on -8.6e-10 instead of 0, which flips the ring for anything sitting exactly on a
# boundary. Matching the writer's own constants removes that class of bug entirely.
DEC_EDGES = np.array([
    -90.0, -87.42857143, -82.28571429, -77.14285714, -72.0, -66.85714286,
    -61.71428571, -56.57142857, -51.42857143, -46.28571429, -41.14285714, -36.0,
    -30.85714286, -25.71428571, -20.57142857, -15.42857143, -10.28571429, -5.142857143,
    0.0,
    5.142857143, 10.28571429, 15.42857143, 20.57142857, 25.71428571, 30.85714286, 36.0,
    41.14285714, 46.28571429, 51.42857143, 56.57142857, 61.71428571, 66.85714286, 72.0,
    77.14285714, 82.28571429, 87.42857143, 90.0,
])
assert len(DEC_EDGES) == N_RINGS + 1

RA_SCALE = 360.0 / ((256 ** 3) - 1)        # 0xFFFFFF is reserved as the marker
DEC_SCALE = 90.0 / ((128 * 256 * 256) - 1)
COLOUR_UNKNOWN = -128
MAG_QUANTISATION_SIGMA = 0.1 / np.sqrt(12.0)   # 0.029 mag

# Common dtype, deliberately close to tycho2.DTYPE so the two are interchangeable.
DTYPE = np.dtype([
    ("ra", "f8"),
    ("dec", "f8"),
    ("V", "f4"),
    ("BV", "f4"),
    ("e_V", "f4"),
])


def bv_from_bp_rp(bp_rp):
    """B-V from Gaia BP-RP, following ASTAP's own route (unit_online_gaia.pas:84-88):
    predict Tycho BT/VT from BP-RP, then B-V = 0.850*(BT-VT). The Gaia G term cancels in
    the difference, so this depends on colour alone. Valid for -0.3 < BP-RP < 3.0.

    Only needed for version-1 databases, where the stored byte is (BP-RP)*10 rather than
    (B-V)*50. The installed V50 is version 2 and does not use this path.
    """
    x = np.asarray(bp_rp, dtype=float)
    bt_minus_vt = (-0.006482 + 0.7865 * x - 0.3631 * x ** 2
                   + 0.93192 * x ** 3 - 0.4843 * x ** 4 + 0.06814 * x ** 5)
    bv = 0.850 * bt_minus_vt
    return np.where((x > -0.3) & (x < 3.0), bv, np.nan)


def ring_of_dec(dec_deg):
    """1-based ring number containing a declination.

    ASTAP's ``area_and_boundaries1476`` walks the boundary table downward testing
    ``dec1 > dec_boundaries[k]`` (unit_star_database.pas:2327+), so a ring spans
    ``(lower, upper]`` -- lower exclusive, upper inclusive. ``side='left'`` reproduces
    that. The distinction only bites for a declination sitting exactly on a boundary
    (e.g. Dec 0.0, a ring edge), but it decides which file the *writer* put the star in,
    so it has to match.
    """
    r = int(np.searchsorted(DEC_EDGES, dec_deg, side="left"))
    return min(max(r, 1), N_RINGS)


def area_name(ra_deg, dec_deg):
    """'RRCC' area label for a position."""
    ring = ring_of_dec(dec_deg)
    n = CELLS_PER_RING[ring - 1]
    cell = int(np.floor((ra_deg % 360.0) * n / 360.0)) + 1
    cell = min(max(cell, 1), n)
    return f"{ring:02d}{cell:02d}"


def areas_for_box(ra_min, ra_max, dec_min, dec_max):
    """Every 'RRCC' area overlapping a RA/Dec box. Handles RA wraparound and the
    single-cell polar caps."""
    out = []
    r_lo = ring_of_dec(max(dec_min, -90.0))
    r_hi = ring_of_dec(min(dec_max, 90.0))
    for ring in range(r_lo, r_hi + 1):
        n = CELLS_PER_RING[ring - 1]
        if n == 1:
            out.append(f"{ring:02d}01")
            continue
        lo = ra_min % 360.0
        hi = ra_max % 360.0
        if (ra_max - ra_min) >= 360.0:
            cells = range(n)
        elif lo <= hi:
            c0 = int(np.floor(lo * n / 360.0))
            c1 = int(np.floor(hi * n / 360.0))
            cells = range(c0, min(c1, n - 1) + 1)
        else:  # wraps through RA = 0
            c0 = int(np.floor(lo * n / 360.0))
            c1 = int(np.floor(hi * n / 360.0))
            cells = list(range(c0, n)) + list(range(0, min(c1, n - 1) + 1))
        for c in cells:
            out.append(f"{ring:02d}{c + 1:02d}")
    return out


def read_area_file(path, mag_limit=None):
    """Decode one .1476/.290 area file into a DTYPE array.

    Vectorised: the marker/differential stream is resolved with
    ``np.maximum.accumulate`` to forward-fill each group's state, so there is no Python
    loop over records.
    """
    with open(path, "rb") as fh:
        header = fh.read(110)
        if len(header) < 110:
            raise ValueError(f"{path}: truncated header")
        version = header[108]
        rec_size = header[109]
        if rec_size == 32:      # space means "default", i.e. the 11-byte HNSKY layout
            rec_size = 11
        if rec_size not in (5, 6):
            raise ValueError(f"{path}: unsupported record size {rec_size}")
        buf = fh.read()

    n = len(buf) // rec_size
    a = np.frombuffer(buf[: n * rec_size], dtype=np.uint8).reshape(n, rec_size)

    ra_raw = (a[:, 0].astype(np.uint32)
              | (a[:, 1].astype(np.uint32) << 8)
              | (a[:, 2].astype(np.uint32) << 16))
    is_marker = ra_raw == 0xFFFFFF

    # forward-fill the state set by the most recent marker
    marker_pos = np.maximum.accumulate(np.where(is_marker, np.arange(n), -1))
    seen = marker_pos >= 0
    dec9_at = np.where(is_marker, a[:, 3].astype(np.int16) - 128, 0)
    mag10_at = np.where(is_marker, a[:, 4].astype(np.int16) - 16, 0)
    safe = np.where(seen, marker_pos, 0)
    cur_dec9 = np.where(seen, dec9_at[safe], 0)
    cur_mag10 = np.where(seen, mag10_at[safe], 0)

    star = (~is_marker) & seen
    if mag_limit is not None:
        star &= cur_mag10 <= int(round(mag_limit * 10))
    if not np.any(star):
        return np.empty(0, dtype=DTYPE), version

    out = np.empty(int(star.sum()), dtype=DTYPE)
    out["ra"] = ra_raw[star] * RA_SCALE
    dec_raw = ((cur_dec9[star].astype(np.int32) << 16)
               | (a[star, 4].astype(np.int32) << 8)
               | a[star, 3].astype(np.int32))
    out["dec"] = dec_raw * DEC_SCALE
    out["V"] = cur_mag10[star] / 10.0
    out["e_V"] = MAG_QUANTISATION_SIGMA

    if rec_size >= 6:
        cb = a[star, 5].astype(np.int8)
        known = cb != COLOUR_UNKNOWN
        if version == 2:
            bv = np.where(known, cb / 50.0, np.nan)
        else:
            bv = np.where(known, bv_from_bp_rp(cb / 10.0), np.nan)
        out["BV"] = bv
    else:
        out["BV"] = np.nan       # D80 and friends carry no colour

    return out, version


def _radec_to_xyz(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    cd = np.cos(dec)
    return np.stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)], axis=-1)


def _chord_to_arcsec(chord):
    return np.degrees(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))) * 3600.0


class GaiaV50:
    """ASTAP local star database, loaded lazily one area at a time.

    ``mag_limit`` defaults to 13.0: the Dwarf sensors do not usefully detect fainter
    than about V = 12-13 at SNR 5, and the cut keeps each cached area around 0.3 MB
    instead of 4 MB. Pass ``None`` for everything (useful for blend checks against very
    faint neighbours, at a real memory cost across many fields).
    """

    def __init__(self, db_path=DEFAULT_DB_PATH, db_name=DEFAULT_DB_NAME,
                 mag_limit=13.0, ext=None):
        self.db_path = Path(db_path)
        self.db_name = db_name
        self.mag_limit = mag_limit
        self.ext = ext or self._detect_ext()
        self.version = None
        self.data = np.empty(0, dtype=DTYPE)
        self._loaded = set()
        if not self.db_path.is_dir():
            raise FileNotFoundError(f"ASTAP database directory not found: {self.db_path}")

    def _detect_ext(self):
        for ext in ("1476", "290"):
            if (self.db_path / f"{self.db_name}_0101.{ext}").exists():
                return ext
        raise FileNotFoundError(
            f"no {self.db_name}_0101.1476 or .290 in {self.db_path}; is the database installed?")

    def __len__(self):
        return len(self.data)

    def _area_path(self, label):
        return self.db_path / f"{self.db_name}_{label}.{self.ext}"

    def ensure_areas(self, labels):
        """Decode and append any of these areas not already loaded. Appending keeps
        previously returned indices into ``self.data`` valid."""
        new = [lab for lab in dict.fromkeys(labels) if lab not in self._loaded]
        chunks = []
        for lab in new:
            p = self._area_path(lab)
            if not p.exists():
                self._loaded.add(lab)      # polar/edge areas may legitimately not exist
                continue
            arr, version = read_area_file(p, self.mag_limit)
            if self.version is None:
                self.version = version
            chunks.append(arr)
            self._loaded.add(lab)
        if chunks:
            self.data = np.concatenate([self.data] + chunks) if len(self.data) \
                else np.concatenate(chunks)
        return len(new)

    def box(self, ra_min, ra_max, dec_min, dec_max):
        """Stars within a RA/Dec box (deg), loading areas as needed. Handles RA wrap."""
        self.ensure_areas(areas_for_box(ra_min, ra_max, dec_min, dec_max))
        d = self.data
        if len(d) == 0:
            return d
        dec_mask = (d["dec"] >= dec_min) & (d["dec"] <= dec_max)
        if ra_min <= ra_max:
            ra_mask = (d["ra"] >= ra_min) & (d["ra"] <= ra_max)
        else:
            ra_mask = (d["ra"] >= ra_min) | (d["ra"] <= ra_max)
        return d[dec_mask & ra_mask]

    def cone(self, ra_deg, dec_deg, radius_deg):
        """Stars within radius_deg of a point, via a padded box then an exact cut."""
        pad = radius_deg / max(np.cos(np.radians(dec_deg)), 1e-6)
        sub = self.box(ra_deg - pad, ra_deg + pad, dec_deg - radius_deg, dec_deg + radius_deg)
        if len(sub) == 0:
            return sub
        center = _radec_to_xyz(np.array([ra_deg]), np.array([dec_deg]))[0]
        sep = _chord_to_arcsec(np.linalg.norm(_radec_to_xyz(sub["ra"], sub["dec"]) - center, axis=1))
        return sub[sep <= radius_deg * 3600.0]

    def match(self, ra_deg, dec_deg, radius_arcsec=8.0, isolation_arcsec=20.0):
        """Match field stars against the catalogue. Same contract as ``tycho2.Tycho2``:

          idx         index into ``self.data``, or -1 if nothing within radius_arcsec
          sep_arcsec  separation to that star (NaN when unmatched)
          nn2_arcsec  separation from the *matched catalogue star* to its own nearest
                      catalogue neighbour (NaN when unmatched) -- compare against
                      isolation_arcsec to flag blends
        """
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
        worst_dec = max(abs(dec_min), abs(dec_max), 1e-6)
        pad_ra = pad / max(np.cos(np.radians(worst_dec)), 1e-6)
        self.ensure_areas(areas_for_box(ra_min - pad_ra, ra_max + pad_ra,
                                         dec_min - pad, dec_max + pad))

        d = self.data
        if len(d) == 0:
            return idx_out, sep_out, nn2_out
        dec_mask = (d["dec"] >= dec_min - pad) & (d["dec"] <= dec_max + pad)
        lo, hi = ra_min - pad_ra, ra_max + pad_ra
        if lo <= hi:
            ra_mask = (d["ra"] >= lo) & (d["ra"] <= hi)
        else:
            ra_mask = (d["ra"] >= lo) | (d["ra"] <= hi)
        sub_idx = np.nonzero(dec_mask & ra_mask)[0]
        if len(sub_idx) == 0:
            return idx_out, sep_out, nn2_out
        sub = d[sub_idx]

        cat_xyz = _radec_to_xyz(sub["ra"], sub["dec"])
        tree = cKDTree(cat_xyz)
        if len(sub) > 1:
            self_chord, _ = tree.query(cat_xyz, k=2)
            cat_nn = _chord_to_arcsec(self_chord[:, 1])
        else:
            cat_nn = np.full(len(sub), np.inf)

        chord, ii = tree.query(_radec_to_xyz(ra_deg, dec_deg), k=1)
        sep = _chord_to_arcsec(chord)
        hit = sep <= radius_arcsec
        idx_out[hit] = sub_idx[ii[hit]]
        sep_out[hit] = sep[hit]
        nn2_out[hit] = cat_nn[ii[hit]]
        return idx_out, sep_out, nn2_out

    def designation(self, i):
        """V50 carries no identifiers, so synthesise a stable coordinate-based one.
        5 decimal places in degrees is ~0.036 arcsec, finer than the 0.077/0.039 arcsec
        storage quantisation, so it is unique and reproducible."""
        rec = self.data[i]
        return f"V50 J{rec['ra']:09.5f}{rec['dec']:+09.5f}"


def _selftest(db_path, db_name, area, mag_lo, mag_hi):
    """Decode one area and cross-check against Tycho-2. This is the gate: if positions
    and magnitudes do not line up here, nothing downstream is meaningful."""
    from tycho2 import Tycho2

    path = Path(db_path) / f"{db_name}_{area}.1476"
    if not path.exists():
        print(f"FAIL: {path} not found", file=sys.stderr)
        return 1
    with open(path, "rb") as fh:
        head = fh.read(110)
    print(f"file    : {path}")
    print(f"header  : {head[:108].decode('ascii', 'replace').strip()}")
    print(f"version : {head[108]}   record_size: {head[109]}")
    print(f"colour  : {'(B-V)*50' if head[108] == 2 else '(BP-RP)*10 -> converted'}")

    arr, version = read_area_file(path, mag_limit=None)
    print(f"\ndecoded : {len(arr)} stars")
    print(f"  RA    {arr['ra'].min():.3f} .. {arr['ra'].max():.3f}")
    print(f"  Dec   {arr['dec'].min():.3f} .. {arr['dec'].max():.3f}")
    print(f"  V     {arr['V'].min():.1f} .. {arr['V'].max():.1f}")
    known = ~np.isnan(arr["BV"])
    if known.any():
        print(f"  B-V   {arr['BV'][known].min():+.2f} .. {arr['BV'][known].max():+.2f}"
              f"   ({100*known.mean():.1f}% known)")

    # geometry check: the decoded extent must match the area's declared ring/cell
    ring = int(area[:2]); cell = int(area[2:])
    n = CELLS_PER_RING[ring - 1]
    exp_dec = (DEC_EDGES[ring - 1], DEC_EDGES[ring])
    exp_ra = ((cell - 1) * 360.0 / n, cell * 360.0 / n)
    ok_geom = (arr["dec"].min() >= exp_dec[0] - 1e-3 and arr["dec"].max() <= exp_dec[1] + 1e-3
               and arr["ra"].min() >= exp_ra[0] - 1e-3 and arr["ra"].max() <= exp_ra[1] + 1e-3)
    print(f"\ngeometry: expected Dec {exp_dec[0]:.3f}..{exp_dec[1]:.3f}, "
          f"RA {exp_ra[0]:.3f}..{exp_ra[1]:.3f}  -> {'OK' if ok_geom else 'MISMATCH'}")

    cat = Tycho2()
    t = cat.box(arr["ra"].min(), arr["ra"].max(), arr["dec"].min(), arr["dec"].max())
    t = t[(t["V"] >= mag_lo) & (t["V"] <= mag_hi) & (~np.isnan(t["BV"]))]
    if len(t) < 20:
        print("not enough Tycho-2 stars here to cross-check", file=sys.stderr)
        return 0 if ok_geom else 1

    tree = cKDTree(_radec_to_xyz(arr["ra"], arr["dec"]))
    chord, i = tree.query(_radec_to_xyz(t["ra"], t["dec"]), k=1)
    sep = _chord_to_arcsec(chord)
    m = sep < 1.0
    dV = arr["V"][i[m]] - t["V"][m]
    dBV = arr["BV"][i[m]] - t["BV"][m]

    def robust(x):
        x = x[np.isfinite(x)]
        return 1.4826 * np.median(np.abs(x - np.median(x)))

    print(f"\ncross-check vs Tycho-2 ({mag_lo}<=V<={mag_hi}): "
          f"{m.sum()} of {len(t)} matched within 1\"")
    print(f"  median separation {np.median(sep[m]):.2f}\"")
    print(f"  V    median {np.median(dV):+.3f}   robust sigma {robust(dV):.3f}")
    print(f"  B-V  median {np.nanmedian(dBV):+.3f}   robust sigma {robust(dBV):.3f}")

    ok_astro = np.median(sep[m]) < 2.0 and m.sum() > 0.8 * len(t)
    ok_phot = abs(np.median(dV)) < 0.15
    print(f"\nastrometry {'OK' if ok_astro else 'FAIL'}   photometry "
          f"{'OK' if ok_phot else 'FAIL'}")
    return 0 if (ok_geom and ok_astro and ok_phot) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-path", default=str(DEFAULT_DB_PATH))
    ap.add_argument("--db-name", default=DEFAULT_DB_NAME, help="v50, d80, ...")
    ap.add_argument("--selftest", action="store_true",
                     help="decode one area and cross-check against Tycho-2")
    ap.add_argument("--area", default="0720", help="RRCC area for --selftest (default: RCen's)")
    ap.add_argument("--mag-lo", type=float, default=7.0)
    ap.add_argument("--mag-hi", type=float, default=11.0)
    ap.add_argument("--cone", nargs=3, type=float, metavar=("RA", "DEC", "RADIUS_DEG"),
                     help="list stars around a position")
    ap.add_argument("--mag-limit", type=float, default=None)
    args = ap.parse_args()

    if args.selftest:
        sys.exit(_selftest(args.db_path, args.db_name, args.area, args.mag_lo, args.mag_hi))

    if args.cone:
        ra, dec, rad = args.cone
        cat = GaiaV50(args.db_path, args.db_name, mag_limit=args.mag_limit)
        sub = cat.cone(ra, dec, rad)
        sub = np.sort(sub, order="V")
        print(f"# {len(sub)} stars within {rad} deg of {ra}, {dec}  "
              f"(area {area_name(ra, dec)}, db {args.db_name} v{cat.version})")
        print(f"{'ra_deg':>11} {'dec_deg':>11} {'V':>6} {'B-V':>7}")
        for r in sub[:200]:
            bv = "  --  " if np.isnan(r["BV"]) else f"{r['BV']:+.2f}"
            print(f"{r['ra']:11.5f} {r['dec']:11.5f} {r['V']:6.1f} {bv:>7}")
        return

    ap.print_help()


if __name__ == "__main__":
    main()
