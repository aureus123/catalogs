#!/usr/bin/env python3
"""Gaia DR3 photometry fetched live from VizieR, as a third reference backend.

This mirrors what ASTAP's GUI does in ``unit_online_gaia.pas``: query VizieR's ASU-text
interface for Gaia DR3 (catalogue I/355), then convert G/BP/RP to Johnson V and B-V with
the DR3 photometric relations. Same URL construction as ``read_stars_online``
(unit_online_gaia.pas:258) and the same polynomials as ``transform_gaia``
(unit_online_gaia.pas:33-109).

Why it might beat the local V50 database:

* **No 0.1 mag quantisation.** VizieR returns Gmag/BPmag/RPmag to 6 decimals, whereas
  V50 stores V in 0.1 mag steps (a 0.029 mag uniform quantisation floor).
* **No B-V cap.** V50's signed colour byte clips at +-2.54, so carbon stars come back as
  unknown. Here B-V is computed on the fly and only limited by the polynomial's validity.
* **Deeper on demand** -- the magnitude limit is a query parameter.

The cost is an internet round trip per field, which is exactly the trade-off this module
exists to measure. Responses are cached under ``estim/cache/online/`` so a field is only
fetched once.

Caveat on epoch: VizieR I/355 positions are at **epoch 2016.0**, while V50's are
propagated to 2030 and Tycho-2's are 2000. Proper motion is fetched and applied here to
bring positions to the image epoch -- ASTAP itself does not do this, so this backend is
slightly *better* positioned than ASTAP's own online path.
"""
import argparse
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CACHE_DIR = HERE / "cache" / "online"
VIZIER_BASE = "http://vizier.u-strasbg.fr/viz-bin/asu-txt"

DTYPE = np.dtype([
    ("ra", "f8"), ("dec", "f8"),
    ("V", "f4"), ("BV", "f4"), ("e_V", "f4"),
    ("G", "f4"), ("BP", "f4"), ("RP", "f4"),
])

# Gaia DR3 photometric relations, transcribed from unit_online_gaia.pas:62-89.
# https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/
#   cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html
_V_COEF = (0.02704, -0.01424, 0.2156, -0.01426)                      # valid -0.5 <= BP-RP <= 5.0
_VT_COEF = (0.01077, 0.0682, 0.2387, -0.02342)                       # valid -0.3 < BP-RP < 3.0
_BT_COEF = (0.004288, 0.8547, -0.1244, 0.9085, -0.4843, 0.06814)

# The DR3 G->V relation has ~0.03-0.05 mag intrinsic scatter (it depends on metallicity
# and reddening, not colour alone), which is the real floor on this backend.
TRANSFORM_SIGMA = 0.046


def transform_gaia(G, BP, RP):
    """G/BP/RP -> (V, B-V), following ASTAP exactly including its quality guard.

    Returns NaN where the transform is not applicable, rather than extrapolating.
    """
    G = np.asarray(G, float); BP = np.asarray(BP, float); RP = np.asarray(RP, float)
    ok = np.isfinite(G) & np.isfinite(BP) & np.isfinite(RP) & (G != 0) & (BP != 0) & (RP != 0)

    # Straylight/blending guard (unit_online_gaia.pas:55-59): if the BP+RP flux is far
    # above G's, BP/RP are unreliable. ASTAP rejects when c > 4 AND G > BP.
    with np.errstate(over="ignore", invalid="ignore"):
        Gf = np.power(10.0, (22.0 - G) / 2.5)
        BPf = np.power(10.0, (22.0 - BP) / 2.5)
        RPf = np.power(10.0, (22.0 - RP) / 2.5)
        c = (BPf + RPf) / Gf
    ok &= ~((c > 4) & (G > BP))

    x = np.where(ok, BP - RP, np.nan)

    v_ok = ok & (x >= -0.5) & (x <= 5.0)
    V = np.where(v_ok, G + _V_COEF[0] + _V_COEF[1]*x + _V_COEF[2]*x**2 + _V_COEF[3]*x**3, np.nan)

    # B = V + 0.850*(Bt - Vt); the G term cancels inside the difference.
    bv_ok = ok & (x > -0.3) & (x < 3.0)
    Vt = _VT_COEF[0] + _VT_COEF[1]*x + _VT_COEF[2]*x**2 + _VT_COEF[3]*x**3
    Bt = (_BT_COEF[0] + _BT_COEF[1]*x + _BT_COEF[2]*x**2 + _BT_COEF[3]*x**3
          + _BT_COEF[4]*x**4 + _BT_COEF[5]*x**5)
    BV = np.where(bv_ok, 0.850 * (Bt - Vt), np.nan)
    return V, BV


def _build_url(ra_deg, dec_deg, box_w_arcsec, box_h_arcsec, mag_limit, with_pm=True):
    """VizieR ASU-text query, same shape as ASTAP's (unit_online_gaia.pas:258) except
    that ``-c.bs`` is given width/height separately instead of a square. The Dwarf field
    is 3.0 x 1.7 deg, so a square box asks for ~2.4x more sky than needed -- and VizieR's
    response time scales with the area searched, not the rows returned."""
    sgn = "%2B" if dec_deg >= 0 else "%2D"
    cols = "RA_ICRS,DE_ICRS,Gmag,BPmag,RPmag"
    if with_pm:
        cols += ",pmRA,pmDE"
    return (f"{VIZIER_BASE}?-source=I/355/Gaiadr3&-out={cols}"
            f"&-c={abs(ra_deg):.10f}{sgn}{abs(dec_deg):.10f}"
            f"&-c.bs={box_w_arcsec:.6f}/{box_h_arcsec:.6f}"
            f"&-out.max=200000&BPmag=<{mag_limit:.2f}")


def _parse_asu_txt(text, with_pm=True):
    """Parse VizieR's fixed-width ASU-text body: skip '#' comments, the column header,
    the units line and the '---' separator, then read whitespace-separated columns."""
    rows = []
    seen_sep = False
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        s = line.strip()
        if not s:
            continue
        if set(s) <= set("- "):
            seen_sep = True
            continue
        if not seen_sep:
            continue
        parts = s.split()
        need = 7 if with_pm else 5
        if len(parts) < 5:
            continue
        try:
            ra = float(parts[0]); dec = float(parts[1])
        except ValueError:
            continue

        def num(i):
            if i >= len(parts):
                return np.nan
            try:
                return float(parts[i])
            except ValueError:
                return np.nan

        rows.append((ra, dec, num(2), num(3), num(4),
                     num(5) if with_pm else 0.0, num(6) if with_pm else 0.0))
    return rows


def fetch(ra_deg, dec_deg, box_w_deg, box_h_deg=None, mag_limit=16.0, epoch=None,
          cache_dir=CACHE_DIR, timeout=300, force=False, verbose=True):
    """Fetch a field from VizieR (cached) and return a DTYPE array.

    ``epoch`` (decimal year) propagates positions from Gaia's 2016.0 using pmRA/pmDE.
    ``box_h_deg`` defaults to ``box_w_deg`` (square, as ASTAP does).
    """
    if box_h_deg is None:
        box_h_deg = box_w_deg
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = (f"{ra_deg:.5f}_{dec_deg:+.5f}_{box_w_deg:.4f}x{box_h_deg:.4f}"
           f"_{mag_limit:.2f}.txt")
    cache_file = cache_dir / key

    if cache_file.exists() and not force:
        text = cache_file.read_text()
        if verbose:
            print(f"  [cache] {cache_file.name}", file=sys.stderr)
    else:
        url = _build_url(ra_deg, dec_deg, box_w_deg * 3600.0, box_h_deg * 3600.0, mag_limit)
        if verbose:
            print(f"  [vizier] BP<{mag_limit} over {box_w_deg:.2f}x{box_h_deg:.2f} deg ...",
                  file=sys.stderr)
        t0 = time.time()
        # A precalibration run is ~100 sequential queries over several hours, so a
        # transient blip must not abort it. Retry with backoff, then give up cleanly.
        last = None
        text = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=timeout) as r:
                    text = r.read().decode("ascii", "replace")
                break
            except (urllib.error.URLError, TimeoutError, OSError) as e:
                last = e
                if attempt < 2:
                    wait = 10 * (attempt + 1)
                    if verbose:
                        print(f"  [vizier] attempt {attempt+1} failed ({e}); "
                              f"retrying in {wait}s", file=sys.stderr)
                    time.sleep(wait)
        if text is None:
            raise RuntimeError(
                f"VizieR query failed after 3 attempts ({last}); no internet connection?")
        if "<html" in text[:200].lower():
            raise RuntimeError("VizieR returned an HTML error page, not data")
        cache_file.write_text(text)
        if verbose:
            print(f"  [vizier] {len(text)} bytes in {time.time()-t0:.1f}s -> {cache_file.name}",
                  file=sys.stderr)

    rows = _parse_asu_txt(text)
    if not rows:
        return np.empty(0, dtype=DTYPE)
    a = np.array(rows, dtype=float)
    ra, dec, G, BP, RP, pmra, pmde = (a[:, i] for i in range(7))

    if epoch is not None:
        dt = epoch - 2016.0
        pmra = np.nan_to_num(pmra); pmde = np.nan_to_num(pmde)
        # pmRA from VizieR is already mu_alpha* (includes cos(dec))
        ra = ra + (pmra / 3.6e6) * dt / np.maximum(np.cos(np.radians(dec)), 1e-6)
        dec = dec + (pmde / 3.6e6) * dt

    V, BV = transform_gaia(G, BP, RP)
    out = np.empty(len(a), dtype=DTYPE)
    out["ra"], out["dec"] = ra, dec
    out["V"], out["BV"] = V, BV
    out["e_V"] = TRANSFORM_SIGMA
    out["G"], out["BP"], out["RP"] = G, BP, RP
    return out[np.isfinite(out["V"])]


def _radec_to_xyz(ra_deg, dec_deg):
    ra = np.radians(ra_deg); dec = np.radians(dec_deg)
    cd = np.cos(dec)
    return np.stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)], axis=-1)


def _chord_to_arcsec(chord):
    return np.degrees(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))) * 3600.0


class GaiaOnline:
    """Same surface as ``tycho2.Tycho2`` / ``gaia_v50.GaiaV50``, backed by VizieR.

    Unlike the local backends this needs the field up front, because the query is
    per-field. Either pass ``preload=(ra, dec, box_deg)`` or just call ``match()`` --
    it fetches whatever the query region needs on first use.
    """

    def __init__(self, mag_limit=16.0, epoch=None, cache_dir=CACHE_DIR, verbose=True):
        self.mag_limit = mag_limit
        self.epoch = epoch
        self.cache_dir = Path(cache_dir)
        self.verbose = verbose
        self.data = np.empty(0, dtype=DTYPE)
        self._regions = []
        self._seen = set()

    def __len__(self):
        return len(self.data)

    def preload(self, ra_deg, dec_deg, box_w_deg, box_h_deg=None):
        if box_h_deg is None:
            box_h_deg = box_w_deg
        for (r, d, w, h) in self._regions:
            if (abs(r - ra_deg) < 1e-6 and abs(d - dec_deg) < 1e-6
                    and w >= box_w_deg - 1e-9 and h >= box_h_deg - 1e-9):
                return
        arr = fetch(ra_deg, dec_deg, box_w_deg, box_h_deg, self.mag_limit, self.epoch,
                    self.cache_dir, verbose=self.verbose)
        arr = self._drop_duplicates(arr)
        self.data = np.concatenate([self.data, arr]) if len(self.data) else arr
        self._regions.append((ra_deg, dec_deg, box_w_deg, box_h_deg))

    def _drop_duplicates(self, arr):
        """Discard stars already held from an earlier query.

        Overlapping fields (ROct x6, SOct x3, YCen x2, SVir x2, and any incidental
        overlap) return the *same* Gaia sources. Without this, ``self.data`` accumulates
        exact duplicates at 0.000" separation, and ``starcat.blend_dmag`` then sees every
        star "blended" with its own copy -- which silently rejected every calibrator in
        the affected fields (ROct5 and ROct6 yielded zero).

        Positions are bit-identical between queries for the same source at the same
        epoch, so rounding to 1e-6 deg (3.6 mas) identifies duplicates exactly without
        ever merging two genuinely distinct stars.
        """
        if len(arr) == 0:
            return arr
        keys = list(zip(np.round(arr["ra"], 6).tolist(),
                        np.round(arr["dec"], 6).tolist()))
        fresh = np.fromiter((k not in self._seen for k in keys), bool, len(keys))
        self._seen.update(k for k, ok in zip(keys, fresh) if ok)
        return arr[fresh]

    def box(self, ra_min, ra_max, dec_min, dec_max):
        cra = (ra_min + ra_max) / 2.0
        cdec = (dec_min + dec_max) / 2.0
        w = (ra_max - ra_min) * np.cos(np.radians(cdec)) * 1.1
        h = (dec_max - dec_min) * 1.1
        self.preload(cra, cdec, w, h)
        d = self.data
        if len(d) == 0:
            return d
        m = ((d["dec"] >= dec_min) & (d["dec"] <= dec_max)
             & (d["ra"] >= ra_min) & (d["ra"] <= ra_max))
        return d[m]

    def cone(self, ra_deg, dec_deg, radius_deg):
        self.preload(ra_deg, dec_deg, radius_deg * 2.2)
        if len(self.data) == 0:
            return self.data
        center = _radec_to_xyz(np.array([ra_deg]), np.array([dec_deg]))[0]
        sep = _chord_to_arcsec(
            np.linalg.norm(_radec_to_xyz(self.data["ra"], self.data["dec"]) - center, axis=1))
        return self.data[sep <= radius_deg * 3600.0]

    def match(self, ra_deg, dec_deg, radius_arcsec=8.0, isolation_arcsec=20.0):
        ra_deg = np.atleast_1d(np.asarray(ra_deg, dtype="f8"))
        dec_deg = np.atleast_1d(np.asarray(dec_deg, dtype="f8"))
        n = len(ra_deg)
        idx_out = np.full(n, -1, dtype="i8")
        sep_out = np.full(n, np.nan)
        nn2_out = np.full(n, np.nan)
        if n == 0:
            return idx_out, sep_out, nn2_out

        cdec = float(np.mean(dec_deg))
        w = (ra_deg.max() - ra_deg.min()) * np.cos(np.radians(cdec)) * 1.15 + 0.05
        h = (dec_deg.max() - dec_deg.min()) * 1.15 + 0.05
        self.preload(float(np.mean(ra_deg)), cdec, w, h)
        if len(self.data) == 0:
            return idx_out, sep_out, nn2_out

        cat_xyz = _radec_to_xyz(self.data["ra"], self.data["dec"])
        tree = cKDTree(cat_xyz)
        if len(self.data) > 1:
            ch, _ = tree.query(cat_xyz, k=2)
            cat_nn = _chord_to_arcsec(ch[:, 1])
        else:
            cat_nn = np.full(len(self.data), np.inf)
        chord, ii = tree.query(_radec_to_xyz(ra_deg, dec_deg), k=1)
        sep = _chord_to_arcsec(chord)
        hit = sep <= radius_arcsec
        idx_out[hit] = ii[hit]
        sep_out[hit] = sep[hit]
        nn2_out[hit] = cat_nn[ii[hit]]
        return idx_out, sep_out, nn2_out

    def designation(self, i):
        rec = self.data[i]
        return f"GaiaDR3 J{rec['ra']:09.5f}{rec['dec']:+09.5f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ra", type=float)
    ap.add_argument("dec", type=float)
    ap.add_argument("box", type=float, help="square box size in degrees")
    ap.add_argument("--mag-limit", type=float, default=16.0)
    ap.add_argument("--epoch", type=float, default=None)
    ap.add_argument("--force", action="store_true", help="ignore the cache")
    args = ap.parse_args()

    arr = fetch(args.ra, args.dec, args.box, args.box, args.mag_limit, args.epoch,
                force=args.force)
    print(f"{len(arr)} stars with a usable V")
    if len(arr):
        bv = arr["BV"][np.isfinite(arr["BV"])]
        print(f"  V    {arr['V'].min():.2f} .. {arr['V'].max():.2f}")
        print(f"  B-V  {bv.min():+.2f} .. {bv.max():+.2f}  ({100*len(bv)/len(arr):.1f}% known)")
        for r in np.sort(arr, order="V")[:5]:
            print(f"   {r['ra']:11.5f} {r['dec']:+11.5f}  V={r['V']:6.2f}  B-V={r['BV']:+.2f}"
                  f"  G={r['G']:.3f} BP={r['BP']:.3f} RP={r['RP']:.3f}")


if __name__ == "__main__":
    main()
