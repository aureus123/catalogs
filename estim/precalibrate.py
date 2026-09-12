#!/usr/bin/env python3
"""Pre-calibration: fit per-scope colour/magnitude transformation coefficients
from many Dwarf 2 / Dwarf 3 fields against Tycho-2, and write
instrument_coeffs.json (schema_version 2, see the design plan).

Fields are discovered under ARCHIVE_ROOT (see below), which holds one subfolder per
scope. The archive is treated as **read-only**: every FITS write (extracted planes,
ASTAP csvs) goes to estim/cache/, never next to the original files.

    python3 precalibrate.py --filter VIS --gain 60 -o instrument_coeffs.json
"""
import argparse
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits

import photometry
import starcat
from starcat import blend_dmag, e_V, get_catalog, isolated


def _astap_version():
    # ASTAP's CLI prints no version string and its Info.plist bundle version
    # (CFBundleShortVersionString "0.1") is a generic Lazarus default, not the
    # actual release -- there is no reliable way to read it programmatically.
    return "unknown"


def _photutils_version():
    try:
        import photutils
        return photutils.__version__
    except Exception:
        return "unknown"
from tycho2 import DEFAULT_NPY as TYCHO2_NPY
from tycho2 import Tycho2

HERE = Path(__file__).resolve().parent
#: Root of the image archive. It is expected to hold one subfolder per scope --
#: ``Dwarf2/`` and ``Dwarf3/`` -- and inside each, one folder per capture named after the
#: field, containing ``<field>.fits`` (the 3-plane cube), usually a ``<field>.wcs`` and a
#: ``shotsInfo.json``.
#:
#: The archive is treated as **strictly read-only**: nothing is ever written back into
#: it. Extracted planes and ASTAP detections go to ``estim/cache/<field>/`` instead, which
#: is what precalibrate.py and benchmark.py reuse between runs.
#:
#: Override with the ``DWARF_ARCHIVE`` environment variable if the archive lives
#: elsewhere.
ARCHIVE_ROOT = Path(os.environ.get("DWARF_ARCHIVE", Path.home() / "Desktop"))
SCOPE_DIRS = ("Dwarf2", "Dwarf3")
CACHE_DIR = HERE / "cache"

CALIB_V_MIN, CALIB_V_MAX = 7.0, 13.0
CALIB_EVT_MAX = 0.10
TYCHO2_MATCH_ARCSEC = 8.0
TYCHO2_ISOLATION_ARCSEC = 20.0
OUR_ISOLATION_ARCSEC = 20.0
SNR_MIN_CALIB = 5.0

# Positions are propagated to roughly the middle of the image archive. Gaia DR3 is at
# epoch 2016.0 and the archive spans 2025-2026, so a single epoch costs at most ~0.1"
# of proper motion for a typical star -- negligible against the 8" match radius.
DEFAULT_EPOCH = 2025.5

V_BANDS = [(7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]
COLOUR_CANDIDATES = ("b-r", "b-g", "g-r")


def _find_cube_fits(folder):
    """Most folders hold the cube at <folder>/<folder>.fits, but a couple of
    dozen older captures (pre-dating that naming convention, e.g. SVir2/SVir.fits,
    V462Lup1/V462-Lup1.fits) don't. Fall back to the single top-level 3-plane
    FITS file if there is exactly one -- folders with many top-level FITS
    (Darks_*/Flats_* calibration-frame dumps) stay excluded since they're
    ambiguous, which is also the correct call for them."""
    expected = folder / f"{folder.name}.fits"
    if expected.exists():
        return expected
    cube_candidates = []
    for p in folder.glob("*.fits"):
        if not p.is_file():
            continue
        try:
            h = fits.getheader(p)
        except Exception:
            continue
        if h.get("NAXIS") == 3 and h.get("NAXIS3") == 3:
            cube_candidates.append(p)
    return cube_candidates[0] if len(cube_candidates) == 1 else None


def discover_fields(filt, gain):
    """Yield (scope, field_name, fits_path) for Dwarf2/Dwarf3 cubes matching
    filt/gain. Only ever reads headers -- never writes into the archive."""
    for root in (ARCHIVE_ROOT / d for d in SCOPE_DIRS):
        if not root.is_dir():
            continue
        for folder in sorted(root.iterdir()):
            if not folder.is_dir():
                continue
            fits_path = _find_cube_fits(folder)
            if fits_path is None:
                continue
            try:
                h = fits.getheader(fits_path)
            except Exception:
                continue
            if h.get("FILTER") != filt:
                continue
            if gain is not None and h.get("GAIN") != gain:
                continue
            if h.get("NAXIS") != 3 or h.get("NAXIS3") != 3:
                continue
            yield photometry.normalize_scope(h.get("TELESCOP")), folder.name, fits_path


def instrumental_colour(name, g, r, b):
    return {"b-r": b - r, "b-g": b - g, "g-r": g - r}[name]


def within_center(values, group):
    out = np.empty(len(values), dtype=float)
    for gid in np.unique(group):
        m = group == gid
        out[m] = values[m] - np.mean(values[m])
    return out


def fit_slope_through_origin(x, y):
    """y = slope * x + eps (both already group-centered). Returns slope, se(slope), resid std."""
    sxx = np.sum(x * x)
    slope = np.sum(x * y) / sxx
    resid = y - slope * x
    dof = max(len(x) - 1, 1)
    resid_std = np.sqrt(np.sum(resid ** 2) / dof)
    se_slope = resid_std / np.sqrt(sxx)
    return slope, se_slope, resid_std, resid


def fit_two_covariate(x1, x2, y):
    """y = a*x1 + b*x2 + eps (all group-centered). Returns (a, b, cov, resid)."""
    X = np.column_stack([x1, x2])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - 2, 1)
    sigma2 = np.sum(resid ** 2) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    return beta[0], beta[1], cov, resid


def measure_one_field(cube_path, cache_root, extraction):
    field_cache = cache_root / cube_path.parent.name
    field_cache.mkdir(parents=True, exist_ok=True)
    return photometry.measure_field(
        cube_path, workdir=field_cache,
        snr_min=extraction["snr_min"], max_stars=extraction["max_stars"],
        aperture_hfd=extraction["aperture_hfd"],
        annulus_in_hfd=extraction["annulus_in_hfd"],
        annulus_out_hfd=extraction["annulus_out_hfd"],
        backend=extraction["backend"], saturation_adu=extraction["saturation_adu"],
        channel_radius_px=extraction["channel_radius_px"],
    )


def collect_calibrators(fp, cat, field_id, image_size):
    idx, sep, nn2 = cat.match(fp.ra, fp.dec, radius_arcsec=TYCHO2_MATCH_ARCSEC,
                               isolation_arcsec=TYCHO2_ISOLATION_ARCSEC)
    matched = idx >= 0
    if not np.any(matched):
        return None
    d = cat.data[idx[matched]]
    mi = np.nonzero(matched)[0]

    with np.errstate(invalid="ignore"):
        flux_ok = (fp.flux_r[mi] > 0) & (fp.flux_g[mi] > 0) & (fp.flux_b[mi] > 0) & \
                  np.isfinite(fp.flux_r[mi]) & np.isfinite(fp.flux_g[mi]) & np.isfinite(fp.flux_b[mi])
    # Blend rejection is a flux-ratio test, not "any neighbour within 20 arcsec".
    # The old rule got stricter as the catalogue deepened (on RCen it cut 7 stars with
    # Tycho-2 but 117 with V50 unlimited), so a better reference perversely shrank the
    # calibrator sample. See starcat.blend_dmag.
    iso = isolated(cat.data, idx, radius_arcsec=TYCHO2_ISOLATION_ARCSEC,
                   dmag=starcat.DEFAULT_BLEND_DMAG)
    ok = (
        (d["V"] >= CALIB_V_MIN) & (d["V"] <= CALIB_V_MAX) &
        (e_V(d) < CALIB_EVT_MAX) & (~np.isnan(d["BV"])) &
        iso[mi] &
        (fp.isolation_arcsec[mi] > OUR_ISOLATION_ARCSEC) &
        (~fp.saturated[mi]) &
        (fp.snr_g[mi] >= SNR_MIN_CALIB) &
        flux_ok
    )
    sel = mi[ok]
    if len(sel) == 0:
        return None
    d = d[ok]

    nx, ny = image_size
    cx, cy = nx / 2.0, ny / 2.0
    r2_norm = ((fp.x[sel] - cx) ** 2 + (fp.y[sel] - cy) ** 2) / (cx ** 2 + cy ** 2)

    return dict(
        g=-2.5 * np.log10(fp.flux_g[sel]),
        r=-2.5 * np.log10(fp.flux_r[sel]),
        b=-2.5 * np.log10(fp.flux_b[sel]),
        V_cat=d["V"].astype(float), BV_cat=d["BV"].astype(float),
        e_VT=np.asarray(e_V(d), dtype=float),
        r2=r2_norm,
        field=np.full(len(sel), field_id),
    )


def fit_scope(scope, fields, cat, extraction, cache_root, verbose):
    tables = []
    used_fields = []
    for i, (name, cube_path) in enumerate(fields):
        # The catalogue lookup is inside the try as well: with the online backend it
        # hits the network, and a single blip must skip one field rather than abort a
        # multi-hour run.
        try:
            fp = measure_one_field(cube_path, cache_root, extraction)
            h = fits.getheader(cube_path)
            tab = collect_calibrators(fp, cat, i, (h["NAXIS1"], h["NAXIS2"]))
        except Exception as e:
            print(f"  [{scope}] {name}: FAILED ({e})", file=sys.stderr)
            if verbose:
                traceback.print_exc()
            continue
        if tab is None or len(tab["g"]) < 5:
            if verbose:
                n = 0 if tab is None else len(tab["g"])
                print(f"  [{scope}] {name}: only {n} calibrators, skipped", file=sys.stderr)
            continue
        tables.append(tab)
        used_fields.append(name)
        if verbose:
            print(f"  [{scope}] {name}: {len(tab['g'])} calibrators", file=sys.stderr)

    if not tables:
        return None

    pooled = {k: np.concatenate([t[k] for t in tables]) for k in tables[0]}
    field_id = pooled["field"]
    n_fields = len(used_fields)
    n_stars = len(pooled["g"])

    # --- colour index selection: pooled fixed-effects fit of instrumental
    # colour vs catalogue B-V, one global slope with a per-field intercept ---
    bv_c = within_center(pooled["BV_cat"], field_id)
    candidates = {}
    best_name, best_sigma_bv = None, np.inf
    best_fit = None
    for name in COLOUR_CANDIDATES:
        c = instrumental_colour(name, pooled["g"], pooled["r"], pooled["b"])
        c_c = within_center(c, field_id)
        slope, se_slope, resid_std, resid = fit_slope_through_origin(bv_c, c_c)
        tbv = 1.0 / slope
        sigma_bv = resid_std * abs(tbv)
        candidates[name] = round(float(tbv), 3)
        if sigma_bv < best_sigma_bv:
            best_sigma_bv = sigma_bv
            best_name = name
            se_tbv = (tbv ** 2) * se_slope
            best_fit = dict(tbv=tbv, se_tbv=se_tbv, sigma_bv=sigma_bv, resid=resid, c=c)

    colour_index = best_name
    Tbv = best_fit["tbv"]

    # --- magnitude transform: (V_cat - g) = zp_field + Tv_bv*BV_cat + k*r2 ---
    y = pooled["V_cat"] - pooled["g"]
    y_c = within_center(y, field_id)
    r2_c = within_center(pooled["r2"], field_id)
    Tv_bv, k_r2, cov, resid_v = fit_two_covariate(bv_c, r2_c, y_c)
    se_tv_bv = np.sqrt(cov[0, 0])
    se_k = np.sqrt(cov[1, 1])

    # --- between-field spread of Tv_bv (per-field OLS vs the pooled formal error) ---
    per_field_tv = []
    per_field_se = []
    for fid in np.unique(field_id):
        m = field_id == fid
        if m.sum() < 10:
            continue
        bv = pooled["BV_cat"][m] - np.mean(pooled["BV_cat"][m])
        yy = y[m] - np.mean(y[m])
        s, se, _, _ = fit_slope_through_origin(bv, yy)
        per_field_tv.append(s)
        per_field_se.append(se)
    if len(per_field_tv) >= 3:
        per_field_tv = np.array(per_field_tv)
        per_field_se = np.array(per_field_se)
        raw_var = np.var(per_field_tv, ddof=1)
        mean_formal_var = np.mean(per_field_se ** 2)
        between_field_sigma = float(np.sqrt(max(raw_var - mean_formal_var, 0.0)))
    else:
        between_field_sigma = None

    # --- error bands, Tycho-2-deducted ---
    v_est_resid = resid_v  # V_cat - g - zp_field - Tv_bv*BV - k*r2, group-centered
    bv_est_resid = best_fit["resid"] * abs(Tbv)  # colour-space residual -> BV-space
    sigma_v_by_band = {}
    sigma_bv_by_band = {}
    for lo, hi in V_BANDS:
        m = (pooled["V_cat"] >= lo) & (pooled["V_cat"] < hi)
        label = f"{lo}-{hi}"
        if m.sum() >= 5:
            rms_v = np.sqrt(np.mean(v_est_resid[m] ** 2))
            med_evt = np.median(pooled["e_VT"][m])
            sigma_v_by_band[label] = round(float(np.sqrt(max(rms_v ** 2 - med_evt ** 2, 0.0))), 3)
            rms_bv = np.sqrt(np.mean(bv_est_resid[m] ** 2))
            sigma_bv_by_band[label] = round(float(np.sqrt(max(rms_bv ** 2 - (med_evt * abs(Tbv)) ** 2, 0.0))), 3)

    valid_bv_range = [0.0, round(float(np.nanpercentile(pooled["BV_cat"], 99)), 2)]

    return {
        "scope": scope, "filter": extraction["filter"], "gain": extraction["gain"],
        # Per-set, not top level: different scopes may legitimately be precalibrated
        # against different catalogues (e.g. Dwarf 3 online, Dwarf 2 offline from V50),
        # and measure.py must be able to pick the right one per scope.
        "reference_catalog": {
            "name": extraction["catalog"],
            "system": starcat.PHOTOMETRIC_SYSTEM[extraction["catalog"]],
            "mag_limit": extraction.get("cat_mag_limit"),
            "epoch": extraction.get("epoch"),
            "requires_network": extraction["catalog"] in starcat.NEEDS_NETWORK,
            "source": {
                "tycho2": "cat/tyc2.txt",
                "gaia_v50": "ASTAP v50 local database (.1476)",
                "gaia_online": "VizieR I/355/Gaiadr3 via asu-txt",
            }[extraction["catalog"]],
        },
        "n_fields": n_fields, "n_stars": int(n_stars),
        "extraction": {
            "detection": {
                "tool": "ASTAP -extract", "astap_version": extraction.get("astap_version", "unknown"),
                "snr_min": extraction["snr_min"], "max_stars": extraction["max_stars"],
                "min_star_size_hfd": extraction.get("min_star_size_hfd", 0.8),
            },
            "photometry": {
                "backend": extraction["backend"], "photutils_version": extraction.get("photutils_version", "unknown"),
                "units": "multiples of the per-star measured HFD",
                "aperture_hfd": extraction["aperture_hfd"],
                "annulus_in_hfd": extraction["annulus_in_hfd"],
                "annulus_out_hfd": extraction["annulus_out_hfd"],
                "method": "exact", "background": "sigma_clipped_median", "sigma_clip": 3.0,
            },
            "saturation_adu": extraction["saturation_adu"],
            "match": {
                "channel_radius_px": extraction["channel_radius_px"],
                "tycho2_radius_arcsec": TYCHO2_MATCH_ARCSEC,
                "isolation_arcsec": TYCHO2_ISOLATION_ARCSEC,
            },
        },
        "colour_index": colour_index,
        "colour_index_candidates": candidates,
        "coefficients": {
            "Tv_bv": {"value": round(float(Tv_bv), 3), "sigma": round(float(se_tv_bv), 3),
                      **({"between_field_sigma": round(between_field_sigma, 3)}
                         if between_field_sigma is not None else {})},
            "Tbv": {"value": round(float(Tbv), 3), "sigma": round(float(best_fit["se_tbv"]), 3)},
            "k_r2": {"value": round(float(k_r2), 3), "sigma": round(float(se_k), 3)},
        },
        "valid_bv_range": valid_bv_range,
        "sigma_V_by_band": sigma_v_by_band,
        "sigma_BV_by_band": sigma_bv_by_band,
        "sweep": extraction.get("sweep_provenance", {
            "objective": "not run -- extraction parameters are the plan's defaults",
            "grid": None, "best_score": None, "astap_baseline_score": None,
        }),
        "fields": used_fields,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--filter", default="VIS")
    ap.add_argument("--gain", type=int, default=60)
    ap.add_argument("--scope", default=None,
                     help="restrict to one TELESCOP value, e.g. 'DWARF 3' or 'DWARF II'")
    ap.add_argument("-o", "--out", default=str(HERE / "instrument_coeffs.json"))
    ap.add_argument("--catalog", choices=starcat.CATALOG_NAMES, default="tycho2",
                     help="reference catalogue. 'gaia_online' REQUIRES AN INTERNET "
                          "CONNECTION (queries VizieR per field, then caches).")
    ap.add_argument("--cat-mag-limit", type=float, default=14.0,
                     help="depth for the Gaia backends (BP for online, V for v50)")
    ap.add_argument("--epoch", type=float, default=DEFAULT_EPOCH,
                     help="epoch to propagate Gaia positions to")
    ap.add_argument("--tycho2", default=str(TYCHO2_NPY))
    ap.add_argument("--max-fields", type=int, default=None,
                     help="limit fields per scope (for quick testing)")
    ap.add_argument("--snr-min", type=float, default=5)
    ap.add_argument("--max-stars", type=int, default=500)
    ap.add_argument("--aperture-hfd", type=float, default=2.0)
    ap.add_argument("--annulus-in-hfd", type=float, default=4.0)
    ap.add_argument("--annulus-out-hfd", type=float, default=6.0)
    ap.add_argument("--backend", choices=["photutils", "astap"], default="photutils")
    ap.add_argument("--saturation-adu", type=float, default=64000)
    ap.add_argument("--channel-radius-px", type=float, default=1.5)
    ap.add_argument("--cache-dir", default=str(CACHE_DIR))
    ap.add_argument("--append", action="store_true",
                     help="merge into an existing -o file instead of overwriting, so "
                          "several (scope, catalogue) schemes coexist and measure.py "
                          "can choose between them")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    extraction = dict(
        filter=args.filter, gain=args.gain, snr_min=args.snr_min, max_stars=args.max_stars,
        aperture_hfd=args.aperture_hfd, annulus_in_hfd=args.annulus_in_hfd,
        annulus_out_hfd=args.annulus_out_hfd, backend=args.backend,
        saturation_adu=args.saturation_adu, channel_radius_px=args.channel_radius_px,
        astap_version=_astap_version(), photutils_version=_photutils_version(),
        catalog=args.catalog, cat_mag_limit=args.cat_mag_limit, epoch=args.epoch,
    )

    if args.catalog in starcat.NEEDS_NETWORK:
        print(f"*** {args.catalog} requires an internet connection: it queries VizieR "
              f"per field.\n*** Responses are cached under estim/cache/online/, so "
              f"re-runs are offline.", file=sys.stderr)
    print(f"loading reference catalogue '{args.catalog}' ...", file=sys.stderr)
    cat = get_catalog(args.catalog, npy_path=args.tycho2, mag_limit=args.cat_mag_limit,
                       epoch=args.epoch, verbose=False)

    print(f"discovering fields (filter={args.filter}, gain={args.gain}) ...", file=sys.stderr)
    by_scope = {}
    for scope, name, path in discover_fields(args.filter, args.gain):
        if args.scope and scope != args.scope:
            continue
        by_scope.setdefault(scope, []).append((name, path))
    if args.scope and not by_scope:
        print(f"error: no fields with TELESCOP={args.scope!r}", file=sys.stderr)
        sys.exit(1)
    for scope in by_scope:
        print(f"  {scope}: {len(by_scope[scope])} candidate fields", file=sys.stderr)
        if args.max_fields:
            by_scope[scope] = by_scope[scope][: args.max_fields]

    cache_root = Path(args.cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)

    sets = []
    t0 = time.time()
    for scope, fields in sorted(by_scope.items()):
        print(f"fitting {scope} ({len(fields)} fields) ...", file=sys.stderr)
        result = fit_scope(scope, fields, cat, extraction, cache_root, args.verbose)
        if result is None:
            print(f"  {scope}: no usable fields, skipped", file=sys.stderr)
            continue
        sets.append(result)
        c = result["coefficients"]
        print(f"  {scope}: colour_index={result['colour_index']} "
              f"Tbv={c['Tbv']['value']}+-{c['Tbv']['sigma']} "
              f"Tv_bv={c['Tv_bv']['value']}+-{c['Tv_bv']['sigma']} "
              f"k_r2={c['k_r2']['value']}+-{c['k_r2']['sigma']} "
              f"n_fields={result['n_fields']} n_stars={result['n_stars']}", file=sys.stderr)

    doc = {
        "schema_version": 4,
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "calibrator_cuts": {
            "V_range": [CALIB_V_MIN, CALIB_V_MAX],
            "e_V_max": CALIB_EVT_MAX,
            "snr_g_min": SNR_MIN_CALIB,
            "blend": {"radius_arcsec": TYCHO2_ISOLATION_ARCSEC,
                      "dmag": starcat.DEFAULT_BLEND_DMAG,
                      "rule": "reject if a catalogue neighbour within radius is brighter "
                              "than V_target + dmag (flux-ratio test, depth-independent)"},
        },
        "sets": sets,
    }

    out_path = Path(args.out)
    if args.append and out_path.exists():
        # Accumulate schemes in one file: a set is identified by
        # (scope, filter, gain, reference catalogue), so re-running the same combination
        # replaces it while a different catalogue or scope is added alongside. This is
        # what lets measure.py offer a choice of precalibration at run time.
        old = json.loads(out_path.read_text())
        merged = []
        new_keys = {(s["scope"], s["filter"], s["gain"],
                     s["reference_catalog"]["name"]) for s in sets}
        for s in old.get("sets", []):
            key = (s["scope"], s["filter"], s["gain"],
                   (s.get("reference_catalog") or {}).get("name", "tycho2"))
            if key in new_keys:
                print(f"  replacing existing set {key}", file=sys.stderr)
                continue
            merged.append(s)
        doc["sets"] = merged + sets
        print(f"  merged: {len(merged)} kept + {len(sets)} new = {len(doc['sets'])} sets",
              file=sys.stderr)

    out_path.write_text(json.dumps(doc, indent=2))
    print(f"wrote {args.out} ({time.time()-t0:.0f}s total)", file=sys.stderr)


if __name__ == "__main__":
    main()
