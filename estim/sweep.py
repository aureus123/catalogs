#!/usr/bin/env python3
"""Extraction hyperparameter sweep for the green-channel photometry pipeline.

ASTAP detection is run once per field per channel, at the lowest snr_min in
the grid; every (aperture_hfd, annulus_in_hfd, snr_min) combination is then
re-measured with photutils over that same detection list (photometry.photometer),
so the grid costs one ASTAP pass per field, not one per combination.

The objective only concerns the green channel / V magnitude pipeline (not the
colour transform, which is fit globally in precalibrate.py from whichever
aperture wins here):

  1. primary: RMS of (V_epoch_a - V_epoch_b) for the same Tycho-2 star seen in
     two independent visits to the same target (ROct x6, SOct x3, YCen x2,
     SVir x2, V462Lup1/2) -- ungameable by catalogue selection, independent
     of Tycho-2's own errors.
  2. secondary: sigma_V(V) = sqrt(a^2 + b^2 * 10^(0.8*(V-10))) fit over the
     pooled per-field zero-point residuals, evaluated at V=10.
  3. hard constraint: completeness >= 90% of a fixed reference calibrator set
     (defined once, at the loosest grid point).

Usage:
    python3 sweep.py --coarse --fields 10 -o sweep_coarse.json
    python3 sweep.py --confirm --coarse-result sweep_coarse.json --top 3 -o sweep_confirm.json
"""
import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

import photometry
from precalibrate import (CACHE_DIR, CALIB_EVT_MAX, CALIB_V_MAX, CALIB_V_MIN,
                           OUR_ISOLATION_ARCSEC, TYCHO2_MATCH_ARCSEC, discover_fields)
from starcat import e_V
from tycho2 import DEFAULT_NPY as TYCHO2_NPY
from tycho2 import Tycho2

# NOTE: this module is still Tycho-2-only. It keys stars by (TYC1, TYC2, TYC3) to pair
# them across epochs, and only Tycho-2 carries those identifiers -- the Gaia backends
# would need a coordinate-based key instead. Not a blocker: the sweep tunes *extraction*
# parameters, which is independent of the reference catalogue.

HERE = Path(__file__).resolve().parent

DEFAULT_GRID = {
    "aperture_hfd": [1.2, 1.6, 2.0, 2.4, 3.0, 3.5],
    "annulus_in_hfd": [3.0, 4.0, 6.0],
    "snr_min": [3, 5, 10, 20],
}
COMPLETENESS_MIN = 0.90
BASELINE = {"aperture_hfd": 2.0, "annulus_in_hfd": 4.0, "snr_min": 3}

REPEAT_GROUPS = {
    "ROct": ["ROct", "ROct2", "ROct3", "ROct4", "ROct5", "ROct6"],
    "SOct": ["SOct", "SOct2", "SOct3"],
    "YCen": ["YCen", "YCen2"],
    "SVir": ["SVir", "SVir2"],
    "V462Lup": ["V462Lup1", "V462Lup2"],
}


def _sigma_clip(x, sigma=3.0, iters=2):
    """A few Tycho-2 stars are themselves variable, and cross-matches
    occasionally land on a blend; clip those before computing scatter so a
    handful of outliers don't dominate the objective."""
    x = np.asarray(x, dtype=float)
    for _ in range(iters):
        if len(x) < 5:
            break
        m, s = np.mean(x), np.std(x)
        keep = np.abs(x - m) <= sigma * s
        if keep.all():
            break
        x = x[keep]
    return x


def pick_representative_fields(all_fields, n):
    """Prefer repeat-group members (both scopes, so the same-target pairs are
    available for the primary objective) then fill up to n with the rest."""
    by_name = {name: (scope, path) for scope, name, path in all_fields}
    chosen = []
    for members in REPEAT_GROUPS.values():
        for m in members:
            if m in by_name and m not in chosen:
                chosen.append(m)
    rest = [name for _, name, _ in all_fields if name not in chosen]
    chosen.extend(rest)
    chosen = chosen[:n] if n else chosen
    return [(by_name[name][0], name, by_name[name][1]) for name in chosen if name in by_name]


def detect_all(fields, cache_root, snr_floor, verbose):
    detections = {}
    image_sizes = {}
    for scope, name, path in fields:
        try:
            field_cache = cache_root / path.parent.name
            field_cache.mkdir(parents=True, exist_ok=True)
            fd = photometry.detect_field(path, workdir=field_cache, snr_floor=snr_floor,
                                          verbose=verbose)
        except Exception as e:
            print(f"  {name}: detect FAILED ({e})", file=sys.stderr)
            continue
        detections[name] = fd
        h = fits.getheader(path)
        image_sizes[name] = (h["NAXIS1"], h["NAXIS2"])
        if verbose:
            print(f"  {name}: detected G={len(fd.det_g)} R={len(fd.det_r)} B={len(fd.det_b)}",
                  file=sys.stderr)
    return detections, image_sizes


def reference_set(detections, cat):
    """Fixed evaluation list: green detections at the loosest grid point that
    match a Tycho-2 calibrator. Completeness is measured against this list,
    never against the candidate's own detections (that would reward keeping
    only bright isolated stars)."""
    ref = {}
    for name, fd in detections.items():
        idx, sep, nn2 = cat.match(*fd.wcs.all_pix2world(fd.det_g["x"], fd.det_g["y"], 1),
                                   radius_arcsec=TYCHO2_MATCH_ARCSEC)
        matched = idx >= 0
        mi = np.nonzero(matched)[0]
        d = cat.data[idx[mi]]
        ok = (d["V"] >= CALIB_V_MIN) & (d["V"] <= CALIB_V_MAX) & (e_V(d) < CALIB_EVT_MAX) & \
             (nn2[mi] > OUR_ISOLATION_ARCSEC)
        keys = set(zip(d["tyc1"][ok].tolist(), d["tyc2"][ok].tolist(), d["tyc3"][ok].tolist()))
        ref[name] = keys
    return ref


def measure_zp_only(fd, image_size, cat, params):
    """Green-only per-field zero point: g + zp ~ V_cat for Tycho-2 calibrators.
    Returns dict tyc_key -> (vg_corrected, V_cat), and the calibrator count."""
    fp = photometry.photometer(
        fd, snr_min=params["snr_min"], max_stars=500,
        aperture_hfd=params["aperture_hfd"], annulus_in_hfd=params["annulus_in_hfd"],
        annulus_out_hfd=params["annulus_in_hfd"] + 2.0,
        backend="photutils", saturation_adu=64000, channel_radius_px=1.5,
    )
    idx, sep, nn2 = cat.match(fp.ra, fp.dec, radius_arcsec=TYCHO2_MATCH_ARCSEC,
                               isolation_arcsec=OUR_ISOLATION_ARCSEC)
    matched = idx >= 0
    mi = np.nonzero(matched)[0]
    d = cat.data[idx[mi]]
    with np.errstate(invalid="ignore"):
        g = -2.5 * np.log10(fp.flux_g[mi])
    ok = (
        (d["V"] >= CALIB_V_MIN) & (d["V"] <= CALIB_V_MAX) & (e_V(d) < CALIB_EVT_MAX) &
        (nn2[mi] > OUR_ISOLATION_ARCSEC) & (fp.isolation_arcsec[mi] > OUR_ISOLATION_ARCSEC) &
        (~fp.saturated[mi]) & np.isfinite(g)
    )
    if ok.sum() < 5:
        return {}, 0
    g_ok = g[ok]
    V_ok = d["V"][ok]
    zp = float(np.mean(V_ok - g_ok))
    out = {}
    for k, rec in enumerate(d[ok]):
        key = (int(rec["tyc1"]), int(rec["tyc2"]), int(rec["tyc3"]))
        out[key] = (g_ok[k] + zp, float(rec["V"]))
    return out, int(ok.sum())


def evaluate_params(detections, image_sizes, cat, params, ref_sets, groups_by_scope):
    per_field = {}
    n_ref_total = 0
    n_ref_hit = 0
    for name, fd in detections.items():
        vgmap, n_calib = measure_zp_only(fd, image_sizes[name], cat, params)
        per_field[name] = vgmap
        ref_keys = ref_sets.get(name, set())
        n_ref_total += len(ref_keys)
        n_ref_hit += len(ref_keys & set(vgmap.keys()))
    completeness = (n_ref_hit / n_ref_total) if n_ref_total else 0.0

    # primary: repeat-pair scatter
    diffs = []
    for group_name, members in groups_by_scope.items():
        present = [m for m in members if m in per_field]
        for a, b in itertools.combinations(present, 2):
            common = set(per_field[a]) & set(per_field[b])
            for key in common:
                va, _ = per_field[a][key]
                vb, _ = per_field[b][key]
                diffs.append(va - vb)
    diffs = _sigma_clip(np.array(diffs))
    repeat_rms = float(np.sqrt(np.mean(diffs ** 2))) if len(diffs) >= 5 else None

    # secondary: pooled sigma_V(V) at V=10
    all_v = []
    all_resid = []
    for name, vgmap in per_field.items():
        for key, (vg, vcat) in vgmap.items():
            all_v.append(vcat)
            all_resid.append(vcat - vg)
    all_v = np.array(all_v)
    all_resid = np.array(all_resid)
    sigma_v10 = None
    if len(all_v) >= 30:
        edges = np.arange(7, 12.01, 1.0)
        centers, rms = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (all_v >= lo) & (all_v < hi)
            if m.sum() >= 5:
                centers.append((lo + hi) / 2.0)
                rms.append(np.sqrt(np.mean(_sigma_clip(all_resid[m]) ** 2)))
        if len(centers) >= 2:
            centers, rms = np.array(centers), np.array(rms)
            weight = 10 ** (0.8 * (centers - 10))
            # sigma^2 = a^2 + b^2 * weight  -> linear least squares in (a^2, b^2)
            A = np.column_stack([np.ones_like(weight), weight])
            coeffs, *_ = np.linalg.lstsq(A, rms ** 2, rcond=None)
            a2, b2 = max(coeffs[0], 0.0), max(coeffs[1], 0.0)
            sigma_v10 = float(np.sqrt(a2 + b2))

    score = repeat_rms if repeat_rms is not None else sigma_v10
    return {
        "params": params, "repeat_rms": repeat_rms, "sigma_v10": sigma_v10,
        "completeness": completeness, "n_repeat_pairs": len(diffs),
        "n_pooled": len(all_v), "score": score,
    }


def run_sweep(fields, grid, cache_root, verbose):
    print(f"loading Tycho-2 ...", file=sys.stderr)
    cat = Tycho2(TYCHO2_NPY)

    snr_floor = min(grid["snr_min"])
    print(f"detecting {len(fields)} fields at snr_floor={snr_floor} ...", file=sys.stderr)
    detections, image_sizes = detect_all(fields, cache_root, snr_floor, verbose)

    groups_by_scope = {}
    for gname, members in REPEAT_GROUPS.items():
        present = [m for m in members if m in detections]
        if len(present) >= 2:
            groups_by_scope[gname] = present
    print(f"repeat groups available: { {k: v for k, v in groups_by_scope.items()} }",
          file=sys.stderr)

    print("building fixed reference calibrator set (baseline detections) ...", file=sys.stderr)
    ref_sets = reference_set(detections, cat)

    combos = [dict(zip(grid.keys(), vals)) for vals in
              itertools.product(*[grid[k] for k in grid])]
    print(f"evaluating {len(combos)} parameter combinations over {len(detections)} fields ...",
          file=sys.stderr)

    results = []
    t0 = time.time()
    for i, params in enumerate(combos):
        r = evaluate_params(detections, image_sizes, cat, params, ref_sets, groups_by_scope)
        results.append(r)
        if verbose or (i + 1) % 20 == 0:
            print(f"  [{i+1}/{len(combos)}] {params} -> repeat_rms={r['repeat_rms']} "
                  f"sigma_v10={r['sigma_v10']} completeness={r['completeness']:.2f} "
                  f"({time.time()-t0:.0f}s)", file=sys.stderr)

    ok = [r for r in results if r["completeness"] >= COMPLETENESS_MIN and r["score"] is not None]
    ranked = sorted(ok, key=lambda r: r["score"])
    if not ranked:
        print("warning: no combination met the completeness constraint; "
              "ranking all combinations instead", file=sys.stderr)
        ranked = sorted([r for r in results if r["score"] is not None], key=lambda r: r["score"])
    return ranked, results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--coarse", action="store_true")
    mode.add_argument("--confirm", action="store_true")
    ap.add_argument("--fields", type=int, default=10, help="(coarse) number of representative fields")
    ap.add_argument("--all-vis", action="store_true", help="(confirm) use every VIS/gain field")
    ap.add_argument("--coarse-result", default=str(HERE / "sweep_coarse.json"))
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--filter", default="VIS")
    ap.add_argument("--gain", type=int, default=60)
    ap.add_argument("--scope", default=None, help="restrict to one TELESCOP value")
    ap.add_argument("--cache-dir", default=str(CACHE_DIR))
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    all_fields = list(discover_fields(args.filter, args.gain))
    if args.scope:
        all_fields = [f for f in all_fields if f[0] == args.scope]
    cache_root = Path(args.cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)

    if args.coarse:
        fields = pick_representative_fields(all_fields, args.fields)
        print(f"coarse sweep over {len(fields)} fields: {[f[1] for f in fields]}",
              file=sys.stderr)
        ranked, all_results = run_sweep(fields, DEFAULT_GRID, cache_root, args.verbose)
        doc = {"mode": "coarse", "fields": [f[1] for f in fields],
               "grid": DEFAULT_GRID, "ranked": ranked}
        Path(args.out).write_text(json.dumps(doc, indent=2, default=str))
        print(f"wrote {args.out}: top candidate {ranked[0] if ranked else None}", file=sys.stderr)

    else:  # confirm
        coarse = json.loads(Path(args.coarse_result).read_text())
        candidates = [r["params"] for r in coarse["ranked"][: args.top]]
        fields = all_fields if args.all_vis else pick_representative_fields(all_fields, args.fields)
        print(f"confirming top {len(candidates)} candidates over {len(fields)} fields ...",
              file=sys.stderr)
        grid = {k: sorted({c[k] for c in candidates}) for k in candidates[0]}
        ranked, all_results = run_sweep(fields, grid, cache_root, args.verbose)
        astap_baseline = next((r for r in all_results if r["params"] == BASELINE), None)
        doc = {"mode": "confirm", "fields": [f[1] for f in fields], "candidates": candidates,
               "ranked": ranked, "astap_baseline": astap_baseline}
        Path(args.out).write_text(json.dumps(doc, indent=2, default=str))
        print(f"wrote {args.out}: winner {ranked[0] if ranked else None}", file=sys.stderr)


if __name__ == "__main__":
    main()
