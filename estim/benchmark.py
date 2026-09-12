#!/usr/bin/env python3
"""Validation suite for the precalibrated coefficients.

    python3 benchmark.py --loo-star --loo-field --repeats --ablation

Reuses precalibrate.py's field cache (estim/cache/) built by the last
`precalibrate.py` run, so nothing here re-runs ASTAP detection from scratch
except the repeat-field and ablation checks, which need a couple of
photometry passes at specific parameter choices.

Not automated here (see estim/README.md "Note on the extracted channel" and
the design plan for the numbers that were obtained manually): the raw-pipeline
vs onboard-cube comparison, which needs a full ASTAP dark/flat/bias stack per
field. --cube-inventory only reports which fields still have a raw/ folder.
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

import photometry
import precalibrate as pc
import sweep
from tycho2 import DEFAULT_NPY as TYCHO2_NPY
from tycho2 import Tycho2

HERE = Path(__file__).resolve().parent
#: Image archive root -- see precalibrate.ARCHIVE_ROOT for the expected layout.
ARCHIVE_ROOT = pc.ARCHIVE_ROOT


def build_pooled(scope, fields, cat, extraction, cache_root, verbose=False):
    """Same as precalibrate.fit_scope's data collection stage, exposed here
    for the benchmark's own analyses (LOO, ablations)."""
    tables = []
    used_fields = []
    field_headers = {}
    for i, (name, cube_path) in enumerate(fields):
        try:
            fp = pc.measure_one_field(cube_path, cache_root, extraction)
        except Exception as e:
            if verbose:
                print(f"  [{scope}] {name}: FAILED ({e})", file=sys.stderr)
            continue
        h = fits.getheader(cube_path)
        tab = pc.collect_calibrators(fp, cat, i, (h["NAXIS1"], h["NAXIS2"]))
        if tab is None or len(tab["g"]) < 5:
            continue
        tables.append(tab)
        used_fields.append(name)
        field_headers[i] = name
    if not tables:
        return None, None
    pooled = {k: np.concatenate([t[k] for t in tables]) for k in tables[0]}
    return pooled, field_headers


def fit_pooled(pooled, colour_index):
    field_id = pooled["field"]
    bv_c = pc.within_center(pooled["BV_cat"], field_id)
    c = pc.instrumental_colour(colour_index, pooled["g"], pooled["r"], pooled["b"])
    c_c = pc.within_center(c, field_id)
    slope, se_slope, resid_std, resid_c = pc.fit_slope_through_origin(bv_c, c_c)
    Tbv = 1.0 / slope

    y = pooled["V_cat"] - pooled["g"]
    y_c = pc.within_center(y, field_id)
    r2_c = pc.within_center(pooled["r2"], field_id)
    Tv_bv, k_r2, cov, resid_v = pc.fit_two_covariate(bv_c, r2_c, y_c)

    return dict(field_id=field_id, bv_c=bv_c, c_c=c_c, slope=slope, resid_c=resid_c, Tbv=Tbv,
                r2_c=r2_c, y_c=y_c, Tv_bv=Tv_bv, k_r2=k_r2, resid_v=resid_v)


def leave_one_out_residuals(x, resid):
    """Exact LOO residuals for y = slope*x (through the origin, already
    group-centered) via the leverage identity resid_loo = resid/(1-h)."""
    h = x ** 2 / np.sum(x ** 2)
    return resid / (1.0 - h)


def leave_one_out_residuals_2cov(X, resid):
    XtX_inv = np.linalg.inv(X.T @ X)
    h = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    return resid / (1.0 - h)


def rms(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x ** 2))) if len(x) else float("nan")


def do_loo_star(scope, pooled, verbose):
    if pooled is None:
        print(f"[{scope}] loo-star: no usable data", file=sys.stderr)
        return
    fit = fit_pooled(pooled, "b-r")
    loo_c = leave_one_out_residuals(fit["bv_c"], fit["resid_c"])
    loo_bv = loo_c * fit["Tbv"]

    X = np.column_stack([fit["bv_c"], fit["r2_c"]])
    loo_v = leave_one_out_residuals_2cov(X, fit["resid_v"])

    print(f"\n[{scope}] leave-one-star-out (n={len(loo_v)}, analytic via regression leverage,\n"
          f"  approximation: ignores the sub-leverage from each field's own group mean):")
    for lo, hi in [(7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]:
        m = (pooled["V_cat"] >= lo) & (pooled["V_cat"] < hi)
        if m.sum() >= 5:
            print(f"    V {lo}-{hi}: n={m.sum():5d}  RMS(V)={rms(loo_v[m]):.3f}  RMS(B-V)={rms(loo_bv[m]):.3f}")
    print(f"    overall: RMS(V)={rms(loo_v):.3f}  RMS(B-V)={rms(loo_bv):.3f}")


def do_loo_field(scope, pooled, n_fields, verbose):
    if pooled is None or n_fields < 4:
        print(f"[{scope}] loo-field: not enough fields", file=sys.stderr)
        return
    field_id = pooled["field"]
    all_ids = np.unique(field_id)
    held_out_v = []
    held_out_bv = []
    for fid in all_ids:
        train = field_id != fid
        test = field_id == fid
        if test.sum() < 5 or train.sum() < 100:
            continue
        sub = {k: pooled[k][train] for k in pooled}
        fit = fit_pooled(sub, "b-r")
        Tbv_loo, Tv_bv_loo, k_loo = fit["Tbv"], fit["Tv_bv"], fit["k_r2"]

        g, r_, b_ = pooled["g"][test], pooled["r"][test], pooled["b"][test]
        BV_cat, V_cat, r2 = pooled["BV_cat"][test], pooled["V_cat"][test], pooled["r2"][test]
        c = pc.instrumental_colour("b-r", g, r_, b_)
        ZPbv = np.mean(BV_cat - Tbv_loo * c)
        zp = np.mean(V_cat - g - Tv_bv_loo * BV_cat - k_loo * r2)
        BVest = ZPbv + Tbv_loo * c
        Vest = g + zp + Tv_bv_loo * BVest + k_loo * r2
        held_out_v.append(V_cat - Vest)
        held_out_bv.append(BV_cat - BVest)

    held_out_v = np.concatenate(held_out_v)
    held_out_bv = np.concatenate(held_out_bv)
    print(f"\n[{scope}] leave-one-field-out (n_fields={len(all_ids)}, n_stars={len(held_out_v)}):")
    print(f"    RMS(V)={rms(held_out_v):.3f}  RMS(B-V)={rms(held_out_bv):.3f}")
    print(f"    (compare to the pooled in-sample sigma_V_by_band in instrument_coeffs.json --"
          f" close agreement means Tbv/Tv_bv/k really do transfer between fields)")


def do_repeats(coeffs, cache_root, verbose):
    cat = Tycho2(TYCHO2_NPY)
    for cset in coeffs["sets"]:
        scope = cset["scope"]
        extraction = {
            "snr_min": cset["extraction"]["detection"]["snr_min"],
            "max_stars": cset["extraction"]["detection"]["max_stars"],
            "aperture_hfd": cset["extraction"]["photometry"]["aperture_hfd"],
            "annulus_in_hfd": cset["extraction"]["photometry"]["annulus_in_hfd"],
        }
        subdir = "Dwarf2" if scope == "DWARF II" else "Dwarf3"
        fields = []
        for group in sweep.REPEAT_GROUPS.values():
            for name in group:
                folder = ARCHIVE_ROOT / subdir / name
                path = pc._find_cube_fits(folder) if folder.is_dir() else None
                if path is None:
                    continue
                try:
                    h = fits.getheader(path)
                except Exception:
                    continue
                if h.get("FILTER") != cset["filter"] or h.get("GAIN") != cset["gain"]:
                    if verbose:
                        print(f"  {name}: filter/gain mismatch "
                              f"({h.get('FILTER')}/{h.get('GAIN')} != "
                              f"{cset['filter']}/{cset['gain']}), excluded from repeats",
                              file=sys.stderr)
                    continue
                fields.append((scope, name, path))
        if not fields:
            continue
        detections, image_sizes = sweep.detect_all(fields, cache_root, extraction["snr_min"], verbose)
        groups = {g: [m for m in members if m in detections]
                  for g, members in sweep.REPEAT_GROUPS.items()}
        groups = {g: m for g, m in groups.items() if len(m) >= 2}
        ref = sweep.reference_set(detections, cat)
        result = sweep.evaluate_params(detections, image_sizes, cat,
                                        {"aperture_hfd": extraction["aperture_hfd"],
                                         "annulus_in_hfd": extraction["annulus_in_hfd"],
                                         "snr_min": extraction["snr_min"]},
                                        ref, groups)
        print(f"\n[{scope}] repeat-field scatter at the precalibrated extraction settings:")
        print(f"    groups: { {g: v for g, v in groups.items()} }")
        print(f"    n_repeat_pairs={result['n_repeat_pairs']}  RMS(V epoch-to-epoch)={result['repeat_rms']}")


def do_ablation_colour_r2(scope, pooled):
    if pooled is None:
        return
    field_id = pooled["field"]
    y = pooled["V_cat"] - pooled["g"]
    y_c = pc.within_center(y, field_id)
    bv_c = pc.within_center(pooled["BV_cat"], field_id)
    r2_c = pc.within_center(pooled["r2"], field_id)

    # full model: Tv_bv * BV + k * r2
    Tv_bv, k_r2, _, resid_full = pc.fit_two_covariate(bv_c, r2_c, y_c)
    # -colour: r2 only
    k_only = np.sum(r2_c * y_c) / np.sum(r2_c ** 2)
    resid_no_colour_term = y_c - k_only * r2_c
    # -r2: BV only (still "with colour", just no radial term)
    tv_only, _, _, resid_no_r2 = pc.fit_slope_through_origin(bv_c, y_c)
    # green-only, no correction at all (the naive comparison)
    resid_naive = y_c

    print(f"\n[{scope}] ablation on the magnitude-transform covariates (n={len(y_c)}):")
    print(f"    full model (Tv_bv*BV + k*r2):        RMS={rms(resid_full):.4f}")
    print(f"    -colour (k*r2 only, no Tv_bv term):   RMS={rms(resid_no_colour_term):.4f}")
    print(f"    -r2 (Tv_bv*BV only):                  RMS={rms(resid_no_r2):.4f}")
    print(f"    naive (zero-point only, no terms):    RMS={rms(resid_naive):.4f}")


def do_ablation_defocus(cache_root, verbose):
    cat = Tycho2(TYCHO2_NPY)
    names = ["V462Lup1", "V462Lup2", "V462LupDefocus"]
    fields = []
    for n in names:
        folder = ARCHIVE_ROOT / "Dwarf3" / n
        p = pc._find_cube_fits(folder) if folder.is_dir() else None
        if p is not None:
            fields.append(("DWARF 3", n, p))
    if len(fields) < 2:
        print("\nV462Lup focused/defocused fields not found, skipping", file=sys.stderr)
        return
    detections, image_sizes = sweep.detect_all(fields, cache_root, 5, verbose)
    present = [n for n in names if n in detections]
    print(f"\nfocus ablation, fields present: {present}")
    for a, b in itertools.combinations(present, 2):
        vg_a, _ = sweep.measure_zp_only(detections[a], image_sizes[a], cat, sweep.BASELINE)
        vg_b, _ = sweep.measure_zp_only(detections[b], image_sizes[b], cat, sweep.BASELINE)
        common = set(vg_a) & set(vg_b)
        if len(common) < 3:
            print(f"    {a} vs {b}: only {len(common)} common calibrators, skipped")
            continue
        diffs = np.array([vg_a[k][0] - vg_b[k][0] for k in common])
        print(f"    {a} vs {b}: n={len(diffs)}  RMS={rms(diffs):.3f}")


def do_cube_inventory():
    n_with_raw, n_without = 0, 0
    for root in (ARCHIVE_ROOT / d for d in pc.SCOPE_DIRS):
        if not root.is_dir():
            continue
        for folder in sorted(root.iterdir()):
            if not folder.is_dir():
                continue
            if (folder / "raw").is_dir():
                n_with_raw += 1
            elif (folder / f"{folder.name}.fits").exists():
                n_without += 1
    print(f"\ncube-vs-raw inventory: {n_with_raw} fields still have a raw/ subfolder, "
          f"{n_without} cube-only.")
    print("  Full raw-pipeline re-measurement (ASTAP dark/flat/bias stack per field) is not "
          "automated here -- see the design plan for the 13-field manual comparison already done "
          "(median RMS 0.117 cube vs 0.125 raw-recalibrated).")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coeffs", default=str(HERE / "instrument_coeffs.json"))
    ap.add_argument("--cache-dir", default=str(HERE / "cache"))
    ap.add_argument("--loo-star", action="store_true")
    ap.add_argument("--loo-field", action="store_true")
    ap.add_argument("--repeats", action="store_true")
    ap.add_argument("--ablation", action="store_true")
    ap.add_argument("--cube-inventory", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    if args.all:
        args.loo_star = args.loo_field = args.repeats = args.ablation = args.cube_inventory = True

    coeffs = json.loads(Path(args.coeffs).read_text())
    cache_root = Path(args.cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)
    cat = Tycho2(TYCHO2_NPY)

    by_scope = {}
    for scope, name, path in pc.discover_fields("VIS", 60):
        by_scope.setdefault(scope, []).append((name, path))

    if args.loo_star or args.loo_field or args.ablation:
        for cset in coeffs["sets"]:
            scope = cset["scope"]
            fields = by_scope.get(scope, [])
            extraction = {
                "filter": cset["filter"], "gain": cset["gain"],
                "snr_min": cset["extraction"]["detection"]["snr_min"],
                "max_stars": cset["extraction"]["detection"]["max_stars"],
                "aperture_hfd": cset["extraction"]["photometry"]["aperture_hfd"],
                "annulus_in_hfd": cset["extraction"]["photometry"]["annulus_in_hfd"],
                "annulus_out_hfd": cset["extraction"]["photometry"]["annulus_out_hfd"],
                "backend": cset["extraction"]["photometry"]["backend"],
                "saturation_adu": cset["extraction"]["saturation_adu"],
                "channel_radius_px": cset["extraction"]["match"]["channel_radius_px"],
            }
            pooled, field_names = build_pooled(scope, fields, cat, extraction, cache_root, args.verbose)
            n_fields = len(field_names) if field_names else 0
            if args.loo_star:
                do_loo_star(scope, pooled, args.verbose)
            if args.loo_field:
                do_loo_field(scope, pooled, n_fields, args.verbose)
            if args.ablation:
                do_ablation_colour_r2(scope, pooled)

    if args.repeats:
        do_repeats(coeffs, cache_root, args.verbose)
    if args.ablation:
        do_ablation_defocus(cache_root, args.verbose)
    if args.cube_inventory:
        do_cube_inventory()


if __name__ == "__main__":
    main()
