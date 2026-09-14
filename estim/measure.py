#!/usr/bin/env python3
"""Estimate V and B-V for every star detected in one Dwarf 2/3 FITS cube.

    python3 measure.py RCen.fits -c instrument_coeffs.json -o RCen.csv

Looks up the (scope, filter, gain) triplet from the FITS header in
instrument_coeffs.json, takes the extraction settings and the colour/
magnitude transform coefficients (Tbv, Tv_bv, k_r2) from there, fits only
the two zero points (zp, ZPbv) that do not transfer between images from
this image's own Tycho-2-matched stars, and applies the model:

    g = -2.5 log10(flux_G)         c = <colour_index> from the JSON
    (B-V)est = ZPbv + Tbv * c
    Vest     = g + zp + Tv_bv * (B-V)est + k_r2 * r2

to every detected star, matched or not.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

import photometry
import starcat
from starcat import e_V, get_catalog, isolated
from tycho2 import DEFAULT_NPY as TYCHO2_NPY

CALIB_V_MIN, CALIB_V_MAX = 7.0, 13.0
CALIB_EVT_MAX = 0.10
OUR_ISOLATION_ARCSEC = 20.0
SNR_MIN_CALIB = 5.0
LOWSNR_THRESHOLD = 10.0
EDGE_MARGIN_PX = 20.0


def instrumental_colour(name, g, r, b):
    return {"b-r": b - r, "b-g": b - g, "g-r": g - r}[name]


def _band_floor(bands, V):
    """Per-star empirical error floor, interpolated from a ``sigma_*_by_band`` table.

    The table is keyed by magnitude band ("8-9": 0.058, ...). Stars outside the tabulated
    range take the nearest band's value rather than zero -- a star fainter than the
    faintest calibrated band is *less* well measured, not perfectly measured, so falling
    back to no floor at all would be exactly backwards.
    """
    V = np.asarray(V, dtype=float)
    if not bands:
        return np.zeros(V.shape)
    centres, values = [], []
    for band, s in sorted(bands.items()):
        try:
            lo, hi = (float(x) for x in band.split("-"))
        except ValueError:
            continue
        centres.append((lo + hi) / 2.0)
        values.append(float(s))
    if not centres:
        return np.zeros(V.shape)
    centres = np.asarray(centres); values = np.asarray(values)
    order = np.argsort(centres)
    out = np.interp(np.nan_to_num(V, nan=float(centres[order][0])),
                    centres[order], values[order])
    return np.where(np.isfinite(V), out, values[order][-1])


def set_catalog(s):
    """Which reference a set was precalibrated against. Schema >= 4 records it per set;
    older files had it top level or not at all (they were all Tycho-2)."""
    return (s.get("reference_catalog") or {}).get("name", "tycho2")


def find_set(coeffs, scope, filt, gain, catalog):
    """Locate the coefficient set to use, given the chosen *fit* catalogue.

    The precalibration is derived from the fit catalogue rather than chosen separately,
    so the two can never land on different photometric systems. Preference order comes
    from ``starcat.PRECAL_PREFERENCE``: Tycho-2 uses only Tycho-2, while either Gaia
    backend prefers the unquantised online precalibration and falls back to V50.

    Returns ``(set, precal_name, available)``.
    """
    scope = photometry.normalize_scope(scope)
    cands = [s for s in coeffs["sets"]
             if photometry.normalize_scope(s["scope"]) == scope
             and s["filter"] == filt and s["gain"] == gain]
    available = sorted({set_catalog(s) for s in cands})
    for want in starcat.PRECAL_PREFERENCE[catalog]:
        for s in cands:
            if set_catalog(s) == want:
                return s, want, available
    return None, None, available


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cube", help="3-channel FITS cube from Dwarf 2/3")
    ap.add_argument("-c", "--coeffs", required=True, help="instrument_coeffs.json")
    ap.add_argument("-o", "--out", required=True, help="output CSV path")
    ap.add_argument("--tycho2", default=str(TYCHO2_NPY))
    ap.add_argument("--catalog", choices=starcat.CATALOG_NAMES, default="gaia_v50",
                     help="fit method: the catalogue used to match this image and fit "
                          "its zero points. The matching precalibration is selected "
                          "automatically -- 'tycho2' uses the Tycho-2 coefficients, "
                          "while 'gaia_v50' and 'gaia_online' both prefer the online "
                          "precalibration and fall back to V50. 'gaia_online' REQUIRES "
                          "AN INTERNET CONNECTION; 'gaia_v50' is the offline default.")
    ap.add_argument("--wcs", default=None,
                     help="use this .wcs instead of the one beside the cube. Needed when "
                          "an archived .wcs is stale -- a wrong plate scale silently "
                          "shifts every position and collapses the catalogue match rate.")
    ap.add_argument("--workdir", help="directory for extracted planes (default: next to cube)")
    ap.add_argument("--force", action="store_true",
                     help="proceed without a matching (scope,filter,gain) set; every row flagged UNCALIBRATED")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    cube_path = Path(args.cube)
    header = fits.getheader(cube_path)
    scope = photometry.normalize_scope(header.get("TELESCOP"))
    gain = header.get("GAIN", -1)
    # Pre-rename firmware writes no FILTER card, so fall back to shotsInfo.json and
    # canonicalise CUT -> VIS / PASS -> Astro.
    filt, filt_raw, filt_src = photometry.resolve_filter(header, cube_path)
    if filt_src == "shotsInfo.json":
        print(f"note: no FILTER in the header; taking filter {filt_raw!r} from "
              f"shotsInfo.json", file=sys.stderr)
    if filt_raw and filt_raw.upper() in photometry.LEGACY_FILTER_NAMES:
        print(
            f"\n*** WARNING: this frame's filter is {filt_raw!r}, the pre-2025-11 firmware\n"
            f"*** name for {filt!r}. Same physical filter, but the older firmware applies\n"
            f"*** different on-board colour processing: it delivers ~22% less instrumental\n"
            f"*** colour separation, so Tbv is 2.21 rather than 1.34 (16 sigma apart).\n"
            f"*** The {filt!r} coefficients being applied here were fitted on post-rename\n"
            f"*** frames, where they give a V residual of 0.068 mag; on {filt_raw!r} frames\n"
            f"*** they give 0.312 mag -- about 4.6x worse.\n"
            f"*** Treat V as approximate and B-V as unreliable. Precalibrating the\n"
            f"*** pre-rename era separately would fix this properly.\n",
            file=sys.stderr)
    nx, ny = header["NAXIS1"], header["NAXIS2"]

    coeffs = json.loads(Path(args.coeffs).read_text())
    cset, precal_name, available = find_set(coeffs, scope, filt, gain, args.catalog)
    if cset is None:
        msg = (f"no coefficient set for (scope={scope!r}, filter={filt!r}, gain={gain!r}) "
               f"usable with --catalog {args.catalog}: wanted a precalibration in "
               f"{list(starcat.PRECAL_PREFERENCE[args.catalog])}"
               + (f", have {available}" if available else ", and this scope has none"))
        if not args.force:
            print(f"error: {msg}. Use --force to proceed uncalibrated.", file=sys.stderr)
            sys.exit(1)
        print(f"warning: {msg}, proceeding uncalibrated (--force)", file=sys.stderr)
    elif precal_name != starcat.PRECAL_PREFERENCE[args.catalog][0]:
        print(f"note: no {starcat.PRECAL_PREFERENCE[args.catalog][0]!r} precalibration "
              f"for this scope, falling back to {precal_name!r}", file=sys.stderr)

    if cset is not None:
        det = cset["extraction"]["detection"]
        pho = cset["extraction"]["photometry"]
        match = cset["extraction"]["match"]
        colour_index = cset["colour_index"]
        Tbv = cset["coefficients"]["Tbv"]["value"]
        sigma_Tbv = cset["coefficients"]["Tbv"]["sigma"]
        Tv_bv = cset["coefficients"]["Tv_bv"]["value"]
        sigma_Tv_bv = cset["coefficients"]["Tv_bv"]["sigma"]
        k_r2 = cset["coefficients"]["k_r2"]["value"]
        valid_bv_range = cset["valid_bv_range"]
        saturation_adu = cset["extraction"]["saturation_adu"]
    else:
        det = {"snr_min": 5, "max_stars": 500}
        pho = {"backend": "photutils", "aperture_hfd": 2.0, "annulus_in_hfd": 4.0, "annulus_out_hfd": 6.0}
        match = {"channel_radius_px": 1.5, "tycho2_radius_arcsec": 8.0, "isolation_arcsec": 20.0}
        colour_index, Tbv, sigma_Tbv, Tv_bv, sigma_Tv_bv, k_r2 = "b-r", 1.5, 0.0, 0.0, 0.0, 0.0
        valid_bv_range = [0.0, 2.0]
        saturation_adu = 64000

    if args.verbose:
        print(f"measuring {cube_path} (scope={scope}, filter={filt}, gain={gain}) ...", file=sys.stderr)
    fp = photometry.measure_field(
        cube_path, workdir=args.workdir,
        snr_min=det["snr_min"], max_stars=det["max_stars"],
        aperture_hfd=pho["aperture_hfd"], annulus_in_hfd=pho["annulus_in_hfd"],
        annulus_out_hfd=pho["annulus_out_hfd"], backend=pho["backend"],
        saturation_adu=saturation_adu, channel_radius_px=match["channel_radius_px"],
        verbose=args.verbose, wcs_path=args.wcs,
    )
    n_extracted = fp.n()

    # The reference catalogue comes from the JSON, not from a default: applying a
    # coefficient set with a different reference than it was fitted against would
    # reintroduce exactly the mismatch the provenance block exists to prevent.
    ref = (cset.get("reference_catalog") if cset else None) \
        or coeffs.get("reference_catalog") or {}
    cat_name = args.catalog
    precal_name = precal_name or ref.get("name", "tycho2")

    # PRECAL_PREFERENCE guarantees precal and fit share a photometric system, so this
    # can only trip if the JSON is hand-edited. Kept as a cheap invariant check: the
    # colour terms genuinely do not transfer between systems (fitted on the same 63
    # Dwarf 3 fields, Tycho-2 gives Tbv 1.671 / Tv_bv -0.039 against Gaia's 1.527 /
    # -0.013), so a mismatch would bias B-V by (dTbv * colour).
    if cset is not None and not starcat.same_system(precal_name, cat_name):
        msg = (f"photometric system mismatch: coefficients precalibrated against "
               f"{precal_name!r} ({starcat.PHOTOMETRIC_SYSTEM.get(precal_name)}) cannot "
               f"be applied while fitting with {cat_name!r} "
               f"({starcat.PHOTOMETRIC_SYSTEM.get(cat_name)}).")
        if not args.force:
            print(f"error: {msg} Use --force to override.", file=sys.stderr)
            sys.exit(1)
        print(f"warning: {msg} Proceeding anyway (--force).", file=sys.stderr)
    if cat_name in starcat.NEEDS_NETWORK:
        print(f"note: catalogue '{cat_name}' needs an internet connection "
              f"(cached fields are reused offline)", file=sys.stderr)
    cat = get_catalog(cat_name, npy_path=args.tycho2,
                       mag_limit=ref.get("mag_limit", 14.0),
                       epoch=ref.get("epoch"), verbose=args.verbose)
    idx, sep, nn2 = cat.match(fp.ra, fp.dec, radius_arcsec=match["tycho2_radius_arcsec"],
                               isolation_arcsec=match["isolation_arcsec"])
    matched = idx >= 0
    n_matched = int(np.sum(matched))
    # A healthy field matches most of its detections. A low rate almost always means the
    # astrometry is wrong -- typically a stale .wcs beside the cube whose plate scale no
    # longer matches the image -- rather than a genuinely sparse field.
    if n_extracted and n_matched / n_extracted < 0.4:
        print(f"warning: only {n_matched}/{n_extracted} detections matched {cat_name}. "
              f"That usually means a bad WCS (check the plate scale in the .wcs beside "
              f"the cube, and re-solve / pass --wcs).", file=sys.stderr)

    with np.errstate(invalid="ignore", divide="ignore"):
        g = -2.5 * np.log10(fp.flux_g)
        r_mag = -2.5 * np.log10(fp.flux_r)
        b_mag = -2.5 * np.log10(fp.flux_b)
    c = instrumental_colour(colour_index, g, r_mag, b_mag)

    cx, cy = nx / 2.0, ny / 2.0
    r2 = ((fp.x - cx) ** 2 + (fp.y - cy) ** 2) / (cx ** 2 + cy ** 2)

    # --- fit this image's zero points from its own Tycho-2-matched calibrators ---
    calib = np.zeros(fp.n(), dtype=bool)
    calib[matched] = True
    d_all = np.full(fp.n(), None, dtype=object)
    tyc_V = np.full(fp.n(), np.nan)
    tyc_BV = np.full(fp.n(), np.nan)
    tyc_name = np.full(fp.n(), "", dtype=object)
    mi = np.nonzero(matched)[0]
    dmatch = cat.data[idx[mi]]
    tyc_V[mi] = dmatch["V"]
    tyc_BV[mi] = dmatch["BV"]
    for j, k in zip(mi, idx[mi]):
        tyc_name[j] = cat.designation(k)

    with np.errstate(invalid="ignore"):
        flux_ok = (fp.flux_r > 0) & (fp.flux_g > 0) & (fp.flux_b > 0) & \
                  np.isfinite(fp.flux_r) & np.isfinite(fp.flux_g) & np.isfinite(fp.flux_b)
    # Flux-ratio blend test rather than "any neighbour within 20 arcsec" -- the latter
    # gets stricter as the catalogue deepens, so a better reference would shrink the
    # calibrator sample. See starcat.blend_dmag.
    blend_cfg = (coeffs.get("calibrator_cuts") or {}).get("blend") or {}
    iso = isolated(cat.data, idx,
                   radius_arcsec=blend_cfg.get("radius_arcsec", match["isolation_arcsec"]),
                   dmag=blend_cfg.get("dmag", starcat.DEFAULT_BLEND_DMAG))
    is_calib = np.zeros(fp.n(), dtype=bool)
    is_calib[mi] = (
        (tyc_V[mi] >= CALIB_V_MIN) & (tyc_V[mi] <= CALIB_V_MAX) &
        (e_V(dmatch) < CALIB_EVT_MAX) & (~np.isnan(tyc_BV[mi])) &
        iso[mi] &
        (fp.isolation_arcsec[mi] > OUR_ISOLATION_ARCSEC) &
        (~fp.saturated[mi]) & (fp.snr_g[mi] >= SNR_MIN_CALIB) & flux_ok[mi]
    )
    n_calib = int(np.sum(is_calib))
    if n_calib < 3:
        print(f"error: only {n_calib} usable calibrators in this field (need >= 3)", file=sys.stderr)
        sys.exit(1)

    zpbv_i = tyc_BV[is_calib] - Tbv * c[is_calib]
    ZPbv = float(np.mean(zpbv_i))
    sigma_ZPbv = float(np.std(zpbv_i, ddof=1) / np.sqrt(n_calib)) if n_calib > 1 else 0.0

    zp_i = tyc_V[is_calib] - g[is_calib] - Tv_bv * tyc_BV[is_calib] - k_r2 * r2[is_calib]
    zp = float(np.mean(zp_i))
    sigma_zp = float(np.std(zp_i, ddof=1) / np.sqrt(n_calib)) if n_calib > 1 else 0.0
    rms_zp = float(np.sqrt(np.mean((zp_i - zp) ** 2)))

    # --- apply the model to every detected star ---
    BVest = ZPbv + Tbv * c
    Vest = g + zp + Tv_bv * BVest + k_r2 * r2

    # Green-only fallback. A star undetected in R or B has no colour index, which
    # previously lost its V entirely -- and that silently discarded exactly the very red
    # targets the pipeline exists to measure (T Pup: bright in green at SNR 852, absent
    # from the blue plane, hence no magnitude at all).
    #
    # The colour term is small enough to survive being guessed: |Tv_bv| is 0.013-0.020,
    # so a whole magnitude of B-V error costs under 0.02 mag. Substituting the field's
    # own mean calibrator colour keeps the estimate unbiased on average -- dropping the
    # term outright would instead assume B-V = 0, biasing every fallback star by
    # Tv_bv * <B-V>, since zp was fitted with the colour term present.
    no_colour = ~np.isfinite(c) & np.isfinite(g)
    n_fallback = int(np.sum(no_colour))
    bv_assumed = float(np.mean(tyc_BV[is_calib]))
    bv_spread = float(np.std(tyc_BV[is_calib])) if n_calib > 1 else 0.5
    sigma_nocolour = np.zeros(fp.n())
    if n_fallback:
        Vest[no_colour] = (g[no_colour] + zp + Tv_bv * bv_assumed + k_r2 * r2[no_colour])
        BVest[no_colour] = np.nan          # V is recoverable, the colour genuinely is not
        sigma_nocolour[no_colour] = abs(Tv_bv) * bv_spread

    sigma_phot_g = 1.0857 / fp.snr_g
    sigma_phot_r = 1.0857 / fp.snr_r
    sigma_phot_b = 1.0857 / fp.snr_b
    BV_err = np.sqrt(Tbv ** 2 * (sigma_phot_b ** 2 + sigma_phot_r ** 2) + sigma_ZPbv ** 2 +
                      (sigma_Tbv * c) ** 2)
    # Saturated stars get an uncertainty covering the measured bias as well as its
    # scatter -- without this, a clipped star reports a photon-noise sigma of ~0.01 mag
    # while sitting up to 0.9 mag faint. See starcat.SATURATION_SIGMA.
    sigma_sat = starcat.saturation_sigma(fp.n_sat_px)
    # BV_err/BVest are NaN for fallback stars, so their colour terms drop out; the
    # substituted-colour uncertainty takes their place via sigma_nocolour.
    colour_var = np.where(no_colour, 0.0, np.nan_to_num((Tv_bv * BV_err) ** 2)
                          + np.nan_to_num((sigma_Tv_bv * BVest) ** 2))
    # Empirical noise floor. Everything above is *modelled* -- photon noise, the
    # zero-point fit, coefficient uncertainty -- and on a bright star those sum to about
    # 0.005 mag, a precision this pipeline has never demonstrated. sigma_V_by_band comes
    # from the leave-one-out residuals measured during precalibration, so it is what the
    # pipeline actually achieves on stars of that brightness, including the systematics
    # no noise model sees: transparency drift, flat-field error, position on the chip.
    # It is consistent with the 0.046 mag night-to-night repeat-field scatter, and with
    # ASTAP's own MERR, which floors on the check star's observed scatter for the same
    # reason (unit_aavso.pas:2511, max(2/SNR, sigma_check)).
    floor_V = _band_floor((cset or {}).get("sigma_V_by_band"), Vest)
    floor_BV = _band_floor((cset or {}).get("sigma_BV_by_band"), Vest)
    V_err = np.sqrt(sigma_phot_g ** 2 + sigma_zp ** 2 + colour_var
                     + sigma_sat ** 2 + sigma_nocolour ** 2 + floor_V ** 2)
    BV_err = np.sqrt(BV_err ** 2 + sigma_sat ** 2 + floor_BV ** 2)

    # --- flags ---
    flags = [[] for _ in range(fp.n())]
    for k in range(fp.n()):
        fl = flags[k]
        if fp.saturated[k]:
            fl.append("SAT")
        if fp.isolation_arcsec[k] < OUR_ISOLATION_ARCSEC:
            fl.append("BLEND")
        if fp.x[k] < EDGE_MARGIN_PX or fp.x[k] > nx - EDGE_MARGIN_PX or \
           fp.y[k] < EDGE_MARGIN_PX or fp.y[k] > ny - EDGE_MARGIN_PX:
            fl.append("EDGE")
        if fp.snr_g[k] < LOWSNR_THRESHOLD:
            fl.append("LOWSNR")
        if not (valid_bv_range[0] <= BVest[k] <= valid_bv_range[1]):
            fl.append("BV_EXTRAP")
        if not flux_ok[k]:
            fl.append("NO_COLOR")      # V from the green-only fallback, B-V unavailable
        if cset is None:
            fl.append("UNCALIBRATED")

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ra_deg", "dec_deg", "x", "y", "flux_r", "flux_g", "flux_b",
                    "snr", "hfd", "n_sat", "V", "V_err", "BV", "BV_err",
                    "ref_id", "ref_V", "ref_BV", "flags"])
        for k in range(fp.n()):
            w.writerow([
                f"{fp.ra[k]:.8f}", f"{fp.dec[k]:.8f}", f"{fp.x[k]:.2f}", f"{fp.y[k]:.2f}",
                f"{fp.flux_r[k]:.1f}" if np.isfinite(fp.flux_r[k]) else "",
                f"{fp.flux_g[k]:.1f}" if np.isfinite(fp.flux_g[k]) else "",
                f"{fp.flux_b[k]:.1f}" if np.isfinite(fp.flux_b[k]) else "",
                f"{fp.snr_g[k]:.0f}", f"{fp.hfd_g[k]:.2f}", f"{fp.n_sat_px[k]:d}",
                f"{Vest[k]:.2f}" if np.isfinite(Vest[k]) else "",
                f"{V_err[k]:.4f}" if np.isfinite(V_err[k]) else "",
                f"{BVest[k]:.2f}" if np.isfinite(BVest[k]) else "",
                f"{BV_err[k]:.4f}" if np.isfinite(BV_err[k]) else "",
                tyc_name[k],
                f"{tyc_V[k]:.2f}" if np.isfinite(tyc_V[k]) else "",
                f"{tyc_BV[k]:.2f}" if np.isfinite(tyc_BV[k]) else "",
                "|".join(flags[k]),
            ])

    print(f"scope={scope} filter={filt} gain={gain} "
          f"{'set=' + cset['scope'] if cset else 'NO SET (--force)'}  "
          f"precal={precal_name} fit={cat_name}", file=sys.stderr)
    print(f"extracted={n_extracted} tycho2_matched={n_matched} calibrators_used={n_calib}",
          file=sys.stderr)
    if n_fallback:
        print(f"green-only fallback for {n_fallback} star(s) missing R or B: assumed "
              f"B-V={bv_assumed:.2f} +-{bv_spread:.2f} (field mean), adding "
              f"{abs(Tv_bv)*bv_spread:.3f} mag to their V_err", file=sys.stderr)
    print(f"colour_index={colour_index} Tbv={Tbv}+-{sigma_Tbv} Tv_bv={Tv_bv}+-{sigma_Tv_bv} k_r2={k_r2}",
          file=sys.stderr)
    print(f"fitted ZPbv={ZPbv:.3f}+-{sigma_ZPbv:.3f}  zp={zp:.3f}+-{sigma_zp:.3f} "
          f"(zero-point fit RMS {rms_zp:.3f} mag)", file=sys.stderr)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
