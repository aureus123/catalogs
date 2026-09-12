#!/usr/bin/env python3
"""Shared photometry backend for precalibrate.py and measure.py.

Given a 3-channel FITS cube from a Dwarf 2/3 smartscope (axis order
plane, y, x; plane 0 = R, 1 = G, 2 = B -- see estim/README.md "Note on the
extracted channel"), this module:

  1. extracts the three planes with wcstools' imextract,
  2. solves (or reuses an existing) astrometric WCS on the green plane only
     -- the three planes share one pixel grid, so one WCS covers all of them,
  3. detects stars independently in each plane with ASTAP's ``-extract``,
  4. cross-matches the R and B detections onto the green (reference) list by
     pixel position,
  5. measures flux with a pluggable aperture stage: ASTAP's own flux (fixed
     aperture 99 / annulus 14, see README) or a re-measurement with
     photutils at an aperture/annulus expressed as multiples of each
     detection's own measured HFD.

The result is per-star raw measurements only (position, per-channel flux,
HFD, SNR, saturation, isolation). Fitting zero points and colour
coefficients is precalibrate.py/measure.py's job, not this module's.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
IMEXTRACT_BIN = Path(os.environ.get(
    "IMEXTRACT_BIN", HERE.parent / "wcstools-3.9.7" / "bin" / "imextract"))
ASTAP_BIN = Path(os.environ.get(
    "ASTAP_BIN", "/Applications/ASTAP.app/Contents/MacOS/astap"))

PLANE_NAMES = {1: "R", 2: "G", 3: "B"}

#: Different Dwarf firmware revisions write the same scope under different TELESCOP
#: strings -- "DWARFII", "DWARF II" and "DWARF 2" are all the same instrument, likewise
#: for the Dwarf 3. Everything downstream keys coefficient sets on (scope, filter, gain),
#: so without canonicalising this the archive splits into groups that are physically
#: identical, and 104 fields end up unmatched against their own coefficients.
SCOPE_ALIASES = {
    "DWARFII": "DWARF 2", "DWARF II": "DWARF 2", "DWARF 2": "DWARF 2", "DWARF2": "DWARF 2",
    "DWARFIII": "DWARF 3", "DWARF III": "DWARF 3", "DWARF 3": "DWARF 3", "DWARF3": "DWARF 3",
}


#: Filter names, canonicalised. The Dwarf firmware renamed its filters around
#: 2025-11: "CUT" became "VIS" and "PASS" became "Astro". Same physical filters --
#: header BAYERPAT, BITPIX and frame size are identical across the rename.
#:
#: **But the photometry is not identical.** Measured on Dwarf 2 gain-60 fields, the
#: pre-rename firmware delivers ~22% less instrumental colour separation for the same
#: real B-V range, giving Tbv = 2.214 +-0.055 against 1.344 +-0.015 after the rename --
#: a 16-sigma difference. That is on-board colour processing (white balance / saturation
#: applied when building the 3-channel cube), not the passband. Applying the new-firmware
#: coefficients to old-firmware frames degrades the V residual from 0.068 to 0.312 mag.
#:
#: So the alias is right for *identifying* the filter and wrong for *transferring*
#: coefficients; measure.py warns loudly when it relies on it.
FILTER_ALIASES = {
    "CUT": "VIS", "VIS": "VIS",
    "PASS": "Astro", "ASTRO": "Astro",
    "DUO-BAND": "Duo-Band", "DUOBAND": "Duo-Band",
}
#: Aliases that indicate a pre-rename (older firmware) capture.
LEGACY_FILTER_NAMES = {"CUT", "PASS"}


def normalize_filter(raw):
    """Canonical filter name, or '?' when unknown. Unrecognised values pass through."""
    if raw is None:
        return "?"
    key = " ".join(str(raw).strip().upper().split())
    if not key or key == "NONE":
        return "?"
    return FILTER_ALIASES.get(key, str(raw).strip())


def resolve_filter(header, cube_path):
    """Filter for a cube, falling back to ``shotsInfo.json`` when the header lacks it.

    Pre-rename firmware writes no FILTER card at all, so the only record of which filter
    was used is the ``ir`` field in the capture's shotsInfo.json.

    Returns ``(canonical, raw, source)`` -- raw is what was actually recorded, so the
    caller can tell whether a legacy alias was involved.
    """
    raw = header.get("FILTER")
    if raw and str(raw).strip().upper() not in ("", "NONE"):
        return normalize_filter(raw), str(raw).strip(), "header"
    si = Path(cube_path).parent / "shotsInfo.json"
    if si.exists():
        try:
            ir = json.loads(si.read_text()).get("ir")
        except Exception:
            ir = None
        if ir:
            return normalize_filter(ir), str(ir).strip(), "shotsInfo.json"
    return "?", None, None


def normalize_scope(raw):
    """Canonical scope name for a TELESCOP/INSTRUME header value.

    Unknown values pass through unchanged (uppercased and whitespace-collapsed) rather
    than being forced into one of the two known scopes -- a third instrument should show
    up as itself, not be silently mislabelled.
    """
    if raw is None:
        return "?"
    key = " ".join(str(raw).strip().upper().split())
    return SCOPE_ALIASES.get(key, key)


def _run(cmd, **kwargs):
    kwargs.setdefault("capture_output", True)
    kwargs.setdefault("text", True)
    r = subprocess.run(cmd, **kwargs)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n"
                            f"stdout: {r.stdout}\nstderr: {r.stderr}")
    return r


def extract_plane(cube_path, plane, out_prefix):
    """Extract one plane (1=R, 2=G, 3=B) of a FITS cube to out_prefix.fits."""
    out_prefix = Path(out_prefix)
    # imextract exits 0 but silently writes nothing when the output directory is
    # missing, which surfaces later as a confusing "did not produce" error. Create it
    # here so the function is self-sufficient rather than relying on every caller.
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    _run([str(IMEXTRACT_BIN), "-v", "-o", str(out_prefix), str(plane), str(cube_path)])
    out_fits = Path(str(out_prefix) + ".fits")
    if not out_fits.exists():
        raise RuntimeError(f"imextract did not produce {out_fits}")
    return out_fits


def resolve_wcs(cube_path, green_fits_path, wcs_path=None):
    """Prefer a WCS already solved for the cube itself (all three planes share
    one pixel grid, so it applies directly) over re-solving the green plane.
    A few older captures keep the FITS under a name that doesn't match its
    folder (e.g. V462Lup1/V462-Lup1.fits) while the .wcs still matches the
    folder name, so check both."""
    if wcs_path:
        p = Path(wcs_path)
        if not p.exists():
            raise RuntimeError(f"--wcs {p} does not exist")
        return p
    cube_path = Path(cube_path)
    candidates = [cube_path.with_suffix(".wcs"),
                  cube_path.parent / f"{cube_path.parent.name}.wcs"]
    for c in candidates:
        if c.exists():
            return c
    return ensure_green_wcs(green_fits_path)


def ensure_green_wcs(green_fits):
    """Solve astrometry for the green plane if a .wcs doesn't already exist."""
    green_fits = Path(green_fits)
    wcs_path = green_fits.with_suffix(".wcs")
    if wcs_path.exists():
        return wcs_path
    _run([
        "solve-field", "--new-fits", "none", "--match", "none", "--solved", "none",
        "--rdls", "none", "--corr", "none", "--scale-units", "arcsecperpix",
        "--scale-low", "2.76", "--scale-high", "2.79", str(green_fits), "--overwrite",
    ])
    if not wcs_path.exists():
        raise RuntimeError(f"solve-field did not produce {wcs_path}")
    return wcs_path


_ASTAP_DTYPE = np.dtype([("x", "f8"), ("y", "f8"), ("hfd", "f8"),
                          ("snr", "f8"), ("flux", "f8")])


def run_astap_extract(fits_path, snr_min=5):
    """Run ASTAP -extract on a single-plane FITS, return a structured array
    with fields x, y, hfd, snr, flux (pixel positions, FITS 1-based, same
    convention as astrometry.net -- verified to agree with a .axy list to
    better than 0.15 px RMS on ROct-TG)."""
    fits_path = Path(fits_path)
    csv_path = fits_path.with_suffix(".csv")
    _run([str(ASTAP_BIN), "-f", str(fits_path), "-extract", str(snr_min), "-log"])
    if not csv_path.exists():
        raise RuntimeError(f"astap -extract did not produce {csv_path}")
    rows = []
    with open(csv_path) as fh:
        next(fh)  # header: x,y,hfd,snr,flux,ra[0..360],dec[0..360]
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            x, y, hfd, snr, flux = (float(v) for v in parts[:5])
            rows.append((x, y, hfd, snr, flux))
    return np.array(rows, dtype=_ASTAP_DTYPE)


def photutils_flux(data2d, x, y, r_aperture, r_in, r_out):
    """Background-subtracted aperture flux at 0-indexed pixel positions (x, y).

    Positions and radii may be per-star arrays (radii scaled by each star's
    own HFD). Background is the sigma-clipped median of the annulus.
    """
    from astropy.stats import sigma_clipped_stats
    from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry

    n = len(x)
    flux = np.full(n, np.nan)
    for i in range(n):
        if not (np.isfinite(x[i]) and np.isfinite(y[i]) and np.isfinite(r_aperture[i])
                and np.isfinite(r_in[i]) and np.isfinite(r_out[i])):
            continue
        pos = (x[i], y[i])
        aperture = CircularAperture(pos, r=r_aperture[i])
        annulus = CircularAnnulus(pos, r_in=r_in[i], r_out=r_out[i])
        mask = annulus.to_mask(method="center")
        ann_data = mask.multiply(data2d)
        if ann_data is None:
            continue
        ann_1d = ann_data[mask.data > 0]
        ann_1d = ann_1d[np.isfinite(ann_1d)]
        if len(ann_1d) == 0:
            continue
        _, bkg_median, _ = sigma_clipped_stats(ann_1d, sigma=3.0)
        phot = aperture_photometry(data2d, aperture, method="exact")
        flux[i] = float(phot["aperture_sum"][0]) - bkg_median * aperture.area
    return flux


@dataclass
class FieldPhotometry:
    scope: str
    filt: str
    gain: float
    x: np.ndarray
    y: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    hfd_r: np.ndarray
    hfd_g: np.ndarray
    hfd_b: np.ndarray
    snr_r: np.ndarray
    snr_g: np.ndarray
    snr_b: np.ndarray
    flux_r: np.ndarray
    flux_g: np.ndarray
    flux_b: np.ndarray
    saturated: np.ndarray
    n_sat_px: np.ndarray
    isolation_arcsec: np.ndarray

    def n(self):
        return len(self.x)


def _cross_match_px(ref_x, ref_y, cand_x, cand_y, cand_field, radius_px):
    """For each reference (x, y), find the nearest candidate within radius_px.
    Returns arrays (matched values of cand_field, or NaN)."""
    n = len(ref_x)
    out = np.full(n, np.nan)
    if len(cand_x) == 0:
        return out
    tree = cKDTree(np.stack([cand_x, cand_y], axis=1))
    d, i = tree.query(np.stack([ref_x, ref_y], axis=1), k=1)
    good = d <= radius_px
    out[good] = cand_field[i[good]]
    return out


def _isolation_arcsec(ra, dec):
    """Distance from each star to its nearest neighbour in the same list."""
    n = len(ra)
    if n < 2:
        return np.full(n, np.inf)
    ra_r = np.radians(ra)
    dec_r = np.radians(dec)
    cd = np.cos(dec_r)
    xyz = np.stack([cd * np.cos(ra_r), cd * np.sin(ra_r), np.sin(dec_r)], axis=1)
    tree = cKDTree(xyz)
    chord, _ = tree.query(xyz, k=2)
    sep = np.degrees(2.0 * np.arcsin(np.clip(chord[:, 1] / 2.0, 0.0, 1.0))) * 3600.0
    return sep


@dataclass
class FieldDetections:
    """Result of the expensive, aperture-independent half of measurement:
    plane extraction, astrometry and ASTAP detection. Re-used by sweep.py to
    try many apertures without re-running ASTAP each time."""
    scope: str
    filt: str
    gain: float
    wcs: object
    cube_data: np.ndarray  # (3, ny, nx) float: 0=R, 1=G, 2=B
    det_r: np.ndarray
    det_g: np.ndarray
    det_b: np.ndarray


def detect_field(cube_path, workdir=None, snr_floor=3, verbose=False, wcs_path=None):
    """Extract planes, resolve WCS, and run ASTAP -extract at a low SNR floor
    so the result can be re-filtered to any snr_min >= snr_floor without
    re-running detection."""
    cube_path = Path(cube_path)
    workdir = Path(workdir) if workdir else cube_path.parent
    stem = cube_path.stem

    header = fits.getheader(cube_path)
    scope = normalize_scope(header.get("TELESCOP"))
    filt = header.get("FILTER", "?")
    gain = header.get("GAIN", -1)

    cube_data = fits.getdata(cube_path).astype(float)  # (3, ny, nx): 0=R, 1=G, 2=B
    if cube_data.ndim != 3 or cube_data.shape[0] != 3:
        raise ValueError(f"{cube_path} is not a 3-plane cube: shape {cube_data.shape}")

    plane_paths = {}
    for plane, name in PLANE_NAMES.items():
        suffix = "TG" if name == "G" else f"T{name}"
        out_prefix = workdir / f"{stem}-{suffix}"
        out_fits = Path(str(out_prefix) + ".fits")
        # Only reuse a cached plane if it really came from *this* cube. The cache is
        # keyed on the file stem, and two different captures can share one (e.g. both
        # Dwarf2/TPup and Dwarf3/TPup exist, with different sensor sizes). Reusing the
        # wrong one silently detects stars on one image and measures apertures on
        # another, which produces plausible-looking but meaningless photometry.
        if out_fits.exists():
            try:
                cached_shape = fits.getdata(out_fits).shape
            except Exception:
                cached_shape = None
            if cached_shape != cube_data.shape[1:]:
                print(f"warning: cached {out_fits.name} is {cached_shape}, but this cube's "
                      f"planes are {cube_data.shape[1:]} -- re-extracting (stale cache "
                      f"from a different image with the same name)", file=sys.stderr)
                out_fits.unlink()
                for stale in (out_fits.with_suffix(".csv"), out_fits.with_suffix(".wcs")):
                    if stale.exists():
                        stale.unlink()
        if not out_fits.exists():
            if verbose:
                print(f"extracting plane {plane} ({name}) -> {out_fits}", file=sys.stderr)
            out_fits = extract_plane(cube_path, plane, out_prefix)
        plane_paths[name] = out_fits

    wcs_file = resolve_wcs(cube_path, plane_paths["G"], wcs_path)
    wcs = WCS(fits.getheader(wcs_file))

    if verbose:
        print(f"running ASTAP -extract (snr_floor={snr_floor}) on R, G, B planes...",
              file=sys.stderr)
    det = {name: run_astap_extract(plane_paths[name], snr_floor) for name in ("R", "G", "B")}

    return FieldDetections(scope=scope, filt=filt, gain=gain, wcs=wcs, cube_data=cube_data,
                            det_r=det["R"], det_g=det["G"], det_b=det["B"])


def photometer(fd, snr_min=5, max_stars=500, min_star_size_hfd=0.0,
               aperture_hfd=2.0, annulus_in_hfd=4.0, annulus_out_hfd=6.0,
               backend="photutils", saturation_adu=64000, channel_radius_px=1.5):
    """Re-measure a FieldDetections at a chosen set of extraction hyperparameters."""
    green = fd.det_g[(fd.det_g["snr"] >= snr_min) & (fd.det_g["hfd"] >= min_star_size_hfd)]
    if len(green) == 0:
        raise RuntimeError("no green-plane detections pass snr_min/min_star_size_hfd")
    if len(green) > max_stars:
        order = np.argsort(green["flux"])[::-1][:max_stars]
        green = green[order]

    ra, dec = fd.wcs.all_pix2world(green["x"], green["y"], 1)

    hfd_r = _cross_match_px(green["x"], green["y"], fd.det_r["x"], fd.det_r["y"],
                             fd.det_r["hfd"], channel_radius_px)
    snr_r = _cross_match_px(green["x"], green["y"], fd.det_r["x"], fd.det_r["y"],
                             fd.det_r["snr"], channel_radius_px)
    flux_r_astap = _cross_match_px(green["x"], green["y"], fd.det_r["x"], fd.det_r["y"],
                                    fd.det_r["flux"], channel_radius_px)
    x_r = _cross_match_px(green["x"], green["y"], fd.det_r["x"], fd.det_r["y"],
                           fd.det_r["x"], channel_radius_px)
    y_r = _cross_match_px(green["x"], green["y"], fd.det_r["x"], fd.det_r["y"],
                           fd.det_r["y"], channel_radius_px)

    hfd_b = _cross_match_px(green["x"], green["y"], fd.det_b["x"], fd.det_b["y"],
                             fd.det_b["hfd"], channel_radius_px)
    snr_b = _cross_match_px(green["x"], green["y"], fd.det_b["x"], fd.det_b["y"],
                             fd.det_b["snr"], channel_radius_px)
    flux_b_astap = _cross_match_px(green["x"], green["y"], fd.det_b["x"], fd.det_b["y"],
                                    fd.det_b["flux"], channel_radius_px)
    x_b = _cross_match_px(green["x"], green["y"], fd.det_b["x"], fd.det_b["y"],
                           fd.det_b["x"], channel_radius_px)
    y_b = _cross_match_px(green["x"], green["y"], fd.det_b["x"], fd.det_b["y"],
                           fd.det_b["y"], channel_radius_px)

    if backend == "astap":
        flux_r = flux_r_astap
        flux_g = green["flux"].astype(float)
        flux_b = flux_b_astap
    elif backend == "photutils":
        r_ap_g = aperture_hfd * green["hfd"]
        r_in_g = annulus_in_hfd * green["hfd"]
        r_out_g = annulus_out_hfd * green["hfd"]
        flux_g = photutils_flux(fd.cube_data[1], green["x"] - 1, green["y"] - 1,
                                 r_ap_g, r_in_g, r_out_g)

        r_ap_r = aperture_hfd * hfd_r
        r_in_r = annulus_in_hfd * hfd_r
        r_out_r = annulus_out_hfd * hfd_r
        flux_r = photutils_flux(fd.cube_data[0], x_r - 1, y_r - 1,
                                 r_ap_r, r_in_r, r_out_r)

        r_ap_b = aperture_hfd * hfd_b
        r_in_b = annulus_in_hfd * hfd_b
        r_out_b = annulus_out_hfd * hfd_b
        flux_b = photutils_flux(fd.cube_data[2], x_b - 1, y_b - 1,
                                 r_ap_b, r_in_b, r_out_b)
    else:
        raise ValueError(f"unknown backend {backend!r}")

    # Count clipped pixels, don't just flag them: the photometric damage scales with how
    # many pixels hit the ceiling. Measured against Gaia over 10 Dwarf 3 fields, the
    # magnitude is biased faint (clipped pixels lose flux) by +0.16 mag at one saturated
    # pixel rising to +0.92 at 50+, so the count is what measure.py needs to inflate
    # V_err sensibly.
    ny, nx = fd.cube_data.shape[1:]
    n_sat_px = np.zeros(len(green), dtype=int)
    xi = np.round(green["x"] - 1).astype(int)
    yi = np.round(green["y"] - 1).astype(int)
    # The search box scales with the star's own measured HFD rather than being a fixed
    # 5x5. A saturated star is bloated -- its flat-topped core spreads out and its
    # centroid drifts -- so a fixed box misses the clipping precisely on the stars that
    # have it. T Pup (HFD 11.0 px against 7.1 typical for its field) has its clipped
    # pixels 3.0-5.1 px from the centroid, while a 5x5 box reaches only 2.8 px
    # diagonally: it counted 0 where an 11x11 box counts 11, so the star reported
    # +-0.01 mag on a measurement biased ~0.19 mag faint.
    half = np.maximum(2, np.ceil(np.nan_to_num(green["hfd"], nan=2.0))).astype(int)
    for k in range(len(green)):
        h = int(half[k])
        x0, x1 = max(0, xi[k] - h), min(nx, xi[k] + h + 1)
        y0, y1 = max(0, yi[k] - h), min(ny, yi[k] + h + 1)
        box = fd.cube_data[:, y0:y1, x0:x1]
        n_sat_px[k] = int(np.count_nonzero(box >= saturation_adu))
    saturated = n_sat_px > 0

    isolation = _isolation_arcsec(ra, dec)

    return FieldPhotometry(
        scope=fd.scope, filt=fd.filt, gain=fd.gain,
        x=green["x"], y=green["y"], ra=ra, dec=dec,
        hfd_r=hfd_r, hfd_g=green["hfd"].astype(float), hfd_b=hfd_b,
        snr_r=snr_r, snr_g=green["snr"].astype(float), snr_b=snr_b,
        flux_r=flux_r, flux_g=flux_g, flux_b=flux_b,
        saturated=saturated, n_sat_px=n_sat_px, isolation_arcsec=isolation,
    )


def measure_field(cube_path, workdir=None, snr_min=5, max_stars=500, min_star_size_hfd=0.0,
                   aperture_hfd=2.0, annulus_in_hfd=4.0, annulus_out_hfd=6.0,
                   backend="photutils", saturation_adu=64000, channel_radius_px=1.5,
                   verbose=False, wcs_path=None):
    fd = detect_field(cube_path, workdir=workdir, snr_floor=snr_min, verbose=verbose,
                       wcs_path=wcs_path)
    return photometer(fd, snr_min=snr_min, max_stars=max_stars, min_star_size_hfd=min_star_size_hfd,
                       aperture_hfd=aperture_hfd, annulus_in_hfd=annulus_in_hfd,
                       annulus_out_hfd=annulus_out_hfd, backend=backend,
                       saturation_adu=saturation_adu, channel_radius_px=channel_radius_px)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cube", help="3-channel FITS cube from Dwarf 2/3")
    ap.add_argument("-o", "--out", help="write per-star CSV here (default: stdout)")
    ap.add_argument("--workdir", help="directory for extracted planes (default: next to cube)")
    ap.add_argument("--snr-min", type=float, default=5)
    ap.add_argument("--max-stars", type=int, default=500)
    ap.add_argument("--aperture-hfd", type=float, default=2.0)
    ap.add_argument("--annulus-in-hfd", type=float, default=4.0)
    ap.add_argument("--annulus-out-hfd", type=float, default=6.0)
    ap.add_argument("--backend", choices=["photutils", "astap"], default="photutils")
    ap.add_argument("--saturation-adu", type=float, default=64000)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    fp = measure_field(args.cube, workdir=args.workdir, snr_min=args.snr_min, max_stars=args.max_stars,
                        aperture_hfd=args.aperture_hfd, annulus_in_hfd=args.annulus_in_hfd,
                        annulus_out_hfd=args.annulus_out_hfd, backend=args.backend,
                        saturation_adu=args.saturation_adu, verbose=args.verbose)

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    w = csv.writer(out)
    w.writerow(["x", "y", "ra", "dec", "hfd_r", "hfd_g", "hfd_b",
                "snr_r", "snr_g", "snr_b", "flux_r", "flux_g", "flux_b",
                "saturated", "isolation_arcsec"])
    for k in range(fp.n()):
        w.writerow([fp.x[k], fp.y[k], fp.ra[k], fp.dec[k],
                    fp.hfd_r[k], fp.hfd_g[k], fp.hfd_b[k],
                    fp.snr_r[k], fp.snr_g[k], fp.snr_b[k],
                    fp.flux_r[k], fp.flux_g[k], fp.flux_b[k],
                    fp.saturated[k], fp.isolation_arcsec[k]])
    if args.out:
        out.close()
    print(f"scope={fp.scope} filter={fp.filt} gain={fp.gain} n_stars={fp.n()}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
