# Variable star estimation

This folder contains examples of FITS images used to estimate the magnitude of variable
stars with the `cross_txt` tool.

The script `estimate.sh` shows a full end-to-end workflow: it extracts the
green channel from a post-processed FITS image produced by a Dwarf 3 sensor, solves the
astrometry, and then calls `cross_txt` in all of its sub-modes
(ensemble fit, one or two comparison stars, PPM or user-supplied V magnitudes). It serves as a reference for how to
invoke the tool in practice.

## Background

The `cross_txt` tool was originally designed to compare magnitude scales of
historical star catalogs against the PPM catalog. The same cross-matching machinery, however, can be applied to an arbitrary set
of stars captured in a CCD image.  In that context, instead of calibrating a
historical catalog, the tool uses the PPM-matched field stars to
build a photometric calibration curve and then applies it to estimate the
magnitude of a user-selected variable star.  Sequence stars with known
magnitudes (from PPM or supplied by the user) anchor the fit, and an optional
control star provides an independent quality check.

Although to obtain a good magnitude estimation, it is suggested to use specialized tools (such as [ASTAP](https://www.hnsky.org/astap.htm)) where green channel is extracted from raw FITS files and then the protometry is performed on the resulting FITS file obtained by stacking and calibrating with darks/flats, the
approach given here allows to handle an already post-processed multi-color FITS file provided by the smartscope and computes a rough estimation of the V magnitude of a star (with better results in the range 8-10).

See Requirements below for the tools and packages needed to run everything
in this folder.

## Requirements

- **[ASTAP](https://www.hnsky.org/astap.htm)** -- star detection (`-extract`)
  and a D80/D50 star database for astrometry fallback. The Tycho-2 pipeline
  invokes the binary directly; on macOS it expects
  `/Applications/ASTAP.app/Contents/MacOS/astap` (override with the
  `ASTAP_BIN` environment variable).
- **[WCSTools](http://tdc-www.harvard.edu/wcstools/)** -- `imextract` to pull
  the R/G/B planes out of the Dwarf's 3-channel FITS cube (expected at
  `wcstools-3.9.7/bin/imextract` next to this repo; override with
  `IMEXTRACT_BIN`).
- **[Astrometry.net](https://astrometry.net/)** -- `solve-field`,
  `wcs-xy2rd`, `tablist`. Only used as a fallback when a field has no
  pre-solved `.wcs`; most fields under `~/Desktop/Dwarf2`/`Dwarf3` already
  have one and skip this step entirely.
- **Python 3** (tested with 3.13) with:
  - `numpy`, `scipy`, `astropy` -- catalogue handling, WCS, matching, fitting
  - `photutils` -- aperture photometry backend (`pip install photutils`)
  - `Pillow` (`PIL`) -- `annotate.py`'s magnitude overlay

None of the above are vendored in this repo; install them before running
`tycho2.py`, `precalibrate.py`, `measure.py`, `sweep.py`, `annotate.py` or
`benchmark.py`.

## Example results

The table below summarises the estimations produced by running both examples
in `estimate.sh` against their respective solved FITS images.  "One comp"
means a constant fit anchored by a single sequence star; "Two comp" means a
linear fit anchored by two sequence stars.  The "Comp mags" column indicates whether the sequence and
control star magnitudes came from the PPM catalog or were supplied by the
user (in this case, from SIMBAD).

| Mode | Comp mags | Variable star | Est. Vmag | Control PPM | Control ref Vmag | Est. control Vmag | Control abs. error |
|------|-----------|---------------|-----------|-------------|------------------|-------------------|--------------------|
| Ensemble  | PPM    | Y Centauri | 8.41 | 263030 | 8.4  | 8.34 | 0.06 |
| One comp  | PPM    | Y Centauri | 7.94 | 263030 | 8.4  | 7.85 | 0.55 |
| Two comp  | PPM    | Y Centauri | 8.34 | 263030 | 8.4  | 8.02 | 0.38 |
| Two comp  | Custom | Y Centauri | 8.21 | 263030 | 8.15 | 7.89 | 0.26 |
| Ensemble  | PPM    | R Octantis | 8.73 | 376400 | 8.7  | 8.77 | 0.07 |
| One comp  | PPM    | R Octantis | 8.29 | 376400 | 8.7  | 8.33 | 0.37 |
| Two comp  | PPM    | R Octantis | 8.45 | 376400 | 8.7  | 8.55 | 0.15 |
| Two comp  | Custom | R Octantis | 8.36 | 376400 | 8.71 | 8.52 | 0.19 |

The ensemble fit is the one to trust.  Its RMSE against the PPM magnitudes of
the field stars is 0.35 mag for Y Centauri (71 matched stars) and 0.23 mag for
R Octantis (68 matched stars).  The two-comparison mode is fragile: when the
two sequence stars happen to have nearly equal instrumental magnitudes the
two-point fit is badly conditioned, and slopes of 2.7 to 3.9 are obtained
instead of the expected value near 1.

## Note on the extracted channel

Until 2026-09 the script passed `1` to `imextract`, believing it selected the
green plane.  It does not: `imextract` numbers planes from 1 (see the
`(nimage - 1) * nbimage` offset in `wcstools-3.9.7/imextract.c`), so `1` is the
first plane of the RGB cube, which is **red**.  Every `*-TG.fits` produced
before that date was in fact red-channel data.  The script now passes `2`.

Measured against Tycho-2 (which, unlike PPM, provides real photometry with
sigma(V_T) of 0.014 to 0.06 mag in the range 8 to 11), the plane matters a
great deal.  Using ASTAP aperture photometry on the R Centauri field:

| plane | colour coefficient vs (B-V) | RMS against Tycho-2 V |
|-------|-----------------------------|-----------------------|
| 1 (red)   | +0.27 | 0.241 |
| 2 (green) | -0.03 | 0.104 |
| 3 (blue)  | -0.34 | 0.226 |

The green plane needs essentially no colour transformation, which is what one
wants for an estimate of V.

The table above, however, barely moved when the channel was corrected: the
control-star errors changed by at most 0.06 mag, and the ensemble RMSE against
PPM did not improve.  That is not a contradiction.  Two other error sources
dominate this pipeline and both are independent of the channel:

* PPM magnitudes are quantised to 0.1 mag and are largely Durchmusterung-era
  visual and photographic estimates, with errors of 0.2 to 0.3 mag.  Fitting
  against them cannot produce a residual smaller than their own scatter.
* The fluxes come from the `flux` column of the Astrometry.net `.axy` table,
  which is `simplexy`'s peak amplitude rather than an integrated aperture
  flux.  On the same field, ASTAP's integrated photometry gives a fitted slope
  of 0.97 against a reference V, whereas the `.axy` flux gives 0.78 - the
  apparent non-linearity that the quadratic fit in `--txt` mode absorbs.

So the estimated magnitudes did shift (Y Cen 8.27 to 8.41, R Oct 8.56 to
8.73), and the green results are the correct ones to use, but reaching the
0.10 mag accuracy that the sensor is capable of needs a better reference
catalogue and real aperture photometry, not just the right plane.

## Tycho-2 / three-channel pipeline

The tools above (`cross_txt` via `estimate.sh`) are kept as-is for
historical-catalogue work and quick single-star estimates against PPM. For a
properly calibrated V and B-V over a whole field, use the newer two-tool
pipeline instead, which cross-matches all three Bayer channels against
Tycho-2 (V = 7-11.5) rather than a single green-only fit against PPM:

```sh
P=~/.venv/bin/python
$P tycho2.py --build                                     # one-time: cat/tyc2.txt -> cat/tyc2_phot.npy
$P gaia_v50.py --selftest                                # verify the ASTAP V50 database decodes
$P precalibrate.py --catalog gaia_online --filter VIS --gain 60 -v -o instrument_coeffs.json
$P measure.py <field>.fits -c instrument_coeffs.json -o <field>.csv
$P annotate.py <field>.csv --image <field>.png -o <field>.mag.png
$P sweep.py --coarse --fields 10 -o sweep_coarse.json    # optional: retune the aperture
$P benchmark.py --all                                    # loo-star/loo-field/repeats/ablation
```

`--catalog` defaults to `tycho2` (no network). See "Reference catalogues" below.

### Reference catalogues

`precalibrate.py` and `measure.py` take `--catalog`, and the choice is recorded in
`instrument_coeffs.json` so a coefficient set is never applied against a different
reference than it was fitted with:

| `--catalog` | network | notes |
|---|---|---|
| `tycho2` | no | `cat/tyc2_phot.npy`. Useful to V ≤ 11.5; σ(V) 0.015–0.09 growing with magnitude. B−V transform valid only below B−V 1.53. |
| `gaia_v50` | no | ASTAP's local V50 database. V quantised to 0.1 mag (a 0.029 mag floor), B−V to 0.02 mag but capped at ±2.54. |
| `gaia_online` | **yes** | VizieR Gaia DR3 — full precision, no colour cap, depth on demand. Cached per field under `estim/cache/online/`, but a *new* field needs a connection. |

Measured head-to-head on RCen (271 detections, 143-star common subset, identical blend
mask), robust residual σ: **0.068** (tycho2), **0.047** (gaia_v50), **0.031**
(gaia_online). The v50→online gap is exactly v50's quantisation floor —
√(0.047² − 0.031²) = 0.035 — since both use the same Gaia source and transform.

Recommended split: **`gaia_online` for precalibration** (one-time, cached, best
coefficients) and **`gaia_v50` for measurement** (fully offline; our own ~0.08 mag
photometry error dominates there anyway). Mixing those two is safe because V50 *is* the
same Gaia transform precomputed; mixing either with `tycho2` is not — their V scales
differ by about 0.04 mag.

Blend rejection is a **flux-ratio** test (`starcat.blend_dmag`): a calibrator is
rejected only if a catalogue neighbour within 20″ is brighter than V + 3. The earlier
"any neighbour within 20″" rule got stricter as the catalogue deepened — on RCen it cut
7 stars with Tycho-2 but 117 with V50 unlimited — so a better reference perversely
shrank the calibrator sample.

### Annotator markers

`annotate.py`'s markers encode two things independently: shape is whether the
star matched Tycho-2 (square = matched, circle = unmatched), and outline
colour is its `measure.py` flag, checked in priority order -- red = `SAT`
(saturated), orange = `BV_EXTRAP` (estimated B-V outside the calibrated
range), green = anything else. The yellow number is V with the decimal point
dropped (`875` = V 8.75), the standard finder-chart convention so a label is
never misread as a star. The full flag list (e.g. a star can be both `SAT`
and `EDGE`) is only in the CSV -- the marker only shows the highest-priority
one.

`precalibrate.py` fits the colour transform (`Tbv`) and magnitude transform
(`Tv_bv`, `k_r2`) once per scope from every VIS/gain-60 field in
`~/Desktop/Dwarf2` and `~/Desktop/Dwarf3` (read-only -- all intermediate
files go to `estim/cache/`), and writes them to `instrument_coeffs.json`.
`measure.py` looks up the (scope, filter, gain) triplet from a new FITS
cube's header, takes those coefficients, and only fits the two per-image
zero points (`zp`, `ZPbv`) that don't transfer between nights. Validated by
`benchmark.py`: leave-one-field-out RMS matches leave-one-star-out RMS to
within 0.003 mag for both scopes, confirming the transform really is a
hardware constant and shouldn't be refit per field. See
`.claude/plans/let-s-brainstorm-about-creating-lovely-rabbit.md` for the
full design rationale.
