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

WCSTools and Astrometry.net packages are required.

## Example results

The table below summarises the estimations produced by running both examples
in `estimate.sh` against their respective solved FITS images.  "One comp"
means a constant fit anchored by a single sequence star; "Two comp" means a
linear fit anchored by two sequence stars.  The "Comp mags" column indicates whether the sequence and
control star magnitudes came from the PPM catalog or were supplied by the
user (in this case, from SIMBAD).

| Mode | Comp mags | Variable star | Est. Vmag | Control PPM | Control ref Vmag | Est. control Vmag | Control abs. error |
|------|-----------|---------------|-----------|-------------|------------------|-------------------|--------------------|
| Ensemble  | PPM    | Y Centauri | 8.27 | 263030 | 8.4  | 8.33 | 0.07 |
| One comp  | PPM    | Y Centauri | 7.82 | 263030 | 8.4  | 7.88 | 0.52 |
| Two comp  | PPM    | Y Centauri | 7.85 | 263030 | 8.4  | 7.99 | 0.41 |
| Two comp  | Custom | Y Centauri | 7.73 | 263030 | 8.15 | 7.87 | 0.28 |
| Ensemble  | PPM    | R Octantis | 8.56 | 376400 | 8.7  | 8.74 | 0.04 |
| One comp  | PPM    | R Octantis | 8.22 | 376400 | 8.7  | 8.39 | 0.31 |
| Two comp  | PPM    | R Octantis | 8.24 | 376400 | 8.7  | 8.57 | 0.13 |
| Two comp  | Custom | R Octantis | 8.01 | 376400 | 8.71 | 8.55 | 0.16 |
