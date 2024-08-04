# Old catalogs comparison

This repository contains some simple algorithms to compare old catalogs in order to detect typo errors during the transcription from original printed catalogs to digital ones, particularly the Durchmusterungs.

The main comparison is given between Cordoba Durchmusterung (CD) and Positions and Proper Motions (PPM) catalogs, but other "experiments" are performed between other catalogs.
The comparison consists of computing the angular distance between the position of a star in both catalogs, and log those cases where it exceeds a given threshold. Also the difference between visual magnitude is logged.
The threshold used for comparing CD and PPM is 2 arcmin, and for magnitude is 1.5 units.
The hypothesis is that, since the cross-identifications between both catalogs were also transcribed from older documents, it provides an independent matching between both catalogs (i.e. the cross-identification does not comes from the comparison between the digital version of the Durchmusterung and the target catalog, but from a man-made identification between the printed version and the target catalog).

### Algorithm:

There are two catalogs, the source (e.g. CD) and the target (e.g PPM).
In first place, all stars from PPM having a cross-identification between both catalogs are considered. Then, coordinates from PPM are precessed
to the epoch of the CD catalog (from J2000 to B1875) and corrected by
proper motion. The resulting coordinates are compared against CD ones.
If the angular distance exceeds a given threshold, a warning is generated, with the data of that star in the CD catalog in the format
used in its printed form. Also, the page of the printed catalog is
estimated.
One can compare these data against the real printed catalog, in order to
find typo errors during the confection of the digital version of CD.  

### Comparison schemes:

- compare_ppm: Compares CD and PPM (from declination -23 to south pole)
- compare_ppm_bd: Compares BD and PPM (only Vol 1, from declination -1 to +19)
- compare_sd: Compares CD and SD (declination -22)
- compare_cpd: Compares CD and CPD, through catalogs 4005 or 4011
- compare_agk: Compares CD and AGK (Cordoba A, B and C, from declination -22 to -37)

### Other experiments:

- compare_gc: Cross-identifies CD and GC (Argentine General Catalog) and compares them
- compare_cd: Logs differences between two digital versions of CD

### Requirements

Sources of WCSTools 3.9.7 from http://tdc-www.harvard.edu/wcstools should be downloaded to wcstools-3.9.7 folder, and compiled.
The routing "wcsconp" is used for transforming coordinates.

