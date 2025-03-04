# Old catalogs comparison

💻 This repository contains some simple algorithms to compare old catalogs in order to detect typo errors during the transcription from original printed catalogs to digital ones, particularly the Durchmusterungs.

✏️ The main comparison is given between Cordoba Durchmusterung (CD) and Positions and Proper Motions (PPM) catalogs, but other "experiments" are performed between other catalogs.
The comparison consists of computing the angular distance between the position of a star in both catalogs, and log those cases where it exceeds a given threshold. Also the difference between visual magnitude is logged.
The threshold used for comparing CD and PPM is 2 arcmin, and for magnitude is 1.5 units.
The hypothesis is that, since the cross-identifications between both catalogs were also transcribed from older documents, it provides an independent matching between both catalogs (i.e. the cross-identification does not comes from the comparison between the digital version of the Durchmusterung and the target catalog, but from a man-made identification between the printed version and the target catalog).

### Algorithm:

📄 The approach considers two catalogs, the source (e.g. CD) and the target (e.g PPM).
In first place, all stars from the target catalog having a cross-identification between both catalogs are considered. Then, coordinates are precessed to the epoch of the source catalog (e.g from J2000 to B1875) and corrected by proper motion.
The resulting coordinates are compared to those from the source catalog.
If the angular distance exceeds a given threshold, a warning is generated, with the data of that star in the source catalog in the format
used in its printed form. Also, the page of the printed catalog is
estimated.
One can compare these data against the real printed catalog, in order to
find typo errors during the confection of the digital version of BD or CD.

### Limitations of the approach:

🛑 On the one hand, only typo mistakes that leads to an excess in the thresholds can be corrected. For instance, a threshold of 0.5 for magnitudes will find some errors in the first digit of them (e.g. if the value reported is 8.8 but the real is 9.3) and will skip some other (e.g. if the real is 9.1), however mistakes in the second digit will not be found. Naturally, one can reduce the threshold to raise the number of hits, but at the expense of increase a lot the number of false-positives (see Experiment 2 in the results folder).

✋ On the other hand, the approach is limited to the cross-identified stars. In the case of CD, from a total of 613778 stars, only 171303 has useful cross-identifications with PPM. That means that there are roughly 72% of stars in the CD digital catalog that the algorithm does not explore.  

### Comparison schemes:

- *compare_ppm*: Compares CD and PPM (from declination -23 to south pole)
- *compare_ppm_bd*: Compares BD and PPM (only Vol 1, from declination -1 to +19)
- *compare_sd*: Compares CD and SD (declination -22), through catalog 4005
- *compare_cpd*: Compares CD and CPD, through catalogs 4005 or 4011
- *compare_agk*: Compares CD and AGK (Cordoba A, B and C, from declination -22 to -37)

### Other experiments:

- *compare_gc*: Cross-identifies CD and GC (Argentine General Catalog) and compares them (here *log_gc.log* was generated with original CD while cross CSV
files were generated with curated CD).
- *compare_cd*: Logs differences between two digital versions of CD
- *gen_tycho2_north* and *gen_tycho2_south*: See README in tycho2 folder
- *find_coord*: Tool to find by coordinates on old catalogs 

### Requirements

🚰 Sources of WCSTools 3.9.7 from http://tdc-www.harvard.edu/wcstools should be downloaded to wcstools-3.9.7 folder, and compiled.
The routine "wcsconp" is used for transforming coordinates.

### Bibliography

- Severin D. E., Sevilla D. J. (2015) Development of a new digital version of "Cordoba Durchmusterung" stellar catalog. Revista Académica Electrónica de la UNR 15 (8), pg. 2250-2260.
https://rephip.unr.edu.ar/server/api/core/bitstreams/0362c1df-b472-4216-99fe-46c7f135921a/content 
- Severin D. E. (2018) Cross-identification of stellar catalogs with multiple stars: Complexity and Resolution. Electronic Notes in Discrete Mathematics 69, pg 29-36.
https://doi.org/10.1016/j.endm.2018.07.005
- Severin, D. E. (2018) Cross-identification between Cordoba Durchmusterung catalog (declinations -22, -23 and -24) and PPMX catalog, Mendeley Data, V1.
 http://dx.doi.org/10.17632/5wwwtv7c8c.1


Enjoy! 🤗

Daniel Severin.-