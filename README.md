# Old catalogs comparison

💻 This repository contains some simple algorithms to compare old catalogs in order to detect typo errors during the transcription from original printed catalogs to digital ones, particularly the Durchmusterungs.
Catalogs are in the [cat](cat) folder.

✏️ The main comparison is given between Cordoba Durchmusterung (CD) and Positions and Proper Motions (PPM) catalogs, but other "experiments" are performed between other catalogs.
The comparison consists of computing the angular distance between the position of a star in both catalogs, and log those cases where it exceeds a given threshold. Also the difference between visual magnitude is logged.
The threshold used for comparing CD and PPM is 2 arcmin, and for magnitude is 1.5 units.
The hypothesis is that, since the cross-identifications between both catalogs were also transcribed from older documents, it provides an independent matching between both catalogs (i.e. the cross-identification does not comes from the comparison between the digital version of the Durchmusterung and the target catalog, but from a man-made identification between the printed version and the target catalog).

### Algorithm

📄 The approach considers two catalogs, the source (e.g. CD) and the target (e.g PPM).
In first place, all stars from the target catalog having a cross-identification between both catalogs are considered. Then, coordinates are precessed to the epoch of the source catalog (e.g from J2000 to B1875) and corrected by proper motion.
The resulting coordinates are compared to those from the source catalog.
If the angular distance exceeds a given threshold, a warning is generated, with the data of that star in the source catalog in the format
used in its printed form. Also, the page of the printed catalog is
estimated.
One can compare these data against the real printed catalog, in order to
find typo errors during the confection of the digital version of BD or CD.
Results and logs are in the [results](results) folder.

### Limitations of the approach

🛑 On the one hand, only typo mistakes that leads to an excess in the thresholds can be corrected. For instance, a threshold of 0.5 for magnitudes will find some errors in the first digit of them (e.g. if the value reported is 8.8 but the real is 9.3) and will skip some other (e.g. if the real is 9.1), however mistakes in the second digit will not be found. Naturally, one can reduce the threshold to raise the number of hits, but at the expense of increase a lot the number of false-positives (see Experiment 2 in the results folder).

✋ On the other hand, the approach is limited to the cross-identified stars. In the case of CD, from a total of 613778 stars, only 171303 has useful cross-identifications with PPM. That means that there are roughly 72% of stars in the CD digital catalog that the algorithm does not explore.  

### Comparison schemes

- *compare_ppm*: Compares CD and PPM (from declination -23 to south pole)
- *compare_ppm_bd*: Compares BD and PPM (only Vol 1, from declination -1 to +19)
- *compare_sd*: Compares CD and SD (declination -22), through catalog 4005
- *compare_cpd*: Compares CD and CPD, through catalogs 4005 or 4011
- *compare_agk*: Compares CD and AGK (Cordoba A, B and C, from declination -22 to -37)

### Other experiments

- *cross_gc*: Cross-identifies curated CD and GC (Argentine General Catalog) and compares them
- *cross_north*: Cross-identified lower hierarchy catalogs, mostly north.
- *cross_south*: Cross-identified lower hierarchy catalogs, mostly south.
- *compare_cd*: Logs differences between two digital versions of CD
- *gen_tycho2_north* and *gen_tycho2_south*: See README in [tycho2](tycho2) folder, also see the [gallery](gallery) folder

### Requirements

🚰 Sources of WCSTools 3.9.7 from http://tdc-www.harvard.edu/wcstools should be downloaded to wcstools-3.9.7 folder, and compiled.
The routine "wcsconp" is used for transforming coordinates.

### Mean accuracy of catalogues

When comparing different old catalogues with PPM, an accuracy can be computed
by assuming that the PPM star, corrected by proper motion, has the real coordinates
in the target epoch (19th. century). Below there are two tables with a row
per catalog: 2nd. column is the epoch of the target catalog, 3rd. column
is the number of stars cross-identified between PPM and the target catalog, and 4th. column reports the accuracy in arcseconds.

In this table, identifications are direct from PPM source:
| Catalog | Epoch | Stars | RSME |
| --- | --- | --- | --- |
| Bonner Durchmusterung (BD, 1st. vol) | 1855 | 63034 | 35.99 |
| Cordoba Durchmusterung (CD) | 1875 | 172068 | 20.48 |
| Cordoba Durchmusterung (CD, 1st. vol) | 1875 | 40832 | 24.00 |

In this table, identifications arise from small angular distances
(below a given threshold) between PPM stars and the target catalog:
| Catalog | Epoch | Stars | RSME |
| --- | --- | --- | --- |
| Uranometría Argentina | 1875 | 7605 | 5.45 |
| Catálogo General Argentino (GC) | 1875 | 31706 | 1.96 |
| Segundo Cat. General Argentino | 1900 | 5357 | 2.84 |
| Resultados XV, pg. 55-74 | 1881 | 739 | 2.52 |
| Resultados XV, pg. 140-164 | 1882 | 1259 | 2.43 |
| Resultados XV, pg. 181-184 | 1883 | 160 | 2.62 |
| Resultados XV, pg. 232-249 | 1884 | 658 | 1.99 |
| Lalande | 1800 | 46914 | 6.88 |
| Taylor | 1835 | 10944 | 4.19 |
| Oeltzen-Argelander (North) | 1842 | 25947 | 4.77 |
| Weiss | 1850 | 18010 | 4.66 |
| Weisse | 1825 | 30933 | 6.50 |
| Stone | 1880 | 12414 | 1.48 |
| USNO 3rd. edition | 1860 | 10606 | 2.77 |
| Gilliss | 1850 | 15792 | 2.87 |
| British Association Catalogue | 1850 | 8165 | 9.41 |

### Wishlist

- Use cross-identification algorithm based on a matching that maximizes likelihood (see bibliography) instead of simple angular distance thresholds.
- Write a list of all double stars from footnotes of Resultados del Observatorio Nacional Argentino, [Vol XVI](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0016&type=SCREEN_THMB) (only declinations -22, -23 and -24 were transcripted at the moment) and perform a cross-identification of that CD volume with a modern catalog.
- Correct typo error of all CD catalog.
- Usually, the supplementary letter of BD/SD/CD designations is ignored. Revise it.
- Refactor code! Consider Jupyter Lab and AstroPy instead of C++ code.

### Bibliography

- Severin D. E., Sevilla D. J. (2015) Development of a new digital version of "Cordoba Durchmusterung" stellar catalog. Revista Académica Electrónica de la UNR 15 (8), pg. 2250-2260.
https://rephip.unr.edu.ar/server/api/core/bitstreams/0362c1df-b472-4216-99fe-46c7f135921a/content 
- Severin D. E. (2018) Cross-identification of stellar catalogs with multiple stars: Complexity and Resolution. Electronic Notes in Discrete Mathematics 69, pg 29-36.
https://doi.org/10.1016/j.endm.2018.07.005
- Severin, D. E. (2018) Cross-identification between Cordoba Durchmusterung catalog (declinations -22, -23 and -24) and PPMX catalog, Mendeley Data, V1.
 http://dx.doi.org/10.17632/5wwwtv7c8c.1


Enjoy! 🤗

Daniel Severin.-