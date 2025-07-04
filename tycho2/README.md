# Cross identifications with Tycho-2 catalogue

Tool *gen_tycho2* generates cross identifications between Tycho-2 and PPM/DM catalogues, however the DM designation is used instead of PPM. In the case designations of other old catalogues were available, the latter are reported. In the case of PPM stars, DM designation are taken from their DM column. In northern hemisphere, stars are reported in the following order of priority:
- Stars from Uranometría Argentina (*UA*)
- Yarnall-Frisby USNO catalog (*U*)
- *BD* designation coming from PPM stars
- All Bonner Durchmusterung stars up to +19 declination (first volume)

In southern hemisphere, we have stars from:
- Uranometría Argentina (*UA*)
- Catálogo General Argentino (*GC*)
- Gould's Zone Catalog designiations coming from other catalogs (*ZC*)
- Oeltzen-Argelander designations coming from Weiss catalog (*OA*)
- Yarnall-Frisby USNO catalog (*U*)
- Gilliss catalog (*G*)
- Segundo Catálogo General Argentino (*G2*)
- *BD* / *SD* / *CD* designation coming from PPM stars
- All Cordoba Durchmusterung stars
- All Cape Photographic Durchmusterung stars (*CPD*)

For the northern hemisphere, there are 987 UA, 4194 U and 211760 BD stars summarizing 216941 identifications.
For the southern hemisphere, there are 7043 UA, 27340 GC, 5213 ZC, 13347 OA, 1360 U, 9091 G, 2593 G2, 6353 BD, 77202 SD, 495897 CD and 79701 CPD stars summarizing  725140 identifications.


### Shortcomings of the current approach

Since all identifications are performed just by angular distance thresholds, misidentifications might happen. On the other hand, some BD/SD stars are
missing since they are only obtained from PPM (in northern hemisphere, it happens from declination +20 and above; in southern hemisphere, it happens from declination -22 and above).
Files *cross_gc.cpp* and *cross_south.cpp* do the cross-identifications between PPM/CD/CPD and UA/GC/ZC/OA/U/G/G2.

Other further improvements (not done yet) might be to use a better
cross-identification algorithm, e.g. that maximizes likelihood between positions
and magnitudes when comparing Tycho-2 and a target catalogue, with Hungarian
algorithm, see this [paper](https://doi.org/10.1016/j.endm.2018.07.005).

### Instructions

Before executing *gen_tycho2*, download Tycho-2 catalog from VizieR (I/259) and concatenate data files into a single *cat/tyc2.txt* file. Do the same for both supplemental catalogs, into a single *cat/tyc2_suppl.txt* file.

Files in this folder:
- cross_tyc2_north.csv = Cross identifications of Tycho-2 stars (northern hemisphere)
- cross_tyc2_south.csv = Cross identifications of Tycho-2 stars (southern hemisphere)
- cross_tyc2_south_*.csv = See next section
- colored.sh = See next section
- plotann.py = Modified Python script (part of Astrometry.net software) that replaces Tycho-2 stars by those generated in the CSV files

In order to use plotann.py, replace it. Then, copy the CSV files here to the target folder where plotann.py will be used. Finally, use plotann with Tycho2 catalogue, e.g.
```
plotann.py image.wcs image.png image.ann.png --no-grid --no-const --tycho2cat=tycho2.kd
```

Example:
![Alt text](C102.png?raw=true "Southern Pleyades")

Some stars that you can see in the example:
| Star (SIMBAD) | Old designation | Tycho-2 | Where to find it |
| --- | --- | --- | --- |
| HD 93540 | 236G Car | 8965-1767-1 | Resultados del Observatorio Nacional Argentino, [Vol I](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0001&type=SCREEN_THMB), pag. 143 |
| HD 308015 | GC 14766 | 8965-137-1 | Resultados del Observatorio Nacional Argentino, [Vol XIV](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0014&type=SCREEN_THMB), pag. 291 |
| HD 93269 | ZC 10h 2923 | 8965-288-1 | Resultados del Observatorio Nacional  Argentino, [Vol VII](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0007&type=SCREEN_THMB), pag. 337 |
| HD 93600 | G 7322 | 8965-547-1 | [A catalogue of 16748 southern stars](https://archive.org/details/catalogueof1674800unitrich/catalogueof1674800unitrich/) deduced by the United States Naval Observatory from the zone observations made at Santiago de Chile, pag. 184 |
| HD 308005 | CPD -62°1747 | 8961-2302-1 | Annals of the Cape Observatory, [Vol. 5](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=AnCap&volume=0005&type=SCREEN_THMB), pag. 357 |

Another example:
![Alt text](51Hya.png?raw=true "51 Hya")

Stars that you can see here:
| Star (SIMBAD) | Old designation | Tycho-2 | Where to find it |
| --- | --- | --- | --- |
| k Hya | k Hya | 6740-785-1 | Resultados del Observatorio Nacional Argentino, [Vol I](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0001&type=SCREEN_THMB), pag. 195 |
| HD 126088 | G2 3591 | 6740-8-1 | Resultados del Observatorio Nacional Argentino, [Vol XX](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0020&type=SCREEN_THMB), pag. 66 |
| HD 125995 | U 6039 | 6740-263-1 | Catalogue of stars observed at the [United States Naval Observatory](https://archive.org/details/cataloguestarsus00unitrich/cataloguestarsus00unitrich/) during the years 1845 to 1877, pag. 153 |
| HD 126349 | OA 13607 | 6740-108-1 | [Argelander's Zonen-Beobachtungen vom 15. bis 31.](https://babel.hathitrust.org/cgi/pt?id=uc1.$b524535&seq=278) Grade südlicher Declination, in mittleren Positionen für 1850.0, pag. 206 |
| CD -27 9807 | CD -27°9807 | 6740-336-1 | Resultados del Observatorio Nacional Argentino, [Vol XVI](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0016&type=SCREEN_THMB), pag. 323 |

### Colored and double CD stars

Based on extra data entered by hand (see folder *cd*), it is possible to flag those CD stars marked as color or double. Three additional identifications were performed:
- cross_tyc2_south_plain.csv = Cross identifications between Tycho-2 and unflagged CD stars
- cross_tyc2_south_color.csv = Cross identifications between Tycho-2 and single colored CD stars
- cross_tyc2_south_dpl.csv = Cross identifications between Tycho-2 and double CD stars (regardless if they are colored or not)

The script *colored.sh* uses *plotann.py* in a way these stars are highlighted with different colors according to their flags: green (unflagged), magenta (double) and red (color).
Currently, only stars of the first volume of CD are considered (declinations -22 to -31) and flags are available only for declinations -22 to -24. The lists are in the folder *cd*.

Example:
![Alt text](NGC4993.png?raw=true "NGC4993")

Here, it is shown that the star CD -22°9774 is colored (psi Hydrae). Also the star CD -22°9808 is double.
