# Cross identifications with Tycho-2 catalogue

Tool *gen_tycho2* generates cross identifications between Tycho-2 and PPM/DM catalogues. Designation are for DM (BD/SD/CD) and other old catalogues. In the case of PPM stars, DM designation are taken from their DM column. In north hemisphere, also BD stars (from -1 to +19) are used. In south hemisphere, also CD stars are used. In the case designations for GC (Catálogo General Argentino), OA (Oeltzen-Argelander catalogue), U (Yarnall-Frisby USNO),
GZC (Gould's Zone catalogue) or G (Gilliss catalogue) are present, they are used instead of CD.

### Shortcomings of the current approach

Since all identifications are performed just by angular distance thresholds, misidentifications might happen. On the other hand, some BD/SD stars are
missing in the south hemisphere; the same happens for the other old
catalogues GC/OA/U/GZC/G, they are reported only if a PPM/CD star is
cross-matched with them. The file *read_old.cpp* does the cross-identifications
between PPM/DM and GC/OA/U/GZC/G.

Other further improvements (not done yet) might be to use a better
cross-identification algorithm, e.g. that maximizes likelihood between positions
and magnitudes when comparing Tycho-2 and a target catalogue, with Hungarian
algorithm, see this [paper](https://doi.org/10.1016/j.endm.2018.07.005). Finally, Tycho-2 supplementary stars are not considered.

### Instructions

Before executing *gen_tycho2*, download Tycho-2 catalog from VizieR (I/259) and concatenate data files into a single *cat/tyc2.txt* file.

Files in this folder:
- cross_tyc2_north.csv = Cross identifications of Tycho-2 stars (northern hemisfere)
- cross_tyc2_south.csv = Cross identifications of Tycho-2 stars (southern hemisfere)
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
| HD 93540 | GC 14764 | 8965-1767-1 | Resultados del Observatorio Nacional Argentino, [Vol XIV](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0014&type=SCREEN_THMB), pag. 291 |
| HD 93269 | GZC 10h 2923 | 8965-288-1 | Resultados del Observatorio Nacional  Argentino, [Vol VII](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0007&type=SCREEN_THMB), pag. 337 |
| HD 93600 | G 7322 | 8965-547-1 | [A catalogue of 16748 southern stars](https://archive.org/details/catalogueof1674800unitrich/catalogueof1674800unitrich/) deduced by the United States Naval Observatory from the zone observations made at Santiago de Chile, pag. 184 |

Another example:
![Alt text](51Hya.png?raw=true "51 Hya")

Stars that you can see here:
| Star (SIMBAD) | Old designation | Tycho-2 | Where to find it |
| --- | --- | --- | --- |
| k Hya | GC 19455 | 6740-785-1 | Resultados del Observatorio Nacional Argentino, [Vol XIV](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0014&type=SCREEN_THMB), pag. 377 |
| HD 125995 | U 6039 | 6740-263-1 | Catalogue of stars observed at the [United States Naval Observatory](https://archive.org/details/cataloguestarsus00unitrich/cataloguestarsus00unitrich/) during the years 1845 to 1877, pag. 153 |
| HD 126349 | OA 13607 | 6740-108-1 | [Argelander's Zonen-Beobachtungen vom 15. bis 31.](https://babel.hathitrust.org/cgi/pt?id=uc1.$b524535&seq=278) Grade südlicher Declination, in mittleren Positionen für 1850.0, pag. 206 |
| CD -27 9807 | CD -27°9807 | 6740-336-1 | Resultados del Observatorio Nacional Argentino, [Vol XVI](https://articles.adsabs.harvard.edu/cgi-bin/iarticle_query?journal=RNAO.&volume=0016&type=SCREEN_THMB), pag. 323 |
