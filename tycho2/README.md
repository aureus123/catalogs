# Cross identifications with Tycho-2 catalogue

Tool *gen_tycho2* generates cross identifications between Tycho-2 and PPM/DM catalogues. Designation are for DM (BD/SD/CD) and other old catalogues. In the case of PPM stars, DM designation are taken from their DM column. In north hemisphere, also BD stars (from -1 to +19) are used. In south hemisphere, also CD stars are used. In the case designations for GC (Catálogo General Argentino), OA (Oeltzen-Argelander catalogue), U (Yarnall-Frisby USNO),
GZC (Gould's Zone catalogue) and G (Gilliss catalogue) are present, they are used instead of CD.

### Shortcomings of the current approach

Since all identifications are performed just by angular distance thresholds, misidentifications might happen. On the other hand, some BD/SD stars are
missing in the south hemisphere since they are exclusively taken from PPM
catalogues. The same applied for the other old catalogues GC/OA/U/GZC/G, they
are reported only if a PPM/DM star is cross-matched with them. The file
*read_old.cpp* does the cross-identifications between PPM/DM and GC/OA/U/GZC/G.

Other further improvements (not done yet) might be to use a better
cross-identification algorithm, e.g. that maximizes likelihood between positions
and magnitudes when comparing Tycho-2 and a target catalogue, with Hungarian
algorithm, see https://doi.org/10.1016/j.endm.2018.07.005. Finally, Tycho-2 supplementary stars are not considered.

### Instructions

Before executing *gen_tycho2*, download Tycho-2 catalog from VizieR (I/259) and concatenate data files into a single *cat/tyc2.txt* file.

Files in this folder:
- cross_tyc2_north.csv = Cross identifications of Tycho-2 stars (northern hemisfere)
- cross_tyc2_south.csv = Cross identifications of Tycho-2 stars (southern hemisfere)
- plotann.py = Modified Python script (part of Astrometry.net software) that replaces Tycho-2 stars by those generated in the CSV files

In order to use plotann.py, replace it. Then, copy the CSV files here to the target folder where plotann.py will be used. Finally, use plotann with Tycho2 catalogue, e.g.
    plotann.py image.wcs image.png image.ann.png --no-grid --no-const --tycho2cat=tycho2.kd

Example:
![Alt text](image.ann.png?raw=true "Southern Pleyades")

