# Experiments and results (by Daniel Severin)

This folder has a copy of different digital catalogs.

### Efficiency of algorithm

An experiment was carried out to analize the efficiency in finding errors in CD digital catalog. Two comparisons are made:

- *Compare files cd_vol1 and I88*: In this case, the first volume of the latest CD catalog was compared against I88, an older one having typo mistakes. By comparing both files (see *compare_cd* tool), 205 differences were found: 202 errors in Position and 3 in Magnitude. Also I88 has 44 unknown codes in the magnitude field (probably also typo errors) that was not taken into account. When comparing cd_vol1 against PPM via *compare_ppm* it logs 84 warnings (72 in Position and 12 in Magnitude), whereas the comparison between I88 and PPM generates 138 warnings (123 in Position and 15 in Magnitude) meaning that using the older catalog adds 51 warnings in Position and 3 in Magnitude. Therefore, the algorithm is able to reach all magnitude differences, but for the position it only hits roughly 25% of the cases (51 out of 202).

- *Compare files cd and 1114*: Here, the latest catalog was compared against an older version from NASA-ADC. By comparing both files, 141 differences in Position were found. When comparing both catalogs against PPM, for "cd" (newer) it yields 121 in Position whereas for "1114" (older) it yields 156. It means that using the older catalog adds 35 warnings. Therefore, the algorithm is able to reach roughly 25% of the cases (35 out of 141).

The thresholds used in both experiments was 2 arcmin for the position and 1.5 for the magnitude.

### Typo errors found in Bonner Durchmusterung digital catalog

By using tool *compare_ppm_bd* with 180 arcsec and 0.4 magnitude units as thresholds, a total of 86 warnings were found (82 in position and 4 in magnitude). By manually comparing those registers with the printed catalog, 13 typo errors in position were found (see file *log_ppm_bd*). It gives a success rate of 15%. If the position threshold is reduced to 150 arcsec, 61 new warnings are generated but almost all of them are false-positives while only 1 typo error is found.

Since the RSME of distance is 36 arcsec, a threshold of 150 arcsec represents roughly 4 standard deviations. In the case of magnitude, the RSME is 0.1 so a threshold of 0.4 represents also 4 standard deviations, although none typo errors were found for that parameter.

Below is the list of corrections:

| Star | Parameter | current | corrected | dist1 | dist2 |
| --- | --- | --- | --- | --- | --- |
| BD +13°2074 | RAm | 12 | 13 | 900 | 866 |
| BD +19°4871 | DEm | 37.9 | 7.9 | 1800 | 1848 |
| BD +4°1318 | DEm | 57.5 | 37.5 | 1200 | 1180 |
| BD +2°1367 | DEm | 50.9 | 55.9 | 300 | 275 |
| BD +1°2078 | DEm | 39.3 | 29.3 | 600 | 597 |
| BD +1°2911 | DEm | 35.3 | 55.3 | 1200 | 1220 |
| BD +9°3277 | DEm | 13.5 | 18.5 | 300 | 310 |
| BD +0°3742 | DEm | 30.2 | 39.2 | 540 | 543 |
| BD +0°3924 | DEm | 53.9 | 58.9 | 300 | 324 |
| BD +3°3794 | DEm | 57.4 | 27.4 | 1800 | 1800 |
| BD +7°4658 | DEm | 8.1 | 38.1 | 1800 | 1783 |
| BD -1°2359 | RAm | 5 | 4 | 900 | 913 |
| BD -0°2683 | DEm | 5.1 | 3.1 | 120 | 222 |
| BD +0°2210 | RAs | 40.8 | 49.8 | 135 | 159 |

where the parameter may be "RAm" or "RAs" for the minutes or seconds of Right Ascension resp. or "DEm" for the minutes of Declination, *dist1* is the distance produced by the difference between the current position and the corrected one, and *dist2* is the angular distance computed by the algorithm between the position in the BD digital catalog and PPM catalog (after precession and correction by proper motion of the PPM star).

Note that *dist1* and *dist2* are similar, thus showing that the cause of the mistake in the digital catalog is typographical.

### Typo errors found in Cordoba Durchmusterung digital catalog

🚧 In construction