# Advances in Correcting Typo Errors in Old Digital Star Catalogs

**Daniel E. Severin**

*Departamento de Matemática, Escuela de Formación Básica, FCEIA, Universidad Nacional de Rosario.*

*Presented at JEHA III, November 2025.*

---

## Abstract

Nineteenth-century astrometric surveys produced a series of large printed star catalogs covering the entire celestial sphere, and these catalogs constitute a fundamental heritage of astronomy. Their digitization, carried out by an international consortium of data centers between the 1970s and 1980s, inevitably introduced typographical (typo) errors due to the sheer volume of data transcribed.
This paper describes a method for detecting and correcting such errors in the digital version of the *Córdoba Durchmusterung* and, to a lesser extent, in the *Bonner Durchmusterung*, *Uranometría Argentina*, *Argentine General* Catalog and other nineteenth-century catalogs. The main technique consists of exploiting pre-existing, independently obtained cross-identifications between the target catalog and a modern reference catalog, comparing positions and magnitudes after precessing to the appropriate epoch, and flagging discrepancies that exceed empirically tuned thresholds.
For the first volume of CD (declinations −22° to −31°), the algorithm flagged 84 candidates, of which 70 were genuine typo corrections. In total, 247 individual corrections have been made across nine digitized catalogs, and curated versions are made publicly available.

**Keywords:** astrometry; historical catalogs; star catalogs; Córdoba Durchmusterung; data quality; cross-identification.

---

## Resumen

Los catálogos astrométricos del siglo XIX fueron publicados en forma impresa y constituyen un patrimonio fundamental de la astronomía. Su digitalización, realizada por un consorcio de centros de datos durante los años 70 y 80, introdujo inevitablemente errores tipográficos (typo) debido al enorme volumen de datos a transcribir.
Este artículo describe un método para detectar y corregir dichos errores en la versión digital del *Córdoba Durchmusterung* y, en menor medida, en el *Bonner Durchmusterung*, *Uranometría Argentina*, *Catálogo General Argentino* y otros catálogos del siglo XIX. La técnica principal consiste en explotar identificaciones cruzadas preexistentes entre el catálogo objetivo y un catálogo de referencia moderno, comparando posiciones y magnitudes tras precesar a la época adecuada y señalando las discrepancias que superan umbrales calibrados empíricamente.
Para el primer tomo del CD (declinaciones −22° a −31°), el algoritmo marcó 84 candidatos, de los cuales 70 resultaron ser correcciones de typo genuinas. En total, se realizaron 247 correcciones individuales en nueve catálogos digitalizados, y las versiones curadas están disponibles públicamente.

**Palabras clave:** astrometría; catálogos históricos; catálogos estelares; Córdoba Durchmusterung; calidad de datos; identificación cruzada.

---

## 1. Introduction

A *typo error* (from *typographical error*) in the context of digital star catalogs refers to a mistake introduced during the manual transcription of data from a printed source to a machine-readable format. Such errors include digit transpositions (e.g., writing 39.1 instead of 29.1 for an arcminute value), swapped rows, and misread characters that arise from the visual similarity between certain printed glyphs (1 and 7, 3 and 8, and so on).

Throughout the nineteenth century (and part of twentieth century), systematic wide-field observations at a small number of observatories produced a family of large *Durchmusterung* (survey) catalogs, recording the positions and visual magnitudes of hundreds of thousands of stars with a precision sufficient for unambiguous identification. Although not designed as precision astrometric tools, these catalogs occupied a unique role in astronomy for over a century: they were comprehensive enough to serve as a standard reference for stellar designations. For instance, they were used to cover faint comparison stars required for variable-star estimation, as they were relatively homogeneous in their magnitude determinations [1]. Even today, many stars are primarily identified by their Durchmusterung designation in databases such as SIMBAD.

Because these catalogs maintained the greatest sky coverage and depth well into the twentieth century, an international consortium of data centers — including the Centre de Données Astronomiques de Strasbourg (CDS) and the NASA Astronomical Data Center (ADC) — undertook their digitization, see e.g. [2]. For the combined *Bonner Durchmusterung*, *Südliche Durchmusterung*, and *Córdoba Durchmusterung*, this process lasted more than a decade, spanning the 1970s to the late 1980s. Given the volume of data — well over one million individual stellar records — typo errors were unavoidable.

Correcting these errors helps to preserve the historical and scientific legacy embedded in these catalogs: the data represent genuine observational records of the nineteenth-century sky that remain irreplaceable for long-baseline studies (e.g. variable-star behavior, proper-motion anomalies, transient events).

This paper extends the preliminary results reported in [3], focusing now on the complete first volume of the *Córdoba Durchmusterung* [4] and on additional catalogs processed since then.

As a byproduct of the comparison process, the tools described here also produce machine-readable cross-identification tables. These tables, summarized in Section 4.7, are provided alongside the corrected catalogs in the public repository and may serve as an independent resource for historical research. Notably, one of these catalogs — Gould's Zone Catalog (ZC), a preliminary record of Córdoba zone observations that was never digitized — is made accessible for the first time in machine-readable form through the cross-identifications inherited from available old catalogs (Section 4.8).

---

## 2. Catalogs

### 2.1 The Durchmusterung Catalogs

The Durchmusterungs form a complementary trio of surveys covering the entire celestial sphere.

The *Bonner Durchmusterung* (BD) was initiated under the direction of F. W. A. Argelander at the Bonn Observatory and catalogued stars down to approximately visual magnitude 9.5 from the north pole to declination −1° [5]. Schönfeld later extended it southward to declination −23° as the *Südliche Durchmusterung* (SD) [6]. Finally, the *Córdoba Durchmusterung* was undertaken at the Observatorio Nacional Argentino (ONA) under the successive directorships of J. M. Thome and other astronomers, extending coverage from declination −22° to the south pole. Publication spanned from 1892 to 1932, and the catalog records about 614,000 stellar entries, with positions given for epoch B1875.0 and magnitudes given in the visual scale of that era. For more historical context, the reader may consult Chapter 15 of [7].

All three catalogs sacrificed astrometric precision in exchange for sky coverage, but their positional uncertainty was sufficient to identify stars without ambiguity. The precision in the first volume of CD is approximately 22.5 arcseconds (1σ) as reported by Thome himself in the proem of that volume [4], consistent with the RMSE of 24 arcseconds derived from comparison with PPM in this work. These catalogs are complemented with the *Cape Photographic Durchmusterung* (CPD) that also covers the southern sky but its magnitudes are photographic, not visual, see Chapter 6 of [8].

| Catalog | Dec. range | Lim. mag | Stars | Epoch | Precision (arcsec) |
|---------|-----------|----------|-------|-------|--------------------|
| Bonner Durchmusterung (BD) | −1° to +89° | ~9.5 | ~325,000 | 1855.0 | ~36 |
| Südliche Durchmusterung (SD) | −23° to −2° | ~9.5 | ~135,000 | 1855.0 |  |
| Córdoba Durchmusterung (CD) | −89° to −22° | ~10 | ~614,000 | 1875.0 | ~24 |
| Cape Photographic Durchmusterung (CPD) | −89° to −18° | ~12.5 | ~455,000 | 1875.0 |  |

*Table 1. Summary of the large nineteenth-century survey catalogs. Precision column gives the RMSE of angular separation relative to PPM (for the first volumes of BD and CD). The CPD limiting magnitude is photographic.*

For BD, we consider as authority the last edition of the printed catalog [5] which was the source base for the digital BD [9]. For CD, we consider the original printed catalog [4] but with those entries overridden by subsequent corrigenda published by the same Observatory: possibly another volume of the same series, such as Vol. 18, proemio pp. 54-56, or an individual paper [10].

This work is mainly focused on finding typo errors in the first volume of CD [11] (around 180,000 stars) and, to a lesser extent, the first volume of BD (around 111,000 stars).

### 2.2 Reference and Comparison Catalogs

The primary reference catalog used for detecting errors in CD is the *Positions and Proper Motions* catalog (PPM), a compilation of positions at J2000.0 for ~469,000 stars [12]. Of these, 172,000 carry a cross-identification with a CD star (about 28% of CD). The PPM positions are sufficiently precise — the star with the largest positional uncertainty (i.e. PPM 266698) among those identified with CD reports $\sigma_{\alpha *}$ = 0.57 arcseconds and $\sigma_{\delta}$ = 0.61 arcseconds when propagated to B1875.0 from their former mean epoch of 1973 — that any residual error from proper-motion propagation is negligible compared to the intrinsic uncertainty of the Durchmusterungs.

Additional comparison catalogs used include the SD and CPD as we will see in the next section, and the *Córdoba* series of *Astronomische Gesellschaft Katalog* (AGK) zones, labelled Cordoba A, B, and C (epoch B1900.0, declinations −22° to −37°) [13]. These catalogs were measured at roughly the same epoch as CD and thus require no proper-motion correction. More about AGK can be consulted in Chapter 17 of [7].

### 2.3 Cross-Identification Catalogs

A crucial component of the method is that the cross-identifications between the Durchmusterungs and comparison catalogs must be *independent* of the digital positions: they must have been derived by humans comparing the printed catalogs, not by comparing the digital files. Besides the cross-identification information present in PPM catalog, two published cross-identification tables meet this requirement and are used here:

- **Catalog IV/5** (*Table of Correspondences BD/CD/CPD*, Jung & Bischoff 1971) [14]: provides correspondences between BD, CD, and CPD for irregular sky patches.

- **Catalog IV/11** (*Correspondences CD/CPD, Zones −18 to −39*, Bonnet) [15]: a more complete table covering CD and CPD for declinations −18° to −39°.

Note that these catalogs, being hand-made, are not exempt from errors. Moreover, some of these errors appear to have survived and are now found in modern databases. For instance, in catalog IV/5, SD -22°5630 is identified
to CD -22°15216 and this identification is also
observed in [SIMBAD](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=HD+202025&submit=submit+id). However, the angular 
distance between both stars is 4356 arcsec. In fact,
SD -22°5630 is the same as HD 202025 which is near CD -22°15266. On the other hand CD
-22°15216 is near TYC 6358-218-1, a different star.

### 2.4 Other Old Catalogs

Besides [13], there are plenty of other old catalogs available at GAVO Data Center. In this work, we tried to curate some of them by comparing their entries with other catalogs via cross-identifications available on the catalog itself. Below, we enumerate them and their peculiarities (except UA, digital catalogs can be downloaded from GAVO Data Center):

- *Uranometría Argentina* (UA): This is a catalog of ~8,500 bright stars down to 7th magnitude for epoch 1875.0 done at the Observatorio Nacional Argentino under the management of B. A. Gould. It covers all the southern sky and 10 degrees of the northern one. It contains Bayer and Flamsteed designations and references to BD, Lacaille, Lalande, Taylor, Oeltzen-Argelander South, Yarnall (2nd Ed.) and Weisse among others. More info in Chapter 5 of [7]. The digital version is available with code V/135A in Vizier.

- *Argentine General Catalog* / *Catálogo General Argentino* (GC): A catalog of ~32,000 stars down to approximately magnitude 9 for epoch 1875.0, published in 1886 as part of the *Resultados del Observatorio Nacional Argentino* series under the direction of B. A. Gould. It fairly covers the southern sky. More info in Chapter 6 of [7]. Although the printed version contains cross-identifications to other catalogs, unfortunately this information was not added in the digital catalog.

- *Histoire Céleste Française* (Lal): A catalog of ~47,000 stars down to approximately magnitude 9, with positions reduced to epoch 1800.0. Observations were made at the Paris Observatory by J. J. de Lalande during the late 18th century; the standard reduced edition was published by F. Baily in 1847. It covers roughly declinations −30° to the north pole. See [16].

- *Cape Catalogue of 12,441 Stars* (St): A catalog of ~12,000 stars for epoch 1880.0, observed at the Royal Observatory, Cape of Good Hope, by E. J. Stone, and published in 1881. It covers the sky accessible from the Cape of Good Hope, from approximately declination +50° to the south pole, although it is focused on the southern sky. In this work, it is also used when Lacaille (L) designations are found, since this catalog fairly covers Lacaille's stars and has better precision. Digital version contains cross-identifications to older catalogs although this information was not exploited in this work.

- *Katalog der Argelander'schen Zonen* (W): A catalog of ~18,000 stars for epoch 1850.0 covering declinations −15° to −31°, based on F. W. A. Argelander's zone observations reduced by E. Weiss at the Vienna Observatory. In this work, it is also used when Oeltzen-Argelander South (OA) designations are found, since this catalog covers most of OA's stars and the latter is not available in digital form. See [16] for more info about OA catalog.

- *Catalogue of 16,748 Southern Stars* (G): A catalog of 16,748 stars for epoch 1850.0 derived from zone observations made by the U.S. Naval Astronomical Expedition at a temporary observatory on Cerro Santa Lucía, Santiago de Chile (1849–1852), under the direction of J. M. Gilliss. The reductions were completed at the U.S. Naval Observatory and published in 1895. It covers about 30 degrees around the south pole. It has references to Lacaille, Stone and GC catalogs among others. It also has references to Gould's Zone Catalog. 

- *Catalogue of Stars observed at the United States Naval Observatory* (U): A catalog of ~11,000 stars for epoch 1860.0 observed at the U.S. Naval Observatory in Washington during 1845–1877 [17]. The 3rd edition was originally compiled by M. Yarnall and later revised by E. Frisby (1889). It irregularly covers declinations from approximately −35° to the north pole. It contains designations to several other catalogs of that era: BD, BAC, Lacaille, Lalande, Taylor, Oeltzen-Argelander North and South, and WB among others. Cordoba catalogs refer to the 2nd edition of this catalog (1845-1871) [16], so number designations may differ up to ~300 units. For instance, star 37G Pyx of UA (HD 75722) corresponds to Y. 3739 (2nd Ed.) which is U. 3817 (3rd Ed.).

- *Zonen-Beobachtungen, nördliche Zonen* (OARN): A catalog of ~26,000 stars for epoch 1842.0 covering declinations +45° to +80°, based on Argelander's northern zone observations and reduced by W. Oeltzen. It is the northern counterpart of the OA South zones used in this work.

- *General Catalogue of the Principal Fixed Stars* (T): A catalog of ~11,000 stars down to approximately magnitude 8 for epoch 1835.0, observed at the Madras Observatory (India) during 1830–1843 under the direction of T. G. Taylor [16]. The digital version contains cross-identifications with Lacaille and GC catalogs, among others. It covers the whole sky.

Other catalogs considered are: *Segundo Catálogo General Argentino* (G2), *Catalogue of Stars of the British Association for the Advancement of Science* (BAC), *Positiones Mediae Stellarum Fixarum* (WB) and some Thome's zone observations in 1881-1884.

---

## 3. Methods

### 3.1 Main Method: Cross-Identification and Threshold Comparison

The core detection strategy exploits the fact that a typo error in a digital catalog manifests as an anomalously large discrepancy between the digitized position (or magnitude) and an independently obtained counterpart. The algorithm operates as follows:

1. For each star in the reference catalog (e.g., PPM) that carries a cross-identification to CD, precess its J2000.0 coordinates and apply proper motion to obtain its expected position at epoch B1875.0 (in the FK4 system, using the WCSTools `wcsconp` routine) [18].
2. Compute the angular distance between the propagated reference position and the position recorded in the digital CD.
3. If this distance exceeds a pre-set threshold $d_{pos}$, log the entry for manual review, together with the entry as it would appear in the printed catalog (including an estimate of the printed page).
4. If the reference catalog provides a visual or Johnson *V* magnitude, compare it with the CD magnitude after applying a polynomial magnitude-scale transformation (see below). If the magnitude difference exceeds a threshold $d_{mag}$, log the entry separately.
5. Logged entries are then examined manually: the digital record is compared against the corresponding entry in the printed catalog, and any genuine digit error is corrected. Also check if the star is mentioned in any corrigenda (in that case, it requires further investigation).

The choice of thresholds represents a trade-off between sensitivity and the number of false positives that must be manually inspected. For CD versus PPM, thresholds of $d_{pos}$ = 120 arcseconds and $d_{mag}$ = 1.5 magnitudes (approximately 5σ in both cases) were found to produce a manageable log while retaining a high success rate.

**Magnitude scale transformation.** The visual magnitudes of the Durchmusterungs and modern Johnson *V* magnitudes follow slightly different scales. For each volume of CD, a quadratic polynomial transformation was fitted by least squares using the cross-identified PPM stars. For the first volume of CD (declinations −22° to −31°), the transformation derived from 11,659 stars is:

$mag_V = −0.157169 + 1.188316 · V − 0.022130 · V^2$

with RMSE = 0.30 magnitudes. Separate fits were performed for each of the five volumes of CD, with RMSEs ranging from 0.28 to 0.35 magnitudes.

**Limitations.** Only errors large enough to exceed the position or magnitude threshold are detected. A typo that shifts a position by less than $d_{pos}$, or a magnitude error smaller than $d_{mag}$, will be missed. Furthermore, the method is limited to the fraction of stars that carry a pre-existing cross-identification: only 28% of CD stars are identified in PPM. Approximately 72% of CD thus remains unexplored by this technique alone.

### 3.2 Alternative Method: Reverse-Engineering of Printed Precession Constants

Certain old catalogs (notably the *Argentine General Catalog*) include, for each star, a pre-computed annual precession constant for right ascension and declination. These constants were computed from the star's own position using the Bessel formulae with Struve/Peters coefficients, see pp. 23-24 of [17]. For instance, for epoch 1875.0 used by GC catalog we have:

$p_{\alpha} = 3.072245 + 1.33695 · sin(\alpha) · tan(\delta) ~~~[s]$

$p_{\delta} = 20.05425 · cos(\alpha) ~~~['']$

where $\alpha$ and $\delta$ are the right ascension and declination of the star. A typo error in the position will therefore cause the reported precession constant to disagree with the value recomputed from the (incorrect) position. Conversely, a typo error in the *precession column itself* — leaving the position unchanged — can be detected directly by comparing the reported constant against the value recomputed from the reported position.

This technique was implemented in GC, *Yarnall-Frisby USNO 3rd ed.* and *Gilliss* catalogs, yielding some candidates in which the discrepancy was unambiguously in the precession column of the printed catalog rather than in the position. For instance, in the printed catalog, GC 1795 reports an annual precession in right ascension of 0.409 s/yr, whereas the value recomputed from its position is $p_{\alpha}$ = 1.409 s/yr — a clear digit dropout. Since the annual precession in declination is fine, it means there is no error in right ascension (because $p_{\delta}$ only depends on $\alpha$) leading to two hypotheses: the aforementioned error is in $p_{\alpha}$, or in $\delta$. Finally, we can discard an error in declination since the distance between GC 1795 and the nearest PPM is less than an arc second (PPM 352586).

Detections with this technique were rather limited.

### 3.3 Isolated-Star Detection

For catalogs without pre-existing cross-identifications, a complementary approach flags stars in the old catalog for which no counterpart is found within a generous angular radius in any of several reference catalogs (PPM, CD, CPD, or the Guide Star Catalog). A star that is truly isolated in the old catalog — with no plausible counterpart nearby — is a candidate for a position typo so large that the star appears displaced far from its true location. This method was applied to GC and yielded additional candidates, although several entries correspond to nebular objects and thus are false positives (e.g., GC 18278 is the star cluster omega Centauri).

---

## 4. Results

### 4.1 Algorithm Efficiency

Before applying the method to detect genuine errors, its sensitivity was evaluated by comparing the current VizieR version (I/114) against two older digital versions of CD with known differences:

**I/114 vs. I/88.** When comparing I/114 against I/88, an earlier digitized version (circa 1982), 205 differences were found: 202 in position and 3 in magnitude, plus 44 entries with unrecognized magnitude codes in I/88. Running both independently against PPM, the older version generated 51 additional position warnings and 3 additional magnitude warnings. The algorithm thus captured roughly 25% of the known positional differences (51 out of 202), but recovered all 3 magnitude differences.

**I/114 vs. 1114.** When comparing against an older NASA-ADC version (code 1114), 141 positional differences were found between them. Running both independently against PPM, the older version generated 35 additional position warnings, meaning the algorithm captured roughly 25% of the known differences (35 out of 141).

In both experiments the positional recovery fraction is roughly 25%, reflecting a fundamental limitation of the approach: the lack of coverage of CD with PPM stars. It does, however, suffice to validate the chosen $d_{pos}$ and $d_{mag}$ thresholds.

### 4.2 Results of the CD×PPM comparison

The comparison of PPM against the first volume of CD (179,804 stars, declinations −22° to −31°) yielded 40,769 cross-identified pairs (23%). Running the detection algorithm with thresholds of 120 arcseconds and 1.5 magnitudes produced:

| Category | Flagged | False positives | True corrections |
|----------|---------|-----------------|-----------------|
| Position errors | 72 | 9 | 63 |
| Magnitude errors | 12 | 5 | 7 |
| **Total** | **84** | **14** | **70** |

*Table 2. Results of the CD×PPM comparison.*

The success rate of 83% (70 genuine corrections out of 84 flagged) compares favorably with the BD comparison (see below). Each flagged entry was examined against the printed catalog [4] to confirm that the digital-to-printed difference is consistent with a typical typographic mistake (e.g., a single digit changed or swapped between consecutive rows).

An independent consistency check is provided by the RMSE of angular separation for the non-flagged pairs: 24.0 arcseconds, in excellent agreement with Thome's own estimate of 22.5 arcseconds mentioned in the first section.

### 4.3 Typo Errors Found in CD

The corrected entries include 78 corrections to right ascension, 40 corrections to declination, and 9 corrections to magnitude (some stars required corrections to more than one parameter). Among the most diagnostic corrections are cases where two consecutive entries in the digital file have their parameters swapped (a signature of a row-level transposition).

It was noticed that the angular discrepancy between the digital position and the PPM counterpart is explained almost exactly by the positional difference between the correct and incorrect printed values, providing strong evidence that the cause is typographic rather than observational. For instance, digital catalog reports 16'4 in declination of CD -24°11535 but it should be 46'4 according to [4]; the difference reported by the angular distance between that star and PPM 263064 is 1772 arcseconds, very near to the 1800 arcsecond difference between 16'4 and 46'4.

The corrected first-volume catalog, as well as the list of differences and corrections is available at the author's repository.

### 4.4 CD Declination −22° via the Southern Durchmusterung

For the declination −22° strip of CD, PPM does not provide identifications because the SD boundary lies at this declination. Instead, the SD was used as the reference, with identifications obtained through catalog IV/5 [14]. Of 3,485 identified pairs (21% of the −22° strip), the algorithm flagged 21 position and 12 magnitude discrepancies. Manual inspection yielded only 4 genuine typo corrections (success rate 12%), notably CD −22°9451 and CD −22°9452, whose values in the digital file appear to have been swapped between the two consecutive entries.

### 4.5 Typo Errors Found in the Bonner Durchmusterung

The comparison of PPM with the first volume of BD (declinations −1° to +19°, 63,034 cross-identified pairs) used more conservative thresholds of 180 arcseconds and 0.4 magnitudes. This yielded 86 flagged entries (82 position, 4 magnitude). Manual inspection found **13 genuine position errors** (15% success rate); reducing the position threshold to 150 arcseconds added 61 more flagged entries but found only 1 additional correction, confirming that the 180-arcsecond threshold already captures most detectable errors at a reasonable false-positive rate.

The curated catalog is also available at the author's repository.

### 4.6 Other Catalogs

Additional tools made by the author extended the methodology to a broader set of nineteenth-century catalogs by performing cross-identifications through angular proximity (in the absence of pre-existing cross-identification tables) against PPM, CD, CPD, and the Guide Star Catalog:

| Catalog | Epoch | Corrections |
|---------|-------|-------------|
| Uranometría Argentina (UA) | 1875 | 48 |
| Catálogo General Argentino (GC) | 1875 | 3 |
| Resultados XV, ONA | 1881–1884 | 7 |
| Yarnall-Frisby USNO (U) | 1860 | 13 |
| Gilliss (G) | 1850 | 22 |
| Weiss (W) | 1850 | 3 |

*Table 3. Summary of corrections applied to other nineteenth-century catalogs.*

Additional experiments comparing CD against CPD (via IV/5 and IV/11) and against the AGK zones were also run, but the resulting logs have not been fully reviewed and are provided as supplementary material in the repository.

### 4.7 Cross-Identification Tables Generated

A byproduct of the comparison workflow is a set of machine-readable cross-identification tables stored as CSV files in the author's repository. For each historical catalog processed, the tool records which stars were successfully matched against PPM, CD, CPD, and — transitively through PPM, which already carries them — SAO and HD identifiers. These pairings were not previously available in digital form for most of the catalogs listed.

Table 4 summarizes the number of cross-identification pairs generated for each source catalog. A dash indicates that no cross-identification file was generated for that combination, typically because the source catalog covers declinations not well matched by CD or CPD (e.g., Lalande and Taylor span largely northern skies).

| Source catalog | ×PPM | ×CD | ×CPD | ×SAO | ×HD |
|----------------|-----:|----:|-----:|-----:|----:|
| CD (all 5 vols.) | 152,177 | — | — | — | — |
| GC | 31,706 | 28,436 | 28,910 | 30,025 | 29,567 |
| GC2 | 5,357 | 5,451 | 5,252 | 4,829 | 4,099 |
| Gilliss | 15,773 | 15,941 | 16,176 | 5,663 | 11,251 |
| Lacaille | 9,615 | 9,606 | 9,599 | 9,603 | 9,501 |
| Lalande | 46,914 | — | — | — | — |
| OA (South) | 18,011 | 9,296 | 13,146 | 17,199 | 15,543 |
| Taylor | 10,944 | — | — | — | — |
| UA | 7,688 | 4,523 | 4,758 | 7,654 | 7,621 |
| Yarnall-Frisby USNO | 10,676 | 4,295 | 4,658 | 10,224 | 9,702 |
| ZC (see next section) | 5,552 | 5,755 | 5,756 | 2,380 | 4,709 |

*Table 4. Number of cross-identification pairs generated per source catalog and target catalog. SAO and HD identifications are obtained transitively through PPM.*

The CD×PPM files (one per volume) are additionally used in the error-detection workflow described in Section 3.1; the remaining files constitute stand-alone cross-identification resources. All files share the same CSV format: `index1, index2, mag, dist`, where `mag` is the source catalog magnitude and `dist` is the angular separation in arcseconds at the comparison epoch.

### 4.8 Cross-Identifications for Gould's Zone Catalog

Gould's Zone Catalog (ZC), compiled by B. A. Gould at the Observatorio Nacional Argentino, represents a record of zone observations made during the same program that ultimately produced the *Córdoba Durchmusterung*. Unlike the CD and the other survey catalogs discussed above, the ZC was never digitized and is therefore not available in any online data center.

However, the digital version of some of Thome's zone observations (*Resultados del Observatorio Nacional Argentino*, Tomo XV, 1881–1884) and the Gilliss catalog include, for many of their entries, a cross-reference to a ZC star. When these entries are in turn cross-identified against PPM, CD, and CPD by the comparison tool, the ZC identifications are automatically inherited, yielding cross-identification lists that connect ZC designations to modern catalog numbers.

In this way, 5,796 (out of 73,160) unique ZC stars were registered, representing the 8% of the printed catalog.

---

## 5. Conclusions

We have described and applied a simple but effective method for detecting and correcting typographical errors in the digital versions of nineteenth-century star catalogs. The key insight is that cross-identifications between old catalogs and a modern reference, when independently derived from the printed sources, can expose errors introduced during digitization without reference to the digital file itself.

The method has important limitations: it can only detect errors that produce discrepancies exceeding the threshold, and it is restricted to the fraction of stars carrying pre-existing cross-identifications. Future work could focus on (i) correcting other volumes of the Durchmusterungs with the available information in digital form (catalogs SAO, PPM, AGK, IV/5 and IV/11), (ii) correcting other old available catalogs, (iii) scanning old catalogs via OCR/AI and comparing them against their digital versions, and (iv) researching systematic errors similar to how it was done with 18th-century catalogs [19].

The curated catalogs and CSV tables are publicly available at https://github.com/aureus123/catalogs, together with the source code of all comparison tools and the full logs of flagged entries.

---

## References

[1] Williams, T. R.; Saladyga, M. (2011). *Advancing Variable Star Astronomy: The Centennial History of the American Association of Variable Star Observers*. Cambridge: Cambridge University Press. ISBN 9780521519120.

[2] Wagner, M. J. (1984). The BD : a progress report. *Bulletin d'Information du CDS*, Vol. 26, p. 87. https://articles.adsabs.harvard.edu/pdf/1984BICDS..26...87W

[3] Severin, D. E.; Sevilla, D. J. (2015). Development of a new digital version of "Cordoba Durchmusterung" stellar catalog. *Revista Académica Electrónica de la UNR*, Vol. 15 (8), pp. 2250–2260. https://rephip.unr.edu.ar/server/api/core/bitstreams/0362c1df-b472-4216-99fe-46c7f135921a/content

[4] Thome, J. M. (1892). Zonas de exploración (Córdoba Durchmusterung): declinación −22° a −32°. *Resultados del Observatorio Nacional Argentino*, Vol. 16. Buenos Aires: Pablo E. Coni e Hijos.

[5] Schmidt, H. (1968). *Bonner Durchmusterung, Nordlicher Teil, Deklinations-Zonen -1° bis +89° Grad Sternverzeichnis, vierte Auflage*. Bonn: Ferd. Duemmlers Verlag.

[6] Schmidt, H. (1967). *Bonner Durchmusterung, Südlicher Teil, Deklinations-Zonen -2° bis -22° Sternverzeichnis, dritte Auflage*. Bonn: Ferd. Dümmlers Verlag.

[7] Minniti, E.; Paolantonio, S. (2024). *Córdoba Estelar: Desde los sueños a la Astrofísica. Historia del Observatorio Nacional Argentino (1871–1942)*. Universidad Nacional de Córdoba. http://www.cordobaestelar.oac.uncor.edu/Cordoba_Estelar_2024.pdf

[8] van der Kruit, P. C. (2014). *Jacobus Cornelius Kapteyn: Born Investigator of the Heavens*. Astrophysics and Space Science Library, Vol. 416. Cham: Springer. ISBN 9783319108759.

[9] Vizier catalog I/122: Bonner Durchmusterung at CDS. https://cdsarc.cds.unistra.fr/viz-bin/cat/I/122

[10] Thome, J. M. (1896). List of additional errors found in the Cordoba catalogues. *The Astronomical Journal*, Vol. 364, p. 30.

[11] Vizier catalog I/114: Córdoba Durchmusterung at CDS. https://cdsarc.cds.unistra.fr/viz-bin/cat/I/114

[12] Röser, S.; Bastian, U. (1991). PPM Star Catalogue. Vizier catalogs I/146, I/193, I/206, I/208, e.g. https://cdsarc.cds.unistra.fr/viz-bin/cat/I/146

[13] Astronomische Gesellschaft Katalog: Cordoba A, B, C at GAVO Data Center in Heidelberg, codes T0986P10, T0987P10, T0988P10. https://dc.g-vo.org/arigfh/katkat/q/form

[14] Jung, J.; Bischoff, M. (1971). *Table of Correspondences BD/CD/CPD*. NASA Astronomical Data Center, catalog IV/5. [Also available at VizieR.](https://cdsarc.cds.unistra.fr/viz-bin/cat/IV/5)

[15] Bonnet, R. M. (circa 1970). *Correspondences CD/CPD, Zones −18 to −39*. NASA Astronomical Data Center, catalog IV/11. [Also available at VizieR.](https://cdsarc.cds.unistra.fr/viz-bin/cat/IV/11)

[16] Paolantonio, S. (2014). Libros con historia III. Available at https://historiadelaastronomia.wordpress.com/documentos/librosconhistoriaiii/

[17] Yarnall, M.; Frisby, E. (1889). *Catalogue of Stars observed at the United States Naval Observatory during the Years 1845 to 1877*, Third Edition. Washington: U.S. Government Printing Office.

[18] Mink, D. J. (1997). WCSTools: Putting Image World Coordinate Systems to Use. *Astronomical Data Analysis Software and Systems VI, A.S.P. Conference Series*, Vol. 125, pp. 249-252. San Francisco: Astronomical Society of the Pacific.

[19] Lequeux, J. (2014). From Flamsteed to Piazzi and Lalande: new standards in 18th century astrometry. *Astronomy & Astrophysics*, Vol. 567, A26. https://doi.org/10.1051/0004-6361/201423899
