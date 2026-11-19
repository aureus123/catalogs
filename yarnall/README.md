# Curiosities of Yarnall-Frisby USNO Catalog

## Cross-Identification File

`yarnall_usno.txt` contains cross-identifications between the **[2nd Edition](https://archive.org/details/catalogueofstars00sand)** of the Yarnall-Frisby catalog (entries prefixed with `Y`) and the **[3rd Edition](https://archive.org/details/cataloguestarsus00unitrich)** (entries prefixed with `U`) obtained from cross-identifications between Uranometría Argentina, Resultados XV and GC with USNO catalog. Each line gives the correspondence between a star's number in the two editions, e.g. `Y 91 = U 92`.

## Zone Files

The files `transit_zone.txt`, `mural_zone.txt`, `meridian_zone.txt`, and `anonymous.txt` list USNO stars that have **no prior identification** in earlier catalogs such as Lalande or Lacaille. These stars can therefore be considered as having been **observed for the first time by the US Naval Observatory**. Each line contains three comma-separated fields:

```
<catalog_id>, <visual_magnitude>, <constellation>
```

The constellation was derived by converting the original B1860.0 coordinates to a modern epoch and using the modern IAU abbreviations.

The three instrument-based files indicate which instrument was used to measure the stars:

| File | Instrument | Stars |
|------|-----------|-------|
| `transit_zone.txt` | Transit Instrument | 149 |
| `mural_zone.txt` | Mural Circle | 393 |
| `meridian_zone.txt` | Meridian Circle | 38 |
| `anonymous.txt` |  | 487 |
| **Total** | | **1,067** |

All 1,067 stars across these four files represent the USNO's original observations, with no cross-reference to any prior catalog.

## Gilliss Expedition

`anon_gil1963.txt` lists, in the same format, the 385 stars of J. M. Gilliss's *Catalogue of 1963 stars* (`cat/gil1963.txt`, epoch 1850.0, observed at Santiago de Chile in 1850-52) that carry neither a Bayer/Flamsteed name nor a cross-reference to an earlier catalogue (Lacaille, B.A.C., Brisbane, etc.). These stars could therefore have been **observed for the first time by the U.S. Naval Astronomical Expedition** led by Gilliss. The constellation was derived from the B1850.0 coordinates, and five stars without a magnitude in the catalogue have an empty second field.

The list comes from the curated digital file, not the printed catalogue, so it doesn't match the printed "Anonymous" entries exactly. Some stars printed as Anonymous carry a Lacaille number in the digital file (e.g. G 3, 80, 121, 140) and are left out here. Some stars whose printed name is a doubtful identification ("1267 B.A.C.?", G 105) or a lost one (G 185 = 2063 Lacaille) appear here. Also, seven stars (G 52, 53, 395, 891, 1156, 1261, 1658) are missing from the digital file.
