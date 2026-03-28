# Curiosities of Yarnall-Frisby USNO Catalog

## Cross-Identification File

`yarnall_usno.txt` (69 entries) contains cross-identifications between the **[2nd Edition](https://archive.org/details/catalogueofstars00sand)** of the Yarnall-Frisby catalog (entries prefixed with `Y`) and the **[3rd Edition](https://archive.org/details/cataloguestarsus00unitrich)** (entries prefixed with `U`) obtained from cross-identifications between Uranometría Argentina and Resultados XV with USNO catalog. Each line gives the correspondence between a star's number in the two editions, e.g. `Y 91 = U 92`.

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
