#!/usr/bin/env python3
"""Common interface over the three photometric reference catalogues.

    from starcat import get_catalog, blend_dmag, e_V
    cat = get_catalog("gaia_online", epoch=2026.3)
    idx, sep, nn2 = cat.match(ra, dec)

Backends, all exposing ``.data`` (a structured array with at least ``ra, dec, V, BV``),
``box()``, ``cone()``, ``match()`` and ``designation(i)``:

===============  ===========  ==========================================================
name             network      notes
===============  ===========  ==========================================================
``tycho2``       no           cat/tyc2_phot.npy. V <= 11.5 useful, sigma(V) 0.015-0.09
                              growing with magnitude. B-V transform valid only B-V < 1.53.
``gaia_v50``     no           ASTAP local database. V quantised to 0.1 mag (0.029 mag
                              floor), B-V to 0.02 mag but capped at +-2.54.
``gaia_online``  **YES**      VizieR Gaia DR3, full precision, no colour cap, depth on
                              demand. Cached per field, but a new field needs network.
===============  ===========  ==========================================================

Measured head-to-head on RCen (271 detections, 143-star common subset, identical blend
mask): robust residual sigma 0.068 (tycho2), 0.047 (gaia_v50), 0.031 (gaia_online). The
v50/online gap is exactly v50's quantisation floor.
"""
import numpy as np
from scipy.spatial import cKDTree

CATALOG_NAMES = ("tycho2", "gaia_v50", "gaia_online")
NEEDS_NETWORK = {"gaia_online"}

#: Which photometric system each backend's V and B-V live on.
#:
#: This governs whether a coefficient set fitted against one catalogue may be *applied*
#: using another. Coefficients are only transferable within a system:
#:
#: * ``gaia_v50`` and ``gaia_online`` are the same construction -- V50 stores exactly the
#:   transform ``gaia_online`` computes on the fly -- so they differ only in precision
#:   and are freely interchangeable.
#: * ``tycho2`` is a different system. Measured over 63 Dwarf 3 fields, the same data
#:   fitted against Tycho-2 vs Gaia gives Tbv 1.671 vs 1.527 and Tv_bv -0.039 vs -0.013.
#:   Those differences *are* the system mismatch: applying one system's colour terms
#:   while fitting zero points against the other biases B-V by (dTbv * colour). The
#:   ~0.04 mag zero-point offset between the systems is harmless by comparison, since
#:   the per-image ``zp`` absorbs it.
PHOTOMETRIC_SYSTEM = {
    "tycho2": "tycho2",
    "gaia_v50": "gaia_dr3",
    "gaia_online": "gaia_dr3",
}


#: Prefix each backend's ``designation()`` puts on the CSV ``ref_id`` column.
#:
#: This is the pipeline's provenance marker: measure.py's CSV is self-describing, so
#: nothing downstream needs a sidecar file to know which catalogue produced it. Keep
#: these in sync with the ``designation()`` methods -- ``test_ref_id_prefixes`` checks it.
#:
#: Note the Gaia forms are **synthesised from coordinates and are not real catalogue
#: identifiers** -- V50 stores no identifiers at all, and the VizieR query does not
#: currently request Gaia's ``Source`` column. They are unique and stable positionally,
#: but not resolvable as designations.
ID_PREFIX = {
    "tycho2": "TYC ",
    "gaia_v50": "V50 J",
    "gaia_online": "GaiaDR3 J",
}


def catalog_from_ref_id(ref_id):
    """Which catalogue produced a ``ref_id``, or None if it is blank/unrecognised."""
    if not ref_id:
        return None
    for name, prefix in ID_PREFIX.items():
        if ref_id.startswith(prefix):
            return name
    return None


def catalog_from_rows(rows, column="ref_id"):
    """Infer the fit catalogue from a measure.py CSV by majority vote over ref_id.

    Rows with no catalogue match are blank, so this ignores them; a vote rather than
    first-hit keeps it robust to a stray malformed row.
    """
    from collections import Counter
    votes = Counter(c for c in (catalog_from_ref_id((r.get(column) or "").strip())
                                for r in rows) if c)
    return votes.most_common(1)[0][0] if votes else None


def same_system(a, b):
    """True when coefficients fitted against catalogue ``a`` may be applied while
    fitting zero points against catalogue ``b``."""
    return PHOTOMETRIC_SYSTEM.get(a) == PHOTOMETRIC_SYSTEM.get(b)


#: Which precalibration to use for a given fit catalogue, most preferred first.
#:
#: The caller picks only the *fit* catalogue; the precalibration follows from it, so the
#: two can never end up on different photometric systems by accident. Within the Gaia
#: system the online precalibration is preferred whenever it exists because it is not
#: quantised (V50 stores V in 0.1 mag steps, a 0.029 mag floor), and V50 is the offline
#: fallback. The coefficients are interchangeable either way -- see PHOTOMETRIC_SYSTEM.
PRECAL_PREFERENCE = {
    "tycho2": ("tycho2",),
    "gaia_v50": ("gaia_online", "gaia_v50"),
    "gaia_online": ("gaia_online", "gaia_v50"),
}

#: A neighbour this many magnitudes fainter than the target contributes
#: 10**(-0.4*3) = 6.3% of its flux. Anything fainter is treated as harmless.
DEFAULT_BLEND_DMAG = 3.0
DEFAULT_BLEND_RADIUS_ARCSEC = 20.0


def get_catalog(name, **kwargs):
    """Instantiate a backend by name. Unknown kwargs are passed through, so callers can
    give e.g. ``epoch=`` (online) or ``mag_limit=`` (v50/online) without branching."""
    name = name.lower()
    if name == "tycho2":
        from tycho2 import Tycho2
        return Tycho2(**{k: v for k, v in kwargs.items() if k in ("npy_path",)})
    if name == "gaia_v50":
        from gaia_v50 import GaiaV50
        allowed = ("db_path", "db_name", "mag_limit", "ext")
        return GaiaV50(**{k: v for k, v in kwargs.items() if k in allowed})
    if name == "gaia_online":
        from gaia_online import GaiaOnline
        allowed = ("mag_limit", "epoch", "cache_dir", "verbose")
        return GaiaOnline(**{k: v for k, v in kwargs.items() if k in allowed})
    raise ValueError(f"unknown catalogue {name!r}; expected one of {CATALOG_NAMES}")


def e_V(data):
    """Per-star magnitude uncertainty, whatever the backend calls it.

    Tycho-2 stores a real measured ``e_VT``; the Gaia backends store a constant floor
    (v50 its 0.029 mag quantisation, online the ~0.046 mag transform scatter).
    """
    names = data.dtype.names
    if "e_V" in names:
        return data["e_V"]
    if "e_VT" in names:
        return data["e_VT"]
    return np.full(len(data), np.nan)


def blend_dmag(data, idx, radius_arcsec=DEFAULT_BLEND_RADIUS_ARCSEC):
    """Contamination metric that does **not** get stricter as the catalogue deepens.

    For each matched catalogue star, returns the magnitude difference to the most
    contaminating catalogue neighbour within ``radius_arcsec``:
    ``min(V_neighbour - V_target)``, or ``+inf`` when it is isolated.

    This replaces the old ``nn2 > 20 arcsec`` rule, which counted *any* neighbour
    regardless of brightness. That rule was depth-coupled and perverse: on RCen it
    removed 7 stars using Tycho-2, 24 using V50 at V<=13 and 117 using V50 unlimited, so
    a better catalogue destroyed the calibrator sample. A V=16 star beside a V=9 star
    contributes ~0.1% of its flux and must not disqualify it.

    ``idx`` is the array returned by ``match()``; entries < 0 (unmatched) yield NaN.
    """
    idx = np.asarray(idx)
    out = np.full(len(idx), np.nan)
    hit = idx >= 0
    if not np.any(hit) or len(data) == 0:
        return out

    ra = np.radians(data["ra"]); dec = np.radians(data["dec"])
    cd = np.cos(dec)
    xyz = np.stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)], axis=1)
    tree = cKDTree(xyz)

    # chord length subtending radius_arcsec
    ang = np.radians(radius_arcsec / 3600.0)
    r_chord = 2.0 * np.sin(ang / 2.0)

    targets = idx[hit]
    V = np.asarray(data["V"], dtype=float)
    neigh = tree.query_ball_point(xyz[targets], r=r_chord)

    worst = np.full(len(targets), np.inf)
    for k, (t, nb) in enumerate(zip(targets, neigh)):
        others = [j for j in nb if j != t]
        if others:
            worst[k] = float(np.min(V[others] - V[t]))
    out[hit] = worst
    return out


#: Extra V uncertainty for a saturated star, keyed on the number of clipped pixels.
#:
#: Measured against Gaia over 12 Dwarf 3 VIS/gain-60 fields (4026 stars, zero point fit
#: on unsaturated calibrators only, so saturated stars are out-of-sample):
#:
#: ===========  =====  ============  ============
#: sat. pixels      n  median bias   robust sigma
#: ===========  =====  ============  ============
#: 0             3941        +0.003         0.064
#: 1               18        +0.167         0.114
#: 2-4             20        +0.213         0.179
#: 5-14            17        +0.314         0.317
#: 15-49           20        +0.776         0.286
#: 50+             10        +0.772         0.759
#: ===========  =====  ============  ============
#:
#: Re-derived after the saturation search box was changed from a fixed 5x5 to one
#: scaled by each star's HFD (photometry.photometer). The relation barely moved, which
#: is expected -- it is physical -- but the *counting* did: bloated saturated stars whose
#: clipped core sits outside a 5x5 box are now detected at all.
#:
#: Clipped pixels lose flux, so the magnitude is biased **faint**. The bias is not
#: corrected -- it depends on the PSF and how far into saturation the star is, so a
#: correction would be fragile -- instead the quoted sigma covers it, as
#: ``sqrt(bias^2 + scatter^2)``. That keeps the true value inside roughly 1 sigma while
#: making it obvious the measurement is degraded. Rederive with the study in the project
#: history if the sub-exposure or gain changes; these numbers are for 10 s subs at gain 60.
SATURATION_SIGMA = (
    (1, 0.20),
    (2, 0.28),
    (5, 0.45),
    (15, 0.83),
    (50, 1.08),
)


def saturation_sigma(n_sat_px):
    """Extra V uncertainty (mag) to add in quadrature for a star with this many
    clipped pixels. Zero when unsaturated."""
    n = np.asarray(n_sat_px)
    out = np.zeros(n.shape, dtype=float)
    for threshold, sigma in SATURATION_SIGMA:
        out = np.where(n >= threshold, sigma, out)
    return out if out.shape else float(out)


def isolated(data, idx, radius_arcsec=DEFAULT_BLEND_RADIUS_ARCSEC,
             dmag=DEFAULT_BLEND_DMAG):
    """Boolean 'clean enough to calibrate with': no catalogue neighbour within
    ``radius_arcsec`` brighter than ``V_target + dmag``. Unmatched stars are False."""
    return blend_dmag(data, idx, radius_arcsec) > dmag
