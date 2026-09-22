#!/usr/bin/env python3
"""Generate the LaTeX source for the star catalogue.

Primary source is the Bright Star Catalogue 5th ed. (book/ybsc5.txt).
Gould designations come from the Uranometria Argentina (cat/ua.txt), joined
on HD.  Positions are the BSC5 ones: the table prints RA to the second and
Dec to a tenth of an arcminute, well inside what BSC5 already carries.

A hand-picked list of deep-sky objects (book/ngc2000.txt) is merged into the
same RA-ordered table, sharing its columns but filling only those that mean
anything for a nebula or a cluster -- see NGC_LIST and the Dso class.

Usage:
    python build_catalog.py                    # 2 pages, 55 rows/page
    python build_catalog.py --pages 0          # whole catalogue
    python build_catalog.py --notes bsc5       # notes from ybsc5.notes.txt
    python build_catalog.py --dso none         # stars only, no deep-sky rows
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
BOOK = ROOT / "book"
CAT = ROOT / "cat"

CROSS = ROOT / "results" / "cross"

# Old-catalogue cross-identifications shown in the notes column, in the order
# they appear in the note.  Each entry renders as "Label (number): mag".
#
# A catalogue is a list of tiers, and a tier is (suffix, files).  Files inside
# one tier are merged -- Oeltzen-Argelander is split across two files (southern
# / northern zones) that reuse each other's numbers but share no HD, so the
# merge is unambiguous.  Tiers are consulted in order: a later tier only
# supplies stars the earlier ones do not carry, and its suffix marks the number
# in print, so Gilliss falls back to the 1963 reduction as "Gilliss (1878*)".
CROSS_SOURCES = [
    ("Gould",       [("",  ["cross_gc_hd.csv"])]),
    ("Yarnall",     [("",  ["cross_usno_hd.csv"])]),
    ("Gilliss",     [("",  ["cross_gilliss_hd.csv"]),
                     ("*", ["cross_gil1963_hd.csv"])]),
    ("Taylor",      [("",  ["cross_taylor_hd.csv"])]),
    ("Oeltzen-Arg", [("",  ["cross_oa_hd.csv", "cross_oarn_hd.csv"])]),
    ("BAC",         [("",  ["cross_bac_hd.csv"])]),
    ("Stone",       [("",  ["cross_stone_hd.csv"])]),
]

# A magnitude of exactly 0 in a cross table means "variable, or no magnitude
# recorded" -- such an entry is dropped rather than printed as 0.0.
CROSS_NO_MAG = 0.0

# BSC5 VarID values that carry no information beyond "this star varies"; the
# asterisk on the magnitude already says that, so they get no note.  Bare
# numbers are NSV (suspected-variable) numbers, too obscure to print.
VAR_UNINFORMATIVE = {"Var", "Var?"}


def var_note(var_id: str) -> str:
    """Leading "Var: ..." note, or "" when the designation says nothing useful."""
    if not var_id or var_id in VAR_UNINFORMATIVE or var_id.isdigit():
        return ""
    return f"Var: {var_id}"

STRUVE = BOOK / "Struve.csv"

# F.G.W. Struve's three series share a numbering sequence -- STF 43 and
# STFA 43 (Albireo) are different stars -- so the appendices keep their
# classical roman numeral.  The sigma is rewritten to $\Sigma$ at escape time.
STRUVE_PREFIX = {"STF": "Σ ", "STFA": "Σ I ", "STFB": "Σ II "}

# ---- visual doubles ------------------------------------------------------
# A stelledoppie export, in the same schema Struve.csv already uses:
#   name;cst;SAO;coord;wds_name;last;obs;pa;sep;m1;m2;d_mag;orb
# cut at m1 <= 9.0 and sep <= 180" over the whole sky.  The reported pair
# always has m1 <= 5.5, but the component count needs a system's *inner* pairs
# too -- beta Mon is a triple only because its B-C row is there -- hence the
# looser magnitude cut.  Do not add a declination filter to the query: the site
# silently reduces "+52" to 0 and returns the southern sky alone.
DOUBLES = BOOK / "doubles.csv"

# A companion is worth reporting when it can actually be split and seen in a
# small telescope: closer than 3" it is not resolved, past 180" it is not a
# double to the eye any more, and fainter than 9.0 it is not there at all.
DBL_MIN_SEP = 3.0
DBL_MAX_SEP = 180.0
DBL_MAX_MAG = 9.0

# Counting components is a laxer question than reporting one -- "triple" says
# something is there even when it is too faint or too tight to be the note's
# subject -- so it runs to 11.0, roughly the limit of the apertures this book
# is for.  The separation bound stays, which is what keeps Proxima (2.2 deg
# from alpha Cen) from making the sky's best pair a triple.
DBL_COUNT_MAG = 11.0

MULT_WORDS = {2: "doble", 3: "triple", 4: "cuádruple",
              5: "quíntuple", 6: "séxtuple"}

BSC5 = BOOK / "ybsc5.txt"
BSC5_NOTES = BOOK / "ybsc5.notes.txt"
UA = CAT / "ua.txt"
MAPS = BOOK / "maps"      # atlas plates produced by gen_maps.py

# ---- deep-sky objects ----------------------------------------------------
# NGC 2000.0 (Sinnott 1988), the machine-readable NGC/IC.  Positions there are
# equinox 2000.0 despite the "B2000" label in the byte description -- M 31 sits
# at 0h42.7m +41d16', which is the J2000 place, not the B1950 one.
NGC2000 = BOOK / "ngc2000.txt"
NGC2000_NAMES = BOOK / "ngc2000.names.txt"

# The objects to include, chosen by hand rather than by any magnitude cut: a
# deep-sky list is a list of things worth pointing a telescope at, which no
# single number selects.  Kept sorted for reading; duplicates are harmless.
NGC_LIST = [
    104, 121, 224, 253, 288, 362, 1068, 1316, 1399, 1535, 1543, 1566, 1574,
    1647, 1851, 1904, 1960, 2064, 2067, 2070, 2071, 2074, 2079, 2080, 2099,
    2100, 2158, 2168, 2169, 2232, 2287, 2298, 2323, 2362, 2392, 2422, 2425,
    2437, 2438, 2447, 2451, 2477, 2478, 2516, 2547, 2548, 2632, 2645, 2660,
    2682, 2808, 2867, 2899, 3114, 3115, 3201, 3228, 3242, 3293, 3532, 3621,
    3766, 3918, 4103, 4361, 4409, 4486, 4594, 4609, 4755, 4833, 5139, 5156,
    5206, 5236, 5253, 5264, 5272, 5408, 5460, 5904, 6025, 6093, 6121, 6231,
    6337, 6369, 6397, 6405, 6475, 6520, 6523, 6530, 6618, 6656, 6705, 6752,
    6809, 6853, 6934, 7009, 7078, 7089, 7213,
]

# NGC 2000.0 type code -> the word that opens the note.
#
# "Nb" is the one code that does not map cleanly: the catalogue defines it as
# "bright emission *or* reflection nebula" and offers nothing that separates
# the two, so it becomes the neutral "nebulosa difusa" -- true of NGC 2064,
# 2067 and 2071 (reflection, the M78 complex) and of NGC 6523 (emission, the
# Lagoon) alike, and contrasting usefully with "nebulosa planetaria".
DSO_TYPE = {
    "Gx":  "Galaxia",
    "OC":  "Cúmulo abierto",
    "Gb":  "Cúmulo globular",
    "Nb":  "Nebulosa difusa",
    "Pl":  "Nebulosa planetaria",
    "C+N": "Cúmulo con nebulosidad",
    "Ast": "Asterismo",
    "Kt":  "Nudo en una galaxia",
    "***": "Estrella triple",
    "D*":  "Estrella doble",
    "*":   "Estrella",
}

# Three objects carry type "-", meaning the RNGC (1973) called them
# nonexistent, so NGC 2000.0 leaves them with no type, magnitude or size.
# Dreyer and modern catalogues both disagree with that verdict, and since the
# list asks for these objects by name they are typed from his description
# rather than printed blank:
#   2478  "cluster"                 open cluster in Puppis
#   2645  "Cl, S, st L and S"       open cluster in Vela
#   4409  "vF, pS, r; = 4420?"      Dreyer's own guess; NGC 4420 is a galaxy
DSO_TYPE_FIX = {2478: "OC", 2645: "OC", 4409: "Gx"}

# Hipparcos parallaxes, borrowed from the atlas star file.  BSC5 carries a
# parallax of its own but it is pre-Hipparcos (1991) and badly wrong at these
# distances: of the 457 stars it places inside 100 ly, 198 are farther -- one
# of them (HR 4511) is really 2568 ly away.  See the log, "Distancia".
BIGSKY = ROOT / "stars.bigksy.0.1.3.mag11.parquet"

# 1 pc = 3.2616 ly; distance in ly = LY_PC / parallax_in_arcsec.
LY_PC = 3.2616

# Colour names from B-V.  Each entry is the upper bound of its bin; the
# boundaries sit on the spectral-class transitions and were checked against
# the Sp. column (azul is 90% B, amarillo 74% G, rojo 72% M).  B-V is the
# *observed* colour, not dereddened: 28 reddened O/B supergiants therefore
# read yellower than they intrinsically are, which is deliberate -- this is a
# visual catalogue, and that is the colour the star shows at the eyepiece.
COLOUR_BINS = [
    (-0.02, "azul"),
    (0.15, "azul-blanco"),
    (0.40, "blanco"),
    (0.60, "blanco-amarillo"),
    (0.85, "amarillo"),
    (1.20, "amarillo-anaranjado"),
    (1.55, "anaranjado"),
    (None, "rojo"),
]

# Spectral type as printed in the "Sp." column: the class letter and its
# numeric grade, nothing else.  Three things the raw BSC5 field forces us to
# handle (the log works through every case, under "La columna Sp."):
#   * a leading lowercase Yerkes luminosity prefix (c, d, g, sd, sg) comes
#     *before* the class -- "gK3" is a K3 giant, not a "gK" anything;
#   * Wolf-Rayet classes are two letters (WC8+O9I -> WC8);
#   * the grade may be fractional (B2.5, O9.7, S3.5) and is kept in full.
# When no grade was determined we print none: a bare "K" (HR 5044, "KIII")
# says the subclass is unknown, which is what the catalogue means.  A trailing
# m/p peculiarity marker is kept instead, since "Am"/"Ap" *is* the class.
SP_RE = re.compile(r"^[a-z]*(W[CNR]|[OBAFGKMRNSC])(\d+(?:\.\d+)?|[mp])?")

# BSC5 3-letter Bayer codes -> LaTeX math.  Omicron has no macro; plain "o".
GREEK = {
    "Alp": r"\alpha",   "Bet": r"\beta",     "Gam": r"\gamma",
    "Del": r"\delta",   "Eps": r"\varepsilon", "Zet": r"\zeta",
    "Eta": r"\eta",     "The": r"\theta",    "Iot": r"\iota",
    "Kap": r"\kappa",   "Lam": r"\lambda",   "Mu":  r"\mu",
    "Nu":  r"\nu",      "Xi":  r"\xi",       "Omi": r"o",
    "Pi":  r"\pi",      "Rho": r"\rho",      "Sig": r"\sigma",
    "Tau": r"\tau",     "Ups": r"\upsilon",  "Phi": r"\varphi",
    "Chi": r"\chi",     "Psi": r"\psi",      "Ome": r"\omega",
}


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

@dataclass
class Star:
    hr: int
    hd: str
    ra_deg: float
    de_deg: float
    vmag: float
    bayer: str          # 3-letter code, e.g. "Del"
    bayer_sup: str      # superscript digit, e.g. "1"
    flamsteed: str
    name_con: str       # constellation of the *designation*, from BSC5
    sao: str = ""       # BSC5 SAO number (bytes 32-37)
    ads: str = ""       # BSC5 ADS number (bytes 45-49); present => double
    var_id: str = ""    # BSC5 VarID (bytes 52-60), whitespace-normalised
    sp_type: str = ""   # BSC5 SpType (bytes 128-147), raw
    bv: str = ""        # BSC5 B-V (bytes 110-114), raw
    dist_ly: float | None = None   # from the Hipparcos parallax, if near enough
    con: str = ""       # IAU constellation from J2000 position
    gould: str = ""
    dbl: "Double | None" = None    # reportable visual companion, if any
    notes: list[str] = field(default_factory=list)


@dataclass
class Double:
    """The companion a double star's note reports, and its system's size."""
    sep: float          # arcsec, last measured
    pa: int             # degrees, last measured
    m2: float           # companion's magnitude
    ncomp: int          # components counted in the system (>= 2)
    struve: str = ""    # Struve designation of the system, if it has one
    wds: str = ""       # the pair's own designation, for the stderr report


@dataclass
class Dso:
    """A deep-sky row.  Shares the table's columns with Star, but only RA,
    Dec, V, Cst., HD/NGC and Notas ever carry anything."""
    ngc: int
    ra_deg: float
    de_deg: float
    vmag: float | None      # integrated magnitude, absent for 15 of them
    phot: bool              # magnitude is photographic (blue), not visual
    typ: str                # NGC 2000.0 type code
    names: list[str] = field(default_factory=list)
    con: str = ""           # IAU constellation from J2000 position
    notes: list[str] = field(default_factory=list)


def _f(s: str):
    s = s.strip()
    return float(s) if s else None


def parse_bsc5(path: Path, vmax: float, decmax: float) -> list[Star]:
    """Read the fixed-width BSC5 records that pass the magnitude/dec cuts."""
    out = []
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if len(line) < 107:
                continue
            vmag = _f(line[102:107])
            if vmag is None or vmag > vmax:
                continue
            # J2000 position, bytes 76-90
            rah, ram, ras = _f(line[75:77]), _f(line[77:79]), _f(line[79:83])
            sign = line[83]
            ded, dem, des = _f(line[84:86]), _f(line[86:88]), _f(line[88:90])
            if None in (rah, ram, ras, ded, dem, des) or sign not in "+-":
                continue  # novae / removed objects have blank positions
            de = ded + dem / 60 + des / 3600
            if sign == "-":
                de = -de
            if de >= decmax:
                continue
            ra = (rah + ram / 60 + ras / 3600) * 15.0

            name = line[4:14]
            out.append(Star(
                hr=int(line[0:4]),
                hd=line[25:31].strip(),
                ra_deg=ra, de_deg=de, vmag=vmag,
                bayer=name[3:6].strip(),
                bayer_sup=name[6:7].strip(),
                flamsteed=name[0:3].strip(),
                name_con=name[7:10].strip(),
                sao=line[31:37].strip(),
                ads=line[44:49].strip(),
                var_id=" ".join(line[51:60].split()),
                sp_type=line[127:147].strip(),
                bv=line[109:114].strip(),
            ))
    out.sort(key=lambda s: s.ra_deg)
    return out


def load_ngc_names(path: Path) -> dict[int, list[str]]:
    """NGC number -> the common names NGC 2000.0 records for it.

    The file's Comment field is either a list of co-designations ("4038-9" for
    the Antennae) or an editorial aside in parentheses.  The parenthesised ones
    are not names of the object at all -- "Beehive cluster (See Praesepe)" and
    "kappa Cru cluster (See Jewel Box)" are pointers to the entry that follows,
    and "Hourglass nebula (Brightest part of NGC 6523)" names a *part* of the
    Lagoon -- so those rows are dropped and the object keeps its real name.

    Names are printed exactly as the file spells them, "omega Cen" and all;
    only a Messier number is renormalised, from "M  31" to "M 31".
    """
    names: dict[int, list[str]] = {}
    if not path.exists():
        return names
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            obj = " ".join(line[0:35].split())
            num, comment = line[36:41].strip(), line[42:70].strip()
            if not obj or not num.isdigit() or comment.startswith("("):
                continue
            names.setdefault(int(num), []).append(obj)
    # Messier first, then the proper names in the file's own (alphabetical)
    # order: the number is the shorter and more useful handle at the eyepiece.
    for lst in names.values():
        lst.sort(key=lambda n: (0 if n.startswith("M ") else 1))
    return names


def parse_ngc2000(path: Path, wanted: list[int], decmax: float
                  ) -> tuple[list[Dso], list[tuple[int, float]]]:
    """Read the requested NGC objects.  Returns (kept, skipped-too-far-north).

    NGC 2000.0 gives RA to a tenth of a minute and Dec to the whole arcminute,
    so a deep-sky row is coarser than the stars around it by roughly 6s and
    30" respectively.  It is still printed in the table's own format: a second
    position style for 104 rows out of 2700 would cost the reader more than the
    spurious final digit does, and no object here is small enough to care.
    """
    want = set(wanted)
    kept: list[Dso] = []
    north: list[tuple[int, float]] = []
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            name = line[0:5].strip()
            if not name.isdigit() or int(name) not in want:
                continue          # IC entries keep their "I" and never match
            n = int(name)
            rah, ram = _f(line[10:12]), _f(line[13:17])
            ded, dem = _f(line[20:22]), _f(line[23:25])
            if None in (rah, ram, ded, dem) or line[19] not in "+-":
                continue
            de = ded + dem / 60.0
            if line[19] == "-":
                de = -de
            if de >= decmax:
                north.append((n, de))
                continue
            kept.append(Dso(
                ngc=n,
                ra_deg=(rah + ram / 60.0) * 15.0,
                de_deg=de,
                vmag=_f(line[40:44]),
                phot=line[44:45] == "p",
                typ=DSO_TYPE_FIX.get(n, line[6:9].strip()),
            ))
    return kept, north


def load_gould(path: Path) -> dict[str, str]:
    """HD -> Gould number, from the Uranometria Argentina.

    Only rows flagged 'G' in column 1 carry a Gould designation; the
    interleaved rows are UA stars without one.
    """
    g: dict[str, str] = {}
    if not path.exists():
        return g
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if not line.startswith("G"):
                continue
            num = line[1:5].strip()
            hd = line[64:71].strip()   # HD is right-aligned, ending at col 71
            if num and hd.isdigit():
                g.setdefault(hd, num)
    return g


def load_struve(path: Path) -> dict[str, str]:
    """SAO -> Struve designation, e.g. "28737" -> "Σ 1744".

    Struve.csv carries one row per component pair, so a system appears several
    times ("STF 1744 AB", "STF 1744 AC"); we key on SAO and keep the lowest
    designation, which is the discovery pair.
    """
    if not path.exists():
        print(f"  warning: {path} not found, no Struve notes", file=sys.stderr)
        return {}
    df = pd.read_csv(path, sep=";")
    ext = df.wds_name.str.extract(r"^(STF[AB]?)\s+(\d+)")
    df = df.assign(pfx=ext[0], num=pd.to_numeric(ext[1], errors="coerce"))
    df = df.dropna(subset=["SAO", "pfx", "num"])

    order = {"STF": 0, "STFA": 1, "STFB": 2}
    out: dict[str, tuple[int, int, str]] = {}
    for sao, pfx, num in zip(df.SAO, df.pfx, df.num):
        key = str(int(sao))
        rank = (order.get(pfx, 9), int(num))
        if key not in out or rank < out[key][:2]:
            out[key] = (*rank, f"{STRUVE_PREFIX[pfx]}{int(num)}")
    return {k: v[2] for k, v in out.items()}


def parse_coord(c: str) -> tuple[float, float]:
    """stelledoppie's "14 39 36 -60 50 02" -> (RA, Dec) in degrees.

    The sign is read from the string, not from int(): a companion at -00 20 24
    would otherwise come out north of the equator.
    """
    p = c.split()
    ra = (int(p[0]) + int(p[1]) / 60 + int(p[2]) / 3600) * 15.0
    de = abs(int(p[3])) + int(p[4]) / 60 + int(p[5]) / 3600
    return ra, -de if p[3].lstrip().startswith("-") else de


def wds_components(name: str) -> str:
    """The component token ending a designation: "STF 1110 AB" -> "AB".

    Taken from the end rather than by parsing the name, because a discoverer
    designation is not one word -- "H 5 102 AB" is Herschel's class-5 number
    102, and "DUN 252" has no component token at all.
    """
    parts = str(name).split()
    if len(parts) > 1 and re.fullmatch(r"[A-Za-z]+(,[A-Za-z]+)?", parts[-1]):
        return parts[-1]
    return ""


def split_comp(tok: str) -> tuple[set[str], set[str]]:
    """Component token -> (letters on the primary side, on the secondary side).

    With a comma the split is given: "AB,C" is the AB pair against C.  Without
    one, the first letter is the primary and the rest the secondary, so "AD" is
    A against D -- getting *that* wrong lets a faint secondary slip past the
    magnitude cut and inflates the count.  Lowercase sub-component suffixes are
    then dropped, so "Aa,Ab" is A against A: a speckle pair is one point of
    light to this book, and contributes no new component.
    """
    tok = tok or "AB"
    pri, _, sec = tok.partition(",")
    if not sec:
        pri, sec = tok[:1], tok[1:]
    return (set(re.sub("[a-z]", "", pri)), set(re.sub("[a-z]", "", sec)))


def load_doubles(path: Path) -> pd.DataFrame:
    """The stelledoppie pairs, RA-sorted, with position and components parsed.

    Rows are *pairs*, and a system is not identifiable from the file: its pairs
    do not even share a position.  Beta Monocerotis keeps its B-C row at the B
    component's own place, 5" from the A-B row's and under a different SAO, and
    alpha Crucis does the same.  So nothing is grouped here -- match_doubles
    assembles a system positionally, around the star it is looking at.
    """
    if not path.exists():
        print(f"  warning: {path} not found, no double notes", file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(path, sep=";").dropna(subset=["coord", "sep", "m1"])
    pos = [parse_coord(c) for c in df.coord]
    df = df.assign(ra=[p[0] for p in pos], de=[p[1] for p in pos],
                   comp=[wds_components(n) for n in df.wds_name])
    return df.sort_values("ra", ignore_index=True)


def match_doubles(stars: list[Star], pairs: pd.DataFrame, struve: dict[str, str],
                  radius_as: float = 20.0) -> list[Star]:
    """Set Star.dbl on the doubles, and return the component rows to delete.

    The system around a star is every pair row within **20"** -- the same
    radius, and for the same reason, as load_distances: BSC5's own position is
    off by that much on the fastest movers.  It is wide enough to pull in the
    rows registered on a *companion's* position (beta Mon's B-C) and so get the
    component count right, and those rows must not be reportable themselves,
    which is what the |V - m1| test separates: a row whose primary side is this
    star, or a row whose primary side is one of its companions.
    """
    if pairs.empty:
        return []
    pra = pairs.ra.to_numpy()
    rad = radius_as / 3600.0
    mates: dict[int, list[tuple[float, float]]] = {}   # HR -> [(sep, m2), ...]

    for s in stars:
        cos_d = max(math.cos(math.radians(s.de_deg)), 0.02)
        dra = rad / cos_d
        idx = slice(np.searchsorted(pra, s.ra_deg - dra),
                    np.searchsorted(pra, s.ra_deg + dra))
        sub = pairs.iloc[idx]
        if sub.empty:
            continue
        d_ra = (sub.ra.to_numpy() - s.ra_deg + 180) % 360 - 180
        sub = sub[np.hypot(d_ra * cos_d, sub.de.to_numpy() - s.de_deg) <= rad]
        if sub.empty:
            continue

        # the star's own pairs: those whose primary side is this star
        own = sub[(sub.m1 - s.vmag).abs() <= 0.5]
        near = own[(own.sep >= DBL_MIN_SEP) & (own.sep <= DBL_MAX_SEP)]
        good = near[near.m2 < DBL_MAX_MAG]
        if good.empty:
            continue
        best = good.loc[good.m2.idxmin()]

        # component count over the whole system, faint companions included
        letters: set[str] = set()
        for r in sub.itertuples(index=False):
            if r.sep > DBL_MAX_SEP:
                continue
            pri, sec = split_comp(r.comp)
            letters |= pri
            if pd.notna(r.m2) and r.m2 <= DBL_COUNT_MAG:
                letters |= sec

        s.dbl = Double(sep=float(best.sep), pa=int(best.pa), m2=float(best.m2),
                       ncomp=max(len(letters), 2),
                       struve=struve.get(s.sao, ""), wds=str(best.wds_name))
        # Only a companion the note could have reported may cost a row.  The
        # window is the reporting one, 3" included: a pair too tight to be
        # written up is also too tight to absorb, which is what keeps both
        # components of xi Scorpii (1.1" apart, while the note talks about the
        # C component 7" away) on the page.
        mates[s.hr] = [(float(r.sep), float(r.m2))
                       for r in near.itertuples(index=False) if pd.notna(r.m2)]

    return _companion_rows(stars, mates)


def _companion_rows(stars: list[Star], mates: dict[int, list[tuple[float, float]]]
                    ) -> list[Star]:
    """The catalogue rows that are components of a system already reported.

    Such a row says nothing its primary's note does not, so it goes.  Position
    alone cannot decide it -- theta1 and theta2 Orionis stand 135" apart and are
    different systems -- so the star's V must also match one of the system's
    companion magnitudes.  Only fainter stars are ever dropped, which settles
    the case where both components are bright enough to have picked up a note
    of their own: beta1 Tucanae keeps the row, beta2 loses it.
    """
    sra = np.array([s.ra_deg for s in stars])
    drop: dict[int, Star] = {}
    for s in stars:
        if s.dbl is None or not mates.get(s.hr):
            continue
        reach = (max(sep for sep, _ in mates[s.hr]) + 10.0) / 3600.0
        cos_d = max(math.cos(math.radians(s.de_deg)), 0.02)
        dra = reach / cos_d
        idx = range(int(np.searchsorted(sra, s.ra_deg - dra)),
                    int(np.searchsorted(sra, s.ra_deg + dra)))
        for t in (stars[i] for i in idx):
            if (t.vmag, t.hr) <= (s.vmag, s.hr):
                continue                        # the primary keeps its row
            d_ra = (t.ra_deg - s.ra_deg + 180) % 360 - 180
            if math.hypot(d_ra * cos_d, t.de_deg - s.de_deg) > reach:
                continue
            hit = [m2 for _, m2 in mates[s.hr] if abs(t.vmag - m2) <= 0.8]
            if hit:
                drop[t.hr] = t
                print(f"  double: HR {t.hr} (V {t.vmag:.2f}) absorbed into "
                      f"HR {s.hr} {s.dbl.wds}", file=sys.stderr)
            else:
                print(f"  double: HR {t.hr} (V {t.vmag:.2f}) is inside "
                      f"HR {s.hr} {s.dbl.wds} but matches no companion "
                      f"magnitude -- kept", file=sys.stderr)
    for t in drop.values():                 # a companion carries no note itself
        t.dbl = None
    return list(drop.values())


def load_cross(max_dist: float | None = None) -> dict[int, list[tuple[str, str, float]]]:
    """HD -> [(label, number, magnitude), ...] from the results/cross tables.

    Entries with magnitude 0 are dropped first: that value means the old
    catalogue recorded the star as variable, or recorded no magnitude at all.
    Where a catalogue still matches the same HD more than once, the brightest
    surviving match wins (nearest breaks a tie), so each old catalogue
    contributes at most one designation per star -- taken from the first tier
    that has the star at all, never from two tiers at once.
    """
    out: dict[int, list[tuple[str, str, float]]] = {}
    for label, tiers in CROSS_SOURCES:
        chosen: dict[int, tuple[str, float]] = {}
        for suffix, files in tiers:
            best: dict[int, tuple[str, float, float]] = {}
            for fn in files:
                p = CROSS / fn
                if not p.exists():
                    print(f"  warning: {p} not found, skipping", file=sys.stderr)
                    continue
                df = pd.read_csv(p)
                df = df[df.mag != CROSS_NO_MAG]
                if max_dist is not None:
                    df = df[df.dist <= max_dist]
                for i1, i2, mag, dist in df.itertuples(index=False):
                    try:
                        hd = int(str(i2).split()[-1])
                        num = str(i1).split()[-1]
                    except (ValueError, IndexError):
                        continue
                    mag, dist = float(mag), float(dist)
                    cur = best.get(hd)
                    if cur is None or (mag, dist) < (cur[1], cur[2]):
                        best[hd] = (num, mag, dist)
            for hd, (num, mag, _) in best.items():
                chosen.setdefault(hd, (num + suffix, mag))
        for hd, (num, mag) in chosen.items():
            out.setdefault(hd, []).append((label, num, mag))
    return out


def fmt_cross(entries: list[tuple[str, str, float]]) -> str:
    return ", ".join(f"{lab} ({num}): {mag:.1f}" for lab, num, mag in entries)


def fit_notes(free: list[str], opt: list[str]) -> str:
    """Wrap the note fragments in the LaTeX macros that fit them to the line.

    `free` is always printed; `opt` is offered in priority order and kept only
    while the accumulated line still fits the notes column.  The measuring is
    done by TeX at the real font and size -- see the preamble.
    """
    return ("\\ntstart"
            + "".join(rf"\ntfree{{{f}}}" for f in free)
            + "".join(rf"\ntopt{{{o}}}" for o in opt)
            + "\\ntend")


def load_notes(path: Path) -> dict[int, list[str]]:
    """HR -> list of remark strings, from the BSC5 notes file."""
    n: dict[int, list[str]] = {}
    if not path.exists():
        return n
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if len(line) < 13:
                continue
            hr = line[1:5].strip()
            if not hr.isdigit():
                continue
            cat = line[7:11].strip()
            txt = line[12:].rstrip()
            if txt:
                n.setdefault(int(hr), []).append(f"{cat} {txt}" if cat else txt)
    return n


# --------------------------------------------------------------------------
# Hipparcos distances
# --------------------------------------------------------------------------

def load_distances(stars: list[Star], max_ly: float, radius_as: float = 20.0) -> int:
    """Fill in Star.dist_ly from the Hipparcos parallaxes in the atlas file.

    BSC5 has no HIP number, so the join is positional -- the same kind of
    cross-match that Tycho-2 once needed, but a far cheaper one: the parquet
    holds 9487 rows down to V=6.5, which is a fraction of a second to read,
    against the 2.5 million lines Tycho-2 cost.

    Two details, both learned the hard way (see the log, section 2.3):

    * the radius is a generous **20"**, because BSC5's own position is off by
      that much for the fastest-moving stars -- alpha Cen sits 7.0" from it and
      61 Cygni 15.1", and a tight 5" radius silently dropped exactly the
      nearest stars, which are the ones this note exists for;
    * among candidates the **brightest** wins, not the nearest.  At 20" the
      nearest neighbour can be an unrelated faint star, whereas the brightest
      inside the circle is the counterpart.  In a close double this picks the
      primary for both components, which is harmless: they share a distance.
    """
    if not BIGSKY.exists():
        print(f"  warning: {BIGSKY} not found, no distance notes", file=sys.stderr)
        return 0
    import pyarrow.parquet as pq
    df = pq.read_table(BIGSKY,
                       columns=["ra", "dec", "magnitude", "parallax_mas"]).to_pandas()
    df = df[df.magnitude <= 6.5].dropna(subset=["ra", "dec", "parallax_mas"])

    tra, tde, tpx, tmg = (df.ra.to_numpy(), df.dec.to_numpy(),
                          df.parallax_mas.to_numpy(), df.magnitude.to_numpy())
    order = np.argsort(tra)
    tra, tde, tpx, tmg = tra[order], tde[order], tpx[order], tmg[order]

    rad = radius_as / 3600.0
    found = 0
    for s in stars:
        # RA window widened by 1/cos(dec); the poles are handled by the clamp
        dra = rad / max(math.cos(math.radians(s.de_deg)), 0.02)
        idx = np.arange(np.searchsorted(tra, s.ra_deg - dra),
                        np.searchsorted(tra, s.ra_deg + dra))
        if idx.size == 0:
            continue
        d_ra = (tra[idx] - s.ra_deg + 180) % 360 - 180
        sep = np.hypot(d_ra * math.cos(math.radians(s.de_deg)), tde[idx] - s.de_deg)
        near = idx[sep <= rad]
        if near.size == 0:
            continue
        px = float(tpx[near[np.argmin(tmg[near])]])   # brightest inside the radius
        if px <= 0:                     # negative parallaxes carry no distance
            continue
        ly = LY_PC / (px / 1000.0)
        if ly <= max_ly:
            s.dist_ly = ly
            found += 1
    return found


# --------------------------------------------------------------------------
# formatting
# --------------------------------------------------------------------------

def ra_hms(deg: float) -> tuple[int, int, int]:
    """Right ascension in whole h, m, s, with the rounding carry applied.

    One second of time is 15" at the equator, so the worst-case error from
    dropping the fraction is 7.5" -- finer than a fifth-magnitude catalogue
    needs, and BSC5 only carries 0.1s anyway.
    """
    h = deg / 15.0
    hh = int(h)
    mm = int((h - hh) * 60)
    ss = int(round((h - hh - mm / 60) * 3600))
    if ss >= 60:
        ss -= 60
        mm += 1
    if mm >= 60:
        mm -= 60
        hh += 1
    if hh >= 24:       # 23h59m59.7s rounds up and wraps to 0h
        hh -= 24
    return hh, mm, ss


def fmt_ra(deg: float, show_h: bool = True, show_m: bool = True) -> str:
    """e.g. 15h 32m 40s, as LaTeX with superscript units.

    A suppressed hour or minute is set as \\hphantom of itself: it prints
    nothing but keeps its width, so the seconds stay in the same column.
    """
    hh, mm, ss = ra_hms(deg)
    h_part = rf"{hh:02d}\ra{{h}}"
    m_part = rf"{mm:02d}\ra{{m}}"
    if not show_h:
        h_part = rf"\hphantom{{{h_part}}}"
    if not show_m:
        m_part = rf"\hphantom{{{m_part}}}"
    return rf"{h_part}{m_part}{ss:02d}\ra{{s}}"


def runs_visible(keys: list, page: int) -> list[bool]:
    """Within each page, keep only the first and last row of every run.

    A run of one or two rows is left untouched; a longer run keeps its ends
    and blanks everything between them.
    """
    vis = [True] * len(keys)
    for start in range(0, len(keys), page):
        end = min(start + page, len(keys))
        i = start
        while i < end:
            j = i
            while j + 1 < end and keys[j + 1] == keys[i]:
                j += 1
            for k in range(i + 1, j):
                vis[k] = False
            i = j + 1
    return vis


def fmt_dec(deg: float) -> str:
    """e.g. -20 deg 17.5', to a tenth of an arcminute (worst case 3" off)."""
    sign = "-" if deg < 0 else "+"
    a = abs(deg)
    dd = int(a)
    tenths = int(round((a - dd) * 600))    # tenths of an arcminute
    if tenths >= 600:
        tenths -= 600
        dd += 1
    whole, frac = divmod(tenths, 10)
    return rf"${sign}${dd:02d}\degr{whole:02d}\farcm{frac:01d}"


def fmt_sp(sp_type: str) -> str:
    """Spectral class and grade from a raw BSC5 SpType, e.g. K0IIIbCN-0.5 -> K0.

    Returns "" if the field is blank or unparseable; no star in the current
    selection hits either case.
    """
    m = SP_RE.match(sp_type)
    if not m:
        return ""
    return m.group(1) + (m.group(2) or "")


def colour_name(bv: str) -> str:
    """Colour word from a raw BSC5 B-V field; "" when the field is blank.

    19 of the 2596 stars have no B-V and so get no colour note, per the rule
    that nothing is printed where the colour cannot be inferred.
    """
    if not bv:
        return ""
    try:
        f = float(bv)
    except ValueError:
        return ""
    for hi, name in COLOUR_BINS:
        if hi is None or f < hi:
            return name
    return ""


def fmt_dist(ly: float) -> str:
    """e.g. "27.5 al".  Hipparcos is precise enough here to earn the decimal."""
    return f"{ly:.1f} al"


def fmt_double(d: Double) -> str:
    """e.g. "doble $\\Sigma$ 1744: d = 14.4$''$, m = 3.9, P = 153$^\\circ$".

    The position angle is printed as the whole degree the source carries: finer
    than that is more than an eyepiece can judge.  The separation keeps its
    tenth below 100", where it still distinguishes one pair from another, and
    loses it above, where it cannot.
    """
    word = MULT_WORDS.get(d.ncomp, "múltiple")
    name = f" {d.struve}" if d.struve else ""
    sep = f"{d.sep:.1f}" if d.sep < 100 else f"{d.sep:.0f}"
    return (tex_escape(f"{word}{name}: d = {sep}") + r"\arcsec"
            + tex_escape(f", m = {d.m2:.1f}, P = {d.pa}") + r"\degr")


def fmt_bayer(s: Star) -> str:
    if not s.bayer:
        return ""
    g = GREEK.get(s.bayer)
    if g is None:
        return tex_escape(s.bayer)
    return f"${g}^{{{s.bayer_sup}}}$" if s.bayer_sup else f"${g}$"


def fmt_desig(s: Star) -> str:
    """Flamsteed if present, else Gould (e.g. 43G)."""
    if s.flamsteed:
        return s.flamsteed
    if s.gould:
        return f"{s.gould}G"
    return ""


_TEX = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
        "Σ": r"$\Sigma$"}


def colour(text: str, name: str) -> str:
    """Wrap already-escaped LaTeX in a colour group."""
    return rf"\textcolor{{{name}}}{{{text}}}"


def tex_escape(t: str) -> str:
    return "".join(_TEX.get(c, c) for c in t)


# --------------------------------------------------------------------------
# LaTeX emission
# --------------------------------------------------------------------------

PREAMBLE = r"""%% twoside makes LaTeX distinguish recto (odd) from verso (even) pages,
%% which is what drives the alternating page-number position.
\documentclass[10pt,twoside]{article}

\usepackage[a4paper,top=14mm,bottom=16mm,left=12mm,right=12mm]{geometry}
\usepackage{longtable}
\usepackage{array}
\usepackage{amsmath}
\usepackage{fancyhdr}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{tikz}

%% ---- table colours ------------------------------------------------------
%% Anything not named here stays black: RA, Bayer, Flamsteed/Gould, the
%% old-catalogue notes, the headers and the rules.
\definecolor{coldec}{HTML}{1F3D7A}   %% declination -- dark blue
\definecolor{colmag}{HTML}{8B1A1A}   %% magnitude -- dark red
\definecolor{colsp}{HTML}{008B8B}    %% spectral class -- turquoise
\definecolor{coldist}{HTML}{4E6E8E}  %% distance -- slate blue
\definecolor{colcon}{HTML}{2E8B57}   %% constellation -- as on the plates
\definecolor{colhd}{HTML}{5B2C87}    %% HD -- dark violet
\definecolor{colsao}{HTML}{7A4B22}   %% SAO -- brown
\definecolor{colvar}{HTML}{B8860B}   %% variability -- dark gold
\definecolor{coldbl}{HTML}{C2185B}   %% double stars -- pink
%% lmodern must precede fontenc: BasicTeX ships no EC bitmap fonts, so plain
%% T1 Computer Modern fails to build the 5pt superscripts (ecrm0500).
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{textcomp}

%% ---- astronomical unit marks -------------------------------------------
%% Units that follow a whole number sit after it, with a thin space.
\newcommand{\ra}[1]{\textsuperscript{\textrm{#1}}\,}
\newcommand{\degr}{\ensuremath{^{\circ}}\,}
%% Units on a fractional number sit *over the decimal point*, the usual
%% A&A/AAS convention: 17\farcm5 sets the prime between the 17 and the 5.
\newcommand{\farcm}{\hbox{$.\!\!^{\prime}$}}
%% The notes column is prose, not a position, so a separation there is written
%% out plainly -- 14.4\arcsec, not 14\farcs4.
\newcommand{\arcsec}{\ensuremath{^{\prime\prime}}}

%% ---- notes that fit themselves to the line ------------------------------
%% Rather than allowing a fixed number of notes per star, we keep as many as
%% the line actually holds.  TeX does the measuring, so it is exact at the
%% real font and size -- no width table on the Python side, no extra pass.
%%
%%   \ntfree{...}  always kept   (colour, distance)
%%   \ntopt{...}   kept only while the accumulated line still fits \notew
%%
%% The first \ntopt that does not fit sets \ntstop, so everything after it is
%% skipped too: what survives is always a *prefix* of the priority order, and
%% a shorter low-priority note can never jump ahead of a dropped one.
%% \ntend unpacks with \unhcopy rather than \usebox, so in the pathological
%% case where the always-kept parts alone overrun the column the line wraps
%% (ugly but complete) instead of spilling into the margin.
\newsavebox{\ntbox}\newsavebox{\ntry}
\newif\ifntfirst \newif\ifntstop
\newcommand{\ntstart}{\sbox{\ntbox}{}\ntfirsttrue\ntstopfalse}
\newcommand{\ntfree}[1]{%
  \sbox{\ntbox}{\usebox{\ntbox}\ifntfirst\else,\space\fi#1}\ntfirstfalse}
\newcommand{\ntopt}[1]{%
  \ifntstop\else
    \sbox{\ntry}{\usebox{\ntbox}\ifntfirst\else,\space\fi#1}%
    \ifdim\wd\ntry>\notew\relax \ntstoptrue
    \else \sbox{\ntbox}{\usebox{\ntry}}\ntfirstfalse\fi
  \fi}
\newcommand{\ntend}{\unhcopy\ntbox}

%% ---- table metrics -----------------------------------------------------
\setlength{\tabcolsep}{3pt}
\renewcommand{\arraystretch}{ARRAYSTRETCH}
\setlength{\LTpre}{0pt}
\setlength{\LTpost}{0pt}
%% longtable centres each chunk by default (\LTleft=\LTright=\fill), and the
%% header block is narrower than the body, so its rules came out inset by a
%% few mm on each side.  Pinning both to 0pt aligns every chunk flush left.
\setlength{\LTleft}{0pt}
\setlength{\LTright}{0pt}

%% Width taken by the seven fixed columns plus their separators.  The notes
%% column absorbs whatever is left, so it re-adapts if the margins change.
%% Tune this one number after the first compile if the table over/underruns.
%% HDEXTRA is the allowance for the wider "HD/NGC" heading, added only when
%% deep-sky rows are present, so --fixedcolw keeps one meaning in both modes.
\newlength{\fixedcolw}\setlength{\fixedcolw}{FIXEDCOLW}
\addtolength{\fixedcolw}{HDEXTRA}
\newlength{\notew}\setlength{\notew}{\dimexpr\textwidth-\fixedcolw\relax}

%% ---- title page ---------------------------------------------------------
%% Times for the title only; the catalogue itself stays Latin Modern.  ptm is
%% the psnfss Times, present in BasicTeX (tgtermes/newtx are not).  Setting the
%% family rather than loading mathptmx keeps the change local to this page, and
%% the digits then come from the same text font as the letters -- which is what
%% makes 2596 and +52 match the surrounding capitals without any math setup.
\newcommand{\titlefont}{\fontfamily{ptm}\selectfont}

%% A printer's flourish: a rule interrupted by a symmetric knot.  #1 = width.
\newcommand{\flourishrule}[1]{%
  \begin{tikzpicture}[baseline=-0.5ex]
    \draw[line width=0.5pt] (-#1/2,0) -- (-7mm,0);
    \draw[line width=0.5pt] (7mm,0) -- (#1/2,0);
    \draw[line width=0.5pt]
      (-7mm,0) .. controls (-5mm,2.4mm) and (-2mm,1.6mm) .. (0,0)
               .. controls (2mm,-1.6mm) and (5mm,-2.4mm) .. (7mm,0);
    \draw[line width=0.5pt]
      (-7mm,0) .. controls (-5mm,-2.4mm) and (-2mm,-1.6mm) .. (0,0)
               .. controls (2mm,1.6mm) and (5mm,2.4mm) .. (7mm,0);
    \fill (0,0) circle (0.75mm);
    \fill (-#1/2,0) circle (0.4mm);
    \fill (#1/2,0) circle (0.4mm);
  \end{tikzpicture}}

%% One corner ornament, drawn for the top-left corner with the corner at (0,0);
%% the border macro mirrors it into the other three.
\newcommand{\cornerorn}{%
  \draw[line width=0.5pt]
    (0,-16mm) .. controls (0,-7mm) and (7mm,0) .. (16mm,0);
  \draw[line width=0.5pt]
    (1.6mm,-13mm) .. controls (2.4mm,-7mm) and (7mm,-2.4mm) .. (13mm,-1.6mm)
                  .. controls (8mm,-3.4mm) and (5mm,-6mm) .. (4.2mm,-10mm)
                  .. controls (3.4mm,-6.4mm) and (2.6mm,-4.6mm) .. (1.6mm,-13mm);
  \fill (2.6mm,-2.6mm) circle (0.55mm);}

%% Double frame plus the four corner ornaments, laid over the page.
\newcommand{\titleborder}{%
  \begin{tikzpicture}[remember picture, overlay]
    \coordinate (tl) at ([shift={( 17mm,-17mm)}]current page.north west);
    \coordinate (tr) at ([shift={(-17mm,-17mm)}]current page.north east);
    \coordinate (bl) at ([shift={( 17mm, 17mm)}]current page.south west);
    \coordinate (br) at ([shift={(-17mm, 17mm)}]current page.south east);
    \draw[line width=1.0pt] (tl) rectangle (br);
    \coordinate (itl) at ([shift={( 2.2mm,-2.2mm)}]tl);
    \coordinate (itr) at ([shift={(-2.2mm,-2.2mm)}]tr);
    \coordinate (ibl) at ([shift={( 2.2mm, 2.2mm)}]bl);
    \coordinate (ibr) at ([shift={(-2.2mm, 2.2mm)}]br);
    \draw[line width=0.4pt] (itl) rectangle (ibr);
    \begin{scope}[shift={(itl)}]                 \cornerorn \end{scope}
    \begin{scope}[shift={(itr)},xscale=-1]       \cornerorn \end{scope}
    \begin{scope}[shift={(ibl)},yscale=-1]       \cornerorn \end{scope}
    \begin{scope}[shift={(ibr)},xscale=-1,yscale=-1] \cornerorn \end{scope}
  \end{tikzpicture}}

%% Carry the border in the page style so it never interacts with the text.
\fancypagestyle{tpborder}{%
  \fancyhf{}%
  \renewcommand{\headrulewidth}{0pt}%
  \renewcommand{\footrulewidth}{0pt}%
  \fancyhead[C]{\titleborder}%
}

%% ---- page numbers ------------------------------------------------------
%% Outer edge of the spread: right on odd (recto) pages, left on even
%% (verso) ones, as a printed and bound book wants them.
\fancypagestyle{catalogo}{%
  \fancyhf{}%
  \fancyfoot[RO,LE]{\thepage}%
  \renewcommand{\headrulewidth}{0pt}%
  \renewcommand{\footrulewidth}{0pt}%
}

\begin{document}

%% ---- front matter, unnumbered ------------------------------------------
\pagestyle{empty}

%% title (recto).  The border is drawn from the page style, not from the text:
%% a tikzpicture placed in the body joins the title's paragraph, and the
%% \fill then lands inside the block instead of above it.
\thispagestyle{tpborder}
\null\vspace*{\fill}
{\titlefont
\begin{center}
{\fontsize{25}{34}\selectfont
 \begin{minipage}{0.74\textwidth}\centering TITLETEXT\end{minipage}}\\[13mm]
\flourishrule{78mm}\\[11mm]
{\fontsize{15}{20}\selectfont AUTHORTEXT}
\end{center}}
\vspace*{\fill}
\newpage

%% deliberately blank (verso), so the catalogue opens on a recto
\null\newpage

%% ---- catalogue, numbered from 1 ----------------------------------------
\pagestyle{catalogo}
\setcounter{page}{1}

FONTSIZE
%% \extrarowheight lifts the text off the rule above it, which the superscript
%% h/m of the RA column would otherwise touch.  \arraystretch opens up both
%% sides of the baseline.
\setlength{\extrarowheight}{EXTRAROWHEIGHT}

\begin{longtable}{@{}l l r l l c l r r p{\notew}@{}}
HEADBLOCK"""

# The header row.  When it lives in \endhead, longtable boxes it separately and
# sizes that box from the header's *own* natural column widths -- so wherever a
# heading is narrower than the data below it (HD, SAO, V ...) the head box comes
# out short, its rules stop early, and the headings drift left of their columns.
# Emitting the row as ordinary table rows puts it in the same alignment pass as
# the data, which fixes both symptoms at once.
HEADER_ROW = (r"\textbf{AR} & \textbf{Dec} & \textbf{V} & \textbf{Sp.} & "
              r"\textbf{Cst.} & \textbf{B.} & \textbf{Fl/G} & \textbf{HDHEAD} & "
              r"\textbf{SAO} & \textbf{Notas} \\")

HEADER_BLOCK = "\\hline\n" + HEADER_ROW + "\n\\hline"

# Only used with --flow, where we do not control the page breaks and so have to
# fall back on longtable repeating the head for us.
HEADBLOCK_FLOW = ("\\hline\n" + HEADER_ROW + "\n\\hline\n\\endfirsthead\n"
                  "\\hline\n" + HEADER_ROW + "\n\\hline\n\\endhead")

POSTAMBLE = r"""\end{longtable}
ATLAS
\end{document}
"""


# Extra width the bold "HD/NGC" heading claims over the widest HD number.
# Measured from the overfull-hbox warning the plain 78/90mm settings produce
# once deep-sky rows are in: 14.42pt = 5.06mm, rounded up for a little slack.
HD_NGC_EXTRA = "5.5mm"

# Space a plate may occupy on the page, in mm (text block minus a little).
ATLAS_MAX_W = 186.0
ATLAS_MAX_H = 264.0


def pdf_bbox(path: Path) -> tuple[float, float]:
    """(width, height) of a PDF's MediaBox, in points."""
    m = re.search(rb"/MediaBox\s*\[([^\]]*)\]", path.read_bytes())
    if not m:
        return (1.0, 1.0)
    a, b, c, d = (float(v) for v in m.group(1).split()[:4])
    return (c - a, d - b)


def atlas_block(maps_dir: Path) -> str:
    """LaTeX for the atlas plates that follow the catalogue.

    Band plates are rotated 90 degrees clockwise so RA reads bottom-to-top and
    Dec left-to-right.  \\rotatebox is used rather than the angle= key of
    \\includegraphics because graphicx scales *before* it rotates, which makes
    width= mean the pre-rotation width and is easy to get wrong.

    The two polar halves are set as a spread: the left-opening half on a verso
    pushed to the spine, the right-opening half on the facing recto.
    """
    if not maps_dir.is_dir():
        return ""
    rel = maps_dir.name
    out = [r"\clearpage"]
    for name in ("map1", "map2", "map3"):
        plate = maps_dir / f"{name}.pdf"
        if not plate.exists():
            continue
        # After the -90 rotation the plate's own width becomes its height on
        # the page, so size from the measured box rather than a fixed number
        # -- the Klein frame changes the aspect and a hard-coded value silently
        # overruns the text width.
        w, h = pdf_bbox(plate)
        long_mm = min(ATLAS_MAX_H, ATLAS_MAX_W * w / h)
        out += [
            r"\begin{center}",
            rf"\rotatebox{{-90}}{{\includegraphics[width={long_mm:.1f}mm]"
            rf"{{{rel}/{name}.pdf}}}}",
            r"\end{center}",
            r"\clearpage",
        ]
    # polar spread: force the left half onto a verso so the circle closes
    # across the gutter rather than across a page turn
    pa, pb = maps_dir / "polar_a.pdf", maps_dir / "polar_b.pdf"
    if pa.exists() and pb.exists():
        def polar_h(path):
            w, h = pdf_bbox(path)
            return min(ATLAS_MAX_H, ATLAS_MAX_W * h / w)
        out += [
            r"\ifodd\value{page}\null\clearpage\fi",
            rf"\noindent\hfill\includegraphics[height={polar_h(pb):.1f}mm]"
            rf"{{{rel}/polar_b.pdf}}",
            r"\clearpage",
            rf"\noindent\includegraphics[height={polar_h(pa):.1f}mm]"
            rf"{{{rel}/polar_a.pdf}}\hfill\mbox{{}}",
            r"\clearpage",
        ]
    return "\n".join(out)


def star_cells(s: Star) -> list[str]:
    return [
        "",                                     # RA, filled in by the caller
        colour(fmt_dec(s.de_deg), "coldec"),
        # \rlap keeps the asterisk out of the cell's measured width, so
        # every magnitude stays flush right and the decimal points line up
        colour(f"{s.vmag:.1f}", "colmag")
        + (r"\rlap{" + colour(r"$^{*}$", "colvar") + "}" if s.var_id else ""),
        colour(fmt_sp(s.sp_type), "colsp") if s.sp_type else "",
        colour(s.con, "colcon"),
        fmt_bayer(s),
        fmt_desig(s),
        colour(s.hd, "colhd") if s.hd else "",
        colour(s.sao, "colsao") if s.sao else "",
        " ".join(s.notes),
    ]


def dso_cells(d: Dso) -> list[str]:
    """The same ten columns, with the six that mean nothing here left empty.

    A dagger on the magnitude marks the photographic (blue) values: NGC 2000.0
    carries whichever of the two it has, and 13 of these objects were never
    measured visually.  It is \\rlap'd for the same reason the variable-star
    asterisk is -- so the digits stay flush right with the stars above.
    """
    if d.vmag is None:
        mag = ""
    else:
        mag = colour(f"{d.vmag:.1f}", "colmag")
        if d.phot:
            mag += r"\rlap{" + colour(r"$^{\dagger}$", "colmag") + "}"
    return [
        "",                                     # RA, filled in by the caller
        colour(fmt_dec(d.de_deg), "coldec"),
        mag,
        "",                                     # Sp.
        colour(d.con, "colcon"),
        "", "",                                 # B., Fl/G
        colour(str(d.ngc), "colhd"),
        "",                                     # SAO
        " ".join(d.notes),
    ]


def emit(stars: list, out: Path, per_page: int, fontsize: str,
         arraystretch: float, notewidth: str, flow: bool = False,
         extrarowheight: str = "2pt", rule_every: int = 5,
         title: str = "", repeat_ra: bool = False,
         atlas: str = "", author: str = "") -> None:
    # The mid-table headers are appended to the body, which never goes through
    # the preamble's placeholder pass, so they need their own substitution.
    # The HD column only earns its second name when deep-sky rows are present.
    # It is not cosmetic: bold "HD/NGC" is wider than any six-digit HD number,
    # so it is the heading, not the data, that then sets the column's width.
    # The measured cost is 14.42pt (5.06mm); HD_NGC_EXTRA rounds that up and is
    # added to \fixedcolw here rather than being left for --fixedcolw to carry,
    # so the same number works whether or not the deep-sky rows are in.
    has_dso = any(isinstance(r, Dso) for r in stars)
    hd_head = "HD/NGC" if has_dso else "HD"
    head_block = HEADER_BLOCK.replace("HDHEAD", hd_head)
    # Ink-saving rule: within a page, an hour or minute that repeats down the
    # column is printed only on the first and last row of the run.  Minutes are
    # keyed on (hour, minute) so a run can never straddle an hour boundary.
    hms = [ra_hms(s.ra_deg) for s in stars]
    if repeat_ra:
        show_h = show_m = [True] * len(stars)
    else:
        page = per_page if (per_page and not flow) else len(stars) or 1
        show_h = runs_visible([h for h, _, _ in hms], page)
        show_m = runs_visible([(h, m) for h, m, _ in hms], page)

    body = []
    for i, s in enumerate(stars):
        cells = dso_cells(s) if isinstance(s, Dso) else star_cells(s)
        cells[0] = fmt_ra(s.ra_deg, show_h[i], show_m[i])
        body.append(" & ".join(cells) + r" \\")
        # a rule every rule_every rows, banding the table instead of boxing
        # every star; always close the table with one
        if (i + 1) % rule_every == 0 or i + 1 == len(stars):
            body.append(r"\hline")
        # force the page break exactly on the requested cadence.  Once the
        # notes column is populated this will overrun the page in note-dense
        # regions; --flow lets longtable break naturally instead.
        if not flow and per_page and (i + 1) % per_page == 0 and i + 1 < len(stars):
            body.append(r"\newpage")
            body.append(head_block)

    post = POSTAMBLE.replace("ATLAS", atlas)
    tex = (PREAMBLE.replace("HEADBLOCK",
                            HEADBLOCK_FLOW if flow else HEADER_BLOCK)
           .replace("ARRAYSTRETCH", f"{arraystretch:.2f}")
           .replace("FONTSIZE", fontsize)
           .replace("FIXEDCOLW", notewidth)
           .replace("HDEXTRA", HD_NGC_EXTRA if has_dso else "0pt")
           .replace("EXTRAROWHEIGHT", extrarowheight)
           .replace("HDHEAD", hd_head)
           .replace("TITLETEXT", title)
           .replace("AUTHORTEXT", author.upper())
           + "\n".join(body) + "\n" + post)
    out.write_text(tex, encoding="utf-8")


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vmax", type=float, default=5.5)
    ap.add_argument("--decmax", type=float, default=52.0)
    ap.add_argument("--per-page", type=int, default=55)
    ap.add_argument("--pages", type=int, default=2,
                    help="number of pages to emit; 0 = the whole catalogue")
    ap.add_argument("--out", type=Path, default=BOOK / "catalog.tex")
    ap.add_argument("--notes", choices=["cross", "bsc5", "none"], default="cross",
                    help="notes source: old-catalogue cross-identifications "
                         "from results/cross (default), the BSC5 remarks, "
                         "or an empty column")
    ap.add_argument("--max-dist", type=float, default=None,
                    help="drop cross-matches farther than this many arcsec")
    ap.add_argument("--max-notes", type=int, default=0,
                    help="optional hard cap on notes per star; 0 (default) "
                         "lets the line width decide, keeping as many as fit "
                         "and dropping the lowest-priority ones first")
    ap.add_argument("--dso", choices=["ngc", "none"], default="ngc",
                    help="merge the hand-picked NGC 2000.0 deep-sky objects "
                         "into the table (default); none leaves the "
                         "catalogue stars-only")
    ap.add_argument("--doubles", choices=["csv", "none"], default="csv",
                    help="report a star's brightest visual companion instead "
                         "of its old-catalogue designations, and absorb the "
                         "companion's own row; none leaves both alone")
    ap.add_argument("--parallax", choices=["bigsky", "none"], default="bigsky",
                    help="source of the distance note; none omits it")
    ap.add_argument("--max-dist-ly", type=float, default=100.0,
                    help="stars farther than this get no distance note")
    ap.add_argument("--note-chars", type=int, default=80,
                    help="truncate demo notes to this many characters; 80 is "
                         "the measured break-even at 50 stars/page, i.e. the "
                         "point where every note still fits on one line")
    ap.add_argument("--fontsize", default=r"\scriptsize")
    ap.add_argument("--arraystretch", type=float, default=1.30)
    ap.add_argument("--atlas", action="store_true",
                    help="append the star-atlas plates from book/maps/")
    ap.add_argument("--repeat-ra", action="store_true",
                    help="print the RA hour/minute on every row instead of "
                         "only at the ends of each repeated run")
    ap.add_argument("--author", default="Daniel Esteban Severin",
                    help="name set below the title-page flourish")
    ap.add_argument("--title", default="",
                    help="override the title-page text (LaTeX)")
    ap.add_argument("--rule-every", type=int, default=5,
                    help="draw a horizontal rule every N rows (1 = one per star)")
    ap.add_argument("--extrarowheight", default="2pt",
                    help="space added above the text in every row, to keep "
                         "the RA superscripts clear of the rule above")
    ap.add_argument("--flow", action="store_true",
                    help="let longtable break pages naturally instead of "
                         "forcing a break every --per-page rows")
    ap.add_argument("--fixedcolw", dest="notewidth", default="90mm",
                    help="width allowance for the nine fixed columns; the "
                         "notes column takes textwidth minus this, less "
                         "HD_NGC_EXTRA when deep-sky rows widen the HD heading")
    args = ap.parse_args()

    stars = parse_bsc5(BSC5, args.vmax, args.decmax)
    # the whole catalogue, before --pages truncates it for a preview build;
    # the title page describes the book, not the excerpt
    total = len(stars)
    print(f"BSC5: {total} stars with V <= {args.vmax} and "
          f"Dec < +{args.decmax}", file=sys.stderr)

    # Visual doubles, before anything else looks at the list: a star that is a
    # component of a system already reported loses its row here.  `total` above
    # deliberately keeps counting it -- the book still covers that star, in its
    # primary's note and with its own dot on the atlas plates.
    struve = load_struve(STRUVE) if args.doubles == "csv" else {}
    if args.doubles == "csv":
        dropped = match_doubles(stars, load_doubles(DOUBLES), struve)
        gone = {s.hr for s in dropped}
        stars = [s for s in stars if s.hr not in gone]
        ndbl = sum(1 for s in stars if s.dbl)
        words: dict[str, int] = {}
        for s in stars:
            if s.dbl:
                w = MULT_WORDS.get(s.dbl.ncomp, "múltiple")
                words[w] = words.get(w, 0) + 1
        print(f"doubles: {ndbl} stars with a companion note "
              f"({', '.join(f'{v} {k}' for k, v in sorted(words.items(), key=lambda kv: -kv[1]))}); "
              f"{len(dropped)} component rows absorbed", file=sys.stderr)

    gould = load_gould(UA)
    for s in stars:
        if s.hd and not s.flamsteed:
            s.gould = gould.get(s.hd, "")
    print(f"Gould designations available: "
          f"{sum(1 for s in stars if s.gould)}", file=sys.stderr)
    n_sp = sum(1 for s in stars if fmt_sp(s.sp_type))
    print(f"spectral classes parsed: {n_sp}/{total}", file=sys.stderr)

    dsos: list[Dso] = []
    if args.dso == "ngc":
        dsos, north = parse_ngc2000(NGC2000, NGC_LIST, args.decmax)
        for n, de in north:
            print(f"  NGC {n} is at Dec {de:+.1f}, north of the "
                  f"{args.decmax:+g} limit -- omitted", file=sys.stderr)
        names = load_ngc_names(NGC2000_NAMES)
        for d in dsos:
            d.names = names.get(d.ngc, [])
        found = {d.ngc for d in dsos} | {n for n, _ in north}
        for n in sorted(set(NGC_LIST) - found):
            print(f"  NGC {n} not found in {NGC2000.name}", file=sys.stderr)
        print(f"deep-sky objects: {len(dsos)} of {len(set(NGC_LIST))} requested",
              file=sys.stderr)
    total_dso = len(dsos)

    # One RA-ordered sequence: a deep-sky object sits among the stars that
    # share its patch of sky, which is how it will be looked up.
    rows = sorted(stars + dsos, key=lambda r: r.ra_deg)
    if args.pages:
        rows = rows[: args.pages * args.per_page]
    stars = [r for r in rows if isinstance(r, Star)]
    dsos = [r for r in rows if isinstance(r, Dso)]

    if args.parallax == "bigsky" and args.notes == "cross":
        n = load_distances(stars, args.max_dist_ly)
        print(f"Hipparcos distances within {args.max_dist_ly:g} ly: "
              f"{n}/{len(stars)}", file=sys.stderr)

    # IAU constellation from the J2000 position
    from astropy.coordinates import SkyCoord, get_constellation
    import astropy.units as u
    coords = SkyCoord([r.ra_deg for r in rows] * u.deg,
                      [r.de_deg for r in rows] * u.deg, frame="icrs")
    for r, c in zip(rows, get_constellation(coords, short_name=True)):
        r.con = c

    # Deep-sky notes: what the object *is*, then what it is called.  The type
    # is never dropped; the names are offered to the same width rule the stars
    # use, so a long one gives way rather than wrapping the row.
    for d in dsos:
        typ = DSO_TYPE.get(d.typ, "")
        d.notes = [fit_notes([tex_escape(typ)] if typ else [],
                             [tex_escape(n) for n in d.names])]
    n_phot = sum(1 for d in dsos if d.phot)
    n_untyped = sum(1 for d in dsos if d.typ not in DSO_TYPE)
    if dsos:
        print(f"deep-sky: {sum(1 for d in dsos if d.vmag is not None)} with a "
              f"magnitude ({n_phot} photographic), "
              f"{sum(1 for d in dsos if d.names)} named, "
              f"{n_untyped} without a type", file=sys.stderr)

    if args.notes == "cross":
        cross = load_cross(args.max_dist)
        struve = struve or load_struve(STRUVE)
        hit = nstruve = ndblstruve = 0
        for s in stars:
            # "free" notes are never dropped: they are short, and they answer
            # what the star *is* rather than what it was once called.
            free = []
            cn = colour_name(s.bv)
            if cn:
                free.append(colour(tex_escape(cn), "colsp"))
            if s.dist_ly is not None:
                free.append(colour(tex_escape(fmt_dist(s.dist_ly)), "coldist"))
            # "opt" notes are offered in priority order and kept while the
            # line holds them; LaTeX decides, not this loop.
            opt = []
            if s.dbl:
                # What the companion is and where to look for it displaces
                # every designation the star would otherwise have carried: at
                # the eyepiece the pair is the fact, and the Struve number --
                # the one designation still worth printing -- rides along
                # inside the note.
                free.append(colour(fmt_double(s.dbl), "coldbl"))
                if s.dbl.struve:
                    ndblstruve += 1
            elif s.hd:
                # A BSC5 ADS number means the star is a known double; if Struve
                # also catalogued it, lead with his designation.
                if s.ads and s.sao:
                    st = struve.get(s.sao)
                    if st:
                        opt.append(colour(tex_escape(st), "coldbl"))
                        nstruve += 1
                # Only a real variable-star designation earns a note; the
                # asterisk in the V column already covers the rest.
                vn = var_note(s.var_id)
                if vn:
                    opt.append(colour(tex_escape(vn), "colvar"))
                # each cross-identification is its own fragment, so they can be
                # dropped one at a time from the low-priority end
                ent = cross.get(int(s.hd), [])
                if args.max_notes:      # optional hard cap, off by default
                    ent = ent[: max(args.max_notes - len(opt), 0)]
                opt += [tex_escape(f"{lab} ({num}): {mag:.1f}")
                        for lab, num, mag in ent]
            if free or opt:
                s.notes = [fit_notes(free, opt)]
                hit += 1
        ndouble = sum(1 for s in stars if s.ads)
        ncol = sum(1 for s in stars if colour_name(s.bv))
        ndist = sum(1 for s in stars if s.dist_ly is not None)
        print(f"stars with a note: {hit}/{len(stars)};  "
              f"ADS doubles {ndouble}, of which Struve {nstruve} on their own"
              f" + {ndblstruve} inside a companion note", file=sys.stderr)
        print(f"colour names: {ncol}/{len(stars)};  "
              f"distances within {args.max_dist_ly:g} ly: {ndist}", file=sys.stderr)
    elif args.notes == "bsc5":
        allnotes = load_notes(BSC5_NOTES)
        for s in stars:
            txt = " ".join(allnotes.get(s.hr, []))
            if len(txt) > args.note_chars:
                txt = txt[: args.note_chars].rsplit(" ", 1)[0] + "..."
            s.notes = [tex_escape(txt)] if txt else []

    # "magnitude 5" is the traditional magnitude class: V <= 5.5 is the whole
    # of the 5th magnitude, so floor() gives the number a reader expects.
    # Upper case, and the degree sign as text (\textdegree) rather than a math
    # superscript so it comes from the same Times face as the digits.
    title = args.title or (
        (f"Catálogo y Atlas de {total} estrellas hasta la quinta magnitud"
         + (f" y {total_dso} objetos de cielo profundo," if total_dso else " y")
         + f" al sur de la declinación {args.decmax:+g}").upper()
        + r"\textdegree{}")

    emit(rows, args.out, args.per_page, args.fontsize,
         args.arraystretch, args.notewidth, args.flow, args.extrarowheight,
         args.rule_every, title, args.repeat_ra,
         atlas_block(MAPS) if args.atlas else "", args.author)

    n_mismatch = sum(1 for s in stars if s.name_con and s.con != s.name_con)
    print(f"wrote {args.out} - {len(stars)} stars + {len(dsos)} deep-sky, "
          f"{math.ceil(len(rows)/args.per_page)} pages", file=sys.stderr)
    print(f"designation/IAU constellation mismatches: {n_mismatch}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
