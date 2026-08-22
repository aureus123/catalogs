#!/usr/bin/env python3
"""Generate the LaTeX source for the star catalogue.

Primary source is the Bright Star Catalogue 5th ed. (book/ybsc5.txt).
Gould designations come from the Uranometria Argentina (cat/ua.txt), joined
on HD.  Positions are optionally refined against Tycho-2 (cat/tyc2.txt), since
BSC5 only carries RA to 0.1s and Dec to 1".

Usage:
    python build_catalog.py                    # 2 pages, 50 stars/page
    python build_catalog.py --pages 0          # whole catalogue
    python build_catalog.py --notes            # fill notes from ybsc5.notes.txt
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
# Oeltzen-Argelander is split across two files (southern / northern zones);
# they reuse each other's numbers but share no HD, so merging is unambiguous.
CROSS_SOURCES = [
    ("Gould",       ["cross_gc_hd.csv"]),
    ("Yarnall",     ["cross_usno_hd.csv"]),
    ("Oeltzen-Arg", ["cross_oa_hd.csv", "cross_oarn_hd.csv"]),
    ("Gilliss",     ["cross_gilliss_hd.csv"]),
    ("BAC",         ["cross_bac_hd.csv"]),
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

BSC5 = BOOK / "ybsc5.txt"
BSC5_NOTES = BOOK / "ybsc5.notes.txt"
UA = CAT / "ua.txt"
TYC2 = CAT / "tyc2.txt"
TYC2_SUPPL = CAT / "tyc2_suppl.txt"
TYC_CACHE = BOOK / ".tyc2_bright.parquet"
MAPS = BOOK / "maps"      # atlas plates produced by gen_maps.py

SUPPL1_N = 17588  # supplement-2 (probably-false stars) starts after this

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
    con: str = ""       # IAU constellation from J2000 position
    gould: str = ""
    notes: list[str] = field(default_factory=list)
    pos_src: str = "bsc5"


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
            ))
    out.sort(key=lambda s: s.ra_deg)
    return out


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


def load_cross(max_dist: float | None = None) -> dict[int, list[tuple[str, str, float]]]:
    """HD -> [(label, number, magnitude), ...] from the results/cross tables.

    Entries with magnitude 0 are dropped first: that value means the old
    catalogue recorded the star as variable, or recorded no magnitude at all.
    Where a catalogue still matches the same HD more than once, the brightest
    surviving match wins (nearest breaks a tie), so each old catalogue
    contributes at most one designation per star.
    """
    out: dict[int, list[tuple[str, str, float]]] = {}
    for label, files in CROSS_SOURCES:
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
            out.setdefault(hd, []).append((label, num, mag))
    return out


def fmt_cross(entries: list[tuple[str, str, float]]) -> str:
    return ", ".join(f"{lab} ({num}): {mag:.1f}" for lab, num, mag in entries)


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
# Tycho-2 positional refinement
# --------------------------------------------------------------------------

def build_tyc_cache(decmax: float, vtmax: float = 7.5) -> pd.DataFrame:
    """Scan Tycho-2 + supplement-1 once, keep the bright end, cache to parquet."""
    rows = []
    with open(TYC2, encoding="latin-1") as fh:
        for line in fh:
            p = line.split("|")
            vt = p[19].strip()
            if not vt or float(vt) > vtmax:
                continue
            ra = p[2].strip() or p[24].strip()
            de = p[3].strip() or p[25].strip()
            if not ra or not de or float(de) >= decmax + 1.0:
                continue
            rows.append((float(ra), float(de), float(vt)))
    with open(TYC2_SUPPL, encoding="latin-1") as fh:
        for i, line in enumerate(fh):
            if i >= SUPPL1_N:
                break
            p = line.split("|")
            vt = p[13].strip()
            if not vt or float(vt) > vtmax:
                continue
            ra, de = p[2].strip(), p[3].strip()
            if not ra or not de or float(de) >= decmax + 1.0:
                continue
            rows.append((float(ra), float(de), float(vt)))
    df = pd.DataFrame(rows, columns=["ra", "de", "vt"]).sort_values("ra")
    df.to_parquet(TYC_CACHE, index=False)
    return df


def refine_positions(stars: list[Star], decmax: float, radius_as: float = 6.0):
    """Replace BSC5 positions with Tycho-2 ones where an unambiguous match exists."""
    if TYC_CACHE.exists():
        tyc = pd.read_parquet(TYC_CACHE)
    else:
        print("  building Tycho-2 bright-star cache (one-time, ~1 min)...",
              file=sys.stderr)
        tyc = build_tyc_cache(decmax)

    tra = tyc.ra.to_numpy()
    tde = tyc.de.to_numpy()
    tvt = tyc.vt.to_numpy()
    order = np.argsort(tra)
    tra, tde, tvt = tra[order], tde[order], tvt[order]

    rad = radius_as / 3600.0
    matched = 0
    for s in stars:
        # RA window, widened by 1/cos(dec); handle wraparound by searching twice
        dra = rad / max(math.cos(math.radians(s.de_deg)), 0.02)
        lo, hi = s.ra_deg - dra, s.ra_deg + dra
        idx = np.arange(np.searchsorted(tra, lo), np.searchsorted(tra, hi))
        if lo < 0:
            idx = np.r_[idx, np.arange(np.searchsorted(tra, lo + 360), len(tra))]
        if hi > 360:
            idx = np.r_[idx, np.arange(0, np.searchsorted(tra, hi - 360))]
        if idx.size == 0:
            continue
        d_de = tde[idx] - s.de_deg
        d_ra = (tra[idx] - s.ra_deg + 180) % 360 - 180
        sep = np.hypot(d_ra * math.cos(math.radians(s.de_deg)), d_de)
        near = idx[sep <= rad]
        if near.size == 0:
            continue
        best = near[np.argmin(tvt[near])]   # brightest inside the radius
        s.ra_deg, s.de_deg = float(tra[best]), float(tde[best])
        s.pos_src = "tyc2"
        matched += 1
    return matched


# --------------------------------------------------------------------------
# formatting
# --------------------------------------------------------------------------

def ra_hms(deg: float) -> tuple[int, int, float]:
    """Right ascension in h, m, s, with the rounding carry already applied."""
    h = deg / 15.0
    hh = int(h)
    mm = int((h - hh) * 60)
    ss = (h - hh - mm / 60) * 3600
    if round(ss, 2) >= 60.0:
        ss -= 60.0
        mm += 1
    if mm >= 60:
        mm -= 60
        hh += 1
    return hh, mm, ss


def fmt_ra(deg: float, show_h: bool = True, show_m: bool = True) -> str:
    """e.g. 15h 32m 40.32s, as LaTeX with superscript units.

    A suppressed hour or minute is set as \\hphantom of itself: it prints
    nothing but keeps its width, so the seconds stay in the same column.
    """
    hh, mm, ss = ra_hms(deg)
    whole, frac = divmod(round(ss * 100), 100)
    h_part = rf"{hh:02d}\ra{{h}}"
    m_part = rf"{mm:02d}\ra{{m}}"
    if not show_h:
        h_part = rf"\hphantom{{{h_part}}}"
    if not show_m:
        m_part = rf"\hphantom{{{m_part}}}"
    return rf"{h_part}{m_part}{whole:02d}\fs{frac:02d}"


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
    """e.g. -20 deg 23' 44.1''."""
    sign = "-" if deg < 0 else "+"
    a = abs(deg)
    dd = int(a)
    mm = int((a - dd) * 60)
    ss = (a - dd - mm / 60) * 3600
    if round(ss, 1) >= 60.0:
        ss -= 60.0
        mm += 1
    if mm >= 60:
        mm -= 60
        dd += 1
    whole, frac = divmod(round(ss * 10), 10)
    return rf"${sign}${dd:02d}\degr{mm:02d}\arcm{whole:02d}\farcs{frac:01d}"


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
\newcommand{\arcm}{\ensuremath{'}\,}
%% Units on a fractional number sit *over the decimal point*, the usual
%% A&A/AAS convention: 35\fs68 sets the "s" between the 35 and the 68.
\newcommand{\fs}{\hbox{$.\!\!^{\mathrm{s}}$}}
\newcommand{\farcs}{\hbox{$.\!\!^{\prime\prime}$}}

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
\newlength{\fixedcolw}\setlength{\fixedcolw}{FIXEDCOLW}
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

\begin{longtable}{@{}l l r l c l rSAOSPEC p{\notew}@{}}
HEADBLOCK"""

# The header row.  When it lives in \endhead, longtable boxes it separately and
# sizes that box from the header's *own* natural column widths -- so wherever a
# heading is narrower than the data below it (HD, SAO, V ...) the head box comes
# out short, its rules stop early, and the headings drift left of their columns.
# Emitting the row as ordinary table rows puts it in the same alignment pass as
# the data, which fixes both symptoms at once.
HEADER_ROW = (r"\textbf{AR (J2000)} & \textbf{Dec (J2000)} & \textbf{V} & "
              r"\textbf{Cst.} & \textbf{B.} & \textbf{Fl/G} & \textbf{HD} &"
              r"SAOHEAD \textbf{Notas} \\")

HEADER_BLOCK = "\\hline\n" + HEADER_ROW + "\n\\hline"

# Only used with --flow, where we do not control the page breaks and so have to
# fall back on longtable repeating the head for us.
HEADBLOCK_FLOW = ("\\hline\n" + HEADER_ROW + "\n\\hline\n\\endfirsthead\n"
                  "\\hline\n" + HEADER_ROW + "\n\\hline\n\\endhead")

POSTAMBLE = r"""\end{longtable}
ATLAS
\end{document}
"""


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


def emit(stars: list[Star], out: Path, per_page: int, fontsize: str,
         arraystretch: float, notewidth: str, flow: bool = False,
         extrarowheight: str = "2pt", rule_every: int = 5,
         sao: bool = False, title: str = "", repeat_ra: bool = False,
         atlas: str = "", author: str = "") -> None:
    # resolve the optional SAO cell once; the mid-table headers are appended to
    # the body, which never goes through the preamble's placeholder pass
    sao_head = r" \textbf{SAO} &" if sao else ""
    head_block = HEADER_BLOCK.replace("SAOHEAD", sao_head)
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
        note = " ".join(s.notes)
        cells = [
            fmt_ra(s.ra_deg, show_h[i], show_m[i]),
            colour(fmt_dec(s.de_deg), "coldec"),
            # \rlap keeps the asterisk out of the cell's measured width, so
            # every magnitude stays flush right and the decimal points line up
            colour(f"{s.vmag:.1f}", "colmag")
            + (r"\rlap{" + colour(r"$^{*}$", "colvar") + "}" if s.var_id else ""),
            colour(s.con, "colcon"),
            fmt_bayer(s),
            fmt_desig(s),
            colour(s.hd, "colhd") if s.hd else "",
            *([colour(s.sao, "colsao") if s.sao else ""] if sao else []),
            note,
        ]
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
           .replace("EXTRAROWHEIGHT", extrarowheight)
           .replace("SAOSPEC", " r" if sao else "")
           .replace("SAOHEAD", r" \textbf{SAO} &" if sao else "")
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
    ap.add_argument("--positions", choices=["bsc5", "tycho2"], default="tycho2",
                    help="tycho2 refines BSC5's 0.1s/1'' positions (default)")
    ap.add_argument("--notes", choices=["cross", "bsc5", "none"], default="cross",
                    help="notes source: old-catalogue cross-identifications "
                         "from results/cross (default), the BSC5 remarks, "
                         "or an empty column")
    ap.add_argument("--max-dist", type=float, default=None,
                    help="drop cross-matches farther than this many arcsec")
    ap.add_argument("--max-notes", type=int, default=3,
                    help="most cross-identifications shown per star; the "
                         "lowest-priority ones are dropped first")
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
    ap.add_argument("--sao", action="store_true",
                    help="add an SAO column to the right of HD")
    ap.add_argument("--rule-every", type=int, default=5,
                    help="draw a horizontal rule every N rows (1 = one per star)")
    ap.add_argument("--extrarowheight", default="2pt",
                    help="space added above the text in every row, to keep "
                         "the RA superscripts clear of the rule above")
    ap.add_argument("--flow", action="store_true",
                    help="let longtable break pages naturally instead of "
                         "forcing a break every --per-page rows")
    ap.add_argument("--fixedcolw", dest="notewidth", default="78mm",
                    help="width allowance for the 7 fixed columns; the notes "
                         "column takes textwidth minus this")
    args = ap.parse_args()

    stars = parse_bsc5(BSC5, args.vmax, args.decmax)
    # the whole catalogue, before --pages truncates it for a preview build;
    # the title page describes the book, not the excerpt
    total = len(stars)
    print(f"BSC5: {total} stars with V <= {args.vmax} and "
          f"Dec < +{args.decmax}", file=sys.stderr)

    gould = load_gould(UA)
    for s in stars:
        if s.hd and not s.flamsteed:
            s.gould = gould.get(s.hd, "")
    print(f"Gould designations available: "
          f"{sum(1 for s in stars if s.gould)}", file=sys.stderr)

    if args.pages:
        stars = stars[: args.pages * args.per_page]

    if args.positions == "tycho2":
        m = refine_positions(stars, args.decmax)
        print(f"positions refined from Tycho-2: {m}/{len(stars)}", file=sys.stderr)

    # IAU constellation from the J2000 position
    from astropy.coordinates import SkyCoord, get_constellation
    import astropy.units as u
    coords = SkyCoord([s.ra_deg for s in stars] * u.deg,
                      [s.de_deg for s in stars] * u.deg, frame="icrs")
    for s, c in zip(stars, get_constellation(coords, short_name=True)):
        s.con = c

    if args.notes == "cross":
        cross = load_cross(args.max_dist)
        struve = load_struve(STRUVE)
        hit = nstruve = 0
        for s in stars:
            if not s.hd:
                continue
            parts = []
            # A BSC5 ADS number means the star is a known double; if Struve
            # also catalogued it, lead with his designation.
            if s.ads and s.sao:
                st = struve.get(s.sao)
                if st:
                    parts.append(colour(tex_escape(st), "coldbl"))
                    nstruve += 1
            # Only a real variable-star designation earns a leading note; the
            # asterisk in the V column already covers the rest.
            vn = var_note(s.var_id)
            if vn:
                parts.append(colour(tex_escape(vn), "colvar"))
            # keep at most --max-notes entries in total, in priority order
            ent = cross.get(int(s.hd), [])[: max(args.max_notes - len(parts), 0)]
            if ent:
                parts.append(tex_escape(fmt_cross(ent)))
            if parts:
                s.notes = [", ".join(parts)]
                hit += 1
        ndouble = sum(1 for s in stars if s.ads)
        print(f"stars with a note: {hit}/{len(stars)};  "
              f"ADS doubles {ndouble}, of which Struve {nstruve}", file=sys.stderr)
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
        (f"Catálogo y Atlas de {total} estrellas hasta la quinta magnitud "
         f"y al sur de la declinación {args.decmax:+g}").upper()
        + r"\textdegree{}")

    emit(stars, args.out, args.per_page, args.fontsize,
         args.arraystretch, args.notewidth, args.flow, args.extrarowheight,
         args.rule_every, args.sao, title, args.repeat_ra,
         atlas_block(MAPS) if args.atlas else "", args.author)

    n_mismatch = sum(1 for s in stars if s.name_con and s.con != s.name_con)
    print(f"wrote {args.out} - {len(stars)} stars, "
          f"{math.ceil(len(stars)/args.per_page)} pages", file=sys.stderr)
    print(f"designation/IAU constellation mismatches: {n_mismatch}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
