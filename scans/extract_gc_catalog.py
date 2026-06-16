#!/usr/bin/env python3
"""
Extract star table pages from the Catalogo General Argentino (RNAO14.pdf)
and convert them to CSV (star number, magnitude, cross-reference).

v2 — column-strip composites:
  1. Render printed page with pdftoppm at 300 DPI (printed N = PDF N+19)
  2. Detect the table's ruled lines (max vertical ink-run length); locate
     the N°, Mag, N°Obs and reference columns
  3. Build two narrow composites: [N° | Mag] and [N° | a | d | Reference],
     split into short segments cut at blank rows
  4. Ask Gemini (one call per composite) for {"num", "mag"/"ref"} pairs
  5. Validate: star numbers consecutive, first number chains to the
     previous page, every magnitude equals cat/gc.txt; retry on failure.
     References accepted when two runs agree (third run breaks ties).
  6. Write rnao14_pageNN.csv

Usage:
  python3 extract_gc_catalog.py --validate           # pages 62-70 vs known CSVs
  python3 extract_gc_catalog.py --pages 71 90        # write CSVs
  python3 extract_gc_catalog.py --pages 71 0         # from 71 until quota runs out
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
from PIL import Image

try:
    from google import genai
except ImportError:
    sys.exit("Install google-genai: pip install google-genai")

# -- configuration ------------------------------------------------------------

PDF         = Path(__file__).parent / "RNAO14.pdf"
GC_TXT      = Path(__file__).parent.parent / "cat" / "gc.txt"
OUT_DIR     = Path(__file__).parent   # write CSVs into scans/
CACHE_DIR   = Path("/tmp/gc_pages")
def _load_api_key():
    """Return the Google API key, kept OUTSIDE this file so it is never
    committed. Reads $GOOGLE_API_KEY if set, otherwise ../../gapikey.txt
    (i.e. the parent directory of the repository, e.g. /Users/<you>/gapikey.txt)."""
    env = os.environ.get("GOOGLE_API_KEY")
    if env:
        return env.strip()
    key_file = Path(__file__).resolve().parent.parent.parent / "gapikey.txt"
    try:
        return key_file.read_text().strip()
    except OSError:
        sys.exit(f"No Google API key: set $GOOGLE_API_KEY or put it in {key_file}")

API_KEY     = _load_api_key()
MODEL       = "models/gemini-3.1-flash-lite"
PAGE_OFFSET = 19      # printed page + 19 = PDF page
CALL_DELAY  = 5       # seconds between API calls (15 RPM limit)
MAX_RUNS    = 3
SEG_HEIGHT  = 700     # max pixel height of a composite segment

MAG_PROMPT = """\
The images are consecutive vertical segments (top to bottom) of a two-column
composite cut out of a 19th-century star catalog table.
Left column: the star number (4 digits, increasing down the page; blank on
continuation rows). Right column: the magnitude printed on the SAME
horizontal line as the star number. A magnitude is an integer ("9"), a
decimal ("7.0"), or an integer followed by a printed fraction glyph;
transcribe the fractions 1/4, 1/2, 3/4 with a space: "8 1/2", "7 1/4",
"9 3/4". Occasionally the cell holds "var." or is empty: use null then.
Ignore any header text. A segment boundary may clip a row; report each star
number exactly once, in order.
Respond with a JSON array only:
[{"num":3018,"mag":"7 1/4"},{"num":3019,"mag":"7.0"}]"""

REF_PROMPT = """\
The images are consecutive vertical segments (top to bottom) of a composite
of five columns cut out of a 19th-century star catalog table. Light gray
horizontal guide lines separate consecutive printed lines. The columns are:
(1) star number (4 digits, increasing down the page; blank on continuation
lines), (2) a date like "73.88" (printed on EVERY line, also on continuation
lines), (3) and (4) two small observation-count integers, (5) a
cross-reference, usually empty. Examples of references: "L. 884", "L. (884)",
"Ll. 5315", "ZC. 1220", "OA. 1846", "B. 425", "WB. 788", "P. 209", "G. 79",
"F. 4", "Y. 1291", "Melb.I. 155". They may carry lowercase annotations
("red", "comes, s.", "sq.", "pr.", "n. sq."): transcribe them verbatim too.
Distinguish "L." from "Ll." carefully. Keep parentheses.
Go line by line (the date column marks every line). For EVERY star number
give the reference printed between the same guide lines, or "" when that
cell is empty. Never move a reference to a neighbouring line. If a
reference sits on a continuation line (no star number), report it as
{"num": null, "ref": "..."}.
Ignore any header text. A segment boundary may clip a row; report each star
number exactly once, in order.
Respond with a JSON array only:
[{"num":3018,"ref":"L. (884)"},{"num":3019,"ref":""}]"""

CACHE_DIR.mkdir(parents=True, exist_ok=True)
client = genai.Client(api_key=API_KEY)

FRACTIONS = {"1/4": "2", "1/2": "5", "3/4": "7"}


class QuotaExhausted(Exception):
    pass


# -- gc.txt magnitudes ----------------------------------------------------------

def load_gc_magnitudes() -> dict[int, float | None]:
    """First-observation magnitudes; None for variables/blank cells."""
    mags: dict[int, float | None] = {}
    for line in open(GC_TXT):
        try:
            n, obs = int(line[0:5]), int(line[5:7])
        except ValueError:
            continue
        if obs != 1:
            continue
        m = line[7:10].strip()
        if m and line[10:11] != "V":
            mags[n] = int(m) / 10 if len(m) == 3 else float(m)
        else:
            mags[n] = None
    return mags


# -- image processing ------------------------------------------------------------

def render_page(printed_page: int) -> Image.Image:
    pdf_page = printed_page + PAGE_OFFSET
    out = CACHE_DIR / f"page_{printed_page:03d}.png"
    if not out.exists():
        subprocess.run(
            ["pdftoppm", "-f", str(pdf_page), "-l", str(pdf_page), "-r", "300",
             "-png", "-singlefile", str(PDF), str(out.with_suffix(""))],
            capture_output=True, check=True,
        )
    return Image.open(out).convert("L")


def _max_runs(d: np.ndarray) -> np.ndarray:
    """Max consecutive-True run length per column of a 2-D bool array."""
    h = d.shape[0]
    idx = np.arange(h)[:, None]
    last_false = np.maximum.accumulate(np.where(~d, idx, -1), axis=0)
    runs = np.where(d, idx - last_false, 0)
    return runs.max(axis=0)


def _group(positions: np.ndarray, tol: int) -> list[tuple[int, int]]:
    groups: list[list[int]] = []
    for p in positions:
        if groups and p - groups[-1][-1] <= tol:
            groups[-1].append(int(p))
        else:
            groups.append([int(p)])
    return [(g[0], g[-1]) for g in groups]


def detect_layout(img: Image.Image) -> dict:
    """Locate table body and the N°, Mag, obs and reference column x-ranges."""
    a = np.asarray(img) < 200
    h, w = a.shape

    vert = _max_runs(a)
    groups = _group(np.where(vert > 140)[0], tol=10)
    if len(groups) < 10:
        raise ValueError(f"only {len(groups)} vertical rule groups")
    mid = lambda g: (g[0] + g[1]) / 2

    # Border pair: the two heavy rules at a plausible table width apart,
    # enclosing as many detected rules as possible (ignores stamps and
    # marginalia on either side of the table).
    best_pair, best_n = None, -1
    for gr in reversed(groups):
        for gl in groups:
            if 1900 <= mid(gr) - mid(gl) <= 2090:
                n = sum(1 for g in groups if gl[0] < g[0] < gr[0])
                if n > best_n:
                    best_pair, best_n = (gl, gr), n
    if best_pair is None:
        raise ValueError("table borders not found")
    gl, gr = best_pair
    width_table = mid(gr) - mid(gl)

    # Snap the 12 internal separators to their template positions
    # (fractions of the table width measured on page 62), synthesizing
    # separators whose rules are too faint to detect.
    TEMPLATE = [0.0618, 0.1735, 0.2132, 0.2696, 0.3735, 0.4421,
                0.5103, 0.6107, 0.6873, 0.7564, 0.7863, 0.8143]
    inner = [g for g in groups if gl[0] < g[0] < gr[0]]
    rules = [gl]
    used = set()
    for frac in TEMPLATE:
        want = mid(gl) + frac * width_table
        cand = [(abs(mid(g) - want), i, g) for i, g in enumerate(inner)
                if i not in used and abs(mid(g) - want) <= 28]
        if cand:
            d, i, g = min(cand)
            used.add(i)
            rules.append(g)
        else:
            x = int(want)
            rules.append((x - 2, x + 2))
    rules.append(gr)
    rules.sort()

    # Horizontal borders: rows with a long contiguous dark run (true rules;
    # footer/footnote text lines have short runs)
    y0 = y1 = None
    hrun = _max_runs(a.T)
    hgroups = _group(np.where(hrun > 0.18 * w)[0], tol=12)
    if len(hgroups) >= 2:
        top = hgroups[0]
        bottoms = [g for g in hgroups if 2000 < g[0] - top[1] < 2620]
        if bottoms:
            y0, y1 = top[1] + 3, bottoms[-1][0] - 3
    if y0 is None:
        # Faint rules: take whichever border rule was detected and fill in
        # the other side from the table geometry (top ~260-320, height
        # ~2550 px on every page of this volume).
        tops = [g for g in hgroups if g[0] < 600]
        bottoms = [g for g in hgroups if g[0] > 2400]
        y0 = tops[0][1] + 3 if tops else 320
        y1 = bottoms[-1][0] - 3 if bottoms else min(y0 + 2550, 3100)

    # White gaps between rule groups (slivers between double rules dropped);
    # crops expand a little into the rule zones so leaning rules never clip
    # cell contents.
    gaps, crops = [], []
    for i in range(len(rules) - 1):
        g = (rules[i][1] + 2, rules[i + 1][0] - 1)
        gaps.append(g)
        x0c = rules[i][1] - 6
        crops.append((x0c, max(rules[i + 1][0] + 7, x0c + 10)))
    if len(gaps) != 13:
        raise ValueError(f"{len(gaps)} column gaps instead of 13: {gaps}")
    width = lambda g: g[1] - g[0]
    crops[2] = (crops[2][0] - 14, crops[2][1])   # '10' can overflow Mag cell
    if not (60 <= width(gaps[0]) <= 260):
        raise ValueError(f"suspicious N column width {gaps[0]}")
    if not (35 <= width(gaps[2]) <= 160):
        raise ValueError(f"suspicious Mag column width {gaps[2]}")
    if not (180 <= width(gaps[12]) <= 520):
        raise ValueError(f"suspicious reference column width {gaps[12]}")
    obs_ok = all(12 <= width(g) <= 95 for g in gaps[10:12])
    return {"y0": y0, "y1": y1, "cols": crops, "obs_ok": obs_ok}


def _row_profile(img: Image.Image, gap, y0, y1) -> np.ndarray:
    a = np.asarray(img.crop((gap[0], y0, gap[1], y1))) < 128
    return a.sum(axis=1).astype(float)


def _best_shift(prof: np.ndarray, lattice: np.ndarray, max_shift=30) -> int:
    """Vertical offset s such that prof shifted up by s best matches lattice."""
    best, bs = -1.0, 0
    n = len(prof)
    for s in range(-max_shift, max_shift + 1):
        a = prof[max(0, s): n + min(0, s)]
        b = lattice[max(0, -s): n + min(0, -s)]
        v = float(np.dot(a, b))
        if v > best:
            best, bs = v, s
    return bs


PART_COLS = {"num": [0], "mag": [2], "fecha": [3], "obs": [10, 11], "ref": [12]}


def compose(img: Image.Image, layout: dict, parts: list[str]) -> Image.Image:
    """
    Horizontally paste the requested column crops, separated by rules.
    The print can be vertically misregistered across the page (a reference
    may sit half a line lower than its star number, and pages can lean), so
    column crops are re-aligned by chaining ink-profile cross-correlation
    through ALL physical columns: every hop is between adjacent columns and
    tightly bounded, which corrects large total drift while making one-row
    jumps impossible.
    """
    y0, y1 = layout["y0"], layout["y1"]
    cols = layout["cols"]
    profiles = [_row_profile(img, g, y0, y1) for g in cols]
    shifts, cum, prev = [0], 0, 0
    for i in range(1, len(cols)):
        if profiles[i].sum() > 0 and profiles[prev].sum() > 0:
            cum += _best_shift(profiles[i], profiles[prev], max_shift=10)
            prev = i
        shifts.append(cum)
    indices = []
    for p in parts:
        if p == "obs" and not layout["obs_ok"]:
            continue
        indices.extend(PART_COLS[p])
    crops = [img.crop((cols[i][0], y0 + shifts[i], cols[i][1], y1 + shifts[i]))
             for i in indices]
    height = y1 - y0
    width = sum(c.width for c in crops) + 3 * (len(crops) - 1)
    out = Image.new("L", (width, height), 255)
    x = 0
    for i, c in enumerate(crops):
        out.paste(c, (x, 0))
        x += c.width
        if i < len(crops) - 1:
            out.paste(0, (x, 0, x + 1, height))
            x += 3
    return _add_row_guides(out)


def _add_row_guides(strip: Image.Image) -> Image.Image:
    """Draw a light line in each blank band between text rows."""
    a = np.asarray(strip) < 128
    row_ink = a.sum(axis=1)
    blanks = _group(np.where(row_ink == 0)[0], tol=1)
    for b0, b1 in blanks:
        if b1 - b0 >= 5:
            y = (b0 + b1) // 2
            strip.paste(190, (0, y, strip.width, y + 1))
    return strip


def split_segments(strip: Image.Image) -> list[Image.Image]:
    """Split a tall strip into segments <= SEG_HEIGHT, cutting at blank rows."""
    a = np.asarray(strip) < 128
    row_ink = a.sum(axis=1)
    segs, top = [], 0
    h = strip.height
    while h - top > SEG_HEIGHT:
        target = top + SEG_HEIGHT
        cut = None
        for off in range(0, 250):
            y = target - off
            if y <= top + 50:
                break
            if row_ink[max(0, y - 4):y + 4].sum() == 0:
                cut = y
                break
        if cut is None:
            cut = target
        segs.append(strip.crop((0, top, strip.width, cut)))
        top = cut
    segs.append(strip.crop((0, top, strip.width, h)))
    return [s.resize((int(s.width * 1.5), int(s.height * 1.5)), Image.LANCZOS)
            for s in segs]


# -- Gemini ----------------------------------------------------------------------

api_calls = 0


def ask_gemini(prompt: str, images: list[Image.Image]) -> list[dict]:
    """One API call with backoff on 503; raises QuotaExhausted on daily 429."""
    global api_calls
    for attempt in range(4):
        time.sleep(CALL_DELAY)
        try:
            api_calls += 1
            raw = client.models.generate_content(
                model=MODEL, contents=[prompt, *images]).text
            text = re.sub(r"^```[a-z]*\n?", "", raw.strip()).rstrip("`").strip()
            rows = json.loads(text)
            if not isinstance(rows, list):
                raise ValueError("not a JSON array")
            return rows
        except Exception as e:
            msg = str(e)
            if "RESOURCE_EXHAUSTED" in msg or "429" in msg:
                if "PerDay" in msg or "per_day" in msg.lower():
                    raise QuotaExhausted(msg[:200])
                print("    429 per-minute, waiting 60s", file=sys.stderr)
                time.sleep(60)
                continue
            if "UNAVAILABLE" in msg or "503" in msg:
                print(f"    503, backing off ({attempt+1})", file=sys.stderr)
                time.sleep(20)
                continue
            raise
    raise RuntimeError("gave up after repeated API errors")


# -- cleaning & validation ---------------------------------------------------------

def convert_mag(raw: str) -> str:
    raw = (raw or "").strip()
    m = re.match(r"^(\d+)\s+(1/4|1/2|3/4)$", raw)
    if m:
        return f"{m.group(1)}.{FRACTIONS[m.group(2)]}"
    return raw


def clean_ref(raw: str) -> str:
    """Drop lowercase annotation tokens, fix OCR slips, normalize spacing."""
    raw = (raw or "").strip()
    tokens = []
    for t in raw.replace(",", " ").split():
        if re.match(r"^[a-z.…'’]+$", t):        # lowercase annotation
            continue
        t = re.sub(r"^L[I1l]\.", "Ll.", t)       # LI./L1. -> Ll.
        t = re.sub(r"^0A\.", "OA.", t)
        tokens.append(t)
    ref = " ".join(tokens)
    ref = re.sub(r"\.(?=\S)", ". ", ref)
    ref = re.sub(r"^Melb\.\s+I\.", "Melb.I.", ref)
    ref = re.sub(r"^O\.\s*A\.", "OA.", ref)
    ref = re.sub(r"^Z\.?\s*C\.", "ZC.", ref)
    ref = re.sub(r"^W\.?\s*B\.", "WB.", ref)
    return ref.strip()


def norm(s: str) -> str:
    return re.sub(r"\s+", "", s)


def check_mags(pairs: dict[int, str], mags: dict[int, float | None],
               anchor: int | None) -> list[str]:
    """Consecutive numbering, chain anchor, magnitudes equal gc.txt."""
    errors = []
    if not pairs:
        return ["no stars extracted"]
    nums = sorted(pairs)
    if nums != list(range(nums[0], nums[-1] + 1)):
        missing = sorted(set(range(nums[0], nums[-1] + 1)) - set(nums))
        errors.append(f"numbers not consecutive, missing {missing}")
    if anchor is not None and nums[0] != anchor:
        errors.append(f"first star {nums[0]} != expected {anchor}")
    if not 35 <= len(nums) <= 80:
        errors.append(f"implausible star count {len(nums)}")
    for n in nums:
        if n not in mags:
            errors.append(f"{n}: not in gc.txt")
            continue
        want = mags[n]
        if want is None:
            continue                      # variable/blank in gc.txt: accept
        try:
            got = float(convert_mag(pairs[n]))
        except ValueError:
            errors.append(f"{n}: unparseable magnitude '{pairs[n]}'")
            continue
        if abs(got - want) > 0.001:
            errors.append(f"{n}: magnitude {got} != gc.txt {want}")
    return errors


# -- per-page pipeline --------------------------------------------------------------

def extract_mags(img, layout, mags, anchor, verbose=True):
    segs = split_segments(compose(img, layout, ["num", "mag"]))
    best, best_err = {}, ["no successful run"]
    for run in range(MAX_RUNS):
        try:
            rows = ask_gemini(MAG_PROMPT, segs)
        except QuotaExhausted:
            raise
        except Exception as e:
            if verbose:
                print(f"    mag run {run+1}: {e}", file=sys.stderr)
            continue
        pairs = {}
        for r in rows:
            if r.get("num") is not None:
                pairs.setdefault(int(r["num"]), (r.get("mag") or "").strip())
        errors = check_mags(pairs, mags, anchor)
        if not best or len(errors) < len(best_err):
            best, best_err = pairs, errors
        if not errors:
            break
        if verbose:
            print(f"    mag run {run+1}: {len(errors)} issue(s): {errors[:3]}",
                  file=sys.stderr)
    return best, best_err


def extract_refs(img, layout, expected, verbose=True):
    segs = split_segments(compose(img, layout, ["num", "fecha", "obs", "ref"]))
    runs: list[dict[int, str]] = []
    notes: list[str] = []
    for run in range(MAX_RUNS):
        try:
            rows = ask_gemini(REF_PROMPT, segs)
        except QuotaExhausted:
            raise
        except Exception as e:
            if verbose:
                print(f"    ref run {run+1}: {e}", file=sys.stderr)
            continue
        pairs, orphans = {}, []
        for r in rows:
            if r.get("num") is not None:
                pairs.setdefault(int(r["num"]), clean_ref(r.get("ref") or ""))
            elif clean_ref(r.get("ref") or ""):
                orphans.append(clean_ref(r["ref"]))
        if orphans:
            notes.append(f"text on unnumbered rows: {orphans}")
        runs.append(pairs)
        if set(pairs) != set(expected):
            if verbose:
                print(f"    ref run {run+1}: star set mismatch "
                      f"(missing {sorted(set(expected)-set(pairs))[:5]}, "
                      f"extra {sorted(set(pairs)-set(expected))[:5]})",
                      file=sys.stderr)
            continue
        if len(runs) >= 2 and \
                all(runs[-1].get(n) == runs[-2].get(n) for n in expected):
            return pairs, notes
    if not runs:
        return {}, notes + ["no successful ref run"]
    final, disagreements = {}, []
    for n in expected:
        votes = Counter(r[n] for r in runs if n in r)
        if not votes:
            final[n] = ""
            disagreements.append(f"{n}: never extracted")
            continue
        ref, cnt = votes.most_common(1)[0]
        final[n] = ref
        if cnt < 2 and len(votes) > 1:
            disagreements.append(f"{n}: no majority {list(votes)}")
    if disagreements:
        notes.append(f"ref disagreements: {disagreements}")
    return final, notes


def process_page(page: int, mags, anchor=None, verbose=True):
    img = render_page(page)
    layout = detect_layout(img)
    mag_pairs, mag_errors = extract_mags(img, layout, mags, anchor, verbose)
    expected = sorted(mag_pairs)
    ref_pairs, ref_notes = extract_refs(img, layout, expected, verbose)
    stars = {n: {"mag": mag_pairs.get(n, ""), "ref": ref_pairs.get(n, "")}
             for n in expected}
    return stars, mag_errors + ref_notes


def write_csv(page: int, stars: dict[int, dict]) -> Path:
    path = OUT_DIR / f"rnao14_page{page}.csv"
    with open(path, "w") as f:
        for n in sorted(stars):
            f.write(f"{n}, {convert_mag(stars[n]['mag'])}, {stars[n]['ref']}\n")
    return path


def read_known(page: int) -> dict[int, tuple[str, str]] | None:
    path = OUT_DIR / f"rnao14_page{page}.csv"
    if not path.exists():
        return None
    known = {}
    for line in path.read_text().splitlines():
        if line.strip():
            n, mag, ref = [p.strip() for p in line.split(",", 2)]
            known[int(n)] = (mag, ref)
    return known


def compare_with_known(page: int, stars: dict[int, dict]) -> list[str]:
    known = read_known(page)
    if known is None:
        return [f"no ground truth for page {page}"]
    diffs = []
    for n in sorted(set(known) | set(stars)):
        if n not in stars:
            diffs.append(f"{n}: missing from extraction")
            continue
        if n not in known:
            diffs.append(f"{n}: extra star not in ground truth")
            continue
        kmag, kref = known[n]
        gmag, gref = convert_mag(stars[n]["mag"]), stars[n]["ref"]
        if gmag != kmag:
            diffs.append(f"{n}: mag '{gmag}' != known '{kmag}'")
        if norm(gref) != norm(kref):
            diffs.append(f"{n}: ref '{gref}' != known '{kref}'")
    return diffs


# -- CLI ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--validate", action="store_true",
                      help="printed pages 62-70, compare against known CSVs (no writes)")
    mode.add_argument("--pages", nargs=2, type=int, metavar=("FIRST", "LAST"),
                      help="process printed pages FIRST..LAST and write CSVs "
                           "(LAST=0: continue until quota runs out)")
    args = parser.parse_args()

    mags = load_gc_magnitudes()
    last_star = max(mags)

    try:
        if args.validate:
            total_stars = total_diffs = 0
            anchor = None
            for page in range(62, 71):
                known = read_known(page)
                anchor = min(known) if known else anchor
                print(f"Page {page} ...", flush=True)
                stars, issues = process_page(page, mags, anchor)
                diffs = compare_with_known(page, stars)
                total_stars += len(stars)
                total_diffs += len(diffs)
                status = "OK" if not diffs and not issues else "ISSUES"
                print(f"  {len(stars)} stars, issues: {len(issues)}, "
                      f"diffs vs ground truth: {len(diffs)}  [{status}]")
                for d in diffs:
                    print(f"    diff {d}")
                for e in issues:
                    print(f"    note {e}")
            acc = 100.0 * (1 - total_diffs / max(1, 2 * total_stars))
            print(f"\nTOTAL: {total_stars} stars, {total_diffs} field diffs "
                  f"(field accuracy ~{acc:.2f}%), {api_calls} API calls")
        else:
            first, last = args.pages
            if last == 0:
                last = 10 ** 6
            prev = read_known(first - 1)
            anchor = max(prev) + 1 if prev else None
            page = first
            failures = 0
            while page <= last:
                if anchor is not None and anchor > last_star:
                    print(f"Catalog exhausted (last star {last_star}).")
                    break
                print(f"Page {page} ...", flush=True)
                try:
                    stars, issues = process_page(page, mags, anchor)
                except QuotaExhausted:
                    raise
                except Exception as e:
                    stars, issues = {}, [f"{type(e).__name__}: {e}"]
                if not stars:
                    failures += 1
                    print(f"  SKIPPED page {page}: {issues}")
                    if failures >= 3:
                        print("  3 consecutive failures — stopping")
                        break
                    anchor = None
                    page += 1
                    continue
                failures = 0
                path = write_csv(page, stars)
                flag = f"  ISSUES: {issues}" if issues else ""
                print(f"  written {path.name} ({len(stars)} stars), "
                      f"{api_calls} API calls so far{flag}", flush=True)
                anchor = max(stars) + 1
                page += 1
    except QuotaExhausted as e:
        print(f"\nQUOTA EXHAUSTED after {api_calls} API calls: {e}")
        sys.exit(2)


if __name__ == "__main__":
    main()
