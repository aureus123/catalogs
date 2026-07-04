#!/usr/bin/env python3
"""
Extract the Dpl./Color footnotes of the Cordoba Durchmusterung (RNAO16.pdf)
into one JSON file per printed page (scans/rnao16_pageNN.json).

Pipeline per printed page (printed N = PDF page N+51):
  1. Render with pdftoppm at 300 DPI
  2. Detect the table boxes from their vertical rules; cut out
       - the page header (declination label, e.g. "-23°")
       - each table's header row (star ranges like "16401-16430")
       - every text band outside the tables (footnotes "Dpl.: Nº 737.",
         "Color: Nº 721, 799.", zone headings "ZONA -23°", printer marks)
  3. One Gemini call transcribes all strips; the code parses declination,
     page star interval and footnote numbers
  4. Cross-checks: footnote numbers inside the page star interval, star
     ranges contiguous, numbers <= last star of the declination (cd_vol1.txt)
  5. Up to 3 runs; accepted when two agree, otherwise majority vote and
     numbers seen in only one run are stored as uncertain ("?")

Output JSON (a list: transition pages hold two declinations):
  [{"page":59,"decl":-23,"first":621,"last":950,"dpl":[737],
    "color":[721,799],"dpl_uncertain":[],"color_uncertain":[],"notes":[]}]

Usage:
  python3 extract_cd_footnotes.py --pages 59 63          # printed pages
  python3 extract_cd_footnotes.py --pages 1 0            # until quota runs out
  python3 extract_cd_footnotes.py --pages 56 56 --force --debug
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
from PIL import Image, ImageDraw, ImageFont

try:
    from google import genai
except ImportError:
    sys.exit("Install google-genai: pip install google-genai")

# -- configuration ------------------------------------------------------------

PDF         = Path(__file__).parent / "RNAO16.pdf"
CD_VOL1     = Path(__file__).parent.parent / "cat" / "cd_vol1.txt"
OUT_DIR     = Path(__file__).parent          # write JSONs into scans/
CACHE_DIR   = Path("/tmp/cd_pages300")

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
PAGE_OFFSET = 51      # printed page + 51 = PDF page
CALL_DELAY  = 5       # seconds between API calls (15 RPM limit)
MAX_RUNS    = 3
DECS        = range(22, 32)

PROMPT = """\
You receive several images cut from one scanned page of the Cordoba
Durchmusterung star catalog (19th century, printed in Spanish).
Image 1 is the page header area. The remaining images are, in top-to-bottom
page order, either:
- a table header strip: star-number ranges, one per column group, like
  "16401-16430" or "1-40",
- a text strip outside the tables: footnotes like "Dpl.: Nº 737, 1348." or
  "Color: Nº 721, 799.", a zone heading like "ZONA -23°", or printer marks.
Each image carries its index printed as "[n]" in the top-left corner of a
white margin: report its content under that image number and do not
transcribe the marker itself.
Transcribe every printed line of black text in each image exactly as
printed, digit by digit (the digits are small; be careful to distinguish
1/4/7, 3/8, 5/6, 0/9). Write the numero sign as "Nº". Ignore the dashed
blue library stamp, handwriting and specks. Keep the lines of one image in
top-to-bottom order.
Respond with JSON only, one entry per image, e.g.:
[{"image":1,"lines":["-23°   - 59 -   1h-2h"]},
 {"image":2,"lines":["621-650   651-680"]},
 {"image":3,"lines":["Dpl.: Nº 737.","Color: Nº 721, 799."]}]"""

CACHE_DIR.mkdir(parents=True, exist_ok=True)
client = genai.Client(api_key=API_KEY)


class QuotaExhausted(Exception):
    pass


# -- last star of each declination (cd_vol1.txt) -------------------------------

def load_max_star() -> dict[int, int]:
    mx: dict[int, int] = {}
    for line in open(CD_VOL1):
        if not line.startswith("CD-"):
            continue
        d, n = int(line[3:5]), int(line[5:10])
        mx[d] = max(mx.get(d, 0), n)
    return mx


MAX_STAR = load_max_star()


# -- image processing ------------------------------------------------------------

def render_page(printed_page: int) -> Image.Image:
    """600 DPI render (Gemini reads the strips at full resolution; the
    layout analysis works on a half-size copy)."""
    pdf_page = printed_page + PAGE_OFFSET
    out = CACHE_DIR / f"page_{printed_page:03d}.png"
    if not out.exists():
        subprocess.run(
            ["pdftoppm", "-f", str(pdf_page), "-l", str(pdf_page), "-r", "600",
             "-png", "-singlefile", str(PDF), str(out.with_suffix(""))],
            capture_output=True, check=True,
        )
    return Image.open(out)


def _group(positions: np.ndarray, tol: int) -> list[tuple[int, int]]:
    groups: list[list[int]] = []
    for p in positions:
        if groups and p - groups[-1][-1] <= tol:
            groups[-1].append(int(p))
        else:
            groups.append([int(p)])
    return [(g[0], g[-1]) for g in groups]


def _rule_mask(dark: np.ndarray, min_len: int) -> np.ndarray:
    """Bool mask of pixels belonging to vertical ink runs >= min_len px
    (table rules)."""
    h, w = dark.shape
    down_all = np.zeros((h, w), np.int32)
    up_all = np.zeros((h, w), np.int32)
    run = np.zeros(w, np.int32)
    for y in range(h):
        run = (run + 1) * dark[y]
        down_all[y] = run
    run = np.zeros(w, np.int32)
    for y in range(h - 1, -1, -1):
        run = (run + 1) * dark[y]
        up_all[y] = run
    return ((down_all + up_all - 1) >= min_len) & dark


def cut_strips(img: Image.Image, debug_tag: str | None = None) -> list[tuple[str, Image.Image]]:
    """Return [(kind, crop)] with kinds 'header', 'tablehead', 'text',
    in top-to-bottom page order (header first). img is the 600 DPI render;
    the analysis runs on a half-size copy."""
    small = img.reduce(2)
    L = np.asarray(small.convert("L"))
    h, w = L.shape
    x0 = 130                       # skip the blue ADS stamp / page edges
    dark = L[:, x0:w - 30] < 200

    rules = _rule_mask(dark, min_len=70)
    # A row belongs to a table when it crosses >= 3 rule columns: the outer
    # frame that encloses the whole page on zone-transition pages only
    # contributes 2 (its left and right border), so the text between the
    # tables (mid-page footnotes, "ZONA -23°" headings) stays outside.
    groups_per_row = (rules[:, 1:] & ~rules[:, :-1]).sum(axis=1) + rules[:, 0]
    in_table = groups_per_row >= 3
    spans = [s for s in _group(np.where(in_table)[0], tol=10) if s[1] - s[0] >= 50]
    if not spans:
        raise ValueError("no table found on page")

    hrules = _rule_mask(dark.T, 60).T      # horizontal rules (table borders)
    text = dark & ~rules & ~hrules         # ink that is not a table rule
    ink = text.sum(axis=1)
    masked = np.zeros(h, bool)
    for s0, s1 in spans:
        masked[max(0, s0 - 6):s1 + 7] = True
    cand = (ink >= 8) & ~masked
    bands = []
    for b0, b1 in _group(np.where(cand)[0], tol=35):
        if b1 - b0 < 12 or b0 < spans[0][0] or b0 > h - 220:
            continue                      # specks, header area, ADS footer
        if ink[b0:b1 + 1].sum() < 150:
            continue
        bands.append((b0, b1))

    items: list[tuple[int, str, tuple, int]] = []
    items.append((0, "header", (0, 0, w, max(10, spans[0][0] - 2)), 0))
    for s0, s1 in spans:
        # split point for the stacked halves: the table rule nearest to the
        # middle, so no star-range number is cut in two
        rule_cols = np.where(rules[s0:min(s0 + 175, s1)].any(axis=0))[0]
        split = x0 + (rule_cols[np.abs(rule_cols - (w // 2 - x0)).argmin()]
                      if len(rule_cols) else w // 2 - x0)
        items.append((s0, "tablehead",
                      (0, max(0, s0 - 4), w, min(h, s0 + 175)), split))
    for b0, b1 in bands:
        # trim the band to its ink so Gemini sees the digits large; when the
        # ink touches the analysis margin the glyphs continue beyond it
        # (left-shifted verso pages clipped "Dpl." to "pl."), so open the
        # crop up to the page edge there
        cols = np.where(text[b0:b1 + 1].any(axis=0))[0]
        bx0 = x0 + cols[0] - 25 if cols[0] >= 25 else 35
        bx1 = x0 + cols[-1] + 25 if cols[-1] < text.shape[1] - 25 else w - 10
        items.append((b0, "text", (bx0, max(0, b0 - 8), bx1, min(h, b1 + 9)), 0))
    items.sort()

    strips = []
    for i, (_, kind, box, split) in enumerate(items):
        crop = img.crop(tuple(2 * v for v in box))
        if kind == "header":
            # big glyphs; halve so Gemini does not downscale the wide strip
            crop = crop.reduce(2)
        elif kind == "tablehead":
            # keep the 600 DPI digits but stack the two halves vertically:
            # Gemini downscales overly wide images and then misreads digits
            cut = 2 * split
            half = max(cut, crop.width - cut)
            stack = Image.new("RGB", (half, crop.height * 2 + 20),
                              (255, 255, 255))
            stack.paste(crop.crop((0, 0, cut, crop.height)), (0, 0))
            stack.paste(crop.crop((cut, 0, crop.width, crop.height)),
                        (0, crop.height + 20))
            crop = stack
        crop = _label(crop, i + 1)
        if debug_tag:
            crop.save(CACHE_DIR / f"debug_{debug_tag}_{i:02d}_{kind}.png")
        strips.append((kind, crop))
    return strips


def _label(crop: Image.Image, idx: int) -> Image.Image:
    """Draw the image index into a white margin on top of the strip, so
    Gemini cannot mix up which content belongs to which image."""
    out = Image.new("RGB", (crop.width, crop.height + 90), (255, 255, 255))
    out.paste(crop, (0, 90))
    ImageDraw.Draw(out).text((12, 8), f"[{idx}]", fill=(0, 0, 0),
                             font=ImageFont.load_default(size=64))
    return out


# -- Gemini ----------------------------------------------------------------------

api_calls = 0


def transcribe(strips: list[tuple[str, Image.Image]],
               scale: float = 1.0) -> list[list[str]]:
    """One API call; returns the lines of each strip (parallel to strips).
    scale != 1 resizes the digit-bearing strips: re-reading the page at a
    slightly different resolution decorrelates the model's misreads, so
    run agreement means much more."""
    global api_calls
    images = [im if kind == "header" or scale == 1.0 else
              im.resize((int(im.width * scale), int(im.height * scale)),
                        Image.LANCZOS)
              for kind, im in strips]
    attempt = net_waits = 0
    while attempt < 4 and net_waits < 60:   # outages tolerated up to ~1 h
        attempt += 1
        time.sleep(CALL_DELAY)
        try:
            api_calls += 1
            raw = client.models.generate_content(
                model=MODEL, contents=[PROMPT, *images]).text
            text = re.sub(r"^```[a-z]*\n?", "", raw.strip()).rstrip("`").strip()
            m = re.search(r"\[.*\]", text, re.S)
            try:
                rows = json.loads(m.group(0) if m else text)
            except json.JSONDecodeError as e:
                print(f"    malformed response ({e}), retrying", file=sys.stderr)
                continue
            out = [[] for _ in strips]
            for r in rows:
                i = int(r.get("image", 0)) - 1
                if 0 <= i < len(strips):
                    out[i] = [str(x) for x in r.get("lines", [])]
            return out
        except Exception as e:
            msg = str(e)
            if "RESOURCE_EXHAUSTED" in msg or "429" in msg:
                if "PerDay" in msg or "per_day" in msg.lower():
                    raise QuotaExhausted(msg[:200])
                print("    429 per-minute, waiting 60s", file=sys.stderr)
                time.sleep(60)
                continue
            if "UNAVAILABLE" in msg or "503" in msg:
                print(f"    503, backing off ({attempt})", file=sys.stderr)
                time.sleep(20)
                continue
            if "Connect" in type(e).__name__ or "unreachable" in msg \
                    or "nodename" in msg or "timed out" in msg.lower():
                # network outage: wait it out without consuming attempts
                attempt -= 1
                net_waits += 1
                print(f"    network error ({msg[:80]}), waiting 60s",
                      file=sys.stderr)
                time.sleep(60)
                continue
            raise
    raise RuntimeError("gave up after repeated API errors")


# -- parsing -----------------------------------------------------------------------

def _new_entry() -> dict:
    return {"first": None, "last": None, "dpl": [], "color": [],
            "dpl_suspect": [], "color_suspect": [], "notes": []}


def _repair(n: int, lo: int, hi: int) -> int | None:
    """A number outside the page star interval usually has one misread
    digit: return the unique single-digit variant inside the interval."""
    s = str(n)
    cands = {int(s[:i] + d + s[i + 1:])
             for i in range(len(s)) for d in "0123456789"
             if d != s[i] and lo <= int(s[:i] + d + s[i + 1:]) <= hi}
    return cands.pop() if len(cands) == 1 else None


def _parse_footnotes(text: str, entry: dict) -> int:
    """Append the numbers after each Dpl/Color keyword to the entry.
    Only the run of digits/commas right after the keyword is taken, so a
    printer mark or page number later in the band cannot leak in.
    "pl"/"olor" cover keywords whose first letter was clipped by the crop.
    Returns the number of footnotes recognized."""
    found = 0
    for m in re.finditer(r"\b(Dpl|Color|pl|olor)\b[\s.:]*(?:N[ºo°]?\.?)?\s*"
                         r"((?:\d+[\s,.]*)+)", text, re.I):
        kind = "dpl" if "pl" in m.group(1).lower() else "color"
        entry[kind].extend(int(n) for n in re.findall(r"\d+", m.group(2)))
        found += 1
    return found


def parse_page(strips: list[tuple[str, Image.Image]],
               lines: list[list[str]], prev_dec: int | None) -> dict[int, dict]:
    """Turn the transcripts into {dec: entry} in page order."""
    header_text = " ".join(lines[0]) if lines else ""
    active = None
    m = re.search(r"ZONA\s*[-—–−]?\s*(\d{1,2})", header_text, re.I)
    if m and int(m.group(1)) in DECS:
        active = int(m.group(1))
    if active is None:
        # some transcripts drop the leading minus of the declination label,
        # so accept a bare "NN°" too (the hour labels carry "h", not "°")
        m = re.search(r"[-—–−]?\s*(\d{2})\s*[°º]", header_text)
        if m and int(m.group(1)) in DECS:
            active = int(m.group(1))
    entries: dict[int, dict] = {}
    if active is None:
        active = prev_dec
        if active is None:
            raise ValueError(f"cannot read declination from header {header_text!r}")
        entries[active] = _new_entry()
        entries[active]["notes"].append(
            f"declination not readable in header {header_text!r}; "
            f"carried over -{active}")
    else:
        entries[active] = _new_entry()

    for (kind, _), lns in zip(strips[1:], lines[1:]):
        text = " ".join(lns)
        if kind == "tablehead":
            src = " ".join(lns)
            ranges = [(int(a), int(b)) for a, b in
                      re.findall(r"(\d{1,5})\s*[-—–−]\s*(\d{1,5})", src)]
            limit = MAX_STAR.get(active, 20000) + 10
            ranges = [(a, b) for a, b in ranges
                      if 1 <= a <= b <= limit and b - a <= 90]
            if not ranges:
                entries[active]["notes"].append(f"unreadable table header {src!r}")
                continue
            # Keep the longest chain of consecutive ranges: stray matches
            # (an hour label "2h-3h" transcribed as "2-3") never chain.
            # Sorted first: the model sometimes reports the two stacked
            # halves of the strip interleaved.
            ranges.sort()
            best, cur = [ranges[0]], [ranges[0]]
            for r in ranges[1:]:
                cur = cur + [r] if r[0] == cur[-1][1] + 1 else [r]
                if len(cur) > len(best):
                    best = cur
            e = entries[active]
            if len(best) < len(ranges):
                e["notes"].append(f"dropped stray ranges: kept {best} "
                                  f"of {ranges}")
            lo, hi = best[0][0], best[-1][1]
            e["first"] = lo if e["first"] is None else min(e["first"], lo)
            e["last"] = hi if e["last"] is None else max(e["last"], hi)
        else:                                   # text band
            m = re.search(r"ZONA\s*[-—–−]?\s*(\d{1,2})", text, re.I)
            if m and int(m.group(1)) in DECS:
                active = int(m.group(1))
                entries.setdefault(active, _new_entry())
                continue
            if not _parse_footnotes(text, entries[active]) \
                    and re.search(r"N[ºo°]\s*\d", text):
                entries[active]["notes"].append(f"unparsed footnote: {text!r}")
    return entries


def validate(entries: dict[int, dict]) -> None:
    """In-place: check numbers against page interval and declination size,
    move numbers that fit another declination on the page."""
    decs = list(entries)
    for d in decs:
        e = entries[d]
        for kind in ("dpl", "color"):
            keep = []
            for n in e[kind]:
                lo, hi = e["first"], e["last"]
                if lo is not None and hi is not None and not lo <= n <= hi:
                    other = next((o for o in decs if o != d
                                  and entries[o]["first"] is not None
                                  and entries[o]["first"] <= n <= entries[o]["last"]),
                                 None)
                    if other is not None:
                        entries[other][kind].append(n)
                        entries[other]["notes"].append(
                            f"{kind} {n} moved here from -{d} (page interval)")
                        continue
                    # repair only against a plausibly complete page interval:
                    # a misread range chain can leave a bogus narrow one
                    fixed = _repair(n, lo, hi) if hi - lo >= 150 else None
                    if fixed is not None:
                        e["notes"].append(f"{kind} {n} outside page interval "
                                          f"{lo}-{hi}: repaired to {fixed}")
                        e[kind + "_suspect"].append(fixed)
                        keep.append(fixed)
                        continue
                    e["notes"].append(
                        f"{kind} {n} outside page interval {lo}-{hi}")
                    e[kind + "_suspect"].append(n)
                if n > MAX_STAR.get(d, 10 ** 6):
                    e["notes"].append(
                        f"{kind} {n} beyond last star {MAX_STAR.get(d)} of -{d}")
                    e[kind + "_suspect"].append(n)
                keep.append(n)
            e[kind] = keep
        for kind in ("dpl", "color"):
            if e[kind] != sorted(e[kind]):
                e["notes"].append(f"{kind} numbers not ascending: {e[kind]}")


def _key(entries: dict[int, dict]) -> str:
    return json.dumps({d: {"first": e["first"], "last": e["last"],
                           "dpl": sorted(set(e["dpl"])),
                           "color": sorted(set(e["color"]))}
                       for d, e in entries.items()}, sort_keys=True)


def merge_runs(runs: list[dict[int, dict]]) -> dict[int, dict]:
    """Majority vote across disagreeing runs; single-vote numbers become
    uncertain."""
    merged: dict[int, dict] = {}
    decs: list[int] = []
    for r in runs:
        for d in r:
            if d not in decs:
                decs.append(d)
    for d in decs:
        present = [r[d] for r in runs if d in r]
        e = _new_entry()
        e["dpl_uncertain"], e["color_uncertain"] = [], []
        for field in ("first", "last"):
            votes = Counter(x[field] for x in present if x[field] is not None)
            if votes:
                e[field] = votes.most_common(1)[0][0]
        for kind in ("dpl", "color"):
            votes = Counter()
            susp = {n for x in present for n in x[kind + "_suspect"]}
            for x in present:
                votes.update(set(x[kind]))
            for n, c in sorted(votes.items()):
                if c >= 2 and n not in susp:
                    e[kind].append(n)
                else:
                    e[kind + "_uncertain"].append(n)
                    if c < 2:
                        e["notes"].append(f"{kind} {n} seen in only 1 of "
                                          f"{len(runs)} runs")
        for x in present:
            for note in x["notes"]:
                if note not in e["notes"]:
                    e["notes"].append(note)
        merged[d] = e
    return merged


# -- per-page pipeline --------------------------------------------------------------

def process_page(page: int, prev_dec: int | None,
                 debug: bool = False) -> list[dict]:
    img = render_page(page)
    strips = cut_strips(img, debug_tag=f"p{page}" if debug else None)
    runs: list[dict[int, dict]] = []
    for run in range(MAX_RUNS):
        lines = transcribe(strips, scale=(1.0, 0.8, 0.9)[run % 3])
        if debug:
            for (kind, _), lns in zip(strips, lines):
                print(f"    [{kind}] {lns}")
        try:
            entries = parse_page(strips, lines, prev_dec)
        except ValueError as e:
            print(f"    run {run+1}: {e}", file=sys.stderr)
            continue
        validate(entries)
        runs.append(entries)
        if len(runs) >= 2 and _key(runs[-1]) == _key(runs[-2]):
            break
    if not runs:
        raise ValueError("no successful run")
    if len(runs) >= 2 and _key(runs[-1]) == _key(runs[-2]):
        final = runs[-1]
        for d, e in final.items():
            for kind in ("dpl", "color"):
                susp = {n for r in runs if d in r
                        for n in r[d][kind + "_suspect"]}
                e[kind + "_uncertain"] = sorted(set(e[kind]) & susp)
                e[kind] = sorted(set(e[kind]) - susp)
    else:
        final = merge_runs(runs)

    out = []
    for d, e in final.items():
        out.append({"page": page, "decl": -d,
                    "first": e["first"], "last": e["last"],
                    "dpl": e["dpl"], "color": e["color"],
                    "dpl_uncertain": e["dpl_uncertain"],
                    "color_uncertain": e["color_uncertain"],
                    "notes": e["notes"]})
    return out


def json_path(page: int) -> Path:
    return OUT_DIR / f"rnao16_page{page}.json"


# -- CLI ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pages", nargs=2, type=int, required=True,
                        metavar=("FIRST", "LAST"),
                        help="printed page range (LAST=0: continue until the "
                             "catalog or the API quota is exhausted)")
    parser.add_argument("--force", action="store_true",
                        help="reprocess pages whose JSON already exists")
    parser.add_argument("--debug", action="store_true",
                        help="save strip crops to the cache dir, print transcripts")
    parser.add_argument("--model", default=None,
                        help="Gemini model override (each model has its own "
                             "free-tier daily quota)")
    args = parser.parse_args()
    if args.model:
        global MODEL
        MODEL = args.model if args.model.startswith("models/") \
            else f"models/{args.model}"

    first, last = args.pages
    if last == 0:
        last = 10 ** 6
    prev_dec: int | None = None
    failures = 0
    try:
        for page in range(first, last + 1):
            if page + PAGE_OFFSET > 655:
                print("End of the PDF reached.")
                break
            path = json_path(page)
            if path.exists() and not args.force:
                data = json.loads(path.read_text())
                prev_dec = -data[-1]["decl"]
                print(f"Page {page}: exists, skipping")
                continue
            print(f"Page {page} ...", flush=True)
            try:
                entries = process_page(page, prev_dec, debug=args.debug)
            except QuotaExhausted:
                raise
            except Exception as e:
                failures += 1
                print(f"  FAILED: {type(e).__name__}: {e}", file=sys.stderr)
                if failures >= 3:
                    print("  3 consecutive failures - stopping", file=sys.stderr)
                    break
                continue
            failures = 0
            path.write_text(json.dumps(entries, indent=1) + "\n")
            prev_dec = -entries[-1]["decl"]
            summary = "; ".join(
                f"dec {e['decl']} [{e['first']}-{e['last']}] "
                f"dpl={e['dpl'] + ['?%d' % n for n in e['dpl_uncertain']]} "
                f"color={e['color'] + ['?%d' % n for n in e['color_uncertain']]}"
                for e in entries)
            flags = "; ".join(n for e in entries for n in e["notes"])
            print(f"  {summary}  ({api_calls} calls)"
                  + (f"  NOTES: {flags}" if flags else ""), flush=True)
    except QuotaExhausted as e:
        print(f"\nQUOTA EXHAUSTED after {api_calls} API calls: {e}")
        sys.exit(2)


if __name__ == "__main__":
    main()
