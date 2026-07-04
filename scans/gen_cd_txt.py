#!/usr/bin/env python3
"""
Aggregate the scans/rnao16_pageNN.json footnote files (written by
extract_cd_footnotes.py) into the cd/dpl_NN.txt and cd/color_NN.txt lists
consumed by gen_tycho2.cpp: one star number per line, ascending, uncertain
numbers marked as "NNNN ?".

Also reports coverage: pages present for each declination, star-interval
continuity between consecutive pages, and numbers beyond the last star of
the declination (cat/cd_vol1.txt).

Usage:
  python3 gen_cd_txt.py --check 22 24          # compare against existing cd files
  python3 gen_cd_txt.py --write 25 31          # write cd/{dpl,color}_NN.txt
  python3 gen_cd_txt.py --write 25 31 --force  # overwrite existing files
"""

import argparse
import json
import re
import sys
from pathlib import Path

SCANS   = Path(__file__).parent
CD_DIR  = Path(__file__).parent.parent / "cd"
CD_VOL1 = Path(__file__).parent.parent / "cat" / "cd_vol1.txt"


def load_max_star() -> dict[int, int]:
    mx: dict[int, int] = {}
    for line in open(CD_VOL1):
        if not line.startswith("CD-"):
            continue
        d, n = int(line[3:5]), int(line[5:10])
        mx[d] = max(mx.get(d, 0), n)
    return mx


def load_json_entries() -> list[dict]:
    entries = []
    for path in SCANS.glob("rnao16_page*.json"):
        m = re.match(r"rnao16_page(\d+)\.json$", path.name)
        if not m:
            continue
        for e in json.loads(path.read_text()):
            e.setdefault("page", int(m.group(1)))
            entries.append(e)
    return sorted(entries, key=lambda e: e["page"])


def collect(entries: list[dict], dec: int, max_star: dict[int, int]):
    """Return ({'dpl': {num: uncertain_bool}, 'color': ...}, warnings)."""
    nums = {"dpl": {}, "color": {}}
    warnings: list[str] = []
    mine = [e for e in entries if e["decl"] == -dec]
    if not mine:
        return nums, [f"no JSON pages found for declination -{dec}"]

    prev = None
    for e in mine:
        if e["first"] is None or e["last"] is None:
            warnings.append(f"page {e['page']}: star interval unknown")
        elif prev is not None and prev["last"] is not None \
                and e["first"] != prev["last"] + 1:
            warnings.append(
                f"pages {prev['page']}->{e['page']}: interval jumps "
                f"{prev['last']} -> {e['first']} (missing/misread pages?)")
        prev = e
        for kind in ("dpl", "color"):
            for n in e[kind]:
                nums[kind][n] = nums[kind].get(n, False)   # certain wins
            for n in e.get(kind + "_uncertain", []):
                nums[kind].setdefault(n, True)
        for note in e.get("notes", []):
            warnings.append(f"page {e['page']}: {note}")

    if mine[0]["first"] != 1:
        warnings.append(f"first page {mine[0]['page']} starts at star "
                        f"{mine[0]['first']}, not 1")
    want = max_star.get(dec)
    if want and mine[-1]["last"] != want:
        warnings.append(f"last page {mine[-1]['page']} ends at star "
                        f"{mine[-1]['last']}, catalog has {want}")
    for kind in ("dpl", "color"):
        for n in nums[kind]:
            if want and n > want:
                warnings.append(f"{kind} {n} beyond last star {want}")
    return nums, warnings


def read_known(path: Path) -> dict[int, bool]:
    """Existing cd file -> {num: has_question_mark}."""
    known: dict[int, bool] = {}
    for line in path.read_text().splitlines():
        m = re.match(r"\s*(\d+)\s*(\?)?\s*$", line)
        if m:
            known[int(m.group(1))] = bool(m.group(2))
    return known


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", nargs=2, type=int, metavar=("DEC1", "DEC2"),
                      help="compare against existing cd/ files (no writes)")
    mode.add_argument("--write", nargs=2, type=int, metavar=("DEC1", "DEC2"),
                      help="write cd/{dpl,color}_NN.txt")
    parser.add_argument("--force", action="store_true",
                        help="overwrite existing cd/ files")
    args = parser.parse_args()

    max_star = load_max_star()
    entries = load_json_entries()
    d1, d2 = args.check or args.write
    all_ok = True

    for dec in range(d1, d2 + 1):
        nums, warnings = collect(entries, dec, max_star)
        print(f"=== declination -{dec} ===")
        for wtext in warnings:
            print(f"  warn: {wtext}")
        for kind in ("dpl", "color"):
            got = nums[kind]
            if args.check:
                path = CD_DIR / f"{kind}_{dec}.txt"
                if not path.exists():
                    print(f"  {kind}: no existing file {path} to compare")
                    continue
                known = read_known(path)
                missing = sorted(set(known) - set(got))
                extra = sorted(set(got) - set(known))
                qdiff = sorted(n for n in set(known) & set(got)
                               if known[n] != got[n])
                status = "OK" if not missing and not extra else "MISMATCH"
                print(f"  {kind}: {len(got)} extracted vs {len(known)} known "
                      f"[{status}]")
                if missing:
                    print(f"    missing (known, not extracted): {missing}")
                if extra:
                    print(f"    extra   (extracted, not known): {extra}")
                if qdiff:
                    print(f"    '?' flag differs: {qdiff}")
                if missing or extra:
                    all_ok = False
            else:
                path = CD_DIR / f"{kind}_{dec}.txt"
                if path.exists() and not args.force:
                    print(f"  {kind}: {path} exists, use --force to overwrite")
                    all_ok = False
                    continue
                CD_DIR.mkdir(exist_ok=True)
                with open(path, "w") as f:
                    for n in sorted(got):
                        f.write(f"{n} ?\n" if got[n] else f"{n}\n")
                print(f"  {kind}: written {path} ({len(got)} entries, "
                      f"{sum(got.values())} uncertain)")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
