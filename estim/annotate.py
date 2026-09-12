#!/usr/bin/env python3
"""Overlay measure.py's V magnitudes on a field PNG.

    python3 annotate.py RCen.csv --image RCen.png -o RCen.mag.png

Labels drop the decimal point (V=5.67 -> "567"), the standard finder-chart
convention so a label is never misread as a star. Never overwrites the
existing astrometry.net *.ann.png files -- always write to a new path.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

FONT_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/Library/Fonts/Arial.ttf",
]
FONT_SIZE = 26

COLOR_NORMAL = (80, 255, 80)
COLOR_SAT = (255, 60, 60)
COLOR_EXTRAP = (255, 165, 0)
COLOR_TEXT = (255, 255, 0)


def load_font(size=FONT_SIZE):
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def read_rows(csv_path):
    with open(csv_path, newline="") as fh:
        return list(csv.DictReader(fh))


def find_image(csv_path, image_arg):
    if image_arg:
        return Path(image_arg)
    stem = Path(csv_path).stem
    folder = Path(csv_path).parent
    for ext in (".png", ".jpg"):
        cand = folder / f"{stem}{ext}"
        if cand.exists():
            return cand
    return None


def render_from_green_fits(fits_path, naxis1, naxis2):
    from astropy.io import fits
    from astropy.visualization import AsinhStretch, MinMaxInterval

    data = fits.getdata(fits_path).astype(float)
    if data.ndim == 3:
        data = data[1] if data.shape[0] == 3 else data[0]
    stretch = AsinhStretch() + MinMaxInterval()
    scaled = stretch(data / max(data.max(), 1.0))
    img8 = (np.clip(scaled, 0, 1) * 255).astype(np.uint8)
    return Image.fromarray(img8, mode="L").convert("RGB")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="measure.py output CSV")
    ap.add_argument("--image", help="background PNG/JPG (default: <csv-stem>.png next to the CSV)")
    ap.add_argument("--green-fits", help="fallback: render this green-plane FITS if no PNG/JPG found")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--max-mag", type=float, default=10.5)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only-matched", "--only-tycho2", dest="only_tycho2",
                     action="store_true",
                     help="label only stars matched to the reference catalogue")
    ap.add_argument("--show", choices=["v", "bv", "both"], default="v")
    args = ap.parse_args()

    rows = read_rows(args.csv)
    for r in rows:
        r["_V"] = float(r["V"]) if r["V"] else None
        r["_BV"] = float(r["BV"]) if r["BV"] else None
        r["_x"] = float(r["x"])
        r["_y"] = float(r["y"])
        r["_flags"] = r["flags"].split("|") if r["flags"] else []
        # 'ref_id' since schema 3 (any reference catalogue); 'tyc2' in older CSVs
        r["_tyc2"] = r.get("ref_id", r.get("tyc2", ""))

    sel = [r for r in rows if r["_V"] is not None and r["_V"] <= args.max_mag]
    if args.only_tycho2:
        sel = [r for r in sel if r["_tyc2"]]
    sel.sort(key=lambda r: r["_V"])
    if args.limit:
        sel = sel[: args.limit]

    image_path = find_image(args.csv, args.image)
    if image_path is not None:
        img = Image.open(image_path).convert("RGB")
    elif args.green_fits:
        naxis1 = naxis2 = None
        img = render_from_green_fits(args.green_fits, naxis1, naxis2)
    else:
        print("error: no PNG/JPG found and --green-fits not given", file=sys.stderr)
        sys.exit(1)

    draw = ImageDraw.Draw(img)
    font = load_font()

    for r in sel:
        # The Dwarf's PNG preview is rendered directly from the pixel array
        # in row order (row 0 = row 0 of the FITS data), not flipped to the
        # bottom-up astronomical convention -- so this is a 1-based -> 0-based
        # conversion only, no vertical mirror.
        x = r["_x"] - 1.0
        y_png = r["_y"] - 1.0

        if "SAT" in r["_flags"]:
            color = COLOR_SAT
        elif "BV_EXTRAP" in r["_flags"]:
            color = COLOR_EXTRAP
        else:
            color = COLOR_NORMAL
        marker = "s" if r["_tyc2"] else "o"  # square = Tycho-2 matched, circle = unmatched

        rad = 10
        if marker == "o":
            draw.ellipse([x - rad, y_png - rad, x + rad, y_png + rad], outline=color, width=2)
        else:
            draw.rectangle([x - rad, y_png - rad, x + rad, y_png + rad], outline=color, width=2)

        if args.show == "v":
            label = f"{r['_V']:.2f}".replace(".", "")
        elif args.show == "bv" and r["_BV"] is not None:
            label = f"{r['_BV']:.2f}".replace(".", "").replace("-", "m")
        elif args.show == "both":
            v_lbl = f"{r['_V']:.2f}".replace(".", "")
            bv_lbl = f"{r['_BV']:.2f}".replace(".", "").replace("-", "m") if r["_BV"] is not None else "??"
            label = f"{v_lbl}/{bv_lbl}"
        else:
            label = f"{r['_V']:.2f}".replace(".", "")

        draw.text((x + rad + 3, y_png - rad - 2), label, fill=COLOR_TEXT, font=font)

    img.save(args.out)
    print(f"wrote {args.out}: {len(sel)} labels (max_mag={args.max_mag})", file=sys.stderr)


if __name__ == "__main__":
    main()
