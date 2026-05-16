#!/usr/bin/env python3
"""
Keep only the best Tycho-2 match for each index2 identifier.
Input columns: index1,index2,mag,dist
Usage: python3 keep_nearest.py <input.csv>
Output: <input.csv>.2
"""

import sys
import csv

def better(current, candidate):
    """Pairwise pick: brightest when both mags are meaningful and differ; else nearest."""
    _, _, _, _, cur_mag, cur_dist = current
    _, _, _, _, cand_mag, cand_dist = candidate
    if cur_mag != 0.0 and cand_mag != 0.0 and cur_mag != cand_mag:
        return candidate if cand_mag < cur_mag else current
    return candidate if cand_dist < cur_dist else current

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.csv>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = input_path + ".2"

    rows = []
    with open(input_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) < 4:
                continue
            index1, index2, mag_str, dist_str = row[0], row[1], row[2], row[3]
            rows.append((index1, index2, mag_str, dist_str, float(mag_str), float(dist_str)))

    # For each index2, reduce candidates pairwise to a single best
    best = {}
    for row in rows:
        index2 = row[1]
        if index2 not in best:
            best[index2] = row
        else:
            best[index2] = better(best[index2], row)

    # Write output preserving original order, one row per index2
    seen = set()
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            index2 = row[1]
            if index2 in seen:
                continue
            if best[index2] is row:
                writer.writerow([row[0], row[1], row[2], row[3]])
                seen.add(index2)

    print(f"Written {len(best)} rows to {output_path}")

if __name__ == "__main__":
    main()
