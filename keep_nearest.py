#!/usr/bin/env python3
"""
Keep only the nearest Tycho-2 match for each index2 identifier.
Usage: python3 keep_nearest.py <input.csv>
Output: <input.csv>.2
"""

import sys
import csv

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
            if len(row) < 3:
                continue
            index1, index2, dist_str = row[0], row[1], row[2]
            rows.append((index1, index2, dist_str, float(dist_str)))

    # For each index2, keep the row with the smallest distance
    best = {}
    for index1, index2, dist_str, dist in rows:
        if index2 not in best or dist < best[index2][3]:
            best[index2] = (index1, index2, dist_str, dist)

    # Write output preserving original order, skipping duplicates
    seen = set()
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for index1, index2, dist_str, dist in rows:
            if index2 in seen:
                continue
            if best[index2][0] == index1:
                writer.writerow([index1, index2, dist_str])
                seen.add(index2)

    print(f"Written {len(best)} rows to {output_path}")

if __name__ == "__main__":
    main()
