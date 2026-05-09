#!/usr/bin/env python3
"""
Maximum-likelihood cross-identification between two star catalogs.

Reads two CSVs in the format produced for likelihood/cat1875 (header
`name,x,y,z,mag` with x,y,z as unit-vector coordinates and mag as visual
magnitude with sentinel >99 for "unknown").

Builds a bipartite graph where each edge (a,b) carries a negative-log-
likelihood cost combining angular distance and magnitude difference,
prunes by hard thresholds, decomposes the graph into connected components,
and solves each component as a minimum-cost assignment (brute force for
small components, Hungarian via scipy for larger ones).

Output schema mirrors results/cross/*.csv:
    index1,index2,mag,dist        (mag column omitted for any row whose
                                   index1 has unknown magnitude)
"""

import sys
from itertools import permutations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

SIGMA_POS_ARCSEC = 30.0
SIGMA_MAG = 0.5
MAX_THETA_ARCSEC = 300.0
MAX_DMAG = 3.0
MISSING_DMAG = 2.0
MAG_MISSING_SENTINEL = 99.0
MAG_ZERO_EPSILON = 0.0001  # |mag| < eps is treated as "no magnitude reported"

ARCSEC_PER_RAD = 180.0 * 3600.0 / np.pi
BRUTE_FORCE_LIMIT = 5

POS_COST_FACTOR = 1.0 / (2.0 * SIGMA_POS_ARCSEC * SIGMA_POS_ARCSEC)
MAG_COST_FACTOR = 1.0 / (2.0 * SIGMA_MAG * SIGMA_MAG)
MISSING_MAG_COST = MISSING_DMAG * MISSING_DMAG * MAG_COST_FACTOR


def mag_unknown(m):
    return m > MAG_MISSING_SENTINEL or abs(m) < MAG_ZERO_EPSILON


def load_catalog(path):
    df = pd.read_csv(path)
    names = df["name"].to_numpy()
    xyz = df[["x", "y", "z"]].to_numpy(dtype=np.float64)
    # Inputs are stored at 8-decimal precision; renormalize so |v|=1 exactly,
    # otherwise arccos(v·v) yields ~17 arcsec of phantom self-distance.
    norms = np.sqrt((xyz * xyz).sum(axis=1, keepdims=True))
    xyz /= norms
    mag = df["mag"].to_numpy(dtype=np.float64)
    return names, xyz, mag


def build_edges(xyz_a, mag_a, xyz_b, mag_b):
    """
    Returns three parallel numpy arrays (i, j, cost) and a parallel
    `theta` array (arcsec) for surviving edges.
    """
    chord_thresh = 2.0 * np.sin(np.deg2rad(MAX_THETA_ARCSEC / 3600.0) / 2.0)
    tree = cKDTree(xyz_b)
    candidates = tree.query_ball_point(xyz_a, r=chord_thresh)

    i_list = []
    j_list = []
    cost_list = []
    theta_list = []

    for i, js in enumerate(candidates):
        if not js:
            continue
        x_a, y_a, z_a = xyz_a[i]
        m_a = mag_a[i]
        for j in js:
            dx = x_a - xyz_b[j, 0]
            dy = y_a - xyz_b[j, 1]
            dz = z_a - xyz_b[j, 2]
            chord = np.sqrt(dx * dx + dy * dy + dz * dz)
            # 2*arcsin(c/2) is numerically stable for small angles, unlike arccos(dot).
            theta_arcsec = 2.0 * np.arcsin(min(chord * 0.5, 1.0)) * ARCSEC_PER_RAD
            if theta_arcsec > MAX_THETA_ARCSEC:
                continue
            m_b = mag_b[j]
            both_known = not mag_unknown(m_a) and not mag_unknown(m_b)
            if both_known and abs(m_a - m_b) > MAX_DMAG:
                continue
            pos = theta_arcsec * theta_arcsec * POS_COST_FACTOR
            if both_known:
                dmag = m_a - m_b
                cost = pos + dmag * dmag * MAG_COST_FACTOR
            else:
                cost = pos + MISSING_MAG_COST
            i_list.append(i)
            j_list.append(j)
            cost_list.append(cost)
            theta_list.append(theta_arcsec)

    return (
        np.asarray(i_list, dtype=np.int64),
        np.asarray(j_list, dtype=np.int64),
        np.asarray(cost_list, dtype=np.float64),
        np.asarray(theta_list, dtype=np.float64),
    )


def find_components(n_a, n_b, edge_i, edge_j):
    """
    Returns an array of length (n_a + n_b) giving the component label for
    each vertex (A vertices first, then B). Isolated vertices get their
    own labels.
    """
    n_total = n_a + n_b
    if edge_i.size == 0:
        return np.arange(n_total, dtype=np.int64)
    rows = np.concatenate([edge_i, n_a + edge_j])
    cols = np.concatenate([n_a + edge_j, edge_i])
    data = np.ones(rows.size, dtype=np.int8)
    graph = csr_matrix((data, (rows, cols)), shape=(n_total, n_total))
    _, labels = connected_components(graph, directed=False)
    return labels


def solve_brute_force(a_idx, b_idx, cost_lookup):
    """
    Enumerate all matchings of a_idx onto subsets of b_idx (or vice versa
    if |b| < |a|). Returns list of (a, b) selected edges.
    """
    if len(a_idx) <= len(b_idx):
        small, large = a_idx, b_idx
        small_is_a = True
    else:
        small, large = b_idx, a_idx
        small_is_a = False

    best_cost = np.inf
    best_match = []
    for perm in permutations(large, len(small)):
        total = 0.0
        edges = []
        feasible = True
        for s, l in zip(small, perm):
            a, b = (s, l) if small_is_a else (l, s)
            c = cost_lookup.get((a, b))
            if c is None:
                feasible = False
                break
            total += c
            edges.append((a, b))
        if feasible and total < best_cost:
            best_cost = total
            best_match = edges
    return best_match


def solve_hungarian(a_idx, b_idx, cost_lookup):
    """
    Build a |A|x|B| cost matrix with a sentinel for absent edges, run the
    Hungarian algorithm, and drop assignments that landed on sentinels.
    """
    sentinel = (
        MAX_THETA_ARCSEC * MAX_THETA_ARCSEC * POS_COST_FACTOR
        + MAX_DMAG * MAX_DMAG * MAG_COST_FACTOR
        + 1.0
    )
    n_a, n_b = len(a_idx), len(b_idx)
    cost = np.full((n_a, n_b), sentinel, dtype=np.float64)
    a_pos = {a: p for p, a in enumerate(a_idx)}
    b_pos = {b: p for p, b in enumerate(b_idx)}
    for (a, b), c in cost_lookup.items():
        if a in a_pos and b in b_pos:
            cost[a_pos[a], b_pos[b]] = c
    row_ind, col_ind = linear_sum_assignment(cost)
    edges = []
    for r, c in zip(row_ind, col_ind):
        if cost[r, c] >= sentinel:
            continue
        edges.append((a_idx[r], b_idx[c]))
    return edges


def solve_components(n_a, n_b, edge_i, edge_j, edge_cost_arr, labels):
    """
    Returns dict (a -> b) of selected matches.
    """
    cost_lookup = {}
    for i, j, c in zip(edge_i, edge_j, edge_cost_arr):
        cost_lookup[(int(i), int(j))] = float(c)

    by_component = {}
    for v in range(n_a + n_b):
        comp = int(labels[v])
        if comp not in by_component:
            by_component[comp] = ([], [])
        if v < n_a:
            by_component[comp][0].append(v)
        else:
            by_component[comp][1].append(v - n_a)

    matches = {}
    sizes = []
    n_brute = 0
    n_hungarian = 0
    for a_idx, b_idx in by_component.values():
        if not a_idx or not b_idx:
            continue
        sizes.append((len(a_idx), len(b_idx)))
        if min(len(a_idx), len(b_idx)) <= BRUTE_FORCE_LIMIT:
            edges = solve_brute_force(a_idx, b_idx, cost_lookup)
            n_brute += 1
        else:
            edges = solve_hungarian(a_idx, b_idx, cost_lookup)
            n_hungarian += 1
        for a, b in edges:
            matches[a] = b

    print(
        f"components: {len(sizes)} non-trivial, "
        f"brute={n_brute}, hungarian={n_hungarian}",
        file=sys.stderr,
    )
    if sizes:
        sizes_arr = np.array(sizes)
        print(
            f"component sizes (|A|,|B|): "
            f"max=({sizes_arr[:, 0].max()},{sizes_arr[:, 1].max()}), "
            f"mean=({sizes_arr[:, 0].mean():.2f},{sizes_arr[:, 1].mean():.2f})",
            file=sys.stderr,
        )
        flat = sizes_arr.sum(axis=1)
        unique, counts = np.unique(flat, return_counts=True)
        hist = ", ".join(f"{u}:{c}" for u, c in zip(unique, counts))
        print(f"|A|+|B| histogram: {hist}", file=sys.stderr)
    return matches


def write_output(path, names_a, mag_a, names_b, matches, edge_lookup_theta):
    with open(path, "w") as f:
        f.write("index1,index2,mag,dist\n")
        for a in sorted(matches):
            b = matches[a]
            theta = edge_lookup_theta[(a, b)]
            m = mag_a[a]
            if mag_unknown(m):
                f.write(f"{names_a[a]},{names_b[b]},{theta:.2f}\n")
            else:
                f.write(f"{names_a[a]},{names_b[b]},{m:.1f},{theta:.2f}\n")


def main():
    if len(sys.argv) != 4:
        print("usage: cross_likelihood.py A.csv B.csv output.csv", file=sys.stderr)
        sys.exit(2)
    path_a, path_b, path_out = sys.argv[1], sys.argv[2], sys.argv[3]

    print(f"loading {path_a}...", file=sys.stderr)
    names_a, xyz_a, mag_a = load_catalog(path_a)
    print(f"loading {path_b}...", file=sys.stderr)
    names_b, xyz_b, mag_b = load_catalog(path_b)
    print(f"|A|={len(names_a)}, |B|={len(names_b)}", file=sys.stderr)

    print("building edges...", file=sys.stderr)
    edge_i, edge_j, edge_cost_arr, edge_theta = build_edges(xyz_a, mag_a, xyz_b, mag_b)
    print(f"edges within thresholds: {edge_i.size}", file=sys.stderr)

    edge_lookup_theta = {
        (int(i), int(j)): float(t) for i, j, t in zip(edge_i, edge_j, edge_theta)
    }

    print("finding connected components...", file=sys.stderr)
    labels = find_components(len(names_a), len(names_b), edge_i, edge_j)

    print("solving per-component matching...", file=sys.stderr)
    matches = solve_components(
        len(names_a), len(names_b), edge_i, edge_j, edge_cost_arr, labels
    )
    print(f"matched pairs: {len(matches)}", file=sys.stderr)

    print(f"writing {path_out}...", file=sys.stderr)
    write_output(path_out, names_a, mag_a, names_b, matches, edge_lookup_theta)
    print("done.", file=sys.stderr)


if __name__ == "__main__":
    main()
