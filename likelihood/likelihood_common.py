"""Shared catalog I/O, calibration, evidence and matching functions. See CD_CROSS.md."""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from dataclasses import dataclass, asdict
from collections import defaultdict
import json
import math
from pathlib import Path
import pandas as pd
from scipy.optimize import least_squares, linear_sum_assignment, minimize_scalar
from scipy.stats import t
from scipy.ndimage import gaussian_filter1d
import re
import warnings
from astropy.coordinates import SkyCoord, FK5, FK4
from astropy.time import Time
from astropy import units as u
from erfa import ErfaWarning
import argparse
import hashlib
import sys
from scipy.optimize import linear_sum_assignment
from collections import defaultdict, Counter
from dataclasses import dataclass
from itertools import combinations
import time
from scipy.special import logsumexp, ndtr
from scipy.stats import t, truncnorm
from scipy.optimize import minimize_scalar
import pickle
import gzip
import importlib.metadata
import copy


# LIKELIHOOD COMMON
ARCSEC_PER_RAD = 180.0 * 3600.0 / np.pi

def normalize(xyz):
    xyz = np.asarray(xyz, dtype=float)
    n = np.linalg.norm(xyz, axis=-1, keepdims=True)
    if not np.isfinite(xyz).all() or np.any(n == 0):
        raise ValueError('Coordinates must be finite, nonzero vectors')
    return xyz / n

def separation(a, b):
    return 2 * np.arcsin(np.clip(np.linalg.norm(np.asarray(a) - b, axis=-1) / 2, 0, 1)) * ARCSEC_PER_RAD

def neighbors(a, b, radius):
    return cKDTree(b).query_ball_point(a, 2 * np.sin(radius / (2 * ARCSEC_PER_RAD)))

def find_components(n_a, n_b, edge_i, edge_j):
    edge_i, edge_j = (np.asarray(edge_i, dtype=int), np.asarray(edge_j, dtype=int))
    if not len(edge_i):
        return np.arange(n_a + n_b)
    rows = np.concatenate([edge_i, n_a + edge_j])
    cols = np.concatenate([n_a + edge_j, edge_i])
    graph = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_a + n_b, n_a + n_b))
    return connected_components(graph, directed=False)[1]


# SINGLE LIKELIHOOD

@dataclass
class Settings:
    radius: float = 300.0
    anchor_radius: float = 120.0
    systematic: bool = True
    photometry: bool = True
    broad_fraction: float = 0.03
    broad_scale_floor: float = 60.0
    density_neighbors: int = 64
    min_anchors: int = 80
    global_margin_limit: int = 50
    q: float | None = None
    zero_missing: bool = True
    variable_sigma: float | None = None

def load_catalog(path, zero_missing=True):
    df = pd.read_csv(path)
    if df.name.isna().any() or df.name.duplicated().any():
        raise ValueError('Identifiers must be nonempty and unique')
    v = normalize(df[['x', 'y', 'z']].to_numpy())
    m = pd.to_numeric(df.mag, errors='coerce').to_numpy(dtype=float)
    m[(m >= 99) | (np.abs(m) < 0.0001 if zero_missing else np.zeros(len(m), bool))] = np.nan
    return (df.name.to_numpy(), v, m)

def load_variability(path, names, default_sigma=None, sidecar=None):
    """Optional variable flag and per-epoch variability STD (mag), not amplitude.

    A sidecar CSV with name[,variability_sigma] marks its listed sources variable.
    Unknown dispersion is NaN and disables photometry for those pairs.
    """
    data = pd.read_csv(path).set_index('name')
    variable = np.zeros(len(names), bool)
    sigma = np.zeros(len(names))
    if 'variable' in data:
        tokens = data.loc[names, 'variable'].fillna('').astype(str).str.lower().str.strip()
        allowed = {'', '0', '0.0', 'false', 'no', '1', '1.0', 'true', 'yes', 'v', 'var', 'variable'}
        if not set(tokens) <= allowed:
            raise ValueError('Unrecognized variable flag')
        variable = tokens.isin(['1', '1.0', 'true', 'yes', 'v', 'var', 'variable']).to_numpy()
    if 'variability_sigma' in data:
        values = pd.to_numeric(data.loc[names, 'variability_sigma'], errors='coerce').to_numpy()
        if np.any(np.isfinite(values) & (values < 0)):
            raise ValueError('Variability sigma must be nonnegative')
        variable |= np.isfinite(values) & (values > 0)
    else:
        values = np.full(len(names), np.nan)
    sigma[variable] = values[variable]
    if sidecar is not None:
        extra = pd.read_csv(sidecar)
        if extra.name.duplicated().any() or not set(extra.name) <= set(names):
            raise ValueError('Duplicate or unknown sidecar identifiers')
        positions = {name: i for i, name in enumerate(names)}
        for row in extra.itertuples():
            i = positions[row.name]
            variable[i] = True
            value = getattr(row, 'variability_sigma', np.nan)
            if np.isfinite(value) and value < 0:
                raise ValueError('Variability sigma must be nonnegative')
            sigma[i] = value
    if default_sigma is not None:
        if not np.isfinite(default_sigma) or default_sigma < 0:
            raise ValueError('Invalid default variability sigma')
        sigma[variable & ~np.isfinite(sigma)] = default_sigma
    return (variable, sigma)

def photometric_scale(base, sigma_a, sigma_b, slope):
    return np.sqrt(base ** 2 + (np.asarray(sigma_a) ** 2 + slope ** 2 * np.asarray(sigma_b) ** 2) / 3)

def basis(v):
    ra = np.arctan2(v[:, 1], v[:, 0])
    dec = np.arcsin(np.clip(v[:, 2], -1, 1))
    east = np.column_stack((-np.sin(ra), np.cos(ra), np.zeros(len(v))))
    north = np.column_stack((-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)))
    features = np.column_stack((np.ones(len(v)), np.sin(ra), np.cos(ra), np.sin(2 * ra), np.cos(2 * ra), np.sin(dec)))
    return (east, north, features)

def residuals(a, b):
    e, n, _ = basis(a)
    dot = np.einsum('ij,ij->i', a, b)
    return np.column_stack((np.arctan2(np.einsum('ij,ij->i', b, e), dot), np.arctan2(np.einsum('ij,ij->i', b, n), dot))) * ARCSEC_PER_RAD

def correct(v, coeff):
    e, n, x = basis(v)
    shift = x @ np.asarray(coeff)
    w = normalize(v + (shift[:, 0, None] * e + shift[:, 1, None] * n) / ARCSEC_PER_RAD)
    return (w, np.linalg.norm(shift, axis=1))

def robust_linear(x, y):
    coef = np.linalg.lstsq(x, y, rcond=None)[0]
    for _ in range(10):
        r = y - x @ coef
        scale = max(0.2, 1.4826 * np.median(np.abs(r - np.median(r))))
        weights = np.minimum(1, 1.5 * scale / np.maximum(np.abs(r), 1e-12))
        coef = np.linalg.solve(x.T * weights @ x + np.eye(x.shape[1]) * 1e-06, x.T * weights @ y)
    return coef

def radial_cdf(r, s):
    return 1 - (1 + np.asarray(r) ** 2 / (4 * s * s)) ** (-2)

def spatial_pdf(r, s, broad_fraction=0.03, broad_scale_floor=60.0):
    broad = max(6 * s, broad_scale_floor)

    def density(scale):
        return (1 + np.asarray(r) ** 2 / (4 * scale * scale)) ** (-3) / (2 * np.pi * scale * scale)
    return (1 - broad_fraction) * density(s) + broad_fraction * density(broad)

def spatial_cdf(r, s, settings):
    b = settings.broad_fraction
    return (1 - b) * radial_cdf(r, s) + b * radial_cdf(r, max(6 * s, settings.broad_scale_floor))

def calibrate_single(av, am, bv, bm, settings, variable_a=None, variable_b=None):
    if len(av) < 3 or len(bv) < 2:
        raise ValueError('Calibration requires at least 3 old and 2 modern records')
    d, ix = cKDTree(bv).query(av, k=2)
    angle = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * ARCSEC_PER_RAD
    back = cKDTree(av).query(bv)[1]
    good = (angle[:, 0] < settings.anchor_radius) & (angle[:, 1] > np.maximum(3 * angle[:, 0], 30)) & (back[ix[:, 0]] == np.arange(len(av)))
    if variable_a is not None:
        good &= ~variable_a
    if variable_b is not None:
        good &= ~variable_b[ix[:, 0]]
    ids = np.flatnonzero(good)
    bi = ix[ids, 0]
    if len(ids) < settings.min_anchors:
        raise ValueError(f'Only {len(ids)} reliable anchors; supply a calibrated --model-in or improve coverage')
    _, _, x = basis(av[ids])
    res = residuals(av[ids], bv[bi])
    ra = np.rad2deg(np.arctan2(av[ids, 1], av[ids, 0])) % 360
    valid = np.floor(ra / 5).astype(int) % 5 == 0
    split = 'RA 5-degree sectors, sector mod 5 == 0'
    if min(valid.sum(), (~valid).sum()) < max(20, len(ids) * 0.1):
        valid = ids % 5 == 0
        split = 'record-index mod 5 == 0 (small footprint)'
    train = ~valid
    coef = np.zeros((6, 2))
    trial = np.column_stack([robust_linear(x[train], res[train, j]) for j in range(2)])
    original = np.linalg.norm(res[valid], axis=1)
    after = np.linalg.norm(res[valid] - x[valid] @ trial, axis=1)
    accepted = settings.systematic and np.median(after) < 0.98 * np.median(original)
    if accepted:
        coef = trial
    rr = np.linalg.norm(res[train] - x[train] @ coef, axis=1)
    sigma = max(0.5, float(np.median(rr) / np.sqrt(4 * (np.sqrt(2) - 1))))
    photo = dict(enabled=False, reason='insufficient or uninformative photometry')
    known = np.isfinite(am[ids]) & np.isfinite(bm[bi])
    pt = train & known
    pv = valid & known
    if settings.photometry and pt.sum() >= 50 and (pv.sum() >= 20):
        xx = np.column_stack((np.ones(len(ids)), np.nan_to_num(bm[bi], nan=9) - 9))
        cc = robust_linear(xx[pt], am[ids][pt])
        cc[1] = np.clip(cc[1], 0.1, 2)
        cc[0] = np.median(am[ids][pt] - cc[1] * (bm[bi][pt] - 9))
        r = am[ids][pt] - xx[pt] @ cc
        sc = max(0.35, float(1.4826 * np.median(np.abs(r - np.median(r)))))
        valr = am[ids][pv] - xx[pv] @ cc
        null = am[ids][pv] - np.median(am[ids][pt])
        photo = dict(enabled=bool(np.mean(np.minimum(np.abs(valr), 2)) < 0.95 * np.mean(np.minimum(np.abs(null), 2))), coefficients=cc.tolist(), scale=sc, n_train=int(pt.sum()), n_validation=int(pv.sum()), validation_mae=float(np.median(np.abs(valr))), validation_null_mae=float(np.median(np.abs(null))), validation_clipped_mean_absolute_error=float(np.mean(np.minimum(np.abs(valr), 2))), validation_null_clipped_mean_absolute_error=float(np.mean(np.minimum(np.abs(null), 2))), distribution='Student t, df=3; scale floor .35 mag')
    values = am[ids[train]]
    values = values[np.isfinite(values)]
    edges = np.arange(-5, 25.25, 0.25)
    hist = np.histogram(values, edges)[0].astype(float) + 0.1
    hist = gaussian_filter1d(hist, 1.5)
    hist /= hist.sum() * 0.25
    valres = np.linalg.norm(res[valid] - x[valid] @ coef, axis=1)
    return dict(sigma=sigma, systematic_coefficients=coef.tolist(), systematic_accepted=bool(accepted), systematic_validation_before_median=float(np.median(original)), systematic_validation_after_median=float(np.median(valres)), systematic_validation_trial_median=float(np.median(after)), residual_validation_quantiles={str(q): float(np.quantile(valres, q)) for q in [0.5, 0.9, 0.99]}, anchor_count=len(ids), train_count=int(train.sum()), validation_count=int(valid.sum()), split=split, anchor_selection='mutual nearest; r<anchor_radius; second>max(3*r,30 arcsec)', selection_warning='Anchor truncation/isolation bias; extreme historical errors are not calibrated from these anchors', photometry=photo, mag_centers=((edges[:-1] + edges[1:]) / 2).tolist(), mag_density=hist.tolist())

def evidence(av, am, bv, bm, model, settings, variability_a=None, variability_b=None):
    corrected, shifts = correct(av, model['systematic_coefficients'])
    tree = cKDTree(bv)
    lists = tree.query_ball_point(corrected, 2 * np.sin(settings.radius / (2 * ARCSEC_PER_RAD)))
    ei = np.repeat(np.arange(len(av)), [len(x) for x in lists])
    ej = np.array([b for js in lists for b in sorted(js)], dtype=int)
    r = separation(corrected[ei], bv[ej])
    raw = separation(av[ei], bv[ej])
    k = min(settings.density_neighbors, len(bv))
    ds = tree.query(corrected, k=k)[0]
    dk = ds[:, -1] if k > 1 else ds
    rk = 2 * np.arcsin(np.clip(dk / 2, 0, 1)) * ARCSEC_PER_RAD
    rho = max(k - 1, 1) / (np.pi * np.maximum(rk, 1) ** 2)
    lr = spatial_pdf(r, model['sigma'], settings.broad_fraction, settings.broad_scale_floor) / rho[ei]
    logphot = np.zeros(len(ei))
    photo = model['photometry']
    if photo['enabled'] and settings.photometry:
        va = np.zeros(len(av)) if variability_a is None else np.asarray(variability_a)
        vb = np.zeros(len(bv)) if variability_b is None else np.asarray(variability_b)
        known = np.isfinite(am[ei]) & np.isfinite(bm[ej]) & np.isfinite(va[ei]) & np.isfinite(vb[ej])
        c = photo['coefficients']
        sc = photometric_scale(photo['scale'], va[ei[known]], vb[ej[known]], c[1])
        predicted = c[0] + c[1] * (bm[ej[known]] - 9)
        numerator = t.pdf((am[ei[known]] - predicted) / sc, df=3) / sc
        denominator = np.interp(am[ei[known]], model['mag_centers'], model['mag_density'], left=1e-05, right=1e-05)
        logphot[known] = np.clip(np.log(numerator / denominator), -math.log(20), math.log(20))
        lr *= np.exp(logphot)
    sums = np.bincount(ei, weights=lr, minlength=len(av))
    F = float(spatial_cdf(settings.radius, model['sigma'], settings))
    if settings.q is not None:
        q = settings.q
    elif 'q' in model:
        q = model['q']
    else:
        objective = lambda q: -float(np.log(1 - q * F + q * sums).sum())
        opt = minimize_scalar(objective, bounds=(0.01, 0.999), method='bounded')
        q = float(opt.x)
    if not 0 < q < 1:
        raise ValueError('Q must lie strictly between 0 and 1')
    null = 1 - q * F
    denom = null + q * sums
    logodds = np.log(q * lr / null)
    frame = pd.DataFrame(dict(a=ei, b=ej, raw_distance=raw, residual_distance=r, log_lr=np.log(lr), log_photometry=logphot, local_probability=q * lr / denom[ei], cost=-logodds))
    return (frame, dict(q=q, search_mass=F, local_null_probability=null / denom, density=rho, systematic_shift=shifts))

def assign(n_a, n_b, edges, margin_limit=50):
    useful = edges[edges.cost < 0]
    labels = find_components(n_a, n_b, useful.a.to_numpy(), useful.b.to_numpy())
    groups = defaultdict(list)
    for e in useful.itertuples():
        groups[int(labels[e.a])].append(e)
    chosen = {}
    margins = {}
    components = []
    for group in groups.values():
        aa = sorted({e.a for e in group})
        bb = sorted({e.b for e in group})
        ai = {a: i for i, a in enumerate(aa)}
        bi = {b: i for i, b in enumerate(bb)}
        mat = np.full((len(aa), len(bb) + len(aa)), np.inf)
        for r in range(len(aa)):
            mat[r, len(bb) + r] = 0
        for e in group:
            mat[ai[e.a], bi[e.b]] = e.cost
        rr, cc = linear_sum_assignment(mat)
        optimum = float(mat[rr, cc].sum())
        components.append((len(aa), len(bb), len(group), optimum))
        for r, c in zip(rr, cc):
            if c >= len(bb):
                continue
            a, b = (aa[r], bb[c])
            chosen[a] = b
            if len(aa) <= margin_limit:
                old = mat[r, c]
                mat[r, c] = np.inf
                ar, ac = linear_sum_assignment(mat)
                margins[a] = max(0.0, float(mat[ar, ac].sum() - optimum))
                mat[r, c] = old
    return (chosen, margins, components)

def run_single(av, am, bv, bm, settings=None, model=None, variable_a=None, variable_b=None, variability_a=None, variability_b=None):
    settings = settings or Settings()
    if settings.radius <= 0 or not 0 <= settings.broad_fraction < 1 or settings.broad_scale_floor <= 0:
        raise ValueError('Invalid settings')
    av, bv = (normalize(av), normalize(bv))
    if len(bv) < 2 or len(av) == 0:
        raise ValueError('Nonempty catalogs and at least two modern objects required')
    for variable, count in [(variable_a, len(av)), (variable_b, len(bv))]:
        if variable is not None and np.asarray(variable).shape != (count,):
            raise ValueError('Variable mask shape mismatch')
    if variability_a is None and variable_a is not None:
        variability_a = np.where(variable_a, np.nan if settings.variable_sigma is None else settings.variable_sigma, 0.0)
    if variability_b is None and variable_b is not None:
        variability_b = np.where(variable_b, np.nan if settings.variable_sigma is None else settings.variable_sigma, 0.0)
    for scatter, count in [(variability_a, len(av)), (variability_b, len(bv))]:
        if scatter is not None and (np.asarray(scatter).shape != (count,) or np.any(np.asarray(scatter) < 0)):
            raise ValueError('Invalid variability sigma array')
    if model is None:
        model = calibrate_single(av, am, bv, bm, settings, variable_a, variable_b)
    if not settings.systematic:
        model = {**model, 'systematic_coefficients': np.zeros((6, 2)).tolist()}
    edges, info = evidence(av, am, bv, bm, model, settings, variability_a, variability_b)
    info['variable_a'] = np.zeros(len(av), bool) if variable_a is None else variable_a
    info['variable_b'] = np.zeros(len(bv), bool) if variable_b is None else variable_b
    chosen, margins, components = assign(len(av), len(bv), edges, settings.global_margin_limit)
    model = {**model, 'q': info['q'], 'settings': asdict(settings), 'search_mass': info['search_mass']}
    return (chosen, edges, info, margins, components, model)

def save_single(output, an, av, am, bn, bv, result):
    chosen, edges, info, margins, components, model = result
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    lookup = {(e.a, e.b): e for e in edges.itertuples()}
    rows = []
    matches = []
    alternatives = defaultdict(list)
    for e in edges.sort_values(['a', 'cost', 'b']).itertuples():
        if len(alternatives[e.a]) < 3:
            alternatives[e.a].append(e)
    for a, name in enumerate(an):
        b = chosen.get(a)
        row = dict(index1=name, index2=bn[b] if b is not None else '', mag=am[a], dist=np.nan, residual_dist=np.nan, local_probability=np.nan, local_null_probability=info['local_null_probability'][a], global_cost_gap=margins.get(a, np.nan), systematic_shift=info['systematic_shift'][a], density_arcsec2=info['density'][a], variable_a=bool(info['variable_a'][a]), variable_b=bool(info['variable_b'][b]) if b is not None else False, status='unmatched')
        if b is not None:
            e = lookup[a, b]
            row.update(dist=e.raw_distance, residual_dist=e.residual_distance, local_probability=e.local_probability, status='review' if margins.get(a, 0) < math.log(10) or e.residual_distance > 5 * model['sigma'] else 'matched')
            matches.append((name, bn[b], am[a], e.raw_distance))
        for rank, e in enumerate(alternatives[a], 1):
            row[f'candidate{rank}'] = bn[e.b]
            row[f'candidate{rank}_cost'] = e.cost
            row[f'candidate{rank}_raw_dist'] = e.raw_distance
        rows.append(row)
    pd.DataFrame(matches, columns=['index1', 'index2', 'mag', 'dist']).to_csv(output, index=False)
    diagnostics = pd.DataFrame(rows)
    diagnostics.to_csv(output.with_suffix('.diagnostics.csv'), index=False)
    diagnostics[diagnostics.status == 'matched'][['index1', 'index2', 'mag', 'dist']].to_csv(output.with_suffix('.secure.csv'), index=False)
    diagnostics[diagnostics.status != 'matched'].to_csv(output.with_suffix('.review.csv'), index=False)
    output.with_suffix('.model.json').write_text(json.dumps(model, indent=2, allow_nan=False) + '\n')
    edges.to_csv(output.with_suffix('.candidates.csv'), index=False)
    return pd.DataFrame(rows)


# CATALOG IO
B1875 = Time('B1875')
J2000 = Time('J2000')
BANDS = np.array([0, 1, 6, 8, 10, 11, 12, 13, 14, 18, 4, 3, 16])

def read_cd(root):
    records = []
    for i, line in enumerate((root / 'cat/cd_vol1_curated.txt').read_text().splitlines()):
        if len(line) != 30 or line[:2] != 'CD':
            raise ValueError(f'Invalid CD record {i + 1}')
        zone, num, suppl = (int(line[2:5]), int(line[5:10]), line[10])
        if not -31 <= zone <= -22:
            raise ValueError(f'CD outside volume I: {line}')
        mag = float(line[11:15])
        ra = 15 * (int(line[15:17]) + int(line[17:19]) / 60 + float(line[19:23]) / 3600)
        dec = -(int(line[24:26]) + float(line[26:30]) / 60)
        records.append(dict(name=f'CD{zone}{num:5d}{suppl}'.rstrip(), zone=zone, num=num, suppl=suppl, mag=mag if mag < 20 else np.nan, mag_code=mag, ra=ra, dec=dec, raw=line, active=suppl != 'D' and (mag < 20 or mag == 30), double=False, color=False))
    df = pd.DataFrame(records)
    if df.name.duplicated().any():
        raise ValueError('Nonunique CD identifiers')
    flag_audit = {}
    for kind, field in [('dpl', 'double'), ('color', 'color')]:
        for zone in range(22, 32):
            p = root / f'cd/{kind}_{zone}.txt'
            entries = [x.strip() for x in p.read_text().splitlines() if x.strip()]
            nums = [int(x.split()[0]) for x in entries]
            uncertain = [int(x.split()[0]) for x in entries if '?' in x]
            if uncertain:
                field_u = field + '_uncertain'
                if field_u not in df:
                    df[field_u] = False
                df.loc[(df.zone == -zone) & df.num.isin(uncertain) & (df.suppl == ' '), field_u] = True
            if len(set(nums)) != len(nums):
                raise ValueError(f'Duplicate flag in {p}')
            mask = (df.zone == -zone) & df.num.isin(nums) & (df.suppl == ' ')
            absent = sorted(set(nums) - set(df.loc[mask, 'num']))
            if absent:
                raise ValueError(f'Orphan flags {p}: {absent}')
            df.loc[mask, field] = True
            flag_audit[p.name] = int(mask.sum())
    return (df, flag_audit)

def xyz(ra, dec):
    a, d = (np.deg2rad(ra), np.deg2rad(dec))
    return np.column_stack((np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)))

def read_ppm(root):
    rows = []
    for s in (root / 'cat/ppm.txt').read_text().splitlines():
        n = int(s[1:7])
        mag = float(s[19:23]) if s[19:23].strip() else np.nan
        ra = 15 * (int(s[27:29]) + int(s[30:32]) / 60 + float(s[33:39]) / 3600)
        dec = (int(s[42:44]) + int(s[45:47]) / 60 + float(s[48:53]) / 3600) * (-1 if s[41] == '-' else 1)
        is_v = s[130:131] == 'V' or 400001 <= n <= 400321
        rows.append((f'PPM {n}', n, ra, dec, float(s[55:62]) * 15 * np.cos(np.deg2rad(dec)), float(s[63:69]), mag, 'PPM_V' if is_v else f"PPM_nonV_{s[130:131].strip() or 'blank'}", s[9:18], s[126:128]))
    return pd.DataFrame(rows, columns=['name', 'ppm', 'ra', 'dec', 'pmra_cosdec', 'pmdec', 'mag', 'band', 'dm', 'flags'])

def ppm_at(ppm, epoch, b1875=False):
    """pmRA source is seconds of time/year, converted to mu_alpha*cos(delta)."""
    c = SkyCoord(ra=ppm.ra.to_numpy() * u.deg, dec=ppm.dec.to_numpy() * u.deg, pm_ra_cosdec=ppm.pmra_cosdec.to_numpy() * u.arcsec / u.yr, pm_dec=ppm.pmdec.to_numpy() * u.arcsec / u.yr, frame=FK5(equinox=J2000), obstime=J2000)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ErfaWarning)
        p = c.apply_space_motion(new_obstime=epoch)
    p = SkyCoord(ra=p.ra, dec=p.dec, frame=FK5(equinox=J2000))
    if b1875:
        p = p.transform_to(FK4(equinox=B1875, obstime=B1875))
    return p.cartesian.xyz.value.T

def gsc_at_1875(gsc):
    c = SkyCoord(ra=gsc.ra.to_numpy() * u.deg, dec=gsc.dec.to_numpy() * u.deg, frame=FK5(equinox=J2000))
    return c.transform_to(FK4(equinox=B1875, obstime=B1875)).cartesian.xyz.value.T

def decode_region(path):
    """Vectorized exact port of gsc/src/decode_c.c encoding version 2."""
    with path.open('rb') as f:
        length = int(f.read(4))
        header = f.read(length - 4).decode('ascii')
        tokens = header.split()
        ver, region, nobj = map(int, tokens[:3])
        if ver != 2:
            raise ValueError(f'Unsupported encoding {ver}: {path}')
        amin, amax, dmin, dmax, off, sra, sdec, spos, smag = map(float, tokens[3:12])
        npl = int(tokens[12])
        plates = tokens[13:13 + npl]
        epochs = tokens[13 + npl:13 + 2 * npl]
        data = f.read()
    if len(data) % 12:
        raise ValueError(f'Truncated GSC file {path}')
    c = np.frombuffer(data, dtype=np.uint8).reshape(-1, 12).astype(np.int64)
    if len(c) != nobj:
        raise ValueError(f'GSC count mismatch {path}: {len(c)} != {nobj}')
    ids = (c[:, 0] & 127) << 7 | c[:, 1] >> 1
    da = (c[:, 1] & 1) << 21 | c[:, 2] << 13 | c[:, 3] << 5 | c[:, 4] >> 3
    dd = (c[:, 4] & 7) << 16 | c[:, 5] << 8 | c[:, 6]
    pl = c[:, 11] & 15
    ba = c[:, 10] >> 1 & 15
    if np.any(pl >= npl) or np.any(ba >= len(BANDS)):
        raise ValueError(f'GSC code {path}')
    epochs = np.array([float(e) if re.fullmatch('\\d{4}\\.\\d+', e) else np.nan for e in epochs])
    if len(epochs) != npl:
        raise ValueError(f'GSC missing epoch list {path}')
    return pd.DataFrame(dict(name=[f'GSC {region:05d}-{i:05d}' for i in ids], ra=(da / sra + amin) % 360, dec=dd / sdec + dmin, poserr=(c[:, 7] << 1 | c[:, 8] >> 7) / spos, mag=(c[:, 9] << 3 | c[:, 10] >> 5) / smag + off, magerr=(c[:, 8] & 127) / smag, band=[f'GSC_{x}' for x in BANDS[ba]], classification=c[:, 11] >> 4 & 7, plate=np.array(plates)[pl], epoch=epochs[pl]))

def read_gsc(root, lower=-34, upper=-19):
    """Use region headers to select a generous J2000 declination halo."""
    parts = []
    files = 0
    for p in sorted((root / 'gsc').glob('*/*.GSC')):
        with p.open('rb') as f:
            n = int(f.read(4))
            h = f.read(n - 4).decode('ascii').split()
        if float(h[6]) < lower or float(h[5]) > upper:
            continue
        files += 1
        parts.append(decode_region(p))
    obs = pd.concat(parts, ignore_index=True)
    obs = obs[(obs.dec >= lower) & (obs.dec <= upper)].copy()
    total = len(obs)
    obs = obs[obs.classification.isin([0, 2])].copy()
    fallback = float(obs.epoch.mean())
    mean_epoch = obs.groupby('name').epoch.mean()
    mean_ra_vec = pd.DataFrame(xyz(obs.ra, obs.dec), index=obs.index).groupby(obs.name).mean()
    v = mean_ra_vec.to_numpy()
    v /= np.linalg.norm(v, axis=1)[:, None]
    coords = pd.DataFrame(dict(name=mean_ra_vec.index, ra=np.rad2deg(np.arctan2(v[:, 1], v[:, 0])) % 360, dec=np.rad2deg(np.arcsin(v[:, 2])))).set_index('name')
    best = obs.sort_values(['name', 'magerr', 'poserr', 'plate']).drop_duplicates('name').set_index('name')
    best[['ra', 'dec']] = coords[['ra', 'dec']]
    best['epoch'] = mean_epoch
    best['epoch_missing'] = best.epoch.isna()
    best['epoch'] = best.epoch.fillna(fallback)
    best['nobs'] = obs.groupby('name').size()
    best = best.reset_index()
    audit = dict(region_files=files, observations_in_halo=total, stellar_observations=len(obs), objects=len(best), epoch_fallback=fallback, epoch_missing=int(best.epoch_missing.sum()), bands={str(k): int(v) for k, v in best.band.value_counts().items()}, epoch_min=float(best.epoch.min()), epoch_max=float(best.epoch.max()))
    return (best, audit)


# PREPARE CATALOGS

def report(message):
    print(message, flush=True)

def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')

def assign_pairs(na, nb, ei, ej, cost, empty_cost):
    labels = find_components(na, nb, ei, ej)
    groups = {}
    for a, b, c in zip(ei, ej, cost):
        groups.setdefault(int(labels[a]), []).append((int(a), int(b), float(c)))
    matches = []
    for edges in groups.values():
        aa = sorted({e[0] for e in edges})
        bb = sorted({e[1] for e in edges})
        ai = {a: k for k, a in enumerate(aa)}
        bi = {b: k for k, b in enumerate(bb)}
        matrix = np.full((len(aa), len(bb) + len(aa)), np.inf)
        for k in range(len(aa)):
            matrix[k, len(bb) + k] = empty_cost
        for a, b, c in edges:
            matrix[ai[a], bi[b]] = c
        rows, cols = linear_sum_assignment(matrix)
        matches.extend(((aa[r], bb[c]) for r, c in zip(rows, cols) if c < len(bb)))
    return matches

def preidentify(ppm, pvec, gsc, gvec, out):
    epoch = float(gsc.epoch.mean())
    pmid = ppm_at(ppm, Time(epoch, format='jyear'))
    motion = np.hypot(ppm.pmra_cosdec, ppm.pmdec).to_numpy()
    radius = 20 + motion * np.max(np.abs(gsc.epoch.to_numpy() - epoch))
    cand = cKDTree(gvec).query_ball_point(pmid, 2 * np.sin(radius / (2 * ARCSEC_PER_RAD)))
    ei = np.repeat(np.arange(len(ppm)), [len(x) for x in cand])
    ej = np.array([j for x in cand for j in x], dtype=int)
    report(f'PPM/GSC: {len(ei)} candidates for epoch-aware distances')
    pred = ppm_at(ppm.iloc[ei], Time(gsc.epoch.to_numpy()[ej], format='jyear'))
    dist = separation(pred, gvec[ej])
    order = np.argsort(dist)
    besta = {}
    bestb = {}
    for k in order:
        besta.setdefault(int(ei[k]), k)
        bestb.setdefault(int(ej[k]), k)
    core = np.array([k for k in besta.values() if bestb[int(ej[k])] == k and dist[k] < 10], dtype=int)
    cut = float(np.clip(np.quantile(dist[core], 0.995), 3, 10))
    keep = dist <= cut
    pairs = assign_pairs(len(ppm), len(gsc), ei[keep], ej[keep], (dist[keep] / cut) ** 2, 1.01)
    lookup = {(int(a), int(b)): float(d) for a, b, d in zip(ei[keep], ej[keep], dist[keep])}
    acount = np.bincount(ei[keep], minlength=len(ppm))
    bcount = np.bincount(ej[keep], minlength=len(gsc))
    rows = [dict(ppm=ppm.name.iloc[a], gsc=gsc.name.iloc[b], dist_arcsec=lookup[a, b], epoch=gsc.epoch.iloc[b], epoch_imputed=bool(gsc.epoch_missing.iloc[b]), ambiguous=bool(acount[a] > 1 or bcount[b] > 1)) for a, b in pairs]
    pd.DataFrame(rows).to_csv(out / 'ppm_gsc.csv', index=False)
    shifted = xyz((gsc.ra + 0.5) % 360, gsc.dec)
    sd, _ = cKDTree(shifted).query(pmid)
    false = int(np.sum(sd < 2 * np.sin(cut / (2 * ARCSEC_PER_RAD))))
    audit = dict(epoch_mean=epoch, epoch_policy='PPM propagated to each GSC effective plate epoch; FK5 J2000 equinox', radius_arcsec=cut, mutual_core_count=len(core), core_quantiles={str(q): float(np.quantile(dist[core], q)) for q in [0.5, 0.9, 0.99, 0.995]}, matches=len(pairs), ambiguous_matches=sum((r['ambiguous'] for r in rows)), shifted_control_matches=false, shift_degrees_ra=0.5)
    return (pairs, audit)

def fit_band(x, y, holdout, existing=False, min_sigma=0.3, degree=2):
    train = ~holdout
    valid = holdout
    if train.sum() < 40 or valid.sum() < 10:
        return None
    keep = train.copy()
    for _ in range(4):
        coef = np.polynomial.polynomial.polyfit(x[keep], y[keep], degree)
        residual = y - np.polynomial.polynomial.polyval(x, coef)
        scale = max(0.1, 1.4826 * np.median(np.abs(residual[keep] - np.median(residual[keep]))))
        keep = train & (np.abs(residual) < 3 * scale)
    fitted = coef.copy()
    if existing:
        coef = np.array([-0.157169, 1.188316, -0.02213])
    r = y - np.polynomial.polynomial.polyval(x, coef)
    rv = r[valid]
    core = rv[np.abs(rv - np.median(rv)) < 3 * max(0.1, 1.4826 * np.median(np.abs(rv - np.median(rv))))]
    sigma = max(min_sigma, float(np.sqrt(np.mean(core ** 2))))
    return dict(coefficients=coef.tolist(), refit_coefficients=fitted.tolist(), n_train=int(keep.sum()), n_validation=int(valid.sum()), validation_rmse=float(np.sqrt(np.mean(rv ** 2))), validation_rmse_refit=float(np.sqrt(np.mean((y[valid] - np.polynomial.polynomial.polyval(x[valid], fitted)) ** 2))), sigma_mag=sigma, mag_training_range=[float(x[keep].min()), float(x[keep].max())], method='existing trig.cpp V polynomial' if existing else f'iterated 3-MAD-clipped degree-{degree} least squares')

def calibrate_photometry(cd, cvec, modern, mvec, out, gsc, pairs):
    report('Photometric and positional calibration')
    dist, index = cKDTree(mvec).query(cvec, k=2)
    theta = 2 * np.arcsin(np.clip(dist / 2, 0, 1)) * ARCSEC_PER_RAD
    reverse = cKDTree(cvec).query(mvec)[1]
    primary = index[:, 0]
    clean = cd.active.to_numpy() & ~cd.double.to_numpy() & ~cd.color.to_numpy() & np.isfinite(cd.mag)
    clean = clean & (theta[:, 0] < 60) & (theta[:, 1] > np.maximum(60, 2 * theta[:, 0])) & (reverse[primary] == np.arange(len(cd)))
    holdout = np.floor(cd.ra.to_numpy() / 5).astype(int) % 5 == 0
    bands = {}
    gsc_by_ppm = {a: b for a, b in pairs}
    phot_rows = []
    for band in sorted(modern.band.unique()):
        if band.startswith('GSC'):
            aa = np.array([a for a in np.flatnonzero(clean) if primary[a] in gsc_by_ppm and gsc.band.iloc[gsc_by_ppm[primary[a]]] == band], dtype=int)
            bb = np.array([gsc_by_ppm[primary[a]] for a in aa], dtype=int)
            x = gsc.mag.to_numpy()[bb]
            y = cd.mag.to_numpy()[aa]
            selected = np.isfinite(x) & (y < 9.5) & (x < 12.5)
            aa = aa[selected]
            x = x[selected]
            y = y[selected]
            validation = holdout[aa]
        else:
            use = clean & (modern.band.to_numpy()[primary] == band) & np.isfinite(modern.mag.to_numpy()[primary])
            aa = np.flatnonzero(use)
            x = modern.mag.to_numpy()[primary[use]]
            y = cd.mag.to_numpy()[use]
            validation = holdout[use]
        result = fit_band(x, y, validation, existing=band == 'PPM_V', min_sigma=0.3 if band == 'PPM_V' else 0.5, degree=1 if band.startswith('GSC') else 2)
        if result:
            bands[band] = result
        phot_rows.extend((dict(cd=cd.name.iloc[a], band=band, raw_mag=float(xx), cd_mag=float(yy), validation=bool(v)) for a, xx, yy, v in zip(aa, x, y, validation)))
        report(f'  {band}: {len(x)} trusted anchors, fitted={result is not None}')
    pd.DataFrame(phot_rows).to_csv(out / 'photometric_anchors.csv', index=False)
    modern['mag_cd'] = np.nan
    modern['sigma_mag'] = 1.5
    modern['photometry_extrapolated'] = False
    for band, cal in bands.items():
        use = modern.band == band
        x = modern.loc[use, 'mag'].to_numpy()
        coef = cal['coefficients']
        lo, hi = cal['mag_training_range']
        xx = np.clip(x, lo, hi)
        converted = np.polynomial.polynomial.polyval(xx, coef)
        deriv = np.polynomial.polynomial.polyval(xx, np.polynomial.polynomial.polyder(coef))
        converted += np.maximum(deriv, 0.5) * (x - xx)
        modern.loc[use, 'mag_cd'] = converted
        extra = (x < lo) | (x > hi)
        modern.loc[use, 'photometry_extrapolated'] = extra
        modern.loc[use, 'sigma_mag'] = cal['sigma_mag'] * np.where(extra, 1.5, 1.0)
    direct = []
    cdmap = {(r.zone, r.num): i for i, r in cd.iterrows() if r.suppl == ' ' and r.active and (not r.double)}
    for b, r in modern[modern.source == 'PPM'].iterrows():
        dm = r.get('dm', '')
        if isinstance(dm, str) and dm[:1] == '-' and dm[1:3].strip().isdigit():
            zone = -int(dm[1:3])
            num = dm[3:8].strip()
            if zone <= -23 and num.isdigit() and ((zone, int(num)) in cdmap):
                a = cdmap[zone, int(num)]
                direct.append(float(separation(cvec[a], mvec[b])))
    direct = np.array(direct)
    core = direct[direct < 120]
    sigma_ppm = float(np.sqrt(np.mean(core ** 2)))
    gt = theta[clean & (modern.source.to_numpy()[primary] == 'GSC'), 0]
    sigma_gsc = max(sigma_ppm, float(np.sqrt(np.mean(gt ** 2))))
    radius = float(np.clip(np.quantile(core, 0.997) * 1.2, 90, 180))
    calibration = dict(bands=bands, sigma_pos_ppm=sigma_ppm, sigma_pos_gsc=sigma_gsc, radius_arcsec=radius, direct_ppm_anchors=len(direct), direct_ppm_over_120=int(np.sum(direct >= 120)), direct_ppm_quantiles={str(q): float(np.quantile(direct, q)) for q in [0.5, 0.9, 0.95, 0.99, 0.995, 0.999]}, position_caveat='PPM RMS uses direct CD anchors truncated at 120 arcsec; GSC anchors limited to 60 arcsec and no PM.', photometry_split='5-degree RA sectors, sector mod 5 == 0 held out', clean_photometry_anchors=int(clean.sum()), missing_photometric_density=0.05, max_dmag=3.0, modern_mag_limit=13.5, p_empty=0.01, separation_mean=34.9, separation_sigma=13.65, contrast_sigma=0.915, pair_radius=80.0, p_single=0.15)
    pd.DataFrame(dict(cd=cd.name[clean], modern=modern.name.to_numpy()[primary[clean]], theta=theta[clean, 0], cd_mag=cd.mag[clean], raw_mag=modern.mag.to_numpy()[primary[clean]], band=modern.band.to_numpy()[primary[clean]], validation=holdout[clean])).to_csv(out / 'calibration_anchors.csv', index=False)
    return (modern, calibration)

def calibrate_doubles(cd, cvec, modern, mvec, cal):
    dbl = np.flatnonzero(cd.double.to_numpy() & cd.active.to_numpy())
    controls = []
    for a in dbl:
        mask = (cd.zone == cd.zone.iloc[a]) & ~cd.double & cd.active & (np.abs(cd.mag - cd.mag.iloc[a]) < 0.5)
        ids = np.flatnonzero(mask.to_numpy())
        da = np.abs((cd.ra.to_numpy()[ids] - cd.ra.iloc[a] + 180) % 360 - 180)
        controls.extend(ids[np.argsort(da)[:3]].tolist())
    sampled = np.concatenate([dbl, np.asarray(controls, dtype=int)])
    cand = neighbors(cvec[sampled], mvec, cal['radius_arcsec'])
    hits = []
    seps = []
    contrasts = []
    mag = modern.mag_cd.to_numpy()
    for ai, js in zip(sampled, cand):
        found = []
        if js:
            ds = separation(cvec[ai], mvec[js])
            primary = int(js[int(np.argmin(ds))])
            if ds.min() < 60 and np.isfinite(mag[primary]):
                for b in js:
                    if b == primary or not np.isfinite(mag[b]):
                        continue
                    s = float(separation(mvec[primary], mvec[b]))
                    dm = abs(mag[primary] - mag[b])
                    if 10 < s < 80 and dm < 2 and (not np.isfinite(cd.mag.iloc[ai]) or abs(cd.mag.iloc[ai] - min(mag[b], mag[primary])) < 1.5):
                        found.append((s, dm))
        hits.append(bool(found))
        if len(found) == 1 and len(hits) <= len(dbl):
            seps.append(found[0][0])
            contrasts.append(found[0][1])
    pdbl = float(np.mean(hits[:len(dbl)]))
    pctrl = float(np.mean(hits[len(dbl):]))
    q = float(np.clip((pdbl - pctrl) / (1 - pctrl), 0.01, 0.99))
    cal['p_single'] = 1 - q
    if len(seps) >= 50:
        cal['separation_mean'] = float(np.mean(seps))
        cal['separation_sigma'] = float(np.std(seps, ddof=1))
        cal['contrast_sigma'] = float(np.sqrt(np.mean(np.asarray(contrasts) ** 2)))
    cal['double_calibration'] = dict(n_double=len(dbl), n_control=len(controls), double_pair_rate=pdbl, control_pair_rate=pctrl, n_unique_pairs=len(seps), p_single_estimate=1 - q, caveat='Empirical excess estimator assumes equal chance-pair rates and perfect detectability inside 10-80 arcsec / contrast<2. Not independently labelled; verify by sensitivity analysis.')
    return cal

def prepare_catalogs(root, out):
    root, out = (Path(root), Path(out))
    out.mkdir(parents=True, exist_ok=True)
    cd, flags = read_cd(root)
    cvec = xyz(cd.ra, cd.dec)
    report(f'CD: {len(cd)} records; {cd.double.sum()} doubles')
    ppm = read_ppm(root)
    all_ppm = len(ppm)
    pvec = ppm_at(ppm, B1875, b1875=True)
    dec = np.rad2deg(np.arcsin(pvec[:, 2]))
    keep = (dec >= -34) & (dec <= -19)
    ppm = ppm.loc[keep].reset_index(drop=True)
    pvec = pvec[keep]
    report(f'PPM: {len(ppm)} in halo of {all_ppm} total')
    gsc, ga = read_gsc(root)
    report(f'GSC: {len(gsc)} unique stellar objects')
    gvec = xyz(gsc.ra, gsc.dec)
    pairs, pa = preidentify(ppm, pvec, gsc, gvec, out)
    used = {b for a, b in pairs}
    aliases = {a: gsc.name.iloc[b] for a, b in pairs}
    ppm['source'] = 'PPM'
    ppm['gsc_alias'] = [aliases.get(i, '') for i in range(len(ppm))]
    rem = gsc.loc[~gsc.index.isin(used)].copy()
    rem['source'] = 'GSC'
    rem['gsc_alias'] = ''
    modern = pd.concat([ppm, rem], ignore_index=True)
    mvec = np.vstack([pvec, gsc_at_1875(rem)])
    modern, cal = calibrate_photometry(cd, cvec, modern, mvec, out, gsc, pairs)
    cal = calibrate_doubles(cd, cvec, modern, mvec, cal)
    for k, col in enumerate(['x', 'y', 'z']):
        cd[col] = cvec[:, k]
        modern[col] = mvec[:, k]
    cd.to_csv(out / 'cd.csv', index=False)
    modern.to_csv(out / 'modern.csv', index=False)
    save_json(out / 'calibration.json', cal)
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [root / 'cat/cd_vol1_curated.txt', root / 'cat/ppm.txt', *sorted((root / 'cd').glob('*.txt'))]}
    audit = dict(cd_records=len(cd), active=int(cd.active.sum()), doubles=int(cd.double.sum()), colors=int(cd.color.sum()), flags=flags, ppm_total=all_ppm, ppm_halo=len(ppm), ppm_bands={str(k): int(v) for k, v in ppm.band.value_counts().items()}, gsc=ga, preidentification=pa, modern_union=len(modern), input_sha256=hashes)
    save_json(out / 'audit.json', audit)
    report(json.dumps(dict(sigma_ppm=cal['sigma_pos_ppm'], sigma_gsc=cal['sigma_pos_gsc'], radius=cal['radius_arcsec'], p_single=cal['p_single'])))


# CROSS DOUBLE LIKELIHOOD

@dataclass(frozen=True)
class Hypothesis:
    members: tuple
    cost: float
    primary: int = -1
    distance: float = math.nan
    pair_separation: float = math.nan

def solve_scip(group, hypotheses, time_limit=120):
    from pyscipopt import Model, quicksum
    model = Model('cd_double_component')
    try:
        model.hideOutput()
        model.setParam('parallel/maxnthreads', 1)
        model.setParam('randomization/randomseedshift', 0)
        model.setParam('limits/time', float(time_limit))
        model.setParam('limits/gap', 0.0)
        model.setParam('limits/absgap', 0.0)
        xs = {(a, k): model.addVar(vtype='B', name=f'x_{a}_{k}') for a in group for k in range(len(hypotheses[a]))}
        model.setObjective(quicksum((hypotheses[a][k].cost * x for (a, k), x in xs.items())), 'minimize')
        for a in group:
            model.addCons(quicksum((xs[a, k] for k in range(len(hypotheses[a])))) == 1)
        use = defaultdict(list)
        for (a, k), x in xs.items():
            for b in hypotheses[a][k].members:
                use[b].append(x)
        for variables in use.values():
            if len(variables) > 1:
                model.addCons(quicksum(variables) <= 1)
        nvars, ncons = (model.getNVars(), model.getNConss())
        model.optimize()
        status = str(model.getStatus())
        if status != 'optimal' or model.getNSols() == 0:
            raise RuntimeError(f'SCIP did not prove optimality: {status}; group size={len(group)}')
        solution = model.getBestSol()
        selected = {}
        for a in group:
            ks = [k for k in range(len(hypotheses[a])) if model.getSolVal(solution, xs[a, k]) > 0.5]
            if len(ks) != 1:
                raise AssertionError('SCIP assignment is not integral')
            selected[a] = ks[0]
        used = [b for a, k in selected.items() for b in hypotheses[a][k].members]
        if len(used) != len(set(used)):
            raise AssertionError('SCIP reused a modern object')
        objective = sum((hypotheses[a][k].cost for a, k in selected.items()))
        if not math.isclose(objective, model.getObjVal(), rel_tol=1e-08, abs_tol=1e-06):
            raise AssertionError('SCIP objective and extracted assignment disagree')
        return (selected, dict(status=status, variables=nvars, constraints=ncons, seconds=float(model.getSolvingTime()), objective=float(objective), gap=float(model.getGap())))
    finally:
        model.freeProb()

def solve_hypotheses(hypotheses, nmodern, time_limit=120, component_labels=None, double_flags=None):
    ei = []
    ej = []
    for a, hs in enumerate(hypotheses):
        for b in sorted({b for h in hs for b in h.members}):
            ei.append(a)
            ej.append(b)
    labels = find_components(len(hypotheses), nmodern, ei, ej) if component_labels is None else component_labels
    groups = defaultdict(list)
    for a in range(len(hypotheses)):
        groups[int(labels[a])].append(a)
    selected = {}
    stats = Counter()
    mip = []
    for group_number, group in enumerate(groups.values(), 1):
        use_scip = bool(np.any(np.asarray(double_flags)[group])) if double_flags is not None else any((len(h.members) == 2 for a in group for h in hypotheses[a]))
        if len(group) == 1 and (not use_scip) and (double_flags is None):
            a = group[0]
            selected[a] = min(range(len(hypotheses[a])), key=lambda k: hypotheses[a][k].cost)
            stats['trivial'] += 1
        elif use_scip:
            if len(group) > 1000:
                print(f'SCIP: component with {len(group)} CD', flush=True)
            chosen, info = solve_scip(group, hypotheses, time_limit)
            selected.update(chosen)
            mip.append(info)
            stats['scip'] += 1
            if len(group) > 1000:
                print(f'SCIP: {info}', flush=True)
        else:
            bs = sorted({h.members[0] for a in group for h in hypotheses[a] if h.members})
            bp = {b: i for i, b in enumerate(bs)}
            matrix = np.full((len(group), len(bs) + len(group)), np.inf)
            lookup = {}
            for row, a in enumerate(group):
                for k, h in enumerate(hypotheses[a]):
                    col = bp[h.members[0]] if h.members else len(bs) + row
                    matrix[row, col] = h.cost
                    lookup[row, col] = k
            rr, cc = linear_sum_assignment(matrix)
            selected.update({group[r]: lookup[r, c] for r, c in zip(rr, cc)})
            stats['hungarian'] += 1
    used = [b for a, k in selected.items() for b in hypotheses[a][k].members]
    if len(used) != len(set(used)) or len(selected) != len(hypotheses):
        raise AssertionError('Invalid global assignment')
    return (selected, dict(components=dict(stats), scip=mip, objective=sum((hypotheses[a][k].cost for a, k in selected.items()))))

def write_result(cd, modern, hypotheses, selected, diagnostics, out, config, stats, catalog_output=None):
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    fixed = []
    for a, r in cd.iterrows():
        h = hypotheses[a][selected[a]]
        members = sorted(h.members, key=lambda b: (modern.source.iloc[b] != 'PPM', b != h.primary, modern.name.iloc[b]))
        names = [modern.name.iloc[b] for b in members] + ['', '']
        row = dict(cd=r['name'], original=r.raw, zone=int(r.zone), number=int(r.num), supplement=r.suppl, double=bool(r.double), color=bool(r.color), double_uncertain=bool(r.get('double_uncertain', False)), color_uncertain=bool(r.get('color_uncertain', False)), type=len(members), id1=names[0], id2=names[1], distance_arcsec=float(separation(r[['x', 'y', 'z']].to_numpy(dtype=float), modern.iloc[members[0]][['x', 'y', 'z']].to_numpy(dtype=float))) if members else math.nan, pair_separation_arcsec=h.pair_separation, local_probability=math.nan if config.get('model') == 'robust_cd_v1' else math.exp(-h.cost), log_odds_vs_empty=-h.cost, astrometric_primary=modern.name.iloc[h.primary] if members else '', candidates_single=diagnostics[a][0], candidates_double=diagnostics[a][1])
        reasons = []
        if not r.active:
            reasons.append('excluded_deleted_or_nonstellar')
        elif not members:
            reasons.append('unmatched')
        if row['double_uncertain'] or row['color_uncertain']:
            reasons.append('uncertain_historical_flag')
        if r.mag_code == 30:
            reasons.append('variable_position_only')
        costs = sorted((x.cost for x in hypotheses[a]))
        row['local_cost_margin'] = costs[1] - costs[0] if len(costs) > 1 else math.nan
        if row['local_cost_margin'] < 1:
            reasons.append('locally_ambiguous')
        if selected[a] != min(range(len(hypotheses[a])), key=lambda k: hypotheses[a][k].cost):
            reasons.append('global_conflict_changes_local_choice')
        for i, b in enumerate(members, 1):
            mr = modern.iloc[b]
            row[f'gsc_alias{i}'] = mr.gsc_alias if isinstance(mr.gsc_alias, str) else ''
            row[f'band{i}'] = mr.band
            row[f'mag_cd{i}'] = mr.mag_cd
            row[f'variable{i}'] = bool(mr.get('variable', False))
            row[f'variable_name{i}'] = mr.get('variable_name', '')
            row[f'variability_sigma_cd{i}'] = mr.get('variability_sigma_cd', 0.0)
            if mr.get('variable', False):
                reasons.append('variable_photometry_broadened')
            if not np.isfinite(mr.mag_cd):
                reasons.append('missing_photometry')
            if mr.photometry_extrapolated:
                reasons.append('photometric_extrapolation')
            if mr.source == 'GSC':
                reasons.append('gsc_no_proper_motion')
            if mr.source == 'GSC' and mr.classification == 2:
                reasons.append('gsc_blend')
            if mr.source == 'GSC' and mr.epoch_missing:
                reasons.append('gsc_epoch_imputed')
            if mr.source == 'PPM' and isinstance(mr.get('flags'), str) and mr['flags'].strip():
                reasons.append('ppm_flag')
        row['review'] = ';'.join(sorted(set(reasons)))
        rows.append(row)
        line = f"{r.raw}{('D' if r.double else ' ')}{('C' if r.color else ' ')}{len(members)}{names[0]:19s}{names[1]:19s}"
        if len(line) != 71:
            raise AssertionError('Fixed width overflow')
        fixed.append(line)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'cd_ppm_gsc.csv', index=False)
    (Path(catalog_output) if catalog_output else Path(__file__).resolve().parent / 'cd_ppm_gsc.txt').write_text('\n'.join(fixed) + '\n', encoding='ascii')
    df[df.review.str.contains('ambiguous|conflict|uncertain|unmatched|blend|missing|extrapolation')].to_csv(out / 'review.csv', index=False)
    stats.update(records=len(df), by_type={str(k): int(v) for k, v in df.type.value_counts().items()}, doubles_by_type={str(k): int(v) for k, v in df[df.double].type.value_counts().items()}, designations_ppm=int(df.id1.str.startswith('PPM').sum() + df.id2.str.startswith('PPM').sum()), designations_gsc=int(df.id1.str.startswith('GSC').sum() + df.id2.str.startswith('GSC').sum()), parameters=config)
    (out / 'summary.json').write_text(json.dumps(stats, indent=2, allow_nan=False) + '\n')
    return df

def load_inputs(prepared):
    cd = pd.read_csv(prepared / 'cd.csv', keep_default_na=False, na_values={'mag': ['']})
    modern = pd.read_csv(prepared / 'modern.csv', low_memory=False)
    config = json.loads((prepared / 'calibration.json').read_text())
    return (cd, modern, config)


# CD EVIDENCE

def pair_area_pdf(r, config):
    lo, hi = (config.get('pair_min', 2.0), config['pair_radius'])
    mu, sig = (config['separation_mean'], config['separation_sigma'])
    if not lo <= r <= hi:
        return 0.0
    return float(truncnorm.pdf(r, (lo - mu) / sig, (hi - mu) / sig, loc=mu, scale=sig) / (2 * np.pi * r))

def make_cd_hypotheses(cd, modern, config, return_components=False):
    av = normalize(cd[['x', 'y', 'z']].to_numpy())
    bv = normalize(modern[['x', 'y', 'z']].to_numpy())
    cv, _ = correct(av, config['systematic_coefficients'])
    tree = cKDTree(bv)
    lists = tree.query_ball_point(cv, 2 * np.sin(config['radius_arcsec'] / (2 * ARCSEC_PER_RAD)))
    k = min(64, len(bv))
    dk = tree.query(cv, k=[k])[0][:, 0]
    rk = 2 * np.arcsin(np.clip(dk / 2, 0, 1)) * ARCSEC_PER_RAD
    rho = (k - 1) / (np.pi * np.maximum(rk, 1) ** 2)
    mag, sm = (modern.mag_cd.to_numpy(), modern.sigma_mag.to_numpy())
    var = modern.variability_sigma_cd.to_numpy()
    source = modern.source.to_numpy()
    am = cd.mag.to_numpy()
    active, doubles, colors = (cd.active.to_numpy(), cd.double.to_numpy(), cd.color.to_numpy())
    bg = np.interp(np.nan_to_num(am, nan=9), config['mag_centers'], config['mag_density'])
    ps = config['p_single']
    records, diagnostics, sums = ([], [], [])
    settings = Settings(broad_fraction=config.get('broad_fraction', 0.03))
    for a, js in enumerate(lists):
        if a and a % 30000 == 0:
            print(f'Evidence: {a}/{len(cd)} CD', flush=True)
        if not active[a]:
            records.append([])
            diagnostics.append((0, 0))
            sums.append(0.0)
            continue
        js = sorted(js)
        rr = separation(cv[a], bv[js]) if js else []
        raw = separation(av[a], bv[js]) if js else []
        js_array = np.asarray(js, dtype=int)
        scales = np.where(source[js_array] == 'PPM', config['sigma_circular_ppm'], config['sigma_circular_gsc'])
        rr = np.asarray(rr)
        broad = np.maximum(6 * scales, 60.0)
        density = lambda s: (1 + rr ** 2 / (4 * s * s)) ** (-3) / (2 * np.pi * s * s)
        lr_array = ((1 - settings.broad_fraction) * density(scales) + settings.broad_fraction * density(broad)) / rho[a]
        if np.isfinite(am[a]):
            known = np.isfinite(mag[js_array])
            kk = js_array[known]
            sc = np.sqrt((sm[kk] ** 2 + var[kk] ** 2 + (0.5 if colors[a] else 0) ** 2) / 3)
            phot = t.pdf((am[a] - mag[kk]) / sc, 3) / sc / max(bg[a], 1e-05)
            lr_array[known] *= np.clip(phot, 0.05, 20.0)
        local = {b: (float(lr), float(r)) for b, lr, r in zip(js, lr_array, raw)}
        rec = [((b,), lr * (ps if doubles[a] else 1.0), b, r, math.nan) for b, (lr, r) in local.items()]
        npairs = 0
        if doubles[a]:
            for b, c in combinations(js, 2):
                sep = float(separation(bv[b], bv[c]))
                area = pair_area_pdf(sep, config)
                if area == 0:
                    continue
                lr = (local[b][0] + local[c][0]) * area / rho[a]
                primary = max((b, c), key=lambda j: local[j][0])
                rec.append(((b, c), (1 - ps) * lr, primary, local[primary][1], sep))
                npairs += 1
        records.append(rec)
        diagnostics.append((len(js), npairs))
        sums.append(sum((x[1] for x in rec)))
    scale = max(config['sigma_circular_ppm'], config['sigma_circular_gsc'])
    F = float(spatial_cdf(config['radius_arcsec'], scale, settings))
    sums = np.asarray(sums)
    mask = active & ~doubles
    q = config.get('q')
    if q is None:
        q = float(minimize_scalar(lambda q: -np.log(1 - q * F + q * sums[mask]).sum(), bounds=(0.01, 0.999), method='bounded').x)
        config['q'] = q
    null = 1 - q * F
    hypotheses = []
    for rec in records:
        hs = [Hypothesis((), 0.0)]
        for members, lr, primary, r, sep in rec:
            cost = -math.log(max(q * lr / null, 1e-300))
            if cost <= 0:
                hs.append(Hypothesis(members, cost, primary, r, sep))
        hypotheses.append(hs)
    config['search_mass_approximation'] = F
    if return_components:
        ei = np.repeat(np.arange(len(cd)), [len(js) if active[a] else 0 for a, js in enumerate(lists)])
        ej = np.asarray([b for a, js in enumerate(lists) if active[a] for b in js], dtype=int)
        labels = find_components(len(cd), len(modern), ei, ej)
        return (hypotheses, diagnostics, labels)
    return (hypotheses, diagnostics)


# RUN CD PIPELINE

def identify_variables(modern, variables, out):
    mv = xyz(modern.ra.to_numpy(), modern.dec.to_numpy())
    vv = xyz(variables.ra_deg.to_numpy(), variables.dec_deg.to_numpy())
    ds, ix = cKDTree(mv).query(vv, k=2)
    ds = 2 * np.arcsin(np.clip(ds / 2, 0, 1)) * ARCSEC_PER_RAD
    modern['variable'] = False
    modern['variable_name'] = ''
    modern['variability_sigma_cd'] = 0.0
    rows = []
    for a, v in variables.iterrows():
        b = int(ix[a, 0])
        dist, second = ds[a]
        if dist > 30:
            continue
        accepted = dist <= 5 and second > max(5, 3 * dist)
        amp = float(v.amplitude_full)
        if v.band == 'V':
            poly = [-0.157169, 1.188316, -0.02213]

            def convert(x):
                z = np.clip(x, 3.9, 11.0)
                return np.polynomial.polynomial.polyval(z, poly) + max(0.5, 1.188316 - 0.04426 * z) * (x - z)
            sc = abs(convert(v.mag_faint) - convert(v.mag_bright)) / math.sqrt(6)
            method = 'V-to-CD full amplitude; independent uniform phases'
        else:
            sc = max(1.5, amp / math.sqrt(6))
            method = 'other-band amplitude proxy, CD sigma floor 1.5; sensitivity required'
        if accepted:
            modern.loc[b, 'variable'] = True
            modern.loc[b, 'variable_name'] = v['name']
            modern.loc[b, 'variability_sigma_cd'] = max(sc, modern.loc[b, 'variability_sigma_cd'])
        rows.append(dict(variable=v['name'], modern=modern.name.iloc[b], distance_arcsec=dist, second_distance_arcsec=second, accepted=accepted, band=v.band, amplitude_full=amp, amplitude_truncated=v.amplitude_truncated, sigma_cd=sc, method=method, source=modern.source.iloc[b], coordinate_epoch=v.coordinate_epoch, modern_epoch=modern.epoch.iloc[b], uncertain_extrema=bool(str(v.bright_uncertain).strip() or str(v.faint_uncertain).strip())))
    pd.DataFrame(rows).to_csv(out / 'variable_identifications.csv', index=False)
    return modern

def prepare_cd_match(base, out):
    cd, modern, config = load_inputs(base / 'prepared')
    variables = pd.read_csv(base / 'variables/variables_sur.csv', keep_default_na=False)
    modern = identify_variables(modern, variables, out)
    keep = (modern.source == 'PPM') | modern.variable | modern.mag.isna() | (modern.mag <= 13.5)
    modern = modern.loc[keep].reset_index(drop=True)
    av, bv = (cd[['x', 'y', 'z']].to_numpy(), modern[['x', 'y', 'z']].to_numpy())
    variable_vectors = bv[modern.variable.to_numpy()]
    nearvar = np.zeros(len(cd), bool)
    if len(variable_vectors):
        nearvar = cKDTree(variable_vectors).query(av)[0] < 2 * np.sin(120 / (2 * ARCSEC_PER_RAD))
    clean = cd.active.to_numpy() & ~cd.double.to_numpy() & ~cd.color.to_numpy() & ~nearvar
    photo = pd.read_csv(base / 'prepared/photometric_anchors.csv')
    photo = photo[photo.cd.isin(cd.loc[clean, 'name'])]
    fits = {}
    for band, rows in photo.groupby('band'):
        fit = fit_band(rows.raw_mag.to_numpy(), rows.cd_mag.to_numpy(), rows.validation.to_numpy(), existing=band == 'PPM_V', min_sigma=0.3 if band == 'PPM_V' else 0.5, degree=1 if band.startswith('GSC') else 2)
        if fit:
            fits[band] = fit
    modern['mag_cd'] = np.nan
    for band, fit in fits.items():
        use = modern.band == band
        raw = modern.loc[use, 'mag'].to_numpy()
        lo, hi = fit['mag_training_range']
        clipped = np.clip(raw, lo, hi)
        coef = fit['coefficients']
        derivative = np.polynomial.polynomial.polyder(coef)
        modern.loc[use, 'mag_cd'] = np.polynomial.polynomial.polyval(clipped, coef) + np.maximum(0.5, np.polynomial.polynomial.polyval(clipped, derivative)) * (raw - clipped)
        modern.loc[use, 'sigma_mag'] = fit['sigma_mag'] * np.where((raw < lo) | (raw > hi), 1.5, 1.0)
        modern.loc[use, 'photometry_extrapolated'] = (raw < lo) | (raw > hi)
    config['bands'] = fits
    cmap = {(r.zone, r.num): i for i, r in cd.iterrows() if r.suppl == ' ' and clean[i]}
    ai, bi = ([], [])
    for b, r in modern[modern.source == 'PPM'].iterrows():
        dm = r.dm
        if isinstance(dm, str) and dm[:1] == '-' and dm[1:3].strip().isdigit() and dm[3:8].strip().isdigit():
            key = (-int(dm[1:3]), int(dm[3:8]))
            if -31 <= key[0] <= -23 and key in cmap:
                ai.append(cmap[key])
                bi.append(b)
    ai, bi = (np.asarray(ai), np.asarray(bi))
    hold = np.floor(cd.ra.to_numpy()[ai] / 5).astype(int) % 5 == 0
    res = residuals(av[ai], bv[bi])
    _, _, features = basis(av[ai])
    train = ~hold & (np.linalg.norm(res, axis=1) < 180)
    coef = np.column_stack([robust_linear(features[train], res[train, j]) for j in range(2)])
    before = np.linalg.norm(res[hold], axis=1)
    after = np.linalg.norm(res[hold] - features[hold] @ coef, axis=1)
    accepted = np.median(after) < 0.98 * np.median(before)
    if not accepted:
        coef = np.zeros((6, 2))
    corrected, _ = correct(av, coef)
    rr = separation(corrected[ai], bv[bi])
    sig = float(np.median(rr[train]) / math.sqrt(4 * (math.sqrt(2) - 1)))
    gm = modern.source.to_numpy() == 'GSC'
    gs = calibrate_single(corrected[clean], cd.mag.to_numpy()[clean], bv[gm], modern.mag_cd.to_numpy()[gm], Settings(photometry=False, systematic=False, anchor_radius=90), variable_b=modern.variable.to_numpy()[gm])
    marginal = calibrate_single(av[clean], cd.mag.to_numpy()[clean], bv, modern.mag_cd.to_numpy(), Settings(photometry=False, systematic=False), variable_b=modern.variable.to_numpy())
    config.update(model='robust_cd_v1', radius_arcsec=360.0, sigma_circular_ppm=sig, sigma_circular_gsc=max(sig, gs['sigma']), systematic_coefficients=coef.tolist(), systematic_accepted=bool(accepted), systematic_validation_before=float(np.median(before)), systematic_validation_after=float(np.median(separation(corrected[ai[hold]], bv[bi[hold]]))), ppm_training_anchors=int(train.sum()), ppm_validation_anchors=int(hold.sum()), mag_centers=marginal['mag_centers'], mag_density=marginal['mag_density'], pair_min=2.0, broad_fraction=0.03, variable_count=int(modern.variable.sum()), p_single_source='provisional excess estimator from prepared data; sensitivity required', pair_model='symmetric positional-primary mixture, truncated radial separation / area; no contrast term', variable_model='full amplitude; independent uniform phases; other-band proxy', modern_population=len(modern), modern_mag_limit=13.5)
    cd.to_csv(out / 'cd.csv', index=False)
    modern.to_csv(out / 'modern.csv', index=False)
    (out / 'calibration.json').write_text(json.dumps(config, indent=2) + '\n')
    return (cd, modern, config)


# VALIDATE CD RUN

def validate_cd_match(out, catalog_output, run_controls=False):
    base = Path(__file__).resolve().parent / 'cd_cross'
    out = Path(out)
    cd, modern, config = load_inputs(out)
    result = pd.read_csv(out / 'cd_ppm_gsc.csv', keep_default_na=False)
    raw = Path(catalog_output).read_text().splitlines()
    assert len(raw) == len(cd) == len(result)
    assert all((len(s) == 71 and s[:30] == r for s, r in zip(raw, cd.raw)))
    used = [x for col in ['id1', 'id2'] for x in result[col] if x]
    assert len(used) == len(set(used))
    assert not ((result.type == 2) & ~cd.double).any()
    assert set(used) <= set(modern.name)
    assert not (result.id2.str.startswith('PPM') & ~result.id1.str.startswith('PPM')).any()
    pre = pd.read_csv(base / 'prepared/ppm_gsc.csv')
    assert not pre.ppm.duplicated().any() and (not pre.gsc.duplicated().any())
    assert not set(pre.gsc) & set(used)
    for line, (_, r) in zip(raw, result.iterrows()):
        assert line[30] == ('D' if r.double else ' ') and line[31] == ('C' if r.color else ' ')
        assert int(line[32]) == r.type and line[33:52].strip() == r.id1 and (line[52:71].strip() == r.id2)
    ref = {}
    for _, r in modern[modern.source == 'PPM'].iterrows():
        dm = r.dm
        if isinstance(dm, str) and dm.startswith('-') and dm[1:3].strip().isdigit() and dm[3:8].strip().isdigit():
            key = (-int(dm[1:3]), int(dm[3:8]))
            if -31 <= key[0] <= -23:
                ref.setdefault(key, set()).add(r['name'])
    rows = []
    for a, r in result.iterrows():
        expected = ref.get((r.zone, r.number), set())
        if not expected or cd.suppl.iloc[a] != ' ' or (not cd.active.iloc[a]):
            continue
        actual = {r.id1, r.id2}
        rows.append(dict(cd=r.cd, zone=r.zone, double=r.double, expected=';'.join(sorted(expected)), id1=r.id1, id2=r.id2, any_agreement=bool(expected & actual), all_agreement=expected <= actual, heldout=int(cd.ra.iloc[a] // 5) % 5 == 0, ppm_reference_count=len(expected)))
    reference = pd.DataFrame(rows)
    reference.to_csv(out / 'ppm_reference_comparison.csv', index=False)
    reference[~reference.any_agreement].to_csv(out / 'ppm_reference_disagreements.csv', index=False)
    zone = []
    for z, g in result.groupby('zone'):
        rg = reference[reference.zone == z]
        zone.append(dict(zone=int(z), cd=len(g), unmatched=int((g.type == 0).sum()), with_ppm=int((g.id1.str.startswith('PPM') | g.id2.str.startswith('PPM')).sum()), gsc_only=int(((g.type > 0) & ~g.id1.str.startswith('PPM') & ~g.id2.str.startswith('PPM')).sum()), pairs=int((g.type == 2).sum()), ppm_reference=len(rg), ppm_reference_agree=int(rg.any_agreement.sum())))
    pd.DataFrame(zone).to_csv(out / 'by_zone.csv', index=False)
    report = dict(integrity_passed=True, records=len(cd), reference_rows=len(reference), reference_agreed=int(reference.any_agreement.sum()), reference_fraction=float(reference.any_agreement.mean()), heldout_reference_fraction=float(reference.loc[reference.heldout, 'any_agreement'].mean()), reference_note='DM labels fit training-sector astrometric scale/offset only; not hard constraints. -22 excluded.', global_probabilities_calibrated=False)
    (out / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(report, flush=True)
    if not run_controls:
        return report
    sample = cd.loc[(np.arange(len(cd)) % 100 == 0) | cd.double.to_numpy()].reset_index(drop=True)
    variants = [('baseline', {}), ('p_single_015', {'p_single': 0.15}), ('p_single_075', {'p_single': 0.75}), ('radius_180', {'radius_arcsec': 180.0}), ('variables_no_broadening', {}), ('shift_ra_05deg', {})]
    sigs = None
    sensitivity = []
    for name, changes in variants:
        cfg = copy.deepcopy(config)
        cfg.update(changes)
        a = sample.copy()
        m = modern
        if name == 'variables_no_broadening':
            m = modern.copy()
            m['variability_sigma_cd'] = 0.0
        if name == 'shift_ra_05deg':
            a[['x', 'y', 'z']] = xyz((a.ra + 0.5) % 360, a.dec)
        hs, diag, labels = make_cd_hypotheses(a, m, cfg, return_components=True)
        chosen, stats = solve_hypotheses(hs, len(m), 600, component_labels=labels, double_flags=a.double.to_numpy())
        sig = [hs[i][chosen[i]].members for i in range(len(a))]
        if sigs is None:
            sigs = sig
        row = dict(variant=name, sample=len(a), assigned=sum((bool(x) for x in sig)), pairs=sum((len(x) == 2 for x in sig)), changed=sum((x != y for x, y in zip(sig, sigs))))
        sensitivity.append(row)
        print(row, flush=True)
        pd.DataFrame(sensitivity).to_csv(out / 'sensitivity.csv', index=False)


def export_cd_catalog(results_csv, output):
    """Export 71 ASCII bytes per row, preserving the original 30 CD bytes."""
    table = pd.read_csv(results_csv, keep_default_na=False)
    lines = []
    for row in table.itertuples(index=False):
        if len(row.original) != 30 or len(row.id1) > 19 or len(row.id2) > 19:
            raise ValueError('Catalog field width exceeded')
        if row.type not in (0, 1, 2) or bool(row.id1) != (row.type > 0) or bool(row.id2) != (row.type == 2):
            raise ValueError('Multiplicity and identifiers disagree')
        line = f'{row.original}{"D" if row.double else " "}{"C" if row.color else " "}{row.type}{row.id1:19s}{row.id2:19s}'
        if len(line.encode('ascii')) != 71:
            raise ValueError('Invalid fixed-width record')
        lines.append(line)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('\n'.join(lines)+'\n', encoding='ascii')
    return len(lines)