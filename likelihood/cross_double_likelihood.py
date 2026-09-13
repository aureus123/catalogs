#!/usr/bin/env python3
"""CD singles/doubles matching with SCIP. See CD_CROSS.md for methods and formats."""
import sys
sys.dont_write_bytecode = True
import argparse
import gzip
import hashlib
import json
import pickle
import time
from pathlib import Path
from importlib.metadata import version
import pandas as pd
from likelihood_common import (
    Hypothesis, load_inputs, prepare_catalogs, prepare_cd_match,
    make_cd_hypotheses, solve_hypotheses, write_result, validate_cd_match,
    export_cd_catalog,
)
# Hypothesis is also exposed here to read the original locally generated checkpoint.


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, default=root/'cd_cross/results')
    ap.add_argument('--output', type=Path, default=root/'cd_ppm_gsc.txt')
    ap.add_argument('--prepare', action='store_true', help='Regenerate PPM/GSC union from repository source catalogs')
    ap.add_argument('--rebuild', action='store_true', help='Recalibrate and rebuild hypotheses from cached source tables')
    ap.add_argument('--build-only', action='store_true', help='Save preparation and hypotheses without solving')
    ap.add_argument('--solve-only', action='store_true', help='Solve the trusted saved checkpoint (default)')
    ap.add_argument('--validate', action='store_true', help='Check existing output and compare PPM CD designations')
    ap.add_argument('--sensitivity', action='store_true', help='Validate output and run frozen-model sensitivity controls')
    ap.add_argument('--export-only', action='store_true', help='Export the existing results CSV in the documented fixed-width format')
    ap.add_argument('--time-limit', type=float, default=600., help='SCIP time limit per component; optimality remains mandatory')
    args = ap.parse_args()
    if args.time_limit <= 0:
        ap.error('--time-limit must be positive')
    if args.solve_only and (args.prepare or args.rebuild or args.build_only):
        ap.error('--solve-only cannot be combined with preparation options')
    args.data.mkdir(parents=True, exist_ok=True)
    if args.export_only:
        export_cd_catalog(args.data/'cd_ppm_gsc.csv', args.output)
        return
    if args.validate or args.sensitivity:
        validate_cd_match(args.data, args.output, run_controls=args.sensitivity)
        return
    start = time.monotonic()
    if args.prepare:
        prepare_catalogs(root.parent, root/'cd_cross/prepared')
    rebuilding = args.prepare or args.rebuild or args.build_only
    if rebuilding:
        cd, modern, config = prepare_cd_match(root/'cd_cross', args.data)
        hs, diag, labels = make_cd_hypotheses(cd, modern, config, return_components=True)
        with gzip.open(args.data/'hypotheses.pkl.gz', 'wb') as f:
            pickle.dump((hs, diag, labels, config), f, protocol=5)
        manifest = {name:hashlib.sha256((args.data/name).read_bytes()).hexdigest()
                    for name in ['cd.csv', 'modern.csv']}
        (args.data/'checkpoint_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        (args.data/'calibration.json').write_text(json.dumps(config, indent=2)+'\n')
    else:
        cd, modern, config = load_inputs(args.data)
        manifest = json.loads((args.data/'checkpoint_manifest.json').read_text())
        for name, digest in manifest.items():
            if hashlib.sha256((args.data/name).read_bytes()).hexdigest() != digest:
                raise ValueError(f'Checkpoint input changed: {name}')
        # This is a trusted local checkpoint, never a downloaded/untrusted pickle.
        with gzip.open(args.data/'hypotheses.pkl.gz', 'rb') as f:
            hs, diag, labels, config = pickle.load(f)
    counts = pd.Series(labels[:len(cd)]).value_counts()
    print(f'{len(cd)} CD; {len(counts)} components; {sum(map(len,hs))} hypotheses', flush=True)
    if args.build_only:
        return
    selected, stats = solve_hypotheses(hs, len(modern), args.time_limit,
                                      component_labels=labels, double_flags=cd.double.to_numpy())
    stats.update(solver='SCIP', pyscipopt_version=version('pyscipopt'),
                 geometric_components=len(counts), largest_component_cd=int(counts.max()),
                 component_radius_arcsec=config['radius_arcsec'], seconds=time.monotonic()-start)
    write_result(cd, modern, hs, selected, diag, args.data, config, stats, catalog_output=args.output)
    (args.data/'calibration.json').write_text(json.dumps(config, indent=2)+'\n')
    print(json.dumps({k:v for k,v in stats.items() if k not in ['parameters','scip']}, indent=2))


if __name__ == '__main__':
    main()
