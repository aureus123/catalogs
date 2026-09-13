#!/usr/bin/env python3
"""Robust one-to-one historical-catalog matching; see README.md.

Old command syntax remains valid: A.csv B.csv output.csv.
Optional variable-star CSVs use catalog identifiers and per-epoch magnitude scatter.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import json
from pathlib import Path
from likelihood_common import Settings,load_catalog as load,load_variability,run_single as run,save_single as save


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('catalog_a',type=Path); ap.add_argument('catalog_b',type=Path); ap.add_argument('output',type=Path)
    ap.add_argument('--radius',type=float,default=None)
    ap.add_argument('--anchor-radius',type=float,default=None); ap.add_argument('--no-systematic',action='store_true')
    ap.add_argument('--no-photometry',action='store_true'); ap.add_argument('--q',type=float)
    ap.add_argument('--broad-fraction',type=float,default=None); ap.add_argument('--model-in',type=Path)
    ap.add_argument('--variables-a',type=Path); ap.add_argument('--variables-b',type=Path)
    ap.add_argument('--variable-sigma',type=float,help='Per-epoch variability standard deviation in mag; unknown by default, hence position-only')
    ap.add_argument('--zero-is-valid',action='store_true'); args=ap.parse_args()
    model=json.loads(args.model_in.read_text()) if args.model_in else None
    cfg=Settings(**model['settings']) if model and 'settings' in model else Settings()
    for key,value in [('radius',args.radius),('anchor_radius',args.anchor_radius),('q',args.q),('broad_fraction',args.broad_fraction),('variable_sigma',args.variable_sigma)]:
        if value is not None: setattr(cfg,key,value)
    if args.no_systematic: cfg.systematic=False
    if args.no_photometry: cfg.photometry=False
    if args.zero_is_valid: cfg.zero_missing=False
    an,av,am=load(args.catalog_a,cfg.zero_missing); bn,bv,bm=load(args.catalog_b,cfg.zero_missing)
    va,sa=load_variability(args.catalog_a,an,cfg.variable_sigma,args.variables_a)
    vb,sb=load_variability(args.catalog_b,bn,cfg.variable_sigma,args.variables_b)
    result=run(av,am,bv,bm,cfg,model,va,vb,sa,sb); df=save(args.output,an,av,am,bn,bv,result)
    print(json.dumps(dict(records=len(df),matched=len(result[0]),sigma=result[-1]['sigma'],q=result[-1]['q'],
                          systematic=result[-1]['systematic_accepted']),indent=2))

if __name__=='__main__': main()
