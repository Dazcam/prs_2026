#!/usr/bin/env python3
"""
Combine genotype PCs, sample annotations, and PRS into a covariate matrix.

Usage:
    python 02_prepare_covariates.py <pcs> <annot> <prs> <output> [n_pcs]
"""

import sys
import argparse
import pandas as pd

sys.path.insert(0, 'scripts')
from lib.covariates import encode_sex, select_pcs, join_covariates, validate_covariates


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pcs',    help='PLINK eigenvec file')
    parser.add_argument('annot',  help='Sample annotations Excel file')
    parser.add_argument('prs',    help='Standardised PRS TSV')
    parser.add_argument('output', help='Output covariate TSV')
    parser.add_argument('--n_pcs', type=int, default=5,
                        help='Number of PCs to retain (default: 5)')
    return parser.parse_args()


def main(args):

    # ---- Load PCs ----
    print(f"Loading genotype PCs from {args.pcs}")
    pcs_raw = pd.read_csv(args.pcs, sep=r'\s+', header=None)
    pcs = select_pcs(pcs_raw, n_pcs=args.n_pcs)
    print(f"  {len(pcs)} samples, {args.n_pcs} PCs retained")

    # ---- Load annotations ----
    print(f"Loading annotations from {args.annot}")
    annot_raw = pd.read_excel(args.annot)
    print(f"  Columns: {list(annot_raw.columns)}")

    annot = annot_raw[['sample ID', 'sex', 'PCW']].copy()
    annot.columns = ['sample', 'sex', 'PCW']
    annot['sample'] = annot['sample'].astype(str)
    annot['sex']    = encode_sex(annot['sex'])

    n_missing_sex = annot['sex'].isna().sum()
    if n_missing_sex > 0:
        print(f"  WARNING: {n_missing_sex} samples with unrecognised sex coding")

    print(f"  {len(annot)} samples in annotations")

    # ---- Load PRS ----
    print(f"Loading PRS from {args.prs}")
    prs = pd.read_csv(args.prs, sep='\t')
    prs['sample'] = prs['sample'].astype(str)
    print(f"  {len(prs)} samples with PRS")

    # ---- Join ----
    print("Joining covariates...")
    covs = join_covariates(prs, pcs, annot)

    summary = validate_covariates(covs)
    print(f"  Samples after join: {summary['n_samples']}")
    print(f"  Columns: {summary['columns']}")
    if summary['missing_values']:
        print(f"  WARNING: Missing values: {summary['missing_values']}")

    covs.to_csv(args.output, sep='\t', index=False)
    print(f"Written to {args.output}")


if __name__ == '__main__':
    main(parse_args())
