#!/usr/bin/env python3
"""
Restrict pseudobulk expression matrix to European samples
defined by the genotype fam file from stage 1.

Input is TensorQTL format:
  #Chr  start  end  TargetID  sample1  sample2  ...

Inputs:
  1: expression file (TensorQTL bed format, optionally gzipped)
  2: PLINK fam file (defines which samples to keep)
  3: output filtered expression TSV
"""

import sys
import argparse
import pandas as pd

sys.path.insert(0, 'scripts')
from lib.expression import (filter_to_common_samples,
                             check_minimum_samples,
                             validate_expression_filter,
                             META_COLS)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('expression', help='TensorQTL format expression bed(.gz)')
    parser.add_argument('fam',        help='PLINK fam file defining target samples')
    parser.add_argument('output',     help='Output filtered expression TSV')
    parser.add_argument('--min_samples', type=int, default=10,
                        help='Minimum samples required (default: 10)')
    return parser.parse_args()


def load_fam_samples(fam_file: str) -> list:
    """
    Extract sample IDs from PLINK fam file.
    Fam file columns: FID IID PAT MAT SEX PHENO
    We use IID (column 2) as the sample identifier.

    Parameters
    ----------
    fam_file : str
        Path to PLINK fam file

    Returns
    -------
    list
        Sample IDs as strings
    """
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    samples = fam['IID'].astype(str).tolist()
    print(f"  {len(samples)} samples in fam file")
    print(f"  First 5: {samples[:5]}")
    return samples


def main(args):

    # ---- Load target samples from fam file ----
    print(f"Loading target samples from {args.fam}")
    keep_samples = load_fam_samples(args.fam)

    # ---- Load expression ----
    print(f"Loading expression from {args.expression}")
    compression = 'gzip' if args.expression.endswith('.gz') else None
    expr = pd.read_csv(args.expression, sep='\t', compression=compression)

    # Handle TensorQTL header which may start with #Chr
    print(f"  Genes: {len(expr)}")
    print(f"  Samples in expression: {len(expr.columns) - len(META_COLS)}")

    # ---- Filter to common samples ----
    print("Filtering to European samples from fam file...")
    expr_filtered = filter_to_common_samples(expr, keep_samples)

    # ---- Production minimum sample guard ----
    n_kept = len([c for c in expr_filtered.columns if c not in META_COLS])
    check_minimum_samples(n_kept, min_samples=args.min_samples)

    # ---- Summary ----
    summary = validate_expression_filter(expr, expr_filtered, keep_samples)
    print(f"  Genes retained:          {summary['n_genes']}")
    print(f"  Samples original:        {summary['n_samples_original']}")
    print(f"  Samples kept:            {summary['n_samples_kept']}")
    print(f"  Samples dropped:         {summary['n_samples_dropped']}")
    print(f"  In fam not in expr:      {summary['n_missing_from_expr']}")
    if summary['missing_from_expr']:
        print(f"  Missing: {summary['missing_from_expr']}")

    expr_filtered.to_csv(args.output, sep='\t', index=False)
    print(f"Written to {args.output}")


if __name__ == '__main__':
    main(parse_args())
