#!/usr/bin/env python3
"""
Prepare covariate matrix for PRS-expression analysis.

Extracts sex, PCW, and genotype PCs from the eQTL covariate file,
joins with PRS scores, and restricts to European samples from fam file.

Expression PCs are deliberately excluded — they would remove the
PRS-expression signal of interest.

Inputs:
  1: eQTL covariate file (covariates x samples format)
  2: PRS scores TSV (sample, PRS, PRS_std)
  3: PLINK fam file (defines target European samples)
  4: output covariate TSV (samples x covariates)
"""

import sys
import argparse
import pandas as pd

sys.path.insert(0, 'scripts')
from lib.covariates import join_covariates, check_minimum_samples, validate_covariates


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('covariates',   help='eQTL covariate file (covariates x samples)')
    parser.add_argument('prs',          help='PRS scores TSV')
    parser.add_argument('fam',          help='PLINK fam file')
    parser.add_argument('output',       help='Output covariate TSV')
    parser.add_argument('--keep_covs',  nargs='+',
                        default=['Sex', 'PCW', 'genPC1', 'genPC2', 'genPC3', 'genPC4'],
                        help='Covariate rows to retain (default: Sex PCW genPC1-4)')
    parser.add_argument('--min_samples', type=int, default=10)
    return parser.parse_args()


def load_covariates(covs_file: str, keep_covs: list) -> pd.DataFrame:
    """
    Load eQTL covariate file (covariates x samples) and transpose
    to samples x covariates format, keeping only requested covariates.

    Parameters
    ----------
    covs_file : str
        Path to covariate file with covariates as rows, samples as columns
    keep_covs : list
        Covariate names to retain

    Returns
    -------
    pd.DataFrame
        Samples x covariates DataFrame with sample column
    """
    covs = pd.read_csv(covs_file, sep='\t', index_col=0)

    print(f"  Covariates available: {list(covs.index)}")
    print(f"  Samples in covariate file: {covs.shape[1]}")

    # Check requested covariates exist
    missing = [c for c in keep_covs if c not in covs.index]
    if missing:
        raise ValueError(f"Requested covariates not found in file: {missing}")

    # Subset to requested covariates and transpose
    covs_sub = covs.loc[keep_covs].T
    covs_sub.index.name = 'sample'
    covs_sub = covs_sub.reset_index()
    covs_sub['sample'] = covs_sub['sample'].astype(str)

    return covs_sub


def load_fam_samples(fam_file: str) -> list:
    """Extract IID column from PLINK fam file as sample list."""
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    return fam['IID'].astype(str).tolist()


def main(args):

    # ---- Load fam sample list ----
    print(f"Loading target samples from {args.fam}")
    euro_samples = load_fam_samples(args.fam)
    print(f"  {len(euro_samples)} European samples")

    # ---- Load and filter covariates ----
    print(f"Loading covariates from {args.covariates}")
    print(f"  Retaining: {args.keep_covs}")
    covs = load_covariates(args.covariates, args.keep_covs)

    # Restrict to European samples
    covs = covs[covs['sample'].isin(euro_samples)]
    print(f"  Samples after Euro filter: {len(covs)}")

    # ---- Load PRS ----
    print(f"Loading PRS from {args.prs}")
    prs = pd.read_csv(args.prs, sep='\t')
    prs['sample'] = prs['sample'].astype(str)
    print(f"  {len(prs)} samples with PRS")

    # ---- Join PRS to covariates ----
    print("Joining PRS to covariates...")
    combined = covs.merge(prs[['sample', 'PRS', 'PRS_std']], on='sample', how='inner')

    check_minimum_samples(combined, min_samples=args.min_samples)

    summary = validate_covariates(combined)
    print(f"  Final samples: {summary['n_samples']}")
    print(f"  Columns: {summary['columns']}")
    if summary['missing_values']:
        print(f"  WARNING missing values: {summary['missing_values']}")

    combined.to_csv(args.output, sep='\t', index=False)
    print(f"Written to {args.output}")


if __name__ == '__main__':
    main(parse_args())
