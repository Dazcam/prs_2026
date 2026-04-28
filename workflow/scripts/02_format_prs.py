#!/usr/bin/env python3
"""
Reformat PLINK .profile output to a clean sample x PRS table.
Standardises PRS to mean 0 SD 1 for interpretable effect sizes.

Input:  PLINK .profile file (FID IID PHENO CNT CNT2 SCORE)
Output: TSV with sample, PRS, PRS_std columns
"""

import sys
import argparse
import pandas as pd

sys.path.insert(0, 'scripts')
from lib.prs import standardise_prs, validate_prs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',  help='PLINK .profile file')
    parser.add_argument('output', help='Output standardised PRS TSV')
    return parser.parse_args()


def main(args):
    print(f"Reading PLINK profile from {args.input}")
    raw = pd.read_csv(args.input, sep=r'\s+')
    print(f"  Columns: {list(raw.columns)}")
    print(f"  Samples: {len(raw)}")

    # Rename to standard format
    prs = pd.DataFrame({
        'sample': raw['IID'].astype(str),
        'PRS':    raw['SCORE']
    })

    print(f"Standardising PRS")
    prs_std = standardise_prs(prs)

    summary = validate_prs(prs_std)
    for k, v in summary.items():
        print(f"  {k}: {v}")

    prs_std.to_csv(args.output, sep='\t', index=False)
    print(f"Written to {args.output}")


if __name__ == '__main__':
    main(parse_args())
