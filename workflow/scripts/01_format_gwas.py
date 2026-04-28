#!/usr/bin/env python3
"""
Reformat PGC3 SCZ GWAS summary statistics for PRS-CS input.
PRS-CS requires: SNP A1 A2 BETA P
- SNP: rsID
- A1: effect allele
- A2: other allele
- BETA: log-OR (used only for direction — p-value drives shrinkage)
- P: p-value
"""

import sys
import pandas as pd

input_file  = sys.argv[1]   # scz_hg38_with_neff.tsv
output_file = sys.argv[2]   # scz_prscs_input.txt

print(f"Reading {input_file}...")
df = pd.read_csv(input_file, sep='\t', low_memory=False)

print(f"Input variants: {len(df)}")
print(f"Columns: {list(df.columns)}")

# ---- QC filters ----

# Remove non-rsIDs
df = df[df['SNP'].str.startswith('rs', na=False)]
print(f"After rsID filter: {len(df)}")

# Remove low imputation quality
df = df[df['IMPINFO'] >= 0.9]
print(f"After IMPINFO >= 0.9 filter: {len(df)}")

# Remove alt chromosomes
df = df[~df['CHR'].astype(str).str.endswith('alt')]
df = df[df['CHR'].astype(str).isin([str(i) for i in range(1, 23)])]
print(f"After chromosome filter: {len(df)}")

# Remove NA in required columns
df = df.dropna(subset=['SNP', 'A1', 'A2', 'BETA', 'PVAL'])
print(f"After NA filter: {len(df)}")

# Remove invalid beta/SE
df = df[df['SE'] > 0]
df = df[df['BETA'].apply(lambda x: abs(x) < 5)]
print(f"After beta/SE filter: {len(df)}")

# Remove duplicates
df = df.drop_duplicates(subset='SNP', keep='first')
print(f"After dedup: {len(df)}")

# ---- Format output ----
out = pd.DataFrame({
    'SNP':  df['SNP'],
    'A1':   df['A1'].str.upper(),
    'A2':   df['A2'].str.upper(),
    'BETA': df['BETA'],
    'P':    df['PVAL']
})

out.to_csv(output_file, sep='\t', index=False)
print(f"Written {len(out)} variants to {output_file}")
