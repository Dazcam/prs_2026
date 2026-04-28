# SCZ Polygenic Risk Score | Developing Human Brain

> **Polygenic risk score analysis of schizophrenia GWAS in cell types of the developing human brain**  
> Cardiff University · Division of Psychological Medicine and Clinical Neurosciences  
> 2026

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0-brightgreen)](https://snakemake.readthedocs.io)
[![Platform](https://img.shields.io/badge/platform-SLURM%20HPC-blue)](https://slurm.schedmd.com/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## What this is

This repository computes polygenic risk scores (PRS) for schizophrenia in 118 European samples with prenatal brain single-nucleus RNA sequencing data, and tests whether genetic risk load correlates with cell-type-specific gene expression across 19 brain cell types.

PRS are derived from the PGC3 schizophrenia GWAS (European subset) using PRS-CS-auto. Association testing between PRS and pseudobulk gene expression is performed per cell type using linear models. Gene set enrichment analysis is then used to identify biological pathways enriched among PRS-associated genes.

## Pipeline overview

The analysis runs as three sequential Snakemake pipelines on a SLURM HPC cluster.

| # | Stage | Tools |
| :---: | :--- | :--- |
| 01 | PRS computation | PRS-CS-auto, PLINK |
| 02 | PRS-expression association | limma |
| 03 | Gene set enrichment analysis | fgsea, MSigDB |

## Repository structure

```bash
.
├── config/                   # Pipeline configuration and SLURM profile
├── resources/
│   ├── containers/           # Singularity container definitions
│   ├── refs/                 # PRS-CS LD reference panels (EUR UKBB)
│   ├── sheets/               # Sample annotations and gene lookup tables
│   └── software/             # PRS-CS source code
├── workflow/
│   ├── Snakefile             # Master workflow entry point
│   ├── rules/                # Snakemake rule files (one per pipeline stage)
│   ├── scripts/              # R and Python analysis scripts
│   │   └── lib/              # Reusable Python modules
│   ├── tests/                # Unit tests for Python modules
│   ├── envs/                 # Conda environment definitions
│   └── snakemake.sh          # SLURM submission wrapper
└── results/                  # Pipeline outputs (not tracked in git)
```

## Inputs

| Input | Source |
| :--- | :--- |
| Genotypes (VCF, hg38) | Prenatal brain cohort (118 European samples) |
| SCZ GWAS summary statistics | PGC3 European subset (`PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz`) |
| PRS-CS LD reference | UK Biobank EUR (`ldblk_ukbb_eur`) |
| Pseudobulk expression | TensorQTL format, quantile normalised, 19 cell types |
| Covariates | PCW, sex, genotype PCs |

## Quickstart

### Prerequisites

| Requirement | Version |
| :--- | :--- |
| Snakemake | ≥ 9.0 |
| Conda / Mamba | any |
| Singularity | ≥ 3.x |
| SLURM | any |
| Python | ≥ 3.12 |

### Setup

```bash
# Clone the repo
git clone https://github.com/Dazcam/prs_2026.git
cd prs_2026/workflow

# Download PRS-CS
git clone https://github.com/getian107/PRScs.git ../resources/software/PRS-CS

# Download LD reference panel (~6.25 GB)
mkdir -p ../resources/refs
wget -P ../resources/refs \
    https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz
tar -zxvf ../resources/refs/ldblk_ukbb_eur.tar.gz -C ../resources/refs/
rm ../resources/refs/ldblk_ukbb_eur.tar.gz

# Download HapMap3+ gene lookup (run once, requires internet)
Rscript scripts/00_get_gene_lookup.R \
    ../resources/sheets/gene_lookup_hg38.tsv
```

### Run

```bash
# Dry run — check the full job graph without executing
snakemake -np

# Submit to SLURM
./snakemake.sh
```

### Tests

```bash
cd workflow
python3 -m pytest tests/ -v
```

## Data access

| Resource | Location |
| :--- | :--- |
| Genotype data | Available on request pending publication |
| Pseudobulk expression | Available on request pending publication |
| SCZ GWAS | [PGC downloads](http://www.med.unc.edu/pgc/results-and-downloads) |

## Citation

> Manuscript in preparation.

## Funding

> This work was funded by Medical Research Council / UK Research and Innovation
> project grant ([MR/Y003756/1](https://gtr.ukri.org/projects?ref=MR%2FY003756%2F1)).

## Licence & copyright

Code in this repository is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).  
© 2026 Cardiff University. See [LICENCE.md](LICENCE.md) for details.
