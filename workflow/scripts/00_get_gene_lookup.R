#--------------------------------------------------------------------------------------
#
#   Fetch Ensembl gene ID to symbol mapping via biomaRt
#   Run once and save to resources/sheets/gene_lookup_hg38.tsv
#
#   Output columns:
#     ensembl_gene_id, hgnc_symbol, chromosome_name,
#     start_position, end_position, gene_biotype
#
#--------------------------------------------------------------------------------------

args     <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]
log_file <- args[2]

# Open log connection directly in R — avoids shell redirect issues in containers
log_con <- file(log_file, open = "wt")
sink(log_con, append = TRUE)
sink(log_con, append = TRUE, type = "message")

suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
})

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush.console()
}

# Create output directory if needed
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Connecting to Ensembl biomaRt (GRCh38)")

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast" # fallback mirrors: "uswest", "asia", "www"
)

log_msg("Fetching gene annotations")

lookup <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype"
  ),
  mart = mart
)

log_msg("Fetched ", nrow(lookup), " records")
head(lookup)

# Basic cleaning
lookup <- lookup |>
  as_tibble() |>
  filter(chromosome_name %in% c(1:22)) |>             # restrict to std chrs
  filter(gene_biotype == "protein_coding") |>	      # restrict to protien coding genes	  
  filter(hgnc_symbol != "") |>                        # restrict to hgnc symbol
  distinct(ensembl_gene_id, .keep_all = TRUE) |>      # retain distinct genes (dedup)
  rename(
    gene       = ensembl_gene_id,
    symbol     = hgnc_symbol,
    chr        = chromosome_name,
    start      = start_position,
    end        = end_position,
    biotype    = gene_biotype
  )

log_msg("NA values after cleaning: ", sum(is.na(lookup)))
log_msg("Duplicate rows after cleaning: ", sum(duplicated(lookup)))
log_msg("Duplicate Ensembl IDs: ", sum(duplicated(lookup$gene)))
log_msg("Duplicate gene symbols: ", sum(duplicated(lookup$symbol)))
log_msg("After cleaning: ", nrow(lookup), " unique Ensembl IDs with symbols")
log_msg("Biotypes represented: ", paste(unique(lookup$biotype), collapse = ", "))

write_tsv(lookup, out_file)
log_msg("Written to ", out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
