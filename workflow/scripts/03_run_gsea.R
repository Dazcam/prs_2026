#--------------------------------------------------------------------------------------
#
#   Gene Set Enrichment Analysis using fgsea
#   Run per cell type on limma results
#
#   Ranking method (toggled via snakemake params):
#     - t_statistic:    signed t-statistic (direction + magnitude) [used]
#     - neg_log10_pval: signed -log10(p-value) (direction + significance)
#
#   Gene sets tested:
#     - GO Biological Process (MSigDB C5 GO:BP)
#
#   Gene sets converted from symbols to Ensembl IDs via lookup table.
#   Optionally restricted to protein-coding genes.
#
#--------------------------------------------------------------------------------------

if (exists("snakemake") && length(snakemake@log) != 0) {
  dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, append = TRUE)
  sink(log_con, append = TRUE, type = "message")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
})

limma_file    <- snakemake@input[["limma"]]
lookup_file   <- snakemake@input[["lookup"]]
out_res       <- snakemake@output[["results"]]
out_rds       <- snakemake@output[["rds"]]
cell_type     <- snakemake@params[["cell_type"]]
filter_coding <- snakemake@params[["filter_coding"]]
min_gs_size   <- snakemake@params[["min_gs_size"]]
max_gs_size   <- snakemake@params[["max_gs_size"]]
rank_method   <- snakemake@params[["rank_method"]]

# Validate rank_method
if (!rank_method %in% c("t_statistic", "neg_log10_pval")) {
  stop(paste("Invalid rank_method:", rank_method,
             "— must be 't_statistic' or 'neg_log10_pval'"))
}

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush.console()
}

log_msg("Cell type:          ", cell_type)
log_msg("Rank method:        ", rank_method)
log_msg("Filter coding only: ", filter_coding)
log_msg("Min gene set size:  ", min_gs_size)
log_msg("Max gene set size:  ", max_gs_size)

#--------------------------------------------------------------------------------------
# Load limma results
#--------------------------------------------------------------------------------------

log_msg("Loading limma results from ", limma_file)
limma <- read_tsv(limma_file, show_col_types = FALSE)
log_msg("Genes in limma results: ", nrow(limma))

#--------------------------------------------------------------------------------------
# Load gene lookup and optionally filter to protein-coding genes
#--------------------------------------------------------------------------------------

log_msg("Loading gene lookup from ", lookup_file)
lookup <- read_tsv(lookup_file, show_col_types = FALSE)
log_msg("Genes in lookup: ", nrow(lookup))

if (filter_coding) {
  coding_genes <- lookup |>
    filter(biotype == "protein_coding") |>
    pull(gene)
  limma <- limma |> filter(gene %in% coding_genes)
  log_msg("Protein-coding genes retained: ", nrow(limma))
} else {
  log_msg("No biotype filter applied — using all genes")
}

#--------------------------------------------------------------------------------------
# Build symbol -> Ensembl mapping for gene set conversion
#--------------------------------------------------------------------------------------

symbol_to_ensembl <- lookup |>
  filter(!is.na(symbol), !is.na(gene), symbol != "") |>
  select(symbol, gene) |>
  distinct()

log_msg("Symbol-to-Ensembl pairs for gene set conversion: ", nrow(symbol_to_ensembl))

#--------------------------------------------------------------------------------------
# Build gene ranking
#--------------------------------------------------------------------------------------

log_msg("Building gene ranking using method: ", rank_method)

ranks <- switch(rank_method,
  t_statistic = limma |>
    mutate(rank = t_stat) |>
    arrange(desc(rank)) |>
    select(gene, rank),
  neg_log10_pval = limma |>
    mutate(rank = sign(beta) * -log10(p_value)) |>
    arrange(desc(rank)) |>
    select(gene, rank)
)

rank_vec        <- ranks$rank
names(rank_vec) <- ranks$gene
rank_vec        <- sort(rank_vec, decreasing = TRUE)

log_msg("Rank vector length: ", length(rank_vec))
log_msg("Rank range: ", round(min(rank_vec), 3), " to ", round(max(rank_vec), 3))

#--------------------------------------------------------------------------------------
# Load gene sets from MSigDB and convert symbols to Ensembl IDs
#--------------------------------------------------------------------------------------

log_msg("Loading gene sets from MSigDB")

convert_gs_to_ensembl <- function(gs_raw, symbol_map) {
  gs_raw |>
    select(gs_name, gene_symbol) |>
    inner_join(symbol_map, by = c("gene_symbol" = "symbol"),
               relationship = "many-to-many") |>
    group_by(gs_name) |>
    summarise(genes = list(unique(gene)), .groups = "drop")
}

gs_gobp_raw <- msigdbr(species = "Homo sapiens",
                        collection = "C5", subcollection = "GO:BP")

gs_gobp     <- convert_gs_to_ensembl(gs_gobp_raw, symbol_to_ensembl)

make_gs_list <- function(gs_df) {
  setNames(gs_df$genes, gs_df$gs_name)
}

gs_list <- list(
  GO_BP     = make_gs_list(gs_gobp)
)

log_msg("Gene set counts after Ensembl conversion:")
for (gs_name in names(gs_list)) {
  log_msg("  ", gs_name, ": ", length(gs_list[[gs_name]]), " gene sets")
}

#--------------------------------------------------------------------------------------
# Run fgsea
#--------------------------------------------------------------------------------------

run_fgsea <- function(rank_vec, gene_sets, gs_category,
                      ct, symbol_map, min_size, max_size) {

  log_msg("Running fgsea: ", gs_category)

  set.seed(42)
  res <- fgsea(
    pathways    = gene_sets,
    stats       = rank_vec,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = 10000
  )

  res |>
    as_tibble() |>
    arrange(pval) |>
    mutate(
      gs_category = gs_category,
      rank_method = rank_method,
      cell_type   = ct,
      leadingEdge_ensembl = sapply(leadingEdge, paste, collapse = ","),
      leadingEdge_symbol  = sapply(leadingEdge, function(ids) {
        syms <- symbol_map$symbol[match(ids, symbol_map$gene)]
        syms[is.na(syms)] <- ids[is.na(syms)]
        paste(syms, collapse = ",")
      })
    ) |>
    rename(
      gene_set   = pathway,
      p_value    = pval,
      fdr        = padj,
      enrichment = NES
    ) |>
    select(cell_type, gs_category, rank_method, gene_set,
           enrichment, p_value, fdr, size,
           leadingEdge_ensembl, leadingEdge_symbol)
}

results_list <- list()

for (gs_cat in names(gs_list)) {
  results_list[[gs_cat]] <- run_fgsea(
    rank_vec, gs_list[[gs_cat]], gs_cat,
    cell_type, symbol_to_ensembl, min_gs_size, max_gs_size
  )
}

results <- bind_rows(results_list)

#--------------------------------------------------------------------------------------
# QC summary
#--------------------------------------------------------------------------------------

log_msg("GSEA QC summary:")
results |>
  group_by(gs_category) |>
  summarise(
    n_tested   = n(),
    n_fdr05    = sum(fdr < 0.05, na.rm = TRUE),
    n_enriched = sum(fdr < 0.05 & enrichment > 0, na.rm = TRUE),
    n_depleted = sum(fdr < 0.05 & enrichment < 0, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  print()

#--------------------------------------------------------------------------------------
# Output
#--------------------------------------------------------------------------------------

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)

log_msg("Writing results to ", out_res)
write_tsv(results, out_res)

log_msg("Saving RDS to ", out_rds)
saveRDS(results_list, out_rds)

log_msg("Done: ", cell_type)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
