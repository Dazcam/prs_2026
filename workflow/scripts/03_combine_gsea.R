#--------------------------------------------------------------------------------------
#
#   Combine GSEA results across all cell types
#   Adds global FDR correction across all cell types and gene set categories
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
})

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush.console()
}

input_files <- unlist(snakemake@input)
out_file    <- snakemake@output[[1]]

log_msg("Combining GSEA results across ", length(input_files), " cell types")

results <- map_dfr(input_files, function(f) {
  log_msg("  Reading: ", f)
  read_tsv(f, show_col_types = FALSE)
})

log_msg("Total gene set tests: ", nrow(results))
log_msg("Cell types: ",           n_distinct(results$cell_type))

# Global FDR across all cell types and gene set categories
results <- results |>
  mutate(fdr_global = p.adjust(p_value, method = "BH"))

# Summary per category and ranking method
log_msg("Summary of significant results (FDR < 0.05):")
results |>
  group_by(cell_type, gs_category, rank_method) |>
  summarise(
    n_tested   = n(),
    n_fdr05    = sum(fdr        < 0.05, na.rm = TRUE),
    n_fdr05_g  = sum(fdr_global < 0.05, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  print(n = Inf)

log_msg("FDR < 0.05 (within cell type): ", sum(results$fdr        < 0.05, na.rm = TRUE))
log_msg("FDR < 0.05 (global):           ", sum(results$fdr_global < 0.05, na.rm = TRUE))

log_msg("Top 10 gene sets across all cell types:")
results |>
  group_by(cell_type, gs_category, rank_method) |>
  summarise(
    n_tested   = n(),
    n_fdr05    = sum(fdr        < 0.05, na.rm = TRUE),
    n_fdr05_g  = sum(fdr_global < 0.05, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  print(n = Inf)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(results, out_file)
log_msg("Written to ", out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
