#--------------------------------------------------------------------------------------
#
#   Combine limma results across all cell types
#   Adds global FDR correction across all cell types and genes
#
#--------------------------------------------------------------------------------------

if (exists("snakemake") && length(snakemake@log) != 0) {
  dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

suppressPackageStartupMessages({
  library(tidyverse)
})

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush.console()
}

input_files <- snakemake@input
out_file    <- snakemake@output[[1]]

log_msg("Combining limma results across ", length(input_files), " cell types")

results <- map_dfr(input_files, function(f) {
  log_msg("  Reading: ", f)
  read_tsv(f, show_col_types = FALSE)
})

log_msg("Total gene tests: ", nrow(results))
log_msg("Cell types: ", n_distinct(results$cell_type))

# Global FDR across all cell types and genes
results <- results |>
  mutate(fdr_global = p.adjust(p_value, method = "BH"))

log_msg("FDR < 0.05 (within cell type): ", sum(results$fdr < 0.05, na.rm = TRUE))
log_msg("FDR < 0.05 (global):           ", sum(results$fdr_global < 0.05, na.rm = TRUE))

log_msg("Top 10 genes across all cell types:")
results |>
  arrange(p_value) |>
  select(cell_type, gene, beta, p_value, fdr, fdr_global) |>
  head(10) |>
  print()

write_tsv(results, out_file)
log_msg("Written to ", out_file)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
