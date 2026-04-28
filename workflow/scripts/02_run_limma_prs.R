#--------------------------------------------------------------------------------------
#
#   PRS-Expression association using limma
#
#   Model: expression ~ PRS_std + Sex + PCW + genPC1 + ... + genPCn
#
#   Args:
#     1: filtered expression TSV (TensorQTL format, genes x samples)
#     2: covariate TSV (samples x covariates, includes PRS_std)
#     3: output results TSV
#     4: output RDS (fitted limma model)
#     5: cell type label (for logging)
#
#--------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
})

if (exists("snakemake") && length(snakemake@log) != 0) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush.console()
}

expr_file  <- snakemake@input[["expr"]]
covs_file  <- snakemake@input[["covs"]]
out_res    <- snakemake@output[["results"]]
out_rds    <- snakemake@output[["rdata"]]
cell_type  <- snakemake@params[["cell_type"]]

log_msg("Cell type: ", cell_type)

#--------------------------------------------------------------------------------------
# Load expression
#--------------------------------------------------------------------------------------

log_msg("Loading expression data from ", expr_file)

expr_raw  <- read_tsv(expr_file, show_col_types = FALSE)
meta_cols <- c("#Chr", "start", "end", "TargetID")

gene_ids <- expr_raw$TargetID
expr_mat <- expr_raw |>
  select(-all_of(meta_cols)) |>
  as.matrix()
rownames(expr_mat) <- gene_ids

log_msg("Genes: ", nrow(expr_mat))
log_msg("Samples in expression: ", ncol(expr_mat))

#--------------------------------------------------------------------------------------
# Load covariates
#--------------------------------------------------------------------------------------

log_msg("Loading covariates from ", covs_file)

covs <- read_tsv(covs_file, show_col_types = FALSE) |>
  mutate(sample = as.character(sample))

log_msg("Samples in covariates: ", nrow(covs))
log_msg("Covariates: ", paste(colnames(covs), collapse = ", "))

# Verify PRS_std is present
if (!"PRS_std" %in% colnames(covs)) {
  stop("PRS_std column not found in covariates file")
}

#--------------------------------------------------------------------------------------
# Align samples
#--------------------------------------------------------------------------------------

log_msg("Aligning samples between expression and covariates")

common_samples <- intersect(colnames(expr_mat), covs$sample)
log_msg("Common samples: ", length(common_samples))

if (length(common_samples) < 10) {
  stop(paste("Only", length(common_samples), "common samples — check sample ID format"))
}

# Subset and align — covariate row order must match expression column order
expr_mat <- expr_mat[, common_samples]
covs     <- covs |>
  filter(sample %in% common_samples) |>
  arrange(match(sample, common_samples))

# Critical check — must be true before proceeding
stopifnot(all(colnames(expr_mat) == covs$sample))
log_msg("Sample alignment verified")

#--------------------------------------------------------------------------------------
# Build design matrix
#--------------------------------------------------------------------------------------

log_msg("Building design matrix")

# Identify genotype PC columns
geno_pc_cols <- grep("^genPC[0-9]+$", colnames(covs), value = TRUE)
log_msg("Genotype PCs included: ", paste(geno_pc_cols, collapse = ", "))

# Build formula dynamically from available columns
covariate_terms <- c("Sex", "PCW", geno_pc_cols)
formula_str     <- paste("~ PRS_std +", paste(covariate_terms, collapse = " + "))
log_msg("Model formula: ", formula_str)

design <- model.matrix(as.formula(formula_str), data = covs)

log_msg("Design matrix: ", nrow(design), " samples x ", ncol(design), " terms")

# Check for rank deficiency
if (qr(design)$rank < ncol(design)) {
  warning("Design matrix is rank deficient — check for collinear covariates")
}

#--------------------------------------------------------------------------------------
# Fit limma model
#--------------------------------------------------------------------------------------

log_msg("Fitting limma model")

# Data are already quantile normalised — use lmFit directly
# trend = TRUE fits a mean-variance trend, robust = TRUE downweights outlier genes
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

#--------------------------------------------------------------------------------------
# Extract PRS results
#--------------------------------------------------------------------------------------

log_msg("Extracting PRS_std coefficient results")

results <- topTable(
  fit,
  coef    = "PRS_std",
  number  = Inf,
  sort.by = "P",
  confint = TRUE
) |>
  as_tibble(rownames = "gene") |>
  rename(
    beta     = logFC,
    CI_low   = CI.L,
    CI_high  = CI.R,
    avg_expr = AveExpr,
    t_stat   = t,
    p_value  = P.Value,
    fdr      = adj.P.Val
  ) |>
  mutate(cell_type = cell_type)

#--------------------------------------------------------------------------------------
# QC
#--------------------------------------------------------------------------------------

log_msg("Running QC checks")

# Genomic inflation factor — should be ~1 for well-calibrated results
chi2   <- qchisq(results$p_value, df = 1, lower.tail = FALSE)
lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1)
log_msg("Genomic inflation lambda: ", round(lambda, 3))
if (lambda > 2) {
  warning(paste("High genomic inflation lambda =", round(lambda, 3),
                "— consider additional covariates"))
}

n_sig_fdr <- sum(results$fdr      < 0.05, na.rm = TRUE)
n_sig_nom <- sum(results$p_value  < 0.05, na.rm = TRUE)

log_msg("Genes tested: ", nrow(results))
log_msg("Nominally significant (p<0.05): ", n_sig_nom,
        " (", round(100 * n_sig_nom / nrow(results), 1), "%)")
log_msg("FDR < 0.05: ", n_sig_fdr)

# Top hits
log_msg("Top 5 genes by p-value:")
print(results |> select(gene, symbol, beta, p_value, fdr) |> head(5))

#--------------------------------------------------------------------------------------
# Output
#--------------------------------------------------------------------------------------

log_msg("Writing results to ", out_res)
write_tsv(results, out_res)

log_msg("Saving fitted model to ", out_rds)
saveRDS(fit, out_rds)

log_msg("Done: ", cell_type)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
