rule format_prs:
    """Reformat PLINK profile output to clean sample x PRS TSV and standardise"""
    input:  "../results/01PRS/prscs_score.profile"
    output: "../results/02PRS_EXPR/prep_input/prs_standardised.tsv"
    resources: time="30:00"
    benchmark: "../results/00LOG/02PRS_EXPR/format_prs.bmk"
    log:    "../results/00LOG/02PRS_EXPR/format_prs.log"
    shell:  "python scripts/02_format_prs.py {input} {output} > {log} 2>&1"

rule filter_expression:
    """Restrict expression to European samples defined by genotype fam file"""
    input:
        expr = "../../eQTL_study_2025/results/05TENSORQTL/prep_input/{cell_type}_quantile.bed",
        fam  = "../results/01PRS/chrALL_final.filt.euro.vcf.fam"
    output: "../results/02PRS_EXPR/prep_input/{cell_type}_exp_filt.tsv"
    resources: time="30:00"
    benchmark: "../results/00LOG/02PRS_EXPR/filter_expression_{cell_type}.bmk"
    log:    "../results/00LOG/02PRS_EXPR/filter_expression_{cell_type}.log"
    shell:
        """
        python scripts/02_filter_expression.py \
            {input.expr} \
            {input.fam} \
            {output} \
        > {log} 2>&1
        """

rule prepare_covariates:
    """
    Extract relevant covariates from eQTL covariate file.
    Keeps sex, PCW, and genotype PCs only — expression PCs are
    not appropriate for PRS-expression analysis as they would
    remove the signal of interest.
    """
    input:
        covs = "../../eQTL_study_2025/results/05TENSORQTL/prep_input/{cell_type}_quantile_genPC_4_expPC_40_split_covariates.txt",
        prs  = "../results/02PRS_EXPR/prep_input/prs_standardised.tsv",
        fam  = "../results/01PRS/chrALL_final.filt.euro.vcf.fam"
    output: "../results/02PRS_EXPR/prep_input/{cell_type}_cov_filt.tsv"
    resources: time="30:00"
    benchmark: "../results/00LOG/02PRS_EXPR/prepare_covariates_{cell_type}.bmk"
    log:    "../results/00LOG/02PRS_EXPR/prepare_covariates_{cell_type}.log"
    shell:
        """
        python scripts/02_prepare_covariates.py \
            {input.covs} \
            {input.prs} \
            {input.fam} \
            {output} \
        > {log} 2>&1
        """

rule run_limma:
    """
    Test association between PRS and gene expression using limma.
    Model: expression ~ PRS_std + Sex + PCW + genPC1..n
    Data are already quantile normalised so voom is not needed.
    """
    input:  expr = rules.filter_expression.output,
            covs = rules.prepare_covariates.output
    output: results = "../results/02PRS_EXPR/limma/{cell_type}_results.tsv",
            rdata   = "../results/02PRS_EXPR/limma/{cell_type}_fit.rds"
    params: cell_type = "{cell_type}"
    resources: mem_mb=16000, time="1:00:00"
    singularity: "../../eQTL_study_2025/resources/containers/seurat5f.sif"
    log:     "../results/00LOG/02PRS_EXPR/limma_{cell_type}.log"
    benchmark: "../results/00LOG/02PRS_EXPR/limma_{cell_type}.bmk"
    script:  "../scripts/02_run_limma_prs.R"

rule combine_limma:
    """
    Combine individual gene limma results across all cell types.
    Adds global FDR correction across all cell types and genes.
    """
    input:
        expand("../results/02PRS_EXPR/limma/{cell_type}_results.tsv",
               cell_type=config['cell_types'])
    output: "../results/02PRS_EXPR/limma/all_cell_types_results.tsv"
    singularity: "../resources/containers/prs_2026.sif"
    resources: time="30:00"
    log:    "../results/00LOG/02PRS_EXPR/combine_limma.log"
    script: "../scripts/02_combine_limma.R"
