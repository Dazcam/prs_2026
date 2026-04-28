#--------------------------------------------------------------------------------------
#
#   Stage 3: Gene Set Enrichment Analysis
#
#--------------------------------------------------------------------------------------

rule run_gsea:
    """
    Run fgsea per cell type on limma results.
    Gene sets: GO Biological Process (MSigDB C5 GO:BP)
    Ranking method toggled via params: 't_statistic' or 'neg_log10_pval'.
    Gene sets converted to Ensembl IDs. Optionally protein-coding only.
    """
    input:   limma  = "../results/02PRS_EXPR/limma/{cell_type}_results.tsv",
             lookup = "../resources/sheets/gene_lookup_hg38.tsv"
    output:  results = "../results/03GSEA/{cell_type}_gsea_results.tsv",
             rds     = "../results/03GSEA/{cell_type}_gsea.rds"
    params:  cell_type     = "{cell_type}",
             filter_coding = True,
             min_gs_size   = 15,
             max_gs_size   = 500,
             rank_method   = "t_statistic"
    resources: mem_mb=16000, time="2:00:00"
    singularity: "../resources/containers/prs_2026.sif"
    log:       "../results/00LOG/03GSEA/run_gsea_{cell_type}.log"
    benchmark: "../results/00LOG/03GSEA/run_gsea_{cell_type}.bmk"
    script: "../scripts/03_run_gsea.R"


rule combine_gsea:
    """
    Combine GSEA results across all cell types into a single table.
    Adds global FDR correction across all cell types and gene set categories.
    """
    input:
        expand("../results/03GSEA/{cell_type}_gsea_results.tsv",
               cell_type=config['cell_types'])
    output: "../results/03GSEA/all_cell_types_gsea.tsv"
    resources: mem_mb=8000, time="0:30:00"
    singularity: "../resources/containers/prs_2026.sif"
    log:    "../results/00LOG/03GSEA/combine_gsea.log"
    benchmark: "../results/00LOG/03GSEA/combine_gsea.bmk"
    script: "../scripts/03_combine_gsea.R"

rule render_gsea_report:
    """
    Render R Markdown report for GSEA results.
    One report covering all cell types with tabbed layout.
    """
    input:  gsea  = "../results/03GSEA/all_cell_types_gsea.tsv",
            limma = "../results/02PRS_EXPR/limma/all_cell_types_results.tsv"
    output: "../results/03GSEA/gsea_report.html"
    params: gsea  = "../../results/03GSEA/all_cell_types_gsea.tsv",
            limma = "../../results/02PRS_EXPR/limma/all_cell_types_results.tsv",
            out_file = "../../results/03GSEA/gsea_report.html"
    log:    "../results/00LOG/03GSEA/render_gsea_report.log"
    singularity: "../resources/containers/prs_2026.sif"
    resources: mem_mb=8000, time="0:30:00"
    shell:
        """
	Rscript -e '
          rmarkdown::render(
            input	= "scripts/03_report_gsea.Rmd",
            output_file = normalizePath("{params.out_file}", mustWork = FALSE),
            params	= list(
              gsea_file  = "{params.gsea}",
              limma_file = "{params.limma}",
              fdr_thr    = 0.05
            )
          )
	' > {log} 2>&1
        """
