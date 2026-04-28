rule extract_euro_samples:
    input:  "../resources/sheets/eQTL_study_2025_sample_anns.xlsx"
    output: "../results/01PRS/euro_samples.txt"
    resources: time="10:00"
    log: "../results/00LOG/01PRS/extract_euro_samples.log"
    script: "../scripts/01_extract_euro_samples.py"

rule restrict_to_euro:
    input: vcf = "../../eQTL_study_2025/results/04GENOTYPES-POST/filtered/chrALL_final.filt.vcf.gz",
           samples = "../results/01PRS/euro_samples.txt"
    output: "../results/01PRS/chrALL_final.filt.euro.vcf.gz"
    envmodules: "BCFtools"
    resources: time="1:00:00"
    log: "../results/00LOG/01PRS/restrict_to_euro.log"
    shell: """
           bcftools view \
             --samples-file {input.samples} \
             --output-type z \
             --output {output} \
             {input.vcf} 2>> {log}
           bcftools index --tbi {output} 2>> {log}
           """

rule vcf_to_plink:
    input:  "../results/01PRS/chrALL_final.filt.euro.vcf.gz"
    output: "../results/01PRS/chrALL_final.filt.euro.vcf.bim"
    singularity: "../resources/containers/prs_2026.sif"
    resources: time="1:00:00"
    params: prefix="../results/01PRS/chrALL_final.filt.euro.vcf"
    log:    "../results/00LOG/01PRS/vcf_to_plink.log"
    shell:  """
            plink \
            --vcf {input} \
            --set-missing-var-ids @:# \
            --double-id \
            --make-bed \
            --out {params.prefix} >> {log}
            """

rule format_gwas:
    """Reformat GWAS sumstats to PRS-CS input format: SNP A1 A2 BETA P"""
    input:  "../results/01PRS/scz_hg38_with_neff.tsv"
    output: "../results/01PRS/scz_hg38_prscs_ready.txt"
    log:    "../results/00LOG/01PRS/format_gwas_prscs.log"
    shell:  "python scripts/01_format_gwas.py {input} {output} > {log} 2>&1"

rule run_prscs_auto:
    """
    Run PRS-CS-auto per chromosome.
    phi=None triggers auto mode which estimates the global shrinkage parameter.
    n_gwas uses NEFF from the GWAS file (58,749 — already effective N).
    Thread control prevents numpy from using all cluster cores.
    """
    input:  gwas    = rules.format_gwas.output,
            bim     = rules.vcf_to_plink.output
    output: "../results/01PRS/prscs/scz_pst_eff_a1_b0.5_phiauto_chr{chr}.txt"
    params: bim_prefix = "../results/01PRS/chrALL_final.filt.euro.vcf",
            ref_dir    = "../resources/refs/ldblk_ukbb_eur",
            prscs_dir  = "../resources/software/PRS-CS",
            out_prefix = "../results/01PRS/prscs/scz",
            n_gwas     = 58749
    resources: mem_mb=40000, time="3:00:00", partition = "htc_genoa", threads = 4
    conda: "../envs/prscs.yaml"
    threads: 4
    benchmark: "../results/00LOG/01PRS/run_prscs_chr{chr}.bmk"
    log:    "../results/00LOG/01PRS/run_prscs_chr{chr}.log"
    shell:
         """
        export MKL_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}

        python {params.prscs_dir}/PRScs.py \
            --ref_dir={params.ref_dir} \
            --bim_prefix={params.bim_prefix} \
            --sst_file={input.gwas} \
            --n_gwas={params.n_gwas} \
            --chrom={wildcards.chr} \
            --out_dir={params.out_prefix} \
            --seed=42 \
            > {log} 2>&1
        """

rule combine_prscs_weights:
    """Concatenate per-chromosome PRS-CS effect size files into one genome-wide file"""
    input:  expand(rules.run_prscs_auto.output, chr=range(1,23))
    output: "../results/01PRS/prscs/scz_prscs_weights_all.txt"
    log:    "../results/00LOG/01PRS/combine_prscs.log"
    benchmark: "../results/00LOG/01PRS/combine_prscs_weights.bmk"
    shell:
        """
        # Add header then concatenate all chromosomes
        echo -e "CHR\tSNP\tBP\tA1\tA2\tBETA" > {output}
        cat {input} >> {output}
        echo "Total variants in weight file: $(wc -l < {output})" > {log}
        """

rule score_prscs:
    """
    Compute per-sample PRS using PLINK --score with PRS-CS posterior effect sizes.
    Column positions: 2=SNP, 4=A1, 7=BETA_EST
    --score-col-nums not needed for plink1.9 — use positional args.
    """
    input:  weights = "../results/01PRS/prscs/scz_prscs_weights_all.txt",
            bim     = "../results/01PRS/chrALL_final.filt.euro.vcf.bim",
    output: "../results/01PRS/prscs_score.profile"
    params: bim_prefix = "../results/01PRS/chrALL_final.filt.euro.vcf",
            out_prefix = "../results/01PRS/prscs_score"
    resources: mem_mb=40000, time="3:00:00", partition = "htc_genoa", threads = 4
    singularity: "../resources/containers/prs_2026.sif"
    threads: 4
    benchmark: "../results/00LOG/01PRS/score_prscs.bmk"
    log:    "logs/01PRS/score_prscs.log"
    shell:
        """
        plink \
            --bfile {params.bim_prefix} \
            --score {input.weights} 2 4 6 header \
            --out {params.out_prefix} \
            --threads {threads} \
        > {log} 2>&1
        """

