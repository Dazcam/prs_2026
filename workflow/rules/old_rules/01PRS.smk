rule get_refs_ldpred2:
  output: zip = "../resources/ldpred2_refs/ldref_hm3_plus.zip",
          rds = "../resources/ldpred2_refs/map_hm3_plus.rds"
  params: gdown = "/shared/home1/c.c1477909/.local/bin/gdown", # requires gdown
          outdir = "../resources/ldpred2_refs/"
  log: "../results/00LOG/00SETUP/get_refs.log"
  shell: """
         {params.gdown} 17dyKGA2PZjMsivlYb_AjDmuZjM1RIGvs -O {params.outdir}
         wget https://ndownloader.figshare.com/files/37802721 -O {params.outdir}map_hm3_plus.rds
         """


rule compute_prs_ldpred2:
    input:  bed = "../results/01PRS/chrALL_final.filt.euro.vcf.bed",
            gwas = "../results/01PRS/scz_hg38_with_neff.tsv",
            hm3 = "../resources/ldpred2_refs/map_hm3_plus.rds"
    output: "../results/01PRS/prs_ldpred2.tsv"
    params: ld_ref_dir = "../resources/ldpred2_refs/",
            geno_prefix = "../results/01PRS/chrALL_final.filt.euro.vcf"
    singularity: "../resources/containers/prs_2026.sif"
    threads: 8
    resources: threads = 8, mem_mb = 80000, time="5:00:00"
    log:    "../results/00LOG/01PRS/compute_prs_ldpred2.log"
    script: "../scripts/prs_compute_ldpred2.R"
