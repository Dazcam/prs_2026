localrules: create_prscs_env, get_gene_lookup

# TODO: update container for get_gene_lookup

rule all:
   input:
      "../resources/refs/prs-cs-refs.done",
      "../resources/software/PRS-CS/get_prs_sc.done",
      "../resources/software/PRS-CS/prscs_env.done",
      "../resources/sheets/gene_lookup_hg38.tsv"

rule get_containers:
  input:  "../config/containers.txt"
  output: touch("../resources/containers/get_containers.done")
  params: gdown = "/shared/home1/c.c1477909/.local/bin/gdown" # requires gdown
  log: "../results/00LOG/00SETUP/get_containers.log"
  shell: "../scripts/pull_containers.sh {input}"

rule get_refs:
  output: touch("../resources/refs/prs-cs-refs.done")
  params: ref_dir = "../resources/refs/"
  log: "../results/00LOG/00SETUP/get_refs.log"
  shell: """
         wget -P {params.ref_dir} https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz
         tar -zxvf {params.ref_dir}ldblk_ukbb_eur.tar.gz -C {params.ref_dir}
         rm {params.ref_dir}ldblk_ukbb_eur.tar.gz
         """  

rule get_prs_cs:
  output: touch("../resources/software/PRS-CS/get_prs_sc.done")
  params: out_dir = "../resources/software/PRS-CS/"
  log: "../results/00LOG/00SETUP/get_prs_cs.log"
  shell: """
         git clone https://github.com/getian107/PRScs.git {params.out_dir}
         """

rule create_prscs_env:
    output: touch("../resources/software/PRS-CS/prscs_env.done")
    log:    "../results/00LOG/00SETUP/create_prscs_env.log"
    shell:
        """
        conda env create -f envs/prscs.yaml > {log} 2>&1
        """

rule get_gene_lookup:
     output: ancient("../resources/sheets/gene_lookup_hg38.tsv")
     singularity: "../resources/containers/r_eqtl.sif"
     log:    "../results/00LOG/00SETUP/get_gene_lookup.log"
     run:
        import os
        script  = os.path.abspath("scripts/00_get_gene_lookup.R")
        out     = os.path.abspath(output[0])
        log_abs = os.path.abspath(log[0])
        os.makedirs(os.path.dirname(log_abs), exist_ok=True)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        shell(f"Rscript {script} {out} {log_abs}")

