# Prep data for analysis

### Genotype data
# Filters:
#   >0.05 minor allele frequency (--maf 0.05)
#   exclude individuals with >10% missing data (--mind 0.1)
#   exclude SNPs with missing rate >10% (--geno 0.1)

rule genotype_data:
    input:
        "data/external/call_method_75/call_method_75_TAIR9.csv"
    output:
        config["bfile"] + ".map",
        config["bfile"] + ".bed",
        config["bfile"] + ".bim",
        config["bfile"] + ".fam"
    params:
        indir = "data/raw/call_method_75_TAIR9",
        outdir = config["bfile"]
    run:
        shell("Rscript src/01_convert_to_ped.R")
        shell("plink --file {params.indir} \
        --recode 12 \
        --maf 0.05 \
        --mind 0.1 \
        --geno 0.1 \
        --make-bed \
        --out {params.outdir}")

### Extract TAIR10 annotations
# Note: TAIR10 has same assembly as TAIR9 with updated annotations
rule tair_ann:
    input:
        config["bfile"] + ".map"
    output:
        "data/processed/gene_list_tair10.txt",
        "data/processed/snp_gene_ids_tair10.txt"
    shell:
        "Rscript src/02_extract_gene_ids.R {input}"

### Calculate SNP weightings
# source: dougspeed.com/get-weightings/
# input: plink .bed/.bim/.fam files for all genomic SNPs
# output: partitions/weights.all and partitions/weights.short
# this is for the union of all SNPs (aka all genomic SNPs)
# note: both union of subsets and weightings include all genomic SNPs

rule snp_weightings:
    input:
        bed = config["bfile"] + ".bed"
    output:
        "data/processed/sections/weights.all",
        "data/processed/sections/weights.predictors",
        "data/processed/sections/weights.short"
    params:
        bfile = config["bfile"],
        outdir = "data/processed/sections"
    run:
        shell("{ldak} --cut-weights {params.outdir} \
        --bfile {params.bfile} \
        --section-kb 0.25")
        shell("{ldak} --calc-weights-all {params.outdir} \
        --bfile {params.bfile}")

################################################################################
# Step 2: run principal component analysis (PCA)
# use as covariates in reml model to account for population structure

# rule pca_input:
#     input:
#         bed = config["bfile"] + ".bed",
#         bim = config["bfile"] + ".bim",
#         fam = config["bfile"] + ".fam"
#     output:
#         "models/pca/ld_pruned.bed"
#     params:
#         bfile = config["bfile"],
#         pruned = "models/pca/prune",
#         outdir = "models/pca/ld_pruned"
#     run:
#         shell("""awk < {input.bim} "{{print \$2}}" > data/processed/all.snps""")
#         shell("plink --bfile {params.bfile} \
#         --extract data/processed/all.snps \
#         --make-founders require-2-missing \
#         --indep-pairwise 10 5 .1 \
#         --out {params.pruned}")
#         shell("plink --bfile {params.bfile} \
#         --extract {params.pruned}.prune.in \
#         --out {params.outdir} \
#         --make-bed")
#
# rule pca_result:
#     input:
#         "models/pca/ld_pruned.bed"
#     output:
#         config["pheno_file"]
#     run:
#         shell("Rscript src/02_population_structure.R")
