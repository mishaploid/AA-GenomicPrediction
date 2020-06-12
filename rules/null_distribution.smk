################################################################################
# Create an empirical null distribution
# AKA the most time consuming step...
# Random SNP sets need to be generated in advance and stored in separate folders

## generate a uniform distribution of SNP sizes
rule null_sampling:
    output:
        "data/interim/null_group_sizes.txt"
    shell:
        "Rscript src/07_null_group_sizes.R"

# sample gene groups
# calls a custom R script with a function to sample all SNPs in a gene
# stops once a target threshold of SNPs is achieved (output from null_sampling)

rule null_gene_groups:
    input:
        "data/processed/snp_gene_ids_tair10.txt"
    output:
        "data/processed/random_sets/null_{null}.txt"
    params:
        num = "{null}"
    shell:
        "Rscript src/08_generate_null_gene_groups.R {params.num}"

### calculate kinships for control sets
### Only interested in partitioning variance
### Ignore weightings and power = 0
# need to improve efficiency of this step...
# consider looping - e.g. for r in range(NULL+49):

rule calc_kins_null:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam"
    output:
        list = expand("data/processed/random_sets/c_{random}/partition.list", random = 5000),
        k1 = expand("data/processed/random_sets/c_{random}/kinships.1.grm.details", random = 5000),
        k2 = expand("data/processed/random_sets/c_{random}/kinships.2.grm.details", random = 5000)
    params:
        bfile = config["bfile"],
        prefix = "data/processed/random_sets/c_"
    run:
        for r in RANDOM:
            shell("{ldak} --cut-kins {params.prefix}{r} \
            --bfile {params.bfile} \
            --partition-number 2 \
            --partition-prefix {params.prefix}{r}/list")
            shell("{ldak} --calc-kins {params.prefix}{r} \
            --bfile {params.bfile} \
            --partition 1 \
            --ignore-weights YES \
            --power 0")
            shell("{ldak} --calc-kins {params.prefix}{r} \
            --bfile {params.bfile} \
            --partition 2 \
            --ignore-weights YES \
            --power 0")

### REML model for control sets

rule null_h2:
    input:
        pheno = config["pheno_file"],
        mgrm = expand("data/processed/random_sets/c_{random}/partition.list", random = 5000),
        covar = config["covar"]
    output:
        out = expand("models/null_h2/c_{random}/{{trait}}.h2.reml", random = RANDOM),
    params:
        prefix = "models/null_h2/c_",
        trait = "{trait}",
        mgrm = "data/processed/random_sets/c_"
    run:
        for r in RANDOM:
            shell("{ldak} --reml {params.prefix}{r}/{params.trait}.h2 \
            --pheno {input.pheno} \
            --pheno-name {params.trait} \
            --mgrm {params.mgrm}{r}/partition.list \
            --covar {input.covar} \
            --constrain YES \
            --reml-iter 500")

rule null_results:
    input:
        expand("models/null_h2/c_{random}/{trait}.h2.reml", random = 5000, trait = TRAIT),
        gblup = "reports/gblup.RData"
    output:
        "reports/lr_null_results.csv"
    run:
        shell("Rscript src/09_summarize_null.R")

### rinse and repeat for each pathway
## sample 1000 random gene groups of a similar size for each pathway
# sample gene groups for each pathway tested

rule null_pathway_sampling:
    input:
        "data/processed/snp_gene_ids_tair10.txt",
        "data/processed/pathways/{pathway}/{pathway}.txt"
    output:
        "data/processed/random_sets_pathways/{pathway}/null_1000.txt"
    params:
        pathway = "{pathway}"
    shell:
        "Rscript src/08a_generate_null_gene_groups_pathways.R {params.pathway}"

### calculate kinships for control sets
### Only interested in partitioning variance
### Ignore weightings and power = 0

rule calc_kins_null_pathways:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam",
        null_snps = "data/processed/random_sets_pathways/{pathway}/null_1000.txt"
    output:
        list = "data/processed/random_sets_pathways/{pathway}/c_1000/partition.list",
        k1 = "data/processed/random_sets_pathways/{pathway}/c_1000/kinships.1.grm.details",
        k2 = "data/processed/random_sets_pathways/{pathway}/c_1000/kinships.2.grm.details"
    params:
        bfile = config["bfile"],
        prefix = "data/processed/random_sets_pathways/{pathway}/c_"
    run:
        for r in range(1,1001):
            shell("{ldak} --cut-kins {params.prefix}{r} \
            --bfile {params.bfile} \
            --partition-number 2 \
            --partition-prefix {params.prefix}{r}/list")
            shell("{ldak} --calc-kins {params.prefix}{r} \
            --bfile {params.bfile} \
            --partition 1 \
            --ignore-weights YES \
            --power 0")
            shell("{ldak} --calc-kins {params.prefix}{r} \
            --bfile {params.bfile} \
            --partition 2 \
            --ignore-weights YES \
            --power 0")

### REML model for control sets

rule null_h2_pathways:
    input:
        pheno = config["pheno_file"],
        mgrm = "data/processed/random_sets_pathways/{pathway}/c_1000/partition.list",
        covar = config["covar"]
    output:
        out = expand("models/null_h2_pathways/{{pathway}}/c_{random}/{{trait}}.h2.reml", random = range(1,1001)),
    params:
        prefix = "models/null_h2_pathways/{pathway}/c_",
        trait = "{trait}",
        mgrm = "data/processed/random_sets_pathways/{pathway}/c_"
    run:
        for r in range(1,1001):
            shell("{ldak} --reml {params.prefix}{r}/{params.trait}.h2 \
            --pheno {input.pheno} \
            --pheno-name {params.trait} \
            --mgrm {params.mgrm}{r}/partition.list \
            --covar {input.covar} \
            --constrain YES \
            --reml-iter 500")
